#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    getwd()
  }
})
OXYGEN_ROOT <- normalizePath(file.path(.script_dir, ".."), mustWork = TRUE)
WORKFLOW_ROOT <- normalizePath(
  file.path(OXYGEN_ROOT, "code", "O2_supply_demand_MAP"),
  mustWork = TRUE
)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R"), local = environment())
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_invitro_utils.R"),
  local = environment(),
  chdir = TRUE
)
rm(.script_dir)

parse_numeric_csv <- function(value, default) {
  if (is.null(value) || !length(value) || !nzchar(trimws(as.character(value[[1L]])))) {
    return(as.numeric(default))
  }
  out <- suppressWarnings(as.numeric(strsplit(as.character(value[[1L]]), ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!is.finite(out))) stop("Grid values must be finite comma-separated numbers.")
  unique(out)
}

first_value <- function(x, default) {
  if (is.null(x) || !length(x)) default else x[[1L]]
}

segment_distribution_metrics <- function(run, segment_id, thresholds = c(70, 80)) {
  hit <- which(vapply(
    run$segment_results,
    function(x) identical(as.character(x$segment$segment_id), as.character(segment_id)),
    logical(1)
  ))
  if (length(hit) != 1L) {
    return(c(mean_N = NA_real_, mass_ge70 = NA_real_, mass_ge80 = NA_real_))
  }
  selected <- run$segment_results[[hit]]$selection$selected_frac
  probs <- pmax(as.numeric(selected), 0)
  total <- sum(probs)
  if (!is.finite(total) || total <= 0) {
    return(c(mean_N = NA_real_, mass_ge70 = NA_real_, mass_ge80 = NA_real_))
  }
  probs <- probs / total
  grid <- as.numeric(run$grid_pre)
  c(
    mean_N = sum(grid * probs),
    mass_ge70 = sum(probs[grid >= thresholds[[1L]]]),
    mass_ge80 = sum(probs[grid >= thresholds[[2L]]])
  )
}

terminal_segment_id <- function(run, lineage_id) {
  matching <- Filter(
    function(x) identical(
      as.character(x$segment$lineage_id),
      as.character(lineage_id)
    ),
    run$segment_results
  )
  if (!length(matching)) return(NA_character_)
  passage_index <- vapply(
    matching,
    function(x) as.integer(x$segment$lineage_passage_index),
    integer(1)
  )
  as.character(matching[[which.max(passage_index)]]$segment$segment_id)
}

buffer_survival <- function(N, smax, beta, n_exp) {
  smax * exp(-beta * (44 / N)^n_exp)
}

main <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- first_value(argv$fit_dir, NULL)
  if (is.null(fit_dir)) {
    stop("Usage: scan_invitro_misseg_buffer_grid.R --fit_dir=/abs/seed --out_dir=/abs/output")
  }
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  fit_result_path <- normalizePath(
    first_value(argv$fit_result, file.path(fit_dir, "fit_result.rds")),
    mustWork = TRUE
  )
  out_dir <- normalizePath(
    first_value(argv$out_dir, file.path(fit_dir, "postfit_misseg_buffer_grid")),
    mustWork = FALSE
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p_misseg_values <- parse_numeric_csv(argv$p_misseg_values, c(0.0075, 0.010, 0.015, 0.020))
  buffer_beta_values <- parse_numeric_csv(argv$buffer_beta_values, c(0.3, 0.5, 0.7, 1.0))
  buffer_n_exp_values <- parse_numeric_csv(argv$buffer_n_exp_values, c(4, 6, 8))
  buffer_smax <- as.numeric(first_value(argv$buffer_smax, 0.98))
  p_wgd <- as.numeric(first_value(argv$p_wgd, 1e-5))
  n_cores <- max(1L, as.integer(first_value(argv$n_cores, 1L)))
  if (!is.finite(buffer_smax) || buffer_smax < 0 || buffer_smax > 1) {
    stop("buffer_smax must lie in [0, 1].")
  }
  if (!is.finite(p_wgd) || p_wgd <= 0) stop("p_wgd must be strictly positive.")

  fit_result <- readRDS(fit_result_path)
  cfg <- fit_result$cfg
  if (is.null(cfg)) stop("fit_result.rds does not contain cfg.")
  base_params <- fit_result$best_params
  if (is.null(base_params)) stop("fit_result.rds does not contain best_params.")
  fit_objects <- ivt_load_fit_objects(
    repo_root = OXYGEN_ROOT,
    fit_objects_dir = normalizePath(fit_result$fit_objects_dir, mustWork = TRUE),
    flow_csv_path = normalizePath(fit_result$flow_density_path, mustWork = TRUE),
    death_data_path = normalizePath(fit_result$death_data_path, mustWork = TRUE)
  )
  death_weight <- as.numeric(first_value(fit_result$death_weight, 1))
  sigma_death_logit <- as.numeric(first_value(fit_result$sigma_death_logit, 0.75))
  death_fraction_eps <- as.numeric(first_value(fit_result$death_fraction_eps, 1e-4))
  passage_time_weight <- as.numeric(first_value(fit_result$passage_time_weight, 0.25))
  passage_time_tolerance_days <- as.numeric(first_value(fit_result$passage_time_tolerance_days, 1))
  passage_time_sigma_days <- as.numeric(first_value(fit_result$passage_time_sigma_days, 1))
  passage_time_df <- as.numeric(first_value(fit_result$passage_time_df, 4))

  grid <- expand.grid(
    p_misseg = p_misseg_values,
    buffer_beta = buffer_beta_values,
    buffer_n_exp = buffer_n_exp_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$grid_id <- sprintf("grid_%03d", seq_len(nrow(grid)))

  evaluate_row <- function(i) {
    row <- grid[i, , drop = FALSE]
    run_params <- base_params
    run_params$p_misseg <- row$p_misseg[[1L]]
    run_params$buffer_smax <- buffer_smax
    run_params$buffer_beta <- row$buffer_beta[[1L]]
    run_params$buffer_n_exp <- row$buffer_n_exp[[1L]]
    run_params$p_wgd <- p_wgd
    comp <- tryCatch(
      ivt_objective_components(
        run_params = run_params,
        fit_objects = fit_objects,
        cfg = cfg,
        fallback_max_passage_days = 14,
        growth_weight = 1,
        ploidy_weight = 1,
        flow_weight = 1,
        death_weight = death_weight,
        passage_time_weight = passage_time_weight,
        passage_time_tolerance_days = passage_time_tolerance_days,
        passage_time_sigma_days = passage_time_sigma_days,
        passage_time_df = passage_time_df,
        sigma_death_logit = sigma_death_logit,
        death_fraction_eps = death_fraction_eps
      ),
      error = function(e) e
    )
    base <- data.frame(
      grid_id = row$grid_id,
      p_misseg = row$p_misseg,
      buffer_smax = buffer_smax,
      buffer_beta = row$buffer_beta,
      buffer_n_exp = row$buffer_n_exp,
      p_wgd = p_wgd,
      survival_N44_m1 = buffer_survival(44, buffer_smax, row$buffer_beta, row$buffer_n_exp),
      survival_N88_m1 = buffer_survival(88, buffer_smax, row$buffer_beta, row$buffer_n_exp),
      stringsAsFactors = FALSE
    )
    if (inherits(comp, "error")) {
      return(cbind(
        base,
        data.frame(
          status = "ERROR",
          error = conditionMessage(comp),
          objective_likelihood = NA_real_,
          growth_loglik = NA_real_,
          ploidy_loglik = NA_real_,
          flow_loglik = NA_real_,
          death_loglik = NA_real_,
          passage_time_loglik = NA_real_,
          protocol_feasibility_status = NA_character_,
          n_insufficient_boundaries = NA_integer_,
          mean_2N_O1_A19 = NA_real_,
          mean_2N_O2_A19 = NA_real_,
          pooled_mean_2N_A19 = NA_real_,
          pooled_mass_ge70_2N_A19 = NA_real_,
          pooled_mass_ge80_2N_A19 = NA_real_,
          mean_4N_O1_A17 = NA_real_,
          mean_4N_O2_A17 = NA_real_,
          pooled_mean_4N_A17 = NA_real_,
          mean_2N_C_A12 = NA_real_,
          mean_4N_C_A12 = NA_real_,
          stringsAsFactors = FALSE
        )
      ))
    }

    m_2n_o1 <- segment_distribution_metrics(comp$run_2N, "2N-O1-A19")
    m_2n_o2 <- segment_distribution_metrics(comp$run_2N, "2N-O2-A19")
    m_4n_o1 <- segment_distribution_metrics(comp$run_4N, "4N-O1-A17")
    m_4n_o2 <- segment_distribution_metrics(comp$run_4N, "4N-O2-A17")
    m_2n_c <- segment_distribution_metrics(
      comp$run_2N,
      terminal_segment_id(comp$run_2N, "C")
    )
    m_4n_c <- segment_distribution_metrics(
      comp$run_4N,
      terminal_segment_id(comp$run_4N, "C")
    )
    cbind(
      base,
      data.frame(
        status = "OK",
        error = "",
        objective_likelihood = comp$objective,
        growth_loglik = comp$growth_loglik,
        ploidy_loglik = comp$ploidy_loglik,
        flow_loglik = comp$flow_loglik,
        death_loglik = comp$death_loglik,
        passage_time_loglik = comp$passage_time_loglik,
        protocol_feasibility_status = comp$protocol_feasibility_status,
        n_insufficient_boundaries = comp$n_insufficient_boundaries,
        mean_2N_O1_A19 = m_2n_o1[["mean_N"]],
        mean_2N_O2_A19 = m_2n_o2[["mean_N"]],
        pooled_mean_2N_A19 = mean(c(m_2n_o1[["mean_N"]], m_2n_o2[["mean_N"]])),
        pooled_mass_ge70_2N_A19 = mean(c(m_2n_o1[["mass_ge70"]], m_2n_o2[["mass_ge70"]])),
        pooled_mass_ge80_2N_A19 = mean(c(m_2n_o1[["mass_ge80"]], m_2n_o2[["mass_ge80"]])),
        mean_4N_O1_A17 = m_4n_o1[["mean_N"]],
        mean_4N_O2_A17 = m_4n_o2[["mean_N"]],
        pooled_mean_4N_A17 = mean(c(m_4n_o1[["mean_N"]], m_4n_o2[["mean_N"]])),
        mean_2N_C_A12 = m_2n_c[["mean_N"]],
        mean_4N_C_A12 = m_4n_c[["mean_N"]],
        stringsAsFactors = FALSE
      )
    )
  }

  indices <- seq_len(nrow(grid))
  rows <- if (n_cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(indices, evaluate_row, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(indices, evaluate_row)
  }
  summary <- dplyr::bind_rows(rows)
  summary$abs_pooled_mean_error_2N_A19 <- abs(summary$pooled_mean_2N_A19 - 77.45)
  summary$abs_pooled_mass_ge70_error_2N_A19 <- abs(summary$pooled_mass_ge70_2N_A19 - 0.625)
  summary <- summary[order(
    summary$status != "OK",
    summary$abs_pooled_mean_error_2N_A19,
    summary$abs_pooled_mass_ge70_error_2N_A19
  ), , drop = FALSE]
  utils::write.table(
    summary,
    file.path(out_dir, "invitro_misseg_buffer_grid_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
  manifest <- data.frame(
    key = c(
      "fit_dir", "fit_result", "out_dir", "p_misseg_values",
      "buffer_beta_values", "buffer_n_exp_values", "buffer_smax", "p_wgd",
      "n_cores", "grid_size", "generated_at"
    ),
    value = c(
      fit_dir,
      fit_result_path,
      normalizePath(out_dir, mustWork = TRUE),
      paste(p_misseg_values, collapse = ","),
      paste(buffer_beta_values, collapse = ","),
      paste(buffer_n_exp_values, collapse = ","),
      buffer_smax,
      p_wgd,
      n_cores,
      nrow(grid),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    manifest,
    file.path(out_dir, "invitro_misseg_buffer_grid_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message("In-vitro missegregation/buffer grid written to: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
