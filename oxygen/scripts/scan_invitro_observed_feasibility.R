#!/usr/bin/env Rscript

# Two-stage observed-data feasibility scan for the in-vitro model.
#
# Stage 1 calibrates the largest contiguous p_mis_base value for which the
# predicted 2N-C-A12 flow mass at ploidy >= 3 stays below the requested
# normoxic limit. Stage 2 samples the calibrated nine-dimensional parameter
# box with a deterministic Latin hypercube and evaluates each candidate once.
# This is a diagnostic scan: it does not change the fitted objective and it
# does not introduce functional soft constraints.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  file_path <- sub("^--file=", "", file_arg)
  file_path <- file_path[
    basename(file_path) == "scan_invitro_observed_feasibility.R"
  ]
  if (length(file_path)) {
    dirname(normalizePath(file_path[[1L]], mustWork = TRUE))
  } else {
    getwd()
  }
})

.configured_project_root <- Sys.getenv("O2SD_PROJECT_ROOT", unset = "")
if (nzchar(.configured_project_root)) {
  PROJECT_ROOT <- normalizePath(.configured_project_root, mustWork = TRUE)
  OXYGEN_ROOT <- normalizePath(file.path(PROJECT_ROOT, "oxygen"), mustWork = TRUE)
} else {
  OXYGEN_ROOT <- normalizePath(file.path(.script_dir, ".."), mustWork = TRUE)
  PROJECT_ROOT <- normalizePath(file.path(OXYGEN_ROOT, ".."), mustWork = TRUE)
}
WORKFLOW_ROOT <- normalizePath(
  file.path(OXYGEN_ROOT, "code", "O2_supply_demand_MAP"),
  mustWork = TRUE
)

source(
  file.path(.script_dir, "invitro_observed_feasibility_utils.R"),
  local = environment()
)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R"), local = environment())
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_invitro_utils.R"),
  local = environment(),
  chdir = TRUE
)

`%||%` <- if (exists("%||%", inherits = TRUE)) {
  get("%||%", inherits = TRUE)
} else {
  function(x, y) if (is.null(x) || !length(x)) y else x
}

.ivt_scan_first <- function(x, default = NULL) {
  if (is.null(x) || !length(x)) default else x[[1L]]
}

.ivt_scan_as_integer <- function(x, name, minimum = 1L) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be one integer >= ", minimum, ".", call. = FALSE)
  }
  value
}

.ivt_scan_as_number <- function(x, name, lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || !is.finite(value) || value < lower || value > upper) {
    stop(
      name, " must be one finite number in [", lower, ", ", upper, "].",
      call. = FALSE
    )
  }
  value
}

.ivt_scan_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(isTRUE(default))
  value <- tolower(trimws(as.character(x[[1L]])))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Expected a boolean value, received: ", x[[1L]], call. = FALSE)
}

.ivt_scan_resolve_input <- function(explicit,
                                    recorded,
                                    fallback,
                                    label,
                                    directory = FALSE) {
  explicit_present <- !is.null(explicit) && length(explicit) &&
    !is.na(explicit[[1L]]) && nzchar(as.character(explicit[[1L]]))
  if (explicit_present) {
    explicit_path <- as.character(explicit[[1L]])
    explicit_exists <- if (directory) dir.exists(explicit_path) else file.exists(explicit_path)
    if (!explicit_exists) {
      stop("Explicit ", label, " does not exist: ", explicit_path, call. = FALSE)
    }
    return(normalizePath(explicit_path, mustWork = TRUE))
  }
  candidates <- unique(Filter(
    function(x) !is.null(x) && length(x) && !is.na(x[[1L]]) && nzchar(as.character(x[[1L]])),
    list(explicit, recorded, fallback)
  ))
  for (candidate in candidates) {
    path <- as.character(candidate[[1L]])
    exists_here <- if (directory) dir.exists(path) else file.exists(path)
    if (exists_here) return(normalizePath(path, mustWork = TRUE))
  }
  rendered <- vapply(candidates, function(x) as.character(x[[1L]]), character(1))
  stop(
    "Unable to resolve ", label, ". Checked: ", paste(rendered, collapse = "; "),
    call. = FALSE
  )
}

.ivt_scan_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
  invisible(path)
}

.ivt_scan_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    line <- system2(sha256sum, path, stdout = TRUE, stderr = FALSE)
    return(strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]])
  }
  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    line <- system2(shasum, c("-a", "256", path), stdout = TRUE, stderr = FALSE)
    return(strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]])
  }
  NA_character_
}

.ivt_scan_git_value <- function(args) {
  tryCatch(
    paste(system2("git", c("-C", PROJECT_ROOT, args), stdout = TRUE, stderr = FALSE), collapse = " | "),
    error = function(e) NA_character_
  )
}

.ivt_scan_file_inventory <- function(role, paths) {
  paths <- sort(unique(normalizePath(paths, mustWork = TRUE)))
  info <- file.info(paths)
  data.frame(
    role = role,
    path = paths,
    size_bytes = as.numeric(info$size),
    sha256 = vapply(paths, .ivt_scan_sha256, character(1)),
    stringsAsFactors = FALSE
  )
}

.ivt_scan_build_run_contract <- function(fit_result_path,
                                         parameter_table_path,
                                         fit_objects_dir,
                                         flow_density_path,
                                         death_data_path,
                                         settings) {
  scan_files <- c(
    file.path(.script_dir, "scan_invitro_observed_feasibility.R"),
    file.path(.script_dir, "invitro_observed_feasibility_utils.R")
  )
  workflow_files <- c(
    list.files(
      file.path(WORKFLOW_ROOT, "model"),
      pattern = "\\.(R|cpp)$",
      full.names = TRUE,
      recursive = TRUE
    ),
    list.files(
      file.path(WORKFLOW_ROOT, "util"),
      pattern = "\\.R$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
  fit_object_files <- list.files(
    fit_objects_dir,
    full.names = TRUE,
    recursive = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  fit_object_files <- fit_object_files[file.info(fit_object_files)$isdir %in% FALSE]
  inputs <- rbind(
    .ivt_scan_file_inventory("scan_code", scan_files),
    .ivt_scan_file_inventory("workflow_code", workflow_files),
    .ivt_scan_file_inventory("fit_result", fit_result_path),
    .ivt_scan_file_inventory("parameter_table", parameter_table_path),
    .ivt_scan_file_inventory("fit_object", fit_object_files),
    .ivt_scan_file_inventory("flow_density", flow_density_path),
    .ivt_scan_file_inventory("death_data", death_data_path)
  )
  settings <- data.frame(
    key = names(settings),
    value = vapply(settings, function(x) paste(as.character(x), collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
  list(schema_version = 1L, inputs = inputs, settings = settings)
}

.ivt_scan_prepare_run_contract <- function(contract, out_dir, resume) {
  rds_path <- file.path(out_dir, "run_contract.rds")
  tsv_path <- file.path(out_dir, "run_contract_inputs.tsv")
  settings_path <- file.path(out_dir, "run_contract_settings.tsv")
  reusable_paths <- c(
    file.path(out_dir, "p_mis_base_calibration.tsv"),
    file.path(out_dir, "p_mis_base_calibration_summary.tsv"),
    file.path(out_dir, "p_mis_base_calibration.rds"),
    file.path(out_dir, "feasibility_design.rds"),
    list.files(file.path(out_dir, "checkpoints"), pattern = "\\.rds$", full.names = TRUE)
  )
  has_reusable <- any(file.exists(reusable_paths))
  if (!isTRUE(resume) && (file.exists(rds_path) || has_reusable)) {
    stop(
      "resume=FALSE requires a new output directory; existing run artifacts were found.",
      call. = FALSE
    )
  }
  if (isTRUE(resume) && file.exists(rds_path)) {
    prior <- readRDS(rds_path)
    if (!identical(prior, contract)) {
      stop(
        "Run contract differs from existing output; refusing to reuse calibration/checkpoints.",
        call. = FALSE
      )
    }
  } else if (isTRUE(resume) && has_reusable) {
    stop(
      "Reusable outputs exist without run_contract.rds; refusing an unverified resume.",
      call. = FALSE
    )
  } else {
    saveRDS(contract, rds_path)
  }
  .ivt_scan_write_tsv(contract$inputs, tsv_path)
  .ivt_scan_write_tsv(contract$settings, settings_path)
  invisible(rds_path)
}

.ivt_scan_objective_args <- function(run_params, context = .ivt_scan_context) {
  list(
    run_params = run_params,
    fit_objects = context$fit_objects,
    cfg = context$cfg,
    fallback_max_passage_days = 14,
    growth_weight = 1,
    ploidy_weight = context$ploidy_weight,
    flow_weight = 1,
    death_weight = context$death_weight,
    passage_time_weight = context$passage_time_weight,
    passage_time_tolerance_days = context$passage_time_tolerance_days,
    passage_time_sigma_days = context$passage_time_sigma_days,
    passage_time_df = context$passage_time_df,
    sigma_death_logit = context$sigma_death_logit,
    death_fraction_eps = context$death_fraction_eps
  )
}

.ivt_scan_metric_summary <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(mean = NA_real_, rmse = NA_real_, median = NA_real_, q05 = NA_real_, min = NA_real_))
  }
  c(
    mean = mean(x),
    rmse = sqrt(mean(x^2)),
    median = stats::median(x),
    q05 = as.numeric(stats::quantile(x, 0.05, names = FALSE, type = 8)),
    min = min(x)
  )
}

.ivt_scan_observed_kary <- function(run, fit_data, segment_id) {
  ids <- vapply(
    run$segment_results,
    function(x) as.character(x$segment$segment_id %||% ""),
    character(1)
  )
  idx <- which(ids == segment_id)
  if (length(idx) != 1L) stop("Expected one segment for observed karyotype: ", segment_id)
  segment <- run$segment_results[[idx]]$segment
  passage_ids <- .ivt_segment_endpoint_data_ids(run, segment)
  values <- unlist(lapply(passage_ids, function(pid) {
    if (!pid %in% names(fit_data)) return(numeric())
    as.numeric(ivt_observed_passage_summary(fit_data[[pid]])$observed_kary)
  }), use.names = FALSE)
  values[is.finite(values)]
}

.ivt_scan_density_mass <- function(df, lower = 3, upper = Inf) {
  grid <- .ivt_scan_density_grid(df)
  keep <- grid$x >= lower & grid$x < upper
  sum(grid$y[keep]) * grid$dx / grid$total
}

.ivt_scan_collect_flow_metrics <- function(overlay, candidate_id, stage) {
  required <- c("segment_id", "passage_id", "sample_name", "series", "ploidy", "density")
  if (!all(required %in% names(overlay))) stop("flow_overlay_df has an incompatible schema.")
  keep_segment <- as.character(overlay$segment_id) == "2N-C-A12" |
    grepl("^2N-O[12]-A(12|18|23)$", as.character(overlay$segment_id))
  use <- overlay[keep_segment, , drop = FALSE]
  keys <- unique(use[c("segment_id", "passage_id", "sample_name", "series")])
  if (!nrow(keys)) return(data.frame())
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    hit <- as.character(use$segment_id) == as.character(key$segment_id[[1L]]) &
      as.character(use$passage_id) == as.character(key$passage_id[[1L]]) &
      as.character(use$sample_name) == as.character(key$sample_name[[1L]]) &
      as.character(use$series) == as.character(key$series[[1L]])
    z <- use[hit, , drop = FALSE]
    shape <- ivt_scan_flow_shape_metrics(z)
    cbind(
      data.frame(
        candidate_id = candidate_id,
        stage = stage,
        segment_id = as.character(key$segment_id[[1L]]),
        passage_id = as.character(key$passage_id[[1L]]),
        sample_name = as.character(key$sample_name[[1L]]),
        series = as.character(key$series[[1L]]),
        high_mass_ge3 = .ivt_scan_density_mass(z, 3, Inf),
        stringsAsFactors = FALSE
      ),
      shape
    )
  })
  do.call(rbind, rows)
}

.ivt_scan_collect_kary_metrics <- function(comp, candidate_id, stage, sigma_kary) {
  segment_ids <- c("2N-O1-A19", "2N-O2-A19", "4N-O1-A17", "4N-O2-A17")
  rows <- lapply(segment_ids, function(segment_id) {
    run <- if (startsWith(segment_id, "2N-")) comp$run_2N else comp$run_4N
    latent <- ivt_scan_segment_distribution_metrics(run, segment_id)
    smoothed <- ivt_scan_segment_distribution_metrics(run, segment_id, sigma_kary = sigma_kary)
    rbind(latent, smoothed)
  })
  out <- do.call(rbind, rows)
  cbind(
    data.frame(candidate_id = candidate_id, stage = stage, stringsAsFactors = FALSE),
    out
  )
}

.ivt_scan_pooled_metric <- function(kary, cohort, field, basis = "latent") {
  keep <- as.character(kary$cohort) == cohort &
    as.character(kary$distribution_basis) == basis
  values <- suppressWarnings(as.numeric(kary[[field]][keep]))
  if (length(values) != 2L || any(!is.finite(values))) return(NA_real_)
  mean(values)
}

.ivt_scan_flow_passage_mean <- function(flow, passage, field, series = "Predicted") {
  keep <- as.character(flow$series) == series &
    grepl(paste0("^2N-O[12]-A", passage, "$"), as.character(flow$segment_id))
  values <- suppressWarnings(as.numeric(flow[[field]][keep]))
  if (!length(values) || any(!is.finite(values))) return(NA_real_)
  mean(values)
}

.ivt_scan_terminal_event_metrics <- function(summary_df, cohort) {
  use <- summary_df[as.character(summary_df$cohort) == cohort, , drop = FALSE]
  if (!nrow(use)) {
    return(c(gross_divisions = NA_real_, hypoxia_deaths = NA_real_, buffer_inflow = NA_real_, nonlive_inflow = NA_real_))
  }
  split_rows <- split(seq_len(nrow(use)), as.character(use$scenario_id))
  terminal <- do.call(rbind, lapply(split_rows, function(idx) {
    use[idx[[which.max(as.numeric(use$lineage_passage_index[idx]))]], , drop = FALSE]
  }))
  sums <- function(name) {
    values <- suppressWarnings(as.numeric(terminal[[name]]))
    if (!length(values) || any(!is.finite(values))) NA_real_ else sum(values)
  }
  c(
    gross_divisions = sums("cumulative_gross_divisions"),
    hypoxia_deaths = sums("cumulative_hypoxia_deaths"),
    buffer_inflow = sums("cumulative_dead_buffer_inflow"),
    nonlive_inflow = sums("cumulative_nonlive_inflow")
  )
}

.ivt_scan_safe_ratio <- function(numerator, denominator) {
  if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) return(NA_real_)
  numerator / denominator
}

.ivt_scan_blank_summary <- function(candidate_id, stage, run_params) {
  scanned <- c(
    "p_mis_base", "p_misseg", "k_o_mis", "gamma_mu", "buffer_smax",
    "buffer_beta", "buffer_n_exp", "O2_crit", "n_O"
  )
  fixed <- c("p_wgd", "gamma_growth", "mu_hp", "sigma_growth")
  values <- c(run_params[scanned], run_params[fixed])
  metrics <- list(
    status = NA_character_, error = NA_character_,
    failed_cohort = NA_character_, failed_scenario_id = NA_character_,
    failed_segment_id = NA_character_, failed_segment_ordinal = NA_real_,
    failed_threshold_target_cells = NA_real_, failed_max_live_cells = NA_real_,
    failed_relative_threshold_shortfall = NA_real_, passage_failure_reason = NA_character_,
    objective = NA_real_, growth_nll = NA_real_, ploidy_nll = NA_real_,
    flow_nll = NA_real_, death_nll = NA_real_, passage_time_nll = NA_real_,
    protocol_feasibility_status = NA_character_, n_scenarios = NA_real_,
    n_summary_rows = NA_real_, n_growth = NA_real_, n_death = NA_real_,
    n_flow_samples = NA_real_, n_ploidy_samples = NA_real_,
    n_insufficient_boundaries = NA_real_, n_growth_missing_pred = NA_real_,
    n_growth_negative_pred = NA_real_, n_missing_predictions = NA_real_,
    n_negative_predictions = NA_real_, all_passage_boundaries_feasible = NA,
    all_passage_executed = NA, all_selected_at_or_after_last_observation = NA,
    all_selected_within_tolerance = NA, n_invalid_distributions = NA_real_,
    control_flow_mass_ge3 = NA_real_, control_flow_observed_mass_ge3 = NA_real_,
    control_flow_mass_ge3_residual = NA_real_,
    pooled_2N_A19_mean_N = NA_real_, pooled_2N_A19_mass_ge70 = NA_real_,
    pooled_2N_A19_smoothed_mean_N = NA_real_,
    pooled_2N_A19_smoothed_mass_ge70 = NA_real_,
    abs_error_2N_A19_mean_N = NA_real_, abs_error_2N_A19_mass_ge70 = NA_real_,
    pooled_4N_A17_mean_N = NA_real_, pooled_4N_A17_smoothed_mean_N = NA_real_,
    abs_error_4N_A17_mean_N = NA_real_,
    flow_2N_low_A12_mass_ge3 = NA_real_, flow_2N_low_A18_mass_ge3 = NA_real_,
    flow_2N_low_A23_mass_ge3 = NA_real_, flow_2N_low_A12_wgd_mass = NA_real_,
    flow_2N_low_A18_wgd_mass = NA_real_, flow_2N_low_A23_wgd_mass = NA_real_,
    observed_flow_2N_low_A12_mass_ge3 = NA_real_, observed_flow_2N_low_A18_mass_ge3 = NA_real_,
    observed_flow_2N_low_A23_mass_ge3 = NA_real_, observed_flow_2N_low_A12_wgd_mass = NA_real_,
    observed_flow_2N_low_A18_wgd_mass = NA_real_, observed_flow_2N_low_A23_wgd_mass = NA_real_,
    flow_2N_low_A12_bimodal_count = NA_real_, flow_2N_low_A18_bimodal_count = NA_real_,
    flow_2N_low_A23_bimodal_count = NA_real_,
    observed_flow_2N_low_A12_bimodal_count = NA_real_,
    observed_flow_2N_low_A18_bimodal_count = NA_real_,
    observed_flow_2N_low_A23_bimodal_count = NA_real_,
    flow_low_o2_mass_ge3_rmse = NA_real_, flow_low_o2_wgd_mass_rmse = NA_real_,
    flow_low_o2_bimodal_count_abs_error = NA_real_,
    growth_log_residual_mean = NA_real_, growth_log_residual_rmse = NA_real_,
    growth_median_predicted_observed_ratio = NA_real_,
    growth_2N_log_rmse = NA_real_, growth_4N_log_rmse = NA_real_,
    growth_low_o2_log_rmse = NA_real_, growth_control_log_rmse = NA_real_,
    death_mean_observed_fraction = NA_real_, death_mean_predicted_fraction = NA_real_,
    death_fraction_bias_pred_minus_obs = NA_real_, death_logit_residual_mean = NA_real_,
    death_logit_residual_rmse = NA_real_, death_2N_logit_rmse = NA_real_,
    death_4N_logit_rmse = NA_real_,
    max_abs_selected_time_residual = NA_real_, n_selected_outside_tolerance = NA_real_,
    median_threshold_time_residual = NA_real_, q05_threshold_time_residual = NA_real_,
    min_threshold_time_residual = NA_real_, fraction_threshold_earlier_than_1_day = NA_real_,
    gross_divisions_2N = NA_real_, hypoxia_deaths_2N = NA_real_,
    dead_buffer_inflow_2N = NA_real_, nonlive_inflow_2N = NA_real_,
    hypoxia_death_to_division_2N = NA_real_, nonlive_to_division_2N = NA_real_,
    gross_divisions_4N = NA_real_, hypoxia_deaths_4N = NA_real_,
    dead_buffer_inflow_4N = NA_real_, nonlive_inflow_4N = NA_real_,
    hypoxia_death_to_division_4N = NA_real_, nonlive_to_division_4N = NA_real_,
    hard_pass = NA
  )
  c(
    list(candidate_id = candidate_id, stage = stage),
    as.list(values),
    metrics
  )
}

.ivt_scan_invalid_distribution_count <- function(run) {
  sum(vapply(run$segment_results, function(x) {
    p <- suppressWarnings(as.numeric(x$selection$selected_frac))
    !length(p) || any(!is.finite(p)) || any(p < 0) || !is.finite(sum(p)) || sum(p) <= 0
  }, logical(1)))
}

.ivt_scan_evaluate_task <- function(task) {
  run_params <- .ivt_scan_context$base_params
  for (name in names(task$params)) run_params[[name]] <- as.numeric(task$params[[name]])
  summary_values <- .ivt_scan_blank_summary(task$candidate_id, task$stage, run_params)

  tryCatch({
    comp <- do.call(ivt_objective_components, .ivt_scan_objective_args(run_params))
    flow <- .ivt_scan_collect_flow_metrics(comp$flow_overlay_df, task$candidate_id, task$stage)
    kary <- .ivt_scan_collect_kary_metrics(
      comp,
      task$candidate_id,
      task$stage,
      sigma_kary = as.numeric(run_params$sigma_kary)
    )

    control_pred <- flow$high_mass_ge3[
      flow$segment_id == "2N-C-A12" & flow$series == "Predicted"
    ]
    control_obs <- flow$high_mass_ge3[
      flow$segment_id == "2N-C-A12" & flow$series == "Observed"
    ]
    control_pred <- if (length(control_pred)) mean(control_pred) else NA_real_
    control_obs <- if (length(control_obs)) mean(control_obs) else NA_real_

    growth <- comp$growth_df
    growth_residual <- suppressWarnings(as.numeric(growth$log_live_cell_residual))
    if (!length(growth_residual) || all(!is.finite(growth_residual))) {
      growth_residual <- log(
        as.numeric(growth$predicted_live_cells_at_observation) /
          as.numeric(growth$observed_live_cells_at_observation)
      )
    }
    growth_stats <- .ivt_scan_metric_summary(growth_residual)
    growth_rmse_group <- function(keep) {
      x <- growth_residual[keep]
      x <- x[is.finite(x)]
      if (!length(x)) NA_real_ else sqrt(mean(x^2))
    }
    predicted_cells <- suppressWarnings(as.numeric(growth$predicted_live_cells_at_observation))
    observed_cells <- suppressWarnings(as.numeric(growth$observed_live_cells_at_observation))
    ratios <- predicted_cells / observed_cells

    death <- comp$death_df
    death_residual <- suppressWarnings(as.numeric(death$logit_residual))
    death_stats <- .ivt_scan_metric_summary(death_residual)
    death_rmse_group <- function(cohort) {
      x <- death_residual[as.character(death$cohort) == cohort]
      x <- x[is.finite(x)]
      if (!length(x)) NA_real_ else sqrt(mean(x^2))
    }

    passage <- comp$passage_time_df
    selected_residual <- suppressWarnings(as.numeric(passage$passage_time_residual_days))
    threshold_residual <- suppressWarnings(as.numeric(comp$summary$threshold_time_residual_days))
    threshold_stats <- .ivt_scan_metric_summary(threshold_residual)

    invalid_distributions <- .ivt_scan_invalid_distribution_count(comp$run_2N) +
      .ivt_scan_invalid_distribution_count(comp$run_4N)
    invalid_growth <- sum(
      !is.finite(predicted_cells) | predicted_cells <= 0 |
        !is.finite(observed_cells) | observed_cells <= 0
    )
    predicted_death <- suppressWarnings(as.numeric(death$predicted_dead_fraction))
    invalid_death <- sum(!is.finite(predicted_death) | predicted_death < 0 | predicted_death > 1)
    selected_day <- suppressWarnings(as.numeric(comp$summary$selected_day))
    last_observation_day <- suppressWarnings(as.numeric(comp$summary$last_observation_day))
    invalid_passage <- sum(!is.finite(selected_day) | !is.finite(last_observation_day))
    n_missing_predictions <- invalid_growth + invalid_death + invalid_passage + invalid_distributions
    n_negative_predictions <- as.numeric(comp$n_growth_negative_pred) +
      sum(is.finite(predicted_cells) & predicted_cells <= 0)

    event_2n <- .ivt_scan_terminal_event_metrics(comp$summary, "2N")
    event_4n <- .ivt_scan_terminal_event_metrics(comp$summary, "4N")
    target <- .ivt_scan_context$targets
    pooled_2n_mean <- .ivt_scan_pooled_metric(kary, "2N", "mean_N", "latent")
    pooled_2n_mass <- .ivt_scan_pooled_metric(kary, "2N", "mass_N_ge_70", "latent")
    pooled_2n_smooth_mean <- .ivt_scan_pooled_metric(kary, "2N", "mean_N", "observation_smoothed")
    pooled_2n_smooth_mass <- .ivt_scan_pooled_metric(kary, "2N", "mass_N_ge_70", "observation_smoothed")
    pooled_4n_mean <- .ivt_scan_pooled_metric(kary, "4N", "mean_N", "latent")
    pooled_4n_smooth_mean <- .ivt_scan_pooled_metric(kary, "4N", "mean_N", "observation_smoothed")
    low_passages <- c(12, 18, 23)
    predicted_high_mass <- vapply(
      low_passages,
      function(passage_id) .ivt_scan_flow_passage_mean(flow, passage_id, "high_mass_ge3", "Predicted"),
      numeric(1)
    )
    observed_high_mass <- vapply(
      low_passages,
      function(passage_id) .ivt_scan_flow_passage_mean(flow, passage_id, "high_mass_ge3", "Observed"),
      numeric(1)
    )
    predicted_wgd_mass <- vapply(
      low_passages,
      function(passage_id) .ivt_scan_flow_passage_mean(flow, passage_id, "wgd_mass", "Predicted"),
      numeric(1)
    )
    observed_wgd_mass <- vapply(
      low_passages,
      function(passage_id) .ivt_scan_flow_passage_mean(flow, passage_id, "wgd_mass", "Observed"),
      numeric(1)
    )
    bimodal_count <- function(passage_id, series) {
      keep <- as.character(flow$series) == series &
        grepl(paste0("^2N-O[12]-A", passage_id, "$"), as.character(flow$segment_id))
      sum(as.logical(flow$has_target_bimodality[keep]), na.rm = TRUE)
    }
    predicted_bimodal <- vapply(low_passages, bimodal_count, numeric(1), series = "Predicted")
    observed_bimodal <- vapply(low_passages, bimodal_count, numeric(1), series = "Observed")
    finite_rmse <- function(predicted, observed) {
      keep <- is.finite(predicted) & is.finite(observed)
      if (!any(keep)) NA_real_ else sqrt(mean((predicted[keep] - observed[keep])^2))
    }

    values <- list(
      status = "OK",
      error = "",
      objective = as.numeric(comp$objective),
      growth_nll = -as.numeric(comp$growth_loglik),
      ploidy_nll = -as.numeric(comp$ploidy_loglik),
      flow_nll = -as.numeric(comp$flow_loglik),
      death_nll = -as.numeric(comp$death_loglik),
      passage_time_nll = -as.numeric(comp$passage_time_loglik),
      protocol_feasibility_status = as.character(comp$protocol_feasibility_status),
      n_scenarios = as.numeric(comp$n_scenarios),
      n_summary_rows = nrow(comp$summary),
      n_growth = nrow(growth),
      n_death = nrow(death),
      n_flow_samples = nrow(comp$flow_df),
      n_ploidy_samples = nrow(comp$ploidy_df),
      n_insufficient_boundaries = as.numeric(comp$n_insufficient_boundaries),
      n_growth_missing_pred = as.numeric(comp$n_growth_missing_pred),
      n_growth_negative_pred = as.numeric(comp$n_growth_negative_pred),
      n_missing_predictions = n_missing_predictions,
      n_negative_predictions = n_negative_predictions,
      all_passage_boundaries_feasible = isTRUE(comp$all_passage_boundaries_feasible),
      all_passage_executed = all(as.logical(comp$summary$passage_executed)),
      all_selected_at_or_after_last_observation = all(selected_day >= last_observation_day),
      all_selected_within_tolerance = all(as.logical(comp$summary$passage_time_within_tolerance)),
      n_invalid_distributions = invalid_distributions,
      control_flow_mass_ge3 = control_pred,
      control_flow_observed_mass_ge3 = control_obs,
      control_flow_mass_ge3_residual = control_pred - control_obs,
      pooled_2N_A19_mean_N = pooled_2n_mean,
      pooled_2N_A19_mass_ge70 = pooled_2n_mass,
      pooled_2N_A19_smoothed_mean_N = pooled_2n_smooth_mean,
      pooled_2N_A19_smoothed_mass_ge70 = pooled_2n_smooth_mass,
      abs_error_2N_A19_mean_N = abs(pooled_2n_mean - target$pooled_2N_A19_mean_N),
      abs_error_2N_A19_mass_ge70 = abs(pooled_2n_mass - target$pooled_2N_A19_mass_ge70),
      pooled_4N_A17_mean_N = pooled_4n_mean,
      pooled_4N_A17_smoothed_mean_N = pooled_4n_smooth_mean,
      abs_error_4N_A17_mean_N = abs(pooled_4n_mean - target$pooled_4N_A17_mean_N),
      flow_2N_low_A12_mass_ge3 = predicted_high_mass[[1L]],
      flow_2N_low_A18_mass_ge3 = predicted_high_mass[[2L]],
      flow_2N_low_A23_mass_ge3 = predicted_high_mass[[3L]],
      flow_2N_low_A12_wgd_mass = predicted_wgd_mass[[1L]],
      flow_2N_low_A18_wgd_mass = predicted_wgd_mass[[2L]],
      flow_2N_low_A23_wgd_mass = predicted_wgd_mass[[3L]],
      observed_flow_2N_low_A12_mass_ge3 = observed_high_mass[[1L]],
      observed_flow_2N_low_A18_mass_ge3 = observed_high_mass[[2L]],
      observed_flow_2N_low_A23_mass_ge3 = observed_high_mass[[3L]],
      observed_flow_2N_low_A12_wgd_mass = observed_wgd_mass[[1L]],
      observed_flow_2N_low_A18_wgd_mass = observed_wgd_mass[[2L]],
      observed_flow_2N_low_A23_wgd_mass = observed_wgd_mass[[3L]],
      flow_2N_low_A12_bimodal_count = predicted_bimodal[[1L]],
      flow_2N_low_A18_bimodal_count = predicted_bimodal[[2L]],
      flow_2N_low_A23_bimodal_count = predicted_bimodal[[3L]],
      observed_flow_2N_low_A12_bimodal_count = observed_bimodal[[1L]],
      observed_flow_2N_low_A18_bimodal_count = observed_bimodal[[2L]],
      observed_flow_2N_low_A23_bimodal_count = observed_bimodal[[3L]],
      flow_low_o2_mass_ge3_rmse = finite_rmse(predicted_high_mass, observed_high_mass),
      flow_low_o2_wgd_mass_rmse = finite_rmse(predicted_wgd_mass, observed_wgd_mass),
      flow_low_o2_bimodal_count_abs_error = sum(abs(predicted_bimodal - observed_bimodal)),
      growth_log_residual_mean = growth_stats[["mean"]],
      growth_log_residual_rmse = growth_stats[["rmse"]],
      growth_median_predicted_observed_ratio = stats::median(ratios[is.finite(ratios)]),
      growth_2N_log_rmse = growth_rmse_group(as.character(growth$cohort) == "2N"),
      growth_4N_log_rmse = growth_rmse_group(as.character(growth$cohort) == "4N"),
      growth_low_o2_log_rmse = growth_rmse_group(as.character(growth$lineage_id) %in% c("O1", "O2")),
      growth_control_log_rmse = growth_rmse_group(as.character(growth$lineage_id) == "C"),
      death_mean_observed_fraction = mean(as.numeric(death$observed_dead_fraction)),
      death_mean_predicted_fraction = mean(predicted_death),
      death_fraction_bias_pred_minus_obs = mean(predicted_death - as.numeric(death$observed_dead_fraction)),
      death_logit_residual_mean = death_stats[["mean"]],
      death_logit_residual_rmse = death_stats[["rmse"]],
      death_2N_logit_rmse = death_rmse_group("2N"),
      death_4N_logit_rmse = death_rmse_group("4N"),
      max_abs_selected_time_residual = max(abs(selected_residual[is.finite(selected_residual)])),
      n_selected_outside_tolerance = sum(!as.logical(passage$passage_time_within_tolerance)),
      median_threshold_time_residual = threshold_stats[["median"]],
      q05_threshold_time_residual = threshold_stats[["q05"]],
      min_threshold_time_residual = threshold_stats[["min"]],
      fraction_threshold_earlier_than_1_day = mean(threshold_residual < -1, na.rm = TRUE),
      gross_divisions_2N = event_2n[["gross_divisions"]],
      hypoxia_deaths_2N = event_2n[["hypoxia_deaths"]],
      dead_buffer_inflow_2N = event_2n[["buffer_inflow"]],
      nonlive_inflow_2N = event_2n[["nonlive_inflow"]],
      hypoxia_death_to_division_2N = .ivt_scan_safe_ratio(event_2n[["hypoxia_deaths"]], event_2n[["gross_divisions"]]),
      nonlive_to_division_2N = .ivt_scan_safe_ratio(event_2n[["nonlive_inflow"]], event_2n[["gross_divisions"]]),
      gross_divisions_4N = event_4n[["gross_divisions"]],
      hypoxia_deaths_4N = event_4n[["hypoxia_deaths"]],
      dead_buffer_inflow_4N = event_4n[["buffer_inflow"]],
      nonlive_inflow_4N = event_4n[["nonlive_inflow"]],
      hypoxia_death_to_division_4N = .ivt_scan_safe_ratio(event_4n[["hypoxia_deaths"]], event_4n[["gross_divisions"]]),
      nonlive_to_division_4N = .ivt_scan_safe_ratio(event_4n[["nonlive_inflow"]], event_4n[["gross_divisions"]])
    )
    for (name in names(values)) summary_values[[name]] <- values[[name]]
    summary_df <- as.data.frame(summary_values, stringsAsFactors = FALSE, check.names = FALSE)
    summary_df$hard_pass <- ivt_scan_hard_gate(summary_df, .ivt_scan_context$control_limit)
    list(summary = summary_df, flow = flow, kary = kary)
  }, invitro_protocol_infeasible = function(e) {
    selection <- e$selection %||% list()
    segment <- e$segment %||% list()
    threshold <- suppressWarnings(as.numeric(selection$threshold_target_cells %||% NA_real_))
    max_live <- suppressWarnings(as.numeric(selection$max_live_cells_in_search %||% NA_real_))
    summary_values$status <- "PROTOCOL_FAIL"
    summary_values$error <- conditionMessage(e)
    summary_values$protocol_feasibility_status <- "FAIL"
    summary_values$failed_cohort <- as.character(segment$cohort %||% NA_character_)
    summary_values$failed_scenario_id <- as.character(segment$scenario_id %||% NA_character_)
    summary_values$failed_segment_id <- as.character(segment$segment_id %||% NA_character_)
    summary_values$failed_segment_ordinal <- as.numeric(e$segment_ordinal %||% NA_real_)
    summary_values$failed_threshold_target_cells <- threshold
    summary_values$failed_max_live_cells <- max_live
    summary_values$failed_relative_threshold_shortfall <- if (
      length(threshold) == 1L && is.finite(threshold) && threshold > 0 &&
        length(max_live) == 1L && is.finite(max_live)
    ) max((threshold - max_live) / threshold, 0) else NA_real_
    summary_values$passage_failure_reason <- as.character(selection$passage_failure_reason %||% NA_character_)
    summary_values$hard_pass <- FALSE
    list(
      summary = as.data.frame(summary_values, stringsAsFactors = FALSE, check.names = FALSE),
      flow = data.frame(),
      kary = data.frame()
    )
  }, error = function(e) {
    summary_values$status <- "SIMULATION_ERROR"
    summary_values$error <- conditionMessage(e)
    summary_values$hard_pass <- FALSE
    list(
      summary = as.data.frame(summary_values, stringsAsFactors = FALSE, check.names = FALSE),
      flow = data.frame(),
      kary = data.frame()
    )
  })
}

.ivt_scan_task <- function(candidate_id, stage, params) {
  list(candidate_id = candidate_id, stage = stage, params = as.list(params))
}

.ivt_scan_rbind <- function(rows, name) {
  pieces <- lapply(rows, `[[`, name)
  pieces <- pieces[vapply(pieces, function(x) is.data.frame(x) && nrow(x), logical(1))]
  if (length(pieces)) return(do.call(rbind, pieces))
  if (identical(name, "flow")) {
    return(data.frame(
      candidate_id = character(), stage = character(), segment_id = character(),
      passage_id = character(), sample_name = character(), series = character(),
      high_mass_ge3 = numeric(), mean_ploidy = numeric(), sd_ploidy = numeric(),
      bridge_mass = numeric(), wgd_mass = numeric(), mass_3_5 = numeric(),
      tail_gt4_5 = numeric(), tail_gt5 = numeric(), low_peak = numeric(),
      high_peak = numeric(), peak_ratio = numeric(), valley_ratio = numeric(),
      high_relative_height = numeric(), has_target_bimodality = logical(),
      stringsAsFactors = FALSE
    ))
  }
  if (identical(name, "kary")) {
    return(data.frame(
      candidate_id = character(), stage = character(), segment_id = character(),
      cohort = character(), lineage_id = character(), scenario_id = character(),
      oxygen_pct = numeric(), distribution_basis = character(), sigma_kary = numeric(),
      mean_N = numeric(), sd_N = numeric(), median_N = numeric(),
      mass_N_ge_70 = numeric(), mass_N_ge_80 = numeric(), mass_N_ge_88 = numeric(),
      n_grid_states = integer(), stringsAsFactors = FALSE
    ))
  }
  data.frame()
}

.ivt_scan_worker_init <- function(cluster, cpp_info) {
  parallel::clusterCall(
    cluster,
    function(wrapper_path, dll_path) {
      Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
      tryCatch({
        if (!file.exists(wrapper_path)) stop("Missing C++ wrapper: ", wrapper_path)
        if (!file.exists(dll_path)) stop("Missing C++ DLL: ", dll_path)
        source(wrapper_path, local = .GlobalEnv)
        required <- c(
          "cpp_o2simps_build_G_for_o2_triplet",
          "cpp_o2simps_simulate_one",
          "cpp_o2simps_objective_components_map"
        )
        missing <- required[!vapply(required, exists, logical(1), envir = .GlobalEnv, mode = "function", inherits = TRUE)]
        if (length(missing)) stop("Wrapper missing function(s): ", paste(missing, collapse = ", "))
        list(ok = TRUE, status = "PASS_WRAPPER_DLL_LOADED", error = NA_character_)
      }, error = function(e) {
        list(ok = FALSE, status = "WORKER_CPP_INIT_ERROR", error = conditionMessage(e))
      })
    },
    cpp_info$wrapper_path,
    cpp_info$path
  )
}

.ivt_scan_start_cluster <- function(n_cores,
                                    probe_task,
                                    audit_path,
                                    max_retries = 5L) {
  n_cores <- .ivt_scan_as_integer(n_cores, "n_cores")
  max_retries <- .ivt_scan_as_integer(max_retries, "max_preflight_retries", minimum = 0L)
  master <- .ivt_scan_evaluate_task(probe_task)
  master_objective <- suppressWarnings(as.numeric(master$summary$objective[[1L]]))
  if (!identical(as.character(master$summary$status[[1L]]), "OK") || !is.finite(master_objective)) {
    stop("Master preflight candidate failed: ", master$summary$error[[1L]], call. = FALSE)
  }
  if (n_cores == 1L) {
    audit <- data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      attempt = 1L, retry_number = 0L, worker = 0L, role = "master",
      ok = TRUE, status = "PASS", objective = master_objective,
      error = NA_character_, stringsAsFactors = FALSE
    )
    .ivt_scan_write_tsv(audit, audit_path)
    return(list(cluster = NULL, master = master, audit = audit))
  }
  if (.Platform$OS.type != "unix" ||
      !exists("makeForkCluster", envir = asNamespace("parallel"), inherits = FALSE)) {
    stop("Parallel scan requires a Unix persistent fork cluster.", call. = FALSE)
  }

  cpp_info <- o2simps_cpp_dll_info()
  all_audit <- list()
  for (attempt in seq_len(max_retries + 1L)) {
    cluster <- NULL
    attempt_rows <- list()
    result <- tryCatch({
      cluster <- parallel::makeForkCluster(
        n_cores,
        outfile = file.path(dirname(audit_path), "worker_output.log")
      )
      init <- .ivt_scan_worker_init(cluster, cpp_info)
      init_rows <- do.call(rbind, lapply(seq_along(init), function(i) {
        x <- init[[i]]
        data.frame(
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
          attempt = attempt, retry_number = attempt - 1L, worker = i,
          role = "worker_init", ok = isTRUE(x$ok), status = as.character(x$status),
          objective = NA_real_, error = as.character(x$error), stringsAsFactors = FALSE
        )
      }))
      attempt_rows[[length(attempt_rows) + 1L]] <- init_rows
      if (!all(init_rows$ok)) stop("At least one worker failed explicit wrapper/DLL initialization.")

      probes <- parallel::clusterCall(cluster, function(task) {
        out <- .ivt_scan_evaluate_task(task)
        list(
          status = as.character(out$summary$status[[1L]]),
          objective = suppressWarnings(as.numeric(out$summary$objective[[1L]])),
          error = as.character(out$summary$error[[1L]])
        )
      }, probe_task)
      probe_rows <- do.call(rbind, lapply(seq_along(probes), function(i) {
        x <- probes[[i]]
        tolerance <- 1e-8 * max(1, abs(master_objective))
        ok <- identical(x$status, "OK") && is.finite(x$objective) &&
          abs(x$objective - master_objective) <= tolerance
        data.frame(
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
          attempt = attempt, retry_number = attempt - 1L, worker = i,
          role = "worker_probe", ok = ok,
          status = if (ok) "PASS_OBJECTIVE_MATCH" else "OBJECTIVE_MISMATCH",
          objective = x$objective, error = x$error, stringsAsFactors = FALSE
        )
      }))
      attempt_rows[[length(attempt_rows) + 1L]] <- probe_rows
      if (!all(probe_rows$ok)) stop("At least one worker failed the objective-consistency probe.")
      list(ok = TRUE, error = NA_character_)
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))

    if (!length(attempt_rows)) {
      attempt_rows[[1L]] <- data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        attempt = attempt, retry_number = attempt - 1L, worker = NA_integer_,
        role = "cluster", ok = FALSE, status = "CLUSTER_ERROR",
        objective = NA_real_, error = result$error, stringsAsFactors = FALSE
      )
    }
    attempt_df <- do.call(rbind, attempt_rows)
    if (!isTRUE(result$ok)) {
      attempt_df <- rbind(attempt_df, data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        attempt = attempt, retry_number = attempt - 1L, worker = NA_integer_,
        role = "attempt", ok = FALSE, status = "PREFLIGHT_FAILED",
        objective = NA_real_, error = result$error, stringsAsFactors = FALSE
      ))
    }
    all_audit[[length(all_audit) + 1L]] <- attempt_df
    .ivt_scan_write_tsv(do.call(rbind, all_audit), audit_path)
    if (isTRUE(result$ok)) {
      return(list(cluster = cluster, master = master, audit = do.call(rbind, all_audit)))
    }
    if (!is.null(cluster)) try(parallel::stopCluster(cluster), silent = TRUE)
    if (attempt <= max_retries) Sys.sleep(min(5, attempt))
  }
  stop(
    "Worker preflight failed after ", max_retries, " retries; see ", audit_path,
    call. = FALSE
  )
}

.ivt_scan_evaluate_tasks <- function(tasks, cluster) {
  if (is.null(cluster)) return(lapply(tasks, .ivt_scan_evaluate_task))
  parallel::parLapplyLB(cluster, tasks, function(task) .ivt_scan_evaluate_task(task))
}

.ivt_scan_calibration_tasks <- function(values, base_params, stage) {
  lapply(seq_along(values), function(i) {
    params <- c(
      p_mis_base = values[[i]],
      p_wgd = base_params$p_wgd,
      gamma_growth = base_params$gamma_growth,
      mu_hp = base_params$mu_hp,
      sigma_growth = base_params$sigma_growth
    )
    .ivt_scan_task(
      sprintf("calibration_%s_%03d", stage, i),
      paste0("calibration_", stage),
      params
    )
  })
}

.ivt_scan_observed_targets <- function(contract_comp, fit_objects) {
  ids_2n <- c("2N-O1-A19", "2N-O2-A19")
  ids_4n <- c("4N-O1-A17", "4N-O2-A17")
  obs_2n <- lapply(ids_2n, function(id) .ivt_scan_observed_kary(contract_comp$run_2N, fit_objects$fit_data, id))
  obs_4n <- lapply(ids_4n, function(id) .ivt_scan_observed_kary(contract_comp$run_4N, fit_objects$fit_data, id))
  pooled_2n <- unlist(obs_2n, use.names = FALSE)
  pooled_4n <- unlist(obs_4n, use.names = FALSE)
  flow <- .ivt_scan_collect_flow_metrics(contract_comp$flow_overlay_df, "observed", "observed")
  control_obs <- flow$high_mass_ge3[flow$segment_id == "2N-C-A12" & flow$series == "Observed"]
  target_list <- list(
    pooled_2N_A19_mean_N = mean(pooled_2n),
    pooled_2N_A19_mass_ge70 = mean(pooled_2n >= 70),
    pooled_4N_A17_mean_N = mean(pooled_4n),
    control_2N_C_A12_flow_mass_ge3 = mean(control_obs)
  )
  table <- data.frame(
    metric = c(
      "2N_O1_A19_kary_mean_N", "2N_O2_A19_kary_mean_N",
      "2N_O1_A19_kary_mass_N_ge70", "2N_O2_A19_kary_mass_N_ge70",
      "pooled_2N_A19_kary_mean_N", "pooled_2N_A19_kary_mass_N_ge70",
      "4N_O1_A17_kary_mean_N", "4N_O2_A17_kary_mean_N",
      "pooled_4N_A17_kary_mean_N", "2N_C_A12_observed_flow_mass_ploidy_ge3"
    ),
    value = c(
      mean(obs_2n[[1L]]), mean(obs_2n[[2L]]),
      mean(obs_2n[[1L]] >= 70), mean(obs_2n[[2L]] >= 70),
      target_list$pooled_2N_A19_mean_N, target_list$pooled_2N_A19_mass_ge70,
      mean(obs_4n[[1L]]), mean(obs_4n[[2L]]),
      target_list$pooled_4N_A17_mean_N, target_list$control_2N_C_A12_flow_mass_ge3
    ),
    n_observations = c(
      length(obs_2n[[1L]]), length(obs_2n[[2L]]),
      length(obs_2n[[1L]]), length(obs_2n[[2L]]),
      length(pooled_2n), length(pooled_2n),
      length(obs_4n[[1L]]), length(obs_4n[[2L]]),
      length(pooled_4n), length(control_obs)
    ),
    stringsAsFactors = FALSE
  )
  list(values = target_list, table = table)
}

.ivt_scan_contract_table <- function(comp, target_table) {
  observed <- c(
    n_scenarios = comp$n_scenarios,
    n_summary_rows = nrow(comp$summary),
    n_growth_observations = nrow(comp$growth_df),
    n_death_observations = nrow(comp$death_df),
    n_passage_time_rows = nrow(comp$passage_time_df),
    n_flow_samples = nrow(comp$flow_df),
    n_karyotype_samples = nrow(comp$ploidy_df),
    pooled_2N_A19_mean_N = target_table$value[target_table$metric == "pooled_2N_A19_kary_mean_N"],
    pooled_2N_A19_mass_N_ge70 = target_table$value[target_table$metric == "pooled_2N_A19_kary_mass_N_ge70"],
    pooled_4N_A17_mean_N = target_table$value[target_table$metric == "pooled_4N_A17_kary_mean_N"]
  )
  expected <- c(6, 114, 219, 90, 114, 20, 12, 77.45, 0.625, 81.475)
  tolerance <- c(rep(0, 7), rep(1e-8, 3))
  data.frame(
    check = names(observed),
    observed = as.numeric(observed),
    expected = expected,
    tolerance = tolerance,
    pass = abs(as.numeric(observed) - expected) <= tolerance,
    stringsAsFactors = FALSE
  )
}

.ivt_scan_finalize <- function(summary_df, flow_df, kary_df, out_dir) {
  summary_df$hard_pass <- ivt_scan_hard_gate(summary_df, .ivt_scan_context$control_limit)
  metric_cols <- c(
    "abs_error_2N_A19_mean_N", "abs_error_2N_A19_mass_ge70",
    "abs_error_4N_A17_mean_N", "flow_low_o2_mass_ge3_rmse",
    "flow_low_o2_wgd_mass_rmse", "flow_low_o2_bimodal_count_abs_error", "flow_nll",
    "growth_log_residual_rmse", "death_logit_residual_rmse"
  )
  summary_df$pareto_front <- FALSE
  pass_idx <- which(summary_df$hard_pass)
  if (length(pass_idx)) {
    summary_df$pareto_front[pass_idx] <- ivt_scan_pareto_front(summary_df[pass_idx, , drop = FALSE], metric_cols)
  }
  for (metric in metric_cols) {
    rank_name <- paste0("rank_", metric)
    summary_df[[rank_name]] <- NA_real_
    if (length(pass_idx)) {
      summary_df[[rank_name]][pass_idx] <- rank(
        as.numeric(summary_df[[metric]][pass_idx]),
        ties.method = "min",
        na.last = "keep"
      )
    }
  }
  order_idx <- order(
    !summary_df$hard_pass,
    !summary_df$pareto_front,
    summary_df$abs_error_2N_A19_mean_N,
    summary_df$abs_error_2N_A19_mass_ge70,
    summary_df$flow_nll,
    na.last = TRUE
  )
  summary_df <- summary_df[order_idx, , drop = FALSE]
  rownames(summary_df) <- NULL
  .ivt_scan_write_tsv(summary_df, file.path(out_dir, "feasibility_summary.tsv"))
  .ivt_scan_write_tsv(
    summary_df[summary_df$hard_pass, , drop = FALSE],
    file.path(out_dir, "feasibility_hard_pass.tsv")
  )
  .ivt_scan_write_tsv(
    summary_df[summary_df$hard_pass & summary_df$pareto_front, , drop = FALSE],
    file.path(out_dir, "feasibility_pareto.tsv")
  )
  failure_cols <- c(
    "candidate_id", "status", "error", "failed_cohort", "failed_scenario_id",
    "failed_segment_id", "failed_segment_ordinal", "failed_threshold_target_cells",
    "failed_max_live_cells", "failed_relative_threshold_shortfall", "passage_failure_reason"
  )
  .ivt_scan_write_tsv(
    summary_df[summary_df$status != "OK", failure_cols, drop = FALSE],
    file.path(out_dir, "feasibility_failures.tsv")
  )
  .ivt_scan_write_tsv(flow_df, file.path(out_dir, "feasibility_flow_metrics.tsv"))
  .ivt_scan_write_tsv(kary_df, file.path(out_dir, "feasibility_kary_metrics.tsv"))
  invisible(summary_df)
}

main <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_result_path <- normalizePath(
    .ivt_scan_first(argv$fit_result, stop("--fit_result is required.", call. = FALSE)),
    mustWork = TRUE
  )
  fit_dir <- dirname(fit_result_path)
  out_dir <- normalizePath(
    .ivt_scan_first(argv$out_dir, stop("--out_dir is required.", call. = FALSE)),
    mustWork = FALSE
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- file.path(out_dir, "checkpoints")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  n_candidates <- .ivt_scan_as_integer(.ivt_scan_first(argv$n_candidates, 5000L), "n_candidates")
  n_cores <- .ivt_scan_as_integer(.ivt_scan_first(argv$n_cores, 1L), "n_cores")
  design_seed <- .ivt_scan_as_integer(.ivt_scan_first(argv$design_seed, 10L), "design_seed", minimum = 0L)
  calibration_points <- .ivt_scan_as_integer(.ivt_scan_first(argv$calibration_points, 31L), "calibration_points", minimum = 3L)
  calibration_refine_points <- .ivt_scan_as_integer(.ivt_scan_first(argv$calibration_refine_points, 20L), "calibration_refine_points", minimum = 0L)
  batch_size <- .ivt_scan_as_integer(.ivt_scan_first(argv$batch_size, max(n_cores * 5L, n_cores)), "batch_size")
  control_limit <- .ivt_scan_as_number(.ivt_scan_first(argv$control_high_mass_limit, 0.05), "control_high_mass_limit", 0, 1)
  max_preflight_retries <- .ivt_scan_as_integer(.ivt_scan_first(argv$max_preflight_retries, 5L), "max_preflight_retries", minimum = 0L)
  resume <- .ivt_scan_as_bool(argv$resume, TRUE)

  fit_result <- readRDS(fit_result_path)
  if (is.null(fit_result$cfg) || is.null(fit_result$best_params)) {
    stop("fit_result must contain cfg and best_params.", call. = FALSE)
  }
  parameter_table_path <- .ivt_scan_resolve_input(
    argv$parameter_table,
    file.path(fit_dir, "parameter_table_input.csv"),
    file.path(OXYGEN_ROOT, "results", "fit_invitro_O2_buffering_500seed", "seed10", "parameter_table_input.csv"),
    "parameter table"
  )
  fit_objects_dir <- .ivt_scan_resolve_input(
    argv$fit_objects_dir,
    fit_result$fit_objects_dir,
    file.path(OXYGEN_ROOT, "ploidyOxygen", "data", "fit_objects"),
    "fit objects directory",
    directory = TRUE
  )
  flow_density_path <- .ivt_scan_resolve_input(
    argv$flow_density_path,
    fit_result$flow_density_path,
    file.path(OXYGEN_ROOT, "data", "g0g1_ploidy_density_grid.csv"),
    "flow-density data"
  )
  death_data_path <- .ivt_scan_resolve_input(
    argv$death_data_path,
    fit_result$death_data_path,
    file.path(PROJECT_ROOT, "data", "InVitroData", "sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"),
    "death-likelihood data"
  )
  parameter_table <- utils::read.csv(parameter_table_path, stringsAsFactors = FALSE, check.names = FALSE)
  p_row <- parameter_table[as.character(parameter_table$param_symbol) == "p_mis_base", , drop = FALSE]
  if (nrow(p_row) != 1L) stop("Parameter table must contain one p_mis_base row.")
  p_lower <- as.numeric(p_row$lower_bound[[1L]])
  p_upper <- as.numeric(p_row$upper_bound[[1L]])
  if (!is.finite(p_lower) || !is.finite(p_upper) || p_lower <= 0 || p_lower > p_upper) {
    stop("Invalid p_mis_base bounds in parameter table.")
  }

  fit_objects <- ivt_load_fit_objects(
    repo_root = OXYGEN_ROOT,
    fit_objects_dir = fit_objects_dir,
    flow_csv_path = flow_density_path,
    death_data_path = death_data_path
  )
  base_params <- fit_result$best_params
  base_params$p_wgd <- .ivt_scan_as_number(.ivt_scan_first(argv$fixed_p_wgd, 0.00165), "fixed_p_wgd", 0, 1)
  base_params$gamma_growth <- .ivt_scan_as_number(.ivt_scan_first(argv$fixed_gamma_growth, 0.24), "fixed_gamma_growth", 0, Inf)
  base_params$mu_hp <- .ivt_scan_as_number(.ivt_scan_first(argv$fixed_mu_hp, 0.00045), "fixed_mu_hp", 0, Inf)
  base_params$sigma_growth <- .ivt_scan_as_number(.ivt_scan_first(argv$fixed_sigma_growth, 0.032), "fixed_sigma_growth", 0, Inf)

  .ivt_scan_context <<- list(
    fit_objects = fit_objects,
    cfg = fit_result$cfg,
    base_params = base_params,
    targets = list(
      pooled_2N_A19_mean_N = 77.45,
      pooled_2N_A19_mass_ge70 = 0.625,
      pooled_4N_A17_mean_N = 81.475
    ),
    control_limit = control_limit,
    ploidy_weight = as.numeric(fit_result$ploidy_weight %||% 1),
    death_weight = as.numeric(fit_result$death_weight %||% 1),
    passage_time_weight = as.numeric(fit_result$passage_time_weight %||% 0.25),
    passage_time_tolerance_days = as.numeric(fit_result$passage_time_tolerance_days %||% 1),
    passage_time_sigma_days = as.numeric(fit_result$passage_time_sigma_days %||% 1),
    passage_time_df = as.numeric(fit_result$passage_time_df %||% 4),
    sigma_death_logit = as.numeric(fit_result$sigma_death_logit %||% 0.75),
    death_fraction_eps = as.numeric(fit_result$death_fraction_eps %||% 1e-4)
  )

  run_contract <- .ivt_scan_build_run_contract(
    fit_result_path = fit_result_path,
    parameter_table_path = parameter_table_path,
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path,
    death_data_path = death_data_path,
    settings = list(
      n_candidates = n_candidates,
      n_cores = n_cores,
      design_seed = design_seed,
      calibration_points = calibration_points,
      calibration_refine_points = calibration_refine_points,
      batch_size = batch_size,
      control_high_mass_limit = control_limit,
      p_mis_base_source_lower = p_lower,
      p_mis_base_source_upper = p_upper,
      fixed_p_wgd = base_params$p_wgd,
      fixed_gamma_growth = base_params$gamma_growth,
      fixed_mu_hp = base_params$mu_hp,
      fixed_sigma_growth = base_params$sigma_growth,
      ploidy_weight = .ivt_scan_context$ploidy_weight,
      death_weight = .ivt_scan_context$death_weight,
      passage_time_weight = .ivt_scan_context$passage_time_weight,
      passage_time_tolerance_days = .ivt_scan_context$passage_time_tolerance_days,
      passage_time_sigma_days = .ivt_scan_context$passage_time_sigma_days,
      passage_time_df = .ivt_scan_context$passage_time_df,
      sigma_death_logit = .ivt_scan_context$sigma_death_logit,
      death_fraction_eps = .ivt_scan_context$death_fraction_eps,
      container_image = Sys.getenv("O2SD_CONTAINER_IMAGE", unset = "")
    )
  )
  .ivt_scan_prepare_run_contract(run_contract, out_dir, resume)

  probe_params <- c(
    p_mis_base = as.numeric(base_params$p_mis_base),
    p_wgd = as.numeric(base_params$p_wgd),
    gamma_growth = as.numeric(base_params$gamma_growth),
    mu_hp = as.numeric(base_params$mu_hp),
    sigma_growth = as.numeric(base_params$sigma_growth)
  )
  probe_task <- .ivt_scan_task("preflight_anchor", "preflight", probe_params)
  cluster_state <- .ivt_scan_start_cluster(
    n_cores = n_cores,
    probe_task = probe_task,
    audit_path = file.path(out_dir, "worker_preflight.tsv"),
    max_retries = max_preflight_retries
  )
  cluster <- cluster_state$cluster
  if (!is.null(cluster)) on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
  contract_result <- cluster_state$master
  contract_comp <- do.call(
    ivt_objective_components,
    .ivt_scan_objective_args(base_params)
  )
  observed_targets <- .ivt_scan_observed_targets(contract_comp, fit_objects)
  .ivt_scan_context$targets <- observed_targets$values
  .ivt_scan_write_tsv(observed_targets$table, file.path(out_dir, "observed_targets.tsv"))
  contract <- .ivt_scan_contract_table(contract_comp, observed_targets$table)
  .ivt_scan_write_tsv(contract, file.path(out_dir, "acceptance_contract.tsv"))
  if (!all(contract$pass)) {
    stop(
      "Observed-data contract failed: ",
      paste(contract$check[!contract$pass], collapse = ", "),
      call. = FALSE
    )
  }

  calibration_path <- file.path(out_dir, "p_mis_base_calibration.tsv")
  calibration_summary_path <- file.path(out_dir, "p_mis_base_calibration_summary.tsv")
  calibration_rds_path <- file.path(out_dir, "p_mis_base_calibration.rds")
  if (resume && file.exists(calibration_rds_path)) {
    calibration_state <- readRDS(calibration_rds_path)
    calibration <- calibration_state$calibration
    cap <- calibration_state$cap
  } else if (resume && (file.exists(calibration_path) || file.exists(calibration_summary_path))) {
    stop(
      "Calibration TSV exists without the exact calibration RDS; refusing an imprecise resume.",
      call. = FALSE
    )
  } else {
    coarse_values <- exp(seq(log(p_lower), log(p_upper), length.out = calibration_points))
    coarse_rows <- .ivt_scan_evaluate_tasks(
      .ivt_scan_calibration_tasks(coarse_values, base_params, "coarse"),
      cluster
    )
    calibration <- .ivt_scan_rbind(coarse_rows, "summary")
    cap <- ivt_scan_select_calibration_cap(calibration, control_limit)
    if (!is.finite(cap$p_mis_base_upper[[1L]])) {
      .ivt_scan_write_tsv(calibration, calibration_path)
      .ivt_scan_write_tsv(cap, calibration_summary_path)
      saveRDS(list(calibration = calibration, cap = cap), calibration_rds_path)
      stop("No p_mis_base value passed the normoxic control calibration.", call. = FALSE)
    }
    if (!isTRUE(cap$all_pass[[1L]]) && calibration_refine_points > 0L) {
      lower_bracket <- cap$p_mis_base_upper[[1L]]
      upper_bracket <- cap$first_failing_p_mis_base[[1L]]
      refine_values <- exp(seq(log(lower_bracket), log(upper_bracket), length.out = calibration_refine_points + 2L))
      refine_values <- refine_values[-c(1L, length(refine_values))]
      refine_rows <- .ivt_scan_evaluate_tasks(
        .ivt_scan_calibration_tasks(refine_values, base_params, "refine"),
        cluster
      )
      calibration <- rbind(calibration, .ivt_scan_rbind(refine_rows, "summary"))
      calibration <- calibration[order(as.numeric(calibration$p_mis_base)), , drop = FALSE]
      cap <- ivt_scan_select_calibration_cap(calibration, control_limit)
    }
    .ivt_scan_write_tsv(calibration, calibration_path)
    .ivt_scan_write_tsv(cap, calibration_summary_path)
    saveRDS(list(calibration = calibration, cap = cap), calibration_rds_path)
  }
  calibrated_upper <- suppressWarnings(as.numeric(cap$p_mis_base_upper[[1L]]))
  if (!is.finite(calibrated_upper) || calibrated_upper < p_lower || calibrated_upper > p_upper) {
    stop("Calibrated p_mis_base upper bound is invalid.", call. = FALSE)
  }

  specs <- ivt_scan_parameter_specs(parameter_table, calibrated_upper)
  design <- ivt_scan_lhs_design(n_candidates, specs, seed = design_seed)
  design$p_wgd <- as.numeric(base_params$p_wgd)
  design$gamma_growth <- as.numeric(base_params$gamma_growth)
  design$mu_hp <- as.numeric(base_params$mu_hp)
  design$sigma_growth <- as.numeric(base_params$sigma_growth)
  .ivt_scan_write_tsv(specs, file.path(out_dir, "feasibility_parameter_specs.tsv"))
  design_rds <- file.path(out_dir, "feasibility_design.rds")
  if (resume && file.exists(design_rds)) {
    prior_design <- readRDS(design_rds)
    if (!identical(prior_design, design)) {
      stop("Existing feasibility design does not match current arguments.", call. = FALSE)
    }
  } else {
    saveRDS(design, design_rds)
  }
  .ivt_scan_write_tsv(design, file.path(out_dir, "feasibility_design.tsv"))

  parameter_names <- c(specs$param_symbol, "p_wgd", "gamma_growth", "mu_hp", "sigma_growth")
  batches <- split(seq_len(nrow(design)), ceiling(seq_len(nrow(design)) / batch_size))
  for (batch_id in seq_along(batches)) {
    checkpoint <- file.path(checkpoint_dir, sprintf("batch_%04d.rds", batch_id))
    if (resume && file.exists(checkpoint)) {
      message("[observed_feasibility] Resume: checkpoint exists for batch ", batch_id, "/", length(batches))
      next
    }
    idx <- batches[[batch_id]]
    tasks <- lapply(idx, function(i) {
      .ivt_scan_task(
        as.character(design$candidate_id[[i]]),
        "feasibility",
        unlist(design[i, parameter_names, drop = FALSE], use.names = TRUE)
      )
    })
    started <- Sys.time()
    rows <- .ivt_scan_evaluate_tasks(tasks, cluster)
    checkpoint_object <- list(
      batch_id = batch_id,
      candidate_ids = as.character(design$candidate_id[idx]),
      summary = .ivt_scan_rbind(rows, "summary"),
      flow = .ivt_scan_rbind(rows, "flow"),
      kary = .ivt_scan_rbind(rows, "kary"),
      started_at = format(started, "%Y-%m-%d %H:%M:%S %Z"),
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
    tmp_checkpoint <- paste0(checkpoint, ".tmp")
    saveRDS(checkpoint_object, tmp_checkpoint)
    if (!file.rename(tmp_checkpoint, checkpoint)) {
      stop("Failed to atomically install checkpoint: ", checkpoint)
    }
    message(
      "[observed_feasibility] Completed batch ", batch_id, "/", length(batches),
      " (", length(idx), " candidates)"
    )
  }

  checkpoint_paths <- file.path(checkpoint_dir, sprintf("batch_%04d.rds", seq_along(batches)))
  missing_checkpoints <- checkpoint_paths[!file.exists(checkpoint_paths)]
  if (length(missing_checkpoints)) {
    stop("Scan incomplete; missing checkpoint(s): ", paste(basename(missing_checkpoints), collapse = ", "))
  }
  checkpoint_objects <- lapply(checkpoint_paths, readRDS)
  checkpoint_membership_ok <- vapply(seq_along(checkpoint_objects), function(i) {
    expected_ids <- as.character(design$candidate_id[batches[[i]]])
    object <- checkpoint_objects[[i]]
    identical(as.integer(object$batch_id), as.integer(i)) &&
      identical(as.character(object$candidate_ids), expected_ids) &&
      identical(as.character(object$summary$candidate_id), expected_ids)
  }, logical(1))
  if (!all(checkpoint_membership_ok)) {
    stop(
      "Checkpoint membership does not match the current design for batch(es): ",
      paste(which(!checkpoint_membership_ok), collapse = ", "),
      call. = FALSE
    )
  }
  summary_df <- do.call(rbind, lapply(checkpoint_objects, `[[`, "summary"))
  flow_df <- do.call(rbind, lapply(checkpoint_objects, `[[`, "flow"))
  kary_df <- do.call(rbind, lapply(checkpoint_objects, `[[`, "kary"))
  if (nrow(summary_df) != n_candidates || anyDuplicated(summary_df$candidate_id) ||
      !identical(sort(as.character(summary_df$candidate_id)), sort(as.character(design$candidate_id)))) {
    stop("Final summary does not contain exactly one row per requested candidate.")
  }
  final_summary <- .ivt_scan_finalize(summary_df, flow_df, kary_df, out_dir)

  checksum_paths <- c(
    normalizePath(file.path(.script_dir, "scan_invitro_observed_feasibility.R"), mustWork = TRUE),
    normalizePath(file.path(.script_dir, "invitro_observed_feasibility_utils.R"), mustWork = TRUE),
    fit_result_path,
    parameter_table_path,
    flow_density_path,
    death_data_path
  )
  input_checksums <- data.frame(
    input = c("scan_script", "scan_helpers", "fit_result", "parameter_table", "flow_density", "death_data"),
    path = checksum_paths,
    sha256 = vapply(
      checksum_paths,
      .ivt_scan_sha256,
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  .ivt_scan_write_tsv(input_checksums, file.path(out_dir, "input_checksums.tsv"))
  manifest <- data.frame(
    key = c(
      "generated_at", "project_root", "git_branch", "git_head", "git_status",
      "fit_result", "parameter_table", "fit_objects_dir", "flow_density_path", "death_data_path",
      "n_candidates", "n_cores", "design_seed", "batch_size", "resume",
      "calibration_points", "calibration_refine_points", "control_high_mass_limit",
      "calibrated_p_mis_base_upper", "p_mis_base_monotone_nondecreasing",
      "p_mis_base_pass_reentry_after_failure", "fixed_p_wgd", "fixed_gamma_growth",
      "fixed_mu_hp", "fixed_sigma_growth", "n_hard_pass", "n_pareto"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), PROJECT_ROOT,
      .ivt_scan_git_value(c("branch", "--show-current")),
      .ivt_scan_git_value(c("rev-parse", "HEAD")),
      .ivt_scan_git_value(c("status", "--short")),
      fit_result_path, parameter_table_path, fit_objects_dir, flow_density_path, death_data_path,
      n_candidates, n_cores, design_seed, batch_size, resume,
      calibration_points, calibration_refine_points, control_limit,
      calibrated_upper, cap$monotone_nondecreasing[[1L]], cap$pass_reentry_after_failure[[1L]],
      base_params$p_wgd, base_params$gamma_growth, base_params$mu_hp, base_params$sigma_growth,
      sum(final_summary$hard_pass), sum(final_summary$hard_pass & final_summary$pareto_front)
    ),
    stringsAsFactors = FALSE
  )
  .ivt_scan_write_tsv(manifest, file.path(out_dir, "manifest.tsv"))
  writeLines(
    c(
      "status=COMPLETE",
      paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste0("n_candidates=", n_candidates),
      paste0("n_hard_pass=", sum(final_summary$hard_pass)),
      paste0("n_pareto=", sum(final_summary$hard_pass & final_summary$pareto_front))
    ),
    file.path(out_dir, "status.log")
  )
  message("[observed_feasibility] COMPLETE: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
