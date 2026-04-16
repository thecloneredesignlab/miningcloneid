#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEoptim))
suppressPackageStartupMessages(library(dplyr))

.o2g_bootstrap_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2g_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_fit_invivo_backend.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "model", "model_O2G_supply_demand_MAP.R"), local = environment())
rm(.o2g_bootstrap_script_dir)

parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_num_vec <- o2sd_as_num_vec
clip <- o2sd_clip

default_parameter_table_path <- function(script_dir = SCRIPT_DIR) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  normalizePath(file.path(workflow_root, "..", "..", "data", "O2G_supply_demand", "parameter_table.csv"), mustWork = FALSE)
}

default_ploidy_data_path <- function(script_dir = SCRIPT_DIR) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  normalizePath(file.path(workflow_root, "..", "..", "data", "ploidy_distribution.csv"), mustWork = FALSE)
}

default_growth_data_path <- function(script_dir = SCRIPT_DIR) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  normalizePath(file.path(workflow_root, "..", "..", "data", "fit_g.Rds"), mustWork = FALSE)
}

default_out_dir <- function(script_dir = SCRIPT_DIR) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  normalizePath(
    file.path(workflow_root, "..", "..", "results", paste0("fit_model_O2G_supply_demand_MAP_invitro_", stamp)),
    mustWork = FALSE
  )
}

load_ploidy_distribution_data <- function(path) {
  if (!file.exists(path)) stop("Missing ploidy distribution data: ", path)
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("id", "ploidy") %in% names(tab))) {
    stop("ploidy_distribution.csv must contain columns: id, ploidy")
  }
  out <- data.frame(
    id = as.character(tab$id),
    P = suppressWarnings(as.numeric(tab$ploidy)),
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(out$P))) stop("ploidy_distribution.csv contains non-numeric ploidy values.")
  out$hypoxia <- grepl("_O", out$id)
  out$ploidy <- ifelse(grepl("4N", out$id), "4N", "2N")
  out$passage <- 0L
  out$passage[grepl("_A7", out$id) & out$hypoxia] <- 7L
  out$passage[grepl("_A17", out$id) | grepl("_A19", out$id)] <- 17L
  out
}

load_growth_data <- function(path) {
  if (!file.exists(path)) stop("Missing growth-rate data: ", path)
  obj <- readRDS(path)
  req <- c("ploidy", "passage", "g", "o2")
  if (!all(req %in% names(obj))) {
    stop("fit_g.Rds must contain columns: ", paste(req, collapse = ", "))
  }
  out <- as.data.frame(obj, stringsAsFactors = FALSE)
  out$ploidy <- as.character(out$ploidy)
  out$passage <- as.integer(out$passage)
  out$g <- as.numeric(out$g)
  out$o2 <- as.numeric(out$o2)
  out <- out[is.finite(out$passage) & is.finite(out$g) & is.finite(out$o2), , drop = FALSE]
  out
}

make_penalty_components <- function(objective = 1e9, reason = "penalty") {
  empty_dists <- data.frame(
    sim_id = character(0),
    passage = integer(0),
    layer = character(0),
    N = integer(0),
    fraction = numeric(0),
    stringsAsFactors = FALSE
  )
  empty_passage_times <- data.frame(
    sim_id = character(0),
    passage = integer(0),
    duration = numeric(0),
    stringsAsFactors = FALSE
  )
  list(
    objective = as.numeric(objective),
    nll_distribution = as.numeric(objective),
    sse_growth = 0.0,
    nll_terms = numeric(0),
    growth_comparison = data.frame(),
    penalty_reason = as.character(reason),
    sim_output = list(all_dists = empty_dists, all_passage_times = empty_passage_times)
  )
}

calculate_distribution_nll <- function(pred_dist, obs_P_values, N_unit = 22L, missing_penalty = 1e6) {
  if (length(obs_P_values) == 0L) return(0)
  if (nrow(pred_dist) == 0L) return(as.numeric(missing_penalty))
  obs_N <- round(as.numeric(obs_P_values) * as.numeric(N_unit))
  obs_counts <- as.data.frame(table(N = obs_N), stringsAsFactors = FALSE)
  obs_counts$N <- as.integer(as.character(obs_counts$N))
  pred_combined <- pred_dist %>%
    group_by(N) %>%
    summarise(fraction = sum(fraction), .groups = "drop")
  if (nrow(pred_combined) == 0L || any(!is.finite(pred_combined$fraction))) {
    return(as.numeric(missing_penalty))
  }
  comp <- merge(obs_counts, pred_combined, by = "N", all.x = TRUE)
  min_prob <- 1e-10
  comp$fraction[is.na(comp$fraction)] <- min_prob
  -sum(comp$Freq * log(comp$fraction + min_prob))
}

calculate_invitro_objective <- function(sim_output, x_data, growth_data, sim_configs, N_unit, pop_growth_factor) {
  missing_penalty <- 1e6
  all_dists <- sim_output$all_dists
  all_passage_times <- sim_output$all_passage_times
  if (!is.data.frame(all_dists) || !is.data.frame(all_passage_times)) {
    return(make_penalty_components(reason = "sim_output_missing_tables"))
  }
  if (nrow(all_dists) == 0L || nrow(all_passage_times) == 0L) {
    return(make_penalty_components(reason = "sim_output_empty"))
  }
  if (any(!is.finite(all_dists$fraction)) || any(!is.finite(all_passage_times$duration))) {
    return(make_penalty_components(reason = "sim_output_non_finite"))
  }
  dist_pairs <- list(
    list(sim_id = "Sim 1: 2N, O2=20.5", passage = 17L, ploidy = "2N", obs_passage = 0L),
    list(sim_id = "Sim 2: 2N, O2=0", passage = 7L, ploidy = "2N", obs_passage = 7L),
    list(sim_id = "Sim 2: 2N, O2=0", passage = 17L, ploidy = "2N", obs_passage = 17L),
    list(sim_id = "Sim 3: 4N, O2=20.5", passage = 17L, ploidy = "4N", obs_passage = 0L),
    list(sim_id = "Sim 4: 4N, O2=0", passage = 7L, ploidy = "4N", obs_passage = 7L),
    list(sim_id = "Sim 4: 4N, O2=0", passage = 17L, ploidy = "4N", obs_passage = 17L)
  )

  nll_terms <- vapply(
    dist_pairs,
    function(item) {
      pred_dist <- subset(all_dists, sim_id == item$sim_id & passage == item$passage)
      obs_P <- x_data$P[x_data$passage == item$obs_passage & x_data$ploidy == item$ploidy]
      calculate_distribution_nll(pred_dist, obs_P, N_unit = N_unit, missing_penalty = missing_penalty)
    },
    numeric(1)
  )
  total_nll <- sum(nll_terms)

  sim_passage_times <- as.data.frame(all_passage_times, stringsAsFactors = FALSE)
  sim_passage_times$g_sim <- log(pop_growth_factor) / sim_passage_times$duration
  mapping <- data.frame(
    sim_id = vapply(sim_configs, `[[`, character(1), "id"),
    ploidy = vapply(sim_configs, `[[`, character(1), "init_ploidy"),
    o2 = c(1, 0, 1, 0),
    stringsAsFactors = FALSE
  )
  comp_g <- sim_passage_times %>%
    inner_join(mapping, by = "sim_id") %>%
    inner_join(growth_data %>% rename(g_obs = g), by = c("ploidy", "o2", "passage")) %>%
    filter(is.finite(g_sim), is.finite(g_obs), duration > 0)
  expected_growth_rows <- mapping %>%
    inner_join(growth_data %>% distinct(ploidy, o2, passage), by = c("ploidy", "o2"))
  n_missing_growth <- max(nrow(expected_growth_rows) - nrow(comp_g), 0L)
  total_growth_sse <- sum((comp_g$g_sim - comp_g$g_obs)^2) + as.numeric(n_missing_growth) * missing_penalty

  list(
    objective = total_nll + total_growth_sse,
    nll_distribution = total_nll,
    sse_growth = total_growth_sse,
    nll_terms = nll_terms,
    growth_comparison = comp_g
  )
}

lock_invitro_wgd_parameters <- function(parameter_table_path, best_run_params, out_path, glucose = TRUE) {
  tab <- read.csv(parameter_table_path, stringsAsFactors = FALSE, check.names = FALSE)
  wgd_names <- if (isTRUE(glucose)) c("p_wgd_max", "O2_wgd") else c("p_wgd")
  for (nm in wgd_names) {
    idx <- match(nm, tab$param_symbol)
    if (is.na(idx)) stop("Missing parameter_table row for ", nm)
    tab$init_value[idx] <- as.numeric(best_run_params[[nm]])
    tab$estimate[idx] <- FALSE
  }
  write.csv(tab, file = out_path, row.names = FALSE, quote = FALSE)
  invisible(out_path)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  parameter_table <- if (!is.null(argv$parameter_table)) argv$parameter_table else default_parameter_table_path()
  x_data_path <- if (!is.null(argv$x_data)) argv$x_data else default_ploidy_data_path()
  growth_data_path <- if (!is.null(argv$growth_data)) argv$growth_data else default_growth_data_path()
  out_dir <- if (!is.null(argv$out_dir)) argv$out_dir else default_out_dir()
  seed <- as.integer(.first_non_null_local(argv$seed, 1L))
  itermax <- as.integer(.first_non_null_local(argv$itermax, 120L))
  NP <- as.integer(.first_non_null_local(argv$NP, 80L))
  DT <<- as.numeric(.first_non_null_local(argv$dt, 0.05))
  POP_GROWTH_FACTOR <<- as.numeric(.first_non_null_local(argv$pop_growth_factor, 10.0))
  report_passages_arg <- .first_non_null_local(argv$report_passages, "7,17")
  REPORT_PASSAGES <<- as.integer(as_num_vec(report_passages_arg))
  if (length(REPORT_PASSAGES) == 0L || any(!is.finite(REPORT_PASSAGES))) REPORT_PASSAGES <<- c(7L, 17L)
  glucose_use <- canonical_glucose_enabled(
    .first_non_null_local(argv$glucose, TRUE),
    default = TRUE
  )
  glucose_dynamic_use <- canonical_glucose_dynamic(
    .first_non_null_local(argv$glucose_dynamic, FALSE),
    default = FALSE
  )
  if (!isTRUE(glucose_use)) glucose_dynamic_use <- FALSE

  if (!file.exists(parameter_table)) stop("parameter_table not found: ", parameter_table)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(parameter_table, file.path(out_dir, "parameter_table_input.csv"), overwrite = TRUE)
  set.seed(seed)

  x_data <- load_ploidy_distribution_data(x_data_path)
  growth_data <- load_growth_data(growth_data_path)

  x <<- x_data
  m_data <<- growth_data
  N_UNIT <<- 22L
  N_MIN <<- 22L
  N_MAX <<- 154L
  grid_pre <<- N_MIN:N_MAX
  grid_post <<- N_MIN:N_MAX
  PASSAGES_TO_RUN <<- as.integer(max(c(growth_data$passage, REPORT_PASSAGES), na.rm = TRUE))
  if (!is.finite(PASSAGES_TO_RUN) || PASSAGES_TO_RUN < max(REPORT_PASSAGES)) {
    PASSAGES_TO_RUN <<- max(REPORT_PASSAGES)
  }
  sim_configs <<- list(
    list(id = "Sim 1: 2N, O2=20.5", O2 = 20.5, init_ploidy = "2N"),
    list(id = "Sim 2: 2N, O2=0", O2 = 0.0, init_ploidy = "2N"),
    list(id = "Sim 3: 4N, O2=20.5", O2 = 20.5, init_ploidy = "4N"),
    list(id = "Sim 4: 4N, O2=0", O2 = 0.0, init_ploidy = "4N")
  )

  cfg_local <- list(
    parameter_table = parameter_table,
    N_UNIT = N_UNIT,
    N_MIN = N_MIN,
    N_MAX = N_MAX,
    DT = DT,
    fit_treatment = FALSE,
    fit_tau_O2 = FALSE,
    dose_zero_only = TRUE,
    paired_only = FALSE,
    truncate_at_treatment = FALSE,
    ploidy_at_harvest = TRUE,
    O2_growth = TRUE,
    o2_burden_feedback = FALSE,
    Crowding = FALSE,
    crowding = "logistic",
    o2_cache_profile = FALSE,
    ploidy_O2_death = "ploidy_related",
    glucose = glucose_use,
    glucose_dynamic = glucose_dynamic_use,
    glucose_ref_mM = as_num(.first_non_null_local(argv$glucose_ref_mM, default_glucose_ref_mM()), default_glucose_ref_mM()),
    glucose_stress_mode = "coupled_to_O2",
    start_with = "chr_number",
    o2_min = 0.5,
    o2_Nref = 1e6,
    init_total_size = 1.0,
    tau_O2 = 0.1,
    tau_O2_init = 0.1
  )
  cfg_local <- normalize_sim_cfg_common(cfg_local, context = "fit")
  param_bundle <- build_transformed_parameter_table(
    path = parameter_table,
    fit_treatment = FALSE,
    fit_tau_O2 = FALSE,
    O2_growth = TRUE,
    glucose = glucose_use,
    glucose_dynamic = glucose_dynamic_use
  )
  cfg_local <- sync_cfg_from_natural_parameter_table(cfg_local, param_bundle$natural)
  cfg_local$use_soft_prior <- FALSE
  cfg_local$glucose <- glucose_use
  cfg_local$glucose_dynamic <- glucose_dynamic_use
  cfg_local$glucose_stress_mode <- resolve_glucose_runtime_mode(
    glucose_dynamic = glucose_dynamic_use,
    glucose_stress_mode = if (isTRUE(glucose_use)) "coupled_to_O2" else "off",
    default_dynamic = glucose_dynamic_use,
    default_static_mode = if (isTRUE(glucose_use)) "coupled_to_O2" else "off"
  )
  cfg <<- cfg_local

  free_names <- if (isTRUE(glucose_use)) c("log10_p_wgd_max", "log10_O2_wgd") else c("log10_p_wgd")
  full_init_t <- param_bundle$optimizer$init
  transformed_estimate <- setNames(
    as.logical(param_bundle$transformed_output$estimate),
    as.character(param_bundle$transformed_output$param_name)
  )
  if (isTRUE(glucose_dynamic_use)) {
    glucose_candidates <- c("log10_G_S0", "log10_kappa_G", "log10_eta_G", "log10_G_c", "log10_tau_G")
    glucose_candidates <- intersect(glucose_candidates, names(full_init_t))
    glucose_candidates <- glucose_candidates[isTRUE(transformed_estimate[glucose_candidates])]
    free_names <- c(free_names, glucose_candidates)
  }
  lower_free <- as.numeric(param_bundle$optimizer$lower[free_names])
  upper_free <- as.numeric(param_bundle$optimizer$upper[free_names])
  names(lower_free) <- free_names
  names(upper_free) <- free_names

  objective_from_free <- function(par_free_t) {
    full_t <- full_init_t
    full_t[free_names] <- as.numeric(par_free_t)
    run_params <- decode_params(full_t, fit_treatment = FALSE, fit_tau_O2 = FALSE, cfg = cfg_local)
    run_params <- normalize_run_params_common(run_params, cfg = cfg_local)
    run_params$glucose <- glucose_use
    run_params$glucose_dynamic <- glucose_dynamic_use
    run_params$glucose_stress_mode <- resolve_glucose_runtime_mode(
      glucose_dynamic = glucose_dynamic_use,
      glucose_stress_mode = if (isTRUE(glucose_use)) "coupled_to_O2" else "off",
      default_dynamic = glucose_dynamic_use,
      default_static_mode = if (isTRUE(glucose_use)) "coupled_to_O2" else "off"
    )
    comp <- tryCatch(
      {
        sim_output <- run_all_sims(run_params)
        comp_local <- calculate_invitro_objective(
          sim_output,
          x_data = x_data,
          growth_data = growth_data,
          sim_configs = sim_configs,
          N_unit = N_UNIT,
          pop_growth_factor = POP_GROWTH_FACTOR
        )
        comp_local$sim_output <- sim_output
        comp_local
      },
      error = function(e) {
        make_penalty_components(reason = paste0("simulation_error: ", conditionMessage(e)))
      }
    )
    comp$run_params <- run_params
    comp$full_t <- full_t
    comp
  }

  objective_value <- function(par_free_t) {
    objective_from_free(par_free_t)$objective
  }

  de_ctrl <- list(
    trace = TRUE,
    itermax = max(itermax, 1L),
    NP = max(NP, 20L),
    parallelType = "none"
  )
  de_fit <- DEoptim::DEoptim(
    fn = objective_value,
    lower = lower_free,
    upper = upper_free,
    control = de_ctrl
  )
  best_free_t <- as.numeric(de_fit$optim$bestmem)
  names(best_free_t) <- free_names

  local_fit <- suppressWarnings(
    optim(
      par = best_free_t,
      fn = objective_value,
      method = "L-BFGS-B",
      lower = lower_free,
      upper = upper_free,
      control = list(maxit = 200)
    )
  )
  if (is.list(local_fit) && is.finite(local_fit$value) && local_fit$value < de_fit$optim$bestval) {
    best_free_t <- as.numeric(local_fit$par)
    names(best_free_t) <- free_names
  }

  best_comp <- objective_from_free(best_free_t)
  best_run_params <- best_comp$run_params
  best_full_t <- best_comp$full_t
  best_sim_output <- best_comp$sim_output
  best_numeric_params <- best_run_params[vapply(best_run_params, is.numeric, logical(1))]
  best_numeric_params <- best_numeric_params[!vapply(best_numeric_params, is.null, logical(1))]

  best_params_df <- data.frame(
    parameter = names(best_numeric_params),
    value = as.numeric(unlist(best_numeric_params)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  write.table(best_params_df, file = file.path(out_dir, "best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  best_transformed_df <- data.frame(
    transformed_parameter = names(best_full_t),
    transformed_value = as.numeric(best_full_t),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  write.table(best_transformed_df, file = file.path(out_dir, "best_params_transformed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  wgd_df <- if (isTRUE(glucose_use)) {
    data.frame(
      parameter = c("p_wgd_max", "O2_wgd"),
      value = c(as.numeric(best_run_params$p_wgd_max), as.numeric(best_run_params$O2_wgd)),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      parameter = c("p_wgd"),
      value = c(as.numeric(best_run_params$p_wgd)),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
  write.table(wgd_df, file = file.path(out_dir, "wgd_induction_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  sim_dists <- as.data.frame(best_sim_output$all_dists, stringsAsFactors = FALSE)
  sim_dists$P <- sim_dists$N / N_UNIT
  write.table(sim_dists, file = file.path(out_dir, "invitro_ploidy_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  passage_times <- as.data.frame(best_sim_output$all_passage_times, stringsAsFactors = FALSE)
  passage_times$g_sim <- log(POP_GROWTH_FACTOR) / passage_times$duration
  mapping <- data.frame(
    sim_id = vapply(sim_configs, `[[`, character(1), "id"),
    ploidy = vapply(sim_configs, `[[`, character(1), "init_ploidy"),
    o2_pct = c(20.5, 0.0, 20.5, 0.0),
    o2 = c(1, 0, 1, 0),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  passage_times <- merge(passage_times, mapping, by = "sim_id", all.x = TRUE)
  write.table(passage_times, file = file.path(out_dir, "invitro_passage_times.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  ploidy_long <- collapse_layers_to_ploidy(sim_dists, N_UNIT = N_UNIT)
  ploidy_summary <- summarize_ploidy_timecourse(ploidy_long)
  write.table(ploidy_summary, file = file.path(out_dir, "invitro_ploidy_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (nrow(best_comp$growth_comparison) > 0L) {
    write.table(best_comp$growth_comparison, file = file.path(out_dir, "invitro_growth_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  summary_df <- data.frame(
    metric = c(
      "objective_total",
      "objective_distribution_nll",
      "objective_growth_sse",
      "seed",
      "itermax",
      "NP",
      "dt",
      "pop_growth_factor",
      "passages_to_run",
      "report_passages",
      "glucose",
      "glucose_dynamic",
      "glucose_stress_mode",
      "glucose_ref_mM"
    ),
    value = c(
      as.character(best_comp$objective),
      as.character(best_comp$nll_distribution),
      as.character(best_comp$sse_growth),
      as.character(seed),
      as.character(itermax),
      as.character(NP),
      as.character(DT),
      as.character(POP_GROWTH_FACTOR),
      as.character(PASSAGES_TO_RUN),
      paste(REPORT_PASSAGES, collapse = ","),
      as.character(glucose_use),
      as.character(glucose_dynamic_use),
      as.character(cfg_local$glucose_stress_mode),
      as.character(cfg_local$glucose_ref_mM)
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  write.table(summary_df, file = file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  locked_table_path <- file.path(out_dir, "parameter_table.invitro_locked.csv")
  lock_invitro_wgd_parameters(parameter_table, best_run_params, locked_table_path, glucose = glucose_use)

  saveRDS(
    list(
      deoptim = de_fit,
      local_optim = local_fit,
      best_components = best_comp,
      best_params = best_run_params,
      free_names = free_names
    ),
    file = file.path(out_dir, "fit_result.rds")
  )

  message("Done. Results written to: ", normalizePath(out_dir, mustWork = FALSE))
  message("Locked parameter table: ", locked_table_path)
  invisible(normalizePath(out_dir, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  main()
}
