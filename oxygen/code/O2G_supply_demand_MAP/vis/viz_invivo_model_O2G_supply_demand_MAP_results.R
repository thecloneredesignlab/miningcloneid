#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Matrix))

.o2sd_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2sd_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
get_script_dir_self <- o2sd_get_script_dir
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool
clip <- o2sd_clip

# -----------------------------------------------------------------------------
# Function: as_num_vec
# Purpose: Convert an input value to the target scalar/vector type with safe defaults.
# Parameters:
#   - x: Input value or vector to process.
#   - default: Fallback value used when the input is NULL or invalid.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
as_num_vec <- function(x, default = numeric(0)) {
  if (is.null(x)) return(as.numeric(default))
  s <- trimws(as.character(x))
  if (!nzchar(s)) return(as.numeric(default))
  parts <- trimws(unlist(strsplit(s, "[,;]", perl = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(as.numeric(default))
  vals <- suppressWarnings(as.numeric(parts))
  if (any(!is.finite(vals))) stop("Invalid numeric vector argument: ", x)
  vals
}

# -----------------------------------------------------------------------------
# Function: find_latest_fit_dir
# Purpose: Discover fit result directories from a root path.
# Parameters:
#   - results_root: Root directory that contains multiple fit result folders.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
find_latest_fit_dir <- function(results_root) {
  dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
  if (length(dirs) == 0) {
    stop("No fit result directories found under: ", results_root)
  }
  dirs[[which.max(file.info(dirs)$mtime)]]
}

# -----------------------------------------------------------------------------
# Function: normalize_cfg_for_viz
# Purpose: Normalize stored fit configuration for visualization-time simulation calls.
# Parameters:
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
normalize_cfg_for_viz <- function(cfg) {
  normalize_sim_cfg_common(cfg, context = "viz")
}

endpoint_mode_label <- function(cfg) {
  mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  if (identical(mode, "chr_number")) "Chromosome Number (N)" else "Ploidy (2N scale)"
}

weighted_mean_series_label <- function(cfg) {
  mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  if (identical(mode, "chr_number")) {
    "Weighted Mean Chromosome Number (N)"
  } else {
    "Weighted Mean Ploidy (P = N / N_UNIT)"
  }
}

functional_state_axis_label <- function(cfg) {
  mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  if (identical(mode, "chr_number")) "Chromosome Number (N)" else "Ploidy"
}

predict_weighted_metric_label <- function(cfg) {
  mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  if (identical(mode, "chr_number")) "Weighted mean chromosome number" else "Weighted mean ploidy"
}

resource_death_language <- function(glucose_use) {
  if (isTRUE(glucose_use)) {
    return(list(
      component_label = "Dead (Resource stress)",
      adjective = "resource-linked",
      figure_phrase = "dead from resource stress",
      report_phrase = "resource-stress-dead"
    ))
  }
  list(
    component_label = "Dead (Hypoxia)",
    adjective = "hypoxia-linked",
    figure_phrase = "dead from hypoxia",
    report_phrase = "hypoxia-dead"
  )
}

resolve_fixed_glucose_o2_min_pct <- function(cfg = NULL, run_params = NULL) {
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(
      if (!is.null(cfg)) cfg$glucose else NULL,
      if (!is.null(run_params)) run_params$glucose else NULL,
      TRUE
    ),
    default = TRUE
  ))
  if (!isTRUE(glucose_use)) return(NA_real_)

  o2_min_use <- suppressWarnings(as.numeric(.first_non_null_local(
    if (!is.null(cfg)) cfg$o2_min else NULL,
    if (!is.null(run_params)) run_params$o2_min else NULL,
    0
  )))
  if (!is.finite(o2_min_use)) o2_min_use <- 0
  as.numeric(clip(o2_min_use, 0, 100))
}

resolve_fixed_glucose_o2_max_pct <- function(cfg = NULL, run_params = NULL) {
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(
      if (!is.null(cfg)) cfg$glucose else NULL,
      if (!is.null(run_params)) run_params$glucose else NULL,
      TRUE
    ),
    default = TRUE
  ))
  if (!isTRUE(glucose_use)) return(NA_real_)

  o2_max_use <- suppressWarnings(as.numeric(.first_non_null_local(
    if (!is.null(cfg)) cfg$o2_max else NULL,
    if (!is.null(run_params)) run_params$o2_max else NULL,
    if (!is.null(cfg)) cfg$o2_S0_upper_bound else NULL,
    if (!is.null(run_params)) run_params$o2_S0_upper_bound else NULL,
    if (!is.null(cfg)) cfg$o2_S0_max else NULL,
    if (!is.null(run_params)) run_params$o2_S0_max else NULL,
    5
  )))
  if (!is.finite(o2_max_use)) o2_max_use <- 5
  as.numeric(clip(o2_max_use, 0, 100))
}

fixed_glucose_diagnostic_specs <- function(cfg = NULL, run_params = NULL) {
  specs <- data.frame(
    suffix = c("gmin", "gmax", "g20"),
    name = c("o2_min", "o2_max", "20"),
    title_label = c("G Fixed at o2_min", "G Fixed at o2_max", "G Fixed at 20"),
    section_title = c("G = O2_min", "G = O2_max", "G = 20"),
    value = c(
      resolve_fixed_glucose_o2_min_pct(cfg = cfg, run_params = run_params),
      resolve_fixed_glucose_o2_max_pct(cfg = cfg, run_params = run_params),
      20
    ),
    stringsAsFactors = FALSE
  )
  specs$value <- as.numeric(clip(suppressWarnings(as.numeric(specs$value)), 0, 100))
  specs[is.finite(specs$value), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# Function: read_run_params
# Purpose: Read fitted parameter table and reconstruct run_params list.
# Parameters:
#   - fit_dir: Directory containing fitted parameters and summary outputs.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
read_run_params <- function(fit_dir, cfg = NULL) {
  p <- file.path(fit_dir, "best_params.tsv")
  if (!file.exists(p)) stop("Missing best_params.tsv in: ", fit_dir)
  tab <- read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain columns: parameter, value")
  }
  vals <- setNames(as.numeric(tab$value), as.character(tab$parameter))
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(cfg$glucose, TRUE),
    default = TRUE
  ))
  if (isTRUE(glucose_use) && !"k_o" %in% names(vals)) {
    # Glucose resource-growth fits omit k_o from output; keep a finite
    # compatibility value for legacy simulation/viz interfaces.
    k_o_compat <- as.numeric(.first_non_null_local(cfg$k_o_init, 50.0))
    if (!is.finite(k_o_compat) || k_o_compat <= 0) k_o_compat <- 50.0
    vals[["k_o"]] <- k_o_compat
  }
  needed_common <- c(
    "lam_min", "lam_max", "k_o", "p_misseg", "k_o_mis",
    "o2_S0", "kappa_O", "eta_o2",
    "alpha_o2", "gamma_growth",
    "mu_hp", "gamma_mu", "O2_crit", "n_O", "k_clear"
  )
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(cfg$misseg_loss_survival, "nullisomy"),
    "nullisomy"
  )
  needed <- c(
    needed_common,
    if (identical(loss_mode, "buffering")) c("buffer_smax", "buffer_beta", "buffer_n_exp") else c("gamma_loss"),
    "p_wgd"
  )
  miss <- setdiff(needed, names(vals))
  if (length(miss) > 0) {
    stop("best_params.tsv missing parameters: ", paste(miss, collapse = ", "))
  }
  out <- as.list(vals[needed])
  extra_init_mult <- grep("^init_mult_", names(vals), value = TRUE)
  if (length(extra_init_mult) > 0L) {
    for (nm in extra_init_mult) out[[nm]] <- vals[[nm]]
  }
  p_mis_base_val <- if ("p_mis_base" %in% names(vals)) vals[["p_mis_base"]] else NULL
  o2_min_val <- if ("o2_min" %in% names(vals)) vals[["o2_min"]] else NULL
  if (!is.null(p_mis_base_val) && is.finite(p_mis_base_val)) out$p_mis_base <- as.numeric(p_mis_base_val)
  if (!is.null(o2_min_val) && is.finite(o2_min_val)) out$o2_min <- as.numeric(o2_min_val)
  if ("rho_2N" %in% names(vals) && is.finite(vals[["rho_2N"]]) && vals[["rho_2N"]] > 0) {
    out$rho_2N <- vals[["rho_2N"]]
  }
  if ("c_vol_2N_mm3" %in% names(vals) && is.finite(vals[["c_vol_2N_mm3"]]) && vals[["c_vol_2N_mm3"]] > 0) {
    out$c_vol_2N_mm3 <- vals[["c_vol_2N_mm3"]]
  }
  fit_treatment_use <- isTRUE(.first_non_null_local(cfg$fit_treatment, FALSE))
  if (fit_treatment_use) {
    miss_tx <- setdiff(c("alpha", "gamma"), names(vals))
    if (length(miss_tx) > 0) {
      stop("best_params.tsv missing treatment parameters while fit_treatment=TRUE: ", paste(miss_tx, collapse = ", "))
    }
    out$alpha <- vals[["alpha"]]
    out$gamma <- vals[["gamma"]]
  } else {
    out$alpha <- if ("alpha" %in% names(vals) && is.finite(vals[["alpha"]])) vals[["alpha"]] else 0
    out$gamma <- if ("gamma" %in% names(vals) && is.finite(vals[["gamma"]])) vals[["gamma"]] else 1
  }
  out$tau_O2 <- if ("tau_O2" %in% names(vals) && is.finite(vals[["tau_O2"]]) && vals[["tau_O2"]] > 0) vals[["tau_O2"]] else NULL
  normalize_run_params_common(out, cfg = cfg)
}

# -----------------------------------------------------------------------------
# Function: simulate_one_full
# Purpose: Run one forward simulation trajectory for a single scenario.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - scenario: Single scenario data object.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - report_dt: Sampling interval for reported trajectories.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
simulate_one_full <- function(run_params, scenario, cfg, report_dt = 1.0) {
  model_core <- build_model_core(cfg = cfg)
  grid_pre <- model_core$grid_pre
  init_state <- if (scenario$cohort == "2N") model_core$init_state_2N else model_core$init_state_4N
  init_mult_name <- .first_non_null_local(scenario$init_mult_param, harvest_init_natural_param_name(scenario$harvest))
  init_mult <- suppressWarnings(as.numeric(.first_non_null_local(run_params[[init_mult_name]], 1.0)))
  if (!is.finite(init_mult) || init_mult <= 0) {
    log_name <- .first_non_null_local(scenario$log_init_mult_param, harvest_init_log_param_name(scenario$harvest))
    init_mult <- exp(as.numeric(.first_non_null_local(run_params[[log_name]], 0.0)))
  }
  if (!is.finite(init_mult) || init_mult <= 0) init_mult <- 1.0
  init_state <- as.numeric(init_state) * init_mult
  sim_end_day <- as.numeric(scenario$sim_end_day)
  full_steps <- 0:as.integer(round(sim_end_day / cfg$DT))
  full_days <- as.numeric(full_steps) * cfg$DT
  keep_days <- sort(unique(c(
    0,
    sim_end_day,
    as.numeric(scenario$obs_days),
    seq(0, sim_end_day, by = report_dt)
  )))
  keep_days <- keep_days[is.finite(keep_days) & keep_days >= 0 & keep_days <= sim_end_day]
  keep_steps <- sort(unique(as.integer(round(keep_days / cfg$DT))))

  obs_steps_cpp <- as.integer(full_steps)
  sim_end_step_cpp <- max(obs_steps_cpp)
  o2_s0_upper_use <- as.numeric(.first_non_null_local(cfg$o2_S0_upper_bound, 5.0))
  if (!is.finite(o2_s0_upper_use) || o2_s0_upper_use <= 0) o2_s0_upper_use <- 5.0
  o2_S0_use <- as.numeric(.first_non_null_local(run_params$o2_S0, cfg$o2_S0_init, 0.5))
  if (!is.finite(o2_S0_use) || o2_S0_use < 0) o2_S0_use <- 0.5
  o2_S0_use <- min(max(o2_S0_use, 0), o2_s0_upper_use)
  kappa_O_use <- as.numeric(.first_non_null_local(run_params$kappa_O, cfg$kappa_O_init, 1.0))
  if (!is.finite(kappa_O_use) || kappa_O_use <= 0) kappa_O_use <- 1.0
  eta_o2_use <- as.numeric(.first_non_null_local(run_params$eta_o2, cfg$eta_o2_init, 1.0))
  if (!is.finite(eta_o2_use) || eta_o2_use <= 0) eta_o2_use <- 1.0
  O2_crit_use <- as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  o2_Nref_use <- as.numeric(.first_non_null_local(cfg$o2_Nref, cfg$init_total_size, 1e6))
  if (!is.finite(o2_Nref_use) || o2_Nref_use <= 0) o2_Nref_use <- 1e6
  o2_min_use <- as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0.5))
  if (!is.finite(o2_min_use) || o2_min_use < 0) o2_min_use <- 0.5
  o2_min_use <- min(max(o2_min_use, 0), o2_s0_upper_use)
  tau_O2_use <- as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 0)
  o2_growth_use <- isTRUE(.first_non_null_local(cfg$O2_growth, TRUE))
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(cfg$glucose, run_params$glucose, TRUE),
    default = TRUE
  ))
  death_language <- resource_death_language(glucose_use)
  loss_mode_use <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(run_params$misseg_loss_survival, cfg$misseg_loss_survival, "nullisomy"),
    "nullisomy"
  )
  buffer_smax_use <- as.numeric(.first_non_null_local(run_params$buffer_smax, cfg$buffer_smax_init, 1.0))
  if (!is.finite(buffer_smax_use)) buffer_smax_use <- 1.0
  buffer_beta_use <- as.numeric(.first_non_null_local(run_params$buffer_beta, cfg$buffer_beta_init, 0.0))
  if (!is.finite(buffer_beta_use)) buffer_beta_use <- 0.0
  buffer_n_exp_use <- as.numeric(.first_non_null_local(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1.0))
  if (!is.finite(buffer_n_exp_use) || buffer_n_exp_use <= 0) buffer_n_exp_use <- 1.0

  sim_cpp <- cpp_o2simps_simulate_one(list(
    init_state = as.numeric(init_state),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = as.integer(obs_steps_cpp),
    sim_end_step = as.integer(sim_end_step_cpp),
    DT = as.numeric(cfg$DT),
    dose = as.numeric(scenario$dose),
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = as.numeric(scenario$treat_day),
    fit_treatment = isTRUE(cfg$fit_treatment),
    alpha = as.numeric(.first_non_null_local(run_params$alpha, 0.0)),
    gamma = as.numeric(.first_non_null_local(run_params$gamma, 1.0)),
    tx_mult_min = as.numeric(cfg$tx_mult_min),
    crowding_enabled = isTRUE(cfg$Crowding),
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_crit = as.numeric(O2_crit_use),
    o2_feedback = isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE)),
    o2_S0 = as.numeric(o2_S0_use),
    kappa_O = as.numeric(kappa_O_use),
    tau_O2 = as.numeric(tau_O2_use),
    o2_Nref = as.numeric(o2_Nref_use),
    o2_min = as.numeric(o2_min_use),
    eta_o2 = as.numeric(eta_o2_use),
    o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01)),
    o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct, 0.005)),
    o2_cache_profile = isTRUE(.first_non_null_local(cfg$o2_cache_profile, FALSE)),
    glucose = isTRUE(glucose_use),
    O2_growth = isTRUE(o2_growth_use),
    lam_min = as.numeric(run_params$lam_min),
    lam_max = as.numeric(run_params$lam_max),
    k_o = as.numeric(run_params$k_o),
    has_p_misseg = !is.null(run_params$p_misseg),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = 0.0,
    p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0.0)),
    p_wgd_max = as.numeric(.first_non_null_local(run_params$p_wgd_max, 0.0)),
    O2_wgd = as.numeric(.first_non_null_local(run_params$O2_wgd, cfg$O2_wgd_init, 0.1)),
    boundary = as.character(.first_non_null_local(run_params$boundary, "drop")),
    eps_tail = as.numeric(1e-8),
    gamma_loss = as.numeric(.first_non_null_local(run_params$gamma_loss, 0.1)),
    misseg_loss_survival = as.character(loss_mode_use),
    buffer_smax = as.numeric(buffer_smax_use),
    buffer_beta = as.numeric(buffer_beta_use),
    buffer_n_exp = as.numeric(buffer_n_exp_use),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = 0.0,
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
    mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)),
    gamma_mu = as.numeric(.first_non_null_local(run_params$gamma_mu, cfg$gamma_mu_init, 1.0)),
    n_O = as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0)),
    ploidy_O2_death = assert_canonical_ploidy_o2_death_mode(
      .first_non_null_local(cfg$ploidy_O2_death, run_params$ploidy_O2_death, "diploid_NULL")
    ),
    start_with = assert_canonical_start_with_mode(
      .first_non_null_local(cfg$start_with, "ploidy")
    ),
    k_clear = as.numeric(.first_non_null_local(run_params$k_clear, cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(burden_floor),
    return_full_trajectory = TRUE
  ))

  live_state_obs <- sim_cpp$live_state_obs
  if (is.null(live_state_obs) || length(live_state_obs) == 0) {
    return(list(burden = data.frame(), ploidy = data.frame()))
  }
  live_state_obs <- as.matrix(live_state_obs)
  if (!identical(dim(live_state_obs), c(length(obs_steps_cpp), length(grid_pre)))) {
    stop(
      "cpp_o2simps_simulate_one returned live_state_obs with unexpected shape: got ",
      paste(dim(live_state_obs), collapse = "x"),
      ", expected ", length(obs_steps_cpp), "x", length(grid_pre)
    )
  }
  live_pop_by_step <- rowSums(live_state_obs, na.rm = TRUE)
  live_frac_obs <- live_state_obs / pmax(live_pop_by_step, 1e-12)
  ploidy_rows_full <- data.frame(
    harvest = scenario$harvest,
    cohort = scenario$cohort,
    dose = scenario$dose,
    day = rep(full_days, each = length(grid_pre)),
    step = rep(obs_steps_cpp, each = length(grid_pre)),
    N = rep(as.numeric(grid_pre), times = length(obs_steps_cpp)),
    fraction = as.numeric(t(live_frac_obs)),
    pop = rep(as.numeric(live_pop_by_step), each = length(grid_pre))
  )

  burden_by_day_full <- data.frame(
    day = full_days,
    step = obs_steps_cpp,
    pred_burden_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_total_obs, sim_cpp$Ntot_obs)),
    pred_burden_live_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_live_obs, sim_cpp$Ntot_obs)),
    pred_burden_dead_hypoxia_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_hypoxia_obs, rep(0, length(obs_steps_cpp)))),
    pred_burden_dead_buffer_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_buffer_obs, rep(0, length(obs_steps_cpp)))),
    pred_burden_dead_total_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_total_obs, rep(0, length(obs_steps_cpp)))),
    pred_burden_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_total_obs, sim_cpp$Vmm3_obs)),
    pred_burden_live_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_live_obs, sim_cpp$Vmm3_obs)),
    pred_burden_dead_hypoxia_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_hypoxia_obs, rep(0, length(obs_steps_cpp)))),
    pred_burden_dead_buffer_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_buffer_obs, rep(0, length(obs_steps_cpp)))),
    pred_burden_dead_total_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_total_obs, rep(0, length(obs_steps_cpp)))),
    pred_o2_target_pct = as.numeric(.first_non_null_local(sim_cpp$O2_target_obs, rep(NA_real_, length(obs_steps_cpp)))),
    pred_o2_pct = as.numeric(.first_non_null_local(sim_cpp$O2_eff_obs, rep(NA_real_, length(obs_steps_cpp)))),
    pred_g_target_pct = as.numeric(.first_non_null_local(sim_cpp$G_target_obs, rep(NA_real_, length(obs_steps_cpp)))),
    pred_g_pct = as.numeric(.first_non_null_local(sim_cpp$G_eff_obs, rep(NA_real_, length(obs_steps_cpp))))
  ) %>% arrange(step)
  burden_by_day_full$pred_o2_lag_gap_pct <- as.numeric(
    burden_by_day_full$pred_o2_target_pct - burden_by_day_full$pred_o2_pct
  )
  burden_by_day <- burden_by_day_full %>% filter(step %in% keep_steps)
  ploidy_rows <- ploidy_rows_full %>% filter(step %in% keep_steps)

  obs_steps <- as.integer(round(as.numeric(scenario$obs_days) / cfg$DT))
  obs_map <- setNames(as.numeric(scenario$obs_burden), as.character(obs_steps))
  burden_by_day$obs_burden <- as.numeric(obs_map[as.character(burden_by_day$step)])

  burden_rows <- burden_by_day %>%
    mutate(
      harvest = scenario$harvest,
      cohort = scenario$cohort,
      dose = scenario$dose,
      treat_day = scenario$treat_day,
      pred_burden = pred_burden_volume_mm3
    ) %>%
    select(
      harvest, cohort, dose, treat_day, step, day,
      pred_burden, pred_burden_volume_mm3, pred_burden_cells,
      pred_burden_live_volume_mm3, pred_burden_dead_hypoxia_volume_mm3,
      pred_burden_dead_buffer_volume_mm3, pred_burden_dead_total_volume_mm3,
      pred_burden_live_cells, pred_burden_dead_hypoxia_cells,
      pred_burden_dead_buffer_cells, pred_burden_dead_total_cells,
      pred_o2_target_pct, pred_o2_pct, pred_o2_lag_gap_pct,
      pred_g_target_pct, pred_g_pct,
      obs_burden
    )

  list(
    burden = as.data.frame(burden_rows),
    ploidy = as.data.frame(select(ploidy_rows, harvest, cohort, dose, day, N, fraction))
  )
}

# -----------------------------------------------------------------------------
# Function: simulate_one_full_horizon
# Purpose: Run one forward simulation trajectory for a single scenario.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - scenario: Single scenario data object.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - horizon_day: Prediction horizon end time in days.
#   - report_dt: Sampling interval for reported trajectories.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
simulate_one_full_horizon <- function(run_params, scenario, cfg, horizon_day, report_dt = 1.0) {
  sc <- scenario
  sc$sim_end_day <- as.numeric(max(horizon_day, 0))
  simulate_one_full(run_params, sc, cfg, report_dt = report_dt)
}

# -----------------------------------------------------------------------------
# Function: normalize_burden_for_plot
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - burden_all: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
normalize_burden_for_plot <- function(burden_all) {
  burden_all %>%
    group_by(harvest, cohort, dose) %>%
    arrange(day, .by_group = TRUE) %>%
    group_modify(function(df, .y) {
      pred_delta <- df$pred_burden - df$pred_burden[[1]]
      pred_scale <- max(abs(pred_delta), na.rm = TRUE)
      df$pred_norm <- if (is.finite(pred_scale) && pred_scale > 0) pred_delta / pred_scale else pred_delta

      obs_vals <- df$obs_burden[is.finite(df$obs_burden)]
      if (length(obs_vals) > 0) {
        obs_delta <- df$obs_burden - obs_vals[[1]]
        obs_scale <- max(abs(obs_delta), na.rm = TRUE)
        df$obs_norm <- if (is.finite(obs_scale) && obs_scale > 0) obs_delta / obs_scale else obs_delta
      } else {
        df$obs_norm <- NA_real_
      }
      df
    }) %>%
    ungroup()
}

# -----------------------------------------------------------------------------
# Function: compute_ploidy_weighted_mean
# Purpose: Compute weighted mean with finite/positive-weight safeguards.
# Parameters:
#   - ploidy_all: Function-specific input argument.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
compute_ploidy_weighted_mean <- function(ploidy_all, cfg) {
  start_with_mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  ploidy_all %>%
    group_by(harvest, cohort, dose, day) %>%
    summarise(
      weighted_mean_N = sum(N * fraction, na.rm = TRUE) / pmax(sum(fraction, na.rm = TRUE), 1e-12),
      .groups = "drop"
    ) %>%
    mutate(
      weighted_mean_ploidy = if (identical(start_with_mode, "chr_number")) {
        weighted_mean_N
      } else {
        weighted_mean_N / cfg$N_UNIT
      },
      weighted_mean_endpoint = weighted_mean_ploidy,
      start_with = start_with_mode
    )
}

# -----------------------------------------------------------------------------
# Function: build_terminal_ploidy_compare_df
# Purpose: Build observed-vs-predicted ploidy distributions at the time points
#   used by the ploidy objective.
# -----------------------------------------------------------------------------
build_terminal_ploidy_compare_df <- function(scenarios, ploidy_all, cfg) {
  start_with_mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  meta_rows <- vector("list", length(scenarios))
  obs_rows <- vector("list", length(scenarios))

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    obs_value <- if (identical(start_with_mode, "chr_number")) {
      as.numeric(sc$chr_number_obs)
    } else {
      obs_z_raw <- as.numeric(sc$ploidy_obs_z)
      obs_N <- map_ploidy_to_N_by_chrlen(
        ploidy_values = obs_z_raw / as.numeric(cfg$N_UNIT),
        N_grid = cfg$N_MIN:cfg$N_MAX,
        chr_lengths_bp = cfg$chr_lengths_bp
      )
      obs_N <- as.integer(clip(obs_N, cfg$N_MIN, cfg$N_MAX))
      as.numeric(obs_N) / as.numeric(cfg$N_UNIT)
    }
    obs_value <- obs_value[is.finite(obs_value)]
    if (length(obs_value) == 0L) next

    target_day <- as.numeric(sc$sim_end_day)
    if (!is.finite(target_day)) next

    meta_rows[[i]] <- data.frame(
      harvest = as.character(sc$harvest),
      cohort = as.character(sc$cohort),
      dose = as.numeric(sc$dose),
      target_day = target_day,
      stringsAsFactors = FALSE
    )

    obs_rows[[i]] <- data.frame(
      harvest = as.character(sc$harvest),
      cohort = as.character(sc$cohort),
      dose = as.numeric(sc$dose),
      target_day = target_day,
      source = "Observed",
      ploidy = as.numeric(obs_value),
      endpoint_value = as.numeric(obs_value),
      endpoint_mode = start_with_mode,
      weight = 1,
      stringsAsFactors = FALSE
    )
  }

  meta_df <- bind_rows(meta_rows)
  obs_df <- bind_rows(obs_rows)
  if (nrow(meta_df) == 0L || nrow(obs_df) == 0L || nrow(ploidy_all) == 0L) {
    return(data.frame())
  }

  pred_df <- ploidy_all %>%
    inner_join(meta_df, by = c("harvest", "cohort", "dose")) %>%
    mutate(day_gap = abs(as.numeric(day) - as.numeric(target_day))) %>%
    group_by(harvest, cohort, dose, target_day) %>%
    filter(day_gap == min(day_gap, na.rm = TRUE)) %>%
    ungroup() %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      target_day = as.numeric(target_day),
      source = "Predicted",
      ploidy = if (identical(start_with_mode, "chr_number")) as.numeric(N) else as.numeric(N) / as.numeric(cfg$N_UNIT),
      endpoint_value = if (identical(start_with_mode, "chr_number")) as.numeric(N) else as.numeric(N) / as.numeric(cfg$N_UNIT),
      endpoint_mode = start_with_mode,
      weight = as.numeric(fraction)
    ) %>%
    filter(is.finite(ploidy), is.finite(weight), weight > 0)

  bind_rows(obs_df, pred_df) %>%
    mutate(
      cohort = factor(as.character(cohort), levels = c("2N", "4N")),
      source = factor(as.character(source), levels = c("Observed", "Predicted")),
      fill_group = factor(
        paste(as.character(cohort), as.character(source)),
        levels = c("2N Observed", "2N Predicted", "4N Observed", "4N Predicted")
      )
    )
}

# -----------------------------------------------------------------------------
# Function: plot_terminal_ploidy_violin_compare
# Purpose: Compare observed and predicted ploidy distributions used in the
#   ploidy objective, grouped by 2N/4N cohort.
# -----------------------------------------------------------------------------
plot_terminal_ploidy_violin_compare <- function(compare_df, fit_dir, out_dir) {
  if (nrow(compare_df) == 0L) return(NULL)
  endpoint_label <- unique(compare_df$endpoint_mode)
  endpoint_label <- if (length(endpoint_label) == 1L && identical(endpoint_label[[1]], "chr_number")) {
    "Chromosome Number (N)"
  } else {
    "Ploidy (2N scale)"
  }

  weighted_quantile_local <- function(x, w, probs) {
    x <- as.numeric(x)
    w <- as.numeric(w)
    keep <- is.finite(x) & is.finite(w) & (w > 0)
    x <- x[keep]
    w <- w[keep]
    if (!length(x)) return(rep(NA_real_, length(probs)))
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    cw <- cumsum(w) / sum(w)
    vapply(
      probs,
      function(p) {
        idx <- which(cw >= p)[1]
        if (!length(idx) || is.na(idx)) x[[length(x)]] else x[[idx]]
      },
      numeric(1)
    )
  }

  box_df <- compare_df %>%
    group_by(cohort, source, fill_group) %>%
    summarise(
      q1 = weighted_quantile_local(endpoint_value, weight, 0.25),
      median = weighted_quantile_local(endpoint_value, weight, 0.5),
      q3 = weighted_quantile_local(endpoint_value, weight, 0.75),
      ymin_raw = min(endpoint_value[weight > 0], na.rm = TRUE),
      ymax_raw = max(endpoint_value[weight > 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      iqr = q3 - q1,
      ymin = pmax(ymin_raw, q1 - 1.5 * iqr),
      ymax = pmin(ymax_raw, q3 + 1.5 * iqr)
    )

  fill_values <- c(
    "2N Observed" = "#8FBF7A",
    "2N Predicted" = "#D7EFCC",
    "4N Observed" = "#F2A9BC",
    "4N Predicted" = "#F9D9E3"
  )

  p <- ggplot(compare_df, aes(x = source, y = endpoint_value, weight = weight, fill = fill_group)) +
    geom_violin(trim = FALSE, scale = "width", quantiles = NULL, color = "grey35", linewidth = 0.35, alpha = 0.95) +
    geom_boxplot(
      data = box_df,
      aes(x = source, ymin = ymin, lower = q1, middle = median, upper = q3, ymax = ymax),
      stat = "identity",
      inherit.aes = FALSE,
      width = 0.18,
      fill = "white",
      color = "black",
      linewidth = 0.45,
      outlier.shape = NA
    ) +
    facet_wrap(~ cohort, nrow = 1, ncol = 2) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    labs(
      title = paste0("Observed vs Predicted ", endpoint_label, " Distributions Used in Endpoint Objective"),
      subtitle = if (identical(endpoint_label, "Chromosome Number (N)")) {
        paste0("fit_dir=", basename(fit_dir))
      } else {
        paste0("Observed ploidy is mapped to the chromosome-count grid used by the objective | fit_dir=", basename(fit_dir))
      },
      x = NULL,
      y = endpoint_label,
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )

  ggsave(file.path(out_dir, "terminal_ploidy_observed_vs_predicted_violin.pdf"), p, width = 6.5, height = 6.5)
  p
}

# -----------------------------------------------------------------------------
# Function: plot_functional_response_curves
# Purpose: Generate and save visualization output for fitted model behavior.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - out_dir: Output directory for generated files and plots.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
plot_functional_response_curves <- function(run_params, cfg, out_dir, fixed_glucose_pct = NULL, fixed_glucose_specs = NULL) {
  start_with_mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  o2_plot_min <- 0
  o2_plot_max <- 5
  o2_grid <- seq(o2_plot_min, o2_plot_max, by = 0.02)
  glucose_heatmap_grid <- seq(0, 100, by = 0.5)
  glucose_heatmap_o2range_grid <- o2_grid
  state_plot_min <- if (identical(start_with_mode, "chr_number")) as.numeric(cfg$N_MIN) else 0
  state_plot_max <- if (identical(start_with_mode, "chr_number")) as.numeric(cfg$N_MAX) else 10
  state_grid_dense <- if (identical(start_with_mode, "chr_number")) {
    seq(state_plot_min, state_plot_max, by = 1)
  } else {
    seq(state_plot_min, state_plot_max, by = 0.05)
  }
  o2_levels_ploidy <- seq(o2_plot_min, o2_plot_max, by = 0.5)
  state_axis_label <- functional_state_axis_label(cfg)
  ref_df <- data.frame(
    cohort = c("2N", "4N"),
    N_ref = as.numeric(c(2 * cfg$N_UNIT, 4 * cfg$N_UNIT)),
    stringsAsFactors = FALSE
  )
  ref_state_mult <- seq(1.5, 5.0, by = 0.5)
  ref_state_label <- ifelse(
    abs(ref_state_mult - round(ref_state_mult)) < 1e-8,
    paste0(as.integer(round(ref_state_mult)), "N"),
    paste0(format(ref_state_mult, trim = TRUE, nsmall = 1), "N")
  )
  ref_state_subtitle <- paste0("Reference states: ", paste(ref_state_label, collapse = ", "))
  ref_df_multi <- data.frame(
    cohort = ref_state_label,
    ploidy_multiple = ref_state_mult,
    N_ref = as.numeric(ref_state_mult * as.numeric(cfg$N_UNIT)),
    stringsAsFactors = FALSE
  )

  O2_crit_use <- as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  o2_growth_use <- isTRUE(.first_non_null_local(cfg$O2_growth, TRUE))
  alpha_o2_use <- pmax(0, as.numeric(run_params$alpha_o2))
  gamma_growth_use <- pmax(as.numeric(run_params$gamma_growth), 1e-12)
  gamma_mu_use <- pmax(as.numeric(run_params$gamma_mu), 1e-12)
  ploidy_O2_death_use <- assert_canonical_ploidy_o2_death_mode(
    .first_non_null_local(cfg$ploidy_O2_death, run_params$ploidy_O2_death, "diploid_NULL")
  )
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(cfg$glucose, run_params$glucose, TRUE),
    default = TRUE
  ))
  death_language <- resource_death_language(glucose_use)
  if (is.null(fixed_glucose_specs) || !is.data.frame(fixed_glucose_specs)) {
    fixed_glucose_pct <- suppressWarnings(as.numeric(fixed_glucose_pct))
    fixed_glucose_pct <- if (length(fixed_glucose_pct) > 0L) fixed_glucose_pct[[1L]] else NA_real_
    fixed_glucose_specs <- data.frame(
      suffix = "gmin",
      name = "o2_min",
      title_label = "G Fixed at o2_min",
      section_title = "G = O2_min",
      value = fixed_glucose_pct,
      stringsAsFactors = FALSE
    )
  }
  fixed_glucose_specs$value <- as.numeric(clip(suppressWarnings(as.numeric(fixed_glucose_specs$value)), 0, 100))
  fixed_glucose_specs <- fixed_glucose_specs[is.finite(fixed_glucose_specs$value), , drop = FALSE]
  has_fixed_glucose <- isTRUE(glucose_use) && nrow(fixed_glucose_specs) > 0L
  n_O <- as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0))
  if (!is.finite(n_O) || n_O < 0) stop("run_params$n_O must be finite and >= 0.")
  mu_hp_use <- pmax(as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)), 0)
  gamma_loss_use <- as.numeric(.first_non_null_local(run_params$gamma_loss, 0.1))
  if (!is.finite(gamma_loss_use) || gamma_loss_use <= 0) gamma_loss_use <- 0.1
  boundary_mode_use <- as.character(.first_non_null_local(run_params$boundary, "drop"))
  p_wgd_use <- as.numeric(.first_non_null_local(run_params$p_wgd, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
  p_wgd_max_use <- as.numeric(.first_non_null_local(run_params$p_wgd_max, 0.0))
  if (!is.finite(p_wgd_max_use)) p_wgd_max_use <- 0.0
  O2_wgd_use <- as.numeric(.first_non_null_local(run_params$O2_wgd, cfg$O2_wgd_init, 0.1))
  if (!is.finite(O2_wgd_use) || O2_wgd_use <= 0) O2_wgd_use <- 1e-12
  loss_mode_use <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(run_params$misseg_loss_survival, cfg$misseg_loss_survival, "nullisomy"),
    "nullisomy"
  )
  buffer_smax_use <- as.numeric(.first_non_null_local(run_params$buffer_smax, cfg$buffer_smax_init, 1.0))
  if (!is.finite(buffer_smax_use)) buffer_smax_use <- 1.0
  buffer_beta_use <- as.numeric(.first_non_null_local(run_params$buffer_beta, cfg$buffer_beta_init, 0.0))
  if (!is.finite(buffer_beta_use)) buffer_beta_use <- 0.0
  buffer_n_exp_use <- as.numeric(.first_non_null_local(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1.0))
  if (!is.finite(buffer_n_exp_use) || buffer_n_exp_use <= 0) buffer_n_exp_use <- 1.0
  N_dip <- 44.0
  n_state <- as.integer(cfg$N_MAX) - as.integer(cfg$N_MIN) + 1L
  .require_cpp_o2simps_fn("cpp_o2simps_build_G_for_o2_triplet")
  functional_rate_cache <- new.env(parent = emptyenv())

  get_functional_rate_curves <- function(O2) {
    O2_use <- .assert_o2_pct(as.numeric(O2), label = "O2")
    key <- sprintf("%.6f", O2_use)
    if (!exists(key, envir = functional_rate_cache, inherits = FALSE)) {
      tri <- cpp_o2simps_build_G_for_o2_triplet(
        O2 = as.numeric(O2_use),
        O2_crit = as.numeric(O2_crit_use),
        N0min = as.integer(cfg$N_MIN),
        N0max = as.integer(cfg$N_MAX),
        N1min = as.integer(cfg$N_MIN),
        N1max = as.integer(cfg$N_MAX),
        lam_min = as.numeric(run_params$lam_min),
        lam_max = as.numeric(run_params$lam_max),
        k_o = as.numeric(run_params$k_o),
        has_p_misseg = !is.null(run_params$p_misseg),
        p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init, 1e-5)),
        p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
        k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
        has_pmis_endpoints = FALSE,
        pmis_O2_0 = 0.0,
        pmis_O2_1 = 0.0,
        p_const = 0.0,
        glucose = isTRUE(glucose_use),
        p_wgd = as.numeric(p_wgd_use),
        p_wgd_max = as.numeric(p_wgd_max_use),
        O2_wgd = as.numeric(O2_wgd_use),
        boundary = as.character(boundary_mode_use),
        eps_tail = as.numeric(1e-8),
        gamma_loss = as.numeric(gamma_loss_use),
        misseg_loss_survival = as.character(loss_mode_use),
        buffer_smax = as.numeric(buffer_smax_use),
        buffer_beta = as.numeric(buffer_beta_use),
        buffer_n_exp = as.numeric(buffer_n_exp_use),
        N_unit = as.integer(cfg$N_UNIT),
        beta_size = 0.0,
        O2_growth = isTRUE(o2_growth_use),
        alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
        gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
        mu_hp = as.numeric(mu_hp_use),
        gamma_mu = as.numeric(gamma_mu_use),
        n_O = as.numeric(n_O),
        ploidy_O2_death = as.character(ploidy_O2_death_use)
      )
      curve_names <- c(
        "dead_buffer_rate",
        "nullisomy_nonviable_rate",
        "boundary_dropped_rate",
        "nullisomy_nonviable_division_prob",
        "nullisomy_nonviable_daughters_per_division"
      )
      curve_list <- setNames(vector("list", length(curve_names)), curve_names)
      for (nm in curve_names) {
        vals <- as.numeric(tri[[nm]])
        if (length(vals) != n_state) {
          stop(
            "cpp_o2simps_build_G_for_o2_triplet returned ",
            nm,
            " with length ",
            length(vals),
            "; expected ",
            n_state,
            "."
          )
        }
        curve_list[[nm]] <- vals
      }
      curve_list$nullisomy_nonviable_daughter_fraction <- pmax(
        pmin(0.5 * curve_list$nullisomy_nonviable_daughters_per_division, 1),
        0
      )
      assign(key, curve_list, envir = functional_rate_cache)
    }
    get(key, envir = functional_rate_cache, inherits = FALSE)
  }

  compute_rate_components <- function(O2, N, G = NULL) {
    O2_vec <- .assert_o2_pct(as.numeric(O2), label = "O2")
    N_vec <- as.numeric(N)
    G_input <- if (is.null(G)) NULL else .assert_o2_pct(as.numeric(G), label = "G")
    n_out <- max(length(O2_vec), length(N_vec), if (is.null(G_input)) 1L else length(G_input))
    if (!(length(O2_vec) %in% c(1L, n_out) &&
          length(N_vec) %in% c(1L, n_out) &&
          (is.null(G_input) || length(G_input) %in% c(1L, n_out)))) {
      stop("O2, G, and N must have compatible lengths in compute_rate_components().")
    }
    O2_vec <- rep_len(O2_vec, n_out)
    N_vec <- rep_len(N_vec, n_out)
    G_vec <- if (is.null(G_input)) O2_vec else rep_len(G_input, n_out)
    same_resource_sweep <- is.null(G_input) || isTRUE(all.equal(
      as.numeric(G_vec),
      as.numeric(O2_vec),
      tolerance = 1e-12,
      check.attributes = FALSE
    ))
    proliferation_rate <- as.numeric(.lambda_eff_of_O2(
      O2 = O2_vec,
      run_params = run_params,
      N = N_vec,
      O2_crit = O2_crit_use,
      O2_growth = o2_growth_use,
      G = G_vec
    ))
    death_rate <- as.numeric(.mu_eff_of_O2(
      O2 = O2_vec,
      run_params = run_params,
      N = N_vec,
      O2_crit = O2_crit_use,
      G = G_vec
    ))
    dead_buffer_rate <- rep(NA_real_, n_out)
    nullisomy_nonviable_rate <- rep(NA_real_, n_out)
    boundary_dropped_rate <- rep(NA_real_, n_out)
    nullisomy_nonviable_division_prob <- rep(NA_real_, n_out)
    nullisomy_nonviable_daughters_per_division <- rep(NA_real_, n_out)
    nullisomy_nonviable_daughter_fraction <- rep(NA_real_, n_out)
    if (isTRUE(same_resource_sweep)) {
      row_groups <- split(
        seq_len(n_out),
        sprintf("%.6f", O2_vec)
      )
      for (row_idx in row_groups) {
        rate_curves <- get_functional_rate_curves(
          O2 = O2_vec[[row_idx[[1]]]]
        )
        state_idx <- as.integer(round(N_vec[row_idx])) - as.integer(cfg$N_MIN) + 1L
        valid_idx <- is.finite(state_idx) & state_idx >= 1L & state_idx <= length(rate_curves$dead_buffer_rate)
        if (any(valid_idx)) {
          idx_use <- state_idx[valid_idx]
          row_use <- row_idx[valid_idx]
          dead_buffer_rate[row_use] <- rate_curves$dead_buffer_rate[idx_use]
          nullisomy_nonviable_rate[row_use] <- rate_curves$nullisomy_nonviable_rate[idx_use]
          boundary_dropped_rate[row_use] <- rate_curves$boundary_dropped_rate[idx_use]
          nullisomy_nonviable_division_prob[row_use] <- rate_curves$nullisomy_nonviable_division_prob[idx_use]
          nullisomy_nonviable_daughters_per_division[row_use] <- rate_curves$nullisomy_nonviable_daughters_per_division[idx_use]
          nullisomy_nonviable_daughter_fraction[row_use] <- rate_curves$nullisomy_nonviable_daughter_fraction[idx_use]
        }
      }
    }
    data.frame(
      O2 = O2_vec,
      G = G_vec,
      N = N_vec,
      proliferation_rate = pmax(as.numeric(proliferation_rate), 0),
      death_rate = pmax(as.numeric(death_rate), 0),
      buffer_death_rate = pmax(as.numeric(dead_buffer_rate), 0),
      buffer_death_per_division = pmax(as.numeric(dead_buffer_rate), 0) / pmax(as.numeric(proliferation_rate), 1e-12),
      nullisomy_nonviable_rate = pmax(as.numeric(nullisomy_nonviable_rate), 0),
      boundary_dropped_rate = pmax(as.numeric(boundary_dropped_rate), 0),
      nullisomy_nonviable_division_prob = pmax(pmin(as.numeric(nullisomy_nonviable_division_prob), 1), 0),
      nullisomy_nonviable_daughters_per_division = pmax(pmin(as.numeric(nullisomy_nonviable_daughters_per_division), 2), 0),
      nullisomy_nonviable_daughter_fraction = pmax(pmin(as.numeric(nullisomy_nonviable_daughter_fraction), 1), 0),
      net_growth_rate = as.numeric(proliferation_rate) - as.numeric(death_rate),
      stringsAsFactors = FALSE
    )
  }

  o2_curve <- dplyr::bind_rows(lapply(seq_len(nrow(ref_df)), function(i) {
    cohort_i <- ref_df$cohort[[i]]
    N_ref <- ref_df$N_ref[[i]]
    ms_rate <- as.numeric(.pmisseg_of_O2(
      O2 = o2_grid,
      run_params = run_params,
      N = N_ref,
      O2_crit = O2_crit_use
    ))
    rate_df <- compute_rate_components(
      O2 = o2_grid,
      N = rep(N_ref, length(o2_grid))
    )
    data.frame(
      oxygen_pct = o2_grid,
      cohort = cohort_i,
      ms_rate = pmax(ms_rate, 0),
      proliferation_rate = rate_df$proliferation_rate,
      death_rate = rate_df$death_rate,
      buffer_death_rate = rate_df$buffer_death_rate,
      buffer_death_per_division = rate_df$buffer_death_per_division,
      nullisomy_nonviable_rate = rate_df$nullisomy_nonviable_rate,
      boundary_dropped_rate = rate_df$boundary_dropped_rate,
      nullisomy_nonviable_division_prob = rate_df$nullisomy_nonviable_division_prob,
      nullisomy_nonviable_daughters_per_division = rate_df$nullisomy_nonviable_daughters_per_division,
      nullisomy_nonviable_daughter_fraction = rate_df$nullisomy_nonviable_daughter_fraction,
      net_growth_rate = rate_df$net_growth_rate,
      N_ref = N_ref,
      row.names = NULL
    )
  }))
  write.table(
    o2_curve,
    file = file.path(out_dir, "functional_curve_oxygen.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  o2_curve_multi <- dplyr::bind_rows(lapply(seq_len(nrow(ref_df_multi)), function(i) {
    cohort_i <- ref_df_multi$cohort[[i]]
    N_ref <- ref_df_multi$N_ref[[i]]
    ploidy_multiple_i <- ref_df_multi$ploidy_multiple[[i]]
    ms_rate <- as.numeric(.pmisseg_of_O2(
      O2 = o2_grid,
      run_params = run_params,
      N = N_ref,
      O2_crit = O2_crit_use
    ))
    rate_df <- compute_rate_components(
      O2 = o2_grid,
      N = rep(N_ref, length(o2_grid))
    )
    data.frame(
      oxygen_pct = o2_grid,
      cohort = cohort_i,
      ploidy_multiple = ploidy_multiple_i,
      N_ref = N_ref,
      ms_rate = pmax(ms_rate, 0),
      proliferation_rate = rate_df$proliferation_rate,
      death_rate = rate_df$death_rate,
      buffer_death_rate = rate_df$buffer_death_rate,
      buffer_death_per_division = rate_df$buffer_death_per_division,
      nullisomy_nonviable_rate = rate_df$nullisomy_nonviable_rate,
      boundary_dropped_rate = rate_df$boundary_dropped_rate,
      nullisomy_nonviable_division_prob = rate_df$nullisomy_nonviable_division_prob,
      nullisomy_nonviable_daughters_per_division = rate_df$nullisomy_nonviable_daughters_per_division,
      nullisomy_nonviable_daughter_fraction = rate_df$nullisomy_nonviable_daughter_fraction,
      net_growth_rate = rate_df$net_growth_rate,
      row.names = NULL
    )
  }))
  write.table(
    o2_curve_multi,
    file = file.path(out_dir, "functional_curve_oxygen_multi_ploidy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  unlink(file.path(out_dir, c(
    "functional_curve_oxygen_gmin.tsv",
    "functional_curve_oxygen_multi_ploidy_gmin.tsv",
    "functional_curve_ploidy_by_o2_gmin.tsv",
    "functional_curve_oxygen_gmax.tsv",
    "functional_curve_oxygen_multi_ploidy_gmax.tsv",
    "functional_curve_ploidy_by_o2_gmax.tsv",
    "functional_curve_oxygen_g18.tsv",
    "functional_curve_oxygen_multi_ploidy_g18.tsv",
    "functional_curve_ploidy_by_o2_g18.tsv",
    "functional_curve_oxygen_g20.tsv",
    "functional_curve_oxygen_multi_ploidy_g20.tsv",
    "functional_curve_ploidy_by_o2_g20.tsv",
    "functional_curve_o2_g_heatmap.tsv",
    "functional_curve_o2_g_heatmap_g0_5.tsv",
    "oxygen_vs_missegregation_rate_gmin.pdf",
    "oxygen_vs_missegregation_rate_multi_ploidy_gmin.pdf",
    "oxygen_vs_proliferation_rate_gmin.pdf",
    "oxygen_vs_death_rate_gmin.pdf",
    "ploidy_vs_proliferation_rate_by_o2_gmin.pdf",
    "ploidy_vs_death_rate_by_o2_gmin.pdf",
    "oxygen_vs_missegregation_rate_gmax.pdf",
    "oxygen_vs_missegregation_rate_multi_ploidy_gmax.pdf",
    "oxygen_vs_proliferation_rate_gmax.pdf",
    "oxygen_vs_death_rate_gmax.pdf",
    "ploidy_vs_proliferation_rate_by_o2_gmax.pdf",
    "ploidy_vs_death_rate_by_o2_gmax.pdf",
    "oxygen_vs_missegregation_rate_g18.pdf",
    "oxygen_vs_missegregation_rate_multi_ploidy_g18.pdf",
    "oxygen_vs_proliferation_rate_g18.pdf",
    "oxygen_vs_death_rate_g18.pdf",
    "ploidy_vs_proliferation_rate_by_o2_g18.pdf",
    "ploidy_vs_death_rate_by_o2_g18.pdf",
    "oxygen_vs_missegregation_rate_g20.pdf",
    "oxygen_vs_missegregation_rate_multi_ploidy_g20.pdf",
    "oxygen_vs_proliferation_rate_g20.pdf",
    "oxygen_vs_death_rate_g20.pdf",
    "ploidy_vs_proliferation_rate_by_o2_g20.pdf",
    "ploidy_vs_death_rate_by_o2_g20.pdf",
    "compare_ploidy_vs_death_rate_by_o2_g20.pdf",
    "compare_ploidy_vs_proliferation_rate_by_o2_g20.pdf",
    "compare_oxygen_vs_missegregation_rate_multi_ploidy_g20.pdf",
    "compare_oxygen_vs_proliferation_rate_g20.pdf",
    "compare_oxygen_vs_death_rate_g20.pdf",
    "o2_g_heatmap_growth.pdf",
    "o2_g_heatmap_death.pdf",
    "o2_g_heatmap_missegregation.pdf",
    "o2_g_heatmap_growth_g0_5.pdf",
    "o2_g_heatmap_death_g0_5.pdf",
    "o2_g_heatmap_missegregation_g0_5.pdf"
  )), force = TRUE)
  fixed_glucose_results <- list()
  if (isTRUE(has_fixed_glucose)) {
    for (spec_i in seq_len(nrow(fixed_glucose_specs))) {
      spec <- fixed_glucose_specs[spec_i, , drop = FALSE]
      suffix <- as.character(spec$suffix[[1]])
      fixed_glucose_pct_i <- as.numeric(spec$value[[1]])
      fixed_glucose_label_i <- paste0(
        "G fixed at ",
        as.character(spec$name[[1]]),
        " = ",
        signif(fixed_glucose_pct_i, 4),
        "%; O2 varies"
      )
      o2_curve_fixed <- dplyr::bind_rows(lapply(seq_len(nrow(ref_df)), function(i) {
        cohort_i <- ref_df$cohort[[i]]
        N_ref <- ref_df$N_ref[[i]]
        G_vec <- rep(fixed_glucose_pct_i, length(o2_grid))
        ms_rate <- as.numeric(.pmisseg_of_O2(
          O2 = o2_grid,
          G = G_vec,
          run_params = run_params,
          N = N_ref,
          O2_crit = O2_crit_use
        ))
        rate_df <- compute_rate_components(
          O2 = o2_grid,
          G = G_vec,
          N = rep(N_ref, length(o2_grid))
        )
        data.frame(
          oxygen_pct = o2_grid,
          glucose_pct = G_vec,
          cohort = cohort_i,
          ms_rate = pmax(ms_rate, 0),
          proliferation_rate = rate_df$proliferation_rate,
          death_rate = rate_df$death_rate,
          net_growth_rate = rate_df$net_growth_rate,
          N_ref = N_ref,
          row.names = NULL
        )
      }))
      write.table(
        o2_curve_fixed,
        file = file.path(out_dir, paste0("functional_curve_oxygen_", suffix, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )

      o2_curve_multi_fixed <- dplyr::bind_rows(lapply(seq_len(nrow(ref_df_multi)), function(i) {
        cohort_i <- ref_df_multi$cohort[[i]]
        N_ref <- ref_df_multi$N_ref[[i]]
        ploidy_multiple_i <- ref_df_multi$ploidy_multiple[[i]]
        G_vec <- rep(fixed_glucose_pct_i, length(o2_grid))
        ms_rate <- as.numeric(.pmisseg_of_O2(
          O2 = o2_grid,
          G = G_vec,
          run_params = run_params,
          N = N_ref,
          O2_crit = O2_crit_use
        ))
        rate_df <- compute_rate_components(
          O2 = o2_grid,
          G = G_vec,
          N = rep(N_ref, length(o2_grid))
        )
        data.frame(
          oxygen_pct = o2_grid,
          glucose_pct = G_vec,
          cohort = cohort_i,
          ploidy_multiple = ploidy_multiple_i,
          N_ref = N_ref,
          ms_rate = pmax(ms_rate, 0),
          proliferation_rate = rate_df$proliferation_rate,
          death_rate = rate_df$death_rate,
          net_growth_rate = rate_df$net_growth_rate,
          row.names = NULL
        )
      }))
      write.table(
        o2_curve_multi_fixed,
        file = file.path(out_dir, paste0("functional_curve_oxygen_multi_ploidy_", suffix, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      fixed_glucose_results[[suffix]] <- list(
        spec = spec,
        label = fixed_glucose_label_i,
        o2_curve = o2_curve_fixed,
        o2_curve_multi = o2_curve_multi_fixed
      )
    }
  }
  o2_g_heatmap <- data.frame()
  o2_g_heatmap_g0_5 <- data.frame()
  if (isTRUE(glucose_use)) {
    compute_o2_g_heatmap <- function(glucose_grid) {
      paired_grid <- expand.grid(
        oxygen_pct = o2_grid,
        glucose_pct = glucose_grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      out <- dplyr::bind_rows(lapply(seq_len(nrow(ref_df)), function(i) {
        cohort_i <- ref_df$cohort[[i]]
        N_ref <- ref_df$N_ref[[i]]
        ms_rate <- as.numeric(.pmisseg_of_O2(
          O2 = paired_grid$oxygen_pct,
          G = paired_grid$glucose_pct,
          run_params = run_params,
          N = N_ref,
          O2_crit = O2_crit_use
        ))
        rate_df <- compute_rate_components(
          O2 = paired_grid$oxygen_pct,
          G = paired_grid$glucose_pct,
          N = rep(N_ref, nrow(paired_grid))
        )
        data.frame(
          oxygen_pct = paired_grid$oxygen_pct,
          glucose_pct = paired_grid$glucose_pct,
          cohort = cohort_i,
          N_ref = N_ref,
          ms_rate = pmax(ms_rate, 0),
          proliferation_rate = rate_df$proliferation_rate,
          death_rate = rate_df$death_rate,
          net_growth_rate = rate_df$net_growth_rate,
          row.names = NULL
        )
      }))
      out$cohort <- factor(out$cohort, levels = ref_df$cohort)
      out
    }
    o2_g_heatmap <- compute_o2_g_heatmap(glucose_heatmap_grid)
    o2_g_heatmap$cohort <- factor(o2_g_heatmap$cohort, levels = ref_df$cohort)
    write.table(
      o2_g_heatmap,
      file = file.path(out_dir, "functional_curve_o2_g_heatmap.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    o2_g_heatmap_g0_5 <- compute_o2_g_heatmap(glucose_heatmap_o2range_grid)
    write.table(
      o2_g_heatmap_g0_5,
      file = file.path(out_dir, "functional_curve_o2_g_heatmap_g0_5.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  multi_colors <- stats::setNames(
    grDevices::hcl.colors(nrow(ref_df_multi), palette = "Dark 3"),
    ref_df_multi$cohort
  )

  unlink(file.path(out_dir, c(
    "functional_curve_glucose.tsv",
    "functional_curve_glucose_multi_ploidy.tsv",
    "glucose_vs_missegregation_rate.pdf",
    "glucose_vs_missegregation_rate_multi_ploidy.pdf",
    "glucose_vs_proliferation_rate.pdf",
    "glucose_vs_death_rate.pdf",
    "glucose_vs_net_growth_rate.pdf"
  )), force = TRUE)

  p_msr_o2 <- ggplot(o2_curve, aes(x = oxygen_pct, y = ms_rate, color = cohort)) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = "Oxygen vs Missegregation Rate",
      subtitle = "MS rate = p_mis_base + p_misseg * mu_eff / (mu_eff + k_o_mis), clamped to [0,1]",
      x = "Oxygen (%)",
      y = "MS rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  p_msr_o2_multi <- ggplot(
    o2_curve_multi,
    aes(x = oxygen_pct, y = ms_rate, color = factor(cohort, levels = ref_df_multi$cohort))
  ) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
    scale_color_manual(values = multi_colors, drop = FALSE) +
    labs(
      title = "Oxygen vs Missegregation Rate Across Reference Ploidy States",
      subtitle = "Reference states: 1.5N, 2N, 2.5N, 3N, 3.5N, 4N, 4.5N, 5N",
      x = "Oxygen (%)",
      y = "MS rate",
      color = "Reference state"
    ) +
    theme_bw(base_size = 11)

  p_msr_death <- ggplot(
    o2_curve_multi,
    aes(x = ms_rate, y = death_rate, color = factor(cohort, levels = ref_df_multi$cohort))
  ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = multi_colors, drop = FALSE) +
    labs(
      title = "Death Rate vs MS Rate Across Reference Ploidy States",
      subtitle = ref_state_subtitle,
      x = "MS rate",
      y = "Death rate",
      color = "Reference state"
    ) +
    theme_bw(base_size = 11)
  p_death_msr <- ggplot(
    o2_curve,
    aes(x = death_rate, y = ms_rate, color = cohort, group = cohort)
  ) +
    geom_path(linewidth = 1) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = "Death Rate vs Missegregation Rate",
      subtitle = "Same oxygen sweep and reference ploidy states as Oxygen vs Missegregation Rate",
      x = "Death rate",
      y = "MS rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  p_msr_buffer_death <- ggplot(
    o2_curve,
    aes(
      x = ms_rate,
      y = buffer_death_rate,
      group = cohort,
      color = oxygen_pct
    )
  ) +
    geom_path(linewidth = 1.1, lineend = "round") +
    facet_wrap(~ cohort, nrow = 1) +
    scale_color_gradient(
      low = "#2C7BB6",
      high = "#F28E2B",
      limits = c(o2_plot_min, o2_plot_max),
      name = "O2 (%)"
    ) +
    labs(
      title = "Buffer-Death Rate vs MS Rate",
      subtitle = "Total dead-buffer inflow = missegregation-linked nonviability + boundary-drop losses",
      x = "MS rate",
      y = "Buffer-death rate"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )
  p_msr_buffer_death_per_division <- ggplot(
    o2_curve_multi,
    aes(
      x = ms_rate,
      y = nullisomy_nonviable_daughters_per_division,
      color = factor(cohort, levels = ref_df_multi$cohort)
    )
  ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = multi_colors, drop = FALSE) +
    labs(
      title = "Expected Nonviable Daughters per Division vs MS Rate",
      subtitle = paste0(
        "Missegregation-linked expected dead daughters per division; excludes boundary-drop losses | ",
        ref_state_subtitle
      ),
      x = "MS rate",
      y = "Expected nonviable daughters / division",
      color = "Reference state"
    ) +
    theme_bw(base_size = 11)
  p_msr_nonviable_daughter_fraction <- ggplot(
    o2_curve_multi,
    aes(
      x = ms_rate,
      y = nullisomy_nonviable_daughter_fraction,
      color = factor(cohort, levels = ref_df_multi$cohort)
    )
  ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = multi_colors, drop = FALSE) +
    labs(
      title = "Nonviable Daughter Fraction vs MS Rate Across Reference Ploidy States",
      subtitle = paste0(
        "Missegregation-linked nonviable daughters / all daughters per division; excludes boundary-drop losses | ",
        ref_state_subtitle
      ),
      x = "MS rate",
      y = "Nonviable daughters / all daughters",
      color = "Reference state"
    ) +
    theme_bw(base_size = 11)
  p_msr_nonviable_division_prob <- ggplot(
    o2_curve_multi,
    aes(
      x = ms_rate,
      y = nullisomy_nonviable_division_prob,
      color = factor(cohort, levels = ref_df_multi$cohort)
    )
  ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = multi_colors, drop = FALSE) +
    labs(
      title = "Capped Nonviable Daughter Burden vs MS Rate Across Reference Ploidy States",
      subtitle = paste0(
        "Legacy capped per-division burden; excludes boundary-drop losses | ",
        ref_state_subtitle
      ),
      x = "MS rate",
      y = "min(expected nonviable daughters / division, 1)",
      color = "Reference state"
    ) +
    theme_bw(base_size = 11)

  p_prolif <- ggplot(
    o2_curve_multi,
    aes(x = oxygen_pct, y = proliferation_rate, color = factor(cohort, levels = ref_df_multi$cohort))
  ) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
    scale_color_manual(values = multi_colors, drop = FALSE) +
    labs(
      title = "Oxygen vs Proliferation Rate Across Reference Ploidy States",
      subtitle = ref_state_subtitle,
      x = "Oxygen (%)",
      y = "Proliferation rate",
      color = "Reference state"
    ) +
    theme_bw(base_size = 11)
  p_death <- ggplot(
    o2_curve_multi,
    aes(x = oxygen_pct, y = death_rate, color = factor(cohort, levels = ref_df_multi$cohort))
  ) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
    scale_color_manual(values = multi_colors, drop = FALSE) +
    labs(
      title = "Oxygen vs Death Rate Across Reference Ploidy States",
      subtitle = ref_state_subtitle,
      x = "Oxygen (%)",
      y = "Death rate",
      color = "Reference state"
    ) +
    theme_bw(base_size = 11)
  if (nrow(o2_g_heatmap) > 0L || nrow(o2_g_heatmap_g0_5) > 0L) {
    make_o2_g_heatmap <- function(data, value_col, title, fill_label, y_breaks, subtitle_suffix) {
      plot_df <- data
      plot_df$value <- as.numeric(plot_df[[value_col]])
      ggplot(plot_df, aes(x = oxygen_pct, y = glucose_pct, fill = value)) +
        geom_raster(interpolate = FALSE) +
        facet_grid(. ~ cohort) +
        scale_x_continuous(
          breaks = seq(o2_plot_min, o2_plot_max, by = 1),
          expand = c(0, 0)
        ) +
        scale_y_continuous(
          breaks = y_breaks,
          expand = c(0, 0)
        ) +
        scale_fill_gradientn(
          colors = grDevices::hcl.colors(256, palette = "Viridis"),
          name = fill_label,
          na.value = "grey90"
        ) +
        labs(
          title = title,
          subtitle = paste0("Paired O2-G grid; panels show 2N and 4N reference states | ", subtitle_suffix),
          x = "Oxygen (%)",
          y = "Glucose (%)"
        ) +
        theme_bw(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95", color = "grey80")
        )
    }
    if (nrow(o2_g_heatmap) > 0L) {
      p_o2_g_growth <- make_o2_g_heatmap(
        o2_g_heatmap,
        "proliferation_rate",
        "O2-G Heatmap: Growth",
        "Growth",
        y_breaks = seq(0, 100, by = 20),
        subtitle_suffix = "G range 0-100%"
      )
      p_o2_g_death <- make_o2_g_heatmap(
        o2_g_heatmap,
        "death_rate",
        "O2-G Heatmap: Death",
        "Death",
        y_breaks = seq(0, 100, by = 20),
        subtitle_suffix = "G range 0-100%"
      )
      p_o2_g_ms <- make_o2_g_heatmap(
        o2_g_heatmap,
        "ms_rate",
        "O2-G Heatmap: MS Rate",
        "MS rate",
        y_breaks = seq(0, 100, by = 20),
        subtitle_suffix = "G range 0-100%"
      )
    }
    if (nrow(o2_g_heatmap_g0_5) > 0L) {
      p_o2_g_growth_g0_5 <- make_o2_g_heatmap(
        o2_g_heatmap_g0_5,
        "proliferation_rate",
        "O2-G Heatmap: Growth (G Range 0-5)",
        "Growth",
        y_breaks = seq(o2_plot_min, o2_plot_max, by = 1),
        subtitle_suffix = "G range 0-5%, matching O2"
      )
      p_o2_g_death_g0_5 <- make_o2_g_heatmap(
        o2_g_heatmap_g0_5,
        "death_rate",
        "O2-G Heatmap: Death (G Range 0-5)",
        "Death",
        y_breaks = seq(o2_plot_min, o2_plot_max, by = 1),
        subtitle_suffix = "G range 0-5%, matching O2"
      )
      p_o2_g_ms_g0_5 <- make_o2_g_heatmap(
        o2_g_heatmap_g0_5,
        "ms_rate",
        "O2-G Heatmap: MS Rate (G Range 0-5)",
        "MS rate",
        y_breaks = seq(o2_plot_min, o2_plot_max, by = 1),
        subtitle_suffix = "G range 0-5%, matching O2"
      )
    }
  }
  p_net <- ggplot(o2_curve, aes(x = oxygen_pct, y = net_growth_rate)) +
    geom_line(linewidth = 1, color = "#2ca02c") +
    facet_wrap(~cohort, ncol = 2) +
    coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
    labs(
      title = "Oxygen vs Net Growth Rate",
      subtitle = paste0("Net rate = proliferation - ", death_language$adjective, " high-ploidy death"),
      x = "Oxygen (%)",
      y = "Net growth rate"
    ) +
    theme_bw(base_size = 11)

  N_states <- seq.int(as.integer(cfg$N_MIN), as.integer(cfg$N_MAX))
  ploidy_grid <- N_states / as.numeric(cfg$N_UNIT)
  gamma_loss_ref <- as.numeric(.first_non_null_local(run_params$gamma_loss, 0.1))
  if (!is.finite(gamma_loss_ref) || gamma_loss_ref <= 0) gamma_loss_ref <- 0.1
  viability <- .loss_survival_nullisomy(
    N_states,
    m_loss = 1L,
    gamma_loss = gamma_loss_ref,
    N_unit = as.integer(cfg$N_UNIT)
  )
  viability_curve <- data.frame(
    N = N_states,
    ploidy = ploidy_grid,
    endpoint_value = if (identical(start_with_mode, "chr_number")) as.numeric(N_states) else ploidy_grid,
    viability_after_ms = pmax(viability, 0),
    row.names = NULL
  )
  write.table(
    viability_curve,
    file = file.path(out_dir, "functional_curve_ploidy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  ploidy_o2_curve <- dplyr::bind_rows(lapply(o2_levels_ploidy, function(o2_level) {
    N_grid <- if (identical(start_with_mode, "chr_number")) {
      as.numeric(state_grid_dense)
    } else {
      as.numeric(state_grid_dense * as.numeric(cfg$N_UNIT))
    }
    rate_df <- compute_rate_components(O2 = rep(o2_level, length(N_grid)), N = N_grid)
    data.frame(
      oxygen_pct = rep(as.numeric(o2_level), length(N_grid)),
      ploidy = if (identical(start_with_mode, "chr_number")) as.numeric(N_grid) else as.numeric(state_grid_dense),
      N = N_grid,
      endpoint_value = if (identical(start_with_mode, "chr_number")) as.numeric(N_grid) else as.numeric(state_grid_dense),
      proliferation_rate = rate_df$proliferation_rate,
      death_rate = rate_df$death_rate,
      net_growth_rate = rate_df$net_growth_rate,
      row.names = NULL
    )
  }))
  write.table(
    ploidy_o2_curve,
    file = file.path(out_dir, "functional_curve_ploidy_by_o2.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  if (length(fixed_glucose_results) > 0L) {
    for (suffix in names(fixed_glucose_results)) {
      fixed_glucose_pct_i <- as.numeric(fixed_glucose_results[[suffix]]$spec$value[[1]])
      ploidy_o2_curve_fixed <- dplyr::bind_rows(lapply(o2_levels_ploidy, function(o2_level) {
        N_grid <- if (identical(start_with_mode, "chr_number")) {
          as.numeric(state_grid_dense)
        } else {
          as.numeric(state_grid_dense * as.numeric(cfg$N_UNIT))
        }
        rate_df <- compute_rate_components(
          O2 = rep(o2_level, length(N_grid)),
          G = rep(fixed_glucose_pct_i, length(N_grid)),
          N = N_grid
        )
        data.frame(
          oxygen_pct = rep(as.numeric(o2_level), length(N_grid)),
          glucose_pct = rep(fixed_glucose_pct_i, length(N_grid)),
          ploidy = if (identical(start_with_mode, "chr_number")) as.numeric(N_grid) else as.numeric(state_grid_dense),
          N = N_grid,
          endpoint_value = if (identical(start_with_mode, "chr_number")) as.numeric(N_grid) else as.numeric(state_grid_dense),
          proliferation_rate = rate_df$proliferation_rate,
          death_rate = rate_df$death_rate,
          net_growth_rate = rate_df$net_growth_rate,
          row.names = NULL
        )
      }))
      write.table(
        ploidy_o2_curve_fixed,
        file = file.path(out_dir, paste0("functional_curve_ploidy_by_o2_", suffix, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      fixed_glucose_results[[suffix]]$ploidy_o2_curve <- ploidy_o2_curve_fixed
    }
  }

  p_viability <- ggplot(viability_curve, aes(x = endpoint_value, y = viability_after_ms)) +
    geom_line(color = "#2ca02c", linewidth = 1) +
    labs(
      title = paste0(state_axis_label, " vs Viability After MS"),
      subtitle = "Nullisomy-risk loss survival for a one-copy loss event",
      x = state_axis_label,
      y = "Viability after MS"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_prolif_o2 <- ggplot(
    ploidy_o2_curve,
    aes(x = endpoint_value, y = proliferation_rate, color = oxygen_pct)
  ) +
    geom_point(shape = 15, size = 1.8, alpha = 0.95) +
    coord_cartesian(xlim = c(state_plot_min, state_plot_max)) +
    scale_color_gradient(
      low = "#2C7BB6",
      high = "#F28E2B",
      limits = c(o2_plot_min, o2_plot_max),
      name = "O2 level"
    ) +
    labs(
      title = paste0(state_axis_label, " vs Proliferation Rate Colored by O2"),
      subtitle = paste0("Dense square markers over ", tolower(state_axis_label), " range, colored by O2 level 0-5"),
      x = state_axis_label,
      y = "Proliferation rate"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_death_o2 <- ggplot(
    ploidy_o2_curve,
    aes(x = endpoint_value, y = death_rate, color = oxygen_pct)
  ) +
    geom_point(shape = 15, size = 1.8, alpha = 0.95) +
    coord_cartesian(xlim = c(state_plot_min, state_plot_max)) +
    scale_color_gradient(
      low = "#2C7BB6",
      high = "#F28E2B",
      limits = c(o2_plot_min, o2_plot_max),
      name = "O2 level"
    ) +
    labs(
      title = paste0(state_axis_label, " vs Death Rate Colored by O2"),
      subtitle = paste0("Dense square markers over ", tolower(state_axis_label), " range, colored by O2 level 0-5"),
      x = state_axis_label,
      y = "Death rate"
    ) +
    theme_bw(base_size = 11)
  if (length(fixed_glucose_results) > 0L) {
    for (suffix in names(fixed_glucose_results)) {
      title_label_i <- as.character(fixed_glucose_results[[suffix]]$spec$title_label[[1]])
      fixed_glucose_label_i <- fixed_glucose_results[[suffix]]$label
      o2_curve_fixed <- fixed_glucose_results[[suffix]]$o2_curve
      o2_curve_multi_fixed <- fixed_glucose_results[[suffix]]$o2_curve_multi
      ploidy_o2_curve_fixed <- fixed_glucose_results[[suffix]]$ploidy_o2_curve
      fixed_glucose_results[[suffix]]$plots <- list(
        ms_o2 = ggplot(o2_curve_fixed, aes(x = oxygen_pct, y = ms_rate, color = cohort)) +
          geom_line(linewidth = 1) +
          coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
          scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
          labs(
            title = paste0("Oxygen vs Missegregation Rate (", title_label_i, ")"),
            subtitle = fixed_glucose_label_i,
            x = "Oxygen (%)",
            y = "MS rate",
            color = "Cohort"
          ) +
          theme_bw(base_size = 11),
        ms_o2_multi = ggplot(
          o2_curve_multi_fixed,
          aes(x = oxygen_pct, y = ms_rate, color = factor(cohort, levels = ref_df_multi$cohort))
        ) +
          geom_line(linewidth = 1) +
          coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
          scale_color_manual(values = multi_colors, drop = FALSE) +
          labs(
            title = paste0("Oxygen vs Missegregation Rate Across Reference Ploidy States (", title_label_i, ")"),
            subtitle = paste0(ref_state_subtitle, " | ", fixed_glucose_label_i),
            x = "Oxygen (%)",
            y = "MS rate",
            color = "Reference state"
          ) +
          theme_bw(base_size = 11),
        proliferation = ggplot(
          o2_curve_multi_fixed,
          aes(x = oxygen_pct, y = proliferation_rate, color = factor(cohort, levels = ref_df_multi$cohort))
        ) +
          geom_line(linewidth = 1) +
          coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
          scale_color_manual(values = multi_colors, drop = FALSE) +
          labs(
            title = paste0("Oxygen vs Proliferation Rate Across Reference Ploidy States (", title_label_i, ")"),
            subtitle = paste0(ref_state_subtitle, " | ", fixed_glucose_label_i),
            x = "Oxygen (%)",
            y = "Proliferation rate",
            color = "Reference state"
          ) +
          theme_bw(base_size = 11),
        death = ggplot(
          o2_curve_multi_fixed,
          aes(x = oxygen_pct, y = death_rate, color = factor(cohort, levels = ref_df_multi$cohort))
        ) +
          geom_line(linewidth = 1) +
          coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
          scale_color_manual(values = multi_colors, drop = FALSE) +
          labs(
            title = paste0("Oxygen vs Death Rate Across Reference Ploidy States (", title_label_i, ")"),
            subtitle = paste0(ref_state_subtitle, " | ", fixed_glucose_label_i),
            x = "Oxygen (%)",
            y = "Death rate",
            color = "Reference state"
          ) +
          theme_bw(base_size = 11),
        ploidy_proliferation = ggplot(
          ploidy_o2_curve_fixed,
          aes(x = endpoint_value, y = proliferation_rate, color = oxygen_pct)
        ) +
          geom_point(shape = 15, size = 1.8, alpha = 0.95) +
          coord_cartesian(xlim = c(state_plot_min, state_plot_max)) +
          scale_color_gradient(
            low = "#2C7BB6",
            high = "#F28E2B",
            limits = c(o2_plot_min, o2_plot_max),
            name = "O2 level"
          ) +
          labs(
            title = paste0(state_axis_label, " vs Proliferation Rate Colored by O2 (", title_label_i, ")"),
            subtitle = paste0(fixed_glucose_label_i, " | Dense square markers over ", tolower(state_axis_label), " range"),
            x = state_axis_label,
            y = "Proliferation rate"
          ) +
          theme_bw(base_size = 11),
        ploidy_death = ggplot(
          ploidy_o2_curve_fixed,
          aes(x = endpoint_value, y = death_rate, color = oxygen_pct)
        ) +
          geom_point(shape = 15, size = 1.8, alpha = 0.95) +
          coord_cartesian(xlim = c(state_plot_min, state_plot_max)) +
          scale_color_gradient(
            low = "#2C7BB6",
            high = "#F28E2B",
            limits = c(o2_plot_min, o2_plot_max),
            name = "O2 level"
          ) +
          labs(
            title = paste0(state_axis_label, " vs Death Rate Colored by O2 (", title_label_i, ")"),
            subtitle = paste0(fixed_glucose_label_i, " | Dense square markers over ", tolower(state_axis_label), " range"),
            x = state_axis_label,
            y = "Death rate"
          ) +
          theme_bw(base_size = 11)
      )
    }
  }
  comparison_plots_g20 <- list()
  if ("g20" %in% names(fixed_glucose_results)) {
    condition_levels <- c("Coupled G=O2", "Fixed G=20")
    coupled_ploidy_o2_curve <- ploidy_o2_curve
    coupled_ploidy_o2_curve$resource_condition <- condition_levels[[1]]
    fixed20_ploidy_o2_curve <- fixed_glucose_results[["g20"]]$ploidy_o2_curve
    fixed20_ploidy_o2_curve$resource_condition <- condition_levels[[2]]
    compare_ploidy_o2_curve <- dplyr::bind_rows(coupled_ploidy_o2_curve, fixed20_ploidy_o2_curve)
    compare_ploidy_o2_curve$resource_condition <- factor(
      compare_ploidy_o2_curve$resource_condition,
      levels = condition_levels
    )

    coupled_o2_curve_multi <- o2_curve_multi
    coupled_o2_curve_multi$resource_condition <- condition_levels[[1]]
    fixed20_o2_curve_multi <- fixed_glucose_results[["g20"]]$o2_curve_multi
    fixed20_o2_curve_multi$resource_condition <- condition_levels[[2]]
    compare_o2_curve_multi <- dplyr::bind_rows(coupled_o2_curve_multi, fixed20_o2_curve_multi)
    compare_o2_curve_multi$resource_condition <- factor(
      compare_o2_curve_multi$resource_condition,
      levels = condition_levels
    )

    comparison_plots_g20$ploidy_death <- ggplot(
      compare_ploidy_o2_curve,
      aes(x = endpoint_value, y = death_rate, color = oxygen_pct)
    ) +
      geom_point(shape = 15, size = 1.8, alpha = 0.95) +
      facet_grid(. ~ resource_condition) +
      coord_cartesian(xlim = c(state_plot_min, state_plot_max)) +
      scale_color_gradient(
        low = "#2C7BB6",
        high = "#F28E2B",
        limits = c(o2_plot_min, o2_plot_max),
        name = "O2 level"
      ) +
      labs(
        title = paste0(state_axis_label, " vs Death Rate by O2: Coupled G vs Fixed G=20"),
        subtitle = paste0("Left: original coupled G=O2; right: fixed G=20 while O2 varies"),
        x = state_axis_label,
        y = "Death rate"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey95", color = "grey80"))

    comparison_plots_g20$ploidy_proliferation <- ggplot(
      compare_ploidy_o2_curve,
      aes(x = endpoint_value, y = proliferation_rate, color = oxygen_pct)
    ) +
      geom_point(shape = 15, size = 1.8, alpha = 0.95) +
      facet_grid(. ~ resource_condition) +
      coord_cartesian(xlim = c(state_plot_min, state_plot_max)) +
      scale_color_gradient(
        low = "#2C7BB6",
        high = "#F28E2B",
        limits = c(o2_plot_min, o2_plot_max),
        name = "O2 level"
      ) +
      labs(
        title = paste0(state_axis_label, " vs Proliferation Rate by O2: Coupled G vs Fixed G=20"),
        subtitle = paste0("Left: original coupled G=O2; right: fixed G=20 while O2 varies"),
        x = state_axis_label,
        y = "Proliferation rate"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey95", color = "grey80"))

    comparison_plots_g20$ms_o2_multi <- ggplot(
      compare_o2_curve_multi,
      aes(x = oxygen_pct, y = ms_rate, color = factor(cohort, levels = ref_df_multi$cohort))
    ) +
      geom_line(linewidth = 1) +
      facet_grid(. ~ resource_condition) +
      coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
      scale_color_manual(values = multi_colors, drop = FALSE) +
      labs(
        title = "Oxygen vs Missegregation Rate: Coupled G vs Fixed G=20",
        subtitle = ref_state_subtitle,
        x = "Oxygen (%)",
        y = "MS rate",
        color = "Reference state"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey95", color = "grey80"))

    comparison_plots_g20$proliferation <- ggplot(
      compare_o2_curve_multi,
      aes(x = oxygen_pct, y = proliferation_rate, color = factor(cohort, levels = ref_df_multi$cohort))
    ) +
      geom_line(linewidth = 1) +
      facet_grid(. ~ resource_condition) +
      coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
      scale_color_manual(values = multi_colors, drop = FALSE) +
      labs(
        title = "Oxygen vs Proliferation Rate: Coupled G vs Fixed G=20",
        subtitle = ref_state_subtitle,
        x = "Oxygen (%)",
        y = "Proliferation rate",
        color = "Reference state"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey95", color = "grey80"))

    comparison_plots_g20$death <- ggplot(
      compare_o2_curve_multi,
      aes(x = oxygen_pct, y = death_rate, color = factor(cohort, levels = ref_df_multi$cohort))
    ) +
      geom_line(linewidth = 1) +
      facet_grid(. ~ resource_condition) +
      coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
      scale_color_manual(values = multi_colors, drop = FALSE) +
      labs(
        title = "Oxygen vs Death Rate: Coupled G vs Fixed G=20",
        subtitle = ref_state_subtitle,
        x = "Oxygen (%)",
        y = "Death rate",
        color = "Reference state"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey95", color = "grey80"))
  }

  ggsave(file.path(out_dir, "oxygen_vs_missegregation_rate.pdf"), p_msr_o2, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_missegregation_rate_multi_ploidy.pdf"), p_msr_o2_multi, width = 10, height = 7)
  ggsave(file.path(out_dir, "ms_rate_vs_death_rate.pdf"), p_msr_death, width = 10, height = 7)
  ggsave(file.path(out_dir, "death_rate_vs_missegregation_rate.pdf"), p_death_msr, width = 10, height = 7)
  ggsave(file.path(out_dir, "ms_rate_vs_buffer_death_rate.pdf"), p_msr_buffer_death, width = 10, height = 7)
  ggsave(file.path(out_dir, "ms_rate_vs_buffer_death_per_division.pdf"), p_msr_buffer_death_per_division, width = 10, height = 7)
  ggsave(file.path(out_dir, "ms_rate_vs_nonviable_daughter_fraction.pdf"), p_msr_nonviable_daughter_fraction, width = 10, height = 7)
  ggsave(file.path(out_dir, "ms_rate_vs_nonviable_division_probability.pdf"), p_msr_nonviable_division_prob, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_proliferation_rate.pdf"), p_prolif, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_death_rate.pdf"), p_death, width = 10, height = 7)
  if (nrow(o2_g_heatmap) > 0L) {
    ggsave(file.path(out_dir, "o2_g_heatmap_growth.pdf"), p_o2_g_growth, width = 10, height = 7)
    ggsave(file.path(out_dir, "o2_g_heatmap_death.pdf"), p_o2_g_death, width = 10, height = 7)
    ggsave(file.path(out_dir, "o2_g_heatmap_missegregation.pdf"), p_o2_g_ms, width = 10, height = 7)
  }
  if (nrow(o2_g_heatmap_g0_5) > 0L) {
    ggsave(file.path(out_dir, "o2_g_heatmap_growth_g0_5.pdf"), p_o2_g_growth_g0_5, width = 10, height = 7)
    ggsave(file.path(out_dir, "o2_g_heatmap_death_g0_5.pdf"), p_o2_g_death_g0_5, width = 10, height = 7)
    ggsave(file.path(out_dir, "o2_g_heatmap_missegregation_g0_5.pdf"), p_o2_g_ms_g0_5, width = 10, height = 7)
  }
  ggsave(file.path(out_dir, "oxygen_vs_net_growth_rate.pdf"), p_net, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_viability_after_ms.pdf"), p_viability, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_proliferation_rate_by_o2.pdf"), p_ploidy_prolif_o2, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_death_rate_by_o2.pdf"), p_ploidy_death_o2, width = 10, height = 7)
  if (length(fixed_glucose_results) > 0L) {
    for (suffix in names(fixed_glucose_results)) {
      plots <- fixed_glucose_results[[suffix]]$plots
      ggsave(file.path(out_dir, paste0("oxygen_vs_missegregation_rate_", suffix, ".pdf")), plots$ms_o2, width = 10, height = 7)
      ggsave(file.path(out_dir, paste0("oxygen_vs_missegregation_rate_multi_ploidy_", suffix, ".pdf")), plots$ms_o2_multi, width = 10, height = 7)
      ggsave(file.path(out_dir, paste0("oxygen_vs_proliferation_rate_", suffix, ".pdf")), plots$proliferation, width = 10, height = 7)
      ggsave(file.path(out_dir, paste0("oxygen_vs_death_rate_", suffix, ".pdf")), plots$death, width = 10, height = 7)
      ggsave(file.path(out_dir, paste0("ploidy_vs_proliferation_rate_by_o2_", suffix, ".pdf")), plots$ploidy_proliferation, width = 10, height = 7)
      ggsave(file.path(out_dir, paste0("ploidy_vs_death_rate_by_o2_", suffix, ".pdf")), plots$ploidy_death, width = 10, height = 7)
    }
  }
  if (length(comparison_plots_g20) > 0L) {
    ggsave(file.path(out_dir, "compare_ploidy_vs_death_rate_by_o2_g20.pdf"), comparison_plots_g20$ploidy_death, width = 14, height = 7)
    ggsave(file.path(out_dir, "compare_ploidy_vs_proliferation_rate_by_o2_g20.pdf"), comparison_plots_g20$ploidy_proliferation, width = 14, height = 7)
    ggsave(file.path(out_dir, "compare_oxygen_vs_missegregation_rate_multi_ploidy_g20.pdf"), comparison_plots_g20$ms_o2_multi, width = 14, height = 7)
    ggsave(file.path(out_dir, "compare_oxygen_vs_proliferation_rate_g20.pdf"), comparison_plots_g20$proliferation, width = 14, height = 7)
    ggsave(file.path(out_dir, "compare_oxygen_vs_death_rate_g20.pdf"), comparison_plots_g20$death, width = 14, height = 7)
  }

  invisible(list(
    p_msr_o2 = p_msr_o2,
    p_msr_o2_multi = p_msr_o2_multi,
    p_msr_death = p_msr_death,
    p_death_msr = p_death_msr,
    p_msr_buffer_death = p_msr_buffer_death,
    p_msr_buffer_death_per_division = p_msr_buffer_death_per_division,
    p_msr_nonviable_daughter_fraction = p_msr_nonviable_daughter_fraction,
    p_msr_nonviable_division_prob = p_msr_nonviable_division_prob,
    p_prolif = p_prolif,
    p_death = p_death,
    p_net = p_net,
    p_viability = p_viability,
    p_ploidy_prolif_o2 = p_ploidy_prolif_o2,
    p_ploidy_death_o2 = p_ploidy_death_o2,
    fixed_glucose_results = fixed_glucose_results,
    comparison_plots_g20 = comparison_plots_g20
  ))
}

# -----------------------------------------------------------------------------
# Function: plot_predict_horizon
# Purpose: Generate and save visualization output for fitted model behavior.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - scenarios: List of scenario-specific observation data and metadata.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - out_dir: Output directory for generated files and plots.
#   - horizon_day: Prediction horizon end time in days.
#   - report_dt: Sampling interval for reported trajectories.
#   - top_n: Maximum number of scenarios selected for detailed plotting.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
extract_ggplot_legend_grob <- function(plot) {
  if (!inherits(plot, "ggplot")) return(NULL)
  plot_grob <- ggplot2::ggplotGrob(plot + theme(legend.position = "right"))
  guide_idx <- which(vapply(plot_grob$grobs, function(g) g$name %||% "", character(1)) == "guide-box")
  if (!length(guide_idx)) return(NULL)
  plot_grob$grobs[[guide_idx[[1]]]]
}

save_aligned_plot_stack_with_shared_legend <- function(plots, path, width = 14, height = 15) {
  plots <- Filter(function(p) inherits(p, "ggplot"), plots)
  if (!length(plots)) return(invisible(NULL))

  legend_grob <- extract_ggplot_legend_grob(plots[[1]])
  plot_grobs <- lapply(plots, function(p) ggplot2::ggplotGrob(p + theme(legend.position = "none")))
  width_lengths <- vapply(plot_grobs, function(g) length(g$widths), integer(1))
  if (length(unique(width_lengths)) == 1L) {
    max_widths <- do.call(grid::unit.pmax, lapply(plot_grobs, function(g) g$widths))
    plot_grobs <- lapply(plot_grobs, function(g) {
      g$widths <- max_widths
      g
    })
  }

  grDevices::pdf(path, width = width, height = height, onefile = TRUE)
  grid::grid.newpage()
  if (is.null(legend_grob)) {
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = length(plot_grobs), ncol = 1)))
    for (i in seq_along(plot_grobs)) {
      grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
      grid::grid.draw(plot_grobs[[i]])
      grid::upViewport()
    }
    grid::upViewport()
  } else {
    lay <- grid::grid.layout(
      nrow = length(plot_grobs),
      ncol = 2,
      widths = grid::unit.c(grid::unit(1, "null"), grid::grobWidth(legend_grob) + grid::unit(0.25, "in"))
    )
    grid::pushViewport(grid::viewport(layout = lay))
    for (i in seq_along(plot_grobs)) {
      grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
      grid::grid.draw(plot_grobs[[i]])
      grid::upViewport()
    }
    grid::pushViewport(grid::viewport(layout.pos.row = seq_along(plot_grobs), layout.pos.col = 2))
    grid::grid.draw(legend_grob)
    grid::upViewport()
    grid::upViewport()
  }
  grDevices::dev.off()
  invisible(path)
}

save_aligned_plot_row_with_shared_legend <- function(plots, path, width = 18, height = 6.5) {
  plots <- Filter(function(p) inherits(p, "ggplot"), plots)
  if (!length(plots)) return(invisible(NULL))

  legend_grob <- extract_ggplot_legend_grob(plots[[1]])
  plot_grobs <- lapply(plots, function(p) ggplot2::ggplotGrob(p + theme(legend.position = "none")))
  height_lengths <- vapply(plot_grobs, function(g) length(g$heights), integer(1))
  if (length(unique(height_lengths)) == 1L) {
    max_heights <- do.call(grid::unit.pmax, lapply(plot_grobs, function(g) g$heights))
    plot_grobs <- lapply(plot_grobs, function(g) {
      g$heights <- max_heights
      g
    })
  }

  grDevices::pdf(path, width = width, height = height, onefile = TRUE)
  grid::grid.newpage()
  if (is.null(legend_grob)) {
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = length(plot_grobs))))
    for (i in seq_along(plot_grobs)) {
      grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
      grid::grid.draw(plot_grobs[[i]])
      grid::upViewport()
    }
    grid::upViewport()
  } else {
    lay <- grid::grid.layout(
      nrow = 1,
      ncol = length(plot_grobs) + 1L,
      widths = grid::unit.c(
        rep(grid::unit(1, "null"), length(plot_grobs)),
        grid::grobWidth(legend_grob) + grid::unit(0.25, "in")
      )
    )
    grid::pushViewport(grid::viewport(layout = lay))
    for (i in seq_along(plot_grobs)) {
      grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
      grid::grid.draw(plot_grobs[[i]])
      grid::upViewport()
    }
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = length(plot_grobs) + 1L))
    grid::grid.draw(legend_grob)
    grid::upViewport()
    grid::upViewport()
  }
  grDevices::dev.off()
  invisible(path)
}

plot_predict_horizon <- function(run_params, scenarios, cfg, out_dir, horizon_day, report_dt = 1.0, top_n = 6L) {
  save_aligned_plot_stack <- function(plots, path, width = 14, height = 7, row_heights = NULL) {
    plots <- Filter(function(p) inherits(p, "ggplot"), plots)
    if (!length(plots)) return(invisible(NULL))
    if (is.null(row_heights)) {
      row_heights <- rep(1, length(plots))
    }
    row_heights <- as.numeric(row_heights)
    if (length(row_heights) != length(plots) || any(!is.finite(row_heights)) || any(row_heights <= 0)) {
      row_heights <- rep(1, length(plots))
    }

    if (requireNamespace("cowplot", quietly = TRUE)) {
      aligned <- cowplot::align_plots(plotlist = plots, align = "v", axis = "lr")
      combined <- cowplot::plot_grid(plotlist = aligned, ncol = 1, rel_heights = row_heights)
      ggsave(path, combined, width = width, height = height, device = grDevices::pdf)
      return(invisible(path))
    }

    grobs <- lapply(plots, ggplot2::ggplotGrob)
    width_lengths <- vapply(grobs, function(g) length(g$widths), integer(1))
    if (length(unique(width_lengths)) == 1L) {
      max_widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
      grobs <- lapply(grobs, function(g) {
        g$widths <- max_widths
        g
      })
    }

    grDevices::pdf(path, width = width, height = height, onefile = TRUE)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      nrow = length(grobs),
      ncol = 1,
      heights = grid::unit(row_heights, "null")
    )))
    for (i in seq_along(grobs)) {
      grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
      grid::grid.draw(grobs[[i]])
      grid::upViewport()
    }
    grid::upViewport()
    grDevices::dev.off()
    invisible(path)
  }

  sim_list <- lapply(scenarios, function(sc) {
    simulate_one_full_horizon(run_params, sc, cfg, horizon_day = horizon_day, report_dt = report_dt)
  })
  burden_all <- bind_rows(lapply(sim_list, `[[`, "burden"))
  ploidy_all <- bind_rows(lapply(sim_list, `[[`, "ploidy"))
  if (nrow(burden_all) == 0 || nrow(ploidy_all) == 0) return(invisible(NULL))
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(cfg$glucose, run_params$glucose, TRUE),
    default = TRUE
  ))
  death_language <- resource_death_language(glucose_use)

  burden_all <- burden_all %>% filter(day <= horizon_day + 1e-9)
  ploidy_all <- ploidy_all %>% filter(day <= horizon_day + 1e-9)
  burden_all <- normalize_burden_for_plot(burden_all)
  ploidy_mean <- compute_ploidy_weighted_mean(ploidy_all, cfg)

  horizon_tag <- paste0("0_", as.integer(round(horizon_day)), "day")
  # Remove deprecated multi-file prediction plot outputs for this horizon to avoid stale files.
  unlink(file.path(out_dir, c(
    paste0("predict_burden_normalized_", horizon_tag, ".pdf"),
    paste0("predict_burden_absolute_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_heatmap_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_top_states_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_weighted_mean_", horizon_tag, ".pdf"),
    paste0("forecast_burden_normalized_", horizon_tag, ".pdf"),
    paste0("forecast_burden_absolute_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_heatmap_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_top_states_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_weighted_mean_", horizon_tag, ".pdf"),
    paste0("predict_g_timecourse_", horizon_tag, ".pdf")
  )), force = TRUE)

  write.table(burden_all, file = file.path(out_dir, paste0("predict_burden_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_all, file = file.path(out_dir, paste0("predict_ploidy_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_mean, file = file.path(out_dir, paste0("predict_ploidy_weighted_mean_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  sample_day_key_from_cols <- function(harvest, cohort, dose, day) {
    paste(
      as.character(harvest),
      as.character(cohort),
      format(as.numeric(dose), trim = TRUE, scientific = FALSE),
      sprintf("%.8f", as.numeric(day)),
      sep = "__"
    )
  }

  burden_plot_df <- burden_all %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      value = as.numeric(pred_burden),
      sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
    )

  endpoint_plot_df <- ploidy_mean %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      value_chr = as.numeric(weighted_mean_N),
      sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
    )

  predict_x_breaks <- unique(as.numeric(seq(0, as.numeric(horizon_day), length.out = 5)))
  predict_x_scale <- function() {
    scale_x_continuous(
      breaks = predict_x_breaks,
      expand = c(0, 0)
    )
  }
  predict_curve_theme <- theme(
    axis.title = element_text(size = 8),
    axis.title.y.right = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.text.y.right = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 8)
  )
  ploidy_n_range <- range(endpoint_plot_df$value_chr, na.rm = TRUE)
  if (!all(is.finite(ploidy_n_range))) {
    ploidy_n_range <- c(as.numeric(cfg$N_MIN), as.numeric(cfg$N_MAX))
  }
  ploidy_n_pad <- diff(ploidy_n_range) * 0.04
  if (!is.finite(ploidy_n_pad) || ploidy_n_pad <= 0) ploidy_n_pad <- 1
  ploidy_n_limits <- c(
    floor(max(as.numeric(cfg$N_MIN), ploidy_n_range[[1]] - ploidy_n_pad)),
    ceiling(min(as.numeric(cfg$N_MAX), ploidy_n_range[[2]] + ploidy_n_pad))
  )
  if (!all(is.finite(ploidy_n_limits)) || ploidy_n_limits[[2]] <= ploidy_n_limits[[1]]) {
    ploidy_n_limits <- c(as.numeric(cfg$N_MIN), as.numeric(cfg$N_MAX))
  }
  ploidy_n_breaks <- pretty(ploidy_n_limits, n = 4)
  ploidy_n_breaks <- ploidy_n_breaks[ploidy_n_breaks >= ploidy_n_limits[[1]] & ploidy_n_breaks <= ploidy_n_limits[[2]]]

  o2_plot_min <- 0
  o2_plot_max <- 5
  o2_plot_df <- if ("pred_o2_pct" %in% names(burden_all)) {
    burden_all %>%
      filter(is.finite(pred_o2_pct)) %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        value = as.numeric(clip(pred_o2_pct, o2_plot_min, o2_plot_max)),
        sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
      )
  } else {
    data.frame()
  }

  live_cells_plot_df <- if ("pred_burden_live_cells" %in% names(burden_all)) {
    burden_all %>%
      filter(is.finite(pred_burden_live_cells)) %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        value = as.numeric(pred_burden_live_cells),
        sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
      )
  } else {
    data.frame()
  }

  chr_density_day_width <- {
    day_vals <- sort(unique(as.numeric(ploidy_all$day[is.finite(ploidy_all$day)])))
    day_step <- diff(day_vals)
    day_step <- day_step[is.finite(day_step) & day_step > 0]
    width_use <- if (length(day_step) > 0L) stats::median(day_step, na.rm = TRUE) else as.numeric(report_dt)
    if (!is.finite(width_use) || width_use <= 0) width_use <- 1
    width_use
  }
  chr_density_n_min <- as.integer(cfg$N_MIN)
  chr_density_n_max <- as.integer(cfg$N_MAX)
  chr_density_bin_width <- 5L
  chr_density_bins <- data.frame(
    chr_bin_lower = seq.int(chr_density_n_min, chr_density_n_max, by = chr_density_bin_width),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      chr_bin_upper = chr_bin_lower + chr_density_bin_width - 1L,
      chr_bin_mid = (as.numeric(chr_bin_lower) + as.numeric(chr_bin_upper)) / 2
    )
  chr_density_df <- ploidy_all %>%
    filter(
      as.character(cohort) %in% c("2N", "4N"),
      is.finite(day),
      is.finite(N),
      is.finite(fraction)
    ) %>%
    transmute(
      cohort = factor(as.character(cohort), levels = c("2N", "4N")),
      day = as.numeric(day),
      N = as.integer(round(as.numeric(N))),
      sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__"),
      fraction = pmax(as.numeric(fraction), 0)
    ) %>%
    filter(N >= chr_density_n_min, N <= chr_density_n_max) %>%
    mutate(
      chr_bin_lower = chr_density_n_min + ((N - chr_density_n_min) %/% chr_density_bin_width) * chr_density_bin_width
    ) %>%
    group_by(cohort, day, sample_id, chr_bin_lower) %>%
    summarise(
      bin_probability = sum(fraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(chr_density_bins, by = "chr_bin_lower") %>%
    group_by(cohort, day, chr_bin_lower, chr_bin_upper, chr_bin_mid) %>%
    summarise(
      density = mean(bin_probability, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    complete(
      cohort,
      day,
      tidyr::nesting(chr_bin_lower, chr_bin_upper, chr_bin_mid),
      fill = list(density = 0)
    ) %>%
    mutate(
      density = pmax(pmin(as.numeric(density), 1), 0)
    )
  write.table(
    chr_density_df,
    file = file.path(out_dir, paste0("predict_chromosome_density_", horizon_tag, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  p_predict_burden <- ggplot(
    burden_plot_df,
    aes(x = day, y = value, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    predict_x_scale() +
    coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = paste0("Predict Curves: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = paste0("Single summary plot (all scenarios overlaid) | fit_dir=", basename(dirname(out_dir)), " | report_dt=", report_dt),
      x = "Day",
      y = "Burden",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank()
    ) +
    predict_curve_theme

  p_predict_endpoint <- ggplot(
    endpoint_plot_df,
    aes(x = day, y = value_chr, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    predict_x_scale() +
    coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
    scale_y_continuous(
      name = "Mean chr. number",
      breaks = ploidy_n_breaks,
      sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Mean ploidy")
    ) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      x = "Day",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank()
    ) +
    predict_curve_theme

  p_predict_chr_density <- if (nrow(chr_density_df) > 0L) {
    ggplot(
      chr_density_df,
      aes(x = day, y = chr_bin_mid, fill = density)
    ) +
      geom_tile(width = chr_density_day_width, height = chr_density_bin_width) +
      facet_grid(cohort ~ .) +
      predict_x_scale() +
      scale_fill_gradientn(
        colors = grDevices::hcl.colors(256, palette = "Viridis"),
        limits = c(0, 1),
        breaks = c(0, 0.5, 1),
        name = "Probability\ndensity",
        na.value = "white"
      ) +
      scale_y_continuous(
        breaks = ploidy_n_breaks,
        expand = c(0, 0),
        sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Ploidy")
      ) +
      coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
      labs(
        x = "Day",
        y = "Chromosome N"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "grey80")
      ) +
      predict_curve_theme
  } else {
    ggplot() +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
      scale_y_continuous(
        breaks = ploidy_n_breaks,
        sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Ploidy")
      ) +
      labs(
        x = "Day",
        y = "Chromosome N"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      predict_curve_theme
  }

  p_predict_o2 <- if (nrow(o2_plot_df) > 0L) {
    ggplot(
      o2_plot_df,
      aes(x = day, y = value, group = sample_id, color = cohort)
    ) +
      geom_line(linewidth = 0.65, alpha = 0.8) +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_y_continuous(
        limits = c(o2_plot_min, o2_plot_max),
        breaks = seq(o2_plot_min, o2_plot_max, by = 1),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
      labs(
        x = "Day",
        y = "O2 (%)",
        color = "Cohort"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      predict_curve_theme
  } else {
    ggplot() +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_y_continuous(
        limits = c(o2_plot_min, o2_plot_max),
        breaks = seq(o2_plot_min, o2_plot_max, by = 1),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      labs(
        x = "Day",
        y = "O2 (%)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      predict_curve_theme
  }

  p_predict_live_cells <- if (nrow(live_cells_plot_df) > 0L) {
    ggplot(
      live_cells_plot_df,
      aes(x = day, y = value, group = sample_id, color = cohort)
    ) +
      geom_line(linewidth = 0.65, alpha = 0.8) +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
      labs(
        x = "Day",
        y = "Live cells",
        color = "Cohort"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      predict_curve_theme
  } else {
    ggplot() +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      labs(
        x = "Day",
        y = "Live cells"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      predict_curve_theme
  }

  predict_curves_pdf <- file.path(out_dir, paste0("predict_curves_", horizon_tag, ".pdf"))
  save_aligned_plot_stack(
    list(p_predict_burden, p_predict_endpoint, p_predict_chr_density, p_predict_o2, p_predict_live_cells),
    predict_curves_pdf,
    width = 12,
    height = 9,
    row_heights = c(1, 1, 1.5, 1, 1)
  )

  p_live_weighted_pms_predict <- NULL
  if ("pred_o2_pct" %in% names(burden_all)) {
    predict_o2_by_sample_day <- burden_all %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__"),
        sample_day_key = sample_day_key_from_cols(harvest, cohort, dose, day),
        o2_pct = as.numeric(clip(pred_o2_pct, o2_plot_min, o2_plot_max))
      )
    predict_ploidy_with_o2 <- ploidy_all %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        N = as.numeric(N),
        fraction = as.numeric(fraction),
        sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__"),
        sample_day_key = sample_day_key_from_cols(harvest, cohort, dose, day)
      ) %>%
      left_join(
        predict_o2_by_sample_day %>% select(sample_day_key, o2_pct),
        by = "sample_day_key"
      ) %>%
      filter(
        cohort %in% c("2N", "4N"),
        is.finite(N),
        is.finite(fraction),
        fraction > 0,
        is.finite(o2_pct)
      ) %>%
      mutate(
        live_state_p_ms = as.numeric(.pmisseg_of_O2(
          O2 = o2_pct,
          run_params = run_params,
          N = N,
          O2_crit = as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0))
        ))
      )
    predict_live_weighted_pms_df <- predict_ploidy_with_o2 %>%
      group_by(harvest, cohort, dose, day, sample_id, o2_pct) %>%
      summarise(
        live_total_fraction = sum(fraction, na.rm = TRUE),
        weighted_p_ms = sum(fraction * live_state_p_ms, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        live_weighted_effective_p_ms = weighted_p_ms / pmax(live_total_fraction, 1e-12),
        cohort = factor(cohort, levels = c("2N", "4N"))
      ) %>%
      arrange(cohort, harvest, dose, day)
    if (nrow(predict_live_weighted_pms_df) > 0L) {
      write.table(
        predict_live_weighted_pms_df,
        file = file.path(out_dir, paste0("predict_live_weighted_pms_", horizon_tag, ".tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      p_live_weighted_pms_predict <- ggplot(
        predict_live_weighted_pms_df,
        aes(x = day, y = live_weighted_effective_p_ms, group = sample_id, color = cohort)
      ) +
        geom_line(linewidth = 0.65, alpha = 0.85) +
        facet_wrap(~ harvest, ncol = 2) +
        coord_cartesian(
          xlim = c(0, horizon_day),
          ylim = c(0, max(predict_live_weighted_pms_df$live_weighted_effective_p_ms, na.rm = TRUE))
        ) +
        scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
        labs(
          title = paste0("Predict Live-Weighted Effective p_ms: 0-", as.integer(round(horizon_day)), " days"),
          subtitle = paste0("Sample-day live-cell effective p_ms under fitted parameters | fit_dir=", basename(dirname(out_dir)), " | report_dt=", report_dt),
          x = "Day",
          y = "Live-weighted effective p_ms",
          color = "Cohort"
        ) +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank())
      ggsave(
        file.path(out_dir, paste0("predict_live_weighted_pms_", horizon_tag, ".pdf")),
        p_live_weighted_pms_predict,
        width = 12,
        height = 9
      )
    }
  }

  burden_decomp_predict <- burden_all %>%
    transmute(
      cohort = factor(as.character(cohort), levels = c("2N", "4N")),
      harvest = as.character(harvest),
      dose = as.numeric(dose),
      day = as.numeric(day),
      burden_live = as.numeric(.first_non_null_local(pred_burden_live_volume_mm3, pred_burden)),
      burden_dead_hypoxia = as.numeric(.first_non_null_local(pred_burden_dead_hypoxia_volume_mm3, 0)),
      burden_dead_buffer = as.numeric(.first_non_null_local(pred_burden_dead_buffer_volume_mm3, 0)),
      burden_dead_total = as.numeric(.first_non_null_local(
        pred_burden_dead_total_volume_mm3,
        .first_non_null_local(pred_burden_dead_hypoxia_volume_mm3, 0) + .first_non_null_local(pred_burden_dead_buffer_volume_mm3, 0)
      )),
      burden_total = as.numeric(pred_burden)
    )

  death_ratio_predict <- burden_decomp_predict %>%
    mutate(
      sample_id = paste(harvest, as.character(cohort), format(dose, trim = TRUE, scientific = FALSE), sep = "__"),
      death_ratio = burden_dead_total / pmax(burden_total, 1e-12)
    ) %>%
    mutate(death_ratio = pmax(0, pmin(1, as.numeric(death_ratio))))
  write.table(
    death_ratio_predict,
    file = file.path(out_dir, paste0("predict_death_ratio_", horizon_tag, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  death_ratio_summary <- death_ratio_predict %>%
    filter(!is.na(cohort)) %>%
    group_by(cohort, day) %>%
    summarise(
      death_ratio = mean(death_ratio, na.rm = TRUE),
      .groups = "drop"
    )

  p_death_ratio_predict <- ggplot(
    death_ratio_predict,
    aes(x = day, y = death_ratio, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.45, alpha = 0.30) +
    geom_line(
      data = death_ratio_summary,
      aes(x = day, y = death_ratio, group = cohort, color = cohort),
      inherit.aes = FALSE,
      linewidth = 1.1,
      alpha = 1.0
    ) +
    coord_cartesian(xlim = c(0, horizon_day), ylim = c(0, 1)) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728"), drop = FALSE) +
    labs(
      title = paste0("Predict Death Ratio: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = "Thin lines = individual scenarios, thick lines = cohort mean",
      x = "Day",
      y = "Death ratio",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, paste0("predict_death_ratio_", horizon_tag, ".pdf")), p_death_ratio_predict, width = 12, height = 7)

  burden_decomp_predict <- burden_decomp_predict %>%
    filter(!is.na(cohort)) %>%
    group_by(cohort, day) %>%
    summarise(
      burden_live = mean(burden_live, na.rm = TRUE),
      burden_dead_hypoxia = mean(burden_dead_hypoxia, na.rm = TRUE),
      burden_dead_buffer = mean(burden_dead_buffer, na.rm = TRUE),
      burden_total = mean(burden_total, na.rm = TRUE),
      .groups = "drop"
    )

  burden_decomp_predict_long <- burden_decomp_predict %>%
    pivot_longer(
      cols = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
      names_to = "component",
      values_to = "value"
    ) %>%
    mutate(
      component = factor(
        component,
        levels = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
        labels = c("Live", death_language$component_label, "Dead (Buffer loss)")
      )
    )

  p_burden_decomp_predict <- ggplot(
    burden_decomp_predict_long,
    aes(x = day, y = value, fill = component, group = component)
  ) +
    geom_area(alpha = 0.55, position = "stack") +
    geom_line(
      data = burden_decomp_predict,
      aes(x = day, y = burden_total),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.65
    ) +
    facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = stats::setNames(c("#1f77b4", "#d62728", "#2ca02c"), c("Live", death_language$component_label, "Dead (Buffer loss)"))) +
    coord_cartesian(xlim = c(0, horizon_day)) +
    labs(
      title = paste0("Predict Burden Live/Dead Decomposition: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = "Cohort-level mean across scenarios (2N top, 4N bottom)",
      x = "Day",
      y = "Tumor burden (mm^3)",
      fill = "Component"
    ) +
    theme_bw(base_size = 11)

  burden_decomp_plot_for_cohort <- function(cohort_use, show_fill_legend = TRUE) {
    row_df <- burden_decomp_predict %>%
      filter(cohort == cohort_use)
    row_long <- burden_decomp_predict_long %>%
      filter(cohort == cohort_use)
    ggplot(
      row_long,
      aes(x = day, y = value, fill = component, group = component)
    ) +
      geom_area(alpha = 0.55, position = "stack") +
      geom_line(
        data = row_df,
        aes(x = day, y = burden_total),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.65
      ) +
      scale_fill_manual(
        values = stats::setNames(
          c("#1f77b4", "#d62728", "#2ca02c"),
          c("Live", death_language$component_label, "Dead (Buffer loss)")
        )
      ) +
      coord_cartesian(xlim = c(0, horizon_day)) +
      labs(
        title = if (identical(as.character(cohort_use), "2N")) {
          paste0("Predict Burden Live/Dead Decomposition: 0-", as.integer(round(horizon_day)), " days")
        } else {
          NULL
        },
        subtitle = paste0(as.character(cohort_use), " cohort mean across scenarios"),
        x = "Day",
        y = "Tumor burden (mm^3)",
        fill = "Component"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (show_fill_legend) "right" else "none"
      )
  }

  save_aligned_plot_stack(
    list(
      burden_decomp_plot_for_cohort("2N", show_fill_legend = TRUE),
      burden_decomp_plot_for_cohort("4N", show_fill_legend = TRUE)
    ),
    file.path(out_dir, paste0("predict_burden_live_dead_decomposition_", horizon_tag, ".pdf")),
    width = 14,
    height = 7
  )

  p_o2_time <- NULL
  if (all(c("pred_o2_target_pct", "pred_o2_pct") %in% names(burden_all))) {
    o2_plot_min <- 0
    o2_plot_max <- 5
    o2_time_df <- burden_all %>%
      filter(is.finite(pred_o2_target_pct), is.finite(pred_o2_pct)) %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        sample_id = paste(as.character(harvest), as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__"),
        o2_target_pct = as.numeric(clip(pred_o2_target_pct, o2_plot_min, o2_plot_max)),
        o2_eff_pct = as.numeric(clip(pred_o2_pct, o2_plot_min, o2_plot_max))
      )
    o2_time_long <- o2_time_df %>%
      pivot_longer(cols = c("o2_target_pct", "o2_eff_pct"), names_to = "o2_series", values_to = "o2_pct") %>%
      mutate(o2_series = factor(o2_series, levels = c("o2_target_pct", "o2_eff_pct"), labels = c("O2_target", "O2_eff")))

    p_o2_time <- ggplot(o2_time_long, aes(x = day, y = o2_pct, color = o2_series, linetype = o2_series, group = interaction(sample_id, o2_series))) +
      geom_line(linewidth = 0.65, alpha = 0.85) +
      facet_wrap(~ harvest, ncol = 2) +
      coord_cartesian(xlim = c(0, horizon_day), ylim = c(o2_plot_min, o2_plot_max)) +
      scale_color_manual(values = c("O2_target" = "#ff7f0e", "O2_eff" = "#1f77b4")) +
      labs(
        title = paste0("Oxygen Evolution Over Time (0-", as.integer(round(horizon_day)), " days)"),
        x = "Day",
        y = "Oxygen (%)",
        color = NULL,
        linetype = NULL
      ) +
      theme_bw(base_size = 11)

    ggsave(file.path(out_dir, paste0("predict_o2_timecourse_", horizon_tag, ".pdf")), p_o2_time, width = 12, height = 9)
  }
  invisible(list(
    p_predict = p_predict_endpoint,
    p_live_weighted_pms_predict = p_live_weighted_pms_predict,
    p_o2_time = p_o2_time,
    p_burden_decomp_predict = p_burden_decomp_predict,
    p_death_ratio_predict = p_death_ratio_predict,
    burden_decomp_predict = burden_decomp_predict,
    burden_decomp_predict_long = burden_decomp_predict_long,
    horizon_day = as.numeric(horizon_day)
  ))
}

plot_predict_burden_live_dead_decomposition_combined <- function(predict_results, out_dir, death_language) {
  predict_results <- Filter(function(x) {
    is.list(x) &&
      is.data.frame(x$burden_decomp_predict) &&
      is.data.frame(x$burden_decomp_predict_long) &&
      is.finite(as.numeric(x$horizon_day))
  }, predict_results)
  if (!length(predict_results)) return(invisible(NULL))

  predict_results <- predict_results[order(vapply(predict_results, function(x) as.numeric(x$horizon_day), numeric(1)))]
  fill_values <- stats::setNames(
    c("#1f77b4", "#d62728", "#2ca02c"),
    c("Live", death_language$component_label, "Dead (Buffer loss)")
  )

  plots <- lapply(predict_results, function(res) {
    horizon_day <- as.numeric(res$horizon_day)
    ggplot(
      res$burden_decomp_predict_long,
      aes(x = day, y = value, fill = component, group = component)
    ) +
      geom_area(alpha = 0.55, position = "stack") +
      geom_line(
        data = res$burden_decomp_predict,
        aes(x = day, y = burden_total),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.65
      ) +
      facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = fill_values, drop = FALSE) +
      coord_cartesian(xlim = c(0, horizon_day)) +
      labs(
        title = paste0("0-", as.integer(round(horizon_day)), " days"),
        subtitle = "Cohort-level mean across scenarios (2N top, 4N bottom)",
        x = "Day",
        y = "Tumor burden (mm^3)",
        fill = "Component"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "right"
      )
  })

  save_aligned_plot_row_with_shared_legend(
    plots,
    file.path(out_dir, "predict_burden_live_dead_decomposition_combined.pdf"),
    width = 18,
    height = 6.5
  )
}

# -----------------------------------------------------------------------------
# Function: find_fit_dirs_under
# Purpose: Discover fit result directories from a root path.
# Parameters:
#   - root_dir: Root directory used for recursive fit-folder discovery.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
find_fit_dirs_under <- function(root_dir) {
  all_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  sub_dirs <- all_dirs[all_dirs != root_dir]

# -----------------------------------------------------------------------------
# Function: is_fit_dir
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - d: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  is_fit_dir <- function(d) {
    file.exists(file.path(d, "fit_config.rds")) &&
      file.exists(file.path(d, "best_params.tsv"))
  }

  # Primary mode: traverse subdirectories under fit_dir.
  fit_sub_dirs <- sub_dirs[vapply(sub_dirs, is_fit_dir, logical(1))]
  if (length(fit_sub_dirs) > 0) {
    return(sort(unique(fit_sub_dirs)))
  }

  # Fallback mode: if no fit subdirectories exist, use fit_dir itself.
  if (is_fit_dir(root_dir)) {
    return(normalizePath(root_dir, mustWork = TRUE))
  }

  character(0)
}

read_fit_summary_value_for_viz <- function(fit_dir, key, default = NA_character_) {
  path <- file.path(fit_dir, "fit_summary.tsv")
  if (!file.exists(path)) return(default)
  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || !all(c("metric", "value") %in% names(tab))) return(default)
  hit <- tab$value[as.character(tab$metric) == as.character(key)]
  if (!length(hit) || is.na(hit[[1]])) return(default)
  trimws(as.character(hit[[1]]))
}

is_joint_fit_dir_for_viz <- function(fit_dir) {
  identical(read_fit_summary_value_for_viz(fit_dir, "fit_mode", default = "fit_invivo"), "fit_joint")
}

cleanup_joint_viz_root <- function(viz_root) {
  dir.create(viz_root, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(viz_root, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (!length(entries)) return(invisible(NULL))
  keep <- basename(entries) %in% c("invivo", "invitro")
  unlink(entries[!keep], recursive = TRUE, force = TRUE)
  invisible(NULL)
}

resolve_invivo_viz_out_dir <- function(fit_dir, out_dir_override = NULL) {
  out_dir_override <- if (is.null(out_dir_override)) NULL else trimws(as.character(out_dir_override[[1]]))
  if (!is.null(out_dir_override) && nzchar(out_dir_override)) {
    return(normalizePath(out_dir_override, mustWork = FALSE))
  }
  viz_root <- file.path(fit_dir, "viz")
  if (is_joint_fit_dir_for_viz(fit_dir)) {
    return(file.path(viz_root, "invivo"))
  }
  viz_root
}

# -----------------------------------------------------------------------------
# Function: run_viz_for_fit_dir
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - fit_dir: Directory containing fitted parameters and summary outputs.
#   - argv: Character vector of command-line arguments in --key=value format.
#   - dt_path: Path to burden observation table (Excel file).
#   - ploidy_path: Path to terminal ploidy table (TSV file).
#   - report_dt: Sampling interval for reported trajectories.
#   - top_n: Maximum number of scenarios selected for detailed plotting.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
run_viz_for_fit_dir <- function(
  fit_dir,
  argv,
  dt_path,
  ploidy_path,
  report_dt,
  top_n,
  out_dir_override = NULL
) {
  joint_fit <- is_joint_fit_dir_for_viz(fit_dir)
  viz_root <- file.path(fit_dir, "viz")
  if (isTRUE(joint_fit)) {
    cleanup_joint_viz_root(viz_root)
  }
  out_dir <- resolve_invivo_viz_out_dir(fit_dir, out_dir_override = out_dir_override)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cfg_path <- file.path(fit_dir, "fit_config.rds")
  if (!file.exists(cfg_path)) stop("Missing fit_config.rds in: ", fit_dir)
  cfg <- readRDS(cfg_path)
  cfg <- normalize_cfg_for_viz(cfg)

  if (!is.null(argv$dose_zero_only)) cfg$dose_zero_only <- as_bool(argv$dose_zero_only, cfg$dose_zero_only)
  if (!is.null(argv$truncate_at_treatment)) cfg$truncate_at_treatment <- as_bool(argv$truncate_at_treatment, cfg$truncate_at_treatment)
  if (!is.null(argv$ploidy_at_harvest)) cfg$ploidy_at_harvest <- as_bool(argv$ploidy_at_harvest, cfg$ploidy_at_harvest)
  if (!is.null(argv$max_scenarios)) cfg$max_scenarios <- as_num(argv$max_scenarios, cfg$max_scenarios)

  run_params <- read_run_params(fit_dir, cfg = cfg)
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(cfg$glucose, run_params$glucose, TRUE),
    default = TRUE
  ))
  death_language <- resource_death_language(glucose_use)
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)

  sim_list <- lapply(scenarios, function(sc) simulate_one_full(run_params, sc, cfg, report_dt = report_dt))
  burden_all <- bind_rows(lapply(sim_list, `[[`, "burden"))
  ploidy_all <- bind_rows(lapply(sim_list, `[[`, "ploidy"))

  if (nrow(burden_all) == 0 || nrow(ploidy_all) == 0) {
    stop("No simulation output generated; check fit/data configuration.")
  }

  burden_all <- normalize_burden_for_plot(burden_all)

  write.table(burden_all, file = file.path(out_dir, "burden_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_all, file = file.path(out_dir, "ploidy_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  fixed_glucose_specs <- fixed_glucose_diagnostic_specs(
    cfg = cfg,
    run_params = run_params
  )
  if (nrow(fixed_glucose_specs) > 0L) {
    message(
      "  Fixed-glucose diagnostics: ",
      paste(
        paste0("G=", fixed_glucose_specs$name, "=", signif(fixed_glucose_specs$value, 6), "%"),
        collapse = "; "
      ),
      "; O2 varies"
    )
  }

  burden_decomp <- burden_all %>%
    mutate(
      burden_live = as.numeric(.first_non_null_local(pred_burden_live_volume_mm3, pred_burden)),
      burden_dead_hypoxia = as.numeric(.first_non_null_local(pred_burden_dead_hypoxia_volume_mm3, 0)),
      burden_dead_buffer = as.numeric(.first_non_null_local(pred_burden_dead_buffer_volume_mm3, 0)),
      burden_dead_total = as.numeric(.first_non_null_local(pred_burden_dead_total_volume_mm3, burden_dead_hypoxia + burden_dead_buffer)),
      burden_total = as.numeric(pred_burden)
    ) %>%
    select(harvest, cohort, dose, day, burden_live, burden_dead_hypoxia, burden_dead_buffer, burden_dead_total, burden_total)
  write.table(burden_decomp, file = file.path(out_dir, "burden_live_dead_decomposition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  ploidy_mean <- compute_ploidy_weighted_mean(ploidy_all, cfg)
  write.table(ploidy_mean, file = file.path(out_dir, "ploidy_weighted_mean_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  terminal_ploidy_compare <- build_terminal_ploidy_compare_df(scenarios = scenarios, ploidy_all = ploidy_all, cfg = cfg)
  if (nrow(terminal_ploidy_compare) > 0) {
    write.table(
      terminal_ploidy_compare,
      file = file.path(out_dir, "terminal_ploidy_observed_vs_predicted.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    plot_terminal_ploidy_violin_compare(
      compare_df = terminal_ploidy_compare,
      fit_dir = fit_dir,
      out_dir = out_dir
    )
  }

  has_o2_lag_cols <- all(c("pred_o2_target_pct", "pred_o2_pct") %in% names(burden_all))
  if (isTRUE(has_o2_lag_cols)) {
    o2_plot_min <- 0
    o2_plot_max <- 5
    o2_lag_df <- burden_all %>%
      filter(is.finite(pred_o2_target_pct), is.finite(pred_o2_pct)) %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        sample_id = paste(as.character(harvest), as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__"),
        o2_target_pct = as.numeric(clip(pred_o2_target_pct, o2_plot_min, o2_plot_max)),
        o2_eff_pct = as.numeric(clip(pred_o2_pct, o2_plot_min, o2_plot_max)),
        o2_lag_gap_pct = as.numeric(clip(pred_o2_target_pct, o2_plot_min, o2_plot_max) - clip(pred_o2_pct, o2_plot_min, o2_plot_max))
      )
    write.table(o2_lag_df, file = file.path(out_dir, "o2_lag_timecourse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    o2_lag_long <- o2_lag_df %>%
      select(harvest, cohort, dose, day, sample_id, o2_target_pct, o2_eff_pct) %>%
      pivot_longer(cols = c("o2_target_pct", "o2_eff_pct"), names_to = "o2_series", values_to = "o2_pct") %>%
      mutate(o2_series = factor(o2_series, levels = c("o2_target_pct", "o2_eff_pct"), labels = c("O2_target", "O2_eff")))

    p_o2_lag <- ggplot(o2_lag_long, aes(x = day, y = o2_pct, color = o2_series, linetype = o2_series, group = interaction(sample_id, o2_series))) +
      geom_line(linewidth = 0.7, alpha = 0.85) +
      facet_wrap(~ harvest, ncol = 2) +
      scale_color_manual(values = c("O2_target" = "#ff7f0e", "O2_eff" = "#1f77b4")) +
      coord_cartesian(ylim = c(o2_plot_min, o2_plot_max)) +
      labs(
        title = "O2 Supply-Demand MAP Model: Oxygen Lag Over Time",
        subtitle = "O2_target (instantaneous) vs O2_eff (lagged state)",
        x = "Day",
        y = "Oxygen (%)",
        color = NULL,
        linetype = NULL
      ) +
      theme_bw(base_size = 11)

    p_o2_lag_gap <- ggplot(o2_lag_df, aes(x = day, y = o2_lag_gap_pct, color = cohort, group = sample_id)) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4, linetype = "dashed") +
      geom_line(linewidth = 0.7, alpha = 0.85) +
      facet_wrap(~ harvest, ncol = 2) +
      scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
      labs(
        title = "O2 Supply-Demand MAP Model: O2 Lag Gap Over Time",
        subtitle = "Lag gap = O2_target - O2_eff",
        x = "Day",
        y = "Lag gap (percentage points)",
        color = "Cohort"
      ) +
      theme_bw(base_size = 11)
  }
  unlink(file.path(out_dir, c(
    "g_lag_timecourse.tsv",
    "g_target_vs_eff_timecourse.pdf",
    "g_lag_gap_timecourse.pdf",
    "predict_burden_vs_g.tsv",
    "predict_burden_vs_g.pdf",
    "glucose_response_4panel.pdf"
  )), force = TRUE)

  p_burden <- ggplot(burden_all, aes(x = day, y = pred_norm)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_line(
      data = burden_all %>% filter(!is.na(obs_norm)),
      aes(y = obs_norm),
      color = "black",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      data = burden_all %>% filter(!is.na(obs_norm)),
      aes(y = obs_norm),
      color = "black",
      size = 1
    ) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      title = "O2 Supply-Demand MAP Model: In Vivo Burden Trajectory (Normalized)",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Day",
      y = "Normalized Burden (delta / max|delta|)"
    ) +
    theme_bw(base_size = 11)

  p_burden_abs <- ggplot(burden_all, aes(x = day, y = pred_burden)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_line(
      data = burden_all %>% filter(!is.na(obs_burden)),
      aes(y = obs_burden),
      color = "black",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      data = burden_all %>% filter(!is.na(obs_burden)),
      aes(y = obs_burden),
      color = "black",
      size = 1
    ) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_y") +
    labs(
      title = "O2 Supply-Demand MAP Model: In Vivo Burden Trajectory (Absolute)",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Day",
      y = "Tumor burden (mm^3)"
    ) +
    theme_bw(base_size = 11)

  rho_2N_min <- as_num(argv$rho_2N_min, 3.2e4)
  rho_2N_max <- as_num(argv$rho_2N_max, 5.6e4)
  if (!is.finite(rho_2N_min) || rho_2N_min <= 0) rho_2N_min <- 3.2e4
  if (!is.finite(rho_2N_max) || rho_2N_max <= 0) rho_2N_max <- 5.6e4
  if (rho_2N_min > rho_2N_max) {
    tmp <- rho_2N_min
    rho_2N_min <- rho_2N_max
    rho_2N_max <- tmp
  }
  rho_2N_mid <- sqrt(rho_2N_min * rho_2N_max)
  pred_cell_col <- if ("pred_burden_cells" %in% names(burden_all)) "pred_burden_cells" else "pred_burden"
  burden_all_real <- burden_all %>%
    mutate(
      pred_burden_cell_number = as.numeric(.data[[pred_cell_col]]),
      obs_burden_cell_number_low = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_min, NA_real_),
      obs_burden_cell_number_mid = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_mid, NA_real_),
      obs_burden_cell_number_high = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_max, NA_real_)
    )

  p_burden_abs_real <- ggplot(burden_all_real, aes(x = day, y = pred_burden_cell_number)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_ribbon(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_low) & !is.na(obs_burden_cell_number_high)),
      aes(x = day, ymin = obs_burden_cell_number_low, ymax = obs_burden_cell_number_high),
      inherit.aes = FALSE,
      fill = "grey50",
      alpha = 0.18
    ) +
    geom_line(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_mid)),
      aes(y = obs_burden_cell_number_mid),
      color = "black",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_mid)),
      aes(y = obs_burden_cell_number_mid),
      color = "black",
      size = 1
    ) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_y") +
    labs(
      title = "O2 Supply-Demand MAP Model: In Vivo Burden Trajectory (Absolute, Real Scale)",
      subtitle = paste0(
        "fit_dir=", basename(fit_dir),
        " | report_dt=", report_dt,
        " | obs burden -> CellNumber using rho_2N range=[",
        signif(rho_2N_min, 4), ", ", signif(rho_2N_max, 4), "] cells/mm^3",
        " (mid=", signif(rho_2N_mid, 4), ")"
      ),
      x = "Day",
      y = "CellNumber (2N-equivalent range)"
    ) +
    theme_bw(base_size = 11)

  burden_decomp_long <- burden_decomp %>%
    pivot_longer(
      cols = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
      names_to = "component",
      values_to = "value"
    ) %>%
    mutate(
      component = factor(
        component,
        levels = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
        labels = c("Live", death_language$component_label, "Dead (Buffer loss)")
      )
    )
  p_burden_decomp <- ggplot(burden_decomp_long, aes(x = day, y = value, fill = component, group = interaction(component, harvest, cohort, dose))) +
    geom_area(alpha = 0.55, position = "stack") +
    geom_line(
      data = burden_decomp,
      aes(x = day, y = burden_total, group = interaction(harvest, cohort, dose)),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.6
    ) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = stats::setNames(c("#1f77b4", "#d62728", "#2ca02c"), c("Live", death_language$component_label, "Dead (Buffer loss)"))) +
    labs(
      title = "O2 Supply-Demand MAP Model: Live/Dead Burden Decomposition",
      subtitle = paste0(
        "Total burden (black) = live + ",
        death_language$figure_phrase,
        " + dead from buffer-derived nonviable offspring"
      ),
      x = "Day",
      y = "Tumor burden (mm^3)",
      fill = "Component"
    ) +
    theme_bw(base_size = 11)

  o2_plot_min <- 0
  o2_plot_max <- 5
  o2_burden_df <- burden_all %>%
    filter(is.finite(pred_burden), is.finite(pred_o2_pct)) %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      burden_mm3 = as.numeric(pred_burden),
      o2_pct = as.numeric(clip(pred_o2_pct, o2_plot_min, o2_plot_max)),
      sample_id = paste(as.character(harvest), as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__")
  )
  write.table(o2_burden_df, file = file.path(out_dir, "predict_burden_vs_o2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sample_day_key_from_cols <- function(harvest, cohort, dose, day) {
    paste(
      as.character(harvest),
      as.character(cohort),
      format(as.numeric(dose), trim = TRUE, scientific = FALSE),
      sprintf("%.8f", as.numeric(day)),
      sep = "__"
    )
  }
  build_supported_surface_local <- function(df, o2_grid, live_fraction_grid) {
    df_use <- df %>%
      filter(
        is.finite(o2_pct),
        is.finite(total_live_fraction),
        is.finite(live_weighted_effective_p_ms)
      ) %>%
      group_by(o2_pct, total_live_fraction) %>%
      summarise(
        live_weighted_effective_p_ms = mean(live_weighted_effective_p_ms, na.rm = TRUE),
        n_points = dplyr::n(),
        .groups = "drop"
      )
    if (nrow(df_use) < 3L || length(unique(df_use$o2_pct)) < 2L || length(unique(df_use$total_live_fraction)) < 2L) {
      return(data.frame())
    }
    x <- as.numeric(df_use$o2_pct)
    y <- as.numeric(df_use$total_live_fraction)
    z <- as.numeric(df_use$live_weighted_effective_p_ms)
    x_min <- min(x, na.rm = TRUE)
    y_min <- min(y, na.rm = TRUE)
    x_scale <- max(diff(range(x, na.rm = TRUE)), 1e-6)
    y_scale <- max(diff(range(y, na.rm = TRUE)), 1e-6)
    obs_xy_scaled <- cbind((x - x_min) / x_scale, (y - y_min) / y_scale)
    obs_dist <- as.matrix(stats::dist(obs_xy_scaled))
    diag(obs_dist) <- Inf
    k_neigh <- min(5L, nrow(obs_xy_scaled) - 1L)
    support_radius <- if (k_neigh >= 1L) {
      kth_dist <- apply(obs_dist, 1, function(v) sort(v, partial = k_neigh)[[k_neigh]])
      max(as.numeric(stats::quantile(kth_dist, probs = 0.9, na.rm = TRUE)) * 1.75, 0.08)
    } else {
      0.12
    }
    kernel_radius <- max(support_radius * 0.85, 0.05)
    grid_df <- expand.grid(
      o2_pct = o2_grid,
      total_live_fraction = live_fraction_grid,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    grid_x_scaled <- (as.numeric(grid_df$o2_pct) - x_min) / x_scale
    grid_y_scaled <- (as.numeric(grid_df$total_live_fraction) - y_min) / y_scale
    dist_mat <- sqrt(
      outer(grid_x_scaled, obs_xy_scaled[, 1], "-")^2 +
        outer(grid_y_scaled, obs_xy_scaled[, 2], "-")^2
    )
    nearest_dist <- apply(dist_mat, 1, min)
    inside_bbox <- (
      as.numeric(grid_df$o2_pct) >= min(x, na.rm = TRUE) &
        as.numeric(grid_df$o2_pct) <= max(x, na.rm = TRUE) &
        as.numeric(grid_df$total_live_fraction) >= min(y, na.rm = TRUE) &
        as.numeric(grid_df$total_live_fraction) <= max(y, na.rm = TRUE)
    )
    support_mask <- inside_bbox & is.finite(nearest_dist) & (nearest_dist <= support_radius)
    weights <- exp(-0.5 * (dist_mat / pmax(kernel_radius, 1e-6))^2)
    weights[dist_mat > support_radius] <- 0
    w_sum <- rowSums(weights)
    z_pred <- as.numeric(weights %*% z) / pmax(w_sum, 1e-12)
    z_pred[!support_mask | !is.finite(z_pred)] <- NA_real_
    grid_df$inside_support <- as.logical(support_mask)
    grid_df$nearest_norm_dist <- as.numeric(nearest_dist)
    grid_df$support_radius_norm <- as.numeric(support_radius)
    grid_df$kernel_radius_norm <- as.numeric(kernel_radius)
    grid_df$live_weighted_effective_p_ms_surface <- as.numeric(z_pred)
    grid_df$n_input_points <- as.integer(nrow(df_use))
    grid_df
  }
  ploidy_with_o2 <- ploidy_all %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      N = as.numeric(N),
      fraction = as.numeric(fraction),
      sample_day_key = sample_day_key_from_cols(harvest, cohort, dose, day)
    ) %>%
    left_join(
      o2_burden_df %>%
        transmute(
          sample_day_key = sample_day_key_from_cols(harvest, cohort, dose, day),
          o2_pct = as.numeric(o2_pct)
        ),
      by = "sample_day_key"
    ) %>%
    filter(
      cohort %in% c("2N", "4N"),
      is.finite(N),
      is.finite(fraction),
      fraction > 0,
      is.finite(o2_pct)
    ) %>%
    mutate(
      live_state_p_ms = as.numeric(.pmisseg_of_O2(
        O2 = o2_pct,
        run_params = run_params,
        N = N,
        O2_crit = as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0))
      ))
    )
  live_weighted_pms_sample_day_df <- ploidy_with_o2 %>%
    group_by(sample_day_key, harvest, cohort, dose, day, o2_pct) %>%
    summarise(
      state_fraction_sum = sum(fraction, na.rm = TRUE),
      weighted_p_ms = sum(fraction * live_state_p_ms, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      live_weighted_effective_p_ms = weighted_p_ms / pmax(state_fraction_sum, 1e-12)
    )
  live_fraction_df <- burden_decomp %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      sample_day_key = sample_day_key_from_cols(harvest, cohort, dose, day),
      burden_live = as.numeric(burden_live),
      burden_total = as.numeric(burden_total),
      total_live_fraction = as.numeric(burden_live) / pmax(as.numeric(burden_total), 1e-12)
    )
  live_state_pms_o2_df <- live_weighted_pms_sample_day_df %>%
    left_join(
      live_fraction_df,
      by = c("sample_day_key", "harvest", "cohort", "dose", "day")
    ) %>%
    mutate(
      cohort = factor(cohort, levels = c("2N", "4N")),
      total_live_fraction = pmax(pmin(as.numeric(total_live_fraction), 1), 0),
      live_weighted_effective_p_ms = pmax(as.numeric(live_weighted_effective_p_ms), 0)
    ) %>%
    filter(
      cohort %in% c("2N", "4N"),
      is.finite(o2_pct),
      is.finite(total_live_fraction),
      is.finite(live_weighted_effective_p_ms)
    ) %>%
    arrange(cohort, harvest, dose, day)
  surface_grid_df <- dplyr::bind_rows(lapply(c("2N", "4N"), function(cohort_i) {
    grid_i <- build_supported_surface_local(
      df = live_state_pms_o2_df %>% filter(as.character(cohort) == cohort_i),
      o2_grid = seq(o2_plot_min, o2_plot_max, by = 0.05),
      live_fraction_grid = seq(0, 1, by = 0.01)
    )
    if (!nrow(grid_i)) return(data.frame())
    grid_i$cohort <- factor(cohort_i, levels = c("2N", "4N"))
    grid_i
  }))
  unlink(file.path(out_dir, c(
    "theoretical_state_pms_vs_o2.tsv",
    "theoretical_state_pms_vs_o2_heatmap.tsv"
  )), force = TRUE)
  write.table(
    live_state_pms_o2_df,
    file = file.path(out_dir, "live_weighted_pms_surface_points.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    live_state_pms_o2_df,
    file = file.path(out_dir, "live_weighted_pms_vs_o2.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    live_state_pms_o2_df,
    file = file.path(out_dir, "live_state_pms_vs_o2.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  o2_bin_width <- 0.05
  live_fraction_bin_width <- 0.01
  write.table(
    surface_grid_df,
    file = file.path(out_dir, "live_weighted_pms_surface_grid.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    surface_grid_df,
    file = file.path(out_dir, "live_weighted_pms_vs_o2_heatmap.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    surface_grid_df,
    file = file.path(out_dir, "live_state_pms_vs_o2_heatmap.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p_burden_vs_o2 <- ggplot(o2_burden_df, aes(x = burden_mm3, y = o2_pct, color = cohort, group = sample_id)) +
    geom_path(linewidth = 0.75, alpha = 0.9) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_x") +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    coord_cartesian(ylim = c(o2_plot_min, o2_plot_max)) +
    labs(
      title = "O2 Supply-Demand MAP Model: Predicted Oxygen vs Burden",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Tumor burden (mm^3)",
      y = "Oxygen (%)",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  surface_plot_df <- surface_grid_df %>%
    filter(inside_support, is.finite(live_weighted_effective_p_ms_surface))
  if (nrow(surface_plot_df) > 0L) {
    p_live_state_pms_o2 <- ggplot(
      surface_plot_df,
      aes(x = o2_pct, y = total_live_fraction, fill = live_weighted_effective_p_ms_surface)
    ) +
      geom_tile(width = o2_bin_width, height = live_fraction_bin_width) +
      facet_wrap(~ cohort, nrow = 1) +
      coord_cartesian(
        xlim = c(o2_plot_min, o2_plot_max),
        ylim = c(0, 1)
      ) +
      scale_fill_gradient(
        low = "#2C7BB6",
        high = "#F28E2B",
        limits = c(
          0,
          max(surface_grid_df$live_weighted_effective_p_ms_surface, na.rm = TRUE)
        ),
        name = "Live-weighted\np_ms"
      ) +
      labs(
        title = "Scenario-Level Surface: Oxygen Supply, Total Live Fraction, and Live-Weighted p_ms",
        subtitle = "Interpolated from model-predicted sample-day points under the fitted seed parameters. Shading is shown only within trajectory-supported regions for the 2N-start and 4N-start cohort scenarios.",
        x = "Oxygen (%)",
        y = "Total live fraction"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "grey80"),
        aspect.ratio = 1
      )
  } else {
    p_live_state_pms_o2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient support for surface interpolation") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      theme_void()
  }

  p_ploidy_heatmap <- ggplot(ploidy_all, aes(x = day, y = N, fill = fraction)) +
    geom_raster(interpolate = FALSE) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(min(ploidy_all$N, na.rm = TRUE), 100)) +
    scale_fill_gradientn(
      colours = viridisLite::viridis(3, option = "C"),
      values = c(0, 0.1, 1),
      limits = c(0, 1)
    ) +
    labs(
      title = if (identical(assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy")), "chr_number")) {
        "O2 Supply-Demand MAP Model: Predicted Chromosome-State Distribution Over Time"
      } else {
        "O2 Supply-Demand MAP Model: Predicted Ploidy Distribution Over Time"
      },
      subtitle = "Heatmap of fraction by chromosome number (N)",
      x = "Day",
      y = "Chromosome Number (N)",
      fill = "Fraction"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey84", linewidth = 0.06),
      panel.grid.minor = element_blank(),
      panel.ontop = TRUE,
      panel.background = element_rect(fill = NA, color = NA)
    )

  top_states <- ploidy_all %>%
    group_by(N) %>%
    summarise(mean_fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_fraction)) %>%
    slice_head(n = top_n) %>%
    pull(N)
  ploidy_top <- ploidy_all %>%
    filter(N %in% top_states) %>%
    mutate(N = factor(N, levels = top_states))

  p_ploidy_lines <- ggplot(ploidy_top, aes(x = day, y = fraction, color = N)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ harvest, ncol = 2) +
    labs(
      title = if (identical(assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy")), "chr_number")) {
        paste0("O2 Supply-Demand MAP Model: Chromosome-State Fractions Over Time (Top ", top_n, " N States)")
      } else {
        paste0("O2 Supply-Demand MAP Model: Ploidy Over Time (Top ", top_n, " N States)")
      },
      x = "Day",
      y = "Fraction",
      color = "N"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_weighted_mean <- ggplot(ploidy_mean, aes(x = day, y = weighted_mean_ploidy)) +
    geom_line(color = "#d62728", linewidth = 0.9) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(min(ploidy_mean$weighted_mean_ploidy, na.rm = TRUE), max(ploidy_mean$weighted_mean_ploidy, na.rm = TRUE))) +
    labs(
      title = paste0("O2 Supply-Demand MAP Model: ", weighted_mean_series_label(cfg), " Over Time"),
      subtitle = "Weighted by predicted viable-state fractions",
      x = "Day",
      y = weighted_mean_series_label(cfg)
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "burden_trend.pdf"), p_burden, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_trend_absolute.pdf"), p_burden_abs, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_trend_absolute(real_scale).pdf"), p_burden_abs_real, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_live_dead_decomposition.pdf"), p_burden_decomp, width = 13, height = 9)
  if (exists("p_o2_lag", inherits = FALSE)) ggsave(file.path(out_dir, "o2_target_vs_eff_timecourse.pdf"), p_o2_lag, width = 13, height = 9)
  if (exists("p_o2_lag_gap", inherits = FALSE)) ggsave(file.path(out_dir, "o2_lag_gap_timecourse.pdf"), p_o2_lag_gap, width = 13, height = 9)
  ggsave(file.path(out_dir, "predict_burden_vs_o2.pdf"), p_burden_vs_o2, width = 13, height = 9)
  ggsave(file.path(out_dir, "oxygen_vs_live_state_pms_colored_by_live_fraction.pdf"), p_live_state_pms_o2, width = 12, height = 6)
  ggsave(file.path(out_dir, "oxygen_livefraction_liveweighted_pms_surface.pdf"), p_live_state_pms_o2, width = 12, height = 6)
  ggsave(file.path(out_dir, "ploidy_heatmap_over_time.pdf"), p_ploidy_heatmap, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_top_states_over_time.pdf"), p_ploidy_lines, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_weighted_mean_over_time.pdf"), p_ploidy_weighted_mean, width = 13, height = 9)
  functional_plots <- plot_functional_response_curves(
    run_params = run_params,
    cfg = cfg,
    out_dir = out_dir,
    fixed_glucose_specs = fixed_glucose_specs
  )
  legend_inside <- function(p, x = 0.98, y = 0.98) {
    p + theme(
      legend.position = c(x, y),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.7), color = "grey70"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.width = grid::unit(0.50, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )
  }
  if (exists("p_o2_lag", inherits = FALSE) &&
      is.list(functional_plots) &&
      all(c("p_msr_buffer_death_per_division", "p_prolif", "p_death") %in% names(functional_plots))) {
    p_o2_panel <- p_o2_lag +
      labs(
        title = "Oxygen Evolution Over Time",
        subtitle = NULL
      )
    p_o2_panel <- legend_inside(p_o2_panel, x = 0.98, y = 0.98)
    p_msr_buffer_death_panel <- legend_inside(functional_plots$p_msr_buffer_death_per_division, x = 0.98, y = 0.98)
    p_prolif_panel <- legend_inside(functional_plots$p_prolif, x = 0.98, y = 0.98)
    p_death_panel <- legend_inside(functional_plots$p_death, x = 0.98, y = 0.98)
    grDevices::pdf(file.path(out_dir, "oxygen_response_4panel.pdf"), width = 18, height = 12, onefile = TRUE)
    grid::grid.newpage()
    lay <- grid::grid.layout(nrow = 2, ncol = 2)
    grid::pushViewport(grid::viewport(layout = lay))
    print(p_o2_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_msr_buffer_death_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(p_prolif_panel, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(p_death_panel, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
    grDevices::dev.off()
  }
  predict_horizons <- as_num_vec(argv$predict_horizons, c(100, 300, 1000))
  predict_horizons <- sort(unique(predict_horizons[is.finite(predict_horizons) & predict_horizons > 0]))
  predict_report_dt <- as_num(argv$predict_report_dt, report_dt)
  if (!is.finite(predict_report_dt) || predict_report_dt <= 0) predict_report_dt <- report_dt
  predict_top_n <- as_int(argv$predict_top_n, top_n)
  if (!is.finite(predict_top_n) || predict_top_n < 1) predict_top_n <- top_n
  do_predict_plots <- as_bool(argv$predict_plots, TRUE)
  overview_predict_horizon_day <- 1000L
  p_predict_for_overview <- NULL
  p_o2_1000_for_overview <- NULL
  p_burden_decomp_predict_for_overview <- NULL
  predict_results <- list()
  unlink(file.path(out_dir, "predict_burden_live_dead_decomposition_combined.pdf"), force = TRUE)

  if (isTRUE(do_predict_plots) && length(predict_horizons) > 0) {
    for (hz in predict_horizons) {
      message("  Predict plots: 0-", hz, " days (report_dt=", predict_report_dt, ")")
      p_hz <- plot_predict_horizon(
        run_params = run_params,
        scenarios = scenarios,
        cfg = cfg,
        out_dir = out_dir,
        horizon_day = hz,
        report_dt = predict_report_dt,
        top_n = predict_top_n
      )
      if (is.list(p_hz)) {
        predict_results[[length(predict_results) + 1L]] <- p_hz
      }
      hz_int <- as.integer(round(hz))
      p_predict_hz <- if (is.list(p_hz)) p_hz$p_predict else NULL
      p_o2_hz <- if (is.list(p_hz)) p_hz$p_o2_time else NULL
      p_burden_decomp_hz <- if (is.list(p_hz)) p_hz$p_burden_decomp_predict else NULL
      if (isTRUE(hz_int == overview_predict_horizon_day) || (is.null(p_predict_for_overview) && isTRUE(hz_int == as.integer(round(max(predict_horizons)))))) {
        p_predict_for_overview <- p_predict_hz
      }
      if (isTRUE(hz_int == overview_predict_horizon_day) || (is.null(p_o2_1000_for_overview) && isTRUE(hz_int == as.integer(round(max(predict_horizons)))))) {
        p_o2_1000_for_overview <- p_o2_hz
      }
      if (isTRUE(hz_int == overview_predict_horizon_day) || (is.null(p_burden_decomp_predict_for_overview) && isTRUE(hz_int == as.integer(round(max(predict_horizons)))))) {
        p_burden_decomp_predict_for_overview <- p_burden_decomp_hz
      }
    }
    plot_predict_burden_live_dead_decomposition_combined(
      predict_results = predict_results,
      out_dir = out_dir,
      death_language = death_language
    )
  }

  if (is.list(functional_plots) &&
      all(c("p_msr_buffer_death_per_division", "p_prolif", "p_death") %in% names(functional_plots)) &&
      inherits(p_burden_abs_real, "ggplot") &&
      inherits(p_burden_decomp, "ggplot") &&
      inherits(p_ploidy_weighted_mean, "ggplot") &&
      inherits(p_o2_1000_for_overview, "ggplot") &&
      inherits(p_predict_for_overview, "ggplot") &&
      inherits(p_burden_decomp_predict_for_overview, "ggplot")) {
    p_row1_left <- p_ploidy_weighted_mean +
      labs(title = paste0(weighted_mean_series_label(cfg), " Over Time"), subtitle = NULL) +
      theme(legend.position = "none")
    p_row1_right <- p_burden_abs_real +
      labs(title = "Burden Trend Absolute (Real Scale)", subtitle = NULL) +
      theme(legend.position = "none")
    p_row2_col1 <- legend_inside(functional_plots$p_msr_buffer_death_per_division, x = 0.98, y = 0.98)
    p_row2_col2 <- legend_inside(functional_plots$p_prolif, x = 0.98, y = 0.98)
    p_row2_col3 <- legend_inside(functional_plots$p_death, x = 0.98, y = 0.98)
    p_row3_left <- p_burden_decomp +
      labs(title = "Burden Live/Dead Decomposition", subtitle = NULL)
    p_row3_left <- legend_inside(p_row3_left, x = 0.98, y = 0.98)
    p_row3_right <- p_o2_1000_for_overview +
      labs(title = "Oxygen Evolution Over Time (0-1000 days)", subtitle = NULL)
    p_row3_right <- legend_inside(p_row3_right, x = 0.98, y = 0.98)
    p_row4_left <- p_predict_for_overview +
      labs(title = "Predict Curves (0-1000 days)", subtitle = NULL)
    p_row4_left <- legend_inside(p_row4_left, x = 0.98, y = 0.98)
    p_row4_right <- p_burden_decomp_predict_for_overview +
      labs(title = "Predict Burden Live/Dead Decomposition (0-1000 days)", subtitle = NULL)
    p_row4_right <- legend_inside(p_row4_right, x = 0.98, y = 0.98)

    grDevices::pdf(file.path(out_dir, "overview_9panel.pdf"), width = 20, height = 30, onefile = TRUE)
    grid::grid.newpage()
    # 4-row layout with variable column count: row1=2, row2=3, row3=2, row4=2.
    # Use a 6-column grid so rows can span full width consistently.
    lay <- grid::grid.layout(nrow = 4, ncol = 6)
    grid::pushViewport(grid::viewport(layout = lay))
    print(p_row1_left,  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3))
    print(p_row1_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 4:6))
    print(p_row2_col1,  vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    print(p_row2_col2,  vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 3:4))
    print(p_row2_col3,  vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 5:6))
    print(p_row3_left,  vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1:3))
    print(p_row3_right, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 4:6))
    print(p_row4_left,  vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:3))
    print(p_row4_right, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 4:6))
    grDevices::dev.off()
  }

  if (isTRUE(joint_fit)) {
    cleanup_joint_viz_root(viz_root)
  }
  normalizePath(out_dir)
}

# -----------------------------------------------------------------------------
# Function: main
# Purpose: Entry point: parse options, run fitting/visualization workflow, and write outputs.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
main <- function() {
  script_dir <- get_script_dir_self()
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  model_path <- file.path(workflow_root, "model", "model_O2G_supply_demand_MAP.R")
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
  source(model_path)

  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  results_root <- normalizePath(file.path(workflow_root, "..", "..", "results"), mustWork = FALSE)
  fit_root <- if (!is.null(argv$fit_dir)) {
    normalizePath(argv$fit_dir, mustWork = TRUE)
  } else {
    normalizePath(find_latest_fit_dir(results_root), mustWork = TRUE)
  }

  report_dt <- as_num(argv$report_dt, 1.0)
  if (!is.finite(report_dt) || report_dt <= 0) stop("report_dt must be > 0")
  top_n <- as_int(argv$top_n, 6L)
  if (!is.finite(top_n) || top_n < 1) stop("top_n must be >= 1")

  data_dir <- if (!is.null(argv$data_dir)) {
    argv$data_dir
  } else {
    normalizePath(file.path(script_dir, "..", "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
  }
  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- resolve_terminal_ploidy_path(data_dir)

  fit_dirs <- find_fit_dirs_under(fit_root)
  if (length(fit_dirs) == 0) {
    stop(
      "No valid fit results found. Need either: ",
      "(a) fit subdirectories under fit_dir, or ",
      "(b) fit_config.rds + best_params.tsv directly under fit_dir. fit_dir=", fit_root
    )
  }

  fit_dirs <- sort(unique(fit_dirs))
  message("Found ", length(fit_dirs), " fit directories under: ", fit_root)
  out_dir_override <- argv$out_dir %||% NULL
  if (!is.null(out_dir_override) && length(fit_dirs) > 1L) {
    stop("--out_dir can only be used when --fit_dir resolves to a single fit directory.")
  }

  n_cores <- as_int(argv$n_cores, 1L)
  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
  n_workers <- as.integer(max(1L, min(length(fit_dirs), n_cores)))

  process_one <- function(i) {
    fit_dir <- fit_dirs[[i]]
    message("[", i, "/", length(fit_dirs), "] Processing: ", fit_dir)
    tryCatch(
      {
        out_dir <- run_viz_for_fit_dir(
          fit_dir = fit_dir,
          argv = argv,
          dt_path = dt_path,
          ploidy_path = ploidy_path,
          report_dt = report_dt,
          top_n = top_n,
          out_dir_override = out_dir_override
        )
        message("  Done: ", out_dir)
        list(ok = TRUE, fit_dir = fit_dir, out_dir = out_dir, error = NA_character_)
      },
      error = function(e) {
        err <- conditionMessage(e)
        message("  Failed: ", err)
        list(ok = FALSE, fit_dir = fit_dir, out_dir = NA_character_, error = err)
      }
    )
  }

  use_mc <- (n_workers > 1L) && (.Platform$OS.type != "windows")
  if (use_mc) {
    message("Visualization parallel mode enabled: workers=", n_workers)
    res <- parallel::mclapply(
      X = seq_along(fit_dirs),
      FUN = process_one,
      mc.cores = n_workers,
      mc.preschedule = FALSE
    )
  } else {
    if (n_workers > 1L && .Platform$OS.type == "windows") {
      message("Windows platform detected; visualization falls back to serial.")
    }
    res <- lapply(seq_along(fit_dirs), process_one)
  }

  ok_idx <- which(vapply(res, function(x) isTRUE(x$ok), logical(1)))
  bad_idx <- setdiff(seq_along(res), ok_idx)
  ok <- vapply(res[ok_idx], function(x) x$fit_dir, character(1))
  failed <- vapply(res[bad_idx], function(x) paste0(x$fit_dir, " :: ", x$error), character(1))

  if (length(ok_idx) == 0) {
    stop("All fit subdirectories failed.")
  }

  message("Finished. Success: ", length(ok), " | Failed: ", length(failed))
  if (length(failed) > 0) {
    message("Failed directories:")
    for (x in failed) message("  - ", x)
  }
}

if (sys.nframe() == 0) {
  main()
}
