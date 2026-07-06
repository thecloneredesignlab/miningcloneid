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
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_parameter_plot.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
get_script_dir_self <- o2sd_get_script_dir
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool
clip <- o2sd_clip

ploidy_fraction_fill_scale <- function(fill_max, name = "Fraction") {
  fill_max <- suppressWarnings(as.numeric(fill_max))
  if (!is.finite(fill_max) || fill_max <= 0) fill_max <- 1
  ggplot2::scale_fill_gradientn(
    colors = c("#f7f7f7", "#2c7fb8", "#ffff33"),
    values = scales::rescale(c(0, 0.05 * fill_max, fill_max), from = c(0, fill_max)),
    limits = c(0, fill_max),
    oob = scales::squish,
    na.value = "white",
    name = name
  )
}

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

make_reference_ploidy_states <- function(cfg, min_multiple = 1.5, max_multiple = 5.0, step = 0.5) {
  n_unit <- suppressWarnings(as.numeric(cfg$N_UNIT))
  if (!is.finite(n_unit) || n_unit <= 0) n_unit <- 22.0
  ref_state_mult <- seq(min_multiple, max_multiple, by = step)
  ref_state_label <- ifelse(
    abs(ref_state_mult - round(ref_state_mult)) < 1e-8,
    paste0(as.integer(round(ref_state_mult)), "N"),
    paste0(format(ref_state_mult, trim = TRUE, nsmall = 1), "N")
  )
  data.frame(
    cohort = ref_state_label,
    ploidy_multiple = ref_state_mult,
    N_ref = as.numeric(ref_state_mult * n_unit),
    stringsAsFactors = FALSE
  )
}

predict_weighted_metric_label <- function(cfg) {
  mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  if (identical(mode, "chr_number")) "Weighted mean chromosome number" else "Weighted mean ploidy"
}

resource_death_language <- function() {
  list(
    live_label = "Viable",
    component_label = "Hypoxia-Origin Dead",
    cin_component_label = "CIN-Associated Dead",
    adjective = "hypoxia-origin",
    figure_phrase = "hypoxia-origin dead",
    report_phrase = "hypoxia-origin dead"
  )
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
  needed_common <- c(
    "lam_max", "p_misseg", "k_o_mis",
    "o2_S0", "kappa_O", "eta_o2",
    "alpha_o2", "gamma_growth",
    "mu_hp", "gamma_mu", "O2_crit", "n_O", "k_clear"
  )
  needed <- c(
    needed_common,
    "buffer_smax", "buffer_beta", "buffer_n_exp",
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
  o2_min_use <- as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0.0))
  if (!is.finite(o2_min_use) || o2_min_use < 0) o2_min_use <- 0.0
  o2_min_use <- min(max(o2_min_use, 0), o2_s0_upper_use)
  tau_O2_use <- as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 0)
  o2_growth_use <- isTRUE(.first_non_null_local(cfg$O2_growth, TRUE))
  death_language <- resource_death_language()
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
    o2_S0_upper_bound = as.numeric(o2_s0_upper_use),
    eta_o2 = as.numeric(eta_o2_use),
    o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01)),
    o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct, 0.005)),
    o2_cache_profile = isTRUE(.first_non_null_local(cfg$o2_cache_profile, FALSE)),
    O2_growth = isTRUE(o2_growth_use),
    lam_max = as.numeric(run_params$lam_max),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0.0)),
    boundary = as.character(.first_non_null_local(run_params$boundary, "drop")),
    eps_tail = as.numeric(1e-8),
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
    pred_o2_pct = as.numeric(.first_non_null_local(sim_cpp$O2_eff_obs, rep(NA_real_, length(obs_steps_cpp))))
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

make_canonical_initial_cohort_scenarios <- function(horizon_day, cfg) {
  horizon_day <- as.numeric(horizon_day)
  if (!is.finite(horizon_day) || horizon_day < 0) horizon_day <- 0
  dt_use <- as.numeric(.first_non_null_local(cfg$DT, 1.0))
  if (!is.finite(dt_use) || dt_use <= 0) dt_use <- 1.0
  lapply(c("2N", "4N"), function(cohort_i) {
    list(
      harvest = paste0("canonical_", cohort_i),
      cohort = cohort_i,
      dose = 0,
      treat_day = horizon_day + dt_use,
      sim_end_day = horizon_day,
      obs_days = numeric(0),
      obs_burden = numeric(0)
    )
  })
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

log10_plot_floor <- function(x, default = 1e-12) {
  x <- suppressWarnings(as.numeric(x))
  positive <- x[is.finite(x) & x > 0]
  if (length(positive) > 0L) {
    floor_use <- min(positive, na.rm = TRUE) / 10
    if (is.finite(floor_use) && floor_use > 0) return(floor_use)
  }
  default <- suppressWarnings(as.numeric(default))
  if (!is.finite(default) || default <= 0) default <- 1e-12
  default
}

floor_for_log10_plot <- function(x, floor) {
  x <- suppressWarnings(as.numeric(x))
  floor <- suppressWarnings(as.numeric(floor))
  if (!is.finite(floor) || floor <= 0) floor <- 1e-12
  out <- x
  finite <- is.finite(out)
  out[finite & out <= floor] <- floor
  out[!finite] <- NA_real_
  out
}

make_burden_decomp_long <- function(burden_decomp, death_language) {
  burden_decomp %>%
    pivot_longer(
      cols = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
      names_to = "component",
      values_to = "value"
    ) %>%
    mutate(
      component = factor(
        component,
        levels = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
        labels = c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
      )
    )
}

make_burden_decomp_ribbon <- function(burden_decomp, death_language, floor) {
  component_levels <- c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
  burden_decomp %>%
    mutate(
      burden_live = pmax(as.numeric(burden_live), 0),
      burden_dead_hypoxia = pmax(as.numeric(burden_dead_hypoxia), 0),
      burden_dead_buffer = pmax(as.numeric(burden_dead_buffer), 0)
    ) %>%
    pivot_longer(
      cols = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
      names_to = "component_raw",
      values_to = "value"
    ) %>%
    mutate(
      component = factor(
        component_raw,
        levels = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
        labels = component_levels
      ),
      component_index = as.integer(component)
    ) %>%
    arrange(cohort, day, component_index) %>%
    group_by(cohort, day) %>%
    mutate(
      ymax_raw = cumsum(value),
      ymin_raw = ymax_raw - value,
      ymin = floor_for_log10_plot(ymin_raw, floor),
      ymax = floor_for_log10_plot(ymax_raw, floor)
    ) %>%
    ungroup() %>%
    mutate(component = factor(component, levels = component_levels))
}

format_original_scale_labels <- function(x) {
  vapply(as.numeric(x), function(z) {
    if (!is.finite(z)) return("")
    if (z != 0 && (abs(z) >= 1e4 || abs(z) < 1e-2)) {
      return(formatC(z, format = "e", digits = 1))
    }
    format(signif(z, 3), trim = TRUE, scientific = FALSE)
  }, character(1))
}

format_log10_axis_labels <- function(x) {
  vapply(log10(as.numeric(x)), function(z) {
    if (!is.finite(z)) return("")
    if (abs(z - round(z)) < 1e-8) return(as.character(as.integer(round(z))))
    format(round(z, 2), trim = TRUE, scientific = FALSE)
  }, character(1))
}

log10_burden_y_scale <- function() {
  scale_y_log10(labels = format_log10_axis_labels)
}

log10_original_breaks <- function(x, floor, n = 4) {
  x <- suppressWarnings(as.numeric(x))
  floor <- suppressWarnings(as.numeric(floor))
  positive <- x[is.finite(x) & x > 0]
  if (!length(positive)) return(floor)
  span <- c(floor, max(positive, na.rm = TRUE))
  span <- span[is.finite(span) & span > 0]
  if (length(span) < 2L) span <- c(min(span), min(span) * 10)
  exponents <- pretty(log10(span), n = n)
  out <- 10^exponents
  out[is.finite(out) & out > 0]
}

cohort_strip_layers <- function(horizon_day, y, height, text_size = 3) {
  horizon_day <- suppressWarnings(as.numeric(horizon_day))
  if (!is.finite(horizon_day) || horizon_day <= 0) horizon_day <- 100
  strip_width <- max(horizon_day * 0.018, 1)
  labels <- c("2N", "4N")
  y_raw <- suppressWarnings(as.numeric(y))
  if (!is.null(names(y))) {
    y_use <- y_raw[match(labels, names(y))]
  } else {
    y_use <- rep(y_raw, length.out = length(labels))
  }
  y_use[!is.finite(y_use)] <- 1
  strip_df <- data.frame(
    cohort = factor(labels, levels = labels),
    x = horizon_day + strip_width / 2,
    y = y_use,
    label = labels,
    stringsAsFactors = FALSE
  )
  list(
    geom_tile(
      data = strip_df[strip_df$label == "2N", , drop = FALSE],
      aes(x = x, y = y),
      inherit.aes = FALSE,
      width = strip_width,
      height = height,
      fill = "#9ecae1",
      color = "grey70",
      linewidth = 0.3
    ),
    geom_tile(
      data = strip_df[strip_df$label == "4N", , drop = FALSE],
      aes(x = x, y = y),
      inherit.aes = FALSE,
      width = strip_width,
      height = height,
      fill = "lightpink",
      color = "grey70",
      linewidth = 0.3
    ),
    geom_text(
      data = strip_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      angle = 270,
      size = text_size,
      color = "black"
    )
  )
}

cohort_facet_grid <- function() {
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Package 'ggh4x' is required for colored 2N/4N facet strips in Predicted plots.", call. = FALSE)
  }
  ggh4x::facet_grid2(
    rows = vars(cohort),
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = c("#9ecae1", "lightpink"),
        colour = "grey70"
      ),
      text_y = ggh4x::elem_list_text(
        colour = "black",
        size = 8
      )
    )
  )
}

make_predict_annotation_track_plot <- function(
  df,
  value_col,
  y_label,
  legend_title,
  day_width,
  horizon_day,
  x_breaks,
  colors,
  transform = "identity",
  limits = NULL,
  breaks = NULL,
  labels = NULL,
  title = NULL,
  subtitle = NULL,
  show_legend = TRUE
) {
  plot_df <- data.frame(
    cohort = factor(as.character(df$cohort), levels = c("2N", "4N")),
    day = as.numeric(df$day),
    value = suppressWarnings(as.numeric(df[[value_col]])),
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[is.finite(plot_df$day), , drop = FALSE]
  plot_df$cohort_y <- ifelse(as.character(plot_df$cohort) == "2N", 2, 1)
  if (identical(transform, "log10")) {
    floor <- log10_plot_floor(plot_df$value, default = 1e-12)
    original_breaks <- log10_original_breaks(plot_df$value, floor = floor, n = 4)
    if (length(original_breaks) > 2L) {
      original_breaks <- original_breaks[c(1L, length(original_breaks))]
    }
    plot_df$value_fill <- log10(floor_for_log10_plot(plot_df$value, floor))
    fill_scale <- scale_fill_gradientn(
      colors = colors,
      breaks = log10(original_breaks),
      labels = format_original_scale_labels(original_breaks),
      na.value = "grey90",
      name = legend_title,
      guide = guide_colorbar(
        barheight = grid::unit(0.20, "in"),
        barwidth = grid::unit(0.12, "in")
      )
    )
  } else {
    plot_df$value_fill <- plot_df$value
    fill_scale <- scale_fill_gradientn(
      colors = colors,
      limits = limits,
      breaks = breaks,
      labels = labels,
      oob = scales::squish,
      na.value = "grey90",
      name = legend_title,
      guide = guide_colorbar(
        barheight = grid::unit(0.20, "in"),
        barwidth = grid::unit(0.12, "in")
      )
    )
  }

  ggplot(plot_df, aes(x = day, y = cohort_y, fill = value_fill)) +
    geom_tile(width = day_width, height = 0.9) +
    cohort_strip_layers(horizon_day = horizon_day, y = c("2N" = 2, "4N" = 1), height = 0.9, text_size = 2.8) +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, limits = c(0.5, 2.5), expand = c(0, 0)) +
    fill_scale +
    coord_cartesian(xlim = c(0, horizon_day), expand = FALSE, clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = y_label
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8),
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 5.5),
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 8),
      plot.margin = margin(1, 10, 1, 1)
    )
}

make_predicted_annotation_legend <- function(annotation_summary, o2_plot_min, o2_plot_max) {
  format_log10_labels <- function(x) {
    vapply(as.numeric(x), function(z) {
      if (!is.finite(z)) return("")
      if (abs(z - round(z)) < 1e-8) return(as.character(as.integer(round(z))))
      format(signif(z, 3), trim = TRUE, scientific = FALSE)
    }, character(1))
  }
  log10_legend_spec <- function(x, default = 1e-12, n = 5) {
    floor <- log10_plot_floor(x, default = default)
    values <- log10(floor_for_log10_plot(x, floor))
    values <- values[is.finite(values)]
    if (!length(values)) values <- log10(c(floor, floor * 10))
    range_use <- range(values, na.rm = TRUE)
    if (!all(is.finite(range_use)) || diff(range_use) <= 0) {
      range_use <- range_use[1] + c(0, 1)
    }
    breaks <- pretty(range_use, n = n)
    breaks <- breaks[is.finite(breaks) & breaks >= range_use[1] & breaks <= range_use[2]]
    if (length(breaks) < 2L) breaks <- range_use
    list(range = range_use, breaks = breaks, labels = format_log10_labels(breaks))
  }
  linear_legend_spec <- function(range_use, n = 5) {
    range_use <- suppressWarnings(as.numeric(range_use))
    range_use <- range_use[is.finite(range_use)]
    if (length(range_use) < 2L) range_use <- c(0, 1)
    range_use <- range(range_use)
    if (diff(range_use) <= 0) range_use <- range_use[1] + c(0, 1)
    breaks <- pretty(range_use, n = n)
    breaks <- breaks[is.finite(breaks) & breaks >= range_use[1] & breaks <= range_use[2]]
    if (length(breaks) < 2L) breaks <- range_use
    list(range = range_use, breaks = breaks, labels = format(signif(breaks, 3), trim = TRUE))
  }
  make_bar_df <- function(x0, x1, colors, n = 80L) {
    data.frame(
      x = seq(x0, x1, length.out = n),
      y = 0.42,
      fill_col = grDevices::colorRampPalette(colors)(n),
      stringsAsFactors = FALSE
    )
  }
  make_tick_df <- function(x0, x1, spec) {
    rng <- spec$range
    pos <- if (diff(rng) > 0) {
      x0 + (spec$breaks - rng[1]) / diff(rng) * (x1 - x0)
    } else {
      rep(mean(c(x0, x1)), length(spec$breaks))
    }
    data.frame(
      x = pos,
      label = spec$labels,
      stringsAsFactors = FALSE
    )
  }
  burden_spec <- log10_legend_spec(annotation_summary$burden, default = 1e-12, n = 5)
  live_spec <- log10_legend_spec(annotation_summary$live_cells, default = 1, n = 5)
  o2_spec <- linear_legend_spec(c(o2_plot_min, o2_plot_max), n = 5)
  burden_ticks <- make_tick_df(0.02, 0.28, burden_spec)
  live_ticks <- make_tick_df(0.37, 0.63, live_spec)
  o2_ticks <- make_tick_df(0.72, 0.98, o2_spec)
  tick_df <- dplyr::bind_rows(
    dplyr::mutate(burden_ticks, track = "burden"),
    dplyr::mutate(live_ticks, track = "live"),
    dplyr::mutate(o2_ticks, track = "o2")
  )
  bars <- dplyr::bind_rows(
    make_bar_df(0.02, 0.28, c("#ffffbf", "#542788")),
    make_bar_df(0.37, 0.63, c("#ffffbf", "#2c7bb6")),
    make_bar_df(0.72, 0.98, c("#f7f7f7", "#9ecae1", "#08519c"))
  )
  ggplot(bars, aes(x = x, y = y, fill = fill_col)) +
    geom_tile(width = 0.004, height = 0.22) +
    scale_fill_identity() +
    geom_segment(
      data = tick_df,
      aes(x = x, xend = x, y = 0.27, yend = 0.33),
      inherit.aes = FALSE,
      linewidth = 0.18,
      color = "grey20"
    ) +
    geom_text(
      data = tick_df,
      aes(x = x, y = 0.12, label = label),
      inherit.aes = FALSE,
      size = 1.9,
      color = "black"
    ) +
    annotate("text", x = 0.02, y = 0.82, label = "log10 Burden (mm^3)", hjust = 0, size = 2.6) +
    annotate("text", x = 0.37, y = 0.82, label = "log10 Viable cells", hjust = 0, size = 2.6) +
    annotate("text", x = 0.72, y = 0.82, label = "Effective O2 (%)", hjust = 0, size = 2.6) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
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

compute_missegregation_probability_timecourse <- function(ploidy_all, burden_all, run_params) {
  empty <- data.frame(
    harvest = character(),
    cohort = character(),
    dose = numeric(),
    day = numeric(),
    sample_id = character(),
    o2_eff_pct = numeric(),
    mean_p_misseg = numeric(),
    weighted_mean_N = numeric(),
    total_fraction = numeric(),
    stringsAsFactors = FALSE
  )
  required_ploidy <- c("harvest", "cohort", "dose", "day", "N", "fraction")
  required_burden <- c("harvest", "cohort", "dose", "day", "pred_o2_pct")
  if (!all(required_ploidy %in% names(ploidy_all)) || !all(required_burden %in% names(burden_all))) {
    return(empty)
  }

  ploidy_df <- ploidy_all %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      day_key = round(as.numeric(day), 8),
      N = as.numeric(N),
      fraction = pmax(as.numeric(fraction), 0)
    ) %>%
    filter(is.finite(day_key), is.finite(N), is.finite(fraction), fraction > 0)

  o2_df <- burden_all %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      day_key = round(as.numeric(day), 8),
      o2_eff_pct = as.numeric(pred_o2_pct)
    ) %>%
    filter(is.finite(day_key), is.finite(o2_eff_pct)) %>%
    distinct(harvest, cohort, dose, day_key, .keep_all = TRUE)

  joined <- inner_join(
    ploidy_df,
    o2_df[, c("harvest", "cohort", "dose", "day_key", "o2_eff_pct"), drop = FALSE],
    by = c("harvest", "cohort", "dose", "day_key")
  )
  if (!nrow(joined)) return(empty)

  joined$p_misseg <- pmax(pmin(as.numeric(.pmisseg_of_O2(
    O2 = joined$o2_eff_pct,
    run_params = run_params,
    N = joined$N
  )), 1), 0)
  joined <- joined[is.finite(joined$p_misseg), , drop = FALSE]
  if (!nrow(joined)) return(empty)

  out <- joined %>%
    group_by(harvest, cohort, dose, day, o2_eff_pct) %>%
    summarise(
      mean_p_misseg = sum(fraction * p_misseg, na.rm = TRUE) / pmax(sum(fraction, na.rm = TRUE), 1e-12),
      weighted_mean_N = sum(fraction * N, na.rm = TRUE) / pmax(sum(fraction, na.rm = TRUE), 1e-12),
      total_fraction = sum(fraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      sample_id = paste(as.character(harvest), as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__")
    ) %>%
    select(harvest, cohort, dose, day, sample_id, o2_eff_pct, mean_p_misseg, weighted_mean_N, total_fraction) %>%
    arrange(harvest, cohort, dose, day)

  as.data.frame(out, stringsAsFactors = FALSE)
}

plot_missegregation_probability_timecourse <- function(ms_timecourse, out_dir, fit_dir = NULL, report_dt = NULL) {
  if (is.null(ms_timecourse) || !is.data.frame(ms_timecourse) || !nrow(ms_timecourse)) {
    return(FALSE)
  }
  plot_df <- ms_timecourse %>%
    mutate(
      day = as.numeric(day),
      mean_p_misseg = as.numeric(mean_p_misseg),
      cohort = as.character(cohort),
      sample_id = as.character(sample_id)
    ) %>%
    filter(is.finite(day), is.finite(mean_p_misseg))
  if (!nrow(plot_df)) return(FALSE)

  subtitle_parts <- character()
  if (!is.null(fit_dir)) subtitle_parts <- c(subtitle_parts, paste0("fit_dir=", basename(fit_dir)))
  if (!is.null(report_dt) && is.finite(as.numeric(report_dt))) {
    subtitle_parts <- c(subtitle_parts, paste0("report_dt=", signif(as.numeric(report_dt), 4)))
  }
  subtitle <- if (length(subtitle_parts)) paste(subtitle_parts, collapse = " | ") else NULL
  p <- ggplot(plot_df, aes(x = day, y = mean_p_misseg, color = cohort, group = sample_id)) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    facet_wrap(~ harvest, ncol = 2) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 3)) +
    labs(
      title = "Resource-Stress Model: Mean Per-Chromosome Missegregation Probability Over Time",
      subtitle = subtitle,
      x = "Day",
      y = "Viable-population-weighted mean per-chromosome missegregation probability",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)

  pdf_path <- file.path(out_dir, "missegregation_probability_over_time.pdf")
  png_path <- file.path(out_dir, "missegregation_probability_over_time.png")
  ggsave(pdf_path, p, width = 13, height = 9)
  ggsave(png_path, p, width = 13, height = 9, dpi = 180, bg = "white")
  TRUE
}

parse_harvest_day_value <- function(x) {
  x_chr <- as.character(x)
  out <- suppressWarnings(as.numeric(x_chr))
  missing <- !is.finite(out)
  if (any(missing)) {
    extracted <- sub("^.*?([0-9]+(?:\\.[0-9]+)?).*$", "\\1", x_chr[missing], perl = TRUE)
    out[missing] <- suppressWarnings(as.numeric(extracted))
  }
  out
}

plot_prediction_horizon_population_average_cin <- function(ms_timecourse,
                                                           out_dir,
                                                           target_days = c(100, 300, 1000),
                                                           day_tol = 1e-6) {
  if (is.null(ms_timecourse) || !is.data.frame(ms_timecourse) || !nrow(ms_timecourse)) {
    return(FALSE)
  }
  required_cols <- c("cohort", "day", "sample_id", "mean_p_misseg")
  if (!all(required_cols %in% names(ms_timecourse))) return(FALSE)

  cohort_levels <- c("2N", "4N")
  cohort_labels <- c("2N" = "2N-derived", "4N" = "4N-derived")
  cohort_colors <- c("2N" = "#1f77b4", "4N" = "#d62728")
  cohort_linetypes <- c("2N" = "solid", "4N" = "solid")

  sample_df <- ms_timecourse %>%
    transmute(
      initial_cohort = as.character(cohort),
      day = as.numeric(day),
      sample_id = as.character(sample_id),
      population_average_cin = as.numeric(mean_p_misseg)
    ) %>%
    filter(
      initial_cohort %in% cohort_levels,
      is.finite(day),
      is.finite(population_average_cin)
    ) %>%
    distinct(initial_cohort, day, sample_id, population_average_cin, .keep_all = TRUE)
  if (!nrow(sample_df)) return(FALSE)

  multi_trajectory <- sample_df %>%
    group_by(initial_cohort) %>%
    summarise(n_trajectories = dplyr::n_distinct(sample_id), .groups = "drop") %>%
    filter(n_trajectories > 1L)
  if (nrow(multi_trajectory)) {
    warning(
      "Population-average CIN plot expected one canonical trajectory per initial cohort; ",
      "using the first trajectory for: ",
      paste(as.character(multi_trajectory$initial_cohort), collapse = ", ")
    )
    sample_df <- sample_df %>%
      arrange(initial_cohort, sample_id, day) %>%
      group_by(initial_cohort) %>%
      mutate(.first_sample_id = dplyr::first(sample_id)) %>%
      ungroup() %>%
      filter(sample_id == .first_sample_id) %>%
      select(-.first_sample_id)
  }

  target_days <- sort(unique(as.numeric(target_days[is.finite(target_days)])))
  all_rows <- list()
  any_saved <- FALSE
  for (target_day in target_days) {
    target_df <- sample_df %>%
      filter(day >= -day_tol, day <= target_day + day_tol)
    if (!nrow(target_df)) next

    plot_df <- target_df
    plot_df$target_day <- target_day
    plot_df$initial_cohort <- factor(plot_df$initial_cohort, levels = cohort_levels)
    plot_df <- plot_df %>%
      mutate(
        cohort_order = match(as.character(initial_cohort), cohort_levels),
        cohort_label = cohort_labels[as.character(initial_cohort)]
      ) %>%
      filter(is.finite(day), is.finite(population_average_cin)) %>%
      arrange(cohort_order, day)
    if (!nrow(plot_df)) next
    all_rows[[length(all_rows) + 1L]] <- plot_df

    target_label <- format(target_day, trim = TRUE, scientific = FALSE)
    target_tag <- gsub("\\.", "p", target_label)
    p <- ggplot(
      plot_df,
      aes(
        x = day,
        y = population_average_cin,
        color = initial_cohort,
        linetype = initial_cohort,
        group = initial_cohort
      )
    ) +
      geom_line(linewidth = 0.9, alpha = 0.95) +
      scale_color_manual(values = cohort_colors, labels = cohort_labels, drop = FALSE) +
      scale_linetype_manual(values = cohort_linetypes, labels = cohort_labels, drop = FALSE) +
      scale_x_continuous(limits = c(0, target_day), breaks = pretty(c(0, target_day), n = 5)) +
      scale_y_continuous(labels = function(x) formatC(x, format = "f", digits = 4)) +
      labs(
        title = paste0("0-", target_label, " Day Population-average CIN rate over time"),
        subtitle = "Canonical 2N-derived and 4N-derived trajectories simulated from the fitted parameters.",
        x = "Day",
        y = "Population-average CIN rate\n(mean per-chromosome missegregation probability)",
        color = "Initial cohort",
        linetype = "Initial cohort"
      ) +
      theme_bw(base_size = 11)

    out_base <- paste0("population_average_cin_by_initial_cohort_day", target_tag)
    ggsave(file.path(out_dir, paste0(out_base, ".pdf")), p, width = 10, height = 7)
    ggsave(file.path(out_dir, paste0(out_base, ".png")), p, width = 10, height = 7, dpi = 180, bg = "white")
    any_saved <- TRUE
  }

  if (length(all_rows)) {
    out_df <- dplyr::bind_rows(all_rows) %>%
      transmute(
        target_day = as.numeric(target_day),
        day = as.numeric(day),
        initial_cohort = as.character(initial_cohort),
        cohort_label = as.character(cohort_label),
        cohort_order = as.integer(cohort_order),
        sample_id = as.character(sample_id),
        population_average_cin = as.numeric(population_average_cin)
      ) %>%
      arrange(target_day, cohort_order, day)
    tsv_path <- file.path(out_dir, "population_average_cin_by_initial_cohort_horizons.tsv")
    if (file.exists(tsv_path)) {
      old_df <- tryCatch(
        utils::read.delim(tsv_path, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (is.data.frame(old_df) && nrow(old_df) > 0L && "target_day" %in% names(old_df)) {
        old_df$target_day <- suppressWarnings(as.numeric(old_df$target_day))
        old_df <- old_df[!(old_df$target_day %in% unique(out_df$target_day)), , drop = FALSE]
        for (nm in setdiff(names(out_df), names(old_df))) {
          old_df[[nm]] <- NA
        }
        old_df <- old_df[, names(out_df), drop = FALSE]
        out_df <- dplyr::bind_rows(old_df, out_df) %>%
          arrange(target_day, cohort_order, day)
      }
    }
    write.table(
      out_df,
      file = tsv_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  any_saved
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
        paste0("Observed ploidy is mapped to the chromosome-number grid used by the objective | fit_dir=", basename(fit_dir))
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

make_o2_crit_adaptive_grid <- function(o2_crit, grid_min = 0, grid_max = 5, base_step = 0.02, dense_n = 220L) {
  grid_min <- as.numeric(grid_min[[1L]])
  grid_max <- as.numeric(grid_max[[1L]])
  base_step <- as.numeric(base_step[[1L]])
  dense_n <- as.integer(dense_n[[1L]])
  if (!is.finite(grid_min) || !is.finite(grid_max)) stop("O2 adaptive grid bounds must be finite.")
  if (grid_max < grid_min) {
    tmp <- grid_min
    grid_min <- grid_max
    grid_max <- tmp
  }
  if (!is.finite(base_step) || base_step <= 0) base_step <- max((grid_max - grid_min) / 250, 1e-6)
  if (!is.finite(dense_n) || dense_n < 2L) dense_n <- 220L
  base_grid <- seq(grid_min, grid_max, by = base_step)
  base_grid <- c(base_grid, grid_min, grid_max)

  crit <- suppressWarnings(as.numeric(o2_crit[[1L]]))
  if (!is.finite(crit) || crit <= 0) {
    return(sort(unique(signif(base_grid[is.finite(base_grid)], 14))))
  }

  near_upper <- min(grid_max, max(base_step, 25 * crit))
  near_zero_grid <- if (near_upper > grid_min) {
    seq(grid_min, near_upper, length.out = dense_n)
  } else {
    numeric(0)
  }
  transition_lower <- max(grid_min, crit * 0.05)
  transition_upper <- min(grid_max, crit * 20)
  transition_grid <- if (transition_upper > transition_lower) {
    seq(transition_lower, transition_upper, length.out = dense_n)
  } else {
    numeric(0)
  }
  log_grid <- crit * 10^seq(-4, 2, length.out = dense_n)
  multiplier_grid <- crit * c(
    0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3,
    0.5, 0.75, 1, 1.25, 1.5, 2, 3, 5, 8, 10, 15, 20, 30, 50, 100
  )

  grid <- c(base_grid, near_zero_grid, transition_grid, log_grid, multiplier_grid)
  grid <- grid[is.finite(grid) & grid >= grid_min & grid <= grid_max]
  sort(unique(signif(grid, 14)))
}

make_o2_crit_reference_levels <- function(o2_crit, grid_min = 0, grid_max = 5, coarse_step = 0.5) {
  grid_min <- as.numeric(grid_min[[1L]])
  grid_max <- as.numeric(grid_max[[1L]])
  coarse_step <- as.numeric(coarse_step[[1L]])
  if (!is.finite(grid_min) || !is.finite(grid_max)) stop("O2 reference level bounds must be finite.")
  if (grid_max < grid_min) {
    tmp <- grid_min
    grid_min <- grid_max
    grid_max <- tmp
  }
  if (!is.finite(coarse_step) || coarse_step <= 0) coarse_step <- max((grid_max - grid_min) / 10, 1e-6)
  coarse_levels <- seq(grid_min, grid_max, by = coarse_step)
  crit <- suppressWarnings(as.numeric(o2_crit[[1L]]))
  crit_levels <- if (is.finite(crit) && crit > 0) {
    crit * c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 25)
  } else {
    numeric(0)
  }
  levels <- c(grid_min, grid_max, coarse_levels, crit_levels)
  levels <- levels[is.finite(levels) & levels >= grid_min & levels <= grid_max]
  sort(unique(signif(levels, 12)))
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
plot_functional_response_curves <- function(run_params, cfg, out_dir, ...) {
  start_with_mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
  o2_plot_min <- 0
  o2_plot_max <- 5
  O2_crit_use <- as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  o2_grid <- make_o2_crit_adaptive_grid(O2_crit_use, o2_plot_min, o2_plot_max, base_step = 0.02, dense_n = 220L)
  state_plot_min <- if (identical(start_with_mode, "chr_number")) as.numeric(cfg$N_MIN) else 0
  state_plot_max <- if (identical(start_with_mode, "chr_number")) as.numeric(cfg$N_MAX) else 10
  state_grid_dense <- if (identical(start_with_mode, "chr_number")) {
    seq(state_plot_min, state_plot_max, by = 1)
  } else {
    seq(state_plot_min, state_plot_max, by = 0.05)
  }
  o2_levels_ploidy <- make_o2_crit_reference_levels(O2_crit_use, o2_plot_min, o2_plot_max, coarse_step = 0.5)
  state_axis_label <- functional_state_axis_label(cfg)
  ref_df <- data.frame(
    cohort = c("2N", "4N"),
    N_ref = as.numeric(c(2 * cfg$N_UNIT, 4 * cfg$N_UNIT)),
    stringsAsFactors = FALSE
  )
  ref_df_multi <- make_reference_ploidy_states(cfg)
  ref_state_subtitle <- paste0("Reference states: ", paste(ref_df_multi$cohort, collapse = ", "))

  o2_growth_use <- isTRUE(.first_non_null_local(cfg$O2_growth, TRUE))
  alpha_o2_use <- pmax(0, as.numeric(run_params$alpha_o2))
  gamma_growth_use <- pmax(as.numeric(run_params$gamma_growth), 1e-12)
  gamma_mu_use <- pmax(as.numeric(run_params$gamma_mu), 1e-12)
  ploidy_O2_death_use <- assert_canonical_ploidy_o2_death_mode(
    .first_non_null_local(cfg$ploidy_O2_death, run_params$ploidy_O2_death, "diploid_NULL")
  )
  n_O <- as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0))
  if (!is.finite(n_O) || n_O < 0) stop("run_params$n_O must be finite and >= 0.")
  mu_hp_use <- pmax(as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)), 0)
  boundary_mode_use <- as.character(.first_non_null_local(run_params$boundary, "drop"))
  p_wgd_use <- as.numeric(.first_non_null_local(run_params$p_wgd, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
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
        lam_max = as.numeric(run_params$lam_max),
        p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init, 1e-5)),
        p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
        k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
        p_wgd = as.numeric(p_wgd_use),
        boundary = as.character(boundary_mode_use),
        eps_tail = as.numeric(1e-8),
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
        "misseg_nonviable_rate",
        "boundary_dropped_rate",
        "misseg_nonviable_division_prob",
        "misseg_nonviable_daughters_per_division"
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
      curve_list$misseg_nonviable_daughter_fraction <- pmax(
        pmin(0.5 * curve_list$misseg_nonviable_daughters_per_division, 1),
        0
      )
      assign(key, curve_list, envir = functional_rate_cache)
    }
    get(key, envir = functional_rate_cache, inherits = FALSE)
  }

  compute_rate_components <- function(O2, N) {
    O2_vec <- .assert_o2_pct(as.numeric(O2), label = "O2")
    N_vec <- as.numeric(N)
    n_out <- max(length(O2_vec), length(N_vec))
    if (!(length(O2_vec) %in% c(1L, n_out) &&
          length(N_vec) %in% c(1L, n_out))) {
      stop("O2 and N must have compatible lengths in compute_rate_components().")
    }
    O2_vec <- rep_len(O2_vec, n_out)
    N_vec <- rep_len(N_vec, n_out)
    proliferation_rate <- as.numeric(.lambda_eff_of_O2(
      O2 = O2_vec,
      run_params = run_params,
      N = N_vec,
      O2_crit = O2_crit_use,
      O2_growth = o2_growth_use
    ))
    death_rate <- as.numeric(.mu_eff_of_O2(
      O2 = O2_vec,
      run_params = run_params,
      N = N_vec,
      O2_crit = O2_crit_use
    ))
    dead_buffer_rate <- rep(NA_real_, n_out)
    misseg_nonviable_rate <- rep(NA_real_, n_out)
    boundary_dropped_rate <- rep(NA_real_, n_out)
    misseg_nonviable_division_prob <- rep(NA_real_, n_out)
    misseg_nonviable_daughters_per_division <- rep(NA_real_, n_out)
    misseg_nonviable_daughter_fraction <- rep(NA_real_, n_out)
    row_groups <- split(seq_len(n_out), sprintf("%.6f", O2_vec))
    for (row_idx in row_groups) {
      rate_curves <- get_functional_rate_curves(O2 = O2_vec[[row_idx[[1]]]])
      state_idx <- as.integer(round(N_vec[row_idx])) - as.integer(cfg$N_MIN) + 1L
      valid_idx <- is.finite(state_idx) & state_idx >= 1L & state_idx <= length(rate_curves$dead_buffer_rate)
      if (any(valid_idx)) {
        idx_use <- state_idx[valid_idx]
        row_use <- row_idx[valid_idx]
        dead_buffer_rate[row_use] <- rate_curves$dead_buffer_rate[idx_use]
        misseg_nonviable_rate[row_use] <- rate_curves$misseg_nonviable_rate[idx_use]
        boundary_dropped_rate[row_use] <- rate_curves$boundary_dropped_rate[idx_use]
        misseg_nonviable_division_prob[row_use] <- rate_curves$misseg_nonviable_division_prob[idx_use]
        misseg_nonviable_daughters_per_division[row_use] <- rate_curves$misseg_nonviable_daughters_per_division[idx_use]
        misseg_nonviable_daughter_fraction[row_use] <- rate_curves$misseg_nonviable_daughter_fraction[idx_use]
      }
    }
    data.frame(
      O2 = O2_vec,
      N = N_vec,
      proliferation_rate = pmax(as.numeric(proliferation_rate), 0),
      death_rate = pmax(as.numeric(death_rate), 0),
      buffer_death_rate = pmax(as.numeric(dead_buffer_rate), 0),
      buffer_death_per_division = pmax(as.numeric(dead_buffer_rate), 0) / pmax(as.numeric(proliferation_rate), 1e-12),
      misseg_nonviable_rate = pmax(as.numeric(misseg_nonviable_rate), 0),
      boundary_dropped_rate = pmax(as.numeric(boundary_dropped_rate), 0),
      misseg_nonviable_division_prob = pmax(pmin(as.numeric(misseg_nonviable_division_prob), 1), 0),
      misseg_nonviable_daughters_per_division = pmax(pmin(as.numeric(misseg_nonviable_daughters_per_division), 2), 0),
      misseg_nonviable_daughter_fraction = pmax(pmin(as.numeric(misseg_nonviable_daughter_fraction), 1), 0),
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
      misseg_nonviable_rate = rate_df$misseg_nonviable_rate,
      boundary_dropped_rate = rate_df$boundary_dropped_rate,
      misseg_nonviable_division_prob = rate_df$misseg_nonviable_division_prob,
      misseg_nonviable_daughters_per_division = rate_df$misseg_nonviable_daughters_per_division,
      misseg_nonviable_daughter_fraction = rate_df$misseg_nonviable_daughter_fraction,
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
      misseg_nonviable_rate = rate_df$misseg_nonviable_rate,
      boundary_dropped_rate = rate_df$boundary_dropped_rate,
      misseg_nonviable_division_prob = rate_df$misseg_nonviable_division_prob,
      misseg_nonviable_daughters_per_division = rate_df$misseg_nonviable_daughters_per_division,
      misseg_nonviable_daughter_fraction = rate_df$misseg_nonviable_daughter_fraction,
      net_growth_rate = rate_df$net_growth_rate,
      row.names = NULL
    )
  }))
  write.table(
    o2_curve_multi,
    file = file.path(out_dir, "functional_curve_oxygen_multi_ploidy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  o2_curve_multi$cohort <- factor(o2_curve_multi$cohort, levels = unique(o2_curve_multi$cohort))
  save_o2_curve_plot <- function(value_col, title, y_label, filename) {
    plot_df <- o2_curve_multi
    plot_df$value <- suppressWarnings(as.numeric(plot_df[[value_col]]))
    plot_df <- plot_df[is.finite(plot_df$oxygen_pct) & is.finite(plot_df$value), , drop = FALSE]
    if (!nrow(plot_df)) return(invisible(FALSE))
    p <- ggplot(plot_df, aes(x = oxygen_pct, y = value, color = cohort)) +
      geom_line(linewidth = 1) +
      labs(
        title = title,
        subtitle = "In vivo fitted rate function across oxygen levels",
        x = "Effective oxygen (%)",
        y = y_label,
        color = "Reference state"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, filename), p, width = 10, height = 7)
    invisible(TRUE)
  }
  save_o2_curve_plot(
    "ms_rate",
    "Effective Oxygen vs Missegregation Rate",
    "Missegregation rate",
    "oxygen_vs_missegregation_rate_multi_ploidy.pdf"
  )
  save_o2_curve_plot(
    "proliferation_rate",
    "Effective Oxygen vs Proliferation Rate",
    "Proliferation rate",
    "oxygen_vs_proliferation_rate.pdf"
  )
  save_o2_curve_plot(
    "death_rate",
    "Effective Oxygen vs Stress-Associated Death Rate",
    "Stress-associated death rate",
    "oxygen_vs_death_rate.pdf"
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
    "functional_curve_oxygen_g0.tsv",
    "functional_curve_oxygen_multi_ploidy_g0.tsv",
    "functional_curve_ploidy_by_o2_g0.tsv",
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
    "oxygen_vs_missegregation_rate_g0.pdf",
    "oxygen_vs_missegregation_rate_multi_ploidy_g0.pdf",
    "oxygen_vs_proliferation_rate_g0.pdf",
    "oxygen_vs_death_rate_g0.pdf",
    "ploidy_vs_proliferation_rate_by_o2_g0.pdf",
    "ploidy_vs_death_rate_by_o2_g0.pdf",
    "compare_ploidy_vs_death_rate_by_o2_g20.pdf",
    "compare_ploidy_vs_proliferation_rate_by_o2_g20.pdf",
    "compare_oxygen_vs_missegregation_rate_multi_ploidy_g20.pdf",
    "compare_oxygen_vs_proliferation_rate_g20.pdf",
    "compare_oxygen_vs_death_rate_g20.pdf",
    "compare_ploidy_vs_death_rate_by_o2_g0.pdf",
    "compare_ploidy_vs_proliferation_rate_by_o2_g0.pdf",
    "compare_oxygen_vs_missegregation_rate_multi_ploidy_g0.pdf",
    "compare_oxygen_vs_proliferation_rate_g0.pdf",
    "compare_oxygen_vs_death_rate_g0.pdf",
    "o2_g_heatmap_growth.pdf",
    "o2_g_heatmap_death.pdf",
    "o2_g_heatmap_missegregation.pdf",
    "o2_g_heatmap_growth_g0_5.pdf",
    "o2_g_heatmap_death_g0_5.pdf",
    "o2_g_heatmap_missegregation_g0_5.pdf",
    "ploidy_vs_proliferation_rate_by_o2.pdf",
    "ploidy_vs_death_rate_by_o2.pdf",
    "oxygen_vs_missegregation_rate.pdf",
    "oxygen_vs_net_growth_rate.pdf",
    "ms_rate_vs_death_rate.pdf",
    "ms_rate_vs_buffer_death_rate.pdf",
    "ms_rate_vs_buffer_death_per_division.pdf",
    "ms_rate_vs_nonviable_division_probability.pdf"
  )), force = TRUE)
  multi_colors <- stats::setNames(
    grDevices::hcl.colors(nrow(ref_df_multi), palette = "Dark 3"),
    ref_df_multi$cohort
  )

  p_death_msr <- ggplot(
    o2_curve,
    aes(x = death_rate, y = ms_rate, color = cohort, group = cohort)
  ) +
    geom_path(linewidth = 1) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = "Stress-Associated Death Rate vs Missegregation Rate",
      subtitle = "Same effective-oxygen sweep and reference chromosome-number states as Effective Oxygen vs Missegregation Rate",
      x = "Stress-associated death rate",
      y = "Missegregation rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  p_msr_nonviable_daughter_fraction <- ggplot(
    o2_curve_multi,
    aes(
      x = ms_rate,
      y = misseg_nonviable_daughter_fraction,
      color = factor(cohort, levels = ref_df_multi$cohort)
    )
  ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = multi_colors, drop = FALSE) +
    labs(
      title = "Nonviable Daughter Fraction vs Missegregation Rate Across Reference Chromosome-Number States",
      subtitle = paste0(
        "Missegregation-linked nonviable daughters / all daughters per division; excludes boundary-drop losses | ",
        ref_state_subtitle
      ),
      x = "Missegregation rate",
      y = "Nonviable daughters / all daughters",
      color = "Reference state"
    ) +
    theme_bw(base_size = 11)
  N_states <- seq.int(as.integer(cfg$N_MIN), as.integer(cfg$N_MAX))
  ploidy_grid <- N_states / as.numeric(cfg$N_UNIT)
  n_chr_use <- if (as.integer(cfg$N_UNIT) > 0L) as.numeric(cfg$N_UNIT) else 22.0
  ratio <- (2.0 * n_chr_use) / pmax(as.numeric(N_states), 1e-12)
  sN <- buffer_smax_use * exp(-buffer_beta_use * pmax(ratio, 0)^buffer_n_exp_use)
  viability <- pmax(pmin(sN, 1), 0)
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

  p_viability <- ggplot(viability_curve, aes(x = endpoint_value, y = viability_after_ms)) +
    geom_line(color = "#2ca02c", linewidth = 1) +
    labs(
      title = paste0(state_axis_label, " vs Post-Missegregation Survival"),
      subtitle = "Ploidy-dependent survival for a one-copy missegregation event",
      x = state_axis_label,
      y = "Post-missegregation survival"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "death_rate_vs_missegregation_rate.pdf"), p_death_msr, width = 10, height = 7)
  ggsave(file.path(out_dir, "ms_rate_vs_nonviable_daughter_fraction.pdf"), p_msr_nonviable_daughter_fraction, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_viability_after_ms.pdf"), p_viability, width = 10, height = 7)
  invisible(list(
    p_death_msr = p_death_msr,
    p_msr_nonviable_daughter_fraction = p_msr_nonviable_daughter_fraction,
    p_viability = p_viability
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

plot_predict_horizon <- function(run_params, scenarios, cfg, out_dir, horizon_day, report_dt = 1.0) {
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
  death_language <- resource_death_language()

  burden_all <- burden_all %>% filter(day <= horizon_day + 1e-9)
  ploidy_all <- ploidy_all %>% filter(day <= horizon_day + 1e-9)
  burden_all <- normalize_burden_for_plot(burden_all)
  ploidy_mean <- compute_ploidy_weighted_mean(ploidy_all, cfg)
  canonical_scenarios <- make_canonical_initial_cohort_scenarios(horizon_day = horizon_day, cfg = cfg)
  canonical_sim_list <- lapply(canonical_scenarios, function(sc) {
    simulate_one_full_horizon(run_params, sc, cfg, horizon_day = horizon_day, report_dt = report_dt)
  })
  canonical_burden_all <- bind_rows(lapply(canonical_sim_list, `[[`, "burden"))
  canonical_ploidy_all <- bind_rows(lapply(canonical_sim_list, `[[`, "ploidy"))
  canonical_misseg_timecourse <- compute_missegregation_probability_timecourse(
    ploidy_all = canonical_ploidy_all,
    burden_all = canonical_burden_all,
    run_params = run_params
  )
  if (nrow(canonical_misseg_timecourse) > 0L) {
    plot_prediction_horizon_population_average_cin(
      ms_timecourse = canonical_misseg_timecourse,
      out_dir = out_dir,
      target_days = as.numeric(horizon_day)
    )
  }

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
    paste0("predict_g_timecourse_", horizon_tag, ".pdf"),
    paste0("predict_live_resource_death_fraction_", horizon_tag, ".pdf"),
    paste0("predict_live_weighted_pms_", horizon_tag, ".pdf"),
    paste0("predict_death_ratio_", horizon_tag, ".pdf"),
    paste0("predict_o2_timecourse_", horizon_tag, ".pdf"),
    paste0("predicted_", horizon_tag, ".pdf")
  )), force = TRUE)

  write.table(burden_all, file = file.path(out_dir, paste0("predict_burden_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_all, file = file.path(out_dir, paste0("predict_ploidy_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ploidy_mean, file = file.path(out_dir, paste0("predict_ploidy_weighted_mean_", horizon_tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  burden_plot_df <- burden_all %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      value = as.numeric(pred_burden),
      sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
    )
  burden_plot_floor <- log10_plot_floor(burden_plot_df$value, default = 1e-12)
  burden_plot_df <- burden_plot_df %>%
    mutate(value_log_plot = floor_for_log10_plot(value, burden_plot_floor))

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
  live_cells_plot_floor <- if (nrow(live_cells_plot_df) > 0L) {
    log10_plot_floor(live_cells_plot_df$value, default = 1)
  } else {
    1
  }
  if (nrow(live_cells_plot_df) > 0L) {
    live_cells_plot_df <- live_cells_plot_df %>%
      mutate(value_log_plot = floor_for_log10_plot(value, live_cells_plot_floor))
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
    aes(x = day, y = value_log_plot, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    predict_x_scale() +
    coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
    scale_y_log10() +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = paste0("Predict Curves: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = paste0("Single summary plot (all scenarios overlaid) | fit_dir=", basename(dirname(out_dir)), " | report_dt=", report_dt),
      x = "Day",
      y = "Burden (log10 scale)",
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
    chr_density_fill_max <- max(chr_density_df$density, na.rm = TRUE)
    if (!is.finite(chr_density_fill_max) || chr_density_fill_max <= 0) {
      chr_density_fill_max <- 1
    }
    ggplot(
      chr_density_df,
      aes(x = day, y = chr_bin_mid, fill = density)
    ) +
      geom_tile(width = chr_density_day_width, height = chr_density_bin_width) +
      facet_grid(cohort ~ .) +
      predict_x_scale() +
      ploidy_fraction_fill_scale(chr_density_fill_max, name = "Probability\ndensity") +
      scale_y_continuous(
        breaks = ploidy_n_breaks,
        expand = c(0, 0),
        sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Ploidy")
      ) +
      coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
      labs(
        x = "Day",
        y = "Chromosome count (N)"
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
        y = "Chromosome count (N)"
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
        y = "Effective O2 (%)",
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
        y = "Effective O2 (%)"
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
      aes(x = day, y = value_log_plot, group = sample_id, color = cohort)
    ) +
      geom_line(linewidth = 0.65, alpha = 0.8) +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_y_log10() +
      scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
      labs(
        x = NULL,
        y = "Viable cells (log10 scale)",
        color = "Cohort"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank()
      ) +
      predict_curve_theme
  } else {
    ggplot() +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_y_log10() +
      labs(
        x = NULL,
        y = "Viable cells (log10 scale)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank()
      ) +
      predict_curve_theme
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

  resource_death_fraction_predict <- burden_decomp_predict %>%
    mutate(
      sample_id = paste(harvest, as.character(cohort), format(dose, trim = TRUE, scientific = FALSE), sep = "__"),
      resource_death_fraction = ifelse(
        burden_dead_total > 0,
        burden_dead_hypoxia / pmax(burden_dead_total, 1e-12),
        0
      ),
      resource_death_fraction = pmax(0, pmin(1, as.numeric(resource_death_fraction)))
    )
  write.table(
    resource_death_fraction_predict,
    file = file.path(out_dir, paste0("predict_resource_death_fraction_", horizon_tag, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  p_resource_death_fraction_predict <- ggplot(
    resource_death_fraction_predict,
    aes(x = day, y = resource_death_fraction, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    predict_x_scale() +
    coord_cartesian(xlim = c(0, horizon_day), ylim = c(0, 1), expand = FALSE) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      x = "Day",
      y = paste0(death_language$adjective, " death / all deaths"),
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank()
    ) +
    predict_curve_theme

  if (identical(as.integer(round(horizon_day)), 1000L)) {
    annotation_summary <- burden_all %>%
      filter(as.character(cohort) %in% c("2N", "4N"), is.finite(day)) %>%
      transmute(
        cohort = factor(as.character(cohort), levels = c("2N", "4N")),
        day = as.numeric(day),
        burden = as.numeric(pred_burden),
        live_cells = as.numeric(.first_non_null_local(pred_burden_live_cells, NA_real_)),
        o2_pct = as.numeric(clip(.first_non_null_local(pred_o2_pct, NA_real_), o2_plot_min, o2_plot_max))
      ) %>%
      group_by(cohort, day) %>%
      summarise(
        burden = mean(burden, na.rm = TRUE),
        live_cells = mean(live_cells, na.rm = TRUE),
        o2_pct = mean(o2_pct, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        burden = ifelse(is.nan(burden), NA_real_, burden),
        live_cells = ifelse(is.nan(live_cells), NA_real_, live_cells),
        o2_pct = ifelse(is.nan(o2_pct), NA_real_, o2_pct)
      )

    p_annotation_title <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0.95,
        hjust = 0,
        vjust = 1,
        label = paste0("Predicted (0-", as.integer(round(horizon_day)), " day)"),
        size = 3.8
      ) +
      annotate(
        "text",
        x = 0,
        y = 0.25,
        hjust = 0,
        vjust = 1,
        label = paste0("Column annotations are cohort-level means; fit_dir=", basename(dirname(out_dir)), " | report_dt=", report_dt),
        size = 2.4
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))

    p_annotation_legend <- make_predicted_annotation_legend(
      annotation_summary = annotation_summary,
      o2_plot_min = o2_plot_min,
      o2_plot_max = o2_plot_max
    )

    p_annotation_burden <- make_predict_annotation_track_plot(
      df = annotation_summary,
      value_col = "burden",
      y_label = "Burden",
      legend_title = "Burden\n(mm^3)",
      day_width = chr_density_day_width,
      horizon_day = horizon_day,
      x_breaks = predict_x_breaks,
      colors = c("#ffffbf", "#542788"),
      transform = "log10",
      show_legend = FALSE
    )
    p_annotation_live <- make_predict_annotation_track_plot(
      df = annotation_summary,
      value_col = "live_cells",
      y_label = "Viable cells",
      legend_title = "Viable cells",
      day_width = chr_density_day_width,
      horizon_day = horizon_day,
      x_breaks = predict_x_breaks,
      colors = c("#ffffbf", "#2c7bb6"),
      transform = "log10",
      show_legend = FALSE
    )
    p_annotation_o2 <- make_predict_annotation_track_plot(
      df = annotation_summary,
      value_col = "o2_pct",
      y_label = "Effective O2 (%)",
      legend_title = "Effective O2 (%)",
      day_width = chr_density_day_width,
      horizon_day = horizon_day,
      x_breaks = predict_x_breaks,
      colors = c("#f7f7f7", "#9ecae1", "#08519c"),
      transform = "identity",
      limits = c(o2_plot_min, o2_plot_max),
      breaks = unique(c(o2_plot_min, o2_plot_max)),
      labels = format(unique(c(o2_plot_min, o2_plot_max)), trim = TRUE),
      show_legend = FALSE
    )

    chr_density_fill_max_annot <- max(chr_density_df$density, na.rm = TRUE)
    if (!is.finite(chr_density_fill_max_annot) || chr_density_fill_max_annot <= 0) {
      chr_density_fill_max_annot <- 1
    }
    p_annotation_chr_density <- ggplot(
      chr_density_df,
      aes(x = day, y = chr_bin_mid, fill = density)
    ) +
      geom_tile(width = chr_density_day_width, height = chr_density_bin_width) +
      cohort_facet_grid() +
      predict_x_scale() +
      ploidy_fraction_fill_scale(chr_density_fill_max_annot, name = "Chromosome\nprobability") +
      scale_y_continuous(
        breaks = ploidy_n_breaks,
        expand = c(0, 0),
        sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Ploidy")
      ) +
      coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
      labs(
        x = NULL,
        y = "Chromosome count (N)"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )

    mean_chr_summary <- endpoint_plot_df %>%
      filter(as.character(cohort) %in% c("2N", "4N"), is.finite(day), is.finite(value_chr)) %>%
      mutate(cohort = factor(as.character(cohort), levels = c("2N", "4N"))) %>%
      group_by(cohort, day) %>%
      summarise(value_chr = mean(value_chr, na.rm = TRUE), .groups = "drop")
    p_annotation_mean_chr <- ggplot(
      mean_chr_summary,
      aes(x = day, y = value_chr, group = cohort, color = cohort)
    ) +
      geom_line(linewidth = 0.85, alpha = 0.95) +
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
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )

    save_aligned_plot_stack(
      list(
        p_annotation_title,
        p_annotation_legend,
        p_annotation_burden,
        p_annotation_live,
        p_annotation_o2,
        p_annotation_chr_density,
        p_annotation_mean_chr
      ),
      file.path(out_dir, paste0("predicted_", horizon_tag, ".pdf")),
      width = 12,
      height = 4.8,
      row_heights = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25)
    )
  }

  save_aligned_plot_stack(
    list(
      p_predict_burden,
      p_predict_live_cells,
      p_resource_death_fraction_predict,
      p_predict_endpoint,
      p_predict_chr_density,
      p_predict_o2
    ),
    file.path(out_dir, paste0("predict_curves_", horizon_tag, ".pdf")),
    width = 12,
    height = 10.8,
    row_heights = c(1, 1, 1, 1, 1.5, 1)
  )

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

  burden_decomp_predict_floor <- log10_plot_floor(
    c(
      burden_decomp_predict$burden_live,
      burden_decomp_predict$burden_dead_hypoxia,
      burden_decomp_predict$burden_dead_buffer,
      burden_decomp_predict$burden_total
    ),
    default = 1e-12
  )
  burden_decomp_predict <- burden_decomp_predict %>%
    mutate(burden_total_log_plot = floor_for_log10_plot(burden_total, burden_decomp_predict_floor))
  burden_decomp_predict_long <- make_burden_decomp_long(burden_decomp_predict, death_language)
  burden_decomp_predict_ribbon <- make_burden_decomp_ribbon(
    burden_decomp_predict,
    death_language,
    burden_decomp_predict_floor
  )

  p_burden_decomp_predict <- ggplot(
    burden_decomp_predict_ribbon,
    aes(x = day, ymin = ymin, ymax = ymax, fill = component, group = component)
  ) +
    geom_ribbon(alpha = 0.55) +
    geom_line(
      data = burden_decomp_predict,
      aes(x = day, y = burden_total_log_plot),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.65
    ) +
    facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = stats::setNames(
      c("#1f77b4", "#d62728", "#2ca02c"),
      c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
    )) +
    log10_burden_y_scale() +
    coord_cartesian(xlim = c(0, horizon_day)) +
    labs(
      title = paste0("Predicted Total Burden Viable/Dead Decomposition: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = "Cohort-level mean across scenarios (2N top, 4N bottom)",
      x = "Day",
      y = "log10 Tumor burden (mm^3)",
      fill = "Component"
    ) +
    theme_bw(base_size = 11)

  burden_decomp_plot_for_cohort <- function(cohort_use, show_fill_legend = TRUE) {
    row_df <- burden_decomp_predict %>%
      filter(cohort == cohort_use)
    row_ribbon <- burden_decomp_predict_ribbon %>%
      filter(cohort == cohort_use)
    ggplot(
      row_ribbon,
      aes(x = day, ymin = ymin, ymax = ymax, fill = component, group = component)
    ) +
      geom_ribbon(alpha = 0.55) +
      geom_line(
        data = row_df,
        aes(x = day, y = burden_total_log_plot),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.65
      ) +
      scale_fill_manual(
        values = stats::setNames(
          c("#1f77b4", "#d62728", "#2ca02c"),
          c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
        )
      ) +
      log10_burden_y_scale() +
      coord_cartesian(xlim = c(0, horizon_day)) +
      labs(
        title = if (identical(as.character(cohort_use), "2N")) {
          paste0("Predicted Total Burden Viable/Dead Decomposition: 0-", as.integer(round(horizon_day)), " days")
        } else {
          NULL
        },
        subtitle = paste0(as.character(cohort_use), " cohort mean across scenarios"),
        x = "Day",
        y = "log10 Tumor burden (mm^3)",
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

  invisible(list(
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
    c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
  )

  plots <- lapply(predict_results, function(res) {
    horizon_day <- as.numeric(res$horizon_day)
    decomp_floor <- log10_plot_floor(
      c(
        res$burden_decomp_predict$burden_live,
        res$burden_decomp_predict$burden_dead_hypoxia,
        res$burden_decomp_predict$burden_dead_buffer,
        res$burden_decomp_predict$burden_total
      ),
      default = 1e-12
    )
    burden_decomp_predict <- res$burden_decomp_predict %>%
      mutate(burden_total_log_plot = floor_for_log10_plot(burden_total, decomp_floor))
    burden_decomp_predict_ribbon <- make_burden_decomp_ribbon(
      burden_decomp_predict,
      death_language,
      decomp_floor
    )
    ggplot(
      burden_decomp_predict_ribbon,
      aes(x = day, ymin = ymin, ymax = ymax, fill = component, group = component)
    ) +
      geom_ribbon(alpha = 0.55) +
      geom_line(
        data = burden_decomp_predict,
        aes(x = day, y = burden_total_log_plot),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.65
      ) +
      facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = fill_values, drop = FALSE) +
      log10_burden_y_scale() +
      coord_cartesian(xlim = c(0, horizon_day)) +
      labs(
        title = paste0("0-", as.integer(round(horizon_day)), " days"),
        subtitle = "Cohort-level mean across scenarios (2N top, 4N bottom)",
        x = "Day",
        y = "log10 Tumor burden (mm^3)",
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
  keep <- basename(entries) %in% c("invivo", "invitro", "joint_parameters")
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
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
run_viz_for_fit_dir <- function(
  fit_dir,
  argv,
  dt_path,
  ploidy_path,
  report_dt,
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
  death_language <- resource_death_language()
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
  misseg_timecourse <- compute_missegregation_probability_timecourse(
    ploidy_all = ploidy_all,
    burden_all = burden_all,
    run_params = run_params
  )
  if (nrow(misseg_timecourse) > 0L) {
    write.table(
      misseg_timecourse,
      file = file.path(out_dir, "missegregation_probability_timecourse.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    plot_missegregation_probability_timecourse(
      ms_timecourse = misseg_timecourse,
      out_dir = out_dir,
      fit_dir = fit_dir,
      report_dt = report_dt
    )
  }
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
      mutate(o2_series = factor(o2_series, levels = c("o2_target_pct", "o2_eff_pct"), labels = c("O2 target", "Effective O2")))

    p_o2_lag <- ggplot(o2_lag_long, aes(x = day, y = o2_pct, color = o2_series, linetype = o2_series, group = interaction(sample_id, o2_series))) +
      geom_line(linewidth = 0.7, alpha = 0.85) +
      facet_wrap(~ harvest, ncol = 2) +
      scale_color_manual(values = c("O2 target" = "#ff7f0e", "Effective O2" = "#1f77b4")) +
      coord_cartesian(ylim = c(o2_plot_min, o2_plot_max)) +
      labs(
        title = "Resource-Stress Model: Effective Oxygen Relaxation Over Time",
        subtitle = "Oxygen supply-demand target vs lagged effective oxygen state",
        x = "Day",
        y = "Effective oxygen (%)",
        color = NULL,
        linetype = NULL
      ) +
      theme_bw(base_size = 11)

  }
  unlink(file.path(out_dir, c(
    "g_lag_timecourse.tsv",
    "g_target_vs_eff_timecourse.pdf",
    "g_lag_gap_timecourse.pdf",
    "predict_burden_vs_g.tsv",
    "predict_burden_vs_g.pdf",
    "burden_trend.pdf",
    "burden_trend_absolute.pdf",
    "o2_lag_gap_timecourse.pdf",
    "oxygen_vs_live_state_pms_colored_by_live_fraction.pdf",
    "oxygen_livefraction_liveweighted_pms_surface.pdf",
    "ploidy_heatmap_over_time.pdf",
    "ploidy_top_states_over_time.pdf",
    "oxygen_response_4panel.pdf",
    "overview_9panel.pdf"
  )), force = TRUE)

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
      title = "Resource-Stress Model: In Vivo Tumor Burden Trajectory",
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
        labels = c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
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
    scale_fill_manual(values = stats::setNames(
      c("#1f77b4", "#d62728", "#2ca02c"),
      c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
    )) +
    labs(
      title = "Resource-Stress Model: Total Tumor Burden Viable/Dead Decomposition",
      subtitle = paste0(
        "Total burden (black) = viable + ",
        death_language$figure_phrase,
        " + CIN-associated dead"
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

  p_burden_vs_o2 <- ggplot(o2_burden_df, aes(x = burden_mm3, y = o2_pct, color = cohort, group = sample_id)) +
    geom_path(linewidth = 0.75, alpha = 0.9) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_x") +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    coord_cartesian(ylim = c(o2_plot_min, o2_plot_max)) +
    labs(
      title = "Resource-Stress Model: Predicted Effective Oxygen vs Tumor Burden",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Tumor burden (mm^3)",
      y = "Effective oxygen (%)",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  p_ploidy_weighted_mean <- ggplot(ploidy_mean, aes(x = day, y = weighted_mean_ploidy)) +
    geom_line(color = "#d62728", linewidth = 0.9) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(min(ploidy_mean$weighted_mean_ploidy, na.rm = TRUE), max(ploidy_mean$weighted_mean_ploidy, na.rm = TRUE))) +
    labs(
      title = paste0("Resource-Stress Model: ", weighted_mean_series_label(cfg), " Over Time"),
      subtitle = "Weighted by predicted viable chromosome-number state fractions",
      x = "Day",
      y = weighted_mean_series_label(cfg)
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "burden_trend_absolute(real_scale).pdf"), p_burden_abs_real, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_live_dead_decomposition.pdf"), p_burden_decomp, width = 13, height = 9)
  if (exists("p_o2_lag", inherits = FALSE)) ggsave(file.path(out_dir, "o2_target_vs_eff_timecourse.pdf"), p_o2_lag, width = 13, height = 9)
  ggsave(file.path(out_dir, "predict_burden_vs_o2.pdf"), p_burden_vs_o2, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_weighted_mean_over_time.pdf"), p_ploidy_weighted_mean, width = 13, height = 9)
  plot_functional_response_curves(
    run_params = run_params,
    cfg = cfg,
    out_dir = out_dir
  )
  predict_horizons <- as_num_vec(argv$predict_horizons, c(100, 300, 1000))
  predict_horizons <- sort(unique(predict_horizons[is.finite(predict_horizons) & predict_horizons > 0]))
  predict_report_dt <- as_num(argv$predict_report_dt, report_dt)
  if (!is.finite(predict_report_dt) || predict_report_dt <= 0) predict_report_dt <- report_dt
  do_predict_plots <- as_bool(argv$predict_plots, TRUE)
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
        report_dt = predict_report_dt
      )
      if (is.list(p_hz)) {
        predict_results[[length(predict_results) + 1L]] <- p_hz
      }
    }
    plot_predict_burden_live_dead_decomposition_combined(
      predict_results = predict_results,
      out_dir = out_dir,
      death_language = death_language
    )
  }

  if (isTRUE(joint_fit)) {
    tryCatch(
      {
        plot_joint_parameter_ratio_figure(fit_dir = fit_dir)
      },
      error = function(e) {
        message("  Joint parameter ratio figure skipped: ", conditionMessage(e))
      }
    )
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
  model_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")
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

  data_dir <- if (!is.null(argv$data_dir)) {
    argv$data_dir
  } else {
    data_dir_candidates <- c(
      file.path(script_dir, "..", "..", "..", "data", "InVivoData_Gemcitabine"),
      file.path(script_dir, "..", "..", "..", "..", "data", "InVivoData_Gemcitabine")
    )
    data_dir_hit <- data_dir_candidates[dir.exists(data_dir_candidates)]
    normalizePath(if (length(data_dir_hit)) data_dir_hit[[1L]] else data_dir_candidates[[1L]], mustWork = FALSE)
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
