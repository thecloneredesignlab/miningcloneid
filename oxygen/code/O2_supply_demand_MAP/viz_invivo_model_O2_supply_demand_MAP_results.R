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
source(file.path(.o2sd_bootstrap_script_dir, "o2_supply_demand_map_shared.R"), local = environment())
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
  needed <- c(
    "lam_min", "lam_max", "k_o", "p_misseg", "k_o_mis",
    "gamma_loss", "p_wgd",
    "o2_S0", "kappa_O", "eta_o2",
    "beta_size", "alpha_o2", "gamma_growth",
    "mu_hp", "gamma_mu", "O2_crit", "n_O", "k_clear"
  )
  miss <- setdiff(needed, names(vals))
  if (length(miss) > 0) {
    stop("best_params.tsv missing parameters: ", paste(miss, collapse = ", "))
  }
  out <- as.list(vals[needed])
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

  sim_cpp <- cpp_o2simps_simulate_one(
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
    boundary = as.character(.first_non_null_local(run_params$boundary, "drop")),
    eps_tail = as.numeric(1e-8),
    gamma_loss = as.numeric(.first_non_null_local(run_params$gamma_loss, 0.1)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = as.numeric(.first_non_null_local(run_params$beta_size, cfg$prior_center_beta_size, default_beta_size_prior_center())),
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
    mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)),
    gamma_mu = as.numeric(.first_non_null_local(run_params$gamma_mu, cfg$gamma_mu_init, 1.0)),
    n_O = as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0)),
    # Keep config mode authoritative in viz re-simulation.
    ploidy_O2_death = assert_canonical_ploidy_o2_death_mode(
      .first_non_null_local(cfg$ploidy_O2_death, run_params$ploidy_O2_death, "diploid_NULL")
    ),
    k_clear = as.numeric(.first_non_null_local(run_params$k_clear, cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(burden_floor),
    return_full_trajectory = TRUE
  )

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
      pred_o2_target_pct, pred_o2_pct, pred_o2_lag_gap_pct, obs_burden
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
  ploidy_all %>%
    group_by(harvest, cohort, dose, day) %>%
    summarise(
      weighted_mean_N = sum(N * fraction, na.rm = TRUE) / pmax(sum(fraction, na.rm = TRUE), 1e-12),
      .groups = "drop"
    ) %>%
    mutate(weighted_mean_ploidy = weighted_mean_N / cfg$N_UNIT)
}

# -----------------------------------------------------------------------------
# Function: build_terminal_ploidy_compare_df
# Purpose: Build observed-vs-predicted ploidy distributions at the time points
#   used by the ploidy objective.
# -----------------------------------------------------------------------------
build_terminal_ploidy_compare_df <- function(scenarios, ploidy_all, cfg) {
  meta_rows <- vector("list", length(scenarios))
  obs_rows <- vector("list", length(scenarios))

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    obs_z <- as.numeric(sc$ploidy_obs_z)
    obs_z <- obs_z[is.finite(obs_z)]
    if (length(obs_z) == 0L) next

    target_day <- as.numeric(sc$sim_end_day)
    if (!is.finite(target_day)) next

    meta_rows[[i]] <- data.frame(
      harvest = as.character(sc$harvest),
      cohort = as.character(sc$cohort),
      dose = as.numeric(sc$dose),
      target_day = target_day,
      stringsAsFactors = FALSE
    )

    obs_N <- map_ploidy_to_N_by_chrlen(
      ploidy_values = obs_z / as.numeric(cfg$N_UNIT),
      N_grid = cfg$N_MIN:cfg$N_MAX,
      chr_lengths_bp = cfg$chr_lengths_bp
    )
    obs_N <- as.integer(clip(obs_N, cfg$N_MIN, cfg$N_MAX))
    obs_N <- obs_N[is.finite(obs_N)]
    if (length(obs_N) == 0L) next

    obs_rows[[i]] <- data.frame(
      harvest = as.character(sc$harvest),
      cohort = as.character(sc$cohort),
      dose = as.numeric(sc$dose),
      target_day = target_day,
      source = "Observed",
      ploidy = as.numeric(obs_N) / as.numeric(cfg$N_UNIT),
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
      ploidy = as.numeric(N) / as.numeric(cfg$N_UNIT),
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
      q1 = weighted_quantile_local(ploidy, weight, 0.25),
      median = weighted_quantile_local(ploidy, weight, 0.5),
      q3 = weighted_quantile_local(ploidy, weight, 0.75),
      ymin_raw = min(ploidy[weight > 0], na.rm = TRUE),
      ymax_raw = max(ploidy[weight > 0], na.rm = TRUE),
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

  p <- ggplot(compare_df, aes(x = source, y = ploidy, weight = weight, fill = fill_group)) +
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
    facet_wrap(~ cohort, nrow = 1) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    labs(
      title = "Observed vs Predicted Ploidy Distributions Used in Ploidy Objective",
      subtitle = paste0(
        "Observed ploidy is mapped to the chromosome-count grid used by the objective | fit_dir=",
        basename(fit_dir)
      ),
      x = NULL,
      y = "Ploidy (2N scale)",
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )

  ggsave(file.path(out_dir, "terminal_ploidy_observed_vs_predicted_violin.pdf"), p, width = 10, height = 6.5)
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
plot_functional_response_curves <- function(run_params, cfg, out_dir) {
  o2_plot_min <- 0
  o2_plot_max <- 5
  o2_grid <- seq(o2_plot_min, o2_plot_max, by = 0.02)
  ploidy_plot_min <- 0
  ploidy_plot_max <- 10
  ploidy_grid_dense <- seq(ploidy_plot_min, ploidy_plot_max, by = 0.05)
  o2_levels_ploidy <- seq(o2_plot_min, o2_plot_max, by = 0.5)
  ref_df <- data.frame(
    cohort = c("2N", "4N"),
    N_ref = as.numeric(c(2 * cfg$N_UNIT, 4 * cfg$N_UNIT)),
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
  n_O <- as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0))
  if (!is.finite(n_O) || n_O < 0) stop("run_params$n_O must be finite and >= 0.")
  mu_hp_use <- pmax(as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)), 0)
  N_dip <- 44.0

  compute_rate_components <- function(O2, N) {
    O2_vec <- as.numeric(O2)
    N_vec <- as.numeric(N)
    n_out <- max(length(O2_vec), length(N_vec))
    if (!(length(O2_vec) %in% c(1L, n_out) && length(N_vec) %in% c(1L, n_out))) {
      stop("O2 and N must have compatible lengths in compute_rate_components().")
    }
    O2_vec <- rep_len(O2_vec, n_out)
    N_vec <- rep_len(N_vec, n_out)
    lam_base <- as.numeric(growth_lambda(
      O2 = O2_vec,
      N = N_vec,
      lam_min = run_params$lam_min,
      lam_max = run_params$lam_max,
      k_o = run_params$k_o
    ))
    o2_c <- pmax(O2_crit_use, 1e-12)
    h_o2 <- (o2_c^n_O) / ((o2_c^n_O) + (pmax(O2_vec, 0)^n_O))
    h_o2 <- pmax(0, pmin(1, h_o2))
    proliferation_rate <- if (isTRUE(o2_growth_use)) {
      denom <- 1 + alpha_o2_use * h_o2 * ((pmax(N_vec, 0) / N_dip)^gamma_growth_use)
      lam_base / pmax(denom, 1e-12)
    } else {
      lam_base
    }
    death_rate <- if (identical(ploidy_O2_death_use, "uniform")) {
      mu_hp_use * h_o2
    } else if (identical(ploidy_O2_death_use, "diploid_NULL")) {
      mu_hp_use * h_o2 * (1 + pmax(N_vec / N_dip - 1, 0)^gamma_mu_use)
    } else {
      mu_hp_use * h_o2 * pmax(N_vec / N_dip, 0)^gamma_mu_use
    }
    data.frame(
      O2 = O2_vec,
      N = N_vec,
      proliferation_rate = pmax(as.numeric(proliferation_rate), 0),
      death_rate = pmax(as.numeric(death_rate), 0),
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
    rate_df <- compute_rate_components(O2 = o2_grid, N = rep(N_ref, length(o2_grid)))
    data.frame(
      oxygen_pct = o2_grid,
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
    o2_curve,
    file = file.path(out_dir, "functional_curve_oxygen.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

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

  p_msr_death <- ggplot(o2_curve, aes(x = death_rate, y = ms_rate, color = cohort)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = "MS Rate vs Death Rate",
      x = "Death rate",
      y = "MS rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)

  p_prolif <- ggplot(o2_curve, aes(x = oxygen_pct, y = proliferation_rate, color = cohort)) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = "Oxygen vs Proliferation Rate",
      x = "Oxygen (%)",
      y = "Proliferation rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  p_death <- ggplot(o2_curve, aes(x = oxygen_pct, y = death_rate, color = cohort)) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#9467bd")) +
    labs(
      title = "Oxygen vs Death Rate",
      x = "Oxygen (%)",
      y = "Death rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  p_net <- ggplot(o2_curve, aes(x = oxygen_pct, y = net_growth_rate)) +
    geom_line(linewidth = 1, color = "#2ca02c") +
    facet_wrap(~cohort, ncol = 2) +
    coord_cartesian(xlim = c(o2_plot_min, o2_plot_max)) +
    labs(
      title = "Oxygen vs Net Growth Rate",
      subtitle = "Net rate = proliferation - hypoxia-linked high-ploidy death",
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
    viability_after_ms = pmax(viability, 0),
    row.names = NULL
  )
  write.table(
    viability_curve,
    file = file.path(out_dir, "functional_curve_ploidy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  ploidy_o2_curve <- dplyr::bind_rows(lapply(o2_levels_ploidy, function(o2_level) {
    N_grid <- as.numeric(ploidy_grid_dense * as.numeric(cfg$N_UNIT))
    rate_df <- compute_rate_components(O2 = rep(o2_level, length(N_grid)), N = N_grid)
    data.frame(
      oxygen_pct = rep(as.numeric(o2_level), length(ploidy_grid_dense)),
      ploidy = as.numeric(ploidy_grid_dense),
      N = N_grid,
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

  p_viability <- ggplot(viability_curve, aes(x = ploidy, y = viability_after_ms)) +
    geom_line(color = "#2ca02c", linewidth = 1) +
    labs(
      title = "Ploidy vs Viability After MS",
      subtitle = "Nullisomy-risk loss survival for a one-copy loss event",
      x = "Ploidy (N / N_UNIT)",
      y = "Viability after MS"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_prolif_o2 <- ggplot(
    ploidy_o2_curve,
    aes(x = ploidy, y = proliferation_rate, color = oxygen_pct)
  ) +
    geom_point(shape = 15, size = 1.8, alpha = 0.95) +
    coord_cartesian(xlim = c(ploidy_plot_min, ploidy_plot_max)) +
    scale_color_gradient(
      low = "#2C7BB6",
      high = "#F28E2B",
      limits = c(o2_plot_min, o2_plot_max),
      name = "O2 level"
    ) +
    labs(
      title = "Ploidy vs Proliferation Rate Colored by O2",
      subtitle = "Dense square markers over ploidy range 0-10, colored by O2 level 0-5",
      x = "Ploidy",
      y = "Proliferation rate"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_death_o2 <- ggplot(
    ploidy_o2_curve,
    aes(x = ploidy, y = death_rate, color = oxygen_pct)
  ) +
    geom_point(shape = 15, size = 1.8, alpha = 0.95) +
    coord_cartesian(xlim = c(ploidy_plot_min, ploidy_plot_max)) +
    scale_color_gradient(
      low = "#2C7BB6",
      high = "#F28E2B",
      limits = c(o2_plot_min, o2_plot_max),
      name = "O2 level"
    ) +
    labs(
      title = "Ploidy vs Death Rate Colored by O2",
      subtitle = "Dense square markers over ploidy range 0-10, colored by O2 level 0-5",
      x = "Ploidy",
      y = "Death rate"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "oxygen_vs_missegregation_rate.pdf"), p_msr_o2, width = 10, height = 7)
  ggsave(file.path(out_dir, "ms_rate_vs_death_rate.pdf"), p_msr_death, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_proliferation_rate.pdf"), p_prolif, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_death_rate.pdf"), p_death, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_net_growth_rate.pdf"), p_net, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_viability_after_ms.pdf"), p_viability, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_proliferation_rate_by_o2.pdf"), p_ploidy_prolif_o2, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_death_rate_by_o2.pdf"), p_ploidy_death_o2, width = 10, height = 7)

  invisible(list(
    p_msr_o2 = p_msr_o2,
    p_msr_death = p_msr_death,
    p_prolif = p_prolif,
    p_death = p_death,
    p_net = p_net,
    p_viability = p_viability,
    p_ploidy_prolif_o2 = p_ploidy_prolif_o2,
    p_ploidy_death_o2 = p_ploidy_death_o2
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
plot_predict_horizon <- function(run_params, scenarios, cfg, out_dir, horizon_day, report_dt = 1.0, top_n = 6L) {
  sim_list <- lapply(scenarios, function(sc) {
    simulate_one_full_horizon(run_params, sc, cfg, horizon_day = horizon_day, report_dt = report_dt)
  })
  burden_all <- bind_rows(lapply(sim_list, `[[`, "burden"))
  ploidy_all <- bind_rows(lapply(sim_list, `[[`, "ploidy"))
  if (nrow(burden_all) == 0 || nrow(ploidy_all) == 0) return(invisible(NULL))

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
    paste0("forecast_ploidy_weighted_mean_", horizon_tag, ".pdf")
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
      metric = "Burden (absolute)",
      value = as.numeric(pred_burden)
    )

  ploidy_plot_df <- ploidy_mean %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      metric = "Weighted mean ploidy",
      value = as.numeric(weighted_mean_ploidy)
    )

  predict_plot_df <- bind_rows(burden_plot_df, ploidy_plot_df) %>%
    mutate(
      sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__"),
      metric = factor(metric, levels = c("Burden (absolute)", "Weighted mean ploidy"))
    )

  p_predict <- ggplot(
    predict_plot_df,
    aes(x = day, y = value, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    coord_cartesian(xlim = c(0, horizon_day)) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = paste0("Predict Curves: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = paste0("Single summary plot (all scenarios overlaid) | fit_dir=", basename(dirname(out_dir)), " | report_dt=", report_dt),
      x = "Day",
      y = NULL,
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_dir, paste0("predict_curves_", horizon_tag, ".pdf")), p_predict, width = 12, height = 11)

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
        labels = c("Live", "Dead (Hypoxia)", "Dead (Buffer loss)")
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
    scale_fill_manual(values = c("Live" = "#1f77b4", "Dead (Hypoxia)" = "#d62728", "Dead (Buffer loss)" = "#2ca02c")) +
    coord_cartesian(xlim = c(0, horizon_day)) +
    labs(
      title = paste0("Predict Burden Live/Dead Decomposition: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = "Cohort-level mean across scenarios (2N top, 4N bottom)",
      x = "Day",
      y = "Tumor burden (mm^3)",
      fill = "Component"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, paste0("predict_burden_live_dead_decomposition_", horizon_tag, ".pdf")), p_burden_decomp_predict, width = 12, height = 11)

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
    p_predict = p_predict,
    p_o2_time = p_o2_time,
    p_burden_decomp_predict = p_burden_decomp_predict,
    p_death_ratio_predict = p_death_ratio_predict
  ))
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
  top_n
) {
  out_dir <- file.path(fit_dir, "viz")
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
        labels = c("Live", "Dead (Hypoxia)", "Dead (Buffer loss)")
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
    scale_fill_manual(values = c("Live" = "#1f77b4", "Dead (Hypoxia)" = "#d62728", "Dead (Buffer loss)" = "#2ca02c")) +
    labs(
      title = "O2 Supply-Demand MAP Model: Live/Dead Burden Decomposition",
      subtitle = "Total burden (black) = live + dead from hypoxia + dead from buffer-derived nonviable offspring",
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
      title = "O2 Supply-Demand MAP Model: Predicted Oxygen vs Burden",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Tumor burden (mm^3)",
      y = "Oxygen (%)",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)

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
      title = "O2 Supply-Demand MAP Model: Predicted Ploidy Distribution Over Time",
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
      title = paste0("O2 Supply-Demand MAP Model: Ploidy Over Time (Top ", top_n, " N States)"),
      x = "Day",
      y = "Fraction",
      color = "N"
    ) +
    theme_bw(base_size = 11)

  p_ploidy_weighted_mean <- ggplot(ploidy_mean, aes(x = day, y = weighted_mean_ploidy)) +
    geom_line(color = "#d62728", linewidth = 0.9) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(min(ploidy_mean$weighted_mean_ploidy, na.rm = TRUE), 5)) +
    labs(
      title = "O2 Supply-Demand MAP Model: Weighted Mean Ploidy Over Time",
      subtitle = "Weighted by predicted ploidy fractions",
      x = "Day",
      y = "Weighted Mean Ploidy (P = N / N_UNIT)"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "burden_trend.pdf"), p_burden, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_trend_absolute.pdf"), p_burden_abs, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_trend_absolute(real_scale).pdf"), p_burden_abs_real, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_live_dead_decomposition.pdf"), p_burden_decomp, width = 13, height = 9)
  if (exists("p_o2_lag", inherits = FALSE)) ggsave(file.path(out_dir, "o2_target_vs_eff_timecourse.pdf"), p_o2_lag, width = 13, height = 9)
  if (exists("p_o2_lag_gap", inherits = FALSE)) ggsave(file.path(out_dir, "o2_lag_gap_timecourse.pdf"), p_o2_lag_gap, width = 13, height = 9)
  ggsave(file.path(out_dir, "predict_burden_vs_o2.pdf"), p_burden_vs_o2, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_heatmap_over_time.pdf"), p_ploidy_heatmap, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_top_states_over_time.pdf"), p_ploidy_lines, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_weighted_mean_over_time.pdf"), p_ploidy_weighted_mean, width = 13, height = 9)
  functional_plots <- plot_functional_response_curves(run_params = run_params, cfg = cfg, out_dir = out_dir)
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
      all(c("p_msr_death", "p_prolif", "p_death") %in% names(functional_plots))) {
    p_o2_panel <- p_o2_lag +
      labs(
        title = "Oxygen Evolution Over Time",
        subtitle = NULL
      )
    p_o2_panel <- legend_inside(p_o2_panel, x = 0.98, y = 0.98)
    p_msr_death_panel <- legend_inside(functional_plots$p_msr_death, x = 0.98, y = 0.98)
    p_prolif_panel <- legend_inside(functional_plots$p_prolif, x = 0.98, y = 0.98)
    p_death_panel <- legend_inside(functional_plots$p_death, x = 0.98, y = 0.98)
    grDevices::pdf(file.path(out_dir, "oxygen_response_4panel.pdf"), width = 18, height = 12, onefile = TRUE)
    grid::grid.newpage()
    lay <- grid::grid.layout(nrow = 2, ncol = 2)
    grid::pushViewport(grid::viewport(layout = lay))
    print(p_o2_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_msr_death_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
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
  }

  if (is.list(functional_plots) &&
      all(c("p_msr_death", "p_prolif", "p_death") %in% names(functional_plots)) &&
      inherits(p_burden_abs_real, "ggplot") &&
      inherits(p_burden_decomp, "ggplot") &&
      inherits(p_ploidy_weighted_mean, "ggplot") &&
      inherits(p_o2_1000_for_overview, "ggplot") &&
      inherits(p_predict_for_overview, "ggplot") &&
      inherits(p_burden_decomp_predict_for_overview, "ggplot")) {
    p_row1_left <- p_ploidy_weighted_mean +
      labs(title = "Ploidy Weighted Mean Over Time", subtitle = NULL) +
      theme(legend.position = "none")
    p_row1_right <- p_burden_abs_real +
      labs(title = "Burden Trend Absolute (Real Scale)", subtitle = NULL) +
      theme(legend.position = "none")
    p_row2_col1 <- legend_inside(functional_plots$p_msr_death, x = 0.98, y = 0.98)
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
  source(file.path(script_dir, "model_O2_supply_demand_MAP.R"))

  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  results_root <- normalizePath(file.path(script_dir, "..", "..", "results"), mustWork = FALSE)
  fit_root <- if (!is.null(argv$fit_dir)) {
    normalizePath(argv$fit_dir, mustWork = TRUE)
  } else {
    normalizePath(find_latest_fit_dir(results_root), mustWork = TRUE)
  }

  if (!is.null(argv$out_dir)) {
    message("Ignoring --out_dir. Outputs are always written to each subdirectory's /viz.")
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
  ploidy_path <- file.path(data_dir, "all_ploidy.tsv")

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
          top_n = top_n
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
