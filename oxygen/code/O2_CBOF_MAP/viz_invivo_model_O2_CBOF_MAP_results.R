#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Matrix))

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

# -----------------------------------------------------------------------------
# Function: get_script_dir_self
# Purpose: Resolve script directory path for robust relative file loading.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
get_script_dir_self <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grepl("^--file=", args)]
  if (length(farg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", farg[[1]])))
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
  cfg$N_UNIT <- as.integer(cfg$N_UNIT %||% 22L)
  cfg$N_MIN <- as.integer(cfg$N_MIN %||% 22L)
  cfg$N_MAX <- as.integer(cfg$N_MAX %||% 154L)
  cfg$DT <- as.numeric(cfg$DT %||% 0.5)
  cfg$O2_fixed <- as.numeric(cfg$O2_fixed %||% 5.0)
  cfg$o2_cap_pct <- as.numeric(cfg$o2_cap_pct %||% 5.0)
  if (!is.finite(cfg$O2_fixed) || cfg$O2_fixed < 0 || cfg$O2_fixed > 100) {
    stop("fit_config O2_fixed must be in percent scale [0, 100].")
  }
  if (!is.finite(cfg$o2_cap_pct) || cfg$o2_cap_pct < 0 || cfg$o2_cap_pct > 100) {
    stop("fit_config o2_cap_pct must be in percent scale [0, 100].")
  }
  cfg$K <- as.numeric(cfg$K %||% 1e12)
  cfg$crowding <- as.character(cfg$crowding %||% "logistic")
  cfg$init_total_size <- as.numeric(cfg$init_total_size %||% 1e6)
  cfg$dose_ref <- as.numeric(cfg$dose_ref %||% 30)
  cfg$tx_mult_min <- as.numeric(cfg$tx_mult_min %||% 0.05)
  cfg$min_pop <- as.numeric(cfg$min_pop %||% 1e-12)
  cfg$dose_zero_only <- isTRUE(cfg$dose_zero_only %||% TRUE)
  cfg$fit_treatment <- isTRUE(cfg$fit_treatment %||% FALSE)
  cfg$max_scenarios <- as.numeric(cfg$max_scenarios %||% Inf)
  cfg$growth_penalty_ploidy <- isTRUE(cfg$growth_penalty_ploidy %||% FALSE)
  cfg$growth_penalty_hypoxia <- isTRUE(cfg$growth_penalty_hypoxia %||% FALSE)
  if (!is.null(cfg$o2_min)) {
    cfg$o2_min <- as.numeric(cfg$o2_min)
    if (!is.finite(cfg$o2_min) || cfg$o2_min < 0 || cfg$o2_min > 100) {
      stop("fit_config o2_min must be in percent scale [0, 100].")
    }
    if (cfg$o2_min > cfg$o2_cap_pct) cfg$o2_min <- cfg$o2_cap_pct
  }

  if (is.null(cfg$truncate_at_treatment)) {
    cfg$truncate_at_treatment <- isTRUE(cfg$pretreat_only %||% FALSE)
  }
  cfg$truncate_at_treatment <- isTRUE(cfg$truncate_at_treatment)

  if (is.null(cfg$ploidy_at_harvest)) {
    cfg$ploidy_at_harvest <- TRUE
  }
  cfg$ploidy_at_harvest <- isTRUE(cfg$ploidy_at_harvest)
  cfg
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
    "beta_buffer", "n_exp", "smax", "p_wgd",
    "K_down", "A_ang", "m_on", "delta_m", "s_on", "s_off",
    "beta_size", "alpha_o2", "o2_ref_pct", "gamma_growth", "mu_hp", "k_clear"
  )
  miss <- setdiff(needed, names(vals))
  if (length(miss) > 0) {
    stop("best_params.tsv missing parameters: ", paste(miss, collapse = ", "))
  }
  out <- as.list(vals[needed])
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
  out$tau_O2 <- if ("tau_O2" %in% names(vals) && is.finite(vals[["tau_O2"]]) && vals[["tau_O2"]] > 0) vals[["tau_O2"]] else as.numeric(.first_non_null_local(cfg$tau_O2, cfg$tau_O2_init, 2.0))
  out$growth_penalty_ploidy <- isTRUE(.first_non_null_local(cfg$growth_penalty_ploidy, FALSE))
  out$growth_penalty_hypoxia <- isTRUE(.first_non_null_local(cfg$growth_penalty_hypoxia, FALSE))
  out
}

# -----------------------------------------------------------------------------
# Function: compute_o2_target_from_burden
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - Ntot: Total predicted cell count (or burden proxy) at current time.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
compute_o2_target_from_burden <- function(Ntot, run_params, cfg) {
  o2_feedback <- isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE))
  o2_cap <- clip(as.numeric(.first_non_null_local(cfg$o2_cap_pct, 5.0)), 0, 100)
  o2_base <- min(clip(as.numeric(.first_non_null_local(cfg$O2_fixed, 5.0)), 0, 100), o2_cap)
  o2_min <- clip(as.numeric(.first_non_null_local(cfg$o2_min, 0.0)), 0, 100)
  h_O2 <- as.numeric(.first_non_null_local(cfg$h_O2, 1.0))
  if (!is.finite(h_O2) || h_O2 <= 0) h_O2 <- 1.0
  o2_logN_eps <- as.numeric(.first_non_null_local(cfg$o2_logN_eps, 1.0))
  if (!is.finite(o2_logN_eps) || o2_logN_eps <= 0) o2_logN_eps <- 1.0

  O2_target <- o2_base
  if (isTRUE(o2_feedback)) {
    O2_target <- as.numeric(.o2_window_supply_from_burden(
      Ntot = Ntot,
      run_params = run_params,
      O2_base = o2_base,
      o2_cap = o2_cap,
      o2_min = o2_min,
      h_down = h_O2,
      o2_logN_eps = o2_logN_eps
    ))
    O2_target <- clip(O2_target, 0, o2_cap)
  }
  as.numeric(O2_target)
}

# -----------------------------------------------------------------------------
# Function: build_g_backend_for_viz
# Purpose: Build a G/get_G backend for visualization when fit-side model_core does
#          not expose matrix generators (e.g., C++ objective path).
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - grid_pre: Integer ploidy grid (pre-WGD layer).
#   - grid_post: Integer ploidy grid (post-WGD layer).
# Returns:
#   list(G=<Matrix>, get_G=<function>)
# -----------------------------------------------------------------------------
build_g_backend_for_viz <- function(run_params, cfg, grid_pre, grid_post) {
  if (!exists(".build_G_with_WGD", mode = "function", inherits = TRUE)) {
    stop("Visualization backend error: .build_G_with_WGD is unavailable.")
  }
  if (!exists("growth_lambda", mode = "function", inherits = TRUE)) {
    stop("Visualization backend error: growth_lambda is unavailable.")
  }

  beta_buffer <- as.numeric(.first_non_null_local(run_params$beta_buffer, run_params$beta, 0.0))
  n_exp <- as.numeric(.first_non_null_local(run_params$n_exp, 1.0))
  smax <- as.numeric(.first_non_null_local(run_params$smax, 1.0))
  p_wgd <- as.numeric(.first_non_null_local(run_params$p_wgd, 0.0))
  if (!is.finite(p_wgd)) p_wgd <- 0.0
  boundary_mode <- "drop"
  N_min <- as.integer(cfg$N_MIN)
  N_max <- as.integer(cfg$N_MAX)
  N_unit <- as.integer(cfg$N_UNIT)
  lam_min <- as.numeric(run_params$lam_min)
  lam_max <- as.numeric(run_params$lam_max)
  k_o <- as.numeric(run_params$k_o)

  G_cache <- new.env(parent = emptyenv())
  build_G_for_O2 <- function(O2_level) {
    O2_level <- clip(as.numeric(O2_level), 0, 100)
    key <- sprintf("%.3f", O2_level)
    if (!exists(key, envir = G_cache, inherits = FALSE)) {
      lambda0 <- as.numeric(growth_lambda(O2 = O2_level, N = grid_pre, lam_min = lam_min, lam_max = lam_max, k_o = k_o))
      lambda1 <- as.numeric(growth_lambda(O2 = O2_level, N = grid_post, lam_min = lam_min, lam_max = lam_max, k_o = k_o))
      if (exists(".pmisseg_of_O2", mode = "function", inherits = TRUE)) {
        p_mis <- as.numeric(.pmisseg_of_O2(O2_level, run_params))
      } else {
        k_o_mis_use <- max(as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)), 1e-12)
        p0 <- as.numeric(.first_non_null_local(run_params$p_misseg, 0.0))
        p_mis <- p0 * (1 - (O2_level / (O2_level + k_o_mis_use)))
      }
      p_mis <- clip(as.numeric(p_mis), 0, 1)
      G <- .build_G_with_WGD(
        N0min = N_min, N0max = N_max, lambda0_vec = lambda0,
        p0_vec = p_mis, wgd_prob_vec = p_wgd,
        N1min = N_min, N1max = N_max, lambda1_vec = lambda1,
        p1_vec = p_mis,
        boundary = boundary_mode, N_unit = N_unit,
        beta_buffer = beta_buffer, n_exp = n_exp, smax = smax
      )
      mu_hp_use <- pmax(as.numeric(.first_non_null_local(run_params$mu_hp, 0.0)), 0)
      if (mu_hp_use > 0) {
        o2_ref_use <- clip(as.numeric(.first_non_null_local(run_params$o2_ref_pct, 0.0)), 0, 100)
        hypoxia_deficit <- pmax(0, o2_ref_use - O2_level)
        d_pre <- pmax(0, as.numeric(grid_pre) / as.numeric(N_unit) - 2)
        d_post <- pmax(0, as.numeric(grid_post) / as.numeric(N_unit) - 2)
        mu_diag <- c(mu_hp_use * d_pre * hypoxia_deficit, mu_hp_use * d_post * hypoxia_deficit)
        if (any(is.finite(mu_diag) & mu_diag > 0)) {
          G <- G - Diagonal(n = nrow(G), x = as.numeric(mu_diag))
        }
      }
      assign(key, G, envir = G_cache)
    }
    get(key, envir = G_cache, inherits = FALSE)
  }

  list(
    G = build_G_for_O2(as.numeric(.first_non_null_local(cfg$O2_fixed, 5.0))),
    get_G = build_G_for_O2
  )
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
  model_core <- build_model_core(run_params, cfg)
  grid_pre <- model_core$grid_pre
  grid_post <- model_core$grid_post
  R0 <- model_core$R0
  R1 <- model_core$R1
  init_state <- if (scenario$cohort == "2N") model_core$init_state_2N else model_core$init_state_4N
  has_explicit_G_backend <- !is.null(model_core$G) && is.function(model_core$get_G)
  if (has_explicit_G_backend) {
    G <- model_core$G
    get_G <- model_core$get_G
  } else {
    g_backend <- build_g_backend_for_viz(run_params = run_params, cfg = cfg, grid_pre = grid_pre, grid_post = grid_post)
    G <- g_backend$G
    get_G <- g_backend$get_G
  }

  sim_end_step <- as.integer(round(scenario$sim_end_day / cfg$DT))
  obs_steps <- as.integer(round(scenario$obs_days / cfg$DT))
  report_steps <- as.integer(round(seq(0, scenario$sim_end_day, by = report_dt) / cfg$DT))
  keep_steps <- sort(unique(c(0L, sim_end_step, obs_steps, report_steps)))
  keep_steps <- keep_steps[keep_steps >= 0L & keep_steps <= sim_end_step]

  obs_map <- setNames(as.numeric(scenario$obs_burden), as.character(obs_steps))

  v <- as.numeric(init_state)
  dose_scaled <- scenario$dose / cfg$dose_ref
  if (!is.finite(dose_scaled) || dose_scaled < 0) dose_scaled <- 0
  vol_by_N <- cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg)
  o2_feedback <- isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE))
  o2_cap <- clip(as.numeric(.first_non_null_local(cfg$o2_cap_pct, 5.0)), 0, 100)
  o2_base <- min(clip(as.numeric(.first_non_null_local(cfg$O2_fixed, 5.0)), 0, 100), o2_cap)
  tau_O2_use <- as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  alpha_tau <- 1 - exp(-cfg$DT / tau_O2_use)
  O2_state <- if (isTRUE(o2_feedback)) {
    compute_o2_target_from_burden(Ntot = sum(v), run_params = run_params, cfg = cfg)
  } else {
    o2_base
  }

  h_down_use <- as.numeric(.first_non_null_local(run_params$h_down, cfg$h_down_init, 1.0))
  if (!is.finite(h_down_use) || h_down_use <= 0) h_down_use <- 1.0
  K_down_use <- as.numeric(.first_non_null_local(run_params$K_down, cfg$K, 1e12))
  if (!is.finite(K_down_use) || K_down_use <= 0) K_down_use <- 1e12
  m_on_use <- as.numeric(.first_non_null_local(run_params$m_on, cfg$o2_m_on_default, 9.0))
  if (!is.finite(m_on_use)) m_on_use <- 9.0
  delta_m_use <- as.numeric(.first_non_null_local(run_params$delta_m, cfg$o2_delta_m_default, 1.0))
  if (!is.finite(delta_m_use) || delta_m_use <= 0) delta_m_use <- 1.0
  m_off_use <- as.numeric(.first_non_null_local(run_params$m_off, m_on_use + delta_m_use))
  if (!is.finite(m_off_use) || m_off_use <= m_on_use) m_off_use <- m_on_use + delta_m_use
  s_on_use <- as.numeric(.first_non_null_local(run_params$s_on, cfg$o2_s_on_default, 0.3))
  if (!is.finite(s_on_use) || s_on_use <= 0) s_on_use <- 0.3
  s_off_use <- as.numeric(.first_non_null_local(run_params$s_off, cfg$o2_s_off_default, 0.3))
  if (!is.finite(s_off_use) || s_off_use <= 0) s_off_use <- 0.3
  p_wgd_use <- as.numeric(.first_non_null_local(run_params$p_wgd, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
  boundary_mode <- as.character(.first_non_null_local(run_params$boundary, "drop"))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 0)
  sim_cpp <- cpp_o2invivo_simulate_one(
    init_state = as.numeric(init_state),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = as.integer(keep_steps),
    sim_end_step = as.integer(sim_end_step),
    DT = as.numeric(cfg$DT),
    dose = as.numeric(scenario$dose),
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = as.numeric(scenario$treat_day),
    fit_treatment = isTRUE(cfg$fit_treatment),
    alpha = as.numeric(.first_non_null_local(run_params$alpha, 0.0)),
    gamma = as.numeric(.first_non_null_local(run_params$gamma, 1.0)),
    tx_mult_min = as.numeric(cfg$tx_mult_min),
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_base = as.numeric(o2_base),
    o2_cap = as.numeric(o2_cap),
    o2_feedback = isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE)),
    o2_min = as.numeric(.first_non_null_local(cfg$o2_min, 0.0)),
    h_O2 = as.numeric(h_down_use),
    K_down = as.numeric(K_down_use),
    A_ang = as.numeric(.first_non_null_local(run_params$A_ang, cfg$o2_A_ang_default, 0.0)),
    m_on = as.numeric(m_on_use),
    m_off = as.numeric(m_off_use),
    s_on = as.numeric(s_on_use),
    s_off = as.numeric(s_off_use),
    tau_O2 = as.numeric(tau_O2_use),
    o2_logN_eps = as.numeric(.first_non_null_local(cfg$o2_logN_eps, 1.0)),
    o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01)),
    o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct, 0.005)),
    o2_cache_profile = isTRUE(.first_non_null_local(cfg$o2_cache_profile, FALSE)),
    lam_min = as.numeric(run_params$lam_min),
    lam_max = as.numeric(run_params$lam_max),
    k_o = as.numeric(run_params$k_o),
    has_p_misseg = !is.null(run_params$p_misseg),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = 0.0,
    p_wgd = as.numeric(p_wgd_use),
    boundary = boundary_mode,
    eps_tail = as.numeric(1e-8),
    beta_buffer = as.numeric(.first_non_null_local(run_params$beta_buffer, 0.0)),
    n_exp = as.numeric(.first_non_null_local(run_params$n_exp, 1.0)),
    smax = as.numeric(.first_non_null_local(run_params$smax, 1.0)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = as.numeric(.first_non_null_local(run_params$beta_size, cfg$prior_center_beta_size, default_beta_size_prior_center())),
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
    o2_ref_pct = as.numeric(.first_non_null_local(run_params$o2_ref_pct, cfg$o2_ref_pct_init, 2.5)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
    growth_penalty_ploidy = isTRUE(.first_non_null_local(cfg$growth_penalty_ploidy, FALSE)),
    growth_penalty_hypoxia = isTRUE(.first_non_null_local(cfg$growth_penalty_hypoxia, FALSE)),
    mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)),
    k_clear = as.numeric(.first_non_null_local(run_params$k_clear, cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(burden_floor)
  )
  step_key <- as.character(keep_steps)
  pred_cells_total_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Ntot_total_obs, sim_cpp$Ntot_obs)), step_key)
  pred_cells_live_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Ntot_live_obs, sim_cpp$Ntot_obs)), step_key)
  pred_cells_dead_h_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_hypoxia_obs, rep(0, length(keep_steps)))), step_key)
  pred_cells_dead_b_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_buffer_obs, rep(0, length(keep_steps)))), step_key)
  pred_cells_dead_t_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_total_obs, pred_cells_dead_h_by_step + pred_cells_dead_b_by_step)), step_key)
  pred_vol_total_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Vmm3_total_obs, sim_cpp$Vmm3_obs)), step_key)
  pred_vol_live_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Vmm3_live_obs, sim_cpp$Vmm3_obs)), step_key)
  pred_vol_dead_h_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_hypoxia_obs, rep(0, length(keep_steps)))), step_key)
  pred_vol_dead_b_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_buffer_obs, rep(0, length(keep_steps)))), step_key)
  pred_vol_dead_t_by_step <- setNames(as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_total_obs, pred_vol_dead_h_by_step + pred_vol_dead_b_by_step)), step_key)

  burden_rows <- vector("list", length(keep_steps))
  ploidy_rows <- vector("list", length(keep_steps))
  k <- 0L

  for (step in 0:sim_end_step) {
    Ntot <- sum(v)
    O2_target <- compute_o2_target_from_burden(Ntot = Ntot, run_params = run_params, cfg = cfg)
    if (isTRUE(o2_feedback)) {
      O2_state <- O2_state + alpha_tau * (O2_target - O2_state)
    } else {
      O2_state <- o2_base
    }
    O2_eff <- clip(as.numeric(O2_state), 0, o2_cap)

    if (step %in% keep_steps) {
      k <- k + 1L
      frac_pre <- v[seq_len(R0)]
      frac_post <- v[R0 + seq_len(R1)]
      frac_N <- frac_pre + frac_post
      if (sum(frac_N) > 0) {
        frac_N <- frac_N / sum(frac_N)
      } else {
        frac_N <- rep(1 / length(frac_N), length(frac_N))
      }

      step_chr <- as.character(step)
      pred_burden_cells <- as.numeric(.first_non_null_local(pred_cells_total_by_step[[step_chr]], sum(v)))
      pred_burden_live_cells <- as.numeric(.first_non_null_local(pred_cells_live_by_step[[step_chr]], sum(v)))
      pred_burden_dead_hypoxia_cells <- as.numeric(.first_non_null_local(pred_cells_dead_h_by_step[[step_chr]], 0))
      pred_burden_dead_buffer_cells <- as.numeric(.first_non_null_local(pred_cells_dead_b_by_step[[step_chr]], 0))
      pred_burden_dead_total_cells <- as.numeric(.first_non_null_local(pred_cells_dead_t_by_step[[step_chr]], pred_burden_dead_hypoxia_cells + pred_burden_dead_buffer_cells))
      pred_burden_vol_mm3 <- as.numeric(.first_non_null_local(pred_vol_total_by_step[[step_chr]], NA_real_))
      if (!is.finite(pred_burden_vol_mm3)) {
        pred_burden_vol_mm3 <- burden_volume_mm3_from_state(
          v = v, grid_pre = grid_pre, R0 = R0, R1 = R1,
          run_params = run_params, cfg = cfg, vol_by_N = vol_by_N
        )
      }
      pred_burden_live_vol_mm3 <- as.numeric(.first_non_null_local(pred_vol_live_by_step[[step_chr]], pred_burden_vol_mm3))
      pred_burden_dead_hypoxia_vol_mm3 <- as.numeric(.first_non_null_local(pred_vol_dead_h_by_step[[step_chr]], 0))
      pred_burden_dead_buffer_vol_mm3 <- as.numeric(.first_non_null_local(pred_vol_dead_b_by_step[[step_chr]], 0))
      pred_burden_dead_total_vol_mm3 <- as.numeric(.first_non_null_local(pred_vol_dead_t_by_step[[step_chr]], pred_burden_dead_hypoxia_vol_mm3 + pred_burden_dead_buffer_vol_mm3))
      burden_rows[[k]] <- data.frame(
        harvest = scenario$harvest,
        cohort = scenario$cohort,
        dose = scenario$dose,
        treat_day = scenario$treat_day,
        step = step,
        day = step * cfg$DT,
        pred_burden = pred_burden_vol_mm3,
        pred_burden_volume_mm3 = pred_burden_vol_mm3,
        pred_burden_cells = pred_burden_cells,
        pred_burden_live_volume_mm3 = pred_burden_live_vol_mm3,
        pred_burden_dead_hypoxia_volume_mm3 = pred_burden_dead_hypoxia_vol_mm3,
        pred_burden_dead_buffer_volume_mm3 = pred_burden_dead_buffer_vol_mm3,
        pred_burden_dead_total_volume_mm3 = pred_burden_dead_total_vol_mm3,
        pred_burden_live_cells = pred_burden_live_cells,
        pred_burden_dead_hypoxia_cells = pred_burden_dead_hypoxia_cells,
        pred_burden_dead_buffer_cells = pred_burden_dead_buffer_cells,
        pred_burden_dead_total_cells = pred_burden_dead_total_cells,
        pred_o2_target_pct = O2_target,
        pred_o2_pct = O2_eff,
        pred_o2_lag_gap_pct = O2_target - O2_eff,
        obs_burden = as.numeric(obs_map[as.character(step)]),
        row.names = NULL
      )
      ploidy_rows[[k]] <- data.frame(
        harvest = scenario$harvest,
        cohort = scenario$cohort,
        dose = scenario$dose,
        day = step * cfg$DT,
        N = as.integer(grid_pre),
        fraction = as.numeric(frac_N),
        row.names = NULL
      )
    }

    if (step >= sim_end_step) break

    t <- step * cfg$DT
    tx_mult <- 1.0
    if (isTRUE(cfg$fit_treatment)) {
      tx_mult <- if (t < scenario$treat_day || scenario$dose <= 0) {
        1.0
      } else {
        exp(-run_params$alpha * (dose_scaled^run_params$gamma))
      }
      tx_mult <- clip(tx_mult, cfg$tx_mult_min, 1.0)
    }

    G_step <- if (is.function(get_G)) get_G(O2_eff) else G
    crowd <- if (cfg$crowding == "logistic") max(0, 1 - Ntot / cfg$K) else exp(-Ntot / cfg$K)
    growth_vec <- as.numeric(G_step %*% v)
    v <- v + cfg$DT * (crowd * tx_mult) * growth_vec
    v[!is.finite(v)] <- 0
    v[v < 0] <- 0
    if (sum(v) <= cfg$min_pop) break
  }

  list(
    burden = bind_rows(burden_rows),
    ploidy = bind_rows(ploidy_rows)
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
  o2_grid <- seq(0, 100, by = 0.2)
  ref_df <- data.frame(
    cohort = c("2N", "4N"),
    N_ref = as.numeric(c(2 * cfg$N_UNIT, 4 * cfg$N_UNIT)),
    stringsAsFactors = FALSE
  )

  o2_curve <- dplyr::bind_rows(lapply(seq_len(nrow(ref_df)), function(i) {
    cohort_i <- ref_df$cohort[[i]]
    N_ref <- ref_df$N_ref[[i]]
    ms_rate <- if (exists(".pmisseg_of_O2", mode = "function", inherits = TRUE)) {
      as.numeric(.pmisseg_of_O2(o2_grid, run_params))
    } else {
      k_o_mis_use <- max(as.numeric(run_params$k_o_mis), 1e-12)
      as.numeric(run_params$p_misseg) * (1 - o2_grid / (o2_grid + k_o_mis_use))
    }
    lam_base <- as.numeric(growth_lambda(
      O2 = o2_grid,
      N = N_ref,
      lam_min = run_params$lam_min,
      lam_max = run_params$lam_max,
      k_o = run_params$k_o
    ))
    d_ref <- pmax(0, as.numeric(N_ref) / as.numeric(cfg$N_UNIT) - 2)
    beta_size_use <- as.numeric(run_params$beta_size)
    alpha_o2_use <- pmax(0, as.numeric(run_params$alpha_o2))
    o2_ref_use <- as.numeric(clip(as.numeric(run_params$o2_ref_pct), 0, 100))
    gamma_growth_use <- pmax(as.numeric(run_params$gamma_growth), 1e-12)
    growth_penalty_ploidy_on <- isTRUE(.first_non_null_local(run_params$growth_penalty_ploidy, cfg$growth_penalty_ploidy, FALSE))
    growth_penalty_hypoxia_on <- isTRUE(.first_non_null_local(run_params$growth_penalty_hypoxia, cfg$growth_penalty_hypoxia, FALSE))
    size_penalty <- if (growth_penalty_ploidy_on) exp(-beta_size_use * (d_ref^gamma_growth_use)) else rep(1, length(o2_grid))
    hypoxia_penalty <- if (growth_penalty_hypoxia_on) {
      1 / (1 + alpha_o2_use * d_ref * pmax(0, o2_ref_use - o2_grid))
    } else {
      rep(1, length(o2_grid))
    }
    prolif_rate <- lam_base * size_penalty * hypoxia_penalty
    mu_hp_use <- pmax(as.numeric(run_params$mu_hp), 0)
    death_rate <- mu_hp_use * d_ref * pmax(0, o2_ref_use - o2_grid)
    net_growth_rate <- prolif_rate - death_rate
    data.frame(
      oxygen_pct = o2_grid,
      cohort = cohort_i,
      ms_rate = pmax(ms_rate, 0),
      proliferation_rate = pmax(prolif_rate, 0),
      death_rate = pmax(death_rate, 0),
      net_growth_rate = net_growth_rate,
      N_ref = N_ref,
      row.names = NULL
    )
  }))
  write.table(
    o2_curve,
    file = file.path(out_dir, "functional_curve_oxygen.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  p_ms <- ggplot(o2_curve, aes(x = oxygen_pct, y = ms_rate, color = cohort)) +
    geom_line(linewidth = 1) +
    labs(
      title = "Oxygen vs Missegregation (MS) Rate",
      x = "Oxygen (%)",
      y = "MS rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)

  p_prolif <- ggplot(o2_curve, aes(x = oxygen_pct, y = proliferation_rate, color = cohort)) +
    geom_line(linewidth = 1) +
    labs(
      title = "Oxygen vs Proliferation Rate",
      subtitle = "Proliferation-only component from fitted lambda_eff(N,O2), split by 2N/4N reference ploidy",
      x = "Oxygen (%)",
      y = "Proliferation rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  p_net <- ggplot(o2_curve, aes(x = oxygen_pct, y = net_growth_rate, color = cohort)) +
    geom_line(linewidth = 1) +
    labs(
      title = "Oxygen vs Net Growth Rate",
      subtitle = "Net rate = proliferation - hypoxia-linked high-ploidy death",
      x = "Oxygen (%)",
      y = "Net growth rate",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11)

  ploidy_grid <- seq(cfg$N_MIN / cfg$N_UNIT, cfg$N_MAX / cfg$N_UNIT, by = 0.02)
  N_grid <- pmax(ploidy_grid * cfg$N_UNIT, 1e-8)
  viability <- as.numeric(run_params$smax) * exp(
    -as.numeric(run_params$beta_buffer) * ((2 * cfg$N_UNIT) / N_grid)^as.numeric(run_params$n_exp)
  )
  viability_curve <- data.frame(
    ploidy = ploidy_grid,
    viability_after_ms = pmax(viability, 0),
    row.names = NULL
  )
  write.table(
    viability_curve,
    file = file.path(out_dir, "functional_curve_ploidy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  p_viability <- ggplot(viability_curve, aes(x = ploidy, y = viability_after_ms)) +
    geom_line(color = "#2ca02c", linewidth = 1) +
    labs(
      title = "Ploidy vs Viability After MS",
      subtitle = "Viability term from fitted buffering functional form",
      x = "Ploidy (N / N_UNIT)",
      y = "Viability after MS"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "oxygen_vs_ms_rate.pdf"), p_ms, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_proliferation_rate.pdf"), p_prolif, width = 10, height = 7)
  ggsave(file.path(out_dir, "oxygen_vs_net_growth_rate.pdf"), p_net, width = 10, height = 7)
  ggsave(file.path(out_dir, "ploidy_vs_viability_after_ms.pdf"), p_viability, width = 10, height = 7)
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
      metric = "Burden (normalized)",
      value = as.numeric(pred_norm)
    ) %>%
    bind_rows(
      burden_all %>%
        transmute(
          harvest = as.character(harvest),
          cohort = as.character(cohort),
          dose = as.numeric(dose),
          day = as.numeric(day),
          metric = "Burden (absolute)",
          value = as.numeric(pred_burden)
        )
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
      metric = factor(metric, levels = c("Burden (normalized)", "Burden (absolute)", "Weighted mean ploidy"))
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

  invisible(NULL)
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

  has_o2_lag_cols <- all(c("pred_o2_target_pct", "pred_o2_pct") %in% names(burden_all))
  if (isTRUE(has_o2_lag_cols)) {
    o2_lag_df <- burden_all %>%
      filter(is.finite(pred_o2_target_pct), is.finite(pred_o2_pct)) %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        sample_id = paste(as.character(harvest), as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__"),
        o2_target_pct = as.numeric(pred_o2_target_pct),
        o2_eff_pct = as.numeric(pred_o2_pct),
        o2_lag_gap_pct = as.numeric(pred_o2_target_pct - pred_o2_pct)
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
      labs(
        title = "O2 CBOF MAP Model: Oxygen Lag Over Time",
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
        title = "O2 CBOF MAP Model: O2 Lag Gap Over Time",
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
      title = "O2 CBOF MAP Model: In Vivo Burden Trajectory (Normalized)",
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
      title = "O2 CBOF MAP Model: In Vivo Burden Trajectory (Absolute)",
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
      title = "O2 CBOF MAP Model: In Vivo Burden Trajectory (Absolute, Real Scale)",
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
      title = "O2 CBOF MAP Model: Live/Dead Burden Decomposition",
      subtitle = "Total burden (black) = live + dead from hypoxia + dead from buffer-derived nonviable offspring",
      x = "Day",
      y = "Tumor burden (mm^3)",
      fill = "Component"
    ) +
    theme_bw(base_size = 11)

  o2_burden_df <- burden_all %>%
    filter(is.finite(pred_burden), is.finite(pred_o2_pct)) %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      burden_mm3 = as.numeric(pred_burden),
      o2_pct = as.numeric(pred_o2_pct),
      sample_id = paste(as.character(harvest), as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__")
    )
  write.table(o2_burden_df, file = file.path(out_dir, "predict_burden_vs_o2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  p_burden_vs_o2 <- ggplot(o2_burden_df, aes(x = burden_mm3, y = o2_pct, color = cohort, group = sample_id)) +
    geom_path(linewidth = 0.75, alpha = 0.9) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_x") +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = "O2 CBOF MAP Model: Predicted Oxygen vs Burden",
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
      title = "O2 CBOF MAP Model: Predicted Ploidy Distribution Over Time",
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
      title = paste0("O2 CBOF MAP Model: Ploidy Over Time (Top ", top_n, " N States)"),
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
      title = "O2 CBOF MAP Model: Weighted Mean Ploidy Over Time",
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
  plot_functional_response_curves(run_params = run_params, cfg = cfg, out_dir = out_dir)

  predict_horizons <- as_num_vec(argv$predict_horizons, c(100, 300, 1000))
  predict_horizons <- sort(unique(predict_horizons[is.finite(predict_horizons) & predict_horizons > 0]))
  predict_report_dt <- as_num(argv$predict_report_dt, report_dt)
  if (!is.finite(predict_report_dt) || predict_report_dt <= 0) predict_report_dt <- report_dt
  predict_top_n <- as_int(argv$predict_top_n, top_n)
  if (!is.finite(predict_top_n) || predict_top_n < 1) predict_top_n <- top_n
  do_predict_plots <- as_bool(argv$predict_plots, TRUE)

  if (isTRUE(do_predict_plots) && length(predict_horizons) > 0) {
    for (hz in predict_horizons) {
      message("  Predict plots: 0-", hz, " days (report_dt=", predict_report_dt, ")")
      plot_predict_horizon(
        run_params = run_params,
        scenarios = scenarios,
        cfg = cfg,
        out_dir = out_dir,
        horizon_day = hz,
        report_dt = predict_report_dt,
        top_n = predict_top_n
      )
    }
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
  source(file.path(script_dir, "fit_invivo_model_O2_CBOF_MAP.R"))
  source(file.path(script_dir, "model_O2_CBOF_MAP.R"))

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
