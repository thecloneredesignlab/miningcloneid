# Fitting-time lineage simulation helpers used by the in-vitro objective.

.ivt_cpp_backend_function <- function(name) {
  # PSOCK workers receive serialized objective closures whose captured sourceCpp
  # functions contain master-process external pointers. Worker initialization
  # loads fresh wrappers into .GlobalEnv, so prefer those process-local bindings.
  if (exists(name, envir = .GlobalEnv, mode = "function", inherits = FALSE)) {
    return(get(name, envir = .GlobalEnv, mode = "function", inherits = FALSE))
  }

  backend_env <- environment(.ivt_cpp_backend_function)
  if (exists(name, envir = backend_env, mode = "function", inherits = TRUE)) {
    return(get(name, envir = backend_env, mode = "function", inherits = TRUE))
  }
  stop("Required process-local C++ backend function is unavailable: ", name)
}

.ivt_cpp_simulate_one <- function(sim_args) {
  .ivt_cpp_backend_function("cpp_o2simps_simulate_one")(sim_args)
}

ivt_set_segment_o2 <- function(target_o2_pct, cfg, run_params) {
  target_o2_use <- as.numeric(target_o2_pct)
  if (!is.finite(target_o2_use)) {
    stop("target_o2_pct must be finite.")
  }
  target_o2_use <- max(0, min(21, target_o2_use))

  sim_cfg <- cfg
  sim_run_params <- run_params
  sim_run_params$o2_S0 <- target_o2_use

  if (!isTRUE(sim_cfg$o2_burden_feedback)) {
    # In the fixed-O2 in vitro path, make the intended constant oxygen explicit.
    sim_cfg$o2_min <- target_o2_use
    sim_run_params$o2_min <- target_o2_use
  }

  list(
    cfg = sim_cfg,
    run_params = sim_run_params,
    target_o2_pct = target_o2_use
  )
}

ivt_build_segment_o2_profile <- function(target_o2_pct, cfg, run_params, burden_grid) {
  o2_setup <- ivt_set_segment_o2(target_o2_pct = target_o2_pct, cfg = cfg, run_params = run_params)
  if (isTRUE(o2_setup$cfg$o2_burden_feedback)) {
    data.frame(
      burden = as.numeric(burden_grid),
      oxygen_target = .o2_supply_demand_from_burden(
        Ntot = burden_grid,
        run_params = o2_setup$run_params,
        o2_Nref = o2_setup$cfg$o2_Nref
      ),
      oxygen_mode = "dynamic feedback",
      target_o2_pct = as.numeric(o2_setup$target_o2_pct),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      burden = as.numeric(burden_grid),
      oxygen_target = rep(as.numeric(o2_setup$target_o2_pct), length(burden_grid)),
      oxygen_mode = paste0("target_O2=", signif(o2_setup$target_o2_pct, 4), "%"),
      target_o2_pct = as.numeric(o2_setup$target_o2_pct),
      stringsAsFactors = FALSE
    )
  }
}

ivt_gaussian_initial_state <- function(grid_pre, mean_N, sd_N) {
  mu <- as.numeric(mean_N)
  sigma <- as.numeric(sd_N)
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
    return(NULL)
  }
  w <- stats::dnorm(as.numeric(grid_pre), mean = mu, sd = sigma)
  w[!is.finite(w)] <- 0
  s <- sum(w)
  if (!is.finite(s) || s <= 0) return(NULL)
  w / s
}

.ivt_initial_state_fraction <- function(cohort, run_params, model_core) {
  cohort <- as.character(cohort)
  if (!cohort %in% c("2N", "4N")) stop("cohort must be '2N' or '4N'.")
  if (identical(cohort, "2N")) {
    gaussian <- ivt_gaussian_initial_state(
      grid_pre = model_core$grid_pre,
      mean_N = .first_non_null_local(run_params$init_mean_2N, NA_real_),
      sd_N = .first_non_null_local(run_params$init_sd_2N, NA_real_)
    )
    base <- if (!is.null(gaussian)) gaussian else as.numeric(model_core$init_state_2N)
  } else {
    gaussian <- ivt_gaussian_initial_state(
      grid_pre = model_core$grid_pre,
      mean_N = .first_non_null_local(run_params$init_mean_4N, NA_real_),
      sd_N = .first_non_null_local(run_params$init_sd_4N, NA_real_)
    )
    base <- if (!is.null(gaussian)) gaussian else as.numeric(model_core$init_state_4N)
  }
  base <- pmax(as.numeric(base), 0)
  total <- sum(base)
  if (!is.finite(total) || total <= 0) {
    return(rep(1 / length(model_core$grid_pre), length(model_core$grid_pre)))
  }
  base / total
}

ivt_run_segment_fixed_o2 <- function(segment,
                                     cfg,
                                     run_params,
                                     model_core,
                                     vol_by_N,
                                     init_state_override = NULL,
                                     init_cells_override = NULL) {
  if (isTRUE(cfg$o2_burden_feedback)) {
    stop(
      "In-vitro passage simulation requires fixed external O2; ",
      "o2_burden_feedback must be FALSE."
    )
  }
  duration_days_use <- suppressWarnings(as.numeric(segment$duration_days))
  obs_days_use <- suppressWarnings(as.numeric(segment$obs_days_local))
  if (length(duration_days_use) != 1L || !is.finite(duration_days_use) ||
      duration_days_use <= 0) {
    stop("segment$duration_days must be one finite positive passage duration.")
  }
  if (length(obs_days_use) == 0L || any(!is.finite(obs_days_use)) ||
      any(obs_days_use < 0) ||
      any(obs_days_use > duration_days_use + 1e-10)) {
    stop("segment$obs_days_local cannot extend beyond the segment search horizon.")
  }
  if (sum(abs(obs_days_use - duration_days_use) <= 1e-10) != 1L) {
    stop("segment$obs_days_local must contain the segment search horizon exactly once.")
  }

  o2_setup <- ivt_set_segment_o2(
    target_o2_pct = segment$oxygen_pct,
    cfg = cfg,
    run_params = run_params
  )
  sim_cfg <- o2_setup$cfg
  rp <- o2_setup$run_params

  if (!is.null(init_state_override)) {
    init_state <- as.numeric(init_state_override)
    if (length(init_state) != length(model_core$grid_pre)) {
      stop("init_state_override length does not match the model chromosome grid.")
    }
    if (any(!is.finite(init_state)) || any(init_state < 0)) {
      stop("init_state_override must contain finite non-negative state values.")
    }
    state_total <- sum(init_state)
    requested_init_cells <- as.numeric(.first_non_null_local(init_cells_override, state_total))
    if (!is.finite(requested_init_cells) || requested_init_cells < 0) {
      stop("init_cells_override must be finite and non-negative.")
    }
    if (!isTRUE(all.equal(state_total, requested_init_cells, tolerance = 1e-12))) {
      stop(
        "init_cells_override cannot rescale an explicit parent state; ",
        "apply the passage-boundary rule before simulation."
      )
    }
    init_cells_use <- state_total
  } else {
    init_cells_use <- as.numeric(.first_non_null_local(
      init_cells_override,
      segment$initial_cells,
      sim_cfg$init_total_size
    ))
    if (!is.finite(init_cells_use) || init_cells_use <= 0) {
      init_cells_use <- as.numeric(sim_cfg$init_total_size)
    }
    init_state <- .ivt_initial_state_fraction(
      cohort = segment$cohort,
      run_params = run_params,
      model_core = model_core
    ) * init_cells_use
  }

  sim <- .ivt_cpp_simulate_one(list(
    init_state = as.numeric(init_state),
    N0min = as.integer(sim_cfg$N_MIN),
    N0max = as.integer(sim_cfg$N_MAX),
    N1min = as.integer(sim_cfg$N_MIN),
    N1max = as.integer(sim_cfg$N_MAX),
    obs_steps = as.integer(round(segment$obs_days_local / sim_cfg$DT)),
    sim_end_step = as.integer(round(segment$duration_days / sim_cfg$DT)),
    DT = as.numeric(sim_cfg$DT),
    dose = 0,
    dose_ref = as.numeric(sim_cfg$dose_ref),
    treat_day = Inf,
    fit_treatment = FALSE,
    alpha = 0,
    gamma = 1,
    tx_mult_min = as.numeric(sim_cfg$tx_mult_min),
    crowding_enabled = isTRUE(sim_cfg$Crowding),
    crowding = as.character(sim_cfg$crowding),
    K = as.numeric(sim_cfg$K),
    min_pop = as.numeric(sim_cfg$min_pop),
    O2_crit = as.numeric(.first_non_null_local(rp$O2_crit, sim_cfg$o2_crit_init, 1.0)),
    o2_feedback = isTRUE(sim_cfg$o2_burden_feedback),
    o2_S0 = as.numeric(rp$o2_S0),
    kappa_O = as.numeric(.first_non_null_local(rp$kappa_O, 1.0)),
    tau_O2 = as.numeric(.first_non_null_local(rp$tau_O2, sim_cfg$tau_O2, sim_cfg$tau_O2_init, 2.0)),
    o2_Nref = as.numeric(sim_cfg$o2_Nref),
    o2_min = as.numeric(.first_non_null_local(rp$o2_min, sim_cfg$o2_min, rp$o2_S0)),
    o2_S0_upper_bound = as.numeric(.first_non_null_local(rp$o2_S0_upper_bound, sim_cfg$o2_S0_upper_bound, 5.0)),
    eta_o2 = as.numeric(.first_non_null_local(rp$eta_o2, 1.0)),
    o2_cache_bin_pct = 0.01,
    o2_cache_hysteresis_pct = 0.005,
    o2_cache_profile = FALSE,
    O2_growth = isTRUE(sim_cfg$O2_growth),
    lam_max = as.numeric(rp$lam_max),
    p_mis_base = as.numeric(.first_non_null_local(rp$p_mis_base, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(rp$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(rp$k_o_mis, 50.0)),
    p_wgd = as.numeric(.first_non_null_local(rp$p_wgd, 0.0)),
    boundary = "drop",
    eps_tail = 1e-8,
    buffer_smax = as.numeric(.first_non_null_local(rp$buffer_smax, sim_cfg$buffer_smax_init, 1.0)),
    buffer_beta = as.numeric(.first_non_null_local(rp$buffer_beta, sim_cfg$buffer_beta_init, 0.0)),
    buffer_n_exp = as.numeric(.first_non_null_local(rp$buffer_n_exp, sim_cfg$buffer_n_exp_init, 1.0)),
    N_unit = as.integer(sim_cfg$N_UNIT),
    beta_size = as.numeric(.first_non_null_local(rp$beta_size, 0.0)),
    alpha_o2 = as.numeric(.first_non_null_local(rp$alpha_o2, sim_cfg$alpha_o2_init, 0.5)),
    gamma_growth = as.numeric(.first_non_null_local(rp$gamma_growth, sim_cfg$gamma_growth_init, 2.0)),
    mu_hp = as.numeric(.first_non_null_local(rp$mu_hp, sim_cfg$mu_hp_init, 1e-3)),
    gamma_mu = as.numeric(.first_non_null_local(rp$gamma_mu, sim_cfg$gamma_mu_init, 1.0)),
    n_O = as.numeric(.first_non_null_local(rp$n_O, sim_cfg$n_O_init, 1.0)),
    ploidy_O2_death = canonical_ploidy_o2_death_mode(
      .first_non_null_local(sim_cfg$ploidy_O2_death, "ploidy_related"),
      "ploidy_related"
    ),
    start_with = assert_canonical_start_with_mode(
      .first_non_null_local(rp$start_with, sim_cfg$start_with, "chr_number")
    ),
    k_clear = as.numeric(.first_non_null_local(rp$k_clear, sim_cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(sim_cfg$burden_log_eps),
    return_full_trajectory = TRUE
  ))
  o2_tolerance <- 1e-10
  o2_target <- suppressWarnings(as.numeric(sim$O2_target_obs))
  o2_eff <- suppressWarnings(as.numeric(sim$O2_eff_obs))
  if (any(!is.finite(o2_target)) ||
      any(!is.finite(o2_eff)) ||
      any(abs(o2_target - o2_setup$target_o2_pct) > o2_tolerance) ||
      any(abs(o2_eff - o2_setup$target_o2_pct) > o2_tolerance)) {
    stop(
      "Fixed-O2 in-vitro segment drifted from its assigned external O2."
    )
  }

  list(
    segment = segment,
    sim = sim,
    initial_state = init_state,
    initial_cells = sum(init_state)
  )
}

ivt_first_threshold_crossing <- function(days,
                                         live_cells,
                                         threshold_target_cells) {
  days_use <- suppressWarnings(as.numeric(days))
  live_use <- suppressWarnings(as.numeric(live_cells))
  if (length(days_use) != length(live_use)) {
    stop("days and live_cells must have the same length.")
  }
  keep <- is.finite(days_use) & is.finite(live_use)
  days_use <- days_use[keep]
  live_use <- live_use[keep]
  if (length(days_use)) {
    ord <- order(days_use)
    days_use <- days_use[ord]
    live_use <- live_use[ord]
  }
  if (anyDuplicated(days_use)) {
    stop("Threshold-crossing time grid must not contain duplicate days.")
  }
  positive_steps <- diff(days_use)
  positive_steps <- positive_steps[is.finite(positive_steps) & positive_steps > 0]
  grid_resolution <- if (length(positive_steps)) {
    max(positive_steps)
  } else {
    NA_real_
  }
  target_use <- suppressWarnings(as.numeric(threshold_target_cells))
  target_valid <- length(target_use) == 1L &&
    is.finite(target_use) &&
    target_use > 0
  empty_result <- list(
    threshold_reached_by_endpoint = FALSE,
    predicted_threshold_crossing_day = NA_real_,
    threshold_time_grid_resolution_days = grid_resolution,
    threshold_crossing_interval_width_days = NA_real_
  )
  if (!target_valid || !length(days_use)) return(empty_result)

  reached_index <- which(live_use >= target_use)
  if (!length(reached_index)) return(empty_result)
  crossing_index <- reached_index[[1L]]
  if (crossing_index == 1L) {
    empty_result$threshold_reached_by_endpoint <- TRUE
    empty_result$predicted_threshold_crossing_day <- days_use[[1L]]
    empty_result$threshold_crossing_interval_width_days <- 0
    return(empty_result)
  }

  left_index <- crossing_index - 1L
  left_day <- days_use[[left_index]]
  right_day <- days_use[[crossing_index]]
  left_cells <- live_use[[left_index]]
  right_cells <- live_use[[crossing_index]]
  interval_width <- right_day - left_day
  crossing_day <- right_day
  if (is.finite(interval_width) &&
      interval_width >= 0 &&
      is.finite(left_cells) &&
      is.finite(right_cells) &&
      left_cells < target_use &&
      right_cells > left_cells) {
    crossing_fraction <- (target_use - left_cells) / (right_cells - left_cells)
    crossing_fraction <- min(max(crossing_fraction, 0), 1)
    crossing_day <- left_day + crossing_fraction * interval_width
  }
  empty_result$threshold_reached_by_endpoint <- TRUE
  empty_result$predicted_threshold_crossing_day <- crossing_day
  empty_result$threshold_crossing_interval_width_days <- interval_width
  empty_result
}

ivt_extract_passage_end_state <- function(sim,
                                         reseed_live_cells,
                                         grid_pre,
                                         target_live_cells = NA_real_,
                                         obs_days_local = NULL,
                                         observed_passage_day = NA_real_,
                                         passage_time_tolerance_days = 1) {
  live_cells <- as.numeric(sim$Ntot_live_obs)
  live_state_mat <- sim$live_state_obs
  obs_n <- length(live_cells)
  if (obs_n == 0L || nrow(live_state_mat) != obs_n) {
    stop("Simulation output does not contain compatible daily trajectories.")
  }
  obs_days_use <- if (is.null(obs_days_local)) {
    seq(0, obs_n - 1L, by = 1)
  } else {
    as.numeric(obs_days_local)
  }
  if (length(obs_days_use) != obs_n) {
    stop("obs_days_local length does not match the number of simulated observations.")
  }

  search_horizon_day <- max(obs_days_use[is.finite(obs_days_use)])
  search_horizon_idx <- which(
    is.finite(obs_days_use) & abs(obs_days_use - search_horizon_day) <= 1e-10
  )
  if (length(search_horizon_idx) != 1L) {
    stop("Simulation output must contain exactly one segment search horizon.")
  }
  search_horizon_idx <- search_horizon_idx[[1]]
  search_horizon_state <- as.numeric(live_state_mat[search_horizon_idx, ])
  search_horizon_live_cells <- sum(search_horizon_state)

  observed_passage_day_use <- suppressWarnings(as.numeric(observed_passage_day))
  if (length(observed_passage_day_use) != 1L ||
      !is.finite(observed_passage_day_use) ||
      observed_passage_day_use < 0) {
    observed_passage_day_use <- search_horizon_day
  }
  last_observation_day_use <- observed_passage_day_use
  last_observation_idx <- which(
    is.finite(obs_days_use) &
      abs(obs_days_use - last_observation_day_use) <= 1e-10
  )
  if (length(last_observation_idx) != 1L) {
    stop("Simulation output must contain the last experimental observation day exactly once.")
  }
  last_observation_idx <- last_observation_idx[[1L]]
  last_observation_live_cells <- as.numeric(live_cells[[last_observation_idx]])
  passage_time_tolerance_days_use <- suppressWarnings(as.numeric(passage_time_tolerance_days))
  if (length(passage_time_tolerance_days_use) != 1L ||
      !is.finite(passage_time_tolerance_days_use) ||
      passage_time_tolerance_days_use < 0) {
    stop("passage_time_tolerance_days must be one finite non-negative value.")
  }

  target_live_cells_use <- suppressWarnings(as.numeric(target_live_cells))
  if (length(target_live_cells_use) != 1L ||
      !is.finite(target_live_cells_use) ||
      target_live_cells_use <= 0) {
    target_live_cells_use <- NA_real_
  }
  required_cells <- suppressWarnings(as.numeric(reseed_live_cells))
  has_boundary <- length(required_cells) == 1L &&
    is.finite(required_cells) && required_cells > 0
  required_cells_use <- if (has_boundary) required_cells else NA_real_
  threshold_values <- c(target_live_cells_use, required_cells_use)
  threshold_values <- threshold_values[is.finite(threshold_values) & threshold_values > 0]
  effective_threshold_cells <- if (length(threshold_values)) max(threshold_values) else NA_real_
  threshold_source_use <- if (is.finite(target_live_cells_use) && has_boundary) {
    if (required_cells_use > target_live_cells_use) {
      "max_observed_final_and_next_initial:next_initial"
    } else {
      "max_observed_final_and_next_initial:observed_final"
    }
  } else if (is.finite(target_live_cells_use)) {
    "observed_final_cells"
  } else if (has_boundary) {
    "next_initial_cells"
  } else {
    "missing"
  }

  positive_day_idx <- which(is.finite(obs_days_use) & obs_days_use > 0)
  candidate_idx <- if (length(positive_day_idx) > 0L) positive_day_idx else seq_len(obs_n)
  passage_window_idx <- candidate_idx[
    obs_days_use[candidate_idx] >= last_observation_day_use - 1e-10 &
      obs_days_use[candidate_idx] <=
        last_observation_day_use + passage_time_tolerance_days_use + 1e-10
  ]
  closest_candidate_idx <- if (length(passage_window_idx)) {
    passage_window_idx
  } else {
    candidate_idx
  }
  closest_idx <- if (is.finite(target_live_cells_use) && target_live_cells_use > 0) {
    ord <- order(
      abs(live_cells[closest_candidate_idx] - target_live_cells_use),
      abs(obs_days_use[closest_candidate_idx] - observed_passage_day_use),
      obs_days_use[closest_candidate_idx],
      closest_candidate_idx
    )
    closest_candidate_idx[[ord[[1]]]]
  } else {
    closest_candidate_idx[[length(closest_candidate_idx)]]
  }
  threshold_crossing <- ivt_first_threshold_crossing(
    days = obs_days_use,
    live_cells = live_cells,
    threshold_target_cells = effective_threshold_cells
  )

  eligible_idx <- if (is.finite(effective_threshold_cells)) {
    passage_window_idx[
      is.finite(live_cells[passage_window_idx]) &
        live_cells[passage_window_idx] >= effective_threshold_cells
    ]
  } else {
    integer()
  }
  max_live_cells_in_search <- if (length(passage_window_idx)) {
    suppressWarnings(max(live_cells[passage_window_idx], na.rm = TRUE))
  } else {
    NA_real_
  }
  if (!is.finite(max_live_cells_in_search)) max_live_cells_in_search <- NA_real_
  passage_executed <- length(eligible_idx) > 0L
  selected_idx <- NA_integer_
  selected_day <- NA_real_
  selected_state <- NULL
  selected_total <- NA_real_
  selected_frac <- rep(NA_real_, length(grid_pre))
  reseeded_state <- NULL
  boundary_scale <- NA_real_
  reseed_mode <- "no_passage_threshold_not_reached"
  passage_failure_reason <- NA_character_
  if (passage_executed) {
    # Passage timing is determined by the closest count from above. A state
    # below the experimental final count is never eligible, even if its
    # absolute count error is smaller.
    eligible_order <- order(
      live_cells[eligible_idx] - effective_threshold_cells,
      abs(obs_days_use[eligible_idx] - observed_passage_day_use),
      obs_days_use[eligible_idx],
      eligible_idx
    )
    selected_idx <- eligible_idx[[eligible_order[[1L]]]]
    selected_day <- as.numeric(obs_days_use[[selected_idx]])
    selected_state <- as.numeric(live_state_mat[selected_idx, ])
    selected_total <- sum(selected_state)
    selected_frac <- selected_state / selected_total
    if (!has_boundary) {
      reseed_mode <- "terminal_no_reseed"
    } else {
      reseed_mode <- "downsample_to_observed_inoculum"
      boundary_scale <- required_cells_use / selected_total
      if (!is.finite(boundary_scale) || boundary_scale < 0 || boundary_scale > 1 + 1e-12) {
        stop("Passage downsampling scale must be finite and no greater than one.")
      }
      boundary_scale <- min(boundary_scale, 1)
      reseeded_state <- selected_state * boundary_scale
    }
  } else {
    passage_failure_reason <- paste0(
      "observed_final_or_next_inoculum_not_reached_after_last_observation; threshold=",
      signif(effective_threshold_cells, 8),
      "; last_observation_day=", signif(last_observation_day_use, 8),
      "; search_horizon_day=", signif(search_horizon_day, 8),
      "; max_live_cells=", signif(max_live_cells_in_search, 8)
    )
  }

  supply_ratio <- if (has_boundary && is.finite(selected_total)) {
    selected_total / required_cells_use
  } else {
    NA_real_
  }
  selected_time_residual <- if (is.finite(selected_day)) {
    selected_day - observed_passage_day_use
  } else {
    NA_real_
  }
  cell_number_after <- if (is.null(reseeded_state)) {
    if (!has_boundary && is.finite(selected_total)) selected_total else NA_real_
  } else {
    sum(reseeded_state)
  }

  list(
    selected_index = selected_idx,
    selected_day = selected_day,
    selected_live_cells = selected_total,
    target_live_cells = if (is.finite(target_live_cells_use) && target_live_cells_use > 0) target_live_cells_use else NA_real_,
    selected_state = selected_state,
    selected_frac = selected_frac,
    endpoint_index = selected_idx,
    endpoint_day = selected_day,
    endpoint_live_cells = selected_total,
    endpoint_state = selected_state,
    endpoint_frac = selected_frac,
    search_horizon_index = search_horizon_idx,
    search_horizon_day = search_horizon_day,
    search_horizon_live_cells = search_horizon_live_cells,
    max_live_cells_in_search = max_live_cells_in_search,
    closest_index_diagnostic = closest_idx,
    closest_day_diagnostic = as.numeric(obs_days_use[[closest_idx]]),
    closest_live_cells_diagnostic = live_cells[[closest_idx]],
    threshold_target_cells = effective_threshold_cells,
    effective_threshold_cells = effective_threshold_cells,
    passage_threshold_cells = target_live_cells_use,
    protocol_threshold_cells = NA_real_,
    observed_final_target_cells = target_live_cells_use,
    threshold_target_source = threshold_source_use,
    threshold_reached_by_endpoint =
      threshold_crossing$threshold_reached_by_endpoint,
    predicted_threshold_crossing_day =
      threshold_crossing$predicted_threshold_crossing_day,
    threshold_time_grid_resolution_days =
      threshold_crossing$threshold_time_grid_resolution_days,
    threshold_crossing_interval_width_days =
      threshold_crossing$threshold_crossing_interval_width_days,
    last_observation_index = last_observation_idx,
    last_observation_day = last_observation_day_use,
    predicted_live_cells_at_observation = last_observation_live_cells,
    observed_passage_day = observed_passage_day_use,
    passage_time_tolerance_days = passage_time_tolerance_days_use,
    passage_time_residual_days = selected_time_residual,
    passage_time_within_tolerance = is.finite(selected_time_residual) &&
      abs(selected_time_residual) <= passage_time_tolerance_days_use + 1e-10,
    selected_day_after_last_observation = is.finite(selected_day) &&
      selected_day >= last_observation_day_use - 1e-10,
    threshold_time_residual_days = if (
      is.finite(threshold_crossing$predicted_threshold_crossing_day)
    ) {
      threshold_crossing$predicted_threshold_crossing_day - observed_passage_day_use
    } else {
      NA_real_
    },
    endpoint_cell_count_residual = if (
      is.finite(target_live_cells_use) &&
        is.finite(selected_total)
    ) {
      selected_total - target_live_cells_use
    } else {
      NA_real_
    },
    cell_count_overshoot = if (is.finite(selected_total) && is.finite(target_live_cells_use)) {
      selected_total - target_live_cells_use
    } else {
      NA_real_
    },
    protocol_threshold_overshoot = if (
      is.finite(selected_total) && is.finite(target_live_cells_use)
    ) {
      selected_total - target_live_cells_use
    } else {
      NA_real_
    },
    passage_executed = passage_executed,
    passage_recorded = passage_executed,
    passage_failure_reason = passage_failure_reason,
    reseed_mode = reseed_mode,
    available_cells = selected_total,
    required_cells = required_cells_use,
    supply_ratio = supply_ratio,
    boundary_scale = boundary_scale,
    cell_number_before = selected_total,
    cell_number_after = cell_number_after,
    reseeded_state = reseeded_state,
    predicted_mean_kary_N = if (passage_executed) {
      ivt_weighted_mean_kary_N(selected_frac, grid_pre = grid_pre)
    } else {
      NA_real_
    }
  )
}

.ivt_stop_protocol_infeasible <- function(segment,
                                          selection,
                                          segment_ordinal,
                                          segment_count) {
  message <- paste0(
    "protocol_infeasible: cohort=", segment$cohort,
    "; scenario=", segment$scenario_id,
    "; segment=", segment$segment_id,
    "; ", selection$passage_failure_reason
  )
  condition <- structure(
    list(
      message = message,
      call = NULL,
      segment = segment,
      selection = selection,
      segment_ordinal = as.integer(segment_ordinal),
      segment_count = as.integer(segment_count)
    ),
    class = c("invitro_protocol_infeasible", "error", "condition")
  )
  stop(condition)
}

ivt_run_lineage <- function(adapter,
                            cfg,
                            run_params,
                            model_core = NULL) {
  if (is.null(model_core)) model_core <- build_model_core(cfg = cfg)
  vol_by_N <- cell_volume_mm3_by_N(model_core$grid_pre, run_params = run_params, cfg = cfg)

  segment_results <- vector("list", length(adapter$segments))
  names(segment_results) <- vapply(adapter$segments, `[[`, character(1), "segment_id")
  if (anyDuplicated(names(segment_results))) {
    stop("In-vitro segment_id values must be unique across scenarios.")
  }

  for (i in seq_along(adapter$segments)) {
    seg <- adapter$segments[[i]]
    parent_key <- seg$parent_segment_id

    parent_res <- NULL
    if (!is.null(parent_key) && is.character(parent_key) && nzchar(parent_key) && !is.na(parent_key)) {
      parent_res <- segment_results[[parent_key]]
      if (is.null(parent_res)) {
        stop("Missing simulated parent segment for ", seg$segment_id, ": ", parent_key)
      }
      if (!identical(as.character(parent_res$segment$scenario_id), as.character(seg$scenario_id)) ||
          !identical(as.character(parent_res$segment$lineage_id), as.character(seg$lineage_id))) {
        stop("Parent state crosses an in-vitro scenario boundary for segment ", seg$segment_id, ".")
      }
    }

    if (!is.null(parent_res)) {
      init_state_override <- as.numeric(parent_res$selection$reseeded_state)
      init_cells_use <- sum(init_state_override)
    } else {
      init_state_override <- NULL
      init_cells_use <- as.numeric(seg$initial_cells)
      if (!is.finite(init_cells_use) || init_cells_use <= 0) {
        init_cells_use <- as.numeric(cfg$init_total_size)
      }
    }

    res <- ivt_run_segment_fixed_o2(
      segment = seg,
      cfg = cfg,
      run_params = run_params,
      model_core = model_core,
      vol_by_N = vol_by_N,
      init_state_override = init_state_override,
      init_cells_override = init_cells_use
    )
    child_idx <- which(vapply(adapter$segments, function(candidate) {
      candidate_parent <- candidate$parent_segment_id
      !is.null(candidate_parent) &&
        length(candidate_parent) == 1L &&
        !is.na(candidate_parent) &&
        identical(as.character(candidate_parent), as.character(seg$segment_id))
    }, logical(1)))
    if (length(child_idx) > 1L) {
      stop("In-vitro scenario is not a linear passage chain at segment ", seg$segment_id, ".")
    }
    next_seg <- if (length(child_idx) == 1L) adapter$segments[[child_idx]] else NULL
    next_init_cells <- if (is.null(next_seg)) NA_real_ else as.numeric(next_seg$initial_cells)
    observed_final_cells <- suppressWarnings(as.numeric(seg$final_cells))
    if (length(observed_final_cells) != 1L ||
        !is.finite(observed_final_cells) ||
        observed_final_cells <= 0) {
      stop(
        "Segment ", seg$segment_id,
        " has no valid observed final live-cell count."
      )
    }
    protocol_threshold_cells <- suppressWarnings(as.numeric(
      seg$protocol_threshold_cells
    ))
    if (length(protocol_threshold_cells) != 1L ||
        !is.finite(protocol_threshold_cells) ||
        protocol_threshold_cells <= 0) {
      stop(
        "Segment ", seg$segment_id,
        " has no valid cohort-level protocol_threshold_cells value."
      )
    }
    protocol_threshold_source <- as.character(seg$protocol_threshold_source)
    if (length(protocol_threshold_source) != 1L ||
        is.na(protocol_threshold_source) ||
        !nzchar(protocol_threshold_source)) {
      stop(
        "Segment ", seg$segment_id,
        " has no valid protocol_threshold_source value."
      )
    }
    segment_time_tolerance <- suppressWarnings(as.numeric(
      .first_non_null_local(
        seg$passage_time_tolerance_days,
        if (is.finite(seg$duration_days) && is.finite(seg$passage_duration)) {
          max(as.numeric(seg$duration_days) - as.numeric(seg$passage_duration), 0)
        } else {
          NULL
        },
        0
      )
    ))
    picked <- ivt_extract_passage_end_state(
      sim = res$sim,
      reseed_live_cells = next_init_cells,
      grid_pre = model_core$grid_pre,
      target_live_cells = observed_final_cells,
      obs_days_local = seg$obs_days_local,
      observed_passage_day = .first_non_null_local(
        seg$last_observation_day,
        seg$observed_passage_day,
        seg$passage_duration
      ),
      passage_time_tolerance_days = segment_time_tolerance
    )
    picked$protocol_threshold_cells <- protocol_threshold_cells
    picked$protocol_threshold_source <- protocol_threshold_source
    picked$passage_threshold_cells <- observed_final_cells
    picked$observed_final_cells <- observed_final_cells
    picked$observed_final_target_cells <- observed_final_cells
    picked$protocol_threshold_overshoot <- if (
      is.finite(picked$selected_live_cells)
    ) {
      picked$selected_live_cells - protocol_threshold_cells
    } else {
      NA_real_
    }
    picked$protocol_threshold_reached_at_selected <-
      is.finite(picked$selected_live_cells) &&
      picked$selected_live_cells >= protocol_threshold_cells
    picked$protocol_threshold_reached_at_observation <-
      is.finite(picked$predicted_live_cells_at_observation) &&
      picked$predicted_live_cells_at_observation >= protocol_threshold_cells
    picked$observed_final_cell_count_residual <- if (
      is.finite(picked$predicted_live_cells_at_observation) &&
        is.finite(observed_final_cells)
    ) {
      picked$predicted_live_cells_at_observation - observed_final_cells
    } else {
      NA_real_
    }
    if (!isTRUE(picked$passage_executed)) {
      .ivt_stop_protocol_infeasible(
        seg,
        picked,
        segment_ordinal = i,
        segment_count = length(adapter$segments)
      )
    }
    res$selection <- picked
    segment_results[[i]] <- res
  }

  list(
    adapter = adapter,
    model_core = model_core,
    grid_pre = model_core$grid_pre,
    segment_results = segment_results,
    initial_observations = adapter$initial_observations,
    landmark_observations = adapter$landmark_observations,
    initial_fraction = .ivt_initial_state_fraction(
      cohort = adapter$cohort,
      run_params = run_params,
      model_core = model_core
    ),
    shared_run_params = run_params
  )
}
