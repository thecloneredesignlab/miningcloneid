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

ivt_run_segment_fixed_o2 <- function(segment,
                                     cfg,
                                     run_params,
                                     model_core,
                                     vol_by_N,
                                     init_state_override = NULL,
                                     init_cells_override = NULL) {
  o2_setup <- ivt_set_segment_o2(
    target_o2_pct = segment$oxygen_pct,
    cfg = cfg,
    run_params = run_params
  )
  sim_cfg <- o2_setup$cfg
  rp <- o2_setup$run_params

  init_state_base <- if (!is.null(init_state_override)) {
    as.numeric(init_state_override)
  } else if (segment$cohort == "2N") {
    init_gauss <- ivt_gaussian_initial_state(
      grid_pre = model_core$grid_pre,
      mean_N = .first_non_null_local(run_params$init_mean_2N, NA_real_),
      sd_N = .first_non_null_local(run_params$init_sd_2N, NA_real_)
    )
    if (!is.null(init_gauss)) init_gauss else as.numeric(model_core$init_state_2N)
  } else {
    init_gauss <- ivt_gaussian_initial_state(
      grid_pre = model_core$grid_pre,
      mean_N = .first_non_null_local(run_params$init_mean_4N, NA_real_),
      sd_N = .first_non_null_local(run_params$init_sd_4N, NA_real_)
    )
    if (!is.null(init_gauss)) init_gauss else as.numeric(model_core$init_state_4N)
  }
  init_cells_use <- as.numeric(.first_non_null_local(
    init_cells_override,
    segment$initial_cells,
    sim_cfg$init_total_size
  ))
  if (!is.finite(init_cells_use) || init_cells_use <= 0) init_cells_use <- as.numeric(sim_cfg$init_total_size)
  init_frac <- if (sum(init_state_base) > 0) init_state_base / sum(init_state_base) else rep(1 / length(init_state_base), length(init_state_base))
  init_state <- init_frac * init_cells_use

  sim <- cpp_o2simps_simulate_one(list(
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

  list(
    segment = segment,
    sim = sim
  )
}

ivt_extract_passage_end_state <- function(sim,
                                         reseed_live_cells,
                                         grid_pre,
                                         target_live_cells = NA_real_,
                                         obs_days_local = NULL) {
  live_cells <- as.numeric(sim$Ntot_live_obs)
  live_state_mat <- sim$live_state_obs
  obs_n <- length(live_cells)
  if (obs_n == 0L || nrow(live_state_mat) != obs_n) {
    stop("Simulation output does not contain compatible daily trajectories.")
  }
  if (!is.finite(reseed_live_cells) || reseed_live_cells <= 0) {
    stop("reseed_live_cells must be positive.")
  }
  obs_days_use <- if (is.null(obs_days_local)) {
    seq(0, obs_n - 1L, by = 1)
  } else {
    as.numeric(obs_days_local)
  }
  if (length(obs_days_use) != obs_n) {
    stop("obs_days_local length does not match the number of simulated observations.")
  }

  target_live_cells_use <- as.numeric(target_live_cells)
  positive_day_idx <- which(is.finite(obs_days_use) & obs_days_use > 0)
  candidate_idx <- if (length(positive_day_idx) > 0L) positive_day_idx else seq_len(obs_n)

  idx <- if (is.finite(target_live_cells_use) && target_live_cells_use > 0) {
    ord <- order(abs(live_cells[candidate_idx] - target_live_cells_use), candidate_idx)
    candidate_idx[[ord[[1]]]]
  } else {
    candidate_idx[[length(candidate_idx)]]
  }
  chosen_state <- as.numeric(live_state_mat[idx, ])
  chosen_total <- sum(chosen_state)
  chosen_frac <- if (is.finite(chosen_total) && chosen_total > 0) {
    chosen_state / chosen_total
  } else {
    rep(1 / length(grid_pre), length(grid_pre))
  }
  reseeded_state <- chosen_frac * reseed_live_cells

  list(
    selected_index = idx,
    selected_day = as.numeric(obs_days_use[[idx]]),
    selected_live_cells = live_cells[[idx]],
    target_live_cells = if (is.finite(target_live_cells_use) && target_live_cells_use > 0) target_live_cells_use else NA_real_,
    selected_frac = chosen_frac,
    reseeded_state = reseeded_state,
    predicted_mean_kary_N = ivt_weighted_mean_kary_N(chosen_frac, grid_pre = grid_pre)
  )
}

ivt_run_lineage <- function(adapter,
                            cfg,
                            run_params,
                            model_core = NULL) {
  if (is.null(model_core)) model_core <- build_model_core(cfg = cfg)
  vol_by_N <- cell_volume_mm3_by_N(model_core$grid_pre, run_params = run_params, cfg = cfg)

  segment_results <- vector("list", length(adapter$segments))
  names(segment_results) <- vapply(adapter$segments, `[[`, character(1), "segment_id")

  for (i in seq_along(adapter$segments)) {
    seg <- adapter$segments[[i]]
    parent_key <- seg$parent_segment_id

    parent_res <- NULL
    if (!is.null(parent_key) && is.character(parent_key) && nzchar(parent_key) && !is.na(parent_key)) {
      parent_res <- segment_results[[parent_key]]
    }

    init_cells_use <- as.numeric(seg$initial_cells)
    if (!is.finite(init_cells_use) || init_cells_use <= 0) {
      init_cells_use <- as.numeric(cfg$init_total_size)
    }

    init_state_override <- if (!is.null(parent_res)) {
      as.numeric(parent_res$selection$selected_frac) * init_cells_use
    } else if (i > 1L) {
      prev_res <- segment_results[[i - 1L]]
      if (!is.null(prev_res)) {
        as.numeric(prev_res$selection$selected_frac) * init_cells_use
      } else {
        NULL
      }
    } else {
      NULL
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
    next_seg <- if (i < length(adapter$segments)) adapter$segments[[i + 1L]] else NULL
    target_live_cells_use <- as.numeric(seg$final_cells)
    if ((!is.finite(target_live_cells_use) || target_live_cells_use <= 0) && !is.null(next_seg)) {
      target_live_cells_use <- as.numeric(next_seg$initial_cells)
    }
    next_init_cells <- as.numeric(seg$initial_cells)
    if (!is.finite(next_init_cells) || next_init_cells <= 0) {
      next_init_cells <- as.numeric(cfg$init_total_size)
    }
    picked <- ivt_extract_passage_end_state(
      sim = res$sim,
      reseed_live_cells = next_init_cells,
      grid_pre = model_core$grid_pre,
      target_live_cells = target_live_cells_use,
      obs_days_local = seg$obs_days_local
    )
    res$selection <- picked
    segment_results[[i]] <- res
  }

  list(
    adapter = adapter,
    model_core = model_core,
    grid_pre = model_core$grid_pre,
    segment_results = segment_results
  )
}
