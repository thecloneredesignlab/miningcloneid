ivt_run_segment_fixed_o2 <- function(segment,
                                     cfg,
                                     run_params,
                                     model_core,
                                     vol_by_N,
                                     init_state_override = NULL,
                                     init_cells_override = NULL) {
  rp <- run_params
  rp$o2_S0 <- as.numeric(segment$oxygen_pct)

  init_state_base <- if (!is.null(init_state_override)) {
    as.numeric(init_state_override)
  } else if (segment$cohort == "2N") {
    as.numeric(model_core$init_state_2N)
  } else {
    as.numeric(model_core$init_state_4N)
  }
  init_cells_use <- as.numeric(.first_non_null_local(
    init_cells_override,
    segment$initial_cells,
    cfg$init_total_size
  ))
  if (!is.finite(init_cells_use) || init_cells_use <= 0) init_cells_use <- as.numeric(cfg$init_total_size)
  init_frac <- if (sum(init_state_base) > 0) init_state_base / sum(init_state_base) else rep(1 / length(init_state_base), length(init_state_base))
  init_state <- init_frac * init_cells_use

  sim <- cpp_o2simps_simulate_one(
    init_state = init_state,
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = as.integer(round(segment$obs_days_local / cfg$DT)),
    sim_end_step = as.integer(round(segment$duration_days / cfg$DT)),
    DT = as.numeric(cfg$DT),
    dose = 0,
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = Inf,
    fit_treatment = FALSE,
    alpha = 0,
    gamma = 1,
    tx_mult_min = as.numeric(cfg$tx_mult_min),
    crowding_enabled = isTRUE(cfg$Crowding),
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_crit = as.numeric(.first_non_null_local(rp$O2_crit, cfg$o2_crit_init, 1.0)),
    o2_feedback = FALSE,
    o2_S0 = as.numeric(segment$oxygen_pct),
    kappa_O = as.numeric(.first_non_null_local(rp$kappa_O, 1.0)),
    tau_O2 = as.numeric(.first_non_null_local(rp$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0)),
    o2_Nref = as.numeric(cfg$o2_Nref),
    o2_min = as.numeric(cfg$o2_min),
    eta_o2 = as.numeric(.first_non_null_local(rp$eta_o2, 1.0)),
    o2_cache_bin_pct = 0.01,
    o2_cache_hysteresis_pct = 0.005,
    o2_cache_profile = FALSE,
    O2_growth = isTRUE(cfg$O2_growth),
    lam_min = as.numeric(rp$lam_min),
    lam_max = as.numeric(rp$lam_max),
    k_o = as.numeric(rp$k_o),
    has_p_misseg = !is.null(rp$p_misseg),
    p_mis_base = as.numeric(.first_non_null_local(rp$p_mis_base, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(rp$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(rp$k_o_mis, 50.0)),
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = 0.0,
    p_wgd = as.numeric(.first_non_null_local(rp$p_wgd, 0.0)),
    boundary = "drop",
    eps_tail = 1e-8,
    gamma_loss = as.numeric(.first_non_null_local(rp$gamma_loss, 0.1)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = as.numeric(.first_non_null_local(rp$beta_size, default_beta_size_prior_center())),
    alpha_o2 = as.numeric(.first_non_null_local(rp$alpha_o2, cfg$alpha_o2_init, 0.5)),
    gamma_growth = as.numeric(.first_non_null_local(rp$gamma_growth, cfg$gamma_growth_init, 2.0)),
    mu_hp = as.numeric(.first_non_null_local(rp$mu_hp, cfg$mu_hp_init, 1e-3)),
    gamma_mu = as.numeric(.first_non_null_local(rp$gamma_mu, cfg$gamma_mu_init, 1.0)),
    n_O = as.numeric(.first_non_null_local(rp$n_O, cfg$n_O_init, 1.0)),
    ploidy_O2_death = canonical_ploidy_o2_death_mode(
      .first_non_null_local(cfg$ploidy_O2_death, "diploid_NULL"),
      "diploid_NULL"
    ),
    k_clear = as.numeric(.first_non_null_local(rp$k_clear, cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(cfg$burden_log_eps),
    return_full_trajectory = TRUE
  )

  list(
    segment = segment,
    sim = sim
  )
}

ivt_extract_passage_end_state <- function(sim, reseed_live_cells, grid_pre) {
  live_cells <- as.numeric(sim$Ntot_live_obs)
  live_state_mat <- sim$live_state_obs
  obs_n <- length(live_cells)
  if (obs_n == 0L || nrow(live_state_mat) != obs_n) {
    stop("Simulation output does not contain compatible daily trajectories.")
  }
  if (!is.finite(reseed_live_cells) || reseed_live_cells <= 0) {
    stop("reseed_live_cells must be positive.")
  }

  idx <- obs_n
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
    selected_day = as.numeric(idx - 1),
    selected_live_cells = live_cells[[idx]],
    selected_frac = chosen_frac,
    reseeded_state = reseeded_state,
    predicted_mean_ploidy = ivt_weighted_mean_ploidy(chosen_frac, grid_pre = grid_pre)
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
    next_init_cells <- as.numeric(seg$initial_cells)
    if (!is.finite(next_init_cells) || next_init_cells <= 0) {
      next_init_cells <- as.numeric(cfg$init_total_size)
    }
    picked <- ivt_extract_passage_end_state(
      sim = res$sim,
      reseed_live_cells = next_init_cells,
      grid_pre = model_core$grid_pre
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
