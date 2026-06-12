make_o2_bound_sim_args <- function(o2_upper = 4.0) {
  list(
    init_state = 100.0,
    N0min = 44L,
    N0max = 44L,
    N1min = 44L,
    N1max = 44L,
    obs_steps = 0:4,
    sim_end_step = 4L,
    DT = 0.25,
    dose = 0.0,
    dose_ref = 30.0,
    treat_day = Inf,
    fit_treatment = FALSE,
    alpha = 0.0,
    gamma = 1.0,
    tx_mult_min = 0.0,
    crowding_enabled = FALSE,
    crowding = "logistic",
    K = 1e12,
    min_pop = 1e-12,
    O2_crit = 1.0,
    o2_feedback = TRUE,
    o2_S0 = 25.0,
    kappa_O = 1.0,
    tau_O2 = 1.0,
    o2_Nref = 1e6,
    o2_min = 20.0,
    o2_S0_upper_bound = as.numeric(o2_upper),
    eta_o2 = 1.0,
    o2_cache_bin_pct = 0.001,
    o2_cache_hysteresis_pct = 0.0,
    o2_cache_profile = FALSE,
    lam_max = 0.1,
    p_mis_base = 0.0,
    p_misseg = 0.0,
    k_o_mis = 1.0,
    p_wgd = 0.0,
    boundary = "drop",
    eps_tail = 0.0,
    buffer_smax = 1.0,
    buffer_beta = 0.0,
    buffer_n_exp = 1.0,
    N_unit = 22L,
    beta_size = 0.0,
    O2_growth = FALSE,
    alpha_o2 = 0.0,
    gamma_growth = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    n_O = 1.0,
    ploidy_O2_death = "ploidy_related",
    start_with = "chr_number",
    k_clear = 0.0,
    vol_by_N = 1.0,
    burden_floor = 1e-12,
    return_full_trajectory = TRUE
  )
}

testthat::test_that("cpp_o2simps_o2_window_supply respects o2_S0_upper_bound", {
  out <- cpp_o2simps_o2_window_supply(
    Ntot = c(0, 1e6, 1e12),
    o2_S0 = 25.0,
    kappa_O = 0.5,
    o2_Nref = 1e6,
    o2_min = 20.0,
    o2_S0_upper_bound = 4.0
  )

  testthat::expect_true(all(is.finite(out)))
  testthat::expect_true(all(out <= 4.0 + 1e-12))
  testthat::expect_true(all(out >= 0.0))
  testthat::expect_equal(as.numeric(out[[1]]), 4.0, tolerance = 1e-12)
})

testthat::test_that("cpp_o2simps_simulate_one respects o2_S0_upper_bound for O2 state and targets", {
  upper <- 3.5
  sim <- cpp_o2simps_simulate_one(make_o2_bound_sim_args(o2_upper = upper))

  testthat::expect_true(all(is.finite(sim$O2_target_obs)))
  testthat::expect_true(all(is.finite(sim$O2_eff_obs)))
  testthat::expect_true(all(as.numeric(sim$O2_target_obs) <= upper + 1e-12))
  testthat::expect_true(all(as.numeric(sim$O2_eff_obs) <= upper + 1e-12))
  testthat::expect_equal(as.numeric(sim$O2_target_obs), rep(upper, length(sim$O2_target_obs)), tolerance = 1e-12)
})

testthat::test_that("R and C++ oxygen supply helpers agree on the active upper bound", {
  run_params <- list(
    o2_S0 = 12.0,
    kappa_O = 0.35,
    o2_min = 1.0,
    o2_S0_upper_bound = 4.0
  )
  Ntot <- c(0, 10, 1e6, 1e9)

  r_out <- .o2_supply_demand_from_burden(Ntot, run_params = run_params, o2_Nref = 1e6)
  cpp_out <- cpp_o2simps_o2_window_supply(
    Ntot = Ntot,
    o2_S0 = run_params$o2_S0,
    kappa_O = run_params$kappa_O,
    o2_Nref = 1e6,
    o2_min = run_params$o2_min,
    o2_S0_upper_bound = run_params$o2_S0_upper_bound
  )

  testthat::expect_equal(r_out, as.numeric(cpp_out), tolerance = 1e-12)
  testthat::expect_true(all(r_out <= run_params$o2_S0_upper_bound + 1e-12))
})

testthat::test_that("invalid C++ oxygen upper bounds are rejected", {
  testthat::expect_error(
    cpp_o2simps_o2_window_supply(
      Ntot = 0,
      o2_S0 = 1.0,
      kappa_O = 1.0,
      o2_Nref = 1e6,
      o2_min = 0.0,
      o2_S0_upper_bound = 0.0
    ),
    "o2_S0_upper_bound"
  )
})
