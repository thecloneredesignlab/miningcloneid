make_minimal_o2simps_sim_args <- function(N = 44L,
                                          init_cells = 100.0,
                                          obs_steps = 0:5,
                                          sim_end_step = max(obs_steps),
                                          p_mis_base = 0.0,
                                          p_wgd = 0.0,
                                          mu_hp = 0.0,
                                          boundary = "drop",
                                          o2_S0_upper_bound = 5.0) {
  list(
    init_state = as.numeric(init_cells),
    N0min = as.integer(N),
    N0max = as.integer(N),
    N1min = as.integer(N),
    N1max = as.integer(N),
    obs_steps = as.integer(obs_steps),
    sim_end_step = as.integer(sim_end_step),
    DT = 0.1,
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
    o2_feedback = FALSE,
    o2_S0 = 1.0,
    kappa_O = 1.0,
    tau_O2 = 2.0,
    o2_Nref = 1e6,
    o2_min = 0.0,
    o2_S0_upper_bound = as.numeric(o2_S0_upper_bound),
    eta_o2 = 1.0,
    o2_cache_bin_pct = 0.001,
    o2_cache_hysteresis_pct = 0.0,
    o2_cache_profile = FALSE,
    lam_max = 0.2,
    p_mis_base = as.numeric(p_mis_base),
    p_misseg = 0.0,
    k_o_mis = 1.0,
    p_wgd = as.numeric(p_wgd),
    boundary = as.character(boundary),
    eps_tail = 0.0,
    buffer_smax = 1.0,
    buffer_beta = 0.0,
    buffer_n_exp = 1.0,
    N_unit = 22L,
    beta_size = 0.0,
    O2_growth = FALSE,
    alpha_o2 = 0.0,
    gamma_growth = 1.0,
    mu_hp = as.numeric(mu_hp),
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

testthat::test_that("no-missegregation no-WGD no-death generator has expected live mass gain", {
  # Division generator semantics: lambda is the parent event rate. In the null
  # regime, one event removes one parent and returns two live daughters on the
  # same state, so the live-state generator column sums to +lambda.
  N <- 44L
  lam <- 0.37
  tri <- cpp_o2simps_build_G_for_o2_triplet(
    O2 = 1.0,
    O2_crit = 1.0,
    N0min = N,
    N0max = N,
    N1min = N,
    N1max = N,
    lam_max = lam,
    p_mis_base = 0.0,
    p_misseg = 0.0,
    k_o_mis = 1.0,
    p_wgd = 0.0,
    boundary = "drop",
    eps_tail = 0.0,
    N_unit = 22L,
    beta_size = 0.0,
    O2_growth = FALSE,
    alpha_o2 = 0.0,
    gamma_growth = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    n_O = 1.0,
    ploidy_O2_death = "ploidy_related"
  )
  G <- triplet_to_sparse(tri)

  testthat::expect_equal(as.numeric(G[1, 1]), lam, tolerance = 1e-12)
  testthat::expect_equal(sum(as.numeric(G[, 1])), lam, tolerance = 1e-12)
  testthat::expect_equal(as.numeric(tri$dead_buffer_rate[[1]]), 0.0, tolerance = 1e-12)
})

testthat::test_that("WGD routing matches current parent-transition accounting", {
  # Current WGD semantics: p_wgd diverts a parent division to one live state at
  # 2N with mass lambda*p_wgd, while ordinary two-daughter division is scaled by
  # (1 - p_wgd). It is not encoded as two WGD daughters.
  N <- 30L
  lam <- 2.0
  p_wgd <- 0.25
  tri <- cpp_o2simps_build_G_for_o2_triplet(
    O2 = 1.0,
    O2_crit = 1.0,
    N0min = N,
    N0max = 2L * N,
    N1min = N,
    N1max = 2L * N,
    lam_max = lam,
    p_mis_base = 0.0,
    p_misseg = 0.0,
    k_o_mis = 1.0,
    p_wgd = p_wgd,
    boundary = "drop",
    eps_tail = 0.0,
    N_unit = 22L,
    beta_size = 0.0,
    O2_growth = FALSE,
    alpha_o2 = 0.0,
    gamma_growth = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    n_O = 1.0,
    ploidy_O2_death = "ploidy_related"
  )
  G <- triplet_to_sparse(tri)
  source_idx <- 1L
  wgd_idx <- N + 1L

  testthat::expect_equal(as.numeric(G[source_idx, source_idx]), lam * (1.0 - 2.0 * p_wgd), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(G[wgd_idx, source_idx]), lam * p_wgd, tolerance = 1e-12)
  testthat::expect_equal(sum(as.numeric(G[, source_idx])), lam * (1.0 - p_wgd), tolerance = 1e-12)
})

testthat::test_that("WGD near upper boundary uses documented drop routing", {
  N <- 40L
  lam <- 1.5
  tri <- cpp_o2simps_build_G_for_o2_triplet(
    O2 = 1.0,
    O2_crit = 1.0,
    N0min = 30L,
    N0max = 60L,
    N1min = 30L,
    N1max = 60L,
    lam_max = lam,
    p_mis_base = 0.0,
    p_misseg = 0.0,
    k_o_mis = 1.0,
    p_wgd = 1.0,
    boundary = "drop",
    eps_tail = 0.0,
    N_unit = 22L,
    beta_size = 0.0,
    O2_growth = FALSE,
    alpha_o2 = 0.0,
    gamma_growth = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    n_O = 1.0,
    ploidy_O2_death = "ploidy_related"
  )
  G <- triplet_to_sparse(tri)
  source_idx <- N - 30L + 1L

  testthat::expect_equal(sum(as.numeric(G[, source_idx])), -lam, tolerance = 1e-12)
  testthat::expect_equal(as.numeric(tri$boundary_dropped_rate[[source_idx]]), lam, tolerance = 1e-12)
  testthat::expect_equal(as.numeric(tri$dead_buffer_rate[[source_idx]]), lam, tolerance = 1e-12)
})

testthat::test_that("mu_hp zero keeps dead-hypoxia compartment identically zero", {
  sim <- cpp_o2simps_simulate_one(make_minimal_o2simps_sim_args(
    p_mis_base = 0.0,
    p_wgd = 0.0,
    mu_hp = 0.0
  ))

  testthat::expect_equal(as.numeric(sim$Ntot_dead_hypoxia_obs), rep(0.0, length(sim$Ntot_dead_hypoxia_obs)), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(sim$Ntot_dead_buffer_obs), rep(0.0, length(sim$Ntot_dead_buffer_obs)), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(sim$Ntot_dead_total_obs), rep(0.0, length(sim$Ntot_dead_total_obs)), tolerance = 1e-12)
  testthat::expect_equal(as.matrix(sim$dead_hypoxia_state_obs), matrix(0.0, nrow = length(sim$Ntot_dead_hypoxia_obs), ncol = 1L), tolerance = 1e-12)
  testthat::expect_equal(as.matrix(sim$dead_buffer_state_obs), matrix(0.0, nrow = length(sim$Ntot_dead_buffer_obs), ncol = 1L), tolerance = 1e-12)
})

testthat::test_that("with mu_hp zero dead growth occurs only through nonviable daughter routing", {
  sim <- cpp_o2simps_simulate_one(make_minimal_o2simps_sim_args(
    p_mis_base = 0.15,
    p_wgd = 0.0,
    mu_hp = 0.0,
    boundary = "drop"
  ))

  dead_h <- as.numeric(sim$Ntot_dead_hypoxia_obs)
  dead_b <- as.numeric(sim$Ntot_dead_buffer_obs)
  dead_total <- as.numeric(sim$Ntot_dead_total_obs)

  testthat::expect_equal(dead_h, rep(0.0, length(dead_h)), tolerance = 1e-12)
  testthat::expect_gt(max(dead_b), 0.0)
  testthat::expect_equal(dead_total, dead_b, tolerance = 1e-12)
  testthat::expect_equal(as.matrix(sim$dead_hypoxia_state_obs), matrix(0.0, nrow = length(dead_h), ncol = 1L), tolerance = 1e-12)
})
