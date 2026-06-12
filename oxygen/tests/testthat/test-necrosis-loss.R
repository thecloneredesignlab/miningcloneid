testthat::test_that("necrosis observations are grouped by mapped harvest", {
  path <- tempfile(fileext = ".csv")
  tab <- data.frame(
    dt_harvest = c("H1", "H1", "H2", "H3", ""),
    mapping_status = c("mapped", "mapped", "unmapped", "mapped", "mapped"),
    percent_necrosis = c("50", "70", "90", "", "40"),
    stringsAsFactors = FALSE
  )
  write.csv(tab, path, row.names = FALSE, quote = TRUE)

  obs <- read_necrosis_observations(path, use_necrosis_loss = TRUE)

  testthat::expect_equal(nrow(obs), 1L)
  testthat::expect_equal(obs$harvest[[1]], "H1")
  testthat::expect_equal(obs$obs_necrosis_fraction[[1]], 0.6, tolerance = 1e-12)
  testthat::expect_equal(obs$n_necrosis_obs[[1]], 2L)
})

testthat::test_that("missing necrosis harvests are skipped as NA observations", {
  obs_tab <- data.frame(
    harvest = "H1",
    obs_necrosis_fraction = 0.6,
    n_necrosis_obs = 2L,
    stringsAsFactors = FALSE
  )

  present <- necrosis_observation_for_harvest(obs_tab, "H1")
  missing <- necrosis_observation_for_harvest(obs_tab, "H2")

  testthat::expect_equal(present$obs_necrosis_fraction, 0.6, tolerance = 1e-12)
  testthat::expect_equal(present$n_necrosis_obs, 2L)
  testthat::expect_true(is.na(missing$obs_necrosis_fraction))
  testthat::expect_equal(missing$n_necrosis_obs, 0L)
})

make_necrosis_objective_inputs <- function(use_necrosis_loss = TRUE,
                                           obs_necrosis_fraction = 0.25,
                                           keep_necrosis = TRUE) {
  scenario_data <- list(
    cohort_code = as.integer(0L),
    dose_vec = 0.0,
    treat_day_vec = Inf,
    obs_steps_list = list(as.integer(c(0L, 10L))),
    sim_end_step_vec = as.integer(10L),
    obs_burden_list = list(c(NA_real_, NA_real_)),
    keep_burden_list = list(c(FALSE, FALSE)),
    ploidy_z_list = list(numeric(0)),
    init_mult_vec = 1.0,
    obs_necrosis_fraction_vec = as.numeric(obs_necrosis_fraction),
    keep_necrosis_vec = as.logical(keep_necrosis),
    necrosis_step_vec = as.integer(10L)
  )
  objective_data <- list(
    mu_by_N = 44.0,
    sigma_burden = 0.1,
    sigma_ploidy = 0.1,
    use_necrosis_loss = isTRUE(use_necrosis_loss),
    sigma_necrosis_logit = 1.0,
    necrosis_fraction_eps = 1e-4
  )
  state_data <- list(
    init_state_2N = 100.0,
    init_state_4N = 100.0,
    N0min = as.integer(44L),
    N0max = as.integer(44L),
    N1min = as.integer(44L),
    N1max = as.integer(44L),
    N_unit = as.integer(22L),
    vol_by_N = 1.0
  )
  sim_args <- list(
    DT = 0.1,
    dose_ref = 30.0,
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
    o2_S0_upper_bound = 5.0,
    eta_o2 = 1.0,
    o2_cache_bin_pct = 0.001,
    o2_cache_hysteresis_pct = 0.0,
    o2_cache_profile = FALSE,
    lam_max = 0.0,
    p_mis_base = 0.0,
    p_misseg = 0.0,
    k_o_mis = 1.0,
    p_wgd = 0.0,
    boundary = "drop",
    eps_tail = 0.0,
    buffer_smax = 1.0,
    buffer_beta = 0.0,
    buffer_n_exp = 1.0,
    beta_size = 0.0,
    alpha_o2 = 0.0,
    gamma_growth = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    n_O = 1.0,
    ploidy_O2_death = "ploidy_related",
    start_with = "chr_number",
    k_clear = 0.0,
    burden_log_eps = 1e-12
  )

  list(
    scenario_data = scenario_data,
    objective_data = objective_data,
    state_data = state_data,
    sim_args = sim_args
  )
}

testthat::test_that("necrosis objective contributes only when enabled and observed", {
  inputs <- make_necrosis_objective_inputs(use_necrosis_loss = TRUE)
  out <- cpp_o2simps_objective_components_map(
    scenario_data = inputs$scenario_data,
    objective_data = inputs$objective_data,
    state_data = inputs$state_data,
    sim_args = inputs$sim_args
  )

  expected <- (qlogis(1e-4) - qlogis(0.25))^2
  testthat::expect_equal(out$n_necrosis, 1L)
  testthat::expect_equal(out$L_n, expected, tolerance = 1e-8)

  disabled_inputs <- make_necrosis_objective_inputs(use_necrosis_loss = FALSE)
  disabled <- cpp_o2simps_objective_components_map(
    scenario_data = disabled_inputs$scenario_data,
    objective_data = disabled_inputs$objective_data,
    state_data = disabled_inputs$state_data,
    sim_args = disabled_inputs$sim_args
  )
  testthat::expect_equal(disabled$n_necrosis, 0L)
  testthat::expect_equal(disabled$L_n, 0, tolerance = 1e-12)

  missing_inputs <- make_necrosis_objective_inputs(
    use_necrosis_loss = TRUE,
    obs_necrosis_fraction = NA_real_,
    keep_necrosis = FALSE
  )
  missing <- cpp_o2simps_objective_components_map(
    scenario_data = missing_inputs$scenario_data,
    objective_data = missing_inputs$objective_data,
    state_data = missing_inputs$state_data,
    sim_args = missing_inputs$sim_args
  )
  testthat::expect_equal(missing$n_necrosis, 0L)
  testthat::expect_equal(missing$L_n, 0, tolerance = 1e-12)
})
