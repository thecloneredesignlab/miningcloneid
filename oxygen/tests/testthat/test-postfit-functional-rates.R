postfit_workflow_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

postfit_seed10_params <- list(
  lam_max = 1.37317637942057,
  p_mis_base = 0.005,
  p_misseg = 0.0278817962094089,
  k_o_mis = 1e-5,
  p_wgd = 0.001,
  boundary = "drop",
  buffer_smax = 1,
  buffer_beta = 1.21454513799383,
  buffer_n_exp = 0.886203818141274,
  alpha_o2 = 4,
  gamma_growth = 3.4698976487065,
  mu_hp = 0.005,
  gamma_mu = 1.2,
  O2_crit = 0.00178169753233183,
  n_O = 3.81806920430669
)

postfit_seed10_cfg <- list(
  start_with = "chr_number",
  O2_growth = TRUE,
  ploidy_O2_death = "ploidy_related",
  N_MIN = 22L,
  N_MAX = 154L,
  N_UNIT = 22L,
  o2_crit_init = 1,
  n_O_init = 1,
  mu_hp_init = 0.001,
  gamma_mu_init = 1,
  buffer_smax_init = 0.8,
  buffer_beta_init = 1,
  buffer_n_exp_init = 1,
  alpha_o2_init = 0.5,
  gamma_growth_init = 2
)

testthat::test_that("in-vitro post-fit net growth includes dead-buffer loss", {
  response_env <- new.env(parent = globalenv())
  sys.source(
    file.path(
      postfit_workflow_root,
      "util",
      "o2_supply_demand_map_invitro_postfit_io_utils.R"
    ),
    envir = response_env
  )
  sys.source(
    file.path(
      postfit_workflow_root,
      "util",
      "o2_supply_demand_map_invitro_postfit_response_utils.R"
    ),
    envir = response_env
  )

  context <- response_env$ivt_sim_build_functional_response_context(
    postfit_seed10_params,
    postfit_seed10_cfg
  )
  rates <- context$compute_rates(o2 = c(0, 0), n = c(44, 88))

  testthat::expect_equal(
    rates$net_growth_rate,
    rates$proliferation_rate - rates$death_rate - rates$buffer_death_rate,
    tolerance = 1e-12
  )
  testthat::expect_true(all(rates$net_growth_rate < 0))
})

testthat::test_that("in-vivo post-fit net growth uses the same live-state balance", {
  response_env <- new.env(parent = globalenv())
  sys.source(
    file.path(
      postfit_workflow_root,
      "simulation",
      "invivo",
      "o2",
      "invivo_o2_simulation.R"
    ),
    envir = response_env
  )
  sys.source(
    file.path(
      postfit_workflow_root,
      "simulation",
      "invivo",
      "functional_response",
      "invivo_functional_response_simulation.R"
    ),
    envir = response_env
  )

  tables <- response_env$generate_invivo_functional_response_tables(
    postfit_seed10_params,
    postfit_seed10_cfg
  )
  rates <- tables$functional_curve_oxygen_multi_ploidy
  rates <- rates[rates$oxygen_pct == 0 & rates$N_ref %in% c(44, 88), ]

  testthat::expect_equal(nrow(rates), 2L)
  testthat::expect_equal(
    rates$net_growth_rate,
    rates$proliferation_rate - rates$death_rate - rates$buffer_death_rate,
    tolerance = 1e-12
  )
  testthat::expect_true(all(rates$net_growth_rate < 0))
})
