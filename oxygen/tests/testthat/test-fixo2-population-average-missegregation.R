load_fixo2_test_env <- function() {
  simulation_path <- file.path(
    repo_info$root,
    "oxygen", "code", "O2_supply_demand_MAP", "simulation", "fix_o2_simulation.R"
  )
  sim_env <- new.env(parent = globalenv())
  source(simulation_path, local = sim_env, chdir = TRUE)
  sim_env
}

testthat::test_that("fixed-O2 population-average missegregation uses the full state distribution", {
  sim_env <- load_fixo2_test_env()
  model_env <- new.env(parent = baseenv())
  model_env$.pmisseg_of_O2 <- function(O2, run_params, N) {
    run_params$p_mis_base + run_params$p_misseg * (N / max(N)) * exp(-O2)
  }
  state <- c(0.2, 0.3, 0.5)
  ngrid <- c(22, 44, 88)
  params <- list(p_mis_base = 0.001, p_misseg = 0.04)
  state_values <- model_env$.pmisseg_of_O2(rep(1, 3), params, ngrid)

  observed <- sim_env$fixo2_population_average_p_misseg(
    state = state,
    ngrid = ngrid,
    O2 = 1,
    run_params = params,
    model_env = model_env
  )

  testthat::expect_equal(observed, sum(state * state_values), tolerance = 1e-12)
  testthat::expect_gte(observed, 0)
  testthat::expect_lte(observed, 1)
  testthat::expect_equal(sum(sim_env$fixo2_normalize_state(c(2, 3, 5))), 1, tolerance = 1e-12)
})

testthat::test_that("fixed-O2 attractor output keeps a stable schema on failure", {
  sim_env <- load_fixo2_test_env()
  model_env <- new.env(parent = baseenv())
  model_env$.pmisseg_of_O2 <- function(O2, run_params, N) rep(0.02, length(N))
  cfg <- list(N_MIN = 22L, N_UNIT = 22L)
  params <- list()

  sim_env$fixo2_fixed_matrix <- function(model_env, cfg, run_params, O2) {
    list(
      M = diag(c(2, 1)),
      G = diag(c(2, 1)),
      mu_all = c(0, 0),
      ngrid = c(22, 23)
    )
  }
  ok <- sim_env$fixo2_dominant_attractor_one("seed1", params, model_env, cfg, 1)
  testthat::expect_equal(ok$population_average_p_misseg, 0.02, tolerance = 1e-12)

  sim_env$fixo2_fixed_matrix <- function(model_env, cfg, run_params, O2) stop("expected failure")
  failed <- sim_env$fixo2_dominant_attractor_one("seed1", params, model_env, cfg, 1)
  testthat::expect_identical(names(failed), names(ok))
  testthat::expect_true(is.na(failed$population_average_p_misseg))
})
