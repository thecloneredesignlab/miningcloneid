if (!exists("pilot_study")) testthat::skip("Requires the isolated culture-pilot validation entry")

testthat::test_that("saved minima survive serialization reset and retain complete vectors", {
  folder <- tempfile("checkpoint_", tmpdir = profile_validation_root)
  p <- ctx$saved$best_params
  testthat::expect_true(o2pilot_checkpoint(p, 5, folder, "run_a"))
  callback <- function(q) o2pilot_checkpoint(p, q, folder, "run_a")
  reset <- unserialize(serialize(callback, NULL))
  testthat::expect_false(reset(6))
  rm(list = ls(.o2pilot_best_cache, all.names = TRUE), envir = .o2pilot_best_cache)
  testthat::expect_false(reset(7))
  testthat::expect_true(reset(4))
  saved <- readRDS(list.files(folder, full.names = TRUE)[1])
  testthat::expect_equal(saved$objective, 4)
  testthat::expect_identical(saved$params, p)
  testthat::expect_error(o2pilot_checkpoint(p, 3, folder, "stale"), "Stale")
  testthat::expect_error(o2pilot_checkpoint(p, Inf, folder, "run_a"), "Nonfinite")
})

testthat::test_that("serialized workers retain their lowest feasible candidate on disk", {
  folder <- tempfile("psock_best_", tmpdir = profile_validation_root)
  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(paths) { for (p in paths) sys.source(p, .GlobalEnv); TRUE }, pilot_modules)
  for (q in c(5, 3, 8)) parallel::clusterCall(cl, function(params, value, folder) {
    o2pilot_checkpoint(params, value, folder, "psock")
  }, ctx$saved$best_params, q, folder)
  candidates <- lapply(list.files(folder, full.names = TRUE), readRDS)
  testthat::expect_length(candidates, 2)
  testthat::expect_true(all(vapply(candidates, function(x) x$objective == 3, logical(1))))
  for (a in candidates) testthat::expect_identical(a$params, ctx$saved$best_params)
})

testthat::test_that("native v2 receives 20 and 19 dimensional starts without changing the objective contract", {
  for (fixed in list(list(), list(p_misseg = 0.01), list(buffer_beta = 1))) {
    toy <- o2prof_load_culture(file.path(project_root, settings$wp0_relative_path), pilot_modules)
    map <- o2prof_map(toy$cfg$parameter_table, toy$saved$best_params, "zero_extended", fixed = fixed)
    directory <- tempfile("native_pilot_toy_", tmpdir = profile_validation_root)
    dir.create(directory)
    o2prof_install_culture_map(toy, map, list(actual_start = toy$saved$best_params), directory)
    toy$backend$ivt_objective_components <- function(run_params, fit_objects, cfg, ...) {
      o2prof_to_free(map, run_params)
      x <- toy$saved$best_components
      x$objective <- (run_params$lam_max - toy$saved$best_params$lam_max)^2
      x$total_loglik <- -x$objective
      x
    }
    cleanup <- o2prof_trace_native_population(directory)
    tryCatch(toy$backend$main(list(parameter_table = toy$cfg$parameter_table,
                                   fit_objects_dir = file.path(toy$wp0_root, "source/inputs/oxygen/ploidyOxygen/data/fit_objects"),
                                   flow_density_path = file.path(toy$wp0_root, "source/inputs/oxygen/data/g0g1_ploidy_density_grid.csv"),
                                   out_dir = directory, seed = 31L, itermax = 2L, itermax_max = 2L,
                                   n_cores = 1L, local_optim_maxit = 2L, auto_viz = FALSE)), finally = cleanup())
    actual <- readRDS(file.path(directory, "native_de_initial_population.rds"))
    testthat::expect_identical(actual, readRDS(file.path(directory, "requested_initial_population.rds")))
    testthat::expect_identical(ncol(actual), if (length(fixed)) 19L else 20L)
    decoded <- o2prof_from_free(map, setNames(actual[1, ], colnames(actual)))
    for (name in names(fixed)) testthat::expect_identical(decoded[[name]], fixed[[name]])
    result <- readRDS(file.path(directory, "fit_result.rds"))
    testthat::expect_identical(result$cfg, toy$cfg)
    testthat::expect_equal(result$best_components$objective, 0, tolerance = 1e-12)
  }
})
