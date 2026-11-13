# Run by the isolated v2 validation entry point, after loading frozen ctx.
if (!exists("profile_validation_root")) testthat::skip("Requires the isolated frozen-v2 validation entry point")
testthat::test_that("integrated production and extended maps reproduce the accepted endpoint", {
  for (domain in c("production", "zero_extended")) {
    map <- o2prof_map(profile_test_table, ctx$saved$best_params, domain)
    full <- o2prof_from_free(map, o2prof_to_free(map, ctx$saved$best_params))
    o2prof_install_culture_map(ctx, map, list(seed144 = full))
    comp <- ctx$backend$ivt_objective_components(full, ctx$objects, ctx$cfg)
    testthat::expect_equal(comp$objective, ctx$saved$best_components$objective, tolerance = 1e-8)
    testthat::expect_equal(o2prof_independent_score(ctx, comp)$objective, comp$objective, tolerance = 1e-8)
    for (cohort in c("2N", "4N")) {
      before <- ctx$saved$best_components[[paste0("run_", cohort)]]
      after <- comp[[paste0("run_", cohort)]]
      testthat::expect_identical(names(after$segment_results), names(before$segment_results))
      for (i in seq_along(after$segment_results)) {
        a <- before$segment_results[[i]]; b <- after$segment_results[[i]]
        testthat::expect_identical(a$selection$selected_day, b$selection$selected_day)
        testthat::expect_identical(a$segment$data_ids, b$segment$data_ids)
        difference <- abs(a$sim$live_state_obs - b$sim$live_state_obs)
        bound <- 1e-8 + 1e-8 * pmax(abs(a$sim$live_state_obs), abs(b$sim$live_state_obs))
        testthat::expect_true(all(difference <= bound))
        testthat::expect_equal(a$selection$selected_frac, b$selection$selected_frac, tolerance = 1e-12)
        if (!is.null(a$parent_boundary_selection)) testthat::expect_identical(a$parent_boundary_selection$selected_day, b$parent_boundary_selection$selected_day)
      }
    }
  }
})

testthat::test_that("exact nulls survive production normalization and kernel dispatch", {
  for (target in c("p_misseg", "buffer_beta")) {
    fixed <- setNames(list(0), target)
    map <- o2prof_map(profile_test_table, ctx$saved$best_params, "zero_extended", fixed = fixed)
    full <- o2prof_project_start(map, ctx$saved$best_params)$params
    o2prof_install_culture_map(ctx, map, list(seed144 = full))
    normalized <- ctx$backend$ivt_optim_par_to_run_params(o2prof_to_free(map, full), ctx$cfg)
    testthat::expect_identical(normalized[[target]], 0)
    testthat::expect_identical(ctx$cfg, ctx$backend$build_invitro_cfg(ctx$cfg$parameter_table))
    testthat::expect_identical(ctx$cfg$prior_center_log10_buffer_beta, ctx$saved$cfg$prior_center_log10_buffer_beta)
  }
  build <- function(O, mu, k, amp = 0) {
    x <- cpp_o2simps_build_G_for_o2_triplet(
      O2 = O, O2_crit = 1, N0min = 22L, N0max = 154L, N1min = 22L, N1max = 154L,
      lam_max = 0.4, p_mis_base = 0.01, p_misseg = amp, k_o_mis = k,
      p_wgd = 0.001, boundary = "drop", eps_tail = 0,
      buffer_smax = 0.8, buffer_beta = 0.5, buffer_n_exp = 1.2,
      N_unit = 22L, O2_growth = FALSE, mu_hp = mu, gamma_mu = 1,
      ploidy_O2_death = "ploidy_related")
    x <- as.matrix(triplet_to_sparse(x)); diag(x) <- 0; x
  }
  reference <- build(20.5, 1e-4, 1e-8)
  testthat::expect_gt(sum(reference), 0)
  for (O in c(0, 1, 20.5)) for (k in c(1e-8, 0.01, 1)) {
    testthat::expect_equal(build(O, 0.1, k), reference, tolerance = 0)
  }
  testthat::expect_gt(max(abs(build(0, 0.1, 1e-4, amp = 0.03) - reference)), 0)
  for (N in c(44L, 88L)) {
    a <- cpp_o2simps_pr_delta_vec(N, p = 0.03, eps_tail = 0, buffer_smax = 0.8, buffer_beta = 0, buffer_n_exp = 1, N_unit = 22L)
    b <- cpp_o2simps_pr_delta_vec(N, p = 0.03, eps_tail = 0, buffer_smax = 0.8, buffer_beta = 0, buffer_n_exp = 9, N_unit = 22L)
    testthat::expect_identical(a$prob, b$prob)
    testthat::expect_equal(shift_weight(a, 0L), 2 * stats::dbinom(0, N, 0.03), tolerance = 1e-12)
    testthat::expect_equal(shift_weight(a, 1L) / stats::dbinom(1, N, 0.03), 0.8, tolerance = 1e-12)
    testthat::expect_equal(shift_weight(a, 2L) / stats::dbinom(2, N, 0.03), 0.8^2, tolerance = 1e-12)
  }
})

testthat::test_that("production C++ initialization and exact zeros reach serialized workers", {
  map <- o2prof_map(profile_test_table, ctx$saved$best_params, "zero_extended", fixed = list(p_misseg = 0, buffer_beta = 0))
  full <- o2prof_project_start(map, ctx$saved$best_params)$params
  o2prof_install_culture_map(ctx, map, list(seed144 = full))
  free <- o2prof_to_free(map, full)
  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  ctx$backend$initialize_invitro_deoptim_workers(cl, ctx$backend$o2simps_cpp_dll_info())
  answer <- parallel::clusterCall(cl, function(convert, free, cfg) {
    p <- convert(free, cfg)
    d <- cpp_o2simps_pr_delta_vec(44L, 0.03, eps_tail = 0, buffer_smax = p$buffer_smax,
                                 buffer_beta = p$buffer_beta, buffer_n_exp = p$buffer_n_exp)
    list(p = p$p_misseg, beta = p$buffer_beta, finite_kernel = all(is.finite(d$prob)))
  }, ctx$backend$ivt_optim_par_to_run_params, free, ctx$cfg)
  for (a in answer) {
    testthat::expect_identical(a$p, 0)
    testthat::expect_identical(a$beta, 0)
    testthat::expect_true(a$finite_kernel)
  }
})

testthat::test_that("native production DEoptim receives the saved population and complete starts", {
  map <- o2prof_map(profile_test_table, ctx$saved$best_params, "zero_extended", fixed = list(p_misseg = 0))
  second <- ctx$saved$best_params; second$init_mean_4N <- 85
  directory <- tempfile("native_toy_fit_", tmpdir = profile_validation_root)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  o2prof_install_culture_map(ctx, map, list(seed144 = ctx$saved$best_params, neighbor = second), directory)
  # Only this one test replaces the simulator with a cheap analytic objective.
  ctx$backend$ivt_objective_components <- function(run_params, fit_objects, cfg, ...) {
    if (!identical(run_params$p_misseg, 0)) stop("Native optimizer lost the exact zero constraint")
    x <- ctx$saved$best_components
    x$objective <- (run_params$lam_max - ctx$saved$best_params$lam_max)^2
    x$total_loglik <- -x$objective
    x
  }
  cleanup <- o2prof_trace_native_population(directory)
  on.exit(cleanup(), add = TRUE)
  ctx$backend$main(list(parameter_table = ctx$cfg$parameter_table,
                        fit_objects_dir = file.path(ctx$wp0_root, "source/inputs/oxygen/ploidyOxygen/data/fit_objects"),
                        flow_density_path = file.path(ctx$wp0_root, "source/inputs/oxygen/data/g0g1_ploidy_density_grid.csv"),
                        out_dir = directory, seed = 13L, itermax = 2L, itermax_max = 2L,
                        n_cores = 1L, local_optim_maxit = 2L, auto_viz = FALSE))
  native <- readRDS(file.path(directory, "native_de_initial_population.rds"))
  requested <- readRDS(file.path(directory, "requested_initial_population.rds"))
  testthat::expect_identical(native, requested)
  testthat::expect_equal(dim(native), c(256L, 18L))
  for (i in 1:2) {
    decoded <- o2prof_from_free(map, setNames(native[i, ], colnames(native)))
    testthat::expect_identical(decoded$p_misseg, 0)
    if (i == 2) testthat::expect_equal(decoded$init_mean_4N, 85)
  }
  result <- readRDS(file.path(directory, "fit_result.rds"))
  testthat::expect_identical(result$best_params$p_misseg, 0)
  testthat::expect_identical(result$cfg, ctx$cfg)
  testthat::expect_equal(result$best_components$objective, 0, tolerance = 1e-12)
})
