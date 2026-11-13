# Pure mapping/orchestration tests; no simulator or production optimization.
if (!exists("o2prof_map", mode = "function")) {
  profile_root <- if (exists("repo_info")) repo_info$root else Sys.getenv("MININGCLONEID_REPO_ROOT")
  profile_map_module <- file.path(profile_root, "oxygen", "code", "O2_supply_demand_MAP", "util", "o2_supply_demand_map_profile.R")
  source(profile_map_module, local = environment())
}
if (!exists("profile_test_table")) {
  profile_root <- if (exists("repo_info")) repo_info$root else Sys.getenv("MININGCLONEID_REPO_ROOT")
  profile_test_table <- read.csv(file.path(profile_root, "oxygen", "data", "O2_supply_demand", "parameter_table_invitro_buffering.csv"), stringsAsFactors = FALSE)
  profile_test_params <- setNames(lapply(profile_test_table$init_value, as.numeric), profile_test_table$param_symbol)
  profile_test_params$ploidy_O2_death <- "ploidy_related"
}
testthat::test_that("production and zero-capable maps reconstruct full parameters", {
  for (domain in c("production", "zero_extended")) {
    map <- o2prof_map(profile_test_table, profile_test_params, domain)
    free <- o2prof_to_free(map, profile_test_params)
    recovered <- o2prof_from_free(map, free)
    testthat::expect_length(free, 20L)
    testthat::expect_equal(unlist(recovered[o2prof_symbols]), unlist(profile_test_params[o2prof_symbols]), tolerance = 1e-12)
    testthat::expect_identical(recovered$o2_S0, profile_test_params$o2_S0)
    testthat::expect_identical(recovered$ploidy_O2_death, profile_test_params$ploidy_O2_death)
  }
  zero <- o2prof_map(profile_test_table, profile_test_params, "zero_extended")
  p <- profile_test_params; p$p_misseg <- 0; p$buffer_beta <- 0
  encoded <- o2prof_to_free(zero, p)
  testthat::expect_identical(unname(encoded[c("log1p_p_misseg", "log1p_buffer_beta")]), c(0, 0))
  testthat::expect_identical(o2prof_from_free(zero, encoded)$p_misseg, 0)
  testthat::expect_identical(o2prof_from_free(zero, encoded)$buffer_beta, 0)
  testthat::expect_error(o2prof_map(profile_test_table, profile_test_params, "production", fixed = list(p_misseg = 0)), "outside")
})

testthat::test_that("exact nulls remove only completely inactive nuisance coordinates", {
  pmap <- o2prof_map(profile_test_table, profile_test_params, "zero_extended", fixed = list(p_misseg = 0))
  bmap <- o2prof_map(profile_test_table, profile_test_params, "zero_extended", fixed = list(buffer_beta = 0))
  for (m in list(pmap, bmap)) {
    full <- o2prof_project_start(m, profile_test_params)$params
    testthat::expect_length(o2prof_to_free(m, full), 18L)
    for (n in c("p_mis_base", "p_wgd", "sigma_growth", "sigma_kary", "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N")) {
      testthat::expect_identical(m$table$status[m$table$parameter == n], "active")
    }
  }
  testthat::expect_identical(pmap$table$status[pmap$table$parameter == "k_o_mis"], "inactive")
  testthat::expect_identical(bmap$table$status[bmap$table$parameter == "buffer_n_exp"], "inactive")
  testthat::expect_identical(bmap$table$status[bmap$table$parameter == "buffer_smax"], "active")
  penalty <- o2prof_map(profile_test_table, profile_test_params, "zero_extended", fixed = list(p_misseg = 0), penalty_active = "k_o_mis")
  testthat::expect_identical(penalty$table$status[penalty$table$parameter == "k_o_mis"], "active")
  testthat::expect_false(penalty$table$data_active[penalty$table$parameter == "k_o_mis"])
  testthat::expect_length(o2prof_to_free(penalty, o2prof_project_start(penalty, profile_test_params)$params), 19L)
  both <- o2prof_map(profile_test_table, profile_test_params, "zero_extended", fixed = list(p_misseg = 0, buffer_beta = 0))
  testthat::expect_equal(sum(both$table$status == "active"), 16L)
})

testthat::test_that("production log-bound roundoff is recorded without relaxing exact targets", {
  map <- o2prof_map(profile_test_table, profile_test_params, "zero_extended", fixed = list(p_misseg = 0))
  start <- profile_test_params
  lower <- map$table$lower[map$table$parameter == "mu_hp"]
  start$mu_hp <- lower * (1 - 2 * .Machine$double.eps)
  projected <- o2prof_project_start(map, start)
  testthat::expect_identical(projected$params$mu_hp, lower)
  testthat::expect_identical(projected$projection$reason[projected$projection$parameter == "mu_hp"], "roundoff_to_bound")
  testthat::expect_identical(projected$params$p_misseg, 0)
  start$mu_hp <- lower / 2
  testthat::expect_error(o2prof_project_start(map, start), "outside")
})

testthat::test_that("fixed zeros survive persistence and real PSOCK serialization", {
  map <- o2prof_map(profile_test_table, profile_test_params, "zero_extended", fixed = list(p_misseg = 0))
  path <- tempfile(fileext = ".rds"); saveRDS(map, path)
  restored <- readRDS(path)
  free <- o2prof_to_free(restored, o2prof_project_start(restored, profile_test_params)$params)
  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  answer <- parallel::clusterCall(cl, function(module, map, free) {
    sys.source(module, envir = .GlobalEnv)
    x <- o2prof_from_free(map, free)
    list(amplitude = x$p_misseg, k = x$k_o_mis)
  }, profile_map_module, restored, free)
  for (x in answer) testthat::expect_identical(x$amplitude, 0)
  unlink(path)
})

testthat::test_that("complete starts reach the initial population without hidden substitution", {
  map <- o2prof_map(profile_test_table, profile_test_params, "zero_extended", fixed = list(buffer_beta = 0))
  second <- profile_test_params; second$lam_max <- second$lam_max * 0.8; second$init_mean_4N <- 85
  set.seed(101)
  pop <- o2prof_initial_population(map, 32L, list(baseline = profile_test_params, neighbor = second))
  for (i in 1:2) testthat::expect_identical(unname(pop$population[i, ]), unname(o2prof_to_free(map, pop$starts[[i]]$params)))
  reconstructed <- o2prof_from_free(map, setNames(pop$population[2, ], colnames(pop$population)))
  testthat::expect_equal(reconstructed$init_mean_4N, 85)
  testthat::expect_identical(reconstructed$buffer_beta, 0)
  testthat::expect_identical(o2prof_initial_population(map, 32L, list(baseline = profile_test_params, neighbor = second), pop$population)$population, pop$population)
  bad <- pop$population; bad[1, 1] <- bad[1, 1] + 0.01
  testthat::expect_error(o2prof_initial_population(map, 32L, list(baseline = profile_test_params, neighbor = second), bad), "declared injected")
  testthat::expect_error(o2prof_initial_population(map, 32L, list(profile_test_params)), "named")
})

testthat::test_that("constrained nuisance optimization matches an analytic profile", {
  for (t in c(0, 1, 2, 3, 4)) {
    fit <- stats::optim(0, function(u) (t + u - 3)^2 + (u - 1)^2, method = "L-BFGS-B", lower = -4, upper = 4)
    testthat::expect_equal(unname(fit$par), 2 - t / 2, tolerance = 1e-6)
    testthat::expect_equal(fit$value, (t - 2)^2 / 2, tolerance = 1e-9)
  }
})

testthat::test_that("missing/nonfinite scored observations cannot disappear from averages", {
  manifest <- data.frame(observation_id = c("a", "b"), growth_scored = c(TRUE, TRUE), karyotype_scored = c(TRUE, FALSE), flow_scored = c(FALSE, TRUE))
  comp <- list(objective = 1, n_growth_missing_pred = 0L,
               growth_df = data.frame(passage_id = c("a", "b"), loglik = c(-1, -2)),
               ploidy_df = data.frame(passage_id = "a", mean_loglik = -1),
               flow_df = data.frame(passage_id = "b", mean_loglik = -1))
  testthat::expect_silent(o2prof_check_samples(comp, manifest))
  missing <- comp; missing$growth_df <- missing$growth_df[1, ]
  testthat::expect_error(o2prof_check_samples(missing, manifest), "sample contract")
  nonfinite <- comp; nonfinite$flow_df$mean_loglik <- NA_real_
  testthat::expect_error(o2prof_check_samples(nonfinite, manifest), "sample contract")
})

testthat::test_that("candidate acceptance and resume reject stale/wrong/interrupted results", {
  map <- o2prof_map(profile_test_table, profile_test_params, "zero_extended", fixed = list(p_misseg = 0))
  candidate <- list(fingerprint = "runA", launch_completed = TRUE, launch_exit_code = 0L,
                    initial_population_hash = "populationA", params = o2prof_project_start(map, profile_test_params)$params,
                    samples_valid = TRUE, objective = 2, reevaluated_objective = 2, budget_limited = TRUE)
  testthat::expect_identical(o2prof_accept_candidate(candidate, "runA", map, "populationA")$search_status, "budget_limited_provisional")
  testthat::expect_error(o2prof_accept_candidate(candidate, "runB", map, "populationA"), "fingerprint")
  testthat::expect_error(o2prof_accept_candidate(candidate, "runA", map, "populationB"), "population")
  wrong <- candidate; wrong$params$p_misseg <- 1e-9
  testthat::expect_error(o2prof_accept_candidate(wrong, "runA", map, "populationA"), "fixed/inactive")
  stopped <- candidate; stopped$launch_completed <- FALSE
  testthat::expect_error(o2prof_accept_candidate(stopped, "runA", map, "populationA"), "Interrupted")
  bad <- candidate; bad$reevaluated_objective <- 3
  testthat::expect_error(o2prof_accept_candidate(bad, "runA", map, "populationA"), "validation")
  path <- tempfile(); dir.create(path)
  jsonlite::write_json(list(fingerprint = "runA"), file.path(path, "run_manifest.json"), auto_unbox = TRUE)
  testthat::expect_silent(o2prof_assert_resume(path, "runA"))
  testthat::expect_error(o2prof_assert_resume(path, "runB"), "fingerprint")
  writeLines("{}", file.path(path, "validated_result.json"))
  testthat::expect_error(o2prof_assert_resume(path, "runA"), "overwritten")
  unlink(path, recursive = TRUE)
})
