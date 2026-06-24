find_o2ipa_repo_root <- function() {
  d <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:8) {
    script <- file.path(d, "oxygen", "code", "O2_supply_demand_MAP", "analysis", "process_fingerprint_utils.R")
    if (file.exists(script)) return(d)
    nd <- dirname(d)
    if (identical(nd, d)) break
    d <- nd
  }
  normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
}

o2ipa_test_script <- function() {
  file.path(find_o2ipa_repo_root(), "oxygen", "code", "O2_supply_demand_MAP", "analysis", "process_fingerprint_utils.R")
}

testthat::test_that("parameter transform follows optimizer scale", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  testthat::expect_equal(o2ipa_transform_parameter_value("lam_max", 0.1), -1)
  testthat::expect_equal(o2ipa_transform_parameter_value("buffer_smax", 0.8), 0.8)
  testthat::expect_true(is.na(o2ipa_transform_parameter_value("p_wgd", 0)))
})

testthat::test_that("alias fallback reads seed summary values", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  summary_row <- data.frame(seed_id = "seed1", value__o2_crit = 0.2, value__alpha_o2 = 1.5, check.names = FALSE)
  hit <- o2ipa_extract_param("seed1", "O2_crit", setNames(numeric(0), character(0)), summary_row, NULL)
  testthat::expect_equal(hit$value, 0.2)
  testthat::expect_equal(hit$alias, "o2_crit")
})

testthat::test_that("zero MAD features are removed from clustering matrix", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  x <- data.frame(
    seed_id = rep(paste0("seed", 1:4), each = 2),
    fingerprint_scope = "static_full",
    module = rep(c("A", "B"), 4),
    feature = rep(c("constant", "variable"), 4),
    raw_value = c(1, 1, 1, 2, 1, 3, 1, 4),
    feature_type = "signed",
    stringsAsFactors = FALSE
  )
  scaled <- o2ipa_scale_long_features(x)
  testthat::expect_false(any(grepl("constant", names(scaled$wide), fixed = TRUE)))
  testthat::expect_true(any(grepl("variable", names(scaled$wide), fixed = TRUE)))
  testthat::expect_true(any(scaled$metadata$zero_variance[scaled$metadata$feature == "constant"]))
})

testthat::test_that("conditional feature missingness is counted across all seeds", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  x <- data.frame(
    seed_id = c(paste0("seed", 1:4), paste0("seed", 1:3), "seed1"),
    fingerprint_scope = "static_full",
    module = "A",
    feature = c(rep("stable", 4), rep("conditional", 3), "mostly_missing"),
    raw_value = c(1, 2, 3, 4, 10, 20, 30, 100),
    feature_type = "signed",
    stringsAsFactors = FALSE
  )
  scaled <- o2ipa_scale_long_features(x, missing_feature_max = 0.5, missing_seed_max = 0.5)
  meta <- scaled$metadata
  testthat::expect_equal(meta$missing_fraction[meta$feature == "conditional"], 0.25)
  testthat::expect_equal(meta$missing_fraction[meta$feature == "mostly_missing"], 0.75)
  testthat::expect_true(meta$retained_for_clustering[meta$feature == "conditional"])
  testthat::expect_false(meta$retained_for_clustering[meta$feature == "mostly_missing"])
})

testthat::test_that("retained conditional values are center-imputed in clustering matrix", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  x <- data.frame(
    seed_id = c(paste0("seed", 1:4), paste0("seed", 1:3)),
    fingerprint_scope = "static_full",
    module = "A",
    feature = c(rep("stable", 4), rep("conditional", 3)),
    raw_value = c(1, 2, 3, 4, 10, 20, 30),
    feature_type = "signed",
    stringsAsFactors = FALSE
  )
  scaled <- o2ipa_scale_long_features(x, missing_feature_max = 0.5, missing_seed_max = 0.5)
  key <- "static_full||A||conditional"
  testthat::expect_true(key %in% names(scaled$wide))
  testthat::expect_equal(nrow(scaled$wide), 4L)
  testthat::expect_true(all(is.finite(scaled$wide[[key]])))
  testthat::expect_equal(scaled$wide[[key]][scaled$wide$seed_id == "seed4"], 0)
})

testthat::test_that("oxygen zero is a finite floor value after transform", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  transformed <- o2ipa_feature_transform(c(0, 1), c("oxygen", "oxygen"))
  testthat::expect_true(all(is.finite(transformed)))
  testthat::expect_lt(transformed[[1]], transformed[[2]])
})

testthat::test_that("module weighting balances retained features by module", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  x <- data.frame(
    seed_id = rep(paste0("seed", 1:4), each = 3),
    fingerprint_scope = "static_full",
    module = rep(c("A", "A", "B"), 4),
    feature = rep(c("a1", "a2", "b1"), 4),
    raw_value = c(1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6),
    feature_type = "signed",
    stringsAsFactors = FALSE
  )
  scaled <- o2ipa_scale_long_features(x)
  weights <- unique(scaled$long[, c("module", "module_weight")])
  testthat::expect_equal(weights$module_weight[weights$module == "A"][[1]], 1 / sqrt(2))
  testthat::expect_equal(weights$module_weight[weights$module == "B"][[1]], 1)
})

testthat::test_that("module centroid table is generated from retained features and membership", {
  script <- file.path(find_o2ipa_repo_root(), "oxygen", "code", "O2_supply_demand_MAP", "analysis", "run_invivo_process_analysis.R")
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  x <- data.frame(
    seed_id = rep(paste0("seed", 1:4), each = 2),
    fingerprint_scope = "static_full",
    module = rep(c("A", "B"), 4),
    feature = rep(c("a1", "b1"), 4),
    raw_value = c(1, 10, 2, 20, 4, 40, 5, 50),
    feature_type = "signed",
    stringsAsFactors = FALSE
  )
  scaled <- o2ipa_scale_long_features(x)
  membership <- data.frame(
    seed_id = paste0("seed", 1:4),
    k = 2L,
    cluster = c(1L, 1L, 2L, 2L),
    stringsAsFactors = FALSE
  )
  cent <- o2ipa_module_centroids(scaled$long, membership, 2L, "static_full")
  testthat::expect_equal(nrow(cent), 4L)
  testthat::expect_equal(sort(unique(cent$module)), c("A", "B"))
  testthat::expect_equal(sort(unique(cent$cluster)), c(1L, 2L))
  testthat::expect_true(all(cent$n_seed == 2L))
  testthat::expect_true(all(cent$n_features == 1L))
})

testthat::test_that("distance matrix is symmetric with zero diagonal", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  wide <- data.frame(seed_id = paste0("seed", 1:3), f1 = c(0, 1, 2), f2 = c(1, 1, 3))
  d <- o2ipa_distance(wide)
  testthat::expect_equal(unname(diag(d)), rep(0, 3))
  testthat::expect_equal(d, t(d))
})

testthat::test_that("clustering is deterministic for fixed seed", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  wide <- data.frame(seed_id = paste0("seed", 1:5), f1 = c(0, 0.1, 5, 5.1, 9), f2 = c(0, 0.1, 5, 5.1, 9))
  d <- o2ipa_distance(wide)
  a <- o2ipa_cluster_distance(d, n_bootstrap = 3, random_seed = 11, feature_wide = wide)
  b <- o2ipa_cluster_distance(d, n_bootstrap = 3, random_seed = 11, feature_wide = wide)
  testthat::expect_equal(a$membership, b$membership)
  testthat::expect_equal(a$recommended_k, b$recommended_k)
})

testthat::test_that("bad seed is logged without stopping discovery", {
  script <- o2ipa_test_script()
  testthat::skip_if_not(file.exists(script))
  source(script, local = TRUE)
  tmp <- tempfile("o2ipa_mock_")
  o2ipa_make_mock_run(tmp, n_seed = 3)
  file.remove(file.path(tmp, "seed2", "best_params.tsv"))
  inputs <- o2ipa_collect_seed_inputs(tmp, objective_source = "total")
  testthat::expect_equal(nrow(inputs$manifest), 3)
  testthat::expect_true(any(!(inputs$manifest$fit_success %in% TRUE)))
  unlink(tmp, recursive = TRUE, force = TRUE)
})
