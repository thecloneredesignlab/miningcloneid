find_o2pr_repo_root <- function() {
  d <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:8) {
    script <- file.path(d, "oxygen", "code", "O2_supply_demand_MAP", "analysis", "process_fingerprints", "ploidy_regime_utils.R")
    if (file.exists(script)) return(d)
    nd <- dirname(d)
    if (identical(nd, d)) break
    d <- nd
  }
  normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
}

source_o2pr <- function() {
  root <- find_o2pr_repo_root()
  source(file.path(root, "oxygen", "code", "O2_supply_demand_MAP", "analysis", "process_fingerprints", "process_fingerprint_utils.R"), local = .GlobalEnv)
  source(file.path(root, "oxygen", "code", "O2_supply_demand_MAP", "analysis", "process_fingerprints", "ploidy_regime_utils.R"), local = .GlobalEnv)
}

testthat::test_that("stable synthetic trajectories are labeled stable_high_chr", {
  source_o2pr()
  day <- 0:1000
  curve <- rbind(
    data.frame(seed_id = "seed1", cohort = "2N", day = day, trajectory_value = rep(46, length(day))),
    data.frame(seed_id = "seed1", cohort = "4N", day = day, trajectory_value = rep(48, length(day)))
  )
  feat <- do.call(rbind, lapply(split(curve, curve$cohort), function(d) {
    o2pr_cohort_features("seed1", d$cohort[[1]], d, c(200, 700), c(900, 1000), 22)
  }))
  status <- data.frame(seed_id = "seed1", trajectory_available = TRUE)
  labels <- o2pr_make_trajectory_labels(curve, feat, status, list(N_MIN = 22, N_UNIT = 22, start_with = "chr_number"), 1000, c(200, 700), c(900, 1000))
  testthat::expect_equal(labels$trajectory_regime, "stable_high_chr")
})

testthat::test_that("collapse synthetic trajectories are labeled late_collapse_to_low_chr", {
  source_o2pr()
  day <- 0:1000
  y <- ifelse(day < 700, 50, 50 - (day - 700) * (30 / 300))
  curve <- rbind(
    data.frame(seed_id = "seed1", cohort = "2N", day = day, trajectory_value = y),
    data.frame(seed_id = "seed1", cohort = "4N", day = day, trajectory_value = y + 0.5)
  )
  feat <- do.call(rbind, lapply(split(curve, curve$cohort), function(d) {
    o2pr_cohort_features("seed1", d$cohort[[1]], d, c(200, 700), c(900, 1000), 22)
  }))
  status <- data.frame(seed_id = "seed1", trajectory_available = TRUE)
  labels <- o2pr_make_trajectory_labels(curve, feat, status, list(N_MIN = 22, N_UNIT = 22, start_with = "chr_number"), 1000, c(200, 700), c(900, 1000))
  testthat::expect_equal(labels$trajectory_regime, "late_collapse_to_low_chr")
})

testthat::test_that("common-grid interpolation does not extrapolate", {
  source_o2pr()
  curve <- rbind(
    data.frame(seed_id = "seed1", cohort = "2N", day = 0:10, trajectory_value = 44),
    data.frame(seed_id = "seed1", cohort = "4N", day = 0:10, trajectory_value = 88)
  )
  mat <- o2pr_common_grid_matrix(curve, 1000)$matrix
  testthat::expect_equal(nrow(mat), 0L)
})

testthat::test_that("process clustering matrix excludes trajectory features", {
  source_o2pr()
  x <- data.frame(
    seed_id = rep(paste0("seed", 1:4), each = 2),
    fingerprint_scope = "combined_full",
    module = rep(c("growth_selection", "death_selection"), 4),
    feature = rep(c("lambda_probe", "mu_probe"), 4),
    raw_value = c(1, 4, 2, 3, 3, 2, 4, 1),
    feature_type = "signed",
    stringsAsFactors = FALSE
  )
  sc <- o2pr_scale_long_features_safe(x)
  testthat::expect_false(any(grepl("trajectory|late_level|drop_score", names(sc$wide), ignore.case = TRUE)))
})

testthat::test_that("safe scaling handles all-zero-MAD features", {
  source_o2pr()
  x <- data.frame(
    seed_id = rep("seed1", 2),
    fingerprint_scope = "combined_full",
    module = c("A", "B"),
    feature = c("f1", "f2"),
    raw_value = c(1, 2),
    feature_type = "signed",
    stringsAsFactors = FALSE
  )
  sc <- o2pr_scale_long_features_safe(x)
  testthat::expect_equal(nrow(sc$wide), 0L)
  testthat::expect_true(all(sc$metadata$zero_variance))
})

testthat::test_that("safe scaling keeps seeds with sparse retained conditional features", {
  source_o2pr()
  x <- data.frame(
    seed_id = c(paste0("seed", 1:4), paste0("seed", 1:3)),
    fingerprint_scope = "combined_full",
    module = "oxygen_feedback",
    feature = c(rep("minimum_O2_floor", 4), rep("log10_burden_at_O2_2pct", 3)),
    raw_value = c(0, 0.1, 1, 2, 10, 20, 30),
    feature_type = c(rep("oxygen", 4), rep("signed", 3)),
    stringsAsFactors = FALSE
  )
  sc <- o2pr_scale_long_features_safe(x, missing_feature_max = 0.5, missing_seed_max = 0.5)
  testthat::expect_equal(nrow(sc$wide), 4L)
  testthat::expect_true(all(is.finite(as.matrix(sc$wide[, setdiff(names(sc$wide), "seed_id"), drop = FALSE]))))
  testthat::expect_equal(sc$metadata$missing_fraction[sc$metadata$feature == "log10_burden_at_O2_2pct"], 0.25)
})
