find_joint_ploidy_repo_root <- function() {
  current <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:8) {
    target <- file.path(current, "oxygen", "code", "O2_supply_demand_MAP", "util", "o2_supply_demand_map_joint_ploidy_category_utils.R")
    if (file.exists(target)) return(current)
    parent <- dirname(current); if (identical(parent, current)) break; current <- parent
  }
  stop("Cannot locate repository root")
}

source_joint_ploidy_utils <- function() {
  root <- find_joint_ploidy_repo_root()
  source(file.path(root, "oxygen", "code", "O2_supply_demand_MAP", "util", "o2_supply_demand_map_joint_ploidy_category_utils.R"), local = .GlobalEnv)
  root
}

synthetic_joint_ploidy_curves <- function() {
  day <- 0:1000
  list(
    day = day,
    CatA = rep(44, length(day)),
    CatB = ifelse(day < 250, 44, ifelse(day <= 450, 44 - (day - 250) * 22 / 200, 22)),
    CatC = ifelse(
      day < 150, 44,
      ifelse(day <= 300, 44 - (day - 150) * 9 / 150, ifelse(day < 650, 35, ifelse(day <= 800, 35 - (day - 650) * 13 / 150, 22)))
    )
  )
}

testthat::test_that("synthetic one-stage, two-stage, and high trajectories map to CatB, CatC, and CatA", {
  source_joint_ploidy_utils()
  curves <- synthetic_joint_ploidy_curves()
  testthat::expect_equal(o2jpc_curve_features(curves$day, curves$CatA)$category, "CatA")
  testthat::expect_equal(o2jpc_curve_features(curves$day, curves$CatB)$category, "CatB")
  catc <- o2jpc_curve_features(curves$day, curves$CatC)
  testthat::expect_equal(catc$category, "CatC")
  testthat::expect_gte(catc$plateau_duration, 60)
  testthat::expect_lte(catc$delta_bic_two_minus_one, -10)
})

testthat::test_that("44 chromosomes is retained in CatA using numerical tolerance", {
  source_joint_ploidy_utils()
  day <- 0:1000
  curve <- rep(44, length(day)); curve[451:550] <- 43.75
  testthat::expect_equal(o2jpc_curve_features(day, curve, high_tolerance = 0.5)$category, "CatA")
  testthat::expect_false(o2jpc_curve_features(day, curve, high_tolerance = 0.1)$high_floor_pass)
})

testthat::test_that("2N/4N disagreement and missing cohorts remain CatU", {
  source_joint_ploidy_utils()
  testthat::expect_equal(o2jpc_seed_consensus(c("CatA", "CatB"))[["category"]], "CatU")
  testthat::expect_equal(o2jpc_seed_consensus("CatA")[["category"]], "CatU")
  testthat::expect_equal(o2jpc_seed_consensus(c("CatC", "CatC"))[["category"]], "CatC")
  testthat::expect_equal(o2jpc_seed_consensus(c(`2N` = "CatB", `4N` = "CatC"))[["category"]], "CatC")
})

testthat::test_that("transition diagnostics are deterministic", {
  source_joint_ploidy_utils()
  curves <- synthetic_joint_ploidy_curves()
  first <- o2jpc_curve_features(curves$day, curves$CatC)
  second <- o2jpc_curve_features(curves$day, curves$CatC)
  testthat::expect_identical(first, second)
})

testthat::test_that("category definition records all prespecified primary thresholds", {
  source_joint_ploidy_utils()
  definition <- o2jpc_category_definition()
  testthat::expect_true(all(c("high_tolerance", "low_endpoint", "plateau_min_days", "two_transition_bic_delta_cutoff") %in% definition$setting))
})
