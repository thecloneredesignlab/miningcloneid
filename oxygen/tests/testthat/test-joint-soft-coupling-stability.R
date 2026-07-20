find_joint_coupling_repo_root <- function() {
  current <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:8) {
    target <- file.path(current, "oxygen", "code", "O2_supply_demand_MAP", "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R")
    if (file.exists(target)) return(current)
    parent <- dirname(current); if (identical(parent, current)) break; current <- parent
  }
  stop("Cannot locate repository root")
}

source_joint_coupling_utils <- function() {
  root <- find_joint_coupling_repo_root()
  source(file.path(root, "oxygen", "code", "O2_supply_demand_MAP", "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"), local = .GlobalEnv)
  root
}

testthat::test_that("ClassA/B/C uses inclusive ClassB boundaries at threshold 1.1", {
  source_joint_coupling_utils()
  values <- c(-1, 0, (1 / 1.1) - 1e-12, 1 / 1.1, 1, 1.1, 1.1 + 1e-12, NA, Inf)
  observed <- as.character(o2jca_ratio_class(values, 1.1))
  testthat::expect_identical(observed, c("Invalid", "Invalid", "ClassA", "ClassB", "ClassB", "ClassB", "ClassC", "Invalid", "Invalid"))
})

testthat::test_that("within-parameter stability distinguishes union, strict intersection, and graded cores", {
  source_joint_coupling_utils()
  data <- data.frame(
    ratio_vivo_to_vitro = c(rep(1, 95), rep(2, 5)),
    log2_ratio_vivo_to_vitro = log2(c(rep(1, 95), rep(2, 5))),
    ratio_class = c(rep("ClassB", 95), rep("ClassC", 5)),
    stringsAsFactors = FALSE
  )
  result <- o2jca_summarize_parameter_group(data, 1.1)
  testthat::expect_true(result$union_ClassB)
  testthat::expect_false(result$intersection_ClassB)
  testthat::expect_true(result$stable90_ClassB)
  testthat::expect_true(result$stable95_ClassB)
  testthat::expect_false(result$intersection_ClassC)
  testthat::expect_equal(result$prop_ClassB, 0.95)
  testthat::expect_equal(result$prop_ratio_gt1, 0.05)
})

testthat::test_that("audited biological process map covers exactly the 14 soft-coupled parameters", {
  source_joint_coupling_utils()
  map <- o2jca_process_map()
  testthat::expect_setequal(map$parameter, o2jca_parameter_levels())
  testthat::expect_equal(nrow(map), 14L)
  testthat::expect_equal(map$primary_process[map$parameter == "p_mis_base"], "CIN/missegregation")
  testthat::expect_equal(map$primary_process[map$parameter == "p_wgd"], "Whole-genome doubling")
  testthat::expect_true(all(map$mapping_status == "audited_from_parameter_table"))
})

testthat::test_that("the named seed1 report input reproduces the stored in-vivo/in-vitro ratio", {
  root <- source_joint_coupling_utils()
  result_root <- file.path(root, "oxygen", "results", "fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540")
  soft_path <- file.path(result_root, "fit_joint_tsne_vi_seed138_C03Sc01_vt_seed10", "extra_results", "joint_soft_coupling_all.tsv")
  testthat::skip_if_not(file.exists(soft_path), "Selective local result mirror is unavailable")
  data <- utils::read.delim(soft_path, check.names = FALSE, stringsAsFactors = FALSE)
  row <- data[data$seed == "seed1" & data$parameter == "O2_crit", , drop = FALSE]
  testthat::expect_equal(nrow(row), 1L)
  testthat::expect_equal(row$ratio_vivo_to_vitro, row$vivo_natural / row$vitro_natural, tolerance = 1e-12)
  testthat::expect_equal(as.character(o2jca_ratio_class(row$ratio_vivo_to_vitro, 1.1)), "ClassA")
})

testthat::test_that("duplicate pair-seed-parameter rows fail the key contract", {
  source_joint_coupling_utils()
  data <- data.frame(pair_id = c("p", "p"), seed = c("seed1", "seed1"), parameter = c("x", "x"))
  testthat::expect_error(o2jca_assert_unique_key(data, c("pair_id", "seed", "parameter"), "fixture"), "duplicate key")
})

testthat::test_that("small objective strata skip empty subsets without failing", {
  source_joint_coupling_utils()
  seeds <- data.frame(seed = paste0("seed", 1:3), objective = c(1, 2, 3))
  strata <- o2jca_objective_strata(seeds)
  testthat::expect_true(nrow(strata) > 0)
  testthat::expect_false("top_10pct" %in% strata$objective_stratum)
  testthat::expect_true(all(c("all_valid", "best_seed") %in% strata$objective_stratum))
})
