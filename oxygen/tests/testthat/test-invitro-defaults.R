testthat::test_that("in vitro default death mode is ploidy-related", {
  env <- new.env(parent = globalenv())
  source(
    file.path(
      repo_info$root,
      "oxygen",
      "code",
      "O2_supply_demand_MAP",
      "util",
      "o2_supply_demand_map_invitro_utils.R"
    ),
    local = env,
    chdir = TRUE
  )

  cfg <- env$ivt_build_default_cfg(
    repo_root = file.path(repo_info$root, "oxygen")
  )

  testthat::expect_identical(cfg$ploidy_O2_death, "ploidy_related")
  testthat::expect_identical(cfg$passage_mode, "org")

  v1_cfg <- env$ivt_build_default_cfg(
    repo_root = file.path(repo_info$root, "oxygen"),
    passage_mode = "V1"
  )
  testthat::expect_identical(v1_cfg$passage_mode, "v1")
  v2_cfg <- env$ivt_build_default_cfg(
    repo_root = file.path(repo_info$root, "oxygen"),
    passage_mode = "V2"
  )
  testthat::expect_identical(v2_cfg$passage_mode, "v2")
  testthat::expect_error(
    env$ivt_build_default_cfg(
      repo_root = file.path(repo_info$root, "oxygen"),
      passage_mode = "future"
    ),
    "passage_mode must be one of: org, v1, v2"
  )
})
