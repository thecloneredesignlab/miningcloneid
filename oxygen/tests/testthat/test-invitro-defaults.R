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
})
