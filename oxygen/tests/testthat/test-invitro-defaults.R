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
  testthat::expect_false("passage_mode" %in% names(cfg))
  testthat::expect_identical(env$INVITRO_PASSAGE_IMPLEMENTATION, "v2")
  testthat::expect_error(
    env$ivt_reject_removed_passage_mode("org", source = "test configuration"),
    "passage_mode has been removed from test configuration"
  )
})
