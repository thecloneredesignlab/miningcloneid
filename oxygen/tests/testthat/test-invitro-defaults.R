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

testthat::test_that("in vitro parameter tables share passage-rate uncertainty", {
  table_names <- c(
    "parameter_table_invitro.csv",
    "parameter_table_invitro_buffering.csv",
    "parameter_table_invitro_wgd_bimodal.csv"
  )
  rows <- lapply(table_names, function(table_name) {
    path <- file.path(
      repo_info$root,
      "oxygen",
      "data",
      "O2_supply_demand",
      table_name
    )
    table <- utils::read.csv(path, stringsAsFactors = FALSE)
    row <- table[table$param_symbol == "sigma_growth", , drop = FALSE]
    testthat::expect_equal(nrow(row), 1L, info = table_name)
    row$parameter_table <- table_name
    row
  })
  rows <- do.call(rbind, rows)

  testthat::expect_equal(rows$init_value, rep(0.2, length(table_names)))
  testthat::expect_equal(rows$lower_bound, rep(0.05, length(table_names)))
  testthat::expect_equal(rows$upper_bound, rep(0.5, length(table_names)))
  testthat::expect_true(all(grepl("passage-average", rows$description, fixed = TRUE)))
  testthat::expect_true(all(grepl("day^-1", rows$description, fixed = TRUE)))
})
