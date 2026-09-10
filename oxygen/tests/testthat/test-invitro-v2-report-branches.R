testthat::test_that("v2 visualization recognizes and isolates O1/O2 branches", {
  viz_script <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "vis",
    "viz_invitro_model_O2_supply_demand_MAP_results.R"
  )
  expression_list <- parse(viz_script)
  env <- new.env(parent = globalenv())
  wanted <- c(
    "derive_invitro_lineage_label",
    "derive_invitro_lineage_passage_index",
    "ensure_invitro_plot_columns",
    "invitro_branch_from_key",
    "invitro_v2_branches",
    "subset_invitro_v2_branch",
    "subset_invitro_v2_control"
  )
  for (expr in expression_list) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && as.character(expr[[2]]) %in% wanted) {
      eval(expr, envir = env)
    }
  }

  keys <- c("20.5", "20.5_2", "20.5_2::O1", "20.5_2::O2")
  testthat::expect_identical(
    env$derive_invitro_lineage_label(keys),
    c("control", "deprived", "O1", "O2")
  )
  testthat::expect_identical(
    env$derive_invitro_lineage_passage_index(keys),
    c(1, 2, 2, 2)
  )

  fixture <- data.frame(
    segment_id = c(
      "20.5",
      "20.5_20.5",
      "20.5_20.5_20.5",
      "20.5_20.5_2::O1",
      "20.5_20.5_2::O2"
    ),
    lineage_label = c("control", "control", "control", "deprived", "deprived"),
    cohort = "2N",
    passage_index = c(1, 2, 3, 3, 3),
    stringsAsFactors = FALSE
  )
  testthat::expect_identical(env$invitro_v2_branches(fixture), c("O1", "O2"))
  o1 <- env$subset_invitro_v2_branch(fixture, "O1")
  testthat::expect_setequal(
    o1$segment_id,
    c("20.5", "20.5_20.5", "20.5_20.5_2::O1")
  )
  testthat::expect_true(all(o1$lineage_label == "O1"))
  testthat::expect_false(any(grepl("::O2$", o1$segment_id)))
  testthat::expect_false("20.5_20.5_20.5" %in% o1$segment_id)
  control <- env$subset_invitro_v2_control(fixture)
  testthat::expect_setequal(
    control$segment_id,
    c("20.5", "20.5_20.5", "20.5_20.5_20.5")
  )
  testthat::expect_true(all(control$lineage_label == "control"))
})

testthat::test_that("v2 report combines only daily trajectories and pairs other branch figures", {
  report_script <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "report",
    "render_invitro_fit_report.R"
  )
  expression_list <- parse(report_script)
  env <- new.env(parent = globalenv())
  wanted <- c(
    "make_figure_spec",
    "optional_figure",
    "branch_figure_pair_or_legacy",
    "branch_figure_triplet_or_legacy",
    "build_invitro_section_specs"
  )
  for (expr in expression_list) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && as.character(expr[[2]]) %in% wanted) {
      eval(expr, envir = env)
    }
  }

  viz_dir <- tempfile("v2-report-assets-")
  dir.create(viz_dir)
  branch_bases <- c(
    "invitro_o2_selected_live_panels",
    "invitro_missegregation_probability_over_passage",
    "invitro_growth_ploidy_burden_composite",
    "invitro_flow_density"
  )
  branch_assets <- unlist(lapply(branch_bases, function(x) paste0(x, "_", c("O1", "O2"))))
  file.create(file.path(viz_dir, paste0(branch_assets, ".png")))
  distribution_assets <- paste0(
    "invitro_distribution_heatmap_",
    c("control", "O1", "O2")
  )
  file.create(file.path(viz_dir, paste0(distribution_assets, ".png")))
  file.create(file.path(viz_dir, "invitro_daily_counts.png"))

  figures <- unlist(lapply(env$build_invitro_section_specs(viz_dir), `[[`, "figures"), recursive = FALSE)
  basenames <- vapply(figures, `[[`, character(1), "basename")
  testthat::expect_true("invitro_daily_counts" %in% basenames)
  testthat::expect_false(any(paste0("invitro_daily_counts_", c("O1", "O2")) %in% basenames))
  testthat::expect_true(all(branch_assets %in% basenames))
  testthat::expect_true(all(distribution_assets %in% basenames))
})
