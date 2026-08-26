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

testthat::test_that("joint report consumes v2 O1/O2/control assets and relationship figures", {
  report_script <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "report",
    "render_fit_report.R"
  )
  expression_list <- parse(report_script)
  env <- new.env(parent = globalenv())
  env$`%||%` <- function(x, y) if (is.null(x)) y else x
  wanted <- c(
    "fresh_png_companion_for_pdf",
    "make_figure_spec",
    "optional_figure",
    "optional_figure_with_layout",
    "joint_branch_figure_pair_or_legacy",
    "joint_branch_figure_triplet_or_legacy",
    "sequential_layout_groups_from_figures",
    "build_invitro_report_section_specs_for_joint"
  )
  for (expr in expression_list) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && as.character(expr[[2]]) %in% wanted) {
      eval(expr, envir = env)
    }
  }

  fit_dir <- tempfile("joint-v2-report-assets-")
  viz_dir <- file.path(fit_dir, "viz", "invitro")
  dir.create(viz_dir, recursive = TRUE)
  branch_assets <- unlist(lapply(
    c(
      "invitro_o2_selected_live_panels",
      "invitro_missegregation_probability_over_passage",
      "invitro_growth_ploidy_burden_composite",
      "invitro_flow_density"
    ),
    function(x) paste0(x, "_", c("O1", "O2"))
  ))
  distribution_assets <- paste0(
    "invitro_distribution_heatmap_",
    c("control", "O1", "O2")
  )
  relationship_assets <- c(
    "ms_rate_vs_nonviable_daughter_fraction",
    "death_rate_vs_missegregation_rate",
    "ploidy_vs_viability_after_ms"
  )
  all_assets <- c(branch_assets, distribution_assets, relationship_assets)
  file.create(file.path(viz_dir, paste0(all_assets, ".pdf")))

  sections <- env$build_invitro_report_section_specs_for_joint(fit_dir)
  figures <- unlist(lapply(sections, `[[`, "figures"), recursive = FALSE)
  filenames <- basename(vapply(figures, `[[`, character(1), "src"))

  testthat::expect_true(all(paste0(branch_assets, ".pdf") %in% filenames))
  testthat::expect_true(all(paste0(distribution_assets, ".pdf") %in% filenames))
  testthat::expect_true(all(paste0(relationship_assets, ".pdf") %in% filenames))
})

testthat::test_that("in-vitro diagnostics build standalone missegregation relationship plots", {
  testthat::skip_if_not_installed("ggplot2")
  plot_utils <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "vis",
    "invitro",
    "o2_supply_demand_map_invitro_diagnostic_plot_utils.R"
  )
  env <- new.env(parent = globalenv())
  sys.source(plot_utils, envir = env)

  oxygen <- rep(seq(0.1, 2, length.out = 5), 2)
  reference <- data.frame(
    oxygen_pct = oxygen,
    cohort = rep(c("2N", "4N"), each = 5),
    death_rate = seq(0.01, 0.1, length.out = 10),
    ms_rate = seq(0.02, 0.2, length.out = 10)
  )
  multi <- data.frame(
    cohort = rep(c("1.5x", "2N", "4N"), each = 5),
    ms_rate = seq(0.02, 0.3, length.out = 15),
    misseg_nonviable_daughter_fraction = seq(0.01, 0.4, length.out = 15)
  )
  viability <- data.frame(
    endpoint_value = seq(20, 100, length.out = 10),
    viability_after_ms = seq(0.7, 0.95, length.out = 10)
  )

  plots <- env$ivt_plot_missegregation_linked_relationships(reference, multi, viability)

  testthat::expect_setequal(
    names(plots),
    c(
      "ms_rate_vs_nonviable_daughter_fraction",
      "death_rate_vs_missegregation_rate",
      "ploidy_vs_viability_after_ms"
    )
  )
  testthat::expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
})
