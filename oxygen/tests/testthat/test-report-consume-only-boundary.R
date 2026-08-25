report_boundary_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP",
  "report"
)

testthat::test_that("report code contains no fitting, simulation, or plotting producer", {
  report_files <- list.files(
    report_boundary_root,
    pattern = "[.](R|Rmd)$",
    full.names = TRUE,
    recursive = TRUE
  )
  testthat::expect_true(length(report_files) >= 5L)
  text <- paste(
    unlist(lapply(report_files, readLines, warn = FALSE)),
    collapse = "\n"
  )
  forbidden <- c(
    "ggplot2::",
    "ggsave\\(",
    "readRDS\\(",
    "fit_result[.]rds",
    "model_O2_supply_demand_MAP",
    "sourceCpp\\(",
    "source\\([^\\n]*(model|optimizer|simulation)",
    "sys[.]source\\([^\\n]*(model|optimizer|simulation)",
    "generate_invivo_invitro_comparison_figures",
    "save_comparison_plot_pair",
    "plot_joint_parameter_ratio_figure"
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("report layer contains forbidden producer dependency:", pattern)
    )
  }
})

testthat::test_that("general fit report consumes materialized visualization scopes", {
  path <- file.path(report_boundary_root, "render_fit_report.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expected_materialized_dirs <- c(
    'file.path\\(fit_dir, "viz"\\)',
    'file.path\\(viz_root, "invivo"\\)',
    'file.path\\(fit_dir, "viz", "invitro"\\)',
    'file.path\\(fit_dir, "viz", "invivo_vs_invitro"\\)',
    'file.path\\(fit_dir, "viz", "joint_parameters"\\)'
  )
  for (pattern in expected_materialized_dirs) {
    testthat::expect_match(text, pattern, perl = TRUE)
  }
  testthat::expect_match(text, "build_invivo_invitro_comparison_section")
  testthat::expect_match(text, "build_parameter_figure_specs")
  testthat::expect_false(grepl("dir[.]create\\(viz_dir", text, perl = TRUE))
})

testthat::test_that("general fit report isolates per-seed render intermediates", {
  path <- file.path(report_boundary_root, "render_fit_report.R")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  testthat::expect_match(
    text,
    'paste0\\("\\.", report_basename, "_intermediates"\\)',
    perl = TRUE
  )
  testthat::expect_equal(
    length(gregexpr("intermediates_dir = intermediates_dir", text, fixed = TRUE)[[1L]]),
    2L
  )
})
