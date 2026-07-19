parameter_stage3_root <- file.path(repo_info$root, "oxygen", "code", "O2_supply_demand_MAP")

parameter_stage3_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

testthat::test_that("parameter-landscape workflow has separate simulation, analysis, vis, report, and runner entrypoints", {
  paths <- c(
    file.path(parameter_stage3_root, "simulation", "parameter_landscape", "generate_parameter_landscape_simulation_tables.R"),
    file.path(parameter_stage3_root, "analysis", "parameter_landscape_clustering", "analyze_parameter_landscape.R"),
    file.path(parameter_stage3_root, "analysis", "parameter_landscape_clustering", "parameter_contribution_analysis.R"),
    file.path(parameter_stage3_root, "vis", "parameter_landscape", "visualize_parameter_landscape.R"),
    file.path(parameter_stage3_root, "report", "parameter_landscape", "clustering_report.R"),
    file.path(parameter_stage3_root, "runner", "parameter_landscape", "run_parameter_landscape.R")
  )
  testthat::expect_true(all(file.exists(paths)))
  testthat::expect_silent(invisible(lapply(paths, parse)))
  testthat::expect_false(grepl("ggplot|png\\(|pdf\\(", parameter_stage3_text(paths[[1L]]), perl = TRUE))
  testthat::expect_false(grepl("simulation/|ggplot|png\\(|pdf\\(", parameter_stage3_text(paths[[2L]]), perl = TRUE))
  testthat::expect_false(grepl("best_params|fit_result|readRDS", parameter_stage3_text(paths[[4L]]), perl = TRUE))
  testthat::expect_false(grepl("ggplot|readRDS|best_params", parameter_stage3_text(paths[[5L]]), perl = TRUE))
})

testthat::test_that("old parameter-landscape entrypoints are short compatibility wrappers", {
  old_dir <- file.path(parameter_stage3_root, "analysis", "best_fit_parameter_feature", "02_parameter_landscape_clustering")
  wrappers <- list.files(old_dir, pattern = "[.]R$", full.names = TRUE)
  testthat::expect_true(length(wrappers) >= 10L)
  for (path in wrappers) {
    testthat::expect_true(length(readLines(path, warn = FALSE)) < 80L, info = path)
    testthat::expect_match(tolower(parameter_stage3_text(path)), "deprecated", fixed = TRUE, info = path)
  }
})

testthat::test_that("parameter-landscape canonical runner resolves all layers from any cwd", {
  runner <- file.path(parameter_stage3_root, "runner", "parameter_landscape", "run_parameter_landscape.R")
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(tempdir())
  output <- system2(file.path(R.home("bin"), "Rscript"), c(runner, "--dry_run=true", "--run_parts=all"), stdout = TRUE, stderr = TRUE)
  testthat::expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  testthat::expect_true(any(grepl("simulation/parameter_landscape", output, fixed = TRUE)))
  testthat::expect_true(any(grepl("analysis/parameter_landscape_clustering", output, fixed = TRUE)))
  testthat::expect_true(any(grepl("vis/parameter_landscape", output, fixed = TRUE)))
  testthat::expect_true(any(grepl("report/parameter_landscape", output, fixed = TRUE)))
})
