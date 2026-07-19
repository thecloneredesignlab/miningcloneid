fixed_o2_stage_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

fixed_o2_read <- function(...) {
  paste(
    readLines(file.path(fixed_o2_stage_root, ...), warn = FALSE),
    collapse = "\n"
  )
}

testthat::test_that("fixed-O2 simulation directory is numerical-producer only", {
  root <- file.path(fixed_o2_stage_root, "simulation", "o2", "fixed_o2")
  paths <- list.files(root, pattern = "[.]R$", full.names = TRUE, recursive = TRUE)
  testthat::expect_setequal(
    basename(paths),
    c(
      "run_fixed_o2_simulation.R",
      "fixed_o2_simulation_utils.R",
      "fixed_o2_numerical_producers.R",
      "build_fixo2_eigen_attractor_features.R"
    )
  )
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  forbidden <- c(
    "ggplot2::",
    "ggsave\\s*\\(",
    "grDevices::pdf\\s*\\(",
    "grDevices::png\\s*\\(",
    "fixed_o2_plot_utils",
    "fixed_o2_report_utils",
    "fixed_o2_legacy_pipeline",
    "[/]vis[/]",
    "[/]report[/]"
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("fixed-O2 simulation contains cross-stage behavior:", pattern)
    )
  }
})

testthat::test_that("fixed-O2 staged analysis consumes tables only", {
  text <- paste(
    fixed_o2_read("analysis", "fixed_o2", "run_fixed_o2_analysis.R"),
    fixed_o2_read("analysis", "fixed_o2", "fixed_o2_analysis_utils.R"),
    sep = "\n"
  )
  forbidden <- c(
    "readRDS\\s*\\(",
    "o2ipa_source_model\\s*\\(",
    "model_O2_supply_demand_MAP",
    "sourceCpp\\s*\\(",
    "best_params",
    "ggplot2::",
    "ggsave\\s*\\(",
    "grDevices::pdf\\s*\\(",
    "write_html_report\\s*\\("
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("fixed-O2 analysis contains producer/render behavior:", pattern)
    )
  }
  testthat::expect_match(text, "fixed_o2_attractors_by_seed.tsv", fixed = TRUE)
  testthat::expect_match(
    text,
    "fixed_o2_counterfactual_summary_by_seed.tsv",
    fixed = TRUE
  )
  testthat::expect_match(text, "agreement_augmented_data.tsv", fixed = TRUE)
})

testthat::test_that("fixed-O2 visualization has no model or producer dependency", {
  text <- paste(
    fixed_o2_read("vis", "fixed_o2", "render_fixed_o2_figures.R"),
    fixed_o2_read("vis", "fixed_o2", "fixed_o2_plot_utils.R"),
    sep = "\n"
  )
  forbidden <- c(
    "readRDS\\s*\\(",
    "o2ipa_source_model\\s*\\(",
    "model_O2_supply_demand_MAP",
    "sourceCpp\\s*\\(",
    "best_params",
    "generate_missing_simulation",
    "run_fixed_o2_simulation[.]R",
    "fixed_o2_analysis_utils[.]R",
    "sys[.]source\\([^\\n]*simulation",
    "source\\([^\\n]*simulation"
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("fixed-O2 visualization contains producer dependency:", pattern)
    )
  }
  testthat::expect_match(text, "fixed_o2_figure_manifest.tsv", fixed = TRUE)
})

testthat::test_that("fixed-O2 report entry is assembly only", {
  text <- paste(
    fixed_o2_read("report", "fixed_o2", "render_fixed_o2_report.R"),
    fixed_o2_read("report", "fixed_o2", "fixed_o2_report_utils.R"),
    sep = "\n"
  )
  forbidden <- c(
    "readRDS\\s*\\(",
    "o2ipa_source_model\\s*\\(",
    "model_O2_supply_demand_MAP",
    "sourceCpp\\s*\\(",
    "ggplot2::",
    "ggsave\\s*\\(",
    "run_fixed_o2_simulation[.]R"
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("fixed-O2 report contains non-assembly behavior:", pattern)
    )
  }
  testthat::expect_match(text, "write_html_report", fixed = TRUE)
})

testthat::test_that("historical fixed-O2 CLIs are thin compatibility wrappers", {
  paths <- c(
    file.path(fixed_o2_stage_root, "analysis", "fixed_o2", "FixO2_invivo.R"),
    file.path(
      fixed_o2_stage_root,
      "analysis",
      "best_fit_parameter_feature",
      "01_fixed_o2",
      "FixO2_invivo.R"
    )
  )
  line_counts <- vapply(paths, function(path) length(readLines(path, warn = FALSE)), integer(1))
  testthat::expect_true(all(line_counts < 100L))
  text <- vapply(paths, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1))
  testthat::expect_true(all(grepl(
    "deprecated compatibility orchestrator",
    tolower(text),
    fixed = TRUE
  )))
  testthat::expect_true(all(grepl('"runner"', text, fixed = TRUE)))
  testthat::expect_match(text[[1L]], "FixO2_invivo_500seed", fixed = TRUE)
  testthat::expect_match(
    text[[2L]],
    "best_fit_parameter_feature",
    fixed = TRUE
  )

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tempdir())
  for (path in paths) {
    env <- new.env(parent = globalenv())
    source(path, local = env, chdir = FALSE)
    testthat::expect_true(exists("fixo2_main", envir = env, inherits = FALSE))
    testthat::expect_identical(names(formals(env$fixo2_main)), "args")
    testthat::expect_false(
      exists("FIXO2_SIMULATION_ENV", envir = env, inherits = FALSE)
    )
  }
  rscript <- file.path(R.home("bin"), "Rscript")
  for (path in paths) {
    output <- suppressWarnings(system2(
      rscript,
      c(path, "--run_parts", "invalid"),
      stdout = TRUE,
      stderr = TRUE
    ))
    testthat::expect_false(is.null(attr(output, "status")))
    testthat::expect_match(
      paste(output, collapse = "\n"),
      "Unknown run_parts value(s): invalid",
      fixed = TRUE
    )
  }
})

testthat::test_that("staged fixed-O2 entrypoints parse and source from any cwd", {
  paths <- c(
    file.path(fixed_o2_stage_root, "analysis", "fixed_o2", "run_fixed_o2_analysis.R"),
    file.path(fixed_o2_stage_root, "vis", "fixed_o2", "render_fixed_o2_figures.R"),
    file.path(fixed_o2_stage_root, "report", "fixed_o2", "render_fixed_o2_report.R")
  )
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tempdir())
  for (path in paths) {
    testthat::expect_silent(parse(path))
    env <- new.env(parent = globalenv())
    testthat::expect_silent(source(path, local = env, chdir = FALSE))
  }
  rscript <- file.path(R.home("bin"), "Rscript")
  for (path in paths) {
    output <- suppressWarnings(system2(
      rscript,
      c(path, "--help"),
      stdout = TRUE,
      stderr = TRUE
    ))
    testthat::expect_null(attr(output, "status"))
    testthat::expect_match(paste(output, collapse = "\n"), "Usage:", fixed = TRUE)
  }
})

testthat::test_that("fixed-O2 analysis mode assignment preserves the legacy contract", {
  env <- new.env(parent = globalenv())
  sys.source(
    file.path(fixed_o2_stage_root, "analysis", "fixed_o2", "run_fixed_o2_analysis.R"),
    envir = env,
    chdir = TRUE
  )
  tab <- data.frame(
    seed_id = rep(c("seed1", "seed2"), each = 2L),
    O2_pct = rep(c(0, 2), times = 2L),
    dominant_mean_ploidy = c(2.2, 2.1, 1.7, 1.6),
    stringsAsFactors = FALSE
  )
  modes <- env$fixo2_assign_attractor_modes(tab)
  testthat::expect_identical(
    unique(modes$trajectory_regime[modes$seed_id == "seed1"]),
    "mode1_attractor_dominant_ploidy_ge_2"
  )
  testthat::expect_identical(
    unique(modes$trajectory_regime[modes$seed_id == "seed2"]),
    "mode2_attractor_dominant_ploidy_lt_2"
  )
  reference <- env$fixo2_reference_mode_table(
    env$fixo2_attractor_mode_table(modes),
    2
  )
  testthat::expect_identical(reference$seed_id, c("seed1", "seed2"))
  testthat::expect_identical(reference$mode_reference_o2_pct, c(2, 2))
})

testthat::test_that("generated fixed-O2 modules contain real header newlines", {
  paths <- c(
    file.path(fixed_o2_stage_root, "simulation", "o2", "fixed_o2", "fixed_o2_numerical_producers.R"),
    file.path(fixed_o2_stage_root, "analysis", "fixed_o2", "fixed_o2_analysis_utils.R"),
    file.path(fixed_o2_stage_root, "vis", "fixed_o2", "fixed_o2_plot_utils.R"),
    file.path(fixed_o2_stage_root, "report", "fixed_o2", "fixed_o2_report_utils.R"),
    file.path(fixed_o2_stage_root, "runner", "fixed_o2", "fixed_o2_legacy_pipeline_functions.R")
  )
  for (path in paths) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl("\\\\n#", text, perl = TRUE), info = path)
  }
})
