stage2_workflow_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

testthat::test_that("in-vitro visualization is a pure materialized-table consumer", {
  path <- file.path(
    stage2_workflow_root,
    "vis",
    "viz_invitro_model_O2_supply_demand_MAP_results.R"
  )
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  forbidden <- c(
    "readRDS\\(",
    "fit_result[.]rds",
    "best_params",
    "model[/]model_O2_supply_demand_MAP",
    "viz_invivo_model_O2_supply_demand_MAP_results",
    "ivt_collect_daily_counts",
    "compute_invitro_missegregation_probability_timecourse"
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, perl = TRUE),
      info = paste("pure vis contains forbidden dependency:", pattern)
    )
  }
  testthat::expect_match(text, 'file.path\\(fit_dir, "simulation", "invitro"\\)')
  testthat::expect_match(text, 'file.path\\(fit_dir, "analysis", "invitro"\\)')
  testthat::expect_match(text, 'file.path\\(fit_dir, "viz", "invitro"\\)')
  testthat::expect_match(text, 'device = headless_png_device', fixed = TRUE)
  testthat::expect_match(text, 'capabilities("cairo")', fixed = TRUE)
})

testthat::test_that("in-vitro simulation and diagnostics producers do not plot", {
  producer_paths <- c(
    list.files(
      file.path(stage2_workflow_root, "simulation", "invitro"),
      pattern = "[.]R$",
      full.names = TRUE,
      recursive = TRUE
    ),
    list.files(
      file.path(stage2_workflow_root, "analysis", "fit_diagnostics"),
      pattern = "[.]R$",
      full.names = TRUE,
      recursive = TRUE
    )
  )
  testthat::expect_true(length(producer_paths) >= 3L)
  text <- paste(unlist(lapply(producer_paths, readLines, warn = FALSE)), collapse = "\n")
  forbidden <- c("ggplot", "ggsave", "[.]pdf", "[.]png", "[.]svg", "source.*[/]vis")
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("producer contains plotting behavior:", pattern)
    )
  }
})

testthat::test_that("stage-2 in-vitro entrypoints and output contracts are stable", {
  simulation_entry <- file.path(
    stage2_workflow_root,
    "simulation",
    "invitro",
    "run_invitro_simulation_outputs.R"
  )
  diagnostics_entry <- file.path(
    stage2_workflow_root,
    "analysis",
    "fit_diagnostics",
    "run_invitro_fit_diagnostics.R"
  )
  vis_entry <- file.path(
    stage2_workflow_root,
    "vis",
    "viz_invitro_model_O2_supply_demand_MAP_results.R"
  )
  testthat::expect_true(all(file.exists(c(simulation_entry, diagnostics_entry, vis_entry))))

  simulation_text <- paste(readLines(simulation_entry, warn = FALSE), collapse = "\n")
  diagnostics_text <- paste(readLines(diagnostics_entry, warn = FALSE), collapse = "\n")
  testthat::expect_match(simulation_text, 'file.path\\(fit_dir, "simulation", "invitro"\\)')
  testthat::expect_match(diagnostics_text, 'file.path\\(fit_dir, "simulation", "invitro"\\)')
  testthat::expect_match(diagnostics_text, 'file.path\\(fit_dir, "analysis", "invitro"\\)')
  testthat::expect_match(diagnostics_text, "invitro_optimizer_population.tsv", fixed = TRUE)
  for (pattern in c("readRDS\\(", "fit_result", "best_params")) {
    testthat::expect_false(
      grepl(pattern, diagnostics_text, perl = TRUE),
      info = paste("diagnostics contains forbidden fitting dependency:", pattern)
    )
  }
  testthat::expect_match(simulation_text, "invitro_simulation_schema.tsv", fixed = TRUE)
  testthat::expect_match(simulation_text, "invitro_simulation_manifest.tsv", fixed = TRUE)
  testthat::expect_match(diagnostics_text, "invitro_fit_diagnostics_schema.tsv", fixed = TRUE)
  testthat::expect_match(diagnostics_text, "invitro_fit_diagnostics_manifest.tsv", fixed = TRUE)
})

testthat::test_that("in-vitro simulation modules keep domain responsibilities separate", {
  simulation_root <- file.path(stage2_workflow_root, "simulation", "invitro")
  paths <- c(
    util = file.path(
      stage2_workflow_root,
      "util",
      "o2_supply_demand_map_invitro_postfit_response_utils.R"
    ),
    population = file.path(
      simulation_root,
      "population",
      "invitro_population_simulation_utils.R"
    ),
    ploidy = file.path(simulation_root, "ploidy", "invitro_ploidy_simulation_utils.R"),
    cin = file.path(simulation_root, "cin", "invitro_cin_simulation_utils.R"),
    o2 = file.path(simulation_root, "o2", "invitro_o2_functional_response_utils.R")
  )
  testthat::expect_true(all(file.exists(paths)))
  text <- lapply(paths, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"))
  testthat::expect_match(text$population, "invitro_daily_counts", fixed = TRUE)
  testthat::expect_match(text$ploidy, "functional_curve_ploidy", fixed = TRUE)
  testthat::expect_match(text$cin, "missegregation", fixed = TRUE)
  testthat::expect_match(text$o2, "functional_curve_oxygen", fixed = TRUE)
  testthat::expect_false(grepl("functional_curve_ploidy", text$o2, fixed = TRUE))
  testthat::expect_false(grepl("functional_curve_oxygen", text$ploidy, fixed = TRUE))
})
