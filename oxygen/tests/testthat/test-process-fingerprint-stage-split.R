process_stage_root <- file.path(
  repo_info$root, "oxygen", "code", "O2_supply_demand_MAP"
)

process_stage_files <- list(
  simulation = file.path(process_stage_root, "simulation", "process_fingerprints", "generate_process_fingerprint_outputs.R"),
  simulation_utils = file.path(process_stage_root, "simulation", "process_fingerprints", "process_fingerprint_simulation_utils.R"),
  event_simulation = file.path(process_stage_root, "simulation", "process_fingerprints", "generate_o2_ploidy_event_inputs.R"),
  process_analysis = file.path(process_stage_root, "analysis", "process_fingerprints", "analyze_invivo_process_fingerprints.R"),
  ploidy_analysis = file.path(process_stage_root, "analysis", "process_fingerprints", "analyze_invivo_ploidy_regimes.R"),
  event_analysis = file.path(process_stage_root, "analysis", "process_fingerprints", "analyze_invivo_o2_ploidy_event_coupling_tables.R"),
  medium_analysis = file.path(process_stage_root, "analysis", "process_fingerprints", "analyze_invivo_medium_o2_window_tables.R"),
  process_utils = file.path(process_stage_root, "analysis", "process_fingerprints", "process_fingerprint_utils.R"),
  ploidy_utils = file.path(process_stage_root, "analysis", "process_fingerprints", "ploidy_regime_utils.R"),
  process_vis = file.path(process_stage_root, "vis", "process_fingerprints", "plot_invivo_process_fingerprints.R"),
  ploidy_vis = file.path(process_stage_root, "vis", "process_fingerprints", "plot_invivo_ploidy_regimes.R"),
  event_vis = file.path(process_stage_root, "vis", "process_fingerprints", "plot_o2_ploidy_event_coupling.R"),
  medium_vis = file.path(process_stage_root, "vis", "process_fingerprints", "plot_medium_o2_windows.R"),
  vis_utils = file.path(process_stage_root, "vis", "process_fingerprints", "process_fingerprint_visualization_utils.R"),
  process_report = file.path(process_stage_root, "report", "process_fingerprints", "render_invivo_process_fingerprint_report.R"),
  ploidy_report = file.path(process_stage_root, "report", "process_fingerprints", "render_invivo_ploidy_regime_report.R"),
  event_report = file.path(process_stage_root, "report", "process_fingerprints", "render_o2_ploidy_event_coupling_report.R"),
  medium_report = file.path(process_stage_root, "report", "process_fingerprints", "render_medium_o2_window_report.R"),
  process_compat = file.path(process_stage_root, "analysis", "process_fingerprints", "run_invivo_process_analysis.R"),
  ploidy_compat = file.path(process_stage_root, "analysis", "process_fingerprints", "run_invivo_ploidy_regime_analysis.R"),
  event_compat = file.path(process_stage_root, "analysis", "process_fingerprints", "analyze_invivo_o2_ploidy_event_coupling.R"),
  medium_compat = file.path(process_stage_root, "analysis", "process_fingerprints", "analyze_invivo_medium_o2_windows.R")
)

read_process_stage_file <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

testthat::test_that("process-fingerprint stage entrypoints exist and parse", {
  testthat::expect_true(all(file.exists(unlist(process_stage_files))))
  for (path in unlist(process_stage_files)) testthat::expect_silent(parse(path))
})

testthat::test_that("simulation producer materializes numeric contracts only", {
  text <- paste(
    read_process_stage_file(process_stage_files$simulation),
    read_process_stage_file(process_stage_files$simulation_utils),
    sep = "\n"
  )
  forbidden <- c(
    "grDevices::pdf", "ggplot", "ggsave", "graphics::plot",
    "cluster_distance", "trajectory_regime", "analysis_manifest.tsv",
    "visualization_manifest.tsv", "report_manifest.tsv"
  )
  for (pattern in forbidden) {
    testthat::expect_false(grepl(pattern, text, fixed = TRUE), info = pattern)
  }
  testthat::expect_match(text, "simulation_manifest.tsv", fixed = TRUE)
  testthat::expect_match(text, "simulation/invivo", fixed = TRUE)

  event_text <- read_process_stage_file(process_stage_files$event_simulation)
  for (pattern in c(
    "grDevices::pdf", "ggplot", "ggsave", "graphics::plot",
    "wilcox.test", "cor.test", "analysis_manifest.tsv",
    "visualization_manifest.tsv", "report_manifest.tsv"
  )) {
    testthat::expect_false(grepl(pattern, event_text, fixed = TRUE), info = pattern)
  }
  testthat::expect_match(event_text, "event_ploidy_timecourses.tsv", fixed = TRUE)
  testthat::expect_match(event_text, "event_o2_timecourses.tsv", fixed = TRUE)
  testthat::expect_match(event_text, "pfs_write_manifest", fixed = TRUE)
})

testthat::test_that("canonical analysis contains no plotting or reporting producer", {
  paths <- c(
    process_stage_files$process_analysis, process_stage_files$ploidy_analysis,
    process_stage_files$event_analysis, process_stage_files$medium_analysis,
    process_stage_files$process_utils, process_stage_files$ploidy_utils
  )
  text <- paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  forbidden <- c(
    "grDevices::pdf", "ggplot", "ggsave", "graphics::plot",
    "graphics::boxplot", "stats::heatmap", "analysis_summary.md"
  )
  for (pattern in forbidden) {
    testthat::expect_false(grepl(pattern, text, fixed = TRUE), info = pattern)
  }
  for (path in c(
    process_stage_files$process_analysis, process_stage_files$ploidy_analysis,
    process_stage_files$event_analysis, process_stage_files$medium_analysis
  )) {
    stage_text <- read_process_stage_file(path)
    testthat::expect_match(stage_text, "--simulation_dir", fixed = TRUE)
    testthat::expect_false(grepl("--run_dir", stage_text, fixed = TRUE), info = path)
    testthat::expect_false(grepl("source_model", stage_text, fixed = TRUE), info = path)
    testthat::expect_false(grepl("best_params", stage_text, fixed = TRUE), info = path)
    testthat::expect_false(grepl("readRDS", stage_text, fixed = TRUE), info = path)
    testthat::expect_match(stage_text, "analysis_manifest.tsv", fixed = TRUE)
  }
})

testthat::test_that("plot entities live in visualization and vis is a strict consumer", {
  analysis_text <- paste(
    read_process_stage_file(process_stage_files$process_utils),
    read_process_stage_file(process_stage_files$ploidy_utils),
    sep = "\n"
  )
  vis_text <- paste(
    read_process_stage_file(process_stage_files$process_vis),
    read_process_stage_file(process_stage_files$ploidy_vis),
    read_process_stage_file(process_stage_files$event_vis),
    read_process_stage_file(process_stage_files$medium_vis),
    read_process_stage_file(process_stage_files$vis_utils),
    sep = "\n"
  )
  testthat::expect_false(grepl("o2ipa_plot_outputs", analysis_text, fixed = TRUE))
  testthat::expect_false(grepl("o2pr_basic_figures", analysis_text, fixed = TRUE))
  testthat::expect_match(vis_text, "pfv_plot_process_outputs", fixed = TRUE)
  testthat::expect_match(vis_text, "pfv_plot_ploidy_regime_outputs", fixed = TRUE)
  testthat::expect_match(vis_text, "run_o2_ploidy_event_visualization", fixed = TRUE)
  testthat::expect_match(vis_text, "run_medium_o2_window_visualization", fixed = TRUE)
  testthat::expect_match(vis_text, "grDevices::pdf", fixed = TRUE)
  for (pattern in c("best_params", "fit_result.rds", "fit_config.rds", "source_model", "optimizer/", "model_O2")) {
    testthat::expect_false(grepl(pattern, vis_text, fixed = TRUE), info = pattern)
  }
  testthat::expect_match(vis_text, "visualization_manifest.tsv", fixed = TRUE)
})

testthat::test_that("report scripts consume manifests without analysis or plotting", {
  text <- paste(
    read_process_stage_file(process_stage_files$process_report),
    read_process_stage_file(process_stage_files$ploidy_report),
    read_process_stage_file(process_stage_files$event_report),
    read_process_stage_file(process_stage_files$medium_report),
    sep = "\n"
  )
  for (pattern in c("source_model", "best_params", "run_process_fingerprint_simulation", "grDevices::pdf", "ggplot", "hclust(")) {
    testthat::expect_false(grepl(pattern, text, fixed = TRUE), info = pattern)
  }
  testthat::expect_match(text, "simulation_manifest.tsv", fixed = TRUE)
  testthat::expect_match(text, "analysis_manifest.tsv", fixed = TRUE)
  testthat::expect_match(text, "visualization_manifest.tsv", fixed = TRUE)
  testthat::expect_match(text, "report_manifest.tsv", fixed = TRUE)
})

testthat::test_that("legacy process entrypoints explicitly orchestrate the four stages", {
  for (path in c(process_stage_files$process_compat, process_stage_files$ploidy_compat)) {
    text <- read_process_stage_file(path)
    positions <- vapply(
      c(
        "run_process_fingerprint_simulation(",
        if (identical(path, process_stage_files$process_compat)) "run_invivo_process_fingerprint_analysis(" else "run_invivo_ploidy_regime_analysis_stage(",
        if (identical(path, process_stage_files$process_compat)) "run_invivo_process_fingerprint_visualization(" else "run_invivo_ploidy_regime_visualization(",
        if (identical(path, process_stage_files$process_compat)) "run_invivo_process_fingerprint_report(" else "run_invivo_ploidy_regime_report("
      ),
      function(marker) regexpr(marker, text, fixed = TRUE)[[1]],
      integer(1)
    )
    testthat::expect_true(all(positions > 0L), info = path)
    testthat::expect_true(all(diff(positions) > 0L), info = path)
  }
  aliases <- file.path(
    process_stage_root, "analysis", "process_fingerprints",
    c(
      "build_invivo_process_fingerprints.R", "cluster_invivo_process_fingerprints.R",
      "build_ploidy_regime_fingerprints.R", "analyze_ploidy_regimes.R"
    )
  )
  testthat::expect_true(all(file.exists(aliases)))

  extra_compat <- list(
    event = list(
      path = process_stage_files$event_compat,
      markers = c(
        "run_o2_ploidy_event_input_simulation(",
        "run_o2_ploidy_event_coupling_analysis(",
        "run_o2_ploidy_event_visualization(",
        "run_o2_ploidy_event_report("
      )
    ),
    medium = list(
      path = process_stage_files$medium_compat,
      markers = c(
        "run_o2_ploidy_event_input_simulation(",
        "run_medium_o2_window_analysis(",
        "run_medium_o2_window_visualization(",
        "run_medium_o2_window_report("
      )
    )
  )
  for (entry in extra_compat) {
    text <- read_process_stage_file(entry$path)
    testthat::expect_match(text, "COMPATIBILITY ORCHESTRATOR", fixed = TRUE)
    positions <- vapply(entry$markers, function(marker) regexpr(marker, text, fixed = TRUE)[[1]], integer(1))
    testthat::expect_true(all(positions > 0L), info = entry$path)
    testthat::expect_true(all(diff(positions) > 0L), info = entry$path)
  }
})

testthat::test_that("consumers fail before creating outputs when upstream manifests are missing", {
  empty <- tempfile("missing_process_stage_")
  dir.create(empty)
  on.exit(unlink(empty, recursive = TRUE, force = TRUE), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  cases <- list(
    list(path = process_stage_files$process_analysis, args = c(paste0("--simulation_dir=", empty), "--analysis_dir=%OUT%")),
    list(path = process_stage_files$ploidy_analysis, args = c(paste0("--simulation_dir=", empty), "--analysis_dir=%OUT%")),
    list(path = process_stage_files$event_analysis, args = c(paste0("--simulation_dir=", empty), "--analysis_dir=%OUT%")),
    list(path = process_stage_files$medium_analysis, args = c(paste0("--simulation_dir=", empty), "--analysis_dir=%OUT%")),
    list(path = process_stage_files$process_vis, args = c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), "--viz_dir=%OUT%")),
    list(path = process_stage_files$ploidy_vis, args = c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), "--viz_dir=%OUT%")),
    list(path = process_stage_files$event_vis, args = c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), "--viz_dir=%OUT%")),
    list(path = process_stage_files$medium_vis, args = c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), "--viz_dir=%OUT%")),
    list(path = process_stage_files$process_report, args = c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), paste0("--viz_dir=", empty), "--report_dir=%OUT%")),
    list(path = process_stage_files$ploidy_report, args = c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), paste0("--viz_dir=", empty), "--report_dir=%OUT%")),
    list(path = process_stage_files$event_report, args = c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), paste0("--viz_dir=", empty), "--report_dir=%OUT%")),
    list(path = process_stage_files$medium_report, args = c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), paste0("--viz_dir=", empty), "--report_dir=%OUT%"))
  )
  for (case in cases) {
    destination <- tempfile("process_stage_output_")
    args <- sub("%OUT%", destination, case$args, fixed = TRUE)
    out <- suppressWarnings(system2(rscript, c(shQuote(case$path), shQuote(args)), stdout = TRUE, stderr = TRUE))
    testthat::expect_false(is.null(attr(out, "status")), info = case$path)
    testthat::expect_false(dir.exists(destination), info = case$path)
  }
})

testthat::test_that("process helpers create only directories owned by their layer", {
  helper_paths <- c(
    file.path(
      process_stage_root,
      "simulation",
      "process_fingerprints",
      "process_fingerprint_simulation_legacy_utils.R"
    ),
    file.path(
      process_stage_root,
      "simulation",
      "process_fingerprints",
      "ploidy_regime_simulation_utils.R"
    ),
    file.path(
      process_stage_root,
      "analysis",
      "process_fingerprints",
      "process_fingerprint_utils.R"
    ),
    file.path(
      process_stage_root,
      "analysis",
      "process_fingerprints",
      "ploidy_regime_utils.R"
    )
  )
  helper_names <- c(
    "o2ipa_mkdirs",
    "o2pr_mkdirs",
    "o2ipa_mkdirs",
    "o2pr_mkdirs"
  )
  for (i in seq_along(helper_paths)) {
    env <- new.env(parent = globalenv())
    source(helper_paths[[i]], local = env)
    out_dir <- tempfile("process_layer_dirs_")
    env[[helper_names[[i]]]](out_dir)
    testthat::expect_true(
      all(dir.exists(file.path(out_dir, c("tables", "cache", "logs")))),
      info = helper_paths[[i]]
    )
    testthat::expect_false(
      dir.exists(file.path(out_dir, "figures")),
      info = helper_paths[[i]]
    )
    testthat::expect_false(
      dir.exists(file.path(out_dir, "report")),
      info = helper_paths[[i]]
    )
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
})
