fit_results_root <- file.path(repo_info$root, "oxygen", "code", "O2_supply_demand_MAP")

fit_results_files <- list(
  simulation = file.path(fit_results_root, "simulation", "fit_results", c(
    "materialize_extra_results_predictions.R", "materialize_invivo_long_ploidy_metrics.R"
  )),
  analysis = file.path(fit_results_root, "analysis", "fit_results", c(
    "analyze_extra_results.R", "extra_results_analysis_utils.R",
    "analyze_joint_sigma_soft_coupled_paired_seeds.R", "analyze_sigma_burden_extra_results.R",
    "select_invivo_best_long_ploidy_seed_from_metrics.R", "select_best_seed_from_summary.R"
  )),
  vis = file.path(fit_results_root, "vis", "fit_results", c(
    "plot_extra_results.R", "plot_extra_results_objective_violin.R",
    "plot_joint_sigma_soft_coupled_paired_seeds.R", "plot_sigma_burden_extra_results.R"
  )),
  report = file.path(fit_results_root, "report", "fit_results", c(
    "render_extra_results_report.R", "render_joint_sigma_soft_coupled_paired_seeds_report.R"
  )),
  runner = file.path(fit_results_root, "runner", "fit_results", c(
    "run_extra_results.R", "run_joint_sigma_soft_coupled_paired_seeds.R",
    "run_sigma_burden_extra_results.R", "select_invivo_best_long_ploidy_seed.R"
  )),
  compatibility = file.path(fit_results_root, "analysis", "fit_results", c(
    "extra_results.R", "extra_results_report.R", "plot_extra_results_objective_violin.R",
    "compare_joint_sigma_soft_coupled_paired_seeds.R", "compare_sigma_burden_extra_results.R",
    "select_invivo_best_long_ploidy_seed.R"
  )),
  util = file.path(fit_results_root, "util", "o2_supply_demand_map_fit_results_utils.R")
)

read_fit_results_code <- function(paths) {
  paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
}

testthat::test_that("fit-results canonical R files exist and parse", {
  paths <- unique(unlist(fit_results_files, use.names = FALSE))
  testthat::expect_true(all(file.exists(paths)))
  for (path in paths) testthat::expect_silent(parse(path))
  python_paths <- file.path(fit_results_root, c(
    "analysis/fit_results/collect_best_separate_fit_reports.py",
    "report/fit_results/collect_best_separate_fit_reports.py"
  ))
  testthat::expect_true(all(file.exists(python_paths)))
  for (path in python_paths) {
    status <- system2("python3", c(
      "-c",
      shQuote("import pathlib, sys; p=sys.argv[1]; compile(pathlib.Path(p).read_text(), p, 'exec')"),
      shQuote(path)
    ))
    testthat::expect_identical(status, 0L, info = path)
  }
})

testthat::test_that("fit-results simulation materializes numeric contracts only", {
  text <- read_fit_results_code(fit_results_files$simulation)
  forbidden <- c(
    "ggplot2::", "grDevices::pdf", "grDevices::png", "graphics::plot",
    "wilcox.test", "cor.test", "analysis_manifest.tsv", "visualization_manifest.tsv",
    "report_manifest.tsv", '"analysis"', '"vis"', '"report"'
  )
  for (pattern in forbidden) testthat::expect_false(grepl(pattern, text, fixed = TRUE), info = pattern)
  testthat::expect_match(text, "simulation_manifest.tsv", fixed = TRUE)
  testthat::expect_match(text, "predict1000_ploidy_seed_day_mean.tsv", fixed = TRUE)
  testthat::expect_match(text, "invivo_long_ploidy_metrics.tsv", fixed = TRUE)
})

testthat::test_that("fit-results analysis is table and statistics only", {
  text <- read_fit_results_code(fit_results_files$analysis)
  patterns <- c(
    ggplot = "(^|[^A-Za-z0-9_.])(ggplot|ggsave)\\s*\\(",
    ggplot_namespace = "ggplot2::",
    base_plot = "(^|[^A-Za-z0-9_.])(plot|lines|points|hist|heatmap|image|boxplot|legend|axis|abline|polygon|segments|matplot)\\s*\\(",
    device = "(grDevices::)?(pdf|png|jpeg|svg)\\s*\\(",
    report = "(rmarkdown|quarto)::render\\s*\\(",
    source_simulation = "(source|sys.source)\\s*\\([^\\n]*simulation",
    invoke_simulation = "(system2|system)\\s*\\([^\\n]*(simulate_|materialize_.*prediction|generate_.*simulation)"
  )
  for (name in names(patterns)) testthat::expect_false(grepl(patterns[[name]], text, perl = TRUE, ignore.case = TRUE), info = name)
  testthat::expect_match(text, "analysis_manifest.tsv", fixed = TRUE)
})

testthat::test_that("fit-results visualization consumes materialized tables only", {
  text <- read_fit_results_code(fit_results_files$vis)
  forbidden <- c(
    "best_params.tsv", "fit_result.rds", "readRDS(", "fit_summary.tsv",
    "model_O2_supply_demand_MAP", '"optimizer"', '"simulation", "fit_results"',
    '"analysis", "fit_results"'
  )
  for (pattern in forbidden) testthat::expect_false(grepl(pattern, text, fixed = TRUE), info = pattern)
  testthat::expect_match(text, "visualization_manifest.tsv", fixed = TRUE)
  testthat::expect_match(text, "analysis_manifest.tsv", fixed = TRUE)
})

testthat::test_that("fit-results reports only assemble upstream artifacts", {
  text <- read_fit_results_code(fit_results_files$report)
  patterns <- c(
    ggplot = "(^|[^A-Za-z0-9_.])(ggplot|ggsave)\\s*\\(",
    ggplot_namespace = "ggplot2::",
    base_plot = "(^|[^A-Za-z0-9_.])(plot|lines|points|hist|heatmap|image|boxplot|legend|axis|abline|polygon|segments|matplot)\\s*\\(",
    device = "(grDevices::)?(pdf|png|jpeg|svg)\\s*\\(",
    read_fit = "readRDS\\s*\\(",
    source_business_layer = "(source|sys.source)\\s*\\([^\\n]*(simulation|analysis|vis)"
  )
  for (name in names(patterns)) testthat::expect_false(grepl(patterns[[name]], text, perl = TRUE, ignore.case = TRUE), info = name)
  testthat::expect_match(text, "report_manifest.tsv", fixed = TRUE)
})

testthat::test_that("fit-results runners preserve canonical stage order", {
  expected <- list(
    run_extra_results.R = c("materialize_extra_results_predictions.R", "analyze_extra_results.R", "plot_extra_results.R", "render_extra_results_report.R"),
    run_joint_sigma_soft_coupled_paired_seeds.R = c("analyze_joint_sigma_soft_coupled_paired_seeds.R", "plot_joint_sigma_soft_coupled_paired_seeds.R", "render_joint_sigma_soft_coupled_paired_seeds_report.R"),
    run_sigma_burden_extra_results.R = c("analyze_sigma_burden_extra_results.R", "plot_sigma_burden_extra_results.R"),
    select_invivo_best_long_ploidy_seed.R = c("materialize_invivo_long_ploidy_metrics.R", "select_invivo_best_long_ploidy_seed_from_metrics.R")
  )
  for (name in names(expected)) {
    path <- file.path(fit_results_root, "runner", "fit_results", name)
    text <- read_fit_results_code(path)
    positions <- vapply(expected[[name]], function(marker) regexpr(marker, text, fixed = TRUE)[[1]], integer(1))
    testthat::expect_true(all(positions > 0L), info = name)
    testthat::expect_true(all(diff(positions) > 0L), info = name)
  }
})

testthat::test_that("fit-results canonical defaults stay inside the active checkout", {
  canonical <- unique(c(
    unlist(fit_results_files[c("simulation", "analysis", "vis", "report", "runner")], use.names = FALSE),
    fit_results_files$util
  ))
  text <- read_fit_results_code(canonical)
  testthat::expect_false(grepl("/Users/", text, fixed = TRUE))
  testthat::expect_false(grepl("Constant_WGD", text, fixed = TRUE))
  testthat::expect_match(
    read_fit_results_code(file.path(fit_results_root, "runner", "fit_results", "run_sigma_burden_extra_results.R")),
    'file.path(WORKFLOW_ROOT, "..", "..", "results")',
    fixed = TRUE
  )
})

testthat::test_that("historical fit-results R paths are explicit compatibility wrappers", {
  for (path in fit_results_files$compatibility) {
    header <- paste(head(readLines(path, warn = FALSE), 20L), collapse = "\n")
    testthat::expect_true(grepl("DEPRECATED COMPATIBILITY", header, fixed = TRUE), info = path)
  }
})

testthat::test_that("long-ploidy simulation and selection preserve the numeric contract", {
  root <- tempfile("fit_results_long_ploidy_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  for (seed in 1:2) {
    seed_dir <- file.path(root, paste0("seed", seed))
    dir.create(file.path(seed_dir, "simulation", "invivo"), recursive = TRUE)
    utils::write.table(
      data.frame(metric = c("objective"), value = c(if (seed == 1) 2 else 1)),
      file.path(seed_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
    )
    value <- if (seed == 1) c(44, 50, 44, 51) else c(44, 42, 44, 43)
    utils::write.table(
      data.frame(cohort = rep(c("2N", "4N"), each = 2), day = rep(c(0, 1000), 2), weighted_mean_N = value),
      file.path(seed_dir, "simulation", "invivo", "predict_ploidy_weighted_mean_0_1000day.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  sim_dir <- file.path(root, "materialized")
  out_tsv <- file.path(root, "selection.tsv")
  out_dir <- file.path(root, "selection.dir")
  rscript <- file.path(R.home("bin"), "Rscript")
  sim_status <- system2(rscript, c(
    shQuote(fit_results_files$simulation[[2]]), shQuote(paste0("--invivo_dir=", root)),
    shQuote(paste0("--simulation_dir=", sim_dir))
  ))
  testthat::expect_identical(sim_status, 0L)
  metrics <- utils::read.delim(file.path(sim_dir, "invivo_long_ploidy_metrics.tsv"), check.names = FALSE)
  testthat::expect_equal(metrics$long_term_weighted_mean_N, c(50, 42))
  analysis_status <- system2(rscript, c(
    shQuote(file.path(fit_results_root, "analysis", "fit_results", "select_invivo_best_long_ploidy_seed_from_metrics.R")),
    shQuote(paste0("--simulation_dir=", sim_dir)), shQuote(paste0("--out_tsv=", out_tsv)),
    shQuote(paste0("--best_dir_file=", out_dir)), "--threshold_N=44"
  ))
  testthat::expect_identical(analysis_status, 0L)
  selected <- utils::read.delim(out_tsv, check.names = FALSE)
  testthat::expect_identical(selected$seed[selected$selected], 1L)
  testthat::expect_equal(selected$long_term_weighted_mean_N, c(50, 42))
  testthat::expect_match(readLines(out_dir), "seed1$", perl = TRUE)
})

testthat::test_that("fit-results consumers fail before creating output without manifests", {
  empty <- tempfile("fit_results_empty_")
  dir.create(empty)
  on.exit(unlink(empty, recursive = TRUE, force = TRUE), add = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  cases <- list(
    list(fit_results_files$analysis[[1]], c(paste0("--run_dir=", empty), paste0("--simulation_dir=", empty), "--analysis_dir=%OUT%")),
    list(fit_results_files$vis[[1]], c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), "--viz_dir=%OUT%")),
    list(fit_results_files$report[[1]], c(paste0("--simulation_dir=", empty), paste0("--analysis_dir=", empty), paste0("--viz_dir=", empty), "--report_dir=%OUT%"))
  )
  for (case in cases) {
    destination <- tempfile("fit_results_uncreated_")
    args <- sub("%OUT%", destination, case[[2]], fixed = TRUE)
    out <- suppressWarnings(system2(rscript, c(shQuote(case[[1]]), shQuote(args)), stdout = TRUE, stderr = TRUE))
    testthat::expect_false(is.null(attr(out, "status")), info = case[[1]])
    testthat::expect_false(dir.exists(destination), info = case[[1]])
  }
})
