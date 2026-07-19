stage2_postfit_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

read_stage2_script <- function(...) {
  paste(
    readLines(file.path(stage2_postfit_root, ...), warn = FALSE),
    collapse = "\n"
  )
}

testthat::test_that("joint comparison visualization only consumes simulation tables", {
  text <- read_stage2_script(
    "vis",
    "joint",
    "plot_invivo_invitro_comparison.R"
  )
  forbidden <- c(
    "best_params",
    "fit_result[.]rds",
    "fit_config[.]rds",
    "model[/]model_O2_supply_demand_MAP",
    "optimizer[/]",
    "cpp_o2",
    "system2\\([^)]*Rscript"
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("joint pure vis contains forbidden dependency:", pattern)
    )
  }
  for (table_name in c(
    "functional_curve_oxygen.tsv",
    "functional_curve_oxygen_multi_ploidy.tsv",
    "functional_curve_ploidy.tsv",
    "functional_curve_ploidy_by_o2.tsv"
  )) {
    testthat::expect_match(text, table_name, fixed = TRUE)
  }
  testthat::expect_match(text, "visualization_manifest.tsv", fixed = TRUE)
})

testthat::test_that("joint parameter visualization only consumes analysis tables", {
  text <- read_stage2_script(
    "vis",
    "joint",
    "plot_joint_parameter_ratios.R"
  )
  forbidden <- c(
    "best_params",
    "fit_result[.]rds",
    "fit_config[.]rds",
    "joint_soft_coupling[.]tsv",
    "model[/]",
    "optimizer[/]"
  )
  for (pattern in forbidden) {
    testthat::expect_false(
      grepl(pattern, text, ignore.case = TRUE, perl = TRUE),
      info = paste("joint parameter vis contains forbidden dependency:", pattern)
    )
  }
  testthat::expect_match(text, "joint_parameter_ratio.tsv", fixed = TRUE)
  testthat::expect_match(
    text,
    "joint_parameter_visualization_manifest.tsv",
    fixed = TRUE
  )
})

testthat::test_that("fit report consumes pre-existing joint figures", {
  text <- read_stage2_script("report", "render_fit_report.R")
  build_start <- regexpr(
    "build_invivo_invitro_comparison_section <- function",
    text,
    fixed = TRUE
  )
  testthat::expect_gt(as.integer(build_start), 0L)
  build_text <- substr(text, build_start, nchar(text))
  build_text <- substr(
    build_text,
    1L,
    regexpr("\nbuild_", build_text, fixed = TRUE)[1L] - 1L
  )
  testthat::expect_false(
    grepl("generate_invivo_invitro_comparison_figures", build_text, fixed = TRUE)
  )
  testthat::expect_match(
    build_text,
    'file.path\\(fit_dir, "viz", "invivo_vs_invitro"\\)'
  )
  testthat::expect_false(
    grepl("plot_joint_parameter_ratio_figure", text, fixed = TRUE)
  )
})

testthat::test_that("post-fit runner preserves producer-to-consumer order", {
  text <- read_stage2_script("runner", "run_postfit_pipeline.R")
  markers <- c(
    '"in-vivo simulation materialization"',
    '"in-vitro simulation materialization"',
    '"in-vitro fit diagnostics"',
    '"joint parameter diagnostics"',
    '"in-vivo visualization"',
    '"in-vitro visualization"',
    '"joint in-vivo/in-vitro visualization"',
    '"joint parameter visualization"',
    '"fit report"'
  )
  positions <- vapply(
    markers,
    function(marker) regexpr(marker, text, fixed = TRUE)[1L],
    integer(1)
  )
  testthat::expect_true(all(positions > 0L))
  testthat::expect_true(all(diff(positions) > 0L))
  testthat::expect_match(text, "--scope must be one of: invivo, invitro, joint")
  testthat::expect_match(text, "dry_run")
})

testthat::test_that("fit backends dispatch one unified post-fit pipeline", {
  backend_files <- c(
    "o2_supply_demand_map_fit_invitro_backend.R",
    "o2_supply_demand_map_fit_invivo_backend.R",
    "o2_supply_demand_map_fit_joint_backend.R"
  )
  texts <- vapply(
    backend_files,
    function(path) read_stage2_script("util", path),
    character(1)
  )
  testthat::expect_true(all(grepl("run_postfit_pipeline.R", texts, fixed = TRUE)))
  testthat::expect_match(texts[[1]], "--scope=invitro", fixed = TRUE)
  testthat::expect_match(texts[[2]], "--scope=invivo", fixed = TRUE)
  testthat::expect_match(texts[[3]], "--scope=joint", fixed = TRUE)
})

testthat::test_that("HPC joint rerender uses the unified post-fit worker", {
  text <- read_stage2_script(
    "hpc",
    "postprocess",
    "submit_rerender_joint_seed_reports.sh"
  )
  testthat::expect_match(text, "runner/run_postfit_pipeline.R", fixed = TRUE)
  testthat::expect_match(text, "--scope=joint", fixed = TRUE)
  testthat::expect_false(
    grepl("INVIVO_VIZ_SCRIPT|INVITRO_VIZ_SCRIPT", text, perl = TRUE)
  )
})

testthat::test_that("shared path helper resolves canonical simulation tables", {
  env <- new.env(parent = baseenv())
  sys.source(
    file.path(
      stage2_postfit_root,
      "util",
      "o2_supply_demand_map_shared.R"
    ),
    envir = env
  )
  fit_dir <- tempfile("o2sd_postfit_path_")
  dir.create(fit_dir)
  on.exit(unlink(fit_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expected <- file.path(
    normalizePath(fit_dir, mustWork = TRUE),
    "simulation",
    "invivo",
    "burden_timecourse.tsv"
  )
  testthat::expect_identical(
    env$o2sd_simulation_table_path(
      fit_dir,
      "burden_timecourse.tsv",
      "invivo",
      must_work = FALSE
    ),
    expected
  )
  testthat::expect_error(
    env$o2sd_simulation_table_path(
      fit_dir,
      "../best_params.tsv",
      "invivo",
      must_work = FALSE
    ),
    "one simulation-table basename",
    fixed = TRUE
  )
})
