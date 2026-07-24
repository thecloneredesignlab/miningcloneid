find_joint_stage_repo_root <- function() {
  current <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:8) {
    target <- file.path(current, "oxygen", "code", "O2_supply_demand_MAP", "README.md")
    if (file.exists(target)) return(current)
    parent <- dirname(current); if (identical(parent, current)) break; current <- parent
  }
  stop("Cannot locate repository root")
}

testthat::test_that("joint coupling stages preserve analysis-vis-report-HPC ownership", {
  root <- find_joint_stage_repo_root()
  workflow <- file.path(root, "oxygen", "code", "O2_supply_demand_MAP")
  files <- list(
    ratio_analysis = list.files(file.path(workflow, "analysis", "joint_soft_coupling_stability"), pattern = "[.]R$", full.names = TRUE),
    ploidy_analysis = list.files(file.path(workflow, "analysis", "joint_ploidy_coupling_association"), pattern = "[.]R$", full.names = TRUE),
    fixed_o2_analysis = list.files(file.path(workflow, "analysis", "joint_fixed_o2_ploidy_classification"), pattern = "[.]R$", full.names = TRUE),
    ratio_vis = list.files(file.path(workflow, "vis", "joint_soft_coupling_stability"), pattern = "[.]R$", full.names = TRUE),
    ploidy_vis = list.files(file.path(workflow, "vis", "joint_ploidy_coupling_association"), pattern = "[.]R$", full.names = TRUE),
    fixed_o2_vis = list.files(file.path(workflow, "vis", "joint_fixed_o2_ploidy_classification"), pattern = "[.]R$", full.names = TRUE)
  )
  read_all <- function(paths) paste(unlist(lapply(paths, readLines, warn = FALSE)), collapse = "\n")
  analysis_text <- read_all(c(files$ratio_analysis, files$ploidy_analysis, files$fixed_o2_analysis))
  vis_text <- read_all(c(files$ratio_vis, files$ploidy_vis, files$fixed_o2_vis))
  testthat::expect_false(grepl("ggplot2::|ggsave\\s*\\(", analysis_text))
  testthat::expect_false(grepl("readRDS\\s*\\(|best_params[.]tsv|fit_result[.]rds", vis_text, ignore.case = TRUE))
  testthat::expect_true(file.exists(file.path(workflow, "runner", "joint_coupling", "run_joint_coupling_pipeline.R")))
  testthat::expect_true(file.exists(file.path(workflow, "hpc", "joint_coupling_analysis", "submit_joint_coupling_analysis.sh")))
  testthat::expect_true(file.exists(file.path(workflow, "report", "joint_coupling", "render_joint_coupling_report.R")))
})

testthat::test_that("joint coupling output contract keeps the fitting tree read-only", {
  root <- find_joint_stage_repo_root()
  workflow <- file.path(root, "oxygen", "code", "O2_supply_demand_MAP")
  util_env <- new.env(parent = baseenv())
  sys.source(file.path(workflow, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"), envir = util_env)
  result_root <- file.path(root, "oxygen", "results", "example_joint_fit")
  output_root <- util_env$o2jca_default_output_root(result_root, workflow)
  testthat::expect_identical(
    output_root,
    file.path(root, "oxygen", "results", "analysis", "joint_coupling", "example_joint_fit")
  )
  testthat::expect_false(util_env$o2jca_path_is_within(output_root, result_root))
  testthat::expect_identical(
    util_env$o2jca_path_is_within(c(file.path(result_root, "analysis"), output_root), result_root),
    c(TRUE, FALSE)
  )
  testthat::expect_error(
    util_env$o2jca_assert_separate_output_root(result_root, file.path(result_root, "analysis")),
    "outside the read-only fitting result_root"
  )

  runner_text <- paste(readLines(file.path(workflow, "runner", "joint_coupling", "run_joint_coupling_pipeline.R"), warn = FALSE), collapse = "\n")
  submit_text <- paste(readLines(file.path(workflow, "hpc", "joint_coupling_analysis", "submit_joint_coupling_analysis.sh"), warn = FALSE), collapse = "\n")
  worker_text <- paste(readLines(file.path(workflow, "hpc", "joint_coupling_analysis", "run_joint_coupling_analysis.sub"), warn = FALSE), collapse = "\n")
  testthat::expect_match(runner_text, "--output_root|args\\$output_root")
  testthat::expect_false(grepl('file[.]path\\(result_root, "(analysis|vis|report|logs)"', runner_text))
  testthat::expect_match(submit_text, "OUTPUT_ROOT")
  testthat::expect_match(submit_text, 'log_dir="\\$\\{OUTPUT_ROOT\\}/logs"')
  testthat::expect_match(worker_text, '--output_root="\\$\\{OUTPUT_ROOT\\}"')
  testthat::expect_match(runner_text, "fixed_o2_source_root")
  testthat::expect_match(submit_text, "FIXED_O2_SOURCE_ROOT")
  testthat::expect_match(worker_text, "fixed_o2_source_root")
})

testthat::test_that("expanded visualization suite covers multiple chart families", {
  root <- find_joint_stage_repo_root()
  workflow <- file.path(root, "oxygen", "code", "O2_supply_demand_MAP")
  vis_files <- c(
    list.files(file.path(workflow, "vis", "joint_soft_coupling_stability"), pattern = "[.]R$", full.names = TRUE),
    list.files(file.path(workflow, "vis", "joint_ploidy_coupling_association"), pattern = "[.]R$", full.names = TRUE),
    list.files(file.path(workflow, "vis", "joint_fixed_o2_ploidy_classification"), pattern = "[.]R$", full.names = TRUE)
  )
  text <- paste(unlist(lapply(vis_files, readLines, warn = FALSE)), collapse = "\n")
  expected_stems <- c(
    "overview_pair_parameter_ratio_map", "within_pair_ratio_distributions",
    "between_pair_all_class_rankings", "stable90_pair_set_membership",
    "ploidy_category_archetypes", "ploidy_representative_trajectories",
    "cat_parameter_vivo_vitro_dumbbells",
    "fixed_o2_regression_smoothed_all_seed_curves_by_pair_and_class",
    "fixed_o2_regression_smoothed_all_pair_curves_by_class_and_cat"
  )
  testthat::expect_true(all(vapply(expected_stems, grepl, logical(1L), x = text, fixed = TRUE)))
  testthat::expect_match(text, "Matrix & Cohort")
  testthat::expect_match(text, "Distribution")
  testthat::expect_match(text, "Uncertainty & Benchmark")
})

testthat::test_that("biological-process facet labels are wrapped and horizontal", {
  root <- find_joint_stage_repo_root()
  script_path <- file.path(
    root, "oxygen", "code", "O2_supply_demand_MAP", "vis",
    "joint_soft_coupling_stability", "plot_soft_coupling_processes.R"
  )
  text <- paste(readLines(script_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "label_wrap_gen\\(width = 18\\)")
  testthat::expect_match(text, "strip[.]text[.]y = ggplot2::element_text\\(angle = 0")
  testthat::expect_match(text, '"process_parameter_class_map"[^\n]+11[.]5, 9')
})
