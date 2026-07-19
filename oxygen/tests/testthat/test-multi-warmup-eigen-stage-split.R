stage_split_find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = FALSE)
  for (depth in 0:10) {
    candidate <- normalizePath(
      file.path(current, paste(rep("..", depth), collapse = "/")),
      mustWork = FALSE
    )
    marker <- file.path(
      candidate,
      "oxygen",
      "code",
      "O2_supply_demand_MAP",
      "util",
      "o2_supply_demand_map_multi_warmup_utils.R"
    )
    if (file.exists(marker)) return(candidate)
  }
  stop("Cannot locate the repository root for stage-split tests.", call. = FALSE)
}

stage_split_repo_root <- if (exists("repo_info", inherits = TRUE)) {
  get("repo_info", inherits = TRUE)$root
} else {
  stage_split_find_repo_root()
}
stage_split_root <- file.path(
  stage_split_repo_root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

stage_split_path <- function(...) file.path(stage_split_root, ...)
stage_split_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

stage_split_leaf_files <- c(
  stage_split_path("util", "o2_supply_demand_map_eigen_attractor_utils.R"),
  stage_split_path("util", "o2_supply_demand_map_multi_warmup_utils.R"),
  stage_split_path(
    "simulation", "o2", "fixed_o2", "eigen_attractor",
    "build_fixo2_eigen_attractor_features.R"
  ),
  stage_split_path(
    "analysis", "fixed_o2_eigen",
    "analyze_fixo2_eigen_attractor_embeddings.R"
  ),
  stage_split_path("analysis", "multi_warmup", "build_multi_warmup_seed_plan_tables.R"),
  stage_split_path("analysis", "multi_warmup", "build_multi_warmup_landscape_tables.R"),
  stage_split_path("analysis", "multi_warmup", "collect_multi_warmup_tables.R"),
  stage_split_path("analysis", "multi_warmup", "analyze_warm_up_joint_results.R"),
  stage_split_path("vis", "fixed_o2_eigen", "render_fixo2_eigen_attractor_figures.R"),
  stage_split_path("vis", "multi_warmup", "render_multi_warmup_seed_plan_figures.R"),
  stage_split_path("vis", "multi_warmup", "render_multi_warmup_landscape_figures.R"),
  stage_split_path("vis", "multi_warmup", "render_multi_warmup_collected_figures.R"),
  stage_split_path("vis", "multi_warmup", "render_multi_warmup_report_figures.R"),
  stage_split_path("vis", "multi_warmup", "render_warm_up_joint_figures.R"),
  stage_split_path("report", "fixed_o2_eigen", "render_fixo2_eigen_attractor_report.R"),
  stage_split_path("report", "multi_warmup", "render_multi_warmup_results_report.R"),
  stage_split_path(
    "hpc", "fixo2_eigen_attractor", "fixo2_eigen_attractor_hpc_task.R"
  )
)

stage_split_legacy_files <- c(
  stage_split_path("analysis", "multi_warmup", "build_multi_warmup_seed_plan.R"),
  stage_split_path(
    "analysis", "multi_warmup",
    "build_multi_warmup_pairs_from_landscape_subclusters.R"
  ),
  stage_split_path("analysis", "multi_warmup", "collect_multi_warmup_results.R"),
  stage_split_path("analysis", "multi_warmup", "multi_warmup_results_report.R"),
  stage_split_path(
    "analysis", "warm_up_joint_fitting_results_extra",
    "warm_up_joint_fitting_results_extra.R"
  ),
  stage_split_path(
    "analysis", "best_fit_parameter_feature",
    "05_FixO2_eigen_attractor_based_clustering",
    "fixo2_eigen_attractor_feature_builder.R"
  ),
  stage_split_path(
    "analysis", "best_fit_parameter_feature",
    "05_FixO2_eigen_attractor_based_clustering",
    "fixo2_eigen_attractor_clustering_runner.R"
  ),
  stage_split_path(
    "analysis", "best_fit_parameter_feature",
    "05_FixO2_eigen_attractor_based_clustering",
    "fixo2_eigen_attractor_report.R"
  )
)

testthat::test_that("multi-warmup and fixed-O2 eigen stage files parse and source from any cwd", {
  testthat::expect_true(all(file.exists(c(stage_split_leaf_files, stage_split_legacy_files))))
  for (path in c(stage_split_leaf_files, stage_split_legacy_files)) {
    testthat::expect_silent(parse(path))
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tempdir())
  for (path in stage_split_leaf_files) {
    env <- new.env(parent = globalenv())
    testthat::expect_silent(source(path, local = env, chdir = FALSE))
  }
  for (path in stage_split_legacy_files) {
    env <- new.env(parent = globalenv())
    testthat::expect_silent(source(path, local = env, chdir = FALSE))
  }
})

testthat::test_that("historical multi-warmup and fixed-O2 eigen entries are thin delegates", {
  counts <- vapply(stage_split_legacy_files, function(path) {
    length(readLines(path, warn = FALSE))
  }, integer(1L))
  testthat::expect_true(all(counts < 80L))
  text <- vapply(stage_split_legacy_files, stage_split_text, character(1L))
  testthat::expect_true(all(grepl("compatibility", text, ignore.case = TRUE)))
  testthat::expect_true(all(grepl("runner|simulation|report", text, ignore.case = TRUE)))
})

testthat::test_that("fixed-O2 eigen numerical producer has no rendering behavior", {
  path <- stage_split_path(
    "simulation", "o2", "fixed_o2", "eigen_attractor",
    "build_fixo2_eigen_attractor_features.R"
  )
  text <- stage_split_text(path)
  forbidden <- c(
    "ggplot2::",
    "(^|[^A-Za-z0-9_.])(ggplot|ggsave)\\s*\\(",
    "(grDevices::)?(pdf|png|jpeg|svg)\\s*\\(",
    "source\\s*\\([^\n]*(analysis|vis|report)"
  )
  for (pattern in forbidden) {
    testthat::expect_false(grepl(pattern, text, perl = TRUE, ignore.case = TRUE), info = pattern)
  }
})

testthat::test_that("canonical multi-warmup and eigen analysis is table-only", {
  paths <- c(
    stage_split_path("analysis", "multi_warmup", "build_multi_warmup_seed_plan_tables.R"),
    stage_split_path("analysis", "multi_warmup", "build_multi_warmup_landscape_tables.R"),
    stage_split_path("analysis", "multi_warmup", "collect_multi_warmup_tables.R"),
    stage_split_path("analysis", "multi_warmup", "analyze_warm_up_joint_results.R"),
    stage_split_path("analysis", "fixed_o2_eigen", "analyze_fixo2_eigen_attractor_embeddings.R")
  )
  text <- paste(vapply(paths, stage_split_text, character(1L)), collapse = "\n")
  forbidden <- c(
    "ggplot2::",
    "(^|[^A-Za-z0-9_.])(ggplot|ggsave)\\s*\\(",
    "(grDevices::)?(pdf|png|jpeg|svg)\\s*\\(",
    "(source|sys[.]source)\\s*\\([^\n]*(simulation|vis)",
    "sourceCpp\\s*\\("
  )
  for (pattern in forbidden) {
    testthat::expect_false(grepl(pattern, text, perl = TRUE, ignore.case = TRUE), info = pattern)
  }
})

testthat::test_that("fixed-O2 eigen status summaries match materialized simulation tables", {
  table_root <- file.path(
    stage_split_repo_root,
    "oxygen", "results", "analysis", "best_fit_parameter_feature",
    "05_FixO2_eigen_attractor_based_clustering", "FixO2EigenAttractors", "Tables"
  )
  wide_path <- file.path(table_root, "invivo_best_fixo2_eigen_attractor_wide.csv")
  status_path <- file.path(table_root, "invivo_best_fixo2_eigen_attractor_status_summary.csv")
  if (!all(file.exists(c(wide_path, status_path)))) {
    testthat::skip("Materialized fixed-O2 eigen regression tables are not available.")
  }

  env <- new.env(parent = globalenv())
  sys.source(
    stage_split_path(
      "simulation", "o2", "fixed_o2", "eigen_attractor",
      "build_fixo2_eigen_attractor_features.R"
    ),
    envir = env,
    chdir = TRUE
  )
  wide <- utils::read.csv(wide_path, check.names = FALSE, stringsAsFactors = FALSE)
  expected <- utils::read.csv(status_path, check.names = FALSE, stringsAsFactors = FALSE)
  testthat::expect_identical(env$fixo2ea_status_summary(wide), expected)
})

testthat::test_that("multi-warmup table split preserves legacy materialized contracts", {
  result_root <- file.path(
    stage_split_repo_root,
    "oxygen", "results", "ver1", "fit_joint_multi_warmup_10seed_k5_sigma1"
  )
  manifest <- file.path(result_root, "multi_warmup_manifest.tsv")
  if (!file.exists(manifest)) {
    testthat::skip("Materialized multi-warmup regression result is not available.")
  }

  out_dir <- tempfile("multi_warmup_stage_split_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(
    stage_split_path("analysis", "multi_warmup", "collect_multi_warmup_tables.R"),
    envir = env,
    chdir = TRUE
  )
  testthat::expect_message(
    env$main(list(
      multi_warmup_root = result_root,
      manifest = manifest,
      out_dir = out_dir
    )),
    "Wrote multi-warmup summary"
  )

  exact_files <- c(
    "multi_warmup_deoptim_iteration_counts.tsv",
    "multi_warmup_final_basin_assignments.tsv",
    "multi_warmup_final_basin_distance_summary.tsv",
    "multi_warmup_final_to_warmup_distance_heatmap.tsv",
    "multi_warmup_final_to_warmup_distance_matrix.tsv",
    "multi_warmup_initial_final_distance_per_seed.tsv",
    "multi_warmup_invivo_only_parameter_long.tsv"
  )
  for (name in exact_files) {
    testthat::expect_identical(
      readLines(file.path(out_dir, name), warn = FALSE),
      readLines(file.path(result_root, name), warn = FALSE),
      info = name
    )
  }

  common_files <- c(
    "multi_warmup_best_seed_summary.tsv",
    "multi_warmup_deoptim_iteration_summary.tsv",
    "multi_warmup_parameter_ratio_long.tsv"
  )
  for (name in common_files) {
    expected <- utils::read.delim(
      file.path(result_root, name),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    actual <- utils::read.delim(
      file.path(out_dir, name),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    shared <- intersect(names(expected), names(actual))
    if ("joint_run_dir" %in% shared) {
      expected$joint_run_dir <- sub(
        "/oxygen/results/fit_joint",
        "/oxygen/results/ver1/fit_joint",
        expected$joint_run_dir,
        fixed = TRUE
      )
    }
    testthat::expect_identical(
      expected[, shared, drop = FALSE],
      actual[, shared, drop = FALSE],
      info = name
    )
  }
})

testthat::test_that("multi-warmup visualization functions render materialized tables", {
  testthat::skip_if_not_installed("ggplot2")
  out_dir <- tempfile("multi_warmup_vis_smoke_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  env <- new.env(parent = globalenv())
  sys.source(
    stage_split_path("vis", "multi_warmup", "render_multi_warmup_collected_figures.R"),
    envir = env,
    chdir = TRUE
  )
  summary <- data.frame(
    warmup_label = c("v01_i01", "v02_i01"),
    objective = c(2, 1),
    objective_invivo = c(1.1, 0.8),
    objective_invitro = c(0.9, 0.7),
    invivo_family = c("C01", "C02"),
    stringsAsFactors = FALSE
  )
  testthat::expect_silent(env$plot_objectives(summary, out_dir))
  testthat::expect_true(file.exists(file.path(out_dir, "multi_warmup_objective_summary.pdf")))
  testthat::expect_true(file.exists(file.path(out_dir, "multi_warmup_invivo_invitro_objective_scatter.pdf")))
})

testthat::test_that("seed-plan visualization preserves its materialized coordinate input", {
  testthat::skip_if_not_installed("ggplot2")
  out_dir <- tempfile("multi_warmup_seed_vis_smoke_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  coordinate_path <- file.path(
    out_dir,
    "joint_soft_coupling_18d_profile_umap_coords.tsv"
  )
  coordinates <- data.frame(
    pair_id = c("v01_i01", "v01_i02", "v02_i01", "v02_i02"),
    invivo_rank = c(1L, 1L, 2L, 2L),
    invitro_rank = c(1L, 2L, 1L, 2L),
    UMAP1 = c(-1, -0.7, 0.8, 1),
    UMAP2 = c(-0.8, 0.9, -0.7, 1),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    coordinates,
    coordinate_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  env <- new.env(parent = globalenv())
  sys.source(
    stage_split_path("vis", "multi_warmup", "render_multi_warmup_seed_plan_figures.R"),
    envir = env,
    chdir = TRUE
  )
  testthat::expect_silent(env$render_multi_warmup_seed_plan_figures(out_dir, umap_seed = 1L))
  testthat::expect_true(file.exists(coordinate_path))
  testthat::expect_true(file.exists(file.path(
    out_dir,
    "joint_soft_coupling_18d_profile_umap_500seed.png"
  )))
})
