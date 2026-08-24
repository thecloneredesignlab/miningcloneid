stage34_root <- file.path(repo_info$root, "oxygen", "code", "O2_supply_demand_MAP")

stage34_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

testthat::test_that("dense-grid numerical production and consumers have directional boundaries", {
  paths <- list(
    simulation = file.path(stage34_root, "simulation", "o2", "dense_grid", "generate_dense_grid_simulation_tables.R"),
    analysis = file.path(stage34_root, "analysis", "dense_grid_monotonicity", "analyze_dense_grid_tables.R"),
    visualization = file.path(stage34_root, "vis", "dense_grid_monotonicity", "render_dense_grid_figures.R"),
    report = file.path(stage34_root, "report", "dense_grid_monotonicity", "render_dense_grid_report.R"),
    runner = file.path(stage34_root, "runner", "dense_grid_monotonicity", "run_dense_grid_monotonicity.R")
  )
  testthat::expect_true(all(file.exists(unlist(paths))))
  testthat::expect_false(grepl("ggplot|png\\(|pdf\\(", stage34_text(paths$simulation), perl = TRUE))
  testthat::expect_false(grepl("simulation/|ggplot|png\\(|pdf\\(", stage34_text(paths$analysis), perl = TRUE))
  testthat::expect_false(grepl("best_params|fit_result|readRDS", stage34_text(paths$visualization), perl = TRUE))
  testthat::expect_false(grepl("ggplot|readRDS|best_params", stage34_text(paths$report), perl = TRUE))

  old_dirs <- c(
    file.path(stage34_root, "analysis", "best_fit_parameter_feature", "03_dense-grid_monotonicity_classification"),
    file.path(stage34_root, "analysis", "dense-grid_monotonicity_classification")
  )
  wrappers <- unlist(lapply(old_dirs, list.files, pattern = "[.]R$", full.names = TRUE))
  testthat::expect_true(length(wrappers) >= 10L)
  for (path in wrappers) {
    testthat::expect_true(length(readLines(path, warn = FALSE)) < 130L, info = path)
    testthat::expect_match(tolower(stage34_text(path)), "deprecated compatibility", fixed = TRUE, info = path)
  }
})

testthat::test_that("dense-grid pointwise classes preserve representative curve semantics", {
  env <- new.env(parent = globalenv())
  env$.dense_analysis_source_override <- file.path(stage34_root, "analysis", "dense_grid_monotonicity", "dense_grid_analysis_utils.R")
  sys.source(env$.dense_analysis_source_override, envir = env, chdir = FALSE)
  o2 <- 0:4
  curves <- rbind(
    data.frame(seed_id = "seed1", O2_pct = o2, dominant_mean_ploidy = 1 + o2 / 5, spectral_gap = 0.1, objective = 1),
    data.frame(seed_id = "seed2", O2_pct = o2, dominant_mean_ploidy = 2 - o2 / 5, spectral_gap = 0.1, objective = 2),
    data.frame(seed_id = "seed3", O2_pct = o2, dominant_mean_ploidy = c(2, 1.5, 1, 1.5, 2), spectral_gap = 0.1, objective = 3),
    data.frame(seed_id = "seed4", O2_pct = o2, dominant_mean_ploidy = c(1, 1.5, 2, 1.5, 1), spectral_gap = 0.1, objective = 4)
  )
  result <- env$dense_classify_attractor_curves(curves, list(flat_range_threshold = "0.01"), smooth = FALSE)
  observed <- stats::setNames(result$by_seed$curve_class, result$by_seed$seed_id)
  testthat::expect_identical(unname(observed[c("seed1", "seed2", "seed3", "seed4")]), c("monotone_increasing", "monotone_decreasing", "u_shaped", "inverted_u_shaped"))
})

testthat::test_that("combined workflow preserves seed-level average-slope contract", {
  env <- new.env(parent = globalenv())
  env$.combined_analysis_source_override <- file.path(stage34_root, "analysis", "combined_parameter_landscape", "prepare_combined_parameter_landscape_tables.R")
  sys.source(env$.combined_analysis_source_override, envir = env, chdir = FALSE)
  curves <- rbind(
    data.frame(seed_id = "seed1", O2_pct = c(0, 1, 2), smoothed_dominant_mean_ploidy = c(1, 2, 3), smooth_step_epsilon = 0.01),
    data.frame(seed_id = "seed2", O2_pct = c(0, 1, 2), smoothed_dominant_mean_ploidy = c(3, 2, 1), smooth_step_epsilon = 0.01)
  )
  annotations <- data.frame(seed_id = c("seed1", "seed2"), curve_class = c("monotone_increasing", "monotone_decreasing"))
  slopes <- env$compute_average_slope_table(curves, annotations)
  testthat::expect_equal(slopes$average_slope, c(1, -1), tolerance = 1e-12)
  testthat::expect_equal(slopes$weighted_mean_interval_slope, slopes$average_slope, tolerance = 1e-12)
  testthat::expect_identical(slopes$curve_class, annotations$curve_class)
})

testthat::test_that("old combined entrypoints are thin compatibility wrappers", {
  old_dirs <- c(
    file.path(stage34_root, "analysis", "best_fit_parameter_feature", "04_combine_parameter_landscape"),
    file.path(stage34_root, "analysis", "best_fit_parameter_feature", "04_combine")
  )
  wrappers <- unlist(lapply(old_dirs, list.files, pattern = "[.]R$", full.names = TRUE))
  testthat::expect_true(length(wrappers) >= 4L)
  for (path in wrappers) {
    testthat::expect_true(length(readLines(path, warn = FALSE)) < 90L, info = path)
    testthat::expect_match(tolower(stage34_text(path)), "deprecated compatibility", fixed = TRUE, info = path)
  }
})

testthat::test_that("combined analysis, visualization, and report run from materialized tables", {
  root <- tempfile("combined_stage34_")
  dense <- file.path(root, "dense", "tables")
  pooled <- file.path(root, "pooled", "PCAs", "Tables", "Sampled")
  out <- file.path(root, "out")
  dir.create(dense, recursive = TRUE); dir.create(pooled, recursive = TRUE)
  curves <- rbind(
    data.frame(seed_id = "seed1", O2_pct = c(0, 1, 2), smoothed_dominant_mean_ploidy = c(1, 2, 3), smooth_step_epsilon = 0.01),
    data.frame(seed_id = "seed2", O2_pct = c(0, 1, 2), smoothed_dominant_mean_ploidy = c(3, 2, 1), smooth_step_epsilon = 0.01)
  )
  classes <- data.frame(seed_id = c("seed1", "seed2"), seed_number = 1:2, curve_class = c("monotone_increasing", "monotone_decreasing"))
  utils::write.table(curves, file.path(dense, "fixed_o2_ploidy_monotonicity_regression_curves.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(classes, file.path(dense, "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  coordinates <- data.frame(PCA1 = c(-1, -0.5, 0.5, 1, -0.8, 0.8), PCA2 = c(-1, -0.5, 0.5, 1, 0.7, -0.7), dataset = c("invivo", "invivo", "invivo", "invivo", "invitro", "invitro"), point_type = c("initial", "initial", "best", "best", "best", "best"), seed = c(1, 2, 1, 2, 1, 2), objective = c(NA, NA, 1, 2, 1, 2))
  utils::write.csv(coordinates, file.path(pooled, "fixture_pca_coordinates.csv"), row.names = FALSE)
  scripts <- c(
    file.path(stage34_root, "analysis", "combined_parameter_landscape", "prepare_combined_parameter_landscape_tables.R"),
    file.path(stage34_root, "vis", "combined_parameter_landscape", "render_combined_parameter_landscape_figures.R"),
    file.path(stage34_root, "report", "combined_parameter_landscape", "render_combined_parameter_landscape_report.R")
  )
  statuses <- c(
    system2(file.path(R.home("bin"), "Rscript"), c(scripts[[1L]], paste0("--dense_grid_dir=", root, "/dense"), paste0("--pooled_root=", root, "/pooled"), paste0("--out_dir=", out), "--reductions=PCAs", "--variants=Sampled")),
    system2(file.path(R.home("bin"), "Rscript"), c(scripts[[2L]], paste0("--out_dir=", out))),
    system2(file.path(R.home("bin"), "Rscript"), c(scripts[[3L]], paste0("--out_dir=", out)))
  )
  testthat::expect_identical(as.integer(statuses), c(0L, 0L, 0L))
  testthat::expect_true(file.exists(file.path(out, "analysis_tables", "combined_embedding_table_manifest.tsv")))
  testthat::expect_true(file.exists(file.path(out, "figures", "combined_embedding_figure_manifest.tsv")))
  testthat::expect_true(file.exists(file.path(out, "report", "index.html")))
})

testthat::test_that("dense-grid HPC submitters use the canonical runner", {
  submitters <- c(
    file.path(stage34_root, "hpc", "dense_grid_monotonicity_classification", "submit_dense_grid_monotonicity_classification.sh"),
    file.path(stage34_root, "hpc", "best_fit_parameter_feature", "submit_best_fit_parameter_feature.sh"),
    file.path(stage34_root, "hpc", "warm_up_joint_fitting_results_extra", "submit_warm_up_joint_curve_array_hpc.sh"),
    file.path(stage34_root, "hpc", "submit", "submit_multi_warmup_joint.sh")
  )
  for (path in submitters) testthat::expect_match(stage34_text(path), "runner/dense_grid_monotonicity/run_dense_grid_monotonicity.R", fixed = TRUE, info = path)

  unified <- file.path(stage34_root, "hpc", "submit", "submit_o2_fit.sh")
  unified_text <- stage34_text(unified)
  testthat::expect_false(
    grepl("runner/dense_grid_monotonicity/run_dense_grid_monotonicity.R", unified_text, fixed = TRUE),
    info = "The only unified joint workflow must not run dense-grid curve filtering."
  )
  testthat::expect_match(unified_text, "primary-cluster workflow", fixed = TRUE)
})
