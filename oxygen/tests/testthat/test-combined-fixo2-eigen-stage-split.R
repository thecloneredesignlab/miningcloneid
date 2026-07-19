combined_fixo2_root <- file.path(repo_info$root, "oxygen", "code", "O2_supply_demand_MAP")

combined_fixo2_files <- list(
  slope = file.path(combined_fixo2_root, "analysis", "combined_fixo2_eigen", "calculate_regression_curve_average_slope.R"),
  prepare = file.path(combined_fixo2_root, "analysis", "combined_fixo2_eigen", "prepare_fixo2_eigen_curve_class_tables.R"),
  vis = file.path(combined_fixo2_root, "vis", "combined_fixo2_eigen", "plot_fixo2_eigen_attractor_embedding_curve_class.R"),
  report = file.path(combined_fixo2_root, "report", "combined_fixo2_eigen", "render_fixo2_eigen_attractor_embedding_curve_class_report.R"),
  runner = file.path(combined_fixo2_root, "runner", "combined_fixo2_eigen", "run_fixo2_eigen_curve_class_pipeline.R"),
  legacy = file.path(
    combined_fixo2_root,
    "analysis", "best_fit_parameter_feature", "06_combine_FixO2_eigen_attractor",
    c(
      "calculate_regression_curve_average_slope.R",
      "plot_fixo2_eigen_attractor_embedding_curve_class.R",
      "render_fixo2_eigen_attractor_embedding_curve_class_report.R"
    )
  )
)

combined_fixo2_code <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  paste(lines, collapse = "\n")
}

combined_fixo2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

testthat::test_that("combined FixO2 eigen-attractor entrypoints exist and parse", {
  paths <- c(unlist(combined_fixo2_files[setdiff(names(combined_fixo2_files), "legacy")]), combined_fixo2_files$legacy)
  testthat::expect_true(all(file.exists(paths)))
  for (path in paths) testthat::expect_silent(parse(path))
})

testthat::test_that("combined FixO2 eigen-attractor responsibilities are separated", {
  analysis_text <- paste(
    combined_fixo2_code(combined_fixo2_files$slope),
    combined_fixo2_code(combined_fixo2_files$prepare),
    sep = "\n"
  )
  vis_text <- combined_fixo2_code(combined_fixo2_files$vis)
  report_text <- combined_fixo2_code(combined_fixo2_files$report)

  testthat::expect_false(grepl("ggplot2::|ggsave\\s*\\(|grDevices::(pdf|png)", analysis_text, perl = TRUE))
  testthat::expect_match(analysis_text, "pooled_embedding_curve_class_analysis_manifest.tsv", fixed = TRUE)
  testthat::expect_match(analysis_text, "annotated_coordinate_table", fixed = TRUE)

  for (forbidden in c(
    "find_latest_class_table", "read_curve_class_map", "read_average_slope_map",
    "add_curve_class", "add_average_slope", "discover_coordinate_tables",
    "--pooled_root", "--dense_grid_dir", "--class_table", "--average_slope_table"
  )) {
    testthat::expect_false(grepl(forbidden, vis_text, fixed = TRUE), info = forbidden)
  }
  testthat::expect_match(vis_text, "annotated_coordinate_table", fixed = TRUE)
  testthat::expect_match(vis_text, "pooled_embedding_curve_class_manifest.tsv", fixed = TRUE)

  testthat::expect_false(grepl("ggplot2::|ggsave\\s*\\(|readRDS\\s*\\(", report_text, perl = TRUE))
  testthat::expect_match(report_text, "pooled_embedding_curve_class_manifest.tsv", fixed = TRUE)
})

testthat::test_that("combined FixO2 eigen-attractor runner preserves stage order", {
  text <- combined_fixo2_code(combined_fixo2_files$runner)
  stages <- c(
    "calculate_regression_curve_average_slope.R",
    "prepare_fixo2_eigen_curve_class_tables.R",
    "plot_fixo2_eigen_attractor_embedding_curve_class.R",
    "render_fixo2_eigen_attractor_embedding_curve_class_report.R"
  )
  positions <- vapply(stages, function(stage) regexpr(stage, text, fixed = TRUE)[[1L]], integer(1))
  testthat::expect_true(all(positions > 0L))
  testthat::expect_true(all(diff(positions) > 0L))
})

testthat::test_that("historical 06 entrypoints are thin compatibility wrappers", {
  for (path in combined_fixo2_files$legacy) {
    lines <- readLines(path, warn = FALSE)
    testthat::expect_true(length(lines) <= 250L, info = path)
    testthat::expect_true(
      grepl("deprecated compatibility", paste(head(lines, 30L), collapse = "\n"), ignore.case = TRUE),
      info = path
    )
  }
})

testthat::test_that("slope and annotated-table analysis preserve the numeric contract on a fixture", {
  root <- tempfile("combined_fixo2_eigen_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  curves <- data.frame(
    seed_id = c("seed1", "seed1", "seed1", "seed2", "seed2"),
    O2_pct = c(0, 1, 2, 0, 2),
    smoothed_dominant_mean_ploidy = c(1, 2, 3, 3, 1),
    smooth_step_epsilon = 0
  )
  by_seed <- data.frame(
    seed_id = c("seed1", "seed2"),
    seed_number = 1:2,
    curve_class = c("monotone_increasing", "monotone_decreasing"),
    objective = c(0.1, 0.2)
  )
  curve_path <- file.path(root, "curves.tsv")
  by_seed_path <- file.path(root, "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv")
  combined_fixo2_write_tsv(curves, curve_path)
  combined_fixo2_write_tsv(by_seed, by_seed_path)
  out_dir <- file.path(root, "out")
  rscript <- file.path(R.home("bin"), "Rscript")
  slope_status <- system2(rscript, c(
    shQuote(combined_fixo2_files$slope),
    shQuote(paste0("--curve_table=", curve_path)),
    shQuote(paste0("--by_seed_table=", by_seed_path)),
    shQuote(paste0("--out_dir=", file.path(out_dir, "tables")))
  ))
  testthat::expect_identical(slope_status, 0L)
  slope_path <- file.path(out_dir, "tables", "fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv")
  slopes <- utils::read.delim(slope_path, check.names = FALSE)
  testthat::expect_equal(slopes$average_slope, c(1, -1), tolerance = 0)
  testthat::expect_identical(slopes$curve_class, by_seed$curve_class)

  pooled_root <- file.path(root, "pooled")
  coordinate_dir <- file.path(pooled_root, "PCAs", "Tables", "BestOnly")
  dir.create(coordinate_dir, recursive = TRUE)
  coordinates <- data.frame(
    dataset = c("invivo", "invivo", "invitro", "invivo"),
    point_type = c("best", "best", "best", "initial"),
    seed = c(1, 2, 1, 1),
    PCA1 = c(0, 1, 2, 3),
    PCA2 = c(3, 2, 1, 0),
    objective = c(0.1, 0.2, 0.3, NA_real_)
  )
  coordinate_path <- file.path(coordinate_dir, "toy_coordinates.csv")
  utils::write.csv(coordinates, coordinate_path, row.names = FALSE, na = "")
  prepare_status <- system2(rscript, c(
    shQuote(combined_fixo2_files$prepare),
    shQuote(paste0("--pooled_root=", pooled_root)),
    shQuote(paste0("--dense_grid_dir=", root)),
    shQuote(paste0("--class_table=", by_seed_path)),
    shQuote(paste0("--average_slope_table=", slope_path)),
    shQuote(paste0("--out_dir=", out_dir)),
    "--reductions=PCAs",
    "--variants=BestOnly"
  ))
  testthat::expect_identical(prepare_status, 0L)
  manifest <- utils::read.delim(
    file.path(out_dir, "tables", "pooled_embedding_curve_class_analysis_manifest.tsv"),
    check.names = FALSE
  )
  testthat::expect_equal(nrow(manifest), 1L)
  annotated <- utils::read.csv(manifest$annotated_coordinate_table[[1L]], check.names = FALSE)
  testthat::expect_identical(annotated$curve_class[1:2], by_seed$curve_class)
  testthat::expect_equal(annotated$average_slope[1:2], c(1, -1), tolerance = 0)
  testthat::expect_true(all(is.na(annotated$curve_class[3:4]) | !nzchar(annotated$curve_class[3:4])))
  testthat::expect_true(file.exists(manifest$best_points_table[[1L]]))
  testthat::expect_true(file.exists(manifest$curve_class_counts_table[[1L]]))
})
