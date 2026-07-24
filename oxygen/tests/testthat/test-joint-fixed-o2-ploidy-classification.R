find_joint_fixed_o2_repo_root <- function() {
  current <- normalizePath(getwd(), mustWork = FALSE)
  for (i in seq_len(8L)) {
    target <- file.path(current, "oxygen", "code", "O2_supply_demand_MAP", "README.md")
    if (file.exists(target)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Cannot locate repository root")
}

write_joint_fixed_o2_fixture <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

testthat::test_that("joint fixed-O2 analysis preserves exact seed, pair, Cat, and grid mappings", {
  testthat::skip_if_not(nzchar(Sys.which("Rscript")), "Rscript is required")
  root <- find_joint_fixed_o2_repo_root()
  script <- file.path(
    root, "oxygen", "code", "O2_supply_demand_MAP", "analysis",
    "joint_fixed_o2_ploidy_classification", "analyze_joint_fixed_o2_ploidy_classes.R"
  )
  fixture <- tempfile("joint-fixed-o2-")
  curve_path <- file.path(fixture, "curves.tsv")
  class_path <- file.path(fixture, "classes.tsv")
  manifest_path <- file.path(fixture, "manifest.tsv")
  category_path <- file.path(fixture, "categories.tsv")
  out_dir <- file.path(fixture, "out")

  seeds <- paste0("seed", seq_len(4L))
  grid <- c(0, 1, 2)
  curves <- expand.grid(seed_id = seeds, O2_pct = grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  curves$dominant_mean_ploidy <- 44 - as.numeric(factor(curves$seed_id)) - curves$O2_pct
  curves$smoothed_dominant_mean_ploidy <- curves$dominant_mean_ploidy + 0.1
  classes <- data.frame(
    seed_id = seeds,
    smooth_curve_class = c("complex_nonmonotone", "complex_nonmonotone", "inverted_u_shaped", "inverted_u_shaped"),
    smooth_classification_rule_version = "loess_persistent_v1",
    stringsAsFactors = FALSE
  )
  manifest <- data.frame(
    synthetic_seed_id = seeds, synthetic_seed_number = seq_len(4L),
    pair_id = rep(c("pair1", "pair2"), each = 2L), joint_seed = rep(seq_len(2L), 2L),
    source_seed_dir = file.path(fixture, "source", seeds), stringsAsFactors = FALSE
  )
  categories <- data.frame(
    pair_id = c("fit_joint_pair1", "fit_joint_pair2"), pair_label = c("Pair 1", "Pair 2"),
    pair_ploidy_category = c("CatA", "CatB"), stringsAsFactors = FALSE
  )
  write_joint_fixed_o2_fixture(curves, curve_path)
  write_joint_fixed_o2_fixture(classes, class_path)
  write_joint_fixed_o2_fixture(manifest, manifest_path)
  write_joint_fixed_o2_fixture(categories, category_path)

  output <- system2(
    Sys.which("Rscript"),
    c(
      shQuote(script), paste0("--curve_table=", shQuote(curve_path)),
      paste0("--class_table=", shQuote(class_path)),
      paste0("--seed_manifest=", shQuote(manifest_path)),
      paste0("--ploidy_category_table=", shQuote(category_path)),
      paste0("--out_dir=", shQuote(out_dir))
    ),
    stdout = TRUE, stderr = TRUE
  )
  testthat::expect_true(is.null(attr(output, "status")), info = paste(output, collapse = "\n"))

  seed_table <- utils::read.delim(file.path(out_dir, "fixed_o2_curve_class_by_seed.tsv"), stringsAsFactors = FALSE)
  pair_table <- utils::read.delim(file.path(out_dir, "fixed_o2_curve_class_summary_by_pair.tsv"), stringsAsFactors = FALSE)
  quality <- utils::read.delim(file.path(out_dir, "fixed_o2_input_quality_summary.tsv"), stringsAsFactors = FALSE)
  manifest_out <- utils::read.delim(file.path(out_dir, "fixed_o2_ploidy_classification_manifest.tsv"), stringsAsFactors = FALSE)
  testthat::expect_equal(nrow(seed_table), 4L)
  testthat::expect_setequal(seed_table$pair_id, categories$pair_id)
  testthat::expect_setequal(seed_table$pair_ploidy_category, categories$pair_ploidy_category)
  testthat::expect_equal(nrow(pair_table), 4L)
  observed_counts <- stats::setNames(pair_table$n_seed, paste(pair_table$pair_id, pair_table$curve_class, sep = "::"))
  testthat::expect_equal(unname(observed_counts[["fit_joint_pair1::complex_nonmonotone"]]), 2L)
  testthat::expect_equal(unname(observed_counts[["fit_joint_pair1::inverted_u_shaped"]]), 0L)
  testthat::expect_equal(unname(observed_counts[["fit_joint_pair2::complex_nonmonotone"]]), 0L)
  testthat::expect_equal(unname(observed_counts[["fit_joint_pair2::inverted_u_shaped"]]), 2L)
  testthat::expect_true(all(quality$n_o2 == 3L))
  testthat::expect_true(all(quality$row_count_pass))
  testthat::expect_true(all(quality$all_seed_grids_complete))
  testthat::expect_equal(nrow(manifest_out), 10L)
})

testthat::test_that("joint fixed-O2 visualization uses Cat only for color", {
  root <- find_joint_fixed_o2_repo_root()
  script <- file.path(
    root, "oxygen", "code", "O2_supply_demand_MAP", "vis",
    "joint_fixed_o2_ploidy_classification", "plot_joint_fixed_o2_ploidy_curves.R"
  )
  text <- paste(readLines(script, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "colour = pair_ploidy_category")
  testthat::expect_match(text, "facet_grid")
  testthat::expect_match(text, "curve_class")
  testthat::expect_false(grepl("colour = curve_class|color = curve_class", text))
})
