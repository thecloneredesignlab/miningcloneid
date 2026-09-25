#!/usr/bin/env Rscript

# Consume existing dense fixed-O2 steady-state curves and their established
# regression-smoothed classifications. This analysis never runs the model.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

curve_class_levels <- function(observed = character()) {
  preferred <- c(
    "approximately_flat", "monotone_increasing", "monotone_decreasing",
    "single_transition_increase_then_plateau",
    "single_transition_decrease_then_plateau", "u_shaped",
    "inverted_u_shaped", "complex_nonmonotone", "insufficient_data"
  )
  c(intersect(preferred, observed), sort(setdiff(observed, preferred)))
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  normalizePath(path, mustWork = TRUE)
}

group_median <- function(data, keys, value, output_name, count_col = "seed_id") {
  groups <- o2jca_group_split(data, keys)
  rows <- lapply(groups, function(group) {
    out <- group[1L, keys, drop = FALSE]
    values <- suppressWarnings(as.numeric(group[[value]]))
    values <- values[is.finite(values)]
    out[[output_name]] <- if (length(values)) stats::median(values) else NA_real_
    out$n_curve <- length(unique(group[[count_col]]))
    out
  })
  do.call(rbind, rows)
}

normalized_entropy_k <- function(counts, k) {
  counts <- as.numeric(counts)
  counts <- counts[is.finite(counts) & counts > 0]
  if (!length(counts) || sum(counts) <= 0) return(NA_real_)
  if (k <= 1L || length(counts) == 1L) return(0)
  p <- counts / sum(counts)
  -sum(p * log(p)) / log(k)
}

required_arg <- function(args, key) {
  value <- args[[key]] %||% ""
  if (!nzchar(value)) stop("--", key, " is required", call. = FALSE)
  normalizePath(path.expand(value), mustWork = TRUE)
}

analyze_joint_fixed_o2 <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  args <- o2jca_parse_args(raw_args)
  curve_path <- required_arg(args, "curve_table")
  class_path <- required_arg(args, "class_table")
  manifest_path <- required_arg(args, "seed_manifest")
  category_path <- required_arg(args, "ploidy_category_table")
  out_dir <- normalizePath(path.expand(args$out_dir %||% stop("--out_dir is required", call. = FALSE)), mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  curves <- o2jca_read_tsv(curve_path)
  classes <- o2jca_read_tsv(class_path)
  manifest <- o2jca_read_tsv(manifest_path)
  categories <- o2jca_read_tsv(category_path)

  o2jca_assert_columns(
    curves,
    c("seed_id", "O2_pct", "dominant_mean_ploidy", "smoothed_dominant_mean_ploidy"),
    "fixed-O2 regression curve table"
  )
  o2jca_assert_columns(
    classes,
    c("seed_id", "smooth_curve_class", "smooth_classification_rule_version"),
    "fixed-O2 by-seed class table"
  )
  o2jca_assert_columns(
    manifest,
    c("synthetic_seed_id", "synthetic_seed_number", "pair_id", "joint_seed", "source_seed_dir"),
    "joint synthetic-seed manifest"
  )
  o2jca_assert_columns(
    categories,
    c("pair_id", "pair_label", "pair_ploidy_category"),
    "pair ploidy-category table"
  )
  o2jca_assert_unique_key(curves, c("seed_id", "O2_pct"), "fixed-O2 regression curve table")
  o2jca_assert_unique_key(classes, "seed_id", "fixed-O2 by-seed class table")
  o2jca_assert_unique_key(manifest, "synthetic_seed_id", "joint synthetic-seed manifest")
  o2jca_assert_unique_key(categories, "pair_id", "pair ploidy-category table")

  manifest$pair_id <- o2jca_normalize_pair_id(manifest$pair_id)
  categories$pair_id <- o2jca_normalize_pair_id(categories$pair_id)
  classes$curve_class <- as.character(classes$smooth_curve_class)
  observed_classes <- curve_class_levels(unique(classes$curve_class))
  classes$curve_class <- factor(classes$curve_class, levels = observed_classes)

  manifest_keep <- unique(c(
    "synthetic_seed_id", "synthetic_seed_number", "pair_id", "joint_seed",
    "source_seed_dir", intersect(c("synthetic_seed_dir"), names(manifest))
  ))
  seed_meta <- merge(
    classes,
    manifest[, manifest_keep, drop = FALSE],
    by.x = "seed_id", by.y = "synthetic_seed_id", all.x = TRUE, sort = FALSE
  )
  if (anyNA(seed_meta$pair_id)) stop("Not every classified seed maps to a pair", call. = FALSE)
  seed_meta <- merge(
    seed_meta,
    categories[, c("pair_id", "pair_label", "pair_ploidy_category"), drop = FALSE],
    by = "pair_id", all.x = TRUE, sort = FALSE
  )
  if (anyNA(seed_meta$pair_ploidy_category)) stop("Not every pair maps to a temporal ploidy Cat", call. = FALSE)
  if (nrow(seed_meta) != nrow(classes)) stop("Seed mapping changed classified seed count", call. = FALSE)
  o2jca_assert_unique_key(seed_meta, "seed_id", "annotated fixed-O2 seed classes")

  curve_keep <- intersect(
    c(
      "seed_id", "O2_pct", "dominant_mean_ploidy", "smoothed_dominant_mean_ploidy",
      "spectral_gap", "dominant_growth_rate", "objective", "trajectory_regime",
      "mode_label", "is_reporting_o2"
    ),
    names(curves)
  )
  curve_meta_keep <- c(
    "seed_id", "synthetic_seed_number", "joint_seed", "pair_id", "pair_label",
    "pair_ploidy_category", "curve_class", "smooth_classification_rule_version",
    intersect(c("monotonicity_reliability", "final_interpretation_class"), names(seed_meta))
  )
  plot_data <- merge(
    curves[, curve_keep, drop = FALSE],
    seed_meta[, curve_meta_keep, drop = FALSE],
    by = "seed_id", all.x = TRUE, sort = FALSE
  )
  if (anyNA(plot_data$pair_id) || anyNA(plot_data$curve_class)) {
    stop("Not every fixed-O2 curve row maps to pair and curve class", call. = FALSE)
  }
  plot_data$curve_class <- factor(as.character(plot_data$curve_class), levels = observed_classes)
  plot_data <- plot_data[order(plot_data$synthetic_seed_number, as.numeric(plot_data$O2_pct)), , drop = FALSE]

  o2_grid <- sort(unique(suppressWarnings(as.numeric(plot_data$O2_pct))))
  if (!length(o2_grid) || any(!is.finite(o2_grid))) stop("Invalid fixed-O2 grid", call. = FALSE)
  grid_key <- paste(format(o2_grid, digits = 15, trim = TRUE), collapse = ",")
  grid_audit <- lapply(o2jca_group_split(plot_data, "seed_id"), function(group) {
    values <- sort(unique(suppressWarnings(as.numeric(group$O2_pct))))
    data.frame(
      seed_id = group$seed_id[[1L]], pair_id = group$pair_id[[1L]],
      n_o2 = length(values), o2_min = min(values), o2_max = max(values),
      grid_complete = identical(paste(format(values, digits = 15, trim = TRUE), collapse = ","), grid_key),
      stringsAsFactors = FALSE
    )
  })
  grid_audit <- do.call(rbind, grid_audit)
  if (!all(grid_audit$grid_complete)) stop("At least one seed has an incomplete fixed-O2 grid", call. = FALSE)

  pair_levels <- unique(categories$pair_id[categories$pair_id %in% seed_meta$pair_id])
  seed_meta$pair_id <- factor(seed_meta$pair_id, levels = pair_levels)
  seed_meta$curve_class <- factor(as.character(seed_meta$curve_class), levels = observed_classes)
  pair_cells <- expand.grid(
    pair_id = pair_levels, curve_class = observed_classes,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  counts <- as.data.frame(table(seed_meta$pair_id, seed_meta$curve_class), stringsAsFactors = FALSE)
  names(counts) <- c("pair_id", "curve_class", "n_seed")
  pair_cells <- merge(pair_cells, counts, by = c("pair_id", "curve_class"), all.x = TRUE, sort = FALSE)
  pair_cells$n_seed[is.na(pair_cells$n_seed)] <- 0L
  pair_totals <- aggregate(n_seed ~ pair_id, pair_cells, sum)
  names(pair_totals)[[2L]] <- "pair_total_seed"
  pair_cells <- merge(pair_cells, pair_totals, by = "pair_id", all.x = TRUE, sort = FALSE)
  pair_cells$fraction_seed <- pair_cells$n_seed / pair_cells$pair_total_seed
  pair_cells <- merge(
    pair_cells,
    categories[, c("pair_id", "pair_label", "pair_ploidy_category"), drop = FALSE],
    by = "pair_id", all.x = TRUE, sort = FALSE
  )
  pair_metrics <- do.call(rbind, lapply(o2jca_group_split(pair_cells, "pair_id"), function(group) {
    winner <- group$curve_class[[which.max(group$n_seed)]]
    ordered <- sort(group$n_seed, decreasing = TRUE)
    data.frame(
      pair_id = group$pair_id[[1L]], dominant_curve_class = winner,
      dominant_fraction = max(group$fraction_seed),
      dominant_margin = if (length(ordered) > 1L) (ordered[[1L]] - ordered[[2L]]) / sum(ordered) else 1,
      normalized_curve_class_entropy = normalized_entropy_k(group$n_seed, length(observed_classes)),
      n_observed_curve_classes = sum(group$n_seed > 0), stringsAsFactors = FALSE
    )
  }))
  pair_cells <- merge(pair_cells, pair_metrics, by = "pair_id", all.x = TRUE, sort = FALSE)
  pair_cells <- pair_cells[order(match(pair_cells$pair_id, pair_levels), match(pair_cells$curve_class, observed_classes)), , drop = FALSE]

  cat_levels <- c("CatA", "CatB", "CatC", "CatU")
  cat_cells <- expand.grid(
    pair_ploidy_category = intersect(cat_levels, unique(pair_cells$pair_ploidy_category)),
    curve_class = observed_classes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  cat_summary <- do.call(rbind, lapply(o2jca_group_split(pair_cells, c("pair_ploidy_category", "curve_class")), function(group) {
    data.frame(
      pair_ploidy_category = group$pair_ploidy_category[[1L]],
      curve_class = group$curve_class[[1L]], n_pair = nrow(group),
      pair_balanced_mean_fraction = mean(group$fraction_seed),
      pair_balanced_median_fraction = stats::median(group$fraction_seed),
      pair_fraction_min = min(group$fraction_seed), pair_fraction_max = max(group$fraction_seed),
      pooled_n_seed = sum(group$n_seed), pooled_total_seed = sum(group$pair_total_seed),
      pooled_fraction_seed = sum(group$n_seed) / sum(group$pair_total_seed),
      stringsAsFactors = FALSE
    )
  }))
  cat_cells <- merge(cat_cells, cat_summary, by = c("pair_ploidy_category", "curve_class"), all.x = TRUE, sort = FALSE)
  cat_cells <- cat_cells[order(match(cat_cells$pair_ploidy_category, cat_levels), match(cat_cells$curve_class, observed_classes)), , drop = FALSE]

  pair_median <- group_median(
    plot_data,
    c("pair_id", "pair_label", "pair_ploidy_category", "curve_class", "O2_pct"),
    "smoothed_dominant_mean_ploidy", "median_smoothed_dominant_mean_ploidy"
  )
  cat_pair_median <- group_median(
    pair_median,
    c("pair_ploidy_category", "curve_class", "O2_pct"),
    "median_smoothed_dominant_mean_ploidy", "pair_balanced_median_smoothed_dominant_mean_ploidy",
    count_col = "pair_id"
  )
  names(cat_pair_median)[names(cat_pair_median) == "n_curve"] <- "n_pair"

  pair_dominant <- merge(
    unique(pair_cells[, c("pair_id", "pair_label", "pair_ploidy_category")]),
    pair_metrics, by = "pair_id", all.x = TRUE, sort = FALSE
  )
  pair_tab <- table(pair_dominant$pair_ploidy_category, pair_dominant$dominant_curve_class)
  seed_tab <- table(seed_meta$pair_ploidy_category, seed_meta$curve_class)
  association <- data.frame(
    analysis_unit = c("pair_dominant_class", "seed_class_pooled_descriptive"),
    n_unit = c(nrow(pair_dominant), nrow(seed_meta)),
    cramers_v = c(o2jca_cramers_v(pair_tab), o2jca_cramers_v(seed_tab)),
    normalized_mutual_information = c(o2jca_normalized_mutual_information(pair_tab), o2jca_normalized_mutual_information(seed_tab)),
    inferential_scope = c(
      paste0(
        "descriptive: ", nrow(pair_dominant), " pairs; Cat pair counts ",
        paste(names(table(pair_dominant$pair_ploidy_category)), as.integer(table(pair_dominant$pair_ploidy_category)), sep = "=", collapse = ", ")
      ),
      "descriptive only: seeds are repeated optimization outcomes nested within pairs"
    ),
    stringsAsFactors = FALSE
  )

  pair_quality <- do.call(rbind, lapply(o2jca_group_split(grid_audit, "pair_id"), function(group) {
    data.frame(
      pair_id = group$pair_id[[1L]], n_seed = nrow(group),
      n_o2 = unique(group$n_o2)[[1L]], o2_min = unique(group$o2_min)[[1L]],
      o2_max = unique(group$o2_max)[[1L]], all_seed_grids_complete = all(group$grid_complete),
      expected_curve_rows = nrow(group) * unique(group$n_o2)[[1L]],
      observed_curve_rows = sum(plot_data$pair_id == group$pair_id[[1L]]),
      stringsAsFactors = FALSE
    )
  }))
  pair_quality$row_count_pass <- pair_quality$expected_curve_rows == pair_quality$observed_curve_rows
  pair_quality <- merge(pair_quality, categories[, c("pair_id", "pair_label", "pair_ploidy_category")], by = "pair_id", all.x = TRUE, sort = FALSE)

  seed_path <- o2jca_write_tsv(seed_meta, file.path(out_dir, "fixed_o2_curve_class_by_seed.tsv"))
  plot_path <- write_tsv_gz(plot_data, file.path(out_dir, "fixed_o2_regression_smoothed_curve_plot_data.tsv.gz"))
  pair_summary_path <- o2jca_write_tsv(pair_cells, file.path(out_dir, "fixed_o2_curve_class_summary_by_pair.tsv"))
  cat_summary_path <- o2jca_write_tsv(cat_cells, file.path(out_dir, "fixed_o2_curve_class_summary_by_cat.tsv"))
  pair_dominant_path <- o2jca_write_tsv(pair_dominant, file.path(out_dir, "fixed_o2_pair_dominant_curve_class.tsv"))
  association_path <- o2jca_write_tsv(association, file.path(out_dir, "fixed_o2_curve_class_association.tsv"))
  pair_median_path <- o2jca_write_tsv(pair_median, file.path(out_dir, "fixed_o2_curve_pair_class_median.tsv"))
  cat_median_path <- o2jca_write_tsv(cat_pair_median, file.path(out_dir, "fixed_o2_curve_cat_pair_balanced_median.tsv"))
  quality_path <- o2jca_write_tsv(pair_quality, file.path(out_dir, "fixed_o2_input_quality_summary.tsv"))
  provenance <- data.frame(
    input_role = c("regression_curves", "regression_by_seed", "synthetic_seed_manifest", "temporal_ploidy_pair_category"),
    path = normalizePath(c(curve_path, class_path, manifest_path, category_path), mustWork = TRUE),
    bytes = as.numeric(file.info(c(curve_path, class_path, manifest_path, category_path))$size),
    stringsAsFactors = FALSE
  )
  provenance_path <- o2jca_write_tsv(provenance, file.path(out_dir, "fixed_o2_input_provenance.tsv"))
  outputs <- c(
    seed_path, plot_path, pair_summary_path, cat_summary_path, pair_dominant_path,
    association_path, pair_median_path, cat_median_path, quality_path, provenance_path
  )
  o2jca_write_manifest(outputs, file.path(out_dir, "fixed_o2_ploidy_classification_manifest.tsv"))
  message(
    "Wrote joint fixed-O2 classification tables: ", nrow(seed_meta), " seeds, ",
    nrow(plot_data), " curve rows, ", length(o2_grid), " O2 points, ",
    length(observed_classes), " observed curve classes"
  )
  invisible(outputs)
}

if (identical(environment(), globalenv())) analyze_joint_fixed_o2()
