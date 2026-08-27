#!/usr/bin/env Rscript

.combined_analysis_dir <- local({
  if (exists(".combined_analysis_source_override", inherits = TRUE)) return(dirname(normalizePath(get(".combined_analysis_source_override", inherits = TRUE), mustWork = TRUE)))
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "prepare_combined_parameter_landscape_tables.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.combined_root <- normalizePath(file.path(.combined_analysis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.combined_root, "util", "o2_supply_demand_map_combined_landscape_utils.R"), local = environment(), chdir = TRUE)

combined_first_finite <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (length(x)) x[[1L]] else default
}

summarize_seed_curve <- function(seed_data) {
  value_col <- if ("smoothed_dominant_mean_ploidy" %in% names(seed_data)) "smoothed_dominant_mean_ploidy" else if ("fitted_value" %in% names(seed_data)) "fitted_value" else "dominant_mean_ploidy"
  epsilon_cols <- grep("^(smooth_)?step_epsilon([.][xy])?$", names(seed_data), value = TRUE)
  x <- suppressWarnings(as.numeric(seed_data$O2_pct))
  y <- suppressWarnings(as.numeric(seed_data[[value_col]]))
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  order_index <- order(x); x <- x[order_index]; y <- y[order_index]
  seed_id <- as.character(seed_data$seed_id[[1L]])
  epsilon <- if (length(epsilon_cols)) combined_first_finite(seed_data[[tail(epsilon_cols, 1L)]], 0) else 0
  dx <- diff(x); dy <- diff(y); valid <- is.finite(dx) & dx > 0 & is.finite(dy)
  dx <- dx[valid]; dy <- dy[valid]; slopes <- dy / dx
  n_intervals <- length(slopes); o2_range <- if (length(x) > 1L) tail(x, 1L) - x[[1L]] else NA_real_
  safe_fraction <- function(condition) if (n_intervals) sum(condition, na.rm = TRUE) / n_intervals else NA_real_
  data.frame(
    seed_id = seed_id, seed_number = combined_seed_number(seed_id), n_o2_points = length(x),
    o2_min = if (length(x)) x[[1L]] else NA_real_, o2_max = if (length(x)) tail(x, 1L) else NA_real_,
    smoothed_ploidy_at_o2_min = if (length(y)) y[[1L]] else NA_real_, smoothed_ploidy_at_o2_max = if (length(y)) tail(y, 1L) else NA_real_,
    smoothed_ploidy_range = if (length(y)) diff(range(y)) else NA_real_, net_smoothed_ploidy_change = if (length(y) > 1L) tail(y, 1L) - y[[1L]] else NA_real_,
    average_slope = if (is.finite(o2_range) && o2_range > 0) (tail(y, 1L) - y[[1L]]) / o2_range else NA_real_,
    mean_interval_slope = if (n_intervals) mean(slopes) else NA_real_, weighted_mean_interval_slope = if (n_intervals) stats::weighted.mean(slopes, dx) else NA_real_,
    median_interval_slope = if (n_intervals) stats::median(slopes) else NA_real_, min_interval_slope = if (n_intervals) min(slopes) else NA_real_, max_interval_slope = if (n_intervals) max(slopes) else NA_real_,
    positive_interval_fraction = safe_fraction(dy > epsilon), negative_interval_fraction = safe_fraction(dy < -epsilon), near_flat_interval_fraction = safe_fraction(abs(dy) <= epsilon),
    n_intervals = n_intervals, smooth_step_epsilon = epsilon, stringsAsFactors = FALSE
  )
}

compute_average_slope_table <- function(curves, by_seed = data.frame()) {
  combined_require_columns(curves, c("seed_id", "O2_pct"), "regression curve table")
  if (!any(c("smoothed_dominant_mean_ploidy", "fitted_value", "dominant_mean_ploidy") %in% names(curves))) stop("Regression curve table has no fitted ploidy column.", call. = FALSE)
  rows <- lapply(split(curves, curves$seed_id), summarize_seed_curve)
  out <- do.call(rbind, rows); row.names(out) <- NULL
  if (nrow(by_seed) && "seed_id" %in% names(by_seed)) {
    keep <- intersect(c("seed_id", "curve_class", "final_interpretation_class", "smooth_curve_class", "pointwise_curve_class", "class_changed", "classification_rule_version", "smooth_classification_rule_version", "objective", "objective_source", "objective_total", "objective_data", "objective_burden", "objective_ploidy", "convergence_status", "monotonicity_reliability"), names(by_seed))
    out <- merge(out, by_seed[!duplicated(by_seed$seed_id), keep, drop = FALSE], by = "seed_id", all.x = TRUE, sort = FALSE)
  }
  out[order(out$seed_number, out$seed_id), , drop = FALSE]
}

write_average_slope_analysis <- function(curve_table, by_seed_table = NULL, out_dir, output_file = NULL, dry_run = FALSE) {
  output_file <- output_file %||% file.path(combined_analysis_dir(out_dir), "fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv")
  if (dry_run) return(invisible(output_file))
  curves <- combined_read_tsv(curve_table)
  by_seed <- if (!is.null(by_seed_table) && file.exists(by_seed_table)) combined_read_tsv(by_seed_table) else data.frame()
  result <- compute_average_slope_table(curves, by_seed)
  combined_write_tsv(result, output_file)
  combined_write_tsv(result, file.path(combined_legacy_table_dir(out_dir), basename(output_file)))
  invisible(output_file)
}

combined_annotation_maps <- function(by_seed, slopes, class_col) {
  if (!class_col %in% names(by_seed)) stop("Class table is missing --class_col=", class_col, call. = FALSE)
  if (!"seed_number" %in% names(by_seed)) by_seed$seed_number <- combined_seed_number(by_seed$seed_id)
  class_rows <- by_seed[is.finite(by_seed$seed_number), c("seed_number", class_col), drop = FALSE]
  class_rows <- class_rows[!duplicated(class_rows$seed_number), , drop = FALSE]
  slope_rows <- slopes[is.finite(slopes$seed_number), c("seed_number", "average_slope"), drop = FALSE]
  list(classes = stats::setNames(as.character(class_rows[[class_col]]), class_rows$seed_number), slopes = stats::setNames(as.numeric(slope_rows$average_slope), slope_rows$seed_number))
}

combined_annotate_coordinates <- function(data, maps) {
  combined_require_columns(data, c("dataset", "point_type", "seed"), "coordinate table")
  data$seed_number <- combined_seed_number(data$seed)
  data$curve_class <- NA_character_; data$average_slope <- NA_real_
  target <- data$dataset == "invivo" & data$point_type == "best" & is.finite(data$seed_number)
  data$curve_class[target] <- unname(maps$classes[as.character(data$seed_number[target])])
  data$average_slope[target] <- suppressWarnings(as.numeric(unname(maps$slopes[as.character(data$seed_number[target])])) )
  data
}

prepare_combined_parameter_landscape_tables <- function(argv = combined_parse_args()) {
  out_dir <- normalizePath(path.expand(argv$out_dir %||% combined_default_out_dir()), mustWork = FALSE)
  dense_root <- normalizePath(path.expand(argv$dense_grid_dir %||% combined_default_dense_root()), mustWork = TRUE)
  curve_table <- combined_find_latest(dense_root, "^fixed_o2_ploidy_monotonicity_regression_curves[.]tsv$", argv$curve_table)
  class_table <- combined_find_latest(dense_root, "^fixed_o2_ploidy_monotonicity_regression_by_seed[.]tsv$", argv$class_table %||% argv$by_seed_table)
  slope_file <- normalizePath(path.expand(argv$output_file %||% file.path(combined_analysis_dir(out_dir), "fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv")), mustWork = FALSE)
  dry_run <- combined_as_bool(argv$dry_run, FALSE)
  write_average_slope_analysis(curve_table, class_table, out_dir, slope_file, dry_run)
  if (combined_as_bool(argv$slope_only %||% argv$slopes_only, FALSE)) return(invisible(slope_file))
  pooled_candidate <- argv$pooled_root %||% if (!is.null(argv$parameter_root)) file.path(argv$parameter_root, "pooled_invivo_invitro") else combined_default_pooled_root()
  pooled_root <- normalizePath(path.expand(pooled_candidate), mustWork = TRUE)
  embeddings <- combined_discover_coordinates(pooled_root, combined_normalize_reductions(argv$reductions), combined_normalize_variants(argv$variants))
  if (!nrow(embeddings)) stop("No coordinate tables under ", pooled_root, call. = FALSE)
  if (dry_run) { print(embeddings); return(invisible(embeddings)) }
  by_seed <- combined_read_tsv(class_table); slopes <- combined_read_tsv(slope_file)
  maps <- combined_annotation_maps(by_seed, slopes, argv$class_col %||% "curve_class")
  rows <- vector("list", nrow(embeddings))
  for (i in seq_len(nrow(embeddings))) {
    item <- embeddings[i, , drop = FALSE]
    annotated <- combined_annotate_coordinates(combined_read_csv(item$coordinate_table), maps)
    rel_dir <- file.path(item$reduction, item$variant)
    annotated_path <- file.path(combined_analysis_dir(out_dir), "coordinates", rel_dir, paste0(item$stub, "_annotated.csv"))
    combined_write_csv(annotated, annotated_path)
    best <- annotated[annotated$point_type == "best", , drop = FALSE]
    combined_write_csv(best, file.path(combined_legacy_table_dir(out_dir), rel_dir, paste0(item$stub, "_best_points_curve_class.csv")))
    counts <- as.data.frame(table(dataset = best$dataset, curve_class = ifelse(is.na(best$curve_class), "not_applicable_or_missing", best$curve_class)), stringsAsFactors = FALSE)
    counts <- counts[counts$Freq > 0L, , drop = FALSE]
    combined_write_tsv(counts, file.path(combined_legacy_table_dir(out_dir), rel_dir, paste0(item$stub, "_curve_class_counts.tsv")))
    rows[[i]] <- data.frame(reduction = item$reduction, variant = item$variant, stub = item$stub, coordinate_table = item$coordinate_table, annotated_table = normalizePath(annotated_path, mustWork = FALSE), original_png = item$original_png, n_rows = nrow(annotated), n_initial = sum(annotated$point_type == "initial"), n_best_invivo = sum(annotated$dataset == "invivo" & annotated$point_type == "best"), n_best_invitro = sum(annotated$dataset == "invitro" & annotated$point_type == "best"), n_invivo_best_with_class = sum(annotated$dataset == "invivo" & annotated$point_type == "best" & !is.na(annotated$curve_class)), n_invivo_best_with_average_slope = sum(annotated$dataset == "invivo" & annotated$point_type == "best" & is.finite(annotated$average_slope)), class_table = class_table, average_slope_table = slope_file, stringsAsFactors = FALSE)
  }
  manifest <- do.call(rbind, rows)
  manifest_path <- file.path(combined_analysis_dir(out_dir), "combined_embedding_table_manifest.tsv")
  combined_write_tsv(manifest, manifest_path)
  message("Wrote combined analysis manifest: ", manifest_path)
  invisible(manifest)
}

main <- prepare_combined_parameter_landscape_tables
if (sys.nframe() == 0L) main()
