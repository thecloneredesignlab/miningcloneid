#!/usr/bin/env Rscript

.dense_analysis_dir <- local({
  if (exists(".dense_analysis_source_override", inherits = TRUE)) return(dirname(normalizePath(get(".dense_analysis_source_override", inherits = TRUE), mustWork = TRUE)))
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "dense_grid_analysis_utils.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.dense_analysis_root <- normalizePath(file.path(.dense_analysis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.dense_analysis_root, "util", "o2_supply_demand_map_dense_grid_utils.R"), local = environment(), chdir = TRUE)
source(file.path(.dense_analysis_root, "util", "o2_supply_demand_map_curve_classification_utils.R"), local = environment(), chdir = TRUE)

dense_analysis_materialized_path <- function(out_dir, part) file.path(dense_grid_simulation_dir(out_dir), if (dense_grid_normalize_part(part) == "monotonicity") "fixed_o2_dense_grid_attractor_curves.tsv" else "fixed_o2_initial_ploidy_trajectories.tsv")

dense_classify_one_curve <- function(curve,
                                     value_col = "dominant_mean_ploidy",
                                     smooth = FALSE,
                                     flat_range_threshold = 0.05,
                                     step_epsilon_abs = 1e-6,
                                     step_epsilon_fraction = 1e-4,
                                     reverse_fraction_tolerance = 0.05,
                                     plateau_min_points = 3L,
                                     smooth_span = 0.20) {
  classifier <- if (smooth) classify_o2_ploidy_curve_smooth else classify_o2_ploidy_curve
  args <- list(
    curve = curve,
    value_col = value_col,
    x_col = "O2_pct",
    id_col = "seed_id",
    flat_range_threshold = flat_range_threshold,
    step_epsilon_abs = step_epsilon_abs,
    step_epsilon_fraction = step_epsilon_fraction,
    reverse_fraction_tolerance = reverse_fraction_tolerance
  )
  if (smooth) args$smooth_span <- smooth_span else args$plateau_min_points <- plateau_min_points
  result <- do.call(classifier, args)
  result$summary$seed_id <- as.character(curve$seed_id[[1L]])
  result$summary$seed_number <- dense_seed_number(curve$seed_id[[1L]])
  result
}

dense_curve_reliability <- function(curve, gap_low = 0.005, gap_caution = 0.01, unreliable_fraction = 0.25, caution_fraction = 0.25) {
  gap <- suppressWarnings(as.numeric(curve$spectral_gap))
  low <- mean(gap < gap_low, na.rm = TRUE)
  caution <- mean(gap < gap_caution, na.rm = TRUE)
  if (is.finite(low) && low >= unreliable_fraction) "unreliable_small_gap" else if ((is.finite(low) && low > 0) || (is.finite(caution) && caution >= caution_fraction)) "caution_small_gap" else "reliable"
}

dense_classify_attractor_curves <- function(curves, argv, smooth = FALSE) {
  reporting_o2 <- sort(unique(dense_as_num_vec(argv$reporting_o2, c(0, 0.1, 0.5, 1, 2, 5))))
  groups <- split(curves, curves$seed_id)
  results <- lapply(groups, function(curve) dense_classify_one_curve(
    curve,
    smooth = smooth,
    flat_range_threshold = dense_as_num(argv$flat_range_threshold, 0.05),
    step_epsilon_abs = dense_as_num(argv$step_epsilon_abs, 1e-6),
    step_epsilon_fraction = dense_as_num(argv$step_epsilon_fraction, 1e-4),
    reverse_fraction_tolerance = dense_as_num(argv$reverse_fraction_tolerance, 0.05),
    plateau_min_points = dense_as_int(argv$plateau_min_points, 3L),
    smooth_span = dense_as_num(argv$smooth_span, 0.20)
  ))
  summaries <- dense_rbind_fill(lapply(results, `[[`, "summary"))
  differences <- dense_rbind_fill(lapply(results, `[[`, "differences"))
  rows <- lapply(seq_along(groups), function(i) {
    curve <- groups[[i]]
    summary <- results[[i]]$summary
    y <- suppressWarnings(as.numeric(curve$dominant_mean_ploidy))
    gap <- suppressWarnings(as.numeric(curve$spectral_gap))
    reliability <- dense_curve_reliability(curve)
    base <- data.frame(
      seed_id = curve$seed_id[[1L]], seed_number = dense_seed_number(curve$seed_id[[1L]]),
      n_o2 = nrow(curve), o2_min = min(curve$O2_pct, na.rm = TRUE), o2_max = max(curve$O2_pct, na.rm = TRUE),
      curve_class = summary$curve_class[[1L]],
      final_interpretation_class = if (reliability == "unreliable_small_gap") "unreliable_small_spectral_gap" else summary$curve_class[[1L]],
      monotonicity_reliability = reliability,
      min_ploidy = min(y, na.rm = TRUE), max_ploidy = max(y, na.rm = TRUE),
      min_spectral_gap = suppressWarnings(min(gap, na.rm = TRUE)), median_spectral_gap = stats::median(gap, na.rm = TRUE),
      fraction_o2_gap_below_0p005 = mean(gap < 0.005, na.rm = TRUE), fraction_o2_gap_below_0p01 = mean(gap < 0.01, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    meta <- summary[, setdiff(names(summary), c("seed_id", "seed_number", "curve_class")), drop = FALSE]
    out <- cbind(base, meta)
    for (o2 in reporting_o2) {
      hit <- which(abs(as.numeric(curve$O2_pct) - o2) < 1e-9)
      slug <- gsub("[.]", "p", format(o2, scientific = FALSE, trim = TRUE))
      out[[paste0("mode_label_o2_", slug)]] <- if (length(hit) && "mode_label" %in% names(curve)) curve$mode_label[[hit[[1L]]]] else NA_character_
      out[[paste0("dominant_mean_ploidy_o2_", slug)]] <- if (length(hit)) y[[hit[[1L]]]] else NA_real_
      out[[paste0("spectral_gap_o2_", slug)]] <- if (length(hit)) gap[[hit[[1L]]]] else NA_real_
    }
    for (column in intersect(c("objective", "delta_objective"), names(curve))) out[[column]] <- curve[[column]][[1L]]
    out
  })
  by_seed <- dense_rbind_fill(rows)
  difference_keep <- intersect(c("seed_id", "O2_pct", "finite_difference_next", "local_slope_sign", "fitted_value", "fitted_difference_next", "fitted_local_slope_sign", "step_epsilon"), names(differences))
  enriched <- merge(curves, differences[, difference_keep, drop = FALSE], by = c("seed_id", "O2_pct"), all.x = TRUE, sort = FALSE)
  enriched <- enriched[order(dense_seed_number(enriched$seed_id), as.numeric(enriched$O2_pct)), , drop = FALSE]
  list(by_seed = by_seed, curves = enriched, raw_summaries = summaries, differences = differences)
}

dense_curve_class_counts <- function(by_seed) {
  out <- as.data.frame(table(by_seed$curve_class), stringsAsFactors = FALSE)
  names(out) <- c("curve_class", "n_seed")
  out <- out[out$n_seed > 0L, , drop = FALSE]
  out$fraction_seed <- out$n_seed / sum(out$n_seed)
  out[order(-out$n_seed, out$curve_class), , drop = FALSE]
}

dense_curve_summary <- function(curves, by_seed, value_col = "dominant_mean_ploidy") {
  joined <- merge(curves, by_seed[, c("seed_id", "curve_class"), drop = FALSE], by = "seed_id", all.x = TRUE, sort = FALSE)
  groups <- split(joined, interaction(joined$curve_class, joined$O2_pct, drop = TRUE))
  dense_rbind_fill(lapply(groups, function(data) {
    values <- suppressWarnings(as.numeric(data[[value_col]]))
    data.frame(curve_class = data$curve_class[[1L]], O2_pct = data$O2_pct[[1L]], n_seed = length(unique(data$seed_id)), median_dominant_mean_ploidy = stats::median(values, na.rm = TRUE), q25_dominant_mean_ploidy = as.numeric(stats::quantile(values, 0.25, na.rm = TRUE, names = FALSE)), q75_dominant_mean_ploidy = as.numeric(stats::quantile(values, 0.75, na.rm = TRUE, names = FALSE)), median_spectral_gap = if ("spectral_gap" %in% names(data)) stats::median(data$spectral_gap, na.rm = TRUE) else NA_real_, stringsAsFactors = FALSE)
  }))
}

dense_representative_seeds <- function(curves, by_seed) {
  dense_rbind_fill(lapply(split(by_seed, by_seed$curve_class), function(seed_table) {
    subset <- curves[curves$seed_id %in% seed_table$seed_id, , drop = FALSE]
    medians <- aggregate(subset$dominant_mean_ploidy, list(O2_pct = subset$O2_pct), stats::median, na.rm = TRUE)
    names(medians)[[2L]] <- "class_median"
    score <- dense_rbind_fill(lapply(split(subset, subset$seed_id), function(curve) {
      joined <- merge(curve[, c("seed_id", "O2_pct", "dominant_mean_ploidy")], medians, by = "O2_pct")
      data.frame(seed_id = curve$seed_id[[1L]], median_curve_rmse = sqrt(mean((joined$dominant_mean_ploidy - joined$class_median)^2, na.rm = TRUE)))
    }))
    hit <- score[which.min(score$median_curve_rmse), , drop = FALSE]
    data.frame(curve_class = seed_table$curve_class[[1L]], representative_seed_id = hit$seed_id, representative_seed_number = dense_seed_number(hit$seed_id), median_curve_rmse = hit$median_curve_rmse, class_n_seed = nrow(seed_table), stringsAsFactors = FALSE)
  }))
}

dense_objective_statistics <- function(by_seed) {
  data <- by_seed[is.finite(suppressWarnings(as.numeric(by_seed$objective))), c("seed_id", "curve_class", "objective"), drop = FALSE]
  if (!nrow(data)) return(list(plot_data = data.frame(), summary = data.frame(), global = data.frame(), pairwise = data.frame()))
  data$objective <- as.numeric(data$objective)
  summary <- dense_rbind_fill(lapply(split(data, data$curve_class), function(group) data.frame(curve_class = group$curve_class[[1L]], n_seed = nrow(group), median_objective = stats::median(group$objective), q25_objective = as.numeric(stats::quantile(group$objective, 0.25)), q75_objective = as.numeric(stats::quantile(group$objective, 0.75)), stringsAsFactors = FALSE)))
  global <- if (length(unique(data$curve_class)) > 1L) {
    test <- stats::kruskal.test(objective ~ curve_class, data = data)
    data.frame(test = "Kruskal-Wallis", statistic = unname(test$statistic), p_value = test$p.value, stringsAsFactors = FALSE)
  } else data.frame(test = "Kruskal-Wallis", statistic = NA_real_, p_value = NA_real_)
  classes <- unique(data$curve_class)
  pairs <- if (length(classes) > 1L) utils::combn(classes, 2L, simplify = FALSE) else list()
  pairwise <- dense_rbind_fill(lapply(pairs, function(pair) {
    x <- data$objective[data$curve_class == pair[[1L]]]; y <- data$objective[data$curve_class == pair[[2L]]]
    test <- tryCatch(stats::wilcox.test(x, y, exact = FALSE), error = function(e) NULL)
    data.frame(class_a = pair[[1L]], class_b = pair[[2L]], p_value = if (is.null(test)) NA_real_ else test$p.value, stringsAsFactors = FALSE)
  }))
  if (nrow(pairwise)) pairwise$fdr <- stats::p.adjust(pairwise$p_value, "BH")
  list(plot_data = data, summary = summary, global = global, pairwise = pairwise)
}

dense_write_analysis_table <- function(data, canonical_dir, legacy_dir, filename, artifact_id, out_dir, source) {
  canonical <- file.path(canonical_dir, filename)
  legacy <- file.path(legacy_dir, filename)
  dense_write_tsv(data, canonical)
  dense_write_tsv(data, legacy)
  dense_record_artifact(out_dir, "analysis", artifact_id, canonical, data, "analyze_dense_grid_tables.R", source)
  canonical
}

dense_materialize_monotonicity_analysis <- function(out_dir, argv) {
  source_path <- dense_analysis_materialized_path(out_dir, "monotonicity")
  if (!file.exists(source_path)) stop("Missing dense-grid simulation table: ", source_path, call. = FALSE)
  curves <- dense_read_tsv(source_path)
  dense_require_columns(curves, c("seed_id", "O2_pct", "dominant_mean_ploidy", "spectral_gap"), source_path)
  pointwise <- dense_classify_attractor_curves(curves, argv, FALSE)
  smooth <- dense_classify_attractor_curves(curves, argv, TRUE)
  # Preserve the established downstream table contract while keeping the
  # canonical classifier's neutral fitted-value naming internally.
  smooth$curves$smoothed_dominant_mean_ploidy <- smooth$curves$fitted_value
  smooth_step_col <- grep("^step_epsilon([.]y)?$", names(smooth$curves), value = TRUE)
  smooth$curves$smooth_step_epsilon <- if (length(smooth_step_col)) smooth$curves[[tail(smooth_step_col, 1L)]] else NA_real_
  smooth$by_seed$smooth_curve_class <- smooth$by_seed$curve_class
  smooth$by_seed$smooth_classification_rule_version <- smooth$by_seed$classification_rule_version
  smooth$by_seed$pointwise_curve_class <- pointwise$by_seed$curve_class[match(smooth$by_seed$seed_id, pointwise$by_seed$seed_id)]
  smooth$by_seed$class_changed <- smooth$by_seed$curve_class != smooth$by_seed$pointwise_curve_class
  canonical_dir <- dense_grid_analysis_dir(out_dir); legacy_dir <- dense_grid_legacy_table_dir(out_dir)
  dense_write_analysis_table(pointwise$by_seed, canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_by_seed.tsv", "dense_grid_monotonicity_by_seed", out_dir, source_path)
  dense_write_analysis_table(pointwise$curves, canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_curves.tsv", "dense_grid_monotonicity_curves", out_dir, source_path)
  dense_write_analysis_table(dense_curve_class_counts(pointwise$by_seed), canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_class_counts.tsv", "dense_grid_monotonicity_class_counts", out_dir, source_path)
  dense_write_analysis_table(dense_curve_summary(pointwise$curves, pointwise$by_seed), canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_class_curve_summary.tsv", "dense_grid_monotonicity_class_summary", out_dir, source_path)
  dense_write_analysis_table(dense_representative_seeds(pointwise$curves, pointwise$by_seed), canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_representative_seeds.tsv", "dense_grid_monotonicity_representatives", out_dir, source_path)
  dense_write_analysis_table(smooth$by_seed, canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv", "dense_grid_smooth_monotonicity_by_seed", out_dir, source_path)
  dense_write_analysis_table(smooth$curves, canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_regression_curves.tsv", "dense_grid_smooth_monotonicity_curves", out_dir, source_path)
  dense_write_analysis_table(dense_curve_class_counts(smooth$by_seed), canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_regression_class_counts.tsv", "dense_grid_smooth_monotonicity_counts", out_dir, source_path)
  dense_write_analysis_table(dense_curve_summary(smooth$curves, smooth$by_seed), canonical_dir, legacy_dir, "fixed_o2_ploidy_monotonicity_regression_class_curve_summary.tsv", "dense_grid_smooth_monotonicity_summary", out_dir, source_path)
  objective <- dense_objective_statistics(smooth$by_seed)
  dense_write_analysis_table(objective$plot_data, canonical_dir, legacy_dir, "fixed_o2_regression_objective_by_curve_class_boxplot_plot_data.tsv", "dense_grid_objective_plot_data", out_dir, source_path)
  dense_write_analysis_table(objective$summary, canonical_dir, legacy_dir, "fixed_o2_regression_objective_by_curve_class_boxplot_summary.tsv", "dense_grid_objective_summary", out_dir, source_path)
  dense_write_analysis_table(objective$global, canonical_dir, legacy_dir, "fixed_o2_regression_objective_by_curve_class_boxplot_global_test.tsv", "dense_grid_objective_global_test", out_dir, source_path)
  dense_write_analysis_table(objective$pairwise, canonical_dir, legacy_dir, "fixed_o2_regression_objective_by_curve_class_boxplot_pairwise_tests.tsv", "dense_grid_objective_pairwise", out_dir, source_path)
  invisible(pointwise)
}

dense_initial_selected_table <- function(raw) {
  selected <- raw
  rename <- c(analytical_mean_N = "mean_N", analytical_mean_ploidy = "mean_ploidy", analytical_fraction_N_le_25 = "fraction_N_le_25", analytical_fraction_N_below_44 = "fraction_N_below_44", analytical_fraction_N_ge_44 = "fraction_N_ge_44", analytical_fraction_N_ge_66 = "fraction_N_ge_66", analytical_fraction_N_ge_88 = "fraction_N_ge_88")
  for (old in names(rename)) if (old %in% names(selected)) names(selected)[names(selected) == old] <- rename[[old]]
  selected$seed_number <- dense_seed_number(selected$seed_id)
  selected
}

dense_initial_delta <- function(selected) {
  keys <- c("seed_id", "seed_number", "O2_pct", "day")
  values <- intersect(c("mean_N", "mean_ploidy", "fraction_N_le_25", "fraction_N_below_44", "fraction_N_ge_44", "fraction_N_ge_66", "fraction_N_ge_88", "used_initial_N", "initial_ploidy"), names(selected))
  two <- selected[abs(as.numeric(selected$initial_ploidy) - 2) < 1e-9, c(keys, values), drop = FALSE]
  four <- selected[abs(as.numeric(selected$initial_ploidy) - 4) < 1e-9, c(keys, values), drop = FALSE]
  names(two)[!names(two) %in% keys] <- paste0(names(two)[!names(two) %in% keys], "_2N")
  names(four)[!names(four) %in% keys] <- paste0(names(four)[!names(four) %in% keys], "_4N")
  out <- merge(two, four, by = keys, all = TRUE, sort = FALSE)
  out$delta_initial <- out$mean_ploidy_4N - out$mean_ploidy_2N
  out$abs_delta_initial <- abs(out$delta_initial)
  dominant <- unique(selected[, intersect(c(keys, "dominant_mean_ploidy", "spectral_gap"), names(selected)), drop = FALSE])
  out <- merge(out, dominant, by = keys, all.x = TRUE, sort = FALSE)
  out$abs_distance_to_dominant_ploidy_2N <- abs(out$mean_ploidy_2N - out$dominant_mean_ploidy)
  out$abs_distance_to_dominant_ploidy_4N <- abs(out$mean_ploidy_4N - out$dominant_mean_ploidy)
  out[order(out$seed_number, out$day, out$O2_pct), , drop = FALSE]
}

dense_initial_curve_classes <- function(selected, argv) {
  groups <- split(selected, interaction(selected$seed_id, selected$day, selected$initial_ploidy, drop = TRUE))
  dense_rbind_fill(lapply(groups, function(curve) {
    result <- dense_classify_one_curve(curve, value_col = "mean_ploidy", flat_range_threshold = dense_as_num(argv$flat_range_threshold, 0.05))
    cbind(data.frame(day = curve$day[[1L]], initial_condition = curve$initial_condition[[1L]], initial_ploidy = curve$initial_ploidy[[1L]], stringsAsFactors = FALSE), result$summary)
  }))
}

dense_initial_convergence <- function(delta) {
  terminal <- max(delta$day, na.rm = TRUE)
  data <- delta[abs(delta$day - terminal) < 1e-9, , drop = FALSE]
  data.frame(scope = "global", group = "all_seed_o2", terminal_day = terminal, n_seed_o2 = nrow(data), n_seed = length(unique(data$seed_id)), n_o2 = length(unique(data$O2_pct)), spectral_gap_median = if ("spectral_gap" %in% names(data)) stats::median(data$spectral_gap, na.rm = TRUE) else NA_real_, spectral_gap_min = if ("spectral_gap" %in% names(data)) min(data$spectral_gap, na.rm = TRUE) else NA_real_, abs_delta_initial_median = stats::median(data$abs_delta_initial, na.rm = TRUE), abs_delta_initial_q95 = as.numeric(stats::quantile(data$abs_delta_initial, 0.95, na.rm = TRUE)), abs_delta_initial_max = max(data$abs_delta_initial, na.rm = TRUE), fraction_abs_delta_le_0p01 = mean(data$abs_delta_initial <= 0.01, na.rm = TRUE), fraction_abs_delta_le_0p05 = mean(data$abs_delta_initial <= 0.05, na.rm = TRUE), stringsAsFactors = FALSE)
}

dense_materialize_initial_analysis <- function(out_dir, argv) {
  source_path <- dense_analysis_materialized_path(out_dir, "initial_ploidy")
  if (!file.exists(source_path)) stop("Missing dense-grid initial-ploidy simulation table: ", source_path, call. = FALSE)
  selected <- dense_initial_selected_table(dense_read_tsv(source_path))
  dense_require_columns(selected, c("seed_id", "day", "initial_ploidy", "O2_pct", "mean_ploidy"), source_path)
  delta <- dense_initial_delta(selected)
  classes <- dense_initial_curve_classes(selected, argv)
  convergence <- dense_initial_convergence(delta)
  canonical_dir <- dense_grid_analysis_dir(out_dir); legacy_dir <- dense_grid_legacy_table_dir(out_dir)
  dense_write_analysis_table(selected, canonical_dir, legacy_dir, "fixed_o2_initial_ploidy_selected_time_curves.tsv", "dense_grid_initial_selected_curves", out_dir, source_path)
  dense_write_analysis_table(delta, canonical_dir, legacy_dir, "fixed_o2_initial_ploidy_delta_by_seed_o2_time.tsv", "dense_grid_initial_delta", out_dir, source_path)
  dense_write_analysis_table(classes, canonical_dir, legacy_dir, "fixed_o2_initial_ploidy_curve_class_by_seed_time.tsv", "dense_grid_initial_classes", out_dir, source_path)
  dense_write_analysis_table(convergence, canonical_dir, legacy_dir, "fixed_o2_initial_ploidy_convergence_summary.tsv", "dense_grid_initial_convergence", out_dir, source_path)
  invisible(list(selected = selected, delta = delta, classes = classes, convergence = convergence))
}
