#!/usr/bin/env Rscript

# Render pair-level Cat x parameter and Cat x ratio-class descriptive figures.
# No within-pair Cat comparison is implied when a pair contains one category.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(SCRIPT_DIR, "..", "joint_soft_coupling_stability", "joint_coupling_visualization_utils.R"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(name, default = NULL) { hit <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default }
  analysis_dir <- get_arg("analysis_dir") %||% stop("--analysis_dir is required", call. = FALSE)
  out_dir <- get_arg("out_dir") %||% file.path(analysis_dir, "figures")
  summary_path <- file.path(analysis_dir, "cat_parameter_pair_balanced_values.tsv")
  pair_path <- file.path(analysis_dir, "pair_level_cat_parameter_summary.tsv")
  cells_path <- file.path(analysis_dir, "cat_ratio_class_pair_balanced_cells.tsv")
  association_path <- file.path(analysis_dir, "cat_ratio_class_pair_level_association.tsv")
  metric_path <- file.path(analysis_dir, "cat_parameter_metric_information_comparison.tsv")
  summary <- o2jcv_read_tsv(summary_path); pair_data <- o2jcv_read_tsv(pair_path); cells <- o2jcv_read_tsv(cells_path); association <- o2jcv_read_tsv(association_path); metrics <- o2jcv_read_tsv(metric_path)
  parameter_levels <- rev(o2jcv_parameter_levels())
  threshold <- unique(pair_data$class_threshold)[[1L]]
  lower_bound <- if ("class_lower_bound" %in% names(pair_data)) unique(pair_data$class_lower_bound)[[1L]] else 1 / threshold
  upper_bound <- if ("class_upper_bound" %in% names(pair_data)) unique(pair_data$class_upper_bound)[[1L]] else threshold

  summary$parameter <- factor(summary$parameter, levels = parameter_levels)
  summary$ploidy_category <- factor(summary$ploidy_category, levels = c("CatA", "CatB", "CatC"))
  p1 <- ggplot2::ggplot(summary, ggplot2::aes(x = pair_balanced_mean_log2_ratio, y = parameter, colour = ploidy_category)) +
    o2jcv_class_band(lower_bound, upper_bound) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = pair_min_log2_ratio, xmax = pair_max_log2_ratio), orientation = "y", width = 0, colour = "#8A9299", linewidth = 0.55) +
    ggplot2::geom_point(size = 2.8) + ggplot2::geom_vline(xintercept = 0, colour = "#343A40", linewidth = 0.35) +
    ggplot2::facet_wrap(~ploidy_category, nrow = 1) + ggplot2::scale_colour_manual(values = o2jcv_category_colors(), guide = "none") +
    ggplot2::labs(x = "Pair-balanced mean log2(in vivo / in vitro), with pair range", y = NULL, title = "Soft-coupling ratios within each ploidy category", subtitle = "Two pairs contribute to each Cat; gold band is ClassB. This is a between-pair descriptive comparison") + o2jcv_theme(8.5, rotate_x = FALSE)
  o2jcv_save_figure(p1, "cat_parameter_pair_level_ratio_intervals", out_dir, c(summary_path), 13.5, 7, family = "Uncertainty & Benchmark", question = "Which parameter ratios are shared by the two pairs assigned to each ploidy category?")

  cells$parameter <- factor(cells$parameter, levels = parameter_levels)
  cells$cell <- paste(cells$ploidy_category, cells$ratio_class, sep = " × ")
  cells$label_colour <- ifelse(cells$pair_balanced_mean_proportion >= 0.65, "white", "#26313A")
  p2 <- ggplot2::ggplot(cells, ggplot2::aes(x = cell, y = parameter, fill = pair_balanced_mean_proportion)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * pair_balanced_mean_proportion), colour = label_colour), size = 2.6) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient(low = "#FFF9E8", high = "#9B7418", limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Ploidy category × ratio class", y = NULL, fill = "Pair-balanced\nseed proportion", title = "CatA/B/C × ClassA/B/C profile", subtitle = "Percentages average available pairs within each Cat and avoid treating seeds as independent pairs") + o2jcv_theme(8)
  o2jcv_save_figure(p2, "cat_ratio_class_pair_balanced_profile", out_dir, c(cells_path), 13, 7, family = "Matrix & Cohort", question = "How does the ClassA/B/C profile change across CatA/B/C pairs?")

  association$parameter <- factor(association$parameter, levels = parameter_levels)
  association_plot <- association[is.finite(association$pair_level_cramers_v), , drop = FALSE]
  p3 <- ggplot2::ggplot(association_plot, ggplot2::aes(x = pair_level_cramers_v, y = parameter)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = pair_level_cramers_v, yend = parameter), colour = "#D8DDE2", linewidth = 1.1) +
    ggplot2::geom_point(shape = 21, fill = "#D4A72C", colour = "#6F5417", size = 3) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(x = "Descriptive pair-level Cramér's V", y = NULL, title = "Pair-level ploidy-category × dominant-class pattern", subtitle = "Only six pairs are available; values are descriptive and cannot separate category effects from pair effects") + o2jcv_theme(9.5, rotate_x = FALSE)
  o2jcv_save_figure(p3, "cat_ratio_class_pair_level_effect_size", out_dir, c(association_path), 9, 6, family = "Relationship", question = "At the pair level, how strongly does ploidy category align with a parameter's dominant ratio class?")

  dumbbell <- rbind(
    data.frame(summary[, c("ploidy_category", "parameter")], source = "In vitro", value = summary$pair_balanced_mean_vitro),
    data.frame(summary[, c("ploidy_category", "parameter")], source = "In vivo", value = summary$pair_balanced_mean_vivo)
  )
  dumbbell$ploidy_category <- factor(dumbbell$ploidy_category, levels = c("CatA", "CatB", "CatC"))
  segment <- summary
  p4 <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = segment, ggplot2::aes(x = pair_balanced_mean_vitro, xend = pair_balanced_mean_vivo, y = ploidy_category, yend = ploidy_category, colour = ploidy_category), linewidth = 0.8) +
    ggplot2::geom_point(data = dumbbell, ggplot2::aes(x = value, y = ploidy_category, shape = source, fill = ploidy_category), size = 2.8, colour = "#26313A", stroke = 0.5) +
    ggplot2::facet_wrap(~parameter, scales = "free_x", ncol = 4) +
    ggplot2::scale_x_log10() + ggplot2::scale_colour_manual(values = o2jcv_category_colors(), guide = "none") +
    ggplot2::scale_fill_manual(values = o2jcv_category_colors(), guide = "none") +
    ggplot2::scale_shape_manual(values = c("In vitro" = 21, "In vivo" = 24)) +
    ggplot2::labs(x = "Pair-balanced natural parameter value (log10 scale)", y = NULL, shape = NULL, title = "In-vivo versus in-vitro parameter values by ploidy category", subtitle = "Dumbbells show whether Cat differences arise from the in-vivo value, the shared in-vitro anchor, or both") + o2jcv_theme(7.5)
  o2jcv_save_figure(p4, "cat_parameter_vivo_vitro_dumbbells", out_dir, c(summary_path), 13, 11, family = "Comparison & Ranking", question = "Which side of the joint fit drives category-specific parameter differences?")

  pair_data$parameter <- factor(pair_data$parameter, levels = parameter_levels)
  pair_data$ploidy_category <- factor(pair_data$ploidy_category, levels = c("CatA", "CatB", "CatC"))
  p5 <- ggplot2::ggplot(pair_data, ggplot2::aes(x = median_log2_ratio, y = parameter, shape = pair_label)) +
    o2jcv_class_band(lower_bound, upper_bound) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = q05_log2_ratio, xmax = q95_log2_ratio), orientation = "y", width = 0, colour = "#A0A7AD", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(fill = ploidy_category), size = 2.5, colour = "#26313A") +
    ggplot2::geom_vline(xintercept = 0, colour = "#343A40", linewidth = 0.3) +
    ggplot2::facet_wrap(~ploidy_category, nrow = 1) +
    ggplot2::scale_fill_manual(values = o2jcv_category_colors(), guide = "none") +
    ggplot2::labs(x = "Pair median log2(in vivo / in vitro), with 5th–95th seed percentiles", y = NULL, shape = "Pair", title = "The two pairs underlying each ploidy category", subtitle = "Pair-level points prevent category averages from hiding disagreement between the two contributing pairs") + o2jcv_theme(8, rotate_x = FALSE)
  o2jcv_save_figure(p5, "cat_parameter_pair_level_detail", out_dir, c(pair_path), 14, 7, family = "Uncertainty & Benchmark", question = "Do the two pairs within each ploidy category agree on parameter direction and spread?")
}

if (identical(environment(), globalenv())) main()
