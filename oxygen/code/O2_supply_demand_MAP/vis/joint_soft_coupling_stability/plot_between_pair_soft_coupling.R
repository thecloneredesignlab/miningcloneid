#!/usr/bin/env Rscript

# Render pair-level ratio forests, class rankings, stable cores, and robustness
# views. All summaries are precomputed by the analysis layer.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(SCRIPT_DIR, "joint_coupling_visualization_utils.R"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(name, default = NULL) { hit <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default }
  analysis_dir <- get_arg("analysis_dir") %||% stop("--analysis_dir is required", call. = FALSE)
  out_dir <- get_arg("out_dir") %||% file.path(analysis_dir, "figures")
  within_path <- file.path(analysis_dir, "within_pair_parameter_stability.tsv")
  between_path <- file.path(analysis_dir, "between_pair_parameter_stability.tsv")
  class_path <- file.path(analysis_dir, "between_pair_class_summary.tsv")
  sensitivity_path <- file.path(analysis_dir, "threshold_sensitivity_pair_balanced.tsv")
  objective_path <- file.path(analysis_dir, "objective_sensitivity_pair_balanced.tsv")
  boundary_path <- file.path(analysis_dir, "within_pair_boundary_projection_summary.tsv")
  membership_path <- file.path(analysis_dir, "class_pair_set_membership.tsv")
  config_path <- file.path(analysis_dir, "analysis_config.tsv")
  appendix_dir <- get_arg("appendix_dir") %||% out_dir
  within <- o2jcv_add_pair_label(o2jcv_read_tsv(within_path))
  sensitivity <- o2jcv_read_tsv(sensitivity_path)
  between <- o2jcv_read_tsv(between_path)
  classes <- o2jcv_read_tsv(class_path)
  objective <- o2jcv_read_tsv(objective_path)
  boundary <- o2jcv_add_pair_label(o2jcv_read_tsv(boundary_path))
  membership <- o2jcv_add_pair_label(o2jcv_read_tsv(membership_path))
  config <- o2jcv_read_tsv(config_path)
  parameter_levels <- rev(o2jcv_parameter_levels())
  within$parameter <- factor(within$parameter, levels = parameter_levels)
  threshold <- suppressWarnings(as.numeric(config$value[config$key == "class_threshold"][[1L]]))
  p1 <- ggplot2::ggplot(within, ggplot2::aes(x = log2_ratio_median, y = parameter)) +
    o2jcv_class_band(threshold) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = log2(ratio_q05), xmax = log2(ratio_q95)), orientation = "y", width = 0, colour = "#71808C", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(fill = dominant_class), shape = 21, size = 2.2, colour = "#26313A", stroke = 0.35) +
    ggplot2::geom_vline(xintercept = 0, colour = "#343A40", linewidth = 0.35) +
    ggplot2::facet_wrap(~pair_label, ncol = 2) +
    ggplot2::scale_fill_manual(values = o2jcv_class_colors(), drop = FALSE) +
    ggplot2::labs(x = "Median log2(in vivo / in vitro), with 5th–95th seed percentiles", y = NULL, fill = "Dominant class", title = "Pair-specific parameter ratio intervals", subtitle = "Gold band marks ClassB; intervals expose both direction and seed-level spread") + o2jcv_theme(8.5, rotate_x = FALSE)
  o2jcv_save_figure(p1, "between_pair_log2_ratio_forest", out_dir, c(within_path, config_path), 11.5, 10.5, family = "Uncertainty & Benchmark", question = "Which parameters remain inside or outside the coupling band across seeds and pairs?")

  between$parameter <- factor(between$parameter, levels = rev(unique(between$parameter[order(between$cross_pair_dominant_fraction, decreasing = TRUE)])))
  between$pair_count_label <- sprintf("%.0f/%.0f pairs", between$cross_pair_dominant_fraction * between$n_pairs, between$n_pairs)
  p2 <- ggplot2::ggplot(between, ggplot2::aes(x = cross_pair_dominant_fraction, y = parameter, colour = cross_pair_dominant_class)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = cross_pair_dominant_fraction, yend = parameter), colour = "#D8DDE2", linewidth = 1.4) +
    ggplot2::geom_point(size = 3.4) +
    ggplot2::geom_text(ggplot2::aes(label = pair_count_label), hjust = -0.12, size = 3, colour = "#26313A") +
    ggplot2::geom_vline(xintercept = c(0.8, 0.9, 0.95), linetype = c(3, 2, 3), colour = "#8A9299", linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = o2jcv_class_colors(), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1.17), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Pairs sharing the cross-pair dominant class", y = NULL, colour = "Dominant class", title = "Between-pair class consensus", subtitle = "Pair is the statistical unit; labels show the exact number out of six pairs") + o2jcv_theme(10, rotate_x = FALSE)
  o2jcv_save_figure(p2, "between_pair_class_consensus", out_dir, c(between_path), 9.5, 6.5, family = "Comparison & Ranking", question = "Which parameters retain the same dominant class across pairs?")

  classes$parameter <- factor(classes$parameter, levels = parameter_levels)
  classes$ratio_class <- factor(classes$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  p3 <- ggplot2::ggplot(classes, ggplot2::aes(x = pair_balanced_mean_proportion, y = parameter, colour = ratio_class)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = bootstrap_ci_lower, xmax = bootstrap_ci_upper), orientation = "y", width = 0, colour = "#89939C", linewidth = 0.45) +
    ggplot2::geom_point(size = 2.8) + ggplot2::facet_wrap(~ratio_class, nrow = 1) +
    ggplot2::scale_colour_manual(values = o2jcv_class_colors(), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Pair-balanced mean seed proportion (95% pair-bootstrap interval)", y = NULL, title = "ClassA, ClassB, and ClassC parameter rankings", subtitle = "Each pair contributes equally, regardless of its number of valid seeds") + o2jcv_theme(8.5, rotate_x = FALSE)
  o2jcv_save_figure(p3, "between_pair_all_class_rankings", out_dir, c(class_path), 13.5, 7, family = "Comparison & Ranking", question = "Which parameters most consistently belong to each class across pairs?")

  stable <- classes[, c("parameter", "ratio_class", "n_pairs", "n_pairs_stable80", "n_pairs_stable90", "n_pairs_stable95", "n_pairs_strict_intersection")]
  stable$parameter <- factor(stable$parameter, levels = parameter_levels)
  stable$ratio_class <- factor(stable$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  n_pair_total <- max(stable$n_pairs, na.rm = TRUE)
  p4 <- ggplot2::ggplot(stable, ggplot2::aes(x = n_pairs_stable90, y = parameter, colour = ratio_class)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = n_pairs_stable90, yend = parameter), colour = "#D8DDE2", linewidth = 1.1) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(n_pairs_stable90, "/", n_pairs)), hjust = -0.25, size = 2.7, colour = "#26313A") +
    ggplot2::facet_wrap(~ratio_class, nrow = 1) +
    ggplot2::scale_colour_manual(values = o2jcv_class_colors(), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, n_pair_total + 0.8), breaks = 0:n_pair_total) +
    ggplot2::labs(x = "Pairs with ≥90% of seeds in the class", y = NULL, title = "Cross-pair stable-90 core coverage", subtitle = sprintf("%d out of %d indicates a class-stable parameter across every available pair", n_pair_total, n_pair_total)) + o2jcv_theme(8.5, rotate_x = FALSE)
  o2jcv_save_figure(p4, "between_pair_stable90_core_coverage", out_dir, c(class_path), 13.5, 7, family = "Comparison & Ranking", question = "How broadly do ClassA/B/C stable cores extend across pairs?")

  sensitivity$class_threshold <- factor(sensitivity$class_threshold, levels = sort(unique(sensitivity$class_threshold)))
  sensitivity$parameter <- factor(sensitivity$parameter, levels = parameter_levels)
  sensitivity$ratio_class <- factor(sensitivity$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  p5 <- ggplot2::ggplot(sensitivity, ggplot2::aes(x = class_threshold, y = parameter, fill = pair_balanced_mean_proportion)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) + ggplot2::facet_wrap(~ratio_class, nrow = 1) +
    ggplot2::scale_fill_gradient(low = "#FFF9E8", high = "#A97616", limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Symmetric ratio threshold", y = NULL, fill = "Pair-balanced\nseed proportion", title = "ClassA/B/C threshold sensitivity", subtitle = "The primary threshold is 1.1; nearby thresholds show whether rankings depend on the cutoff") + o2jcv_theme(8)
  o2jcv_save_figure(p5, "all_class_threshold_sensitivity", appendix_dir, c(sensitivity_path), 13.5, 7, family = "Uncertainty & Benchmark", question = "Do parameter classifications persist under nearby ratio thresholds?")

  objective$order <- match(objective$objective_stratum, c("all_valid", "delta_le_2", "delta_le_5", "delta_le_10", "top_10pct", "top_25pct", "top_50pct", "best_seed"))
  objective$objective_stratum <- factor(objective$objective_stratum, levels = c("all_valid", "delta_le_2", "delta_le_5", "delta_le_10", "top_10pct", "top_25pct", "top_50pct", "best_seed"))
  objective$parameter <- factor(objective$parameter, levels = parameter_levels)
  objective$ratio_class <- factor(objective$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  p6 <- ggplot2::ggplot(objective, ggplot2::aes(x = objective_stratum, y = parameter, fill = pair_balanced_mean_proportion)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) + ggplot2::facet_wrap(~ratio_class, nrow = 1) +
    ggplot2::scale_fill_gradient(low = "#F5F7F8", high = "#4B708A", limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Objective-quality seed subset", y = NULL, fill = "Pair-balanced\nseed proportion", title = "Objective-stratum robustness", subtitle = "Compare all seeds with near-optimal, top-percentile, and best-seed subsets") + o2jcv_theme(7.8)
  o2jcv_save_figure(p6, "all_class_objective_strata", appendix_dir, c(objective_path), 15, 7, family = "Uncertainty & Benchmark", question = "Are class patterns preserved among better-fitting seeds?")

  boundary_long <- rbind(
    data.frame(boundary[, c("pair_id", "pair_label", "parameter")], diagnostic = "Projection applied", rate = boundary$projection_applied_rate),
    data.frame(boundary[, c("pair_id", "pair_label", "parameter")], diagnostic = "At least one active bound", rate = boundary$any_bound_active_rate),
    data.frame(boundary[, c("pair_id", "pair_label", "parameter")], diagnostic = "Feasible before projection", rate = boundary$feasible_before_projection_rate)
  )
  boundary_long <- boundary_long[is.finite(boundary_long$rate), , drop = FALSE]
  if (nrow(boundary_long)) {
    boundary_long$parameter <- factor(boundary_long$parameter, levels = parameter_levels)
    p7 <- ggplot2::ggplot(boundary_long, ggplot2::aes(x = pair_label, y = parameter, fill = rate)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.2) + ggplot2::facet_wrap(~diagnostic, nrow = 1) +
      ggplot2::scale_fill_gradient(low = "#F5F7F8", high = "#D56A27", limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Seed rate", title = "Boundary and feasibility diagnostics", subtitle = "High rates flag classifications that may be influenced by projection or active bounds") + o2jcv_theme(7.8)
    o2jcv_save_figure(p7, "boundary_projection_sensitivity", appendix_dir, c(boundary_path), 15, 7, family = "Uncertainty & Benchmark", question = "Could optimizer feasibility or active bounds explain apparent class stability?")
  }

  membership <- membership[membership$criterion == "stable90", , drop = FALSE]
  membership$parameter <- factor(membership$parameter, levels = o2jcv_parameter_levels())
  membership$pair_label <- factor(membership$pair_label, levels = rev(sort(unique(membership$pair_label))))
  membership$ratio_class <- factor(membership$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  p8 <- ggplot2::ggplot(membership, ggplot2::aes(x = parameter, y = pair_label)) +
    ggplot2::geom_tile(fill = "#F4F5F5", colour = "white", linewidth = 0.3) +
    ggplot2::geom_point(data = membership[membership$included, , drop = FALSE], shape = 21, fill = "#D4A72C", colour = "#6F5417", size = 3) +
    ggplot2::facet_wrap(~ratio_class, nrow = 1) +
    ggplot2::labs(x = NULL, y = NULL, title = "Stable-90 pair-set membership", subtitle = "UpSet-style matrix: a filled dot means ≥90% of seeds in that pair belong to the displayed class") + o2jcv_theme(8)
  o2jcv_save_figure(p8, "stable90_pair_set_membership", appendix_dir, c(membership_path), 15, 5.5, family = "Matrix & Cohort", question = "Which exact pair sets define each parameter's stable-90 core?")
}

if (identical(environment(), globalenv())) main()
