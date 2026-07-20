#!/usr/bin/env Rscript

# Render pair assignments, category archetypes, representative trajectories,
# diagnostic margins, and threshold robustness from analysis tables.

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
  composition_path <- file.path(analysis_dir, "ploidy_category_summary.tsv")
  trajectory_path <- file.path(analysis_dir, "ploidy_trajectory_category_summary.tsv")
  feature_path <- file.path(analysis_dir, "ploidy_cohort_features.tsv")
  pair_path <- file.path(analysis_dir, "ploidy_pair_category_assignment.tsv")
  archetype_path <- file.path(analysis_dir, "ploidy_category_archetype_summary.tsv")
  representative_path <- file.path(analysis_dir, "ploidy_representative_trajectory_plot_data.tsv")
  confidence_path <- file.path(analysis_dir, "ploidy_classification_confidence.tsv")
  sensitivity_path <- file.path(analysis_dir, "ploidy_category_sensitivity_pair_summary.tsv")
  appendix_dir <- get_arg("appendix_dir") %||% out_dir
  composition <- o2jcv_read_tsv(composition_path)
  trajectory <- o2jcv_add_pair_label(o2jcv_read_tsv(trajectory_path))
  features <- o2jcv_add_pair_label(o2jcv_read_tsv(feature_path))
  pairs <- o2jcv_read_tsv(pair_path)
  archetypes <- o2jcv_read_tsv(archetype_path)
  representative <- o2jcv_add_pair_label(o2jcv_read_tsv(representative_path))
  confidence <- o2jcv_add_pair_label(o2jcv_read_tsv(confidence_path))
  sensitivity <- o2jcv_read_tsv(sensitivity_path)

  pairs$pair_ploidy_category <- factor(pairs$pair_ploidy_category, levels = c("CatA", "CatB", "CatC", "CatU"))
  pairs$pair_label <- factor(pairs$pair_label, levels = rev(pairs$pair_label[order(match(pairs$pair_ploidy_category, c("CatA", "CatB", "CatC", "CatU")), pairs$pair_label)]))
  pairs$label <- sprintf("%s  •  %.0f%% of %d seeds", pairs$pair_ploidy_category, 100 * pairs$dominant_fraction, pairs$n_seed)
  p1 <- ggplot2::ggplot(pairs, ggplot2::aes(x = dominant_fraction, y = pair_label, colour = pair_ploidy_category)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = dominant_fraction, yend = pair_label), colour = "#D8DDE2", linewidth = 1.5) +
    ggplot2::geom_point(size = 4) + ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.05, colour = "#26313A", size = 3.3) +
    ggplot2::scale_colour_manual(values = o2jcv_category_colors(), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 1.35), breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Seeds assigned to the pair's dominant category", y = NULL, title = "Ploidy-trajectory category assigned to each pair", subtitle = "Every pair contains one category only; therefore within-pair CatA/B/C comparisons are not estimable") + o2jcv_theme(10, rotate_x = FALSE)
  o2jcv_save_figure(p1, "ploidy_pair_category_assignment", out_dir, c(pair_path), 10.5, 5, family = "Comparison & Ranking", question = "Which ploidy category characterizes each pair, and is there within-pair category variation?")

  archetypes$ploidy_category <- factor(archetypes$ploidy_category, levels = c("CatA", "CatB", "CatC", "CatU"))
  p_archetype <- ggplot2::ggplot(archetypes, ggplot2::aes(x = day, y = pair_balanced_median, colour = ploidy_category, fill = ploidy_category)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pair_q25, ymax = pair_q75), alpha = 0.16, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::facet_wrap(~cohort, ncol = 1) +
    ggplot2::geom_hline(yintercept = 22, linetype = 3, colour = "#6E7780", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = 44, linetype = 2, colour = "#6E7780", linewidth = 0.35) +
    ggplot2::scale_colour_manual(values = o2jcv_category_colors(), drop = FALSE) +
    ggplot2::scale_fill_manual(values = o2jcv_category_colors(), drop = FALSE) +
    ggplot2::labs(x = "Day", y = "Mean chromosome number", colour = "Ploidy category", fill = "Ploidy category", title = "CatA/B/C ploidy-trajectory archetypes", subtitle = "Pair-balanced category curves; ribbons span the 25th–75th percentile of pair medians") + o2jcv_theme(9, rotate_x = FALSE)
  o2jcv_save_figure(p_archetype, "ploidy_category_archetypes", out_dir, c(archetype_path), 10, 7, family = "Trend", question = "What trajectory shape defines CatA, CatB, and CatC across pairs?")

  trajectory$ploidy_category <- factor(trajectory$ploidy_category, levels = c("CatA", "CatB", "CatC", "CatU"))
  p2 <- ggplot2::ggplot(trajectory, ggplot2::aes(x = day, y = median_ploidy, colour = ploidy_category, fill = ploidy_category)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25_ploidy, ymax = q75_ploidy), alpha = 0.18, colour = NA) +
    ggplot2::geom_line(linewidth = 0.55) + ggplot2::facet_grid(cohort ~ pair_label) +
    ggplot2::geom_hline(yintercept = 22, linetype = 3, colour = "#666666") +
    ggplot2::geom_hline(yintercept = 44, linetype = 2, colour = "#666666") +
    ggplot2::scale_colour_manual(values = o2jcv_category_colors(), drop = FALSE) +
    ggplot2::scale_fill_manual(values = o2jcv_category_colors(), drop = FALSE) +
    ggplot2::labs(x = "Day", y = "Mean chromosome number (median and IQR across seeds)", colour = "Category", fill = "Category", title = "Pair-specific 0–1000 day ploidy trajectories", subtitle = "Each panel summarizes 500 seeds for the indicated initial-ploidy cohort") + o2jcv_theme(8, rotate_x = FALSE)
  o2jcv_save_figure(p2, "ploidy_pair_trajectories", out_dir, c(trajectory_path), 15, 7, family = "Trend", question = "How reproducible is each category trajectory within its assigned pair?")

  representative$ploidy_category <- factor(representative$ploidy_category, levels = c("CatA", "CatB", "CatC", "CatU"))
  annotations <- unique(representative[, c("pair_label", "cohort", "first_drop_start_day", "second_drop_start_day")])
  p_rep <- ggplot2::ggplot(representative, ggplot2::aes(x = day, y = ploidy_value, colour = ploidy_category)) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_vline(data = annotations[is.finite(annotations$first_drop_start_day), , drop = FALSE], ggplot2::aes(xintercept = first_drop_start_day), linetype = 2, colour = "#56616B", linewidth = 0.35) +
    ggplot2::geom_vline(data = annotations[is.finite(annotations$second_drop_start_day), , drop = FALSE], ggplot2::aes(xintercept = second_drop_start_day), linetype = 3, colour = "#56616B", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = 22, linetype = 3, colour = "#9AA2A9", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 44, linetype = 2, colour = "#9AA2A9", linewidth = 0.3) +
    ggplot2::facet_grid(cohort ~ pair_label) + ggplot2::scale_colour_manual(values = o2jcv_category_colors(), guide = "none") +
    ggplot2::labs(x = "Day", y = "Mean chromosome number", title = "Representative trajectory and detected drop starts", subtitle = "One deterministic representative seed per pair; dashed/dotted vertical lines mark first/second detected drop starts") + o2jcv_theme(7.8, rotate_x = FALSE)
  o2jcv_save_figure(p_rep, "ploidy_representative_trajectories", out_dir, c(representative_path), 15, 7, family = "Trend", question = "What does a typical classified trajectory look like, including detected transitions?")

  features$category <- factor(features$category, levels = c("CatA", "CatB", "CatC", "CatU"))
  p3 <- ggplot2::ggplot(features, ggplot2::aes(x = plateau_duration, y = delta_bic_two_minus_one, colour = category, shape = cohort)) +
    ggplot2::geom_hline(yintercept = -10, linetype = 2, colour = "#666666") +
    ggplot2::geom_vline(xintercept = 60, linetype = 2, colour = "#666666") +
    ggplot2::geom_point(alpha = 0.35, size = 1.1) + ggplot2::facet_wrap(~pair_label, ncol = 2) +
    ggplot2::scale_colour_manual(values = o2jcv_category_colors(), drop = FALSE) +
    ggplot2::labs(x = "Longest middle plateau (days)", y = "BIC(two transitions) − BIC(one transition)", colour = "Cohort category", shape = "Initial cohort", title = "CatB/C transition diagnostics", subtitle = "CatC evidence lies right of 60 days and below ΔBIC = −10; points are pair × seed × cohort") + o2jcv_theme(8, rotate_x = FALSE)
  o2jcv_save_figure(p3, "ploidy_category_transition_diagnostics", out_dir, c(feature_path), 11, 8, family = "Relationship", question = "Do plateau duration and change-point evidence separate CatB from CatC?")

  confidence_long <- rbind(
    data.frame(confidence[, c("pair_id", "pair_label", "seed", "cohort", "category")], diagnostic = "High-floor margin (chr)", margin = confidence$margin_high_floor_chr),
    data.frame(confidence[, c("pair_id", "pair_label", "seed", "cohort", "category")], diagnostic = "Low-endpoint margin (chr)", margin = confidence$margin_low_endpoint_chr),
    data.frame(confidence[, c("pair_id", "pair_label", "seed", "cohort", "category")], diagnostic = "Plateau margin (days)", margin = confidence$margin_plateau_days),
    data.frame(confidence[, c("pair_id", "pair_label", "seed", "cohort", "category")], diagnostic = "Two-transition BIC margin", margin = confidence$margin_two_transition_bic)
  )
  confidence_long$category <- factor(confidence_long$category, levels = c("CatA", "CatB", "CatC", "CatU"))
  p4 <- ggplot2::ggplot(confidence_long, ggplot2::aes(x = margin, y = category, fill = cohort)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#56616B", linetype = 2, linewidth = 0.35) +
    ggplot2::geom_boxplot(outlier.alpha = 0.12, width = 0.65, linewidth = 0.35) +
    ggplot2::facet_wrap(~diagnostic, scales = "free_x", ncol = 2) +
    ggplot2::scale_fill_manual(values = c("2N" = "#8FB7C9", "4N" = "#D9A56C")) +
    ggplot2::labs(x = "Signed distance from primary cutoff", y = NULL, fill = "Initial cohort", title = "Classification margins relative to primary cutoffs", subtitle = "Values near zero flag threshold-edge classifications; the biological meaning of the sign differs by diagnostic") + o2jcv_theme(8.5, rotate_x = FALSE)
  o2jcv_save_figure(p4, "ploidy_classification_margins", appendix_dir, c(confidence_path), 11, 7, family = "Uncertainty & Benchmark", question = "How far are classified trajectories from the category thresholds?")

  sensitivity$scenario_category <- factor(sensitivity$scenario_category, levels = c("CatA", "CatB", "CatC", "CatU"))
  sensitivity$pair_label <- factor(sensitivity$pair_label, levels = levels(pairs$pair_label))
  p5 <- ggplot2::ggplot(sensitivity, ggplot2::aes(x = proportion_scenarios, y = pair_label, fill = scenario_category)) +
    ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(data = sensitivity[sensitivity$scenario_category == sensitivity$primary_category, , drop = FALSE], ggplot2::aes(x = 1.02, label = sprintf("primary retained in %.0f%%", 100 * overall_primary_agreement_rate)), hjust = 0, colour = "#26313A", size = 3) +
    ggplot2::scale_fill_manual(values = o2jcv_category_colors(), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1.38), breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Share of 81 prespecified threshold scenarios", y = NULL, fill = "Scenario category", title = "Ploidy-category threshold robustness", subtitle = "Stacked bars show how often each pair maps to CatA/B/C/U across the full sensitivity grid") + o2jcv_theme(9, rotate_x = FALSE)
  o2jcv_save_figure(p5, "ploidy_category_threshold_sensitivity", appendix_dir, c(sensitivity_path), 11, 5.5, family = "Composition", question = "Does each pair retain its primary category under reasonable cutoff changes?")
}

if (identical(environment(), globalenv())) main()
