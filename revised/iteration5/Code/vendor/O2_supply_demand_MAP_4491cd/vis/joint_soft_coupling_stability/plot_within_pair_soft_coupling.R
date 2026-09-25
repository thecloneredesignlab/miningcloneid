#!/usr/bin/env Rscript

# Render within-pair class composition, seed distributions, stability cores,
# and direction specifically among ClassB seeds.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(SCRIPT_DIR, "joint_coupling_visualization_utils.R"))

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(name, default = NULL) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default
  }
  analysis_dir <- get_arg("analysis_dir") %||% stop("--analysis_dir is required", call. = FALSE)
  out_dir <- get_arg("out_dir") %||% file.path(analysis_dir, "figures")
  class_path <- file.path(analysis_dir, "within_pair_class_summary.tsv")
  stability_path <- file.path(analysis_dir, "within_pair_parameter_stability.tsv")
  master_path <- file.path(analysis_dir, "soft_coupling_master_long.tsv")
  direction_path <- file.path(analysis_dir, "within_pair_classB_direction_summary.tsv")
  membership_path <- file.path(analysis_dir, "within_pair_class_membership.tsv")
  config_path <- file.path(analysis_dir, "analysis_config.tsv")
  classes <- o2jcv_add_pair_label(o2jcv_read_tsv(class_path))
  stability <- o2jcv_add_pair_label(o2jcv_read_tsv(stability_path))
  master <- o2jcv_add_pair_label(o2jcv_read_tsv(master_path))
  direction_data <- o2jcv_add_pair_label(o2jcv_read_tsv(direction_path))
  membership <- o2jcv_add_pair_label(o2jcv_read_tsv(membership_path))
  config <- o2jcv_read_tsv(config_path)
  class_spec <- o2jcv_classification_spec(config)
  parameter_levels <- rev(o2jcv_parameter_levels())
  classes$parameter <- factor(classes$parameter, levels = parameter_levels)
  classes$ratio_class <- factor(classes$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  p1 <- ggplot2::ggplot(classes, ggplot2::aes(x = proportion, y = parameter, fill = ratio_class)) +
    ggplot2::geom_col(width = 0.78, colour = "white", linewidth = 0.2) + ggplot2::facet_wrap(~pair_label, ncol = 2) +
    ggplot2::scale_fill_manual(values = o2jcv_class_colors(), drop = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    ggplot2::labs(x = "Proportion of available seeds", y = NULL, fill = "Ratio class", title = "Within-pair ClassA/B/C composition", subtitle = paste0("ClassB includes ", o2jcv_classb_label(class_spec), "; each horizontal bar sums to 100%")) + o2jcv_theme(9, rotate_x = FALSE)
  o2jcv_save_figure(p1, "within_pair_class_composition", out_dir, c(class_path), 11, 10.5, family = "Composition", question = "How are 500 seeds distributed across ClassA/B/C within each pair?")

  stability$parameter <- factor(stability$parameter, levels = parameter_levels)
  p2 <- ggplot2::ggplot(stability, ggplot2::aes(x = pair_label, y = parameter, fill = normalized_entropy)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", normalized_entropy)), size = 2.5, colour = "#26313A") +
    ggplot2::scale_fill_gradient(low = "#FFF9E8", high = "#B27A18", limits = c(0, 1), na.value = "grey90") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Normalized\nentropy", title = "Within-pair class instability", subtitle = "0 means one class across all seeds; 1 means maximally mixed A/B/C") + o2jcv_theme(9)
  o2jcv_save_figure(p2, "within_pair_class_entropy", out_dir, c(stability_path), 9.5, 6.5, family = "Matrix & Cohort", question = "Which parameter classifications are sensitive to seed initialization?")

  direction <- rbind(
    data.frame(direction_data[, c("pair_id", "pair_label", "parameter", "n_ClassB")], direction = "ClassB ratio < 1", proportion = direction_data$ClassB_prop_ratio_lt1),
    data.frame(direction_data[, c("pair_id", "pair_label", "parameter", "n_ClassB")], direction = "ClassB ratio = 1", proportion = direction_data$ClassB_prop_ratio_eq1),
    data.frame(direction_data[, c("pair_id", "pair_label", "parameter", "n_ClassB")], direction = "ClassB ratio > 1", proportion = direction_data$ClassB_prop_ratio_gt1)
  )
  direction <- direction[direction$n_ClassB > 0, , drop = FALSE]
  direction$parameter <- factor(direction$parameter, levels = parameter_levels)
  p3 <- ggplot2::ggplot(direction, ggplot2::aes(x = proportion, y = parameter, fill = direction)) +
    ggplot2::geom_col(width = 0.78, colour = "white", linewidth = 0.2) + ggplot2::facet_wrap(~pair_label, ncol = 2) +
    ggplot2::scale_fill_manual(values = c("ClassB ratio < 1" = "#2F6B9A", "ClassB ratio = 1" = "#F2EEE1", "ClassB ratio > 1" = "#D56A27")) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(x = "Direction among ClassB seeds", y = NULL, fill = NULL, title = "Direction inside the ClassB coupling band", subtitle = "Only seeds already classified as ClassB are included; parameters with zero ClassB seeds are omitted") + o2jcv_theme(9, rotate_x = FALSE)
  o2jcv_save_figure(p3, "within_pair_classB_direction", out_dir, c(direction_path), 11, 10.5, family = "Composition", question = "When a parameter is coupled, is the in-vivo value consistently above or below its in-vitro value?")

  master$parameter <- factor(master$parameter, levels = parameter_levels)
  p4 <- ggplot2::ggplot(master, ggplot2::aes(x = log2_ratio_vivo_to_vitro, y = parameter)) +
    o2jcv_class_band(class_spec$lower, class_spec$upper) +
    ggplot2::geom_violin(fill = "#9CB8CC", colour = "#46677F", linewidth = 0.25, scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "#26313A", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#343A40", linewidth = 0.35) +
    ggplot2::facet_wrap(~pair_label, ncol = 2) +
    ggplot2::labs(x = "log2(in vivo / in vitro)", y = NULL, title = "Seed-level ratio distributions within each pair", subtitle = "Gold band is ClassB; violin shape shows all available seeds and the box marks median and IQR") + o2jcv_theme(8.5, rotate_x = FALSE)
  o2jcv_save_figure(p4, "within_pair_ratio_distributions", out_dir, c(master_path, config_path), 11.5, 11, family = "Distribution", question = "Are class assignments driven by compact or boundary-spanning seed distributions?")

  membership$parameter <- factor(membership$parameter, levels = parameter_levels)
  membership$ratio_class <- factor(membership$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  membership$stability_level <- factor(membership$stability_level, levels = c("Absent", "Union only", "Stable >=80%", "Stable >=90%", "Stable >=95%", "Strict 100%"))
  stability_colors <- c("Absent" = "#F4F4F2", "Union only" = "#F5E9C5", "Stable >=80%" = "#E8CA78", "Stable >=90%" = "#D4A72C", "Stable >=95%" = "#B78720", "Strict 100%" = "#7B5B18")
  membership$label_colour <- ifelse(membership$stability_level %in% c("Stable >=95%", "Strict 100%"), "white", "#26313A")
  p5 <- ggplot2::ggplot(membership, ggplot2::aes(x = pair_label, y = parameter, fill = stability_level)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * proportion), colour = label_colour), size = 2.1) +
    ggplot2::scale_colour_identity() +
    ggplot2::facet_grid(. ~ ratio_class) +
    ggplot2::scale_fill_manual(values = stability_colors, drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Seed stability", title = "Union, graded core, and strict intersection by class", subtitle = "Labels give seed proportion; darkest cells are classified identically in all 500 seeds") + o2jcv_theme(8)
  o2jcv_save_figure(p5, "within_pair_class_stability_levels", out_dir, c(membership_path), 14, 7.2, family = "Matrix & Cohort", question = "Which parameters belong to each class as a union, graded core, or strict intersection?")
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
if (identical(environment(), globalenv())) main()
