#!/usr/bin/env Rscript

# Visualization-only overview of pair-level soft-coupling ratios and ploidy
# categories. All statistics are read from materialized analysis tables.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(SCRIPT_DIR, "joint_coupling_visualization_utils.R"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(name, default = NULL) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default
  }
  ratio_dir <- get_arg("ratio_analysis_dir") %||% stop("--ratio_analysis_dir is required", call. = FALSE)
  ploidy_dir <- get_arg("ploidy_analysis_dir") %||% stop("--ploidy_analysis_dir is required", call. = FALSE)
  out_dir <- get_arg("out_dir") %||% stop("--out_dir is required", call. = FALSE)
  stability_path <- file.path(ratio_dir, "within_pair_parameter_stability.tsv")
  config_path <- file.path(ratio_dir, "analysis_config.tsv")
  pair_path <- file.path(ploidy_dir, "ploidy_pair_category_assignment.tsv")
  stability <- o2jcv_add_pair_label(o2jcv_read_tsv(stability_path))
  pairs <- o2jcv_read_tsv(pair_path)
  class_spec <- o2jcv_classification_spec(o2jcv_read_tsv(config_path))
  stability$ploidy_category <- pairs$pair_ploidy_category[match(stability$pair_id, pairs$pair_id)]
  stability$pair_display <- paste0(stability$pair_label, "  [", stability$ploidy_category, "]")
  pair_order <- unique(stability[order(match(stability$ploidy_category, c("CatA", "CatB", "CatC")), stability$pair_label), "pair_display"])
  stability$pair_display <- factor(stability$pair_display, levels = pair_order)
  stability$parameter <- factor(stability$parameter, levels = rev(o2jcv_parameter_levels()))

  bound <- max(abs(stability$log2_ratio_median), na.rm = TRUE)
  stability$ratio_label_colour <- ifelse(abs(stability$log2_ratio_median) >= 0.55 * bound, "white", "#26313A")
  p1 <- ggplot2::ggplot(stability, ggplot2::aes(x = pair_display, y = parameter, fill = log2_ratio_median)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", log2_ratio_median), colour = ratio_label_colour), size = 2.6) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient2(low = "#2F6B9A", mid = "#FFFDF7", high = "#D56A27", midpoint = 0, limits = c(-bound, bound)) +
    ggplot2::labs(
      x = NULL, y = NULL, fill = "Median log2\n(vivo/vitro)",
      title = "Pair-level in-vivo / in-vitro ratio map",
      subtitle = "Each cell is the median across 500 fitting seeds; pair labels include the pair-level ploidy category",
      caption = "Negative values indicate lower in vivo values; positive values indicate higher in vivo values."
    ) + o2jcv_theme(9)
  o2jcv_save_figure(
    p1, "overview_pair_parameter_ratio_map", out_dir, c(stability_path, pair_path), 10.5, 7.2,
    family = "Matrix & Cohort", question = "Which parameters differ in direction and magnitude across pairs?"
  )

  stability$class_letter <- sub("Class", "", stability$dominant_class)
  stability$consensus_label <- sprintf("%s\n%.0f%%", stability$class_letter, 100 * stability$dominant_fraction)
  stability$class_label_colour <- ifelse(stability$dominant_class == "ClassB", "#26313A", "white")
  p2 <- ggplot2::ggplot(stability, ggplot2::aes(x = pair_display, y = parameter, fill = dominant_class)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = consensus_label, colour = class_label_colour), size = 2.7, lineheight = 0.9) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_manual(values = o2jcv_class_colors(), drop = FALSE) +
    ggplot2::labs(
      x = NULL, y = NULL, fill = "Dominant class",
      title = "Dominant coupling class and seed agreement",
      subtitle = paste0("Cell labels show Class A/B/C and seed agreement; ClassB uses ", o2jcv_classb_label(class_spec))
    ) + o2jcv_theme(9)
  o2jcv_save_figure(
    p2, "overview_dominant_class_matrix", out_dir, c(stability_path, pair_path, config_path), 10.5, 7.2,
    family = "Matrix & Cohort", question = "Which classes are stable within each pair and how do they align with ploidy category?"
  )
}

if (identical(environment(), globalenv())) main()
