#!/usr/bin/env Rscript

# Render all-seed fixed-O2 steady-state ploidy curves. Curve class is a facet;
# the separate temporal ploidy Cat is the only color encoding.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "vis", "joint_soft_coupling_stability", "joint_coupling_visualization_utils.R"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[gsub("-", "_", bits[[1L]], fixed = TRUE)]] <- if (length(bits) > 1L) paste(bits[-1L], collapse = "=") else "TRUE"
  }
  out
}

pretty_curve_class <- function(x) {
  labels <- c(
    complex_nonmonotone = "Complex nonmonotone",
    inverted_u_shaped = "Inverted U-shaped",
    monotone_decreasing = "Monotone decreasing",
    monotone_increasing = "Monotone increasing"
  )
  x <- as.character(x)
  fallback <- gsub("_", " ", x, fixed = TRUE)
  fallback <- paste0(toupper(substring(fallback, 1L, 1L)), substring(fallback, 2L))
  ifelse(x %in% names(labels), unname(labels[x]), fallback)
}

args <- parse_args()
analysis_dir <- normalizePath(args$analysis_dir %||% stop("--analysis_dir is required", call. = FALSE), mustWork = TRUE)
out_dir <- normalizePath(args$out_dir %||% stop("--out_dir is required", call. = FALSE), mustWork = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)

plot_path <- file.path(analysis_dir, "fixed_o2_regression_smoothed_curve_plot_data.tsv.gz")
pair_median_path <- file.path(analysis_dir, "fixed_o2_curve_pair_class_median.tsv")
cat_median_path <- file.path(analysis_dir, "fixed_o2_curve_cat_pair_balanced_median.tsv")
pair_summary_path <- file.path(analysis_dir, "fixed_o2_curve_class_summary_by_pair.tsv")
plot_data <- o2jcv_read_tsv(plot_path)
pair_median <- o2jcv_read_tsv(pair_median_path)
cat_median <- o2jcv_read_tsv(cat_median_path)
pair_summary <- o2jcv_read_tsv(pair_summary_path)

required_plot <- c("seed_id", "pair_id", "pair_label", "pair_ploidy_category", "curve_class", "O2_pct", "smoothed_dominant_mean_ploidy")
missing_plot <- setdiff(required_plot, names(plot_data))
if (length(missing_plot)) stop("Fixed-O2 plot table is missing: ", paste(missing_plot, collapse = ", "), call. = FALSE)
cat_order <- c("CatA", "CatB", "CatC", "CatU")
cats <- cat_order[cat_order %in% unique(as.character(plot_data$pair_ploidy_category))]
classes <- unique(as.character(plot_data$curve_class))
classes <- classes[!is.na(classes)]
pair_order <- unique(as.character(pair_summary$pair_id))
pair_lookup <- unique(pair_summary[, c("pair_id", "pair_label", "pair_ploidy_category"), drop = FALSE])
pair_lookup <- pair_lookup[match(pair_order, pair_lookup$pair_id), , drop = FALSE]
pair_lookup$pair_panel <- paste0(pair_lookup$pair_label, " — ", pair_lookup$pair_ploidy_category)

decorate <- function(data) {
  data$pair_ploidy_category <- factor(as.character(data$pair_ploidy_category), levels = cats)
  data$curve_class <- factor(as.character(data$curve_class), levels = classes, labels = pretty_curve_class(classes))
  data$pair_panel <- pair_lookup$pair_panel[match(data$pair_id, pair_lookup$pair_id)]
  data$pair_panel <- factor(data$pair_panel, levels = pair_lookup$pair_panel)
  data
}
plot_data <- decorate(plot_data)
pair_median <- decorate(pair_median)
pair_summary <- decorate(pair_summary)
cat_median$pair_ploidy_category <- factor(as.character(cat_median$pair_ploidy_category), levels = cats)
cat_median$curve_class <- factor(as.character(cat_median$curve_class), levels = classes, labels = pretty_curve_class(classes))

panel_n <- pair_summary[, c("pair_panel", "curve_class", "n_seed"), drop = FALSE]
panel_n$O2_pct <- min(plot_data$O2_pct, na.rm = TRUE)
y_range <- range(plot_data$smoothed_dominant_mean_ploidy, finite = TRUE)
panel_n$label_y <- y_range[[2L]] - 0.02 * diff(y_range)
panel_n$label <- paste0("n=", panel_n$n_seed)

cat_colors <- o2jcv_category_colors()[cats]
pair_figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = O2_pct, y = smoothed_dominant_mean_ploidy,
    group = seed_id, colour = pair_ploidy_category
  )
) +
  ggplot2::geom_line(alpha = 0.10, linewidth = 0.16, na.rm = TRUE) +
  ggplot2::geom_line(
    data = pair_median,
    ggplot2::aes(x = O2_pct, y = median_smoothed_dominant_mean_ploidy, group = interaction(pair_id, curve_class)),
    inherit.aes = FALSE, colour = "#111827", linewidth = 0.7, na.rm = TRUE
  ) +
  ggplot2::geom_text(
    data = panel_n,
    ggplot2::aes(x = O2_pct, y = label_y, label = label),
    inherit.aes = FALSE, hjust = -0.05, vjust = 1, size = 2.4, colour = "#374151"
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(pair_panel), cols = ggplot2::vars(curve_class), drop = FALSE) +
  ggplot2::scale_colour_manual(values = cat_colors, drop = TRUE, na.value = "#8A9299") +
  ggplot2::scale_x_continuous(breaks = 0:5, expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
  ggplot2::labs(
    title = "Fixed-O2 steady-state ploidy curves by pair and curve class",
    subtitle = "All 500 best-fit seeds per pair; curve class is faceted and color identifies the pair's temporal ploidy Cat",
    x = expression("Fixed " * O[2] * " (%)"),
    y = "Regression-smoothed steady-state dominant mean ploidy",
    colour = "Temporal ploidy Cat",
    caption = "Thin lines: individual joint-fit seeds. Black line: pair-by-curve-class seed median."
  ) +
  o2jcv_theme(base_size = 9, rotate_x = FALSE) +
  ggplot2::theme(
    strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
    panel.grid.major.y = ggplot2::element_line(colour = "#EEF1F4", linewidth = 0.2),
    legend.position = "top"
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, linewidth = 1.1)))

o2jcv_save_figure(
  pair_figure,
  "fixed_o2_regression_smoothed_all_seed_curves_by_pair_and_class",
  out_dir,
  c(plot_path, pair_median_path, pair_summary_path),
  width = 15.5, height = 17, dpi = 300,
  family = "Trend / faceted multi-series line",
  question = "Within each pair, which regression-smoothed fixed-O2 steady-state ploidy curve classes occur and how reproducible are their shapes across fitting seeds?"
)

class_n <- aggregate(seed_id ~ curve_class, plot_data, function(x) length(unique(x)))
class_labels <- setNames(
  paste0(as.character(class_n$curve_class), " (n=", class_n$seed_id, ")"),
  as.character(class_n$curve_class)
)
pooled_figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = O2_pct, y = smoothed_dominant_mean_ploidy,
    group = seed_id, colour = pair_ploidy_category
  )
) +
  ggplot2::geom_line(alpha = 0.075, linewidth = 0.16, na.rm = TRUE) +
  ggplot2::geom_line(
    data = cat_median,
    ggplot2::aes(
      x = O2_pct,
      y = pair_balanced_median_smoothed_dominant_mean_ploidy,
      group = pair_ploidy_category, colour = pair_ploidy_category
    ),
    inherit.aes = FALSE, linewidth = 1.05, na.rm = TRUE
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(curve_class), ncol = 2, scales = "fixed",
    labeller = ggplot2::as_labeller(class_labels)
  ) +
  ggplot2::scale_colour_manual(values = cat_colors, drop = TRUE, na.value = "#8A9299") +
  ggplot2::scale_x_continuous(breaks = 0:5, expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
  ggplot2::labs(
    title = "Fixed-O2 steady-state ploidy curves across all pairs",
    subtitle = "All 3,000 best-fit seeds; pair identity is not encoded and color identifies only the temporal ploidy Cat",
    x = expression("Fixed " * O[2] * " (%)"),
    y = "Regression-smoothed steady-state dominant mean ploidy",
    colour = "Temporal ploidy Cat",
    caption = "Thin lines: individual seeds. Thick colored lines: pair-balanced Cat medians within each fixed-O2 curve class."
  ) +
  o2jcv_theme(base_size = 10, rotate_x = FALSE) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(colour = "#EEF1F4", linewidth = 0.2),
    legend.position = "top"
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, linewidth = 1.2)))

o2jcv_save_figure(
  pooled_figure,
  "fixed_o2_regression_smoothed_all_pair_curves_by_class_and_cat",
  out_dir,
  c(plot_path, cat_median_path),
  width = 13.5, height = 7.5, dpi = 300,
  family = "Trend / faceted multi-series line",
  question = "Across all pairs, do temporal ploidy Cats occupy different fixed-O2 steady-state ploidy response shapes within each curve class?"
)

message("Wrote two joint fixed-O2 steady-state ploidy figures to: ", out_dir)
