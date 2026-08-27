#!/usr/bin/env Rscript

.dense_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "render_dense_grid_figures.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.dense_vis_root <- normalizePath(file.path(.dense_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.dense_vis_root, "util", "o2_supply_demand_map_dense_grid_utils.R"), local = environment(), chdir = TRUE)

dense_vis_require <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required for dense-grid figures.", call. = FALSE)
}
dense_vis_palette <- function(classes) stats::setNames(grDevices::hcl.colors(max(length(classes), 1L), "Dark 3"), classes)
dense_vis_save <- function(figure, figure_dir, stub, width = 8, height = 5.5) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf <- file.path(figure_dir, paste0(stub, ".pdf")); png <- file.path(figure_dir, paste0(stub, ".png"))
  ggplot2::ggsave(pdf, figure, width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(png, figure, width = width, height = height, dpi = 220)
  data.frame(figure_id = stub, pdf = pdf, png = png, stringsAsFactors = FALSE)
}
dense_vis_read <- function(out_dir, filename) {
  path <- file.path(dense_grid_analysis_dir(out_dir), filename)
  if (!file.exists(path)) stop("Missing materialized dense-grid analysis table: ", path, call. = FALSE)
  dense_read_tsv(path)
}

dense_render_monotonicity <- function(out_dir) {
  curves <- dense_vis_read(out_dir, "fixed_o2_ploidy_monotonicity_curves.tsv")
  by_seed <- dense_vis_read(out_dir, "fixed_o2_ploidy_monotonicity_by_seed.tsv")
  summary <- dense_vis_read(out_dir, "fixed_o2_ploidy_monotonicity_class_curve_summary.tsv")
  curves <- merge(curves, by_seed[, c("seed_id", "curve_class"), drop = FALSE], by = "seed_id", all.x = TRUE, sort = FALSE)
  classes <- sort(unique(curves$curve_class)); colors <- dense_vis_palette(classes); figure_dir <- dense_grid_figure_dir(out_dir)
  figures <- list()
  p <- ggplot2::ggplot(curves, ggplot2::aes(x = O2_pct, y = dominant_mean_ploidy, group = seed_id, color = curve_class)) + ggplot2::geom_line(alpha = 0.18, linewidth = 0.25) + ggplot2::facet_wrap(~curve_class, scales = "free_y") + ggplot2::scale_color_manual(values = colors) + ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "none") + ggplot2::labs(x = expression(O[2] * " (%)"), y = "Dominant mean ploidy")
  figures[[1L]] <- dense_vis_save(p, figure_dir, "fixed_o2_dominant_ploidy_all_seed_curves_by_class", 10, 7)
  p <- ggplot2::ggplot(summary, ggplot2::aes(x = O2_pct, y = median_dominant_mean_ploidy, color = curve_class, fill = curve_class)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = q25_dominant_mean_ploidy, ymax = q75_dominant_mean_ploidy), alpha = 0.18, color = NA) + ggplot2::geom_line(linewidth = 0.8) + ggplot2::scale_color_manual(values = colors) + ggplot2::scale_fill_manual(values = colors) + ggplot2::theme_bw(base_size = 11) + ggplot2::labs(x = expression(O[2] * " (%)"), y = "Dominant mean ploidy")
  figures[[2L]] <- dense_vis_save(p, figure_dir, "fixed_o2_dominant_ploidy_median_iqr_by_class")
  heat <- curves[order(curves$curve_class, dense_seed_number(curves$seed_id)), , drop = FALSE]
  heat$seed_order <- match(heat$seed_id, unique(heat$seed_id))
  p <- ggplot2::ggplot(heat, ggplot2::aes(x = O2_pct, y = seed_order, fill = dominant_mean_ploidy)) + ggplot2::geom_raster() + ggplot2::scale_fill_viridis_c() + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(x = expression(O[2] * " (%)"), y = "Seeds ordered by curve class", fill = "Ploidy")
  figures[[3L]] <- dense_vis_save(p, figure_dir, "fixed_o2_dominant_ploidy_heatmap_ordered_by_class", 9, 7)
  p <- ggplot2::ggplot(by_seed, ggplot2::aes(x = ploidy_range, y = min_spectral_gap, color = curve_class)) + ggplot2::geom_point(alpha = 0.7) + ggplot2::scale_color_manual(values = colors) + ggplot2::scale_y_continuous(trans = "log10") + ggplot2::theme_bw(base_size = 11) + ggplot2::labs(x = "Ploidy range", y = "Minimum spectral gap")
  figures[[4L]] <- dense_vis_save(p, figure_dir, "fixed_o2_min_spectral_gap_vs_ploidy_range")
  objective_path <- file.path(dense_grid_analysis_dir(out_dir), "fixed_o2_regression_objective_by_curve_class_boxplot_plot_data.tsv")
  if (file.exists(objective_path)) {
    objective <- dense_read_tsv(objective_path)
    if (nrow(objective)) {
      p <- ggplot2::ggplot(objective, ggplot2::aes(x = curve_class, y = objective, fill = curve_class)) + ggplot2::geom_boxplot(outlier.alpha = 0.25) + ggplot2::scale_fill_manual(values = colors) + ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none") + ggplot2::labs(x = NULL, y = "Objective")
      figures[[5L]] <- dense_vis_save(p, figure_dir, "fixed_o2_regression_objective_by_curve_class_boxplot", 9, 5.5)
    }
  }
  dense_rbind_fill(figures)
}

dense_render_initial_ploidy <- function(out_dir) {
  curves <- dense_vis_read(out_dir, "fixed_o2_initial_ploidy_selected_time_curves.tsv")
  delta <- dense_vis_read(out_dir, "fixed_o2_initial_ploidy_delta_by_seed_o2_time.tsv")
  classes <- dense_vis_read(out_dir, "fixed_o2_initial_ploidy_curve_class_by_seed_time.tsv")
  figure_dir <- dense_grid_figure_dir(out_dir); figures <- list()
  summary <- aggregate(curves$mean_ploidy, list(day = curves$day, initial_ploidy = curves$initial_ploidy, O2_pct = curves$O2_pct), stats::median, na.rm = TRUE)
  names(summary)[[4L]] <- "median_mean_ploidy"
  p <- ggplot2::ggplot(summary, ggplot2::aes(x = O2_pct, y = median_mean_ploidy, color = factor(initial_ploidy), group = initial_ploidy)) + ggplot2::geom_line(linewidth = 0.7) + ggplot2::facet_wrap(~day) + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(x = expression(O[2] * " (%)"), y = "Median mean ploidy", color = "Initial ploidy")
  figures[[1L]] <- dense_vis_save(p, figure_dir, "median_mean_ploidy_curves_by_selected_time", 10, 7)
  terminal <- max(delta$day, na.rm = TRUE); terminal_data <- delta[abs(delta$day - terminal) < 1e-9, , drop = FALSE]
  p <- ggplot2::ggplot(terminal_data, ggplot2::aes(x = factor(O2_pct), y = abs_delta_initial)) + ggplot2::geom_boxplot(outlier.alpha = 0.15, fill = "#56B4E9") + ggplot2::theme_bw(base_size = 9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 70, hjust = 1)) + ggplot2::labs(x = expression(O[2] * " (%)"), y = "Absolute 2N-vs-4N ploidy difference")
  figures[[2L]] <- dense_vis_save(p, figure_dir, paste0("day", format(terminal, trim = TRUE), "_convergence_distance_boxplot"), 11, 5.5)
  heat <- terminal_data[order(dense_seed_number(terminal_data$seed_id)), , drop = FALSE]; heat$seed_order <- match(heat$seed_id, unique(heat$seed_id))
  p <- ggplot2::ggplot(heat, ggplot2::aes(x = O2_pct, y = seed_order, fill = abs_delta_initial)) + ggplot2::geom_raster() + ggplot2::scale_fill_viridis_c() + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(x = expression(O[2] * " (%)"), y = "Seed", fill = "Absolute delta")
  figures[[3L]] <- dense_vis_save(p, figure_dir, paste0("day", format(terminal, trim = TRUE), "_delta_initial_heatmap"), 9, 7)
  class_counts <- as.data.frame(table(classes$day, classes$curve_class), stringsAsFactors = FALSE); names(class_counts) <- c("day", "curve_class", "n")
  p <- ggplot2::ggplot(class_counts, ggplot2::aes(x = factor(day), y = n, fill = curve_class)) + ggplot2::geom_col(position = "fill") + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(x = "Day", y = "Fraction", fill = "Curve class")
  figures[[4L]] <- dense_vis_save(p, figure_dir, "curve_class_consistency_by_time", 9, 5.5)
  dense_rbind_fill(figures)
}

dense_visualization_main <- function(argv = dense_parse_args()) {
  dense_vis_require()
  part <- dense_grid_normalize_part(argv$part %||% argv$workflow_part)
  out_dir <- normalizePath(path.expand(argv$out_dir %||% argv$output_root %||% dense_grid_default_out_dir(part)), mustWork = FALSE)
  manifest <- if (part == "monotonicity") dense_render_monotonicity(out_dir) else dense_render_initial_ploidy(out_dir)
  dense_write_tsv(manifest, file.path(dense_grid_figure_dir(out_dir), "figure_manifest.tsv"))
  dense_record_artifact(out_dir, "visualization", paste0("dense_grid_", part, "_figures"), file.path(dense_grid_figure_dir(out_dir), "figure_manifest.tsv"), manifest, "render_dense_grid_figures.R", dense_grid_analysis_dir(out_dir))
  message("Dense-grid visualization complete: ", out_dir)
}

if (sys.nframe() == 0L) dense_visualization_main()
