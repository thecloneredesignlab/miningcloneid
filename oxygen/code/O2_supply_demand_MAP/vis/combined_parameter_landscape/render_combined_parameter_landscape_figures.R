#!/usr/bin/env Rscript

.combined_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "render_combined_parameter_landscape_figures.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.combined_vis_root <- normalizePath(file.path(.combined_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.combined_vis_root, "util", "o2_supply_demand_map_combined_landscape_utils.R"), local = environment(), chdir = TRUE)

combined_class_palette <- function(classes) {
  base <- c(complex_nonmonotone = "#6A5ACD", inverted_u_shaped = "#D55E00", monotone_increasing = "#009E73", monotone_decreasing = "#0072B2", single_transition_increase_then_plateau = "#CC79A7", single_transition_decrease_then_plateau = "#56B4E9", u_shaped = "#E69F00", approximately_flat = "#666666", insufficient_data = "#999999", missing_curve_class = "#BDBDBD")
  missing <- setdiff(classes, names(base))
  if (length(missing)) base <- c(base, stats::setNames(grDevices::hcl.colors(length(missing), "Dark 3"), missing))
  base[classes]
}

combined_plot_theme <- function() ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "right", axis.line = ggplot2::element_line(linewidth = 0.35))

combined_curve_plot <- function(data, xy, reduction) {
  data$.x <- suppressWarnings(as.numeric(data[[xy[[1L]]]])); data$.y <- suppressWarnings(as.numeric(data[[xy[[2L]]]]))
  initial <- data[data$point_type == "initial", , drop = FALSE]
  invitro <- data[data$dataset == "invitro" & data$point_type == "best", , drop = FALSE]
  invivo <- data[data$dataset == "invivo" & data$point_type == "best", , drop = FALSE]
  invivo$curve_class_plot <- ifelse(is.na(invivo$curve_class) | !nzchar(invivo$curve_class), "missing_curve_class", invivo$curve_class)
  classes <- sort(unique(invivo$curve_class_plot))
  ggplot2::ggplot() +
    ggplot2::geom_point(data = initial, ggplot2::aes(.x, .y), color = "#D0D0D0", size = 0.25, alpha = 0.45) +
    ggplot2::geom_point(data = invitro, ggplot2::aes(.x, .y), color = "#D73027", shape = 17, size = 1.15, alpha = 0.9) +
    ggplot2::geom_point(data = invivo, ggplot2::aes(.x, .y, color = curve_class_plot), size = 1.15, alpha = 0.92) +
    ggplot2::scale_color_manual(values = combined_class_palette(classes), labels = gsub("_", " ", classes, fixed = TRUE), name = "O2-ploidy curve class") +
    ggplot2::coord_equal() + ggplot2::labs(x = xy[[1L]], y = xy[[2L]], title = paste(reduction, "curve classes")) + combined_plot_theme()
}

combined_slope_plot <- function(data, xy, reduction) {
  data$.x <- suppressWarnings(as.numeric(data[[xy[[1L]]]])); data$.y <- suppressWarnings(as.numeric(data[[xy[[2L]]]]))
  background <- data[data$point_type == "best", , drop = FALSE]
  invivo <- data[data$dataset == "invivo" & data$point_type == "best" & is.finite(data$average_slope), , drop = FALSE]
  limit <- suppressWarnings(max(abs(invivo$average_slope), na.rm = TRUE)); if (!is.finite(limit) || limit <= 0) limit <- 1
  ggplot2::ggplot() +
    ggplot2::geom_point(data = background, ggplot2::aes(.x, .y), color = "#D0D0D0", size = 0.8, alpha = 0.45) +
    ggplot2::geom_point(data = invivo, ggplot2::aes(.x, .y, color = average_slope), size = 1.25, alpha = 0.95) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-limit, limit), name = "Average slope") +
    ggplot2::coord_equal() + ggplot2::labs(x = xy[[1L]], y = xy[[2L]], title = paste(reduction, "average slope")) + combined_plot_theme()
}

combined_save_pair <- function(plot, stem, width = 7.4, height = 6.2) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  png <- paste0(stem, ".png"); pdf <- paste0(stem, ".pdf")
  ggplot2::ggsave(png, plot, width = width, height = height, units = "in", dpi = 300)
  ggplot2::ggsave(pdf, plot, width = width, height = height, units = "in")
  c(png = normalizePath(png, mustWork = FALSE), pdf = normalizePath(pdf, mustWork = FALSE))
}

render_combined_parameter_landscape_figures <- function(argv = combined_parse_args()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.", call. = FALSE)
  out_dir <- normalizePath(path.expand(argv$out_dir %||% combined_default_out_dir()), mustWork = TRUE)
  manifest_path <- normalizePath(path.expand(argv$analysis_manifest %||% file.path(combined_analysis_dir(out_dir), "combined_embedding_table_manifest.tsv")), mustWork = TRUE)
  manifest <- combined_read_tsv(manifest_path)
  combined_require_columns(manifest, c("reduction", "variant", "stub", "annotated_table"), manifest_path)
  if (combined_as_bool(argv$dry_run, FALSE)) { print(manifest); return(invisible(manifest)) }
  rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    item <- manifest[i, , drop = FALSE]; data <- combined_read_csv(item$annotated_table)
    xy <- combined_coordinate_columns(item$reduction, names(data))
    figure_dir <- file.path(combined_figure_dir(out_dir), item$reduction, item$variant)
    curve <- combined_save_pair(combined_curve_plot(data, xy, item$reduction), file.path(figure_dir, paste0(item$stub, "_curve_class")))
    invivo_curve <- combined_save_pair(combined_curve_plot(data[data$point_type == "best", , drop = FALSE], xy, item$reduction), file.path(figure_dir, paste0(item$stub, "_invivo_best_curve_class")))
    slope <- combined_save_pair(combined_slope_plot(data, xy, item$reduction), file.path(figure_dir, paste0(item$stub, "_invivo_best_average_slope")))
    rows[[i]] <- cbind(item, data.frame(curve_png = curve[["png"]], curve_pdf = curve[["pdf"]], invivo_best_curve_png = invivo_curve[["png"]], invivo_best_curve_pdf = invivo_curve[["pdf"]], average_slope_png = slope[["png"]], average_slope_pdf = slope[["pdf"]], average_slope_highlight_3x3_png = slope[["png"]], average_slope_highlight_3x3_pdf = slope[["pdf"]], stringsAsFactors = FALSE))
  }
  output <- do.call(rbind, rows)
  combined_write_tsv(output, file.path(combined_figure_dir(out_dir), "combined_embedding_figure_manifest.tsv"))
  combined_write_tsv(output, file.path(combined_legacy_table_dir(out_dir), "pooled_embedding_curve_class_manifest.tsv"))
  message("Rendered ", nrow(output), " combined embedding figure sets.")
  invisible(output)
}

main <- render_combined_parameter_landscape_figures
if (sys.nframe() == 0L) main()
