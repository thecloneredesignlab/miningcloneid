#!/usr/bin/env Rscript

# Canonical visualization consumer for the combined FixO2 eigen-attractor workflow.

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

WORKFLOW_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
UTIL_DIR <- file.path(WORKFLOW_DIR, "util")
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_path_utils.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_cli_utils.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_combined_fixo2_eigen_utils.R"))

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript plot_fixo2_eigen_attractor_embedding_curve_class.R [options]",
      "",
      "Purpose:",
      "  Render FixO2 eigen-attractor plots from analysis-prepared annotated coordinate tables.",
      "",
      "Inputs:",
      "  --analysis_manifest=FILE  Analysis manifest containing annotated_coordinate_table paths.",
      "",
      "Outputs:",
      "  --out_dir=DIR           Output root for curve-class figures and tables.",
      "  --reductions=LIST       Optional reduction subset. Default: all manifest rows.",
      "  --variants=LIST         Optional variant subset. Default: all manifest rows.",
      "  --dry_run=TRUE|FALSE    Print planned inputs/outputs without writing files.",
      sep = "\n"
    ),
    "\n"
  )
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}

read_csv_plain <- o2cfe_read_csv_plain
read_tsv_plain <- o2cfe_read_tsv_plain
write_tsv <- o2cfe_write_tsv
normalize_reductions <- o2cfe_normalize_reductions
normalize_variants <- o2cfe_normalize_variants

curve_class_colors <- function(classes) {
  pal <- c(
    complex_nonmonotone = "#6A5ACD",
    inverted_u_shaped = "#D55E00",
    monotone_increasing = "#009E73",
    monotone_decreasing = "#0072B2",
    single_transition_increase_then_plateau = "#CC79A7",
    single_transition_decrease_then_plateau = "#56B4E9",
    u_shaped = "#E69F00",
    approximately_flat = "#666666",
    insufficient_data = "#999999",
    missing_curve_class = "#999999"
  )
  missing <- setdiff(classes, names(pal))
  if (length(missing)) {
    pal <- c(pal, stats::setNames(grDevices::rainbow(length(missing)), missing))
  }
  pal[classes]
}

curve_class_order <- function(classes) {
  base <- c(
    "complex_nonmonotone",
    "inverted_u_shaped",
    "monotone_increasing",
    "u_shaped",
    "approximately_flat",
    "single_transition_increase_then_plateau",
    "single_transition_decrease_then_plateau",
    "monotone_decreasing",
    "insufficient_data",
    "missing_curve_class"
  )
  classes <- unique(as.character(classes))
  classes <- classes[!is.na(classes) & nzchar(classes)]
  c(intersect(base, classes), sort(setdiff(classes, base)))
}

coord_spec <- function(reduction, data_names) {
  candidates <- switch(
    reduction,
    PCAs = list(c("PCA1", "PCA2")),
    UMAPs = list(c("UMAP1", "UMAP2")),
    TSNEs = list(c("tSNE1", "tSNE2"), c("TSNE1", "TSNE2")),
    stop("Unknown reduction: ", reduction, call. = FALSE)
  )
  for (candidate in candidates) {
    if (all(candidate %in% data_names)) {
      labels <- switch(
        reduction,
        PCAs = c("PCA 1", "PCA 2"),
        UMAPs = c("UMAP 1", "UMAP 2"),
        TSNEs = c("t-SNE 1", "t-SNE 2")
      )
      return(list(columns = candidate, labels = labels))
    }
  }
  stop(
    "Could not find coordinate columns for ",
    reduction,
    ". Available columns: ",
    paste(data_names, collapse = ", "),
    call. = FALSE
  )
}

square_limits <- function(plot_data, coord_names) {
  x <- suppressWarnings(as.numeric(plot_data[[coord_names[[1L]]]]))
  y <- suppressWarnings(as.numeric(plot_data[[coord_names[[2L]]]]))
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) stop("No finite embedding coordinates.", call. = FALSE)
  x_range <- range(x)
  y_range <- range(y)
  x_center <- mean(x_range)
  y_center <- mean(y_range)
  span <- max(diff(x_range), diff(y_range))
  if (!is.finite(span) || span <= 0) span <- 1
  pad <- span * 0.03
  span <- span + 2 * pad
  list(
    xlim = c(x_center - span / 2, x_center + span / 2),
    ylim = c(y_center - span / 2, y_center + span / 2)
  )
}

class_labels <- function(classes) {
  labels <- gsub("_", " ", classes, fixed = TRUE)
  labels[classes == "missing_curve_class"] <- "missing curve class"
  stats::setNames(labels, classes)
}

prepare_plot_data <- function(plot_data, coord_names) {
  plot_data$.embedding_x <- suppressWarnings(as.numeric(plot_data[[coord_names[[1L]]]]))
  plot_data$.embedding_y <- suppressWarnings(as.numeric(plot_data[[coord_names[[2L]]]]))
  plot_data$dataset <- as.character(plot_data$dataset)
  plot_data$point_type <- as.character(plot_data$point_type)
  plot_data
}

initial_points <- function(plot_data) {
  initial <- plot_data[plot_data$point_type == "initial", , drop = FALSE]
  initial$initial_group <- ifelse(initial$dataset == "invivo", "invivo_initial", "invitro_initial")
  initial$initial_group <- factor(initial$initial_group, levels = c("invivo_initial", "invitro_initial"))
  initial
}

initial_legend_points <- function(lims) {
  data.frame(
    .embedding_x = rep(mean(lims$xlim), 2L),
    .embedding_y = rep(mean(lims$ylim), 2L),
    initial_group = factor(c("invivo_initial", "invitro_initial"), levels = c("invivo_initial", "invitro_initial"))
  )
}

common_embedding_theme <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      legend.box = "vertical",
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

build_objective_plot <- function(plot_data,
                                 coord_names,
                                 axis_labels,
                                 lims,
                                 initial_size = 0.22,
                                 best_size = 1.25,
                                 initial_invivo_color = "#7FA2FF",
                                 initial_invitro_color = "#FF4FBF",
                                 initial_invivo_alpha = 0.48,
                                 initial_invitro_alpha = 0.45) {
  plot_data <- prepare_plot_data(plot_data, coord_names)
  initial <- initial_points(plot_data)
  initial_invivo <- initial[initial$dataset == "invivo", , drop = FALSE]
  initial_invitro <- initial[initial$dataset == "invitro", , drop = FALSE]
  initial_legend <- initial_legend_points(lims)
  best_invivo <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invivo", , drop = FALSE]
  best_invitro <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invitro", , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = initial_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      color = initial_invivo_color,
      shape = 16,
      size = initial_size,
      alpha = initial_invivo_alpha,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      color = initial_invitro_color,
      shape = 17,
      size = initial_size,
      alpha = initial_invitro_alpha,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_legend,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = initial_group, shape = initial_group),
      color = "transparent",
      alpha = 0,
      size = 0.01,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Initial",
      values = c(invivo_initial = initial_invivo_color, invitro_initial = initial_invitro_color),
      labels = c(invivo_initial = "in vivo", invitro_initial = "in vitro"),
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          alpha = 1,
          size = 3,
          shape = c(21, 24),
          color = "transparent"
        )
      )
    ) +
    ggplot2::scale_shape_manual(
      name = "Initial",
      values = c(invivo_initial = 21, invitro_initial = 24),
      labels = c(invivo_initial = "in vivo", invitro_initial = "in vitro"),
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::geom_point(
      data = best_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = objective),
      shape = 16,
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_gradient(
      name = "in vivo\nobjective",
      low = "#2C7BB6",
      high = "#FDE725",
      guide = ggplot2::guide_colorbar(order = 2, barheight = ggplot2::unit(26, "mm"), display = "rectangles")
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      data = best_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = objective),
      shape = 17,
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_gradient(
      name = "in vitro\nobjective",
      low = "#1A9850",
      high = "#D73027",
      guide = ggplot2::guide_colorbar(order = 3, barheight = ggplot2::unit(26, "mm"), display = "rectangles")
    ) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(x = axis_labels[[1L]], y = axis_labels[[2L]]) +
    common_embedding_theme()
}

build_curve_class_plot <- function(plot_data,
                                   reduction,
                                   coord_names,
                                   axis_labels,
                                   lims,
                                   initial_size = 0.22,
                                   best_size = 1.25,
                                   initial_invivo_color = "#7FA2FF",
                                   initial_invitro_color = "#FF4FBF",
                                   initial_invivo_alpha = 0.48,
                                   initial_invitro_alpha = 0.45,
                                   invitro_best_color = "#D73027") {
  plot_data <- prepare_plot_data(plot_data, coord_names)
  initial <- initial_points(plot_data)
  initial_invivo <- initial[initial$dataset == "invivo", , drop = FALSE]
  initial_invitro <- initial[initial$dataset == "invitro", , drop = FALSE]
  initial_legend <- initial_legend_points(lims)
  best_invivo <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invivo", , drop = FALSE]
  best_invitro <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invitro", , drop = FALSE]

  best_invivo$curve_class_for_plot <- ifelse(
    is.na(best_invivo$curve_class) | !nzchar(best_invivo$curve_class),
    "missing_curve_class",
    best_invivo$curve_class
  )
  class_order <- sort(setdiff(unique(best_invivo$curve_class_for_plot), "missing_curve_class"))
  if ("missing_curve_class" %in% best_invivo$curve_class_for_plot) {
    class_order <- c(class_order, "missing_curve_class")
  }
  best_invivo$curve_class_for_plot <- factor(best_invivo$curve_class_for_plot, levels = class_order)
  best_fit_legend <- data.frame(
    .embedding_x = rep(mean(lims$xlim), 2L),
    .embedding_y = rep(mean(lims$ylim), 2L),
    best_fit_group = factor(c("invivo_best", "invitro_best"), levels = c("invivo_best", "invitro_best"))
  )
  class_values <- curve_class_colors(class_order)
  labels <- class_labels(class_order)

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = initial_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      color = initial_invivo_color,
      shape = 16,
      size = initial_size,
      alpha = initial_invivo_alpha,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      color = initial_invitro_color,
      shape = 17,
      size = initial_size,
      alpha = initial_invitro_alpha,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_legend,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = initial_group, shape = initial_group),
      color = "transparent",
      alpha = 0,
      size = 0.01,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Initial",
      values = c(invivo_initial = initial_invivo_color, invitro_initial = initial_invitro_color),
      labels = c(invivo_initial = "in vivo", invitro_initial = "in vitro"),
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          alpha = 1,
          size = 3,
          shape = c(21, 24),
          color = "transparent"
        )
      )
    ) +
    ggplot2::scale_shape_manual(
      name = "Initial",
      values = c(invivo_initial = 21, invitro_initial = 24),
      labels = c(invivo_initial = "in vivo", invitro_initial = "in vitro"),
      drop = FALSE,
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
      data = best_fit_legend,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = best_fit_group),
      shape = 21,
      color = "transparent",
      alpha = 0,
      size = 0.01,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Best fits",
      values = c(invivo_best = "#555555", invitro_best = invitro_best_color),
      labels = c(invivo_best = "in vivo", invitro_best = "in vitro"),
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          shape = c(16, 17),
          color = c("#555555", invitro_best_color),
          fill = c("#555555", invitro_best_color),
          alpha = 1,
          size = 3
        )
      )
    ) +
    ggplot2::geom_point(
      data = best_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 17,
      color = invitro_best_color,
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = best_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = curve_class_for_plot),
      shape = 16,
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = expression(O[2] * "-Ploidy Curve Class"),
      values = class_values,
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 3,
        override.aes = list(shape = 16, alpha = 1, size = 3)
      )
    ) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(
      x = axis_labels[[1L]],
      y = axis_labels[[2L]]
    ) +
    common_embedding_theme()
}

prepare_invivo_best_class_data <- function(plot_data, class_order = NULL) {
  best_invivo <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invivo", , drop = FALSE]
  best_invivo$curve_class_for_plot <- ifelse(
    is.na(best_invivo$curve_class) | !nzchar(best_invivo$curve_class),
    "missing_curve_class",
    best_invivo$curve_class
  )
  if (is.null(class_order)) class_order <- curve_class_order(best_invivo$curve_class_for_plot)
  best_invivo$curve_class_for_plot <- factor(best_invivo$curve_class_for_plot, levels = class_order)
  best_invivo
}

sparse_axis_breaks <- function(lim, n = 3L) {
  breaks <- pretty(lim, n = n)
  breaks[breaks >= min(lim) & breaks <= max(lim)]
}

sparse_axis_scales <- function(lims, n = 3L) {
  list(
    ggplot2::scale_x_continuous(breaks = sparse_axis_breaks(lims$xlim, n = n)),
    ggplot2::scale_y_continuous(breaks = sparse_axis_breaks(lims$ylim, n = n))
  )
}

finite_range <- function(x, fallback = c(-1, 1)) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(fallback)
  lim <- range(x)
  if (!all(is.finite(lim))) return(fallback)
  if (diff(lim) <= 0) {
    pad <- max(abs(lim[[1L]]) * 0.05, 0.01)
    lim <- lim + c(-pad, pad)
  }
  lim
}

compact_embedding_theme <- function() {
  common_embedding_theme() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 9.8, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 5.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.2),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

build_invivo_best_curve_class_plot <- function(plot_data,
                                               coord_names,
                                               axis_labels,
                                               lims,
                                               best_size = 1.25,
                                               class_order = NULL) {
  plot_data <- prepare_plot_data(plot_data, coord_names)
  best_invivo <- prepare_invivo_best_class_data(plot_data, class_order = class_order)
  class_order <- levels(best_invivo$curve_class_for_plot)
  class_values <- curve_class_colors(class_order)
  labels <- class_labels(class_order)

  ggplot2::ggplot(best_invivo, ggplot2::aes(x = .embedding_x, y = .embedding_y, color = curve_class_for_plot)) +
    ggplot2::geom_point(shape = 16, alpha = 0.95, size = best_size, stroke = 0, show.legend = TRUE) +
    ggplot2::scale_color_manual(
      name = expression(O[2] * "-Ploidy Curve Class"),
      values = class_values,
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(override.aes = list(shape = 16, alpha = 1, size = 3))
    ) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(x = axis_labels[[1L]], y = axis_labels[[2L]]) +
    common_embedding_theme()
}

build_curve_class_highlight_panel <- function(plot_data,
                                              coord_names,
                                              lims,
                                              class_order,
                                              highlight_class,
                                              highlight_size = 1.25,
                                              background_size = 1.05,
                                              grey_color = "#C9C9C9") {
  plot_data <- prepare_plot_data(plot_data, coord_names)
  best_invivo <- prepare_invivo_best_class_data(plot_data, class_order = class_order)
  other <- best_invivo[as.character(best_invivo$curve_class_for_plot) != highlight_class, , drop = FALSE]
  current <- best_invivo[as.character(best_invivo$curve_class_for_plot) == highlight_class, , drop = FALSE]
  class_values <- curve_class_colors(class_order)

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = other,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 16,
      color = grey_color,
      alpha = 0.45,
      size = background_size,
      stroke = 0
    ) +
    ggplot2::geom_point(
      data = current,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 16,
      color = unname(class_values[[highlight_class]]),
      alpha = 0.95,
      size = highlight_size,
      stroke = 0
    ) +
    sparse_axis_scales(lims, n = 2L) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(title = unname(class_labels(highlight_class))) +
    compact_embedding_theme()
}

build_curve_class_legend_panel <- function(class_order) {
  d <- data.frame(
    .embedding_x = 0,
    .embedding_y = seq_along(class_order),
    curve_class_for_plot = factor(class_order, levels = class_order)
  )
  ggplot2::ggplot(d, ggplot2::aes(x = .embedding_x, y = .embedding_y, color = curve_class_for_plot)) +
    ggplot2::geom_point(alpha = 0, size = 0.01, show.legend = TRUE) +
    ggplot2::scale_color_manual(
      name = expression(O[2] * "-Ploidy Curve Class"),
      values = curve_class_colors(class_order),
      labels = class_labels(class_order),
      drop = FALSE,
      guide = ggplot2::guide_legend(
        override.aes = list(alpha = 1, size = 2.5, shape = 16),
        keyheight = ggplot2::unit(3.4, "mm"),
        keywidth = ggplot2::unit(3.4, "mm")
      )
    ) +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      legend.position = c(0.5, 0.5),
      legend.justification = "center",
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 6.2),
      legend.key.size = ggplot2::unit(3.4, "mm"),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

build_curve_class_highlight_3x3 <- function(plot_data,
                                            coord_names,
                                            lims,
                                            class_order,
                                            highlight_size = 1.25,
                                            background_size = 1.05) {
  panels <- lapply(class_order, function(cls) {
    build_curve_class_highlight_panel(
      plot_data = plot_data,
      coord_names = coord_names,
      lims = lims,
      class_order = class_order,
      highlight_class = cls,
      highlight_size = highlight_size,
      background_size = background_size
    )
  })
  panels <- c(panels, list(build_curve_class_legend_panel(class_order)))
  patchwork::wrap_plots(panels, ncol = 3)
}

build_average_slope_highlight_panel <- function(plot_data,
                                                coord_names,
                                                lims,
                                                class_order,
                                                highlight_class,
                                                slope_limits,
                                                highlight_size = 1.25,
                                                background_size = 1.05,
                                                show_legend = FALSE,
                                                grey_color = "#C9C9C9") {
  plot_data <- prepare_plot_data(plot_data, coord_names)
  best_invivo <- prepare_invivo_best_class_data(plot_data, class_order = class_order)
  best_invivo$average_slope <- suppressWarnings(as.numeric(best_invivo$average_slope))
  other <- best_invivo[as.character(best_invivo$curve_class_for_plot) != highlight_class, , drop = FALSE]
  current <- best_invivo[as.character(best_invivo$curve_class_for_plot) == highlight_class, , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = other,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 16,
      color = grey_color,
      alpha = 0.45,
      size = background_size,
      stroke = 0
    ) +
    ggplot2::geom_point(
      data = current,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = average_slope),
      shape = 21,
      color = "transparent",
      alpha = 0.95,
      size = highlight_size,
      stroke = 0
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = slope_limits,
      oob = scales::squish,
      name = "Average\nslope",
      guide = if (isTRUE(show_legend)) {
        ggplot2::guide_colorbar(barheight = ggplot2::unit(30, "mm"), display = "rectangles")
      } else {
        "none"
      }
    ) +
    sparse_axis_scales(lims, n = 2L) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(title = unname(class_labels(highlight_class))) +
    compact_embedding_theme() +
    ggplot2::theme(
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7)
    )
}

build_average_slope_legend_panel <- function(slope_limits) {
  d <- data.frame(
    .embedding_x = 0,
    .embedding_y = 0,
    average_slope = mean(slope_limits)
  )
  ggplot2::ggplot(d, ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = average_slope)) +
    ggplot2::geom_point(shape = 21, color = "transparent", alpha = 0, size = 0.01, show.legend = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = slope_limits,
      oob = scales::squish,
      name = "Average\nslope",
      guide = ggplot2::guide_colorbar(
        barheight = ggplot2::unit(30, "mm"),
        barwidth = ggplot2::unit(3.8, "mm"),
        display = "rectangles"
      )
    ) +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      legend.position = c(0.5, 0.5),
      legend.justification = "center",
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 6.2),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

build_average_slope_highlight_3x3 <- function(plot_data,
                                              coord_names,
                                              lims,
                                              class_order,
                                              slope_limits,
                                              highlight_size = 1.25,
                                              background_size = 1.05) {
  panels <- lapply(class_order, function(cls) {
    build_average_slope_highlight_panel(
      plot_data = plot_data,
      coord_names = coord_names,
      lims = lims,
      class_order = class_order,
      highlight_class = cls,
      slope_limits = slope_limits,
      highlight_size = highlight_size,
      background_size = background_size
    )
  })
  panels <- c(panels, list(build_average_slope_legend_panel(slope_limits)))
  patchwork::wrap_plots(panels, ncol = 3)
}

build_monotone_increasing_slope_panel <- function(plot_data,
                                                  coord_names,
                                                  lims,
                                                  class_order,
                                                  slope_limits,
                                                  highlight_size = 2.35,
                                                  background_size = 1.15,
                                                  grey_color = "#C9C9C9") {
  plot_data <- prepare_plot_data(plot_data, coord_names)
  best_invivo <- prepare_invivo_best_class_data(plot_data, class_order = class_order)
  best_invivo$average_slope <- suppressWarnings(as.numeric(best_invivo$average_slope))
  current <- best_invivo[as.character(best_invivo$curve_class_for_plot) == "monotone_increasing", , drop = FALSE]
  other <- best_invivo[as.character(best_invivo$curve_class_for_plot) != "monotone_increasing", , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = other,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 16,
      color = grey_color,
      alpha = 0.45,
      size = background_size,
      stroke = 0
    ) +
    ggplot2::geom_point(
      data = current,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = average_slope),
      shape = 21,
      color = "transparent",
      alpha = 0.95,
      size = highlight_size,
      stroke = 0
    ) +
    ggplot2::scale_fill_gradient(
      low = "#F7F7F7",
      high = "#B2182B",
      limits = slope_limits,
      oob = scales::squish,
      name = "Average\nslope",
      guide = ggplot2::guide_colorbar(barheight = ggplot2::unit(32, "mm"), display = "rectangles")
    ) +
    sparse_axis_scales(lims, n = 2L) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(title = "monotone increasing average slope") +
    compact_embedding_theme() +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7)
    )
}

build_monotone_increasing_combined <- function(plot_data,
                                               coord_names,
                                               lims,
                                               class_order,
                                               slope_limits,
                                               highlight_size = 2.35,
                                               background_size = 1.15) {
  class_panel <- build_curve_class_highlight_panel(
    plot_data = plot_data,
    coord_names = coord_names,
    lims = lims,
    class_order = class_order,
    highlight_class = "monotone_increasing",
    highlight_size = highlight_size,
    background_size = background_size
  ) +
    ggplot2::labs(title = "monotone increasing curve class")
  slope_panel <- build_monotone_increasing_slope_panel(
    plot_data = plot_data,
    coord_names = coord_names,
    lims = lims,
    class_order = class_order,
    slope_limits = slope_limits,
    highlight_size = highlight_size,
    background_size = background_size
  )
  class_panel + slope_panel + patchwork::plot_layout(ncol = 2, widths = c(1, 1.08))
}

save_plot_pair <- function(plot, out_stem, width = 7.4, height = 6.5, dpi = 300) {
  png_path <- paste0(out_stem, ".png")
  pdf_path <- paste0(out_stem, ".pdf")
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  c(png = png_path, pdf = pdf_path)
}

process_embedding <- function(row, class_map, slope_map, average_slope_table, out_dir) {
  if (!"annotated_coordinate_table" %in% names(row)) {
    stop("Analysis manifest is missing annotated_coordinate_table.", call. = FALSE)
  }
  plot_data <- read_csv_plain(row$annotated_coordinate_table)
  required_annotations <- c("seed_number", "curve_class", "average_slope")
  missing_annotations <- setdiff(required_annotations, names(plot_data))
  if (length(missing_annotations)) {
    stop(
      "Annotated coordinate table is missing analysis columns: ",
      paste(missing_annotations, collapse = ", "),
      call. = FALSE
    )
  }
  spec <- coord_spec(row$reduction, names(plot_data))
  lims <- square_limits(plot_data, spec$columns)
  best_invivo_idx <- plot_data$dataset == "invivo" & plot_data$point_type == "best"
  best_invivo_lims <- if (any(best_invivo_idx)) {
    square_limits(plot_data[best_invivo_idx, , drop = FALSE], spec$columns)
  } else {
    lims
  }
  best_size <- 1.25
  highlight_best_size <- 2.35
  highlight_background_size <- 1.15
  initial_size <- if (identical(row$variant, "Sampled")) 0.8 * best_size else 0.22
  n_invivo_best_missing <- sum(
    plot_data$dataset == "invivo" &
      plot_data$point_type == "best" &
      (is.na(plot_data$curve_class) | !nzchar(plot_data$curve_class))
  )
  n_invivo_best_missing_slope <- sum(
    plot_data$dataset == "invivo" &
      plot_data$point_type == "best" &
      !is.finite(plot_data$average_slope)
  )
  best_invivo_for_order <- prepare_invivo_best_class_data(
    prepare_plot_data(plot_data, spec$columns),
    class_order = NULL
  )
  class_order <- curve_class_order(best_invivo_for_order$curve_class_for_plot)
  slope_values <- suppressWarnings(as.numeric(best_invivo_for_order$average_slope))
  max_abs_slope <- suppressWarnings(max(abs(slope_values[is.finite(slope_values)]), na.rm = TRUE))
  if (!is.finite(max_abs_slope) || max_abs_slope <= 0) max_abs_slope <- 1
  slope_limits <- c(-max_abs_slope, max_abs_slope)
  monotone_increasing_slope_values <- suppressWarnings(as.numeric(
    best_invivo_for_order$average_slope[as.character(best_invivo_for_order$curve_class_for_plot) == "monotone_increasing"]
  ))
  monotone_increasing_slope_limits <- finite_range(monotone_increasing_slope_values, fallback = slope_limits)

  original_plot <- build_objective_plot(
    plot_data = plot_data,
    coord_names = spec$columns,
    axis_labels = spec$labels,
    lims = lims,
    initial_size = initial_size,
    best_size = best_size
  )
  curve_plot <- build_curve_class_plot(
    plot_data = plot_data,
    reduction = row$reduction,
    coord_names = spec$columns,
    axis_labels = spec$labels,
    lims = lims,
    initial_size = initial_size,
    best_size = best_size
  )
  invivo_best_curve_plot <- build_invivo_best_curve_class_plot(
    plot_data = plot_data,
    coord_names = spec$columns,
    axis_labels = spec$labels,
    lims = best_invivo_lims,
    best_size = best_size,
    class_order = class_order
  )
  class_highlight_plot <- build_curve_class_highlight_3x3(
    plot_data = plot_data,
    coord_names = spec$columns,
    lims = best_invivo_lims,
    class_order = class_order,
    highlight_size = highlight_best_size,
    background_size = highlight_background_size
  )
  slope_highlight_plot <- build_average_slope_highlight_3x3(
    plot_data = plot_data,
    coord_names = spec$columns,
    lims = best_invivo_lims,
    class_order = class_order,
    slope_limits = slope_limits,
    highlight_size = highlight_best_size,
    background_size = highlight_background_size
  )
  monotone_increasing_plot <- build_monotone_increasing_combined(
    plot_data = plot_data,
    coord_names = spec$columns,
    lims = best_invivo_lims,
    class_order = class_order,
    slope_limits = monotone_increasing_slope_limits,
    highlight_size = highlight_best_size,
    background_size = highlight_background_size
  )
  combined_plot <- original_plot + invivo_best_curve_plot + patchwork::plot_layout(ncol = 2, widths = c(1, 1))
  extended_combined_plot <- original_plot + invivo_best_curve_plot + class_highlight_plot + slope_highlight_plot +
    patchwork::plot_layout(ncol = 2, widths = c(1, 1), heights = c(1, 1))

  figure_dir <- file.path(out_dir, "figures", row$reduction, row$variant)
  curve_stem <- file.path(figure_dir, paste0(row$stub, "_curve_class"))
  invivo_best_curve_stem <- file.path(figure_dir, paste0(row$stub, "_invivo_best_curve_class"))
  class_highlight_stem <- file.path(figure_dir, paste0(row$stub, "_invivo_best_curve_class_highlight_3x3"))
  slope_highlight_stem <- file.path(figure_dir, paste0(row$stub, "_invivo_best_average_slope_highlight_3x3"))
  monotone_increasing_stem <- file.path(figure_dir, paste0(row$stub, "_monotone_increasing_curve_class_vs_average_slope_1x2"))
  combined_stem <- file.path(figure_dir, paste0(row$stub, "_original_vs_curve_class"))
  extended_combined_stem <- file.path(figure_dir, paste0(row$stub, "_original_vs_invivo_best_curve_class_highlight_slope_2x2"))
  plot_paths <- save_plot_pair(curve_plot, curve_stem)
  invivo_best_curve_paths <- save_plot_pair(invivo_best_curve_plot, invivo_best_curve_stem)
  class_highlight_paths <- save_plot_pair(class_highlight_plot, class_highlight_stem, width = 10.5, height = 10.5)
  slope_highlight_paths <- save_plot_pair(slope_highlight_plot, slope_highlight_stem, width = 10.5, height = 10.5)
  monotone_increasing_paths <- save_plot_pair(monotone_increasing_plot, monotone_increasing_stem, width = 12.8, height = 5.4)
  combined_paths <- save_plot_pair(combined_plot, combined_stem, width = 14.8, height = 6.5)
  extended_combined_paths <- save_plot_pair(extended_combined_plot, extended_combined_stem, width = 18, height = 18)

  data.frame(
    reduction = row$reduction,
    variant = row$variant,
    stub = row$stub,
    coordinate_table = row$coordinate_table,
    original_png = row$original_png,
    curve_png = unname(plot_paths[["png"]]),
    curve_pdf = unname(plot_paths[["pdf"]]),
    invivo_best_curve_png = unname(invivo_best_curve_paths[["png"]]),
    invivo_best_curve_pdf = unname(invivo_best_curve_paths[["pdf"]]),
    curve_class_highlight_3x3_png = unname(class_highlight_paths[["png"]]),
    curve_class_highlight_3x3_pdf = unname(class_highlight_paths[["pdf"]]),
    average_slope_highlight_3x3_png = unname(slope_highlight_paths[["png"]]),
    average_slope_highlight_3x3_pdf = unname(slope_highlight_paths[["pdf"]]),
    monotone_increasing_combined_png = unname(monotone_increasing_paths[["png"]]),
    monotone_increasing_combined_pdf = unname(monotone_increasing_paths[["pdf"]]),
    combined_png = unname(combined_paths[["png"]]),
    combined_pdf = unname(combined_paths[["pdf"]]),
    extended_combined_png = unname(extended_combined_paths[["png"]]),
    extended_combined_pdf = unname(extended_combined_paths[["pdf"]]),
    n_rows = nrow(plot_data),
    n_initial = sum(plot_data$point_type == "initial"),
    n_best_invivo = sum(plot_data$dataset == "invivo" & plot_data$point_type == "best"),
    n_best_invitro = sum(plot_data$dataset == "invitro" & plot_data$point_type == "best"),
    n_invivo_best_with_class = sum(plot_data$dataset == "invivo" & plot_data$point_type == "best" & !is.na(plot_data$curve_class) & nzchar(plot_data$curve_class)),
    n_invivo_best_missing_class = n_invivo_best_missing,
    n_invivo_best_with_average_slope = sum(plot_data$dataset == "invivo" & plot_data$point_type == "best" & is.finite(plot_data$average_slope)),
    n_invivo_best_missing_average_slope = n_invivo_best_missing_slope,
    average_slope_min = suppressWarnings(min(plot_data$average_slope[plot_data$dataset == "invivo" & plot_data$point_type == "best"], na.rm = TRUE)),
    average_slope_max = suppressWarnings(max(plot_data$average_slope[plot_data$dataset == "invivo" & plot_data$point_type == "best"], na.rm = TRUE)),
    monotone_increasing_n = sum(
      plot_data$dataset == "invivo" &
        plot_data$point_type == "best" &
        plot_data$curve_class == "monotone_increasing"
    ),
    monotone_increasing_average_slope_min = monotone_increasing_slope_limits[[1L]],
    monotone_increasing_average_slope_max = monotone_increasing_slope_limits[[2L]],
    average_slope_table = if ("average_slope_table" %in% names(row)) {
      row$average_slope_table
    } else {
      average_slope_table
    },
    stringsAsFactors = FALSE
  )
}

default_analysis_manifest <- function(out_dir) {
  file.path(out_dir, "tables", "pooled_embedding_curve_class_analysis_manifest.tsv")
}

main <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  argv <- bpf_parse_args(raw_args)
  if (bpf_as_bool(argv$help %||% argv$h, FALSE)) {
    usage()
    return(invisible(NULL))
  }

  require_package("ggplot2")
  require_package("ggnewscale")
  require_package("patchwork")
  require_package("scales")

  repo_root <- bpf_repo_root(SCRIPT_DIR)
  out_dir <- bpf_resolve_repo_path(
    argv$out_dir %||% file.path(bpf_combine_fixo2_eigen_attractor_result_dir(repo_root), "pooled_embedding_curve_class"),
    repo_root = repo_root,
    mustWork = FALSE
  )
  dry_run <- bpf_as_bool(argv$dry_run, FALSE)
  analysis_manifest <- bpf_resolve_repo_path(
    argv$analysis_manifest %||% default_analysis_manifest(out_dir),
    repo_root = repo_root,
    mustWork = TRUE
  )
  reductions <- normalize_reductions(argv$reductions)
  variants <- normalize_variants(argv$variants)
  embeddings <- read_tsv_plain(analysis_manifest)
  required_manifest_columns <- c(
    "reduction", "variant", "stub", "coordinate_table", "annotated_coordinate_table",
    "original_png", "average_slope_table", "class_table", "class_col"
  )
  missing_manifest_columns <- setdiff(required_manifest_columns, names(embeddings))
  if (length(missing_manifest_columns)) {
    stop(
      "Analysis manifest is missing required columns: ",
      paste(missing_manifest_columns, collapse = ", "),
      call. = FALSE
    )
  }
  embeddings <- embeddings[
    embeddings$reduction %in% reductions & embeddings$variant %in% variants,
    ,
    drop = FALSE
  ]
  if (!nrow(embeddings)) {
    stop("No matching analysis-prepared embedding rows in: ", analysis_manifest, call. = FALSE)
  }

  message("Analysis manifest: ", analysis_manifest)
  message("Output root: ", out_dir)
  message("Annotated embedding tables: ", nrow(embeddings))

  if (isTRUE(dry_run)) {
    print(
      embeddings[, c("reduction", "variant", "annotated_coordinate_table", "original_png")],
      row.names = FALSE
    )
    return(invisible(embeddings))
  }

  manifest_rows <- vector("list", nrow(embeddings))
  for (i in seq_len(nrow(embeddings))) {
    row <- embeddings[i, , drop = FALSE]
    message("[", i, "/", nrow(embeddings), "] ", row$reduction, " ", row$variant, " ", row$stub)
    manifest_rows[[i]] <- process_embedding(
      row,
      class_map = character(),
      slope_map = numeric(),
      average_slope_table = row$average_slope_table,
      out_dir = out_dir
    )
  }
  manifest <- do.call(rbind, manifest_rows)
  manifest$class_table <- embeddings$class_table
  manifest$class_col <- embeddings$class_col
  manifest$average_slope_table <- embeddings$average_slope_table
  manifest_path <- file.path(out_dir, "tables", "pooled_embedding_curve_class_manifest.tsv")
  write_tsv(manifest, manifest_path)
  message("Wrote manifest: ", manifest_path)
  invisible(manifest)
}

if (identical(environment(), globalenv())) {
  main()
}
