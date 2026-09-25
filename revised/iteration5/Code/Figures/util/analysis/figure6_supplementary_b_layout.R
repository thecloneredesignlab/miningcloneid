#!/usr/bin/env Rscript

# Shared B-only layout for Supplementary Figures 6-14 through 6-16.
# Each heatmap panel has an exact 38 mm by 38 mm plotting region.

options(stringsAsFactors = FALSE, warn = 1)

f6sb_ploidy_fill <- function() ggplot2::scale_fill_gradientn(
  colours = c("#2166AC", "#FFFFBF", "#B2182B"), trans = "log10",
  limits = c(1, 7), breaks = c(1, 1.5, 2, 3, 4, 6),
  name = "Mean ploidy (log colors)", na.value = "#D9D9D9"
)

f6sb_growth_limits <- function(objects, day_limits, o2_limits) {
  values <- unlist(lapply(objects, function(object) {
    day <- which(object$day_values >= day_limits[[1L]] &
      object$day_values <= day_limits[[2L]])
    oxygen <- which(object$o2_values >= o2_limits[[1L]] - 1e-12 &
      object$o2_values <= o2_limits[[2L]] + 1e-12)
    as.vector(object$mean_net_growth_rate[, day, oxygen, , drop = FALSE])
  }), use.names = FALSE)
  values <- values[is.finite(values)]
  if (!length(values)) stop("No finite population net-growth values are available.")
  observed <- range(values)
  if (observed[[1L]] < 0 && observed[[2L]] > 0) {
    bound <- max(abs(observed))
    c(-bound, bound)
  } else observed
}

f6sb_growth_fill <- function(
    limits, growth_scale = c("linear", "signed_log")
) {
  growth_scale <- match.arg(growth_scale)
  if (limits[[1L]] < 0 && limits[[2L]] > 0) {
    transformation <- if (identical(growth_scale, "signed_log")) {
      scales::pseudo_log_trans(sigma = 0.01, base = 10)
    } else {
      "identity"
    }
    scale_breaks <- if (identical(growth_scale, "signed_log")) {
      c(-1, -0.1, 0, 0.1, 1)
    } else {
      scales::breaks_pretty(n = 5)
    }
    scale_labels <- if (identical(growth_scale, "signed_log")) {
      c("-1.0", "-0.1", "0", "0.1", "1.0")
    } else {
      ggplot2::waiver()
    }
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = limits, trans = transformation,
      breaks = scale_breaks, labels = scale_labels,
      name = expression("Population net-live growth rate (day"^{-1}*")"),
      na.value = "#D9D9D9"
    )
  } else {
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "#FFFFBF", "#B2182B"), limits = limits,
      breaks = scales::breaks_pretty(n = 5),
      name = expression("Population net-live growth rate (day"^{-1}*")"),
      na.value = "#D9D9D9"
    )
  }
}

f6sb_extract <- function(
    object, metric, initial_ploidy, p_misseg, day_limits, o2_limits
) {
  values <- object[[metric]]
  expected <- c(
    length(object$initial_ploidy), length(object$day_values),
    length(object$o2_values), length(object$p_misseg)
  )
  if (is.null(values) || !identical(as.integer(dim(values)), as.integer(expected))) {
    stop("Panel metric has unexpected dimensions: ", metric)
  }
  initial <- match(initial_ploidy, object$initial_ploidy)
  p <- match(sprintf("%.12f", p_misseg), sprintf("%.12f", object$p_misseg))
  day <- which(object$day_values >= day_limits[[1L]] &
    object$day_values <= day_limits[[2L]])
  oxygen <- which(object$o2_values >= o2_limits[[1L]] - 1e-12 &
    object$o2_values <= o2_limits[[2L]] + 1e-12)
  if (is.na(initial) || is.na(p) || !length(day) || !length(oxygen)) {
    stop("Requested supplementary Figure 6 slice is absent from its panel object.")
  }
  matrix_values <- values[initial, day, oxygen, p, drop = TRUE]
  if (!is.matrix(matrix_values)) {
    matrix_values <- matrix(matrix_values, nrow = length(day), ncol = length(oxygen))
  }
  data.frame(
    day = rep(object$day_values[day], times = length(oxygen)),
    O2_pct = rep(object$o2_values[oxygen], each = length(day)),
    value = as.vector(matrix_values), stringsAsFactors = FALSE
  )
}

f6sb_tile_theme <- function() {
  ggplot2::theme_classic(base_size = 9, base_family = "Helvetica") +
    ggplot2::theme(
      legend.position = "none", aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(size = 8.2, colour = "#333333"),
      axis.text.y = ggplot2::element_text(size = 8.2, colour = "#333333"),
      axis.ticks = ggplot2::element_line(linewidth = .25),
      panel.border = ggplot2::element_rect(
        fill = NA, colour = "#555555", linewidth = .22
      ),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
}

f6sb_parts <- function(plot) {
  grob <- ggplot2::ggplotGrob(plot)
  get <- function(name) grob$grobs[[match(name, grob$layout$name)]]
  list(panel = get("panel"), x = get("axis-b"), y = get("axis-l"))
}

f6sb_build <- function(
    objects, metric = c("mean_ploidy", "mean_net_growth_rate"),
    day_limits, o2_limits, title, panel_mm = 38,
    growth_scale = c("linear", "signed_log")
) {
  metric <- match.arg(metric)
  growth_scale <- match.arg(growth_scale)
  if (identical(metric, "mean_ploidy") && !identical(growth_scale, "linear")) {
    stop("Signed-log colors are defined only for population net-growth panels.")
  }
  f6r_require_packages(c("ggplot2", "scales"))
  context <- c("in vivo", "in vivo", "in vitro", "in vitro")
  initial <- c(4, 2, 4, 2)
  family <- c("C01", "C02")
  context_color <- c("in vivo" = "#0072B2", "in vitro" = "#CC79A7")
  family_color <- c(C01 = "#C99700", C02 = "#6A3D9A")
  p_values <- f6ft_p_values()
  if (length(p_values) != 5L) stop("B-only layout requires five fixed p_misseg values.")
  gap <- .55
  cluster_gap <- 10
  row_gap <- 2
  row_top <- 24 + c(
    0, panel_mm + row_gap,
    2 * panel_mm + row_gap + 5,
    3 * panel_mm + 2 * row_gap + 5
  )
  panel_x <- 23 + (0:9) * (panel_mm + gap) +
    c(rep(0, 5), rep(cluster_gap - gap, 5))
  width <- tail(panel_x, 1) + panel_mm + 5
  bottom <- tail(row_top, 1) + panel_mm
  height <- bottom + 29
  children <- list()
  geometry <- list()

  add <- function(grob, x, y, w, h, clip = "off") {
    children[[length(children) + 1L]] <<- grid::grobTree(
      grob,
      vp = grid::viewport(
        x = grid::unit(x, "mm"), y = grid::unit(height - y, "mm"),
        width = grid::unit(w, "mm"), height = grid::unit(h, "mm"),
        just = c("left", "top"), clip = clip
      )
    )
  }
  text <- function(
      label, x, y, w, h, size = 7, bold = FALSE, rotate = 0,
      just = "centre", colour = "#222222"
  ) {
    add(grid::textGrob(
      label, x = if (just == "left") 0 else .5, just = just, rot = rotate,
      gp = grid::gpar(
        fontfamily = "Helvetica", fontsize = size, col = colour,
        fontface = if (bold) "bold" else "plain"
      )
    ), x, y, w, h)
  }
  strip <- function(label, colour, x, y, w, h, rotate = 0) {
    add(grid::rectGrob(gp = grid::gpar(
      fill = colour, col = "#BEBEBE", lwd = .4
    )), x, y, w, h)
    text(label, x, y, w, h, size = 8.2, bold = TRUE, rotate = rotate,
      colour = "#FFFFFF")
  }
  tile <- function(parts, x, y, label, show_y, show_x) {
    add(parts$panel, x, y, panel_mm, panel_mm, "on")
    if (show_y) add(parts$y, x - 6.5, y, 6.5, panel_mm)
    if (show_x) add(parts$x, x, y + panel_mm, panel_mm, 7)
    geometry[[length(geometry) + 1L]] <<- data.frame(
      panel = label, x_mm = x, y_mm = y,
      width_mm = panel_mm, height_mm = panel_mm,
      stringsAsFactors = FALSE
    )
  }

  fill_limits <- if (identical(metric, "mean_ploidy")) {
    c(1, 7)
  } else {
    f6sb_growth_limits(objects, day_limits, o2_limits)
  }
  fill_scale <- if (identical(metric, "mean_ploidy")) {
    f6sb_ploidy_fill()
  } else {
    f6sb_growth_fill(fill_limits, growth_scale = growth_scale)
  }
  text(title, 5, 1, width - 10, 6, 9.5, TRUE, just = "left")
  text(expression("Fixed " * p[misseg]), panel_x[[1L]], 7,
    width - panel_x[[1L]] - 5, 5, 9, TRUE)
  for (cluster in seq_along(family)) {
    columns <- if (cluster == 1L) 1:5 else 6:10
    strip(
      family[[cluster]], family_color[[family[[cluster]]]],
      panel_x[columns[[1L]]], 14,
      5 * panel_mm + 4 * gap, 4
    )
    for (p in seq_along(columns)) {
      add(grid::rectGrob(gp = grid::gpar(
        fill = "#F2F2F2", col = "#BEBEBE", lwd = .4
      )), panel_x[columns[[p]]], 18.3, panel_mm, 4.8)
      text(
        f6ft_format_p(p_values[[p]]), panel_x[columns[[p]]], 18.3,
        panel_mm, 4.8, size = 8.2, bold = TRUE
      )
    }
  }
  for (group in 1:2) {
    first_row <- if (group == 1L) 1L else 3L
    group_height <- 2 * panel_mm + row_gap
    strip(
      context[[first_row]], context_color[[context[[first_row]]]],
      5, row_top[[first_row]], 3.8, group_height, 90
    )
    text(
      "Fixed oxygen (%)", 9.2, row_top[[first_row]], 4,
      group_height, 8.8, rotate = 90
    )
  }

  for (row in 1:4) {
    add(grid::rectGrob(gp = grid::gpar(
      fill = "#F2F2F2", col = "#BEBEBE", lwd = .4
    )), 13.3, row_top[[row]], 4.2, panel_mm)
    text(
      paste0(initial[[row]], "N"), 13.3, row_top[[row]],
      4.2, panel_mm, size = 8.2, bold = TRUE, rotate = 90
    )
    for (cluster in seq_along(family)) for (p in seq_along(p_values)) {
      key <- paste(context[[row]], family[[cluster]], sep = "|")
      data <- f6sb_extract(
        objects[[key]], metric, initial[[row]], p_values[[p]],
        day_limits, o2_limits
      )
      last_day <- day_limits[[2L]]
      labels <- if (p == 1L) c("0", "") else if (p == length(p_values)) {
        c("", as.character(last_day))
      } else c("", "")
      plot <- ggplot2::ggplot(
        data, ggplot2::aes(x = day, y = O2_pct, fill = value)
      ) +
        ggplot2::geom_raster(interpolate = FALSE) + fill_scale +
        ggplot2::scale_x_continuous(
          breaks = c(day_limits[[1L]], day_limits[[2L]]), labels = labels,
          expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
          breaks = seq(o2_limits[[1L]], o2_limits[[2L]], by = 1),
          expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(
          xlim = day_limits, ylim = o2_limits, expand = FALSE
        ) +
        ggplot2::labs(x = NULL, y = NULL) + f6sb_tile_theme()
      if (row <= 2L) plot <- plot + ggplot2::geom_hline(
        yintercept = .5, colour = "#7A7A7A",
        linetype = "dashed", linewidth = .28
      )
      column <- (cluster - 1L) * 5L + p
      tile(
        f6sb_parts(plot), panel_x[[column]], row_top[[row]],
        paste0("row", row, "_column", column),
        show_y = column == 1L, show_x = row == 4L
      )
    }
  }
  text(
    "Experimental time (days)", panel_x[[1L]], bottom + 7.5,
    width - panel_x[[1L]] - 5, 5, 9
  )

  legend_data <- data.frame(
    x = seq(fill_limits[[1L]], fill_limits[[2L]], length.out = 101L), y = 1
  )
  legend_plot <- ggplot2::ggplot(
    legend_data, ggplot2::aes(x = x, y = y, fill = x)
  ) + ggplot2::geom_tile() + fill_scale +
    ggplot2::theme_void(base_family = "Helvetica") +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7)
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      title.position = "top", title.hjust = .5,
      barwidth = grid::unit(64, "mm"), barheight = grid::unit(3, "mm")
    ))
  legend <- ggplot2::ggplotGrob(legend_plot)
  index <- which(legend$layout$name == "guide-box-bottom")
  if (!length(index)) index <- which(legend$layout$name == "guide-box")
  add(legend$grobs[[index[[1L]]]], (width - 90) / 2, bottom + 13.5, 90, 12)
  list(
    plot = do.call(grid::grobTree, children), width = width / 25.4,
    height = height / 25.4, geometry = do.call(rbind, geometry),
    fill_limits = fill_limits
  )
}

f6sb_objects_full_range <- function(paths) {
  override <- trimws(Sys.getenv("FIGURE6_FULL_RANGE_SOURCE_RUN_ROOT", ""))
  run <- if (!nzchar(override)) {
    f6g_paths(paths)
  } else {
    run_root <- normalizePath(override, mustWork = TRUE)
    allowed_root <- normalizePath(
      file.path(paths$root, "data", "Figures"), mustWork = TRUE
    )
    if (!startsWith(run_root, paste0(allowed_root, .Platform$file.sep))) {
      stop("Full-range source override must remain inside iteration4/data/Figures.")
    }
    list(run_id = basename(run_root), run_root = run_root)
  }
  list(
    "in vivo|C01" = f6g_read_panel(run, "in vivo", "continuous", "C01"),
    "in vivo|C02" = f6g_read_panel(run, "in vivo", "continuous", "C02"),
    "in vitro|C01" = f6g_read_panel(run, "in vitro", "passage", "C01"),
    "in vitro|C02" = f6g_read_panel(run, "in vitro", "passage", "C02")
  )
}

f6sb_objects_net_growth <- function(paths) {
  run <- f6ng_paths(paths)
  list(
    "in vivo|C01" = f6ng_read_panel(run, "in vivo", "C01"),
    "in vivo|C02" = f6ng_read_panel(run, "in vivo", "C02"),
    "in vitro|C01" = f6ng_read_panel(run, "in vitro", "C01"),
    "in vitro|C02" = f6ng_read_panel(run, "in vitro", "C02")
  )
}

f6sb_draw <- function(
    workspace_root = f6r_find_workspace_root(), supplement,
    filename, metric, day_limits, o2_limits, title,
    growth_scale = c("linear", "signed_log")
) {
  growth_scale <- match.arg(growth_scale)
  paths <- f6r_paths(workspace_root)
  objects <- if (identical(metric, "mean_ploidy")) {
    f6sb_objects_full_range(paths)
  } else {
    f6sb_objects_net_growth(paths)
  }
  layout <- f6sb_build(
    objects, metric = metric, day_limits = day_limits,
    o2_limits = o2_limits, title = title, growth_scale = growth_scale
  )
  geometry <- layout$geometry
  stopifnot(
    nrow(geometry) == 40L,
    all(geometry$width_mm == geometry$height_mm),
    length(unique(geometry$width_mm)) == 1L,
    length(unique(geometry$x_mm)) == 10L,
    length(unique(geometry$y_mm)) == 4L
  )
  directory <- file.path(paths$root, "data", "Figures", supplement)
  rendered <- file.path(directory, "rendered")
  dir.create(rendered, recursive = TRUE, showWarnings = FALSE)
  base_path <- file.path(rendered, filename)
  pdf_path <- paste0(base_path, ".pdf")
  png_path <- paste0(base_path, ".png")

  # grid grobs can acquire device-specific viewport state while being drawn.
  # Build a fresh layout for the second device so PDF and PNG cannot alter one
  # another's strips, panels, or clipping regions.
  ggplot2::ggsave(
    pdf_path, plot = layout$plot, width = layout$width,
    height = layout$height, units = "in",
    device = function(file, width, height, ...) grDevices::cairo_pdf(
      filename = file, width = width, height = height, family = "Helvetica", ...
    ),
    bg = "white", limitsize = FALSE
  )
  png_layout <- f6sb_build(
    objects, metric = metric, day_limits = day_limits,
    o2_limits = o2_limits, title = title, growth_scale = growth_scale
  )
  stopifnot(
    identical(layout$geometry, png_layout$geometry),
    identical(layout$fill_limits, png_layout$fill_limits)
  )
  ggplot2::ggsave(
    png_path, plot = png_layout$plot, width = png_layout$width,
    height = png_layout$height, units = "in", dpi = 300,
    device = "png", type = "cairo", bg = "white", limitsize = FALSE
  )
  output <- c(
    png = normalizePath(png_path, mustWork = TRUE),
    pdf = normalizePath(pdf_path, mustWork = TRUE)
  )
  published <- f6g_publish_plot(output, paths, filename)
  geometry_path <- f6ft_atomic_write_tsv(
    geometry, file.path(directory, paste0(filename, "_panel_geometry.tsv"))
  )
  validation <- f6g_render_hash_validation(
    output, published,
    file.path(directory, paste0(filename, "_render_validation.tsv")),
    data.frame(
      check = c(
        "forty_panels", "all_plotting_regions_square_38mm",
        "four_rows", "ten_columns", "requested_day_range",
        "requested_oxygen_range", "requested_growth_color_scale"
      ),
      observed = c(
        nrow(geometry) == 40L,
        all(geometry$width_mm == 38 & geometry$height_mm == 38),
        length(unique(geometry$y_mm)) == 4L,
        length(unique(geometry$x_mm)) == 10L,
        all(vapply(objects, function(x) all(day_limits %in% x$day_values), logical(1L))),
        all(vapply(objects, function(x) {
          min(x$o2_values) <= o2_limits[[1L]] + 1e-12 &&
            max(x$o2_values) >= o2_limits[[2L]] - 1e-12
        }, logical(1L))),
        growth_scale
      ),
      expected = c(rep(TRUE, 6L), growth_scale), stringsAsFactors = FALSE
    )
  )
  invisible(list(
    output = output, published = published, validation = validation,
    geometry = geometry_path, fill_limits = layout$fill_limits
  ))
}
