#!/usr/bin/env Rscript

# Five-initial-ploidy, full oxygen/time layouts for the signed-log net-growth
# extension. In vivo and in vitro are rendered separately while sharing one
# persisted color contract.

options(stringsAsFactors = FALSE, warn = 1)

f6ngfr_objects <- function(paths, context) {
  run <- f6ng_paths(paths)
  objects <- setNames(lapply(f6ft_family_levels(), function(family) {
    f6ng_read_panel(run, context, family)
  }), f6ft_family_levels())
  contract_path <- file.path(
    run$run_root, "net_growth_full_range_color_contract.tsv"
  )
  f6r_require_files(contract_path, "full-range net-growth color contract")
  contract <- f6r_read_tsv(contract_path)
  if (nrow(contract) != 1L ||
      !identical(as.character(contract$scale[[1L]]), "signed_pseudo_log10")) {
    stop("Malformed full-range net-growth color contract.")
  }
  limits <- as.numeric(contract[1L, c(
    "displayed_minimum_per_day", "displayed_maximum_per_day"
  )])
  if (any(!is.finite(limits)) || limits[[1L]] >= 0 || limits[[2L]] <= 0 ||
      abs(limits[[1L]] + limits[[2L]]) > 1e-12) {
    stop("Full-range signed-log color limits must be finite and symmetric.")
  }
  list(run = run, objects = objects, limits = limits, contract = contract)
}

f6ngfr_build <- function(
    objects, context, fill_limits, panel_mm = 38
) {
  f6r_require_packages(c("ggplot2", "scales"))
  stopifnot(
    context %in% c("in vivo", "in vitro"),
    identical(names(objects), f6ft_family_levels())
  )
  family <- f6ft_family_levels()
  initial <- rev(f6ft_initial_ploidy())
  p_values <- f6ft_p_values()
  context_color <- c("in vivo" = "#0072B2", "in vitro" = "#CC79A7")
  family_color <- c(C01 = "#C99700", C02 = "#6A3D9A")
  day_limits <- c(0, 10000)
  o2_limits <- if (identical(context, "in vivo")) c(0, 5) else c(0, 20)
  o2_breaks <- if (identical(context, "in vivo")) 0:5 else seq(0, 20, by = 5)

  gap <- .55
  cluster_gap <- 28
  row_gap <- 2
  row_top <- 24 + (seq_along(initial) - 1L) * (panel_mm + row_gap)
  panel_x <- 23 + (0:9) * (panel_mm + gap) +
    c(rep(0, 5), rep(cluster_gap - gap, 5))
  group_x <- panel_x[c(1L, 6L)]
  group_span <- 5 * panel_mm + 4 * gap
  width <- tail(panel_x, 1L) + panel_mm + 5
  bottom <- tail(row_top, 1L) + panel_mm
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
    text(
      label, x, y, w, h, size = 8.2, bold = TRUE, rotate = rotate,
      colour = "#FFFFFF"
    )
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

  fill_scale <- f6sb_growth_fill(fill_limits, growth_scale = "signed_log")
  for (cluster in seq_along(family)) {
    text(paste0(LETTERS[[cluster]], ". ", family[[cluster]], " net growth rate"),
         group_x[[cluster]] - 18, 1, group_span + 18, 6,
         9.5, TRUE, just = "left")
    text("Fixed p_misseg", group_x[[cluster]], 7,
         group_span, 5, 9, TRUE)
  }
  for (cluster in seq_along(family)) {
    columns <- if (cluster == 1L) 1:5 else 6:10
    strip(
      family[[cluster]], family_color[[family[[cluster]]]],
      panel_x[columns[[1L]]], 14, 5 * panel_mm + 4 * gap, 4
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
  grid_height <- length(initial) * panel_mm + (length(initial) - 1L) * row_gap
  for (cluster in seq_along(family)) {
    strip(context, context_color[[context]],
          group_x[[cluster]] - 18.5, row_top[[1L]], 3.8,
          grid_height, 90)
    text("Fixed oxygen (%)", group_x[[cluster]] - 14.7,
         row_top[[1L]], 5, grid_height, 8.8, rotate = 90)
  }

  for (row in seq_along(initial)) {
    for (cluster in seq_along(family)) {
      add(grid::rectGrob(gp = grid::gpar(
        fill = "#F2F2F2", col = "#BEBEBE", lwd = .4
      )), group_x[[cluster]] - 9.7, row_top[[row]], 4.2, panel_mm)
      text(paste0(initial[[row]], "N"), group_x[[cluster]] - 9.7,
           row_top[[row]], 4.2, panel_mm,
           size = 8.2, bold = TRUE, rotate = 90)
    }
    for (cluster in seq_along(family)) for (p in seq_along(p_values)) {
      data <- f6sb_extract(
        objects[[family[[cluster]]]], "mean_net_growth_rate",
        initial[[row]], p_values[[p]], day_limits, o2_limits
      )
      labels <- if (p == 1L) c("0", "") else if (p == length(p_values)) {
        c("", "10000")
      } else c("", "")
      plot <- ggplot2::ggplot(
        data, ggplot2::aes(x = day, y = O2_pct, fill = value)
      ) +
        ggplot2::geom_raster(interpolate = FALSE) + fill_scale +
        ggplot2::geom_hline(
          yintercept = .5, colour = "#7A7A7A",
          linetype = "dashed", linewidth = .28
        ) +
        ggplot2::scale_x_continuous(
          breaks = day_limits, labels = labels, expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
          breaks = o2_breaks, expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(
          xlim = day_limits, ylim = o2_limits, expand = FALSE
        ) +
        ggplot2::labs(x = NULL, y = NULL) + f6sb_tile_theme()
      column <- (cluster - 1L) * 5L + p
      tile(
        f6sb_parts(plot), panel_x[[column]], row_top[[row]],
        paste0(LETTERS[[cluster]], row, "_", p),
        show_y = column == 1L || column == 6L,
        show_x = row == length(initial)
      )
    }
  }
  for (cluster in seq_along(family)) {
    text("Experimental time (days)", group_x[[cluster]], bottom + 7.5,
         group_span, 5, 9)
  }

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
  legend_index <- which(legend$layout$name == "guide-box-bottom")
  if (!length(legend_index)) legend_index <- which(legend$layout$name == "guide-box")
  add(
    legend$grobs[[legend_index[[1L]]]], (width - 90) / 2,
    bottom + 13.5, 90, 12
  )
  list(
    plot = do.call(grid::grobTree, children), width = width / 25.4,
    height = height / 25.4, geometry = do.call(rbind, geometry),
    fill_limits = fill_limits
  )
}

f6ngfr_draw <- function(
    workspace_root = f6r_find_workspace_root(), context,
    supplement, filename
) {
  paths <- f6r_paths(workspace_root)
  bundle <- f6ngfr_objects(paths, context)
  layout <- f6ngfr_build(
    bundle$objects, context, bundle$limits
  )
  geometry <- layout$geometry
  stopifnot(
    nrow(geometry) == 50L,
    all(geometry$width_mm == geometry$height_mm),
    all(geometry$width_mm == 38),
    length(unique(geometry$x_mm)) == 10L,
    length(unique(geometry$y_mm)) == 5L,
    sum(startsWith(geometry$panel, "A")) == 25L,
    sum(startsWith(geometry$panel, "B")) == 25L
  )
  directory <- file.path(paths$root, "data", "Figures", supplement)
  rendered <- file.path(directory, "rendered")
  dir.create(rendered, recursive = TRUE, showWarnings = FALSE)
  base_path <- file.path(rendered, filename)
  pdf_path <- paste0(base_path, ".pdf")
  png_path <- paste0(base_path, ".png")

  ggplot2::ggsave(
    pdf_path, plot = layout$plot, width = layout$width,
    height = layout$height, units = "in",
    device = function(file, width, height, ...) grDevices::cairo_pdf(
      filename = file, width = width, height = height, family = "Helvetica", ...
    ),
    bg = "white", limitsize = FALSE
  )
  png_layout <- f6ngfr_build(
    bundle$objects, context, bundle$limits
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
  expected_o2 <- if (identical(context, "in vivo")) "0:5" else "0:20"
  validation <- f6g_render_hash_validation(
    output, published,
    file.path(directory, paste0(filename, "_render_validation.tsv")),
    data.frame(
      check = c(
        "fifty_panels", "all_plotting_regions_square_38mm",
        "five_initial_ploidies", "ten_columns", "twenty_five_A_twenty_five_B",
        "A_B_gap_at_least_25mm", "day_range_0_10000",
        "oxygen_range", "signed_log_shared_color_contract"
      ),
      observed = c(
        nrow(geometry) == 50L,
        all(geometry$width_mm == 38 & geometry$height_mm == 38),
        length(unique(geometry$y_mm)) == 5L,
        length(unique(geometry$x_mm)) == 10L,
        sum(startsWith(geometry$panel, "A")) == 25L &&
          sum(startsWith(geometry$panel, "B")) == 25L,
        min(geometry$x_mm[startsWith(geometry$panel, "B")]) -
          max(geometry$x_mm[startsWith(geometry$panel, "A")] +
                geometry$width_mm[startsWith(geometry$panel, "A")]) >= 25,
        all(vapply(bundle$objects, function(x) {
          identical(range(x$day_values), c(0L, 10000L))
        }, logical(1L))),
        paste(range(bundle$objects[[1L]]$o2_values), collapse = ":"),
        bundle$contract$scale[[1L]]
      ),
      expected = c(
        rep(TRUE, 7L), expected_o2, "signed_pseudo_log10"
      ), stringsAsFactors = FALSE
    )
  )
  invisible(list(
    output = output, published = published, validation = validation,
    geometry = geometry_path, fill_limits = layout$fill_limits
  ))
}
