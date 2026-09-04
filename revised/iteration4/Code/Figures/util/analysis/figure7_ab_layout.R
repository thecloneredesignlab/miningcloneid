# Figure 7 publication layout: physical panel geometry, not outer plot boxes,
# is shared between the four steady-state panels and forty finite-time panels.
f7ab_fill <- function() ggplot2::scale_fill_gradientn(
  colours = c("#2166AC", "#FFFFBF", "#B2182B"), trans = "log10",
  limits = c(1, 7), breaks = c(1, 1.5, 2, 3, 4, 6),
  name = "Mean ploidy (log colors)", na.value = "#D9D9D9")

f7ab_surface_base <- function(paths, oxygen = c(0, 2)) {
  plot <- f7x_main_surface_plot(paths)
  plot$layers <- lapply(plot$layers, function(layer) {
    copied <- ggplot2::ggproto(NULL, layer)
    mapping <- layer$mapping
    old_x <- mapping$x
    mapping$x <- mapping$y
    mapping$y <- old_x
    copied$mapping <- mapping
    copied
  })
  suppressMessages(plot + f7ab_fill() +
    ggplot2::scale_x_continuous(
      breaks = log10(c(0.005, 0.01, 0.05, 0.1, 0.5)),
      labels = c("0.005", "0.01", "0.05", "0.10", "0.50")) +
    ggplot2::scale_y_continuous(breaks = if (oxygen[2] == 2) seq(0, 2, .5) else 0:5) +
    ggplot2::coord_cartesian(xlim = log10(c(.005, .5)), ylim = oxygen, expand = FALSE) +
    ggplot2::labs(x = expression(p[misseg]), y = "Fixed oxygen (%)"))
}

f7ab_tile_theme <- function() ggplot2::theme_classic(base_size = 7.7, base_family = "Helvetica") +
  ggplot2::theme(legend.position = "none", aspect.ratio = 1,
    axis.text = ggplot2::element_text(size = 5.7, colour = "#333333"),
    axis.ticks = ggplot2::element_line(linewidth = .25),
    panel.border = ggplot2::element_rect(fill = NA, colour = "#555555", linewidth = .22),
    plot.margin = ggplot2::margin(0, 0, 0, 0))

f7ab_parts <- function(plot) {
  grob <- ggplot2::ggplotGrob(plot)
  get <- function(name) grob$grobs[[match(name, grob$layout$name)]]
  list(panel = get("panel"), x = get("axis-b"), y = get("axis-l"))
}

f7ab_build <- function(paths, objects, panel_mm = 38) {
  base <- f7ab_surface_base(paths)
  context <- c("in vivo", "in vivo", "in vitro", "in vitro")
  family <- c("C01", "C02", "C01", "C02")
  initial <- c(4, 2, 4, 2)
  context_color <- c("in vivo" = "#0072B2", "in vitro" = "#CC79A7")
  family_color <- c(C01 = "#C99700", C02 = "#6A3D9A")
  gap <- .55
  row_gap <- 2.0
  row_top <- 24 + c(0, panel_mm + row_gap, 2 * panel_mm + row_gap + 5,
                    3 * panel_mm + 2 * row_gap + 5)
  a_x <- 23
  b_x <- a_x + panel_mm + 27
  b_columns <- b_x + (0:9) * (panel_mm + gap) + c(rep(0, 5), rep(5 - gap, 5))
  width <- tail(b_columns, 1) + panel_mm + 5
  bottom <- tail(row_top, 1) + panel_mm
  height <- bottom + 37
  children <- list()
  geometry <- list()
  add <- function(grob, x, y, w, h, clip = "off") {
    children[[length(children) + 1L]] <<- grid::grobTree(grob,
      vp = grid::viewport(x = grid::unit(x, "mm"), y = grid::unit(height - y, "mm"),
        width = grid::unit(w, "mm"), height = grid::unit(h, "mm"),
        just = c("left", "top"), clip = clip))
  }
  text <- function(label, x, y, w, h, size = 7, bold = FALSE, rotate = 0, just = "centre") {
    add(grid::textGrob(label, x = if (just == "left") 0 else .5,
      just = just, rot = rotate,
      gp = grid::gpar(fontfamily = "Helvetica", fontsize = size,
                     fontface = if (bold) "bold" else "plain")), x, y, w, h)
  }
  strip <- function(label, color, x, y, w, h, rotate = 0) {
    add(grid::rectGrob(gp = grid::gpar(fill = color, col = "#BEBEBE", lwd = .4)), x, y, w, h)
    text(label, x, y, w, h, size = 7.2, bold = TRUE, rotate = rotate)
  }
  tile <- function(parts, x, y, label, show_y = FALSE, show_x = FALSE) {
    add(parts$panel, x, y, panel_mm, panel_mm, "on")
    if (show_y) add(parts$y, x - 5.5, y, 5.5, panel_mm)
    if (show_x) add(parts$x, x, y + panel_mm, panel_mm, 4)
    geometry[[length(geometry) + 1L]] <<- data.frame(
      panel = label, x_mm = x, y_mm = y, width_mm = panel_mm, height_mm = panel_mm)
  }
  text("A. Steady-state ploidy", 5, 1, a_x + panel_mm - 5, 6, 9.5, TRUE, just = "left")
  text("B. Finite-time ploidy", b_x - 18, 1, width - b_x, 6, 9.5, TRUE, just = "left")
  text("Fixed p_misseg", b_x, 7, width - b_x - 5, 5, 8, TRUE)
  for (group in 1:2) {
    start <- if (group == 1) 1L else 3L
    context_height <- 2 * panel_mm + row_gap
    strip(context[start], context_color[[context[start]]], 5, row_top[start], 3.8, context_height, 90)
    text("Fixed oxygen (%)", 12.8, row_top[start], 4.2, context_height, 7.5, rotate = 90)
    strip(context[start], context_color[[context[start]]], b_x - 18.5, row_top[start], 3.8, context_height, 90)
    text("Fixed oxygen (%)", b_x - 14.7, row_top[start], 5, context_height, 7.5, rotate = 90)
  }
  for (column_group in 1:2) {
    index <- if (column_group == 1) 1:5 else 6:10
    strip(c("C01", "C02")[column_group], family_color[column_group],
      b_columns[index[1]], 14, 5 * panel_mm + 4 * gap, 4)
    for (j in seq_along(index)) strip(f7ft_format_p(f7ft_p_values()[j]), "#F2F2F2",
      b_columns[index[j]], 18.3, panel_mm, 4.8)
  }
  for (row in 1:4) {
    a <- base
    a$layers <- lapply(base$layers, function(layer) {
      copy <- ggplot2::ggproto(NULL, layer)
      d <- layer$data
      copy$data <- d[as.character(d$model_context) == context[row] &
                      as.character(d$display_label) == family[row], , drop = FALSE]
      copy
    })
    a <- a + ggplot2::facet_null() + ggplot2::labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) + f7ab_tile_theme()
    tile(f7ab_parts(a), a_x, row_top[row], paste0("A", row), TRUE, row == 4)
    strip(family[row], family_color[[family[row]]], 8.8, row_top[row], 4, panel_mm, 90)
    strip(paste0(initial[row], "N"), "#F2F2F2", b_x - 9.7, row_top[row], 4.2, panel_mm, 90)
    for (cluster in 1:2) {
      key <- if (row <= 2) c("C", "D")[cluster] else c("E", "F")[cluster]
      data <- f7ft_panel_long_subset(objects[[key]], initial[row], c(0, 2), c(0, 1000), 1L)
      for (p in 1:5) {
        z <- data[abs(data$p_misseg - f7ft_p_values()[p]) < 1e-12, ]
        b <- ggplot2::ggplot(z, ggplot2::aes(day, O2_pct, fill = mean_ploidy)) +
          ggplot2::geom_raster(interpolate = FALSE) + f7ab_fill() +
          ggplot2::scale_x_continuous(breaks = c(0, 500, 1000),
            labels = c(if (p == 1) "0" else "", "500", "1000"), expand = c(0, 0)) +
          ggplot2::scale_y_continuous(breaks = seq(0, 2, .5), expand = c(0, 0)) +
          ggplot2::coord_cartesian(xlim = c(0, 1000), ylim = c(0, 2), expand = FALSE) +
          ggplot2::labs(x = NULL, y = NULL) + f7ab_tile_theme()
        if (row <= 2) b <- b + ggplot2::geom_hline(yintercept = .5,
          colour = "#7A7A7A", linetype = "dashed", linewidth = .28)
        col <- (cluster - 1L) * 5L + p
        tile(f7ab_parts(b), b_columns[col], row_top[row],
             paste0("B", row, "_", col), col == 1, row == 4)
      }
    }
  }
  text(expression(p[misseg]), a_x, bottom + 4, panel_mm, 5, 8)
  text("Experimental time (days)", b_x, bottom + 4, width - b_x - 5, 5, 8)
  legend_plot <- ggplot2::ggplot(data.frame(x = 1:7, y = 1, p = 1:7),
    ggplot2::aes(x, y, fill = p)) + ggplot2::geom_tile() + f7ab_fill() +
    ggplot2::theme_void(base_family = "Helvetica") +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(size = 8),
                   legend.text = ggplot2::element_text(size = 7)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
      title.hjust = .5, barwidth = grid::unit(58, "mm"), barheight = grid::unit(3, "mm")))
  legend <- ggplot2::ggplotGrob(legend_plot)
  index <- which(legend$layout$name == "guide-box-bottom")
  if (!length(index)) index <- which(legend$layout$name == "guide-box")
  add(legend$grobs[[index[1]]], (width - 80) / 2, bottom + 10, 80, 12)
  text("A overlays: black solid = fitted p_misseg mean; black dashed = population-averaged p_mis(N,O2) mean; white dotted boundary / hatching = weak spectral gap.",
    5, bottom + 23, width - 10, 4, 6.8)
  text("Arithmetic means across the same 50 q10 optimizer endpoints. A: steady state. B: continuous in vivo; nearest-target segment selection in vitro. Gray dashed: oxygen 0.5%.",
    5, bottom + 27, width - 10, 4, 6.8)
  add(grid::rectGrob(gp = grid::gpar(fill = "#D9D9D9", col = "#888888", lwd = .4)),
      width / 2 - 99, bottom + 32, 3, 3)
  text("Gray fill: passage protocol infeasible for at least one of 50 endpoints; no survivor-only mean is shown.",
    width / 2 - 94, bottom + 31.5, 205, 4, 6.8, just = "left")
  list(plot = do.call(grid::grobTree, children), width = width / 25.4, height = height / 25.4,
       geometry = do.call(rbind, geometry))
}

f7ft_draw_main <- function(workspace_root = f7r_find_workspace_root()) {
  paths <- f7r_paths(workspace_root)
  run_paths <- f7g_paths(paths)
  objects <- list(C = f7g_read_panel(run_paths, "in vivo", "continuous", "C01"),
    D = f7g_read_panel(run_paths, "in vivo", "continuous", "C02"),
    E = f7g_read_panel(run_paths, "in vitro", "passage", "C01"),
    F = f7g_read_panel(run_paths, "in vitro", "passage", "C02"))
  layout <- f7ab_build(paths, objects)
  geometry <- layout$geometry
  stopifnot(nrow(geometry) == 44L, all(geometry$width_mm == geometry$height_mm),
    length(unique(geometry$width_mm)) == 1L, length(unique(geometry$y_mm)) == 4L)
  dir.create(run_paths$rendered, recursive = TRUE, showWarnings = FALSE)
  output <- f7r_save_plot(layout$plot, file.path(run_paths$rendered, "assembled_fig7"),
    width = layout$width, height = layout$height, dpi = 300)
  published <- f7g_publish_plot(output, paths, "assembled_fig7")
  f7ft_atomic_write_tsv(geometry, file.path(run_paths$run_root, "figure7_ab_panel_geometry.tsv"))
  validation <- f7g_render_hash_validation(output, published,
    file.path(run_paths$run_root, "figure7_full_range_render_validation.tsv"),
    data.frame(check = c("same_square_panel_size", "four_aligned_rows", "shared_color_limits_1_7"),
      observed = c(length(unique(geometry$width_mm)) == 1L && all(geometry$width_mm == geometry$height_mm),
        length(unique(geometry$y_mm)) == 4L,
        isTRUE(all.equal(f7ab_fill()$get_transformation()$inverse(f7ab_fill()$limits),
                         c(1, 7)))), expected = TRUE))
  invisible(list(output = output, published = published, validation = validation))
}

f7ft_draw_supp7_8 <- function(workspace_root = f7r_find_workspace_root()) {
  paths <- f7r_paths(workspace_root)
  plot <- f7ab_surface_base(paths, c(0, 5)) + ggplot2::labs(
    title = "Steady-state ploidy response: full oxygen range", subtitle = "Arithmetic means across 50 q10 optimizer endpoints") +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  plot$facet$params$free$x <- TRUE
  plot <- plot + ggh4x::facetted_pos_scales(x = list(
    display_label == "C01" ~ ggplot2::scale_x_continuous(
      breaks = log10(c(.005, .01, .05, .1, .5)),
      labels = c("0.005", "0.01", "0.05", "0.10", "0.50")),
    TRUE ~ ggplot2::scale_x_continuous(
      breaks = log10(c(.005, .01, .05, .1, .5)),
      labels = c("", "0.01", "0.05", "0.10", "0.50"))
  ))
  directory <- file.path(paths$root, "data", "Figures", "Supp_Figure7_8")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  output <- f7r_save_plot(plot, file.path(directory, "supp_fig7-8_steady_state_full_oxygen_range"),
    width = 18.2, height = 6.2)
  published <- f7g_publish_plot(output, paths, "supp_fig7-8_steady_state_full_oxygen_range")
  f7g_render_hash_validation(output, published,
    file.path(directory, "supp_fig7-8_steady_state_full_oxygen_range_render_validation.tsv"))
}
