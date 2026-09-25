# Figure 6 publication layout: physical panel geometry, not outer plot boxes,
# is shared between the four steady-state panels and the two finite-time panels.
# The same finite-time layout is reused for the C02-only Supplementary Figure 6-17.
f6ab_fill <- function() ggplot2::scale_fill_gradientn(
  colours = c("#2166AC", "#FFFFBF", "#B2182B"), trans = "log10",
  limits = c(1, 7), breaks = c(1, 1.5, 2, 3, 4, 6),
  name = "Mean ploidy (log colors)", na.value = "#D9D9D9")

f6ab_surface_base <- function(paths, oxygen = c(0, 2)) {
  plot <- f6x_main_surface_plot(paths)
  plot$layers <- lapply(plot$layers, function(layer) {
    copied <- ggplot2::ggproto(NULL, layer)
    mapping <- layer$mapping
    old_x <- mapping$x
    mapping$x <- mapping$y
    mapping$y <- old_x
    copied$mapping <- mapping
    copied
  })
  suppressMessages(plot + f6ab_fill() +
    ggplot2::scale_x_continuous(
      breaks = log10(c(0.005, 0.01, 0.05, 0.1, 0.5)),
      labels = c("0.005", "0.01", "0.05", "0.10", "0.50")) +
    ggplot2::scale_y_continuous(breaks = if (oxygen[2] == 2) seq(0, 2, .5) else 0:5) +
    ggplot2::coord_cartesian(xlim = log10(c(.005, .5)), ylim = oxygen, expand = FALSE) +
    ggplot2::labs(x = expression(p[misseg]), y = "Fixed oxygen (%)"))
}

f6ab_tile_theme <- function(x_angle = 0) ggplot2::theme_classic(base_size = 9, base_family = "Helvetica") +
  ggplot2::theme(legend.position = "none", aspect.ratio = 1,
    axis.text.x = ggplot2::element_text(
      size = 8.2, colour = "#333333", angle = x_angle,
      hjust = if (x_angle == 0) 0.5 else 1, vjust = if (x_angle == 0) 1 else 1
    ),
    axis.text.y = ggplot2::element_text(size = 8.2, colour = "#333333"),
    axis.ticks = ggplot2::element_line(linewidth = .25),
    panel.border = ggplot2::element_rect(fill = NA, colour = "#555555", linewidth = .22),
    plot.margin = ggplot2::margin(0, 0, 0, 0))

f6ab_parts <- function(plot) {
  grob <- ggplot2::ggplotGrob(plot)
  get <- function(name) grob$grobs[[match(name, grob$layout$name)]]
  list(panel = get("panel"), x = get("axis-b"), y = get("axis-l"))
}

f6ab_growth_limits <- function(paths) {
  run <- f6ng_paths(paths)
  contract_path <- file.path(run$run_root, "net_growth_full_range_color_contract.tsv")
  f6r_require_files(contract_path, "full-range net-growth color contract")
  contract <- f6r_read_tsv(contract_path)
  if (nrow(contract) != 1L ||
      !identical(as.character(contract$scale[[1L]]), "signed_pseudo_log10")) {
    stop("Malformed full-range net-growth color contract.")
  }
  limits <- as.numeric(contract[1L, c(
    "displayed_minimum_per_day", "displayed_maximum_per_day")])
  if (any(!is.finite(limits)) || limits[[1L]] >= 0 || limits[[2L]] <= 0 ||
      abs(sum(limits)) > 1e-12) {
    stop("Full-range net-growth color limits must be finite and symmetric.")
  }
  limits
}

f6ab_panel_coverage <- function(objects, mode = c("all", "oxygen", "time")) {
  mode <- match.arg(mode)
  all(vapply(objects, function(object) {
    oxygen_ok <- all(sprintf("%.12f", c(0, 5)) %in%
      sprintf("%.12f", object$o2_values))
    time_ok <- all(c(0, 1000) %in% object$day_values)
    if (mode == "oxygen") return(oxygen_ok)
    if (mode == "time") return(time_ok)
    all(c(2, 4) %in% object$initial_ploidy) &&
      all(sprintf("%.12f", f6ft_p_values()) %in%
            sprintf("%.12f", object$p_misseg)) &&
      time_ok && oxygen_ok
  }, logical(1L)))
}

f6ab_build <- function(paths, ploidy_objects, growth_objects, family_b = "C01",
                       growth_limits, panel_mm = 38, show_a = TRUE) {
  stopifnot(family_b %in% c("C01", "C02"),
            identical(names(ploidy_objects), c("in vivo", "in vitro")),
            identical(names(growth_objects), c("in vivo", "in vitro")),
            f6ab_panel_coverage(c(ploidy_objects, growth_objects)))
  base <- if (show_a) f6ab_surface_base(paths, c(0, 5)) else NULL
  context <- c("in vivo", "in vivo", "in vitro", "in vitro")
  family <- c("C01", "C02", "C01", "C02")
  initial <- c(4, 2, 4, 2)
  context_color <- c("in vivo" = "#0072B2", "in vitro" = "#CC79A7")
  family_color <- c(C01 = "#C99700", C02 = "#6A3D9A")
  gap <- .55
  cluster_gap <- 28
  row_gap <- 2.0
  row_top <- 24 + c(0, panel_mm + row_gap, 2 * panel_mm + row_gap + 5,
                    3 * panel_mm + 2 * row_gap + 5)
  a_x <- 23
  b_x <- if (show_a) a_x + panel_mm + 27 else 23
  b_columns <- b_x + (0:9) * (panel_mm + gap) +
    c(rep(0, 5), rep(cluster_gap - gap, 5))
  c_x <- b_columns[[6L]]
  group_span <- 5 * panel_mm + 4 * gap
  width <- tail(b_columns, 1) + panel_mm + 5
  bottom <- tail(row_top, 1) + panel_mm
  height <- bottom + 29
  children <- list()
  geometry <- list()
  add <- function(grob, x, y, w, h, clip = "off") {
    children[[length(children) + 1L]] <<- grid::grobTree(grob,
      vp = grid::viewport(x = grid::unit(x, "mm"), y = grid::unit(height - y, "mm"),
        width = grid::unit(w, "mm"), height = grid::unit(h, "mm"),
        just = c("left", "top"), clip = clip))
  }
  text <- function(label, x, y, w, h, size = 7, bold = FALSE, rotate = 0,
                   just = "centre", colour = "#222222") {
    if (length(label) != 1L) stop("Figure 6 text must be drawn as a whole label.")
    add(grid::textGrob(label, x = if (just == "left") 0 else .5,
      just = just, rot = rotate,
      gp = grid::gpar(fontfamily = "Helvetica", fontsize = size, col = colour,
                     fontface = if (bold) "bold" else "plain")), x, y, w, h)
  }
  strip <- function(label, color, x, y, w, h, rotate = 0,
                    text_colour = "#222222") {
    add(grid::rectGrob(gp = grid::gpar(fill = color, col = "#BEBEBE", lwd = .4)), x, y, w, h)
    text(label, x, y, w, h, size = 8.2, bold = TRUE, rotate = rotate,
         colour = text_colour)
  }
  tile <- function(parts, x, y, label, show_y = FALSE, show_x = FALSE) {
    add(parts$panel, x, y, panel_mm, panel_mm, "on")
    if (show_y) add(parts$y, x - 6.5, y, 6.5, panel_mm)
    if (show_x) add(parts$x, x, y + panel_mm, panel_mm, 7)
    geometry[[length(geometry) + 1L]] <<- data.frame(
      panel = label, x_mm = x, y_mm = y, width_mm = panel_mm, height_mm = panel_mm)
  }
  if (show_a) text("A. Steady-state ploidy", 5, 1, a_x + panel_mm - 5, 6,
                   9.5, TRUE, just = "left")
  text(if (show_a) "B. Finite-time mean ploidy" else
         "A. Finite-time mean ploidy", b_x - 18, 1,
       group_span + 18, 6, 9.5, TRUE, just = "left")
  text(if (show_a) "C. Finite-time net growth rate" else
         "B. Finite-time net growth rate", c_x - 18, 1,
       group_span + 18, 6, 9.5, TRUE, just = "left")
  text("Fixed p_misseg", b_x, 7, group_span, 5, 9, TRUE)
  text("Fixed p_misseg", c_x, 7, group_span, 5, 9, TRUE)
  for (group in 1:2) {
    start <- if (group == 1) 1L else 3L
    context_height <- 2 * panel_mm + row_gap
    if (show_a) {
      strip(context[start], context_color[[context[start]]], 5, row_top[start],
            3.8, context_height, 90, text_colour = "#FFFFFF")
      text("Fixed oxygen (%)", 12.8, row_top[start], 3.7,
           context_height, 8.8, rotate = 90)
    }
    strip(context[start], context_color[[context[start]]], b_x - 18.5, row_top[start], 3.8, context_height, 90,
          text_colour = "#FFFFFF")
    text("Fixed oxygen (%)", b_x - 14.7, row_top[start], 5, context_height, 8.8, rotate = 90)
    strip(context[start], context_color[[context[start]]], c_x - 18.5,
          row_top[start], 3.8, context_height, 90,
          text_colour = "#FFFFFF")
    text("Fixed oxygen (%)", c_x - 14.7, row_top[start], 5,
         context_height, 8.8, rotate = 90)
  }
  for (column_group in 1:2) {
    index <- if (column_group == 1) 1:5 else 6:10
    strip(family_b,
      family_color[[family_b]],
      b_columns[index[1]], 14, 5 * panel_mm + 4 * gap, 4,
      text_colour = "#FFFFFF")
    for (j in seq_along(index)) strip(f6ft_format_p(f6ft_p_values()[j]), "#F2F2F2",
      b_columns[index[j]], 18.3, panel_mm, 4.8)
  }
  for (row in 1:4) {
    if (show_a) {
      a <- base
      a$layers <- lapply(base$layers, function(layer) {
        copy <- ggplot2::ggproto(NULL, layer)
        d <- layer$data
        copy$data <- d[as.character(d$model_context) == context[row] &
                        as.character(d$display_label) == family[row], , drop = FALSE]
        copy
      })
      a <- a + ggplot2::facet_null() + ggplot2::labs(
        title = NULL, subtitle = NULL, x = NULL, y = NULL) + f6ab_tile_theme(45)
      tile(f6ab_parts(a), a_x, row_top[row], paste0("A", row), TRUE, row == 4)
      strip(family[row], family_color[[family[row]]], 8.8, row_top[row],
            4, panel_mm, 90, text_colour = "#FFFFFF")
    }
    strip(paste0(initial[row], "N"), "#F2F2F2", b_x - 9.7, row_top[row], 4.2, panel_mm, 90)
    strip(paste0(initial[row], "N"), "#F2F2F2",
          c_x - 9.7, row_top[row], 4.2, panel_mm, 90)
    for (metric_group in 1:2) {
      metric <- c("mean_ploidy", "mean_net_growth_rate")[[metric_group]]
      object <- if (metric_group == 1L) ploidy_objects[[context[[row]]]] else
        growth_objects[[context[[row]]]]
      for (p in 1:5) {
        z <- f6sb_extract(object, metric, initial[row], f6ft_p_values()[p],
                          c(0, 1000), c(0, 5))
        fill_scale <- if (metric_group == 1L) f6ab_fill() else
          f6sb_growth_fill(growth_limits, growth_scale = "signed_log")
        b <- ggplot2::ggplot(z, ggplot2::aes(day, O2_pct, fill = value)) +
          ggplot2::geom_raster(interpolate = FALSE) + fill_scale +
          ggplot2::scale_x_continuous(
            breaks = c(0, 1000),
            labels = if (p == 1) c("0", "") else if (p == 5) c("", "1000") else c("", ""),
            expand = c(0, 0)
          ) +
          ggplot2::scale_y_continuous(breaks = 0:5, expand = c(0, 0)) +
          ggplot2::coord_cartesian(xlim = c(0, 1000), ylim = c(0, 5), expand = FALSE) +
          ggplot2::labs(x = NULL, y = NULL) + f6ab_tile_theme()
        if (row <= 2) b <- b + ggplot2::geom_hline(yintercept = .5,
          colour = "#7A7A7A", linetype = "dashed", linewidth = .28)
        col <- (metric_group - 1L) * 5L + p
        panel_id <- if (metric_group == 1L) {
          if (show_a) "B" else "A"
        } else {
          if (show_a) "C" else "B"
        }
        tile(f6ab_parts(b), b_columns[col], row_top[row],
             paste0(panel_id, row, "_", p),
             col == 1L || col == 6L, row == 4L)
      }
    }
  }
  if (show_a) text("p_misseg", a_x, bottom + 7.5,
                   panel_mm, 5, 9)
  text("Experimental time (days)", b_x, bottom + 7.5,
       group_span, 5, 9)
  text("Experimental time (days)", c_x, bottom + 7.5,
       group_span, 5, 9)
  add_legend <- function(values, fill_scale, x) {
    legend_plot <- ggplot2::ggplot(data.frame(value = values, y = 1),
      ggplot2::aes(value, y, fill = value)) + ggplot2::geom_tile() +
      fill_scale + ggplot2::theme_void(base_family = "Helvetica") +
      ggplot2::theme(legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 8),
        legend.text = ggplot2::element_text(size = 7)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(
        title.position = "top", title.hjust = .5,
        barwidth = grid::unit(58, "mm"), barheight = grid::unit(3, "mm")))
    legend <- ggplot2::ggplotGrob(legend_plot)
    index <- which(legend$layout$name == "guide-box-bottom")
    if (!length(index)) index <- which(legend$layout$name == "guide-box")
    add(legend$grobs[[index[[1L]]]], x, bottom + 13.5, 84, 12)
  }
  add_legend(1:7, f6ab_fill(), b_columns[[1L]] + (group_span - 84) / 2)
  add_legend(seq(growth_limits[[1L]], growth_limits[[2L]], length.out = 101L),
    f6sb_growth_fill(growth_limits, growth_scale = "signed_log"),
    b_columns[[6L]] + (group_span - 84) / 2)
  list(plot = do.call(grid::grobTree, children), width = width / 25.4, height = height / 25.4,
       geometry = do.call(rbind, geometry))
}

f6ft_draw_main <- function(workspace_root = f6r_find_workspace_root()) {
  paths <- f6r_paths(workspace_root)
  run_paths <- f6g_paths(paths)
  growth_run <- f6ng_paths(paths)
  ploidy <- list(
    "in vivo" = f6g_read_panel(run_paths, "in vivo", "continuous", "C01"),
    "in vitro" = f6g_read_panel(run_paths, "in vitro", "passage", "C01"))
  growth <- list(
    "in vivo" = f6ng_read_panel(growth_run, "in vivo", "C01"),
    "in vitro" = f6ng_read_panel(growth_run, "in vitro", "C01"))
  growth_limits <- f6ab_growth_limits(paths)
  layout <- f6ab_build(paths, ploidy, growth, growth_limits = growth_limits)
  geometry <- layout$geometry
  stopifnot(nrow(geometry) == 44L, all(geometry$width_mm == geometry$height_mm),
    length(unique(geometry$width_mm)) == 1L, length(unique(geometry$y_mm)) == 4L,
    sum(startsWith(geometry$panel, "A")) == 4L,
    sum(startsWith(geometry$panel, "B")) == 20L,
    sum(startsWith(geometry$panel, "C")) == 20L)
  dir.create(run_paths$rendered, recursive = TRUE, showWarnings = FALSE)
  output <- f6r_save_plot(layout$plot, file.path(run_paths$rendered, "assembled_fig6"),
    width = layout$width, height = layout$height, dpi = 300)
  published <- f6g_publish_plot(output, paths, "assembled_fig6")
  f6ft_atomic_write_tsv(geometry, file.path(run_paths$run_root, "figure6_ab_panel_geometry.tsv"))
  validation <- f6g_render_hash_validation(output, published,
    file.path(run_paths$run_root, "figure6_ab_render_validation.tsv"),
    data.frame(check = c("same_square_panel_size", "four_aligned_rows",
                         "four_A_twenty_B_twenty_C", "B_C_gap_at_least_25mm",
                         "shared_ploidy_color_limits_1_7", "signed_log_growth_colors",
                         "oxygen_range_0_5", "time_range_0_1000"),
      observed = c(length(unique(geometry$width_mm)) == 1L && all(geometry$width_mm == geometry$height_mm),
        length(unique(geometry$y_mm)) == 4L,
        sum(startsWith(geometry$panel, "A")) == 4L &&
          sum(startsWith(geometry$panel, "B")) == 20L &&
          sum(startsWith(geometry$panel, "C")) == 20L,
        min(geometry$x_mm[startsWith(geometry$panel, "C")]) -
          max(geometry$x_mm[startsWith(geometry$panel, "B")] +
                geometry$width_mm[startsWith(geometry$panel, "B")]) >= 25,
        isTRUE(all.equal(f6ab_fill()$get_transformation()$inverse(f6ab_fill()$limits),
                         c(1, 7))),
        growth_limits[[1L]] < 0 && growth_limits[[2L]] > 0,
        f6ab_panel_coverage(c(ploidy, growth), "oxygen"),
        f6ab_panel_coverage(c(ploidy, growth), "time")), expected = TRUE))
  invisible(list(output = output, published = published, validation = validation))
}

f6ft_draw_supp6_17 <- function(workspace_root = f6r_find_workspace_root()) {
  paths <- f6r_paths(workspace_root)
  run_paths <- f6g_paths(paths)
  growth_run <- f6ng_paths(paths)
  ploidy <- list(
    "in vivo" = f6g_read_panel(run_paths, "in vivo", "continuous", "C02"),
    "in vitro" = f6g_read_panel(run_paths, "in vitro", "passage", "C02"))
  growth <- list(
    "in vivo" = f6ng_read_panel(growth_run, "in vivo", "C02"),
    "in vitro" = f6ng_read_panel(growth_run, "in vitro", "C02"))
  growth_limits <- f6ab_growth_limits(paths)
  layout <- f6ab_build(paths, ploidy, growth, family_b = "C02",
    growth_limits = growth_limits, show_a = FALSE)
  geometry <- layout$geometry
  stopifnot(nrow(geometry) == 40L,
    all(geometry$width_mm == 38 & geometry$height_mm == 38),
    length(unique(geometry$x_mm)) == 10L,
    length(unique(geometry$y_mm)) == 4L,
    sum(startsWith(geometry$panel, "A")) == 20L,
    sum(startsWith(geometry$panel, "B")) == 20L)
  directory <- file.path(paths$root, "data", "Figures", "Supp_Figure6_17")
  rendered <- file.path(directory, "rendered")
  dir.create(rendered, recursive = TRUE, showWarnings = FALSE)
  name <- "supp_fig6-17_c02_finite_time_ploidy_and_net_growth_signed_log"
  output <- f6r_save_plot(layout$plot, file.path(rendered, name),
    width = layout$width, height = layout$height, dpi = 300)
  published <- f6g_publish_plot(output, paths, name)
  f6ft_atomic_write_tsv(geometry,
    file.path(directory, paste0(name, "_panel_geometry.tsv")))
  validation <- f6g_render_hash_validation(output, published,
    file.path(directory, paste0(name, "_render_validation.tsv")),
    data.frame(check = c("forty_panels", "square_38mm", "ten_columns",
                         "four_rows", "twenty_A_twenty_B", "A_B_gap_at_least_25mm",
                         "oxygen_range_0_5", "time_range_0_1000",
                         "signed_log_growth_colors"),
      observed = c(nrow(geometry) == 40L,
                   all(geometry$width_mm == 38 & geometry$height_mm == 38),
                   length(unique(geometry$x_mm)) == 10L,
                   length(unique(geometry$y_mm)) == 4L,
                   sum(startsWith(geometry$panel, "A")) == 20L &&
                     sum(startsWith(geometry$panel, "B")) == 20L,
                   min(geometry$x_mm[startsWith(geometry$panel, "B")]) -
                     max(geometry$x_mm[startsWith(geometry$panel, "A")] +
                           geometry$width_mm[startsWith(geometry$panel, "A")]) >= 25,
                   f6ab_panel_coverage(c(ploidy, growth), "oxygen"),
                   f6ab_panel_coverage(c(ploidy, growth), "time"),
                   growth_limits[[1L]] < 0 && growth_limits[[2L]] > 0),
      expected = rep(TRUE, 9L)))
  invisible(list(output = output, published = published, validation = validation))
}

f6ft_draw_supp6_8 <- function(workspace_root = f6r_find_workspace_root()) {
  paths <- f6r_paths(workspace_root)
  plot <- f6ab_surface_base(paths, c(0, 5)) + ggplot2::labs(
    title = "Steady-state ploidy response: full oxygen range", subtitle = NULL,
    caption = NULL) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
      panel.spacing.x = grid::unit(c(7, 10, 7), "mm"))
  directory <- file.path(paths$root, "data", "Figures", "Supp_Figure6_8")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  output <- f6r_save_plot(plot, file.path(directory, "supp_fig6-8_steady_state_full_oxygen_range"),
    width = 18.2, height = 6.2)
  published <- f6g_publish_plot(output, paths, "supp_fig6-8_steady_state_full_oxygen_range")
  f6g_render_hash_validation(output, published,
    file.path(directory, "supp_fig6-8_steady_state_full_oxygen_range_render_validation.tsv"))
}
