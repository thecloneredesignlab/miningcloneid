#!/usr/bin/env Rscript

# Rendering utilities for Figure 7 finite-time q10 panels and diagnostics.

options(stringsAsFactors = FALSE, warn = 1)

f7ft_read_panel_object <- function(run_paths, panel_letter) {
  path <- file.path(
    run_paths$run_root,
    paste0("finite_time_panel_", tolower(panel_letter), ".rds")
  )
  f7r_require_files(path, paste0("Figure 7", panel_letter, " finite-time data"))
  object <- readRDS(path)
  if (!identical(object$profile, f7ft_profile()) ||
      !identical(object$panel_letter, panel_letter)) {
    stop("Unexpected finite-time panel object: ", path)
  }
  object
}

f7ft_panel_long <- function(object) {
  values <- object$mean_ploidy
  expected_dim <- c(
    length(object$initial_ploidy), length(object$day_values),
    length(object$o2_values), length(object$p_misseg)
  )
  if (!identical(as.integer(dim(values)), as.integer(expected_dim))) {
    stop("Finite-time panel array dimensions do not match its coordinate grids.")
  }
  n_initial <- length(object$initial_ploidy)
  n_day <- length(object$day_values)
  n_o2 <- length(object$o2_values)
  n_p <- length(object$p_misseg)
  data.frame(
    initial_ploidy = rep(
      object$initial_ploidy, times = n_day * n_o2 * n_p
    ),
    day = rep(
      rep(object$day_values, each = n_initial), times = n_o2 * n_p
    ),
    O2_pct = rep(
      rep(object$o2_values, each = n_initial * n_day), times = n_p
    ),
    p_misseg = rep(
      object$p_misseg, each = n_initial * n_day * n_o2
    ),
    mean_ploidy = as.vector(values),
    stringsAsFactors = FALSE
  )
}

f7ft_format_p <- function(x) {
  ifelse(abs(x - 0.005) < 1e-12, "0.005", sprintf("%.2f", x))
}

f7ft_panel_title <- function(object) {
  if (!is.null(object$title_override) && nzchar(object$title_override)) {
    return(object$title_override)
  }
  paste0(
    object$panel_letter, ". ", object$pair_label, " - ", object$model_context
  )
}

f7ft_initial_row_levels <- function(initial_ploidy = f7ft_initial_ploidy()) {
  paste0(rev(as.numeric(initial_ploidy)), "N")
}

f7ft_finite_time_panel_plot <- function(
    object, show_legend = TRUE,
    initial_ploidy_display = object$initial_ploidy,
    o2_limits = range(object$o2_values),
    o2_breaks = NULL,
    reference_o2 = NULL
) {
  f7r_require_packages(c("ggplot2", "scales"))
  data <- f7ft_panel_long(object)
  initial_ploidy_display <- as.numeric(initial_ploidy_display)
  o2_limits <- as.numeric(o2_limits)
  if (length(o2_limits) != 2L || any(!is.finite(o2_limits)) ||
      o2_limits[[1L]] >= o2_limits[[2L]]) {
    stop("o2_limits must contain two increasing finite values.")
  }
  if (!all(initial_ploidy_display %in% object$initial_ploidy)) {
    stop("Requested initial-ploidy rows are absent from the panel object.")
  }
  data <- data[
    data$initial_ploidy %in% initial_ploidy_display &
      data$O2_pct >= o2_limits[[1L]] - 1e-12 &
      data$O2_pct <= o2_limits[[2L]] + 1e-12,
    , drop = FALSE
  ]
  data$initial_label <- factor(
    paste0(data$initial_ploidy, "N"),
    levels = f7ft_initial_row_levels(initial_ploidy_display)
  )
  data$p_label <- factor(
    f7ft_format_p(data$p_misseg),
    levels = f7ft_format_p(f7ft_p_values())
  )
  is_passage <- isTRUE(grepl(
    "^passage_constrained_expm", as.character(object$propagation_mode)
  ))
  x_breaks <- if (!is.null(object$x_breaks)) {
    as.numeric(object$x_breaks)
  } else if (is_passage) {
    unique(pretty(range(object$day_values), n = 4))
  } else {
    c(0, 500, 1000)
  }
  if (is.null(o2_breaks)) {
    o2_breaks <- if (isTRUE(all.equal(o2_limits, c(0, 2)))) {
      c(0, 0.5, 1, 1.5, 2)
    } else if (isTRUE(all.equal(o2_limits, c(0, 20)))) {
      c(0, 5, 10, 15, 20)
    } else {
      unique(pretty(o2_limits, n = 3))
    }
  }
  subtitle <- if (!is.null(object$subtitle_override)) {
    object$subtitle_override
  } else if (is_passage) {
    paste0(
      "Rows (bottom to top): initial ploidy ",
      min(initial_ploidy_display), "N to ", max(initial_ploidy_display), "N; ",
      "columns: fixed p_misseg parameter; ",
      "fitted passage timing, target-cell state selection, and reseeding"
    )
  } else {
    paste0(
      "Rows (bottom to top): initial ploidy ",
      min(initial_ploidy_display), "N to ", max(initial_ploidy_display), "N; ",
      "columns: fixed p_misseg parameter; ",
      "daily expm propagation"
    )
  }
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = day, y = O2_pct, fill = mean_ploidy)
  ) +
    ggplot2::geom_raster(interpolate = FALSE)
  if (!is.null(reference_o2)) {
    reference_o2 <- as.numeric(reference_o2)
    if (length(reference_o2) != 1L || !is.finite(reference_o2) ||
        reference_o2 < o2_limits[[1L]] || reference_o2 > o2_limits[[2L]]) {
      stop("reference_o2 must be one finite value inside o2_limits.")
    }
    plot <- plot + ggplot2::geom_hline(
      yintercept = reference_o2,
      colour = "#7A7A7A", linetype = "dashed", linewidth = 0.28
    )
  }
  plot +
    ggplot2::facet_grid(initial_label ~ p_label, switch = "y") +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "#FFFFBF", "#B2182B"),
      trans = "log10", limits = c(1, 7),
      breaks = c(1, 1.5, 2, 3, 4, 6),
      name = paste0(
        "Mean finite-time ploidy\n",
        "across 50 q10 optimizer endpoints\n(log colors)"
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks, expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = o2_breaks, expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(ylim = o2_limits, expand = FALSE) +
    ggplot2::labs(
      title = f7ft_panel_title(object),
      subtitle = subtitle,
      x = if (is_passage) "Cumulative experimental time (day)" else "Time (day)",
      y = if (is_passage) "Fixed oxygen across passages (%)" else "Fixed oxygen (%)"
    ) +
    ggplot2::theme_classic(base_size = 7.7, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.2),
      plot.subtitle = ggplot2::element_text(size = 6.7, colour = "#555555"),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", colour = "#BEBEBE", linewidth = 0.25
      ),
      strip.text.x = ggplot2::element_text(face = "bold", size = 7.2),
      strip.text.y = ggplot2::element_text(face = "bold", size = 7.2),
      strip.placement = "outside",
      panel.spacing = grid::unit(0.55, "mm"),
      panel.border = ggplot2::element_rect(
        colour = "#555555", fill = NA, linewidth = 0.22
      ),
      aspect.ratio = 1,
      axis.text = ggplot2::element_text(size = 5.7, colour = "#333333"),
      axis.title = ggplot2::element_text(size = 7.5),
      legend.position = if (isTRUE(show_legend)) "bottom" else "none",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 6.8),
      legend.text = ggplot2::element_text(size = 6.2),
      legend.key.width = grid::unit(13, "mm"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        direction = "horizontal", title.position = "top",
        title.hjust = 0.5, barwidth = grid::unit(46, "mm"),
        barheight = grid::unit(2.6, "mm")
      )
    )
}

f7ft_figure7_validation <- function(
    paths, run_paths, output, published, panel_objects, image_info,
    display_specs
) {
  panels <- names(panel_objects)
  checks <- data.frame(
    check = c(
      "main_panel_count", "finite_time_panel_count",
      "finite_time_sources_have_5_initial_ploidies",
      "finite_time_display_has_3_initial_ploidies",
      "finite_time_display_rows_top_to_bottom_4N_to_2N",
      "finite_time_each_has_5_p_values",
      "finite_time_each_has_201_o2_values",
      "finite_time_each_has_21_o2_values_from_0_to_0p5",
      "finite_time_display_o2_is_0_to_2",
      "invivo_reference_line_is_o2_0p5",
      "invitro_has_no_reference_line",
      "invivo_has_1001_days",
      "invitro_uses_passage_experimental_days",
      "invitro_has_six_lineage_schedules",
      "finite_time_weights_equal_50",
      "png_width_px", "png_height_px",
      "figures_png_md5_match", "figures_pdf_md5_match",
      "manuscript_png_md5_match", "manuscript_pdf_md5_match"
    ),
    observed = c(
      6L, length(panels),
      all(vapply(panel_objects, function(x) length(x$initial_ploidy) == 5L, logical(1L))),
      all(vapply(display_specs, function(x) identical(x$initial_ploidy_display, 2:4), logical(1L))),
      identical(f7ft_initial_row_levels(2:4), paste0(4:2, "N")),
      all(vapply(panel_objects, function(x) length(x$p_misseg) == 5L, logical(1L))),
      all(vapply(panel_objects, function(x) length(x$o2_values) == 201L, logical(1L))),
      all(vapply(panel_objects, function(x) {
        sum(x$o2_values >= 0 & x$o2_values <= 0.5 + 1e-12) == 21L
      }, logical(1L))),
      all(vapply(display_specs, function(x) isTRUE(all.equal(x$o2_limits, c(0, 2))), logical(1L))),
      all(vapply(display_specs[c("C", "D")], function(x) identical(x$reference_o2, 0.5), logical(1L))),
      all(vapply(display_specs[c("E", "F")], function(x) is.null(x$reference_o2), logical(1L))),
      all(vapply(panel_objects[c("C", "D")], function(x) length(x$day_values) == 1001L, logical(1L))),
      all(vapply(panel_objects[c("E", "F")], function(x) {
        min(x$day_values) == 0 && max(x$day_values) < 1000 &&
          identical(as.integer(x$day_values), 0:max(x$day_values))
      }, logical(1L))),
      all(vapply(panel_objects[c("E", "F")], function(x) identical(x$n_lineage_schedule, 6L), logical(1L))),
      all(vapply(panel_objects, function(x) all(x$optimizer_endpoint_weight == 50L), logical(1L))),
      image_info$width[[1L]], image_info$height[[1L]],
      f7r_md5(output[["png"]]) == f7r_md5(published[["figures_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["figures_pdf"]]),
      f7r_md5(output[["png"]]) == f7r_md5(published[["manuscript_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(
      6L, 4L, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE, TRUE, TRUE, 5040L, 6480L, TRUE, TRUE, TRUE, TRUE
    ),
    stringsAsFactors = FALSE
  )
  checks$passed <- as.character(checks$observed) == as.character(checks$expected)
  path <- f7ft_atomic_write_tsv(
    checks, file.path(run_paths$run_root, "figure7_finite_time_render_validation.tsv")
  )
  if (!all(checks$passed)) stop("Rendered Figure 7 finite-time validation failed.")
  path
}

f7ft_draw_main <- function(workspace_root = f7r_find_workspace_root()) {
  f7r_require_packages(c(
    "ggplot2", "patchwork", "magick", "isoband", "scales"
  ))
  paths <- f7r_paths(workspace_root)
  base_run_paths <- f7ft_paths(paths, run_id = NULL, create = FALSE)
  run_paths <- f7p_paths(paths, run_id = NULL, create = TRUE)
  panel_objects <- list(
    C = f7ft_read_panel_object(base_run_paths, "C"),
    D = f7ft_read_panel_object(base_run_paths, "D"),
    E = f7p_read_panel_object(run_paths, "E"),
    F = f7p_read_panel_object(run_paths, "F")
  )
  p_a <- f7x_main_surface_plot(paths) +
    ggplot2::theme(legend.position = "bottom")
  p_b <- f7x_main_inverse_plot(paths) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center"
    ) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(
        nrow = 1L, byrow = TRUE, order = 1L, title.position = "left"
      ),
      fill = ggplot2::guide_colourbar(
        barwidth = grid::unit(52, "mm"),
        barheight = grid::unit(3.2, "mm"),
        order = 2L, title.position = "left"
      )
    )
  display_specs <- list(
    C = list(initial_ploidy_display = 2:4, o2_limits = c(0, 2), reference_o2 = 0.5),
    D = list(initial_ploidy_display = 2:4, o2_limits = c(0, 2), reference_o2 = 0.5),
    E = list(initial_ploidy_display = 2:4, o2_limits = c(0, 2), reference_o2 = NULL),
    F = list(initial_ploidy_display = 2:4, o2_limits = c(0, 2), reference_o2 = NULL)
  )
  finite_plots <- lapply(names(panel_objects), function(letter) {
    do.call(
      f7ft_finite_time_panel_plot,
      c(list(object = panel_objects[[letter]]), display_specs[[letter]])
    )
  })
  names(finite_plots) <- names(panel_objects)

  panel_paths <- c()
  for (letter in names(finite_plots)) {
    panel_paths <- c(
      panel_paths,
      stats::setNames(
        f7r_save_plot(
          finite_plots[[letter]],
          file.path(
            run_paths$panels,
            paste0("figure7", tolower(letter), "_finite_time_q10")
          ),
          width = 8.2, height = 6.5
        ),
        paste0(letter, c("_png", "_pdf"))
      )
    )
  }

  top <- p_a + p_b + patchwork::plot_layout(widths = c(1, 1))
  finite_block <- (
    (finite_plots$C + finite_plots$D) /
      (finite_plots$E + finite_plots$F) +
      patchwork::plot_layout(guides = "collect", heights = c(1, 1))
  ) & ggplot2::theme(legend.position = "bottom")
  combined <- (top / finite_block +
    patchwork::plot_layout(heights = c(1, 2.05)) +
    patchwork::plot_annotation(
      caption = paste0(
        "C-D show continuous fixed-environment expm solutions; E-F add the ",
        "in-vitro fitting passage clock, target-cell state selection, and ",
        "reseeding. All panels average the same q10 50 optimizer endpoints ",
        "used for A-B, retaining repeated endpoint multiplicity."
      )
    )) &
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        size = 7, colour = "#555555", hjust = 0
      )
    )

  output <- f7r_save_plot(
    combined, file.path(run_paths$rendered, "assembled_fig7"),
    width = 16.8, height = 21.6, dpi = 300
  )
  published <- c(
    figures_png = f7r_publish(
      output[["png"]], file.path(paths$figures, "assembled_fig7.png")
    ),
    figures_pdf = f7r_publish(
      output[["pdf"]], file.path(paths$figures, "assembled_fig7.pdf")
    ),
    manuscript_png = f7r_publish(
      output[["png"]], file.path(paths$manuscript_figures, "assembled_fig7.png")
    ),
    manuscript_pdf = f7r_publish(
      output[["pdf"]], file.path(paths$manuscript_figures, "assembled_fig7.pdf")
    )
  )
  info <- magick::image_info(magick::image_read(output[["png"]]))
  validation <- f7ft_figure7_validation(
    paths, run_paths, output, published, panel_objects, info, display_specs
  )
  invisible(list(
    plot = combined, output = output, published = published,
    panel_paths = panel_paths, validation = validation
  ))
}

f7ft_draw_supp7_8 <- function(workspace_root = f7r_find_workspace_root()) {
  f7r_require_packages(c("ggplot2", "patchwork", "magick", "scales"))
  paths <- f7r_paths(workspace_root)
  base_run_paths <- f7ft_paths(paths, run_id = NULL, create = FALSE)
  output_dir <- file.path(paths$root, "data", "Figures", "Supp_Figure7_8")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  panel_a <- f7ft_read_panel_object(base_run_paths, "E")
  panel_b <- f7ft_read_panel_object(base_run_paths, "F")
  panel_a$panel_letter <- "A"
  panel_b$panel_letter <- "B"
  panel_a$title_override <- "A. C01 - in vitro, no passaging"
  panel_b$title_override <- "B. C02 - in vitro, no passaging"
  common_subtitle <- paste0(
    "Rows (bottom to top): initial ploidy 2N to 6N; columns: fixed p_misseg ",
    "parameter; continuous fixed-environment expm propagation"
  )
  panel_a$subtitle_override <- common_subtitle
  panel_b$subtitle_override <- common_subtitle
  panel_a$propagation_mode <- "continuous_no_passaging"
  panel_b$propagation_mode <- "continuous_no_passaging"

  combined <- f7ft_finite_time_panel_plot(panel_a) +
    f7ft_finite_time_panel_plot(panel_b) +
    patchwork::plot_layout(guides = "collect", widths = c(1, 1)) +
    patchwork::plot_annotation(
      caption = paste0(
        "Archived continuous-culture counterfactuals formerly shown as Figure 7E-F. ",
        "Oxygen and p_misseg remain fixed for 1,000 days and no ",
        "passage boundary or reseeding operation is applied."
      )
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(size = 7, colour = "#555555", hjust = 0)
    )
  output <- f7r_save_plot(
    combined,
    file.path(output_dir, "supp_fig7-8_continuous_invitro_no_passaging"),
    width = 16.8, height = 10.2, dpi = 300
  )
  published <- c(
    figures_png = f7r_publish(
      output[["png"]],
      file.path(paths$figures, "supp_fig7-8_continuous_invitro_no_passaging.png")
    ),
    figures_pdf = f7r_publish(
      output[["pdf"]],
      file.path(paths$figures, "supp_fig7-8_continuous_invitro_no_passaging.pdf")
    ),
    manuscript_png = f7r_publish(
      output[["png"]],
      file.path(paths$manuscript_figures, "supp_fig7-8_continuous_invitro_no_passaging.png")
    ),
    manuscript_pdf = f7r_publish(
      output[["pdf"]],
      file.path(paths$manuscript_figures, "supp_fig7-8_continuous_invitro_no_passaging.pdf")
    )
  )
  info <- magick::image_info(magick::image_read(output[["png"]]))
  checks <- data.frame(
    check = c(
      "panel_count", "source_profile", "source_days", "source_o2_values",
      "source_p_values", "source_optimizer_weight", "png_width_px",
      "png_height_px", "figures_png_md5_match", "figures_pdf_md5_match",
      "manuscript_png_md5_match", "manuscript_pdf_md5_match"
    ),
    observed = c(
      2L,
      identical(panel_a$profile, f7ft_profile()) && identical(panel_b$profile, f7ft_profile()),
      length(panel_a$day_values) == 1001L && length(panel_b$day_values) == 1001L,
      length(panel_a$o2_values) == 201L && length(panel_b$o2_values) == 201L,
      length(panel_a$p_misseg) == 5L && length(panel_b$p_misseg) == 5L,
      all(panel_a$optimizer_endpoint_weight == 50L) && all(panel_b$optimizer_endpoint_weight == 50L),
      info$width[[1L]], info$height[[1L]],
      f7r_md5(output[["png"]]) == f7r_md5(published[["figures_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["figures_pdf"]]),
      f7r_md5(output[["png"]]) == f7r_md5(published[["manuscript_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(
      2L, TRUE, TRUE, TRUE, TRUE, TRUE, 5040L, 3060L,
      TRUE, TRUE, TRUE, TRUE
    ),
    stringsAsFactors = FALSE
  )
  checks$passed <- as.character(checks$observed) == as.character(checks$expected)
  validation <- f7ft_atomic_write_tsv(
    checks, file.path(output_dir, "supp_fig7-8_render_validation.tsv")
  )
  if (!all(checks$passed)) stop("Supplementary Figure 7-8 render validation failed.")
  invisible(list(plot = combined, output = output, published = published, validation = validation))
}

f7ft_draw_full_pair_supp <- function(
    paths, panel_a, panel_b, output_dir_name, filename, caption,
    title_a, title_b, subtitle, expected_profile,
    expected_days, expected_o2, expected_mode
) {
  f7r_require_packages(c("ggplot2", "patchwork", "magick", "scales"))
  output_dir <- file.path(paths$root, "data", "Figures", output_dir_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  panel_a$panel_letter <- "A"
  panel_b$panel_letter <- "B"
  panel_a$title_override <- title_a
  panel_b$title_override <- title_b
  panel_a$subtitle_override <- subtitle
  panel_b$subtitle_override <- subtitle
  panel_a$x_breaks <- if (max(expected_days) == 1000) c(0, 500, 1000) else NULL
  panel_b$x_breaks <- panel_a$x_breaks
  o2_breaks <- if (max(expected_o2) == 20) c(0, 5, 10, 15, 20) else c(0, 2.5, 5)

  plot_a <- f7ft_finite_time_panel_plot(
    panel_a, initial_ploidy_display = 2:6,
    o2_limits = range(expected_o2), o2_breaks = o2_breaks
  )
  plot_b <- f7ft_finite_time_panel_plot(
    panel_b, initial_ploidy_display = 2:6,
    o2_limits = range(expected_o2), o2_breaks = o2_breaks
  )
  combined <- (plot_a + plot_b +
    patchwork::plot_layout(guides = "collect", widths = c(1, 1)) +
    patchwork::plot_annotation(caption = caption)) &
    ggplot2::theme(
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(
        size = 7, colour = "#555555", hjust = 0
      )
    )
  output <- f7r_save_plot(
    combined, file.path(output_dir, filename),
    width = 16.8, height = 10.2, dpi = 300
  )
  published <- c(
    figures_png = f7r_publish(
      output[["png"]], file.path(paths$figures, paste0(filename, ".png"))
    ),
    figures_pdf = f7r_publish(
      output[["pdf"]], file.path(paths$figures, paste0(filename, ".pdf"))
    ),
    manuscript_png = f7r_publish(
      output[["png"]], file.path(paths$manuscript_figures, paste0(filename, ".png"))
    ),
    manuscript_pdf = f7r_publish(
      output[["pdf"]], file.path(paths$manuscript_figures, paste0(filename, ".pdf"))
    )
  )
  info <- magick::image_info(magick::image_read(output[["png"]]))
  objects <- list(panel_a, panel_b)
  mode_ok <- if (is.null(expected_mode)) {
    all(vapply(objects, function(x) is.null(x$propagation_mode), logical(1L)))
  } else {
    all(vapply(
      objects,
      function(x) identical(x$propagation_mode, expected_mode), logical(1L)
    ))
  }
  checks <- data.frame(
    check = c(
      "panel_count", "source_profile", "propagation_mode",
      "initial_ploidy_grid", "day_grid", "oxygen_grid",
      "p_misseg_parameter_grid", "optimizer_weight", "png_width_px",
      "png_height_px", "figures_png_md5_match", "figures_pdf_md5_match",
      "manuscript_png_md5_match", "manuscript_pdf_md5_match"
    ),
    observed = c(
      length(objects),
      all(vapply(objects, function(x) identical(x$profile, expected_profile), logical(1L))),
      mode_ok,
      all(vapply(objects, function(x) isTRUE(all.equal(as.numeric(x$initial_ploidy), as.numeric(2:6))), logical(1L))),
      all(vapply(objects, function(x) identical(as.numeric(x$day_values), as.numeric(expected_days)), logical(1L))),
      all(vapply(objects, function(x) isTRUE(all.equal(as.numeric(x$o2_values), as.numeric(expected_o2))), logical(1L))),
      all(vapply(objects, function(x) isTRUE(all.equal(as.numeric(x$p_misseg), f7ft_p_values())), logical(1L))),
      all(vapply(objects, function(x) all(x$optimizer_endpoint_weight == 50L), logical(1L))),
      info$width[[1L]], info$height[[1L]],
      f7r_md5(output[["png"]]) == f7r_md5(published[["figures_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["figures_pdf"]]),
      f7r_md5(output[["png"]]) == f7r_md5(published[["manuscript_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(
      2L, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      5040L, 3060L, TRUE, TRUE, TRUE, TRUE
    ), stringsAsFactors = FALSE
  )
  checks$passed <- as.character(checks$observed) == as.character(checks$expected)
  validation <- f7ft_atomic_write_tsv(
    checks, file.path(output_dir, paste0(filename, "_render_validation.tsv"))
  )
  if (!all(checks$passed)) stop("Supplementary render validation failed: ", filename)
  invisible(list(
    plot = combined, output = output, published = published,
    validation = validation
  ))
}

f7ft_draw_supp7_9 <- function(workspace_root = f7r_find_workspace_root()) {
  paths <- f7r_paths(workspace_root)
  run_paths <- f7ft_paths(paths, run_id = NULL, create = FALSE)
  f7ft_draw_full_pair_supp(
    paths,
    f7ft_read_panel_object(run_paths, "C"),
    f7ft_read_panel_object(run_paths, "D"),
    output_dir_name = "Supp_Figure7_9",
    filename = "supp_fig7-9_invivo_finite_time_full_grid",
    caption = paste0(
      "Full-grid in-vivo finite-time counterfactuals archived from Figure 7C-D. ",
      "Oxygen and p_misseg remain fixed during continuous ",
      "daily expm propagation."
    ),
    title_a = "A. C01 - in vivo",
    title_b = "B. C02 - in vivo",
    subtitle = paste0(
      "Rows (bottom to top): initial ploidy 2N to 6N; columns: fixed ",
      "p_misseg parameter; daily expm propagation"
    ),
    expected_profile = f7ft_profile(), expected_days = 0:1000,
    expected_o2 = f7ft_o2_values(), expected_mode = NULL
  )
}

f7ft_draw_supp7_10 <- function(workspace_root = f7r_find_workspace_root()) {
  paths <- f7r_paths(workspace_root)
  run_paths <- f7p_paths(paths, run_id = NULL, create = FALSE)
  panel_e <- f7p_read_panel_object(run_paths, "E")
  panel_f <- f7p_read_panel_object(run_paths, "F")
  f7ft_draw_full_pair_supp(
    paths, panel_e, panel_f,
    output_dir_name = "Supp_Figure7_10",
    filename = "supp_fig7-10_invitro_passage_full_grid",
    caption = paste0(
      "Full-grid in-vitro passage-aware counterfactuals archived from ",
      "Figure 7E-F. Experimental time follows the six fitted lineage schedules ",
      "without extrapolation beyond the recorded passage histories."
    ),
    title_a = "A. C01 - in vitro with passaging",
    title_b = "B. C02 - in vitro with passaging",
    subtitle = paste0(
      "Rows (bottom to top): initial ploidy 2N to 6N; columns: fixed p_misseg ",
      "parameter; fitted passage selection and reseeding"
    ),
    expected_profile = f7p_profile(), expected_days = panel_e$day_values,
    expected_o2 = f7ft_o2_values(), expected_mode = "passage_constrained_expm"
  )
}

f7ft_draw_supp7_11 <- function(workspace_root = f7r_find_workspace_root()) {
  paths <- f7r_paths(workspace_root)
  run_paths <- f7e_paths(paths, mode = "passage", run_id = NULL, create = FALSE)
  f7ft_draw_full_pair_supp(
    paths,
    f7e_read_panel_object(run_paths, "E", mode = "passage"),
    f7e_read_panel_object(run_paths, "F", mode = "passage"),
    output_dir_name = "Supp_Figure7_11",
    filename = "supp_fig7-11_invitro_passage_1000d_extended_o2",
    caption = paste0(
      "Passage-aware in-vitro counterfactuals extended to 1,000 days and ",
      "0-20% oxygen. Each lineage cyclically repeats its complete fitted ",
      "passage schedule after the recorded history ends."
    ),
    title_a = "A. C01 - in vitro with repeated passaging",
    title_b = "B. C02 - in vitro with repeated passaging",
    subtitle = paste0(
      "Rows (bottom to top): initial ploidy 2N to 6N; columns: fixed p_misseg ",
      "parameter; repeated fitted passage schedule"
    ),
    expected_profile = f7e_profile("passage"), expected_days = 0:1000,
    expected_o2 = f7e_o2_values(),
    expected_mode = "passage_constrained_expm_repeated_schedule"
  )
}

f7ft_draw_supp7_12 <- function(workspace_root = f7r_find_workspace_root()) {
  paths <- f7r_paths(workspace_root)
  run_paths <- f7e_paths(paths, mode = "continuous", run_id = NULL, create = FALSE)
  f7ft_draw_full_pair_supp(
    paths,
    f7e_read_panel_object(run_paths, "E", mode = "continuous"),
    f7e_read_panel_object(run_paths, "F", mode = "continuous"),
    output_dir_name = "Supp_Figure7_12",
    filename = "supp_fig7-12_continuous_invitro_1000d_extended_o2",
    caption = paste0(
      "Continuous-culture in-vitro counterfactuals extended to 0-20% oxygen. ",
      "Oxygen and p_misseg remain fixed for 1,000 days, with ",
      "no passage boundary or reseeding operation."
    ),
    title_a = "A. C01 - in vitro, no passaging",
    title_b = "B. C02 - in vitro, no passaging",
    subtitle = paste0(
      "Rows (bottom to top): initial ploidy 2N to 6N; columns: fixed p_misseg ",
      "parameter; continuous fixed-environment expm propagation"
    ),
    expected_profile = f7e_profile("continuous"), expected_days = 0:1000,
    expected_o2 = f7e_o2_values(),
    expected_mode = "continuous_no_passaging_extended_o2"
  )
}

f7ft_read_diagnostics <- function(run_paths) {
  path <- file.path(run_paths$run_root, "finite_time_diagnostic_rows.rds")
  f7r_require_files(path, "Figure 7 finite-time diagnostic rows")
  data <- readRDS(path)
  required <- c(
    "panel_letter", "model_context", "pair_label", "endpoint_multiplicity_q10",
    "day", "eigen_mean_ploidy", "expm_mean_ploidy", "euler_mean_ploidy",
    "steady_mean_ploidy", "spectral_gap", "eigenvector_condition_number"
  )
  if (!all(required %in% names(data))) {
    stop("Finite-time diagnostic table lacks required fields.")
  }
  data$context_label <- factor(
    paste0(data$panel_letter, ". ", data$pair_label, " — ", data$model_context),
    levels = c(
      "C. C01 — in vivo", "D. C02 — in vivo",
      "E. C01 — in vitro", "F. C02 — in vitro"
    )
  )
  data
}

f7ft_weighted_metrics_for_plot <- function(data, x_column, y_column, groups) {
  split_data <- split(data, interaction(data[groups], drop = TRUE, lex.order = TRUE))
  rows <- lapply(split_data, function(z) {
    total_rows <- nrow(z)
    total_weight_all <- sum(
      z$endpoint_multiplicity_q10[is.finite(z$endpoint_multiplicity_q10)],
      na.rm = TRUE
    )
    keep <- is.finite(z[[x_column]]) & is.finite(z[[y_column]]) &
      is.finite(z$endpoint_multiplicity_q10) & z$endpoint_multiplicity_q10 > 0
    z <- z[keep, , drop = FALSE]
    if (!nrow(z)) stop("No finite matched diagnostic rows for plotting metrics.")
    residual <- z[[y_column]] - z[[x_column]]
    weight <- z$endpoint_multiplicity_q10
    total <- sum(weight)
    data.frame(
      z[1L, groups, drop = FALSE],
      bias = sum(weight * residual) / total,
      rmse = sqrt(sum(weight * residual^2) / total),
      q95 = f7ft_weighted_quantile(abs(residual), weight, 0.95),
      maximum = max(abs(residual), na.rm = TRUE),
      weighted_n = total,
      weighted_n_total = total_weight_all,
      finite_row_fraction = nrow(z) / total_rows,
      weighted_coverage = total / total_weight_all,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$annotation <- sprintf(
    "valid %.1f%%\nbias %.2g\nRMSE %.2g\nq95 %.2g\nmax %.2g",
    100 * out$weighted_coverage,
    out$bias, out$rmse, out$q95, out$maximum
  )
  out
}

f7ft_calibration_theme <- function() {
  ggplot2::theme_classic(base_size = 8.2, base_family = "Helvetica") +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "#F1F1F1", colour = "#BEBEBE", linewidth = 0.3
      ),
      strip.text = ggplot2::element_text(face = "bold", size = 7.5),
      panel.border = ggplot2::element_rect(
        fill = NA, colour = "#555555", linewidth = 0.25
      ),
      plot.title = ggplot2::element_text(face = "bold", size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 7.2, colour = "#555555"),
      legend.position = "bottom"
    )
}

f7ft_calibration_scatter <- function(
    data, x_column, y_column, row_column = NULL, title, subtitle,
    x_label, y_label
) {
  groups <- c(if (!is.null(row_column)) row_column, "context_label")
  annotations <- f7ft_weighted_metrics_for_plot(
    data, x_column, y_column, groups
  )
  plot_data <- data[
    is.finite(data[[x_column]]) & is.finite(data[[y_column]]) &
      is.finite(data$endpoint_multiplicity_q10) &
      data$endpoint_multiplicity_q10 > 0,
    , drop = FALSE
  ]
  mapping <- ggplot2::aes(
    x = .data[[x_column]], y = .data[[y_column]],
    weight = endpoint_multiplicity_q10
  )
  p <- ggplot2::ggplot(plot_data, mapping) +
    ggplot2::geom_bin_2d(bins = 42) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1, colour = "#222222",
      linewidth = 0.35, linetype = "dashed"
    ) +
    ggplot2::geom_text(
      data = annotations,
      ggplot2::aes(x = 1.08, y = 6.92, label = annotation),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      size = 2.05, colour = "#222222", lineheight = 0.92
    ) +
    ggplot2::coord_equal(xlim = c(1, 7), ylim = c(1, 7), expand = FALSE) +
    ggplot2::scale_fill_viridis_c(
      option = "C", trans = "log10", name = "Weighted\nrow count"
    ) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = x_label, y = y_label
    ) +
    f7ft_calibration_theme()
  if (is.null(row_column)) {
    p + ggplot2::facet_grid(. ~ context_label)
  } else {
    p + ggplot2::facet_grid(
      stats::as.formula(paste(row_column, "~ context_label"))
    )
  }
}

f7ft_publish_supp <- function(output, paths, filename) {
  c(
    figures_png = f7r_publish(
      output[["png"]], file.path(paths$figures, paste0(filename, ".png"))
    ),
    figures_pdf = f7r_publish(
      output[["pdf"]], file.path(paths$figures, paste0(filename, ".pdf"))
    ),
    manuscript_png = f7r_publish(
      output[["png"]], file.path(paths$manuscript_figures, paste0(filename, ".png"))
    ),
    manuscript_pdf = f7r_publish(
      output[["pdf"]], file.path(paths$manuscript_figures, paste0(filename, ".pdf"))
    )
  )
}

f7ft_validate_supp_render <- function(
    data_dir, output, published, expected_width, expected_height, figure_id
) {
  f7r_require_packages("magick")
  info <- magick::image_info(magick::image_read(output[["png"]]))
  checks <- data.frame(
    check = c(
      "png_width_px", "png_height_px", "figures_png_md5_match",
      "figures_pdf_md5_match", "manuscript_png_md5_match",
      "manuscript_pdf_md5_match"
    ),
    observed = c(
      info$width[[1L]], info$height[[1L]],
      f7r_md5(output[["png"]]) == f7r_md5(published[["figures_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["figures_pdf"]]),
      f7r_md5(output[["png"]]) == f7r_md5(published[["manuscript_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(expected_width, expected_height, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  checks$passed <- as.character(checks$observed) == as.character(checks$expected)
  path <- f7ft_atomic_write_tsv(
    checks, file.path(data_dir, paste0(figure_id, "_render_validation.tsv"))
  )
  if (!all(checks$passed)) stop(figure_id, " render validation failed.")
  path
}

f7ft_draw_supplement_7_5 <- function(workspace_root = f7r_find_workspace_root()) {
  f7r_require_packages(c("ggplot2", "magick", "scales"))
  paths <- f7r_paths(workspace_root)
  run_paths <- f7ft_paths(paths, run_id = NULL, create = TRUE)
  data <- f7ft_read_diagnostics(run_paths)
  long <- rbind(
    transform(
      data,
      method = "Full eigen finite-time",
      analytical_mean_ploidy = eigen_mean_ploidy
    ),
    transform(
      data,
      method = "Expm finite-time",
      analytical_mean_ploidy = expm_mean_ploidy
    )
  )
  long$method <- factor(
    long$method, levels = c("Full eigen finite-time", "Expm finite-time")
  )
  plot <- f7ft_calibration_scatter(
    long, "analytical_mean_ploidy", "euler_mean_ploidy",
    row_column = "method",
    title = "Finite-time analytical solutions versus Euler numerical integration",
    subtitle = paste0(
      "Same q10 50 optimizer endpoints as Figure 7A-B; density and metrics restore ",
      "optimizer-seed multiplicity"
    ),
    x_label = "Analytical solution mean ploidy",
    y_label = "Euler numerical-integration mean ploidy"
  )
  filename <- "supp_fig7-5_finite_time_eigen_expm_vs_euler"
  output <- f7r_save_plot(
    plot, file.path(run_paths$supp7_5, filename),
    width = 14.0, height = 7.8, dpi = 300
  )
  published <- f7ft_publish_supp(output, paths, filename)
  validation <- f7ft_validate_supp_render(
    run_paths$supp7_5, output, published, 4200L, 2340L, "supp_fig7-5"
  )
  invisible(list(plot = plot, output = output, published = published,
                 validation = validation))
}

f7ft_draw_supplement_7_6 <- function(workspace_root = f7r_find_workspace_root()) {
  f7r_require_packages(c("ggplot2", "patchwork", "magick", "scales"))
  paths <- f7r_paths(workspace_root)
  run_paths <- f7ft_paths(paths, run_id = NULL, create = TRUE)
  data <- f7ft_read_diagnostics(run_paths)
  scatter <- f7ft_calibration_scatter(
    data, "eigen_mean_ploidy", "expm_mean_ploidy",
    title = "A. Full-eigen versus expm finite-time solutions",
    subtitle = "Identity line marks exact agreement",
    x_label = "Full-eigen finite-time mean ploidy",
    y_label = "Expm finite-time mean ploidy"
  )
  data$method_mean <- (data$eigen_mean_ploidy + data$expm_mean_ploidy) / 2
  data$residual <- data$expm_mean_ploidy - data$eigen_mean_ploidy
  residual_metrics <- f7ft_weighted_metrics_for_plot(
    data, "eigen_mean_ploidy", "expm_mean_ploidy", "context_label"
  )
  residual_data <- data[
    is.finite(data$method_mean) & is.finite(data$residual) &
      is.finite(data$endpoint_multiplicity_q10) &
      data$endpoint_multiplicity_q10 > 0,
    , drop = FALSE
  ]
  residual <- ggplot2::ggplot(
    residual_data,
    ggplot2::aes(
      x = method_mean, y = residual, weight = endpoint_multiplicity_q10
    )
  ) +
    ggplot2::geom_bin_2d(bins = 42) +
    ggplot2::geom_hline(
      data = residual_metrics, ggplot2::aes(yintercept = bias),
      colour = "#B2182B", linewidth = 0.42
    ) +
    ggplot2::geom_hline(
      yintercept = 0, colour = "#222222", linewidth = 0.30,
      linetype = "dashed"
    ) +
    ggplot2::facet_grid(. ~ context_label, scales = "free_y") +
    ggplot2::scale_fill_viridis_c(
      option = "C", trans = "log10", name = "Weighted\nrow count"
    ) +
    ggplot2::labs(
      title = "B. Expm minus full-eigen residuals",
      subtitle = paste0(
        "Red line: weighted bias; conditioning and spectral-gap values are ",
        "retained in the source table"
      ),
      x = "Mean of expm and full-eigen ploidy",
      y = "Expm − full-eigen mean ploidy"
    ) +
    f7ft_calibration_theme()
  combined <- scatter / residual + patchwork::plot_layout(heights = c(1.12, 0.88))
  filename <- "supp_fig7-6_eigen_expm_agreement"
  output <- f7r_save_plot(
    combined, file.path(run_paths$supp7_6, filename),
    width = 14.0, height = 9.0, dpi = 300
  )
  published <- f7ft_publish_supp(output, paths, filename)
  validation <- f7ft_validate_supp_render(
    run_paths$supp7_6, output, published, 4200L, 2700L, "supp_fig7-6"
  )
  invisible(list(plot = combined, output = output, published = published,
                 validation = validation))
}

f7ft_draw_supplement_7_7 <- function(workspace_root = f7r_find_workspace_root()) {
  f7r_require_packages(c("ggplot2", "magick", "scales"))
  paths <- f7r_paths(workspace_root)
  run_paths <- f7ft_paths(paths, run_id = NULL, create = TRUE)
  data <- f7ft_read_diagnostics(run_paths)
  data <- data[data$day %in% c(25, 100, 500, 1000), , drop = FALSE]
  data$day_label <- factor(
    paste0("Day ", data$day),
    levels = paste0("Day ", c(25, 100, 500, 1000))
  )
  annotations <- f7ft_weighted_metrics_for_plot(
    data, "steady_mean_ploidy", "expm_mean_ploidy",
    c("context_label", "day_label")
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = steady_mean_ploidy, y = expm_mean_ploidy,
      weight = endpoint_multiplicity_q10
    )
  ) +
    ggplot2::geom_bin_2d(bins = 38) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1, colour = "#222222",
      linewidth = 0.35, linetype = "dashed"
    ) +
    ggplot2::geom_text(
      data = annotations,
      ggplot2::aes(x = 1.08, y = 6.92, label = annotation),
      inherit.aes = FALSE, hjust = 0, vjust = 1,
      size = 1.85, lineheight = 0.90
    ) +
    ggplot2::facet_grid(context_label ~ day_label) +
    ggplot2::coord_equal(xlim = c(1, 7), ylim = c(1, 7), expand = FALSE) +
    ggplot2::scale_fill_viridis_c(
      option = "C", trans = "log10", name = "Weighted\nrow count"
    ) +
    ggplot2::labs(
      title = "Finite-time expm solutions approach the fixed-environment attractor",
      subtitle = paste0(
        "The same q10 endpoint ensemble is shown at four time horizons; identity ",
        "marks equality with the dominant-eigenvector steady state"
      ),
      x = "Dominant-eigenvector steady-state mean ploidy",
      y = "Expm finite-time mean ploidy"
    ) +
    f7ft_calibration_theme()
  filename <- "supp_fig7-7_finite_time_vs_steady_attractor"
  output <- f7r_save_plot(
    plot, file.path(run_paths$supp7_7, filename),
    width = 13.6, height = 11.8, dpi = 300
  )
  published <- f7ft_publish_supp(output, paths, filename)
  validation <- f7ft_validate_supp_render(
    run_paths$supp7_7, output, published, 4080L, 3540L, "supp_fig7-7"
  )
  invisible(list(plot = plot, output = output, published = published,
                 validation = validation))
}
