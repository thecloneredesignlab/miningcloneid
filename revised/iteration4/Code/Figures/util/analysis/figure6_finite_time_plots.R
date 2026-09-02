#!/usr/bin/env Rscript

# Rendering utilities for Figure 6 finite-time q10 panels and diagnostics.

options(stringsAsFactors = FALSE, warn = 1)

f6ft_read_panel_object <- function(run_paths, panel_letter) {
  path <- file.path(
    run_paths$run_root,
    paste0("finite_time_panel_", tolower(panel_letter), ".rds")
  )
  f6r_require_files(path, paste0("Figure 6", panel_letter, " finite-time data"))
  object <- readRDS(path)
  if (!identical(object$profile, f6ft_profile()) ||
      !identical(object$panel_letter, panel_letter)) {
    stop("Unexpected finite-time panel object: ", path)
  }
  object
}

f6ft_panel_long <- function(object) {
  values <- object$mean_ploidy
  expected_dim <- c(
    length(object$initial_ploidy), length(object$day_values),
    length(object$o2_values), length(object$effective_p_misseg)
  )
  if (!identical(as.integer(dim(values)), as.integer(expected_dim))) {
    stop("Finite-time panel array dimensions do not match its coordinate grids.")
  }
  n_initial <- length(object$initial_ploidy)
  n_day <- length(object$day_values)
  n_o2 <- length(object$o2_values)
  n_p <- length(object$effective_p_misseg)
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
    effective_p_misseg = rep(
      object$effective_p_misseg, each = n_initial * n_day * n_o2
    ),
    mean_ploidy = as.vector(values),
    stringsAsFactors = FALSE
  )
}

f6ft_format_p <- function(x) {
  ifelse(abs(x - 0.005) < 1e-12, "0.005", sprintf("%.2f", x))
}

f6ft_panel_title <- function(object) {
  paste0(
    object$panel_letter, ". ", object$pair_label, " — ", object$model_context
  )
}

f6ft_finite_time_panel_plot <- function(object, show_legend = TRUE) {
  f6r_require_packages(c("ggplot2", "scales"))
  data <- f6ft_panel_long(object)
  data$initial_label <- factor(
    paste0(data$initial_ploidy, "N"),
    levels = paste0(f6ft_initial_ploidy(), "N")
  )
  data$p_label <- factor(
    f6ft_format_p(data$effective_p_misseg),
    levels = f6ft_format_p(f6ft_p_values())
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = day, y = O2_pct, fill = mean_ploidy)
  ) +
    ggplot2::geom_raster(interpolate = FALSE) +
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
      breaks = c(0, 500, 1000), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 2.5, 5), expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = f6ft_panel_title(object),
      subtitle = paste0(
        "Rows: initial ploidy; columns: fixed effective missegregation probability; ",
        "daily expm propagation"
      ),
      x = "Time (day)", y = "Fixed oxygen (%)"
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

f6ft_figure6_validation <- function(
    paths, run_paths, output, published, panel_objects, image_info
) {
  panels <- names(panel_objects)
  checks <- data.frame(
    check = c(
      "main_panel_count", "finite_time_panel_count",
      "finite_time_each_has_5_initial_ploidies",
      "finite_time_each_has_5_p_values",
      "finite_time_each_has_201_o2_values",
      "finite_time_each_has_1001_days",
      "finite_time_weights_equal_50",
      "png_width_px", "png_height_px",
      "figures_png_md5_match", "figures_pdf_md5_match",
      "manuscript_png_md5_match", "manuscript_pdf_md5_match"
    ),
    observed = c(
      6L, length(panels),
      all(vapply(panel_objects, function(x) length(x$initial_ploidy) == 5L, logical(1L))),
      all(vapply(panel_objects, function(x) length(x$effective_p_misseg) == 5L, logical(1L))),
      all(vapply(panel_objects, function(x) length(x$o2_values) == 201L, logical(1L))),
      all(vapply(panel_objects, function(x) length(x$day_values) == 1001L, logical(1L))),
      all(vapply(panel_objects, function(x) all(x$optimizer_endpoint_weight == 50L), logical(1L))),
      image_info$width[[1L]], image_info$height[[1L]],
      f6r_md5(output[["png"]]) == f6r_md5(published[["figures_png"]]),
      f6r_md5(output[["pdf"]]) == f6r_md5(published[["figures_pdf"]]),
      f6r_md5(output[["png"]]) == f6r_md5(published[["manuscript_png"]]),
      f6r_md5(output[["pdf"]]) == f6r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(
      6L, 4L, TRUE, TRUE, TRUE, TRUE, TRUE,
      5040L, 7560L, TRUE, TRUE, TRUE, TRUE
    ),
    stringsAsFactors = FALSE
  )
  checks$passed <- as.character(checks$observed) == as.character(checks$expected)
  path <- f6ft_atomic_write_tsv(
    checks, file.path(run_paths$run_root, "figure6_finite_time_render_validation.tsv")
  )
  if (!all(checks$passed)) stop("Rendered Figure 6 finite-time validation failed.")
  path
}

f6ft_draw_main <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c(
    "ggplot2", "patchwork", "magick", "isoband", "scales"
  ))
  paths <- f6r_paths(workspace_root)
  run_paths <- f6ft_paths(paths, run_id = NULL, create = TRUE)
  panel_objects <- stats::setNames(
    lapply(c("C", "D", "E", "F"), function(letter) {
      f6ft_read_panel_object(run_paths, letter)
    }),
    c("C", "D", "E", "F")
  )
  p_a <- f6x_main_surface_plot(paths) +
    ggplot2::theme(legend.position = "bottom")
  p_b <- f6x_main_inverse_plot(paths) +
    ggplot2::theme(legend.position = "bottom")
  finite_plots <- lapply(panel_objects, f6ft_finite_time_panel_plot)

  panel_paths <- c()
  for (letter in names(finite_plots)) {
    panel_paths <- c(
      panel_paths,
      stats::setNames(
        f6r_save_plot(
          finite_plots[[letter]],
          file.path(
            run_paths$panels,
            paste0("figure6", tolower(letter), "_finite_time_q10")
          ),
          width = 8.2, height = 8.8
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
  combined <- top / finite_block +
    patchwork::plot_layout(heights = c(1, 2.05)) +
    patchwork::plot_annotation(
      caption = paste0(
        "C-F show finite-time expm solutions averaged across the same q10 50 ",
        "optimizer endpoints used for A-B; repeated exact parameter endpoints ",
        "retain their optimizer-seed multiplicity."
      )
    ) &
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        size = 7, colour = "#555555", hjust = 0
      )
    )

  output <- f6r_save_plot(
    combined, file.path(run_paths$rendered, "assembled_fig6"),
    width = 16.8, height = 25.2, dpi = 300
  )
  published <- c(
    figures_png = f6r_publish(
      output[["png"]], file.path(paths$figures, "assembled_fig6.png")
    ),
    figures_pdf = f6r_publish(
      output[["pdf"]], file.path(paths$figures, "assembled_fig6.pdf")
    ),
    manuscript_png = f6r_publish(
      output[["png"]], file.path(paths$manuscript_figures, "assembled_fig6.png")
    ),
    manuscript_pdf = f6r_publish(
      output[["pdf"]], file.path(paths$manuscript_figures, "assembled_fig6.pdf")
    )
  )
  info <- magick::image_info(magick::image_read(output[["png"]]))
  validation <- f6ft_figure6_validation(
    paths, run_paths, output, published, panel_objects, info
  )
  invisible(list(
    plot = combined, output = output, published = published,
    panel_paths = panel_paths, validation = validation
  ))
}

f6ft_read_diagnostics <- function(run_paths) {
  path <- file.path(run_paths$run_root, "finite_time_diagnostic_rows.rds")
  f6r_require_files(path, "Figure 6 finite-time diagnostic rows")
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

f6ft_weighted_metrics_for_plot <- function(data, x_column, y_column, groups) {
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
      q95 = f6ft_weighted_quantile(abs(residual), weight, 0.95),
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

f6ft_calibration_theme <- function() {
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

f6ft_calibration_scatter <- function(
    data, x_column, y_column, row_column = NULL, title, subtitle,
    x_label, y_label
) {
  groups <- c(if (!is.null(row_column)) row_column, "context_label")
  annotations <- f6ft_weighted_metrics_for_plot(
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
    f6ft_calibration_theme()
  if (is.null(row_column)) {
    p + ggplot2::facet_grid(. ~ context_label)
  } else {
    p + ggplot2::facet_grid(
      stats::as.formula(paste(row_column, "~ context_label"))
    )
  }
}

f6ft_publish_supp <- function(output, paths, filename) {
  c(
    figures_png = f6r_publish(
      output[["png"]], file.path(paths$figures, paste0(filename, ".png"))
    ),
    figures_pdf = f6r_publish(
      output[["pdf"]], file.path(paths$figures, paste0(filename, ".pdf"))
    ),
    manuscript_png = f6r_publish(
      output[["png"]], file.path(paths$manuscript_figures, paste0(filename, ".png"))
    ),
    manuscript_pdf = f6r_publish(
      output[["pdf"]], file.path(paths$manuscript_figures, paste0(filename, ".pdf"))
    )
  )
}

f6ft_validate_supp_render <- function(
    data_dir, output, published, expected_width, expected_height, figure_id
) {
  f6r_require_packages("magick")
  info <- magick::image_info(magick::image_read(output[["png"]]))
  checks <- data.frame(
    check = c(
      "png_width_px", "png_height_px", "figures_png_md5_match",
      "figures_pdf_md5_match", "manuscript_png_md5_match",
      "manuscript_pdf_md5_match"
    ),
    observed = c(
      info$width[[1L]], info$height[[1L]],
      f6r_md5(output[["png"]]) == f6r_md5(published[["figures_png"]]),
      f6r_md5(output[["pdf"]]) == f6r_md5(published[["figures_pdf"]]),
      f6r_md5(output[["png"]]) == f6r_md5(published[["manuscript_png"]]),
      f6r_md5(output[["pdf"]]) == f6r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(expected_width, expected_height, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  checks$passed <- as.character(checks$observed) == as.character(checks$expected)
  path <- f6ft_atomic_write_tsv(
    checks, file.path(data_dir, paste0(figure_id, "_render_validation.tsv"))
  )
  if (!all(checks$passed)) stop(figure_id, " render validation failed.")
  path
}

f6ft_draw_supplement_6_5 <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c("ggplot2", "magick", "scales"))
  paths <- f6r_paths(workspace_root)
  run_paths <- f6ft_paths(paths, run_id = NULL, create = TRUE)
  data <- f6ft_read_diagnostics(run_paths)
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
  plot <- f6ft_calibration_scatter(
    long, "analytical_mean_ploidy", "euler_mean_ploidy",
    row_column = "method",
    title = "Finite-time analytical solutions versus Euler numerical integration",
    subtitle = paste0(
      "Same q10 50 optimizer endpoints as Figure 6A-B; density and metrics restore ",
      "optimizer-seed multiplicity"
    ),
    x_label = "Analytical solution mean ploidy",
    y_label = "Euler numerical-integration mean ploidy"
  )
  filename <- "supp_fig6-5_finite_time_eigen_expm_vs_euler"
  output <- f6r_save_plot(
    plot, file.path(run_paths$supp6_5, filename),
    width = 14.0, height = 7.8, dpi = 300
  )
  published <- f6ft_publish_supp(output, paths, filename)
  validation <- f6ft_validate_supp_render(
    run_paths$supp6_5, output, published, 4200L, 2340L, "supp_fig6-5"
  )
  invisible(list(plot = plot, output = output, published = published,
                 validation = validation))
}

f6ft_draw_supplement_6_6 <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c("ggplot2", "patchwork", "magick", "scales"))
  paths <- f6r_paths(workspace_root)
  run_paths <- f6ft_paths(paths, run_id = NULL, create = TRUE)
  data <- f6ft_read_diagnostics(run_paths)
  scatter <- f6ft_calibration_scatter(
    data, "eigen_mean_ploidy", "expm_mean_ploidy",
    title = "A. Full-eigen versus expm finite-time solutions",
    subtitle = "Identity line marks exact agreement",
    x_label = "Full-eigen finite-time mean ploidy",
    y_label = "Expm finite-time mean ploidy"
  )
  data$method_mean <- (data$eigen_mean_ploidy + data$expm_mean_ploidy) / 2
  data$residual <- data$expm_mean_ploidy - data$eigen_mean_ploidy
  residual_metrics <- f6ft_weighted_metrics_for_plot(
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
    f6ft_calibration_theme()
  combined <- scatter / residual + patchwork::plot_layout(heights = c(1.12, 0.88))
  filename <- "supp_fig6-6_eigen_expm_agreement"
  output <- f6r_save_plot(
    combined, file.path(run_paths$supp6_6, filename),
    width = 14.0, height = 9.0, dpi = 300
  )
  published <- f6ft_publish_supp(output, paths, filename)
  validation <- f6ft_validate_supp_render(
    run_paths$supp6_6, output, published, 4200L, 2700L, "supp_fig6-6"
  )
  invisible(list(plot = combined, output = output, published = published,
                 validation = validation))
}

f6ft_draw_supplement_6_7 <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c("ggplot2", "magick", "scales"))
  paths <- f6r_paths(workspace_root)
  run_paths <- f6ft_paths(paths, run_id = NULL, create = TRUE)
  data <- f6ft_read_diagnostics(run_paths)
  data <- data[data$day %in% c(25, 100, 500, 1000), , drop = FALSE]
  data$day_label <- factor(
    paste0("Day ", data$day),
    levels = paste0("Day ", c(25, 100, 500, 1000))
  )
  annotations <- f6ft_weighted_metrics_for_plot(
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
    f6ft_calibration_theme()
  filename <- "supp_fig6-7_finite_time_vs_steady_attractor"
  output <- f6r_save_plot(
    plot, file.path(run_paths$supp6_7, filename),
    width = 13.6, height = 11.8, dpi = 300
  )
  published <- f6ft_publish_supp(output, paths, filename)
  validation <- f6ft_validate_supp_render(
    run_paths$supp6_7, output, published, 4080L, 3540L, "supp_fig6-7"
  )
  invisible(list(plot = plot, output = output, published = published,
                 validation = validation))
}
