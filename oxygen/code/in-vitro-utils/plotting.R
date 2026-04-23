ivt_plot_lineage_counts <- function(summary_df) {
  count_df <- summary_df |>
    dplyr::distinct(passage_index, oxygen_pct, predicted_live_cells, selected_day)

  ggplot2::ggplot(count_df, ggplot2::aes(passage_index, predicted_live_cells, color = factor(oxygen_pct))) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      title = "Predicted live cells at the selected propagation day",
      x = "Passage index",
      y = "Predicted live cells",
      color = "Oxygen (%)"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

ivt_plot_lineage_ploidy <- function(summary_df,
                                    quantile_df = NULL,
                                    observed_kary_df = NULL,
                                    primary_label = "Loaded best",
                                    quantile_alpha = 0.5) {
  if (is.null(quantile_df)) {
    stop("quantile_df is required for ivt_plot_lineage_ploidy().")
  }
  if (is.null(observed_kary_df)) {
    stop("observed_kary_df is required for ivt_plot_lineage_ploidy().")
  }

  ploidy_lines <- quantile_df |>
    dplyr::distinct(
      dplyr::across(dplyr::any_of(c("cohort", "lineage_label", "lineage_terminal_key", "lineage_passage_index", "passage_index", "oxygen_pct", "quantile_prob", "predicted_quantile_kary_N")))
    ) |>
    dplyr::filter(is.finite(predicted_quantile_kary_N)) |>
    dplyr::mutate(fit = primary_label)
  ploidy_obs <- observed_kary_df |>
    dplyr::filter(is.finite(observed_kary_N))

  x_var <- if ("lineage_passage_index" %in% names(ploidy_lines)) "lineage_passage_index" else "passage_index"
  facet_var <- if ("lineage_label" %in% names(ploidy_lines)) "lineage_label" else "lineage_terminal_key"
  group_var <- if ("lineage_terminal_key" %in% names(ploidy_lines)) "lineage_terminal_key" else x_var
  facet_formula <- if (all(c("cohort", facet_var) %in% names(ploidy_lines))) {
    stats::as.formula(paste("cohort ~", facet_var))
  } else {
    NULL
  }

  fit_levels <- unique(ploidy_lines$fit)
  fit_palette <- stats::setNames(c("#1b9e77", "#377eb8", "#4e79a7", "#e15759", "#f28e2b")[seq_along(fit_levels)], fit_levels)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = ploidy_lines,
      ggplot2::aes(.data[[x_var]], predicted_quantile_kary_N, color = fit, group = interaction(fit, quantile_prob, .data[[group_var]])),
      linewidth = 0.5,
      alpha = quantile_alpha
    ) +
    ggplot2::geom_point(
      data = ploidy_obs,
      ggplot2::aes(.data[[x_var]], observed_kary_N),
      size = 0.8,
      alpha = 0.55,
      color = "#d95f02",
      position = ggplot2::position_jitter(width = 0.12, height = 0)
    ) +
    ggplot2::scale_color_manual(values = fit_palette) +
    ggplot2::labs(
      title = "Predicted chromosome-count quantiles versus observed cells by passage",
      x = "Passage index",
      y = "Chromosome count (N)",
      color = "Predicted fit"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.null(facet_formula)) {
    p <- p + ggplot2::facet_grid(facet_formula, scales = "free_x", space = "free_x")
  }

  p
}

ivt_plot_daily_counts <- function(count_df) {
  chosen <- count_df |>
    dplyr::filter(.data$day == .data$selected_day)

  ggplot2::ggplot(count_df, ggplot2::aes(day, live_cells, group = passage_index, color = factor(oxygen_pct))) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(data = chosen, size = 2) +
    ggplot2::facet_wrap(~ passage_index, scales = "free_y") +
    ggplot2::labs(
      title = "Daily live-cell trajectories within each passage",
      x = "Day within passage",
      y = "Live cells",
      color = "Oxygen (%)"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

ivt_plot_lineage_growth <- function(summary_df,
                                    comparison_summary_df = NULL,
                                    primary_label = "Loaded best",
                                    comparison_label = "Comparison") {
  pred <- summary_df |>
    dplyr::distinct(
      dplyr::across(dplyr::any_of(c("cohort", "lineage_label", "lineage_terminal_key", "lineage_passage_index", "passage_index", "oxygen_pct", "predicted_growth_rate")))
    ) |>
    dplyr::mutate(fit = primary_label)
  obs <- summary_df |>
    dplyr::filter(is.finite(observed_growth))
  pred_comp <- if (is.null(comparison_summary_df)) {
    NULL
  } else {
    comparison_summary_df |>
      dplyr::distinct(
        dplyr::across(dplyr::any_of(c("cohort", "lineage_label", "lineage_terminal_key", "lineage_passage_index", "passage_index", "oxygen_pct", "predicted_growth_rate")))
      ) |>
      dplyr::mutate(fit = comparison_label)
  }
  pred_lines <- dplyr::bind_rows(pred, pred_comp)
  x_var <- if ("lineage_passage_index" %in% names(pred_lines)) "lineage_passage_index" else "passage_index"
  facet_var <- if ("lineage_label" %in% names(pred_lines)) "lineage_label" else "lineage_terminal_key"
  group_var <- if ("lineage_terminal_key" %in% names(pred_lines)) "lineage_terminal_key" else x_var
  facet_formula <- if (all(c("cohort", facet_var) %in% names(pred_lines))) {
    stats::as.formula(paste("cohort ~", facet_var))
  } else {
    NULL
  }

  fit_levels <- unique(pred_lines$fit)
  fit_palette <- stats::setNames(c("#1b9e77", "#377eb8", "#4e79a7", "#e15759", "#f28e2b")[seq_along(fit_levels)], fit_levels)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = pred_lines,
      ggplot2::aes(.data[[x_var]], predicted_growth_rate, color = fit, group = interaction(fit, .data[[group_var]])),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = pred_lines,
      ggplot2::aes(.data[[x_var]], predicted_growth_rate, color = fit),
      size = 2
    ) +
    ggplot2::geom_point(
      data = obs,
      ggplot2::aes(.data[[x_var]], observed_growth),
      size = 2,
      alpha = 0.8,
      color = "#7570b3",
      position = ggplot2::position_jitter(width = 0.08, height = 0)
    ) +
    ggplot2::scale_color_manual(values = fit_palette) +
    ggplot2::labs(
      title = "Predicted versus observed growth rate by passage",
      x = "Passage index",
      y = "Growth rate",
      color = "Predicted fit"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.null(facet_formula)) {
    p <- p + ggplot2::facet_grid(facet_formula, scales = "free_x", space = "free_x")
  }

  p
}

ivt_plot_distribution_heatmap <- function(dist_df, max_N = 110) {
  plot_df <- dist_df |>
    dplyr::filter(.data$N <= max_N, .data$fraction > 0)
  y_var <- if ("lineage_passage_index" %in% names(plot_df)) "lineage_passage_index" else "passage_index"

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(N, .data[[y_var]], fill = fraction)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(
      title = "Predicted ploidy distribution across passages",
      x = "Chromosome-count state N",
      y = "Passage index",
      fill = "Fraction"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  facet_var <- if ("lineage_label" %in% names(plot_df)) "lineage_label" else "lineage_terminal_key"
  if (all(c("cohort", facet_var) %in% names(plot_df))) {
    p <- p + ggplot2::facet_grid(
      stats::as.formula(paste("cohort ~", facet_var)),
      scales = "free_y",
      space = "free_y"
    )
  }

  p
}

ivt_plot_lineage_flow_density <- function(flow_overlay_df,
                                          max_facets = 20L) {
  plot_df <- flow_overlay_df |>
    dplyr::filter(is.finite(ploidy), is.finite(density), density >= 0)
  if (!nrow(plot_df)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No matched flow-density observations available"))
  }

  plot_df <- plot_df |>
    dplyr::mutate(
      lineage_display = dplyr::coalesce(.data$lineage_label, .data$lineage_terminal_key, .data$sample_name, .data$passage_id),
      passage_display = dplyr::coalesce(.data$lineage_passage_index, .data$passage_index),
      facet_label = paste0(lineage_display, " | p", passage_display)
    )

  facet_levels <- plot_df |>
    dplyr::distinct(.data$cohort, .data$facet_label, .data$passage_display) |>
    dplyr::arrange(.data$cohort, .data$passage_display, .data$facet_label) |>
    dplyr::slice_head(n = max(1L, as.integer(max_facets))) |>
    dplyr::pull(.data$facet_label)

  plot_df <- plot_df |>
    dplyr::filter(.data$facet_label %in% facet_levels)

  ggplot2::ggplot(plot_df, ggplot2::aes(ploidy, density, color = series)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_grid(cohort ~ facet_label, scales = "free_y", space = "free_x") +
    ggplot2::scale_color_manual(values = c("Observed" = "#d95f02", "Predicted" = "#1b9e77")) +
    ggplot2::labs(
      title = "Observed versus predicted G0/G1 ploidy-density curves",
      x = "Ploidy",
      y = "Density",
      color = "Series"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(size = 6)
    )
}
