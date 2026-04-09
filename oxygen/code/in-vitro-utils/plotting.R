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

ivt_plot_lineage_ploidy <- function(summary_df, comparison_summary_df = NULL) {
  ploidy_pred <- summary_df |>
    dplyr::distinct(
      dplyr::across(dplyr::any_of(c("cohort", "lineage_label", "lineage_terminal_key", "lineage_passage_index", "passage_index", "oxygen_pct", "predicted_mean_kary_N")))
    ) |>
    dplyr::mutate(fit = "Loaded best")
  ploidy_obs <- summary_df |>
    dplyr::filter(is.finite(observed_mean_kary_N))
  ploidy_comp <- if (is.null(comparison_summary_df)) {
    NULL
  } else {
    comparison_summary_df |>
      dplyr::distinct(
        dplyr::across(dplyr::any_of(c("cohort", "lineage_label", "lineage_terminal_key", "lineage_passage_index", "passage_index", "oxygen_pct", "predicted_mean_kary_N")))
      ) |>
      dplyr::mutate(fit = "Optimized")
  }
  ploidy_lines <- dplyr::bind_rows(ploidy_pred, ploidy_comp)
  x_var <- if ("lineage_passage_index" %in% names(ploidy_lines)) "lineage_passage_index" else "passage_index"
  facet_var <- if ("lineage_label" %in% names(ploidy_lines)) "lineage_label" else "lineage_terminal_key"
  group_var <- if ("lineage_terminal_key" %in% names(ploidy_lines)) "lineage_terminal_key" else x_var
  facet_formula <- if (all(c("cohort", facet_var) %in% names(ploidy_lines))) {
    stats::as.formula(paste("cohort ~", facet_var))
  } else {
    NULL
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = ploidy_lines,
      ggplot2::aes(.data[[x_var]], predicted_mean_kary_N, color = fit, group = interaction(fit, .data[[group_var]])),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = ploidy_lines,
      ggplot2::aes(.data[[x_var]], predicted_mean_kary_N, color = fit),
      size = 2
    ) +
    ggplot2::geom_point(
      data = ploidy_obs,
      ggplot2::aes(.data[[x_var]], observed_mean_kary_N),
      size = 2,
      alpha = 0.8,
      color = "#d95f02",
      position = ggplot2::position_jitter(width = 0.08, height = 0)
    ) +
    ggplot2::scale_color_manual(values = c("Loaded best" = "#1b9e77", "Optimized" = "#377eb8")) +
    ggplot2::labs(
      title = "Predicted versus observed mean chromosome count by passage",
      x = "Passage index",
      y = "Mean chromosome count (N)",
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

ivt_plot_lineage_growth <- function(summary_df, comparison_summary_df = NULL) {
  pred <- summary_df |>
    dplyr::distinct(
      dplyr::across(dplyr::any_of(c("cohort", "lineage_label", "lineage_terminal_key", "lineage_passage_index", "passage_index", "oxygen_pct", "predicted_growth_rate")))
    ) |>
    dplyr::mutate(fit = "Loaded best")
  obs <- summary_df |>
    dplyr::filter(is.finite(observed_growth))
  pred_comp <- if (is.null(comparison_summary_df)) {
    NULL
  } else {
    comparison_summary_df |>
      dplyr::distinct(
        dplyr::across(dplyr::any_of(c("cohort", "lineage_label", "lineage_terminal_key", "lineage_passage_index", "passage_index", "oxygen_pct", "predicted_growth_rate")))
      ) |>
      dplyr::mutate(fit = "Optimized")
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
    ggplot2::scale_color_manual(values = c("Loaded best" = "#1b9e77", "Optimized" = "#377eb8")) +
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
