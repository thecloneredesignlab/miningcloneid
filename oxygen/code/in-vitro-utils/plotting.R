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
  if (is.null(observed_kary_df) && (is.null(summary_df) || !"observed_mean_kary_N" %in% names(summary_df))) {
    stop("observed_kary_df or summary_df$observed_mean_kary_N is required for ivt_plot_lineage_ploidy().")
  }

  ploidy_lines <- quantile_df |>
    dplyr::distinct(
      dplyr::across(dplyr::any_of(c(
        "cohort",
        "lineage_label",
        "lineage_terminal_key",
        "lineage_passage_index",
        "passage_index",
        "oxygen_pct",
        "quantile_prob",
        "predicted_quantile_kary_N"
      )))
    ) |>
    dplyr::filter(is.finite(quantile_prob), is.finite(predicted_quantile_kary_N)) |>
    dplyr::mutate(
      predicted_kary_N = predicted_quantile_kary_N,
      fit = primary_label
    )
  if (!nrow(ploidy_lines) && !is.null(summary_df) && "predicted_mean_kary_N" %in% names(summary_df)) {
    ploidy_lines <- summary_df |>
      dplyr::distinct(
        dplyr::across(dplyr::any_of(c("cohort", "lineage_label", "lineage_terminal_key", "lineage_passage_index", "passage_index", "oxygen_pct", "predicted_mean_kary_N")))
      ) |>
      dplyr::filter(is.finite(predicted_mean_kary_N)) |>
      dplyr::mutate(
        quantile_prob = 0.5,
        predicted_kary_N = predicted_mean_kary_N,
        fit = primary_label
      )
  }
  x_var <- if ("lineage_passage_index" %in% names(ploidy_lines)) "lineage_passage_index" else "passage_index"
  has_lineage_label <- "lineage_label" %in% names(ploidy_lines) &&
    any(!is.na(ploidy_lines$lineage_label)) &&
    (!"cohort" %in% names(ploidy_lines) ||
      any(as.character(ploidy_lines$lineage_label) != as.character(ploidy_lines$cohort), na.rm = TRUE))
  facet_var <- if (has_lineage_label) "lineage_label" else NULL
  facet_formula <- if (!is.null(facet_var) && all(c("cohort", facet_var) %in% names(ploidy_lines))) {
    stats::as.formula(paste("cohort ~", facet_var))
  } else {
    NULL
  }
  facet_wrap_formula <- if (is.null(facet_formula) && "cohort" %in% names(ploidy_lines)) {
    stats::as.formula("~ cohort")
  } else {
    NULL
  }

  line_group_vars <- c("fit", intersect("cohort", names(ploidy_lines)))
  if (!is.null(facet_var)) {
    line_group_vars <- c(line_group_vars, facet_var)
  }
  line_group_vars <- c(line_group_vars, "quantile_prob")
  line_group_vars <- unique(line_group_vars)
  ploidy_lines <- ploidy_lines |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(line_group_vars, x_var)))) |>
    dplyr::summarise(predicted_kary_N = mean(predicted_kary_N, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(is.finite(predicted_kary_N))
  ploidy_lines$line_group <- do.call(
    interaction,
    c(ploidy_lines[line_group_vars], list(drop = TRUE, lex.order = TRUE))
  )

  ploidy_obs <- if (!is.null(observed_kary_df) && "observed_kary_N" %in% names(observed_kary_df)) {
    observed_kary_df |>
      dplyr::filter(is.finite(observed_kary_N)) |>
      dplyr::mutate(observed_plot_kary_N = observed_kary_N)
  } else {
    summary_df |>
      dplyr::filter(is.finite(observed_mean_kary_N)) |>
      dplyr::mutate(observed_plot_kary_N = observed_mean_kary_N)
  }

  fit_levels <- unique(ploidy_lines$fit)
  fit_palette <- stats::setNames(c("#1b9e77", "#377eb8", "#4e79a7", "#e15759", "#f28e2b")[seq_along(fit_levels)], fit_levels)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = ploidy_lines,
      ggplot2::aes(.data[[x_var]], predicted_kary_N, color = fit, group = line_group),
      linewidth = 0.7,
      alpha = quantile_alpha
    ) +
    ggplot2::geom_point(
      data = ploidy_obs,
      ggplot2::aes(.data[[x_var]], observed_plot_kary_N),
      size = 1.8,
      alpha = 0.8,
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
  } else if (!is.null(facet_wrap_formula)) {
    p <- p + ggplot2::facet_wrap(facet_wrap_formula, scales = "free_x")
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
  has_lineage_label <- "lineage_label" %in% names(pred_lines) &&
    any(!is.na(pred_lines$lineage_label)) &&
    (!"cohort" %in% names(pred_lines) ||
      any(as.character(pred_lines$lineage_label) != as.character(pred_lines$cohort), na.rm = TRUE))
  facet_var <- if (has_lineage_label) "lineage_label" else NULL
  facet_formula <- if (!is.null(facet_var) && all(c("cohort", facet_var) %in% names(pred_lines))) {
    stats::as.formula(paste("cohort ~", facet_var))
  } else {
    NULL
  }
  facet_wrap_formula <- if (is.null(facet_formula) && "cohort" %in% names(pred_lines)) {
    stats::as.formula("~ cohort")
  } else {
    NULL
  }

  line_group_vars <- c("fit", intersect("cohort", names(pred_lines)))
  if (!is.null(facet_var)) {
    line_group_vars <- c(line_group_vars, facet_var)
  }
  if (!is.null(facet_var) && "lineage_terminal_key" %in% names(pred_lines)) {
    group_check_vars <- c(intersect("cohort", names(pred_lines)), facet_var, "lineage_terminal_key")
    group_check <- pred_lines |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_check_vars))) |>
      dplyr::summarise(n_x = dplyr::n_distinct(.data[[x_var]]), .groups = "drop")
    if (any(group_check$n_x > 1L)) {
      line_group_vars <- c(line_group_vars, "lineage_terminal_key")
    }
  }
  line_group_vars <- unique(line_group_vars)
  pred_lines <- pred_lines |>
    dplyr::filter(is.finite(predicted_growth_rate)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(line_group_vars, x_var)))) |>
    dplyr::summarise(predicted_growth_rate = mean(predicted_growth_rate, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(is.finite(predicted_growth_rate))
  pred_lines$line_group <- do.call(
    interaction,
    c(pred_lines[line_group_vars], list(drop = TRUE, lex.order = TRUE))
  )

  fit_levels <- unique(pred_lines$fit)
  fit_palette <- stats::setNames(c("#1b9e77", "#377eb8", "#4e79a7", "#e15759", "#f28e2b")[seq_along(fit_levels)], fit_levels)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = pred_lines,
      ggplot2::aes(.data[[x_var]], predicted_growth_rate, color = fit, group = line_group),
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
  } else if (!is.null(facet_wrap_formula)) {
    p <- p + ggplot2::facet_wrap(facet_wrap_formula, scales = "free_x")
  }

  p
}

ivt_plot_distribution_heatmap <- function(dist_df, max_N = 110) {
  plot_df <- dist_df |>
    dplyr::filter(is.finite(.data$N), is.finite(.data$fraction), .data$N <= max_N, .data$fraction >= 0)
  x_var <- if ("lineage_passage_index" %in% names(plot_df)) "lineage_passage_index" else "passage_index"
  plot_df <- plot_df |>
    dplyr::filter(is.finite(.data[[x_var]])) |>
    dplyr::mutate(N = as.integer(round(.data$N)))
  if (!nrow(plot_df)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No predicted ploidy distribution rows"))
  }

  has_lineage_label <- "lineage_label" %in% names(plot_df) &&
    any(!is.na(plot_df$lineage_label)) &&
    (!"cohort" %in% names(plot_df) ||
      any(as.character(plot_df$lineage_label) != as.character(plot_df$cohort), na.rm = TRUE))
  facet_var <- if (has_lineage_label) "lineage_label" else NULL
  facet_vars <- intersect(c("cohort", facet_var), names(plot_df))
  group_vars <- unique(c(facet_vars, x_var, "N"))
  plot_df <- plot_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(fraction = mean(.data$fraction, na.rm = TRUE), .groups = "drop")

  n_seq <- seq.int(
    max(0L, floor(min(plot_df$N, na.rm = TRUE))),
    min(as.integer(max_N), ceiling(max(plot_df$N, na.rm = TRUE)))
  )
  split_key <- if (length(facet_vars)) {
    do.call(interaction, c(plot_df[facet_vars], list(drop = TRUE, lex.order = TRUE)))
  } else {
    factor(rep("all", nrow(plot_df)))
  }
  completed <- lapply(split(plot_df, split_key), function(df_group) {
    x_vals <- sort(unique(df_group[[x_var]]))
    grid <- expand.grid(
      x_value = x_vals,
      N = n_seq,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    names(grid)[names(grid) == "x_value"] <- x_var
    for (var in facet_vars) {
      grid[[var]] <- df_group[[var]][[1]]
    }
    merged <- merge(grid, df_group, by = unique(c(facet_vars, x_var, "N")), all.x = TRUE, sort = FALSE)
    merged$fraction[is.na(merged$fraction)] <- 0
    merged
  })
  plot_df <- do.call(rbind, completed)
  plot_df$x_plot <- factor(
    as.character(plot_df[[x_var]]),
    levels = as.character(sort(unique(plot_df[[x_var]])))
  )
  fraction_fill_max <- max(plot_df$fraction, na.rm = TRUE)
  if (!is.finite(fraction_fill_max) || fraction_fill_max <= 0) {
    fraction_fill_max <- 1
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x_plot, N, fill = fraction)) +
    ggplot2::geom_tile(width = 0.98, height = 1) +
    ggplot2::scale_y_continuous(limits = c(min(n_seq) - 0.5, max(n_seq) + 0.5), expand = c(0, 0)) +
    ggplot2::scale_x_discrete(drop = TRUE, expand = c(0, 0)) +
    ggplot2::scale_fill_gradientn(
      colors = c("#f7f7f7", "#2c7fb8", "#ffff33"),
      values = scales::rescale(c(0, 0.05, 1)),
      limits = c(0, fraction_fill_max)
    ) +
    ggplot2::labs(
      title = "Predicted ploidy distribution across passages",
      x = "Passage index",
      y = "Chromosome-count state N",
      fill = "Fraction"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.null(facet_var) && all(c("cohort", facet_var) %in% names(plot_df))) {
    p <- p + ggplot2::facet_grid(
      stats::as.formula(paste("cohort ~", facet_var)),
      scales = "free_x",
      space = "free_x"
    )
  } else if ("cohort" %in% names(plot_df)) {
    p <- p + ggplot2::facet_wrap(~cohort, scales = "free_x")
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
