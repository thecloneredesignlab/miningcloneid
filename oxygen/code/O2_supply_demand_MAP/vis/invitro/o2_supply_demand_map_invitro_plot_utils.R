# Pure in-vitro ggplot constructors shared by in-vitro visualization entrypoints.

.ivt_plot_utils_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[
    basename(frame_files) == "o2_supply_demand_map_invitro_plot_utils.R"
  ]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(
  file.path(.ivt_plot_utils_dir, "..", "o2_supply_demand_map_common_plot_utils.R"),
  local = TRUE,
  chdir = TRUE
)
rm(.ivt_plot_utils_dir)

ivt_ploidy_fraction_fill_scale <- o2sd_ploidy_fraction_fill_scale

ivt_plot_lineage_counts <- function(summary_df) {
  count_df <- summary_df |>
    dplyr::distinct(dplyr::across(dplyr::any_of(c(
      "cohort", "lineage_id", "lineage_label", "scenario_id",
      "passage_index", "oxygen_pct", "predicted_live_cells", "selected_day"
    ))))
  if (!"lineage_label" %in% names(count_df)) count_df$lineage_label <- "lineage"
  if (!"cohort" %in% names(count_df)) count_df$cohort <- "cohort"

  ggplot2::ggplot(
    count_df,
    ggplot2::aes(
      passage_index,
      predicted_live_cells,
      color = factor(oxygen_pct),
      group = interaction(cohort, lineage_label)
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(cohort),
      cols = ggplot2::vars(lineage_label),
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::labs(
      title = "Predicted live cells at the selected passage time",
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
  if (!"lineage_label" %in% names(count_df)) {
    count_df$lineage_label <- if ("cohort" %in% names(count_df)) as.character(count_df$cohort) else "lineage"
  }
  if (!"lineage_passage_index" %in% names(count_df)) {
    count_df$lineage_passage_index <- count_df$passage_index
  }
  count_df$day_num <- suppressWarnings(as.numeric(count_df$day))
  count_df$live_cells_num <- suppressWarnings(as.numeric(count_df$live_cells))
  count_df$lineage_passage_index_num <- suppressWarnings(as.numeric(count_df$lineage_passage_index))
  count_df$passage_index_num <- if ("passage_index" %in% names(count_df)) {
    suppressWarnings(as.numeric(count_df$passage_index))
  } else {
    count_df$lineage_passage_index_num
  }
  count_df <- count_df[
    is.finite(count_df$day_num) &
      is.finite(count_df$live_cells_num) &
      count_df$live_cells_num > 0 &
      is.finite(count_df$lineage_passage_index_num) &
      is.finite(count_df$passage_index_num),
    ,
    drop = FALSE
  ]
  if (!nrow(count_df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Daily live-cell trajectories unavailable")
    )
  }
  cohort_values <- if ("cohort" %in% names(count_df)) unique(as.character(count_df$cohort)) else character(0)
  preferred_cohort <- c("2N", "4N")
  cohort_levels <- c(preferred_cohort[preferred_cohort %in% cohort_values], sort(setdiff(cohort_values, preferred_cohort)))
  if (length(cohort_levels)) {
    count_df$cohort <- factor(as.character(count_df$cohort), levels = cohort_levels)
  }
  lineage_values <- unique(as.character(count_df$lineage_label))
  preferred_lineage <- c("C", "O1", "O2", "control", "deprived")
  lineage_levels <- c(preferred_lineage[preferred_lineage %in% lineage_values], sort(setdiff(lineage_values, preferred_lineage)))
  count_df$lineage_label <- factor(as.character(count_df$lineage_label), levels = lineage_levels)
  oxygen_numeric <- suppressWarnings(as.numeric(count_df$oxygen_pct))
  oxygen_levels <- unique(as.character(count_df$oxygen_pct[order(oxygen_numeric, na.last = TRUE)]))
  oxygen_levels <- oxygen_levels[nzchar(oxygen_levels)]
  count_df$oxygen_factor <- factor(as.character(count_df$oxygen_pct), levels = oxygen_levels)
  oxygen_base_palette <- c(
    "0" = "#F8766D",
    "0.1" = "#C49A00",
    "0.2" = "#7CAE00",
    "0.3" = "#00BA38",
    "0.5" = "#00BFC4",
    "1" = "#00A9FF",
    "2" = "#C77CFF",
    "20.5" = "#FF61C3"
  )
  missing_oxygen <- setdiff(oxygen_levels, names(oxygen_base_palette))
  if (length(missing_oxygen)) {
    extra_cols <- grDevices::hcl.colors(length(missing_oxygen), palette = "Dark 3")
    names(extra_cols) <- missing_oxygen
    oxygen_base_palette <- c(oxygen_base_palette, extra_cols)
  }
  oxygen_palette <- oxygen_base_palette[oxygen_levels]

  group_vars <- intersect(
    c("cohort", "lineage_label", "lineage_terminal_key", "segment_id", "passage_index", "lineage_passage_index"),
    names(count_df)
  )
  if (length(group_vars)) {
    count_df$.daily_group <- do.call(
      interaction,
      c(count_df[group_vars], list(drop = TRUE, lex.order = TRUE))
    )
  } else {
    count_df$.daily_group <- seq_len(nrow(count_df))
  }
  duration_df <- stats::aggregate(
    day_num ~ .daily_group,
    data = count_df,
    FUN = function(x) max(x, na.rm = TRUE)
  )
  names(duration_df)[names(duration_df) == "day_num"] <- "duration_days"
  duration_df$duration_days[!is.finite(duration_df$duration_days) | duration_df$duration_days <= 0] <- 1
  count_df <- merge(count_df, duration_df, by = ".daily_group", all.x = TRUE, sort = FALSE)

  chosen <- count_df |>
    dplyr::filter(.data$day_num == suppressWarnings(as.numeric(.data$selected_day)))

  compact_count_label <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    out <- format(signif(x, 3), trim = TRUE, scientific = FALSE)
    big <- is.finite(x) & abs(x) >= 1e6
    mid <- is.finite(x) & abs(x) >= 1e3 & abs(x) < 1e6
    out[big] <- paste0(format(signif(x[big] / 1e6, 3), trim = TRUE), "M")
    out[mid] <- paste0(format(signif(x[mid] / 1e3, 3), trim = TRUE), "k")
    out
  }

  add_passage_label <- function(df) {
    if (is.null(df) || !nrow(df)) return(df)
    label_map <- unique(df[, c("lineage_passage_index_num", "passage_index_num"), drop = FALSE])
    label_map <- label_map[order(label_map$lineage_passage_index_num, label_map$passage_index_num), , drop = FALSE]
    label_map$key <- paste(label_map$lineage_passage_index_num, label_map$passage_index_num, sep = "__")
    label_map$label <- paste0(
      format(label_map$lineage_passage_index_num, trim = TRUE, scientific = FALSE),
      " (p",
      format(label_map$passage_index_num, trim = TRUE, scientific = FALSE),
      ")"
    )
    row_key <- paste(df$lineage_passage_index_num, df$passage_index_num, sep = "__")
    df$lineage_passage_label <- factor(
      label_map$label[match(row_key, label_map$key)],
      levels = label_map$label
    )
    df
  }

  make_nested_daily_panel <- function(cohort_value, lineage_value, show_legend = FALSE) {
    panel_df <- count_df[
      as.character(count_df$cohort) == as.character(cohort_value) &
        as.character(count_df$lineage_label) == as.character(lineage_value),
      ,
      drop = FALSE
    ]
    panel_chosen <- chosen[
      as.character(chosen$cohort) == as.character(cohort_value) &
        as.character(chosen$lineage_label) == as.character(lineage_value),
      ,
      drop = FALSE
    ]
    if (!nrow(panel_df)) {
      return(
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste(as.character(cohort_value), as.character(lineage_value), sep = " / "))
      )
    }
    panel_df <- add_passage_label(panel_df)
    panel_chosen <- add_passage_label(panel_chosen)
    legend_seed <- data.frame(
      day_num = min(panel_df$day_num, na.rm = TRUE),
      live_cells_num = min(panel_df$live_cells_num, na.rm = TRUE),
      oxygen_factor = factor(oxygen_levels, levels = oxygen_levels)
    )
    n_passages <- length(unique(panel_df$lineage_passage_label))
    facet_cols <- if (n_passages > 18L) 5L else if (n_passages > 10L) 4L else 3L
    ggplot2::ggplot(
      panel_df,
      ggplot2::aes(
        x = .data$day_num,
        y = .data$live_cells_num,
        group = .data$.daily_group,
        color = .data$oxygen_factor
      )
    ) +
      ggplot2::geom_line(linewidth = 0.62, alpha = 0.9) +
      ggplot2::geom_point(data = panel_chosen, size = 1.25, alpha = 0.9) +
      ggplot2::geom_point(
        data = legend_seed,
        ggplot2::aes(x = .data$day_num, y = .data$live_cells_num, color = .data$oxygen_factor),
        inherit.aes = FALSE,
        alpha = 0,
        size = 0.01,
        show.legend = TRUE
      ) +
      ggplot2::facet_wrap(~ lineage_passage_label, scales = "free_y", ncol = facet_cols) +
      ggplot2::scale_color_manual(values = oxygen_palette, limits = oxygen_levels, drop = FALSE) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, linewidth = 0.8, size = 2))) +
      ggplot2::scale_y_continuous(
        breaks = function(limits) pretty(limits, n = 3),
        labels = compact_count_label
      ) +
      ggplot2::labs(
        title = paste(as.character(cohort_value), as.character(lineage_value), sep = " / "),
        x = "Day within passage",
        y = "Live cells",
        color = "Oxygen (%)"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        legend.position = if (isTRUE(show_legend)) "bottom" else "none",
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        strip.text = ggplot2::element_text(face = "bold", size = 8),
        axis.title = ggplot2::element_text(size = 9),
        axis.text = ggplot2::element_text(size = 7),
        panel.spacing = grid::unit(0.45, "lines"),
        plot.margin = grid::unit(c(0.05, 0.08, 0.05, 0.08), "in")
      )
  }

  if (requireNamespace("patchwork", quietly = TRUE) && length(cohort_levels) && length(lineage_levels)) {
    panel_specs <- expand.grid(
      cohort = cohort_levels,
      lineage = lineage_levels,
      stringsAsFactors = FALSE
    )
    panel_specs <- panel_specs[order(match(panel_specs$cohort, cohort_levels), match(panel_specs$lineage, lineage_levels)), , drop = FALSE]
    panel_specs$show_legend <- seq_len(nrow(panel_specs)) == nrow(panel_specs)
    panel_plots <- Map(
      make_nested_daily_panel,
      panel_specs$cohort,
      panel_specs$lineage,
      panel_specs$show_legend
    )
    lineage_widths <- vapply(lineage_levels, function(lineage_value) {
      max(vapply(cohort_levels, function(cohort_value) {
        rows <- count_df[
          as.character(count_df$cohort) == as.character(cohort_value) &
            as.character(count_df$lineage_label) == as.character(lineage_value),
          ,
          drop = FALSE
        ]
        length(unique(rows$lineage_passage_index_num))
      }, integer(1)), 1L)
    }, numeric(1))
    return(
      patchwork::wrap_plots(
        panel_plots,
        ncol = length(lineage_levels),
        widths = pmax(lineage_widths, 1)
      ) +
        patchwork::plot_annotation(title = "Daily live-cell trajectories within each passage")
    )
  }

  count_df <- add_passage_label(count_df)
  chosen <- add_passage_label(chosen)
  legend_seed <- data.frame(
    day_num = min(count_df$day_num, na.rm = TRUE),
    live_cells_num = min(count_df$live_cells_num, na.rm = TRUE),
    oxygen_factor = factor(oxygen_levels, levels = oxygen_levels)
  )
  ggplot2::ggplot(
    count_df,
    ggplot2::aes(.data$day_num, .data$live_cells_num, group = .data$.daily_group, color = .data$oxygen_factor)
  ) +
    ggplot2::geom_line(linewidth = 0.62, alpha = 0.9) +
    ggplot2::geom_point(data = chosen, size = 1.25, alpha = 0.9) +
    ggplot2::geom_point(
      data = legend_seed,
      ggplot2::aes(x = .data$day_num, y = .data$live_cells_num, color = .data$oxygen_factor),
      inherit.aes = FALSE,
      alpha = 0,
      size = 0.01,
      show.legend = TRUE
    ) +
    ggplot2::facet_grid(cohort + lineage_label ~ lineage_passage_label, scales = "free_y") +
    ggplot2::scale_color_manual(values = oxygen_palette, limits = oxygen_levels, drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, linewidth = 0.8, size = 2))) +
    ggplot2::scale_y_continuous(breaks = function(limits) pretty(limits, n = 3), labels = compact_count_label) +
    ggplot2::labs(
      title = "Daily live-cell trajectories within each passage",
      x = "Day within passage",
      y = "Live cells",
      color = "Oxygen (%)"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
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
    ivt_ploidy_fraction_fill_scale(fraction_fill_max, name = "Fraction") +
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
