# Pure ggplot constructors for materialized in-vitro simulation/analysis tables.
# These functions do not read or write files and do not source model code.

ivt_plot_identifiability_diagnostics <- function(variance_df, loadings_df) {
  if (is.null(variance_df) || is.null(loadings_df) ||
      !all(c("direction", "variance", "weakest_rank") %in% names(variance_df)) ||
      !all(c("direction", "parameter", "loading") %in% names(loadings_df))) {
    return(NULL)
  }
  variance_df$variance <- suppressWarnings(as.numeric(variance_df$variance))
  variance_df$weakest_rank <- suppressWarnings(as.numeric(variance_df$weakest_rank))
  variance_df <- variance_df[
    is.finite(variance_df$variance) & is.finite(variance_df$weakest_rank),
    ,
    drop = FALSE
  ]
  variance_df <- variance_df[order(variance_df$weakest_rank), , drop = FALSE]
  variance_df <- utils::head(variance_df, min(5L, nrow(variance_df)))
  loadings_df$loading <- suppressWarnings(as.numeric(loadings_df$loading))
  loadings_df <- loadings_df[is.finite(loadings_df$loading), , drop = FALSE]
  if (!nrow(variance_df) || !nrow(loadings_df)) return(NULL)
  variance_df$direction <- factor(
    variance_df$direction,
    levels = variance_df$direction
  )
  loadings_df$parameter <- factor(
    loadings_df$parameter,
    levels = rev(unique(loadings_df$parameter))
  )
  p_variance <- ggplot2::ggplot(
    variance_df,
    ggplot2::aes(.data$direction, .data$variance)
  ) +
    ggplot2::geom_col(fill = "#4B6F8A", width = 0.72) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = "Smallest optimizer-population variances",
      x = NULL,
      y = "Variance (log10)"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  p_loadings <- ggplot2::ggplot(
    loadings_df,
    ggplot2::aes(.data$parameter, .data$loading, fill = .data$loading > 0)
  ) +
    ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~direction, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#2c7fb8", "FALSE" = "#d95f02")) +
    ggplot2::labs(
      title = "Dominant loadings in weakest optimizer directions",
      x = NULL,
      y = "Loading"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(
      (p_variance / p_loadings) +
        patchwork::plot_layout(heights = c(1, 2.6)) +
        patchwork::plot_annotation(title = "In-vitro identifiability diagnostics")
    )
  }
  p_loadings
}

ivt_plot_missegregation_probability_timecourse <- function(timecourse_df) {
  required <- c(
    "cohort", "lineage_label", "lineage_passage_index", "oxygen_pct", "mean_p_misseg"
  )
  if (is.null(timecourse_df) || !all(required %in% names(timecourse_df))) return(NULL)
  timecourse_df$lineage_passage_index <- suppressWarnings(as.numeric(timecourse_df$lineage_passage_index))
  timecourse_df$oxygen_pct <- suppressWarnings(as.numeric(timecourse_df$oxygen_pct))
  timecourse_df$mean_p_misseg <- suppressWarnings(as.numeric(timecourse_df$mean_p_misseg))
  timecourse_df <- timecourse_df[
    is.finite(timecourse_df$lineage_passage_index) &
      is.finite(timecourse_df$mean_p_misseg),
    ,
    drop = FALSE
  ]
  if (!nrow(timecourse_df)) return(NULL)
  ggplot2::ggplot(
    timecourse_df,
    ggplot2::aes(.data$lineage_passage_index, .data$mean_p_misseg)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = interaction(.data$cohort, .data$lineage_label)),
      color = "grey45",
      linewidth = 0.55,
      alpha = 0.8
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$oxygen_pct), size = 2.2) +
    ggplot2::facet_grid(cohort ~ lineage_label, scales = "free_x") +
    ggplot2::scale_colour_viridis_c(option = "B", na.value = "grey50") +
    ggplot2::labs(
      title = "Mean per-chromosome missegregation probability over passage",
      x = "Lineage passage index",
      y = "Mean missegregation probability",
      colour = "Fixed oxygen (%)"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

ivt_plot_functional_response_diagnostics <- function(o2_curve_df, viability_df, ploidy_o2_df) {
  if (is.null(o2_curve_df) ||
      !all(c("oxygen_pct", "cohort", "ms_rate", "proliferation_rate", "death_rate") %in% names(o2_curve_df))) {
    return(NULL)
  }
  numeric_columns <- intersect(
    c("oxygen_pct", "ms_rate", "proliferation_rate", "death_rate"),
    names(o2_curve_df)
  )
  for (name in numeric_columns) o2_curve_df[[name]] <- suppressWarnings(as.numeric(o2_curve_df[[name]]))
  make_o2_plot <- function(value, title, y_label) {
    df <- o2_curve_df[is.finite(o2_curve_df$oxygen_pct) & is.finite(o2_curve_df[[value]]), , drop = FALSE]
    ggplot2::ggplot(df, ggplot2::aes(.data$oxygen_pct, .data[[value]], color = .data$cohort)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::labs(title = title, x = "Effective oxygen (%)", y = y_label, color = "Reference state") +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(legend.position = "bottom")
  }
  plots <- list(
    make_o2_plot("ms_rate", "Oxygen vs missegregation rate", "Missegregation rate"),
    make_o2_plot("proliferation_rate", "Oxygen vs proliferation rate", "Proliferation rate"),
    make_o2_plot("death_rate", "Oxygen vs stress-associated death", "Death rate")
  )
  if (!is.null(viability_df) && all(c("endpoint_value", "viability_after_ms") %in% names(viability_df))) {
    viability_df$endpoint_value <- suppressWarnings(as.numeric(viability_df$endpoint_value))
    viability_df$viability_after_ms <- suppressWarnings(as.numeric(viability_df$viability_after_ms))
    viability_df <- viability_df[is.finite(viability_df$endpoint_value) & is.finite(viability_df$viability_after_ms), , drop = FALSE]
    if (nrow(viability_df)) {
      plots[[length(plots) + 1L]] <- ggplot2::ggplot(
        viability_df,
        ggplot2::aes(.data$endpoint_value, .data$viability_after_ms)
      ) +
        ggplot2::geom_line(color = "#2ca02c", linewidth = 0.8) +
        ggplot2::labs(title = "Post-missegregation survival", x = "State", y = "Survival") +
        ggplot2::theme_minimal(base_size = 9)
    }
  }
  if (!is.null(ploidy_o2_df) &&
      all(c("endpoint_value", "oxygen_pct", "net_growth_rate") %in% names(ploidy_o2_df))) {
    for (name in c("endpoint_value", "oxygen_pct", "net_growth_rate")) {
      ploidy_o2_df[[name]] <- suppressWarnings(as.numeric(ploidy_o2_df[[name]]))
    }
    ploidy_o2_df <- ploidy_o2_df[
      is.finite(ploidy_o2_df$endpoint_value) & is.finite(ploidy_o2_df$oxygen_pct) &
        is.finite(ploidy_o2_df$net_growth_rate),
      ,
      drop = FALSE
    ]
    if (nrow(ploidy_o2_df)) {
      plots[[length(plots) + 1L]] <- ggplot2::ggplot(
        ploidy_o2_df,
        ggplot2::aes(.data$endpoint_value, .data$net_growth_rate, color = .data$oxygen_pct, group = .data$oxygen_pct)
      ) +
        ggplot2::geom_line(linewidth = 0.65) +
        ggplot2::scale_color_viridis_c(option = "B") +
        ggplot2::labs(title = "State vs net growth by oxygen", x = "State", y = "Net growth", color = "Oxygen (%)") +
        ggplot2::theme_minimal(base_size = 9)
    }
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(
      patchwork::wrap_plots(plots, ncol = 3, guides = "collect") +
        patchwork::plot_annotation(title = "In-vitro fitted functional responses")
    )
  }
  plots[[1L]]
}

ivt_plot_missegregation_linked_relationships <- function(o2_curve_df,
                                                          o2_curve_multi_df,
                                                          viability_df) {
  plots <- list(
    ms_rate_vs_nonviable_daughter_fraction = NULL,
    death_rate_vs_missegregation_rate = NULL,
    ploidy_vs_viability_after_ms = NULL
  )

  if (!is.null(o2_curve_multi_df) && all(c(
    "cohort", "ms_rate", "misseg_nonviable_daughter_fraction"
  ) %in% names(o2_curve_multi_df))) {
    df <- o2_curve_multi_df
    df$ms_rate <- suppressWarnings(as.numeric(df$ms_rate))
    df$misseg_nonviable_daughter_fraction <- suppressWarnings(as.numeric(
      df$misseg_nonviable_daughter_fraction
    ))
    df <- df[
      is.finite(df$ms_rate) & is.finite(df$misseg_nonviable_daughter_fraction),
      ,
      drop = FALSE
    ]
    if (nrow(df)) {
      cohort_levels <- unique(as.character(df$cohort))
      df$cohort <- factor(as.character(df$cohort), levels = cohort_levels)
      cohort_colors <- stats::setNames(
        grDevices::hcl.colors(length(cohort_levels), palette = "Dark 3"),
        cohort_levels
      )
      plots$ms_rate_vs_nonviable_daughter_fraction <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          .data$ms_rate,
          .data$misseg_nonviable_daughter_fraction,
          color = .data$cohort,
          group = .data$cohort
        )
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_color_manual(values = cohort_colors, drop = FALSE) +
        ggplot2::labs(
          title = "Nonviable Daughter Fraction vs Missegregation Rate",
          subtitle = "In-vitro fitted relationship across reference chromosome-number states",
          x = "Missegregation rate",
          y = "Nonviable daughters / all daughters",
          color = "Reference state"
        ) +
        ggplot2::theme_bw(base_size = 11)
    }
  }

  if (!is.null(o2_curve_df) && all(c(
    "cohort", "death_rate", "ms_rate"
  ) %in% names(o2_curve_df))) {
    df <- o2_curve_df
    df$death_rate <- suppressWarnings(as.numeric(df$death_rate))
    df$ms_rate <- suppressWarnings(as.numeric(df$ms_rate))
    df <- df[is.finite(df$death_rate) & is.finite(df$ms_rate), , drop = FALSE]
    if (nrow(df)) {
      df$cohort <- factor(as.character(df$cohort), levels = unique(as.character(df$cohort)))
      plots$death_rate_vs_missegregation_rate <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          .data$death_rate,
          .data$ms_rate,
          color = .data$cohort,
          group = .data$cohort
        )
      ) +
        ggplot2::geom_path(linewidth = 1) +
        ggplot2::scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
        ggplot2::labs(
          title = "Stress-Associated Death Rate vs Missegregation Rate",
          subtitle = "In-vitro fitted relationship across the effective-oxygen sweep",
          x = "Stress-associated death rate",
          y = "Missegregation rate",
          color = "Reference state"
        ) +
        ggplot2::theme_bw(base_size = 11)
    }
  }

  if (!is.null(viability_df) && all(c(
    "endpoint_value", "viability_after_ms"
  ) %in% names(viability_df))) {
    df <- viability_df
    df$endpoint_value <- suppressWarnings(as.numeric(df$endpoint_value))
    df$viability_after_ms <- suppressWarnings(as.numeric(df$viability_after_ms))
    df <- df[
      is.finite(df$endpoint_value) & is.finite(df$viability_after_ms),
      ,
      drop = FALSE
    ]
    if (nrow(df)) {
      plots$ploidy_vs_viability_after_ms <- ggplot2::ggplot(
        df,
        ggplot2::aes(.data$endpoint_value, .data$viability_after_ms)
      ) +
        ggplot2::geom_line(color = "#2ca02c", linewidth = 1) +
        ggplot2::labs(
          title = "Chromosome Number vs Post-Missegregation Survival",
          subtitle = "In-vitro fitted survival for a one-copy missegregation event",
          x = "Chromosome number",
          y = "Post-missegregation survival"
        ) +
        ggplot2::theme_bw(base_size = 11)
    }
  }

  plots
}
