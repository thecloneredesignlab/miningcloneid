.joint_parameter_plot_utils_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[
    basename(frame_files) == "o2_supply_demand_map_joint_parameter_plot_utils.R"
  ]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(
  file.path(
    .joint_parameter_plot_utils_dir,
    "..",
    "..",
    "util",
    "o2_supply_demand_map_joint_parameter_utils.R"
  ),
  local = TRUE,
  chdir = TRUE
)
rm(.joint_parameter_plot_utils_dir)

joint_parameter_plot_mechanism_levels <- joint_parameter_ratio_mechanism_levels
joint_parameter_plot_mechanism_colors <- joint_parameter_ratio_mechanism_colors
joint_parameter_plot_direction_colors <- joint_parameter_ratio_direction_colors

joint_parameter_plot_value_label <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return("")
  if (abs(x) >= 0.01 && abs(x) < 1000) return(sprintf("%.2fx", x))
  paste0(formatC(x, format = "e", digits = 2), "x")
}

joint_parameter_plot_log2_axis_label <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return("")
  if (abs(x - round(x)) < 1e-8) return(as.character(as.integer(round(x))))
  sprintf("%.2f", x)
}

joint_parameter_plot_log2_natural_label <- function(log2_value, ratio_value) {
  log2_value <- suppressWarnings(as.numeric(log2_value))
  if (!is.finite(log2_value)) return("")
  paste0(
    sprintf("%.2f", log2_value),
    "(",
    joint_parameter_plot_value_label(ratio_value),
    ")"
  )
}

joint_parameter_ratio_plot <- function(plot_df, fit_label = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for joint parameter visualization.")
  }
  if (!is.data.frame(plot_df) || !nrow(plot_df)) {
    stop("joint parameter analysis table has no rows.", call. = FALSE)
  }
  required <- c("parameter", "ratio_vivo_to_vitro", "mechanism_class")
  missing <- setdiff(required, names(plot_df))
  if (length(missing)) {
    stop(
      "joint parameter analysis table is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  plot_df$ratio_vivo_to_vitro <- suppressWarnings(
    as.numeric(plot_df$ratio_vivo_to_vitro)
  )
  plot_df <- plot_df[
    is.finite(plot_df$ratio_vivo_to_vitro) &
      plot_df$ratio_vivo_to_vitro > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) {
    stop("joint parameter analysis table has no finite positive ratios.")
  }
  plot_df <- plot_df[order(plot_df$ratio_vivo_to_vitro), , drop = FALSE]
  plot_df$log2_fold_change <- log2(plot_df$ratio_vivo_to_vitro)
  plot_df$mechanism_class <- factor(
    plot_df$mechanism_class,
    levels = joint_parameter_plot_mechanism_levels()
  )

  n <- nrow(plot_df)
  plot_df$y <- rev(seq_len(n))
  plot_df$segment_x <- pmin(0, plot_df$log2_fold_change)
  plot_df$segment_xend <- pmax(0, plot_df$log2_fold_change)
  plot_df$ratio_label <- mapply(
    joint_parameter_plot_log2_natural_label,
    plot_df$log2_fold_change,
    plot_df$ratio_vivo_to_vitro,
    USE.NAMES = FALSE
  )
  plot_df$direction <- ifelse(
    plot_df$ratio_vivo_to_vitro < 1 / 1.10,
    "Higher in vitro",
    ifelse(
      plot_df$ratio_vivo_to_vitro > 1.10,
      "Higher in vivo",
      "Near 1x"
    )
  )
  plot_df$direction <- factor(
    plot_df$direction,
    levels = c("Higher in vitro", "Near 1x", "Higher in vivo")
  )
  plot_df$ratio_label_x <- ifelse(
    plot_df$log2_fold_change < 0,
    plot_df$log2_fold_change - 0.08,
    plot_df$log2_fold_change + 0.08
  )
  plot_df$ratio_hjust <- ifelse(plot_df$log2_fold_change < 0, 1, 0)

  limit_power <- max(
    6,
    ceiling(max(abs(plot_df$log2_fold_change), na.rm = TRUE))
  )
  limit_power <- min(max(limit_power, 1), 20)
  x_breaks <- seq(-limit_power, limit_power, by = 1)
  x_limits <- c(-limit_power - 0.75, limit_power + 0.75)
  parameter_label_x <- -limit_power + 0.10
  subtitle <- paste(
    "Point/segment color = direction of difference",
    "Label color = mechanism class",
    sep = " | "
  )
  if (!is.null(fit_label) &&
      length(fit_label) &&
      nzchar(as.character(fit_label[[1]]))) {
    subtitle <- paste0(as.character(fit_label[[1]]), " | ", subtitle)
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(y = y)) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "black",
      linetype = "22",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = segment_x,
        xend = segment_xend,
        yend = y,
        color = direction
      ),
      linewidth = 1.35,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = log2_fold_change, color = direction),
      size = 4.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = ratio_label_x,
        label = ratio_label,
        color = direction,
        hjust = ratio_hjust
      ),
      size = 4.6,
      show.legend = FALSE
    ) +
    ggplot2::geom_label(
      ggplot2::aes(
        x = parameter_label_x,
        label = parameter,
        fill = mechanism_class
      ),
      color = "black",
      fontface = "bold",
      size = 4.9,
      linewidth = 0.25,
      label.padding = grid::unit(0.16, "lines"),
      label.r = grid::unit(0.11, "lines"),
      show.legend = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = vapply(
        x_breaks,
        joint_parameter_plot_log2_axis_label,
        character(1)
      ),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.55, n + 0.75),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_color_manual(
      name = "Direction",
      values = joint_parameter_plot_direction_colors(),
      breaks = c("Higher in vitro", "Higher in vivo"),
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = joint_parameter_plot_mechanism_colors(),
      breaks = joint_parameter_plot_mechanism_levels(),
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(xlim = x_limits, clip = "off") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(label = "a", color = "black", size = 4.4)
      ),
      color = ggplot2::guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(size = 4.5, linewidth = 1.5)
      )
    ) +
    ggplot2::labs(
      title = "In Vivo vs In Vitro Parameter Ratios",
      subtitle = subtitle,
      x = "log2 fold change (vivo / vitro)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(
        color = "#E6E6E6",
        linewidth = 0.7
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "#555555", size = 12),
      axis.title.x = ggplot2::element_text(size = 15),
      plot.title = ggplot2::element_text(
        size = 23,
        face = "plain",
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        size = 14,
        hjust = 0,
        color = "black"
      ),
      legend.position = "top",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.spacing.x = grid::unit(0.30, "cm"),
      legend.margin = ggplot2::margin(0, 0, 6, 0),
      plot.margin = ggplot2::margin(10, 22, 10, 10)
    )
}
