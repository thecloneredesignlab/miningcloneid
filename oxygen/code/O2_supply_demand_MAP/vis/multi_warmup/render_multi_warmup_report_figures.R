#!/usr/bin/env Rscript

# Pure visualization for the multi-warmup HTML report. All input tables are
# materialized by analysis stages before this module runs.

.o2mw_report_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_multi_warmup_report_figures.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_report_vis_workflow_root <- normalizePath(file.path(.o2mw_report_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_report_vis_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

plot_integrated_tradeoff_by_invivo_warmup <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  summary_path <- file.path(root, "multi_warmup_integrated_joint_run", "extra_results", "seed_summary.tsv")
  summary_df <- read_table_optional(summary_path)
  required <- c("seed", "objective_invivo", "objective_invitro")
  if (is.null(summary_df) || !is.data.frame(summary_df) || !all(required %in% names(summary_df))) {
    return(invisible(NULL))
  }
  plot_df <- summary_df[, required, drop = FALSE]
  plot_df$objective_invivo <- num_col(plot_df$objective_invivo)
  plot_df$objective_invitro <- num_col(plot_df$objective_invitro)
  plot_df$warmup_label <- sub("__.*$", "", as.character(plot_df$seed))
  plot_df$invivo_warmup <- toupper(sub("_.*$", "", plot_df$warmup_label))
  plot_df <- plot_df[
    is.finite(plot_df$objective_invivo) &
      is.finite(plot_df$objective_invitro) &
      nzchar(plot_df$invivo_warmup),
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))

  manifest <- read_table_optional(file.path(root, "multi_warmup_manifest.tsv"))
  if (is.data.frame(manifest) && nrow(manifest) && "warmup_label" %in% names(manifest)) {
    warmup_levels <- unique(toupper(sub("_.*$", "", as.character(manifest$warmup_label))))
    warmup_levels <- warmup_levels[nzchar(warmup_levels)]
  } else {
    warmup_levels <- unique(plot_df$invivo_warmup)
  }
  plot_df$invivo_warmup <- factor(plot_df$invivo_warmup, levels = warmup_levels)

  plot_data_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.tsv")
  utils::write.table(plot_df, plot_data_path, sep = "\t", quote = FALSE, row.names = FALSE)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = objective_invivo, y = objective_invitro, color = invivo_warmup)
  ) +
    ggplot2::geom_point(size = 1, alpha = 0.9) +
    ggplot2::labs(
      title = "Joint Objective Tradeoff by In Vivo Warm-Up Family",
      subtitle = "Each point is one integrated joint seed; color shows Vxx warm-up source.",
      x = "In vivo objective",
      y = "In vitro objective",
      color = "In vivo warm-up"
    ) +
    scale_color_invivo_warmup(warmup_levels) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  pdf_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.pdf")
  png_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.png")
  ggplot2::ggsave(pdf_path, p, width = 9, height = 7, bg = "white")
  ggplot2::ggsave(png_path, p, width = 9, height = 7, dpi = 220, bg = "white")
  invisible(pdf_path)
}

invivo_warmup_palette <- family_palette

scale_color_invivo_warmup <- function(levels, name = "In vivo warm-up", breaks = levels, drop = FALSE) {
  ggplot2::scale_color_manual(
    values = invivo_warmup_palette(levels),
    limits = levels,
    breaks = breaks,
    drop = drop,
    name = name
  )
}

save_plot_pdf_png <- function(plot, root, stem, width, height, dpi = 220) {
  pdf_path <- file.path(root, paste0(stem, ".pdf"))
  png_path <- file.path(root, paste0(stem, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(c(pdf = pdf_path, png = png_path))
}

point_shape_values <- function(n) {
  base <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 5, 6, 9, 10, 12, 13, 14)
  rep(base, length.out = n)
}

square_limits <- function(x, y, pad_frac = 0.08) {
  x <- num_col(x)
  y <- num_col(y)
  xr <- range(x[is.finite(x)], na.rm = TRUE)
  yr <- range(y[is.finite(y)], na.rm = TRUE)
  cx <- mean(xr)
  cy <- mean(yr)
  span <- max(diff(xr), diff(yr))
  if (!is.finite(span) || span <= 0) span <- 1
  half <- span * (0.5 + pad_frac)
  list(x = c(cx - half, cx + half), y = c(cy - half, cy + half))
}

wrap_label <- function(x, width = 68L) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

plot_part1_objective_figures <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  summary_df <- read_table_optional(file.path(root, "multi_warmup_best_seed_summary.tsv"))
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df)) return(invisible(NULL))
  if (!all(c("warmup_label", "objective", "invivo_family") %in% names(summary_df))) return(invisible(NULL))
  summary_df$objective <- num_col(summary_df$objective)
  family_levels <- unique(as.character(summary_df$invivo_family))
  family_levels <- family_levels[nzchar(family_levels)]
  summary_df$invivo_family <- factor(as.character(summary_df$invivo_family), levels = family_levels)

  p1 <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = stats::reorder(warmup_label, objective), y = objective, color = invivo_family)
  ) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::coord_flip() +
    scale_color_invivo_warmup(family_levels, name = "In vivo family") +
    ggplot2::labs(
      title = "Multi-Warm-Up Best Joint Objective",
      x = "Warm-up pair",
      y = "Best total objective"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  save_plot_pdf_png(p1, root, "multi_warmup_objective_summary", width = 8, height = 5)

  if (all(c("objective_invivo", "objective_invitro") %in% names(summary_df))) {
    summary_df$objective_invivo <- num_col(summary_df$objective_invivo)
    summary_df$objective_invitro <- num_col(summary_df$objective_invitro)
    plot_df <- summary_df[is.finite(summary_df$objective_invivo) & is.finite(summary_df$objective_invitro), , drop = FALSE]
    if (nrow(plot_df)) {
      p2 <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(objective_invivo, objective_invitro, color = invivo_family, label = warmup_label)
      ) +
        ggplot2::geom_point(size = 2.8) +
        scale_color_invivo_warmup(family_levels, name = "In vivo family") +
        ggplot2::labs(
          title = "In Vivo vs In Vitro Objective by Warm-Up Pair",
          x = "Best seed in vivo objective",
          y = "Best seed in vitro objective"
        ) +
        ggplot2::theme_minimal(base_size = 11)
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p2 <- p2 + ggrepel::geom_text_repel(size = 3, max.overlaps = Inf)
      } else {
        p2 <- p2 + ggplot2::geom_text(size = 3, vjust = -0.6)
      }
      save_plot_pdf_png(p2, root, "multi_warmup_invivo_invitro_objective_scatter", width = 7, height = 5)
    }
  }
  invisible(NULL)
}

plot_part1_cluster_umap <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  coords <- read_table_optional(file.path(root, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_coords.tsv"))
  required <- c("pair_id", "invivo_cluster", "invitro_rank", "UMAP1", "UMAP2")
  if (is.null(coords) || !is.data.frame(coords) || !all(required %in% names(coords))) return(invisible(NULL))
  plot_df <- coords
  cluster_levels <- unique(as.character(plot_df$invivo_cluster))
  cluster_levels <- cluster_levels[nzchar(cluster_levels)]
  plot_df$invivo_cluster <- factor(as.character(plot_df$invivo_cluster), levels = cluster_levels)
  top_n <- max(suppressWarnings(as.integer(plot_df$invitro_rank)), na.rm = TRUE)
  if (!is.finite(top_n) || top_n <= 0L) top_n <- length(unique(plot_df$invitro_rank))
  plot_df$invitro_rank_factor <- factor(plot_df$invitro_rank, levels = seq_len(top_n))
  if (!("is_invivo_representative" %in% names(plot_df))) plot_df$is_invivo_representative <- FALSE
  plot_df$is_invivo_representative <- tolower(as.character(plot_df$is_invivo_representative)) %in% c("true", "t", "1", "yes")
  limits <- square_limits(plot_df$UMAP1, plot_df$UMAP2)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(UMAP1, UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = invivo_cluster, shape = invitro_rank_factor), size = 3.0, stroke = 0.9) +
    scale_color_invivo_warmup(cluster_levels, name = "In Vivo Cluster") +
    ggplot2::scale_shape_manual(values = point_shape_values(top_n), name = "In Vitro Rank") +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3.2)),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::labs(
      title = "18D Warm-Start Profile UMAP by In Vivo Cluster",
      subtitle = wrap_label(
        "Same 18D paired warm-start profile UMAP as the standalone view, colored by profile-space in vivo cluster. Red labels mark selected representative in vivo ranks.",
        width = 68L
      ),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    ggplot2::coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE, clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical"
    )
  label_nonrep <- plot_df[!plot_df$is_invivo_representative, , drop = FALSE]
  label_rep <- plot_df[plot_df$is_invivo_representative, , drop = FALSE]
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    if (nrow(label_nonrep)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_nonrep,
        ggplot2::aes(label = pair_id),
        size = 2.1,
        color = "grey15",
        max.overlaps = Inf,
        seed = 1L
      )
    }
    if (nrow(label_rep)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_rep,
        ggplot2::aes(label = pair_id),
        size = 2.35,
        color = "#d62728",
        fontface = "bold",
        max.overlaps = Inf,
        seed = 1L
      )
    }
  } else {
    if (nrow(label_nonrep)) {
      p <- p + ggplot2::geom_text(data = label_nonrep, ggplot2::aes(label = pair_id), size = 2.1, color = "grey15", vjust = -0.6)
    }
    if (nrow(label_rep)) {
      p <- p + ggplot2::geom_text(data = label_rep, ggplot2::aes(label = pair_id), size = 2.35, color = "#d62728", fontface = "bold", vjust = -0.6)
    }
  }
  save_plot_pdf_png(p, root, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_500seed", width = 7, height = 7)
  invisible(NULL)
}

plot_part1_color_aligned_figures <- function(root) {
  plot_part1_cluster_umap(root)
  plot_part1_objective_figures(root)
  invisible(NULL)
}

mw_boundary_axis_config <- function(x,
                                    near_thresh,
                                    x_scale = c("relative", "log10_original"),
                                    raw_value = NULL,
                                    raw_lower = NULL,
                                    raw_upper = NULL) {
  x_scale <- match.arg(x_scale)
  if (identical(x_scale, "relative")) {
    return(list(
      axis_type = "relative",
      x = x,
      lower_limit = 0,
      upper_limit = 1,
      vlines = c(0, near_thresh, 0.5, 1 - near_thresh, 1),
      vline_colors = c("grey50", "grey70", "grey82", "grey70", "grey50"),
      vline_linetypes = c("solid", "dashed", "dotted", "dashed", "solid"),
      lower_rect = c(0, near_thresh),
      upper_rect = c(1 - near_thresh, 1),
      scale = ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, near_thresh, 0.5, 1 - near_thresh, 1)),
      x_label = "Relative position in transformed fit range",
      subtitle_note = ""
    ))
  }

  raw_value <- suppressWarnings(as.numeric(raw_value))
  raw_lower <- suppressWarnings(as.numeric(raw_lower))
  raw_upper <- suppressWarnings(as.numeric(raw_upper))
  positive_x <- c(raw_value, raw_lower, raw_upper)
  positive_x <- positive_x[is.finite(positive_x) & positive_x > 0]
  log_floor <- min(positive_x, na.rm = TRUE) / 10
  if (!is.finite(log_floor) || log_floor <= 0) log_floor <- 1e-6
  lower_plot <- ifelse(is.finite(raw_lower), pmax(raw_lower, log_floor), NA_real_)
  upper_plot <- ifelse(is.finite(raw_upper), pmax(raw_upper, log_floor), NA_real_)
  value_plot <- ifelse(is.finite(raw_value), pmax(raw_value, log_floor), NA_real_)
  axis_upper <- max(c(upper_plot, value_plot), na.rm = TRUE)
  if (!is.finite(axis_upper) || axis_upper <= log_floor) axis_upper <- log_floor * 10
  breaks <- 10^seq(floor(log10(log_floor)), ceiling(log10(axis_upper)), by = 1)
  breaks <- breaks[breaks >= log_floor & breaks <= axis_upper]
  needs_floor_label <- any(
    is.finite(c(raw_value, raw_lower, raw_upper)) & c(raw_value, raw_lower, raw_upper) <= 0,
    na.rm = TRUE
  )
  labels <- function(vals) {
    vapply(
      vals,
      function(v) {
        if (isTRUE(needs_floor_label) && isTRUE(all.equal(v, log_floor, tolerance = 1e-12))) return("0")
        formatC(v, format = "fg", digits = 4)
      },
      character(1)
    )
  }
  list(
    axis_type = "log10_original",
    x = value_plot,
    lower_limit = log_floor,
    upper_limit = axis_upper,
    lower_plot = lower_plot,
    upper_plot = upper_plot,
    scale = ggplot2::scale_x_log10(limits = c(log_floor, axis_upper), breaks = breaks, labels = labels),
    x_label = "Original parameter value (log10 scale)",
    subtitle_note = if (isTRUE(needs_floor_label)) " Values <=0 shown at log floor." else ""
  )
}

mw_top_ranked_seeds <- function(summary_df, n = 3L, pred1000_gate = FALSE) {
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df) || !all(c("seed", "objective") %in% names(summary_df))) {
    return(character(0))
  }
  keep <- is.finite(num_col(summary_df$objective))
  if (isTRUE(pred1000_gate) && "pred1000_both_gt44" %in% names(summary_df)) {
    gate <- summary_df$pred1000_both_gt44
    if (!is.logical(gate)) gate <- tolower(as.character(gate)) %in% c("true", "t", "1", "yes")
    keep <- keep & isTRUE(TRUE) & gate
  }
  ranked <- summary_df[keep, , drop = FALSE]
  if (!nrow(ranked)) return(character(0))
  ranked$objective_num <- num_col(ranked$objective)
  ranked <- ranked[order(ranked$objective_num), , drop = FALSE]
  head(as.character(ranked$seed), n)
}

plot_multi_warmup_boundary_forest <- function(root,
                                              stem,
                                              title_suffix = NULL,
                                              pred1000_top3_only = FALSE,
                                              x_scale = c("relative", "log10_original"),
                                              near_thresh = 0.05) {
  x_scale <- match.arg(x_scale)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  extra_dir <- integrated_extra_results_dir(root)
  long_df <- read_table_optional(file.path(extra_dir, "parameter_boundary_long.tsv"))
  summary_df <- read_table_optional(file.path(extra_dir, "seed_summary.tsv"))
  required <- c(
    "seed",
    "active_in_fit",
    "rel_pos_plot",
    "rel_dist_to_nearest",
    "param_prototype",
    "prototype_value",
    "prototype_lower_bound",
    "prototype_upper_bound"
  )
  if (is.null(long_df) || !is.data.frame(long_df) || !all(required %in% names(long_df))) return(invisible(NULL))
  plot_df <- long_df[as.logical(long_df$active_in_fit) & is.finite(num_col(long_df$rel_pos_plot)), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  top3_seeds <- mw_top_ranked_seeds(summary_df, n = 3L, pred1000_gate = pred1000_top3_only)
  plot_df <- add_invivo_warmup_from_seed(plot_df, root)
  plot_df$rel_pos_plot <- num_col(plot_df$rel_pos_plot)
  plot_df$rel_dist_to_nearest <- num_col(plot_df$rel_dist_to_nearest)
  plot_df$prototype_value <- num_col(plot_df$prototype_value)
  plot_df$prototype_lower_bound <- num_col(plot_df$prototype_lower_bound)
  plot_df$prototype_upper_bound <- num_col(plot_df$prototype_upper_bound)
  axis_cfg <- mw_boundary_axis_config(
    plot_df$rel_pos_plot,
    near_thresh = near_thresh,
    x_scale = x_scale,
    raw_value = plot_df$prototype_value,
    raw_lower = plot_df$prototype_lower_bound,
    raw_upper = plot_df$prototype_upper_bound
  )
  plot_df$boundary_x_plot <- axis_cfg$x
  if (identical(axis_cfg$axis_type, "log10_original")) {
    plot_df$boundary_x_lower <- axis_cfg$lower_plot
    plot_df$boundary_x_upper <- axis_cfg$upper_plot
    plot_df <- plot_df[
      is.finite(plot_df$boundary_x_plot) &
        is.finite(plot_df$boundary_x_lower) &
        is.finite(plot_df$boundary_x_upper) &
        plot_df$boundary_x_upper > plot_df$boundary_x_lower,
      ,
      drop = FALSE
    ]
    if (!nrow(plot_df)) return(invisible(NULL))
  }
  param_rank <- tapply(plot_df$rel_dist_to_nearest, plot_df$param_prototype, min, na.rm = TRUE)
  param_levels <- names(sort(param_rank, decreasing = FALSE))
  plot_df$param_prototype <- factor(plot_df$param_prototype, levels = rev(param_levels))
  point_pos <- ggplot2::position_jitter(height = 0.14, width = 0, seed = 1L)
  title_detail <- if (!is.null(title_suffix) && nzchar(title_suffix)) paste0(" (", title_suffix, ")") else ""
  top3_label <- if (length(top3_seeds)) {
    if (isTRUE(pred1000_top3_only)) {
      " Black symbols mark eligible objective top 3."
    } else {
      " Black symbols mark objective top 3."
    }
  } else if (isTRUE(pred1000_top3_only)) {
    "No seeds met the 2N/4N 1000-day prediction gate."
  } else {
    ""
  }
  color_levels <- levels(plot_df$invivo_warmup)
  color_breaks <- color_levels
  point_size <- 2.0
  point_alpha <- 0.62
  top_df <- plot_df[as.character(plot_df$seed) %in% top3_seeds, , drop = FALSE]
  if (nrow(top_df)) {
    marker_levels <- c(
      if (length(top3_seeds) >= 1L) paste0("Top 1: ", top3_seeds[[1]], " (*)"),
      if (length(top3_seeds) >= 2L) paste0("Top 2: ", top3_seeds[[2]], " (triangle)"),
      if (length(top3_seeds) >= 3L) paste0("Top 3: ", top3_seeds[[3]], " (black dot)")
    )
    names(marker_levels) <- top3_seeds
    top_df$top3_marker <- factor(marker_levels[as.character(top_df$seed)], levels = unname(marker_levels))
    top_shape_values <- c(8, 17, 16)[seq_along(marker_levels)]
    names(top_shape_values) <- unname(marker_levels)
  } else {
    top_shape_values <- numeric(0)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = boundary_x_plot, y = param_prototype))
  if (identical(axis_cfg$axis_type, "relative")) {
    ref_df <- unique(plot_df["param_prototype"])
    ref_df$boundary_x_start <- axis_cfg$lower_limit
    ref_df$boundary_x_end <- axis_cfg$upper_limit
    p <- p +
      ggplot2::annotate("rect", xmin = axis_cfg$lower_rect[[1]], xmax = axis_cfg$lower_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#fddbc7", alpha = 0.28) +
      ggplot2::annotate("rect", xmin = axis_cfg$upper_rect[[1]], xmax = axis_cfg$upper_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#d1e5f0", alpha = 0.28) +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = boundary_x_start, xend = boundary_x_end, y = param_prototype, yend = param_prototype),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      ) +
      ggplot2::geom_vline(xintercept = axis_cfg$vlines, color = axis_cfg$vline_colors, linetype = axis_cfg$vline_linetypes, linewidth = 0.35)
  } else {
    ref_df <- unique(plot_df[c("param_prototype", "boundary_x_lower", "boundary_x_upper")])
    p <- p +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = boundary_x_lower, xend = boundary_x_upper, y = param_prototype, yend = param_prototype),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      )
  }

  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(color = invivo_warmup),
      shape = 16,
      size = point_size,
      alpha = point_alpha,
      position = point_pos
    ) +
    axis_cfg$scale +
    ggplot2::labs(
      title = paste0("Parameter Positions Within Fitted Bounds", title_detail, ": integrated multi-warm-up"),
      subtitle = paste0(
        if (identical(axis_cfg$axis_type, "relative")) {
          paste0("0=lower, 1=upper; shaded zones are within ", sprintf("%.0f", 100 * near_thresh), "% of a bound. ")
        } else {
          "Horizontal lines span original bounds. "
        },
        "Points are colored by Vxx warm-up family.",
        top3_label,
        axis_cfg$subtitle_note
      ),
      x = axis_cfg$x_label,
      y = NULL,
      color = "In vivo warm-up"
    ) +
    scale_color_invivo_warmup(color_levels, breaks = color_breaks) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )
  if (nrow(top_df)) {
    p <- p +
      ggplot2::geom_point(
        data = top_df,
        ggplot2::aes(shape = top3_marker),
        color = "black",
        size = 3.0,
        position = point_pos
      ) +
      ggplot2::scale_shape_manual(values = top_shape_values, drop = FALSE) +
      ggplot2::labs(shape = "Top seed marker") +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1), order = 1),
        shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 2)
      )
  }
  save_plot_pdf_png(p, root, stem, width = 12, height = 5.2)
}

plot_multi_warmup_paired_boundary_scatter <- function(root,
                                                      stem,
                                                      title_suffix = NULL,
                                                      pred1000_top3_only = FALSE,
                                                      x_scale = c("relative", "log10_original"),
                                                      near_thresh = 0.05) {
  x_scale <- match.arg(x_scale)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  extra_dir <- integrated_extra_results_dir(root)
  soft_df <- read_table_optional(file.path(extra_dir, "joint_soft_coupling_all.tsv"))
  summary_df <- read_table_optional(file.path(extra_dir, "seed_summary.tsv"))
  required <- c(
    "seed",
    "parameter",
    "vivo_transformed",
    "vitro_transformed",
    "vivo_natural",
    "vitro_natural"
  )
  if (is.null(soft_df) || !is.data.frame(soft_df) || !nrow(soft_df) || !all(required %in% names(soft_df))) {
    return(invisible(NULL))
  }
  lower_t_col <- if ("joint_union_lower_transformed" %in% names(soft_df)) "joint_union_lower_transformed" else "center_lower_transformed"
  upper_t_col <- if ("joint_union_upper_transformed" %in% names(soft_df)) "joint_union_upper_transformed" else "center_upper_transformed"
  lower_nat_col <- if ("joint_union_lower_bound" %in% names(soft_df)) "joint_union_lower_bound" else "center_lower_bound"
  upper_nat_col <- if ("joint_union_upper_bound" %in% names(soft_df)) "joint_union_upper_bound" else "center_upper_bound"
  bound_cols <- c(lower_t_col, upper_t_col, lower_nat_col, upper_nat_col)
  if (!all(bound_cols %in% names(soft_df))) return(invisible(NULL))

  make_context_df <- function(context, value_t_col, value_nat_col) {
    data.frame(
      seed = as.character(soft_df$seed),
      parameter = as.character(soft_df$parameter),
      context = context,
      value_transformed = num_col(soft_df[[value_t_col]]),
      value_natural = num_col(soft_df[[value_nat_col]]),
      bound_lower_transformed = num_col(soft_df[[lower_t_col]]),
      bound_upper_transformed = num_col(soft_df[[upper_t_col]]),
      bound_lower_natural = num_col(soft_df[[lower_nat_col]]),
      bound_upper_natural = num_col(soft_df[[upper_nat_col]]),
      stringsAsFactors = FALSE
    )
  }
  plot_df <- rbind(
    make_context_df("in vivo", "vivo_transformed", "vivo_natural"),
    make_context_df("in vitro", "vitro_transformed", "vitro_natural")
  )
  width_t <- plot_df$bound_upper_transformed - plot_df$bound_lower_transformed
  plot_df$rel_pos_in_range <- (plot_df$value_transformed - plot_df$bound_lower_transformed) / width_t
  plot_df$rel_pos_plot <- pmin(pmax(plot_df$rel_pos_in_range, 0), 1)
  plot_df$rel_dist_to_lower <- (plot_df$value_transformed - plot_df$bound_lower_transformed) / width_t
  plot_df$rel_dist_to_upper <- (plot_df$bound_upper_transformed - plot_df$value_transformed) / width_t
  plot_df$rel_dist_to_nearest <- pmin(plot_df$rel_dist_to_lower, plot_df$rel_dist_to_upper)
  plot_df <- plot_df[
    is.finite(plot_df$rel_pos_plot) &
      is.finite(plot_df$rel_dist_to_nearest) &
      is.finite(width_t) &
      width_t > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df <- add_invivo_warmup_from_seed(plot_df, root)
  family_levels <- levels(plot_df$invivo_warmup)
  top3_seeds <- mw_top_ranked_seeds(summary_df, n = 3L, pred1000_gate = pred1000_top3_only)

  axis_cfg <- mw_boundary_axis_config(
    plot_df$rel_pos_plot,
    near_thresh = near_thresh,
    x_scale = x_scale,
    raw_value = plot_df$value_natural,
    raw_lower = plot_df$bound_lower_natural,
    raw_upper = plot_df$bound_upper_natural
  )
  plot_df$boundary_x_plot <- axis_cfg$x
  if (identical(axis_cfg$axis_type, "log10_original")) {
    plot_df$boundary_x_lower <- axis_cfg$lower_plot
    plot_df$boundary_x_upper <- axis_cfg$upper_plot
    plot_df <- plot_df[
      is.finite(plot_df$boundary_x_plot) &
        is.finite(plot_df$boundary_x_lower) &
        is.finite(plot_df$boundary_x_upper) &
        plot_df$boundary_x_upper > plot_df$boundary_x_lower,
      ,
      drop = FALSE
    ]
    if (!nrow(plot_df)) return(invisible(NULL))
  }

  plot_df$parameter_label <- as.character(plot_df$parameter)
  param_rank <- tapply(plot_df$rel_dist_to_nearest, plot_df$parameter_label, min, na.rm = TRUE)
  param_levels <- names(sort(param_rank, decreasing = FALSE))
  plot_df$parameter <- factor(plot_df$parameter_label, levels = rev(param_levels))
  y_breaks <- seq_along(levels(plot_df$parameter))
  plot_df$y_base <- as.numeric(plot_df$parameter)
  plot_df$context <- factor(plot_df$context, levels = c("in vivo", "in vitro"))
  pair_key <- paste(plot_df$seed, plot_df$parameter_label, sep = "\r")
  pair_hash <- vapply(
    pair_key,
    function(key) {
      ints <- utf8ToInt(key)
      if (!length(ints)) return(0)
      sum((seq_along(ints) * ints) %% 997)
    },
    numeric(1)
  )
  plot_df$pair_jitter <- ((pair_hash %% 101) / 100 - 0.5) * 0.08
  plot_df$context_offset <- ifelse(as.character(plot_df$context) == "in vivo", 0.18, -0.18)
  plot_df$y_plot <- plot_df$y_base + plot_df$context_offset + plot_df$pair_jitter
  ref_df <- if (identical(axis_cfg$axis_type, "log10_original")) {
    unique(plot_df[c("parameter", "y_base", "boundary_x_lower", "boundary_x_upper")])
  } else {
    ref <- unique(plot_df[c("parameter", "y_base")])
    ref$boundary_x_start <- axis_cfg$lower_limit
    ref$boundary_x_end <- axis_cfg$upper_limit
    ref
  }

  plot_df$seed_marker <- "Other seeds"
  if (length(top3_seeds) >= 1L) plot_df$seed_marker[plot_df$seed == top3_seeds[[1]]] <- paste0("Top 1: ", top3_seeds[[1]], " (*)")
  if (length(top3_seeds) >= 2L) plot_df$seed_marker[plot_df$seed == top3_seeds[[2]]] <- paste0("Top 2: ", top3_seeds[[2]], " (triangle)")
  if (length(top3_seeds) >= 3L) plot_df$seed_marker[plot_df$seed == top3_seeds[[3]]] <- paste0("Top 3: ", top3_seeds[[3]], " (black dot)")
  other_df <- plot_df[plot_df$seed_marker == "Other seeds", , drop = FALSE]
  top_df <- plot_df[plot_df$seed_marker != "Other seeds", , drop = FALSE]
  other_vivo_df <- other_df[as.character(other_df$context) == "in vivo", , drop = FALSE]
  other_vitro_df <- other_df[as.character(other_df$context) == "in vitro", , drop = FALSE]
  top_breaks <- c(
    if (length(top3_seeds) >= 1L) paste0("Top 1: ", top3_seeds[[1]], " (*)"),
    if (length(top3_seeds) >= 2L) paste0("Top 2: ", top3_seeds[[2]], " (triangle)"),
    if (length(top3_seeds) >= 3L) paste0("Top 3: ", top3_seeds[[3]], " (black dot)")
  )
  shape_values <- c(
    if (length(top3_seeds) >= 1L) setNames(8, paste0("Top 1: ", top3_seeds[[1]], " (*)")),
    if (length(top3_seeds) >= 2L) setNames(17, paste0("Top 2: ", top3_seeds[[2]], " (triangle)")),
    if (length(top3_seeds) >= 3L) setNames(16, paste0("Top 3: ", top3_seeds[[3]], " (black dot)"))
  )
  top3_label <- if (length(top3_seeds)) {
    if (isTRUE(pred1000_top3_only)) {
      paste0("Black symbols mark eligible objective top 3: ", paste(top3_seeds, collapse = "; "), ".")
    } else {
      paste0("Black symbols mark objective top 3: ", paste(top3_seeds, collapse = "; "), ".")
    }
  } else if (isTRUE(pred1000_top3_only)) {
    "No seeds met the 2N/4N 1000-day prediction gate."
  } else {
    ""
  }

  vivo_pair_df <- plot_df[
    as.character(plot_df$context) == "in vivo",
    c("seed", "parameter_label", "boundary_x_plot", "y_plot", "seed_marker"),
    drop = FALSE
  ]
  vitro_pair_df <- plot_df[
    as.character(plot_df$context) == "in vitro",
    c("seed", "parameter_label", "boundary_x_plot", "y_plot", "seed_marker"),
    drop = FALSE
  ]
  pair_df <- merge(
    vivo_pair_df,
    vitro_pair_df,
    by = c("seed", "parameter_label"),
    suffixes = c("_vivo", "_vitro"),
    all = FALSE,
    sort = FALSE
  )
  pair_df$is_top_seed <- pair_df$seed %in% top3_seeds
  other_pair_df <- pair_df[!pair_df$is_top_seed, , drop = FALSE]
  top_pair_df <- pair_df[pair_df$is_top_seed, , drop = FALSE]

  title_detail <- if (!is.null(title_suffix) && nzchar(title_suffix)) paste0(" (", title_suffix, ")") else ""
  subtitle_line1 <- if (identical(axis_cfg$axis_type, "relative")) {
    paste0(
      "0 = joint union lower bound, 1 = joint union upper bound; shaded zones are within ",
      sprintf("%.0f", 100 * near_thresh),
      "% of a bound."
    )
  } else {
    "Horizontal lines span natural joint union lower-to-upper bounds."
  }
  subtitle_text <- paste(
    c(
      subtitle_line1,
      "In vivo points are colored by Vxx warm-up family; in vitro points are blue; lines connect paired seed-parameter values.",
      top3_label,
      trimws(axis_cfg$subtitle_note)
    )[nzchar(c(
      subtitle_line1,
      "In vivo points are colored by Vxx warm-up family; in vitro points are blue; lines connect paired seed-parameter values.",
      top3_label,
      trimws(axis_cfg$subtitle_note)
    ))],
    collapse = " "
  )
  subtitle_text <- paste(strwrap(subtitle_text, width = 125L), collapse = "\n")

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = boundary_x_plot, y = y_plot))
  if (identical(axis_cfg$axis_type, "relative")) {
    p <- p +
      ggplot2::annotate("rect", xmin = axis_cfg$lower_rect[[1]], xmax = axis_cfg$lower_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#fddbc7", alpha = 0.28) +
      ggplot2::annotate("rect", xmin = axis_cfg$upper_rect[[1]], xmax = axis_cfg$upper_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#d1e5f0", alpha = 0.28) +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = boundary_x_start, xend = boundary_x_end, y = y_base, yend = y_base),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      ) +
      ggplot2::geom_vline(xintercept = axis_cfg$vlines, color = axis_cfg$vline_colors, linetype = axis_cfg$vline_linetypes, linewidth = 0.35)
  } else {
    p <- p +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = boundary_x_lower, xend = boundary_x_upper, y = y_base, yend = y_base),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      )
  }
  if (nrow(other_pair_df)) {
    p <- p +
      ggplot2::geom_segment(
        data = other_pair_df,
        ggplot2::aes(
          x = boundary_x_plot_vivo,
          xend = boundary_x_plot_vitro,
          y = y_plot_vivo,
          yend = y_plot_vitro
        ),
        inherit.aes = FALSE,
        color = "grey55",
        alpha = 0.12,
        linewidth = 0.18
      )
  }
  if (nrow(top_pair_df)) {
    p <- p +
      ggplot2::geom_segment(
        data = top_pair_df,
        ggplot2::aes(
          x = boundary_x_plot_vivo,
          xend = boundary_x_plot_vitro,
          y = y_plot_vivo,
          yend = y_plot_vitro
        ),
        inherit.aes = FALSE,
        color = "grey25",
        alpha = 0.70,
        linewidth = 0.45
      )
  }

  p <- p +
    ggplot2::geom_point(
      data = other_vitro_df,
      shape = 16,
      size = 2.1,
      color = "#0072B2",
      alpha = 0.42
    ) +
    ggplot2::geom_point(
      data = other_vivo_df,
      ggplot2::aes(color = invivo_warmup),
      shape = 16,
      size = 2.1,
      alpha = 0.62
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = levels(plot_df$parameter),
      expand = ggplot2::expansion(add = 0.45)
    ) +
    axis_cfg$scale +
    ggplot2::labs(
      title = paste0("Joint Soft-Coupled In Vivo/In Vitro Paired Parameter Positions", title_detail),
      subtitle = subtitle_text,
      x = axis_cfg$x_label,
      y = NULL,
      color = "In vivo warm-up"
    ) +
    scale_color_invivo_warmup(family_levels) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1), order = 1)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )

  if (nrow(top_df)) {
    p <- p +
      ggplot2::geom_point(
        data = top_df,
        ggplot2::aes(shape = seed_marker),
        size = 3.0,
        color = "black"
      ) +
      ggplot2::scale_shape_manual(values = shape_values, breaks = top_breaks, drop = FALSE) +
      ggplot2::labs(shape = if (isTRUE(pred1000_top3_only)) "Top 3 eligible seeds" else "Objective top 3 seeds") +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1), order = 1),
        shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 2)
      )
  }
  save_plot_pdf_png(p, root, stem, width = 12, height = 6.2)
}

plot_multi_warmup_seed_trajectory <- function(root,
                                              stem,
                                              data_file,
                                              summary_file,
                                              cohort_value,
                                              y_col,
                                              summary_transform = identity,
                                              y_transform = identity,
                                              title,
                                              subtitle,
                                              y_label,
                                              add_ploidy_axis = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  extra_dir <- integrated_extra_results_dir(root)
  seed_df <- read_table_optional(file.path(extra_dir, data_file))
  mean_df <- read_table_optional(file.path(extra_dir, summary_file))
  if (is.null(seed_df) || !is.data.frame(seed_df) || !all(c("seed", "cohort", "day", y_col) %in% names(seed_df))) return(invisible(NULL))
  if (is.null(mean_df) || !is.data.frame(mean_df) || !all(c("cohort", "day", "mean_value", "median_value") %in% names(mean_df))) return(invisible(NULL))
  seed_df <- seed_df[as.character(seed_df$cohort) == cohort_value, , drop = FALSE]
  mean_df <- mean_df[as.character(mean_df$cohort) == cohort_value, , drop = FALSE]
  if (!nrow(seed_df) || !nrow(mean_df)) return(invisible(NULL))
  seed_df <- add_invivo_warmup_from_seed(seed_df, root)
  seed_df$day <- num_col(seed_df$day)
  seed_df$plot_value <- y_transform(num_col(seed_df[[y_col]]))
  mean_df$day <- num_col(mean_df$day)
  mean_df$plot_mean <- summary_transform(num_col(mean_df$mean_value))
  mean_df$plot_median <- summary_transform(num_col(mean_df$median_value))
  seed_df <- seed_df[is.finite(seed_df$day) & is.finite(seed_df$plot_value), , drop = FALSE]
  mean_df <- mean_df[is.finite(mean_df$day) & is.finite(mean_df$plot_mean), , drop = FALSE]
  if (!nrow(seed_df) || !nrow(mean_df)) return(invisible(NULL))
  summary_line_color <- if (identical(cohort_value, "2N")) "#2166AC" else "#B2182B"
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = seed_df,
      ggplot2::aes(x = day, y = plot_value, group = seed, color = invivo_warmup),
      linewidth = 0.28,
      alpha = 0.58
    ) +
    ggplot2::geom_line(
      data = mean_df,
      ggplot2::aes(x = day, y = plot_mean),
      color = summary_line_color,
      linewidth = 1.1
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Day",
      y = y_label,
      color = "In vivo warm-up"
    ) +
    scale_color_invivo_warmup(levels(seed_df$invivo_warmup)) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.84, alpha = 1))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
  if (isTRUE(add_ploidy_axis)) {
    p <- p +
      ggplot2::geom_line(
        data = mean_df[is.finite(mean_df$plot_median), , drop = FALSE],
        ggplot2::aes(x = day, y = plot_median),
        color = summary_line_color,
        linewidth = 0.75,
        linetype = "dashed"
      ) +
      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ . / 22, name = "Ploidy"))
  }
  save_plot_pdf_png(p, root, stem, width = 12, height = 6)
}

plot_integrated_colored_seed_figures <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_multi_warmup_boundary_forest(
    root,
    stem = "multi_warmup_integrated_parameter_boundary_forest_by_invivo_warmup",
    x_scale = "relative"
  )
  plot_multi_warmup_paired_boundary_scatter(
    root,
    stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_by_invivo_warmup",
    x_scale = "relative"
  )
  plot_multi_warmup_boundary_forest(
    root,
    stem = "multi_warmup_integrated_parameter_boundary_forest_pred1000_gt44_top3_by_invivo_warmup",
    title_suffix = "Pred1000 > 44 Top 3",
    pred1000_top3_only = TRUE,
    x_scale = "relative"
  )
  plot_multi_warmup_paired_boundary_scatter(
    root,
    stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_pred1000_gt44_top3_by_invivo_warmup",
    title_suffix = "Pred1000 > 44 Top 3",
    pred1000_top3_only = TRUE,
    x_scale = "relative"
  )
  plot_multi_warmup_boundary_forest(
    root,
    stem = "multi_warmup_integrated_parameter_boundary_forest_log_x_by_invivo_warmup",
    title_suffix = "Log X Scale",
    x_scale = "log10_original"
  )
  plot_multi_warmup_paired_boundary_scatter(
    root,
    stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_log_x_by_invivo_warmup",
    title_suffix = "Log X Scale",
    x_scale = "log10_original"
  )
  plot_multi_warmup_boundary_forest(
    root,
    stem = "multi_warmup_integrated_parameter_boundary_forest_pred1000_gt44_top3_log_x_by_invivo_warmup",
    title_suffix = "Pred1000 > 44 Top 3, Log X Scale",
    pred1000_top3_only = TRUE,
    x_scale = "log10_original"
  )
  plot_multi_warmup_paired_boundary_scatter(
    root,
    stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_pred1000_gt44_top3_log_x_by_invivo_warmup",
    title_suffix = "Pred1000 > 44 Top 3, Log X Scale",
    pred1000_top3_only = TRUE,
    x_scale = "log10_original"
  )
  plot_multi_warmup_seed_trajectory(
    root,
    stem = "multi_warmup_integrated_burden_total_log_seed_mean_2N_by_invivo_warmup",
    data_file = "predict1000_burden_total_seed_day_mean.tsv",
    summary_file = "predict1000_burden_total_mean_ci.tsv",
    cohort_value = "2N",
    y_col = "burden_value",
    y_transform = function(x) log10(x),
    summary_transform = function(x) log10(x),
    title = "1000-day Total Burden Prediction by In Vivo Warm-Up Family: 2N",
    subtitle = "Colored hairlines = individual seed trajectories by Vxx warm-up source; solid line = cross-seed mean.",
    y_label = "log10 total tumor burden (mm^3)"
  )
  plot_multi_warmup_seed_trajectory(
    root,
    stem = "multi_warmup_integrated_burden_total_log_seed_mean_4N_by_invivo_warmup",
    data_file = "predict1000_burden_total_seed_day_mean.tsv",
    summary_file = "predict1000_burden_total_mean_ci.tsv",
    cohort_value = "4N",
    y_col = "burden_value",
    y_transform = function(x) log10(x),
    summary_transform = function(x) log10(x),
    title = "1000-day Total Burden Prediction by In Vivo Warm-Up Family: 4N",
    subtitle = "Colored hairlines = individual seed trajectories by Vxx warm-up source; solid line = cross-seed mean.",
    y_label = "log10 total tumor burden (mm^3)"
  )
  plot_multi_warmup_seed_trajectory(
    root,
    stem = "multi_warmup_integrated_ploidy_seed_curves_2N_by_invivo_warmup",
    data_file = "predict1000_ploidy_seed_day_mean.tsv",
    summary_file = "predict1000_ploidy_mean_ci.tsv",
    cohort_value = "2N",
    y_col = "ploidy_value",
    title = "1000-day Ploidy Seed Trajectories by In Vivo Warm-Up Family: 2N",
    subtitle = "Colored hairlines = individual seed trajectories by Vxx warm-up source; solid line = cross-seed mean; dashed line = median.",
    y_label = "Mean chromosome count",
    add_ploidy_axis = TRUE
  )
  plot_multi_warmup_seed_trajectory(
    root,
    stem = "multi_warmup_integrated_ploidy_seed_curves_4N_by_invivo_warmup",
    data_file = "predict1000_ploidy_seed_day_mean.tsv",
    summary_file = "predict1000_ploidy_mean_ci.tsv",
    cohort_value = "4N",
    y_col = "ploidy_value",
    title = "1000-day Ploidy Seed Trajectories by In Vivo Warm-Up Family: 4N",
    subtitle = "Colored hairlines = individual seed trajectories by Vxx warm-up source; solid line = cross-seed mean; dashed line = median.",
    y_label = "Mean chromosome count",
    add_ploidy_axis = TRUE
  )
  invisible(NULL)
}


render_multi_warmup_report_figures <- function(root) {
  root <- normalizePath(path.expand(root), mustWork = TRUE)
  plot_part1_color_aligned_figures(root)
  plot_integrated_tradeoff_by_invivo_warmup(root)
  plot_integrated_colored_seed_figures(root)
  invisible(TRUE)
}
