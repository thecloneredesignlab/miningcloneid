#!/usr/bin/env Rscript

# Pure visualization of materialized multi-warmup summary and diagnostic tables.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to render multi-warmup figures.", call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))

.o2mw_collect_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_multi_warmup_collected_figures.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_collect_workflow_root <- normalizePath(file.path(.o2mw_collect_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_collect_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

plot_objectives <- function(summary_df, out_dir) {
  if (!nrow(summary_df)) return(invisible(NULL))
  family_levels <- unique(as.character(summary_df$invivo_family))
  family_levels <- family_levels[nzchar(family_levels)]
  summary_df$invivo_family <- factor(as.character(summary_df$invivo_family), levels = family_levels)
  p1 <- ggplot(summary_df, aes(x = reorder(warmup_label, objective), y = objective, color = invivo_family)) +
    geom_point(size = 2.8) +
    coord_flip() +
    scale_color_manual(values = family_palette(family_levels), drop = FALSE) +
    labs(title = "Multi-Warm-Up Best Joint Objective", x = "Warm-up pair", y = "Best total objective", color = "In vivo family") +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "multi_warmup_objective_summary.pdf"), p1, width = 8, height = 5, bg = "white")

  if (all(c("objective_invivo", "objective_invitro") %in% names(summary_df))) {
    p2 <- ggplot(summary_df, aes(objective_invivo, objective_invitro, color = invivo_family, label = warmup_label)) +
      geom_point(size = 2.8) +
      scale_color_manual(values = family_palette(family_levels), drop = FALSE) +
      labs(title = "In Vivo vs In Vitro Objective by Warm-Up Pair", x = "Best seed in vivo objective", y = "Best seed in vitro objective", color = "In vivo family") +
      theme_minimal(base_size = 11)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p2 <- p2 + ggrepel::geom_text_repel(size = 3, max.overlaps = Inf)
    } else {
      p2 <- p2 + geom_text(size = 3, vjust = -0.6)
    }
    ggsave(file.path(out_dir, "multi_warmup_invivo_invitro_objective_scatter.pdf"), p2, width = 7, height = 5, bg = "white")
  }
}

plot_deoptim_iteration_histograms <- function(iter_df, out_dir) {
  if (is.null(iter_df) || !is.data.frame(iter_df) || !nrow(iter_df)) return(invisible(NULL))
  if (!all(c("warmup_label", "optimizer_iter_completed", "seed_frequency") %in% names(iter_df))) {
    stop("DEoptim plot requires the materialized iteration-count table.", call. = FALSE)
  }
  iter_df$optimizer_iter_completed <- num_col(iter_df$optimizer_iter_completed)
  iter_df <- iter_df[is.finite(iter_df$optimizer_iter_completed), , drop = FALSE]
  if (!nrow(iter_df)) return(invisible(NULL))
  pair_levels <- unique(iter_df$warmup_label)
  iter_df$warmup_label <- factor(iter_df$warmup_label, levels = pair_levels)
  counts <- iter_df
  n_pairs <- length(pair_levels)
  ncol_use <- min(3L, max(1L, n_pairs))
  height_use <- max(4.8, 2.4 * ceiling(n_pairs / ncol_use))
  width_use <- max(8, 3.8 * ncol_use)
  iter_breaks <- sort(unique(counts$optimizer_iter_completed))
  bar_width <- if (length(iter_breaks) <= 1L) 0.35 else 0.85
  x_scale <- if (length(iter_breaks) <= 12L) {
    scale_x_continuous(breaks = iter_breaks)
  } else {
    scale_x_continuous(n.breaks = 5)
  }
  p <- ggplot(counts, aes(x = optimizer_iter_completed, y = seed_frequency)) +
    geom_col(width = bar_width, fill = "#2b6cb0", alpha = 0.82) +
    facet_wrap(~warmup_label, scales = "free_y", ncol = ncol_use) +
    x_scale +
    labs(
      title = "DEoptim Iteration Frequency by Warm-Up Pair",
      x = "Completed DEoptim iterations",
      y = "Seed frequency"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
  ggsave(file.path(out_dir, "multi_warmup_deoptim_iteration_histograms.pdf"), p, width = width_use, height = height_use, bg = "white")
}

plot_invivo_only_warmstart_parameters <- function(invivo_only_long, out_dir) {
  if (is.null(invivo_only_long) || !is.data.frame(invivo_only_long) || !nrow(invivo_only_long)) return(invisible(NULL))
  plot_df <- invivo_only_long[is.finite(invivo_only_long$scaled_log10_value), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  n_seed <- length(unique(plot_df$invivo_label))
  height_use <- max(4.5, 0.42 * n_seed + 2.2)
  p <- ggplot(plot_df, aes(parameter, invivo_label, fill = scaled_log10_value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = display_value, color = label_color), size = 2.7) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#d33f3f", midpoint = 0, name = "scaled\nlog10 value") +
    scale_color_identity() +
    labs(
      title = "In Vivo-Only Warm-Start Parameters",
      subtitle = "Source in vivo seed values for parameters not represented in paired soft-coupling ratios.",
      x = "Parameter",
      y = "In vivo warm-start seed"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  ggsave(file.path(out_dir, "multi_warmup_invivo_only_parameter_heatmap.pdf"), p, width = 7.5, height = height_use, bg = "white")
}

plot_parameter_ratios <- function(ratio_long, out_dir) {
  if (is.null(ratio_long) || !nrow(ratio_long)) return(invisible(NULL))
  if (!"final_log10_ratio_vivo_to_vitro" %in% names(ratio_long)) {
    stop("Parameter-ratio plot requires final_log10_ratio_vivo_to_vitro from analysis.", call. = FALSE)
  }
  ratio_long <- ratio_long[is.finite(ratio_long$final_log10_ratio_vivo_to_vitro), , drop = FALSE]
  if (!nrow(ratio_long)) return(invisible(NULL))
  p <- ggplot(ratio_long, aes(parameter, warmup_label, fill = final_log10_ratio_vivo_to_vitro)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#d33f3f", midpoint = 0, name = "log10 ratio") +
    labs(title = "Best-Seed Final In Vivo / In Vitro Parameter Ratios", x = "Parameter", y = "Warm-up pair") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "multi_warmup_final_parameter_ratio_heatmap.pdf"), p, width = 9, height = 5.5, bg = "white")
}

plot_initial_final_distance_diagnostics <- function(per_seed, heatmap_long, out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (is.data.frame(per_seed) && nrow(per_seed)) {
    pair_levels <- if ("warmup_order" %in% names(per_seed)) {
      unique(as.character(per_seed$warmup_label[order(per_seed$warmup_order)]))
    } else {
      unique(as.character(per_seed$warmup_label))
    }
    per_seed$warmup_label <- factor(as.character(per_seed$warmup_label), levels = pair_levels)
    per_seed$final_basin_id <- factor(as.character(per_seed$final_basin_id), levels = unique(as.character(per_seed$final_basin_id)))
    p <- ggplot(per_seed, aes(x = warmup_label, y = distance_to_own_warmup, color = final_basin_id)) +
      geom_jitter(width = 0.18, height = 0, size = 2.2, alpha = 0.82) +
      labs(
        title = "Initial-to-Final 18D Profile Distance by Joint Seed",
        subtitle = "Distance is computed in scaled 18D warm-start profile space.",
        x = "Warm-up pair",
        y = "Distance to own warm-up profile",
        color = "Final basin"
      ) +
      theme_minimal(base_size = 11) +
      guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3, alpha = 1))) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom"
      )
    square_size <- 7.2
    ggsave(file.path(out_dir, "multi_warmup_initial_to_final_distance.pdf"), p, width = square_size, height = square_size, bg = "white")
    ggsave(file.path(out_dir, "multi_warmup_initial_to_final_distance.png"), p, width = square_size, height = square_size, dpi = 220, bg = "white")
  }

  if (is.data.frame(heatmap_long) && nrow(heatmap_long)) {
    heatmap_long$final_representative <- factor(
      as.character(heatmap_long$final_representative),
      levels = rev(unique(as.character(heatmap_long$final_representative)))
    )
    heatmap_long$warmup_label <- factor(as.character(heatmap_long$warmup_label), levels = unique(as.character(heatmap_long$warmup_label)))
    heatmap_long$distance_label <- format_heatmap_value(heatmap_long$distance)
    p <- ggplot(heatmap_long, aes(x = warmup_label, y = final_representative, fill = distance)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = distance_label), size = 2.8) +
      scale_fill_gradient(low = "#f7fbff", high = "#4a1486", guide = "none") +
      labs(
        title = "Final Basin Representative vs Warm-Up 18D Profile Distance",
        subtitle = "Rows are best-objective representatives from each final basin;\ncolumns are initial warm-up profiles.",
        x = "Warm-up profile",
        y = "Final basin representative"
      ) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
    square_size <- 7.2
    ggsave(file.path(out_dir, "multi_warmup_final_to_warmup_distance_heatmap.pdf"), p, width = square_size, height = square_size, bg = "white")
    ggsave(file.path(out_dir, "multi_warmup_final_to_warmup_distance_heatmap.png"), p, width = square_size, height = square_size, dpi = 220, bg = "white")
  }
}

cleanup_removed_basin_pca_files <- function(out_dir) {
  unlink(
    file.path(
      out_dir,
      c("multi_warmup_final_parameter_basin_pca.pdf", "multi_warmup_final_parameter_pca_coords.tsv")
    ),
    force = TRUE
  )
}


render_multi_warmup_collected_figures <- function(out_dir) {
  out_dir <- normalizePath(path.expand(out_dir), mustWork = TRUE)
  cleanup_removed_basin_pca_files(out_dir)
  summary_df <- read_table_optional(file.path(out_dir, "multi_warmup_best_seed_summary.tsv"))
  ratio_long <- read_table_optional(file.path(out_dir, "multi_warmup_parameter_ratio_long.tsv"))
  iteration_counts <- read_table_optional(file.path(out_dir, "multi_warmup_deoptim_iteration_counts.tsv"))
  invivo_only_long <- read_table_optional(file.path(out_dir, "multi_warmup_invivo_only_parameter_long.tsv"))
  per_seed <- read_table_optional(file.path(out_dir, "multi_warmup_initial_final_distance_per_seed.tsv"))
  heatmap_long <- read_table_optional(file.path(out_dir, "multi_warmup_final_to_warmup_distance_heatmap.tsv"))
  plot_objectives(summary_df, out_dir)
  plot_deoptim_iteration_histograms(iteration_counts, out_dir)
  plot_parameter_ratios(ratio_long, out_dir)
  plot_invivo_only_warmstart_parameters(invivo_only_long, out_dir)
  plot_initial_final_distance_diagnostics(per_seed, heatmap_long, out_dir)
  invisible(TRUE)
}
