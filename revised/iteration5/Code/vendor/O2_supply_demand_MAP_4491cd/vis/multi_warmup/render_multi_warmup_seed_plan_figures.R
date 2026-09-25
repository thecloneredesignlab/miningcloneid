#!/usr/bin/env Rscript

# Pure figure rendering for materialized multi-warmup seed-plan tables.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to render multi-warmup seed-plan figures.", call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))

.o2mw_seed_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_multi_warmup_seed_plan_figures.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_seed_workflow_root <- normalizePath(file.path(.o2mw_seed_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_seed_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

square_limits <- function(x, y, padding_frac = 0.08) {
  xr <- range(x, finite = TRUE)
  yr <- range(y, finite = TRUE)
  span <- max(diff(xr), diff(yr), 1e-6)
  span <- span * (1 + padding_frac)
  xmid <- mean(xr)
  ymid <- mean(yr)
  list(
    x = xmid + c(-0.5, 0.5) * span,
    y = ymid + c(-0.5, 0.5) * span
  )
}

wrap_label <- function(x, width = 86L) {
  paste(strwrap(x, width = width), collapse = "\n")
}

clear_removed_ratio_umap_outputs <- function(out_dir) {
  stale <- file.path(
    out_dir,
    c(
      "cross_paired_top10_ratio_matrix.tsv",
      "cross_paired_top10_umap_coords.tsv",
      "cross_paired_top10_umap_coords_by_invivo_cluster.tsv",
      "joint_soft_coupling_ratio_umap_500seed.pdf",
      "joint_soft_coupling_ratio_umap_500seed.png",
      "joint_soft_coupling_ratio_umap_by_invivo_cluster_500seed.pdf",
      "joint_soft_coupling_ratio_umap_by_invivo_cluster_500seed.png",
      "joint_soft_coupling_invivo_cluster_umap.pdf",
      "joint_soft_coupling_invivo_cluster_umap.png",
      "joint_soft_coupling_invivo_cluster_umap_coords.tsv"
    )
  )
  unlink(stale[file.exists(stale)], force = TRUE)
}

point_shape_values <- function(top_n) {
  values <- rep(c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4), length.out = top_n)
  names(values) <- as.character(seq_len(top_n))
  values
}

plot_paired_profile_umap <- function(coords, out_dir, invivo_top_n, invitro_top_n, umap_seed) {
  required <- c("pair_id", "invivo_rank", "invitro_rank", "UMAP1", "UMAP2")
  missing <- setdiff(required, names(coords))
  if (length(missing)) stop("18D UMAP coords are missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

  plot_df <- coords
  plot_df$invitro_rank_factor <- factor(plot_df$invitro_rank, levels = seq_len(invitro_top_n))
  limits <- square_limits(plot_df$UMAP1, plot_df$UMAP2)
  p <- ggplot(plot_df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = invivo_rank, shape = invitro_rank_factor), size = 2.8, stroke = 0.9) +
    scale_color_gradient(low = "#2b6cb0", high = "#d33f3f", name = "In Vivo Rank") +
    scale_shape_manual(values = point_shape_values(invitro_top_n), name = "In Vitro Rank") +
    labs(
      title = "18D Warm-Start Profile UMAP",
      subtitle = wrap_label(
        paste0(nrow(plot_df), " source top", invivo_top_n, " x top", invitro_top_n, " pairings; scaled 18D warm-start profile features."),
        width = 68L
      ),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE, clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(aes(label = pair_id), size = 2.1, color = "grey15", max.overlaps = Inf, seed = umap_seed)
  } else {
    p <- p + geom_text(aes(label = pair_id), size = 2.1, vjust = -0.6)
  }

  ggsave(file.path(out_dir, "joint_soft_coupling_18d_profile_umap_500seed.pdf"), p, width = 7, height = 7, bg = "white")
  ggsave(file.path(out_dir, "joint_soft_coupling_18d_profile_umap_500seed.png"), p, width = 7, height = 7, dpi = 220, bg = "white")
  plot_df
}

plot_paired_profile_umap_by_invivo_cluster <- function(coords, invivo_means, invivo_reps, out_dir, invitro_top_n, umap_seed) {
  required <- c("pair_id", "invivo_rank", "invitro_rank", "UMAP1", "UMAP2")
  missing <- setdiff(required, names(coords))
  if (length(missing)) stop("18D UMAP coords are missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  missing_means <- setdiff(c("rank", "cluster"), names(invivo_means))
  if (length(missing_means)) stop("In vivo cluster table is missing columns: ", paste(missing_means, collapse = ", "), call. = FALSE)

  plot_df <- coords
  cluster_idx <- match(plot_df$invivo_rank, invivo_means$rank)
  if (any(is.na(cluster_idx))) stop("Could not map paired UMAP in vivo ranks to cluster assignments.", call. = FALSE)
  cluster_levels <- sort(unique(as.integer(invivo_means$cluster)))
  plot_df$invivo_cluster <- factor(
    sprintf("C%02d", as.integer(invivo_means$cluster[cluster_idx])),
    levels = sprintf("C%02d", cluster_levels)
  )
  rep_ranks <- if (!is.null(invivo_reps) && "rank" %in% names(invivo_reps)) invivo_reps$rank else integer(0)
  plot_df$is_invivo_representative <- plot_df$invivo_rank %in% rep_ranks
  plot_df$invitro_rank_factor <- factor(plot_df$invitro_rank, levels = seq_len(invitro_top_n))

  cluster_colors <- family_palette(sprintf("C%02d", cluster_levels))
  limits <- square_limits(plot_df$UMAP1, plot_df$UMAP2)

  p <- ggplot(plot_df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = invivo_cluster, shape = invitro_rank_factor), size = 3.0, stroke = 0.9) +
    scale_color_manual(values = cluster_colors, name = "In Vivo Cluster") +
    scale_shape_manual(values = point_shape_values(invitro_top_n), name = "In Vitro Rank") +
    guides(
      color = guide_legend(nrow = 1, override.aes = list(size = 3.2)),
      shape = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    labs(
      title = "18D Warm-Start Profile UMAP by In Vivo Cluster",
      subtitle = wrap_label(
        "Same 18D paired warm-start profile UMAP as the standalone view, colored by profile-space in vivo cluster. Red labels mark selected representative in vivo ranks.",
        width = 68L
      ),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE, clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical"
    )
  label_nonrep <- plot_df[!plot_df$is_invivo_representative, , drop = FALSE]
  label_rep <- plot_df[plot_df$is_invivo_representative, , drop = FALSE]
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    if (nrow(label_nonrep)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_nonrep,
        aes(label = pair_id),
        size = 2.1,
        color = "grey15",
        max.overlaps = Inf,
        seed = umap_seed
      )
    }
    if (nrow(label_rep)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_rep,
        aes(label = pair_id),
        size = 2.35,
        color = "#d62728",
        fontface = "bold",
        max.overlaps = Inf,
        seed = umap_seed
      )
    }
  } else {
    if (nrow(label_nonrep)) {
      p <- p + geom_text(data = label_nonrep, aes(label = pair_id), size = 2.1, color = "grey15", vjust = -0.6)
    }
    if (nrow(label_rep)) {
      p <- p + geom_text(data = label_rep, aes(label = pair_id), size = 2.35, color = "#d62728", fontface = "bold", vjust = -0.6)
    }
  }

  write_tsv(plot_df, file.path(out_dir, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_coords.tsv"))
  ggsave(file.path(out_dir, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_500seed.pdf"), p, width = 7, height = 7, bg = "white")
  ggsave(file.path(out_dir, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_500seed.png"), p, width = 7, height = 7, dpi = 220, bg = "white")
  plot_df
}

plot_single_side_umap <- function(coords, side, out_dir, umap_seed) {
  side <- match.arg(side, c("invivo", "invitro"))
  side_title <- if (identical(side, "invivo")) "In Vivo" else "In Vitro"
  rank_prefix <- if (identical(side, "invivo")) "V" else "I"
  cluster_prefix <- if (identical(side, "invivo")) "C" else "I"
  plot_df <- coords
  plot_df$rank_label <- sprintf("%s%02d", rank_prefix, as.integer(plot_df$rank))
  cluster_levels <- sort(unique(as.integer(plot_df$cluster)))
  plot_df$cluster_label <- factor(
    sprintf("%s%02d", cluster_prefix, as.integer(plot_df$cluster)),
    levels = sprintf("%s%02d", cluster_prefix, cluster_levels)
  )
  cluster_colors <- family_palette(levels(plot_df$cluster_label))
  limits <- square_limits(plot_df$UMAP1, plot_df$UMAP2)

  p <- ggplot(plot_df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = cluster_label), size = 3.0, stroke = 0.9) +
    scale_color_manual(values = cluster_colors, name = paste(side_title, "Cluster")) +
    labs(
      title = paste(side_title, "Warm-Start Profile UMAP"),
      subtitle = wrap_label(
        paste0(nrow(plot_df), " source ", side_title, " seeds clustered directly from optimizer-scale best-parameter features."),
        width = 68L
      ),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE, clip = "off") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(aes(label = rank_label), size = 2.1, color = "grey15", max.overlaps = Inf, seed = umap_seed)
  } else {
    p <- p + geom_text(aes(label = rank_label), size = 2.1, vjust = -0.6)
  }

  ggsave(file.path(out_dir, paste0("multi_warmup_", side, "_only_profile_umap.pdf")), p, width = 7, height = 7, bg = "white")
  ggsave(file.path(out_dir, paste0("multi_warmup_", side, "_only_profile_umap.png")), p, width = 7, height = 7, dpi = 220, bg = "white")
  plot_df
}


render_multi_warmup_seed_plan_figures <- function(out_dir, umap_seed = 1L) {
  out_dir <- normalizePath(path.expand(out_dir), mustWork = TRUE)
  clear_removed_ratio_umap_outputs(out_dir)
  paired_path <- file.path(out_dir, "joint_soft_coupling_18d_profile_umap_coords.tsv")
  invivo_summary_path <- file.path(out_dir, "multi_warmup_invivo_cluster_summary.tsv")
  invivo_reps_path <- file.path(out_dir, "multi_warmup_invivo_representatives.tsv")
  if (file.exists(paired_path)) {
    coords <- read_tsv(paired_path)
    plot_paired_profile_umap(
      coords, out_dir,
      invivo_top_n = length(unique(coords$invivo_rank)),
      invitro_top_n = length(unique(coords$invitro_rank)),
      umap_seed = umap_seed
    )
    if (file.exists(invivo_summary_path) && file.exists(invivo_reps_path)) {
      plot_paired_profile_umap_by_invivo_cluster(
        coords, read_tsv(invivo_summary_path), read_tsv(invivo_reps_path),
        out_dir, invitro_top_n = length(unique(coords$invitro_rank)),
        umap_seed = umap_seed
      )
    }
  }
  for (side in c("invivo", "invitro")) {
    path <- file.path(out_dir, paste0("multi_warmup_", side, "_only_profile_umap_coords.tsv"))
    if (file.exists(path)) plot_single_side_umap(read_tsv(path), side, out_dir, umap_seed)
  }
  invisible(TRUE)
}
