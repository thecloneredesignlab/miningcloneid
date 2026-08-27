#!/usr/bin/env Rscript

# Pure visualization consumer for warm-up joint embedding and curve-class tables.

.o2mw_warm_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_warm_up_joint_figures.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_warm_vis_root <- normalizePath(file.path(.o2mw_warm_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_warm_vis_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

make_pair_palette <- function(pair_ids) {
  pair_ids <- sort(unique(as.character(pair_ids)))
  cols <- grDevices::hcl.colors(length(pair_ids), palette = "Dark 3")
  names(cols) <- pair_ids
  cols
}

blend_hex_colors <- function(cols, grey = "grey50", color_weight = 0.5) {
  if (!length(cols)) return(cols)
  color_weight <- max(0, min(1, color_weight))
  rgb_cols <- grDevices::col2rgb(cols)
  rgb_grey <- grDevices::col2rgb(grey)
  blended <- round(color_weight * rgb_cols + (1 - color_weight) * as.numeric(rgb_grey))
  out <- grDevices::rgb(blended[1, ], blended[2, ], blended[3, ], maxColorValue = 255)
  names(out) <- names(cols)
  out
}

embedding_plot <- function(coords, reduction, x_col, y_col, out_prefix, param_label = "parameters") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed; skipping plot: ", out_prefix, call. = FALSE)
    return(invisible(NULL))
  }
  d <- coords
  d$point_type <- factor(d$point_type, levels = c("initial", "best"))
  pair_levels <- sort(unique(as.character(d$pair_id)))
  d$pair_id <- factor(d$pair_id, levels = pair_levels)
  pair_cols <- make_pair_palette(pair_levels)
  initial_cols <- blend_hex_colors(pair_cols, grey = "grey50", color_weight = 0.5)
  initial <- d[d$point_type == "initial", , drop = FALSE]
  best <- d[d$point_type == "best", , drop = FALSE]
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = initial,
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]),
      color = unname(initial_cols[as.character(initial$pair_id)]),
      alpha = 0.5, size = 0.36, stroke = 0
    ) +
    ggplot2::geom_point(
      data = best,
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]),
      shape = 8, color = "white", alpha = 0.95, size = 3.0, stroke = 0.9
    ) +
    ggplot2::geom_point(
      data = best,
      ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = pair_id),
      shape = 8, alpha = 0.95, size = 2.2, stroke = 0.65
    ) +
    ggplot2::scale_color_manual(values = pair_cols, drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::labs(
      title = paste0("Joint in vivo ", reduction),
      subtitle = paste0("Regenerated initial population vs best; ", param_label),
      x = x_col, y = y_col, color = "Pair"
    )
  ggplot2::ggsave(paste0(out_prefix, ".pdf"), p, width = 7.2, height = 7.2, useDingbats = FALSE)
  ggplot2::ggsave(paste0(out_prefix, ".png"), p, width = 7.2, height = 7.2, dpi = 300)

  p_best <- ggplot2::ggplot(best, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = pair_id)) +
    ggplot2::geom_point(shape = 8, color = "white", alpha = 0.95, size = 3.0, stroke = 0.9) +
    ggplot2::geom_point(shape = 8, alpha = 0.95, size = 2.2, stroke = 0.65) +
    ggplot2::scale_color_manual(values = pair_cols, drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::labs(
      title = paste0("Joint in vivo ", reduction, ": best only"),
      subtitle = param_label,
      x = x_col, y = y_col, color = "Pair"
    )
  ggplot2::ggsave(paste0(out_prefix, "_best_only.pdf"), p_best, width = 7.2, height = 7.2, useDingbats = FALSE)
  ggplot2::ggsave(paste0(out_prefix, "_best_only.png"), p_best, width = 7.2, height = 7.2, dpi = 300)
  invisible(NULL)
}

plot_pair_curve_summary <- function(summary_df, out_root) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !nrow(summary_df)) return(invisible(NULL))
  fig_dir <- file.path(out_root, "figures", "curve_classification")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = pair_id, y = fraction, fill = curve_class)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = "Pair", y = "Fraction of joint seeds", fill = "Curve class", title = "Joint in vivo curve class composition by pair")
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_curve_class_fraction_by_pair.pdf"), p, width = 8.4, height = 4.8, useDingbats = FALSE)
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_curve_class_fraction_by_pair.png"), p, width = 8.4, height = 4.8, dpi = 300)

  p2 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = curve_class, y = pair_id, fill = enrichment)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b83232", midpoint = 1, na.value = "grey90") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = "Curve class", y = "Pair", fill = "Observed / expected", title = "Curve class enrichment by pair")
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_curve_class_enrichment_heatmap.pdf"), p2, width = 8.4, height = 4.8, useDingbats = FALSE)
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_curve_class_enrichment_heatmap.png"), p2, width = 8.4, height = 4.8, dpi = 300)
  invisible(NULL)
}

plot_embedding_by_curve_class <- function(master, out_root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  emb_dir <- file.path(out_root, "embeddings")
  fig_dir <- file.path(out_root, "figures", "embedding_curve_class")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  coord_paths <- list.files(emb_dir, pattern = "_coordinates\\.tsv$", full.names = TRUE, recursive = TRUE)
  if (!length(coord_paths)) return(invisible(NULL))
  for (path in coord_paths) {
    coords <- read_tsv(path)
    best <- coords[coords$point_type == "best", , drop = FALSE]
    if (!nrow(best)) next
    best <- merge(
      best,
      master[, c("synthetic_seed_id", "curve_class", "smooth_curve_class", "pair_id"), drop = FALSE],
      by = c("synthetic_seed_id", "pair_id"),
      all.x = TRUE,
      sort = FALSE
    )
    coord_cols <- intersect(c("PCA1", "PCA2", "UMAP1", "UMAP2", "tSNE1", "tSNE2"), names(best))
    if (length(coord_cols) < 2L) next
    x_col <- coord_cols[[1L]]
    y_col <- coord_cols[[2L]]
    emb_base <- paste0(normalizePath(emb_dir, mustWork = FALSE), "/")
    path_norm <- normalizePath(path, mustWork = FALSE)
    rel_path <- if (startsWith(path_norm, emb_base)) substring(path_norm, nchar(emb_base) + 1L) else basename(path_norm)
    stub <- safe_id(sub("_coordinates\\.tsv$", "", gsub("/", "_", rel_path, fixed = TRUE)))
    p <- ggplot2::ggplot(best, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], color = curve_class)) +
      ggplot2::geom_point(shape = 8, alpha = 0.95, size = 2.2, stroke = 0.65) +
      ggplot2::coord_fixed() +
      ggplot2::facet_wrap(~ pair_id) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::labs(x = x_col, y = y_col, color = "Curve class", title = paste0(stub, ": best points by curve class"))
    ggplot2::ggsave(file.path(fig_dir, paste0(stub, "_best_by_curve_class.pdf")), p, width = 8.2, height = 8.2, useDingbats = FALSE)
    ggplot2::ggsave(file.path(fig_dir, paste0(stub, "_best_by_curve_class.png")), p, width = 8.2, height = 8.2, dpi = 300)
  }
  invisible(NULL)
}


plot_smoothed_curves_by_pair_and_class <- function(out_root, master) {
  curves_path <- file.path(
    out_root, "curve_classification", "dense-grid_monotonicity_regression_classification",
    "tables", "fixed_o2_ploidy_monotonicity_regression_curves.tsv"
  )
  if (!file.exists(curves_path) || !nrow(master) || !requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  curves <- read_tsv(curves_path)
  required <- c("synthetic_seed_id", "pair_id", "joint_seed", "curve_class")
  if (!all(required %in% names(master))) return(invisible(NULL))
  curves <- merge(
    curves, master[, required, drop = FALSE],
    by.x = "seed_id", by.y = "synthetic_seed_id", all.x = TRUE, sort = FALSE
  )
  fig_dir <- file.path(out_root, "figures", "curve_classification")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  figure <- ggplot2::ggplot(
    curves,
    ggplot2::aes(x = O2_pct, y = smoothed_dominant_mean_ploidy, group = seed_id, color = .data$curve_class)
  ) +
    ggplot2::geom_line(alpha = 0.22, linewidth = 0.25) +
    ggplot2::facet_wrap(~ pair_id) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(
      x = "O2 (%)", y = "Smoothed dominant mean ploidy", color = "Curve class",
      title = "Joint in vivo fixed-O2 curves by pair"
    )
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_smoothed_curves_by_pair_and_class.pdf"), figure, width = 10, height = 6.5, useDingbats = FALSE)
  ggplot2::ggsave(file.path(fig_dir, "joint_invivo_smoothed_curves_by_pair_and_class.png"), figure, width = 10, height = 6.5, dpi = 300)
  invisible(TRUE)
}

render_warm_up_joint_figures <- function(out_root) {
  out_root <- normalizePath(path.expand(out_root), mustWork = TRUE)
  embedding_tables <- list.files(
    file.path(out_root, "embeddings"), recursive = TRUE, full.names = TRUE,
    pattern = "_coordinates\\.tsv$"
  )
  for (path in embedding_tables) {
    coords <- read_tsv(path)
    pair <- if (all(c("PCA1", "PCA2") %in% names(coords))) c("PCA1", "PCA2") else
      if (all(c("UMAP1", "UMAP2") %in% names(coords))) c("UMAP1", "UMAP2") else
        if (all(c("tSNE1", "tSNE2") %in% names(coords))) c("tSNE1", "tSNE2") else character()
    if (!length(pair)) next
    token <- basename(dirname(path))
    out_prefix <- file.path(out_root, "figures", "embeddings", token, sub("_coordinates.tsv$", "", basename(path)))
    embedding_plot(coords, toupper(sub(".*_(pca|umap|tsne).*", "\\1", basename(path))),
                   pair[[1L]], pair[[2L]], out_prefix, param_label = token)
  }
  summary_path <- file.path(out_root, "tables", "pair_curve_class_summary.tsv")
  master_path <- file.path(out_root, "tables", "joint_best_master_table.tsv")
  summary_df <- if (file.exists(summary_path)) read_tsv(summary_path) else data.frame()
  master <- if (file.exists(master_path)) read_tsv(master_path) else data.frame()
  plot_pair_curve_summary(summary_df, out_root)
  plot_embedding_by_curve_class(master, out_root)
  plot_smoothed_curves_by_pair_and_class(out_root, master)
  invisible(TRUE)
}
