#!/usr/bin/env Rscript

# Visualization-only helpers for the post-fit process-fingerprint workflows.
# This module consumes materialized simulation/analysis tables and deliberately
# contains no fit discovery, model loading, simulation, or statistical analysis.

.pfv_utils_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[
    basename(frame_files) == "process_fingerprint_visualization_utils.R"
  ]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(
  file.path(
    .pfv_utils_dir,
    "..",
    "..",
    "util",
    "o2_supply_demand_map_process_fingerprint_utils.R"
  ),
  local = TRUE,
  chdir = TRUE
)
rm(.pfv_utils_dir)

`%||%` <- o2ipa_null_coalesce
pfv_parse_args <- o2ipa_parse_args
pfv_as_chr <- o2ipa_as_chr
pfv_read_tsv <- o2ipa_read_tsv
pfv_write_tsv <- o2ipa_write_tsv

pfv_upper_vec <- function(x) {
  if (is.null(x) || nrow(x) < 2L) return(numeric())
  x[upper.tri(x)]
}

pfv_plot_process_outputs <- function(out_dir, d_parameter, d_static, d_static18,
                                     d_dynamic, static_scaled, dynamic_scaled,
                                     static_cluster, manifest, static_long) {
  fig <- function(name) file.path(out_dir, "figures", paste0(name, ".pdf"))
  scatter <- function(a, b, xlab, ylab, name) {
    if (is.null(a) || is.null(b)) return(invisible(NULL))
    common <- intersect(rownames(a), rownames(b))
    if (length(common) < 3L) return(invisible(NULL))
    av <- pfv_upper_vec(a[common, common, drop = FALSE])
    bv <- pfv_upper_vec(b[common, common, drop = FALSE])
    grDevices::pdf(fig(name), width = 6, height = 5)
    graphics::plot(av, bv, pch = 16, xlab = xlab, ylab = ylab, main = name)
    grDevices::dev.off()
  }
  scatter(d_parameter, d_static, "Parameter distance", "Static process distance", "parameter_vs_static_process_distance")
  scatter(d_parameter, d_static18, "Parameter distance", "Static 18-only distance", "parameter_vs_static_18only_distance")
  scatter(d_static, d_dynamic, "Static process distance", "Dynamic phenotype distance", "static_process_vs_dynamic_phenotype_distance")
  scatter(d_static, d_static18, "Static full distance", "Static 18-only distance", "static_full_vs_static_18only_distance")

  if (is.data.frame(static_scaled) && nrow(static_scaled) >= 2L && ncol(static_scaled) >= 3L) {
    mat <- as.matrix(static_scaled[, setdiff(names(static_scaled), "seed_id"), drop = FALSE])
    rownames(mat) <- static_scaled$seed_id
    mat <- mat[rowSums(!is.finite(mat)) == 0, colSums(!is.finite(mat)) == 0, drop = FALSE]
    if (nrow(mat) >= 2L && ncol(mat) >= 2L) {
      grDevices::pdf(fig("static_process_feature_heatmap"), width = 9, height = 7)
      stats::heatmap(mat, scale = "none", margins = c(5, 5), main = "Static process features")
      grDevices::dev.off()
      grDevices::pdf(fig("static_process_hierarchical_dendrogram"), width = 7, height = 5)
      graphics::plot(stats::hclust(stats::dist(mat), method = "ward.D2"), main = "Static process dendrogram", xlab = "")
      grDevices::dev.off()
      pc <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
      grDevices::pdf(fig("static_process_pca"), width = 6, height = 5)
      graphics::plot(pc$x[, 1], pc$x[, 2], pch = 16, xlab = "PC1", ylab = "PC2", main = "Static process PCA")
      graphics::text(pc$x[, 1], pc$x[, 2], labels = rownames(mat), pos = 3, cex = 0.6)
      grDevices::dev.off()
    }
  }
  if (is.data.frame(dynamic_scaled) && nrow(dynamic_scaled) >= 2L && ncol(dynamic_scaled) >= 3L) {
    mat <- as.matrix(dynamic_scaled[, setdiff(names(dynamic_scaled), "seed_id"), drop = FALSE])
    rownames(mat) <- dynamic_scaled$seed_id
    mat <- mat[rowSums(!is.finite(mat)) == 0, colSums(!is.finite(mat)) == 0, drop = FALSE]
    if (nrow(mat) >= 2L && ncol(mat) >= 2L) {
      grDevices::pdf(fig("dynamic_phenotype_heatmap"), width = 8, height = 6)
      stats::heatmap(mat, scale = "none", margins = c(5, 5), main = "Dynamic phenotype features")
      grDevices::dev.off()
      pc <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
      grDevices::pdf(fig("dynamic_phenotype_pca"), width = 6, height = 5)
      graphics::plot(pc$x[, 1], pc$x[, 2], pch = 16, xlab = "PC1", ylab = "PC2", main = "Dynamic phenotype PCA")
      grDevices::dev.off()
    }
  }
  if (!is.null(static_cluster$consensus) && all(is.finite(static_cluster$consensus))) {
    grDevices::pdf(fig("bootstrap_consensus_matrix"), width = 7, height = 6)
    graphics::image(seq_len(nrow(static_cluster$consensus)), seq_len(ncol(static_cluster$consensus)), static_cluster$consensus,
                    xlab = "", ylab = "", main = "Bootstrap consensus")
    grDevices::dev.off()
  }
  if (nrow(static_cluster$membership)) {
    rec <- static_cluster$recommended_k
    memb <- static_cluster$membership[static_cluster$membership$k == rec, , drop = FALSE]
    mf <- merge(manifest, memb, by = "seed_id", all.x = FALSE)
    if (nrow(mf)) {
      grDevices::pdf(fig("objective_by_cluster"), width = 6, height = 5)
      graphics::boxplot(objective ~ cluster, data = mf, xlab = "Cluster", ylab = "Objective", main = "Objective by cluster")
      grDevices::dev.off()
      grDevices::pdf(fig("boundary_risk_by_cluster"), width = 6, height = 5)
      graphics::barplot(tapply(as.numeric(mf$boundary_risk), mf$cluster, mean, na.rm = TRUE),
                        ylab = "Boundary risk fraction", xlab = "Cluster", main = "Boundary risk by cluster")
      grDevices::dev.off()
    }
  }
  if (nrow(static_long) && nrow(static_cluster$medoids)) {
    med <- static_cluster$medoids$medoid_seed
    curve_df <- static_long[
      static_long$seed_id %in% med &
        static_long$module %in% c("hypoxia_sensing", "proliferation", "death", "CIN_missegregation") &
        grepl("_O2_|stress_at_O2", static_long$feature),
      , drop = FALSE
    ]
    if (nrow(curve_df)) {
      grDevices::pdf(fig("cluster_medoids_response_curves"), width = 8, height = 6)
      graphics::par(mfrow = c(2, 2))
      modules <- unique(curve_df$module)
      for (mod in modules[seq_len(min(4L, length(modules)))]) {
        dd <- curve_df[curve_df$module == mod, , drop = FALSE]
        graphics::stripchart(raw_value ~ seed_id, data = dd, vertical = TRUE, method = "jitter", pch = 16, main = mod, ylab = "value")
      }
      grDevices::dev.off()
    }
  }
  invisible(out_dir)
}

pfv_plot_ploidy_regime_outputs <- function(traj, scaled, labels, out_dir) {
  fig <- function(name) file.path(out_dir, "figures", paste0(name, ".pdf"))
  if (is.data.frame(traj$curves) && nrow(traj$curves)) {
    d <- merge(traj$curves, labels[, c("seed_id", "trajectory_regime")], by = "seed_id", all.x = TRUE)
    grDevices::pdf(fig("all_trajectories_by_regime"), width = 9, height = 6)
    cols <- c(stable_high_chr = "#1b9e77", late_collapse_to_low_chr = "#d95f02", ambiguous = "grey60")
    graphics::plot(NA, xlim = range(d$day, na.rm = TRUE), ylim = range(d$trajectory_value, na.rm = TRUE),
                   xlab = "Day", ylab = "Weighted mean chromosome number", main = "Trajectories by regime")
    for (key in split(seq_len(nrow(d)), paste(d$seed_id, d$cohort))) {
      x <- d[key, ]
      col <- cols[x$trajectory_regime[[1]]]
      if (!length(col) || is.na(col)) col <- "grey60"
      graphics::lines(x$day, x$trajectory_value, col = grDevices::adjustcolor(col, alpha.f = 0.35))
    }
    graphics::legend("topright", legend = names(cols), col = cols, lty = 1, bty = "n")
    grDevices::dev.off()
  }
  if (!is.null(scaled$combined_full) && nrow(scaled$combined_full$wide) >= 2L && ncol(scaled$combined_full$wide) >= 3L) {
    mat <- as.matrix(scaled$combined_full$wide[, setdiff(names(scaled$combined_full$wide), "seed_id"), drop = FALSE])
    rownames(mat) <- scaled$combined_full$wide$seed_id
    mat <- mat[rowSums(!is.finite(mat)) == 0, colSums(!is.finite(mat)) == 0, drop = FALSE]
    if (nrow(mat) >= 2L && ncol(mat) >= 2L) {
      grDevices::pdf(fig("process_PCA"), width = 6, height = 5)
      pc <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
      lab <- labels$trajectory_regime[match(rownames(mat), labels$seed_id)]
      col <- c(stable_high_chr = "#1b9e77", late_collapse_to_low_chr = "#d95f02", ambiguous = "grey60")[lab]
      graphics::plot(pc$x[, 1], pc$x[, 2], pch = 16, col = col, xlab = "PC1", ylab = "PC2", main = "Process PCA")
      grDevices::dev.off()
      grDevices::pdf(fig("intrinsic_process_fingerprint_heatmap"), width = 9, height = 7)
      stats::heatmap(mat, scale = "none", margins = c(5, 5), main = "Combined process fingerprint")
      grDevices::dev.off()
      if (nrow(mat) >= 3L) {
        grDevices::pdf(fig("process_dendrogram"), width = 7, height = 5)
        graphics::plot(stats::hclust(stats::dist(mat), method = "ward.D2"), main = "Process dendrogram", xlab = "")
        grDevices::dev.off()
      }
    }
  }
  invisible(out_dir)
}
