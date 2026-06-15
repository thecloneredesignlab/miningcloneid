#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos > 0L) out[[substr(kv, 1L, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv)) else out[[kv]] <- TRUE
  }
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) y else x
as_chr <- function(x, default = "") as.character((x %||% default)[[1]])

read_table_optional <- function(path, sep = "\t") {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}

num_col <- function(x) suppressWarnings(as.numeric(x))

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

choose_k <- function(x, k_max = 6L) {
  n <- nrow(x)
  if (n <= 2L) return(1L)
  d <- stats::dist(x)
  hc <- stats::hclust(d, method = "ward.D2")
  max_k <- max(2L, min(k_max, n - 1L))
  if (!requireNamespace("cluster", quietly = TRUE)) return(min(2L, max_k))
  scores <- sapply(2L:max_k, function(k) {
    cl <- stats::cutree(hc, k = k)
    mean(cluster::silhouette(cl, d)[, "sil_width"], na.rm = TRUE)
  })
  as.integer((2L:max_k)[which.max(scores)])
}

best_seed_row <- function(seed_summary) {
  if (is.null(seed_summary) || !nrow(seed_summary) || !("objective" %in% names(seed_summary))) return(NULL)
  seed_summary$objective <- num_col(seed_summary$objective)
  seed_summary <- seed_summary[is.finite(seed_summary$objective), , drop = FALSE]
  if (!nrow(seed_summary)) return(NULL)
  seed_summary[order(seed_summary$objective, seed_summary$seed), , drop = FALSE][1L, , drop = FALSE]
}

collect_parameter_ratios <- function(soft_tab, best_seed, warmup_label) {
  if (is.null(soft_tab) || !nrow(soft_tab)) return(NULL)
  if (!all(c("seed", "parameter") %in% names(soft_tab))) return(NULL)
  soft_tab <- soft_tab[as.character(soft_tab$seed) == as.character(best_seed), , drop = FALSE]
  if (!nrow(soft_tab)) return(NULL)
  keep <- intersect(
    c(
      "parameter", "ratio_vivo_to_vitro", "log10_ratio_vivo_to_vitro",
      "fold_change_vivo_to_vitro", "penalty_paid", "vivo_natural", "vitro_natural",
      "vivo_transformed", "vitro_transformed"
    ),
    names(soft_tab)
  )
  out <- soft_tab[, keep, drop = FALSE]
  out$warmup_label <- warmup_label
  out$best_joint_seed <- best_seed
  out
}

plot_objectives <- function(summary_df, out_dir) {
  if (!nrow(summary_df)) return(invisible(NULL))
  p1 <- ggplot(summary_df, aes(x = reorder(warmup_label, objective), y = objective, color = invivo_family)) +
    geom_point(size = 2.8) +
    coord_flip() +
    labs(title = "Multi-Warm-Up Best Joint Objective", x = "Warm-up pair", y = "Best total objective", color = "In vivo family") +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "multi_warmup_objective_summary.pdf"), p1, width = 8, height = 5, bg = "white")

  if (all(c("objective_invivo", "objective_invitro") %in% names(summary_df))) {
    p2 <- ggplot(summary_df, aes(objective_invivo, objective_invitro, color = invivo_family, label = warmup_label)) +
      geom_point(size = 2.8) +
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

plot_parameter_ratios <- function(ratio_long, out_dir) {
  if (is.null(ratio_long) || !nrow(ratio_long) || !("log10_ratio_vivo_to_vitro" %in% names(ratio_long))) return(invisible(NULL))
  ratio_long$log10_ratio_vivo_to_vitro <- num_col(ratio_long$log10_ratio_vivo_to_vitro)
  ratio_long <- ratio_long[is.finite(ratio_long$log10_ratio_vivo_to_vitro), , drop = FALSE]
  if (!nrow(ratio_long)) return(invisible(NULL))
  p <- ggplot(ratio_long, aes(parameter, warmup_label, fill = log10_ratio_vivo_to_vitro)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#d33f3f", midpoint = 0, name = "log10 ratio") +
    labs(title = "Best-Seed Final In Vivo / In Vitro Parameter Ratios", x = "Parameter", y = "Warm-up pair") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "multi_warmup_final_parameter_ratio_heatmap.pdf"), p, width = 9, height = 5.5, bg = "white")
}

assign_final_basins <- function(ratio_long, out_dir) {
  if (is.null(ratio_long) || !nrow(ratio_long) || !("log10_ratio_vivo_to_vitro" %in% names(ratio_long))) return(NULL)
  ratio_long$log10_ratio_vivo_to_vitro <- num_col(ratio_long$log10_ratio_vivo_to_vitro)
  ratio_long <- ratio_long[is.finite(ratio_long$log10_ratio_vivo_to_vitro), , drop = FALSE]
  if (!nrow(ratio_long)) return(NULL)
  wide <- reshape(
    ratio_long[, c("warmup_label", "parameter", "log10_ratio_vivo_to_vitro")],
    idvar = "warmup_label",
    timevar = "parameter",
    direction = "wide"
  )
  labels <- wide$warmup_label
  x <- as.matrix(wide[, setdiff(names(wide), "warmup_label"), drop = FALSE])
  complete <- stats::complete.cases(x)
  x <- x[complete, , drop = FALSE]
  labels <- labels[complete]
  if (!nrow(x)) return(NULL)
  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  k <- choose_k(x_scaled)
  basin <- if (k <= 1L || nrow(x_scaled) <= 2L) rep(1L, nrow(x_scaled)) else stats::cutree(stats::hclust(stats::dist(x_scaled), method = "ward.D2"), k = k)
  out <- data.frame(warmup_label = labels, final_basin_id = paste0("B", sprintf("%02d", basin)), stringsAsFactors = FALSE)
  if (nrow(x_scaled) >= 2L && ncol(x_scaled) >= 2L) {
    pc <- stats::prcomp(x_scaled, center = FALSE, scale. = FALSE)
    coords <- data.frame(warmup_label = labels, final_basin_id = out$final_basin_id, PC1 = pc$x[, 1], PC2 = pc$x[, 2], stringsAsFactors = FALSE)
    write_tsv(coords, file.path(out_dir, "multi_warmup_final_parameter_pca_coords.tsv"))
    p <- ggplot(coords, aes(PC1, PC2, color = final_basin_id, label = warmup_label)) +
      geom_point(size = 3) +
      labs(title = "Final Parameter-Ratio Basin PCA", color = "Basin") +
      theme_minimal(base_size = 11)
    if (requireNamespace("ggrepel", quietly = TRUE)) p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) else p <- p + geom_text(size = 3, vjust = -0.6)
    ggsave(file.path(out_dir, "multi_warmup_final_parameter_basin_pca.pdf"), p, width = 7, height = 5, bg = "white")
  }
  out
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root <- normalizePath(as_chr(argv$multi_warmup_root, as_chr(argv$out_dir)), mustWork = TRUE)
  manifest_path <- normalizePath(as_chr(argv$manifest, file.path(root, "multi_warmup_manifest.tsv")), mustWork = TRUE)
  out_dir <- normalizePath(as_chr(argv$out_dir, root), mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- utils::read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
  rows <- list()
  ratio_rows <- list()
  for (i in seq_len(nrow(manifest))) {
    run_dir <- file.path(root, manifest$joint_run_prefix[[i]])
    seed_summary <- read_table_optional(file.path(run_dir, "extra_results", "seed_summary.tsv"))
    best <- best_seed_row(seed_summary)
    if (is.null(best)) next
    base_cols <- intersect(
      c(
        "seed", "objective", "objective_invivo", "objective_invitro", "objective_soft_coupling",
        "objective_burden", "objective_ploidy", "objective_necrosis",
        "boundary_penalty_active", "min_rel_dist_active", "pred1000_2N", "pred1000_4N"
      ),
      names(best)
    )
    row <- cbind(
      manifest[i, , drop = FALSE],
      data.frame(joint_run_dir = run_dir, stringsAsFactors = FALSE),
      best[, base_cols, drop = FALSE]
    )
    names(row)[names(row) == "seed"] <- "best_joint_seed"
    rows[[length(rows) + 1L]] <- row
    soft <- read_table_optional(file.path(run_dir, "extra_results", "joint_soft_coupling_all.tsv"))
    ratio_tab <- collect_parameter_ratios(soft, best$seed[[1]], manifest$warmup_label[[i]])
    if (!is.null(ratio_tab)) ratio_rows[[length(ratio_rows) + 1L]] <- ratio_tab
  }
  summary_df <- if (length(rows)) do.call(rbind, rows) else data.frame()
  ratio_long <- if (length(ratio_rows)) do.call(rbind, ratio_rows) else data.frame()
  basin_df <- assign_final_basins(ratio_long, out_dir)
  if (!is.null(basin_df) && nrow(summary_df)) {
    summary_df <- merge(summary_df, basin_df, by = "warmup_label", all.x = TRUE, sort = FALSE)
  }

  if (nrow(summary_df)) {
    for (col in intersect(c("objective", "objective_invivo", "objective_invitro", "objective_soft_coupling", "objective_burden", "objective_ploidy", "objective_necrosis", "boundary_penalty_active", "min_rel_dist_active"), names(summary_df))) {
      summary_df[[col]] <- num_col(summary_df[[col]])
    }
  }
  write_tsv(summary_df, file.path(out_dir, "multi_warmup_best_seed_summary.tsv"))
  write_tsv(ratio_long, file.path(out_dir, "multi_warmup_parameter_ratio_long.tsv"))
  if (!is.null(basin_df)) write_tsv(basin_df, file.path(out_dir, "multi_warmup_final_basin_assignments.tsv"))
  plot_objectives(summary_df, out_dir)
  plot_parameter_ratios(ratio_long, out_dir)
  message("Wrote multi-warmup summary: ", file.path(out_dir, "multi_warmup_best_seed_summary.tsv"))
}

if (identical(environment(), globalenv())) {
  main()
}
