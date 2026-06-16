#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos > 0L) {
      out[[substr(kv, 1L, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv))
    } else {
      out[[kv]] <- TRUE
    }
  }
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) y else x

as_chr <- function(x, default = "") as.character((x %||% default)[[1]])
as_int <- function(x, default) {
  out <- suppressWarnings(as.integer(as_chr(x, as.character(default))))
  if (!is.finite(out)) default else out
}
as_bool <- function(x, default = FALSE) {
  val <- tolower(trimws(as_chr(x, if (default) "TRUE" else "FALSE")))
  val %in% c("true", "t", "1", "yes", "y", "on")
}

sanitize_label <- function(x, fallback = "warmup") {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", trimws(as.character(x)))
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) fallback else x
}

seed_order_key <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(sub("^seed", "", x)))
  out[!is.finite(out)] <- Inf
  out
}

seed_integer <- function(x) {
  out <- seed_order_key(x)
  out[!is.finite(out)] <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", x[!is.finite(out)])))
  as.integer(out)
}

ratio_parameters <- c(
  "O2_crit", "alpha_o2", "mu_hp", "p_misseg", "k_o_mis",
  "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O",
  "gamma_growth", "lam_max", "p_mis_base", "p_wgd", "gamma_mu"
)
ratio_feature_cols <- paste0("log10_ratio_", ratio_parameters)
invivo_only_oxygen_parameters <- c("o2_S0", "kappa_O", "eta_o2", "rho_2N")
invivo_only_oxygen_feature_cols <- paste0("log10_invivo_", invivo_only_oxygen_parameters)
invivo_profile_feature_cols <- c(ratio_feature_cols, invivo_only_oxygen_feature_cols)

read_metric_map <- function(path) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("metric", "value") %in% names(tab))) stop("Missing metric/value columns: ", path, call. = FALSE)
  setNames(tab$value, tab$metric)
}

metric_value <- function(vals, key, default = NA_character_) {
  if (is.null(vals) || !(key %in% names(vals))) return(default)
  vals[[key]]
}

find_seed_dirs <- function(run_dir) {
  dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "fit_summary.tsv")) & file.exists(file.path(dirs, "best_params.tsv"))]
  dirs[order(seed_order_key(basename(dirs)), basename(dirs))]
}

read_seed_summary <- function(seed_dir) {
  vals <- read_metric_map(file.path(seed_dir, "fit_summary.tsv"))
  objective <- suppressWarnings(as.numeric(metric_value(vals, "objective", metric_value(vals, "objective_total", NA))))
  data.frame(
    seed = basename(seed_dir),
    seed_int = seed_integer(basename(seed_dir)),
    seed_dir = normalizePath(seed_dir, mustWork = TRUE),
    objective = objective,
    stringsAsFactors = FALSE
  )
}

read_best_params <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("parameter", "value") %in% names(tab))) stop("Missing parameter/value columns: ", path, call. = FALSE)
  vals <- suppressWarnings(as.numeric(tab$value))
  names(vals) <- as.character(tab$parameter)
  vals
}

read_top_seed_table <- function(run_dir, top_n, label) {
  seed_dirs <- find_seed_dirs(run_dir)
  if (!length(seed_dirs)) stop("No valid seed dirs found for ", label, ": ", run_dir, call. = FALSE)
  tab <- do.call(rbind, lapply(seed_dirs, read_seed_summary))
  tab <- tab[is.finite(tab$objective), , drop = FALSE]
  if (nrow(tab) < top_n) {
    stop(label, " source run has only ", nrow(tab), " finite-objective seeds; need at least ", top_n, ".", call. = FALSE)
  }
  tab <- tab[order(tab$objective, tab$seed_int, tab$seed), , drop = FALSE]
  tab <- utils::head(tab, top_n)
  tab$rank <- seq_len(nrow(tab))
  row.names(tab) <- NULL
  tab
}

params_for_seed <- function(seed_dir) {
  vals <- read_best_params(seed_dir)
  out <- suppressWarnings(as.numeric(vals[ratio_parameters]))
  if (any(!is.finite(out))) {
    missing <- ratio_parameters[!is.finite(out)]
    stop("Missing/non-finite ratio parameters in ", seed_dir, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out
}

invivo_only_oxygen_features_for_seed <- function(seed_dir) {
  vals <- read_best_params(seed_dir)
  raw <- suppressWarnings(as.numeric(vals[invivo_only_oxygen_parameters]))
  if (any(!is.finite(raw) | raw <= 0)) {
    missing <- invivo_only_oxygen_parameters[!is.finite(raw) | raw <= 0]
    stop(
      "Missing/non-positive in vivo-only oxygen parameters in ",
      seed_dir, ": ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  out <- log10(raw)
  names(out) <- invivo_only_oxygen_feature_cols
  out
}

build_pair_tables <- function(invivo_top, invitro_top) {
  invivo_params <- lapply(invivo_top$seed_dir, params_for_seed)
  invivo_only_features <- lapply(invivo_top$seed_dir, invivo_only_oxygen_features_for_seed)
  invitro_params <- lapply(invitro_top$seed_dir, params_for_seed)
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(invivo_top))) {
    for (j in seq_len(nrow(invitro_top))) {
      ratio_vals <- invivo_params[[i]] / invitro_params[[j]]
      log_vals <- log10(ratio_vals)
      if (any(!is.finite(log_vals))) next
      row <- data.frame(
        pair_id = sprintf("V%02d-I%02d", invivo_top$rank[[i]], invitro_top$rank[[j]]),
        invivo_rank = invivo_top$rank[[i]],
        invitro_rank = invitro_top$rank[[j]],
        invivo_seed = invivo_top$seed_int[[i]],
        invitro_seed = invitro_top$seed_int[[j]],
        invivo_seed_dir = invivo_top$seed_dir[[i]],
        invitro_seed_dir = invitro_top$seed_dir[[j]],
        invivo_objective = invivo_top$objective[[i]],
        invitro_objective = invitro_top$objective[[j]],
        stringsAsFactors = FALSE
      )
      for (p in ratio_parameters) row[[paste0("ratio_", p)]] <- ratio_vals[[match(p, ratio_parameters)]]
      for (p in ratio_parameters) row[[paste0("log10_ratio_", p)]] <- log_vals[[match(p, ratio_parameters)]]
      for (col in invivo_only_oxygen_feature_cols) row[[col]] <- invivo_only_features[[i]][[col]]
      rows[[idx]] <- row
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

choose_k <- function(x, k_spec = "auto", k_max = 8L) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  n <- nrow(x)
  if (n <= 2L) return(1L)
  if (!identical(k_spec, "auto")) {
    k <- suppressWarnings(as.integer(k_spec))
    return(max(1L, min(k, n)))
  }
  max_k <- max(2L, min(k_max, n - 1L))
  d <- stats::dist(x)
  hc <- stats::hclust(d, method = "ward.D2")
  if (!requireNamespace("cluster", quietly = TRUE)) return(min(5L, max_k))
  scores <- sapply(2L:max_k, function(k) {
    cl <- stats::cutree(hc, k = k)
    mean(cluster::silhouette(cl, d)[, "sil_width"], na.rm = TRUE)
  })
  as.integer((2L:max_k)[which.max(scores)])
}

umap_embedding <- function(x, n_neighbors, umap_seed) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  if (nrow(x) <= 1L) {
    emb <- matrix(c(0, 0), nrow = 1L)
    colnames(emb) <- c("UMAP1", "UMAP2")
    return(emb)
  }
  if (nrow(x) < 3L) {
    emb <- stats::cmdscale(stats::dist(x), k = 2L)
    emb <- as.matrix(emb)
    if (ncol(emb) < 2L) emb <- cbind(emb[, 1L], 0)
    colnames(emb) <- c("UMAP1", "UMAP2")
    return(emb)
  }
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Package 'uwot' is required for UMAP.", call. = FALSE)
  set.seed(umap_seed)
  umap_args <- list(
    X = x,
    n_neighbors = min(n_neighbors, nrow(x) - 1L),
    min_dist = 0.1,
    metric = "euclidean",
    n_components = 2L,
    n_threads = 1L,
    ret_model = FALSE,
    verbose = FALSE
  )
  if ("n_sgd_threads" %in% names(formals(uwot::umap))) {
    umap_args$n_sgd_threads <- 1L
  }
  emb <- do.call(uwot::umap, umap_args)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb
}

cluster_means <- function(pairs,
                          rank_col,
                          seed_col,
                          objective_col,
                          seed_dir_col,
                          k_spec,
                          k_max,
                          extra_feature_fn = NULL,
                          extra_feature_cols = character()) {
  log_cols <- ratio_feature_cols
  keys <- unique(pairs[, c(rank_col, seed_col, objective_col, seed_dir_col), drop = FALSE])
  names(keys) <- c("rank", "seed", "objective", "seed_dir")
  means <- aggregate(pairs[, log_cols], list(rank = pairs[[rank_col]]), mean)
  means <- merge(keys, means, by = "rank", all.x = TRUE, sort = FALSE)
  means <- means[order(means$rank), , drop = FALSE]
  if (!is.null(extra_feature_fn)) {
    extra <- do.call(rbind, lapply(means$seed_dir, extra_feature_fn))
    extra <- as.data.frame(extra, check.names = FALSE, stringsAsFactors = FALSE)
    missing_extra <- setdiff(extra_feature_cols, names(extra))
    if (length(missing_extra)) {
      stop("Extra feature function did not return columns: ", paste(missing_extra, collapse = ", "), call. = FALSE)
    }
    for (col in extra_feature_cols) means[[col]] <- suppressWarnings(as.numeric(extra[[col]]))
  }
  feature_cols <- c(log_cols, extra_feature_cols)
  x <- as.matrix(means[, feature_cols, drop = FALSE])
  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  cluster_x <- x_scaled
  k <- choose_k(cluster_x, k_spec = k_spec, k_max = k_max)
  cl <- if (k <= 1L) rep(1L, nrow(means)) else stats::cutree(stats::hclust(stats::dist(cluster_x), method = "ward.D2"), k = k)
  means$cluster <- as.integer(cl)
  means$cluster_space <- "profile"
  means$cluster_feature_count <- length(feature_cols)
  means$cluster_feature_set <- paste(feature_cols, collapse = ",")
  means
}

representatives_from_clusters <- function(means, feature_cols = ratio_feature_cols) {
  missing <- setdiff(feature_cols, names(means))
  if (length(missing)) stop("Representative feature columns are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x_all <- scale(as.matrix(means[, feature_cols, drop = FALSE]))
  x_all[!is.finite(x_all)] <- 0
  means$.representative_row_id <- seq_len(nrow(means))
  reps <- do.call(rbind, lapply(split(means, means$cluster), function(df) {
    x <- x_all[df$.representative_row_id, , drop = FALSE]
    centroid <- colMeans(x)
    df$centroid_distance <- sqrt(rowSums(sweep(x, 2L, centroid, "-")^2))
    df <- df[order(df$centroid_distance, df$objective, df$rank), , drop = FALSE]
    df[1L, , drop = FALSE]
  }))
  reps$.representative_row_id <- NULL
  reps[order(reps$cluster), , drop = FALSE]
}

parse_rank_list <- function(spec, available_ranks) {
  spec <- trimws(spec)
  if (!nzchar(spec)) return(integer(0))
  vals <- suppressWarnings(as.integer(trimws(strsplit(spec, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  unique(vals[vals %in% available_ranks])
}

auto_phase2_invitro_ranks <- function(invitro_means) {
  log_cols <- ratio_feature_cols
  x <- as.matrix(invitro_means[, log_cols, drop = FALSE])
  global <- colMeans(x)
  dist_global <- sqrt(rowSums(sweep(x, 2L, global, "-")^2))
  outlier_rank <- invitro_means$rank[[which.max(dist_global)]]
  middle_candidates <- invitro_means[!invitro_means$rank %in% c(1L, outlier_rank), , drop = FALSE]
  if (nrow(middle_candidates)) {
    middle_x <- as.matrix(middle_candidates[, log_cols, drop = FALSE])
    middle_rank <- middle_candidates$rank[[which.min(sqrt(rowSums(sweep(middle_x, 2L, global, "-")^2)))]]
  } else {
    middle_rank <- integer(0)
  }
  unique(c(1L, middle_rank, outlier_rank))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

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

paired_profile_umap_coords <- function(pairs, out_dir, umap_seed) {
  feature_cols <- invivo_profile_feature_cols
  required <- c(
    "pair_id", "invivo_rank", "invitro_rank", "invivo_seed", "invitro_seed",
    "invivo_seed_dir", "invitro_seed_dir", "invivo_objective", "invitro_objective",
    feature_cols
  )
  missing <- setdiff(required, names(pairs))
  if (length(missing)) stop("Paired profile table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

  x <- as.matrix(pairs[, feature_cols, drop = FALSE])
  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  emb <- umap_embedding(x_scaled, n_neighbors = min(15L, nrow(x_scaled) - 1L), umap_seed = umap_seed)
  coords <- cbind(
    pairs[, c(
      "pair_id", "invivo_rank", "invitro_rank", "invivo_seed", "invitro_seed",
      "invivo_seed_dir", "invitro_seed_dir", "invivo_objective", "invitro_objective"
    )],
    data.frame(UMAP1 = emb[, 1L], UMAP2 = emb[, 2L], stringsAsFactors = FALSE)
  )
  write_tsv(coords, file.path(out_dir, "joint_soft_coupling_18d_profile_umap_coords.tsv"))
  coords
}

plot_paired_profile_umap <- function(coords, out_dir, top_n, umap_seed) {
  required <- c("pair_id", "invivo_rank", "invitro_rank", "UMAP1", "UMAP2")
  missing <- setdiff(required, names(coords))
  if (length(missing)) stop("18D UMAP coords are missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

  plot_df <- coords
  plot_df$invitro_rank_factor <- factor(plot_df$invitro_rank, levels = seq_len(top_n))
  limits <- square_limits(plot_df$UMAP1, plot_df$UMAP2)
  p <- ggplot(plot_df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = invivo_rank, shape = invitro_rank_factor), size = 2.8, stroke = 0.9) +
    scale_color_gradient(low = "#2b6cb0", high = "#d33f3f", name = "In Vivo Rank") +
    scale_shape_manual(values = point_shape_values(top_n), name = "In Vitro Rank") +
    labs(
      title = "18D Warm-Start Profile UMAP",
      subtitle = wrap_label(
        paste0(nrow(plot_df), " source top", top_n, " x top", top_n, " pairings; scaled 18D warm-start profile features."),
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

plot_paired_profile_umap_by_invivo_cluster <- function(coords, invivo_means, invivo_reps, out_dir, top_n, umap_seed) {
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
  plot_df$invitro_rank_factor <- factor(plot_df$invitro_rank, levels = seq_len(top_n))

  cluster_colors <- rep(
    c("#2b6cb0", "#d33f3f", "#2f855a", "#805ad5", "#dd6b20", "#319795", "#d53f8c", "#718096"),
    length.out = length(cluster_levels)
  )
  names(cluster_colors) <- sprintf("C%02d", cluster_levels)
  limits <- square_limits(plot_df$UMAP1, plot_df$UMAP2)

  p <- ggplot(plot_df, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = invivo_cluster, shape = invitro_rank_factor), size = 3.0, stroke = 0.9) +
    scale_color_manual(values = cluster_colors, name = "In Vivo Cluster") +
    scale_shape_manual(values = point_shape_values(top_n), name = "In Vitro Rank") +
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

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  invivo_run_dir <- normalizePath(as_chr(argv$invivo_run_dir), mustWork = TRUE)
  invitro_run_dir <- normalizePath(as_chr(argv$invitro_run_dir), mustWork = TRUE)
  out_dir <- normalizePath(as_chr(argv$out_dir), mustWork = FALSE)
  if (!nzchar(out_dir)) stop("--out_dir is required.", call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  clear_removed_ratio_umap_outputs(out_dir)

  top_n <- as_int(argv$top_n, 10L)
  umap_seed <- as_int(argv$umap_seed %||% argv$multi_warmup_umap_seed, 1L)
  invivo_k <- as_chr(argv$invivo_k, "auto")
  invitro_k <- as_chr(argv$invitro_k, "auto")
  invivo_cluster_space <- tolower(as_chr(argv$invivo_cluster_space, "profile"))
  if (!identical(invivo_cluster_space, "profile")) {
    stop("--invivo_cluster_space=umap has been removed; multi-warmup representative selection now always uses the 18D profile space.", call. = FALSE)
  }
  include_phase2 <- as_bool(argv$include_phase2, FALSE)
  anchor_ranks <- parse_rank_list(as_chr(argv$invitro_anchor_ranks, "1"), seq_len(top_n))

  invivo_top <- read_top_seed_table(invivo_run_dir, top_n, "in vivo")
  invitro_top <- read_top_seed_table(invitro_run_dir, top_n, "in vitro")
  pairs <- build_pair_tables(invivo_top, invitro_top)

  invivo_means <- cluster_means(
    pairs, "invivo_rank", "invivo_seed", "invivo_objective", "invivo_seed_dir",
    invivo_k,
    k_max = min(8L, top_n - 1L),
    extra_feature_fn = invivo_only_oxygen_features_for_seed,
    extra_feature_cols = invivo_only_oxygen_feature_cols
  )
  invitro_means <- cluster_means(
    pairs, "invitro_rank", "invitro_seed", "invitro_objective", "invitro_seed_dir",
    invitro_k, k_max = min(3L, top_n - 1L)
  )
  invivo_reps <- representatives_from_clusters(invivo_means, feature_cols = invivo_profile_feature_cols)
  write_tsv(invivo_means, file.path(out_dir, "multi_warmup_invivo_cluster_summary.tsv"))
  write_tsv(invitro_means, file.path(out_dir, "multi_warmup_invitro_cluster_summary.tsv"))
  write_tsv(invivo_reps, file.path(out_dir, "multi_warmup_invivo_representatives.tsv"))
  write_tsv(
    data.frame(
      feature = invivo_profile_feature_cols,
      source = c(rep("mean_paired_log10_invivo_over_invitro_ratio", length(ratio_feature_cols)), rep("invivo_only_log10_best_param", length(invivo_only_oxygen_feature_cols))),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "multi_warmup_invivo_feature_columns.tsv")
  )
  paired_coords <- paired_profile_umap_coords(pairs, out_dir = out_dir, umap_seed = umap_seed)
  plot_paired_profile_umap(paired_coords, out_dir = out_dir, top_n = top_n, umap_seed = umap_seed)
  plot_paired_profile_umap_by_invivo_cluster(paired_coords, invivo_means, invivo_reps, out_dir = out_dir, top_n = top_n, umap_seed = umap_seed)

  if (!length(anchor_ranks)) anchor_ranks <- 1L
  if (include_phase2) {
    phase2_spec <- as_chr(argv$phase2_invitro_anchor_ranks, "auto")
    phase2_ranks <- if (identical(phase2_spec, "auto")) auto_phase2_invitro_ranks(invitro_means) else parse_rank_list(phase2_spec, seq_len(top_n))
    anchor_ranks <- unique(c(anchor_ranks, phase2_ranks))
  }
  anchors <- invitro_means[invitro_means$rank %in% anchor_ranks, , drop = FALSE]
  anchors <- anchors[order(match(anchors$rank, anchor_ranks)), , drop = FALSE]

  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(invivo_reps))) {
    for (j in seq_len(nrow(anchors))) {
      warmup_label <- sprintf("v%02d_i%02d", invivo_reps$rank[[i]], anchors$rank[[j]])
      rows[[idx]] <- data.frame(
        warmup_label = warmup_label,
        phase = if (anchors$rank[[j]] == 1L) "phase1" else "phase2",
        invivo_family = sprintf("C%02d", invivo_reps$cluster[[i]]),
        invivo_rank = invivo_reps$rank[[i]],
        invivo_seed = invivo_reps$seed[[i]],
        invivo_seed_dir = invivo_reps$seed_dir[[i]],
        invitro_family = sprintf("I%02d", anchors$cluster[[j]]),
        invitro_rank = anchors$rank[[j]],
        invitro_seed = anchors$seed[[j]],
        invitro_seed_dir = anchors$seed_dir[[j]],
        selection_reason = "invivo_nearest_18d_profile_cluster_centroid",
        joint_run_prefix = paste0("fit_joint_", warmup_label),
        joint_soft_coupling_parameters_table = file.path(out_dir, "joint_soft_coupling_tables", paste0("joint_soft_coupling_parameters_table__", warmup_label, ".csv")),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  manifest <- do.call(rbind, rows)
  write_tsv(manifest, file.path(out_dir, "multi_warmup_manifest.tsv"))
  message("Wrote multi-warmup manifest: ", file.path(out_dir, "multi_warmup_manifest.tsv"))
  message("Selected ", length(unique(manifest$invivo_family)), " in vivo clusters and ", nrow(manifest), " warm-up pairs.")
}

if (identical(environment(), globalenv())) {
  main()
}
