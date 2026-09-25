#!/usr/bin/env Rscript

# Table-only fixed-O2 eigen-attractor embedding and clustering analysis.
# Input contract: materialized *_fixo2_eigen_attractor_wide.csv tables.
# Output contract: coordinate, preprocessing, variance, cluster, and silhouette
# tables.  Figure rendering is owned by vis/fixed_o2_eigen/.

.fixo2ea_analysis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "analyze_fixo2_eigen_attractor_embeddings.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  workflow_frames <- frames[
    grepl("/O2_supply_demand_MAP/", frames, fixed = TRUE)
  ]
  if (length(workflow_frames)) {
    root <- dirname(workflow_frames[[length(workflow_frames)]])
    while (!identical(basename(root), "O2_supply_demand_MAP")) {
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    if (identical(basename(root), "O2_supply_demand_MAP")) {
      return(file.path(root, "analysis", "fixed_o2_eigen"))
    }
  }
  normalizePath(getwd(), mustWork = FALSE)
})
.fixo2ea_workflow_root <- normalizePath(file.path(.fixo2ea_analysis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.fixo2ea_workflow_root, "util", "o2_supply_demand_map_eigen_attractor_utils.R"), local = environment())

fixo2ea_wide_path <- function(result_root, dataset, point_type) {
  file.path(result_root, "FixO2EigenAttractors", "Tables",
            paste0(normalize_dataset(dataset), "_", point_type, "_fixo2_eigen_attractor_wide.csv"))
}

fixo2ea_read_wide <- function(result_root, dataset, point_type) {
  path <- fixo2ea_wide_path(result_root, dataset, point_type)
  if (!file.exists(path)) stop("Missing FixO2 eigen wide table: ", path, call. = FALSE)
  df <- read_csv_plain(path)
  if (!"seed" %in% names(df)) stop("Wide table is missing seed column: ", path, call. = FALSE)
  df$seed <- suppressWarnings(as.integer(df$seed))
  if (!"dataset" %in% names(df)) df$dataset <- normalize_dataset(dataset)
  if (!"point_type" %in% names(df)) df$point_type <- point_type
  if (!"source_group" %in% names(df)) df$source_group <- paste0(normalize_dataset(dataset), "_", point_type)
  if (!"objective" %in% names(df)) df$objective <- NA_real_
  df
}

fixo2ea_feature_columns_from_tables <- function(tables) {
  columns <- lapply(tables, function(df) grep("^fixo2_eigen_ploidy_o2_", names(df), value = TRUE))
  shared <- Reduce(intersect, columns)
  if (!length(shared)) stop("No shared FixO2 eigen ploidy feature columns found.", call. = FALSE)
  shared
}

fixo2ea_metadata <- function(df) {
  keep <- intersect(c(
    "dataset", "point_type", "source_group", "point_id", "seed", "source_row_id",
    "initial_index", "objective", "objective_ploidy", "objective_burden",
    "objective_growth", "objective_flow", "n_o2", "n_o2_ok", "n_o2_failed",
    "min_dominant_mean_ploidy", "max_dominant_mean_ploidy", "min_spectral_gap",
    "median_spectral_gap"
  ), names(df))
  out <- df[, keep, drop = FALSE]
  if (!"objective" %in% names(out)) out$objective <- NA_real_
  out$objective <- suppressWarnings(as.numeric(out$objective))
  out
}

fixo2ea_scale_features <- function(feature_df) {
  matrix <- as.matrix(data.frame(lapply(feature_df, as.numeric), check.names = FALSE))
  if (!nrow(matrix) || !ncol(matrix)) stop("Empty FixO2 eigen feature matrix.", call. = FALSE)
  for (column in seq_len(ncol(matrix))) {
    finite <- is.finite(matrix[, column])
    replacement <- if (any(finite)) stats::median(matrix[finite, column]) else 0
    matrix[!finite, column] <- replacement
  }
  center <- colMeans(matrix)
  scale <- apply(matrix, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  scaled <- sweep(sweep(matrix, 2L, center, "-"), 2L, scale, "/")
  list(
    matrix = scaled,
    metadata = data.frame(feature = colnames(matrix), center = center, scale = scale,
                          preprocessing = "zscore", stringsAsFactors = FALSE)
  )
}

fixo2ea_run_embedding <- function(matrix,
                                  reduction,
                                  umap_seed = 123L,
                                  n_neighbors = 80L,
                                  min_dist = 0.1,
                                  n_threads = 1L,
                                  tsne_seed = 123L,
                                  tsne_perplexity = 30,
                                  tsne_theta = 0.5,
                                  tsne_max_iter = 1000L) {
  reduction <- normalize_reduction(reduction)
  coord_names <- reduction_coordinate_names(reduction)
  if (reduction == "pca") {
    fit <- stats::prcomp(matrix, center = FALSE, scale. = FALSE)
    coordinates <- fit$x[, seq_len(min(2L, ncol(fit$x))), drop = FALSE]
    if (ncol(coordinates) == 1L) coordinates <- cbind(coordinates, 0)
    colnames(coordinates) <- coord_names
    variance <- fit$sdev^2
    variance_table <- data.frame(
      PC = paste0("PC", seq_along(variance)),
      variance = variance,
      variance_fraction = variance / sum(variance),
      cumulative_variance_fraction = cumsum(variance) / sum(variance),
      stringsAsFactors = FALSE
    )
    return(list(coordinates = coordinates, variance = variance_table))
  }
  if (nrow(matrix) < 4L) stop(reduction, " requires at least four rows.", call. = FALSE)
  if (reduction == "umap") {
    if (!requireNamespace("uwot", quietly = TRUE)) stop("Package 'uwot' is required for UMAP.", call. = FALSE)
    set.seed(as.integer(umap_seed))
    coordinates <- uwot::umap(
      matrix,
      n_neighbors = max(2L, min(as.integer(n_neighbors), nrow(matrix) - 1L)),
      min_dist = as.numeric(min_dist),
      n_threads = max(1L, as.integer(n_threads)),
      ret_model = FALSE,
      verbose = FALSE
    )
    colnames(coordinates) <- coord_names
    return(list(coordinates = coordinates, variance = NULL))
  }
  if (!requireNamespace("Rtsne", quietly = TRUE)) stop("Package 'Rtsne' is required for t-SNE.", call. = FALSE)
  set.seed(as.integer(tsne_seed))
  perplexity <- min(as.numeric(tsne_perplexity), max(1, floor((nrow(matrix) - 1) / 3)))
  fit <- Rtsne::Rtsne(matrix, dims = 2L, perplexity = perplexity,
                      theta = as.numeric(tsne_theta), max_iter = as.integer(tsne_max_iter),
                      check_duplicates = FALSE, pca = FALSE, verbose = FALSE)
  coordinates <- fit$Y
  colnames(coordinates) <- coord_names
  list(coordinates = coordinates, variance = NULL)
}

fixo2ea_silhouette_mean <- function(matrix, cluster) {
  if (length(unique(cluster)) < 2L || !requireNamespace("cluster", quietly = TRUE)) return(NA_real_)
  mean(cluster::silhouette(cluster, stats::dist(matrix))[, "sil_width"])
}

fixo2ea_cluster_assignments <- function(matrix,
                                        seed = 123L,
                                        k_min = 2L,
                                        k_max = 8L,
                                        sample_n = 5000L) {
  matrix <- as.matrix(matrix)
  n <- nrow(matrix)
  if (n < 3L) {
    return(list(cluster = rep.int(1L, n), selected_k = 1L,
                silhouette = data.frame(k = 1L, average_silhouette = NA_real_, sample_n = n)))
  }
  set.seed(as.integer(seed))
  selected <- if (n > sample_n) sort(sample.int(n, sample_n)) else seq_len(n)
  basis <- matrix[selected, , drop = FALSE]
  candidates <- seq.int(max(2L, as.integer(k_min)), min(as.integer(k_max), nrow(basis) - 1L))
  scores <- lapply(candidates, function(k) {
    fit <- stats::kmeans(basis, centers = k, nstart = 20L, iter.max = 100L)
    data.frame(k = k, average_silhouette = fixo2ea_silhouette_mean(basis, fit$cluster),
               sample_n = nrow(basis), stringsAsFactors = FALSE)
  })
  score_table <- do.call(rbind, scores)
  rank_score <- ifelse(is.finite(score_table$average_silhouette), score_table$average_silhouette, -Inf)
  selected_k <- score_table$k[[which.max(rank_score)]]
  set.seed(as.integer(seed))
  fit <- stats::kmeans(matrix, centers = selected_k, nstart = 30L, iter.max = 100L)
  list(cluster = fit$cluster, selected_k = selected_k, silhouette = score_table)
}

fixo2ea_reduction_dirs <- function(root, dataset_name, reduction, variant) {
  base <- file.path(root, dataset_name, reduction_dir_name(reduction))
  list(
    tables = file.path(base, "Tables", variant),
    figures = file.path(base, "Figures", variant),
    tables_wc = file.path(base, "TablesWclusters", variant),
    figures_wc = file.path(base, "FiguresWclusters", variant)
  )
}

fixo2ea_coordinate_cluster_label <- function(reduction) {
  paste0("By", switch(normalize_reduction(reduction), pca = "PCA", umap = "UMAP", tsne = "TSNE"), "Coordinates")
}

fixo2ea_write_cluster_tables <- function(coordinates,
                                         feature_matrix,
                                         table_prefix,
                                         tables_wclusters_dir,
                                         reduction,
                                         cluster_seed,
                                         cluster_k_min,
                                         cluster_k_max,
                                         cluster_silhouette_sample_n) {
  methods <- list(ByInputFeatures = feature_matrix)
  methods[[fixo2ea_coordinate_cluster_label(reduction)]] <- as.matrix(coordinates[, reduction_coordinate_names(reduction), drop = FALSE])
  for (method in names(methods)) {
    result <- fixo2ea_cluster_assignments(
      methods[[method]], seed = cluster_seed,
      k_min = cluster_k_min, k_max = cluster_k_max,
      sample_n = cluster_silhouette_sample_n
    )
    clustered <- coordinates
    clustered$cluster_num <- result$cluster
    clustered$cluster_id <- paste0("C", result$cluster)
    clustered$cluster_source <- method
    clustered$cluster_k <- result$selected_k
    out_dir <- file.path(tables_wclusters_dir, method)
    write_csv(clustered, file.path(out_dir, paste0(basename(table_prefix), "_coordinates.csv")))
    write_csv(result$silhouette, file.path(out_dir, paste0(basename(table_prefix), "_silhouette.csv")))
  }
}

fixo2ea_run_reduction_block <- function(feature_df,
                                        initial_meta,
                                        best_meta,
                                        dataset_name,
                                        result_root,
                                        variant,
                                        output_stem,
                                        reductions = c("pca", "umap", "tsne"),
                                        run_clustered = TRUE,
                                        umap_seed = 123L,
                                        n_neighbors = 80L,
                                        min_dist = 0.1,
                                        n_threads = 1L,
                                        tsne_seed = 123L,
                                        tsne_perplexity = 30,
                                        tsne_theta = 0.5,
                                        tsne_max_iter = 1000L,
                                        cluster_seed = 123L,
                                        cluster_k_min = 2L,
                                        cluster_k_max = 8L,
                                        cluster_silhouette_sample_n = 5000L,
                                        ...) {
  prepared <- fixo2ea_scale_features(feature_df)
  metadata <- rbind_fill_plain(list(initial_meta, best_meta))
  if (nrow(metadata) != nrow(prepared$matrix)) stop("Feature and metadata row counts differ.", call. = FALSE)
  for (reduction in unique(vapply(reductions, normalize_reduction, character(1L)))) {
    dirs <- fixo2ea_reduction_dirs(result_root, dataset_name, reduction, variant)
    dir.create(dirs$tables, recursive = TRUE, showWarnings = FALSE)
    embedding <- fixo2ea_run_embedding(
      prepared$matrix, reduction,
      umap_seed = umap_seed, n_neighbors = n_neighbors, min_dist = min_dist, n_threads = n_threads,
      tsne_seed = tsne_seed, tsne_perplexity = tsne_perplexity,
      tsne_theta = tsne_theta, tsne_max_iter = tsne_max_iter
    )
    coordinate_table <- cbind(metadata, as.data.frame(embedding$coordinates, check.names = FALSE))
    stem <- paste0(output_stem, "_", reduction)
    table_prefix <- file.path(dirs$tables, stem)
    write_csv(coordinate_table, paste0(table_prefix, "_coordinates.csv"))
    write_csv(prepared$metadata, paste0(table_prefix, "_preprocessing_metadata.csv"))
    if (!is.null(embedding$variance)) write_csv(embedding$variance, paste0(table_prefix, "_pca_variance.csv"))
    if (isTRUE(run_clustered)) {
      fixo2ea_write_cluster_tables(
        coordinate_table, prepared$matrix, table_prefix, dirs$tables_wc, reduction,
        cluster_seed, cluster_k_min, cluster_k_max, cluster_silhouette_sample_n
      )
    }
  }
  invisible(TRUE)
}

fixo2ea_run_dataset_reductions <- function(dataset,
                                           result_root = fixo2ea_default_result_root(),
                                           reductions = c("pca", "umap", "tsne"),
                                           run_initial_vs_best = TRUE,
                                           run_best_only = TRUE,
                                           ...) {
  dataset <- normalize_dataset(dataset)
  best <- fixo2ea_read_wide(result_root, dataset, "best")
  initial <- fixo2ea_read_wide(result_root, dataset, "initial")
  features <- fixo2ea_feature_columns_from_tables(list(best, initial))
  if (isTRUE(run_initial_vs_best)) {
    fixo2ea_run_reduction_block(
      rbind(initial[, features, drop = FALSE], best[, features, drop = FALSE]),
      fixo2ea_metadata(initial), fixo2ea_metadata(best), dataset, result_root, "Full",
      paste0(dataset, "_full_fixo2_eigen_attractor_initial_vs_best"), reductions, ...
    )
  }
  if (isTRUE(run_best_only)) {
    fixo2ea_run_reduction_block(
      best[, features, drop = FALSE], fixo2ea_metadata(initial)[0L, , drop = FALSE],
      fixo2ea_metadata(best), dataset, result_root, "BestOnly",
      paste0(dataset, "_best_fixo2_eigen_attractor"), reductions, ...
    )
  }
  invisible(TRUE)
}

fixo2ea_run_pooled_reductions <- function(result_root = fixo2ea_default_result_root(),
                                          reductions = c("pca", "umap", "tsne"),
                                          run_initial_vs_best = TRUE,
                                          run_best_only = TRUE,
                                          ...) {
  best <- rbind_fill_plain(list(
    fixo2ea_read_wide(result_root, "invivo", "best"),
    fixo2ea_read_wide(result_root, "invitro", "best")
  ))
  initial <- rbind_fill_plain(list(
    fixo2ea_read_wide(result_root, "invivo", "initial"),
    fixo2ea_read_wide(result_root, "invitro", "initial")
  ))
  features <- fixo2ea_feature_columns_from_tables(list(best, initial))
  if (isTRUE(run_initial_vs_best)) {
    fixo2ea_run_reduction_block(
      rbind(initial[, features, drop = FALSE], best[, features, drop = FALSE]),
      fixo2ea_metadata(initial), fixo2ea_metadata(best), "pooled_invivo_invitro", result_root, "Full",
      "pooled_invivo_invitro_full_fixo2_eigen_attractor_initial_vs_best", reductions, ...
    )
  }
  if (isTRUE(run_best_only)) {
    fixo2ea_run_reduction_block(
      best[, features, drop = FALSE], fixo2ea_metadata(initial)[0L, , drop = FALSE],
      fixo2ea_metadata(best), "pooled_invivo_invitro", result_root, "BestOnly",
      "pooled_invivo_invitro_best_fixo2_eigen_attractor", reductions, ...
    )
  }
  invisible(TRUE)
}

fixo2ea_run_all_reductions <- function(result_root = fixo2ea_default_result_root(),
                                       datasets = c("invivo", "invitro"),
                                       reductions = c("pca", "umap", "tsne"),
                                       run_single = TRUE,
                                       run_pooled = TRUE,
                                       run_initial_vs_best = TRUE,
                                       run_best_only = TRUE,
                                       ...) {
  if (isTRUE(run_single)) {
    for (dataset in datasets) fixo2ea_run_dataset_reductions(
      dataset, result_root, reductions, run_initial_vs_best, run_best_only, ...
    )
  }
  if (isTRUE(run_pooled)) fixo2ea_run_pooled_reductions(
    result_root, reductions, run_initial_vs_best, run_best_only, ...
  )
  invisible(TRUE)
}
