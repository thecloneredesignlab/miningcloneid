#!/usr/bin/env Rscript

.o2pl_analysis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "parameter_landscape_analysis_utils.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2pl_analysis_workflow_root <- normalizePath(file.path(.o2pl_analysis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_analysis_workflow_root, "util", "o2_supply_demand_map_parameter_landscape_io_utils.R"), local = environment(), chdir = TRUE)

o2pl_parameter_set <- umap_parameter_set
o2pl_log10_parameter_set <- umap_log10_parameter_set

# Legacy warm-up API adapters retained for callers that still use the former
# old function signatures while delegating all numerical work to this layer.
transform_umap_features <- function(df, params, log10_params) {
  o2pl_require_columns(df, params, "materialized parameter feature table")
  out <- as.data.frame(lapply(params, function(parameter) {
    value <- suppressWarnings(as.numeric(df[[parameter]]))
    if (parameter %in% log10_params) {
      if (any(!is.finite(value) | value <= 0)) stop("Cannot log10-transform parameter: ", parameter, call. = FALSE)
      value <- log10(value)
    }
    value
  }), check.names = FALSE)
  names(out) <- params
  out
}

pooled_umap_parameter_set <- function() intersect(umap_parameter_set("invivo"), umap_parameter_set("invitro"))
pooled_umap_log10_parameter_set <- function() intersect(
  pooled_umap_parameter_set(),
  union(umap_log10_parameter_set("invivo"), umap_log10_parameter_set("invitro"))
)

prepare_feature_matrix <- function(feature_df, preprocess_mode = "zscore", pooled = FALSE) {
  result <- o2pl_scale_features(feature_df, preprocess_mode)
  list(mat = result$matrix, metadata = result$metadata)
}

run_pca_embedding <- function(feature_mat,
                              label = "features",
                              variance_path = NULL,
                              variance_figure_prefix = NULL,
                              center = FALSE) {
  matrix <- as.matrix(feature_mat)
  fit <- stats::prcomp(matrix, center = isTRUE(center), scale. = FALSE)
  variance <- fit$sdev^2
  if (!is.null(variance_path)) {
    write_csv(data.frame(
      PC = paste0("PC", seq_along(variance)),
      variance = variance,
      variance_fraction = variance / sum(variance),
      cumulative_variance_fraction = cumsum(variance) / sum(variance),
      stringsAsFactors = FALSE
    ), variance_path)
  }
  if (ncol(fit$x) < 2L) stop("PCA embedding requires at least two input dimensions for ", label, call. = FALSE)
  out <- fit$x[, 1:2, drop = FALSE]
  colnames(out) <- c("PCA1", "PCA2")
  out
}

run_umap_embedding <- function(feature_mat,
                               label = "features",
                               umap_seed = 123L,
                               n_neighbors = 80L,
                               min_dist = 0.1,
                               n_threads = 1L) {
  o2pl_run_reduction(feature_mat, "umap", umap_seed, n_neighbors, min_dist, n_threads)
}

run_tsne_embedding <- function(feature_mat,
                               label = "features",
                               tsne_seed = 123L,
                               perplexity = 30,
                               theta = 0.5,
                               max_iter = 1000L,
                               num_threads = 1L) {
  o2pl_run_reduction(
    feature_mat, "tsne", tsne_seed,
    n_threads = num_threads,
    perplexity = perplexity,
    theta = theta,
    max_iter = max_iter
  )
}

o2pl_materialized_feature_path <- function(root_dir, dataset, kind = c("fitted", "initial", "growth_turnover")) {
  dataset <- o2pl_normalize_dataset(dataset)
  kind <- match.arg(kind)
  filename <- switch(kind,
    fitted = "fitted_parameter_features.tsv",
    initial = "deoptim_initial_population.tsv",
    growth_turnover = "growth_turnover_100d.tsv"
  )
  file.path(o2pl_simulation_tables_dir(root_dir, dataset), filename)
}

o2pl_read_materialized_features <- function(root_dir, dataset, kind = c("fitted", "initial", "growth_turnover")) {
  path <- o2pl_materialized_feature_path(root_dir, dataset, match.arg(kind))
  if (!file.exists(path)) stop("Missing materialized parameter-landscape input: ", path, call. = FALSE)
  read_tsv(path)
}

o2pl_transform_features <- function(df, dataset, parameters = o2pl_parameter_set(dataset)) {
  parameters <- intersect(parameters, names(df))
  if (length(parameters) < 2L) stop("Fewer than two fitted parameter columns are available.", call. = FALSE)
  log_parameters <- o2pl_log10_parameter_set(dataset)
  out <- as.data.frame(lapply(parameters, function(parameter) {
    values <- suppressWarnings(as.numeric(df[[parameter]]))
    if (parameter %in% log_parameters) {
      if (any(!is.finite(values) | values <= 0)) stop("Cannot log10-transform ", parameter, call. = FALSE)
      values <- log10(values)
    }
    values
  }), check.names = FALSE)
  names(out) <- parameters
  if (any(!is.finite(as.matrix(out)))) stop("Feature matrix contains non-finite values.", call. = FALSE)
  out
}

o2pl_scale_features <- function(features, method = "zscore") {
  method <- tolower(trimws(method %||% "zscore"))
  matrix <- as.matrix(features)
  storage.mode(matrix) <- "double"
  if (method %in% c("zscore", "standard", "standardized")) {
    center <- colMeans(matrix)
    scale <- apply(matrix, 2L, stats::sd)
    scale[!is.finite(scale) | scale <= 0] <- 1
    result <- sweep(sweep(matrix, 2L, center, `-`), 2L, scale, `/`)
    metadata <- data.frame(parameter = colnames(matrix), preprocessing = "zscore", center = center, scale = scale)
  } else if (method %in% c("prior_unit", "context_prior_unit", "common_prior_unit", "unit")) {
    lower <- apply(matrix, 2L, min)
    upper <- apply(matrix, 2L, max)
    span <- upper - lower
    span[!is.finite(span) | span <= 0] <- 1
    result <- sweep(sweep(matrix, 2L, lower, `-`), 2L, span, `/`)
    metadata <- data.frame(parameter = colnames(matrix), preprocessing = "materialized_range_unit", lower = lower, upper = upper)
  } else {
    stop("Unsupported preprocessing method: ", method, call. = FALSE)
  }
  list(matrix = result, metadata = metadata)
}

o2pl_run_reduction <- function(feature_matrix,
                               reduction = "umap",
                               seed = 123L,
                               n_neighbors = 80L,
                               min_dist = 0.1,
                               n_threads = 1L,
                               perplexity = 30,
                               theta = 0.5,
                               max_iter = 1000L) {
  reduction <- o2pl_normalize_reduction(reduction)
  feature_matrix <- as.matrix(feature_matrix)
  if (nrow(feature_matrix) < 2L) stop("Reduction requires at least two rows.", call. = FALSE)
  set.seed(as.integer(seed))
  if (reduction == "pca") {
    fit <- stats::prcomp(feature_matrix, center = FALSE, scale. = FALSE)
    coordinates <- fit$x[, seq_len(min(2L, ncol(fit$x))), drop = FALSE]
    if (ncol(coordinates) == 1L) coordinates <- cbind(coordinates, 0)
  } else if (reduction == "umap" && requireNamespace("uwot", quietly = TRUE) && nrow(feature_matrix) >= 4L) {
    coordinates <- uwot::umap(
      feature_matrix,
      n_neighbors = min(max(2L, as.integer(n_neighbors)), nrow(feature_matrix) - 1L),
      min_dist = as.numeric(min_dist),
      n_components = 2L,
      n_threads = max(1L, as.integer(n_threads)),
      seed = as.integer(seed),
      verbose = FALSE
    )
  } else if (reduction == "tsne" && requireNamespace("Rtsne", quietly = TRUE) && nrow(feature_matrix) >= 4L) {
    coordinates <- Rtsne::Rtsne(
      feature_matrix,
      dims = 2L,
      perplexity = min(as.numeric(perplexity), max(1, floor((nrow(feature_matrix) - 1L) / 3L))),
      theta = as.numeric(theta),
      max_iter = as.integer(max_iter),
      pca = FALSE,
      check_duplicates = FALSE,
      num_threads = max(1L, as.integer(n_threads)),
      verbose = FALSE
    )$Y
  } else {
    fit <- stats::prcomp(feature_matrix, center = FALSE, scale. = FALSE)
    coordinates <- fit$x[, seq_len(min(2L, ncol(fit$x))), drop = FALSE]
    if (ncol(coordinates) == 1L) coordinates <- cbind(coordinates, 0)
    attr(coordinates, "fallback") <- paste0(reduction, "_to_pca")
  }
  colnames(coordinates) <- reduction_coordinate_names(reduction)
  coordinates
}

o2pl_silhouette_score <- function(matrix, groups) {
  groups <- as.integer(groups)
  if (length(unique(groups)) < 2L || any(table(groups) < 1L)) return(NA_real_)
  distances <- as.matrix(stats::dist(matrix))
  scores <- vapply(seq_along(groups), function(i) {
    same <- which(groups == groups[[i]] & seq_along(groups) != i)
    other <- split(which(groups != groups[[i]]), groups[groups != groups[[i]]])
    a <- if (length(same)) mean(distances[i, same]) else 0
    b <- min(vapply(other, function(index) mean(distances[i, index]), numeric(1)))
    if (max(a, b) <= 0) 0 else (b - a) / max(a, b)
  }, numeric(1))
  mean(scores)
}

o2pl_cluster_coordinates <- function(coordinates, seed = 123L, k_min = 2L, k_max = 8L) {
  coordinates <- as.matrix(coordinates)
  if (nrow(coordinates) < 3L) {
    return(list(cluster = rep(1L, nrow(coordinates)), selected_k = 1L, diagnostics = data.frame(k = 1L, silhouette = NA_real_)))
  }
  k_values <- seq.int(max(2L, as.integer(k_min)), min(as.integer(k_max), nrow(coordinates) - 1L))
  diagnostics <- lapply(k_values, function(k) {
    set.seed(as.integer(seed) + k)
    fit <- stats::kmeans(coordinates, centers = k, nstart = 25L, iter.max = 200L)
    data.frame(k = k, silhouette = o2pl_silhouette_score(coordinates, fit$cluster), stringsAsFactors = FALSE)
  })
  diagnostics <- do.call(rbind, diagnostics)
  selected_k <- diagnostics$k[[which.max(ifelse(is.finite(diagnostics$silhouette), diagnostics$silhouette, -Inf))]]
  set.seed(as.integer(seed) + selected_k)
  fit <- stats::kmeans(coordinates, centers = selected_k, nstart = 50L, iter.max = 300L)
  centers <- aggregate(coordinates, list(cluster = fit$cluster), mean)
  order_index <- order(centers[[2L]], centers[[3L]])
  labels <- stats::setNames(seq_along(order_index), centers$cluster[order_index])
  list(cluster = unname(labels[as.character(fit$cluster)]), selected_k = selected_k, diagnostics = diagnostics)
}

o2pl_embedding_rows <- function(initial, fitted, coordinates, dataset, reduction, preprocessing) {
  combined <- rbind_fill_plain(list(initial, fitted))
  parameter_columns <- intersect(o2pl_parameter_set(dataset), names(combined))
  out <- data.frame(
    dataset = dataset,
    point_type = c(rep("initial", nrow(initial)), rep("fitted", nrow(fitted))),
    seed = suppressWarnings(as.integer(combined$seed)),
    objective = if ("objective" %in% names(combined)) suppressWarnings(as.numeric(combined$objective)) else NA_real_,
    reduction = reduction,
    preprocessing = preprocessing,
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(coordinates, check.names = FALSE), combined[, parameter_columns, drop = FALSE])
  row.names(out) <- NULL
  out
}

o2pl_analysis_output_prefix <- function(dataset, reduction, variant = "combined", preprocessing = "zscore") {
  suffix <- reduction_file_suffix(reduction)
  core <- if (variant == "combined") paste(dataset, "deoptim", "initial", "vs", "best", suffix, sep = "_") else paste(dataset, "best", "parameters", suffix, sep = "_")
  if (preprocessing == "zscore") core else paste(core, gsub("[^A-Za-z0-9]+", "_", preprocessing), sep = "_")
}

o2pl_materialize_embedding <- function(root_dir,
                                       dataset,
                                       reduction = "umap",
                                       preprocessing = "zscore",
                                       seed = 123L,
                                       n_neighbors = 80L,
                                       min_dist = 0.1,
                                       n_threads = 1L,
                                       cluster_seed = 123L,
                                       cluster_k_min = 2L,
                                       cluster_k_max = 8L) {
  dataset <- o2pl_normalize_dataset(dataset)
  reduction <- o2pl_normalize_reduction(reduction)
  initial <- o2pl_read_materialized_features(root_dir, dataset, "initial")
  fitted <- o2pl_read_materialized_features(root_dir, dataset, "fitted")
  parameters <- o2pl_parameter_set(dataset)
  feature_rows <- rbind(o2pl_transform_features(initial, dataset, parameters), o2pl_transform_features(fitted, dataset, parameters))
  prepared <- o2pl_scale_features(feature_rows, preprocessing)
  coordinates <- o2pl_run_reduction(prepared$matrix, reduction, seed, n_neighbors, min_dist, n_threads)
  rows <- o2pl_embedding_rows(initial, fitted, coordinates, dataset, reduction, preprocessing)
  fitted_index <- which(rows$point_type == "fitted")
  clustering <- o2pl_cluster_coordinates(coordinates[fitted_index, , drop = FALSE], cluster_seed, cluster_k_min, cluster_k_max)
  rows$cluster_id <- NA_character_
  rows$cluster_id[fitted_index] <- paste0(if (dataset == "invivo") "vi" else "vt", "_C", sprintf("%02d", clustering$cluster))
  canonical_dir <- file.path(o2pl_analysis_tables_dir(root_dir, dataset), reduction, preprocessing)
  canonical <- file.path(canonical_dir, "embedding_coordinates.tsv")
  cluster_table <- file.path(canonical_dir, "cluster_assignments.tsv")
  diagnostics <- file.path(canonical_dir, "clustering_diagnostics.tsv")
  metadata <- file.path(canonical_dir, "preprocessing_metadata.tsv")
  schema <- file.path(canonical_dir, "embedding_coordinates.schema.tsv")
  write_tsv_plain(rows, canonical)
  write_tsv_plain(rows[fitted_index, c("dataset", "seed", "objective", reduction_coordinate_names(reduction), "cluster_id"), drop = FALSE], cluster_table)
  write_tsv_plain(clustering$diagnostics, diagnostics)
  write_tsv_plain(prepared$metadata, metadata)
  o2pl_write_schema(rows, paste(dataset, reduction, preprocessing, "embedding", sep = "_"), schema)
  legacy_dir <- file.path(paper_reduction_tables_dir(dataset, reduction, root_dir), "Full")
  prefix <- o2pl_analysis_output_prefix(dataset, reduction, "combined", preprocessing)
  legacy <- file.path(legacy_dir, paste0(prefix, "_coordinates.csv"))
  legacy_rows <- rows
  legacy_rows$point_type[legacy_rows$point_type == "fitted"] <- "best"
  write_csv(legacy_rows, legacy)
  legacy_cluster_dir <- file.path(paper_reduction_tables_wclusters_dir(dataset, reduction, root_dir), reduction_coordinate_cluster_dir(reduction), "Full")
  write_csv(legacy_rows, file.path(legacy_cluster_dir, paste0(prefix, "_coordinates_clustered.csv")))
  write_csv(clustering$diagnostics, file.path(legacy_cluster_dir, paste0(prefix, "_silhouette_scores.csv")))
  o2pl_record_artifact(root_dir, "analysis", paste(dataset, reduction, preprocessing, "embedding", sep = "_"), canonical, rows, schema, "analyze_parameter_landscape.R", dataset, paste(o2pl_materialized_feature_path(root_dir, dataset, "initial"), o2pl_materialized_feature_path(root_dir, dataset, "fitted"), sep = ";"))
  invisible(list(coordinates = canonical, clusters = cluster_table, diagnostics = diagnostics, legacy_coordinates = legacy))
}

o2pl_materialize_pooled_embedding <- function(root_dir,
                                              reduction = "umap",
                                              preprocessing = "zscore",
                                              seed = 123L,
                                              n_neighbors = 80L,
                                              min_dist = 0.1,
                                              n_threads = 1L,
                                              cluster_seed = 123L,
                                              cluster_k_min = 2L,
                                              cluster_k_max = 8L) {
  tables <- lapply(c("invivo", "invitro"), function(dataset) {
    initial <- o2pl_read_materialized_features(root_dir, dataset, "initial")
    fitted <- o2pl_read_materialized_features(root_dir, dataset, "fitted")
    parameters <- o2pl_parameter_set(dataset)
    common <- intersect(o2pl_parameter_set("invivo"), o2pl_parameter_set("invitro"))
    selected <- if (length(common) >= 2L) common else parameters
    transformed <- rbind(o2pl_transform_features(initial, dataset, selected), o2pl_transform_features(fitted, dataset, selected))
    rows <- data.frame(dataset = dataset, point_type = c(rep("initial", nrow(initial)), rep("fitted", nrow(fitted))), seed = c(initial$seed, fitted$seed), objective = c(rep(NA_real_, nrow(initial)), fitted$objective), stringsAsFactors = FALSE)
    list(features = transformed, rows = rows)
  })
  common_parameters <- Reduce(intersect, lapply(tables, function(x) names(x$features)))
  features <- do.call(rbind, lapply(tables, function(x) x$features[, common_parameters, drop = FALSE]))
  rows <- do.call(rbind, lapply(tables, `[[`, "rows"))
  prepared <- o2pl_scale_features(features, preprocessing)
  coordinates <- o2pl_run_reduction(prepared$matrix, reduction, seed, n_neighbors, min_dist, n_threads)
  rows <- cbind(rows, as.data.frame(coordinates, check.names = FALSE), features)
  fitted_index <- which(rows$point_type == "fitted")
  clustering <- o2pl_cluster_coordinates(coordinates[fitted_index, , drop = FALSE], cluster_seed, cluster_k_min, cluster_k_max)
  rows$cluster_id <- NA_character_
  rows$cluster_id[fitted_index] <- paste0("pool_C", sprintf("%02d", clustering$cluster))
  canonical_dir <- file.path(o2pl_analysis_tables_dir(root_dir, "pooled"), reduction, preprocessing)
  canonical <- file.path(canonical_dir, "embedding_coordinates.tsv")
  schema <- file.path(canonical_dir, "embedding_coordinates.schema.tsv")
  write_tsv_plain(rows, canonical)
  write_tsv_plain(clustering$diagnostics, file.path(canonical_dir, "clustering_diagnostics.tsv"))
  write_tsv_plain(prepared$metadata, file.path(canonical_dir, "preprocessing_metadata.tsv"))
  o2pl_write_schema(rows, paste("pooled", reduction, preprocessing, "embedding", sep = "_"), schema)
  legacy_dir <- file.path(paper_pooled_reduction_tables_dir(reduction, root_dir), "Full")
  prefix <- paste("pooled", "invivo", "invitro", "initial", "vs", "best", reduction_file_suffix(reduction), sep = "_")
  if (preprocessing != "zscore") prefix <- paste(prefix, preprocessing, sep = "_")
  legacy <- file.path(legacy_dir, paste0(prefix, "_coordinates.csv"))
  legacy_rows <- rows
  legacy_rows$point_type[legacy_rows$point_type == "fitted"] <- "best"
  write_csv(legacy_rows, legacy)
  o2pl_record_artifact(root_dir, "analysis", paste("pooled", reduction, preprocessing, "embedding", sep = "_"), canonical, rows, schema, "analyze_parameter_landscape.R", "pooled", "materialized invivo and invitro feature tables")
  invisible(list(coordinates = canonical, legacy_coordinates = legacy))
}
