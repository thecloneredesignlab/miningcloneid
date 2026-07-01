#!/usr/bin/env Rscript

runner_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

source(file.path(runner_dir, "fixo2_eigen_attractor_feature_builder.R"))

fixo2ea_wide_path <- function(result_root, dataset, point_type) {
  file.path(
    result_root,
    "FixO2EigenAttractors", "Tables",
    paste0(normalize_dataset(dataset), "_", point_type, "_fixo2_eigen_attractor_wide.csv")
  )
}

fixo2ea_read_wide <- function(result_root, dataset, point_type) {
  path <- fixo2ea_wide_path(result_root, dataset, point_type)
  if (!file.exists(path)) stop("Missing FixO2 eigen wide table: ", path, call. = FALSE)
  df <- read_csv_plain(path)
  if (!"seed" %in% names(df)) stop("Wide table is missing seed column: ", path, call. = FALSE)
  df$seed <- as.integer(df$seed)
  if (!"dataset" %in% names(df)) df$dataset <- normalize_dataset(dataset)
  if (!"point_type" %in% names(df)) df$point_type <- point_type
  if (!"source_group" %in% names(df)) df$source_group <- paste0(normalize_dataset(dataset), "_", point_type)
  if (!"objective" %in% names(df)) df$objective <- NA_real_
  df
}

fixo2ea_feature_columns_from_tables <- function(tables) {
  cols <- lapply(tables, function(df) grep("^fixo2_eigen_ploidy_o2_", names(df), value = TRUE))
  shared <- Reduce(intersect, cols)
  if (!length(shared)) stop("No shared FixO2 eigen ploidy feature columns found.", call. = FALSE)
  shared
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

fixo2ea_metadata <- function(df) {
  keep <- intersect(
    c(
      "dataset", "point_type", "source_group", "point_id", "seed",
      "source_row_id", "initial_index", "objective", "objective_ploidy",
      "objective_burden", "objective_growth", "objective_flow",
      "n_o2", "n_o2_ok", "n_o2_failed", "min_dominant_mean_ploidy",
      "max_dominant_mean_ploidy", "min_spectral_gap", "median_spectral_gap"
    ),
    names(df)
  )
  out <- df[, keep, drop = FALSE]
  if (!"objective" %in% names(out)) out$objective <- NA_real_
  out$objective <- suppressWarnings(as.numeric(out$objective))
  out
}

fixo2ea_run_reduction_block <- function(feature_df,
                                        initial_meta,
                                        best_meta,
                                        dataset_name,
                                        result_root,
                                        variant,
                                        output_stem,
                                        reductions = c("pca", "umap", "tsne"),
                                        preprocess_mode = "zscore",
                                        pooled = FALSE,
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
                                        initial_size = 0.22,
                                        best_size = 1.25) {
  reductions <- unique(vapply(reductions, normalize_reduction, character(1L)))
  prepared <- prepare_feature_matrix(feature_df, preprocess_mode = preprocess_mode, pooled = pooled)
  for (reduction in reductions) {
    dirs <- fixo2ea_reduction_dirs(result_root, dataset_name, reduction, variant)
    for (path in unlist(dirs, use.names = FALSE)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    out_prefix <- paste0(output_stem, "_", reduction_file_suffix(reduction))
    figure_prefix <- file.path(dirs$figures, out_prefix)
    table_prefix <- file.path(dirs$tables, out_prefix)
    emb <- run_reduction_embedding(
      prepared$mat,
      reduction = reduction,
      label = paste(dataset_name, variant, "FixO2 eigen", reduction),
      table_prefix = table_prefix,
      figure_prefix = figure_prefix,
      preprocess_mode = preprocess_mode,
      umap_seed = umap_seed,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      n_threads = n_threads,
      tsne_seed = tsne_seed,
      tsne_perplexity = tsne_perplexity,
      tsne_theta = tsne_theta,
      tsne_max_iter = tsne_max_iter
    )
    if (isTRUE(pooled)) {
      write_pooled_reduction_outputs(
        reduction = reduction,
        emb = emb,
        initial_df = initial_meta,
        best_df = best_meta,
        reduction_label = paste(dataset_name, variant, reduction, sep = "_"),
        feature_mat = prepared$mat,
        feature_metadata = prepared$metadata,
        prior_metadata = NULL,
        figure_prefix = figure_prefix,
        table_prefix = table_prefix,
        distance_table_path = file.path(dirs$tables, paste0(out_prefix, "_best_pair_distances.csv")),
        initial_size = initial_size,
        best_size = best_size,
        figures_wclusters_dir = if (isTRUE(run_clustered)) dirs$figures_wc else NULL,
        tables_wclusters_dir = if (isTRUE(run_clustered)) dirs$tables_wc else NULL,
        cluster_seed = cluster_seed,
        cluster_k_min = cluster_k_min,
        cluster_k_max = cluster_k_max,
        cluster_silhouette_sample_n = cluster_silhouette_sample_n
      )
    } else {
      write_reduction_outputs(
        reduction = reduction,
        emb = emb,
        initial_df = initial_meta,
        best_df = best_meta,
        reduction_label = paste(dataset_name, variant, reduction, sep = "_"),
        feature_mat = prepared$mat,
        feature_metadata = prepared$metadata,
        prior_metadata = NULL,
        figure_prefix = figure_prefix,
        table_prefix = table_prefix,
        initial_size = initial_size,
        best_size = best_size,
        shape_by_pred = FALSE,
        figures_wclusters_dir = if (isTRUE(run_clustered)) dirs$figures_wc else NULL,
        tables_wclusters_dir = if (isTRUE(run_clustered)) dirs$tables_wc else NULL,
        cluster_seed = cluster_seed,
        cluster_k_min = cluster_k_min,
        cluster_k_max = cluster_k_max,
        cluster_silhouette_sample_n = cluster_silhouette_sample_n
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
  feature_cols <- fixo2ea_feature_columns_from_tables(list(best, initial))
  best_meta <- fixo2ea_metadata(best)
  initial_meta <- fixo2ea_metadata(initial)
  if (isTRUE(run_initial_vs_best)) {
    feature_df <- rbind(initial[, feature_cols, drop = FALSE], best[, feature_cols, drop = FALSE])
    fixo2ea_run_reduction_block(
      feature_df = feature_df,
      initial_meta = initial_meta,
      best_meta = best_meta,
      dataset_name = dataset,
      result_root = result_root,
      variant = "Full",
      output_stem = paste0(dataset, "_full_fixo2_eigen_attractor_initial_vs_best"),
      reductions = reductions,
      pooled = FALSE,
      ...
    )
  }
  if (isTRUE(run_best_only)) {
    empty_initial <- initial_meta[0L, , drop = FALSE]
    fixo2ea_run_reduction_block(
      feature_df = best[, feature_cols, drop = FALSE],
      initial_meta = empty_initial,
      best_meta = best_meta,
      dataset_name = dataset,
      result_root = result_root,
      variant = "BestOnly",
      output_stem = paste0(dataset, "_best_fixo2_eigen_attractor"),
      reductions = reductions,
      pooled = FALSE,
      ...
    )
  }
}

fixo2ea_run_pooled_reductions <- function(result_root = fixo2ea_default_result_root(),
                                          reductions = c("pca", "umap", "tsne"),
                                          run_initial_vs_best = TRUE,
                                          run_best_only = TRUE,
                                          ...) {
  invivo_best <- fixo2ea_read_wide(result_root, "invivo", "best")
  invitro_best <- fixo2ea_read_wide(result_root, "invitro", "best")
  invivo_initial <- fixo2ea_read_wide(result_root, "invivo", "initial")
  invitro_initial <- fixo2ea_read_wide(result_root, "invitro", "initial")
  feature_cols <- fixo2ea_feature_columns_from_tables(list(invivo_best, invitro_best, invivo_initial, invitro_initial))
  best <- rbind_fill_plain(list(invivo_best, invitro_best))
  initial <- rbind_fill_plain(list(invivo_initial, invitro_initial))
  best_meta <- fixo2ea_metadata(best)
  initial_meta <- fixo2ea_metadata(initial)
  if (isTRUE(run_initial_vs_best)) {
    feature_df <- rbind(initial[, feature_cols, drop = FALSE], best[, feature_cols, drop = FALSE])
    fixo2ea_run_reduction_block(
      feature_df = feature_df,
      initial_meta = initial_meta,
      best_meta = best_meta,
      dataset_name = "pooled_invivo_invitro",
      result_root = result_root,
      variant = "Full",
      output_stem = "pooled_invivo_invitro_full_fixo2_eigen_attractor_initial_vs_best",
      reductions = reductions,
      pooled = TRUE,
      ...
    )
  }
  if (isTRUE(run_best_only)) {
    empty_initial <- initial_meta[0L, , drop = FALSE]
    fixo2ea_run_reduction_block(
      feature_df = best[, feature_cols, drop = FALSE],
      initial_meta = empty_initial,
      best_meta = best_meta,
      dataset_name = "pooled_invivo_invitro",
      result_root = result_root,
      variant = "BestOnly",
      output_stem = "pooled_invivo_invitro_best_fixo2_eigen_attractor",
      reductions = reductions,
      pooled = TRUE,
      ...
    )
  }
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
    for (dataset in datasets) {
      fixo2ea_run_dataset_reductions(
        dataset = dataset,
        result_root = result_root,
        reductions = reductions,
        run_initial_vs_best = run_initial_vs_best,
        run_best_only = run_best_only,
        ...
      )
    }
  }
  if (isTRUE(run_pooled)) {
    fixo2ea_run_pooled_reductions(
      result_root = result_root,
      reductions = reductions,
      run_initial_vs_best = run_initial_vs_best,
      run_best_only = run_best_only,
      ...
    )
  }
  invisible(TRUE)
}

fixo2ea_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  result_root <- normalizePath(path.expand(argv$result_root %||% fixo2ea_default_result_root()), mustWork = FALSE)
  source_root <- normalizePath(path.expand(argv$source_root %||% fixo2ea_default_parameter_source_root()), mustWork = FALSE)
  run_parts <- tolower(as_char_vec(argv$run_parts, c("features", "reductions", "report")))
  datasets <- as_char_vec(argv$datasets, c("invivo", "invitro"))
  reductions <- as_char_vec(argv$reductions, c("pca", "umap", "tsne"))
  o2_values <- fixo2ea_o2_grid(
    o2_values = argv$o2_values,
    o2_min = argv$o2_min %||% 0,
    o2_max = argv$o2_max %||% 5,
    o2_n = argv$o2_n %||% 201
  )
  if ("features" %in% run_parts) {
    fixo2ea_build_feature_tables(
      result_root = result_root,
      source_root = source_root,
      invivo_input = argv$invivo_input %||% default_dataset_input_dir("invivo"),
      invitro_input = argv$invitro_input %||% default_dataset_input_dir("invitro"),
      invivo_best_csv = argv$invivo_best_csv %||% fixo2ea_default_best_csv("invivo", source_root),
      invitro_best_csv = argv$invitro_best_csv %||% fixo2ea_default_best_csv("invitro", source_root),
      invivo_initial_csv = argv$invivo_initial_csv %||% fixo2ea_default_initial_csv("invivo", source_root),
      invitro_initial_csv = argv$invitro_initial_csv %||% fixo2ea_default_initial_csv("invitro", source_root),
      o2_values = o2_values,
      datasets = datasets,
      point_types = as_char_vec(argv$point_types, c("best", "initial")),
      seeds = argv$seeds,
      max_seeds = as_int(argv$max_seeds, NA_integer_),
      n_workers = as_int(argv$n_workers, 1L),
      chunk_size = as_int(argv$chunk_size, 250L),
      write_initial_long = as_bool(argv$write_initial_long, FALSE),
      force = as_bool(argv$force, FALSE)
    )
  }
  if ("reductions" %in% run_parts) {
    fixo2ea_run_all_reductions(
      result_root = result_root,
      datasets = datasets,
      reductions = reductions,
      run_single = as_bool(argv$run_single, TRUE),
      run_pooled = as_bool(argv$run_pooled, TRUE),
      run_initial_vs_best = as_bool(argv$run_initial_vs_best, TRUE),
      run_best_only = as_bool(argv$run_best_only, TRUE),
      run_clustered = as_bool(argv$run_clustered, TRUE),
      umap_seed = as_int(argv$umap_seed, 123L),
      n_neighbors = as_int(argv$n_neighbors, 80L),
      min_dist = as_num(argv$min_dist, 0.1),
      n_threads = as_int(argv$n_threads, 1L),
      tsne_seed = as_int(argv$tsne_seed, 123L),
      tsne_perplexity = as_num(argv$tsne_perplexity, 30),
      tsne_theta = as_num(argv$tsne_theta, 0.5),
      tsne_max_iter = as_int(argv$tsne_max_iter, 1000L),
      cluster_seed = as_int(argv$cluster_seed, 123L),
      cluster_k_min = as_int(argv$cluster_k_min, 2L),
      cluster_k_max = as_int(argv$cluster_k_max, 8L),
      cluster_silhouette_sample_n = as_int(argv$cluster_silhouette_sample_n, 5000L)
    )
  }
  if ("report" %in% run_parts) {
    report_script <- file.path(runner_dir, "fixo2_eigen_attractor_report.R")
    if (file.exists(report_script)) {
      source(report_script)
      fixo2ea_write_report(result_root = result_root)
    }
  }
  invisible(TRUE)
}

fixo2ea_main()
