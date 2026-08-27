#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

source(file.path(local_script_dir(), "parameter_landscape_utils.R"))

argv <- parse_args(commandArgs(trailingOnly = TRUE))
analysis_part <- tolower(trimws(as.character(argv$analysis_part %||% argv$part %||% "invivo_reductions")))
analysis_part <- switch(
  analysis_part,
  invivo_table = "invivo_tables",
  invivo_tables = "invivo_tables",
  invitro_table = "invitro_tables",
  invitro_tables = "invitro_tables",
  invivo_reduction = "invivo_reductions",
  invivo_reductions = "invivo_reductions",
  invitro_reduction = "invitro_reductions",
  invitro_reductions = "invitro_reductions",
  pooled_reduction = "pooled_reductions",
  pooled_reductions = "pooled_reductions",
  stop("Unknown analysis_part: ", analysis_part, call. = FALSE)
)

if (identical(analysis_part, "invivo_tables")) {
  dataset <- "invivo"
  root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
  tables_dir <- argv$tables_dir %||% paper_tables_dir(dataset, root_dir = root_dir)
  seed_tables_dir <- argv$seed_tables_dir %||% argv$seed_parameter_tables_dir %||% paper_seed_parameter_tables_dir(root_dir = root_dir)
  best_csv <- argv$best_csv %||% paper_best_params_csv(dataset, root_dir = root_dir)
  write_best <- as_bool(argv$write_best, TRUE)
  write_initial <- as_bool(argv$write_initial, TRUE)
  write_growth_turnover <- as_bool(argv$write_growth_turnover, TRUE)
  write_ploidy_ratios <- as_bool(argv$write_ploidy_ratios, FALSE)
  write_modes <- as_bool(argv$write_modes, TRUE)

  if (!write_best && !write_initial && !write_growth_turnover && !write_ploidy_ratios && !write_modes) {
    stop("Nothing to write: enable write_best, write_initial, write_growth_turnover, write_ploidy_ratios, or write_modes.")
  }
  if (write_best || write_initial) {
    paper_generate_umap_tables(
      dataset = dataset,
      input_dir = argv$input_dir %||% default_dataset_input_dir(dataset),
      root_dir = root_dir,
      tables_dir = seed_tables_dir,
      max_seeds = as_int(argv$max_seeds, NA_integer_),
      write_best = write_best,
      write_initial = write_initial
    )
  }
  if (write_growth_turnover) {
    paper_generate_invivo_growth_turnover_table(
      input_dir = argv$input_dir %||% default_dataset_input_dir(dataset),
      best_csv = best_csv,
      output_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
      data_dir = argv$data_dir %||% default_invivo_data_dir(),
      max_seeds = as_int(argv$max_seeds, NA_integer_),
      horizon_day = as_num(argv$growth_turnover_horizon_day, 100),
      report_dt = as_num(argv$growth_turnover_report_dt, 1.0)
    )
  } else if (write_ploidy_ratios) {
    paper_update_invivo_growth_turnover_ploidy_ratios(
      input_dir = argv$input_dir %||% default_dataset_input_dir(dataset),
      growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
      target_day = as_num(argv$pred1000_target_day, 1000)
    )
  }
  if (write_modes) {
    paper_generate_invivo_fixed_o2_mode_tables(
      input_dir = argv$input_dir %||% default_dataset_input_dir(dataset),
      mode_tables_dir = argv$mode_tables_dir %||% paper_fixo2_mode_tables_dir(root_dir = root_dir),
      best_csv = best_csv,
      augmented_best_csv = argv$augmented_best_csv %||% file.path(tables_dir, "invivo_best_params_by_seed_with_fixed_o2_mode.csv"),
      attractor_o2_grid = as_num_vec(argv$attractor_o2_grid, paper_default_attractor_o2_grid()),
      summary_o2 = as_num_vec(argv$summary_o2, paper_default_mode_summary_o2()),
      mode_reference_o2 = as_num(argv$mode_reference_o2, 2),
      mode_reference_o2_values = argv$mode_reference_o2_values,
      max_seeds = as_int(argv$max_seeds, NA_integer_),
      n_workers = as_int(argv$n_workers, 1L),
      write_augmented_best = as_bool(argv$write_augmented_best, TRUE),
      overwrite_modes = as_bool(argv$overwrite_modes %||% argv$force_modes, FALSE)
    )
  }
  quit(save = "no", status = 0)
}

if (identical(analysis_part, "invitro_tables")) {
  dataset <- "invitro"
  root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
  seed_tables_dir <- argv$seed_tables_dir %||% argv$seed_parameter_tables_dir %||% paper_seed_parameter_tables_dir(root_dir = root_dir)
  paper_generate_umap_tables(
    dataset = dataset,
    input_dir = argv$input_dir %||% default_dataset_input_dir(dataset),
    root_dir = root_dir,
    tables_dir = seed_tables_dir,
    max_seeds = as_int(argv$max_seeds, NA_integer_),
    write_best = as_bool(argv$write_best, TRUE),
    write_initial = as_bool(argv$write_initial, TRUE)
  )
  quit(save = "no", status = 0)
}

if (identical(analysis_part, "invitro_reductions")) {
  dataset <- "invitro"
  root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
  tables_dir <- argv$tables_dir %||% paper_tables_dir(dataset, root_dir = root_dir)
  figures_dir <- argv$figures_dir %||% paper_figures_dir(dataset, root_dir = root_dir)
  support_tables_dir <- argv$support_tables_dir %||% tables_dir
  figures_wclusters_dir <- argv$figures_wclusters_dir %||% paper_figures_wclusters_dir(dataset, root_dir = root_dir)
  tables_wclusters_dir <- argv$tables_wclusters_dir %||% paper_tables_wclusters_dir(dataset, root_dir = root_dir)
  paper_generate_umap_figures(
    dataset = dataset,
    root_dir = root_dir,
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    pca_tables_dir = argv$pca_tables_dir %||% paper_reduction_tables_dir(dataset, "pca", root_dir = root_dir),
    pca_figures_dir = argv$pca_figures_dir %||% paper_reduction_figures_dir(dataset, "pca", root_dir = root_dir),
    pca_figures_wclusters_dir = argv$pca_figures_wclusters_dir %||% paper_reduction_figures_wclusters_dir(dataset, "pca", root_dir = root_dir),
    pca_tables_wclusters_dir = argv$pca_tables_wclusters_dir %||% paper_reduction_tables_wclusters_dir(dataset, "pca", root_dir = root_dir),
    tsne_tables_dir = argv$tsne_tables_dir %||% paper_reduction_tables_dir(dataset, "tsne", root_dir = root_dir),
    tsne_figures_dir = argv$tsne_figures_dir %||% paper_reduction_figures_dir(dataset, "tsne", root_dir = root_dir),
    tsne_figures_wclusters_dir = argv$tsne_figures_wclusters_dir %||% paper_reduction_figures_wclusters_dir(dataset, "tsne", root_dir = root_dir),
    tsne_tables_wclusters_dir = argv$tsne_tables_wclusters_dir %||% paper_reduction_tables_wclusters_dir(dataset, "tsne", root_dir = root_dir),
    initial_csv = argv$initial_csv %||% paper_initial_population_csv(dataset, root_dir = root_dir),
    best_csv = argv$best_csv %||% paper_best_params_csv(dataset, root_dir = root_dir),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    output_prefix = argv$output_prefix %||% "invitro_deoptim_initial_vs_best_umap",
    best_output_prefix = argv$best_output_prefix %||% "invitro_best_params_umap",
    run_combined = as_bool(argv$run_combined, TRUE),
    run_sampled = as_bool(argv$run_sampled, TRUE),
    run_best_only = as_bool(argv$run_best_only, TRUE),
    run_clustered_umaps = as_bool(argv$run_clustered_umaps, TRUE),
    reductions = as_char_vec(argv$reductions, c("umap")),
    preprocess_modes = as_char_vec(argv$preprocess_modes, c("zscore")),
    run_pca_umap = as_bool(argv$run_pca_umap, TRUE),
    run_full_tsne = as_bool(argv$run_full_tsne, FALSE),
    drop_parameter_table_initial = as_bool(argv$drop_parameter_table_initial, FALSE),
    shape_by_pred = as_bool(argv$shape_by_pred, FALSE),
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    tsne_seed = as_int(argv$tsne_seed, 123L),
    tsne_perplexity = as_num(argv$tsne_perplexity, 30),
    tsne_theta = as_num(argv$tsne_theta, 0.5),
    tsne_max_iter = as_int(argv$tsne_max_iter, 1000L),
    sample_initial_n = as_int(argv$sample_initial_n, 500L),
    sample_initial_seed = as_int(argv$sample_initial_seed, 123L),
    sample_initial_by_seed = as_bool(argv$sample_initial_by_seed, TRUE),
    cluster_seed = as_int(argv$cluster_seed, 123L),
    cluster_k_min = as_int(argv$cluster_k_min, 2L),
    cluster_k_max = as_int(argv$cluster_k_max, 8L),
    cluster_silhouette_sample_n = as_int(argv$cluster_silhouette_sample_n, 5000L)
  )
  quit(save = "no", status = 0)
}

if (identical(analysis_part, "pooled_reductions")) {
  root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
  paper_generate_pooled_invivo_invitro_umap_figures(
    root_dir = root_dir,
    tables_dir = argv$tables_dir %||% paper_pooled_tables_dir(root_dir = root_dir),
    figures_dir = argv$figures_dir %||% paper_pooled_figures_dir(root_dir = root_dir),
    figures_wclusters_dir = argv$figures_wclusters_dir %||% paper_pooled_figures_wclusters_dir(root_dir = root_dir),
    tables_wclusters_dir = argv$tables_wclusters_dir %||% paper_pooled_tables_wclusters_dir(root_dir = root_dir),
    pca_tables_dir = argv$pca_tables_dir %||% paper_pooled_reduction_tables_dir("pca", root_dir = root_dir),
    pca_figures_dir = argv$pca_figures_dir %||% paper_pooled_reduction_figures_dir("pca", root_dir = root_dir),
    pca_figures_wclusters_dir = argv$pca_figures_wclusters_dir %||% paper_pooled_reduction_figures_wclusters_dir("pca", root_dir = root_dir),
    pca_tables_wclusters_dir = argv$pca_tables_wclusters_dir %||% paper_pooled_reduction_tables_wclusters_dir("pca", root_dir = root_dir),
    tsne_tables_dir = argv$tsne_tables_dir %||% paper_pooled_reduction_tables_dir("tsne", root_dir = root_dir),
    tsne_figures_dir = argv$tsne_figures_dir %||% paper_pooled_reduction_figures_dir("tsne", root_dir = root_dir),
    tsne_figures_wclusters_dir = argv$tsne_figures_wclusters_dir %||% paper_pooled_reduction_figures_wclusters_dir("tsne", root_dir = root_dir),
    tsne_tables_wclusters_dir = argv$tsne_tables_wclusters_dir %||% paper_pooled_reduction_tables_wclusters_dir("tsne", root_dir = root_dir),
    invivo_best_csv = argv$invivo_best_csv %||% paper_best_params_csv("invivo", root_dir = root_dir),
    invivo_initial_csv = argv$invivo_initial_csv %||% paper_initial_population_csv("invivo", root_dir = root_dir),
    invitro_best_csv = argv$invitro_best_csv %||% paper_best_params_csv("invitro", root_dir = root_dir),
    invitro_initial_csv = argv$invitro_initial_csv %||% paper_initial_population_csv("invitro", root_dir = root_dir),
    invivo_objective_seed_dir = argv$invivo_objective_seed_dir %||% default_dataset_input_dir("invivo"),
    invitro_objective_seed_dir = argv$invitro_objective_seed_dir %||% default_dataset_input_dir("invitro"),
    output_prefix = argv$output_prefix %||% "pooled_invivo_invitro_initial_vs_best_umap",
    run_full = as_bool(argv$run_full, TRUE),
    run_sampled = as_bool(argv$run_sampled, TRUE),
    run_clustered_umaps = as_bool(argv$run_clustered_umaps, TRUE),
    reductions = as_char_vec(argv$reductions, c("umap")),
    preprocess_modes = as_char_vec(argv$preprocess_modes, c("zscore")),
    run_pca_umap = as_bool(argv$run_pca_umap, TRUE),
    drop_parameter_table_initial = as_bool(argv$drop_parameter_table_initial, TRUE),
    drop_invitro_parameter_table_initial = as_bool(argv$drop_invitro_parameter_table_initial, TRUE),
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    tsne_seed = as_int(argv$tsne_seed, 123L),
    tsne_perplexity = as_num(argv$tsne_perplexity, 30),
    tsne_theta = as_num(argv$tsne_theta, 0.5),
    tsne_max_iter = as_int(argv$tsne_max_iter, 1000L),
    sample_initial_n = as_int(argv$sample_initial_n, 500L),
    sample_initial_seed = as_int(argv$sample_initial_seed, 123L),
    sample_initial_by_seed = as_bool(argv$sample_initial_by_seed, TRUE),
    cluster_seed = as_int(argv$cluster_seed, 123L),
    cluster_k_min = as_int(argv$cluster_k_min, 2L),
    cluster_k_max = as_int(argv$cluster_k_max, 8L),
    cluster_silhouette_sample_n = as_int(argv$cluster_silhouette_sample_n, 5000L),
    initial_size = as_num(argv$initial_size, 0.22),
    sampled_initial_size = as_num(argv$sampled_initial_size, NA_real_),
    best_size = as_num(argv$best_size, 1.25)
  )
  quit(save = "no", status = 0)
}

dataset <- "invivo"
root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
tables_dir <- argv$tables_dir %||% paper_tables_dir(dataset, root_dir = root_dir)
figures_dir <- argv$figures_dir %||% paper_figures_dir(dataset, root_dir = root_dir)
support_tables_dir <- argv$support_tables_dir %||% tables_dir
figures_wclusters_dir <- argv$figures_wclusters_dir %||% paper_figures_wclusters_dir(dataset, root_dir = root_dir)
tables_wclusters_dir <- argv$tables_wclusters_dir %||% paper_tables_wclusters_dir(dataset, root_dir = root_dir)
reductions <- as_char_vec(argv$reductions, c("umap"))
preprocess_modes <- as_char_vec(argv$preprocess_modes, c("zscore"))
run_combined <- as_bool(argv$run_combined, TRUE)
run_sampled <- as_bool(argv$run_sampled, TRUE)
run_best_only <- as_bool(argv$run_best_only, TRUE)
run_clustered_umaps <- as_bool(argv$run_clustered_umaps, TRUE)
run_growth_turnover_best_only <- as_bool(argv$run_growth_turnover_best_only, TRUE)
run_growth_turnover_o2_best_only <- as_bool(argv$run_growth_turnover_o2_best_only, TRUE)
run_growth_turnover_o2_ploidy_ratio_best_only <- as_bool(argv$run_growth_turnover_o2_ploidy_ratio_best_only, TRUE)
run_growth_turnover_misseg_best_only <- as_bool(argv$run_growth_turnover_misseg_best_only, TRUE)
run_growth_turnover_misseg_o2_best_only <- as_bool(argv$run_growth_turnover_misseg_o2_best_only, TRUE)
run_growth_turnover_misseg_o2_ploidy_ratio_best_only <- as_bool(argv$run_growth_turnover_misseg_o2_ploidy_ratio_best_only, TRUE)
run_growth_turnover_misseg_o2_attractor_dominant_ploidy_best_only <- as_bool(argv$run_growth_turnover_misseg_o2_attractor_dominant_ploidy_best_only, TRUE)
run_growth_turnover_boxplot <- as_bool(argv$run_growth_turnover_boxplot, TRUE)
run_growth_turnover_misseg_boxplot <- as_bool(argv$run_growth_turnover_misseg_boxplot, TRUE)
shape_by_mode <- as_bool(argv$shape_by_mode %||% argv$shape_by_pred, TRUE)
mode_tables_dir <- argv$mode_tables_dir %||% paper_fixo2_mode_tables_dir(root_dir = root_dir)
mode_reference_o2 <- as_num(argv$mode_reference_o2, 2)
attractor_feature_o2_values <- as_num_vec(argv$attractor_feature_o2_values, paper_default_mode_summary_o2())
cluster_seed <- as_int(argv$cluster_seed, 123L)
cluster_k_min <- as_int(argv$cluster_k_min, 2L)
cluster_k_max <- as_int(argv$cluster_k_max, 8L)
cluster_silhouette_sample_n <- as_int(argv$cluster_silhouette_sample_n, 5000L)

if (!run_combined && !run_sampled && !run_best_only && !run_growth_turnover_best_only &&
    !run_growth_turnover_o2_best_only &&
    !run_growth_turnover_o2_ploidy_ratio_best_only &&
    !run_growth_turnover_misseg_best_only &&
    !run_growth_turnover_misseg_o2_best_only &&
    !run_growth_turnover_misseg_o2_ploidy_ratio_best_only &&
    !run_growth_turnover_misseg_o2_attractor_dominant_ploidy_best_only &&
    !run_growth_turnover_boxplot && !run_growth_turnover_misseg_boxplot) {
  stop("Nothing to plot: enable one of run_combined, run_sampled, run_best_only, run_growth_turnover_best_only, run_growth_turnover_o2_best_only, run_growth_turnover_o2_ploidy_ratio_best_only, run_growth_turnover_misseg_best_only, run_growth_turnover_misseg_o2_best_only, run_growth_turnover_misseg_o2_ploidy_ratio_best_only, run_growth_turnover_misseg_o2_attractor_dominant_ploidy_best_only, run_growth_turnover_boxplot, or run_growth_turnover_misseg_boxplot.")
}

if (run_combined || run_sampled || run_best_only) {
  paper_generate_umap_figures(
    dataset = dataset,
    root_dir = root_dir,
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    pca_tables_dir = argv$pca_tables_dir %||% paper_reduction_tables_dir(dataset, "pca", root_dir = root_dir),
    pca_figures_dir = argv$pca_figures_dir %||% paper_reduction_figures_dir(dataset, "pca", root_dir = root_dir),
    pca_figures_wclusters_dir = argv$pca_figures_wclusters_dir %||% paper_reduction_figures_wclusters_dir(dataset, "pca", root_dir = root_dir),
    pca_tables_wclusters_dir = argv$pca_tables_wclusters_dir %||% paper_reduction_tables_wclusters_dir(dataset, "pca", root_dir = root_dir),
    tsne_tables_dir = argv$tsne_tables_dir %||% paper_reduction_tables_dir(dataset, "tsne", root_dir = root_dir),
    tsne_figures_dir = argv$tsne_figures_dir %||% paper_reduction_figures_dir(dataset, "tsne", root_dir = root_dir),
    tsne_figures_wclusters_dir = argv$tsne_figures_wclusters_dir %||% paper_reduction_figures_wclusters_dir(dataset, "tsne", root_dir = root_dir),
    tsne_tables_wclusters_dir = argv$tsne_tables_wclusters_dir %||% paper_reduction_tables_wclusters_dir(dataset, "tsne", root_dir = root_dir),
    initial_csv = argv$initial_csv %||% paper_initial_population_csv(dataset, root_dir = root_dir),
    best_csv = argv$best_csv %||% paper_best_params_csv(dataset, root_dir = root_dir),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2,
    output_prefix = argv$output_prefix %||% "invivo_deoptim_initial_vs_best_umap",
    best_output_prefix = argv$best_output_prefix %||% "invivo_best_params_umap",
    run_combined = run_combined,
    run_sampled = run_sampled,
    run_best_only = run_best_only,
    run_clustered_umaps = run_clustered_umaps,
    reductions = reductions,
    preprocess_modes = preprocess_modes,
    run_pca_umap = as_bool(argv$run_pca_umap, TRUE),
    run_full_tsne = as_bool(argv$run_full_tsne, FALSE),
    drop_parameter_table_initial = as_bool(argv$drop_parameter_table_initial, TRUE),
    shape_by_pred = shape_by_mode,
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    tsne_seed = as_int(argv$tsne_seed, 123L),
    tsne_perplexity = as_num(argv$tsne_perplexity, 30),
    tsne_theta = as_num(argv$tsne_theta, 0.5),
    tsne_max_iter = as_int(argv$tsne_max_iter, 1000L),
    sample_initial_n = as_int(argv$sample_initial_n, 500L),
    sample_initial_seed = as_int(argv$sample_initial_seed, 123L),
    sample_initial_by_seed = as_bool(argv$sample_initial_by_seed, TRUE),
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n
  )
}

if (run_growth_turnover_best_only) {
  paper_generate_invivo_growth_turnover_umap_figures(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2,
    output_prefix = argv$growth_turnover_output_prefix %||% "invivo_best_params_growth_turnover_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = shape_by_mode,
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    best_size = as_num(argv$best_size, 1.6),
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n,
    rate_version = "base"
  )
}

if (run_growth_turnover_o2_best_only) {
  paper_generate_invivo_growth_turnover_umap_figures(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2,
    output_prefix = argv$growth_turnover_o2_output_prefix %||% "invivo_best_params_growth_turnover_O2_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = shape_by_mode,
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    best_size = as_num(argv$best_size, 1.6),
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n,
    rate_version = "base_with_o2"
  )
}

if (run_growth_turnover_o2_ploidy_ratio_best_only) {
  paper_generate_invivo_growth_turnover_umap_figures(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2,
    output_prefix = argv$growth_turnover_o2_ploidy_ratio_output_prefix %||% "invivo_best_params_growth_turnover_O2_ploidy_ratio_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = shape_by_mode,
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    best_size = as_num(argv$best_size, 1.6),
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n,
    rate_version = "base_with_o2_ploidy_ratio"
  )
}

if (run_growth_turnover_misseg_best_only) {
  paper_generate_invivo_growth_turnover_umap_figures(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2,
    output_prefix = argv$growth_turnover_misseg_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = shape_by_mode,
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    best_size = as_num(argv$best_size, 1.6),
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n,
    rate_version = "with_misseg_death"
  )
}

if (run_growth_turnover_misseg_o2_best_only) {
  paper_generate_invivo_growth_turnover_umap_figures(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2,
    output_prefix = argv$growth_turnover_misseg_o2_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_O2_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = shape_by_mode,
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    best_size = as_num(argv$best_size, 1.6),
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n,
    rate_version = "with_misseg_death_with_o2"
  )
}

if (run_growth_turnover_misseg_o2_ploidy_ratio_best_only) {
  paper_generate_invivo_growth_turnover_umap_figures(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2,
    output_prefix = argv$growth_turnover_misseg_o2_ploidy_ratio_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_O2_ploidy_ratio_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = shape_by_mode,
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    best_size = as_num(argv$best_size, 1.6),
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n,
    rate_version = "with_misseg_death_with_o2_ploidy_ratio"
  )
}

if (run_growth_turnover_misseg_o2_attractor_dominant_ploidy_best_only) {
  paper_generate_invivo_growth_turnover_umap_figures(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    mode_tables_dir = mode_tables_dir,
    mode_reference_o2 = mode_reference_o2,
    attractor_feature_o2_values = attractor_feature_o2_values,
    output_prefix = argv$growth_turnover_misseg_o2_attractor_dominant_ploidy_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_O2_attractor_dominant_ploidy_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = shape_by_mode,
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
    best_size = as_num(argv$best_size, 1.6),
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n,
    rate_version = "with_misseg_death_with_o2_attractor_dominant_ploidy"
  )
}

if (run_growth_turnover_boxplot) {
  paper_generate_invivo_growth_turnover_boxplot(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    best_csv = argv$best_csv %||% paper_best_params_csv(dataset, root_dir = root_dir),
    output_prefix = argv$growth_turnover_boxplot_output_prefix %||% "invivo_best_params_growth_turnover_100d_paired_boxplot",
    jitter_seed = as_int(argv$growth_turnover_boxplot_jitter_seed, 123L),
    jitter_width = as_num(argv$growth_turnover_boxplot_jitter_width, 0.22),
    rate_version = "base"
  )
}

if (run_growth_turnover_misseg_boxplot) {
  paper_generate_invivo_growth_turnover_boxplot(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    best_csv = argv$best_csv %||% paper_best_params_csv(dataset, root_dir = root_dir),
    output_prefix = argv$growth_turnover_misseg_boxplot_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_100d_paired_boxplot",
    jitter_seed = as_int(argv$growth_turnover_boxplot_jitter_seed, 123L),
    jitter_width = as_num(argv$growth_turnover_boxplot_jitter_width, 0.22),
    rate_version = "with_misseg_death"
  )
}
