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
dataset <- "invivo"
root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
tables_dir <- argv$tables_dir %||% paper_tables_dir(dataset, root_dir = root_dir)
figures_dir <- argv$figures_dir %||% paper_figures_dir(dataset, root_dir = root_dir)
support_tables_dir <- argv$support_tables_dir %||% tables_dir
figures_wclusters_dir <- argv$figures_wclusters_dir %||% paper_figures_wclusters_dir(dataset, root_dir = root_dir)
tables_wclusters_dir <- argv$tables_wclusters_dir %||% paper_tables_wclusters_dir(dataset, root_dir = root_dir)
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
run_growth_turnover_boxplot <- as_bool(argv$run_growth_turnover_boxplot, TRUE)
run_growth_turnover_misseg_boxplot <- as_bool(argv$run_growth_turnover_misseg_boxplot, TRUE)
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
    !run_growth_turnover_boxplot && !run_growth_turnover_misseg_boxplot) {
  stop("Nothing to plot: enable one of run_combined, run_sampled, run_best_only, run_growth_turnover_best_only, run_growth_turnover_o2_best_only, run_growth_turnover_o2_ploidy_ratio_best_only, run_growth_turnover_misseg_best_only, run_growth_turnover_misseg_o2_best_only, run_growth_turnover_misseg_o2_ploidy_ratio_best_only, run_growth_turnover_boxplot, or run_growth_turnover_misseg_boxplot.")
}

if (run_combined || run_sampled || run_best_only) {
  paper_generate_umap_figures(
    dataset = dataset,
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    support_tables_dir = support_tables_dir,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    initial_csv = argv$initial_csv %||% file.path(tables_dir, "invivo_deoptim_initial_population.csv"),
    best_csv = argv$best_csv %||% file.path(tables_dir, "invivo_best_params_by_seed.csv"),
    objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
    output_prefix = argv$output_prefix %||% "invivo_deoptim_initial_vs_best_umap",
    best_output_prefix = argv$best_output_prefix %||% "invivo_best_params_umap",
    run_combined = run_combined,
    run_sampled = run_sampled,
    run_best_only = run_best_only,
    run_clustered_umaps = run_clustered_umaps,
    drop_parameter_table_initial = as_bool(argv$drop_parameter_table_initial, TRUE),
    shape_by_pred = as_bool(argv$shape_by_pred, TRUE),
    umap_seed = as_int(argv$umap_seed, 123L),
    n_neighbors = as_int(argv$n_neighbors, 80L),
    min_dist = as_num(argv$min_dist, 0.1),
    n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
    pca_n = as_int(argv$pca_n, 10L),
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
    output_prefix = argv$growth_turnover_output_prefix %||% "invivo_best_params_growth_turnover_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = as_bool(argv$shape_by_pred, TRUE),
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
    output_prefix = argv$growth_turnover_o2_output_prefix %||% "invivo_best_params_growth_turnover_O2_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = as_bool(argv$shape_by_pred, TRUE),
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
    output_prefix = argv$growth_turnover_o2_ploidy_ratio_output_prefix %||% "invivo_best_params_growth_turnover_O2_ploidy_ratio_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = as_bool(argv$shape_by_pred, TRUE),
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
    output_prefix = argv$growth_turnover_misseg_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = as_bool(argv$shape_by_pred, TRUE),
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
    output_prefix = argv$growth_turnover_misseg_o2_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_O2_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = as_bool(argv$shape_by_pred, TRUE),
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
    output_prefix = argv$growth_turnover_misseg_o2_ploidy_ratio_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_O2_ploidy_ratio_100d_umap",
    run_clustered_umaps = run_clustered_umaps,
    shape_by_pred = as_bool(argv$shape_by_pred, TRUE),
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

if (run_growth_turnover_boxplot) {
  paper_generate_invivo_growth_turnover_boxplot(
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    growth_turnover_csv = argv$growth_turnover_csv %||% file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
    best_csv = argv$best_csv %||% file.path(tables_dir, "invivo_best_params_by_seed.csv"),
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
    best_csv = argv$best_csv %||% file.path(tables_dir, "invivo_best_params_by_seed.csv"),
    output_prefix = argv$growth_turnover_misseg_boxplot_output_prefix %||% "invivo_best_params_growth_turnover_with_misseg_death_100d_paired_boxplot",
    jitter_seed = as_int(argv$growth_turnover_boxplot_jitter_seed, 123L),
    jitter_width = as_num(argv$growth_turnover_boxplot_jitter_width, 0.22),
    rate_version = "with_misseg_death"
  )
}
