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
root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)
reductions <- as_char_vec(argv$reductions, c("umap"))
preprocess_modes <- as_char_vec(argv$preprocess_modes, c("zscore"))

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
  invivo_best_csv = argv$invivo_best_csv %||% file.path(paper_tables_dir("invivo", root_dir = root_dir), "invivo_best_params_by_seed.csv"),
  invivo_initial_csv = argv$invivo_initial_csv %||% file.path(paper_tables_dir("invivo", root_dir = root_dir), "invivo_deoptim_initial_population.csv"),
  invitro_best_csv = argv$invitro_best_csv %||% file.path(paper_tables_dir("invitro", root_dir = root_dir), "invitro_best_params_by_seed.csv"),
  invitro_initial_csv = argv$invitro_initial_csv %||% file.path(paper_tables_dir("invitro", root_dir = root_dir), "invitro_deoptim_initial_population.csv"),
  invivo_objective_seed_dir = argv$invivo_objective_seed_dir %||% default_dataset_input_dir("invivo"),
  invitro_objective_seed_dir = argv$invitro_objective_seed_dir %||% default_dataset_input_dir("invitro"),
  output_prefix = argv$output_prefix %||% "pooled_invivo_invitro_initial_vs_best_umap",
  run_full = as_bool(argv$run_full, TRUE),
  run_sampled = as_bool(argv$run_sampled, TRUE),
  run_clustered_umaps = as_bool(argv$run_clustered_umaps, TRUE),
  reductions = reductions,
  preprocess_modes = preprocess_modes,
  run_pca_umap = as_bool(argv$run_pca_umap, TRUE),
  run_full_tsne = as_bool(argv$run_full_tsne, FALSE),
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
