#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

source(file.path(local_script_dir(), "util.R"))

argv <- parse_args(commandArgs(trailingOnly = TRUE))
dataset <- "invitro"
root_dir <- normalizePath(path.expand(argv$result_root %||% default_paper_figures_tables_dir()), mustWork = FALSE)
tables_dir <- argv$tables_dir %||% paper_tables_dir(dataset, root_dir = root_dir)
figures_dir <- argv$figures_dir %||% paper_figures_dir(dataset, root_dir = root_dir)
support_tables_dir <- argv$support_tables_dir %||% tables_dir

paper_generate_umap_figures(
  dataset = dataset,
  tables_dir = tables_dir,
  figures_dir = figures_dir,
  support_tables_dir = support_tables_dir,
  initial_csv = argv$initial_csv %||% file.path(tables_dir, "invitro_deoptim_initial_population.csv"),
  best_csv = argv$best_csv %||% file.path(tables_dir, "invitro_best_params_by_seed.csv"),
  objective_seed_dir = argv$objective_seed_dir %||% default_dataset_input_dir(dataset),
  output_prefix = argv$output_prefix %||% "invitro_deoptim_initial_vs_best_umap",
  best_output_prefix = argv$best_output_prefix %||% "invitro_best_params_umap",
  run_combined = as_bool(argv$run_combined, TRUE),
  run_sampled = as_bool(argv$run_sampled, TRUE),
  run_best_only = as_bool(argv$run_best_only, TRUE),
  drop_parameter_table_initial = as_bool(argv$drop_parameter_table_initial, FALSE),
  shape_by_pred = as_bool(argv$shape_by_pred, FALSE),
  umap_seed = as_int(argv$umap_seed, 123L),
  n_neighbors = as_int(argv$n_neighbors, 80L),
  min_dist = as_num(argv$min_dist, 0.1),
  n_threads = as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L))),
  pca_n = as_int(argv$pca_n, 10L),
  sample_initial_n = as_int(argv$sample_initial_n, 500L),
  sample_initial_seed = as_int(argv$sample_initial_seed, 123L),
  sample_initial_by_seed = as_bool(argv$sample_initial_by_seed, TRUE)
)
