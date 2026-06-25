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
dataset <- "invitro"
root_dir <- normalizePath(path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()), mustWork = FALSE)

paper_generate_umap_tables(
  dataset = dataset,
  input_dir = argv$input_dir %||% default_dataset_input_dir(dataset),
  tables_dir = argv$tables_dir %||% paper_tables_dir(dataset, root_dir = root_dir),
  max_seeds = as_int(argv$max_seeds, NA_integer_),
  write_best = as_bool(argv$write_best, TRUE),
  write_initial = as_bool(argv$write_initial, TRUE)
)
