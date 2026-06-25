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
    tables_dir = tables_dir,
    max_seeds = as_int(argv$max_seeds, NA_integer_),
    write_best = write_best,
    write_initial = write_initial
  )
}

if (write_growth_turnover) {
  paper_generate_invivo_growth_turnover_table(
    input_dir = argv$input_dir %||% default_dataset_input_dir(dataset),
    best_csv = argv$best_csv %||% file.path(tables_dir, "invivo_best_params_by_seed.csv"),
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
    best_csv = argv$best_csv %||% file.path(tables_dir, "invivo_best_params_by_seed.csv"),
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
