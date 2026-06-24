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
dataset <- "invivo"
root_dir <- normalizePath(path.expand(argv$result_root %||% default_paper_figures_tables_dir()), mustWork = FALSE)
tables_dir <- argv$tables_dir %||% paper_tables_dir(dataset, root_dir = root_dir)
write_best <- as_bool(argv$write_best, TRUE)
write_initial <- as_bool(argv$write_initial, TRUE)
write_growth_turnover <- as_bool(argv$write_growth_turnover, TRUE)
write_ploidy_ratios <- as_bool(argv$write_ploidy_ratios, FALSE)

if (!write_best && !write_initial && !write_growth_turnover && !write_ploidy_ratios) {
  stop("Nothing to write: enable write_best, write_initial, write_growth_turnover, or write_ploidy_ratios.")
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
