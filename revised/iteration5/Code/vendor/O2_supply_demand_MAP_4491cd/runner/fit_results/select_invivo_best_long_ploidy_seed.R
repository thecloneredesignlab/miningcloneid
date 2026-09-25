#!/usr/bin/env Rscript

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)

run_long_ploidy_seed_selection <- function(argv = o2fr_parse_args()) {
  invivo_dir <- normalizePath(argv$invivo_dir, mustWork = TRUE)
  simulation_dir <- normalizePath(o2fr_null_coalesce(argv$simulation_dir, file.path(invivo_dir, "simulation", "fit_results", "long_ploidy")), mustWork = FALSE)
  out_tsv <- normalizePath(o2fr_null_coalesce(argv$out_tsv, file.path(invivo_dir, "best_long_ploidy_gt2_seed.tsv")), mustWork = FALSE)
  best_dir_file <- normalizePath(o2fr_null_coalesce(argv$best_dir_file, file.path(invivo_dir, "best_long_ploidy_gt2_seed.dir")), mustWork = FALSE)
  sim_args <- c(
    paste0("--invivo_dir=", invivo_dir), paste0("--simulation_dir=", simulation_dir),
    paste0("--horizon=", o2fr_null_coalesce(argv$horizon, "1000")), paste0("--cohort=", o2fr_null_coalesce(argv$cohort, "2N")),
    paste0("--dose=", o2fr_null_coalesce(argv$dose, "ALL"))
  )
  o2fr_run_rscript_stage(file.path(WORKFLOW_ROOT, "simulation", "fit_results", "materialize_invivo_long_ploidy_metrics.R"), sim_args, "long-ploidy simulation")
  analysis_args <- c(
    paste0("--simulation_dir=", simulation_dir), paste0("--out_tsv=", out_tsv),
    paste0("--best_dir_file=", best_dir_file), paste0("--threshold_N=", o2fr_null_coalesce(argv$threshold_N, "44"))
  )
  o2fr_run_rscript_stage(file.path(WORKFLOW_ROOT, "analysis", "fit_results", "select_invivo_best_long_ploidy_seed_from_metrics.R"), analysis_args, "long-ploidy selection analysis")
  invisible(out_tsv)
}

if (sys.nframe() == 0L) run_long_ploidy_seed_selection()
