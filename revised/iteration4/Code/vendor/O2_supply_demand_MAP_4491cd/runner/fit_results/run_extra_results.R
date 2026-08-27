#!/usr/bin/env Rscript

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)

run_extra_results_pipeline <- function(argv = o2fr_parse_args()) {
  run_dir <- normalizePath(argv$run_dir, mustWork = TRUE)
  out_dir <- normalizePath(o2fr_null_coalesce(argv$out_dir, file.path(run_dir, "extra_results")), mustWork = FALSE)
  simulation_dir <- normalizePath(o2fr_null_coalesce(argv$simulation_dir, out_dir), mustWork = FALSE)
  analysis_dir <- normalizePath(o2fr_null_coalesce(argv$analysis_dir, out_dir), mustWork = FALSE)
  viz_dir <- normalizePath(o2fr_null_coalesce(argv$viz_dir, out_dir), mustWork = FALSE)
  report_dir <- normalizePath(o2fr_null_coalesce(argv$report_dir, out_dir), mustWork = FALSE)
  o2fr_run_rscript_stage(
    file.path(WORKFLOW_ROOT, "simulation", "fit_results", "materialize_extra_results_predictions.R"),
    c(paste0("--run_dir=", run_dir), paste0("--simulation_dir=", simulation_dir)), "extra-results simulation"
  )
  o2fr_run_rscript_stage(
    file.path(WORKFLOW_ROOT, "analysis", "fit_results", "analyze_extra_results.R"),
    c(
      paste0("--run_dir=", run_dir), paste0("--simulation_dir=", simulation_dir), paste0("--analysis_dir=", analysis_dir),
      paste0("--near_thresh=", o2fr_null_coalesce(argv$near_thresh, "0.05")),
      paste0("--allow_partial_seed_dirs=", o2fr_null_coalesce(argv$allow_partial_seed_dirs, "FALSE"))
    ), "extra-results analysis"
  )
  o2fr_run_rscript_stage(
    file.path(WORKFLOW_ROOT, "vis", "fit_results", "plot_extra_results.R"),
    c(paste0("--simulation_dir=", simulation_dir), paste0("--analysis_dir=", analysis_dir), paste0("--viz_dir=", viz_dir)), "extra-results visualization"
  )
  o2fr_run_rscript_stage(
    file.path(WORKFLOW_ROOT, "report", "fit_results", "render_extra_results_report.R"),
    c(
      paste0("--simulation_dir=", simulation_dir), paste0("--analysis_dir=", analysis_dir), paste0("--viz_dir=", viz_dir),
      paste0("--report_dir=", report_dir), paste0("--out_path=", file.path(report_dir, "extra_results_report.html"))
    ), "extra-results report"
  )
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_extra_results_pipeline()
