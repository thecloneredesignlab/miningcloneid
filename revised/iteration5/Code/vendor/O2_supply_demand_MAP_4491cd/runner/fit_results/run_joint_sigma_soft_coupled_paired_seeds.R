#!/usr/bin/env Rscript

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)

run_joint_sigma_fit_results <- function(argv = o2fr_parse_args()) {
  out_dir <- normalizePath(o2fr_null_coalesce(argv$out_dir, argv$results_dir), mustWork = FALSE)
  common <- c(
    if (!is.null(argv$results_dir)) paste0("--results_dir=", argv$results_dir),
    if (!is.null(argv$run_dirs)) paste0("--run_dirs=", argv$run_dirs),
    paste0("--out_dir=", out_dir)
  )
  o2fr_run_rscript_stage(file.path(WORKFLOW_ROOT, "analysis", "fit_results", "analyze_joint_sigma_soft_coupled_paired_seeds.R"), common, "joint-sigma analysis")
  o2fr_run_rscript_stage(file.path(WORKFLOW_ROOT, "vis", "fit_results", "plot_joint_sigma_soft_coupled_paired_seeds.R"), c(paste0("--analysis_dir=", out_dir), paste0("--viz_dir=", out_dir)), "joint-sigma visualization")
  o2fr_run_rscript_stage(file.path(WORKFLOW_ROOT, "report", "fit_results", "render_joint_sigma_soft_coupled_paired_seeds_report.R"), c(paste0("--analysis_dir=", out_dir), paste0("--viz_dir=", out_dir), paste0("--report_dir=", out_dir)), "joint-sigma report")
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_joint_sigma_fit_results()
