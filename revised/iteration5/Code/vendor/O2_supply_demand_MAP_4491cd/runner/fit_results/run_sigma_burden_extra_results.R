#!/usr/bin/env Rscript

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)

run_sigma_burden_fit_results <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  results_root <- normalizePath(file.path(WORKFLOW_ROOT, "..", "..", "results"), mustWork = FALSE)
  default_template <- file.path(
    results_root,
    "fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma}"
  )
  default_out_dir <- file.path(results_root, "comp")
  run_template <- o2sd_null_coalesce(argv$run_dir_template, default_template)
  sigma_caps <- o2sd_null_coalesce(argv$sigma_caps, "0.05,0.15,0.3,0.6")
  out_dir <- normalizePath(o2sd_null_coalesce(argv$out_dir, default_out_dir), mustWork = FALSE)
  analysis_args <- c(
    paste0("--run_dir_template=", run_template),
    paste0("--sigma_caps=", sigma_caps),
    paste0("--out_dir=", out_dir)
  )
  o2fr_run_rscript_stage(file.path(WORKFLOW_ROOT, "analysis", "fit_results", "analyze_sigma_burden_extra_results.R"), analysis_args, "sigma-burden analysis")
  o2fr_run_rscript_stage(file.path(WORKFLOW_ROOT, "vis", "fit_results", "plot_sigma_burden_extra_results.R"), c(paste0("--analysis_dir=", out_dir), paste0("--viz_dir=", out_dir)), "sigma-burden visualization")
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_sigma_burden_fit_results()
