#!/usr/bin/env Rscript

# COMPATIBILITY ORCHESTRATOR.  Canonical work is split across simulation,
# analysis, visualization, and report entrypoints.

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "simulation", "process_fingerprints", "generate_o2_ploidy_event_inputs.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "analyze_invivo_medium_o2_window_tables.R"), local = TRUE)
source(file.path(WORKFLOW_ROOT, "vis", "process_fingerprints", "plot_medium_o2_windows.R"), local = TRUE)
source(file.path(WORKFLOW_ROOT, "report", "process_fingerprints", "render_medium_o2_window_report.R"), local = TRUE)

main <- function(argv = o2ipa_parse_args()) {
  run_dir <- o2ipa_as_chr(argv$run_dir)
  if (!nzchar(run_dir) || !dir.exists(run_dir)) stop("Missing or invalid --run_dir.")
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  out_dir <- o2ipa_as_chr(argv$out_dir, file.path(dirname(run_dir), "analysis", "invivo_medium_o2_windows_500seed"))
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  message("Compatibility entry: orchestrating event simulation -> medium-O2 analysis -> visualization -> report.")
  sim_args <- argv
  sim_args$run_dir <- run_dir
  sim_args$simulation_dir <- out_dir
  run_o2_ploidy_event_input_simulation(sim_args)
  analysis_args <- argv
  analysis_args$simulation_dir <- out_dir
  analysis_args$analysis_dir <- out_dir
  run_medium_o2_window_analysis(analysis_args)
  viz_args <- list(simulation_dir = out_dir, analysis_dir = out_dir, viz_dir = out_dir)
  run_medium_o2_window_visualization(viz_args)
  report_args <- list(simulation_dir = out_dir, analysis_dir = out_dir, viz_dir = out_dir, report_dir = file.path(out_dir, "report"))
  run_medium_o2_window_report(report_args)
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
