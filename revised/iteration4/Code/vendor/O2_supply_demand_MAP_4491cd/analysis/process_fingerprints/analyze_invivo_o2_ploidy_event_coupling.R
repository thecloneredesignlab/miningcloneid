#!/usr/bin/env Rscript

# COMPATIBILITY ORCHESTRATOR.  Canonical work is split across simulation,
# analysis, visualization, and report entrypoints.

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "simulation", "process_fingerprints", "generate_o2_ploidy_event_inputs.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "analyze_invivo_o2_ploidy_event_coupling_tables.R"), local = TRUE)
source(file.path(WORKFLOW_ROOT, "vis", "process_fingerprints", "plot_o2_ploidy_event_coupling.R"), local = TRUE)
source(file.path(WORKFLOW_ROOT, "report", "process_fingerprints", "render_o2_ploidy_event_coupling_report.R"), local = TRUE)

main <- function(argv = parse_args()) {
  run_dir <- as_chr(argv$run_dir)
  if (!nzchar(run_dir) || !dir.exists(run_dir)) stop("Missing or invalid --run_dir.")
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  out_dir <- as_chr(argv$out_dir, file.path(dirname(run_dir), "analysis", "invivo_o2_ploidy_event_coupling_500seed"))
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  message("Compatibility entry: orchestrating event simulation -> analysis -> visualization -> report.")
  sim_args <- argv
  sim_args$run_dir <- run_dir
  sim_args$simulation_dir <- out_dir
  run_o2_ploidy_event_input_simulation(sim_args)
  analysis_args <- argv
  analysis_args$simulation_dir <- out_dir
  analysis_args$analysis_dir <- out_dir
  run_o2_ploidy_event_coupling_analysis(analysis_args)
  viz_args <- list(simulation_dir = out_dir, analysis_dir = out_dir, viz_dir = out_dir)
  run_o2_ploidy_event_visualization(viz_args)
  report_args <- list(simulation_dir = out_dir, analysis_dir = out_dir, viz_dir = out_dir, report_dir = file.path(out_dir, "report"))
  run_o2_ploidy_event_report(report_args)
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
