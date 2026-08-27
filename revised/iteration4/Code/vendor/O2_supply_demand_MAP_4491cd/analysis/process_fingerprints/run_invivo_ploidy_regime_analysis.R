#!/usr/bin/env Rscript

# Compatibility orchestrator for the former monolithic ploidy-regime entry.
# Each canonical stage consumes only the preceding stage's materialized files.

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
      if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
    }, character(1)))
    own <- frame_files[
      basename(frame_files) == "run_invivo_ploidy_regime_analysis.R"
    ]
    if (length(own)) {
      dirname(own[[length(own)]])
    } else if (length(frame_files)) {
      dirname(frame_files[[length(frame_files)]])
    } else {
      getwd()
    }
  }
})
locate_process_workflow_root <- function(starts) {
  for (start in unique(starts)) {
    current <- normalizePath(start, mustWork = FALSE)
    for (depth in 0:10) {
      candidates <- c(
        current,
        file.path(current, "oxygen", "code", "O2_supply_demand_MAP"),
        file.path(current, "code", "O2_supply_demand_MAP")
      )
      hits <- candidates[file.exists(file.path(
        candidates,
        "analysis",
        "process_fingerprints",
        "process_fingerprint_utils.R"
      ))]
      if (length(hits)) return(normalizePath(hits[[1L]], mustWork = TRUE))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate the O2_supply_demand_MAP workflow root.")
}
WORKFLOW_ROOT <- locate_process_workflow_root(c(SCRIPT_DIR, getwd()))
SCRIPT_DIR <- file.path(WORKFLOW_ROOT, "analysis", "process_fingerprints")
rm(locate_process_workflow_root)

source(file.path(WORKFLOW_ROOT, "simulation", "process_fingerprints", "generate_process_fingerprint_outputs.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "analyze_invivo_ploidy_regimes.R"), local = TRUE)
source(file.path(WORKFLOW_ROOT, "vis", "process_fingerprints", "plot_invivo_ploidy_regimes.R"), local = TRUE)
source(file.path(WORKFLOW_ROOT, "report", "process_fingerprints", "render_invivo_ploidy_regime_report.R"), local = TRUE)

main <- function(argv = o2ipa_parse_args()) {
  run_dir <- o2ipa_as_chr(argv$run_dir)
  if (!nzchar(run_dir) || !dir.exists(run_dir)) stop("Missing or invalid --run_dir.")
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  out_dir <- o2ipa_as_chr(argv$out_dir, file.path(run_dir, "analysis", "ploidy_regimes"))
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("Compatibility entry: orchestrating simulation -> ploidy analysis -> visualization -> report.")

  sim_args <- argv
  sim_args$run_dir <- run_dir
  sim_args$simulation_dir <- out_dir
  run_process_fingerprint_simulation(sim_args)

  analysis_args <- argv
  analysis_args$simulation_dir <- out_dir
  analysis_args$analysis_dir <- out_dir
  run_invivo_ploidy_regime_analysis_stage(analysis_args)

  viz_args <- list(simulation_dir = out_dir, analysis_dir = out_dir, viz_dir = out_dir)
  run_invivo_ploidy_regime_visualization(viz_args)

  report_dir <- file.path(out_dir, "report")
  report_args <- list(simulation_dir = out_dir, analysis_dir = out_dir, viz_dir = out_dir, report_dir = report_dir)
  run_invivo_ploidy_regime_report(report_args)
  message("Completed compatibility ploidy-regime workflow: ", out_dir)
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
