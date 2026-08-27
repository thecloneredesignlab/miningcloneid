#!/usr/bin/env Rscript

# DEPRECATED COMPATIBILITY ORCHESTRATOR.
# Canonical execution is runner/fit_results/run_extra_results.R and invokes
# simulation -> analysis -> visualization -> report stages in that order.
.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "runner", "fit_results", "run_extra_results.R"), local = TRUE)

if (sys.nframe() == 0L) run_extra_results_pipeline()
