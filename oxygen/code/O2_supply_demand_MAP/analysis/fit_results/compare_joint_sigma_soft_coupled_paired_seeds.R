#!/usr/bin/env Rscript

# DEPRECATED COMPATIBILITY ORCHESTRATOR.
# Canonical stages live under analysis/fit_results, vis/fit_results,
# report/fit_results, and runner/fit_results.

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "runner", "fit_results", "run_joint_sigma_soft_coupled_paired_seeds.R"), local = TRUE)

if (sys.nframe() == 0L) run_joint_sigma_fit_results()
