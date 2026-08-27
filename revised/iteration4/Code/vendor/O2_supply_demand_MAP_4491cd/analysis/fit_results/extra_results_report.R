#!/usr/bin/env Rscript

# DEPRECATED COMPATIBILITY WRAPPER. Report implementation moved to report/fit_results.
.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "report", "fit_results", "render_extra_results_report.R"), local = TRUE)

if (sys.nframe() == 0L) main()
