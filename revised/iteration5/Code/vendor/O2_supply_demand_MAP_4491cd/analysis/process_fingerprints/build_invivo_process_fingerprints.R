#!/usr/bin/env Rscript

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})

message("build_invivo_process_fingerprints.R is provided for compatibility; run_invivo_process_analysis.R executes the full cached fingerprint pipeline.")
source(file.path(script_dir, "run_invivo_process_analysis.R"), local = TRUE)
main()

