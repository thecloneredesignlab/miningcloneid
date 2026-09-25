#!/usr/bin/env Rscript

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})

message("analyze_ploidy_regimes.R is a compatibility entry; run_invivo_ploidy_regime_analysis.R performs labeling, process fingerprints, clustering, concordance, and reporting.")
source(file.path(script_dir, "run_invivo_ploidy_regime_analysis.R"), local = TRUE)
main()

