#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})

source(file.path(script_dir, "data_Figure5.R"))
figure5_build_local_solution_ensemble(WORKSPACE_ROOT)
