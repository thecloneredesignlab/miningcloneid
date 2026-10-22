#!/usr/bin/env Rscript
# One current entry point: raw full joint fits -> simulation -> data -> figures.
.script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
.analysis <- file.path(dirname(.script),"util/analysis")
for(.module in c("figure78_pipeline.R","figure78_model.R","figure7_simulation.R"))
  source(file.path(.analysis,.module),local=globalenv())
fp_cli(7L,.script)
