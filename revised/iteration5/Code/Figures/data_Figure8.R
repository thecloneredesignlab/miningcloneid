#!/usr/bin/env Rscript
# Data generation uses the same complete current pipeline, not a validator stub.
.data_script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
source(file.path(dirname(.data_script),"run_Figure8.R"),local=globalenv())
