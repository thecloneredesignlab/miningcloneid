#!/usr/bin/env Rscript
# Compatibility filename: current explicit interfaces only.
script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
source(file.path(dirname(script),"util/analysis/figure78_compat.R"))
f78_compat_cli(script,figure=7L,role="validation")
