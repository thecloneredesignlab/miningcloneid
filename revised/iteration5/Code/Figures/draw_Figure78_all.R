#!/usr/bin/env Rscript
# Render existing validated data only; no archived fallback or model execution.
script <- normalizePath(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1]))
source(file.path(dirname(script),"util/analysis/figure78_render_only.R"))
f78_render_cli(script,default_scope="all")
