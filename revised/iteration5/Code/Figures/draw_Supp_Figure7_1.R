#!/usr/bin/env Rscript
# Canonical retained-data entry; no version selection, simulation, or archive fallback.
script <- normalizePath(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1]))
source(file.path(dirname(script),"util/analysis/figure78_render_only.R"))
f78_render_cli(script,default_scope="all",default_products=c("supp_fig7_1"))
