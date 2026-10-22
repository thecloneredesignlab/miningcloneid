#!/usr/bin/env Rscript
# Canonical retained-data entry; no version selection, simulation, or archive fallback.
script <- normalizePath(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1]))
source(file.path(dirname(script),"util/analysis/figure78_render_only.R"))
f78_render_cli(script,default_scope="all",default_products=c("figure7_ver2","supp_fig7_1_ver2","supp_fig7_2_ver2","supp_fig7_3_ver2","supp_fig7_4_ver2","supp_fig7_5_ver2","supp_fig7_6_ver2","supp_fig7_7_ver2"))
