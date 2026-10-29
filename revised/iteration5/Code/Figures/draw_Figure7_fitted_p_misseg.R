#!/usr/bin/env Rscript
# Render only the authenticated fitted-p products; never simulate.
.script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
source(file.path(dirname(.script),"util/analysis/figure7_fitted_p_misseg.R"),local=globalenv())
source(file.path(.f7fit_dir,"figure7_fitted_plotting.R"),local=globalenv())
f7fit_plot(f7fit_init(.script,f7fit_args()))
