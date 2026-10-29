#!/usr/bin/env Rscript
# Fitted-parameter condition adapter to the unchanged Figure7 numerical library.
.script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
source(file.path(dirname(.script),"util/analysis/figure7_fitted_p_misseg.R"),local=globalenv())
f7fit_cli(.script)
