#!/usr/bin/env Rscript
root <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT", getwd()), mustWork=TRUE)
analysis <- file.path(root,"Code","Figures","util","analysis")
source(file.path(analysis,"figure7_robustness.R"))
source(file.path(analysis,"figure7_finite_time_q10.R"))
source(file.path(analysis,"figure7_stochastic_passage.R"))
source(file.path(analysis,"figure7_stochastic_validation.R"))
directory <- file.path(root,"audit","figure7_stochastic_20260904","local_validation")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path(analysis,"figure7_full_range_propagator.cpp"),
  cacheDir=file.path(directory,"rcpp"))
Rcpp::sourceCpp(file.path(analysis,"figure7_stochastic_propagator.cpp"),
  cacheDir=file.path(directory,"rcpp"))
f7s_validate(list(),list(run_root=directory))
