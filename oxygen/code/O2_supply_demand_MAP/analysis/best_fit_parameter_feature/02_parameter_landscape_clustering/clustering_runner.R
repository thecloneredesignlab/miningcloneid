#!/usr/bin/env Rscript
# Deprecated compatibility wrapper; canonical implementation is in runner/parameter_landscape.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])), "deprecated_forward.R"))
bpf_parameter_landscape_forward("clustering_runner.R")
