#!/usr/bin/env Rscript
# Deprecated compatibility wrapper for the historical misspelled entrypoint.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])), "deprecated_forward.R"))
bpf_parameter_landscape_forward("full_data_in_vivo_clustring.R")
