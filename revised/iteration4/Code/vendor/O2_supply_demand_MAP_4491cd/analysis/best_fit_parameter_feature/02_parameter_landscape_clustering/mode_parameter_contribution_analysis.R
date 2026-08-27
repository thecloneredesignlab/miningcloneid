#!/usr/bin/env Rscript
# Deprecated compatibility wrapper; canonical implementation consumes materialized simulation tables.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])), "deprecated_forward.R"))
bpf_parameter_landscape_forward("mode_parameter_contribution_analysis.R")
