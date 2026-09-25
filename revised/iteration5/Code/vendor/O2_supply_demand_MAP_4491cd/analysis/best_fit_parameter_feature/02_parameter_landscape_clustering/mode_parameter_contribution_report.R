#!/usr/bin/env Rscript
# Deprecated compatibility wrapper; report assembly is canonical under report/parameter_landscape.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])), "deprecated_forward.R"))
bpf_parameter_landscape_forward("mode_parameter_contribution_report.R")
