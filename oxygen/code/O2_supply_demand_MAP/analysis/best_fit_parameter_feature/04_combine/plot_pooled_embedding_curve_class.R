#!/usr/bin/env Rscript
# Deprecated compatibility mirror; the canonical workflow owns separate analysis and vis stages.
.legacy_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)
.helper <- file.path(dirname(.legacy_file), "..", "04_combine_parameter_landscape", "deprecated_forward.R")
source(.helper, local = environment(), chdir = TRUE)
if (sys.nframe() == 0L) deprecated_combined_direct("run", "--run_report=false")
