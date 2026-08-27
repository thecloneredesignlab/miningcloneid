#!/usr/bin/env Rscript
# Deprecated compatibility entrypoint; joins live in analysis and graphics live in vis.
.legacy_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)
source(file.path(dirname(.legacy_file), "deprecated_forward.R"), local = environment())
if (sys.nframe() == 0L) deprecated_combined_direct("run", "--run_report=false")
