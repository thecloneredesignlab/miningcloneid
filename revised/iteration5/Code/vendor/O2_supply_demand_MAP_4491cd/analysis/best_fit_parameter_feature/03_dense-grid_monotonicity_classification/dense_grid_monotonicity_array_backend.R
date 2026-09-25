#!/usr/bin/env Rscript
# Deprecated compatibility entrypoint; HPC numerical work now belongs to simulation/o2/dense_grid.
.legacy_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)
source(file.path(dirname(.legacy_file), "deprecated_forward.R"), local = environment())
if (sys.nframe() == 0L) deprecated_dense_direct("run", "monotonicity")
