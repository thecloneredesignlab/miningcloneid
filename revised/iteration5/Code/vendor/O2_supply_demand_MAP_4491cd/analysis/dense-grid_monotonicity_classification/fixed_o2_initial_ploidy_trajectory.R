#!/usr/bin/env Rscript
# Deprecated compatibility entrypoint for initial-ploidy simulation and analysis.
.legacy_file <- local({ frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1))); own <- frames[basename(frames) == "fixed_o2_initial_ploidy_trajectory.R"]; if (length(own)) own[[length(own)]] else normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE) })
source(file.path(dirname(.legacy_file), "deprecated_forward.R"), local = environment())
generate_outputs <- function(args = list()) deprecated_dense_run(args, "run", "initial_ploidy")
main <- generate_outputs
if (sys.nframe() == 0L) deprecated_dense_direct("run", "initial_ploidy")
