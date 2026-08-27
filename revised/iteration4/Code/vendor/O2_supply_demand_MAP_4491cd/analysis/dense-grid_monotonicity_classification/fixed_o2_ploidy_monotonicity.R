#!/usr/bin/env Rscript
# Deprecated compatibility entrypoint for simulation plus curve analysis.
.legacy_file <- local({ frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1))); own <- frames[basename(frames) == "fixed_o2_ploidy_monotonicity.R"]; if (length(own)) own[[length(own)]] else normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE) })
source(file.path(dirname(.legacy_file), "deprecated_forward.R"), local = environment())
generate_outputs <- function(args = list()) deprecated_dense_run(args, "run", "monotonicity")
main <- generate_outputs
if (sys.nframe() == 0L) deprecated_dense_direct("run", "monotonicity")
