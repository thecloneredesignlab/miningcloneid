#!/usr/bin/env Rscript
# Deprecated compatibility entrypoint for report assembly.
.legacy_file <- local({ frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1))); own <- frames[basename(frames) == "dense_grid_monotonicity_report.R"]; if (length(own)) own[[length(own)]] else normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE) })
source(file.path(dirname(.legacy_file), "deprecated_forward.R"), local = environment())
write_dense_grid_monotonicity_report <- function(result_dir, ...) deprecated_dense_run(list(out_dir = result_dir), "report", "monotonicity")
main <- write_dense_grid_monotonicity_report
if (sys.nframe() == 0L) deprecated_dense_direct("report", "monotonicity")
