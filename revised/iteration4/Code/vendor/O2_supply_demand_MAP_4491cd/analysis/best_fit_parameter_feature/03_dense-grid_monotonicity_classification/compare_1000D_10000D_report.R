#!/usr/bin/env Rscript
# Deprecated compatibility entrypoint; report assembly now belongs to report/dense_grid_monotonicity.
.legacy_file <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "compare_1000D_10000D_report.R"]
  if (length(own)) own[[length(own)]] else normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE)
})
source(file.path(dirname(.legacy_file), "deprecated_forward.R"), local = environment())
if (sys.nframe() == 0L) deprecated_dense_direct("report", "initial_ploidy")
