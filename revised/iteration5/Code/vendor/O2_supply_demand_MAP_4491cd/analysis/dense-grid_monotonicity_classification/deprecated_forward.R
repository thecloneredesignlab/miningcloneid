#!/usr/bin/env Rscript

# Deprecated compatibility loader for the historical mirror directory.
.deprecated_dense_file <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "deprecated_forward.R"]
  if (length(own)) own[[length(own)]] else ""
})
.deprecated_dense_dir <- dirname(.deprecated_dense_file)
.deprecated_dense_root <- normalizePath(file.path(.deprecated_dense_dir, "..", ".."), mustWork = TRUE)
source(file.path(.deprecated_dense_root, "analysis", "best_fit_parameter_feature", "03_dense-grid_monotonicity_classification", "deprecated_forward.R"), local = environment(), chdir = TRUE)
