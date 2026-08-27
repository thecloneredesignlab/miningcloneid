#!/usr/bin/env Rscript

# Deprecated compatibility entrypoint.  The runner now orchestrates the
# simulation, analysis, and visualization layers without mixing them here.
.o2pl_compat_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else x$ofile, character(1)))
  own <- frames[basename(frames) == "clustering_analysis.R"]
  if (length(own)) dirname(normalizePath(own[[length(own)]], mustWork = FALSE)) else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = FALSE))
  }
})
.o2pl_root <- normalizePath(file.path(.o2pl_compat_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_root, "runner", "parameter_landscape", "run_parameter_landscape.R"), local = environment(), chdir = TRUE)
.o2pl_argv <- parse_args(commandArgs(trailingOnly = TRUE))
.o2pl_argv$run_parts <- .o2pl_argv$analysis_part %||% .o2pl_argv$part %||% "invivo_reductions"
o2pl_parameter_landscape_runner_main(.o2pl_argv)
