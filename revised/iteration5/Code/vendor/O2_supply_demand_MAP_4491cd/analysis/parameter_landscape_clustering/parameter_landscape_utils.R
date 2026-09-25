#!/usr/bin/env Rscript

# Historical loader retained for callers that source this path. The canonical
# analysis-only API is kept in parameter_landscape_analysis_utils.R.
.o2pl_compat_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "parameter_landscape_utils.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(.o2pl_compat_dir, "parameter_landscape_analysis_utils.R"), local = environment(), chdir = TRUE)
rm(.o2pl_compat_dir)
