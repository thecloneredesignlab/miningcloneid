#!/usr/bin/env Rscript

# Deprecated name retained for compatibility; the canonical analysis consumes
# only simulation-layer materialized tables.
.o2pl_compat_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else x$ofile, character(1)))
  own <- frames[basename(frames) == "mode_parameter_contribution_analysis.R"]
  if (length(own)) dirname(normalizePath(own[[length(own)]], mustWork = FALSE)) else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = FALSE))
  }
})
source(file.path(.o2pl_compat_dir, "parameter_contribution_analysis.R"), local = environment(), chdir = TRUE)
o2pl_contribution_main(parse_args(commandArgs(trailingOnly = TRUE)), forced_target = "mode")
