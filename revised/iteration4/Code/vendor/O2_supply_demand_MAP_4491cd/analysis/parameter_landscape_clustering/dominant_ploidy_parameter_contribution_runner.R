#!/usr/bin/env Rscript

# Deprecated compatibility entrypoint; orchestration belongs in runner/.
.o2pl_compat_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else x$ofile, character(1)))
  own <- frames[basename(frames) == "dominant_ploidy_parameter_contribution_runner.R"]
  if (length(own)) dirname(normalizePath(own[[length(own)]], mustWork = FALSE)) else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = FALSE))
  }
})
.o2pl_root <- normalizePath(file.path(.o2pl_compat_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_root, "runner", "parameter_landscape", "run_parameter_landscape.R"), local = environment(), chdir = TRUE)
.o2pl_argv <- parse_args(commandArgs(trailingOnly = TRUE))
.o2pl_argv$run_parts <- "dominant_ploidy_contribution"
o2pl_parameter_landscape_runner_main(.o2pl_argv)
