#!/usr/bin/env Rscript

# Deprecated compatibility orchestrator for the historical landscape-pair CLI.

.o2mw_landscape_legacy_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "build_multi_warmup_pairs_from_landscape_subclusters.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_landscape_legacy_root <- normalizePath(file.path(.o2mw_landscape_legacy_dir, "..", ".."), mustWork = TRUE)
sys.source(file.path(.o2mw_landscape_legacy_root, "runner", "multi_warmup", "run_multi_warmup_landscape_pipeline.R"),
           envir = environment(), chdir = TRUE)
if (sys.nframe() == 0L) run_multi_warmup_landscape_pipeline()
