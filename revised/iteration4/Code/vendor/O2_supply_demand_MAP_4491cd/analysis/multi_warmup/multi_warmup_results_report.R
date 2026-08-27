#!/usr/bin/env Rscript

# Deprecated compatibility orchestrator for the historical multi-warmup report.

.o2mw_report_legacy_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "multi_warmup_results_report.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_report_legacy_root <- normalizePath(file.path(.o2mw_report_legacy_dir, "..", ".."), mustWork = TRUE)
sys.source(file.path(.o2mw_report_legacy_root, "runner", "multi_warmup", "render_multi_warmup_results.R"),
           envir = environment(), chdir = TRUE)
if (sys.nframe() == 0L) render_multi_warmup_results()
