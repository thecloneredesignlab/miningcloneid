#!/usr/bin/env Rscript

# Deprecated compatibility orchestrator for the historical seed-plan CLI.

.o2mw_seed_legacy_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "build_multi_warmup_seed_plan.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_seed_workflow_root <- normalizePath(file.path(.o2mw_seed_legacy_dir, "..", ".."), mustWork = TRUE)
sys.source(file.path(.o2mw_seed_workflow_root, "runner", "multi_warmup", "run_multi_warmup_seed_plan.R"),
           envir = environment(), chdir = TRUE)
if (sys.nframe() == 0L) run_multi_warmup_seed_plan()
