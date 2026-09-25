#!/usr/bin/env Rscript

# Deprecated compatibility orchestrator.  The historical all-in-one entrypoint
# now delegates to runner/fixed_o2_eigen/, which sequences simulation, analysis,
# visualization, and report stages through materialized file contracts.

.fixo2ea_legacy_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "fixo2_eigen_attractor_clustering_runner.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
.fixo2ea_legacy_workflow_root <- normalizePath(file.path(.fixo2ea_legacy_runner_dir, "..", "..", ".."), mustWork = TRUE)
.fixo2ea_pipeline <- file.path(
  .fixo2ea_legacy_workflow_root, "runner", "fixed_o2_eigen",
  "run_fixo2_eigen_attractor_pipeline.R"
)
sys.source(.fixo2ea_pipeline, envir = environment(), chdir = TRUE)
if (sys.nframe() == 0L) fixo2ea_main()
