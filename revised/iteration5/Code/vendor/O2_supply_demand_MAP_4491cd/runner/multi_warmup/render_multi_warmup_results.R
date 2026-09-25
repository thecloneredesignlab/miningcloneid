#!/usr/bin/env Rscript

# Visualization-plus-report orchestrator for already materialized multi-warmup
# analysis products.

.o2mw_report_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_multi_warmup_results.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  workflow_frames <- frames[
    grepl("/O2_supply_demand_MAP/", frames, fixed = TRUE)
  ]
  if (length(workflow_frames)) {
    root <- dirname(workflow_frames[[length(workflow_frames)]])
    while (!identical(basename(root), "O2_supply_demand_MAP")) {
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    if (identical(basename(root), "O2_supply_demand_MAP")) {
      return(file.path(root, "runner", "multi_warmup"))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_report_workflow_root <- normalizePath(file.path(.o2mw_report_runner_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_report_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

.o2mw_report_vis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_report_workflow_root, "vis", "multi_warmup", "render_multi_warmup_report_figures.R"),
           envir = .o2mw_report_vis, chdir = TRUE)
.o2mw_report_env <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_report_workflow_root, "report", "multi_warmup", "render_multi_warmup_results_report.R"),
           envir = .o2mw_report_env, chdir = TRUE)

render_multi_warmup_results <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root <- normalizePath(path.expand(as_chr(argv$multi_warmup_root, as_chr(argv$out_dir))), mustWork = TRUE)
  .o2mw_report_vis$render_multi_warmup_report_figures(root)
  .o2mw_report_env$main(argv)
  invisible(TRUE)
}

main <- render_multi_warmup_results
if (sys.nframe() == 0L) render_multi_warmup_results()
