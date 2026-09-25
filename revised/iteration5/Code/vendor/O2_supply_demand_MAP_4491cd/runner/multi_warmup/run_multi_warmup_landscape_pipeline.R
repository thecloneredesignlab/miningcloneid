#!/usr/bin/env Rscript

# Cross-stage orchestrator for landscape-based multi-warmup pair selection.

.o2mw_landscape_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "run_multi_warmup_landscape_pipeline.R"]
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
.o2mw_landscape_workflow_root <- normalizePath(file.path(.o2mw_landscape_runner_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_landscape_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

.o2mw_landscape_analysis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_landscape_workflow_root, "analysis", "multi_warmup", "build_multi_warmup_landscape_tables.R"),
           envir = .o2mw_landscape_analysis, chdir = TRUE)
.o2mw_landscape_vis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_landscape_workflow_root, "vis", "multi_warmup", "render_multi_warmup_landscape_figures.R"),
           envir = .o2mw_landscape_vis, chdir = TRUE)

run_multi_warmup_landscape_pipeline <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  .o2mw_landscape_analysis$dispatch_main(argv)
  stage <- tolower(gsub("-", "_", as_chr(argv$stage %||% argv$mode, "build_pairs"), fixed = TRUE))
  render_stages <- c("build_pairs", "pair", "prepare_landscape", "prepare")
  if (stage %in% render_stages && as_bool(argv$render_figures, TRUE)) {
    out_dir <- normalizePath(path.expand(as_chr(argv$out_dir)), mustWork = TRUE)
    analysis_root <- normalizePath(
      path.expand(as_chr(argv$analysis_root, .o2mw_landscape_analysis$landscape_analysis_root(out_dir))),
      mustWork = TRUE
    )
    .o2mw_landscape_vis$render_multi_warmup_landscape_figures(analysis_root)
  }
  invisible(TRUE)
}

dispatch_main <- run_multi_warmup_landscape_pipeline
if (sys.nframe() == 0L) run_multi_warmup_landscape_pipeline()
