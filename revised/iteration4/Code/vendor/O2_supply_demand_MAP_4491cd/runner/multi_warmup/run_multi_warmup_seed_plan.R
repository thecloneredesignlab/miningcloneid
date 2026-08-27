#!/usr/bin/env Rscript

# Cross-stage seed-plan orchestrator: table analysis followed by pure rendering.

.o2mw_seed_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "run_multi_warmup_seed_plan.R"]
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
.o2mw_seed_workflow_root <- normalizePath(file.path(.o2mw_seed_runner_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_seed_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

.o2mw_seed_analysis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_seed_workflow_root, "analysis", "multi_warmup", "build_multi_warmup_seed_plan_tables.R"),
           envir = .o2mw_seed_analysis, chdir = TRUE)
.o2mw_seed_vis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_seed_workflow_root, "vis", "multi_warmup", "render_multi_warmup_seed_plan_figures.R"),
           envir = .o2mw_seed_vis, chdir = TRUE)

run_multi_warmup_seed_plan <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  .o2mw_seed_analysis$main(argv)
  out_dir <- normalizePath(path.expand(as_chr(argv$out_dir)), mustWork = TRUE)
  .o2mw_seed_vis$render_multi_warmup_seed_plan_figures(
    out_dir,
    umap_seed = as_int(argv$umap_seed %||% argv$multi_warmup_umap_seed, 1L)
  )
  invisible(TRUE)
}

main <- run_multi_warmup_seed_plan
if (sys.nframe() == 0L) run_multi_warmup_seed_plan()
