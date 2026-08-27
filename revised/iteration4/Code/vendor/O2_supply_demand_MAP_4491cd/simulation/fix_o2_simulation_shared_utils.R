#!/usr/bin/env Rscript

# Backward-compatible loader for the canonical fixed-O2 utilities.
.legacy_fixed_o2_wrapper_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frame_files[basename(frame_files) == "fix_o2_simulation_shared_utils.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
    if (basename(path) == "fix_o2_simulation_shared_utils.R") return(dirname(path))
  }
  normalizePath(getwd(), mustWork = FALSE)
})
.legacy_fixed_o2_canonical <- normalizePath(
  file.path(.legacy_fixed_o2_wrapper_dir, "..", "util", "o2_supply_demand_map_fixed_o2_utils.R"),
  mustWork = TRUE
)
source(.legacy_fixed_o2_canonical, local = environment(), chdir = TRUE)
rm(
  bpf_default_fixed_o2_grid,
  bpf_default_dense_attractor_o2_grid,
  bpf_o2_key,
  bpf_o2_slug,
  fixed_o2_shared_utils_dir
)
rm(.legacy_fixed_o2_wrapper_dir, .legacy_fixed_o2_canonical)
