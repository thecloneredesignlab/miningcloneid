#!/usr/bin/env Rscript

# Backward-compatible analysis loader for the canonical fixed-O2 utilities.
.bpf_fixed_o2_wrapper_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frame_files[basename(frame_files) == "fixed_o2_shared_utils.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
    if (basename(path) == "fixed_o2_shared_utils.R") return(dirname(path))
  }
  normalizePath(getwd(), mustWork = FALSE)
})
.bpf_fixed_o2_canonical <- normalizePath(
  file.path(.bpf_fixed_o2_wrapper_dir, "..", "..", "..", "util", "o2_supply_demand_map_fixed_o2_utils.R"),
  mustWork = TRUE
)
source(.bpf_fixed_o2_canonical, local = environment(), chdir = TRUE)
rm(.bpf_fixed_o2_wrapper_dir, .bpf_fixed_o2_canonical)
