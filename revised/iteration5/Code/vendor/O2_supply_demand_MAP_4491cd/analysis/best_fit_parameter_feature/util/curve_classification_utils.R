#!/usr/bin/env Rscript

# Backward-compatible loader. Canonical implementation lives in workflow util/.
.curve_wrapper_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frame_files[basename(frame_files) == "curve_classification_utils.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
    if (basename(path) == "curve_classification_utils.R") return(dirname(path))
  }
  normalizePath(getwd(), mustWork = FALSE)
})
.curve_canonical <- normalizePath(
  file.path(.curve_wrapper_dir, "..", "..", "..", "util", "o2_supply_demand_map_curve_classification_utils.R"),
  mustWork = TRUE
)
source(.curve_canonical, local = environment(), chdir = TRUE)
rm(.curve_wrapper_dir, .curve_canonical)
