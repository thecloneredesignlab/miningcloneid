#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
source(file.path(SCRIPT_DIR, "..", "util", "curve_classification_utils.R"))
