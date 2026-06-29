#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
source(file.path(SCRIPT_DIR, "dominant_ploidy_parameter_contribution_analysis.R"))

argv <- parse_args(commandArgs(trailingOnly = TRUE))
argv$mode_contribution_target <- "dominant_ploidy"
run_dominant_ploidy_parameter_contribution_analysis(argv)
