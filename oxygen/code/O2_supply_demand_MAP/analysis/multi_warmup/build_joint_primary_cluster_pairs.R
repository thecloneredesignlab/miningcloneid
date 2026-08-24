#!/usr/bin/env Rscript

# Canonical joint-pair entrypoint: in-vivo primary clusters to one in-vitro anchor.

.o2_joint_primary_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)))
  }
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "build_joint_primary_cluster_pairs.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  normalizePath(getwd(), mustWork = TRUE)
})
.o2_joint_primary_analysis <- new.env(parent = globalenv())
source(
  file.path(.o2_joint_primary_dir, "build_multi_warmup_landscape_tables.R"),
  local = .o2_joint_primary_analysis,
  chdir = TRUE
)
dispatch_main <- .o2_joint_primary_analysis$dispatch_main
if (sys.nframe() == 0L) dispatch_main()
