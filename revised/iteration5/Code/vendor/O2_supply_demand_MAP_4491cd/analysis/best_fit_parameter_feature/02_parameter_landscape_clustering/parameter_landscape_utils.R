#!/usr/bin/env Rscript
# Deprecated compatibility loader for callers that source the former duplicate utility.
Sys.setenv(O2PL_DEFAULT_RESULT_LAYOUT = "legacy_bpf")
.o2pl_frame_file <- sys.frame(1)$ofile
if (is.null(.o2pl_frame_file) || !nzchar(.o2pl_frame_file)) {
  .o2pl_frame_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
}
.o2pl_here <- dirname(normalizePath(.o2pl_frame_file, mustWork = FALSE))
.o2pl_root <- normalizePath(file.path(.o2pl_here, "..", "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_root, "analysis", "parameter_landscape_clustering", "parameter_landscape_utils.R"), local = parent.frame(), chdir = TRUE)
