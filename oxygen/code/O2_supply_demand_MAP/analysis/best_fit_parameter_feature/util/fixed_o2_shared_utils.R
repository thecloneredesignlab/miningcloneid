#!/usr/bin/env Rscript

bpf_default_fixed_o2_grid <- function() {
  c(0, 0.1, 0.5, 1, 2, 5)
}

bpf_default_dense_attractor_o2_grid <- function() {
  seq(0, 5, by = 0.05)
}

bpf_o2_key <- function(x) {
  vapply(x, function(xx) format(signif(as.numeric(xx), 12), scientific = FALSE, trim = TRUE), character(1))
}

bpf_o2_slug <- function(x) {
  key <- bpf_o2_key(x)
  key <- gsub("-", "minus", key, fixed = TRUE)
  key <- gsub("[^0-9A-Za-z]+", "p", key)
  key <- gsub("^p+|p+$", "", key)
  ifelse(nzchar(key), key, "NA")
}

fixed_o2_shared_utils_dir <- function() {
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files)) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

FIXO2_SIMULATION_SHARED_UTILS <- normalizePath(
  file.path(fixed_o2_shared_utils_dir(), "..", "..", "..", "simulation", "fix_o2_simulation_shared_utils.R"),
  mustWork = TRUE
)
source(FIXO2_SIMULATION_SHARED_UTILS, local = environment())
rm(FIXO2_SIMULATION_SHARED_UTILS)
