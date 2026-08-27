#!/usr/bin/env Rscript

# Deprecated compatibility helper for the former duplicate workflow tree.
bpf_parameter_landscape_forward <- function(target = basename(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))) {
  here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])), mustWork = TRUE)
  root <- normalizePath(file.path(here, "..", "..", ".."), mustWork = TRUE)
  canonical <- file.path(root, "analysis", "parameter_landscape_clustering", target)
  if (!file.exists(canonical)) stop("Missing canonical parameter-landscape entrypoint: ", canonical, call. = FALSE)
  env <- c(paste0("O2PL_DEFAULT_RESULT_LAYOUT=legacy_bpf"), paste0("O2PL_WORKFLOW_ROOT=", root))
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", canonical, commandArgs(trailingOnly = TRUE)), env = env)
  if (!identical(status, 0L)) stop("Canonical parameter-landscape entrypoint failed with status ", status, call. = FALSE)
  invisible(status)
}
