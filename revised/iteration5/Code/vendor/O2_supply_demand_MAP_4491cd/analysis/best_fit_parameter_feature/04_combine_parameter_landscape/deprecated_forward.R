#!/usr/bin/env Rscript

# Deprecated compatibility helper for the pre-stage-3 combined workflow.
.deprecated_combined_file <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "deprecated_forward.R"]
  if (length(own)) own[[length(own)]] else ""
})
deprecated_combined_root <- function() normalizePath(file.path(dirname(.deprecated_combined_file), "..", "..", ".."), mustWork = TRUE)
deprecated_combined_runner <- function() file.path(deprecated_combined_root(), "runner", "combined_parameter_landscape", "run_combined_parameter_landscape.R")
deprecated_combined_analysis <- function() file.path(deprecated_combined_root(), "analysis", "combined_parameter_landscape", "prepare_combined_parameter_landscape_tables.R")
deprecated_combined_exec <- function(raw, mode = "run", extra = character()) {
  if (!any(grepl("^--out[_-]dir=", raw))) {
    embedding <- sub("^--embedding[_-]dir=", "", grep("^--embedding[_-]dir=", raw, value = TRUE))
    if (length(embedding)) raw <- c(raw, paste0("--out_dir=", embedding[[1L]]))
  }
  if (!any(grepl("^--mode=", raw))) raw <- c(raw, paste0("--mode=", mode))
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", deprecated_combined_runner(), raw, extra))
  if (!identical(status, 0L)) stop("Canonical combined runner failed with status ", status, call. = FALSE)
  invisible(TRUE)
}
deprecated_combined_direct <- function(mode = "run", extra = character()) deprecated_combined_exec(commandArgs(trailingOnly = TRUE), mode, extra)
