#!/usr/bin/env Rscript

# Deprecated compatibility wrapper.  The canonical cross-workflow runner moved
# to runner/best_fit_parameter_feature/runner.R.
.bpf_old_runner <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE)
.bpf_root <- normalizePath(file.path(dirname(.bpf_old_runner), "..", ".."), mustWork = TRUE)
.bpf_runner <- file.path(.bpf_root, "runner", "best_fit_parameter_feature", "runner.R")
.bpf_status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", .bpf_runner, commandArgs(trailingOnly = TRUE)))
if (!identical(.bpf_status, 0L)) stop("Canonical best-fit-parameter-feature runner failed with status ", .bpf_status, call. = FALSE)
