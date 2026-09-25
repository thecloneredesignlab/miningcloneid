#!/usr/bin/env Rscript

# Deprecated compatibility entrypoint; report assembly lives under report/.
.o2pl_compat_dir <- local({
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = FALSE))
})
.o2pl_root <- normalizePath(file.path(.o2pl_compat_dir, "..", ".."), mustWork = TRUE)
.o2pl_report <- file.path(.o2pl_root, "report", "parameter_landscape", "mode_parameter_contribution_report.R")
.o2pl_status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", .o2pl_report, commandArgs(trailingOnly = TRUE)))
if (!identical(.o2pl_status, 0L)) stop("Canonical mode report failed with status ", .o2pl_status, call. = FALSE)
