#!/usr/bin/env Rscript

Sys.setenv(
  KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("This entry point must be run with Rscript.")
  dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
})
source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure6_context_extension.R"))
source(file.path(script_dir, "util", "analysis", "figure6_finite_time_q10.R"))

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = character()) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}
as_flag <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
f6ft_data(
  workspace_root = workspace_root,
  n_core = as.integer(option_value("n-core", "1")),
  run_id = option_value("run-id", f6ft_resolve_run_id()),
  smoke = as_flag(option_value("smoke", "FALSE")),
  publish_current = as_flag(option_value("publish-current", "TRUE")),
  compute_diagnostics = as_flag(option_value("compute-diagnostics", "TRUE"))
)
