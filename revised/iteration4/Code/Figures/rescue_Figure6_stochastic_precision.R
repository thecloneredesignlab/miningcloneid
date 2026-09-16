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
source(file.path(script_dir, "util", "analysis", "figure6_invitro_passage_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_full_range_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_precision_rescue.R"))

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = character()) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}
as_flag <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")

base_run <- option_value("base-run")
run_id <- option_value("run-id")
if (!nzchar(base_run) || !nzchar(run_id)) {
  stop("Required arguments: --base-run=ABSOLUTE_PATH --run-id=ID")
}
workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
f6pr_data(
  workspace_root = workspace_root,
  base_run_root = base_run,
  run_id = run_id,
  n_core = as.integer(option_value("n-core", "16")),
  replicates = as.integer(option_value("replicates", "200")),
  master_seed = as.integer(option_value("master-seed", "20260907")),
  publish_current = as_flag(option_value("publish-current", "TRUE"))
)
