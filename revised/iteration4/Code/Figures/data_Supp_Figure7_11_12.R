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
source(file.path(script_dir, "util", "analysis", "figure7_robustness.R"))
source(file.path(script_dir, "util", "analysis", "figure7_context_extension.R"))
source(file.path(script_dir, "util", "analysis", "figure7_finite_time_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure7_invitro_passage_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure7_extended_time_o2.R"))

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = character()) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}
as_flag <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

mode <- option_value("mode", "all")
if (!mode %in% c(f7e_modes(), "all")) {
  stop("--mode must be passage, continuous, or all.")
}
n_core <- as.integer(option_value("n-core", "1"))
run_id <- option_value("run-id", format(Sys.time(), "%Y%m%d_%H%M%S"))
smoke <- as_flag(option_value("smoke", "FALSE"))
publish_current <- as_flag(option_value("publish-current", "TRUE"))
workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

modes <- if (mode == "all") f7e_modes() else mode
for (selected_mode in modes) {
  f7e_data(
    workspace_root = workspace_root, mode = selected_mode,
    n_core = n_core, run_id = paste0(run_id, "_", selected_mode),
    smoke = smoke, publish_current = publish_current && !smoke
  )
}
