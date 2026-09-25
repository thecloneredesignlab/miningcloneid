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
source(file.path(script_dir, "util", "analysis", "figure6_stochastic_passage.R"))
source(file.path(script_dir, "util", "analysis", "figure6_invitro_passage_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_full_range_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_net_growth_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_net_growth_full_range_q10.R"))
source(file.path(script_dir, "util", "analysis", "figure6_growth_permissive.R"))

workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1L]])
}
phase <- option_value("phase", "all")
smoke <- tolower(option_value("smoke", "false")) %in% c("1", "true", "yes")
if (!phase %in% c("all", "prepare", "task", "finalize")) {
  stop("Use --phase=all, --phase=prepare, --phase=task, or --phase=finalize.")
}
if (phase == "all") {
  result <- f6gp_run(workspace_root, n_core = as.integer(option_value("n-core", "1")),
                     smoke = smoke)
  message("Wrote day-1000 positive-growth curves: ", result$output$curve)
} else if (phase == "prepare") {
  result <- f6gp_prepare(workspace_root, smoke = smoke)
  message("Prepared tasks: ", result$output$task)
} else if (phase == "task") {
  setup <- f6gp_prepare(workspace_root, smoke = smoke)
  task_index <- as.integer(option_value("task-index", NA_character_))
  if (length(task_index) != 1L || is.na(task_index) || task_index < 1L ||
      task_index > nrow(setup$tasks)) stop("Invalid --task-index.")
  f6r_require_packages(c("Matrix", "Rcpp"))
  f6r_load_response_engine(setup$paths)
  f6ng_load_propagator(setup$paths)
  objective_bundle <- f6r_objective_selection(setup$paths)
  contexts <- lapply(
    unique(setup$source$endpoints$pair_id), f6r_pair_model_context,
    selected = objective_bundle$selected, paths = setup$paths
  )
  names(contexts) <- unique(setup$source$endpoints$pair_id)
  qc <- f6gp_compute_task(
    setup$tasks[task_index, , drop = FALSE], setup$source,
    objective_bundle, contexts, setup$paths, setup$fingerprint,
    smoke = smoke
  )
  message("Completed task ", task_index, ": positive counts ",
          qc$n_positive_2N, ", ", qc$n_positive_4N)
} else {
  result <- f6gp_finalize(workspace_root, smoke = smoke)
  message("Wrote validated curves: ", result$output$curve)
}
