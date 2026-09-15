#!/usr/bin/env Rscript

Sys.setenv(
  KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)
options(stringsAsFactors = FALSE, warn = 1)

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
source(file.path(script_dir, "util", "analysis", "figure6_net_growth_q10.R"))
source(file.path(
  script_dir, "util", "analysis", "figure6_net_growth_full_range_q10.R"
))

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = character()) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}

mode <- option_value("mode", "task")
if (!mode %in% c("prewarm", "task")) {
  stop("--mode must be prewarm or task.")
}
workspace_root <- normalizePath(
  option_value("workspace-root", Sys.getenv("FIGURE_WORKSPACE_ROOT")),
  mustWork = TRUE
)
run_id <- f6ft_sanitize_run_id(option_value("run-id"))
if (!nzchar(run_id)) stop("--run-id is required.")

paths <- f6r_paths(workspace_root)
run_paths <- f6ng_paths(paths, run_id = run_id, create = FALSE)
f6r_load_response_engine(paths)
f6g_load_propagator(paths)
propagator <- f6ng_load_propagator(paths)
source_bundle <- f6ng_source_bundle(paths)
objective_bundle <- f6r_objective_selection(paths)
contexts <- lapply(
  unique(source_bundle$endpoints$pair_id), f6r_pair_model_context,
  selected = objective_bundle$selected, paths = paths
)
names(contexts) <- unique(source_bundle$endpoints$pair_id)
fingerprint <- f6ng_fingerprint(
  paths, source_bundle, propagator, smoke = FALSE
)

manifest_path <- file.path(run_paths$run_root, "net_growth_task_manifest.tsv")
f6r_require_files(manifest_path, "full-range net-growth task manifest")
manifest <- f6r_read_tsv(manifest_path)
required <- c("task_id", "cache_path")
if (nrow(manifest) != 4020L || !all(required %in% names(manifest))) {
  stop("Unexpected full-range task manifest: ", manifest_path)
}

if (identical(mode, "prewarm")) {
  cat("figure6_net_growth_array_prewarm_ok\n")
  cat("run_id=", run_id, "\n", sep = "")
  cat("task_count=", nrow(manifest), "\n", sep = "")
  cat("fingerprint=", fingerprint, "\n", sep = "")
  quit(save = "no", status = 0L)
}

task_file <- normalizePath(option_value("task-file"), mustWork = TRUE)
array_index <- suppressWarnings(as.integer(option_value(
  "array-index", Sys.getenv("SLURM_ARRAY_TASK_ID")
)))
if (is.na(array_index) || array_index < 1L) {
  stop("A positive --array-index or SLURM_ARRAY_TASK_ID is required.")
}
selected <- f6r_read_tsv(task_file)
if (!all(c("array_index", "task_id") %in% names(selected))) {
  stop("Malformed missing-task table: ", task_file)
}
row <- selected[selected$array_index == array_index, , drop = FALSE]
if (nrow(row) != 1L) {
  stop("Array index does not resolve exactly one missing task: ", array_index)
}
task <- manifest[manifest$task_id == row$task_id[[1L]], , drop = FALSE]
if (nrow(task) != 1L) {
  stop("Task id does not resolve exactly one manifest row: ", row$task_id[[1L]])
}
expected_cache <- normalizePath(run_paths$cache, mustWork = TRUE)
observed_cache <- normalizePath(dirname(task$cache_path[[1L]]), mustWork = TRUE)
if (!identical(expected_cache, observed_cache)) {
  stop("Task cache escaped the selected run root: ", task$cache_path[[1L]])
}

started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
qc <- f6ng_compute_task(
  task, source_bundle, objective_bundle, contexts, paths, run_paths,
  fingerprint, smoke = FALSE
)
if (!isTRUE(qc$passed[[1L]])) {
  stop("Task completed but failed QC: ", task$task_id[[1L]])
}

audit_dir <- option_value("audit-dir", "")
if (nzchar(audit_dir)) {
  result_dir <- file.path(audit_dir, "task_results")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  result <- data.frame(
    array_index = array_index,
    task_id = task$task_id[[1L]],
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = ""),
    host = unname(Sys.info()[["nodename"]]),
    started_at = started_at,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    cache_path = normalizePath(task$cache_path[[1L]], mustWork = TRUE),
    passed = TRUE,
    stringsAsFactors = FALSE
  )
  f6ft_atomic_write_tsv(
    result,
    file.path(result_dir, sprintf("array_%04d_%s.tsv", array_index, task$task_id[[1L]]))
  )
}
cat(
  "figure6_net_growth_array_task_complete array_index=", array_index,
  " task_id=", task$task_id[[1L]], "\n", sep = ""
)
