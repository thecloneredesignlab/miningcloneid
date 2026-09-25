#!/usr/bin/env Rscript

# Copy only the eight diagnostic traces that are outside the rescued condition.
# Each source/destination pair is checked byte-for-byte and recorded.

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
base_run <- option_value("base-run")
run_id <- option_value("run-id")
if (!nzchar(base_run) || !nzchar(run_id)) {
  stop("Required arguments: --base-run=ABSOLUTE_PATH --run-id=ID")
}

workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
paths <- f6r_paths(workspace_root)
base_run <- f6pr_assert_base_root(paths, base_run)
run <- f6g_paths(paths, run_id = run_id, create = FALSE)
run_root <- normalizePath(run$run_root, mustWork = TRUE)
base_fingerprint <- readRDS(file.path(
  base_run, "task_cache", "full_range_task_0337.rds"
))$fingerprint
trace_files <- sort(list.files(
  base_run, pattern = "^stochastic_trace_.*[.]rds$", full.names = TRUE
))
if (length(trace_files) != 8L) {
  stop("Expected exactly eight immutable stochastic diagnostic traces.")
}

records <- lapply(trace_files, function(source_path) {
  object <- readRDS(source_path)
  if (!identical(object$fingerprint, base_fingerprint) ||
      !object$endpoint$pair_label[[1L]] %in% f6ft_family_levels() ||
      !object$O2_pct %in% c(0.5, 20) ||
      !object$p_misseg %in% c(0.005, 0.3) ||
      (abs(object$O2_pct - 0.7) < 1e-12 &&
       abs(object$p_misseg - 0.01) < 1e-12)) {
    stop("Trace is not outside the rescued grid condition: ", source_path)
  }
  destination <- file.path(run_root, basename(source_path))
  source_md5 <- unname(tools::md5sum(source_path))
  if (!file.exists(destination) ||
      !identical(unname(tools::md5sum(destination)), source_md5)) {
    temporary <- paste0(destination, ".tmp-", Sys.getpid())
    on.exit(unlink(temporary), add = TRUE)
    if (!file.copy(source_path, temporary, overwrite = TRUE) ||
        !identical(unname(tools::md5sum(temporary)), source_md5) ||
        !file.rename(temporary, destination)) {
      stop("Failed atomic trace reuse: ", basename(source_path))
    }
  }
  data.frame(
    file = basename(source_path), pair_label = object$endpoint$pair_label[[1L]],
    p_misseg = object$p_misseg, O2_pct = object$O2_pct,
    endpoint_group = object$endpoint$endpoint_group[[1L]],
    source_run = basename(base_run), source_md5 = source_md5,
    destination_md5 = unname(tools::md5sum(destination)),
    outside_rescued_condition = TRUE,
    passed = identical(unname(tools::md5sum(destination)), source_md5),
    stringsAsFactors = FALSE
  )
})
manifest <- do.call(rbind, records)
if (!all(manifest$passed) ||
    !identical(sort(unique(manifest$pair_label)), f6ft_family_levels()) ||
    !isTRUE(all.equal(sort(unique(manifest$O2_pct)), c(0.5, 20))) ||
    !isTRUE(all.equal(sort(unique(manifest$p_misseg)), c(0.005, 0.3)))) {
  stop("Diagnostic trace reuse contract failed.")
}
f6ft_atomic_write_tsv(
  manifest, file.path(run_root, "stochastic_trace_reuse_manifest.tsv")
)
message("Reused and verified eight unaffected stochastic diagnostic traces.")
