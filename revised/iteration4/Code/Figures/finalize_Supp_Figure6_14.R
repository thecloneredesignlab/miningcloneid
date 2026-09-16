#!/usr/bin/env Rscript

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

args <- commandArgs(trailingOnly = TRUE)
hit <- args[startsWith(args, "--run-id=")]
if (length(hit) != 1L) stop("Exactly one --run-id=ID argument is required.")
run_id <- sub("^--run-id=", "", hit[[1L]])
workspace_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
f6ng_finalize_existing(workspace_root, run_id)
