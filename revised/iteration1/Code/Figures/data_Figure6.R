#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

data_Figure6 <- function(n_core = 8L, rebuild = FALSE) {
  data_dir <- file.path(DATA_ROOT, "Figure6")
  precompute_figure_dir <- file.path(data_dir, "precompute_panels")
  precompute_output_dir <- file.path(data_dir, "precompute_assembled")
  run_oxygen_response_pipeline(
    args = c(
      paste0("--joint_root=", JOINT_RESULT_ROOT),
      paste0("--workspace_dir=", data_dir),
      paste0("--data_dir=", data_dir),
      paste0("--figure_dir=", precompute_figure_dir),
      paste0("--output_dir=", precompute_output_dir),
      paste0("--workers=", as.integer(n_core)),
      paste0("--rebuild=", if (isTRUE(rebuild)) "TRUE" else "FALSE")
    ),
    env = paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT)
  )
  required <- file.path(data_dir, c(
    "response_class_smoothed_curves.tsv",
    "response_class_curve_class_by_seed.tsv",
    "response_class_class_counts.tsv",
    "pair_surface_best_joint_seed_summary.tsv",
    "pair_surface_best_trajectory.tsv",
    "pair_surface_complete_surface.tsv",
    "selection_diagnostic_selection_data.tsv",
    "response_workflow_validation.tsv",
    "response_workflow_source_file_provenance.tsv"
  ))
  require_files(required, "Figure 6 analysis intermediate")
  provenance <- utils::read.delim(
    file.path(data_dir, "response_workflow_source_file_provenance.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  scientific <- provenance[
    startsWith(
      provenance$source_file,
      paste0(normalizePath(JOINT_RESULT_ROOT), .Platform$file.sep)
    ),
    ,
    drop = FALSE
  ]
  pair_dirs <- sort(list.dirs(
    JOINT_RESULT_ROOT, recursive = FALSE, full.names = TRUE
  ))
  pair_dirs <- pair_dirs[grepl("^fit_joint_", basename(pair_dirs))]
  selection_inputs <- unlist(lapply(pair_dirs, function(pair_dir) {
    seed_dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = TRUE)
    seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
    file.path(seed_dirs, "fit_summary.tsv")
  }), use.names = FALSE)
  require_files(
    selection_inputs, "Figure 6 joint winner-selection input"
  )
  contract <- data.frame(
    role = c(
      scientific$role,
      rep("Figure 6B joint winner-selection fit summary", length(selection_inputs))
    ),
    source = c(scientific$source_file, selection_inputs),
    local_file = rep(
      NA_character_, nrow(scientific) + length(selection_inputs)
    ),
    source_md5 = c(
      scientific$md5, unname(tools::md5sum(selection_inputs))
    ),
    local_md5 = rep(
      NA_character_, nrow(scientific) + length(selection_inputs)
    ),
    stringsAsFactors = FALSE
  )
  write_data_contract("Figure6", contract)
  invisible(contract)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  n_core_arg <- sub("^--n-core=", "", args[grepl("^--n-core=", args)])
  rebuild_arg <- sub("^--rebuild=", "", args[grepl("^--rebuild=", args)])
  data_Figure6(
    n_core = if (length(n_core_arg)) as.integer(n_core_arg[[1L]]) else 8L,
    rebuild = if (length(rebuild_arg)) {
      as_boolean(rebuild_arg[[1L]])
    } else {
      FALSE
    }
  )
}
