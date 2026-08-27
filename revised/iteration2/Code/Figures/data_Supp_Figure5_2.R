#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))
source(resolve_utility_file("analysis/optimizer_diagnostics.R"))

data_Supp_Figure5_2 <- function() {
  output_root <- file.path(DATA_ROOT, "Supp_Figure5_2")
  selection_path <- file.path(
    DATA_ROOT, "Figure5", "figure5_frozen_inputs",
    "selected_results.tsv"
  )
  require_files(selection_path, "Figure 5 selected-results table")
  selection <- utils::read.delim(
    selection_path, check.names = FALSE, stringsAsFactors = FALSE
  )
  joint_selection <- selection[
    selection$record_type == "joint_pair_best", , drop = FALSE
  ]
  if (nrow(joint_selection) != 6L) {
    stop("Supplementary Figure 5-2 requires six selected joint winners.")
  }

  joint_source_files <- character()
  boundary_source_files <- character()
  for (i in seq_len(nrow(joint_selection))) {
    row <- joint_selection[i, , drop = FALSE]
    pair_id <- paste0("fit_joint_", row$warmup_label[[1L]])
    pair_dir <- file.path(JOINT_RESULT_ROOT, pair_id)
    metrics <- optimizer_collect_seed_metrics(pair_dir, "objective")
    selected_summary <- optimizer_write_bundle(
      metrics,
      file.path(output_root, "joint", row$warmup_label[[1L]]),
      "joint fitting",
      row$selected_seed[[1L]],
      pair_dir
    )
    joint_source_files <- c(joint_source_files, metrics$source_file)
    boundary_source_files <- c(
      boundary_source_files,
      attr(selected_summary, "source_files")
    )
  }

  sources <- unique(c(joint_source_files, boundary_source_files))
  contract <- data.frame(
    role = "joint optimizer diagnostic input",
    source = sources,
    local_file = NA_character_,
    source_md5 = unname(tools::md5sum(sources)),
    local_md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  write_data_contract("Supp_Figure5_2", contract)
  invisible(TRUE)
}

if (sys.nframe() == 0L) data_Supp_Figure5_2()
