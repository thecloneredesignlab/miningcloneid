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
    DATA_ROOT, "Figure5", "figure5f_selected_pair_inputs.tsv"
  )
  require_files(selection_path, "Figure 5 selected primary-family pairs")
  selection <- utils::read.delim(
    selection_path, check.names = FALSE, stringsAsFactors = FALSE
  )
  required_selection <- c(
    "family", "warmup_label", "selected_seed", "invitro_seed",
    "selected_for_figure5f"
  )
  if (!all(required_selection %in% names(selection))) {
    stop("Figure 5 selected-pair table lacks required fields.")
  }
  joint_selection <- selection[
    as.logical(selection$selected_for_figure5f), , drop = FALSE
  ]
  family_order <- c("C01", "C02", "C03")
  joint_selection <- joint_selection[
    order(match(joint_selection$family, family_order)), , drop = FALSE
  ]
  expected_selection <- c(
    C01 = "tsne_vi_seed366_C01Sc01_vt_seed10",
    C02 = "tsne_vi_seed25_C02Sc01_vt_seed10",
    C03 = "tsne_vi_seed311_C03Sc02_vt_seed10"
  )
  observed_selection <- setNames(
    joint_selection$warmup_label, joint_selection$family
  )
  if (nrow(joint_selection) != 3L ||
      !identical(joint_selection$family, family_order) ||
      !identical(observed_selection[family_order], expected_selection) ||
      any(joint_selection$invitro_seed != "seed10")) {
    stop(
      "Supplementary Figure 5-2 selected pairs do not match the approved ",
      "C01/C02/C03 primary-family inputs."
    )
  }

  joint_output_root <- file.path(output_root, "joint")
  dir.create(joint_output_root, recursive = TRUE, showWarnings = FALSE)
  existing_bundles <- list.dirs(
    joint_output_root, recursive = FALSE, full.names = TRUE
  )
  obsolete_bundles <- existing_bundles[
    !basename(existing_bundles) %in% joint_selection$warmup_label
  ]
  if (length(obsolete_bundles)) {
    unlink(obsolete_bundles, recursive = TRUE, force = TRUE)
    if (any(dir.exists(obsolete_bundles))) {
      stop("Could not remove obsolete secondary-group diagnostic bundles.")
    }
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
    role = c(
      "selected primary-family pair table",
      rep("joint optimizer diagnostic input", length(sources))
    ),
    source = c(selection_path, sources),
    local_file = NA_character_,
    source_md5 = c(
      unname(tools::md5sum(selection_path)),
      unname(tools::md5sum(sources))
    ),
    local_md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  write_data_contract("Supp_Figure5_2", contract)
  invisible(TRUE)
}

if (sys.nframe() == 0L) data_Supp_Figure5_2()
