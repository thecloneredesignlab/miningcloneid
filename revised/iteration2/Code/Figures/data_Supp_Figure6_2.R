#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "analysis", "figure6_robustness.R"))

data_Supp_Figure6_2 <- function() {
  workspace_root <- normalizePath(
    file.path(script_dir, "..", ".."), mustWork = TRUE
  )
  paths <- f6r_paths(workspace_root)
  required <- file.path(paths$figure6, c(
    "cluster_k_selection.tsv",
    "cluster_k_resample_selection_frequency.tsv",
    "cluster_stability.tsv",
    "joint_seed_acceptance.tsv",
    "joint_seed_claim_robustness.tsv",
    "joint_seed_cutoff_consistency.tsv",
    "joint_unique_parameter_endpoint_counts.tsv",
    "joint_seed_vs_unique_parameter_robustness.tsv",
    "joint_seed_vs_unique_parameter_surface_comparison.tsv",
    "joint_seed_vs_unique_parameter_trajectory_comparison.tsv",
    "joint_seed_response_delta_summary.tsv",
    "joint_seed_robustness_audit_summary.tsv"
  ))
  f6r_require_files(required, "Supplementary Figure 6-2 analytical input")
  dir.create(paths$supp6_2, recursive = TRUE, showWarnings = FALSE)
  contract <- data.frame(
    figure_id = "Supp_Figure6_2",
    input_file = normalizePath(required, mustWork = TRUE),
    size_bytes = file.info(required)$size,
    md5 = unname(tools::md5sum(required)),
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(contract, file.path(paths$supp6_2, "data_contract.tsv"))
  invisible(required)
}

if (sys.nframe() == 0L) data_Supp_Figure6_2()
