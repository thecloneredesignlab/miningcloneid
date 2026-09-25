#!/usr/bin/env Rscript

# Map parameter-level ClassA/B/C stability onto audited biological processes.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

summarize_soft_coupling_processes <- function(out_dir) {
  class_summary <- o2jca_read_tsv(file.path(out_dir, "between_pair_class_summary.tsv"))
  stability <- o2jca_read_tsv(file.path(out_dir, "between_pair_parameter_stability.tsv"))
  process_map <- o2jca_read_tsv(file.path(out_dir, "parameter_process_map.tsv"))
  class_mapped <- merge(class_summary, process_map, by = "parameter", all.x = TRUE, sort = FALSE)
  stability_mapped <- merge(stability, process_map, by = "parameter", all.x = TRUE, sort = FALSE)
  groups <- o2jca_group_split(class_mapped, c("primary_process", "ratio_class"))
  process_summary <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      primary_process = group$primary_process[[1L]], ratio_class = group$ratio_class[[1L]],
      n_parameters = nrow(group),
      mean_pair_balanced_proportion = mean(group$pair_balanced_mean_proportion, na.rm = TRUE),
      median_pair_balanced_proportion = stats::median(group$pair_balanced_mean_proportion, na.rm = TRUE),
      n_parameters_all_pairs_stable90 = sum(group$n_pairs_stable90 == group$n_pairs, na.rm = TRUE),
      n_parameters_all_pairs_strict = sum(group$n_pairs_strict_intersection == group$n_pairs, na.rm = TRUE),
      parameters = paste(group$parameter[order(-group$pair_balanced_mean_proportion)], collapse = ","),
      stringsAsFactors = FALSE
    )
  }))
  report_summary <- stability_mapped[order(-stability_mapped$cross_pair_dominant_fraction, stability_mapped$parameter), c(
    "parameter", "primary_process", "cross_pair_dominant_class", "cross_pair_dominant_fraction",
    "all_pairs_same_dominant_class", "pair_median_direction_consistency", "intraclass_correlation",
    "shared_invitro_anchor", "invitro_anchor_seed"
  )]
  report_summary$interpretation_scope <- "model-derived stability across in-vivo warm-up pairs conditional on a shared in-vitro seed10 anchor"
  files <- c(
    o2jca_write_tsv(class_mapped, file.path(out_dir, "parameter_class_process_summary.tsv")),
    o2jca_write_tsv(stability_mapped, file.path(out_dir, "parameter_process_consensus.tsv")),
    o2jca_write_tsv(process_summary, file.path(out_dir, "process_class_summary.tsv")),
    o2jca_write_tsv(report_summary, paste0(out_dir, "/soft_coupling_report_summary.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "process_analysis_manifest.tsv"))
  invisible(process_summary)
}

main <- function() {
  args <- o2jca_parse_args(); out_dir <- args$out_dir %||% stop("--out_dir is required", call. = FALSE)
  summarize_soft_coupling_processes(out_dir)
}

if (identical(environment(), globalenv())) main()
