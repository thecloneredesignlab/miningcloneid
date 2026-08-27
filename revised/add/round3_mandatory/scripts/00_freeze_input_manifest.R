#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)
out_root <- file.path(repo_root, "revised", "add", "round3_mandatory")
dir.create(file.path(out_root, "provenance"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), recursive = TRUE, showWarnings = FALSE)

paths <- c(
  "revised/iteration1/manuscript/ltee_hypoxia_model.tex",
  "revised/SuppFiles/round3/deep-research-report-7.md",
  "revised/SuppFiles/round3/LTEE_revision_review_package/01_core_results_evidence_index.md",
  "revised/SuppFiles/round3/LTEE_revision_review_package/02_scientific_expansions.md",
  "revised/SuppFiles/round3/LTEE_revision_review_package/03_revised_manuscript_source.zip",
  "revised/SuppFiles/round3/LTEE_revision_review_package/06_revision_diff.md",
  "revised/SuppFiles/round3/LTEE_deliverables/LTEE_core_supporting_results.md",
  "revised/SuppFiles/round3/LTEE_deliverables/LTEE_expansion_notes.md",
  "revised/SuppFiles/round3/LTEE_deliverables/LTEE_revision_diff.md",
  "revised/iteration1/data/Figures/Figure3/source_seed10/invitro_observed_kary.tsv",
  "revised/iteration1/data/Figures/Figure3/source_seed10/invitro_distribution_summary.tsv",
  "revised/iteration1/data/Figures/Figure3/source_seed10/invitro_growth_loglik.tsv",
  "revised/iteration1/data/Figures/Figure3/source_seed10/invitro_flow_overlay.tsv",
  "revised/iteration1/data/Figures/Figure3/source_seed10/best_params.tsv",
  "revised/iteration1/data/Figures/Figure5/figure5_frozen_inputs/selected_results.tsv",
  "revised/iteration1/data/Figures/Figure4/fixed_o2_analysis/data/validation_summary.tsv"
)

abs_paths <- file.path(repo_root, paths)
present <- file.exists(abs_paths)
hashes <- rep(NA_character_, length(paths))
sizes <- rep(NA_real_, length(paths))
modified <- rep(NA_character_, length(paths))
if (any(present)) {
  hashes[present] <- unname(tools::md5sum(abs_paths[present]))
  info <- file.info(abs_paths[present])
  sizes[present] <- info$size
  modified[present] <- format(info$mtime, "%Y-%m-%dT%H:%M:%S%z")
}

manifest <- data.frame(
  relative_path = paths,
  exists = present,
  bytes = sizes,
  md5 = hashes,
  modified_local = modified,
  stringsAsFactors = FALSE
)
write.table(
  manifest,
  file.path(out_root, "provenance", "input_manifest.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(out_root, "provenance", "R_sessionInfo.txt"))
message("Wrote manifest for ", nrow(manifest), " inputs; present=", sum(present))
