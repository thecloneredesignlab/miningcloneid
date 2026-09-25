#!/usr/bin/env Rscript

# Report-only ploidy-regime summary.

.ploidy_report_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_invivo_ploidy_regime_report.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.ploidy_report_root <- normalizePath(file.path(.ploidy_report_dir, "..", ".."), mustWork = FALSE)
source(
  file.path(.ploidy_report_root, "util", "o2_supply_demand_map_process_fingerprint_utils.R"),
  local = environment()
)
rm(.ploidy_report_dir, .ploidy_report_root)
parse_args <- o2ipa_parse_args
`%||%` <- o2ipa_null_coalesce
read_tsv <- o2ipa_read_tsv

run_invivo_ploidy_regime_report <- function(argv = parse_args()) {
  simulation_dir <- argv$simulation_dir
  analysis_dir <- argv$analysis_dir
  viz_dir <- argv$viz_dir
  report_dir <- argv$report_dir %||% argv$out_dir
  required <- c(
    file.path(simulation_dir, "simulation_manifest.tsv"),
    file.path(analysis_dir, "analysis_manifest.tsv"),
    file.path(viz_dir, "visualization_manifest.tsv")
  )
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Ploidy-regime report requires completed upstream manifests: ", paste(missing, collapse = ", "))
  if (is.null(report_dir) || !nzchar(report_dir)) stop("Missing --report_dir (or --out_dir).")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  labels <- read_tsv(file.path(analysis_dir, "tables", "trajectory_regime_labels.tsv"))
  evidence <- read_tsv(file.path(analysis_dir, "tables", "final_evidence_classification.tsv"))
  figures <- read_tsv(file.path(viz_dir, "visualization_manifest.tsv"))
  counts <- table(labels$trajectory_regime, useNA = "ifany")
  lines <- c(
    "# In vivo ploidy-regime report", "",
    paste0("- Analyzed seeds: ", nrow(labels)),
    paste0("- Figures: ", nrow(figures)), "", "## Regime counts", "",
    paste0("- ", names(counts), ": ", as.integer(counts)), "", "## Evidence classification", "",
    paste0("- objective cutoff ", evidence$objective_cutoff, ": ", evidence$evidence_classification)
  )
  path <- file.path(report_dir, "analysis_summary.md")
  writeLines(lines, path)
  manifest <- data.frame(artifact = "analysis_summary", relative_path = "analysis_summary.md", role = "report", path = normalizePath(path, mustWork = TRUE), exists = TRUE)
  utils::write.table(manifest, file.path(report_dir, "report_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(report_dir)
}

if (sys.nframe() == 0L) run_invivo_ploidy_regime_report()
