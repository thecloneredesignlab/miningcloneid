#!/usr/bin/env Rscript

# Report-only process-fingerprint summary.  All tables and figures must already
# exist; this script performs no simulation, clustering, or plotting.

.process_report_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_invivo_process_fingerprint_report.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.process_report_root <- normalizePath(file.path(.process_report_dir, "..", ".."), mustWork = FALSE)
source(
  file.path(.process_report_root, "util", "o2_supply_demand_map_process_fingerprint_utils.R"),
  local = environment()
)
rm(.process_report_dir, .process_report_root)
parse_args <- o2ipa_parse_args
`%||%` <- o2ipa_null_coalesce
read_tsv <- o2ipa_read_tsv

run_invivo_process_fingerprint_report <- function(argv = parse_args()) {
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
  if (length(missing)) stop("Process report requires completed upstream manifests: ", paste(missing, collapse = ", "))
  if (is.null(report_dir) || !nzchar(report_dir)) stop("Missing --report_dir (or --out_dir).")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  seed <- read_tsv(file.path(analysis_dir, "tables", "seed_manifest.tsv"))
  diag <- read_tsv(file.path(analysis_dir, "tables", "cluster_k_diagnostics.tsv"))
  figures <- read_tsv(file.path(viz_dir, "visualization_manifest.tsv"))
  lines <- c(
    "# In vivo process-fingerprint report",
    "",
    paste0("- Valid seeds: ", sum(seed$fit_success %in% TRUE, na.rm = TRUE)),
    paste0("- Analysis tables: ", nrow(read_tsv(file.path(analysis_dir, "analysis_manifest.tsv")))),
    paste0("- Figures: ", nrow(figures)),
    "",
    "## Cluster diagnostics",
    "",
    paste0("- ", diag$cluster_source, ": k=", diag$k, ", silhouette=", signif(diag$mean_silhouette, 4))
  )
  path <- file.path(report_dir, "analysis_summary.md")
  writeLines(lines, path)
  manifest <- data.frame(artifact = "analysis_summary", relative_path = "analysis_summary.md", role = "report", path = normalizePath(path, mustWork = TRUE), exists = TRUE)
  utils::write.table(manifest, file.path(report_dir, "report_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(report_dir)
}

if (sys.nframe() == 0L) run_invivo_process_fingerprint_report()
