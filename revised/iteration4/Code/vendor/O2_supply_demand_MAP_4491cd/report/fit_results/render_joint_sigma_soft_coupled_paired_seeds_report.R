#!/usr/bin/env Rscript

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_html_utils.R"), local = TRUE)
parse_args <- o2fr_parse_args
`%||%` <- o2fr_null_coalesce
escape_html <- o2sd_html_escape_minimal

table_html <- function(tab, n = 20L) {
  if (!nrow(tab)) return("<p>No rows.</p>")
  tab <- utils::head(tab, n)
  head <- paste0("<tr>", paste0("<th>", escape_html(names(tab)), "</th>", collapse = ""), "</tr>")
  rows <- apply(tab, 1L, function(row) paste0("<tr>", paste0("<td>", escape_html(row), "</td>", collapse = ""), "</tr>"))
  paste0("<table>", head, paste(rows, collapse = "\n"), "</table>")
}

main <- function(argv = parse_args()) {
  analysis_dir <- normalizePath(argv$analysis_dir %||% argv$out_dir, mustWork = TRUE)
  viz_dir <- normalizePath(argv$viz_dir %||% argv$out_dir %||% analysis_dir, mustWork = TRUE)
  report_dir <- normalizePath(argv$report_dir %||% argv$out_dir %||% analysis_dir, mustWork = FALSE)
  required <- c(file.path(analysis_dir, "analysis_manifest.tsv"), file.path(viz_dir, "visualization_manifest.tsv"))
  if (any(!file.exists(required))) stop("Joint-sigma report requires completed analysis and visualization manifests.", call. = FALSE)
  prefix <- "joint_sigma_soft_coupled_paired_seed_comparison"
  summary_path <- file.path(analysis_dir, paste0(prefix, "_summary.tsv"))
  contrast_path <- file.path(analysis_dir, paste0(prefix, "_two_sigma_contrast.tsv"))
  if (!file.exists(summary_path)) stop("Missing joint-sigma summary: ", summary_path, call. = FALSE)
  summary <- utils::read.delim(summary_path, check.names = FALSE, stringsAsFactors = FALSE)
  contrast <- if (file.exists(contrast_path) && file.info(contrast_path)$size > 1L) utils::read.delim(contrast_path, check.names = FALSE, stringsAsFactors = FALSE) else data.frame()
  figures <- utils::read.delim(file.path(viz_dir, "visualization_manifest.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
  pngs <- as.character(figures$file[grepl("[.]png$", figures$file)])
  image_blocks <- vapply(pngs, function(name) {
    path <- file.path(viz_dir, name)
    if (!file.exists(path)) return("")
    src <- if (requireNamespace("base64enc", quietly = TRUE)) paste0("data:image/png;base64,", base64enc::base64encode(path)) else normalizePath(path, mustWork = TRUE)
    paste0("<figure><img src=\"", src, "\"><figcaption>", escape_html(name), "</figcaption></figure>")
  }, character(1))
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Joint sigma paired-seed comparison</title>",
    "<style>body{font-family:Arial;margin:2rem;max-width:1500px}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #ccc;padding:4px}img{max-width:100%}figure{margin:2rem 0}</style></head><body>",
    "<h1>Joint sigma soft-coupled paired-seed comparison</h1>",
    paste0("<p>Rows in paired summary: ", nrow(summary), "; unique seeds: ", length(unique(summary$seed)), ".</p>"),
    "<h2>Paired summary preview</h2>", table_html(summary),
    "<h2>Two-sigma contrast preview</h2>", table_html(contrast),
    "<h2>Figures</h2>", image_blocks, "</body></html>"
  )
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(report_dir, paste0(prefix, ".html"))
  writeLines(html, out_path)
  o2sd_inject_report_image_lightbox(out_path)
  utils::write.table(data.frame(stage = "report", file = basename(out_path), stringsAsFactors = FALSE), file.path(report_dir, "report_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote joint-sigma report: ", out_path)
  invisible(out_path)
}

if (sys.nframe() == 0L) main()
