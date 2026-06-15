#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos > 0L) out[[substr(kv, 1L, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv)) else out[[kv]] <- TRUE
  }
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) y else x
as_chr <- function(x, default = "") as.character((x %||% default)[[1]])

escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

read_table_optional <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}

table_html <- function(tab, max_rows = 50L) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return('<p class="empty">No data available.</p>')
  tab <- utils::head(tab, max_rows)
  header <- paste0("<tr>", paste(sprintf("<th>%s</th>", escape_html(names(tab))), collapse = ""), "</tr>")
  rows <- apply(tab, 1L, function(row) paste0("<tr>", paste(sprintf("<td>%s</td>", escape_html(row)), collapse = ""), "</tr>"))
  paste0('<table class="report-table">', header, paste(rows, collapse = ""), "</table>")
}

figure_html <- function(root, filename, title, caption) {
  path <- file.path(root, filename)
  if (!file.exists(path)) return("")
  paste0(
    '<section class="card">',
    '<h2>', escape_html(title), '</h2>',
    '<p>', escape_html(caption), '</p>',
    '<object class="figure-object" data="', escape_html(filename), '" type="application/pdf">',
    '<a href="', escape_html(filename), '">Open ', escape_html(filename), '</a>',
    '</object>',
    '</section>'
  )
}

status_counts <- function(summary_df, manifest) {
  total <- if (is.data.frame(manifest)) nrow(manifest) else 0L
  completed <- if (is.data.frame(summary_df)) nrow(summary_df) else 0L
  data.frame(
    metric = c("planned_pairs", "completed_pairs", "missing_or_failed_pairs"),
    value = c(total, completed, max(total - completed, 0L)),
    stringsAsFactors = FALSE
  )
}

decision_text <- function(summary_df) {
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df)) {
    return("No completed warm-up pair summaries were available at report time.")
  }
  if ("final_basin_id" %in% names(summary_df)) {
    basins <- unique(na.omit(summary_df$final_basin_id))
    if (length(basins) <= 1L) {
      return("All completed warm-up pairs assigned to one final parameter basin; this supports warm-start stability for the completed runs.")
    }
    return(paste0("Completed warm-up pairs assigned to ", length(basins), " final parameter basins; compare objective and biological diagnostics before selecting a final result."))
  }
  "Final basin assignments were not available; use objective and diagnostic summaries for the completed runs."
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root <- normalizePath(as_chr(argv$multi_warmup_root, as_chr(argv$out_dir)), mustWork = TRUE)
  out_path <- normalizePath(as_chr(argv$out, file.path(root, "multi-warm-up_results.html")), mustWork = FALSE)
  manifest <- read_table_optional(file.path(root, "multi_warmup_manifest.tsv"))
  invivo_clusters <- read_table_optional(file.path(root, "multi_warmup_invivo_cluster_summary.tsv"))
  summary_df <- read_table_optional(file.path(root, "multi_warmup_best_seed_summary.tsv"))
  basins <- read_table_optional(file.path(root, "multi_warmup_final_basin_assignments.tsv"))
  jobs <- read_table_optional(file.path(root, "multi_warmup_jobs.tsv"))

  html <- paste0(
    '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    '<title>Multi-Warm-Up Results</title>',
    '<style>',
    'body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f7fa;color:#1b2a38;}',
    '.shell{max-width:1500px;margin:0 auto;padding:24px;}',
    '.hero,.card{background:#fff;border:1px solid #d6dde6;border-radius:10px;box-shadow:0 8px 22px rgba(0,0,0,.05);padding:16px;margin-bottom:18px;}',
    'h1{margin:0 0 8px;font-size:28px;}h2{margin:0 0 8px;font-size:18px;}p{color:#425365;line-height:1.45;}',
    '.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px;}',
    '.figure-object{width:100%;min-height:620px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-table{width:100%;border-collapse:collapse;font-size:12px;}.report-table th,.report-table td{padding:7px 9px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}',
    '.report-table th{background:#f7f9fb;font-weight:700;}.empty{font-style:italic;color:#657789;}code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}',
    '@media(max-width:900px){.grid{grid-template-columns:1fr}.shell{padding:14px}.figure-object{min-height:420px;}}',
    '</style></head><body><main class="shell">',
    '<section class="hero"><h1>Multi-Warm-Up Results</h1>',
    '<p><strong>Root:</strong> <code>', escape_html(root), '</code><br/>',
    '<strong>Generated:</strong> ', escape_html(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), '<br/>',
    '<strong>Decision summary:</strong> ', escape_html(decision_text(summary_df)), '</p></section>',
    '<section class="card"><h2>Run Status</h2>', table_html(status_counts(summary_df, manifest), max_rows = 20), '</section>',
    figure_html(root, "joint_soft_coupling_ratio_umap_500seed.pdf", "Joint Soft Coupling Ratio UMAP", "Pre-fitting source top-seed ratio UMAP used to choose multi-warm-up in vivo representatives."),
    '<div class="grid">',
    figure_html(root, "multi_warmup_objective_summary.pdf", "Objective Summary", "Best total objective for each completed warm-up pair."),
    figure_html(root, "multi_warmup_invivo_invitro_objective_scatter.pdf", "In Vivo vs In Vitro Objective", "Best-seed in vivo and in vitro objective components by warm-up pair."),
    figure_html(root, "multi_warmup_final_parameter_ratio_heatmap.pdf", "Final Parameter Ratios", "Best-seed final in vivo/in vitro parameter ratios."),
    figure_html(root, "multi_warmup_final_parameter_basin_pca.pdf", "Final Basin PCA", "PCA of best-seed final parameter-ratio vectors."),
    '</div>',
    '<section class="card"><h2>Warm-Up Manifest</h2>', table_html(manifest, max_rows = 200), '</section>',
    '<section class="card"><h2>In Vivo Cluster Summary</h2>', table_html(invivo_clusters, max_rows = 200), '</section>',
    '<section class="card"><h2>Best Seed Summary</h2>', table_html(summary_df, max_rows = 200), '</section>',
    '<section class="card"><h2>Final Basin Assignments</h2>', table_html(basins, max_rows = 200), '</section>',
    '<section class="card"><h2>HPC/Run Jobs</h2>', table_html(jobs, max_rows = 300), '</section>',
    '</main></body></html>'
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_path, useBytes = TRUE)
  message("Wrote multi-warm-up report: ", out_path)
}

if (identical(environment(), globalenv())) {
  main()
}
