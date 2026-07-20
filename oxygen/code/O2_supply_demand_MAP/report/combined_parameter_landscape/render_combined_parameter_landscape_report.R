#!/usr/bin/env Rscript

.combined_report_script_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "render_combined_parameter_landscape_report.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.combined_report_root <- normalizePath(file.path(.combined_report_script_dir, "..", ".."), mustWork = TRUE)
source(file.path(.combined_report_root, "util", "o2_supply_demand_map_combined_landscape_utils.R"), local = environment(), chdir = TRUE)
source(file.path(.combined_report_root, "util", "o2_supply_demand_map_html_utils.R"), local = environment(), chdir = TRUE)

combined_html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE); x <- gsub(">", "&gt;", x, fixed = TRUE); gsub('"', "&quot;", x, fixed = TRUE)
}
combined_html_table <- function(data, limit = 12L) {
  if (!nrow(data)) return("<p>No rows.</p>")
  data <- utils::head(data, limit)
  header <- paste0("<tr>", paste0("<th>", combined_html_escape(names(data)), "</th>", collapse = ""), "</tr>")
  body <- apply(data, 1L, function(row) paste0("<tr>", paste0("<td>", combined_html_escape(row), "</td>", collapse = ""), "</tr>"))
  paste0("<table>", header, paste(body, collapse = ""), "</table>")
}
combined_report_relative <- function(path, report_dir) {
  path <- normalizePath(path, mustWork = FALSE); report_dir <- normalizePath(report_dir, mustWork = FALSE)
  prefix <- paste0(report_dir, .Platform$file.sep)
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}

render_combined_parameter_landscape_report <- function(argv = combined_parse_args()) {
  out_dir <- normalizePath(path.expand(argv$out_dir %||% combined_default_out_dir()), mustWork = TRUE)
  report_dir <- combined_report_dir(out_dir)
  output_html <- normalizePath(path.expand(argv$output_html %||% file.path(report_dir, "index.html")), mustWork = FALSE)
  figure_manifest <- file.path(combined_figure_dir(out_dir), "combined_embedding_figure_manifest.tsv")
  analysis_manifest <- file.path(combined_analysis_dir(out_dir), "combined_embedding_table_manifest.tsv")
  slope_table <- file.path(combined_analysis_dir(out_dir), "fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv")
  if (!file.exists(analysis_manifest)) stop("Missing combined analysis manifest: ", analysis_manifest, call. = FALSE)
  analysis <- combined_read_tsv(analysis_manifest)
  figures <- if (file.exists(figure_manifest)) combined_read_tsv(figure_manifest) else data.frame()
  slopes <- if (file.exists(slope_table)) combined_read_tsv(slope_table) else data.frame()
  if (combined_as_bool(argv$dry_run, FALSE)) return(invisible(output_html))
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  cards <- character()
  if (nrow(figures) && "curve_png" %in% names(figures)) for (i in seq_len(nrow(figures))) {
    src <- combined_report_relative(figures$curve_png[[i]], dirname(output_html))
    title <- paste(figures$reduction[[i]], figures$variant[[i]], figures$stub[[i]])
    cards <- c(cards, paste0("<section><h3>", combined_html_escape(title), "</h3><img src=\"", combined_html_escape(src), "\" alt=\"", combined_html_escape(title), "\"></section>"))
  }
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Combined parameter landscape</title><style>",
    "body{font-family:Arial,sans-serif;margin:2rem;color:#222}h1,h2{color:#17365d}section{margin:1.5rem 0}img{max-width:100%;height:auto;border:1px solid #ddd}table{border-collapse:collapse;font-size:.82rem;display:block;overflow:auto}th,td{border:1px solid #ccc;padding:.3rem .45rem;white-space:nowrap}th{background:#eef3f8}",
    "</style></head><body><h1>Combined parameter-landscape and O2-ploidy analysis</h1>",
    "<p>This report assembles only materialized analysis tables and visualization files.</p>",
    "<h2>Embedding inventory</h2>", combined_html_table(analysis),
    "<h2>Average-slope summary</h2>", combined_html_table(slopes),
    "<h2>Figures</h2>", paste(cards, collapse = "\n"), "</body></html>"
  )
  writeLines(html, output_html, useBytes = TRUE)
  o2sd_inject_report_image_lightbox(output_html)
  legacy <- normalizePath(path.expand(argv$legacy_output_html %||% file.path(out_dir, "pooled_embedding_curve_class_report.html")), mustWork = FALSE)
  if (!identical(legacy, output_html)) file.copy(output_html, legacy, overwrite = TRUE)
  message("Wrote combined parameter-landscape report: ", output_html)
  invisible(output_html)
}

main <- render_combined_parameter_landscape_report
if (sys.nframe() == 0L) main()
