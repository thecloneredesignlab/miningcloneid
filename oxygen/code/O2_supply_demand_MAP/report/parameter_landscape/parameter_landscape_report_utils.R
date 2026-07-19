#!/usr/bin/env Rscript

# Presentation-only HTML helpers.  This file never derives scientific values;
# it formats already materialized tables and references already rendered plots.
.o2pl_report_utils_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[
    basename(frame_files) == "parameter_landscape_report_utils.R"
  ]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(
  file.path(
    .o2pl_report_utils_dir,
    "..",
    "..",
    "util",
    "o2_supply_demand_map_html_utils.R"
  ),
  local = TRUE,
  chdir = TRUE
)
rm(.o2pl_report_utils_dir)

o2pl_report_escape <- o2sd_html_escape_standard

o2pl_report_table <- function(data, max_rows = 100L) {
  if (!is.data.frame(data) || !nrow(data)) return("<p class='empty'>No materialized rows were available.</p>")
  data <- head(data, as.integer(max_rows))
  header <- paste0("<th>", o2pl_report_escape(names(data)), "</th>", collapse = "")
  rows <- apply(data, 1L, function(row) paste0("<tr>", paste0("<td>", o2pl_report_escape(row), "</td>", collapse = ""), "</tr>"))
  paste0("<div class='table-wrap'><table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>")
}

o2pl_report_file_table <- function(paths, output_html) {
  paths <- unique(normalizePath(paths[file.exists(paths)], mustWork = FALSE))
  if (!length(paths)) return(data.frame())
  data.frame(
    file = basename(paths),
    path = vapply(paths, function(path) {
      relative <- tryCatch(file.path(".", substring(path, nchar(dirname(output_html)) + 2L)), error = function(e) path)
      if (startsWith(path, paste0(normalizePath(dirname(output_html), mustWork = FALSE), "/"))) relative else path
    }, character(1)),
    bytes = as.numeric(file.info(paths)$size),
    stringsAsFactors = FALSE
  )
}

o2pl_report_image_grid <- function(paths, output_html) {
  paths <- unique(normalizePath(paths[file.exists(paths)], mustWork = FALSE))
  if (!length(paths)) return("<p class='empty'>No rendered figures were available.</p>")
  cards <- vapply(paths, function(path) {
    base <- normalizePath(dirname(output_html), mustWork = FALSE)
    src <- if (startsWith(path, paste0(base, "/"))) substring(path, nchar(base) + 2L) else path
    paste0("<figure><img src='", o2pl_report_escape(src), "' alt='", o2pl_report_escape(basename(path)), "'><figcaption>", o2pl_report_escape(basename(path)), "</figcaption></figure>")
  }, character(1))
  paste0("<div class='figure-grid'>", paste(cards, collapse = ""), "</div>")
}

o2pl_write_report_html <- function(title, sections, output_html) {
  dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
  section_html <- paste0(vapply(sections, function(section) {
    paste0("<section><h2>", o2pl_report_escape(section$title), "</h2>", section$body, "</section>")
  }, character(1)), collapse = "")
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>", o2pl_report_escape(title), "</title><style>",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:1200px;margin:0 auto;padding:28px;color:#1f2937}",
    "h1,h2{color:#111827}section{margin:2rem 0}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.86rem}",
    "th,td{border:1px solid #d1d5db;padding:.42rem;text-align:left;white-space:nowrap}th{background:#f3f4f6}.figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:1rem}",
    "figure{margin:0;border:1px solid #d1d5db;padding:.7rem}img{width:100%;height:auto}figcaption{font-size:.78rem;color:#4b5563}.empty{color:#6b7280;font-style:italic}",
    "</style></head><body><h1>", o2pl_report_escape(title), "</h1>", section_html, "</body></html>"
  )
  writeLines(html, output_html, useBytes = TRUE)
  invisible(output_html)
}
