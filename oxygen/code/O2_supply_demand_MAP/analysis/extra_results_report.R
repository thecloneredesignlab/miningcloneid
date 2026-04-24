#!/usr/bin/env Rscript

.o2sd_bootstrap_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2sd_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

escape_html <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

make_figure_spec <- function(extra_results_dir, filename, title, legend) {
  path <- file.path(extra_results_dir, filename)
  if (!file.exists(path)) {
    stop("Missing required figure for extra_results report: ", path)
  }
  list(
    filename = filename,
    path = normalizePath(path, mustWork = TRUE),
    title = title,
    legend = legend
  )
}

build_figure_specs <- function(extra_results_dir) {
  list(
    make_figure_spec(
      extra_results_dir,
      "objective_vs_boundary_risk.pdf",
      "Objective vs Boundary Risk",
      "Objective score against minimum relative distance to the nearest fitted parameter bound."
    ),
    make_figure_spec(
      extra_results_dir,
      "objective_components_violin.pdf",
      "Objective Components Violin",
      "Across-seed distributions of total objective, burden objective, and ploidy objective."
    ),
    make_figure_spec(
      extra_results_dir,
      "parameter_boundary_forest.pdf",
      "Parameter Boundary Forest",
      "Relative fitted positions of active parameters within their transformed bounds across seeds."
    ),
    make_figure_spec(
      extra_results_dir,
      "parameter_boundary_forest_pred1000_gt44_top3.pdf",
      "Parameter Boundary Forest (Pred1000 > 44 Top 3)",
      "Boundary forest restricted to the recommended top 3 seeds among runs with both 2N and 4N 1000-day predictions above 44."
    )
  )
}

infer_run_label <- function(extra_results_dir) {
  base <- basename(extra_results_dir)
  if (identical(base, "extra_results")) {
    return(basename(dirname(extra_results_dir)))
  }
  base
}

report_magick_available <- function() {
  if (identical(Sys.getenv("O2SD_REPORT_FORCE_NO_MAGICK", unset = ""), "TRUE")) {
    return(FALSE)
  }
  requireNamespace("magick", quietly = TRUE)
}

report_gs_available <- function() {
  nzchar(Sys.which("gs"))
}

report_base64enc_available <- function() {
  requireNamespace("base64enc", quietly = TRUE)
}

file_to_data_uri <- function(path, mime) {
  if (report_base64enc_available()) {
    return(base64enc::dataURI(file = path, mime = mime))
  }
  base64_bin <- Sys.which("base64")
  if (nzchar(base64_bin)) {
    enc <- tryCatch(
      suppressWarnings(system2(base64_bin, c("-w", "0", path), stdout = TRUE, stderr = TRUE)),
      error = function(e) character()
    )
    if (!length(enc)) {
      enc <- tryCatch(
        suppressWarnings(system2(base64_bin, path, stdout = TRUE, stderr = TRUE)),
        error = function(e) character()
      )
    }
    if (length(enc) > 0L) {
      return(sprintf("data:%s;base64,%s", mime, paste(enc, collapse = "")))
    }
  }
  stop(
    "HTML report fallback requires either the R package 'base64enc' or a system 'base64' command ",
    "when 'magick' is unavailable."
  )
}

render_pdf_preview_png_gs <- function(src_pdf, dest_png, density = 180) {
  gs_bin <- Sys.which("gs")
  if (!nzchar(gs_bin)) {
    stop("Ghostscript ('gs') was requested for PDF preview rendering but is not available in PATH.")
  }
  src_pdf_use <- normalizePath(src_pdf, mustWork = TRUE)
  dest_png_use <- normalizePath(dest_png, mustWork = FALSE)
  density_use <- suppressWarnings(as.integer(density))
  if (!is.finite(density_use) || density_use <= 0L) density_use <- 180L
  args <- c(
    "-dSAFER",
    "-dBATCH",
    "-dNOPAUSE",
    "-sDEVICE=pngalpha",
    sprintf("-r%d", density_use),
    sprintf("-sOutputFile=%s", shQuote(dest_png_use)),
    shQuote(src_pdf_use)
  )
  out <- suppressWarnings(system2(gs_bin, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop(
      "Ghostscript failed while rendering PDF preview for ", src_pdf, ": ",
      paste(out, collapse = "\n")
    )
  }
  if (!file.exists(dest_png_use)) {
    stop("Ghostscript did not create expected PNG preview: ", dest_png_use)
  }
  normalizePath(dest_png_use, mustWork = TRUE)
}

pdf_to_data_uri <- function(pdf_path, density = 180) {
  if (!report_magick_available() && !report_gs_available()) {
    return(file_to_data_uri(pdf_path, mime = "application/pdf"))
  }
  png_path <- tempfile("o2sd_extra_results_", fileext = ".png")
  on.exit(unlink(png_path, force = TRUE), add = TRUE)
  if (report_magick_available()) {
    img <- magick::image_read(pdf_path, density = density)
    if (length(img) > 1L) {
      img <- img[1]
    }
    magick::image_write(img, path = png_path, format = "png")
  } else {
    render_pdf_preview_png_gs(pdf_path, png_path, density = density)
  }
  file_to_data_uri(png_path, mime = "image/png")
}

figure_media_html <- function(fig) {
  data_uri <- pdf_to_data_uri(fig$path)
  if (report_magick_available() || report_gs_available()) {
    sprintf(
      '<div class="report-figure"><img src="%s" alt="%s" class="report-figure-image"/></div>',
      data_uri,
      escape_html(fig$title)
    )
  } else {
    sprintf(
      paste0(
        '<div class="report-figure">',
        '<object data="%s" type="application/pdf" class="report-figure-object">',
        '<div class="report-figure-fallback"><a href="%s">Open PDF figure</a></div>',
        '</object></div>'
      ),
      data_uri,
      escape_html(fig$filename)
    )
  }
}

build_report_html <- function(extra_results_dir, figure_specs) {
  run_label <- infer_run_label(extra_results_dir)
  nav_items <- vapply(seq_along(figure_specs), function(i) {
    fig <- figure_specs[[i]]
    sprintf(
      '<li class="report-nav-item"><a class="report-nav-link" href="#figure-%d">Figure %d %s</a></li>',
      i,
      i,
      escape_html(fig$title)
    )
  }, character(1))

  figure_blocks <- vapply(seq_along(figure_specs), function(i) {
    fig <- figure_specs[[i]]
    sprintf(
      paste0(
        '<section class="report-section" id="figure-%d">',
        '%s',
        '<h2 class="report-figure-title">Figure %d. %s</h2>',
        '<p class="report-figure-legend">%s</p>',
        '<p class="report-figure-file"><code>%s</code></p>',
        '</section>'
      ),
      i,
      figure_media_html(fig),
      i,
      escape_html(fig$title),
      escape_html(fig$legend),
      escape_html(fig$filename)
    )
  }, character(1))

  paste0(
    '<!DOCTYPE html>',
    '<html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    '<title>Extra Results Report</title>',
    '<style>',
    'body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f7fa;color:#1b2a38;}',
    '.report-shell{display:flex;gap:28px;max-width:1600px;margin:0 auto;padding:24px;}',
    '.report-sidebar{position:sticky;top:24px;align-self:flex-start;width:280px;border:1px solid #d6dde6;border-radius:12px;background:#f7f9fb;box-shadow:0 10px 28px rgba(0,0,0,0.08);overflow:hidden;}',
    '.report-sidebar-header{padding:14px;background:linear-gradient(180deg,#1f3348 0%,#284662 100%);color:#fff;}',
    '.report-kicker{font-size:11px;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;opacity:0.78;}',
    '.report-title{margin-top:4px;font-size:18px;font-weight:700;line-height:1.15;}',
    '.report-subtitle{margin-top:4px;font-size:12px;opacity:0.85;}',
    '.report-nav{padding:10px 8px 12px 8px;}',
    '.report-nav-list{margin:0;padding:0;list-style:none;}',
    '.report-nav-item{margin:4px 0;}',
    '.report-nav-link{display:block;padding:10px 12px;border-radius:8px;text-decoration:none;color:#17324c;font-size:14px;font-weight:600;line-height:1.35;}',
    '.report-nav-link:hover{background:rgba(47,110,164,0.08);}',
    '.report-main{flex:1;min-width:0;max-width:1120px;}',
    '.report-hero{margin-bottom:20px;padding:18px 20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}',
    '.report-hero h1{margin:0 0 6px 0;font-size:28px;line-height:1.15;}',
    '.report-meta{margin:0;color:#516274;font-size:14px;}',
    '.report-section{margin-bottom:36px;padding:20px;border:1px solid #d6dde6;border-radius:12px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}',
    '.report-figure{margin:0 0 14px 0;}',
    '.report-figure-image{display:block;width:100%;max-width:100%;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-figure-object{width:100%;min-height:680px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-figure-fallback{padding:18px;text-align:center;}.report-figure-fallback a{color:#2f6ea4;font-weight:600;text-decoration:none;}',
    '.report-figure-title{margin:0 0 8px 0;font-size:22px;line-height:1.2;}',
    '.report-figure-legend{margin:0 0 8px 0;font-size:14px;line-height:1.6;color:#425365;}',
    '.report-figure-file{margin:0;color:#5f7082;font-size:12px;}',
    'code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}',
    '@media (max-width: 991px){.report-shell{display:block;padding:16px;}.report-sidebar{position:relative;top:auto;width:auto;margin-bottom:16px;}.report-main{max-width:none;}}',
    '</style></head><body>',
    '<div class="report-shell">',
    '<aside class="report-sidebar" aria-label="Figure navigation">',
    '<div class="report-sidebar-header">',
    '<div class="report-kicker">Navigation</div>',
    '<div class="report-title">Extra Results</div>',
    '<div class="report-subtitle">Figure guide for ', escape_html(run_label), '</div>',
    '</div>',
    '<nav class="report-nav"><ul class="report-nav-list">', paste(nav_items, collapse = ""), '</ul></nav>',
    '</aside>',
    '<main class="report-main">',
    '<section class="report-hero">',
    '<h1>Extra Results Report</h1>',
    '<p class="report-meta"><strong>Run:</strong> ', escape_html(run_label), '<br/>',
    '<strong>Source directory:</strong> <code>', escape_html(extra_results_dir), '</code><br/>',
    '<strong>Generated at:</strong> ', escape_html(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), '</p>',
    '</section>',
    paste(figure_blocks, collapse = ""),
    '</main></div></body></html>'
  )
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  extra_results_dir <- argv$extra_results_dir
  if (is.null(extra_results_dir) || !nzchar(trimws(extra_results_dir))) {
    stop("Usage: Rscript extra_results_report.R --extra_results_dir=/path/to/extra_results [--out_path=/path/to/extra_results_report.html]")
  }
  extra_results_dir <- normalizePath(extra_results_dir, mustWork = TRUE)
  out_path <- normalizePath(
    argv$out_path %||% file.path(extra_results_dir, "extra_results_report.html"),
    mustWork = FALSE
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  figure_specs <- build_figure_specs(extra_results_dir)
  html <- build_report_html(extra_results_dir = extra_results_dir, figure_specs = figure_specs)
  writeLines(html, con = out_path, useBytes = TRUE)

  message("Wrote extra results report: ", out_path)
}

if (sys.nframe() == 0) {
  main()
}
