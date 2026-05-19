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
REPORT_SCRIPT_DIR <- normalizePath(.o2sd_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(REPORT_SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)
source(file.path(WORKFLOW_ROOT, "analysis", "plot_joint_qc_diagnostics.R"), local = environment())

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

report_magick_available <- function() {
  if (identical(Sys.getenv("O2G_REPORT_FORCE_NO_MAGICK", unset = ""), "TRUE")) {
    return(FALSE)
  }
  requireNamespace("magick", quietly = TRUE)
}

report_gs_available <- function() {
  nzchar(Sys.which("gs"))
}

file_to_data_uri <- function(path, mime) {
  if (requireNamespace("base64enc", quietly = TRUE)) {
    return(base64enc::dataURI(file = path, mime = mime))
  }
  base64_bin <- Sys.which("base64")
  if (nzchar(base64_bin)) {
    attempts <- list(c("-w", "0", path), path, c("-i", path))
    for (args in attempts) {
      out <- tryCatch(
        suppressWarnings(system2(base64_bin, args = args, stdout = TRUE, stderr = TRUE)),
        error = function(e) character()
      )
      status <- attr(out, "status")
      if ((is.null(status) || identical(status, 0L)) && length(out) > 0L) {
        enc <- gsub("[[:space:]]+", "", paste(out, collapse = ""))
        if (nzchar(enc)) {
          return(sprintf("data:%s;base64,%s", mime, enc))
        }
      }
    }
  }
  stop("HTML report requires either the R package 'base64enc' or a system 'base64' command.")
}

render_pdf_preview_png_gs <- function(src_pdf, dest_png, density = 180) {
  gs_bin <- Sys.which("gs")
  if (!nzchar(gs_bin)) {
    stop("Ghostscript ('gs') was requested for PDF preview rendering but is not available in PATH.")
  }
  args <- c(
    "-dSAFER",
    "-dBATCH",
    "-dNOPAUSE",
    "-sDEVICE=pngalpha",
    sprintf("-r%d", as.integer(density)),
    sprintf("-sOutputFile=%s", shQuote(normalizePath(dest_png, mustWork = FALSE))),
    shQuote(normalizePath(src_pdf, mustWork = TRUE))
  )
  out <- suppressWarnings(system2(gs_bin, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("Ghostscript failed while rendering PDF preview for ", src_pdf, ": ", paste(out, collapse = "\n"))
  }
  if (!file.exists(dest_png)) {
    stop("Ghostscript did not create expected PNG preview: ", dest_png)
  }
  normalizePath(dest_png, mustWork = TRUE)
}

pdf_to_data_uri <- function(pdf_path, density = 180) {
  if (!report_magick_available() && !report_gs_available()) {
    return(file_to_data_uri(pdf_path, mime = "application/pdf"))
  }
  png_path <- tempfile("joint_qc_diagnostics_", fileext = ".png")
  on.exit(unlink(png_path, force = TRUE), add = TRUE)
  if (report_magick_available()) {
    img <- magick::image_read(pdf_path, density = density)
    if (length(img) > 1L) img <- img[1]
    magick::image_write(img, path = png_path, format = "png")
  } else {
    render_pdf_preview_png_gs(pdf_path, png_path, density = density)
  }
  file_to_data_uri(png_path, mime = "image/png")
}

figure_media_html <- function(pdf_path, title) {
  data_uri <- pdf_to_data_uri(pdf_path)
  if (report_magick_available() || report_gs_available()) {
    return(sprintf(
      '<div class="figure-media"><img src="%s" alt="%s"/></div>',
      data_uri,
      escape_html(title)
    ))
  }
  sprintf(
    paste0(
      '<div class="figure-media">',
      '<object data="%s" type="application/pdf" class="figure-object">',
      '<div class="figure-fallback"><a href="%s">Open PDF figure</a></div>',
      '</object></div>'
    ),
    data_uri,
    escape_html(basename(pdf_path))
  )
}

build_figure_section <- function(index, title, legend, pdf_path) {
  paste0(
    '<section class="figure-section" id="figure-', index, '">',
    figure_media_html(pdf_path, title),
    '<h2>Figure ', index, '. ', escape_html(title), '</h2>',
    '<p class="legend">', escape_html(legend), '</p>',
    '<p class="file"><code>', escape_html(basename(pdf_path)), '</code></p>',
    '</section>'
  )
}

build_joint_qc_report_html <- function(run_label, seed_summary_path, figure_paths) {
  figure_specs <- list(
    list(
      title = "QC vs Joint Objective Components",
      legend = "Seed-level qc against total, in vivo, and in vitro objectives. Each panel reports a Spearman correlation; lower objective is better.",
      path = figure_paths[["objective"]]
    ),
    list(
      title = "QC vs Key Joint Parameters",
      legend = "Seed-level qc against p_misseg, mu_hp, gamma_growth, and buffer_n_exp fitted values. Each panel reports a Spearman correlation.",
      path = figure_paths[["parameter"]]
    )
  )
  nav <- vapply(seq_along(figure_specs), function(i) {
    sprintf(
      '<li><a href="#figure-%d">Figure %d %s</a></li>',
      i,
      i,
      escape_html(figure_specs[[i]]$title)
    )
  }, character(1))
  figures <- vapply(seq_along(figure_specs), function(i) {
    build_figure_section(
      index = i,
      title = figure_specs[[i]]$title,
      legend = figure_specs[[i]]$legend,
      pdf_path = figure_specs[[i]]$path
    )
  }, character(1))

  paste0(
    '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    '<title>Joint QC Diagnostics</title>',
    '<style>',
    'body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f7fa;color:#1b2a38;}',
    '.shell{max-width:1280px;margin:0 auto;padding:24px;}',
    '.hero,.figure-section{border:1px solid #d6dde6;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,0.05);}',
    '.hero{padding:18px 20px;margin-bottom:18px;}',
    '.hero h1{margin:0 0 6px 0;font-size:28px;line-height:1.15;}',
    '.meta{margin:0;color:#516274;font-size:13px;line-height:1.45;}',
    '.nav{margin:14px 0 0 0;padding:0;list-style:none;display:flex;gap:10px;flex-wrap:wrap;}',
    '.nav a{display:block;padding:7px 10px;border-radius:7px;background:#eef3f8;color:#17324c;text-decoration:none;font-size:13px;font-weight:600;}',
    '.figure-section{padding:14px;margin-bottom:18px;}',
    '.figure-media{margin:0 0 9px 0;}',
    '.figure-media img{display:block;width:100%;max-width:100%;height:auto;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.figure-object{display:block;width:100%;min-height:720px;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.figure-fallback{padding:18px;text-align:center;}',
    '.figure-fallback a{color:#2f6ea4;font-weight:600;text-decoration:none;}',
    'h2{margin:0 0 5px 0;font-size:16px;line-height:1.22;}',
    '.legend{margin:0 0 5px 0;font-size:12px;line-height:1.35;color:#425365;}',
    '.file{margin:0;color:#5f7082;font-size:11px;}',
    'code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}',
    '</style></head><body><main class="shell">',
    '<section class="hero">',
    '<h1>Joint QC Diagnostics</h1>',
    '<p class="meta"><strong>Run:</strong> ', escape_html(run_label), '<br/>',
    '<strong>Seed summary:</strong> <code>', escape_html(seed_summary_path), '</code><br/>',
    '<strong>Generated at:</strong> ', escape_html(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), '</p>',
    '<ul class="nav">', paste(nav, collapse = ""), '</ul>',
    '</section>',
    paste(figures, collapse = ""),
    '</main></body></html>'
  )
}

default_out_dir <- function(seed_summary_path = NULL, extra_results_dir = NULL) {
  if (!is.null(extra_results_dir) && nzchar(trimws(extra_results_dir))) {
    return(file.path(extra_results_dir, "joint_qc_diagnostics"))
  }
  file.path(dirname(seed_summary_path), "joint_qc_diagnostics")
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  seed_summary_path <- resolve_seed_summary_path(
    seed_summary_path = argv$seed_summary,
    extra_results_dir = argv$extra_results_dir
  )
  out_dir <- argv$out_dir %||% default_out_dir(
    seed_summary_path = seed_summary_path,
    extra_results_dir = argv$extra_results_dir
  )
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_label <- argv$run_label %||% infer_run_label_from_paths(
    seed_summary_path = seed_summary_path,
    extra_results_dir = argv$extra_results_dir
  )
  figure_paths <- plot_joint_qc_diagnostics(
    seed_summary_path = seed_summary_path,
    extra_results_dir = argv$extra_results_dir,
    out_dir = out_dir,
    run_label = run_label
  )
  out_path <- normalizePath(
    argv$out_path %||% file.path(out_dir, "joint_qc_diagnostics_report.html"),
    mustWork = FALSE
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  html <- build_joint_qc_report_html(
    run_label = run_label,
    seed_summary_path = seed_summary_path,
    figure_paths = figure_paths
  )
  writeLines(html, con = out_path, useBytes = TRUE)
  message("Wrote joint QC diagnostics report: ", out_path)
}

if (sys.nframe() == 0) {
  main()
}
