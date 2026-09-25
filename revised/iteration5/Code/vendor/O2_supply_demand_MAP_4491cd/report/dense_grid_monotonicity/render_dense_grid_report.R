#!/usr/bin/env Rscript

.dense_report_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "render_dense_grid_report.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.dense_report_root <- normalizePath(file.path(.dense_report_dir, "..", ".."), mustWork = TRUE)
source(file.path(.dense_report_root, "util", "o2_supply_demand_map_dense_grid_utils.R"), local = environment(), chdir = TRUE)
source(file.path(.dense_report_root, "report", "parameter_landscape", "parameter_landscape_report_utils.R"), local = environment(), chdir = TRUE)

dense_report_main <- function(argv = dense_parse_args()) {
  part <- dense_grid_normalize_part(argv$part %||% argv$workflow_part)
  out_dir <- normalizePath(path.expand(argv$out_dir %||% argv$output_root %||% dense_grid_default_out_dir(part)), mustWork = FALSE)
  output_html <- normalizePath(path.expand(argv$output_html %||% file.path(dense_grid_report_dir(out_dir), "index.html")), mustWork = FALSE)
  analysis_manifest_path <- file.path(out_dir, "analysis_artifact_manifest.tsv")
  figure_manifest_path <- file.path(dense_grid_figure_dir(out_dir), "figure_manifest.tsv")
  analysis_manifest <- if (file.exists(analysis_manifest_path)) dense_read_tsv(analysis_manifest_path) else data.frame()
  figure_manifest <- if (file.exists(figure_manifest_path)) dense_read_tsv(figure_manifest_path) else data.frame()
  tables <- list.files(dense_grid_analysis_dir(out_dir), pattern = "[.]tsv$", full.names = TRUE)
  figures <- list.files(dense_grid_figure_dir(out_dir), pattern = "[.]png$", full.names = TRUE, recursive = TRUE)
  sections <- list(
    list(title = "Analysis artifacts", body = o2pl_report_table(analysis_manifest, 200L)),
    list(title = "Materialized analysis tables", body = o2pl_report_table(o2pl_report_file_table(tables, output_html), 200L)),
    list(title = "Visualization artifacts", body = o2pl_report_table(figure_manifest, 200L)),
    list(title = "Rendered figures", body = o2pl_report_image_grid(figures, output_html))
  )
  title <- if (part == "monotonicity") "Dense fixed-O2 ploidy monotonicity report" else "Dense fixed-O2 initial-ploidy trajectory report"
  o2pl_write_report_html(title, sections, output_html)
  dense_record_artifact(out_dir, "report", paste0("dense_grid_", part, "_report"), output_html, producer = "render_dense_grid_report.R", source = analysis_manifest_path)
  message("Wrote dense-grid report: ", output_html)
}

if (sys.nframe() == 0L) dense_report_main()
