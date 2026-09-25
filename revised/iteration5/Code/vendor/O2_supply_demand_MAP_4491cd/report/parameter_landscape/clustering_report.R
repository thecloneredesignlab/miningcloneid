#!/usr/bin/env Rscript

.o2pl_report_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "clustering_report.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.o2pl_root <- normalizePath(file.path(.o2pl_report_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_root, "util", "o2_supply_demand_map_parameter_landscape_io_utils.R"), local = environment(), chdir = TRUE)
source(file.path(.o2pl_report_dir, "parameter_landscape_report_utils.R"), local = environment(), chdir = TRUE)

o2pl_clustering_report_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  result_root <- normalizePath(path.expand(argv$result_root %||% o2pl_default_result_root()), mustWork = FALSE)
  output_html <- normalizePath(path.expand(argv$output_html %||% file.path(result_root, "parameter_landscape_clustering_report.html")), mustWork = FALSE)
  manifest_path <- o2pl_manifest_path(result_root, "analysis")
  manifest <- if (file.exists(manifest_path)) read_tsv(manifest_path) else data.frame()
  coordinate_paths <- list.files(file.path(result_root, "analysis_tables"), pattern = "embedding_coordinates[.]tsv$", recursive = TRUE, full.names = TRUE)
  coordinate_index <- o2pl_report_file_table(coordinate_paths, output_html)
  figure_paths <- list.files(file.path(result_root, "figures"), pattern = "[.]png$", recursive = TRUE, full.names = TRUE)
  sections <- list(
    list(title = "Analysis artifact manifest", body = o2pl_report_table(manifest)),
    list(title = "Materialized coordinate tables", body = o2pl_report_table(coordinate_index)),
    list(title = "Rendered embeddings", body = o2pl_report_image_grid(figure_paths, output_html))
  )
  o2pl_write_report_html("O2 parameter-landscape clustering report", sections, output_html)
  o2pl_record_artifact(result_root, "report", "parameter_landscape_clustering_report", output_html, producer = "clustering_report.R", source = manifest_path)
  message("Wrote parameter-landscape clustering report: ", output_html)
}

if (sys.nframe() == 0L) o2pl_clustering_report_main()
