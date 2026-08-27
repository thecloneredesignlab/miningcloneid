#!/usr/bin/env Rscript

.o2pl_report_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "dominant_ploidy_parameter_contribution_report.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.o2pl_root <- normalizePath(file.path(.o2pl_report_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_root, "util", "o2_supply_demand_map_parameter_landscape_io_utils.R"), local = environment(), chdir = TRUE)
source(file.path(.o2pl_report_dir, "parameter_landscape_report_utils.R"), local = environment(), chdir = TRUE)

o2pl_dominant_ploidy_report_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  result_root <- normalizePath(path.expand(argv$result_root %||% o2pl_default_result_root()), mustWork = FALSE)
  contribution_dir <- normalizePath(path.expand(argv$mode_contribution_dir %||% argv$output_dir %||% file.path(result_root, "dominant_ploidy_parameter_contribution")), mustWork = FALSE)
  output_html <- normalizePath(path.expand(argv$output_html %||% argv$report_html %||% file.path(contribution_dir, "dominant_ploidy_parameter_contribution_report.html")), mustWork = FALSE)
  index_path <- file.path(contribution_dir, "dominant_ploidy_parameter_contribution_index.csv")
  features_path <- file.path(contribution_dir, "dominant_ploidy_parameter_top_features_across_reference_o2.csv")
  index <- if (file.exists(index_path)) read_csv_plain(index_path) else data.frame()
  features <- if (file.exists(features_path)) read_csv_plain(features_path) else data.frame()
  figures <- list.files(contribution_dir, pattern = "[.]png$", recursive = TRUE, full.names = TRUE)
  o2pl_write_report_html("Fixed-O2 dominant-ploidy parameter contribution report", list(
    list(title = "Reference-O2 runs", body = o2pl_report_table(index)),
    list(title = "Ranked materialized feature contributions", body = o2pl_report_table(features, 120L)),
    list(title = "Rendered contribution figures", body = o2pl_report_image_grid(figures, output_html))
  ), output_html)
  o2pl_record_artifact(result_root, "report", "dominant_ploidy_parameter_contribution_report", output_html, producer = "dominant_ploidy_parameter_contribution_report.R", source = index_path)
  message("Wrote dominant-ploidy parameter contribution report: ", output_html)
}

if (sys.nframe() == 0L) o2pl_dominant_ploidy_report_main()
