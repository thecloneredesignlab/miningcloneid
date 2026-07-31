#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
source(file.path(SCRIPT_DIR, "dominant_ploidy_parameter_contribution_analysis.R"))

append_arg <- function(args, key, value) {
  if (is.null(value) || length(value) == 0L || all(is.na(value))) return(args)
  value <- as.character(value[[1L]])
  if (!nzchar(value)) return(args)
  c(args, paste0("--", key, "=", value))
}

run_report_script <- function(label, script, args) {
  script_path <- file.path(SCRIPT_DIR, script)
  if (!file.exists(script_path)) stop("Missing report script: ", script_path)
  rscript <- file.path(R.home("bin"), "Rscript")
  cmd <- c("--vanilla", script_path, args)
  message("")
  message("[", format(Sys.time(), "%F %T"), "] ", label)
  message("Command: ", paste(shQuote(c(rscript, cmd)), collapse = " "))
  status <- system2(rscript, cmd)
  if (!identical(status, 0L)) {
    stop("Report script failed with status ", status, ": ", script)
  }
  invisible(status)
}

render_dominant_ploidy_parameter_contribution_report <- function(argv) {
  result_root <- normalizePath(
    path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()),
    mustWork = FALSE
  )
  contribution_dir <- normalizePath(
    path.expand(argv$output_dir %||% mode_contribution_default_output_dir(result_root, "dominant_ploidy")),
    mustWork = FALSE
  )
  output_html <- argv$report_html %||% argv$dominant_ploidy_report_html %||% file.path(contribution_dir, "dominant_ploidy_parameter_contribution_report.html")
  args <- character()
  args <- append_arg(args, "result_root", result_root)
  args <- append_arg(args, "mode_contribution_dir", contribution_dir)
  args <- append_arg(args, "output_html", output_html)
  args <- append_arg(args, "embed_assets", argv$embed_report_assets %||% argv$embed_assets %||% "TRUE")
  args <- append_arg(args, "top_n", argv$report_top_n)
  run_report_script("Render dominant ploidy parameter contribution HTML report", "dominant_ploidy_parameter_contribution_report.R", args)
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
argv$mode_contribution_target <- "dominant_ploidy"
run_dominant_ploidy_parameter_contribution_analysis(argv)
if (as_bool(argv$render_html_report %||% argv$run_report %||% "TRUE", TRUE)) {
  render_dominant_ploidy_parameter_contribution_report(argv)
} else {
  message("Skipping dominant ploidy parameter contribution HTML report because --render_html_report=FALSE")
}
