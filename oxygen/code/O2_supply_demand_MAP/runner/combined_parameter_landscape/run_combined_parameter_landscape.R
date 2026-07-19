#!/usr/bin/env Rscript

.combined_runner_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
.combined_runner_root <- normalizePath(file.path(.combined_runner_dir, "..", ".."), mustWork = TRUE)
source(file.path(.combined_runner_root, "util", "o2_supply_demand_map_combined_landscape_utils.R"), local = environment(), chdir = TRUE)

combined_runner_child <- function(label, script, raw, dry_run) {
  command <- c("--vanilla", normalizePath(script, mustWork = TRUE), raw)
  message("[", label, "] ", paste(shQuote(c(file.path(R.home("bin"), "Rscript"), command)), collapse = " "))
  if (dry_run) return(invisible(0L))
  status <- system2(file.path(R.home("bin"), "Rscript"), command)
  if (!identical(status, 0L)) stop(label, " failed with status ", status, call. = FALSE)
  invisible(status)
}

run_combined_parameter_landscape <- function(argv = combined_parse_args(), raw = commandArgs(trailingOnly = TRUE)) {
  mode <- tolower(gsub("-", "_", argv$mode %||% argv$action %||% "run", fixed = TRUE))
  dry_run <- combined_as_bool(argv$dry_run, FALSE)
  analysis <- file.path(.combined_runner_root, "analysis", "combined_parameter_landscape", "prepare_combined_parameter_landscape_tables.R")
  visualization <- file.path(.combined_runner_root, "vis", "combined_parameter_landscape", "render_combined_parameter_landscape_figures.R")
  report <- file.path(.combined_runner_root, "report", "combined_parameter_landscape", "render_combined_parameter_landscape_report.R")
  if (mode %in% c("run", "all", "analyze", "analysis")) combined_runner_child("combined analysis", analysis, raw, dry_run)
  if (mode %in% c("run", "all", "visualize", "visualise") && combined_as_bool(argv$generate_figures, TRUE)) combined_runner_child("combined visualization", visualization, raw, dry_run)
  if (mode %in% c("run", "all", "report") && combined_as_bool(argv$run_report, TRUE)) combined_runner_child("combined report", report, raw, dry_run)
  if (!mode %in% c("run", "all", "analyze", "analysis", "visualize", "visualise", "report")) stop("Unknown combined mode: ", mode, call. = FALSE)
  invisible(TRUE)
}

main <- run_combined_parameter_landscape
if (sys.nframe() == 0L) main()
