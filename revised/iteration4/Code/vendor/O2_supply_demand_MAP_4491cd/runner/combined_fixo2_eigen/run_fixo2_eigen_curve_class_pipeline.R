#!/usr/bin/env Rscript

# Sequential orchestrator: slope analysis -> embedding-table analysis ->
# consume-only visualization -> consume-only HTML report.

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})
WORKFLOW_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
UTIL_DIR <- file.path(WORKFLOW_DIR, "util")
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_path_utils.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_cli_utils.R"))

run_r_stage <- function(label, script, args = character(), dry_run = FALSE) {
  cmd <- c(normalizePath(script, mustWork = TRUE), args)
  message("[", label, "] Rscript ", paste(shQuote(cmd), collapse = " "))
  if (isTRUE(dry_run)) return(invisible(0L))
  status <- system2("Rscript", args = vapply(cmd, shQuote, character(1)))
  if (!identical(as.integer(status), 0L)) {
    stop(label, " failed with exit status ", status, call. = FALSE)
  }
  invisible(status)
}

forward_if_present <- function(argv, names) {
  out <- character()
  for (name in names) {
    value <- argv[[name]]
    if (!is.null(value) && length(value) && !is.na(value[[1L]]) && nzchar(as.character(value[[1L]]))) {
      out <- c(out, paste0("--", name, "=", as.character(value[[1L]])))
    }
  }
  out
}

run_fixo2_eigen_curve_class_pipeline <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  argv <- bpf_parse_args(raw_args)
  repo_root <- bpf_repo_root(SCRIPT_DIR)
  out_dir <- bpf_resolve_repo_path(
    argv$out_dir %||% file.path(
      bpf_combine_fixo2_eigen_attractor_result_dir(repo_root),
      "pooled_embedding_curve_class"
    ),
    repo_root
  )
  classification_dir <- bpf_resolve_repo_path(
    argv$classification_dir %||% file.path(
      bpf_dense_grid_result_root(repo_root),
      "dense-grid_monotonicity_classification"
    ),
    repo_root
  )
  output_html <- bpf_resolve_repo_path(
    argv$output_html %||% file.path(
      out_dir,
      "fixo2_eigen_attractor_embedding_curve_class_report.html"
    ),
    repo_root
  )
  dry_run <- bpf_as_bool(argv$dry_run, FALSE)
  run_report <- bpf_as_bool(argv$run_report, TRUE)
  scripts <- list(
    slope = file.path(WORKFLOW_DIR, "analysis", "combined_fixo2_eigen", "calculate_regression_curve_average_slope.R"),
    prepare = file.path(WORKFLOW_DIR, "analysis", "combined_fixo2_eigen", "prepare_fixo2_eigen_curve_class_tables.R"),
    vis = file.path(WORKFLOW_DIR, "vis", "combined_fixo2_eigen", "plot_fixo2_eigen_attractor_embedding_curve_class.R"),
    report = file.path(WORKFLOW_DIR, "report", "combined_fixo2_eigen", "render_fixo2_eigen_attractor_embedding_curve_class_report.R")
  )
  common <- c(paste0("--out_dir=", out_dir), paste0("--dry_run=", dry_run))
  slope_args <- c(common, forward_if_present(argv, c("curve_table", "by_seed_table", "output_file")))
  prepare_args <- c(
    common,
    forward_if_present(
      argv,
      c("pooled_root", "dense_grid_dir", "class_table", "class_col", "average_slope_table", "reductions", "variants")
    )
  )
  vis_args <- c(common, forward_if_present(argv, c("analysis_manifest", "reductions", "variants")))
  report_args <- c(
    paste0("--classification_dir=", classification_dir),
    paste0("--embedding_dir=", out_dir),
    paste0("--output_html=", output_html),
    paste0("--dry_run=", dry_run)
  )
  run_r_stage("slope analysis", scripts$slope, slope_args, dry_run)
  run_r_stage("embedding table analysis", scripts$prepare, prepare_args, dry_run)
  run_r_stage("visualization", scripts$vis, vis_args, dry_run)
  if (isTRUE(run_report)) {
    run_r_stage("report", scripts$report, report_args, dry_run)
  }
  invisible(output_html)
}

main <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  run_fixo2_eigen_curve_class_pipeline(raw_args)
}

if (identical(environment(), globalenv())) {
  main()
}
