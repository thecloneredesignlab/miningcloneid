#!/usr/bin/env Rscript

.event_report_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_o2_ploidy_event_coupling_report.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
})
.event_report_root <- normalizePath(file.path(.event_report_dir, "..", ".."), mustWork = FALSE)
source(
  file.path(.event_report_root, "util", "o2_supply_demand_map_report_utils.R"),
  local = environment()
)
rm(.event_report_dir, .event_report_root)

parse_args <- o2sd_report_parse_equals_args
`%||%` <- o2sd_report_null_coalesce_simple
safe_read <- o2sd_report_safe_read_delim

run_o2_ploidy_event_report <- function(argv = parse_args()) {
  simulation_dir <- argv$simulation_dir
  analysis_dir <- argv$analysis_dir
  viz_dir <- argv$viz_dir
  report_dir <- argv$report_dir %||% argv$out_dir
  required <- c(
    file.path(simulation_dir, "simulation_manifest.tsv"),
    file.path(analysis_dir, "analysis_manifest.tsv"),
    file.path(viz_dir, "visualization_manifest.tsv")
  )
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Event report requires completed upstream manifests: ", paste(missing, collapse = ", "))
  if (is.null(report_dir) || !nzchar(report_dir)) stop("Missing --report_dir (or --out_dir).")
  summary <- safe_read(file.path(analysis_dir, "tables", "seed_event_summary.tsv"))
  coupling <- safe_read(file.path(analysis_dir, "tables", "event_coupling_diagnostics.tsv"))
  objective <- safe_read(file.path(analysis_dir, "tables", "objective_and_seed_level_regime_tests.tsv"))
  parameters <- safe_read(file.path(analysis_dir, "tables", "parameter_regime_association.tsv"))
  oxygen <- safe_read(file.path(analysis_dir, "tables", "o2_feature_regime_association.tsv"))
  figures <- safe_read(file.path(viz_dir, "visualization_manifest.tsv"))
  counts <- table(summary$trajectory_regime, useNA = "ifany")
  finite_top <- function(x) {
    if (!nrow(x) || !"BH_adjusted_p_value" %in% names(x)) return(x)
    x <- x[is.finite(x$BH_adjusted_p_value), , drop = FALSE]
    utils::head(x[order(x$BH_adjusted_p_value), , drop = FALSE], 15L)
  }
  lines <- c(
    "# In vivo O2-ploidy event coupling analysis", "",
    paste0("- analyzed seeds: ", nrow(summary)),
    paste0("- figures: ", nrow(figures)), "", "## Regime counts", "",
    paste0("- ", names(counts), ": ", as.integer(counts)), "",
    "## Interpretation boundary", "",
    "This analysis uses materialized fitted trajectories and O2 timecourses. It evaluates mechanistic consistency, not biological causality.",
    "", "## Event-coupling summary", "",
    utils::capture.output(print(coupling, row.names = FALSE)),
    "", "## Objective comparison", "",
    utils::capture.output(print(objective, row.names = FALSE)),
    "", "## Top parameter associations", "",
    utils::capture.output(print(finite_top(parameters), row.names = FALSE)),
    "", "## Top O2 associations", "",
    utils::capture.output(print(finite_top(oxygen), row.names = FALSE))
  )
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(report_dir, "analysis_summary.md")
  writeLines(lines, path)
  manifest <- data.frame(artifact = "analysis_summary", relative_path = "analysis_summary.md", role = "report", path = normalizePath(path, mustWork = TRUE), exists = TRUE)
  utils::write.table(manifest, file.path(report_dir, "report_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(report_dir)
}

if (sys.nframe() == 0L) run_o2_ploidy_event_report()
