#!/usr/bin/env Rscript

# Derive factorial endpoint tables from pre-existing simulation trajectories.

suppressPackageStartupMessages(library(dplyr))

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_perturbation_utils.R"), local = environment())
rm(.script_dir)

make_factorial_endpoint_summary <- function(design_out, burden_all, ploidy_summary) {
  endpoint_burden <- burden_all %>%
    group_by(scenario_id, experiment) %>%
    arrange(day, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    rename(
      endpoint_p_mis_base = p_mis_base,
      endpoint_p_wgd = p_wgd,
      endpoint_o2_S0 = o2_S0
    ) %>%
    select(
      scenario_id, experiment, day, segment, local_day, step,
      endpoint_p_mis_base, endpoint_p_wgd, endpoint_o2_S0,
      starts_with("pred_burden"), starts_with("pred_log10"), starts_with("pred_o2")
    )
  endpoint_ploidy <- ploidy_summary %>%
    group_by(scenario_id, experiment) %>%
    arrange(day, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    select(
      scenario_id, experiment, day, live_cells, mean_N, mean_ploidy,
      frac_ploidy_lt2, frac_ploidy_2to3, frac_ploidy_3to4,
      frac_ploidy_ge4, frac_ploidy_ge5
    )
  design_out %>%
    left_join(endpoint_burden, by = c("scenario_id", "experiment")) %>%
    left_join(endpoint_ploidy, by = c("scenario_id", "experiment", "day"))
}

run_factorial_interaction_analysis <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  mode <- normalize_interaction_mode(argv$mode %||% "triggered_drug")
  fit_dir <- resolve_path_value(argv$fit_dir %||% argv$run_dir, getwd())
  if (!is.null(fit_dir)) fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  simulation_dir <- resolve_path_value(argv$simulation_dir %||% argv$input_dir, getwd())
  analysis_dir <- resolve_path_value(argv$analysis_dir %||% argv$out_dir, getwd())
  stem <- interaction_output_stem(mode)
  if (is.null(simulation_dir)) {
    if (is.null(fit_dir)) stop("Provide --simulation_dir or --fit_dir.")
    simulation_dir <- file.path(fit_dir, "simulation", "interactions", stem)
  }
  if (is.null(analysis_dir)) {
    if (is.null(fit_dir)) stop("Provide --analysis_dir when --fit_dir is omitted.")
    analysis_dir <- file.path(fit_dir, "analysis", "interactions", stem)
  }
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  validate_artifact_manifest(
    simulation_dir,
    "simulation_manifest.tsv",
    c("interaction_design.tsv", "burden_timecourse.tsv", "ploidy_summary_timecourse.tsv"),
    "Factorial interaction analysis"
  )
  design <- read_required_tsv(file.path(simulation_dir, "interaction_design.tsv"))
  burden <- read_required_tsv(file.path(simulation_dir, "burden_timecourse.tsv"))
  ploidy <- read_required_tsv(file.path(simulation_dir, "ploidy_summary_timecourse.tsv"))
  endpoint <- make_factorial_endpoint_summary(design, burden, ploidy)

  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(endpoint, file.path(analysis_dir, "endpoint_summary.tsv"))
  write_artifact_manifest(
    analysis_dir,
    list(artifact_manifest_row("endpoint_summary", "endpoint_summary.tsv", "analysis_table", endpoint)),
    "analysis_manifest.tsv"
  )
  message("Done. Wrote analysis outputs to: ", normalizePath(analysis_dir, mustWork = FALSE))
  invisible(analysis_dir)
}

if (sys.nframe() == 0L) {
  run_factorial_interaction_analysis()
}
