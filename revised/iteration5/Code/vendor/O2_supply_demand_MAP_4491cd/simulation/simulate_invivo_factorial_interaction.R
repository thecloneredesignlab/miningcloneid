#!/usr/bin/env Rscript

# Compatibility pipeline wrapper.  New callers should invoke the dedicated
# simulation, analysis, and visualization entrypoints directly.

.wrapper_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.wrapper_dir, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "simulation", "perturbation", "generate_factorial_interaction_outputs.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "analysis", "interactions", "analyze_factorial_interactions.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "vis", "interactions", "plot_factorial_interactions.R"), local = environment())
rm(.wrapper_dir)

run_legacy_factorial_interaction_pipeline <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  warning(
    "Legacy factorial entrypoint: use simulation/perturbation, analysis/interactions, and vis/interactions entrypoints for staged execution.",
    call. = FALSE
  )
  fit_dir <- resolve_path_value(argv$fit_dir %||% argv$run_dir, getwd())
  if (is.null(fit_dir)) stop("Missing required argument: --fit_dir=/path/to/seed_dir")
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  mode <- normalize_interaction_mode(argv$mode %||% "triggered_drug")
  out_dir <- resolve_path_value(argv$out_dir %||% argv$simulation_dir, getwd())
  if (is.null(out_dir)) out_dir <- file.path(fit_dir, "simulation", interaction_output_stem(mode))

  sim_argv <- argv
  sim_argv$fit_dir <- fit_dir
  sim_argv$mode <- mode
  sim_argv$simulation_dir <- out_dir
  run_factorial_interaction_simulation(sim_argv)
  run_factorial_interaction_analysis(list(mode = mode, simulation_dir = out_dir, analysis_dir = out_dir))
  if (as_bool(argv$make_plots, TRUE)) {
    run_factorial_interaction_visualization(list(
      mode = mode, simulation_dir = out_dir, analysis_dir = out_dir, out_dir = out_dir
    ))
  }
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  run_legacy_factorial_interaction_pipeline()
}
