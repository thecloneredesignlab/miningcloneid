#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Figure6 <- function(n_core = 8L) {
  data_dir <- file.path(DATA_ROOT, "Figure6")
  panel_dir <- file.path(data_dir, "panels")
  render_output_dir <- file.path(data_dir, "rendered")
  require_files(file.path(data_dir, c(
    "response_class_smoothed_curves.tsv",
    "response_class_curve_class_by_seed.tsv",
    "pair_surface_best_joint_seed_summary.tsv",
    "pair_surface_best_trajectory.tsv",
    "pair_surface_complete_surface.tsv"
  )), "Figure 6 analysis cache")
  run_oxygen_response_pipeline(
    args = c(
      paste0("--joint_root=", JOINT_RESULT_ROOT),
      paste0("--workspace_dir=", data_dir),
      paste0("--data_dir=", data_dir),
      paste0("--figure_dir=", panel_dir),
      paste0("--output_dir=", render_output_dir),
      paste0("--workers=", as.integer(n_core)),
      "--rebuild=FALSE"
    ),
    env = paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT)
  )
  stage_output_file(
    file.path(render_output_dir, "assembled_oxygen_response.png"),
    "assembled_fig6.png"
  )
  stage_output_file(
    file.path(render_output_dir, "assembled_oxygen_response.pdf"),
    "assembled_fig6.pdf"
  )
  require_files(
    file.path(OUTPUT_ROOT, c("assembled_fig6.png", "assembled_fig6.pdf")),
    "Figure 6 output"
  )
  invisible(file.path(OUTPUT_ROOT, "assembled_fig6.png"))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  n_core_arg <- sub("^--n-core=", "", args[grepl("^--n-core=", args)])
  draw_Figure6(
    n_core = if (length(n_core_arg)) as.integer(n_core_arg[[1L]]) else 8L
  )
}
