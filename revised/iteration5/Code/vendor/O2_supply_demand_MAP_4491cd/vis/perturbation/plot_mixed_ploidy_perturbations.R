#!/usr/bin/env Rscript

# Pure visualization consumer for mixed-ploidy perturbation outputs.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

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

plot_mixed_outputs <- function(burden_all, ploidy_summary, design, out_dir) {
  burden_plot <- burden_all %>%
    mutate(
      varied_value_label = format(signif(varied_value, 4), trim = TRUE),
      initial_burden_label = paste0("init burden=", format(signif(initial_burden_cells, 3), scientific = TRUE)),
      trigger_label = ifelse(
        is.finite(trigger_burden_cells),
        paste0("trigger burden=", format(signif(trigger_burden_cells, 3), scientific = TRUE)),
        initial_burden_label
      )
    )
  ploidy_plot <- ploidy_summary %>%
    left_join(design, by = c("scenario_id", "experiment")) %>%
    mutate(
      varied_value_label = format(signif(varied_value, 4), trim = TRUE),
      initial_burden_label = paste0("init burden=", format(signif(initial_burden_cells, 3), scientific = TRUE)),
      trigger_label = ifelse(
        is.finite(trigger_burden_cells),
        paste0("trigger burden=", format(signif(trigger_burden_cells, 3), scientific = TRUE)),
        initial_burden_label
      )
    )

  static_burden <- burden_plot %>% filter(experiment %in% c("o2_initial_supply", "p_wgd_sweep"))
  if (nrow(static_burden)) {
    p <- ggplot(static_burden, aes(day, pred_log10_burden_cells, color = varied_value_label, group = scenario_id)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      facet_grid(experiment ~ initial_burden_label, scales = "free_y") +
      labs(x = "Day", y = "log10 Burden (total cells)", color = "Value") +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "burden_static_sweeps.pdf"), p, width = 14, height = 7)
  }
  static_ploidy <- ploidy_plot %>% filter(experiment %in% c("o2_initial_supply", "p_wgd_sweep"))
  if (nrow(static_ploidy)) {
    p <- ggplot(static_ploidy, aes(day, mean_ploidy, color = varied_value_label, group = scenario_id)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      facet_grid(experiment ~ initial_burden_label, scales = "free_y") +
      labs(x = "Day", y = "Mean ploidy", color = "Value") +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "mean_ploidy_static_sweeps.pdf"), p, width = 14, height = 7)
  }
  tx_burden <- burden_plot %>% filter(experiment == "pmiss_triggered_treatment")
  if (nrow(tx_burden)) {
    p <- ggplot(tx_burden, aes(day, pred_log10_burden_cells, color = varied_value_label, group = scenario_id)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      geom_vline(aes(xintercept = trigger_day), linetype = "dashed", linewidth = 0.25, alpha = 0.4) +
      facet_wrap(~ trigger_label, scales = "free_y") +
      labs(x = "Day", y = "log10 Burden (total cells)", color = "Post p_mis_base") +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "burden_triggered_pmiss_treatment.pdf"), p, width = 12, height = 8)
  }
  tx_ploidy <- ploidy_plot %>% filter(experiment == "pmiss_triggered_treatment")
  if (nrow(tx_ploidy)) {
    p <- ggplot(tx_ploidy, aes(day, mean_ploidy, color = varied_value_label, group = scenario_id)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      geom_vline(aes(xintercept = trigger_day), linetype = "dashed", linewidth = 0.25, alpha = 0.4) +
      facet_wrap(~ trigger_label, scales = "free_y") +
      labs(x = "Day", y = "Mean ploidy", color = "Post p_mis_base") +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "mean_ploidy_triggered_pmiss_treatment.pdf"), p, width = 12, height = 8)
  }
}

run_mixed_ploidy_perturbation_visualization <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- resolve_path_value(argv$fit_dir %||% argv$run_dir, getwd())
  if (!is.null(fit_dir)) fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  simulation_dir <- resolve_path_value(argv$simulation_dir %||% argv$input_dir, getwd())
  analysis_dir <- resolve_path_value(argv$analysis_dir, getwd())
  out_dir <- resolve_path_value(argv$viz_dir %||% argv$out_dir, getwd())
  if (is.null(simulation_dir)) {
    if (is.null(fit_dir)) stop("Provide --simulation_dir or --fit_dir.")
    simulation_dir <- file.path(fit_dir, "simulation", "perturbation", "mixed_ploidy")
  }
  if (is.null(analysis_dir)) {
    if (is.null(fit_dir)) stop("Provide --analysis_dir when --fit_dir is omitted.")
    analysis_dir <- file.path(fit_dir, "analysis", "perturbation", "mixed_ploidy")
  }
  if (is.null(out_dir)) {
    if (is.null(fit_dir)) stop("Provide --out_dir when --fit_dir is omitted.")
    out_dir <- file.path(fit_dir, "viz", "perturbation", "mixed_ploidy")
  }
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  validate_artifact_manifest(
    simulation_dir,
    "simulation_manifest.tsv",
    c("simulation_design.tsv", "burden_timecourse.tsv", "ploidy_summary_timecourse.tsv"),
    "Mixed-ploidy perturbation visualization"
  )
  validate_artifact_manifest(
    analysis_dir, "analysis_manifest.tsv", "endpoint_summary.tsv",
    "Mixed-ploidy perturbation visualization"
  )

  design <- read_required_tsv(file.path(simulation_dir, "simulation_design.tsv"))
  burden <- read_required_tsv(file.path(simulation_dir, "burden_timecourse.tsv"))
  ploidy <- read_required_tsv(file.path(simulation_dir, "ploidy_summary_timecourse.tsv"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_mixed_outputs(burden, ploidy, design, out_dir)
  figures <- list.files(out_dir, pattern = "\\.pdf$", full.names = FALSE)
  if (!length(figures)) stop("Visualization did not create any PDF figures.")
  rows <- lapply(figures, function(rel) data.frame(
    artifact = tools::file_path_sans_ext(rel), relative_path = rel, role = "figure",
    rows = NA_integer_, columns = NA_integer_, stringsAsFactors = FALSE
  ))
  write_artifact_manifest(out_dir, rows, "visualization_manifest.tsv")
  message("Done. Wrote visualization outputs to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  run_mixed_ploidy_perturbation_visualization()
}
