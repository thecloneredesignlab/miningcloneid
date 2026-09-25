#!/usr/bin/env Rscript

# Generate post-fit mixed-ploidy perturbation trajectories.  This producer
# writes numerical simulation artifacts only; analysis and plotting are separate
# consumers under analysis/perturbation and vis/perturbation.

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
load_o2sd_perturbation_model(environment())

run_mixed_ploidy_perturbation_simulation <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- resolve_path_value(argv$fit_dir %||% argv$run_dir, getwd())
  if (is.null(fit_dir)) stop("Missing required argument: --fit_dir=/path/to/seed_dir")
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  cfg <- prepare_sim_cfg(read_fit_config(fit_dir), argv)
  run_params <- read_run_params(fit_dir, cfg)

  out_dir <- resolve_path_value(argv$simulation_dir %||% argv$out_dir, getwd())
  if (is.null(out_dir)) out_dir <- file.path(fit_dir, "simulation", "perturbation", "mixed_ploidy")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  initial_burden_cells <- as_num(argv$initial_burden_cells, as_num(argv$initial_live_cells, 1e5))
  initial_burden_values <- parse_num_list(argv$initial_burden_values, c(1e5, 2.5e5, 5e5, 1e6, 2e6))
  trigger_burden_values <- parse_num_list(argv$trigger_burden_values, c(5e5, 1e6, 2e6, 5e6, 1e7))
  o2_values <- parse_num_list(argv$o2_values, c(0.5, 1, 2, 3, 4, 5))
  pmis_base_values <- parse_num_list(argv$pmis_base_values, c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1))
  p_wgd_values <- parse_num_list(argv$p_wgd_values, as.numeric(run_params$p_wgd) * c(1e-2, 1e-1, 1, 10, 100))
  if (as_bool(argv$smoke, FALSE)) {
    initial_burden_values <- head(initial_burden_values, 1)
    trigger_burden_values <- head(trigger_burden_values, 1)
    o2_values <- head(o2_values, 1)
    pmis_base_values <- head(pmis_base_values, 1)
    p_wgd_values <- head(p_wgd_values, 1)
  }

  horizon_day <- as_num(argv$horizon_day, 1000)
  report_dt <- as_num(argv$report_dt, 1.0)
  trigger_check_dt <- as_num(argv$trigger_check_dt, report_dt)
  init_mean <- as_num(argv$initial_ploidy_mean, 3.0)
  init_sd <- as_num(argv$initial_ploidy_sd, 0.4)
  init_min <- as_num(argv$initial_ploidy_min, 1.5)
  init_max <- as_num(argv$initial_ploidy_max, 6.0)

  design <- build_design_rows(
    run_params, initial_burden_values, o2_values, p_wgd_values,
    pmis_base_values, trigger_burden_values, initial_burden_cells
  )
  design$initial_ploidy_mean <- init_mean
  design$initial_ploidy_sd <- init_sd
  design$initial_ploidy_min <- init_min
  design$initial_ploidy_max <- init_max
  design$fit_dir <- fit_dir

  message("Running ", nrow(design), " mixed-ploidy scenarios. Output: ", out_dir)
  burden_rows <- vector("list", nrow(design))
  ploidy_rows <- vector("list", nrow(design))
  design_rows <- vector("list", nrow(design))
  for (i in seq_len(nrow(design))) {
    dr <- design[i, , drop = FALSE]
    message("[", i, "/", nrow(design), "] ", dr$scenario_id)
    init_state <- make_continuous_ploidy_state(
      cfg, as.numeric(dr$initial_burden_cells), init_mean, init_sd, init_min, init_max
    )
    if (identical(as.character(dr$experiment), "pmiss_triggered_treatment")) {
      res <- run_triggered_pmiss_scenario(
        run_params, cfg, init_state, horizon_day, report_dt, trigger_check_dt, dr
      )
      design_rows[[i]] <- res$design
      burden_rows[[i]] <- res$burden
      ploidy_rows[[i]] <- res$ploidy
    } else {
      res <- run_static_scenario(run_params, cfg, init_state, horizon_day, report_dt, dr)
      design_rows[[i]] <- res$burden[1, names(dr), drop = FALSE]
      burden_rows[[i]] <- res$burden
      ploidy_rows[[i]] <- res$ploidy
    }
  }

  design_out <- bind_rows(design_rows)
  burden_all <- bind_rows(burden_rows)
  ploidy_all <- bind_rows(ploidy_rows)
  ploidy_summary <- summarise_ploidy_timecourse(ploidy_all)

  write_tsv(design_out, file.path(out_dir, "simulation_design.tsv"))
  write_tsv(burden_all, file.path(out_dir, "burden_timecourse.tsv"))
  write_tsv(ploidy_summary, file.path(out_dir, "ploidy_summary_timecourse.tsv"))
  write_tsv_gz(ploidy_all, file.path(out_dir, "ploidy_distribution_timecourse.tsv.gz"))
  write_artifact_manifest(
    out_dir,
    list(
      artifact_manifest_row("design", "simulation_design.tsv", "simulation_design", design_out),
      artifact_manifest_row("burden", "burden_timecourse.tsv", "simulated_trajectory", burden_all),
      artifact_manifest_row("ploidy_summary", "ploidy_summary_timecourse.tsv", "simulated_trajectory", ploidy_summary),
      artifact_manifest_row("ploidy_distribution", "ploidy_distribution_timecourse.tsv.gz", "simulated_state_distribution", ploidy_all)
    ),
    "simulation_manifest.tsv"
  )
  message("Done. Wrote simulation outputs to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  run_mixed_ploidy_perturbation_simulation()
}
