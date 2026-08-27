#!/usr/bin/env Rscript

# Materialize all in-vivo simulation products from an existing fitted parameter
# set.  This entry point intentionally produces tables only; visualization is a
# separate consumer of the manifest written here.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Matrix))

.bootstrap_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own_frame <- frame_files[basename(frame_files) == "generate_invivo_simulation_outputs.R"]
  if (length(own_frame)) return(dirname(own_frame[[length(own_frame)]]))
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
})

SCRIPT_DIR <- normalizePath(.bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())

model_path <- file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R")
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
source(model_path, local = environment())
source(file.path(SCRIPT_DIR, "o2_supply_demand_map_invivo_simulation_utils.R"), local = environment())
rm(.bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_bool <- o2sd_as_bool
clip <- o2sd_clip

.resolve_data_dir <- function(argv) {
  if (!is.null(argv$data_dir)) {
    return(normalizePath(argv$data_dir, mustWork = TRUE))
  }
  candidates <- c(
    file.path(WORKFLOW_ROOT, "..", "..", "data", "InVivoData_Gemcitabine"),
    file.path(WORKFLOW_ROOT, "..", "..", "..", "data", "InVivoData_Gemcitabine")
  )
  hit <- candidates[dir.exists(candidates)]
  normalizePath(if (length(hit)) hit[[1L]] else candidates[[1L]], mustWork = FALSE)
}

.apply_simulation_overrides <- function(cfg, argv) {
  if (!is.null(argv$dose_zero_only)) {
    cfg$dose_zero_only <- as_bool(argv$dose_zero_only, cfg$dose_zero_only)
  }
  if (!is.null(argv$truncate_at_treatment)) {
    cfg$truncate_at_treatment <- as_bool(argv$truncate_at_treatment, cfg$truncate_at_treatment)
  }
  if (!is.null(argv$ploidy_at_harvest)) {
    cfg$ploidy_at_harvest <- as_bool(argv$ploidy_at_harvest, cfg$ploidy_at_harvest)
  }
  if (!is.null(argv$max_scenarios)) {
    cfg$max_scenarios <- as_num(argv$max_scenarios, cfg$max_scenarios)
  }
  cfg
}

.manifest_value <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  paste(as.character(x), collapse = ",")
}

.standardize_observations <- function(observed_tables) {
  burden <- observed_tables$burden_timecourse %>%
    filter(is.finite(obs_burden)) %>%
    transmute(
      observation_type = "tumor_burden",
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      value = as.numeric(obs_burden),
      weight = 1,
      unit = "mm3"
    )
  terminal <- observed_tables$terminal_ploidy_observed_vs_predicted
  if (nrow(terminal)) {
    terminal <- terminal %>%
      filter(as.character(source) == "Observed") %>%
      transmute(
        observation_type = "terminal_ploidy",
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(target_day),
        value = as.numeric(endpoint_value),
        weight = as.numeric(weight),
        unit = as.character(endpoint_mode)
      )
  } else {
    terminal <- burden[0, , drop = FALSE]
  }
  bind_rows(burden, terminal)
}

generate_invivo_simulation_outputs <- function(
  fit_dir,
  out_dir = file.path(fit_dir, "simulation", "invivo"),
  data_dir = NULL,
  report_dt = 1,
  predict_horizons = c(100, 300, 1000),
  predict_report_dt = report_dt,
  argv = list()
) {
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  cfg_path <- file.path(fit_dir, "fit_config.rds")
  params_path <- file.path(fit_dir, "best_params.tsv")
  missing_fit_inputs <- c(cfg_path, params_path)[!file.exists(c(cfg_path, params_path))]
  if (length(missing_fit_inputs)) {
    stop("Missing fitted-input file(s): ", paste(missing_fit_inputs, collapse = ", "))
  }

  report_dt <- as.numeric(report_dt)
  predict_report_dt <- as.numeric(predict_report_dt)
  if (!is.finite(report_dt) || report_dt <= 0) stop("report_dt must be > 0")
  if (!is.finite(predict_report_dt) || predict_report_dt <= 0) {
    stop("predict_report_dt must be > 0")
  }
  predict_horizons <- sort(unique(as.numeric(predict_horizons)))
  predict_horizons <- predict_horizons[is.finite(predict_horizons) & predict_horizons > 0]

  if (is.null(data_dir)) data_dir <- .resolve_data_dir(argv)
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- resolve_terminal_ploidy_path(data_dir)
  if (!file.exists(dt_path)) stop("Missing in-vivo burden workbook: ", dt_path)
  if (!file.exists(ploidy_path)) stop("Missing in-vivo terminal ploidy table: ", ploidy_path)

  cfg <- normalize_cfg_for_simulation(readRDS(cfg_path))
  cfg <- .apply_simulation_overrides(cfg, argv)
  run_params <- read_run_params(fit_dir, cfg = cfg)
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  written <- list()
  roles <- character()

  write_contract <- function(filename, value, role) {
    value <- as.data.frame(value, stringsAsFactors = FALSE)
    path <- file.path(out_dir, filename)
    write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE)
    written[[filename]] <<- value
    roles[[filename]] <<- role
    invisible(path)
  }

  observed <- generate_invivo_observed_tables(
    run_params = run_params,
    scenarios = scenarios,
    cfg = cfg,
    report_dt = report_dt
  )
  observed_names <- c(
    burden_timecourse = "burden_timecourse.tsv",
    ploidy_timecourse = "ploidy_timecourse.tsv",
    burden_live_dead_decomposition = "burden_live_dead_decomposition.tsv",
    ploidy_weighted_mean_timecourse = "ploidy_weighted_mean_timecourse.tsv",
    missegregation_probability_timecourse = "missegregation_probability_timecourse.tsv",
    terminal_ploidy_observed_vs_predicted = "terminal_ploidy_observed_vs_predicted.tsv",
    o2_lag_timecourse = "o2_lag_timecourse.tsv",
    predict_burden_vs_o2 = "predict_burden_vs_o2.tsv"
  )
  for (nm in names(observed_names)) {
    write_contract(observed_names[[nm]], observed[[nm]], "observed_fit_simulation")
  }

  functional <- generate_invivo_functional_response_tables(run_params, cfg)
  for (nm in names(functional)) {
    write_contract(paste0(nm, ".tsv"), functional[[nm]], "functional_response_simulation")
  }

  cin_rows <- list()
  for (horizon_day in predict_horizons) {
    horizon <- generate_invivo_horizon_tables(
      run_params = run_params,
      scenarios = scenarios,
      cfg = cfg,
      horizon_day = horizon_day,
      report_dt = predict_report_dt
    )
    tag <- paste0("0_", as.integer(round(horizon_day)), "day")
    write_contract(paste0("predict_burden_", tag, ".tsv"), horizon$burden, "forecast_simulation")
    write_contract(paste0("predict_ploidy_", tag, ".tsv"), horizon$ploidy, "forecast_simulation")
    write_contract(
      paste0("predict_ploidy_weighted_mean_", tag, ".tsv"),
      horizon$ploidy_weighted_mean,
      "forecast_simulation"
    )
    write_contract(
      paste0("predict_chromosome_density_", tag, ".tsv"),
      horizon$chromosome_density,
      "forecast_simulation"
    )
    write_contract(paste0("predict_death_ratio_", tag, ".tsv"), horizon$death_ratio, "forecast_simulation")
    write_contract(
      paste0("predict_resource_death_fraction_", tag, ".tsv"),
      horizon$resource_death_fraction,
      "forecast_simulation"
    )
    cin_rows[[tag]] <- horizon$population_average_cin
  }
  write_contract(
    "population_average_cin_by_initial_cohort_horizons.tsv",
    bind_rows(cin_rows),
    "forecast_cin_simulation"
  )

  write_contract("observations.tsv", .standardize_observations(observed), "observed_input_contract")
  params_used <- data.frame(
    parameter = names(run_params),
    value = vapply(run_params, .manifest_value, character(1)),
    stringsAsFactors = FALSE
  )
  write_contract("parameters_used.tsv", params_used, "fitted_parameter_provenance")

  schema <- bind_rows(lapply(names(written), function(filename) {
    value <- written[[filename]]
    data.frame(
      file = filename,
      role = unname(roles[[filename]]),
      n_rows = nrow(value),
      columns = paste(names(value), collapse = ","),
      stringsAsFactors = FALSE
    )
  }))
  schema_path <- file.path(out_dir, "output_schema.tsv")
  write.table(schema, schema_path, sep = "\t", quote = FALSE, row.names = FALSE)

  manifest <- data.frame(
    key = c(
      "schema_version", "status", "fit_dir", "simulation_dir", "data_dir",
      "report_dt", "predict_report_dt", "predict_horizons", "start_with",
      "endpoint_mode", "N_UNIT", "N_MIN", "N_MAX", "DT", "rho_2N_min",
      "rho_2N_max", "table_count", "created_at"
    ),
    value = c(
      "o2sd-invivo-simulation-v1", "complete", fit_dir, out_dir, data_dir,
      .manifest_value(report_dt), .manifest_value(predict_report_dt),
      .manifest_value(predict_horizons), .manifest_value(cfg$start_with),
      .manifest_value(cfg$start_with), .manifest_value(cfg$N_UNIT),
      .manifest_value(cfg$N_MIN), .manifest_value(cfg$N_MAX), .manifest_value(cfg$DT),
      .manifest_value(as_num(argv$rho_2N_min, 3.2e4)),
      .manifest_value(as_num(argv$rho_2N_max, 5.6e4)),
      as.character(nrow(schema)), format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(out_dir, "simulation_manifest.tsv")
  write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(
    fit_dir = fit_dir,
    simulation_dir = out_dir,
    manifest = manifest_path,
    schema = schema_path,
    tables = schema$file
  ))
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  fit_dir <- argv$fit_dir %||% argv$run_dir
  if (is.null(fit_dir) || !nzchar(trimws(fit_dir))) {
    stop("Provide --fit_dir=<fit-result-directory> (alias: --run_dir).")
  }
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  out_dir <- argv$out_dir %||% argv$simulation_dir %||% file.path(fit_dir, "simulation", "invivo")
  report_dt <- as_num(argv$report_dt, 1)
  predict_report_dt <- as_num(argv$predict_report_dt, report_dt)
  predict_horizons <- o2sd_as_num_vec(argv$predict_horizons, c(100, 300, 1000))

  result <- generate_invivo_simulation_outputs(
    fit_dir = fit_dir,
    out_dir = out_dir,
    data_dir = argv$data_dir %||% NULL,
    report_dt = report_dt,
    predict_horizons = predict_horizons,
    predict_report_dt = predict_report_dt,
    argv = argv
  )
  message("In-vivo simulation tables written: ", result$simulation_dir)
  message("Manifest: ", result$manifest)
}

if (sys.nframe() == 0L) main()
