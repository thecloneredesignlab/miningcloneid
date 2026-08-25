#!/usr/bin/env Rscript

# Pure visualization consumer for materialized in-vivo simulation products.
# It never performs simulation and never reads fitting or raw observation inputs.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

if (!nzchar(Sys.getenv("DISPLAY")) && isTRUE(capabilities("cairo"))) {
  options(bitmapType = "cairo")
}

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
  own_frame <- frame_files[basename(frame_files) == "viz_invivo_model_O2_supply_demand_MAP_results.R"]
  if (length(own_frame)) return(dirname(own_frame[[length(own_frame)]]))
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
})

SCRIPT_DIR <- normalizePath(.bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(SCRIPT_DIR, "invivo", "o2_supply_demand_map_invivo_plot_utils.R"), local = environment())
rm(.bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_bool <- o2sd_as_bool
clip <- o2sd_clip

.read_key_value_manifest <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing simulation manifest: ", path, ". ",
      "Run the in-vivo simulation producer before visualization."
    )
  }
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(tab), c("key", "value"))) {
    stop("Invalid simulation manifest (expected columns key,value): ", path)
  }
  if (anyDuplicated(tab$key)) stop("Invalid simulation manifest (duplicate keys): ", path)
  stats::setNames(as.character(tab$value), as.character(tab$key))
}

.manifest_num <- function(manifest, key, default = NA_real_) {
  value <- suppressWarnings(as.numeric(manifest[[key]] %||% default))
  if (!length(value) || !is.finite(value[[1L]])) as.numeric(default) else value[[1L]]
}

.manifest_num_vec <- function(manifest, key, default = numeric()) {
  o2sd_as_num_vec(manifest[[key]] %||% NULL, default)
}

.required_simulation_tables <- function(horizons) {
  base <- c(
    "burden_timecourse.tsv",
    "ploidy_timecourse.tsv",
    "burden_live_dead_decomposition.tsv",
    "ploidy_weighted_mean_timecourse.tsv",
    "missegregation_probability_timecourse.tsv",
    "terminal_ploidy_observed_vs_predicted.tsv",
    "o2_lag_timecourse.tsv",
    "predict_burden_vs_o2.tsv",
    "functional_curve_oxygen.tsv",
    "functional_curve_oxygen_multi_ploidy.tsv",
    "functional_curve_ploidy.tsv",
    "functional_curve_ploidy_by_o2.tsv",
    "population_average_cin_by_initial_cohort_horizons.tsv"
  )
  horizon_files <- unlist(lapply(horizons, function(day) {
    tag <- paste0("0_", as.integer(round(day)), "day")
    paste0(
      c(
        "predict_burden_", "predict_ploidy_", "predict_ploidy_weighted_mean_",
        "predict_chromosome_density_", "predict_death_ratio_",
        "predict_resource_death_fraction_"
      ),
      tag,
      ".tsv"
    )
  }), use.names = FALSE)
  c(base, horizon_files)
}

.validate_simulation_contract <- function(simulation_dir, manifest, horizons) {
  required_manifest <- c(
    "schema_version", "status", "fit_dir", "report_dt", "predict_report_dt",
    "predict_horizons", "start_with", "N_UNIT", "N_MIN", "N_MAX"
  )
  missing_keys <- setdiff(required_manifest, names(manifest))
  if (length(missing_keys)) {
    stop("Simulation manifest is missing key(s): ", paste(missing_keys, collapse = ", "))
  }
  if (!identical(manifest[["status"]], "complete")) {
    stop("Simulation manifest status is not complete: ", manifest[["status"]])
  }
  schema_path <- file.path(simulation_dir, "output_schema.tsv")
  if (!file.exists(schema_path)) stop("Missing simulation output schema: ", schema_path)
  schema <- utils::read.delim(schema_path, check.names = FALSE, stringsAsFactors = FALSE)
  schema_cols <- c("file", "role", "n_rows", "columns")
  if (!all(schema_cols %in% names(schema))) {
    stop("Invalid simulation output schema: ", schema_path)
  }
  required_files <- .required_simulation_tables(horizons)
  missing_files <- required_files[!file.exists(file.path(simulation_dir, required_files))]
  if (length(missing_files)) {
    stop(
      "Simulation output is incomplete; missing table(s): ",
      paste(missing_files, collapse = ", "), ". Run the simulation producer first."
    )
  }
  undeclared <- setdiff(required_files, schema$file)
  if (length(undeclared)) {
    stop("Simulation schema does not declare required table(s): ", paste(undeclared, collapse = ", "))
  }
  invisible(schema)
}

.cfg_from_manifest <- function(manifest) {
  fit_dir <- normalizePath(manifest[["fit_dir"]], mustWork = FALSE)
  list(
    fit_dir = fit_dir,
    fit_label = basename(fit_dir),
    report_dt = .manifest_num(manifest, "report_dt", 1),
    predict_report_dt = .manifest_num(manifest, "predict_report_dt", 1),
    start_with = manifest[["start_with"]] %||% "ploidy",
    N_UNIT = .manifest_num(manifest, "N_UNIT", 22),
    N_MIN = .manifest_num(manifest, "N_MIN", 1),
    N_MAX = .manifest_num(manifest, "N_MAX", 200),
    DT = .manifest_num(manifest, "DT", 1),
    rho_2N_min = .manifest_num(manifest, "rho_2N_min", 3.2e4),
    rho_2N_max = .manifest_num(manifest, "rho_2N_max", 5.6e4)
  )
}

.write_viz_manifest <- function(out_dir) {
  files <- list.files(out_dir, pattern = "[.](pdf|png)$", full.names = TRUE)
  files <- sort(files)
  info <- file.info(files)
  manifest <- data.frame(
    file = basename(files),
    type = tolower(tools::file_ext(files)),
    bytes = as.numeric(info$size),
    stringsAsFactors = FALSE
  )
  path <- file.path(out_dir, "viz_manifest.tsv")
  write.table(manifest, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

render_invivo_visualizations <- function(simulation_dir, out_dir, manifest_path, horizons = NULL) {
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  manifest <- .read_key_value_manifest(manifest_path)
  available_horizons <- .manifest_num_vec(manifest, "predict_horizons", numeric())
  available_horizons <- sort(unique(available_horizons[is.finite(available_horizons) & available_horizons > 0]))
  if (is.null(horizons)) horizons <- available_horizons
  horizons <- sort(unique(as.numeric(horizons)))
  horizons <- horizons[is.finite(horizons) & horizons > 0]
  unavailable <- setdiff(horizons, available_horizons)
  if (length(unavailable)) {
    stop("Requested horizon(s) are absent from the simulation manifest: ", paste(unavailable, collapse = ", "))
  }
  .validate_simulation_contract(simulation_dir, manifest, horizons)
  cfg <- .cfg_from_manifest(manifest)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  stale <- list.files(out_dir, pattern = "[.](pdf|png)$", full.names = TRUE)
  if (length(stale)) unlink(stale, force = TRUE)
  unlink(file.path(out_dir, "viz_manifest.tsv"), force = TRUE)

  render_invivo_observed_plots(simulation_dir, cfg, out_dir)
  plot_functional_response_curves(simulation_dir, cfg, out_dir)
  predict_results <- lapply(horizons, function(day) {
    plot_predict_horizon(
      simulation_dir = simulation_dir,
      cfg = cfg,
      out_dir = out_dir,
      horizon_day = day,
      report_dt = cfg$predict_report_dt
    )
  })
  plot_predict_burden_live_dead_decomposition_combined(
    predict_results,
    out_dir,
    resource_death_language()
  )
  manifest_out <- .write_viz_manifest(out_dir)
  invisible(list(out_dir = out_dir, manifest = manifest_out))
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  fit_dir <- argv$fit_dir %||% argv$run_dir %||% NULL
  manifest_arg <- argv$manifest %||% argv$simulation_manifest %||% NULL
  simulation_dir <- argv$simulation_dir %||% NULL

  if (!is.null(manifest_arg)) {
    manifest_path <- normalizePath(manifest_arg, mustWork = TRUE)
    if (is.null(simulation_dir)) simulation_dir <- dirname(manifest_path)
  } else {
    if (is.null(simulation_dir)) {
      if (is.null(fit_dir)) {
        stop("Provide --fit_dir/--run_dir, --simulation_dir, or --manifest/--simulation_manifest.")
      }
      simulation_dir <- file.path(fit_dir, "simulation", "invivo")
    }
    manifest_path <- file.path(simulation_dir, "simulation_manifest.tsv")
  }
  manifest <- .read_key_value_manifest(manifest_path)
  if (is.null(fit_dir)) fit_dir <- manifest[["fit_dir"]]
  out_dir <- argv$out_dir %||% file.path(fit_dir, "viz", "invivo")
  horizons <- if (!is.null(argv$predict_horizons)) {
    o2sd_as_num_vec(argv$predict_horizons, numeric())
  } else {
    NULL
  }
  if (!as_bool(argv$predict_plots, TRUE)) horizons <- numeric()

  result <- render_invivo_visualizations(
    simulation_dir = simulation_dir,
    out_dir = out_dir,
    manifest_path = manifest_path,
    horizons = horizons
  )
  message("In-vivo visualization written: ", result$out_dir)
  message("Visualization manifest: ", result$manifest)
}

if (sys.nframe() == 0L) main()
