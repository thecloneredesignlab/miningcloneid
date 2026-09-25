#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = TRUE)
OXYGEN_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", ".."), mustWork = TRUE)

source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R"), local = environment())
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_invitro_utils.R"),
  local = environment(),
  chdir = TRUE
)
module_paths <- c(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_invitro_postfit_io_utils.R"),
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_invitro_postfit_response_utils.R"),
  file.path(.script_dir, "population", "invitro_population_simulation_utils.R"),
  file.path(.script_dir, "ploidy", "invitro_ploidy_simulation_utils.R"),
  file.path(.script_dir, "cin", "invitro_cin_simulation_utils.R"),
  file.path(.script_dir, "o2", "invitro_o2_functional_response_utils.R")
)
missing_modules <- module_paths[!file.exists(module_paths)]
if (length(missing_modules)) stop("Missing in-vitro simulation module(s): ", paste(missing_modules, collapse = ", "))
for (module_path in module_paths) {
  sys.source(module_path, envir = environment(), chdir = TRUE)
}
rm(module_path, module_paths, missing_modules)
rm(.script_dir)

main <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir %||% stop(
    "Usage: run_invitro_simulation_outputs.R --fit_dir=/abs/seed [--simulation_dir=/abs/output]",
    call. = FALSE
  )
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  simulation_dir <- argv$simulation_dir %||% argv$out_dir %||%
    file.path(fit_dir, "simulation", "invitro")
  dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)

  fit_result_path <- argv$fit_result %||% file.path(fit_dir, "fit_result.rds")
  if (!file.exists(fit_result_path)) stop("Missing fit_result.rds: ", fit_result_path)
  fit_result_path <- normalizePath(fit_result_path, mustWork = TRUE)
  best_params_path <- argv$best_params %||% file.path(fit_dir, "best_params.tsv")
  if (!file.exists(best_params_path)) best_params_path <- NULL

  fit_result <- readRDS(fit_result_path)
  payload <- ivt_sim_extract_payload(
    fit_result,
    fit_dir = fit_dir,
    best_params_path = best_params_path
  )
  tables <- c(
    ivt_sim_collect_population_tables(payload$components),
    ivt_sim_collect_ploidy_tables(payload$components)
  )
  tables$invitro_missegregation_probability_timecourse <-
    ivt_sim_compute_missegregation_timecourse(
      tables$invitro_distribution_summary,
      payload$run_params
    )
  tables$invitro_optimizer_population <-
    ivt_sim_extract_optimizer_population(fit_result)
  functional_context <- ivt_sim_build_functional_response_context(
    payload$run_params,
    payload$cfg
  )
  tables <- c(
    tables,
    ivt_sim_compute_o2_response_tables(
      payload$run_params,
      payload$cfg,
      context = functional_context
    ),
    ivt_sim_compute_ploidy_response_tables(
      payload$run_params,
      payload$cfg,
      context = functional_context
    )
  )

  for (name in names(tables)) {
    ivt_sim_write_tsv(tables[[name]], file.path(simulation_dir, paste0(name, ".tsv")))
  }
  schema <- ivt_sim_schema(tables)
  ivt_sim_write_tsv(schema, file.path(simulation_dir, "invitro_simulation_schema.tsv"))
  manifest <- data.frame(
    key = c(
      "fit_dir", "simulation_dir", "fit_result", "best_params", "generated_at",
      "model_r", "model_cpp", "table_count", paste0("rows.", names(tables))
    ),
    value = c(
      fit_dir,
      simulation_dir,
      fit_result_path,
      if (is.null(best_params_path)) "" else normalizePath(best_params_path, mustWork = TRUE),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      normalizePath(file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R"), mustWork = TRUE),
      normalizePath(file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.cpp"), mustWork = TRUE),
      as.character(length(tables)),
      as.character(vapply(tables, nrow, integer(1)))
    ),
    stringsAsFactors = FALSE
  )
  ivt_sim_write_tsv(manifest, file.path(simulation_dir, "invitro_simulation_manifest.tsv"))
  message("In-vitro simulation outputs written to: ", simulation_dir)
  invisible(simulation_dir)
}

if (sys.nframe() == 0L) main()
