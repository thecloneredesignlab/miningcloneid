#!/usr/bin/env Rscript

# Materialize raw post-fit process/ploidy/O2/CIN numerical artifacts.  This
# producer reads fit inputs and each seed's simulation/invivo contract, runs
# model probes, and writes no classifications, clusters, figures, or reports.

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
SIMULATION_UTIL_DIR <- .script_dir
source(file.path(SIMULATION_UTIL_DIR, "process_fingerprint_simulation_legacy_utils.R"), local = TRUE)
source(file.path(SIMULATION_UTIL_DIR, "ploidy_regime_simulation_utils.R"), local = TRUE)
source(file.path(.script_dir, "process_fingerprint_simulation_utils.R"), local = TRUE)

pfs_parse_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  values <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  values[is.finite(values)]
}

run_process_fingerprint_simulation <- function(argv = o2ipa_parse_args()) {
  run_dir <- o2ipa_as_chr(argv$run_dir)
  if (!nzchar(run_dir) || !dir.exists(run_dir)) stop("Missing or invalid --run_dir.")
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  out_dir <- o2ipa_as_chr(argv$simulation_dir %||% argv$out_dir, file.path(run_dir, "simulation", "process_fingerprints"))
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  o2ipa_mkdirs(out_dir)

  objective_source <- o2ipa_as_chr(argv$objective_source, "auto")
  max_seeds <- o2ipa_as_int(argv$max_seeds, NA_integer_)
  o2_grid <- pfs_parse_num_vec(argv$o2_grid, c(0.05, 0.1, 0.25, 0.5, 1, 2, 5))
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = objective_source)
  if (!is.na(max_seeds) && max_seeds > 0L) {
    valid <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
    valid <- valid[order(valid$objective, o2ipa_seed_number(valid$seed_id)), , drop = FALSE]
    keep <- valid$seed_id[seq_len(min(max_seeds, nrow(valid)))]
    inputs$manifest <- inputs$manifest[inputs$manifest$seed_id %in% keep, , drop = FALSE]
    inputs$params_long <- inputs$params_long[inputs$params_long$seed_id %in% keep, , drop = FALSE]
  }
  if (!nrow(inputs$manifest)) stop("No seeds were discovered in run_dir: ", run_dir)
  pfs_require_invivo_contract(inputs$manifest)

  cfg <- o2pr_first_seed_cfg(inputs$manifest)
  cfg_meta <- o2pr_cfg_metadata(cfg)
  raw_traj <- pfs_collect_raw_trajectories(inputs, cfg, SIMULATION_UTIL_DIR)
  if (any(!(raw_traj$status$trajectory_available %in% TRUE))) {
    bad <- raw_traj$status$seed_id[!(raw_traj$status$trajectory_available %in% TRUE)]
    stop("Failed to read required ploidy trajectories for: ", paste(bad, collapse = ", "))
  }
  dynamic <- pfs_dynamic_from_seed_simulations(inputs, raw_traj$curves)

  message("Sourcing the O2 model for raw process probes.")
  model_env <- o2ipa_source_model(SIMULATION_UTIL_DIR)
  static <- pfs_build_static_outputs(inputs, model_env)
  process <- o2pr_build_process_outputs(inputs, out_dir, model_env, cfg, o2_grid, force = TRUE)

  params <- inputs$params_long
  params$module <- vapply(params$parameter, o2ipa_param_module, character(1))
  o2ipa_write_tsv(inputs$manifest, file.path(out_dir, "tables", "seed_input_manifest.tsv"))
  o2ipa_write_tsv(params, file.path(out_dir, "tables", "parameter_values_long.tsv"))
  o2ipa_write_tsv(inputs$boundary_long %||% data.frame(), file.path(out_dir, "tables", "parameter_boundary_long.tsv"))
  o2ipa_write_tsv(cfg_meta, file.path(out_dir, "tables", "fit_config_schema_summary.tsv"))
  o2ipa_write_tsv(raw_traj$status, file.path(out_dir, "tables", "trajectory_status.tsv"))
  o2ipa_write_tsv(raw_traj$curves, file.path(out_dir, "tables", "trajectory_curves.tsv"))
  o2ipa_write_tsv(raw_traj$density, file.path(out_dir, "tables", "trajectory_density.tsv"))
  o2ipa_write_tsv(dynamic$burden, file.path(out_dir, "tables", "burden_timecourses.tsv"))
  o2ipa_write_tsv(dynamic$long, file.path(out_dir, "tables", "process_fingerprint_dynamic_long.tsv"))
  o2ipa_write_tsv(dynamic$unavailable, file.path(out_dir, "tables", "unavailable_dynamic_features.tsv"))
  o2ipa_write_tsv(static$full_long, file.path(out_dir, "tables", "process_fingerprint_static_full_long.tsv"))
  o2ipa_write_tsv(static$only18_long, file.path(out_dir, "tables", "process_fingerprint_static_18only_long.tsv"))
  o2ipa_write_tsv(static$reference, file.path(out_dir, "tables", "static_reference_parameter_set.tsv"))
  o2ipa_write_tsv(data.frame(
    key = c("run_dir", "objective_source", "max_seeds", "o2_grid"),
    value = c(run_dir, objective_source, max_seeds, paste(o2_grid, collapse = ",")),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "tables", "simulation_run_arguments.tsv"))

  files <- list.files(file.path(out_dir, "tables"), full.names = FALSE)
  files <- file.path("tables", files)
  pfs_write_manifest(out_dir, files)
  message("Completed process-fingerprint simulation: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  run_process_fingerprint_simulation()
}
