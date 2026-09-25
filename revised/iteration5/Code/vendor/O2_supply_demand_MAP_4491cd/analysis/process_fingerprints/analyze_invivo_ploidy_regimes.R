#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    frame_files <- Filter(
      nzchar,
      vapply(sys.frames(), function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      }, character(1))
    )
    if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
  }
})

source(file.path(SCRIPT_DIR, "process_fingerprint_utils.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "ploidy_regime_utils.R"), local = TRUE)

options(error = function() {
  traceback(2)
  quit(status = 1)
})

o2pr_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals[is.finite(vals)]
}

o2pr_update_manifest_for_ploidy <- function(manifest, traj_status, objective_deltas) {
  manifest$selected_objective <- manifest$objective
  manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  config_file <- if ("config_file" %in% names(manifest)) as.character(manifest$config_file) else rep(NA_character_, nrow(manifest))
  manifest$fit_config_available <- !is.na(config_file) & nzchar(config_file) & file.exists(config_file)
  manifest <- merge(manifest, traj_status, by = "seed_id", all.x = TRUE, sort = FALSE)
  if (!"trajectory_available" %in% names(manifest)) manifest$trajectory_available <- FALSE
  for (d in objective_deltas) {
    manifest[[paste0("objective_delta_le_", gsub("[^0-9]+", "p", d))]] <- manifest$fit_success %in% TRUE & manifest$delta_objective <= d
  }
  manifest
}

o2pr_seed_qc_tables <- function(manifest, inputs, out_dir, objective_deltas) {
  o2pr_write_tsv(manifest, file.path(out_dir, "tables", "seed_manifest.tsv"))
  qc <- data.frame(
    metric = c("discovered_seeds", "valid_seeds", "trajectory_available_seeds", paste0("n_delta_le_", objective_deltas)),
    value = c(
      nrow(manifest),
      sum(manifest$fit_success %in% TRUE, na.rm = TRUE),
      sum(manifest$trajectory_available %in% TRUE, na.rm = TRUE),
      vapply(objective_deltas, function(d) sum(manifest$fit_success %in% TRUE & manifest$delta_objective <= d, na.rm = TRUE), integer(1))
    ),
    stringsAsFactors = FALSE
  )
  o2pr_write_tsv(qc, file.path(out_dir, "tables", "seed_qc_summary.tsv"))
  excl <- manifest[!(manifest$fit_success %in% TRUE) | !(manifest$trajectory_available %in% TRUE), c("seed_id", "failure_reason", "trajectory_status"), drop = FALSE]
  o2pr_write_tsv(excl, file.path(out_dir, "tables", "seed_exclusion_log.tsv"))
  if (!is.null(inputs$boundary_long) && nrow(inputs$boundary_long)) {
    o2pr_write_tsv(inputs$boundary_long, file.path(out_dir, "tables", "parameter_boundary_long.tsv"))
  } else {
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "parameter_boundary_long.tsv"))
  }
}

o2pr_objective_by_regime <- function(manifest, labels, out_dir) {
  d <- merge(manifest, labels[, c("seed_id", "trajectory_regime")], by = "seed_id", all.x = TRUE)
  rows <- lapply(split(d, d$trajectory_regime), function(x) {
    data.frame(
      trajectory_regime = x$trajectory_regime[[1]],
      n_seed = nrow(x),
      selected_objective_median = stats::median(x$selected_objective, na.rm = TRUE),
      selected_objective_min = min(x$selected_objective, na.rm = TRUE),
      selected_objective_max = max(x$selected_objective, na.rm = TRUE),
      delta_objective_median = stats::median(x$delta_objective, na.rm = TRUE),
      boundary_risk_fraction = mean(x$boundary_risk %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  o2pr_write_tsv(do.call(rbind, rows), file.path(out_dir, "tables", "trajectory_regime_objective_summary.tsv"))
}

o2pr_write_fixed_summaries <- function(proc, out_dir) {
  intr <- proc$intrinsic_full
  if (!is.data.frame(intr) || !nrow(intr)) {
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "fixed_O2_attractor_summary.tsv"))
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "fixed_O2_simulation_summary.tsv"))
    return(invisible(NULL))
  }
  attractor <- intr[intr$module == "attractor", , drop = FALSE]
  o2pr_write_tsv(attractor, file.path(out_dir, "tables", "fixed_O2_attractor_summary.tsv"))
  sim_like <- intr[intr$module %in% c("growth_selection", "death_selection", "CIN", "buffering", "WGD"), , drop = FALSE]
  sim_like$note <- "core pipeline records fixed-O2 rate/kernel probes; long-run standardized time simulation is deferred unless explicitly extended"
  o2pr_write_tsv(sim_like, file.path(out_dir, "tables", "fixed_O2_simulation_summary.tsv"))
}

o2pr_continuous_axis_tables <- function(labels, scaled, out_dir) {
  if (is.null(scaled$combined_full) || !nrow(scaled$combined_full$wide)) {
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "continuous_axis_diagnostics.tsv"))
    return(invisible(NULL))
  }
  wide <- scaled$combined_full$wide
  mat <- as.matrix(wide[, setdiff(names(wide), "seed_id"), drop = FALSE])
  rownames(mat) <- wide$seed_id
  mat <- mat[rowSums(!is.finite(mat)) == 0, colSums(!is.finite(mat)) == 0, drop = FALSE]
  if (nrow(mat) < 3L || ncol(mat) < 2L) {
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "continuous_axis_diagnostics.tsv"))
    return(invisible(NULL))
  }
  pc <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
  lab <- labels[match(rownames(mat), labels$seed_id), , drop = FALSE]
  out <- data.frame(
    seed_id = rownames(mat),
    trajectory_regime = lab$trajectory_regime,
    trajectory_drop_score = lab$trajectory_drop_score,
    process_PC1 = pc$x[, 1],
    process_PC2 = if (ncol(pc$x) >= 2L) pc$x[, 2] else NA_real_,
    stringsAsFactors = FALSE
  )
  rows <- lapply(c("process_PC1", "process_PC2"), function(v) {
    ok <- is.finite(out$trajectory_drop_score) & is.finite(out[[v]])
    data.frame(
      relationship = paste0("trajectory_drop_score_vs_", v),
      n_seed = sum(ok),
      spearman_rho = if (sum(ok) >= 3L) suppressWarnings(stats::cor(out$trajectory_drop_score[ok], out[[v]][ok], method = "spearman")) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  o2pr_write_tsv(out, file.path(out_dir, "tables", "process_pc_scores_with_trajectory.tsv"))
  o2pr_write_tsv(do.call(rbind, rows), file.path(out_dir, "tables", "continuous_axis_diagnostics.tsv"))
}

o2pr_require_simulation_tables <- function(simulation_dir, required) {
  manifest_path <- file.path(simulation_dir, "simulation_manifest.tsv")
  if (!file.exists(manifest_path)) {
    stop("Missing simulation manifest: ", manifest_path, ". Run the process-fingerprint simulation producer first.")
  }
  manifest <- o2pr_read_tsv(manifest_path)
  if (!"relative_path" %in% names(manifest)) stop("Invalid simulation manifest: missing relative_path.")
  missing_entries <- setdiff(required, manifest$relative_path)
  missing_files <- required[!file.exists(file.path(simulation_dir, required))]
  if (length(missing_entries) || length(missing_files)) {
    stop(
      "Ploidy-regime analysis requires complete pre-existing simulation artifacts. Missing manifest entries: ",
      paste(missing_entries, collapse = ", "), "; missing files: ", paste(missing_files, collapse = ", ")
    )
  }
  invisible(manifest)
}

o2pr_cfg_from_metadata <- function(tab) {
  values <- setNames(as.character(tab$value), as.character(tab$metric))
  list(
    N_MIN = suppressWarnings(as.numeric(values[["Nmin"]])),
    N_MAX = suppressWarnings(as.numeric(values[["Nmax"]])),
    N_UNIT = suppressWarnings(as.numeric(values[["N_UNIT"]])),
    DT = suppressWarnings(as.numeric(values[["DT"]])),
    start_with = values[["start_with"]],
    boundary = values[["boundary_default"]]
  )
}

o2pr_analyze_trajectories <- function(curves, status, cfg, trajectory_end_day, mid_window, late_window) {
  nmin <- as.numeric(cfg$N_MIN %||% 22)
  long_features <- list()
  if (nrow(curves)) {
    for (idx in split(seq_len(nrow(curves)), paste(curves$seed_id, curves$cohort, sep = "||"))) {
      d <- curves[idx, , drop = FALSE]
      long_features[[length(long_features) + 1L]] <- o2pr_cohort_features(
        d$seed_id[[1]], d$cohort[[1]], d, mid_window, late_window, nmin
      )
    }
  }
  features_long <- if (length(long_features)) do.call(rbind, long_features) else data.frame()
  labels <- o2pr_make_trajectory_labels(
    curves, features_long, status, cfg, trajectory_end_day, mid_window, late_window
  )
  list(
    curves = curves,
    status = status,
    features_long = features_long,
    labels = labels,
    features_wide = o2pr_features_wide(features_long, labels)
  )
}

o2pr_write_analysis_manifest <- function(out_dir) {
  files <- list.files(file.path(out_dir, "tables"), full.names = FALSE)
  rows <- lapply(files, function(file) {
    path <- file.path(out_dir, "tables", file)
    tab <- tryCatch(o2pr_read_tsv(path), error = function(e) NULL)
    data.frame(
      artifact = tools::file_path_sans_ext(file), relative_path = file.path("tables", file),
      role = "analysis_table", rows = if (is.data.frame(tab)) nrow(tab) else NA_integer_,
      columns = if (is.data.frame(tab)) ncol(tab) else NA_integer_,
      path = normalizePath(path, mustWork = TRUE), exists = TRUE,
      stringsAsFactors = FALSE
    )
  })
  o2pr_write_tsv(do.call(rbind, rows), file.path(out_dir, "analysis_manifest.tsv"))
}

run_invivo_ploidy_regime_analysis_stage <- function(argv = o2ipa_parse_args()) {
  simulation_dir <- o2ipa_as_chr(argv$simulation_dir)
  if (!nzchar(simulation_dir) || !dir.exists(simulation_dir)) {
    stop("Missing or invalid --simulation_dir. Run simulation/process_fingerprints/generate_process_fingerprint_outputs.R first.")
  }
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  out_dir <- o2ipa_as_chr(argv$analysis_dir %||% argv$out_dir)
  if (!nzchar(out_dir)) stop("Missing required --analysis_dir (or --out_dir).")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)

  required <- file.path("tables", c(
    "seed_input_manifest.tsv", "parameter_values_long.tsv", "parameter_boundary_long.tsv",
    "fit_config_schema_summary.tsv", "trajectory_status.tsv", "trajectory_curves.tsv",
    "trajectory_density.tsv", "process_intrinsic_full_long.tsv",
    "process_intrinsic_18only_long.tsv", "process_exposure_full_long.tsv",
    "process_exposure_18only_long.tsv"
  ))
  o2pr_require_simulation_tables(simulation_dir, required)
  o2pr_mkdirs(out_dir)
  read_sim <- function(name) {
    path <- file.path(simulation_dir, "tables", name)
    if (!is.finite(file.info(path)$size) || file.info(path)$size <= 1L) return(data.frame())
    tryCatch(o2pr_read_tsv(path), error = function(e) data.frame())
  }
  inputs <- list(
    manifest = read_sim("seed_input_manifest.tsv"),
    params_long = read_sim("parameter_values_long.tsv"),
    boundary_long = read_sim("parameter_boundary_long.tsv")
  )
  cfg_meta <- read_sim("fit_config_schema_summary.tsv")
  cfg <- o2pr_cfg_from_metadata(cfg_meta)
  curves <- read_sim("trajectory_curves.tsv")
  status <- read_sim("trajectory_status.tsv")
  density <- read_sim("trajectory_density.tsv")
  proc <- list(
    intrinsic_full = read_sim("process_intrinsic_full_long.tsv"),
    intrinsic_18 = read_sim("process_intrinsic_18only_long.tsv"),
    exposure_full = read_sim("process_exposure_full_long.tsv"),
    exposure_18 = read_sim("process_exposure_18only_long.tsv")
  )

  objective_deltas <- o2pr_as_num_vec(argv$objective_deltas, c(2, 5, 10))
  if (!length(objective_deltas)) objective_deltas <- c(2, 5, 10)
  trajectory_end_day <- o2ipa_as_num(argv$trajectory_end_day, 1000)
  mid_window <- o2pr_as_num_vec(argv$mid_window, c(200, 700))
  late_window <- o2pr_as_num_vec(argv$late_window, c(900, 1000))
  n_bootstrap <- o2ipa_as_int(argv$n_bootstrap, 200L)
  n_permutations <- o2ipa_as_int(argv$n_permutations, 1000L)
  random_seed <- o2ipa_as_int(argv$random_seed, 20260623L)
  analysis_level <- o2ipa_as_chr(argv$analysis_level, "core")
  module_factorial <- o2ipa_as_bool(argv$module_factorial, FALSE)
  n_local_perturbations <- o2ipa_as_int(argv$n_local_perturbations, 0L)
  local_perturbation_sd <- o2ipa_as_num(argv$local_perturbation_sd, NA_real_)
  set.seed(random_seed)
  o2pr_write_tsv(data.frame(
    argument = c("simulation_dir", "objective_deltas", "trajectory_end_day", "mid_window", "late_window", "n_bootstrap", "n_permutations", "random_seed"),
    value = c(simulation_dir, paste(objective_deltas, collapse = ","), trajectory_end_day,
              paste(mid_window, collapse = ","), paste(late_window, collapse = ","),
              n_bootstrap, n_permutations, random_seed),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "tables", "analysis_run_arguments.tsv"))

  traj <- o2pr_analyze_trajectories(curves, status, cfg, trajectory_end_day, mid_window, late_window)
  traj$density <- density
  inputs$manifest <- o2pr_update_manifest_for_ploidy(inputs$manifest, traj$status, objective_deltas)
  o2pr_seed_qc_tables(inputs$manifest, inputs, out_dir, objective_deltas)
  o2pr_write_tsv(cfg_meta, file.path(out_dir, "tables", "fit_config_schema_summary.tsv"))
  o2pr_write_tsv(traj$curves, file.path(out_dir, "tables", "trajectory_curves_common_input.tsv"))
  o2pr_write_tsv(traj$features_long, file.path(out_dir, "tables", "trajectory_features_long.tsv"))
  o2pr_write_tsv(traj$features_wide, file.path(out_dir, "tables", "trajectory_features_wide.tsv"))
  o2pr_write_tsv(traj$labels, file.path(out_dir, "tables", "trajectory_regime_labels.tsv"))
  mat <- o2pr_common_grid_matrix(traj$curves, trajectory_end_day)$matrix
  stability <- o2pr_bootstrap_trajectory_stability(mat, n_bootstrap = min(n_bootstrap, 200L), random_seed = random_seed)
  o2pr_write_tsv(stability, file.path(out_dir, "tables", "trajectory_cluster_stability.tsv"))
  o2pr_objective_by_regime(inputs$manifest, traj$labels, out_dir)

  param <- o2pr_parameter_outputs(inputs, out_dir)
  o2pr_write_tsv(param$metadata, file.path(out_dir, "tables", "parameter_alias_transform_metadata.tsv"))
  o2pr_write_fixed_summaries(proc, out_dir)
  o2pr_oxygen_switch_diagnostics(proc, out_dir)
  scaled <- o2pr_scale_and_write_process(proc, out_dir)
  clusters <- o2pr_cluster_process(scaled, out_dir, n_bootstrap, random_seed)
  o2pr_concordance_outputs(traj$labels, clusters, scaled, out_dir, n_permutations, random_seed)
  o2pr_continuous_axis_tables(traj$labels, scaled, out_dir)
  o2pr_boundary_diagnostics(traj, cfg, out_dir)
  o2pr_placeholder_advanced(out_dir, analysis_level, module_factorial, n_local_perturbations, local_perturbation_sd)
  evidence <- o2pr_evidence_summary(out_dir, inputs$manifest, traj$labels, clusters, objective_deltas)
  o2pr_write_tsv(evidence, file.path(out_dir, "tables", "final_evidence_classification.tsv"))
  o2pr_write_analysis_manifest(out_dir)
  message("Completed ploidy-regime analysis: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  run_invivo_ploidy_regime_analysis_stage()
}
