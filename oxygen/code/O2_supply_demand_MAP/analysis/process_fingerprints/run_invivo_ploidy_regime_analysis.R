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
  manifest$fit_config_available <- !is.na(manifest$config_file) & file.exists(manifest$config_file)
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

main <- function(argv = o2ipa_parse_args()) {
  repo_root <- o2ipa_repo_root(SCRIPT_DIR)
  run_dir <- normalizePath(o2ipa_as_chr(argv$run_dir), mustWork = FALSE)
  out_dir <- normalizePath(o2ipa_as_chr(argv$out_dir, file.path(repo_root, "oxygen", "results", "analysis")), mustWork = FALSE)
  if (!dir.exists(run_dir)) stop("run_dir does not exist: ", run_dir)
  o2pr_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "run_invivo_ploidy_regime_analysis.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)
  objective_source <- o2ipa_as_chr(argv$objective_source, "auto")
  objective_deltas <- o2pr_as_num_vec(argv$objective_deltas, c(2, 5, 10))
  if (!length(objective_deltas)) objective_deltas <- c(2, 5, 10)
  main_objective_delta <- o2ipa_as_num(argv$main_objective_delta, max(objective_deltas, na.rm = TRUE))
  trajectory_end_day <- o2ipa_as_num(argv$trajectory_end_day, 1000)
  mid_window <- o2pr_as_num_vec(argv$mid_window, c(200, 700))
  late_window <- o2pr_as_num_vec(argv$late_window, c(900, 1000))
  o2_grid <- o2pr_as_num_vec(argv$o2_grid, c(0.05, 0.1, 0.25, 0.5, 1, 2, 5))
  n_bootstrap <- o2ipa_as_int(argv$n_bootstrap, 200L)
  n_permutations <- o2ipa_as_int(argv$n_permutations, 1000L)
  workers <- o2ipa_as_int(argv$workers, 1L)
  random_seed <- o2ipa_as_int(argv$random_seed, 20260623L)
  max_seeds <- o2ipa_as_int(argv$max_seeds, NA_integer_)
  force <- o2ipa_as_bool(argv$force, FALSE)
  generate_viz <- o2ipa_as_bool(argv$generate_viz, FALSE)
  analysis_level <- o2ipa_as_chr(argv$analysis_level, "core")
  module_factorial <- o2ipa_as_bool(argv$module_factorial, FALSE)
  n_local_perturbations <- o2ipa_as_int(argv$n_local_perturbations, 0L)
  local_perturbation_sd <- o2ipa_as_num(argv$local_perturbation_sd, NA_real_)
  set.seed(random_seed)
  run_args <- data.frame(
    argument = c(
      "run_dir", "out_dir", "objective_source", "objective_deltas", "main_objective_delta",
      "trajectory_end_day", "mid_window", "late_window", "o2_grid", "n_bootstrap",
      "n_permutations", "workers", "max_seeds", "force", "generate_viz", "analysis_level",
      "module_factorial", "n_local_perturbations", "local_perturbation_sd", "random_seed"
    ),
    value = c(
      run_dir, out_dir, objective_source, paste(objective_deltas, collapse = ","),
      main_objective_delta, trajectory_end_day, paste(mid_window, collapse = ","),
      paste(late_window, collapse = ","), paste(o2_grid, collapse = ","), n_bootstrap,
      n_permutations, workers, max_seeds, force, generate_viz, analysis_level,
      module_factorial, n_local_perturbations, local_perturbation_sd, random_seed
    ),
    note = c(
      rep("", 11),
      "accepted for CLI compatibility; current implementation is single-process",
      rep("", 4),
      "accepted for CLI compatibility; factorial swaps are unavailable unless simulator/evaluator wiring is added",
      "accepted for CLI compatibility; local perturbations are unavailable unless simulator/evaluator wiring is added",
      "", ""
    ),
    stringsAsFactors = FALSE
  )
  o2pr_write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))

  message("Collecting seed inputs.")
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = objective_source)
  if (!is.na(max_seeds) && max_seeds > 0L) {
    valid <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
    valid_with_dirs <- valid[!is.na(valid$seed_dir) & nzchar(valid$seed_dir) & dir.exists(valid$seed_dir), , drop = FALSE]
    if (nrow(valid_with_dirs) > 0L) valid <- valid_with_dirs
    valid <- valid[order(valid$objective, o2ipa_seed_number(valid$seed_id)), , drop = FALSE]
    keep <- valid$seed_id[seq_len(min(max_seeds, nrow(valid)))]
    inputs$manifest <- inputs$manifest[inputs$manifest$seed_id %in% keep, , drop = FALSE]
    inputs$params_long <- inputs$params_long[inputs$params_long$seed_id %in% keep, , drop = FALSE]
    message("Using max_seeds subset: ", paste(keep, collapse = ", "))
  }

  cfg <- o2pr_first_seed_cfg(inputs$manifest)
  cfg_meta <- o2pr_cfg_metadata(cfg)
  o2pr_write_tsv(cfg_meta, file.path(out_dir, "tables", "fit_config_schema_summary.tsv"))

  message("Reading trajectories and assigning regimes.")
  traj <- o2pr_collect_trajectories(
    inputs = inputs,
    script_dir = SCRIPT_DIR,
    trajectory_end_day = trajectory_end_day,
    mid_window = mid_window,
    late_window = late_window,
    generate_viz = generate_viz,
    cfg = cfg
  )
  inputs$manifest <- o2pr_update_manifest_for_ploidy(inputs$manifest, traj$status, objective_deltas)
  o2pr_seed_qc_tables(inputs$manifest, inputs, out_dir, objective_deltas)
  o2pr_write_tsv(traj$curves, file.path(out_dir, "tables", "trajectory_curves_common_input.tsv"))
  o2pr_write_tsv(traj$features_long, file.path(out_dir, "tables", "trajectory_features_long.tsv"))
  o2pr_write_tsv(traj$features_wide, file.path(out_dir, "tables", "trajectory_features_wide.tsv"))
  o2pr_write_tsv(traj$labels, file.path(out_dir, "tables", "trajectory_regime_labels.tsv"))
  mat <- o2pr_common_grid_matrix(traj$curves, trajectory_end_day)$matrix
  stab <- o2pr_bootstrap_trajectory_stability(mat, n_bootstrap = min(n_bootstrap, 200L), random_seed = random_seed)
  o2pr_write_tsv(stab, file.path(out_dir, "tables", "trajectory_cluster_stability.tsv"))
  o2pr_objective_by_regime(inputs$manifest, traj$labels, out_dir)

  message("Writing parameter tables.")
  param <- o2pr_parameter_outputs(inputs, out_dir)
  o2pr_write_tsv(param$metadata, file.path(out_dir, "tables", "parameter_alias_transform_metadata.tsv"))

  message("Sourcing model and building process fingerprints.")
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  proc <- o2pr_build_process_outputs(inputs, out_dir, model_env, cfg, o2_grid, force = force)
  o2pr_write_fixed_summaries(proc, out_dir)
  o2pr_oxygen_switch_diagnostics(proc, out_dir)
  scaled <- o2pr_scale_and_write_process(proc, out_dir)

  message("Clustering process fingerprints.")
  clusters <- o2pr_cluster_process(scaled, out_dir, n_bootstrap = n_bootstrap, random_seed = random_seed)

  message("Computing trajectory-process concordance.")
  o2pr_concordance_outputs(traj$labels, clusters, scaled, out_dir, n_permutations = n_permutations, random_seed = random_seed)
  o2pr_continuous_axis_tables(traj$labels, scaled, out_dir)

  message("Writing boundary and advanced diagnostic tables.")
  o2pr_boundary_diagnostics(traj, cfg, out_dir)
  o2pr_placeholder_advanced(
    out_dir,
    analysis_level = analysis_level,
    module_factorial = module_factorial,
    n_local_perturbations = n_local_perturbations,
    local_perturbation_sd = local_perturbation_sd
  )
  evidence <- o2pr_evidence_summary(out_dir, inputs$manifest, traj$labels, clusters, objective_deltas = objective_deltas)
  o2pr_write_tsv(evidence, file.path(out_dir, "tables", "final_evidence_classification.tsv"))

  message("Writing figures and report.")
  o2pr_basic_figures(traj, scaled, clusters, traj$labels, out_dir)
  o2pr_report(
    out_dir = out_dir,
    run_dir = run_dir,
    repo_root = repo_root,
    manifest = inputs$manifest,
    cfg_meta = cfg_meta,
    labels = traj$labels,
    unavailable = c(
      "Objective recomputation for interpolation/module swaps was not invoked by the core pipeline.",
      "Alternative boundary and lower-grid sensitivity simulations are recorded as unavailable unless a dedicated safe simulation API is wired."
    )
  )

  message("Completed ploidy regime analysis: ", out_dir)
}

if (sys.nframe() == 0L) {
  main()
}
