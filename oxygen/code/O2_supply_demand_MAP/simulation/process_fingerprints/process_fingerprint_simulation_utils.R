#!/usr/bin/env Rscript

# Simulation-layer helpers for post-fit process fingerprints.  These functions
# may read fit artifacts and the already-materialized simulation/invivo contract;
# they do not cluster, classify, plot, or report.

pfs_write_manifest <- function(out_dir, files, role = "simulated_numeric_artifact") {
  files <- unique(as.character(files))
  paths <- file.path(out_dir, files)
  missing <- files[!file.exists(paths)]
  if (length(missing)) stop("Cannot write simulation manifest; missing artifacts: ", paste(missing, collapse = ", "))
  rows <- lapply(seq_along(files), function(i) {
    tab <- tryCatch(o2ipa_read_tsv(paths[[i]]), error = function(e) NULL)
    data.frame(
      artifact = tools::file_path_sans_ext(gsub("[/\\\\]", "_", files[[i]])),
      relative_path = files[[i]],
      role = role,
      rows = if (is.data.frame(tab)) nrow(tab) else NA_integer_,
      columns = if (is.data.frame(tab)) ncol(tab) else NA_integer_,
      path = normalizePath(paths[[i]], mustWork = TRUE),
      exists = TRUE,
      stringsAsFactors = FALSE
    )
  })
  o2ipa_write_tsv(do.call(rbind, rows), file.path(out_dir, "simulation_manifest.tsv"))
}

pfs_require_invivo_contract <- function(manifest) {
  failures <- list()
  for (i in seq_len(nrow(manifest))) {
    if (!(manifest$fit_success[[i]] %in% TRUE)) next
    seed <- manifest$seed_id[[i]]
    seed_dir <- manifest$seed_dir[[i]]
    sim_dir <- file.path(seed_dir, "simulation", "invivo")
    required <- c("simulation_manifest.tsv")
    trajectory <- o2pr_find_trajectory_path(seed_dir)
    burden_candidates <- file.path(sim_dir, c(
      "predict_burden_0_1000day.tsv", "predict_burden_0_300day.tsv",
      "burden_timecourse.tsv", "predict_burden_vs_o2.tsv"
    ))
    burden <- burden_candidates[file.exists(burden_candidates)]
    missing <- required[!file.exists(file.path(sim_dir, required))]
    if (is.na(trajectory)) missing <- c(missing, "ploidy trajectory table")
    if (!length(burden)) missing <- c(missing, "burden/O2 trajectory table")
    if (length(missing)) {
      failures[[seed]] <- paste0(seed, ": ", paste(missing, collapse = ", "))
    }
  }
  if (length(failures)) {
    stop(
      "Process-fingerprint simulation requires pre-existing simulation/invivo artifacts for every valid seed: ",
      paste(unlist(failures, use.names = FALSE), collapse = "; "),
      ". Run simulation/invivo/generate_invivo_simulation_outputs.R first."
    )
  }
  invisible(TRUE)
}

pfs_collect_raw_trajectories <- function(inputs, cfg, script_dir) {
  curves <- list()
  density <- list()
  status <- list()
  for (i in seq_len(nrow(inputs$manifest))) {
    seed <- inputs$manifest$seed_id[[i]]
    seed_dir <- inputs$manifest$seed_dir[[i]]
    tr <- o2pr_read_trajectory_seed(seed, seed_dir, cfg, script_dir, generate_viz = FALSE)
    if (nrow(tr$curve)) curves[[seed]] <- tr$curve
    den <- o2pr_read_density_seed(seed, seed_dir)
    if (nrow(den)) density[[seed]] <- den
    status[[seed]] <- data.frame(
      seed_id = seed,
      trajectory_available = identical(tr$status, "ok"),
      trajectory_status = tr$status,
      trajectory_source_file = tr$source_file,
      stringsAsFactors = FALSE
    )
  }
  list(
    curves = if (length(curves)) do.call(rbind, curves) else data.frame(),
    density = if (length(density)) do.call(rbind, density) else data.frame(),
    status = do.call(rbind, status)
  )
}

pfs_read_seed_burden <- function(seed_id, seed_dir) {
  sim_dir <- file.path(seed_dir, "simulation", "invivo")
  candidates <- file.path(sim_dir, c(
    "predict_burden_0_1000day.tsv", "predict_burden_0_300day.tsv",
    "burden_timecourse.tsv", "predict_burden_vs_o2.tsv"
  ))
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) return(data.frame())
  tab <- o2ipa_read_tsv(hits[[1]])
  day_col <- intersect(c("day", "report_day"), names(tab))
  burden_col <- intersect(c("pred_burden_cells", "pred_burden_total_cells", "burden_value", "pred_burden_volume_mm3"), names(tab))
  if (!length(day_col) || !length(burden_col)) return(data.frame())
  cohort <- if ("cohort" %in% names(tab)) as.character(tab$cohort) else "all"
  out <- data.frame(
    seed_id = seed_id,
    cohort = cohort,
    day = suppressWarnings(as.numeric(tab[[day_col[[1]]]])),
    burden_value = suppressWarnings(as.numeric(tab[[burden_col[[1]]]])),
    stringsAsFactors = FALSE
  )
  o2_col <- intersect(c("pred_o2_pct", "O2_eff", "oxygen_pct"), names(tab))
  out$pred_o2_pct <- if (length(o2_col)) suppressWarnings(as.numeric(tab[[o2_col[[1]]]])) else NA_real_
  out <- out[is.finite(out$day) & is.finite(out$burden_value), , drop = FALSE]
  out$source_file <- normalizePath(hits[[1]], mustWork = TRUE)
  out
}

pfs_dynamic_from_seed_simulations <- function(inputs, curves) {
  rows <- list()
  burden_long <- list()
  for (i in seq_len(nrow(inputs$manifest))) {
    seed <- inputs$manifest$seed_id[[i]]
    b <- pfs_read_seed_burden(seed, inputs$manifest$seed_dir[[i]])
    if (nrow(b)) burden_long[[seed]] <- b
    for (idx in split(seq_len(nrow(b)), b$cohort)) {
      d <- b[idx, , drop = FALSE]
      d <- d[order(d$day), , drop = FALSE]
      prefix <- paste0("cohort_", d$cohort[[1]], "_")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "terminal_burden"), tail(d$burden_value, 1), "positive", "simulation_invivo")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "burden_AUC"), o2ipa_auc(d$day, d$burden_value), "positive", "simulation_invivo")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "maximum_burden"), max(d$burden_value, na.rm = TRUE), "positive", "simulation_invivo")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "time_to_burden_double"), o2ipa_crossing(d$day, d$burden_value, 2 * d$burden_value[[1]])$value, "time", "simulation_invivo")
    }
  }
  if (nrow(curves)) {
    for (idx in split(seq_len(nrow(curves)), paste(curves$seed_id, curves$cohort, sep = "||"))) {
      d <- curves[idx, , drop = FALSE]
      d <- d[order(d$day), , drop = FALSE]
      prefix <- paste0("cohort_", d$cohort[[1]], "_")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(d$seed_id[[1]], "dynamic", "dynamic_ploidy", paste0(prefix, "terminal_mean_chromosome"), tail(d$trajectory_value, 1), "positive", "simulation_invivo")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(d$seed_id[[1]], "dynamic", "dynamic_ploidy", paste0(prefix, "maximum_mean_chromosome"), max(d$trajectory_value, na.rm = TRUE), "positive", "simulation_invivo")
    }
  }
  unavailable <- data.frame(
    feature = c("ploidy_entropy_AUC", "terminal_ploidy_entropy", "WGD_or_high_ploidy_lineage_fraction"),
    reason = "not_available_in_materialized_simulation_outputs",
    stringsAsFactors = FALSE
  )
  list(
    long = if (length(rows)) do.call(rbind, rows) else data.frame(),
    burden = if (length(burden_long)) do.call(rbind, burden_long) else data.frame(),
    unavailable = unavailable
  )
}

pfs_build_static_outputs <- function(inputs, model_env) {
  manifest <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
  if (!nrow(manifest)) stop("No valid seeds after input QC.")
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  best_seed <- manifest$seed_id[which.min(manifest$objective)]
  target <- o2ipa_target_params()
  ref <- as.numeric(param_mat[best_seed, , drop = TRUE])
  names(ref) <- colnames(param_mat)
  ref_full <- c(ref, o2_min = 0, o2_S0_upper_bound = 5, o2_Nref = 1e6)
  full_rows <- list()
  only18_rows <- list()
  for (seed in manifest$seed_id) {
    pvec <- as.numeric(param_mat[seed, , drop = TRUE]); names(pvec) <- colnames(param_mat)
    full_params <- c(pvec, o2_min = 0, o2_S0_upper_bound = 5, o2_Nref = 1e6)
    only18_params <- ref_full; only18_params[target] <- pvec[target]
    full_rows[[seed]] <- o2ipa_build_static_one(seed, full_params, model_env, scope = "static_full")
    only18_rows[[seed]] <- o2ipa_build_static_one(seed, only18_params, model_env, scope = "static_18only")
  }
  list(
    full_long = do.call(rbind, full_rows),
    only18_long = do.call(rbind, only18_rows),
    reference = data.frame(parameter = names(ref_full), value = as.numeric(ref_full), reference_seed = best_seed, stringsAsFactors = FALSE)
  )
}
