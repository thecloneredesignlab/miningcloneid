#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files)) {
    dirname(frame_files[[length(frame_files)]])
  } else {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
  }
})
source(file.path(SCRIPT_DIR, "process_fingerprint_utils.R"), local = TRUE)

o2ipa_bind_rows_fill <- function(rows) {
  rows <- rows[vapply(rows, is.data.frame, logical(1))]
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  if (!length(cols)) return(data.frame())
  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (m in missing) x[[m]] <- rep(NA, nrow(x))
    x[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

o2ipa_add_cluster_source <- function(x, source) {
  if (!is.data.frame(x)) x <- data.frame()
  x$cluster_source <- rep(source, nrow(x))
  x
}

o2ipa_empty_module_centroids <- function() {
  data.frame(
    cluster_source = character(),
    k = integer(),
    cluster = integer(),
    module = character(),
    n_seed = integer(),
    n_features = integer(),
    mean_scaled_value = numeric(),
    mean_clustering_value = numeric(),
    stringsAsFactors = FALSE
  )
}

o2ipa_module_centroids <- function(scaled_long, membership, recommended_k, cluster_source = "static_full") {
  needed_long <- c("seed_id", "module", "feature_key", "retained_for_clustering", "scaled_value", "clustering_value")
  needed_membership <- c("seed_id", "k", "cluster")
  if (!is.data.frame(scaled_long) || !is.data.frame(membership) || is.na(recommended_k)) {
    return(o2ipa_empty_module_centroids())
  }
  if (!all(needed_long %in% names(scaled_long)) || !all(needed_membership %in% names(membership))) {
    return(o2ipa_empty_module_centroids())
  }

  retained_flag <- scaled_long$retained_for_clustering %in% c(TRUE, "TRUE", "true", "1", 1)
  d <- scaled_long[retained_flag, needed_long, drop = FALSE]
  memb <- membership[membership$k == recommended_k, c("seed_id", "cluster"), drop = FALSE]
  if (!nrow(d) || !nrow(memb)) return(o2ipa_empty_module_centroids())

  d$seed_id <- as.character(d$seed_id)
  memb$seed_id <- as.character(memb$seed_id)
  memb$cluster <- as.integer(memb$cluster)
  sf <- merge(d, memb, by = "seed_id", all = FALSE)
  if (!nrow(sf)) return(o2ipa_empty_module_centroids())

  mean_scaled <- aggregate(sf$scaled_value, by = list(cluster = sf$cluster, module = sf$module), FUN = mean, na.rm = TRUE)
  names(mean_scaled)[3] <- "mean_scaled_value"
  mean_clustering <- aggregate(sf$clustering_value, by = list(cluster = sf$cluster, module = sf$module), FUN = mean, na.rm = TRUE)
  names(mean_clustering)[3] <- "mean_clustering_value"
  n_seed <- aggregate(sf$seed_id, by = list(cluster = sf$cluster, module = sf$module), FUN = function(x) length(unique(x)))
  names(n_seed)[3] <- "n_seed"
  n_features <- aggregate(sf$feature_key, by = list(cluster = sf$cluster, module = sf$module), FUN = function(x) length(unique(x)))
  names(n_features)[3] <- "n_features"

  out <- Reduce(function(x, y) merge(x, y, by = c("cluster", "module"), all = TRUE), list(mean_scaled, mean_clustering, n_seed, n_features))
  out$cluster_source <- cluster_source
  out$k <- as.integer(recommended_k)
  out <- out[, c("cluster_source", "k", "cluster", "module", "n_seed", "n_features", "mean_scaled_value", "mean_clustering_value")]
  out[order(out$cluster, out$module), , drop = FALSE]
}

o2ipa_parameter_scaled_distance <- function(param_long) {
  p <- param_long
  p$fingerprint_scope <- "parameter"
  p$feature <- p$parameter
  p$feature_type <- "signed"
  p$raw_value <- p$transformed_value
  p$module <- vapply(p$parameter, o2ipa_param_module, character(1))
  o2ipa_scale_long_features(p[, c("seed_id", "fingerprint_scope", "module", "feature", "raw_value", "feature_type")], missing_feature_max = 0, missing_seed_max = 0)
}

o2ipa_seed_qc_outputs <- function(inputs, out_dir, objective_deltas) {
  manifest <- inputs$manifest
  manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  o2ipa_write_tsv(manifest, file.path(out_dir, "tables", "seed_manifest.tsv"))
  qc <- data.frame(
    metric = c("discovered_seeds", "valid_seeds", paste0("n_delta_le_", objective_deltas)),
    value = c(nrow(manifest), sum(manifest$fit_success %in% TRUE), vapply(objective_deltas, function(d) sum(manifest$fit_success %in% TRUE & manifest$delta_objective <= d, na.rm = TRUE), integer(1))),
    stringsAsFactors = FALSE
  )
  o2ipa_write_tsv(qc, file.path(out_dir, "tables", "seed_qc_summary.tsv"))
  excl <- manifest[!(manifest$fit_success %in% TRUE), c("seed_id", "failure_reason"), drop = FALSE]
  o2ipa_write_tsv(excl, file.path(out_dir, "tables", "seed_exclusion_log.tsv"))
  if (!is.null(inputs$boundary_long) && nrow(inputs$boundary_long)) {
    o2ipa_write_tsv(inputs$boundary_long, file.path(out_dir, "tables", "parameter_boundary_long.tsv"))
  } else {
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", "parameter_boundary_long.tsv"))
  }
  manifest
}

o2ipa_cluster_qc <- function(manifest, cluster, out_dir, prefix = "static_full") {
  if (!nrow(cluster$membership) || is.na(cluster$recommended_k)) {
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", "cluster_qc_summary.tsv"))
    return(invisible(NULL))
  }
  memb <- cluster$membership[cluster$membership$k == cluster$recommended_k, , drop = FALSE]
  df <- merge(manifest, memb, by = "seed_id", all.x = FALSE)
  rows <- lapply(sort(unique(df$cluster)), function(cl) {
    d <- df[df$cluster == cl, , drop = FALSE]
    data.frame(
      cluster_source = prefix,
      k = cluster$recommended_k,
      cluster = cl,
      n_seed = nrow(d),
      objective_median = stats::median(d$objective, na.rm = TRUE),
      objective_iqr = stats::IQR(d$objective, na.rm = TRUE),
      objective_min = min(d$objective, na.rm = TRUE),
      objective_max = max(d$objective, na.rm = TRUE),
      delta_objective_median = stats::median(d$delta_objective, na.rm = TRUE),
      boundary_risk_fraction = mean(d$boundary_risk %in% TRUE, na.rm = TRUE),
      convergence_failure_fraction = mean(!(d$fit_success %in% TRUE), na.rm = TRUE),
      visualization_failure_fraction = mean(!(d$visualization_available %in% TRUE), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  qc <- do.call(rbind, rows)
  o2ipa_write_tsv(qc, file.path(out_dir, "tables", "cluster_qc_summary.tsv"))
  o2ipa_write_tsv(qc[, c("cluster_source", "k", "cluster", "boundary_risk_fraction")], file.path(out_dir, "tables", "cluster_boundary_enrichment.tsv"))
  o2ipa_write_tsv(qc[, c("cluster_source", "k", "cluster", "objective_median", "objective_iqr", "objective_min", "objective_max")], file.path(out_dir, "tables", "cluster_objective_summary.tsv"))
}

o2ipa_leave_one_module_out <- function(static_scaled_long, recommended_k, base_membership, out_dir) {
  if (is.na(recommended_k) || !nrow(static_scaled_long)) {
    out <- data.frame()
    o2ipa_write_tsv(out, file.path(out_dir, "tables", "leave_one_module_out_cluster_stability.tsv"))
    return(out)
  }
  modules <- sort(unique(static_scaled_long$module[static_scaled_long$retained_for_clustering %in% TRUE]))
  base <- base_membership[base_membership$k == recommended_k, , drop = FALSE]
  rows <- lapply(modules, function(mod) {
    d <- static_scaled_long[static_scaled_long$module != mod, , drop = FALSE]
    wide <- o2ipa_long_to_wide(d[d$retained_for_clustering %in% TRUE, c("seed_id", "feature_key", "clustering_value")], "feature_key", "clustering_value")
    dm <- o2ipa_distance(wide)
    if (is.null(dm) || nrow(dm) < recommended_k + 1L) {
      return(data.frame(left_out_module = mod, pairwise_assignment_concordance = NA_real_, status = "insufficient_features", stringsAsFactors = FALSE))
    }
    hc <- stats::hclust(stats::as.dist(dm), method = "ward.D2")
    cl <- stats::cutree(hc, k = recommended_k)
    common <- intersect(base$seed_id, names(cl))
    if (length(common) < 3L) {
      agree <- NA_real_
    } else {
      idx <- utils::combn(common, 2)
      bcl <- setNames(base$cluster, base$seed_id)
      agree <- mean((bcl[idx[1, ]] == bcl[idx[2, ]]) == (cl[idx[1, ]] == cl[idx[2, ]]))
    }
    data.frame(left_out_module = mod, pairwise_assignment_concordance = agree, status = "ok", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  o2ipa_write_tsv(out, file.path(out_dir, "tables", "leave_one_module_out_cluster_stability.tsv"))
  out
}

o2ipa_concordance <- function(memb_a, memb_b, k_a, k_b, name) {
  a <- memb_a[memb_a$k == k_a, , drop = FALSE]
  b <- memb_b[memb_b$k == k_b, , drop = FALSE]
  common <- intersect(a$seed_id, b$seed_id)
  if (length(common) < 3L) {
    return(data.frame(comparison = name, n_seed = length(common), pairwise_concordance = NA_real_, stringsAsFactors = FALSE))
  }
  ca <- setNames(a$cluster, a$seed_id)[common]
  cb <- setNames(b$cluster, b$seed_id)[common]
  idx <- utils::combn(common, 2)
  data.frame(comparison = name, n_seed = length(common), pairwise_concordance = mean((ca[idx[1, ]] == ca[idx[2, ]]) == (cb[idx[1, ]] == cb[idx[2, ]])), stringsAsFactors = FALSE)
}

o2ipa_require_simulation_tables <- function(simulation_dir, required) {
  manifest_path <- file.path(simulation_dir, "simulation_manifest.tsv")
  if (!file.exists(manifest_path)) {
    stop("Missing simulation manifest: ", manifest_path, ". Run the process-fingerprint simulation producer first.")
  }
  manifest <- o2ipa_read_tsv(manifest_path)
  if (!"relative_path" %in% names(manifest)) stop("Invalid simulation manifest: missing relative_path.")
  missing_entries <- setdiff(required, manifest$relative_path)
  missing_files <- required[!file.exists(file.path(simulation_dir, required))]
  if (length(missing_entries) || length(missing_files)) {
    stop(
      "Process-fingerprint analysis requires complete pre-existing simulation artifacts. Missing manifest entries: ",
      paste(missing_entries, collapse = ", "), "; missing files: ", paste(missing_files, collapse = ", ")
    )
  }
  invisible(manifest)
}

o2ipa_parameter_outputs_from_simulation <- function(params_long, out_dir) {
  params_long$module <- vapply(params_long$parameter, o2ipa_param_module, character(1))
  params_long$transformed_value <- mapply(
    o2ipa_transform_parameter_value, params_long$parameter, params_long$value
  )
  o2ipa_write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))
  raw <- o2ipa_params_wide(params_long, "value")
  raw_out <- data.frame(seed_id = rownames(raw), raw, check.names = FALSE)
  rownames(raw_out) <- NULL
  trans <- o2ipa_params_wide(params_long, "transformed_value")
  trans_out <- data.frame(seed_id = rownames(trans), trans, check.names = FALSE)
  rownames(trans_out) <- NULL
  metadata <- o2ipa_transform_metadata()
  o2ipa_write_tsv(raw_out, file.path(out_dir, "tables", "parameter_matrix_raw.tsv"))
  o2ipa_write_tsv(trans_out, file.path(out_dir, "tables", "parameter_matrix_transformed.tsv"))
  o2ipa_write_tsv(metadata, file.path(out_dir, "tables", "parameter_transform_metadata.tsv"))
  list(raw = raw_out, transformed = trans_out, metadata = metadata, long = params_long)
}

o2ipa_write_analysis_manifest <- function(out_dir) {
  files <- list.files(file.path(out_dir, "tables"), full.names = FALSE)
  rows <- lapply(files, function(file) {
    path <- file.path(out_dir, "tables", file)
    tab <- tryCatch(o2ipa_read_tsv(path), error = function(e) NULL)
    data.frame(
      artifact = tools::file_path_sans_ext(file),
      relative_path = file.path("tables", file),
      role = "analysis_table",
      rows = if (is.data.frame(tab)) nrow(tab) else NA_integer_,
      columns = if (is.data.frame(tab)) ncol(tab) else NA_integer_,
      path = normalizePath(path, mustWork = TRUE),
      exists = TRUE,
      stringsAsFactors = FALSE
    )
  })
  o2ipa_write_tsv(do.call(rbind, rows), file.path(out_dir, "analysis_manifest.tsv"))
}

run_invivo_process_fingerprint_analysis <- function(argv = o2ipa_parse_args()) {
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
    "process_fingerprint_static_full_long.tsv", "process_fingerprint_static_18only_long.tsv",
    "process_fingerprint_dynamic_long.tsv", "unavailable_dynamic_features.tsv"
  ))
  o2ipa_require_simulation_tables(simulation_dir, required)
  o2ipa_mkdirs(out_dir)
  read_sim <- function(name) {
    path <- file.path(simulation_dir, "tables", name)
    if (!is.finite(file.info(path)$size) || file.info(path)$size <= 1L) return(data.frame())
    tryCatch(o2ipa_read_tsv(path), error = function(e) data.frame())
  }
  inputs <- list(
    manifest = read_sim("seed_input_manifest.tsv"),
    params_long = read_sim("parameter_values_long.tsv"),
    boundary_long = read_sim("parameter_boundary_long.tsv")
  )
  static_full <- read_sim("process_fingerprint_static_full_long.tsv")
  static_18 <- read_sim("process_fingerprint_static_18only_long.tsv")
  dynamic_long <- read_sim("process_fingerprint_dynamic_long.tsv")
  unavailable <- read_sim("unavailable_dynamic_features.tsv")

  objective_deltas <- suppressWarnings(as.numeric(o2ipa_split_csv(argv$objective_deltas, c("2", "5", "10"))))
  objective_deltas <- objective_deltas[is.finite(objective_deltas)]
  n_bootstrap <- o2ipa_as_int(argv$n_bootstrap, 200L)
  random_seed <- o2ipa_as_int(argv$random_seed, 20260623L)
  inputs$manifest <- o2ipa_seed_qc_outputs(inputs, out_dir, objective_deltas)

  param <- o2ipa_parameter_outputs_from_simulation(inputs$params_long, out_dir)
  param_scaled <- o2ipa_parameter_scaled_distance(param$long)
  o2ipa_write_tsv(param_scaled$metadata, file.path(out_dir, "tables", "parameter_feature_metadata.tsv"))
  d_parameter <- o2ipa_distance(param_scaled$wide)
  o2ipa_write_distance(d_parameter, "parameter", out_dir)

  static_full_scaled <- o2ipa_scale_long_features(static_full, missing_feature_max = 0, missing_seed_max = 0)
  static_18_scaled <- o2ipa_scale_long_features(static_18, missing_feature_max = 0, missing_seed_max = 0)
  o2ipa_write_tsv(static_full_scaled$long, file.path(out_dir, "tables", "process_fingerprint_static_full_scaled_long.tsv"))
  o2ipa_write_tsv(static_18_scaled$long, file.path(out_dir, "tables", "process_fingerprint_static_18only_scaled_long.tsv"))
  o2ipa_write_tsv(static_full_scaled$wide, file.path(out_dir, "tables", "process_fingerprint_static_scaled.tsv"))
  o2ipa_write_tsv(static_18_scaled$wide, file.path(out_dir, "tables", "process_fingerprint_static_18only_scaled.tsv"))
  o2ipa_write_tsv(static_full_scaled$metadata, file.path(out_dir, "tables", "process_feature_metadata.tsv"))
  o2ipa_write_tsv(static_full_scaled$missingness, file.path(out_dir, "tables", "process_feature_missingness.tsv"))
  d_static <- o2ipa_distance(static_full_scaled$wide)
  d_static18 <- o2ipa_distance(static_18_scaled$wide)
  o2ipa_write_distance(d_static, "static_full", out_dir)
  o2ipa_write_distance(d_static18, "static_18only", out_dir)

  dyn_raw_wide <- if (nrow(dynamic_long)) {
    o2ipa_long_to_wide(transform(dynamic_long, feature_key = paste(module, feature, sep = "||")), "feature_key", "raw_value")
  } else data.frame()
  o2ipa_write_tsv(dynamic_long, file.path(out_dir, "tables", "process_fingerprint_dynamic_long.tsv"))
  o2ipa_write_tsv(dyn_raw_wide, file.path(out_dir, "tables", "process_fingerprint_dynamic_wide_raw.tsv"))
  o2ipa_write_tsv(unavailable, file.path(out_dir, "tables", "unavailable_dynamic_features.tsv"))
  dynamic_scaled <- o2ipa_scale_long_features(dynamic_long, missing_feature_max = 0, missing_seed_max = 0)
  o2ipa_write_tsv(dynamic_scaled$wide, file.path(out_dir, "tables", "process_fingerprint_dynamic_scaled.tsv"))
  d_dynamic <- o2ipa_distance(dynamic_scaled$wide)
  o2ipa_write_distance(d_dynamic, "dynamic", out_dir)

  dlist <- Filter(Negate(is.null), list(parameter = d_parameter, static_full = d_static, static_18only = d_static18, dynamic = d_dynamic))
  dist_cor <- if (length(dlist) >= 2L) o2ipa_distance_correlations(dlist) else data.frame()
  o2ipa_write_tsv(dist_cor, file.path(out_dir, "tables", "distance_matrix_correlations.tsv"))
  pair_diag <- o2ipa_pairwise_diagnostics(d_parameter, d_static, d_static18, d_dynamic, inputs$manifest)
  if (is.list(pair_diag)) {
    o2ipa_write_tsv(pair_diag$pairs, file.path(out_dir, "tables", "pairwise_regime_diagnostics.tsv"))
    o2ipa_write_tsv(pair_diag$thresholds, file.path(out_dir, "tables", "distance_thresholds.tsv"))
  } else {
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", "pairwise_regime_diagnostics.tsv"))
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", "distance_thresholds.tsv"))
  }

  cl_static <- o2ipa_cluster_distance(d_static, n_bootstrap, random_seed, static_full_scaled$wide, static_full_scaled$metadata)
  cl_static18 <- o2ipa_cluster_distance(d_static18, n_bootstrap, random_seed, static_18_scaled$wide, static_18_scaled$metadata)
  cl_dynamic <- o2ipa_cluster_distance(d_dynamic, n_bootstrap, random_seed, dynamic_scaled$wide, dynamic_scaled$metadata)
  diagnostics <- o2ipa_bind_rows_fill(list(
    o2ipa_add_cluster_source(cl_static$diagnostics, "static_full"),
    o2ipa_add_cluster_source(cl_static18$diagnostics, "static_18only"),
    o2ipa_add_cluster_source(cl_dynamic$diagnostics, "dynamic")
  ))
  o2ipa_write_tsv(diagnostics, file.path(out_dir, "tables", "cluster_k_diagnostics.tsv"))
  o2ipa_write_tsv(cl_static$membership, file.path(out_dir, "tables", "cluster_membership_static_full.tsv"))
  o2ipa_write_tsv(cl_static18$membership, file.path(out_dir, "tables", "cluster_membership_static_18only.tsv"))
  o2ipa_write_tsv(cl_dynamic$membership, file.path(out_dir, "tables", "cluster_membership_dynamic.tsv"))
  o2ipa_write_tsv(cl_static$stability, file.path(out_dir, "tables", "cluster_stability.tsv"))
  medoids <- o2ipa_bind_rows_fill(list(
    o2ipa_add_cluster_source(cl_static$medoids, "static_full"),
    o2ipa_add_cluster_source(cl_static18$medoids, "static_18only"),
    o2ipa_add_cluster_source(cl_dynamic$medoids, "dynamic")
  ))
  o2ipa_write_tsv(medoids, file.path(out_dir, "tables", "cluster_medoids.tsv"))
  if (!is.null(cl_static$consensus)) {
    o2ipa_write_tsv(data.frame(seed_id = rownames(cl_static$consensus), cl_static$consensus, check.names = FALSE), file.path(out_dir, "tables", "cluster_consensus_matrix.tsv"))
  } else o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", "cluster_consensus_matrix.tsv"))

  o2ipa_write_tsv(o2ipa_concordance(cl_static$membership, cl_static18$membership, cl_static$recommended_k, cl_static18$recommended_k, "static_full_vs_static_18only"), file.path(out_dir, "tables", "cluster_concordance_full_vs_18only.tsv"))
  o2ipa_write_tsv(o2ipa_concordance(cl_static$membership, cl_dynamic$membership, cl_static$recommended_k, cl_dynamic$recommended_k, "static_vs_dynamic"), file.path(out_dir, "tables", "cluster_concordance_static_vs_dynamic.tsv"))
  o2ipa_cluster_qc(inputs$manifest, cl_static, out_dir, "static_full")
  o2ipa_leave_one_module_out(static_full_scaled$long, cl_static$recommended_k, cl_static$membership, out_dir)
  o2ipa_write_tsv(o2ipa_module_centroids(static_full_scaled$long, cl_static$membership, cl_static$recommended_k, "static_full"), file.path(out_dir, "tables", "cluster_process_module_centroid_profile.tsv"))
  o2ipa_write_analysis_manifest(out_dir)
  message("Completed process-fingerprint analysis: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  run_invivo_process_fingerprint_analysis()
}
