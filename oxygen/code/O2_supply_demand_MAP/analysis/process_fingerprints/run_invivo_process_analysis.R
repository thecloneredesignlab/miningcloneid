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

o2ipa_collect_parameter_tables <- function(manifest) {
  paths <- unique(na.omit(manifest$seed_dir))
  out <- list()
  for (d in paths[seq_len(min(10L, length(paths)))]) {
    path <- file.path(d, "parameter_table.csv")
    if (!file.exists(path)) next
    tab <- tryCatch(o2ipa_read_csv_or_tsv(path), error = function(e) NULL)
    if (is.null(tab)) next
    tab$source_file <- normalizePath(path, mustWork = FALSE)
    out[[length(out) + 1L]] <- tab
  }
  out
}

o2ipa_parameter_outputs <- function(inputs, out_dir) {
  params_long <- inputs$params_long
  params_long$module <- vapply(params_long$parameter, o2ipa_param_module, character(1))
  params_long$transformed_value <- mapply(o2ipa_transform_parameter_value, params_long$parameter, params_long$value)
  o2ipa_write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))

  raw <- o2ipa_params_wide(params_long, "value")
  raw_out <- data.frame(seed_id = rownames(raw), raw, check.names = FALSE)
  rownames(raw_out) <- NULL
  o2ipa_write_tsv(raw_out, file.path(out_dir, "tables", "parameter_matrix_raw.tsv"))

  trans <- o2ipa_params_wide(params_long, "transformed_value")
  trans_out <- data.frame(seed_id = rownames(trans), trans, check.names = FALSE)
  rownames(trans_out) <- NULL
  o2ipa_write_tsv(trans_out, file.path(out_dir, "tables", "parameter_matrix_transformed.tsv"))

  param_tables <- o2ipa_collect_parameter_tables(inputs$manifest)
  meta <- o2ipa_transform_metadata(param_tables)
  o2ipa_write_tsv(meta, file.path(out_dir, "tables", "parameter_transform_metadata.tsv"))
  list(raw = raw_out, transformed = trans_out, metadata = meta, long = params_long)
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

o2ipa_build_static_outputs <- function(inputs, out_dir, model_env, force = FALSE) {
  manifest <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
  if (!nrow(manifest)) stop("No valid seeds after input QC.")
  params_long <- inputs$params_long
  seeds <- manifest$seed_id
  best_seed <- manifest$seed_id[which.min(manifest$objective)]
  target <- o2ipa_target_params()
  all_param_names <- unique(params_long$parameter)
  param_mat <- o2ipa_params_wide(params_long, "value")
  ref <- as.numeric(param_mat[best_seed, , drop = TRUE])
  names(ref) <- colnames(param_mat)
  ref_full <- ref
  ref_full <- c(ref_full, o2_min = 0, o2_S0_upper_bound = 5, o2_Nref = 1e6)
  ref_df <- data.frame(parameter = names(ref_full), value = as.numeric(ref_full), reference_seed = best_seed, stringsAsFactors = FALSE)
  o2ipa_write_tsv(ref_df, file.path(out_dir, "tables", "reference_parameter_set.tsv"))

  full_rows <- list()
  only18_rows <- list()
  for (seed in seeds) {
    cache_file <- file.path(out_dir, "cache", paste0(seed, "_static_fingerprint.rds"))
    if (file.exists(cache_file) && !isTRUE(force)) {
      cached <- readRDS(cache_file)
      full_rows[[seed]] <- cached$full
      only18_rows[[seed]] <- cached$only18
      next
    }
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    full_params <- c(pvec, o2_min = 0, o2_S0_upper_bound = 5, o2_Nref = 1e6)
    only18_params <- ref_full
    only18_params[target] <- pvec[target]
    full <- o2ipa_build_static_one(seed, full_params, model_env, scope = "static_full")
    only18 <- o2ipa_build_static_one(seed, only18_params, model_env, scope = "static_18only")
    saveRDS(list(full = full, only18 = only18), cache_file)
    full_rows[[seed]] <- full
    only18_rows[[seed]] <- only18
  }
  full_long <- do.call(rbind, full_rows)
  only18_long <- do.call(rbind, only18_rows)
  o2ipa_write_tsv(full_long, file.path(out_dir, "tables", "process_fingerprint_static_full_long.tsv"))
  o2ipa_write_tsv(only18_long, file.path(out_dir, "tables", "process_fingerprint_static_18only_long.tsv"))
  list(full_long = full_long, only18_long = only18_long, reference_seed = best_seed)
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

main <- function(argv = o2ipa_parse_args()) {
  script_dir <- SCRIPT_DIR
  run_dir <- normalizePath(o2ipa_as_chr(argv$run_dir), mustWork = FALSE)
  out_dir <- normalizePath(o2ipa_as_chr(argv$out_dir, o2ipa_default_out_dir(script_dir)), mustWork = FALSE)
  if (!dir.exists(run_dir)) stop("run_dir does not exist: ", run_dir)
  o2ipa_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "run_invivo_process_analysis.log")
  sink(log_file, split = TRUE)
  on.exit(sink(), add = TRUE)
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)

  objective_source <- o2ipa_as_chr(argv$objective_source, "auto")
  objective_deltas <- suppressWarnings(as.numeric(o2ipa_split_csv(argv$objective_deltas, c("2", "5", "10"))))
  objective_deltas <- objective_deltas[is.finite(objective_deltas)]
  main_objective_delta <- o2ipa_as_num(argv$main_objective_delta, 10)
  n_bootstrap <- o2ipa_as_int(argv$n_bootstrap, 200L)
  random_seed <- o2ipa_as_int(argv$random_seed, 20260623L)
  max_seeds <- o2ipa_as_int(argv$max_seeds, NA_integer_)
  force <- o2ipa_as_bool(argv$force, FALSE)

  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = objective_source)
  manifest <- o2ipa_seed_qc_outputs(inputs, out_dir, objective_deltas)
  inputs$manifest <- manifest
  if (!is.na(max_seeds) && max_seeds > 0L) {
    valid <- manifest[manifest$fit_success %in% TRUE, , drop = FALSE]
    valid <- valid[order(valid$objective, o2ipa_seed_number(valid$seed_id)), , drop = FALSE]
    keep <- valid$seed_id[seq_len(min(max_seeds, nrow(valid)))]
    inputs$manifest <- manifest[manifest$seed_id %in% keep, , drop = FALSE]
    inputs$params_long <- inputs$params_long[inputs$params_long$seed_id %in% keep, , drop = FALSE]
    message("Using max_seeds subset: ", paste(keep, collapse = ", "))
  }

  param <- o2ipa_parameter_outputs(inputs, out_dir)
  param_scaled <- o2ipa_parameter_scaled_distance(param$long)
  o2ipa_write_tsv(param_scaled$metadata, file.path(out_dir, "tables", "parameter_feature_metadata.tsv"))
  d_parameter <- o2ipa_distance(param_scaled$wide)
  o2ipa_write_distance(d_parameter, "parameter", out_dir)

  message("Sourcing model helpers.")
  model_env <- o2ipa_source_model(script_dir)
  static <- o2ipa_build_static_outputs(inputs, out_dir, model_env, force = force)

  static_full_scaled <- o2ipa_scale_long_features(static$full_long, missing_feature_max = 0, missing_seed_max = 0)
  static_18_scaled <- o2ipa_scale_long_features(static$only18_long, missing_feature_max = 0, missing_seed_max = 0)
  o2ipa_write_tsv(static_full_scaled$long, file.path(out_dir, "tables", "process_fingerprint_static_full_scaled_long.tsv"))
  o2ipa_write_tsv(static_18_scaled$long, file.path(out_dir, "tables", "process_fingerprint_static_18only_scaled_long.tsv"))
  o2ipa_write_tsv(static_full_scaled$wide, file.path(out_dir, "tables", "process_fingerprint_static_scaled.tsv"))
  o2ipa_write_tsv(static_18_scaled$wide, file.path(out_dir, "tables", "process_fingerprint_static_18only_scaled.tsv"))
  feature_meta <- static_full_scaled$metadata
  o2ipa_write_tsv(feature_meta, file.path(out_dir, "tables", "process_feature_metadata.tsv"))
  o2ipa_write_tsv(static_full_scaled$missingness, file.path(out_dir, "tables", "process_feature_missingness.tsv"))

  d_static <- o2ipa_distance(static_full_scaled$wide)
  d_static18 <- o2ipa_distance(static_18_scaled$wide)
  o2ipa_write_distance(d_static, "static_full", out_dir)
  o2ipa_write_distance(d_static18, "static_18only", out_dir)

  dynamic <- o2ipa_dynamic_fingerprints(run_dir, inputs$manifest$seed_id)
  o2ipa_write_tsv(dynamic$long, file.path(out_dir, "tables", "process_fingerprint_dynamic_long.tsv"))
  dyn_raw_wide <- if (nrow(dynamic$long)) o2ipa_long_to_wide(transform(dynamic$long, feature_key = paste(module, feature, sep = "||")), "feature_key", "raw_value") else data.frame()
  o2ipa_write_tsv(dyn_raw_wide, file.path(out_dir, "tables", "process_fingerprint_dynamic_wide_raw.tsv"))
  o2ipa_write_tsv(dynamic$unavailable, file.path(out_dir, "tables", "unavailable_dynamic_features.tsv"))
  dynamic_scaled <- o2ipa_scale_long_features(dynamic$long, missing_feature_max = 0, missing_seed_max = 0)
  o2ipa_write_tsv(dynamic_scaled$wide, file.path(out_dir, "tables", "process_fingerprint_dynamic_scaled.tsv"))
  d_dynamic <- o2ipa_distance(dynamic_scaled$wide)
  o2ipa_write_distance(d_dynamic, "dynamic", out_dir)

  dlist <- Filter(Negate(is.null), list(parameter = d_parameter, static_full = d_static, static_18only = d_static18, dynamic = d_dynamic))
  dist_cor <- if (length(dlist) >= 2L) o2ipa_distance_correlations(dlist) else data.frame()
  o2ipa_write_tsv(dist_cor, file.path(out_dir, "tables", "distance_matrix_correlations.tsv"))

  diag <- o2ipa_pairwise_diagnostics(d_parameter, d_static, d_static18, d_dynamic, inputs$manifest)
  if (is.list(diag)) {
    o2ipa_write_tsv(diag$pairs, file.path(out_dir, "tables", "pairwise_regime_diagnostics.tsv"))
    o2ipa_write_tsv(diag$thresholds, file.path(out_dir, "tables", "distance_thresholds.tsv"))
  } else {
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", "pairwise_regime_diagnostics.tsv"))
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", "distance_thresholds.tsv"))
  }

  cl_static <- o2ipa_cluster_distance(d_static, n_bootstrap = n_bootstrap, random_seed = random_seed, feature_wide = static_full_scaled$wide, feature_meta = static_full_scaled$metadata)
  cl_static18 <- o2ipa_cluster_distance(d_static18, n_bootstrap = n_bootstrap, random_seed = random_seed, feature_wide = static_18_scaled$wide, feature_meta = static_18_scaled$metadata)
  cl_dynamic <- o2ipa_cluster_distance(d_dynamic, n_bootstrap = n_bootstrap, random_seed = random_seed, feature_wide = dynamic_scaled$wide, feature_meta = dynamic_scaled$metadata)
  diag_all <- o2ipa_bind_rows_fill(list(
    o2ipa_add_cluster_source(cl_static$diagnostics, "static_full"),
    o2ipa_add_cluster_source(cl_static18$diagnostics, "static_18only"),
    o2ipa_add_cluster_source(cl_dynamic$diagnostics, "dynamic")
  ))
  o2ipa_write_tsv(diag_all, file.path(out_dir, "tables", "cluster_k_diagnostics.tsv"))
  o2ipa_write_tsv(cl_static$membership, file.path(out_dir, "tables", "cluster_membership_static_full.tsv"))
  o2ipa_write_tsv(cl_static18$membership, file.path(out_dir, "tables", "cluster_membership_static_18only.tsv"))
  o2ipa_write_tsv(cl_dynamic$membership, file.path(out_dir, "tables", "cluster_membership_dynamic.tsv"))
  o2ipa_write_tsv(cl_static$stability, file.path(out_dir, "tables", "cluster_stability.tsv"))
  medoid_all <- o2ipa_bind_rows_fill(list(
    o2ipa_add_cluster_source(cl_static$medoids, "static_full"),
    o2ipa_add_cluster_source(cl_static18$medoids, "static_18only"),
    o2ipa_add_cluster_source(cl_dynamic$medoids, "dynamic")
  ))
  o2ipa_write_tsv(medoid_all, file.path(out_dir, "tables", "cluster_medoids.tsv"))
  if (!is.null(cl_static$consensus)) o2ipa_write_tsv(data.frame(seed_id = rownames(cl_static$consensus), cl_static$consensus, check.names = FALSE), file.path(out_dir, "tables", "cluster_consensus_matrix.tsv")) else o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", "cluster_consensus_matrix.tsv"))

  conc18 <- o2ipa_concordance(cl_static$membership, cl_static18$membership, cl_static$recommended_k, cl_static18$recommended_k, "static_full_vs_static_18only")
  concdyn <- o2ipa_concordance(cl_static$membership, cl_dynamic$membership, cl_static$recommended_k, cl_dynamic$recommended_k, "static_vs_dynamic")
  o2ipa_write_tsv(conc18, file.path(out_dir, "tables", "cluster_concordance_full_vs_18only.tsv"))
  o2ipa_write_tsv(concdyn, file.path(out_dir, "tables", "cluster_concordance_static_vs_dynamic.tsv"))

  o2ipa_cluster_qc(inputs$manifest, cl_static, out_dir, "static_full")
  lomo <- o2ipa_leave_one_module_out(static_full_scaled$long, cl_static$recommended_k, cl_static$membership, out_dir)
  module_centroids <- o2ipa_module_centroids(static_full_scaled$long, cl_static$membership, cl_static$recommended_k, "static_full")
  o2ipa_write_tsv(module_centroids, file.path(out_dir, "tables", "cluster_process_module_centroid_profile.tsv"))

  o2ipa_plot_outputs(out_dir, d_parameter, d_static, d_static18, d_dynamic, static_full_scaled$wide, dynamic_scaled$wide, cl_static, inputs$manifest, static$full_long)

  limitations <- c(
    "Static process probes use current repository model helpers and C++ buffering backend where available.",
    "Dynamic phenotype features are limited to existing prediction tables; unavailable O2 and entropy metrics are listed explicitly.",
    paste0("The main objective delta argument was recorded as ", main_objective_delta, "; seeds were not silently dropped solely for boundary risk.")
  )
  o2ipa_report(out_dir, run_dir, inputs$manifest, param$metadata, feature_meta, dist_cor, cl_static, dynamic$unavailable, limitations)
  message("Completed analysis: ", out_dir)
}

if (sys.nframe() == 0L) {
  main()
}
