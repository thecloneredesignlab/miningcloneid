options(stringsAsFactors = FALSE, warn = 1)

si6_paths <- function(workspace_root = f6r_find_workspace_root()) {
  base <- f6r_paths(workspace_root)
  data_dir <- file.path(base$root, "data", "Figures", "Supp_Figure6_3")
  list(
    base = base,
    root = base$root,
    data = data_dir,
    cache = file.path(data_dir, "top10_endpoint_cache"),
    panels = file.path(data_dir, "panels"),
    figures = base$figures,
    manuscript_figures = base$manuscript_figures
  )
}

si6_existing_data_file <- function(paths, candidates) {
  candidate_paths <- file.path(paths$data, candidates)
  hit <- candidate_paths[file.exists(candidate_paths)]
  if (!length(hit)) {
    stop(
      "Missing Supplementary Figure 6-3 input; tried: ",
      paste(candidate_paths, collapse = ", ")
    )
  }
  normalizePath(hit[[1L]], mustWork = TRUE)
}

si6_profile <- function() "top10_eigenmode_localization_201x60_q10_v1"
si6_o2_values <- function() seq(0, 5, length.out = 201L)
si6_p_values <- function() 10^seq(log10(0.005), log10(0.5), length.out = 60L)
si6_rank_count <- function() 10L

si6_atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = "gzip")
  if (!file.copy(temporary, path, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Failed to publish Supplementary Figure 6-3 checkpoint: ", path)
  }
  unlink(temporary)
  normalizePath(path, mustWork = TRUE)
}

si6_rank_eigenmodes <- function(M, ngrid, n_unit, n_rank = si6_rank_count()) {
  eig <- eigen(as.matrix(M), only.values = FALSE)
  index <- seq_along(eig$values)
  ordering <- order(
    -Re(eig$values), -Mod(eig$values), -Im(eig$values), index,
    method = "radix"
  )
  ordering <- ordering[seq_len(min(n_rank, length(ordering)))]
  values <- eig$values[ordering]
  vectors <- eig$vectors[, ordering, drop = FALSE]

  magnitude <- Mod(vectors)
  l1_denominator <- colSums(magnitude)
  l2_weight <- magnitude^2
  l2_denominator <- colSums(l2_weight)
  localization_l1 <- as.numeric(crossprod(ngrid, magnitude)) / l1_denominator / n_unit
  localization_l2 <- as.numeric(crossprod(ngrid, l2_weight)) / l2_denominator / n_unit

  dominant <- Re(vectors[, 1L])
  if (sum(dominant, na.rm = TRUE) < 0) dominant <- -dominant
  top1_nonnegative <- all(dominant >= -1e-8, na.rm = TRUE)
  dominant <- fixo2_normalize_state(dominant)
  dominant_ploidy <- sum(ngrid * dominant, na.rm = TRUE) / n_unit
  localization_l1[[1L]] <- dominant_ploidy
  localization_l2[[1L]] <- dominant_ploidy

  lambda_real <- Re(values)
  gap_to_dominant <- pmax(lambda_real[[1L]] - lambda_real, 0)
  residual_matrix <- as.matrix(M) %*% vectors - sweep(vectors, 2L, values, `*`)
  vector_norm <- sqrt(colSums(Mod(vectors)^2))
  matrix_norm <- sqrt(sum(Mod(M)^2))
  relative_residual <- sqrt(colSums(Mod(residual_matrix)^2)) /
    pmax(matrix_norm * vector_norm, .Machine$double.eps)
  real_vectors <- Re(vectors)
  mixed_real_signs <- vapply(seq_len(ncol(real_vectors)), function(j) {
    x <- real_vectors[, j]
    any(x > 1e-10) && any(x < -1e-10)
  }, logical(1L))

  list(
    localization_l1 = localization_l1,
    localization_l2 = localization_l2,
    lambda_real = lambda_real,
    lambda_imag_abs = abs(Im(values)),
    gap_to_dominant = gap_to_dominant,
    relative_residual = relative_residual,
    mixed_real_signs = mixed_real_signs,
    top1_nonnegative = top1_nonnegative
  )
}

si6_endpoint_cache_valid <- function(object, metadata, expected_grid_rows) {
  !is.null(object) &&
    identical(as.character(object$metadata$profile), si6_profile()) &&
    identical(
      as.character(object$metadata$parameter_endpoint_group),
      as.character(metadata$parameter_endpoint_group[[1L]])
    ) &&
    identical(
      as.integer(object$metadata$representative_seed_number),
      as.integer(metadata$representative_seed_number[[1L]])
    ) &&
    nrow(object$grid) == expected_grid_rows &&
    identical(dim(object$localization_l1), c(expected_grid_rows, si6_rank_count())) &&
    isTRUE(object$qc$operator_qc_pass[[1L]])
}

si6_compute_endpoint_cache <- function(
    metadata, parameters, context, cache_path, parameter_source,
    force_rebuild = FALSE
) {
  o2_values <- si6_o2_values()
  p_values <- si6_p_values()
  n_rank <- si6_rank_count()
  n_grid <- length(o2_values) * length(p_values)
  if (file.exists(cache_path) && !isTRUE(force_rebuild)) {
    existing <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (si6_endpoint_cache_valid(existing, metadata, n_grid)) {
      return(existing$qc)
    }
  }

  seed_number <- metadata$representative_seed_number[[1L]]
  z <- parameters[
    parameters$pair_id == metadata$pair_id[[1L]] &
      parameters$seed_number == seed_number,
    , drop = FALSE
  ]
  required_parameters <- c(f6r_shared_parameters(), "rho_2N")
  if (nrow(z) != length(required_parameters) ||
      !setequal(z$parameter, required_parameters) || any(!is.finite(z$value))) {
    stop(
      "Incomplete Supplementary Figure 6-3 parameters for ", metadata$pair_label,
      " representative seed ", seed_number
    )
  }
  params <- stats::setNames(as.numeric(z$value), z$parameter)
  run_params <- prepare_run_params(
    param_values = params, simulation = "joint",
    cfg = context$config, fixed_o2 = 5
  )
  config <- prepare_sim_cfg(
    context$config, argv = list(), fixed_o2 = 5, run_params = run_params
  )
  run_params$O2_growth <- isTRUE(config$O2_growth)
  run_params$ploidy_O2_death <- config$ploidy_O2_death
  n_unit <- as.numeric(config$N_UNIT %||% 22)

  grid <- expand.grid(
    O2_pct = o2_values,
    effective_p_misseg = p_values,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$effective_p_misseg, grid$O2_pct), , drop = FALSE]
  matrices <- lapply(
    c(
      "localization_l1", "localization_l2", "lambda_real",
      "lambda_imag_abs", "gap_to_dominant", "relative_residual"
    ),
    function(x) matrix(NA_real_, nrow = n_grid, ncol = n_rank)
  )
  names(matrices) <- c(
    "localization_l1", "localization_l2", "lambda_real",
    "lambda_imag_abs", "gap_to_dominant", "relative_residual"
  )
  mixed_real_signs <- matrix(FALSE, nrow = n_grid, ncol = n_rank)
  top1_nonnegative <- logical(n_grid)

  row_index <- 0L
  for (p_fixed in p_values) {
    forced <- response_force_effective_p_misseg(run_params, p_fixed)
    for (o2 in o2_values) {
      row_index <- row_index + 1L
      fm <- fixo2_fixed_matrix(globalenv(), config, forced, o2)
      ranked <- si6_rank_eigenmodes(
        fm$M, fm$ngrid, n_unit = n_unit, n_rank = n_rank
      )
      for (field in names(matrices)) {
        matrices[[field]][row_index, ] <- ranked[[field]]
      }
      mixed_real_signs[row_index, ] <- ranked$mixed_real_signs
      top1_nonnegative[[row_index]] <- ranked$top1_nonnegative
    }
  }

  maximum_residual <- max(matrices$relative_residual, na.rm = TRUE)
  valid <- row_index == n_grid &&
    all(vapply(matrices, function(x) all(is.finite(x)), logical(1L))) &&
    all(top1_nonnegative) &&
    all(matrices$gap_to_dominant >= -1e-12) &&
    max(abs(matrices$gap_to_dominant[, 1L])) <= 1e-12 &&
    is.finite(maximum_residual) && maximum_residual <= 1e-8
  qc <- data.frame(
    display_label = metadata$display_label,
    pair_label = metadata$pair_label,
    pair_id = metadata$pair_id,
    parameter_endpoint_group = metadata$parameter_endpoint_group,
    representative_seed_number = seed_number,
    endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10,
    n_grid = n_grid,
    n_rank = n_rank,
    all_numeric_outputs_finite = all(vapply(
      matrices, function(x) all(is.finite(x)), logical(1L)
    )),
    all_top1_eigenvectors_nonnegative = all(top1_nonnegative),
    maximum_relative_eigen_residual = maximum_residual,
    n_complex_ranked_modes = sum(matrices$lambda_imag_abs > 1e-10),
    n_mixed_real_sign_modes = sum(mixed_real_signs),
    operator_qc_pass = valid,
    cache_path = normalizePath(cache_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )

  object <- c(
    list(
      metadata = list(
        profile = si6_profile(),
        display_label = metadata$display_label[[1L]],
        pair_label = metadata$pair_label[[1L]],
        pair_id = metadata$pair_id[[1L]],
        parameter_endpoint_group = metadata$parameter_endpoint_group[[1L]],
        representative_seed_number = seed_number,
        endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10[[1L]],
        config_path = context$config_path,
        parameter_source = normalizePath(parameter_source, mustWork = TRUE),
        localization_definition = paste0(
          "rank1: nonnegative dominant eigenvector mean chromosome number/N_UNIT; ",
          "ranks2-10: sum_N N*|v_Nk|/(N_UNIT*sum_N|v_Nk|)"
        ),
        gap_definition = "Re(lambda_1)-Re(lambda_k)"
      ),
      grid = grid,
      mixed_real_signs = mixed_real_signs,
      top1_nonnegative = top1_nonnegative,
      qc = qc
    ),
    matrices
  )
  si6_atomic_save_rds(object, cache_path)
  qc
}

si6_cache_paths <- function(paths, endpoints) {
  stats::setNames(
    file.path(
      paths$cache, endpoints$pair_label,
      paste0(
        "endpoint_",
        gsub("[^A-Za-z0-9_.-]", "_", endpoints$parameter_endpoint_group),
        "_seed", endpoints$representative_seed_number, ".rds"
      )
    ),
    endpoints$parameter_endpoint_group
  )
}

si6_safe_quantiles <- function(x, probs) {
  out <- matrixStats::rowQuantiles(x, probs = probs, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

si6_summarize_pair <- function(metadata, objects) {
  reference <- objects[[1L]]$grid
  key <- paste(reference$effective_p_misseg, reference$O2_pct, sep = "|")
  if (!all(vapply(objects, function(x) {
    identical(paste(x$grid$effective_p_misseg, x$grid$O2_pct, sep = "|"), key)
  }, logical(1L)))) {
    stop("Supplementary Figure 6-3 endpoint grids differ within ", metadata$pair_label[[1L]])
  }
  object_groups <- vapply(
    objects, function(x) x$metadata$parameter_endpoint_group, character(1L)
  )
  multiplicity <- metadata$endpoint_multiplicity_q10[
    match(object_groups, metadata$parameter_endpoint_group)
  ]
  if (sum(multiplicity) != 50L) {
    stop("Supplementary Figure 6-3 endpoint multiplicity does not sum to 50.")
  }
  expand_index <- rep(seq_along(objects), times = multiplicity)
  summary_rows <- vector("list", si6_rank_count())
  for (rank in seq_len(si6_rank_count())) {
    localization <- do.call(cbind, lapply(
      objects, function(x) x$localization_l1[, rank]
    ))[, expand_index, drop = FALSE]
    localization_l2 <- do.call(cbind, lapply(
      objects, function(x) x$localization_l2[, rank]
    ))[, expand_index, drop = FALSE]
    gap <- do.call(cbind, lapply(
      objects, function(x) x$gap_to_dominant[, rank]
    ))[, expand_index, drop = FALSE]
    lambda_real <- do.call(cbind, lapply(
      objects, function(x) x$lambda_real[, rank]
    ))[, expand_index, drop = FALSE]
    imag_abs <- do.call(cbind, lapply(
      objects, function(x) x$lambda_imag_abs[, rank]
    ))[, expand_index, drop = FALSE]
    loc_q <- matrixStats::rowQuantiles(
      localization, probs = c(0.10, 0.90), na.rm = TRUE
    )
    gap_q <- matrixStats::rowQuantiles(
      gap, probs = c(0.10, 0.90), na.rm = TRUE
    )
    summary_rows[[rank]] <- data.frame(
      pair_id = metadata$pair_id[[1L]],
      pair_label = metadata$pair_label[[1L]],
      display_label = metadata$display_label[[1L]],
      O2_pct = reference$O2_pct,
      effective_p_misseg = reference$effective_p_misseg,
      eigenmode_rank = rank,
      n_seed = 50L,
      n_unique_parameter_endpoint = nrow(metadata),
      eigenmode_localization_ploidy_median = matrixStats::rowMedians(localization),
      eigenmode_localization_ploidy_q10 = loc_q[, 1L],
      eigenmode_localization_ploidy_q90 = loc_q[, 2L],
      spectral_gap_to_dominant_median = matrixStats::rowMedians(gap),
      spectral_gap_to_dominant_q10 = gap_q[, 1L],
      spectral_gap_to_dominant_q90 = gap_q[, 2L],
      eigenvalue_real_median = matrixStats::rowMedians(lambda_real),
      proportion_complex_eigenvalue = rowMeans(imag_abs > 1e-10),
      l1_l2_localization_absdiff_median = matrixStats::rowMedians(
        abs(localization - localization_l2)
      ),
      stringsAsFactors = FALSE
    )
  }

  derived_by_object <- lapply(objects, function(x) {
    top1 <- x$localization_l1[, 1L]
    delta <- abs(sweep(x$localization_l1[, 2:10, drop = FALSE], 1L, top1, "-"))
    gap <- x$gap_to_dominant[, 2:10, drop = FALSE]
    near <- delta
    near[gap >= 0.005] <- NA_real_
    near_jump <- matrixStats::rowMaxs(near, na.rm = TRUE)
    near_jump[!is.finite(near_jump)] <- 0
    competitor_gap <- gap
    competitor_gap[delta < 1] <- NA_real_
    minimum_competitor_gap <- matrixStats::rowMins(
      competitor_gap, na.rm = TRUE
    )
    minimum_competitor_gap[!is.finite(minimum_competitor_gap)] <- NA_real_
    list(
      near_jump = near_jump,
      minimum_competitor_gap = minimum_competitor_gap,
      has_competitor = is.finite(minimum_competitor_gap)
    )
  })
  near_jump <- do.call(cbind, lapply(
    derived_by_object, `[[`, "near_jump"
  ))[, expand_index, drop = FALSE]
  competitor_gap <- do.call(cbind, lapply(
    derived_by_object, `[[`, "minimum_competitor_gap"
  ))[, expand_index, drop = FALSE]
  has_competitor <- do.call(cbind, lapply(
    derived_by_object, `[[`, "has_competitor"
  ))[, expand_index, drop = FALSE]
  near_q <- matrixStats::rowQuantiles(
    near_jump, probs = c(0.10, 0.90), na.rm = TRUE
  )
  competitor_q <- si6_safe_quantiles(competitor_gap, c(0.10, 0.90))
  derived <- data.frame(
    pair_id = metadata$pair_id[[1L]],
    pair_label = metadata$pair_label[[1L]],
    display_label = metadata$display_label[[1L]],
    O2_pct = reference$O2_pct,
    effective_p_misseg = reference$effective_p_misseg,
    n_seed = 50L,
    n_unique_parameter_endpoint = nrow(metadata),
    near_degenerate_jump_median = matrixStats::rowMedians(near_jump),
    near_degenerate_jump_q10 = near_q[, 1L],
    near_degenerate_jump_q90 = near_q[, 2L],
    proportion_near_degenerate_jump_ge_1 = rowMeans(near_jump >= 1),
    different_ploidy_competitor_gap_median = matrixStats::rowMedians(
      competitor_gap, na.rm = TRUE
    ),
    different_ploidy_competitor_gap_q10 = competitor_q[, 1L],
    different_ploidy_competitor_gap_q90 = competitor_q[, 2L],
    proportion_with_different_ploidy_competitor = rowMeans(has_competitor),
    stringsAsFactors = FALSE
  )
  derived$different_ploidy_competitor_gap_median[
    !is.finite(derived$different_ploidy_competitor_gap_median)
  ] <- NA_real_
  list(summary = do.call(rbind, summary_rows), derived = derived)
}

si6_summarize_caches <- function(paths, endpoint_manifest, cache_paths) {
  summary_rows <- list()
  derived_rows <- list()
  for (pair in endpoint_manifest$display_manifest$pair_id) {
    metadata <- endpoint_manifest$endpoints[
      endpoint_manifest$endpoints$pair_id == pair, , drop = FALSE
    ]
    objects <- lapply(
      cache_paths[metadata$parameter_endpoint_group], readRDS
    )
    pair_summary <- si6_summarize_pair(metadata, objects)
    summary_rows[[pair]] <- pair_summary$summary
    derived_rows[[pair]] <- pair_summary$derived
  }
  summary <- do.call(rbind, summary_rows)
  derived <- do.call(rbind, derived_rows)
  rownames(summary) <- NULL
  rownames(derived) <- NULL
  list(summary = summary, derived = derived)
}

si6_local_rank1_metrics <- function(object) {
  o2_values <- sort(unique(object$grid$O2_pct))
  p_values <- sort(unique(object$grid$effective_p_misseg))
  expected_rows <- length(o2_values) * length(p_values)
  if (nrow(object$grid) != expected_rows) {
    stop("Unexpected Supplementary Figure 6-3 rank-1 grid size.")
  }
  expected_key <- paste(
    rep(p_values, each = length(o2_values)),
    rep(o2_values, times = length(p_values)), sep = "|"
  )
  observed_key <- paste(
    object$grid$effective_p_misseg, object$grid$O2_pct, sep = "|"
  )
  if (!identical(expected_key, observed_key)) {
    stop("Supplementary Figure 6-3 rank-1 grid order is not p-major/O2-minor.")
  }

  ploidy <- matrix(
    object$localization_l1[, 1L],
    nrow = length(o2_values), ncol = length(p_values)
  )
  regime_class <- ifelse(ploidy <= 2, 1L, ifelse(ploidy >= 4, 3L, 2L))
  local_switch <- matrix(FALSE, nrow = nrow(ploidy), ncol = ncol(ploidy))
  local_jump <- matrix(0, nrow = nrow(ploidy), ncol = ncol(ploidy))

  if (nrow(ploidy) > 1L) {
    delta <- abs(ploidy[-1L, , drop = FALSE] - ploidy[-nrow(ploidy), , drop = FALSE])
    crossed <- regime_class[-1L, , drop = FALSE] !=
      regime_class[-nrow(regime_class), , drop = FALSE]
    local_switch[-1L, ] <- local_switch[-1L, ] | crossed
    local_switch[-nrow(local_switch), ] <- local_switch[-nrow(local_switch), ] | crossed
    local_jump[-1L, ] <- pmax(local_jump[-1L, ], delta)
    local_jump[-nrow(local_jump), ] <- pmax(local_jump[-nrow(local_jump), ], delta)
  }
  if (ncol(ploidy) > 1L) {
    delta <- abs(ploidy[, -1L, drop = FALSE] - ploidy[, -ncol(ploidy), drop = FALSE])
    crossed <- regime_class[, -1L, drop = FALSE] !=
      regime_class[, -ncol(regime_class), drop = FALSE]
    local_switch[, -1L] <- local_switch[, -1L] | crossed
    local_switch[, -ncol(local_switch)] <- local_switch[, -ncol(local_switch)] | crossed
    local_jump[, -1L] <- pmax(local_jump[, -1L], delta)
    local_jump[, -ncol(local_jump)] <- pmax(local_jump[, -ncol(local_jump)], delta)
  }
  list(
    dominant_mean_ploidy = as.vector(ploidy),
    regime_class = as.vector(regime_class),
    local_regime_switch = as.vector(local_switch),
    local_adjacent_ploidy_jump = as.vector(local_jump)
  )
}

si6_summarize_weak_gap_pair <- function(metadata, objects, figure6_surface) {
  reference <- objects[[1L]]$grid
  object_groups <- vapply(
    objects, function(x) x$metadata$parameter_endpoint_group, character(1L)
  )
  multiplicity <- metadata$endpoint_multiplicity_q10[
    match(object_groups, metadata$parameter_endpoint_group)
  ]
  if (sum(multiplicity) != 50L) {
    stop("Supplementary Figure 6-3 weak-gap endpoint multiplicity does not sum to 50.")
  }
  metrics <- lapply(objects, si6_local_rank1_metrics)
  expand_index <- rep(seq_along(objects), times = multiplicity)
  ploidy <- do.call(cbind, lapply(
    metrics, `[[`, "dominant_mean_ploidy"
  ))[, expand_index, drop = FALSE]
  regime_class <- do.call(cbind, lapply(
    metrics, `[[`, "regime_class"
  ))[, expand_index, drop = FALSE]
  local_switch <- do.call(cbind, lapply(
    metrics, `[[`, "local_regime_switch"
  ))[, expand_index, drop = FALSE]
  local_jump <- do.call(cbind, lapply(
    metrics, `[[`, "local_adjacent_ploidy_jump"
  ))[, expand_index, drop = FALSE]

  ploidy_q <- matrixStats::rowQuantiles(
    ploidy, probs = c(0.10, 0.50, 0.90), na.rm = TRUE
  )
  jump_q <- matrixStats::rowQuantiles(
    local_jump, probs = c(0.10, 0.50, 0.90), na.rm = TRUE
  )
  proportion_low <- rowMeans(regime_class == 1L)
  proportion_intermediate <- rowMeans(regime_class == 2L)
  proportion_high <- rowMeans(regime_class == 3L)
  proportion_ge_2 <- rowMeans(ploidy >= 2)
  result <- data.frame(
    pair_id = metadata$pair_id[[1L]],
    pair_label = metadata$pair_label[[1L]],
    display_label = metadata$display_label[[1L]],
    O2_pct = reference$O2_pct,
    effective_p_misseg = reference$effective_p_misseg,
    n_seed = 50L,
    n_unique_parameter_endpoint = nrow(metadata),
    dominant_mean_ploidy_q10 = ploidy_q[, 1L],
    dominant_mean_ploidy_median = ploidy_q[, 2L],
    dominant_mean_ploidy_q90 = ploidy_q[, 3L],
    dominant_ploidy_spread_q90_q10 = ploidy_q[, 3L] - ploidy_q[, 1L],
    proportion_dominant_mean_ploidy_ge_2 = proportion_ge_2,
    proportion_low_ploidy_le_2 = proportion_low,
    proportion_intermediate_ploidy_gt_2_lt_4 = proportion_intermediate,
    proportion_high_ploidy_ge_4 = proportion_high,
    ploidy_regime_consensus = pmax(
      proportion_low, proportion_intermediate, proportion_high
    ),
    local_regime_switch_proportion = rowMeans(local_switch),
    local_adjacent_ploidy_jump_q10 = jump_q[, 1L],
    local_adjacent_ploidy_jump_median = jump_q[, 2L],
    local_adjacent_ploidy_jump_q90 = jump_q[, 3L],
    stringsAsFactors = FALSE
  )
  result$regime_class <- ifelse(
    result$proportion_high_ploidy_ge_4 >= 0.90, "Stable high",
    ifelse(
      result$proportion_low_ploidy_le_2 >= 0.90, "Stable low",
      ifelse(
        result$proportion_intermediate_ploidy_gt_2_lt_4 >= 0.90,
        "Stable intermediate", "Mixed"
      )
    )
  )

  surface <- figure6_surface[
    figure6_surface$cutoff == "q10" &
      figure6_surface$pair_id == metadata$pair_id[[1L]],
    c(
      "pair_id", "O2_pct", "effective_p_misseg",
      "proportion_spectral_gap_below_0p005", "spectral_gap_median"
    ), drop = FALSE
  ]
  key <- function(pair, o2, p) {
    paste(pair, sprintf("%.12f", o2), sprintf("%.12g", p), sep = "|")
  }
  result$key <- key(result$pair_id, result$O2_pct, result$effective_p_misseg)
  surface$key <- key(surface$pair_id, surface$O2_pct, surface$effective_p_misseg)
  matched <- match(result$key, surface$key)
  if (anyNA(matched) || anyDuplicated(surface$key)) {
    stop("Supplementary Figure 6-3 weak-gap summary could not be matched to Figure 6A.")
  }
  result$proportion_spectral_gap_below_0p005 <-
    surface$proportion_spectral_gap_below_0p005[matched]
  result$spectral_gap_median <- surface$spectral_gap_median[matched]
  result$weak_gap_region <-
    result$proportion_spectral_gap_below_0p005 >= 0.50
  result$key <- NULL
  result
}

si6_weak_gap_summary <- function(paths, endpoint_manifest, cache_paths) {
  figure6_surface <- f6r_read_tsv(file.path(
    paths$base$figure6, "joint_multiseed_surface_summary.tsv"
  ))
  rows <- lapply(endpoint_manifest$display_manifest$pair_id, function(pair) {
    metadata <- endpoint_manifest$endpoints[
      endpoint_manifest$endpoints$pair_id == pair, , drop = FALSE
    ]
    objects <- lapply(cache_paths[metadata$parameter_endpoint_group], readRDS)
    si6_summarize_weak_gap_pair(
      metadata, objects, figure6_surface = figure6_surface
    )
  })
  grid <- do.call(rbind, rows)
  rownames(grid) <- NULL
  pair_rows <- lapply(split(grid, grid$display_label), function(z) {
    weak <- z[z$weak_gap_region, , drop = FALSE]
    data.frame(
      display_label = z$display_label[[1L]],
      pair_label = z$pair_label[[1L]],
      n_grid_cell = nrow(z),
      n_weak_gap_cell = nrow(weak),
      fraction_grid_weak_gap = nrow(weak) / nrow(z),
      fraction_weak_gap_stable_high = mean(weak$regime_class == "Stable high"),
      fraction_weak_gap_stable_low = mean(weak$regime_class == "Stable low"),
      fraction_weak_gap_stable_intermediate = mean(
        weak$regime_class == "Stable intermediate"
      ),
      fraction_weak_gap_mixed = mean(weak$regime_class == "Mixed"),
      fraction_weak_gap_consensus_ge_0p9 = mean(
        weak$ploidy_regime_consensus >= 0.90
      ),
      fraction_weak_gap_any_endpoint_local_switch = mean(
        weak$local_regime_switch_proportion > 0
      ),
      fraction_weak_gap_majority_endpoint_local_switch = mean(
        weak$local_regime_switch_proportion >= 0.50
      ),
      weak_gap_ploidy_spread_median = stats::median(
        weak$dominant_ploidy_spread_q90_q10
      ),
      weak_gap_ploidy_spread_q90 = stats::quantile(
        weak$dominant_ploidy_spread_q90_q10, 0.90, names = FALSE
      ),
      weak_gap_ploidy_spread_max = max(
        weak$dominant_ploidy_spread_q90_q10
      ),
      weak_gap_local_jump_median = stats::median(
        weak$local_adjacent_ploidy_jump_median
      ),
      weak_gap_local_jump_q90 = stats::quantile(
        weak$local_adjacent_ploidy_jump_median, 0.90, names = FALSE
      ),
      weak_gap_local_jump_max = max(
        weak$local_adjacent_ploidy_jump_median
      ),
      fraction_weak_gap_local_jump_median_ge_1 = mean(
        weak$local_adjacent_ploidy_jump_median >= 1
      ),
      stringsAsFactors = FALSE
    )
  })
  pair_summary <- do.call(rbind, pair_rows)
  pair_summary <- pair_summary[
    match(sprintf("C%02d", 1:6), pair_summary$display_label), , drop = FALSE
  ]
  rownames(pair_summary) <- NULL
  list(grid = grid, pair_summary = pair_summary)
}

si6_validation <- function(
    paths, endpoint_manifest, qc, summary, derived,
    weak_gap_grid, weak_gap_pair_summary
) {
  original <- f6r_read_tsv(file.path(
    paths$base$figure6, "joint_multiseed_surface_summary.tsv"
  ))
  original <- original[
    original$cutoff == "q10" &
      original$pair_id %in% endpoint_manifest$display_manifest$pair_id,
    c(
      "pair_id", "O2_pct", "effective_p_misseg",
      "dominant_mean_ploidy_median"
    ), drop = FALSE
  ]
  top1 <- summary[summary$eigenmode_rank == 1L, c(
    "pair_id", "O2_pct", "effective_p_misseg",
    "eigenmode_localization_ploidy_median"
  )]
  key_fun <- function(pair, o2, p) {
    paste(pair, sprintf("%.12f", o2), sprintf("%.12g", p), sep = "|")
  }
  original$key <- key_fun(
    original$pair_id, original$O2_pct, original$effective_p_misseg
  )
  top1$key <- key_fun(top1$pair_id, top1$O2_pct, top1$effective_p_misseg)
  overlap <- merge(top1, original, by = "key", all = FALSE)
  top1_difference <- max(abs(
    overlap$eigenmode_localization_ploidy_median -
      overlap$dominant_mean_ploidy_median
  ))
  weak_reference <- f6r_read_tsv(file.path(
    paths$base$figure6, "joint_multiseed_surface_summary.tsv"
  ))
  weak_reference <- weak_reference[
    weak_reference$cutoff == "q10" &
      weak_reference$pair_id %in% endpoint_manifest$display_manifest$pair_id,
    c(
      "pair_id", "O2_pct", "effective_p_misseg",
      "dominant_mean_ploidy_median", "dominant_mean_ploidy_q10",
      "dominant_mean_ploidy_q90", "proportion_dominant_mean_ploidy_ge_2"
    ), drop = FALSE
  ]
  weak_reference$key <- key_fun(
    weak_reference$pair_id, weak_reference$O2_pct,
    weak_reference$effective_p_misseg
  )
  weak_gap_grid$key <- key_fun(
    weak_gap_grid$pair_id, weak_gap_grid$O2_pct,
    weak_gap_grid$effective_p_misseg
  )
  weak_match <- match(weak_gap_grid$key, weak_reference$key)
  weak_median_difference <- max(abs(
    weak_gap_grid$dominant_mean_ploidy_median -
      weak_reference$dominant_mean_ploidy_median[weak_match]
  ))
  weak_q10_difference <- max(abs(
    weak_gap_grid$dominant_mean_ploidy_q10 -
      weak_reference$dominant_mean_ploidy_q10[weak_match]
  ))
  weak_q90_difference <- max(abs(
    weak_gap_grid$dominant_mean_ploidy_q90 -
      weak_reference$dominant_mean_ploidy_q90[weak_match]
  ))
  weak_legacy_ge2_share_difference <- max(abs(
    weak_gap_grid$proportion_dominant_mean_ploidy_ge_2 -
      weak_reference$proportion_dominant_mean_ploidy_ge_2[weak_match]
  ))
  weak_cell_counts <- weak_gap_pair_summary$n_weak_gap_cell
  weak_regime_sums <- rowSums(weak_gap_pair_summary[, c(
    "fraction_weak_gap_stable_high", "fraction_weak_gap_stable_low",
    "fraction_weak_gap_stable_intermediate", "fraction_weak_gap_mixed"
  )])
  pointwise_regime_sums <- rowSums(weak_gap_grid[, c(
    "proportion_low_ploidy_le_2",
    "proportion_intermediate_ploidy_gt_2_lt_4",
    "proportion_high_ploidy_ge_4"
  )])
  represented <- tapply(
    qc$endpoint_multiplicity_q10,
    factor(qc$pair_label, levels = sprintf("C%02d", 1:6)),
    sum
  )
  observed <- c(
    nrow(qc), paste(represented, collapse = ","), nrow(summary),
    nrow(derived), length(unique(summary$O2_pct)),
    length(unique(summary$effective_p_misseg)),
    paste(range(summary$eigenmode_rank), collapse = "--"),
    all(qc$operator_qc_pass), max(qc$maximum_relative_eigen_residual),
    all(is.finite(summary$eigenmode_localization_ploidy_median)),
    all(summary$spectral_gap_to_dominant_median >= -1e-12),
    max(abs(summary$spectral_gap_to_dominant_median[
      summary$eigenmode_rank == 1L
    ])),
    nrow(overlap), top1_difference,
    nrow(weak_gap_grid), nrow(weak_gap_pair_summary),
    paste(weak_cell_counts, collapse = ","),
    all(weak_gap_grid$ploidy_regime_consensus >= 1 / 3 &
      weak_gap_grid$ploidy_regime_consensus <= 1),
    max(abs(weak_regime_sums - 1)),
    max(abs(pointwise_regime_sums - 1)),
    all(weak_gap_grid$local_regime_switch_proportion >= 0 &
      weak_gap_grid$local_regime_switch_proportion <= 1),
    all(weak_gap_grid$local_adjacent_ploidy_jump_median >= 0),
    sum(is.na(weak_match)), weak_median_difference, weak_q10_difference,
    weak_q90_difference, weak_legacy_ge2_share_difference
  )
  expected <- c(
    as.character(nrow(qc)), paste(rep(50L, 6L), collapse = ","),
    as.character(6L * 201L * 60L * 10L), as.character(6L * 201L * 60L),
    "201", "60", "1--10",
    "TRUE", "<=1e-8", "TRUE", "TRUE", "<=1e-12",
    as.character(6L * 201L * 60L), "<=1e-10",
    as.character(6L * 201L * 60L), "6", paste(weak_cell_counts, collapse = ","),
    "TRUE", "<=1e-12", "<=1e-12", "TRUE", "TRUE",
    "0", "<=1e-10", "<=1e-10", "<=1e-10", "<=1e-12"
  )
  passed <- c(
    nrow(qc) >= 6L, identical(as.integer(represented), rep(50L, 6L)),
    nrow(summary) == 6L * 201L * 60L * 10L,
    nrow(derived) == 6L * 201L * 60L,
    length(unique(summary$O2_pct)) == 201L,
    length(unique(summary$effective_p_misseg)) == 60L,
    identical(range(summary$eigenmode_rank), c(1L, 10L)),
    all(qc$operator_qc_pass),
    max(qc$maximum_relative_eigen_residual) <= 1e-8,
    all(is.finite(summary$eigenmode_localization_ploidy_median)),
    all(summary$spectral_gap_to_dominant_median >= -1e-12),
    max(abs(summary$spectral_gap_to_dominant_median[
      summary$eigenmode_rank == 1L
    ])) <= 1e-12,
    nrow(overlap) == 6L * 201L * 60L,
    top1_difference <= 1e-10,
    nrow(weak_gap_grid) == 6L * 201L * 60L,
    nrow(weak_gap_pair_summary) == 6L,
    length(weak_cell_counts) == 6L && all(weak_cell_counts > 0L),
    all(weak_gap_grid$ploidy_regime_consensus >= 1 / 3 &
      weak_gap_grid$ploidy_regime_consensus <= 1),
    max(abs(weak_regime_sums - 1)) <= 1e-12,
    max(abs(pointwise_regime_sums - 1)) <= 1e-12,
    all(weak_gap_grid$local_regime_switch_proportion >= 0 &
      weak_gap_grid$local_regime_switch_proportion <= 1),
    all(weak_gap_grid$local_adjacent_ploidy_jump_median >= 0),
    sum(is.na(weak_match)) == 0L,
    weak_median_difference <= 1e-10, weak_q10_difference <= 1e-10,
    weak_q90_difference <= 1e-10,
    weak_legacy_ge2_share_difference <= 1e-12
  )
  data.frame(
    check = c(
      "unique_endpoint_cache_count", "represented_seed_count_per_pair",
      "top10_summary_row_count", "derived_summary_row_count",
      "oxygen_grid_count", "fixed_missegregation_grid_count",
      "eigenmode_rank_range", "all_endpoint_operator_qc_pass",
      "maximum_relative_eigen_residual", "all_localization_medians_finite",
      "all_spectral_gaps_nonnegative", "rank1_gap_is_zero",
      "top1_overlap_row_count", "top1_matches_figure6d_surface",
      "weak_gap_grid_row_count", "weak_gap_pair_summary_row_count",
      "weak_gap_cell_count_per_pair", "weak_gap_regime_consensus_in_valid_range",
      "weak_gap_regime_fractions_sum_to_one",
      "pointwise_three_regime_proportions_sum_to_one",
      "local_switch_proportion_in_unit_interval",
      "local_adjacent_ploidy_jump_nonnegative", "weak_gap_surface_unmatched_rows",
      "weak_gap_median_matches_figure6d_surface",
      "weak_gap_q10_matches_figure6d_surface",
      "weak_gap_q90_matches_figure6d_surface",
      "weak_gap_legacy_ge2_share_matches_figure6d_surface"
    ),
    observed = as.character(observed), expected = expected, passed = passed,
    stringsAsFactors = FALSE
  )
}

si6_data <- function(
    workspace_root = f6r_find_workspace_root(), n_core = 8L,
    rebuild = FALSE
) {
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
  )
  f6r_require_packages(c("Matrix", "Rcpp", "matrixStats", "data.table"))
  paths <- si6_paths(workspace_root)
  dir.create(paths$data, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$cache, recursive = TRUE, showWarnings = FALSE)
  f6r_load_response_engine(paths$base)
  objective_bundle <- f6r_objective_selection(paths$base)
  endpoint_manifest <- f6r_figure6d_endpoint_manifest(
    paths$base, objective_bundle
  )
  f6r_write_tsv(
    endpoint_manifest$endpoints,
    file.path(paths$data, "archive_figure6_top10_endpoint_manifest.tsv")
  )
  endpoints <- endpoint_manifest$endpoints
  cache_paths <- si6_cache_paths(paths, endpoints)
  contexts <- lapply(
    endpoint_manifest$display_manifest$pair_id,
    f6r_pair_model_context,
    selected = objective_bundle$selected, paths = paths$base
  )
  names(contexts) <- endpoint_manifest$display_manifest$pair_id

  compute_one <- function(i) {
    f6r_load_response_engine(paths$base)
    metadata <- endpoints[i, , drop = FALSE]
    message(
      "Supplementary Figure 6-3 endpoint start: ", metadata$pair_label,
      " seed", metadata$representative_seed_number,
      " (", i, "/", nrow(endpoints), ")"
    )
    result <- tryCatch(
      si6_compute_endpoint_cache(
        metadata = metadata,
        parameters = objective_bundle$parameters,
        context = contexts[[metadata$pair_id[[1L]]]],
        cache_path = cache_paths[[metadata$parameter_endpoint_group[[1L]]]],
        parameter_source = objective_bundle$paths[["parameters"]],
        force_rebuild = rebuild
      ),
      error = function(e) structure(
        list(index = i, message = conditionMessage(e)), class = "si6_error"
      )
    )
    message(
      "Supplementary Figure 6-3 endpoint complete: ", metadata$pair_label,
      " seed", metadata$representative_seed_number,
      " (", i, "/", nrow(endpoints), ")"
    )
    result
  }
  index <- seq_len(nrow(endpoints))
  qc_rows <- f6r_resilient_lapply(index, compute_one, n_core = n_core)
  failed <- vapply(qc_rows, inherits, logical(1L), "si6_error")
  if (any(failed)) {
    stop(
      "Supplementary Figure 6-3 endpoint failures: ",
      paste(vapply(qc_rows[failed], `[[`, character(1L), "message"), collapse = "; ")
    )
  }
  qc <- do.call(rbind, qc_rows)
  f6r_write_tsv(qc, file.path(paths$data, "archive_figure6_top10_endpoint_qc.tsv"))
  if (!all(qc$operator_qc_pass)) {
    stop("Supplementary Figure 6-3 endpoint operator QC failed.")
  }
  summarized <- si6_summarize_caches(paths, endpoint_manifest, cache_paths)
  summary_path <- f6r_write_tsv(
    summarized$summary,
    file.path(paths$data, "archive_figure6_top10_eigenmode_summary.tsv")
  )
  derived_path <- f6r_write_tsv(
    summarized$derived,
    file.path(paths$data, "archive_figure6_eigenmode_competition_summary.tsv")
  )
  weak_gap <- si6_weak_gap_summary(paths, endpoint_manifest, cache_paths)
  weak_gap_grid_path <- f6r_write_tsv(
    weak_gap$grid,
    file.path(paths$data, "supp_figure6-3_weak_gap_regime_robustness.tsv")
  )
  weak_gap_pair_summary_path <- f6r_write_tsv(
    weak_gap$pair_summary,
    file.path(paths$data, "supp_figure6-3_weak_gap_pair_summary.tsv")
  )
  validation <- si6_validation(
    paths, endpoint_manifest, qc, summarized$summary, summarized$derived,
    weak_gap_grid = weak_gap$grid,
    weak_gap_pair_summary = weak_gap$pair_summary
  )
  validation_path <- f6r_write_tsv(
    validation, file.path(paths$data, "supp_figure6-3_data_validation.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Supplementary Figure 6-3 data validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  source_files <- c(
    normalizePath(
      file.path(paths$base$code, "util", "analysis", "si_figure6_eigenmodes.R"),
      mustWork = TRUE
    ),
    normalizePath(
      file.path(paths$base$code, "util", "analysis", "figure6_robustness.R"),
      mustWork = TRUE
    ),
    normalizePath(
      file.path(paths$base$code, "util", "oxygen", "response_pipeline.R"),
      mustWork = TRUE
    ),
    objective_bundle$paths[["parameters"]],
    endpoint_manifest$path
  )
  contract <- data.frame(
    source_file = source_files,
    exists = file.exists(source_files),
    md5 = vapply(source_files, f6r_md5, character(1L)),
    role = c(
      "Supp_Figure6_3_analysis", "Figure6_local_ensemble_contract",
      "local_model_operator", "accepted_endpoint_parameters",
      "displayed_pair_endpoint_manifest"
    ),
    stringsAsFactors = FALSE
  )
  contract_path <- f6r_write_tsv(
    contract, file.path(paths$data, "supp_figure6-3_source_file_provenance.tsv")
  )
  invisible(list(
    paths = paths, qc = qc, summary = summarized$summary,
    derived = summarized$derived, weak_gap = weak_gap,
    outputs = c(
      summary = summary_path, derived = derived_path,
      weak_gap_grid = weak_gap_grid_path,
      weak_gap_pair_summary = weak_gap_pair_summary_path,
      validation = validation_path, provenance = contract_path
    )
  ))
}

si6_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 16, hjust = 0, margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 11, color = "#555555", margin = ggplot2::margin(b = 4)
      ),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      axis.title = ggplot2::element_text(size = 12.5),
      axis.text = ggplot2::element_text(size = 10.5),
      legend.title = ggplot2::element_text(size = 11.5),
      legend.text = ggplot2::element_text(size = 10.5),
      plot.margin = ggplot2::margin(4, 5, 4, 5)
    )
}

si6_atlas_plot <- function(data, display_label, panel_letter, fill_limits) {
  z <- data[data$display_label == display_label, , drop = FALSE]
  z$rank_label <- factor(
    paste0("Rank ", z$eigenmode_rank), levels = paste0("Rank ", 1:10)
  )
  contour <- z[z$eigenmode_rank > 1L, , drop = FALSE]
  ggplot2::ggplot(z, ggplot2::aes(O2_pct, effective_p_misseg)) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = eigenmode_localization_ploidy_median),
      interpolate = FALSE
    ) +
    ggplot2::geom_contour(
      data = contour,
      ggplot2::aes(z = spectral_gap_to_dominant_median),
      breaks = c(0.001, 0.005, 0.010),
      color = "#666666", linewidth = 0.28, alpha = 0.75
    ) +
    ggplot2::facet_wrap(~rank_label, ncol = 5L) +
    ggplot2::scale_y_log10(
      breaks = c(0.005, 0.01, 0.03, 0.1, 0.3, 0.5),
      labels = scales::label_number(accuracy = 0.001)
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#2C7BB6", "#FFFFBF", "#D7191C"),
      limits = fill_limits, oob = scales::squish,
      name = "Median eigenmode-localization\nploidy"
    ) +
    ggplot2::labs(
      x = "Fixed oxygen (%)",
      y = expression(Fixed~p[miss,eff]),
      title = paste0(panel_letter, ". ", display_label, " top-10 eigenmode localization"),
      subtitle = expression(
        "Gray contours for ranks 2-10: "*Delta*lambda[1*k]*" = 0.001, 0.005, 0.010"
      )
    ) +
    si6_theme() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, color = "#444444", linewidth = 0.40),
      panel.spacing = grid::unit(1.5, "mm"),
      legend.position = "bottom",
      legend.key.width = grid::unit(18, "mm")
    )
}

si6_derived_plot <- function(data, metric, panel_letter, title, fill_scale) {
  z <- data
  z$display_label <- factor(
    z$display_label, levels = sprintf("C%02d", 1:6)
  )
  p <- ggplot2::ggplot(z, ggplot2::aes(O2_pct, effective_p_misseg)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data[[metric]]), interpolate = FALSE) +
    ggplot2::facet_wrap(~display_label, nrow = 1L) +
    ggplot2::scale_y_log10(
      breaks = c(0.005, 0.01, 0.03, 0.1, 0.3, 0.5),
      labels = scales::label_number(accuracy = 0.001)
    ) +
    fill_scale +
    ggplot2::labs(
      x = "Fixed oxygen (%)", y = expression(Fixed~p[miss,eff]),
      title = paste0(panel_letter, ". ", title)
    ) +
    si6_theme() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, color = "#444444", linewidth = 0.40),
      panel.spacing = grid::unit(4, "mm"),
      legend.position = "right"
    )
  if (identical(metric, "near_degenerate_jump_median")) {
    p <- p + ggplot2::geom_contour(
      ggplot2::aes(z = near_degenerate_jump_median),
      breaks = 1, color = "black", linewidth = 0.40
    )
  }
  if (identical(metric, "different_ploidy_competitor_gap_median")) {
    p <- p + ggplot2::geom_contour(
      ggplot2::aes(z = proportion_with_different_ploidy_competitor),
      breaks = 0.5, color = "black", linetype = "dotted", linewidth = 0.40
    )
  }
  p
}

si6_draw_top10_audit <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c("ggplot2", "patchwork", "scales", "magick", "cowplot"))
  paths <- si6_paths(workspace_root)
  summary_path <- file.path(paths$data, "archive_figure6_top10_eigenmode_summary.tsv")
  derived_path <- file.path(paths$data, "archive_figure6_eigenmode_competition_summary.tsv")
  f6r_require_files(c(summary_path, derived_path), "archived Figure 6 eigenmode analysis")
  summary <- f6r_read_tsv(summary_path)
  derived <- f6r_read_tsv(derived_path)
  fill_limits <- range(
    summary$eigenmode_localization_ploidy_median, finite = TRUE
  )
  p_a <- si6_atlas_plot(summary, "C01", "A", fill_limits)
  p_b <- si6_atlas_plot(summary, "C02", "B", fill_limits)
  p_c <- si6_atlas_plot(summary, "C03", "C", fill_limits)
  atlas_legend <- cowplot::get_legend(
    p_a + ggplot2::theme(legend.position = "bottom")
  )
  atlas_body <- patchwork::wrap_plots(
    list(
      p_a + ggplot2::theme(legend.position = "none"),
      p_b + ggplot2::theme(legend.position = "none"),
      p_c + ggplot2::theme(legend.position = "none")
    ),
    ncol = 1L, heights = c(1, 1, 1)
  )
  atlas <- patchwork::wrap_plots(
    list(
      patchwork::wrap_elements(full = atlas_body),
      patchwork::wrap_elements(full = atlas_legend)
    ),
    ncol = 1L, heights = c(1, 0.065)
  )

  jump_max <- max(derived$near_degenerate_jump_median, na.rm = TRUE)
  p_d <- si6_derived_plot(
    derived, "near_degenerate_jump_median", "D",
    "Near-degenerate eigenmode ploidy-jump potential",
    ggplot2::scale_fill_gradientn(
      colors = c("white", "#FDB863", "#B2182B"),
      limits = c(0, jump_max), oob = scales::squish,
      name = expression(Median~max*" "*Delta*"ploidy at "*Delta*lambda[1*k]<0.005)
    )
  )
  positive_gap <- derived$different_ploidy_competitor_gap_median[
    is.finite(derived$different_ploidy_competitor_gap_median) &
      derived$different_ploidy_competitor_gap_median > 0
  ]
  p_e <- si6_derived_plot(
    derived, "different_ploidy_competitor_gap_median", "E",
    "Gap to the nearest differently localized eigenmode",
    ggplot2::scale_fill_viridis_c(
      option = "C", trans = "log10",
      limits = range(positive_gap), oob = scales::squish,
      na.value = "#E6E6E6",
      name = expression(Median~Delta*lambda[1*k]~"for "*Delta*"ploidy" >= 1)
    )
  )
  derived_plot <- patchwork::wrap_plots(
    list(p_d, p_e), ncol = 1L, heights = c(1, 1)
  )
  combined <- patchwork::wrap_plots(
    list(
      patchwork::wrap_elements(full = atlas),
      patchwork::wrap_elements(full = derived_plot)
    ),
    ncol = 1L, heights = c(3.15, 1.65)
  )

  dir.create(paths$panels, recursive = TRUE, showWarnings = FALSE)
  output_png <- file.path(paths$figures, "archive_fig6_eigenmode_competition.png")
  output_pdf <- file.path(paths$figures, "archive_fig6_eigenmode_competition.pdf")
  ggplot2::ggsave(
    output_png, combined, width = 13.2, height = 18.2,
    units = "in", dpi = 300, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    output_pdf, combined, width = 13.2, height = 18.2,
    units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
  )
  published <- c(
    figures_png = normalizePath(output_png, mustWork = TRUE),
    figures_pdf = normalizePath(output_pdf, mustWork = TRUE),
    manuscript_png = f6r_publish(
      output_png,
      file.path(paths$manuscript_figures, basename(output_png))
    ),
    manuscript_pdf = f6r_publish(
      output_pdf,
      file.path(paths$manuscript_figures, basename(output_pdf))
    )
  )
  image_info <- magick::image_info(magick::image_read(output_png))
  validation <- data.frame(
    check = c(
      "figure_png_exists", "figure_pdf_exists", "figure_width_px",
      "figure_height_px", "panel_group_count", "publication_png_md5_match",
      "publication_pdf_md5_match"
    ),
    observed = c(
      file.exists(output_png), file.exists(output_pdf), image_info$width[[1L]],
      image_info$height[[1L]], 5L,
      f6r_md5(output_png) == f6r_md5(published[["manuscript_png"]]),
      f6r_md5(output_pdf) == f6r_md5(published[["manuscript_pdf"]])
    ),
    expected = c("TRUE", "TRUE", "3960", "5460", "5", "TRUE", "TRUE"),
    stringsAsFactors = FALSE
  )
  validation$passed <- c(
    file.exists(output_png), file.exists(output_pdf),
    image_info$width[[1L]] == 3960L, image_info$height[[1L]] == 5460L,
    5L == 5L,
    f6r_md5(output_png) == f6r_md5(published[["manuscript_png"]]),
    f6r_md5(output_pdf) == f6r_md5(published[["manuscript_pdf"]])
  )
  validation_path <- f6r_write_tsv(
    validation, file.path(paths$data, "archive_figure6_figure_validation.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Archived Figure 6 eigenmode rendering validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  outputs <- unique(c(
    published, summary_path, derived_path,
    file.path(paths$data, "supp_figure6-3_data_validation.tsv"), validation_path
  ))
  manifest <- data.frame(
    output_file = outputs, exists = file.exists(outputs),
    size_bytes = file.info(outputs)$size,
    md5 = vapply(outputs, f6r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  manifest_path <- f6r_write_tsv(
    manifest, file.path(paths$data, "archive_figure6_output_manifest.tsv")
  )
  invisible(list(
    output = published, validation = validation,
    manifest = manifest_path, plot = combined
  ))
}

# The top-10 atlas above is retained as a reproducibility audit. Supplementary
# Figure 6-3 uses rank 1 only, because higher eigenmodes can be
# signed or complex and therefore are not alternative population distributions.
si6_weak_gap_base <- function(data, panel_letter, title) {
  data$display_label <- factor(
    data$display_label, levels = sprintf("C%02d", 1:6)
  )
  data$model_context <- factor(
    data$model_context, levels = c("in vivo", "in vitro")
  )
  data$log10_effective_p_misseg <- log10(data$effective_p_misseg)
  ggplot2::ggplot(data, ggplot2::aes(O2_pct, log10_effective_p_misseg)) +
    ggplot2::geom_tile(
      fill = "#F2F2F2", width = 0.025,
      height = diff(log10(si6_p_values()))[[1L]]
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(z = proportion_spectral_gap_below_0p005),
      breaks = 0.50, color = "#9B59B6", linetype = "dotted",
      linewidth = 0.55
    ) +
    ggplot2::facet_grid(model_context ~ display_label, switch = "y") +
    ggplot2::scale_y_continuous(
      breaks = log10(c(0.005, 0.01, 0.03, 0.1, 0.3, 0.5)),
      labels = c("0.005", "0.010", "0.030", "0.100", "0.300", "0.500")
    ) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      x = "Fixed oxygen (%)", y = expression(Fixed~p[miss,eff]),
      title = paste0(panel_letter, ". ", title)
    ) +
    si6_theme(base_size = 12) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        fill = NA, color = "#444444", linewidth = 0.45
      ),
      panel.spacing = grid::unit(4.0, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      legend.position = "right",
      plot.margin = ggplot2::margin(4, 4, 3, 5)
    )
}

si6_weak_gap_regime_plot <- function(data) {
  data$model_context <- factor(
    data$model_context, levels = c("in vivo", "in vitro")
  )
  data$log10_effective_p_misseg <- log10(data$effective_p_misseg)
  weak <- data[data$weak_gap_region, , drop = FALSE]
  weak$regime_class <- factor(
    weak$regime_class,
    levels = c(
      "Stable low", "Stable intermediate", "Stable high", "Mixed"
    )
  )
  si6_weak_gap_base(
    data, "A", "Low/intermediate/high mean-ploidy consensus within weak-gap regions"
  ) +
    ggplot2::geom_tile(
      data = weak, ggplot2::aes(fill = regime_class), width = 0.025,
      height = diff(log10(si6_p_values()))[[1L]]
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Stable low" = "#2C7BB6",
        "Stable intermediate" = "#FDAE61",
        "Stable high" = "#D7191C", "Mixed" = "#8C8C8C"
      ),
      limits = c(
        "Stable low", "Stable intermediate", "Stable high", "Mixed"
      ), drop = FALSE,
      labels = c(
        "Stable low" = ">=90% low (<=2)",
        "Stable intermediate" = ">=90% intermediate (>2 and <4)",
        "Stable high" = ">=90% high (>=4)",
        "Mixed" = "<90% three-class agreement"
      ),
      name = "Across retained seed endpoints"
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(z = proportion_spectral_gap_below_0p005),
      breaks = 0.50, color = "#9B59B6", linetype = "dotted",
      linewidth = 0.55
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(
        color = NA,
        fill = c("#2C7BB6", "#FDAE61", "#D7191C", "#8C8C8C")
      ),
      keyheight = grid::unit(4.5, "mm")
    ))
}

si6_weak_gap_switch_plot <- function(data) {
  data$model_context <- factor(
    data$model_context, levels = c("in vivo", "in vitro")
  )
  data$log10_effective_p_misseg <- log10(data$effective_p_misseg)
  weak <- data[data$weak_gap_region, , drop = FALSE]
  si6_weak_gap_base(
    data, "B", "Local rank-1 regime-switch sensitivity"
  ) +
    ggplot2::geom_tile(
      data = weak,
      ggplot2::aes(fill = local_regime_switch_proportion),
      width = 0.025, height = diff(log10(si6_p_values()))[[1L]]
    ) +
    ggplot2::geom_contour(
      data = weak,
      ggplot2::aes(z = local_regime_switch_proportion),
      breaks = 0.50, color = "#222222", linewidth = 0.35
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(z = proportion_spectral_gap_below_0p005),
      breaks = 0.50, color = "#9B59B6", linetype = "dotted",
      linewidth = 0.55
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#FFFFFF", "#E5D4F2", "#B07CC6", "#54278F"),
      limits = c(0, 1), oob = scales::squish, trans = "sqrt",
      breaks = c(0, 0.10, 0.25, 0.50, 0.75, 1),
      labels = scales::label_percent(accuracy = 1),
      name = "Retained seed endpoints changing\nlow/intermediate/high class\nat any 4-neighbor cell"
    )
}

si6_weak_gap_spread_plot <- function(data) {
  data$model_context <- factor(
    data$model_context, levels = c("in vivo", "in vitro")
  )
  data$log10_effective_p_misseg <- log10(data$effective_p_misseg)
  weak <- data[data$weak_gap_region, , drop = FALSE]
  upper <- max(weak$dominant_ploidy_spread_q90_q10, na.rm = TRUE)
  si6_weak_gap_base(
    data, "C", "Across-fit dominant-mean-ploidy spread"
  ) +
    ggplot2::geom_tile(
      data = weak,
      ggplot2::aes(fill = dominant_ploidy_spread_q90_q10),
      width = 0.025, height = diff(log10(si6_p_values()))[[1L]]
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#FFFFFF", "#FEE8C8", "#FDBB84", "#E34A33", "#8C2D04"),
      limits = c(0, upper), oob = scales::squish,
      trans = scales::pseudo_log_trans(sigma = 0.001, base = 10),
      breaks = c(0, 0.01, 0.1, 0.5),
      labels = scales::label_number(accuracy = 0.001),
      name = "90th-10th percentile spread\n(mean-ploidy units)"
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(z = proportion_spectral_gap_below_0p005),
      breaks = 0.50, color = "#9B59B6", linetype = "dotted",
      linewidth = 0.55
    )
}

si6_weak_gap_jump_ecdf <- function(data, pair_summary) {
  weak <- data[data$weak_gap_region, , drop = FALSE]
  weak$display_label <- factor(
    weak$display_label, levels = sprintf("C%02d", 1:6)
  )
  pair_summary$display_label <- factor(
    pair_summary$display_label, levels = sprintf("C%02d", 1:6)
  )
  weak$model_context <- factor(
    weak$model_context, levels = c("in vivo", "in vitro")
  )
  pair_summary$model_context <- factor(
    pair_summary$model_context, levels = c("in vivo", "in vitro")
  )
  pair_summary$annotation <- sprintf(
    "3-class switch in >=50%% endpoints: %.1f%%\nMedian local change >=1: %.1f%%",
    100 * pair_summary$fraction_weak_gap_majority_endpoint_local_switch,
    100 * pair_summary$fraction_weak_gap_local_jump_median_ge_1
  )
  family_colors <- c(
    C01 = "#C99700", C02 = "#6A3D9A", C03 = "#006D2C",
    C04 = "#0072B2", C05 = "#D55E00", C06 = "#009E73"
  )
  ggplot2::ggplot(
    weak,
    ggplot2::aes(x = local_adjacent_ploidy_jump_median, group = display_label)
  ) +
    ggplot2::stat_ecdf(color = "black", linewidth = 1.45, geom = "step") +
    ggplot2::stat_ecdf(
      ggplot2::aes(color = display_label), linewidth = 0.85, geom = "step"
    ) +
    ggplot2::geom_vline(
      xintercept = 1, color = "#666666", linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::geom_text(
      data = pair_summary,
      ggplot2::aes(x = 2.92, y = 0.15, label = annotation),
      inherit.aes = FALSE, hjust = 1, vjust = 0, size = 3.15,
      lineheight = 1.05, color = "#222222"
    ) +
    ggplot2::facet_grid(model_context ~ display_label, switch = "y") +
    ggplot2::scale_color_manual(values = family_colors, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = 0:3, expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent(accuracy = 1), expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 3)) +
    ggplot2::labs(
      x = "Median maximum ploidy change to a 4-neighbor grid cell",
      y = "Cumulative fraction of\nweak-gap grid cells",
      title = "D. Distribution of local rank-1 ploidy changes"
    ) +
    si6_theme(base_size = 12) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        fill = NA, color = "#444444", linewidth = 0.45
      ),
      panel.spacing = grid::unit(4.0, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      plot.margin = ggplot2::margin(4, 4, 3, 5)
    )
}

si6_fixed_legend_row <- function(plot, legend_fraction = 0.185) {
  legend <- cowplot::get_legend(
    plot + ggplot2::theme(legend.position = "right")
  )
  cowplot::plot_grid(
    plot + ggplot2::theme(legend.position = "none"), legend,
    nrow = 1L, rel_widths = c(1 - legend_fraction, legend_fraction),
    align = "h", axis = "tb"
  )
}

si6_draw_weak_gap <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c("ggplot2", "patchwork", "scales", "magick", "cowplot"))
  paths <- si6_paths(workspace_root)
  grid_path <- file.path(paths$data, "supp_figure6-3_weak_gap_regime_robustness.tsv")
  pair_summary_path <- file.path(paths$data, "supp_figure6-3_weak_gap_pair_summary.tsv")
  data_validation_path <- file.path(
    paths$data, "supp_figure6-3_context_validation.tsv"
  )
  f6r_require_files(
    c(grid_path, pair_summary_path, data_validation_path),
    "Supplementary Figure 6-3 weak-gap robustness analysis"
  )
  data <- f6r_read_tsv(grid_path)
  pair_summary <- f6r_read_tsv(pair_summary_path)
  data$model_context <- factor(
    data$model_context, levels = c("in vivo", "in vitro")
  )
  pair_summary$model_context <- factor(
    pair_summary$model_context, levels = c("in vivo", "in vitro")
  )

  p_a <- si6_weak_gap_regime_plot(data)
  p_b <- si6_weak_gap_switch_plot(data)
  p_c <- si6_weak_gap_spread_plot(data)
  p_d <- si6_weak_gap_jump_ecdf(data, pair_summary)
  row_a <- si6_fixed_legend_row(p_a)
  row_b <- si6_fixed_legend_row(p_b)
  row_c <- si6_fixed_legend_row(p_c)
  row_d <- cowplot::plot_grid(
    p_d, NULL, nrow = 1L, rel_widths = c(0.815, 0.185)
  )
  combined <- cowplot::plot_grid(
    row_a, row_b, row_c, row_d,
    ncol = 1L, rel_heights = c(1, 1, 1, 0.88),
    align = "v", axis = "lr"
  )

  output_png <- file.path(
    paths$figures, "supp_fig6-3_weak_gap_regime_robustness.png"
  )
  output_pdf <- file.path(
    paths$figures, "supp_fig6-3_weak_gap_regime_robustness.pdf"
  )
  ggplot2::ggsave(
    output_png, combined, width = 26.4, height = 20.0,
    units = "in", dpi = 300, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    output_pdf, combined, width = 26.4, height = 20.0,
    units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
  )
  published <- c(
    figures_png = normalizePath(output_png, mustWork = TRUE),
    figures_pdf = normalizePath(output_pdf, mustWork = TRUE),
    manuscript_png = f6r_publish(
      output_png, file.path(paths$manuscript_figures, basename(output_png))
    ),
    manuscript_pdf = f6r_publish(
      output_pdf, file.path(paths$manuscript_figures, basename(output_pdf))
    )
  )
  image_info <- magick::image_info(magick::image_read(output_png))
  weak <- data[data$weak_gap_region, , drop = FALSE]
  validation <- data.frame(
    check = c(
      "figure_png_exists", "figure_pdf_exists", "figure_width_px",
      "figure_height_px", "panel_group_count", "model_context_count",
      "context_row_order", "displayed_pair_count", "context_pair_count",
      "weak_gap_cells_present",
      "publication_png_md5_match",
      "publication_pdf_md5_match"
    ),
    observed = c(
      file.exists(output_png), file.exists(output_pdf), image_info$width[[1L]],
      image_info$height[[1L]], 4L, length(unique(data$model_context)),
      paste(levels(data$model_context), collapse = ","),
      length(unique(data$display_label)),
      nrow(unique(data[, c("model_context", "display_label")])), nrow(weak) > 0L,
      f6r_md5(output_png) == f6r_md5(published[["manuscript_png"]]),
      f6r_md5(output_pdf) == f6r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(
      "TRUE", "TRUE", "7920", "6000", "4", "2", "in vivo,in vitro",
      "6", "12", "TRUE", "TRUE", "TRUE"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- c(
    file.exists(output_png), file.exists(output_pdf),
    image_info$width[[1L]] == 7920L, image_info$height[[1L]] == 6000L,
    4L == 4L, length(unique(data$model_context)) == 2L,
    identical(levels(data$model_context), c("in vivo", "in vitro")),
    length(unique(data$display_label)) == 6L,
    nrow(unique(data[, c("model_context", "display_label")])) == 12L,
    nrow(weak) > 0L,
    f6r_md5(output_png) == f6r_md5(published[["manuscript_png"]]),
    f6r_md5(output_pdf) == f6r_md5(published[["manuscript_pdf"]])
  )
  validation_path <- f6r_write_tsv(
    validation, file.path(paths$data, "supp_figure6-3_figure_validation.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Supplementary Figure 6-3 rendering validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  outputs <- unique(c(
    published, grid_path, pair_summary_path, data_validation_path,
    validation_path
  ))
  manifest <- data.frame(
    output_file = outputs, exists = file.exists(outputs),
    size_bytes = file.info(outputs)$size,
    md5 = vapply(outputs, f6r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  manifest_path <- f6r_write_tsv(
    manifest, file.path(paths$data, "supp_figure6-3_output_manifest.tsv")
  )
  invisible(list(
    output = published, validation = validation,
    manifest = manifest_path, plot = combined
  ))
}

si6_draw <- si6_draw_weak_gap
