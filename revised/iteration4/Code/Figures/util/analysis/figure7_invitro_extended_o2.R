#!/usr/bin/env Rscript

# Supplementary Figure 7-4: in-vitro fixed-oxygen response diagnostics on an
# extended 0--20% O2 grid. This workflow deliberately writes to its own data
# root and never overwrites the canonical 0--5% Figure 7 products.

options(stringsAsFactors = FALSE, warn = 1)

s64_o2_values <- function() seq(0, 20, by = 0.1)
s64_p_values <- function() f7r_figure7d_p_values()
s64_cluster_labels <- function() f7r_family_levels()
s64_pair_labels <- s64_cluster_labels
s64_profile <- function() "figure7_fixed_pmisseg_invitro_extended_o2_0to20_step0p1_v1"

s64_endpoint_manifest <- function(paths) {
  manifest_path <- file.path(
    paths$base$figure7, "figure7_invitro_dense_endpoint_manifest.tsv"
  )
  endpoints <- f7r_read_tsv(manifest_path)
  endpoints <- endpoints[
    endpoints$pair_label %in% s64_pair_labels() &
      endpoints$endpoint_multiplicity_q10 > 0, , drop = FALSE
  ]
  endpoints <- endpoints[order(
    match(endpoints$display_label, s64_cluster_labels()),
    endpoints$representative_objective_rank,
    endpoints$representative_seed_number
  ), , drop = FALSE]
  if (nrow(endpoints) < f7r_family_count() ||
      !identical(sort(unique(endpoints$display_label)), s64_cluster_labels()) ||
      !identical(
        as.integer(tapply(
          endpoints$endpoint_multiplicity_q10,
          factor(endpoints$display_label, levels = s64_cluster_labels()), sum
        )),
        rep(50L, f7r_family_count())
      )) {
    stop("The extended-O2 endpoint manifest must represent 50 q10 seeds in each primary family.")
  }
  endpoints
}

s64_paths <- function(workspace_root = f7r_find_workspace_root()) {
  base <- f7r_paths(workspace_root)
  data <- file.path(base$root, "data", "Figures", "Supp_Figure7_4")
  list(
    base = base,
    data = data,
    separate_cache = file.path(data, "cache", "separate_invitro"),
    joint_cache = file.path(data, "cache", "joint_q10_surface"),
    dense_cache = file.path(data, "cache", "joint_q10_dense"),
    panels = file.path(data, "panels"),
    rendered = file.path(data, "rendered"),
    output_base = "supp_fig7-4_extended_invitro_o2_response"
  )
}

s64_analysis_paths <- function(paths) {
  out <- paths$base
  out$figure7 <- paths$data
  out
}

s64_inverse_paths <- function(paths) {
  out <- s64_analysis_paths(paths)
  out$root <- paths$data
  out
}

s64_require_grid <- function(values = s64_o2_values()) {
  expected <- seq(0, 20, by = 0.1)
  if (!identical(as.numeric(values), as.numeric(expected)) ||
      length(values) != 201L) {
    stop("Supplementary Figure 7-4 requires 201 O2 values from 0 to 20 by 0.1.")
  }
  invisible(values)
}

s64_context_with_extended_o2 <- function(context) {
  context$config$o2_S0_upper_bound <- 20
  context
}

s64_model_signature <- function(paths) {
  source_files <- f7x_model_source_files(paths$base)
  f7r_require_files(source_files, "Supplementary Figure 7-4 model source")
  paste(vapply(source_files, f7r_md5, character(1L)), collapse = "|")
}

s64_write_provenance <- function(paths) {
  source_files <- f7x_model_source_files(paths$base)
  grid <- s64_o2_values()
  provenance <- data.frame(
    artifact = c(rep("model source", length(source_files)), "oxygen grid"),
    source_file = c(normalizePath(source_files, mustWork = TRUE), NA_character_),
    md5 = c(vapply(source_files, f7r_md5, character(1L)), NA_character_),
    model_context = "in vitro",
    oxygen_min_pct = min(grid),
    oxygen_max_pct = max(grid),
    oxygen_interval_pct = unique(diff(grid))[[1L]],
    oxygen_count = length(grid),
    interpretation = paste0(
      "Post-fit extended-range model diagnostic; values above 5% O2 are ",
      "extrapolations of existing fitted parameter vectors."
    ),
    stringsAsFactors = FALSE
  )
  f7r_write_tsv(
    provenance,
    file.path(paths$data, "supp_figure7-4_source_provenance.tsv")
  )
}

s64_objective_bundle_from_frozen <- function(paths) {
  acceptance_path <- file.path(
    paths$base$figure7, "joint_seed_acceptance.tsv"
  )
  parameter_path <- file.path(
    paths$base$figure7, "joint_acceptable_seed_invitro_parameters.tsv"
  )
  f7r_require_files(
    c(acceptance_path, parameter_path),
    "frozen Figure 7 q10 selection and in-vitro parameter inputs"
  )
  objectives <- f7r_read_tsv(acceptance_path)
  objectives <- objectives[
    objectives$pair_label %in% s64_pair_labels(), , drop = FALSE
  ]
  selected <- f7r_selected_results(paths$base)
  selected <- selected[selected$pair_label %in% s64_pair_labels(), , drop = FALSE]
  parameters <- f7r_read_tsv(parameter_path)
  parameters <- parameters[
    parameters$pair_id %in% selected$pair_id, , drop = FALSE
  ]
  expected_parameter_rows <-
    f7r_family_count() * 100L * length(f7r_shared_parameters())
  if (nrow(objectives) != f7r_family_count() * 500L ||
      any(table(objectives$pair_label) != 500L) ||
      any(tapply(objectives$eligible_q10, objectives$pair_label, sum) != 50L) ||
      any(!objectives$hard_qc_pass) ||
      nrow(selected) != f7r_family_count() ||
      nrow(parameters) != expected_parameter_rows ||
      any(table(interaction(
        parameters$pair_id, parameters$seed_number, drop = TRUE
      )) != length(f7r_shared_parameters())) ||
      any(!is.finite(parameters$value))) {
    stop("Frozen primary-family Figure 7 objective bundle failed integrity checks.")
  }
  list(
    objectives = objectives, selected = selected,
    parameters_invitro = parameters,
    paths = c(
      objectives = normalizePath(acceptance_path, mustWork = TRUE),
      parameters_invitro = normalizePath(parameter_path, mustWork = TRUE)
    )
  )
}

s64_cache_valid <- function(
    object, profile, model_signature, expected_rows, source_md5 = NULL
) {
  if (is.null(object) || is.null(object$metadata) || is.null(object$qc)) {
    return(FALSE)
  }
  cache_table <- if (!is.null(object$surface)) object$surface else object$curve
  valid <- identical(object$metadata$profile, profile) &&
    identical(object$metadata$model_signature, model_signature) &&
    nrow(cache_table) == expected_rows &&
    isTRUE(object$qc$operator_qc_pass[[1L]])
  if (!is.null(source_md5)) {
    valid <- valid && identical(object$metadata$parameter_source_md5, source_md5)
  }
  valid
}

s64_prepare_run <- function(values, config, simulation = "invitro") {
  config$o2_S0_upper_bound <- 20
  run_params <- prepare_run_params(
    values, simulation = simulation, cfg = config, fixed_o2 = 20
  )
  sim_cfg <- prepare_sim_cfg(
    config, argv = list(), fixed_o2 = 20, run_params = run_params
  )
  run_params$O2_growth <- isTRUE(sim_cfg$O2_growth)
  run_params$ploidy_O2_death <- sim_cfg$ploidy_O2_death
  list(run_params = run_params, config = sim_cfg)
}

s64_separate_seed_cache <- function(
    seed_number, endpoints, cache_path, model_signature,
    parameter_source_md5, rebuild = FALSE
) {
  grid <- s64_require_grid()
  profile <- paste0(s64_profile(), "_separate")
  if (file.exists(cache_path) && !isTRUE(rebuild)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (s64_cache_valid(
      cached, profile, model_signature, length(grid), parameter_source_md5
    )) return(cached$qc)
  }
  z <- endpoints[endpoints$seed_number == seed_number, , drop = FALSE]
  if (nrow(z) != length(f7r_shared_parameters()) ||
      !setequal(z$parameter, f7r_shared_parameters())) {
    stop("Incomplete separate in-vitro endpoint for seed ", seed_number)
  }
  values <- stats::setNames(as.numeric(z$value), z$parameter)
  prepared <- s64_prepare_run(
    values, config = list(o2_S0_upper_bound = 20), simulation = "invitro"
  )
  rows <- lapply(grid, function(o2) {
    result <- fixo2_dominant_attractor_one(
      seed_id = paste0("seed", seed_number),
      run_params = prepared$run_params,
      model_env = globalenv(), cfg = prepared$config, O2 = o2
    )
    data.frame(
      seed_id = paste0("seed", seed_number), seed_number = seed_number,
      O2_pct = as.numeric(result$O2_pct[[1L]]),
      population_average_p_misseg =
        as.numeric(result$population_average_p_misseg[[1L]]),
      dominant_mean_ploidy = as.numeric(result$dominant_mean_ploidy[[1L]]),
      spectral_gap = as.numeric(result$spectral_gap[[1L]]),
      dominant_growth_rate = as.numeric(result$dominant_growth_rate[[1L]]),
      status = as.character(result$status[[1L]]),
      stringsAsFactors = FALSE
    )
  })
  curve <- do.call(rbind, rows)
  oxygen_error <- max(abs(curve$O2_pct - grid))
  finite <- all(is.finite(as.matrix(curve[, c(
    "population_average_p_misseg", "dominant_mean_ploidy",
    "spectral_gap", "dominant_growth_rate"
  )])))
  valid <- nrow(curve) == length(grid) && all(curve$status == "ok") &&
    finite && is.finite(oxygen_error) && oxygen_error <= 1e-12
  qc <- data.frame(
    seed_number = seed_number, n_o2 = nrow(curve),
    minimum_o2 = min(curve$O2_pct), maximum_o2 = max(curve$O2_pct),
    oxygen_interval = unique(round(diff(curve$O2_pct), 12))[[1L]],
    maximum_requested_o2_error = oxygen_error,
    all_status_ok = all(curve$status == "ok"),
    all_numeric_outputs_finite = finite,
    operator_qc_pass = valid,
    cache_path = normalizePath(cache_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  f7x_atomic_save_rds(
    list(
      metadata = list(
        profile = profile, seed_number = seed_number,
        objective = unique(z$objective), model_signature = model_signature,
        parameter_source_md5 = parameter_source_md5,
        oxygen_values = grid
      ),
      curve = curve, qc = qc
    ),
    cache_path
  )
  qc$cache_path <- normalizePath(cache_path, mustWork = TRUE)
  qc
}

s64_compute_separate <- function(paths, n_core = 8L, rebuild = FALSE) {
  f7r_require_packages(c("Matrix", "Rcpp", "data.table"))
  analysis_paths <- s64_analysis_paths(paths)
  f7r_load_response_engine(analysis_paths)
  endpoint_bundle <- f7x_separate_invitro_endpoint_table(analysis_paths)
  endpoints <- endpoint_bundle$data
  parameter_source_md5 <- f7r_md5(endpoint_bundle$path)
  model_signature <- s64_model_signature(paths)
  dir.create(paths$separate_cache, recursive = TRUE, showWarnings = FALSE)
  seeds <- sort(unique(endpoints$seed_number))
  compute_one <- function(seed_number) {
    f7r_load_response_engine(analysis_paths)
    if (seed_number %% 25L == 0L) message("Extended separate in-vitro seed ", seed_number, "/500")
    tryCatch(
      s64_separate_seed_cache(
        seed_number = seed_number, endpoints = endpoints,
        cache_path = file.path(paths$separate_cache, paste0("seed", seed_number, ".rds")),
        model_signature = model_signature,
        parameter_source_md5 = parameter_source_md5,
        rebuild = rebuild
      ),
      error = function(e) structure(
        list(seed_number = seed_number, message = conditionMessage(e)),
        class = "s64_error"
      )
    )
  }
  result <- f7r_resilient_lapply(seeds, compute_one, n_core = n_core)
  failed <- vapply(result, inherits, logical(1L), "s64_error")
  if (any(failed)) stop(
    "Extended separate in-vitro endpoint failures: ",
    paste(vapply(result[failed], `[[`, character(1L), "message"), collapse = "; ")
  )
  qc <- do.call(rbind, result)
  if (nrow(qc) != 500L || !all(qc$operator_qc_pass)) {
    stop("Extended separate in-vitro endpoint QC failed.")
  }
  objects <- lapply(qc$cache_path, readRDS)
  curves <- do.call(rbind, lapply(objects, `[[`, "curve"))
  curves <- curves[order(curves$seed_number, curves$O2_pct), , drop = FALSE]
  if (nrow(curves) != 500L * 201L || anyDuplicated(
    curves[, c("seed_number", "O2_pct")]
  )) stop("Extended separate in-vitro grid is incomplete or duplicated.")
  classified <- response_classify_all_curves(curves)
  by_seed <- classified$summary
  smooth <- classified$smooth
  segments <- classified$segments
  objective_map <- unique(endpoints[, c("seed_number", "objective")])
  by_seed <- merge(by_seed, objective_map, by = "seed_number", all.x = TRUE)
  gap_rows <- lapply(split(curves$spectral_gap, curves$seed_number), function(gap) {
    data.frame(
      fraction_gap_below_0p005 = mean(gap < 0.005),
      fraction_gap_below_0p01 = mean(gap < 0.01),
      any_gap_below_0p005 = any(gap < 0.005),
      stringsAsFactors = FALSE
    )
  })
  gap_summary <- do.call(rbind, gap_rows)
  gap_summary$seed_number <- as.integer(rownames(gap_summary))
  rownames(gap_summary) <- NULL
  gap_summary$spectral_gap_class <- ifelse(
    gap_summary$fraction_gap_below_0p005 >= 0.25, "unreliable",
    ifelse(
      gap_summary$any_gap_below_0p005 |
        gap_summary$fraction_gap_below_0p01 >= 0.10,
      "caution", "reliable"
    )
  )
  by_seed <- merge(by_seed, gap_summary, by = "seed_number", all.x = TRUE)
  by_seed$context <- "in vitro"
  smooth$context <- "in vitro"
  if (nrow(segments)) segments$context <- "in vitro"
  class_counts <- data.frame(
    smooth_curve_class = response_curve_class_order,
    n_seed = as.integer(table(factor(
      by_seed$smooth_curve_class, levels = response_curve_class_order
    ))), stringsAsFactors = FALSE
  )
  class_counts$fraction_seed <- class_counts$n_seed / 500
  outputs <- c(
    curves = f7r_write_tsv(curves, file.path(paths$data, "separate_invitro_raw_curves.tsv")),
    smooth = f7r_write_tsv(smooth, file.path(paths$data, "separate_invitro_smoothed_curves.tsv")),
    by_seed = f7r_write_tsv(by_seed, file.path(paths$data, "separate_invitro_curve_class_by_seed.tsv")),
    segments = f7r_write_tsv(segments, file.path(paths$data, "separate_invitro_persistent_segments.tsv")),
    counts = f7r_write_tsv(class_counts, file.path(paths$data, "separate_invitro_class_counts.tsv")),
    qc = f7r_write_tsv(qc, file.path(paths$data, "separate_invitro_operator_qc.tsv"))
  )
  list(
    curves = curves, smooth = smooth, by_seed = by_seed,
    segments = segments, class_counts = class_counts, qc = qc,
    paths = outputs
  )
}

s64_extended_endpoint_cache <- function(
    metadata, parameters, context, base_cache_path, cache_path,
    parameter_source, p_profile = c("standard", "dense"),
    model_signature, rebuild = FALSE
) {
  p_profile <- match.arg(p_profile)
  o2_values <- s64_require_grid()
  extended_o2 <- o2_values[o2_values > 5]
  profile <- paste0(s64_profile(), "_joint_", p_profile)
  parameter_source_md5 <- f7r_md5(parameter_source)
  base_cache_md5 <- f7r_md5(base_cache_path)
  f7r_require_files(base_cache_path, "canonical 0--5 in-vitro endpoint cache")
  base <- readRDS(base_cache_path)
  p_values <- sort(unique(base$surface$p_misseg))
  expected_p <- if (p_profile == "standard") 60L else 496L
  expected_rows <- length(o2_values) * expected_p
  if (length(p_values) != expected_p) {
    stop("Unexpected ", p_profile, " p_misseg grid in ", base_cache_path)
  }
  if (file.exists(cache_path) && !isTRUE(rebuild)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (s64_cache_valid(
      cached, profile, model_signature, expected_rows, parameter_source_md5
    ) && identical(cached$metadata$base_cache_md5, base_cache_md5) &&
        (p_profile != "standard" || nrow(cached$trajectory) == 201L)) {
      return(cached$qc)
    }
  }

  seed_number <- metadata$representative_seed_number[[1L]]
  z <- parameters[
    parameters$pair_id == metadata$pair_id[[1L]] &
      parameters$seed_number == seed_number, , drop = FALSE
  ]
  required <- f7r_shared_parameters()
  if (nrow(z) != length(required) || !setequal(z$parameter, required) ||
      any(!is.finite(z$value))) {
    stop("Incomplete extended in-vitro parameters for ", metadata$pair_label)
  }
  params <- stats::setNames(as.numeric(z$value), z$parameter)
  prepared <- s64_prepare_run(
    params, config = context$config, simulation = "invitro"
  )
  seed_id <- paste0("seed", seed_number)

  base_o2 <- o2_values[o2_values <= 5]
  base_keep <- vapply(
    base$surface$O2_pct,
    function(x) any(abs(x - base_o2) <= 1e-12), logical(1L)
  )
  base_surface <- base$surface[base_keep, , drop = FALSE]
  required_columns <- c(
    "O2_pct", "p_misseg", "forced_p_misseg", "status",
    "population_average_p_misseg", "dominant_mean_ploidy",
    "spectral_gap", "dominant_growth_rate"
  )
  if (!all(required_columns %in% names(base_surface)) ||
      nrow(base_surface) != length(base_o2) * expected_p) {
    stop("Canonical base cache cannot provide the 0--5 overlap grid.")
  }
  base_surface <- base_surface[, required_columns, drop = FALSE]
  extension_rows <- lapply(p_values, function(p_fixed) {
    forced <- figure7_force_p_misseg(
      prepared$run_params, p_fixed
    )
    rows <- lapply(extended_o2, function(o2) {
      result <- fixo2_dominant_attractor_one(
        seed_id = seed_id, run_params = forced, model_env = globalenv(),
        cfg = prepared$config, O2 = o2
      )
      data.frame(
        O2_pct = as.numeric(result$O2_pct[[1L]]),
        p_misseg = p_fixed,
        forced_p_misseg = as.numeric(forced$p_misseg),
        status = as.character(result$status[[1L]]),
        population_average_p_misseg =
          as.numeric(result$population_average_p_misseg[[1L]]),
        dominant_mean_ploidy =
          as.numeric(result$dominant_mean_ploidy[[1L]]),
        spectral_gap = as.numeric(result$spectral_gap[[1L]]),
        dominant_growth_rate =
          as.numeric(result$dominant_growth_rate[[1L]]),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  extension <- do.call(rbind, extension_rows)
  surface <- rbind(base_surface, extension)
  surface <- surface[order(
    surface$p_misseg, surface$O2_pct
  ), , drop = FALSE]
  surface$model_context <- "in vitro"
  trajectory <- NULL
  if (p_profile == "standard") {
    trajectory_columns <- c(
      "O2_pct", "status", "population_average_p_misseg",
      "dominant_mean_ploidy", "spectral_gap", "dominant_growth_rate",
      "fitted_p_misseg", "fitted_p_mis_base", "fitted_k_o_mis"
    )
    if (is.null(base$trajectory) ||
        !all(trajectory_columns %in% names(base$trajectory))) {
      stop("Canonical standard cache lacks its natural-p trajectory.")
    }
    trajectory_base <- base$trajectory[vapply(
      base$trajectory$O2_pct,
      function(x) any(abs(x - base_o2) <= 1e-12), logical(1L)
    ), trajectory_columns, drop = FALSE]
    trajectory_extension <- do.call(rbind, lapply(extended_o2, function(o2) {
      result <- fixo2_dominant_attractor_one(
        seed_id = seed_id, run_params = prepared$run_params,
        model_env = globalenv(), cfg = prepared$config, O2 = o2
      )
      data.frame(
        O2_pct = as.numeric(result$O2_pct[[1L]]),
        status = as.character(result$status[[1L]]),
        population_average_p_misseg =
          as.numeric(result$population_average_p_misseg[[1L]]),
        dominant_mean_ploidy =
          as.numeric(result$dominant_mean_ploidy[[1L]]),
        spectral_gap = as.numeric(result$spectral_gap[[1L]]),
        dominant_growth_rate =
          as.numeric(result$dominant_growth_rate[[1L]]),
        fitted_p_misseg = as.numeric(prepared$run_params$p_misseg),
        fitted_p_mis_base = as.numeric(prepared$run_params$p_mis_base),
        fitted_k_o_mis = as.numeric(prepared$run_params$k_o_mis),
        stringsAsFactors = FALSE
      )
    }))
    trajectory <- rbind(trajectory_base, trajectory_extension)
    trajectory <- trajectory[order(trajectory$O2_pct), , drop = FALSE]
    trajectory$model_context <- "in vitro"
  }
  oxygen_error <- max(abs(
    surface$O2_pct - rep(o2_values, times = expected_p)
  ))
  p_error <- max(abs(surface$forced_p_misseg - surface$p_misseg), na.rm = TRUE)
  formula_qc <- figure7_p_misseg_formula_qc(prepared$run_params, p_values)
  finite <- all(is.finite(as.matrix(surface[, c(
    "population_average_p_misseg", "dominant_mean_ploidy",
    "spectral_gap", "dominant_growth_rate"
  )])))
  trajectory_valid <- p_profile != "standard" ||
    (nrow(trajectory) == 201L && all(trajectory$status == "ok") &&
       all(is.finite(as.matrix(trajectory[, c(
         "population_average_p_misseg", "dominant_mean_ploidy",
         "spectral_gap", "dominant_growth_rate"
       )]))))
  valid <- nrow(surface) == expected_rows && all(surface$status == "ok") &&
    finite && oxygen_error <= 1e-12 && p_error <= 1e-8 &&
    isTRUE(formula_qc$passed) && trajectory_valid
  qc <- data.frame(
    display_label = metadata$display_label,
    pair_label = metadata$pair_label, pair_id = metadata$pair_id,
    model_context = "in vitro",
    parameter_endpoint_group = metadata$parameter_endpoint_group,
    representative_seed_number = seed_number,
    endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10,
    p_profile = p_profile, n_fixed_p = expected_p,
    maximum_direct_formula_error = formula_qc$maximum_direct_formula_error,
    n_o2_per_fixed_p = length(o2_values), n_surface = nrow(surface),
    minimum_o2 = min(surface$O2_pct), maximum_o2 = max(surface$O2_pct),
    oxygen_interval = unique(round(diff(sort(unique(surface$O2_pct))), 12))[[1L]],
    maximum_requested_o2_error = oxygen_error,
    all_status_ok = all(surface$status == "ok"),
    all_numeric_outputs_finite = finite,
    max_abs_forced_minus_requested_p_misseg = p_error,
    trajectory_qc_pass = trajectory_valid, operator_qc_pass = valid,
    cache_path = normalizePath(cache_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  f7x_atomic_save_rds(
    list(
      metadata = list(
        profile = profile, model_context = "in vitro",
        display_label = metadata$display_label[[1L]],
        pair_label = metadata$pair_label[[1L]],
        pair_id = metadata$pair_id[[1L]],
        parameter_endpoint_group = metadata$parameter_endpoint_group[[1L]],
        representative_seed_number = seed_number,
        endpoint_multiplicity_q10 = metadata$endpoint_multiplicity_q10[[1L]],
        model_signature = model_signature,
        parameter_source_md5 = parameter_source_md5,
        base_cache_md5 = base_cache_md5,
        oxygen_values = o2_values, p_misseg_values = p_values,
        interpretation = paste0(
          "O2 <= 5 copied from the canonical validated endpoint cache; ",
          "O2 > 5 calculated with o2_S0_upper_bound = 20."
        )
      ),
      surface = surface, trajectory = trajectory, qc = qc
    ),
    cache_path
  )
  qc$cache_path <- normalizePath(cache_path, mustWork = TRUE)
  qc
}

s64_summarize_joint_trajectory <- function(manifest, cache_paths) {
  f7r_require_packages("matrixStats")
  rows <- lapply(manifest$display_manifest$pair_id, function(pair) {
    metadata <- manifest$endpoints[
      manifest$endpoints$pair_id == pair, , drop = FALSE
    ]
    objects <- lapply(
      cache_paths[metadata$parameter_endpoint_group], readRDS
    )
    reference <- objects[[1L]]$trajectory
    if (is.null(reference)) stop("Missing extended natural-p trajectory.")
    same_grid <- vapply(objects, function(x) {
      identical(x$trajectory$O2_pct, reference$O2_pct)
    }, logical(1L))
    if (!all(same_grid)) stop("Extended trajectory grids do not match.")
    object_groups <- vapply(
      objects,
      function(x) x$metadata$parameter_endpoint_group, character(1L)
    )
    multiplicity <- metadata$endpoint_multiplicity_q10[
      match(object_groups, metadata$parameter_endpoint_group)
    ]
    if (sum(multiplicity) != 50L) stop("Trajectory weights do not sum to 50.")
    expand_index <- rep(seq_along(objects), times = multiplicity)
    summarize_field <- function(field) {
      matrix_values <- do.call(cbind, lapply(
        objects, function(x) x$trajectory[[field]]
      ))[, expand_index, drop = FALSE]
      quantiles <- matrixStats::rowQuantiles(
        matrix_values, probs = c(0.10, 0.50, 0.90)
      )
      stats::setNames(
        as.data.frame(quantiles), paste0(field, c("_q10", "_median", "_q90"))
      )
    }
    out <- data.frame(
      pair_id = pair, pair_label = metadata$pair_label[[1L]],
      display_label = metadata$display_label[[1L]],
      O2_pct = reference$O2_pct, n_seed = 50L,
      n_unique_parameter_endpoint = nrow(metadata), stringsAsFactors = FALSE
    )
    for (field in c(
      "fitted_p_misseg", "fitted_p_mis_base", "fitted_k_o_mis",
      "population_average_p_misseg", "dominant_mean_ploidy",
      "spectral_gap", "dominant_growth_rate"
    )) {
      out <- cbind(out, summarize_field(field))
    }
    out$model_context <- "in vitro"
    out
  })
  out <- do.call(rbind, rows)
  out <- out[order(
    match(out$display_label, s64_cluster_labels()), out$O2_pct
  ), , drop = FALSE]
  rownames(out) <- NULL
  out
}

s64_compute_joint_profile <- function(
    paths, objective_bundle, p_profile = c("standard", "dense"),
    n_core = 8L, rebuild = FALSE
) {
  p_profile <- match.arg(p_profile)
  f7r_require_packages(c("Matrix", "Rcpp", "matrixStats", "data.table"))
  analysis_paths <- s64_analysis_paths(paths)
  manifest <- f7x_joint_context_endpoint_manifest(
    analysis_paths, objective_bundle, cutoff = "q10", displayed_only = TRUE,
    output_name = paste0("joint_invitro_q10_", p_profile, "_endpoint_manifest.tsv")
  )
  endpoints <- manifest$endpoints
  if (nrow(endpoints) < f7r_family_count() ||
      !identical(sort(unique(endpoints$display_label)), s64_cluster_labels()) ||
      any(tapply(
        endpoints$endpoint_multiplicity_q10, endpoints$display_label, sum
      ) != 50L)) {
    stop("Extended joint profile requires validated q10 endpoints representing 50 seeds in each primary family.")
  }
  contexts <- lapply(
    manifest$display_manifest$pair_id, f7r_pair_model_context,
    selected = objective_bundle$selected, paths = paths$base
  )
  names(contexts) <- manifest$display_manifest$pair_id
  cache_root <- if (p_profile == "standard") paths$joint_cache else paths$dense_cache
  base_root <- file.path(
    paths$base$figure7,
    if (p_profile == "standard") {
      "multiseed_endpoint_cache_invitro"
    } else {
      "figure7_invitro_dense_endpoint_cache"
    }
  )
  cache_paths <- stats::setNames(
    file.path(
      cache_root, endpoints$pair_label,
      paste0("endpoint_", endpoints$parameter_endpoint_group, ".rds")
    ), endpoints$parameter_endpoint_group
  )
  base_paths <- stats::setNames(
    file.path(
      base_root, endpoints$pair_label,
      paste0("endpoint_", endpoints$parameter_endpoint_group, ".rds")
    ), endpoints$parameter_endpoint_group
  )
  f7r_require_files(unname(base_paths), paste0("canonical ", p_profile, " caches"))
  model_signature <- s64_model_signature(paths)
  parameter_source <- objective_bundle$paths[["parameters_invitro"]]
  compute_one <- function(i) {
    f7r_load_response_engine(analysis_paths)
    z <- endpoints[i, , drop = FALSE]
    if (i %% 10L == 0L || i == 1L) {
      message("Extended ", p_profile, " endpoint ", i, "/", nrow(endpoints))
    }
    tryCatch(
      s64_extended_endpoint_cache(
        metadata = z, parameters = objective_bundle$parameters_invitro,
        context = contexts[[z$pair_id[[1L]]]],
        base_cache_path = base_paths[[z$parameter_endpoint_group[[1L]]]],
        cache_path = cache_paths[[z$parameter_endpoint_group[[1L]]]],
        parameter_source = parameter_source, p_profile = p_profile,
        model_signature = model_signature, rebuild = rebuild
      ),
      error = function(e) structure(
        list(index = i, message = conditionMessage(e)), class = "s64_error"
      )
    )
  }
  result <- f7r_resilient_lapply(
    seq_len(nrow(endpoints)), compute_one, n_core = n_core
  )
  failed <- vapply(result, inherits, logical(1L), "s64_error")
  if (any(failed)) stop(
    "Extended ", p_profile, " endpoint failures: ",
    paste(vapply(result[failed], `[[`, character(1L), "message"), collapse = "; ")
  )
  qc <- do.call(rbind, result)
  if (!all(qc$operator_qc_pass)) stop("Extended joint cache QC failed.")
  summary <- f7r_figure7d_summarize_dense_caches(
    analysis_paths, manifest, cache_paths
  )
  summary$model_context <- "in vitro"
  trajectory <- NULL
  if (p_profile == "standard") {
    trajectory <- s64_summarize_joint_trajectory(manifest, cache_paths)
  }
  summary_path <- f7r_write_tsv(
    summary,
    file.path(paths$data, paste0("joint_invitro_", p_profile, "_surface.tsv"))
  )
  qc_path <- f7r_write_tsv(
    qc, file.path(paths$data, paste0("joint_invitro_", p_profile, "_operator_qc.tsv"))
  )
  trajectory_path <- if (!is.null(trajectory)) {
    f7r_write_tsv(
      trajectory,
      file.path(paths$data, "joint_invitro_standard_trajectory.tsv")
    )
  } else {
    character()
  }
  list(
    manifest = manifest, cache_paths = cache_paths, qc = qc,
    summary = summary, trajectory = trajectory,
    paths = c(summary = summary_path, qc = qc_path,
              trajectory = trajectory_path)
  )
}

s64_profile_cache_complete <- function(paths, p_profile) {
  cache_root <- if (p_profile == "standard") paths$joint_cache else paths$dense_cache
  files <- list.files(
    cache_root, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE
  )
  expected_endpoint_count <- nrow(s64_endpoint_manifest(paths))
  if (length(files) != expected_endpoint_count) return(FALSE)
  objects <- lapply(files, function(path) {
    tryCatch(readRDS(path), error = function(e) NULL)
  })
  expected_rows <- 201L * if (p_profile == "standard") 60L else 496L
  expected_profile <- paste0(s64_profile(), "_joint_", p_profile)
  model_signature <- s64_model_signature(paths)
  parameter_source_md5 <- f7r_md5(file.path(
    paths$base$figure7, "joint_acceptable_seed_invitro_parameters.tsv"
  ))
  base_root <- file.path(
    paths$base$figure7,
    if (p_profile == "standard") {
      "multiseed_endpoint_cache_invitro"
    } else {
      "figure7_invitro_dense_endpoint_cache"
    }
  )
  all(vapply(seq_along(objects), function(index) {
    object <- objects[[index]]
    if (is.null(object) || is.null(object$metadata) || is.null(object$surface) ||
        is.null(object$qc)) return(FALSE)
    base_path <- file.path(
      base_root, object$metadata$pair_label,
      paste0("endpoint_", object$metadata$parameter_endpoint_group, ".rds")
    )
    file.exists(base_path) && nrow(object$surface) == expected_rows &&
      identical(object$metadata$profile, expected_profile) &&
      identical(object$metadata$model_signature, model_signature) &&
      identical(object$metadata$parameter_source_md5, parameter_source_md5) &&
      identical(object$metadata$base_cache_md5, f7r_md5(base_path)) &&
      isTRUE(object$qc$operator_qc_pass[[1L]])
  }, logical(1L)))
}

s64_run_endpoint_workers <- function(
    paths, p_profile, n_core = 8L, rebuild = FALSE
) {
  if (!isTRUE(rebuild) && s64_profile_cache_complete(paths, p_profile)) {
    message("Extended ", p_profile, " endpoint cache is complete; workers skipped.")
    return(invisible(TRUE))
  }
  worker <- file.path(
    paths$base$code, "worker_Supp_Figure7_4_endpoint.R"
  )
  f7r_require_files(worker, "Supplementary Figure 7-4 endpoint worker")
  endpoint_count <- nrow(s64_endpoint_manifest(paths))
  workers <- max(1L, min(as.integer(n_core), endpoint_count))
  requested_plan <- tolower(trimws(Sys.getenv(
    "FIGURE7_FUTURE_PLAN", unset = "multisession"
  )))
  if (identical(requested_plan, "multicore")) {
    # On the Linux compute node, load the external model once and fork from
    # that process.  Starting one independent R session per endpoint makes all
    # workers contend for the external sourceCpp lock and can time out before
    # numerical evaluation begins.
    f7r_load_response_engine(paths$base)
    bundle <- s64_objective_bundle_from_frozen(paths)
    endpoints <- s64_endpoint_manifest(paths)
    contexts <- lapply(
      unique(endpoints$pair_id), f7r_pair_model_context,
      selected = bundle$selected, paths = paths$base
    )
    names(contexts) <- unique(endpoints$pair_id)
    base_root <- file.path(
      paths$base$figure7,
      if (p_profile == "standard") {
        "multiseed_endpoint_cache_invitro"
      } else {
        "figure7_invitro_dense_endpoint_cache"
      }
    )
    cache_root <- if (p_profile == "standard") {
      paths$joint_cache
    } else {
      paths$dense_cache
    }
    model_signature <- s64_model_signature(paths)
    parameter_source <- bundle$paths[["parameters_invitro"]]
    compute_one <- function(index) {
      z <- endpoints[index, , drop = FALSE]
      base_cache_path <- file.path(
        base_root, z$pair_label[[1L]],
        paste0("endpoint_", z$parameter_endpoint_group[[1L]], ".rds")
      )
      cache_path <- file.path(
        cache_root, z$pair_label[[1L]],
        paste0("endpoint_", z$parameter_endpoint_group[[1L]], ".rds")
      )
      qc <- s64_extended_endpoint_cache(
        metadata = z, parameters = bundle$parameters_invitro,
        context = contexts[[z$pair_id[[1L]]]],
        base_cache_path = base_cache_path, cache_path = cache_path,
        parameter_source = parameter_source, p_profile = p_profile,
        model_signature = model_signature, rebuild = rebuild
      )
      if (!isTRUE(qc$operator_qc_pass[[1L]])) {
        stop("Endpoint QC failed for index ", index, ".")
      }
      message(
        p_profile, " endpoint ", index, "/", endpoint_count,
        " complete: ", z$display_label[[1L]], " ",
        z$parameter_endpoint_group[[1L]]
      )
      qc
    }
    message(
      "Launching ", workers, " fork workers for extended ",
      p_profile, " endpoints."
    )
    results <- f7r_resilient_lapply(
      seq_len(endpoint_count), compute_one, n_core = workers
    )
    if (any(vapply(
      results, function(x) inherits(x, "try-error"), logical(1L)
    )) || !s64_profile_cache_complete(paths, p_profile)) {
      stop("Forked extended ", p_profile, " endpoint workers failed QC.")
    }
    return(invisible(TRUE))
  }
  command <- paste(
    "seq 1", endpoint_count, "| xargs -P", workers, "-I{}",
    shQuote(file.path(R.home("bin"), "Rscript")), shQuote(worker),
    paste0("--profile=", p_profile), "--index={}",
    paste0("--rebuild=", if (isTRUE(rebuild)) "true" else "false")
  )
  message(
    "Launching ", workers, " independent R workers for extended ",
    p_profile, " endpoints."
  )
  status <- system(command)
  if (!identical(as.integer(status), 0L) ||
      !s64_profile_cache_complete(paths, p_profile)) {
    stop("Independent extended ", p_profile, " endpoint workers failed QC.")
  }
  invisible(TRUE)
}

s64_weak_gap_pair_summary <- function(grid) {
  rows <- lapply(split(grid, grid$display_label), function(z) {
    weak <- z[z$weak_gap_region, , drop = FALSE]
    data.frame(
      display_label = z$display_label[[1L]], pair_label = z$pair_label[[1L]],
      n_grid_cell = nrow(z), n_weak_gap_cell = nrow(weak),
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
      weak_gap_ploidy_spread_max = max(weak$dominant_ploidy_spread_q90_q10),
      weak_gap_local_jump_median = stats::median(
        weak$local_adjacent_ploidy_jump_median
      ),
      weak_gap_local_jump_q90 = stats::quantile(
        weak$local_adjacent_ploidy_jump_median, 0.90, names = FALSE
      ),
      weak_gap_local_jump_max = max(weak$local_adjacent_ploidy_jump_median),
      fraction_weak_gap_local_jump_median_ge_1 = mean(
        weak$local_adjacent_ploidy_jump_median >= 1
      ),
      model_context = "in vitro", stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[match(s64_cluster_labels(), out$display_label), , drop = FALSE]
}

s64_compute_weak_gap <- function(paths, standard) {
  f7r_require_packages("matrixStats")
  rows <- lapply(standard$manifest$display_manifest$pair_id, function(pair) {
    metadata <- standard$manifest$endpoints[
      standard$manifest$endpoints$pair_id == pair, , drop = FALSE
    ]
    raw <- lapply(
      standard$cache_paths[metadata$parameter_endpoint_group], readRDS
    )
    objects <- Map(
      f7x_si7_rank1_object, raw, metadata$parameter_endpoint_group
    )
    si7_summarize_weak_gap_pair(
      metadata, objects, figure7_surface = standard$summary
    )
  })
  grid <- do.call(rbind, rows)
  rownames(grid) <- NULL
  grid$model_context <- "in vitro"
  pair_summary <- s64_weak_gap_pair_summary(grid)
  paths_out <- c(
    grid = f7r_write_tsv(
      grid, file.path(paths$data, "weak_gap_regime_robustness.tsv")
    ),
    pair_summary = f7r_write_tsv(
      pair_summary, file.path(paths$data, "weak_gap_pair_summary.tsv")
    )
  )
  list(grid = grid, pair_summary = pair_summary, paths = paths_out)
}

s64_key <- function(data, columns) {
  formatted <- lapply(data[, columns, drop = FALSE], function(x) {
    if (is.numeric(x)) sprintf("%.12g", x) else as.character(x)
  })
  do.call(paste, c(formatted, sep = "|"))
}

s64_overlap_metric <- function(new, old, keys, values, label) {
  new <- new[new$O2_pct <= 5 + 1e-12, , drop = FALSE]
  old_key <- s64_key(old, keys)
  new_key <- s64_key(new, keys)
  matched <- match(new_key, old_key)
  if (anyNA(matched) || anyDuplicated(old_key)) {
    stop("Cannot match 0--5 overlap rows for ", label)
  }
  deltas <- vapply(values, function(field) {
    a <- as.numeric(new[[field]])
    b <- as.numeric(old[[field]][matched])
    finite <- is.finite(a) & is.finite(b)
    if (!any(finite)) 0 else max(abs(a[finite] - b[finite]))
  }, numeric(1L))
  data.frame(
    comparison = label, n_overlap = nrow(new),
    field = values, maximum_absolute_difference = deltas,
    tolerance = 1e-8, passed = deltas <= 1e-8,
    stringsAsFactors = FALSE
  )
}

s64_validate <- function(paths, separate, standard, dense, inverse, weak_gap) {
  base <- paths$base$figure7
  original_separate <- f7r_read_tsv(file.path(
    base, "response_class_invitro_raw_curves.tsv"
  ))
  original_separate <- original_separate[vapply(
    original_separate$O2_pct,
    function(x) any(abs(x - seq(0, 5, by = 0.1)) <= 1e-12), logical(1L)
  ), , drop = FALSE]
  original_standard <- f7r_read_tsv(file.path(
    base, "joint_multiseed_surface_summary_invitro.tsv"
  ))
  original_standard <- original_standard[
    original_standard$cutoff == "q10" &
      original_standard$pair_label %in% s64_pair_labels() &
      vapply(
        original_standard$O2_pct,
        function(x) any(abs(x - seq(0, 5, by = 0.1)) <= 1e-12), logical(1L)
      ), , drop = FALSE
  ]
  original_dense <- f7r_read_tsv(file.path(
    base, "figure7_invitro_fixed_p_curve_family.tsv"
  ))
  original_dense <- original_dense[
    vapply(
      original_dense$O2_pct,
      function(x) any(abs(x - seq(0, 5, by = 0.1)) <= 1e-12), logical(1L)
    ), , drop = FALSE
  ]
  original_inverse <- f7r_read_tsv(file.path(
    base, "figure7_invitro_inverse_response_summary.tsv"
  ))
  original_inverse <- original_inverse[
    vapply(
      original_inverse$O2_pct,
      function(x) any(abs(x - seq(0, 5, by = 0.1)) <= 1e-12), logical(1L)
    ), , drop = FALSE
  ]

  overlap <- rbind(
    s64_overlap_metric(
      separate$curves, original_separate,
      c("seed_number", "O2_pct"),
      c("population_average_p_misseg", "dominant_mean_ploidy", "spectral_gap"),
      "separate in-vitro curves"
    ),
    s64_overlap_metric(
      standard$summary, original_standard,
      c("pair_id", "O2_pct", "p_misseg"),
      c("dominant_mean_ploidy_median", "spectral_gap_median"),
      "joint 60-level response surface"
    ),
    s64_overlap_metric(
      dense$summary, original_dense,
      c("pair_id", "O2_pct", "p_misseg"),
      c("dominant_mean_ploidy_median", "spectral_gap_median"),
      "joint 496-level dense surface"
    ),
    s64_overlap_metric(
      inverse$summary, original_inverse,
      c("pair_id", "O2_pct", "target_ploidy"),
      c("fraction_any_solution", "fraction_unique_solution",
        "fraction_multiple_solutions", "p_display"),
      "inverse response"
    )
  )
  overlap_path <- f7r_write_tsv(
    overlap, file.path(paths$data, "zero_to_five_overlap_validation.tsv")
  )

  grid <- s64_o2_values()
  cluster_check <- function(x) {
    identical(sort(unique(as.character(x$display_label))), s64_cluster_labels())
  }
  checks <- data.frame(
    check = c(
      "oxygen_grid_count", "oxygen_grid_range", "oxygen_grid_interval",
      "separate_seed_count", "separate_row_count", "separate_operator_qc",
      "standard_cluster_set", "standard_row_count", "standard_q10_seed_weights",
      "standard_operator_qc", "dense_cluster_set", "dense_row_count",
      "dense_q10_seed_weights", "dense_operator_qc", "inverse_row_count",
      "weak_gap_cluster_set", "weak_gap_row_count", "overlap_reproduction"
    ),
    observed = c(
      length(grid), paste(range(grid), collapse = ","),
      paste(unique(round(diff(grid), 12)), collapse = ","),
      length(unique(separate$curves$seed_number)), nrow(separate$curves),
      all(separate$qc$operator_qc_pass), cluster_check(standard$summary),
      nrow(standard$summary),
      paste(tapply(
        standard$manifest$endpoints$endpoint_multiplicity_q10,
        standard$manifest$endpoints$display_label, sum
      ), collapse = ","),
      all(standard$qc$operator_qc_pass), cluster_check(dense$summary),
      nrow(dense$summary),
      paste(tapply(
        dense$manifest$endpoints$endpoint_multiplicity_q10,
        dense$manifest$endpoints$display_label, sum
      ), collapse = ","),
      all(dense$qc$operator_qc_pass), nrow(inverse$summary),
      cluster_check(weak_gap$grid), nrow(weak_gap$grid), all(overlap$passed)
    ),
    expected = c(
      201, "0,20", "0.1", 500, 500L * 201L, TRUE, TRUE,
      f7r_family_count() * 60L * 201L,
      paste(rep(50L, f7r_family_count()), collapse = ","), TRUE, TRUE,
      f7r_family_count() * 496L * 201L,
      paste(rep(50L, f7r_family_count()), collapse = ","), TRUE,
      f7r_family_count() * 201L * 241L, TRUE,
      f7r_family_count() * 60L * 201L, TRUE
    ),
    stringsAsFactors = FALSE
  )
  checks$passed <- as.character(checks$observed) == as.character(checks$expected)
  validation_path <- f7r_write_tsv(
    checks, file.path(paths$data, "supp_figure7-4_validation.tsv")
  )
  if (!all(checks$passed) || !all(overlap$passed)) {
    stop(
      "Supplementary Figure 7-4 validation failed: ",
      paste(checks$check[!checks$passed], collapse = ", ")
    )
  }
  list(checks = checks, overlap = overlap, paths = c(
    validation = validation_path, overlap = overlap_path
  ))
}

s64_write_data_manifest <- function(paths, named_paths) {
  files <- unname(named_paths[file.exists(named_paths)])
  manifest <- data.frame(
    role = names(named_paths)[file.exists(named_paths)],
    path = normalizePath(files, mustWork = TRUE),
    md5 = vapply(files, f7r_md5, character(1L)),
    bytes = file.info(files)$size,
    stringsAsFactors = FALSE
  )
  f7r_write_tsv(
    manifest, file.path(paths$data, "supp_figure7-4_data_manifest.tsv")
  )
}

s64_data <- function(
    workspace_root = f7r_find_workspace_root(), n_core = 8L,
    rebuild = FALSE
) {
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
  )
  paths <- s64_paths(workspace_root)
  dir.create(paths$data, recursive = TRUE, showWarnings = FALSE)
  f7r_load_response_engine(paths$base)
  provenance_path <- s64_write_provenance(paths)
  message("Supplementary Figure 7-4: separate-in-vitro 500-seed curves")
  separate <- s64_compute_separate(
    paths, n_core = n_core, rebuild = rebuild
  )
  objective_bundle <- s64_objective_bundle_from_frozen(paths)
  message("Supplementary Figure 7-4: q10 60-level response surfaces")
  s64_run_endpoint_workers(
    paths, "standard", n_core = n_core, rebuild = rebuild
  )
  standard <- s64_compute_joint_profile(
    paths, objective_bundle, p_profile = "standard",
    n_core = 1L, rebuild = FALSE
  )
  message("Supplementary Figure 7-4: q10 496-level dense surfaces")
  s64_run_endpoint_workers(
    paths, "dense", n_core = n_core, rebuild = rebuild
  )
  dense <- s64_compute_joint_profile(
    paths, objective_bundle, p_profile = "dense",
    n_core = 1L, rebuild = FALSE
  )
  message("Supplementary Figure 7-4: inverse response")
  inverse <- f7r_inverse_panel_data(
    s64_inverse_paths(paths), rebuild = rebuild, n_core = n_core,
    dense_qc_path = dense$paths[["qc"]],
    output_prefix = "joint_invitro_extended_o2",
    model_context = "in vitro"
  )
  message("Supplementary Figure 7-4: weak-gap ensemble diagnostics")
  weak_gap <- s64_compute_weak_gap(paths, standard)
  validation <- s64_validate(
    paths, separate, standard, dense, inverse, weak_gap
  )
  manifest_path <- s64_write_data_manifest(paths, c(
    provenance = provenance_path, separate$paths, standard$paths,
    dense$paths, inverse$paths, weak_gap$paths, validation$paths
  ))
  invisible(list(
    paths = paths, separate = separate, standard = standard, dense = dense,
    inverse = inverse, weak_gap = weak_gap, validation = validation,
    manifest = manifest_path
  ))
}

s64_theme <- function(base_size = 8.2) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 10.2, hjust = 0, margin = ggplot2::margin(b = 3)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 7.0, colour = "#555555", margin = ggplot2::margin(b = 3)
      ),
      axis.title = ggplot2::element_text(size = 8.0),
      axis.text = ggplot2::element_text(size = 7.0, colour = "#222222"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 8.2),
      panel.border = ggplot2::element_rect(
        fill = NA, colour = "#444444", linewidth = 0.35
      ),
      panel.spacing = grid::unit(2.2, "mm"),
      legend.title = ggplot2::element_text(size = 7.1),
      legend.text = ggplot2::element_text(size = 6.7),
      legend.key.height = grid::unit(3.5, "mm"),
      plot.margin = ggplot2::margin(4, 4, 3, 4)
    )
}

s64_extended_marker <- function(indexed = FALSE) {
  boundary <- if (isTRUE(indexed)) 50 else 5
  breaks <- if (isTRUE(indexed)) c(0, 50, 100, 150, 200) else c(0, 5, 10, 15, 20)
  list(
    ggplot2::geom_vline(
      xintercept = boundary, colour = "#777777", linetype = "dashed",
      linewidth = 0.28
    ),
    ggplot2::scale_x_continuous(
      breaks = breaks, labels = c(0, 5, 10, 15, 20),
      expand = ggplot2::expansion(mult = c(0, 0))
    )
  )
}

s64_add_o2_index <- function(data) {
  data$O2_index <- round(10 * data$O2_pct, 8)
  data
}

s64_panel_a <- function(paths) {
  smooth <- f7r_read_tsv(file.path(
    paths$data, "separate_invitro_smoothed_curves.tsv"
  ))
  classes <- f7r_read_tsv(file.path(
    paths$data, "separate_invitro_curve_class_by_seed.tsv"
  ))
  merged <- merge(
    smooth,
    classes[, c("seed_id", "smooth_curve_class"), drop = FALSE],
    by = "seed_id", all.x = TRUE, sort = FALSE
  )
  counts <- table(factor(
    classes$smooth_curve_class, levels = response_curve_class_order
  ))
  class_labels <- paste0(
    unname(response_curve_class_labels[response_curve_class_order]),
    " (n=", as.integer(counts), ")"
  )
  names(class_labels) <- response_curve_class_order
  merged$class_panel <- factor(
    class_labels[merged$smooth_curve_class], levels = class_labels
  )
  merged$curve_color <- factor(
    merged$smooth_curve_class, levels = response_curve_class_order
  )
  median_curve <- stats::aggregate(
    smoothed_dominant_mean_ploidy ~ class_panel + O2_pct,
    merged, stats::median, na.rm = TRUE
  )
  ggplot2::ggplot(
    merged,
    ggplot2::aes(
      O2_pct, smoothed_dominant_mean_ploidy,
      group = seed_id, colour = curve_color
    )
  ) +
    ggplot2::geom_line(linewidth = 0.17, alpha = 0.24) +
    ggplot2::geom_line(
      data = median_curve,
      ggplot2::aes(
        O2_pct, smoothed_dominant_mean_ploidy, group = class_panel
      ),
      inherit.aes = FALSE, colour = "#111111", linewidth = 0.62
    ) +
    s64_extended_marker() +
    ggplot2::facet_grid(
      class_panel ~ ., scales = "free_y", switch = "y", drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = response_curve_class_colors, drop = FALSE
    ) +
    ggplot2::labs(
      title = "A. Separate-fit in vitro oxygen-ploidy response classes",
      subtitle = paste0(
        "500 fitted seeds; black, class median; dashed line, fitted O2 boundary"
      ),
      x = "Fixed oxygen (%)", y = "Smoothed dominant mean ploidy"
    ) +
    s64_theme(base_size = 8.0) +
    ggplot2::theme(
      legend.position = "none", strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(
        angle = 0, hjust = 1, size = 6.7, margin = ggplot2::margin(r = 3)
      ),
      panel.spacing.y = grid::unit(1.25, "mm")
    )
}

s64_standard_hatch <- function(surface) {
  surface$log10_p_misseg <- log10(surface$p_misseg)
  rows <- lapply(s64_cluster_labels(), function(label) {
    z <- surface[surface$display_label == label, , drop = FALSE]
    h <- f7r_weak_gap_hatch_data(z)
    if (!nrow(h)) return(NULL)
    h$display_label <- label
    h
  })
  rows <- Filter(Negate(is.null), rows)
  out <- if (length(rows)) do.call(rbind, rows) else data.frame(
    O2_pct = numeric(), log10_p_misseg = numeric(),
    hatch_group = integer(), display_label = character(),
    stringsAsFactors = FALSE
  )
  out$display_label <- factor(out$display_label, levels = s64_cluster_labels())
  out
}

s64_panel_b <- function(paths) {
  surface <- f7r_read_tsv(file.path(
    paths$data, "joint_invitro_standard_surface.tsv"
  ))
  trajectory <- f7r_read_tsv(file.path(
    paths$data, "joint_invitro_standard_trajectory.tsv"
  ))
  surface$display_label <- factor(
    surface$display_label, levels = s64_cluster_labels()
  )
  trajectory$display_label <- factor(
    trajectory$display_label, levels = s64_cluster_labels()
  )
  surface$log10_p_misseg <- log10(surface$p_misseg)
  surface <- s64_add_o2_index(surface)
  trajectory <- s64_add_o2_index(trajectory)
  trajectory$log10_p_median <- log10(trajectory$fitted_p_misseg_median)
  trajectory$log10_p_q10 <- log10(trajectory$fitted_p_misseg_q10)
  trajectory$log10_p_q90 <- log10(trajectory$fitted_p_misseg_q90)
  hatch <- s64_standard_hatch(surface)
  hatch <- s64_add_o2_index(hatch)
  fill_limits <- range(c(surface$dominant_mean_ploidy_median, 1, 7), na.rm = TRUE)
  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = surface,
      ggplot2::aes(
        O2_index, log10_p_misseg,
        fill = dominant_mean_ploidy_median
      )
    ) +
    ggplot2::geom_path(
      data = hatch,
      ggplot2::aes(
        O2_index, log10_p_misseg, group = hatch_group
      ),
      colour = "#9B59B6", linewidth = 0.16, alpha = 0.72
    ) +
    ggplot2::geom_ribbon(
      data = trajectory,
      ggplot2::aes(
        O2_index, ymin = log10_p_q10, ymax = log10_p_q90
      ),
      inherit.aes = FALSE, fill = "#E0E0E0", alpha = 0.62
    ) +
    ggplot2::geom_path(
      data = trajectory, ggplot2::aes(O2_index, log10_p_median),
      inherit.aes = FALSE, colour = "#111111", linewidth = 0.45
    ) +
    s64_extended_marker(indexed = TRUE) +
    ggplot2::facet_grid(. ~ display_label) +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "#FFFFBF", "#B2182B"),
      trans = "log10", limits = fill_limits,
      breaks = c(1, 1.5, 2, 3, 4, 6),
      name = "Median dominant\nmean ploidy\n(log colors)"
    ) +
    ggplot2::scale_y_continuous(
      breaks = log10(c(0.005, 0.01, 0.05, 0.10, 0.50)),
      labels = c("0.005", "0.01", "0.05", "0.10", "0.50")
    ) +
    ggplot2::coord_cartesian(ylim = log10(c(0.005, 0.5)), expand = FALSE) +
    ggplot2::labs(
      title = "B. Oxygen-missegregation-ploidy response surface",
      subtitle = paste0(
        "q10 ensemble; purple, weak-gap region; black and gray, fitted p_misseg median and 10-90%"
      ),
      x = "Fixed oxygen (%)", y = expression(p[misseg])
    ) +
    s64_theme()
}

s64_inverse_hatch <- function(inverse) {
  rows <- lapply(s64_cluster_labels(), function(label) {
    z <- inverse[inverse$display_label == label, , drop = FALSE]
    h <- f7r_inverse_multivalue_hatch_data(z)
    if (!nrow(h)) return(NULL)
    h$display_label <- label
    h
  })
  rows <- Filter(Negate(is.null), rows)
  out <- if (length(rows)) do.call(rbind, rows) else data.frame(
    O2_pct = numeric(), target_ploidy = numeric(), hatch_group = integer(),
    display_label = character(), stringsAsFactors = FALSE
  )
  out$display_label <- factor(out$display_label, levels = s64_cluster_labels())
  out
}

s64_panel_c <- function(paths) {
  inverse <- f7r_read_tsv(file.path(
    paths$data, "joint_invitro_extended_o2_inverse_response_summary.tsv"
  ))
  dense <- f7r_read_tsv(file.path(
    paths$data, "joint_invitro_dense_surface.tsv"
  ))
  inverse$display_label <- factor(
    inverse$display_label, levels = s64_cluster_labels()
  )
  dense$display_label <- factor(dense$display_label, levels = s64_cluster_labels())
  inverse <- s64_add_o2_index(inverse)
  dense <- s64_add_o2_index(dense)
  highlighted_values <- c(0.01, 0.10, 0.20, 0.30)
  highlighted <- dense[vapply(
    dense$p_misseg,
    function(x) any(abs(x - highlighted_values) <= 1e-12), logical(1L)
  ), , drop = FALSE]
  highlighted$reference_label <- factor(
    sprintf("%.2f", highlighted$p_misseg),
    levels = c("0.01", "0.10", "0.20", "0.30")
  )
  mean_curve <- stats::aggregate(
    dominant_mean_ploidy_median ~ display_label + pair_id + O2_pct,
    dense, mean
  )
  mean_curve$reference_label <- factor(
    "Mean across 496 fixed p_misseg values",
    levels = c(
      "0.01", "0.10", "0.20", "0.30",
      "Mean across 496 fixed p_misseg values"
    )
  )
  mean_curve <- s64_add_o2_index(mean_curve)
  hatch <- s64_inverse_hatch(inverse)
  hatch <- s64_add_o2_index(hatch)
  reference_levels <- levels(mean_curve$reference_label)
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = inverse,
      ggplot2::aes(O2_index, target_ploidy, fill = p_display)
    ) +
    ggplot2::geom_hline(
      yintercept = c(2, 4), colour = "#666666", linewidth = 0.22,
      linetype = "longdash"
    ) +
    ggplot2::geom_path(
      data = hatch,
      ggplot2::aes(O2_index, target_ploidy, group = hatch_group),
      colour = "#7B3294", linewidth = 0.16, alpha = 0.72
    ) +
    ggplot2::geom_path(
      data = highlighted,
      ggplot2::aes(
        O2_index, dominant_mean_ploidy_median,
        group = reference_label, linetype = reference_label
      ),
      colour = "#111111", linewidth = 0.47
    ) +
    ggplot2::geom_path(
      data = mean_curve,
      ggplot2::aes(
        O2_index, dominant_mean_ploidy_median,
        group = reference_label, linetype = reference_label
      ),
      colour = "#D62728", linewidth = 0.62
    ) +
    s64_extended_marker(indexed = TRUE) +
    ggplot2::facet_grid(. ~ display_label) +
    ggplot2::scale_fill_viridis_c(
      option = "D", trans = "log10", limits = c(0.005, 0.5),
      breaks = c(0.005, 0.01, 0.05, 0.10, 0.50),
      na.value = "#EFEFEF", name = "Median required\np_misseg\n(log colors)"
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "0.01" = "solid", "0.10" = "F28282",
        "0.20" = "dotdash", "0.30" = "dotted",
        "Mean across 496 fixed p_misseg values" = "solid"
      ),
      breaks = reference_levels,
      labels = c(
        "p_misseg = 0.01", "p_misseg = 0.10",
        "p_misseg = 0.20", "p_misseg = 0.30",
        "Mean across 496 fixed p_misseg values"
      ), name = "Reference curves"
    ) +
    ggplot2::coord_cartesian(ylim = c(1, 7), expand = FALSE) +
    ggplot2::labs(
      title = "C. p_misseg required for target ploidy",
      subtitle = paste0(
        "Gray, no stable unique inverse; purple, multiple inverse solutions"
      ),
      x = "Fixed oxygen (%)", y = "Target dominant mean ploidy"
    ) +
    s64_theme()
}

s64_weak_gap_base_plot <- function(data, title) {
  data$display_label <- factor(data$display_label, levels = s64_cluster_labels())
  data$log10_p_misseg <- log10(data$p_misseg)
  data <- s64_add_o2_index(data)
  ggplot2::ggplot(data, ggplot2::aes(O2_index, log10_p_misseg)) +
    ggplot2::geom_raster(fill = "#F2F2F2") +
    s64_extended_marker(indexed = TRUE) +
    ggplot2::facet_grid(. ~ display_label) +
    ggplot2::scale_y_continuous(
      breaks = log10(c(0.005, 0.01, 0.05, 0.10, 0.50)),
      labels = c("0.005", "0.01", "0.05", "0.10", "0.50")
    ) +
    ggplot2::coord_cartesian(ylim = log10(c(0.005, 0.5)), expand = FALSE) +
    ggplot2::labs(
      title = title, x = "Fixed oxygen (%)",
      y = expression(p[misseg])
    ) +
    s64_theme()
}

s64_panel_d <- function(paths) {
  data <- f7r_read_tsv(file.path(paths$data, "weak_gap_regime_robustness.tsv"))
  weak <- data[data$weak_gap_region, , drop = FALSE]
  weak$log10_p_misseg <- log10(weak$p_misseg)
  weak <- s64_add_o2_index(weak)
  boundary <- s64_standard_hatch(data)
  boundary <- s64_add_o2_index(boundary)
  weak$regime_class <- factor(
    weak$regime_class,
    levels = c("Stable low", "Stable intermediate", "Stable high", "Mixed")
  )
  s64_weak_gap_base_plot(
    data, "D. Ploidy consensus within weak-gap regions"
  ) +
    ggplot2::geom_tile(
      data = weak, ggplot2::aes(fill = regime_class),
      width = 1,
      height = stats::median(diff(log10(sort(unique(data$p_misseg)))))
    ) +
    ggplot2::geom_path(
      data = boundary,
      ggplot2::aes(
        O2_index, log10_p_misseg, group = hatch_group
      ),
      inherit.aes = FALSE, colour = "#9B59B6", linewidth = 0.18,
      alpha = 0.75
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Stable low" = "#2C7BB6", "Stable intermediate" = "#FDAE61",
        "Stable high" = "#D7191C", "Mixed" = "#8C8C8C"
      ),
      labels = c(
        "Stable low" = ">=90% low (<=2)",
        "Stable intermediate" = ">=90% intermediate (>2 and <4)",
        "Stable high" = ">=90% high (>=4)",
        "Mixed" = "<90% three-class agreement"
      ), name = "Across retained seed endpoints", drop = FALSE
    )
}

s64_panel_e <- function(paths) {
  data <- f7r_read_tsv(file.path(paths$data, "weak_gap_regime_robustness.tsv"))
  weak <- data[data$weak_gap_region, , drop = FALSE]
  weak$log10_p_misseg <- log10(weak$p_misseg)
  weak <- s64_add_o2_index(weak)
  boundary <- s64_standard_hatch(data)
  boundary <- s64_add_o2_index(boundary)
  s64_weak_gap_base_plot(
    data, "E. Local rank-1 regime-switch sensitivity"
  ) +
    ggplot2::geom_tile(
      data = weak, ggplot2::aes(fill = local_regime_switch_proportion),
      width = 1,
      height = stats::median(diff(log10(sort(unique(data$p_misseg)))))
    ) +
    ggplot2::geom_path(
      data = boundary,
      ggplot2::aes(
        O2_index, log10_p_misseg, group = hatch_group
      ),
      inherit.aes = FALSE, colour = "#9B59B6", linewidth = 0.18,
      alpha = 0.75
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFFFFF", "#E5D4F2", "#B07CC6", "#54278F"),
      limits = c(0, 1), oob = scales::squish, trans = "sqrt",
      labels = scales::label_percent(accuracy = 1),
      name = "Endpoints changing class\nat any 4-neighbor cell"
    )
}

s64_panel_f <- function(paths) {
  data <- f7r_read_tsv(file.path(paths$data, "weak_gap_regime_robustness.tsv"))
  weak <- data[data$weak_gap_region, , drop = FALSE]
  weak$log10_p_misseg <- log10(weak$p_misseg)
  weak <- s64_add_o2_index(weak)
  boundary <- s64_standard_hatch(data)
  boundary <- s64_add_o2_index(boundary)
  upper <- max(weak$dominant_ploidy_spread_q90_q10, na.rm = TRUE)
  s64_weak_gap_base_plot(
    data, "F. Across-fit dominant-mean-ploidy spread"
  ) +
    ggplot2::geom_tile(
      data = weak,
      ggplot2::aes(fill = dominant_ploidy_spread_q90_q10),
      width = 1,
      height = stats::median(diff(log10(sort(unique(data$p_misseg)))))
    ) +
    ggplot2::geom_path(
      data = boundary,
      ggplot2::aes(
        O2_index, log10_p_misseg, group = hatch_group
      ),
      inherit.aes = FALSE, colour = "#9B59B6", linewidth = 0.18,
      alpha = 0.75
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFFFFF", "#FEE8C8", "#FDBB84", "#E34A33", "#8C2D04"),
      limits = c(0, upper), oob = scales::squish,
      trans = scales::pseudo_log_trans(sigma = 0.001, base = 10),
      breaks = c(0, 0.01, 0.1, 0.5),
      labels = scales::label_number(accuracy = 0.001),
      name = "90th-10th percentile spread\n(mean-ploidy units)"
    )
}

s64_save_plot <- function(plot, png_path, pdf_path, width, height) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    png_path, plot = plot, width = width, height = height,
    units = "in", dpi = 300, device = "png", type = "cairo",
    bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    pdf_path, plot = plot, width = width, height = height,
    units = "in", device = grDevices::cairo_pdf, bg = "white",
    limitsize = FALSE
  )
  c(png = normalizePath(png_path), pdf = normalizePath(pdf_path))
}

s64_draw <- function(workspace_root = f7r_find_workspace_root()) {
  font_cache <- file.path(tempdir(), "supp-figure7-4-fontconfig-cache")
  dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = font_cache)
  f7r_require_packages(c(
    "ggplot2", "patchwork", "scales", "viridisLite", "isoband", "magick"
  ))
  paths <- s64_paths(workspace_root)
  f7r_load_response_engine(paths$base)
  required <- file.path(paths$data, c(
    "separate_invitro_smoothed_curves.tsv",
    "separate_invitro_curve_class_by_seed.tsv",
    "joint_invitro_standard_surface.tsv",
    "joint_invitro_standard_trajectory.tsv",
    "joint_invitro_dense_surface.tsv",
    "joint_invitro_extended_o2_inverse_response_summary.tsv",
    "weak_gap_regime_robustness.tsv",
    "supp_figure7-4_validation.tsv"
  ))
  f7r_require_files(required, "validated Supplementary Figure 7-4 data")
  validation <- f7r_read_tsv(file.path(
    paths$data, "supp_figure7-4_validation.tsv"
  ))
  if (!all(validation$passed)) stop("Refusing to draw from non-passing data.")
  dir.create(paths$panels, recursive = TRUE, showWarnings = FALSE)
  plots <- list(
    A = s64_panel_a(paths), B = s64_panel_b(paths), C = s64_panel_c(paths),
    D = s64_panel_d(paths), E = s64_panel_e(paths), F = s64_panel_f(paths)
  )
  panel_sizes <- list(
    A = c(5.3, 16.0), B = c(24.4, 3.0), C = c(24.4, 3.25),
    D = c(24.4, 2.75), E = c(24.4, 2.75), F = c(24.4, 2.75)
  )
  panel_outputs <- unlist(lapply(names(plots), function(label) {
    size <- panel_sizes[[label]]
    out <- s64_save_plot(
      plots[[label]],
      file.path(paths$panels, paste0("supp_fig7-4", tolower(label), ".png")),
      file.path(paths$panels, paste0("supp_fig7-4", tolower(label), ".pdf")),
      width = size[[1L]], height = size[[2L]]
    )
    stats::setNames(out, paste0(label, "_", names(out)))
  }))
  right <- plots$B / plots$C / plots$D / plots$E / plots$F +
    patchwork::plot_layout(heights = c(1.05, 1.12, 0.96, 0.96, 0.96))
  assembled <- plots$A | right
  assembled <- assembled + patchwork::plot_layout(widths = c(0.18, 0.82)) +
    patchwork::plot_annotation(
      title = "Supplementary Figure 7-4. Extended-range in vitro oxygen-ploidy response",
      subtitle = paste0(
        "O2 = 0-20% at 0.1% intervals. Values above 5% are post-fit model extrapolations."
      )
    )
  output_png <- file.path(paths$base$figures, paste0(paths$output_base, ".png"))
  output_pdf <- file.path(paths$base$figures, paste0(paths$output_base, ".pdf"))
  output <- s64_save_plot(
    assembled, output_png, output_pdf, width = 30.0, height = 16.0
  )
  published <- c(
    manuscript_png = f7r_publish(
      output[["png"]],
      file.path(paths$base$manuscript_figures, basename(output[["png"]]))
    ),
    manuscript_pdf = f7r_publish(
      output[["pdf"]],
      file.path(paths$base$manuscript_figures, basename(output[["pdf"]]))
    )
  )
  image <- magick::image_read(output[["png"]])
  info <- magick::image_info(image)
  draw_validation <- data.frame(
    check = c(
      "png_width", "png_height", "pdf_nontrivial",
      "manuscript_png_identical", "manuscript_pdf_identical",
      "panel_file_count"
    ),
    observed = c(
      info$width[[1L]], info$height[[1L]], file.info(output[["pdf"]])$size > 10000,
      f7r_md5(output[["png"]]) == f7r_md5(published[["manuscript_png"]]),
      f7r_md5(output[["pdf"]]) == f7r_md5(published[["manuscript_pdf"]]),
      length(panel_outputs)
    ),
    expected = c(9000, 4800, TRUE, TRUE, TRUE, 12),
    stringsAsFactors = FALSE
  )
  draw_validation$passed <- as.character(draw_validation$observed) ==
    as.character(draw_validation$expected)
  draw_validation_path <- f7r_write_tsv(
    draw_validation, file.path(paths$data, "supp_figure7-4_draw_validation.tsv")
  )
  if (!all(draw_validation$passed)) {
    stop("Supplementary Figure 7-4 rendering validation failed.")
  }
  output_manifest <- data.frame(
    role = c("assembled_png", "assembled_pdf", names(panel_outputs)),
    path = c(output, panel_outputs),
    md5 = vapply(c(output, panel_outputs), f7r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  output_manifest_path <- f7r_write_tsv(
    output_manifest, file.path(paths$data, "supp_figure7-4_output_manifest.tsv")
  )
  invisible(list(
    output = output, published = published, panels = panel_outputs,
    validation = draw_validation, validation_path = draw_validation_path,
    manifest = output_manifest_path, plot = assembled
  ))
}
