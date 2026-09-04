#!/usr/bin/env Rscript

# Full-range data for Figure 7B. In vivo is continuous only. In vitro passage
# uses fitting-compatible nearest-target selection within fixed culture windows,
# independently of continuous culture. Validated continuous caches may be reused.

options(stringsAsFactors = FALSE, warn = 1)

f7g_profile <- function() "figure7_fixed_pmisseg_full_range_q10_v3_feasible_segments"
f7g_contexts <- function() c("in vivo", "in vitro")
f7g_modes <- function(context) {
  if (identical(context, "in vivo")) "continuous" else c("continuous", "passage")
}
f7g_days <- function(smoke = FALSE) {
  if (isTRUE(smoke)) 0:10 else 0:10000
}
f7g_o2 <- function(context, smoke = FALSE) {
  if (identical(context, "in vivo")) {
    if (isTRUE(smoke)) c(0, 2.5, 5) else seq(0, 5, length.out = 201L)
  } else if (identical(context, "in vitro")) {
    if (isTRUE(smoke)) c(0, 10, 20) else seq(0, 20, by = 0.1)
  } else {
    stop("Unknown Figure 7 model context: ", context)
  }
}

f7g_paths <- function(paths, run_id = NULL, create = FALSE) {
  base <- file.path(paths$figure7, "finite_time_full_q10_runs")
  current <- file.path(paths$figure7, "finite_time_full_q10_current.tsv")
  if (is.null(run_id)) {
    f7r_require_files(current, "current Figure 7 full-range run pointer")
    pointer <- f7r_read_tsv(current)
    required <- c("run_id", "relative_run_path", "profile", "fingerprint")
    if (nrow(pointer) != 1L || !all(required %in% names(pointer))) {
      stop("Malformed Figure 7 full-range pointer: ", current)
    }
    run_id <- f7ft_sanitize_run_id(pointer$run_id[[1L]])
    run_root <- normalizePath(
      file.path(paths$figure7, pointer$relative_run_path[[1L]]),
      mustWork = TRUE
    )
  } else {
    run_id <- f7ft_sanitize_run_id(run_id)
    run_root <- file.path(base, run_id)
  }
  out <- list(
    run_id = run_id, base = base, run_root = run_root,
    cache = file.path(run_root, "task_cache"),
    panels = file.path(run_root, "panels"),
    rendered = file.path(run_root, "rendered"), current = current
  )
  if (isTRUE(create)) {
    invisible(lapply(
      unname(unlist(out[c("base", "run_root", "cache", "panels", "rendered")])),
      dir.create, recursive = TRUE, showWarnings = FALSE
    ))
  }
  out
}

f7g_panel_filename <- function(context, mode, family) {
  context_key <- if (identical(context, "in vivo")) "invivo" else "invitro"
  paste0("full_range_panel_", context_key, "_", mode, "_", tolower(family), ".rds")
}

f7g_read_panel <- function(run_paths, context, mode, family) {
  if (!context %in% f7g_contexts() || !mode %in% f7g_modes(context) ||
      !family %in% f7ft_family_levels()) {
    stop("Invalid Figure 7 full-range panel request.")
  }
  path <- file.path(
    run_paths$run_root, f7g_panel_filename(context, mode, family)
  )
  f7r_require_files(path, paste(context, mode, family, "full-range panel"))
  object <- readRDS(path)
  expected <- c(
    identical(object$profile, f7g_profile()),
    identical(object$model_context, context),
    identical(object$propagation_mode, mode),
    identical(object$pair_label, family)
  )
  if (!all(expected)) stop("Unexpected full-range panel object: ", path)
  object
}

f7g_load_propagator <- function(paths) {
  f7r_require_packages("Rcpp")
  cpp <- file.path(
    paths$code, "util", "analysis", "figure7_full_range_propagator.cpp"
  )
  f7r_require_files(cpp, "Figure 7 full-range C++ propagator")
  cache_dir <- file.path(paths$figure7, ".rcpp_cache_full_range")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(cpp, rebuild = FALSE, cacheDir = cache_dir, showOutput = FALSE)
  if (!exists("f7g_propagate_operator_cpp", mode = "function", inherits = TRUE)) {
    stop("Figure 7 full-range C++ propagator did not load.")
  }
  normalizePath(cpp, mustWork = TRUE)
}

f7g_geometric_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(NA_real_)
  exp(mean(log(x)))
}

f7g_canonical_passage_rule <- function(schedule_bundle, run_paths) {
  raw <- schedule_bundle$schedule
  valid <- raw[
    is.finite(raw$initial_cells) & raw$initial_cells > 0 &
      is.finite(raw$target_live_cells) & raw$target_live_cells > raw$initial_cells,
    , drop = FALSE
  ]
  if (length(unique(valid$lineage)) != 6L) {
    stop("Canonical passage rule requires valid maintenance passages in six lineages.")
  }
  valid$expansion_ratio <- valid$target_live_cells / valid$initial_cells
  split_valid <- split(valid, valid$lineage, drop = TRUE)
  lineage <- do.call(rbind, lapply(split_valid, function(z) data.frame(
    cohort = z$cohort[[1L]], lineage = z$lineage[[1L]],
    n_valid_passage = nrow(z),
    seed_cells_geometric_mean = f7g_geometric_mean(z$initial_cells),
    expansion_ratio_geometric_mean = f7g_geometric_mean(z$expansion_ratio),
    passage_duration_median_day = stats::median(z$passage_duration),
    stringsAsFactors = FALSE
  )))
  rownames(lineage) <- NULL
  seed <- f7g_geometric_mean(lineage$seed_cells_geometric_mean)
  ratio <- f7g_geometric_mean(lineage$expansion_ratio_geometric_mean)
  rule <- data.frame(
    rule_id = "six_lineage_equal_weight_nearest_target_feasible_segments_v3",
    n_source_lineage = nrow(lineage),
    n_source_maintenance_passage = nrow(valid),
    seed_cells = seed,
    expansion_ratio = ratio,
    target_live_cells = seed * ratio,
    representative_passage_duration_day =
      stats::median(lineage$passage_duration_median_day),
    trigger = "nearest_absolute_live_cell_target_within_positive_segment_days",
    tie_break = "earliest_candidate_day",
    time_axis = "accumulated_segment_duration_not_selected_model_day",
    inoculum_constraint = "downsample_only_or_protocol_infeasible",
    ensemble_policy = "mean_requires_all_50_endpoints_no_survivor_renormalization",
    reseeding = "composition_neutral_scalar_reset_to_seed_cells",
    stringsAsFactors = FALSE
  )
  if (nrow(lineage) != 6L || !all(is.finite(unlist(rule[c(
      "seed_cells", "expansion_ratio", "target_live_cells"
    )]))) || rule$expansion_ratio[[1L]] <= 1) {
    stop("Canonical passage-rule validation failed.")
  }
  paths <- c(
    raw = f7ft_atomic_write_tsv(
      raw, file.path(run_paths$run_root, "source_six_lineage_passage_schedule.tsv")
    ),
    valid = f7ft_atomic_write_tsv(
      valid, file.path(run_paths$run_root, "canonical_passage_source_rows.tsv")
    ),
    lineage = f7ft_atomic_write_tsv(
      lineage, file.path(run_paths$run_root, "canonical_passage_lineage_summary.tsv")
    ),
    rule = f7ft_atomic_write_tsv(
      rule, file.path(run_paths$run_root, "canonical_passage_rule.tsv")
    )
  )
  list(rule = rule, lineage = lineage, valid = valid, paths = paths,
       source_fit_result = schedule_bundle$fit_result_path)
}

f7g_task_manifest <- function(
    endpoint_manifest, run_paths, smoke = FALSE, o2_chunk_size = 10L
) {
  endpoints <- endpoint_manifest$endpoints
  p_values <- if (isTRUE(smoke)) f7ft_p_values()[1:2] else f7ft_p_values()
  rows <- list()
  index <- 0L
  for (context in f7g_contexts()) {
    oxygen <- f7g_o2(context, smoke)
    chunks <- split(
      seq_along(oxygen),
      ceiling(seq_along(oxygen) / max(1L, as.integer(o2_chunk_size)))
    )
    for (family in f7ft_family_levels()) {
      selected <- endpoints[
        endpoints$model_context == context & endpoints$pair_label == family,
        , drop = FALSE
      ]
      selected <- selected[order(
        selected$representative_objective_rank,
        selected$representative_seed_number
      ), , drop = FALSE]
      if (isTRUE(smoke)) selected <- selected[1L, , drop = FALSE]
      for (p_value in p_values) for (chunk_index in seq_along(chunks)) {
        index <- index + 1L
        oxygen_index <- chunks[[chunk_index]]
        rows[[index]] <- data.frame(
          task_id = sprintf("G%04d", index), model_context = context,
          pair_id = selected$pair_id[[1L]], pair_label = family,
          p_misseg = p_value, o2_chunk_index = chunk_index,
          o2_index_start = min(oxygen_index), o2_index_end = max(oxygen_index),
          endpoint_indices = paste(selected$endpoint_index, collapse = ","),
          n_unique_endpoint = nrow(selected),
          represented_optimizer_endpoint =
            sum(selected$endpoint_multiplicity_q10),
          cache_path = file.path(
            run_paths$cache, sprintf("full_range_task_%04d.rds", index)
          ), stringsAsFactors = FALSE
        )
      }
    }
  }
  tasks <- do.call(rbind, rows)
  f7ft_atomic_write_tsv(
    tasks, file.path(run_paths$run_root, "full_range_task_manifest.tsv")
  )
  tasks
}

f7g_fingerprint <- function(
    paths, endpoint_manifest, passage_bundle, propagator_path, smoke
) {
  inputs <- c(
    endpoint_manifest$paths[c("endpoints", "expanded")],
    passage_bundle$paths[c("valid", "lineage", "rule")], propagator_path,
    file.path(paths$code, "util", "analysis", "figure7_full_range_q10.R"),
    file.path(paths$oxygen_code, "util", "o2_supply_demand_map_invitro_lineage_simulation_utils.R")
  )
  paste(
    f7g_profile(), paste0("model=", f7r_model_source_fingerprint(paths)),
    paste0("inputs=", paste(unname(tools::md5sum(inputs)), collapse = ":")),
    paste0("initial=", paste(f7ft_initial_ploidy(), collapse = ",")),
    paste0("p=", paste(f7ft_p_values(), collapse = ",")),
    paste0("invivo_o2=", paste(range(f7g_o2("in vivo", smoke)), collapse = ":"),
           ":", length(f7g_o2("in vivo", smoke))),
    paste0("invitro_o2=", paste(range(f7g_o2("in vitro", smoke)), collapse = ":"),
           ":", length(f7g_o2("in vitro", smoke))),
    paste0("day=", paste(range(f7g_days(smoke)), collapse = ":"),
           ":", length(f7g_days(smoke))), sep = "|"
  )
}

f7g_compute_task <- function(
    task, endpoint_manifest, objective_bundle, contexts, paths, run_paths,
    passage_bundle, fingerprint, smoke = FALSE
) {
  f7r_load_response_engine(paths)
  indices <- as.integer(strsplit(
    task$endpoint_indices[[1L]], ",", fixed = TRUE
  )[[1L]])
  endpoints <- endpoint_manifest$endpoints[
    match(indices, endpoint_manifest$endpoints$endpoint_index), , drop = FALSE
  ]
  if (anyNA(endpoints$endpoint_index)) stop("Full-range endpoint lookup failed.")
  context <- task$model_context[[1L]]
  oxygen_all <- f7g_o2(context, smoke)
  oxygen_index <- seq.int(
    as.integer(task$o2_index_start[[1L]]),
    as.integer(task$o2_index_end[[1L]])
  )
  oxygen <- oxygen_all[oxygen_index]
  days <- f7g_days(smoke)
  p_value <- as.numeric(task$p_misseg[[1L]])
  weighted_sum <- array(
    0, dim = c(length(f7ft_initial_ploidy()), length(days), length(oxygen)),
    dimnames = list(
      initial_ploidy = paste0(f7ft_initial_ploidy(), "N"),
      day = as.character(days), O2_pct = as.character(oxygen)
    )
  )
  passage_weighted_sum <- if (context == "in vitro") weighted_sum else NULL
  passage_feasible_weight <- if (context == "in vitro") weighted_sum else NULL
  reused <- NULL
  reuse_root <- Sys.getenv("FIGURE7_REUSE_CONTINUOUS_RUN", "")
  if (nzchar(reuse_root) && !smoke) {
    reuse_path <- file.path(reuse_root, "task_cache", basename(task$cache_path[[1L]]))
    reused <- readRDS(reuse_path)
    keys <- setdiff(names(task), "cache_path")
    model_hash <- function(x) strsplit(x, "|", fixed = TRUE)[[1L]][2L]
    stopifnot(isTRUE(all.equal(task[keys], reused$task[keys], check.attributes = FALSE)),
      identical(model_hash(fingerprint), model_hash(reused$fingerprint)),
      identical(days, reused$day_values), isTRUE(all.equal(oxygen, reused$o2_values)),
      all(is.finite(reused$weighted_mean)), isTRUE(reused$qc$passed))
    # Compare endpoint parameter identities/multiplicities, not just row numbers.
    old_endpoints <- f7r_read_tsv(file.path(reuse_root, "q10_unique_endpoint_manifest.tsv"))
    endpoint_key <- function(x) paste(x$model_context, x$pair_label, x$endpoint_group)
    matched <- match(endpoint_key(endpoints), endpoint_key(old_endpoints))
    if (anyNA(matched) || !all(endpoints$endpoint_multiplicity_q10 ==
        old_endpoints$endpoint_multiplicity_q10[matched]) ||
        !all(endpoints$parameter_signature == old_endpoints$parameter_signature[matched])) {
      stop("Continuous reuse endpoint identities differ.")
    }
    if (context == "in vivo") {
      reused$profile <- f7g_profile()
      reused$fingerprint <- fingerprint
      reused$task <- task
      reused$qc$cache_path <- task$cache_path[[1L]]
      reused$qc$continuous_source <- reuse_path
      reused$continuous_source_md5 <- unname(tools::md5sum(reuse_path))
      f7ft_atomic_save_rds(reused, task$cache_path[[1L]], compress = FALSE)
      return(reused$qc)
    }
  }
  passage_rows <- list()
  passage_index <- 0L
  maximum_override_error <- 0
  maximum_formula_error <- 0
  maximum_day0_error <- 0
  for (endpoint_index in seq_len(nrow(endpoints))) {
    endpoint <- endpoints[endpoint_index, , drop = FALSE]
    prepared <- f7ft_prepare_endpoint(endpoint, objective_bundle, contexts)
    forced <- figure7_force_p_misseg(prepared$run_params, p_value)
    formula_qc <- figure7_p_misseg_formula_qc(prepared$run_params, p_value)
    maximum_override_error <- max(
      maximum_override_error, abs(as.numeric(forced$p_misseg) - p_value)
    )
    maximum_formula_error <- max(
      maximum_formula_error, formula_qc$maximum_direct_formula_error
    )
    weight <- as.integer(endpoint$endpoint_multiplicity_q10[[1L]])
    for (o2_index in seq_along(oxygen)) {
      fixed <- fixo2_fixed_matrix(
        model_env = globalenv(), cfg = prepared$config,
        run_params = forced, O2 = oxygen[[o2_index]]
      )
      n_unit <- prepared$config$N_UNIT %||% 22L
      initial_state <- f7ft_initial_matrix(
        fixed$ngrid, n_unit, f7ft_initial_ploidy()
      )
      step <- as.matrix(Matrix::expm(fixed$M))
      response <- if (is.null(reused)) f7g_propagate_operator_cpp(
        step, initial_state,
        as.numeric(fixed$ngrid) / as.numeric(n_unit), max(days),
        log(passage_bundle$rule$seed_cells[[1L]]),
        log(passage_bundle$rule$target_live_cells[[1L]]),
        FALSE
      ) else NULL
      trajectory <- if (is.null(reused)) as.matrix(response$mean_ploidy) else
        reused$weighted_mean[, , o2_index]
      maximum_day0_error <- max(
        maximum_day0_error,
        max(abs(trajectory[, 1L] - f7ft_initial_ploidy()))
      )
      if (is.null(reused)) weighted_sum[, , o2_index] <-
        weighted_sum[, , o2_index] + weight * trajectory
      if (identical(context, "in vitro")) {
        response <- f7g_propagate_segments_cpp(
          step, initial_state, as.numeric(fixed$ngrid) / as.numeric(n_unit),
          max(days), log(passage_bundle$rule$seed_cells[[1L]]),
          log(passage_bundle$rule$target_live_cells[[1L]]),
          as.integer(passage_bundle$rule$representative_passage_duration_day[[1L]])
        )
        supported <- is.finite(response$mean_ploidy)
        expected_missing <- outer(response$protocol_failure_day, days,
          function(failure, day) !is.na(failure) & day >= failure)
        stopifnot(identical(!supported, expected_missing))
        contribution <- response$mean_ploidy
        contribution[!supported] <- 0
        passage_weighted_sum[, , o2_index] <-
          passage_weighted_sum[, , o2_index] + weight * contribution
        passage_feasible_weight[, , o2_index] <-
          passage_feasible_weight[, , o2_index] + weight * supported
        maximum_day0_error <- max(maximum_day0_error,
          max(abs(response$mean_ploidy[, 1L] - f7ft_initial_ploidy())))
        passage_index <- passage_index + 1L
        passage_rows[[passage_index]] <- data.frame(
          pair_label = endpoint$pair_label[[1L]],
          endpoint_group = endpoint$endpoint_group[[1L]],
          endpoint_multiplicity_q10 = weight, p_misseg = p_value,
          O2_pct = oxygen[[o2_index]],
          initial_ploidy = f7ft_initial_ploidy(),
          passage_count = as.integer(response$passage_count),
          first_passage_day = as.integer(response$first_passage_day),
          last_passage_day = as.integer(response$last_passage_day),
          no_crossing = as.logical(response$no_crossing),
          protocol_failure_day = as.integer(response$protocol_failure_day),
          maximum_pre_post_mean_error =
            as.numeric(response$maximum_pre_post_mean_error),
          maximum_boundary_mean_jump = as.numeric(response$maximum_boundary_mean_jump),
          selected_model_day_sum = as.numeric(response$selected_model_day_sum),
          earlier_than_segment_end_count = as.integer(response$earlier_than_segment_end_count),
          selected_relative_target_distance_sum =
            as.numeric(response$selected_relative_target_distance_sum),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  represented <- sum(endpoints$endpoint_multiplicity_q10)
  weighted_mean <- if (is.null(reused)) weighted_sum / represented else reused$weighted_mean
  passage_weighted_mean <- if (context == "in vitro") passage_weighted_sum / represented else NULL
  if (context == "in vitro") {
    passage_weighted_mean[passage_feasible_weight < represented] <- NA_real_
  }
  missing_explained <- context != "in vitro" ||
    (all(is.finite(passage_weighted_mean) | is.na(passage_weighted_mean)) &&
     identical(is.na(passage_weighted_mean), passage_feasible_weight < represented))
  qc <- data.frame(
    task_id = task$task_id[[1L]], model_context = context,
    pair_label = task$pair_label[[1L]], p_misseg = p_value,
    o2_index_start = min(oxygen_index), o2_index_end = max(oxygen_index),
    n_unique_endpoint = nrow(endpoints),
    represented_optimizer_endpoint = represented,
    n_operator = nrow(endpoints) * length(oxygen),
    n_o2 = length(oxygen), n_day = length(days),
    maximum_day0_abs_error = maximum_day0_error,
    maximum_p_misseg_override_error = maximum_override_error,
    maximum_direct_formula_error = maximum_formula_error,
    maximum_passage_pre_post_mean_error = if (length(passage_rows)) {
      max(vapply(
        passage_rows, function(x) max(x$maximum_pre_post_mean_error), numeric(1L)
      ))
    } else 0,
    total_passage_events = if (length(passage_rows)) {
      sum(vapply(passage_rows, function(x) sum(x$passage_count), numeric(1L)))
    } else 0,
    all_finite = all(is.finite(weighted_mean)) && all(is.finite(passage_weighted_mean)),
    protocol_mask_valid = missing_explained,
    passed = all(is.finite(weighted_mean)) && missing_explained && maximum_day0_error <= 1e-10 &&
      maximum_override_error <= 1e-12 && maximum_formula_error <= 1e-12,
    continuous_source = if (is.null(reused)) "independent_continuous_propagation" else reuse_path,
    cache_path = task$cache_path[[1L]], stringsAsFactors = FALSE
  )
  f7ft_atomic_save_rds(list(
    profile = f7g_profile(), fingerprint = fingerprint, task = task,
    oxygen_index = oxygen_index, o2_values = oxygen, day_values = days,
    represented_optimizer_endpoint = represented,
    weighted_mean = weighted_mean,
    passage_weighted_mean = passage_weighted_mean,
    passage_feasible_weight = if (is.null(passage_feasible_weight)) NULL else {
      storage.mode(passage_feasible_weight) <- "integer"; passage_feasible_weight
    },
    continuous_source_md5 = if (is.null(reused)) NA_character_ else unname(tools::md5sum(reuse_path)),
    passage = if (length(passage_rows)) do.call(rbind, passage_rows) else NULL,
    qc = qc
  ), task$cache_path[[1L]], compress = FALSE)
  qc
}

f7g_aggregate_passage <- function(passage_rows, run_paths) {
  f7r_require_packages("data.table")
  rows <- data.table::rbindlist(passage_rows, use.names = TRUE, fill = TRUE)
  data.table::setDT(rows)
  summary <- rows[, .(
    optimizer_weight = sum(endpoint_multiplicity_q10),
    weighted_mean_passage_count = sum(
      endpoint_multiplicity_q10 * passage_count
    ) / sum(endpoint_multiplicity_q10),
    weighted_no_crossing_fraction = sum(
      endpoint_multiplicity_q10 * as.numeric(no_crossing)
    ) / sum(endpoint_multiplicity_q10),
    earliest_first_passage_day = if (all(is.na(first_passage_day))) {
      NA_integer_
    } else min(first_passage_day, na.rm = TRUE),
    latest_last_passage_day = if (all(is.na(last_passage_day))) {
      NA_integer_
    } else max(last_passage_day, na.rm = TRUE),
    maximum_pre_post_mean_error = max(maximum_pre_post_mean_error),
    maximum_boundary_mean_jump = max(maximum_boundary_mean_jump),
    weighted_mean_selected_model_day = sum(endpoint_multiplicity_q10 * selected_model_day_sum) /
      sum(endpoint_multiplicity_q10 * passage_count),
    weighted_earlier_selection_fraction = sum(endpoint_multiplicity_q10 * earlier_than_segment_end_count) /
      sum(endpoint_multiplicity_q10 * passage_count),
    weighted_protocol_infeasible_fraction = sum(endpoint_multiplicity_q10 * !is.na(protocol_failure_day)) /
      sum(endpoint_multiplicity_q10),
    earliest_protocol_failure_day = if (all(is.na(protocol_failure_day))) NA_integer_ else
      min(protocol_failure_day, na.rm = TRUE)
  ), by = .(pair_label, p_misseg, O2_pct, initial_ploidy)]
  full_path <- f7ft_atomic_save_rds(
    as.data.frame(rows),
    file.path(run_paths$run_root, "canonical_passage_endpoint_summary.rds"),
    compress = "gzip"
  )
  summary_path <- f7ft_atomic_write_tsv(
    as.data.frame(summary),
    file.path(run_paths$run_root, "canonical_passage_grid_summary.tsv")
  )
  list(full = full_path, summary = summary_path, data = summary)
}

f7g_aggregate <- function(
    tasks, run_paths, fingerprint, passage_bundle, smoke = FALSE
) {
  task_qc <- as.data.frame(data.table::rbindlist(
    lapply(tasks$cache_path, function(path) readRDS(path)$qc), fill = TRUE))
  qc_path <- f7ft_atomic_write_tsv(
    task_qc, file.path(run_paths$run_root, "full_range_task_validation.tsv")
  )
  if (!all(task_qc$passed)) stop("One or more full-range tasks failed validation.")

  panel_validation <- list()
  panel_index <- 0L
  passage_rows <- list()
  passage_row_index <- 0L
  for (context in f7g_contexts()) for (family in f7ft_family_levels()) {
    oxygen <- f7g_o2(context, smoke)
    days <- f7g_days(smoke)
    p_values <- if (isTRUE(smoke)) f7ft_p_values()[1:2] else f7ft_p_values()
    values <- array(
      NA_real_, dim = c(
        length(f7ft_initial_ploidy()), length(days), length(oxygen),
        length(p_values)
      ), dimnames = list(
        initial_ploidy = paste0(f7ft_initial_ploidy(), "N"),
        day = as.character(days), O2_pct = as.character(oxygen),
        p_misseg = format(p_values, scientific = FALSE, trim = TRUE)
      )
    )
    selected_tasks <- tasks[
      tasks$model_context == context & tasks$pair_label == family,
      , drop = FALSE
    ]
    passage_values <- if (context == "in vitro") values else NULL
    feasible_weight <- if (context == "in vitro") {
      weights <- values; storage.mode(weights) <- "integer"; weights
    } else NULL
    for (task_index in seq_len(nrow(selected_tasks))) {
      task <- selected_tasks[task_index, , drop = FALSE]
      cached <- readRDS(task$cache_path[[1L]])
      if (!identical(cached$fingerprint, fingerprint)) {
        stop("Full-range task fingerprint mismatch.")
      }
      p_index <- match(as.numeric(task$p_misseg[[1L]]), p_values)
      values[, , cached$oxygen_index, p_index] <- cached$weighted_mean
      if (context == "in vitro") {
        if (is.null(cached$passage_weighted_mean)) stop("Missing independently computed passage data.")
        passage_values[, , cached$oxygen_index, p_index] <- cached$passage_weighted_mean
        feasible_weight[, , cached$oxygen_index, p_index] <- cached$passage_feasible_weight
      }
      if (!is.null(cached$passage)) {
        passage_row_index <- passage_row_index + 1L
        passage_rows[[passage_row_index]] <- cached$passage
      }
    }
    if (any(!is.finite(values))) {
      stop("Incomplete full-range aggregate for ", context, " ", family)
    }
    for (mode in f7g_modes(context)) {
      mode_values <- if (mode == "passage") passage_values else values
      if (mode == "continuous" && any(!is.finite(mode_values))) stop("Incomplete continuous aggregate.")
      represented <- if (isTRUE(smoke)) unique(selected_tasks$represented_optimizer_endpoint)[[1L]] else 50L
      mask_valid <- mode != "passage" || (!anyNA(feasible_weight) &&
        all(feasible_weight >= 0 & feasible_weight <= represented) &&
        identical(is.na(mode_values), feasible_weight < represented))
      if (!mask_valid) stop("Unexplained missing passage values.")
      object <- list(
        profile = f7g_profile(), fingerprint = fingerprint,
        pair_label = family,
        panel_letter = if (context == "in vivo") {
          if (family == "C01") "C" else "D"
        } else {
          if (family == "C01") "E" else "F"
        },
        model_context = context, propagation_mode = mode,
        initial_ploidy = f7ft_initial_ploidy(), day_values = days,
        o2_values = oxygen, p_misseg = p_values,
        optimizer_endpoint_weight = rep(
          if (isTRUE(smoke)) unique(selected_tasks$represented_optimizer_endpoint)[[1L]] else 50L,
          length(p_values)
        ),
        canonical_passage_rule = if (mode == "passage") passage_bundle$rule else NULL,
        protocol_feasible_optimizer_weight = if (mode == "passage") feasible_weight else NULL,
        mean_ploidy = mode_values
      )
      path <- file.path(
        run_paths$run_root, f7g_panel_filename(context, mode, family)
      )
      f7ft_atomic_save_rds(object, path, compress = "gzip")
      day0 <- mode_values[, 1L, , , drop = FALSE]
      expected <- array(
        f7ft_initial_ploidy(), dim = dim(day0), dimnames = dimnames(day0)
      )
      panel_index <- panel_index + 1L
      panel_validation[[panel_index]] <- data.frame(
        model_context = context, propagation_mode = mode,
        pair_label = family, n_initial_ploidy = length(object$initial_ploidy),
        n_day = length(days), n_o2 = length(oxygen),
        n_p_misseg = length(p_values),
        maximum_day0_abs_error = max(abs(day0 - expected)),
        minimum_mean_ploidy = min(mode_values, na.rm = TRUE), maximum_mean_ploidy = max(mode_values, na.rm = TRUE),
        all_finite = all(is.finite(mode_values)),
        fraction_protocol_infeasible = mean(is.na(mode_values)),
        protocol_mask_valid = mask_valid,
        panel_path = normalizePath(path, mustWork = TRUE),
        passed = max(abs(day0 - expected)) <= 1e-10 &&
          min(mode_values, na.rm = TRUE) >= 1 - 1e-8 && max(mode_values, na.rm = TRUE) <= 7 + 1e-8 && mask_valid,
        stringsAsFactors = FALSE
      )
    }
  }
  validation <- do.call(rbind, panel_validation)
  expected_modes <- data.frame(
    model_context = c("in vivo", "in vitro", "in vitro"),
    propagation_mode = c("continuous", "continuous", "passage")
  )
  observed_modes <- unique(validation[c("model_context", "propagation_mode")])
  mode_ok <- nrow(observed_modes) == nrow(expected_modes) &&
    all(paste(expected_modes$model_context, expected_modes$propagation_mode) %in%
      paste(observed_modes$model_context, observed_modes$propagation_mode)) &&
    !any(validation$model_context == "in vivo" &
      validation$propagation_mode == "passage")
  validation$mode_contract_passed <- mode_ok
  panel_path <- f7ft_atomic_write_tsv(
    validation, file.path(run_paths$run_root, "full_range_panel_validation.tsv")
  )
  if (!all(validation$passed) || !mode_ok || nrow(validation) != 6L) {
    stop("Figure 7 full-range panel validation failed.")
  }
  passage <- f7g_aggregate_passage(passage_rows, run_paths)
  equivalence <- do.call(rbind, lapply(f7ft_family_levels(), function(family) {
    continuous <- f7g_read_panel(run_paths, "in vitro", "continuous", family)
    segmented <- f7g_read_panel(run_paths, "in vitro", "passage", family)
    delta <- abs(segmented$mean_ploidy - continuous$mean_ploidy)
    data.frame(comparison = "invitro_passage_vs_continuous_mean_ploidy",
      pair_label = family, maximum_absolute_delta = max(delta, na.rm = TRUE),
      mean_absolute_delta = mean(delta, na.rm = TRUE), fraction_different_gt_1e_10 = mean(delta > 1e-10, na.rm = TRUE),
      fraction_protocol_infeasible = mean(is.na(segmented$mean_ploidy)),
      reason = "measured_difference_of_independent_propagation_modes_no_required_sign_or_size",
      passed = all(is.finite(delta) | is.na(segmented$mean_ploidy)) &&
        identical(is.na(segmented$mean_ploidy), segmented$protocol_feasible_optimizer_weight <
          segmented$optimizer_endpoint_weight[[1L]]),
      stringsAsFactors = FALSE)
  }))
  equivalence_path <- f7ft_atomic_write_tsv(
    equivalence,
    file.path(run_paths$run_root, "passage_vs_continuous_validation.tsv")
  )
  list(
    task_qc = qc_path, panels = panel_path, passage = passage,
    equivalence = equivalence_path
  )
}

f7g_publish_current <- function(paths, run_paths, fingerprint) {
  pointer <- data.frame(
    run_id = run_paths$run_id,
    relative_run_path = file.path("finite_time_full_q10_runs", run_paths$run_id),
    profile = f7g_profile(), fingerprint = fingerprint,
    published_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  f7ft_atomic_write_tsv(pointer, run_paths$current)
}

f7g_chart_contract <- function(run_paths) {
  contract <- data.frame(
    artifact = c(
      "Figure 7B", "Supplementary Figure 7-10",
      "Supplementary Figure 7-11", "Supplementary Figure 7-12"
    ),
    model_context = c("in vivo and in vitro", "in vivo", "in vitro", "in vitro"),
    propagation_mode = c(
      "continuous in vivo; passage in vitro", "continuous",
      "continuous", "passage"
    ),
    displayed_time_day = c("0:1000", "0:10000", "0:10000", "0:10000"),
    displayed_oxygen_pct = c("0:2", "0:5", "0:20", "0:20"),
    displayed_initial_ploidy = c("2,4", "2,3,4,5,6", "2,3,4,5,6", "2,3,4,5,6"),
    data_source = "one full-range q10 run; drawing-only subsets",
    primary_encoding = paste0(
      "experimental day on x; fixed oxygen on y; arithmetic-mean finite-time ",
      "ploidy as log-scaled color; columns are fixed p_misseg"
    ), stringsAsFactors = FALSE
  )
  f7ft_atomic_write_tsv(
    contract, file.path(run_paths$run_root, "full_range_chart_contract.tsv")
  )
}

f7g_data <- function(
    workspace_root = f7r_find_workspace_root(), n_core = 1L,
    run_id = f7ft_resolve_run_id(), smoke = FALSE,
    publish_current = !isTRUE(smoke), o2_chunk_size = 10L
) {
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  f7r_require_packages(c("Matrix", "Rcpp", "data.table", "future", "future.apply"))
  paths <- f7r_paths(workspace_root)
  run_paths <- f7g_paths(paths, run_id = run_id, create = TRUE)
  if (length(list.files(run_paths$cache, all.files = TRUE, no.. = TRUE))) {
    stop("Fresh Figure 7 full-range cache is not empty: ", run_paths$cache)
  }
  f7r_load_response_engine(paths)
  propagator <- f7g_load_propagator(paths)
  objective_bundle <- f7r_objective_selection(paths)
  endpoint_manifest <- f7ft_build_endpoint_manifest(
    paths, objective_bundle, run_paths
  )
  f7g_chart_contract(run_paths)
  schedule <- f7p_extract_schedule(run_paths, paths)
  passage_bundle <- f7g_canonical_passage_rule(schedule, run_paths)
  fingerprint <- f7g_fingerprint(
    paths, endpoint_manifest, passage_bundle, propagator, smoke
  )
  contexts <- lapply(
    unique(endpoint_manifest$endpoints$pair_id),
    f7r_pair_model_context, selected = objective_bundle$selected, paths = paths
  )
  names(contexts) <- unique(endpoint_manifest$endpoints$pair_id)
  source(file.path(paths$code, "util", "analysis", "figure7_segment_validation.R"))
  validation_cases <- list()
  for (family in f7ft_family_levels()) {
    endpoint <- endpoint_manifest$endpoints[
      endpoint_manifest$endpoints$model_context == "in vitro" &
        endpoint_manifest$endpoints$pair_label == family, , drop = FALSE][1L, ]
    prepared <- f7ft_prepare_endpoint(endpoint, objective_bundle, contexts)
    for (oxygen in c(0, 0.5, 2, 20)) for (p in c(0.005, 0.3)) {
      fixed <- fixo2_fixed_matrix(globalenv(), prepared$config,
        figure7_force_p_misseg(prepared$run_params, p), O2 = oxygen)
      unit <- prepared$config$N_UNIT %||% 22L
      validation_cases[[length(validation_cases) + 1L]] <- list(
        name = paste(family, oxygen, p, sep = "_"),
        step = as.matrix(Matrix::expm(fixed$M)),
        initial = f7ft_initial_matrix(fixed$ngrid, unit, f7ft_initial_ploidy()),
        ploidy = as.numeric(fixed$ngrid) / unit, days = 25L,
        seed = passage_bundle$rule$seed_cells[[1L]],
        target = passage_bundle$rule$target_live_cells[[1L]],
        duration = as.integer(passage_bundle$rule$representative_passage_duration_day[[1L]]))
    }
  }
  f7g_validate_segments(paths, run_paths, validation_cases)
  tasks <- f7g_task_manifest(
    endpoint_manifest, run_paths, smoke = smoke,
    o2_chunk_size = o2_chunk_size
  )
  execution_order <- order(
    tasks$o2_chunk_index, tasks$p_misseg,
    tasks$model_context, tasks$pair_label
  )
  task_list <- split(tasks[execution_order, , drop = FALSE], seq_len(nrow(tasks)))
  compute_one <- function(task) tryCatch(
    f7g_compute_task(
      task, endpoint_manifest, objective_bundle, contexts, paths, run_paths,
      passage_bundle, fingerprint, smoke = smoke
    ), error = function(e) structure(
      list(task_id = task$task_id[[1L]], message = conditionMessage(e)),
      class = "f7g_error"
    )
  )
  message(
    "Figure 7 full-range q10: ", nrow(tasks), " tasks, workers=",
    min(as.integer(n_core), nrow(tasks)), ", run_id=", run_paths$run_id
  )
  results <- f7ft_parallel_lapply(task_list, compute_one, n_core = n_core)
  failed <- vapply(results, inherits, logical(1L), "f7g_error")
  if (any(failed)) stop(
    "Figure 7 full-range task failures: ",
    paste(vapply(results[failed], `[[`, character(1L), "message"), collapse = "; ")
  )
  aggregation <- f7g_aggregate(
    tasks, run_paths, fingerprint, passage_bundle, smoke = smoke
  )
  if (!identical(fingerprint, f7g_fingerprint(paths, endpoint_manifest, passage_bundle, propagator, smoke))) {
    stop("Model, selector, or computation source changed while the run was active.")
  }
  provenance <- data.frame(
    key = c(
      "profile", "run_id", "fingerprint", "workspace_root",
      "model_code_root", "model_source_fingerprint", "joint_result_root",
      "n_core", "smoke", "optimizer_endpoint_per_panel", "day_grid",
      "invivo_oxygen_grid", "invitro_oxygen_grid", "invivo_modes",
      "invitro_modes", "endpoint_aggregation", "passage_source_fit_result"
    ), value = c(
      f7g_profile(), run_paths$run_id, fingerprint,
      normalizePath(paths$root, mustWork = TRUE),
      normalizePath(paths$oxygen_code, mustWork = TRUE),
      f7r_model_source_fingerprint(paths), f7p_joint_result_root(),
      as.character(n_core), as.character(smoke),
      if (isTRUE(smoke)) "one unique endpoint per panel" else "50",
      paste(range(f7g_days(smoke)), collapse = ":"),
      paste0(paste(range(f7g_o2("in vivo", smoke)), collapse = ":"),
             ":", length(f7g_o2("in vivo", smoke))),
      paste0(paste(range(f7g_o2("in vitro", smoke)), collapse = ":"),
             ":", length(f7g_o2("in vitro", smoke))),
      paste(f7g_modes("in vivo"), collapse = ","),
      paste(f7g_modes("in vitro"), collapse = ","),
      "arithmetic mean with optimizer-endpoint multiplicity restored",
      passage_bundle$source_fit_result
    ), stringsAsFactors = FALSE
  )
  provenance_path <- f7ft_atomic_write_tsv(
    provenance, file.path(run_paths$run_root, "full_range_run_provenance.tsv")
  )
  if (isTRUE(publish_current)) f7g_publish_current(paths, run_paths, fingerprint)
  invisible(list(
    paths = run_paths, endpoint_manifest = endpoint_manifest,
    passage = passage_bundle, tasks = tasks, aggregation = aggregation,
    provenance = provenance_path
  ))
}
