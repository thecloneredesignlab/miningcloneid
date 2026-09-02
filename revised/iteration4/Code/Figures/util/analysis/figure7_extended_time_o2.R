#!/usr/bin/env Rscript

# Extended 0-20% oxygen, 0-1000 day in-vitro q10 response surfaces for
# Supplementary Figures 7-11 (repeated fitted passaging) and 7-12
# (continuous culture without passaging).

options(stringsAsFactors = FALSE, warn = 1)

f7e_modes <- function() c("passage", "continuous")

f7e_profile <- function(mode) {
  mode <- match.arg(mode, f7e_modes())
  paste0("figure7_fixed_pmisseg_invitro_", mode, "_1000d_o2_0_20_q10_v1")
}

f7e_o2_values <- function(smoke = FALSE) {
  if (isTRUE(smoke)) c(0, 0.5, 20) else seq(0, 20, by = 0.1)
}

f7e_day_values <- function(smoke = FALSE) {
  if (isTRUE(smoke)) 0:10 else 0:1000
}

f7e_paths <- function(paths, mode, run_id = NULL, create = FALSE) {
  mode <- match.arg(mode, f7e_modes())
  base <- file.path(paths$figure7, paste0("invitro_", mode, "_extended_q10_runs"))
  current <- file.path(
    paths$figure7, paste0("invitro_", mode, "_extended_q10_current.tsv")
  )
  if (is.null(run_id)) {
    f7r_require_files(current, paste0("current Figure 7 ", mode, " extended run pointer"))
    pointer <- f7r_read_tsv(current)
    required <- c("run_id", "relative_run_path", "profile", "fingerprint")
    if (nrow(pointer) != 1L || !all(required %in% names(pointer))) {
      stop("Malformed Figure 7 extended run pointer: ", current)
    }
    run_id <- f7ft_sanitize_run_id(pointer$run_id[[1L]])
    run_root <- normalizePath(
      file.path(paths$figure7, pointer$relative_run_path[[1L]]), mustWork = TRUE
    )
  } else {
    run_id <- f7ft_sanitize_run_id(run_id)
    run_root <- file.path(base, run_id)
  }
  out <- list(
    mode = mode, run_id = run_id, base = base, run_root = run_root,
    cache = file.path(run_root, "task_cache"), current = current
  )
  if (isTRUE(create)) {
    invisible(lapply(
      unname(unlist(out[c("base", "run_root", "cache")])),
      dir.create, recursive = TRUE, showWarnings = FALSE
    ))
  }
  out
}

f7e_read_panel_object <- function(run_paths, panel_letter, mode = run_paths$mode) {
  mode <- match.arg(mode, f7e_modes())
  path <- file.path(
    run_paths$run_root,
    paste0("extended_time_panel_", tolower(panel_letter), ".rds")
  )
  f7r_require_files(path, paste0("Figure 7 ", mode, " extended panel"))
  object <- readRDS(path)
  if (!identical(object$profile, f7e_profile(mode)) ||
      !identical(object$panel_letter, panel_letter)) {
    stop("Unexpected Figure 7 extended panel object: ", path)
  }
  object
}

f7e_expand_schedule <- function(schedule_bundle, max_day = 1000L) {
  max_day <- as.integer(max_day)
  if (!is.finite(max_day) || max_day < 1L) {
    stop("Extended passage max_day must be a positive integer.")
  }
  cycle_days <- as.integer(schedule_bundle$max_experimental_day) + 1L
  if (!is.finite(cycle_days) || cycle_days < 1L) {
    stop("Recorded passage schedule has no positive experimental duration.")
  }
  rows <- list()
  row_index <- 0L
  for (lineage in unique(schedule_bundle$schedule$lineage)) {
    template <- schedule_bundle$schedule[
      schedule_bundle$schedule$lineage == lineage, , drop = FALSE
    ]
    template <- template[order(template$passage_index), , drop = FALSE]
    passage_index <- 0L
    for (cycle_index in seq_len(ceiling((max_day + 1L) / cycle_days))) {
      cycle_start <- (cycle_index - 1L) * cycle_days
      for (template_index in seq_len(nrow(template))) {
        segment <- template[template_index, , drop = FALSE]
        segment_start <- cycle_start +
          as.integer(round(segment$experimental_day_start[[1L]]))
        if (segment_start >= max_day) next
        duration <- as.integer(round(segment$passage_duration[[1L]]))
        passage_index <- passage_index + 1L
        row_index <- row_index + 1L
        segment$cycle_index <- cycle_index
        segment$template_passage_index <- as.integer(segment$passage_index[[1L]])
        segment$passage_index <- passage_index
        segment$experimental_day_start <- segment_start
        segment$experimental_day_end <- segment_start + duration
        rows[[row_index]] <- segment
      }
    }
  }
  expanded <- do.call(rbind, rows)
  rownames(expanded) <- NULL
  if (length(unique(expanded$lineage)) != 6L ||
      any(expanded$passage_duration <= 0)) {
    stop("Repeated passage schedule failed coverage validation.")
  }
  list(
    schedule = expanded,
    source_schedule = schedule_bundle$schedule,
    fit_result_path = schedule_bundle$fit_result_path,
    schedule_path = NULL,
    source_schedule_path = NULL,
    max_experimental_day = max_day,
    repeat_policy = paste0(
      "repeat the complete six-lineage experiment in ", cycle_days,
      "-day blocks; preserve inactive intervals and carry terminal composition"
    )
  )
}

f7e_schedule_coverage <- function(schedule_bundle) {
  schedule <- schedule_bundle$schedule
  lineage_table <- unique(schedule[c("cohort", "lineage")])
  do.call(rbind, lapply(0:schedule_bundle$max_experimental_day, function(day) {
    active <- vapply(lineage_table$lineage, function(lineage) {
      rows <- schedule[schedule$lineage == lineage, , drop = FALSE]
      any(rows$experimental_day_start <= day & rows$experimental_day_end >= day)
    }, logical(1L))
    data.frame(
      day = day,
      active_lineages_2N = sum(active & lineage_table$cohort == "2N"),
      active_lineages_4N = sum(active & lineage_table$cohort == "4N"),
      active_lineages_total = sum(active),
      active_cohorts = sum(vapply(c("2N", "4N"), function(cohort) {
        any(active & lineage_table$cohort == cohort)
      }, logical(1L))),
      stringsAsFactors = FALSE
    )
  }))
}

f7e_task_manifest <- function(
    endpoints, run_paths, smoke = FALSE, chunk_size = 2L
) {
  endpoints <- endpoints[endpoints$model_context == "in vitro", , drop = FALSE]
  if (isTRUE(smoke)) {
    endpoints <- do.call(rbind, lapply(
      split(endpoints, endpoints$pair_label), function(x) x[1L, , drop = FALSE]
    ))
    p_values <- f7ft_p_values()[1:2]
  } else {
    p_values <- f7ft_p_values()
  }
  rows <- list()
  index <- 0L
  for (family in f7ft_family_levels()) {
    family_endpoints <- endpoints[endpoints$pair_label == family, , drop = FALSE]
    family_endpoints <- family_endpoints[order(
      family_endpoints$representative_objective_rank,
      family_endpoints$representative_seed_number
    ), , drop = FALSE]
    chunks <- split(
      seq_len(nrow(family_endpoints)),
      ceiling(seq_len(nrow(family_endpoints)) / chunk_size)
    )
    for (p_value in p_values) for (chunk_index in seq_along(chunks)) {
      selected <- family_endpoints[chunks[[chunk_index]], , drop = FALSE]
      index <- index + 1L
      rows[[index]] <- data.frame(
        task_id = sprintf("X%04d", index), mode = run_paths$mode,
        pair_id = family_endpoints$pair_id[[1L]], pair_label = family,
        panel_letter = if (family == "C01") "E" else "F",
        p_misseg = p_value, chunk_index = chunk_index,
        endpoint_indices = paste(selected$endpoint_index, collapse = ","),
        n_unique_endpoint = nrow(selected),
        represented_optimizer_endpoint = sum(selected$endpoint_multiplicity_q10),
        cache_path = file.path(
          run_paths$cache, sprintf("extended_task_%04d.rds", index)
        ), stringsAsFactors = FALSE
      )
    }
  }
  tasks <- do.call(rbind, rows)
  f7ft_atomic_write_tsv(
    tasks, file.path(run_paths$run_root, "extended_task_manifest.tsv")
  )
  tasks
}

f7e_lineage_trajectory <- function(
    step, ploidy_grid, ngrid, n_unit, initial_ploidy, lineage_schedule,
    default_initial_cells, max_day
) {
  values <- rep(NA_real_, max_day + 1L)
  state_fraction <- rep(0, length(ngrid))
  initial_index <- match(initial_ploidy * n_unit, ngrid)
  if (is.na(initial_index)) stop("Initial ploidy is outside the model grid.")
  state_fraction[[initial_index]] <- 1
  selected_days <- integer()
  current_day <- 0L

  for (passage_index in seq_len(nrow(lineage_schedule))) {
    if (current_day >= max_day) break
    segment <- lineage_schedule[passage_index, , drop = FALSE]
    segment_start <- as.integer(round(segment$experimental_day_start[[1L]]))
    if (segment_start > current_day) current_day <- segment_start
    if (current_day >= max_day) break
    recorded_duration <- as.integer(round(segment$passage_duration[[1L]]))
    duration <- min(recorded_duration, max_day - current_day)
    initial_cells <- as.numeric(segment$initial_cells[[1L]])
    if (!is.finite(initial_cells) || initial_cells <= 0) {
      initial_cells <- default_initial_cells
    }
    state <- state_fraction * initial_cells
    states <- matrix(NA_real_, nrow = length(ngrid), ncol = duration + 1L)
    states[, 1L] <- state
    if (duration > 0L) for (day_index in seq_len(duration)) {
      state <- f7p_clean_state(step %*% state)
      if (anyNA(state)) stop("Non-finite repeated-passage state.")
      states[, day_index + 1L] <- state
    }
    means <- apply(states, 2L, f7p_state_mean, ploidy_grid = ploidy_grid)
    indices <- current_day + 0:duration + 1L
    values[indices] <- means
    if (duration < recorded_duration) {
      current_day <- current_day + duration
      break
    }
    live_cells <- colSums(states)
    target <- as.numeric(segment$target_live_cells[[1L]])
    candidates <- 2:ncol(states)
    selected_column <- if (is.finite(target) && target > 0) {
      candidates[[order(abs(live_cells[candidates] - target), candidates)[[1L]]]]
    } else {
      ncol(states)
    }
    selected_days <- c(selected_days, selected_column - 1L)
    selected_state <- states[, selected_column]
    selected_total <- sum(selected_state)
    if (!is.finite(selected_total) || selected_total <= 0) {
      stop("Repeated passage selection returned an empty state.")
    }
    state_fraction <- selected_state / selected_total
    values[[current_day + duration + 1L]] <- f7p_state_mean(
      state_fraction, ploidy_grid
    )
    current_day <- current_day + duration
  }
  list(values = values, selected_days = selected_days)
}

f7e_passage_operator_response <- function(
    M, ngrid, n_unit, schedule_bundle, default_initial_cells, day_values
) {
  if (!identical(as.integer(day_values), 0:max(day_values))) {
    stop("Extended passage calculation requires a complete daily grid.")
  }
  step <- as.matrix(Matrix::expm(M))
  ploidy_grid <- as.numeric(ngrid) / as.numeric(n_unit)
  response <- matrix(
    NA_real_, nrow = length(f7ft_initial_ploidy()), ncol = length(day_values),
    dimnames = list(
      paste0(f7ft_initial_ploidy(), "N"), as.character(day_values)
    )
  )
  selection <- c(
    n_selection = 0, n_selected_before_duration = 0,
    minimum_selected_day = Inf, maximum_selected_day = -Inf
  )
  for (initial_index in seq_along(f7ft_initial_ploidy())) {
    cohort_trajectories <- list()
    for (cohort in c("2N", "4N")) {
      lineages <- unique(schedule_bundle$schedule$lineage[
        schedule_bundle$schedule$cohort == cohort
      ])
      lineage_values <- matrix(
        NA_real_, nrow = length(lineages), ncol = length(day_values)
      )
      for (lineage_index in seq_along(lineages)) {
        lineage_schedule <- schedule_bundle$schedule[
          schedule_bundle$schedule$lineage == lineages[[lineage_index]],
          , drop = FALSE
        ]
        lineage_schedule <- lineage_schedule[
          order(lineage_schedule$passage_index), , drop = FALSE
        ]
        trajectory <- f7e_lineage_trajectory(
          step, ploidy_grid, ngrid, n_unit,
          f7ft_initial_ploidy()[[initial_index]], lineage_schedule,
          default_initial_cells, max(day_values)
        )
        lineage_values[lineage_index, ] <- trajectory$values
        completed <- lineage_schedule$experimental_day_end <= max(day_values)
        completed_duration <- as.integer(round(
          lineage_schedule$passage_duration[completed]
        ))
        selection[["n_selection"]] <- selection[["n_selection"]] +
          length(trajectory$selected_days)
        selection[["n_selected_before_duration"]] <-
          selection[["n_selected_before_duration"]] +
          sum(trajectory$selected_days < completed_duration)
        if (length(trajectory$selected_days)) {
          selection[["minimum_selected_day"]] <- min(
            selection[["minimum_selected_day"]], trajectory$selected_days
          )
          selection[["maximum_selected_day"]] <- max(
            selection[["maximum_selected_day"]], trajectory$selected_days
          )
        }
      }
      cohort_trajectories[[cohort]] <- apply(
        lineage_values, 2L,
        function(x) if (all(!is.finite(x))) NA_real_ else mean(x[is.finite(x)])
      )
    }
    response[initial_index, ] <- apply(
      do.call(rbind, cohort_trajectories), 2L,
      function(x) if (all(!is.finite(x))) NA_real_ else mean(x[is.finite(x)])
    )
  }
  list(mean_ploidy = response, selection = selection)
}

f7e_compute_task <- function(
    task, endpoints, objective_bundle, contexts, paths, run_paths,
    schedule_bundle, fingerprint, smoke = FALSE
) {
  f7r_load_response_engine(paths)
  indices <- as.integer(strsplit(
    task$endpoint_indices[[1L]], ",", fixed = TRUE
  )[[1L]])
  selected <- endpoints[match(indices, endpoints$endpoint_index), , drop = FALSE]
  if (anyNA(selected$endpoint_index)) stop("Extended endpoint lookup failed.")
  o2_values <- f7e_o2_values(smoke)
  day_values <- f7e_day_values(smoke)
  weighted_sum <- array(
    0,
    dim = c(length(f7ft_initial_ploidy()), length(day_values), length(o2_values)),
    dimnames = list(
      initial_ploidy = paste0(f7ft_initial_ploidy(), "N"),
      day = as.character(day_values), O2_pct = as.character(o2_values)
    )
  )
  p_value <- as.numeric(task$p_misseg[[1L]])
  operator_count <- 0L
  maximum_override_error <- 0
  maximum_formula_error <- 0
  selection_total <- c(
    n_selection = 0, n_selected_before_duration = 0,
    minimum_selected_day = Inf, maximum_selected_day = -Inf
  )
  for (endpoint_index in seq_len(nrow(selected))) {
    endpoint <- selected[endpoint_index, , drop = FALSE]
    prepared <- f7ft_prepare_endpoint(endpoint, objective_bundle, contexts)
    forced <- figure7_force_p_misseg(prepared$run_params, p_value)
    formula_qc <- figure7_p_misseg_formula_qc(
      prepared$run_params, p_value
    )
    maximum_override_error <- max(
      maximum_override_error, abs(as.numeric(forced$p_misseg) - p_value)
    )
    maximum_formula_error <- max(
      maximum_formula_error, formula_qc$maximum_direct_formula_error
    )
    weight <- as.integer(endpoint$endpoint_multiplicity_q10[[1L]])
    for (o2_index in seq_along(o2_values)) {
      fixed <- fixo2_fixed_matrix(
        model_env = globalenv(), cfg = prepared$config,
        run_params = forced, O2 = o2_values[[o2_index]]
      )
      if (run_paths$mode == "passage") {
        result <- f7e_passage_operator_response(
          fixed$M, fixed$ngrid, prepared$config$N_UNIT %||% 22L,
          schedule_bundle, prepared$config$init_total_size %||% 1e6,
          day_values
        )
        trajectory <- result$mean_ploidy
        selection_total[c("n_selection", "n_selected_before_duration")] <-
          selection_total[c("n_selection", "n_selected_before_duration")] +
          result$selection[c("n_selection", "n_selected_before_duration")]
        if (is.finite(result$selection[["minimum_selected_day"]])) {
          selection_total[["minimum_selected_day"]] <- min(
            selection_total[["minimum_selected_day"]],
            result$selection[["minimum_selected_day"]]
          )
          selection_total[["maximum_selected_day"]] <- max(
            selection_total[["maximum_selected_day"]],
            result$selection[["maximum_selected_day"]]
          )
        }
      } else {
        trajectory <- f7ft_expm_daily_mean_ploidy(
          fixed$M, fixed$ngrid, prepared$config$N_UNIT %||% 22L,
          day_values
        )
      }
      weighted_sum[, , o2_index] <- weighted_sum[, , o2_index] +
        weight * trajectory
      operator_count <- operator_count + 1L
    }
  }
  represented <- sum(selected$endpoint_multiplicity_q10)
  day0 <- weighted_sum[, 1L, , drop = FALSE] / represented
  expected <- array(
    f7ft_initial_ploidy(), dim = dim(day0), dimnames = dimnames(day0)
  )
  qc <- data.frame(
    task_id = task$task_id[[1L]], mode = run_paths$mode,
    pair_label = task$pair_label[[1L]], p_misseg = p_value,
    n_unique_endpoint = nrow(selected),
    represented_optimizer_endpoint = represented,
    n_operator = operator_count, n_o2 = length(o2_values),
    n_day = length(day_values),
    n_passage_selection = selection_total[["n_selection"]],
    n_selected_before_duration = selection_total[["n_selected_before_duration"]],
    minimum_selected_day = selection_total[["minimum_selected_day"]],
    maximum_selected_day = selection_total[["maximum_selected_day"]],
    maximum_day0_abs_error = max(abs(day0 - expected), na.rm = TRUE),
    maximum_p_misseg_override_error = maximum_override_error,
    maximum_direct_formula_error = maximum_formula_error,
    all_finite = all(is.finite(weighted_sum)),
    passed = all(is.finite(weighted_sum)) &&
      max(abs(day0 - expected), na.rm = TRUE) <= 1e-10 &&
      maximum_override_error <= 1e-12 && maximum_formula_error <= 1e-12,
    cache_path = task$cache_path[[1L]], stringsAsFactors = FALSE
  )
  f7ft_atomic_save_rds(list(
    profile = f7e_profile(run_paths$mode), fingerprint = fingerprint,
    task = task, endpoints = selected, o2_values = o2_values,
    day_values = day_values, weighted_sum = weighted_sum, qc = qc
  ), task$cache_path[[1L]], compress = FALSE)
  qc
}

f7e_aggregate <- function(
    tasks, run_paths, fingerprint, smoke = FALSE
) {
  objects <- lapply(tasks$cache_path, readRDS)
  if (any(!vapply(
    objects, function(x) identical(x$fingerprint, fingerprint), logical(1L)
  ))) stop("Extended task caches do not share one fingerprint.")
  task_qc <- do.call(rbind, lapply(objects, `[[`, "qc"))
  f7ft_atomic_write_tsv(
    task_qc, file.path(run_paths$run_root, "extended_task_validation.tsv")
  )
  if (!all(task_qc$passed)) stop("One or more extended tasks failed validation.")

  panel_objects <- list()
  for (family in f7ft_family_levels()) {
    letter <- if (family == "C01") "E" else "F"
    p_values <- sort(unique(tasks$p_misseg[tasks$pair_label == family]))
    arrays <- vector("list", length(p_values))
    weights <- integer(length(p_values))
    for (p_index in seq_along(p_values)) {
      keep <- which(
        tasks$pair_label == family &
          abs(tasks$p_misseg - p_values[[p_index]]) < 1e-12
      )
      selected <- objects[keep]
      arrays[[p_index]] <- Reduce(`+`, lapply(selected, `[[`, "weighted_sum"))
      weights[[p_index]] <- sum(vapply(
        selected,
        function(x) sum(x$endpoints$endpoint_multiplicity_q10), integer(1L)
      ))
      expected_weight <- if (isTRUE(smoke)) weights[[p_index]] else 50L
      if (weights[[p_index]] != expected_weight) {
        stop("Unexpected extended optimizer weight for panel ", letter)
      }
      arrays[[p_index]] <- arrays[[p_index]] / weights[[p_index]]
    }
    values <- simplify2array(arrays)
    dimnames(values)[[4L]] <- format(p_values, scientific = FALSE, trim = TRUE)
    object <- list(
      profile = f7e_profile(run_paths$mode), fingerprint = fingerprint,
      panel_letter = letter, pair_label = family,
      model_context = if (run_paths$mode == "passage") {
        "in vitro with repeated passaging"
      } else {
        "in vitro, no passaging"
      },
      propagation_mode = if (run_paths$mode == "passage") {
        "passage_constrained_expm_repeated_schedule"
      } else {
        "continuous_no_passaging_extended_o2"
      },
      initial_ploidy = f7ft_initial_ploidy(),
      day_values = f7e_day_values(smoke), o2_values = f7e_o2_values(smoke),
      p_misseg = p_values,
      optimizer_endpoint_weight = weights,
      n_lineage_schedule = if (run_paths$mode == "passage") 6L else 0L,
      mean_ploidy = values
    )
    f7ft_atomic_save_rds(
      object,
      file.path(
        run_paths$run_root,
        paste0("extended_time_panel_", tolower(letter), ".rds")
      ), compress = "gzip"
    )
    panel_objects[[letter]] <- object
  }
  validation <- do.call(rbind, lapply(panel_objects, function(object) {
    day0 <- object$mean_ploidy[, 1L, , , drop = FALSE]
    expected <- array(
      object$initial_ploidy, dim = dim(day0), dimnames = dimnames(day0)
    )
    data.frame(
      panel_letter = object$panel_letter, mode = run_paths$mode,
      n_initial_ploidy = length(object$initial_ploidy),
      n_day = length(object$day_values), n_o2 = length(object$o2_values),
      n_o2_0_to_0p5 = sum(
        object$o2_values >= 0 & object$o2_values <= 0.5 + 1e-12
      ),
      n_p_misseg = length(object$p_misseg),
      optimizer_weight = paste(object$optimizer_endpoint_weight, collapse = ","),
      minimum_mean_ploidy = min(object$mean_ploidy),
      maximum_mean_ploidy = max(object$mean_ploidy),
      maximum_day0_abs_error = max(abs(day0 - expected)),
      all_finite = all(is.finite(object$mean_ploidy)),
      passed = all(is.finite(object$mean_ploidy)) &&
        min(object$mean_ploidy) >= 1 - 1e-8 &&
        max(object$mean_ploidy) <= 7 + 1e-8 &&
        max(abs(day0 - expected)) <= 1e-10 &&
        length(object$initial_ploidy) == 5L &&
        length(object$day_values) == length(f7e_day_values(smoke)) &&
        isTRUE(all.equal(object$o2_values, f7e_o2_values(smoke))) &&
        (isTRUE(smoke) || all(object$optimizer_endpoint_weight == 50L)),
      stringsAsFactors = FALSE
    )
  }))
  f7ft_atomic_write_tsv(
    validation, file.path(run_paths$run_root, "extended_panel_validation.tsv")
  )
  if (nrow(validation) != 2L || !all(validation$passed)) {
    stop("Extended panel validation failed.")
  }
  list(panels = panel_objects, task_qc = task_qc, validation = validation)
}

f7e_compare_existing <- function(paths, run_paths, panels) {
  if (run_paths$mode == "passage") {
    source_paths <- f7p_paths(paths, run_id = NULL, create = FALSE)
    source_reader <- function(letter) f7p_read_panel_object(source_paths, letter)
  } else {
    source_paths <- f7ft_paths(paths, run_id = NULL, create = FALSE)
    source_reader <- function(letter) {
      path <- file.path(
        source_paths$run_root,
        paste0("finite_time_panel_", tolower(letter), ".rds")
      )
      f7r_require_files(
        path, paste0("Figure 7", letter, " continuous finite-time data")
      )
      object <- readRDS(path)
      if (!identical(object$profile, f7ft_profile()) ||
          !identical(object$panel_letter, letter)) {
        stop("Unexpected continuous finite-time panel object: ", path)
      }
      object
    }
  }
  rows <- lapply(names(panels), function(letter) {
    extended <- panels[[letter]]
    source <- source_reader(letter)
    day_index <- match(source$day_values, extended$day_values)
    o2_index <- match(
      round(source$o2_values, 10), round(extended$o2_values, 10)
    )
    keep <- which(!is.na(o2_index))
    if (anyNA(day_index) || !length(keep)) {
      stop("No common existing grid for extended comparison panel ", letter)
    }
    source_values <- source$mean_ploidy[, , keep, , drop = FALSE]
    extended_values <- extended$mean_ploidy[
      , day_index, o2_index[keep], , drop = FALSE
    ]
    delta <- extended_values - source_values
    data.frame(
      panel_letter = letter, mode = run_paths$mode,
      n_compared = length(delta), n_common_o2 = length(keep),
      comparison_day_min = min(source$day_values),
      comparison_day_max = max(source$day_values),
      maximum_absolute_delta = max(abs(delta)),
      mean_absolute_delta = mean(abs(delta)),
      passed = max(abs(delta)) <= 1e-10,
      stringsAsFactors = FALSE
    )
  })
  comparison <- do.call(rbind, rows)
  path <- f7ft_atomic_write_tsv(
    comparison,
    file.path(run_paths$run_root, "extended_vs_existing_validation.tsv")
  )
  if (!all(comparison$passed)) {
    stop("Extended common-grid results do not reproduce existing panels.")
  }
  path
}

f7e_fingerprint <- function(
    paths, base_run_paths, run_paths, schedule_bundle = NULL, smoke = FALSE
) {
  inputs <- c(
    file.path(base_run_paths$run_root, "q10_unique_endpoint_manifest.tsv"),
    file.path(base_run_paths$run_root, "q10_optimizer_seed_manifest.tsv"),
    f7p_external_semantics_files(paths)
  )
  if (!is.null(schedule_bundle)) {
    inputs <- c(
      inputs, schedule_bundle$source_schedule_path,
      schedule_bundle$schedule_path, schedule_bundle$fit_result_path
    )
  }
  paste(
    f7e_profile(run_paths$mode),
    paste0("fixed_model=", f7r_model_source_fingerprint(paths)),
    paste0("inputs=", paste(unname(tools::md5sum(inputs)), collapse = ":")),
    paste0("p=", paste(f7ft_p_values(), collapse = ",")),
    paste0("o2=", paste(f7e_o2_values(smoke), collapse = ",")),
    paste0("day=", paste(range(f7e_day_values(smoke)), collapse = ":")),
    sep = "|"
  )
}

f7e_publish_current <- function(paths, run_paths, fingerprint) {
  pointer <- data.frame(
    run_id = run_paths$run_id,
    relative_run_path = file.path(
      basename(run_paths$base), run_paths$run_id
    ),
    profile = f7e_profile(run_paths$mode), fingerprint = fingerprint,
    published_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  f7ft_atomic_write_tsv(pointer, run_paths$current)
}

f7e_data <- function(
    workspace_root = f7r_find_workspace_root(), mode = f7e_modes(),
    n_core = 1L, run_id = format(Sys.time(), "%Y%m%d_%H%M%S"),
    smoke = FALSE, publish_current = !isTRUE(smoke)
) {
  mode <- match.arg(mode, f7e_modes())
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  f7r_require_packages(c("Matrix", "Rcpp", "data.table"))
  paths <- f7r_paths(workspace_root)
  base_run_paths <- f7ft_paths(paths, run_id = NULL, create = FALSE)
  run_paths <- f7e_paths(paths, mode, run_id = run_id, create = TRUE)
  if (length(list.files(run_paths$cache, all.files = TRUE, no.. = TRUE))) {
    stop("Fresh extended run cache is not empty: ", run_paths$cache)
  }
  f7r_load_response_engine(paths)
  objective_bundle <- f7r_objective_selection(paths)
  endpoints <- f7r_read_tsv(file.path(
    base_run_paths$run_root, "q10_unique_endpoint_manifest.tsv"
  ))
  schedule_bundle <- NULL
  if (mode == "passage") {
    recorded <- f7p_extract_schedule(base_run_paths, paths)
    schedule_bundle <- f7e_expand_schedule(
      recorded, max_day = max(f7e_day_values(smoke))
    )
    schedule_bundle$source_schedule_path <- f7ft_atomic_write_tsv(
      schedule_bundle$source_schedule,
      file.path(run_paths$run_root, "recorded_passage_schedule.tsv")
    )
    schedule_bundle$schedule_path <- f7ft_atomic_write_tsv(
      schedule_bundle$schedule,
      file.path(run_paths$run_root, "repeated_passage_schedule.tsv")
    )
    f7ft_atomic_write_tsv(
      f7e_schedule_coverage(schedule_bundle),
      file.path(run_paths$run_root, "repeated_passage_time_coverage.tsv")
    )
  }
  fingerprint <- f7e_fingerprint(
    paths, base_run_paths, run_paths, schedule_bundle, smoke
  )
  contexts <- lapply(
    unique(endpoints$pair_id), f7r_pair_model_context,
    selected = objective_bundle$selected, paths = paths
  )
  names(contexts) <- unique(endpoints$pair_id)
  tasks <- f7e_task_manifest(
    endpoints, run_paths, smoke = smoke, chunk_size = 2L
  )
  task_list <- split(tasks, seq_len(nrow(tasks)))
  compute_one <- function(task) tryCatch(
    f7e_compute_task(
      task, endpoints, objective_bundle, contexts, paths, run_paths,
      schedule_bundle, fingerprint, smoke
    ),
    error = function(e) structure(
      list(task_id = task$task_id[[1L]], message = conditionMessage(e)),
      class = "f7e_error"
    )
  )
  message(
    "Figure 7 ", mode, " extended q10: ", nrow(tasks),
    " tasks, workers=", min(as.integer(n_core), nrow(tasks)),
    ", run_id=", run_paths$run_id
  )
  results <- f7ft_parallel_lapply(task_list, compute_one, n_core = n_core)
  failed <- vapply(results, inherits, logical(1L), "f7e_error")
  if (any(failed)) stop(
    "Figure 7 extended task failures: ",
    paste(vapply(results[failed], `[[`, character(1L), "message"), collapse = "; ")
  )
  aggregation <- f7e_aggregate(tasks, run_paths, fingerprint, smoke)
  comparison <- NULL
  if (!isTRUE(smoke)) {
    comparison <- f7e_compare_existing(paths, run_paths, aggregation$panels)
  }
  provenance <- data.frame(
    key = c(
      "profile", "run_id", "mode", "fingerprint", "workspace_root",
      "model_code_root", "joint_result_root", "base_continuous_run_id",
      "n_core", "smoke", "day_grid", "oxygen_grid", "aggregation",
      "passage_repeat_policy"
    ),
    value = c(
      f7e_profile(mode), run_paths$run_id, mode, fingerprint,
      normalizePath(paths$root, mustWork = TRUE),
      normalizePath(paths$oxygen_code, mustWork = TRUE),
      f7p_joint_result_root(), base_run_paths$run_id,
      as.character(n_core), as.character(smoke),
      paste(range(f7e_day_values(smoke)), collapse = ":"),
      paste0(
        paste(range(f7e_o2_values(smoke)), collapse = ":"),
        "; n=", length(f7e_o2_values(smoke))
      ),
      "optimizer multiplicity; equal lineage within cohort; equal cohort",
      if (mode == "passage") schedule_bundle$repeat_policy else "not applicable"
    ), stringsAsFactors = FALSE
  )
  f7ft_atomic_write_tsv(
    provenance, file.path(run_paths$run_root, "extended_run_provenance.tsv")
  )
  if (isTRUE(publish_current)) {
    f7e_publish_current(paths, run_paths, fingerprint)
  }
  invisible(list(
    paths = run_paths, tasks = tasks, aggregation = aggregation,
    comparison = comparison, schedule = schedule_bundle,
    fingerprint = fingerprint
  ))
}
