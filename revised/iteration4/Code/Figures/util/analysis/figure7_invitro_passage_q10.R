#!/usr/bin/env Rscript

# Passage-constrained finite-time q10 response surfaces for Figure 7E-F.
#
# The fixed-O2 transition operator is evaluated from the required external
# O2_supply_demand_MAP checkout.  Passage durations, target-cell state
# selection, and reseeding schedules are recovered from the latest joint-fit
# RDS objects.  Experimental time advances by the recorded passage duration;
# the state transferred at a boundary is the daily state closest to the
# recorded target live-cell count, exactly as in the fitting lineage code.

options(stringsAsFactors = FALSE, warn = 1)

f7p_profile <- function() "figure7_fixed_pmisseg_invitro_passage_q10_v1"

f7p_joint_result_root <- function() {
  root <- trimws(Sys.getenv(
    "FIGURE_JOINT_RESULT_ROOT",
    unset = paste0(
      "/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/",
      "fit_joint_unified_global_invitro_500seed_all_xxlarge_r442_exact_20260828_145253"
    )
  ))
  normalizePath(root, mustWork = TRUE)
}

f7p_paths <- function(paths, run_id = NULL, create = FALSE) {
  base <- file.path(paths$figure7, "invitro_passage_q10_runs")
  current <- file.path(paths$figure7, "invitro_passage_q10_current.tsv")
  if (is.null(run_id)) {
    f7r_require_files(current, "current Figure 7 in-vitro passage run pointer")
    pointer <- f7r_read_tsv(current)
    if (nrow(pointer) != 1L || !all(c("run_id", "relative_run_path") %in% names(pointer))) {
      stop("Malformed Figure 7 in-vitro passage pointer: ", current)
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
    run_id = run_id,
    base = base,
    run_root = run_root,
    cache = file.path(run_root, "task_cache"),
    panels = file.path(run_root, "panels"),
    rendered = file.path(run_root, "rendered"),
    current = current,
    supp7_8 = file.path(paths$root, "data", "Figures", "Supp_Figure7_8")
  )
  if (isTRUE(create)) {
    invisible(lapply(
      unname(unlist(out[c("base", "run_root", "cache", "panels", "rendered", "supp7_8")])),
      dir.create, recursive = TRUE, showWarnings = FALSE
    ))
  }
  out
}

f7p_read_panel_object <- function(run_paths, panel_letter) {
  path <- file.path(
    run_paths$run_root,
    paste0("passage_time_panel_", tolower(panel_letter), ".rds")
  )
  f7r_require_files(path, paste0("Figure 7", panel_letter, " passage data"))
  object <- readRDS(path)
  if (!identical(object$profile, f7p_profile()) ||
      !identical(object$panel_letter, panel_letter)) {
    stop("Unexpected passage panel object: ", path)
  }
  object
}

f7p_external_semantics_files <- function(paths) {
  files <- file.path(
    paths$oxygen_code,
    c(
      "util/o2_supply_demand_map_invitro_lineage_utils.R",
      "util/o2_supply_demand_map_invitro_lineage_simulation_utils.R",
      "util/o2_supply_demand_map_invitro_objective_utils.R"
    )
  )
  f7r_require_files(files, "external in-vitro passage implementation")
  normalizePath(files, mustWork = TRUE)
}

f7p_trace_path <- function(ids, parents, terminal) {
  path <- character()
  current <- terminal
  repeat {
    path <- c(current, path)
    index <- match(current, ids)
    if (is.na(index)) stop("Missing lineage segment: ", current)
    parent <- parents[[index]]
    if (is.na(parent) || !nzchar(parent)) break
    current <- parent
  }
  path
}

f7p_terminal_group <- function(x) {
  if (!grepl("::", x, fixed = TRUE)) return("control")
  if (grepl("::O1$", x)) return("O1")
  if (grepl("::O2$", x)) return("O2")
  NA_character_
}

f7p_extract_schedule <- function(base_run_paths, paths) {
  manifest <- f7r_read_tsv(file.path(
    base_run_paths$run_root, "q10_unique_endpoint_manifest.tsv"
  ))
  source_endpoint <- manifest[
    manifest$model_context == "in vitro" & manifest$pair_label == "C01",
    , drop = FALSE
  ][1L, , drop = FALSE]
  fit_result_path <- file.path(
    f7p_joint_result_root(), source_endpoint$pair_id[[1L]],
    source_endpoint$representative_seed[[1L]], "fit_result.rds"
  )
  f7r_require_files(fit_result_path, "joint fit result carrying passage adapters")
  fit <- readRDS(fit_result_path)
  invitro <- fit$best_components$invitro
  if (is.null(invitro$run_2N) || is.null(invitro$run_4N)) {
    stop("Joint fit result does not contain fitted in-vitro lineage runs.")
  }

  output <- list()
  output_index <- 0L
  for (cohort in c("2N", "4N")) {
    run <- invitro[[paste0("run_", cohort)]]
    segments <- run$adapter$segments
    ids <- vapply(segments, `[[`, character(1L), "segment_id")
    parents <- vapply(segments, function(segment) {
      parent <- segment$parent_segment_id
      if (is.null(parent) || !length(parent)) NA_character_ else as.character(parent[[1L]])
    }, character(1L))
    terminals <- ids[!ids %in% stats::na.omit(parents)]
    terminal_table <- data.frame(
      terminal = terminals,
      lineage_group = vapply(terminals, f7p_terminal_group, character(1L)),
      stringsAsFactors = FALSE
    )
    terminal_table$path <- lapply(
      terminal_table$terminal,
      f7p_trace_path, ids = ids, parents = parents
    )
    terminal_table$total_duration <- vapply(terminal_table$path, function(path) {
      sum(vapply(segments[match(path, ids)], `[[`, numeric(1L), "duration_days"))
    }, numeric(1L))
    selected_terminals <- do.call(rbind, lapply(
      c("control", "O1", "O2"),
      function(group) {
        candidates <- terminal_table[terminal_table$lineage_group == group, , drop = FALSE]
        if (!nrow(candidates)) stop("Missing ", cohort, " ", group, " terminal lineage.")
        candidates[order(-candidates$total_duration, candidates$terminal), , drop = FALSE][1L, ]
      }
    ))
    for (terminal_index in seq_len(nrow(selected_terminals))) {
      terminal <- selected_terminals$terminal[[terminal_index]]
      group <- selected_terminals$lineage_group[[terminal_index]]
      path <- selected_terminals$path[[terminal_index]]
      cumulative <- 0
      for (passage_index in seq_along(path)) {
        segment_index <- match(path[[passage_index]], ids)
        segment <- segments[[segment_index]]
        selection <- run$segment_results[[segment_index]]$selection
        duration <- as.numeric(segment$duration_days)
        output_index <- output_index + 1L
        output[[output_index]] <- data.frame(
          cohort = cohort,
          lineage = paste(cohort, group, sep = "_"),
          lineage_group = group,
          terminal_segment_id = terminal,
          passage_index = passage_index,
          segment_id = segment$segment_id,
          parent_segment_id = if (is.null(segment$parent_segment_id)) NA_character_ else as.character(segment$parent_segment_id),
          source_oxygen_pct = as.numeric(segment$oxygen_pct),
          passage_duration = duration,
          experimental_day_start = cumulative,
          experimental_day_end = cumulative + duration,
          initial_cells = as.numeric(segment$initial_cells),
          final_cells = as.numeric(segment$final_cells),
          target_live_cells = as.numeric(selection$target_live_cells),
          fitted_selected_day = as.numeric(selection$selected_day),
          stringsAsFactors = FALSE
        )
        cumulative <- cumulative + duration
      }
    }
  }
  schedule <- do.call(rbind, output)
  rownames(schedule) <- NULL
  if (length(unique(schedule$lineage)) != 6L ||
      !identical(sort(unique(schedule$cohort)), c("2N", "4N")) ||
      !all(schedule$passage_duration > 0) ||
      any(abs(schedule$passage_duration - round(schedule$passage_duration)) > 1e-9)) {
    stop("Recovered passage schedule failed the six-lineage duration contract.")
  }
  list(
    schedule = schedule,
    fit_result_path = normalizePath(fit_result_path, mustWork = TRUE),
    schedule_path = NULL,
    max_experimental_day = as.integer(max(schedule$experimental_day_end))
  )
}

f7p_schedule_coverage <- function(schedule_bundle) {
  schedule <- schedule_bundle$schedule
  days <- 0:schedule_bundle$max_experimental_day
  rows <- lapply(days, function(day) {
    terminal <- stats::aggregate(
      experimental_day_end ~ cohort + lineage, schedule, max
    )
    active <- terminal$experimental_day_end >= day
    data.frame(
      day = day,
      active_lineages_2N = sum(active & terminal$cohort == "2N"),
      active_lineages_4N = sum(active & terminal$cohort == "4N"),
      active_lineages_total = sum(active),
      active_cohorts = sum(vapply(c("2N", "4N"), function(cohort) {
        any(active & terminal$cohort == cohort)
      }, logical(1L))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

f7p_fingerprint <- function(paths, base_run_paths, schedule_bundle) {
  semantics <- f7p_external_semantics_files(paths)
  inputs <- c(
    file.path(base_run_paths$run_root, "q10_unique_endpoint_manifest.tsv"),
    file.path(base_run_paths$run_root, "q10_optimizer_seed_manifest.tsv"),
    schedule_bundle$schedule_path,
    schedule_bundle$fit_result_path,
    semantics
  )
  paste(
    f7p_profile(),
    paste0("fixed_model=", f7r_model_source_fingerprint(paths)),
    paste0("inputs=", paste(unname(tools::md5sum(inputs)), collapse = ":")),
    paste0("p=", paste(f7ft_p_values(), collapse = ",")),
    paste0("o2=", paste(range(f7ft_o2_values()), collapse = ":"), ":", length(f7ft_o2_values())),
    paste0("experimental_day=0:", schedule_bundle$max_experimental_day),
    sep = "|"
  )
}

f7p_task_manifest <- function(endpoints, run_paths, smoke = FALSE, chunk_size = 2L) {
  endpoints <- endpoints[endpoints$model_context == "in vitro", , drop = FALSE]
  if (isTRUE(smoke)) {
    endpoints <- do.call(rbind, lapply(split(endpoints, endpoints$pair_label), function(x) x[1L, , drop = FALSE]))
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
    chunks <- split(seq_len(nrow(family_endpoints)), ceiling(seq_len(nrow(family_endpoints)) / chunk_size))
    for (p_value in p_values) for (chunk_index in seq_along(chunks)) {
      selected <- family_endpoints[chunks[[chunk_index]], , drop = FALSE]
      index <- index + 1L
      rows[[index]] <- data.frame(
        task_id = sprintf("P%04d", index),
        pair_id = family_endpoints$pair_id[[1L]],
        pair_label = family,
        panel_letter = if (family == "C01") "E" else "F",
        p_misseg = p_value,
        chunk_index = chunk_index,
        endpoint_indices = paste(selected$endpoint_index, collapse = ","),
        n_unique_endpoint = nrow(selected),
        represented_optimizer_endpoint = sum(selected$endpoint_multiplicity_q10),
        cache_path = file.path(run_paths$cache, sprintf("passage_task_%04d.rds", index)),
        stringsAsFactors = FALSE
      )
    }
  }
  tasks <- do.call(rbind, rows)
  f7ft_atomic_write_tsv(tasks, file.path(run_paths$run_root, "passage_task_manifest.tsv"))
  tasks
}

f7p_clean_state <- function(x) {
  x <- Re(as.numeric(x))
  x[x < 0 & x > -1e-8] <- 0
  x[!is.finite(x) | x < 0] <- NA_real_
  x
}

f7p_state_mean <- function(state, ploidy_grid) {
  total <- sum(state)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  sum(state * ploidy_grid) / total
}

f7p_lineage_trajectory <- function(
    step, ploidy_grid, ngrid, n_unit, initial_ploidy, lineage_schedule,
    default_initial_cells, max_day = NULL
) {
  schedule_end <- as.integer(max(lineage_schedule$experimental_day_end))
  if (is.null(max_day)) max_day <- schedule_end
  max_day <- as.integer(max_day)
  if (!is.finite(max_day) || max_day < 0L || max_day > schedule_end) {
    stop("max_day must be between zero and the lineage schedule end.")
  }
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
    recorded_duration <- as.integer(round(segment$passage_duration[[1L]]))
    duration <- min(recorded_duration, max_day - current_day)
    initial_cells <- as.numeric(segment$initial_cells[[1L]])
    if (!is.finite(initial_cells) || initial_cells <= 0) initial_cells <- default_initial_cells
    state <- state_fraction * initial_cells
    states <- matrix(NA_real_, nrow = length(ngrid), ncol = duration + 1L)
    states[, 1L] <- state
    if (duration > 0L) for (day_index in seq_len(duration)) {
      state <- f7p_clean_state(step %*% state)
      if (anyNA(state)) stop("Non-finite passage state.")
      states[, day_index + 1L] <- state
    }
    live_cells <- colSums(states)
    if (duration < recorded_duration) {
      current_day <- current_day + duration
      break
    }
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
    if (!is.finite(selected_total) || selected_total <= 0) stop("Passage selection returned an empty state.")
    state_fraction <- selected_state / selected_total

    means <- apply(states, 2L, f7p_state_mean, ploidy_grid = ploidy_grid)
    indices <- current_day + 0:duration + 1L
    values[indices] <- means
    values[[current_day + duration + 1L]] <- f7p_state_mean(state_fraction, ploidy_grid)
    current_day <- current_day + duration
  }
  list(values = values, selected_days = selected_days)
}

f7p_operator_response <- function(
    M, ngrid, n_unit, schedule_bundle, default_initial_cells,
    day_values = NULL
) {
  step <- as.matrix(Matrix::expm(M))
  ploidy_grid <- as.numeric(ngrid) / as.numeric(n_unit)
  days <- if (is.null(day_values)) {
    0:schedule_bundle$max_experimental_day
  } else {
    as.integer(day_values)
  }
  if (!identical(days, 0:max(days))) {
    stop("Passage-aware propagation requires a complete integer day grid from zero.")
  }
  response <- matrix(
    NA_real_, nrow = length(f7ft_initial_ploidy()), ncol = length(days),
    dimnames = list(paste0(f7ft_initial_ploidy(), "N"), as.character(days))
  )
  n_selection <- 0L
  n_selected_before_duration <- 0L
  minimum_selected_day <- Inf
  maximum_selected_day <- -Inf

  for (initial_index in seq_along(f7ft_initial_ploidy())) {
    cohort_trajectories <- list()
    for (cohort in c("2N", "4N")) {
      lineages <- unique(schedule_bundle$schedule$lineage[
        schedule_bundle$schedule$cohort == cohort
      ])
      lineage_values <- matrix(NA_real_, nrow = length(lineages), ncol = length(days))
      for (lineage_index in seq_along(lineages)) {
        lineage_schedule <- schedule_bundle$schedule[
          schedule_bundle$schedule$lineage == lineages[[lineage_index]],
          , drop = FALSE
        ]
        lineage_schedule <- lineage_schedule[order(lineage_schedule$passage_index), , drop = FALSE]
        trajectory <- f7p_lineage_trajectory(
          step = step, ploidy_grid = ploidy_grid, ngrid = ngrid,
          n_unit = n_unit,
          initial_ploidy = f7ft_initial_ploidy()[[initial_index]],
          lineage_schedule = lineage_schedule,
          default_initial_cells = default_initial_cells,
          max_day = min(
            max(days), as.integer(max(lineage_schedule$experimental_day_end))
          )
        )
        lineage_values[lineage_index, seq_along(trajectory$values)] <- trajectory$values
        duration <- as.integer(round(lineage_schedule$passage_duration))
        n_selection <- n_selection + length(trajectory$selected_days)
        n_selected_before_duration <- n_selected_before_duration +
          sum(trajectory$selected_days < duration)
        minimum_selected_day <- min(minimum_selected_day, trajectory$selected_days)
        maximum_selected_day <- max(maximum_selected_day, trajectory$selected_days)
      }
      cohort_trajectories[[cohort]] <- apply(lineage_values, 2L, function(x) {
        if (all(!is.finite(x))) NA_real_ else mean(x[is.finite(x)])
      })
    }
    cohort_matrix <- do.call(rbind, cohort_trajectories)
    response[initial_index, ] <- apply(cohort_matrix, 2L, function(x) {
      if (all(!is.finite(x))) NA_real_ else mean(x[is.finite(x)])
    })
  }
  list(
    mean_ploidy = response,
    selection = c(
      n_selection = n_selection,
      n_selected_before_duration = n_selected_before_duration,
      minimum_selected_day = minimum_selected_day,
      maximum_selected_day = maximum_selected_day
    )
  )
}

f7p_compute_task <- function(
    task, endpoints, objective_bundle, contexts, paths, run_paths,
    schedule_bundle, fingerprint, smoke = FALSE
) {
  f7r_load_response_engine(paths)
  indices <- as.integer(strsplit(task$endpoint_indices[[1L]], ",", fixed = TRUE)[[1L]])
  selected_endpoints <- endpoints[match(indices, endpoints$endpoint_index), , drop = FALSE]
  if (anyNA(selected_endpoints$endpoint_index)) stop("Passage task endpoint lookup failed.")
  o2_values <- f7ft_o2_values(smoke)
  days <- 0:schedule_bundle$max_experimental_day
  weighted_sum <- array(
    0,
    dim = c(length(f7ft_initial_ploidy()), length(days), length(o2_values)),
    dimnames = list(
      initial_ploidy = paste0(f7ft_initial_ploidy(), "N"),
      day = as.character(days), O2_pct = as.character(o2_values)
    )
  )
  p_value <- as.numeric(task$p_misseg[[1L]])
  selection_total <- c(
    n_selection = 0, n_selected_before_duration = 0,
    minimum_selected_day = Inf, maximum_selected_day = -Inf
  )
  operator_count <- 0L
  maximum_override_error <- 0
  maximum_formula_error <- 0
  for (endpoint_index in seq_len(nrow(selected_endpoints))) {
    endpoint <- selected_endpoints[endpoint_index, , drop = FALSE]
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
      response <- f7p_operator_response(
        fixed$M, fixed$ngrid, prepared$config$N_UNIT %||% 22L,
        schedule_bundle, prepared$config$init_total_size %||% 1e6
      )
      weighted_sum[, , o2_index] <- weighted_sum[, , o2_index] +
        weight * response$mean_ploidy
      selection_total[c("n_selection", "n_selected_before_duration")] <-
        selection_total[c("n_selection", "n_selected_before_duration")] +
        response$selection[c("n_selection", "n_selected_before_duration")]
      selection_total[["minimum_selected_day"]] <- min(
        selection_total[["minimum_selected_day"]],
        response$selection[["minimum_selected_day"]]
      )
      selection_total[["maximum_selected_day"]] <- max(
        selection_total[["maximum_selected_day"]],
        response$selection[["maximum_selected_day"]]
      )
      operator_count <- operator_count + 1L
    }
  }
  day0 <- weighted_sum[, 1L, , drop = FALSE] /
    sum(selected_endpoints$endpoint_multiplicity_q10)
  expected <- array(
    f7ft_initial_ploidy(), dim = dim(day0), dimnames = dimnames(day0)
  )
  qc <- data.frame(
    task_id = task$task_id[[1L]], pair_label = task$pair_label[[1L]],
    p_misseg = p_value,
    n_unique_endpoint = nrow(selected_endpoints),
    represented_optimizer_endpoint = sum(selected_endpoints$endpoint_multiplicity_q10),
    n_operator = operator_count,
    n_o2 = length(o2_values), n_day = length(days),
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
    profile = f7p_profile(), fingerprint = fingerprint, task = task,
    endpoints = selected_endpoints, o2_values = o2_values,
    day_values = days, weighted_sum = weighted_sum, qc = qc
  ), task$cache_path[[1L]], compress = FALSE)
  qc
}

f7p_aggregate <- function(
    tasks, endpoints, run_paths, schedule_bundle, fingerprint, smoke = FALSE
) {
  objects <- lapply(tasks$cache_path, readRDS)
  if (any(!vapply(objects, function(x) identical(x$fingerprint, fingerprint), logical(1L)))) {
    stop("Passage task caches do not share the current fingerprint.")
  }
  task_qc <- do.call(rbind, lapply(objects, `[[`, "qc"))
  f7ft_atomic_write_tsv(task_qc, file.path(run_paths$run_root, "passage_task_validation.tsv"))
  if (!all(task_qc$passed)) stop("One or more passage tasks failed validation.")

  panel_objects <- list()
  for (family in f7ft_family_levels()) {
    letter <- if (family == "C01") "E" else "F"
    p_values <- sort(unique(tasks$p_misseg[tasks$pair_label == family]))
    arrays <- vector("list", length(p_values))
    weights <- integer(length(p_values))
    for (p_index in seq_along(p_values)) {
      keep <- which(tasks$pair_label == family &
        abs(tasks$p_misseg - p_values[[p_index]]) < 1e-12)
      selected <- objects[keep]
      arrays[[p_index]] <- Reduce(`+`, lapply(selected, `[[`, "weighted_sum"))
      weights[[p_index]] <- sum(vapply(
        selected,
        function(x) sum(x$endpoints$endpoint_multiplicity_q10), integer(1L)
      ))
      expected_weight <- if (isTRUE(smoke)) weights[[p_index]] else 50L
      if (weights[[p_index]] != expected_weight) {
        stop("Unexpected optimizer endpoint weight for passage panel ", letter)
      }
      arrays[[p_index]] <- arrays[[p_index]] / weights[[p_index]]
    }
    values <- simplify2array(arrays)
    dimnames(values)[[4L]] <- format(p_values, scientific = FALSE, trim = TRUE)
    object <- list(
      profile = f7p_profile(), fingerprint = fingerprint,
      panel_letter = letter, pair_label = family,
      model_context = "in vitro with passaging",
      propagation_mode = "passage_constrained_expm",
      initial_ploidy = f7ft_initial_ploidy(),
      day_values = 0:schedule_bundle$max_experimental_day,
      o2_values = f7ft_o2_values(smoke),
      p_misseg = p_values,
      optimizer_endpoint_weight = weights,
      n_lineage_schedule = 6L,
      mean_ploidy = values
    )
    f7ft_atomic_save_rds(
      object,
      file.path(run_paths$run_root, paste0("passage_time_panel_", tolower(letter), ".rds")),
      compress = "gzip"
    )
    panel_objects[[letter]] <- object
  }
  validation <- do.call(rbind, lapply(panel_objects, function(object) {
    day0 <- object$mean_ploidy[, 1L, , , drop = FALSE]
    expected <- array(object$initial_ploidy, dim = dim(day0), dimnames = dimnames(day0))
    data.frame(
      panel_letter = object$panel_letter,
      n_initial_ploidy = length(object$initial_ploidy),
      n_day = length(object$day_values), n_o2 = length(object$o2_values),
      n_p_misseg = length(object$p_misseg),
      n_lineage_schedule = object$n_lineage_schedule,
      minimum_mean_ploidy = min(object$mean_ploidy),
      maximum_mean_ploidy = max(object$mean_ploidy),
      maximum_day0_abs_error = max(abs(day0 - expected)),
      all_finite = all(is.finite(object$mean_ploidy)),
      passed = all(is.finite(object$mean_ploidy)) &&
        min(object$mean_ploidy) >= 1 - 1e-8 &&
        max(object$mean_ploidy) <= 7 + 1e-8 &&
        max(abs(day0 - expected)) <= 1e-10 &&
        length(object$initial_ploidy) == 5L &&
        object$n_lineage_schedule == 6L,
      stringsAsFactors = FALSE
    )
  }))
  f7ft_atomic_write_tsv(validation, file.path(run_paths$run_root, "passage_panel_validation.tsv"))
  if (nrow(validation) != 2L || !all(validation$passed)) {
    stop("Figure 7 passage panel validation failed.")
  }
  list(panels = panel_objects, task_qc = task_qc, validation = validation)
}

f7p_compare_continuous <- function(base_run_paths, panel_objects, run_paths) {
  rows <- lapply(names(panel_objects), function(letter) {
    passage <- panel_objects[[letter]]
    continuous_path <- file.path(
      base_run_paths$run_root,
      paste0("finite_time_panel_", tolower(letter), ".rds")
    )
    f7r_require_files(
      continuous_path,
      paste0("continuous Figure 7", letter, " comparison panel")
    )
    continuous <- readRDS(continuous_path)
    if (!identical(continuous$profile, f7ft_profile()) ||
        !identical(continuous$panel_letter, letter)) {
      stop("Unexpected continuous comparison panel: ", continuous_path)
    }
    day_index <- match(passage$day_values, continuous$day_values)
    if (anyNA(day_index) ||
        !identical(passage$initial_ploidy, continuous$initial_ploidy) ||
        !identical(passage$o2_values, continuous$o2_values) ||
        !isTRUE(all.equal(passage$p_misseg, continuous$p_misseg))) {
      stop("Continuous and passage-aware panel grids are not comparable for ", letter)
    }
    reference <- continuous$mean_ploidy[, day_index, , , drop = FALSE]
    delta <- passage$mean_ploidy - reference
    data.frame(
      panel_letter = letter,
      n_compared = length(delta),
      mean_delta = mean(delta),
      mean_absolute_delta = mean(abs(delta)),
      q95_absolute_delta = unname(stats::quantile(abs(delta), 0.95)),
      maximum_absolute_delta = max(abs(delta)),
      final_day = max(passage$day_values),
      final_day_mean_absolute_delta = mean(abs(delta[, dim(delta)[2L], , ])),
      stringsAsFactors = FALSE
    )
  })
  f7ft_atomic_write_tsv(
    do.call(rbind, rows),
    file.path(run_paths$run_root, "passage_vs_continuous_comparison.tsv")
  )
}

f7p_publish_current <- function(paths, run_paths, fingerprint) {
  pointer <- data.frame(
    run_id = run_paths$run_id,
    relative_run_path = file.path("invitro_passage_q10_runs", run_paths$run_id),
    profile = f7p_profile(), fingerprint = fingerprint,
    published_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  f7ft_atomic_write_tsv(pointer, run_paths$current)
}

f7p_data <- function(
    workspace_root = f7r_find_workspace_root(), n_core = 1L,
    run_id = f7ft_resolve_run_id(), smoke = FALSE,
    publish_current = !isTRUE(smoke)
) {
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  f7r_require_packages(c("Matrix", "Rcpp", "data.table"))
  paths <- f7r_paths(workspace_root)
  base_run_paths <- f7ft_paths(paths, run_id = NULL, create = FALSE)
  run_paths <- f7p_paths(paths, run_id = run_id, create = TRUE)
  if (length(list.files(run_paths$cache, all.files = TRUE, no.. = TRUE))) {
    stop("Fresh Figure 7 passage run cache is not empty: ", run_paths$cache)
  }
  f7r_load_response_engine(paths)
  objective_bundle <- f7r_objective_selection(paths)
  endpoints <- f7r_read_tsv(file.path(
    base_run_paths$run_root, "q10_unique_endpoint_manifest.tsv"
  ))
  schedule_bundle <- f7p_extract_schedule(base_run_paths, paths)
  schedule_bundle$schedule_path <- f7ft_atomic_write_tsv(
    schedule_bundle$schedule,
    file.path(run_paths$run_root, "passage_schedule.tsv")
  )
  coverage <- f7p_schedule_coverage(schedule_bundle)
  f7ft_atomic_write_tsv(coverage, file.path(run_paths$run_root, "passage_time_coverage.tsv"))
  fingerprint <- f7p_fingerprint(paths, base_run_paths, schedule_bundle)
  contexts <- lapply(
    unique(endpoints$pair_id), f7r_pair_model_context,
    selected = objective_bundle$selected, paths = paths
  )
  names(contexts) <- unique(endpoints$pair_id)
  tasks <- f7p_task_manifest(endpoints, run_paths, smoke = smoke, chunk_size = 2L)
  task_list <- split(tasks, seq_len(nrow(tasks)))
  compute_one <- function(task) tryCatch(
    f7p_compute_task(
      task, endpoints, objective_bundle, contexts, paths, run_paths,
      schedule_bundle, fingerprint, smoke = smoke
    ),
    error = function(e) structure(
      list(task_id = task$task_id[[1L]], message = conditionMessage(e)),
      class = "f7p_error"
    )
  )
  message(
    "Figure 7 passage q10: ", nrow(tasks), " tasks, workers=",
    min(as.integer(n_core), nrow(tasks)), ", run_id=", run_paths$run_id
  )
  results <- f7ft_parallel_lapply(task_list, compute_one, n_core = n_core)
  failed <- vapply(results, inherits, logical(1L), "f7p_error")
  if (any(failed)) stop(
    "Figure 7 passage task failures: ",
    paste(vapply(results[failed], `[[`, character(1L), "message"), collapse = "; ")
  )
  aggregation <- f7p_aggregate(
    tasks, endpoints, run_paths, schedule_bundle, fingerprint, smoke = smoke
  )
  comparison <- NULL
  if (!isTRUE(smoke)) {
    comparison <- f7p_compare_continuous(base_run_paths, aggregation$panels, run_paths)
  }
  provenance <- data.frame(
    key = c(
      "profile", "run_id", "fingerprint", "workspace_root",
      "model_code_root", "joint_result_root", "passage_source_fit_result",
      "base_continuous_run_id", "n_core", "smoke", "aggregation",
      "experimental_clock", "boundary_state", "oxygen_counterfactual"
    ),
    value = c(
      f7p_profile(), run_paths$run_id, fingerprint,
      normalizePath(paths$root, mustWork = TRUE),
      normalizePath(paths$oxygen_code, mustWork = TRUE),
      f7p_joint_result_root(), schedule_bundle$fit_result_path,
      base_run_paths$run_id, as.character(n_core), as.character(smoke),
      "optimizer multiplicity; equal lineage within cohort; equal covered cohort",
      "cumulative recorded passage_duration; no extrapolation past lineage end",
      "closest positive daily state to target_live_cells, then reseed by next passage initial_cells",
      "fixed counterfactual O2 value across every passage"
    ), stringsAsFactors = FALSE
  )
  f7ft_atomic_write_tsv(
    provenance, file.path(run_paths$run_root, "passage_run_provenance.tsv")
  )
  if (isTRUE(publish_current)) f7p_publish_current(paths, run_paths, fingerprint)
  invisible(list(
    paths = run_paths, base_paths = base_run_paths, schedule = schedule_bundle,
    tasks = tasks, aggregation = aggregation, comparison = comparison
  ))
}
