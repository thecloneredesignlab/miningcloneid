# Independent R propagation and the actual external fitting selector are the
# oracle. Production BLAS propagation must agree in state, day selection, and
# daily mean, including experiment-time boundaries (not just final means).
f7g_reference_segments <- function(step, initial, ploidy, n_day, seed, target,
                                  duration, selector) {
  state <- sweep(initial, 2L, colSums(initial), "/") * seed
  output <- matrix(NA_real_, ncol(initial), n_day + 1L)
  output[, 1L] <- colSums(state * ploidy) / colSums(state)
  picked <- matrix(NA_integer_, ncol(initial), n_day %/% duration)
  active <- rep(TRUE, ncol(initial))
  failure_day <- rep(NA_integer_, ncol(initial))
  offset <- 0L
  while (offset < n_day) {
    length <- min(duration, n_day - offset)
    candidates <- array(0, c(nrow(state), ncol(state), length + 1L))
    candidates[, , 1L] <- state
    for (day in seq_len(length)) {
      state <- step %*% state
      candidates[, , day + 1L] <- state
      output[, offset + day + 1L] <- colSums(state * ploidy) / colSums(state)
    }
    if (length == duration) for (col in which(active)) {
      trajectory <- t(candidates[, col, ])
      selection <- selector(
        sim = list(Ntot_live_obs = rowSums(trajectory), live_state_obs = trajectory),
        reseed_live_cells = seed, grid_pre = ploidy,
        target_live_cells = target, obs_days_local = 0:length
      )
      if (!isTRUE(selection$passage_executed)) {
        active[col] <- FALSE
        failure_day[col] <- offset + length
        state[, col] <- 0
        output[col, offset + length + 1L] <- NA_real_
        next
      }
      state[, col] <- selection$selected_frac * seed
      picked[col, offset %/% duration + 1L] <- as.integer(selection$selected_day)
      output[col, offset + length + 1L] <- sum(selection$selected_frac * ploidy)
    }
    offset <- offset + length
  }
  list(mean_ploidy = output, final_state = sweep(state, 2L, colSums(state), "/"),
       selected_day_trace = picked, protocol_failure_day = failure_day)
}

f7g_validate_segments <- function(paths, run_paths, model_cases = list()) {
  helper <- Sys.getenv("FIGURE7_SELECTOR_REFERENCE", file.path(paths$oxygen_code,
    "util", "o2_supply_demand_map_invitro_lineage_simulation_utils.R"))
  oracle <- new.env(parent = globalenv())
  sys.source(helper, envir = oracle)
  if (!exists("ivt_select_segment_observation", envir = oracle, inherits = FALSE)) {
    stop("The selected reference lacks the current HPC inoculum-supply constraint.")
  }
  # This scalar return is not used in the comparison; selection and reseeding
  # above are executed by the unmodified fitting function.
  oracle$ivt_weighted_mean_kary_N <- function(frac, grid_pre) sum(frac * grid_pre) / sum(frac)
  initial <- diag(2)
  coupled <- as.matrix(Matrix::expm(matrix(c(0.45, 0.12, 0.03, 0.1), 2L)))
  cases <- c(list(
    list(name = "unreachable_target", step = coupled, initial = initial,
         ploidy = c(2, 4), seed = 100, target = 1e9, days = 25L),
    list(name = "early_target_selection", step = coupled, initial = initial,
         ploidy = c(2, 4), seed = 100, target = 230, days = 25L),
    list(name = "declining_population", step = diag(c(0.7, 0.8)), initial = initial,
         ploidy = c(2, 4), seed = 100, target = 500, days = 25L),
    list(name = "fallback_to_smallest_eligible_population",
         step = matrix(c(0, .5, 100, 0), 2), initial = initial,
         ploidy = c(2, 4), seed = 100, target = 150, days = 25L),
    list(name = "exact_tie_earliest_positive_day", step = diag(2), initial = initial,
         ploidy = c(2, 4), seed = 100, target = 200, days = 25L),
    list(name = "unfinished_last_segment", step = coupled, initial = initial,
         ploidy = c(2, 4), seed = 100, target = 230, days = 13L)
  ), model_cases)
  rows <- lapply(cases, function(z) {
    duration <- z$duration %||% 5L
    actual <- f7g_propagate_segments_cpp(z$step, z$initial, z$ploidy,
      z$days, log(z$seed), log(z$target), duration, TRUE)
    expected <- f7g_reference_segments(z$step, z$initial, z$ploidy,
      z$days, z$seed, z$target, duration, oracle$ivt_extract_passage_end_state)
    finite_error <- function(a, b) {
      d <- abs(a - b)
      if (any(is.finite(d))) max(d[is.finite(d)]) else 0
    }
    mean_error <- finite_error(actual$mean_ploidy, expected$mean_ploidy)
    state_error <- finite_error(actual$final_state, expected$final_state)
    missing_match <- identical(is.na(actual$mean_ploidy), is.na(expected$mean_ploidy)) &&
      identical(is.na(actual$final_state), is.na(expected$final_state))
    failure_match <- identical(actual$protocol_failure_day, expected$protocol_failure_day)
    selection_match <- identical(actual$selected_day_trace, expected$selected_day_trace)
    data.frame(case = z$name, mean_error = mean_error, state_error = state_error,
      selection_match = selection_match,
      infeasible_mask_match = missing_match, failure_day_match = failure_match,
      n_failed_initial_conditions = sum(!is.na(actual$protocol_failure_day)),
      minimum_selected_day = if (all(is.na(actual$selected_day_trace))) NA_integer_ else min(actual$selected_day_trace, na.rm = TRUE),
      maximum_selected_day = if (all(is.na(actual$selected_day_trace))) NA_integer_ else max(actual$selected_day_trace, na.rm = TRUE),
      external_selector_md5 = unname(tools::md5sum(helper)),
      passed = mean_error < 1e-10 && state_error < 1e-10 && selection_match && missing_match && failure_match)
  })
  result <- do.call(rbind, rows)
  f7ft_atomic_write_tsv(result, file.path(run_paths$run_root, "segment_selector_validation.tsv"))
  if (!all(result$passed)) stop("Fitting-selector comparison failed: ",
                              paste(result$case[!result$passed], collapse = ", "))
  result
}
