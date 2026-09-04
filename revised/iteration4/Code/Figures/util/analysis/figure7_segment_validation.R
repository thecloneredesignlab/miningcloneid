# Independent R propagation and the actual external fitting selector are the
# oracle. Production BLAS propagation must agree in state, day selection, and
# daily mean, including experiment-time boundaries (not just final means).
f7g_reference_segments <- function(step, initial, ploidy, n_day, seed, target,
                                  duration, selector) {
  state <- sweep(initial, 2L, colSums(initial), "/") * seed
  output <- matrix(NA_real_, ncol(initial), n_day + 1L)
  output[, 1L] <- colSums(state * ploidy) / colSums(state)
  picked <- matrix(NA_integer_, ncol(initial), n_day %/% duration)
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
    if (length == duration) for (col in seq_len(ncol(initial))) {
      trajectory <- t(candidates[, col, ])
      selection <- selector(
        sim = list(Ntot_live_obs = rowSums(trajectory), live_state_obs = trajectory),
        reseed_live_cells = seed, grid_pre = ploidy,
        target_live_cells = target, obs_days_local = 0:length
      )
      state[, col] <- selection$selected_frac * seed
      picked[col, offset %/% duration + 1L] <- as.integer(selection$selected_day)
      output[col, offset + length + 1L] <- sum(selection$selected_frac * ploidy)
    }
    offset <- offset + length
  }
  list(mean_ploidy = output, final_state = sweep(state, 2L, colSums(state), "/"),
       selected_day_trace = picked)
}

f7g_validate_segments <- function(paths, run_paths, model_cases = list()) {
  helper <- file.path(paths$oxygen_code, "util", "o2_supply_demand_map_invitro_lineage_simulation_utils.R")
  oracle <- new.env(parent = globalenv())
  sys.source(helper, envir = oracle)
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
    mean_error <- max(abs(actual$mean_ploidy - expected$mean_ploidy))
    state_error <- max(abs(actual$final_state - expected$final_state))
    selection_match <- identical(actual$selected_day_trace, expected$selected_day_trace)
    data.frame(case = z$name, mean_error = mean_error, state_error = state_error,
      selection_match = selection_match,
      minimum_selected_day = min(actual$selected_day_trace),
      maximum_selected_day = max(actual$selected_day_trace),
      external_selector_md5 = unname(tools::md5sum(helper)),
      passed = mean_error < 1e-10 && state_error < 1e-10 && selection_match)
  })
  result <- do.call(rbind, rows)
  f7ft_atomic_write_tsv(result, file.path(run_paths$run_root, "segment_selector_validation.tsv"))
  if (!all(result$passed)) stop("Fitting-selector comparison failed: ",
                              paste(result$case[!result$passed], collapse = ", "))
  result
}
