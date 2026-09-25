#!/usr/bin/env Rscript

# Side-effect-free utilities for deterministic CatA/CatB/CatC/CatU
# classification of materialized 0-1000 day mean chromosome trajectories.

o2jpc_smooth_curve <- function(value, window = 21L) {
  value <- suppressWarnings(as.numeric(value))
  window <- max(3L, as.integer(window))
  if (window %% 2L == 0L) window <- window + 1L
  if (length(value) < window) return(value)
  smoothed <- as.numeric(stats::filter(value, rep(1 / window, window), sides = 2L))
  smoothed[!is.finite(smoothed)] <- value[!is.finite(smoothed)]
  smoothed
}

o2jpc_crossing_day <- function(day, value, threshold) {
  hit <- which(is.finite(value) & value <= threshold)
  if (length(hit)) day[[hit[[1L]]]] else NA_real_
}

o2jpc_segment_sse <- function(csum, csum2, start, end) {
  n <- end - start + 1L
  total <- csum[[end + 1L]] - csum[[start]]
  total2 <- csum2[[end + 1L]] - csum2[[start]]
  max(total2 - total^2 / n, 0)
}

o2jpc_transition_bic <- function(day, value, step = 20L, min_segment = 4L) {
  keep <- !duplicated(round(day / step)) & is.finite(day) & is.finite(value)
  y <- value[keep]
  x <- day[keep]
  n <- length(y)
  if (n < 3L * min_segment) {
    return(c(one_transition_bic = NA, two_transition_bic = NA, delta_bic_two_minus_one = NA, one_change_day = NA, two_change_day_1 = NA, two_change_day_2 = NA))
  }
  csum <- c(0, cumsum(y)); csum2 <- c(0, cumsum(y^2))
  one_candidates <- min_segment:(n - min_segment)
  one_sse <- vapply(one_candidates, function(k) {
    o2jpc_segment_sse(csum, csum2, 1L, k) + o2jpc_segment_sse(csum, csum2, k + 1L, n)
  }, numeric(1L))
  one_index <- which.min(one_sse)
  one_cut <- one_candidates[[one_index]]
  best_two_sse <- Inf; best_cuts <- c(NA_integer_, NA_integer_)
  for (first in min_segment:(n - 2L * min_segment)) {
    for (second in (first + min_segment):(n - min_segment)) {
      sse <- o2jpc_segment_sse(csum, csum2, 1L, first) +
        o2jpc_segment_sse(csum, csum2, first + 1L, second) +
        o2jpc_segment_sse(csum, csum2, second + 1L, n)
      if (sse < best_two_sse) {
        best_two_sse <- sse
        best_cuts <- c(first, second)
      }
    }
  }
  one_bic <- n * log(max(one_sse[[one_index]] / n, .Machine$double.eps)) + 3L * log(n)
  two_bic <- n * log(max(best_two_sse / n, .Machine$double.eps)) + 5L * log(n)
  c(
    one_transition_bic = one_bic,
    two_transition_bic = two_bic,
    delta_bic_two_minus_one = two_bic - one_bic,
    one_change_day = x[[one_cut]],
    two_change_day_1 = x[[best_cuts[[1L]]]],
    two_change_day_2 = x[[best_cuts[[2L]]]]
  )
}

o2jpc_detect_drop_episodes <- function(
    day,
    value,
    slope_window = 21L,
    slope_threshold = 0.025,
    min_drop = 3,
    min_duration = 5,
    merge_gap = 20) {
  n <- length(day)
  if (n < 3L) return(data.frame())
  smooth <- o2jpc_smooth_curve(value, slope_window)
  half <- max(1L, floor(slope_window / 2L))
  slope <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    left <- max(1L, i - half); right <- min(n, i + half)
    elapsed <- day[[right]] - day[[left]]
    if (elapsed > 0) slope[[i]] <- (smooth[[right]] - smooth[[left]]) / elapsed
  }
  active <- is.finite(slope) & slope <= -abs(slope_threshold)
  runs <- rle(active)
  ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1L
  segments <- data.frame(start = starts[runs$values], end = ends[runs$values])
  if (!nrow(segments)) return(data.frame())
  merged <- list(); current <- segments[1L, , drop = FALSE]
  if (nrow(segments) > 1L) {
    for (i in 2:nrow(segments)) {
      gap <- day[[segments$start[[i]]]] - day[[current$end[[1L]]]]
      if (gap <= merge_gap) {
        current$end <- segments$end[[i]]
      } else {
        merged[[length(merged) + 1L]] <- current
        current <- segments[i, , drop = FALSE]
      }
    }
  }
  merged[[length(merged) + 1L]] <- current
  segments <- do.call(rbind, merged)
  rows <- lapply(seq_len(nrow(segments)), function(i) {
    start <- segments$start[[i]]; end <- segments$end[[i]]
    local <- smooth[start:end]
    data.frame(
      episode = i,
      start_index = start,
      end_index = end,
      start_day = day[[start]],
      end_day = day[[end]],
      duration = day[[end]] - day[[start]],
      drop = max(local, na.rm = TRUE) - min(local, na.rm = TRUE),
      mean_slope = mean(slope[start:end], na.rm = TRUE),
      min_slope = min(slope[start:end], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[out$drop >= min_drop & out$duration >= min_duration, , drop = FALSE]
  if (nrow(out)) out$episode <- seq_len(nrow(out))
  out
}

o2jpc_longest_plateau_between_episodes <- function(day, value, episodes, slope_limit = 0.02) {
  if (nrow(episodes) < 2L) return(c(duration = 0, slope = NA, amplitude = NA))
  best <- c(duration = 0, slope = NA, amplitude = NA)
  for (i in seq_len(nrow(episodes) - 1L)) {
    start <- episodes$end_index[[i]] + 1L
    end <- episodes$start_index[[i + 1L]] - 1L
    if (end <= start) next
    x <- day[start:end]; y <- value[start:end]
    fit_slope <- if (length(x) > 1L) unname(stats::coef(stats::lm(y ~ x))[[2L]]) else NA_real_
    duration <- max(x) - min(x)
    amplitude <- max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
    if (is.finite(fit_slope) && abs(fit_slope) <= slope_limit && duration > best[["duration"]]) {
      best <- c(duration = duration, slope = fit_slope, amplitude = amplitude)
    }
  }
  best
}

o2jpc_curve_features <- function(
    day,
    value,
    high_floor = 44,
    high_tolerance = 0.5,
    low_endpoint = 30,
    terminal_window = 50,
    slope_threshold = 0.025,
    plateau_min_days = 60,
    plateau_slope_limit = 0.02,
    bic_delta_cutoff = -10) {
  day <- suppressWarnings(as.numeric(day)); value <- suppressWarnings(as.numeric(value))
  ok <- is.finite(day) & is.finite(value)
  day <- day[ok]; value <- value[ok]
  ord <- order(day); day <- day[ord]; value <- value[ord]
  if (length(day) < 20L) {
    return(data.frame(category = "CatU", category_reason = "insufficient_points", n_points = length(day), stringsAsFactors = FALSE))
  }
  smooth <- o2jpc_smooth_curve(value, 21L)
  terminal_keep <- day >= max(day) - terminal_window
  terminal <- stats::median(smooth[terminal_keep], na.rm = TRUE)
  episodes <- o2jpc_detect_drop_episodes(day, smooth, slope_threshold = slope_threshold)
  plateau <- o2jpc_longest_plateau_between_episodes(day, smooth, episodes, plateau_slope_limit)
  transition <- o2jpc_transition_bic(day, smooth)
  diff_value <- diff(smooth) / diff(day)
  high_pass <- min(smooth, na.rm = TRUE) >= high_floor - high_tolerance
  low_pass <- terminal <= low_endpoint
  two_stage_evidence <- nrow(episodes) >= 2L && plateau[["duration"]] >= plateau_min_days
  two_stage_bic <- is.finite(transition[["delta_bic_two_minus_one"]]) && transition[["delta_bic_two_minus_one"]] <= bic_delta_cutoff
  category <- "CatU"; reason <- "does_not_meet_prespecified_rules"
  if (high_pass) {
    category <- "CatA"; reason <- "trajectory_remains_at_or_above_high_floor_with_tolerance"
  } else if (low_pass && two_stage_evidence && two_stage_bic) {
    category <- "CatC"; reason <- "low_endpoint_with_two_drop_episodes_and_sustained_plateau"
  } else if (low_pass && nrow(episodes) >= 1L && !two_stage_evidence) {
    category <- "CatB"; reason <- "low_endpoint_with_single_effective_drop_and_no_qualifying_plateau"
  }
  data.frame(
    category = category,
    category_reason = reason,
    n_points = length(day),
    day_min = min(day), day_max = max(day),
    initial_value = smooth[[1L]], terminal_median = terminal,
    minimum_value = min(smooth, na.rm = TRUE), maximum_value = max(smooth, na.rm = TRUE),
    crossing_day_44 = o2jpc_crossing_day(day, smooth, 44),
    crossing_day_35 = o2jpc_crossing_day(day, smooth, 35),
    crossing_day_30 = o2jpc_crossing_day(day, smooth, 30),
    crossing_day_25 = o2jpc_crossing_day(day, smooth, 25),
    max_negative_slope = if (length(diff_value)) min(diff_value, na.rm = TRUE) else NA_real_,
    max_negative_slope_day = if (length(diff_value)) day[[which.min(diff_value)]] else NA_real_,
    n_drop_episodes = nrow(episodes),
    first_drop_start_day = if (nrow(episodes)) episodes$start_day[[1L]] else NA_real_,
    first_drop_end_day = if (nrow(episodes)) episodes$end_day[[1L]] else NA_real_,
    first_drop_size = if (nrow(episodes)) episodes$drop[[1L]] else NA_real_,
    second_drop_start_day = if (nrow(episodes) >= 2L) episodes$start_day[[2L]] else NA_real_,
    second_drop_end_day = if (nrow(episodes) >= 2L) episodes$end_day[[2L]] else NA_real_,
    second_drop_size = if (nrow(episodes) >= 2L) episodes$drop[[2L]] else NA_real_,
    plateau_duration = plateau[["duration"]], plateau_slope = plateau[["slope"]], plateau_amplitude = plateau[["amplitude"]],
    one_transition_bic = transition[["one_transition_bic"]],
    two_transition_bic = transition[["two_transition_bic"]],
    delta_bic_two_minus_one = transition[["delta_bic_two_minus_one"]],
    one_change_day = transition[["one_change_day"]],
    two_change_day_1 = transition[["two_change_day_1"]],
    two_change_day_2 = transition[["two_change_day_2"]],
    high_floor_pass = high_pass, low_endpoint_pass = low_pass,
    two_stage_episode_pass = two_stage_evidence, two_stage_bic_pass = two_stage_bic,
    stringsAsFactors = FALSE
  )
}

o2jpc_seed_consensus <- function(cohort_categories) {
  cohort_categories <- as.character(cohort_categories)
  if (length(cohort_categories) != 2L || any(!cohort_categories %in% c("CatA", "CatB", "CatC"))) {
    return(c(category = "CatU", reason = "missing_uncertain_or_nonconcordant_cohort"))
  }
  if (length(unique(cohort_categories)) == 1L) {
    return(c(category = cohort_categories[[1L]], reason = "2N_and_4N_concordant"))
  }
  if (setequal(cohort_categories, c("CatB", "CatC"))) {
    return(c(
      category = "CatC",
      reason = "complementary_two_stage_pattern_2N_starts_on_the_4N_middle_plateau"
    ))
  }
  c(category = "CatU", reason = "2N_and_4N_disagree")
}

o2jpc_category_definition <- function(
    high_floor = 44,
    high_tolerance = 0.5,
    low_endpoint = 30,
    terminal_window = 50,
    slope_threshold = 0.025,
    plateau_min_days = 60,
    plateau_slope_limit = 0.02,
    bic_delta_cutoff = -10) {
  data.frame(
    setting = c(
      "analysis_window_days", "high_floor", "high_tolerance", "low_endpoint",
      "terminal_window_days", "rolling_slope_threshold_chr_per_day",
      "plateau_min_days", "plateau_abs_slope_limit_chr_per_day", "two_transition_bic_delta_cutoff",
      "cohort_CatA_rule", "cohort_CatB_rule", "cohort_CatC_rule",
      "seed_CatA_rule", "seed_CatB_rule", "seed_CatC_rule", "seed_CatU_rule"
    ),
    value = c(
      "0-1000", high_floor, high_tolerance, low_endpoint, terminal_window,
      slope_threshold, plateau_min_days, plateau_slope_limit, bic_delta_cutoff,
      "trajectory remains >= high_floor-high_tolerance",
      "trajectory ends <= low_endpoint with >=1 effective drop and no qualifying middle plateau",
      "trajectory ends <= low_endpoint with >=2 effective drops, sustained plateau, and two-transition BIC support",
      "both cohorts are cohort CatA",
      "both cohorts are cohort CatB",
      "both cohorts are cohort CatC, or one is CatC and the other CatB because 2N starts on the 4N middle plateau",
      "missing/edge/unclassified cohort or incompatible A-versus-low pattern"
    ),
    stringsAsFactors = FALSE
  )
}
