# Canonical in-vitro fitting and exported-table summary helpers.

ivt_weighted_mean_kary_N <- function(weights, grid_pre) {
  w <- as.numeric(weights)
  if (length(w) != length(grid_pre)) stop("weights and grid_pre length mismatch.")
  s <- sum(w)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  sum((w / s) * as.numeric(grid_pre))
}

ivt_weighted_quantile_kary_N <- function(weights, grid_pre, probs) {
  w <- as.numeric(weights)
  grid <- as.numeric(grid_pre)
  probs <- as.numeric(probs)
  if (length(w) != length(grid)) stop("weights and grid_pre length mismatch.")
  if (!length(probs)) return(numeric(0))

  keep <- is.finite(w) & is.finite(grid) & w > 0
  if (!any(keep)) return(rep(NA_real_, length(probs)))

  w <- w[keep]
  grid <- grid[keep]
  ord <- order(grid)
  w <- w[ord]
  grid <- grid[ord]
  cdf <- cumsum(w / sum(w))
  vapply(pmin(pmax(probs, 0), 1), function(p) {
    grid[which(cdf >= p)[[1]]]
  }, numeric(1))
}

ivt_log_growth_rate <- function(initial_cells, final_cells, duration_days, eps = 1e-12) {
  ni <- as.numeric(initial_cells)
  nf <- as.numeric(final_cells)
  dt <- as.numeric(duration_days)
  if (!is.finite(ni) || !is.finite(nf) || !is.finite(dt) || ni <= 0 || nf <= 0 || dt <= 0) {
    return(NA_real_)
  }
  (log(pmax(nf, eps)) - log(pmax(ni, eps))) / dt
}

.ivt_log2_population_change <- function(initial_cells, final_cells) {
  ni <- suppressWarnings(as.numeric(initial_cells))
  nf <- suppressWarnings(as.numeric(final_cells))
  if (length(ni) != 1L ||
      length(nf) != 1L ||
      !is.finite(ni) ||
      !is.finite(nf) ||
      ni <= 0 ||
      nf <= 0) {
    return(NA_real_)
  }
  log2(nf / ni)
}

.ivt_nonnegative_ratio <- function(numerator, denominator) {
  numerator <- suppressWarnings(as.numeric(numerator))
  denominator <- suppressWarnings(as.numeric(denominator))
  if (length(numerator) != 1L ||
      length(denominator) != 1L ||
      !is.finite(numerator) ||
      !is.finite(denominator) ||
      numerator < 0 ||
      denominator <= 0) {
    return(NA_real_)
  }
  numerator / denominator
}

.ivt_first_scalar <- function(..., default = NA) {
  for (value in list(...)) {
    if (is.null(value) || !length(value)) next
    candidate <- value[[1L]]
    if (length(candidate) != 1L || is.na(candidate)) next
    if (is.character(candidate) && !nzchar(candidate)) next
    return(candidate)
  }
  default
}

.ivt_legacy_lineage_label <- function(segment) {
  key <- as.character(.ivt_first_scalar(
    segment$lineage_terminal_key,
    segment$segment_id,
    default = ""
  ))
  if (!nzchar(key)) return("lineage")
  parts <- strsplit(key, "_", fixed = TRUE)[[1L]]
  if (length(parts) && all(trimws(parts) == "20.5")) "control" else "deprived"
}

.ivt_legacy_endpoint_state <- function(seg_res, endpoint_index, grid_pre = NULL) {
  if (!is.list(seg_res$sim) || is.null(seg_res$sim$live_state_obs)) return(NULL)
  raw_state <- seg_res$sim$live_state_obs
  live_state <- if (is.null(dim(raw_state))) {
    matrix(as.numeric(raw_state), nrow = 1L)
  } else {
    as.matrix(raw_state)
  }
  if (!nrow(live_state) || !ncol(live_state)) return(NULL)

  endpoint_index <- suppressWarnings(as.integer(endpoint_index))
  if (length(endpoint_index) != 1L ||
      !is.finite(endpoint_index) ||
      endpoint_index < 1L ||
      endpoint_index > nrow(live_state)) {
    endpoint_index <- nrow(live_state)
  }
  endpoint_state <- suppressWarnings(as.numeric(live_state[endpoint_index, ]))
  if (any(!is.finite(endpoint_state)) || any(endpoint_state < 0)) return(NULL)
  endpoint_total <- sum(endpoint_state)
  endpoint_frac <- if (is.finite(endpoint_total) && endpoint_total > 0) {
    endpoint_state / endpoint_total
  } else {
    rep(0, length(endpoint_state))
  }
  predicted_mean_kary_N <- NA_real_
  grid_pre <- suppressWarnings(as.numeric(grid_pre))
  if (length(grid_pre) == length(endpoint_frac) && all(is.finite(grid_pre))) {
    predicted_mean_kary_N <- ivt_weighted_mean_kary_N(
      endpoint_frac,
      grid_pre = grid_pre
    )
  }

  list(
    endpoint_index = endpoint_index,
    endpoint_state = endpoint_state,
    endpoint_live_cells = endpoint_total,
    endpoint_frac = endpoint_frac,
    predicted_mean_kary_N = predicted_mean_kary_N
  )
}

.ivt_normalize_segment_result <- function(seg_res, grid_pre = NULL) {
  if (!is.list(seg_res) || !is.list(seg_res$segment)) {
    stop("Each segment result must contain a segment list.")
  }
  seg <- seg_res$segment
  selection <- if (is.list(seg_res$selection)) seg_res$selection else list()
  segment_id <- as.character(.ivt_first_scalar(seg$segment_id, default = "legacy-segment"))
  cohort <- as.character(.ivt_first_scalar(seg$cohort, default = "cohort"))
  lineage_label <- as.character(.ivt_first_scalar(
    seg$lineage_label,
    .ivt_legacy_lineage_label(seg),
    default = "lineage"
  ))
  lineage_id <- as.character(.ivt_first_scalar(seg$lineage_id, lineage_label, default = "lineage"))
  lineage_group <- as.character(.ivt_first_scalar(
    seg$lineage_group,
    if (lineage_label %in% c("C", "control")) "control" else "deprived",
    default = "lineage"
  ))
  scenario_id <- as.character(.ivt_first_scalar(
    seg$scenario_id,
    paste(cohort, lineage_id, sep = "-"),
    default = paste(cohort, "lineage", sep = "-")
  ))
  passage_index <- as.integer(.ivt_first_scalar(seg$passage_index, default = 1L))
  lineage_passage_index <- as.integer(.ivt_first_scalar(
    seg$lineage_passage_index,
    if (!is.null(seg$depth) && length(seg$depth) && is.finite(as.numeric(seg$depth[[1L]]))) {
      as.integer(seg$depth[[1L]]) + 1L
    } else {
      passage_index
    },
    default = passage_index
  ))
  obs_days <- suppressWarnings(as.numeric(seg$obs_days_local))
  obs_days <- obs_days[is.finite(obs_days)]
  observed_endpoint <- if (length(obs_days)) max(obs_days) else NA_real_
  passage_duration <- as.numeric(.ivt_first_scalar(
    seg$passage_duration,
    seg$duration_days,
    observed_endpoint,
    default = NA_real_
  ))
  endpoint_day <- as.numeric(.ivt_first_scalar(
    selection$endpoint_day,
    seg$endpoint_day,
    passage_duration,
    observed_endpoint,
    default = NA_real_
  ))
  legacy_selected_day <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$selected_day,
    default = NA_real_
  )))
  legacy_selected_index <- suppressWarnings(as.integer(.ivt_first_scalar(
    selection$selected_index,
    default = NA_integer_
  )))
  legacy_selected_live_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$selected_live_cells,
    default = NA_real_
  )))
  data_ids <- as.character(seg$data_ids)
  data_ids <- data_ids[!is.na(data_ids) & nzchar(data_ids)]
  passage_id <- as.character(.ivt_first_scalar(
    seg$passage_id,
    if (length(data_ids)) paste(data_ids, collapse = ";") else segment_id,
    default = segment_id
  ))
  endpoint_index <- as.integer(.ivt_first_scalar(
    selection$endpoint_index,
    selection$selected_index,
    if (length(obs_days)) length(obs_days) else NA_integer_,
    default = NA_integer_
  ))
  legacy_endpoint <- if (is.null(selection$endpoint_state) &&
                         is.null(selection$selected_state)) {
    .ivt_legacy_endpoint_state(
      seg_res = seg_res,
      endpoint_index = endpoint_index,
      grid_pre = grid_pre
    )
  } else {
    NULL
  }
  endpoint_was_derived <- !is.null(legacy_endpoint)
  if (endpoint_was_derived) {
    endpoint_index <- legacy_endpoint$endpoint_index
  }
  endpoint_live_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$endpoint_live_cells,
    if (endpoint_was_derived) legacy_endpoint$endpoint_live_cells else NULL,
    if (is.list(seg_res$sim)) tail(seg_res$sim$Ntot_live_obs, 1L) else NULL,
    default = NA_real_
  )))
  selected_day <- as.numeric(.ivt_first_scalar(
    legacy_selected_day,
    endpoint_day,
    default = NA_real_
  ))
  closest_day_diagnostic <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$closest_day_diagnostic,
    NA_real_,
    default = NA_real_
  )))

  seg$segment_id <- segment_id
  seg$parent_segment_id <- as.character(.ivt_first_scalar(seg$parent_segment_id, default = NA_character_))
  seg$cohort <- cohort
  seg$lineage_id <- lineage_id
  seg$lineage_group <- lineage_group
  seg$lineage_label <- lineage_label
  seg$scenario_id <- scenario_id
  seg$lineage_terminal_key <- as.character(.ivt_first_scalar(
    seg$lineage_terminal_key,
    segment_id,
    default = segment_id
  ))
  seg$passage_index <- passage_index
  seg$lineage_passage_index <- lineage_passage_index
  seg$passage_id <- passage_id
  seg$duration_days <- as.numeric(.ivt_first_scalar(seg$duration_days, passage_duration, default = passage_duration))
  seg$passage_duration <- passage_duration
  seg$observed_passage_day <- suppressWarnings(as.numeric(.ivt_first_scalar(
    seg$observed_passage_day,
    passage_duration,
    default = NA_real_
  )))
  seg$last_observation_day <- suppressWarnings(as.numeric(.ivt_first_scalar(
    seg$last_observation_day,
    seg$observed_passage_day,
    passage_duration,
    default = NA_real_
  )))
  seg$search_horizon_day <- suppressWarnings(as.numeric(.ivt_first_scalar(
    seg$search_horizon_day,
    seg$duration_days,
    observed_endpoint,
    default = NA_real_
  )))
  seg$passage_time_tolerance_days <- suppressWarnings(as.numeric(.ivt_first_scalar(
    seg$passage_time_tolerance_days,
    if (is.finite(seg$search_horizon_day) && is.finite(seg$observed_passage_day)) {
      max(seg$search_horizon_day - seg$observed_passage_day, 0)
    } else {
      NULL
    },
    default = NA_real_
  )))
  seg$endpoint_day <- selected_day
  seg$oxygen_pct <- suppressWarnings(as.numeric(.ivt_first_scalar(seg$oxygen_pct, default = NA_real_)))
  fallback_n_days <- if (is.finite(endpoint_index) && endpoint_index > 0L) endpoint_index else 1L
  seg$obs_days_local <- if (length(obs_days)) obs_days else seq_len(fallback_n_days) - 1L

  selection$endpoint_index <- endpoint_index
  selection$endpoint_day <- selected_day
  selection$endpoint_live_cells <- endpoint_live_cells
  if (endpoint_was_derived) {
    selection$endpoint_state <- legacy_endpoint$endpoint_state
    selection$endpoint_frac <- legacy_endpoint$endpoint_frac
    selection$selected_index <- endpoint_index
    selection$selected_live_cells <- endpoint_live_cells
    selection$selected_frac <- legacy_endpoint$endpoint_frac
    if (is.finite(legacy_endpoint$predicted_mean_kary_N)) {
      selection$predicted_mean_kary_N <- legacy_endpoint$predicted_mean_kary_N
    }
  }
  selection$selected_day <- selected_day
  selection$selected_index <- as.integer(.ivt_first_scalar(
    selection$selected_index,
    endpoint_index,
    default = NA_integer_
  ))
  selection$selected_live_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$selected_live_cells,
    endpoint_live_cells,
    default = NA_real_
  )))
  selection$search_horizon_day <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$search_horizon_day,
    seg$search_horizon_day,
    default = NA_real_
  )))
  selection$search_horizon_live_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$search_horizon_live_cells,
    if (is.list(seg_res$sim)) tail(seg_res$sim$Ntot_live_obs, 1L) else NULL,
    default = NA_real_
  )))
  selection$max_live_cells_in_search <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$max_live_cells_in_search,
    if (is.list(seg_res$sim)) max(seg_res$sim$Ntot_live_obs, na.rm = TRUE) else NULL,
    default = NA_real_
  )))
  selection$closest_index_diagnostic <- as.integer(.ivt_first_scalar(
    selection$closest_index_diagnostic,
    if (is.finite(closest_day_diagnostic)) legacy_selected_index else NULL,
    default = NA_integer_
  ))
  selection$closest_day_diagnostic <- closest_day_diagnostic
  selection$closest_live_cells_diagnostic <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$closest_live_cells_diagnostic,
    if (is.finite(closest_day_diagnostic)) legacy_selected_live_cells else NULL,
    default = NA_real_
  )))
  selection$target_live_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$target_live_cells,
    seg$final_cells,
    seg$protocol_threshold_cells,
    default = NA_real_
  )))
  selection$protocol_threshold_cells <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$protocol_threshold_cells,
      seg$protocol_threshold_cells,
      selection$target_live_cells,
      default = NA_real_
    )
  ))
  selection$protocol_threshold_source <- as.character(.ivt_first_scalar(
    selection$protocol_threshold_source,
    seg$protocol_threshold_source,
    default = "legacy_observed_final_cells"
  ))
  selection$observed_final_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$observed_final_cells,
    seg$final_cells,
    default = NA_real_
  )))
  selection$observed_final_target_cells <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$observed_final_target_cells,
      selection$observed_final_cells,
      selection$target_live_cells,
      default = NA_real_
    )
  ))
  selection$passage_recorded <- as.logical(.ivt_first_scalar(
    selection$passage_recorded,
    length(data_ids) > 0L,
    default = FALSE
  ))
  selection$passage_executed <- as.logical(.ivt_first_scalar(
    selection$passage_executed,
    selection$passage_recorded,
    default = FALSE
  ))
  selection$passage_failure_reason <- as.character(.ivt_first_scalar(
    selection$passage_failure_reason,
    default = NA_character_
  ))
  selection$reseed_mode <- as.character(.ivt_first_scalar(
    selection$reseed_mode,
    default = "legacy_unavailable"
  ))
  selection$available_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$available_cells,
    endpoint_live_cells,
    default = NA_real_
  )))
  selection$required_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$required_cells,
    default = NA_real_
  )))
  selection$supply_ratio <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$supply_ratio,
    default = NA_real_
  )))
  selection$boundary_scale <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$boundary_scale,
    default = NA_real_
  )))
  selection$cell_number_before <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$cell_number_before,
    endpoint_live_cells,
    default = NA_real_
  )))
  selection$cell_number_after <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$cell_number_after,
    default = NA_real_
  )))
  threshold_target_cells <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$threshold_target_cells,
    selection$target_live_cells,
    seg$final_cells,
    selection$protocol_threshold_cells,
    default = NA_real_
  )))
  if (!is.finite(threshold_target_cells) || threshold_target_cells <= 0) {
    threshold_target_cells <- NA_real_
  }
  threshold_target_source <- as.character(.ivt_first_scalar(
    selection$threshold_target_source,
    if (is.finite(threshold_target_cells)) "legacy_observed_final_cells" else "missing",
    default = "missing"
  ))
  threshold_crossing <- NULL
  if (exists("ivt_first_threshold_crossing", mode = "function", inherits = TRUE) &&
      is.list(seg_res$sim) &&
      length(seg_res$sim$Ntot_live_obs) == length(seg$obs_days_local)) {
    threshold_crossing <- ivt_first_threshold_crossing(
      days = seg$obs_days_local,
      live_cells = seg_res$sim$Ntot_live_obs,
      threshold_target_cells = threshold_target_cells
    )
  }
  selection$threshold_target_cells <- threshold_target_cells
  selection$effective_threshold_cells <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$effective_threshold_cells,
      threshold_target_cells,
      default = NA_real_
    )
  ))
  selection$threshold_target_source <- threshold_target_source
  selection$threshold_reached_by_endpoint <- as.logical(.ivt_first_scalar(
    selection$threshold_reached_by_endpoint,
    if (is.list(threshold_crossing)) {
      threshold_crossing$threshold_reached_by_endpoint
    } else {
      FALSE
    },
    default = FALSE
  ))
  selection$predicted_threshold_crossing_day <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$predicted_threshold_crossing_day,
      if (is.list(threshold_crossing)) {
        threshold_crossing$predicted_threshold_crossing_day
      } else {
        NULL
      },
      default = NA_real_
    )
  ))
  selection$threshold_time_grid_resolution_days <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$threshold_time_grid_resolution_days,
      if (is.list(threshold_crossing)) {
        threshold_crossing$threshold_time_grid_resolution_days
      } else {
        NULL
      },
      default = NA_real_
    )
  ))
  selection$threshold_crossing_interval_width_days <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$threshold_crossing_interval_width_days,
      if (is.list(threshold_crossing)) {
        threshold_crossing$threshold_crossing_interval_width_days
      } else {
        NULL
      },
      default = NA_real_
    )
  ))
  selection$observed_passage_day <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$observed_passage_day,
    seg$observed_passage_day,
    passage_duration,
    default = NA_real_
  )))
  selection$last_observation_day <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$last_observation_day,
    seg$last_observation_day,
    selection$observed_passage_day,
    passage_duration,
    default = NA_real_
  )))
  matched_last_observation_idx <- which(
    is.finite(seg$obs_days_local) &
      is.finite(selection$last_observation_day) &
      abs(seg$obs_days_local - selection$last_observation_day) <= 1e-10
  )
  fallback_last_observation_idx <- if (length(matched_last_observation_idx) == 1L) {
    matched_last_observation_idx[[1L]]
  } else {
    NA_integer_
  }
  selection$last_observation_index <- suppressWarnings(as.integer(.ivt_first_scalar(
    selection$last_observation_index,
    fallback_last_observation_idx,
    default = NA_integer_
  )))
  fallback_observation_live_cells <- if (
    is.list(seg_res$sim) &&
      is.finite(selection$last_observation_index) &&
      selection$last_observation_index >= 1L &&
      selection$last_observation_index <= length(seg_res$sim$Ntot_live_obs)
  ) {
    seg_res$sim$Ntot_live_obs[[selection$last_observation_index]]
  } else {
    NA_real_
  }
  selection$predicted_live_cells_at_observation <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$predicted_live_cells_at_observation,
      fallback_observation_live_cells,
      default = NA_real_
    )
  ))
  selection$selected_day_after_last_observation <- as.logical(.ivt_first_scalar(
    selection$selected_day_after_last_observation,
    if (is.finite(selected_day) && is.finite(selection$last_observation_day)) {
      selected_day >= selection$last_observation_day - 1e-10
    } else {
      NULL
    },
    default = FALSE
  ))
  selection$passage_time_tolerance_days <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$passage_time_tolerance_days,
    seg$passage_time_tolerance_days,
    default = NA_real_
  )))
  selection$passage_time_residual_days <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$passage_time_residual_days,
    if (is.finite(selected_day) && is.finite(selection$observed_passage_day)) {
      selected_day - selection$observed_passage_day
    } else {
      NULL
    },
    default = NA_real_
  )))
  selection$passage_time_within_tolerance <- as.logical(.ivt_first_scalar(
    selection$passage_time_within_tolerance,
    if (is.finite(selection$passage_time_residual_days) &&
        is.finite(selection$passage_time_tolerance_days)) {
      abs(selection$passage_time_residual_days) <=
        selection$passage_time_tolerance_days + 1e-10
    } else {
      NULL
    },
    default = FALSE
  ))
  selection$threshold_time_residual_days <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$threshold_time_residual_days,
      if (is.finite(selection$predicted_threshold_crossing_day) &&
          is.finite(selection$observed_passage_day)) {
        selection$predicted_threshold_crossing_day -
          selection$observed_passage_day
      } else {
        NULL
      },
      default = NA_real_
    )
  ))
  selection$endpoint_cell_count_residual <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$endpoint_cell_count_residual,
      if (is.finite(endpoint_live_cells) &&
          is.finite(threshold_target_cells)) {
        endpoint_live_cells - threshold_target_cells
      } else {
        NULL
      },
      default = NA_real_
    )
  ))
  selection$cell_count_overshoot <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$cell_count_overshoot,
    selection$endpoint_cell_count_residual,
    default = NA_real_
  )))
  selection$protocol_threshold_overshoot <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$protocol_threshold_overshoot,
      if (is.finite(endpoint_live_cells) &&
          is.finite(selection$protocol_threshold_cells)) {
        endpoint_live_cells - selection$protocol_threshold_cells
      } else {
        NULL
      },
      default = NA_real_
    )
  ))
  selection$protocol_threshold_reached_at_selected <- as.logical(
    .ivt_first_scalar(
      selection$protocol_threshold_reached_at_selected,
      if (is.finite(selection$selected_live_cells) &&
          is.finite(selection$protocol_threshold_cells)) {
        selection$selected_live_cells >= selection$protocol_threshold_cells
      } else {
        NULL
      },
      default = FALSE
    )
  )
  selection$protocol_threshold_reached_at_observation <- as.logical(
    .ivt_first_scalar(
      selection$protocol_threshold_reached_at_observation,
      if (is.finite(selection$predicted_live_cells_at_observation) &&
          is.finite(selection$protocol_threshold_cells)) {
        selection$predicted_live_cells_at_observation >=
          selection$protocol_threshold_cells
      } else {
        NULL
      },
      default = FALSE
    )
  )
  selection$observed_final_cell_count_residual <- suppressWarnings(as.numeric(
    .ivt_first_scalar(
      selection$observed_final_cell_count_residual,
      if (is.finite(selection$predicted_live_cells_at_observation) &&
          is.finite(selection$observed_final_cells)) {
        selection$predicted_live_cells_at_observation -
          selection$observed_final_cells
      } else {
        NULL
      },
      default = NA_real_
    )
  ))

  seg_res$segment <- seg
  seg_res$selection <- selection
  seg_res
}

ivt_collect_lineage_summary <- function(run, fit_data) {
  out <- do.call(rbind, lapply(run$segment_results, function(seg_res) {
    seg_res <- .ivt_normalize_segment_result(seg_res, grid_pre = run$grid_pre)
    seg <- seg_res$segment
    pred_mean_kary_N <- as.numeric(seg_res$selection$predicted_mean_kary_N)
    sim_live <- as.numeric(seg_res$sim$Ntot_live_obs)
    selected_idx <- as.integer(seg_res$selection$selected_index)
    if (!is.finite(selected_idx) || selected_idx < 1L || selected_idx > length(sim_live)) {
      stop("Selected passage index is unavailable for segment ", seg$segment_id, ".")
    }
    last_observation_idx <- as.integer(seg_res$selection$last_observation_index)
    if (!is.finite(last_observation_idx) ||
        last_observation_idx < 1L ||
        last_observation_idx > length(sim_live)) {
      stop("Last experimental observation index is unavailable for segment ", seg$segment_id, ".")
    }
    sim_terminal_scalar <- function(terminal_name, obs_name) {
      obs_value <- suppressWarnings(as.numeric(seg_res$sim[[obs_name]]))
      if (length(obs_value) >= selected_idx &&
          is.finite(obs_value[[selected_idx]])) {
        return(obs_value[[selected_idx]])
      }
      terminal_value <- suppressWarnings(as.numeric(seg_res$sim[[terminal_name]]))
      if (selected_idx == length(sim_live) &&
          length(terminal_value) == 1L && is.finite(terminal_value)) {
        return(terminal_value)
      }
      NA_real_
    }
    predicted_gross_division_events <- sim_terminal_scalar(
      "cumulative_gross_divisions_terminal",
      "cumulative_gross_divisions_obs"
    )
    predicted_cumulative_hypoxia_deaths <- sim_terminal_scalar(
      "cumulative_hypoxia_deaths_terminal",
      "cumulative_hypoxia_deaths_obs"
    )
    predicted_cumulative_dead_buffer_inflow <- sim_terminal_scalar(
      "cumulative_dead_buffer_inflow_terminal",
      "cumulative_dead_buffer_inflow_obs"
    )
    predicted_cumulative_nonlive_inflow <- sim_terminal_scalar(
      "cumulative_nonlive_inflow_terminal",
      "cumulative_nonlive_inflow_obs"
    )
    if (!is.finite(predicted_cumulative_nonlive_inflow) &&
        is.finite(predicted_cumulative_hypoxia_deaths) &&
        is.finite(predicted_cumulative_dead_buffer_inflow)) {
      predicted_cumulative_nonlive_inflow <-
        predicted_cumulative_hypoxia_deaths +
        predicted_cumulative_dead_buffer_inflow
    }
    predicted_initial_cells <- sim_live[[1]]
    predicted_passage_live_cells <- sim_live[[selected_idx]]
    predicted_live_cells_at_observation <- sim_live[[last_observation_idx]]
    predicted_net_gain <- predicted_passage_live_cells - predicted_initial_cells
    predicted_net_population_doublings <- .ivt_log2_population_change(
      predicted_initial_cells,
      predicted_live_cells_at_observation
    )
    pred_growth <- ivt_log_growth_rate(
      initial_cells = predicted_initial_cells,
      final_cells = predicted_live_cells_at_observation,
      duration_days = seg_res$selection$last_observation_day
    )
    do.call(rbind, lapply(seg$data_ids, function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      endpoint_local_obs <- .ivt_endpoint_local_observed_growth(
        fit_data[[pid]]
      )
      observed_net_population_doublings <- .ivt_log2_population_change(
        obs$initial_cells,
        obs$final_cells
      )
      observed_minimum_division_events <- if (
        is.finite(obs$initial_cells) &&
          is.finite(obs$final_cells)
      ) {
        max(obs$final_cells - obs$initial_cells, 0)
      } else {
        NA_real_
      }
      data.frame(
        segment_id = seg$segment_id,
        parent_segment_id = if (is.null(seg$parent_segment_id)) NA_character_ else as.character(seg$parent_segment_id),
        cohort = seg$cohort,
        lineage_id = seg$lineage_id,
        lineage_group = seg$lineage_group,
        lineage_label = seg$lineage_label,
        scenario_id = seg$scenario_id,
        lineage_terminal_key = seg$lineage_terminal_key,
        passage_index = seg$passage_index,
        lineage_passage_index = seg$lineage_passage_index,
        oxygen_pct = seg$oxygen_pct,
        duration_days = seg$duration_days,
        passage_duration = seg$passage_duration,
        observed_passage_day = seg_res$selection$observed_passage_day,
        last_observation_day = seg_res$selection$last_observation_day,
        search_horizon_day = seg_res$selection$search_horizon_day,
        search_horizon_live_cells =
          seg_res$selection$search_horizon_live_cells,
        max_live_cells_in_search = seg_res$selection$max_live_cells_in_search,
        passage_time_tolerance_days =
          seg_res$selection$passage_time_tolerance_days,
        endpoint_day = seg_res$selection$selected_day,
        initial_cells = seg$initial_cells,
        final_cells = seg$final_cells,
        selected_day = seg_res$selection$selected_day,
        selected_live_cells = seg_res$selection$selected_live_cells,
        selected_day_after_last_observation =
          seg_res$selection$selected_day_after_last_observation,
        passage_time_residual_days =
          seg_res$selection$passage_time_residual_days,
        passage_time_within_tolerance =
          seg_res$selection$passage_time_within_tolerance,
        closest_day_diagnostic = seg_res$selection$closest_day_diagnostic,
        closest_live_cells_diagnostic = seg_res$selection$closest_live_cells_diagnostic,
        target_live_cells = seg_res$selection$target_live_cells,
        protocol_threshold_cells =
          seg_res$selection$protocol_threshold_cells,
        protocol_threshold_source =
          seg_res$selection$protocol_threshold_source,
        protocol_culture_vessel = as.character(.ivt_first_scalar(
          seg$protocol_culture_vessel,
          default = NA_character_
        )),
        protocol_growth_area_cm2 = suppressWarnings(as.numeric(.ivt_first_scalar(
          seg$protocol_growth_area_cm2,
          default = NA_real_
        ))),
        protocol_target_confluence = suppressWarnings(as.numeric(.ivt_first_scalar(
          seg$protocol_target_confluence,
          default = NA_real_
        ))),
        threshold_target_cells =
          seg_res$selection$threshold_target_cells,
        effective_threshold_cells =
          seg_res$selection$effective_threshold_cells,
        threshold_target_source =
          seg_res$selection$threshold_target_source,
        threshold_reached_by_endpoint =
          seg_res$selection$threshold_reached_by_endpoint,
        predicted_threshold_crossing_day =
          seg_res$selection$predicted_threshold_crossing_day,
        threshold_time_residual_days =
          seg_res$selection$threshold_time_residual_days,
        endpoint_cell_count_residual =
          seg_res$selection$endpoint_cell_count_residual,
        cell_count_overshoot = seg_res$selection$cell_count_overshoot,
        protocol_threshold_overshoot =
          seg_res$selection$protocol_threshold_overshoot,
        protocol_threshold_reached_at_selected =
          seg_res$selection$protocol_threshold_reached_at_selected,
        protocol_threshold_reached_at_observation =
          seg_res$selection$protocol_threshold_reached_at_observation,
        observed_final_cell_count_residual =
          seg_res$selection$observed_final_cell_count_residual,
        threshold_time_grid_resolution_days =
          seg_res$selection$threshold_time_grid_resolution_days,
        threshold_crossing_interval_width_days =
          seg_res$selection$threshold_crossing_interval_width_days,
        passage_id = pid,
        predicted_initial_cells = predicted_initial_cells,
        predicted_final_cells = predicted_live_cells_at_observation,
        predicted_passage_live_cells = predicted_passage_live_cells,
        observed_initial_cells = obs$initial_cells,
        observed_final_cells = obs$final_cells,
        predicted_live_cells = predicted_live_cells_at_observation,
        predicted_live_cells_at_observation =
          predicted_live_cells_at_observation,
        observed_live_cells_at_observation = obs$final_cells,
        measurement_day_cell_count_residual = if (
          is.finite(predicted_live_cells_at_observation) &&
            is.finite(obs$final_cells)
        ) {
          predicted_live_cells_at_observation - obs$final_cells
        } else {
          NA_real_
        },
        log_live_cell_residual = if (
          is.finite(predicted_live_cells_at_observation) &&
            predicted_live_cells_at_observation > 0 &&
            is.finite(obs$final_cells) &&
            obs$final_cells > 0
        ) {
          log(predicted_live_cells_at_observation / obs$final_cells)
        } else {
          NA_real_
        },
        passage_recorded = seg_res$selection$passage_recorded,
        passage_executed = seg_res$selection$passage_executed,
        passage_failure_reason = seg_res$selection$passage_failure_reason,
        reseed_mode = seg_res$selection$reseed_mode,
        insufficient_boundary = identical(
          as.character(seg_res$selection$reseed_mode),
          "carry_forward_insufficient"
        ),
        boundary_feasible = if (identical(
          as.character(seg_res$selection$reseed_mode),
          "terminal_no_reseed"
        )) {
          NA
        } else {
          !identical(
            as.character(seg_res$selection$reseed_mode),
            "carry_forward_insufficient"
          )
        },
        available_cells = seg_res$selection$available_cells,
        required_cells = seg_res$selection$required_cells,
        supply_ratio = seg_res$selection$supply_ratio,
        boundary_scale = seg_res$selection$boundary_scale,
        cell_number_before = seg_res$selection$cell_number_before,
        cell_number_after = seg_res$selection$cell_number_after,
        predicted_growth = pred_growth,
        predicted_growth_rate = pred_growth,
        predicted_endpoint_instantaneous_net_growth_rate =
          suppressWarnings(as.numeric(.ivt_first_scalar(
            seg_res$selection$predicted_endpoint_instantaneous_net_growth_rate,
            default = NA_real_
          ))),
        predicted_endpoint_division_linked_live_rate =
          suppressWarnings(as.numeric(.ivt_first_scalar(
            seg_res$selection$predicted_endpoint_division_linked_live_rate,
            default = NA_real_
          ))),
        predicted_endpoint_hypoxia_death_rate =
          suppressWarnings(as.numeric(.ivt_first_scalar(
            seg_res$selection$predicted_endpoint_hypoxia_death_rate,
            default = NA_real_
          ))),
        predicted_endpoint_crowding_multiplier =
          suppressWarnings(as.numeric(.ivt_first_scalar(
            seg_res$selection$predicted_endpoint_crowding_multiplier,
            default = NA_real_
          ))),
        observed_endpoint_local_net_growth_rate = endpoint_local_obs$rate,
        observed_endpoint_local_interval_start_day =
          endpoint_local_obs$interval_start_day,
        observed_endpoint_local_interval_end_day =
          endpoint_local_obs$interval_end_day,
        observed_endpoint_local_n_positive_timepoints =
          endpoint_local_obs$n_positive_timepoints,
        observed_net_population_doublings =
          observed_net_population_doublings,
        predicted_net_population_doublings =
          predicted_net_population_doublings,
        observed_minimum_division_events =
          observed_minimum_division_events,
        predicted_gross_division_events =
          predicted_gross_division_events,
        predicted_cumulative_hypoxia_deaths =
          predicted_cumulative_hypoxia_deaths,
        predicted_cumulative_dead_buffer_inflow =
          predicted_cumulative_dead_buffer_inflow,
        predicted_cumulative_nonlive_inflow =
          predicted_cumulative_nonlive_inflow,
        predicted_divisions_per_initial_cell = .ivt_nonnegative_ratio(
          predicted_gross_division_events,
          predicted_initial_cells
        ),
        predicted_nonlive_inflow_to_division_ratio =
          .ivt_nonnegative_ratio(
            predicted_cumulative_nonlive_inflow,
            predicted_gross_division_events
          ),
        predicted_gross_division_to_net_gain_ratio = if (
          is.finite(predicted_net_gain) &&
            predicted_net_gain > 0
        ) {
          predicted_gross_division_events / predicted_net_gain
        } else {
          NA_real_
        },
        dead_buffer_inflow_definition = as.character(.ivt_first_scalar(
          seg_res$sim$cumulative_dead_buffer_inflow_definition,
          default = paste(
            "missegregation-linked nonviable daughters plus",
            "grid boundary-routed loss"
          )
        )),
        predicted_mean_kary_N = pred_mean_kary_N,
        observed_growth = obs$observed_growth,
        observed_mean_kary_N = obs$observed_mean_kary_N,
        observed_n_kary = obs$observed_n_kary,
        observed_n_flow = obs$observed_n_flow,
        stringsAsFactors = FALSE
      )
    }))
  }))
  if (!is.null(out) && nrow(out) > 0L) {
    out <- out[order(
      match(out$cohort, c("2N", "4N")),
      match(out$lineage_id, c("C", "O1", "O2")),
      out$lineage_passage_index
    ), , drop = FALSE]
    out$cumulative_time <- ave(
      out$selected_day,
      out$scenario_id,
      FUN = cumsum
    )
    out$cumulative_experimental_time <- ave(
      out$passage_duration,
      out$scenario_id,
      FUN = cumsum
    )
    out$cumulative_observed_net_population_doublings <- ave(
      out$observed_net_population_doublings,
      out$scenario_id,
      FUN = cumsum
    )
    out$cumulative_predicted_net_population_doublings <- ave(
      out$predicted_net_population_doublings,
      out$scenario_id,
      FUN = cumsum
    )
    out$cumulative_gross_divisions <- ave(
      out$predicted_gross_division_events,
      out$scenario_id,
      FUN = cumsum
    )
    out$cumulative_hypoxia_deaths <- ave(
      out$predicted_cumulative_hypoxia_deaths,
      out$scenario_id,
      FUN = cumsum
    )
    out$cumulative_dead_buffer_inflow <- ave(
      out$predicted_cumulative_dead_buffer_inflow,
      out$scenario_id,
      FUN = cumsum
    )
    out$cumulative_nonlive_inflow <- ave(
      out$predicted_cumulative_nonlive_inflow,
      out$scenario_id,
      FUN = cumsum
    )
    n_insufficient_boundaries <- sum(
      out$reseed_mode == "carry_forward_insufficient",
      na.rm = TRUE
    )
    out$n_insufficient_boundaries <- n_insufficient_boundaries
    out$all_passage_boundaries_feasible <-
      n_insufficient_boundaries == 0L
    out$protocol_feasibility_status <- if (
      n_insufficient_boundaries == 0L
    ) {
      "PASS"
    } else {
      "FAIL"
    }
    rownames(out) <- NULL
  }
  out
}

ivt_collect_observed_kary_summary <- function(run, fit_data) {
  out <- lapply(run$segment_results, function(seg_res) {
    seg_res <- .ivt_normalize_segment_result(seg_res, grid_pre = run$grid_pre)
    seg <- seg_res$segment
    dplyr::bind_rows(lapply(.ivt_segment_endpoint_data_ids(run, seg), function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_kary <- obs$observed_kary
      observed_kary <- observed_kary[is.finite(observed_kary)]
      if (!length(observed_kary)) return(NULL)
      data.frame(
        segment_id = seg$segment_id,
        cohort = seg$cohort,
        lineage_id = seg$lineage_id,
        lineage_group = seg$lineage_group,
        lineage_label = seg$lineage_label,
        scenario_id = seg$scenario_id,
        lineage_terminal_key = seg$lineage_terminal_key,
        passage_index = seg$passage_index,
        lineage_passage_index = seg$lineage_passage_index,
        oxygen_pct = seg$oxygen_pct,
        passage_id = pid,
        cell_index = seq_along(observed_kary),
        observed_kary_N = observed_kary,
        stringsAsFactors = FALSE
      )
    }))
  })
  initial_out <- lapply(run$initial_observations, function(record) {
    obs <- record$observed
    observed_kary <- obs$observed_kary
    observed_kary <- observed_kary[is.finite(observed_kary)]
    if (!length(observed_kary)) return(NULL)
    data.frame(
      segment_id = record$segment_id,
      cohort = record$cohort,
      lineage_id = record$lineage_id,
      lineage_group = record$lineage_group,
      lineage_label = record$lineage_id,
      scenario_id = record$scenario_id,
      lineage_terminal_key = record$scenario_id,
      passage_index = 0L,
      lineage_passage_index = 0L,
      oxygen_pct = record$oxygen_pct,
      passage_id = record$passage_id,
      cell_index = seq_along(observed_kary),
      observed_kary_N = observed_kary,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(c(out, initial_out))
}

ivt_collect_observed_flow_summary <- function(run, fit_data) {
  out <- lapply(run$segment_results, function(seg_res) {
    seg_res <- .ivt_normalize_segment_result(seg_res, grid_pre = run$grid_pre)
    seg <- seg_res$segment
    dplyr::bind_rows(lapply(.ivt_segment_endpoint_data_ids(run, seg), function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_flow <- obs$observed_flow
      if (is.null(observed_flow) || !nrow(observed_flow)) return(NULL)
      data.frame(
        segment_id = seg$segment_id,
        cohort = seg$cohort,
        lineage_id = seg$lineage_id,
        lineage_group = seg$lineage_group,
        lineage_label = seg$lineage_label,
        scenario_id = seg$scenario_id,
        lineage_terminal_key = seg$lineage_terminal_key,
        passage_index = seg$passage_index,
        lineage_passage_index = seg$lineage_passage_index,
        oxygen_pct = seg$oxygen_pct,
        passage_id = pid,
        sample_name = obs$observed_flow_sample_name,
        grid_index = observed_flow$grid_index,
        ploidy = observed_flow$ploidy,
        observed_density = observed_flow$density,
        observed_log_density = observed_flow$log_density,
        stringsAsFactors = FALSE
      )
    }))
  })
  initial_out <- lapply(run$initial_observations, function(record) {
    obs <- record$observed
    observed_flow <- obs$observed_flow
    if (is.null(observed_flow) || !nrow(observed_flow)) return(NULL)
    data.frame(
      segment_id = record$segment_id,
      cohort = record$cohort,
      lineage_id = record$lineage_id,
      lineage_group = record$lineage_group,
      lineage_label = record$lineage_id,
      scenario_id = record$scenario_id,
      lineage_terminal_key = record$scenario_id,
      passage_index = 0L,
      lineage_passage_index = 0L,
      oxygen_pct = record$oxygen_pct,
      passage_id = record$passage_id,
      sample_name = obs$observed_flow_sample_name,
      grid_index = observed_flow$grid_index,
      ploidy = observed_flow$ploidy,
      observed_density = observed_flow$density,
      observed_log_density = observed_flow$log_density,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(c(out, initial_out))
}

ivt_collect_distribution_summary <- function(run) {
  do.call(rbind, lapply(run$segment_results, function(seg_res) {
    seg_res <- .ivt_normalize_segment_result(seg_res, grid_pre = run$grid_pre)
    frac <- as.numeric(seg_res$selection$selected_frac)
    data.frame(
      segment_id = seg_res$segment$segment_id,
      cohort = seg_res$segment$cohort,
      lineage_id = seg_res$segment$lineage_id,
      lineage_group = seg_res$segment$lineage_group,
      lineage_label = seg_res$segment$lineage_label,
      scenario_id = seg_res$segment$scenario_id,
      lineage_terminal_key = seg_res$segment$lineage_terminal_key,
      passage_index = seg_res$segment$passage_index,
      lineage_passage_index = seg_res$segment$lineage_passage_index,
      oxygen_pct = seg_res$segment$oxygen_pct,
      passage_id = seg_res$segment$passage_id,
      passage_duration = seg_res$segment$passage_duration,
      endpoint_day = seg_res$selection$endpoint_day,
      selected_day = seg_res$selection$selected_day,
      closest_day_diagnostic = seg_res$selection$closest_day_diagnostic,
      N = run$grid_pre,
      fraction = frac,
      stringsAsFactors = FALSE
    )
  }))
}

ivt_collect_distribution_quantiles <- function(run, probs) {
  probs <- as.numeric(probs)
  out <- lapply(run$segment_results, function(seg_res) {
    seg_res <- .ivt_normalize_segment_result(seg_res, grid_pre = run$grid_pre)
    frac <- as.numeric(seg_res$selection$selected_frac)
    quantiles <- ivt_weighted_quantile_kary_N(frac, run$grid_pre, probs)
    data.frame(
      segment_id = seg_res$segment$segment_id,
      cohort = seg_res$segment$cohort,
      lineage_id = seg_res$segment$lineage_id,
      lineage_group = seg_res$segment$lineage_group,
      lineage_label = seg_res$segment$lineage_label,
      scenario_id = seg_res$segment$scenario_id,
      lineage_terminal_key = seg_res$segment$lineage_terminal_key,
      passage_index = seg_res$segment$passage_index,
      lineage_passage_index = seg_res$segment$lineage_passage_index,
      oxygen_pct = seg_res$segment$oxygen_pct,
      passage_id = seg_res$segment$passage_id,
      passage_duration = seg_res$segment$passage_duration,
      endpoint_day = seg_res$selection$endpoint_day,
      selected_day = seg_res$selection$selected_day,
      closest_day_diagnostic = seg_res$selection$closest_day_diagnostic,
      quantile_prob = probs,
      predicted_quantile_kary_N = quantiles,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(out)
}

ivt_collect_daily_counts <- function(run) {
  do.call(rbind, lapply(run$segment_results, function(seg_res) {
    seg_res <- .ivt_normalize_segment_result(seg_res, grid_pre = run$grid_pre)
    seg <- seg_res$segment
    days <- seg$obs_days_local
    n_days <- length(days)
    sim_vec <- function(name, default = NA_real_) {
      vals <- seg_res$sim[[name]]
      if (is.null(vals)) return(rep_len(default, n_days))
      vals <- as.numeric(vals)
      if (length(vals) == n_days) return(vals)
      rep_len(vals, n_days)
    }
    live_cells <- sim_vec("Ntot_live_obs")
    dead_hypoxia_cells <- sim_vec("Ntot_dead_hypoxia_obs", 0)
    dead_buffer_cells <- sim_vec("Ntot_dead_buffer_obs", 0)
    dead_total_cells <- sim_vec("Ntot_dead_total_obs", dead_hypoxia_cells + dead_buffer_cells)
    total_cells <- sim_vec("Ntot_total_obs", live_cells + dead_total_cells)
    burden_live <- sim_vec("Vmm3_live_obs", live_cells)
    burden_dead_hypoxia <- sim_vec("Vmm3_dead_hypoxia_obs", 0)
    burden_dead_buffer <- sim_vec("Vmm3_dead_buffer_obs", 0)
    burden_dead_total <- sim_vec("Vmm3_dead_total_obs", burden_dead_hypoxia + burden_dead_buffer)
    burden_total <- sim_vec("Vmm3_total_obs", burden_live + burden_dead_total)
    o2_target <- sim_vec("O2_target_obs")
    o2_eff <- sim_vec("O2_eff_obs")
    cumulative_gross_divisions <- sim_vec(
      "cumulative_gross_divisions_obs"
    )
    cumulative_hypoxia_deaths <- sim_vec(
      "cumulative_hypoxia_deaths_obs"
    )
    cumulative_dead_buffer_inflow <- sim_vec(
      "cumulative_dead_buffer_inflow_obs"
    )
    cumulative_nonlive_inflow <- sim_vec(
      "cumulative_nonlive_inflow_obs",
      cumulative_hypoxia_deaths + cumulative_dead_buffer_inflow
    )
    data.frame(
      segment_id = seg$segment_id,
      cohort = seg$cohort,
      lineage_id = seg$lineage_id,
      lineage_group = seg$lineage_group,
      lineage_label = seg$lineage_label,
      scenario_id = seg$scenario_id,
      lineage_terminal_key = seg$lineage_terminal_key,
      passage_index = seg$passage_index,
      lineage_passage_index = seg$lineage_passage_index,
      oxygen_pct = seg$oxygen_pct,
      passage_id = seg$passage_id,
      passage_duration = seg$passage_duration,
      observed_passage_day = seg_res$selection$observed_passage_day,
      last_observation_day = seg_res$selection$last_observation_day,
      search_horizon_day = seg_res$selection$search_horizon_day,
      search_horizon_live_cells =
        seg_res$selection$search_horizon_live_cells,
      max_live_cells_in_search = seg_res$selection$max_live_cells_in_search,
      endpoint_day = seg_res$selection$endpoint_day,
      day = days,
      live_cells = live_cells,
      dead_hypoxia_cells = dead_hypoxia_cells,
      dead_buffer_cells = dead_buffer_cells,
      dead_total_cells = dead_total_cells,
      total_cells = total_cells,
      burden_live = burden_live,
      burden_dead_hypoxia = burden_dead_hypoxia,
      burden_dead_buffer = burden_dead_buffer,
      burden_dead_total = burden_dead_total,
      burden_total = burden_total,
      assigned_o2 = seg$oxygen_pct,
      o2_target = o2_target,
      o2_eff = o2_eff,
      cumulative_gross_divisions = cumulative_gross_divisions,
      cumulative_hypoxia_deaths = cumulative_hypoxia_deaths,
      cumulative_dead_buffer_inflow =
        cumulative_dead_buffer_inflow,
      cumulative_nonlive_inflow = cumulative_nonlive_inflow,
      selected_day = seg_res$selection$selected_day,
      selected_live_cells = seg_res$selection$selected_live_cells,
      selected_day_after_last_observation =
        seg_res$selection$selected_day_after_last_observation,
      passage_time_tolerance_days =
        seg_res$selection$passage_time_tolerance_days,
      passage_time_residual_days =
        seg_res$selection$passage_time_residual_days,
      passage_time_within_tolerance =
        seg_res$selection$passage_time_within_tolerance,
      closest_day_diagnostic = seg_res$selection$closest_day_diagnostic,
      target_live_cells = seg_res$selection$target_live_cells,
      protocol_threshold_cells =
        seg_res$selection$protocol_threshold_cells,
      protocol_threshold_source =
        seg_res$selection$protocol_threshold_source,
      observed_final_cells = seg_res$selection$observed_final_cells,
      threshold_target_cells =
        seg_res$selection$threshold_target_cells,
      effective_threshold_cells =
        seg_res$selection$effective_threshold_cells,
      threshold_target_source =
        seg_res$selection$threshold_target_source,
      threshold_reached_by_endpoint =
        seg_res$selection$threshold_reached_by_endpoint,
      predicted_threshold_crossing_day =
        seg_res$selection$predicted_threshold_crossing_day,
      threshold_time_residual_days =
        seg_res$selection$threshold_time_residual_days,
      endpoint_cell_count_residual =
        seg_res$selection$endpoint_cell_count_residual,
      cell_count_overshoot = seg_res$selection$cell_count_overshoot,
      protocol_threshold_overshoot =
        seg_res$selection$protocol_threshold_overshoot,
      protocol_threshold_reached_at_selected =
        seg_res$selection$protocol_threshold_reached_at_selected,
      protocol_threshold_reached_at_observation =
        seg_res$selection$protocol_threshold_reached_at_observation,
      observed_final_cell_count_residual =
        seg_res$selection$observed_final_cell_count_residual,
      threshold_time_grid_resolution_days =
        seg_res$selection$threshold_time_grid_resolution_days,
      threshold_crossing_interval_width_days =
        seg_res$selection$threshold_crossing_interval_width_days,
      passage_recorded = seg_res$selection$passage_recorded,
      passage_executed = seg_res$selection$passage_executed,
      passage_failure_reason = seg_res$selection$passage_failure_reason,
      reseed_mode = seg_res$selection$reseed_mode,
      available_cells = seg_res$selection$available_cells,
      required_cells = seg_res$selection$required_cells,
      supply_ratio = seg_res$selection$supply_ratio,
      boundary_scale = seg_res$selection$boundary_scale,
      stringsAsFactors = FALSE
    )
  }))
}

.ivt_select_existing_columns <- function(df, columns) {
  if (!is.data.frame(df)) return(data.frame())
  df[, intersect(columns, names(df)), drop = FALSE]
}

.ivt_collect_division_death_diagnostics <- function(summary_df) {
  .ivt_select_existing_columns(summary_df, c(
    "cohort", "lineage_id", "scenario_id", "passage_id",
    "passage_index", "lineage_passage_index", "oxygen_pct",
    "passage_duration", "observed_passage_day", "last_observation_day",
    "search_horizon_day",
    "search_horizon_live_cells", "max_live_cells_in_search",
    "endpoint_day", "selected_day", "selected_live_cells",
    "selected_day_after_last_observation",
    "passage_time_tolerance_days",
    "passage_time_residual_days", "passage_time_within_tolerance",
    "predicted_initial_cells", "predicted_final_cells",
    "predicted_passage_live_cells",
    "observed_initial_cells", "observed_final_cells",
    "predicted_live_cells_at_observation",
    "observed_live_cells_at_observation",
    "measurement_day_cell_count_residual", "log_live_cell_residual",
    "observed_net_population_doublings",
    "predicted_net_population_doublings",
    "observed_minimum_division_events",
    "predicted_gross_division_events",
    "predicted_cumulative_hypoxia_deaths",
    "predicted_cumulative_dead_buffer_inflow",
    "predicted_cumulative_nonlive_inflow",
    "predicted_divisions_per_initial_cell",
    "predicted_nonlive_inflow_to_division_ratio",
    "predicted_gross_division_to_net_gain_ratio",
    "cumulative_experimental_time",
    "cumulative_observed_net_population_doublings",
    "cumulative_predicted_net_population_doublings",
    "cumulative_gross_divisions",
    "cumulative_hypoxia_deaths",
    "cumulative_dead_buffer_inflow",
    "cumulative_nonlive_inflow",
    "dead_buffer_inflow_definition"
  ))
}

.ivt_collect_protocol_feasibility <- function(summary_df) {
  if (!is.data.frame(summary_df)) return(data.frame())
  out <- summary_df
  if (!"insufficient_boundary" %in% names(out)) {
    out$insufficient_boundary <-
      out$reseed_mode == "carry_forward_insufficient"
  }
  n_insufficient <- sum(out$insufficient_boundary, na.rm = TRUE)
  out$n_insufficient_boundaries <- n_insufficient
  out$all_passage_boundaries_feasible <- n_insufficient == 0L
  out$protocol_feasibility_status <- if (
    n_insufficient == 0L
  ) {
    "PASS"
  } else {
    "FAIL"
  }
  out$protocol_boundary_status <- ifelse(
    out$reseed_mode == "terminal_no_reseed",
    "TERMINAL_NO_RESEED",
    ifelse(out$insufficient_boundary, "FAIL", "PASS")
  )
  .ivt_select_existing_columns(out, c(
    "cohort", "lineage_id", "scenario_id", "passage_id",
    "passage_index", "lineage_passage_index", "passage_duration",
    "observed_passage_day", "last_observation_day",
    "search_horizon_day", "endpoint_day",
    "search_horizon_live_cells", "max_live_cells_in_search",
    "selected_day", "selected_live_cells", "passage_executed",
    "selected_day_after_last_observation",
    "passage_failure_reason",
    "protocol_threshold_cells", "protocol_threshold_source",
    "protocol_culture_vessel", "protocol_growth_area_cm2",
    "protocol_target_confluence", "observed_final_cells",
    "effective_threshold_cells", "threshold_target_source",
    "protocol_threshold_overshoot",
    "protocol_threshold_reached_at_selected",
    "protocol_threshold_reached_at_observation",
    "observed_final_cell_count_residual",
    "reseed_mode", "protocol_boundary_status",
    "insufficient_boundary", "boundary_feasible",
    "available_cells", "required_cells", "supply_ratio",
    "boundary_scale", "cell_number_before", "cell_number_after",
    "n_insufficient_boundaries",
    "all_passage_boundaries_feasible",
    "protocol_feasibility_status"
  ))
}

.ivt_collect_threshold_crossing_diagnostics <- function(summary_df) {
  .ivt_select_existing_columns(summary_df, c(
    "cohort", "lineage_id", "scenario_id", "passage_id",
    "passage_index", "lineage_passage_index",
    "passage_duration", "observed_passage_day", "last_observation_day",
    "search_horizon_day",
    "search_horizon_live_cells", "max_live_cells_in_search",
    "endpoint_day", "selected_day", "selected_live_cells",
    "selected_day_after_last_observation",
    "passage_time_tolerance_days",
    "passage_time_residual_days", "passage_time_within_tolerance",
    "protocol_threshold_cells", "protocol_threshold_source",
    "protocol_culture_vessel", "protocol_growth_area_cm2",
    "protocol_target_confluence", "observed_final_cells",
    "threshold_target_cells", "threshold_target_source",
    "effective_threshold_cells",
    "threshold_reached_by_endpoint",
    "predicted_threshold_crossing_day", "observed_passage_day",
    "threshold_time_residual_days", "endpoint_cell_count_residual",
    "cell_count_overshoot", "protocol_threshold_overshoot",
    "protocol_threshold_reached_at_selected",
    "protocol_threshold_reached_at_observation",
    "observed_final_cell_count_residual",
    "threshold_time_grid_resolution_days",
    "threshold_crossing_interval_width_days",
    "closest_day_diagnostic", "closest_live_cells_diagnostic"
  ))
}

.ivt_attach_endpoint_instantaneous_growth <- function(summary_df,
                                                      components) {
  if (!is.data.frame(summary_df)) return(summary_df)
  rate_columns <- c(
    "predicted_endpoint_instantaneous_net_growth_rate",
    "predicted_endpoint_division_linked_live_rate",
    "predicted_endpoint_hypoxia_death_rate",
    "predicted_endpoint_crowding_multiplier"
  )
  for (column in rate_columns) {
    if (!column %in% names(summary_df)) {
      summary_df[[column]] <- rep(NA_real_, nrow(summary_df))
    }
  }
  if (!nrow(summary_df) ||
      !"segment_id" %in% names(summary_df) || !is.list(components)) {
    return(summary_df)
  }

  rate_rows <- lapply(c("run_2N", "run_4N"), function(run_name) {
    run <- components[[run_name]]
    if (!is.list(run) || !length(run$segment_results)) return(NULL)
    cfg <- .first_non_null_local(components$cfg, run$simulation_cfg)
    run_params <- .first_non_null_local(
      components$run_params,
      run$shared_run_params
    )
    model_core <- run$model_core
    if (!is.list(cfg) || !is.list(run_params) || !is.list(model_core)) {
      return(NULL)
    }
    rate_cache <- new.env(parent = emptyenv())

    rows <- lapply(run$segment_results, function(seg_res) {
      seg_res <- .ivt_normalize_segment_result(
        seg_res,
        grid_pre = run$grid_pre
      )
      selected_index <- suppressWarnings(as.integer(
        seg_res$selection$selected_index
      ))
      simulated_o2 <- suppressWarnings(as.numeric(seg_res$sim$O2_eff_obs))
      live_state <- suppressWarnings(as.numeric(
        seg_res$selection$selected_state
      ))
      if ((!length(live_state) || any(!is.finite(live_state))) &&
          is.matrix(seg_res$sim$live_state_obs) &&
          length(selected_index) == 1L && is.finite(selected_index) &&
          selected_index >= 1L &&
          selected_index <= nrow(seg_res$sim$live_state_obs)) {
        live_state <- as.numeric(
          seg_res$sim$live_state_obs[selected_index, , drop = TRUE]
        )
      }
      if (length(selected_index) != 1L || !is.finite(selected_index) ||
          selected_index < 1L || selected_index > length(simulated_o2) ||
          !is.finite(simulated_o2[[selected_index]]) ||
          length(live_state) != length(run$grid_pre) ||
          any(!is.finite(live_state))) {
        return(NULL)
      }
      rates <- .ivt_instantaneous_net_growth_at_state(
        live_state = live_state,
        oxygen_pct = simulated_o2[[selected_index]],
        cfg = cfg,
        run_params = run_params,
        model_core = model_core,
        rate_cache = rate_cache
      )
      data.frame(
        segment_id = as.character(seg_res$segment$segment_id),
        predicted_endpoint_instantaneous_net_growth_rate =
          rates$net_growth_rate,
        predicted_endpoint_division_linked_live_rate =
          rates$division_linked_live_rate,
        predicted_endpoint_hypoxia_death_rate = rates$hypoxia_death_rate,
        predicted_endpoint_crowding_multiplier = rates$crowding_multiplier,
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(rows)
  })
  rate_df <- dplyr::bind_rows(rate_rows)
  if (!nrow(rate_df)) return(summary_df)
  if (anyDuplicated(rate_df$segment_id)) {
    stop("Endpoint instantaneous growth diagnostics require unique segment IDs.")
  }
  index <- match(as.character(summary_df$segment_id), rate_df$segment_id)
  matched <- !is.na(index)
  for (column in rate_columns) {
    summary_df[[column]][matched] <- rate_df[[column]][index[matched]]
  }
  summary_df
}

ivt_collect_postfit_tables <- function(components) {
  if (!is.list(components) ||
      !is.list(components$run_2N) ||
      !is.list(components$run_4N)) {
    stop("In-vitro components must contain run_2N and run_4N.")
  }
  summary_df <- as.data.frame(components$summary, stringsAsFactors = FALSE)
  summary_df <- .ivt_attach_endpoint_instantaneous_growth(
    summary_df,
    components
  )
  insufficient_boundary <- if (
    "insufficient_boundary" %in% names(summary_df)
  ) {
    as.logical(summary_df$insufficient_boundary)
  } else {
    summary_df$reseed_mode == "carry_forward_insufficient"
  }
  n_insufficient <- sum(insufficient_boundary, na.rm = TRUE)
  all_boundaries_feasible <- n_insufficient == 0L
  feasibility_status <- if (all_boundaries_feasible) "PASS" else "FAIL"
  summary_df$n_insufficient_boundaries <- n_insufficient
  summary_df$all_passage_boundaries_feasible <- all_boundaries_feasible
  summary_df$protocol_feasibility_status <- feasibility_status
  growth_df <- components$growth_df
  if (is.data.frame(growth_df)) {
    growth_df$n_insufficient_boundaries <- n_insufficient
    growth_df$all_passage_boundaries_feasible <- all_boundaries_feasible
    growth_df$protocol_feasibility_status <- feasibility_status
  }
  growth_count_diagnostics_df <- components$growth_count_diagnostics_df
  if (!is.data.frame(growth_count_diagnostics_df)) {
    growth_count_diagnostics_df <- data.frame()
  }
  if (nrow(growth_count_diagnostics_df) > 0L) {
    growth_count_diagnostics_df$n_insufficient_boundaries <- n_insufficient
    growth_count_diagnostics_df$all_passage_boundaries_feasible <-
      all_boundaries_feasible
    growth_count_diagnostics_df$protocol_feasibility_status <-
      feasibility_status
  }
  passage_audit <- .ivt_select_existing_columns(summary_df, c(
    "cohort", "lineage_id", "scenario_id", "passage_id",
    "passage_index", "lineage_passage_index",
    "passage_duration", "observed_passage_day", "last_observation_day",
    "search_horizon_day",
    "search_horizon_live_cells", "max_live_cells_in_search",
    "endpoint_day", "selected_day", "selected_live_cells",
    "selected_day_after_last_observation",
    "passage_time_tolerance_days",
    "passage_time_residual_days", "passage_time_within_tolerance",
    "closest_day_diagnostic",
    "predicted_initial_cells", "predicted_final_cells",
    "predicted_passage_live_cells",
    "observed_initial_cells", "observed_final_cells",
    "predicted_live_cells_at_observation",
    "observed_live_cells_at_observation",
    "measurement_day_cell_count_residual", "log_live_cell_residual",
    "predicted_growth", "predicted_growth_rate", "observed_growth",
    "predicted_endpoint_instantaneous_net_growth_rate",
    "predicted_endpoint_division_linked_live_rate",
    "predicted_endpoint_hypoxia_death_rate",
    "predicted_endpoint_crowding_multiplier",
    "observed_endpoint_local_net_growth_rate",
    "observed_endpoint_local_interval_start_day",
    "observed_endpoint_local_interval_end_day",
    "observed_endpoint_local_n_positive_timepoints",
    "observed_net_population_doublings",
    "predicted_net_population_doublings",
    "observed_minimum_division_events",
    "predicted_gross_division_events",
    "predicted_cumulative_hypoxia_deaths",
    "predicted_cumulative_dead_buffer_inflow",
    "predicted_cumulative_nonlive_inflow",
    "passage_recorded", "passage_executed", "passage_failure_reason",
    "reseed_mode", "insufficient_boundary",
    "boundary_feasible", "available_cells", "required_cells",
    "supply_ratio", "boundary_scale", "cell_number_before",
    "cell_number_after", "cumulative_time",
    "cumulative_experimental_time",
    "n_insufficient_boundaries",
    "all_passage_boundaries_feasible",
    "protocol_feasibility_status",
    "protocol_threshold_cells", "protocol_threshold_source",
    "protocol_culture_vessel", "protocol_growth_area_cm2",
    "protocol_target_confluence",
    "threshold_target_cells", "threshold_target_source",
    "effective_threshold_cells",
    "threshold_reached_by_endpoint",
    "predicted_threshold_crossing_day", "observed_passage_day",
    "threshold_time_residual_days", "endpoint_cell_count_residual",
    "cell_count_overshoot", "protocol_threshold_overshoot",
    "protocol_threshold_reached_at_selected",
    "protocol_threshold_reached_at_observation",
    "observed_final_cell_count_residual",
    "threshold_time_grid_resolution_days",
    "threshold_crossing_interval_width_days"
  ))
  list(
    invitro_lineage_summary = summary_df,
    invitro_passage_audit = passage_audit,
    invitro_growth_loglik = growth_df,
    invitro_growth_count_diagnostics = growth_count_diagnostics_df,
    invitro_passage_time_loglik = if (is.data.frame(components$passage_time_df)) components$passage_time_df else data.frame(),
    invitro_death_loglik = if (is.data.frame(components$death_df)) components$death_df else data.frame(),
    invitro_daily_counts = dplyr::bind_rows(
      ivt_collect_daily_counts(components$run_2N),
      ivt_collect_daily_counts(components$run_4N)
    ),
    invitro_division_death_diagnostics =
      .ivt_collect_division_death_diagnostics(summary_df),
    invitro_protocol_feasibility =
      .ivt_collect_protocol_feasibility(summary_df),
    invitro_threshold_crossing_diagnostics =
      .ivt_collect_threshold_crossing_diagnostics(summary_df)
  )
}
