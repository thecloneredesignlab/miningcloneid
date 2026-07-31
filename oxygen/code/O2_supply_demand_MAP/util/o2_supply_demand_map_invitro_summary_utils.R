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
  has_fixed_endpoint_fields <- !is.null(selection$endpoint_day) && length(selection$endpoint_day)
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
    if (length(obs_days)) length(obs_days) else NA_integer_,
    default = NA_integer_
  ))
  legacy_endpoint <- if (!has_fixed_endpoint_fields) {
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
  selected_day <- if (has_fixed_endpoint_fields || endpoint_was_derived) {
    as.numeric(.ivt_first_scalar(selection$selected_day, endpoint_day, default = endpoint_day))
  } else {
    as.numeric(.ivt_first_scalar(legacy_selected_day, endpoint_day, default = endpoint_day))
  }
  if (endpoint_was_derived) selected_day <- endpoint_day
  closest_day_diagnostic <- suppressWarnings(as.numeric(.ivt_first_scalar(
    selection$closest_day_diagnostic,
    if (is.finite(legacy_selected_day) &&
        is.finite(endpoint_day) &&
        abs(legacy_selected_day - endpoint_day) > 1e-10) {
      legacy_selected_day
    } else {
      NA_real_
    },
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
  seg$endpoint_day <- endpoint_day
  seg$oxygen_pct <- suppressWarnings(as.numeric(.ivt_first_scalar(seg$oxygen_pct, default = NA_real_)))
  fallback_n_days <- if (is.finite(endpoint_index) && endpoint_index > 0L) endpoint_index else 1L
  seg$obs_days_local <- if (length(obs_days)) obs_days else seq_len(fallback_n_days) - 1L

  selection$endpoint_index <- endpoint_index
  selection$endpoint_day <- endpoint_day
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
    default = NA_real_
  )))
  selection$passage_recorded <- as.logical(.ivt_first_scalar(
    selection$passage_recorded,
    length(data_ids) > 0L,
    default = FALSE
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
    endpoint_idx <- as.integer(seg_res$selection$endpoint_index)
    if (!is.finite(endpoint_idx) || endpoint_idx < 1L || endpoint_idx > length(sim_live)) {
      endpoint_idx <- length(sim_live)
    }
    pred_growth <- ivt_log_growth_rate(
      initial_cells = sim_live[[1]],
      final_cells = sim_live[[endpoint_idx]],
      duration_days = seg$passage_duration
    )
    do.call(rbind, lapply(seg$data_ids, function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
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
        endpoint_day = seg_res$selection$endpoint_day,
        initial_cells = seg$initial_cells,
        final_cells = seg$final_cells,
        selected_day = seg_res$selection$selected_day,
        closest_day_diagnostic = seg_res$selection$closest_day_diagnostic,
        closest_live_cells_diagnostic = seg_res$selection$closest_live_cells_diagnostic,
        target_live_cells = seg_res$selection$target_live_cells,
        passage_id = pid,
        predicted_initial_cells = sim_live[[1]],
        predicted_final_cells = sim_live[[endpoint_idx]],
        observed_initial_cells = obs$initial_cells,
        observed_final_cells = obs$final_cells,
        predicted_live_cells = sim_live[[endpoint_idx]],
        passage_recorded = seg_res$selection$passage_recorded,
        reseed_mode = seg_res$selection$reseed_mode,
        available_cells = seg_res$selection$available_cells,
        required_cells = seg_res$selection$required_cells,
        supply_ratio = seg_res$selection$supply_ratio,
        boundary_scale = seg_res$selection$boundary_scale,
        cell_number_before = seg_res$selection$cell_number_before,
        cell_number_after = seg_res$selection$cell_number_after,
        predicted_growth = pred_growth,
        predicted_growth_rate = pred_growth,
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
      out$passage_duration,
      out$scenario_id,
      FUN = cumsum
    )
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
      selected_day = seg_res$selection$selected_day,
      closest_day_diagnostic = seg_res$selection$closest_day_diagnostic,
      target_live_cells = seg_res$selection$target_live_cells,
      passage_recorded = seg_res$selection$passage_recorded,
      reseed_mode = seg_res$selection$reseed_mode,
      available_cells = seg_res$selection$available_cells,
      required_cells = seg_res$selection$required_cells,
      supply_ratio = seg_res$selection$supply_ratio,
      boundary_scale = seg_res$selection$boundary_scale,
      stringsAsFactors = FALSE
    )
  }))
}
