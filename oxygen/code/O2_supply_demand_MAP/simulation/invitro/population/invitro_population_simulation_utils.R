# Population-count outputs derived from fitted in-vitro simulations.

ivt_sim_endpoint_summary_map <- function(run) {
  if (!is.list(run) || !is.list(run$segment_results) || !length(run$segment_results)) {
    return(data.frame())
  }
  dplyr::bind_rows(lapply(run$segment_results, function(seg_res) {
    seg_res <- .ivt_normalize_segment_result(
      seg_res,
      grid_pre = run$grid_pre
    )
    seg <- seg_res$segment
    sim_live <- suppressWarnings(as.numeric(seg_res$sim$Ntot_live_obs))
    endpoint_index <- suppressWarnings(as.integer(seg_res$selection$endpoint_index))
    if (length(endpoint_index) != 1L ||
        !is.finite(endpoint_index) ||
        endpoint_index < 1L ||
        endpoint_index > length(sim_live)) {
      endpoint_index <- length(sim_live)
    }
    predicted_initial_cells <- if (length(sim_live)) sim_live[[1L]] else NA_real_
    predicted_final_cells <- if (length(sim_live)) sim_live[[endpoint_index]] else NA_real_
    sim_terminal_scalar <- function(terminal_name, obs_name) {
      terminal_value <- suppressWarnings(as.numeric(seg_res$sim[[terminal_name]]))
      if (length(terminal_value) == 1L && is.finite(terminal_value)) {
        return(terminal_value)
      }
      obs_value <- suppressWarnings(as.numeric(seg_res$sim[[obs_name]]))
      if (length(obs_value) >= endpoint_index &&
          is.finite(obs_value[[endpoint_index]])) {
        return(obs_value[[endpoint_index]])
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
    predicted_net_gain <- predicted_final_cells - predicted_initial_cells
    predicted_growth <- ivt_log_growth_rate(
      initial_cells = predicted_initial_cells,
      final_cells = predicted_final_cells,
      duration_days = seg$passage_duration
    )
    data.frame(
      cohort = seg$cohort,
      segment_id = seg$segment_id,
      duration_days = seg$duration_days,
      passage_duration = seg$passage_duration,
      endpoint_day = seg_res$selection$endpoint_day,
      selected_day = seg_res$selection$selected_day,
      closest_day_diagnostic = seg_res$selection$closest_day_diagnostic,
      closest_live_cells_diagnostic = seg_res$selection$closest_live_cells_diagnostic,
      threshold_target_cells =
        seg_res$selection$threshold_target_cells,
      threshold_target_source =
        seg_res$selection$threshold_target_source,
      threshold_reached_by_endpoint =
        seg_res$selection$threshold_reached_by_endpoint,
      predicted_threshold_crossing_day =
        seg_res$selection$predicted_threshold_crossing_day,
      observed_passage_day =
        seg_res$selection$observed_passage_day,
      threshold_time_residual_days =
        seg_res$selection$threshold_time_residual_days,
      endpoint_cell_count_residual =
        seg_res$selection$endpoint_cell_count_residual,
      threshold_time_grid_resolution_days =
        seg_res$selection$threshold_time_grid_resolution_days,
      threshold_crossing_interval_width_days =
        seg_res$selection$threshold_crossing_interval_width_days,
      predicted_initial_cells = predicted_initial_cells,
      predicted_final_cells = predicted_final_cells,
      predicted_live_cells = predicted_final_cells,
      predicted_growth = predicted_growth,
      predicted_growth_rate = predicted_growth,
      observed_net_population_doublings = .ivt_log2_population_change(
        seg$initial_cells,
        seg$final_cells
      ),
      predicted_net_population_doublings =
        .ivt_log2_population_change(
          predicted_initial_cells,
          predicted_final_cells
        ),
      observed_minimum_division_events = if (
        is.finite(seg$initial_cells) &&
          is.finite(seg$final_cells)
      ) {
        max(seg$final_cells - seg$initial_cells, 0)
      } else {
        NA_real_
      },
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
      predicted_mean_kary_N = suppressWarnings(as.numeric(.ivt_first_scalar(
        seg_res$selection$predicted_mean_kary_N,
        default = NA_real_
      ))),
      passage_recorded = seg_res$selection$passage_recorded,
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
      stringsAsFactors = FALSE
    )
  }))
}

ivt_sim_refresh_lineage_summary <- function(components) {
  summary_df <- as.data.frame(components$summary, stringsAsFactors = FALSE)
  if (!nrow(summary_df)) return(summary_df)
  summary_df <- ivt_sim_normalize_lineage_columns(summary_df)
  endpoint_map <- dplyr::bind_rows(
    ivt_sim_endpoint_summary_map(components$run_2N),
    ivt_sim_endpoint_summary_map(components$run_4N)
  )
  if (!nrow(endpoint_map) ||
      !all(c("cohort", "segment_id") %in% names(summary_df))) {
    return(summary_df)
  }
  endpoint_key <- paste(endpoint_map$cohort, endpoint_map$segment_id, sep = "\r")
  if (anyDuplicated(endpoint_key)) {
    stop("Endpoint summary map contains duplicate cohort/segment keys.")
  }
  summary_key <- paste(summary_df$cohort, summary_df$segment_id, sep = "\r")
  endpoint_index <- match(summary_key, endpoint_key)
  matched <- !is.na(endpoint_index)
  value_columns <- setdiff(
    names(endpoint_map),
    c("cohort", "segment_id")
  )
  for (column in value_columns) {
    if (!column %in% names(summary_df)) {
      prototype <- endpoint_map[[column]]
      summary_df[[column]] <- if (is.character(prototype)) {
        rep(NA_character_, nrow(summary_df))
      } else if (is.logical(prototype)) {
        rep(NA, nrow(summary_df))
      } else {
        rep(NA_real_, nrow(summary_df))
      }
    }
    summary_df[[column]][matched] <- endpoint_map[[column]][endpoint_index[matched]]
  }
  if ("scenario_id" %in% names(summary_df)) {
    cumulative_columns <- c(
      cumulative_experimental_time = "passage_duration",
      cumulative_observed_net_population_doublings =
        "observed_net_population_doublings",
      cumulative_predicted_net_population_doublings =
        "predicted_net_population_doublings",
      cumulative_gross_divisions =
        "predicted_gross_division_events",
      cumulative_hypoxia_deaths =
        "predicted_cumulative_hypoxia_deaths",
      cumulative_dead_buffer_inflow =
        "predicted_cumulative_dead_buffer_inflow",
      cumulative_nonlive_inflow =
        "predicted_cumulative_nonlive_inflow"
    )
    for (output_column in names(cumulative_columns)) {
      input_column <- cumulative_columns[[output_column]]
      if (input_column %in% names(summary_df)) {
        summary_df[[output_column]] <- ave(
          suppressWarnings(as.numeric(summary_df[[input_column]])),
          summary_df$scenario_id,
          FUN = cumsum
        )
      }
    }
    if ("cumulative_experimental_time" %in% names(summary_df)) {
      summary_df$cumulative_time <-
        summary_df$cumulative_experimental_time
    }
  }
  if ("reseed_mode" %in% names(summary_df)) {
    n_insufficient <- sum(
      summary_df$reseed_mode == "carry_forward_insufficient",
      na.rm = TRUE
    )
    summary_df$n_insufficient_boundaries <- n_insufficient
    summary_df$all_passage_boundaries_feasible <-
      n_insufficient == 0L
    summary_df$protocol_feasibility_status <- if (
      n_insufficient == 0L
    ) {
      "PASS"
    } else {
      "FAIL"
    }
  }
  summary_df
}

ivt_sim_collect_population_tables <- function(components) {
  run_2n <- components$run_2N
  run_4n <- components$run_4N
  if (!is.list(run_2n) || !is.list(run_4n)) {
    stop("In-vitro best_components must contain run_2N and run_4N.")
  }
  summary_df <- ivt_sim_refresh_lineage_summary(components)
  list(
    invitro_lineage_summary = summary_df,
    invitro_daily_counts = dplyr::bind_rows(
      ivt_collect_daily_counts(run_2n),
      ivt_collect_daily_counts(run_4n)
    ),
    invitro_division_death_diagnostics =
      .ivt_collect_division_death_diagnostics(summary_df),
    invitro_protocol_feasibility =
      .ivt_collect_protocol_feasibility(summary_df),
    invitro_threshold_crossing_diagnostics =
      .ivt_collect_threshold_crossing_diagnostics(summary_df)
  )
}
