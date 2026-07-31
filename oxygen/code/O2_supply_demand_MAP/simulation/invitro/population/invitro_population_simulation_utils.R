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
      predicted_initial_cells = predicted_initial_cells,
      predicted_final_cells = predicted_final_cells,
      predicted_live_cells = predicted_final_cells,
      predicted_growth = predicted_growth,
      predicted_growth_rate = predicted_growth,
      predicted_mean_kary_N = suppressWarnings(as.numeric(.ivt_first_scalar(
        seg_res$selection$predicted_mean_kary_N,
        default = NA_real_
      ))),
      passage_recorded = seg_res$selection$passage_recorded,
      reseed_mode = seg_res$selection$reseed_mode,
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
  summary_df
}

ivt_sim_collect_population_tables <- function(components) {
  run_2n <- components$run_2N
  run_4n <- components$run_4N
  if (!is.list(run_2n) || !is.list(run_4n)) {
    stop("In-vitro best_components must contain run_2N and run_4N.")
  }
  list(
    invitro_lineage_summary = ivt_sim_refresh_lineage_summary(components),
    invitro_daily_counts = dplyr::bind_rows(
      ivt_collect_daily_counts(run_2n),
      ivt_collect_daily_counts(run_4n)
    )
  )
}
