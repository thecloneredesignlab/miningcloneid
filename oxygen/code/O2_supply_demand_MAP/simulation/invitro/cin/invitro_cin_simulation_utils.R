# CIN and missegregation-derived in-vitro simulation tables.

ivt_sim_compute_missegregation_timecourse <- function(dist_df, run_params) {
  empty <- data.frame(
    cohort = character(),
    lineage_id = character(),
    lineage_label = character(),
    scenario_id = character(),
    passage_id = character(),
    lineage_terminal_key = character(),
    segment_id = character(),
    parent_segment_id = character(),
    passage_index = numeric(),
    lineage_passage_index = numeric(),
    oxygen_pct = numeric(),
    endpoint_day = numeric(),
    selected_day = numeric(),
    mean_p_misseg = numeric(),
    weighted_mean_N = numeric(),
    total_fraction = numeric(),
    stringsAsFactors = FALSE
  )
  required <- c("cohort", "passage_index", "oxygen_pct", "N", "fraction")
  if (is.null(dist_df) || !is.data.frame(dist_df) || !all(required %in% names(dist_df))) {
    return(empty)
  }
  if (!exists(".pmisseg_of_O2", mode = "function", inherits = TRUE)) return(empty)

  df <- ivt_sim_normalize_lineage_columns(dist_df)
  numeric_columns <- c(
    "passage_index", "lineage_passage_index", "oxygen_pct", "selected_day", "N", "fraction"
  )
  for (name in numeric_columns) df[[name]] <- suppressWarnings(as.numeric(df[[name]]))
  df$fraction <- pmax(df$fraction, 0)
  df <- df[
    is.finite(df$oxygen_pct) & is.finite(df$N) &
      is.finite(df$fraction) & df$fraction > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(df)) return(empty)
  df$p_misseg <- pmax(pmin(as.numeric(.pmisseg_of_O2(
    O2 = df$oxygen_pct,
    run_params = run_params,
    N = df$N
  )), 1), 0)
  df <- df[is.finite(df$p_misseg), , drop = FALSE]
  if (!nrow(df)) return(empty)

  out <- df |>
    dplyr::group_by(
      .data$cohort,
      .data$lineage_id,
      .data$lineage_label,
      .data$scenario_id,
      .data$passage_id,
      .data$lineage_terminal_key,
      .data$segment_id,
      .data$parent_segment_id,
      .data$passage_index,
      .data$lineage_passage_index,
      .data$oxygen_pct,
      .data$endpoint_day,
      .data$selected_day
    ) |>
    dplyr::summarise(
      mean_p_misseg = sum(.data$fraction * .data$p_misseg, na.rm = TRUE) /
        pmax(sum(.data$fraction, na.rm = TRUE), 1e-12),
      weighted_mean_N = sum(.data$fraction * .data$N, na.rm = TRUE) /
        pmax(sum(.data$fraction, na.rm = TRUE), 1e-12),
      total_fraction = sum(.data$fraction, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$cohort,
      .data$lineage_label,
      .data$lineage_passage_index,
      .data$passage_index,
      .data$segment_id
    )
  as.data.frame(out, stringsAsFactors = FALSE)
}
