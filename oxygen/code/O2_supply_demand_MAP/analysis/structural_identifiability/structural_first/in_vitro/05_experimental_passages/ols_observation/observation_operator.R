# The observation operator consumes a trajectory supplied by the passage model.
# It has no default values for the number or detection of dead cells.

calendar_keys <- function(rows) {
  paste(rows$cohort, rows$model_segment_id, rows$calendar_day, sep = "\t")
}

build_ols_calendar_map <- function(audited_counts, dt_days) {
  stopifnot(length(dt_days) == 1L, is.finite(dt_days), dt_days > 0,
            nrow(audited_counts) > 0L, !anyDuplicated(audited_counts$record_id))
  groups <- split(audited_counts, audited_counts$passage_id)
  mapping <- lapply(groups, function(rows) {
    rows <- rows[order(rows$ols_calendar_day), ]
    day <- rows$ols_calendar_day
    n <- length(day)
    day_sum <- sum(day)
    denominator <- n * sum(day^2) - day_sum^2
    ticks <- day / dt_days
    stopifnot(n >= 2L, day[[1L]] == 0, !anyDuplicated(day),
              all(day == rows$archived_model_day), all(day == round(day)),
              length(unique(rows$model_segment_id)) == 1L,
              length(unique(rows$cohort)) == 1L, denominator > 0,
              all(abs(ticks - round(ticks)) <= 1e-10))
    data.frame(
      record_id = rows$record_id,
      passage_id = rows$passage_id,
      cohort = rows$cohort,
      oxygen_pct = rows$oxygen_pct,
      model_segment_id = rows$model_segment_id,
      calendar_day = day,
      euler_step = as.integer(round(ticks)),
      within_passage_count_role = rows$within_passage_count_role,
      observed_after_seed144_selected_day = rows$ols_day_vs_seed144_selected_day == "after",
      ols_weight_numerator = n * day - day_sum,
      ols_weight_denominator = denominator
    )
  })
  mapping <- do.call(rbind, mapping)
  rownames(mapping) <- NULL
  stopifnot(!anyDuplicated(calendar_keys(mapping)))
  mapping
}

apply_ols_calendar_map <- function(mapping, count_trajectory) {
  required <- c("cohort", "model_segment_id", "calendar_day", "visible_count")
  if (!all(required %in% names(count_trajectory))) {
    stop("Count trajectory needs cohort, model_segment_id, calendar_day and visible_count.")
  }
  keys <- calendar_keys(count_trajectory)
  if (anyNA(count_trajectory[required[1:3]]) || anyDuplicated(keys)) {
    stop("Count trajectory must identify each cohort/segment/day at most once.")
  }
  index <- match(calendar_keys(mapping), keys)
  if (anyNA(index)) {
    stop("Missing trajectory sample for count record: ",
         mapping$record_id[which(is.na(index))[[1L]]])
  }
  count <- count_trajectory$visible_count[index]
  if (!is.numeric(count) || any(!is.finite(count)) || any(count < 0)) {
    stop("Visible counts at all requested sample times must be finite and nonnegative.")
  }
  weights <- mapping$ols_weight_numerator / mapping$ols_weight_denominator
  slopes <- vapply(
    split(weights * log(pmax(1, count)), mapping$passage_id),
    sum, numeric(1)
  )
  data.frame(passage_id = names(slopes), predicted_ols_per_day = unname(slopes))
}

predict_ols_from_trajectory <- function(mapping, state_trajectory) {
  states <- c("live_cells", "dead_hypoxia_cells", "dead_mitotic_biological_cells")
  fractions <- c("hypoxia_visible_fraction", "mitotic_visible_fraction")
  required <- c("cohort", "model_segment_id", "calendar_day", states, fractions)
  missing <- setdiff(required, names(state_trajectory))
  if (length(missing)) {
    stop("State trajectory needs explicit fields: ", paste(missing, collapse = ", "))
  }
  valid_states <- vapply(state_trajectory[states], function(x) {
    is.numeric(x) && all(is.finite(x) & x >= 0)
  }, logical(1))
  valid_fractions <- vapply(state_trajectory[fractions], function(x) {
    is.numeric(x) && all(is.finite(x) & x >= 0 & x <= 1)
  }, logical(1))
  if (!all(valid_states) || !all(valid_fractions)) {
    stop("State stocks must be finite and nonnegative and visibility fractions must be in [0, 1].")
  }
  state_trajectory$visible_count <- with(
    state_trajectory,
    live_cells + hypoxia_visible_fraction * dead_hypoxia_cells +
      mitotic_visible_fraction * dead_mitotic_biological_cells
  )
  apply_ols_calendar_map(mapping, state_trajectory)
}
