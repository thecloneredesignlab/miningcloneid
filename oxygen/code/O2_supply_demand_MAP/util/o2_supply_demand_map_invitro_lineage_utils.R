# Canonical in-vitro observation and lineage mapping helpers.

ivt_observed_passage_summary <- function(fit_entry) {
  g_val <- suppressWarnings(as.numeric(fit_entry$g))
  if (length(g_val) == 0L || !is.finite(g_val)) g_val <- NA_real_
  duration_val <- suppressWarnings(as.numeric(fit_entry$passage_duration))
  if (length(duration_val) == 0L || !is.finite(duration_val)) duration_val <- NA_real_
  initial_cells_val <- suppressWarnings(as.numeric(fit_entry$initial_cells))
  if (length(initial_cells_val) == 0L || !is.finite(initial_cells_val)) initial_cells_val <- NA_real_
  final_cells_val <- suppressWarnings(as.numeric(fit_entry$final_cells))
  if (length(final_cells_val) == 0L || !is.finite(final_cells_val)) final_cells_val <- NA_real_
  kary <- suppressWarnings(as.numeric(fit_entry$kary))
  kary <- kary[is.finite(kary)]
  flow_entry <- fit_entry$flow
  observed_flow <- NULL
  observed_flow_n <- 0L
  flow_sample_name <- NA_character_
  if (is.list(flow_entry) && !is.data.frame(flow_entry) && !is.null(flow_entry$g0g1_ploidy_density)) {
    observed_flow <- flow_entry$g0g1_ploidy_density
    observed_flow_n <- nrow(observed_flow)
    if (!is.null(flow_entry$sample_name_canonical)) {
      flow_sample_name <- as.character(flow_entry$sample_name_canonical)
    }
  }
  list(
    observed_growth = g_val,
    passage_duration = duration_val,
    initial_cells = initial_cells_val,
    final_cells = final_cells_val,
    observed_mean_kary_N = if (length(kary) > 0L) mean(kary) else NA_real_,
    observed_n_kary = length(kary),
    observed_kary = kary,
    observed_n_flow = observed_flow_n,
    observed_flow = observed_flow,
    observed_flow_sample_name = flow_sample_name
  )
}

ivt_segment_median_value <- function(observed_list, field, default = NA_real_) {
  vals <- vapply(observed_list, function(x) {
    v <- suppressWarnings(as.numeric(x[[field]]))
    if (length(v) == 0L || !is.finite(v)) NA_real_ else v
  }, numeric(1))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(as.numeric(default))
  stats::median(vals)
}

ivt_nested_observed_median <- function(observed_nested_list, field, default = NA_real_) {
  vals <- unlist(lapply(observed_nested_list, function(obs_list) {
    vapply(obs_list, function(x) {
      v <- suppressWarnings(as.numeric(x[[field]]))
      if (length(v) == 0L || !is.finite(v)) NA_real_ else v
    }, numeric(1))
  }), use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(as.numeric(default))
  stats::median(vals)
}

.ivt_parse_lineage_passage_id <- function(passage_id, cohort = NULL) {
  passage_id <- as.character(passage_id)
  hit <- regexec("_([24]N)_(C|O1|O2)_A([0-9]+)_seed$", passage_id, perl = TRUE)
  parts <- regmatches(passage_id, hit)[[1]]
  if (length(parts) != 4L) return(NULL)
  parsed_cohort <- parts[[2]]
  if (!is.null(cohort) && !identical(parsed_cohort, as.character(cohort))) return(NULL)
  list(
    cohort = parsed_cohort,
    lineage_id = parts[[3]],
    passage_number = as.integer(parts[[4]])
  )
}

.ivt_observation_has_likelihood_data <- function(obs) {
  is.finite(suppressWarnings(as.numeric(obs$observed_growth))) ||
    length(obs$observed_kary) > 0L ||
    as.integer(obs$observed_n_flow) > 0L
}

.ivt_segment_endpoint_data_ids <- function(run, segment) {
  formal_ids <- as.character(segment$data_ids)
  landmark_records <- run$landmark_observations
  if (is.null(landmark_records) || length(landmark_records) == 0L) {
    return(formal_ids)
  }

  target_segment_id <- as.character(segment$segment_id)
  is_target <- vapply(landmark_records, function(record) {
    record_target <- record$target_segment_id
    length(record_target) == 1L &&
      !is.na(record_target) &&
      identical(as.character(record_target), target_segment_id)
  }, logical(1))
  landmark_ids <- vapply(
    landmark_records[is_target],
    `[[`,
    character(1),
    "passage_id"
  )
  endpoint_ids <- c(formal_ids, landmark_ids)
  if (anyDuplicated(endpoint_ids)) {
    stop("Endpoint observation appears more than once for segment ", target_segment_id, ".")
  }
  endpoint_ids
}

.ivt_build_independent_scenario_adapter <- function(jobs,
                                                    fit_data,
                                                    cohort,
                                                    obs_days_local = NULL,
                                                    passage_time_tolerance_days = 1) {
  cohort <- as.character(cohort)
  if (!cohort %in% c("2N", "4N")) {
    stop("cohort must be '2N' or '4N'.")
  }
  passage_time_tolerance_days <- suppressWarnings(as.numeric(passage_time_tolerance_days))
  if (length(passage_time_tolerance_days) != 1L ||
      !is.finite(passage_time_tolerance_days) ||
      passage_time_tolerance_days < 0) {
    stop("passage_time_tolerance_days must be one finite non-negative value.")
  }
  required_job_columns <- c("sim_key", "parent_key", "oxygen", "depth", "data_ids")
  missing_job_columns <- setdiff(required_job_columns, names(jobs))
  if (length(missing_job_columns) > 0L) {
    stop("jobs is missing required columns: ", paste(missing_job_columns, collapse = ", "))
  }

  formal_records <- list()
  initial_observations <- list()
  landmark_observations <- list()
  seen_data_ids <- character()
  formal_n <- 0L
  initial_n <- 0L
  landmark_n <- 0L

  for (i in seq_len(nrow(jobs))) {
    ids <- as.character(jobs$data_ids[[i]])
    for (passage_id in ids) {
      if (passage_id %in% seen_data_ids) {
        stop("data_id appears in more than one in-vitro job: ", passage_id)
      }
      seen_data_ids <- c(seen_data_ids, passage_id)
      if (is.null(fit_data[[passage_id]])) {
        stop("Missing fit_data entry for in-vitro data_id: ", passage_id)
      }
      obs <- ivt_observed_passage_summary(fit_data[[passage_id]])
      obs$passage_id <- passage_id
      identity <- .ivt_parse_lineage_passage_id(passage_id, cohort = cohort)
      is_formal_passage <- !is.null(identity) &&
        is.finite(obs$passage_duration) && obs$passage_duration > 0 &&
        is.finite(obs$initial_cells) && obs$initial_cells > 0 &&
        is.finite(obs$final_cells) && obs$final_cells > 0

      record <- list(
        passage_id = passage_id,
        observed = obs,
        source_sim_key = as.character(jobs$sim_key[[i]]),
        source_parent_key = as.character(jobs$parent_key[[i]]),
        oxygen_pct = as.numeric(jobs$oxygen[[i]]),
        source_depth = as.integer(jobs$depth[[i]])
      )
      if (is_formal_passage) {
        formal_n <- formal_n + 1L
        record$lineage_id <- identity$lineage_id
        record$passage_number <- identity$passage_number
        formal_records[[formal_n]] <- record
      } else if (.ivt_observation_has_likelihood_data(obs)) {
        if (!is.finite(record$source_depth)) {
          stop("Non-formal observation has no finite source depth: ", passage_id)
        }
        if (record$source_depth == 0L) {
          initial_n <- initial_n + 1L
          record$cohort <- cohort
          record$lineage_id <- "INITIAL"
          record$lineage_group <- "initial"
          record$scenario_id <- paste(cohort, "INITIAL", sep = "-")
          record$segment_id <- paste(cohort, "INITIAL", initial_n, sep = "-")
          record$passage_index <- 0L
          record$lineage_passage_index <- 0L
          initial_observations[[initial_n]] <- record
        } else {
          landmark_n <- landmark_n + 1L
          landmark_observations[[landmark_n]] <- record
        }
      }
    }
  }

  if (length(formal_records) == 0L) {
    stop("No formal C/O1/O2 passage records were found for cohort ", cohort, ".")
  }

  lineage_order <- c("C", "O1", "O2")
  segments <- list()
  segment_n <- 0L
  scenario_ids <- character()
  for (lineage_id in lineage_order) {
    records <- Filter(
      function(x) identical(as.character(x$lineage_id), lineage_id),
      formal_records
    )
    if (length(records) == 0L) next
    passage_numbers <- vapply(records, `[[`, integer(1), "passage_number")
    records <- records[order(passage_numbers)]
    passage_numbers <- vapply(records, `[[`, integer(1), "passage_number")
    if (anyDuplicated(passage_numbers)) {
      stop("Duplicate passage number in ", cohort, "-", lineage_id, " scenario.")
    }
    if (!identical(passage_numbers, seq_len(length(records)))) {
      stop(
        "Passage sequence must be contiguous from A1 for ", cohort, "-", lineage_id,
        "; found A", paste(passage_numbers, collapse = ", A"), "."
      )
    }

    scenario_id <- paste(cohort, lineage_id, sep = "-")
    scenario_ids <- c(scenario_ids, scenario_id)
    terminal_segment_id <- paste0(scenario_id, "-A", max(passage_numbers))
    previous_segment_id <- NA_character_
    for (j in seq_along(records)) {
      record <- records[[j]]
      obs <- record$observed
      observed_passage_day <- as.numeric(obs$passage_duration)
      search_horizon_day <- observed_passage_day + passage_time_tolerance_days
      local_days <- if (is.null(obs_days_local)) {
        sort(unique(c(
          seq(0, search_horizon_day, by = 1),
          observed_passage_day,
          search_horizon_day
        )))
      } else {
        custom_days <- suppressWarnings(as.numeric(obs_days_local))
        custom_days <- custom_days[
          is.finite(custom_days) &
            custom_days >= 0 &
            custom_days <= search_horizon_day
        ]
        sort(unique(c(0, custom_days, observed_passage_day, search_horizon_day)))
      }
      segment_id <- paste0(scenario_id, "-A", record$passage_number)
      segment_n <- segment_n + 1L
      segments[[segment_n]] <- list(
        segment_id = segment_id,
        parent_segment_id = previous_segment_id,
        source_sim_key = record$source_sim_key,
        source_parent_key = record$source_parent_key,
        cohort = cohort,
        lineage_id = lineage_id,
        lineage_group = if (identical(lineage_id, "C")) "control" else "deprived",
        lineage_label = lineage_id,
        scenario_id = scenario_id,
        lineage_terminal_key = terminal_segment_id,
        oxygen_pct = record$oxygen_pct,
        passage_id = record$passage_id,
        passage_number = record$passage_number,
        passage_index = j,
        lineage_passage_index = j,
        duration_days = search_horizon_day,
        passage_duration = observed_passage_day,
        observed_passage_day = observed_passage_day,
        last_observation_day = observed_passage_day,
        search_horizon_day = search_horizon_day,
        passage_time_tolerance_days = passage_time_tolerance_days,
        endpoint_day = search_horizon_day,
        initial_cells = as.numeric(obs$initial_cells),
        final_cells = as.numeric(obs$final_cells),
        obs_days_local = local_days,
        depth = record$source_depth,
        data_ids = record$passage_id,
        observed = list(obs)
      )
      previous_segment_id <- segment_id
    }
  }

  if (length(landmark_observations) > 0L) {
    for (i in seq_along(landmark_observations)) {
      record <- landmark_observations[[i]]
      target_idx <- which(vapply(segments, function(segment) {
        identical(
          as.character(segment$source_sim_key),
          as.character(record$source_sim_key)
        )
      }, logical(1)))
      if (length(target_idx) != 1L) {
        stop(
          "Non-root landmark observation ", record$passage_id,
          " must resolve to exactly one formal source segment; found ",
          length(target_idx), "."
        )
      }
      target <- segments[[target_idx[[1]]]]
      record$cohort <- target$cohort
      record$lineage_id <- target$lineage_id
      record$lineage_group <- target$lineage_group
      record$lineage_label <- target$lineage_label
      record$scenario_id <- target$scenario_id
      record$segment_id <- target$segment_id
      record$target_segment_id <- target$segment_id
      record$lineage_terminal_key <- target$lineage_terminal_key
      record$passage_index <- target$passage_index
      record$lineage_passage_index <- target$lineage_passage_index
      record$oxygen_pct <- target$oxygen_pct
      landmark_observations[[i]] <- record
    }
  }

  list(
    cohort = cohort,
    terminal_sim_key = NA_character_,
    n_segments = length(segments),
    n_scenarios = length(scenario_ids),
    scenario_ids = scenario_ids,
    segments = segments,
    initial_observations = initial_observations,
    landmark_observations = landmark_observations
  )
}

ivt_choose_demo_terminal_key <- function(jobs, fit_data, min_kary = 1L) {
  parent_keys <- unique(stats::na.omit(jobs$parent_key))
  terminal <- jobs[!jobs$sim_key %in% parent_keys, , drop = FALSE]
  if (nrow(terminal) == 0L) stop("No terminal lineage keys found.")
  score_kary <- vapply(terminal$data_ids, function(ids) {
    sum(vapply(ids, function(id) {
      k <- fit_data[[id]]$kary
      !is.null(k) && any(is.finite(k))
    }, logical(1)))
  }, integer(1))
  score_g <- vapply(terminal$data_ids, function(ids) {
    sum(vapply(ids, function(id) is.finite(suppressWarnings(as.numeric(fit_data[[id]]$g))), logical(1)))
  }, integer(1))
  ord <- order(-(score_kary >= min_kary), -terminal$depth, -score_kary, -score_g)
  terminal$sim_key[[ord[[1]]]]
}

ivt_trace_lineage <- function(jobs, terminal_sim_key) {
  idx <- match(terminal_sim_key, jobs$sim_key)
  if (is.na(idx)) stop("terminal_sim_key not found in jobs: ", terminal_sim_key)
  path_idx <- integer(0)
  repeat {
    path_idx <- c(idx, path_idx)
    parent_key <- jobs$parent_key[[idx]]
    if (is.na(parent_key) || !nzchar(parent_key)) break
    idx <- match(parent_key, jobs$sim_key)
    if (is.na(idx)) stop("Missing parent key in jobs: ", parent_key)
  }
  jobs[path_idx, , drop = FALSE]
}

ivt_list_terminal_keys <- function(jobs) {
  parent_keys <- unique(stats::na.omit(jobs$parent_key))
  terminal <- jobs[!jobs$sim_key %in% parent_keys, , drop = FALSE]
  terminal$sim_key
}

ivt_build_lineage_adapter <- function(jobs,
                                      fit_data,
                                      terminal_sim_key,
                                      cohort,
                                      max_segment_days = 14,
                                      obs_days_local = NULL,
                                      passage_time_tolerance_days = 1) {
  adapter <- .ivt_build_independent_scenario_adapter(
    jobs = jobs,
    fit_data = fit_data,
    cohort = cohort,
    obs_days_local = obs_days_local,
    passage_time_tolerance_days = passage_time_tolerance_days
  )
  lineage_jobs <- ivt_trace_lineage(jobs, terminal_sim_key)
  path_keys <- as.character(lineage_jobs$sim_key)
  adapter$segments <- Filter(
    function(seg) as.character(seg$source_sim_key) %in% path_keys,
    adapter$segments
  )
  adapter$initial_observations <- Filter(
    function(obs) as.character(obs$source_sim_key) %in% path_keys,
    adapter$initial_observations
  )
  retained_segment_ids <- vapply(
    adapter$segments,
    `[[`,
    character(1),
    "segment_id"
  )
  adapter$landmark_observations <- Filter(
    function(obs) as.character(obs$target_segment_id) %in% retained_segment_ids,
    adapter$landmark_observations
  )
  adapter$terminal_sim_key <- terminal_sim_key
  adapter$n_segments <- length(adapter$segments)
  adapter$scenario_ids <- unique(vapply(
    adapter$segments,
    `[[`,
    character(1),
    "scenario_id"
  ))
  adapter$n_scenarios <- length(adapter$scenario_ids)
  adapter
}

ivt_make_passage_map <- function(adapter) {
  do.call(rbind, lapply(adapter$segments, function(seg) {
    do.call(rbind, lapply(seg$observed, function(obs) {
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
        observed_passage_day = seg$observed_passage_day,
        search_horizon_day = seg$search_horizon_day,
        passage_time_tolerance_days = seg$passage_time_tolerance_days,
        endpoint_day = seg$endpoint_day,
        initial_cells = seg$initial_cells,
        final_cells = seg$final_cells,
        passage_id = obs$passage_id,
        observed_growth = obs$observed_growth,
        observed_passage_duration = obs$passage_duration,
        observed_initial_cells = obs$initial_cells,
        observed_final_cells = obs$final_cells,
        observed_mean_kary_N = obs$observed_mean_kary_N,
        observed_n_kary = obs$observed_n_kary,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

ivt_terminal_path_map <- function(jobs, cohort) {
  records <- list()
  record_n <- 0L
  for (i in seq_len(nrow(jobs))) {
    for (passage_id in as.character(jobs$data_ids[[i]])) {
      identity <- .ivt_parse_lineage_passage_id(passage_id, cohort = cohort)
      if (is.null(identity)) next
      record_n <- record_n + 1L
      records[[record_n]] <- data.frame(
        cohort = as.character(cohort),
        lineage_id = identity$lineage_id,
        passage_number = identity$passage_number,
        passage_id = passage_id,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(records) == 0L) return(data.frame())
  records <- do.call(rbind, records)
  records <- records[order(match(records$lineage_id, c("C", "O1", "O2")), records$passage_number), , drop = FALSE]
  records$scenario_id <- paste(records$cohort, records$lineage_id, sep = "-")
  records$segment_id <- paste0(records$scenario_id, "-A", records$passage_number)
  records$parent_segment_id <- ave(
    records$segment_id,
    records$scenario_id,
    FUN = function(x) c(NA_character_, head(x, -1L))
  )
  terminal_by_scenario <- tapply(records$segment_id, records$scenario_id, function(x) tail(x, 1L))
  records$lineage_terminal_key <- unname(terminal_by_scenario[records$scenario_id])
  records$lineage_label <- records$lineage_id
  records$lineage_group <- ifelse(records$lineage_id == "C", "control", "deprived")
  records$lineage_passage_index <- records$passage_number
  records[, c(
    "cohort", "lineage_id", "lineage_group", "lineage_label", "scenario_id",
    "lineage_terminal_key", "segment_id", "parent_segment_id",
    "lineage_passage_index", "passage_id"
  ), drop = FALSE]
}
