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
                                      obs_days_local = NULL) {
  lineage_jobs <- ivt_trace_lineage(jobs, terminal_sim_key)
  segments <- lapply(seq_len(nrow(lineage_jobs)), function(i) {
    job <- lineage_jobs[i, , drop = FALSE]
    ids <- job$data_ids[[1]]
    obs <- lapply(ids, function(id) {
      out <- ivt_observed_passage_summary(fit_data[[id]])
      out$passage_id <- id
      out
    })
    duration_use <- ivt_segment_median_value(obs, "passage_duration", default = max_segment_days)
    if (!is.finite(duration_use) || duration_use <= 0) duration_use <- as.numeric(max_segment_days)
    initial_cells_use <- ivt_segment_median_value(obs, "initial_cells", default = NA_real_)
    final_cells_use <- ivt_segment_median_value(obs, "final_cells", default = NA_real_)
    local_days <- if (is.null(obs_days_local)) {
      seq(0, duration_use, by = 1)
    } else {
      sort(unique(as.numeric(obs_days_local)))
    }
    list(
      segment_id = job$sim_key[[1]],
      parent_segment_id = job$parent_key[[1]],
      cohort = cohort,
      oxygen_pct = as.numeric(job$oxygen[[1]]),
      duration_days = duration_use,
      initial_cells = initial_cells_use,
      final_cells = final_cells_use,
      obs_days_local = local_days,
      passage_index = i,
      depth = as.integer(job$depth[[1]]),
      data_ids = ids,
      observed = obs
    )
  })

  list(
    cohort = cohort,
    terminal_sim_key = terminal_sim_key,
    n_segments = length(segments),
    segments = segments
  )
}

ivt_make_passage_map <- function(adapter) {
  do.call(rbind, lapply(adapter$segments, function(seg) {
    do.call(rbind, lapply(seg$observed, function(obs) {
      data.frame(
        segment_id = seg$segment_id,
        parent_segment_id = if (is.null(seg$parent_segment_id)) NA_character_ else as.character(seg$parent_segment_id),
        passage_index = seg$passage_index,
        oxygen_pct = seg$oxygen_pct,
        duration_days = seg$duration_days,
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
  terminal_keys <- ivt_list_terminal_keys(jobs)
  do.call(rbind, lapply(terminal_keys, function(term_key) {
    alias <- if (grepl("^[0-9.]+(_[0-9.]+)*$", term_key) &&
      all(strsplit(term_key, "_", fixed = TRUE)[[1]] == "20.5")) {
      "control"
    } else {
      "deprived"
    }
    path <- ivt_trace_lineage(jobs, term_key)
    data.frame(
      cohort = cohort,
      lineage_terminal_key = term_key,
      lineage_label = alias,
      segment_id = path$sim_key,
      lineage_passage_index = seq_len(nrow(path)),
      stringsAsFactors = FALSE
    )
  }))
}
