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

ivt_collect_lineage_summary <- function(run, fit_data) {
  do.call(rbind, lapply(run$segment_results, function(seg_res) {
    seg <- seg_res$segment
    pred_mean_kary_N <- as.numeric(seg_res$selection$predicted_mean_kary_N)
    sim_live <- as.numeric(seg_res$sim$Ntot_live_obs)
    selected_idx <- as.integer(seg_res$selection$selected_index)
    if (!is.finite(selected_idx) || selected_idx < 1L || selected_idx > length(sim_live)) {
      selected_idx <- length(sim_live)
    }
    pred_growth <- ivt_log_growth_rate(
      initial_cells = sim_live[[1]],
      final_cells = sim_live[[selected_idx]],
      duration_days = seg_res$selection$selected_day
    )
    do.call(rbind, lapply(seg$data_ids, function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      data.frame(
        segment_id = seg$segment_id,
        parent_segment_id = if (is.null(seg$parent_segment_id)) NA_character_ else as.character(seg$parent_segment_id),
        cohort = seg$cohort,
        passage_index = seg$passage_index,
        oxygen_pct = seg$oxygen_pct,
        duration_days = seg$duration_days,
        initial_cells = seg$initial_cells,
        final_cells = seg$final_cells,
        selected_day = seg_res$selection$selected_day,
        target_live_cells = seg_res$selection$target_live_cells,
        passage_id = pid,
        predicted_live_cells = seg_res$selection$selected_live_cells,
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
}

ivt_collect_observed_kary_summary <- function(run, fit_data) {
  out <- lapply(run$segment_results, function(seg_res) {
    seg <- seg_res$segment
    dplyr::bind_rows(lapply(seg$data_ids, function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_kary <- obs$observed_kary
      observed_kary <- observed_kary[is.finite(observed_kary)]
      if (!length(observed_kary)) return(NULL)
      data.frame(
        segment_id = seg$segment_id,
        cohort = seg$cohort,
        passage_index = seg$passage_index,
        oxygen_pct = seg$oxygen_pct,
        passage_id = pid,
        cell_index = seq_along(observed_kary),
        observed_kary_N = observed_kary,
        stringsAsFactors = FALSE
      )
    }))
  })

  dplyr::bind_rows(out)
}

ivt_collect_observed_flow_summary <- function(run, fit_data) {
  out <- lapply(run$segment_results, function(seg_res) {
    seg <- seg_res$segment
    dplyr::bind_rows(lapply(seg$data_ids, function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_flow <- obs$observed_flow
      if (is.null(observed_flow) || !nrow(observed_flow)) return(NULL)
      data.frame(
        segment_id = seg$segment_id,
        cohort = seg$cohort,
        passage_index = seg$passage_index,
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

  dplyr::bind_rows(out)
}

ivt_collect_distribution_summary <- function(run) {
  do.call(rbind, lapply(run$segment_results, function(seg_res) {
    frac <- as.numeric(seg_res$selection$selected_frac)
    data.frame(
      segment_id = seg_res$segment$segment_id,
      cohort = seg_res$segment$cohort,
      passage_index = seg_res$segment$passage_index,
      oxygen_pct = seg_res$segment$oxygen_pct,
      selected_day = seg_res$selection$selected_day,
      N = run$grid_pre,
      fraction = frac,
      stringsAsFactors = FALSE
    )
  }))
}

ivt_collect_distribution_quantiles <- function(run, probs) {
  probs <- as.numeric(probs)
  out <- lapply(run$segment_results, function(seg_res) {
    frac <- as.numeric(seg_res$selection$selected_frac)
    quantiles <- ivt_weighted_quantile_kary_N(frac, run$grid_pre, probs)
    data.frame(
      segment_id = seg_res$segment$segment_id,
      cohort = seg_res$segment$cohort,
      passage_index = seg_res$segment$passage_index,
      oxygen_pct = seg_res$segment$oxygen_pct,
      selected_day = seg_res$selection$selected_day,
      quantile_prob = probs,
      predicted_quantile_kary_N = quantiles,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(out)
}

ivt_collect_daily_counts <- function(run) {
  do.call(rbind, lapply(run$segment_results, function(seg_res) {
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
      passage_index = seg$passage_index,
      oxygen_pct = seg$oxygen_pct,
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
      target_live_cells = seg_res$selection$target_live_cells,
      stringsAsFactors = FALSE
    )
  }))
}
