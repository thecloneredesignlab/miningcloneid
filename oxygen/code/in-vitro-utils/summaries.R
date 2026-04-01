ivt_weighted_mean_ploidy <- function(weights, grid_pre, chr_lengths_bp = NULL) {
  w <- as.numeric(weights)
  if (length(w) != length(grid_pre)) stop("weights and grid_pre length mismatch.")
  s <- sum(w)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  p <- weighted_ploidy_from_total_N(as.numeric(grid_pre), chr_lengths_bp = chr_lengths_bp)
  sum((w / s) * p)
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
    pred_mean_ploidy <- as.numeric(seg_res$selection$predicted_mean_ploidy)
    sim_live <- as.numeric(seg_res$sim$Ntot_live_obs)
    pred_growth <- ivt_log_growth_rate(
      initial_cells = sim_live[[1]],
      final_cells = sim_live[[length(sim_live)]],
      duration_days = seg$duration_days
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
        passage_id = pid,
        predicted_live_cells = seg_res$selection$selected_live_cells,
        predicted_growth_rate = pred_growth,
        predicted_mean_ploidy = pred_mean_ploidy,
        observed_growth = obs$observed_growth,
        observed_mean_ploidy = obs$observed_mean_ploidy,
        observed_n_kary = obs$observed_n_kary,
        stringsAsFactors = FALSE
      )
    }))
  }))
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

ivt_collect_daily_counts <- function(run) {
  do.call(rbind, lapply(run$segment_results, function(seg_res) {
    seg <- seg_res$segment
    days <- seg$obs_days_local
    data.frame(
      segment_id = seg$segment_id,
      cohort = seg$cohort,
      passage_index = seg$passage_index,
      oxygen_pct = seg$oxygen_pct,
      day = days,
      live_cells = as.numeric(seg_res$sim$Ntot_live_obs),
      selected_day = seg_res$selection$selected_day,
      stringsAsFactors = FALSE
    )
  }))
}
