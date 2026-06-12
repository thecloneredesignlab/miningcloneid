#!/usr/bin/env Rscript

# Factorial in vivo interaction simulations:
# 1) triggered_drug: pre-treatment uses fitted best parameters except O2_S0;
#    post-treatment switches to requested p_mis_base and p_wgd values when
#    total burden reaches a trigger threshold, then follows 500 days by default.
# 2) untreated_factorial: requested O2_S0, p_mis_base, and p_wgd values are
#    active from simulation day 0, with no treatment trigger.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

.interaction_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
})

source(file.path(.interaction_script_dir, "simulate_invivo_mixed_ploidy_perturbations.R"), local = environment())
rm(.interaction_script_dir)

format_value_label <- function(x, digits = 4) {
  format(signif(as.numeric(x), digits), trim = TRUE, scientific = TRUE)
}

format_compact_value_label <- function(x, digits = 3) {
  vapply(as.numeric(x), function(z) {
    if (!is.finite(z)) return(as.character(z))
    use_scientific <- z != 0 && abs(z) < 1e-3
    format(signif(z, digits), trim = TRUE, scientific = use_scientific)
  }, character(1))
}

factor_labels <- function(x, digits = 4) {
  factor(format_value_label(x, digits = digits), levels = format_value_label(sort(unique(x)), digits = digits))
}

make_interaction_design <- function(run_params,
                                    o2_values,
                                    pmis_base_values,
                                    p_wgd_values,
                                    trigger_burden_values,
                                    initial_burden_cells,
                                    experiment_name = "triggered_drug",
                                    scenario_prefix = "drug") {
  baseline_p_mis_base <- as.numeric(run_params$p_mis_base)
  baseline_p_wgd <- as.numeric(run_params$p_wgd)
  expand.grid(
    o2_S0 = as.numeric(o2_values),
    p_wgd_post = as.numeric(p_wgd_values),
    trigger_burden_cells = as.numeric(trigger_burden_values),
    p_mis_base_post = as.numeric(pmis_base_values),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      experiment = experiment_name,
      scenario_id = safe_id(
        scenario_prefix,
        "o2", num_tag(o2_S0),
        "pwgd", num_tag(p_wgd_post),
        "trigger", num_tag(trigger_burden_cells),
        "pmiss", num_tag(p_mis_base_post)
      ),
      initial_burden_cells = as.numeric(initial_burden_cells),
      p_mis_base_pre = baseline_p_mis_base,
      p_wgd_pre = baseline_p_wgd,
      p_wgd = p_wgd_post,
      varied_parameter = "factorial",
      varied_value = NA_real_,
      status = "planned",
      trigger_day = NA_real_,
      actual_trigger_burden_cells = NA_real_,
      post_treatment_duration_day = NA_real_,
      sim_end_day = NA_real_
    ) %>%
    select(
      experiment, scenario_id, varied_parameter, varied_value,
      initial_burden_cells, trigger_burden_cells, trigger_day,
      actual_trigger_burden_cells, o2_S0, p_mis_base_pre,
      p_mis_base_post, p_wgd_pre, p_wgd_post, p_wgd,
      post_treatment_duration_day, sim_end_day, status
    )
}

make_untreated_design <- function(run_params,
                                  o2_values,
                                  pmis_base_values,
                                  p_wgd_values,
                                  initial_burden_cells,
                                  horizon_day) {
  expand.grid(
    o2_S0 = as.numeric(o2_values),
    p_wgd_post = as.numeric(p_wgd_values),
    p_mis_base_post = as.numeric(pmis_base_values),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      experiment = "untreated_factorial",
      scenario_id = safe_id(
        "untreated",
        "o2", num_tag(o2_S0),
        "pwgd", num_tag(p_wgd_post),
        "pmiss", num_tag(p_mis_base_post)
      ),
      initial_burden_cells = as.numeric(initial_burden_cells),
      trigger_burden_cells = NA_real_,
      trigger_day = NA_real_,
      actual_trigger_burden_cells = NA_real_,
      p_mis_base_pre = as.numeric(run_params$p_mis_base),
      p_wgd_pre = as.numeric(run_params$p_wgd),
      p_wgd = p_wgd_post,
      varied_parameter = "factorial",
      varied_value = NA_real_,
      status = "untreated",
      post_treatment_duration_day = NA_real_,
      sim_end_day = as.numeric(horizon_day)
    ) %>%
    select(
      experiment, scenario_id, varied_parameter, varied_value,
      initial_burden_cells, trigger_burden_cells, trigger_day,
      actual_trigger_burden_cells, o2_S0, p_mis_base_pre,
      p_mis_base_post, p_wgd_pre, p_wgd_post, p_wgd,
      post_treatment_duration_day, sim_end_day, status
    )
}

make_exact_ploidy_mixture_state <- function(cfg,
                                            total_live_cells,
                                            ploidy_values,
                                            fractions) {
  ploidy_values <- as.numeric(ploidy_values)
  fractions <- as.numeric(fractions)
  if (length(ploidy_values) != length(fractions)) {
    stop("ploidy_values and fractions must have the same length.")
  }
  if (!length(ploidy_values) || any(!is.finite(ploidy_values)) || any(!is.finite(fractions))) {
    stop("Initial ploidy mixture values must be finite.")
  }
  if (any(fractions < 0) || sum(fractions) <= 0) {
    stop("Initial ploidy mixture fractions must be non-negative and have positive total mass.")
  }

  grid_pre <- cfg$N_MIN:cfg$N_MAX
  target_N <- as.integer(round(ploidy_values * as.numeric(cfg$N_UNIT)))
  if (any(target_N < cfg$N_MIN | target_N > cfg$N_MAX)) {
    stop("Initial ploidy mixture maps outside cfg N range.")
  }

  fractions <- fractions / sum(fractions)
  x <- rep(0, length(grid_pre))
  names(x) <- as.character(grid_pre)
  for (i in seq_along(target_N)) {
    x[as.character(target_N[[i]])] <- x[as.character(target_N[[i]])] +
      as.numeric(total_live_cells) * fractions[[i]]
  }
  x
}

make_initial_state_for_mode <- function(cfg,
                                        total_live_cells,
                                        initial_state_mode,
                                        init_mean,
                                        init_sd,
                                        init_min,
                                        init_max) {
  initial_state_mode <- tolower(as.character(initial_state_mode))
  if (initial_state_mode %in% c("2n4n", "2n_4n", "two_four_equal", "two_four_50_50")) {
    return(make_exact_ploidy_mixture_state(
      cfg = cfg,
      total_live_cells = total_live_cells,
      ploidy_values = c(2, 4),
      fractions = c(0.5, 0.5)
    ))
  }

  if (initial_state_mode %in% c("continuous", "normal", "truncated_normal")) {
    return(make_continuous_ploidy_state(
      cfg = cfg,
      total_live_cells = total_live_cells,
      mean_ploidy = init_mean,
      sd_ploidy = init_sd,
      min_ploidy = init_min,
      max_ploidy = init_max
    ))
  }

  stop("Unknown initial_state_mode=", initial_state_mode, ". Use continuous or 2n4n.")
}

simulate_cpp_light <- function(init_state,
                               run_params,
                               cfg,
                               horizon_day,
                               report_dt,
                               init_dead_hypoxia_state = NULL,
                               init_dead_buffer_state = NULL,
                               start_day = 0,
                               scenario_id = "scenario",
                               experiment = "factorial_interaction",
                               segment = "single",
                               live_min_cells = 1e-6) {
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  horizon_day <- as.numeric(horizon_day)
  report_dt <- as.numeric(report_dt)
  if (!is.finite(horizon_day) || horizon_day < 0) stop("horizon_day must be finite and >= 0.")
  if (!is.finite(report_dt) || report_dt <= 0) stop("report_dt must be finite and > 0.")

  step_max <- as.integer(round(horizon_day / cfg$DT))
  keep_steps <- unique(as.integer(round(seq(0, horizon_day, by = report_dt) / cfg$DT)))
  keep_steps <- sort(unique(c(0L, keep_steps, step_max)))
  keep_steps <- keep_steps[keep_steps >= 0L & keep_steps <= step_max]
  local_days <- as.numeric(keep_steps) * as.numeric(cfg$DT)
  day <- as.numeric(start_day) + local_days
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg))

  o2_s0_upper <- as.numeric(.first_non_null_local(cfg$o2_S0_upper_bound, 5.0))
  o2_S0_use <- as.numeric(.first_non_null_local(run_params$o2_S0, cfg$o2_S0_init, 0.5))
  o2_S0_use <- min(max(o2_S0_use, 0), o2_s0_upper)

  sim_cpp <- cpp_o2simps_simulate_one(list(
    init_state = as.numeric(init_state),
    init_dead_hypoxia_state = as.numeric(init_dead_hypoxia_state %||% rep(0, length(grid_pre))),
    init_dead_buffer_state = as.numeric(init_dead_buffer_state %||% rep(0, length(grid_pre))),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = as.integer(keep_steps),
    sim_end_step = as.integer(step_max),
    DT = as.numeric(cfg$DT),
    dose = 0.0,
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = as.numeric(horizon_day + 1),
    fit_treatment = FALSE,
    alpha = 0.0,
    gamma = 1.0,
    tx_mult_min = as.numeric(cfg$tx_mult_min),
    crowding_enabled = isTRUE(cfg$Crowding),
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_crit = as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0)),
    o2_feedback = isTRUE(cfg$o2_burden_feedback),
    o2_S0 = as.numeric(o2_S0_use),
    kappa_O = as.numeric(.first_non_null_local(run_params$kappa_O, cfg$kappa_O_init, 1.0)),
    tau_O2 = as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0)),
    o2_Nref = as.numeric(.first_non_null_local(cfg$o2_Nref, cfg$init_total_size, 1e6)),
    o2_min = as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0.0)),
    o2_S0_upper_bound = as.numeric(o2_s0_upper),
    eta_o2 = as.numeric(.first_non_null_local(run_params$eta_o2, cfg$eta_o2_init, 1.0)),
    o2_cache_bin_pct = as.numeric(cfg$o2_cache_bin_pct),
    o2_cache_hysteresis_pct = as.numeric(cfg$o2_cache_hysteresis_pct),
    o2_cache_profile = FALSE,
    O2_growth = isTRUE(cfg$O2_growth),
    lam_max = as.numeric(run_params$lam_max),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0.0)),
    boundary = as.character(.first_non_null_local(run_params$boundary, "drop")),
    eps_tail = 1e-8,
    buffer_smax = as.numeric(.first_non_null_local(run_params$buffer_smax, cfg$buffer_smax_init, 0.8)),
    buffer_beta = as.numeric(.first_non_null_local(run_params$buffer_beta, cfg$buffer_beta_init, 1.0)),
    buffer_n_exp = as.numeric(.first_non_null_local(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1.0)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = 0.0,
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
    mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)),
    gamma_mu = as.numeric(.first_non_null_local(run_params$gamma_mu, cfg$gamma_mu_init, 1.0)),
    n_O = as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0)),
    ploidy_O2_death = as.character(cfg$ploidy_O2_death),
    start_with = as.character(cfg$start_with),
    k_clear = as.numeric(.first_non_null_local(run_params$k_clear, cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(cfg$burden_log_eps),
    return_full_trajectory = TRUE
  ))

  live_state <- as.matrix(sim_cpp$live_state_obs)
  dead_hypoxia_state <- as.matrix(sim_cpp$dead_hypoxia_state_obs)
  dead_buffer_state <- as.matrix(sim_cpp$dead_buffer_state_obs)
  if (!identical(dim(live_state), c(length(keep_steps), length(grid_pre)))) {
    stop("Unexpected live_state_obs shape returned by cpp_o2simps_simulate_one.")
  }

  burden <- data.frame(
    scenario_id = scenario_id,
    experiment = experiment,
    segment = segment,
    local_day = local_days,
    day = day,
    step = as.integer(keep_steps),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, 1e-5)),
    p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0.0)),
    o2_S0 = as.numeric(o2_S0_use),
    pred_burden_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_total_obs, sim_cpp$Ntot_obs)),
    pred_burden_live_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_live_obs, sim_cpp$Ntot_obs)),
    pred_burden_dead_hypoxia_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_hypoxia_obs, rep(0, length(keep_steps)))),
    pred_burden_dead_buffer_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_buffer_obs, rep(0, length(keep_steps)))),
    pred_burden_dead_total_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_total_obs, rep(0, length(keep_steps)))),
    pred_burden_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_total_obs, sim_cpp$Vmm3_obs)),
    pred_burden_live_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_live_obs, sim_cpp$Vmm3_obs)),
    pred_burden_dead_total_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_total_obs, rep(0, length(keep_steps)))),
    pred_o2_target_pct = as.numeric(.first_non_null_local(sim_cpp$O2_target_obs, rep(NA_real_, length(keep_steps)))),
    pred_o2_pct = as.numeric(.first_non_null_local(sim_cpp$O2_eff_obs, rep(NA_real_, length(keep_steps)))),
    stringsAsFactors = FALSE
  )
  burden$pred_log10_burden_cells <- log10(pmax(burden$pred_burden_cells, 1e-300))
  burden$pred_log10_burden_live_cells <- log10(pmax(burden$pred_burden_live_cells, 1e-300))
  burden$pred_log10_burden_dead_total_cells <- log10(pmax(burden$pred_burden_dead_total_cells, 1e-300))
  burden$pred_log10_burden_volume_mm3 <- log10(pmax(burden$pred_burden_volume_mm3, 1e-300))
  burden$pred_o2_lag_gap_pct <- burden$pred_o2_target_pct - burden$pred_o2_pct

  live_pop <- rowSums(live_state, na.rm = TRUE)
  frac_sum <- rep(NA_real_, length(live_pop))
  mean_N <- rep(NA_real_, length(live_pop))
  mean_ploidy <- rep(NA_real_, length(live_pop))
  frac_lt2 <- rep(NA_real_, length(live_pop))
  frac_2to3 <- rep(NA_real_, length(live_pop))
  frac_3to4 <- rep(NA_real_, length(live_pop))
  frac_ge4 <- rep(NA_real_, length(live_pop))
  frac_ge5 <- rep(NA_real_, length(live_pop))
  ploidy_grid <- as.numeric(grid_pre) / as.numeric(cfg$N_UNIT)
  has_live <- is.finite(live_pop) & live_pop > live_min_cells
  for (i in which(has_live)) {
    frac <- as.numeric(live_state[i, ]) / live_pop[[i]]
    s <- sum(frac, na.rm = TRUE)
    frac_sum[[i]] <- s
    if (is.finite(s) && s > 1e-9) {
      mean_N[[i]] <- sum(as.numeric(grid_pre) * frac, na.rm = TRUE) / s
      mean_ploidy[[i]] <- sum(ploidy_grid * frac, na.rm = TRUE) / s
      frac_lt2[[i]] <- sum(frac[ploidy_grid < 2], na.rm = TRUE)
      frac_2to3[[i]] <- sum(frac[ploidy_grid >= 2 & ploidy_grid < 3], na.rm = TRUE)
      frac_3to4[[i]] <- sum(frac[ploidy_grid >= 3 & ploidy_grid < 4], na.rm = TRUE)
      frac_ge4[[i]] <- sum(frac[ploidy_grid >= 4], na.rm = TRUE)
      frac_ge5[[i]] <- sum(frac[ploidy_grid >= 5], na.rm = TRUE)
    }
  }

  ploidy_summary <- data.frame(
    scenario_id = scenario_id,
    experiment = experiment,
    segment = segment,
    local_day = local_days,
    day = day,
    live_cells = live_pop,
    mean_N = mean_N,
    mean_ploidy = mean_ploidy,
    frac_ploidy_lt2 = frac_lt2,
    frac_ploidy_2to3 = frac_2to3,
    frac_ploidy_3to4 = frac_3to4,
    frac_ploidy_ge4 = frac_ge4,
    frac_ploidy_ge5 = frac_ge5,
    stringsAsFactors = FALSE
  )

  list(
    burden = burden,
    ploidy_summary = ploidy_summary,
    live_state = live_state,
    dead_hypoxia_state = dead_hypoxia_state,
    dead_buffer_state = dead_buffer_state
  )
}

annotate_light_timecourses <- function(sim, design_row) {
  sim$burden$scenario_id <- as.character(design_row$scenario_id)
  sim$burden$experiment <- as.character(design_row$experiment)
  sim$ploidy_summary$scenario_id <- as.character(design_row$scenario_id)
  sim$ploidy_summary$experiment <- as.character(design_row$experiment)
  add_cols <- setdiff(names(design_row), names(sim$burden))
  burden <- bind_cols(as.data.frame(design_row[rep(1, nrow(sim$burden)), add_cols, drop = FALSE]), sim$burden)
  add_cols_ploidy <- setdiff(names(design_row), names(sim$ploidy_summary))
  ploidy <- bind_cols(as.data.frame(design_row[rep(1, nrow(sim$ploidy_summary)), add_cols_ploidy, drop = FALSE]), sim$ploidy_summary)
  list(burden = burden, ploidy_summary = ploidy)
}

simulate_intermittent_post <- function(init_state,
                                       init_dead_hypoxia_state,
                                       init_dead_buffer_state,
                                       run_params_on,
                                       run_params_off,
                                       cfg,
                                       start_day,
                                       followup_day,
                                       on_day = 7,
                                       off_day = 21,
                                       report_dt = 1,
                                       scenario_id,
                                       experiment,
                                       live_min_cells = 1e-6) {
  current_live <- as.numeric(init_state)
  current_dead_hypoxia <- as.numeric(init_dead_hypoxia_state)
  current_dead_buffer <- as.numeric(init_dead_buffer_state)
  current_day <- as.numeric(start_day)
  remaining <- as.numeric(followup_day)
  cycle <- 1L
  segment_index <- 1L
  burden_rows <- list()
  ploidy_rows <- list()

  run_segment <- function(rp, duration, segment_name, treatment_on, cycle_index, keep_start) {
    seg <- simulate_cpp_light(
      init_state = current_live,
      init_dead_hypoxia_state = current_dead_hypoxia,
      init_dead_buffer_state = current_dead_buffer,
      run_params = rp,
      cfg = cfg,
      horizon_day = duration,
      report_dt = report_dt,
      start_day = current_day,
      scenario_id = scenario_id,
      experiment = experiment,
      segment = segment_name,
      live_min_cells = live_min_cells
    )
    current_live <<- as.numeric(seg$live_state[nrow(seg$live_state), ])
    current_dead_hypoxia <<- as.numeric(seg$dead_hypoxia_state[nrow(seg$dead_hypoxia_state), ])
    current_dead_buffer <<- as.numeric(seg$dead_buffer_state[nrow(seg$dead_buffer_state), ])
    current_day <<- current_day + as.numeric(duration)

    keep <- rep(TRUE, nrow(seg$burden))
    if (!keep_start) keep <- seg$burden$local_day > 1e-9
    seg$burden <- seg$burden[keep, , drop = FALSE]
    seg$burden$cycle <- cycle_index
    seg$burden$treatment_on <- treatment_on

    keep_p <- rep(TRUE, nrow(seg$ploidy_summary))
    if (!keep_start) keep_p <- seg$ploidy_summary$local_day > 1e-9
    seg$ploidy_summary <- seg$ploidy_summary[keep_p, , drop = FALSE]
    seg$ploidy_summary$cycle <- cycle_index
    seg$ploidy_summary$treatment_on <- treatment_on
    seg
  }

  while (remaining > 1e-9) {
    duration_on <- min(as.numeric(on_day), remaining)
    if (duration_on > 1e-9) {
      seg <- run_segment(
        rp = run_params_on,
        duration = duration_on,
        segment_name = "drug_on",
        treatment_on = TRUE,
        cycle_index = cycle,
        keep_start = segment_index == 1L
      )
      burden_rows[[segment_index]] <- seg$burden
      ploidy_rows[[segment_index]] <- seg$ploidy_summary
      segment_index <- segment_index + 1L
      remaining <- remaining - duration_on
    }

    duration_off <- min(as.numeric(off_day), remaining)
    if (duration_off > 1e-9) {
      seg <- run_segment(
        rp = run_params_off,
        duration = duration_off,
        segment_name = "drug_off",
        treatment_on = FALSE,
        cycle_index = cycle,
        keep_start = FALSE
      )
      burden_rows[[segment_index]] <- seg$burden
      ploidy_rows[[segment_index]] <- seg$ploidy_summary
      segment_index <- segment_index + 1L
      remaining <- remaining - duration_off
    }
    cycle <- cycle + 1L
  }

  list(
    burden = bind_rows(burden_rows),
    ploidy_summary = bind_rows(ploidy_rows),
    live_state = current_live,
    dead_hypoxia_state = current_dead_hypoxia,
    dead_buffer_state = current_dead_buffer
  )
}

make_endpoint_summary <- function(design_out, burden_all, ploidy_summary) {
  endpoint_burden <- burden_all %>%
    group_by(scenario_id, experiment) %>%
    arrange(day, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    rename(
      endpoint_p_mis_base = p_mis_base,
      endpoint_p_wgd = p_wgd,
      endpoint_o2_S0 = o2_S0
    ) %>%
    select(
      scenario_id, experiment, day, segment, local_day, step,
      endpoint_p_mis_base, endpoint_p_wgd, endpoint_o2_S0,
      starts_with("pred_burden"),
      starts_with("pred_log10"),
      starts_with("pred_o2")
    )

  endpoint_ploidy <- ploidy_summary %>%
    group_by(scenario_id, experiment) %>%
    arrange(day, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    select(
      scenario_id, experiment, day,
      live_cells, mean_N, mean_ploidy,
      frac_ploidy_lt2, frac_ploidy_2to3, frac_ploidy_3to4,
      frac_ploidy_ge4, frac_ploidy_ge5
    )

  design_out %>%
    left_join(endpoint_burden, by = c("scenario_id", "experiment")) %>%
    left_join(endpoint_ploidy, by = c("scenario_id", "experiment", "day"))
}

add_plot_labels <- function(df, run_params, condition_pmiss_prefix = "post p_miss=") {
  if ("p_wgd_post" %in% names(df)) {
    df$.plot_p_wgd <- as.numeric(df$p_wgd_post)
  } else {
    df$.plot_p_wgd <- as.numeric(df$p_wgd)
  }
  if ("p_mis_base_post" %in% names(df)) {
    df$.plot_p_mis_base <- as.numeric(df$p_mis_base_post)
  } else {
    df$.plot_p_mis_base <- as.numeric(df$p_mis_base)
  }

  baseline_pwgd <- as.numeric(run_params$p_wgd)
  pwgd_values <- sort(unique(df$.plot_p_wgd))
  pwgd_ratio <- pwgd_values / baseline_pwgd
  pwgd_map <- setNames(
    paste0(
      format_compact_value_label(pwgd_ratio),
      "x fitted p_wgd (",
      format_compact_value_label(pwgd_values),
      ")"
    ),
    format_value_label(pwgd_values)
  )
  out <- df %>%
    mutate(
      o2_label = factor(
        paste0("O2=", format_compact_value_label(o2_S0)),
        levels = paste0("O2=", format_compact_value_label(sort(unique(o2_S0))))
      ),
      trigger_label = factor(
        paste0("trigger=", format(signif(trigger_burden_cells, 3), scientific = TRUE)),
        levels = paste0("trigger=", format(signif(sort(unique(trigger_burden_cells)), 3), scientific = TRUE))
      ),
      pmiss_label = factor_labels(.plot_p_mis_base),
      pwgd_label_raw = format_value_label(.plot_p_wgd),
      pwgd_label = factor(unname(pwgd_map[pwgd_label_raw]), levels = unname(pwgd_map[format_value_label(pwgd_values)]))
    )
  condition_grid <- expand.grid(
    pmiss_label = levels(out$pmiss_label),
    pwgd_label = levels(out$pwgd_label),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  condition_levels <- paste(condition_grid$pwgd_label, paste0(condition_pmiss_prefix, condition_grid$pmiss_label), sep = " / ")
  out %>%
    mutate(
      condition_label = factor(paste(pwgd_label, paste0(condition_pmiss_prefix, pmiss_label), sep = " / "), levels = condition_levels),
      condition_index = as.integer(condition_label)
    )
}

scale_fill_log10_burden <- function(name, limits = NULL) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#FFF7BC", "#FEC44F", "#FB6A4A", "#CC2C7A", "#4A1486"),
    limits = limits,
    oob = scales::squish,
    na.value = "grey85",
    name = name
  )
}

scale_fill_mean_ploidy <- function(name) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#F3E8FF", "#C084FC", "#7C3AED", "#1D4ED8"),
    values = scales::rescale(c(0, 1.5, 3, 4.5, 6)),
    limits = c(0, 6),
    oob = scales::squish,
    na.value = "grey85",
    name = name
  )
}

reset_generated_dir <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

complete_timecourse_for_plot <- function(plot_df, value_col, fill_value = 0, max_day = NULL) {
  if (is.null(max_day)) {
    max_day <- ceiling(max(plot_df$day, na.rm = TRUE))
  } else {
    max_day <- ceiling(as.numeric(max_day))
  }
  if (!is.finite(max_day) || max_day < 0) max_day <- 0
  scenario_ids <- unique(as.character(plot_df$scenario_id))
  meta_cols <- c(
    "scenario_id", "experiment", "trigger_burden_cells", "trigger_day",
    "actual_trigger_burden_cells", "o2_S0", "p_mis_base_pre",
    "p_mis_base_post", "p_wgd_pre", "p_wgd_post", "p_wgd", "status", "o2_label", "trigger_label",
    "pmiss_label", "pwgd_label", "condition_label", "condition_index", "vary_label", "vary_index",
    "unit_o2_label", "unit_o2_index", "unit_vary_label", "unit_vary_index"
  )
  meta_cols <- intersect(meta_cols, names(plot_df))
  meta <- plot_df[!duplicated(as.character(plot_df$scenario_id)), meta_cols, drop = FALSE]
  meta$scenario_id <- as.character(meta$scenario_id)
  values <- plot_df[, c("scenario_id", "day", value_col), drop = FALSE]
  values$scenario_id <- as.character(values$scenario_id)
  values$.observed_value <- TRUE
  grid <- expand.grid(
    scenario_id = scenario_ids,
    day = seq(0, max_day, by = 1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  out <- merge(grid, meta, by = "scenario_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, values, by = c("scenario_id", "day"), all.x = TRUE, sort = FALSE)
  missing_from_grid <- is.na(out$.observed_value)
  missing_value <- is.na(out[[value_col]])
  if (is.numeric(out[[value_col]])) {
    missing_value <- missing_value | !is.finite(out[[value_col]])
  }
  out[[value_col]][missing_from_grid | missing_value] <- fill_value
  out$.observed_value <- NULL
  out
}

make_treatment_line_segments <- function(endpoint_plot,
                                         row_index_col = "condition_index",
                                         show_treatment_lines = TRUE) {
  if (!isTRUE(show_treatment_lines) || !nrow(endpoint_plot)) return(endpoint_plot[FALSE, , drop = FALSE])
  if (!row_index_col %in% names(endpoint_plot)) return(endpoint_plot[FALSE, , drop = FALSE])
  endpoint_plot <- endpoint_plot[is.finite(endpoint_plot$trigger_day), , drop = FALSE]
  if (!nrow(endpoint_plot)) return(endpoint_plot)

  rows <- vector("list", nrow(endpoint_plot))
  for (i in seq_len(nrow(endpoint_plot))) {
    dr <- endpoint_plot[i, , drop = FALSE]
    treatment_days <- as.numeric(dr$trigger_day)
    is_intermittent <- "intermittent_on_day" %in% names(dr) &&
      "intermittent_off_day" %in% names(dr) &&
      "post_treatment_duration_day" %in% names(dr) &&
      is.finite(as.numeric(dr$intermittent_on_day)) &&
      is.finite(as.numeric(dr$intermittent_off_day)) &&
      is.finite(as.numeric(dr$post_treatment_duration_day)) &&
      as.numeric(dr$intermittent_on_day) > 0 &&
      as.numeric(dr$intermittent_off_day) >= 0

    if (isTRUE(is_intermittent)) {
      period <- as.numeric(dr$intermittent_on_day) + as.numeric(dr$intermittent_off_day)
      if (period > 0) {
        n_starts <- max(1L, as.integer(ceiling((as.numeric(dr$post_treatment_duration_day) - 1e-9) / period)))
        treatment_days <- as.numeric(dr$trigger_day) + seq(0, by = period, length.out = n_starts)
      }
    }

    line_rows <- dr[rep(1, length(treatment_days)), , drop = FALSE]
    line_rows$treatment_day <- treatment_days
    line_rows$y0 <- as.numeric(line_rows[[row_index_col]]) - 0.48
    line_rows$y1 <- as.numeric(line_rows[[row_index_col]]) + 0.48
    rows[[i]] <- line_rows
  }

  bind_rows(rows)
}

plot_one_split_timecourse <- function(plot_df,
                                      endpoint_plot,
                                      fixed_col,
                                      fixed_value,
                                      vary_col,
                                      value_col,
                                      out_pdf,
                                      include_trigger = TRUE,
                                      show_treatment_lines = TRUE,
                                      fill_kind = c("burden", "ploidy"),
                                      y_axis_label = NULL) {
  fill_kind <- match.arg(fill_kind)
  fixed_chr <- as.character(fixed_value)
  sub <- plot_df[as.character(plot_df[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(invisible(FALSE))

  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  sub$vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)
  sub$vary_index <- as.integer(sub$vary_label)
  sub_full <- complete_timecourse_for_plot(sub, value_col, fill_value = 0)

  endpoint_sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  endpoint_sub$vary_label <- factor(as.character(endpoint_sub[[vary_col]]), levels = vary_levels)
  endpoint_sub$vary_index <- as.integer(endpoint_sub$vary_label)
  treatment_lines <- make_treatment_line_segments(
    endpoint_sub,
    row_index_col = "vary_index",
    show_treatment_lines = show_treatment_lines
  )

  facet_layer <- if (include_trigger) facet_grid(trigger_label ~ o2_label) else facet_grid(. ~ o2_label)
  y_axis_label <- y_axis_label %||% vary_col
  p <- ggplot(sub_full, aes(day, vary_index, fill = .data[[value_col]])) +
    geom_raster() +
    facet_layer +
    scale_y_continuous(breaks = seq_along(vary_levels), labels = vary_levels, expand = c(0, 0)) +
    labs(x = "Day", y = y_axis_label) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.ticks.y = element_line(linewidth = 0.2),
      strip.text = element_text(size = 7)
    )
  if (identical(fill_kind, "burden")) {
    p <- p + scale_fill_log10_burden("log10 burden")
  } else {
    p <- p + scale_fill_mean_ploidy("mean ploidy")
  }
  if (nrow(treatment_lines) > 0L) {
    p <- p + geom_segment(
      data = treatment_lines,
      aes(x = treatment_day, xend = treatment_day, y = y0, yend = y1),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.25
    )
  }

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  plot_height <- if (include_trigger) 10 else 4.5
  ggsave(out_pdf, p, width = 20, height = plot_height)
  invisible(TRUE)
}

plot_split_timecourse_heatmaps <- function(endpoint_summary,
                                           burden_all,
                                           ploidy_summary,
                                           run_params,
                                           out_dir,
                                           include_trigger = TRUE,
                                           show_treatment_lines = TRUE,
                                           condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  burden_plot <- add_plot_labels(burden_all, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  ploidy_plot <- add_plot_labels(ploidy_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)

  split_root <- file.path(out_dir, "timecourse_split")
  reset_generated_dir(split_root)

  for (pmiss in levels(burden_plot$pmiss_label)) {
    sub_dir <- file.path(split_root, "fixed_p_mis_base", safe_id("p_miss", pmiss))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(data.frame(fixed_parameter = "p_mis_base", fixed_value = pmiss), file.path(sub_dir, "fixed_value.tsv"))
    plot_one_split_timecourse(
      burden_plot, endpoint_plot,
      fixed_col = "pmiss_label", fixed_value = pmiss,
      vary_col = "pwgd_label", value_col = "pred_log10_burden_cells",
      out_pdf = file.path(sub_dir, "timecourse_log10_burden_by_p_wgd.pdf"),
      include_trigger = include_trigger,
      show_treatment_lines = show_treatment_lines,
      fill_kind = "burden",
      y_axis_label = "p_wgd relative to fitted (actual p_wgd)"
    )
    plot_one_split_timecourse(
      ploidy_plot, endpoint_plot,
      fixed_col = "pmiss_label", fixed_value = pmiss,
      vary_col = "pwgd_label", value_col = "mean_ploidy",
      out_pdf = file.path(sub_dir, "timecourse_mean_ploidy_by_p_wgd.pdf"),
      include_trigger = include_trigger,
      show_treatment_lines = show_treatment_lines,
      fill_kind = "ploidy",
      y_axis_label = "p_wgd relative to fitted (actual p_wgd)"
    )
  }

  for (pwgd in levels(burden_plot$pwgd_label)) {
    sub_dir <- file.path(split_root, "fixed_p_wgd", safe_id("p_wgd", pwgd))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(data.frame(fixed_parameter = "p_wgd", fixed_value = pwgd), file.path(sub_dir, "fixed_value.tsv"))
    plot_one_split_timecourse(
      burden_plot, endpoint_plot,
      fixed_col = "pwgd_label", fixed_value = pwgd,
      vary_col = "pmiss_label", value_col = "pred_log10_burden_cells",
      out_pdf = file.path(sub_dir, "timecourse_log10_burden_by_p_mis_base.pdf"),
      include_trigger = include_trigger,
      show_treatment_lines = show_treatment_lines,
      fill_kind = "burden",
      y_axis_label = "p_mis_base"
    )
    plot_one_split_timecourse(
      ploidy_plot, endpoint_plot,
      fixed_col = "pwgd_label", fixed_value = pwgd,
      vary_col = "pmiss_label", value_col = "mean_ploidy",
      out_pdf = file.path(sub_dir, "timecourse_mean_ploidy_by_p_mis_base.pdf"),
      include_trigger = include_trigger,
      show_treatment_lines = show_treatment_lines,
      fill_kind = "ploidy",
      y_axis_label = "p_mis_base"
    )
  }
  invisible(split_root)
}

plot_one_split_endpoint <- function(endpoint_plot,
                                    fixed_col,
                                    fixed_value,
                                    vary_col,
                                    value_col,
                                    out_pdf,
                                    include_trigger = TRUE,
                                    fill_kind = c("burden", "ploidy"),
                                    y_axis_label = NULL) {
  fill_kind <- match.arg(fill_kind)
  fixed_chr <- as.character(fixed_value)
  sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(invisible(FALSE))

  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  sub$vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)

  y_axis_label <- y_axis_label %||% vary_col
  p <- ggplot(sub, aes(o2_label, vary_label, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.15) +
    labs(x = "O2_S0", y = y_axis_label) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 7)
    )
  if (include_trigger) {
    p <- p + facet_grid(trigger_label ~ .)
  }
  if (identical(fill_kind, "burden")) {
    p <- p + scale_fill_log10_burden("Endpoint\nlog10 burden")
  } else {
    p <- p + scale_fill_mean_ploidy("Endpoint\nmean ploidy")
  }

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  plot_height <- if (include_trigger) 8 else 4.2
  ggsave(out_pdf, p, width = 8, height = plot_height)
  invisible(TRUE)
}

plot_split_endpoint_heatmaps <- function(endpoint_summary,
                                         run_params,
                                         out_dir,
                                         include_trigger = TRUE,
                                         condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  split_root <- file.path(out_dir, "endpoint_split")
  reset_generated_dir(split_root)

  for (pmiss in levels(endpoint_plot$pmiss_label)) {
    sub_dir <- file.path(split_root, "fixed_p_mis_base", safe_id("p_miss", pmiss))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(data.frame(fixed_parameter = "p_mis_base", fixed_value = pmiss), file.path(sub_dir, "fixed_value.tsv"))
    plot_one_split_endpoint(
      endpoint_plot,
      fixed_col = "pmiss_label", fixed_value = pmiss,
      vary_col = "pwgd_label", value_col = "pred_log10_burden_cells",
      out_pdf = file.path(sub_dir, "endpoint_log10_burden_by_o2.pdf"),
      include_trigger = include_trigger,
      fill_kind = "burden",
      y_axis_label = "p_wgd relative to fitted (actual p_wgd)"
    )
    plot_one_split_endpoint(
      endpoint_plot,
      fixed_col = "pmiss_label", fixed_value = pmiss,
      vary_col = "pwgd_label", value_col = "mean_ploidy",
      out_pdf = file.path(sub_dir, "endpoint_mean_ploidy_by_o2.pdf"),
      include_trigger = include_trigger,
      fill_kind = "ploidy",
      y_axis_label = "p_wgd relative to fitted (actual p_wgd)"
    )
  }

  for (pwgd in levels(endpoint_plot$pwgd_label)) {
    sub_dir <- file.path(split_root, "fixed_p_wgd", safe_id("p_wgd", pwgd))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(data.frame(fixed_parameter = "p_wgd", fixed_value = pwgd), file.path(sub_dir, "fixed_value.tsv"))
    plot_one_split_endpoint(
      endpoint_plot,
      fixed_col = "pwgd_label", fixed_value = pwgd,
      vary_col = "pmiss_label", value_col = "pred_log10_burden_cells",
      out_pdf = file.path(sub_dir, "endpoint_log10_burden_by_o2.pdf"),
      include_trigger = include_trigger,
      fill_kind = "burden",
      y_axis_label = "p_mis_base"
    )
    plot_one_split_endpoint(
      endpoint_plot,
      fixed_col = "pwgd_label", fixed_value = pwgd,
      vary_col = "pmiss_label", value_col = "mean_ploidy",
      out_pdf = file.path(sub_dir, "endpoint_mean_ploidy_by_o2.pdf"),
      include_trigger = include_trigger,
      fill_kind = "ploidy",
      y_axis_label = "p_mis_base"
    )
  }
  invisible(split_root)
}

plot_endpoint_unit_metric <- function(sub,
                                      vary_col,
                                      value_col,
                                      fill_kind = c("burden", "ploidy"),
                                      y_axis_label = NULL,
                                      show_y_axis = TRUE,
                                      burden_limits = NULL) {
  fill_kind <- match.arg(fill_kind)
  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  sub$vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)
  y_axis_label <- y_axis_label %||% vary_col

  p <- ggplot(sub, aes(o2_label, vary_label, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.12) +
    coord_fixed(ratio = 1) +
    labs(x = "O2_S0", y = if (show_y_axis) y_axis_label else NULL) +
    theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 7, hjust = 0.5)
    )
  if (!show_y_axis) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  if (identical(fill_kind, "burden")) {
    p <- p + scale_fill_log10_burden("Endpoint\nlog10 burden", limits = burden_limits)
  } else {
    p <- p + scale_fill_mean_ploidy("Endpoint\nmean ploidy")
  }
  p
}

make_endpoint_unit_plot <- function(endpoint_plot,
                                    fixed_col,
                                    fixed_value,
                                    vary_col,
                                    fixed_axis_label,
                                    vary_axis_label,
                                    burden_limits = NULL) {
  fixed_chr <- as.character(fixed_value)
  sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  unit_title <- paste(strwrap(paste0(fixed_axis_label, " = ", fixed_chr), width = 34), collapse = "\n")

  p_burden <- plot_endpoint_unit_metric(
    sub,
    vary_col = vary_col,
    value_col = "pred_log10_burden_cells",
    fill_kind = "burden",
    y_axis_label = vary_axis_label,
    show_y_axis = TRUE,
    burden_limits = burden_limits
  ) +
    ggtitle(paste("burden", unit_title, sep = "\n"))
  p_ploidy <- plot_endpoint_unit_metric(
    sub,
    vary_col = vary_col,
    value_col = "mean_ploidy",
    fill_kind = "ploidy",
    y_axis_label = vary_axis_label,
    show_y_axis = FALSE,
    burden_limits = burden_limits
  ) +
    ggtitle("ploidy")

  wrap_plots(p_burden, p_ploidy, ncol = 2)
}

make_endpoint_unit_grid_plot <- function(endpoint_plot,
                                         fixed_col,
                                         vary_col,
                                         fixed_axis_label,
                                         vary_axis_label,
                                         fixed_levels,
                                         burden_limits = NULL,
                                         page_title = NULL) {
  unit_plots <- lapply(fixed_levels, function(fixed_value) {
    make_endpoint_unit_plot(
      endpoint_plot = endpoint_plot,
      fixed_col = fixed_col,
      fixed_value = fixed_value,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      burden_limits = burden_limits
    )
  })
  unit_plots <- Filter(Negate(is.null), unit_plots)
  if (!length(unit_plots)) return(NULL)

  wrap_plots(unit_plots, ncol = 3, guides = "collect") +
    plot_layout(guides = "collect") +
    plot_annotation(title = page_title)
}

save_endpoint_unit_grid_pages <- function(endpoint_plot,
                                          fixed_col,
                                          vary_col,
                                          fixed_axis_label,
                                          vary_axis_label,
                                          fixed_levels,
                                          out_pdf,
                                          include_trigger = TRUE) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  burden_limits <- range(endpoint_plot$pred_log10_burden_cells, na.rm = TRUE)
  if (any(!is.finite(burden_limits))) burden_limits <- NULL

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  grDevices::pdf(out_pdf, width = 22, height = 16, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_group in page_groups) {
    page_df <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_title <- if (identical(page_group, "all")) {
      paste0(fixed_axis_label, " endpoint unit grid")
    } else {
      paste0(fixed_axis_label, " endpoint unit grid / ", page_group)
    }
    page_plot <- make_endpoint_unit_grid_plot(
      endpoint_plot = page_df,
      fixed_col = fixed_col,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      fixed_levels = fixed_levels,
      burden_limits = burden_limits,
      page_title = page_title
    )
    if (!is.null(page_plot)) print(page_plot)
  }
  invisible(out_pdf)
}

save_endpoint_unit_individual_plots <- function(endpoint_plot,
                                                fixed_col,
                                                vary_col,
                                                fixed_axis_label,
                                                vary_axis_label,
                                                fixed_levels,
                                                out_dir,
                                                include_trigger = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  burden_limits <- range(endpoint_plot$pred_log10_burden_cells, na.rm = TRUE)
  if (any(!is.finite(burden_limits))) burden_limits <- NULL

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  manifest_rows <- list()
  row_i <- 1L
  for (page_group in page_groups) {
    page_df <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_dir <- file.path(out_dir, safe_id(page_group))
    dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)

    for (fixed_value in fixed_levels) {
      unit_plot <- make_endpoint_unit_plot(
        endpoint_plot = page_df,
        fixed_col = fixed_col,
        fixed_value = fixed_value,
        vary_col = vary_col,
        fixed_axis_label = fixed_axis_label,
        vary_axis_label = vary_axis_label,
        burden_limits = burden_limits
      )
      if (is.null(unit_plot)) next

      if (!identical(page_group, "all")) {
        unit_plot <- unit_plot + plot_annotation(title = page_group)
      }

      unit_pdf <- file.path(page_dir, paste0(safe_id("unit", fixed_axis_label, fixed_value), ".pdf"))
      ggsave(unit_pdf, unit_plot, width = 7.4, height = if (identical(page_group, "all")) 4.4 else 4.7)
      manifest_rows[[row_i]] <- data.frame(
        trigger = page_group,
        fixed_parameter = fixed_axis_label,
        fixed_value = fixed_value,
        file = unit_pdf,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }

  if (length(manifest_rows)) {
    write_tsv(bind_rows(manifest_rows), file.path(out_dir, "unit_files.tsv"))
  }
  invisible(out_dir)
}

map_gradient_colors <- function(x,
                                colors,
                                limits = NULL,
                                values = NULL,
                                na_color = "grey85") {
  x <- as.numeric(x)
  out <- rep(na_color, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)

  if (is.null(limits) || length(limits) != 2L || any(!is.finite(limits))) {
    limits <- range(x[ok], na.rm = TRUE)
  }
  limits <- as.numeric(limits)
  if (!is.finite(diff(limits)) || diff(limits) <= 0) {
    limits <- limits[[1]] + c(-0.5, 0.5)
  }

  scaled <- (x - limits[[1]]) / diff(limits)
  scaled <- pmin(1, pmax(0, scaled))
  pal <- scales::gradient_n_pal(colors, values = values)
  out[ok] <- pal(scaled[ok])
  out
}

make_triangle_polygon_rows <- function(tile_df,
                                       x_col,
                                       y_col,
                                       fill_col,
                                       metric_label,
                                       upper = TRUE) {
  pieces <- vector("list", nrow(tile_df))
  for (i in seq_len(nrow(tile_df))) {
    x <- as.numeric(tile_df[[x_col]][[i]])
    y <- as.numeric(tile_df[[y_col]][[i]])
    if (isTRUE(upper)) {
      px <- c(x - 0.5, x + 0.5, x - 0.5)
      py <- c(y + 0.5, y + 0.5, y - 0.5)
    } else {
      px <- c(x + 0.5, x + 0.5, x - 0.5)
      py <- c(y + 0.5, y - 0.5, y - 0.5)
    }
    pieces[[i]] <- data.frame(
      triangle_id = paste(metric_label, i, sep = "_"),
      x = px,
      y = py,
      fill_color = tile_df[[fill_col]][[i]],
      stringsAsFactors = FALSE
    )
  }
  bind_rows(pieces)
}

make_endpoint_triangle_data <- function(sub,
                                        vary_col,
                                        burden_limits = NULL) {
  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  o2_levels <- levels(droplevels(sub$o2_label))
  if (is.null(o2_levels) || !length(o2_levels)) o2_levels <- sort(unique(as.character(sub$o2_label)))

  sub$vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)
  sub$o2_label <- factor(as.character(sub$o2_label), levels = o2_levels)
  sub$.tile_x <- as.integer(sub$o2_label)
  sub$.tile_y <- as.integer(sub$vary_label)
  sub$.burden_fill <- map_gradient_colors(
    sub$pred_log10_burden_cells,
    colors = c("#FFFFFF", "#FFF7BC", "#FEC44F", "#FB6A4A", "#CC2C7A", "#4A1486"),
    limits = burden_limits
  )
  sub$.ploidy_fill <- map_gradient_colors(
    sub$mean_ploidy,
    colors = c("#FFFFFF", "#F3E8FF", "#C084FC", "#7C3AED", "#1D4ED8"),
    limits = c(0, 6),
    values = scales::rescale(c(0, 1.5, 3, 4.5, 6))
  )

  upper <- make_triangle_polygon_rows(sub, ".tile_x", ".tile_y", ".burden_fill", "burden", upper = TRUE)
  lower <- make_triangle_polygon_rows(sub, ".tile_x", ".tile_y", ".ploidy_fill", "ploidy", upper = FALSE)
  list(
    polygons = bind_rows(upper, lower),
    o2_levels = o2_levels,
    vary_levels = vary_levels
  )
}

plot_gradient_legend_identity <- function(colors,
                                          limits,
                                          title,
                                          values = NULL,
                                          breaks = NULL) {
  if (is.null(limits) || length(limits) != 2L || any(!is.finite(limits)) || diff(limits) <= 0) {
    limits <- c(0, 1)
  }
  vals <- seq(limits[[1]], limits[[2]], length.out = 160)
  df <- data.frame(
    x = 1,
    y = vals,
    fill_color = map_gradient_colors(vals, colors = colors, limits = limits, values = values),
    stringsAsFactors = FALSE
  )
  breaks <- breaks %||% pretty(limits, n = 4)
  ggplot(df, aes(x, y, fill = fill_color)) +
    geom_raster() +
    scale_fill_identity() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = breaks, position = "right", expand = c(0, 0)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 7) +
    theme(
      axis.text.y = element_text(size = 6),
      axis.ticks.y = element_line(linewidth = 0.2),
      plot.title = element_text(size = 7, hjust = 0)
    )
}

make_triangle_legend_plot <- function(burden_limits = NULL) {
  burden_legend <- plot_gradient_legend_identity(
    colors = c("#FFFFFF", "#FFF7BC", "#FEC44F", "#FB6A4A", "#CC2C7A", "#4A1486"),
    limits = burden_limits,
    title = "upper\nlog10 burden"
  )
  ploidy_legend <- plot_gradient_legend_identity(
    colors = c("#FFFFFF", "#F3E8FF", "#C084FC", "#7C3AED", "#1D4ED8"),
    limits = c(0, 6),
    values = scales::rescale(c(0, 1.5, 3, 4.5, 6)),
    breaks = 0:6,
    title = "lower\nmean ploidy"
  )
  wrap_plots(burden_legend, ploidy_legend, ncol = 1)
}

make_endpoint_triangle_unit_plot <- function(endpoint_plot,
                                             fixed_col,
                                             fixed_value,
                                             vary_col,
                                             fixed_axis_label,
                                             vary_axis_label,
                                             burden_limits = NULL) {
  fixed_chr <- as.character(fixed_value)
  sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  tri <- make_endpoint_triangle_data(sub, vary_col = vary_col, burden_limits = burden_limits)
  unit_title <- paste(strwrap(paste0(fixed_axis_label, " = ", fixed_chr), width = 34), collapse = "\n")

  ggplot(tri$polygons, aes(x, y, group = triangle_id, fill = fill_color)) +
    geom_polygon(color = "white", linewidth = 0.08) +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = seq_along(tri$o2_levels),
      labels = tri$o2_levels,
      limits = c(0.5, length(tri$o2_levels) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(tri$vary_levels),
      labels = tri$vary_levels,
      limits = c(0.5, length(tri$vary_levels) + 0.5),
      expand = c(0, 0)
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(x = "O2_S0", y = vary_axis_label, title = unit_title) +
    theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      plot.title = element_text(size = 7, hjust = 0.5)
    )
}

make_endpoint_triangle_unit_grid_plot <- function(endpoint_plot,
                                                  fixed_col,
                                                  vary_col,
                                                  fixed_axis_label,
                                                  vary_axis_label,
                                                  fixed_levels,
                                                  burden_limits = NULL,
                                                  page_title = NULL) {
  unit_plots <- lapply(fixed_levels, function(fixed_value) {
    make_endpoint_triangle_unit_plot(
      endpoint_plot = endpoint_plot,
      fixed_col = fixed_col,
      fixed_value = fixed_value,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      burden_limits = burden_limits
    )
  })
  unit_plots <- Filter(Negate(is.null), unit_plots)
  if (!length(unit_plots)) return(NULL)

  grid_plot <- wrap_plots(unit_plots, ncol = 3) +
    plot_annotation(title = page_title)
  (grid_plot | make_triangle_legend_plot(burden_limits)) +
    plot_layout(widths = c(1, 0.08))
}

save_endpoint_triangle_unit_grid_pages <- function(endpoint_plot,
                                                   fixed_col,
                                                   vary_col,
                                                   fixed_axis_label,
                                                   vary_axis_label,
                                                   fixed_levels,
                                                   out_pdf,
                                                   include_trigger = TRUE) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  burden_limits <- range(endpoint_plot$pred_log10_burden_cells, na.rm = TRUE)
  if (any(!is.finite(burden_limits))) burden_limits <- NULL

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  grDevices::pdf(out_pdf, width = 16, height = 22, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_group in page_groups) {
    page_df <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_title <- if (identical(page_group, "all")) {
      paste0(fixed_axis_label, " endpoint triangle unit grid")
    } else {
      paste0(fixed_axis_label, " endpoint triangle unit grid / ", page_group)
    }
    page_plot <- make_endpoint_triangle_unit_grid_plot(
      endpoint_plot = page_df,
      fixed_col = fixed_col,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      fixed_levels = fixed_levels,
      burden_limits = burden_limits,
      page_title = page_title
    )
    if (!is.null(page_plot)) print(page_plot)
  }
  invisible(out_pdf)
}

save_endpoint_triangle_unit_individual_plots <- function(endpoint_plot,
                                                         fixed_col,
                                                         vary_col,
                                                         fixed_axis_label,
                                                         vary_axis_label,
                                                         fixed_levels,
                                                         out_dir,
                                                         include_trigger = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  burden_limits <- range(endpoint_plot$pred_log10_burden_cells, na.rm = TRUE)
  if (any(!is.finite(burden_limits))) burden_limits <- NULL

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  manifest_rows <- list()
  row_i <- 1L
  for (page_group in page_groups) {
    page_df <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_dir <- file.path(out_dir, safe_id(page_group))
    dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)

    for (fixed_value in fixed_levels) {
      unit_plot <- make_endpoint_triangle_unit_plot(
        endpoint_plot = page_df,
        fixed_col = fixed_col,
        fixed_value = fixed_value,
        vary_col = vary_col,
        fixed_axis_label = fixed_axis_label,
        vary_axis_label = vary_axis_label,
        burden_limits = burden_limits
      )
      if (is.null(unit_plot)) next

      unit_plot <- (unit_plot | make_triangle_legend_plot(burden_limits)) +
        plot_layout(widths = c(1, 0.14))
      if (!identical(page_group, "all")) {
        unit_plot <- unit_plot + plot_annotation(title = page_group)
      }

      unit_pdf <- file.path(page_dir, paste0(safe_id("unit_triangle", fixed_axis_label, fixed_value), ".pdf"))
      ggsave(unit_pdf, unit_plot, width = 6.2, height = if (identical(page_group, "all")) 6.2 else 6.5)
      manifest_rows[[row_i]] <- data.frame(
        trigger = page_group,
        fixed_parameter = fixed_axis_label,
        fixed_value = fixed_value,
        file = unit_pdf,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }

  if (length(manifest_rows)) {
    write_tsv(bind_rows(manifest_rows), file.path(out_dir, "unit_files.tsv"))
  }
  invisible(out_dir)
}

make_timecourse_unit_treatment_lines <- function(endpoint_sub,
                                                 max_day,
                                                 cell_width = 2,
                                                 cell_height = 1,
                                                 show_treatment_lines = TRUE) {
  if (!isTRUE(show_treatment_lines) || !nrow(endpoint_sub)) return(endpoint_sub[FALSE, , drop = FALSE])
  required_cols <- c("trigger_day", "unit_o2_index", "unit_vary_index")
  if (!all(required_cols %in% names(endpoint_sub))) return(endpoint_sub[FALSE, , drop = FALSE])

  max_day <- as.numeric(max_day)
  if (!is.finite(max_day) || max_day < 0) max_day <- 0
  denom_day <- max(1, max_day)
  endpoint_sub <- endpoint_sub[is.finite(as.numeric(endpoint_sub$trigger_day)), , drop = FALSE]
  if (!nrow(endpoint_sub)) return(endpoint_sub)

  rows <- list()
  row_i <- 1L
  for (i in seq_len(nrow(endpoint_sub))) {
    dr <- endpoint_sub[i, , drop = FALSE]
    treatment_days <- as.numeric(dr$trigger_day)
    is_intermittent <- "intermittent_on_day" %in% names(dr) &&
      "intermittent_off_day" %in% names(dr) &&
      "post_treatment_duration_day" %in% names(dr) &&
      is.finite(as.numeric(dr$intermittent_on_day)) &&
      is.finite(as.numeric(dr$intermittent_off_day)) &&
      is.finite(as.numeric(dr$post_treatment_duration_day)) &&
      as.numeric(dr$intermittent_on_day) > 0 &&
      as.numeric(dr$intermittent_off_day) >= 0

    if (isTRUE(is_intermittent)) {
      period <- as.numeric(dr$intermittent_on_day) + as.numeric(dr$intermittent_off_day)
      if (period > 0) {
        n_starts <- max(1L, as.integer(ceiling((as.numeric(dr$post_treatment_duration_day) - 1e-9) / period)))
        treatment_days <- as.numeric(dr$trigger_day) + seq(0, by = period, length.out = n_starts)
      }
    }

    treatment_days <- treatment_days[is.finite(treatment_days) & treatment_days >= 0 & treatment_days <= max_day]
    if (!length(treatment_days)) next

    line_rows <- dr[rep(1, length(treatment_days)), , drop = FALSE]
    line_rows$treatment_day <- treatment_days
    line_rows$x0 <- (as.numeric(line_rows$unit_o2_index) - 1) * cell_width +
      (pmin(max_day, pmax(0, treatment_days)) / denom_day) * cell_width
    line_rows$x1 <- line_rows$x0
    line_rows$y0 <- (as.numeric(line_rows$unit_vary_index) - 1) * cell_height
    line_rows$y1 <- as.numeric(line_rows$unit_vary_index) * cell_height
    rows[[row_i]] <- line_rows
    row_i <- row_i + 1L
  }

  if (!length(rows)) return(endpoint_sub[FALSE, , drop = FALSE])
  bind_rows(rows)
}

plot_timecourse_unit_metric <- function(plot_df,
                                        endpoint_plot,
                                        fixed_col,
                                        fixed_value,
                                        vary_col,
                                        value_col,
                                        fill_kind = c("burden", "ploidy"),
                                        y_axis_label = NULL,
                                        show_y_axis = TRUE,
                                        max_day,
                                        show_treatment_lines = TRUE,
                                        burden_limits = NULL,
                                        cell_width = 2,
                                        cell_height = 1) {
  fill_kind <- match.arg(fill_kind)
  fixed_chr <- as.character(fixed_value)
  sub <- plot_df[as.character(plot_df[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(NULL)

  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  o2_levels <- levels(droplevels(sub$o2_label))
  if (is.null(o2_levels) || !length(o2_levels)) o2_levels <- sort(unique(as.character(sub$o2_label)))

  sub$unit_vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)
  sub$unit_vary_index <- as.integer(sub$unit_vary_label)
  sub$unit_o2_label <- factor(as.character(sub$o2_label), levels = o2_levels)
  sub$unit_o2_index <- as.integer(sub$unit_o2_label)

  max_day_i <- ceiling(as.numeric(max_day))
  if (!is.finite(max_day_i) || max_day_i < 0) max_day_i <- 0
  day_count <- max_day_i + 1L
  sub_full <- complete_timecourse_for_plot(sub, value_col, fill_value = 0, max_day = max_day_i)
  sub_full$.x <- (as.numeric(sub_full$unit_o2_index) - 1) * cell_width +
    (as.numeric(sub_full$day) + 0.5) * (cell_width / day_count)
  sub_full$.y <- (as.numeric(sub_full$unit_vary_index) - 0.5) * cell_height

  endpoint_sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (nrow(endpoint_sub)) {
    endpoint_sub$unit_vary_label <- factor(as.character(endpoint_sub[[vary_col]]), levels = vary_levels)
    endpoint_sub$unit_vary_index <- as.integer(endpoint_sub$unit_vary_label)
    endpoint_sub$unit_o2_label <- factor(as.character(endpoint_sub$o2_label), levels = o2_levels)
    endpoint_sub$unit_o2_index <- as.integer(endpoint_sub$unit_o2_label)
  }
  treatment_lines <- make_timecourse_unit_treatment_lines(
    endpoint_sub,
    max_day = max_day_i,
    cell_width = cell_width,
    cell_height = cell_height,
    show_treatment_lines = show_treatment_lines
  )

  border <- expand.grid(
    unit_o2_index = seq_along(o2_levels),
    unit_vary_index = seq_along(vary_levels),
    KEEP.OUT.ATTRS = FALSE
  )
  border$xmin <- (border$unit_o2_index - 1) * cell_width
  border$xmax <- border$unit_o2_index * cell_width
  border$ymin <- (border$unit_vary_index - 1) * cell_height
  border$ymax <- border$unit_vary_index * cell_height

  y_axis_label <- y_axis_label %||% vary_col
  p <- ggplot(sub_full, aes(.x, .y, fill = .data[[value_col]])) +
    geom_tile(width = cell_width / day_count, height = cell_height) +
    geom_rect(
      data = border,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      color = "grey80",
      linewidth = 0.08
    ) +
    scale_x_continuous(
      breaks = (seq_along(o2_levels) - 0.5) * cell_width,
      labels = o2_levels,
      limits = c(0, length(o2_levels) * cell_width),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = (seq_along(vary_levels) - 0.5) * cell_height,
      labels = vary_levels,
      limits = c(0, length(vary_levels) * cell_height),
      expand = c(0, 0)
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(x = "O2_S0", y = if (show_y_axis) y_axis_label else NULL) +
    theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 7, hjust = 0.5)
    )
  if (!show_y_axis) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  if (identical(fill_kind, "burden")) {
    p <- p + scale_fill_log10_burden("Timecourse\nlog10 burden", limits = burden_limits)
  } else {
    p <- p + scale_fill_mean_ploidy("Timecourse\nmean ploidy")
  }
  if (nrow(treatment_lines) > 0L) {
    p <- p + geom_segment(
      data = treatment_lines,
      aes(x = x0, xend = x1, y = y0, yend = y1),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.18
    )
  }
  p
}

make_timecourse_unit_plot <- function(burden_plot,
                                      ploidy_plot,
                                      endpoint_plot,
                                      fixed_col,
                                      fixed_value,
                                      vary_col,
                                      fixed_axis_label,
                                      vary_axis_label,
                                      max_day,
                                      show_treatment_lines = TRUE,
                                      burden_limits = NULL) {
  fixed_chr <- as.character(fixed_value)
  if (!any(as.character(burden_plot[[fixed_col]]) == fixed_chr)) return(NULL)
  unit_title <- paste(strwrap(paste0(fixed_axis_label, " = ", fixed_chr), width = 34), collapse = "\n")

  p_burden <- plot_timecourse_unit_metric(
    plot_df = burden_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = fixed_col,
    fixed_value = fixed_value,
    vary_col = vary_col,
    value_col = "pred_log10_burden_cells",
    fill_kind = "burden",
    y_axis_label = vary_axis_label,
    show_y_axis = TRUE,
    max_day = max_day,
    show_treatment_lines = show_treatment_lines,
    burden_limits = burden_limits
  ) +
    ggtitle(paste("burden", unit_title, sep = "\n"))
  p_ploidy <- plot_timecourse_unit_metric(
    plot_df = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = fixed_col,
    fixed_value = fixed_value,
    vary_col = vary_col,
    value_col = "mean_ploidy",
    fill_kind = "ploidy",
    y_axis_label = vary_axis_label,
    show_y_axis = FALSE,
    max_day = max_day,
    show_treatment_lines = show_treatment_lines,
    burden_limits = burden_limits
  ) +
    ggtitle("ploidy")

  wrap_plots(p_burden, p_ploidy, ncol = 2)
}

make_timecourse_unit_grid_plot <- function(burden_plot,
                                           ploidy_plot,
                                           endpoint_plot,
                                           fixed_col,
                                           vary_col,
                                           fixed_axis_label,
                                           vary_axis_label,
                                           fixed_levels,
                                           max_day,
                                           show_treatment_lines = TRUE,
                                           burden_limits = NULL,
                                           page_title = NULL) {
  unit_plots <- lapply(fixed_levels, function(fixed_value) {
    make_timecourse_unit_plot(
      burden_plot = burden_plot,
      ploidy_plot = ploidy_plot,
      endpoint_plot = endpoint_plot,
      fixed_col = fixed_col,
      fixed_value = fixed_value,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      max_day = max_day,
      show_treatment_lines = show_treatment_lines,
      burden_limits = burden_limits
    )
  })
  unit_plots <- Filter(Negate(is.null), unit_plots)
  if (!length(unit_plots)) return(NULL)

  wrap_plots(unit_plots, ncol = 3, guides = "collect") +
    plot_layout(guides = "collect") +
    plot_annotation(title = page_title)
}

save_timecourse_unit_grid_pages <- function(burden_plot,
                                            ploidy_plot,
                                            endpoint_plot,
                                            fixed_col,
                                            vary_col,
                                            fixed_axis_label,
                                            vary_axis_label,
                                            fixed_levels,
                                            out_pdf,
                                            max_day,
                                            include_trigger = TRUE,
                                            show_treatment_lines = TRUE) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  burden_max <- max(burden_plot$pred_log10_burden_cells, na.rm = TRUE)
  burden_limits <- if (is.finite(burden_max) && burden_max > 0) c(0, burden_max) else c(0, 1)

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  grDevices::pdf(out_pdf, width = 26, height = 18, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_group in page_groups) {
    page_endpoint <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_burden <- if (identical(page_group, "all")) {
      burden_plot
    } else {
      burden_plot[as.character(burden_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_ploidy <- if (identical(page_group, "all")) {
      ploidy_plot
    } else {
      ploidy_plot[as.character(ploidy_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_title <- if (identical(page_group, "all")) {
      paste0(fixed_axis_label, " timecourse unit grid / day 0-", max_day)
    } else {
      paste0(fixed_axis_label, " timecourse unit grid / ", page_group, " / day 0-", max_day)
    }
    page_plot <- make_timecourse_unit_grid_plot(
      burden_plot = page_burden,
      ploidy_plot = page_ploidy,
      endpoint_plot = page_endpoint,
      fixed_col = fixed_col,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      fixed_levels = fixed_levels,
      max_day = max_day,
      show_treatment_lines = show_treatment_lines,
      burden_limits = burden_limits,
      page_title = page_title
    )
    if (!is.null(page_plot)) print(page_plot)
  }
  invisible(out_pdf)
}

save_timecourse_unit_individual_plots <- function(burden_plot,
                                                  ploidy_plot,
                                                  endpoint_plot,
                                                  fixed_col,
                                                  vary_col,
                                                  fixed_axis_label,
                                                  vary_axis_label,
                                                  fixed_levels,
                                                  out_dir,
                                                  max_day,
                                                  include_trigger = TRUE,
                                                  show_treatment_lines = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  burden_max <- max(burden_plot$pred_log10_burden_cells, na.rm = TRUE)
  burden_limits <- if (is.finite(burden_max) && burden_max > 0) c(0, burden_max) else c(0, 1)

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  manifest_rows <- list()
  row_i <- 1L
  for (page_group in page_groups) {
    page_endpoint <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_burden <- if (identical(page_group, "all")) {
      burden_plot
    } else {
      burden_plot[as.character(burden_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_ploidy <- if (identical(page_group, "all")) {
      ploidy_plot
    } else {
      ploidy_plot[as.character(ploidy_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_dir <- file.path(out_dir, safe_id(page_group))
    dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)

    for (fixed_value in fixed_levels) {
      unit_plot <- make_timecourse_unit_plot(
        burden_plot = page_burden,
        ploidy_plot = page_ploidy,
        endpoint_plot = page_endpoint,
        fixed_col = fixed_col,
        fixed_value = fixed_value,
        vary_col = vary_col,
        fixed_axis_label = fixed_axis_label,
        vary_axis_label = vary_axis_label,
        max_day = max_day,
        show_treatment_lines = show_treatment_lines,
        burden_limits = burden_limits
      )
      if (is.null(unit_plot)) next

      if (!identical(page_group, "all")) {
        unit_plot <- unit_plot + plot_annotation(title = page_group)
      }

      unit_pdf <- file.path(page_dir, paste0(safe_id("unit_timecourse", fixed_axis_label, fixed_value), ".pdf"))
      ggsave(unit_pdf, unit_plot, width = 12, height = if (identical(page_group, "all")) 5.2 else 5.5)
      manifest_rows[[row_i]] <- data.frame(
        trigger = page_group,
        fixed_parameter = fixed_axis_label,
        fixed_value = fixed_value,
        max_day = max_day,
        file = unit_pdf,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }

  if (length(manifest_rows)) {
    write_tsv(bind_rows(manifest_rows), file.path(out_dir, "unit_files.tsv"))
  }
  invisible(out_dir)
}

plot_timecourse_unit_grid_heatmaps <- function(endpoint_summary,
                                               burden_all,
                                               ploidy_summary,
                                               run_params,
                                               out_dir,
                                               include_trigger = TRUE,
                                               show_treatment_lines = TRUE,
                                               condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  burden_plot <- add_plot_labels(burden_all, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  ploidy_plot <- add_plot_labels(ploidy_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  max_day <- ceiling(max(c(burden_plot$day, ploidy_plot$day), na.rm = TRUE))
  if (!is.finite(max_day) || max_day < 0) max_day <- 0

  unit_root <- file.path(out_dir, "timecourse_unit_grid")
  reset_generated_dir(unit_root)
  write_tsv(data.frame(time_start_day = 0, time_end_day = max_day), file.path(unit_root, "time_window.tsv"))

  fixed_pmiss_dir <- file.path(unit_root, "fixed_p_mis_base")
  dir.create(fixed_pmiss_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    data.frame(fixed_parameter = "p_mis_base", fixed_value = levels(endpoint_plot$pmiss_label)),
    file.path(fixed_pmiss_dir, "fixed_values.tsv")
  )
  save_timecourse_unit_grid_pages(
    burden_plot = burden_plot,
    ploidy_plot = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_pdf = file.path(fixed_pmiss_dir, "timecourse_burden_ploidy_by_o2_grid.pdf"),
    max_day = max_day,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines
  )
  save_timecourse_unit_individual_plots(
    burden_plot = burden_plot,
    ploidy_plot = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_dir = file.path(fixed_pmiss_dir, "units"),
    max_day = max_day,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines
  )

  fixed_pwgd_dir <- file.path(unit_root, "fixed_p_wgd")
  dir.create(fixed_pwgd_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    data.frame(fixed_parameter = "p_wgd", fixed_value = levels(endpoint_plot$pwgd_label)),
    file.path(fixed_pwgd_dir, "fixed_values.tsv")
  )
  save_timecourse_unit_grid_pages(
    burden_plot = burden_plot,
    ploidy_plot = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_pdf = file.path(fixed_pwgd_dir, "timecourse_burden_ploidy_by_o2_grid.pdf"),
    max_day = max_day,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines
  )
  save_timecourse_unit_individual_plots(
    burden_plot = burden_plot,
    ploidy_plot = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_dir = file.path(fixed_pwgd_dir, "units"),
    max_day = max_day,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines
  )
  invisible(unit_root)
}

plot_endpoint_unit_grid_heatmaps <- function(endpoint_summary,
                                             run_params,
                                             out_dir,
                                             include_trigger = TRUE,
                                             condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  unit_root <- file.path(out_dir, "endpoint_unit_grid")
  reset_generated_dir(unit_root)

  fixed_pmiss_dir <- file.path(unit_root, "fixed_p_mis_base")
  dir.create(fixed_pmiss_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    data.frame(fixed_parameter = "p_mis_base", fixed_value = levels(endpoint_plot$pmiss_label)),
    file.path(fixed_pmiss_dir, "fixed_values.tsv")
  )
  save_endpoint_unit_grid_pages(
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_pdf = file.path(fixed_pmiss_dir, "endpoint_burden_ploidy_by_o2_grid.pdf"),
    include_trigger = include_trigger
  )
  save_endpoint_triangle_unit_grid_pages(
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_pdf = file.path(fixed_pmiss_dir, "endpoint_burden_ploidy_triangle_by_o2_grid.pdf"),
    include_trigger = include_trigger
  )
  save_endpoint_unit_individual_plots(
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_dir = file.path(fixed_pmiss_dir, "units"),
    include_trigger = include_trigger
  )
  save_endpoint_triangle_unit_individual_plots(
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_dir = file.path(fixed_pmiss_dir, "triangle_units"),
    include_trigger = include_trigger
  )

  fixed_pwgd_dir <- file.path(unit_root, "fixed_p_wgd")
  dir.create(fixed_pwgd_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    data.frame(fixed_parameter = "p_wgd", fixed_value = levels(endpoint_plot$pwgd_label)),
    file.path(fixed_pwgd_dir, "fixed_values.tsv")
  )
  save_endpoint_unit_grid_pages(
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_pdf = file.path(fixed_pwgd_dir, "endpoint_burden_ploidy_by_o2_grid.pdf"),
    include_trigger = include_trigger
  )
  save_endpoint_triangle_unit_grid_pages(
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_pdf = file.path(fixed_pwgd_dir, "endpoint_burden_ploidy_triangle_by_o2_grid.pdf"),
    include_trigger = include_trigger
  )
  save_endpoint_unit_individual_plots(
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_dir = file.path(fixed_pwgd_dir, "units"),
    include_trigger = include_trigger
  )
  save_endpoint_triangle_unit_individual_plots(
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_dir = file.path(fixed_pwgd_dir, "triangle_units"),
    include_trigger = include_trigger
  )
  invisible(unit_root)
}

plot_interaction_heatmaps <- function(endpoint_summary,
                                      burden_all,
                                      ploidy_summary,
                                      run_params,
                                      out_dir,
                                      include_trigger = TRUE,
                                      show_treatment_lines = TRUE,
                                      pmiss_axis_label = "post p_mis_base",
                                      condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  ref_group_cols <- if (include_trigger) c("o2_S0", "trigger_burden_cells") else "o2_S0"
  ref <- endpoint_plot %>%
    group_by(across(all_of(ref_group_cols))) %>%
    filter(
      .plot_p_mis_base == min(.plot_p_mis_base, na.rm = TRUE),
      abs(.plot_p_wgd - as.numeric(run_params$p_wgd)) == min(abs(.plot_p_wgd - as.numeric(run_params$p_wgd)), na.rm = TRUE)
    ) %>%
    slice(1) %>%
    ungroup() %>%
    select(
      all_of(ref_group_cols),
      ref_log10_burden = pred_log10_burden_cells,
      ref_mean_ploidy = mean_ploidy
    )
  endpoint_plot <- endpoint_plot %>%
    left_join(ref, by = ref_group_cols) %>%
    mutate(
      delta_log10_burden = pred_log10_burden_cells - ref_log10_burden,
      delta_mean_ploidy = mean_ploidy - ref_mean_ploidy
    )

  heat_theme <- theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 7)
    )
  facet_layer <- if (include_trigger) facet_grid(trigger_label ~ o2_label) else facet_grid(. ~ o2_label)
  endpoint_height <- if (include_trigger) 12 else 3.6
  timecourse_height <- if (include_trigger) 13 else 5.2

  p <- ggplot(endpoint_plot, aes(pmiss_label, pwgd_label, fill = pred_log10_burden_cells)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_layer +
    scale_fill_log10_burden("Endpoint\nlog10 burden") +
    labs(x = pmiss_axis_label, y = "p_wgd relative to fitted (actual p_wgd)") +
    heat_theme
  ggsave(file.path(out_dir, "endpoint_log10_burden_heatmap.pdf"), p, width = 18, height = endpoint_height)

  p <- ggplot(endpoint_plot, aes(pmiss_label, pwgd_label, fill = mean_ploidy)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_layer +
    scale_fill_mean_ploidy("Endpoint\nmean ploidy") +
    labs(x = pmiss_axis_label, y = "p_wgd relative to fitted (actual p_wgd)", fill = "Endpoint\nmean ploidy") +
    heat_theme
  ggsave(file.path(out_dir, "endpoint_mean_ploidy_heatmap.pdf"), p, width = 18, height = endpoint_height)

  p <- ggplot(endpoint_plot, aes(pmiss_label, pwgd_label, fill = delta_log10_burden)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_layer +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") +
    labs(x = pmiss_axis_label, y = "p_wgd relative to fitted (actual p_wgd)", fill = "Delta log10\nburden") +
    heat_theme
  ggsave(file.path(out_dir, "endpoint_delta_log10_burden_heatmap.pdf"), p, width = 18, height = endpoint_height)

  p <- ggplot(endpoint_plot, aes(pmiss_label, pwgd_label, fill = delta_mean_ploidy)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_layer +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") +
    labs(x = pmiss_axis_label, y = "p_wgd relative to fitted (actual p_wgd)", fill = "Delta\nmean ploidy") +
    heat_theme
  ggsave(file.path(out_dir, "endpoint_delta_mean_ploidy_heatmap.pdf"), p, width = 18, height = endpoint_height)

  burden_plot <- add_plot_labels(burden_all, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  burden_plot_full <- complete_timecourse_for_plot(burden_plot, "pred_log10_burden_cells", fill_value = 0)
  treatment_lines <- make_treatment_line_segments(
    add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix),
    row_index_col = "condition_index",
    show_treatment_lines = show_treatment_lines
  )
  condition_breaks <- seq_along(levels(burden_plot$condition_label))
  condition_labels <- levels(burden_plot$condition_label)
  p <- ggplot(burden_plot_full, aes(day, condition_index, fill = pred_log10_burden_cells)) +
    geom_raster() +
    facet_layer +
    scale_fill_log10_burden("log10 burden") +
    scale_y_continuous(breaks = condition_breaks, labels = condition_labels, expand = c(0, 0)) +
    labs(x = "Day", y = "p_wgd relative to fitted (actual p_wgd)") +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 5),
      axis.ticks.y = element_line(linewidth = 0.2),
      strip.text = element_text(size = 7)
    )
  if (nrow(treatment_lines) > 0L) {
    p <- p + geom_segment(
      data = treatment_lines,
      aes(x = treatment_day, xend = treatment_day, y = y0, yend = y1),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.25
    )
  }
  ggsave(file.path(out_dir, "timecourse_log10_burden_heatmap.pdf"), p, width = 20, height = timecourse_height)

  ploidy_plot <- add_plot_labels(ploidy_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  ploidy_plot_full <- complete_timecourse_for_plot(ploidy_plot, "mean_ploidy", fill_value = 0)
  treatment_lines <- make_treatment_line_segments(
    add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix),
    row_index_col = "condition_index",
    show_treatment_lines = show_treatment_lines
  )
  condition_breaks <- seq_along(levels(ploidy_plot$condition_label))
  condition_labels <- levels(ploidy_plot$condition_label)
  p <- ggplot(ploidy_plot_full, aes(day, condition_index, fill = mean_ploidy)) +
    geom_raster() +
    facet_layer +
    scale_fill_mean_ploidy("mean ploidy") +
    scale_y_continuous(breaks = condition_breaks, labels = condition_labels, expand = c(0, 0)) +
    labs(x = "Day", y = "p_wgd relative to fitted (actual p_wgd)", fill = "mean ploidy") +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 5),
      axis.ticks.y = element_line(linewidth = 0.2),
      strip.text = element_text(size = 7)
    )
  if (nrow(treatment_lines) > 0L) {
    p <- p + geom_segment(
      data = treatment_lines,
      aes(x = treatment_day, xend = treatment_day, y = y0, yend = y1),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.25
    )
  }
  ggsave(file.path(out_dir, "timecourse_mean_ploidy_heatmap.pdf"), p, width = 20, height = timecourse_height)

  plot_split_timecourse_heatmaps(
    endpoint_summary = endpoint_summary,
    burden_all = burden_all,
    ploidy_summary = ploidy_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines,
    condition_pmiss_prefix = condition_pmiss_prefix
  )
  plot_split_endpoint_heatmaps(
    endpoint_summary = endpoint_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    condition_pmiss_prefix = condition_pmiss_prefix
  )
  plot_endpoint_unit_grid_heatmaps(
    endpoint_summary = endpoint_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    condition_pmiss_prefix = condition_pmiss_prefix
  )
  plot_timecourse_unit_grid_heatmaps(
    endpoint_summary = endpoint_summary,
    burden_all = burden_all,
    ploidy_summary = ploidy_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines,
    condition_pmiss_prefix = condition_pmiss_prefix
  )
}

run_factorial_interaction <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- resolve_path_value(argv$fit_dir, getwd())
  if (is.null(fit_dir)) stop("Missing required argument: --fit_dir=/path/to/seed_dir")
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  cfg <- prepare_sim_cfg(read_fit_config(fit_dir), argv)
  run_params <- read_run_params(fit_dir, cfg)

  mode <- as.character(.first_non_null_local(argv$mode, "triggered_drug"))
  mode <- tolower(mode)
  if (mode %in% c("drug", "triggered", "triggered-drug")) mode <- "triggered_drug"
  if (mode %in% c("untreated", "factorial", "untreated-factorial")) mode <- "untreated_factorial"
  if (mode %in% c("intermittent", "intermittent-drug", "cycle", "cyclic")) mode <- "intermittent_drug"
  if (!mode %in% c("triggered_drug", "untreated_factorial", "intermittent_drug")) {
    stop("Unknown --mode=", mode, ". Use triggered_drug, untreated_factorial, or intermittent_drug.")
  }

  out_dir <- resolve_path_value(argv$out_dir, getwd())
  if (is.null(out_dir)) {
    out_dir <- file.path(
      fit_dir,
      "simulation",
      switch(
        mode,
        triggered_drug = "invivo_triggered_drug_fixed_pre_1000post",
        untreated_factorial = "invivo_untreated_factorial_2N4N_1e6_100day",
        intermittent_drug = "invivo_intermittent_drug_7on21off_156week"
      )
    )
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  default_initial_burden_cells <- if (identical(mode, "untreated_factorial")) 1e6 else 1e5
  initial_burden_cells <- as_num(argv$initial_burden_cells, as_num(argv$initial_live_cells, default_initial_burden_cells))
  o2_values <- parse_num_list(argv$o2_values, c(0.5, 1, 2, 3, 4, 5))
  sweep_ratios <- c(1e-2, 5e-2, 1e-1, 5e-1, 1, 5, 10, 50, 100)
  pmis_base_values <- parse_num_list(argv$pmis_base_values, 5e-3 * sweep_ratios)
  p_wgd_values <- parse_num_list(argv$p_wgd_values, as.numeric(run_params$p_wgd) * sweep_ratios)
  trigger_burden_values <- parse_num_list(argv$trigger_burden_values, c(5e5, 1e6, 2e6, 5e6, 1e7))
  if (as_bool(argv$smoke, FALSE)) {
    o2_values <- head(o2_values, 2)
    pmis_base_values <- head(pmis_base_values, 2)
    p_wgd_values <- head(p_wgd_values, 2)
    trigger_burden_values <- head(trigger_burden_values, 2)
  }

  pre_horizon_day <- as_num(argv$pre_horizon_day, 1000)
  post_treatment_day <- as_num(argv$post_treatment_day, 1000)
  untreated_horizon_day <- as_num(argv$horizon_day, if (identical(mode, "untreated_factorial")) 100 else 1000)
  intermittent_followup_day <- as_num(argv$intermittent_followup_day, as_num(argv$intermittent_followup_weeks, 156) * 7)
  intermittent_on_day <- as_num(argv$intermittent_on_day, as_num(argv$drug_on_day, 7))
  intermittent_off_day <- as_num(argv$intermittent_off_day, as_num(argv$drug_off_day, 21))
  report_dt <- as_num(argv$report_dt, 1.0)
  trigger_check_dt <- as_num(argv$trigger_check_dt, report_dt)
  make_plots <- as_bool(argv$make_plots, TRUE)
  live_min_cells <- as_num(argv$live_min_cells, 1e-6)

  init_mean <- as_num(argv$initial_ploidy_mean, 3.0)
  init_sd <- as_num(argv$initial_ploidy_sd, 0.4)
  init_min <- as_num(argv$initial_ploidy_min, 1.5)
  init_max <- as_num(argv$initial_ploidy_max, 6.0)
  initial_state_mode <- as.character(.first_non_null_local(
    argv$initial_state_mode,
    if (identical(mode, "untreated_factorial")) "2n4n" else "continuous"
  ))

  design <- if (mode %in% c("triggered_drug", "intermittent_drug")) {
    make_interaction_design(
      run_params = run_params,
      o2_values = o2_values,
      pmis_base_values = pmis_base_values,
      p_wgd_values = p_wgd_values,
      trigger_burden_values = trigger_burden_values,
      initial_burden_cells = initial_burden_cells,
      experiment_name = mode,
      scenario_prefix = if (identical(mode, "intermittent_drug")) "intermittent" else "drug"
    )
  } else {
    make_untreated_design(
      run_params = run_params,
      o2_values = o2_values,
      pmis_base_values = pmis_base_values,
      p_wgd_values = p_wgd_values,
      initial_burden_cells = initial_burden_cells,
      horizon_day = untreated_horizon_day
    )
  }
  design$initial_ploidy_mean <- init_mean
  design$initial_ploidy_sd <- init_sd
  design$initial_ploidy_min <- init_min
  design$initial_ploidy_max <- init_max
  design$initial_state_mode <- initial_state_mode
  design$fit_dir <- fit_dir

  message("Running ", mode, " design: ", nrow(design), " scenarios. Output: ", out_dir)
  message("Fitted baseline p_mis_base: ", signif(as.numeric(run_params$p_mis_base), 6))
  message("Fitted baseline p_wgd: ", signif(as.numeric(run_params$p_wgd), 6))
  message("Requested p_mis_base values: ", paste(signif(sort(unique(pmis_base_values)), 6), collapse = ", "))
  message("Requested p_wgd values: ", paste(signif(sort(unique(p_wgd_values)), 6), collapse = ", "))
  if (mode %in% c("triggered_drug", "intermittent_drug")) {
    message("Pre-treatment uses fitted p_mis_base and fitted p_wgd; only O2_S0 and trigger burden change before treatment.")
    if (identical(mode, "intermittent_drug")) {
      message(
        "Intermittent post-treatment: ", intermittent_on_day, " days on / ",
        intermittent_off_day, " days off for ", intermittent_followup_day,
        " days after first treatment."
      )
    } else {
      message("Pre-treatment search horizon: ", pre_horizon_day, " days; post-treatment follow-up: ", post_treatment_day, " days.")
    }
    message("Pre-treatment trajectories are cached by O2_S0 only.")
  } else {
    message("Untreated mode uses requested O2_S0, p_mis_base, and p_wgd from day 0 for ", untreated_horizon_day, " days.")
  }

  init_state <- make_initial_state_for_mode(
    cfg = cfg,
    total_live_cells = as.numeric(initial_burden_cells),
    initial_state_mode = initial_state_mode,
    init_mean = init_mean,
    init_sd = init_sd,
    init_min = init_min,
    init_max = init_max
  )

  pre_cache <- list()
  burden_rows <- vector("list", nrow(design))
  ploidy_rows <- vector("list", nrow(design))
  design_rows <- vector("list", nrow(design))

  if (identical(mode, "untreated_factorial")) {
    for (i in seq_len(nrow(design))) {
      dr <- design[i, , drop = FALSE]
      message("[", i, "/", nrow(design), "] ", dr$scenario_id)
      rp <- run_params
      rp$o2_S0 <- as.numeric(dr$o2_S0)
      rp$p_wgd <- as.numeric(dr$p_wgd_post)
      rp$p_mis_base <- as.numeric(dr$p_mis_base_post)
      sim <- simulate_cpp_light(
        init_state = init_state,
        run_params = rp,
        cfg = cfg,
        horizon_day = untreated_horizon_day,
        report_dt = report_dt,
        scenario_id = dr$scenario_id,
        experiment = dr$experiment,
        segment = "untreated",
        live_min_cells = live_min_cells
      )
      sim_ann <- annotate_light_timecourses(sim, dr)
      burden_rows[[i]] <- sim_ann$burden
      ploidy_rows[[i]] <- sim_ann$ploidy_summary
      design_rows[[i]] <- dr
    }

    design_out <- bind_rows(design_rows)
    burden_all <- bind_rows(burden_rows)
    ploidy_summary <- bind_rows(ploidy_rows)
    endpoint_summary <- make_endpoint_summary(design_out, burden_all, ploidy_summary)

    write_tsv(design_out, file.path(out_dir, "interaction_design.tsv"))
    write_tsv(burden_all, file.path(out_dir, "burden_timecourse.tsv"))
    write_tsv(ploidy_summary, file.path(out_dir, "ploidy_summary_timecourse.tsv"))
    write_tsv(endpoint_summary, file.path(out_dir, "endpoint_summary.tsv"))

    if (make_plots) {
      plot_interaction_heatmaps(
        endpoint_summary, burden_all, ploidy_summary, run_params, out_dir,
        include_trigger = FALSE,
        show_treatment_lines = FALSE,
        pmiss_axis_label = "p_mis_base",
        condition_pmiss_prefix = "p_miss="
      )
    }

    message("Done. Wrote outputs to: ", normalizePath(out_dir, mustWork = FALSE))
    return(invisible(out_dir))
  }

  pre_grid <- unique(design[, c("o2_S0"), drop = FALSE])
  for (i in seq_len(nrow(pre_grid))) {
    key <- format_value_label(pre_grid$o2_S0[[i]])
    rp <- run_params
    rp$o2_S0 <- as.numeric(pre_grid$o2_S0[[i]])
    rp$p_wgd <- as.numeric(run_params$p_wgd)
    rp$p_mis_base <- as.numeric(run_params$p_mis_base)
    message(
      "[pre ", i, "/", nrow(pre_grid), "] O2=", signif(rp$o2_S0, 6),
      " p_wgd=", signif(rp$p_wgd, 6),
      " p_mis_base=", signif(rp$p_mis_base, 6)
    )
    pre_cache[[key]] <- simulate_cpp_light(
      init_state = init_state,
      run_params = rp,
      cfg = cfg,
      horizon_day = pre_horizon_day,
      report_dt = trigger_check_dt,
      scenario_id = paste0("trajectory_", i),
      experiment = "factorial_interaction_cache",
      segment = "pre",
      live_min_cells = live_min_cells
    )
  }

  for (i in seq_len(nrow(design))) {
    dr <- design[i, , drop = FALSE]
    key <- format_value_label(dr$o2_S0)
    pre <- pre_cache[[key]]
    message("[", i, "/", nrow(design), "] ", dr$scenario_id)

    hit_idx <- which(pre$burden$pred_burden_cells >= as.numeric(dr$trigger_burden_cells))
    if (!length(hit_idx)) {
      dr$status <- "trigger_not_reached"
      dr$trigger_day <- NA_real_
      dr$actual_trigger_burden_cells <- NA_real_
      dr$post_treatment_duration_day <- NA_real_
      dr$sim_end_day <- pre_horizon_day
      pre_ann <- annotate_light_timecourses(pre, dr)
      burden_rows[[i]] <- pre_ann$burden
      ploidy_rows[[i]] <- pre_ann$ploidy_summary
      design_rows[[i]] <- dr
      next
    }

    hit <- hit_idx[[1]]
    trigger_day <- as.numeric(pre$burden$day[[hit]])
    dr$status <- if (identical(mode, "intermittent_drug")) "intermittent_treated" else "treated"
    dr$trigger_day <- trigger_day
    dr$actual_trigger_burden_cells <- as.numeric(pre$burden$pred_burden_cells[[hit]])
    followup_day <- if (identical(mode, "intermittent_drug")) intermittent_followup_day else post_treatment_day
    dr$post_treatment_duration_day <- followup_day
    dr$sim_end_day <- trigger_day + followup_day
    if (identical(mode, "intermittent_drug")) {
      dr$intermittent_on_day <- intermittent_on_day
      dr$intermittent_off_day <- intermittent_off_day
    }

    rp_post <- run_params
    rp_post$o2_S0 <- as.numeric(dr$o2_S0)
    rp_post$p_wgd <- as.numeric(dr$p_wgd_post)
    rp_post$p_mis_base <- as.numeric(dr$p_mis_base_post)
    if (identical(mode, "intermittent_drug")) {
      rp_off <- run_params
      rp_off$o2_S0 <- as.numeric(dr$o2_S0)
      rp_off$p_wgd <- as.numeric(run_params$p_wgd)
      rp_off$p_mis_base <- as.numeric(run_params$p_mis_base)
      post <- simulate_intermittent_post(
        init_state = as.numeric(pre$live_state[hit, ]),
        init_dead_hypoxia_state = as.numeric(pre$dead_hypoxia_state[hit, ]),
        init_dead_buffer_state = as.numeric(pre$dead_buffer_state[hit, ]),
        run_params_on = rp_post,
        run_params_off = rp_off,
        cfg = cfg,
        start_day = trigger_day,
        followup_day = intermittent_followup_day,
        on_day = intermittent_on_day,
        off_day = intermittent_off_day,
        report_dt = report_dt,
        scenario_id = dr$scenario_id,
        experiment = dr$experiment,
        live_min_cells = live_min_cells
      )
    } else {
      post <- simulate_cpp_light(
        init_state = as.numeric(pre$live_state[hit, ]),
        init_dead_hypoxia_state = as.numeric(pre$dead_hypoxia_state[hit, ]),
        init_dead_buffer_state = as.numeric(pre$dead_buffer_state[hit, ]),
        run_params = rp_post,
        cfg = cfg,
        horizon_day = post_treatment_day,
        report_dt = report_dt,
        start_day = trigger_day,
        scenario_id = dr$scenario_id,
        experiment = dr$experiment,
        segment = "post",
        live_min_cells = live_min_cells
      )
    }

    pre_keep <- pre
    pre_keep$burden <- pre_keep$burden %>% filter(day < trigger_day - 1e-9)
    pre_keep$ploidy_summary <- pre_keep$ploidy_summary %>% filter(day < trigger_day - 1e-9)
    pre_ann <- annotate_light_timecourses(pre_keep, dr)
    post_ann <- annotate_light_timecourses(post, dr)

    burden_rows[[i]] <- bind_rows(pre_ann$burden, post_ann$burden)
    ploidy_rows[[i]] <- bind_rows(pre_ann$ploidy_summary, post_ann$ploidy_summary)
    design_rows[[i]] <- dr
  }

  design_out <- bind_rows(design_rows)
  burden_all <- bind_rows(burden_rows)
  ploidy_summary <- bind_rows(ploidy_rows)
  endpoint_summary <- make_endpoint_summary(design_out, burden_all, ploidy_summary)

  write_tsv(design_out, file.path(out_dir, "interaction_design.tsv"))
  write_tsv(burden_all, file.path(out_dir, "burden_timecourse.tsv"))
  write_tsv(ploidy_summary, file.path(out_dir, "ploidy_summary_timecourse.tsv"))
  write_tsv(endpoint_summary, file.path(out_dir, "endpoint_summary.tsv"))

  if (make_plots) {
    plot_interaction_heatmaps(
      endpoint_summary, burden_all, ploidy_summary, run_params, out_dir,
      include_trigger = TRUE,
      show_treatment_lines = TRUE,
      pmiss_axis_label = "post p_mis_base",
      condition_pmiss_prefix = "post p_miss="
    )
  }

  message("Done. Wrote outputs to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(out_dir)
}

if (sys.nframe() == 0) {
  run_factorial_interaction()
}
