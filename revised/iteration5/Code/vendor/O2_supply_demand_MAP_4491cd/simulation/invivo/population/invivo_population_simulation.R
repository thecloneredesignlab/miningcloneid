# In-vivo population simulation products. This module is data-only.

read_run_params <-
function(fit_dir, cfg = NULL) {
    p <- file.path(fit_dir, "best_params.tsv")
    if (!file.exists(p))
        stop("Missing best_params.tsv in: ", fit_dir)
    tab <- read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("parameter", "value") %in% names(tab))) {
        stop("best_params.tsv must contain columns: parameter, value")
    }
    vals <- setNames(as.numeric(tab$value), as.character(tab$parameter))
    needed_common <- c("lam_max", "p_misseg", "k_o_mis", "o2_S0", "kappa_O", "eta_o2", "alpha_o2", "gamma_growth", "mu_hp",
        "gamma_mu", "O2_crit", "n_O", "k_clear")
    needed <- c(needed_common, "buffer_smax", "buffer_beta", "buffer_n_exp", "p_wgd")
    miss <- setdiff(needed, names(vals))
    if (length(miss) > 0) {
        stop("best_params.tsv missing parameters: ", paste(miss, collapse = ", "))
    }
    out <- as.list(vals[needed])
    extra_init_mult <- grep("^init_mult_", names(vals), value = TRUE)
    if (length(extra_init_mult) > 0L) {
        for (nm in extra_init_mult) out[[nm]] <- vals[[nm]]
    }
    p_mis_base_val <- if ("p_mis_base" %in% names(vals))
        vals[["p_mis_base"]]
    else NULL
    o2_min_val <- if ("o2_min" %in% names(vals))
        vals[["o2_min"]]
    else NULL
    if (!is.null(p_mis_base_val) && is.finite(p_mis_base_val))
        out$p_mis_base <- as.numeric(p_mis_base_val)
    if (!is.null(o2_min_val) && is.finite(o2_min_val))
        out$o2_min <- as.numeric(o2_min_val)
    if ("rho_2N" %in% names(vals) && is.finite(vals[["rho_2N"]]) && vals[["rho_2N"]] > 0) {
        out$rho_2N <- vals[["rho_2N"]]
    }
    if ("c_vol_2N_mm3" %in% names(vals) && is.finite(vals[["c_vol_2N_mm3"]]) && vals[["c_vol_2N_mm3"]] > 0) {
        out$c_vol_2N_mm3 <- vals[["c_vol_2N_mm3"]]
    }
    fit_treatment_use <- isTRUE(.first_non_null_local(cfg$fit_treatment, FALSE))
    if (fit_treatment_use) {
        miss_tx <- setdiff(c("alpha", "gamma"), names(vals))
        if (length(miss_tx) > 0) {
            stop("best_params.tsv missing treatment parameters while fit_treatment=TRUE: ", paste(miss_tx, collapse = ", "))
        }
        out$alpha <- vals[["alpha"]]
        out$gamma <- vals[["gamma"]]
    }
    else {
        out$alpha <- if ("alpha" %in% names(vals) && is.finite(vals[["alpha"]]))
            vals[["alpha"]]
        else 0
        out$gamma <- if ("gamma" %in% names(vals) && is.finite(vals[["gamma"]]))
            vals[["gamma"]]
        else 1
    }
    out$tau_O2 <- if ("tau_O2" %in% names(vals) && is.finite(vals[["tau_O2"]]) && vals[["tau_O2"]] > 0)
        vals[["tau_O2"]]
    else NULL
    normalize_run_params_common(out, cfg = cfg)
}

simulate_one_full <-
function(run_params, scenario, cfg, report_dt = 1) {
    model_core <- build_model_core(cfg = cfg)
    grid_pre <- model_core$grid_pre
    init_state <- if (scenario$cohort == "2N")
        model_core$init_state_2N
    else model_core$init_state_4N
    init_mult_name <- .first_non_null_local(scenario$init_mult_param, harvest_init_natural_param_name(scenario$harvest))
    init_mult <- suppressWarnings(as.numeric(.first_non_null_local(run_params[[init_mult_name]], 1)))
    if (!is.finite(init_mult) || init_mult <= 0) {
        log_name <- .first_non_null_local(scenario$log_init_mult_param, harvest_init_log_param_name(scenario$harvest))
        init_mult <- exp(as.numeric(.first_non_null_local(run_params[[log_name]], 0)))
    }
    if (!is.finite(init_mult) || init_mult <= 0)
        init_mult <- 1
    init_state <- as.numeric(init_state) * init_mult
    sim_end_day <- as.numeric(scenario$sim_end_day)
    full_steps <- 0:as.integer(round(sim_end_day/cfg$DT))
    full_days <- as.numeric(full_steps) * cfg$DT
    keep_days <- sort(unique(c(0, sim_end_day, as.numeric(scenario$obs_days), seq(0, sim_end_day, by = report_dt))))
    keep_days <- keep_days[is.finite(keep_days) & keep_days >= 0 & keep_days <= sim_end_day]
    keep_steps <- sort(unique(as.integer(round(keep_days/cfg$DT))))
    obs_steps_cpp <- as.integer(full_steps)
    sim_end_step_cpp <- max(obs_steps_cpp)
    o2_s0_upper_use <- as.numeric(.first_non_null_local(cfg$o2_S0_upper_bound, 5))
    if (!is.finite(o2_s0_upper_use) || o2_s0_upper_use <= 0)
        o2_s0_upper_use <- 5
    o2_S0_use <- as.numeric(.first_non_null_local(run_params$o2_S0, cfg$o2_S0_init, 0.5))
    if (!is.finite(o2_S0_use) || o2_S0_use < 0)
        o2_S0_use <- 0.5
    o2_S0_use <- min(max(o2_S0_use, 0), o2_s0_upper_use)
    kappa_O_use <- as.numeric(.first_non_null_local(run_params$kappa_O, cfg$kappa_O_init, 1))
    if (!is.finite(kappa_O_use) || kappa_O_use <= 0)
        kappa_O_use <- 1
    eta_o2_use <- as.numeric(.first_non_null_local(run_params$eta_o2, cfg$eta_o2_init, 1))
    if (!is.finite(eta_o2_use) || eta_o2_use <= 0)
        eta_o2_use <- 1
    O2_crit_use <- as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1))
    if (!is.finite(O2_crit_use) || O2_crit_use < 0)
        O2_crit_use <- 1
    o2_Nref_use <- as.numeric(.first_non_null_local(cfg$o2_Nref, cfg$init_total_size, 1e+06))
    if (!is.finite(o2_Nref_use) || o2_Nref_use <= 0)
        o2_Nref_use <- 1e+06
    o2_min_use <- as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0))
    if (!is.finite(o2_min_use) || o2_min_use < 0)
        o2_min_use <- 0
    o2_min_use <- min(max(o2_min_use, 0), o2_s0_upper_use)
    tau_O2_use <- as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2))
    if (!is.finite(tau_O2_use) || tau_O2_use <= 0)
        tau_O2_use <- 2
    vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg))
    burden_floor <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 0)
    o2_growth_use <- isTRUE(.first_non_null_local(cfg$O2_growth, TRUE))
    buffer_smax_use <- as.numeric(.first_non_null_local(run_params$buffer_smax, cfg$buffer_smax_init, 1))
    if (!is.finite(buffer_smax_use))
        buffer_smax_use <- 1
    buffer_beta_use <- as.numeric(.first_non_null_local(run_params$buffer_beta, cfg$buffer_beta_init, 0))
    if (!is.finite(buffer_beta_use))
        buffer_beta_use <- 0
    buffer_n_exp_use <- as.numeric(.first_non_null_local(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1))
    if (!is.finite(buffer_n_exp_use) || buffer_n_exp_use <= 0)
        buffer_n_exp_use <- 1
    sim_cpp <- cpp_o2simps_simulate_one(list(init_state = as.numeric(init_state), N0min = as.integer(cfg$N_MIN), N0max = as.integer(cfg$N_MAX),
        N1min = as.integer(cfg$N_MIN), N1max = as.integer(cfg$N_MAX), obs_steps = as.integer(obs_steps_cpp), sim_end_step = as.integer(sim_end_step_cpp),
        DT = as.numeric(cfg$DT), dose = as.numeric(scenario$dose), dose_ref = as.numeric(cfg$dose_ref), treat_day = as.numeric(scenario$treat_day),
        fit_treatment = isTRUE(cfg$fit_treatment), alpha = as.numeric(.first_non_null_local(run_params$alpha, 0)), gamma = as.numeric(.first_non_null_local(run_params$gamma,
            1)), tx_mult_min = as.numeric(cfg$tx_mult_min), crowding_enabled = isTRUE(cfg$Crowding), crowding = as.character(cfg$crowding),
        K = as.numeric(cfg$K), min_pop = as.numeric(cfg$min_pop), O2_crit = as.numeric(O2_crit_use), o2_feedback = isTRUE(.first_non_null_local(cfg$o2_burden_feedback,
            TRUE)), o2_S0 = as.numeric(o2_S0_use), kappa_O = as.numeric(kappa_O_use), tau_O2 = as.numeric(tau_O2_use), o2_Nref = as.numeric(o2_Nref_use),
        o2_min = as.numeric(o2_min_use), o2_S0_upper_bound = as.numeric(o2_s0_upper_use), eta_o2 = as.numeric(eta_o2_use),
        o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01)), o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct,
            0.005)), o2_cache_profile = isTRUE(.first_non_null_local(cfg$o2_cache_profile, FALSE)), O2_growth = isTRUE(o2_growth_use),
        lam_max = as.numeric(run_params$lam_max), p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base,
            cfg$p_mis_base_init, 1e-05)), p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0)), k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis,
            50)), p_wgd = as.numeric(.first_non_null_local(run_params$p_wgd, 0)), boundary = as.character(.first_non_null_local(run_params$boundary,
            "drop")), eps_tail = as.numeric(1e-08), buffer_smax = as.numeric(buffer_smax_use), buffer_beta = as.numeric(buffer_beta_use),
        buffer_n_exp = as.numeric(buffer_n_exp_use), N_unit = as.integer(cfg$N_UNIT), beta_size = 0, alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2,
            cfg$alpha_o2_init, 0.5)), gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init,
            2)), mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 0.001)), gamma_mu = as.numeric(.first_non_null_local(run_params$gamma_mu,
            cfg$gamma_mu_init, 1)), n_O = as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1)), ploidy_O2_death = assert_canonical_ploidy_o2_death_mode(.first_non_null_local(cfg$ploidy_O2_death,
            run_params$ploidy_O2_death, "diploid_NULL")), start_with = assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with,
            "ploidy")), k_clear = as.numeric(.first_non_null_local(run_params$k_clear, cfg$k_clear_init, 0.001)), vol_by_N = as.numeric(vol_by_N),
        burden_floor = as.numeric(burden_floor), return_full_trajectory = TRUE))
    live_state_obs <- sim_cpp$live_state_obs
    if (is.null(live_state_obs) || length(live_state_obs) == 0) {
        return(list(burden = data.frame(), ploidy = data.frame()))
    }
    live_state_obs <- as.matrix(live_state_obs)
    if (!identical(dim(live_state_obs), c(length(obs_steps_cpp), length(grid_pre)))) {
        stop("cpp_o2simps_simulate_one returned live_state_obs with unexpected shape: got ", paste(dim(live_state_obs), collapse = "x"),
            ", expected ", length(obs_steps_cpp), "x", length(grid_pre))
    }
    live_pop_by_step <- rowSums(live_state_obs, na.rm = TRUE)
    live_frac_obs <- live_state_obs/pmax(live_pop_by_step, 1e-12)
    ploidy_rows_full <- data.frame(harvest = scenario$harvest, cohort = scenario$cohort, dose = scenario$dose, day = rep(full_days,
        each = length(grid_pre)), step = rep(obs_steps_cpp, each = length(grid_pre)), N = rep(as.numeric(grid_pre), times = length(obs_steps_cpp)),
        fraction = as.numeric(t(live_frac_obs)), pop = rep(as.numeric(live_pop_by_step), each = length(grid_pre)))
    burden_by_day_full <- data.frame(day = full_days, step = obs_steps_cpp, pred_burden_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_total_obs,
        sim_cpp$Ntot_obs)), pred_burden_live_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_live_obs, sim_cpp$Ntot_obs)),
        pred_burden_dead_hypoxia_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_hypoxia_obs, rep(0, length(obs_steps_cpp)))),
        pred_burden_dead_buffer_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_buffer_obs, rep(0, length(obs_steps_cpp)))),
        pred_burden_dead_total_cells = as.numeric(.first_non_null_local(sim_cpp$Ntot_dead_total_obs, rep(0, length(obs_steps_cpp)))),
        pred_burden_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_total_obs, sim_cpp$Vmm3_obs)), pred_burden_live_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_live_obs,
            sim_cpp$Vmm3_obs)), pred_burden_dead_hypoxia_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_hypoxia_obs,
            rep(0, length(obs_steps_cpp)))), pred_burden_dead_buffer_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_buffer_obs,
            rep(0, length(obs_steps_cpp)))), pred_burden_dead_total_volume_mm3 = as.numeric(.first_non_null_local(sim_cpp$Vmm3_dead_total_obs,
            rep(0, length(obs_steps_cpp)))), pred_o2_target_pct = as.numeric(.first_non_null_local(sim_cpp$O2_target_obs,
            rep(NA_real_, length(obs_steps_cpp)))), pred_o2_pct = as.numeric(.first_non_null_local(sim_cpp$O2_eff_obs, rep(NA_real_,
            length(obs_steps_cpp))))) %>% arrange(step)
    burden_by_day_full$pred_o2_lag_gap_pct <- as.numeric(burden_by_day_full$pred_o2_target_pct - burden_by_day_full$pred_o2_pct)
    burden_by_day <- burden_by_day_full %>% filter(step %in% keep_steps)
    ploidy_rows <- ploidy_rows_full %>% filter(step %in% keep_steps)
    obs_steps <- as.integer(round(as.numeric(scenario$obs_days)/cfg$DT))
    obs_map <- setNames(as.numeric(scenario$obs_burden), as.character(obs_steps))
    burden_by_day$obs_burden <- as.numeric(obs_map[as.character(burden_by_day$step)])
    burden_rows <- burden_by_day %>% mutate(harvest = scenario$harvest, cohort = scenario$cohort, dose = scenario$dose, treat_day = scenario$treat_day,
        pred_burden = pred_burden_volume_mm3) %>% select(harvest, cohort, dose, treat_day, step, day, pred_burden, pred_burden_volume_mm3,
        pred_burden_cells, pred_burden_live_volume_mm3, pred_burden_dead_hypoxia_volume_mm3, pred_burden_dead_buffer_volume_mm3,
        pred_burden_dead_total_volume_mm3, pred_burden_live_cells, pred_burden_dead_hypoxia_cells, pred_burden_dead_buffer_cells,
        pred_burden_dead_total_cells, pred_o2_target_pct, pred_o2_pct, pred_o2_lag_gap_pct, obs_burden)
    list(burden = as.data.frame(burden_rows), ploidy = as.data.frame(select(ploidy_rows, harvest, cohort, dose, day, N, fraction)))
}

simulate_one_full_horizon <-
function(run_params, scenario, cfg, horizon_day, report_dt = 1) {
    sc <- scenario
    sc$sim_end_day <- as.numeric(max(horizon_day, 0))
    simulate_one_full(run_params, sc, cfg, report_dt = report_dt)
}

make_canonical_initial_cohort_scenarios <-
function(horizon_day, cfg) {
    horizon_day <- as.numeric(horizon_day)
    if (!is.finite(horizon_day) || horizon_day < 0)
        horizon_day <- 0
    dt_use <- as.numeric(.first_non_null_local(cfg$DT, 1))
    if (!is.finite(dt_use) || dt_use <= 0)
        dt_use <- 1
    lapply(c("2N", "4N"), function(cohort_i) {
        list(harvest = paste0("canonical_", cohort_i), cohort = cohort_i, dose = 0, treat_day = horizon_day + dt_use, sim_end_day = horizon_day,
            obs_days = numeric(0), obs_burden = numeric(0))
    })
}

add_legacy_normalized_burden_columns <-
function(burden_all) {
    burden_all %>% group_by(harvest, cohort, dose) %>% arrange(day, .by_group = TRUE) %>% group_modify(function(df, .y) {
        pred_delta <- df$pred_burden - df$pred_burden[[1]]
        pred_scale <- max(abs(pred_delta), na.rm = TRUE)
        df$pred_norm <- if (is.finite(pred_scale) && pred_scale > 0)
            pred_delta/pred_scale
        else pred_delta
        obs_vals <- df$obs_burden[is.finite(df$obs_burden)]
        if (length(obs_vals) > 0) {
            obs_delta <- df$obs_burden - obs_vals[[1]]
            obs_scale <- max(abs(obs_delta), na.rm = TRUE)
            df$obs_norm <- if (is.finite(obs_scale) && obs_scale > 0)
                obs_delta/obs_scale
            else obs_delta
        }
        else {
            df$obs_norm <- NA_real_
        }
        df
    }) %>% ungroup()
}

generate_invivo_horizon_tables <-
function(run_params, scenarios, cfg, horizon_day, report_dt = 1) {
    sim_list <- lapply(scenarios, function(sc) {
        simulate_one_full_horizon(run_params, sc, cfg, horizon_day = horizon_day, report_dt = report_dt)
    })
    burden_all <- bind_rows(lapply(sim_list, `[[`, "burden")) %>% filter(day <= horizon_day + 1e-09)
    ploidy_all <- bind_rows(lapply(sim_list, `[[`, "ploidy")) %>% filter(day <= horizon_day + 1e-09)
    if (!nrow(burden_all) || !nrow(ploidy_all))
        stop("No prediction output for horizon ", horizon_day)
    burden_all <- add_legacy_normalized_burden_columns(burden_all)
    ploidy_mean <- compute_ploidy_weighted_mean(ploidy_all, cfg)
    chr_density_n_min <- as.integer(cfg$N_MIN)
    chr_density_n_max <- as.integer(cfg$N_MAX)
    chr_density_bin_width <- 5L
    chr_density_bins <- data.frame(chr_bin_lower = seq.int(chr_density_n_min, chr_density_n_max, by = chr_density_bin_width),
        stringsAsFactors = FALSE) %>% mutate(chr_bin_upper = chr_bin_lower + chr_density_bin_width - 1L, chr_bin_mid = (as.numeric(chr_bin_lower) +
        as.numeric(chr_bin_upper))/2)
    chr_density_df <- ploidy_all %>% filter(as.character(cohort) %in% c("2N", "4N"), is.finite(day), is.finite(N), is.finite(fraction)) %>%
        transmute(cohort = factor(as.character(cohort), levels = c("2N", "4N")), day = as.numeric(day), N = as.integer(round(as.numeric(N))),
            sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__"), fraction = pmax(as.numeric(fraction),
                0)) %>% filter(N >= chr_density_n_min, N <= chr_density_n_max) %>% mutate(chr_bin_lower = chr_density_n_min +
        ((N - chr_density_n_min)%/%chr_density_bin_width) * chr_density_bin_width) %>% group_by(cohort, day, sample_id, chr_bin_lower) %>%
        summarise(bin_probability = sum(fraction, na.rm = TRUE), .groups = "drop") %>% left_join(chr_density_bins, by = "chr_bin_lower") %>%
        group_by(cohort, day, chr_bin_lower, chr_bin_upper, chr_bin_mid) %>% summarise(density = mean(bin_probability, na.rm = TRUE),
        .groups = "drop") %>% complete(cohort, day, tidyr::nesting(chr_bin_lower, chr_bin_upper, chr_bin_mid), fill = list(density = 0)) %>%
        mutate(density = pmax(pmin(as.numeric(density), 1), 0))
    burden_decomp <- burden_all %>% transmute(cohort = factor(as.character(cohort), levels = c("2N", "4N")), harvest = as.character(harvest),
        dose = as.numeric(dose), day = as.numeric(day), burden_live = as.numeric(.first_non_null_local(pred_burden_live_volume_mm3,
            pred_burden)), burden_dead_hypoxia = as.numeric(.first_non_null_local(pred_burden_dead_hypoxia_volume_mm3, 0)),
        burden_dead_buffer = as.numeric(.first_non_null_local(pred_burden_dead_buffer_volume_mm3, 0)), burden_dead_total = as.numeric(.first_non_null_local(pred_burden_dead_total_volume_mm3,
            .first_non_null_local(pred_burden_dead_hypoxia_volume_mm3, 0) + .first_non_null_local(pred_burden_dead_buffer_volume_mm3,
                0))), burden_total = as.numeric(pred_burden))
    death_ratio <- burden_decomp %>% mutate(sample_id = paste(harvest, as.character(cohort), format(dose, trim = TRUE, scientific = FALSE),
        sep = "__"), death_ratio = pmax(0, pmin(1, burden_dead_total/pmax(burden_total, 1e-12))))
    resource_death_fraction <- burden_decomp %>% mutate(sample_id = paste(harvest, as.character(cohort), format(dose, trim = TRUE,
        scientific = FALSE), sep = "__"), resource_death_fraction = ifelse(burden_dead_total > 0, burden_dead_hypoxia/pmax(burden_dead_total,
        1e-12), 0), resource_death_fraction = pmax(0, pmin(1, as.numeric(resource_death_fraction))))
    canonical_scenarios <- make_canonical_initial_cohort_scenarios(horizon_day, cfg)
    canonical_sim <- lapply(canonical_scenarios, function(sc) {
        simulate_one_full_horizon(run_params, sc, cfg, horizon_day = horizon_day, report_dt = report_dt)
    })
    canonical_burden <- bind_rows(lapply(canonical_sim, `[[`, "burden"))
    canonical_ploidy <- bind_rows(lapply(canonical_sim, `[[`, "ploidy"))
    canonical_ms <- compute_missegregation_probability_timecourse(canonical_ploidy, canonical_burden, run_params)
    cin_rows <- generate_population_average_cin_rows(canonical_ms, target_day = horizon_day)
    list(burden = burden_all, ploidy = ploidy_all, ploidy_weighted_mean = ploidy_mean, chromosome_density = chr_density_df,
        death_ratio = death_ratio, resource_death_fraction = resource_death_fraction, population_average_cin = cin_rows)
}

generate_invivo_observed_tables <-
function(run_params, scenarios, cfg, report_dt = 1) {
    sim_list <- lapply(scenarios, function(sc) simulate_one_full(run_params, sc, cfg, report_dt = report_dt))
    burden_all <- bind_rows(lapply(sim_list, `[[`, "burden"))
    ploidy_all <- bind_rows(lapply(sim_list, `[[`, "ploidy"))
    if (!nrow(burden_all) || !nrow(ploidy_all))
        stop("No simulation output generated; check fit/data configuration.")
    burden_all <- add_legacy_normalized_burden_columns(burden_all)
    burden_decomp <- burden_all %>% mutate(burden_live = as.numeric(.first_non_null_local(pred_burden_live_volume_mm3, pred_burden)),
        burden_dead_hypoxia = as.numeric(.first_non_null_local(pred_burden_dead_hypoxia_volume_mm3, 0)), burden_dead_buffer = as.numeric(.first_non_null_local(pred_burden_dead_buffer_volume_mm3,
            0)), burden_dead_total = as.numeric(.first_non_null_local(pred_burden_dead_total_volume_mm3, burden_dead_hypoxia +
            burden_dead_buffer)), burden_total = as.numeric(pred_burden)) %>% select(harvest, cohort, dose, day, burden_live,
        burden_dead_hypoxia, burden_dead_buffer, burden_dead_total, burden_total)
    ploidy_mean <- compute_ploidy_weighted_mean(ploidy_all, cfg)
    misseg <- compute_missegregation_probability_timecourse(ploidy_all, burden_all, run_params)
    terminal <- build_terminal_ploidy_compare_df(scenarios, ploidy_all, cfg)
    o2_lag <- burden_all %>% filter(is.finite(pred_o2_target_pct), is.finite(pred_o2_pct)) %>% transmute(harvest = as.character(harvest),
        cohort = as.character(cohort), dose = as.numeric(dose), day = as.numeric(day), sample_id = paste(as.character(harvest),
            as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__"), o2_target_pct = as.numeric(clip(pred_o2_target_pct,
            0, 5)), o2_eff_pct = as.numeric(clip(pred_o2_pct, 0, 5)), o2_lag_gap_pct = as.numeric(clip(pred_o2_target_pct,
            0, 5) - clip(pred_o2_pct, 0, 5)))
    burden_vs_o2 <- burden_all %>% filter(is.finite(pred_burden), is.finite(pred_o2_pct)) %>% transmute(harvest = as.character(harvest),
        cohort = as.character(cohort), dose = as.numeric(dose), day = as.numeric(day), burden_mm3 = as.numeric(pred_burden),
        o2_pct = as.numeric(clip(pred_o2_pct, 0, 5)), sample_id = paste(as.character(harvest), as.character(cohort), format(as.numeric(dose),
            trim = TRUE, scientific = FALSE), sep = "__"))
    list(burden_timecourse = burden_all, ploidy_timecourse = ploidy_all, burden_live_dead_decomposition = burden_decomp,
        ploidy_weighted_mean_timecourse = ploidy_mean, missegregation_probability_timecourse = misseg, terminal_ploidy_observed_vs_predicted = terminal,
        o2_lag_timecourse = o2_lag, predict_burden_vs_o2 = burden_vs_o2)
}
