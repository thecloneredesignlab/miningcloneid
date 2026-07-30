# In-vivo functional_response simulation products. This module is data-only.

make_reference_ploidy_states <-
function(cfg, min_multiple = 1.5, max_multiple = 5, step = 0.5) {
    n_unit <- suppressWarnings(as.numeric(cfg$N_UNIT))
    if (!is.finite(n_unit) || n_unit <= 0)
        n_unit <- 22
    ref_state_mult <- seq(min_multiple, max_multiple, by = step)
    ref_state_label <- ifelse(abs(ref_state_mult - round(ref_state_mult)) < 1e-08, paste0(as.integer(round(ref_state_mult)),
        "N"), paste0(format(ref_state_mult, trim = TRUE, nsmall = 1), "N"))
    data.frame(cohort = ref_state_label, ploidy_multiple = ref_state_mult, N_ref = as.numeric(ref_state_mult * n_unit), stringsAsFactors = FALSE)
}

generate_invivo_functional_response_tables <-
function(run_params, cfg) {
    start_with_mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
    o2_plot_min <- 0
    o2_plot_max <- 5
    O2_crit_use <- as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1))
    if (!is.finite(O2_crit_use) || O2_crit_use < 0)
        O2_crit_use <- 1
    o2_grid <- make_o2_crit_adaptive_grid(O2_crit_use, o2_plot_min, o2_plot_max, base_step = 0.02, dense_n = 220L)
    state_plot_min <- if (identical(start_with_mode, "chr_number"))
        as.numeric(cfg$N_MIN)
    else 0
    state_plot_max <- if (identical(start_with_mode, "chr_number"))
        as.numeric(cfg$N_MAX)
    else 10
    state_grid_dense <- if (identical(start_with_mode, "chr_number")) {
        seq(state_plot_min, state_plot_max, by = 1)
    }
    else {
        seq(state_plot_min, state_plot_max, by = 0.05)
    }
    o2_levels_ploidy <- make_o2_crit_reference_levels(O2_crit_use, o2_plot_min, o2_plot_max, coarse_step = 0.5)
    ref_df <- data.frame(cohort = c("2N", "4N"), N_ref = as.numeric(c(2 * cfg$N_UNIT, 4 * cfg$N_UNIT)), stringsAsFactors = FALSE)
    ref_df_multi <- make_reference_ploidy_states(cfg)
    o2_growth_use <- isTRUE(.first_non_null_local(cfg$O2_growth, TRUE))
    gamma_mu_use <- pmax(as.numeric(run_params$gamma_mu), 1e-12)
    ploidy_O2_death_use <- assert_canonical_ploidy_o2_death_mode(.first_non_null_local(cfg$ploidy_O2_death, run_params$ploidy_O2_death,
        "diploid_NULL"))
    n_O <- as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1))
    if (!is.finite(n_O) || n_O < 0)
        stop("run_params$n_O must be finite and >= 0.")
    mu_hp_use <- pmax(as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 0.001)), 0)
    boundary_mode_use <- as.character(.first_non_null_local(run_params$boundary, "drop"))
    p_wgd_use <- as.numeric(.first_non_null_local(run_params$p_wgd, 0))
    if (!is.finite(p_wgd_use))
        p_wgd_use <- 0
    buffer_smax_use <- as.numeric(.first_non_null_local(run_params$buffer_smax, cfg$buffer_smax_init, 1))
    if (!is.finite(buffer_smax_use))
        buffer_smax_use <- 1
    buffer_beta_use <- as.numeric(.first_non_null_local(run_params$buffer_beta, cfg$buffer_beta_init, 0))
    if (!is.finite(buffer_beta_use))
        buffer_beta_use <- 0
    buffer_n_exp_use <- as.numeric(.first_non_null_local(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1))
    if (!is.finite(buffer_n_exp_use) || buffer_n_exp_use <= 0)
        buffer_n_exp_use <- 1
    n_state <- as.integer(cfg$N_MAX) - as.integer(cfg$N_MIN) + 1L
    .require_cpp_o2simps_fn("cpp_o2simps_build_G_for_o2_triplet")
    functional_rate_cache <- new.env(parent = emptyenv())
    get_functional_rate_curves <- function(O2) {
        O2_use <- .assert_o2_pct(as.numeric(O2), label = "O2")
        key <- sprintf("%.6f", O2_use)
        if (!exists(key, envir = functional_rate_cache, inherits = FALSE)) {
            tri <- cpp_o2simps_build_G_for_o2_triplet(O2 = as.numeric(O2_use), O2_crit = as.numeric(O2_crit_use), N0min = as.integer(cfg$N_MIN),
                N0max = as.integer(cfg$N_MAX), N1min = as.integer(cfg$N_MIN), N1max = as.integer(cfg$N_MAX), lam_max = as.numeric(run_params$lam_max),
                p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init,
                  1e-05)), p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0)), k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis,
                  50)), p_wgd = as.numeric(p_wgd_use), boundary = as.character(boundary_mode_use), eps_tail = 1e-08, buffer_smax = as.numeric(buffer_smax_use),
                buffer_beta = as.numeric(buffer_beta_use), buffer_n_exp = as.numeric(buffer_n_exp_use), N_unit = as.integer(cfg$N_UNIT),
                beta_size = 0, O2_growth = isTRUE(o2_growth_use), alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2,
                  cfg$alpha_o2_init, 0.5)), gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init,
                  2)), mu_hp = as.numeric(mu_hp_use), gamma_mu = as.numeric(gamma_mu_use), n_O = as.numeric(n_O), ploidy_O2_death = as.character(ploidy_O2_death_use))
            curve_names <- c("dead_buffer_rate", "misseg_nonviable_rate", "boundary_dropped_rate", "misseg_nonviable_division_prob",
                "misseg_nonviable_daughters_per_division")
            curves <- setNames(vector("list", length(curve_names)), curve_names)
            for (nm in curve_names) {
                vals <- as.numeric(tri[[nm]])
                if (length(vals) != n_state)
                  stop("Unexpected functional rate length for ", nm, ": ", length(vals))
                curves[[nm]] <- vals
            }
            curves$misseg_nonviable_daughter_fraction <- pmax(pmin(0.5 * curves$misseg_nonviable_daughters_per_division,
                1), 0)
            assign(key, curves, envir = functional_rate_cache)
        }
        get(key, envir = functional_rate_cache, inherits = FALSE)
    }
    compute_rate_components <- function(O2, N) {
        O2_vec <- .assert_o2_pct(as.numeric(O2), label = "O2")
        N_vec <- as.numeric(N)
        n_out <- max(length(O2_vec), length(N_vec))
        if (!(length(O2_vec) %in% c(1L, n_out) && length(N_vec) %in% c(1L, n_out))) {
            stop("O2 and N must have compatible lengths in compute_rate_components().")
        }
        O2_vec <- rep_len(O2_vec, n_out)
        N_vec <- rep_len(N_vec, n_out)
        proliferation_rate <- as.numeric(.lambda_eff_of_O2(O2 = O2_vec, run_params = run_params, N = N_vec, O2_crit = O2_crit_use,
            O2_growth = o2_growth_use))
        death_rate <- as.numeric(.mu_eff_of_O2(O2 = O2_vec, run_params = run_params, N = N_vec, O2_crit = O2_crit_use))
        curve_cols <- c("dead_buffer_rate", "misseg_nonviable_rate", "boundary_dropped_rate", "misseg_nonviable_division_prob",
            "misseg_nonviable_daughters_per_division", "misseg_nonviable_daughter_fraction")
        values <- setNames(lapply(curve_cols, function(x) rep(NA_real_, n_out)), curve_cols)
        row_groups <- split(seq_len(n_out), sprintf("%.6f", O2_vec))
        for (row_idx in row_groups) {
            curves <- get_functional_rate_curves(O2_vec[[row_idx[[1L]]]])
            state_idx <- as.integer(round(N_vec[row_idx])) - as.integer(cfg$N_MIN) + 1L
            valid <- is.finite(state_idx) & state_idx >= 1L & state_idx <= n_state
            if (any(valid)) {
                for (nm in curve_cols) values[[nm]][row_idx[valid]] <- curves[[nm]][state_idx[valid]]
            }
        }
        data.frame(O2 = O2_vec, N = N_vec, proliferation_rate = pmax(proliferation_rate, 0), death_rate = pmax(death_rate,
            0), buffer_death_rate = pmax(values$dead_buffer_rate, 0), buffer_death_per_division = pmax(values$dead_buffer_rate,
            0)/pmax(proliferation_rate, 1e-12), misseg_nonviable_rate = pmax(values$misseg_nonviable_rate, 0), boundary_dropped_rate = pmax(values$boundary_dropped_rate,
            0), misseg_nonviable_division_prob = pmax(pmin(values$misseg_nonviable_division_prob, 1), 0), misseg_nonviable_daughters_per_division = pmax(pmin(values$misseg_nonviable_daughters_per_division,
            2), 0), misseg_nonviable_daughter_fraction = pmax(pmin(values$misseg_nonviable_daughter_fraction, 1), 0), net_growth_rate = proliferation_rate -
            death_rate, stringsAsFactors = FALSE)
    }
    build_o2_curve <- function(refs, include_multiple = FALSE) {
        dplyr::bind_rows(lapply(seq_len(nrow(refs)), function(i) {
            N_ref <- refs$N_ref[[i]]
            ms_rate <- as.numeric(.pmisseg_of_O2(O2 = o2_grid, run_params = run_params, N = N_ref, O2_crit = O2_crit_use))
            rate_df <- compute_rate_components(o2_grid, rep(N_ref, length(o2_grid)))
            prefix <- data.frame(oxygen_pct = o2_grid, cohort = refs$cohort[[i]], stringsAsFactors = FALSE)
            if (isTRUE(include_multiple)) {
                prefix$ploidy_multiple <- refs$ploidy_multiple[[i]]
                prefix$N_ref <- N_ref
            }
            suffix <- data.frame(ms_rate = pmax(ms_rate, 0), proliferation_rate = rate_df$proliferation_rate, death_rate = rate_df$death_rate,
                buffer_death_rate = rate_df$buffer_death_rate, buffer_death_per_division = rate_df$buffer_death_per_division,
                misseg_nonviable_rate = rate_df$misseg_nonviable_rate, boundary_dropped_rate = rate_df$boundary_dropped_rate,
                misseg_nonviable_division_prob = rate_df$misseg_nonviable_division_prob, misseg_nonviable_daughters_per_division = rate_df$misseg_nonviable_daughters_per_division,
                misseg_nonviable_daughter_fraction = rate_df$misseg_nonviable_daughter_fraction, net_growth_rate = rate_df$net_growth_rate,
                stringsAsFactors = FALSE)
            out <- cbind(prefix, suffix)
            if (!isTRUE(include_multiple))
                out$N_ref <- N_ref
            out
        }))
    }
    o2_curve <- build_o2_curve(ref_df, include_multiple = FALSE)
    o2_curve_multi <- build_o2_curve(ref_df_multi, include_multiple = TRUE)
    N_states <- seq.int(as.integer(cfg$N_MIN), as.integer(cfg$N_MAX))
    ploidy_grid <- N_states/as.numeric(cfg$N_UNIT)
    ratio <- (2 * as.numeric(cfg$N_UNIT))/pmax(as.numeric(N_states), 1e-12)
    viability <- pmax(pmin(buffer_smax_use * exp(-buffer_beta_use * pmax(ratio, 0)^buffer_n_exp_use), 1), 0)
    viability_curve <- data.frame(N = N_states, ploidy = ploidy_grid, endpoint_value = if (identical(start_with_mode, "chr_number"))
        as.numeric(N_states)
    else ploidy_grid, viability_after_ms = viability, row.names = NULL)
    ploidy_o2_curve <- dplyr::bind_rows(lapply(o2_levels_ploidy, function(o2_level) {
        N_grid <- if (identical(start_with_mode, "chr_number"))
            state_grid_dense
        else state_grid_dense * cfg$N_UNIT
        rate_df <- compute_rate_components(rep(o2_level, length(N_grid)), N_grid)
        data.frame(oxygen_pct = rep(as.numeric(o2_level), length(N_grid)), ploidy = if (identical(start_with_mode, "chr_number"))
            as.numeric(N_grid)
        else as.numeric(state_grid_dense), N = N_grid, endpoint_value = if (identical(start_with_mode, "chr_number"))
            as.numeric(N_grid)
        else as.numeric(state_grid_dense), proliferation_rate = rate_df$proliferation_rate, death_rate = rate_df$death_rate,
            net_growth_rate = rate_df$net_growth_rate, row.names = NULL)
    }))
    list(functional_curve_oxygen = o2_curve, functional_curve_oxygen_multi_ploidy = o2_curve_multi, functional_curve_ploidy = viability_curve,
        functional_curve_ploidy_by_o2 = ploidy_o2_curve)
}
