#!/usr/bin/env Rscript

# Simulation-specific compatibility helpers. Shared parsing, transformation,
# distance, and clustering functions are loaded from the canonical util module.
process_fingerprint_simulation_utils_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files)) {
    dirname(frame_files[[length(frame_files)]])
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})
process_fingerprint_simulation_workflow_root <- normalizePath(
  file.path(process_fingerprint_simulation_utils_dir, "..", ".."),
  mustWork = FALSE
)
source(
  file.path(
    process_fingerprint_simulation_workflow_root,
    "util",
    "o2_supply_demand_map_process_fingerprint_utils.R"
  ),
  local = TRUE
)
source(
  file.path(
    process_fingerprint_simulation_workflow_root,
    "util",
    "o2_supply_demand_map_postfit_input_utils.R"
  ),
  local = TRUE
)
rm(
  process_fingerprint_simulation_utils_dir,
  process_fingerprint_simulation_workflow_root
)


o2ipa_find_workflow_root <- function(path = o2ipa_script_dir()) {
  cur <- normalizePath(path, mustWork = FALSE)
  if (file.exists(cur) && !dir.exists(cur)) cur <- dirname(cur)
  for (i in seq_len(8L)) {
    if (
      file.exists(file.path(cur, "util", "o2_supply_demand_map_shared.R")) &&
        file.exists(file.path(cur, "model", "model_O2_supply_demand_MAP.R"))
    ) {
      return(normalizePath(cur, mustWork = FALSE))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(file.path(path, "..", ".."), mustWork = FALSE)
}











o2ipa_extract_all_params <- function(seed_id, seed_dir, summary_tab, matrix_tab) {
  best_path <- if (!is.na(seed_dir) && nzchar(seed_dir)) file.path(seed_dir, "best_params.tsv") else NA_character_
  best_vals <- if (!is.na(best_path) && file.exists(best_path)) o2ipa_read_best_params(best_path)$values else setNames(numeric(0), character(0))
  summary_row <- if (!is.null(summary_tab)) summary_tab[summary_tab$seed_id == seed_id, , drop = FALSE] else NULL
  if (!is.null(summary_row) && nrow(summary_row) > 1L) summary_row <- summary_row[1, , drop = FALSE]
  matrix_row <- if (!is.null(matrix_tab)) matrix_tab[matrix_tab$seed_id == seed_id, , drop = FALSE] else NULL
  if (!is.null(matrix_row) && nrow(matrix_row) > 1L) matrix_row <- matrix_row[1, , drop = FALSE]
  params <- o2ipa_target_params()
  rows <- lapply(params, function(p) {
    hit <- o2ipa_extract_param(seed_id, p, best_vals, summary_row, matrix_row)
    data.frame(
      seed_id = seed_id,
      parameter = p,
      value = hit$value,
      parameter_source = hit$source,
      matched_alias = hit$alias,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


o2ipa_collect_seed_inputs <- function(run_dir, objective_source = "auto") {
  manifest0 <- o2ipa_discover_seeds(run_dir)
  summary_tab <- o2ipa_read_extra_summary(run_dir)
  matrix_tab <- o2ipa_read_param_matrix(run_dir)
  boundary_path <- o2ipa_find_extra(run_dir, "parameter_boundary_long.tsv")
  boundary_long <- if (!is.na(boundary_path)) o2ipa_read_tsv(boundary_path) else NULL
  if (!is.null(boundary_long) && "seed" %in% names(boundary_long)) {
    boundary_long$seed_id <- o2ipa_norm_seed(boundary_long$seed)
  }

  param_rows <- list()
  manifest_rows <- vector("list", nrow(manifest0))
  for (i in seq_len(nrow(manifest0))) {
    seed_id <- manifest0$seed_id[[i]]
    seed_dir <- manifest0$seed_dir[[i]]
    if (is.na(seed_dir)) seed_dir <- ""
    fit_summary_path <- if (nzchar(seed_dir)) file.path(seed_dir, "fit_summary.tsv") else NA_character_
    summary_vals <- if (!is.na(fit_summary_path) && file.exists(fit_summary_path)) {
      o2ipa_metric_map(fit_summary_path)
    } else {
      row <- if (!is.null(summary_tab)) summary_tab[summary_tab$seed_id == seed_id, , drop = FALSE] else NULL
      if (!is.null(row) && nrow(row)) as.list(row[1, , drop = TRUE]) else list()
    }
    obj <- o2ipa_choose_objective(summary_vals, objective_source = objective_source)
    params_long <- o2ipa_extract_all_params(seed_id, seed_dir, summary_tab, matrix_tab)
    param_rows[[i]] <- params_long
    missing_params <- params_long$parameter[!is.finite(params_long$value)]

    bseed <- if (!is.null(boundary_long)) boundary_long[boundary_long$seed_id == seed_id, , drop = FALSE] else NULL
    target_boundary <- if (!is.null(bseed) && nrow(bseed)) {
      pname <- if ("param_prototype" %in% names(bseed)) bseed$param_prototype else if ("parameter" %in% names(bseed)) bseed$parameter else bseed$param_name
      bseed[pname %in% o2ipa_target_params(), , drop = FALSE]
    } else {
      data.frame()
    }
    boundary_risk <- FALSE
    n_near <- 0L
    if (nrow(target_boundary)) {
      status_col <- if ("bound_status" %in% names(target_boundary)) "bound_status" else if ("status" %in% names(target_boundary)) "status" else NA_character_
      if (!is.na(status_col)) {
        boundary_risk <- any(!is.na(target_boundary[[status_col]]) & !target_boundary[[status_col]] %in% c("interior", ""))
        n_near <- sum(!is.na(target_boundary[[status_col]]) & !target_boundary[[status_col]] %in% c("interior", ""))
      }
    }
    fit_success <- length(missing_params) == 0L && is.finite(obj$value)
    convergence <- o2ipa_chr_from_map(summary_vals, "deoptim_stop_reason")
    if (is.na(convergence)) convergence <- o2ipa_chr_from_map(summary_vals, "optimizer_local_convergence")
    failure_parts <- character(0)
    if (!is.finite(obj$value)) failure_parts <- c(failure_parts, "missing_objective")
    if (length(missing_params)) failure_parts <- c(failure_parts, paste0("missing_params:", paste(missing_params, collapse = ",")))
    manifest_rows[[i]] <- data.frame(
      seed_id = seed_id,
      seed_dir = if (nzchar(seed_dir)) normalizePath(seed_dir, mustWork = FALSE) else NA_character_,
      fit_success = fit_success,
      convergence_status = convergence,
      objective = obj$value,
      objective_source = obj$source,
      objective_total = o2ipa_num_from_map(summary_vals, "objective_total"),
      objective_data = o2ipa_num_from_map(summary_vals, "objective_data"),
      objective_burden = o2ipa_num_from_map(summary_vals, "objective_burden"),
      objective_ploidy = o2ipa_num_from_map(summary_vals, "objective_ploidy"),
      burden_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_burden_neg2loglik_raw"),
      ploidy_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_ploidy_neg2loglik_raw"),
      runtime = o2ipa_num_from_map(summary_vals, "runtime"),
      parameter_file = if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir, "best_params.tsv"))) file.path(seed_dir, "best_params.tsv") else NA_character_,
      config_file = if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir, "fit_config.rds"))) file.path(seed_dir, "fit_config.rds") else NA_character_,
      visualization_available = !is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir, "viz_status.log")),
      failure_reason = if (fit_success) NA_character_ else paste(failure_parts, collapse = ";"),
      boundary_risk = boundary_risk,
      number_of_parameters_near_boundary = n_near,
      stringsAsFactors = FALSE
    )
  }
  manifest <- do.call(rbind, manifest_rows)
  params_long <- do.call(rbind, param_rows)
  list(manifest = manifest, params_long = params_long, boundary_long = boundary_long, summary_tab = summary_tab, matrix_tab = matrix_tab)
}




o2ipa_crossing <- function(x, y, target, decreasing = FALSE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(list(value = NA_real_, status = "insufficient_grid"))
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  if (target < min(y) || target > max(y)) return(list(value = NA_real_, status = "no_crossing_in_grid"))
  if (isTRUE(decreasing)) {
    ord2 <- order(y)
    return(list(value = as.numeric(stats::approx(y[ord2], x[ord2], xout = target, ties = "ordered")$y), status = "crossed"))
  }
  list(value = as.numeric(stats::approx(y, x, xout = target, ties = "ordered")$y), status = "crossed")
}

o2ipa_feature_row <- function(seed_id, scope, module, feature, value, feature_type = "positive", source = "model_probe") {
  data.frame(
    seed_id = seed_id,
    fingerprint_scope = scope,
    module = module,
    feature = feature,
    raw_value = as.numeric(value),
    feature_type = feature_type,
    source = source,
    stringsAsFactors = FALSE
  )
}

o2ipa_feature_name_num <- function(prefix, value) {
  paste0(prefix, gsub("[^A-Za-z0-9]+", "p", format(value, trim = TRUE, scientific = FALSE)))
}

o2ipa_build_static_one <- function(seed_id, run_params, model_env, scope = "full") {
  rp <- as.list(run_params)
  names(rp) <- names(run_params)
  rp$ploidy_O2_death <- o2ipa_null_coalesce(rp$ploidy_O2_death, "ploidy_related")
  rp$O2_growth <- o2ipa_null_coalesce(rp$O2_growth, TRUE)
  rp$o2_min <- as.numeric(o2ipa_null_coalesce(rp$o2_min, 0))
  rp$o2_S0_upper_bound <- as.numeric(o2ipa_null_coalesce(rp$o2_S0_upper_bound, 5))
  rp$o2_Nref <- as.numeric(o2ipa_null_coalesce(rp$o2_Nref, 1e6))
  if (!is.finite(rp$o2_S0_upper_bound) || rp$o2_S0_upper_bound <= 0) rp$o2_S0_upper_bound <- 5
  o2_grid <- sort(unique(pmin(c(0.05, 0.1, 0.25, 0.5, 1, 2, 5), rp$o2_S0_upper_bound)))
  state_grid <- c(44, 66, 88, 132)
  state_grid <- state_grid[state_grid >= 22 & state_grid <= 154]
  p_mis_grid <- c(1e-4, 1e-3, 1e-2, 5e-2)
  demand_grid <- c(0, 1e3, 1e4, 1e5, 1e6, 1e7)
  burden_grid <- c(0.01, 0.1, 1, 10, 100)
  rows <- list()
  add <- function(module, feature, value, feature_type = "positive", source = "model_probe") {
    rows[[length(rows) + 1L]] <<- o2ipa_feature_row(seed_id, scope, module, feature, value, feature_type, source)
  }

  stress <- o2ipa_call_model(model_env, ".resource_stress_of_O2", O2 = o2_grid, run_params = rp)
  for (i in seq_along(o2_grid)) add("hypoxia_sensing", o2ipa_feature_name_num("stress_at_O2_", o2_grid[[i]]), stress[[i]], "fraction")
  add("hypoxia_sensing", "stress_AUC", o2ipa_auc(o2_grid, stress), "positive")
  dense_o2 <- seq(min(o2_grid), max(o2_grid), length.out = 400L)
  dense_stress <- o2ipa_call_model(model_env, ".resource_stress_of_O2", O2 = dense_o2, run_params = rp)
  for (thr in c(0.1, 0.5, 0.9)) {
    cr <- o2ipa_crossing(dense_o2, dense_stress, thr, decreasing = TRUE)
    add("hypoxia_sensing", paste0("O2_at_stress_", as.integer(thr * 100), "pct"), cr$value, "threshold")
    add("hypoxia_sensing", paste0("crossing_status_stress_", as.integer(thr * 100), "pct_crossed"), ifelse(cr$status == "crossed", 1, 0), "binary")
  }
  cr10 <- o2ipa_crossing(dense_o2, dense_stress, 0.1, decreasing = TRUE)
  cr90 <- o2ipa_crossing(dense_o2, dense_stress, 0.9, decreasing = TRUE)
  add("hypoxia_sensing", "stress_transition_width", if (is.finite(cr10$value) && is.finite(cr90$value)) abs(cr10$value - cr90$value) else NA_real_, "positive")

  grid <- expand.grid(O2 = o2_grid, N = state_grid)
  lam <- o2ipa_call_model(model_env, ".lambda_eff_of_O2", O2 = grid$O2, run_params = rp, N = grid$N, O2_growth = isTRUE(rp$O2_growth))
  mu <- o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = grid$O2, run_params = rp, N = grid$N)
  pm <- o2ipa_call_model(model_env, ".pmisseg_of_O2", O2 = grid$O2, run_params = rp, N = grid$N)
  net <- lam - mu
  for (i in seq_len(nrow(grid))) {
    suffix <- paste0("_O2_", gsub("[^0-9]+", "p", format(grid$O2[[i]], trim = TRUE, scientific = FALSE)), "_N_", grid$N[[i]])
    add("proliferation", paste0("lambda", suffix), lam[[i]], "positive")
    add("death", paste0("mu", suffix), mu[[i]], "positive")
    add("net_growth", paste0("net_growth", suffix), net[[i]], "signed")
    add("CIN_missegregation", paste0("p_misseg", suffix), pm[[i]], "probability")
  }
  for (n in state_grid) {
    idx <- grid$N == n
    add("proliferation", paste0("proliferation_AUC_over_O2_N_", n), o2ipa_auc(grid$O2[idx], lam[idx]), "positive")
    add("death", paste0("death_AUC_over_O2_N_", n), o2ipa_auc(grid$O2[idx], mu[idx]), "positive")
    add("net_growth", paste0("positive_net_growth_AUC_N_", n), o2ipa_auc(grid$O2[idx], pmax(net[idx], 0)), "positive")
    add("net_growth", paste0("negative_net_growth_AUC_N_", n), o2ipa_auc(grid$O2[idx], pmax(-net[idx], 0)), "positive")
    cr <- o2ipa_crossing(grid$O2[idx], net[idx], 0, decreasing = FALSE)
    add("net_growth", paste0("critical_O2_net_growth_zero_N_", n), cr$value, "threshold")
  }
  dip_hi <- lam[grid$N == 44 & grid$O2 == max(o2_grid)][[1]]
  dip_lo <- lam[grid$N == 44 & grid$O2 == min(o2_grid)][[1]]
  add("proliferation", "low_O2_growth_suppression_ratio_diploid", if (is.finite(dip_hi) && dip_hi != 0) dip_lo / dip_hi else NA_real_, "ratio")
  high_hi <- lam[grid$N == max(state_grid) & grid$O2 == max(o2_grid)][[1]]
  add("proliferation", "high_ploidy_vs_diploid_growth_ratio_high_O2", if (is.finite(dip_hi) && dip_hi != 0) high_hi / dip_hi else NA_real_, "ratio")
  mu_hi <- mu[grid$N == 44 & grid$O2 == max(o2_grid)][[1]]
  mu_lo <- mu[grid$N == 44 & grid$O2 == min(o2_grid)][[1]]
  add("death", "low_O2_death_induction_ratio_diploid", if (is.finite(mu_hi) && mu_hi != 0) mu_lo / mu_hi else NA_real_, "ratio")
  high_mu <- mu[grid$N == max(state_grid) & grid$O2 == max(o2_grid)][[1]]
  add("death", "high_ploidy_vs_diploid_death_ratio_high_O2", if (is.finite(mu_hi) && mu_hi != 0) high_mu / mu_hi else NA_real_, "ratio")
  pm_hi <- pm[grid$N == 44 & grid$O2 == max(o2_grid)][[1]]
  pm_lo <- pm[grid$N == 44 & grid$O2 == min(o2_grid)][[1]]
  add("CIN_missegregation", "hypoxia_induced_misseg_increment_diploid", pm_lo - pm_hi, "probability")
  add("CIN_missegregation", "low_O2_to_high_O2_p_misseg_ratio_diploid", if (is.finite(pm_hi) && pm_hi != 0) pm_lo / pm_hi else NA_real_, "ratio")
  add("CIN_missegregation", "expected_missegregation_events_per_100_divisions_low_O2_diploid", 100 * pm_lo, "positive")

  for (n in state_grid) {
    for (p in p_mis_grid) {
      delta <- o2ipa_call_model(model_env, ".pr_delta_vec",
        N = n, p = p, eps_tail = 1e-8, N_unit = 22L,
        buffer_smax = as.numeric(rp$buffer_smax),
        buffer_beta = as.numeric(rp$buffer_beta),
        buffer_n_exp = as.numeric(rp$buffer_n_exp)
      )
      shifts <- suppressWarnings(as.numeric(names(delta)))
      prob <- as.numeric(delta)
      md <- as.numeric(attr(delta, "mass_dropped"))
      suffix <- paste0("_N_", n, "_p_", gsub("[^0-9]+", "p", format(p, trim = TRUE, scientific = FALSE)))
      add("aneuploidy_buffering", paste0("surviving_transition_mass", suffix), sum(prob, na.rm = TRUE), "positive")
      add("aneuploidy_buffering", paste0("mass_dropped", suffix), md, "fraction")
      add("aneuploidy_buffering", paste0("mean_absolute_chromosome_shift_after_buffering", suffix), sum(abs(shifts) * prob, na.rm = TRUE), "positive")
      mu_shift <- sum(shifts * prob, na.rm = TRUE)
      add("aneuploidy_buffering", paste0("transition_variance", suffix), sum(((shifts - mu_shift)^2) * prob, na.rm = TRUE), "positive")
      add("aneuploidy_buffering", paste0("large_shift_survival_fraction", suffix), sum(prob[abs(shifts) >= 3], na.rm = TRUE), "fraction")
    }
  }

  pw <- as.numeric(rp$p_wgd)
  add("WGD", "p_wgd_per_division", pw, "probability")
  add("WGD", "expected_WGD_events_per_100_divisions", 100 * pw, "positive")
  add("WGD", "probability_at_least_one_WGD_in_100_divisions", 1 - (1 - min(max(pw, 0), 1))^100, "probability")

  o2_eff <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = demand_grid, run_params = rp, o2_Nref = as.numeric(rp$o2_Nref))
  for (i in seq_along(demand_grid)) add("O2_supply_demand", o2ipa_feature_name_num("O2_at_effective_demand_", demand_grid[[i]]), o2_eff[[i]], "oxygen")
  for (n in state_grid) {
    demand <- burden_grid * as.numeric(rp$rho_2N) * ((n / 44)^as.numeric(rp$eta_o2))
    o2_obs <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = demand, run_params = rp, o2_Nref = as.numeric(rp$o2_Nref))
    for (i in seq_along(burden_grid)) {
      add("O2_supply_demand", paste0(o2ipa_feature_name_num("O2_at_observed_burden_", burden_grid[[i]]), "_N_", n), o2_obs[[i]], "oxygen")
      add("O2_supply_demand", paste0(o2ipa_feature_name_num("effective_demand_at_observed_burden_", burden_grid[[i]]), "_N_", n), demand[[i]], "positive")
    }
    for (target_o2 in c(2, 1, 0.5)) {
      dense_b <- 10^seq(log10(min(burden_grid)), log10(max(burden_grid)), length.out = 200L)
      dense_d <- dense_b * as.numeric(rp$rho_2N) * ((n / 44)^as.numeric(rp$eta_o2))
      dense_o <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = dense_d, run_params = rp, o2_Nref = as.numeric(rp$o2_Nref))
      cr <- o2ipa_crossing(log10(dense_b), dense_o, target_o2, decreasing = TRUE)
      add("O2_supply_demand", paste0("log_burden_at_O2_", gsub("[^0-9]+", "p", target_o2), "pct_N_", n), cr$value, "signed")
    }
  }
  add("O2_supply_demand", "minimum_O2_floor", min(o2_eff, na.rm = TRUE), "oxygen")
  add("O2_supply_demand", "oxygen_depletion_AUC_effective_demand", o2ipa_auc(log10(demand_grid + 1), pmax(max(o2_eff, na.rm = TRUE) - o2_eff, 0)), "positive")
  do.call(rbind, rows)
}


o2ipa_dynamic_fingerprints <- function(run_dir, seed_ids) {
  burden_path <- o2ipa_find_extra(run_dir, "predict1000_burden_total_seed_day_mean.tsv")
  ploidy_path <- o2ipa_find_extra(run_dir, "predict1000_ploidy_seed_day_mean.tsv")
  rows <- list()
  unavailable <- c("O2_min", "O2_terminal", "O2_AUC", "time_below_O2_2pct", "time_below_O2_1pct", "time_below_O2_0.5pct", "terminal_high_ploidy_fraction", "maximum_high_ploidy_fraction", "ploidy_entropy_AUC", "terminal_ploidy_entropy", "WGD_or_high_ploidy_lineage_fraction")
  if (!is.na(burden_path)) {
    tab <- o2ipa_read_tsv(burden_path)
    tab$seed_id <- o2ipa_norm_seed(tab$seed)
    tab <- tab[tab$seed_id %in% seed_ids, , drop = FALSE]
    for (key in split(seq_len(nrow(tab)), paste(tab$seed_id, tab$cohort, sep = "||"))) {
      d <- tab[key, , drop = FALSE]
      d <- d[order(d$day), , drop = FALSE]
      seed <- d$seed_id[[1]]
      cohort <- d$cohort[[1]]
      y <- as.numeric(d$burden_value)
      x <- as.numeric(d$day)
      prefix <- paste0("cohort_", cohort, "_")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "terminal_burden"), tail(y, 1), "positive", "existing_prediction")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "burden_AUC"), o2ipa_auc(x, y), "positive", "existing_prediction")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "maximum_burden"), max(y, na.rm = TRUE), "positive", "existing_prediction")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "time_to_burden_double"), o2ipa_crossing(x, y, 2 * y[[1]])$value, "time", "existing_prediction")
    }
  }
  if (!is.na(ploidy_path)) {
    tab <- o2ipa_read_tsv(ploidy_path)
    tab$seed_id <- o2ipa_norm_seed(tab$seed)
    tab <- tab[tab$seed_id %in% seed_ids, , drop = FALSE]
    for (key in split(seq_len(nrow(tab)), paste(tab$seed_id, tab$cohort, sep = "||"))) {
      d <- tab[key, , drop = FALSE]
      d <- d[order(d$day), , drop = FALSE]
      seed <- d$seed_id[[1]]
      cohort <- d$cohort[[1]]
      y <- as.numeric(d$ploidy_value)
      x <- as.numeric(d$day)
      prefix <- paste0("cohort_", cohort, "_")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_ploidy", paste0(prefix, "terminal_mean_chromosome"), tail(y, 1), "positive", "existing_prediction")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_ploidy", paste0(prefix, "maximum_mean_chromosome"), max(y, na.rm = TRUE), "positive", "existing_prediction")
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(long = out, unavailable = data.frame(feature = unavailable, reason = "not_available_in_existing_prediction_outputs", stringsAsFactors = FALSE))
}

o2ipa_make_mock_run <- function(run_dir, n_seed = 5L) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  params <- o2ipa_target_params()
  base <- c(
    O2_crit = 0.2, alpha_o2 = 1.2, gamma_growth = 3, mu_hp = 0.05, p_misseg = 0.1,
    k_o_mis = 0.003, buffer_smax = 0.9, buffer_beta = 0.2, buffer_n_exp = 4,
    n_O = 1.5, lam_max = 0.22, p_mis_base = 1e-4, p_wgd = 1e-4, gamma_mu = 2.5,
    o2_S0 = 3.5, kappa_O = 0.8, eta_o2 = 1.1, rho_2N = 1e5
  )
  set.seed(123)
  for (i in seq_len(n_seed)) {
    sdir <- file.path(run_dir, paste0("seed", i))
    dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
    mult <- exp(stats::rnorm(length(params), 0, 0.12))
    vals <- base[params] * mult
    vals[c("gamma_growth", "gamma_mu", "n_O", "buffer_smax")] <- base[c("gamma_growth", "gamma_mu", "n_O", "buffer_smax")] + stats::rnorm(4, 0, 0.05)
    vals["buffer_smax"] <- min(max(vals["buffer_smax"], 0.05), 0.99)
    best <- data.frame(parameter = names(vals), value = as.numeric(vals), stringsAsFactors = FALSE)
    best <- rbind(best, data.frame(parameter = c("o2_min", "o2_S0_upper_bound", "o2_Nref"), value = c(0, 5, 1e6)))
    o2ipa_write_tsv(best, file.path(sdir, "best_params.tsv"))
    fs <- data.frame(metric = c("objective", "objective_data", "objective_burden", "objective_ploidy", "objective_burden_neg2loglik_raw", "objective_ploidy_neg2loglik_raw", "deoptim_stop_reason"), value = c(10 + i / 10, 9 + i / 10, 4 + i / 20, 5 + i / 20, 1000 + i, 2000 + i, "mock_converged"))
    o2ipa_write_tsv(fs, file.path(sdir, "fit_summary.tsv"))
    pt <- data.frame(param_name = paste0(ifelse(params %in% c("gamma_growth", "gamma_mu", "n_O", "buffer_smax"), "", "log10_"), params), estimate = TRUE, init_value = 0, lower_bound = -10, upper_bound = 10, param_prototype = params, prototype_init_value = base[params], prototype_lower_bound = base[params] / 10, prototype_upper_bound = base[params] * 10, source = "mock", note = "mock", stringsAsFactors = FALSE)
    utils::write.csv(pt, file.path(sdir, "parameter_table.csv"), row.names = FALSE)
  }
  extra <- file.path(run_dir, "extra_results")
  dir.create(extra, recursive = TRUE, showWarnings = FALSE)
  days <- 0:20
  burden <- do.call(rbind, lapply(seq_len(n_seed), function(i) {
    do.call(rbind, lapply(c("2N", "4N"), function(cohort) data.frame(seed = paste0("seed", i), cohort = cohort, day = days, burden_value = (10 + i) * exp(0.05 * days) * ifelse(cohort == "4N", 1.1, 1), stringsAsFactors = FALSE)))
  }))
  ploidy <- do.call(rbind, lapply(seq_len(n_seed), function(i) {
    do.call(rbind, lapply(c("2N", "4N"), function(cohort) data.frame(seed = paste0("seed", i), cohort = cohort, day = days, ploidy_value = ifelse(cohort == "4N", 88, 44) + i * 0.2 + sin(days / 5), stringsAsFactors = FALSE)))
  }))
  o2ipa_write_tsv(burden, file.path(extra, "predict1000_burden_total_seed_day_mean.tsv"))
  o2ipa_write_tsv(ploidy, file.path(extra, "predict1000_ploidy_seed_day_mean.tsv"))
  invisible(run_dir)
}
