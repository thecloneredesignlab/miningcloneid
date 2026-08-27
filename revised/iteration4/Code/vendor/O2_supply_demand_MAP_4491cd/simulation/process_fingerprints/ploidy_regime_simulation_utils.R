#!/usr/bin/env Rscript

# Simulation-specific helpers. Shared ploidy-regime functions are loaded from
# the canonical util module so analysis and simulation expose identical APIs.

ploidy_regime_simulation_utils_dir <- local({
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
ploidy_regime_simulation_workflow_root <- normalizePath(
  file.path(ploidy_regime_simulation_utils_dir, "..", ".."),
  mustWork = FALSE
)
source(
  file.path(
    ploidy_regime_simulation_workflow_root,
    "util",
    "o2_supply_demand_map_ploidy_regime_utils.R"
  ),
  local = TRUE
)
source(
  file.path(
    ploidy_regime_simulation_workflow_root,
    "util",
    "o2_supply_demand_map_postfit_probe_utils.R"
  ),
  local = TRUE
)
rm(
  ploidy_regime_simulation_utils_dir,
  ploidy_regime_simulation_workflow_root
)

o2pr_source_process_utils <- function(script_dir) {
  source(file.path(script_dir, "process_fingerprint_simulation_legacy_utils.R"), local = TRUE)
}

o2pr_git_value <- function(args, repo_root) {
  out <- tryCatch(
    system2("git", c("-C", repo_root, args), stdout = TRUE, stderr = TRUE),
    error = function(e) NA_character_
  )
  if (!length(out)) NA_character_ else paste(out, collapse = "\n")
}



o2pr_generate_viz_if_requested <- function(seed_dir, script_dir, generate_viz = FALSE) {
  if (!isTRUE(generate_viz)) return(invisible(FALSE))
  stop(
    "Analysis no longer generates visualization or simulation outputs. ",
    "Run simulation/invivo/generate_invivo_simulation_outputs.R for ",
    seed_dir,
    " before process-fingerprint analysis.",
    call. = FALSE
  )
}

o2pr_find_trajectory_path <- function(seed_dir) {
  candidates <- file.path(seed_dir, "simulation", "invivo", c(
    "predict_ploidy_weighted_mean_0_1000day.tsv",
    "ploidy_weighted_mean_timecourse.tsv",
    "predict_ploidy_weighted_mean_0_300day.tsv",
    "predict_ploidy_0_1000day.tsv",
    "ploidy_timecourse.tsv",
    "predict_ploidy_0_300day.tsv"
  ))
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1]] else NA_character_
}

o2pr_value_col_for_trajectory <- function(tab) {
  for (nm in c("weighted_mean_N", "weighted_mean_endpoint", "ploidy_value", "weighted_mean_ploidy")) {
    if (nm %in% names(tab)) return(nm)
  }
  NA_character_
}

o2pr_weighted_from_density <- function(tab) {
  needed <- c("cohort", "day", "N", "fraction")
  if (!all(needed %in% names(tab))) return(NULL)
  tab$day <- suppressWarnings(as.numeric(tab$day))
  tab$N <- suppressWarnings(as.numeric(tab$N))
  tab$fraction <- suppressWarnings(as.numeric(tab$fraction))
  tab <- tab[is.finite(tab$day) & is.finite(tab$N) & is.finite(tab$fraction), , drop = FALSE]
  if (!nrow(tab)) return(NULL)
  key <- paste(tab$cohort, tab$day, sep = "\r")
  rows <- lapply(split(seq_len(nrow(tab)), key), function(idx) {
    d <- tab[idx, , drop = FALSE]
    denom <- sum(d$fraction, na.rm = TRUE)
    data.frame(
      cohort = as.character(d$cohort[[1]]),
      day = as.numeric(d$day[[1]]),
      value = if (is.finite(denom) && denom > 0) sum(d$N * d$fraction, na.rm = TRUE) / denom else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

o2pr_read_trajectory_seed <- function(seed_id, seed_dir, cfg, script_dir, generate_viz = FALSE) {
  if (is.na(seed_dir) || !nzchar(seed_dir) || !dir.exists(seed_dir)) {
    return(list(curve = data.frame(), status = "missing_seed_dir", source_file = NA_character_))
  }
  path <- o2pr_find_trajectory_path(seed_dir)
  if (is.na(path)) {
    o2pr_generate_viz_if_requested(seed_dir, script_dir, generate_viz)
    path <- o2pr_find_trajectory_path(seed_dir)
  }
  if (is.na(path)) {
    return(list(curve = data.frame(), status = "missing_trajectory_table", source_file = NA_character_))
  }
  tab <- tryCatch(o2pr_read_tsv(path), error = function(e) NULL)
  if (is.null(tab) || !nrow(tab) || !"cohort" %in% names(tab) || !"day" %in% names(tab)) {
    return(list(curve = data.frame(), status = "bad_trajectory_table", source_file = path))
  }
  vcol <- o2pr_value_col_for_trajectory(tab)
  if (is.na(vcol)) {
    dens <- o2pr_weighted_from_density(tab)
    if (is.null(dens)) return(list(curve = data.frame(), status = "no_weighted_value_column", source_file = path))
    tab2 <- dens
  } else {
    tab2 <- data.frame(
      cohort = as.character(tab$cohort),
      day = suppressWarnings(as.numeric(tab$day)),
      value = suppressWarnings(as.numeric(tab[[vcol]])),
      stringsAsFactors = FALSE
    )
    if (identical(vcol, "weighted_mean_ploidy") && !identical(as.character(cfg$start_with %||% "ploidy"), "chr_number")) {
      tab2$value <- tab2$value * as.numeric(cfg$N_UNIT %||% 22)
    }
  }
  tab2 <- tab2[tab2$cohort %in% c("2N", "4N") & is.finite(tab2$day) & is.finite(tab2$value), , drop = FALSE]
  if (!nrow(tab2)) return(list(curve = data.frame(), status = "no_valid_2N_4N_rows", source_file = path))
  agg <- aggregate(tab2$value, by = list(cohort = tab2$cohort, day = tab2$day), FUN = mean, na.rm = TRUE)
  names(agg)[3] <- "trajectory_value"
  agg$seed_id <- seed_id
  agg$source_file <- normalizePath(path, mustWork = FALSE)
  agg <- agg[, c("seed_id", "cohort", "day", "trajectory_value", "source_file")]
  list(curve = agg[order(agg$cohort, agg$day), , drop = FALSE], status = "ok", source_file = path)
}

o2pr_read_density_seed <- function(seed_id, seed_dir) {
  if (is.na(seed_dir) || !nzchar(seed_dir)) return(data.frame())
  path <- file.path(
    seed_dir,
    "simulation",
    "invivo",
    "ploidy_timecourse.tsv"
  )
  if (!file.exists(path)) {
    path <- file.path(
      seed_dir,
      "simulation",
      "invivo",
      "predict_ploidy_0_1000day.tsv"
    )
  }
  if (!file.exists(path)) return(data.frame())
  tab <- tryCatch(o2pr_read_tsv(path), error = function(e) data.frame())
  needed <- c("cohort", "day", "N", "fraction")
  if (!all(needed %in% names(tab))) return(data.frame())
  tab$seed_id <- seed_id
  tab[, c("seed_id", needed), drop = FALSE]
}

o2pr_collect_trajectories <- function(inputs, script_dir, trajectory_end_day, mid_window, late_window, generate_viz = FALSE, cfg = NULL) {
  cfg <- cfg %||% o2pr_first_seed_cfg(inputs$manifest)
  curves <- list()
  status <- list()
  density <- list()
  manifest <- inputs$manifest
  for (i in seq_len(nrow(manifest))) {
    seed <- manifest$seed_id[[i]]
    seed_dir <- manifest$seed_dir[[i]]
    tr <- o2pr_read_trajectory_seed(seed, seed_dir, cfg, script_dir, generate_viz)
    if (nrow(tr$curve)) curves[[seed]] <- tr$curve
    density[[seed]] <- o2pr_read_density_seed(seed, seed_dir)
    status[[seed]] <- data.frame(
      seed_id = seed,
      trajectory_available = identical(tr$status, "ok"),
      trajectory_status = tr$status,
      trajectory_source_file = tr$source_file,
      stringsAsFactors = FALSE
    )
  }
  curve_all <- if (length(curves)) do.call(rbind, curves) else data.frame()
  status_all <- do.call(rbind, status)
  density_all <- do.call(rbind, density)
  nmin <- as.numeric(cfg$N_MIN %||% 22)
  long_features <- list()
  if (nrow(curve_all)) {
    for (key in split(seq_len(nrow(curve_all)), paste(curve_all$seed_id, curve_all$cohort, sep = "||"))) {
      d <- curve_all[key, , drop = FALSE]
      if (!nrow(d)) next
      long_features[[length(long_features) + 1L]] <- o2pr_cohort_features(
        d$seed_id[[1]], d$cohort[[1]], d, mid_window, late_window, nmin
      )
    }
  }
  feat_long <- if (length(long_features)) do.call(rbind, long_features) else data.frame()
  labels <- o2pr_make_trajectory_labels(curve_all, feat_long, status_all, cfg, trajectory_end_day, mid_window, late_window)
  wide <- o2pr_features_wide(feat_long, labels)
  list(curves = curve_all, density = density_all, status = status_all, features_long = feat_long, features_wide = wide, labels = labels)
}



o2pr_build_process_one <- function(seed_id, run_params, model_env, cfg, o2_grid, scope = "full") {
  rows <- list()
  add <- function(module, feature, value, type = "positive", source = "model_probe") {
    rows[[length(rows) + 1L]] <<- o2ipa_feature_row(seed_id, scope, module, feature, value, type, source)
  }
  states <- as.integer(pmin(pmax(c(22, 44, 66, 88), cfg$N_MIN %||% 22), cfg$N_MAX %||% 154))
  states <- sort(unique(states))
  state_note <- data.frame(requested_state = c(22, 44, 66, 88), used_state = as.integer(pmin(pmax(c(22, 44, 66, 88), cfg$N_MIN %||% 22), cfg$N_MAX %||% 154)))
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22), as.integer(cfg$N_MAX %||% 154))
  for (O2 in o2_grid) {
    lam <- o2ipa_call_model(model_env, ".lambda_eff_of_O2", O2 = rep(O2, length(states)), run_params = run_params, N = states, O2_growth = isTRUE(run_params$O2_growth %||% TRUE))
    mu <- o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(states)), run_params = run_params, N = states)
    pm <- o2ipa_call_model(model_env, ".pmisseg_of_O2", O2 = rep(O2, length(states)), run_params = run_params, N = states)
    stress <- o2ipa_call_model(model_env, ".resource_stress_of_O2", O2 = O2, run_params = run_params)
    for (i in seq_along(states)) {
      suffix <- paste0("_O2_", gsub("[^0-9]+", "p", O2), "_N_", states[[i]])
      add("hypoxia_sensing", paste0("hypoxia_stress", suffix), stress, "fraction")
      add("growth_selection", paste0("lambda", suffix), lam[[i]], "positive")
      add("death_selection", paste0("mu", suffix), mu[[i]], "positive")
      add("CIN", paste0("p_misseg", suffix), pm[[i]], "probability")
      add("CIN", paste0("missegregation_flux", suffix), lam[[i]] * pm[[i]], "positive")
      add("WGD", paste0("WGD_flux", suffix), lam[[i]] * as.numeric(run_params$p_wgd %||% 0), "positive")
    }
    G <- o2pr_build_G(model_env, cfg, run_params, O2)
    mu_all <- as.numeric(o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(ngrid)), run_params = run_params, N = ngrid))
    M <- G - Matrix::Diagonal(x = mu_all)
    eff <- vapply(states, function(N) {
      col <- as.integer(N - (cfg$N_MIN %||% 22) + 1L)
      if (col < 1L || col > ncol(G)) return(NA_real_)
      sum(G[, col]) - mu_all[[col]]
    }, numeric(1))
    names(eff) <- states
    if (all(c("22", "44") %in% names(eff))) add("attractor", paste0("selection_22_vs_44_O2_", gsub("[^0-9]+", "p", O2)), eff[["22"]] - eff[["44"]], "signed")
    if (all(c("44", "88") %in% names(eff))) add("attractor", paste0("selection_44_vs_88_O2_", gsub("[^0-9]+", "p", O2)), eff[["44"]] - eff[["88"]], "signed")
    eig <- tryCatch(eigen(as.matrix(M), only.values = FALSE), error = function(e) NULL)
    if (!is.null(eig)) {
      idx <- which.max(Re(eig$values))
      v <- Re(eig$vectors[, idx])
      if (sum(v, na.rm = TRUE) < 0) v <- -v
      nonneg <- all(v >= -1e-8, na.rm = TRUE)
      v <- pmax(v, 0)
      if (sum(v) > 0) v <- v / sum(v)
      lambda1 <- Re(eig$values[[idx]])
      lambda2 <- sort(Re(eig$values), decreasing = TRUE)[min(2L, length(eig$values))]
      add("attractor", paste0("dominant_mean_N_O2_", gsub("[^0-9]+", "p", O2)), sum(ngrid * v), "positive", "eigen")
      add("attractor", paste0("dominant_fraction_N_le_25_O2_", gsub("[^0-9]+", "p", O2)), sum(v[ngrid <= 25]), "fraction", "eigen")
      add("attractor", paste0("dominant_fraction_N_below_44_O2_", gsub("[^0-9]+", "p", O2)), sum(v[ngrid < 44]), "fraction", "eigen")
      add("attractor", paste0("dominant_growth_rate_O2_", gsub("[^0-9]+", "p", O2)), lambda1, "signed", "eigen")
      add("attractor", paste0("spectral_gap_O2_", gsub("[^0-9]+", "p", O2)), lambda1 - lambda2, "signed", "eigen")
      add("attractor", paste0("eigenvector_nonnegative_O2_", gsub("[^0-9]+", "p", O2)), as.numeric(nonneg), "binary", "eigen")
    }
    tri <- attr(G, "triplet")
    if (!is.null(tri$boundary_dropped_rate)) {
      for (N in states) {
        col <- as.integer(N - (cfg$N_MIN %||% 22) + 1L)
        add("buffering", paste0("transition_mass_dropped_at_N_", N, "_O2_", gsub("[^0-9]+", "p", O2)), tri$boundary_dropped_rate[[col]], "positive")
      }
    }
  }
  for (N in states) {
    pN <- as.numeric(o2ipa_call_model(model_env, ".pmisseg_of_O2", O2 = min(o2_grid), run_params = run_params, N = N))
    delta <- o2ipa_call_model(model_env, ".pr_delta_vec",
      N = N, p = pN, eps_tail = 1e-8, N_unit = as.integer(cfg$N_UNIT %||% 22),
      buffer_smax = as.numeric(run_params$buffer_smax %||% 1),
      buffer_beta = as.numeric(run_params$buffer_beta %||% 0),
      buffer_n_exp = as.numeric(run_params$buffer_n_exp %||% 1)
    )
    shifts <- suppressWarnings(as.numeric(names(delta)))
    prob <- as.numeric(delta)
    add("buffering", paste0("surviving_transition_mass_N_", N), sum(prob, na.rm = TRUE), "positive")
    add("buffering", paste0("mass_dropped_N_", N), as.numeric(attr(delta, "mass_dropped") %||% NA), "fraction")
    add("buffering", paste0("expected_delta_N_N_", N), sum(shifts * prob, na.rm = TRUE), "signed")
    add("buffering", paste0("expected_absolute_delta_N_N_", N), sum(abs(shifts) * prob, na.rm = TRUE), "positive")
    add("buffering", paste0("probability_surviving_below_source_N_", N), sum(prob[shifts < 0], na.rm = TRUE), "fraction")
    add("buffering", paste0("probability_surviving_above_source_N_", N), sum(prob[shifts > 0], na.rm = TRUE), "fraction")
    add("buffering", paste0("probability_entering_N_below_44_from_N_", N), sum(prob[N + shifts < 44], na.rm = TRUE), "fraction")
    add("buffering", paste0("probability_entering_N_le_25_from_N_", N), sum(prob[N + shifts <= 25], na.rm = TRUE), "fraction")
  }
  out <- do.call(rbind, rows)
  attr(out, "state_note") <- state_note
  out
}

o2pr_o2_exposure_one <- function(seed_id, run_params, model_env, cfg, burden_grid, realized_o2 = NULL, switch_o2 = NA_real_, scope = "full") {
  rows <- list()
  add <- function(module, feature, value, type = "positive", source = "model_probe") {
    rows[[length(rows) + 1L]] <<- o2ipa_feature_row(seed_id, scope, module, feature, value, type, source)
  }
  demand_grid <- burden_grid * as.numeric(run_params$rho_2N %||% 1) * ((44 / 44)^as.numeric(run_params$eta_o2 %||% 1))
  o2 <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = demand_grid, run_params = run_params, o2_Nref = as.numeric(cfg$o2_Nref %||% 1e6))
  for (i in seq_along(burden_grid)) {
    suffix <- paste0("_burden_", gsub("[^0-9]+", "p", burden_grid[[i]]))
    add("oxygen_feedback", paste0("O2_at_common", suffix), o2[[i]], "oxygen")
    add("oxygen_feedback", paste0("effective_demand_at_common", suffix), demand_grid[[i]], "positive")
  }
  for (target in c(2, 1, 0.5)) {
    dense_b <- 10^seq(log10(min(burden_grid)), log10(max(burden_grid)), length.out = 300L)
    dense_d <- dense_b * as.numeric(run_params$rho_2N %||% 1) * ((44 / 44)^as.numeric(run_params$eta_o2 %||% 1))
    dense_o <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = dense_d, run_params = run_params, o2_Nref = as.numeric(cfg$o2_Nref %||% 1e6))
    cr <- o2ipa_crossing(log10(dense_b), dense_o, target, decreasing = TRUE)
    target_key <- gsub("[^0-9]+", "p", target)
    add("oxygen_feedback", paste0("burden_reaches_O2_", target_key, "pct"), as.numeric(is.finite(cr$value)), "binary")
    add("oxygen_feedback", paste0("log10_burden_at_O2_", target_key, "pct"), cr$value, "signed")
  }
  add("oxygen_feedback", "dO2_dlogBurden_median", stats::median(diff(o2) / diff(log10(demand_grid + 1)), na.rm = TRUE), "signed")
  add("oxygen_feedback", "minimum_O2_floor", min(o2, na.rm = TRUE), "oxygen")
  if (!is.null(realized_o2) && nrow(realized_o2)) {
    r <- realized_o2[realized_o2$seed_id == seed_id, , drop = FALSE]
    if (nrow(r)) {
      add("oxygen_feedback", "minimum_realized_O2", min(r$pred_o2_pct, na.rm = TRUE), "oxygen", "realized_prediction")
      add("oxygen_feedback", "terminal_O2", tail(r$pred_o2_pct[order(r$day)], 1), "oxygen", "realized_prediction")
      add("oxygen_feedback", "O2_AUC", o2pr_auc(r$day, r$pred_o2_pct), "positive", "realized_prediction")
      for (target in c(2, 1, 0.5)) {
        rr <- r[order(r$day), , drop = FALSE]
        hit <- rr$day[rr$pred_o2_pct < target]
        first_hit <- if (length(hit)) hit[[1]] else NA_real_
        target_key <- gsub("[^0-9]+", "p", target)
        add("oxygen_feedback", paste0("ever_below_O2_", target_key, "pct"), as.numeric(is.finite(first_hit)), "binary", "realized_prediction")
        add("oxygen_feedback", paste0("time_below_O2_", target_key, "pct"), first_hit, "time", "realized_prediction")
      }
      min_realized_o2 <- min(r$pred_o2_pct, na.rm = TRUE)
      switch_available <- is.finite(switch_o2)
      add("oxygen_feedback", "switch_available", as.numeric(switch_available), "binary", "switch_diagnostic")
      add("oxygen_feedback", "switch_margin", if (switch_available && is.finite(min_realized_o2)) min_realized_o2 - switch_o2 else NA_real_, "signed", "switch_diagnostic")
      add("oxygen_feedback", "switch_crossed", as.numeric(switch_available && is.finite(min_realized_o2) && min_realized_o2 <= switch_o2), "binary", "switch_diagnostic")
    }
  }
  do.call(rbind, rows)
}

o2pr_read_realized_o2 <- function(inputs) {
  rows <- list()
  for (i in seq_len(nrow(inputs$manifest))) {
    seed <- inputs$manifest$seed_id[[i]]
    d <- inputs$manifest$seed_dir[[i]]
    p <- file.path(
      d,
      "simulation",
      "invivo",
      "predict_burden_0_1000day.tsv"
    )
    if (!file.exists(p)) {
      p <- file.path(
        d,
        "simulation",
        "invivo",
        "predict_burden_vs_o2.tsv"
      )
    }
    if (!file.exists(p)) next
    tab <- tryCatch(o2pr_read_tsv(p), error = function(e) data.frame())
    if (!all(c("day", "pred_o2_pct") %in% names(tab))) next
    tab$seed_id <- seed
    rows[[seed]] <- tab[, intersect(c("seed_id", "cohort", "day", "pred_o2_pct", "pred_o2_target_pct"), names(tab)), drop = FALSE]
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

o2pr_build_process_outputs <- function(inputs, out_dir, model_env, cfg, o2_grid, force = FALSE) {
  manifest <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  best_seed <- manifest$seed_id[which.min(manifest$objective)]
  ref <- as.numeric(param_mat[best_seed, , drop = TRUE])
  names(ref) <- colnames(param_mat)
  target <- o2ipa_target_params()
  real_o2 <- o2pr_read_realized_o2(inputs)
  burden_grid <- c(0.01, 0.1, 1, 10, 100)
  intr_full <- list()
  intr_18 <- list()
  exp_full <- list()
  exp_18 <- list()
  state_note <- NULL
  for (seed in manifest$seed_id) {
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    full_rp <- o2pr_run_params_from_vec(pvec, cfg)
    only_rp <- ref
    only_rp[target] <- pvec[target]
    only_rp <- o2pr_run_params_from_vec(only_rp, cfg)
    full_intr <- o2pr_build_process_one(seed, full_rp, model_env, cfg, o2_grid, "intrinsic_full")
    only_intr <- o2pr_build_process_one(seed, only_rp, model_env, cfg, o2_grid, "intrinsic_18only")
    if (is.null(state_note)) state_note <- attr(full_intr, "state_note")
    intr_full[[seed]] <- full_intr
    intr_18[[seed]] <- only_intr
    switch_o2 <- o2pr_switch_from_intrinsic(full_intr)
    exp_full[[seed]] <- o2pr_o2_exposure_one(seed, full_rp, model_env, cfg, burden_grid, real_o2, switch_o2, "exposure_full")
    exp_18[[seed]] <- o2pr_o2_exposure_one(seed, only_rp, model_env, cfg, burden_grid, real_o2, switch_o2, "exposure_18only")
  }
  intrinsic_full <- do.call(rbind, intr_full)
  intrinsic_18 <- do.call(rbind, intr_18)
  exposure_full <- do.call(rbind, exp_full)
  exposure_18 <- do.call(rbind, exp_18)
  o2pr_write_tsv(intrinsic_full, file.path(out_dir, "tables", "process_intrinsic_full_long.tsv"))
  o2pr_write_tsv(intrinsic_18, file.path(out_dir, "tables", "process_intrinsic_18only_long.tsv"))
  o2pr_write_tsv(exposure_full, file.path(out_dir, "tables", "process_exposure_full_long.tsv"))
  o2pr_write_tsv(exposure_18, file.path(out_dir, "tables", "process_exposure_18only_long.tsv"))
  o2pr_write_tsv(intrinsic_full, file.path(out_dir, "tables", "process_intrinsic_long.tsv"))
  o2pr_write_tsv(exposure_full, file.path(out_dir, "tables", "process_exposure_long.tsv"))
  o2pr_write_tsv(state_note, file.path(out_dir, "tables", "state_anchor_mapping.tsv"))
  ref_df <- data.frame(parameter = names(ref), value = as.numeric(ref), reference_seed = best_seed, stringsAsFactors = FALSE)
  o2pr_write_tsv(ref_df, file.path(out_dir, "tables", "reference_parameter_set.tsv"))
  list(
    intrinsic_full = intrinsic_full,
    intrinsic_18 = intrinsic_18,
    exposure_full = exposure_full,
    exposure_18 = exposure_18,
    realized_o2 = real_o2
  )
}

o2pr_switch_from_intrinsic <- function(long_df) {
  d <- long_df[long_df$feature_type == "signed" & grepl("^selection_22_vs_44_O2_", long_df$feature), , drop = FALSE]
  if (nrow(d) < 2L) return(NA_real_)
  o2 <- suppressWarnings(as.numeric(sub("^.*_O2_", "", gsub("p", ".", d$feature))))
  y <- d$raw_value
  ok <- is.finite(o2) & is.finite(y)
  o2 <- o2[ok]
  y <- y[ok]
  if (length(o2) < 2L || min(y) > 0 || max(y) < 0) return(NA_real_)
  ord <- order(o2)
  o2 <- o2[ord]
  y <- y[ord]
  as.numeric(stats::approx(y, o2, xout = 0, ties = "ordered")$y)
}
