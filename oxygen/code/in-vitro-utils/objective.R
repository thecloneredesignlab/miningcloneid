ivt_build_job_table_adapter <- function(jobs,
                                        fit_data,
                                        cohort,
                                        fallback_max_passage_days = 14) {
  all_observed <- lapply(seq_len(nrow(jobs)), function(i) {
    ids <- jobs$data_ids[[i]]
    lapply(ids, function(id) {
      out <- ivt_observed_passage_summary(fit_data[[id]])
      out$passage_id <- id
      out
    })
  })
  cohort_initial_median <- ivt_nested_observed_median(all_observed, "initial_cells", default = NA_real_)
  cohort_final_median <- ivt_nested_observed_median(all_observed, "final_cells", default = NA_real_)

  segments <- lapply(seq_len(nrow(jobs)), function(i) {
    job <- jobs[i, , drop = FALSE]
    ids <- job$data_ids[[1]]
    obs <- lapply(ids, function(id) {
      out <- ivt_observed_passage_summary(fit_data[[id]])
      out$passage_id <- id
      out
    })
    duration_use <- ivt_segment_median_value(obs, "passage_duration", default = fallback_max_passage_days)
    if (!is.finite(duration_use) || duration_use <= 0) duration_use <- as.numeric(fallback_max_passage_days)
    initial_cells_use <- as.numeric(cohort_initial_median)
    final_cells_use <- as.numeric(cohort_final_median)
    local_days <- sort(unique(c(seq(0, duration_use, by = 1), duration_use)))
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
    terminal_sim_key = NA_character_,
    n_segments = length(segments),
    segments = segments
  )
}

ivt_optimizer_spec <- function(cfg) {
  natural_tab <- read.csv(cfg$parameter_table, stringsAsFactors = FALSE)
  if (!"use_invitro_fit" %in% names(natural_tab)) {
    stop("In vitro parameter table must contain a 'use_invitro_fit' column: ", cfg$parameter_table)
  }

  expected_fit_symbols <- c(
    "lam_min", "lam_max", "k_o", "p_misseg", "k_o_mis", "gamma_loss",
    "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
    "O2_crit", "n_O", "p_mis_base"
  )
  flagged_fit_symbols <- natural_tab$param_symbol[which(tolower(trimws(as.character(natural_tab$use_invitro_fit))) %in% c("true", "t", "1", "yes", "y"))]
  if (!setequal(flagged_fit_symbols, expected_fit_symbols)) {
    stop(
      "In vitro parameter table 'use_invitro_fit' flags do not match the current optimizer definition. ",
      "Expected: ", paste(sort(expected_fit_symbols), collapse = ", "),
      ". Found: ", paste(sort(flagged_fit_symbols), collapse = ", ")
    )
  }

  row_for <- function(symbol) {
    out <- natural_tab[natural_tab$param_symbol == symbol, , drop = FALSE]
    if (nrow(out) != 1L) stop("Missing parameter row: ", symbol)
    out
  }

  positive_bound <- function(symbol, slot = c("init", "lower", "upper")) {
    slot <- match.arg(slot)
    row <- row_for(symbol)
    value <- switch(
      slot,
      init = row$init_value[[1]],
      lower = row$lower_bound[[1]],
      upper = row$upper_bound[[1]]
    )
    value <- as.numeric(value)
    if (!is.finite(value) || value <= 0) {
      stop("Parameter '", symbol, "' must have a strictly positive ", slot, " value for log-scale optimization.")
    }
    value
  }

  data.frame(
    param_name = c(
      "log10_lam_min", "log10_lam_max", "log10_k_o", "log10_p_misseg",
      "log10_k_o_mis", "log10_gamma_loss", "log10_p_wgd",
      "log10_alpha_o2", "gamma_growth", "log10_mu_hp",
      "gamma_mu", "log10_O2_crit", "n_O", "log10_p_mis_base"
    ),
    lower = c(
      log10(positive_bound("lam_min", "lower")),
      log10(positive_bound("lam_max", "lower")),
      log10(positive_bound("k_o", "lower")),
      log10(positive_bound("p_misseg", "lower")),
      log10(positive_bound("k_o_mis", "lower")),
      log10(positive_bound("gamma_loss", "lower")),
      log10(positive_bound("p_wgd", "lower")),
      log10(positive_bound("alpha_o2", "lower")),
      row_for("gamma_growth")$lower_bound[[1]],
      log10(positive_bound("mu_hp", "lower")),
      row_for("gamma_mu")$lower_bound[[1]],
      log10(positive_bound("O2_crit", "lower")),
      row_for("n_O")$lower_bound[[1]],
      log10(positive_bound("p_mis_base", "lower"))
    ),
    upper = c(
      log10(positive_bound("lam_min", "upper")),
      log10(positive_bound("lam_max", "upper")),
      log10(positive_bound("k_o", "upper")),
      log10(positive_bound("p_misseg", "upper")),
      log10(positive_bound("k_o_mis", "upper")),
      log10(positive_bound("gamma_loss", "upper")),
      log10(positive_bound("p_wgd", "upper")),
      log10(positive_bound("alpha_o2", "upper")),
      row_for("gamma_growth")$upper_bound[[1]],
      log10(positive_bound("mu_hp", "upper")),
      row_for("gamma_mu")$upper_bound[[1]],
      log10(positive_bound("O2_crit", "upper")),
      row_for("n_O")$upper_bound[[1]],
      log10(positive_bound("p_mis_base", "upper"))
    ),
    init = c(
      log10(positive_bound("lam_min", "init")),
      log10(positive_bound("lam_max", "init")),
      log10(positive_bound("k_o", "init")),
      log10(positive_bound("p_misseg", "init")),
      log10(positive_bound("k_o_mis", "init")),
      log10(positive_bound("gamma_loss", "init")),
      log10(positive_bound("p_wgd", "init")),
      log10(positive_bound("alpha_o2", "init")),
      row_for("gamma_growth")$init_value[[1]],
      log10(positive_bound("mu_hp", "init")),
      row_for("gamma_mu")$init_value[[1]],
      log10(positive_bound("O2_crit", "init")),
      row_for("n_O")$init_value[[1]],
      log10(positive_bound("p_mis_base", "init"))
    ),
    stringsAsFactors = FALSE
  )
}

ivt_run_params_to_optim_par <- function(run_params, cfg) {
  spec <- ivt_optimizer_spec(cfg)
  par_t <- setNames(spec$init, spec$param_name)

  lam_min <- as.numeric(run_params$lam_min)
  lam_max <- as.numeric(run_params$lam_max)
  if (!is.finite(lam_min) || lam_min <= 0) {
    stop("run_params$lam_min must be strictly positive for log-scale optimization.")
  }
  if (!is.finite(lam_max) || lam_max <= 0) {
    stop("run_params$lam_max must be strictly positive for log-scale optimization.")
  }

  require_positive <- function(value, name) {
    value <- as.numeric(value)
    if (!is.finite(value) || value <= 0) {
      stop("run_params$", name, " must be strictly positive for log-scale optimization.")
    }
    value
  }

  par_t[["log10_lam_min"]] <- log10(lam_min)
  par_t[["log10_lam_max"]] <- log10(lam_max)
  par_t[["log10_k_o"]] <- log10(require_positive(run_params$k_o, "k_o"))
  par_t[["log10_p_misseg"]] <- log10(require_positive(run_params$p_misseg, "p_misseg"))
  par_t[["log10_k_o_mis"]] <- log10(require_positive(run_params$k_o_mis, "k_o_mis"))
  par_t[["log10_gamma_loss"]] <- log10(require_positive(run_params$gamma_loss, "gamma_loss"))
  par_t[["log10_p_wgd"]] <- log10(require_positive(run_params$p_wgd, "p_wgd"))
  par_t[["log10_alpha_o2"]] <- log10(require_positive(run_params$alpha_o2, "alpha_o2"))
  par_t[["gamma_growth"]] <- as.numeric(run_params$gamma_growth)
  par_t[["log10_mu_hp"]] <- log10(require_positive(run_params$mu_hp, "mu_hp"))
  par_t[["gamma_mu"]] <- as.numeric(run_params$gamma_mu)
  par_t[["log10_O2_crit"]] <- log10(require_positive(run_params$O2_crit, "O2_crit"))
  par_t[["n_O"]] <- as.numeric(run_params$n_O)
  par_t[["log10_p_mis_base"]] <- log10(require_positive(run_params$p_mis_base, "p_mis_base"))

  pmin(pmax(par_t, spec$lower), spec$upper)
}

ivt_optim_par_to_run_params <- function(par_t, cfg) {
  par_t <- as.numeric(par_t)
  spec <- ivt_optimizer_spec(cfg)
  if (length(par_t) != nrow(spec)) {
    stop("Optimization parameter length mismatch.")
  }
  names(par_t) <- spec$param_name

  run_params <- ivt_load_default_run_params(cfg)
  run_params$lam_min <- 10^par_t[["log10_lam_min"]]
  run_params$lam_max <- 10^par_t[["log10_lam_max"]]
  run_params$k_o <- 10^par_t[["log10_k_o"]]
  run_params$p_misseg <- 10^par_t[["log10_p_misseg"]]
  run_params$k_o_mis <- 10^par_t[["log10_k_o_mis"]]
  run_params$gamma_loss <- 10^par_t[["log10_gamma_loss"]]
  run_params$p_wgd <- 10^par_t[["log10_p_wgd"]]
  run_params$alpha_o2 <- 10^par_t[["log10_alpha_o2"]]
  run_params$gamma_growth <- par_t[["gamma_growth"]]
  run_params$mu_hp <- 10^par_t[["log10_mu_hp"]]
  run_params$gamma_mu <- par_t[["gamma_mu"]]
  run_params$O2_crit <- 10^par_t[["log10_O2_crit"]]
  run_params$n_O <- par_t[["n_O"]]
  run_params$p_mis_base <- 10^par_t[["log10_p_mis_base"]]

  normalize_run_params_common(run_params, cfg = cfg)
}

ivt_objective_components <- function(run_params,
                                     fit_objects,
                                     cfg,
                                     fallback_max_passage_days = 14,
                                     growth_weight = 1,
                                     ploidy_weight = 1) {
  model_core <- build_model_core(cfg = cfg)

  adapter_2N <- ivt_build_job_table_adapter(
    jobs = fit_objects$jobs_2N,
    fit_data = fit_objects$fit_data,
    cohort = "2N",
    fallback_max_passage_days = fallback_max_passage_days
  )
  adapter_4N <- ivt_build_job_table_adapter(
    jobs = fit_objects$jobs_4N,
    fit_data = fit_objects$fit_data,
    cohort = "4N",
    fallback_max_passage_days = fallback_max_passage_days
  )

  run_2N <- ivt_run_lineage(adapter_2N, cfg = cfg, run_params = run_params, model_core = model_core)
  run_4N <- ivt_run_lineage(adapter_4N, cfg = cfg, run_params = run_params, model_core = model_core)

  sum_2N <- ivt_collect_lineage_summary(run_2N, fit_objects$fit_data)
  sum_4N <- ivt_collect_lineage_summary(run_4N, fit_objects$fit_data)
  summary_df <- dplyr::bind_rows(sum_2N, sum_4N)

  growth_df <- summary_df[is.finite(summary_df$observed_growth) & is.finite(summary_df$predicted_growth_rate), , drop = FALSE]
  kary_df <- summary_df[is.finite(summary_df$observed_mean_kary_N) & is.finite(summary_df$predicted_mean_kary_N), , drop = FALSE]

  growth_sse <- if (nrow(growth_df) > 0L) {
    mean((growth_df$predicted_growth_rate - growth_df$observed_growth)^2)
  } else {
    0
  }
  kary_sse <- if (nrow(kary_df) > 0L) {
    mean((kary_df$predicted_mean_kary_N - kary_df$observed_mean_kary_N)^2)
  } else {
    0
  }

  total <- as.numeric(growth_weight) * growth_sse + as.numeric(ploidy_weight) * kary_sse

  list(
    objective = total,
    growth_sse = growth_sse,
    kary_sse = kary_sse,
    n_growth = nrow(growth_df),
    n_kary = nrow(kary_df),
    summary = summary_df,
    run_2N = run_2N,
    run_4N = run_4N
  )
}

ivt_objective <- function(run_params,
                          fit_objects,
                          cfg,
                          fallback_max_passage_days = 14,
                          growth_weight = 1,
                          ploidy_weight = 1) {
  ivt_objective_components(
    run_params = run_params,
    fit_objects = fit_objects,
    cfg = cfg,
    fallback_max_passage_days = fallback_max_passage_days,
    growth_weight = growth_weight,
    ploidy_weight = ploidy_weight
  )$objective
}
