ivt_build_job_table_adapter <- function(jobs,
                                        fit_data,
                                        cohort,
                                        fallback_max_passage_days = 14) {
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
    initial_cells_use <- ivt_segment_median_value(obs, "initial_cells", default = NA_real_)
    final_cells_use <- ivt_segment_median_value(obs, "final_cells", default = NA_real_)
    list(
      segment_id = job$sim_key[[1]],
      parent_segment_id = job$parent_key[[1]],
      cohort = cohort,
      oxygen_pct = as.numeric(job$oxygen[[1]]),
      duration_days = duration_use,
      initial_cells = initial_cells_use,
      final_cells = final_cells_use,
      obs_days_local = c(0, duration_use),
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

  row_for <- function(symbol) {
    out <- natural_tab[natural_tab$param_symbol == symbol, , drop = FALSE]
    if (nrow(out) != 1L) stop("Missing parameter row: ", symbol)
    out
  }

  delta_lam_bound <- function(slot = c("init", "lower", "upper")) {
    slot <- match.arg(slot)
    lam_min <- row_for("lam_min")
    lam_max <- row_for("lam_max")
    gap <- switch(
      slot,
      init = lam_max$init_value[[1]] - lam_min$init_value[[1]],
      lower = lam_max$lower_bound[[1]] - lam_min$upper_bound[[1]],
      upper = lam_max$upper_bound[[1]] - lam_min$lower_bound[[1]]
    )
    log(max(as.numeric(gap), 1e-8))
  }

  data.frame(
    param_name = c(
      "log10_lam_min", "delta_lam", "log10_k_o", "log10_p_misseg",
      "log10_k_o_mis", "log10_gamma_loss", "log10_p_wgd",
      "log10_alpha_o2", "gamma_growth", "log10_mu_hp",
      "gamma_mu", "log10_O2_crit", "n_O", "log10_p_mis_base"
    ),
    lower = c(
      log10(row_for("lam_min")$lower_bound[[1]]),
      delta_lam_bound("lower"),
      log10(row_for("k_o")$lower_bound[[1]]),
      log10(row_for("p_misseg")$lower_bound[[1]]),
      log10(row_for("k_o_mis")$lower_bound[[1]]),
      log10(row_for("gamma_loss")$lower_bound[[1]]),
      log10(row_for("p_wgd")$lower_bound[[1]]),
      log10(row_for("alpha_o2")$lower_bound[[1]]),
      row_for("gamma_growth")$lower_bound[[1]],
      log10(row_for("mu_hp")$lower_bound[[1]]),
      row_for("gamma_mu")$lower_bound[[1]],
      log10(max(row_for("O2_crit")$lower_bound[[1]], 1e-8)),
      row_for("n_O")$lower_bound[[1]],
      log10(row_for("p_mis_base")$lower_bound[[1]])
    ),
    upper = c(
      log10(row_for("lam_min")$upper_bound[[1]]),
      delta_lam_bound("upper"),
      log10(row_for("k_o")$upper_bound[[1]]),
      log10(row_for("p_misseg")$upper_bound[[1]]),
      log10(row_for("k_o_mis")$upper_bound[[1]]),
      log10(row_for("gamma_loss")$upper_bound[[1]]),
      log10(row_for("p_wgd")$upper_bound[[1]]),
      log10(row_for("alpha_o2")$upper_bound[[1]]),
      row_for("gamma_growth")$upper_bound[[1]],
      log10(row_for("mu_hp")$upper_bound[[1]]),
      row_for("gamma_mu")$upper_bound[[1]],
      log10(max(row_for("O2_crit")$upper_bound[[1]], 1e-8)),
      row_for("n_O")$upper_bound[[1]],
      log10(row_for("p_mis_base")$upper_bound[[1]])
    ),
    init = c(
      log10(row_for("lam_min")$init_value[[1]]),
      delta_lam_bound("init"),
      log10(row_for("k_o")$init_value[[1]]),
      log10(row_for("p_misseg")$init_value[[1]]),
      log10(row_for("k_o_mis")$init_value[[1]]),
      log10(row_for("gamma_loss")$init_value[[1]]),
      log10(row_for("p_wgd")$init_value[[1]]),
      log10(row_for("alpha_o2")$init_value[[1]]),
      row_for("gamma_growth")$init_value[[1]],
      log10(row_for("mu_hp")$init_value[[1]]),
      row_for("gamma_mu")$init_value[[1]],
      log10(max(row_for("O2_crit")$init_value[[1]], 1e-8)),
      row_for("n_O")$init_value[[1]],
      log10(row_for("p_mis_base")$init_value[[1]])
    ),
    stringsAsFactors = FALSE
  )
}

ivt_run_params_to_optim_par <- function(run_params, cfg) {
  spec <- ivt_optimizer_spec(cfg)
  par_t <- setNames(spec$init, spec$param_name)

  lam_min <- as.numeric(run_params$lam_min)
  lam_max <- as.numeric(run_params$lam_max)
  delta_lam <- max(lam_max - lam_min, 1e-8)

  par_t[["log10_lam_min"]] <- log10(max(lam_min, 1e-8))
  par_t[["delta_lam"]] <- log(delta_lam)
  par_t[["log10_k_o"]] <- log10(max(as.numeric(run_params$k_o), 1e-8))
  par_t[["log10_p_misseg"]] <- log10(max(as.numeric(run_params$p_misseg), 1e-12))
  par_t[["log10_k_o_mis"]] <- log10(max(as.numeric(run_params$k_o_mis), 1e-12))
  par_t[["log10_gamma_loss"]] <- log10(max(as.numeric(run_params$gamma_loss), 1e-12))
  par_t[["log10_p_wgd"]] <- log10(max(as.numeric(run_params$p_wgd), 1e-12))
  par_t[["log10_alpha_o2"]] <- log10(max(as.numeric(run_params$alpha_o2), 1e-12))
  par_t[["gamma_growth"]] <- as.numeric(run_params$gamma_growth)
  par_t[["log10_mu_hp"]] <- log10(max(as.numeric(run_params$mu_hp), 1e-12))
  par_t[["gamma_mu"]] <- as.numeric(run_params$gamma_mu)
  par_t[["log10_O2_crit"]] <- log10(max(as.numeric(run_params$O2_crit), 1e-8))
  par_t[["n_O"]] <- as.numeric(run_params$n_O)
  par_t[["log10_p_mis_base"]] <- log10(max(as.numeric(run_params$p_mis_base), 1e-12))

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
  lam_min <- 10^par_t[["log10_lam_min"]]
  run_params$lam_min <- lam_min
  run_params$lam_max <- lam_min + exp(par_t[["delta_lam"]])
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
  ploidy_df <- summary_df[is.finite(summary_df$observed_mean_ploidy) & is.finite(summary_df$predicted_mean_ploidy), , drop = FALSE]

  growth_sse <- if (nrow(growth_df) > 0L) {
    mean((growth_df$predicted_growth_rate - growth_df$observed_growth)^2)
  } else {
    0
  }
  ploidy_sse <- if (nrow(ploidy_df) > 0L) {
    mean((ploidy_df$predicted_mean_ploidy - ploidy_df$observed_mean_ploidy)^2)
  } else {
    0
  }

  total <- as.numeric(growth_weight) * growth_sse + as.numeric(ploidy_weight) * ploidy_sse

  list(
    objective = total,
    growth_sse = growth_sse,
    ploidy_sse = ploidy_sse,
    n_growth = nrow(growth_df),
    n_ploidy = nrow(ploidy_df),
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
