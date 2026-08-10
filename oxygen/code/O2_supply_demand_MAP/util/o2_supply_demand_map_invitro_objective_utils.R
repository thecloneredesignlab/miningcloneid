# Canonical in-vitro objective, likelihood, and parameter-transform helpers.

ivt_build_job_table_adapter <- function(jobs,
                                        fit_data,
                                        cohort,
                                        fallback_max_passage_days = 14,
                                        passage_time_tolerance_days = 1) {
  .ivt_build_independent_scenario_adapter(
    jobs = jobs,
    fit_data = fit_data,
    cohort = cohort,
    obs_days_local = NULL,
    passage_time_tolerance_days = passage_time_tolerance_days
  )
}

ivt_optimizer_spec <- function(cfg) {
  natural_tab <- read.csv(cfg$parameter_table, stringsAsFactors = FALSE)
  if (!"estimate" %in% names(natural_tab)) {
    stop("In vitro parameter table must contain an 'estimate' column: ", cfg$parameter_table)
  }
  if (!"use_invitro_fit" %in% names(natural_tab)) {
    stop("In vitro parameter table must contain a 'use_invitro_fit' column: ", cfg$parameter_table)
  }
  expected_fit_symbols <- c(
    "lam_max", "p_misseg", "k_o_mis",
    "buffer_smax", "buffer_beta", "buffer_n_exp",
    "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
    "O2_crit", "n_O", "p_mis_base", "sigma_growth", "sigma_kary",
    "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N"
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

  loss_param_name <- c("buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp")
  loss_lower <- c(
    as.numeric(row_for("buffer_smax")$lower_bound[[1]]),
    log10(positive_bound("buffer_beta", "lower")),
    log10(positive_bound("buffer_n_exp", "lower"))
  )
  loss_upper <- c(
    as.numeric(row_for("buffer_smax")$upper_bound[[1]]),
    log10(positive_bound("buffer_beta", "upper")),
    log10(positive_bound("buffer_n_exp", "upper"))
  )
  loss_init <- c(
    as.numeric(row_for("buffer_smax")$init_value[[1]]),
    log10(positive_bound("buffer_beta", "init")),
    log10(positive_bound("buffer_n_exp", "init"))
  )

  natural_name <- c(
    "lam_max", "p_misseg", "k_o_mis",
    "buffer_smax", "buffer_beta", "buffer_n_exp", "p_wgd",
    "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
    "O2_crit", "n_O", "p_mis_base", "sigma_growth", "sigma_kary",
    "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N"
  )
  estimate_flag <- vapply(natural_name, function(symbol) {
    value <- tolower(trimws(as.character(row_for(symbol)$estimate[[1]])))
    value %in% c("true", "t", "1", "yes", "y")
  }, logical(1))

  data.frame(
    param_name = c(
      "log10_lam_max", "log10_p_misseg",
      "log10_k_o_mis", loss_param_name, "log10_p_wgd",
      "log10_alpha_o2", "gamma_growth", "log10_mu_hp",
      "gamma_mu", "log10_O2_crit", "n_O", "log10_p_mis_base",
      "log10_sigma_growth", "log10_sigma_kary", "init_mean_2N", "log10_init_sd_2N",
      "init_mean_4N", "log10_init_sd_4N"
    ),
    lower = c(
      log10(positive_bound("lam_max", "lower")),
      log10(positive_bound("p_misseg", "lower")),
      log10(positive_bound("k_o_mis", "lower")),
      loss_lower,
      log10(positive_bound("p_wgd", "lower")),
      log10(positive_bound("alpha_o2", "lower")),
      row_for("gamma_growth")$lower_bound[[1]],
      log10(positive_bound("mu_hp", "lower")),
      row_for("gamma_mu")$lower_bound[[1]],
      log10(positive_bound("O2_crit", "lower")),
      row_for("n_O")$lower_bound[[1]],
      log10(positive_bound("p_mis_base", "lower")),
      log10(positive_bound("sigma_growth", "lower")),
      log10(positive_bound("sigma_kary", "lower")),
      row_for("init_mean_2N")$lower_bound[[1]],
      log10(positive_bound("init_sd_2N", "lower")),
      row_for("init_mean_4N")$lower_bound[[1]],
      log10(positive_bound("init_sd_4N", "lower"))
    ),
    upper = c(
      log10(positive_bound("lam_max", "upper")),
      log10(positive_bound("p_misseg", "upper")),
      log10(positive_bound("k_o_mis", "upper")),
      loss_upper,
      log10(positive_bound("p_wgd", "upper")),
      log10(positive_bound("alpha_o2", "upper")),
      row_for("gamma_growth")$upper_bound[[1]],
      log10(positive_bound("mu_hp", "upper")),
      row_for("gamma_mu")$upper_bound[[1]],
      log10(positive_bound("O2_crit", "upper")),
      row_for("n_O")$upper_bound[[1]],
      log10(positive_bound("p_mis_base", "upper")),
      log10(positive_bound("sigma_growth", "upper")),
      log10(positive_bound("sigma_kary", "upper")),
      row_for("init_mean_2N")$upper_bound[[1]],
      log10(positive_bound("init_sd_2N", "upper")),
      row_for("init_mean_4N")$upper_bound[[1]],
      log10(positive_bound("init_sd_4N", "upper"))
    ),
    init = c(
      log10(positive_bound("lam_max", "init")),
      log10(positive_bound("p_misseg", "init")),
      log10(positive_bound("k_o_mis", "init")),
      loss_init,
      log10(positive_bound("p_wgd", "init")),
      log10(positive_bound("alpha_o2", "init")),
      row_for("gamma_growth")$init_value[[1]],
      log10(positive_bound("mu_hp", "init")),
      row_for("gamma_mu")$init_value[[1]],
      log10(positive_bound("O2_crit", "init")),
      row_for("n_O")$init_value[[1]],
      log10(positive_bound("p_mis_base", "init")),
      log10(positive_bound("sigma_growth", "init")),
      log10(positive_bound("sigma_kary", "init")),
      row_for("init_mean_2N")$init_value[[1]],
      log10(positive_bound("init_sd_2N", "init")),
      row_for("init_mean_4N")$init_value[[1]],
      log10(positive_bound("init_sd_4N", "init"))
    ),
    natural_name = natural_name,
    estimate = estimate_flag,
    stringsAsFactors = FALSE
  )
}

.ivt_buffer_soft_prior <- function(run_params,
                                  weight = 0,
                                  center_smax = 0.98,
                                  sd_smax = 0.10,
                                  center_log10_beta = log10(0.7),
                                  sd_log10_beta = 0.30,
                                  center_log10_n_exp = log10(5),
                                  sd_log10_n_exp = 0.30) {
  scalar_finite <- function(value, name) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) {
      stop(name, " must be one finite numeric value.")
    }
    value
  }
  positive_scalar <- function(value, name) {
    value <- scalar_finite(value, name)
    if (value <= 0) stop(name, " must be strictly positive.")
    value
  }

  weight_use <- scalar_finite(weight, "buffer prior weight")
  if (weight_use < 0) stop("buffer prior weight must be non-negative.")
  center_smax_use <- scalar_finite(center_smax, "buffer_smax prior center")
  if (center_smax_use < 0 || center_smax_use > 1) {
    stop("buffer_smax prior center must lie in [0, 1].")
  }
  sd_smax_use <- positive_scalar(sd_smax, "buffer_smax prior SD")
  center_log10_beta_use <- scalar_finite(
    center_log10_beta,
    "log10(buffer_beta) prior center"
  )
  sd_log10_beta_use <- positive_scalar(
    sd_log10_beta,
    "log10(buffer_beta) prior SD"
  )
  center_log10_n_exp_use <- scalar_finite(
    center_log10_n_exp,
    "log10(buffer_n_exp) prior center"
  )
  sd_log10_n_exp_use <- positive_scalar(
    sd_log10_n_exp,
    "log10(buffer_n_exp) prior SD"
  )

  value_smax <- scalar_finite(run_params$buffer_smax, "run_params$buffer_smax")
  value_beta <- positive_scalar(run_params$buffer_beta, "run_params$buffer_beta")
  value_n_exp <- positive_scalar(run_params$buffer_n_exp, "run_params$buffer_n_exp")
  transformed_value <- c(
    buffer_smax = value_smax,
    log10_buffer_beta = log10(value_beta),
    log10_buffer_n_exp = log10(value_n_exp)
  )
  center <- c(
    buffer_smax = center_smax_use,
    log10_buffer_beta = center_log10_beta_use,
    log10_buffer_n_exp = center_log10_n_exp_use
  )
  prior_sd <- c(
    buffer_smax = sd_smax_use,
    log10_buffer_beta = sd_log10_beta_use,
    log10_buffer_n_exp = sd_log10_n_exp_use
  )
  z_score <- (transformed_value - center) / prior_sd
  raw_term <- 0.5 * z_score^2
  raw_total <- sum(raw_term)
  weighted_total <- weight_use * raw_total
  terms <- data.frame(
    parameter = names(transformed_value),
    transformed_value = as.numeric(transformed_value),
    prior_center = as.numeric(center),
    prior_sd = as.numeric(prior_sd),
    z_score = as.numeric(z_score),
    raw_penalty = as.numeric(raw_term),
    prior_weight = rep(weight_use, length(raw_term)),
    weighted_penalty = as.numeric(weight_use * raw_term),
    stringsAsFactors = FALSE
  )
  list(
    weight = weight_use,
    raw_penalty = raw_total,
    weighted_penalty = weighted_total,
    terms = terms
  )
}

ivt_run_params_to_optim_par <- function(run_params, cfg) {
  spec <- ivt_optimizer_spec(cfg)
  par_t <- setNames(spec$init, spec$param_name)

  lam_max <- as.numeric(run_params$lam_max)
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

  par_t[["log10_lam_max"]] <- log10(lam_max)
  par_t[["log10_p_misseg"]] <- log10(require_positive(run_params$p_misseg, "p_misseg"))
  par_t[["log10_k_o_mis"]] <- log10(require_positive(run_params$k_o_mis, "k_o_mis"))
  if ("buffer_smax" %in% names(par_t)) {
    smax <- as.numeric(run_params$buffer_smax)
    if (!is.finite(smax) || smax < 0 || smax > 1) {
      stop("run_params$buffer_smax must be in [0,1].")
    }
    par_t[["buffer_smax"]] <- smax
  }
  if ("log10_buffer_beta" %in% names(par_t)) {
    par_t[["log10_buffer_beta"]] <- log10(require_positive(run_params$buffer_beta, "buffer_beta"))
  }
  if ("log10_buffer_n_exp" %in% names(par_t)) {
    par_t[["log10_buffer_n_exp"]] <- log10(require_positive(run_params$buffer_n_exp, "buffer_n_exp"))
  }
  par_t[["log10_p_wgd"]] <- log10(require_positive(run_params$p_wgd, "p_wgd"))
  par_t[["log10_alpha_o2"]] <- log10(require_positive(run_params$alpha_o2, "alpha_o2"))
  par_t[["gamma_growth"]] <- as.numeric(run_params$gamma_growth)
  par_t[["log10_mu_hp"]] <- log10(require_positive(run_params$mu_hp, "mu_hp"))
  par_t[["gamma_mu"]] <- as.numeric(run_params$gamma_mu)
  par_t[["log10_O2_crit"]] <- log10(require_positive(run_params$O2_crit, "O2_crit"))
  par_t[["n_O"]] <- as.numeric(run_params$n_O)
  par_t[["log10_p_mis_base"]] <- log10(require_positive(run_params$p_mis_base, "p_mis_base"))
  par_t[["log10_sigma_growth"]] <- log10(require_positive(run_params$sigma_growth, "sigma_growth"))
  par_t[["log10_sigma_kary"]] <- log10(require_positive(run_params$sigma_kary, "sigma_kary"))
  par_t[["init_mean_2N"]] <- as.numeric(run_params$init_mean_2N)
  par_t[["log10_init_sd_2N"]] <- log10(require_positive(run_params$init_sd_2N, "init_sd_2N"))
  par_t[["init_mean_4N"]] <- as.numeric(run_params$init_mean_4N)
  par_t[["log10_init_sd_4N"]] <- log10(require_positive(run_params$init_sd_4N, "init_sd_4N"))

  pmin(pmax(par_t, spec$lower), spec$upper)
}

ivt_optim_par_to_run_params <- function(par_t, cfg) {
  spec <- ivt_optimizer_spec(cfg)
  spec_init <- setNames(as.numeric(spec$init), spec$param_name)

  if (is.null(names(par_t))) {
    par_t <- as.numeric(par_t)
    if (length(par_t) != nrow(spec)) {
      stop("Optimization parameter length mismatch.")
    }
    names(par_t) <- spec$param_name
  } else {
    par_names <- names(par_t)
    par_t <- as.numeric(par_t)
    names(par_t) <- par_names
    extra_names <- setdiff(par_names, spec$param_name)
    if (length(extra_names) > 0L) {
      stop("Unknown optimization parameters: ", paste(extra_names, collapse = ", "))
    }
    spec_init[par_names] <- par_t
    par_t <- spec_init
  }

  run_params <- ivt_load_default_run_params(cfg)
  run_params$lam_max <- 10^par_t[["log10_lam_max"]]
  run_params$p_misseg <- 10^par_t[["log10_p_misseg"]]
  run_params$k_o_mis <- 10^par_t[["log10_k_o_mis"]]
  if ("buffer_smax" %in% names(par_t)) run_params$buffer_smax <- par_t[["buffer_smax"]]
  if ("log10_buffer_beta" %in% names(par_t)) run_params$buffer_beta <- 10^par_t[["log10_buffer_beta"]]
  if ("log10_buffer_n_exp" %in% names(par_t)) run_params$buffer_n_exp <- 10^par_t[["log10_buffer_n_exp"]]
  run_params$p_wgd <- 10^par_t[["log10_p_wgd"]]
  run_params$alpha_o2 <- 10^par_t[["log10_alpha_o2"]]
  run_params$gamma_growth <- par_t[["gamma_growth"]]
  run_params$mu_hp <- 10^par_t[["log10_mu_hp"]]
  run_params$gamma_mu <- par_t[["gamma_mu"]]
  run_params$O2_crit <- 10^par_t[["log10_O2_crit"]]
  run_params$n_O <- par_t[["n_O"]]
  run_params$p_mis_base <- 10^par_t[["log10_p_mis_base"]]
  run_params$sigma_growth <- 10^par_t[["log10_sigma_growth"]]
  run_params$sigma_kary <- 10^par_t[["log10_sigma_kary"]]
  run_params$init_mean_2N <- par_t[["init_mean_2N"]]
  run_params$init_sd_2N <- 10^par_t[["log10_init_sd_2N"]]
  run_params$init_mean_4N <- par_t[["init_mean_4N"]]
  run_params$init_sd_4N <- 10^par_t[["log10_init_sd_4N"]]

  normalize_run_params_common(run_params, cfg = cfg)
}

.ivt_growth_measurement_summary <- function(summary_df,
                                            fit_data,
                                            daily_counts,
                                            day_tolerance = 1e-8) {
  required_summary <- c(
    "passage_id", "segment_id", "last_observation_day", "selected_day",
    "observed_live_cells_at_observation", "observed_initial_cells",
    "predicted_initial_cells"
  )
  missing_summary <- setdiff(required_summary, names(summary_df))
  if (length(missing_summary) > 0L) {
    stop(
      "Growth passage summary is missing columns: ",
      paste(missing_summary, collapse = ", ")
    )
  }
  required_daily <- c("segment_id", "passage_id", "day", "live_cells")
  missing_daily <- setdiff(required_daily, names(daily_counts))
  if (length(missing_daily) > 0L) {
    stop(
      "Daily model counts are missing growth-likelihood columns: ",
      paste(missing_daily, collapse = ", ")
    )
  }
  tolerance_use <- suppressWarnings(as.numeric(day_tolerance))
  if (length(tolerance_use) != 1L || !is.finite(tolerance_use) || tolerance_use < 0) {
    stop("day_tolerance must be one finite non-negative value.")
  }

  rows <- lapply(seq_len(nrow(summary_df)), function(i) {
    passage_row <- summary_df[i, , drop = FALSE]
    passage_id <- as.character(passage_row$passage_id[[1]])
    segment_id <- as.character(passage_row$segment_id[[1]])
    entry <- fit_data[[passage_id]]
    timecourse <- if (is.null(entry)) NULL else entry$growth_timecourse
    if (!is.data.frame(timecourse)) {
      timecourse <- data.frame(
        growth_observation_id = paste0(passage_id, "__endpoint"),
        observation_day = suppressWarnings(as.numeric(
          passage_row$last_observation_day[[1]]
        )),
        observed_live_cells = suppressWarnings(as.numeric(
          passage_row$observed_live_cells_at_observation[[1]]
        )),
        growth_data_source = "summary_endpoint_fallback",
        stringsAsFactors = FALSE
      )
    }
    required_timecourse <- c(
      "growth_observation_id", "observation_day", "observed_live_cells",
      "growth_data_source"
    )
    missing_timecourse <- setdiff(required_timecourse, names(timecourse))
    if (length(missing_timecourse) > 0L) {
      stop(
        "Growth timecourse for passage ", passage_id, " is missing columns: ",
        paste(missing_timecourse, collapse = ", ")
      )
    }
    timecourse$observation_day <- suppressWarnings(as.numeric(timecourse$observation_day))
    timecourse$observed_live_cells <- suppressWarnings(as.numeric(
      timecourse$observed_live_cells
    ))
    timecourse <- timecourse[
      is.finite(timecourse$observation_day) & timecourse$observation_day > 0 &
        is.finite(timecourse$observed_live_cells) & timecourse$observed_live_cells > 0,
      ,
      drop = FALSE
    ]
    if (nrow(timecourse) == 0L) return(NULL)
    timecourse <- timecourse[
      order(timecourse$observation_day, timecourse$growth_observation_id),
      ,
      drop = FALSE
    ]
    last_observation_day <- suppressWarnings(as.numeric(
      passage_row$last_observation_day[[1]]
    ))
    if (!is.finite(last_observation_day) ||
        abs(max(timecourse$observation_day) - last_observation_day) > tolerance_use) {
      stop(
        "Growth timecourse endpoint does not match last_observation_day for passage ",
        passage_id, "."
      )
    }
    selected_day <- suppressWarnings(as.numeric(passage_row$selected_day[[1]]))
    if (!is.finite(selected_day) ||
        selected_day + tolerance_use < max(timecourse$observation_day)) {
      stop(
        "Passage ", passage_id,
        " was selected before its last growth observation day."
      )
    }

    predicted <- vapply(timecourse$observation_day, function(day_value) {
      hit <- which(
        as.character(daily_counts$segment_id) == segment_id &
          as.character(daily_counts$passage_id) == passage_id &
          abs(suppressWarnings(as.numeric(daily_counts$day)) - day_value) <= tolerance_use
      )
      if (length(hit) != 1L) {
        stop(
          "Expected exactly one model count for passage ", passage_id,
          " at day ", format(day_value, trim = TRUE),
          "; found ", length(hit), "."
        )
      }
      suppressWarnings(as.numeric(daily_counts$live_cells[[hit]]))
    }, numeric(1))
    if (any(!is.finite(predicted)) || any(predicted <= 0)) {
      stop("Model live-cell predictions must be positive for passage ", passage_id, ".")
    }

    out <- passage_row[rep(1L, nrow(timecourse)), , drop = FALSE]
    row.names(out) <- NULL
    out$growth_observation_id <- as.character(timecourse$growth_observation_id)
    out$growth_data_source <- as.character(timecourse$growth_data_source)
    out$observation_day <- timecourse$observation_day
    out$is_last_observation <-
      abs(out$observation_day - last_observation_day) <= tolerance_use
    out$n_growth_timepoints_in_passage <- nrow(timecourse)
    out$observed_live_cells_at_observation <- timecourse$observed_live_cells
    out$predicted_live_cells_at_observation <- predicted
    out$measurement_day_cell_count_residual <-
      predicted - timecourse$observed_live_cells
    out$log_live_cell_residual <-
      log(predicted / timecourse$observed_live_cells)
    observed_initial_cells <- suppressWarnings(as.numeric(
      passage_row$observed_initial_cells[[1]]
    ))
    predicted_initial_cells <- suppressWarnings(as.numeric(
      passage_row$predicted_initial_cells[[1]]
    ))
    if (!is.finite(observed_initial_cells) || observed_initial_cells <= 0 ||
        !is.finite(predicted_initial_cells) || predicted_initial_cells <= 0) {
      stop("Growth passage ", passage_id, " has invalid initial live-cell counts.")
    }
    out$observed_log_population_change <-
      log(out$observed_live_cells_at_observation / observed_initial_cells)
    out$predicted_log_population_change <-
      log(out$predicted_live_cells_at_observation / predicted_initial_cells)
    out$observed_timepoint_average_growth_rate <-
      out$observed_log_population_change / out$observation_day
    out$predicted_timepoint_average_growth_rate <-
      out$predicted_log_population_change / out$observation_day
    out$growth_objective_role <-
      "input_to_passage_rate_estimator_not_independent_likelihood_unit"
    out
  })
  dplyr::bind_rows(rows)
}

.ivt_growth_passage_loglik_df <- function(growth_df) {
  if (is.null(growth_df) || !is.data.frame(growth_df) || nrow(growth_df) == 0L) {
    return(data.frame())
  }
  required <- c(
    "passage_id", "cohort", "lineage_id", "scenario_id",
    "n_growth_timepoints", "loglik"
  )
  missing <- setdiff(required, names(growth_df))
  if (length(missing) > 0L) {
    stop("Growth passage likelihood is missing columns: ", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(as.character(growth_df$passage_id))) {
    stop("Growth likelihood must contain exactly one row per passage.")
  }
  values <- suppressWarnings(as.numeric(growth_df$loglik))
  if (any(!is.finite(values))) {
    stop("Growth passage log-likelihood values must be finite.")
  }
  data.frame(
    passage_id = as.character(growth_df$passage_id),
    cohort = as.character(growth_df$cohort),
    lineage_id = as.character(growth_df$lineage_id),
    scenario_id = as.character(growth_df$scenario_id),
    n_growth_timepoints = as.integer(growth_df$n_growth_timepoints),
    growth_loglik_sum = values,
    mean_loglik = values,
    stringsAsFactors = FALSE
  )
}

ivt_growth_loglik_df <- function(summary_df, sigma_growth) {
  sigma_use <- as.numeric(sigma_growth)
  if (!is.finite(sigma_use) || sigma_use <= 0) {
    stop("sigma_growth must be strictly positive.")
  }

  required <- c(
    "passage_id", "cohort", "lineage_id", "scenario_id",
    "observation_day", "observed_initial_cells", "predicted_initial_cells",
    "observed_live_cells_at_observation",
    "predicted_live_cells_at_observation"
  )
  missing <- setdiff(required, names(summary_df))
  if (length(missing)) {
    stop(
      "Growth summary is missing passage-rate columns: ",
      paste(missing, collapse = ", ")
    )
  }
  observed_cells <- suppressWarnings(as.numeric(
    summary_df$observed_live_cells_at_observation
  ))
  predicted_cells <- suppressWarnings(as.numeric(
    summary_df$predicted_live_cells_at_observation
  ))
  observation_days <- suppressWarnings(as.numeric(summary_df$observation_day))
  observed_initial <- suppressWarnings(as.numeric(summary_df$observed_initial_cells))
  predicted_initial <- suppressWarnings(as.numeric(summary_df$predicted_initial_cells))
  keep <- is.finite(observed_cells) & observed_cells > 0 &
    is.finite(predicted_cells) & predicted_cells > 0 &
    is.finite(observation_days) & observation_days > 0 &
    is.finite(observed_initial) & observed_initial > 0 &
    is.finite(predicted_initial) & predicted_initial > 0
  work <- summary_df[keep, , drop = FALSE]
  if (nrow(work) == 0L) {
    out <- summary_df[FALSE, , drop = FALSE]
    out$n_growth_timepoints <- integer(0)
    out$observed_passage_growth_rate <- numeric(0)
    out$predicted_passage_growth_rate <- numeric(0)
    out$growth_rate_residual <- numeric(0)
    out$sigma_growth_rate <- numeric(0)
    out$growth_rate_estimator <- character(0)
    out$growth_likelihood_scale <- character(0)
    out$loglik <- numeric(0)
    return(out)
  }

  passage_ids <- as.character(work$passage_id)
  passage_groups <- split(
    seq_len(nrow(work)),
    factor(passage_ids, levels = unique(passage_ids))
  )
  point_columns <- c(
    "growth_observation_id", "growth_data_source", "observation_day",
    "is_last_observation", "n_growth_timepoints_in_passage",
    "observed_live_cells_at_observation",
    "predicted_live_cells_at_observation",
    "measurement_day_cell_count_residual", "log_live_cell_residual",
    "observed_log_population_change", "predicted_log_population_change",
    "observed_timepoint_average_growth_rate",
    "predicted_timepoint_average_growth_rate", "growth_objective_role"
  )
  passage_rows <- lapply(passage_groups, function(idx) {
    passage <- work[idx, , drop = FALSE]
    identity_columns <- c("passage_id", "cohort", "lineage_id", "scenario_id")
    identities <- lapply(identity_columns, function(column) {
      unique(as.character(passage[[column]]))
    })
    if (any(lengths(identities) != 1L)) {
      stop("Growth timepoints have inconsistent passage identity fields.")
    }
    unique_observed_initial <- unique(suppressWarnings(as.numeric(
      passage$observed_initial_cells
    )))
    unique_predicted_initial <- unique(suppressWarnings(as.numeric(
      passage$predicted_initial_cells
    )))
    if (length(unique_observed_initial) != 1L ||
        length(unique_predicted_initial) != 1L) {
      stop("Growth timepoints have inconsistent initial cell counts for passage ", identities[[1]], ".")
    }
    days <- suppressWarnings(as.numeric(passage$observation_day))
    denominator <- sum(days^2)
    if (!is.finite(denominator) || denominator <= 0) {
      stop("Growth passage ", identities[[1]], " has no positive time support.")
    }
    observed_log_change <- log(
      suppressWarnings(as.numeric(passage$observed_live_cells_at_observation)) /
        unique_observed_initial
    )
    predicted_log_change <- log(
      suppressWarnings(as.numeric(passage$predicted_live_cells_at_observation)) /
        unique_predicted_initial
    )
    observed_rate <- sum(days * observed_log_change) / denominator
    predicted_rate <- sum(days * predicted_log_change) / denominator
    retained_columns <- setdiff(names(passage), point_columns)
    out <- passage[1L, retained_columns, drop = FALSE]
    out$n_growth_timepoints <- nrow(passage)
    out$growth_time_min_day <- min(days)
    out$growth_time_max_day <- max(days)
    out$growth_time_sum_squares <- denominator
    out$observed_passage_growth_rate <- observed_rate
    out$predicted_passage_growth_rate <- predicted_rate
    out$observed_growth <- observed_rate
    out$predicted_growth <- predicted_rate
    out$predicted_growth_rate <- predicted_rate
    out$growth_rate_residual <- observed_rate - predicted_rate
    out
  })
  growth_df <- dplyr::bind_rows(passage_rows)
  growth_df$sigma_growth_rate <- rep(sigma_use, nrow(growth_df))
  growth_df$growth_rate_estimator <-
    "zero_intercept_ols_log_fold_change_on_observation_day"
  growth_df$growth_likelihood_scale <-
    "passage_average_log_growth_rate_per_day"
  growth_df$loglik <- stats::dnorm(
    x = growth_df$observed_passage_growth_rate,
    mean = growth_df$predicted_passage_growth_rate,
    sd = growth_df$sigma_growth_rate,
    log = TRUE
  )
  growth_df
}

ivt_passage_time_loglik_df <- function(summary_df,
                                       passage_time_tolerance_days = 1,
                                       passage_time_sigma_days = 1,
                                       passage_time_df = 4) {
  tolerance_use <- suppressWarnings(as.numeric(passage_time_tolerance_days))
  sigma_use <- suppressWarnings(as.numeric(passage_time_sigma_days))
  df_use <- suppressWarnings(as.numeric(passage_time_df))
  if (length(tolerance_use) != 1L || !is.finite(tolerance_use) || tolerance_use < 0) {
    stop("passage_time_tolerance_days must be one finite non-negative value.")
  }
  if (length(sigma_use) != 1L || !is.finite(sigma_use) || sigma_use <= 0) {
    stop("passage_time_sigma_days must be one finite strictly positive value.")
  }
  if (length(df_use) != 1L || !is.finite(df_use) || df_use <= 0) {
    stop("passage_time_df must be one finite strictly positive value.")
  }
  required <- c(
    "passage_id", "cohort", "lineage_id", "scenario_id",
    "selected_day", "observed_passage_day"
  )
  missing <- setdiff(required, names(summary_df))
  if (length(missing)) {
    stop("Passage-time likelihood summary is missing columns: ", paste(missing, collapse = ", "))
  }
  out <- summary_df[
    is.finite(summary_df$selected_day) & is.finite(summary_df$observed_passage_day),
    ,
    drop = FALSE
  ]
  if (!nrow(out)) {
    out$passage_time_residual_days <- numeric()
    out$passage_time_excess_days <- numeric()
    out$passage_time_loglik <- numeric()
    out$passage_time_nll <- numeric()
    return(out)
  }
  out$passage_time_tolerance_days <- tolerance_use
  out$passage_time_sigma_days <- sigma_use
  out$passage_time_df <- df_use
  out$passage_time_residual_days <- out$selected_day - out$observed_passage_day
  out$passage_time_excess_days <- pmax(
    abs(out$passage_time_residual_days) - tolerance_use,
    0
  )
  standardized_excess <- out$passage_time_excess_days / sigma_use
  out$passage_time_loglik <- stats::dt(standardized_excess, df = df_use, log = TRUE) -
    stats::dt(0, df = df_use, log = TRUE)
  out$passage_time_nll <- -out$passage_time_loglik
  out$loglik <- out$passage_time_loglik
  out
}

ivt_death_loglik_df <- function(run_2N,
                                run_4N,
                                death_data,
                                sigma_death_logit = 0.75,
                                death_fraction_eps = 1e-4,
                                day_tolerance = 1e-8) {
  empty <- data.frame(
    observation_id = character(),
    cohort = character(),
    lineage_id = character(),
    scenario_id = character(),
    segment_id = character(),
    passage_id = character(),
    passage_index = integer(),
    lineage_passage_index = integer(),
    likelihood_observation_day = numeric(),
    model_day = numeric(),
    dead_count = numeric(),
    eligible_denominator = numeric(),
    observed_dead_fraction = numeric(),
    predicted_dead_fraction = numeric(),
    observed_dead_logit = numeric(),
    predicted_dead_logit = numeric(),
    logit_residual = numeric(),
    loglik = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(death_data) || !is.data.frame(death_data) || nrow(death_data) == 0L) {
    return(empty)
  }
  sigma_use <- as.numeric(sigma_death_logit)
  eps_use <- as.numeric(death_fraction_eps)
  tolerance_use <- as.numeric(day_tolerance)
  if (length(sigma_use) != 1L || !is.finite(sigma_use) || sigma_use <= 0) {
    stop("sigma_death_logit must be one finite strictly positive value.")
  }
  if (length(eps_use) != 1L || !is.finite(eps_use) || eps_use <= 0 || eps_use >= 0.5) {
    stop("death_fraction_eps must be one finite value strictly between 0 and 0.5.")
  }
  if (length(tolerance_use) != 1L || !is.finite(tolerance_use) || tolerance_use < 0) {
    stop("day_tolerance must be one finite non-negative value.")
  }

  daily <- dplyr::bind_rows(
    ivt_collect_daily_counts(run_2N),
    ivt_collect_daily_counts(run_4N)
  )
  selected_passages <- dplyr::bind_rows(lapply(
    c(run_2N$segment_results, run_4N$segment_results),
    function(seg_res) {
      data.frame(
        passage_id = as.character(seg_res$segment$passage_id),
        segment_id = as.character(seg_res$segment$segment_id),
        selected_day = as.numeric(seg_res$selection$selected_day),
        stringsAsFactors = FALSE
      )
    }
  ))
  required_daily <- c(
    "cohort", "lineage_id", "scenario_id", "segment_id", "passage_id",
    "passage_index", "lineage_passage_index", "day", "dead_total_cells",
    "total_cells"
  )
  missing_daily <- setdiff(required_daily, names(daily))
  if (length(missing_daily)) {
    stop("In vitro daily predictions are missing death-likelihood columns: ", paste(missing_daily, collapse = ", "))
  }

  matched <- lapply(seq_len(nrow(death_data)), function(i) {
    obs <- death_data[i, , drop = FALSE]
    selected_hit <- which(
      as.character(selected_passages$passage_id) == as.character(obs$model_passage_id[[1]]) &
        as.character(selected_passages$segment_id) == as.character(obs$model_segment_id[[1]])
    )
    if (length(selected_hit) != 1L) {
      stop(
        "Death observation ", obs$observation_id[[1]],
        " must match exactly one selected model passage/segment; found ",
        length(selected_hit), "."
      )
    }
    selected_day <- selected_passages$selected_day[[selected_hit[[1L]]]]
    hits <- which(
      as.character(daily$passage_id) == as.character(obs$model_passage_id[[1]]) &
        as.character(daily$segment_id) == as.character(obs$model_segment_id[[1]]) &
        abs(as.numeric(daily$day) - selected_day) <= tolerance_use
    )
    if (length(hits) != 1L) {
      stop(
        "Death observation ", obs$observation_id[[1]],
        " must match exactly one selected model passage/segment/day; found ", length(hits), "."
      )
    }
    pred <- daily[hits[[1]], , drop = FALSE]
    identity_fields <- c("cohort", "lineage_id", "scenario_id")
    mismatched <- identity_fields[vapply(identity_fields, function(field) {
      !identical(as.character(obs[[field]][[1]]), as.character(pred[[field]][[1]]))
    }, logical(1))]
    if (length(mismatched)) {
      stop(
        "Death observation ", obs$observation_id[[1]],
        " disagrees with the model mapping for: ", paste(mismatched, collapse = ", "), "."
      )
    }
    total_cells <- as.numeric(pred$total_cells[[1]])
    dead_cells <- as.numeric(pred$dead_total_cells[[1]])
    if (!is.finite(total_cells) || total_cells <= 0 ||
        !is.finite(dead_cells) || dead_cells < 0 || dead_cells > total_cells * (1 + 1e-10)) {
      stop("Invalid predicted live/dead cell totals for death observation ", obs$observation_id[[1]], ".")
    }
    data.frame(
      observation_id = as.character(obs$observation_id[[1]]),
      cohort = as.character(pred$cohort[[1]]),
      lineage_id = as.character(pred$lineage_id[[1]]),
      scenario_id = as.character(pred$scenario_id[[1]]),
      segment_id = as.character(pred$segment_id[[1]]),
      passage_id = as.character(pred$passage_id[[1]]),
      passage_index = as.integer(pred$passage_index[[1]]),
      lineage_passage_index = as.integer(pred$lineage_passage_index[[1]]),
      likelihood_observation_day = as.numeric(obs$likelihood_observation_day[[1]]),
      model_day = as.numeric(pred$day[[1]]),
      dead_count = as.numeric(obs$dead_count[[1]]),
      eligible_denominator = as.numeric(obs$eligible_denominator[[1]]),
      observed_dead_fraction = as.numeric(obs$observed_dead_fraction[[1]]),
      predicted_dead_fraction = min(max(dead_cells / total_cells, 0), 1),
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(matched)
  clamp_probability <- function(x) pmin(pmax(as.numeric(x), eps_use), 1 - eps_use)
  out$observed_dead_logit <- stats::qlogis(clamp_probability(out$observed_dead_fraction))
  out$predicted_dead_logit <- stats::qlogis(clamp_probability(out$predicted_dead_fraction))
  out$logit_residual <- out$observed_dead_logit - out$predicted_dead_logit
  out$loglik <- stats::dnorm(
    out$observed_dead_logit,
    mean = out$predicted_dead_logit,
    sd = sigma_use,
    log = TRUE
  )
  out$sigma_death_logit <- sigma_use
  out$death_fraction_eps <- eps_use
  out
}

.ivt_assert_unique_likelihood_units <- function(df,
                                                modality,
                                                unit_id_col = "passage_id") {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(invisible(TRUE))
  if (!unit_id_col %in% names(df)) {
    stop(modality, " likelihood table is missing ", unit_id_col, ".")
  }
  ids <- as.character(df[[unit_id_col]])
  duplicate_ids <- unique(ids[duplicated(ids)])
  if (length(duplicate_ids) > 0L) {
    stop(
      modality, " observations enter the objective more than once: ",
      paste(duplicate_ids, collapse = ", ")
    )
  }
  invisible(TRUE)
}

.ivt_hierarchical_loglik <- function(df, value_col, modality) {
  empty_hierarchy <- data.frame(
    modality = character(),
    aggregation_level = character(),
    cohort = character(),
    lineage_id = character(),
    scenario_id = character(),
    n_units = integer(),
    mean_loglik = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(list(value = 0, lineage = empty_hierarchy, cohort = empty_hierarchy, hierarchy = empty_hierarchy))
  }
  required <- c("cohort", "lineage_id", "scenario_id", value_col)
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop(modality, " likelihood hierarchy is missing columns: ", paste(missing, collapse = ", "))
  }
  work <- data.frame(
    cohort = as.character(df$cohort),
    lineage_id = as.character(df$lineage_id),
    scenario_id = as.character(df$scenario_id),
    unit_loglik = suppressWarnings(as.numeric(df[[value_col]])),
    stringsAsFactors = FALSE
  )
  work <- work[
    is.finite(work$unit_loglik) &
      nzchar(work$cohort) &
      nzchar(work$lineage_id) &
      nzchar(work$scenario_id),
    ,
    drop = FALSE
  ]
  if (nrow(work) == 0L) {
    return(list(value = 0, lineage = empty_hierarchy, cohort = empty_hierarchy, hierarchy = empty_hierarchy))
  }

  lineage_mean <- stats::aggregate(
    unit_loglik ~ cohort + lineage_id + scenario_id,
    data = work,
    FUN = mean
  )
  lineage_n <- stats::aggregate(
    unit_loglik ~ cohort + lineage_id + scenario_id,
    data = work,
    FUN = length
  )
  names(lineage_mean)[names(lineage_mean) == "unit_loglik"] <- "mean_loglik"
  names(lineage_n)[names(lineage_n) == "unit_loglik"] <- "n_units"
  lineage <- merge(
    lineage_mean,
    lineage_n,
    by = c("cohort", "lineage_id", "scenario_id"),
    sort = FALSE
  )
  lineage$modality <- modality
  lineage$aggregation_level <- "lineage"
  lineage <- lineage[, names(empty_hierarchy), drop = FALSE]

  cohort_mean <- stats::aggregate(
    mean_loglik ~ cohort,
    data = lineage,
    FUN = mean
  )
  cohort_n <- stats::aggregate(
    scenario_id ~ cohort,
    data = lineage,
    FUN = length
  )
  names(cohort_n)[names(cohort_n) == "scenario_id"] <- "n_units"
  cohort_df <- merge(cohort_mean, cohort_n, by = "cohort", sort = FALSE)
  cohort_df$modality <- modality
  cohort_df$aggregation_level <- "cohort"
  cohort_df$lineage_id <- NA_character_
  cohort_df$scenario_id <- NA_character_
  cohort_df <- cohort_df[, names(empty_hierarchy), drop = FALSE]

  value <- mean(cohort_df$mean_loglik)
  global_df <- data.frame(
    modality = modality,
    aggregation_level = "global",
    cohort = NA_character_,
    lineage_id = NA_character_,
    scenario_id = NA_character_,
    n_units = nrow(cohort_df),
    mean_loglik = value,
    stringsAsFactors = FALSE
  )
  list(
    value = value,
    lineage = lineage,
    cohort = cohort_df,
    hierarchy = dplyr::bind_rows(lineage, cohort_df, global_df)
  )
}

ivt_smooth_kary_distribution <- function(grid, probs, sigma_kary) {
  sigma_use <- as.numeric(sigma_kary)
  if (!is.finite(sigma_use) || sigma_use <= 0) {
    stop("sigma_kary must be strictly positive.")
  }

  grid <- as.numeric(grid)
  probs <- as.numeric(probs)
  probs <- pmax(probs, 0)
  prob_sum <- sum(probs)
  if (!is.finite(prob_sum) || prob_sum <= 0) {
    probs <- rep(1 / length(grid), length(grid))
  } else {
    probs <- probs / prob_sum
  }

  kernel_mat <- outer(
    X = grid,
    Y = grid,
    FUN = function(obs_n, true_n) stats::dnorm(obs_n, mean = true_n, sd = sigma_use)
  )
  kernel_colsum <- colSums(kernel_mat)
  kernel_colsum[!is.finite(kernel_colsum) | kernel_colsum <= 0] <- 1
  kernel_mat <- sweep(kernel_mat, 2, kernel_colsum, "/", check.margin = FALSE)

  as.numeric(kernel_mat %*% probs)
}

ivt_default_flow_kernel_sd_ploidy <- function(run, fit_data) {
  observed_grids <- unlist(lapply(fit_data, function(entry) {
    obs <- ivt_observed_passage_summary(entry)$observed_flow
    if (is.null(obs) || nrow(obs) < 2L) return(numeric(0))
    unique(as.numeric(obs$ploidy))
  }), use.names = FALSE)
  observed_grids <- sort(unique(observed_grids[is.finite(observed_grids)]))
  step <- if (length(observed_grids) >= 2L) stats::median(diff(observed_grids)) else NA_real_

  default_sd <- max(
    if (is.finite(step) && step > 0) 2 * step else 0,
    1 / 22
  )
  as.numeric(default_sd)
}

ivt_predicted_flow_density <- function(observed_grid,
                                       run_grid_pre,
                                       selected_frac,
                                       sigma_flow_ploidy,
                                       n_unit) {
  observed_grid <- as.numeric(observed_grid)
  model_grid <- as.numeric(run_grid_pre) / as.numeric(n_unit)
  probs <- as.numeric(selected_frac)
  probs <- pmax(probs, 0)
  prob_sum <- sum(probs)
  if (!is.finite(prob_sum) || prob_sum <= 0) {
    probs <- rep(1 / length(model_grid), length(model_grid))
  } else {
    probs <- probs / prob_sum
  }

  kernel_sd <- as.numeric(sigma_flow_ploidy)
  if (!is.finite(kernel_sd) || kernel_sd <= 0) {
    stop("sigma_flow_ploidy must be strictly positive.")
  }

  density_vals <- as.numeric(vapply(observed_grid, function(obs_p) {
    sum(probs * stats::dnorm(obs_p, mean = model_grid, sd = kernel_sd))
  }, numeric(1)))

  if (length(observed_grid) >= 2L) {
    step <- stats::median(diff(sort(unique(observed_grid))))
  } else {
    step <- 1
  }
  norm_const <- sum(density_vals) * step
  if (!is.finite(norm_const) || norm_const <= 0) {
    density_vals <- rep(1 / (length(observed_grid) * step), length(observed_grid))
  } else {
    density_vals <- density_vals / norm_const
  }

  density_vals
}

ivt_ploidy_loglik_df <- function(run, fit_data, sigma_kary, prob_floor = 1e-12) {
  prob_floor <- as.numeric(prob_floor)
  if (!is.finite(prob_floor) || prob_floor <= 0 || prob_floor >= 1) {
    stop("prob_floor must lie in (0, 1).")
  }

  out <- lapply(run$segment_results, function(seg_res) {
    seg <- seg_res$segment
    probs <- as.numeric(seg_res$selection$selected_frac)
    probs <- pmax(probs, 0)
    prob_sum <- sum(probs)
    if (!is.finite(prob_sum) || prob_sum <= 0) {
      probs <- rep(1 / length(run$grid_pre), length(run$grid_pre))
    } else {
      probs <- probs / prob_sum
    }
    probs <- ivt_smooth_kary_distribution(
      grid = run$grid_pre,
      probs = probs,
      sigma_kary = sigma_kary
    )
    names(probs) <- as.character(run$grid_pre)

    do.call(rbind, lapply(.ivt_segment_endpoint_data_ids(run, seg), function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_kary <- obs$observed_kary
      observed_kary <- observed_kary[is.finite(observed_kary)]
      if (length(observed_kary) == 0L) {
        return(data.frame(
          segment_id = seg$segment_id,
          cohort = seg$cohort,
          lineage_id = seg$lineage_id,
          lineage_group = seg$lineage_group,
          lineage_label = seg$lineage_label,
          scenario_id = seg$scenario_id,
          lineage_terminal_key = seg$lineage_terminal_key,
          passage_index = seg$passage_index,
          lineage_passage_index = seg$lineage_passage_index,
          oxygen_pct = seg$oxygen_pct,
          passage_id = pid,
          n_cells = 0L,
          mean_loglik = NA_real_,
          total_loglik = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      prob_lookup <- probs[as.character(observed_kary)]
      prob_lookup[is.na(prob_lookup)] <- prob_floor
      prob_lookup <- pmax(prob_lookup, prob_floor)
      cell_loglik <- log(prob_lookup)

      data.frame(
        segment_id = seg$segment_id,
        cohort = seg$cohort,
        lineage_id = seg$lineage_id,
        lineage_group = seg$lineage_group,
        lineage_label = seg$lineage_label,
        scenario_id = seg$scenario_id,
        lineage_terminal_key = seg$lineage_terminal_key,
        passage_index = seg$passage_index,
        lineage_passage_index = seg$lineage_passage_index,
        oxygen_pct = seg$oxygen_pct,
        passage_id = pid,
        n_cells = length(observed_kary),
        mean_loglik = mean(cell_loglik),
        total_loglik = sum(cell_loglik),
        stringsAsFactors = FALSE
      )
    }))
  })

  initial_probs <- as.numeric(run$initial_fraction)
  initial_probs <- pmax(initial_probs, 0)
  initial_prob_sum <- sum(initial_probs)
  if (!is.finite(initial_prob_sum) || initial_prob_sum <= 0) {
    initial_probs <- rep(1 / length(run$grid_pre), length(run$grid_pre))
  } else {
    initial_probs <- initial_probs / initial_prob_sum
  }
  initial_probs <- ivt_smooth_kary_distribution(
    grid = run$grid_pre,
    probs = initial_probs,
    sigma_kary = sigma_kary
  )
  names(initial_probs) <- as.character(run$grid_pre)
  initial_out <- lapply(run$initial_observations, function(record) {
    obs <- record$observed
    observed_kary <- obs$observed_kary
    observed_kary <- observed_kary[is.finite(observed_kary)]
    if (length(observed_kary) == 0L) return(NULL)
    prob_lookup <- initial_probs[as.character(observed_kary)]
    prob_lookup[is.na(prob_lookup)] <- prob_floor
    cell_loglik <- log(pmax(prob_lookup, prob_floor))
    data.frame(
      segment_id = record$segment_id,
      cohort = record$cohort,
      lineage_id = record$lineage_id,
      lineage_group = record$lineage_group,
      lineage_label = record$lineage_id,
      scenario_id = record$scenario_id,
      lineage_terminal_key = record$scenario_id,
      passage_index = 0L,
      lineage_passage_index = 0L,
      oxygen_pct = record$oxygen_pct,
      passage_id = record$passage_id,
      n_cells = length(observed_kary),
      mean_loglik = mean(cell_loglik),
      total_loglik = sum(cell_loglik),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(c(out, initial_out))
}

ivt_flow_loglik_df <- function(run,
                               fit_data,
                               n_unit,
                               sigma_flow_ploidy,
                               density_floor = 1e-12) {
  density_floor <- as.numeric(density_floor)
  if (!is.finite(density_floor) || density_floor <= 0) {
    stop("density_floor must be strictly positive.")
  }

  out <- lapply(run$segment_results, function(seg_res) {
    seg <- seg_res$segment
    do.call(rbind, lapply(.ivt_segment_endpoint_data_ids(run, seg), function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_flow <- obs$observed_flow
      if (is.null(observed_flow) || nrow(observed_flow) == 0L) {
        return(data.frame(
          segment_id = seg$segment_id,
          cohort = seg$cohort,
          lineage_id = seg$lineage_id,
          lineage_group = seg$lineage_group,
          lineage_label = seg$lineage_label,
          scenario_id = seg$scenario_id,
          lineage_terminal_key = seg$lineage_terminal_key,
          passage_index = seg$passage_index,
          lineage_passage_index = seg$lineage_passage_index,
          oxygen_pct = seg$oxygen_pct,
          passage_id = pid,
          sample_name = NA_character_,
          n_grid = 0L,
          mean_loglik = NA_real_,
          total_loglik = NA_real_,
          sigma_flow_ploidy = as.numeric(sigma_flow_ploidy),
          stringsAsFactors = FALSE
        ))
      }

      obs_grid <- as.numeric(observed_flow$ploidy)
      obs_density <- as.numeric(observed_flow$density)
      keep <- is.finite(obs_grid) & is.finite(obs_density) & obs_density >= 0
      obs_grid <- obs_grid[keep]
      obs_density <- obs_density[keep]
      if (!length(obs_grid)) {
        return(data.frame(
          segment_id = seg$segment_id,
          cohort = seg$cohort,
          lineage_id = seg$lineage_id,
          lineage_group = seg$lineage_group,
          lineage_label = seg$lineage_label,
          scenario_id = seg$scenario_id,
          lineage_terminal_key = seg$lineage_terminal_key,
          passage_index = seg$passage_index,
          lineage_passage_index = seg$lineage_passage_index,
          oxygen_pct = seg$oxygen_pct,
          passage_id = pid,
          sample_name = obs$observed_flow_sample_name,
          n_grid = 0L,
          mean_loglik = NA_real_,
          total_loglik = NA_real_,
          sigma_flow_ploidy = as.numeric(sigma_flow_ploidy),
          stringsAsFactors = FALSE
        ))
      }

      if (length(obs_grid) >= 2L) {
        step <- stats::median(diff(sort(unique(obs_grid))))
      } else {
        step <- 1
      }
      obs_mass <- obs_density * step
      obs_mass_sum <- sum(obs_mass)
      if (!is.finite(obs_mass_sum) || obs_mass_sum <= 0) {
        obs_mass <- rep(1 / length(obs_grid), length(obs_grid))
      } else {
        obs_mass <- obs_mass / obs_mass_sum
      }

      pred_density <- ivt_predicted_flow_density(
        observed_grid = obs_grid,
        run_grid_pre = run$grid_pre,
        selected_frac = seg_res$selection$selected_frac,
        sigma_flow_ploidy = sigma_flow_ploidy,
        n_unit = n_unit
      )
      pred_density <- pmax(pred_density, density_floor)
      point_loglik <- log(pred_density)
      sample_loglik <- sum(obs_mass * point_loglik)

      data.frame(
        segment_id = seg$segment_id,
        cohort = seg$cohort,
        lineage_id = seg$lineage_id,
        lineage_group = seg$lineage_group,
        lineage_label = seg$lineage_label,
        scenario_id = seg$scenario_id,
        lineage_terminal_key = seg$lineage_terminal_key,
        passage_index = seg$passage_index,
        lineage_passage_index = seg$lineage_passage_index,
        oxygen_pct = seg$oxygen_pct,
        passage_id = pid,
        sample_name = obs$observed_flow_sample_name,
        n_grid = length(obs_grid),
        mean_loglik = sample_loglik,
        total_loglik = sample_loglik,
        sigma_flow_ploidy = as.numeric(sigma_flow_ploidy),
        stringsAsFactors = FALSE
      )
    }))
  })

  dplyr::bind_rows(out)
}

ivt_flow_overlay_df <- function(run,
                                fit_data,
                                n_unit,
                                sigma_flow_ploidy) {
  out <- lapply(run$segment_results, function(seg_res) {
    seg <- seg_res$segment
    dplyr::bind_rows(lapply(.ivt_segment_endpoint_data_ids(run, seg), function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_flow <- obs$observed_flow
      if (is.null(observed_flow) || nrow(observed_flow) == 0L) return(NULL)

      obs_grid <- as.numeric(observed_flow$ploidy)
      pred_density <- ivt_predicted_flow_density(
        observed_grid = obs_grid,
        run_grid_pre = run$grid_pre,
        selected_frac = seg_res$selection$selected_frac,
        sigma_flow_ploidy = sigma_flow_ploidy,
        n_unit = n_unit
      )

      dplyr::bind_rows(
        data.frame(
          segment_id = seg$segment_id,
          cohort = seg$cohort,
          lineage_id = seg$lineage_id,
          lineage_group = seg$lineage_group,
          lineage_label = seg$lineage_label,
          scenario_id = seg$scenario_id,
          lineage_terminal_key = seg$lineage_terminal_key,
          passage_index = seg$passage_index,
          lineage_passage_index = seg$lineage_passage_index,
          oxygen_pct = seg$oxygen_pct,
          passage_id = pid,
          sample_name = obs$observed_flow_sample_name,
          grid_index = observed_flow$grid_index,
          ploidy = obs_grid,
          density = as.numeric(observed_flow$density),
          series = "Observed",
          stringsAsFactors = FALSE
        ),
        data.frame(
          segment_id = seg$segment_id,
          cohort = seg$cohort,
          lineage_id = seg$lineage_id,
          lineage_group = seg$lineage_group,
          lineage_label = seg$lineage_label,
          scenario_id = seg$scenario_id,
          lineage_terminal_key = seg$lineage_terminal_key,
          passage_index = seg$passage_index,
          lineage_passage_index = seg$lineage_passage_index,
          oxygen_pct = seg$oxygen_pct,
          passage_id = pid,
          sample_name = obs$observed_flow_sample_name,
          grid_index = observed_flow$grid_index,
          ploidy = obs_grid,
          density = pred_density,
          series = "Predicted",
          stringsAsFactors = FALSE
        )
      )
    }))
  })

  dplyr::bind_rows(out)
}

ivt_objective_components <- function(run_params,
                                     fit_objects,
                                     cfg,
                                     fallback_max_passage_days = 14,
                                     growth_weight = 1,
                                     ploidy_weight = 1,
                                     flow_weight = 1,
                                     death_weight = 1,
                                     passage_time_weight = 0.25,
                                     passage_time_tolerance_days = 1,
                                     passage_time_sigma_days = 1,
                                     passage_time_df = 4,
                                     sigma_death_logit = 0.75,
                                     death_fraction_eps = 1e-4,
                                     ploidy_prob_floor = 1e-12,
                                     flow_density_floor = 1e-12,
                                     flow_kernel_sd_ploidy = NULL) {
  model_core <- build_model_core(cfg = cfg)

  adapter_2N <- ivt_build_job_table_adapter(
    jobs = fit_objects$jobs_2N,
    fit_data = fit_objects$fit_data,
    cohort = "2N",
    fallback_max_passage_days = fallback_max_passage_days,
    passage_time_tolerance_days = passage_time_tolerance_days
  )
  adapter_4N <- ivt_build_job_table_adapter(
    jobs = fit_objects$jobs_4N,
    fit_data = fit_objects$fit_data,
    cohort = "4N",
    fallback_max_passage_days = fallback_max_passage_days,
    passage_time_tolerance_days = passage_time_tolerance_days
  )

  run_2N <- ivt_run_lineage(adapter_2N, cfg = cfg, run_params = run_params, model_core = model_core)
  run_4N <- ivt_run_lineage(adapter_4N, cfg = cfg, run_params = run_params, model_core = model_core)

  sum_2N <- ivt_collect_lineage_summary(run_2N, fit_objects$fit_data)
  sum_4N <- ivt_collect_lineage_summary(run_4N, fit_objects$fit_data)
  summary_df <- dplyr::bind_rows(sum_2N, sum_4N)
  n_insufficient_boundaries <- sum(
    summary_df$reseed_mode == "carry_forward_insufficient",
    na.rm = TRUE
  )
  all_passage_boundaries_feasible <- n_insufficient_boundaries == 0L
  protocol_feasibility_status <- if (
    all_passage_boundaries_feasible
  ) {
    "PASS"
  } else {
    "FAIL"
  }
  summary_df$n_insufficient_boundaries <- n_insufficient_boundaries
  summary_df$all_passage_boundaries_feasible <-
    all_passage_boundaries_feasible
  summary_df$protocol_feasibility_status <-
    protocol_feasibility_status
  .ivt_assert_unique_likelihood_units(summary_df, "growth/passage summary")

  sigma_growth <- as.numeric(run_params$sigma_growth)
  sigma_kary <- as.numeric(run_params$sigma_kary)
  sigma_flow_ploidy <- if (is.null(flow_kernel_sd_ploidy)) {
    ivt_default_flow_kernel_sd_ploidy(run = run_2N, fit_data = fit_objects$fit_data)
  } else {
    as.numeric(flow_kernel_sd_ploidy)
  }
  daily_counts <- dplyr::bind_rows(
    ivt_collect_daily_counts(run_2N),
    ivt_collect_daily_counts(run_4N)
  )
  growth_count_diagnostics_df <- .ivt_growth_measurement_summary(
    summary_df = summary_df,
    fit_data = fit_objects$fit_data,
    daily_counts = daily_counts
  )
  growth_df <- ivt_growth_loglik_df(
    summary_df = growth_count_diagnostics_df,
    sigma_growth = sigma_growth
  )
  growth_passage_df <- .ivt_growth_passage_loglik_df(growth_df)
  rate_index <- match(
    as.character(summary_df$passage_id),
    as.character(growth_df$passage_id)
  )
  if (anyNA(rate_index)) {
    stop("Every simulated passage must have one passage-average growth-rate likelihood row.")
  }
  summary_df$observed_growth <-
    growth_df$observed_passage_growth_rate[rate_index]
  summary_df$predicted_growth <-
    growth_df$predicted_passage_growth_rate[rate_index]
  summary_df$predicted_growth_rate <-
    growth_df$predicted_passage_growth_rate[rate_index]
  passage_time_df_table <- ivt_passage_time_loglik_df(
    summary_df = summary_df,
    passage_time_tolerance_days = passage_time_tolerance_days,
    passage_time_sigma_days = passage_time_sigma_days,
    passage_time_df = passage_time_df
  )
  observed_growth_rates <- suppressWarnings(as.numeric(
    growth_df$observed_passage_growth_rate
  ))
  predicted_growth_rates <- suppressWarnings(as.numeric(
    growth_df$predicted_passage_growth_rate
  ))
  n_growth_observed <- sum(is.finite(observed_growth_rates))
  n_growth_missing_pred <- sum(
    is.finite(observed_growth_rates) & !is.finite(predicted_growth_rates)
  )
  n_growth_negative_pred <- sum(
    is.finite(predicted_growth_rates) & predicted_growth_rates < 0
  )
  ploidy_df <- dplyr::bind_rows(
    ivt_ploidy_loglik_df(run = run_2N, fit_data = fit_objects$fit_data, sigma_kary = sigma_kary, prob_floor = ploidy_prob_floor),
    ivt_ploidy_loglik_df(run = run_4N, fit_data = fit_objects$fit_data, sigma_kary = sigma_kary, prob_floor = ploidy_prob_floor)
  )
  ploidy_df <- ploidy_df[is.finite(ploidy_df$mean_loglik), , drop = FALSE]
  flow_df <- dplyr::bind_rows(
    ivt_flow_loglik_df(run = run_2N, fit_data = fit_objects$fit_data, n_unit = cfg$N_UNIT, sigma_flow_ploidy = sigma_flow_ploidy, density_floor = flow_density_floor),
    ivt_flow_loglik_df(run = run_4N, fit_data = fit_objects$fit_data, n_unit = cfg$N_UNIT, sigma_flow_ploidy = sigma_flow_ploidy, density_floor = flow_density_floor)
  )
  flow_df <- flow_df[is.finite(flow_df$mean_loglik), , drop = FALSE]
  flow_overlay_df <- dplyr::bind_rows(
    ivt_flow_overlay_df(run = run_2N, fit_data = fit_objects$fit_data, n_unit = cfg$N_UNIT, sigma_flow_ploidy = sigma_flow_ploidy),
    ivt_flow_overlay_df(run = run_4N, fit_data = fit_objects$fit_data, n_unit = cfg$N_UNIT, sigma_flow_ploidy = sigma_flow_ploidy)
  )
  death_weight_use <- as.numeric(death_weight)
  if (length(death_weight_use) != 1L || !is.finite(death_weight_use) || death_weight_use < 0) {
    stop("death_weight must be one finite non-negative value.")
  }
  death_df <- if (death_weight_use == 0) {
    ivt_death_loglik_df(run_2N, run_4N, data.frame())
  } else {
    if (is.null(fit_objects$death_data) || !is.data.frame(fit_objects$death_data) ||
        nrow(fit_objects$death_data) == 0L) {
      stop("death_weight is positive but fit_objects$death_data has no enabled observations.")
    }
    ivt_death_loglik_df(
      run_2N = run_2N,
      run_4N = run_4N,
      death_data = fit_objects$death_data,
      sigma_death_logit = sigma_death_logit,
      death_fraction_eps = death_fraction_eps
    )
  }

  .ivt_assert_unique_likelihood_units(
    growth_count_diagnostics_df,
    "growth count diagnostic",
    unit_id_col = "growth_observation_id"
  )
  .ivt_assert_unique_likelihood_units(growth_df, "growth")
  .ivt_assert_unique_likelihood_units(growth_passage_df, "growth passage")
  .ivt_assert_unique_likelihood_units(ploidy_df, "karyotype")
  .ivt_assert_unique_likelihood_units(flow_df, "flow")
  .ivt_assert_unique_likelihood_units(death_df, "death")
  .ivt_assert_unique_likelihood_units(passage_time_df_table, "passage_time")
  growth_aggregation <- .ivt_hierarchical_loglik(
    growth_passage_df,
    value_col = "mean_loglik",
    modality = "growth"
  )
  ploidy_aggregation <- .ivt_hierarchical_loglik(
    ploidy_df,
    value_col = "mean_loglik",
    modality = "karyotype"
  )
  flow_aggregation <- .ivt_hierarchical_loglik(
    flow_df,
    value_col = "mean_loglik",
    modality = "flow"
  )
  death_aggregation <- .ivt_hierarchical_loglik(
    death_df,
    value_col = "loglik",
    modality = "death"
  )
  passage_time_aggregation <- .ivt_hierarchical_loglik(
    passage_time_df_table,
    value_col = "passage_time_loglik",
    modality = "passage_time"
  )

  growth_loglik_sum <- if (nrow(growth_df) > 0L) sum(growth_df$loglik) else 0
  growth_passage_loglik_sum <- if (nrow(growth_passage_df) > 0L) {
    sum(growth_passage_df$mean_loglik)
  } else {
    0
  }
  ploidy_loglik_sum <- if (nrow(ploidy_df) > 0L) sum(ploidy_df$mean_loglik) else 0
  flow_loglik_sum <- if (nrow(flow_df) > 0L) sum(flow_df$mean_loglik) else 0
  death_loglik_sum <- if (nrow(death_df) > 0L) sum(death_df$loglik) else 0
  passage_time_loglik_sum <- if (nrow(passage_time_df_table) > 0L) {
    sum(passage_time_df_table$passage_time_loglik)
  } else {
    0
  }
  growth_loglik <- growth_aggregation$value
  ploidy_loglik <- ploidy_aggregation$value
  flow_loglik <- flow_aggregation$value
  death_loglik <- death_aggregation$value
  passage_time_loglik <- passage_time_aggregation$value
  total_loglik <- as.numeric(growth_weight) * growth_loglik +
    as.numeric(ploidy_weight) * ploidy_loglik +
    as.numeric(flow_weight) * flow_loglik
  if (death_weight_use != 0) {
    total_loglik <- total_loglik + death_weight_use * death_loglik
  }
  passage_time_weight_use <- as.numeric(passage_time_weight)
  if (length(passage_time_weight_use) != 1L ||
      !is.finite(passage_time_weight_use) || passage_time_weight_use < 0) {
    stop("passage_time_weight must be one finite non-negative value.")
  }
  if (passage_time_weight_use != 0) {
    total_loglik <- total_loglik + passage_time_weight_use * passage_time_loglik
  }
  total <- -total_loglik
  finite_boundary_scales <- suppressWarnings(as.numeric(summary_df$boundary_scale))
  finite_boundary_scales <- finite_boundary_scales[is.finite(finite_boundary_scales)]
  list(
    objective = total,
    total_loglik = total_loglik,
    growth_loglik = growth_loglik,
    ploidy_loglik = ploidy_loglik,
    flow_loglik = flow_loglik,
    death_loglik = death_loglik,
    passage_time_loglik = passage_time_loglik,
    growth_loglik_sum = growth_loglik_sum,
    growth_passage_loglik_sum = growth_passage_loglik_sum,
    ploidy_loglik_sum = ploidy_loglik_sum,
    flow_loglik_sum = flow_loglik_sum,
    death_loglik_sum = death_loglik_sum,
    passage_time_loglik_sum = passage_time_loglik_sum,
    sigma_growth = sigma_growth,
    sigma_kary = sigma_kary,
    sigma_flow_ploidy = sigma_flow_ploidy,
    sigma_death_logit = as.numeric(sigma_death_logit),
    death_fraction_eps = as.numeric(death_fraction_eps),
    death_weight = death_weight_use,
    passage_time_weight = passage_time_weight_use,
    passage_time_tolerance_days = as.numeric(passage_time_tolerance_days),
    passage_time_sigma_days = as.numeric(passage_time_sigma_days),
    passage_time_df_value = as.numeric(passage_time_df),
    n_growth = nrow(growth_df),
    n_growth_timepoints = nrow(growth_count_diagnostics_df),
    n_growth_passages = nrow(growth_passage_df),
    n_growth_observed = n_growth_observed,
    n_growth_missing_pred = n_growth_missing_pred,
    n_growth_negative_pred = n_growth_negative_pred,
    n_ploidy_passages = nrow(ploidy_df),
    n_kary_cells = sum(ploidy_df$n_cells),
    n_flow_passages = nrow(flow_df),
    n_flow_samples = nrow(flow_df),
    n_death_passages = nrow(death_df),
    n_passage_time_observations = nrow(passage_time_df_table),
    n_scenarios = length(unique(summary_df$scenario_id)),
    n_insufficient_boundaries = n_insufficient_boundaries,
    all_passage_boundaries_feasible =
      all_passage_boundaries_feasible,
    protocol_feasibility_status = protocol_feasibility_status,
    max_boundary_scale = if (length(finite_boundary_scales)) {
      max(finite_boundary_scales)
    } else {
      NA_real_
    },
    summary = summary_df,
    growth_df = growth_df,
    growth_count_diagnostics_df = growth_count_diagnostics_df,
    growth_passage_df = growth_passage_df,
    ploidy_df = ploidy_df,
    flow_df = flow_df,
    death_df = death_df,
    passage_time_df = passage_time_df_table,
    flow_overlay_df = flow_overlay_df,
    objective_hierarchy = dplyr::bind_rows(
      growth_aggregation$hierarchy,
      ploidy_aggregation$hierarchy,
      flow_aggregation$hierarchy,
      death_aggregation$hierarchy,
      passage_time_aggregation$hierarchy
    ),
    growth_lineage_loglik = growth_aggregation$lineage,
    growth_cohort_loglik = growth_aggregation$cohort,
    ploidy_lineage_loglik = ploidy_aggregation$lineage,
    ploidy_cohort_loglik = ploidy_aggregation$cohort,
    flow_lineage_loglik = flow_aggregation$lineage,
    flow_cohort_loglik = flow_aggregation$cohort,
    death_lineage_loglik = death_aggregation$lineage,
    death_cohort_loglik = death_aggregation$cohort,
    passage_time_lineage_loglik = passage_time_aggregation$lineage,
    passage_time_cohort_loglik = passage_time_aggregation$cohort,
    run_2N = run_2N,
    run_4N = run_4N
  )
}

ivt_objective <- function(run_params,
                          fit_objects,
                          cfg,
                          fallback_max_passage_days = 14,
                          growth_weight = 1,
                          ploidy_weight = 1,
                          flow_weight = 1,
                          death_weight = 1,
                          passage_time_weight = 0.25,
                          passage_time_tolerance_days = 1,
                          passage_time_sigma_days = 1,
                          passage_time_df = 4,
                          sigma_death_logit = 0.75,
                          death_fraction_eps = 1e-4) {
  ivt_objective_components(
    run_params = run_params,
    fit_objects = fit_objects,
    cfg = cfg,
    fallback_max_passage_days = fallback_max_passage_days,
    growth_weight = growth_weight,
    ploidy_weight = ploidy_weight,
    flow_weight = flow_weight,
    death_weight = death_weight,
    passage_time_weight = passage_time_weight,
    passage_time_tolerance_days = passage_time_tolerance_days,
    passage_time_sigma_days = passage_time_sigma_days,
    passage_time_df = passage_time_df,
    sigma_death_logit = sigma_death_logit,
    death_fraction_eps = death_fraction_eps
  )$objective
}
