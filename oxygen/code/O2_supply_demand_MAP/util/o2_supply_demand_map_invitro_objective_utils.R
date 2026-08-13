# Canonical in-vitro objective, likelihood, and parameter-transform helpers.

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
    stringsAsFactors = FALSE
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

ivt_growth_loglik_df <- function(summary_df, sigma_growth) {
  sigma_use <- as.numeric(sigma_growth)
  if (!is.finite(sigma_use) || sigma_use <= 0) {
    stop("sigma_growth must be strictly positive.")
  }

  growth_df <- summary_df[is.finite(summary_df$observed_growth) & is.finite(summary_df$predicted_growth_rate), , drop = FALSE]
  if (nrow(growth_df) == 0L) {
    growth_df$loglik <- numeric(0)
    return(growth_df)
  }

  growth_df$loglik <- stats::dnorm(
    x = growth_df$observed_growth,
    mean = growth_df$predicted_growth_rate,
    sd = sigma_use,
    log = TRUE
  )
  growth_df
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

    do.call(rbind, lapply(seg$data_ids, function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_kary <- obs$observed_kary
      observed_kary <- observed_kary[is.finite(observed_kary)]
      if (length(observed_kary) == 0L) {
        return(data.frame(
          segment_id = seg$segment_id,
          cohort = seg$cohort,
          passage_index = seg$passage_index,
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
        passage_index = seg$passage_index,
        oxygen_pct = seg$oxygen_pct,
        passage_id = pid,
        n_cells = length(observed_kary),
        mean_loglik = mean(cell_loglik),
        total_loglik = sum(cell_loglik),
        stringsAsFactors = FALSE
      )
    }))
  })

  dplyr::bind_rows(out)
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
    do.call(rbind, lapply(seg$data_ids, function(pid) {
      obs <- ivt_observed_passage_summary(fit_data[[pid]])
      observed_flow <- obs$observed_flow
      if (is.null(observed_flow) || nrow(observed_flow) == 0L) {
        return(data.frame(
          segment_id = seg$segment_id,
          cohort = seg$cohort,
          passage_index = seg$passage_index,
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
          passage_index = seg$passage_index,
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
        passage_index = seg$passage_index,
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
    dplyr::bind_rows(lapply(seg$data_ids, function(pid) {
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
          passage_index = seg$passage_index,
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
          passage_index = seg$passage_index,
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

ivt_objective_components <- local({
  ivt_empty_death_loglik_df <- function() {
    data.frame(
      observation_id = character(),
      model_passage_id = character(),
      segment_id = character(),
      matched_old_segment_id = character(),
      cohort = character(),
      lineage_id = character(),
      selected_index = integer(),
      selected_day = numeric(),
      observed_dead_count = numeric(),
      observed_live_count = numeric(),
      eligible_denominator = numeric(),
      predicted_dead_hypoxia_count = numeric(),
      predicted_dead_buffer_count = numeric(),
      predicted_dead_count = numeric(),
      predicted_live_count = numeric(),
      observed_dead_fraction = numeric(),
      predicted_dead_fraction = numeric(),
      observed_dead_logit = numeric(),
      predicted_dead_logit = numeric(),
      logit_residual = numeric(),
      sigma_death_logit = numeric(),
      death_fraction_eps = numeric(),
      loglik = numeric(),
      stringsAsFactors = FALSE
    )
  }

  ivt_death_loglik_df <- function(run_2N,
                                  run_4N,
                                  death_data,
                                  sigma_death_logit = 0.75,
                                  death_fraction_eps = 1e-4) {
    if (is.null(death_data) || !is.data.frame(death_data) || nrow(death_data) == 0L) {
      return(ivt_empty_death_loglik_df())
    }
    if (nrow(death_data) != 90L) {
      stop("Death likelihood requires exactly 90 validated observation rows.")
    }
    required_cols <- c(
      "observation_id", "model_passage_id", "cohort", "lineage_id",
      "dead_count", "observed_live_count", "eligible_denominator",
      "observed_dead_fraction"
    )
    missing_cols <- setdiff(required_cols, names(death_data))
    if (length(missing_cols) > 0L) {
      stop("Validated Death data are missing: ", paste(missing_cols, collapse = ", "))
    }
    sigma_use <- as.numeric(sigma_death_logit)
    eps_use <- as.numeric(death_fraction_eps)
    if (length(sigma_use) != 1L || !is.finite(sigma_use) || sigma_use <= 0) {
      stop("sigma_death_logit must be one finite strictly positive value.")
    }
    if (length(eps_use) != 1L || !is.finite(eps_use) || eps_use <= 0 || eps_use >= 0.5) {
      stop("death_fraction_eps must be one finite value strictly between 0 and 0.5.")
    }

    collect_run_map <- function(run) {
      rows <- lapply(run$segment_results, function(seg_res) {
        seg <- seg_res$segment
        passage_ids <- as.character(seg$data_ids)
        if (!length(passage_ids)) return(NULL)

        selected_index <- as.integer(seg_res$selection$selected_index)
        selected_day <- as.numeric(seg_res$selection$selected_day)
        live_cells <- as.numeric(seg_res$sim$Ntot_live_obs)
        dead_hypoxia <- as.numeric(seg_res$sim$Ntot_dead_hypoxia_obs)
        dead_buffer <- as.numeric(seg_res$sim$Ntot_dead_buffer_obs)
        trajectory_lengths <- c(length(live_cells), length(dead_hypoxia), length(dead_buffer))
        if (length(unique(trajectory_lengths)) != 1L || trajectory_lengths[[1]] == 0L) {
          stop(
            "Live/dead stock trajectories are missing or have incompatible lengths for segment ",
            seg$segment_id, "."
          )
        }
        if (length(selected_index) != 1L || is.na(selected_index) ||
            selected_index < 1L || selected_index > trajectory_lengths[[1]]) {
          stop("Invalid selected_index for Death likelihood segment ", seg$segment_id, ".")
        }
        if (length(selected_day) != 1L || !is.finite(selected_day)) {
          stop("Invalid selected_day for Death likelihood segment ", seg$segment_id, ".")
        }
        if (length(seg$obs_days_local) >= selected_index &&
            !isTRUE(all.equal(
              selected_day,
              as.numeric(seg$obs_days_local[[selected_index]]),
              tolerance = 1e-12
            ))) {
          stop("selected_index and selected_day disagree for Death likelihood segment ", seg$segment_id, ".")
        }

        predicted_live <- live_cells[[selected_index]]
        predicted_dead_hypoxia <- dead_hypoxia[[selected_index]]
        predicted_dead_buffer <- dead_buffer[[selected_index]]
        predicted_dead <- predicted_dead_hypoxia + predicted_dead_buffer
        selected_live <- as.numeric(seg_res$selection$selected_live_cells)
        if (!isTRUE(all.equal(predicted_live, selected_live, tolerance = 1e-12))) {
          stop(
            "Death likelihood live prediction does not use the existing selected_index for segment ",
            seg$segment_id, "."
          )
        }
        predicted_counts <- c(
          predicted_live, predicted_dead_hypoxia, predicted_dead_buffer, predicted_dead
        )
        if (any(!is.finite(predicted_counts)) || predicted_live <= 0 ||
            any(predicted_counts[-1L] < 0)) {
          stop("Death likelihood encountered invalid selected live/dead stock for segment ", seg$segment_id, ".")
        }

        data.frame(
          model_passage_id = passage_ids,
          matched_old_segment_id = rep(as.character(seg$segment_id), length(passage_ids)),
          cohort = rep(as.character(seg$cohort), length(passage_ids)),
          selected_index = rep(selected_index, length(passage_ids)),
          selected_day = rep(selected_day, length(passage_ids)),
          predicted_live_count = rep(predicted_live, length(passage_ids)),
          predicted_dead_hypoxia_count = rep(predicted_dead_hypoxia, length(passage_ids)),
          predicted_dead_buffer_count = rep(predicted_dead_buffer, length(passage_ids)),
          predicted_dead_count = rep(predicted_dead, length(passage_ids)),
          stringsAsFactors = FALSE
        )
      })
      dplyr::bind_rows(rows)
    }

    passage_map <- dplyr::bind_rows(collect_run_map(run_2N), collect_run_map(run_4N))
    passage_match_counts <- table(passage_map$model_passage_id)
    matched_counts <- unname(passage_match_counts[death_data$model_passage_id])
    if (any(is.na(matched_counts)) || any(matched_counts != 1L)) {
      stop("Every Death observation must match exactly one old adapter segment by model_passage_id.")
    }
    hit <- match(death_data$model_passage_id, passage_map$model_passage_id)
    matched <- passage_map[hit, , drop = FALSE]
    if (any(death_data$cohort != matched$cohort)) {
      stop("Death observation cohort disagrees with its matched old adapter segment.")
    }

    segment_key <- paste(matched$cohort, matched$matched_old_segment_id, sep = "::")
    segment_counts <- table(segment_key)
    if (length(segment_counts) != 45L || any(segment_counts != 2L)) {
      stop("Death likelihood mapping must yield 45 old segments with exactly two observations each.")
    }
    lineage_by_segment <- split(as.character(death_data$lineage_id), segment_key)
    lineage_ok <- vapply(
      lineage_by_segment,
      function(x) identical(sort(x), c("O1", "O2")),
      logical(1)
    )
    if (any(!lineage_ok)) {
      stop("Each Death likelihood segment must contain one O1 and one O2 observation.")
    }

    observed_dead_fraction <- as.numeric(death_data$observed_dead_fraction)
    predicted_dead_fraction <- matched$predicted_dead_count /
      (matched$predicted_live_count + matched$predicted_dead_count)
    clamp_probability <- function(x) pmin(pmax(as.numeric(x), eps_use), 1 - eps_use)
    observed_dead_logit <- stats::qlogis(clamp_probability(observed_dead_fraction))
    predicted_dead_logit <- stats::qlogis(clamp_probability(predicted_dead_fraction))
    row_loglik <- stats::dnorm(
      observed_dead_logit,
      mean = predicted_dead_logit,
      sd = sigma_use,
      log = TRUE
    )
    if (any(!is.finite(observed_dead_logit)) || any(!is.finite(predicted_dead_logit)) ||
        any(!is.finite(row_loglik))) {
      stop("Death likelihood produced a non-finite logit or row log-likelihood.")
    }

    data.frame(
      observation_id = as.character(death_data$observation_id),
      model_passage_id = as.character(death_data$model_passage_id),
      segment_id = matched$matched_old_segment_id,
      matched_old_segment_id = matched$matched_old_segment_id,
      cohort = as.character(death_data$cohort),
      lineage_id = as.character(death_data$lineage_id),
      selected_index = matched$selected_index,
      selected_day = matched$selected_day,
      observed_dead_count = as.numeric(death_data$dead_count),
      observed_live_count = as.numeric(death_data$observed_live_count),
      eligible_denominator = as.numeric(death_data$eligible_denominator),
      predicted_dead_hypoxia_count = matched$predicted_dead_hypoxia_count,
      predicted_dead_buffer_count = matched$predicted_dead_buffer_count,
      predicted_dead_count = matched$predicted_dead_count,
      predicted_live_count = matched$predicted_live_count,
      observed_dead_fraction = observed_dead_fraction,
      predicted_dead_fraction = predicted_dead_fraction,
      observed_dead_logit = observed_dead_logit,
      predicted_dead_logit = predicted_dead_logit,
      logit_residual = observed_dead_logit - predicted_dead_logit,
      sigma_death_logit = rep(sigma_use, length(row_loglik)),
      death_fraction_eps = rep(eps_use, length(row_loglik)),
      loglik = row_loglik,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  function(run_params,
           fit_objects,
           cfg,
           fallback_max_passage_days = 14,
           growth_weight = 1,
           ploidy_weight = 1,
           flow_weight = 1,
           ploidy_prob_floor = 1e-12,
           flow_density_floor = 1e-12,
           flow_kernel_sd_ploidy = NULL) {
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

  sigma_growth <- as.numeric(run_params$sigma_growth)
  sigma_kary <- as.numeric(run_params$sigma_kary)
  sigma_flow_ploidy <- if (is.null(flow_kernel_sd_ploidy)) {
    ivt_default_flow_kernel_sd_ploidy(run = run_2N, fit_data = fit_objects$fit_data)
  } else {
    as.numeric(flow_kernel_sd_ploidy)
  }
  growth_df <- ivt_growth_loglik_df(summary_df = summary_df, sigma_growth = sigma_growth)
  n_growth_observed <- sum(is.finite(summary_df$observed_growth))
  n_growth_missing_pred <- sum(is.finite(summary_df$observed_growth) & !is.finite(summary_df$predicted_growth_rate))
  n_growth_negative_pred <- sum(is.finite(summary_df$predicted_growth_rate) & summary_df$predicted_growth_rate < 0)
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
  death_enabled <- isTRUE(fit_objects$death_enabled)
  sigma_death_logit <- 0.75
  death_fraction_eps <- 1e-4
  death_weight <- if (death_enabled) 1.0 else 0.0
  death_df <- if (death_enabled) {
    if (is.null(fit_objects$death_data) || !is.data.frame(fit_objects$death_data) ||
        nrow(fit_objects$death_data) == 0L) {
      stop("Death likelihood is enabled but fit_objects$death_data has no validated observations.")
    }
    ivt_death_loglik_df(
      run_2N = run_2N,
      run_4N = run_4N,
      death_data = fit_objects$death_data,
      sigma_death_logit = sigma_death_logit,
      death_fraction_eps = death_fraction_eps
    )
  } else {
    ivt_empty_death_loglik_df()
  }

  growth_loglik_sum <- if (nrow(growth_df) > 0L) sum(growth_df$loglik) else 0
  ploidy_loglik_sum <- if (nrow(ploidy_df) > 0L) sum(ploidy_df$mean_loglik) else 0
  flow_loglik_sum <- if (nrow(flow_df) > 0L) sum(flow_df$mean_loglik) else 0
  death_loglik_sum <- if (nrow(death_df) > 0L) sum(death_df$loglik) else 0
  growth_loglik <- if (nrow(growth_df) > 0L) growth_loglik_sum / nrow(growth_df) else 0
  ploidy_loglik <- if (nrow(ploidy_df) > 0L) ploidy_loglik_sum / nrow(ploidy_df) else 0
  flow_loglik <- if (nrow(flow_df) > 0L) flow_loglik_sum / nrow(flow_df) else 0
  death_loglik <- if (nrow(death_df) > 0L) death_loglik_sum / nrow(death_df) else 0
  total_loglik <- as.numeric(growth_weight) * growth_loglik +
    as.numeric(ploidy_weight) * ploidy_loglik +
    as.numeric(flow_weight) * flow_loglik +
    death_weight * death_loglik
  total <- -total_loglik

  death_data_path <- if (death_enabled && length(fit_objects$death_data_path)) {
    as.character(fit_objects$death_data_path[[1]])
  } else {
    NA_character_
  }
  death_data_md5 <- if (death_enabled && length(fit_objects$death_data_md5)) {
    as.character(fit_objects$death_data_md5[[1]])
  } else {
    NA_character_
  }
  death_data_n_file_rows <- if (death_enabled && length(fit_objects$death_data_n_file_rows)) {
    as.integer(fit_objects$death_data_n_file_rows[[1]])
  } else {
    0L
  }

  list(
    objective = total,
    total_loglik = total_loglik,
    growth_loglik = growth_loglik,
    ploidy_loglik = ploidy_loglik,
    flow_loglik = flow_loglik,
    death_loglik = death_loglik,
    growth_loglik_sum = growth_loglik_sum,
    ploidy_loglik_sum = ploidy_loglik_sum,
    flow_loglik_sum = flow_loglik_sum,
    death_loglik_sum = death_loglik_sum,
    sigma_growth = sigma_growth,
    sigma_kary = sigma_kary,
    sigma_flow_ploidy = sigma_flow_ploidy,
    sigma_death_logit = sigma_death_logit,
    death_fraction_eps = death_fraction_eps,
    death_weight = death_weight,
    n_growth = nrow(growth_df),
    n_growth_observed = n_growth_observed,
    n_growth_missing_pred = n_growth_missing_pred,
    n_growth_negative_pred = n_growth_negative_pred,
    n_ploidy_passages = nrow(ploidy_df),
    n_kary_cells = sum(ploidy_df$n_cells),
    n_flow_passages = nrow(flow_df),
    n_flow_samples = nrow(flow_df),
    n_death_observations = nrow(death_df),
    death_data_path = death_data_path,
    death_data_md5 = death_data_md5,
    death_data_n_file_rows = death_data_n_file_rows,
    summary = summary_df,
    growth_df = growth_df,
    ploidy_df = ploidy_df,
    flow_df = flow_df,
    death_df = death_df,
    flow_overlay_df = flow_overlay_df,
    run_2N = run_2N,
    run_4N = run_4N
  )
  }
})

ivt_objective <- function(run_params,
                          fit_objects,
                          cfg,
                          fallback_max_passage_days = 14,
                          growth_weight = 1,
                          ploidy_weight = 1,
                          flow_weight = 1) {
  ivt_objective_components(
    run_params = run_params,
    fit_objects = fit_objects,
    cfg = cfg,
    fallback_max_passage_days = fallback_max_passage_days,
    growth_weight = growth_weight,
    ploidy_weight = ploidy_weight,
    flow_weight = flow_weight
  )$objective
}
