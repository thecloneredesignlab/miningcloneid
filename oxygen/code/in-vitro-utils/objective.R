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
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(cfg$misseg_loss_survival, "nullisomy"),
    "nullisomy"
  )

  expected_fit_symbols <- c(
    "lam_min", "lam_max", "k_o", "p_misseg", "k_o_mis",
    if (identical(loss_mode, "nullisomy")) "gamma_loss" else c("buffer_smax", "buffer_beta", "buffer_n_exp"),
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

  loss_param_name <- if (identical(loss_mode, "nullisomy")) {
    "log10_gamma_loss"
  } else {
    c("buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp")
  }
  loss_lower <- if (identical(loss_mode, "nullisomy")) {
    log10(positive_bound("gamma_loss", "lower"))
  } else {
    c(
      as.numeric(row_for("buffer_smax")$lower_bound[[1]]),
      log10(positive_bound("buffer_beta", "lower")),
      log10(positive_bound("buffer_n_exp", "lower"))
    )
  }
  loss_upper <- if (identical(loss_mode, "nullisomy")) {
    log10(positive_bound("gamma_loss", "upper"))
  } else {
    c(
      as.numeric(row_for("buffer_smax")$upper_bound[[1]]),
      log10(positive_bound("buffer_beta", "upper")),
      log10(positive_bound("buffer_n_exp", "upper"))
    )
  }
  loss_init <- if (identical(loss_mode, "nullisomy")) {
    log10(positive_bound("gamma_loss", "init"))
  } else {
    c(
      as.numeric(row_for("buffer_smax")$init_value[[1]]),
      log10(positive_bound("buffer_beta", "init")),
      log10(positive_bound("buffer_n_exp", "init"))
    )
  }

  data.frame(
    param_name = c(
      "log10_lam_min", "log10_lam_max", "log10_k_o", "logit_p_misseg",
      "log10_k_o_mis", loss_param_name, "logit_p_wgd",
      "log10_alpha_o2", "gamma_growth", "log10_mu_hp",
      "gamma_mu", "log10_O2_crit", "n_O", "log10_p_mis_base",
      "log10_sigma_growth", "log10_sigma_kary", "init_mean_2N", "log10_init_sd_2N",
      "init_mean_4N", "log10_init_sd_4N"
    ),
    lower = c(
      log10(positive_bound("lam_min", "lower")),
      log10(positive_bound("lam_max", "lower")),
      log10(positive_bound("k_o", "lower")),
      qlogis(as.numeric(row_for("p_misseg")$lower_bound[[1]])),
      log10(positive_bound("k_o_mis", "lower")),
      loss_lower,
      qlogis(as.numeric(row_for("p_wgd")$lower_bound[[1]])),
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
      log10(positive_bound("lam_min", "upper")),
      log10(positive_bound("lam_max", "upper")),
      log10(positive_bound("k_o", "upper")),
      qlogis(as.numeric(row_for("p_misseg")$upper_bound[[1]])),
      log10(positive_bound("k_o_mis", "upper")),
      loss_upper,
      qlogis(as.numeric(row_for("p_wgd")$upper_bound[[1]])),
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
      log10(positive_bound("lam_min", "init")),
      log10(positive_bound("lam_max", "init")),
      log10(positive_bound("k_o", "init")),
      qlogis(as.numeric(row_for("p_misseg")$init_value[[1]])),
      log10(positive_bound("k_o_mis", "init")),
      loss_init,
      qlogis(as.numeric(row_for("p_wgd")$init_value[[1]])),
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
  p_misseg <- as.numeric(run_params$p_misseg)
  if (!is.finite(p_misseg) || p_misseg <= 0 || p_misseg >= 1) {
    stop("run_params$p_misseg must lie in (0,1) for logit-scale optimization.")
  }
  par_t[["logit_p_misseg"]] <- qlogis(p_misseg)
  par_t[["log10_k_o_mis"]] <- log10(require_positive(run_params$k_o_mis, "k_o_mis"))
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(cfg$misseg_loss_survival, run_params$misseg_loss_survival, "nullisomy"),
    "nullisomy"
  )
  if (identical(loss_mode, "nullisomy") && "log10_gamma_loss" %in% names(par_t)) {
    par_t[["log10_gamma_loss"]] <- log10(require_positive(run_params$gamma_loss, "gamma_loss"))
  }
  if (identical(loss_mode, "buffering")) {
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
  }
  p_wgd <- as.numeric(run_params$p_wgd)
  if (!is.finite(p_wgd) || p_wgd <= 0 || p_wgd >= 1) {
    stop("run_params$p_wgd must lie in (0,1) for logit-scale optimization.")
  }
  par_t[["logit_p_wgd"]] <- qlogis(p_wgd)
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
  run_params$lam_min <- 10^par_t[["log10_lam_min"]]
  run_params$lam_max <- 10^par_t[["log10_lam_max"]]
  run_params$k_o <- 10^par_t[["log10_k_o"]]
  run_params$p_misseg <- plogis(par_t[["logit_p_misseg"]])
  run_params$k_o_mis <- 10^par_t[["log10_k_o_mis"]]
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(cfg$misseg_loss_survival, run_params$misseg_loss_survival, "nullisomy"),
    "nullisomy"
  )
  run_params$misseg_loss_survival <- loss_mode
  if (identical(loss_mode, "nullisomy") && "log10_gamma_loss" %in% names(par_t)) {
    run_params$gamma_loss <- 10^par_t[["log10_gamma_loss"]]
  }
  if (identical(loss_mode, "buffering")) {
    if ("buffer_smax" %in% names(par_t)) run_params$buffer_smax <- par_t[["buffer_smax"]]
    if ("log10_buffer_beta" %in% names(par_t)) run_params$buffer_beta <- 10^par_t[["log10_buffer_beta"]]
    if ("log10_buffer_n_exp" %in% names(par_t)) run_params$buffer_n_exp <- 10^par_t[["log10_buffer_n_exp"]]
  }
  run_params$p_wgd <- plogis(par_t[["logit_p_wgd"]])
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

ivt_objective_components <- function(run_params,
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

  growth_loglik_sum <- if (nrow(growth_df) > 0L) sum(growth_df$loglik) else 0
  ploidy_loglik_sum <- if (nrow(ploidy_df) > 0L) sum(ploidy_df$mean_loglik) else 0
  flow_loglik_sum <- if (nrow(flow_df) > 0L) sum(flow_df$mean_loglik) else 0
  growth_loglik <- if (nrow(growth_df) > 0L) growth_loglik_sum / nrow(growth_df) else 0
  ploidy_loglik <- if (nrow(ploidy_df) > 0L) ploidy_loglik_sum / nrow(ploidy_df) else 0
  flow_loglik <- if (nrow(flow_df) > 0L) flow_loglik_sum / nrow(flow_df) else 0
  total_loglik <- as.numeric(growth_weight) * growth_loglik +
    as.numeric(ploidy_weight) * ploidy_loglik +
    as.numeric(flow_weight) * flow_loglik
  total <- -total_loglik

  list(
    objective = total,
    total_loglik = total_loglik,
    growth_loglik = growth_loglik,
    ploidy_loglik = ploidy_loglik,
    flow_loglik = flow_loglik,
    growth_loglik_sum = growth_loglik_sum,
    ploidy_loglik_sum = ploidy_loglik_sum,
    flow_loglik_sum = flow_loglik_sum,
    sigma_growth = sigma_growth,
    sigma_kary = sigma_kary,
    sigma_flow_ploidy = sigma_flow_ploidy,
    n_growth = nrow(growth_df),
    n_growth_observed = n_growth_observed,
    n_growth_missing_pred = n_growth_missing_pred,
    n_growth_negative_pred = n_growth_negative_pred,
    n_ploidy_passages = nrow(ploidy_df),
    n_kary_cells = sum(ploidy_df$n_cells),
    n_flow_passages = nrow(flow_df),
    n_flow_samples = nrow(flow_df),
    summary = summary_df,
    growth_df = growth_df,
    ploidy_df = ploidy_df,
    flow_df = flow_df,
    flow_overlay_df = flow_overlay_df,
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
