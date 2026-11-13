# Postprocessing only. Canonical rates and daughter kernels come from frozen v2.

pilotdiag_offspring_identity <- function(N, p, s) {
  if (any(!is.finite(c(N, p, s))) || any(N < 1 | N != round(N)) || any(p < 0 | p > 1) || any(s < 0 | s > 1))
    stop("Invalid binomial-survival diagnostic inputs")
  z <- 1 - p + p * s
  data.frame(viability_factor = z,
    effective_error_probability_conditional_on_viability = ifelse(z > 0, p * s / z, NA_real_),
    biological_viable_daughters_per_ordinary_division = 2 * z^N)
}

pilotdiag_rates <- function(ctx, params, cfg, oxygen) {
  grid <- seq.int(cfg$N_MIN, cfg$N_MAX)
  backend <- ctx$backend
  lambda <- backend$.lambda_eff_of_O2(oxygen, params, N = grid, O2_growth = cfg$O2_growth)
  mu <- backend$.mu_eff_of_O2(oxygen, params, N = grid)
  probability <- backend$.pmisseg_of_O2(oxygen, params, N = grid)
  wgd <- params$p_wgd
  delta <- get("cpp_o2simps_pr_delta_vec", envir = backend, inherits = TRUE)
  ledger <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    N <- grid[i]; m <- 0:N
    pn <- stats::dbinom(m, N, probability[i])
    s <- params$buffer_smax * exp(-params$buffer_beta * (2 * cfg$N_UNIT / N)^params$buffer_n_exp)
    survival <- s^m; survival[1] <- 1
    keep <- pn >= 1e-8 & pn > 0
    viable <- pn * survival * keep
    loss <- m > 0 & N - m >= cfg$N_MIN
    gain <- m > 0 & N + m <= cfg$N_MAX
    ordinary_live <- sum(viable[loss]) + sum(viable[gain]) + 2 * viable[1]
    biological <- 2 * sum(pn * (1 - survival))
    pruning <- 2 * sum(pn * survival * !keep)
    outside_loss <- sum(viable[m > 0 & N - m < cfg$N_MIN])
    outside_gain <- sum(viable[m > 0 & N + m > cfg$N_MAX])
    cpp <- delta(N, probability[i], eps_tail = 1e-8, buffer_smax = params$buffer_smax,
                 buffer_beta = params$buffer_beta, buffer_n_exp = params$buffer_n_exp, N_unit = cfg$N_UNIT)
    cpp_live <- sum(cpp$prob[N + cpp$ts >= cfg$N_MIN & N + cpp$ts <= cfg$N_MAX])
    if (abs(cpp_live - ordinary_live) > 1e-10 ||
        abs(2 * cpp$mass_dropped - biological - pruning) > 1e-10 ||
        abs(ordinary_live + biological + pruning + outside_loss + outside_gain - 2) > 1e-10)
      stop("Ordinary daughter ledger does not match the production kernel")
    rate <- lambda[i] * (1 - wgd)
    ledger[[i]] <- data.frame(N = N, oxygen_pct = oxygen, lambda = lambda[i], mu = mu[i], p_per_chromosome = probability[i],
      division = lambda[i], ordinary_generated_loss = rate * sum(pn[m > 0]), ordinary_generated_gain = rate * sum(pn[m > 0]),
      ordinary_retained_loss = rate * sum(viable[loss]), ordinary_retained_gain = rate * sum(viable[gain]),
      ordinary_unchanged = rate * 2 * viable[1], ordinary_biological_nonviability = rate * biological,
      numerical_tail_pruning = rate * pruning, ordinary_boundary_loss = rate * outside_loss,
      ordinary_boundary_gain = rate * outside_gain, wgd_retained = lambda[i] * wgd * as.numeric(2 * N <= cfg$N_MAX),
      wgd_boundary = lambda[i] * wgd * as.numeric(2 * N > cfg$N_MAX))
  }
  ledger <- do.call(rbind, ledger)
  tri <- get("cpp_o2simps_build_G_for_o2_triplet", envir = backend, inherits = TRUE)(
    O2 = oxygen, O2_crit = params$O2_crit, N0min = cfg$N_MIN, N0max = cfg$N_MAX,
    N1min = cfg$N_MIN, N1max = cfg$N_MAX, lam_max = params$lam_max,
    p_mis_base = params$p_mis_base, p_misseg = params$p_misseg, k_o_mis = params$k_o_mis,
    p_wgd = wgd, boundary = "drop", eps_tail = 1e-8, buffer_smax = params$buffer_smax,
    buffer_beta = params$buffer_beta, buffer_n_exp = params$buffer_n_exp, N_unit = cfg$N_UNIT,
    beta_size = params$beta_size, O2_growth = cfg$O2_growth, alpha_o2 = params$alpha_o2,
    gamma_growth = params$gamma_growth, mu_hp = params$mu_hp, gamma_mu = params$gamma_mu,
    n_O = params$n_O, ploidy_O2_death = cfg$ploidy_O2_death)
  G <- Matrix::sparseMatrix(i = tri$i, j = tri$j, x = tri$x, dims = c(length(grid), length(grid)))
  expected <- ledger$ordinary_retained_loss + ledger$ordinary_retained_gain + ledger$ordinary_unchanged + ledger$wgd_retained - ledger$division
  dead <- ledger$ordinary_biological_nonviability + ledger$numerical_tail_pruning + ledger$ordinary_boundary_loss + ledger$ordinary_boundary_gain + ledger$wgd_boundary
  if (max(abs(as.numeric(Matrix::colSums(G)) - expected)) > 1e-9 || max(abs(tri$dead_buffer_rate - dead)) > 1e-9)
    stop("Generator/dead-buffer ledger mismatch")
  list(rates = ledger, generator = G)
}

pilotdiag_dense_flux <- function(ctx, fit, cohort, result, kernel) {
  cfg <- fit$cfg; params <- fit$best_params; segment <- result$segment
  day <- result$selection$selected_day
  steps <- as.integer(round(day / cfg$DT))
  dense_segment <- segment
  dense_segment$duration_days <- day
  dense_segment$obs_days_local <- (0:steps) * cfg$DT
  core <- fit$best_components[[paste0("run_", cohort)]]$model_core
  vol <- ctx$backend$cell_volume_mm3_by_N(core$grid_pre, params, cfg)
  original_start <- result$sim$live_state_obs[1, ]
  dense <- ctx$backend$ivt_run_segment_fixed_o2(dense_segment, cfg, params, core, vol,
                                               init_state_override = original_start,
                                               init_cells_override = sum(original_start))$sim
  X <- dense$live_state_obs
  original_rows <- which(segment$obs_days_local <= day + 1e-12)
  dense_rows <- as.integer(round(segment$obs_days_local[original_rows] / cfg$DT)) + 1L
  difference <- abs(X[dense_rows, , drop = FALSE] - result$sim$live_state_obs[original_rows, , drop = FALSE])
  bound <- 1e-8 + 1e-8 * abs(result$sim$live_state_obs[original_rows, , drop = FALSE])
  if (any(difference > bound)) stop("Dense postprocessing changed the saved trajectory")
  if (any(abs(dense$O2_eff_obs - segment$oxygen_pct) > 1e-8)) stop("Unexpected oxygen history in fixed-oxygen flux audit")
  total <- rowSums(X)
  crowd <- if (!isTRUE(cfg$Crowding)) rep(1, nrow(X)) else if (cfg$crowding == "logistic") pmax(0, 1 - total / cfg$K) else exp(-total / cfg$K)
  rate_names <- setdiff(names(kernel$rates), c("N", "oxygen_pct", "lambda", "mu", "p_per_chromosome"))
  flux <- as.data.frame(as.matrix(X %*% as.matrix(kernel$rates[, rate_names, drop = FALSE])) * crowd)
  flux$direct_death <- as.numeric(X %*% kernel$rates$mu)
  divisions <- pmax(as.numeric(X %*% kernel$rates$lambda), 1e-300)
  chromosome_divisions <- pmax(as.numeric(X %*% (kernel$rates$lambda * kernel$rates$N)), 1e-300)
  errors <- as.numeric(X %*% (kernel$rates$lambda * kernel$rates$N * kernel$rates$p_per_chromosome))
  flux$cell_weighted_mean_p_per_chromosome <- as.numeric(X %*% kernel$rates$p_per_chromosome) / pmax(total, 1e-300)
  flux$chromosome_weighted_mean_p_per_chromosome <- as.numeric(X %*% (kernel$rates$N * kernel$rates$p_per_chromosome)) / pmax(as.numeric(X %*% kernel$rates$N), 1e-300)
  flux$division_weighted_mean_p_per_chromosome <- as.numeric(X %*% (kernel$rates$lambda * kernel$rates$p_per_chromosome)) / divisions
  flux$division_chromosome_weighted_mean_p <- errors / chromosome_divisions
  flux$mean_errors_per_ordinary_division <- errors / divisions
  flux$mean_chromosomes_per_division <- chromosome_divisions / divisions
  if (any(abs(flux$mean_errors_per_ordinary_division - flux$division_chromosome_weighted_mean_p * flux$mean_chromosomes_per_division) > 1e-10))
    stop("Division-weighted error-summary identity failed")
  correction <- numeric(steps)
  maximum_euler_difference <- 0
  for (j in seq_len(steps)) {
    raw <- X[j, ] + cfg$DT * (crowd[j] * as.numeric(kernel$generator %*% X[j, ]) - kernel$rates$mu * X[j, ])
    expected <- pmax(raw, 0)
    delta <- abs(expected - X[j + 1L, ])
    if (any(delta > 1e-8 + 1e-8 * abs(X[j + 1L, ]))) stop("Dense Euler update does not match the event ledger")
    maximum_euler_difference <- max(maximum_euler_difference, delta)
    correction[j] <- sum(expected - raw)
  }
  integrated_names <- c(rate_names, "direct_death")
  integrated <- if (steps) colSums(flux[seq_len(steps), integrated_names, drop = FALSE]) * cfg$DT else setNames(rep(0, length(integrated_names)), integrated_names)
  balance <- integrated[["ordinary_retained_loss"]] + integrated[["ordinary_retained_gain"]] + integrated[["ordinary_unchanged"]] +
    integrated[["wgd_retained"]] - integrated[["division"]] - integrated[["direct_death"]] + sum(correction)
  residual <- total[length(total)] - total[1] - balance
  if (abs(residual) > 1e-7 + 1e-8 * max(total)) stop("Integrated live-cell event balance failed")
  flux$time_days <- dense_segment$obs_days_local
  flux$live_cells <- total
  flux$cohort <- cohort
  flux$segment_id <- segment$segment_id
  flux$oxygen_pct <- segment$oxygen_pct
  list(flux = flux, summary = data.frame(cohort = cohort, segment_id = segment$segment_id,
      selected_observation_day = day, oxygen_pct = segment$oxygen_pct, initial_live_cells = total[1],
      final_live_cells = total[length(total)], as.list(integrated), numerical_nonnegative_correction_cells = sum(correction),
      live_balance_residual_cells = residual, maximum_dense_Euler_state_difference = maximum_euler_difference),
       daily_replay_maximum_difference = max(difference))
}
