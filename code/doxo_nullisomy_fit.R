#!/usr/bin/env Rscript

# Exploratory standalone nullisomy-only fitter.
# This script reuses the active O2_supply_demand_MAP backend, but it is not yet
# integrated into the main optimizer/config path used by the primary invivo
# fitting workflow.

suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

repo_root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
workflow_root <- file.path(repo_root, "oxygen", "code", "O2_supply_demand_MAP")
model_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")
shared_path <- file.path(workflow_root, "util", "o2_supply_demand_map_shared.R")
if (!file.exists(model_path)) stop("Missing active model backend: ", model_path)
if (!file.exists(shared_path)) stop("Missing active shared utils: ", shared_path)
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
source(shared_path, local = TRUE)
source(model_path, local = TRUE)

parse_args <- function(argv) {
  out <- list()
  if (length(argv) == 0L) return(out)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a)
    a <- substring(a, 3L)
    if (grepl("=", a, fixed = TRUE)) {
      kv <- strsplit(a, "=", fixed = TRUE)[[1]]
      out[[kv[[1]]]] <- paste(kv[-1], collapse = "=")
      i <- i + 1L
    } else {
      key <- a
      if (i == length(argv)) stop("Missing value for --", key)
      out[[key]] <- argv[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v)) default else v
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) default else v
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

default_data_path <- function() {
  file.path(
    repo_root,
    "data",
    "DrugResponseData_PloidyJumps",
    "DrugResponseData",
    "SNU668_doxorubicin_response_and_karyotype_ploidy.json"
  )
}

timestamp_string <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

make_output_dir <- function(out_dir = NULL) {
  path <- if (!is.null(out_dir) && nzchar(out_dir)) {
    out_dir
  } else {
    file.path(getwd(), "doxo_outputs", paste0("doxo_nullisomy_", timestamp_string()))
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

load_doxo_data <- function(json_path) {
  if (!file.exists(json_path)) stop("JSON not found: ", json_path)
  dat <- jsonlite::fromJSON(json_path, simplifyDataFrame = TRUE)
  need <- c(
    "drug_response_summary_by_dose",
    "drug_response_replicate_level",
    "cell_ploidy",
    "sample_ploidy_summary"
  )
  miss <- setdiff(need, names(dat))
  if (length(miss) > 0L) stop("JSON missing fields: ", paste(miss, collapse = ", "))
  dat
}

default_cfg <- function() {
  list(
    mode = "doxorubicin_nullisomy",
    N_UNIT = 22L,
    N_MIN = 22L,
    N_MAX = 154L,
    DT = 0.25,
    start_with = "chr_number",
    chr_lengths_bp = default_chr_lengths_bp_1to22(),
    Crowding = FALSE,
    crowding = "logistic",
    K = 1e12,
    min_pop = 1e-12,
    fit_treatment = FALSE,
    dose_ref = 1.0,
    tx_mult_min = 1.0,
    o2_burden_feedback = FALSE,
    O2_growth = FALSE,
    o2_cache_bin_pct = 0.01,
    o2_cache_hysteresis_pct = 0.005,
    o2_cache_profile = FALSE,
    p_mis_base = 1e-5,
    gamma_loss = 1.0,
    p_wgd = 1e-8,
    rho_2N = 42332.020717,
    o2_fixed = 5.0,
    o2_S0 = 5.0,
    o2_min = 5.0,
    kappa_O = 1.0,
    tau_O2 = 1.0,
    o2_Nref = 1e6,
    k_o = 1.0,
    alpha_o2 = 0.5,
    gamma_growth = 2.0,
    O2_crit = 1.0,
    n_O = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    k_clear = 1e-3,
    burden_log_eps = 1e-12,
    sigma_log_response_floor = 0.05,
    assay_days_fixed = 4.0,
    seeds_default = c(1L, 2L, 3L, 4L, 5L)
  )
}

hill_pmis <- function(dose_uM, p_mis_base, p_mis_doxo_max, EC50, hill_n) {
  dose <- pmax(as.numeric(dose_uM), 0)
  p_base <- max(as.numeric(p_mis_base), 0)
  p_amp <- max(as.numeric(p_mis_doxo_max), 0)
  ec50 <- max(as.numeric(EC50), 1e-12)
  n <- max(as.numeric(hill_n), 1e-12)
  frac <- dose^n / (ec50^n + dose^n)
  clip01(p_base + p_amp * frac)
}

expected_mis_copies <- function(N, p) {
  as.numeric(N) * as.numeric(p)
}

sample_initial_state <- function(cell_ploidy, sample_name, cfg) {
  grid <- seq.int(cfg$N_MIN, cfg$N_MAX)
  sub <- cell_ploidy[cell_ploidy$sample == sample_name, , drop = FALSE]
  vals <- suppressWarnings(as.numeric(sub$chromosome_total))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    vals <- rep(46, 50L)
  }
  vals <- pmin(pmax(round(vals), cfg$N_MIN), cfg$N_MAX)
  tab <- table(factor(vals, levels = grid))
  prob <- as.numeric(tab)
  prob <- prob / sum(prob)
  names(prob) <- as.character(grid)
  prob
}

make_run_params <- function(theta, cfg) {
  lambda_eff <- max(as.numeric(theta$lambda_eff), 1e-6)
  list(
    lam_min = lambda_eff,
    lam_max = lambda_eff,
    k_o = cfg$k_o,
    o2_min = cfg$o2_min,
    p_mis_base = cfg$p_mis_base,
    p_misseg = NULL,
    k_o_mis = 1.0,
    gamma_loss = cfg$gamma_loss,
    p_wgd = cfg$p_wgd,
    o2_S0 = cfg$o2_S0,
    kappa_O = cfg$kappa_O,
    eta_o2 = 1.0,
    rho_2N = cfg$rho_2N,
    alpha_o2 = cfg$alpha_o2,
    gamma_growth = cfg$gamma_growth,
    mu_hp = 0.0,
    gamma_mu = cfg$gamma_mu,
    O2_crit = cfg$O2_crit,
    n_O = cfg$n_O,
    k_clear = cfg$k_clear,
    tau_O2 = cfg$tau_O2,
    alpha = 0.0,
    gamma = 1.0,
    ploidy_O2_death = "uniform",
    start_with = "chr_number",
    boundary = "drop"
  )
}

simulate_one_doxo <- function(init_state_prob, dose_uM, theta, cfg) {
  grid <- seq.int(cfg$N_MIN, cfg$N_MAX)
  init_total <- cfg$o2_Nref
  init_state <- as.numeric(init_state_prob) * init_total
  run_params <- make_run_params(theta, cfg)
  p_const <- hill_pmis(
    dose_uM = dose_uM,
    p_mis_base = cfg$p_mis_base,
    p_mis_doxo_max = theta$p_mis_doxo_max,
    EC50 = theta$EC50,
    hill_n = theta$hill_n
  )
  obs_steps <- as.integer(round(cfg$assay_days_fixed / cfg$DT))
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid, run_params = run_params, cfg = cfg))
  sim <- cpp_o2simps_simulate_one(
    init_state = as.numeric(init_state),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = as.integer(obs_steps),
    sim_end_step = as.integer(obs_steps),
    DT = as.numeric(cfg$DT),
    dose = as.numeric(dose_uM),
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = 0.0,
    fit_treatment = FALSE,
    alpha = 0.0,
    gamma = 1.0,
    tx_mult_min = 1.0,
    crowding_enabled = FALSE,
    crowding = "logistic",
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_crit = as.numeric(cfg$O2_crit),
    o2_feedback = FALSE,
    o2_S0 = as.numeric(cfg$o2_S0),
    kappa_O = as.numeric(cfg$kappa_O),
    tau_O2 = as.numeric(cfg$tau_O2),
    o2_Nref = as.numeric(cfg$o2_Nref),
    o2_min = as.numeric(cfg$o2_min),
    eta_o2 = 1.0,
    o2_cache_bin_pct = as.numeric(cfg$o2_cache_bin_pct),
    o2_cache_hysteresis_pct = as.numeric(cfg$o2_cache_hysteresis_pct),
    o2_cache_profile = FALSE,
    O2_growth = FALSE,
    lam_min = as.numeric(run_params$lam_min),
    lam_max = as.numeric(run_params$lam_max),
    k_o = as.numeric(run_params$k_o),
    has_p_misseg = FALSE,
    p_mis_base = as.numeric(cfg$p_mis_base),
    p_misseg = 0.0,
    k_o_mis = 1.0,
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = as.numeric(p_const),
    p_wgd = as.numeric(cfg$p_wgd),
    boundary = "drop",
    eps_tail = 1e-8,
    gamma_loss = as.numeric(cfg$gamma_loss),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = 0.0,
    alpha_o2 = as.numeric(cfg$alpha_o2),
    gamma_growth = as.numeric(cfg$gamma_growth),
    mu_hp = 0.0,
    gamma_mu = as.numeric(cfg$gamma_mu),
    n_O = as.numeric(cfg$n_O),
    ploidy_O2_death = "uniform",
    start_with = "chr_number",
    k_clear = as.numeric(cfg$k_clear),
    vol_by_N = vol_by_N,
    burden_floor = as.numeric(cfg$burden_log_eps),
    return_full_trajectory = FALSE
  )
  list(
    p_const = p_const,
    live_burden = as.numeric(sim$Vmm3_live_obs[[1]]),
    total_burden = as.numeric(sim$Vmm3_total_obs[[1]]),
    dead_hypoxia = as.numeric(sim$Vmm3_dead_hypoxia_obs[[1]]),
    dead_buffer = as.numeric(sim$Vmm3_dead_buffer_obs[[1]]),
    live_fraction = as.numeric(sim$frac_N_live),
    sim = sim
  )
}

predict_dataset <- function(theta, data_obj, cfg) {
  samples <- unique(as.character(data_obj$drug_response_summary_by_dose$sample))
  init_states <- list(
    "SNU-668_P1_A19kT_harvest" = sample_initial_state(data_obj$cell_ploidy, "SNU-668_P1_A19kT_harvest", cfg),
    "SNU-668_r2_GFP_VT_A7_harvesT3" = sample_initial_state(data_obj$cell_ploidy, "SNU-668_r2_GFP_VT_A7_harvesT3", cfg)
  )
  rows <- vector("list", length = 0L)
  for (sample_name in samples) {
    sub <- data_obj$drug_response_summary_by_dose[
      data_obj$drug_response_summary_by_dose$sample == sample_name,
      ,
      drop = FALSE
    ]
    sub <- sub[order(sub$dose_uM), , drop = FALSE]
    sims <- lapply(sub$dose_uM, function(dose) simulate_one_doxo(init_states[[sample_name]], dose, theta, cfg))
    base_live <- sims[[which.min(abs(sub$dose_uM - 0))]]$live_burden
    base_live <- max(base_live, cfg$burden_log_eps)
    for (i in seq_len(nrow(sub))) {
      rows[[length(rows) + 1L]] <- data.frame(
        sample = sample_name,
        dose_uM = sub$dose_uM[[i]],
        observed = sub$mean_normalized_to_0uM[[i]],
        observed_sd = sub$sd_normalized_to_0uM[[i]],
        predicted = sims[[i]]$live_burden / base_live,
        p_const = sims[[i]]$p_const,
        dead_hypoxia = sims[[i]]$dead_hypoxia,
        dead_buffer = sims[[i]]$dead_buffer,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

response_sigma <- function(obs, obs_sd, sigma_floor = 0.05) {
  obs_use <- pmax(as.numeric(obs), 1e-8)
  sd_use <- as.numeric(obs_sd)
  sigma <- ifelse(is.finite(sd_use) & sd_use > 0, sd_use / obs_use, sigma_floor)
  pmax(sigma, sigma_floor)
}

objective_from_table <- function(pred_tab, cfg) {
  obs <- pmax(as.numeric(pred_tab$observed), cfg$burden_log_eps)
  pred <- pmax(as.numeric(pred_tab$predicted), cfg$burden_log_eps)
  sigma <- response_sigma(obs, pred_tab$observed_sd, cfg$sigma_log_response_floor)
  z <- (log(pred) - log(obs)) / sigma
  0.5 * sum(z^2 + log(2 * pi * sigma^2), na.rm = TRUE)
}

theta_from_par <- function(par, cfg) {
  p_headroom <- max(1 - cfg$p_mis_base - 1e-8, 1e-8)
  list(
    p_mis_doxo_max = p_headroom * plogis(par[[1]]),
    EC50 = 10^par[[2]],
    hill_n = 10^par[[3]],
    lambda_eff = 10^par[[4]]
  )
}

par_from_theta <- function(theta, cfg) {
  p_headroom <- max(1 - cfg$p_mis_base - 1e-8, 1e-8)
  frac <- clip01(theta$p_mis_doxo_max / p_headroom)
  c(
    qlogis(pmin(pmax(frac, 1e-8), 1 - 1e-8)),
    log10(theta$EC50),
    log10(theta$hill_n),
    log10(theta$lambda_eff)
  )
}

default_theta <- function(cfg) {
  list(
    p_mis_doxo_max = 0.05,
    EC50 = 0.2,
    hill_n = 1.5,
    lambda_eff = 0.2
  )
}

fit_one_seed <- function(seed, data_obj, cfg) {
  set.seed(seed)
  base <- default_theta(cfg)
  jittered <- list(
    p_mis_doxo_max = clip01(base$p_mis_doxo_max * 10^runif(1, -0.5, 0.5)),
    EC50 = base$EC50 * 10^runif(1, -0.75, 0.75),
    hill_n = base$hill_n * 10^runif(1, -0.5, 0.5),
    lambda_eff = base$lambda_eff * 10^runif(1, -0.4, 0.4)
  )
  init_par <- par_from_theta(jittered, cfg)
  lower <- c(qlogis(1e-8), log10(1e-4), log10(0.2), log10(1e-3))
  upper <- c(qlogis(1 - 1e-8), log10(100), log10(8), log10(2))
  fn <- function(par) {
    theta <- theta_from_par(par, cfg)
    pred <- predict_dataset(theta, data_obj, cfg)
    objective_from_table(pred, cfg)
  }
  opt <- optim(
    par = init_par,
    fn = fn,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 400)
  )
  theta_hat <- theta_from_par(opt$par, cfg)
  pred_hat <- predict_dataset(theta_hat, data_obj, cfg)
  list(
    seed = seed,
    objective = opt$value,
    convergence = opt$convergence,
    message = opt$message,
    theta = theta_hat,
    predictions = pred_hat
  )
}

diagnostic_curves <- function(cfg, out_dir) {
  doses <- c(0, exp(seq(log(1e-4), log(10), length.out = 200L)))
  doses_plot <- pmax(doses, 1e-4)
  p_base <- cfg$p_mis_base
  pmaxs <- c(0.01, 0.05, 0.15)
  ec50s <- c(0.05, 0.2, 1.0)
  hill_ns <- c(1, 2)
  curve_rows <- vector("list", length = 0L)
  for (pmax in pmaxs) {
    for (ec50 in ec50s) {
      for (hn in hill_ns) {
        curve_rows[[length(curve_rows) + 1L]] <- data.frame(
          dose_uM = doses,
          dose_uM_plot = doses_plot,
          p_mis = hill_pmis(doses, p_base, pmax, ec50, hn),
          spec = sprintf("pmax=%.3f EC50=%.3f n=%.1f", pmax, ec50, hn),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  curve_tab <- bind_rows(curve_rows)
  g1 <- ggplot(curve_tab, aes(dose_uM_plot, p_mis, color = spec)) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() +
    labs(title = "p_mis vs doxorubicin dose", x = "Dose (uM)", y = "Per-copy p_mis")
  ggsave(file.path(out_dir, "diag_pmis_vs_dose.png"), g1, width = 8, height = 5, dpi = 150)

  states <- c(46, 53, 90, 95)
  expected_tab <- bind_rows(lapply(states, function(N) {
    data.frame(
      dose_uM = doses,
      dose_uM_plot = doses_plot,
      expected_mis = expected_mis_copies(N, hill_pmis(doses, p_base, 0.05, 0.2, 2)),
      N = as.factor(N),
      stringsAsFactors = FALSE
    )
  }))
  g2 <- ggplot(expected_tab, aes(dose_uM_plot, expected_mis, color = N)) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() +
    labs(title = "Expected missegregated copies per division", x = "Dose (uM)", y = "N * p_mis(C)")
  ggsave(file.path(out_dir, "diag_expected_mis_copies.png"), g2, width = 8, height = 5, dpi = 150)

  p_grid <- exp(seq(log(1e-5), log(0.5), length.out = 120L))
  nulli_tab <- bind_rows(lapply(states, function(N) {
    dead_prob <- vapply(p_grid, function(p) {
      attr(.pr_delta_vec(N, p, gamma_loss = cfg$gamma_loss, N_unit = cfg$N_UNIT), "mass_dropped")
    }, numeric(1))
    data.frame(p_mis = p_grid, dead_prob = dead_prob, N = as.factor(N), stringsAsFactors = FALSE)
  }))
  g3 <- ggplot(nulli_tab, aes(p_mis, dead_prob, color = N)) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() +
    labs(title = "Nullisomy/nonviable daughter mass vs p_mis", x = "Per-copy p_mis", y = "Dropped daughter mass")
  ggsave(file.path(out_dir, "diag_nullisomy_vs_pmis.png"), g3, width = 8, height = 5, dpi = 150)

  invisible(list(curve_tab = curve_tab, expected_tab = expected_tab, nulli_tab = nulli_tab))
}

run_tests <- function(data_obj, cfg, out_dir) {
  doses <- c(0, 0.01, 0.1, 1.0, 10.0)
  p <- hill_pmis(doses, cfg$p_mis_base, 0.05, 0.2, 2.0)
  stopifnot(all(diff(p) >= -1e-12))
  stopifnot(abs(p[[1]] - cfg$p_mis_base) < 1e-12)
  stopifnot(all(p <= cfg$p_mis_base + 0.05 + 1e-12))

  p_a <- hill_pmis(c(0.01, 1.0), cfg$p_mis_base, 0.05, 0.2, 2.0)
  stopifnot(expected_mis_copies(46, p_a[[2]]) > expected_mis_copies(46, p_a[[1]]))
  stopifnot(expected_mis_copies(95, 0.02) > expected_mis_copies(46, 0.02))

  surv_low <- .loss_survival_nullisomy(46, m_loss = 1, gamma_loss = cfg$gamma_loss, N_unit = cfg$N_UNIT)
  surv_high <- .loss_survival_nullisomy(95, m_loss = 1, gamma_loss = cfg$gamma_loss, N_unit = cfg$N_UNIT)
  stopifnot(surv_high >= surv_low)

  theta <- default_theta(cfg)
  init_r2 <- sample_initial_state(data_obj$cell_ploidy, "SNU-668_r2_GFP_VT_A7_harvesT3", cfg)
  sim <- simulate_one_doxo(init_r2, dose_uM = 1.0, theta = theta, cfg = cfg)
  stopifnot(abs(sim$dead_hypoxia) < 1e-12)

  writeLines("All doxorubicin-nullisomy checks passed.", con = file.path(out_dir, "tests_ok.txt"))
  invisible(TRUE)
}

write_fit_outputs <- function(best_fit, data_obj, cfg, out_dir) {
  pred <- best_fit$predictions
  residual <- transform(pred, residual_log = log(pmax(predicted, cfg$burden_log_eps)) - log(pmax(observed, cfg$burden_log_eps)))
  pred$dose_uM_plot <- pmax(pred$dose_uM, 1e-4)
  residual$dose_uM_plot <- pmax(residual$dose_uM, 1e-4)
  params <- data.frame(
    parameter = c("p_mis_base", "p_mis_doxo_max", "EC50_uM", "hill_n", "lambda_eff", "assay_days_fixed", "objective"),
    value = c(cfg$p_mis_base, best_fit$theta$p_mis_doxo_max, best_fit$theta$EC50, best_fit$theta$hill_n, best_fit$theta$lambda_eff, cfg$assay_days_fixed, best_fit$objective),
    stringsAsFactors = FALSE
  )
  utils::write.table(params, file.path(out_dir, "best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(residual, file.path(out_dir, "predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  g_fit <- ggplot(pred, aes(dose_uM_plot, observed, color = sample)) +
    geom_point(size = 2) +
    geom_line(aes(y = predicted), linewidth = 0.9) +
    geom_errorbar(aes(ymin = pmax(observed - observed_sd, 1e-6), ymax = observed + observed_sd), width = 0.02, alpha = 0.6) +
    scale_x_log10() +
    labs(title = "Observed vs predicted normalized doxorubicin response", x = "Dose (uM)", y = "Normalized viable burden")
  ggsave(file.path(out_dir, "fit_observed_vs_predicted.png"), g_fit, width = 8, height = 5, dpi = 150)

  g_res <- ggplot(residual, aes(dose_uM_plot, residual_log, color = sample)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 2) +
    scale_x_log10() +
    labs(title = "Residuals by sample and dose", x = "Dose (uM)", y = "log(pred) - log(obs)")
  ggsave(file.path(out_dir, "fit_residuals.png"), g_res, width = 8, height = 5, dpi = 150)

  p_curve <- data.frame(
    dose_uM = sort(unique(pred$dose_uM)),
    dose_uM_plot = pmax(sort(unique(pred$dose_uM)), 1e-4),
    p_mis = hill_pmis(sort(unique(pred$dose_uM)), cfg$p_mis_base, best_fit$theta$p_mis_doxo_max, best_fit$theta$EC50, best_fit$theta$hill_n)
  )
  utils::write.table(p_curve, file.path(out_dir, "pmis_curve.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  g_p <- ggplot(p_curve, aes(dose_uM_plot, p_mis)) +
    geom_line(linewidth = 1.0, color = "#b44f00") +
    geom_point(size = 2, color = "#b44f00") +
    scale_x_log10() +
    labs(title = "Inferred p_mis(C) curve", x = "Dose (uM)", y = "Per-copy p_mis")
  ggsave(file.path(out_dir, "fit_pmis_curve.png"), g_p, width = 7, height = 4.5, dpi = 150)

  nulli_states <- c(46, 53, 90, 95)
  nulli_curve <- bind_rows(lapply(nulli_states, function(N) {
    pvals <- hill_pmis(sort(unique(pred$dose_uM)), cfg$p_mis_base, best_fit$theta$p_mis_doxo_max, best_fit$theta$EC50, best_fit$theta$hill_n)
    dose_vals <- sort(unique(pred$dose_uM))
    dead_prob <- vapply(pvals, function(p) {
      attr(.pr_delta_vec(N, p, gamma_loss = cfg$gamma_loss, N_unit = cfg$N_UNIT), "mass_dropped")
    }, numeric(1))
    data.frame(
      dose_uM = dose_vals,
      dose_uM_plot = pmax(dose_vals, 1e-4),
      dead_prob = dead_prob,
      N = as.factor(N),
      stringsAsFactors = FALSE
    )
  }))
  utils::write.table(nulli_curve, file.path(out_dir, "nullisomy_prob_by_dose.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  g_n <- ggplot(nulli_curve, aes(dose_uM_plot, dead_prob, color = N)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    labs(title = "Nullisomy/nonviable daughter probability across dose", x = "Dose (uM)", y = "Dropped daughter mass")
  ggsave(file.path(out_dir, "fit_nullisomy_prob_by_dose.png"), g_n, width = 8, height = 5, dpi = 150)

  by_sample_rmse <- residual %>%
    group_by(sample) %>%
    summarise(
      rmse_log = sqrt(mean(residual_log^2, na.rm = TRUE)),
      .groups = "drop"
    )
  high_dose_obs <- pred %>%
    group_by(sample) %>%
    filter(dose_uM == max(dose_uM, na.rm = TRUE)) %>%
    summarise(observed = observed[[1]], predicted = predicted[[1]], .groups = "drop")
  more_sensitive_sample <- high_dose_obs$sample[[which.min(high_dose_obs$observed)]]
  plausibility_note <- if (best_fit$theta$p_mis_doxo_max > 0.3) {
    "The fitted p_mis_doxo_max is very high and likely biologically implausible for a per-copy missegregation probability."
  } else {
    "The fitted p_mis_doxo_max stays below an obviously catastrophic per-copy regime."
  }
  adequacy_note <- if (max(by_sample_rmse$rmse_log, na.rm = TRUE) > 0.5) {
    "Shared-parameter nullisomy-only death does not adequately match both dose-response curves."
  } else {
    "Shared-parameter nullisomy-only death gives a tolerable first-pass fit."
  }

  report_lines <- c(
    "Doxorubicin-nullisomy fit report",
    "",
    paste("Objective:", signif(best_fit$objective, 6)),
    paste("p_mis_base:", signif(cfg$p_mis_base, 6)),
    paste("p_mis_doxo_max:", signif(best_fit$theta$p_mis_doxo_max, 6)),
    paste("EC50_uM:", signif(best_fit$theta$EC50, 6)),
    paste("Hill_n:", signif(best_fit$theta$hill_n, 6)),
    paste("lambda_eff:", signif(best_fit$theta$lambda_eff, 6)),
    paste("assay_days_fixed:", signif(cfg$assay_days_fixed, 6)),
    "",
    "Adequacy summary:",
    adequacy_note,
    plausibility_note,
    paste("Most sensitive observed sample at max dose:", more_sensitive_sample),
    "This fitter uses shared doxorubicin-nullisomy parameters across both SNU-668 samples.",
    "",
    "Sample-wise RMSE on log normalized viability:",
    paste(sprintf("%s: %.4f", by_sample_rmse$sample, by_sample_rmse$rmse_log), collapse = "; "),
    "",
    "Interpretation:",
    "This fit uses the active O2_supply_demand_MAP C++ missegregation/nullisomy kernel with p_const driven by the doxorubicin Hill law.",
    "Direct hypoxia death is set to zero by forcing mu_hp = 0 and treatment scaling is disabled with fit_treatment = FALSE.",
    "Predicted response is normalized viable burden at dose C divided by predicted viable burden at dose 0 for the same sample.",
    "If this mode remains inadequate, the next minimal extension is a direct doxorubicin-linked death or growth-suppression term in addition to nullisomy-only death."
  )
  writeLines(report_lines, con = file.path(out_dir, "technical_report.txt"))
}

run_fit <- function(data_obj, cfg, out_dir, seeds) {
  fits <- lapply(seeds, function(seed) fit_one_seed(seed, data_obj, cfg))
  summary_tab <- bind_rows(lapply(fits, function(x) {
    data.frame(
      seed = x$seed,
      objective = x$objective,
      convergence = x$convergence,
      p_mis_doxo_max = x$theta$p_mis_doxo_max,
      EC50_uM = x$theta$EC50,
      hill_n = x$theta$hill_n,
      lambda_eff = x$theta$lambda_eff,
      assay_days_fixed = cfg$assay_days_fixed,
      stringsAsFactors = FALSE
    )
  }))
  utils::write.table(summary_tab, file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  best_idx <- which.min(summary_tab$objective)
  best_fit <- fits[[best_idx]]
  write_fit_outputs(best_fit, data_obj, cfg, out_dir)
  invisible(best_fit)
}

main <- function(argv) {
  args <- parse_args(argv)
  mode <- tolower(trimws(as.character(if (is.null(args$mode)) "fit" else args$mode)))
  cfg <- default_cfg()
  if (!is.null(args$p_mis_base)) cfg$p_mis_base <- as_num(args$p_mis_base, cfg$p_mis_base)
  if (!is.null(args$gamma_loss)) cfg$gamma_loss <- as_num(args$gamma_loss, cfg$gamma_loss)
  if (!is.null(args$dt)) cfg$DT <- as_num(args$dt, cfg$DT)
  json_path <- if (!is.null(args$json)) args$json else default_data_path()
  out_dir <- make_output_dir(args$out_dir)
  data_obj <- load_doxo_data(json_path)
  diagnostic_curves(cfg, out_dir)
  if (mode == "diagnostics") return(invisible(TRUE))
  if (mode == "test") return(invisible(run_tests(data_obj, cfg, out_dir)))

  seeds <- if (!is.null(args$seeds_csv)) {
    as.integer(strsplit(args$seeds_csv, ",", fixed = TRUE)[[1]])
  } else {
    cfg$seeds_default
  }
  best_fit <- run_fit(data_obj, cfg, out_dir, seeds)
  message("Wrote doxorubicin-nullisomy outputs to: ", out_dir)
  invisible(best_fit)
}

if (sys.nframe() == 0L) {
  main(commandArgs(trailingOnly = TRUE))
}
