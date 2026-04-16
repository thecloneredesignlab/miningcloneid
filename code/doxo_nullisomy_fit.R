#!/usr/bin/env Rscript

# Exploratory standalone nullisomy-only fitter.
# This script reuses the active O2_supply_demand_MAP backend, but it is not yet
# integrated into the main optimizer/config path used by the primary invivo
# fitting workflow.

suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

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
    p_mis_cap = 0.95,
    p_mis_implausible_threshold = 0.3,
    N_2N_ref = 44,
    N_ref_2N_state = 46,
    N_ref_high_state = 95,
    high_ploidy_threshold = 90,
    assay_days_fixed = 4.0,
    optim_maxit = 12,
    fit_models_default = c("shared_per_copy", "continuous_ploidy_amplified", "categorical_high_ploidy"),
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

canonical_model_variant <- function(x) {
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("shared", "shared_per_copy", "null")) return("shared_per_copy")
  if (s %in% c("continuous", "continuous_ploidy_amplified", "alpha")) return("continuous_ploidy_amplified")
  if (s %in% c("categorical", "categorical_high_ploidy", "a4n")) return("categorical_high_ploidy")
  stop("Unsupported model variant: ", x)
}

model_variant_label <- function(x) {
  x <- canonical_model_variant(x)
  switch(
    x,
    shared_per_copy = "Shared per-copy",
    continuous_ploidy_amplified = "Continuous ploidy-amplified",
    categorical_high_ploidy = "Categorical high-ploidy amplified"
  )
}

variant_param_count <- function(variant) {
  variant <- canonical_model_variant(variant)
  if (identical(variant, "shared_per_copy")) 4L else 5L
}

theta_from_par <- function(par, cfg, variant) {
  variant <- canonical_model_variant(variant)
  p_headroom <- max(cfg$p_mis_cap - cfg$p_mis_base - 1e-8, 1e-8)
  out <- list(
    p_mis_doxo_max = p_headroom * plogis(par[[1]]),
    EC50 = 10^par[[2]],
    hill_n = 10^par[[3]],
    lambda_eff = 10^par[[4]]
  )
  if (identical(variant, "continuous_ploidy_amplified")) {
    out$alpha <- as.numeric(par[[5]])
  } else if (identical(variant, "categorical_high_ploidy")) {
    out$a_4N <- 10^par[[5]]
  }
  out
}

par_from_theta <- function(theta, cfg, variant) {
  variant <- canonical_model_variant(variant)
  p_headroom <- max(cfg$p_mis_cap - cfg$p_mis_base - 1e-8, 1e-8)
  frac <- clip01(theta$p_mis_doxo_max / p_headroom)
  out <- c(
    qlogis(pmin(pmax(frac, 1e-8), 1 - 1e-8)),
    log10(theta$EC50),
    log10(theta$hill_n),
    log10(theta$lambda_eff)
  )
  if (identical(variant, "continuous_ploidy_amplified")) {
    out <- c(out, as.numeric(theta$alpha))
  } else if (identical(variant, "categorical_high_ploidy")) {
    out <- c(out, log10(theta$a_4N))
  }
  out
}

default_theta <- function(cfg, variant) {
  variant <- canonical_model_variant(variant)
  out <- list(
    p_mis_doxo_max = 0.05,
    EC50 = 0.2,
    hill_n = 1.5,
    lambda_eff = 0.2
  )
  if (identical(variant, "continuous_ploidy_amplified")) {
    out$alpha <- 0.0
  } else if (identical(variant, "categorical_high_ploidy")) {
    out$a_4N <- 1.0
  }
  out
}

fit_bounds <- function(variant) {
  variant <- canonical_model_variant(variant)
  lower <- c(qlogis(1e-8), log10(1e-4), log10(0.2), log10(1e-3))
  upper <- c(qlogis(1 - 1e-8), log10(100), log10(8), log10(2))
  if (identical(variant, "continuous_ploidy_amplified")) {
    lower <- c(lower, -2.0)
    upper <- c(upper, 2.0)
  } else if (identical(variant, "categorical_high_ploidy")) {
    lower <- c(lower, log10(0.25))
    upper <- c(upper, log10(8.0))
  }
  list(lower = lower, upper = upper)
}

effective_pmis_by_state <- function(dose_uM, N, theta, cfg, variant) {
  variant <- canonical_model_variant(variant)
  N_use <- as.numeric(N)
  p_base <- hill_pmis(
    dose_uM = dose_uM,
    p_mis_base = cfg$p_mis_base,
    p_mis_doxo_max = theta$p_mis_doxo_max,
    EC50 = theta$EC50,
    hill_n = theta$hill_n
  )
  if (identical(variant, "shared_per_copy")) {
    mult <- rep(1.0, length(N_use))
  } else if (identical(variant, "continuous_ploidy_amplified")) {
    mult <- (pmax(N_use, 1) / cfg$N_2N_ref)^theta$alpha
  } else {
    mult <- ifelse(N_use >= cfg$high_ploidy_threshold, theta$a_4N, 1.0)
  }
  pmin(pmax(p_base * mult, 0), cfg$p_mis_cap)
}

build_state_generator <- function(theta, dose_uM, cfg, variant) {
  grid <- seq.int(cfg$N_MIN, cfg$N_MAX)
  p_vec <- effective_pmis_by_state(dose_uM, grid, theta, cfg, variant)
  G <- .build_G_with_WGD(
    N0min = cfg$N_MIN,
    N0max = cfg$N_MAX,
    lambda0_vec = rep(theta$lambda_eff, length(grid)),
    p0_vec = p_vec,
    wgd_prob_vec = rep(cfg$p_wgd, length(grid)),
    boundary = "drop",
    eps_tail = 1e-8,
    gamma_loss = cfg$gamma_loss,
    N_unit = cfg$N_UNIT
  )
  list(G = G, p_vec = p_vec, grid = grid)
}

simulate_one_doxo <- function(init_state_prob, dose_uM, theta, cfg, variant) {
  grid <- seq.int(cfg$N_MIN, cfg$N_MAX)
  init_total <- cfg$o2_Nref
  init_state <- as.numeric(init_state_prob) * init_total
  run_params <- make_run_params(theta, cfg)
  obs_steps <- as.integer(round(cfg$assay_days_fixed / cfg$DT))
  vol_by_N <- as.numeric(cell_volume_mm3_by_N(grid, run_params = run_params, cfg = cfg))
  built <- build_state_generator(theta, dose_uM, cfg, variant)
  x_live <- as.numeric(init_state)
  for (step in seq_len(obs_steps)) {
    x_live <- step_dt(built$G, x_live, cfg$DT, steps = 1L, normalize = FALSE)
    x_live[!is.finite(x_live) | x_live < 0] <- 0
  }
  p_ref_2N <- effective_pmis_by_state(dose_uM, cfg$N_ref_2N_state, theta, cfg, variant)[[1]]
  p_ref_high <- effective_pmis_by_state(dose_uM, cfg$N_ref_high_state, theta, cfg, variant)[[1]]
  list(
    p_const = hill_pmis(dose_uM, cfg$p_mis_base, theta$p_mis_doxo_max, theta$EC50, theta$hill_n),
    p_ref_2N = p_ref_2N,
    p_ref_high = p_ref_high,
    expected_mis_2N = expected_mis_copies(cfg$N_ref_2N_state, p_ref_2N),
    expected_mis_high = expected_mis_copies(cfg$N_ref_high_state, p_ref_high),
    live_burden = sum(x_live * vol_by_N),
    frac_N = x_live / max(sum(x_live), cfg$burden_log_eps)
  )
}

predict_dataset <- function(theta, data_obj, cfg, variant) {
  variant <- canonical_model_variant(variant)
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
    sims <- lapply(sub$dose_uM, function(dose) simulate_one_doxo(init_states[[sample_name]], dose, theta, cfg, variant))
    base_live <- sims[[which.min(abs(sub$dose_uM - 0))]]$live_burden
    base_live <- max(base_live, cfg$burden_log_eps)
    for (i in seq_len(nrow(sub))) {
      rows[[length(rows) + 1L]] <- data.frame(
        model = variant,
        model_label = model_variant_label(variant),
        sample = sample_name,
        dose_uM = sub$dose_uM[[i]],
        observed = sub$mean_normalized_to_0uM[[i]],
        observed_sd = sub$sd_normalized_to_0uM[[i]],
        predicted = sims[[i]]$live_burden / base_live,
        p_const = sims[[i]]$p_const,
        p_eff_2N = sims[[i]]$p_ref_2N,
        p_eff_high = sims[[i]]$p_ref_high,
        expected_mis_2N = sims[[i]]$expected_mis_2N,
        expected_mis_high = sims[[i]]$expected_mis_high,
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

fit_one_seed <- function(seed, data_obj, cfg, variant) {
  variant <- canonical_model_variant(variant)
  set.seed(seed)
  base <- default_theta(cfg, variant)
  jittered <- modifyList(base, list(
    p_mis_doxo_max = clip01(base$p_mis_doxo_max * 10^runif(1, -0.5, 0.5)),
    EC50 = base$EC50 * 10^runif(1, -0.75, 0.75),
    hill_n = base$hill_n * 10^runif(1, -0.5, 0.5),
    lambda_eff = base$lambda_eff * 10^runif(1, -0.4, 0.4)
  ))
  if (identical(variant, "continuous_ploidy_amplified")) {
    jittered$alpha <- runif(1, -0.2, 0.2)
  } else if (identical(variant, "categorical_high_ploidy")) {
    jittered$a_4N <- 10^runif(1, -0.2, 0.2)
  }
  init_par <- par_from_theta(jittered, cfg, variant)
  bounds <- fit_bounds(variant)
  fn <- function(par) {
    theta <- theta_from_par(par, cfg, variant)
    pred <- predict_dataset(theta, data_obj, cfg, variant)
    objective_from_table(pred, cfg)
  }
  opt <- optim(
    par = init_par,
    fn = fn,
    method = "L-BFGS-B",
    lower = bounds$lower,
    upper = bounds$upper,
    control = list(maxit = cfg$optim_maxit)
  )
  theta_hat <- theta_from_par(opt$par, cfg, variant)
  pred_hat <- predict_dataset(theta_hat, data_obj, cfg, variant)
  list(
    model = variant,
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

  theta <- default_theta(cfg, "shared_per_copy")
  init_r2 <- sample_initial_state(data_obj$cell_ploidy, "SNU-668_r2_GFP_VT_A7_harvesT3", cfg)
  sim_shared <- simulate_one_doxo(init_r2, dose_uM = 1.0, theta = theta, cfg = cfg, variant = "shared_per_copy")
  theta_alpha <- modifyList(default_theta(cfg, "continuous_ploidy_amplified"), list(alpha = 0.5))
  sim_alpha <- simulate_one_doxo(init_r2, dose_uM = 1.0, theta = theta_alpha, cfg = cfg, variant = "continuous_ploidy_amplified")
  stopifnot(sim_alpha$p_ref_high >= sim_shared$p_ref_high)

  writeLines("All doxorubicin-nullisomy checks passed.", con = file.path(out_dir, "tests_ok.txt"))
  invisible(TRUE)
}

theta_to_named_rows <- function(theta, variant) {
  variant <- canonical_model_variant(variant)
  rows <- data.frame(
    parameter = c("p_mis_doxo_max", "EC50_uM", "hill_n", "lambda_eff"),
    value = c(theta$p_mis_doxo_max, theta$EC50, theta$hill_n, theta$lambda_eff),
    stringsAsFactors = FALSE
  )
  if (identical(variant, "continuous_ploidy_amplified")) {
    rows <- bind_rows(rows, data.frame(parameter = "alpha", value = theta$alpha, stringsAsFactors = FALSE))
  } else if (identical(variant, "categorical_high_ploidy")) {
    rows <- bind_rows(rows, data.frame(parameter = "a_4N", value = theta$a_4N, stringsAsFactors = FALSE))
  }
  rows
}

implausibility_note_for_predictions <- function(pred_tab, cfg) {
  max_p <- max(c(pred_tab$p_eff_2N, pred_tab$p_eff_high), na.rm = TRUE)
  if (!is.finite(max_p)) return("Could not evaluate p_mis plausibility.")
  if (max_p > cfg$p_mis_implausible_threshold) {
    paste0("Effective per-copy p_mis reaches ", signif(max_p, 4), ", which is likely biologically implausible.")
  } else {
    paste0("Effective per-copy p_mis stays below ", signif(cfg$p_mis_implausible_threshold, 4), ".")
  }
}

assess_identifiability <- function(fit_rows, variant) {
  variant <- canonical_model_variant(variant)
  if (identical(variant, "shared_per_copy")) return("Not applicable.")
  conv <- fit_rows[fit_rows$convergence == 0, , drop = FALSE]
  if (nrow(conv) < 2L) return("Insufficient converged seeds to assess identifiability.")
  target <- if (identical(variant, "continuous_ploidy_amplified")) conv$alpha else conv$a_4N
  target <- as.numeric(target)
  if (any(!is.finite(target))) return("Amplification parameter unavailable.")
  if (sd(target) > max(0.15, 0.5 * abs(mean(target)))) {
    return("Weakly identifiable across seeds.")
  }
  "Reasonably stable across seeds."
}

write_fit_outputs <- function(best_fits, fit_summary, data_obj, cfg, out_dir) {
  pred_all <- bind_rows(lapply(best_fits, `[[`, "predictions"))
  pred_all$dose_uM_plot <- pmax(pred_all$dose_uM, 1e-4)
  residual <- transform(
    pred_all,
    residual_log = log(pmax(predicted, cfg$burden_log_eps)) - log(pmax(observed, cfg$burden_log_eps))
  )

  params_all <- bind_rows(lapply(best_fits, function(fit) {
    bind_rows(
      data.frame(parameter = "p_mis_base", value = cfg$p_mis_base, stringsAsFactors = FALSE),
      theta_to_named_rows(fit$theta, fit$model),
      data.frame(parameter = "assay_days_fixed", value = cfg$assay_days_fixed, stringsAsFactors = FALSE),
      data.frame(parameter = "objective", value = fit$objective, stringsAsFactors = FALSE)
    ) %>%
      mutate(model = fit$model, model_label = model_variant_label(fit$model), .before = 1)
  }))
  utils::write.table(params_all, file.path(out_dir, "best_params_by_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(residual, file.path(out_dir, "predictions_all_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  fit_summary$aic <- 2 * fit_summary$objective + 2 * fit_summary$n_params
  utils::write.table(fit_summary, file.path(out_dir, "model_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  g_fit <- ggplot(pred_all, aes(dose_uM_plot, observed, color = sample)) +
    geom_point(size = 1.8) +
    geom_line(aes(y = predicted), linewidth = 0.9) +
    geom_errorbar(aes(ymin = pmax(observed - observed_sd, 1e-6), ymax = observed + observed_sd), width = 0.02, alpha = 0.5) +
    scale_x_log10() +
    facet_wrap(~ model_label) +
    labs(title = "Observed vs predicted normalized doxorubicin response", x = "Dose (uM)", y = "Normalized viable burden")
  ggsave(file.path(out_dir, "fit_observed_vs_predicted_all_models.png"), g_fit, width = 11, height = 7, dpi = 150)

  g_res <- ggplot(residual, aes(dose_uM_plot, residual_log, color = sample)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 1.8) +
    scale_x_log10() +
    facet_wrap(~ model_label) +
    labs(title = "Residuals by model, sample, and dose", x = "Dose (uM)", y = "log(pred) - log(obs)")
  ggsave(file.path(out_dir, "fit_residuals_all_models.png"), g_res, width = 11, height = 7, dpi = 150)

  p_curve <- pred_all %>%
    distinct(model, model_label, dose_uM, p_eff_2N, p_eff_high) %>%
    mutate(dose_uM_plot = pmax(dose_uM, 1e-4)) %>%
    tidyr::pivot_longer(c(p_eff_2N, p_eff_high), names_to = "state", values_to = "p_eff")
  utils::write.table(p_curve, file.path(out_dir, "effective_pmis_by_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  g_p <- ggplot(p_curve, aes(dose_uM_plot, p_eff, color = state)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    facet_wrap(~ model_label) +
    labs(title = "Effective per-copy p_mis by model", x = "Dose (uM)", y = "Per-copy p_mis")
  ggsave(file.path(out_dir, "effective_pmis_by_model.png"), g_p, width = 11, height = 7, dpi = 150)

  mis_curve <- pred_all %>%
    distinct(model, model_label, dose_uM, expected_mis_2N, expected_mis_high) %>%
    mutate(dose_uM_plot = pmax(dose_uM, 1e-4)) %>%
    tidyr::pivot_longer(c(expected_mis_2N, expected_mis_high), names_to = "state", values_to = "expected_mis")
  utils::write.table(mis_curve, file.path(out_dir, "expected_missegregated_chromosomes_by_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  g_m <- ggplot(mis_curve, aes(dose_uM_plot, expected_mis, color = state)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    facet_wrap(~ model_label) +
    labs(title = "Expected missegregated chromosomes per division", x = "Dose (uM)", y = "Expected missegregated chromosomes")
  ggsave(file.path(out_dir, "expected_missegregated_chromosomes_by_model.png"), g_m, width = 11, height = 7, dpi = 150)

  report_lines <- c("Doxorubicin-nullisomy model comparison", "")
  for (fit in best_fits) {
    variant <- fit$model
    fit_tab <- fit_summary[fit_summary$model == variant, , drop = FALSE]
    pred_tab <- fit$predictions
    report_lines <- c(
      report_lines,
      paste0("[", model_variant_label(variant), "]"),
      paste("Objective:", signif(fit$objective, 6)),
      paste("AIC:", signif(fit_tab$aic[[1]], 6)),
      paste("Parameters:", paste(paste(theta_to_named_rows(fit$theta, variant)$parameter, signif(theta_to_named_rows(fit$theta, variant)$value, 6), sep = "="), collapse = ", ")),
      paste("Identifiability:", assess_identifiability(fit_summary[fit_summary$model == variant, , drop = FALSE], variant)),
      implausibility_note_for_predictions(pred_tab, cfg),
      ""
    )
  }
  best_variant <- fit_summary$model[[which.min(fit_summary$aic)]]
  report_lines <- c(
    report_lines,
    paste("Best AIC model:", model_variant_label(best_variant)),
    "Amplification was tested without adding any direct death term."
  )
  writeLines(report_lines, con = file.path(out_dir, "technical_report.txt"))
}

run_fit <- function(data_obj, cfg, out_dir, seeds, variants = NULL) {
  if (is.null(variants)) variants <- cfg$fit_models_default
  variants <- vapply(variants, canonical_model_variant, character(1))
  all_fits <- unlist(
    lapply(variants, function(variant) lapply(seeds, function(seed) fit_one_seed(seed, data_obj, cfg, variant))),
    recursive = FALSE
  )
  fit_summary <- bind_rows(lapply(all_fits, function(x) {
    row <- data.frame(
      model = x$model,
      model_label = model_variant_label(x$model),
      seed = x$seed,
      objective = x$objective,
      convergence = x$convergence,
      n_params = variant_param_count(x$model),
      p_mis_doxo_max = x$theta$p_mis_doxo_max,
      EC50_uM = x$theta$EC50,
      hill_n = x$theta$hill_n,
      lambda_eff = x$theta$lambda_eff,
      assay_days_fixed = cfg$assay_days_fixed,
      stringsAsFactors = FALSE
    )
    if (!is.null(x$theta$alpha)) row$alpha <- x$theta$alpha
    if (!is.null(x$theta$a_4N)) row$a_4N <- x$theta$a_4N
    row
  }))
  utils::write.table(fit_summary, file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  best_fits <- lapply(variants, function(variant) {
    idx <- which(fit_summary$model == variant)
    all_fits[[idx[which.min(fit_summary$objective[idx])]]]
  })
  names(best_fits) <- variants
  write_fit_outputs(best_fits, fit_summary, data_obj, cfg, out_dir)
  invisible(list(best_fits = best_fits, fit_summary = fit_summary))
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
