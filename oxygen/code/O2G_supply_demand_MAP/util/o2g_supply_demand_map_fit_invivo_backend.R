#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readxl))

.o2sd_bootstrap_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
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
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2sd_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool
clip <- o2sd_clip
get_script_dir <- o2sd_get_script_dir

# -----------------------------------------------------------------------------
# Function: require_cli_args
# Purpose: Enforce that required runtime parameters are explicitly provided.
# Parameters:
#   - argv: Parsed CLI key-value list.
#   - keys: Character vector of required argument keys.
# Returns:
#   NULL (errors on missing keys).
# -----------------------------------------------------------------------------
require_cli_args <- function(argv, keys) {
  missing <- keys[vapply(
    keys,
    function(k) {
      v <- argv[[k]]
      is.null(v) || !nzchar(trimws(as.character(v)))
    },
    logical(1)
  )]
  if (length(missing) > 0L) {
    stop(
      "Missing required CLI args (must be provided via YAML runner): ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Function: .sample_uniform_box
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - n: Requested number of samples, points, or steps.
#   - lower: Lower bounds for optimization variables.
#   - upper: Upper bounds for optimization variables.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.sample_uniform_box <- function(n, lower, upper) {
  n <- as.integer(n)
  d <- length(lower)
  if (is.na(n) || n <= 0L || d <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = d))
  }
  u <- matrix(stats::runif(n * d), nrow = n, ncol = d)
  span <- as.numeric(upper - lower)
  out <- sweep(u, 2, span, `*`)
  out <- sweep(out, 2, as.numeric(lower), `+`)
  colnames(out) <- names(lower)
  out
}

# -----------------------------------------------------------------------------
# Function: .sample_truncnorm_box
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - n: Requested number of samples, points, or steps.
#   - center: Function-specific input argument.
#   - lower: Lower bounds for optimization variables.
#   - upper: Upper bounds for optimization variables.
#   - sigma_frac: Relative local sampling width around warm start in transformed scale.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.sample_truncnorm_box <- function(n, center, lower, upper, sigma_frac = 0.1) {
  n <- as.integer(n)
  d <- length(lower)
  if (is.na(n) || n <= 0L || d <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = d))
  }
  sigma_frac <- as.numeric(sigma_frac)
  if (!is.finite(sigma_frac) || sigma_frac <= 0) sigma_frac <- 0.1

  sd_vec <- pmax((as.numeric(upper) - as.numeric(lower)) * sigma_frac, 1e-12)
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  out <- sweep(z, 2, sd_vec, `*`)
  out <- sweep(out, 2, as.numeric(center), `+`)
  out <- sweep(out, 2, as.numeric(lower), pmax)
  out <- sweep(out, 2, as.numeric(upper), pmin)
  colnames(out) <- names(lower)
  out
}

prob_to_logit <- function(p, param_symbol = "probability", slot_label = "value") {
  p_num <- as.numeric(p)
  if (!is.finite(p_num) || p_num <= 0 || p_num >= 1) {
    stop(param_symbol, " ", slot_label, " must be strictly within (0,1) for logit transform.")
  }
  qlogis(p_num)
}

logit_to_prob <- function(z) {
  p <- plogis(as.numeric(z))
  pmin(pmax(p, 1e-12), 1 - 1e-12)
}

# -----------------------------------------------------------------------------
# Function: .build_de_initialpop
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - np: Population size for DEoptim.
#   - lower: Lower bounds for optimization variables.
#   - upper: Upper bounds for optimization variables.
#   - init_use: Optional warm-start parameter vector used to seed optimization.
#   - mode: Initialization or aggregation mode selector.
#   - uniform_frac: Fraction of DE population initialized from global uniform sampling.
#   - sigma_frac: Relative local sampling width around warm start in transformed scale.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_de_initialpop <- function(np, lower, upper, init_use = NULL, mode = "hybrid", uniform_frac = 0.3, sigma_frac = 0.1) {
  np <- as.integer(np)
  if (is.na(np) || np < 1L) stop("DEoptim NP must be >= 1.")
  mode <- tolower(trimws(as.character(.first_non_null_local(mode, "hybrid"))))
  if (!mode %in% c("hybrid", "uniform")) {
    stop("de_init_mode must be one of: hybrid, uniform.")
  }
  uniform_frac <- as.numeric(uniform_frac)
  if (!is.finite(uniform_frac) || uniform_frac < 0 || uniform_frac > 1) {
    stop("de_init_uniform_frac must be in [0,1].")
  }
  sigma_frac <- as.numeric(sigma_frac)
  if (!is.finite(sigma_frac) || sigma_frac <= 0) {
    stop("de_init_sigma_frac must be > 0.")
  }

  pop <- .sample_uniform_box(np, lower, upper)
  warm_start_used <- !is.null(init_use)
  if (!warm_start_used) {
    return(list(
      pop = pop,
      mode_effective = "uniform_no_warm_start",
      warm_start_used = FALSE,
      n_local = 0L,
      n_uniform = as.integer(np)
    ))
  }

  init_vec <- as.numeric(init_use[names(lower)])
  names(init_vec) <- names(lower)
  if (any(!is.finite(init_vec))) stop("Warm-start vector has missing/non-finite values.")
  init_vec <- clip(init_vec, lower, upper)
  pop[1, ] <- init_vec

  rem <- as.integer(np - 1L)
  if (rem <= 0L) {
    return(list(
      pop = pop,
      mode_effective = mode,
      warm_start_used = TRUE,
      n_local = 0L,
      n_uniform = 0L
    ))
  }

  if (mode == "uniform") {
    if (rem > 0L) {
      pop[2:np, ] <- .sample_uniform_box(rem, lower, upper)
    }
    return(list(
      pop = pop,
      mode_effective = "uniform",
      warm_start_used = TRUE,
      n_local = 0L,
      n_uniform = rem
    ))
  }

  n_uniform <- as.integer(round(rem * uniform_frac))
  n_uniform <- max(0L, min(rem, n_uniform))
  n_local <- as.integer(rem - n_uniform)
  cursor <- 2L
  if (n_local > 0L) {
    idx <- cursor:(cursor + n_local - 1L)
    pop[idx, ] <- .sample_truncnorm_box(
      n = n_local,
      center = init_vec,
      lower = lower,
      upper = upper,
      sigma_frac = sigma_frac
    )
    cursor <- cursor + n_local
  }
  if (n_uniform > 0L) {
    idx <- cursor:(cursor + n_uniform - 1L)
    pop[idx, ] <- .sample_uniform_box(n_uniform, lower, upper)
  }
  if (rem > 1L) {
    ord <- sample.int(rem, size = rem, replace = FALSE)
    pop[2:np, ] <- pop[1L + ord, , drop = FALSE]
  }
  list(
    pop = pop,
    mode_effective = "hybrid",
    warm_start_used = TRUE,
    n_local = n_local,
    n_uniform = n_uniform
  )
}

.first_non_null_local <- o2sd_first_non_null

# -----------------------------------------------------------------------------
# Function: resolve_glucose_settings_local
# Purpose: Resolve runtime glucose mode plus parameter defaults from cfg/run_params.
# -----------------------------------------------------------------------------
resolve_glucose_settings_local <- function(run_params = NULL, cfg = NULL) {
  rp <- if (is.null(run_params)) list() else as.list(run_params)
  cfg_use <- if (is.null(cfg)) list() else cfg
  glucose_defaults <- default_glucose_pct_scale()
  glucose_use <- canonical_glucose_enabled(
    .first_non_null_local(rp$glucose, cfg_use$glucose, TRUE),
    default = TRUE
  )

  glucose_dynamic_use <- canonical_glucose_dynamic(
    .first_non_null_local(rp$glucose_dynamic, cfg_use$glucose_dynamic, FALSE),
    default = FALSE
  )
  if (!isTRUE(glucose_use)) {
    glucose_dynamic_use <- FALSE
  }
  glucose_stress_mode_use <- resolve_glucose_runtime_mode(
    glucose_dynamic = glucose_dynamic_use,
    glucose_stress_mode = if (isTRUE(glucose_use)) {
      .first_non_null_local(rp$glucose_stress_mode, cfg_use$glucose_stress_mode, "coupled_to_O2")
    } else {
      "off"
    },
    default_dynamic = glucose_dynamic_use,
    default_static_mode = if (isTRUE(glucose_use)) "coupled_to_O2" else "off"
  )

  glucose_ref_mM_use <- normalize_glucose_ref_mM(
    .first_non_null_local(rp$glucose_ref_mM, cfg_use$glucose_ref_mM, default_glucose_ref_mM())
  )
  G_S0_use <- as.numeric(.first_non_null_local(rp$G_S0, cfg_use$G_S0_init, glucose_defaults$G_S0))
  kappa_G_use <- as.numeric(.first_non_null_local(rp$kappa_G, cfg_use$kappa_G_init, glucose_defaults$kappa_G))
  eta_G_use <- as.numeric(.first_non_null_local(rp$eta_G, cfg_use$eta_G_init, glucose_defaults$eta_G))
  G_c_use <- as.numeric(.first_non_null_local(rp$G_c, cfg_use$G_c_init, glucose_defaults$G_c))
  tau_G_use <- as.numeric(.first_non_null_local(rp$tau_G, cfg_use$tau_G, cfg_use$tau_G_init, cfg_use$tau_O2, cfg_use$tau_O2_init, 0.1))

  if (!is.finite(G_S0_use) || G_S0_use <= 0) G_S0_use <- glucose_defaults$G_S0
  if (!is.finite(kappa_G_use) || kappa_G_use <= 0) kappa_G_use <- glucose_defaults$kappa_G
  if (!is.finite(eta_G_use) || eta_G_use <= 0) eta_G_use <- glucose_defaults$eta_G
  if (!is.finite(G_c_use) || G_c_use <= 0) G_c_use <- glucose_defaults$G_c
  if (!is.finite(tau_G_use) || tau_G_use <= 0) tau_G_use <- 0.1

  list(
    glucose = glucose_use,
    glucose_dynamic = glucose_dynamic_use,
    glucose_stress_mode = glucose_stress_mode_use,
    glucose_ref_mM = glucose_ref_mM_use,
    G_S0 = G_S0_use,
    kappa_G = kappa_G_use,
    eta_G = eta_G_use,
    G_c = G_c_use,
    tau_G = tau_G_use
  )
}

resolve_harvest_init_settings_local <- function(run_params = NULL, cfg = NULL, harvest_ids = NULL) {
  rp <- if (is.null(run_params)) list() else as.list(run_params)
  cfg_use <- if (is.null(cfg)) list() else cfg
  enabled <- o2sd_as_bool_scalar(
    .first_non_null_local(rp$harvest_init_multiplier, cfg_use$harvest_init_multiplier, FALSE),
    FALSE
  )
  ids <- .first_non_null_local(harvest_ids, cfg_use$harvest_param_ids, character(0))
  ids <- as.character(ids)
  ids <- ids[nzchar(ids)]
  ids <- unique(ids)
  list(
    enabled = isTRUE(enabled),
    harvest_ids = ids,
    log_param_names = setNames(vapply(ids, harvest_init_log_param_name, character(1)), ids),
    natural_param_names = setNames(vapply(ids, harvest_init_natural_param_name, character(1)), ids),
    prior_center = as.numeric(.first_non_null_local(cfg_use$prior_center_log_init_mult, 0.0)),
    prior_sd = as.numeric(.first_non_null_local(cfg_use$prior_sd_log_init_mult, 0.35)),
    lower = as.numeric(.first_non_null_local(cfg_use$log_init_mult_lower, -1.0)),
    upper = as.numeric(.first_non_null_local(cfg_use$log_init_mult_upper, 1.0))
  )
}

scenario_init_multiplier_local <- function(run_params, scenario, cfg) {
  hs <- resolve_harvest_init_settings_local(run_params = run_params, cfg = cfg)
  if (!isTRUE(hs$enabled)) return(1.0)
  natural_name <- .first_non_null_local(scenario$init_mult_param, harvest_init_natural_param_name(scenario$harvest))
  log_name <- .first_non_null_local(scenario$log_init_mult_param, harvest_init_log_param_name(scenario$harvest))
  mult <- suppressWarnings(as.numeric(.first_non_null_local(run_params[[natural_name]], NA_real_)))
  if (!is.finite(mult) || mult <= 0) {
    log_mult <- suppressWarnings(as.numeric(.first_non_null_local(run_params[[log_name]], 0.0)))
    if (!is.finite(log_mult)) log_mult <- 0.0
    mult <- exp(log_mult)
  }
  if (!is.finite(mult) || mult <= 0) mult <- 1.0
  mult
}

# -----------------------------------------------------------------------------
# Function: get_param_names
# Purpose: Return ordered parameter names used in transformed optimization vectors.
# Parameters:
#   - fit_treatment: Logical flag indicating whether treatment-effect parameters are optimized.
#   - fit_tau_O2: Logical flag indicating whether tau_O2 is estimated instead of fixed.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
get_param_names <- function(fit_treatment = TRUE,
                            fit_tau_O2 = FALSE,
                            glucose = TRUE,
                            glucose_dynamic = FALSE,
                            misseg_loss_survival = "nullisomy",
                            harvest_init_multiplier = FALSE,
                            harvest_ids = NULL) {
  glucose_use <- isTRUE(glucose)
  glucose_dynamic_use <- isTRUE(glucose_dynamic) && glucose_use
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(misseg_loss_survival, "nullisomy"),
    default = "nullisomy"
  )
  harvest_init_use <- isTRUE(harvest_init_multiplier)
  harvest_ids_use <- unique(as.character(.first_non_null_local(harvest_ids, character(0))))
  harvest_ids_use <- harvest_ids_use[nzchar(harvest_ids_use)]
  nm <- c(
    "log10_lam_min",
    "delta_lam",
    "log10_p_mis_base",
    "logit_p_misseg",
    "log10_k_o_mis"
  )
  if (identical(loss_mode, "nullisomy")) {
    nm <- c(nm, "log10_gamma_loss")
  } else {
    nm <- c(nm, "buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp")
  }
  if (!glucose_use || !glucose_dynamic_use) {
    nm <- append(nm, "log10_k_o", after = 2L)
  }
  nm <- c(nm, "logit_p_wgd")
  if (glucose_use && glucose_dynamic_use) {
    nm <- c(nm, "log10_G_S0", "log10_kappa_G", "log10_eta_G", "log10_G_c", "log10_tau_G")
  }
  nm <- c(
    nm,
    "log10_o2_S0",
    "log10_kappa_O",
    "log10_eta_o2",
    "log10_rho_2N",
    "log10_alpha_o2",
    "gamma_growth",
    "log10_mu_hp",
    "gamma_mu",
    "log10_O2_crit",
    "n_O",
    "log10_k_clear",
    "log10_sigma_burden"
  )
  if (isTRUE(fit_tau_O2)) {
    nm <- c(nm, "log10_tau_O2")
  }
  if (isTRUE(fit_treatment)) {
    nm <- c(nm, "log10_alpha", "gamma")
  }
  if (harvest_init_use && length(harvest_ids_use) > 0L) {
    nm <- c(nm, vapply(harvest_ids_use, harvest_init_log_param_name, character(1)))
  }
  nm
}

# -----------------------------------------------------------------------------
# Function: compute_soft_prior_penalty
# Purpose: Compute optional soft-prior penalty term from transformed parameters.
# Parameters:
#   - par_transformed: Parameter vector in transformed optimization scale.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
compute_soft_prior_penalty <- function(par_transformed, cfg) {
  if (!isTRUE(.first_non_null_local(cfg$use_soft_prior, FALSE))) {
    return(list(L_prior_raw = 0, n_terms = 0, terms = numeric(0)))
  }
  lambda <- as.numeric(.first_non_null_local(cfg$lambda_prior, 0))
  if (!is.finite(lambda) || lambda <= 0) {
    return(list(L_prior_raw = 0, n_terms = 0, terms = numeric(0)))
  }

  p <- as.numeric(par_transformed)
  p_names <- names(par_transformed)
  if (is.null(p_names) || length(p_names) != length(p)) {
    p_names <- get_param_names(
      fit_treatment = isTRUE(cfg$fit_treatment),
      fit_tau_O2 = isTRUE(.first_non_null_local(cfg$fit_tau_O2, FALSE)),
      glucose = isTRUE(.first_non_null_local(cfg$glucose, TRUE)),
      glucose_dynamic = isTRUE(.first_non_null_local(cfg$glucose_dynamic, FALSE)),
      misseg_loss_survival = .first_non_null_local(cfg$misseg_loss_survival, "nullisomy"),
      harvest_init_multiplier = isTRUE(.first_non_null_local(cfg$harvest_init_multiplier, FALSE)),
      harvest_ids = .first_non_null_local(cfg$harvest_param_ids, character(0))
    )
    if (length(p_names) != length(p)) p_names <- rep("", length(p))
  }
  names(p) <- p_names

  centers <- c(
    log10_k_o = as.numeric(cfg$prior_center_log10_k_o),
    log10_kappa_O = as.numeric(cfg$prior_center_log10_kappa_O),
    log10_o2_S0 = as.numeric(cfg$prior_center_log10_o2_S0),
    log10_eta_o2 = as.numeric(cfg$prior_center_log10_eta_o2),
    log10_gamma_loss = as.numeric(cfg$prior_center_log10_gamma_loss),
    buffer_smax = as.numeric(.first_non_null_local(cfg$prior_center_buffer_smax, cfg$buffer_smax_init, 0.8)),
    log10_buffer_beta = as.numeric(.first_non_null_local(cfg$prior_center_log10_buffer_beta, log10(.first_non_null_local(cfg$buffer_beta_init, 1.0)))),
    log10_buffer_n_exp = as.numeric(.first_non_null_local(cfg$prior_center_log10_buffer_n_exp, log10(.first_non_null_local(cfg$buffer_n_exp_init, 1.0)))),
    log10_rho_2N = as.numeric(cfg$prior_center_log10_rho_2N),
    log10_mu_hp = as.numeric(cfg$prior_center_log10_mu_hp),
    gamma_mu = as.numeric(cfg$prior_center_gamma_mu),
    log10_k_clear = as.numeric(cfg$prior_center_log10_k_clear)
  )
  sds <- c(
    log10_k_o = as.numeric(cfg$prior_sd_log10_k_o),
    log10_kappa_O = as.numeric(cfg$prior_sd_log10_kappa_O),
    log10_o2_S0 = as.numeric(cfg$prior_sd_log10_o2_S0),
    log10_eta_o2 = as.numeric(cfg$prior_sd_log10_eta_o2),
    log10_gamma_loss = as.numeric(cfg$prior_sd_log10_gamma_loss),
    buffer_smax = as.numeric(.first_non_null_local(cfg$prior_sd_buffer_smax, 0.25)),
    log10_buffer_beta = as.numeric(.first_non_null_local(cfg$prior_sd_log10_buffer_beta, 0.75)),
    log10_buffer_n_exp = as.numeric(.first_non_null_local(cfg$prior_sd_log10_buffer_n_exp, 0.75)),
    log10_rho_2N = as.numeric(cfg$prior_sd_log10_rho_2N),
    log10_mu_hp = as.numeric(cfg$prior_sd_log10_mu_hp),
    gamma_mu = as.numeric(cfg$prior_sd_gamma_mu),
    log10_k_clear = as.numeric(cfg$prior_sd_log10_k_clear)
  )

  shared <- intersect(intersect(names(p), names(centers)), names(sds))
  if (length(shared) == 0L) {
    return(list(L_prior_raw = 0, n_terms = 0, terms = numeric(0)))
  }

  terms <- setNames(numeric(0), character(0))
  for (nm in shared) {
    mu <- centers[[nm]]
    sdv <- sds[[nm]]
    pv <- p[[nm]]
    if (!is.finite(mu) || !is.finite(sdv) || sdv <= 0 || !is.finite(pv)) next
    z <- (pv - mu) / sdv
    terms[[nm]] <- 0.5 * z^2
  }
  harvest_settings <- resolve_harvest_init_settings_local(cfg = cfg)
  if (isTRUE(harvest_settings$enabled) && is.finite(harvest_settings$prior_sd) && harvest_settings$prior_sd > 0) {
    init_names <- grep("^log_init_mult_", names(p), value = TRUE)
    for (nm in init_names) {
      pv <- p[[nm]]
      if (!is.finite(pv)) next
      z <- (pv - harvest_settings$prior_center) / harvest_settings$prior_sd
      terms[[nm]] <- 0.5 * z^2
    }
  }
  if (length(terms) == 0L) {
    return(list(L_prior_raw = 0, n_terms = 0, terms = numeric(0)))
  }
  list(L_prior_raw = sum(terms), n_terms = length(terms), terms = terms)
}

# -----------------------------------------------------------------------------
# Function: burden_volume_mm3_from_state
# Purpose: Convert state vector to predicted burden volume in mm^3.
# Parameters:
#   - v: Function-specific input argument.
#   - grid_pre: Ploidy grid.
#   - R0: Number of ploidy states.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - vol_by_N: Optional precomputed per-state cell volume lookup.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
burden_volume_mm3_from_state <- function(v, grid_pre, R0, run_params, cfg, vol_by_N = NULL) {
  if (length(v) != R0) stop("State length mismatch in burden_volume_mm3_from_state.")
  if (is.null(vol_by_N)) {
    vol_by_N <- cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg)
  }
  counts_N <- as.numeric(v)
  sum(as.numeric(counts_N) * as.numeric(vol_by_N), na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# Function: default_n_cores
# Purpose: Infer a safe default worker count from available CPU cores.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
default_n_cores <- function() {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(n) || is.na(n)) {
    n <- suppressWarnings(parallel::detectCores())
  }
  if (!is.finite(n) || is.na(n)) return(1L)
  as.integer(max(1L, n - 1L))
}

# -----------------------------------------------------------------------------
# Function: map_scenarios_parallel
# Purpose: Apply a function to scenarios with optional parallel execution.
# Parameters:
#   - scenarios: List of scenario-specific observation data and metadata.
#   - n_cores: Requested number of CPU workers.
#   - label: Text label used for logging and progress messages.
#   - fn: Function applied to each scenario or element.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
map_scenarios_parallel <- function(scenarios, n_cores = 1L, label = "predict", fn) {
  n <- length(scenarios)
  if (n == 0L) return(list())
  n_cores_use <- suppressWarnings(as.integer(n_cores))
  if (!is.finite(n_cores_use) || n_cores_use < 1L) n_cores_use <- 1L
  workers <- as.integer(max(1L, min(n_cores_use, n)))
  use_mc <- (workers > 1L) && (.Platform$OS.type != "windows")

# -----------------------------------------------------------------------------
# Function: runner
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - i: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  runner <- function(i) {
    tryCatch(
      fn(scenarios[[i]], i),
      error = function(e) structure(conditionMessage(e), class = "try-error")
    )
  }

  if (use_mc) {
    message("[", label, "] scenario-parallel enabled: workers=", workers, ", scenarios=", n)
    out <- parallel::mclapply(
      X = seq_len(n),
      FUN = runner,
      mc.cores = workers,
      mc.preschedule = FALSE
    )
  } else {
    if (workers > 1L && .Platform$OS.type == "windows") {
      message("[", label, "] windows detected; using serial prediction.")
    }
    out <- lapply(seq_len(n), runner)
  }

  bad <- which(vapply(out, inherits, logical(1), what = "try-error"))
  if (length(bad) > 0L) {
    stop(
      "[", label, "] prediction failed in ", length(bad), "/", n,
      " scenarios. First error: ", as.character(out[[bad[[1]]]])
    )
  }
  out
}

# -----------------------------------------------------------------------------
# Function: decode_params
# Purpose: Decode transformed optimization parameters into natural-scale model parameters.
# Parameters:
#   - par_transformed: Parameter vector in transformed optimization scale.
#   - fit_treatment: Logical flag indicating whether treatment-effect parameters are optimized.
#   - fit_tau_O2: Logical flag indicating whether tau_O2 is estimated instead of fixed.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
decode_params <- function(par_transformed, fit_treatment = TRUE, fit_tau_O2 = FALSE, cfg = NULL) {
  glucose_use <- isTRUE(.first_non_null_local(if (!is.null(cfg)) cfg$glucose else NULL, TRUE))
  glucose_dynamic_use <- glucose_use && isTRUE(.first_non_null_local(if (!is.null(cfg)) cfg$glucose_dynamic else NULL, FALSE))
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(if (!is.null(cfg)) cfg$misseg_loss_survival else NULL, "nullisomy"),
    default = "nullisomy"
  )
  harvest_settings <- resolve_harvest_init_settings_local(cfg = cfg)
  names(par_transformed) <- get_param_names(
    fit_treatment = fit_treatment,
    fit_tau_O2 = fit_tau_O2,
    glucose = glucose_use,
    glucose_dynamic = glucose_dynamic_use,
    misseg_loss_survival = loss_mode,
    harvest_init_multiplier = harvest_settings$enabled,
    harvest_ids = harvest_settings$harvest_ids
  )
  p_mis_base_use <- as.numeric(.first_non_null_local(
    if ("log10_p_mis_base" %in% names(par_transformed)) 10^par_transformed["log10_p_mis_base"] else NULL,
    if (!is.null(cfg)) cfg$p_mis_base else NULL,
    if (!is.null(cfg)) cfg$p_mis_base_init else NULL,
    1e-5
  ))
  if (!is.finite(p_mis_base_use) || p_mis_base_use < 0) p_mis_base_use <- 1e-5
  p_mis_base_use <- clip(p_mis_base_use, 0, 1)
  lam_min <- 10^par_transformed["log10_lam_min"]
  lam_max <- lam_min + exp(par_transformed["delta_lam"])
  tau_O2 <- as.numeric(.first_non_null_local(
    if (isTRUE(fit_tau_O2) && "log10_tau_O2" %in% names(par_transformed)) 10^par_transformed["log10_tau_O2"] else NULL,
    if (!is.null(cfg)) cfg$tau_O2 else NULL,
    if (!is.null(cfg)) cfg$tau_O2_init else NULL,
    2.0
  ))
  if (!is.finite(tau_O2) || tau_O2 <= 0) tau_O2 <- 2.0
  n_O <- as.numeric(.first_non_null_local(par_transformed["n_O"], if (!is.null(cfg)) cfg$n_O_init else NULL, 1.0))
  if (!is.finite(n_O) || n_O < 0) n_O <- 1.0
  O2_crit <- as.numeric(.first_non_null_local(
    if ("log10_O2_crit" %in% names(par_transformed)) 10^par_transformed["log10_O2_crit"] else NULL,
    if (!is.null(cfg)) cfg$o2_crit_init else NULL,
    1.0
  ))
  if (!is.finite(O2_crit) || O2_crit < 0) O2_crit <- 1.0
  glucose_settings <- resolve_glucose_settings_local(run_params = list(), cfg = cfg)
  G_S0 <- as.numeric(.first_non_null_local(
    if (glucose_dynamic_use && "log10_G_S0" %in% names(par_transformed)) 10^par_transformed["log10_G_S0"] else NULL,
    glucose_settings$G_S0
  ))
  kappa_G <- as.numeric(.first_non_null_local(
    if (glucose_dynamic_use && "log10_kappa_G" %in% names(par_transformed)) 10^par_transformed["log10_kappa_G"] else NULL,
    glucose_settings$kappa_G
  ))
  eta_G <- as.numeric(.first_non_null_local(
    if (glucose_dynamic_use && "log10_eta_G" %in% names(par_transformed)) 10^par_transformed["log10_eta_G"] else NULL,
    glucose_settings$eta_G
  ))
  G_c <- as.numeric(.first_non_null_local(
    if (glucose_dynamic_use && "log10_G_c" %in% names(par_transformed)) 10^par_transformed["log10_G_c"] else NULL,
    glucose_settings$G_c
  ))
  tau_G <- as.numeric(.first_non_null_local(
    if (glucose_dynamic_use && "log10_tau_G" %in% names(par_transformed)) 10^par_transformed["log10_tau_G"] else NULL,
    glucose_settings$tau_G
  ))
  if (!is.finite(G_S0) || G_S0 <= 0) G_S0 <- glucose_settings$G_S0
  if (!is.finite(kappa_G) || kappa_G <= 0) kappa_G <- glucose_settings$kappa_G
  if (!is.finite(eta_G) || eta_G <= 0) eta_G <- glucose_settings$eta_G
  if (!is.finite(G_c) || G_c <= 0) G_c <- glucose_settings$G_c
  if (!is.finite(tau_G) || tau_G <= 0) tau_G <- glucose_settings$tau_G
  k_o_use <- as.numeric(.first_non_null_local(
    if ("log10_k_o" %in% names(par_transformed)) 10^par_transformed["log10_k_o"] else NULL,
    if (!is.null(cfg)) cfg$k_o_init else NULL,
    50.0
  ))
  if (!is.finite(k_o_use) || k_o_use <= 0) k_o_use <- 50.0
  gamma_loss_use <- as.numeric(.first_non_null_local(
    if ("log10_gamma_loss" %in% names(par_transformed)) 10^par_transformed["log10_gamma_loss"] else NULL,
    if (!is.null(cfg)) cfg$gamma_loss_init else NULL,
    1.0
  ))
  if (!is.finite(gamma_loss_use) || gamma_loss_use <= 0) gamma_loss_use <- 1.0
  buffer_smax_use <- as.numeric(.first_non_null_local(
    if ("buffer_smax" %in% names(par_transformed)) par_transformed["buffer_smax"] else NULL,
    if (!is.null(cfg)) cfg$buffer_smax_init else NULL,
    0.8
  ))
  if (!is.finite(buffer_smax_use)) buffer_smax_use <- 0.8
  buffer_smax_use <- clip(buffer_smax_use, 0, 1)
  buffer_beta_use <- as.numeric(.first_non_null_local(
    if ("log10_buffer_beta" %in% names(par_transformed)) 10^par_transformed["log10_buffer_beta"] else NULL,
    if (!is.null(cfg)) cfg$buffer_beta_init else NULL,
    1.0
  ))
  if (!is.finite(buffer_beta_use) || buffer_beta_use < 0) buffer_beta_use <- 1.0
  buffer_n_exp_use <- as.numeric(.first_non_null_local(
    if ("log10_buffer_n_exp" %in% names(par_transformed)) 10^par_transformed["log10_buffer_n_exp"] else NULL,
    if (!is.null(cfg)) cfg$buffer_n_exp_init else NULL,
    1.0
  ))
  if (!is.finite(buffer_n_exp_use) || buffer_n_exp_use < 0) buffer_n_exp_use <- 1.0
  p_wgd_use <- as.numeric(.first_non_null_local(
    if ("logit_p_wgd" %in% names(par_transformed)) logit_to_prob(par_transformed["logit_p_wgd"]) else NULL,
    if ("log10_p_wgd" %in% names(par_transformed)) 10^par_transformed["log10_p_wgd"] else NULL,
    if (!is.null(cfg)) cfg$p_wgd_init else NULL,
    1e-4
  ))
  if (!is.finite(p_wgd_use) || p_wgd_use < 0) p_wgd_use <- 0.0
  p_wgd_max_use <- as.numeric(.first_non_null_local(
    if ("logit_p_wgd_max" %in% names(par_transformed)) logit_to_prob(par_transformed["logit_p_wgd_max"]) else NULL,
    if ("log10_p_wgd_max" %in% names(par_transformed)) 10^par_transformed["log10_p_wgd_max"] else NULL,
    if (!is.null(cfg)) cfg$p_wgd_max_init else NULL,
    1e-3
  ))
  if (!is.finite(p_wgd_max_use) || p_wgd_max_use < 0) p_wgd_max_use <- 0.0
  O2_wgd_use <- as.numeric(.first_non_null_local(
    if ("log10_O2_wgd" %in% names(par_transformed)) 10^par_transformed["log10_O2_wgd"] else NULL,
    if (!is.null(cfg)) cfg$O2_wgd_init else NULL,
    0.1
  ))
  if (!is.finite(O2_wgd_use) || O2_wgd_use <= 0) O2_wgd_use <- 1e-12
  out <- list(
    lam_min = lam_min,
    lam_max = lam_max,
    k_o = k_o_use,
    o2_min = as.numeric(.first_non_null_local(
      if (!is.null(cfg)) cfg$o2_min else NULL,
      0.5
    )),
    p_mis_base = p_mis_base_use,
    p_misseg = as.numeric(.first_non_null_local(
      if ("logit_p_misseg" %in% names(par_transformed)) logit_to_prob(par_transformed["logit_p_misseg"]) else NULL,
      if ("log10_p_misseg" %in% names(par_transformed)) 10^par_transformed["log10_p_misseg"] else NULL,
      1e-4
    )),
    k_o_mis = 10^par_transformed["log10_k_o_mis"],
    gamma_loss = gamma_loss_use,
    misseg_loss_survival = loss_mode,
    buffer_smax = buffer_smax_use,
    buffer_beta = buffer_beta_use,
    buffer_n_exp = buffer_n_exp_use,
    p_wgd = p_wgd_use,
    p_wgd_max = p_wgd_max_use,
    O2_wgd = O2_wgd_use,
    o2_S0 = 10^par_transformed["log10_o2_S0"],
    kappa_O = 10^par_transformed["log10_kappa_O"],
    eta_o2 = 10^par_transformed["log10_eta_o2"],
    G_S0 = G_S0,
    kappa_G = kappa_G,
    eta_G = eta_G,
    G_c = G_c,
    tau_G = tau_G,
    glucose_ref_mM = glucose_settings$glucose_ref_mM,
    rho_2N = 10^par_transformed["log10_rho_2N"],
    alpha_o2 = 10^par_transformed["log10_alpha_o2"],
    gamma_growth = par_transformed["gamma_growth"],
    mu_hp = 10^par_transformed["log10_mu_hp"],
    gamma_mu = par_transformed["gamma_mu"],
    O2_crit = O2_crit,
    n_O = n_O,
    k_clear = 10^par_transformed["log10_k_clear"],
    sigma_burden = 10^par_transformed["log10_sigma_burden"],
    tau_O2 = tau_O2,
    c_vol_2N_eff_mm3 = 10^-par_transformed["log10_rho_2N"],
    ratio_4N_2N = 1.0,
    glucose = glucose_use,
    harvest_init_multiplier = harvest_settings$enabled,
    glucose_dynamic = glucose_dynamic_use,
    glucose_stress_mode = resolve_glucose_runtime_mode(
      glucose_dynamic = glucose_dynamic_use,
      glucose_stress_mode = if (isTRUE(glucose_use)) {
        .first_non_null_local(if (!is.null(cfg)) cfg$glucose_stress_mode else NULL, "coupled_to_O2")
      } else {
        "off"
      },
      default_dynamic = glucose_dynamic_use,
      default_static_mode = if (isTRUE(glucose_use)) "coupled_to_O2" else "off"
    ),
    alpha = if (isTRUE(fit_treatment)) 10^par_transformed["log10_alpha"] else 0,
    gamma = if (isTRUE(fit_treatment)) par_transformed["gamma"] else 1
  )
  if (isTRUE(harvest_settings$enabled) && length(harvest_settings$harvest_ids) > 0L) {
    for (harvest_id in harvest_settings$harvest_ids) {
      log_name <- harvest_settings$log_param_names[[harvest_id]]
      natural_name <- harvest_settings$natural_param_names[[harvest_id]]
      log_val <- as.numeric(.first_non_null_local(par_transformed[[log_name]], harvest_settings$prior_center))
      if (!is.finite(log_val)) log_val <- harvest_settings$prior_center
      out[[natural_name]] <- exp(log_val)
    }
  }
  out
}

# -----------------------------------------------------------------------------
# Function: encode_params
# Purpose: Encode natural-scale parameters into transformed optimization scale.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - fit_treatment: Logical flag indicating whether treatment-effect parameters are optimized.
#   - fit_tau_O2: Logical flag indicating whether tau_O2 is estimated instead of fixed.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
encode_params <- function(run_params, fit_treatment = TRUE, fit_tau_O2 = FALSE, cfg = NULL) {
  rp <- as.list(run_params)
  glucose_use <- isTRUE(.first_non_null_local(
    rp$glucose,
    if (!is.null(cfg)) cfg$glucose else NULL,
    TRUE
  ))
  glucose_dynamic_use <- isTRUE(.first_non_null_local(
    rp$glucose_dynamic,
    if (!is.null(cfg)) cfg$glucose_dynamic else NULL,
    FALSE
  )) && glucose_use
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(rp$misseg_loss_survival, if (!is.null(cfg)) cfg$misseg_loss_survival else NULL, "nullisomy"),
    default = "nullisomy"
  )
  harvest_settings <- resolve_harvest_init_settings_local(run_params = rp, cfg = cfg)
# -----------------------------------------------------------------------------
# Function: getv
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - keys: Function-specific input argument.
#   - default: Fallback value used when the input is NULL or invalid.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  getv <- function(keys, default = NA_real_) {
    for (k in keys) {
      v <- rp[[k]]
      if (!is.null(v)) {
        vv <- suppressWarnings(as.numeric(v))
        if (is.finite(vv)) return(vv)
      }
    }
    default
  }
# -----------------------------------------------------------------------------
# Function: need_pos
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
#   - nm: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  need_pos <- function(x, nm) {
    if (!is.finite(x) || x <= 0) stop("Warm-start parameter must be > 0: ", nm)
    x
  }

  lam_min_v <- need_pos(getv(c("lam_min")), "lam_min")
  lam_max_v <- need_pos(getv(c("lam_max"), default = lam_min_v), "lam_max")
  lam_gap_v <- lam_max_v - lam_min_v
  if (!is.finite(lam_gap_v) || lam_gap_v <= 0) {
    lam_gap_v <- 1e-8
    lam_max_v <- lam_min_v + lam_gap_v
  }
  k_o_v <- need_pos(
    getv(c("k_o"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$k_o_init else NULL, 50))),
    "k_o"
  )
  p_misseg_v <- need_pos(getv(c("p_misseg"), default = 1e-4), "p_misseg")
  p_mis_base_v <- need_pos(
    getv(c("p_mis_base"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$p_mis_base else NULL, 1e-5))),
    "p_mis_base"
  )
  k_o_mis_v <- need_pos(getv(c("k_o_mis"), default = 50), "k_o_mis")
  gamma_loss_v <- need_pos(getv(c("gamma_loss"), default = 0.1), "gamma_loss")
  buffer_smax_v <- getv(c("buffer_smax"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$buffer_smax_init else NULL, 0.8)))
  if (!is.finite(buffer_smax_v) || buffer_smax_v < 0 || buffer_smax_v > 1) {
    stop("Warm-start parameter must be in [0,1]: buffer_smax")
  }
  buffer_beta_v <- getv(c("buffer_beta"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$buffer_beta_init else NULL, 1.0)))
  if (!is.finite(buffer_beta_v) || buffer_beta_v <= 0) {
    stop("Warm-start parameter must be > 0: buffer_beta")
  }
  buffer_n_exp_v <- getv(c("buffer_n_exp"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$buffer_n_exp_init else NULL, 1.0)))
  if (!is.finite(buffer_n_exp_v) || buffer_n_exp_v <= 0) {
    stop("Warm-start parameter must be > 0: buffer_n_exp")
  }
  p_wgd_v <- need_pos(
    getv(c("p_wgd"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$p_wgd_init else NULL, 1e-4))),
    "p_wgd"
  )
  o2_s0_upper_v <- as.numeric(.first_non_null_local(
    if (!is.null(cfg)) cfg$o2_S0_upper_bound else NULL,
    read_o2_S0_natural_upper_bound_common(if (!is.null(cfg)) cfg$parameter_table else NULL, fallback = 5.0),
    5.0
  ))
  if (!is.finite(o2_s0_upper_v) || o2_s0_upper_v <= 0) o2_s0_upper_v <- 5.0
  o2_init_v <- need_pos(
    getv(c("o2_S0"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_S0_init else NULL, 0.5))),
    "o2_S0"
  )
  o2_init_v <- min(max(o2_init_v, 1e-6), pmax(1e-6, o2_s0_upper_v))
  kappa_O_v <- need_pos(
    getv(c("kappa_O"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$kappa_O_init else NULL, 1.0))),
    "kappa_O"
  )
  eta_o2_v <- need_pos(
    getv(c("eta_o2"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$eta_o2_init else NULL, 1.0))),
    "eta_o2"
  )
  G_S0_v <- need_pos(
    getv(c("G_S0"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$G_S0_init else NULL, default_glucose_pct_scale()$G_S0))),
    "G_S0"
  )
  kappa_G_v <- need_pos(
    getv(c("kappa_G"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$kappa_G_init else NULL, default_glucose_pct_scale()$kappa_G))),
    "kappa_G"
  )
  eta_G_v <- need_pos(
    getv(c("eta_G"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$eta_G_init else NULL, default_glucose_pct_scale()$eta_G))),
    "eta_G"
  )
  G_c_v <- need_pos(
    getv(c("G_c"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$G_c_init else NULL, default_glucose_pct_scale()$G_c))),
    "G_c"
  )
  tau_G_v <- need_pos(
    getv(c("tau_G"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$tau_G_init else NULL, if (!is.null(cfg)) cfg$tau_G else NULL, if (!is.null(cfg)) cfg$tau_O2_init else NULL, 0.1))),
    "tau_G"
  )
  rho_2N_v <- getv(c("rho_2N"), default = NA_real_)
  rho_2N_v <- need_pos(
    if (is.finite(rho_2N_v) && rho_2N_v > 0) rho_2N_v else default_rho_2N_prior_center(cfg),
    "rho_2N"
  )
  alpha_o2_v <- need_pos(getv(c("alpha_o2"), default = 0.5), "alpha_o2")
  gamma_growth_v <- getv(c("gamma_growth"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$gamma_growth_init else NULL, 2.0)))
  if (!is.finite(gamma_growth_v) || gamma_growth_v <= 0) {
    stop("Warm-start parameter must be > 0: gamma_growth")
  }
  gamma_mu_v <- getv(c("gamma_mu"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$gamma_mu_init else NULL, 1.0)))
  if (!is.finite(gamma_mu_v) || gamma_mu_v <= 0) {
    stop("Warm-start parameter must be > 0: gamma_mu")
  }
  O2_crit_v <- getv(c("O2_crit"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_crit_init else NULL, 1.0)))
  if (!is.finite(O2_crit_v) || O2_crit_v < 0) {
    stop("Warm-start parameter must be >= 0: O2_crit")
  }
  if (O2_crit_v <= 0) O2_crit_v <- 1e-6
  n_O_v <- getv(c("n_O"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$n_O_init else NULL, 1.0)))
  if (!is.finite(n_O_v) || n_O_v < 0) {
    stop("Warm-start parameter must be >= 0: n_O")
  }
  mu_hp_v <- need_pos(
    getv(c("mu_hp"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$mu_hp_init else NULL, 1e-3))),
    "mu_hp"
  )
  k_clear_v <- need_pos(
    getv(c("k_clear"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$k_clear_init else NULL, 1e-3))),
    "k_clear"
  )
  sigma_burden_v <- need_pos(
    getv(c("sigma_burden"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$sigma_burden else NULL, 0.35))),
    "sigma_burden"
  )
  tau_O2_v <- need_pos(
    getv(c("tau_O2"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$tau_O2_init else NULL, 2.0))),
    "tau_O2"
  )

  out <- c(
    log10_lam_min = log10(lam_min_v),
    delta_lam = log(lam_gap_v)
  )
  if (!glucose_use || !isTRUE(glucose_dynamic_use)) {
    out <- c(out, log10_k_o = log10(k_o_v))
  }
  out <- c(
    out,
    log10_p_mis_base = log10(p_mis_base_v),
    logit_p_misseg = prob_to_logit(p_misseg_v, "p_misseg", "warm_start"),
    log10_k_o_mis = log10(k_o_mis_v)
  )
  if (identical(loss_mode, "nullisomy")) {
    out <- c(out, log10_gamma_loss = log10(gamma_loss_v))
  } else {
    out <- c(
      out,
      buffer_smax = buffer_smax_v,
      log10_buffer_beta = log10(buffer_beta_v),
      log10_buffer_n_exp = log10(buffer_n_exp_v)
    )
  }
  out <- c(out, logit_p_wgd = prob_to_logit(p_wgd_v, "p_wgd", "warm_start"))
  if (glucose_use && isTRUE(glucose_dynamic_use)) {
    out <- c(
      out,
      log10_G_S0 = log10(G_S0_v),
      log10_kappa_G = log10(kappa_G_v),
      log10_eta_G = log10(eta_G_v),
      log10_G_c = log10(G_c_v),
      log10_tau_G = log10(tau_G_v)
    )
  }
  out <- c(
    out,
    log10_o2_S0 = log10(o2_init_v),
    log10_kappa_O = log10(kappa_O_v),
    log10_eta_o2 = log10(eta_o2_v),
    log10_rho_2N = log10(rho_2N_v),
    log10_alpha_o2 = log10(alpha_o2_v),
    gamma_growth = gamma_growth_v,
    log10_mu_hp = log10(mu_hp_v),
    gamma_mu = gamma_mu_v,
    log10_O2_crit = log10(O2_crit_v),
    n_O = n_O_v,
    log10_k_clear = log10(k_clear_v),
    log10_sigma_burden = log10(sigma_burden_v)
  )
  if (isTRUE(fit_tau_O2)) {
    out <- c(out, log10_tau_O2 = log10(tau_O2_v))
  }

  if (isTRUE(fit_treatment)) {
    alphav <- need_pos(getv(c("alpha")), "alpha")
    gammav <- getv(c("gamma"))
    out <- c(out, log10_alpha = log10(alphav), gamma = gammav)
  }
  if (isTRUE(harvest_settings$enabled) && length(harvest_settings$harvest_ids) > 0L) {
    for (harvest_id in harvest_settings$harvest_ids) {
      natural_name <- harvest_settings$natural_param_names[[harvest_id]]
      log_name <- harvest_settings$log_param_names[[harvest_id]]
      mult_v <- suppressWarnings(as.numeric(rp[[natural_name]]))
      if (!is.finite(mult_v) || mult_v <= 0) {
        mult_v <- exp(as.numeric(.first_non_null_local(rp[[log_name]], harvest_settings$prior_center)))
      }
      mult_v <- need_pos(mult_v, natural_name)
      out[[log_name]] <- log(mult_v)
    }
  }
  out
}

# -----------------------------------------------------------------------------
# Function: read_init_params_t
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - init_path: Function-specific input argument.
#   - bounds: List containing optimization lower/upper bounds.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
read_init_params_t <- function(init_path, bounds, cfg) {
  if (!file.exists(init_path)) stop("init_params_tsv not found: ", init_path)
  tab <- read.delim(init_path, check.names = FALSE, stringsAsFactors = FALSE)
  full_names <- names(bounds$lower)

  out <- NULL
  if (all(c("transformed_parameter", "transformed_value") %in% names(tab))) {
    vals <- setNames(as.numeric(tab$transformed_value), as.character(tab$transformed_parameter))
    missing_names <- setdiff(full_names, names(vals))
    if ("log10_alpha_o2" %in% missing_names) {
      vals[["log10_alpha_o2"]] <- log10(as.numeric(.first_non_null_local(cfg$alpha_o2_init, 0.5)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("gamma_growth" %in% missing_names) {
      vals[["gamma_growth"]] <- as.numeric(.first_non_null_local(cfg$gamma_growth_init, 2.0))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_p_mis_base" %in% missing_names) {
      vals[["log10_p_mis_base"]] <- log10(as.numeric(.first_non_null_local(cfg$p_mis_base, cfg$p_mis_base_init, 1e-5)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("logit_p_wgd_max" %in% missing_names) {
      if ("log10_p_wgd_max" %in% names(vals) && is.finite(vals[["log10_p_wgd_max"]])) {
        vals[["logit_p_wgd_max"]] <- prob_to_logit(10^vals[["log10_p_wgd_max"]], "p_wgd_max", "legacy_init_params_tsv")
      } else {
        vals[["logit_p_wgd_max"]] <- prob_to_logit(as.numeric(.first_non_null_local(cfg$p_wgd_max_init, 1e-3)), "p_wgd_max", "init")
      }
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("logit_p_wgd" %in% missing_names) {
      if ("log10_p_wgd" %in% names(vals) && is.finite(vals[["log10_p_wgd"]])) {
        vals[["logit_p_wgd"]] <- prob_to_logit(10^vals[["log10_p_wgd"]], "p_wgd", "legacy_init_params_tsv")
      } else {
        vals[["logit_p_wgd"]] <- prob_to_logit(as.numeric(.first_non_null_local(cfg$p_wgd_init, 1e-4)), "p_wgd", "init")
      }
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("logit_p_misseg" %in% missing_names) {
      if ("log10_p_misseg" %in% names(vals) && is.finite(vals[["log10_p_misseg"]])) {
        vals[["logit_p_misseg"]] <- prob_to_logit(10^vals[["log10_p_misseg"]], "p_misseg", "legacy_init_params_tsv")
      } else {
        vals[["logit_p_misseg"]] <- prob_to_logit(1e-4, "p_misseg", "init")
      }
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_O2_wgd" %in% missing_names) {
      vals[["log10_O2_wgd"]] <- log10(as.numeric(.first_non_null_local(cfg$O2_wgd_init, 0.1)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_G_S0" %in% missing_names) {
      vals[["log10_G_S0"]] <- log10(as.numeric(.first_non_null_local(cfg$G_S0_init, default_glucose_pct_scale()$G_S0)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_kappa_G" %in% missing_names) {
      vals[["log10_kappa_G"]] <- log10(as.numeric(.first_non_null_local(cfg$kappa_G_init, default_glucose_pct_scale()$kappa_G)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_eta_G" %in% missing_names) {
      vals[["log10_eta_G"]] <- log10(as.numeric(.first_non_null_local(cfg$eta_G_init, default_glucose_pct_scale()$eta_G)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_G_c" %in% missing_names) {
      vals[["log10_G_c"]] <- log10(as.numeric(.first_non_null_local(cfg$G_c_init, default_glucose_pct_scale()$G_c)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_tau_G" %in% missing_names) {
      vals[["log10_tau_G"]] <- log10(as.numeric(.first_non_null_local(cfg$tau_G_init, cfg$tau_O2_init, 0.1)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_mu_hp" %in% missing_names) {
      vals[["log10_mu_hp"]] <- log10(as.numeric(.first_non_null_local(cfg$mu_hp_init, 1e-3)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("gamma_mu" %in% missing_names) {
      vals[["gamma_mu"]] <- as.numeric(.first_non_null_local(cfg$gamma_mu_init, 1.0))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("n_O" %in% missing_names) {
      vals[["n_O"]] <- as.numeric(.first_non_null_local(cfg$n_O_init, 1.0))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_k_clear" %in% missing_names) {
      vals[["log10_k_clear"]] <- log10(as.numeric(.first_non_null_local(cfg$k_clear_init, 1e-3)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_sigma_burden" %in% missing_names) {
      vals[["log10_sigma_burden"]] <- log10(as.numeric(.first_non_null_local(cfg$sigma_burden, 0.35)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_tau_O2" %in% missing_names) {
      vals[["log10_tau_O2"]] <- log10(as.numeric(.first_non_null_local(cfg$tau_O2_init, 2.0)))
      missing_names <- setdiff(full_names, names(vals))
    }
    init_mult_missing <- grep("^log_init_mult_", missing_names, value = TRUE)
    if (length(init_mult_missing) > 0L) {
      vals[init_mult_missing] <- as.numeric(.first_non_null_local(cfg$prior_center_log_init_mult, 0.0))
      missing_names <- setdiff(full_names, names(vals))
    }
    if (length(missing_names) > 0) {
      stop(
        "init_params_tsv transformed table missing required parameters: ",
        paste(missing_names, collapse = ", ")
      )
    }
    common <- intersect(full_names, names(vals))
    out <- (bounds$lower + bounds$upper) / 2
    names(out) <- full_names
    out[common] <- vals[common]
  }
  if (is.null(out) && all(c("parameter", "value") %in% names(tab))) {
    vals <- setNames(as.numeric(tab$value), as.character(tab$parameter))
    if (all(full_names %in% names(vals))) {
      out <- vals[full_names]
    } else {
      out <- encode_params(vals, fit_treatment = cfg$fit_treatment, fit_tau_O2 = cfg$fit_tau_O2, cfg = cfg)
      out <- out[full_names]
    }
  }
  if (is.null(out) && nrow(tab) >= 1) {
    row1 <- suppressWarnings(as.numeric(tab[1, , drop = TRUE]))
    names(row1) <- names(tab)
    if (all(full_names %in% names(row1))) out <- row1[full_names]
  }
  if (is.null(out)) {
    stop(
      "Could not parse init_params_tsv. Supported formats: ",
      "(parameter,value), (transformed_parameter,transformed_value), or one-row transformed table."
    )
  }
  if (any(!is.finite(out))) stop("init_params_tsv contains non-finite warm-start values.")
  out <- as.numeric(out)
  names(out) <- full_names
  clipped <- clip(out, bounds$lower, bounds$upper)
  if (any(clipped != out)) {
    message("Warm-start values clipped to parameter bounds for: ", paste(full_names[clipped != out], collapse = ", "))
  }
  clipped
}

# -----------------------------------------------------------------------------
# Function: make_bounds
# Purpose: Build transformed parameter bounds for optimization.
# Parameters:
#   - fit_treatment: Logical flag indicating whether treatment-effect parameters are optimized.
#   - fit_tau_O2: Logical flag indicating whether tau_O2 is estimated instead of fixed.
#   - rho_2N_min: Function-specific input argument.
#   - rho_2N_max: Function-specific input argument.
#   - o2_S0_min: Function-specific input argument.
#   - o2_S0_max: Function-specific input argument.
#   - kappa_O_min: Function-specific input argument.
#   - kappa_O_max: Function-specific input argument.
#   - tau_O2_min: Lower bound of tau_O2 when estimated.
#   - tau_O2_max: Upper bound of tau_O2 when estimated.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
make_bounds <- function(fit_treatment = TRUE,
                        fit_tau_O2 = FALSE,
                        glucose = TRUE,
                        glucose_dynamic = FALSE,
                        misseg_loss_survival = "nullisomy",
                        harvest_init_multiplier = FALSE,
                        harvest_ids = NULL,
                        log_init_mult_lower = -1.0,
                        log_init_mult_upper = 1.0,
                        rho_2N_min = 3.2e4, rho_2N_max = 5.6e4,
                        p_mis_base_min = 1e-6, p_mis_base_max = 0.1,
                        o2_S0_min = 1e-3, o2_S0_max = 4.9,
                        kappa_O_min = 1e-3, kappa_O_max = 1e2,
                        eta_o2_min = 1e-3, eta_o2_max = 1e1,
                        G_S0_min = 1.0, G_S0_max = 100.0,
                        kappa_G_min = 1.0, kappa_G_max = 100.0,
                        eta_G_min = 1e-3, eta_G_max = 1e1,
                        G_c_min = 1.0, G_c_max = 100.0,
                        tau_G_min = 1e-3, tau_G_max = 1e3,
                        alpha_o2_min = 1e-2, alpha_o2_max = 10,
                        gamma_growth_min = 2.0, gamma_growth_max = 2.0,
                        gamma_mu_min = 0.3, gamma_mu_max = 3.0,
                        O2_crit_min = 1e-6, O2_crit_max = 2.5,
                        n_O_min = 0.0, n_O_max = 5.0,
                        mu_hp_min = 1e-8, mu_hp_max = 1.0,
                        k_clear_min = 1e-8, k_clear_max = 1.0,
                        sigma_burden_min = 0.05, sigma_burden_max = 1.0,
                        tau_O2_min = 1e-3, tau_O2_max = 1e3) {
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(misseg_loss_survival, "nullisomy"),
    default = "nullisomy"
  )
  rho_2N_min <- as.numeric(rho_2N_min)
  rho_2N_max <- as.numeric(rho_2N_max)
  if (!is.finite(rho_2N_min) || rho_2N_min <= 0) rho_2N_min <- 3.2e4
  if (!is.finite(rho_2N_max) || rho_2N_max <= 0) rho_2N_max <- 5.6e4
  if (rho_2N_min > rho_2N_max) {
    tmp <- rho_2N_min
    rho_2N_min <- rho_2N_max
    rho_2N_max <- tmp
  }
  p_mis_base_min <- as.numeric(p_mis_base_min)
  p_mis_base_max <- as.numeric(p_mis_base_max)
  if (!is.finite(p_mis_base_min) || p_mis_base_min <= 0) p_mis_base_min <- 1e-6
  if (!is.finite(p_mis_base_max) || p_mis_base_max <= 0) p_mis_base_max <- 0.1
  if (p_mis_base_min > p_mis_base_max) {
    tmp <- p_mis_base_min
    p_mis_base_min <- p_mis_base_max
    p_mis_base_max <- tmp
  }
  o2_S0_min <- as.numeric(o2_S0_min)
  o2_S0_max <- as.numeric(o2_S0_max)
  if (!is.finite(o2_S0_min) || o2_S0_min <= 0) o2_S0_min <- 1e-3
  if (!is.finite(o2_S0_max) || o2_S0_max <= 0) o2_S0_max <- 4.9
  if (o2_S0_min > o2_S0_max) {
    tmp <- o2_S0_min
    o2_S0_min <- o2_S0_max
    o2_S0_max <- tmp
  }
  kappa_O_min <- as.numeric(kappa_O_min)
  kappa_O_max <- as.numeric(kappa_O_max)
  if (!is.finite(kappa_O_min) || kappa_O_min <= 0) kappa_O_min <- 1e-3
  if (!is.finite(kappa_O_max) || kappa_O_max <= 0) kappa_O_max <- 1e2
  if (kappa_O_min > kappa_O_max) {
    tmp <- kappa_O_min
    kappa_O_min <- kappa_O_max
    kappa_O_max <- tmp
  }
  eta_o2_min <- as.numeric(eta_o2_min)
  eta_o2_max <- as.numeric(eta_o2_max)
  if (!is.finite(eta_o2_min) || eta_o2_min <= 0) eta_o2_min <- 1e-3
  if (!is.finite(eta_o2_max) || eta_o2_max <= 0) eta_o2_max <- 1e1
  if (eta_o2_min > eta_o2_max) {
    tmp <- eta_o2_min
    eta_o2_min <- eta_o2_max
    eta_o2_max <- tmp
  }
  G_S0_min <- as.numeric(G_S0_min)
  G_S0_max <- as.numeric(G_S0_max)
  if (!is.finite(G_S0_min) || G_S0_min <= 0) G_S0_min <- 1.0
  if (!is.finite(G_S0_max) || G_S0_max <= 0) G_S0_max <- 100.0
  if (G_S0_min > G_S0_max) {
    tmp <- G_S0_min
    G_S0_min <- G_S0_max
    G_S0_max <- tmp
  }
  kappa_G_min <- as.numeric(kappa_G_min)
  kappa_G_max <- as.numeric(kappa_G_max)
  if (!is.finite(kappa_G_min) || kappa_G_min <= 0) kappa_G_min <- 1.0
  if (!is.finite(kappa_G_max) || kappa_G_max <= 0) kappa_G_max <- 100.0
  if (kappa_G_min > kappa_G_max) {
    tmp <- kappa_G_min
    kappa_G_min <- kappa_G_max
    kappa_G_max <- tmp
  }
  eta_G_min <- as.numeric(eta_G_min)
  eta_G_max <- as.numeric(eta_G_max)
  if (!is.finite(eta_G_min) || eta_G_min <= 0) eta_G_min <- 1e-3
  if (!is.finite(eta_G_max) || eta_G_max <= 0) eta_G_max <- 1e1
  if (eta_G_min > eta_G_max) {
    tmp <- eta_G_min
    eta_G_min <- eta_G_max
    eta_G_max <- tmp
  }
  G_c_min <- as.numeric(G_c_min)
  G_c_max <- as.numeric(G_c_max)
  if (!is.finite(G_c_min) || G_c_min <= 0) G_c_min <- 1.0
  if (!is.finite(G_c_max) || G_c_max <= 0) G_c_max <- 100.0
  if (G_c_min > G_c_max) {
    tmp <- G_c_min
    G_c_min <- G_c_max
    G_c_max <- tmp
  }
  tau_G_min <- as.numeric(tau_G_min)
  tau_G_max <- as.numeric(tau_G_max)
  if (!is.finite(tau_G_min) || tau_G_min <= 0) tau_G_min <- 1e-3
  if (!is.finite(tau_G_max) || tau_G_max <= 0) tau_G_max <- 1e3
  if (tau_G_min > tau_G_max) {
    tmp <- tau_G_min
    tau_G_min <- tau_G_max
    tau_G_max <- tmp
  }
  alpha_o2_min <- as.numeric(alpha_o2_min)
  alpha_o2_max <- as.numeric(alpha_o2_max)
  if (!is.finite(alpha_o2_min) || alpha_o2_min <= 0) alpha_o2_min <- 1e-2
  if (!is.finite(alpha_o2_max) || alpha_o2_max <= 0) alpha_o2_max <- 10
  if (alpha_o2_min > alpha_o2_max) {
    tmp <- alpha_o2_min
    alpha_o2_min <- alpha_o2_max
    alpha_o2_max <- tmp
  }
  gamma_growth_min <- as.numeric(gamma_growth_min)
  gamma_growth_max <- as.numeric(gamma_growth_max)
  if (!is.finite(gamma_growth_min) || gamma_growth_min <= 0) gamma_growth_min <- 2.0
  if (!is.finite(gamma_growth_max) || gamma_growth_max <= 0) gamma_growth_max <- 2.0
  if (gamma_growth_min > gamma_growth_max) {
    tmp <- gamma_growth_min
    gamma_growth_min <- gamma_growth_max
    gamma_growth_max <- tmp
  }
  gamma_mu_min <- as.numeric(gamma_mu_min)
  gamma_mu_max <- as.numeric(gamma_mu_max)
  if (!is.finite(gamma_mu_min) || gamma_mu_min <= 0) gamma_mu_min <- 0.3
  if (!is.finite(gamma_mu_max) || gamma_mu_max <= 0) gamma_mu_max <- 3.0
  if (gamma_mu_min > gamma_mu_max) {
    tmp <- gamma_mu_min
    gamma_mu_min <- gamma_mu_max
    gamma_mu_max <- tmp
  }
  O2_crit_min <- as.numeric(O2_crit_min)
  O2_crit_max <- as.numeric(O2_crit_max)
  if (!is.finite(O2_crit_min) || O2_crit_min <= 0) O2_crit_min <- 1e-6
  if (!is.finite(O2_crit_max) || O2_crit_max <= 0) O2_crit_max <- 2.5
  if (O2_crit_min > O2_crit_max) {
    tmp <- O2_crit_min
    O2_crit_min <- O2_crit_max
    O2_crit_max <- tmp
  }
  n_O_min <- as.numeric(n_O_min)
  n_O_max <- as.numeric(n_O_max)
  if (!is.finite(n_O_min) || n_O_min < 0) n_O_min <- 0.0
  if (!is.finite(n_O_max) || n_O_max < 0) n_O_max <- 5.0
  if (n_O_min > n_O_max) {
    tmp <- n_O_min
    n_O_min <- n_O_max
    n_O_max <- tmp
  }
  mu_hp_min <- as.numeric(mu_hp_min)
  mu_hp_max <- as.numeric(mu_hp_max)
  if (!is.finite(mu_hp_min) || mu_hp_min <= 0) mu_hp_min <- 1e-8
  if (!is.finite(mu_hp_max) || mu_hp_max <= 0) mu_hp_max <- 1.0
  if (mu_hp_min > mu_hp_max) {
    tmp <- mu_hp_min
    mu_hp_min <- mu_hp_max
    mu_hp_max <- tmp
  }
  k_clear_min <- as.numeric(k_clear_min)
  k_clear_max <- as.numeric(k_clear_max)
  if (!is.finite(k_clear_min) || k_clear_min <= 0) k_clear_min <- 1e-8
  if (!is.finite(k_clear_max) || k_clear_max <= 0) k_clear_max <- 1.0
  if (k_clear_min > k_clear_max) {
    tmp <- k_clear_min
    k_clear_min <- k_clear_max
    k_clear_max <- tmp
  }
  sigma_burden_min <- as.numeric(sigma_burden_min)
  sigma_burden_max <- as.numeric(sigma_burden_max)
  if (!is.finite(sigma_burden_min) || sigma_burden_min <= 0) sigma_burden_min <- 0.05
  if (!is.finite(sigma_burden_max) || sigma_burden_max <= 0) sigma_burden_max <- 1.0
  if (sigma_burden_min > sigma_burden_max) {
    tmp <- sigma_burden_min
    sigma_burden_min <- sigma_burden_max
    sigma_burden_max <- tmp
  }
  tau_O2_min <- as.numeric(tau_O2_min)
  tau_O2_max <- as.numeric(tau_O2_max)
  if (!is.finite(tau_O2_min) || tau_O2_min <= 0) tau_O2_min <- 1e-3
  if (!is.finite(tau_O2_max) || tau_O2_max <= 0) tau_O2_max <- 1e3
  if (tau_O2_min > tau_O2_max) {
    tmp <- tau_O2_min
    tau_O2_min <- tau_O2_max
    tau_O2_max <- tmp
  }
  lower <- c(
    log10_lam_min = log10(1e-3),
    delta_lam = log(1e-8),
    log10_k_o = log10(1e-1),
    log10_p_mis_base = log10(p_mis_base_min),
    logit_p_misseg = prob_to_logit(1e-8, "p_misseg", "lower_bound"),
    log10_k_o_mis = log10(1e-1),
    log10_gamma_loss = log10(5e-3),
    buffer_smax = 0.0,
    log10_buffer_beta = log10(1e-2),
    log10_buffer_n_exp = log10(1e-1),
    logit_p_wgd = prob_to_logit(1e-8, "p_wgd", "lower_bound"),
    logit_p_wgd_max = prob_to_logit(1e-8, "p_wgd_max", "lower_bound"),
    log10_O2_wgd = log10(1e-6),
    log10_o2_S0 = log10(o2_S0_min),
    log10_kappa_O = log10(kappa_O_min),
    log10_eta_o2 = log10(eta_o2_min),
    log10_G_S0 = log10(G_S0_min),
    log10_kappa_G = log10(kappa_G_min),
    log10_eta_G = log10(eta_G_min),
    log10_G_c = log10(G_c_min),
    log10_tau_G = log10(tau_G_min),
    log10_rho_2N = log10(rho_2N_min),
    log10_alpha_o2 = log10(alpha_o2_min),
    gamma_growth = gamma_growth_min,
    log10_mu_hp = log10(mu_hp_min),
    gamma_mu = gamma_mu_min,
    log10_O2_crit = log10(O2_crit_min),
    n_O = n_O_min,
    log10_k_clear = log10(k_clear_min),
    log10_sigma_burden = log10(sigma_burden_min)
  )
  if (isTRUE(glucose) && isTRUE(glucose_dynamic)) {
    lower <- lower[setdiff(names(lower), "log10_k_o")]
  }
  if (identical(loss_mode, "nullisomy")) {
    lower <- lower[setdiff(names(lower), c("buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp"))]
  } else {
    lower <- lower[setdiff(names(lower), "log10_gamma_loss")]
  }
  lower <- lower[setdiff(names(lower), c("logit_p_wgd_max", "log10_O2_wgd"))]
  if (!isTRUE(glucose)) {
    lower <- lower[setdiff(names(lower), c("log10_G_S0", "log10_kappa_G", "log10_eta_G", "log10_G_c", "log10_tau_G"))]
  } else if (!isTRUE(glucose_dynamic)) {
    lower <- lower[setdiff(
      names(lower),
      c("log10_G_S0", "log10_kappa_G", "log10_eta_G", "log10_G_c", "log10_tau_G")
    )]
  }
  upper <- c(
    log10_lam_min = log10(5),
    delta_lam = log(5),
    log10_k_o = log10(1e4),
    log10_p_mis_base = log10(p_mis_base_max),
    logit_p_misseg = prob_to_logit(0.08, "p_misseg", "upper_bound"),
    log10_k_o_mis = log10(1e4),
    log10_gamma_loss = log10(0.5),
    buffer_smax = 1.0,
    log10_buffer_beta = log10(10.0),
    log10_buffer_n_exp = log10(10.0),
    logit_p_wgd = prob_to_logit(1e-1, "p_wgd", "upper_bound"),
    logit_p_wgd_max = prob_to_logit(1e-1, "p_wgd_max", "upper_bound"),
    log10_O2_wgd = log10(1.0),
    log10_o2_S0 = log10(o2_S0_max),
    log10_kappa_O = log10(kappa_O_max),
    log10_eta_o2 = log10(eta_o2_max),
    log10_G_S0 = log10(G_S0_max),
    log10_kappa_G = log10(kappa_G_max),
    log10_eta_G = log10(eta_G_max),
    log10_G_c = log10(G_c_max),
    log10_tau_G = log10(tau_G_max),
    log10_rho_2N = log10(rho_2N_max),
    log10_alpha_o2 = log10(alpha_o2_max),
    gamma_growth = gamma_growth_max,
    log10_mu_hp = log10(mu_hp_max),
    gamma_mu = gamma_mu_max,
    log10_O2_crit = log10(O2_crit_max),
    n_O = n_O_max,
    log10_k_clear = log10(k_clear_max),
    log10_sigma_burden = log10(sigma_burden_max)
  )
  if (isTRUE(glucose) && isTRUE(glucose_dynamic)) {
    upper <- upper[setdiff(names(upper), "log10_k_o")]
  }
  if (identical(loss_mode, "nullisomy")) {
    upper <- upper[setdiff(names(upper), c("buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp"))]
  } else {
    upper <- upper[setdiff(names(upper), "log10_gamma_loss")]
  }
  upper <- upper[setdiff(names(upper), c("logit_p_wgd_max", "log10_O2_wgd"))]
  if (!isTRUE(glucose)) {
    upper <- upper[setdiff(names(upper), c("log10_G_S0", "log10_kappa_G", "log10_eta_G", "log10_G_c", "log10_tau_G"))]
  } else if (!isTRUE(glucose_dynamic)) {
    upper <- upper[setdiff(
      names(upper),
      c("log10_G_S0", "log10_kappa_G", "log10_eta_G", "log10_G_c", "log10_tau_G")
    )]
  }

  if (isTRUE(fit_tau_O2)) {
    lower <- c(lower, log10_tau_O2 = log10(tau_O2_min))
    upper <- c(upper, log10_tau_O2 = log10(tau_O2_max))
  }

  if (isTRUE(fit_treatment)) {
    lower <- c(
      lower,
      log10_alpha = log10(1e-4),
      gamma = 0.2
    )
    upper <- c(
      upper,
      log10_alpha = log10(5),
      gamma = 4.0
    )
  }

  harvest_ids_use <- unique(as.character(.first_non_null_local(harvest_ids, character(0))))
  harvest_ids_use <- harvest_ids_use[nzchar(harvest_ids_use)]
  if (isTRUE(harvest_init_multiplier) && length(harvest_ids_use) > 0L) {
    lower <- c(lower, setNames(rep(as.numeric(log_init_mult_lower), length(harvest_ids_use)), vapply(harvest_ids_use, harvest_init_log_param_name, character(1))))
    upper <- c(upper, setNames(rep(as.numeric(log_init_mult_upper), length(harvest_ids_use)), vapply(harvest_ids_use, harvest_init_log_param_name, character(1))))
  }

  list(lower = lower, upper = upper)
}

# -----------------------------------------------------------------------------
# Function: parameter_table_specs
# Purpose: Declare the natural-scale input schema and transformed output schema.
# -----------------------------------------------------------------------------
parameter_table_specs <- function() {
  data.frame(
    param_symbol = c(
      "lam_min",
      "lam_max",
      "k_o",
      "p_mis_base",
      "p_misseg",
      "k_o_mis",
      "gamma_loss",
      "buffer_smax",
      "buffer_beta",
      "buffer_n_exp",
      "p_wgd",
      "p_wgd_max",
      "O2_wgd",
      "o2_S0",
      "kappa_O",
      "eta_o2",
      "G_S0",
      "kappa_G",
      "eta_G",
      "G_c",
      "tau_G",
      "rho_2N",
      "alpha_o2",
      "gamma_growth",
      "mu_hp",
      "gamma_mu",
      "O2_crit",
      "n_O",
      "tau_O2",
      "k_clear",
      "sigma_burden",
      "alpha",
      "gamma"
    ),
    param_name = c(
      "log10_lam_min",
      "delta_lam",
      "log10_k_o",
      "log10_p_mis_base",
      "logit_p_misseg",
      "log10_k_o_mis",
      "log10_gamma_loss",
      "buffer_smax",
      "log10_buffer_beta",
      "log10_buffer_n_exp",
      "logit_p_wgd",
      "logit_p_wgd_max",
      "log10_O2_wgd",
      "log10_o2_S0",
      "log10_kappa_O",
      "log10_eta_o2",
      "log10_G_S0",
      "log10_kappa_G",
      "log10_eta_G",
      "log10_G_c",
      "log10_tau_G",
      "log10_rho_2N",
      "log10_alpha_o2",
      "gamma_growth",
      "log10_mu_hp",
      "gamma_mu",
      "log10_O2_crit",
      "n_O",
      "log10_tau_O2",
      "log10_k_clear",
      "log10_sigma_burden",
      "log10_alpha",
      "gamma"
    ),
    transform = c(
      "log10",
      "delta_lam",
      "log10",
      "log10",
      "logit01",
      "log10",
      "log10",
      "identity",
      "log10",
      "log10",
      "logit01",
      "logit01",
      "log10",
      "log10",
      "log10",
      "log10",
      "log10",
      "log10",
      "log10",
      "log10",
      "log10",
      "log10",
      "log10",
      "identity",
      "log10",
      "identity",
      "log10_nonnegative",
      "identity",
      "log10",
      "log10",
      "log10",
      "log10",
      "identity"
    ),
    output_when = c(
      "always",
      "always",
      "glucose_off_or_static",
      "always",
      "always",
      "always",
      "nullisomy_only",
      "buffering_only",
      "buffering_only",
      "buffering_only",
      "always",
      "legacy_inert",
      "legacy_inert",
      "always",
      "always",
      "always",
      "glucose_on_dynamic",
      "glucose_on_dynamic",
      "glucose_on_dynamic",
      "glucose_on_dynamic",
      "glucose_on_dynamic",
      "always",
      "always",
      "always",
      "always",
      "always",
      "always",
      "always",
      "always",
      "always",
      "always",
      "fit_treatment",
      "fit_treatment"
    ),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Function: legacy_wgd_optional_parameter_defaults
# Purpose: Provide inert compatibility rows for obsolete oxygen-triggered WGD fields.
# -----------------------------------------------------------------------------
legacy_wgd_optional_parameter_defaults <- function() {
  data.frame(
    param_symbol = c("p_wgd_max", "O2_wgd"),
    estimate = c(FALSE, FALSE),
    init_value = c(1e-3, 0.2),
    lower_bound = c(1e-6, 0.01),
    upper_bound = c(5e-2, 1.0),
    source = rep("legacy_inert", 2L),
    description = c(
      "legacy oxygen-triggered WGD maximum retained for backward compatibility only",
      "legacy oxygen-triggered WGD threshold retained for backward compatibility only"
    ),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Function: glucose_off_optional_parameter_defaults
# Purpose: Provide compatibility rows for legacy O2-only parameter tables.
# -----------------------------------------------------------------------------
glucose_off_optional_parameter_defaults <- function() {
  glucose_defaults <- default_glucose_pct_scale()
  data.frame(
    param_symbol = c("buffer_smax", "buffer_beta", "buffer_n_exp", "p_wgd_max", "O2_wgd", "G_S0", "kappa_G", "eta_G", "G_c", "tau_G"),
    estimate = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    init_value = c(0.8, 1.0, 1.0, 1e-3, 0.2, glucose_defaults$G_S0, glucose_defaults$kappa_G, glucose_defaults$eta_G, glucose_defaults$G_c, 0.1),
    lower_bound = c(0.0, 0.01, 0.1, 1e-6, 0.01, 70.0, 1.0, 1.0, 10.0, 0.01),
    upper_bound = c(1.0, 10.0, 10.0, 5e-2, 1.0, 100.0, 30.0, 2.0, 60.0, 30.0),
    source = rep("glucose_off_compat", 10L),
    description = c(
      "compatibility default for buffering survival maximum when glucose=FALSE",
      "compatibility default for buffering survival slope when glucose=FALSE",
      "compatibility default for buffering survival exponent when glucose=FALSE",
      "legacy inert WGD maximum retained for backward-compatible parameter tables",
      "legacy inert WGD threshold retained for backward-compatible parameter tables",
      "compatibility default for glucose baseline when glucose=FALSE",
      "compatibility default for glucose drop coefficient when glucose=FALSE",
      "compatibility default for glucose demand exponent when glucose=FALSE",
      "compatibility default for glucose stress threshold when glucose=FALSE",
      "compatibility default for glucose relaxation time constant when glucose=FALSE"
    ),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Function: nullisomy_optional_parameter_defaults
# Purpose: Provide compatibility rows for legacy parameter tables that predate
#   the buffering survival family when buffering is not active.
# -----------------------------------------------------------------------------
nullisomy_optional_parameter_defaults <- function() {
  data.frame(
    param_symbol = c("buffer_smax", "buffer_beta", "buffer_n_exp"),
    estimate = c(FALSE, FALSE, FALSE),
    init_value = c(0.8, 1.0, 1.0),
    lower_bound = c(0.0, 0.01, 0.1),
    upper_bound = c(1.0, 10.0, 10.0),
    source = rep("nullisomy_compat", 3L),
    description = c(
      "compatibility default for buffering survival maximum when misseg_loss_survival=nullisomy",
      "compatibility default for buffering survival slope when misseg_loss_survival=nullisomy",
      "compatibility default for buffering survival exponent when misseg_loss_survival=nullisomy"
    ),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Function: read_parameter_table_natural
# Purpose: Read the natural-scale parameter input table.
# -----------------------------------------------------------------------------
read_parameter_table_natural <- function(path,
                                         glucose = TRUE,
                                         glucose_dynamic = FALSE,
                                         misseg_loss_survival = "nullisomy") {
  if (!file.exists(path)) stop("Parameter table CSV not found: ", path)
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)

  if (!("description" %in% names(tab)) && ("note" %in% names(tab))) {
    tab$description <- tab$note
  }
  if (!("source" %in% names(tab))) {
    tab$source <- "parameter_table"
  }

  req_cols <- c("param_symbol", "estimate", "init_value", "lower_bound", "upper_bound", "description")
  miss <- setdiff(req_cols, names(tab))
  if (length(miss) > 0L) {
    stop("Parameter table missing required columns: ", paste(miss, collapse = ", "))
  }

  tab$param_symbol <- trimws(as.character(tab$param_symbol))
  tab <- tab[nzchar(tab$param_symbol), , drop = FALSE]
  if (anyDuplicated(tab$param_symbol)) {
    dup <- unique(tab$param_symbol[duplicated(tab$param_symbol)])
    stop("Duplicate param_symbol in parameter table: ", paste(dup, collapse = ", "))
  }

  req_symbols <- unique(parameter_table_specs()$param_symbol)
  missing_symbols <- setdiff(req_symbols, tab$param_symbol)
  glucose_use <- isTRUE(glucose)
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(misseg_loss_survival, "nullisomy"),
    default = "nullisomy"
  )
  if (length(missing_symbols) > 0L) {
    compat_rows <- legacy_wgd_optional_parameter_defaults()
    compat_symbols <- intersect(missing_symbols, compat_rows$param_symbol)
    if (length(compat_symbols) > 0L) {
      tab <- bind_rows(
        tab,
        compat_rows[match(compat_symbols, compat_rows$param_symbol), , drop = FALSE]
      )
      missing_symbols <- setdiff(req_symbols, tab$param_symbol)
    }
  }
  if (!glucose_use && length(missing_symbols) > 0L) {
    compat_rows <- glucose_off_optional_parameter_defaults()
    compat_symbols <- intersect(missing_symbols, compat_rows$param_symbol)
    if (length(compat_symbols) > 0L) {
      tab <- bind_rows(
        tab,
        compat_rows[match(compat_symbols, compat_rows$param_symbol), , drop = FALSE]
      )
      missing_symbols <- setdiff(req_symbols, tab$param_symbol)
    }
  }
  if (identical(loss_mode, "nullisomy") && length(missing_symbols) > 0L) {
    compat_rows <- nullisomy_optional_parameter_defaults()
    compat_symbols <- intersect(missing_symbols, compat_rows$param_symbol)
    if (length(compat_symbols) > 0L) {
      tab <- bind_rows(
        tab,
        compat_rows[match(compat_symbols, compat_rows$param_symbol), , drop = FALSE]
      )
      missing_symbols <- setdiff(req_symbols, tab$param_symbol)
    }
  }
  if (length(missing_symbols) > 0L) {
    stop("Parameter table missing required rows: ", paste(missing_symbols, collapse = ", "))
  }

  tab <- tab[match(req_symbols, tab$param_symbol), , drop = FALSE]
  tab$estimate <- vapply(tab$estimate, function(x) as_bool(x, FALSE), logical(1))
  tab$init_value <- suppressWarnings(as.numeric(tab$init_value))
  tab$lower_bound <- suppressWarnings(as.numeric(tab$lower_bound))
  tab$upper_bound <- suppressWarnings(as.numeric(tab$upper_bound))
  tab$source <- as.character(tab$source)
  tab$description <- as.character(tab$description)

  for (col in c("init_value", "lower_bound", "upper_bound")) {
    bad <- !is.finite(tab[[col]])
    if (any(bad)) {
      stop(
        "Non-finite ", col, " in parameter table for: ",
        paste(tab$param_symbol[bad], collapse = ", ")
      )
    }
  }
  if (any(tab$lower_bound > tab$upper_bound)) {
    bad <- tab$param_symbol[tab$lower_bound > tab$upper_bound]
    stop("lower_bound > upper_bound in parameter table for: ", paste(bad, collapse = ", "))
  }
  if (any(tab$init_value < tab$lower_bound | tab$init_value > tab$upper_bound)) {
    bad <- tab$param_symbol[tab$init_value < tab$lower_bound | tab$init_value > tab$upper_bound]
    stop("init_value outside [lower_bound, upper_bound] for: ", paste(bad, collapse = ", "))
  }

  positive_required <- c(
    "lam_min", "lam_max", "k_o", "p_mis_base", "p_misseg", "k_o_mis", "gamma_loss",
    "buffer_beta", "buffer_n_exp", "p_wgd", "p_wgd_max", "O2_wgd", "o2_S0", "kappa_O", "eta_o2", "G_S0", "kappa_G", "eta_G", "G_c", "tau_G", "rho_2N", "alpha_o2",
    "gamma_growth", "mu_hp", "gamma_mu", "tau_O2", "k_clear", "sigma_burden",
    "alpha", "gamma"
  )
  nonnegative_allowed <- c("O2_crit", "n_O")
  unit_interval_required <- c("p_misseg", "p_wgd", "p_wgd_max")
  unit_interval_closed <- c("buffer_smax")

  pos_bad <- tab$param_symbol %in% positive_required &
    (tab$init_value <= 0 | tab$lower_bound <= 0 | tab$upper_bound <= 0)
  if (any(pos_bad)) {
    stop(
      "Positive-only parameters must have init/lower/upper > 0: ",
      paste(tab$param_symbol[pos_bad], collapse = ", ")
    )
  }
  nonneg_bad <- tab$param_symbol %in% nonnegative_allowed &
    (tab$init_value < 0 | tab$lower_bound < 0 | tab$upper_bound < 0)
  if (any(nonneg_bad)) {
    stop(
      "Non-negative parameters must have init/lower/upper >= 0: ",
      paste(tab$param_symbol[nonneg_bad], collapse = ", ")
    )
  }
  unit_bad <- tab$param_symbol %in% unit_interval_required &
    (tab$init_value <= 0 | tab$init_value >= 1 |
       tab$lower_bound <= 0 | tab$lower_bound >= 1 |
       tab$upper_bound <= 0 | tab$upper_bound >= 1)
  if (any(unit_bad)) {
    stop(
      "Probability parameters must have init/lower/upper strictly within (0,1): ",
      paste(tab$param_symbol[unit_bad], collapse = ", ")
    )
  }
  unit_closed_bad <- tab$param_symbol %in% unit_interval_closed &
    (tab$init_value < 0 | tab$init_value > 1 |
       tab$lower_bound < 0 | tab$lower_bound > 1 |
       tab$upper_bound < 0 | tab$upper_bound > 1)
  if (any(unit_closed_bad)) {
    stop(
      "Buffer survival parameters must have init/lower/upper within [0,1]: ",
      paste(tab$param_symbol[unit_closed_bad], collapse = ", ")
    )
  }

  lam_min_row <- tab[tab$param_symbol == "lam_min", , drop = FALSE]
  lam_max_row <- tab[tab$param_symbol == "lam_max", , drop = FALSE]
  if (lam_max_row$init_value[[1]] <= lam_min_row$init_value[[1]]) {
    stop("Parameter table requires lam_max init_value > lam_min init_value.")
  }
  if (lam_max_row$upper_bound[[1]] <= lam_min_row$lower_bound[[1]]) {
    stop("Parameter table requires feasible lam_max > lam_min bounds.")
  }

  tab
}

# -----------------------------------------------------------------------------
# Function: transform_param_slot
# Purpose: Convert one natural-scale slot value to the transformed optimizer scale.
# -----------------------------------------------------------------------------
transform_param_slot <- function(value, transform, param_symbol, slot_label) {
  if (!is.finite(value)) {
    stop("Non-finite ", slot_label, " for parameter ", param_symbol)
  }
  if (identical(transform, "identity")) {
    return(as.numeric(value))
  }
  if (identical(transform, "log10")) {
    if (value <= 0) stop(param_symbol, " ", slot_label, " must be > 0 for log10 transform.")
    return(log10(value))
  }
  if (identical(transform, "log10_nonnegative")) {
    if (value < 0) stop(param_symbol, " ", slot_label, " must be >= 0 for log10 transform.")
    return(log10(max(value, 1e-6)))
  }
  if (identical(transform, "logit01")) {
    return(prob_to_logit(value, param_symbol = param_symbol, slot_label = slot_label))
  }
  stop("Unsupported transform type: ", transform)
}

# -----------------------------------------------------------------------------
# Function: transform_delta_lam_slot
# Purpose: Convert natural lam_min/lam_max slots to the transformed delta_lam slot.
# -----------------------------------------------------------------------------
transform_delta_lam_slot <- function(tab, slot = c("init", "lower", "upper")) {
  slot <- match.arg(slot)
  lam_min <- tab[tab$param_symbol == "lam_min", , drop = FALSE]
  lam_max <- tab[tab$param_symbol == "lam_max", , drop = FALSE]
  if (nrow(lam_min) != 1L || nrow(lam_max) != 1L) {
    stop("lam_min and lam_max must both be present in parameter_table.")
  }

  gap <- switch(
    slot,
    init = lam_max$init_value[[1]] - lam_min$init_value[[1]],
    lower = lam_max$lower_bound[[1]] - lam_min$upper_bound[[1]],
    upper = lam_max$upper_bound[[1]] - lam_min$lower_bound[[1]]
  )
  if (identical(slot, "init") && gap <= 0) {
    stop("lam_max init_value must be > lam_min init_value.")
  }
  gap <- max(as.numeric(gap), 1e-8)
  log(gap)
}

# -----------------------------------------------------------------------------
# Function: build_transformed_parameter_table
# Purpose: Build the transformed optimizer/output table from the natural input table.
# -----------------------------------------------------------------------------
build_transformed_parameter_table <- function(path,
                                              fit_treatment = FALSE,
                                              fit_tau_O2 = FALSE,
                                              O2_growth = TRUE,
                                              glucose = TRUE,
                                              glucose_dynamic = FALSE,
                                              misseg_loss_survival = "nullisomy",
                                              harvest_init_multiplier = FALSE,
                                              harvest_ids = NULL,
                                              prior_center_log_init_mult = 0.0,
                                              log_init_mult_lower = -1.0,
                                              log_init_mult_upper = 1.0) {
  natural_tab <- read_parameter_table_natural(
    path,
    glucose = isTRUE(glucose),
    glucose_dynamic = isTRUE(glucose_dynamic),
    misseg_loss_survival = misseg_loss_survival
  )
  specs <- parameter_table_specs()
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(misseg_loss_survival, "nullisomy"),
    default = "nullisomy"
  )
  include_row <- specs$output_when == "always" |
    (specs$output_when == "fit_treatment" & isTRUE(fit_treatment)) |
    (specs$output_when == "glucose_off" & !isTRUE(glucose)) |
    (specs$output_when == "glucose_on" & isTRUE(glucose)) |
    (specs$output_when == "glucose_off_or_static" & (!isTRUE(glucose) || !isTRUE(glucose_dynamic))) |
    (specs$output_when == "glucose_on_dynamic" & isTRUE(glucose) & isTRUE(glucose_dynamic)) |
    (specs$output_when == "nullisomy_only" & identical(loss_mode, "nullisomy")) |
    (specs$output_when == "buffering_only" & identical(loss_mode, "buffering"))
  specs_out <- specs[include_row, , drop = FALSE]

  natural_row <- function(symbol) {
    natural_tab[natural_tab$param_symbol == symbol, , drop = FALSE]
  }

  out_rows <- lapply(seq_len(nrow(specs_out)), function(i) {
    spec <- specs_out[i, , drop = FALSE]
    row <- natural_row(spec$param_symbol[[1]])
    estimate_effective <- isTRUE(row$estimate[[1]])
    if (!isTRUE(O2_growth) && spec$param_symbol[[1]] %in% c("alpha_o2", "gamma_growth")) {
      estimate_effective <- FALSE
    }
    if (!isTRUE(glucose_dynamic) && spec$param_symbol[[1]] %in% c("G_S0", "kappa_G", "eta_G", "G_c", "tau_G")) {
      estimate_effective <- FALSE
    }
    if (!isTRUE(glucose) && spec$param_symbol[[1]] %in% c("p_wgd_max", "O2_wgd", "G_S0", "kappa_G", "eta_G", "G_c", "tau_G")) {
      estimate_effective <- FALSE
    }
    if (identical(loss_mode, "nullisomy") && spec$param_symbol[[1]] %in% c("buffer_smax", "buffer_beta", "buffer_n_exp")) {
      estimate_effective <- FALSE
    }
    if (identical(loss_mode, "buffering") && spec$param_symbol[[1]] %in% c("gamma_loss")) {
      estimate_effective <- FALSE
    }
    init_t <- if (spec$transform[[1]] == "delta_lam") {
      transform_delta_lam_slot(natural_tab, "init")
    } else {
      transform_param_slot(row$init_value[[1]], spec$transform[[1]], spec$param_symbol[[1]], "init_value")
    }
    lower_t <- if (spec$transform[[1]] == "delta_lam") {
      transform_delta_lam_slot(natural_tab, "lower")
    } else {
      transform_param_slot(row$lower_bound[[1]], spec$transform[[1]], spec$param_symbol[[1]], "lower_bound")
    }
    upper_t <- if (spec$transform[[1]] == "delta_lam") {
      transform_delta_lam_slot(natural_tab, "upper")
    } else {
      transform_param_slot(row$upper_bound[[1]], spec$transform[[1]], spec$param_symbol[[1]], "upper_bound")
    }
    if (lower_t > upper_t) {
      stop("Transformed lower_bound > upper_bound for parameter ", spec$param_name[[1]])
    }

    data.frame(
      param_name = spec$param_name[[1]],
      estimate = estimate_effective,
      init_value = clip(init_t, lower_t, upper_t),
      lower_bound = lower_t,
      upper_bound = upper_t,
      param_prototype = spec$param_symbol[[1]],
      prototype_init_value = row$init_value[[1]],
      prototype_lower_bound = row$lower_bound[[1]],
      prototype_upper_bound = row$upper_bound[[1]],
      source = row$source[[1]],
      note = row$description[[1]],
      stringsAsFactors = FALSE
    )
  })
  transformed_tab <- bind_rows(out_rows)

  harvest_ids_use <- unique(as.character(.first_non_null_local(harvest_ids, character(0))))
  harvest_ids_use <- harvest_ids_use[nzchar(harvest_ids_use)]
  if (isTRUE(harvest_init_multiplier) && length(harvest_ids_use) > 0L) {
    init_rows <- lapply(harvest_ids_use, function(harvest_id) {
      param_name <- harvest_init_log_param_name(harvest_id)
      data.frame(
        param_name = param_name,
        estimate = TRUE,
        init_value = as.numeric(prior_center_log_init_mult),
        lower_bound = as.numeric(log_init_mult_lower),
        upper_bound = as.numeric(log_init_mult_upper),
        param_prototype = "harvest_init_multiplier",
        prototype_init_value = as.numeric(prior_center_log_init_mult),
        prototype_lower_bound = as.numeric(log_init_mult_lower),
        prototype_upper_bound = as.numeric(log_init_mult_upper),
        source = "config",
        note = paste0("harvest-specific log initial-size multiplier for ", harvest_id),
        stringsAsFactors = FALSE
      )
    })
    transformed_tab <- bind_rows(transformed_tab, bind_rows(init_rows))
  }

  full_names <- get_param_names(
    fit_treatment = isTRUE(fit_treatment),
    fit_tau_O2 = isTRUE(fit_tau_O2),
    glucose = isTRUE(glucose),
    glucose_dynamic = isTRUE(glucose_dynamic),
    misseg_loss_survival = loss_mode,
    harvest_init_multiplier = isTRUE(harvest_init_multiplier),
    harvest_ids = harvest_ids_use
  )
  missing_names <- setdiff(full_names, transformed_tab$param_name)
  if (length(missing_names) > 0L) {
    stop("Transformed parameter table missing optimizer rows: ", paste(missing_names, collapse = ", "))
  }
  optim_tab <- transformed_tab[match(full_names, transformed_tab$param_name), , drop = FALSE]

  init <- as.numeric(optim_tab$init_value)
  lower <- as.numeric(optim_tab$lower_bound)
  upper <- as.numeric(optim_tab$upper_bound)
  fixed_in_optimizer <- !vapply(optim_tab$estimate, isTRUE, logical(1))
  lower[fixed_in_optimizer] <- init[fixed_in_optimizer]
  upper[fixed_in_optimizer] <- init[fixed_in_optimizer]
  names(init) <- full_names
  names(lower) <- full_names
  names(upper) <- full_names

  list(
    natural = natural_tab,
    transformed_output = transformed_tab,
    optimizer = list(init = init, lower = lower, upper = upper)
  )
}

# -----------------------------------------------------------------------------
# Function: sync_cfg_from_natural_parameter_table
# Purpose: Copy natural-scale init/bounds from the input table onto cfg fields.
# -----------------------------------------------------------------------------
sync_cfg_from_natural_parameter_table <- function(cfg, natural_tab) {
  slot_val <- function(symbol, slot = c("init", "lower", "upper")) {
    slot <- match.arg(slot)
    row <- natural_tab[natural_tab$param_symbol == symbol, , drop = FALSE]
    if (nrow(row) != 1L) stop("Missing parameter_table row for ", symbol)
    switch(
      slot,
      init = as.numeric(row$init_value[[1]]),
      lower = as.numeric(row$lower_bound[[1]]),
      upper = as.numeric(row$upper_bound[[1]])
    )
  }

  cfg$parameter_table_natural <- natural_tab
  cfg$k_o_init <- slot_val("k_o", "init")
  cfg$k_o_min <- slot_val("k_o", "lower")
  cfg$k_o_max <- slot_val("k_o", "upper")
  cfg$p_mis_base <- slot_val("p_mis_base", "init")
  cfg$gamma_loss_init <- slot_val("gamma_loss", "init")
  cfg$gamma_loss_min <- slot_val("gamma_loss", "lower")
  cfg$gamma_loss_max <- slot_val("gamma_loss", "upper")
  cfg$buffer_smax_init <- slot_val("buffer_smax", "init")
  cfg$buffer_smax_min <- slot_val("buffer_smax", "lower")
  cfg$buffer_smax_max <- slot_val("buffer_smax", "upper")
  cfg$buffer_beta_init <- slot_val("buffer_beta", "init")
  cfg$buffer_beta_min <- slot_val("buffer_beta", "lower")
  cfg$buffer_beta_max <- slot_val("buffer_beta", "upper")
  cfg$buffer_n_exp_init <- slot_val("buffer_n_exp", "init")
  cfg$buffer_n_exp_min <- slot_val("buffer_n_exp", "lower")
  cfg$buffer_n_exp_max <- slot_val("buffer_n_exp", "upper")
  cfg$p_wgd_init <- slot_val("p_wgd", "init")
  cfg$p_wgd_min <- slot_val("p_wgd", "lower")
  cfg$p_wgd_max <- slot_val("p_wgd", "upper")
  cfg$p_wgd_max_init <- slot_val("p_wgd_max", "init")
  cfg$p_wgd_max_min <- slot_val("p_wgd_max", "lower")
  cfg$p_wgd_max_max <- slot_val("p_wgd_max", "upper")
  cfg$O2_wgd_init <- slot_val("O2_wgd", "init")
  cfg$O2_wgd_min <- slot_val("O2_wgd", "lower")
  cfg$O2_wgd_max <- slot_val("O2_wgd", "upper")

  cfg$o2_S0_init <- slot_val("o2_S0", "init")
  cfg$o2_S0_min <- slot_val("o2_S0", "lower")
  cfg$o2_S0_max <- slot_val("o2_S0", "upper")
  cfg$o2_S0_upper_bound <- cfg$o2_S0_max

  cfg$kappa_O_init <- slot_val("kappa_O", "init")
  cfg$kappa_O_min <- slot_val("kappa_O", "lower")
  cfg$kappa_O_max <- slot_val("kappa_O", "upper")

  cfg$eta_o2_init <- slot_val("eta_o2", "init")
  cfg$eta_o2_min <- slot_val("eta_o2", "lower")
  cfg$eta_o2_max <- slot_val("eta_o2", "upper")

  cfg$G_S0_init <- slot_val("G_S0", "init")
  cfg$G_S0_min <- slot_val("G_S0", "lower")
  cfg$G_S0_max <- slot_val("G_S0", "upper")

  cfg$kappa_G_init <- slot_val("kappa_G", "init")
  cfg$kappa_G_min <- slot_val("kappa_G", "lower")
  cfg$kappa_G_max <- slot_val("kappa_G", "upper")

  cfg$eta_G_init <- slot_val("eta_G", "init")
  cfg$eta_G_min <- slot_val("eta_G", "lower")
  cfg$eta_G_max <- slot_val("eta_G", "upper")

  cfg$G_c_init <- slot_val("G_c", "init")
  cfg$G_c_min <- slot_val("G_c", "lower")
  cfg$G_c_max <- slot_val("G_c", "upper")

  cfg$tau_G <- slot_val("tau_G", "init")
  cfg$tau_G_init <- slot_val("tau_G", "init")
  cfg$tau_G_min <- slot_val("tau_G", "lower")
  cfg$tau_G_max <- slot_val("tau_G", "upper")

  cfg$rho_2N_init <- slot_val("rho_2N", "init")
  cfg$rho_2N_min <- slot_val("rho_2N", "lower")
  cfg$rho_2N_max <- slot_val("rho_2N", "upper")

  cfg$alpha_o2_init <- slot_val("alpha_o2", "init")
  cfg$alpha_o2_min <- slot_val("alpha_o2", "lower")
  cfg$alpha_o2_max <- slot_val("alpha_o2", "upper")

  cfg$gamma_growth_init <- slot_val("gamma_growth", "init")
  cfg$gamma_growth_min <- slot_val("gamma_growth", "lower")
  cfg$gamma_growth_max <- slot_val("gamma_growth", "upper")

  cfg$mu_hp_init <- slot_val("mu_hp", "init")
  cfg$mu_hp_min <- slot_val("mu_hp", "lower")
  cfg$mu_hp_max <- slot_val("mu_hp", "upper")

  cfg$gamma_mu_init <- slot_val("gamma_mu", "init")
  cfg$gamma_mu_min <- slot_val("gamma_mu", "lower")
  cfg$gamma_mu_max <- slot_val("gamma_mu", "upper")

  cfg$o2_crit_init <- slot_val("O2_crit", "init")
  cfg$o2_crit_min <- slot_val("O2_crit", "lower")
  cfg$o2_crit_max <- slot_val("O2_crit", "upper")

  cfg$n_O_init <- slot_val("n_O", "init")
  cfg$n_O_min <- slot_val("n_O", "lower")
  cfg$n_O_max <- slot_val("n_O", "upper")

  cfg$k_clear_init <- slot_val("k_clear", "init")
  cfg$k_clear_min <- slot_val("k_clear", "lower")
  cfg$k_clear_max <- slot_val("k_clear", "upper")

  cfg$tau_O2 <- slot_val("tau_O2", "init")
  cfg$tau_O2_init <- slot_val("tau_O2", "init")
  cfg$tau_O2_min <- slot_val("tau_O2", "lower")
  cfg$tau_O2_max <- slot_val("tau_O2", "upper")

  cfg$sigma_burden <- slot_val("sigma_burden", "init")
  cfg$sigma_burden_min <- slot_val("sigma_burden", "lower")
  cfg$sigma_burden_max <- slot_val("sigma_burden", "upper")

  cfg
}

# -----------------------------------------------------------------------------
# Function: sanitize_cfg_for_persistence
# Purpose: Drop runtime-only parameter-table expansions before writing cfg to disk.
# -----------------------------------------------------------------------------
sanitize_cfg_for_persistence <- function(cfg) {
  cfg_out <- cfg
  cfg_out$parameter_table_natural <- NULL
  cfg_out
}

# -----------------------------------------------------------------------------
# Function: finalize_prior_defaults
# Purpose: Fill prior centers that used to depend on YAML model-parameter numerics.
# -----------------------------------------------------------------------------
finalize_prior_defaults <- function(cfg) {
  if (!is.finite(cfg$prior_center_log10_k_o)) cfg$prior_center_log10_k_o <- log10(cfg$k_o_init)
  if (!is.finite(cfg$prior_center_log10_kappa_O)) cfg$prior_center_log10_kappa_O <- log10(cfg$kappa_O_init)
  if (!is.finite(cfg$prior_center_log10_o2_S0)) cfg$prior_center_log10_o2_S0 <- log10(cfg$o2_S0_init)
  if (!is.finite(cfg$prior_center_log10_eta_o2)) cfg$prior_center_log10_eta_o2 <- log10(cfg$eta_o2_init)
  if (!is.finite(cfg$prior_center_log10_gamma_loss)) cfg$prior_center_log10_gamma_loss <- log10(cfg$gamma_loss_init)
  if (!is.finite(cfg$prior_center_buffer_smax)) cfg$prior_center_buffer_smax <- cfg$buffer_smax_init
  if (!is.finite(cfg$prior_center_log10_buffer_beta)) cfg$prior_center_log10_buffer_beta <- log10(cfg$buffer_beta_init)
  if (!is.finite(cfg$prior_center_log10_buffer_n_exp)) cfg$prior_center_log10_buffer_n_exp <- log10(cfg$buffer_n_exp_init)
  if (!is.finite(cfg$prior_center_log10_rho_2N)) {
    cfg$prior_center_log10_rho_2N <- log10(sqrt(cfg$rho_2N_min * cfg$rho_2N_max))
  }
  if (!is.finite(cfg$prior_center_log10_mu_hp)) cfg$prior_center_log10_mu_hp <- log10(cfg$mu_hp_init)
  if (!is.finite(cfg$prior_center_gamma_mu)) cfg$prior_center_gamma_mu <- cfg$gamma_mu_init
  if (!is.finite(cfg$prior_center_log10_k_clear)) cfg$prior_center_log10_k_clear <- log10(cfg$k_clear_init)
  cfg
}

# -----------------------------------------------------------------------------
# Function: simulate_one
# Purpose: Run one forward simulation trajectory for a single scenario.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - scenario: Single scenario data object.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - model_core: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
simulate_one <- function(run_params, scenario, cfg, model_core = NULL) {
  if (is.null(model_core)) {
    if (!is.null(cfg$model_core)) {
      model_core <- cfg$model_core
    } else {
      model_core <- build_model_core(cfg = cfg)
    }
  }

  grid_pre <- model_core$grid_pre
  R0 <- model_core$R0
  init_state <- if (scenario$cohort == "2N") model_core$init_state_2N else model_core$init_state_4N
  init_mult <- scenario_init_multiplier_local(run_params = run_params, scenario = scenario, cfg = cfg)
  init_state <- as.numeric(init_state) * as.numeric(init_mult)

  obs_steps <- as.integer(round(scenario$obs_days / cfg$DT))
  sim_end_step <- as.integer(round(scenario$sim_end_day / cfg$DT))
  vol_by_N <- cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg)

  if (!exists("cpp_o2simps_simulate_one", mode = "function", inherits = TRUE)) {
    stop("Required C++ function missing: cpp_o2simps_simulate_one")
  }

  o2_s0_upper_use <- as.numeric(.first_non_null_local(cfg$o2_S0_upper_bound, 5.0))
  if (!is.finite(o2_s0_upper_use) || o2_s0_upper_use <= 0) o2_s0_upper_use <- 5.0
  o2_S0_use <- as.numeric(.first_non_null_local(run_params$o2_S0, cfg$o2_S0_init, 0.5))
  if (!is.finite(o2_S0_use) || o2_S0_use < 0) o2_S0_use <- 0.5
  o2_S0_use <- min(max(o2_S0_use, 0), o2_s0_upper_use)
  kappa_O_use <- as.numeric(.first_non_null_local(run_params$kappa_O, cfg$kappa_O_init, 1.0))
  if (!is.finite(kappa_O_use) || kappa_O_use <= 0) kappa_O_use <- 1.0
  eta_o2_use <- as.numeric(.first_non_null_local(run_params$eta_o2, cfg$eta_o2_init, 1.0))
  if (!is.finite(eta_o2_use) || eta_o2_use <= 0) eta_o2_use <- 1.0
  O2_crit_use <- as.numeric(.first_non_null_local(run_params$O2_crit, cfg$o2_crit_init, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  tau_O2_use <- as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  o2_Nref_use <- as.numeric(.first_non_null_local(cfg$o2_Nref, cfg$init_total_size, 1e6))
  if (!is.finite(o2_Nref_use) || o2_Nref_use <= 0) o2_Nref_use <- 1e6
  o2_min_use <- as.numeric(.first_non_null_local(run_params$o2_min, cfg$o2_min, 0.5))
  if (!is.finite(o2_min_use) || o2_min_use < 0) o2_min_use <- 0.5
  o2_min_use <- min(max(o2_min_use, 0), o2_s0_upper_use)
  p_wgd_use <- as.numeric(.first_non_null_local(run_params$p_wgd, cfg$p_wgd_init, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
  p_wgd_max_use <- as.numeric(.first_non_null_local(run_params$p_wgd_max, cfg$p_wgd_max_init, 0.0))
  if (!is.finite(p_wgd_max_use)) p_wgd_max_use <- 0.0
  O2_wgd_use <- as.numeric(.first_non_null_local(run_params$O2_wgd, cfg$O2_wgd_init, 0.1))
  if (!is.finite(O2_wgd_use) || O2_wgd_use <= 0) O2_wgd_use <- 1e-12
  glucose_settings <- resolve_glucose_settings_local(run_params = run_params, cfg = cfg)
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(run_params$misseg_loss_survival, cfg$misseg_loss_survival, "nullisomy"),
    "nullisomy"
  )
  boundary_mode <- as.character(.first_non_null_local(run_params$boundary, "drop"))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 0)

  sim <- cpp_o2simps_simulate_one(list(
    init_state = as.numeric(init_state),
    N0min = as.integer(cfg$N_MIN),
    N0max = as.integer(cfg$N_MAX),
    N1min = as.integer(cfg$N_MIN),
    N1max = as.integer(cfg$N_MAX),
    obs_steps = as.integer(obs_steps),
    sim_end_step = as.integer(sim_end_step),
    DT = as.numeric(cfg$DT),
    dose = as.numeric(scenario$dose),
    dose_ref = as.numeric(cfg$dose_ref),
    treat_day = as.numeric(scenario$treat_day),
    fit_treatment = isTRUE(cfg$fit_treatment),
    alpha = as.numeric(.first_non_null_local(run_params$alpha, 0.0)),
    gamma = as.numeric(.first_non_null_local(run_params$gamma, 1.0)),
    tx_mult_min = as.numeric(cfg$tx_mult_min),
    crowding_enabled = isTRUE(cfg$Crowding),
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_crit = as.numeric(O2_crit_use),
    o2_feedback = isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE)),
    o2_S0 = as.numeric(o2_S0_use),
    kappa_O = as.numeric(kappa_O_use),
    tau_O2 = as.numeric(tau_O2_use),
    o2_Nref = as.numeric(o2_Nref_use),
    o2_min = as.numeric(o2_min_use),
    eta_o2 = as.numeric(eta_o2_use),
    o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01)),
    o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct, 0.005)),
    o2_cache_profile = isTRUE(.first_non_null_local(cfg$o2_cache_profile, FALSE)),
    G_S0 = as.numeric(glucose_settings$G_S0),
    kappa_G = as.numeric(glucose_settings$kappa_G),
    tau_G = as.numeric(glucose_settings$tau_G),
    G_c = as.numeric(glucose_settings$G_c),
    eta_G = as.numeric(glucose_settings$eta_G),
    glucose = isTRUE(glucose_settings$glucose),
    glucose_dynamic = isTRUE(glucose_settings$glucose_dynamic),
    O2_growth = isTRUE(.first_non_null_local(cfg$O2_growth, TRUE)),
    lam_min = as.numeric(run_params$lam_min),
    lam_max = as.numeric(run_params$lam_max),
    k_o = as.numeric(run_params$k_o),
    has_p_misseg = !is.null(run_params$p_misseg),
    p_mis_base = as.numeric(.first_non_null_local(run_params$p_mis_base, cfg$p_mis_base, cfg$p_mis_base_init, 1e-5)),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = 0.0,
    p_wgd = as.numeric(p_wgd_use),
    p_wgd_max = as.numeric(p_wgd_max_use),
    O2_wgd = as.numeric(O2_wgd_use),
    boundary = boundary_mode,
    eps_tail = as.numeric(1e-8),
    gamma_loss = as.numeric(.first_non_null_local(run_params$gamma_loss, 0.1)),
    misseg_loss_survival = as.character(loss_mode),
    buffer_smax = as.numeric(.first_non_null_local(run_params$buffer_smax, cfg$buffer_smax_init, 0.8)),
    buffer_beta = as.numeric(.first_non_null_local(run_params$buffer_beta, cfg$buffer_beta_init, 1.0)),
    buffer_n_exp = as.numeric(.first_non_null_local(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1.0)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = 0.0,
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
    mu_hp = as.numeric(.first_non_null_local(run_params$mu_hp, cfg$mu_hp_init, 1e-3)),
    gamma_mu = as.numeric(.first_non_null_local(run_params$gamma_mu, cfg$gamma_mu_init, 1.0)),
    n_O = as.numeric(.first_non_null_local(run_params$n_O, cfg$n_O_init, 1.0)),
    # Config mode is authoritative for objective/simulation calls.
    ploidy_O2_death = canonical_ploidy_o2_death_mode(
      .first_non_null_local(cfg$ploidy_O2_death, "diploid_NULL"),
      "diploid_NULL"
    ),
    glucose_stress_mode = as.character(glucose_settings$glucose_stress_mode),
    start_with = assert_canonical_start_with_mode(
      .first_non_null_local(cfg$start_with, "ploidy")
    ),
    k_clear = as.numeric(.first_non_null_local(run_params$k_clear, cfg$k_clear_init, 1e-3)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(burden_floor),
    return_full_trajectory = FALSE
  ))

  frac_N <- as.numeric(.first_non_null_local(sim$frac_N_live, sim$frac_N))
  names(frac_N) <- as.character(grid_pre)

  list(
    Ntot_obs = as.numeric(.first_non_null_local(sim$Ntot_total_obs, sim$Ntot_obs)),
    Vmm3_obs = as.numeric(.first_non_null_local(sim$Vmm3_total_obs, sim$Vmm3_obs)),
    frac_N = frac_N,
    Ntot_live_obs = as.numeric(.first_non_null_local(sim$Ntot_live_obs, sim$Ntot_obs)),
    Ntot_dead_hypoxia_obs = as.numeric(.first_non_null_local(sim$Ntot_dead_hypoxia_obs, rep(0, length(obs_steps)))),
    Ntot_dead_buffer_obs = as.numeric(.first_non_null_local(sim$Ntot_dead_buffer_obs, rep(0, length(obs_steps)))),
    Ntot_dead_total_obs = as.numeric(.first_non_null_local(sim$Ntot_dead_total_obs, rep(0, length(obs_steps)))),
    Ntot_total_obs = as.numeric(.first_non_null_local(sim$Ntot_total_obs, sim$Ntot_obs)),
    Vmm3_live_obs = as.numeric(.first_non_null_local(sim$Vmm3_live_obs, sim$Vmm3_obs)),
    Vmm3_dead_hypoxia_obs = as.numeric(.first_non_null_local(sim$Vmm3_dead_hypoxia_obs, rep(0, length(obs_steps)))),
    Vmm3_dead_buffer_obs = as.numeric(.first_non_null_local(sim$Vmm3_dead_buffer_obs, rep(0, length(obs_steps)))),
    Vmm3_dead_total_obs = as.numeric(.first_non_null_local(sim$Vmm3_dead_total_obs, rep(0, length(obs_steps)))),
    Vmm3_total_obs = as.numeric(.first_non_null_local(sim$Vmm3_total_obs, sim$Vmm3_obs))
  )
}

# -----------------------------------------------------------------------------
# Function: evaluate_objective_components_raw
# Purpose: Compute raw loss components before optional aggregation and scaling.
# Parameters:
#   - par_transformed: Parameter vector in transformed optimization scale.
#   - scenarios: List of scenario-specific observation data and metadata.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
evaluate_objective_components_raw <- function(par_transformed, scenarios, cfg) {
  cfg_eval <- cfg
  rp <- decode_params(
    par_transformed,
    fit_treatment = cfg_eval$fit_treatment,
    fit_tau_O2 = cfg_eval$fit_tau_O2,
    cfg = cfg_eval
  )
  model_core <- if (!is.null(cfg_eval$model_core)) cfg_eval$model_core else build_model_core(cfg = cfg_eval)
  scenario_cpp <- if (!is.null(cfg_eval$scenario_cpp)) cfg_eval$scenario_cpp else prepare_cpp_scenarios(scenarios, cfg_eval)
  init_mult_vec <- vapply(
    scenarios,
    function(sc) scenario_init_multiplier_local(run_params = rp, scenario = sc, cfg = cfg_eval),
    numeric(1)
  )
  vol_by_N <- cell_volume_mm3_by_N(model_core$grid_pre, run_params = rp, cfg = cfg_eval)

  o2_s0_upper_use <- as.numeric(.first_non_null_local(cfg_eval$o2_S0_upper_bound, 5.0))
  if (!is.finite(o2_s0_upper_use) || o2_s0_upper_use <= 0) o2_s0_upper_use <- 5.0
  o2_S0_use <- as.numeric(.first_non_null_local(rp$o2_S0, cfg_eval$o2_S0_init, 0.5))
  if (!is.finite(o2_S0_use) || o2_S0_use < 0) o2_S0_use <- 0.5
  o2_S0_use <- min(max(o2_S0_use, 0), o2_s0_upper_use)
  kappa_O_use <- as.numeric(.first_non_null_local(rp$kappa_O, cfg_eval$kappa_O_init, 1.0))
  if (!is.finite(kappa_O_use) || kappa_O_use <= 0) kappa_O_use <- 1.0
  eta_o2_use <- as.numeric(.first_non_null_local(rp$eta_o2, cfg_eval$eta_o2_init, 1.0))
  if (!is.finite(eta_o2_use) || eta_o2_use <= 0) eta_o2_use <- 1.0
  O2_crit_use <- as.numeric(.first_non_null_local(rp$O2_crit, cfg_eval$o2_crit_init, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  tau_O2_use <- as.numeric(.first_non_null_local(rp$tau_O2, cfg_eval$tau_O2, cfg_eval$tau_O2_init, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  o2_Nref_use <- as.numeric(.first_non_null_local(cfg_eval$o2_Nref, cfg_eval$init_total_size, 1e6))
  if (!is.finite(o2_Nref_use) || o2_Nref_use <= 0) o2_Nref_use <- 1e6
  o2_min_use <- as.numeric(.first_non_null_local(rp$o2_min, cfg_eval$o2_min, 0.5))
  if (!is.finite(o2_min_use) || o2_min_use < 0) o2_min_use <- 0.5
  o2_min_use <- min(max(o2_min_use, 0), o2_s0_upper_use)
  p_wgd_use <- as.numeric(.first_non_null_local(rp$p_wgd, cfg_eval$p_wgd_init, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
  p_wgd_max_use <- as.numeric(.first_non_null_local(rp$p_wgd_max, cfg_eval$p_wgd_max_init, 0.0))
  if (!is.finite(p_wgd_max_use)) p_wgd_max_use <- 0.0
  O2_wgd_use <- as.numeric(.first_non_null_local(rp$O2_wgd, cfg_eval$O2_wgd_init, 0.1))
  if (!is.finite(O2_wgd_use) || O2_wgd_use <= 0) O2_wgd_use <- 1e-12
  glucose_settings <- resolve_glucose_settings_local(run_params = rp, cfg = cfg_eval)
  loss_mode <- canonical_misseg_loss_survival_mode(
    .first_non_null_local(rp$misseg_loss_survival, cfg_eval$misseg_loss_survival, "nullisomy"),
    "nullisomy"
  )
  boundary_mode <- as.character(.first_non_null_local(rp$boundary, "drop"))
  sigma_burden_use <- as.numeric(.first_non_null_local(rp$sigma_burden, cfg_eval$sigma_burden, 0.35))
  if (!is.finite(sigma_burden_use) || sigma_burden_use <= 0) sigma_burden_use <- 0.35
  sigma_ploidy_use <- as.numeric(.first_non_null_local(cfg_eval$sigma_ploidy, 0.08))
  if (!is.finite(sigma_ploidy_use) || sigma_ploidy_use <= 0) sigma_ploidy_use <- 0.08
  mu_hp_use <- as.numeric(.first_non_null_local(rp$mu_hp, cfg_eval$mu_hp_init, 1e-3))
  if (!is.finite(mu_hp_use) || mu_hp_use < 0) mu_hp_use <- 0.0
  alpha_o2_use <- as.numeric(.first_non_null_local(rp$alpha_o2, cfg_eval$alpha_o2_init, 0.5))
  if (!is.finite(alpha_o2_use)) alpha_o2_use <- as.numeric(.first_non_null_local(cfg_eval$alpha_o2_init, 0.5))
  if (!isTRUE(.first_non_null_local(cfg_eval$O2_growth, TRUE))) {
    alpha_o2_use <- -abs(alpha_o2_use)
  }
  gamma_mu_use <- as.numeric(.first_non_null_local(rp$gamma_mu, cfg_eval$gamma_mu_init, 1.0))
  if (!is.finite(gamma_mu_use) || gamma_mu_use <= 0) gamma_mu_use <- 1.0
  k_clear_use <- as.numeric(.first_non_null_local(rp$k_clear, cfg_eval$k_clear_init, 1e-3))
  if (!is.finite(k_clear_use) || k_clear_use < 0) k_clear_use <- 0.0
  start_with_mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg_eval$start_with, "ploidy"))
  mu_by_N <- if (identical(start_with_mode, "chr_number")) {
    as.numeric(model_core$grid_pre)
  } else {
    vapply(
      model_core$grid_pre,
      function(n) weighted_ploidy_from_total_N(n, chr_lengths_bp = cfg_eval$chr_lengths_bp) * cfg_eval$N_UNIT,
      numeric(1)
    )
  }

  comp <- cpp_o2simps_objective_components_map(
    scenario_data = list(
      cohort_code = as.integer(scenario_cpp$cohort_code),
      dose_vec = as.numeric(scenario_cpp$dose),
      treat_day_vec = as.numeric(scenario_cpp$treat_day),
      obs_steps_list = scenario_cpp$obs_steps,
      sim_end_step_vec = as.integer(scenario_cpp$sim_end_step),
      obs_burden_list = scenario_cpp$obs_burden,
      keep_burden_list = scenario_cpp$keep_burden,
      ploidy_z_list = scenario_cpp$ploidy_z,
      init_mult_vec = as.numeric(init_mult_vec)
    ),
    objective_data = list(
      mu_by_N = as.numeric(mu_by_N),
      sigma_burden = as.numeric(sigma_burden_use),
      sigma_ploidy = as.numeric(sigma_ploidy_use)
    ),
    state_data = list(
      init_state_2N = as.numeric(model_core$init_state_2N),
      init_state_4N = as.numeric(model_core$init_state_4N),
      N0min = as.integer(cfg_eval$N_MIN),
      N0max = as.integer(cfg_eval$N_MAX),
      N1min = as.integer(cfg_eval$N_MIN),
      N1max = as.integer(cfg_eval$N_MAX),
      N_unit = as.integer(cfg_eval$N_UNIT),
      vol_by_N = as.numeric(vol_by_N)
    ),
    sim_args = list(
      DT = as.numeric(cfg_eval$DT),
      dose_ref = as.numeric(cfg_eval$dose_ref),
      fit_treatment = isTRUE(cfg_eval$fit_treatment),
      alpha = as.numeric(.first_non_null_local(rp$alpha, 0.0)),
      gamma = as.numeric(.first_non_null_local(rp$gamma, 1.0)),
      tx_mult_min = as.numeric(cfg_eval$tx_mult_min),
      crowding_enabled = isTRUE(cfg_eval$Crowding),
      crowding = as.character(cfg_eval$crowding),
      K = as.numeric(cfg_eval$K),
      min_pop = as.numeric(cfg_eval$min_pop),
      O2_crit = as.numeric(O2_crit_use),
      o2_feedback = isTRUE(.first_non_null_local(cfg_eval$o2_burden_feedback, TRUE)),
      o2_S0 = as.numeric(o2_S0_use),
      kappa_O = as.numeric(kappa_O_use),
      tau_O2 = as.numeric(tau_O2_use),
      o2_Nref = as.numeric(o2_Nref_use),
      o2_min = as.numeric(o2_min_use),
      eta_o2 = as.numeric(eta_o2_use),
      o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg_eval$o2_cache_bin_pct, 0.01)),
      o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg_eval$o2_cache_hysteresis_pct, 0.005)),
      o2_cache_profile = isTRUE(.first_non_null_local(cfg_eval$o2_cache_profile, FALSE)),
      G_S0 = as.numeric(glucose_settings$G_S0),
      kappa_G = as.numeric(glucose_settings$kappa_G),
      tau_G = as.numeric(glucose_settings$tau_G),
      G_c = as.numeric(glucose_settings$G_c),
      eta_G = as.numeric(glucose_settings$eta_G),
      glucose = isTRUE(glucose_settings$glucose),
      glucose_dynamic = isTRUE(glucose_settings$glucose_dynamic),
      lam_min = as.numeric(rp$lam_min),
      lam_max = as.numeric(rp$lam_max),
      k_o = as.numeric(rp$k_o),
      has_p_misseg = !is.null(rp$p_misseg),
      p_mis_base = as.numeric(.first_non_null_local(rp$p_mis_base, cfg_eval$p_mis_base, cfg_eval$p_mis_base_init, 1e-5)),
      p_misseg = as.numeric(.first_non_null_local(rp$p_misseg, 0.0)),
      k_o_mis = as.numeric(.first_non_null_local(rp$k_o_mis, 50.0)),
      has_pmis_endpoints = FALSE,
      pmis_O2_0 = 0.0,
      pmis_O2_1 = 0.0,
      p_const = 0.0,
      p_wgd = as.numeric(p_wgd_use),
      p_wgd_max = as.numeric(p_wgd_max_use),
      O2_wgd = as.numeric(O2_wgd_use),
      boundary = boundary_mode,
      eps_tail = as.numeric(1e-8),
      gamma_loss = as.numeric(.first_non_null_local(rp$gamma_loss, 0.1)),
      misseg_loss_survival = as.character(loss_mode),
      buffer_smax = as.numeric(.first_non_null_local(rp$buffer_smax, cfg_eval$buffer_smax_init, 0.8)),
      buffer_beta = as.numeric(.first_non_null_local(rp$buffer_beta, cfg_eval$buffer_beta_init, 1.0)),
      buffer_n_exp = as.numeric(.first_non_null_local(rp$buffer_n_exp, cfg_eval$buffer_n_exp_init, 1.0)),
      beta_size = 0.0,
      alpha_o2 = as.numeric(alpha_o2_use),
      gamma_growth = as.numeric(.first_non_null_local(rp$gamma_growth, cfg_eval$gamma_growth_init, 2.0)),
      mu_hp = as.numeric(mu_hp_use),
      gamma_mu = as.numeric(gamma_mu_use),
      n_O = as.numeric(.first_non_null_local(rp$n_O, cfg_eval$n_O_init, 1.0)),
      ploidy_O2_death = canonical_ploidy_o2_death_mode(
        .first_non_null_local(cfg_eval$ploidy_O2_death, "diploid_NULL"),
        "diploid_NULL"
      ),
      glucose_stress_mode = as.character(glucose_settings$glucose_stress_mode),
      start_with = start_with_mode,
      k_clear = as.numeric(k_clear_use),
      burden_log_eps = as.numeric(.first_non_null_local(cfg_eval$burden_log_eps, 1e-12))
    )
  )

  L_b <- as.numeric(comp$L_b)
  L_p <- as.numeric(comp$L_p)
  if (!is.finite(L_b)) L_b <- 0
  if (!is.finite(L_p)) L_p <- 0
  cache_g_build <- as.integer(.first_non_null_local(comp$cache_g_build, 0L))
  cache_g_hit <- as.integer(.first_non_null_local(comp$cache_g_hit, 0L))
  cache_g_hysteresis <- as.integer(.first_non_null_local(comp$cache_g_hysteresis, 0L))
  if (!is.finite(cache_g_build)) cache_g_build <- 0L
  if (!is.finite(cache_g_hit)) cache_g_hit <- 0L
  if (!is.finite(cache_g_hysteresis)) cache_g_hysteresis <- 0L
  if (isTRUE(.first_non_null_local(cfg_eval$o2_cache_profile, FALSE)) &&
      isTRUE(.first_non_null_local(cfg_eval$trace_obj, FALSE))) {
    cache_total <- cache_g_build + cache_g_hit
    cache_hit_rate <- if (cache_total > 0L) cache_g_hit / cache_total else NA_real_
    message(
      "[o2_cache] build=", cache_g_build,
      ", hit=", cache_g_hit,
      ", hysteresis=", cache_g_hysteresis,
      ", hit_rate=", if (is.finite(cache_hit_rate)) signif(cache_hit_rate, 4) else "NA"
    )
  }
  list(
    L_b = L_b,
    L_p = L_p,
    burden_nll_total = as.numeric(.first_non_null_local(comp$burden_nll_total, 0)),
    ploidy_nll_total = as.numeric(.first_non_null_local(comp$ploidy_nll_total, 0)),
    objective_burden_neg2loglik_raw = as.numeric(.first_non_null_local(comp$objective_burden_neg2loglik_raw, 0)),
    objective_ploidy_neg2loglik_raw = as.numeric(.first_non_null_local(comp$objective_ploidy_neg2loglik_raw, 0)),
    n_burden = as.integer(.first_non_null_local(comp$n_burden, 0L)),
    n_burden_obs_total = as.integer(.first_non_null_local(comp$n_burden_obs_total, 0L)),
    n_ploidy = as.integer(.first_non_null_local(comp$n_ploidy, 0L)),
    n_ploidy_obs_total = as.integer(.first_non_null_local(comp$n_ploidy_obs_total, 0L)),
    n_ploidy_2N = as.integer(.first_non_null_local(comp$n_ploidy_2N, 0L)),
    n_ploidy_4N = as.integer(.first_non_null_local(comp$n_ploidy_4N, 0L)),
    cache_g_build = cache_g_build,
    cache_g_hit = cache_g_hit,
    cache_g_hysteresis = cache_g_hysteresis
  )
}

# -----------------------------------------------------------------------------
# Function: evaluate_objective_components
# Purpose: Compute full objective decomposition used for reporting and optimization.
# Parameters:
#   - par_transformed: Parameter vector in transformed optimization scale.
#   - scenarios: List of scenario-specific observation data and metadata.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
evaluate_objective_components <- function(par_transformed, scenarios, cfg) {
  raw <- evaluate_objective_components_raw(par_transformed, scenarios = scenarios, cfg = cfg)
  L_b <- raw$L_b
  L_p <- raw$L_p

  L_data <- L_b + L_p
  prior <- compute_soft_prior_penalty(par_transformed, cfg = cfg)
  lambda_prior <- as.numeric(.first_non_null_local(cfg$lambda_prior, 0))
  if (!is.finite(lambda_prior) || lambda_prior < 0) lambda_prior <- 0
  L_prior_raw <- as.numeric(.first_non_null_local(prior$L_prior_raw, 0))
  if (!is.finite(L_prior_raw) || L_prior_raw < 0) L_prior_raw <- 0
  L_prior <- if (isTRUE(.first_non_null_local(cfg$use_soft_prior, FALSE))) lambda_prior * L_prior_raw else 0
  if (!is.finite(L_prior) || L_prior < 0) L_prior <- 0
  L <- L_data + L_prior
  if (!is.finite(L)) L <- 1e9

  if (cfg$trace_obj) {
    cat(sprintf(
      "L=%.6f (data=%.6f, prior=%.6f; burden_norm=%.6f, ploidy_norm=%.6f)\n",
      L, L_data, L_prior, L_b, L_p
    ))
  }
  list(
    L = L,
    L_data = L_data,
    L_prior = L_prior,
    L_prior_raw = L_prior_raw,
    L_b = L_b,
    L_p = L_p,
    burden_nll_total = as.numeric(.first_non_null_local(raw$burden_nll_total, 0)),
    ploidy_nll_total = as.numeric(.first_non_null_local(raw$ploidy_nll_total, 0)),
    objective_burden_neg2loglik_raw = as.numeric(.first_non_null_local(raw$objective_burden_neg2loglik_raw, 0)),
    objective_ploidy_neg2loglik_raw = as.numeric(.first_non_null_local(raw$objective_ploidy_neg2loglik_raw, 0)),
    n_burden = raw$n_burden,
    n_burden_obs_total = raw$n_burden_obs_total,
    n_ploidy = raw$n_ploidy,
    n_ploidy_obs_total = raw$n_ploidy_obs_total,
    n_ploidy_2N = raw$n_ploidy_2N,
    n_ploidy_4N = raw$n_ploidy_4N,
    cache_g_build = raw$cache_g_build,
    cache_g_hit = raw$cache_g_hit,
    cache_g_hysteresis = raw$cache_g_hysteresis,
    n_prior_terms = as.integer(.first_non_null_local(prior$n_terms, 0))
  )
}

# -----------------------------------------------------------------------------
# Function: evaluate_objective
# Purpose: Compute scalar objective value minimized by the optimizer.
# Parameters:
#   - par_transformed: Parameter vector in transformed optimization scale.
#   - scenarios: List of scenario-specific observation data and metadata.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
evaluate_objective <- function(par_transformed, scenarios, cfg) {
  evaluate_objective_components(par_transformed, scenarios = scenarios, cfg = cfg)$L
}

# -----------------------------------------------------------------------------
# Function: run_optimizer
# Purpose: Execute configured optimizer (DEoptim/optim) and return the best parameters.
# Parameters:
#   - objective_fn: Objective function returning scalar loss for optimizer.
#   - lower: Lower bounds for optimization variables.
#   - upper: Upper bounds for optimization variables.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - argv: Character vector of command-line arguments in --key=value format.
#   - stage_label: Text label for the current fitting stage.
#   - init_par: Optional warm-start transformed parameter vector.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
run_optimizer <- function(objective_fn, lower, upper, cfg, argv, stage_label = "fit", init_par = NULL, checkpoint_fn = NULL) {
  n_cores <- as.integer(max(1L, ifelse(is.finite(cfg$n_cores), cfg$n_cores, 1L)))
  use_deoptim <- isTRUE(cfg$use_deoptim)
  deoptim_parallel <- isTRUE(cfg$deoptim_parallel)
  interrupted <- FALSE
  iter_completed <- NA_integer_
  iter_target <- NA_integer_
# -----------------------------------------------------------------------------
# Function: fmt_int
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  fmt_int <- function(x) {
    if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
    as.character(as.integer(x))
  }
# -----------------------------------------------------------------------------
# Function: cluster_workers
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - cl: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  cluster_workers <- function(cl) {
    n <- tryCatch(length(cl), error = function(e) NA_integer_)
    if (is.na(n) || !is.finite(n)) return(NA_integer_)
    as.integer(max(0L, n))
  }
# -----------------------------------------------------------------------------
# Function: log_parallel_backend
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - solver: Function-specific input argument.
#   - requested: Function-specific input argument.
#   - started: Function-specific input argument.
#   - active: Function-specific input argument.
#   - extra: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  log_parallel_backend <- function(solver, requested, started, active, extra = NULL) {
    msg <- paste0(
      "[", stage_label, "] ", solver, " backend: requested_workers=", fmt_int(requested),
      ", started_workers=", fmt_int(started),
      ", active_workers=", fmt_int(active)
    )
    if (!is.null(extra) && nzchar(extra)) {
      msg <- paste0(msg, ", ", extra)
    }
    message(msg)
  }
# -----------------------------------------------------------------------------
# Function: init_cluster_workers
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - cl: Function-specific input argument.
#   - objective_fn: Objective function returning scalar loss for optimizer.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - stage_label: Text label for the current fitting stage.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  init_cluster_workers <- function(cl, objective_fn, cfg, stage_label) {
    if (is.null(cfg$model_path) || !nzchar(cfg$model_path) || !file.exists(cfg$model_path)) {
      stop("[", stage_label, "] Missing model_path for worker initialization.")
    }
    required_cpp <- c("cpp_o2simps_build_G_for_o2_triplet", "cpp_o2simps_simulate_one", "cpp_o2simps_objective_components_map")
    init_modes <- parallel::clusterCall(
      cl,
      function(model_path, wrapper_path, dll_path, required_cpp, stage_label) {
        Sys.setenv(
          OMP_NUM_THREADS = "1",
          MKL_NUM_THREADS = "1",
          OPENBLAS_NUM_THREADS = "1",
          MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path)
        )

        has_required <- function() {
          missing_cpp <- required_cpp[!vapply(required_cpp, exists, logical(1), mode = "function", inherits = TRUE)]
          length(missing_cpp) == 0L
        }
        has_expected_formals <- function() {
          if (!exists("cpp_o2simps_simulate_one", mode = "function", inherits = TRUE) ||
              !exists("cpp_o2simps_objective_components_map", mode = "function", inherits = TRUE)) {
            return(FALSE)
          }
          build_formals <- names(formals(get("cpp_o2simps_build_G_for_o2_triplet", mode = "function", inherits = TRUE)))
          sim_formals <- names(formals(get("cpp_o2simps_simulate_one", mode = "function", inherits = TRUE)))
          obj_formals <- names(formals(get("cpp_o2simps_objective_components_map", mode = "function", inherits = TRUE)))
          all(c("glucose", "p_wgd", "glucose_dynamic", "G_c", "G") %in% build_formals) &&
            ("sim_args" %in% sim_formals) &&
            all(c("scenario_data", "objective_data", "state_data", "sim_args") %in% obj_formals)
        }
        has_usable_backend <- function() {
          has_required() && has_expected_formals()
        }

        load_mode <- "unknown"
        last_err <- NULL
        loaded_ok <- FALSE

        # Path 1: use pre-generated wrapper/DLL from main process if both are visible.
        if (nzchar(wrapper_path) && nzchar(dll_path) && file.exists(wrapper_path) && file.exists(dll_path)) {
          max_attempts <- 20L
          for (attempt in seq_len(max_attempts)) {
            loaded_ok <- tryCatch({
              source(wrapper_path, local = .GlobalEnv)
              if (!has_usable_backend()) {
                stop("Wrapper backend missing required functions/formals after source(wrapper_path).")
              }
              TRUE
            }, error = function(e) {
              last_err <<- conditionMessage(e)
              FALSE
            })
            if (loaded_ok) {
              load_mode <- "shared_wrapper"
              break
            }
            if (!is.null(last_err) && grepl("unable to load shared object|cannot open shared object file|failed to map segment", last_err, ignore.case = TRUE)) {
              tryCatch(suppressWarnings(dyn.load(dll_path)), error = function(e2) {
                last_err <<- paste0(last_err, " | dyn.load retry failed: ", conditionMessage(e2))
              })
            }
            Sys.sleep(min(0.05 * attempt, 0.5))
          }
        }

        # Path 2: if shared wrapper path is unavailable on workers, compile/load in worker.
        if (!loaded_ok) {
          loaded_ok <- tryCatch({
            source(model_path, local = .GlobalEnv)
            if (!has_usable_backend()) {
              stop("Model backend missing required functions/formals after source(model_path).")
            }
            TRUE
          }, error = function(e) {
            if (is.null(last_err) || !nzchar(last_err)) {
              last_err <<- conditionMessage(e)
            } else {
              last_err <<- paste0(last_err, " | ", conditionMessage(e))
            }
            FALSE
          })
          if (loaded_ok) load_mode <- "worker_source_model"
        }

        if (!loaded_ok) {
          stop(
            "[", stage_label, "] Worker failed C++ initialization. ",
            "model_path=", model_path,
            "; wrapper_path=", wrapper_path,
            "; dll_path=", dll_path,
            "; last_error=", as.character(last_err)
          )
        }
        load_mode
      },
      cfg$model_path,
      cfg$cpp_wrapper_path,
      cfg$cpp_dll_path,
      required_cpp,
      stage_label
    )
    init_mode_tab <- sort(table(unlist(init_modes, use.names = FALSE)), decreasing = TRUE)
    message(
      "[", stage_label, "] worker C++ init modes: ",
      paste(paste0(names(init_mode_tab), "=", as.integer(init_mode_tab)), collapse = ", ")
    )
    parallel::clusterCall(cl, function(required_cpp) {
      missing <- required_cpp[!vapply(required_cpp, exists, logical(1), mode = "function", inherits = TRUE)]
      if (length(missing) > 0L) {
        stop("Worker missing required C++ wrapper functions: ", paste(missing, collapse = ", "))
      }
      build_formals <- names(formals(cpp_o2simps_build_G_for_o2_triplet))
      sim_formals <- names(formals(cpp_o2simps_simulate_one))
      obj_formals <- names(formals(cpp_o2simps_objective_components_map))
      if (!all(c("glucose", "p_wgd", "glucose_dynamic", "G_c", "G") %in% build_formals) ||
          !("sim_args" %in% sim_formals) ||
          !all(c("scenario_data", "objective_data", "state_data", "sim_args") %in% obj_formals)) {
        stop("Worker C++ wrappers are stale: required packed objective formals are missing.")
      }
      NULL
    }, required_cpp)
    export_global <- c(
      "evaluate_objective",
      "evaluate_objective_components",
      "evaluate_objective_components_raw",
      "build_model_core",
      "prepare_cpp_scenarios",
      "simulate_one",
      "make_init_state",
      "get_param_names",
      "decode_params",
      "compute_soft_prior_penalty",
      "as_num",
      "o2sd_first_non_null",
      "o2sd_as_bool_scalar",
      "canonical_ploidy_o2_death_mode",
      "canonical_start_with_mode",
      "assert_canonical_start_with_mode",
      "clip",
      ".first_non_null_local",
      "default_rho_2N_prior_bounds",
      "default_rho_2N_prior_center",
      "default_chr_lengths_bp_1to22",
      "normalize_chr_lengths_bp_1to22",
      "weighted_ploidy_from_total_N",
      "map_ploidy_to_N_by_chrlen",
      "cell_volume_mm3_by_ploidy",
      "cell_volume_mm3_by_chr_number",
      "cell_volume_mm3_by_N",
      "burden_volume_mm3_from_state"
    )
    export_env <- environment()
    export_global <- export_global[vapply(
      export_global,
      function(nm) exists(nm, envir = export_env, inherits = TRUE),
      logical(1)
    )]
    if (length(export_global) > 0) {
      parallel::clusterExport(cl, varlist = export_global, envir = export_env)
    }
    parallel::clusterExport(cl, varlist = c("objective_fn"), envir = environment())
    invisible(TRUE)
  }
  init_use <- NULL
  if (!is.null(init_par)) {
    init_use <- as.numeric(init_par[names(lower)])
    names(init_use) <- names(lower)
    if (any(!is.finite(init_use))) stop("[", stage_label, "] warm-start vector has missing/non-finite values.")
    init_use <- clip(init_use, lower, upper)
    message("[", stage_label, "] Using warm start for ", length(init_use), " parameters.")
  }
  emit_checkpoint <- function(best_par_t, best_val, iter_completed, iter_target, interrupted_flag = FALSE) {
    if (!is.function(checkpoint_fn) || is.null(best_par_t) || length(best_par_t) == 0L) return(invisible(NULL))
    try(
      checkpoint_fn(
        best_par_t = as.numeric(best_par_t),
        best_val = as.numeric(best_val),
        iter_completed = as.integer(iter_completed),
        iter_target = as.integer(iter_target),
        interrupted = isTRUE(interrupted_flag),
        stage_label = stage_label
      ),
      silent = TRUE
    )
    invisible(NULL)
  }

  deoptim_available <- requireNamespace("DEoptim", quietly = TRUE)
  deoptim_mode_ok <- (n_cores == 1L || deoptim_parallel)
  has_deoptim <- use_deoptim && deoptim_available && deoptim_mode_ok

  if (use_deoptim && !deoptim_available) {
    stop("[", stage_label, "] use_deoptim=TRUE but DEoptim is not available in this R environment.")
  }
  if (use_deoptim && !deoptim_mode_ok) {
    stop(
      "[", stage_label, "] use_deoptim=TRUE requires n_cores=1 or deoptim_parallel=TRUE; got n_cores=",
      n_cores, ", deoptim_parallel=", deoptim_parallel, "."
    )
  }

  if (has_deoptim) {
    de_requested <- n_cores
    de_started <- if (n_cores > 1L) NA_integer_ else 1L
    de_active <- 1L
    iter_target <- as.integer(cfg$itermax)
    de_init_plan <- .build_de_initialpop(
      np = cfg$NP,
      lower = lower,
      upper = upper,
      init_use = init_use,
      mode = .first_non_null_local(cfg$de_init_mode, "hybrid"),
      uniform_frac = .first_non_null_local(cfg$de_init_uniform_frac, 0.3),
      sigma_frac = .first_non_null_local(cfg$de_init_sigma_frac, 0.1)
    )
    de_extra <- paste0(
      "NP=", cfg$NP, ", itermax=", cfg$itermax,
      ", reltol=", signif(cfg$de_reltol, 6),
      ", steptol=", cfg$de_steptol,
      ", execution_mode=single_call",
      ", init_mode=", de_init_plan$mode_effective,
      ", warm_start=", if (isTRUE(de_init_plan$warm_start_used)) "TRUE" else "FALSE",
      ", init_local=", as.integer(de_init_plan$n_local),
      ", init_uniform=", as.integer(de_init_plan$n_uniform)
    )
    message(
      "[", stage_label, "] Starting DEoptim with itermax=", cfg$itermax,
      ", NP=", cfg$NP,
      ", reltol=", signif(cfg$de_reltol, 6),
      ", steptol=", cfg$de_steptol,
      ", execution_mode=single_call",
      ", n_cores=", n_cores
    )
    de_ctrl <- list(
      trace = TRUE,
      NP = cfg$NP,
      strategy = 2,
      reltol = cfg$de_reltol,
      steptol = cfg$de_steptol
    )
    with_realtime_deoptim_trace <- function(expr) {
      filter_cmd <- paste(
        "perl",
        "-ne",
        shQuote(
          "BEGIN { $| = 1; } s/\\s*bestmemit:.*$//; if (/Iteration:|bestvalit/) { print; }",
          type = "sh"
        )
      )
      sink_depth <- sink.number(type = "output")
      con <- tryCatch(pipe(filter_cmd, open = "w"), error = function(e) NULL)
      if (is.null(con)) {
        warning("Failed to start realtime DEoptim trace filter; falling back to unfiltered DEoptim trace.", call. = FALSE)
        return(force(expr))
      }
      on.exit({
        while (sink.number(type = "output") > sink_depth) sink(NULL, type = "output")
        try(close(con), silent = TRUE)
      }, add = TRUE)
      sink(con, type = "output")
      force(expr)
    }
    current_pop <- de_init_plan$pop
    best_par <- as.numeric(current_pop[1, ])
    names(best_par) <- names(lower)
    best_val <- suppressWarnings(tryCatch(as.numeric(objective_fn(best_par)), error = function(e) Inf))
    iter_completed <- 0L
    chunk_log <- list()
    emit_checkpoint(
      best_par_t = best_par,
      best_val = best_val,
      iter_completed = iter_completed,
      iter_target = iter_target,
      interrupted_flag = FALSE
    )
    if (n_cores > 1L) {
      message("[", stage_label, "] DEoptim parallel requested with n_cores=", n_cores, ".")
      cl <- tryCatch(
        parallel::makePSOCKcluster(n_cores),
        error = function(e) {
          stop("[", stage_label, "] Could not start parallel workers for DEoptim: ", conditionMessage(e))
        }
      )
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      de_started <- cluster_workers(cl)
      message(
        "[", stage_label, "] Initializing DEoptim workers C++ backend (workers=",
        de_started, "). This can take minutes on first compile."
      )
      t_init_workers <- Sys.time()
      tryCatch(
        init_cluster_workers(cl, objective_fn = objective_fn, cfg = cfg, stage_label = stage_label),
        error = function(e) {
          stop("[", stage_label, "] Failed to initialize DEoptim workers: ", conditionMessage(e))
        }
      )
      dt_init_workers <- as.numeric(difftime(Sys.time(), t_init_workers, units = "secs"))
      message("[", stage_label, "] Worker initialization completed in ", sprintf("%.1f", dt_init_workers), "s.")
      # IMPORTANT: when passing an explicit cluster, do NOT set
      # parallelType='parallel', otherwise DEoptim may try stopCluster(cl)
      # on an internal symbol 'cl' that does not exist in this code path.
      de_ctrl$cluster <- cl
      de_active <- ifelse(is.na(de_started) || de_started < 1L, 1L, de_started)
    } else {
      message("[", stage_label, "] DEoptim running in serial mode (n_cores=1).")
    }
    log_parallel_backend(
      solver = "DEoptim",
      requested = de_requested,
      started = de_started,
      active = de_active,
      extra = de_extra
    )
    de_ctrl$itermax <- iter_target
    de_ctrl$initialpop <- current_pop
    message(
      "[", stage_label, "] DEoptim single-call start: iter 1-", iter_target, "/", iter_target
    )
    chunk_res <- tryCatch(
      with_realtime_deoptim_trace(
        DEoptim::DEoptim(
          fn = objective_fn,
          lower = lower,
          upper = upper,
          control = de_ctrl
        )
      ),
      interrupt = function(e) {
        interrupted <<- TRUE
        NULL
      },
      error = function(e) {
        stop("[", stage_label, "] DEoptim failed: ", conditionMessage(e))
      }
    )
    if (isTRUE(interrupted) || is.null(chunk_res)) {
      emit_checkpoint(
        best_par_t = best_par,
        best_val = best_val,
        iter_completed = iter_completed,
        iter_target = iter_target,
        interrupted_flag = TRUE
      )
      message(
        "[", stage_label, "] Interrupt detected. Returning best-so-far from ",
        iter_completed, "/", iter_target, " completed DEoptim iterations."
      )
    } else {
      if (!is.null(chunk_res$member$pop) && is.matrix(chunk_res$member$pop)) {
        current_pop <- chunk_res$member$pop
      }
      iter_completed <- as.integer(.first_non_null_local(
        if (!is.null(chunk_res$optim$iter)) as.integer(chunk_res$optim$iter) else NULL,
        if (!is.null(chunk_res$member$bestmemit) && (is.matrix(chunk_res$member$bestmemit) || is.data.frame(chunk_res$member$bestmemit))) {
          as.integer(nrow(chunk_res$member$bestmemit))
        } else {
          NULL
        },
        iter_target
      ))
      if (!is.finite(iter_completed) || iter_completed < 0L) iter_completed <- as.integer(iter_target)
      if (is.finite(iter_target) && iter_completed > iter_target) iter_completed <- as.integer(iter_target)
      chunk_best_val <- suppressWarnings(as.numeric(.first_non_null_local(chunk_res$optim$bestval, Inf)))
      chunk_best_par <- as.numeric(.first_non_null_local(chunk_res$optim$bestmem, best_par))
      names(chunk_best_par) <- names(lower)
      if (is.finite(chunk_best_val) && chunk_best_val < best_val) {
        best_val <- chunk_best_val
        best_par <- chunk_best_par
      }
      chunk_log[[1L]] <- data.frame(
        chunk = 1L,
        iter_completed = as.integer(iter_completed),
        bestval = as.numeric(best_val),
        interrupted = FALSE,
        stringsAsFactors = FALSE
      )
      emit_checkpoint(
        best_par_t = best_par,
        best_val = best_val,
        iter_completed = iter_completed,
        iter_target = iter_target,
        interrupted_flag = FALSE
      )
    }
    if (!is.finite(best_val)) {
      best_val <- suppressWarnings(tryCatch(as.numeric(objective_fn(best_par)), error = function(e) Inf))
    }
    de_best_par <- best_par
    de_best_val <- best_val
    local_maxit <- as_int(.first_non_null_local(argv$local_optim_maxit, argv$optim_maxit, 200L), 200L)
    if (!is.finite(local_maxit) || is.na(local_maxit) || local_maxit < 1L) local_maxit <- 200L
    local_attempted <- FALSE
    local_accepted <- FALSE
    local_fit <- NULL
    local_best_val <- NA_real_
    local_convergence <- NA_integer_
    local_message <- NA_character_
    if (!isTRUE(interrupted) && is.finite(de_best_val)) {
      local_attempted <- TRUE
      message("[", stage_label, "] Starting L-BFGS-B local refinement from DEoptim best; maxit=", local_maxit, ".")
      local_fit <- tryCatch(
        suppressWarnings(
          optim(
            par = clip(best_par, lower, upper),
            fn = objective_fn,
            method = "L-BFGS-B",
            lower = lower,
            upper = upper,
            control = list(maxit = local_maxit)
          )
        ),
        error = function(e) {
          warning("[", stage_label, "] L-BFGS-B local refinement failed: ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (is.list(local_fit)) {
        local_best_val <- suppressWarnings(as.numeric(local_fit$value))
        local_convergence <- suppressWarnings(as.integer(local_fit$convergence))
        local_message <- as.character(.first_non_null_local(local_fit$message, NA_character_))
        if (is.finite(local_best_val) && local_best_val < de_best_val) {
          best_par <- as.numeric(local_fit$par)
          names(best_par) <- names(lower)
          best_par <- clip(best_par, lower, upper)
          best_val <- local_best_val
          local_accepted <- TRUE
          message(
            "[", stage_label, "] L-BFGS-B local refinement improved objective: ",
            signif(de_best_val, 8), " -> ", signif(best_val, 8), "."
          )
        } else {
          message(
            "[", stage_label, "] L-BFGS-B local refinement did not improve objective; keeping DEoptim best."
          )
        }
      }
    } else if (isTRUE(interrupted)) {
      message("[", stage_label, "] Skipping L-BFGS-B local refinement because DEoptim was interrupted.")
    }
    emit_checkpoint(
      best_par_t = best_par,
      best_val = best_val,
      iter_completed = iter_completed,
      iter_target = iter_target,
      interrupted_flag = isTRUE(interrupted)
    )
    de_method <- if (de_active > 1L) "DEoptim_parallel" else "DEoptim_serial"
    optimizer_method <- if (isTRUE(local_attempted)) {
      paste0(de_method, "_plus_LBFGSB_serial")
    } else {
      de_method
    }
    optim_res <- list(
      optim = list(bestmem = best_par, bestval = best_val),
      method = optimizer_method,
      deoptim = list(bestmem = de_best_par, bestval = de_best_val),
      local_refinement = list(
        method = "L-BFGS-B",
        attempted = isTRUE(local_attempted),
        accepted = isTRUE(local_accepted),
        maxit = as.integer(local_maxit),
        bestval = as.numeric(local_best_val),
        convergence = as.integer(local_convergence),
        message = local_message
      ),
      de_chunk_log = bind_rows(chunk_log),
      de_chunk_info = list(
        iter_chunk = as.integer(iter_target),
        iter_completed = as.integer(iter_completed),
        iter_target = as.integer(iter_target),
        interrupted = isTRUE(interrupted),
        execution_mode = "single_call"
      ),
      parallel_info = list(
      requested_workers = de_requested,
      started_workers = de_started,
      active_workers = de_active
      )
    )
    emit_checkpoint(
      best_par_t = best_par,
      best_val = best_val,
      iter_completed = iter_completed,
      iter_target = iter_target,
      interrupted_flag = isTRUE(interrupted)
    )
  }

  if (!use_deoptim) {
    n_starts <- as_int(argv$n_starts, 20L)
    maxit <- as_int(argv$optim_maxit, max(200L, cfg$itermax * 50L))
    trace_optim <- isTRUE(cfg$optim_trace)
    trace_every <- as.integer(max(1L, ifelse(is.finite(cfg$optim_trace_every), cfg$optim_trace_every, 1L)))
    opt_requested <- n_cores
    opt_started <- if (n_cores > 1L) NA_integer_ else 1L
    opt_active <- 1L
    message(
      "[", stage_label, "] Using multi-start optim (L-BFGS-B). starts=",
      n_starts, ", maxit=", maxit, ", n_cores=", n_cores
    )
    mid <- (lower + upper) / 2
    starts <- vector("list", n_starts)
    for (s in seq_len(n_starts)) {
      starts[[s]] <- if (s == 1L && !is.null(init_use)) {
        init_use
      } else if (s == 1L) {
        mid
      } else {
        stats::runif(length(lower), min = lower, max = upper)
      }
    }

# -----------------------------------------------------------------------------
# Function: worker_fit
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - p0: Function-specific input argument.
#   - objective_fn: Objective function returning scalar loss for optimizer.
#   - lower: Lower bounds for optimization variables.
#   - upper: Upper bounds for optimization variables.
#   - maxit: Maximum number of optimizer or IRLS iterations.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
    worker_fit <- function(p0, objective_fn, lower, upper, maxit) {
      fit <- optim(
        par = p0,
        fn = objective_fn,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = maxit)
      )
      list(par = fit$par, value = fit$value, convergence = fit$convergence)
    }

# -----------------------------------------------------------------------------
# Function: summarize_runs
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - run_log: Function-specific input argument.
#   - live_label: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
    summarize_runs <- function(run_log, live_label = "finished") {
      best_val <- Inf
      best_par <- NULL
      for (i in seq_along(run_log)) {
        fit <- run_log[[i]]
        if (is.finite(fit$value) && fit$value < best_val) {
          best_val <- fit$value
          best_par <- fit$par
        }
        if (trace_optim && (i %% trace_every == 0L || i == length(run_log))) {
          message(
            "[", stage_label, "] start ", i, "/", length(run_log), " ", live_label,
            ": val=", signif(fit$value, 6), ", best=", signif(best_val, 6)
          )
        }
      }
      list(best_val = best_val, best_par = best_par)
    }

# -----------------------------------------------------------------------------
# Function: run_starts_serial
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - starts: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
    run_starts_serial <- function(starts) {
      run_log <- vector("list", length(starts))
      for (i in seq_along(starts)) {
        run_log[[i]] <- worker_fit(
          starts[[i]],
          objective_fn = objective_fn,
          lower = lower,
          upper = upper,
          maxit = maxit
        )
      }
      stats <- summarize_runs(run_log, live_label = "finished")
      list(run_log = run_log, best_val = stats$best_val, best_par = stats$best_par)
    }

    if (n_cores > 1L) {
      effective_cap <- as.integer(max(1L, min(n_cores, length(starts))))
      if (length(starts) < n_cores) {
        message(
          "[", stage_label, "] Requested n_cores=", n_cores,
          " but only ", length(starts),
          " multi-start tasks exist; active workers will be capped by tasks."
        )
      }
      message(
        "[", stage_label, "] optim parallel requested with n_cores=", n_cores,
        ", n_starts=", length(starts),
        ", effective_worker_cap=", effective_cap, "."
      )
      cl <- tryCatch(
        parallel::makePSOCKcluster(n_cores),
        error = function(e) {
          stop("[", stage_label, "] Could not start parallel workers for optim: ", conditionMessage(e))
        }
      )
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
      opt_started <- cluster_workers(cl)
      tryCatch(
        init_cluster_workers(cl, objective_fn = objective_fn, cfg = cfg, stage_label = stage_label),
        error = function(e) {
          stop("[", stage_label, "] Failed to initialize optim workers: ", conditionMessage(e))
        }
      )
      opt_active <- as.integer(max(1L, min(length(starts), ifelse(is.na(opt_started), n_cores, opt_started))))
      if (trace_optim) {
        message("[", stage_label, "] Parallel backend does not stream per-start logs; progress will be reported after results are collected.")
      }
      run_log <- parallel::parLapplyLB(
        cl,
        starts,
        worker_fit,
        objective_fn = objective_fn,
        lower = lower,
        upper = upper,
        maxit = maxit
      )
      stats <- summarize_runs(run_log, live_label = "collected")
      best_val <- stats$best_val
      best_par <- stats$best_par
    } else {
      message("[", stage_label, "] optim running in serial mode (n_cores=1).")
      serial_out <- run_starts_serial(starts)
      run_log <- serial_out$run_log
      best_val <- serial_out$best_val
      best_par <- serial_out$best_par
    }
    log_parallel_backend(
      solver = "optim",
      requested = opt_requested,
      started = opt_started,
      active = opt_active,
      extra = paste0("n_starts=", length(starts), ", maxit=", maxit)
    )
    optim_res <- list(
      optim = list(bestmem = best_par, bestval = best_val),
      method = if (opt_active > 1L) "optim_L-BFGS-B_multistart_parallel" else "optim_L-BFGS-B_multistart",
      runs = run_log,
      parallel_info = list(
        requested_workers = opt_requested,
        started_workers = opt_started,
        active_workers = opt_active,
        n_starts = length(starts)
      )
    )
  }

  best_par <- as.numeric(best_par)
  names(best_par) <- names(lower)
  list(
    best_par = best_par,
    optim_res = optim_res,
    interrupted = isTRUE(interrupted),
    iter_completed = iter_completed,
    iter_target = iter_target
  )
}

# -----------------------------------------------------------------------------
# Function: collect_predictions
# Purpose: Generate fitted burden/ploidy predictions across all scenarios.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - scenarios: List of scenario-specific observation data and metadata.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
collect_predictions <- function(run_params, scenarios, cfg) {
  model_core <- if (!is.null(cfg$model_core)) cfg$model_core else build_model_core(cfg = cfg)
# -----------------------------------------------------------------------------
# Function: pred_one
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - sc: Function-specific input argument.
#   - i: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  pred_one <- function(sc, i) {
    sim <- simulate_one(run_params, sc, cfg, model_core = model_core)

    obs <- as.numeric(sc$obs_burden)
    pred_pop <- as.numeric(sim$Ntot_obs)
    pred_pop_live <- as.numeric(.first_non_null_local(sim$Ntot_live_obs, rep(NA_real_, length(pred_pop))))
    pred_pop_dead_hypoxia <- as.numeric(.first_non_null_local(sim$Ntot_dead_hypoxia_obs, rep(0, length(pred_pop))))
    pred_pop_dead_buffer <- as.numeric(.first_non_null_local(sim$Ntot_dead_buffer_obs, rep(0, length(pred_pop))))
    pred_pop_dead_total <- as.numeric(.first_non_null_local(sim$Ntot_dead_total_obs, pred_pop_dead_hypoxia + pred_pop_dead_buffer))
    pred_vol <- as.numeric(sim$Vmm3_obs)
    pred_vol_live <- as.numeric(.first_non_null_local(sim$Vmm3_live_obs, rep(NA_real_, length(pred_vol))))
    pred_vol_dead_hypoxia <- as.numeric(.first_non_null_local(sim$Vmm3_dead_hypoxia_obs, rep(0, length(pred_vol))))
    pred_vol_dead_buffer <- as.numeric(.first_non_null_local(sim$Vmm3_dead_buffer_obs, rep(0, length(pred_vol))))
    pred_vol_dead_total <- as.numeric(.first_non_null_local(sim$Vmm3_dead_total_obs, pred_vol_dead_hypoxia + pred_vol_dead_buffer))
    obs_delta <- obs - obs[1]
    pred_delta <- pred_vol - pred_vol[1]
    s_obs <- max(abs(obs_delta), na.rm = TRUE)
    s_pred <- max(abs(pred_delta), na.rm = TRUE)
    obs_norm <- if (s_obs > 0) obs_delta / s_obs else obs_delta
    pred_norm <- if (s_pred > 0) pred_delta / s_pred else pred_delta
    log_eps <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 1e-15)

    burden_df <- data.frame(
      harvest = sc$harvest,
      cohort = sc$cohort,
      dose = sc$dose,
      treat_day = sc$treat_day,
      day = sc$obs_days,
      obs_burden = obs,
      pred_pop = pred_pop,
      pred_pop_live = pred_pop_live,
      pred_pop_dead_hypoxia = pred_pop_dead_hypoxia,
      pred_pop_dead_buffer = pred_pop_dead_buffer,
      pred_pop_dead_total = pred_pop_dead_total,
      pred_burden_volume_mm3 = pred_vol,
      pred_burden_live_volume_mm3 = pred_vol_live,
      pred_burden_dead_hypoxia_volume_mm3 = pred_vol_dead_hypoxia,
      pred_burden_dead_buffer_volume_mm3 = pred_vol_dead_buffer,
      pred_burden_dead_total_volume_mm3 = pred_vol_dead_total,
      obs_log_burden = ifelse(is.finite(obs) & obs >= 0, log(pmax(obs, log_eps)), NA_real_),
      pred_log_burden = ifelse(is.finite(pred_vol) & pred_vol >= 0, log(pmax(pred_vol, log_eps)), NA_real_),
      obs_norm = obs_norm,
      pred_norm = pred_norm
    )

    ploidy_df <- NULL
    if (length(sc$ploidy_obs_z) > 0) {
      obs_N <- if (identical(assert_canonical_start_with_mode(cfg$start_with), "chr_number")) {
        suppressWarnings(as.numeric(sc$chr_number_obs))
      } else {
        map_ploidy_to_N_by_chrlen(
          ploidy_values = as.numeric(sc$ploidy_obs_z) / as.numeric(cfg$N_UNIT),
          N_grid = cfg$N_MIN:cfg$N_MAX,
          chr_lengths_bp = cfg$chr_lengths_bp
        )
      }
      obs_N <- as.integer(round(clip(obs_N, cfg$N_MIN, cfg$N_MAX)))
      obs_N <- obs_N[is.finite(obs_N)]
      obs_tab <- table(obs_N)
      ploidy_df <- data.frame(
        harvest = sc$harvest,
        cohort = sc$cohort,
        dose = sc$dose,
        N = as.integer(names(sim$frac_N)),
        pred_fraction = as.numeric(sim$frac_N)
      )
      ploidy_df$obs_count <- as.numeric(obs_tab[as.character(ploidy_df$N)])
      ploidy_df$obs_count[is.na(ploidy_df$obs_count)] <- 0
    }

    list(burden = burden_df, ploidy = ploidy_df)
  }

  pred_rows <- map_scenarios_parallel(
    scenarios = scenarios,
    n_cores = as.integer(.first_non_null_local(cfg$predict_n_cores, 1L)),
    label = "predict",
    fn = pred_one
  )

  list(
    burden = bind_rows(lapply(pred_rows, function(x) x$burden)),
    ploidy = bind_rows(Filter(Negate(is.null), lapply(pred_rows, function(x) x$ploidy)))
  )
}

# -----------------------------------------------------------------------------
# Function: main_fit_single_seed
# Purpose: Run a single-seed fit with fully resolved CLI/runtime arguments.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
main_fit_single_seed <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  # Accept legacy uppercase key used by some YAML files.
  if (is.null(argv$dt) && !is.null(argv$DT)) argv$dt <- argv$DT
  # Accept --parameters as a higher-priority alias for the parameter table path.
  if (!is.null(argv$parameters) && nzchar(trimws(as.character(argv$parameters)))) {
    argv$parameter_table <- argv$parameters
  }
  script_dir <- get_script_dir()
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  source(file.path(workflow_root, "util", "o2g_supply_demand_map_common_semantics.R"), local = environment())
  default_glucose_use <- canonical_glucose_enabled(.first_non_null_local(argv$glucose, TRUE), TRUE)
  default_parameter_table <- default_o2g_parameter_table_path_common(
    script_dir = script_dir,
    glucose = default_glucose_use,
    must_exist = TRUE
  )
  if (is.null(argv$parameter_table) || !nzchar(trimws(as.character(argv$parameter_table)))) {
    argv$parameter_table <- default_parameter_table
  }
  require_cli_args(argv, c(
    "data_dir", "n_cores", "use_deoptim", "deoptim_parallel",
    "fit_treatment", "dose_zero_only", "paired_only", "truncate_at_treatment", "ploidy_at_harvest",
    "itermax", "NP", "n_starts", "optim_maxit",
    "sigma_ploidy", "burden_log_eps", "burden_exclude_day0",
    "use_soft_prior", "lambda_prior",
    "fit_tau_O2",
    "start_with",
    "ploidy_O2_death", "o2_Nref", "o2_min",
    "parameter_table",
    "Crowding", "K", "crowding",
    "N_UNIT", "N_MIN", "N_MAX", "dt", "o2_burden_feedback",
    "o2_cache_bin_pct", "o2_cache_hysteresis_pct", "o2_cache_profile",
    "init_total_size", "dose_ref", "tx_mult_min", "min_pop",
    "prior_center_log10_k_o", "prior_sd_log10_k_o",
    "prior_center_log10_kappa_O", "prior_sd_log10_kappa_O",
    "prior_center_log10_o2_S0", "prior_sd_log10_o2_S0",
    "prior_center_log10_eta_o2", "prior_sd_log10_eta_o2",
    "prior_center_log10_gamma_loss", "prior_sd_log10_gamma_loss",
    "prior_center_log10_rho_2N", "prior_sd_log10_rho_2N",
    "prior_center_log10_mu_hp", "prior_sd_log10_mu_hp",
    "prior_center_gamma_mu", "prior_sd_gamma_mu",
    "prior_center_log10_k_clear", "prior_sd_log10_k_clear",
    "optim_trace", "optim_trace_every", "trace_obj",
    "de_init_mode", "de_init_uniform_frac", "de_init_sigma_frac", "de_reltol", "de_steptol",
    "predict_n_cores", "max_scenarios", "seed"
  ))
  model_path <- file.path(workflow_root, "model", "model_O2G_supply_demand_MAP.R")
  if (!file.exists(model_path)) stop("Cannot find model_O2G_supply_demand_MAP.R at ", model_path)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
  source(model_path)
  required_cpp_fit <- c("cpp_o2simps_build_G_for_o2_triplet", "cpp_o2simps_simulate_one", "cpp_o2simps_objective_components_map")
  missing_cpp_fit <- required_cpp_fit[!vapply(required_cpp_fit, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing_cpp_fit) > 0L) {
    stop("Required C++ symbols missing for fit path: ", paste(missing_cpp_fit, collapse = ", "))
  }
  cpp_dll <- tryCatch(
    o2simps_cpp_dll_info(),
    error = function(e) stop("Failed to resolve compiled O2G_supply_demand_MAP DLL info: ", conditionMessage(e))
  )

  default_data_dir <- normalizePath(file.path(script_dir, "..", "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
  data_dir <- if (!is.null(argv$data_dir)) argv$data_dir else default_data_dir
  truncate_at_treatment <- as_bool(argv$truncate_at_treatment, FALSE)
  n_cores_arg <- as_int(argv$n_cores, NA_integer_)
  n_cores_use <- if (is.finite(n_cores_arg)) n_cores_arg else default_n_cores()
  o2_S0_upper_arg <- read_o2_S0_natural_upper_bound_common(argv$parameter_table, fallback = NA_real_)
  if (!is.finite(o2_S0_upper_arg) || o2_S0_upper_arg <= 0) {
    stop("Failed to read natural-scale o2_S0 upper bound from parameter_table (param_prototype == 'o2_S0').")
  }

  cfg <- list(
    # model constants
    model_path = model_path,
    cpp_dll_name = as.character(cpp_dll$name),
    cpp_dll_path = as.character(cpp_dll$path),
    cpp_wrapper_path = as.character(cpp_dll$wrapper_path),
    N_UNIT = as_int(argv$N_UNIT, 22L),
    chr_lengths_bp = default_chr_lengths_bp_1to22(),
    N_MIN = as_int(argv$N_MIN, 22L),
    N_MAX = as_int(argv$N_MAX, 154L),
    DT = as_num(argv$dt, 0.5),
    o2_S0_upper_bound = o2_S0_upper_arg,
    ploidy_O2_death = canonical_ploidy_o2_death_mode(argv$ploidy_O2_death, "diploid_NULL"),
    misseg_loss_survival = canonical_misseg_loss_survival_mode(.first_non_null_local(argv$misseg_loss_survival, "nullisomy"), "nullisomy"),
    glucose = canonical_glucose_enabled(.first_non_null_local(argv$glucose, TRUE), TRUE),
    glucose_dynamic = canonical_glucose_dynamic(.first_non_null_local(argv$glucose_dynamic, FALSE), FALSE),
    glucose_ref_mM = as_num(argv$glucose_ref_mM, default_glucose_ref_mM()),
    start_with = canonical_start_with_mode(argv$start_with, "ploidy"),
    o2_burden_feedback = as_bool(argv$o2_burden_feedback, TRUE),
    O2_growth = as_bool(argv$O2_growth, TRUE),
    o2_cache_bin_pct = as_num(argv$o2_cache_bin_pct, 0.01),
    o2_cache_hysteresis_pct = as_num(argv$o2_cache_hysteresis_pct, 0.005),
    o2_cache_profile = as_bool(argv$o2_cache_profile, FALSE),
    o2_Nref = as_num(argv$o2_Nref, as_num(argv$init_total_size, 1e6)),
    o2_min = as_num(argv$o2_min, 0.5),
    Crowding = o2sd_as_bool_scalar(argv$Crowding, TRUE),
    fit_tau_O2 = as_bool(argv$fit_tau_O2, FALSE),
    parameter_table = if (!is.null(argv$parameter_table)) argv$parameter_table else default_parameter_table,
    K = as_num(argv$K, 1e12),
    crowding = if (!is.null(argv$crowding)) argv$crowding else "logistic",
    init_total_size = as_num(argv$init_total_size, 1e6),
    dose_ref = as_num(argv$dose_ref, 30),
    tx_mult_min = as_num(argv$tx_mult_min, 0.05),
    min_pop = as_num(argv$min_pop, 1e-12),
    # objective settings (MAP likelihood)
    sigma_ploidy = as_num(argv$sigma_ploidy, 0.08),
    burden_log_eps = as_num(argv$burden_log_eps, 1e-12),
    burden_exclude_day0 = as_bool(argv$burden_exclude_day0, TRUE),
    use_soft_prior = as_bool(argv$use_soft_prior, TRUE),
    lambda_prior = as_num(argv$lambda_prior, 0.1),
    harvest_init_multiplier = as_bool(argv$harvest_init_multiplier, FALSE),
    prior_center_log_init_mult = as_num(argv$prior_center_log_init_mult, 0.0),
    prior_sd_log_init_mult = as_num(argv$prior_sd_log_init_mult, 0.35),
    log_init_mult_lower = as_num(argv$log_init_mult_lower, -1.0),
    log_init_mult_upper = as_num(argv$log_init_mult_upper, 1.0),
    prior_center_log10_k_o = as_num(argv$prior_center_log10_k_o, log10(50)),
    prior_sd_log10_k_o = as_num(argv$prior_sd_log10_k_o, 1.0),
    prior_center_log10_kappa_O = as_num(argv$prior_center_log10_kappa_O, NA_real_),
    prior_sd_log10_kappa_O = as_num(argv$prior_sd_log10_kappa_O, 1.0),
    prior_center_log10_o2_S0 = as_num(argv$prior_center_log10_o2_S0, NA_real_),
    prior_sd_log10_o2_S0 = as_num(argv$prior_sd_log10_o2_S0, 0.5),
    prior_center_log10_eta_o2 = as_num(argv$prior_center_log10_eta_o2, NA_real_),
    prior_sd_log10_eta_o2 = as_num(argv$prior_sd_log10_eta_o2, 0.5),
    prior_center_log10_gamma_loss = as_num(argv$prior_center_log10_gamma_loss, log10(0.1)),
    prior_sd_log10_gamma_loss = as_num(argv$prior_sd_log10_gamma_loss, 0.5),
    prior_center_buffer_smax = as_num(argv$prior_center_buffer_smax, 0.8),
    prior_sd_buffer_smax = as_num(argv$prior_sd_buffer_smax, 0.25),
    prior_center_log10_buffer_beta = as_num(argv$prior_center_log10_buffer_beta, 0.0),
    prior_sd_log10_buffer_beta = as_num(argv$prior_sd_log10_buffer_beta, 0.75),
    prior_center_log10_buffer_n_exp = as_num(argv$prior_center_log10_buffer_n_exp, 0.0),
    prior_sd_log10_buffer_n_exp = as_num(argv$prior_sd_log10_buffer_n_exp, 0.75),
    prior_center_log10_rho_2N = as_num(argv$prior_center_log10_rho_2N, NA_real_),
    prior_sd_log10_rho_2N = as_num(argv$prior_sd_log10_rho_2N, 0.35),
    prior_center_log10_mu_hp = as_num(argv$prior_center_log10_mu_hp, NA_real_),
    prior_sd_log10_mu_hp = as_num(argv$prior_sd_log10_mu_hp, 1.0),
    prior_center_gamma_mu = as_num(argv$prior_center_gamma_mu, NA_real_),
    prior_sd_gamma_mu = as_num(argv$prior_sd_gamma_mu, 0.5),
    prior_center_log10_k_clear = as_num(argv$prior_center_log10_k_clear, NA_real_),
    prior_sd_log10_k_clear = as_num(argv$prior_sd_log10_k_clear, 1.0),
    optim_trace = as_bool(argv$optim_trace, TRUE),
    optim_trace_every = as_int(argv$optim_trace_every, 1L),
    trace_obj = as_bool(argv$trace_obj, FALSE),
    # fitting options
    fit_treatment = as_bool(argv$fit_treatment, FALSE),
    dose_zero_only = as_bool(argv$dose_zero_only, TRUE),
    paired_only = as_bool(argv$paired_only, TRUE),
    truncate_at_treatment = truncate_at_treatment,
    ploidy_at_harvest = as_bool(argv$ploidy_at_harvest, TRUE),
    use_deoptim = as_bool(argv$use_deoptim, TRUE),
    deoptim_parallel = as_bool(argv$deoptim_parallel, FALSE),
    de_init_mode = tolower(trimws(as.character(.first_non_null_local(argv$de_init_mode, "hybrid")))),
    de_init_uniform_frac = as_num(argv$de_init_uniform_frac, 0.3),
    de_init_sigma_frac = as_num(argv$de_init_sigma_frac, 0.1),
    de_reltol = as_num(argv$de_reltol, 1e-4),
    de_steptol = as_int(argv$de_steptol, 25L),
    itermax = as_int(argv$itermax, 40L),
    NP = as_int(argv$NP, 80L),
    n_cores = n_cores_use,
    # Keep post-fit prediction stable: use serial by default unless explicitly requested.
    predict_n_cores = as_int(argv$predict_n_cores, 1L),
    seed = as_int(argv$seed, 1L),
    max_scenarios = as_num(argv$max_scenarios, Inf)
  )

  cfg <- normalize_sim_cfg_common(cfg, context = "fit")
  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- resolve_terminal_ploidy_path(data_dir)
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)
  cfg$harvest_param_ids <- unique(vapply(scenarios, function(sc) as.character(sc$harvest), character(1)))
  param_bundle <- build_transformed_parameter_table(
    path = cfg$parameter_table,
    fit_treatment = cfg$fit_treatment,
    fit_tau_O2 = cfg$fit_tau_O2,
    O2_growth = cfg$O2_growth,
    glucose = cfg$glucose,
    glucose_dynamic = cfg$glucose_dynamic,
    misseg_loss_survival = cfg$misseg_loss_survival,
    harvest_init_multiplier = cfg$harvest_init_multiplier,
    harvest_ids = cfg$harvest_param_ids,
    prior_center_log_init_mult = cfg$prior_center_log_init_mult,
    log_init_mult_lower = cfg$log_init_mult_lower,
    log_init_mult_upper = cfg$log_init_mult_upper
  )
  cfg <- sync_cfg_from_natural_parameter_table(cfg, param_bundle$natural)
  cfg <- finalize_prior_defaults(cfg)

  if (!is.finite(cfg$o2_S0_upper_bound) || cfg$o2_S0_upper_bound <= 0) {
    stop("Failed to read natural-scale o2_S0 upper bound from parameter_table.")
  }
  if (!is.finite(cfg$o2_min) || cfg$o2_min < 0) cfg$o2_min <- 0.5
  cfg$o2_min <- min(max(cfg$o2_min, 0), cfg$o2_S0_upper_bound)

  if (!is.finite(cfg$K) || cfg$K <= 0) stop("K must be > 0")
  if (!cfg$crowding %in% c("logistic", "gompertz")) stop("crowding must be logistic or gompertz")
  if (cfg$DT <= 0) stop("dt must be > 0")
  if (cfg$N_MAX < cfg$N_MIN) stop("N_MAX must be >= N_MIN")
  if (!is.finite(cfg$o2_cache_bin_pct) || cfg$o2_cache_bin_pct <= 0 || cfg$o2_cache_bin_pct > 100) {
    stop("o2_cache_bin_pct must be in (0,100].")
  }
  if (!is.finite(cfg$o2_cache_hysteresis_pct) || cfg$o2_cache_hysteresis_pct < 0 || cfg$o2_cache_hysteresis_pct > 100) {
    stop("o2_cache_hysteresis_pct must be in [0,100].")
  }
  if (!is.finite(cfg$o2_Nref) || cfg$o2_Nref <= 0) stop("o2_Nref must be > 0")
  if (!is.finite(cfg$p_mis_base) || cfg$p_mis_base < 0 || cfg$p_mis_base > 1) stop("p_mis_base must be in [0,1]")
  if (!is.finite(cfg$o2_S0_init) || cfg$o2_S0_init <= 0 || cfg$o2_S0_init > cfg$o2_S0_upper_bound) {
    stop("o2_S0_init must be in (0, o2_S0_upper_bound]")
  }
  if (!is.finite(cfg$o2_S0_min) || cfg$o2_S0_min <= 0) stop("o2_S0_min must be > 0")
  if (!is.finite(cfg$o2_S0_max) || cfg$o2_S0_max <= 0) stop("o2_S0_max must be > 0")
  if (cfg$o2_S0_max > cfg$o2_S0_upper_bound) stop("o2_S0_max must be <= o2_S0_upper_bound")
  if (cfg$o2_S0_max < cfg$o2_S0_min) stop("o2_S0_max must be >= o2_S0_min")
  if (cfg$o2_min > cfg$o2_S0_upper_bound) stop("o2_min must be <= o2_S0_upper_bound")
  if (!is.finite(cfg$kappa_O_init) || cfg$kappa_O_init <= 0) stop("kappa_O_init must be > 0")
  if (!is.finite(cfg$kappa_O_min) || cfg$kappa_O_min <= 0) stop("kappa_O_min must be > 0")
  if (!is.finite(cfg$kappa_O_max) || cfg$kappa_O_max <= 0) stop("kappa_O_max must be > 0")
  if (cfg$kappa_O_max < cfg$kappa_O_min) stop("kappa_O_max must be >= kappa_O_min")
  if (!is.finite(cfg$eta_o2_init) || cfg$eta_o2_init <= 0) stop("eta_o2_init must be > 0")
  if (!is.finite(cfg$eta_o2_min) || cfg$eta_o2_min <= 0) stop("eta_o2_min must be > 0")
  if (!is.finite(cfg$eta_o2_max) || cfg$eta_o2_max <= 0) stop("eta_o2_max must be > 0")
  if (cfg$eta_o2_max < cfg$eta_o2_min) stop("eta_o2_max must be >= eta_o2_min")
  if (!is.finite(cfg$o2_crit_init) || cfg$o2_crit_init < 0) stop("o2_crit_init must be >= 0")
  if (!is.finite(cfg$o2_crit_min) || cfg$o2_crit_min < 0) stop("o2_crit_min must be >= 0")
  if (!is.finite(cfg$o2_crit_max) || cfg$o2_crit_max <= 0) stop("o2_crit_max must be > 0")
  if (cfg$o2_crit_max < cfg$o2_crit_min) stop("o2_crit_max must be >= o2_crit_min")
  if (!is.finite(cfg$n_O_init) || cfg$n_O_init < 0) stop("n_O_init must be >= 0")
  if (!is.finite(cfg$n_O_min) || cfg$n_O_min < 0) stop("n_O_min must be >= 0")
  if (!is.finite(cfg$n_O_max) || cfg$n_O_max < 0) stop("n_O_max must be >= 0")
  if (cfg$n_O_max < cfg$n_O_min) stop("n_O_max must be >= n_O_min")
  if (!is.finite(cfg$alpha_o2_init) || cfg$alpha_o2_init <= 0) stop("alpha_o2_init must be > 0")
  if (!is.finite(cfg$alpha_o2_min) || cfg$alpha_o2_min <= 0) stop("alpha_o2_min must be > 0")
  if (!is.finite(cfg$alpha_o2_max) || cfg$alpha_o2_max <= 0) stop("alpha_o2_max must be > 0")
  if (cfg$alpha_o2_max < cfg$alpha_o2_min) stop("alpha_o2_max must be >= alpha_o2_min")
  if (!is.finite(cfg$gamma_growth_init) || cfg$gamma_growth_init <= 0) stop("gamma_growth_init must be > 0")
  if (!is.finite(cfg$gamma_growth_min) || cfg$gamma_growth_min <= 0) stop("gamma_growth_min must be > 0")
  if (!is.finite(cfg$gamma_growth_max) || cfg$gamma_growth_max <= 0) stop("gamma_growth_max must be > 0")
  if (cfg$gamma_growth_max < cfg$gamma_growth_min) stop("gamma_growth_max must be >= gamma_growth_min")
  if (!is.finite(cfg$gamma_mu_init) || cfg$gamma_mu_init <= 0) stop("gamma_mu_init must be > 0")
  if (!is.finite(cfg$gamma_mu_min) || cfg$gamma_mu_min <= 0) stop("gamma_mu_min must be > 0")
  if (!is.finite(cfg$gamma_mu_max) || cfg$gamma_mu_max <= 0) stop("gamma_mu_max must be > 0")
  if (cfg$gamma_mu_max < cfg$gamma_mu_min) stop("gamma_mu_max must be >= gamma_mu_min")
  if (!is.finite(cfg$mu_hp_init) || cfg$mu_hp_init <= 0) stop("mu_hp_init must be > 0")
  if (!is.finite(cfg$mu_hp_min) || cfg$mu_hp_min <= 0) stop("mu_hp_min must be > 0")
  if (!is.finite(cfg$mu_hp_max) || cfg$mu_hp_max <= 0) stop("mu_hp_max must be > 0")
  if (cfg$mu_hp_max < cfg$mu_hp_min) stop("mu_hp_max must be >= mu_hp_min")
  if (!is.finite(cfg$k_clear_init) || cfg$k_clear_init <= 0) stop("k_clear_init must be > 0")
  if (!is.finite(cfg$k_clear_min) || cfg$k_clear_min <= 0) stop("k_clear_min must be > 0")
  if (!is.finite(cfg$k_clear_max) || cfg$k_clear_max <= 0) stop("k_clear_max must be > 0")
  if (cfg$k_clear_max < cfg$k_clear_min) stop("k_clear_max must be >= k_clear_min")
  if (!is.finite(cfg$tau_O2_init) || cfg$tau_O2_init <= 0) stop("tau_O2_init must be > 0")
  if (!is.finite(cfg$tau_O2_min) || cfg$tau_O2_min <= 0) stop("tau_O2_min must be > 0")
  if (!is.finite(cfg$tau_O2_max) || cfg$tau_O2_max <= 0) stop("tau_O2_max must be > 0")
  if (cfg$tau_O2_max < cfg$tau_O2_min) stop("tau_O2_max must be >= tau_O2_min")
  if (!isTRUE(cfg$fit_tau_O2)) {
    if (!is.finite(cfg$tau_O2) || cfg$tau_O2 <= 0) stop("tau_O2 must be > 0 when provided.")
  } else {
    cfg$tau_O2 <- as.numeric(.first_non_null_local(cfg$tau_O2, cfg$tau_O2_init))
    if (!is.finite(cfg$tau_O2) || cfg$tau_O2 <= 0) {
      stop("tau_O2 must be > 0 when fit_tau_O2=TRUE (used as init/default only).")
    }
  }
  if (cfg$itermax < 1) stop("itermax must be >= 1")
  if (cfg$n_cores < 1) stop("n_cores must be >= 1")
  if (cfg$optim_trace_every < 1) stop("optim_trace_every must be >= 1")
  if (!cfg$de_init_mode %in% c("hybrid", "uniform")) stop("de_init_mode must be one of: hybrid, uniform")
  if (!is.finite(cfg$de_init_uniform_frac) || cfg$de_init_uniform_frac < 0 || cfg$de_init_uniform_frac > 1) {
    stop("de_init_uniform_frac must be in [0,1]")
  }
  if (!is.finite(cfg$de_init_sigma_frac) || cfg$de_init_sigma_frac <= 0) {
    stop("de_init_sigma_frac must be > 0")
  }
  if (!is.finite(cfg$de_reltol) || cfg$de_reltol <= 0) stop("de_reltol must be > 0")
  if (is.na(cfg$de_steptol) || cfg$de_steptol < 1L) stop("de_steptol must be >= 1")
  if (is.null(cfg$cpp_dll_name) || !nzchar(cfg$cpp_dll_name)) stop("cpp_dll_name must be set")
  if (is.null(cfg$cpp_dll_path) || !nzchar(cfg$cpp_dll_path) || !file.exists(cfg$cpp_dll_path)) {
    stop("cpp_dll_path must exist and be readable: ", as.character(cfg$cpp_dll_path))
  }
  if (is.null(cfg$cpp_wrapper_path) || !nzchar(cfg$cpp_wrapper_path) || !file.exists(cfg$cpp_wrapper_path)) {
    stop("cpp_wrapper_path must exist and be readable: ", as.character(cfg$cpp_wrapper_path))
  }
  if (!is.finite(cfg$burden_log_eps) || cfg$burden_log_eps <= 0) stop("burden_log_eps must be > 0")
  if (!is.finite(cfg$sigma_burden) || cfg$sigma_burden <= 0) stop("sigma_burden must be > 0")
  if (!is.finite(cfg$sigma_burden_min) || cfg$sigma_burden_min <= 0) stop("sigma_burden_min must be > 0")
  if (!is.finite(cfg$sigma_burden_max) || cfg$sigma_burden_max <= 0) stop("sigma_burden_max must be > 0")
  if (cfg$sigma_burden_max < cfg$sigma_burden_min) stop("sigma_burden_max must be >= sigma_burden_min")
  if (!is.finite(cfg$sigma_ploidy) || cfg$sigma_ploidy <= 0) stop("sigma_ploidy must be > 0")
  if (!is.finite(cfg$rho_2N_min) || cfg$rho_2N_min <= 0) stop("rho_2N_min must be > 0")
  if (!is.finite(cfg$rho_2N_max) || cfg$rho_2N_max <= 0) stop("rho_2N_max must be > 0")
  if (cfg$rho_2N_max < cfg$rho_2N_min) stop("rho_2N_max must be >= rho_2N_min")
  if (!is.finite(cfg$lambda_prior) || cfg$lambda_prior < 0) stop("lambda_prior must be >= 0")
  if (!is.finite(cfg$prior_sd_log_init_mult) || cfg$prior_sd_log_init_mult <= 0) stop("prior_sd_log_init_mult must be > 0")
  if (!is.finite(cfg$log_init_mult_lower) || !is.finite(cfg$log_init_mult_upper)) stop("log_init_mult bounds must be finite")
  if (cfg$log_init_mult_upper < cfg$log_init_mult_lower) stop("log_init_mult_upper must be >= log_init_mult_lower")
  if ((!isTRUE(cfg$glucose) || !isTRUE(cfg$glucose_dynamic)) &&
      (!is.finite(cfg$prior_sd_log10_k_o) || cfg$prior_sd_log10_k_o <= 0)) {
    stop("prior_sd_log10_k_o must be > 0")
  }
  if (!is.finite(cfg$prior_sd_log10_kappa_O) || cfg$prior_sd_log10_kappa_O <= 0) stop("prior_sd_log10_kappa_O must be > 0")
  if (!is.finite(cfg$prior_sd_log10_o2_S0) || cfg$prior_sd_log10_o2_S0 <= 0) stop("prior_sd_log10_o2_S0 must be > 0")
  if (!is.finite(cfg$prior_sd_log10_eta_o2) || cfg$prior_sd_log10_eta_o2 <= 0) stop("prior_sd_log10_eta_o2 must be > 0")
  if (!is.finite(cfg$prior_sd_log10_gamma_loss) || cfg$prior_sd_log10_gamma_loss <= 0) stop("prior_sd_log10_gamma_loss must be > 0")
  if (!is.finite(cfg$prior_sd_log10_rho_2N) || cfg$prior_sd_log10_rho_2N <= 0) stop("prior_sd_log10_rho_2N must be > 0")
  if (!is.finite(cfg$prior_sd_log10_mu_hp) || cfg$prior_sd_log10_mu_hp <= 0) stop("prior_sd_log10_mu_hp must be > 0")
  if (!is.finite(cfg$prior_sd_gamma_mu) || cfg$prior_sd_gamma_mu <= 0) stop("prior_sd_gamma_mu must be > 0")
  if (!is.finite(cfg$prior_sd_log10_k_clear) || cfg$prior_sd_log10_k_clear <= 0) stop("prior_sd_log10_k_clear must be > 0")
  if (!cfg$use_deoptim && cfg$deoptim_parallel) stop("deoptim_parallel=TRUE requires use_deoptim=TRUE")

  append_ts_out_dir <- as_bool(argv$append_timestamp_out_dir, FALSE)
  ts_format <- if (!is.null(argv$timestamp_format)) argv$timestamp_format else "%Y%m%d_%H%M%S"
  run_stamp <- format(Sys.time(), ts_format)
  out_dir <- if (!is.null(argv$out_dir)) {
    if (append_ts_out_dir) {
      paste0(argv$out_dir, "_", run_stamp)
    } else {
      argv$out_dir
    }
  } else {
    file.path(workflow_root, "..", "..", "results", paste0("fit_model_O2G_supply_demand_MAP_", run_stamp))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- file.path(out_dir, "checkpoints")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    param_bundle$transformed_output,
    file = file.path(out_dir, "parameter_table.csv"),
    sep = ",",
    quote = FALSE,
    row.names = FALSE
  )
  if (grepl("^seed[0-9]+$", basename(out_dir))) {
    write.table(
      param_bundle$transformed_output,
      file = file.path(dirname(out_dir), "parameter_table.csv"),
      sep = ",",
      quote = FALSE,
      row.names = FALSE
    )
  }

  cfg$model_core <- build_model_core(cfg = cfg)
  cfg$scenario_cpp <- prepare_cpp_scenarios(scenarios, cfg)

  full_names <- get_param_names(
    fit_treatment = isTRUE(cfg$fit_treatment),
    fit_tau_O2 = isTRUE(cfg$fit_tau_O2),
    glucose = isTRUE(cfg$glucose),
    glucose_dynamic = isTRUE(cfg$glucose_dynamic),
    harvest_init_multiplier = isTRUE(cfg$harvest_init_multiplier),
    harvest_ids = cfg$harvest_param_ids
  )
  bounds <- list(
    lower = param_bundle$optimizer$lower,
    upper = param_bundle$optimizer$upper
  )
  default_par_t <- param_bundle$optimizer$init
  init_params_tsv <- if (!is.null(argv$init_params_tsv)) argv$init_params_tsv else NULL
  warm_start_t <- if (!is.null(init_params_tsv)) {
    read_init_params_t(init_params_tsv, bounds = bounds, cfg = cfg)
  } else {
    NULL
  }
  message(
    "MAP likelihood objective enabled: burden=lognormal NLL (per-tumor mean), ",
    "endpoint=continuous single-cell mixture NLL (2N/4N group-balanced mean, 0.5/0.5), ",
    "practical weighting default (equal tumor weighting); sigma_burden is estimated (init=", signif(cfg$sigma_burden, 6), ")",
    ", sigma_ploidy=", signif(cfg$sigma_ploidy, 6),
    ", start_with=", cfg$start_with
  )
  if (isTRUE(cfg$use_soft_prior) && cfg$lambda_prior > 0) {
    if (isTRUE(cfg$glucose_dynamic)) {
      message(
        "Soft prior enabled: lambda_prior=", signif(cfg$lambda_prior, 6),
        "; centers(log10_kappa_O, log10_o2_S0, log10_eta_o2, log10_gamma_loss, log10_rho_2N, log10_mu_hp, gamma_mu, log10_k_clear)=(",
        signif(cfg$prior_center_log10_kappa_O, 6), ", ",
        signif(cfg$prior_center_log10_o2_S0, 6), ", ",
        signif(cfg$prior_center_log10_eta_o2, 6), ", ",
        signif(cfg$prior_center_log10_gamma_loss, 6), ", ",
        signif(cfg$prior_center_log10_rho_2N, 6), ", ",
        signif(cfg$prior_center_log10_mu_hp, 6), ", ",
        signif(cfg$prior_center_gamma_mu, 6), ", ",
        signif(cfg$prior_center_log10_k_clear, 6), ")"
      )
    } else {
      message(
        "Soft prior enabled: lambda_prior=", signif(cfg$lambda_prior, 6),
        "; centers(log10_k_o, log10_kappa_O, log10_o2_S0, log10_eta_o2, log10_gamma_loss, log10_rho_2N, log10_mu_hp, gamma_mu, log10_k_clear)=(",
        signif(cfg$prior_center_log10_k_o, 6), ", ",
        signif(cfg$prior_center_log10_kappa_O, 6), ", ",
        signif(cfg$prior_center_log10_o2_S0, 6), ", ",
        signif(cfg$prior_center_log10_eta_o2, 6), ", ",
        signif(cfg$prior_center_log10_gamma_loss, 6), ", ",
        signif(cfg$prior_center_log10_rho_2N, 6), ", ",
        signif(cfg$prior_center_log10_mu_hp, 6), ", ",
        signif(cfg$prior_center_gamma_mu, 6), ", ",
        signif(cfg$prior_center_log10_k_clear, 6), ")"
      )
    }
    if (isTRUE(cfg$harvest_init_multiplier) && length(cfg$harvest_param_ids) > 0L) {
      message(
        "Harvest-specific initial-size multipliers enabled: n=",
        length(cfg$harvest_param_ids),
        ", prior N(",
        signif(cfg$prior_center_log_init_mult, 6), ", ",
        signif(cfg$prior_sd_log_init_mult, 6),
        "^2), bounds=[",
        signif(cfg$log_init_mult_lower, 6), ", ",
        signif(cfg$log_init_mult_upper, 6), "]."
      )
    }
  } else {
    message("Soft prior disabled.")
  }
  if (isTRUE(cfg$use_deoptim)) {
    message(
      "DEoptim init: mode=", cfg$de_init_mode,
      ", uniform_frac=", signif(cfg$de_init_uniform_frac, 6),
      ", sigma_frac=", signif(cfg$de_init_sigma_frac, 6),
      ", reltol=", signif(cfg$de_reltol, 6),
      ", steptol=", cfg$de_steptol,
      " (no warm-start => full uniform population)"
    )
  }
  message(
    "Growth module: ",
    if (isTRUE(cfg$O2_growth)) {
      "oxygen/ploidy growth penalty enabled"
    } else {
      "O2_growth=FALSE, using lambda_eff=lambda_base"
    },
    "; death mode=", cfg$ploidy_O2_death, "."
  )
  message(
    "Burden observation model enabled: log-normal likelihood on V(mm^3), ",
    "V_pred = sum_n n_n * (1/rho_2N), volume independent of chromosome number/ploidy, ",
    "rho_2N_range=[", signif(cfg$rho_2N_min, 6), ", ", signif(cfg$rho_2N_max, 6), "] cells/mm^3",
    "; burden_exclude_day0=", if (isTRUE(cfg$burden_exclude_day0)) "TRUE" else "FALSE"
  )
  message(
    "O2 mode: ",
    if (isTRUE(cfg$o2_burden_feedback)) "dynamic feedback" else "fixed",
    ", o2_S0_upper_bound=", signif(cfg$o2_S0_upper_bound, 6),
    ", o2_Nref=", signif(cfg$o2_Nref, 6),
    ", tau_O2_mode=", if (isTRUE(cfg$fit_tau_O2)) "fit" else "fixed",
    ", tau_O2=", if (isTRUE(cfg$fit_tau_O2)) {
      paste0("init=", signif(cfg$tau_O2_init, 6), ",range=[", signif(cfg$tau_O2_min, 6), ",", signif(cfg$tau_O2_max, 6), "]")
    } else {
      signif(cfg$tau_O2, 6)
    },
    ", o2_S0_init=", signif(cfg$o2_S0_init, 6),
    ", kappa_O_init=", signif(cfg$kappa_O_init, 6),
    ", eta_o2_init=", signif(cfg$eta_o2_init, 6),
    ", O2_crit_init=", signif(cfg$o2_crit_init, 6),
    ", n_O_init=", signif(cfg$n_O_init, 6),
    if (isTRUE(cfg$o2_burden_feedback)) {
      if (isTRUE(cfg$fit_tau_O2)) {
        "; fitted O2 params: o2_S0, kappa_O, eta_o2, O2_crit, tau_O2"
      } else {
        "; fitted O2 params: o2_S0, kappa_O, eta_o2, O2_crit"
      }
    } else {
      "; O2 fixed path still uses fitted O2_crit for hypoxia response."
    },
    "; death params: mu_hp init=", signif(cfg$mu_hp_init, 6),
    ", gamma_mu init=", signif(cfg$gamma_mu_init, 6),
    ", ploidy_O2_death=", as.character(cfg$ploidy_O2_death),
    ", k_clear init=", signif(cfg$k_clear_init, 6)
  )
  message(
    "O2 cache: bin_pct=", signif(cfg$o2_cache_bin_pct, 6),
    ", hysteresis_pct=", signif(cfg$o2_cache_hysteresis_pct, 6),
    ", profile=", if (isTRUE(cfg$o2_cache_profile)) "TRUE" else "FALSE"
  )
  set.seed(cfg$seed)
  initial_par_t <- if (!is.null(warm_start_t)) warm_start_t else default_par_t
  pass_cfg <- cfg
  pass_cfg$de_init_mode <- "uniform"
  checkpoint_writer <- function(best_par_t, best_val, iter_completed, iter_target, interrupted, stage_label) {
    best_par_t <- as.numeric(best_par_t)
    names(best_par_t) <- full_names
    status_df <- data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      stage = as.character(stage_label),
      bestval = as.numeric(best_val),
      iter_completed = as.integer(iter_completed),
      iter_target = as.integer(iter_target),
      interrupted = isTRUE(interrupted),
      stringsAsFactors = FALSE
    )
    write.table(
      status_df,
      file = file.path(checkpoint_dir, "optimizer_status_latest.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(
      data.frame(parameter = full_names, value = as.numeric(best_par_t[full_names]), row.names = NULL),
      file = file.path(checkpoint_dir, "best_params_transformed_latest.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    best_par_natural <- tryCatch(
      decode_params(best_par_t, fit_treatment = cfg$fit_treatment, fit_tau_O2 = cfg$fit_tau_O2, cfg = cfg),
      error = function(e) NULL
    )
    if (!is.null(best_par_natural)) {
      best_par_natural <- filter_family_specific_run_params_for_output_common(
        best_par_natural,
        glucose = cfg$glucose,
        glucose_dynamic = cfg$glucose_dynamic,
        misseg_loss_survival = cfg$misseg_loss_survival
      )
      best_par_natural_num <- best_par_natural[vapply(best_par_natural, is.numeric, logical(1))]
      best_par_natural_num <- best_par_natural_num[!vapply(best_par_natural_num, is.null, logical(1))]
      write.table(
        data.frame(
          parameter = names(best_par_natural_num),
          value = as.numeric(unlist(best_par_natural_num)),
          row.names = NULL
        ),
        file = file.path(checkpoint_dir, "best_params_latest.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
    saveRDS(
      list(
        timestamp = status_df$timestamp[[1]],
        stage = as.character(stage_label),
        best_par_t = best_par_t[full_names],
        best_val = as.numeric(best_val),
        iter_completed = as.integer(iter_completed),
        iter_target = as.integer(iter_target),
        interrupted = isTRUE(interrupted)
      ),
      file = file.path(checkpoint_dir, "optimizer_checkpoint_latest.rds")
    )
    invisible(NULL)
  }
  message("[single_stage] pass1 initialization: candidate[1]=initial_value, candidates[2:NP]=uniform(lower,upper).")
  objective_fn <- function(par) evaluate_objective(par, scenarios = scenarios, cfg = pass_cfg)
  single_fit <- run_optimizer(
    objective_fn = objective_fn,
    lower = bounds$lower,
    upper = bounds$upper,
    cfg = pass_cfg,
    argv = argv,
    stage_label = "single_stage",
    init_par = initial_par_t,
    checkpoint_fn = checkpoint_writer
  )
  if (isTRUE(single_fit$interrupted)) {
    message(
      "[single_stage] Optimization interrupted by user; writing best-so-far result from completed iterations: ",
      as.integer(.first_non_null_local(single_fit$iter_completed, 0L)),
      "/", as.integer(.first_non_null_local(single_fit$iter_target, cfg$itermax))
    )
  }
  best_par_t <- single_fit$best_par
  pass_comp <- evaluate_objective_components(best_par_t, scenarios = scenarios, cfg = pass_cfg)
  single_pass_log <- list(list(
    pass = 1L,
    objective = pass_comp$L,
    objective_data = pass_comp$L_data,
    objective_prior_raw = pass_comp$L_prior_raw,
    objective_prior = pass_comp$L_prior,
    objective_burden = pass_comp$L_b,
    objective_ploidy = pass_comp$L_p,
    objective_burden_neg2loglik_raw = pass_comp$objective_burden_neg2loglik_raw,
    objective_ploidy_neg2loglik_raw = pass_comp$objective_ploidy_neg2loglik_raw,
    optim = single_fit$optim_res
  ))
  optim_res <- list(
    mode = "single_stage",
    passes = single_pass_log,
    optim = single_fit$optim_res$optim
  )

  best_par <- decode_params(
    best_par_t,
    fit_treatment = cfg$fit_treatment,
    fit_tau_O2 = cfg$fit_tau_O2,
    cfg = cfg
  )
  preds <- collect_predictions(best_par, scenarios, cfg)
  final_comp <- evaluate_objective_components(best_par_t, scenarios = scenarios, cfg = cfg)
  final_obj <- final_comp$L

  best_par_num <- best_par[vapply(best_par, is.numeric, logical(1))]
  best_par_num <- filter_family_specific_run_params_for_output_common(
    best_par_num,
    glucose = cfg$glucose,
    glucose_dynamic = cfg$glucose_dynamic,
    misseg_loss_survival = cfg$misseg_loss_survival
  )
  best_par_num <- best_par_num[!vapply(best_par_num, is.null, logical(1))]
  params_df <- data.frame(
    parameter = names(best_par_num),
    value = as.numeric(unlist(best_par_num)),
    row.names = NULL
  )
  write.table(params_df, file = file.path(out_dir, "best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  stage_map <- data.frame(
    transformed_parameter = full_names,
    optimized_in = rep("single_stage", length(full_names)),
    transformed_value = as.numeric(best_par_t[full_names]),
    row.names = NULL
  )
  write.table(stage_map, file = file.path(out_dir, "fit_parameter_stages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  pass_df <- bind_rows(lapply(single_pass_log, function(x) {
    data.frame(
      pass = as.integer(x$pass),
      objective = as.numeric(x$objective),
      objective_data = as.numeric(x$objective_data),
      objective_prior_raw = as.numeric(x$objective_prior_raw),
      objective_prior = as.numeric(x$objective_prior),
      objective_burden = as.numeric(x$objective_burden),
      objective_ploidy = as.numeric(x$objective_ploidy),
      objective_burden_neg2loglik_raw = as.numeric(x$objective_burden_neg2loglik_raw),
      objective_ploidy_neg2loglik_raw = as.numeric(x$objective_ploidy_neg2loglik_raw),
      row.names = NULL
    )
  }))
  write.table(pass_df, file = file.path(out_dir, "single_stage_pass_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  n_ploidy_scenarios <- sum(vapply(scenarios, function(s) length(s$endpoint_obs_z) > 0, logical(1)))
  n_ploidy_loss_2N <- as.character(.first_non_null_local(final_comp$n_ploidy_2N, NA_integer_))
  n_ploidy_loss_4N <- as.character(.first_non_null_local(final_comp$n_ploidy_4N, NA_integer_))
  fit_mode <- if (!is.null(optim_res$mode)) as.character(optim_res$mode) else "single_stage"
  optimizer_method <- as.character(.first_non_null_local(single_fit$optim_res$method, NA_character_))
  local_refinement <- single_fit$optim_res$local_refinement
  if (is.null(local_refinement)) local_refinement <- list()
  optimizer_deoptim_objective <- as.numeric(.first_non_null_local(single_fit$optim_res$deoptim$bestval, NA_real_))
  optimizer_local_objective <- as.numeric(.first_non_null_local(local_refinement$bestval, NA_real_))
  optimizer_local_attempted <- isTRUE(.first_non_null_local(local_refinement$attempted, FALSE))
  optimizer_local_accepted <- isTRUE(.first_non_null_local(local_refinement$accepted, FALSE))
  optimizer_local_convergence <- as.integer(.first_non_null_local(local_refinement$convergence, NA_integer_))
  optimizer_local_maxit <- as.integer(.first_non_null_local(local_refinement$maxit, NA_integer_))
  deoptim_stop_iteration <- as.integer(.first_non_null_local(single_fit$iter_completed, NA_integer_))
  deoptim_iter_target <- as.integer(.first_non_null_local(single_fit$iter_target, NA_integer_))
  deoptim_stop_reason <- if (!isTRUE(cfg$use_deoptim)) {
    "not_deoptim"
  } else if (isTRUE(single_fit$interrupted)) {
    "user_interrupt"
  } else if (is.finite(deoptim_stop_iteration) && is.finite(deoptim_iter_target) && deoptim_stop_iteration < deoptim_iter_target) {
    "early_stop_reltol_or_steptol"
  } else if (is.finite(deoptim_stop_iteration) && is.finite(deoptim_iter_target) && deoptim_stop_iteration >= deoptim_iter_target) {
    "itermax_reached"
  } else {
    NA_character_
  }
  summary_df <- data.frame(
    metric = c(
      "fit_mode",
      "optimizer_method",
      "optimizer_deoptim_objective",
      "optimizer_local_objective",
      "optimizer_local_attempted",
      "optimizer_local_accepted",
      "optimizer_local_convergence",
      "optimizer_local_maxit",
      "objective",
      "objective_data",
      "objective_prior_raw",
      "objective_prior",
      "objective_burden",
      "objective_ploidy",
      "objective_burden_neg2loglik_raw",
      "objective_ploidy_neg2loglik_raw",
      "burden_exclude_day0",
      "sigma_burden",
      "sigma_ploidy",
      "use_soft_prior",
      "lambda_prior",
      "prior_center_log10_k_o",
      "prior_sd_log10_k_o",
      "prior_center_log10_kappa_O",
      "prior_sd_log10_kappa_O",
      "prior_center_log10_o2_S0",
      "prior_sd_log10_o2_S0",
      "prior_center_log10_eta_o2",
      "prior_sd_log10_eta_o2",
      "misseg_loss_survival",
      "gamma_loss_init",
      "gamma_loss_min",
      "gamma_loss_max",
      "prior_center_log10_gamma_loss",
      "prior_sd_log10_gamma_loss",
      "buffer_smax_init",
      "buffer_smax_min",
      "buffer_smax_max",
      "buffer_beta_init",
      "buffer_beta_min",
      "buffer_beta_max",
      "buffer_n_exp_init",
      "buffer_n_exp_min",
      "buffer_n_exp_max",
      "prior_center_buffer_smax",
      "prior_sd_buffer_smax",
      "prior_center_log10_buffer_beta",
      "prior_sd_log10_buffer_beta",
      "prior_center_log10_buffer_n_exp",
      "prior_sd_log10_buffer_n_exp",
      "prior_center_log10_rho_2N",
      "prior_sd_log10_rho_2N",
      "prior_center_log10_mu_hp",
      "prior_sd_log10_mu_hp",
      "prior_center_gamma_mu",
      "prior_sd_gamma_mu",
      "prior_center_log10_k_clear",
      "prior_sd_log10_k_clear",
      "n_scenarios",
      "n_ploidy_scenarios",
      "n_ploidy_loss_2N_tumors",
      "n_ploidy_loss_4N_tumors",
      "itermax",
      "NP",
      "n_cores",
      "optim_trace",
      "optim_trace_every",
      "use_deoptim",
      "deoptim_parallel",
      "de_init_mode",
      "de_init_uniform_frac",
      "de_init_sigma_frac",
      "de_reltol",
      "de_steptol",
      "optimizer_interrupted",
      "optimizer_iter_completed",
      "optimizer_iter_target",
      "deoptim_stop_iteration",
      "deoptim_iter_target",
      "deoptim_stop_reason",
      "fit_treatment",
      "o2_burden_feedback",
      "O2_growth",
      "o2_cache_bin_pct",
      "o2_cache_hysteresis_pct",
      "o2_cache_profile",
      "o2_S0_upper_bound",
      "ploidy_O2_death",
      "glucose",
      "harvest_init_multiplier",
      "n_harvest_init_params",
      "glucose_dynamic",
      "glucose_stress_mode",
      "glucose_ref_mM",
      "prior_center_log_init_mult",
      "prior_sd_log_init_mult",
      "log_init_mult_lower",
      "log_init_mult_upper",
      "start_with",
      "o2_Nref",
      "o2_S0_init",
      "o2_S0_min",
      "o2_S0_max",
      "kappa_O_init",
      "kappa_O_min",
      "kappa_O_max",
      "eta_o2_init",
      "eta_o2_min",
      "eta_o2_max",
      "G_S0_init",
      "G_S0_min",
      "G_S0_max",
      "kappa_G_init",
      "kappa_G_min",
      "kappa_G_max",
      "eta_G_init",
      "eta_G_min",
      "eta_G_max",
      "G_c_init",
      "G_c_min",
      "G_c_max",
      "tau_G",
      "tau_G_init",
      "tau_G_min",
      "tau_G_max",
      "o2_crit_init",
      "o2_crit_min",
      "o2_crit_max",
      "n_O_init",
      "n_O_min",
      "n_O_max",
      "alpha_o2_init",
      "alpha_o2_min",
      "alpha_o2_max",
      "gamma_growth_init",
      "gamma_growth_min",
      "gamma_growth_max",
      "mu_hp_init",
      "mu_hp_min",
      "mu_hp_max",
      "gamma_mu_init",
      "gamma_mu_min",
      "gamma_mu_max",
      "p_wgd_init",
      "p_wgd_min",
      "p_wgd_max",
      "p_wgd_max_init",
      "p_wgd_max_min",
      "p_wgd_max_max",
      "O2_wgd_init",
      "O2_wgd_min",
      "O2_wgd_max",
      "k_clear_init",
      "k_clear_min",
      "k_clear_max",
      "fit_tau_O2",
      "tau_O2",
      "tau_O2_init",
      "tau_O2_min",
      "tau_O2_max",
      "final_cache_g_build",
      "final_cache_g_hit",
      "final_cache_g_hysteresis",
      "paired_only",
      "init_params_tsv",
      "append_timestamp_out_dir",
      "timestamp_format",
      "dose_zero_only",
      "truncate_at_treatment",
      "ploidy_at_harvest"
    ),
    value = c(
      fit_mode,
      optimizer_method,
      as.character(optimizer_deoptim_objective),
      as.character(optimizer_local_objective),
      as.character(optimizer_local_attempted),
      as.character(optimizer_local_accepted),
      as.character(optimizer_local_convergence),
      as.character(optimizer_local_maxit),
      as.character(final_obj),
      as.character(final_comp$L_data),
      as.character(final_comp$L_prior_raw),
      as.character(final_comp$L_prior),
      as.character(final_comp$L_b),
      as.character(final_comp$L_p),
      as.character(final_comp$objective_burden_neg2loglik_raw),
      as.character(final_comp$objective_ploidy_neg2loglik_raw),
      as.character(cfg$burden_exclude_day0),
      as.character(.first_non_null_local(best_par[["sigma_burden"]], cfg$sigma_burden)),
      as.character(cfg$sigma_ploidy),
      as.character(cfg$use_soft_prior),
      as.character(cfg$lambda_prior),
      as.character(cfg$prior_center_log10_k_o),
      as.character(cfg$prior_sd_log10_k_o),
      as.character(cfg$prior_center_log10_kappa_O),
      as.character(cfg$prior_sd_log10_kappa_O),
      as.character(cfg$prior_center_log10_o2_S0),
      as.character(cfg$prior_sd_log10_o2_S0),
      as.character(cfg$prior_center_log10_eta_o2),
      as.character(cfg$prior_sd_log10_eta_o2),
      as.character(cfg$misseg_loss_survival),
      as.character(cfg$gamma_loss_init),
      as.character(cfg$gamma_loss_min),
      as.character(cfg$gamma_loss_max),
      as.character(cfg$prior_center_log10_gamma_loss),
      as.character(cfg$prior_sd_log10_gamma_loss),
      as.character(cfg$buffer_smax_init),
      as.character(cfg$buffer_smax_min),
      as.character(cfg$buffer_smax_max),
      as.character(cfg$buffer_beta_init),
      as.character(cfg$buffer_beta_min),
      as.character(cfg$buffer_beta_max),
      as.character(cfg$buffer_n_exp_init),
      as.character(cfg$buffer_n_exp_min),
      as.character(cfg$buffer_n_exp_max),
      as.character(cfg$prior_center_buffer_smax),
      as.character(cfg$prior_sd_buffer_smax),
      as.character(cfg$prior_center_log10_buffer_beta),
      as.character(cfg$prior_sd_log10_buffer_beta),
      as.character(cfg$prior_center_log10_buffer_n_exp),
      as.character(cfg$prior_sd_log10_buffer_n_exp),
      as.character(cfg$prior_center_log10_rho_2N),
      as.character(cfg$prior_sd_log10_rho_2N),
      as.character(cfg$prior_center_log10_mu_hp),
      as.character(cfg$prior_sd_log10_mu_hp),
      as.character(cfg$prior_center_gamma_mu),
      as.character(cfg$prior_sd_gamma_mu),
      as.character(cfg$prior_center_log10_k_clear),
      as.character(cfg$prior_sd_log10_k_clear),
      as.character(length(scenarios)),
      as.character(n_ploidy_scenarios),
      n_ploidy_loss_2N,
      n_ploidy_loss_4N,
      as.character(cfg$itermax),
      as.character(cfg$NP),
      as.character(cfg$n_cores),
      as.character(cfg$optim_trace),
      as.character(cfg$optim_trace_every),
      as.character(cfg$use_deoptim),
      as.character(cfg$deoptim_parallel),
      as.character(cfg$de_init_mode),
      as.character(cfg$de_init_uniform_frac),
      as.character(cfg$de_init_sigma_frac),
      as.character(cfg$de_reltol),
      as.character(cfg$de_steptol),
      as.character(isTRUE(single_fit$interrupted)),
      as.character(.first_non_null_local(single_fit$iter_completed, NA_integer_)),
      as.character(.first_non_null_local(single_fit$iter_target, NA_integer_)),
      as.character(deoptim_stop_iteration),
      as.character(deoptim_iter_target),
      as.character(deoptim_stop_reason),
      as.character(cfg$fit_treatment),
      as.character(cfg$o2_burden_feedback),
      as.character(cfg$O2_growth),
      as.character(cfg$o2_cache_bin_pct),
      as.character(cfg$o2_cache_hysteresis_pct),
      as.character(cfg$o2_cache_profile),
      as.character(cfg$o2_S0_upper_bound),
      as.character(cfg$ploidy_O2_death),
      as.character(cfg$glucose),
      as.character(cfg$harvest_init_multiplier),
      as.character(length(cfg$harvest_param_ids)),
      as.character(cfg$glucose_dynamic),
      as.character(cfg$glucose_stress_mode),
      as.character(cfg$glucose_ref_mM),
      as.character(cfg$prior_center_log_init_mult),
      as.character(cfg$prior_sd_log_init_mult),
      as.character(cfg$log_init_mult_lower),
      as.character(cfg$log_init_mult_upper),
      as.character(cfg$start_with),
      as.character(cfg$o2_Nref),
      as.character(cfg$o2_S0_init),
      as.character(cfg$o2_S0_min),
      as.character(cfg$o2_S0_max),
      as.character(cfg$kappa_O_init),
      as.character(cfg$kappa_O_min),
      as.character(cfg$kappa_O_max),
      as.character(cfg$eta_o2_init),
      as.character(cfg$eta_o2_min),
      as.character(cfg$eta_o2_max),
      as.character(cfg$G_S0_init),
      as.character(cfg$G_S0_min),
      as.character(cfg$G_S0_max),
      as.character(cfg$kappa_G_init),
      as.character(cfg$kappa_G_min),
      as.character(cfg$kappa_G_max),
      as.character(cfg$eta_G_init),
      as.character(cfg$eta_G_min),
      as.character(cfg$eta_G_max),
      as.character(cfg$G_c_init),
      as.character(cfg$G_c_min),
      as.character(cfg$G_c_max),
      as.character(cfg$tau_G),
      as.character(cfg$tau_G_init),
      as.character(cfg$tau_G_min),
      as.character(cfg$tau_G_max),
      as.character(cfg$o2_crit_init),
      as.character(cfg$o2_crit_min),
      as.character(cfg$o2_crit_max),
      as.character(cfg$n_O_init),
      as.character(cfg$n_O_min),
      as.character(cfg$n_O_max),
      as.character(cfg$alpha_o2_init),
      as.character(cfg$alpha_o2_min),
      as.character(cfg$alpha_o2_max),
      as.character(cfg$gamma_growth_init),
      as.character(cfg$gamma_growth_min),
      as.character(cfg$gamma_growth_max),
      as.character(cfg$mu_hp_init),
      as.character(cfg$mu_hp_min),
      as.character(cfg$mu_hp_max),
      as.character(cfg$gamma_mu_init),
      as.character(cfg$gamma_mu_min),
      as.character(cfg$gamma_mu_max),
      as.character(cfg$p_wgd_init),
      as.character(cfg$p_wgd_min),
      as.character(cfg$p_wgd_max),
      as.character(cfg$p_wgd_max_init),
      as.character(cfg$p_wgd_max_min),
      as.character(cfg$p_wgd_max_max),
      as.character(cfg$O2_wgd_init),
      as.character(cfg$O2_wgd_min),
      as.character(cfg$O2_wgd_max),
      as.character(cfg$k_clear_init),
      as.character(cfg$k_clear_min),
      as.character(cfg$k_clear_max),
      as.character(cfg$fit_tau_O2),
      as.character(cfg$tau_O2),
      as.character(cfg$tau_O2_init),
      as.character(cfg$tau_O2_min),
      as.character(cfg$tau_O2_max),
      as.character(.first_non_null_local(final_comp$cache_g_build, NA_integer_)),
      as.character(.first_non_null_local(final_comp$cache_g_hit, NA_integer_)),
      as.character(.first_non_null_local(final_comp$cache_g_hysteresis, NA_integer_)),
      as.character(cfg$paired_only),
      as.character(if (is.null(init_params_tsv)) NA_character_ else normalizePath(init_params_tsv, mustWork = FALSE)),
      as.character(append_ts_out_dir),
      as.character(ts_format),
      as.character(cfg$dose_zero_only),
      as.character(cfg$truncate_at_treatment),
      as.character(cfg$ploidy_at_harvest)
    )
  )
  summary_df <- filter_fit_summary_metrics_for_output_common(
    summary_df,
    glucose = cfg$glucose,
    glucose_dynamic = cfg$glucose_dynamic,
    misseg_loss_survival = cfg$misseg_loss_survival
  )
  write.table(summary_df, file = file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  write.table(preds$burden, file = file.path(out_dir, "burden_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$ploidy, file = file.path(out_dir, "terminal_ploidy_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(optim_res, file = file.path(out_dir, "deoptim_result.rds"))
  saveRDS(sanitize_cfg_for_persistence(cfg), file = file.path(out_dir, "fit_config.rds"))

  message("Done. Results written to: ", normalizePath(out_dir))
  message("Best objective: ", signif(final_obj, 6))
}

.runner_default_config_path <- function(script_dir = get_script_dir()) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  normalizePath(file.path(workflow_root, "..", "..", "config", "O2G_supply_demand.yaml"), mustWork = FALSE)
}

.runner_default_parameter_table_path <- function(script_dir = get_script_dir(), glucose = TRUE, must_exist = FALSE) {
  default_o2g_parameter_table_path_common(
    script_dir = script_dir,
    glucose = glucose,
    must_exist = must_exist
  )
}

.runner_cli_string <- function(x) {
  if (is.null(x) || !length(x)) return(NULL)
  if (is.list(x) && !is.object(x)) {
    x <- unlist(x, recursive = TRUE, use.names = FALSE)
  }
  if (!length(x)) return(NULL)
  if (length(x) == 1L && is.atomic(x) && is.na(x[[1]])) return(NULL)
  if (is.logical(x)) {
    return(if (isTRUE(x[[1]])) "TRUE" else "FALSE")
  }
  paste(as.character(x), collapse = ",")
}

.runner_resolve_path <- function(path_value, base_dir) {
  p <- .runner_cli_string(path_value)
  if (is.null(p)) return(NULL)
  p <- trimws(p)
  if (!nzchar(p)) return(p)
  if (startsWith(p, "~")) {
    return(normalizePath(path.expand(p), mustWork = FALSE))
  }
  if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) {
    return(normalizePath(p, mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, p), mustWork = FALSE)
}

.runner_read_yaml_config <- function(config_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for --mode=run but is not installed.")
  }
  cfg <- yaml::read_yaml(config_path)
  if (is.null(cfg)) cfg <- list()
  if (!is.list(cfg) || is.null(names(cfg))) {
    stop("Config must be a named YAML mapping: ", config_path)
  }
  cfg
}

.runner_resolve_config <- function(argv, script_dir = get_script_dir(), caller_wd = getwd()) {
  default_config <- .runner_default_config_path(script_dir)
  config_path <- .runner_resolve_path(argv$config, caller_wd)
  if (is.null(config_path) || !nzchar(trimws(config_path))) {
    config_path <- default_config
  }
  if (!file.exists(config_path)) {
    stop("Config not found: ", config_path)
  }
  config_dir <- dirname(config_path)
  yaml_cfg <- .runner_read_yaml_config(config_path)

  cli_cfg <- argv
  cli_cfg$config <- NULL
  cli_cfg$mode <- NULL

  if (is.null(cli_cfg$dt) && !is.null(cli_cfg$DT)) cli_cfg$dt <- cli_cfg$DT
  if (!is.null(cli_cfg$parameters) && nzchar(trimws(as.character(cli_cfg$parameters)))) {
    cli_cfg$parameter_table <- cli_cfg$parameters
  }

  cfg <- yaml_cfg
  for (nm in names(cli_cfg)) {
    cfg[[nm]] <- cli_cfg[[nm]]
  }
  if (is.null(cfg$dt) && !is.null(cfg$DT)) cfg$dt <- cfg$DT

  path_keys <- c(
    "out_root", "data_dir", "seeds_file", "parameter_table", "parameters",
    "init_params_tsv", "joint_invivo_best_dir", "joint_invitro_best_dir",
    "joint_init_candidates_tsv"
  )
  for (key in path_keys) {
    if (is.null(cfg[[key]])) next
    base_dir <- if (!is.null(cli_cfg[[key]])) caller_wd else config_dir
    cfg[[key]] <- .runner_resolve_path(cfg[[key]], base_dir)
  }

  if (!is.null(cfg$parameters) && nzchar(trimws(as.character(cfg$parameters)))) {
    cfg$parameter_table <- cfg$parameters
  }
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null_local(cfg$glucose, TRUE),
    default = TRUE
  ))
  if (is.null(cfg$parameter_table) || !nzchar(trimws(as.character(cfg$parameter_table)))) {
    cfg$parameter_table <- .runner_default_parameter_table_path(
      script_dir = script_dir,
      glucose = glucose_use,
      must_exist = TRUE
    )
  }

  list(
    script_dir = script_dir,
    caller_wd = caller_wd,
    config_path = normalizePath(config_path, mustWork = FALSE),
    config_dir = normalizePath(config_dir, mustWork = FALSE),
    yaml_cfg = yaml_cfg,
    cli_cfg = cli_cfg,
    cfg = cfg
  )
}

.runner_require_nonempty <- function(name, value) {
  txt <- .runner_cli_string(value)
  if (is.null(txt) || !nzchar(trimws(txt))) {
    stop("Required value '", name, "' is empty. Set it in YAML or CLI.")
  }
  invisible(txt)
}

.runner_parse_seed_values <- function(raw) {
  txt <- .runner_cli_string(raw)
  if (is.null(txt) || !nzchar(trimws(txt))) return(integer(0))
  parts <- trimws(unlist(strsplit(txt, "[,]", fixed = FALSE)))
  parts <- parts[nzchar(parts)]
  parts <- parts[!tolower(parts) %in% c("seed", "seeds")]
  nums <- suppressWarnings(as.integer(parts))
  nums[is.finite(nums)]
}

.runner_read_seeds_from_file <- function(path) {
  p <- .runner_cli_string(path)
  if (is.null(p) || !nzchar(trimws(p)) || !file.exists(p)) return(integer(0))
  ln <- readLines(p, warn = FALSE)
  ln <- gsub("\r", "", ln, fixed = TRUE)
  ln <- ln[nzchar(trimws(ln))]
  if (length(ln) == 0L) return(integer(0))
  .runner_parse_seed_values(paste(ln, collapse = ","))
}

.runner_resolve_seeds <- function(parsed_cfg) {
  cfg <- parsed_cfg$cfg
  cli_cfg <- parsed_cfg$cli_cfg

  cli_has_seeds_csv <- !is.null(cli_cfg$seeds_csv)
  cli_has_seeds_file <- !is.null(cli_cfg$seeds_file)

  seeds_from_file <- .runner_read_seeds_from_file(cfg$seeds_file)
  cli_seed_values <- if (cli_has_seeds_csv) .runner_parse_seed_values(cli_cfg$seeds_csv) else integer(0)
  yaml_seed_values <- .runner_parse_seed_values(cfg$seeds_csv)

  if (cli_has_seeds_csv && length(cli_seed_values) > 0L) {
    seeds <- cli_seed_values
    source <- "arg:--seeds_csv"
  } else if (cli_has_seeds_file) {
    seeds <- seeds_from_file
    source <- paste0("arg:--seeds_file(", cfg$seeds_file, ")")
  } else if (length(seeds_from_file) > 0L) {
    seeds <- seeds_from_file
    source <- paste0("file:", cfg$seeds_file)
  } else {
    seeds <- yaml_seed_values
    source <- "yaml:seeds_csv"
  }

  if (length(seeds) == 0L) {
    stop("No seeds found. Provide seeds_file or --seeds_csv.")
  }

  list(
    seeds = as.integer(seeds),
    seeds_csv = paste(as.integer(seeds), collapse = ","),
    seed_source = source
  )
}

.runner_build_fit_base_args <- function(cfg) {
  run_only_keys <- c(
    "config", "mode",
    "run_prefix", "out_root", "data_dir",
    "seeds_file", "seeds_csv",
    "append_run_prefix_timestamp", "run_prefix_timestamp_format",
    "auto_viz", "viz_report_dt", "viz_top_n",
    "seed", "out_dir",
    "append_timestamp_out_dir", "timestamp_format",
    "parameters", "DT"
  )
  fit_cfg <- cfg
  for (key in run_only_keys) fit_cfg[[key]] <- NULL

  args <- character(0)
  for (key in names(fit_cfg)) {
    val <- .runner_cli_string(fit_cfg[[key]])
    if (is.null(val) || !nzchar(trimws(val))) next
    args <- c(args, paste0("--", key, "=", val))
  }

  list(args = args, values = fit_cfg)
}

.runner_write_config_snapshots <- function(run_dir, config_path, resolved_cfg, fit_arg_values, seed_plan, fit_script, viz_script, report_script) {
  config_input_snapshot <- file.path(run_dir, "config.input.yaml")
  config_resolved_snapshot <- file.path(run_dir, "config.resolved.yaml")

  ok <- file.copy(config_path, config_input_snapshot, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to write config input snapshot: ", config_input_snapshot)
  }

  resolved_snapshot <- list(
    config_source = normalizePath(config_path, mustWork = FALSE),
    run_prefix = .runner_cli_string(resolved_cfg$run_prefix),
    out_root = .runner_cli_string(resolved_cfg$out_root),
    data_dir = .runner_cli_string(resolved_cfg$data_dir),
    seeds_file = .runner_cli_string(resolved_cfg$seeds_file),
    seeds_csv = .runner_cli_string(resolved_cfg$seeds_csv),
    seeds_use = seed_plan$seeds_csv,
    seed_source = seed_plan$seed_source,
    append_run_prefix_timestamp = as_bool(resolved_cfg$append_run_prefix_timestamp, FALSE),
    run_prefix_timestamp_format = .runner_cli_string(o2sd_first_non_null(resolved_cfg$run_prefix_timestamp_format, "%Y%m%d_%H%M%S")),
    auto_viz = as_bool(resolved_cfg$auto_viz, TRUE),
    viz_report_dt = as_num(resolved_cfg$viz_report_dt, 1),
    viz_top_n = as_int(resolved_cfg$viz_top_n, 6L),
    run_dir = normalizePath(run_dir, mustWork = FALSE),
    fit_script = normalizePath(fit_script, mustWork = FALSE),
    viz_script = normalizePath(viz_script, mustWork = FALSE),
    report_script = normalizePath(report_script, mustWork = FALSE),
    fit_args = fit_arg_values
  )
  yaml::write_yaml(resolved_snapshot, config_resolved_snapshot)

  list(
    input = config_input_snapshot,
    resolved = config_resolved_snapshot
  )
}

.runner_copy_parameter_table_snapshot <- function(parameter_table_path, dest_dir) {
  p <- .runner_cli_string(parameter_table_path)
  if (is.null(p) || !nzchar(trimws(p)) || !file.exists(p)) return(FALSE)
  file.copy(p, file.path(dest_dir, "parameter_table_input.csv"), overwrite = TRUE)
}

.runner_exec_to_log <- function(command, args, log_path, run_log_path = NULL) {
  shell_quote <- function(x) shQuote(as.character(x), type = "sh")
  cmd_txt <- paste(c(shell_quote(command), vapply(args, shell_quote, character(1))), collapse = " ")
  pipeline <- paste0("set -o pipefail; ", cmd_txt, " 2>&1 | tee ", shell_quote(log_path))
  if (!is.null(run_log_path) && nzchar(trimws(as.character(run_log_path)))) {
    pipeline <- paste0(pipeline, " | tee -a ", shell_quote(run_log_path))
  }
  shell_cmd <- paste("/bin/bash", "-lc", shQuote(pipeline, type = "sh"))
  status <- system(shell_cmd, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE)
  if (is.null(status)) 0L else as.integer(status)
}

.runner_stop_with_log_tail <- function(label, log_path, status) {
  tail_lines <- tryCatch(utils::tail(readLines(log_path, warn = FALSE), 20L), error = function(e) character(0))
  detail <- if (length(tail_lines) > 0L) {
    paste0("\nLast log lines:\n", paste(tail_lines, collapse = "\n"))
  } else {
    ""
  }
  stop(label, " failed with exit status ", status, ". See ", log_path, detail)
}

main_run_from_config <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  script_dir <- get_script_dir()
  parsed <- .runner_resolve_config(argv = argv, script_dir = script_dir, caller_wd = getwd())
  cfg <- parsed$cfg

  ignored_keys <- intersect(names(parsed$cli_cfg), c("seed", "out_dir", "append_timestamp_out_dir", "timestamp_format"))
  if (length(ignored_keys) > 0L) {
    warning(
      "--mode=run ignores CLI keys managed by the runner: ",
      paste(sort(unique(ignored_keys)), collapse = ", "),
      call. = FALSE
    )
  }

  .runner_require_nonempty("run_prefix", cfg$run_prefix)
  .runner_require_nonempty("out_root", cfg$out_root)
  .runner_require_nonempty("data_dir", cfg$data_dir)

  seed_plan <- .runner_resolve_seeds(parsed)
  append_ts <- as_bool(cfg$append_run_prefix_timestamp, FALSE)
  ts_format <- .runner_cli_string(o2sd_first_non_null(cfg$run_prefix_timestamp_format, "%Y%m%d_%H%M%S"))
  run_prefix <- .runner_cli_string(cfg$run_prefix)
  if (append_ts) {
    run_prefix <- paste0(run_prefix, "_", format(Sys.time(), ts_format))
  }

  out_root <- normalizePath(.runner_cli_string(cfg$out_root), mustWork = FALSE)
  data_dir <- normalizePath(.runner_cli_string(cfg$data_dir), mustWork = FALSE)
  auto_viz <- as_bool(cfg$auto_viz, TRUE)
  viz_report_dt <- as_num(cfg$viz_report_dt, 1)
  viz_top_n <- as_int(cfg$viz_top_n, 6L)

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  run_dir <- file.path(out_root, run_prefix)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  run_log <- file.path(run_dir, "run_status.log")
  if (!file.exists(run_log)) file.create(run_log)
  log_line <- function(...) {
    line <- paste0("[", format(Sys.time(), "%F %T"), "] ", paste0(..., collapse = ""))
    cat(line, "\n", sep = "")
    cat(line, "\n", sep = "", file = run_log, append = TRUE)
  }

  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  fit_script <- normalizePath(file.path(script_dir, "fit_model_O2G_supply_demand_MAP.R"), mustWork = FALSE)
  viz_script <- normalizePath(file.path(workflow_root, "vis", "viz_invivo_model_O2G_supply_demand_MAP_results.R"), mustWork = FALSE)
  report_script <- normalizePath(file.path(workflow_root, "report", "render_fit_report.R"), mustWork = FALSE)
  fit_base <- .runner_build_fit_base_args(cfg)
  snapshots <- .runner_write_config_snapshots(
    run_dir = run_dir,
    config_path = parsed$config_path,
    resolved_cfg = cfg,
    fit_arg_values = fit_base$values,
    seed_plan = seed_plan,
    fit_script = fit_script,
    viz_script = viz_script,
    report_script = report_script
  )

  has_parameter_snapshot <- isTRUE(.runner_copy_parameter_table_snapshot(cfg$parameter_table, run_dir))

  log_line("Running O2G_supply_demand_MAP")
  log_line("Config: ", parsed$config_path)
  log_line("Config input snapshot: ", snapshots$input)
  log_line("Config resolved snapshot: ", snapshots$resolved)
  log_line("Fit script: ", fit_script)
  log_line("Viz script: ", viz_script)
  log_line("Report script: ", report_script)
  log_line("Data dir: ", data_dir)
  log_line("Seeds: ", seed_plan$seeds_csv, " (", seed_plan$seed_source, ")")
  log_line("Run dir: ", run_dir)
  log_line("Run log: ", run_log)
  if (has_parameter_snapshot) {
    log_line("Parameter table input snapshot: ", file.path(run_dir, "parameter_table_input.csv"))
  } else {
    log_line("Parameter table snapshot: (missing --parameter_table or file not found)")
  }
  log_line("Run prefix timestamp suffix: ", append_ts, " (format=", ts_format, ")")
  log_line("Auto viz/report: ", auto_viz, " (report_dt=", viz_report_dt, ", top_n=", viz_top_n, ")")

  for (seed in seed_plan$seeds) {
    if (!is.finite(seed)) next
    seed <- as.integer(seed)
    seed_dir <- file.path(run_dir, paste0("seed", seed))
    dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)
    .runner_copy_parameter_table_snapshot(cfg$parameter_table, seed_dir)

    fit_log <- file.path(seed_dir, "fit_status.log")
    viz_log <- file.path(seed_dir, "viz_status.log")
    report_log <- file.path(seed_dir, "report_status.log")
    fit_args <- c(
      fit_script,
      "--fit_invivo",
      "--mode=fit_seed",
      paste0("--seed=", seed),
      paste0("--out_dir=", seed_dir),
      paste0("--data_dir=", data_dir),
      fit_base$args
    )

    log_line("seed=", seed, ": start")
    log_line("seed=", seed, ": fit_log=", fit_log)
    log_line("Fit command: Rscript ", paste(fit_args, collapse = " "))
    fit_status <- .runner_exec_to_log("Rscript", fit_args, fit_log, run_log_path = run_log)
    if (!identical(fit_status, 0L)) {
      .runner_stop_with_log_tail(paste0("seed=", seed, " fit"), fit_log, fit_status)
    }
    log_line("seed=", seed, ": done")

    if (auto_viz) {
      viz_args <- c(
        viz_script,
        paste0("--fit_dir=", seed_dir),
        paste0("--data_dir=", data_dir),
        paste0("--report_dt=", viz_report_dt),
        paste0("--top_n=", viz_top_n),
        "--n_cores=1"
      )
      log_line("seed=", seed, ": viz start")
      log_line("seed=", seed, ": viz_log=", viz_log)
      log_line("Viz command: Rscript ", paste(viz_args, collapse = " "))
      viz_status <- .runner_exec_to_log("Rscript", viz_args, viz_log, run_log_path = run_log)
      if (!identical(viz_status, 0L)) {
        .runner_stop_with_log_tail(paste0("seed=", seed, " viz"), viz_log, viz_status)
      }
      log_line("seed=", seed, ": viz done")

      report_args <- c(
        report_script,
        paste0("--fit_dir=", seed_dir)
      )
      log_line("seed=", seed, ": report start")
      log_line("seed=", seed, ": report_log=", report_log)
      log_line("Report command: Rscript ", paste(report_args, collapse = " "))
      report_status <- .runner_exec_to_log("Rscript", report_args, report_log, run_log_path = run_log)
      if (!identical(report_status, 0L)) {
        .runner_stop_with_log_tail(paste0("seed=", seed, " report"), report_log, report_status)
      }
      log_line("seed=", seed, ": report done")
    }
  }

  message("All done. Run directory: ", normalizePath(run_dir, mustWork = FALSE))
  invisible(normalizePath(run_dir, mustWork = FALSE))
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  mode_raw <- .runner_cli_string(argv$mode)
  if (is.null(mode_raw) || !nzchar(trimws(mode_raw))) {
    inferred_run_mode <- !is.null(argv$config)
    mode <- if (isTRUE(inferred_run_mode)) "run" else "fit_seed"
  } else {
    mode <- tolower(trimws(mode_raw))
  }

  if (mode %in% c("run", "runner")) {
    return(main_run_from_config(argv))
  }
  if (mode %in% c("fit_seed", "fit", "single_seed")) {
    return(main_fit_single_seed(argv))
  }

  stop("Unsupported mode: ", mode, ". Expected one of: run, fit_seed.")
}

if (sys.nframe() == 0) {
  main()
}
