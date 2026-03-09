#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readxl))

# -----------------------------------------------------------------------------
# Function: parse_args
# Purpose: Parse command-line arguments into a structured options object.
# Parameters:
#   - argv: Character vector of command-line arguments in --key=value format.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
parse_args <- function(argv) {
  out <- list()
  if (length(argv) == 0) return(out)
  for (a in argv) {
    if (!startsWith(a, "--")) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

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
# Function: as_num
# Purpose: Convert an input value to the target scalar/vector type with safe defaults.
# Parameters:
#   - x: Input value or vector to process.
#   - default: Fallback value used when the input is NULL or invalid.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.numeric(x))
}

# -----------------------------------------------------------------------------
# Function: as_int
# Purpose: Convert an input value to the target scalar/vector type with safe defaults.
# Parameters:
#   - x: Input value or vector to process.
#   - default: Fallback value used when the input is NULL or invalid.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.integer(x))
}

# -----------------------------------------------------------------------------
# Function: as_bool
# Purpose: Convert an input value to the target scalar/vector type with safe defaults.
# Parameters:
#   - x: Input value or vector to process.
#   - default: Fallback value used when the input is NULL or invalid.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

# -----------------------------------------------------------------------------
# Function: clip
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
#   - lo: Function-specific input argument.
#   - hi: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

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

# -----------------------------------------------------------------------------
# Function: .first_non_null_local
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - ...: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.first_non_null_local <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
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
get_param_names <- function(fit_treatment = TRUE, fit_tau_O2 = FALSE) {
  nm <- c(
    "log10_lam_min",
    "delta_lam",
    "log10_k_o",
    "log10_p_misseg",
    "log10_k_o_mis",
    "beta_buffer",
    "log10_n_exp",
    "log10_smax",
    "log10_p_wgd",
    "log10_o2_init_pct",
    "log10_o2_rate",
    "log10_o2_shape_v",
    "log10_rho_2N",
    "beta_size",
    "log10_alpha_o2",
    "o2_ref_pct",
    "gamma_growth"
  )
  if (isTRUE(fit_tau_O2)) {
    nm <- c(nm, "log10_tau_O2")
  }
  if (isTRUE(fit_treatment)) {
    nm <- c(nm, "log10_alpha", "gamma")
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
      fit_tau_O2 = isTRUE(.first_non_null_local(cfg$fit_tau_O2, FALSE))
    )
    if (length(p_names) != length(p)) p_names <- rep("", length(p))
  }
  names(p) <- p_names

  centers <- c(
    log10_k_o = as.numeric(cfg$prior_center_log10_k_o),
    log10_o2_rate = as.numeric(cfg$prior_center_log10_o2_rate),
    log10_o2_init_pct = as.numeric(cfg$prior_center_log10_o2_init_pct),
    log10_o2_shape_v = as.numeric(cfg$prior_center_log10_o2_shape_v),
    beta_size = as.numeric(cfg$prior_center_beta_size),
    log10_n_exp = as.numeric(cfg$prior_center_log10_n_exp),
    log10_rho_2N = as.numeric(cfg$prior_center_log10_rho_2N)
  )
  sds <- c(
    log10_k_o = as.numeric(cfg$prior_sd_log10_k_o),
    log10_o2_rate = as.numeric(cfg$prior_sd_log10_o2_rate),
    log10_o2_init_pct = as.numeric(cfg$prior_sd_log10_o2_init_pct),
    log10_o2_shape_v = as.numeric(cfg$prior_sd_log10_o2_shape_v),
    beta_size = as.numeric(cfg$prior_sd_beta_size),
    log10_n_exp = as.numeric(cfg$prior_sd_log10_n_exp),
    log10_rho_2N = as.numeric(cfg$prior_sd_log10_rho_2N)
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
  if (length(terms) == 0L) {
    return(list(L_prior_raw = 0, n_terms = 0, terms = numeric(0)))
  }
  list(L_prior_raw = sum(terms), n_terms = length(terms), terms = terms)
}

# -----------------------------------------------------------------------------
# Function: default_beta_size_prior_center
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
default_beta_size_prior_center <- function() {
  log(1.5) / log(2)
}

# -----------------------------------------------------------------------------
# Function: default_rho_2N_prior_bounds
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
default_rho_2N_prior_bounds <- function(cfg = NULL) {
  lo <- as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$rho_2N_min else NULL, 3.2e4))
  hi <- as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$rho_2N_max else NULL, 5.6e4))
  if (!is.finite(lo) || lo <= 0) lo <- 3.2e4
  if (!is.finite(hi) || hi <= 0) hi <- 5.6e4
  if (lo > hi) {
    tmp <- lo
    lo <- hi
    hi <- tmp
  }
  c(rho_2N_min = lo, rho_2N_max = hi)
}

# -----------------------------------------------------------------------------
# Function: default_rho_2N_prior_center
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
default_rho_2N_prior_center <- function(cfg = NULL) {
  b <- default_rho_2N_prior_bounds(cfg)
  sqrt(b[["rho_2N_min"]] * b[["rho_2N_max"]])
}

# -----------------------------------------------------------------------------
# Function: cell_volume_mm3_by_ploidy
# Purpose: Map ploidy/state values to effective cell volume under the observation model.
# Parameters:
#   - ploidy: Function-specific input argument.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
cell_volume_mm3_by_ploidy <- function(ploidy, run_params, cfg) {
  p <- pmax(as.numeric(ploidy), 1e-8)
  rho_2N <- suppressWarnings(as.numeric(run_params$rho_2N))
  rho_2N <- if (length(rho_2N) > 0) rho_2N[[1]] else NA_real_
  if (is.na(rho_2N) || !is.finite(rho_2N) || rho_2N <= 0) rho_2N <- default_rho_2N_prior_center(cfg)
  beta_size <- as.numeric(.first_non_null_local(run_params$beta_size, default_beta_size_prior_center()))
  if (!is.finite(beta_size)) beta_size <- default_beta_size_prior_center()
  (1 / rho_2N) * (p / 2)^beta_size
}

# -----------------------------------------------------------------------------
# Function: cell_volume_mm3_by_N
# Purpose: Map ploidy/state values to effective cell volume under the observation model.
# Parameters:
#   - N: Ploidy state value or chromosome-copy count.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
cell_volume_mm3_by_N <- function(N, run_params, cfg) {
  p_weighted <- weighted_ploidy_from_total_N(
    N_total = as.numeric(N),
    chr_lengths_bp = cfg$chr_lengths_bp
  )
  cell_volume_mm3_by_ploidy(p_weighted, run_params = run_params, cfg = cfg)
}

# -----------------------------------------------------------------------------
# Function: burden_volume_mm3_from_state
# Purpose: Convert state vector to predicted burden volume in mm^3.
# Parameters:
#   - v: Function-specific input argument.
#   - grid_pre: Ploidy grid for pre-missegregation compartment.
#   - R0: Number of pre-missegregation states.
#   - R1: Number of post-missegregation states.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
#   - vol_by_N: Optional precomputed per-state cell volume lookup.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
burden_volume_mm3_from_state <- function(v, grid_pre, R0, R1, run_params, cfg, vol_by_N = NULL) {
  if (length(v) != (R0 + R1)) stop("State length mismatch in burden_volume_mm3_from_state.")
  if (is.null(vol_by_N)) {
    vol_by_N <- cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg)
  }
  counts_pre <- v[seq_len(R0)]
  counts_post <- v[R0 + seq_len(R1)]
  counts_N <- counts_pre + counts_post
  sum(as.numeric(counts_N) * as.numeric(vol_by_N), na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# Function: get_script_dir
# Purpose: Resolve script directory path for robust relative file loading.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grepl("^--file=", args)]
  if (length(farg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", farg[[1]])))
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
  names(par_transformed) <- get_param_names(
    fit_treatment = fit_treatment,
    fit_tau_O2 = fit_tau_O2
  )
  lam_min <- 10^par_transformed["log10_lam_min"]
  lam_max <- lam_min + exp(par_transformed["delta_lam"])
  tau_O2 <- as.numeric(.first_non_null_local(
    if (isTRUE(fit_tau_O2) && "log10_tau_O2" %in% names(par_transformed)) 10^par_transformed["log10_tau_O2"] else NULL,
    if (!is.null(cfg)) cfg$tau_O2 else NULL,
    if (!is.null(cfg)) cfg$tau_O2_init else NULL,
    2.0
  ))
  if (!is.finite(tau_O2) || tau_O2 <= 0) tau_O2 <- 2.0
  list(
    lam_min = lam_min,
    lam_max = lam_max,
    k_o = 10^par_transformed["log10_k_o"],
    p_misseg = 10^par_transformed["log10_p_misseg"],
    k_o_mis = 10^par_transformed["log10_k_o_mis"],
    beta_buffer = par_transformed["beta_buffer"],
    n_exp = 10^par_transformed["log10_n_exp"],
    smax = 10^par_transformed["log10_smax"],
    p_wgd = 10^par_transformed["log10_p_wgd"],
    o2_init_pct = 10^par_transformed["log10_o2_init_pct"],
    o2_rate = 10^par_transformed["log10_o2_rate"],
    o2_shape_v = 10^par_transformed["log10_o2_shape_v"],
    rho_2N = 10^par_transformed["log10_rho_2N"],
    beta_size = par_transformed["beta_size"],
    alpha_o2 = 10^par_transformed["log10_alpha_o2"],
    o2_ref_pct = par_transformed["o2_ref_pct"],
    gamma_growth = par_transformed["gamma_growth"],
    tau_O2 = tau_O2,
    c_vol_2N_eff_mm3 = 10^-par_transformed["log10_rho_2N"],
    ratio_4N_2N = 2^par_transformed["beta_size"],
    alpha = if (isTRUE(fit_treatment)) 10^par_transformed["log10_alpha"] else 0,
    gamma = if (isTRUE(fit_treatment)) par_transformed["gamma"] else 1
  )
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
  k_o_v <- need_pos(getv(c("k_o"), default = 50), "k_o")
  p_misseg_v <- need_pos(getv(c("p_misseg"), default = 1e-4), "p_misseg")
  k_o_mis_v <- need_pos(getv(c("k_o_mis"), default = 50), "k_o_mis")
  beta_buffer_v <- getv(c("beta_buffer"), default = 0.0)
  if (!is.finite(beta_buffer_v)) stop("Warm-start parameter must be finite: beta_buffer")
  n_exp_v <- need_pos(getv(c("n_exp"), default = 1.0), "n_exp")
  smax_v <- need_pos(getv(c("smax"), default = 1.0), "smax")
  p_wgd_v <- need_pos(getv(c("p_wgd"), default = 1e-6), "p_wgd")
  o2_cap_v <- as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_cap_pct else NULL, 5.0))
  if (!is.finite(o2_cap_v) || o2_cap_v <= 0 || o2_cap_v > 100) o2_cap_v <- 5.0
  o2_init_v <- need_pos(
    getv(c("o2_init_pct"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_init_pct_init else NULL, 0.5))),
    "o2_init_pct"
  )
  o2_init_v <- min(max(o2_init_v, 1e-6), o2_cap_v - 1e-6)
  o2_rate_v <- need_pos(
    getv(c("o2_rate"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_rate_init else NULL, 1.0))),
    "o2_rate"
  )
  o2_shape_v <- need_pos(
    getv(c("o2_shape_v"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_shape_v_init else NULL, 1.0))),
    "o2_shape_v"
  )
  rho_2N_v <- getv(c("rho_2N"), default = NA_real_)
  rho_2N_v <- need_pos(
    if (is.finite(rho_2N_v) && rho_2N_v > 0) rho_2N_v else default_rho_2N_prior_center(cfg),
    "rho_2N"
  )
  beta_size_v <- getv(c("beta_size"), default = default_beta_size_prior_center())
  if (!is.finite(beta_size_v)) stop("Warm-start parameter must be finite: beta_size")
  alpha_o2_v <- need_pos(getv(c("alpha_o2"), default = 0.5), "alpha_o2")
  o2_ref_pct_v <- getv(c("o2_ref_pct"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_ref_pct_init else NULL, 2.5)))
  if (!is.finite(o2_ref_pct_v) || o2_ref_pct_v < 0 || o2_ref_pct_v > 100) {
    stop("Warm-start parameter must be finite and in [0,100]: o2_ref_pct")
  }
  gamma_growth_v <- getv(c("gamma_growth"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$gamma_growth_init else NULL, 2.0)))
  if (!is.finite(gamma_growth_v) || gamma_growth_v <= 0) {
    stop("Warm-start parameter must be > 0: gamma_growth")
  }
  tau_O2_v <- need_pos(
    getv(c("tau_O2"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$tau_O2_init else NULL, 2.0))),
    "tau_O2"
  )

  out <- c(
    log10_lam_min = log10(lam_min_v),
    delta_lam = log(lam_gap_v),
    log10_k_o = log10(k_o_v),
    log10_p_misseg = log10(p_misseg_v),
    log10_k_o_mis = log10(k_o_mis_v),
    beta_buffer = beta_buffer_v,
    log10_n_exp = log10(n_exp_v),
    log10_smax = log10(smax_v),
    log10_p_wgd = log10(p_wgd_v),
    log10_o2_init_pct = log10(o2_init_v),
    log10_o2_rate = log10(o2_rate_v),
    log10_o2_shape_v = log10(o2_shape_v),
    log10_rho_2N = log10(rho_2N_v),
    beta_size = beta_size_v,
    log10_alpha_o2 = log10(alpha_o2_v),
    o2_ref_pct = o2_ref_pct_v,
    gamma_growth = gamma_growth_v
  )
  if (isTRUE(fit_tau_O2)) {
    out <- c(out, log10_tau_O2 = log10(tau_O2_v))
  }

  if (isTRUE(fit_treatment)) {
    alphav <- need_pos(getv(c("alpha")), "alpha")
    gammav <- getv(c("gamma"))
    out <- c(out, log10_alpha = log10(alphav), gamma = gammav)
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
    if ("log10_o2_shape_v" %in% missing_names) {
      vals[["log10_o2_shape_v"]] <- log10(as.numeric(.first_non_null_local(cfg$o2_shape_v_init, 1.0)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_alpha_o2" %in% missing_names) {
      vals[["log10_alpha_o2"]] <- log10(as.numeric(.first_non_null_local(cfg$alpha_o2_init, 0.5)))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("o2_ref_pct" %in% missing_names) {
      vals[["o2_ref_pct"]] <- as.numeric(.first_non_null_local(cfg$o2_ref_pct_init, 2.5))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("gamma_growth" %in% missing_names) {
      vals[["gamma_growth"]] <- as.numeric(.first_non_null_local(cfg$gamma_growth_init, 2.0))
      missing_names <- setdiff(full_names, names(vals))
    }
    if ("log10_tau_O2" %in% missing_names) {
      vals[["log10_tau_O2"]] <- log10(as.numeric(.first_non_null_local(cfg$tau_O2_init, 2.0)))
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
#   - o2_init_pct_min: Function-specific input argument.
#   - o2_init_pct_max: Function-specific input argument.
#   - o2_rate_min: Function-specific input argument.
#   - o2_rate_max: Function-specific input argument.
#   - o2_shape_v_min: Function-specific input argument.
#   - o2_shape_v_max: Function-specific input argument.
#   - tau_O2_min: Lower bound of tau_O2 when estimated.
#   - tau_O2_max: Upper bound of tau_O2 when estimated.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
make_bounds <- function(fit_treatment = TRUE,
                        fit_tau_O2 = FALSE,
                        rho_2N_min = 3.2e4, rho_2N_max = 5.6e4,
                        o2_init_pct_min = 1e-3, o2_init_pct_max = 4.9,
                        o2_rate_min = 1e-3, o2_rate_max = 1e2,
                        o2_shape_v_min = 1e-2, o2_shape_v_max = 20,
                        alpha_o2_min = 1e-2, alpha_o2_max = 10,
                        o2_ref_pct_min = 0, o2_ref_pct_max = 5,
                        gamma_growth_min = 2.0, gamma_growth_max = 2.0,
                        tau_O2_min = 1e-3, tau_O2_max = 1e3) {
  rho_2N_min <- as.numeric(rho_2N_min)
  rho_2N_max <- as.numeric(rho_2N_max)
  if (!is.finite(rho_2N_min) || rho_2N_min <= 0) rho_2N_min <- 3.2e4
  if (!is.finite(rho_2N_max) || rho_2N_max <= 0) rho_2N_max <- 5.6e4
  if (rho_2N_min > rho_2N_max) {
    tmp <- rho_2N_min
    rho_2N_min <- rho_2N_max
    rho_2N_max <- tmp
  }
  o2_init_pct_min <- as.numeric(o2_init_pct_min)
  o2_init_pct_max <- as.numeric(o2_init_pct_max)
  if (!is.finite(o2_init_pct_min) || o2_init_pct_min <= 0) o2_init_pct_min <- 1e-3
  if (!is.finite(o2_init_pct_max) || o2_init_pct_max <= 0) o2_init_pct_max <- 4.9
  if (o2_init_pct_min > o2_init_pct_max) {
    tmp <- o2_init_pct_min
    o2_init_pct_min <- o2_init_pct_max
    o2_init_pct_max <- tmp
  }
  o2_rate_min <- as.numeric(o2_rate_min)
  o2_rate_max <- as.numeric(o2_rate_max)
  if (!is.finite(o2_rate_min) || o2_rate_min <= 0) o2_rate_min <- 1e-3
  if (!is.finite(o2_rate_max) || o2_rate_max <= 0) o2_rate_max <- 1e2
  if (o2_rate_min > o2_rate_max) {
    tmp <- o2_rate_min
    o2_rate_min <- o2_rate_max
    o2_rate_max <- tmp
  }
  o2_shape_v_min <- as.numeric(o2_shape_v_min)
  o2_shape_v_max <- as.numeric(o2_shape_v_max)
  if (!is.finite(o2_shape_v_min) || o2_shape_v_min <= 0) o2_shape_v_min <- 1e-2
  if (!is.finite(o2_shape_v_max) || o2_shape_v_max <= 0) o2_shape_v_max <- 20
  if (o2_shape_v_min > o2_shape_v_max) {
    tmp <- o2_shape_v_min
    o2_shape_v_min <- o2_shape_v_max
    o2_shape_v_max <- tmp
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
  o2_ref_pct_min <- as.numeric(o2_ref_pct_min)
  o2_ref_pct_max <- as.numeric(o2_ref_pct_max)
  if (!is.finite(o2_ref_pct_min) || o2_ref_pct_min < 0) o2_ref_pct_min <- 0
  if (!is.finite(o2_ref_pct_max) || o2_ref_pct_max < 0) o2_ref_pct_max <- 5
  if (o2_ref_pct_min > o2_ref_pct_max) {
    tmp <- o2_ref_pct_min
    o2_ref_pct_min <- o2_ref_pct_max
    o2_ref_pct_max <- tmp
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
    log10_p_misseg = log10(1e-8),
    log10_k_o_mis = log10(1e-1),
    beta_buffer = 0.0,
    log10_n_exp = log10(1e-1),
    log10_smax = log10(1e-4),
    log10_p_wgd = log10(1e-8),
    log10_o2_init_pct = log10(o2_init_pct_min),
    log10_o2_rate = log10(o2_rate_min),
    log10_o2_shape_v = log10(o2_shape_v_min),
    log10_rho_2N = log10(rho_2N_min),
    beta_size = 0.2,
    log10_alpha_o2 = log10(alpha_o2_min),
    o2_ref_pct = o2_ref_pct_min,
    gamma_growth = gamma_growth_min
  )
  upper <- c(
    log10_lam_min = log10(5),
    delta_lam = log(5),
    log10_k_o = log10(1e4),
    log10_p_misseg = log10(0.08),
    log10_k_o_mis = log10(1e4),
    beta_buffer = 10.0,
    log10_n_exp = log10(5),
    log10_smax = log10(0.9),
    log10_p_wgd = log10(1e-1),
    log10_o2_init_pct = log10(o2_init_pct_max),
    log10_o2_rate = log10(o2_rate_max),
    log10_o2_shape_v = log10(o2_shape_v_max),
    log10_rho_2N = log10(rho_2N_max),
    beta_size = 1.2,
    log10_alpha_o2 = log10(alpha_o2_max),
    o2_ref_pct = o2_ref_pct_max,
    gamma_growth = gamma_growth_max
  )

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

  list(lower = lower, upper = upper)
}

# -----------------------------------------------------------------------------
# Function: read_parameter_table_transformed
# Purpose: Read transformed-parameter initialization and bounds from CSV.
# Parameters:
#   - path: Path to CSV file containing transformed parameter metadata.
#   - param_names: Ordered transformed parameter names required by optimizer.
# Returns:
#   List with named vectors: init, lower, upper.
# -----------------------------------------------------------------------------
read_parameter_table_transformed <- function(path, param_names) {
  if (!file.exists(path)) stop("Parameter table CSV not found: ", path)
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
  req_cols <- c("param_name", "init_value", "lower_bound", "upper_bound")
  miss <- setdiff(req_cols, names(tab))
  if (length(miss) > 0L) {
    stop("Parameter table missing required columns: ", paste(miss, collapse = ", "))
  }
  tab$param_name <- trimws(as.character(tab$param_name))
  tab <- tab[nzchar(tab$param_name), , drop = FALSE]
  if (anyDuplicated(tab$param_name)) {
    dup <- unique(tab$param_name[duplicated(tab$param_name)])
    stop("Duplicate param_name in parameter table: ", paste(dup, collapse = ", "))
  }
  missing_params <- setdiff(param_names, tab$param_name)
  if (length(missing_params) > 0L) {
    stop("Parameter table missing required rows: ", paste(missing_params, collapse = ", "))
  }
  tab <- tab[match(param_names, tab$param_name), , drop = FALSE]

  init <- suppressWarnings(as.numeric(tab$init_value))
  lower <- suppressWarnings(as.numeric(tab$lower_bound))
  upper <- suppressWarnings(as.numeric(tab$upper_bound))
  if (any(!is.finite(init))) {
    bad <- tab$param_name[!is.finite(init)]
    stop("Non-finite init_value in parameter table for: ", paste(bad, collapse = ", "))
  }
  if (any(!is.finite(lower))) {
    bad <- tab$param_name[!is.finite(lower)]
    stop("Non-finite lower_bound in parameter table for: ", paste(bad, collapse = ", "))
  }
  if (any(!is.finite(upper))) {
    bad <- tab$param_name[!is.finite(upper)]
    stop("Non-finite upper_bound in parameter table for: ", paste(bad, collapse = ", "))
  }
  if (any(lower > upper)) {
    bad <- tab$param_name[lower > upper]
    stop("lower_bound > upper_bound in parameter table for: ", paste(bad, collapse = ", "))
  }

  init <- pmin(pmax(init, lower), upper)
  names(init) <- param_names
  names(lower) <- param_names
  names(upper) <- param_names
  list(init = init, lower = lower, upper = upper)
}

# -----------------------------------------------------------------------------
# Function: prepare_data
# Purpose: Load raw experiment inputs and assemble per-scenario fitting data.
# Parameters:
#   - dt_path: Path to burden observation table (Excel file).
#   - ploidy_path: Path to terminal ploidy table (TSV file).
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
prepare_data <- function(dt_path, ploidy_path, cfg) {
  if (!file.exists(dt_path)) stop("Tumor-burden xlsx not found: ", dt_path)
  if (!file.exists(ploidy_path)) stop("Ploidy tsv not found: ", ploidy_path)

  dt <- readxl::read_excel(dt_path)
  required <- c("harvest", "Dose", "Day of 1st treatment")
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in tumor-burden sheet: ", paste(missing_cols, collapse = ", "))
  }

  day_cols <- grep("^Day_", names(dt), value = TRUE)
  if (length(day_cols) == 0) stop("No Day_* columns found in tumor-burden sheet.")
  day_vals <- as.numeric(sub("^Day_", "", day_cols))

  day_num_df <- lapply(day_cols, function(col) suppressWarnings(as.numeric(dt[[col]])))
  names(day_num_df) <- day_cols
  day_num_df <- as.data.frame(day_num_df, stringsAsFactors = FALSE)

  pl <- read.delim(ploidy_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("file", "ploidy") %in% names(pl))) {
    stop("Ploidy table must include columns: file, ploidy")
  }
  pl$harvest <- sub(".sps.cbs", "", pl$file, fixed = TRUE)
  pl_by_harvest <- split(pl$ploidy, pl$harvest)

  scenarios <- vector("list", nrow(dt))
  keep <- logical(nrow(dt))
  n_ploidy_scaled_by22 <- 0L
  for (i in seq_len(nrow(dt))) {
    h <- as.character(dt$harvest[[i]])
    if (!nzchar(h)) next

    cohort <- if (grepl("2N", h, fixed = TRUE)) "2N" else if (grepl("4N", h, fixed = TRUE)) "4N" else NA_character_
    if (!is.finite(as_num(dt$Dose[[i]], NA_real_))) next
    dose <- as_num(dt$Dose[[i]], NA_real_)
    if (isTRUE(cfg$dose_zero_only) && dose != 0) next
    treat_day <- as_num(dt[["Day of 1st treatment"]][[i]], Inf)
    if (!is.finite(treat_day)) treat_day <- Inf

    y <- as.numeric(day_num_df[i, ])
    idx <- which(is.finite(y))
    if (length(idx) < 2) next

    full_days <- day_vals[idx]
    full_burden <- y[idx]

    # Dataset doc says missing are trailing NAs; enforce to avoid ambiguous rows.
    if (any(diff(idx) > 1)) next

    if (isTRUE(cfg$truncate_at_treatment)) {
      keep_pre <- full_days <= treat_day
      obs_days <- full_days[keep_pre]
      obs_burden <- full_burden[keep_pre]
    } else {
      obs_days <- full_days
      obs_burden <- full_burden
    }
    if (length(obs_days) < 2) next

    obs_pl <- pl_by_harvest[[h]]
    if (is.null(obs_pl)) {
      obs_z <- numeric(0)
    } else {
      obs_raw <- as.numeric(obs_pl)
      obs_raw <- obs_raw[is.finite(obs_raw)]
      if (length(obs_raw) == 0L) {
        obs_z <- numeric(0)
      } else {
        # Use chromosome-count scale everywhere in fitting:
        # observed ploidy (2N-scale) -> chromosome-count scale.
        obs_z <- obs_raw * 22
        n_ploidy_scaled_by22 <- n_ploidy_scaled_by22 + 1L
      }
    }

    scenarios[[i]] <- list(
      harvest = h,
      cohort = cohort,
      dose = dose,
      treat_day = treat_day,
      obs_days = obs_days,
      obs_burden = obs_burden,
      sim_end_day = if (isTRUE(cfg$ploidy_at_harvest)) max(full_days) else max(obs_days),
      harvest_day = max(full_days),
      ploidy_obs_z = obs_z
    )
    keep[[i]] <- TRUE
  }

  scenarios <- scenarios[keep]
  if (length(scenarios) == 0) stop("No valid scenarios after preprocessing.")
  paired_only <- isTRUE(.first_non_null_local(cfg$paired_only, TRUE))
  n_before_pair_filter <- length(scenarios)
  n_ploidy_before_pair_filter <- sum(vapply(scenarios, function(s) length(s$ploidy_obs_z) > 0, logical(1)))

  if (paired_only) {
    scenarios <- scenarios[vapply(scenarios, function(s) length(s$ploidy_obs_z) > 0, logical(1))]
    if (length(scenarios) == 0) {
      stop("paired_only=TRUE but no scenarios have both burden and terminal ploidy data.")
    }
    if (length(scenarios) < n_before_pair_filter) {
      message(
        "paired_only=TRUE: retained ", length(scenarios), "/", n_before_pair_filter,
        " scenarios with both burden+ploidy (dropped ", n_before_pair_filter - length(scenarios), ")."
      )
    }
  }

  if (is.finite(cfg$max_scenarios) && cfg$max_scenarios > 0) {
    scenarios <- scenarios[seq_len(min(length(scenarios), as.integer(cfg$max_scenarios)))]
  }

  matched_ploidy <- sum(vapply(scenarios, function(s) length(s$ploidy_obs_z) > 0, logical(1)))
  message(
    "Prepared scenarios: ", length(scenarios),
    " (with terminal ploidy: ", matched_ploidy,
    "; paired_only=", paired_only,
    "; pre_pair_filter_ploidy=", n_ploidy_before_pair_filter, "/", n_before_pair_filter, ")"
  )
  message("Ploidy observation scaling: chromosome-count mode enabled (obs_z = raw_ploidy * 22). scaled_rows=", n_ploidy_scaled_by22)
  scenarios
}

# -----------------------------------------------------------------------------
# Function: prepare_cpp_scenarios
# Purpose: Convert R scenario objects into C++-friendly arrays and indices.
# Parameters:
#   - scenarios: List of scenario-specific observation data and metadata.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
prepare_cpp_scenarios <- function(scenarios, cfg) {
  n <- length(scenarios)
  cohort_code <- integer(n)
  dose_vec <- numeric(n)
  treat_day_vec <- numeric(n)
  obs_steps_list <- vector("list", n)
  sim_end_step_vec <- integer(n)
  obs_burden_list <- vector("list", n)
  keep_burden_list <- vector("list", n)
  ploidy_z_list <- vector("list", n)

  for (i in seq_len(n)) {
    sc <- scenarios[[i]]
    cohort_code[[i]] <- if (identical(sc$cohort, "2N")) 0L else 1L
    dose_vec[[i]] <- as.numeric(sc$dose)
    treat_day_vec[[i]] <- as.numeric(sc$treat_day)
    obs_steps_list[[i]] <- as.integer(round(as.numeric(sc$obs_days) / cfg$DT))
    sim_end_step_vec[[i]] <- as.integer(round(as.numeric(sc$sim_end_day) / cfg$DT))
    obs_burden <- as.numeric(sc$obs_burden)
    obs_burden_list[[i]] <- obs_burden

    day_obs <- as.numeric(sc$obs_days)
    keep_day <- rep(TRUE, length(obs_burden))
    if (isTRUE(.first_non_null_local(cfg$burden_exclude_day0, TRUE)) &&
        length(day_obs) == length(obs_burden)) {
      keep_day <- is.finite(day_obs) & (day_obs > 0)
    }
    keep_burden_list[[i]] <- as.logical(keep_day)

    z <- as.numeric(sc$ploidy_obs_z)
    ploidy_z_list[[i]] <- z[is.finite(z)]
  }

  list(
    cohort_code = cohort_code,
    dose = dose_vec,
    treat_day = treat_day_vec,
    obs_steps = obs_steps_list,
    sim_end_step = sim_end_step_vec,
    obs_burden = obs_burden_list,
    keep_burden = keep_burden_list,
    ploidy_z = ploidy_z_list
  )
}

# -----------------------------------------------------------------------------
# Function: build_model_core
# Purpose: Construct model core objects used repeatedly during simulation/evaluation.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - cfg: Configuration list controlling model options, bounds, and optimization settings.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
build_model_core <- function(run_params = NULL, cfg) {
  grid_pre <- cfg$N_MIN:cfg$N_MAX
  grid_post <- cfg$N_MIN:cfg$N_MAX
  R0 <- length(grid_pre)
  R1 <- length(grid_post)

  init_state_2N <- make_init_state(
    grid_pre = grid_pre,
    grid_post = grid_post,
    ploidy = 2,
    layer = "pre",
    N_UNIT = cfg$N_UNIT,
    total_size = cfg$init_total_size,
    chr_lengths_bp = cfg$chr_lengths_bp
  )
  init_state_4N <- make_init_state(
    grid_pre = grid_pre,
    grid_post = grid_post,
    ploidy = 4,
    layer = "post",
    N_UNIT = cfg$N_UNIT,
    total_size = cfg$init_total_size,
    chr_lengths_bp = cfg$chr_lengths_bp
  )

  list(
    grid_pre = grid_pre,
    grid_post = grid_post,
    R0 = R0,
    R1 = R1,
    init_state_2N = init_state_2N,
    init_state_4N = init_state_4N
  )
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
  R1 <- model_core$R1
  init_state <- if (scenario$cohort == "2N") model_core$init_state_2N else model_core$init_state_4N

  obs_steps <- as.integer(round(scenario$obs_days / cfg$DT))
  sim_end_step <- as.integer(round(scenario$sim_end_day / cfg$DT))
  vol_by_N <- cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg)

  if (!exists("cpp_o2simps_simulate_one", mode = "function", inherits = TRUE)) {
    stop("Required C++ function missing: cpp_o2simps_simulate_one")
  }

  o2_cap_use <- as.numeric(.first_non_null_local(cfg$o2_cap_pct, 5.0))
  if (!is.finite(o2_cap_use) || o2_cap_use <= 0 || o2_cap_use > 100) o2_cap_use <- 5.0
  o2_curve_type_use <- as.character(.first_non_null_local(cfg$o2_curve_type, "gompertz"))
  o2_init_use <- as.numeric(.first_non_null_local(run_params$o2_init_pct, cfg$o2_init_pct_init, 0.5))
  if (!is.finite(o2_init_use) || o2_init_use <= 0) o2_init_use <- 0.5
  o2_init_use <- min(max(o2_init_use, 1e-6), o2_cap_use - 1e-6)
  o2_rate_use <- as.numeric(.first_non_null_local(run_params$o2_rate, cfg$o2_rate_init, 1.0))
  if (!is.finite(o2_rate_use) || o2_rate_use <= 0) o2_rate_use <- 1.0
  o2_shape_v_use <- as.numeric(.first_non_null_local(run_params$o2_shape_v, cfg$o2_shape_v_init, 1.0))
  if (!is.finite(o2_shape_v_use) || o2_shape_v_use <= 0) o2_shape_v_use <- 1.0
  tau_O2_use <- as.numeric(.first_non_null_local(run_params$tau_O2, cfg$tau_O2, cfg$tau_O2_init, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  o2_anchor_use <- as.numeric(.first_non_null_local(cfg$o2_anchor_N, cfg$init_total_size, 1e6))
  if (!is.finite(o2_anchor_use) || o2_anchor_use < 0) o2_anchor_use <- 1e6
  p_wgd_use <- as.numeric(.first_non_null_local(run_params$p_wgd, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
  boundary_mode <- as.character(.first_non_null_local(run_params$boundary, "drop"))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 0)

  sim <- cpp_o2simps_simulate_one(
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
    crowding = as.character(cfg$crowding),
    K = as.numeric(cfg$K),
    min_pop = as.numeric(cfg$min_pop),
    O2_cap = as.numeric(o2_cap_use),
    o2_feedback = isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE)),
    o2_curve_type = as.character(o2_curve_type_use),
    o2_init = as.numeric(o2_init_use),
    o2_rate = as.numeric(o2_rate_use),
    o2_shape_v = as.numeric(o2_shape_v_use),
    tau_O2 = as.numeric(tau_O2_use),
    o2_anchor_N = as.numeric(o2_anchor_use),
    o2_logN_eps = as.numeric(.first_non_null_local(cfg$o2_logN_eps, 1.0)),
    o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg$o2_cache_bin_pct, 0.01)),
    o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg$o2_cache_hysteresis_pct, 0.005)),
    o2_cache_profile = isTRUE(.first_non_null_local(cfg$o2_cache_profile, FALSE)),
    lam_min = as.numeric(run_params$lam_min),
    lam_max = as.numeric(run_params$lam_max),
    k_o = as.numeric(run_params$k_o),
    has_p_misseg = !is.null(run_params$p_misseg),
    p_misseg = as.numeric(.first_non_null_local(run_params$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(run_params$k_o_mis, 50.0)),
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = 0.0,
    p_wgd = as.numeric(p_wgd_use),
    boundary = boundary_mode,
    eps_tail = as.numeric(1e-8),
    beta_buffer = as.numeric(.first_non_null_local(run_params$beta_buffer, 0.0)),
    n_exp = as.numeric(.first_non_null_local(run_params$n_exp, 1.0)),
    smax = as.numeric(.first_non_null_local(run_params$smax, 1.0)),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = as.numeric(.first_non_null_local(run_params$beta_size, cfg$prior_center_beta_size, default_beta_size_prior_center())),
    alpha_o2 = as.numeric(.first_non_null_local(run_params$alpha_o2, cfg$alpha_o2_init, 0.5)),
    o2_ref_pct = as.numeric(.first_non_null_local(run_params$o2_ref_pct, cfg$o2_ref_pct_init, 2.5)),
    gamma_growth = as.numeric(.first_non_null_local(run_params$gamma_growth, cfg$gamma_growth_init, 2.0)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(burden_floor)
  )

  frac_N <- as.numeric(sim$frac_N)
  names(frac_N) <- as.character(grid_pre)

  list(
    Ntot_obs = as.numeric(sim$Ntot_obs),
    Vmm3_obs = as.numeric(sim$Vmm3_obs),
    frac_N = frac_N
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
  vol_by_N <- cell_volume_mm3_by_N(model_core$grid_pre, run_params = rp, cfg = cfg_eval)

  o2_cap_use <- as.numeric(.first_non_null_local(cfg_eval$o2_cap_pct, 5.0))
  if (!is.finite(o2_cap_use) || o2_cap_use <= 0 || o2_cap_use > 100) o2_cap_use <- 5.0
  o2_curve_type_use <- as.character(.first_non_null_local(cfg_eval$o2_curve_type, "gompertz"))
  o2_init_use <- as.numeric(.first_non_null_local(rp$o2_init_pct, cfg_eval$o2_init_pct_init, 0.5))
  if (!is.finite(o2_init_use) || o2_init_use <= 0) o2_init_use <- 0.5
  o2_init_use <- min(max(o2_init_use, 1e-6), o2_cap_use - 1e-6)
  o2_rate_use <- as.numeric(.first_non_null_local(rp$o2_rate, cfg_eval$o2_rate_init, 1.0))
  if (!is.finite(o2_rate_use) || o2_rate_use <= 0) o2_rate_use <- 1.0
  o2_shape_v_use <- as.numeric(.first_non_null_local(rp$o2_shape_v, cfg_eval$o2_shape_v_init, 1.0))
  if (!is.finite(o2_shape_v_use) || o2_shape_v_use <= 0) o2_shape_v_use <- 1.0
  tau_O2_use <- as.numeric(.first_non_null_local(rp$tau_O2, cfg_eval$tau_O2, cfg_eval$tau_O2_init, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  o2_anchor_use <- as.numeric(.first_non_null_local(cfg_eval$o2_anchor_N, cfg_eval$init_total_size, 1e6))
  if (!is.finite(o2_anchor_use) || o2_anchor_use < 0) o2_anchor_use <- 1e6
  p_wgd_use <- as.numeric(.first_non_null_local(rp$p_wgd, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
  boundary_mode <- as.character(.first_non_null_local(rp$boundary, "drop"))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg_eval$burden_log_eps, 1e-12)), 0)
  sigma_burden_use <- as.numeric(.first_non_null_local(cfg_eval$sigma_burden, 0.35))
  if (!is.finite(sigma_burden_use) || sigma_burden_use <= 0) sigma_burden_use <- 0.35
  sigma_ploidy_use <- as.numeric(.first_non_null_local(cfg_eval$sigma_ploidy, 0.08))
  if (!is.finite(sigma_ploidy_use) || sigma_ploidy_use <= 0) sigma_ploidy_use <- 0.08
  mu_by_N <- vapply(
    model_core$grid_pre,
    function(n) weighted_ploidy_from_total_N(n, chr_lengths_bp = cfg_eval$chr_lengths_bp) * cfg_eval$N_UNIT,
    numeric(1)
  )

  comp <- cpp_o2simps_objective_components_map(
    cohort_code = as.integer(scenario_cpp$cohort_code),
    dose_vec = as.numeric(scenario_cpp$dose),
    treat_day_vec = as.numeric(scenario_cpp$treat_day),
    obs_steps_list = scenario_cpp$obs_steps,
    sim_end_step_vec = as.integer(scenario_cpp$sim_end_step),
    obs_burden_list = scenario_cpp$obs_burden,
    keep_burden_list = scenario_cpp$keep_burden,
    ploidy_z_list = scenario_cpp$ploidy_z,
    mu_by_N = as.numeric(mu_by_N),
    sigma_burden = as.numeric(sigma_burden_use),
    sigma_ploidy = as.numeric(sigma_ploidy_use),
    init_state_2N = as.numeric(model_core$init_state_2N),
    init_state_4N = as.numeric(model_core$init_state_4N),
    N0min = as.integer(cfg_eval$N_MIN),
    N0max = as.integer(cfg_eval$N_MAX),
    N1min = as.integer(cfg_eval$N_MIN),
    N1max = as.integer(cfg_eval$N_MAX),
    DT = as.numeric(cfg_eval$DT),
    dose_ref = as.numeric(cfg_eval$dose_ref),
    fit_treatment = isTRUE(cfg_eval$fit_treatment),
    alpha = as.numeric(.first_non_null_local(rp$alpha, 0.0)),
    gamma = as.numeric(.first_non_null_local(rp$gamma, 1.0)),
    tx_mult_min = as.numeric(cfg_eval$tx_mult_min),
    crowding = as.character(cfg_eval$crowding),
    K = as.numeric(cfg_eval$K),
    min_pop = as.numeric(cfg_eval$min_pop),
    O2_cap = as.numeric(o2_cap_use),
    o2_feedback = isTRUE(.first_non_null_local(cfg_eval$o2_burden_feedback, TRUE)),
    o2_curve_type = as.character(o2_curve_type_use),
    o2_init = as.numeric(o2_init_use),
    o2_rate = as.numeric(o2_rate_use),
    o2_shape_v = as.numeric(o2_shape_v_use),
    tau_O2 = as.numeric(tau_O2_use),
    o2_anchor_N = as.numeric(o2_anchor_use),
    o2_logN_eps = as.numeric(.first_non_null_local(cfg_eval$o2_logN_eps, 1.0)),
    o2_cache_bin_pct = as.numeric(.first_non_null_local(cfg_eval$o2_cache_bin_pct, 0.01)),
    o2_cache_hysteresis_pct = as.numeric(.first_non_null_local(cfg_eval$o2_cache_hysteresis_pct, 0.005)),
    o2_cache_profile = isTRUE(.first_non_null_local(cfg_eval$o2_cache_profile, FALSE)),
    lam_min = as.numeric(rp$lam_min),
    lam_max = as.numeric(rp$lam_max),
    k_o = as.numeric(rp$k_o),
    has_p_misseg = !is.null(rp$p_misseg),
    p_misseg = as.numeric(.first_non_null_local(rp$p_misseg, 0.0)),
    k_o_mis = as.numeric(.first_non_null_local(rp$k_o_mis, 50.0)),
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = 0.0,
    p_wgd = as.numeric(p_wgd_use),
    boundary = boundary_mode,
    eps_tail = as.numeric(1e-8),
    beta_buffer = as.numeric(.first_non_null_local(rp$beta_buffer, 0.0)),
    n_exp = as.numeric(.first_non_null_local(rp$n_exp, 1.0)),
    smax = as.numeric(.first_non_null_local(rp$smax, 1.0)),
    N_unit = as.integer(cfg_eval$N_UNIT),
    beta_size = as.numeric(.first_non_null_local(rp$beta_size, cfg_eval$prior_center_beta_size, default_beta_size_prior_center())),
    alpha_o2 = as.numeric(.first_non_null_local(rp$alpha_o2, cfg_eval$alpha_o2_init, 0.5)),
    o2_ref_pct = as.numeric(.first_non_null_local(rp$o2_ref_pct, cfg_eval$o2_ref_pct_init, 2.5)),
    gamma_growth = as.numeric(.first_non_null_local(rp$gamma_growth, cfg_eval$gamma_growth_init, 2.0)),
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(burden_floor),
    burden_log_eps = as.numeric(.first_non_null_local(cfg_eval$burden_log_eps, 1e-12))
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
    n_burden = as.integer(.first_non_null_local(comp$n_burden, 0L)),
    n_ploidy = as.integer(.first_non_null_local(comp$n_ploidy, 0L)),
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
    n_burden = raw$n_burden,
    n_ploidy = raw$n_ploidy,
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
run_optimizer <- function(objective_fn, lower, upper, cfg, argv, stage_label = "fit", init_par = NULL) {
  n_cores <- as.integer(max(1L, ifelse(is.finite(cfg$n_cores), cfg$n_cores, 1L)))
  use_deoptim <- isTRUE(cfg$use_deoptim)
  deoptim_parallel <- isTRUE(cfg$deoptim_parallel)
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

        load_mode <- "unknown"
        last_err <- NULL
        loaded_ok <- FALSE

        # Path 1: use pre-generated wrapper/DLL from main process if both are visible.
        if (nzchar(wrapper_path) && nzchar(dll_path) && file.exists(wrapper_path) && file.exists(dll_path)) {
          max_attempts <- 20L
          for (attempt in seq_len(max_attempts)) {
            loaded_ok <- tryCatch({
              source(wrapper_path, local = .GlobalEnv)
              if (!has_required()) stop("Missing required C++ wrappers after source(wrapper_path).")
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
            if (!has_required()) stop("Missing required C++ wrappers after source(model_path).")
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
      "clip",
      ".first_non_null_local",
      "default_beta_size_prior_center",
      "default_rho_2N_prior_bounds",
      "default_rho_2N_prior_center",
      "default_chr_lengths_bp_1to22",
      "normalize_chr_lengths_bp_1to22",
      "weighted_ploidy_from_total_N",
      "map_ploidy_to_N_by_chrlen",
      "cell_volume_mm3_by_ploidy",
      "cell_volume_mm3_by_N",
      "burden_volume_mm3_from_state"
    )
    export_global <- export_global[export_global %in% ls(.GlobalEnv, all.names = TRUE)]
    if (length(export_global) > 0) {
      parallel::clusterExport(cl, varlist = export_global, envir = .GlobalEnv)
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
      ", n_cores=", n_cores
    )
    de_ctrl <- list(
      trace = TRUE,
      itermax = cfg$itermax,
      NP = cfg$NP,
      strategy = 2,
      reltol = cfg$de_reltol,
      steptol = cfg$de_steptol
    )
    de_ctrl$initialpop <- de_init_plan$pop
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
    optim_res <- tryCatch(
      DEoptim::DEoptim(
        fn = objective_fn,
        lower = lower,
        upper = upper,
        control = de_ctrl
      ),
      error = function(e) {
        stop("[", stage_label, "] DEoptim failed: ", conditionMessage(e))
      }
    )
    optim_res$method <- if (de_active > 1L) "DEoptim_parallel" else "DEoptim_serial"
    optim_res$parallel_info <- list(
      requested_workers = de_requested,
      started_workers = de_started,
      active_workers = de_active
    )
    best_par <- optim_res$optim$bestmem
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
  list(best_par = best_par, optim_res = optim_res)
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
    pred_vol <- as.numeric(sim$Vmm3_obs)
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
      pred_burden_volume_mm3 = pred_vol,
      obs_log_burden = ifelse(is.finite(obs) & obs >= 0, log(pmax(obs, log_eps)), NA_real_),
      pred_log_burden = ifelse(is.finite(pred_vol) & pred_vol >= 0, log(pmax(pred_vol, log_eps)), NA_real_),
      obs_norm = obs_norm,
      pred_norm = pred_norm
    )

    ploidy_df <- NULL
    if (length(sc$ploidy_obs_z) > 0) {
      obs_N <- map_ploidy_to_N_by_chrlen(
        ploidy_values = as.numeric(sc$ploidy_obs_z) / as.numeric(cfg$N_UNIT),
        N_grid = cfg$N_MIN:cfg$N_MAX,
        chr_lengths_bp = cfg$chr_lengths_bp
      )
      obs_N <- as.integer(clip(obs_N, cfg$N_MIN, cfg$N_MAX))
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
# Function: main
# Purpose: Entry point: parse options, run fitting/visualization workflow, and write outputs.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  require_cli_args(argv, c(
    "data_dir", "n_cores", "use_deoptim", "deoptim_parallel",
    "fit_treatment", "dose_zero_only", "paired_only", "truncate_at_treatment", "ploidy_at_harvest",
    "itermax", "NP", "n_starts", "optim_maxit",
    "sigma_burden", "sigma_ploidy", "burden_log_eps", "burden_exclude_day0",
    "use_soft_prior", "lambda_prior",
    "tau_O2", "tau_O2_init", "tau_O2_min", "tau_O2_max",
    "o2_curve_type", "o2_cap_pct", "o2_anchor_N",
    "o2_init_pct_init", "o2_init_pct_min", "o2_init_pct_max",
    "o2_rate_init", "o2_rate_min", "o2_rate_max",
    "o2_shape_v_init", "o2_shape_v_min", "o2_shape_v_max",
    "alpha_o2_init", "alpha_o2_min", "alpha_o2_max",
    "o2_ref_pct_init", "o2_ref_pct_min", "o2_ref_pct_max",
    "gamma_growth_init", "gamma_growth_min", "gamma_growth_max",
    "parameter_table",
    "N_UNIT", "N_MIN", "N_MAX", "dt", "o2_burden_feedback",
    "o2_logn_eps", "o2_cache_bin_pct", "o2_cache_hysteresis_pct", "o2_cache_profile",
    "K", "crowding", "init_total_size", "dose_ref", "tx_mult_min", "min_pop",
    "rho_2N_min", "rho_2N_max",
    "prior_center_log10_k_o", "prior_sd_log10_k_o",
    "prior_center_log10_o2_rate", "prior_sd_log10_o2_rate",
    "prior_center_log10_o2_init_pct", "prior_sd_log10_o2_init_pct",
    "prior_center_log10_o2_shape_v", "prior_sd_log10_o2_shape_v",
    "prior_center_beta_size", "prior_sd_beta_size",
    "prior_center_log10_n_exp", "prior_sd_log10_n_exp",
    "prior_center_log10_rho_2N", "prior_sd_log10_rho_2N",
    "optim_trace", "optim_trace_every", "trace_obj",
    "de_init_mode", "de_init_uniform_frac", "de_init_sigma_frac", "de_reltol", "de_steptol",
    "predict_n_cores", "max_scenarios", "seed"
  ))
  script_dir <- get_script_dir()

  model_path <- file.path(script_dir, "model_O2_GLF_MAP.R")
  if (!file.exists(model_path)) stop("Cannot find model_O2_GLF_MAP.R at ", model_path)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = script_dir)
  source(model_path)
  required_cpp_fit <- c("cpp_o2simps_build_G_for_o2_triplet", "cpp_o2simps_simulate_one", "cpp_o2simps_objective_components_map")
  missing_cpp_fit <- required_cpp_fit[!vapply(required_cpp_fit, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing_cpp_fit) > 0L) {
    stop("Required C++ symbols missing for fit path: ", paste(missing_cpp_fit, collapse = ", "))
  }
  cpp_dll <- tryCatch(
    o2simps_cpp_dll_info(),
    error = function(e) stop("Failed to resolve compiled O2_GLF_MAP DLL info: ", conditionMessage(e))
  )

  default_data_dir <- normalizePath(file.path(script_dir, "..", "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
  default_parameter_table <- normalizePath(file.path(script_dir, "..", "..", "data", "O2_GLF", "parameter_table.csv"), mustWork = FALSE)
  data_dir <- if (!is.null(argv$data_dir)) argv$data_dir else default_data_dir
  truncate_at_treatment <- as_bool(argv$truncate_at_treatment, FALSE)
  n_cores_arg <- as_int(argv$n_cores, NA_integer_)
  n_cores_use <- if (is.finite(n_cores_arg)) n_cores_arg else default_n_cores()
  o2_cap_arg <- as_num(argv$o2_cap_pct, 5.0)
  if (!is.finite(o2_cap_arg) || o2_cap_arg <= 0 || o2_cap_arg > 100) o2_cap_arg <- 5.0

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
    o2_cap_pct = o2_cap_arg,
    o2_curve_type = tolower(trimws(as.character(.first_non_null_local(argv$o2_curve_type, "gompertz")))),
    o2_burden_feedback = as_bool(argv$o2_burden_feedback, TRUE),
    o2_logN_eps = as_num(argv$o2_logn_eps, 1.0),
    o2_cache_bin_pct = as_num(argv$o2_cache_bin_pct, 0.01),
    o2_cache_hysteresis_pct = as_num(argv$o2_cache_hysteresis_pct, 0.005),
    o2_cache_profile = as_bool(argv$o2_cache_profile, FALSE),
    o2_anchor_N = as_num(argv$o2_anchor_N, as_num(argv$init_total_size, 1e6)),
    o2_init_pct_init = as_num(argv$o2_init_pct_init, min(0.5, o2_cap_arg * 0.5)),
    o2_init_pct_min = as_num(argv$o2_init_pct_min, 1e-3),
    o2_init_pct_max = as_num(argv$o2_init_pct_max, max(1e-3, o2_cap_arg - 1e-3)),
    o2_rate_init = as_num(argv$o2_rate_init, 1.0),
    o2_rate_min = as_num(argv$o2_rate_min, 1e-3),
    o2_rate_max = as_num(argv$o2_rate_max, 1e2),
    o2_shape_v_init = as_num(argv$o2_shape_v_init, 1.0),
    o2_shape_v_min = as_num(argv$o2_shape_v_min, 1e-2),
    o2_shape_v_max = as_num(argv$o2_shape_v_max, 20),
    alpha_o2_init = as_num(argv$alpha_o2_init, 0.5),
    alpha_o2_min = as_num(argv$alpha_o2_min, 1e-2),
    alpha_o2_max = as_num(argv$alpha_o2_max, 10),
    o2_ref_pct_init = as_num(argv$o2_ref_pct_init, 2.5),
    o2_ref_pct_min = as_num(argv$o2_ref_pct_min, 0),
    o2_ref_pct_max = as_num(argv$o2_ref_pct_max, 5),
    gamma_growth_init = as_num(argv$gamma_growth_init, 2.0),
    gamma_growth_min = as_num(argv$gamma_growth_min, 2.0),
    gamma_growth_max = as_num(argv$gamma_growth_max, 2.0),
    tau_O2 = as_num(argv$tau_O2, NA_real_),
    tau_O2_init = as_num(argv$tau_O2_init, 2.0),
    tau_O2_min = as_num(argv$tau_O2_min, 1e-3),
    tau_O2_max = as_num(argv$tau_O2_max, 1e3),
    fit_tau_O2 = !is.finite(as_num(argv$tau_O2, NA_real_)),
    parameter_table = if (!is.null(argv$parameter_table)) argv$parameter_table else default_parameter_table,
    K = as_num(argv$K, 1e12),
    crowding = if (!is.null(argv$crowding)) argv$crowding else "logistic",
    init_total_size = as_num(argv$init_total_size, 1e6),
    dose_ref = as_num(argv$dose_ref, 30),
    tx_mult_min = as_num(argv$tx_mult_min, 0.05),
    min_pop = as_num(argv$min_pop, 1e-12),
    # objective settings (MAP likelihood)
    sigma_burden = as_num(argv$sigma_burden, 0.35),
    sigma_ploidy = as_num(argv$sigma_ploidy, 0.08),
    burden_log_eps = as_num(argv$burden_log_eps, 1e-12),
    burden_exclude_day0 = as_bool(argv$burden_exclude_day0, TRUE),
    rho_2N_min = as_num(argv$rho_2N_min, 3.2e4),
    rho_2N_max = as_num(argv$rho_2N_max, 5.6e4),
    use_soft_prior = as_bool(argv$use_soft_prior, TRUE),
    lambda_prior = as_num(argv$lambda_prior, 0.1),
    prior_center_log10_k_o = as_num(argv$prior_center_log10_k_o, log10(50)),
    prior_sd_log10_k_o = as_num(argv$prior_sd_log10_k_o, 1.0),
    prior_center_log10_o2_rate = as_num(argv$prior_center_log10_o2_rate, log10(as_num(argv$o2_rate_init, 1.0))),
    prior_sd_log10_o2_rate = as_num(argv$prior_sd_log10_o2_rate, 1.0),
    prior_center_log10_o2_init_pct = as_num(argv$prior_center_log10_o2_init_pct, log10(as_num(argv$o2_init_pct_init, 0.5))),
    prior_sd_log10_o2_init_pct = as_num(argv$prior_sd_log10_o2_init_pct, 0.5),
    prior_center_log10_o2_shape_v = as_num(argv$prior_center_log10_o2_shape_v, log10(as_num(argv$o2_shape_v_init, 1.0))),
    prior_sd_log10_o2_shape_v = as_num(argv$prior_sd_log10_o2_shape_v, 0.5),
    prior_center_beta_size = as_num(argv$prior_center_beta_size, default_beta_size_prior_center()),
    prior_sd_beta_size = as_num(argv$prior_sd_beta_size, 0.5),
    prior_center_log10_n_exp = as_num(argv$prior_center_log10_n_exp, 0.0),
    prior_sd_log10_n_exp = as_num(argv$prior_sd_log10_n_exp, 0.5),
    prior_center_log10_rho_2N = as_num(argv$prior_center_log10_rho_2N, log10(sqrt(as_num(argv$rho_2N_min, 3.2e4) * as_num(argv$rho_2N_max, 5.6e4)))),
    prior_sd_log10_rho_2N = as_num(argv$prior_sd_log10_rho_2N, 0.35),
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
    de_reltol = as_num(argv$de_reltol, 1e-3),
    de_steptol = as_int(argv$de_steptol, 25L),
    itermax = as_int(argv$itermax, 40L),
    NP = as_int(argv$NP, 80L),
    n_cores = n_cores_use,
    # Keep post-fit prediction stable: use serial by default unless explicitly requested.
    predict_n_cores = as_int(argv$predict_n_cores, 1L),
    seed = as_int(argv$seed, 1L),
    max_scenarios = as_num(argv$max_scenarios, Inf)
  )

  if (!cfg$crowding %in% c("logistic", "gompertz")) stop("crowding must be logistic or gompertz")
  if (cfg$DT <= 0) stop("dt must be > 0")
  if (cfg$N_MAX < cfg$N_MIN) stop("N_MAX must be >= N_MIN")
  if (!is.finite(cfg$o2_logN_eps) || cfg$o2_logN_eps <= 0) stop("o2_logN_eps must be > 0")
  if (!is.finite(cfg$o2_cache_bin_pct) || cfg$o2_cache_bin_pct <= 0 || cfg$o2_cache_bin_pct > 100) {
    stop("o2_cache_bin_pct must be in (0,100].")
  }
  if (!is.finite(cfg$o2_cache_hysteresis_pct) || cfg$o2_cache_hysteresis_pct < 0 || cfg$o2_cache_hysteresis_pct > 100) {
    stop("o2_cache_hysteresis_pct must be in [0,100].")
  }
  if (!cfg$o2_curve_type %in% c("gompertz", "glogistic")) stop("o2_curve_type must be gompertz or glogistic")
  if (!is.finite(cfg$o2_cap_pct) || cfg$o2_cap_pct <= 0 || cfg$o2_cap_pct > 100) stop("o2_cap_pct must be in (0,100]")
  if (!is.finite(cfg$o2_anchor_N) || cfg$o2_anchor_N < 0) stop("o2_anchor_N must be >= 0")
  if (!is.finite(cfg$o2_init_pct_init) || cfg$o2_init_pct_init <= 0 || cfg$o2_init_pct_init >= cfg$o2_cap_pct) stop("o2_init_pct_init must be in (0,o2_cap_pct)")
  if (!is.finite(cfg$o2_init_pct_min) || cfg$o2_init_pct_min <= 0) stop("o2_init_pct_min must be > 0")
  if (!is.finite(cfg$o2_init_pct_max) || cfg$o2_init_pct_max <= 0) stop("o2_init_pct_max must be > 0")
  if (cfg$o2_init_pct_max >= cfg$o2_cap_pct) stop("o2_init_pct_max must be < o2_cap_pct")
  if (cfg$o2_init_pct_max < cfg$o2_init_pct_min) stop("o2_init_pct_max must be >= o2_init_pct_min")
  if (!is.finite(cfg$o2_rate_init) || cfg$o2_rate_init <= 0) stop("o2_rate_init must be > 0")
  if (!is.finite(cfg$o2_rate_min) || cfg$o2_rate_min <= 0) stop("o2_rate_min must be > 0")
  if (!is.finite(cfg$o2_rate_max) || cfg$o2_rate_max <= 0) stop("o2_rate_max must be > 0")
  if (cfg$o2_rate_max < cfg$o2_rate_min) stop("o2_rate_max must be >= o2_rate_min")
  if (!is.finite(cfg$o2_shape_v_init) || cfg$o2_shape_v_init <= 0) stop("o2_shape_v_init must be > 0")
  if (!is.finite(cfg$o2_shape_v_min) || cfg$o2_shape_v_min <= 0) stop("o2_shape_v_min must be > 0")
  if (!is.finite(cfg$o2_shape_v_max) || cfg$o2_shape_v_max <= 0) stop("o2_shape_v_max must be > 0")
  if (cfg$o2_shape_v_max < cfg$o2_shape_v_min) stop("o2_shape_v_max must be >= o2_shape_v_min")
  if (!is.finite(cfg$alpha_o2_init) || cfg$alpha_o2_init <= 0) stop("alpha_o2_init must be > 0")
  if (!is.finite(cfg$alpha_o2_min) || cfg$alpha_o2_min <= 0) stop("alpha_o2_min must be > 0")
  if (!is.finite(cfg$alpha_o2_max) || cfg$alpha_o2_max <= 0) stop("alpha_o2_max must be > 0")
  if (cfg$alpha_o2_max < cfg$alpha_o2_min) stop("alpha_o2_max must be >= alpha_o2_min")
  if (!is.finite(cfg$o2_ref_pct_init) || cfg$o2_ref_pct_init < 0 || cfg$o2_ref_pct_init > 100) stop("o2_ref_pct_init must be in [0,100]")
  if (!is.finite(cfg$o2_ref_pct_min) || cfg$o2_ref_pct_min < 0 || cfg$o2_ref_pct_min > 100) stop("o2_ref_pct_min must be in [0,100]")
  if (!is.finite(cfg$o2_ref_pct_max) || cfg$o2_ref_pct_max < 0 || cfg$o2_ref_pct_max > 100) stop("o2_ref_pct_max must be in [0,100]")
  if (cfg$o2_ref_pct_max < cfg$o2_ref_pct_min) stop("o2_ref_pct_max must be >= o2_ref_pct_min")
  if (!is.finite(cfg$gamma_growth_init) || cfg$gamma_growth_init <= 0) stop("gamma_growth_init must be > 0")
  if (!is.finite(cfg$gamma_growth_min) || cfg$gamma_growth_min <= 0) stop("gamma_growth_min must be > 0")
  if (!is.finite(cfg$gamma_growth_max) || cfg$gamma_growth_max <= 0) stop("gamma_growth_max must be > 0")
  if (cfg$gamma_growth_max < cfg$gamma_growth_min) stop("gamma_growth_max must be >= gamma_growth_min")
  if (!is.finite(cfg$tau_O2_init) || cfg$tau_O2_init <= 0) stop("tau_O2_init must be > 0")
  if (!is.finite(cfg$tau_O2_min) || cfg$tau_O2_min <= 0) stop("tau_O2_min must be > 0")
  if (!is.finite(cfg$tau_O2_max) || cfg$tau_O2_max <= 0) stop("tau_O2_max must be > 0")
  if (cfg$tau_O2_max < cfg$tau_O2_min) stop("tau_O2_max must be >= tau_O2_min")
  if (!isTRUE(cfg$fit_tau_O2)) {
    if (!is.finite(cfg$tau_O2) || cfg$tau_O2 <= 0) stop("tau_O2 must be > 0 when provided.")
  } else {
    cfg$tau_O2 <- NA_real_
  }
  if (cfg$o2_curve_type != "glogistic") {
    v_fix <- max(cfg$o2_shape_v_init, 1e-6)
    cfg$o2_shape_v_min <- v_fix * (1 - 1e-6)
    cfg$o2_shape_v_max <- v_fix * (1 + 1e-6)
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
  if (!is.finite(cfg$sigma_ploidy) || cfg$sigma_ploidy <= 0) stop("sigma_ploidy must be > 0")
  if (!is.finite(cfg$rho_2N_min) || cfg$rho_2N_min <= 0) stop("rho_2N_min must be > 0")
  if (!is.finite(cfg$rho_2N_max) || cfg$rho_2N_max <= 0) stop("rho_2N_max must be > 0")
  if (cfg$rho_2N_max < cfg$rho_2N_min) stop("rho_2N_max must be >= rho_2N_min")
  if (!is.finite(cfg$lambda_prior) || cfg$lambda_prior < 0) stop("lambda_prior must be >= 0")
  if (!is.finite(cfg$prior_sd_log10_k_o) || cfg$prior_sd_log10_k_o <= 0) stop("prior_sd_log10_k_o must be > 0")
  if (!is.finite(cfg$prior_sd_log10_o2_rate) || cfg$prior_sd_log10_o2_rate <= 0) stop("prior_sd_log10_o2_rate must be > 0")
  if (!is.finite(cfg$prior_sd_log10_o2_init_pct) || cfg$prior_sd_log10_o2_init_pct <= 0) stop("prior_sd_log10_o2_init_pct must be > 0")
  if (!is.finite(cfg$prior_sd_log10_o2_shape_v) || cfg$prior_sd_log10_o2_shape_v <= 0) stop("prior_sd_log10_o2_shape_v must be > 0")
  if (!is.finite(cfg$prior_sd_beta_size) || cfg$prior_sd_beta_size <= 0) stop("prior_sd_beta_size must be > 0")
  if (!is.finite(cfg$prior_sd_log10_n_exp) || cfg$prior_sd_log10_n_exp <= 0) stop("prior_sd_log10_n_exp must be > 0")
  if (!is.finite(cfg$prior_sd_log10_rho_2N) || cfg$prior_sd_log10_rho_2N <= 0) stop("prior_sd_log10_rho_2N must be > 0")
  if (!cfg$use_deoptim && cfg$deoptim_parallel) stop("deoptim_parallel=TRUE requires use_deoptim=TRUE")

  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- file.path(data_dir, "all_ploidy.tsv")
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)
  cfg$model_core <- build_model_core(cfg = cfg)
  cfg$scenario_cpp <- prepare_cpp_scenarios(scenarios, cfg)

  full_names <- get_param_names(
    fit_treatment = isTRUE(cfg$fit_treatment),
    fit_tau_O2 = isTRUE(cfg$fit_tau_O2)
  )
  param_table <- read_parameter_table_transformed(cfg$parameter_table, full_names)
  bounds <- list(lower = param_table$lower, upper = param_table$upper)
  default_par_t <- param_table$init
  init_params_tsv <- if (!is.null(argv$init_params_tsv)) argv$init_params_tsv else NULL
  warm_start_t <- if (!is.null(init_params_tsv)) {
    read_init_params_t(init_params_tsv, bounds = bounds, cfg = cfg)
  } else {
    NULL
  }
  message(
    "MAP likelihood objective enabled: burden=lognormal NLL (per-tumor mean), ",
    "ploidy=continuous single-cell mixture NLL (2N/4N group-balanced mean, 0.5/0.5), ",
    "practical weighting default (equal tumor weighting); sigma_burden=", signif(cfg$sigma_burden, 6),
    ", sigma_ploidy=", signif(cfg$sigma_ploidy, 6)
  )
  if (isTRUE(cfg$use_soft_prior) && cfg$lambda_prior > 0) {
    message(
      "Soft prior enabled: lambda_prior=", signif(cfg$lambda_prior, 6),
      "; centers(log10_k_o, log10_o2_rate, log10_o2_init_pct, log10_o2_shape_v, beta_size, log10_n_exp, log10_rho_2N)=(",
      signif(cfg$prior_center_log10_k_o, 6), ", ",
      signif(cfg$prior_center_log10_o2_rate, 6), ", ",
      signif(cfg$prior_center_log10_o2_init_pct, 6), ", ",
      signif(cfg$prior_center_log10_o2_shape_v, 6), ", ",
      signif(cfg$prior_center_beta_size, 6), ", ",
      signif(cfg$prior_center_log10_n_exp, 6), ", ",
      signif(cfg$prior_center_log10_rho_2N, 6), ")"
    )
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
    "Burden observation model enabled: log-normal likelihood on V(mm^3), ",
    "V_pred = sum_n n_n * [(1/rho_2N) * (P/2)^beta_size], ",
    "rho_2N_range=[", signif(cfg$rho_2N_min, 6), ", ", signif(cfg$rho_2N_max, 6), "] cells/mm^3",
    "; burden_exclude_day0=", if (isTRUE(cfg$burden_exclude_day0)) "TRUE" else "FALSE"
  )
  message(
    "O2 mode: ",
    if (isTRUE(cfg$o2_burden_feedback)) "dynamic feedback" else "fixed",
    "; o2_curve_type=", cfg$o2_curve_type,
    ", o2_cap_pct=", signif(cfg$o2_cap_pct, 6),
    ", tau_O2_mode=", if (isTRUE(cfg$fit_tau_O2)) "fit" else "fixed",
    ", tau_O2=", if (isTRUE(cfg$fit_tau_O2)) {
      paste0("init=", signif(cfg$tau_O2_init, 6), ",range=[", signif(cfg$tau_O2_min, 6), ",", signif(cfg$tau_O2_max, 6), "]")
    } else {
      signif(cfg$tau_O2, 6)
    },
    ", o2_init_pct_init=", signif(cfg$o2_init_pct_init, 6),
    ", o2_rate_init=", signif(cfg$o2_rate_init, 6),
    ", o2_shape_v_init=", signif(cfg$o2_shape_v_init, 6),
    if (isTRUE(cfg$o2_burden_feedback)) {
      "; fitted O2 params: o2_init_pct, o2_rate, o2_shape_v"
    } else {
      "; O2 fixed at o2_cap_pct."
    }
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
  message("[single_stage] pass1 initialization: candidate[1]=initial_value, candidates[2:NP]=uniform(lower,upper).")
  objective_fn <- function(par) evaluate_objective(par, scenarios = scenarios, cfg = pass_cfg)
  single_fit <- run_optimizer(
    objective_fn = objective_fn,
    lower = bounds$lower,
    upper = bounds$upper,
    cfg = pass_cfg,
    argv = argv,
    stage_label = "single_stage",
    init_par = initial_par_t
  )
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
    file.path(script_dir, "..", "..", "results", paste0("fit_invivo_model_O2_GLF_MAP_", run_stamp))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  params_df <- data.frame(
    parameter = names(best_par),
    value = as.numeric(best_par),
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
      row.names = NULL
    )
  }))
  write.table(pass_df, file = file.path(out_dir, "single_stage_pass_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  n_ploidy_scenarios <- sum(vapply(scenarios, function(s) length(s$ploidy_obs_z) > 0, logical(1)))
  fit_mode <- if (!is.null(optim_res$mode)) as.character(optim_res$mode) else "single_stage"
  summary_df <- data.frame(
    metric = c(
      "fit_mode",
      "objective",
      "objective_data",
      "objective_prior_raw",
      "objective_prior",
      "objective_burden",
      "objective_ploidy",
      "burden_exclude_day0",
      "sigma_burden",
      "sigma_ploidy",
      "use_soft_prior",
      "lambda_prior",
      "prior_center_log10_k_o",
      "prior_sd_log10_k_o",
      "prior_center_log10_o2_rate",
      "prior_sd_log10_o2_rate",
      "prior_center_log10_o2_init_pct",
      "prior_sd_log10_o2_init_pct",
      "prior_center_log10_o2_shape_v",
      "prior_sd_log10_o2_shape_v",
      "prior_center_beta_size",
      "prior_sd_beta_size",
      "prior_center_log10_n_exp",
      "prior_sd_log10_n_exp",
      "prior_center_log10_rho_2N",
      "prior_sd_log10_rho_2N",
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
      "fit_treatment",
      "o2_burden_feedback",
      "o2_logN_eps",
      "o2_cache_bin_pct",
      "o2_cache_hysteresis_pct",
      "o2_cache_profile",
      "o2_curve_type",
      "o2_cap_pct",
      "o2_anchor_N",
      "o2_init_pct_init",
      "o2_init_pct_min",
      "o2_init_pct_max",
      "o2_rate_init",
      "o2_rate_min",
      "o2_rate_max",
      "o2_shape_v_init",
      "o2_shape_v_min",
      "o2_shape_v_max",
      "alpha_o2_init",
      "alpha_o2_min",
      "alpha_o2_max",
      "o2_ref_pct_init",
      "o2_ref_pct_min",
      "o2_ref_pct_max",
      "gamma_growth_init",
      "gamma_growth_min",
      "gamma_growth_max",
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
      as.character(final_obj),
      as.character(final_comp$L_data),
      as.character(final_comp$L_prior_raw),
      as.character(final_comp$L_prior),
      as.character(final_comp$L_b),
      as.character(final_comp$L_p),
      as.character(cfg$burden_exclude_day0),
      as.character(cfg$sigma_burden),
      as.character(cfg$sigma_ploidy),
      as.character(cfg$use_soft_prior),
      as.character(cfg$lambda_prior),
      as.character(cfg$prior_center_log10_k_o),
      as.character(cfg$prior_sd_log10_k_o),
      as.character(cfg$prior_center_log10_o2_rate),
      as.character(cfg$prior_sd_log10_o2_rate),
      as.character(cfg$prior_center_log10_o2_init_pct),
      as.character(cfg$prior_sd_log10_o2_init_pct),
      as.character(cfg$prior_center_log10_o2_shape_v),
      as.character(cfg$prior_sd_log10_o2_shape_v),
      as.character(cfg$prior_center_beta_size),
      as.character(cfg$prior_sd_beta_size),
      as.character(cfg$prior_center_log10_n_exp),
      as.character(cfg$prior_sd_log10_n_exp),
      as.character(cfg$prior_center_log10_rho_2N),
      as.character(cfg$prior_sd_log10_rho_2N),
      as.character(length(scenarios)),
      as.character(n_ploidy_scenarios),
      as.character(.first_non_null_local(final_comp$n_ploidy_2N, NA_integer_)),
      as.character(.first_non_null_local(final_comp$n_ploidy_4N, NA_integer_)),
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
      as.character(cfg$fit_treatment),
      as.character(cfg$o2_burden_feedback),
      as.character(cfg$o2_logN_eps),
      as.character(cfg$o2_cache_bin_pct),
      as.character(cfg$o2_cache_hysteresis_pct),
      as.character(cfg$o2_cache_profile),
      as.character(cfg$o2_curve_type),
      as.character(cfg$o2_cap_pct),
      as.character(cfg$o2_anchor_N),
      as.character(cfg$o2_init_pct_init),
      as.character(cfg$o2_init_pct_min),
      as.character(cfg$o2_init_pct_max),
      as.character(cfg$o2_rate_init),
      as.character(cfg$o2_rate_min),
      as.character(cfg$o2_rate_max),
      as.character(cfg$o2_shape_v_init),
      as.character(cfg$o2_shape_v_min),
      as.character(cfg$o2_shape_v_max),
      as.character(cfg$alpha_o2_init),
      as.character(cfg$alpha_o2_min),
      as.character(cfg$alpha_o2_max),
      as.character(cfg$o2_ref_pct_init),
      as.character(cfg$o2_ref_pct_min),
      as.character(cfg$o2_ref_pct_max),
      as.character(cfg$gamma_growth_init),
      as.character(cfg$gamma_growth_min),
      as.character(cfg$gamma_growth_max),
      as.character(cfg$fit_tau_O2),
      as.character(if (isTRUE(cfg$fit_tau_O2)) NA_real_ else cfg$tau_O2),
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
  write.table(summary_df, file = file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  write.table(preds$burden, file = file.path(out_dir, "burden_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(preds$ploidy, file = file.path(out_dir, "terminal_ploidy_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(optim_res, file = file.path(out_dir, "deoptim_result.rds"))
  saveRDS(cfg, file = file.path(out_dir, "fit_config.rds"))

  message("Done. Results written to: ", normalizePath(out_dir))
  message("Best objective: ", signif(final_obj, 6))
  message("Best parameters:")
  print(params_df)
}

if (sys.nframe() == 0) {
  setwd(get_script_dir())
  main()
}
