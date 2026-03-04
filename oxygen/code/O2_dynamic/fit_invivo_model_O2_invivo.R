#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readxl))

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

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.numeric(x))
}

as_num_vec <- function(x, default = NA_real_) {
  if (is.null(x)) return(as.numeric(default))
  s <- trimws(as.character(x))
  if (!nzchar(s)) return(as.numeric(default))
  parts <- unlist(strsplit(s, "[,;]", perl = TRUE))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(as.numeric(default))
  vals <- suppressWarnings(as.numeric(parts))
  if (any(!is.finite(vals))) {
    stop("Invalid numeric vector argument: ", x)
  }
  vals
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.integer(x))
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

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

.first_non_null_local <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

huber_mean <- function(r, k = 0.1) {
  a <- abs(r)
  mean(ifelse(a <= k, 0.5 * r^2, k * (a - 0.5 * k)))
}

weighted_mean_safe <- function(x, w = NULL) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(0)
  if (is.null(w)) return(mean(x))
  w <- as.numeric(w)
  if (length(w) != length(x)) return(mean(x))
  keep <- is.finite(x) & is.finite(w) & (w > 0)
  if (!any(keep)) return(mean(x[is.finite(x)]))
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_median_safe <- function(x, w = NULL) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(0)
  if (is.null(w)) return(stats::median(x))
  w <- as.numeric(w)
  if (length(w) != length(x)) return(stats::median(x))
  keep <- is.finite(x) & is.finite(w) & (w > 0)
  if (!any(keep)) return(stats::median(x[is.finite(x)]))
  xk <- x[keep]
  wk <- w[keep]
  ord <- order(xk)
  xk <- xk[ord]
  wk <- wk[ord]
  cw <- cumsum(wk) / sum(wk)
  xk[which(cw >= 0.5)[1]]
}

robust_huber_location <- function(x, w = NULL, k = 1.5, maxit = 50L, tol = 1e-8) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(0)
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w)
    if (length(w) != length(x)) w <- rep(1, length(x))
  }
  keep <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0L) return(0)
  if (length(x) == 1L) return(x[[1]])
  if (!is.finite(k) || k <= 0) k <- 1.5

  mu <- weighted_median_safe(x, w)
  s <- stats::mad(x, center = mu, constant = 1, na.rm = TRUE)
  if (!is.finite(s) || s <= 1e-12) return(weighted_mean_safe(x, w))

  for (i in seq_len(max(1L, as.integer(maxit)))) {
    u <- (x - mu) / (k * s)
    psi_over_u <- ifelse(abs(u) <= 1, 1, 1 / pmax(abs(u), 1e-12))
    ww <- w * psi_over_u
    sw <- sum(ww)
    if (!is.finite(sw) || sw <= 0) break
    mu_new <- sum(ww * x) / sw
    if (!is.finite(mu_new)) break
    if (abs(mu_new - mu) <= tol) {
      mu <- mu_new
      break
    }
    mu <- mu_new
  }
  mu
}

normalize_agg_method <- function(x, default = "huber") {
  m <- tolower(trimws(as.character(.first_non_null_local(x, default))))
  if (!m %in% c("mean", "median", "huber")) m <- default
  m
}

aggregate_scenario_losses <- function(losses, weights = NULL, method = "huber", huber_k = 1.5) {
  x <- as.numeric(losses)
  if (length(x) == 0L) return(0)
  keep <- is.finite(x)
  x <- x[keep]
  if (length(x) == 0L) return(0)

  w <- NULL
  if (!is.null(weights)) {
    wv <- as.numeric(weights)
    if (length(wv) == length(losses)) {
      wv <- wv[keep]
      keep_w <- is.finite(wv) & (wv > 0)
      if (any(keep_w)) {
        w <- wv
        w[!keep_w] <- 0
      }
    }
  }
  method <- normalize_agg_method(method, default = "huber")
  if (method == "mean") return(weighted_mean_safe(x, w))
  if (method == "median") return(weighted_median_safe(x, w))
  robust_huber_location(x, w = w, k = huber_k)
}

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
    p_names <- if (isTRUE(cfg$fit_treatment) && length(p) == 20L) {
      c(
        "log10_lam_min", "delta_lam", "log10_k_o", "log10_p_misseg",
        "log10_k_o_mis", "beta_buffer", "log10_n_exp", "log10_smax",
        "log10_p_wgd", "log10_K_down", "log10_h_down", "A_ang", "m_on",
        "log10_delta_m", "log10_s_on", "log10_s_off", "log10_rho_2N",
        "beta_size", "log10_alpha", "gamma"
      )
    } else if (length(p) == 18L) {
      c(
        "log10_lam_min", "delta_lam", "log10_k_o", "log10_p_misseg",
        "log10_k_o_mis", "beta_buffer", "log10_n_exp", "log10_smax",
        "log10_p_wgd", "log10_K_down", "log10_h_down", "A_ang", "m_on",
        "log10_delta_m", "log10_s_on", "log10_s_off", "log10_rho_2N",
        "beta_size"
      )
    } else {
      rep("", length(p))
    }
  }
  names(p) <- p_names

  centers <- c(
    log10_k_o = as.numeric(cfg$prior_center_log10_k_o),
    log10_K_down = as.numeric(cfg$prior_center_log10_K_down),
    log10_h_down = as.numeric(cfg$prior_center_log10_h_down),
    beta_size = as.numeric(cfg$prior_center_beta_size),
    log10_n_exp = as.numeric(cfg$prior_center_log10_n_exp),
    log10_rho_2N = as.numeric(cfg$prior_center_log10_rho_2N)
  )
  sds <- c(
    log10_k_o = as.numeric(cfg$prior_sd_log10_k_o),
    log10_K_down = as.numeric(cfg$prior_sd_log10_K_down),
    log10_h_down = as.numeric(cfg$prior_sd_log10_h_down),
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

default_beta_size_prior_center <- function() {
  log(1.5) / log(2)
}

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

default_rho_2N_prior_center <- function(cfg = NULL) {
  b <- default_rho_2N_prior_bounds(cfg)
  sqrt(b[["rho_2N_min"]] * b[["rho_2N_max"]])
}

cell_volume_mm3_by_ploidy <- function(ploidy, run_params, cfg) {
  p <- pmax(as.numeric(ploidy), 1e-8)
  rho_2N <- suppressWarnings(as.numeric(run_params$rho_2N))
  rho_2N <- if (length(rho_2N) > 0) rho_2N[[1]] else NA_real_
  if (is.na(rho_2N) || !is.finite(rho_2N) || rho_2N <= 0) rho_2N <- default_rho_2N_prior_center(cfg)
  beta_size <- as.numeric(.first_non_null_local(run_params$beta_size, default_beta_size_prior_center()))
  if (!is.finite(beta_size)) beta_size <- default_beta_size_prior_center()
  (1 / rho_2N) * (p / 2)^beta_size
}

cell_volume_mm3_by_N <- function(N, run_params, cfg) {
  cell_volume_mm3_by_ploidy(as.numeric(N) / as.numeric(cfg$N_UNIT), run_params = run_params, cfg = cfg)
}

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

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grepl("^--file=", args)]
  if (length(farg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", farg[[1]])))
}

default_n_cores <- function() {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(n) || is.na(n)) {
    n <- suppressWarnings(parallel::detectCores())
  }
  if (!is.finite(n) || is.na(n)) return(1L)
  as.integer(max(1L, n - 1L))
}

map_scenarios_parallel <- function(scenarios, n_cores = 1L, label = "predict", fn) {
  n <- length(scenarios)
  if (n == 0L) return(list())
  n_cores_use <- suppressWarnings(as.integer(n_cores))
  if (!is.finite(n_cores_use) || n_cores_use < 1L) n_cores_use <- 1L
  workers <- as.integer(max(1L, min(n_cores_use, n)))
  use_mc <- (workers > 1L) && (.Platform$OS.type != "windows")

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
      message("[", label, "] windows detected; using serial fallback for prediction.")
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

make_weight_schedule <- function(w_burden_vec, w_ploidy_vec) {
  wb <- as.numeric(w_burden_vec)
  wp <- as.numeric(w_ploidy_vec)
  if (length(wb) == 0 || length(wp) == 0) {
    stop("w_burden and w_ploidy must contain at least one value.")
  }
  nb <- length(wb)
  np <- length(wp)
  n <- max(nb, np)
  if (!(nb == np || nb == 1L || np == 1L)) {
    stop("w_burden and w_ploidy must have the same length, or one of them must have length 1.")
  }
  if (nb == 1L) wb <- rep(wb, n)
  if (np == 1L) wp <- rep(wp, n)
  if (any(!is.finite(wb)) || any(!is.finite(wp))) {
    stop("w_burden/w_ploidy schedules contain non-finite values.")
  }
  list(w_burden = wb, w_ploidy = wp, n_pass = n)
}

resolve_loss_scales <- function(cfg, default_par_t, warm_start_t, scenarios) {
  if (!isTRUE(cfg$loss_rescale)) {
    cfg$loss_scale_burden <- 1.0
    cfg$loss_scale_ploidy <- 1.0
    cfg$loss_scale_source <- "disabled"
    return(cfg)
  }

  sb <- cfg$loss_scale_burden
  sp <- cfg$loss_scale_ploidy
  has_sb <- is.finite(sb) && sb > 0
  has_sp <- is.finite(sp) && sp > 0

  if (has_sb && has_sp) {
    cfg$loss_scale_source <- "manual"
    return(cfg)
  }

  ref_par_t <- if (!is.null(warm_start_t)) warm_start_t else default_par_t
  ref_comp <- evaluate_objective_components_raw(ref_par_t, scenarios = scenarios, cfg = cfg)

  if (!has_sb) {
    sb <- max(ref_comp$L_b, cfg$loss_scale_eps)
    has_sb <- TRUE
  }
  if (!has_sp) {
    sp <- max(ref_comp$L_p, cfg$loss_scale_eps)
    has_sp <- TRUE
  }

  if (!(has_sb && has_sp)) {
    stop("Could not determine valid loss scales for burden/ploidy.")
  }

  cfg$loss_scale_burden <- sb
  cfg$loss_scale_ploidy <- sp
  cfg$loss_scale_source <- if (!is.null(warm_start_t)) "auto_warm_start" else "auto_midpoint"
  cfg
}

decode_params <- function(par_transformed, fit_treatment = TRUE) {
  if (isTRUE(fit_treatment)) {
    names(par_transformed) <- c(
      "log10_lam_min",
      "delta_lam",
      "log10_k_o",
      "log10_p_misseg",
      "log10_k_o_mis",
      "beta_buffer",
      "log10_n_exp",
      "log10_smax",
      "log10_p_wgd",
      "log10_K_down",
      "log10_h_down",
      "A_ang",
      "m_on",
      "log10_delta_m",
      "log10_s_on",
      "log10_s_off",
      "log10_rho_2N",
      "beta_size",
      "log10_alpha",
      "gamma"
    )
    lam_min <- 10^par_transformed["log10_lam_min"]
    lam_max <- lam_min + exp(par_transformed["delta_lam"])
    h_down <- 10^par_transformed["log10_h_down"]
    delta_m <- 10^par_transformed["log10_delta_m"]
    m_on <- par_transformed["m_on"]
    return(list(
      lam_min = lam_min,
      lam_max = lam_max,
      k_o = 10^par_transformed["log10_k_o"],
      p_misseg = 10^par_transformed["log10_p_misseg"],
      k_o_mis = 10^par_transformed["log10_k_o_mis"],
      beta_buffer = par_transformed["beta_buffer"],
      n_exp = 10^par_transformed["log10_n_exp"],
      smax = 10^par_transformed["log10_smax"],
      p_wgd = 10^par_transformed["log10_p_wgd"],
      K_down = 10^par_transformed["log10_K_down"],
      h_down = h_down,
      A_ang = par_transformed["A_ang"],
      m_on = m_on,
      delta_m = delta_m,
      m_off = m_on + delta_m,
      s_on = 10^par_transformed["log10_s_on"],
      s_off = 10^par_transformed["log10_s_off"],
      rho_2N = 10^par_transformed["log10_rho_2N"],
      beta_size = par_transformed["beta_size"],
      c_vol_2N_eff_mm3 = 10^-par_transformed["log10_rho_2N"],
      ratio_4N_2N = 2^par_transformed["beta_size"],
      alpha = 10^par_transformed["log10_alpha"],
      gamma = par_transformed["gamma"]
    ))
  }

  names(par_transformed) <- c(
    "log10_lam_min",
    "delta_lam",
    "log10_k_o",
    "log10_p_misseg",
    "log10_k_o_mis",
    "beta_buffer",
    "log10_n_exp",
    "log10_smax",
    "log10_p_wgd",
    "log10_K_down",
    "log10_h_down",
    "A_ang",
    "m_on",
    "log10_delta_m",
    "log10_s_on",
    "log10_s_off",
    "log10_rho_2N",
    "beta_size"
  )
  lam_min <- 10^par_transformed["log10_lam_min"]
  lam_max <- lam_min + exp(par_transformed["delta_lam"])
  h_down <- 10^par_transformed["log10_h_down"]
  delta_m <- 10^par_transformed["log10_delta_m"]
  m_on <- par_transformed["m_on"]
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
    K_down = 10^par_transformed["log10_K_down"],
    h_down = h_down,
    A_ang = par_transformed["A_ang"],
    m_on = m_on,
    delta_m = delta_m,
    m_off = m_on + delta_m,
    s_on = 10^par_transformed["log10_s_on"],
    s_off = 10^par_transformed["log10_s_off"],
    rho_2N = 10^par_transformed["log10_rho_2N"],
    beta_size = par_transformed["beta_size"],
    c_vol_2N_eff_mm3 = 10^-par_transformed["log10_rho_2N"],
    ratio_4N_2N = 2^par_transformed["beta_size"],
    alpha = 0,
    gamma = 1
  )
}

encode_params <- function(run_params, fit_treatment = TRUE, cfg = NULL) {
  rp <- as.list(run_params)
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
  K_down_v <- need_pos(
    getv(c("K_down"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$K else NULL, 1e12))),
    "K_down"
  )
  h_down_v <- need_pos(
    getv(c("h_down"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$h_down_init else NULL, 1.0))),
    "h_down"
  )
  A_ang_v <- getv(c("A_ang"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_A_ang_default else NULL, 25)))
  if (!is.finite(A_ang_v)) A_ang_v <- as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_A_ang_default else NULL, 25))
  A_ang_v <- clip(A_ang_v, 0, 100)
  m_on_v <- getv(c("m_on"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_m_on_default else NULL, 9.0)))
  if (!is.finite(m_on_v)) m_on_v <- as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_m_on_default else NULL, 9.0))
  m_off_v <- getv(c("m_off"), default = NA_real_)
  delta_m_v <- getv(c("delta_m", "m_off_minus_m_on"), default = NA_real_)
  if ((!is.finite(delta_m_v) || delta_m_v <= 0) && is.finite(m_off_v) && is.finite(m_on_v)) {
    delta_m_v <- m_off_v - m_on_v
  }
  if (!is.finite(delta_m_v) || delta_m_v <= 0) {
    delta_m_v <- as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_delta_m_default else NULL, 1.0))
  }
  delta_m_v <- max(delta_m_v, 1e-3)
  s_on_v <- need_pos(
    getv(c("s_on"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_s_on_default else NULL, 0.3))),
    "s_on"
  )
  s_off_v <- need_pos(
    getv(c("s_off"), default = as.numeric(.first_non_null_local(if (!is.null(cfg)) cfg$o2_s_off_default else NULL, 0.3))),
    "s_off"
  )
  rho_2N_v <- getv(c("rho_2N"), default = NA_real_)
  rho_2N_v <- need_pos(
    if (is.finite(rho_2N_v) && rho_2N_v > 0) rho_2N_v else default_rho_2N_prior_center(cfg),
    "rho_2N"
  )
  beta_size_v <- getv(c("beta_size"), default = default_beta_size_prior_center())
  if (!is.finite(beta_size_v)) stop("Warm-start parameter must be finite: beta_size")

  if (isTRUE(fit_treatment)) {
    alphav <- need_pos(getv(c("alpha")), "alpha")
    gammav <- getv(c("gamma"))
    return(c(
      log10_lam_min = log10(lam_min_v),
      delta_lam = log(lam_gap_v),
      log10_k_o = log10(k_o_v),
      log10_p_misseg = log10(p_misseg_v),
      log10_k_o_mis = log10(k_o_mis_v),
      beta_buffer = beta_buffer_v,
      log10_n_exp = log10(n_exp_v),
      log10_smax = log10(smax_v),
      log10_p_wgd = log10(p_wgd_v),
      log10_K_down = log10(K_down_v),
      log10_h_down = log10(h_down_v),
      A_ang = A_ang_v,
      m_on = m_on_v,
      log10_delta_m = log10(delta_m_v),
      log10_s_on = log10(s_on_v),
      log10_s_off = log10(s_off_v),
      log10_rho_2N = log10(rho_2N_v),
      beta_size = beta_size_v,
      log10_alpha = log10(alphav),
      gamma = gammav
    ))
  }

  c(
    log10_lam_min = log10(lam_min_v),
    delta_lam = log(lam_gap_v),
    log10_k_o = log10(k_o_v),
    log10_p_misseg = log10(p_misseg_v),
    log10_k_o_mis = log10(k_o_mis_v),
    beta_buffer = beta_buffer_v,
    log10_n_exp = log10(n_exp_v),
    log10_smax = log10(smax_v),
    log10_p_wgd = log10(p_wgd_v),
    log10_K_down = log10(K_down_v),
    log10_h_down = log10(h_down_v),
    A_ang = A_ang_v,
    m_on = m_on_v,
    log10_delta_m = log10(delta_m_v),
    log10_s_on = log10(s_on_v),
    log10_s_off = log10(s_off_v),
    log10_rho_2N = log10(rho_2N_v),
    beta_size = beta_size_v
  )
}

read_init_params_t <- function(init_path, bounds, cfg) {
  if (!file.exists(init_path)) stop("init_params_tsv not found: ", init_path)
  tab <- read.delim(init_path, check.names = FALSE, stringsAsFactors = FALSE)
  full_names <- names(bounds$lower)

  out <- NULL
  if (all(c("transformed_parameter", "transformed_value") %in% names(tab))) {
    vals <- setNames(as.numeric(tab$transformed_value), as.character(tab$transformed_parameter))
    missing_names <- setdiff(full_names, names(vals))
    if ("log10_h_down" %in% missing_names) {
      vals[["log10_h_down"]] <- log10(as.numeric(.first_non_null_local(cfg$h_down_init, 1.0)))
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
      out <- encode_params(vals, fit_treatment = cfg$fit_treatment, cfg = cfg)
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

make_bounds <- function(fit_treatment = TRUE,
                        rho_2N_min = 3.2e4, rho_2N_max = 5.6e4,
                        h_down_min = 0.2, h_down_max = 5.0) {
  # Richard-aligned parameterization directly fits p_misseg and k_o_mis.
  rho_2N_min <- as.numeric(rho_2N_min)
  rho_2N_max <- as.numeric(rho_2N_max)
  if (!is.finite(rho_2N_min) || rho_2N_min <= 0) rho_2N_min <- 3.2e4
  if (!is.finite(rho_2N_max) || rho_2N_max <= 0) rho_2N_max <- 5.6e4
  if (rho_2N_min > rho_2N_max) {
    tmp <- rho_2N_min
    rho_2N_min <- rho_2N_max
    rho_2N_max <- tmp
  }
  h_down_min <- as.numeric(h_down_min)
  h_down_max <- as.numeric(h_down_max)
  if (!is.finite(h_down_min) || h_down_min <= 0) h_down_min <- 0.2
  if (!is.finite(h_down_max) || h_down_max <= 0) h_down_max <- 5.0
  if (h_down_min > h_down_max) {
    tmp <- h_down_min
    h_down_min <- h_down_max
    h_down_max <- tmp
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
    log10_K_down = log10(1e4),
    log10_h_down = log10(h_down_min),
    A_ang = 0.0,
    m_on = 6.0,
    log10_delta_m = log10(1e-2),
    log10_s_on = log10(1e-2),
    log10_s_off = log10(1e-2),
    log10_rho_2N = log10(rho_2N_min),
    beta_size = 0.0
  )
  upper <- c(
    log10_lam_min = log10(5),
    delta_lam = log(5),
    log10_k_o = log10(1e4),
    log10_p_misseg = log10(5e-1),
    log10_k_o_mis = log10(1e4),
    beta_buffer = 10.0,
    log10_n_exp = log10(5),
    log10_smax = log10(1),
    log10_p_wgd = log10(1e-1),
    log10_K_down = log10(1e14),
    log10_h_down = log10(h_down_max),
    A_ang = 100.0,
    m_on = 14.0,
    log10_delta_m = log10(4.0),
    log10_s_on = log10(2.0),
    log10_s_off = log10(2.0),
    log10_rho_2N = log10(rho_2N_max),
    beta_size = 2.0
  )

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
      obs_N <- integer(0)
    } else {
      obs_N <- round(as.numeric(obs_pl) * cfg$N_UNIT)
      obs_N <- as.integer(clip(obs_N, cfg$N_MIN, cfg$N_MAX))
      obs_N <- obs_N[is.finite(obs_N)]
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
      ploidy_obs_N = obs_N
    )
    keep[[i]] <- TRUE
  }

  scenarios <- scenarios[keep]
  if (length(scenarios) == 0) stop("No valid scenarios after preprocessing.")
  paired_only <- isTRUE(.first_non_null_local(cfg$paired_only, TRUE))
  n_before_pair_filter <- length(scenarios)
  n_ploidy_before_pair_filter <- sum(vapply(scenarios, function(s) length(s$ploidy_obs_N) > 0, logical(1)))

  if (paired_only) {
    scenarios <- scenarios[vapply(scenarios, function(s) length(s$ploidy_obs_N) > 0, logical(1))]
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

  matched_ploidy <- sum(vapply(scenarios, function(s) length(s$ploidy_obs_N) > 0, logical(1)))
  message(
    "Prepared scenarios: ", length(scenarios),
    " (with terminal ploidy: ", matched_ploidy,
    "; paired_only=", paired_only,
    "; pre_pair_filter_ploidy=", n_ploidy_before_pair_filter, "/", n_before_pair_filter, ")"
  )
  scenarios
}

prepare_cpp_scenarios <- function(scenarios, cfg) {
  n <- length(scenarios)
  cohort_code <- integer(n)
  dose_vec <- numeric(n)
  treat_day_vec <- numeric(n)
  obs_steps_list <- vector("list", n)
  sim_end_step_vec <- integer(n)
  obs_burden_list <- vector("list", n)
  keep_burden_list <- vector("list", n)
  ploidy_idx_list <- vector("list", n)
  ploidy_count_list <- vector("list", n)

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

    if (length(sc$ploidy_obs_N) > 0) {
      obs_tab <- table(as.character(sc$ploidy_obs_N))
      idx <- as.integer(names(obs_tab)) - as.integer(cfg$N_MIN) + 1L
      cnt <- as.numeric(obs_tab)
      keep_idx <- is.finite(idx) & is.finite(cnt) & (cnt > 0) & (idx >= 1L) & (idx <= (cfg$N_MAX - cfg$N_MIN + 1L))
      ploidy_idx_list[[i]] <- as.integer(idx[keep_idx])
      ploidy_count_list[[i]] <- as.numeric(cnt[keep_idx])
    } else {
      ploidy_idx_list[[i]] <- integer(0)
      ploidy_count_list[[i]] <- numeric(0)
    }
  }

  list(
    cohort_code = cohort_code,
    dose = dose_vec,
    treat_day = treat_day_vec,
    obs_steps = obs_steps_list,
    sim_end_step = sim_end_step_vec,
    obs_burden = obs_burden_list,
    keep_burden = keep_burden_list,
    ploidy_idx = ploidy_idx_list,
    ploidy_count = ploidy_count_list
  )
}

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
    total_size = cfg$init_total_size
  )
  init_state_4N <- make_init_state(
    grid_pre = grid_pre,
    grid_post = grid_post,
    ploidy = 4,
    layer = "post",
    N_UNIT = cfg$N_UNIT,
    total_size = cfg$init_total_size
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

simulate_one <- function(run_params, scenario, cfg, model_core = NULL) {
  if (is.null(model_core)) {
    if (!is.null(cfg$model_core)) {
      model_core <- cfg$model_core
    } else {
      model_core <- build_model_core(run_params, cfg)
    }
  }

  grid_pre <- model_core$grid_pre
  R0 <- model_core$R0
  R1 <- model_core$R1
  init_state <- if (scenario$cohort == "2N") model_core$init_state_2N else model_core$init_state_4N

  obs_steps <- as.integer(round(scenario$obs_days / cfg$DT))
  sim_end_step <- as.integer(round(scenario$sim_end_day / cfg$DT))
  vol_by_N <- cell_volume_mm3_by_N(grid_pre, run_params = run_params, cfg = cfg)

  if (!exists("cpp_o2invivo_simulate_one", mode = "function", inherits = TRUE)) {
    stop("Required C++ function missing: cpp_o2invivo_simulate_one")
  }

  h_down_use <- as.numeric(.first_non_null_local(run_params$h_down, cfg$h_down_init, 1.0))
  if (!is.finite(h_down_use) || h_down_use <= 0) h_down_use <- 1.0
  K_down_use <- as.numeric(.first_non_null_local(run_params$K_down, cfg$K, 1e12))
  if (!is.finite(K_down_use) || K_down_use <= 0) K_down_use <- 1e12
  m_on_use <- as.numeric(.first_non_null_local(run_params$m_on, cfg$o2_m_on_default, 9.0))
  if (!is.finite(m_on_use)) m_on_use <- 9.0
  delta_m_use <- as.numeric(.first_non_null_local(run_params$delta_m, cfg$o2_delta_m_default, 1.0))
  if (!is.finite(delta_m_use) || delta_m_use <= 0) delta_m_use <- 1.0
  m_off_use <- as.numeric(.first_non_null_local(run_params$m_off, m_on_use + delta_m_use))
  if (!is.finite(m_off_use) || m_off_use <= m_on_use) m_off_use <- m_on_use + delta_m_use
  s_on_use <- as.numeric(.first_non_null_local(run_params$s_on, cfg$o2_s_on_default, 0.3))
  if (!is.finite(s_on_use) || s_on_use <= 0) s_on_use <- 0.3
  s_off_use <- as.numeric(.first_non_null_local(run_params$s_off, cfg$o2_s_off_default, 0.3))
  if (!is.finite(s_off_use) || s_off_use <= 0) s_off_use <- 0.3
  p_wgd_use <- as.numeric(.first_non_null_local(run_params$p_wgd, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
  boundary_mode <- as.character(.first_non_null_local(run_params$boundary, "drop"))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 0)

  sim <- cpp_o2invivo_simulate_one(
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
    O2_base = as.numeric(.first_non_null_local(cfg$O2_fixed, 5.0)),
    o2_feedback = isTRUE(.first_non_null_local(cfg$o2_burden_feedback, TRUE)),
    o2_min = as.numeric(.first_non_null_local(cfg$o2_min, 0.0)),
    h_O2 = as.numeric(h_down_use),
    K_down = as.numeric(K_down_use),
    A_ang = as.numeric(.first_non_null_local(run_params$A_ang, cfg$o2_A_ang_default, 0.0)),
    m_on = as.numeric(m_on_use),
    m_off = as.numeric(m_off_use),
    s_on = as.numeric(s_on_use),
    s_off = as.numeric(s_off_use),
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

evaluate_objective_components_raw <- function(par_transformed, scenarios, cfg) {
  cfg_eval <- cfg
  rp <- decode_params(
    par_transformed,
    fit_treatment = cfg_eval$fit_treatment
  )
  model_core <- if (!is.null(cfg_eval$model_core)) cfg_eval$model_core else build_model_core(cfg = cfg_eval)
  scenario_cpp <- if (!is.null(cfg_eval$scenario_cpp)) cfg_eval$scenario_cpp else prepare_cpp_scenarios(scenarios, cfg_eval)
  vol_by_N <- cell_volume_mm3_by_N(model_core$grid_pre, run_params = rp, cfg = cfg_eval)

  h_down_use <- as.numeric(.first_non_null_local(rp$h_down, cfg_eval$h_down_init, 1.0))
  if (!is.finite(h_down_use) || h_down_use <= 0) h_down_use <- 1.0
  K_down_use <- as.numeric(.first_non_null_local(rp$K_down, cfg_eval$K, 1e12))
  if (!is.finite(K_down_use) || K_down_use <= 0) K_down_use <- 1e12
  m_on_use <- as.numeric(.first_non_null_local(rp$m_on, cfg_eval$o2_m_on_default, 9.0))
  if (!is.finite(m_on_use)) m_on_use <- 9.0
  delta_m_use <- as.numeric(.first_non_null_local(rp$delta_m, cfg_eval$o2_delta_m_default, 1.0))
  if (!is.finite(delta_m_use) || delta_m_use <= 0) delta_m_use <- 1.0
  m_off_use <- as.numeric(.first_non_null_local(rp$m_off, m_on_use + delta_m_use))
  if (!is.finite(m_off_use) || m_off_use <= m_on_use) m_off_use <- m_on_use + delta_m_use
  s_on_use <- as.numeric(.first_non_null_local(rp$s_on, cfg_eval$o2_s_on_default, 0.3))
  if (!is.finite(s_on_use) || s_on_use <= 0) s_on_use <- 0.3
  s_off_use <- as.numeric(.first_non_null_local(rp$s_off, cfg_eval$o2_s_off_default, 0.3))
  if (!is.finite(s_off_use) || s_off_use <= 0) s_off_use <- 0.3
  p_wgd_use <- as.numeric(.first_non_null_local(rp$p_wgd, 0.0))
  if (!is.finite(p_wgd_use)) p_wgd_use <- 0.0
  boundary_mode <- as.character(.first_non_null_local(rp$boundary, "drop"))
  burden_floor <- pmax(as.numeric(.first_non_null_local(cfg_eval$burden_log_eps, 1e-12)), 0)

  comp <- cpp_o2invivo_objective_components(
    cohort_code = as.integer(scenario_cpp$cohort_code),
    dose_vec = as.numeric(scenario_cpp$dose),
    treat_day_vec = as.numeric(scenario_cpp$treat_day),
    obs_steps_list = scenario_cpp$obs_steps,
    sim_end_step_vec = as.integer(scenario_cpp$sim_end_step),
    obs_burden_list = scenario_cpp$obs_burden,
    keep_burden_list = scenario_cpp$keep_burden,
    ploidy_idx_list = scenario_cpp$ploidy_idx,
    ploidy_count_list = scenario_cpp$ploidy_count,
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
    O2_base = as.numeric(.first_non_null_local(cfg_eval$O2_fixed, 5.0)),
    o2_feedback = isTRUE(.first_non_null_local(cfg_eval$o2_burden_feedback, TRUE)),
    o2_min = as.numeric(.first_non_null_local(cfg_eval$o2_min, 0.0)),
    h_O2 = as.numeric(h_down_use),
    K_down = as.numeric(K_down_use),
    A_ang = as.numeric(.first_non_null_local(rp$A_ang, cfg_eval$o2_A_ang_default, 0.0)),
    m_on = as.numeric(m_on_use),
    m_off = as.numeric(m_off_use),
    s_on = as.numeric(s_on_use),
    s_off = as.numeric(s_off_use),
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
    vol_by_N = as.numeric(vol_by_N),
    burden_floor = as.numeric(burden_floor),
    burden_log_eps = as.numeric(.first_non_null_local(cfg_eval$burden_log_eps, 1e-12)),
    huber_k_burden_log = as.numeric(.first_non_null_local(cfg_eval$huber_k_burden_log, cfg_eval$huber_k, 0.1)),
    eps_prob = as.numeric(cfg_eval$eps_prob),
    agg_burden = as.character(normalize_agg_method(cfg_eval$scenario_agg_burden, default = "huber")),
    agg_ploidy = as.character(normalize_agg_method(cfg_eval$scenario_agg_ploidy, default = "huber")),
    scenario_weight_burden = isTRUE(.first_non_null_local(cfg_eval$scenario_weight_burden, TRUE)),
    scenario_weight_ploidy = isTRUE(.first_non_null_local(cfg_eval$scenario_weight_ploidy, TRUE)),
    scenario_agg_huber_k = as.numeric(.first_non_null_local(cfg_eval$scenario_agg_huber_k, 1.5))
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
  agg_b <- normalize_agg_method(cfg_eval$scenario_agg_burden, default = "huber")
  agg_p <- normalize_agg_method(cfg_eval$scenario_agg_ploidy, default = "huber")
  list(
    L_b = L_b,
    L_p = L_p,
    n_burden = as.integer(.first_non_null_local(comp$n_burden, 0L)),
    n_ploidy = as.integer(.first_non_null_local(comp$n_ploidy, 0L)),
    cache_g_build = cache_g_build,
    cache_g_hit = cache_g_hit,
    cache_g_hysteresis = cache_g_hysteresis,
    agg_burden = agg_b,
    agg_ploidy = agg_p
  )
}

evaluate_objective_components <- function(par_transformed, scenarios, cfg) {
  raw <- evaluate_objective_components_raw(par_transformed, scenarios = scenarios, cfg = cfg)
  L_b <- raw$L_b
  L_p <- raw$L_p

  scale_b <- if (isTRUE(cfg$loss_rescale)) cfg$loss_scale_burden else 1.0
  scale_p <- if (isTRUE(cfg$loss_rescale)) cfg$loss_scale_ploidy else 1.0
  if (!is.finite(scale_b) || scale_b <= 0) scale_b <- 1.0
  if (!is.finite(scale_p) || scale_p <= 0) scale_p <- 1.0

  L_b_scaled <- L_b / scale_b
  L_p_scaled <- L_p / scale_p
  L_data <- cfg$w_burden * L_b_scaled + cfg$w_ploidy * L_p_scaled
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
      "L=%.6f (data=%.6f, prior=%.6f; burden=%.6f, ploidy=%.6f; scaled burden=%.6f, scaled ploidy=%.6f)\n",
      L, L_data, L_prior, L_b, L_p, L_b_scaled, L_p_scaled
    ))
  }
  list(
    L = L,
    L_data = L_data,
    L_prior = L_prior,
    L_prior_raw = L_prior_raw,
    L_b = L_b,
    L_p = L_p,
    L_b_scaled = L_b_scaled,
    L_p_scaled = L_p_scaled,
    scale_b = scale_b,
    scale_p = scale_p,
    n_burden = raw$n_burden,
    n_ploidy = raw$n_ploidy,
    cache_g_build = raw$cache_g_build,
    cache_g_hit = raw$cache_g_hit,
    cache_g_hysteresis = raw$cache_g_hysteresis,
    agg_burden = raw$agg_burden,
    agg_ploidy = raw$agg_ploidy,
    n_prior_terms = as.integer(.first_non_null_local(prior$n_terms, 0))
  )
}

evaluate_objective <- function(par_transformed, scenarios, cfg) {
  evaluate_objective_components(par_transformed, scenarios = scenarios, cfg = cfg)$L
}

run_optimizer <- function(objective_fn, lower, upper, cfg, argv, stage_label = "fit", init_par = NULL) {
  n_cores <- as.integer(max(1L, ifelse(is.finite(cfg$n_cores), cfg$n_cores, 1L)))
  use_deoptim <- isTRUE(cfg$use_deoptim)
  deoptim_parallel <- isTRUE(cfg$deoptim_parallel)
  fmt_int <- function(x) {
    if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
    as.character(as.integer(x))
  }
  cluster_workers <- function(cl) {
    n <- tryCatch(length(cl), error = function(e) NA_integer_)
    if (is.na(n) || !is.finite(n)) return(NA_integer_)
    as.integer(max(0L, n))
  }
  log_parallel_backend <- function(solver, requested, started, active, fallback, reason = NULL, extra = NULL) {
    msg <- paste0(
      "[", stage_label, "] ", solver, " backend: requested_workers=", fmt_int(requested),
      ", started_workers=", fmt_int(started),
      ", active_workers=", fmt_int(active),
      ", fallback=", if (isTRUE(fallback)) "TRUE" else "FALSE"
    )
    if (isTRUE(fallback) && !is.null(reason) && nzchar(reason)) {
      msg <- paste0(msg, ", reason=", reason)
    }
    if (!is.null(extra) && nzchar(extra)) {
      msg <- paste0(msg, ", ", extra)
    }
    message(msg)
  }
  init_cluster_workers <- function(cl, objective_fn, cfg, stage_label) {
    if (is.null(cfg$cpp_wrapper_path) || !nzchar(cfg$cpp_wrapper_path) || !file.exists(cfg$cpp_wrapper_path)) {
      stop("[", stage_label, "] Missing sourceCpp wrapper path for worker initialization.")
    }
    parallel::clusterCall(
      cl,
      function(wrapper_path) {
        source(wrapper_path, local = .GlobalEnv)
        NULL
      },
      cfg$cpp_wrapper_path
    )
    parallel::clusterCall(cl, function() {
      required <- c("cpp_o2invivo_build_G_for_o2_triplet", "cpp_o2invivo_simulate_one", "cpp_o2invivo_objective_components")
      missing <- required[!vapply(required, exists, logical(1), mode = "function", inherits = TRUE)]
      if (length(missing) > 0L) {
        stop("Worker missing required C++ wrapper functions: ", paste(missing, collapse = ", "))
      }
      NULL
    })
    parallel::clusterCall(cl, function() {
      Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
      NULL
    })
    if (!is.null(cfg$model_path) && nzchar(cfg$model_path) && file.exists(cfg$model_path)) {
      parallel::clusterCall(cl, function(path) {
        Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(path))
        NULL
      }, cfg$model_path)
    }
    export_global <- c(
      "evaluate_objective",
      "evaluate_objective_components",
      "evaluate_objective_components_raw",
      "build_model_core",
      "prepare_cpp_scenarios",
      "simulate_one",
      "make_init_state",
      "decode_params",
      "huber_mean",
      "weighted_mean_safe",
      "weighted_median_safe",
      "robust_huber_location",
      "normalize_agg_method",
      "aggregate_scenario_losses",
      "compute_soft_prior_penalty",
      "as_num",
      "clip",
      ".first_non_null_local",
      "default_beta_size_prior_center",
      "default_rho_2N_prior_bounds",
      "default_rho_2N_prior_center",
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
    de_fallback <- FALSE
    de_fallback_reason <- ""
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
      tryCatch(
        init_cluster_workers(cl, objective_fn = objective_fn, cfg = cfg, stage_label = stage_label),
        error = function(e) {
          stop("[", stage_label, "] Failed to initialize DEoptim workers: ", conditionMessage(e))
        }
      )
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
      fallback = de_fallback,
      reason = de_fallback_reason,
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
      active_workers = de_active,
      fallback = de_fallback,
      fallback_reason = if (isTRUE(de_fallback)) de_fallback_reason else ""
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
    opt_fallback <- FALSE
    opt_fallback_reason <- ""
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
      fallback = opt_fallback,
      reason = opt_fallback_reason,
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
        fallback = opt_fallback,
        fallback_reason = if (isTRUE(opt_fallback)) opt_fallback_reason else "",
        n_starts = length(starts)
      )
    )
  }

  best_par <- as.numeric(best_par)
  names(best_par) <- names(lower)
  list(best_par = best_par, optim_res = optim_res)
}

run_subset_fit <- function(vary_names, base_par, bounds, scenarios, cfg, argv, stage_label = "fit", init_par = NULL) {
  stopifnot(all(vary_names %in% names(base_par)))
  lower_sub <- bounds$lower[vary_names]
  upper_sub <- bounds$upper[vary_names]
  init_sub <- NULL
  if (!is.null(init_par)) {
    init_sub <- init_par[vary_names]
  }
  objective_subset <- function(par_sub) {
    full <- base_par
    full[vary_names] <- par_sub
    evaluate_objective(full, scenarios = scenarios, cfg = cfg)
  }
  opt <- run_optimizer(
    objective_fn = objective_subset,
    lower = lower_sub,
    upper = upper_sub,
    cfg = cfg,
    argv = argv,
    stage_label = stage_label,
    init_par = init_sub
  )
  full_best <- base_par
  full_best[vary_names] <- opt$best_par[vary_names]
  list(best_par = full_best, optim_res = opt$optim_res)
}

collect_predictions <- function(run_params, scenarios, cfg) {
  model_core <- if (!is.null(cfg$model_core)) cfg$model_core else build_model_core(cfg = cfg)
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
    if (length(sc$ploidy_obs_N) > 0) {
      obs_tab <- table(sc$ploidy_obs_N)
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
    n_cores = cfg$n_cores,
    label = "predict",
    fn = pred_one
  )

  list(
    burden = bind_rows(lapply(pred_rows, function(x) x$burden)),
    ploidy = bind_rows(Filter(Negate(is.null), lapply(pred_rows, function(x) x$ploidy)))
  )
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  script_dir <- get_script_dir()

  model_path <- file.path(script_dir, "model_O2_invivo.R")
  if (!file.exists(model_path)) stop("Cannot find model_O2_invivo.R at ", model_path)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = script_dir)
  source(model_path)
  required_cpp_fit <- c("cpp_o2invivo_build_G_for_o2_triplet", "cpp_o2invivo_simulate_one", "cpp_o2invivo_objective_components")
  missing_cpp_fit <- required_cpp_fit[!vapply(required_cpp_fit, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing_cpp_fit) > 0L) {
    stop("Required C++ symbols missing for fit path: ", paste(missing_cpp_fit, collapse = ", "))
  }
  cpp_dll <- tryCatch(
    o2invivo_cpp_dll_info(),
    error = function(e) stop("Failed to resolve compiled O2_invivo DLL info: ", conditionMessage(e))
  )

  default_data_dir <- normalizePath(file.path(script_dir, "..", "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
  data_dir <- if (!is.null(argv$data_dir)) argv$data_dir else default_data_dir
  truncate_at_treatment <- as_bool(argv$truncate_at_treatment, FALSE)
  w_burden_vec <- as_num_vec(argv$w_burden, 1.0)
  w_ploidy_vec <- as_num_vec(argv$w_ploidy, 1.0)
  weight_schedule <- make_weight_schedule(w_burden_vec, w_ploidy_vec)
  n_cores_arg <- as_int(argv$n_cores, NA_integer_)
  n_cores_use <- if (is.finite(n_cores_arg)) n_cores_arg else default_n_cores()

  cfg <- list(
    # model constants
    model_path = model_path,
    cpp_dll_name = as.character(cpp_dll$name),
    cpp_dll_path = as.character(cpp_dll$path),
    cpp_wrapper_path = as.character(cpp_dll$wrapper_path),
    N_UNIT = as_int(argv$N_UNIT, 22L),
    N_MIN = as_int(argv$N_MIN, 22L),
    N_MAX = as_int(argv$N_MAX, 154L),
    DT = as_num(argv$dt, 0.5),
    O2_fixed = as_num(argv$O2, 5.0),
    o2_burden_feedback = as_bool(argv$o2_burden_feedback, TRUE),
    o2_logN_eps = as_num(argv$o2_logn_eps, 1.0),
    o2_cache_bin_pct = as_num(argv$o2_cache_bin_pct, 0.01),
    o2_cache_hysteresis_pct = as_num(argv$o2_cache_hysteresis_pct, 0.005),
    o2_cache_profile = as_bool(argv$o2_cache_profile, FALSE),
    o2_A_ang_default = as_num(argv$o2_A_ang_default, 25.0),
    o2_m_on_default = as_num(argv$o2_m_on_default, 9.0),
    o2_delta_m_default = as_num(argv$o2_delta_m_default, 1.0),
    o2_s_on_default = as_num(argv$o2_s_on_default, 0.3),
    o2_s_off_default = as_num(argv$o2_s_off_default, 0.3),
    o2_min = as_num(argv$o2_min, 0.0),
    h_down_init = as_num(argv$h_down_init, 1.0),
    h_down_min = as_num(argv$h_down_min, 0.2),
    h_down_max = as_num(argv$h_down_max, 5.0),
    K = as_num(argv$K, 1e12),
    crowding = if (!is.null(argv$crowding)) argv$crowding else "logistic",
    init_total_size = as_num(argv$init_total_size, 1e6),
    dose_ref = as_num(argv$dose_ref, 30),
    tx_mult_min = as_num(argv$tx_mult_min, 0.05),
    min_pop = as_num(argv$min_pop, 1e-12),
    # objective settings
    huber_k = as_num(argv$huber_k, 0.1),
    huber_k_burden_log = as_num(argv$huber_k_burden_log, as_num(argv$huber_k, 0.1)),
    burden_log_eps = as_num(argv$burden_log_eps, 1e-12),
    burden_exclude_day0 = as_bool(argv$burden_exclude_day0, TRUE),
    scenario_agg_burden = normalize_agg_method(.first_non_null_local(argv$scenario_agg_burden, argv$scenario_agg, "huber"), default = "huber"),
    scenario_agg_ploidy = normalize_agg_method(.first_non_null_local(argv$scenario_agg_ploidy, argv$scenario_agg, "huber"), default = "huber"),
    scenario_agg_huber_k = as_num(argv$scenario_agg_huber_k, 1.5),
    scenario_weight_burden = as_bool(argv$scenario_weight_burden, TRUE),
    scenario_weight_ploidy = as_bool(argv$scenario_weight_ploidy, TRUE),
    rho_2N_min = as_num(argv$rho_2N_min, 3.2e4),
    rho_2N_max = as_num(argv$rho_2N_max, 5.6e4),
    use_soft_prior = as_bool(argv$use_soft_prior, TRUE),
    lambda_prior = as_num(argv$lambda_prior, 0.1),
    prior_center_log10_k_o = as_num(argv$prior_center_log10_k_o, log10(50)),
    prior_sd_log10_k_o = as_num(argv$prior_sd_log10_k_o, 1.0),
    prior_center_log10_K_down = as_num(argv$prior_center_log10_K_down, log10(1e12)),
    prior_sd_log10_K_down = as_num(argv$prior_sd_log10_K_down, 1.0),
    prior_center_log10_h_down = as_num(argv$prior_center_log10_h_down, log10(as_num(argv$h_down_init, 1.0))),
    prior_sd_log10_h_down = as_num(argv$prior_sd_log10_h_down, 0.5),
    prior_center_beta_size = as_num(argv$prior_center_beta_size, default_beta_size_prior_center()),
    prior_sd_beta_size = as_num(argv$prior_sd_beta_size, 0.5),
    prior_center_log10_n_exp = as_num(argv$prior_center_log10_n_exp, 0.0),
    prior_sd_log10_n_exp = as_num(argv$prior_sd_log10_n_exp, 0.5),
    prior_center_log10_rho_2N = as_num(argv$prior_center_log10_rho_2N, log10(sqrt(as_num(argv$rho_2N_min, 3.2e4) * as_num(argv$rho_2N_max, 5.6e4)))),
    prior_sd_log10_rho_2N = as_num(argv$prior_sd_log10_rho_2N, 0.35),
    w_burden = weight_schedule$w_burden[[length(weight_schedule$w_burden)]],
    w_ploidy = weight_schedule$w_ploidy[[length(weight_schedule$w_ploidy)]],
    w_burden_schedule = weight_schedule$w_burden,
    w_ploidy_schedule = weight_schedule$w_ploidy,
    n_weight_passes = weight_schedule$n_pass,
    loss_rescale = as_bool(argv$loss_rescale, FALSE),
    loss_scale_burden = as_num(argv$loss_scale_burden, NA_real_),
    loss_scale_ploidy = as_num(argv$loss_scale_ploidy, NA_real_),
    loss_scale_eps = as_num(argv$loss_scale_eps, 1e-8),
    loss_scale_source = "unset",
    optim_trace = as_bool(argv$optim_trace, TRUE),
    optim_trace_every = as_int(argv$optim_trace_every, 1L),
    eps_prob = as_num(argv$eps_prob, 1e-12),
    trace_obj = as_bool(argv$trace_obj, FALSE),
    # fitting options
    fit_treatment = as_bool(argv$fit_treatment, FALSE),
    dose_zero_only = as_bool(argv$dose_zero_only, TRUE),
    paired_only = as_bool(argv$paired_only, TRUE),
    truncate_at_treatment = truncate_at_treatment,
    ploidy_at_harvest = as_bool(argv$ploidy_at_harvest, TRUE),
    two_stage = as_bool(argv$two_stage, FALSE),
    stage1_w_burden = as_num(argv$stage1_w_burden, 1.0),
    stage1_w_ploidy = as_num(argv$stage1_w_ploidy, 0.0),
    stage2_w_burden = as_num(argv$stage2_w_burden, 0.0),
    stage2_w_ploidy = as_num(argv$stage2_w_ploidy, 1.0),
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
  if (!is.finite(cfg$o2_A_ang_default) || cfg$o2_A_ang_default < 0 || cfg$o2_A_ang_default > 100) stop("o2_A_ang_default must be in [0,100]")
  if (!is.finite(cfg$o2_delta_m_default) || cfg$o2_delta_m_default <= 0) stop("o2_delta_m_default must be > 0")
  if (!is.finite(cfg$o2_s_on_default) || cfg$o2_s_on_default <= 0) stop("o2_s_on_default must be > 0")
  if (!is.finite(cfg$o2_s_off_default) || cfg$o2_s_off_default <= 0) stop("o2_s_off_default must be > 0")
  if (!is.finite(cfg$O2_fixed) || cfg$O2_fixed < 0 || cfg$O2_fixed > 100) {
    stop("O2 must be in percent scale [0, 100].")
  }
  if (!is.finite(cfg$o2_min) || cfg$o2_min < 0 || cfg$o2_min > 100) {
    stop("o2_min must be in percent scale [0, 100].")
  }
  if (!is.finite(cfg$h_down_init) || cfg$h_down_init <= 0) stop("h_down_init must be > 0")
  if (!is.finite(cfg$h_down_min) || cfg$h_down_min <= 0) stop("h_down_min must be > 0")
  if (!is.finite(cfg$h_down_max) || cfg$h_down_max <= 0) stop("h_down_max must be > 0")
  if (cfg$h_down_max < cfg$h_down_min) stop("h_down_max must be >= h_down_min")
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
  if (!is.finite(cfg$loss_scale_eps) || cfg$loss_scale_eps <= 0) stop("loss_scale_eps must be > 0")
  if (!is.finite(cfg$burden_log_eps) || cfg$burden_log_eps <= 0) stop("burden_log_eps must be > 0")
  if (!cfg$scenario_agg_burden %in% c("mean", "median", "huber")) stop("scenario_agg_burden must be one of: mean, median, huber")
  if (!cfg$scenario_agg_ploidy %in% c("mean", "median", "huber")) stop("scenario_agg_ploidy must be one of: mean, median, huber")
  if (!is.finite(cfg$scenario_agg_huber_k) || cfg$scenario_agg_huber_k <= 0) stop("scenario_agg_huber_k must be > 0")
  if (!is.finite(cfg$rho_2N_min) || cfg$rho_2N_min <= 0) stop("rho_2N_min must be > 0")
  if (!is.finite(cfg$rho_2N_max) || cfg$rho_2N_max <= 0) stop("rho_2N_max must be > 0")
  if (cfg$rho_2N_max < cfg$rho_2N_min) stop("rho_2N_max must be >= rho_2N_min")
  if (!is.finite(cfg$lambda_prior) || cfg$lambda_prior < 0) stop("lambda_prior must be >= 0")
  if (!is.finite(cfg$prior_sd_log10_k_o) || cfg$prior_sd_log10_k_o <= 0) stop("prior_sd_log10_k_o must be > 0")
  if (!is.finite(cfg$prior_sd_log10_K_down) || cfg$prior_sd_log10_K_down <= 0) stop("prior_sd_log10_K_down must be > 0")
  if (!is.finite(cfg$prior_sd_log10_h_down) || cfg$prior_sd_log10_h_down <= 0) stop("prior_sd_log10_h_down must be > 0")
  if (!is.finite(cfg$prior_sd_beta_size) || cfg$prior_sd_beta_size <= 0) stop("prior_sd_beta_size must be > 0")
  if (!is.finite(cfg$prior_sd_log10_n_exp) || cfg$prior_sd_log10_n_exp <= 0) stop("prior_sd_log10_n_exp must be > 0")
  if (!is.finite(cfg$prior_sd_log10_rho_2N) || cfg$prior_sd_log10_rho_2N <= 0) stop("prior_sd_log10_rho_2N must be > 0")
  if (!is.na(cfg$loss_scale_burden) && (!is.finite(cfg$loss_scale_burden) || cfg$loss_scale_burden <= 0)) stop("loss_scale_burden must be > 0")
  if (!is.na(cfg$loss_scale_ploidy) && (!is.finite(cfg$loss_scale_ploidy) || cfg$loss_scale_ploidy <= 0)) stop("loss_scale_ploidy must be > 0")
  if (!cfg$use_deoptim && cfg$deoptim_parallel) stop("deoptim_parallel=TRUE requires use_deoptim=TRUE")
  if (isTRUE(cfg$two_stage) && cfg$n_weight_passes > 1L) {
    stop("Weight arrays for w_burden/w_ploidy are supported only when two_stage=FALSE.")
  }

  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- file.path(data_dir, "all_ploidy.tsv")
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)
  cfg$model_core <- build_model_core(cfg = cfg)
  cfg$scenario_cpp <- prepare_cpp_scenarios(scenarios, cfg)

  bounds <- make_bounds(
    fit_treatment = cfg$fit_treatment,
    rho_2N_min = cfg$rho_2N_min,
    rho_2N_max = cfg$rho_2N_max,
    h_down_min = cfg$h_down_min,
    h_down_max = cfg$h_down_max
  )
  full_names <- names(bounds$lower)
  default_par_t <- (bounds$lower + bounds$upper) / 2
  names(default_par_t) <- full_names
  if ("log10_K_down" %in% full_names) default_par_t[["log10_K_down"]] <- log10(as.numeric(.first_non_null_local(cfg$K, 1e12)))
  if ("log10_h_down" %in% full_names) default_par_t[["log10_h_down"]] <- log10(as.numeric(.first_non_null_local(cfg$h_down_init, 1.0)))
  if ("A_ang" %in% full_names) default_par_t[["A_ang"]] <- clip(cfg$o2_A_ang_default, 0, 100)
  if ("m_on" %in% full_names) default_par_t[["m_on"]] <- cfg$o2_m_on_default
  if ("log10_delta_m" %in% full_names) default_par_t[["log10_delta_m"]] <- log10(cfg$o2_delta_m_default)
  if ("log10_s_on" %in% full_names) default_par_t[["log10_s_on"]] <- log10(cfg$o2_s_on_default)
  if ("log10_s_off" %in% full_names) default_par_t[["log10_s_off"]] <- log10(cfg$o2_s_off_default)
  if ("log10_rho_2N" %in% full_names) default_par_t[["log10_rho_2N"]] <- log10(default_rho_2N_prior_center(cfg))
  if ("beta_size" %in% full_names) default_par_t[["beta_size"]] <- default_beta_size_prior_center()
  init_params_tsv <- if (!is.null(argv$init_params_tsv)) argv$init_params_tsv else NULL
  warm_start_t <- if (!is.null(init_params_tsv)) {
    read_init_params_t(init_params_tsv, bounds = bounds, cfg = cfg)
  } else {
    NULL
  }
  cfg <- resolve_loss_scales(cfg, default_par_t = default_par_t, warm_start_t = warm_start_t, scenarios = scenarios)
  if (isTRUE(cfg$loss_rescale)) {
    message(
      "Loss rescaling enabled. scale_burden=", signif(cfg$loss_scale_burden, 6),
      ", scale_ploidy=", signif(cfg$loss_scale_ploidy, 6),
      " (source=", cfg$loss_scale_source, ")"
    )
  }
  message(
    "Scenario aggregation: burden=", cfg$scenario_agg_burden,
    " (weighted=", if (isTRUE(cfg$scenario_weight_burden)) "TRUE" else "FALSE", "), ",
    "ploidy=", cfg$scenario_agg_ploidy,
    " (weighted=", if (isTRUE(cfg$scenario_weight_ploidy)) "TRUE" else "FALSE", "), ",
    "huber_k=", signif(cfg$scenario_agg_huber_k, 6)
  )
  if (isTRUE(cfg$use_soft_prior) && cfg$lambda_prior > 0) {
    message(
      "Soft prior enabled: lambda_prior=", signif(cfg$lambda_prior, 6),
      "; centers(log10_k_o, log10_K_down, log10_h_down, beta_size, log10_n_exp, log10_rho_2N)=(",
      signif(cfg$prior_center_log10_k_o, 6), ", ",
      signif(cfg$prior_center_log10_K_down, 6), ", ",
      signif(cfg$prior_center_log10_h_down, 6), ", ",
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
    "Burden observation model enabled: log-volume Huber on V(mm^3), ",
    "V_pred = sum_n n_n * [(1/rho_2N) * (P/2)^beta_size], ",
    "rho_2N_range=[", signif(cfg$rho_2N_min, 6), ", ", signif(cfg$rho_2N_max, 6), "] cells/mm^3",
    "; burden_exclude_day0=", if (isTRUE(cfg$burden_exclude_day0)) "TRUE" else "FALSE"
  )
  message(
    "O2 mode: ",
    if (isTRUE(cfg$o2_burden_feedback)) "dynamic feedback" else "fixed",
    "; O2_base(%)=", signif(cfg$O2_fixed, 6),
    ", o2_min=", signif(cfg$o2_min, 6),
    ", h_down_init=", signif(cfg$h_down_init, 6),
    if (isTRUE(cfg$o2_burden_feedback)) {
      "; window params fitted: K_down,h_down,A_ang,m_on,delta_m,s_on,s_off"
    } else {
      "; O2 window inactive."
    }
  )
  message(
    "O2 cache: bin_pct=", signif(cfg$o2_cache_bin_pct, 6),
    ", hysteresis_pct=", signif(cfg$o2_cache_hysteresis_pct, 6),
    ", profile=", if (isTRUE(cfg$o2_cache_profile)) "TRUE" else "FALSE"
  )
  set.seed(cfg$seed)
  stage1_comp <- NULL
  stage2_comp <- NULL
  stage1_best_par_t <- NULL
  stage1_vary <- character(0)
  stage2_vary <- character(0)
  single_pass_log <- NULL

  if (isTRUE(cfg$two_stage)) {
    stage1_vary <- intersect(full_names, c("log10_lam_min", "delta_lam", "log10_k_o"))
    stage2_vary <- setdiff(full_names, stage1_vary)
    if (length(stage1_vary) == 0 || length(stage2_vary) == 0) {
      stop("two_stage mode requires both growth and ploidy parameter sets to be present.")
    }

    stage1_cfg <- cfg
    stage1_cfg$w_burden <- cfg$stage1_w_burden
    stage1_cfg$w_ploidy <- cfg$stage1_w_ploidy
    stage2_cfg <- cfg
    stage2_cfg$w_burden <- cfg$stage2_w_burden
    stage2_cfg$w_ploidy <- cfg$stage2_w_ploidy

    message(
      "Two-stage fit enabled. Stage1(growth): ", paste(stage1_vary, collapse = ", "),
      "; Stage2(ploidy): ", paste(stage2_vary, collapse = ", ")
    )
    message(
      "Stage1 weights: w_burden=", stage1_cfg$w_burden,
      ", w_ploidy=", stage1_cfg$w_ploidy,
      "; Stage2 weights: w_burden=", stage2_cfg$w_burden,
      ", w_ploidy=", stage2_cfg$w_ploidy
    )

    stage1_fit <- run_subset_fit(
      vary_names = stage1_vary,
      base_par = default_par_t,
      bounds = bounds,
      scenarios = scenarios,
      cfg = stage1_cfg,
      argv = argv,
      stage_label = "stage1_growth",
      init_par = warm_start_t
    )
    stage1_best_par_t <- stage1_fit$best_par
    stage1_comp <- evaluate_objective_components(stage1_best_par_t, scenarios = scenarios, cfg = stage1_cfg)

    stage2_fit <- run_subset_fit(
      vary_names = stage2_vary,
      base_par = stage1_best_par_t,
      bounds = bounds,
      scenarios = scenarios,
      cfg = stage2_cfg,
      argv = argv,
      stage_label = "stage2_ploidy",
      init_par = warm_start_t
    )
    best_par_t <- stage2_fit$best_par
    stage2_comp <- evaluate_objective_components(best_par_t, scenarios = scenarios, cfg = stage2_cfg)

    optim_res <- list(
      mode = "two_stage",
      stage1 = stage1_fit$optim_res,
      stage2 = stage2_fit$optim_res,
      stage1_param_names = stage1_vary,
      stage2_param_names = stage2_vary,
      optim = list(bestmem = best_par_t, bestval = stage2_comp$L)
    )
  } else {
    n_pass <- cfg$n_weight_passes
    pass_best <- warm_start_t
    pass_logs <- vector("list", n_pass)
    for (pass_i in seq_len(n_pass)) {
      pass_cfg <- cfg
      pass_cfg$w_burden <- cfg$w_burden_schedule[[pass_i]]
      pass_cfg$w_ploidy <- cfg$w_ploidy_schedule[[pass_i]]
      pass_label <- if (n_pass == 1L) "single_stage" else paste0("single_stage_pass", pass_i)
      message(
        "[", pass_label, "] weights: w_burden=", pass_cfg$w_burden,
        ", w_ploidy=", pass_cfg$w_ploidy
      )
      objective_fn <- function(par) evaluate_objective(par, scenarios = scenarios, cfg = pass_cfg)
      single_fit <- run_optimizer(
        objective_fn = objective_fn,
        lower = bounds$lower,
        upper = bounds$upper,
        cfg = pass_cfg,
        argv = argv,
        stage_label = pass_label,
        init_par = pass_best
      )
      pass_best <- single_fit$best_par
      pass_comp <- evaluate_objective_components(pass_best, scenarios = scenarios, cfg = pass_cfg)
      pass_logs[[pass_i]] <- list(
        pass = pass_i,
        w_burden = pass_cfg$w_burden,
        w_ploidy = pass_cfg$w_ploidy,
        objective = pass_comp$L,
        objective_data = pass_comp$L_data,
        objective_prior_raw = pass_comp$L_prior_raw,
        objective_prior = pass_comp$L_prior,
        objective_burden = pass_comp$L_b,
        objective_ploidy = pass_comp$L_p,
        objective_burden_scaled = pass_comp$L_b_scaled,
        objective_ploidy_scaled = pass_comp$L_p_scaled,
        optim = single_fit$optim_res
      )
    }
    best_par_t <- pass_best
    single_pass_log <- pass_logs
    optim_res <- list(
      mode = if (cfg$n_weight_passes > 1L) "single_stage_weight_schedule" else "single_stage",
      passes = pass_logs,
      optim = pass_logs[[length(pass_logs)]]$optim$optim
    )
  }

  best_par <- decode_params(
    best_par_t,
    fit_treatment = cfg$fit_treatment
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
    file.path(script_dir, "..", "..", "results", paste0("fit_invivo_model_O2_invivo_", run_stamp))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  params_df <- data.frame(
    parameter = names(best_par),
    value = as.numeric(best_par),
    row.names = NULL
  )
  write.table(params_df, file = file.path(out_dir, "best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (isTRUE(cfg$two_stage)) {
    stage1_params <- decode_params(
      stage1_best_par_t,
      fit_treatment = cfg$fit_treatment
    )
    stage1_df <- data.frame(
      parameter = names(stage1_params),
      value = as.numeric(stage1_params),
      row.names = NULL
    )
    write.table(stage1_df, file = file.path(out_dir, "stage1_best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  stage_map <- data.frame(
    transformed_parameter = full_names,
    optimized_in = if (isTRUE(cfg$two_stage)) {
      ifelse(full_names %in% stage1_vary, "stage1_growth", "stage2_ploidy")
    } else {
      rep(if (cfg$n_weight_passes > 1L) "single_stage_weight_schedule" else "single_stage", length(full_names))
    },
    transformed_value = as.numeric(best_par_t[full_names]),
    row.names = NULL
  )
  write.table(stage_map, file = file.path(out_dir, "fit_parameter_stages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (!is.null(single_pass_log)) {
    pass_df <- bind_rows(lapply(single_pass_log, function(x) {
      data.frame(
        pass = as.integer(x$pass),
        w_burden = as.numeric(x$w_burden),
        w_ploidy = as.numeric(x$w_ploidy),
        objective = as.numeric(x$objective),
        objective_data = as.numeric(x$objective_data),
        objective_prior_raw = as.numeric(x$objective_prior_raw),
        objective_prior = as.numeric(x$objective_prior),
        objective_burden = as.numeric(x$objective_burden),
        objective_ploidy = as.numeric(x$objective_ploidy),
        objective_burden_scaled = as.numeric(x$objective_burden_scaled),
        objective_ploidy_scaled = as.numeric(x$objective_ploidy_scaled),
        row.names = NULL
      )
    }))
    write.table(pass_df, file = file.path(out_dir, "single_stage_pass_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  n_ploidy_scenarios <- sum(vapply(scenarios, function(s) length(s$ploidy_obs_N) > 0, logical(1)))
  fit_mode <- if (!is.null(optim_res$mode)) as.character(optim_res$mode) else if (isTRUE(cfg$two_stage)) "two_stage" else "single_stage"
  summary_df <- data.frame(
    metric = c(
      "fit_mode",
      "weight_passes",
      "w_burden_schedule",
      "w_ploidy_schedule",
      "objective",
      "objective_data",
      "objective_prior_raw",
      "objective_prior",
      "objective_burden",
      "objective_ploidy",
      "objective_burden_scaled",
      "objective_ploidy_scaled",
      "stage1_objective",
      "stage1_objective_data",
      "stage1_objective_prior_raw",
      "stage1_objective_prior",
      "stage1_objective_burden",
      "stage1_objective_ploidy",
      "stage1_objective_burden_scaled",
      "stage1_objective_ploidy_scaled",
      "stage2_objective",
      "stage2_objective_data",
      "stage2_objective_prior_raw",
      "stage2_objective_prior",
      "stage2_objective_burden",
      "stage2_objective_ploidy",
      "stage2_objective_burden_scaled",
      "stage2_objective_ploidy_scaled",
      "stage1_w_burden",
      "stage1_w_ploidy",
      "stage2_w_burden",
      "stage2_w_ploidy",
      "loss_rescale",
      "loss_scale_burden",
      "loss_scale_ploidy",
      "loss_scale_source",
      "loss_scale_eps",
      "burden_exclude_day0",
      "scenario_agg_burden",
      "scenario_agg_ploidy",
      "scenario_agg_huber_k",
      "scenario_weight_burden",
      "scenario_weight_ploidy",
      "use_soft_prior",
      "lambda_prior",
      "prior_center_log10_k_o",
      "prior_sd_log10_k_o",
      "prior_center_log10_K_down",
      "prior_sd_log10_K_down",
      "prior_center_log10_h_down",
      "prior_sd_log10_h_down",
      "prior_center_beta_size",
      "prior_sd_beta_size",
      "prior_center_log10_n_exp",
      "prior_sd_log10_n_exp",
      "prior_center_log10_rho_2N",
      "prior_sd_log10_rho_2N",
      "n_scenarios",
      "n_ploidy_scenarios",
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
      "o2_min",
      "h_down_init",
      "h_down_min",
      "h_down_max",
      "O2_fixed",
      "final_cache_g_build",
      "final_cache_g_hit",
      "final_cache_g_hysteresis",
      "paired_only",
      "init_params_tsv",
      "append_timestamp_out_dir",
      "timestamp_format",
      "dose_zero_only",
      "truncate_at_treatment",
      "ploidy_at_harvest",
      "two_stage"
    ),
    value = c(
      fit_mode,
      as.character(cfg$n_weight_passes),
      paste(cfg$w_burden_schedule, collapse = ","),
      paste(cfg$w_ploidy_schedule, collapse = ","),
      as.character(final_obj),
      as.character(final_comp$L_data),
      as.character(final_comp$L_prior_raw),
      as.character(final_comp$L_prior),
      as.character(final_comp$L_b),
      as.character(final_comp$L_p),
      as.character(final_comp$L_b_scaled),
      as.character(final_comp$L_p_scaled),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_data else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_prior_raw else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_prior else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_b else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_p else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_b_scaled else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_p_scaled else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L_data else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L_prior_raw else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L_prior else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L_b else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L_p else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L_b_scaled else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L_p_scaled else NA_real_),
      as.character(cfg$stage1_w_burden),
      as.character(cfg$stage1_w_ploidy),
      as.character(cfg$stage2_w_burden),
      as.character(cfg$stage2_w_ploidy),
      as.character(cfg$loss_rescale),
      as.character(cfg$loss_scale_burden),
      as.character(cfg$loss_scale_ploidy),
      as.character(cfg$loss_scale_source),
      as.character(cfg$loss_scale_eps),
      as.character(cfg$burden_exclude_day0),
      as.character(cfg$scenario_agg_burden),
      as.character(cfg$scenario_agg_ploidy),
      as.character(cfg$scenario_agg_huber_k),
      as.character(cfg$scenario_weight_burden),
      as.character(cfg$scenario_weight_ploidy),
      as.character(cfg$use_soft_prior),
      as.character(cfg$lambda_prior),
      as.character(cfg$prior_center_log10_k_o),
      as.character(cfg$prior_sd_log10_k_o),
      as.character(cfg$prior_center_log10_K_down),
      as.character(cfg$prior_sd_log10_K_down),
      as.character(cfg$prior_center_log10_h_down),
      as.character(cfg$prior_sd_log10_h_down),
      as.character(cfg$prior_center_beta_size),
      as.character(cfg$prior_sd_beta_size),
      as.character(cfg$prior_center_log10_n_exp),
      as.character(cfg$prior_sd_log10_n_exp),
      as.character(cfg$prior_center_log10_rho_2N),
      as.character(cfg$prior_sd_log10_rho_2N),
      as.character(length(scenarios)),
      as.character(n_ploidy_scenarios),
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
      as.character(cfg$o2_min),
      as.character(cfg$h_down_init),
      as.character(cfg$h_down_min),
      as.character(cfg$h_down_max),
      as.character(cfg$O2_fixed),
      as.character(.first_non_null_local(final_comp$cache_g_build, NA_integer_)),
      as.character(.first_non_null_local(final_comp$cache_g_hit, NA_integer_)),
      as.character(.first_non_null_local(final_comp$cache_g_hysteresis, NA_integer_)),
      as.character(cfg$paired_only),
      as.character(if (is.null(init_params_tsv)) NA_character_ else normalizePath(init_params_tsv, mustWork = FALSE)),
      as.character(append_ts_out_dir),
      as.character(ts_format),
      as.character(cfg$dose_zero_only),
      as.character(cfg$truncate_at_treatment),
      as.character(cfg$ploidy_at_harvest),
      as.character(cfg$two_stage)
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
