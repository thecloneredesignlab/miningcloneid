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

decode_params <- function(par_transformed, fit_full_pmis = FALSE, fit_treatment = TRUE) {
  if (isTRUE(fit_treatment)) {
    names(par_transformed) <- c(
      "log10_lam_min",
      "log10_lam_max",
      "log10_k_o",
      "log10_p_misseg",
      "log10_k_o_mis",
      "beta_buffer",
      "log10_n_exp",
      "log10_smax",
      "log10_p_wgd",
      "log10_alpha",
      "gamma"
    )
    return(list(
      lam_min = 10^par_transformed["log10_lam_min"],
      lam_max = 10^par_transformed["log10_lam_max"],
      k_o = 10^par_transformed["log10_k_o"],
      p_misseg = 10^par_transformed["log10_p_misseg"],
      k_o_mis = 10^par_transformed["log10_k_o_mis"],
      beta_buffer = par_transformed["beta_buffer"],
      n_exp = 10^par_transformed["log10_n_exp"],
      smax = 10^par_transformed["log10_smax"],
      p_wgd = 10^par_transformed["log10_p_wgd"],
      alpha = 10^par_transformed["log10_alpha"],
      gamma = par_transformed["gamma"]
    ))
  }

  names(par_transformed) <- c(
    "log10_lam_min",
    "log10_lam_max",
    "log10_k_o",
    "log10_p_misseg",
    "log10_k_o_mis",
    "beta_buffer",
    "log10_n_exp",
    "log10_smax",
    "log10_p_wgd"
  )
  list(
    lam_min = 10^par_transformed["log10_lam_min"],
    lam_max = 10^par_transformed["log10_lam_max"],
    k_o = 10^par_transformed["log10_k_o"],
    p_misseg = 10^par_transformed["log10_p_misseg"],
    k_o_mis = 10^par_transformed["log10_k_o_mis"],
    beta_buffer = par_transformed["beta_buffer"],
    n_exp = 10^par_transformed["log10_n_exp"],
    smax = 10^par_transformed["log10_smax"],
    p_wgd = 10^par_transformed["log10_p_wgd"],
    alpha = 0,
    gamma = 1
  )
}

encode_params <- function(run_params, fit_full_pmis = FALSE, fit_treatment = TRUE) {
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

  lam_min_v <- need_pos(getv(c("lam_min", "R")), "lam_min")
  lam_max_v <- need_pos(getv(c("lam_max", "R"), default = lam_min_v), "lam_max")
  k_o_v <- need_pos(getv(c("k_o"), default = 0.5), "k_o")
  p_misseg_v <- need_pos(getv(c("p_misseg", "pmis_O2_0", "pmis_O2_1", "pmis", "p_mis"), default = 1e-4), "p_misseg")
  k_o_mis_v <- need_pos(getv(c("k_o_mis"), default = 0.5), "k_o_mis")
  beta_buffer_v <- getv(c("beta_buffer", "beta"), default = 0.0)
  if (!is.finite(beta_buffer_v)) stop("Warm-start parameter must be finite: beta_buffer")
  n_exp_v <- need_pos(getv(c("n_exp"), default = 1.0), "n_exp")
  smax_v <- need_pos(getv(c("smax"), default = 1.0), "smax")
  p_wgd_v <- need_pos(getv(c("p_wgd", "pwgd"), default = 1e-6), "p_wgd")

  if (isTRUE(fit_treatment)) {
    alphav <- need_pos(getv(c("alpha")), "alpha")
    gammav <- getv(c("gamma"))
    return(c(
      log10_lam_min = log10(lam_min_v),
      log10_lam_max = log10(lam_max_v),
      log10_k_o = log10(k_o_v),
      log10_p_misseg = log10(p_misseg_v),
      log10_k_o_mis = log10(k_o_mis_v),
      beta_buffer = beta_buffer_v,
      log10_n_exp = log10(n_exp_v),
      log10_smax = log10(smax_v),
      log10_p_wgd = log10(p_wgd_v),
      log10_alpha = log10(alphav),
      gamma = gammav
    ))
  }

  c(
    log10_lam_min = log10(lam_min_v),
    log10_lam_max = log10(lam_max_v),
    log10_k_o = log10(k_o_v),
    log10_p_misseg = log10(p_misseg_v),
    log10_k_o_mis = log10(k_o_mis_v),
    beta_buffer = beta_buffer_v,
    log10_n_exp = log10(n_exp_v),
    log10_smax = log10(smax_v),
    log10_p_wgd = log10(p_wgd_v)
  )
}

read_init_params_t <- function(init_path, bounds, cfg) {
  if (!file.exists(init_path)) stop("init_params_tsv not found: ", init_path)
  tab <- read.delim(init_path, check.names = FALSE, stringsAsFactors = FALSE)
  full_names <- names(bounds$lower)

  out <- NULL
  if (all(c("transformed_parameter", "transformed_value") %in% names(tab))) {
    vals <- setNames(as.numeric(tab$transformed_value), as.character(tab$transformed_parameter))
    if (all(full_names %in% names(vals))) out <- vals[full_names]
  }
  if (is.null(out) && all(c("parameter", "value") %in% names(tab))) {
    vals <- setNames(as.numeric(tab$value), as.character(tab$parameter))
    if (all(full_names %in% names(vals))) {
      out <- vals[full_names]
    } else {
      out <- encode_params(
        vals,
        fit_full_pmis = cfg$fit_full_pmis,
        fit_treatment = cfg$fit_treatment
      )
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

make_bounds <- function(fit_full_pmis = FALSE, fit_treatment = TRUE) {
  # fit_full_pmis is retained only for CLI/backward compatibility in this script.
  # Richard-aligned parameterization always fits p_misseg and k_o_mis directly.
  lower <- c(
    log10_lam_min = log10(1e-3),
    log10_lam_max = log10(1e-3),
    log10_k_o = log10(1e-3),
    log10_p_misseg = log10(1e-8),
    log10_k_o_mis = log10(1e-3),
    beta_buffer = 0.0,
    log10_n_exp = log10(1e-1),
    log10_smax = log10(1e-4),
    log10_p_wgd = log10(1e-8)
  )
  upper <- c(
    log10_lam_min = log10(5),
    log10_lam_max = log10(5),
    log10_k_o = log10(1e2),
    log10_p_misseg = log10(5e-1),
    log10_k_o_mis = log10(1e2),
    beta_buffer = 10.0,
    log10_n_exp = log10(5),
    log10_smax = log10(1),
    log10_p_wgd = log10(1e-1)
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

build_model_core <- function(run_params, cfg) {
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

  lam_min <- as.numeric(run_params$lam_min)
  lam_max <- as.numeric(run_params$lam_max)
  k_o <- as.numeric(run_params$k_o)
  p_misseg <- as.numeric(run_params$p_misseg)
  k_o_mis <- as.numeric(run_params$k_o_mis)
  beta_buffer <- as.numeric(.first_non_null_local(run_params$beta_buffer, run_params$beta, 0))
  n_exp <- as.numeric(.first_non_null_local(run_params$n_exp, 1.0))
  smax <- as.numeric(.first_non_null_local(run_params$smax, 1.0))
  p_wgd <- as.numeric(.first_non_null_local(run_params$p_wgd, run_params$pwgd, 0))

  lambda0 <- growth_lambda(
    cfg$O2_fixed, grid_pre,
    lam_min = lam_min, lam_max = lam_max, k_o = k_o,
    R = lam_min, beta = beta_buffer, eta = 0, N_unit = cfg$N_UNIT
  )
  lambda1 <- growth_lambda(
    cfg$O2_fixed, grid_post,
    lam_min = lam_min, lam_max = lam_max, k_o = k_o,
    R = lam_min, beta = beta_buffer, eta = 0, N_unit = cfg$N_UNIT
  )

  if (exists(".pmisseg_of_O2", mode = "function")) {
    p_mis <- as.numeric(.pmisseg_of_O2(cfg$O2_fixed, run_params))
  } else {
    k_o_mis <- max(k_o_mis, 1e-12)
    p_mis <- p_misseg * (1 - (cfg$O2_fixed / (cfg$O2_fixed + k_o_mis)))
  }
  p_mis <- clip(p_mis, 0, 1)

  G <- .build_G_with_WGD(
    N0min = cfg$N_MIN, N0max = cfg$N_MAX, lambda0_vec = lambda0,
    p0_vec = p_mis, wgd_prob_vec = p_wgd,
    N1min = cfg$N_MIN, N1max = cfg$N_MAX, lambda1_vec = lambda1,
    p1_vec = p_mis,
    boundary = "drop",
    N_unit = cfg$N_UNIT,
    beta_buffer = beta_buffer,
    n_exp = n_exp,
    smax = smax
  )

  list(
    grid_pre = grid_pre,
    grid_post = grid_post,
    R0 = R0,
    R1 = R1,
    init_state_2N = init_state_2N,
    init_state_4N = init_state_4N,
    G = G,
    I = Diagonal(n = length(init_state_2N))
  )
}

simulate_one <- function(run_params, scenario, cfg, model_core = NULL) {
  if (is.null(model_core)) {
    model_core <- build_model_core(run_params, cfg)
  }

  grid_pre <- model_core$grid_pre
  R0 <- model_core$R0
  R1 <- model_core$R1
  init_state <- if (scenario$cohort == "2N") model_core$init_state_2N else model_core$init_state_4N
  G <- model_core$G
  I <- model_core$I

  obs_steps <- as.integer(round(scenario$obs_days / cfg$DT))
  sim_end_step <- as.integer(round(scenario$sim_end_day / cfg$DT))
  step_unique <- sort(unique(obs_steps))
  final_step <- max(sim_end_step, max(step_unique))
  step_to_idx <- setNames(seq_along(step_unique), as.character(step_unique))
  Ntot_at_step <- rep(NA_real_, length(step_unique))

  v <- as.numeric(init_state)
  dose_scaled <- scenario$dose / cfg$dose_ref
  if (!is.finite(dose_scaled) || dose_scaled < 0) dose_scaled <- 0

  for (step in 0:final_step) {
    key <- as.character(step)
    i_obs <- step_to_idx[key]
    if (!is.na(i_obs)) Ntot_at_step[[i_obs]] <- sum(v)
    if (step >= final_step) break

    t <- step * cfg$DT
    tx_mult <- 1.0
    if (isTRUE(cfg$fit_treatment)) {
      tx_mult <- if (t < scenario$treat_day || scenario$dose <= 0) {
        1.0
      } else {
        exp(-run_params$alpha * (dose_scaled^run_params$gamma))
      }
      tx_mult <- clip(tx_mult, cfg$tx_mult_min, 1.0)
    }

    Ntot <- sum(v)
    crowd <- if (cfg$crowding == "logistic") max(0, 1 - Ntot / cfg$K) else exp(-Ntot / cfg$K)
    A <- I + cfg$DT * (crowd * tx_mult * G)
    v <- as.numeric(A %*% v)
    v[!is.finite(v)] <- 0
    v[v < 0] <- 0
    if (sum(v) <= cfg$min_pop) break
  }

  Ntot_obs <- Ntot_at_step[match(obs_steps, step_unique)]
  if (any(!is.finite(Ntot_obs))) {
    # Fallback in case integration stopped early
    Ntot_obs[!is.finite(Ntot_obs)] <- cfg$min_pop
  }

  frac_pre <- v[seq_len(R0)]
  frac_post <- v[R0 + seq_len(R1)]
  frac_N <- frac_pre + frac_post
  if (sum(frac_N) > 0) {
    frac_N <- frac_N / sum(frac_N)
  } else {
    frac_N <- rep(1 / length(frac_N), length(frac_N))
  }
  names(frac_N) <- as.character(grid_pre)

  list(
    Ntot_obs = as.numeric(Ntot_obs),
    frac_N = frac_N
  )
}

evaluate_objective_components_raw <- function(par_transformed, scenarios, cfg) {
  rp <- decode_params(
    par_transformed,
    fit_full_pmis = cfg$fit_full_pmis,
    fit_treatment = cfg$fit_treatment
  )
  model_core <- build_model_core(rp, cfg)
  burden_losses <- numeric(0)
  ploidy_losses <- numeric(0)

  for (sc in scenarios) {
    sim <- simulate_one(rp, sc, cfg, model_core = model_core)

    obs <- sc$obs_burden
    pred <- sim$Ntot_obs
    if (length(obs) >= 2 && length(pred) == length(obs) && all(is.finite(pred))) {
      obs_delta <- obs - obs[1]
      pred_delta <- pred - pred[1]
      s_obs <- max(abs(obs_delta), na.rm = TRUE)
      s_pred <- max(abs(pred_delta), na.rm = TRUE)
      if (s_obs > 0 && s_pred > 0) {
        r <- (obs_delta / s_obs) - (pred_delta / s_pred)
        burden_losses <- c(burden_losses, huber_mean(r, k = cfg$huber_k))
      }
    }

    if (length(sc$ploidy_obs_N) > 0) {
      obs_tab <- table(as.character(sc$ploidy_obs_N))
      p <- sim$frac_N[names(obs_tab)]
      p[is.na(p)] <- cfg$eps_prob
      nll <- -sum(as.numeric(obs_tab) * log(p + cfg$eps_prob)) / sum(obs_tab)
      ploidy_losses <- c(ploidy_losses, nll)
    }
  }

  L_b <- if (length(burden_losses) > 0) mean(burden_losses) else 0
  L_p <- if (length(ploidy_losses) > 0) mean(ploidy_losses) else 0
  list(L_b = L_b, L_p = L_p)
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
  L <- cfg$w_burden * L_b_scaled + cfg$w_ploidy * L_p_scaled
  if (!is.finite(L)) L <- 1e9

  if (cfg$trace_obj) {
    cat(sprintf(
      "L=%.6f (burden=%.6f, ploidy=%.6f; scaled burden=%.6f, scaled ploidy=%.6f)\n",
      L, L_b, L_p, L_b_scaled, L_p_scaled
    ))
  }
  list(
    L = L,
    L_b = L_b,
    L_p = L_p,
    L_b_scaled = L_b_scaled,
    L_p_scaled = L_p_scaled,
    scale_b = scale_b,
    scale_p = scale_p
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
    if (!is.null(cfg$model_path) && nzchar(cfg$model_path) && file.exists(cfg$model_path)) {
      parallel::clusterCall(cl, function(path) {
        Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(path))
        source(path)
        NULL
      }, cfg$model_path)
    }
    export_global <- c(
      "evaluate_objective",
      "evaluate_objective_components",
      "evaluate_objective_components_raw",
      "build_model_core",
      "simulate_one",
      "decode_params",
      "huber_mean",
      "clip",
      ".first_non_null_local"
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
    de_extra <- paste0("NP=", cfg$NP, ", itermax=", cfg$itermax)
    message(
      "[", stage_label, "] Starting DEoptim with itermax=", cfg$itermax,
      ", NP=", cfg$NP,
      ", n_cores=", n_cores
    )
    de_ctrl <- list(
      trace = TRUE,
      itermax = cfg$itermax,
      NP = cfg$NP,
      strategy = 2
    )
    if (!is.null(init_use)) {
      initpop <- matrix(
        stats::runif(cfg$NP * length(lower), min = lower, max = upper),
        nrow = cfg$NP, byrow = TRUE
      )
      initpop[1, ] <- init_use
      colnames(initpop) <- names(lower)
      de_ctrl$initialpop <- initpop
    }
    if (n_cores > 1L) {
      message("[", stage_label, "] DEoptim parallel requested with n_cores=", n_cores, ".")
      cl <- tryCatch(
        parallel::makePSOCKcluster(n_cores),
        error = function(e) {
          message("[", stage_label, "] Could not start parallel workers for DEoptim: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(cl)) {
        de_started <- cluster_workers(cl)
        init_ok <- tryCatch(
          init_cluster_workers(cl, objective_fn = objective_fn, cfg = cfg, stage_label = stage_label),
          error = function(e) {
            message("[", stage_label, "] Failed to initialize DEoptim workers: ", conditionMessage(e))
            FALSE
          }
        )
        if (isTRUE(init_ok)) {
          on.exit(parallel::stopCluster(cl), add = TRUE)
          # IMPORTANT: when passing an explicit cluster, do NOT set
          # parallelType='parallel', otherwise DEoptim may try stopCluster(cl)
          # on an internal symbol 'cl' that does not exist in this code path.
          de_ctrl$cluster <- cl
          de_active <- ifelse(is.na(de_started) || de_started < 1L, 1L, de_started)
        } else {
          try(parallel::stopCluster(cl), silent = TRUE)
          de_fallback <- TRUE
          de_fallback_reason <- "worker_init_failed"
          if (is.na(de_started)) de_started <- 0L
          de_active <- 1L
        }
      } else {
        de_started <- 0L
        de_fallback <- TRUE
        de_fallback_reason <- "cluster_start_failed"
        de_active <- 1L
      }
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
          message("[", stage_label, "] Could not start parallel workers for optim fallback: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(cl)) {
        opt_started <- cluster_workers(cl)
        init_ok <- tryCatch(
          init_cluster_workers(cl, objective_fn = objective_fn, cfg = cfg, stage_label = stage_label),
          error = function(e) {
            message("[", stage_label, "] Failed to initialize optim workers: ", conditionMessage(e))
            FALSE
          }
        )
        if (isTRUE(init_ok)) {
          on.exit(parallel::stopCluster(cl), add = TRUE)
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
          try(parallel::stopCluster(cl), silent = TRUE)
          opt_fallback <- TRUE
          opt_fallback_reason <- "worker_init_failed"
          if (is.na(opt_started)) opt_started <- 0L
          opt_active <- 1L
          serial_out <- run_starts_serial(starts)
          run_log <- serial_out$run_log
          best_val <- serial_out$best_val
          best_par <- serial_out$best_par
        }
      } else {
        opt_started <- 0L
        opt_fallback <- TRUE
        opt_fallback_reason <- "cluster_start_failed"
        opt_active <- 1L
        serial_out <- run_starts_serial(starts)
        run_log <- serial_out$run_log
        best_val <- serial_out$best_val
        best_par <- serial_out$best_par
      }
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
  model_core <- build_model_core(run_params, cfg)
  burden_rows <- list()
  ploidy_rows <- list()

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    sim <- simulate_one(run_params, sc, cfg, model_core = model_core)

    obs <- sc$obs_burden
    pred <- sim$Ntot_obs
    obs_delta <- obs - obs[1]
    pred_delta <- pred - pred[1]
    s_obs <- max(abs(obs_delta), na.rm = TRUE)
    s_pred <- max(abs(pred_delta), na.rm = TRUE)
    obs_norm <- if (s_obs > 0) obs_delta / s_obs else obs_delta
    pred_norm <- if (s_pred > 0) pred_delta / s_pred else pred_delta

    burden_rows[[length(burden_rows) + 1]] <- data.frame(
      harvest = sc$harvest,
      cohort = sc$cohort,
      dose = sc$dose,
      treat_day = sc$treat_day,
      day = sc$obs_days,
      obs_burden = obs,
      pred_pop = pred,
      obs_norm = obs_norm,
      pred_norm = pred_norm
    )

    if (length(sc$ploidy_obs_N) > 0) {
      obs_tab <- table(sc$ploidy_obs_N)
      dfp <- data.frame(
        harvest = sc$harvest,
        cohort = sc$cohort,
        dose = sc$dose,
        N = as.integer(names(sim$frac_N)),
        pred_fraction = as.numeric(sim$frac_N)
      )
      dfp$obs_count <- as.numeric(obs_tab[as.character(dfp$N)])
      dfp$obs_count[is.na(dfp$obs_count)] <- 0
      ploidy_rows[[length(ploidy_rows) + 1]] <- dfp
    }
  }

  list(
    burden = bind_rows(burden_rows),
    ploidy = bind_rows(ploidy_rows)
  )
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  script_dir <- get_script_dir()

  model_path <- file.path(script_dir, "model_buffering_align_with_Richard.R")
  if (!file.exists(model_path)) stop("Cannot find model_buffering_align_with_Richard.R at ", model_path)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = script_dir)
  source(model_path)

  default_data_dir <- normalizePath(file.path(script_dir, "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
  data_dir <- if (!is.null(argv$data_dir)) argv$data_dir else default_data_dir
  truncate_at_treatment <- if (!is.null(argv$truncate_at_treatment)) {
    as_bool(argv$truncate_at_treatment, FALSE)
  } else {
    # backward compatibility: old flag name
    as_bool(argv$pretreat_only, FALSE)
  }
  w_burden_vec <- as_num_vec(argv$w_burden, 1.0)
  w_ploidy_vec <- as_num_vec(argv$w_ploidy, 1.0)
  weight_schedule <- make_weight_schedule(w_burden_vec, w_ploidy_vec)
  n_cores_arg <- as_int(argv$n_cores, NA_integer_)
  n_cores_use <- if (is.finite(n_cores_arg)) n_cores_arg else default_n_cores()

  cfg <- list(
    # model constants
    model_path = model_path,
    N_UNIT = as_int(argv$N_UNIT, 22L),
    N_MIN = as_int(argv$N_MIN, 22L),
    N_MAX = as_int(argv$N_MAX, 154L),
    DT = as_num(argv$dt, 0.5),
    O2_fixed = as_num(argv$O2, 1.0),
    K = as_num(argv$K, 1e12),
    crowding = if (!is.null(argv$crowding)) argv$crowding else "logistic",
    init_total_size = as_num(argv$init_total_size, 1e6),
    dose_ref = as_num(argv$dose_ref, 30),
    tx_mult_min = as_num(argv$tx_mult_min, 0.05),
    min_pop = as_num(argv$min_pop, 1e-12),
    # objective settings
    huber_k = as_num(argv$huber_k, 0.1),
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
    fit_full_pmis = as_bool(argv$fit_full_pmis, FALSE),
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
    itermax = as_int(argv$itermax, 40L),
    NP = as_int(argv$NP, 80L),
    n_cores = n_cores_use,
    seed = as_int(argv$seed, 1L),
    max_scenarios = as_num(argv$max_scenarios, Inf)
  )

  if (!cfg$crowding %in% c("logistic", "gompertz")) stop("crowding must be logistic or gompertz")
  if (cfg$DT <= 0) stop("dt must be > 0")
  if (cfg$N_MAX < cfg$N_MIN) stop("N_MAX must be >= N_MIN")
  if (cfg$itermax < 1) stop("itermax must be >= 1")
  if (cfg$n_cores < 1) stop("n_cores must be >= 1")
  if (cfg$optim_trace_every < 1) stop("optim_trace_every must be >= 1")
  if (!is.finite(cfg$loss_scale_eps) || cfg$loss_scale_eps <= 0) stop("loss_scale_eps must be > 0")
  if (!is.na(cfg$loss_scale_burden) && (!is.finite(cfg$loss_scale_burden) || cfg$loss_scale_burden <= 0)) stop("loss_scale_burden must be > 0")
  if (!is.na(cfg$loss_scale_ploidy) && (!is.finite(cfg$loss_scale_ploidy) || cfg$loss_scale_ploidy <= 0)) stop("loss_scale_ploidy must be > 0")
  if (!cfg$use_deoptim && cfg$deoptim_parallel) stop("deoptim_parallel=TRUE requires use_deoptim=TRUE")
  if (isTRUE(cfg$two_stage) && cfg$n_weight_passes > 1L) {
    stop("Weight arrays for w_burden/w_ploidy are supported only when two_stage=FALSE.")
  }

  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- file.path(data_dir, "all_ploidy.tsv")
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)

  bounds <- make_bounds(
    fit_full_pmis = cfg$fit_full_pmis,
    fit_treatment = cfg$fit_treatment
  )
  full_names <- names(bounds$lower)
  default_par_t <- (bounds$lower + bounds$upper) / 2
  names(default_par_t) <- full_names
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
  set.seed(cfg$seed)
  stage1_comp <- NULL
  stage2_comp <- NULL
  stage1_best_par_t <- NULL
  stage1_vary <- character(0)
  stage2_vary <- character(0)
  single_pass_log <- NULL

  if (isTRUE(cfg$fit_full_pmis)) {
    message("fit_full_pmis is ignored in this script; Richard parameterization directly fits p_misseg and k_o_mis.")
  }

  if (isTRUE(cfg$two_stage)) {
    stage1_vary <- intersect(full_names, c("log10_lam_min", "log10_lam_max", "log10_k_o"))
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
    fit_full_pmis = cfg$fit_full_pmis,
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
    file.path(script_dir, "..", "results", paste0("fit_invivo_model_buffering_align_with_Richard_", run_stamp))
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
      fit_full_pmis = cfg$fit_full_pmis,
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
      "objective_burden",
      "objective_ploidy",
      "objective_burden_scaled",
      "objective_ploidy_scaled",
      "stage1_objective",
      "stage1_objective_burden",
      "stage1_objective_ploidy",
      "stage1_objective_burden_scaled",
      "stage1_objective_ploidy_scaled",
      "stage2_objective",
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
      "n_scenarios",
      "n_ploidy_scenarios",
      "itermax",
      "NP",
      "n_cores",
      "optim_trace",
      "optim_trace_every",
      "use_deoptim",
      "deoptim_parallel",
      "fit_full_pmis",
      "fit_treatment",
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
      as.character(final_comp$L_b),
      as.character(final_comp$L_p),
      as.character(final_comp$L_b_scaled),
      as.character(final_comp$L_p_scaled),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_b else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_p else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_b_scaled else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage1_comp$L_p_scaled else NA_real_),
      as.character(if (isTRUE(cfg$two_stage)) stage2_comp$L else NA_real_),
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
      as.character(length(scenarios)),
      as.character(n_ploidy_scenarios),
      as.character(cfg$itermax),
      as.character(cfg$NP),
      as.character(cfg$n_cores),
      as.character(cfg$optim_trace),
      as.character(cfg$optim_trace_every),
      as.character(cfg$use_deoptim),
      as.character(cfg$deoptim_parallel),
      as.character(cfg$fit_full_pmis),
      as.character(cfg$fit_treatment),
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
