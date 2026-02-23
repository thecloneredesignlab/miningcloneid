suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

# ----------------------------------------------------------------------------
# Align miningcloneid oxygen model to Richard's buffering.R logic.
# - Keep outer API from model_functions_ploidy_buffer.R
# - Override core dynamics (lambda, misseg delta, B/G construction)
# - Enforce defaults requested:
#   * WGD branch offspring weight = +1 (not +2)
#   * boundary default = "drop"
# ----------------------------------------------------------------------------

.resolve_align_script_dir <- function() {
  env_dir <- Sys.getenv("MININGCLONEID_OXYGEN_CODE_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = FALSE))
  }

  frs <- sys.frames()
  for (i in rev(seq_along(frs))) {
    ofile <- frs[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
  }

  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grepl("^--file=", args)]
  if (length(farg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", farg[[1]]), mustWork = FALSE)))
  }

  getwd()
}

.ALIGN_MODEL_DIR <- .resolve_align_script_dir()
.BASE_MODEL_PATH <- file.path(.ALIGN_MODEL_DIR, "model_functions_ploidy_buffer.R")
if (!file.exists(.BASE_MODEL_PATH)) {
  stop("Cannot find base model file: ", .BASE_MODEL_PATH)
}

# Load full base API first, then override core dynamics below.
source(.BASE_MODEL_PATH)

.init_cpp_richard_backend <- local({
  initialized <- FALSE
  available <- FALSE
  function() {
    if (initialized) return(available)
    initialized <<- TRUE

    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      available <<- FALSE
      return(available)
    }

    # Dedicated Richard backend (do not use model_ploidy_buffer_cpp.cpp here).
    cpp_path <- file.path(.ALIGN_MODEL_DIR, "model_buffering_align_with_Richard.cpp")
    if (!file.exists(cpp_path)) {
      available <<- FALSE
      return(available)
    }

    cache_dir <- file.path(.ALIGN_MODEL_DIR, ".rcpp_cache_richard")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    available <<- isTRUE(tryCatch({
      Rcpp::sourceCpp(
        file = cpp_path,
        rebuild = FALSE,
        showOutput = FALSE,
        verbose = FALSE,
        cacheDir = cache_dir
      )
      exists("cpp_richard_pr_delta_vec", mode = "function", inherits = TRUE) &&
        exists("cpp_richard_build_B_total_triplet", mode = "function", inherits = TRUE) &&
        exists("cpp_richard_build_B_WGD_triplet", mode = "function", inherits = TRUE)
    }, error = function(e) {
      message("Richard C++ backend unavailable; fallback to R implementation: ", conditionMessage(e))
      FALSE
    }))

    available
  }
})

.USE_CPP_RICHARD_BACKEND <- .init_cpp_richard_backend()

.first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

.clip01 <- function(x) pmin(pmax(x, 0), 1)

.get_pwgd <- function(run_params) {
  as.numeric(.first_non_null(run_params$p_wgd, run_params$pwgd, 0))
}

# Richard buffering.R-style lambda: O2-only saturating rate.
# Keep legacy args in signature for compatibility with existing callers.
growth_lambda <- function(O2, N, R = 1.0, beta = 0.35, N_unit = 22L, eta = 0.01,
                          lam_min = NULL, lam_max = NULL, k_o = NULL) {
  lam_min_use <- as.numeric(.first_non_null(lam_min, R, 1.0))
  lam_max_use <- as.numeric(.first_non_null(lam_max, R, lam_min_use))
  k_o_use <- as.numeric(.first_non_null(k_o, 0.5))
  k_o_use <- max(k_o_use, 1e-12)

  frac <- O2 / (O2 + k_o_use)
  lam <- lam_min_use + (lam_max_use - lam_min_use) * frac
  rep(pmax(lam, 0), length(N))
}

# Richard buffering.R-style O2-dependent missegregation.
# Backward-compatible fallback: old log-linear interpolation between pmis_O2_0/1.
.pmisseg_of_O2 <- function(O2, run_params) {
  if (!is.null(run_params$p_misseg)) {
    p0 <- as.numeric(run_params$p_misseg)
    k_o_mis <- as.numeric(.first_non_null(run_params$k_o_mis, 0.5))
    k_o_mis <- max(k_o_mis, 1e-12)
    p <- p0 * (1 - (O2 / (O2 + k_o_mis)))
    return(.clip01(p))
  }

  if (!is.null(run_params$pmis_O2_0) && !is.null(run_params$pmis_O2_1)) {
    p0 <- as.numeric(run_params$pmis_O2_0)
    p1 <- as.numeric(run_params$pmis_O2_1)
    if (p0 <= 0 || p1 <= 0) return(0)
    logp <- (1 - O2) * log10(p0) + O2 * log10(p1)
    return(.clip01(10^logp))
  }

  .clip01(as.numeric(.first_non_null(run_params$pmis, run_params$p_mis, 0)))
}

# Richard buffering.R delta weight formula.
.pr_delta_vec <- function(N, p, eps_tail = 1e-8, mr_lethality = 0.9,
                          beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L) {
  if (isTRUE(.USE_CPP_RICHARD_BACKEND) &&
      exists("cpp_richard_pr_delta_vec", mode = "function", inherits = TRUE)) {
    res <- cpp_richard_pr_delta_vec(
      as.integer(N),
      as.numeric(p),
      eps_tail = as.numeric(eps_tail),
      beta_buffer = as.numeric(beta_buffer),
      n_exp = as.numeric(n_exp),
      smax = as.numeric(smax),
      N_unit = as.integer(N_unit)
    )
    out <- as.numeric(res$prob)
    names(out) <- as.character(res$ts)
    attr(out, "mass_dropped") <- as.numeric(res$mass_dropped)
    return(out)
  }

  if (p <= 0 || N <= 0) {
    out <- c("0" = 1)
    attr(out, "mass_dropped") <- 0
    return(out)
  }

  sd <- sqrt(N * p)
  if (sd == 0) {
    out <- c("0" = 1)
    attr(out, "mass_dropped") <- 0
    return(out)
  }

  sN <- smax * exp(-beta_buffer * ((2 * N_unit) / N)^n_exp)
  z <- 9.0
  T <- min(N, max(0L, ceiling(z * sd)))
  ts <- (-T):T
  out <- numeric(length(ts))

  for (idx in seq_along(ts)) {
    t <- ts[idx]
    ks <- seq.int(abs(t), N, by = 2)
    if (!length(ks)) next
    out[idx] <- sum(
      stats::dbinom(ks, N, p) *
        stats::dbinom((ks + t) / 2, ks, 0.5) *
        (sN^ks)
    )
  }

  names(out) <- as.character(ts)
  attr(out, "mass_dropped") <- max(0, 1 - sum(out))
  out
}

.build_B_total <- function(Nmin, Nmax, p_vec, mr_lethality = 0.9,
                           boundary = c("drop", "absorb_minmax"),
                           eps_tail = 1e-8, return_sparse = TRUE,
                           beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L) {
  boundary <- match.arg(boundary)
  R <- Nmax - Nmin + 1L
  if (length(p_vec) == 1L) p_vec <- rep(p_vec, R)

  if (isTRUE(.USE_CPP_RICHARD_BACKEND) &&
      exists("cpp_richard_build_B_total_triplet", mode = "function", inherits = TRUE)) {
    tri <- cpp_richard_build_B_total_triplet(
      as.integer(Nmin),
      as.integer(Nmax),
      as.numeric(p_vec),
      boundary = boundary,
      eps_tail = as.numeric(eps_tail),
      beta_buffer = as.numeric(beta_buffer),
      n_exp = as.numeric(n_exp),
      smax = as.numeric(smax),
      N_unit = as.integer(N_unit)
    )
    B <- sparseMatrix(
      i = as.integer(tri$i),
      j = as.integer(tri$j),
      x = as.numeric(tri$x),
      dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
      repr = "C"
    )
    return(if (isTRUE(return_sparse)) B else as.matrix(B))
  }

  ii <- integer(0)
  jj <- integer(0)
  xx <- numeric(0)

  for (col in seq_len(R)) {
    N <- Nmin + col - 1L
    pN <- .clip01(p_vec[col])
    prDelta <- .pr_delta_vec(
      N, pN, eps_tail = eps_tail, mr_lethality = mr_lethality,
      beta_buffer = beta_buffer, n_exp = n_exp, smax = smax, N_unit = N_unit
    )
    ts <- as.integer(names(prDelta))
    pr <- as.numeric(prDelta)

    for (k in seq_along(ts)) {
      t <- ts[k]
      w <- pr[k]
      if (w == 0) next

      if (t == 0L) {
        Np <- N
        val <- 2 * w
        if (Np < Nmin || Np > Nmax) {
          if (boundary == "absorb_minmax") {
            Np2 <- max(min(Np, Nmax), Nmin)
            ii <- c(ii, Np2 - Nmin + 1L)
            jj <- c(jj, col)
            xx <- c(xx, val)
          }
        } else {
          ii <- c(ii, Np - Nmin + 1L)
          jj <- c(jj, col)
          xx <- c(xx, val)
        }
      } else {
        for (Np in c(N + t, N - t)) {
          if (Np < Nmin || Np > Nmax) {
            if (boundary == "absorb_minmax") {
              Np2 <- max(min(Np, Nmax), Nmin)
              ii <- c(ii, Np2 - Nmin + 1L)
              jj <- c(jj, col)
              xx <- c(xx, w)
            }
          } else {
            ii <- c(ii, Np - Nmin + 1L)
            jj <- c(jj, col)
            xx <- c(xx, w)
          }
        }
      }
    }
  }

  B <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(R, R), repr = "C")
  if (isTRUE(return_sparse)) B else as.matrix(B)
}

.build_B_WGD <- function(N0min, N0max, N1min, N1max,
                         boundary = c("drop", "absorb_minmax"),
                         return_sparse = TRUE) {
  boundary <- match.arg(boundary)
  R0 <- N0max - N0min + 1L
  R1 <- N1max - N1min + 1L

  if (isTRUE(.USE_CPP_RICHARD_BACKEND) &&
      exists("cpp_richard_build_B_WGD_triplet", mode = "function", inherits = TRUE)) {
    tri <- cpp_richard_build_B_WGD_triplet(
      as.integer(N0min),
      as.integer(N0max),
      as.integer(N1min),
      as.integer(N1max),
      boundary = boundary,
      wgd_value = 1.0
    )
    B <- sparseMatrix(
      i = as.integer(tri$i),
      j = as.integer(tri$j),
      x = as.numeric(tri$x),
      dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
      repr = "C"
    )
    return(if (isTRUE(return_sparse)) B else as.matrix(B))
  }

  ii <- integer(0)
  jj <- integer(0)
  xx <- numeric(0)

  for (col in seq_len(R0)) {
    N0 <- N0min + col - 1L
    Np <- 2L * N0
    val <- 1.0 # Strictly align with buffering.R requested behavior.

    if (Np < N1min || Np > N1max) {
      if (boundary == "absorb_minmax") {
        Np2 <- max(min(Np, N1max), N1min)
        ii <- c(ii, Np2 - N1min + 1L)
        jj <- c(jj, col)
        xx <- c(xx, val)
      }
    } else {
      ii <- c(ii, Np - N1min + 1L)
      jj <- c(jj, col)
      xx <- c(xx, val)
    }
  }

  B <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(R1, R0), repr = "C")
  if (isTRUE(return_sparse)) B else as.matrix(B)
}

.build_G_with_WGD <- function(
    N0min, N0max, lambda0_vec, p0_vec, wgd_prob_vec,
    N1min, N1max, lambda1_vec, p1_vec,
    mr_lethality0 = 0.9, mr_lethality1 = 0.9,
    mr_buffer_by_ploidy = TRUE, N_unit = 22L, P_low = 2.0, P_high = 4.0,
    boundary = "drop", eps_tail = 1e-8,
    beta_buffer = 0.0, n_exp = 1.0, smax = 1.0
) {
  R0 <- N0max - N0min + 1L
  R1 <- N1max - N1min + 1L
  if (length(lambda0_vec) == 1L) lambda0_vec <- rep(lambda0_vec, R0)
  if (length(p0_vec) == 1L) p0_vec <- rep(p0_vec, R0)
  if (length(wgd_prob_vec) == 1L) wgd_prob_vec <- rep(wgd_prob_vec, R0)
  if (length(lambda1_vec) == 1L) lambda1_vec <- rep(lambda1_vec, R1)
  if (length(p1_vec) == 1L) p1_vec <- rep(p1_vec, R1)
  wgd_prob_vec <- .clip01(wgd_prob_vec)

  B0 <- .build_B_total(
    N0min, N0max, p_vec = p0_vec,
    boundary = boundary, eps_tail = eps_tail,
    beta_buffer = beta_buffer, n_exp = n_exp, smax = smax, N_unit = N_unit
  )
  B1 <- .build_B_total(
    N1min, N1max, p_vec = p1_vec,
    boundary = boundary, eps_tail = eps_tail,
    beta_buffer = beta_buffer, n_exp = n_exp, smax = smax, N_unit = N_unit
  )
  BW <- .build_B_WGD(N0min, N0max, N1min, N1max, boundary = boundary)

  L0 <- Diagonal(x = lambda0_vec)
  L1 <- Diagonal(x = lambda1_vec)
  S0 <- Diagonal(x = (1 - wgd_prob_vec))
  SW <- Diagonal(x = wgd_prob_vec)

  UL <- (B0 %*% S0) %*% L0 - L0
  LR <- (B1 %*% L1) - L1
  LL <- (BW %*% SW) %*% L0
  UR <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(R0, R1))

  rbind(cbind(UL, UR), cbind(LL, LR))
}

run_all_sims <- function(run_params) {
  all_results_list <- list()
  passage_times <- list()

  init_P_2N <- x$P[x$passage == 0 & x$ploidy == "2N"]
  init_P_4N <- x$P[x$passage == 0 & x$ploidy == "4N"]

  beta_buffer <- as.numeric(.first_non_null(run_params$beta_buffer, run_params$beta, 0.0))
  n_exp <- as.numeric(.first_non_null(run_params$n_exp, 1.0))
  smax <- as.numeric(.first_non_null(run_params$smax, 1.0))
  boundary_mode <- as.character(.first_non_null(run_params$boundary, "drop"))
  pwgd_val <- .get_pwgd(run_params)

  for (sim in sim_configs) {
    O2_LEVEL <- sim$O2
    pmis_const <- .pmisseg_of_O2(O2_LEVEL, run_params)

    x0_pre <- rep(0, R0)
    names(x0_pre) <- grid_pre
    x0_post <- rep(0, R1)
    names(x0_post) <- grid_post

    if (sim$init_layer == "pre") {
      init_P_values <- if (sim$init_ploidy == "2N") init_P_2N else init_P_4N
      x0_pre <- create_initial_dist(init_P_values, grid_pre, N_UNIT)
    } else {
      init_P_values <- if (sim$init_ploidy == "2N") init_P_2N else init_P_4N
      x0_post <- create_initial_dist(init_P_values, grid_post, N_UNIT)
    }

    x_current <- c(x0_pre, x0_post)
    x_current <- x_current / sum(x_current)

    lambda0_vec <- growth_lambda(
      O2_LEVEL, grid_pre,
      R = run_params$R,
      beta = run_params$beta,
      eta = run_params$eta,
      N_unit = N_UNIT,
      lam_min = run_params$lam_min,
      lam_max = run_params$lam_max,
      k_o = run_params$k_o
    )
    lambda1_vec <- growth_lambda(
      O2_LEVEL, grid_post,
      R = run_params$R,
      beta = run_params$beta,
      eta = run_params$eta,
      N_unit = N_UNIT,
      lam_min = run_params$lam_min,
      lam_max = run_params$lam_max,
      k_o = run_params$k_o
    )

    G <- .build_G_with_WGD(
      N0min = N_MIN, N0max = N_MAX,
      lambda0_vec = lambda0_vec,
      p0_vec = pmis_const,
      wgd_prob_vec = pwgd_val,
      N1min = N_MIN, N1max = N_MAX,
      lambda1_vec = lambda1_vec,
      p1_vec = pmis_const,
      boundary = boundary_mode,
      N_unit = N_UNIT,
      beta_buffer = beta_buffer,
      n_exp = n_exp,
      smax = smax
    )

    sim_passage_times <- numeric(PASSAGES_TO_RUN)

    for (p in 1:PASSAGES_TO_RUN) {
      pop_start <- sum(x_current)
      pop_target <- pop_start * POP_GROWTH_FACTOR
      time_in_passage <- 0.0

      while (sum(x_current) < pop_target) {
        x_current <- step_dt(G, x_current, DT, 1L)
        time_in_passage <- time_in_passage + DT
        if (sum(x_current) < pop_start * 1e-3 || time_in_passage > 1000) {
          break
        }
      }
      sim_passage_times[p] <- time_in_passage

      if (p %in% REPORT_PASSAGES) {
        pop_total <- sum(x_current)
        dist_df <- data.frame(
          sim_id = sim$id,
          passage = p,
          layer = c(rep("pre", R0), rep("post", R1)),
          N = c(grid_pre, grid_post),
          fraction = x_current / pop_total
        )
        all_results_list[[length(all_results_list) + 1]] <- dist_df
      }
      x_current <- x_current / sum(x_current) * pop_start
    }

    passage_times[[sim$id]] <- data.frame(
      sim_id = sim$id,
      passage = 1:PASSAGES_TO_RUN,
      duration = sim_passage_times
    )
  }

  all_dists <- do.call(rbind, all_results_list)
  all_passage_times <- do.call(rbind, passage_times)

  list(all_dists = all_dists, all_passage_times = all_passage_times)
}

run_in_vivo_crowd <- function(run_params,
                              O2_schedule = list(c(t0 = 0, t1 = Inf, O2 = 1.0)),
                              T_end = 28, sample_days = c(0, 7, 14, 21, 28),
                              N_UNIT = 22L, DT = 0.1,
                              K = 1e9, crowding = c("logistic", "gompertz"),
                              grid_pre = get("grid_pre", inherits = TRUE),
                              grid_post = get("grid_post", inherits = TRUE),
                              init_state) {
  crowding <- match.arg(crowding)

  R0 <- length(grid_pre)
  R1 <- length(grid_post)
  N0min <- min(grid_pre)
  N0max <- max(grid_pre)
  N1min <- min(grid_post)
  N1max <- max(grid_post)

  stopifnot(length(init_state) == (R0 + R1))
  v <- as.numeric(init_state)

  beta_buffer <- as.numeric(.first_non_null(run_params$beta_buffer, run_params$beta, 0.0))
  n_exp <- as.numeric(.first_non_null(run_params$n_exp, 1.0))
  smax <- as.numeric(.first_non_null(run_params$smax, 1.0))
  boundary_mode <- as.character(.first_non_null(run_params$boundary, "drop"))
  pwgd_val <- .get_pwgd(run_params)
  o2_burden_feedback <- isTRUE(.first_non_null(run_params$o2_burden_feedback, TRUE))
  o2_min <- as.numeric(.first_non_null(run_params$o2_min, 0.0))
  K_O2 <- as.numeric(.first_non_null(run_params$K_O2, K))
  h_O2 <- as.numeric(.first_non_null(run_params$h_O2, 1.0))

  get_O2 <- function(t) {
    for (seg in O2_schedule) {
      if (t >= seg["t0"] && t < seg["t1"]) return(as.numeric(seg["O2"]))
    }
    as.numeric(O2_schedule[[length(O2_schedule)]]["O2"])
  }

  apply_O2_feedback <- function(O2_base, Ntot) {
    O2_base <- .clip01(as.numeric(O2_base))
    if (!o2_burden_feedback) return(O2_base)

    o2_floor <- .clip01(o2_min)
    K_O2_use <- max(K_O2, 1e-12)
    h_O2_use <- max(h_O2, 1e-8)
    O2_eff <- o2_floor + (O2_base - o2_floor) / (1 + (Ntot / K_O2_use)^h_O2_use)
    .clip01(O2_eff)
  }

  G_cache <- new.env(parent = emptyenv())
  build_G_for_O2 <- function(O2) {
    key <- sprintf("%.3f", O2)
    if (!exists(key, envir = G_cache)) {
      lambda0 <- growth_lambda(
        O2, grid_pre, R = run_params$R, beta = run_params$beta, eta = run_params$eta,
        N_unit = N_UNIT, lam_min = run_params$lam_min, lam_max = run_params$lam_max, k_o = run_params$k_o
      )
      lambda1 <- growth_lambda(
        O2, grid_post, R = run_params$R, beta = run_params$beta, eta = run_params$eta,
        N_unit = N_UNIT, lam_min = run_params$lam_min, lam_max = run_params$lam_max, k_o = run_params$k_o
      )
      p_mis <- .pmisseg_of_O2(O2, run_params)
      G <- .build_G_with_WGD(
        N0min, N0max, lambda0, p0_vec = p_mis, wgd_prob_vec = pwgd_val,
        N1min, N1max, lambda1, p1_vec = p_mis,
        boundary = boundary_mode, N_unit = N_UNIT,
        beta_buffer = beta_buffer, n_exp = n_exp, smax = smax
      )
      assign(key, G, envir = G_cache)
    }
    get(key, envir = G_cache)
  }

  crowd <- function(Ntot) {
    if (crowding == "logistic") return(max(0, 1 - Ntot / K))
    exp(-Ntot / K)
  }

  I <- Diagonal(n = length(v))
  times <- seq(0, T_end, by = DT)
  snapshots <- list()
  size_trace <- data.frame(day = 0, Ntot = sum(v))

  for (t in times) {
    if (t %in% sample_days) {
      snapshots[[as.character(t)]] <- data.frame(
        day = t,
        layer = c(rep("pre", R0), rep("post", R1)),
        N = c(grid_pre, grid_post),
        fraction = v / sum(v),
        pop = sum(v)
      )
    }
    if (t >= T_end) break
    Ntot <- sum(v)
    O2t <- apply_O2_feedback(get_O2(t), Ntot)
    G <- build_G_for_O2(O2t)
    cfac <- crowd(Ntot)
    v <- as.numeric((I + DT * (cfac * G)) %*% v)
    size_trace <- rbind(size_trace, data.frame(day = t + DT, Ntot = sum(v)))
    if (sum(v) <= 1e-9) break
  }

  list(
    all_dists = do.call(rbind, snapshots),
    tumor_size = size_trace
  )
}

plot_misseg_interp <- function(par, o2_ref = 20.5) {
  O2 <- seq(0, 1, length.out = 401)
  if (!is.null(par$p_misseg)) {
    p0 <- as.numeric(par$p_misseg)
    k_o_mis <- as.numeric(.first_non_null(par$k_o_mis, 0.5))
    k_o_mis <- max(k_o_mis, 1e-12)
    p <- p0 * (1 - O2 / (O2 + k_o_mis))
  } else {
    p0 <- as.numeric(.first_non_null(par$pmis_O2_0, par$pmis, 1e-4))
    p1 <- as.numeric(.first_non_null(par$pmis_O2_1, par$pmis, p0))
    if (p0 <= 0 || p1 <= 0) {
      p <- rep(0, length(O2))
    } else {
      p <- exp((1 - O2) * log(p0) + O2 * log(p1))
    }
  }
  df <- data.frame(O2 = O2, O2_pct = O2 * o2_ref, p = pmax(p, 0))
  ggplot(df, aes(O2_pct, p)) +
    geom_line(linewidth = 1, color = "black") +
    geom_point(data = df[c(1, nrow(df)), ], size = 2, color = "red") +
    labs(x = "Oxygen (%)", y = "Missegregation rate") +
    theme_bw()
}
