suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

# ----------------------------------------------------------------------------
# Align miningcloneid oxygen model to Richard's buffering.R logic.
# model_O2_invivo extension:
# - Keep karyotype dynamics identical to Richard-aligned model.
# - Replace burden->O2 feedback with a windowed angiogenesis term:
#   baseline decline + mid-burden recovery + late decline.
# - Self-contained model file (no dependency on model_functions_ploidy_buffer.R)
# - Core dynamics (lambda, misseg delta, B/G construction) are defined here
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

.init_cpp_o2invivo_backend <- local({
  initialized <- FALSE
  available <- FALSE

  acquire_dir_lock <- function(lock_dir, timeout_sec = 300, poll_sec = 0.1) {
    start <- Sys.time()
    repeat {
      if (dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)) {
        return(TRUE)
      }

      # Recover from stale lock directories left by crashed workers.
      if (dir.exists(lock_dir)) {
        info <- suppressWarnings(file.info(lock_dir))
        mtime <- info$mtime[[1]]
        if (is.finite(as.numeric(mtime))) {
          age <- as.numeric(difftime(Sys.time(), mtime, units = "secs"))
          if (is.finite(age) && age > timeout_sec) {
            unlink(lock_dir, recursive = TRUE, force = TRUE)
            Sys.sleep(0.02)
            next
          }
        }
      }

      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      if (is.finite(elapsed) && elapsed >= timeout_sec) return(FALSE)
      Sys.sleep(poll_sec + stats::runif(1L, min = 0, max = 0.05))
    }
  }

  function() {
    if (initialized) return(available)
    initialized <<- TRUE

    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      stop("Rcpp package is required for model_O2_invivo but is not installed.")
    }

    # Dedicated O2 invivo backend (do not use Richard/shared cpp here).
    cpp_path <- file.path(.ALIGN_MODEL_DIR, "model_O2_invivo.cpp")
    if (!file.exists(cpp_path)) {
      stop("Cannot find required C++ backend file: ", cpp_path)
    }

    cache_root <- file.path(.ALIGN_MODEL_DIR, ".rcpp_cache_o2_invivo")
    cache_dir <- file.path(cache_root, "shared")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    lock_dir <- file.path(cache_root, ".sourcecpp_lock")
    lock_ok <- acquire_dir_lock(lock_dir, timeout_sec = 300, poll_sec = 0.1)
    if (!isTRUE(lock_ok)) {
      stop("Timed out waiting for sourceCpp lock: ", lock_dir)
    }
    on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)

    tryCatch({
      Rcpp::sourceCpp(
        file = cpp_path,
        rebuild = FALSE,
        showOutput = FALSE,
        verbose = FALSE,
        cacheDir = cache_dir
      )
    }, error = function(e) {
      stop("Failed to compile/load model_O2_invivo.cpp: ", conditionMessage(e))
    })

    required_fns <- c(
      "cpp_o2invivo_pr_delta_vec",
      "cpp_o2invivo_build_B_total_triplet",
      "cpp_o2invivo_build_B_WGD_triplet",
      "cpp_o2invivo_o2_window_supply",
      "cpp_o2invivo_build_G_for_o2_triplet",
      "cpp_o2invivo_simulate_one"
    )
    missing_fns <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function", inherits = TRUE)]
    if (length(missing_fns) > 0L) {
      stop(
        "model_O2_invivo C++ backend loaded but required symbols are missing: ",
        paste(missing_fns, collapse = ", ")
      )
    }

    available <<- TRUE

    available
  }
})

.USE_CPP_O2INVIVO_BACKEND <- .init_cpp_o2invivo_backend()

.first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

# Euler stepper for generator matrix dynamics.
step_dt <- function(G, x, dt, steps = 1L, normalize = FALSE) {
  I <- Diagonal(n = nrow(G))
  A <- I + dt * G
  v <- as.numeric(x)
  for (k in seq_len(steps)) {
    v <- as.numeric(A %*% v)
    if (normalize) v <- v / sum(v)
  }
  v
}

# Build a normalized initial distribution on an integer N grid from ploidy values.
create_initial_dist <- function(ploidy_values, N_grid, N_unit = 22L) {
  N_values <- round(as.numeric(ploidy_values) * as.numeric(N_unit))
  N_counts <- table(N_values)
  N_fracs <- as.numeric(N_counts) / sum(N_counts)
  names(N_fracs) <- names(N_counts)
  x_vec <- rep(0, length(N_grid))
  names(x_vec) <- N_grid
  valid_names <- names(N_fracs)[names(N_fracs) %in% names(x_vec)]
  x_vec[valid_names] <- N_fracs[valid_names]
  x_vec
}

# Build explicit initial state vectors.
# - Either pass ploidy + layer, or custom init_pre/init_post fractions.
# - Returns absolute cell counts with total mass = total_size.
make_init_state <- function(grid_pre, grid_post,
                            ploidy = c(2, 4), layer = c("pre", "post"),
                            init_pre = NULL, init_post = NULL,
                            N_UNIT = 22L, total_size = 1e6) {
  layer <- match.arg(layer)
  ploidy <- match.arg(as.character(ploidy), choices = c("2", "4"))
  Pnum <- as.numeric(ploidy)

  R0 <- length(grid_pre)
  R1 <- length(grid_post)
  x_pre <- rep(0, R0)
  names(x_pre) <- grid_pre
  x_post <- rep(0, R1)
  names(x_post) <- grid_post

  if (!is.null(init_pre) || !is.null(init_post)) {
    if (!is.null(init_pre)) x_pre[names(init_pre)] <- as.numeric(init_pre)
    if (!is.null(init_post)) x_post[names(init_post)] <- as.numeric(init_post)
  } else {
    N_delta <- as.integer(Pnum * N_UNIT)
    if (layer == "pre") {
      stopifnot(N_delta %in% grid_pre)
      x_pre[as.character(N_delta)] <- 1
    } else {
      stopifnot(N_delta %in% grid_post)
      x_post[as.character(N_delta)] <- 1
    }
  }

  s <- sum(x_pre) + sum(x_post)
  if (s <= 0) stop("Init mass is zero.")
  c(x_pre, x_post) / s * total_size
}

.clip01 <- function(x) pmin(pmax(x, 0), 1)
.clip_o2pct <- function(x) pmin(pmax(x, 0), 100)
.sigmoid01 <- function(z) 1 / (1 + exp(-z))
.assert_o2_pct <- function(x, label = "O2") {
  x_num <- as.numeric(x)
  if (any(!is.finite(x_num))) stop(label, " must be finite.")
  if (any(x_num < 0 | x_num > 100)) stop(label, " must be in percent scale [0, 100].")
  x_num
}

.require_cpp_o2invivo_fn <- function(fn_name) {
  if (!(isTRUE(.USE_CPP_O2INVIVO_BACKEND) &&
        exists(fn_name, mode = "function", inherits = TRUE))) {
    stop("Required C++ backend function is unavailable: ", fn_name)
  }
}

# O2(N): baseline decline plus transient angiogenesis window on log10(N).
# O2_eff = o2_min + (O2_base-o2_min)/(1 + (N/K_down)^h_down) + A_ang * W(log10(N))
# W(x) = sig((x-m_on)/s_on) * (1 - sig((x-m_off)/s_off))
.o2_window_supply_from_burden <- function(Ntot,
                                          run_params,
                                          O2_base = 5.0,
                                          o2_min = 0.0,
                                          h_down = 1.0,
                                          o2_logN_eps = 1.0) {
  Ntot_use <- pmax(as.numeric(Ntot), 0)
  O2_base_use <- .assert_o2_pct(O2_base, label = "O2_base")
  o2_floor <- .assert_o2_pct(o2_min, label = "o2_min")
  h_use <- as.numeric(h_down)
  if (!is.finite(h_use) || h_use <= 0) h_use <- 1.0

  K_down <- as.numeric(.first_non_null(run_params$K_down, 1e12))
  if (!is.finite(K_down) || K_down <= 0) K_down <- 1e12
  A_ang <- .clip_o2pct(as.numeric(.first_non_null(run_params$A_ang, 0.0)))
  m_on <- as.numeric(.first_non_null(run_params$m_on, 9.0))
  delta_m <- as.numeric(.first_non_null(run_params$delta_m, 1.0))
  if (!is.finite(delta_m) || delta_m <= 0) delta_m <- 1.0
  m_off <- as.numeric(.first_non_null(run_params$m_off, m_on + delta_m))
  if (!is.finite(m_off) || m_off <= m_on) m_off <- m_on + delta_m
  s_on <- as.numeric(.first_non_null(run_params$s_on, 0.3))
  s_off <- as.numeric(.first_non_null(run_params$s_off, 0.3))
  if (!is.finite(s_on) || s_on <= 0) s_on <- 0.3
  if (!is.finite(s_off) || s_off <= 0) s_off <- 0.3

  eps_use <- as.numeric(o2_logN_eps)
  if (!is.finite(eps_use) || eps_use <= 0) eps_use <- 1.0

  .require_cpp_o2invivo_fn("cpp_o2invivo_o2_window_supply")
  return(as.numeric(cpp_o2invivo_o2_window_supply(
    Ntot = as.numeric(Ntot_use),
    O2_base = as.numeric(O2_base_use),
    o2_min = as.numeric(o2_floor),
    h_down = as.numeric(h_use),
    K_down = as.numeric(K_down),
    A_ang = as.numeric(A_ang),
    m_on = as.numeric(m_on),
    m_off = as.numeric(m_off),
    s_on = as.numeric(s_on),
    s_off = as.numeric(s_off),
    o2_logN_eps = as.numeric(eps_use)
  )))
}

.get_pwgd <- function(run_params) {
  as.numeric(.first_non_null(run_params$p_wgd, 0))
}

# Richard buffering.R-style lambda: O2-only saturating rate.
growth_lambda <- function(O2, N, lam_min, lam_max, k_o) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  lam_min_use <- as.numeric(lam_min)
  lam_max_use <- as.numeric(lam_max)
  k_o_use <- as.numeric(k_o)
  if (!is.finite(lam_min_use)) stop("lam_min must be finite.")
  if (!is.finite(lam_max_use)) stop("lam_max must be finite.")
  if (!is.finite(k_o_use) || k_o_use <= 0) stop("k_o must be > 0.")
  k_o_use <- max(k_o_use, 1e-12)

  frac <- O2_use / (O2_use + k_o_use)
  lam <- lam_min_use + (lam_max_use - lam_min_use) * frac
  rep(pmax(lam, 0), length(N))
}

# Richard buffering.R-style O2-dependent missegregation.
# Endpoint interpolation branch uses O2 in percent scale [0, 100].
.pmisseg_of_O2 <- function(O2, run_params) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  p0 <- as.numeric(run_params$p_misseg)
  k_o_mis <- as.numeric(run_params$k_o_mis)
  if (!is.finite(p0) || p0 < 0) stop("run_params$p_misseg must be finite and >= 0.")
  if (!is.finite(k_o_mis) || k_o_mis <= 0) stop("run_params$k_o_mis must be > 0.")
  p <- p0 * (1 - (O2_use / (O2_use + k_o_mis)))
  .clip01(p)
}

# Richard buffering.R delta weight formula.
.pr_delta_vec <- function(N, p, eps_tail = 1e-8, mr_lethality = 0.9,
                          beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L) {
  .require_cpp_o2invivo_fn("cpp_o2invivo_pr_delta_vec")
  res <- cpp_o2invivo_pr_delta_vec(
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

.build_B_total <- function(Nmin, Nmax, p_vec, mr_lethality = 0.9,
                           boundary = c("drop", "absorb_minmax"),
                           eps_tail = 1e-8, return_sparse = TRUE,
                           beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L) {
  boundary <- match.arg(boundary)
  R <- Nmax - Nmin + 1L
  if (length(p_vec) == 1L) p_vec <- rep(p_vec, R)

  .require_cpp_o2invivo_fn("cpp_o2invivo_build_B_total_triplet")
  tri <- cpp_o2invivo_build_B_total_triplet(
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

.build_B_WGD <- function(N0min, N0max, N1min, N1max,
                         boundary = c("drop", "absorb_minmax"),
                         return_sparse = TRUE) {
  boundary <- match.arg(boundary)
  R0 <- N0max - N0min + 1L
  R1 <- N1max - N1min + 1L

  .require_cpp_o2invivo_fn("cpp_o2invivo_build_B_WGD_triplet")
  tri <- cpp_o2invivo_build_B_WGD_triplet(
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

  beta_buffer <- as.numeric(.first_non_null(run_params$beta_buffer, 0.0))
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
      lam_min = run_params$lam_min,
      lam_max = run_params$lam_max,
      k_o = run_params$k_o
    )
    lambda1_vec <- growth_lambda(
      O2_LEVEL, grid_post,
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
                              O2_schedule = list(c(t0 = 0, t1 = Inf, O2 = 5.0)),
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

  beta_buffer <- as.numeric(.first_non_null(run_params$beta_buffer, 0.0))
  n_exp <- as.numeric(.first_non_null(run_params$n_exp, 1.0))
  smax <- as.numeric(.first_non_null(run_params$smax, 1.0))
  boundary_mode <- as.character(.first_non_null(run_params$boundary, "drop"))
  pwgd_val <- .get_pwgd(run_params)
  o2_burden_feedback <- isTRUE(.first_non_null(run_params$o2_burden_feedback, TRUE))
  o2_min <- .assert_o2_pct(as.numeric(.first_non_null(run_params$o2_min, 0.0)), label = "o2_min")
  h_O2 <- as.numeric(.first_non_null(run_params$h_down, 1.0))

  get_O2 <- function(t) {
    for (seg in O2_schedule) {
      if (t >= seg["t0"] && t < seg["t1"]) return(as.numeric(seg["O2"]))
    }
    as.numeric(O2_schedule[[length(O2_schedule)]]["O2"])
  }

  apply_O2_feedback <- function(O2_base, Ntot) {
    O2_base <- .assert_o2_pct(as.numeric(O2_base), label = "O2_schedule value")
    if (!o2_burden_feedback) return(O2_base)
    .o2_window_supply_from_burden(
      Ntot = Ntot,
      run_params = run_params,
      O2_base = O2_base,
      o2_min = o2_min,
      h_down = h_O2,
      o2_logN_eps = 1.0
    )
  }

  G_cache <- new.env(parent = emptyenv())
  build_G_for_O2 <- function(O2) {
    O2_use <- .assert_o2_pct(as.numeric(O2), label = "O2")
    key <- sprintf("%.3f", O2_use)
    if (!exists(key, envir = G_cache, inherits = FALSE)) {
      .require_cpp_o2invivo_fn("cpp_o2invivo_build_G_for_o2_triplet")

      lam_min_use <- as.numeric(.first_non_null(run_params$lam_min, 1.0))
      lam_max_use <- as.numeric(.first_non_null(run_params$lam_max, lam_min_use))
      k_o_use <- as.numeric(.first_non_null(run_params$k_o, 50.0))
      has_p_misseg <- !is.null(run_params$p_misseg)

      tri <- cpp_o2invivo_build_G_for_o2_triplet(
        O2 = as.numeric(O2_use),
        N0min = as.integer(N0min),
        N0max = as.integer(N0max),
        N1min = as.integer(N1min),
        N1max = as.integer(N1max),
        lam_min = as.numeric(lam_min_use),
        lam_max = as.numeric(lam_max_use),
        k_o = as.numeric(k_o_use),
        has_p_misseg = isTRUE(has_p_misseg),
        p_misseg = as.numeric(.first_non_null(run_params$p_misseg, 0.0)),
        k_o_mis = as.numeric(.first_non_null(run_params$k_o_mis, 50.0)),
        has_pmis_endpoints = FALSE,
        pmis_O2_0 = 0.0,
        pmis_O2_1 = 0.0,
        p_const = 0.0,
        p_wgd = as.numeric(pwgd_val),
        boundary = as.character(boundary_mode),
        eps_tail = as.numeric(1e-8),
        beta_buffer = as.numeric(beta_buffer),
        n_exp = as.numeric(n_exp),
        smax = as.numeric(smax),
        N_unit = as.integer(N_UNIT)
      )
      G <- sparseMatrix(
        i = as.integer(tri$i),
        j = as.integer(tri$j),
        x = as.numeric(tri$x),
        dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
        repr = "C"
      )
      assign(key, G, envir = G_cache)
    }
    get(key, envir = G_cache, inherits = FALSE)
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
  O2 <- seq(0, 100, length.out = 401)
  p0 <- as.numeric(par$p_misseg)
  k_o_mis <- as.numeric(.first_non_null(par$k_o_mis, 50.0))
  if (!is.finite(p0) || p0 < 0) stop("par$p_misseg must be finite and >= 0.")
  if (!is.finite(k_o_mis) || k_o_mis <= 0) stop("par$k_o_mis must be > 0.")
  p <- p0 * (1 - O2 / (O2 + k_o_mis))
  df <- data.frame(O2_pct = O2, p = pmax(p, 0))
  ggplot(df, aes(O2_pct, p)) +
    geom_line(linewidth = 1, color = "black") +
    geom_point(data = df[c(1, nrow(df)), ], size = 2, color = "red") +
    labs(x = "Oxygen (%)", y = "Missegregation rate") +
    theme_bw()
}
