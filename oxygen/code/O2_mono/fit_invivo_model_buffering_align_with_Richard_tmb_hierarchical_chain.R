#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

get_script_dir_self <- function() {
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
  if (length(farg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", farg[[1]]), mustWork = FALSE))
}

script_dir <- get_script_dir_self()
model_script <- file.path(script_dir, "model_buffering_align_with_Richard.R")
fit_script <- file.path(script_dir, "fit_invivo_model_buffering_align_with_Richard.R")
if (!file.exists(model_script)) stop("Cannot find model script: ", model_script)
if (!file.exists(fit_script)) stop("Cannot find fit script: ", fit_script)

Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = script_dir)
source(model_script)
source(fit_script)

as_char_vec <- function(x, default = character(0)) {
  if (is.null(x)) return(default)
  s <- trimws(as.character(x))
  if (!nzchar(s)) return(default)
  parts <- trimws(unlist(strsplit(s, "[,;]", perl = TRUE)))
  parts[nzchar(parts)]
}

sanitize_tag <- function(x) {
  s <- trimws(as.character(x %||% "NA"))
  s <- gsub(" ", "", s, fixed = TRUE)
  s <- gsub(",", "", s, fixed = TRUE)
  s <- gsub("\\.", "p", s)
  s <- gsub("-", "m", s, fixed = TRUE)
  s
}

weighted_median <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  keep <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[[which(cw >= 0.5)[1]]]
}

weighted_mad <- function(x, w, center = NULL) {
  if (!is.finite(center)) center <- weighted_median(x, w)
  ad <- abs(as.numeric(x) - center)
  m <- weighted_median(ad, w)
  if (!is.finite(m)) return(NA_real_)
  1.4826 * m
}

robust_pool_local <- function(local_mat, sample_weight, tau_floor = 0.05) {
  stopifnot(is.matrix(local_mat))
  p <- ncol(local_mat)
  mu <- numeric(p)
  tau <- numeric(p)
  names(mu) <- colnames(local_mat)
  names(tau) <- colnames(local_mat)
  for (j in seq_len(p)) {
    x <- local_mat[, j]
    med <- weighted_median(x, sample_weight)
    if (!is.finite(med)) med <- stats::median(x, na.rm = TRUE)
    md <- weighted_mad(x, sample_weight, center = med)
    if (!is.finite(md) || md < tau_floor) md <- tau_floor
    mu[[j]] <- med
    tau[[j]] <- md
  }
  list(mu = mu, tau = tau)
}

.deoptim_optimize_box <- function(objective_fn,
                                  lower,
                                  upper,
                                  itermax = 120L,
                                  NP = 80L,
                                  init_vec = NULL,
                                  alt_init = NULL,
                                  trace = TRUE,
                                  seed = NULL,
                                  label = "deoptim",
                                  deoptim_parallel = FALSE,
                                  n_cores = 1L,
                                  script_path_for_workers = NULL,
                                  cluster_port = NULL,
                                  strict_parallel = FALSE) {
  if (!requireNamespace("DEoptim", quietly = TRUE)) {
    stop("[", label, "] DEoptim requested but package DEoptim is not available.")
  }
  if (!is.null(seed) && is.finite(seed)) set.seed(as.integer(seed))

  nm <- names(lower)
  if (is.null(nm) || any(!nzchar(nm))) {
    nm <- names(upper)
  }
  if (is.null(nm) || length(nm) != length(lower) || any(!nzchar(nm))) {
    nm <- paste0("p", seq_along(lower))
  }
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  names(lower) <- nm
  names(upper) <- nm

  itermax_use <- as.integer(max(1L, itermax))
  NP_use <- as.integer(max(4L, NP))
  n_cores_use <- as.integer(max(1L, n_cores))
  do_parallel <- isTRUE(deoptim_parallel) && n_cores_use > 1L
  requested_workers <- if (do_parallel) as.integer(max(1L, min(n_cores_use, NP_use))) else 1L
  started_workers <- if (do_parallel) NA_integer_ else 1L
  active_workers <- 1L
  cluster_port_used <- NA_integer_
  fallback <- FALSE
  fallback_reason <- ""

  fmt_int <- function(x) {
    if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
    as.character(as.integer(x))
  }
  cluster_workers <- function(cl) {
    n <- tryCatch(length(cl), error = function(e) NA_integer_)
    if (is.na(n) || !is.finite(n)) return(NA_integer_)
    as.integer(max(0L, n))
  }

  ctrl <- list(
    trace = isTRUE(trace),
    itermax = itermax_use,
    NP = NP_use,
    strategy = 2
  )

  if (!is.null(init_vec)) {
    align_init <- function(v) {
      if (is.null(v)) return(rep(NA_real_, length(lower)))
      vv <- as.numeric(v)
      names(vv) <- names(v) %||% rep("", length(vv))
      out <- rep(NA_real_, length(lower))
      names(out) <- names(lower)
      if (!is.null(names(v)) && any(nzchar(names(v)))) {
        idx <- match(names(lower), names(v))
        good <- is.finite(idx)
        out[good] <- vv[idx[good]]
      } else {
        ncopy <- min(length(out), length(vv))
        out[seq_len(ncopy)] <- vv[seq_len(ncopy)]
      }
      out[!is.finite(out)] <- (lower[!is.finite(out)] + upper[!is.finite(out)]) / 2
      out
    }

    d <- length(lower)
    initpop <- matrix(
      stats::runif(NP_use * d, min = lower, max = upper),
      nrow = NP_use,
      byrow = TRUE
    )
    init1 <- clip(align_init(init_vec), lower, upper)
    initpop[1, ] <- init1
    if (!is.null(alt_init) && NP_use >= 2L) {
      init2 <- clip(align_init(alt_init), lower, upper)
      initpop[2, ] <- init2
    }
    colnames(initpop) <- names(lower)
    ctrl$initialpop <- initpop
  }

  if (do_parallel) {
    message("[", label, "] DEoptim parallel requested with n_cores=", requested_workers, ".")
    seed_use <- if (!is.null(seed) && is.finite(seed)) abs(as.integer(seed)) else 0L
    label_hash <- if (!is.null(label) && nzchar(as.character(label))) {
      sum(utf8ToInt(as.character(label)))
    } else {
      0L
    }
    pid_use <- abs(as.integer(Sys.getpid()))
    calc_port <- function(k) {
      raw <- (as.numeric(pid_use) * 131) + (as.numeric(seed_use) * 17) + as.numeric(label_hash) + (as.numeric(k) * 977)
      as.integer(11000L + (raw %% 50000L))
    }
    port_candidates <- integer(0)
    if (!is.null(cluster_port) && is.finite(cluster_port)) {
      port_candidates <- c(port_candidates, as.integer(cluster_port))
    }
    for (k in seq_len(8L)) {
      port_candidates <- c(port_candidates, calc_port(k))
    }
    port_candidates <- unique(port_candidates[is.finite(port_candidates)])
    if (length(port_candidates) == 0L) {
      port_candidates <- 11000L + seq_len(8L)
    }

    cl <- NULL
    cluster_start_last_err <- ""
    for (p in port_candidates) {
      cl <- tryCatch(
        parallel::makePSOCKcluster(requested_workers, port = as.integer(p)),
        error = function(e) {
          cluster_start_last_err <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(cl)) {
        cluster_port_used <- as.integer(p)
        break
      }
    }
    if (is.null(cl)) {
      message(
        "[", label, "] Could not start parallel workers for DEoptim after ",
        length(port_candidates), " port attempts",
        if (nzchar(cluster_start_last_err)) paste0(": ", cluster_start_last_err) else "."
      )
    }
    if (!is.null(cl)) {
      started_workers <- cluster_workers(cl)
      init_ok <- tryCatch(
        {
          if (!is.null(script_path_for_workers) && nzchar(script_path_for_workers) && file.exists(script_path_for_workers)) {
            parallel::clusterCall(cl, function(path) {
              Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(path))
              source(path)
              NULL
            }, script_path_for_workers)
          }
          parallel::clusterExport(cl, varlist = c("objective_fn"), envir = environment())
          TRUE
        },
        error = function(e) {
          message("[", label, "] Failed to initialize DEoptim workers: ", conditionMessage(e))
          FALSE
        }
      )
      if (isTRUE(init_ok)) {
        on.exit(parallel::stopCluster(cl), add = TRUE)
        ctrl$cluster <- cl
        active_workers <- ifelse(is.na(started_workers) || started_workers < 1L, 1L, started_workers)
        if (isTRUE(strict_parallel) && active_workers <= 1L) {
          try(parallel::stopCluster(cl), silent = TRUE)
          stop(
            "[", label, "] strict parallel mode: n_cores>1 but active_workers=", fmt_int(active_workers),
            ". Aborting instead of fallback."
          )
        }
      } else {
        try(parallel::stopCluster(cl), silent = TRUE)
        fallback <- TRUE
        fallback_reason <- "worker_init_failed"
        if (is.na(started_workers)) started_workers <- 0L
        active_workers <- 1L
        if (isTRUE(strict_parallel)) {
          stop(
            "[", label, "] strict parallel mode: worker initialization failed; ",
            "requested_workers=", fmt_int(requested_workers),
            ", started_workers=", fmt_int(started_workers),
            ". Aborting instead of fallback."
          )
        }
      }
    } else {
      started_workers <- 0L
      fallback <- TRUE
      fallback_reason <- "cluster_start_failed"
      active_workers <- 1L
      if (isTRUE(strict_parallel)) {
        stop(
          "[", label, "] strict parallel mode: cluster start failed; ",
          "requested_workers=", fmt_int(requested_workers),
          ", started_workers=0",
          if (nzchar(cluster_start_last_err)) paste0(", last_error=", cluster_start_last_err) else "",
          ". Aborting instead of fallback."
        )
      }
    }
  } else {
    message("[", label, "] DEoptim running in serial mode (n_cores=", n_cores_use, ").")
  }

  backend_msg <- paste0(
    "[", label, "] DEoptim backend: requested_workers=", fmt_int(requested_workers),
    ", started_workers=", fmt_int(started_workers),
    ", active_workers=", fmt_int(active_workers),
    ", fallback=", if (isTRUE(fallback)) "TRUE" else "FALSE",
    ", port=", fmt_int(cluster_port_used),
    ", NP=", NP_use,
    ", itermax=", itermax_use
  )
  if (isTRUE(fallback) && nzchar(fallback_reason)) {
    backend_msg <- paste0(backend_msg, ", reason=", fallback_reason)
  }
  message(backend_msg)

  # DEoptim workers may pass unnamed numeric vectors in parallel mode.
  # Normalize to named parameter vectors and harden objective against transient eval errors.
  fn_names <- names(lower)
  fn_fails <- 0L
  objective_safe <- function(x) {
    xx <- suppressWarnings(as.numeric(x))
    if (length(xx) != length(fn_names)) {
      fn_fails <<- fn_fails + 1L
      return(1e30)
    }
    names(xx) <- fn_names
    v <- tryCatch(
      objective_fn(xx),
      error = function(e) {
        fn_fails <<- fn_fails + 1L
        NA_real_
      }
    )
    if (!is.finite(v)) {
      fn_fails <<- fn_fails + 1L
      return(1e30)
    }
    as.numeric(v)
  }

  fit <- tryCatch(
    DEoptim::DEoptim(
      fn = objective_safe,
      lower = lower,
      upper = upper,
      control = ctrl
    ),
    error = function(e) {
      stop("[", label, "] DEoptim failed: ", conditionMessage(e))
    }
  )

  best <- clip(as.numeric(fit$optim$bestmem), lower, upper)
  names(best) <- names(lower)
  list(
    par = best,
    value = as.numeric(fit$optim$bestval),
    convergence = 0L,
    method = if (active_workers > 1L) "DEoptim_parallel" else "DEoptim_serial",
    itermax = itermax_use,
    NP = NP_use,
    parallel_info = list(
      requested_workers = requested_workers,
      started_workers = started_workers,
      active_workers = active_workers,
      cluster_port = cluster_port_used,
      fallback = fallback,
      fallback_reason = if (isTRUE(fallback)) fallback_reason else "",
      objective_fail_count = as.integer(fn_fails)
    )
  )
}

.init_tmb_hier_backend <- local({
  dll_name <- NULL
  function(script_dir, rebuild = FALSE) {
    if (!requireNamespace("TMB", quietly = TRUE)) {
      stop("TMB package is required but not available.")
    }

    cpp_base <- "model_buffering_align_with_Richard_tmb_hierarchical"
    cpp_file <- file.path(script_dir, paste0(cpp_base, ".cpp"))
    if (!file.exists(cpp_file)) stop("Missing TMB C++ file: ", cpp_file)

    dll_file <- TMB::dynlib(file.path(script_dir, cpp_base))
    if (isTRUE(rebuild) || !file.exists(dll_file)) {
      message("Compiling TMB model: ", basename(cpp_file))
      TMB::compile(cpp_file, safeunload = TRUE)
    }

    dll_name_now <- basename(tools::file_path_sans_ext(dll_file))
    if (!dll_name_now %in% names(getLoadedDLLs())) {
      dyn.load(dll_file)
    }
    dll_name <<- dll_name_now
    dll_name
  }
})

pool_locals_with_tmb <- function(local_hat_mat,
                                 sample_weight,
                                 dll_name,
                                 init_mu,
                                 init_tau,
                                 tau_min = 1e-3,
                                 log_tau_prior_sd = 2.0,
                                 maxit = 200) {
  p <- ncol(local_hat_mat)
  n <- nrow(local_hat_mat)
  if (length(init_mu) != p) stop("init_mu length mismatch")
  if (length(init_tau) != p) stop("init_tau length mismatch")

  init_log_tau <- log(pmax(as.numeric(init_tau) - tau_min, 1e-6))
  names(init_log_tau) <- colnames(local_hat_mat)

  data <- list(
    theta_hat = unname(as.matrix(local_hat_mat)),
    sample_weight = as.numeric(sample_weight),
    tau_min = as.numeric(tau_min),
    log_tau_prior_sd = as.numeric(log_tau_prior_sd)
  )
  parameters <- list(
    mu = as.numeric(init_mu[colnames(local_hat_mat)]),
    log_tau = as.numeric(init_log_tau[colnames(local_hat_mat)]),
    theta_shrunk = unname(as.matrix(local_hat_mat))
  )

  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    DLL = dll_name,
    silent = TRUE
  )

  fallback <- function(msg = "tmb_pool_failed") {
    mu_fb <- as.numeric(init_mu[colnames(local_hat_mat)])
    tau_fb <- pmax(as.numeric(init_tau[colnames(local_hat_mat)]), tau_min)
    theta_fb <- as.matrix(local_hat_mat)
    rownames(theta_fb) <- rownames(local_hat_mat)
    colnames(theta_fb) <- colnames(local_hat_mat)
    names(mu_fb) <- colnames(local_hat_mat)
    names(tau_fb) <- colnames(local_hat_mat)
    list(
      mu = mu_fb,
      tau = tau_fb,
      theta_shrunk = theta_fb,
      objective = NA_real_,
      convergence = 999L,
      message = msg
    )
  }

  opt <- tryCatch(
    nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      control = list(iter.max = as.integer(maxit), eval.max = as.integer(maxit) * 2L)
    ),
    error = function(e) {
      warning("TMB pooling failed; fallback to robust pool: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(opt) || !is.finite(opt$objective)) {
    return(fallback("nlminb_failed"))
  }

  rep <- tryCatch(obj$report(), error = function(e) NULL)
  if (is.null(rep)) {
    warning("TMB report failed; fallback to robust pool.")
    return(fallback("report_failed"))
  }

  mu <- as.numeric(rep$mu)
  tau <- as.numeric(rep$tau)
  theta_shrunk <- as.matrix(rep$theta_shrunk)
  ok <- all(is.finite(mu)) && all(is.finite(tau)) && all(is.finite(theta_shrunk))
  if (!ok) {
    warning("TMB pooling produced non-finite values; fallback to robust pool.")
    return(fallback("nonfinite_report"))
  }

  rownames(theta_shrunk) <- rownames(local_hat_mat)
  colnames(theta_shrunk) <- colnames(local_hat_mat)
  names(mu) <- colnames(local_hat_mat)
  names(tau) <- colnames(local_hat_mat)

  list(
    mu = mu,
    tau = tau,
    theta_shrunk = theta_shrunk,
    objective = opt$objective,
    convergence = opt$convergence,
    message = opt$message
  )
}

compose_full_t <- function(global_t, local_t, full_names, global_names, local_names) {
  out <- numeric(length(full_names))
  names(out) <- full_names
  out[global_names] <- as.numeric(global_t[global_names])
  out[local_names] <- as.numeric(local_t[local_names])
  out
}

scenario_id <- function(sc) {
  paste(sc$harvest, sc$cohort, sc$dose, sep = "__")
}

scenario_data_weight <- function(sc) {
  n_b <- length(sc$obs_days)
  n_p <- length(sc$ploidy_obs_N)
  w <- n_b + n_p
  if (!is.finite(w) || w <= 0) w <- 1
  as.numeric(w)
}

.eval_scenario_components <- function(full_t, sc, cfg_step) {
  raw <- evaluate_objective_components_raw(full_t, scenarios = list(sc), cfg = cfg_step)

  scale_b <- if (isTRUE(cfg_step$loss_rescale)) cfg_step$loss_scale_burden else 1.0
  scale_p <- if (isTRUE(cfg_step$loss_rescale)) cfg_step$loss_scale_ploidy else 1.0
  if (!is.finite(scale_b) || scale_b <= 0) scale_b <- 1.0
  if (!is.finite(scale_p) || scale_p <= 0) scale_p <- 1.0

  L_b_scaled <- raw$L_b / scale_b
  L_p_scaled <- raw$L_p / scale_p
  L_data <- cfg_step$w_burden * L_b_scaled + cfg_step$w_ploidy * L_p_scaled

  list(
    L_b = raw$L_b,
    L_p = raw$L_p,
    L_b_scaled = L_b_scaled,
    L_p_scaled = L_p_scaled,
    L_data = L_data
  )
}

fit_local_single <- function(sc,
                             global_t,
                             local_start,
                             mu_local,
                             tau_local,
                             cfg_step,
                             full_names,
                             global_names,
                             local_names,
                             local_lower,
                             local_upper,
                             lambda_shrink = 1.0,
                             tau_floor = 0.05,
                             use_deoptim = TRUE,
                             deoptim_itermax = 120L,
                             deoptim_np = 80L,
                             deoptim_trace = TRUE,
                             inner_cores = 1L,
                             script_path_for_workers = NULL,
                             cluster_port = NULL,
                             n_starts = 6,
                             maxit = 2500,
                             seed = 1) {
  set.seed(seed)
  local_start <- clip(local_start, local_lower, local_upper)
  mu_local <- clip(mu_local, local_lower, local_upper)
  tau_use <- pmax(as.numeric(tau_local), tau_floor)

  penalty_fn <- function(par_local) {
    if (lambda_shrink <= 0) return(0)
    z <- (par_local - mu_local) / tau_use
    lambda_shrink * sum(z * z)
  }

  objective_local <- function(par_local) {
    full_t <- compose_full_t(global_t, par_local, full_names, global_names, local_names)
    comp <- .eval_scenario_components(full_t, sc, cfg_step)
    comp$L_data + penalty_fn(par_local)
  }

  if (isTRUE(use_deoptim)) {
    de_fit <- .deoptim_optimize_box(
      objective_fn = objective_local,
      lower = local_lower,
      upper = local_upper,
      itermax = deoptim_itermax,
      NP = deoptim_np,
      init_vec = local_start,
      alt_init = mu_local,
      trace = deoptim_trace,
      seed = seed,
      label = paste0("local:", scenario_id(sc)),
      deoptim_parallel = as.integer(max(1L, inner_cores)) > 1L,
      n_cores = as.integer(max(1L, inner_cores)),
      script_path_for_workers = script_path_for_workers,
      cluster_port = cluster_port,
      strict_parallel = as.integer(max(1L, inner_cores)) > 1L
    )
    best_par <- de_fit$par
    best_conv <- de_fit$convergence
    best_method <- de_fit$method
  } else {
    starts <- vector("list", n_starts)
    starts[[1]] <- local_start
    if (n_starts >= 2) starts[[2]] <- mu_local
    if (n_starts >= 3) {
      for (k in 3:n_starts) {
        starts[[k]] <- stats::runif(length(local_lower), min = local_lower, max = local_upper)
        names(starts[[k]]) <- names(local_lower)
      }
    }

    best <- list(par = local_start, value = Inf, conv = 999)
    for (k in seq_along(starts)) {
      fit <- tryCatch(
        optim(
          par = starts[[k]],
          fn = objective_local,
          method = "L-BFGS-B",
          lower = local_lower,
          upper = local_upper,
          control = list(maxit = as.integer(maxit))
        ),
        error = function(e) NULL
      )
      if (!is.null(fit) && is.finite(fit$value) && fit$value < best$value) {
        best <- list(par = fit$par, value = fit$value, conv = fit$convergence)
      }
    }
    best_par <- clip(best$par, local_lower, local_upper)
    best_conv <- best$conv
    best_method <- "optim_multistart"
  }

  full_best <- compose_full_t(global_t, best_par, full_names, global_names, local_names)
  comp_best <- .eval_scenario_components(full_best, sc, cfg_step)
  pen_best <- penalty_fn(best_par)

  list(
    local_t = best_par,
    objective_total = comp_best$L_data + pen_best,
    objective_data = comp_best$L_data,
    objective_penalty = pen_best,
    L_b = comp_best$L_b,
    L_p = comp_best$L_p,
    L_b_scaled = comp_best$L_b_scaled,
    L_p_scaled = comp_best$L_p_scaled,
    convergence = best_conv,
    solver = best_method
  )
}

fit_locals_parallel <- function(scenarios,
                                global_t,
                                local_mat_start,
                                mu_local,
                                tau_local,
                                cfg_step,
                                full_names,
                                global_names,
                                local_names,
                                local_lower,
                                local_upper,
                                lambda_shrink,
                                tau_floor,
                                n_starts_local,
                                maxit_local,
                                use_deoptim_local,
                                deoptim_itermax_local,
                                deoptim_np_local,
                                deoptim_trace,
                                n_cores,
                                script_path_for_workers = NULL,
                                seed_base = 1) {
  n <- length(scenarios)
  n_cores_use <- as.integer(max(1L, n_cores))
  workers <- as.integer(max(1L, min(n_cores_use, n)))
  inner_alloc <- rep(1L, n)
  if (n_cores_use >= n) {
    base <- as.integer(n_cores_use %/% n)
    rem <- as.integer(n_cores_use %% n)
    inner_alloc <- rep(base, n)
    if (rem > 0L) inner_alloc[seq_len(rem)] <- inner_alloc[seq_len(rem)] + 1L
    inner_alloc <- pmax(1L, inner_alloc)
  }

  mode_tag <- if (n_cores_use < n) "scenario_only" else "dual_level"
  alloc_tab <- table(inner_alloc)
  alloc_tag <- paste0(names(alloc_tab), "x", as.integer(alloc_tab), collapse = ",")

  message(
    "[local] fitting ", n, " scenarios with workers=", workers,
    ", mode=", mode_tag,
    ", total_cores=", n_cores_use,
    ", inner_alloc=[", alloc_tag, "]",
    if (isTRUE(use_deoptim_local)) {
      paste0(", solver=DEoptim, itermax=", deoptim_itermax_local, ", NP=", deoptim_np_local)
    } else {
      paste0(", solver=optim_multistart, n_starts=", n_starts_local, ", maxit=", maxit_local)
    }
  )
  if (!isTRUE(use_deoptim_local) && n_cores_use >= n) {
    message("[local] NOTE: use_deoptim_local=FALSE, inner per-scenario cores are not used by optim_multistart.")
  }

  idx <- seq_len(n)
  sid <- vapply(scenarios, scenario_id, character(1))
  run_one <- function(i) {
    sc <- scenarios[[i]]
    cores_i <- as.integer(inner_alloc[[i]])
    sid_i <- sid[[i]]
    seed_i <- as.integer(seed_base + i)
    sid_hash <- if (!is.null(sid_i) && nzchar(sid_i)) sum(utf8ToInt(sid_i)) else 0L
    pid_use <- abs(as.integer(Sys.getpid()))
    cluster_port_i <- if (cores_i > 1L) {
      raw <- (as.numeric(seed_i) * 97) + (as.numeric(i) * 193) + (as.numeric(sid_hash) * 17) + (as.numeric(pid_use) * 131)
      as.integer(11000L + (raw %% 50000L))
    } else {
      NA_integer_
    }
    tryCatch(
      {
        fit <- fit_local_single(
          sc = sc,
          global_t = global_t,
          local_start = local_mat_start[i, ],
          mu_local = mu_local,
          tau_local = tau_local,
          cfg_step = cfg_step,
          full_names = full_names,
          global_names = global_names,
          local_names = local_names,
          local_lower = local_lower,
          local_upper = local_upper,
          lambda_shrink = lambda_shrink,
          tau_floor = tau_floor,
          use_deoptim = use_deoptim_local,
          deoptim_itermax = deoptim_itermax_local,
          deoptim_np = deoptim_np_local,
          deoptim_trace = deoptim_trace,
          inner_cores = cores_i,
          script_path_for_workers = script_path_for_workers,
          cluster_port = cluster_port_i,
          n_starts = n_starts_local,
          maxit = maxit_local,
          seed = seed_i
        )
        list(ok = TRUE, sid = sid_i, inner_cores = cores_i, cluster_port = cluster_port_i, fit = fit, error = "")
      },
      error = function(e) {
        list(ok = FALSE, sid = sid_i, inner_cores = cores_i, cluster_port = cluster_port_i, fit = NULL, error = conditionMessage(e))
      }
    )
  }

  out <- NULL
  if (workers > 1L && .Platform$OS.type != "windows") {
    out <- parallel::mclapply(idx, run_one, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    if (workers > 1L && .Platform$OS.type == "windows") {
      message("[local] windows detected; using serial fallback (mclapply unavailable).")
    }
    out <- lapply(idx, run_one)
  }
  out <- lapply(seq_along(out), function(i) {
    x <- out[[i]]
    if (inherits(x, "try-error")) {
      list(
        ok = FALSE,
        sid = sid[[i]],
        inner_cores = as.integer(inner_alloc[[i]]),
        cluster_port = NA_integer_,
        fit = NULL,
        error = as.character(x)
      )
    } else {
      x
    }
  })

  failed <- which(!vapply(out, function(x) isTRUE(x$ok), logical(1)))
  if (length(failed) > 0L) {
    err_lines <- vapply(
      failed,
      function(i) {
        paste0(
          out[[i]]$sid,
          " (inner_cores=", as.integer(out[[i]]$inner_cores),
          ", cluster_port=", as.integer(out[[i]]$cluster_port %||% NA_integer_),
          "): ",
          out[[i]]$error
        )
      },
      character(1)
    )
    stop(
      "[local] E-step barrier failed: ", length(failed), "/", length(out),
      " scenarios failed before global step.\n",
      paste(err_lines, collapse = "\n")
    )
  }
  message("[local] E-step barrier passed: ", length(out), "/", length(out), " scenarios completed.")

  local_mat <- matrix(NA_real_, nrow = n, ncol = length(local_names))
  rownames(local_mat) <- sid
  colnames(local_mat) <- local_names

  loss_df <- vector("list", n)
  for (i in idx) {
    fit_i <- out[[i]]$fit
    local_mat[i, ] <- as.numeric(fit_i$local_t[local_names])
    sc <- scenarios[[i]]
    loss_df[[i]] <- data.frame(
      scenario_id = scenario_id(sc),
      harvest = sc$harvest,
      cohort = sc$cohort,
      dose = sc$dose,
      objective_total = fit_i$objective_total,
      objective_data = fit_i$objective_data,
      objective_penalty = fit_i$objective_penalty,
      objective_burden = fit_i$L_b,
      objective_ploidy = fit_i$L_p,
      objective_burden_scaled = fit_i$L_b_scaled,
      objective_ploidy_scaled = fit_i$L_p_scaled,
      convergence = fit_i$convergence,
      solver = as.character(fit_i$solver %||% ""),
      inner_cores = as.integer(out[[i]]$inner_cores),
      cluster_port = as.integer(out[[i]]$cluster_port %||% NA_integer_),
      n_obs_burden = length(sc$obs_days),
      n_obs_ploidy = length(sc$ploidy_obs_N),
      data_weight = scenario_data_weight(sc),
      row.names = NULL
    )
  }

  list(local_mat = local_mat, per_sample_loss = bind_rows(loss_df))
}

evaluate_global_with_locals <- function(global_t,
                                        local_mat,
                                        scenarios,
                                        cfg_step,
                                        full_names,
                                        global_names,
                                        local_names,
                                        mu_local = NULL,
                                        tau_local = NULL,
                                        lambda_shrink = 0,
                                        tau_floor = 0.05) {
  n <- length(scenarios)
  rows <- vector("list", n)
  data_losses <- numeric(n)
  data_b <- numeric(n)
  data_p <- numeric(n)
  data_bs <- numeric(n)
  data_ps <- numeric(n)
  pen <- numeric(n)

  for (i in seq_len(n)) {
    sc <- scenarios[[i]]
    local_i <- local_mat[i, ]
    full_t <- compose_full_t(global_t, local_i, full_names, global_names, local_names)
    comp <- .eval_scenario_components(full_t, sc, cfg_step)

    if (!is.null(mu_local) && !is.null(tau_local) && lambda_shrink > 0) {
      tau_use <- pmax(as.numeric(tau_local), tau_floor)
      z <- (as.numeric(local_i[local_names]) - as.numeric(mu_local[local_names])) / tau_use
      pen[[i]] <- lambda_shrink * sum(z * z)
    } else {
      pen[[i]] <- 0
    }

    data_losses[[i]] <- comp$L_data
    data_b[[i]] <- comp$L_b
    data_p[[i]] <- comp$L_p
    data_bs[[i]] <- comp$L_b_scaled
    data_ps[[i]] <- comp$L_p_scaled

    rows[[i]] <- data.frame(
      scenario_id = scenario_id(sc),
      harvest = sc$harvest,
      cohort = sc$cohort,
      dose = sc$dose,
      objective_data = comp$L_data,
      objective_penalty = pen[[i]],
      objective_total = comp$L_data + pen[[i]],
      objective_burden = comp$L_b,
      objective_ploidy = comp$L_p,
      objective_burden_scaled = comp$L_b_scaled,
      objective_ploidy_scaled = comp$L_p_scaled,
      n_obs_burden = length(sc$obs_days),
      n_obs_ploidy = length(sc$ploidy_obs_N),
      data_weight = scenario_data_weight(sc),
      row.names = NULL
    )
  }

  per_sample <- bind_rows(rows)
  list(
    objective_data = mean(data_losses),
    objective_penalized = mean(data_losses + pen),
    objective_burden = mean(data_b),
    objective_ploidy = mean(data_p),
    objective_burden_scaled = mean(data_bs),
    objective_ploidy_scaled = mean(data_ps),
    per_sample = per_sample
  )
}

optimize_global_shared <- function(global_start,
                                   local_mat,
                                   scenarios,
                                   cfg_step,
                                   bounds,
                                   full_names,
                                   global_names,
                                   local_names,
                                   use_deoptim = TRUE,
                                   deoptim_itermax = 180L,
                                   deoptim_np = 160L,
                                   deoptim_trace = TRUE,
                                   deoptim_parallel = TRUE,
                                   n_cores = 1L,
                                   script_path_for_workers = NULL,
                                   n_starts = 12,
                                   maxit = 6000,
                                   seed = 1) {
  set.seed(seed)

  lower <- bounds$lower[global_names]
  upper <- bounds$upper[global_names]
  mid <- (lower + upper) / 2
  start0 <- clip(global_start[global_names], lower, upper)

  objective_global <- function(g_sub) {
    g <- start0
    g[names(g_sub)] <- g_sub
    eval <- evaluate_global_with_locals(
      global_t = g,
      local_mat = local_mat,
      scenarios = scenarios,
      cfg_step = cfg_step,
      full_names = full_names,
      global_names = global_names,
      local_names = local_names,
      mu_local = NULL,
      tau_local = NULL,
      lambda_shrink = 0
    )
    eval$objective_data
  }

  if (isTRUE(use_deoptim)) {
    message(
      "[global] Starting DEoptim with itermax=", as.integer(max(1L, deoptim_itermax)),
      ", NP=", as.integer(max(4L, deoptim_np)),
      ", n_cores=", as.integer(max(1L, n_cores)),
      ", parallel=", if (isTRUE(deoptim_parallel)) "TRUE" else "FALSE"
    )
    de_fit <- .deoptim_optimize_box(
      objective_fn = objective_global,
      lower = lower,
      upper = upper,
      itermax = deoptim_itermax,
      NP = deoptim_np,
      init_vec = start0,
      alt_init = mid,
      trace = deoptim_trace,
      seed = seed,
      label = "global",
      deoptim_parallel = deoptim_parallel,
      n_cores = n_cores,
      script_path_for_workers = script_path_for_workers,
      strict_parallel = isTRUE(deoptim_parallel) && as.integer(max(1L, n_cores)) > 1L
    )
    best_global <- de_fit$par
    best_val <- de_fit$value
    best_conv <- de_fit$convergence
    best_method <- de_fit$method
    best_parallel_info <- de_fit$parallel_info %||% NULL
  } else {
    starts <- vector("list", n_starts)
    starts[[1]] <- start0
    if (n_starts >= 2) starts[[2]] <- mid
    if (n_starts >= 3) {
      for (k in 3:n_starts) {
        starts[[k]] <- stats::runif(length(lower), min = lower, max = upper)
        names(starts[[k]]) <- names(lower)
      }
    }

    best <- list(par = start0, value = Inf, conv = 999)
    for (k in seq_along(starts)) {
      fit <- tryCatch(
        optim(
          par = starts[[k]],
          fn = objective_global,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          control = list(maxit = as.integer(maxit))
        ),
        error = function(e) NULL
      )
      if (!is.null(fit) && is.finite(fit$value) && fit$value < best$value) {
        best <- list(par = fit$par, value = fit$value, conv = fit$convergence)
      }
    }
    best_global <- clip(best$par, lower, upper)
    names(best_global) <- names(lower)
    best_val <- best$value
    best_conv <- best$conv
    best_method <- "optim_multistart"
    best_parallel_info <- NULL
  }

  list(
    best_global = best_global,
    objective = best_val,
    convergence = best_conv,
    solver = best_method,
    parallel_info = best_parallel_info
  )
}

collect_predictions_hierarchical <- function(global_t,
                                             local_mat,
                                             scenarios,
                                             cfg,
                                             full_names,
                                             global_names,
                                             local_names) {
  pred_one <- function(sc, i) {
    full_t <- compose_full_t(global_t, local_mat[i, ], full_names, global_names, local_names)
    rp <- decode_params(
      full_t,
      fit_full_pmis = cfg$fit_full_pmis,
      fit_treatment = cfg$fit_treatment
    )
    sim <- simulate_one(rp, sc, cfg, model_core = NULL)

    obs <- as.numeric(sc$obs_burden)
    pred_pop <- as.numeric(sim$Ntot_obs)
    pred_vol <- as.numeric(sim$Vmm3_obs)
    obs_delta <- obs - obs[1]
    pred_delta <- pred_vol - pred_vol[1]
    s_obs <- max(abs(obs_delta), na.rm = TRUE)
    s_pred <- max(abs(pred_delta), na.rm = TRUE)
    log_eps <- pmax(as.numeric(.first_non_null_local(cfg$burden_log_eps, 1e-12)), 1e-15)

    burden_df <- data.frame(
      scenario_id = scenario_id(sc),
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
      obs_norm = if (s_obs > 0) obs_delta / s_obs else obs_delta,
      pred_norm = if (s_pred > 0) pred_delta / s_pred else pred_delta,
      row.names = NULL
    )

    obs_tab <- table(sc$ploidy_obs_N)
    ploidy_df <- data.frame(
      scenario_id = scenario_id(sc),
      harvest = sc$harvest,
      cohort = sc$cohort,
      dose = sc$dose,
      N = as.integer(names(sim$frac_N)),
      pred_fraction = as.numeric(sim$frac_N),
      row.names = NULL
    )
    ploidy_df$obs_count <- as.numeric(obs_tab[as.character(ploidy_df$N)])
    ploidy_df$obs_count[is.na(ploidy_df$obs_count)] <- 0

    list(burden = burden_df, ploidy = ploidy_df)
  }

  pred_rows <- map_scenarios_parallel(
    scenarios = scenarios,
    n_cores = cfg$n_cores,
    label = "predict_hier",
    fn = pred_one
  )

  list(
    burden = bind_rows(lapply(pred_rows, function(x) x$burden)),
    ploidy = bind_rows(lapply(pred_rows, function(x) x$ploidy))
  )
}

build_sample_theta_tables <- function(global_t,
                                      local_mat,
                                      scenarios,
                                      cfg,
                                      full_names,
                                      global_names,
                                      local_names) {
  rows_t <- list()
  rows_nat <- list()

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    sid <- scenario_id(sc)
    full_t <- compose_full_t(global_t, local_mat[i, ], full_names, global_names, local_names)
    rp <- decode_params(full_t, fit_full_pmis = cfg$fit_full_pmis, fit_treatment = cfg$fit_treatment)

    base_meta <- list(
      scenario_id = sid,
      harvest = sc$harvest,
      cohort = sc$cohort,
      dose = sc$dose,
      n_obs_burden = length(sc$obs_days),
      n_obs_ploidy = length(sc$ploidy_obs_N),
      data_weight = scenario_data_weight(sc)
    )

    row_t <- as.list(full_t)
    row_t <- c(base_meta, row_t)
    rows_t[[i]] <- as.data.frame(row_t, check.names = FALSE)

    row_nat <- c(base_meta, rp)
    rows_nat[[i]] <- as.data.frame(row_nat, check.names = FALSE)
  }

  list(
    transformed = bind_rows(rows_t),
    natural = bind_rows(rows_nat)
  )
}

select_best_step <- function(chain_df, rule = "min_objective_data") {
  if (nrow(chain_df) == 0) stop("No rows in chain summary.")
  rule <- as.character(rule %||% "min_objective_data")

  if (rule == "min_objective_data") {
    ord <- order(chain_df$objective_data, chain_df$objective_penalized)
    idx <- ord[[1]]
    reason <- "Minimize objective_data (scaled weighted burden/ploidy mean over scenarios)."
  } else if (rule == "balanced_rank") {
    rb <- rank(chain_df$objective_burden_scaled, ties.method = "average")
    rp <- rank(chain_df$objective_ploidy_scaled, ties.method = "average")
    score <- rb + rp
    ord <- order(score, chain_df$objective_data)
    idx <- ord[[1]]
    reason <- "Minimize rank(objective_burden_scaled)+rank(objective_ploidy_scaled), tie-break by objective_data."
  } else {
    stop("Unknown select_rule: ", rule)
  }

  list(index = idx, rule = rule, reason = reason)
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  default_data_dir <- normalizePath(file.path(script_dir, "..", "..", "..", "data", "InVivoData_Gemcitabine"), mustWork = FALSE)
  data_dir <- if (!is.null(argv$data_dir)) argv$data_dir else default_data_dir

  truncate_at_treatment <- if (!is.null(argv$truncate_at_treatment)) {
    as_bool(argv$truncate_at_treatment, FALSE)
  } else {
    as_bool(argv$pretreat_only, FALSE)
  }

  w_burden_vec <- as_num_vec(argv$w_burden, 1.0)
  w_ploidy_vec <- as_num_vec(argv$w_ploidy, 1.0)
  weight_schedule <- make_weight_schedule(w_burden_vec, w_ploidy_vec)

  n_cores_arg <- as_int(argv$n_cores, NA_integer_)
  n_cores_use <- if (is.finite(n_cores_arg)) n_cores_arg else default_n_cores()

  cfg <- list(
    model_path = model_script,
    N_UNIT = as_int(argv$N_UNIT, 22L),
    N_MIN = as_int(argv$N_MIN, 22L),
    N_MAX = as_int(argv$N_MAX, 154L),
    DT = as_num(argv$dt, 0.5),
    O2_fixed = as_num(argv$O2, 1.0),
    o2_burden_feedback = if (!is.null(argv$o2_burden_feedback)) {
      as_bool(argv$o2_burden_feedback, TRUE)
    } else {
      as_bool(argv$o2_dynamic, TRUE)
    },
    o2_min = as_num(argv$o2_min, 0.0),
    h_O2 = as_num(argv$h_O2, 1.0),
    K = as_num(argv$K, 1e12),
    crowding = if (!is.null(argv$crowding)) argv$crowding else "logistic",
    init_total_size = as_num(argv$init_total_size, 1e6),
    dose_ref = as_num(argv$dose_ref, 30),
    tx_mult_min = as_num(argv$tx_mult_min, 0.05),
    min_pop = as_num(argv$min_pop, 1e-12),
    huber_k = as_num(argv$huber_k, 0.1),
    huber_k_burden_log = as_num(argv$huber_k_burden_log, as_num(argv$huber_k, 0.1)),
    burden_log_eps = as_num(argv$burden_log_eps, 1e-12),
    rho_2N_min = as_num(argv$rho_2N_min, 3.2e4),
    rho_2N_max = as_num(argv$rho_2N_max, 5.6e4),
    legacy_c_vol_2N_mm3 = as_num(argv$c_vol_2N_mm3, 4.19e-06),
    w_burden = weight_schedule$w_burden[[1]],
    w_ploidy = weight_schedule$w_ploidy[[1]],
    w_burden_schedule = weight_schedule$w_burden,
    w_ploidy_schedule = weight_schedule$w_ploidy,
    n_weight_passes = weight_schedule$n_pass,
    loss_rescale = as_bool(argv$loss_rescale, TRUE),
    loss_scale_burden = as_num(argv$loss_scale_burden, NA_real_),
    loss_scale_ploidy = as_num(argv$loss_scale_ploidy, NA_real_),
    loss_scale_eps = as_num(argv$loss_scale_eps, 1e-8),
    loss_scale_source = "unset",
    optim_trace = as_bool(argv$optim_trace, TRUE),
    optim_trace_every = as_int(argv$optim_trace_every, 1L),
    eps_prob = as_num(argv$eps_prob, 1e-12),
    trace_obj = as_bool(argv$trace_obj, FALSE),
    fit_full_pmis = as_bool(argv$fit_full_pmis, FALSE),
    fit_treatment = as_bool(argv$fit_treatment, FALSE),
    dose_zero_only = as_bool(argv$dose_zero_only, TRUE),
    paired_only = as_bool(argv$paired_only, TRUE),
    truncate_at_treatment = truncate_at_treatment,
    ploidy_at_harvest = as_bool(argv$ploidy_at_harvest, TRUE),
    two_stage = FALSE,
    stage1_w_burden = 1.0,
    stage1_w_ploidy = 0.0,
    stage2_w_burden = 0.0,
    stage2_w_ploidy = 1.0,
    use_deoptim = FALSE,
    deoptim_parallel = FALSE,
    itermax = as_int(argv$itermax, 80L),
    NP = as_int(argv$NP, 120L),
    n_cores = n_cores_use,
    seed = as_int(argv$seed, 1L),
    max_scenarios = as_num(argv$max_scenarios, Inf)
  )

  n_alt_iter <- as_int(argv$n_alt_iter, 3L)
  n_starts_local <- as_int(argv$n_starts_local, 6L)
  n_starts_global <- as_int(argv$n_starts_global, 12L)
  maxit_local <- as_int(argv$maxit_local, 2500L)
  maxit_global <- as_int(argv$maxit_global, 6000L)
  use_deoptim_local <- as_bool(argv$use_deoptim_local, TRUE)
  use_deoptim_global <- as_bool(argv$use_deoptim_global, TRUE)
  deoptim_itermax_local <- as_int(argv$deoptim_itermax_local, 120L)
  deoptim_np_local <- as_int(argv$deoptim_np_local, 80L)
  deoptim_itermax_global <- as_int(argv$deoptim_itermax_global, 180L)
  deoptim_np_global <- as_int(argv$deoptim_np_global, 160L)
  deoptim_trace <- as_bool(argv$deoptim_trace, TRUE)
  lambda_shrink <- as_num(argv$lambda_shrink, 1.0)
  tau_floor <- as_num(argv$tau_floor, 0.05)
  tau_min <- as_num(argv$tmb_tau_min, 1e-3)
  tmb_log_tau_prior_sd <- as_num(argv$tmb_log_tau_prior_sd, 2.0)
  tmb_maxit <- as_int(argv$tmb_maxit, 200L)
  select_rule <- as.character(argv$select_rule %||% "min_objective_data")
  tmb_rebuild <- as_bool(argv$tmb_rebuild, FALSE)

  if (cfg$DT <= 0) stop("dt must be > 0")
  if (cfg$N_MAX < cfg$N_MIN) stop("N_MAX must be >= N_MIN")
  if (!is.finite(cfg$o2_min)) stop("o2_min must be finite")
  cfg$o2_min <- clip(cfg$o2_min, 0, 1)
  if (!is.finite(cfg$h_O2) || cfg$h_O2 <= 0) stop("h_O2 must be > 0")
  if (!is.finite(cfg$burden_log_eps) || cfg$burden_log_eps <= 0) stop("burden_log_eps must be > 0")
  if (!is.finite(cfg$rho_2N_min) || cfg$rho_2N_min <= 0) stop("rho_2N_min must be > 0")
  if (!is.finite(cfg$rho_2N_max) || cfg$rho_2N_max <= 0) stop("rho_2N_max must be > 0")
  if (cfg$rho_2N_max < cfg$rho_2N_min) stop("rho_2N_max must be >= rho_2N_min")
  if (!is.finite(cfg$legacy_c_vol_2N_mm3) || cfg$legacy_c_vol_2N_mm3 <= 0) stop("c_vol_2N_mm3 (legacy warm-start baseline) must be > 0")
  if (cfg$n_cores < 1) stop("n_cores must be >= 1")
  if (n_alt_iter < 1) stop("n_alt_iter must be >= 1")
  if (n_starts_local < 1 || n_starts_global < 1) stop("n_starts_local/n_starts_global must be >= 1")
  if (maxit_local < 1 || maxit_global < 1) stop("maxit_local/maxit_global must be >= 1")
  if (deoptim_itermax_local < 1 || deoptim_itermax_global < 1) stop("deoptim_itermax_local/deoptim_itermax_global must be >= 1")
  if (deoptim_np_local < 4 || deoptim_np_global < 4) stop("deoptim_np_local/deoptim_np_global must be >= 4")
  if (!is.finite(lambda_shrink) || lambda_shrink < 0) stop("lambda_shrink must be >= 0")
  if (!is.finite(tau_floor) || tau_floor <= 0) stop("tau_floor must be > 0")
  if ((use_deoptim_local || use_deoptim_global) && !requireNamespace("DEoptim", quietly = TRUE)) {
    stop("DEoptim requested (local/global) but package DEoptim is not available.")
  }

  set.seed(cfg$seed)

  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- file.path(data_dir, "all_ploidy.tsv")
  scenarios <- prepare_data(dt_path, ploidy_path, cfg)

  bounds <- make_bounds(
    fit_full_pmis = cfg$fit_full_pmis,
    fit_treatment = cfg$fit_treatment,
    rho_2N_min = cfg$rho_2N_min,
    rho_2N_max = cfg$rho_2N_max
  )
  full_names <- names(bounds$lower)
  default_par_t <- (bounds$lower + bounds$upper) / 2
  names(default_par_t) <- full_names
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
    "Burden observation model enabled: log-volume Huber on V(mm^3), ",
    "V_pred = sum_n n_n * [(1/rho_2N) * (P/2)^beta_size], ",
    "rho_2N_range=[", signif(cfg$rho_2N_min, 6), ", ", signif(cfg$rho_2N_max, 6), "] cells/mm^3"
  )
  message(
    "O2 mode: ",
    if (isTRUE(cfg$o2_burden_feedback)) "dynamic feedback" else "fixed",
    "; O2_base=", signif(cfg$O2_fixed, 6),
    ", o2_min=", signif(cfg$o2_min, 6),
    ", h_O2=", signif(cfg$h_O2, 6),
    if (isTRUE(cfg$o2_burden_feedback)) "; K_O2 is fitted." else "; K_O2 is inactive."
  )

  default_local <- intersect(full_names, c("log10_lam_min", "log10_lam_max", "log10_p_misseg", "log10_p_wgd"))
  local_from_cli <- as_char_vec(argv$local_params, character(0))
  local_names <- if (length(local_from_cli) > 0) local_from_cli else default_local
  if (!all(local_names %in% full_names)) {
    stop("local_params contains unknown transformed parameters: ", paste(setdiff(local_names, full_names), collapse = ", "))
  }
  if (length(local_names) == 0) stop("No local parameters selected.")

  global_names <- setdiff(full_names, local_names)
  if (length(global_names) == 0) stop("No global parameters left after local/global split.")

  message("Local transformed parameters: ", paste(local_names, collapse = ", "))
  message("Global transformed parameters: ", paste(global_names, collapse = ", "))
  message(
    "Solvers: local=", if (use_deoptim_local) "DEoptim" else "optim_multistart",
    ", global=", if (use_deoptim_global) "DEoptim" else "optim_multistart",
    ", deoptim_trace=", deoptim_trace,
    ", global_deoptim_parallel=TRUE",
    ", n_cores=", cfg$n_cores
  )

  sample_ids <- make.unique(vapply(scenarios, scenario_id, character(1)))
  sample_weight <- vapply(scenarios, scenario_data_weight, numeric(1))
  self_script_path <- normalizePath(
    file.path(script_dir, "fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.R"),
    mustWork = FALSE
  )

  base_t <- if (!is.null(warm_start_t)) warm_start_t else default_par_t
  global_t_cur <- base_t[global_names]
  local_mat_cur <- matrix(rep(base_t[local_names], each = length(scenarios)), nrow = length(scenarios), byrow = FALSE)
  rownames(local_mat_cur) <- sample_ids
  colnames(local_mat_cur) <- local_names

  robust_init <- robust_pool_local(local_mat_cur, sample_weight, tau_floor = tau_floor)
  mu_local_cur <- robust_init$mu
  tau_local_cur <- robust_init$tau

  dll_name <- .init_tmb_hier_backend(script_dir, rebuild = tmb_rebuild)
  message("Using TMB DLL: ", dll_name)

  ts_format <- as.character(argv$timestamp_format %||% "%Y%m%d_%H%M%S")
  run_stamp <- format(Sys.time(), ts_format)
  out_dir <- if (!is.null(argv$out_dir)) {
    if (as_bool(argv$append_timestamp_out_dir, FALSE)) paste0(argv$out_dir, "_", run_stamp) else argv$out_dir
  } else {
    file.path(script_dir, "..", "..", "results", paste0("fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain_", run_stamp))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(cfg, file = file.path(out_dir, "fit_config.rds"))

  chain_rows <- list()

  for (pass_i in seq_len(cfg$n_weight_passes)) {
    wb <- cfg$w_burden_schedule[[pass_i]]
    wp <- cfg$w_ploidy_schedule[[pass_i]]
    pass_cfg <- cfg
    pass_cfg$w_burden <- wb
    pass_cfg$w_ploidy <- wp

    step_dir <- file.path(
      out_dir,
      sprintf("step%02d_wb%s_wp%s", pass_i, sanitize_tag(wb), sanitize_tag(wp))
    )
    dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)

    message("[step", pass_i, "] weights: w_burden=", wb, ", w_ploidy=", wp)

    alt_rows <- list()
    for (alt_i in seq_len(n_alt_iter)) {
      message("[step", pass_i, "][iter", alt_i, "] local E-step")
      local_fit <- fit_locals_parallel(
        scenarios = scenarios,
        global_t = global_t_cur,
        local_mat_start = local_mat_cur,
        mu_local = mu_local_cur,
        tau_local = tau_local_cur,
        cfg_step = pass_cfg,
        full_names = full_names,
        global_names = global_names,
        local_names = local_names,
        local_lower = bounds$lower[local_names],
        local_upper = bounds$upper[local_names],
        lambda_shrink = lambda_shrink,
        tau_floor = tau_floor,
        n_starts_local = n_starts_local,
        maxit_local = maxit_local,
        use_deoptim_local = use_deoptim_local,
        deoptim_itermax_local = deoptim_itermax_local,
        deoptim_np_local = deoptim_np_local,
        deoptim_trace = deoptim_trace,
        n_cores = cfg$n_cores,
        script_path_for_workers = self_script_path,
        seed_base = cfg$seed + pass_i * 10000L + alt_i * 1000L
      )

      local_hat_mat <- local_fit$local_mat
      robust_pool <- robust_pool_local(local_hat_mat, sample_weight, tau_floor = tau_floor)

      message("[step", pass_i, "][iter", alt_i, "] TMB pooling")
      tmb_pool <- pool_locals_with_tmb(
        local_hat_mat = local_hat_mat,
        sample_weight = sample_weight,
        dll_name = dll_name,
        init_mu = robust_pool$mu,
        init_tau = robust_pool$tau,
        tau_min = tau_min,
        log_tau_prior_sd = tmb_log_tau_prior_sd,
        maxit = tmb_maxit
      )

      local_mat_cur <- tmb_pool$theta_shrunk
      mu_local_cur <- tmb_pool$mu
      tau_local_cur <- pmax(tmb_pool$tau, tau_floor)

      message("[step", pass_i, "][iter", alt_i, "] global M-step")
      global_fit <- optimize_global_shared(
        global_start = global_t_cur,
        local_mat = local_mat_cur,
        scenarios = scenarios,
        cfg_step = pass_cfg,
        bounds = bounds,
        full_names = full_names,
        global_names = global_names,
        local_names = local_names,
        use_deoptim = use_deoptim_global,
        deoptim_itermax = deoptim_itermax_global,
        deoptim_np = deoptim_np_global,
        deoptim_trace = deoptim_trace,
        deoptim_parallel = TRUE,
        n_cores = cfg$n_cores,
        script_path_for_workers = self_script_path,
        n_starts = n_starts_global,
        maxit = maxit_global,
        seed = cfg$seed + pass_i * 10000L + alt_i * 2000L + 123L
      )

      global_t_cur <- global_fit$best_global

      eval_now <- evaluate_global_with_locals(
        global_t = global_t_cur,
        local_mat = local_mat_cur,
        scenarios = scenarios,
        cfg_step = pass_cfg,
        full_names = full_names,
        global_names = global_names,
        local_names = local_names,
        mu_local = mu_local_cur,
        tau_local = tau_local_cur,
        lambda_shrink = lambda_shrink,
        tau_floor = tau_floor
      )

      g_par <- global_fit$parallel_info %||% list(
        requested_workers = as.integer(max(1L, cfg$n_cores)),
        started_workers = if (cfg$n_cores > 1L) NA_integer_ else 1L,
        active_workers = 1L,
        cluster_port = NA_integer_,
        fallback = FALSE,
        fallback_reason = "",
        objective_fail_count = NA_integer_
      )

      alt_rows[[alt_i]] <- data.frame(
        pass = pass_i,
        alt_iter = alt_i,
        w_burden = wb,
        w_ploidy = wp,
        objective_data = eval_now$objective_data,
        objective_penalized = eval_now$objective_penalized,
        objective_burden = eval_now$objective_burden,
        objective_ploidy = eval_now$objective_ploidy,
        objective_burden_scaled = eval_now$objective_burden_scaled,
        objective_ploidy_scaled = eval_now$objective_ploidy_scaled,
        tmb_objective = tmb_pool$objective,
        tmb_convergence = tmb_pool$convergence,
        local_solver = if (use_deoptim_local) "DEoptim" else "optim_multistart",
        global_solver = as.character(global_fit$solver %||% if (use_deoptim_global) "DEoptim" else "optim_multistart"),
        global_convergence = global_fit$convergence,
        global_requested_workers = as.integer(g_par$requested_workers %||% NA_integer_),
        global_started_workers = as.integer(g_par$started_workers %||% NA_integer_),
        global_active_workers = as.integer(g_par$active_workers %||% NA_integer_),
        global_cluster_port = as.integer(g_par$cluster_port %||% NA_integer_),
        global_fallback = as.logical(g_par$fallback %||% FALSE),
        global_objective_fail_count = as.integer(g_par$objective_fail_count %||% NA_integer_),
        global_fallback_reason = as.character(g_par$fallback_reason %||% ""),
        row.names = NULL
      )
    }

    alt_df <- bind_rows(alt_rows)
    write.table(alt_df, file = file.path(step_dir, "alt_iter_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    eval_final <- evaluate_global_with_locals(
      global_t = global_t_cur,
      local_mat = local_mat_cur,
      scenarios = scenarios,
      cfg_step = pass_cfg,
      full_names = full_names,
      global_names = global_names,
      local_names = local_names,
      mu_local = mu_local_cur,
      tau_local = tau_local_cur,
      lambda_shrink = lambda_shrink,
      tau_floor = tau_floor
    )

    theta_tables <- build_sample_theta_tables(
      global_t = global_t_cur,
      local_mat = local_mat_cur,
      scenarios = scenarios,
      cfg = cfg,
      full_names = full_names,
      global_names = global_names,
      local_names = local_names
    )
    write.table(theta_tables$natural, file = file.path(step_dir, "per_sample_theta_i.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(theta_tables$transformed, file = file.path(step_dir, "per_sample_theta_i_transformed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(eval_final$per_sample, file = file.path(step_dir, "per_sample_loss.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    robust_final <- robust_pool_local(local_mat_cur, sample_weight, tau_floor = tau_floor)
    theta0_robust_t <- compose_full_t(global_t_cur, robust_final$mu, full_names, global_names, local_names)
    theta0_robust <- decode_params(theta0_robust_t, fit_full_pmis = cfg$fit_full_pmis, fit_treatment = cfg$fit_treatment)
    theta0_df <- bind_rows(
      data.frame(
        space = "natural",
        parameter = names(theta0_robust),
        value = as.numeric(theta0_robust),
        row.names = NULL
      ),
      data.frame(
        space = "transformed",
        parameter = full_names,
        value = as.numeric(theta0_robust_t[full_names]),
        row.names = NULL
      )
    )
    write.table(theta0_df, file = file.path(step_dir, "theta0_robust.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    shared_global_t <- numeric(length(full_names))
    names(shared_global_t) <- full_names
    shared_global_t[global_names] <- global_t_cur
    shared_global_t[local_names] <- robust_final$mu
    shared_global <- decode_params(shared_global_t, fit_full_pmis = cfg$fit_full_pmis, fit_treatment = cfg$fit_treatment)
    shared_df <- data.frame(parameter = names(shared_global), value = as.numeric(shared_global), row.names = NULL)
    write.table(shared_df, file = file.path(step_dir, "global_best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    preds <- collect_predictions_hierarchical(
      global_t = global_t_cur,
      local_mat = local_mat_cur,
      scenarios = scenarios,
      cfg = pass_cfg,
      full_names = full_names,
      global_names = global_names,
      local_names = local_names
    )
    write.table(preds$burden, file = file.path(step_dir, "global_burden_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(preds$ploidy, file = file.path(step_dir, "global_terminal_ploidy_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    g_par_final <- global_fit$parallel_info %||% list(
      requested_workers = as.integer(max(1L, cfg$n_cores)),
      started_workers = if (cfg$n_cores > 1L) NA_integer_ else 1L,
      active_workers = 1L,
      cluster_port = NA_integer_,
      fallback = FALSE,
      fallback_reason = "",
      objective_fail_count = NA_integer_
    )

    fit_summary <- data.frame(
      metric = c(
        "pass", "w_burden", "w_ploidy", "n_alt_iter", "objective_data", "objective_penalized",
        "objective_burden", "objective_ploidy", "objective_burden_scaled", "objective_ploidy_scaled",
        "n_scenarios", "n_cores", "lambda_shrink", "tau_floor", "tmb_tau_min", "select_rule",
        "local_solver", "global_solver", "deoptim_itermax_local", "deoptim_np_local",
        "deoptim_itermax_global", "deoptim_np_global", "deoptim_trace",
        "global_deoptim_parallel", "global_requested_workers", "global_started_workers",
        "global_active_workers", "global_cluster_port", "global_fallback",
        "global_objective_fail_count", "global_fallback_reason"
      ),
      value = c(
        as.character(pass_i), as.character(wb), as.character(wp), as.character(n_alt_iter),
        as.character(eval_final$objective_data), as.character(eval_final$objective_penalized),
        as.character(eval_final$objective_burden), as.character(eval_final$objective_ploidy),
        as.character(eval_final$objective_burden_scaled), as.character(eval_final$objective_ploidy_scaled),
        as.character(length(scenarios)), as.character(cfg$n_cores), as.character(lambda_shrink),
        as.character(tau_floor), as.character(tau_min), as.character(select_rule),
        as.character(if (use_deoptim_local) "DEoptim" else "optim_multistart"),
        as.character(if (use_deoptim_global) "DEoptim" else "optim_multistart"),
        as.character(deoptim_itermax_local), as.character(deoptim_np_local),
        as.character(deoptim_itermax_global), as.character(deoptim_np_global),
        as.character(deoptim_trace),
        as.character(TRUE),
        as.character(g_par_final$requested_workers %||% NA_integer_),
        as.character(g_par_final$started_workers %||% NA_integer_),
        as.character(g_par_final$active_workers %||% NA_integer_),
        as.character(g_par_final$cluster_port %||% NA_integer_),
        as.character(g_par_final$fallback %||% FALSE),
        as.character(g_par_final$objective_fail_count %||% NA_integer_),
        as.character(g_par_final$fallback_reason %||% "")
      ),
      row.names = NULL
    )
    write.table(fit_summary, file = file.path(step_dir, "global_fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    chain_rows[[pass_i]] <- data.frame(
      pass = pass_i,
      w_burden = wb,
      w_ploidy = wp,
      objective_data = eval_final$objective_data,
      objective_penalized = eval_final$objective_penalized,
      objective_burden = eval_final$objective_burden,
      objective_ploidy = eval_final$objective_ploidy,
      objective_burden_scaled = eval_final$objective_burden_scaled,
      objective_ploidy_scaled = eval_final$objective_ploidy_scaled,
      local_solver = if (use_deoptim_local) "DEoptim" else "optim_multistart",
      global_solver = if (use_deoptim_global) "DEoptim" else "optim_multistart",
      global_requested_workers = as.integer(g_par_final$requested_workers %||% NA_integer_),
      global_started_workers = as.integer(g_par_final$started_workers %||% NA_integer_),
      global_active_workers = as.integer(g_par_final$active_workers %||% NA_integer_),
      global_cluster_port = as.integer(g_par_final$cluster_port %||% NA_integer_),
      global_fallback = as.logical(g_par_final$fallback %||% FALSE),
      global_objective_fail_count = as.integer(g_par_final$objective_fail_count %||% NA_integer_),
      global_fallback_reason = as.character(g_par_final$fallback_reason %||% ""),
      n_scenarios = length(scenarios),
      step_dir = normalizePath(step_dir),
      row.names = NULL
    )
  }

  chain_df <- bind_rows(chain_rows)
  write.table(chain_df, file = file.path(out_dir, "chain_global_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  sel <- select_best_step(chain_df, rule = select_rule)
  selected <- chain_df[sel$index, , drop = FALSE]
  selected$select_rule <- sel$rule
  selected$select_reason <- sel$reason
  write.table(selected, file = file.path(out_dir, "selected_best_step.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  ranked <- chain_df %>% arrange(objective_data, objective_penalized)
  write.table(ranked, file = file.path(out_dir, "chain_global_summary_ranked.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  run_cfg <- data.frame(
    key = c(
      "local_params", "global_params", "n_alt_iter", "n_starts_local", "n_starts_global",
      "maxit_local", "maxit_global", "lambda_shrink", "tau_floor", "tmb_tau_min",
      "tmb_log_tau_prior_sd", "tmb_maxit", "paired_only", "fit_treatment",
      "o2_burden_feedback", "o2_min", "h_O2", "O2_fixed",
      "use_deoptim_local", "use_deoptim_global", "deoptim_itermax_local", "deoptim_np_local",
      "deoptim_itermax_global", "deoptim_np_global", "deoptim_trace", "global_deoptim_parallel"
    ),
    value = c(
      paste(local_names, collapse = ","), paste(global_names, collapse = ","), as.character(n_alt_iter),
      as.character(n_starts_local), as.character(n_starts_global), as.character(maxit_local), as.character(maxit_global),
      as.character(lambda_shrink), as.character(tau_floor), as.character(tau_min),
      as.character(tmb_log_tau_prior_sd), as.character(tmb_maxit), as.character(cfg$paired_only), as.character(cfg$fit_treatment),
      as.character(cfg$o2_burden_feedback), as.character(cfg$o2_min), as.character(cfg$h_O2), as.character(cfg$O2_fixed),
      as.character(use_deoptim_local), as.character(use_deoptim_global), as.character(deoptim_itermax_local), as.character(deoptim_np_local),
      as.character(deoptim_itermax_global), as.character(deoptim_np_global), as.character(deoptim_trace), as.character(TRUE)
    ),
    row.names = NULL
  )
  write.table(run_cfg, file = file.path(out_dir, "run_config_extra.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  message("Done. Output directory: ", normalizePath(out_dir))
  message("Best step (", sel$rule, "): pass=", selected$pass[[1]],
          " w=", selected$w_burden[[1]], ",", selected$w_ploidy[[1]],
          " objective_data=", signif(selected$objective_data[[1]], 6))
}

if (sys.nframe() == 0) {
  main()
}
