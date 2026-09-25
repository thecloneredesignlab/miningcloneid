# Fixed-O2 numerical producers extracted from the legacy monolith.
# This module owns model-based attractors, analytical trajectories,
# Euler/expm/eigen validation, missing simulation task generation, and
# materialized simulation-trajectory readers.

fo2_dominant_attractor_one <- function(seed_id, run_params, model_env, cfg, O2) {
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22L), as.integer(cfg$N_MAX %||% 154L))
  G <- o2pr_build_G(model_env, cfg, run_params, O2)
  mu_all <- as.numeric(o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(ngrid)), run_params = run_params, N = ngrid))
  M <- G - Matrix::Diagonal(x = mu_all)
  eig <- tryCatch(eigen(as.matrix(M), only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) {
    return(data.frame(
      seed_id = seed_id,
      O2_pct = O2,
      status = "eigen_failed",
      dominant_mean_N = NA_real_,
      dominant_mean_ploidy = NA_real_,
      dominant_fraction_N_le_25 = NA_real_,
      dominant_fraction_N_below_44 = NA_real_,
      dominant_fraction_N_ge_44 = NA_real_,
      dominant_growth_rate = NA_real_,
      spectral_gap = NA_real_,
      eigenvector_nonnegative = NA,
      selection_22_vs_44 = NA_real_,
      selection_44_vs_88 = NA_real_,
      eff_growth_N22 = NA_real_,
      eff_growth_N44 = NA_real_,
      eff_growth_N88 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  idx <- which.max(Re(eig$values))
  v <- Re(eig$vectors[, idx])
  if (sum(v, na.rm = TRUE) < 0) v <- -v
  nonneg <- all(v >= -1e-8, na.rm = TRUE)
  v <- pmax(v, 0)
  status <- "ok"
  if (!sum(v, na.rm = TRUE) > 0) {
    v <- rep(NA_real_, length(ngrid))
    status <- "empty_dominant_vector_after_truncation"
  } else {
    v <- v / sum(v, na.rm = TRUE)
  }
  lambda1 <- Re(eig$values[[idx]])
  lambda2 <- sort(Re(eig$values), decreasing = TRUE)[min(2L, length(eig$values))]
  eff <- vapply(c(22L, 44L, 88L), function(N) {
    col <- as.integer(N - (cfg$N_MIN %||% 22L) + 1L)
    if (col < 1L || col > ncol(G)) return(NA_real_)
    sum(G[, col]) - mu_all[[col]]
  }, numeric(1))
  names(eff) <- c("22", "44", "88")
  data.frame(
    seed_id = seed_id,
    O2_pct = O2,
    status = status,
    dominant_mean_N = sum(ngrid * v, na.rm = TRUE),
    dominant_mean_ploidy = sum(ngrid * v, na.rm = TRUE) / as.numeric(cfg$N_UNIT %||% 22),
    dominant_fraction_N_le_25 = sum(v[ngrid <= 25], na.rm = TRUE),
    dominant_fraction_N_below_44 = sum(v[ngrid < 44], na.rm = TRUE),
    dominant_fraction_N_ge_44 = sum(v[ngrid >= 44], na.rm = TRUE),
    dominant_growth_rate = lambda1,
    spectral_gap = lambda1 - lambda2,
    eigenvector_nonnegative = nonneg,
    selection_22_vs_44 = eff[["22"]] - eff[["44"]],
    selection_44_vs_88 = eff[["44"]] - eff[["88"]],
    eff_growth_N22 = eff[["22"]],
    eff_growth_N44 = eff[["44"]],
    eff_growth_N88 = eff[["88"]],
    stringsAsFactors = FALSE
  )
}


cf2_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  fixo2_fixed_matrix(model_env, cfg, run_params, O2)
}


cf2_init_vector <- function(ngrid, init_N) {
  fixo2_init_vector(ngrid, init_N, n_unit = 22)
}


cf2_normalize_state <- function(x) {
  fixo2_normalize_state(x)
}


cf2_eigen_trajectory <- function(M, ngrid, init, time_grid, n_unit) {
  fixo2_eigen_trajectory(M, ngrid, init, time_grid, n_unit)
}


cf2_dominant_one <- function(M, ngrid, n_unit) {
  fixo2_dominant_one(M, ngrid, n_unit)
}


fixo2_simulation_mode_objective_metadata <- function(analysis_dir, run_dir, o2_grid,
                                                     mode_reference_o2, seed_ids = NULL,
                                                     n_workers = 1L) {
  dat <- read_seed_objectives(
    analysis_dir = analysis_dir,
    fit_dir = run_dir,
    o2_values = o2_grid,
    seed_ids = seed_ids,
    n_workers = n_workers,
    mode_reference_o2 = mode_reference_o2
  )
  if (!nrow(dat)) return(data.frame())
  dat$seed_id <- o2ipa_norm_seed(dat$seed_id)
  if (!"objective" %in% names(dat)) dat$objective <- NA_real_
  dat$objective <- suppressWarnings(as.numeric(dat$objective))
  if (!"objective_source" %in% names(dat)) dat$objective_source <- NA_character_
  if (!"mode_label" %in% names(dat)) dat$mode_label <- NA_character_
  if (!"trajectory_regime" %in% names(dat)) dat$trajectory_regime <- NA_character_
  dat <- dat[!duplicated(dat$seed_id), , drop = FALSE]
  if (!is.null(seed_ids) && length(seed_ids)) dat <- dat[dat$seed_id %in% o2ipa_norm_seed(seed_ids), , drop = FALSE]
  dat
}


vf2_default_time_grid <- function() {
  sort(unique(c(seq(0, 100, by = 1), 125, 150, 175, 200, 250, 300, 400, 500, 700, 1000)))
}


vf2_normalize_state <- function(x) {
  x <- as.numeric(Re(x))
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  x <- pmax(x, 0)
  s <- sum(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  x / s
}


vf2_init_vector <- function(ngrid, init_N) {
  idx <- which.min(abs(ngrid - init_N))
  v <- numeric(length(ngrid))
  v[[idx]] <- 1
  list(vector = v, used_N = ngrid[[idx]], used_ploidy = ngrid[[idx]] / 22)
}


vf2_state_metrics <- function(state, ngrid, n_unit) {
  w <- as.numeric(Re(state))
  w[!is.finite(w)] <- 0
  w <- pmax(w, 0)
  total <- sum(w, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(data.frame(
      mean_N = NA_real_,
      mean_ploidy = NA_real_,
      sd_ploidy = NA_real_,
      fraction_N_le_25 = NA_real_,
      fraction_N_below_44 = NA_real_,
      fraction_N_ge_44 = NA_real_,
      fraction_N_ge_66 = NA_real_,
      fraction_N_ge_88 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  w <- w / total
  ploidy_grid <- ngrid / n_unit
  mean_N <- sum(ngrid * w, na.rm = TRUE)
  mean_ploidy <- sum(ploidy_grid * w, na.rm = TRUE)
  second_ploidy <- sum(ploidy_grid^2 * w, na.rm = TRUE)
  sd_ploidy <- sqrt(max(0, second_ploidy - mean_ploidy^2))
  data.frame(
    mean_N = mean_N,
    mean_ploidy = mean_ploidy,
    sd_ploidy = sd_ploidy,
    fraction_N_le_25 = sum(w[ngrid <= 25], na.rm = TRUE),
    fraction_N_below_44 = sum(w[ngrid < 44], na.rm = TRUE),
    fraction_N_ge_44 = sum(w[ngrid >= 44], na.rm = TRUE),
    fraction_N_ge_66 = sum(w[ngrid >= 66], na.rm = TRUE),
    fraction_N_ge_88 = sum(w[ngrid >= 88], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}


vf2_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22L), as.integer(cfg$N_MAX %||% 154L))
  G <- o2pr_build_G(model_env, cfg, run_params, O2)
  mu_all <- as.numeric(o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(ngrid)), run_params = run_params, N = ngrid))
  M <- G - Matrix::Diagonal(x = mu_all)
  list(M = M, G = G, mu_all = mu_all, ngrid = ngrid)
}


vf2_eigen_states <- function(M, init, time_grid) {
  Mdense <- as.matrix(M)
  eig <- tryCatch(eigen(Mdense, only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) stop("eigen decomposition failed")
  coef <- tryCatch(solve(eig$vectors, init), error = function(e) NULL)
  if (is.null(coef)) stop("eigen coefficient solve failed")
  lambda_ref <- max(Re(eig$values), na.rm = TRUE)
  states <- vector("list", length(time_grid))
  for (i in seq_along(time_grid)) {
    tt <- time_grid[[i]]
    weights <- exp((eig$values - lambda_ref) * tt) * coef
    states[[i]] <- vf2_normalize_state(eig$vectors %*% weights)
  }
  list(states = states, lambda_ref = lambda_ref)
}


vf2_expm_states <- function(M, init, time_grid) {
  time_grid <- sort(unique(as.numeric(time_grid)))
  if (!length(time_grid) || time_grid[[1]] < -1e-12) stop("time_grid must start at a non-negative time")
  x <- as.numeric(init)
  t_now <- 0
  states <- vector("list", length(time_grid))
  expm_cache <- new.env(parent = emptyenv())
  get_step_expm <- function(delta) {
    key <- format(signif(delta, 15), scientific = FALSE, trim = TRUE)
    if (!exists(key, envir = expm_cache, inherits = FALSE)) {
      assign(key, Matrix::expm(M * delta), envir = expm_cache)
    }
    get(key, envir = expm_cache, inherits = FALSE)
  }
  for (i in seq_along(time_grid)) {
    target <- time_grid[[i]]
    if (target < t_now - 1e-12) stop("time_grid must be sorted")
    delta <- target - t_now
    if (delta > 1e-12) {
      x <- as.numeric(get_step_expm(delta) %*% x)
      scale <- max(abs(x), na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) {
        x[] <- NA_real_
      } else {
        x <- x / scale
      }
      t_now <- target
    }
    states[[i]] <- vf2_normalize_state(x)
  }
  states
}


vf2_euler_states <- function(M, init, time_grid, dt) {
  if (!is.finite(dt) || dt <= 0) stop("dt must be positive")
  time_grid <- sort(unique(as.numeric(time_grid)))
  x <- as.numeric(init)
  t_now <- 0
  states <- vector("list", length(time_grid))
  for (i in seq_along(time_grid)) {
    target <- time_grid[[i]]
    if (target < t_now - 1e-12) stop("time_grid must be sorted")
    while (t_now < target - 1e-12) {
      h <- min(dt, target - t_now)
      x <- x + h * as.numeric(M %*% x)
      scale <- max(abs(x), na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) {
        x[] <- NA_real_
        break
      }
      x <- x / scale
      t_now <- t_now + h
    }
    states[[i]] <- vf2_normalize_state(x)
  }
  states
}


vf2_compare_states <- function(eigen_states, expm_states, euler_states, ngrid, n_unit, time_grid) {
  rows <- vector("list", length(time_grid))
  for (i in seq_along(time_grid)) {
    e_state <- eigen_states[[i]]
    x_state <- expm_states[[i]]
    u_state <- euler_states[[i]]
    em <- vf2_state_metrics(e_state, ngrid, n_unit)
    xm <- vf2_state_metrics(x_state, ngrid, n_unit)
    um <- vf2_state_metrics(u_state, ngrid, n_unit)
    rows[[i]] <- data.frame(
      day = time_grid[[i]],
      eigen_mean_N = em$mean_N,
      expm_mean_N = xm$mean_N,
      euler_mean_N = um$mean_N,
      diff_euler_expm_mean_N = um$mean_N - xm$mean_N,
      abs_diff_euler_expm_mean_N = abs(um$mean_N - xm$mean_N),
      diff_eigen_expm_mean_N = em$mean_N - xm$mean_N,
      abs_diff_eigen_expm_mean_N = abs(em$mean_N - xm$mean_N),
      eigen_mean_ploidy = em$mean_ploidy,
      expm_mean_ploidy = xm$mean_ploidy,
      euler_mean_ploidy = um$mean_ploidy,
      diff_euler_expm_mean_ploidy = um$mean_ploidy - xm$mean_ploidy,
      abs_diff_euler_expm_mean_ploidy = abs(um$mean_ploidy - xm$mean_ploidy),
      diff_eigen_expm_mean_ploidy = em$mean_ploidy - xm$mean_ploidy,
      abs_diff_eigen_expm_mean_ploidy = abs(em$mean_ploidy - xm$mean_ploidy),
      eigen_sd_ploidy = em$sd_ploidy,
      expm_sd_ploidy = xm$sd_ploidy,
      euler_sd_ploidy = um$sd_ploidy,
      diff_euler_expm_sd_ploidy = um$sd_ploidy - xm$sd_ploidy,
      abs_diff_euler_expm_sd_ploidy = abs(um$sd_ploidy - xm$sd_ploidy),
      diff_eigen_expm_sd_ploidy = em$sd_ploidy - xm$sd_ploidy,
      abs_diff_eigen_expm_sd_ploidy = abs(em$sd_ploidy - xm$sd_ploidy),
      eigen_fraction_N_le_25 = em$fraction_N_le_25,
      expm_fraction_N_le_25 = xm$fraction_N_le_25,
      euler_fraction_N_le_25 = um$fraction_N_le_25,
      eigen_fraction_N_below_44 = em$fraction_N_below_44,
      expm_fraction_N_below_44 = xm$fraction_N_below_44,
      euler_fraction_N_below_44 = um$fraction_N_below_44,
      eigen_fraction_N_ge_44 = em$fraction_N_ge_44,
      expm_fraction_N_ge_44 = xm$fraction_N_ge_44,
      euler_fraction_N_ge_44 = um$fraction_N_ge_44,
      eigen_fraction_N_ge_66 = em$fraction_N_ge_66,
      expm_fraction_N_ge_66 = xm$fraction_N_ge_66,
      euler_fraction_N_ge_66 = um$fraction_N_ge_66,
      eigen_fraction_N_ge_88 = em$fraction_N_ge_88,
      expm_fraction_N_ge_88 = xm$fraction_N_ge_88,
      euler_fraction_N_ge_88 = um$fraction_N_ge_88,
      max_abs_state_diff_euler_expm = max(abs(u_state - x_state), na.rm = TRUE),
      l1_state_diff_euler_expm = sum(abs(u_state - x_state), na.rm = TRUE),
      max_abs_state_diff_eigen_expm = max(abs(e_state - x_state), na.rm = TRUE),
      l1_state_diff_eigen_expm = sum(abs(e_state - x_state), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


vf2_error_summary_one <- function(comp) {
  terminal <- comp[which.max(comp$day), , drop = FALSE]
  data.frame(
    n_timepoint = nrow(comp),
    max_abs_diff_euler_expm_mean_ploidy = max(comp$abs_diff_euler_expm_mean_ploidy, na.rm = TRUE),
    terminal_abs_diff_euler_expm_mean_ploidy = terminal$abs_diff_euler_expm_mean_ploidy[[1]],
    max_abs_diff_eigen_expm_mean_ploidy = max(comp$abs_diff_eigen_expm_mean_ploidy, na.rm = TRUE),
    terminal_abs_diff_eigen_expm_mean_ploidy = terminal$abs_diff_eigen_expm_mean_ploidy[[1]],
    max_abs_diff_euler_expm_sd_ploidy = max(comp$abs_diff_euler_expm_sd_ploidy, na.rm = TRUE),
    terminal_abs_diff_euler_expm_sd_ploidy = terminal$abs_diff_euler_expm_sd_ploidy[[1]],
    max_abs_diff_eigen_expm_sd_ploidy = max(comp$abs_diff_eigen_expm_sd_ploidy, na.rm = TRUE),
    terminal_abs_diff_eigen_expm_sd_ploidy = terminal$abs_diff_eigen_expm_sd_ploidy[[1]],
    max_abs_diff_euler_expm_mean_N = max(comp$abs_diff_euler_expm_mean_N, na.rm = TRUE),
    terminal_abs_diff_euler_expm_mean_N = terminal$abs_diff_euler_expm_mean_N[[1]],
    max_abs_state_diff_euler_expm = max(comp$max_abs_state_diff_euler_expm, na.rm = TRUE),
    terminal_max_abs_state_diff_euler_expm = terminal$max_abs_state_diff_euler_expm[[1]],
    max_l1_state_diff_euler_expm = max(comp$l1_state_diff_euler_expm, na.rm = TRUE),
    terminal_l1_state_diff_euler_expm = terminal$l1_state_diff_euler_expm[[1]],
    stringsAsFactors = FALSE
  )
}


vf2_read_simulation_state_mean <- function(path) {
  if (!file.exists(path)) stop("Simulation state file not found: ", path)
  header_con <- gzfile(path, open = "rt")
  hdr <- strsplit(readLines(header_con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
  close(header_con)
  needed <- c("simulation_id", "day", "N", "ploidy", "status", "cell_count")
  missing <- setdiff(needed, hdr)
  if (length(missing)) stop("Simulation state file missing column(s): ", paste(missing, collapse = ", "), " in ", path)
  classes <- rep("NULL", length(hdr))
  names(classes) <- hdr
  classes["simulation_id"] <- "integer"
  classes["day"] <- "numeric"
  classes["N"] <- "integer"
  classes["ploidy"] <- "numeric"
  classes["status"] <- "character"
  classes["cell_count"] <- "numeric"
  data_con <- gzfile(path, open = "rt")
  on.exit(close(data_con), add = TRUE)
  tab <- utils::read.delim(data_con, check.names = FALSE, stringsAsFactors = FALSE, colClasses = unname(classes))
  live <- tab[tab$status == "live", , drop = FALSE]
  if (!nrow(live)) stop("No live rows found in simulation state file: ", path)
  split_live <- split(live, live$day)
  rows <- lapply(split_live, function(sub) {
    total <- sum(sub$cell_count, na.rm = TRUE)
    mean_ploidy <- sum(sub$ploidy * sub$cell_count, na.rm = TRUE) / total
    second_ploidy <- sum(sub$ploidy^2 * sub$cell_count, na.rm = TRUE) / total
    data.frame(
      day = as.numeric(sub$day[[1]]),
      simulation_mean_N = sum(sub$N * sub$cell_count, na.rm = TRUE) / total,
      simulation_mean_ploidy = mean_ploidy,
      simulation_sd_ploidy = sqrt(max(0, second_ploidy - mean_ploidy^2)),
      simulation_live_cells = total,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$day), , drop = FALSE]
}


vf2_read_simulation_trajectories <- function(sim_root, simulation_mode, seed_ids, seed_label_map, seed_mode_map, o2_grid, init_specs, simulation_ids) {
  rows <- list()
  k <- 0L
  for (seed_id in seed_ids) {
    for (O2 in o2_grid) {
      o2_tag <- vf2_num_path_tag(O2)
      for (j in seq_len(nrow(init_specs))) {
        ploidy_tag <- vf2_num_path_tag(init_specs$initial_ploidy[[j]])
        for (sim_id in simulation_ids) {
          state_path <- file.path(
            sim_root,
            simulation_mode,
            paste0("O2_", o2_tag),
            paste0("ploidy", ploidy_tag),
            seed_id,
            paste0("sim", sim_id),
            "state_trajectory.tsv.gz"
          )
          sim <- vf2_read_simulation_state_mean(state_path)
          sim$seed_id <- seed_id
          sim$seed_label <- seed_label_map[[seed_id]]
          sim$seed_mode <- seed_mode_map[[seed_id]]
          sim$O2_pct <- O2
          sim$initial_condition <- init_specs$initial_condition[[j]]
          sim$initial_ploidy <- init_specs$initial_ploidy[[j]]
          sim$simulation_id <- sim_id
          sim$state_file <- state_path
          k <- k + 1L
          rows[[k]] <- sim[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy",
            "simulation_id", "day", "simulation_mean_N", "simulation_mean_ploidy", "simulation_sd_ploidy",
            "simulation_live_cells", "state_file"
          )]
        }
      }
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}


vf2_simulation_summary <- function(sim_traj) {
  if (is.null(sim_traj) || !nrow(sim_traj)) return(data.frame())
  keys <- c("seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy", "day")
  split_key <- interaction(sim_traj[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(sim_traj, split_key), function(sub) {
    data.frame(
      seed_id = sub$seed_id[[1]],
      seed_mode = sub$seed_mode[[1]],
      seed_label = sub$seed_label[[1]],
      O2_pct = sub$O2_pct[[1]],
      initial_condition = sub$initial_condition[[1]],
      initial_ploidy = sub$initial_ploidy[[1]],
      day = sub$day[[1]],
      simulation_n = length(unique(sub$simulation_id)),
      simulation_mean_mean_N = mean(sub$simulation_mean_N, na.rm = TRUE),
      simulation_sd_mean_N = stats::sd(sub$simulation_mean_N, na.rm = TRUE),
      simulation_min_mean_N = min(sub$simulation_mean_N, na.rm = TRUE),
      simulation_max_mean_N = max(sub$simulation_mean_N, na.rm = TRUE),
      simulation_mean_mean_ploidy = mean(sub$simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_mean_ploidy = stats::sd(sub$simulation_mean_ploidy, na.rm = TRUE),
      simulation_min_mean_ploidy = min(sub$simulation_mean_ploidy, na.rm = TRUE),
      simulation_max_mean_ploidy = max(sub$simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(sub$simulation_sd_ploidy, na.rm = TRUE),
      simulation_sd_sd_ploidy = stats::sd(sub$simulation_sd_ploidy, na.rm = TRUE),
      simulation_min_sd_ploidy = min(sub$simulation_sd_ploidy, na.rm = TRUE),
      simulation_max_sd_ploidy = max(sub$simulation_sd_ploidy, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$seed_id, out$O2_pct, out$initial_condition, out$day), , drop = FALSE]
}


vf2_compare_solution_to_simulation <- function(solution_traj, sim_summary, solution_name) {
  if (is.null(sim_summary) || !nrow(sim_summary)) return(data.frame())
  mean_N_col <- paste0(solution_name, "_mean_N")
  mean_ploidy_col <- paste0(solution_name, "_mean_ploidy")
  sd_ploidy_col <- paste0(solution_name, "_sd_ploidy")
  if (!all(c(mean_N_col, mean_ploidy_col, sd_ploidy_col) %in% names(solution_traj))) {
    stop("solution_traj missing columns for method: ", solution_name)
  }
  sol <- solution_traj[, c(
    "seed_id", "seed_label", "O2_pct", "initial_condition", "day",
    mean_N_col, mean_ploidy_col, sd_ploidy_col
  ), drop = FALSE]
  names(sol)[names(sol) == mean_N_col] <- "analytical_mean_N"
  names(sol)[names(sol) == mean_ploidy_col] <- "analytical_mean_ploidy"
  names(sol)[names(sol) == sd_ploidy_col] <- "analytical_sd_ploidy"
  comp <- merge(
    sim_summary,
    sol,
    by = c("seed_id", "seed_label", "O2_pct", "initial_condition", "day"),
    all = FALSE,
    sort = FALSE
  )
  comp$solution_method <- solution_name
  comp$diff_simulation_analytical_mean_ploidy <- comp$simulation_mean_mean_ploidy - comp$analytical_mean_ploidy
  comp$abs_diff_simulation_analytical_mean_ploidy <- abs(comp$diff_simulation_analytical_mean_ploidy)
  comp$diff_simulation_analytical_sd_ploidy <- comp$simulation_mean_sd_ploidy - comp$analytical_sd_ploidy
  comp$abs_diff_simulation_analytical_sd_ploidy <- abs(comp$diff_simulation_analytical_sd_ploidy)
  comp$diff_simulation_analytical_mean_N <- comp$simulation_mean_mean_N - comp$analytical_mean_N
  comp$abs_diff_simulation_analytical_mean_N <- abs(comp$diff_simulation_analytical_mean_N)
  comp
}


vf2_simulation_comparison_summary <- function(comp) {
  if (is.null(comp) || !nrow(comp)) return(data.frame())
  split_key <- interaction(comp$seed_id, comp$O2_pct, comp$initial_condition, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(comp, split_key), function(sub) {
    terminal <- sub[which.max(sub$day), , drop = FALSE]
    sd_vals <- sub$simulation_sd_mean_ploidy
    max_sd <- if (all(is.na(sd_vals))) 0 else max(sd_vals, na.rm = TRUE)
    data.frame(
      seed_id = sub$seed_id[[1]],
      seed_label = sub$seed_label[[1]],
      O2_pct = sub$O2_pct[[1]],
      initial_condition = sub$initial_condition[[1]],
      solution_method = sub$solution_method[[1]],
      n_timepoint = nrow(sub),
      simulation_n = max(sub$simulation_n, na.rm = TRUE),
      max_abs_diff_simulation_analytical_mean_ploidy = max(sub$abs_diff_simulation_analytical_mean_ploidy, na.rm = TRUE),
      terminal_abs_diff_simulation_analytical_mean_ploidy = terminal$abs_diff_simulation_analytical_mean_ploidy[[1]],
      max_abs_diff_simulation_analytical_sd_ploidy = max(sub$abs_diff_simulation_analytical_sd_ploidy, na.rm = TRUE),
      terminal_abs_diff_simulation_analytical_sd_ploidy = terminal$abs_diff_simulation_analytical_sd_ploidy[[1]],
      max_abs_diff_simulation_analytical_mean_N = max(sub$abs_diff_simulation_analytical_mean_N, na.rm = TRUE),
      terminal_abs_diff_simulation_analytical_mean_N = terminal$abs_diff_simulation_analytical_mean_N[[1]],
      max_simulation_sd_mean_ploidy = max_sd,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


fixo2_simulation_state_path <- function(simulation_dir, simulation_mode, O2, initial_ploidy, seed_id, simulation_id) {
  file.path(
    simulation_dir,
    simulation_mode,
    paste0("O2_", vf2_num_path_tag(O2)),
    paste0("ploidy", vf2_num_path_tag(initial_ploidy)),
    seed_id,
    paste0("sim", simulation_id),
    "state_trajectory.tsv.gz"
  )
}


fixo2_expected_simulation_tasks <- function(run_dir, simulation_dir, simulation_mode, seed_ids, o2_grid, init_specs, simulation_ids) {
  rows <- list()
  k <- 0L
  for (seed_id in seed_ids) {
    seed_dir <- file.path(run_dir, seed_id)
    for (O2 in o2_grid) {
      for (j in seq_len(nrow(init_specs))) {
        initial_ploidy <- init_specs$initial_ploidy[[j]]
        leaf_dir <- dirname(dirname(fixo2_simulation_state_path(simulation_dir, simulation_mode, O2, initial_ploidy, seed_id, 1L)))
        for (sim_id in simulation_ids) {
          state_file <- fixo2_simulation_state_path(simulation_dir, simulation_mode, O2, initial_ploidy, seed_id, sim_id)
          k <- k + 1L
          rows[[k]] <- data.frame(
            task_id = k,
            seed_id = seed_id,
            seed_dir = seed_dir,
            O2_pct = O2,
            initial_ploidy = initial_ploidy,
            initial_condition = init_specs$initial_condition[[j]],
            simulation_id = as.integer(sim_id),
            output_dir = dirname(state_file),
            state_file = state_file,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    info <- file.info(out$state_file)
    out$complete_before <- file.exists(out$state_file) & !is.na(info$size) & info$size > 0
  }
  out
}


fixo2_simulation_forward_args <- function(args) {
  keys <- c("dt", "report_dt", "initial_cells", "joint_scope", "Crowding", "O2_growth", "start_with", "ploidy_O2_death")
  out <- character(0)
  for (key in keys) {
    if (!is.null(args[[key]]) && length(args[[key]]) > 0L) {
      out <- c(out, paste0("--", key, "=", as.character(args[[key]][[1]])))
    }
  }
  if (!is.null(args$save_every_days) && length(args$save_every_days) > 0L) {
    out <- c(out, paste0("--save_every_days=", as.character(args$save_every_days[[1]])))
  } else if (!is.null(args$time_step_days) && length(args$time_step_days) > 0L) {
    out <- c(out, paste0("--save_every_days=", as.character(args$time_step_days[[1]])))
  }
  out
}


fixo2_simulation_n_core <- function(args, task_count = Inf) {
  raw <- args$simulation_n_core %||% args$simulation_workers %||% args$n_core
  if (is.null(raw) || !length(raw)) raw <- 1L
  n_core <- o2ipa_as_int(raw, 1L)
  if (!is.finite(n_core) || is.na(n_core) || n_core < 1L) n_core <- 1L
  if (is.finite(task_count)) n_core <- min(as.integer(n_core), as.integer(task_count))
  max(1L, as.integer(n_core))
}


fixo2_simulation_worker_threads <- function(args) {
  raw <- args$simulation_worker_threads %||% args$simulation_threads_per_worker
  if (is.null(raw) || !length(raw)) raw <- 1L
  n_threads <- o2ipa_as_int(raw, 1L)
  if (!is.finite(n_threads) || is.na(n_threads) || n_threads < 1L) n_threads <- 1L
  as.integer(n_threads)
}


fixo2_simulation_worker_env <- function(args) {
  n_threads <- fixo2_simulation_worker_threads(args)
  c(
    paste0(
      c(
        "OMP_NUM_THREADS",
        "OMP_THREAD_LIMIT",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS"
      ),
      "=",
      n_threads
    ),
    "OMP_WAIT_POLICY=PASSIVE"
  )
}


fixo2_task_log_completed <- function(log_file) {
  if (!file.exists(log_file)) return(NA)
  lines <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) character(0))
  if (!length(lines)) return(FALSE)
  any(grepl("^Done\\.$", lines))
}


fixo2_run_missing_simulation_task <- function(idx, tasks, rscript, sim_script, simulation_mode, time_days, extra_args, out_dir, worker_env) {
  task <- tasks[idx, , drop = FALSE]
  log_file <- file.path(out_dir, "logs", sprintf("fix_o2_simulation_task_%04d.log", task$task_id[[1]]))
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  cmd_args <- c(
    sim_script,
    paste0("--fit_dir=", task$seed_dir[[1]]),
    paste0("--simulation=", simulation_mode),
    paste0("--initial_ploidy=", task$initial_ploidy[[1]]),
    paste0("--time_days=", time_days),
    "--n_sim=1",
    paste0("--simulation_id=", task$simulation_id[[1]]),
    paste0("--o2=", task$O2_pct[[1]]),
    paste0("--seed_id=", task$seed_id[[1]]),
    paste0("--out_dir=", task$output_dir[[1]]),
    extra_args
  )
  error_message <- ""
  status <- tryCatch(
    system2(rscript, cmd_args, stdout = log_file, stderr = log_file, env = worker_env),
    error = function(e) {
      error_message <<- conditionMessage(e)
      cat("ERROR launching fixed-O2 simulation task: ", error_message, "\n", file = log_file, append = TRUE, sep = "")
      1L
    }
  )
  status <- if (is.null(status)) 0L else as.integer(status)
  data.frame(
    task_index = as.integer(idx),
    task_id = as.integer(task$task_id[[1]]),
    run_status = status,
    generated = TRUE,
    log_file = log_file,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}


fixo2_ensure_simulation <- function(args, run_dir, simulation_dir, simulation_mode, seed_ids, o2_grid, init_specs, simulation_ids, out_dir) {
  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
  tasks <- fixo2_expected_simulation_tasks(run_dir, simulation_dir, simulation_mode, seed_ids, o2_grid, init_specs, simulation_ids)
  if (!nrow(tasks)) return(tasks)

  generate_missing <- o2ipa_as_bool(args$generate_missing_simulation, TRUE)
  time_days <- o2ipa_as_num(args$time_days, 1000)
  sim_script <- normalizePath(file.path(SCRIPT_DIR, "..", "..", "simulation", "fix_o2_simulation.R"), mustWork = FALSE)
  if (!file.exists(sim_script)) stop("Simulation script was not found: ", sim_script)
  rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
  tasks$run_status <- NA_integer_
  tasks$generated <- FALSE
  tasks$log_file <- file.path(out_dir, "logs", sprintf("fix_o2_simulation_task_%04d.log", tasks$task_id))
  tasks$log_completed_before <- vapply(tasks$log_file, fixo2_task_log_completed, logical(1))
  tasks$error_message <- ""
  tasks$complete_before <- tasks$complete_before & (is.na(tasks$log_completed_before) | tasks$log_completed_before)

  missing <- which(!tasks$complete_before)
  if (length(missing) && !generate_missing) {
    fixo2_write_tsv(tasks, file.path(out_dir, "tables", "fixed_o2_simulation_file_status.tsv"))
    stop("Missing ", length(missing), " fixed-O2 simulation file(s); set --generate_missing_simulation=TRUE to create them.")
  }

  if (length(missing)) {
    n_core <- fixo2_simulation_n_core(args, length(missing))
    worker_threads <- fixo2_simulation_worker_threads(args)
    worker_env <- fixo2_simulation_worker_env(args)
    message("Generating missing fixed-O2 simulation task(s): ", length(missing))
    message("simulation_n_core: ", n_core)
    message("simulation_worker_threads: ", worker_threads)
    extra_args <- fixo2_simulation_forward_args(args)
    bad_seed_dir <- missing[!dir.exists(tasks$seed_dir[missing])]
    if (length(bad_seed_dir)) {
      stop("Seed directory does not exist for simulation task(s): ", paste(tasks$task_id[bad_seed_dir], collapse = ","))
    }
    runner <- function(idx) {
      tryCatch(
        fixo2_run_missing_simulation_task(
          idx = idx,
          tasks = tasks,
          rscript = rscript,
          sim_script = sim_script,
          simulation_mode = simulation_mode,
          time_days = time_days,
          extra_args = extra_args,
          out_dir = out_dir,
          worker_env = worker_env
        ),
        error = function(e) {
          log_file <- file.path(out_dir, "logs", sprintf("fix_o2_simulation_task_%04d.log", tasks$task_id[[idx]]))
          cat("ERROR running fixed-O2 simulation task: ", conditionMessage(e), "\n", file = log_file, append = TRUE, sep = "")
          data.frame(
            task_index = as.integer(idx),
            task_id = as.integer(tasks$task_id[[idx]]),
            run_status = 1L,
            generated = FALSE,
            log_file = log_file,
            error_message = conditionMessage(e),
            stringsAsFactors = FALSE
          )
        }
      )
    }
    if (n_core > 1L && .Platform$OS.type != "windows") {
      run_results <- parallel::mclapply(missing, runner, mc.cores = n_core, mc.preschedule = FALSE)
    } else {
      if (n_core > 1L) warning("Parallel fixed-O2 simulation generation is not supported on Windows; falling back to sequential.")
      run_results <- lapply(missing, runner)
    }
    run_results <- do.call(rbind, run_results)
    tasks$run_status[run_results$task_index] <- run_results$run_status
    tasks$generated[run_results$task_index] <- run_results$generated
    tasks$log_file[run_results$task_index] <- run_results$log_file
    tasks$error_message[run_results$task_index] <- run_results$error_message
    failed <- run_results[run_results$run_status != 0L | nzchar(run_results$error_message), , drop = FALSE]
    if (nrow(failed)) {
      fixo2_write_tsv(tasks, file.path(out_dir, "tables", "fixed_o2_simulation_file_status.tsv"))
      stop(
        "Fixed-O2 simulation task(s) failed: task_id=",
        paste(failed$task_id, collapse = ","),
        "; log=",
        paste(failed$log_file, collapse = ",")
      )
    }
  } else {
    message("All requested fixed-O2 simulation files already exist; skipping simulation generation.")
  }

  info_after <- file.info(tasks$state_file)
  tasks$complete_after <- file.exists(tasks$state_file) & !is.na(info_after$size) & info_after$size > 0
  fixo2_write_tsv(tasks, file.path(out_dir, "tables", "fixed_o2_simulation_file_status.tsv"))
  if (any(!tasks$complete_after)) {
    failed <- tasks[!tasks$complete_after, , drop = FALSE]
    stop("Missing fixed-O2 simulation output(s) after generation: ", paste(failed$task_id, collapse = ","))
  }
  invisible(tasks)
}


analytical_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  fm <- fixo2_fixed_matrix(model_env, cfg, run_params, O2)
  list(M = fm$M, ngrid = fm$ngrid)
}


analytical_init_vector <- function(ngrid, init_N) {
  fixo2_init_vector(ngrid, init_N, n_unit = 22)
}


analytical_normalize_state <- function(x) {
  fixo2_normalize_state(x)
}


analytical_state_metrics <- function(state, ngrid, n_unit) {
  w <- analytical_normalize_state(state)
  if (all(is.na(w))) {
    return(data.frame(
      analytical_mean_N = NA_real_,
      analytical_mean_ploidy = NA_real_,
      analytical_sd_ploidy = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  ploidy_grid <- ngrid / n_unit
  mean_N <- sum(ngrid * w, na.rm = TRUE)
  mean_ploidy <- sum(ploidy_grid * w, na.rm = TRUE)
  second_ploidy <- sum(ploidy_grid^2 * w, na.rm = TRUE)
  data.frame(
    analytical_mean_N = mean_N,
    analytical_mean_ploidy = mean_ploidy,
    analytical_sd_ploidy = sqrt(max(0, second_ploidy - mean_ploidy^2)),
    stringsAsFactors = FALSE
  )
}


analytical_eigen_states <- function(M, init, time_grid) {
  fixo2_eigen_states(M, init, time_grid)
}


analytical_expm_states <- function(M, init, time_grid) {
  fixo2_expm_states(M, init, time_grid)
}


generate_seed_analytical_trajectories <- function(seed_id, inputs, param_mat, model_env,
                                                  time_points, o2_values, initial_ploidy_values,
                                                  analytical_methods) {
  manifest <- inputs$manifest[inputs$manifest$seed_id == seed_id, , drop = FALSE]
  if (!nrow(manifest) || !seed_id %in% rownames(param_mat)) return(data.frame())
  cfg <- o2pr_first_seed_cfg(manifest)
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
  names(pvec) <- colnames(param_mat)
  run_params <- o2pr_run_params_from_vec(pvec, cfg)
  init_specs <- data.frame(
    initial_condition = paste0("init_", format(initial_ploidy_values, scientific = FALSE, trim = TRUE), "N"),
    initial_ploidy = initial_ploidy_values,
    requested_N = initial_ploidy_values * n_unit,
    stringsAsFactors = FALSE
  )

  rows <- list()
  k <- 0L
  for (O2 in o2_values) {
    fm <- analytical_fixed_matrix(model_env, cfg, run_params, O2)
    for (j in seq_len(nrow(init_specs))) {
      init <- analytical_init_vector(fm$ngrid, init_specs$requested_N[[j]])
      states_by_method <- list()
      if ("eigen" %in% analytical_methods) {
        states_by_method$eigen <- tryCatch(analytical_eigen_states(fm$M, init$vector, time_points), error = function(e) NULL)
      }
      if ("expm" %in% analytical_methods) {
        states_by_method$expm <- tryCatch(analytical_expm_states(fm$M, init$vector, time_points), error = function(e) NULL)
      }
      for (method in names(states_by_method)) {
        states <- states_by_method[[method]]
        if (is.null(states)) next
        for (i in seq_along(time_points)) {
          met <- analytical_state_metrics(states[[i]], fm$ngrid, n_unit)
          k <- k + 1L
          rows[[k]] <- data.frame(
            seed_id = seed_id,
            O2_pct = O2,
            O2_key = num_key(O2),
            initial_condition = init_specs$initial_condition[[j]],
            initial_ploidy = init_specs$initial_ploidy[[j]],
            requested_initial_N = init_specs$requested_N[[j]],
            used_initial_N = init$used_N,
            day = time_points[[i]],
            day_key = num_key(time_points[[i]]),
            analytical_method = method,
            analytical_method_label = method_label(method),
            met,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame()
  out
}


generate_analytical_trajectories <- function(run_dir, time_points, o2_values, initial_ploidy_values,
                                             analytical_methods, n_workers = 1L, seed_ids = NULL) {
  if (!dir.exists(run_dir)) stop("run_dir does not exist and is required to generate analytical trajectories: ", run_dir)
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  seeds <- if (is.null(seed_ids) || !length(seed_ids)) rownames(param_mat) else intersect(seed_ids, rownames(param_mat))
  if (!length(seeds)) stop("No seed parameters were found for analytical trajectory generation.")
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, length(seeds))
  message("Generating analytical trajectories from fitted parameters: ", length(seeds), " seeds, methods=", paste(analytical_methods, collapse = ","), ", workers=", n_workers)
  worker <- function(seed_id) {
    generate_seed_analytical_trajectories(
      seed_id = seed_id,
      inputs = inputs,
      param_mat = param_mat,
      model_env = model_env,
      time_points = time_points,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      analytical_methods = analytical_methods
    )
  }
  rows <- if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
    parallel::mclapply(seeds, worker, mc.cores = n_workers)
  } else {
    lapply(seeds, worker)
  }
  out <- do.call(rbind, rows[vapply(rows, nrow, integer(1)) > 0L])
  if (is.null(out)) out <- data.frame()
  if (nrow(out)) out$analytical_source_run_dir <- normalizePath(run_dir, mustWork = FALSE)
  out
}


filter_analytical_trajectories <- function(analytical, time_points, o2_values, initial_ploidy_values, analytical_methods, seed_ids = NULL) {
  if (!nrow(analytical)) return(analytical)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "analytical_mean_N", "analytical_mean_ploidy", "analytical_sd_ploidy"),
    names(analytical)
  )
  for (col in numeric_cols) analytical[[col]] <- suppressWarnings(as.numeric(analytical[[col]]))
  analytical$analytical_method <- vapply(analytical$analytical_method, method_slug, character(1))
  out <- filter_by_values(analytical, "day", time_points)
  out <- filter_by_values(out, "O2_pct", o2_values)
  out <- filter_by_values(out, "initial_ploidy", initial_ploidy_values)
  out <- out[out$analytical_method %in% analytical_methods, , drop = FALSE]
  if (!is.null(seed_ids) && length(seed_ids)) out <- out[out$seed_id %in% seed_ids, , drop = FALSE]
  out$O2_key <- num_key(out$O2_pct)
  out$day_key <- num_key(out$day)
  out
}


expected_analytical_seed_ids <- function(run_dir, seed_ids = NULL) {
  seed_ids <- normalize_seed_ids(seed_ids)
  if (length(seed_ids)) return(seed_ids)
  if (is.null(run_dir) || !nzchar(run_dir) || !dir.exists(run_dir)) return(character())
  seeds <- tryCatch(o2ipa_discover_seeds(run_dir)$seed_id, error = function(e) character())
  normalize_seed_ids(seeds)
}


analytical_cache_missing_keys <- function(analytical, time_points, o2_values, initial_ploidy_values,
                                          analytical_methods, expected_seed_ids = character(), expected_run_dir = NULL) {
  if (!nrow(analytical)) return("all")
  required_cols <- c("seed_id", "O2_pct", "initial_ploidy", "day", "analytical_method", "analytical_mean_ploidy")
  missing_cols <- setdiff(required_cols, names(analytical))
  if (length(missing_cols)) return(paste0("missing column: ", paste(missing_cols, collapse = ",")))

  if (!is.null(expected_run_dir) && nzchar(expected_run_dir)) {
    expected_run_dir <- normalizePath(expected_run_dir, mustWork = FALSE)
    if (!"analytical_source_run_dir" %in% names(analytical)) return("missing analytical source run_dir")
    observed_run_dirs <- unique(normalizePath(as.character(analytical$analytical_source_run_dir), mustWork = FALSE))
    observed_ok <- observed_run_dirs %in% expected_run_dir |
      path_results_suffix(observed_run_dirs) %in% path_results_suffix(expected_run_dir)
    if (!length(observed_run_dirs) || any(is.na(observed_run_dirs)) || !all(observed_ok)) {
      return(paste0("analytical source run_dir mismatch: ", paste(head(observed_run_dirs, 3L), collapse = ",")))
    }
  }

  seeds <- normalize_seed_ids(expected_seed_ids)
  if (!length(seeds)) seeds <- normalize_seed_ids(unique(analytical$seed_id))
  if (!length(seeds)) return("no seeds")

  cache <- data.frame(
    seed_id = normalize_seed_ids(analytical$seed_id),
    O2_key = num_key(analytical$O2_pct),
    initial_ploidy_key = num_key(analytical$initial_ploidy),
    day_key = num_key(analytical$day),
    analytical_method = vapply(analytical$analytical_method, method_slug, character(1)),
    analytical_mean_ploidy = suppressWarnings(as.numeric(analytical$analytical_mean_ploidy)),
    stringsAsFactors = FALSE
  )
  available_keys <- unique(paste(cache$seed_id, cache$O2_key, cache$initial_ploidy_key, cache$day_key, cache$analytical_method, sep = "\r"))
  required <- expand.grid(
    seed_id = seeds,
    O2_key = num_key(o2_values),
    initial_ploidy_key = num_key(initial_ploidy_values),
    day_key = num_key(time_points),
    analytical_method = analytical_methods,
    stringsAsFactors = FALSE
  )
  required_keys <- paste(required$seed_id, required$O2_key, required$initial_ploidy_key, required$day_key, required$analytical_method, sep = "\r")
  missing <- setdiff(required_keys, available_keys)
  if (!length(missing)) return(character())
  head(missing, 5L)
}


generate_fixo2_attractor_mode_table <- function(run_dir, o2_values, seed_ids = NULL, n_workers = 1L) {
  if (is.null(run_dir) || !nzchar(run_dir) || !dir.exists(run_dir)) {
    stop("run_dir is required to generate FixO2 attractor mode table: ", run_dir)
  }
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  seeds <- if (is.null(seed_ids) || !length(seed_ids)) rownames(param_mat) else intersect(normalize_seed_ids(seed_ids), rownames(param_mat))
  if (!length(seeds)) stop("No seed parameters were found for FixO2 attractor mode generation.")
  cfg <- o2pr_first_seed_cfg(inputs$manifest)
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  o2_values <- sort(unique(as.numeric(o2_values)))
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, length(seeds))
  message("Generating FixO2 attractor mode table: ", length(seeds), " seeds, ", length(o2_values), " O2 values, workers=", n_workers)
  worker <- function(seed) {
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    rows <- lapply(o2_values, function(O2) fo2_dominant_attractor_one(seed, run_params, model_env, cfg, O2))
    do.call(rbind, rows)
  }
  rows <- if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
    parallel::mclapply(seeds, worker, mc.cores = n_workers)
  } else {
    lapply(seeds, worker)
  }
  attractors <- do.call(rbind, rows[vapply(rows, nrow, integer(1)) > 0L])
  if (is.null(attractors)) attractors <- data.frame()
  if (!nrow(attractors)) return(data.frame())
  manifest <- inputs$manifest
  manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  attractors <- merge(
    attractors,
    manifest[, intersect(c("seed_id", "objective", "delta_objective"), names(manifest)), drop = FALSE],
    by = "seed_id",
    all.x = TRUE,
    sort = FALSE
  )
  attractors <- fixo2_assign_attractor_modes(attractors, "dominant_mean_ploidy")
  fixo2_attractor_mode_table(attractors)
}


read_seed_objectives <- function(analysis_dir, fit_dir = NULL, o2_values = NULL, seed_ids = NULL, n_workers = 1L, mode_reference_o2 = 2) {
  mode_reference_path <- file.path(
    analysis_dir,
    "attractors",
    "tables",
    "fixed_o2_attractor_mode_reference_by_seed.tsv"
  )
  mode_seed_path <- file.path(
    analysis_dir,
    "attractors",
    "tables",
    "fixed_o2_attractor_mode_by_seed.tsv"
  )
  mode_o2_path <- file.path(
    analysis_dir,
    "attractors",
    "tables",
    "fixed_o2_attractor_mode_by_seed_o2.tsv"
  )
  mode_tab <- data.frame(seed_id = character(), trajectory_regime = character(), mode_label = character(), stringsAsFactors = FALSE)
  if (file.exists(mode_reference_path) || file.exists(mode_seed_path)) {
    mode_path <- if (file.exists(mode_reference_path)) mode_reference_path else mode_seed_path
    tab <- read_tsv(mode_path)
    cols <- intersect(c(
      "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
      "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
      "mode_reference_dominant_mean_ploidy", "mode_reference_status",
      "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
    ), names(tab))
    if ("seed_id" %in% cols) {
      mode_tab <- tab[, cols, drop = FALSE]
      mode_tab <- mode_tab[!duplicated(mode_tab$seed_id), , drop = FALSE]
      if ("mode_reference_o2_pct" %in% names(mode_tab)) {
        ref_vals <- suppressWarnings(as.numeric(mode_tab$mode_reference_o2_pct))
        if (any(is.finite(ref_vals)) && !any(abs(ref_vals - mode_reference_o2) < 1e-9, na.rm = TRUE)) {
          warning("Existing reference mode table was generated for a different mode_reference_o2; regenerating if fit_dir is available.")
          mode_tab <- data.frame(seed_id = character(), trajectory_regime = character(), mode_label = character(), stringsAsFactors = FALSE)
        }
      } else {
        warning("Existing seed-level mode table does not record mode_reference_o2; regenerating if fit_dir is available.")
        mode_tab <- data.frame(seed_id = character(), trajectory_regime = character(), mode_label = character(), stringsAsFactors = FALSE)
      }
    }
  } else if (file.exists(mode_o2_path)) {
    tab <- read_tsv(mode_o2_path)
    if ("O2_pct" %in% names(tab)) {
      reference_tab <- tryCatch(
        fixo2_reference_mode_table(tab, mode_reference_o2),
        error = function(e) {
          warning(conditionMessage(e))
          data.frame()
        }
      )
      if (nrow(reference_tab)) {
        cols <- intersect(c(
          "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
          "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
          "mode_reference_dominant_mean_ploidy", "mode_reference_status",
          "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
        ), names(reference_tab))
        mode_tab <- reference_tab[, cols, drop = FALSE]
        dir.create(dirname(mode_reference_path), recursive = TRUE, showWarnings = FALSE)
        write_tsv(reference_tab, mode_reference_path)
        write_tsv(reference_tab, mode_seed_path)
      }
    }
  }
  if (!nrow(mode_tab) && !is.null(fit_dir) && nzchar(fit_dir) && dir.exists(fit_dir) && !is.null(o2_values) && length(o2_values)) {
    mode_by_seed_o2 <- generate_fixo2_attractor_mode_table(
      run_dir = fit_dir,
      o2_values = sort(unique(c(o2_values, mode_reference_o2))),
      seed_ids = seed_ids,
      n_workers = n_workers
    )
    mode_reference_tab <- fixo2_reference_mode_table(mode_by_seed_o2, mode_reference_o2)
    mode_cols <- intersect(c(
      "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
      "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
      "mode_reference_dominant_mean_ploidy", "mode_reference_status",
      "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
    ), names(mode_reference_tab))
    mode_tab <- mode_reference_tab[, mode_cols, drop = FALSE]
    dir.create(dirname(mode_o2_path), recursive = TRUE, showWarnings = FALSE)
    write_tsv(mode_by_seed_o2, mode_o2_path)
    write_tsv(mode_reference_tab, mode_reference_path)
    write_tsv(mode_reference_tab, mode_seed_path)
  }

  if (is.null(fit_dir) || !nzchar(fit_dir) || !dir.exists(fit_dir)) {
    warning("Final seed objective values require fit_dir, but fit_dir is unavailable. Objective-colored plots will have missing objective values.")
    mode_tab$objective <- NA_real_
    mode_tab$objective_source <- NA_character_
    return(mode_tab)
  }

  read_summary_objectives <- function(fit_dir) {
    summary_path <- file.path(fit_dir, "extra_results", "seed_summary.tsv")
    if (!file.exists(summary_path)) return(data.frame())
    tab <- tryCatch(read_tsv(summary_path), error = function(e) data.frame())
    if (!nrow(tab)) return(data.frame())
    if (!"seed_id" %in% names(tab)) {
      if ("seed" %in% names(tab)) {
        tab$seed_id <- normalize_seed_ids(tab$seed)
      } else {
        return(data.frame())
      }
    }
    final_cols <- intersect(
      c("objective_total", "objective", "optimizer_local_objective", "optimizer_deoptim_objective", "objective_data"),
      names(tab)
    )
    if (!length(final_cols)) return(data.frame())
    objective <- rep(NA_real_, nrow(tab))
    source <- rep(NA_character_, nrow(tab))
    for (col in final_cols) {
      vals <- suppressWarnings(as.numeric(tab[[col]]))
      fill <- !is.finite(objective) & is.finite(vals)
      objective[fill] <- vals[fill]
      source[fill] <- col
    }
    out <- data.frame(
      seed_id = normalize_seed_ids(tab$seed_id),
      objective = objective,
      objective_source = source,
      stringsAsFactors = FALSE
    )
    out <- out[is.finite(out$objective), , drop = FALSE]
    out[!duplicated(out$seed_id), , drop = FALSE]
  }

  objectives <- read_summary_objectives(fit_dir)
  if (nrow(objectives)) {
    out <- if (nrow(mode_tab)) merge(mode_tab, objectives, by = "seed_id", all = TRUE) else objectives
    return(out)
  }

  seed_dirs <- list.dirs(fit_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("seed[0-9]+$", basename(seed_dirs))]
  rows <- lapply(seed_dirs, function(seed_dir) {
    candidates <- c(
      file.path(seed_dir, "fit_summary.tsv"),
      file.path(seed_dir, "best_params.tsv"),
      file.path(seed_dir, "best_parameters.tsv")
    )
    hits <- candidates[file.exists(candidates)]
    path <- if (length(hits)) hits[[1]] else NA_character_
    if (is.na(path)) return(NULL)
    tab <- tryCatch(read_tsv(path), error = function(e) data.frame())
    if (!nrow(tab)) return(NULL)
    if (all(c("metric", "value") %in% names(tab))) {
      vals <- as.list(tab$value)
      names(vals) <- tab$metric
      objective_cols <- c("objective_total", "objective", "optimizer_local_objective", "optimizer_deoptim_objective", "objective_data")
      hit <- objective_cols[vapply(objective_cols, function(col) {
        is.finite(suppressWarnings(as.numeric(vals[[col]])))
      }, logical(1))]
      if (!length(hit)) return(NULL)
      objective <- suppressWarnings(as.numeric(vals[[hit[[1]]]]))
      source <- hit[[1]]
    } else {
      objective_cols <- intersect(
        c("objective_total", "objective", "optimizer_local_objective", "optimizer_deoptim_objective", "objective_data"),
        names(tab)
      )
      if (!length(objective_cols)) return(NULL)
      objective <- NA_real_
      source <- NA_character_
      for (col in objective_cols) {
        val <- suppressWarnings(as.numeric(tab[[col]][[1]]))
        if (is.finite(val)) {
          objective <- val
          source <- col
          break
        }
      }
      if (!is.finite(objective)) return(NULL)
    }
    data.frame(
      seed_id = basename(seed_dir),
      objective = objective,
      objective_source = source,
      stringsAsFactors = FALSE
    )
  })
  objectives <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(objectives)) objectives <- data.frame(seed_id = character(), objective = numeric(), objective_source = character(), stringsAsFactors = FALSE)
  if (nrow(mode_tab)) merge(mode_tab, objectives, by = "seed_id", all = TRUE) else objectives
}


task_table_path <- function(simulation_dir, simulation_mode) {
  file.path(simulation_dir, simulation_mode, "task_list.tsv")
}


build_task_table_from_paths <- function(simulation_dir, simulation_mode, seed_ids, o2_values, initial_ploidy_values, simulation_ids) {
  rows <- list()
  k <- 0L
  for (seed_id in seed_ids) {
    for (O2 in o2_values) {
      for (initial_ploidy in initial_ploidy_values) {
        for (simulation_id in simulation_ids) {
          k <- k + 1L
          output_dir <- file.path(
            simulation_dir,
            simulation_mode,
            paste0("O2_", num_path_tag(O2)),
            paste0("ploidy", num_path_tag(initial_ploidy)),
            seed_id,
            paste0("sim", simulation_id)
          )
          rows[[k]] <- data.frame(
            task_id = k,
            seed_id = seed_id,
            fixed_o2_pct = O2,
            initial_ploidy = initial_ploidy,
            initial_condition = initial_condition_from_ploidy(initial_ploidy),
            simulation_id = simulation_id,
            output_dir = output_dir,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}


read_simulation_tasks <- function(simulation_dir, simulation_mode, analytical, o2_values, initial_ploidy_values, simulation_ids) {
  path <- task_table_path(simulation_dir, simulation_mode)
  if (file.exists(path)) {
    tasks <- read_tsv(path)
    if (!"fixed_o2_pct" %in% names(tasks) && "O2_pct" %in% names(tasks)) {
      names(tasks)[names(tasks) == "O2_pct"] <- "fixed_o2_pct"
    }
    required <- c("seed_id", "fixed_o2_pct", "initial_ploidy", "simulation_id", "output_dir")
    missing <- setdiff(required, names(tasks))
    if (length(missing)) stop("Simulation task table is missing column(s): ", paste(missing, collapse = ", "))
    if (!"initial_condition" %in% names(tasks)) {
      tasks$initial_condition <- initial_condition_from_ploidy(tasks$initial_ploidy)
    }
  } else {
    seed_ids <- sort(unique(analytical$seed_id))
    tasks <- build_task_table_from_paths(
      simulation_dir = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      simulation_ids = simulation_ids
    )
  }
  tasks$fixed_o2_pct <- as.numeric(tasks$fixed_o2_pct)
  tasks$initial_ploidy <- as.numeric(tasks$initial_ploidy)
  tasks$simulation_id <- as.integer(tasks$simulation_id)
  expected_output_dir <- file.path(
    simulation_dir,
    simulation_mode,
    paste0("O2_", vapply(tasks$fixed_o2_pct, num_path_tag, character(1))),
    paste0("ploidy", vapply(tasks$initial_ploidy, num_path_tag, character(1))),
    tasks$seed_id,
    paste0("sim", tasks$simulation_id)
  )
  use_expected <- !dir.exists(tasks$output_dir) & dir.exists(expected_output_dir)
  tasks$output_dir[use_expected] <- expected_output_dir[use_expected]
  tasks <- filter_by_values(tasks, "fixed_o2_pct", o2_values)
  tasks <- filter_by_values(tasks, "initial_ploidy", initial_ploidy_values)
  tasks <- tasks[tasks$simulation_id %in% simulation_ids, , drop = FALSE]
  tasks <- tasks[tasks$seed_id %in% unique(analytical$seed_id), , drop = FALSE]
  tasks$state_file <- file.path(tasks$output_dir, "state_trajectory.tsv.gz")
  tasks
}


read_state_metrics_awk <- function(path, time_points = NULL) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    return(data.frame())
  }
  gzip <- Sys.which("gzip")
  awk <- Sys.which("awk")
  if (!nzchar(gzip) || !nzchar(awk)) {
    stop("Reading state trajectories requires gzip and awk on PATH.")
  }
  all_days <- is.null(time_points) || !length(time_points)
  days_arg <- if (all_days) "" else paste(num_key(time_points), collapse = ",")
  awk_script <- paste(
    'BEGIN {',
    '  keep_all = all_days + 0;',
    '  if (!keep_all) {',
    '    split(days, d, ",");',
    '    for (i in d) keep[sprintf("%.10g", d[i] + 0)] = 1;',
    '  }',
    '  OFS = "\t";',
    '}',
    'NR == 1 {',
    '  for (i = 1; i <= NF; i++) idx[$i] = i;',
    '  next;',
    '}',
    '{',
    '  day = sprintf("%.10g", $(idx["day"]) + 0);',
    '  if ((keep_all || (day in keep)) && $(idx["status"]) == "live") {',
    '    c = $(idx["cell_count"]) + 0;',
    '    n = $(idx["N"]) + 0;',
    '    p = $(idx["ploidy"]) + 0;',
    '    sumc[day] += c;',
    '    sumn[day] += n * c;',
    '    sump[day] += p * c;',
    '    sump2[day] += p * p * c;',
    '  }',
    '}',
    'END {',
    '  for (day in sumc) {',
    '    if (sumc[day] > 0) {',
    '      meanp = sump[day] / sumc[day];',
    '      varp = sump2[day] / sumc[day] - meanp * meanp;',
    '      if (varp < 0 && varp > -1e-9) varp = 0;',
    '      print day, sumn[day] / sumc[day], meanp, sqrt(varp), sumc[day];',
    '    }',
    '  }',
    '}',
    sep = "\n"
  )
  cmd <- paste(
    shQuote(gzip), "-cd", shQuote(path), "|",
    shQuote(awk),
    "-v", shQuote(paste0("days=", days_arg)),
    "-v", shQuote(paste0("all_days=", if (all_days) "1" else "0")),
    shQuote(awk_script)
  )
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
  if (!length(out)) return(data.frame())
  con <- textConnection(out)
  on.exit(close(con), add = TRUE)
  tab <- utils::read.table(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(tab) <- c("day", "simulation_mean_N", "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells")
  tab$day <- as.numeric(tab$day)
  tab <- tab[order(tab$day), , drop = FALSE]
  tab
}


read_simulation_metrics_chunk <- function(tasks, idx, time_points = NULL, progress_every = 100L, worker_label = NULL, total_tasks = nrow(tasks)) {
  rows <- vector("list", length(idx))
  missing <- 0L
  for (j in seq_along(idx)) {
    i <- idx[[j]]
    if (progress_every > 0L && (j == 1L || j %% progress_every == 0L || j == length(idx))) {
      label <- if (is.null(worker_label) || !nzchar(worker_label)) "" else paste0(worker_label, ": ")
      message(label, "Reading simulation state metrics: ", j, "/", length(idx), " (task ", i, "/", total_tasks, ")")
    }
    task <- tasks[i, , drop = FALSE]
    metric <- read_state_metrics_awk(task$state_file[[1]], time_points)
    if (!nrow(metric)) {
      missing <- missing + 1L
      next
    }
    metric$seed_id <- task$seed_id[[1]]
    metric$O2_pct <- task$fixed_o2_pct[[1]]
    metric$O2_key <- num_key(metric$O2_pct)
    metric$initial_ploidy <- task$initial_ploidy[[1]]
    metric$initial_condition <- task$initial_condition[[1]]
    metric$simulation_id <- task$simulation_id[[1]]
    metric$day_key <- num_key(metric$day)
    rows[[j]] <- metric
  }
  if (missing) {
    label <- if (is.null(worker_label) || !nzchar(worker_label)) "" else paste0(worker_label, ": ")
    warning(label, "Missing or unreadable simulation state files: ", missing)
  }
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) out <- data.frame()
  out
}


read_simulation_metrics <- function(tasks, time_points = NULL, progress_every = 100L, n_workers = 1L) {
  n_workers <- suppressWarnings(as.integer(n_workers[[1]]))
  if (!is.finite(n_workers) || is.na(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, nrow(tasks))
  if (n_workers <= 1L || !identical(.Platform$OS.type, "unix")) {
    return(read_simulation_metrics_chunk(
      tasks = tasks,
      idx = seq_len(nrow(tasks)),
      time_points = time_points,
      progress_every = progress_every,
      total_tasks = nrow(tasks)
    ))
  }

  message("Reading simulation state metrics with ", n_workers, " workers: ", nrow(tasks), " tasks")
  idx <- seq_len(nrow(tasks))
  chunks <- split(idx, ((idx - 1L) %% n_workers) + 1L)
  res <- parallel::mclapply(
    seq_along(chunks),
    function(worker_id) {
      read_simulation_metrics_chunk(
        tasks = tasks,
        idx = chunks[[worker_id]],
        time_points = time_points,
        progress_every = progress_every,
        worker_label = sprintf("worker %02d", worker_id),
        total_tasks = nrow(tasks)
      )
    },
    mc.cores = n_workers
  )
  out <- do.call(rbind, res[vapply(res, nrow, integer(1)) > 0L])
  if (is.null(out)) out <- data.frame()
  out
}


filter_simulation_metrics <- function(sim_metrics, time_points, o2_values, initial_ploidy_values, simulation_ids, seed_ids = NULL) {
  if (!nrow(sim_metrics)) return(sim_metrics)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "simulation_id", "simulation_mean_N",
      "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells"),
    names(sim_metrics)
  )
  for (col in numeric_cols) sim_metrics[[col]] <- suppressWarnings(as.numeric(sim_metrics[[col]]))
  out <- filter_by_values(sim_metrics, "day", time_points)
  out <- filter_by_values(out, "O2_pct", o2_values)
  out <- filter_by_values(out, "initial_ploidy", initial_ploidy_values)
  if ("simulation_id" %in% names(out)) out <- out[out$simulation_id %in% simulation_ids, , drop = FALSE]
  if (!is.null(seed_ids) && length(seed_ids) && "seed_id" %in% names(out)) {
    out <- out[out$seed_id %in% seed_ids, , drop = FALSE]
  }
  out$O2_key <- num_key(out$O2_pct)
  out$day_key <- num_key(out$day)
  out
}
