#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    frame_files <- Filter(
      nzchar,
      vapply(sys.frames(), function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      }, character(1))
    )
    if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
  }
})

source(file.path(SCRIPT_DIR, "process_fingerprint_utils.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "ploidy_regime_utils.R"), local = TRUE)

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required")
})

options(error = function() {
  traceback(2)
  quit(status = 1)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

vf2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

vf2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

vf2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

vf2_as_seed_vec <- function(x, default) {
  vals <- o2ipa_split_csv(x, default)
  vals <- unique(o2ipa_norm_seed(vals))
  if (length(vals)) vals else default
}

vf2_seed_modes <- function(seed_ids, mode_arg = NULL) {
  default_modes <- c(seed322 = "mode1", seed25 = "mode2")
  modes <- o2ipa_split_csv(mode_arg, character())
  if (length(modes) == length(seed_ids)) {
    mode_map <- setNames(modes, seed_ids)
  } else {
    mode_map <- default_modes
  }
  out <- rep("", length(seed_ids))
  names(out) <- seed_ids
  has_mode <- seed_ids %in% names(mode_map) & nzchar(mode_map[seed_ids])
  out[has_mode] <- unname(mode_map[seed_ids[has_mode]])
  out
}

vf2_seed_labels <- function(seed_ids, mode_arg = NULL) {
  seed_modes <- vf2_seed_modes(seed_ids, mode_arg)
  labels <- seed_ids
  has_mode <- nzchar(seed_modes)
  labels[has_mode] <- paste0(seed_ids[has_mode], " (", seed_modes[has_mode], ")")
  labels
}

vf2_as_int_vec <- function(x, default) {
  vals <- suppressWarnings(as.integer(o2ipa_split_csv(x, as.character(default))))
  vals <- vals[!is.na(vals) & vals > 0L]
  if (length(vals)) unique(vals) else default
}

vf2_num_path_tag <- function(x) {
  val <- suppressWarnings(as.numeric(x))
  if (length(val) != 1L || !is.finite(val)) return("NA")
  txt <- format(val, scientific = FALSE, trim = TRUE)
  txt <- sub("^-", "m", txt)
  txt <- gsub("\\.", "p", txt)
  txt <- gsub("[^A-Za-z0-9]+", "", txt)
  if (!nzchar(txt)) "NA" else txt
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

vf2_plot <- function(traj, out_path, plot_dt) {
  d <- traj[abs(traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  seed_use <- unique(d$seed_id)
  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  eigen_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(d$seed_label[d$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init <- d[d$initial_condition == init, , drop = FALSE]
    yr <- range(c(sub_init$expm_mean_ploidy, sub_init$euler_mean_ploidy), na.rm = TRUE)
    xr <- range(sub_init$day, na.rm = TRUE)
    for (seed in seed_use) {
      sub_seed <- sub_init[sub_init$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Day",
        ylab = "Mean ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed[sub_seed$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub$expm_mean_ploidy, col = eigen_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub <- sub_seed[sub_seed$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub$euler_mean_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.6, lty = 2)
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), "Expm analytical", paste0("Euler dt=", plot_dt)),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    "Fixed-O2 expm analytical vs Euler integration",
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}

vf2_plot_solution_vs_simulation <- function(solution_traj, simulation_traj, out_path, plot_dt, solution_name, solution_label, title) {
  d <- solution_traj[abs(solution_traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d) || is.null(simulation_traj) || !nrow(simulation_traj)) return(invisible(FALSE))
  solution_col <- paste0(solution_name, "_mean_ploidy")
  if (!solution_col %in% names(d)) stop("solution_traj missing column: ", solution_col)
  seed_use <- unique(d$seed_id)
  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  analytical_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(d$seed_label[d$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init_sol <- d[d$initial_condition == init, , drop = FALSE]
    sub_init_sim <- simulation_traj[simulation_traj$initial_condition == init, , drop = FALSE]
    yr <- range(c(sub_init_sol[[solution_col]], sub_init_sim$simulation_mean_ploidy), na.rm = TRUE)
    xr <- range(c(sub_init_sol$day, sub_init_sim$day), na.rm = TRUE)
    for (seed in seed_use) {
      sub_seed_sol <- sub_init_sol[sub_init_sol$seed_id == seed, , drop = FALSE]
      sub_seed_sim <- sub_init_sim[sub_init_sim$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Day",
        ylab = "Mean ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed_sol[sub_seed_sol$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub[[solution_col]], col = analytical_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub_o2_sim <- sub_seed_sim[sub_seed_sim$O2_pct == O2, , drop = FALSE]
        reps <- split(sub_o2_sim, sub_o2_sim$simulation_id)
        for (rep_dat in reps) {
          rep_dat <- rep_dat[order(rep_dat$day), , drop = FALSE]
          lines(rep_dat$day, rep_dat$simulation_mean_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.2, lty = 2)
        }
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), solution_label, "Simulation reps"),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    title,
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}

vf2_range_pad <- function(x, frac = 0.04) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  xr <- range(x)
  span <- diff(xr)
  if (!is.finite(span) || span <= 0) {
    pad <- max(abs(xr[[1]]) * frac, 0.05)
  } else {
    pad <- span * frac
  }
  xr + c(-pad, pad)
}

vf2_solution_phase_plane <- function(solution_traj, plot_dt, solution_name, seed_mode_map) {
  d <- solution_traj[abs(solution_traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  mean_col <- paste0(solution_name, "_mean_ploidy")
  sd_col <- paste0(solution_name, "_sd_ploidy")
  if (!all(c(mean_col, sd_col) %in% names(d))) stop("solution_traj missing phase-plane columns for method: ", solution_name)
  out <- data.frame(
    seed_id = d$seed_id,
    seed_mode = unname(seed_mode_map[d$seed_id]),
    seed_label = d$seed_label,
    initial_ploidy = d$initial_ploidy,
    initial_condition = d$initial_condition,
    o2_pct = d$O2_pct,
    simulation_id = NA_integer_,
    solution_method = solution_name,
    time_days = d$day,
    mean_ploidy = d[[mean_col]],
    sd_ploidy = d[[sd_col]],
    dt = d$dt,
    stringsAsFactors = FALSE
  )
  seed_rank <- match(out$seed_id, names(seed_mode_map))
  out[order(seed_rank, out$initial_ploidy, out$o2_pct, out$time_days), , drop = FALSE]
}

vf2_simulation_phase_plane <- function(simulation_traj, seed_mode_map) {
  if (is.null(simulation_traj) || !nrow(simulation_traj)) return(data.frame())
  out <- data.frame(
    seed_id = simulation_traj$seed_id,
    seed_mode = if ("seed_mode" %in% names(simulation_traj)) simulation_traj$seed_mode else unname(seed_mode_map[simulation_traj$seed_id]),
    seed_label = simulation_traj$seed_label,
    initial_ploidy = simulation_traj$initial_ploidy,
    initial_condition = simulation_traj$initial_condition,
    o2_pct = simulation_traj$O2_pct,
    simulation_id = simulation_traj$simulation_id,
    solution_method = "simulation",
    time_days = simulation_traj$day,
    mean_ploidy = simulation_traj$simulation_mean_ploidy,
    sd_ploidy = simulation_traj$simulation_sd_ploidy,
    dt = NA_real_,
    stringsAsFactors = FALSE
  )
  seed_rank <- match(out$seed_id, names(seed_mode_map))
  out[order(seed_rank, out$initial_ploidy, out$o2_pct, out$simulation_id, out$time_days), , drop = FALSE]
}

vf2_plot_phase_plane_solution_vs_simulation <- function(solution_phase, simulation_phase, out_path, solution_label, title) {
  if (is.null(solution_phase) || !nrow(solution_phase) || is.null(simulation_phase) || !nrow(simulation_phase)) {
    return(invisible(FALSE))
  }
  seed_use <- unique(solution_phase$seed_id)
  o2_use <- sort(unique(solution_phase$o2_pct))
  init_use <- unique(solution_phase$initial_condition)
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  analytical_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(solution_phase$seed_label[solution_phase$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4.4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init_sol <- solution_phase[solution_phase$initial_condition == init, , drop = FALSE]
    sub_init_sim <- simulation_phase[simulation_phase$initial_condition == init, , drop = FALSE]
    xr <- vf2_range_pad(c(sub_init_sol$mean_ploidy, sub_init_sim$mean_ploidy))
    yr <- vf2_range_pad(c(sub_init_sol$sd_ploidy, sub_init_sim$sd_ploidy))
    for (seed in seed_use) {
      sub_seed_sol <- sub_init_sol[sub_init_sol$seed_id == seed, , drop = FALSE]
      sub_seed_sim <- sub_init_sim[sub_init_sim$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Mean ploidy",
        ylab = "SD ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed_sol[sub_seed_sol$o2_pct == O2, , drop = FALSE]
        sub <- sub[order(sub$time_days), , drop = FALSE]
        lines(sub$mean_ploidy, sub$sd_ploidy, col = analytical_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub_o2_sim <- sub_seed_sim[sub_seed_sim$o2_pct == O2, , drop = FALSE]
        reps <- split(sub_o2_sim, sub_o2_sim$simulation_id)
        for (rep_dat in reps) {
          rep_dat <- rep_dat[order(rep_dat$time_days), , drop = FALSE]
          lines(rep_dat$mean_ploidy, rep_dat$sd_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.1, lty = 2)
        }
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), solution_label, "Simulation reps"),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    title,
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}

main <- function() {
  args <- o2ipa_parse_args()
  repo_root <- o2ipa_repo_root(SCRIPT_DIR)
  default_hpc_run_dir <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed"
  local_run_dir <- file.path(repo_root, "oxygen", "results", "fit_invivo_O2_buffering_500seed")
  default_run_dir <- if (dir.exists(default_hpc_run_dir)) default_hpc_run_dir else local_run_dir
  run_dir <- normalizePath(o2ipa_as_chr(args$run_dir, default_run_dir), mustWork = FALSE)
  out_dir <- normalizePath(o2ipa_as_chr(args$out_dir, file.path(repo_root, "oxygen", "results", "analysis", "fixed_o2_counterfactual_eigen_validation_seed1")), mustWork = FALSE)
  simulation_dir <- normalizePath(o2ipa_as_chr(args$simulation_dir, file.path(repo_root, "oxygen", "results", "O2_fixed_simulation")), mustWork = FALSE)
  simulation_mode <- o2ipa_as_chr(args$simulation_mode, "invivo")
  simulation_ids <- vf2_as_int_vec(args$simulation_ids, seq_len(10L))
  include_simulation <- o2ipa_as_bool(args$include_simulation, dir.exists(simulation_dir))
  seed_ids <- vf2_as_seed_vec(args$seed_ids %||% args$seed_id, "seed1")
  seed_mode_arg <- args$seed_modes %||% args$seed_labels
  seed_modes <- vf2_seed_modes(seed_ids, seed_mode_arg)
  seed_labels <- vf2_seed_labels(seed_ids, seed_mode_arg)
  seed_mode_map <- stats::setNames(seed_modes, seed_ids)
  seed_label_map <- stats::setNames(seed_labels, seed_ids)
  o2_grid <- sort(unique(vf2_as_num_vec(args$o2_grid, c(0, 0.5, 1, 2, 5))))
  default_time_grid <- if (include_simulation) {
    seq(0, o2ipa_as_num(args$time_days, 1000), by = o2ipa_as_num(args$time_step_days, 1))
  } else {
    vf2_default_time_grid()
  }
  time_grid <- sort(unique(vf2_as_num_vec(args$time_grid, default_time_grid)))
  dt_grid <- sort(unique(vf2_as_num_vec(args$dt_grid %||% args$dt, c(0.05))))
  plot_dt <- o2ipa_as_num(args$plot_dt, min(dt_grid))

  vf2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "validate_fixed_o2_counterfactual_eigen_vs_simulation.log")
  sink(log_file, split = TRUE)
  on.exit(sink(), add = TRUE)
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)
  message("simulation_dir: ", simulation_dir)
  message("include_simulation: ", include_simulation)
  if (include_simulation) message("simulation_ids: ", paste(simulation_ids, collapse = ","))
  message("seed_ids: ", paste(seed_ids, collapse = ","))
  message("seed_modes: ", paste(seed_modes, collapse = ","))
  message("seed_labels: ", paste(seed_labels, collapse = ","))
  message("O2 grid: ", paste(o2_grid, collapse = ","))
  message("dt grid: ", paste(dt_grid, collapse = ","))

  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  init_specs <- data.frame(
    initial_condition = c("init_2N", "init_4N"),
    requested_N = c(44, 88),
    initial_ploidy = c(2, 4),
    stringsAsFactors = FALSE
  )

  trajectory_rows <- list()
  state_error_rows <- list()
  summary_rows <- list()
  tc <- 0L
  ec <- 0L
  sc <- 0L
  for (seed_id in seed_ids) {
    message("Processing seed=", seed_id, " label=", seed_label_map[[seed_id]])
    manifest <- inputs$manifest[inputs$manifest$seed_id == seed_id, , drop = FALSE]
    if (!nrow(manifest)) stop("seed_id not found in run_dir: ", seed_id)
    params_long <- inputs$params_long[inputs$params_long$seed_id == seed_id, , drop = FALSE]
    if (!nrow(params_long)) stop("No parameters found for ", seed_id)
    cfg <- o2pr_first_seed_cfg(manifest)
    n_unit <- as.numeric(cfg$N_UNIT %||% 22)
    if (!seed_id %in% rownames(param_mat)) stop("seed_id missing from parameter matrix: ", seed_id)
    pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    for (O2 in o2_grid) {
      message("  O2=", O2)
      fm <- vf2_fixed_matrix(model_env, cfg, run_params, O2)
      for (j in seq_len(nrow(init_specs))) {
        init <- vf2_init_vector(fm$ngrid, init_specs$requested_N[[j]])
        eig <- vf2_eigen_states(fm$M, init$vector, time_grid)
        expm <- vf2_expm_states(fm$M, init$vector, time_grid)
        for (dt in dt_grid) {
          message("    ", init_specs$initial_condition[[j]], ", dt=", dt)
          eu <- vf2_euler_states(fm$M, init$vector, time_grid, dt)
          comp <- vf2_compare_states(eig$states, expm, eu, fm$ngrid, n_unit, time_grid)
          comp$seed_id <- seed_id
          comp$seed_mode <- seed_mode_map[[seed_id]]
          comp$seed_label <- seed_label_map[[seed_id]]
          comp$O2_pct <- O2
          comp$initial_condition <- init_specs$initial_condition[[j]]
          comp$requested_initial_N <- init_specs$requested_N[[j]]
          comp$initial_ploidy <- init_specs$initial_ploidy[[j]]
          comp$used_initial_N <- init$used_N
          comp$dt <- dt
          comp$dominant_growth_rate <- eig$lambda_ref
          comp <- comp[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy", "requested_initial_N", "used_initial_N",
            "dt", "day", "dominant_growth_rate",
            "eigen_mean_N", "expm_mean_N", "euler_mean_N",
            "diff_euler_expm_mean_N", "abs_diff_euler_expm_mean_N",
            "diff_eigen_expm_mean_N", "abs_diff_eigen_expm_mean_N",
            "eigen_mean_ploidy", "expm_mean_ploidy", "euler_mean_ploidy",
            "diff_euler_expm_mean_ploidy", "abs_diff_euler_expm_mean_ploidy",
            "diff_eigen_expm_mean_ploidy", "abs_diff_eigen_expm_mean_ploidy",
            "eigen_sd_ploidy", "expm_sd_ploidy", "euler_sd_ploidy",
            "diff_euler_expm_sd_ploidy", "abs_diff_euler_expm_sd_ploidy",
            "diff_eigen_expm_sd_ploidy", "abs_diff_eigen_expm_sd_ploidy",
            "eigen_fraction_N_le_25", "expm_fraction_N_le_25", "euler_fraction_N_le_25",
            "eigen_fraction_N_below_44", "expm_fraction_N_below_44", "euler_fraction_N_below_44",
            "eigen_fraction_N_ge_44", "expm_fraction_N_ge_44", "euler_fraction_N_ge_44",
            "eigen_fraction_N_ge_66", "expm_fraction_N_ge_66", "euler_fraction_N_ge_66",
            "eigen_fraction_N_ge_88", "expm_fraction_N_ge_88", "euler_fraction_N_ge_88",
            "max_abs_state_diff_euler_expm", "l1_state_diff_euler_expm",
            "max_abs_state_diff_eigen_expm", "l1_state_diff_eigen_expm"
          )]
          tc <- tc + 1L
          trajectory_rows[[tc]] <- comp
          serr <- comp[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "dt", "day",
            "max_abs_state_diff_euler_expm", "l1_state_diff_euler_expm",
            "max_abs_state_diff_eigen_expm", "l1_state_diff_eigen_expm"
          ), drop = FALSE]
          ec <- ec + 1L
          state_error_rows[[ec]] <- serr
          sm <- vf2_error_summary_one(comp)
          sm$seed_id <- seed_id
          sm$seed_mode <- seed_mode_map[[seed_id]]
          sm$seed_label <- seed_label_map[[seed_id]]
          sm$O2_pct <- O2
          sm$initial_condition <- init_specs$initial_condition[[j]]
          sm$initial_ploidy <- init_specs$initial_ploidy[[j]]
          sm$requested_initial_N <- init_specs$requested_N[[j]]
          sm$used_initial_N <- init$used_N
          sm$dt <- dt
          sm$dominant_growth_rate <- eig$lambda_ref
          sm <- sm[, c(
            "seed_id", "seed_mode", "seed_label", "O2_pct", "initial_condition", "initial_ploidy", "requested_initial_N", "used_initial_N",
            "dt", "dominant_growth_rate", "n_timepoint",
            "max_abs_diff_euler_expm_mean_ploidy", "terminal_abs_diff_euler_expm_mean_ploidy",
            "max_abs_diff_eigen_expm_mean_ploidy", "terminal_abs_diff_eigen_expm_mean_ploidy",
            "max_abs_diff_euler_expm_sd_ploidy", "terminal_abs_diff_euler_expm_sd_ploidy",
            "max_abs_diff_eigen_expm_sd_ploidy", "terminal_abs_diff_eigen_expm_sd_ploidy",
            "max_abs_diff_euler_expm_mean_N", "terminal_abs_diff_euler_expm_mean_N",
            "max_abs_state_diff_euler_expm", "terminal_max_abs_state_diff_euler_expm",
            "max_l1_state_diff_euler_expm", "terminal_l1_state_diff_euler_expm"
          )]
          sc <- sc + 1L
          summary_rows[[sc]] <- sm
        }
      }
    }
  }

  trajectories <- do.call(rbind, trajectory_rows)
  state_errors <- do.call(rbind, state_error_rows)
  summaries <- do.call(rbind, summary_rows)
  summaries <- summaries[order(match(summaries$seed_id, seed_ids), summaries$O2_pct, summaries$initial_condition, summaries$dt), , drop = FALSE]
  seed_slug <- paste(seed_ids, collapse = "_")

  simulation_traj <- data.frame()
  simulation_summary <- data.frame()
  expm_simulation_comparison <- data.frame()
  expm_simulation_comparison_summary <- data.frame()
  eigen_simulation_comparison <- data.frame()
  eigen_simulation_comparison_summary <- data.frame()
  expm_phase_plane <- data.frame()
  eigen_phase_plane <- data.frame()
  simulation_phase_plane <- data.frame()
  if (include_simulation) {
    message("Reading fixed-O2 simulation trajectories")
    simulation_traj <- vf2_read_simulation_trajectories(
      sim_root = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      seed_label_map = seed_label_map,
      seed_mode_map = seed_mode_map,
      o2_grid = o2_grid,
      init_specs = init_specs,
      simulation_ids = simulation_ids
    )
    simulation_summary <- vf2_simulation_summary(simulation_traj)
    expm_simulation_comparison <- vf2_compare_solution_to_simulation(trajectories, simulation_summary, "expm")
    expm_simulation_comparison_summary <- vf2_simulation_comparison_summary(expm_simulation_comparison)
    eigen_simulation_comparison <- vf2_compare_solution_to_simulation(trajectories, simulation_summary, "eigen")
    eigen_simulation_comparison_summary <- vf2_simulation_comparison_summary(eigen_simulation_comparison)
  }
  expm_phase_plane <- vf2_solution_phase_plane(trajectories, plot_dt, "expm", seed_mode_map)
  eigen_phase_plane <- vf2_solution_phase_plane(trajectories, plot_dt, "eigen", seed_mode_map)
  if (include_simulation) {
    simulation_phase_plane <- vf2_simulation_phase_plane(simulation_traj, seed_mode_map)
  }

  run_args <- data.frame(
    argument = c(
      "run_dir", "out_dir", "simulation_dir", "simulation_mode", "include_simulation",
      "simulation_ids", "seed_ids", "seed_modes", "seed_labels", "o2_grid", "time_grid", "dt_grid", "plot_dt"
    ),
    value = c(
      run_dir, out_dir, simulation_dir, simulation_mode, as.character(include_simulation),
      paste(simulation_ids, collapse = ","), paste(seed_ids, collapse = ","), paste(seed_modes, collapse = ","), paste(seed_labels, collapse = ","),
      paste(o2_grid, collapse = ","), paste(time_grid, collapse = ","), paste(dt_grid, collapse = ","), as.character(plot_dt)
    ),
    stringsAsFactors = FALSE
  )
  vf2_write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  vf2_write_tsv(trajectories, file.path(out_dir, "tables", "eigen_vs_euler_trajectories.tsv"))
  vf2_write_tsv(state_errors, file.path(out_dir, "tables", "eigen_vs_euler_state_errors.tsv"))
  vf2_write_tsv(summaries, file.path(out_dir, "tables", "eigen_vs_euler_error_summary.tsv"))
  vf2_write_tsv(expm_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_expm_analytical.tsv"))
  vf2_write_tsv(eigen_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_eigen_analytical.tsv"))
  if (include_simulation) {
    vf2_write_tsv(simulation_traj, file.path(out_dir, "tables", "simulation_replicate_mean_ploidy_trajectories.tsv"))
    vf2_write_tsv(simulation_summary, file.path(out_dir, "tables", "simulation_summary_mean_ploidy_trajectories.tsv"))
    vf2_write_tsv(expm_simulation_comparison, file.path(out_dir, "tables", "expm_vs_simulation_mean_ploidy_comparison.tsv"))
    vf2_write_tsv(expm_simulation_comparison_summary, file.path(out_dir, "tables", "expm_vs_simulation_error_summary.tsv"))
    vf2_write_tsv(eigen_simulation_comparison, file.path(out_dir, "tables", "eigen_vs_simulation_mean_ploidy_comparison.tsv"))
    vf2_write_tsv(eigen_simulation_comparison_summary, file.path(out_dir, "tables", "eigen_vs_simulation_error_summary.tsv"))
    vf2_write_tsv(simulation_phase_plane, file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_simulation_replicates.tsv"))
    vf2_write_tsv(rbind(expm_phase_plane, simulation_phase_plane), file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_expm_vs_simulation.tsv"))
    vf2_write_tsv(rbind(eigen_phase_plane, simulation_phase_plane), file.path(out_dir, "tables", "phase_plane_mean_sd_ploidy_eigen_vs_simulation.tsv"))
  }
  vf2_plot(
    trajectories,
    file.path(out_dir, "figures", sprintf("expm_vs_euler_mean_ploidy_%s.pdf", seed_slug)),
    plot_dt = plot_dt
  )
  if (include_simulation) {
    vf2_plot_solution_vs_simulation(
      trajectories,
      simulation_traj,
      file.path(out_dir, "figures", sprintf("expm_vs_simulation_replicates_mean_ploidy_%s.pdf", seed_slug)),
      plot_dt = plot_dt,
      solution_name = "expm",
      solution_label = "Expm analytical",
      title = "Fixed-O2 expm analytical vs simulation replicates"
    )
    vf2_plot_solution_vs_simulation(
      trajectories,
      simulation_traj,
      file.path(out_dir, "figures", sprintf("eigen_vs_simulation_replicates_mean_ploidy_%s.pdf", seed_slug)),
      plot_dt = plot_dt,
      solution_name = "eigen",
      solution_label = "Eigen analytical",
      title = "Fixed-O2 eigen analytical vs simulation replicates"
    )
    vf2_plot_phase_plane_solution_vs_simulation(
      expm_phase_plane,
      simulation_phase_plane,
      file.path(out_dir, "figures", sprintf("expm_vs_simulation_phase_plane_mean_sd_ploidy_%s.pdf", seed_slug)),
      solution_label = "Expm analytical",
      title = "Fixed-O2 phase plane: expm analytical vs simulation replicates"
    )
    vf2_plot_phase_plane_solution_vs_simulation(
      eigen_phase_plane,
      simulation_phase_plane,
      file.path(out_dir, "figures", sprintf("eigen_vs_simulation_phase_plane_mean_sd_ploidy_%s.pdf", seed_slug)),
      solution_label = "Eigen analytical",
      title = "Fixed-O2 phase plane: eigen analytical vs simulation replicates"
    )
  }
  message("Completed validation: ", out_dir)
}

if (identical(environment(), globalenv())) {
  main()
}
