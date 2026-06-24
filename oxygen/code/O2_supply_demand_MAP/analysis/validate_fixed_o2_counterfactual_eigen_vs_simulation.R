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

vf2_seed_labels <- function(seed_ids, mode_arg = NULL) {
  default_modes <- c(seed322 = "mode1", seed25 = "mode2")
  modes <- o2ipa_split_csv(mode_arg, character())
  if (length(modes) == length(seed_ids)) {
    mode_map <- setNames(modes, seed_ids)
  } else {
    mode_map <- default_modes
  }
  labels <- seed_ids
  has_mode <- seed_ids %in% names(mode_map) & nzchar(mode_map[seed_ids])
  labels[has_mode] <- paste0(seed_ids[has_mode], " (", mode_map[seed_ids[has_mode]], ")")
  labels
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
  data.frame(
    mean_N = sum(ngrid * state, na.rm = TRUE),
    mean_ploidy = sum(ngrid * state, na.rm = TRUE) / n_unit,
    fraction_N_le_25 = sum(state[ngrid <= 25], na.rm = TRUE),
    fraction_N_below_44 = sum(state[ngrid < 44], na.rm = TRUE),
    fraction_N_ge_44 = sum(state[ngrid >= 44], na.rm = TRUE),
    fraction_N_ge_66 = sum(state[ngrid >= 66], na.rm = TRUE),
    fraction_N_ge_88 = sum(state[ngrid >= 88], na.rm = TRUE),
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

vf2_compare_states <- function(eigen_states, euler_states, ngrid, n_unit, time_grid) {
  rows <- vector("list", length(time_grid))
  for (i in seq_along(time_grid)) {
    e_state <- eigen_states[[i]]
    u_state <- euler_states[[i]]
    em <- vf2_state_metrics(e_state, ngrid, n_unit)
    um <- vf2_state_metrics(u_state, ngrid, n_unit)
    rows[[i]] <- data.frame(
      day = time_grid[[i]],
      eigen_mean_N = em$mean_N,
      euler_mean_N = um$mean_N,
      diff_mean_N = um$mean_N - em$mean_N,
      abs_diff_mean_N = abs(um$mean_N - em$mean_N),
      eigen_mean_ploidy = em$mean_ploidy,
      euler_mean_ploidy = um$mean_ploidy,
      diff_mean_ploidy = um$mean_ploidy - em$mean_ploidy,
      abs_diff_mean_ploidy = abs(um$mean_ploidy - em$mean_ploidy),
      eigen_fraction_N_le_25 = em$fraction_N_le_25,
      euler_fraction_N_le_25 = um$fraction_N_le_25,
      eigen_fraction_N_below_44 = em$fraction_N_below_44,
      euler_fraction_N_below_44 = um$fraction_N_below_44,
      eigen_fraction_N_ge_44 = em$fraction_N_ge_44,
      euler_fraction_N_ge_44 = um$fraction_N_ge_44,
      eigen_fraction_N_ge_66 = em$fraction_N_ge_66,
      euler_fraction_N_ge_66 = um$fraction_N_ge_66,
      eigen_fraction_N_ge_88 = em$fraction_N_ge_88,
      euler_fraction_N_ge_88 = um$fraction_N_ge_88,
      max_abs_state_diff = max(abs(u_state - e_state), na.rm = TRUE),
      l1_state_diff = sum(abs(u_state - e_state), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

vf2_error_summary_one <- function(comp) {
  terminal <- comp[which.max(comp$day), , drop = FALSE]
  data.frame(
    n_timepoint = nrow(comp),
    max_abs_diff_mean_ploidy = max(comp$abs_diff_mean_ploidy, na.rm = TRUE),
    terminal_abs_diff_mean_ploidy = terminal$abs_diff_mean_ploidy[[1]],
    max_abs_diff_mean_N = max(comp$abs_diff_mean_N, na.rm = TRUE),
    terminal_abs_diff_mean_N = terminal$abs_diff_mean_N[[1]],
    max_abs_state_diff = max(comp$max_abs_state_diff, na.rm = TRUE),
    terminal_max_abs_state_diff = terminal$max_abs_state_diff[[1]],
    max_l1_state_diff = max(comp$l1_state_diff, na.rm = TRUE),
    terminal_l1_state_diff = terminal$l1_state_diff[[1]],
    stringsAsFactors = FALSE
  )
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
    yr <- range(c(sub_init$eigen_mean_ploidy, sub_init$euler_mean_ploidy), na.rm = TRUE)
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
        lines(sub$day, sub$eigen_mean_ploidy, col = eigen_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub <- sub_seed[sub_seed$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub$euler_mean_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.6, lty = 2)
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), "Eigen analytical", paste0("Euler dt=", plot_dt)),
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
    "Fixed-O2 eigen vs numerical simulation",
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
  seed_ids <- vf2_as_seed_vec(args$seed_ids %||% args$seed_id, "seed1")
  seed_labels <- vf2_seed_labels(seed_ids, args$seed_modes %||% args$seed_labels)
  seed_label_map <- stats::setNames(seed_labels, seed_ids)
  o2_grid <- sort(unique(vf2_as_num_vec(args$o2_grid, c(0, 0.5, 1, 2, 5))))
  time_grid <- sort(unique(vf2_as_num_vec(args$time_grid, vf2_default_time_grid())))
  dt_grid <- sort(unique(vf2_as_num_vec(args$dt_grid %||% args$dt, c(0.05))))
  plot_dt <- o2ipa_as_num(args$plot_dt, min(dt_grid))

  vf2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "validate_fixed_o2_counterfactual_eigen_vs_simulation.log")
  sink(log_file, split = TRUE)
  on.exit(sink(), add = TRUE)
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)
  message("seed_ids: ", paste(seed_ids, collapse = ","))
  message("seed_labels: ", paste(seed_labels, collapse = ","))
  message("O2 grid: ", paste(o2_grid, collapse = ","))
  message("dt grid: ", paste(dt_grid, collapse = ","))

  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  init_specs <- data.frame(
    initial_condition = c("init_2N", "init_4N"),
    requested_N = c(44, 88),
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
        for (dt in dt_grid) {
          message("    ", init_specs$initial_condition[[j]], ", dt=", dt)
          eu <- vf2_euler_states(fm$M, init$vector, time_grid, dt)
          comp <- vf2_compare_states(eig$states, eu, fm$ngrid, n_unit, time_grid)
          comp$seed_id <- seed_id
          comp$seed_label <- seed_label_map[[seed_id]]
          comp$O2_pct <- O2
          comp$initial_condition <- init_specs$initial_condition[[j]]
          comp$requested_initial_N <- init_specs$requested_N[[j]]
          comp$used_initial_N <- init$used_N
          comp$dt <- dt
          comp$dominant_growth_rate <- eig$lambda_ref
          comp <- comp[, c(
            "seed_id", "seed_label", "O2_pct", "initial_condition", "requested_initial_N", "used_initial_N",
            "dt", "day", "dominant_growth_rate",
            "eigen_mean_N", "euler_mean_N", "diff_mean_N", "abs_diff_mean_N",
            "eigen_mean_ploidy", "euler_mean_ploidy", "diff_mean_ploidy", "abs_diff_mean_ploidy",
            "eigen_fraction_N_le_25", "euler_fraction_N_le_25",
            "eigen_fraction_N_below_44", "euler_fraction_N_below_44",
            "eigen_fraction_N_ge_44", "euler_fraction_N_ge_44",
            "eigen_fraction_N_ge_66", "euler_fraction_N_ge_66",
            "eigen_fraction_N_ge_88", "euler_fraction_N_ge_88",
            "max_abs_state_diff", "l1_state_diff"
          )]
          tc <- tc + 1L
          trajectory_rows[[tc]] <- comp
          serr <- comp[, c("seed_id", "seed_label", "O2_pct", "initial_condition", "dt", "day", "max_abs_state_diff", "l1_state_diff"), drop = FALSE]
          ec <- ec + 1L
          state_error_rows[[ec]] <- serr
          sm <- vf2_error_summary_one(comp)
          sm$seed_id <- seed_id
          sm$seed_label <- seed_label_map[[seed_id]]
          sm$O2_pct <- O2
          sm$initial_condition <- init_specs$initial_condition[[j]]
          sm$requested_initial_N <- init_specs$requested_N[[j]]
          sm$used_initial_N <- init$used_N
          sm$dt <- dt
          sm$dominant_growth_rate <- eig$lambda_ref
          sm <- sm[, c(
            "seed_id", "seed_label", "O2_pct", "initial_condition", "requested_initial_N", "used_initial_N",
            "dt", "dominant_growth_rate", "n_timepoint",
            "max_abs_diff_mean_ploidy", "terminal_abs_diff_mean_ploidy",
            "max_abs_diff_mean_N", "terminal_abs_diff_mean_N",
            "max_abs_state_diff", "terminal_max_abs_state_diff",
            "max_l1_state_diff", "terminal_l1_state_diff"
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

  run_args <- data.frame(
    argument = c("run_dir", "out_dir", "seed_ids", "seed_labels", "o2_grid", "time_grid", "dt_grid", "plot_dt"),
    value = c(run_dir, out_dir, paste(seed_ids, collapse = ","), paste(seed_labels, collapse = ","), paste(o2_grid, collapse = ","), paste(time_grid, collapse = ","), paste(dt_grid, collapse = ","), as.character(plot_dt)),
    stringsAsFactors = FALSE
  )
  vf2_write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  vf2_write_tsv(trajectories, file.path(out_dir, "tables", "eigen_vs_euler_trajectories.tsv"))
  vf2_write_tsv(state_errors, file.path(out_dir, "tables", "eigen_vs_euler_state_errors.tsv"))
  vf2_write_tsv(summaries, file.path(out_dir, "tables", "eigen_vs_euler_error_summary.tsv"))
  vf2_plot(
    trajectories,
    file.path(out_dir, "figures", sprintf("eigen_vs_euler_mean_ploidy_%s.pdf", seed_slug)),
    plot_dt = plot_dt
  )
  message("Completed validation: ", out_dir)
}

if (identical(environment(), globalenv())) {
  main()
}
