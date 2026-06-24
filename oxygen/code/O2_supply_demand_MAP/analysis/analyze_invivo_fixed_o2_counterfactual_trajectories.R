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

cf2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

cf2_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

cf2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs", "report")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

cf2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

cf2_default_time_grid <- function() {
  sort(unique(c(seq(0, 100, by = 1), 125, 150, 175, 200, 250, 300, 400, 500, 700, 1000)))
}

cf2_load_labels <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(data.frame(seed_id = character()))
  tab <- cf2_read_tsv(path)
  if (!"seed_id" %in% names(tab)) stop("label_file must contain seed_id: ", path)
  keep <- intersect(c(
    "seed_id", "trajectory_regime", "mode_label", "objective", "delta_objective",
    "mean_terminal_ploidy", "mean_late_drop_amplitude",
    "2N__terminal_median_ploidy", "4N__terminal_median_ploidy",
    "2N__late_drop_amplitude", "4N__late_drop_amplitude"
  ), names(tab))
  tab[, keep, drop = FALSE]
}

cf2_fixed_matrix <- function(model_env, cfg, run_params, O2) {
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22L), as.integer(cfg$N_MAX %||% 154L))
  G <- o2pr_build_G(model_env, cfg, run_params, O2)
  mu_all <- as.numeric(o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(ngrid)), run_params = run_params, N = ngrid))
  M <- as.matrix(G - Matrix::Diagonal(x = mu_all))
  list(M = M, G = G, mu_all = mu_all, ngrid = ngrid)
}

cf2_init_vector <- function(ngrid, init_N) {
  idx <- which.min(abs(ngrid - init_N))
  v <- numeric(length(ngrid))
  v[[idx]] <- 1
  list(vector = v, used_N = ngrid[[idx]], used_ploidy = ngrid[[idx]] / 22)
}

cf2_normalize_state <- function(x) {
  x <- Re(x)
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  x <- pmax(x, 0)
  s <- sum(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  x / s
}

cf2_eigen_trajectory <- function(M, ngrid, init, time_grid, n_unit) {
  eig <- tryCatch(eigen(M, only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) return(list(status = "eigen_failed", trajectory = data.frame()))
  coef <- tryCatch(solve(eig$vectors, init), error = function(e) NULL)
  if (is.null(coef)) return(list(status = "eigen_solve_failed", trajectory = data.frame()))
  lambda_ref <- max(Re(eig$values), na.rm = TRUE)
  rows <- lapply(time_grid, function(tt) {
    weights <- exp((eig$values - lambda_ref) * tt) * coef
    state <- cf2_normalize_state(eig$vectors %*% weights)
    data.frame(
      day = tt,
      mean_N = sum(ngrid * state, na.rm = TRUE),
      mean_ploidy = sum(ngrid * state, na.rm = TRUE) / n_unit,
      fraction_N_le_25 = sum(state[ngrid <= 25], na.rm = TRUE),
      fraction_N_below_44 = sum(state[ngrid < 44], na.rm = TRUE),
      fraction_N_ge_44 = sum(state[ngrid >= 44], na.rm = TRUE),
      fraction_N_ge_66 = sum(state[ngrid >= 66], na.rm = TRUE),
      fraction_N_ge_88 = sum(state[ngrid >= 88], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (any(!is.finite(out$mean_ploidy))) {
    return(list(status = "nonfinite_state", trajectory = out))
  }
  list(status = "ok", trajectory = out)
}

cf2_crossing_time <- function(day, value, threshold, direction = "down") {
  ok <- is.finite(day) & is.finite(value)
  day <- day[ok]
  value <- value[ok]
  if (length(day) < 2L) return(NA_real_)
  ord <- order(day)
  day <- day[ord]
  value <- value[ord]
  hit <- if (identical(direction, "down")) which(value <= threshold) else which(value >= threshold)
  if (!length(hit)) return(NA_real_)
  i <- hit[[1]]
  if (i == 1L) return(day[[1]])
  x0 <- day[[i - 1L]]
  x1 <- day[[i]]
  y0 <- value[[i - 1L]]
  y1 <- value[[i]]
  if (!is.finite(y1 - y0) || abs(y1 - y0) < 1e-12) return(x1)
  x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)
}

cf2_summarize_trajectory <- function(traj, max_time) {
  if (!nrow(traj)) {
    return(data.frame(
      terminal_mean_ploidy = NA_real_,
      terminal_fraction_N_le_25 = NA_real_,
      terminal_fraction_N_below_44 = NA_real_,
      min_mean_ploidy = NA_real_,
      day_min_mean_ploidy = NA_real_,
      time_crossing_ploidy_2p0_down = NA_real_,
      time_crossing_ploidy_1p8_down = NA_real_,
      time_crossing_ploidy_1p5_down = NA_real_,
      time_crossing_ploidy_1p3_down = NA_real_,
      time_crossing_ploidy_1p1_down = NA_real_,
      time_crossing_ploidy_1p5_down_censored = max_time + 1,
      reaches_ploidy_1p5 = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  last <- traj[which.max(traj$day), , drop = FALSE]
  i_min <- which.min(traj$mean_ploidy)
  t15 <- cf2_crossing_time(traj$day, traj$mean_ploidy, 1.5, "down")
  data.frame(
    terminal_mean_ploidy = last$mean_ploidy[[1]],
    terminal_fraction_N_le_25 = last$fraction_N_le_25[[1]],
    terminal_fraction_N_below_44 = last$fraction_N_below_44[[1]],
    min_mean_ploidy = traj$mean_ploidy[[i_min]],
    day_min_mean_ploidy = traj$day[[i_min]],
    time_crossing_ploidy_2p0_down = cf2_crossing_time(traj$day, traj$mean_ploidy, 2.0, "down"),
    time_crossing_ploidy_1p8_down = cf2_crossing_time(traj$day, traj$mean_ploidy, 1.8, "down"),
    time_crossing_ploidy_1p5_down = t15,
    time_crossing_ploidy_1p3_down = cf2_crossing_time(traj$day, traj$mean_ploidy, 1.3, "down"),
    time_crossing_ploidy_1p1_down = cf2_crossing_time(traj$day, traj$mean_ploidy, 1.1, "down"),
    time_crossing_ploidy_1p5_down_censored = if (is.finite(t15)) t15 else max_time + 1,
    reaches_ploidy_1p5 = is.finite(t15),
    stringsAsFactors = FALSE
  )
}

cf2_dominant_one <- function(M, ngrid, n_unit) {
  eig <- tryCatch(eigen(M, only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) {
    return(data.frame(dominant_mean_ploidy = NA_real_, dominant_fraction_N_le_25 = NA_real_, dominant_growth_rate = NA_real_, spectral_gap = NA_real_))
  }
  idx <- which.max(Re(eig$values))
  v <- Re(eig$vectors[, idx])
  if (sum(v, na.rm = TRUE) < 0) v <- -v
  v <- cf2_normalize_state(v)
  lambda1 <- Re(eig$values[[idx]])
  lambda2 <- sort(Re(eig$values), decreasing = TRUE)[min(2L, length(eig$values))]
  data.frame(
    dominant_mean_ploidy = sum(ngrid * v, na.rm = TRUE) / n_unit,
    dominant_fraction_N_le_25 = sum(v[ngrid <= 25], na.rm = TRUE),
    dominant_growth_rate = lambda1,
    spectral_gap = lambda1 - lambda2,
    stringsAsFactors = FALSE
  )
}

cf2_wilcox_row <- function(tab, value_col, group_col = "trajectory_regime",
                           a = "mode1_ploidy_stable", b = "mode2_second_ploidy_collapse",
                           feature_name = value_col) {
  d <- tab[tab[[group_col]] %in% c(a, b) & is.finite(tab[[value_col]]), , drop = FALSE]
  if (!nrow(d) || length(unique(d[[group_col]])) < 2L) {
    return(data.frame(feature = feature_name, n_mode1 = 0L, n_mode2 = 0L, median_mode1 = NA_real_, mean_mode1 = NA_real_, median_mode2 = NA_real_, mean_mode2 = NA_real_, median_diff_mode2_minus_mode1 = NA_real_, mean_diff_mode2_minus_mode1 = NA_real_, wilcox_p_value = NA_real_))
  }
  x <- d[[value_col]][d[[group_col]] == a]
  y <- d[[value_col]][d[[group_col]] == b]
  p <- if (length(x) >= 1L && length(y) >= 1L) suppressWarnings(stats::wilcox.test(x, y)$p.value) else NA_real_
  data.frame(
    feature = feature_name,
    n_mode1 = length(x),
    n_mode2 = length(y),
    median_mode1 = stats::median(x, na.rm = TRUE),
    mean_mode1 = mean(x, na.rm = TRUE),
    median_mode2 = stats::median(y, na.rm = TRUE),
    mean_mode2 = mean(y, na.rm = TRUE),
    median_diff_mode2_minus_mode1 = stats::median(y, na.rm = TRUE) - stats::median(x, na.rm = TRUE),
    mean_diff_mode2_minus_mode1 = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
    wilcox_p_value = p,
    stringsAsFactors = FALSE
  )
}

cf2_regime_summary <- function(summary_by_seed) {
  rows <- list()
  counter <- 0L
  regs <- c("mode1_ploidy_stable", "mode2_second_ploidy_collapse", "ambiguous")
  for (O2 in sort(unique(summary_by_seed$O2_pct))) {
    for (init in unique(summary_by_seed$initial_condition)) {
      for (reg in regs) {
        d <- summary_by_seed[summary_by_seed$O2_pct == O2 & summary_by_seed$initial_condition == init & summary_by_seed$trajectory_regime == reg, , drop = FALSE]
        if (!nrow(d)) next
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          O2_pct = O2,
          initial_condition = init,
          trajectory_regime = reg,
          n_seed = nrow(d),
          terminal_mean_ploidy_median = stats::median(d$terminal_mean_ploidy, na.rm = TRUE),
          terminal_mean_ploidy_mean = mean(d$terminal_mean_ploidy, na.rm = TRUE),
          terminal_mean_ploidy_q25 = as.numeric(stats::quantile(d$terminal_mean_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
          terminal_mean_ploidy_q75 = as.numeric(stats::quantile(d$terminal_mean_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
          reaches_ploidy_1p5_fraction = mean(d$reaches_ploidy_1p5 %in% TRUE, na.rm = TRUE),
          time_crossing_1p5_censored_median = stats::median(d$time_crossing_ploidy_1p5_down_censored, na.rm = TRUE),
          time_crossing_1p5_censored_mean = mean(d$time_crossing_ploidy_1p5_down_censored, na.rm = TRUE),
          terminal_minus_dominant_median = stats::median(d$terminal_minus_dominant_ploidy, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

cf2_regime_tests <- function(summary_by_seed) {
  metrics <- c(
    "terminal_mean_ploidy", "terminal_fraction_N_le_25", "terminal_fraction_N_below_44",
    "min_mean_ploidy", "time_crossing_ploidy_1p5_down_censored",
    "time_crossing_ploidy_1p3_down", "terminal_minus_dominant_ploidy"
  )
  rows <- list()
  counter <- 0L
  for (O2 in sort(unique(summary_by_seed$O2_pct))) {
    for (init in unique(summary_by_seed$initial_condition)) {
      d <- summary_by_seed[summary_by_seed$O2_pct == O2 & summary_by_seed$initial_condition == init, , drop = FALSE]
      for (metric in metrics) {
        counter <- counter + 1L
        r <- cf2_wilcox_row(d, metric, feature_name = paste(metric, "O2", O2, init, sep = "__"))
        r$O2_pct <- O2
        r$initial_condition <- init
        r$metric <- metric
        rows[[counter]] <- r
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) out$BH_adjusted_p_value <- stats::p.adjust(out$wilcox_p_value, method = "BH")
  out[order(out$BH_adjusted_p_value, out$O2_pct, out$initial_condition, out$metric), , drop = FALSE]
}

cf2_parameter_correlations <- function(summary_by_seed, params_long) {
  if (!nrow(summary_by_seed) || !nrow(params_long)) return(data.frame())
  params <- params_long
  params$transformed_value <- mapply(o2ipa_transform_parameter_value, params$parameter, params$value)
  params_wide <- o2ipa_params_wide(params, value_col = "transformed_value")
  params_wide$seed_id <- rownames(params_wide)
  feature_cols <- c(
    "terminal_mean_ploidy", "terminal_fraction_N_le_25", "terminal_fraction_N_below_44",
    "time_crossing_ploidy_1p5_down_censored", "terminal_minus_dominant_ploidy"
  )
  merged <- merge(summary_by_seed[, c("seed_id", "trajectory_regime", "O2_pct", "initial_condition", feature_cols), drop = FALSE], params_wide, by = "seed_id", all.x = TRUE)
  scopes <- list(
    all_seeds = merged,
    mode1_mode2_only = merged[merged$trajectory_regime %in% c("mode1_ploidy_stable", "mode2_second_ploidy_collapse"), , drop = FALSE]
  )
  rows <- list()
  counter <- 0L
  for (scope_name in names(scopes)) {
    ds <- scopes[[scope_name]]
    for (O2 in sort(unique(ds$O2_pct))) {
      for (init in unique(ds$initial_condition)) {
        dd <- ds[ds$O2_pct == O2 & ds$initial_condition == init, , drop = FALSE]
        for (feature in feature_cols) {
          y <- dd[[feature]]
          for (param in o2ipa_target_params()) {
            x <- dd[[param]]
            ok <- is.finite(x) & is.finite(y)
            if (sum(ok) < 5L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) next
            ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
            counter <- counter + 1L
            rows[[counter]] <- data.frame(
              correlation_scope = scope_name,
              O2_pct = O2,
              initial_condition = init,
              feature = feature,
              parameter = param,
              parameter_transform = o2ipa_optimizer_transform(param),
              n = sum(ok),
              spearman_rho = unname(ct$estimate),
              abs_spearman_rho = abs(unname(ct$estimate)),
              p_value = ct$p.value,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$BH_adjusted_p_value <- ave(out$p_value, out$correlation_scope, FUN = function(p) stats::p.adjust(p, method = "BH"))
  out[order(out$correlation_scope, -out$abs_spearman_rho, out$BH_adjusted_p_value), , drop = FALSE]
}

cf2_plot <- function(trajectory, summary_by_seed, fig_dir) {
  d <- trajectory[trajectory$trajectory_regime %in% c("mode1_ploidy_stable", "mode2_second_ploidy_collapse"), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  median_pdf_open <- FALSE
  grDevices::pdf(file.path(fig_dir, "fixed_o2_counterfactual_median_trajectories.pdf"), width = 10, height = 8, onefile = TRUE, bg = "white")
  median_pdf_open <- TRUE
  on.exit(if (isTRUE(median_pdf_open)) grDevices::dev.off(), add = TRUE)
  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02")
  line_cols <- setNames(grDevices::adjustcolor(unname(cols), alpha.f = 0.22), names(cols))
  reg_labels <- c(
    mode1_ploidy_stable = "Mode 1 ploidy stable",
    mode2_second_ploidy_collapse = "Mode 2 second ploidy collapse"
  )
  for (O2 in o2_use) {
    d_o2 <- d[d$O2_pct == O2, , drop = FALSE]
    par(mfrow = c(length(init_use), 1), mar = c(4, 5, 3, 1), oma = c(0, 0, 3, 0))
    for (init in init_use) {
      plot(NA, xlim = range(d$day), ylim = range(d$mean_ploidy, na.rm = TRUE),
           xlab = "Day", ylab = "Mean ploidy", main = paste("Initial condition:", init))
      abline(h = c(1.5, 2), col = "grey80", lty = 2)
      for (reg in names(cols)) {
        sub <- d_o2[d_o2$initial_condition == init & d_o2$trajectory_regime == reg, , drop = FALSE]
        if (!nrow(sub)) next
        seed_traces <- split(sub, sub$seed_id)
        for (trace in seed_traces) {
          trace <- trace[order(trace$day), , drop = FALSE]
          lines(trace$day, trace$mean_ploidy, col = line_cols[[reg]], lwd = 0.25)
        }
      }
      legend("topright", legend = unname(reg_labels[names(cols)]), col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.9)
    }
    mtext(paste0("Fixed-O2 seed trajectories (O2 = ", O2, "%)"), outer = TRUE, cex = 1.2, font = 2)
  }
  grDevices::dev.off()
  median_pdf_open <- FALSE

  grDevices::pdf(file.path(fig_dir, "fixed_o2_counterfactual_terminal_boxplots.pdf"), width = 10, height = 7, bg = "white")
  d2 <- summary_by_seed[summary_by_seed$trajectory_regime %in% c("mode1_ploidy_stable", "mode2_second_ploidy_collapse"), , drop = FALSE]
  boxplot(terminal_mean_ploidy ~ interaction(O2_pct, initial_condition, trajectory_regime, drop = TRUE),
          data = d2, las = 2, ylab = "Terminal mean ploidy", main = "Fixed-O2 counterfactual terminal ploidy")
  grDevices::dev.off()
  invisible(TRUE)
}

cf2_report <- function(out_dir, args, regime_summary, tests, correlations) {
  key_summary <- regime_summary[regime_summary$trajectory_regime %in% c("mode1_ploidy_stable", "mode2_second_ploidy_collapse") &
                                  regime_summary$O2_pct %in% c(1, 2, 5), , drop = FALSE]
  key_tests <- tests[tests$metric %in% c("terminal_mean_ploidy", "time_crossing_ploidy_1p5_down_censored") &
                       tests$O2_pct %in% c(1, 2, 5), , drop = FALSE]
  key_tests <- key_tests[order(key_tests$BH_adjusted_p_value), , drop = FALSE]
  top_corr <- data.frame()
  if (is.data.frame(correlations) && nrow(correlations) &&
      all(c("correlation_scope", "O2_pct", "abs_spearman_rho", "BH_adjusted_p_value") %in% names(correlations))) {
    top_corr <- correlations[correlations$correlation_scope == "mode1_mode2_only" & correlations$O2_pct %in% c(1, 2, 5), , drop = FALSE]
    if (nrow(top_corr)) {
      top_corr <- head(top_corr[order(-top_corr$abs_spearman_rho, top_corr$BH_adjusted_p_value), , drop = FALSE], 30L)
    }
  }
  lines <- c(
    "# In vivo fixed-O2 trajectory counterfactual analysis",
    "",
    paste0("- run_dir: `", o2ipa_as_chr(args$run_dir, ""), "`"),
    paste0("- label_file: `", o2ipa_as_chr(args$label_file, ""), "`"),
    paste0("- O2 grid: `", o2ipa_as_chr(args$o2_grid, "0,0.5,1,2,5"), "`"),
    paste0("- time grid: default dense early grid unless `--time_grid` was supplied"),
    "",
    "## Key high-O2 terminal summaries",
    "",
    paste(capture.output(print(key_summary, row.names = FALSE)), collapse = "\n"),
    "",
    "## Key mode1-vs-mode2 tests",
    "",
    paste(capture.output(print(head(key_tests, 30L), row.names = FALSE)), collapse = "\n"),
    "",
    "## Top mode1/mode2 parameter correlations at O2=1/2/5",
    "",
    paste(capture.output(print(top_corr, row.names = FALSE)), collapse = "\n"),
    "",
    "## Notes",
    "",
    "- This counterfactual uses the same fixed-O2 matrix as the attractor analysis.",
    "- Initial conditions are point masses at N=44 (`init_2N`) and N=88 (`init_4N`) unless clipped by the configured N grid.",
    "- State vectors are propagated through the exact eigen-decomposition of the fixed linear system and normalized to composition at each time point.",
    "- `time_crossing_ploidy_1p5_down_censored` is set to max_time + 1 when the trajectory never crosses 1.5."
  )
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}

main <- function() {
  args <- o2ipa_parse_args()
  repo_root <- o2ipa_repo_root(SCRIPT_DIR)
  run_dir <- normalizePath(o2ipa_as_chr(args$run_dir, file.path(repo_root, "oxygen", "results", "fit_invivo_O2_buffering_500seed")), mustWork = FALSE)
  label_file <- normalizePath(o2ipa_as_chr(args$label_file, file.path(repo_root, "oxygen", "results", "analysis", "invivo_o2_ploidy_event_coupling_500seed", "tables", "seed_event_summary.tsv")), mustWork = FALSE)
  out_dir <- normalizePath(o2ipa_as_chr(args$out_dir, file.path(repo_root, "oxygen", "results", "analysis", "invivo_fixed_o2_counterfactual_trajectories_500seed")), mustWork = FALSE)
  o2_grid <- sort(unique(cf2_as_num_vec(args$o2_grid, c(0, 0.5, 1, 2, 5))))
  time_grid <- sort(unique(cf2_as_num_vec(args$time_grid, cf2_default_time_grid())))
  max_seeds <- o2ipa_as_int(args$max_seeds, 0L)
  generate_figures <- o2ipa_as_bool(args$generate_figures, TRUE)
  max_time <- max(time_grid, na.rm = TRUE)

  cf2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "analyze_invivo_fixed_o2_counterfactual_trajectories.log")
  sink(log_file, split = TRUE)
  on.exit(sink(), add = TRUE)
  message("run_dir: ", run_dir)
  message("label_file: ", label_file)
  message("out_dir: ", out_dir)
  message("O2 grid: ", paste(o2_grid, collapse = ","))
  message("time grid length: ", length(time_grid), "; max_time=", max_time)

  seed_inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  manifest <- seed_inputs$manifest
  params_long <- seed_inputs$params_long
  if (max_seeds > 0L && nrow(manifest) > max_seeds) {
    manifest <- manifest[seq_len(max_seeds), , drop = FALSE]
    params_long <- params_long[params_long$seed_id %in% manifest$seed_id, , drop = FALSE]
  }
  labels <- cf2_load_labels(label_file)
  if (nrow(labels)) {
    manifest <- merge(manifest, labels, by = "seed_id", all.x = TRUE, sort = FALSE)
  }
  cfg <- o2pr_first_seed_cfg(manifest)
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  model_env <- o2ipa_source_model(SCRIPT_DIR)
  param_mat <- o2ipa_params_wide(params_long, "value")

  traj_rows <- list()
  summary_rows <- list()
  matrix_rows <- list()
  counter_traj <- 0L
  counter_summary <- 0L
  counter_matrix <- 0L
  init_specs <- data.frame(
    initial_condition = c("init_2N", "init_4N"),
    requested_N = c(44, 88),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(manifest))) {
    seed <- manifest$seed_id[[i]]
    if (i == 1L || i %% 25L == 0L) message("Processing ", seed, " (", i, "/", nrow(manifest), ")")
    if (!seed %in% rownames(param_mat)) next
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    for (O2 in o2_grid) {
      fm <- tryCatch(cf2_fixed_matrix(model_env, cfg, run_params, O2), error = function(e) NULL)
      if (is.null(fm)) next
      dom <- cf2_dominant_one(fm$M, fm$ngrid, n_unit)
      counter_matrix <- counter_matrix + 1L
      matrix_rows[[counter_matrix]] <- data.frame(
        seed_id = seed,
        O2_pct = O2,
        dom,
        stringsAsFactors = FALSE
      )
      for (j in seq_len(nrow(init_specs))) {
        init <- cf2_init_vector(fm$ngrid, init_specs$requested_N[[j]])
        sim <- cf2_eigen_trajectory(fm$M, fm$ngrid, init$vector, time_grid, n_unit)
        tr <- sim$trajectory
        if (nrow(tr)) {
          tr$seed_id <- seed
          tr$O2_pct <- O2
          tr$initial_condition <- init_specs$initial_condition[[j]]
          tr$requested_initial_N <- init_specs$requested_N[[j]]
          tr$used_initial_N <- init$used_N
          tr$status <- sim$status
          tr$trajectory_regime <- manifest$trajectory_regime[match(seed, manifest$seed_id)]
          tr$mode_label <- manifest$mode_label[match(seed, manifest$seed_id)]
          counter_traj <- counter_traj + 1L
          traj_rows[[counter_traj]] <- tr[, c(
            "seed_id", "trajectory_regime", "mode_label", "O2_pct", "initial_condition",
            "requested_initial_N", "used_initial_N", "status", "day",
            "mean_N", "mean_ploidy", "fraction_N_le_25", "fraction_N_below_44",
            "fraction_N_ge_44", "fraction_N_ge_66", "fraction_N_ge_88"
          )]
        }
        sm <- cf2_summarize_trajectory(tr, max_time)
        counter_summary <- counter_summary + 1L
        summary_rows[[counter_summary]] <- data.frame(
          seed_id = seed,
          trajectory_regime = manifest$trajectory_regime[match(seed, manifest$seed_id)],
          mode_label = manifest$mode_label[match(seed, manifest$seed_id)],
          O2_pct = O2,
          initial_condition = init_specs$initial_condition[[j]],
          requested_initial_N = init_specs$requested_N[[j]],
          used_initial_N = init$used_N,
          status = sim$status,
          sm,
          dom,
          terminal_minus_dominant_ploidy = sm$terminal_mean_ploidy - dom$dominant_mean_ploidy,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  trajectory <- if (length(traj_rows)) do.call(rbind, traj_rows) else data.frame()
  summary_by_seed <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
  dominant_by_seed <- if (length(matrix_rows)) do.call(rbind, matrix_rows) else data.frame()
  regime_summary <- cf2_regime_summary(summary_by_seed)
  tests <- cf2_regime_tests(summary_by_seed)
  correlations <- cf2_parameter_correlations(summary_by_seed, params_long)

  cf2_write_tsv(data.frame(
    argument = c("run_dir", "label_file", "out_dir", "o2_grid", "time_grid", "max_seeds", "generate_figures"),
    value = c(run_dir, label_file, out_dir, paste(o2_grid, collapse = ","), paste(time_grid, collapse = ","), as.character(max_seeds), as.character(generate_figures)),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  cf2_write_tsv(trajectory, file.path(out_dir, "tables", "fixed_o2_counterfactual_trajectories.tsv"))
  cf2_write_tsv(summary_by_seed, file.path(out_dir, "tables", "fixed_o2_counterfactual_summary_by_seed.tsv"))
  cf2_write_tsv(dominant_by_seed, file.path(out_dir, "tables", "fixed_o2_counterfactual_dominant_by_seed.tsv"))
  cf2_write_tsv(regime_summary, file.path(out_dir, "tables", "fixed_o2_counterfactual_regime_summary.tsv"))
  cf2_write_tsv(tests, file.path(out_dir, "tables", "fixed_o2_counterfactual_regime_tests.tsv"))
  cf2_write_tsv(correlations, file.path(out_dir, "tables", "fixed_o2_counterfactual_parameter_correlations.tsv"))
  cf2_write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))
  if (isTRUE(generate_figures)) cf2_plot(trajectory, summary_by_seed, file.path(out_dir, "figures"))
  cf2_report(out_dir, args, regime_summary, tests, correlations)
  message("Completed fixed-O2 counterfactual trajectory analysis: ", out_dir)
}

if (identical(environment(), globalenv())) {
  main()
}
