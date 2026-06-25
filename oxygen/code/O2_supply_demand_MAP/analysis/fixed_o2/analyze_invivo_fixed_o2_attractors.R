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

source(file.path(SCRIPT_DIR, "..", "process_fingerprints", "process_fingerprint_utils.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "..", "process_fingerprints", "ploidy_regime_utils.R"), local = TRUE)

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required")
})

options(error = function() {
  traceback(2)
  quit(status = 1)
})

fo2_as_num_vec <- function(x, default) {
  txt <- o2ipa_as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

fo2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

fo2_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

fo2_mkdirs <- function(out_dir) {
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs", "report")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

fo2_metric_map <- function(path) {
  if (!file.exists(path)) return(list())
  tab <- tryCatch(fo2_read_tsv(path), error = function(e) data.frame())
  if (!all(c("metric", "value") %in% names(tab))) return(list())
  vals <- as.list(tab$value)
  names(vals) <- tab$metric
  vals
}

fo2_map_num <- function(map, keys) {
  for (key in keys) {
    val <- suppressWarnings(as.numeric(map[[key]]))
    if (length(val) && is.finite(val[[1]])) return(val[[1]])
  }
  NA_real_
}

fo2_map_chr <- function(map, keys) {
  for (key in keys) {
    val <- map[[key]]
    if (!is.null(val) && length(val) && !is.na(val[[1]]) && nzchar(as.character(val[[1]]))) return(as.character(val[[1]]))
  }
  NA_character_
}

fo2_seed_manifest_extras <- function(manifest) {
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    seed <- manifest$seed_id[[i]]
    seed_dir <- manifest$seed_dir[[i]]
    m <- if (!is.na(seed_dir) && nzchar(seed_dir)) fo2_metric_map(file.path(seed_dir, "fit_summary.tsv")) else list()
    data.frame(
      seed_id = seed,
      optimizer_iter_completed = fo2_map_num(m, c("optimizer_iter_completed", "deoptim_iter_completed", "itermax")),
      optimizer_stop_reason = fo2_map_chr(m, c("deoptim_stop_reason", "optimizer_stop_reason", "optimizer_local_convergence")),
      optimizer_local_objective = fo2_map_num(m, c("optimizer_local_objective")),
      optimizer_deoptim_objective = fo2_map_num(m, c("optimizer_deoptim_objective")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

fo2_load_labels <- function(path, manifest) {
  if (nzchar(path) && file.exists(path)) {
    tab <- fo2_read_tsv(path)
    if (!"seed_id" %in% names(tab)) stop("label_file must contain seed_id: ", path)
    keep <- intersect(c(
      "seed_id", "trajectory_regime", "mode_label", "mean_terminal_ploidy",
      "mean_late_drop_amplitude", "4N__terminal_median_ploidy",
      "4N__late_drop_amplitude", "4N__time_o2_near_floor"
    ), names(tab))
    return(tab[, keep, drop = FALSE])
  }
  data.frame(
    seed_id = manifest$seed_id,
    trajectory_regime = NA_character_,
    mode_label = NA_character_,
    stringsAsFactors = FALSE
  )
}

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

fo2_wilcox_row <- function(tab, value_col, group_col = "trajectory_regime",
                           a = "mode1_ploidy_stable", b = "mode2_second_ploidy_collapse",
                           feature_name = value_col) {
  d <- tab[tab[[group_col]] %in% c(a, b) & is.finite(tab[[value_col]]), , drop = FALSE]
  if (!nrow(d) || length(unique(d[[group_col]])) < 2L) {
    return(data.frame(feature = feature_name, n_mode1 = 0L, n_mode2 = 0L, median_mode1 = NA_real_, median_mode2 = NA_real_, median_diff_mode2_minus_mode1 = NA_real_, wilcox_p_value = NA_real_))
  }
  x <- d[[value_col]][d[[group_col]] == a]
  y <- d[[value_col]][d[[group_col]] == b]
  p <- if (length(x) >= 1L && length(y) >= 1L) suppressWarnings(stats::wilcox.test(x, y)$p.value) else NA_real_
  data.frame(
    feature = feature_name,
    n_mode1 = length(x),
    n_mode2 = length(y),
    median_mode1 = stats::median(x, na.rm = TRUE),
    median_mode2 = stats::median(y, na.rm = TRUE),
    median_diff_mode2_minus_mode1 = stats::median(y, na.rm = TRUE) - stats::median(x, na.rm = TRUE),
    wilcox_p_value = p,
    stringsAsFactors = FALSE
  )
}

fo2_attractor_regime_tests <- function(attractors) {
  metrics <- c(
    "dominant_mean_ploidy", "dominant_fraction_N_le_25",
    "dominant_fraction_N_below_44", "dominant_growth_rate",
    "spectral_gap", "selection_22_vs_44", "selection_44_vs_88"
  )
  rows <- list()
  for (o2 in sort(unique(attractors$O2_pct))) {
    d <- attractors[attractors$O2_pct == o2, , drop = FALSE]
    for (m in metrics) {
      r <- fo2_wilcox_row(d, m, feature_name = paste0(m, "_O2_", gsub("[^0-9]+", "p", o2)))
      r$O2_pct <- o2
      r$metric <- m
      rows[[paste(o2, m, sep = "||")]] <- r
    }
  }
  out <- do.call(rbind, rows)
  out$BH_adjusted_p_value <- stats::p.adjust(out$wilcox_p_value, method = "BH")
  out[, c("O2_pct", "metric", "feature", "n_mode1", "n_mode2", "median_mode1", "median_mode2", "median_diff_mode2_minus_mode1", "wilcox_p_value", "BH_adjusted_p_value")]
}

fo2_parameter_correlations <- function(attractors, params_long) {
  if (!nrow(params_long)) return(data.frame())
  params_wide <- o2ipa_params_wide(params_long, "value")
  params_wide$seed_id <- rownames(params_wide)
  rows <- list()
  for (o2 in sort(unique(attractors$O2_pct))) {
    d0 <- attractors[attractors$O2_pct == o2, , drop = FALSE]
    d <- merge(d0, params_wide, by = "seed_id", all.x = TRUE, sort = FALSE)
    for (p in setdiff(names(params_wide), "seed_id")) {
      x <- suppressWarnings(as.numeric(d[[p]]))
      xtr <- if (all(x > 0, na.rm = TRUE)) log10(x) else x
      for (metric in c("dominant_mean_ploidy", "dominant_fraction_N_le_25", "selection_22_vs_44", "selection_44_vs_88")) {
        y <- suppressWarnings(as.numeric(d[[metric]]))
        ok <- is.finite(xtr) & is.finite(y)
        rows[[paste(o2, p, metric, sep = "||")]] <- data.frame(
          O2_pct = o2,
          parameter = p,
          metric = metric,
          n = sum(ok),
          spearman_rho = if (sum(ok) >= 3L) suppressWarnings(stats::cor(xtr[ok], y[ok], method = "spearman")) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$abs_spearman_rho <- abs(out$spearman_rho)
  out
}

fo2_mode_seed_stack_table <- function(attractors) {
  regimes <- c("mode1_ploidy_stable", "mode2_second_ploidy_collapse")
  ok_status <- if ("status" %in% names(attractors)) attractors$status == "ok" else rep(TRUE, nrow(attractors))
  d <- attractors[
    ok_status &
      attractors$trajectory_regime %in% regimes &
      is.finite(attractors$dominant_mean_ploidy),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(data.frame())
  d$seed_number <- o2ipa_seed_number(d$seed_id)
  d$stack_panel <- ifelse(
    d$trajectory_regime == "mode1_ploidy_stable",
    "mode1",
    "mode2"
  )
  keep <- intersect(c(
    "stack_panel", "seed_id", "seed_number", "trajectory_regime", "mode_label",
    "O2_pct", "dominant_mean_ploidy", "dominant_mean_N",
    "dominant_fraction_N_le_25", "dominant_fraction_N_below_44",
    "dominant_fraction_N_ge_44", "dominant_growth_rate", "spectral_gap",
    "selection_22_vs_44", "selection_44_vs_88",
    "objective", "delta_objective", "mean_terminal_ploidy", "mean_late_drop_amplitude"
  ), names(d))
  d <- d[, keep, drop = FALSE]
  d$trajectory_regime <- factor(d$trajectory_regime, levels = regimes)
  d <- d[order(d$trajectory_regime, d$seed_number, d$O2_pct), , drop = FALSE]
  d$trajectory_regime <- as.character(d$trajectory_regime)
  d
}

fo2_write_mode_comparison_tables <- function(attractors, out_dir) {
  regimes <- c(mode1 = "mode1_ploidy_stable", mode2 = "mode2_second_ploidy_collapse")
  ok_status <- if ("status" %in% names(attractors)) attractors$status == "ok" else rep(TRUE, nrow(attractors))
  d <- attractors[
    ok_status &
      attractors$trajectory_regime %in% unname(regimes) &
      is.finite(attractors$dominant_mean_ploidy),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(invisible(FALSE))
  o2_vals <- sort(unique(d$O2_pct))
  rows <- lapply(o2_vals, function(O2) {
    row <- data.frame(O2_pct = O2, stringsAsFactors = FALSE)
    for (nm in names(regimes)) {
      vals <- d$dominant_mean_ploidy[d$O2_pct == O2 & d$trajectory_regime == regimes[[nm]]]
      row[[paste0(nm, "_n")]] <- length(vals)
      row[[paste0(nm, "_median")]] <- stats::median(vals, na.rm = TRUE)
      row[[paste0(nm, "_mean")]] <- mean(vals, na.rm = TRUE)
      row[[paste0(nm, "_q25")]] <- as.numeric(stats::quantile(vals, 0.25, na.rm = TRUE, names = FALSE))
      row[[paste0(nm, "_q75")]] <- as.numeric(stats::quantile(vals, 0.75, na.rm = TRUE, names = FALSE))
    }
    row
  })
  summary <- do.call(rbind, rows)
  fo2_write_tsv(summary, file.path(out_dir, "tables", "dominant_mean_ploidy_summary_mode1_mode2.tsv"))

  median_table <- data.frame(
    O2 = formatC(summary$O2_pct, format = "f", digits = 2),
    mode1 = round(summary$mode1_median, 3),
    mode2 = round(summary$mode2_median, 3),
    stringsAsFactors = FALSE
  )
  fo2_write_tsv(median_table, file.path(out_dir, "tables", "dominant_mean_ploidy_median_mode1_mode2.tsv"))

  display <- summary
  names(display)[names(display) == "O2_pct"] <- "O2"
  display$O2 <- formatC(display$O2, format = "f", digits = 2)
  numeric_cols <- setdiff(names(display), c("O2", "mode1_n", "mode2_n"))
  for (col in numeric_cols) display[[col]] <- round(display[[col]], 3)
  fo2_write_tsv(display, file.path(out_dir, "tables", "dominant_mean_ploidy_summary_mode1_mode2_display.tsv"))
  invisible(TRUE)
}

fo2_finite_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
}

fo2_finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

fo2_finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

fo2_spectral_gap_by_seed <- function(attractors) {
  ok_status <- if ("status" %in% names(attractors)) attractors$status == "ok" else rep(TRUE, nrow(attractors))
  d <- attractors[
    ok_status &
      is.finite(attractors$dominant_growth_rate) &
      is.finite(attractors$spectral_gap),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(data.frame())
  d$seed_number <- o2ipa_seed_number(d$seed_id)
  d$lambda1 <- d$dominant_growth_rate
  d$lambda2 <- d$dominant_growth_rate - d$spectral_gap
  d$relative_spectral_gap <- d$spectral_gap / pmax(abs(d$lambda1), .Machine$double.eps)
  d$relax_time_days <- ifelse(d$spectral_gap > 0, 1 / d$spectral_gap, NA_real_)
  d$time_to_10x_days <- ifelse(d$spectral_gap > 0, log(10) / d$spectral_gap, NA_real_)
  d$time_to_100x_days <- ifelse(d$spectral_gap > 0, log(100) / d$spectral_gap, NA_real_)
  d$log10_advantage_1000d <- d$spectral_gap * 1000 / log(10)
  d$dominance_class <- ifelse(
    !is.finite(d$spectral_gap) | d$spectral_gap <= 0, "non_positive",
    ifelse(d$spectral_gap < 0.001, "very_weak",
           ifelse(d$spectral_gap < 0.005, "weak",
                  ifelse(d$spectral_gap < 0.01, "moderate", "strong")))
  )
  d$dominance_class <- factor(d$dominance_class, levels = c("non_positive", "very_weak", "weak", "moderate", "strong"))
  keep <- intersect(c(
    "seed_id", "seed_number", "O2_pct", "trajectory_regime", "mode_label",
    "dominant_mean_ploidy", "lambda1", "lambda2", "spectral_gap",
    "relative_spectral_gap", "relax_time_days", "time_to_10x_days",
    "time_to_100x_days", "log10_advantage_1000d", "dominance_class",
    "dominant_fraction_N_le_25", "dominant_fraction_N_below_44",
    "objective", "delta_objective", "mean_terminal_ploidy", "mean_late_drop_amplitude"
  ), names(d))
  d <- d[, keep, drop = FALSE]
  d <- d[order(d$trajectory_regime, d$seed_number, d$O2_pct), , drop = FALSE]
  d$dominance_class <- as.character(d$dominance_class)
  d
}

fo2_spectral_gap_summary <- function(gap_by_seed) {
  if (!nrow(gap_by_seed)) return(data.frame())
  rows <- lapply(split(gap_by_seed, list(gap_by_seed$O2_pct, gap_by_seed$trajectory_regime), drop = TRUE), function(d) {
    data.frame(
      O2_pct = d$O2_pct[[1]],
      trajectory_regime = d$trajectory_regime[[1]],
      n_seed = length(unique(d$seed_id)),
      dominant_mean_ploidy_median = fo2_finite_median(d$dominant_mean_ploidy),
      spectral_gap_median = fo2_finite_median(d$spectral_gap),
      spectral_gap_q25 = fo2_finite_quantile(d$spectral_gap, 0.25),
      spectral_gap_q75 = fo2_finite_quantile(d$spectral_gap, 0.75),
      spectral_gap_min = suppressWarnings(min(d$spectral_gap, na.rm = TRUE)),
      relative_gap_median = fo2_finite_median(d$relative_spectral_gap),
      relax_time_days_median = fo2_finite_median(d$relax_time_days),
      time_to_10x_days_median = fo2_finite_median(d$time_to_10x_days),
      time_to_10x_days_q25 = fo2_finite_quantile(d$time_to_10x_days, 0.25),
      time_to_10x_days_q75 = fo2_finite_quantile(d$time_to_10x_days, 0.75),
      time_to_100x_days_median = fo2_finite_median(d$time_to_100x_days),
      log10_advantage_1000d_median = fo2_finite_median(d$log10_advantage_1000d),
      log10_advantage_1000d_q25 = fo2_finite_quantile(d$log10_advantage_1000d, 0.25),
      log10_advantage_1000d_q75 = fo2_finite_quantile(d$log10_advantage_1000d, 0.75),
      fraction_gap_lt_0p001 = mean(d$spectral_gap < 0.001, na.rm = TRUE),
      fraction_gap_lt_0p005 = mean(d$spectral_gap < 0.005, na.rm = TRUE),
      fraction_gap_lt_0p01 = mean(d$spectral_gap < 0.01, na.rm = TRUE),
      fraction_gap_ge_0p005 = mean(d$spectral_gap >= 0.005, na.rm = TRUE),
      fraction_gap_ge_0p01 = mean(d$spectral_gap >= 0.01, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$O2_pct, out$trajectory_regime), , drop = FALSE]
}

fo2_ploidy_gap_reliability_composite_table <- function(gap_by_seed) {
  regimes <- c("mode1_ploidy_stable", "mode2_second_ploidy_collapse")
  d <- gap_by_seed[
    gap_by_seed$trajectory_regime %in% regimes &
      is.finite(gap_by_seed$dominant_mean_ploidy) &
      is.finite(gap_by_seed$spectral_gap),
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(data.frame())
  keep <- intersect(c(
    "seed_id", "seed_number", "trajectory_regime", "mode_label", "O2_pct",
    "dominant_mean_ploidy", "lambda1", "lambda2", "spectral_gap",
    "relative_spectral_gap", "relax_time_days", "time_to_10x_days",
    "time_to_100x_days", "log10_advantage_1000d", "dominance_class"
  ), names(d))
  d <- d[, keep, drop = FALSE]
  d$trajectory_regime <- factor(d$trajectory_regime, levels = regimes)
  d <- d[order(d$trajectory_regime, d$seed_number, d$O2_pct), , drop = FALSE]
  d$trajectory_regime <- as.character(d$trajectory_regime)
  d
}

fo2_shade_uncertain_o2 <- function(summary, reg, y_rng, o2_vals) {
  if (!nrow(summary)) return(invisible(FALSE))
  sub <- summary[
    summary$trajectory_regime == reg &
      (
        (!is.na(summary$spectral_gap_median) & summary$spectral_gap_median < 0.005) |
          (!is.na(summary$time_to_10x_days_median) & summary$time_to_10x_days_median > 500)
      ),
    ,
    drop = FALSE
  ]
  if (!nrow(sub)) return(invisible(FALSE))
  step <- suppressWarnings(min(diff(sort(unique(o2_vals))), na.rm = TRUE))
  if (!is.finite(step) || step <= 0) step <- 0.05
  half_step <- step / 2
  rect(
    xleft = pmax(min(o2_vals, na.rm = TRUE), sub$O2_pct - half_step),
    ybottom = y_rng[[1]],
    xright = pmin(max(o2_vals, na.rm = TRUE), sub$O2_pct + half_step),
    ytop = y_rng[[2]],
    col = grDevices::adjustcolor("grey80", alpha.f = 0.35),
    border = NA
  )
  invisible(TRUE)
}

fo2_plot_ploidy_gap_reliability_composite <- function(gap_by_seed, summary, fig_dir) {
  regimes <- c("mode1_ploidy_stable", "mode2_second_ploidy_collapse")
  d <- fo2_ploidy_gap_reliability_composite_table(gap_by_seed)
  d <- d[d$trajectory_regime %in% regimes, , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  d_gap <- d[d$spectral_gap > 0, , drop = FALSE]
  d_time <- d[is.finite(d$time_to_10x_days) & d$time_to_10x_days > 0, , drop = FALSE]
  if (!nrow(d_gap) || !nrow(d_time)) return(invisible(FALSE))

  panel_names <- c(
    mode1_ploidy_stable = "Mode 1: ploidy stable",
    mode2_second_ploidy_collapse = "Mode 2: 2nd ploidy collapse"
  )
  cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02")
  o2_vals <- sort(unique(d$O2_pct))
  ploidy_rng <- range(d$dominant_mean_ploidy, 1, 1.5, 2, 2.5, na.rm = TRUE)
  ploidy_pad <- diff(ploidy_rng) * 0.04
  if (!is.finite(ploidy_pad) || ploidy_pad == 0) ploidy_pad <- 0.1
  ploidy_rng <- ploidy_rng + c(-ploidy_pad, ploidy_pad)
  gap_rng <- range(d_gap$spectral_gap, 0.001, 0.005, 0.01, na.rm = TRUE)
  time_rng <- range(d_time$time_to_10x_days, 100, 500, 1000, na.rm = TRUE)

  grDevices::pdf(file.path(fig_dir, "fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.pdf"), width = 11, height = 10)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2), mar = c(4.1, 4.4, 2.1, 0.8), oma = c(0.6, 0.4, 2, 0.4))
  seed_lwd <- 0.95
  for (metric_row in c("ploidy", "gap", "time10x")) {
    for (j in seq_along(regimes)) {
      reg <- regimes[[j]]
      sub <- d[d$trajectory_regime == reg, , drop = FALSE]
      seed_ids <- unique(sub$seed_id[order(sub$seed_number)])
      seed_col <- grDevices::adjustcolor(cols[[reg]], alpha.f = 0.18)
      main_label <- if (metric_row == "ploidy") panel_names[[reg]] else ""
      x_axis <- if (metric_row == "time10x") "s" else "n"
      y_axis <- if (j == 1L) "s" else "n"
      if (metric_row == "ploidy") {
        plot(NA, xlim = range(o2_vals), ylim = ploidy_rng,
             xlab = "", ylab = if (j == 1L) "Dominant mean ploidy" else "",
             main = main_label, xaxt = x_axis, yaxt = y_axis)
        fo2_shade_uncertain_o2(summary, reg, ploidy_rng, o2_vals)
        for (seed in seed_ids) {
          ss <- sub[sub$seed_id == seed, , drop = FALSE]
          ss <- ss[order(ss$O2_pct), , drop = FALSE]
          lines(ss$O2_pct, ss$dominant_mean_ploidy, col = seed_col, lwd = seed_lwd)
        }
        abline(h = c(1, 1.5, 2, 2.5), col = "grey62", lty = 3, lwd = 1.4)
        if (j == 1L) {
          legend("topright", legend = c("uncertain O2 band", "single seed", "reference line"),
                 fill = c(grDevices::adjustcolor("grey80", alpha.f = 0.35), NA, NA),
                 border = c(NA, NA, NA), lty = c(NA, 1, 3),
                 col = c(NA, seed_col, "grey62"), lwd = c(NA, seed_lwd, 1.4),
                 bty = "n", cex = 0.75)
        }
      } else if (metric_row == "gap") {
        plot(NA, xlim = range(o2_vals), ylim = gap_rng, log = "y",
             xlab = "", ylab = if (j == 1L) "Spectral gap" else "",
             main = main_label, xaxt = x_axis, yaxt = y_axis)
        fo2_shade_uncertain_o2(summary, reg, gap_rng, o2_vals)
        for (seed in seed_ids) {
          ss <- sub[sub$seed_id == seed & sub$spectral_gap > 0, , drop = FALSE]
          ss <- ss[order(ss$O2_pct), , drop = FALSE]
          lines(ss$O2_pct, ss$spectral_gap, col = seed_col, lwd = seed_lwd)
        }
        abline(h = c(0.001, 0.005, 0.01), col = "grey55", lty = 3, lwd = 1.4)
      } else {
        plot(NA, xlim = range(o2_vals), ylim = time_rng, log = "y",
             xlab = "", ylab = if (j == 1L) "Days to 10x dominance" else "",
             main = main_label, xaxt = x_axis, yaxt = y_axis)
        fo2_shade_uncertain_o2(summary, reg, time_rng, o2_vals)
        for (seed in seed_ids) {
          ss <- sub[sub$seed_id == seed & is.finite(sub$time_to_10x_days) & sub$time_to_10x_days > 0, , drop = FALSE]
          ss <- ss[order(ss$O2_pct), , drop = FALSE]
          lines(ss$O2_pct, ss$time_to_10x_days, col = seed_col, lwd = seed_lwd)
        }
        abline(h = c(100, 500, 1000), col = "grey55", lty = 3, lwd = 1.4)
        mtext("Fixed O2 pct", side = 1, line = 2.7, cex = 0.9)
      }
    }
  }
  mtext("Fixed-O2 ploidy attractor and spectral-gap reliability", outer = TRUE, cex = 1.1, font = 2)
  par(oldpar)
  grDevices::dev.off()
  invisible(TRUE)
}

fo2_ploidy_gap_reliability_violin_table <- function(gap_by_seed,
                                                    o2_values = c(0, 0.1, 0.2, 0.5, 0.75, 1, 2, 3, 4, 5)) {
  regimes <- c(mode1 = "mode1_ploidy_stable", mode2 = "mode2_second_ploidy_collapse")
  d <- fo2_ploidy_gap_reliability_composite_table(gap_by_seed)
  if (!nrow(d)) return(data.frame())
  target <- round(o2_values, 8)
  key <- round(d$O2_pct, 8)
  keep <- key %in% target & d$trajectory_regime %in% unname(regimes)
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  d$O2_index <- match(round(d$O2_pct, 8), target)
  d$O2_label <- c("0", "0.1", "0.2", "0.5", "0.75", "1", "2", "3", "4", "5")[d$O2_index]
  d$mode_group <- ifelse(d$trajectory_regime == regimes[["mode1"]], "mode1", "mode2")
  d <- d[order(d$O2_index, d$mode_group, d$seed_number), , drop = FALSE]
  keep_cols <- intersect(c(
    "O2_index", "O2_label", "O2_pct", "mode_group", "seed_id", "seed_number",
    "trajectory_regime", "mode_label", "dominant_mean_ploidy", "spectral_gap",
    "time_to_10x_days", "time_to_100x_days", "log10_advantage_1000d",
    "dominance_class", "lambda1", "lambda2", "relative_spectral_gap",
    "relax_time_days"
  ), names(d))
  d[, keep_cols, drop = FALSE]
}

fo2_draw_violin_box <- function(y, x, width, fill_col, border_col = "grey25") {
  y <- y[is.finite(y)]
  if (!length(y)) return(invisible(FALSE))
  if (length(unique(y)) >= 2L) {
    dens <- stats::density(y, n = 128)
    scale <- if (max(dens$y, na.rm = TRUE) > 0) dens$y / max(dens$y, na.rm = TRUE) else dens$y
    graphics::polygon(
      c(x - scale * width, rev(x + scale * width)),
      c(dens$x, rev(dens$x)),
      col = fill_col,
      border = border_col,
      lwd = 0.7
    )
  } else {
    graphics::segments(x - width * 0.55, y[[1]], x + width * 0.55, y[[1]], col = border_col, lwd = 1.2)
  }
  graphics::boxplot(
    y, at = x, add = TRUE, axes = FALSE, outline = FALSE,
    boxwex = width * 0.55, col = grDevices::adjustcolor("white", alpha.f = 0.75),
    border = border_col, staplewex = 0.55
  )
  invisible(TRUE)
}

fo2_axis_ticks_log10 <- function(vals, rng) {
  vals <- vals[is.finite(vals) & vals > 0]
  vals <- vals[log10(vals) >= rng[[1]] & log10(vals) <= rng[[2]]]
  vals
}

fo2_plot_ploidy_gap_reliability_violin <- function(gap_by_seed, fig_dir,
                                                   o2_values = c(0, 0.1, 0.2, 0.5, 0.75, 1, 2, 3, 4, 5)) {
  d <- fo2_ploidy_gap_reliability_violin_table(gap_by_seed, o2_values = o2_values)
  if (!nrow(d)) return(invisible(FALSE))
  o2_labels <- c("0", "0.1", "0.2", "0.5", "0.75", "1", "2", "3", "4", "5")
  x_centers <- seq_along(o2_labels)
  offsets <- c(mode1 = -0.18, mode2 = 0.18)
  cols <- c(mode1 = "#1b9e77", mode2 = "#d95f02")
  fill_cols <- stats::setNames(grDevices::adjustcolor(cols, alpha.f = 0.45), names(cols))
  violin_width <- 0.16

  ploidy_rng <- range(d$dominant_mean_ploidy, c(1, 1.5, 2, 2.5), na.rm = TRUE)
  ploidy_pad <- diff(ploidy_rng) * 0.05
  if (!is.finite(ploidy_pad) || ploidy_pad == 0) ploidy_pad <- 0.1
  ploidy_rng <- ploidy_rng + c(-ploidy_pad, ploidy_pad)

  gap_vals <- d$spectral_gap[d$spectral_gap > 0]
  time_vals <- d$time_to_10x_days[is.finite(d$time_to_10x_days) & d$time_to_10x_days > 0]
  gap_rng <- range(log10(gap_vals), log10(c(0.001, 0.005, 0.01)), na.rm = TRUE)
  time_rng <- range(log10(time_vals), log10(c(100, 500, 1000)), na.rm = TRUE)

  grDevices::pdf(file.path(fig_dir, "fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.pdf"), width = 11, height = 9.5)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(3, 1), mar = c(2.5, 5, 2.1, 1), oma = c(3.3, 0.5, 1.5, 0.5))

  panel_specs <- list(
    list(metric = "dominant_mean_ploidy", ylab = "Dominant mean ploidy", ylim = ploidy_rng, log10 = FALSE,
         ref = c(1, 1.5, 2, 2.5), ref_label = "ploidy reference"),
    list(metric = "spectral_gap", ylab = "Spectral gap", ylim = gap_rng, log10 = TRUE,
         ref = c(0.001, 0.005, 0.01), ref_label = "gap threshold"),
    list(metric = "time_to_10x_days", ylab = "Days to 10x dominance", ylim = time_rng, log10 = TRUE,
         ref = c(100, 500, 1000), ref_label = "time reference")
  )

  for (i in seq_along(panel_specs)) {
    spec <- panel_specs[[i]]
    xaxt <- if (i == length(panel_specs)) "s" else "n"
    plot(NA, xlim = c(0.5, length(o2_labels) + 0.5), ylim = spec$ylim,
         axes = FALSE, xlab = "", ylab = spec$ylab, main = "")
    if (spec$log10) {
      ticks <- if (spec$metric == "spectral_gap") {
        fo2_axis_ticks_log10(c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1), spec$ylim)
      } else {
        fo2_axis_ticks_log10(c(10, 30, 100, 300, 500, 1000, 3000, 10000), spec$ylim)
      }
      axis(2, at = log10(ticks), labels = format(ticks, trim = TRUE, scientific = FALSE), las = 1)
    } else {
      axis(2, las = 1)
    }
    if (xaxt == "s") axis(1, at = x_centers, labels = o2_labels)
    box()
    for (idx in x_centers) {
      for (mode in names(offsets)) {
        vals <- d[d$O2_index == idx & d$mode_group == mode, spec$metric]
        vals <- suppressWarnings(as.numeric(vals))
        if (spec$log10) vals <- log10(vals[is.finite(vals) & vals > 0])
        fo2_draw_violin_box(vals, idx + offsets[[mode]], violin_width, fill_cols[[mode]], cols[[mode]])
      }
    }
    ref_y <- if (spec$log10) log10(spec$ref) else spec$ref
    graphics::abline(h = ref_y, col = "grey45", lty = 3, lwd = 1.5)
    if (i == 1L) {
      legend("topright",
             legend = c("mode1", "mode2", spec$ref_label),
             fill = c(fill_cols[["mode1"]], fill_cols[["mode2"]], NA),
             border = c(cols[["mode1"]], cols[["mode2"]], NA),
             lty = c(NA, NA, 3), col = c(NA, NA, "grey45"),
             lwd = c(NA, NA, 1.5), bty = "n", cex = 0.85)
    }
  }
  mtext("Fixed O2 pct", side = 1, outer = TRUE, line = 1.8, cex = 0.95)
  mtext("Fixed-O2 ploidy attractor and reliability distributions", side = 3, outer = TRUE, line = 0.3, cex = 1.1, font = 2)
  par(oldpar)
  grDevices::dev.off()
  invisible(TRUE)
}

fo2_plot_spectral_gap_outputs <- function(gap_by_seed, summary, fig_dir) {
  if (!nrow(gap_by_seed)) return(invisible(FALSE))
  regimes <- c("mode1_ploidy_stable", "mode2_second_ploidy_collapse")
  panel_names <- c(
    mode1_ploidy_stable = "Mode 1: ploidy stable",
    mode2_second_ploidy_collapse = "Mode 2: 2nd ploidy collapse"
  )
  cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02")
  d <- gap_by_seed[gap_by_seed$trajectory_regime %in% regimes & gap_by_seed$spectral_gap > 0, , drop = FALSE]
  if (nrow(d)) {
    y_rng <- range(d$spectral_gap, 0.001, 0.005, 0.01, na.rm = TRUE)
    grDevices::pdf(file.path(fig_dir, "fixed_o2_spectral_gap_all_seed_stack_mode1_mode2.pdf"), width = 8.5, height = 8)
    oldpar <- par(no.readonly = TRUE)
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
    for (reg in regimes) {
      sub <- d[d$trajectory_regime == reg, , drop = FALSE]
      plot(NA, xlim = range(d$O2_pct), ylim = y_rng, log = "y",
           xlab = "Fixed O2 pct", ylab = "Spectral gap",
           main = panel_names[[reg]])
      abline(h = c(0.001, 0.005, 0.01), col = "grey80", lty = 3)
      seed_ids <- unique(sub$seed_id[order(sub$seed_number)])
      seed_col <- grDevices::adjustcolor(cols[[reg]], alpha.f = 0.18)
      for (seed in seed_ids) {
        ss <- sub[sub$seed_id == seed, , drop = FALSE]
        ss <- ss[order(ss$O2_pct), , drop = FALSE]
        lines(ss$O2_pct, ss$spectral_gap, col = seed_col, lwd = 0.7)
      }
      legend("topright", legend = c("0.001 / 0.005 / 0.01 thresholds", "single seed"),
             col = c("grey80", seed_col), lty = c(3, 1), lwd = c(1, 1), bty = "n")
    }
    par(oldpar)
    grDevices::dev.off()
  }

  d10 <- gap_by_seed[gap_by_seed$trajectory_regime %in% regimes & is.finite(gap_by_seed$time_to_10x_days) & gap_by_seed$time_to_10x_days > 0, , drop = FALSE]
  if (nrow(d10)) {
    y_rng <- range(d10$time_to_10x_days, 10, 100, 1000, na.rm = TRUE)
    grDevices::pdf(file.path(fig_dir, "fixed_o2_time_to_10x_all_seed_stack_mode1_mode2.pdf"), width = 8.5, height = 8)
    oldpar <- par(no.readonly = TRUE)
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
    for (reg in regimes) {
      sub <- d10[d10$trajectory_regime == reg, , drop = FALSE]
      plot(NA, xlim = range(d10$O2_pct), ylim = y_rng, log = "y",
           xlab = "Fixed O2 pct", ylab = "Days for dominant mode to reach 10x",
           main = panel_names[[reg]])
      abline(h = c(10, 100, 1000), col = "grey80", lty = 3)
      seed_ids <- unique(sub$seed_id[order(sub$seed_number)])
      seed_col <- grDevices::adjustcolor(cols[[reg]], alpha.f = 0.18)
      for (seed in seed_ids) {
        ss <- sub[sub$seed_id == seed, , drop = FALSE]
        ss <- ss[order(ss$O2_pct), , drop = FALSE]
        lines(ss$O2_pct, ss$time_to_10x_days, col = seed_col, lwd = 0.7)
      }
      legend("topright", legend = c("10 / 100 / 1000 days", "single seed"),
             col = c("grey80", seed_col), lty = c(3, 1), lwd = c(1, 1), bty = "n")
    }
    par(oldpar)
    grDevices::dev.off()
  }

  s <- summary[summary$trajectory_regime %in% regimes, , drop = FALSE]
  if (nrow(s)) {
    grDevices::pdf(file.path(fig_dir, "fixed_o2_gap_reliability_fraction_mode1_mode2.pdf"), width = 8.5, height = 6)
    oldpar <- par(no.readonly = TRUE)
    plot(NA, xlim = range(s$O2_pct), ylim = c(0, 1),
         xlab = "Fixed O2 pct", ylab = "Fraction of seeds",
         main = "Spectral gap reliability fractions")
    ltys <- c(fraction_gap_ge_0p005 = 2, fraction_gap_ge_0p01 = 1)
    for (reg in regimes) {
      sub <- s[s$trajectory_regime == reg, , drop = FALSE]
      sub <- sub[order(sub$O2_pct), , drop = FALSE]
      lines(sub$O2_pct, sub$fraction_gap_ge_0p005, col = cols[[reg]], lty = ltys[["fraction_gap_ge_0p005"]], lwd = 2)
      lines(sub$O2_pct, sub$fraction_gap_ge_0p01, col = cols[[reg]], lty = ltys[["fraction_gap_ge_0p01"]], lwd = 2)
    }
    legend("bottomright",
           legend = c("mode1 gap >= 0.005", "mode1 gap >= 0.01", "mode2 gap >= 0.005", "mode2 gap >= 0.01"),
           col = c(cols[["mode1_ploidy_stable"]], cols[["mode1_ploidy_stable"]], cols[["mode2_second_ploidy_collapse"]], cols[["mode2_second_ploidy_collapse"]]),
           lty = c(2, 1, 2, 1), lwd = 2, bty = "n")
    par(oldpar)
    grDevices::dev.off()
  }
  invisible(TRUE)
}

fo2_write_spectral_gap_outputs <- function(attractors, out_dir, generate_figures = TRUE) {
  gap_by_seed <- fo2_spectral_gap_by_seed(attractors)
  gap_summary <- fo2_spectral_gap_summary(gap_by_seed)
  fo2_write_tsv(gap_by_seed, file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_by_seed.tsv"))
  fo2_write_tsv(gap_summary, file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_regime_summary.tsv"))
  fo2_write_tsv(fo2_ploidy_gap_reliability_composite_table(gap_by_seed), file.path(out_dir, "tables", "fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv"))
  fo2_write_tsv(fo2_ploidy_gap_reliability_violin_table(gap_by_seed), file.path(out_dir, "tables", "fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv"))
  if (isTRUE(generate_figures)) {
    dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
    fo2_plot_spectral_gap_outputs(gap_by_seed, gap_summary, file.path(out_dir, "figures"))
    fo2_plot_ploidy_gap_reliability_composite(gap_by_seed, gap_summary, file.path(out_dir, "figures"))
    fo2_plot_ploidy_gap_reliability_violin(gap_by_seed, file.path(out_dir, "figures"))
  }
  invisible(list(by_seed = gap_by_seed, summary = gap_summary))
}

fo2_plot_mode_seed_stack <- function(attractors, fig_dir) {
  d <- fo2_mode_seed_stack_table(attractors)
  if (!nrow(d)) return(invisible(FALSE))
  regimes <- c("mode1_ploidy_stable", "mode2_second_ploidy_collapse")
  panel_names <- c(
    mode1_ploidy_stable = "Mode 1: ploidy stable",
    mode2_second_ploidy_collapse = "Mode 2: 2nd ploidy collapse"
  )
  cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02")
  x_vals <- sort(unique(d$O2_pct))
  y_rng <- range(d$dominant_mean_ploidy, na.rm = TRUE)
  y_pad <- diff(y_rng) * 0.04
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.1
  y_rng <- y_rng + c(-y_pad, y_pad)

  grDevices::pdf(file.path(fig_dir, "fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.pdf"), width = 8.5, height = 8)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 1), mar = c(4, 5, 3, 1))
  for (reg in regimes) {
    sub <- d[d$trajectory_regime == reg, , drop = FALSE]
    plot(NA,
         xlim = range(x_vals, na.rm = TRUE), ylim = y_rng, xaxt = "n",
         xlab = "Fixed O2 pct", ylab = "Dominant mean ploidy",
         main = panel_names[[reg]])
    axis(1, at = x_vals, labels = format(x_vals, trim = TRUE, scientific = FALSE))
    abline(h = c(1, 1.5, 2, 2.5), col = "grey85", lty = 3)
    seed_ids <- unique(sub$seed_id[order(sub$seed_number)])
    seed_col <- grDevices::adjustcolor(cols[[reg]], alpha.f = 0.18)
    for (seed in seed_ids) {
      ss <- sub[sub$seed_id == seed, , drop = FALSE]
      ss <- ss[order(ss$O2_pct), , drop = FALSE]
      lines(ss$O2_pct, ss$dominant_mean_ploidy, col = seed_col, lwd = 0.7)
    }
    legend("topright", legend = "single seed", col = seed_col, lwd = 1, bty = "n")
  }
  invisible(TRUE)
}

fo2_plot_outputs <- function(attractors, out_dir) {
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02", ambiguous = "grey60")
  d <- attractors[attractors$status == "ok" & is.finite(attractors$dominant_mean_ploidy), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  grDevices::pdf(file.path(fig_dir, "fixed_o2_dominant_ploidy_by_regime.pdf"), width = 8, height = 6)
  plot(NA, xlim = range(d$O2_pct, na.rm = TRUE), ylim = range(d$dominant_mean_ploidy, na.rm = TRUE),
       xlab = "Fixed O2 pct", ylab = "Dominant mean ploidy", main = "Fixed-O2 attractor by regime")
  for (reg in unique(d$trajectory_regime)) {
    sub <- d[d$trajectory_regime == reg, , drop = FALSE]
    if (!nrow(sub)) next
    by_o2 <- aggregate(sub$dominant_mean_ploidy, by = list(O2_pct = sub$O2_pct), FUN = median, na.rm = TRUE)
    lines(by_o2$O2_pct, by_o2$x, col = cols[reg] %||% "grey60", lwd = 2, type = "b", pch = 16)
  }
  legend("topright", legend = names(cols), col = cols, lwd = 2, pch = 16, bty = "n")
  grDevices::dev.off()

  low <- d[d$O2_pct %in% sort(unique(d$O2_pct))[seq_len(min(4L, length(unique(d$O2_pct))))], , drop = FALSE]
  if (nrow(low)) {
    grDevices::pdf(file.path(fig_dir, "low_o2_dominant_ploidy_boxplot.pdf"), width = 9, height = 6)
    boxplot(dominant_mean_ploidy ~ interaction(O2_pct, trajectory_regime, drop = TRUE), data = low, las = 2,
            ylab = "Dominant mean ploidy", main = "Low fixed-O2 attractor distribution")
    grDevices::dev.off()
  }
  fo2_plot_mode_seed_stack(attractors, fig_dir)
  invisible(TRUE)
}

fo2_write_report <- function(out_dir, run_dir, label_file, attractors, regime_summary, tests, correlations) {
  counts <- as.data.frame(table(attractors$trajectory_regime[!duplicated(attractors$seed_id)]), stringsAsFactors = FALSE)
  names(counts) <- c("trajectory_regime", "n_seed")
  low_o2 <- regime_summary[regime_summary$O2_pct %in% c(0, 0.01, 0.05, 0.1), , drop = FALSE]
  top_tests <- tests[is.finite(tests$BH_adjusted_p_value), , drop = FALSE]
  top_tests <- top_tests[order(top_tests$BH_adjusted_p_value, top_tests$wilcox_p_value), , drop = FALSE]
  top_corr <- correlations[is.finite(correlations$abs_spearman_rho), , drop = FALSE]
  top_corr <- top_corr[order(-top_corr$abs_spearman_rho), , drop = FALSE]
  gap_summary_path <- file.path(out_dir, "tables", "fixed_o2_attractor_spectral_gap_regime_summary.tsv")
  gap_summary <- if (file.exists(gap_summary_path)) fo2_read_tsv(gap_summary_path) else data.frame()
  key_gap <- gap_summary[
    nrow(gap_summary) > 0 &
      gap_summary$trajectory_regime %in% c("mode1_ploidy_stable", "mode2_second_ploidy_collapse") &
      gap_summary$O2_pct %in% c(0, 0.5, 1, 2, 5),
    ,
    drop = FALSE
  ]
  key_gap <- key_gap[, intersect(c(
    "O2_pct", "trajectory_regime", "n_seed", "spectral_gap_median",
    "time_to_10x_days_median", "time_to_100x_days_median",
    "log10_advantage_1000d_median", "fraction_gap_ge_0p005",
    "fraction_gap_ge_0p01"
  ), names(key_gap)), drop = FALSE]
  lines <- c(
    "# In vivo fixed-O2 ploidy attractor analysis",
    "",
    paste0("- run_dir: ", run_dir),
    paste0("- label_file: ", label_file),
    paste0("- analyzed seeds: ", length(unique(attractors$seed_id))),
    "",
    "## Regime counts",
    "",
    paste(utils::capture.output(print(counts, row.names = FALSE)), collapse = "\n"),
    "",
    "## Low-O2 Attractor Summary",
    "",
    paste(utils::capture.output(print(low_o2, row.names = FALSE)), collapse = "\n"),
    "",
    "## Strongest Mode1-vs-Mode2 Attractor Differences",
    "",
    paste(utils::capture.output(print(utils::head(top_tests, 20), row.names = FALSE)), collapse = "\n"),
    "",
    "## Strongest Parameter-Attractor Correlations",
    "",
    paste(utils::capture.output(print(utils::head(top_corr, 30), row.names = FALSE)), collapse = "\n"),
    "",
    "## Seed-Level Mode1/Mode2 Stack Outputs",
    "",
    "- Table: `tables/fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv`",
    "- Figure: `figures/fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.pdf`",
    "",
    "## Spectral Gap Reliability Outputs",
    "",
    "- Seed table: `tables/fixed_o2_attractor_spectral_gap_by_seed.tsv`",
    "- Regime summary: `tables/fixed_o2_attractor_spectral_gap_regime_summary.tsv`",
    "- Figures: `figures/fixed_o2_spectral_gap_all_seed_stack_mode1_mode2.pdf`, `figures/fixed_o2_time_to_10x_all_seed_stack_mode1_mode2.pdf`, `figures/fixed_o2_gap_reliability_fraction_mode1_mode2.pdf`",
    "- Composite table: `tables/fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.tsv`",
    "- Composite figure: `figures/fixed_o2_ploidy_gap_reliability_composite_mode1_mode2.pdf`",
    "- Violin table: `tables/fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.tsv`",
    "- Violin figure: `figures/fixed_o2_ploidy_gap_reliability_violin_mode1_mode2.pdf`",
    "",
    paste(utils::capture.output(print(key_gap, row.names = FALSE)), collapse = "\n"),
    "",
    "## Interpretation Boundary",
    "",
    "These are fixed-O2 attractors from the fitted model parameters. They test model-internal mechanism consistency under standardized O2; they do not by themselves prove biological causality."
  )
  dir.create(file.path(out_dir, "report"), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}

main <- function(argv = o2ipa_parse_args()) {
  repo_root <- o2ipa_repo_root(SCRIPT_DIR)
  run_dir <- normalizePath(o2ipa_as_chr(argv$run_dir), mustWork = FALSE)
  out_dir <- normalizePath(o2ipa_as_chr(argv$out_dir, file.path(repo_root, "oxygen", "results", "analysis", "invivo_fixed_o2_attractors_500seed")), mustWork = FALSE)
  label_file <- normalizePath(o2ipa_as_chr(argv$label_file, ""), mustWork = FALSE)
  if (!dir.exists(run_dir)) stop("run_dir does not exist: ", run_dir)
  fo2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "analyze_invivo_fixed_o2_attractors.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  objective_source <- o2ipa_as_chr(argv$objective_source, "auto")
  o2_grid <- sort(unique(fo2_as_num_vec(argv$o2_grid, seq(0, 5, by = 0.05))))
  max_seeds <- o2ipa_as_int(argv$max_seeds, NA_integer_)
  generate_figures <- o2ipa_as_bool(argv$generate_figures, TRUE)
  random_seed <- o2ipa_as_int(argv$random_seed, 20260623L)
  set.seed(random_seed)

  run_args <- data.frame(
    argument = c("run_dir", "out_dir", "label_file", "objective_source", "o2_grid", "max_seeds", "generate_figures", "random_seed"),
    value = c(run_dir, out_dir, label_file, objective_source, paste(o2_grid, collapse = ","), max_seeds, generate_figures, random_seed),
    stringsAsFactors = FALSE
  )
  fo2_write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)
  message("label_file: ", label_file)

  message("Collecting seed inputs.")
  inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = objective_source)
  if (!is.na(max_seeds) && max_seeds > 0L) {
    valid <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
    valid <- valid[order(o2ipa_seed_number(valid$seed_id)), , drop = FALSE]
    keep <- valid$seed_id[seq_len(min(max_seeds, nrow(valid)))]
    inputs$manifest <- inputs$manifest[inputs$manifest$seed_id %in% keep, , drop = FALSE]
    inputs$params_long <- inputs$params_long[inputs$params_long$seed_id %in% keep, , drop = FALSE]
    message("Using max_seeds subset: ", paste(keep, collapse = ", "))
  }
  extras <- fo2_seed_manifest_extras(inputs$manifest)
  manifest <- merge(inputs$manifest, extras, by = "seed_id", all.x = TRUE, sort = FALSE)
  manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  labels <- fo2_load_labels(label_file, manifest)
  manifest <- merge(manifest, labels, by = "seed_id", all.x = TRUE, sort = FALSE)
  fo2_write_tsv(manifest, file.path(out_dir, "tables", "seed_manifest_with_labels.tsv"))
  fo2_write_tsv(inputs$params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))

  cfg <- o2pr_first_seed_cfg(inputs$manifest)
  cfg_meta <- o2pr_cfg_metadata(cfg)
  fo2_write_tsv(cfg_meta, file.path(out_dir, "tables", "fit_config_schema_summary.tsv"))
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  model_env <- o2ipa_source_model(SCRIPT_DIR)

  message("Computing fixed-O2 attractors.")
  rows <- list()
  valid_manifest <- manifest[manifest$fit_success %in% TRUE & manifest$seed_id %in% rownames(param_mat), , drop = FALSE]
  valid_manifest <- valid_manifest[order(o2ipa_seed_number(valid_manifest$seed_id)), , drop = FALSE]
  counter <- 0L
  for (i in seq_len(nrow(valid_manifest))) {
    seed <- valid_manifest$seed_id[[i]]
    if (i %% 25L == 0L || i == 1L) message("Processing ", seed, " (", i, "/", nrow(valid_manifest), ")")
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- o2pr_run_params_from_vec(pvec, cfg)
    for (O2 in o2_grid) {
      counter <- counter + 1L
      rows[[counter]] <- fo2_dominant_attractor_one(seed, run_params, model_env, cfg, O2)
    }
  }
  attractors <- do.call(rbind, rows)
  attractors <- merge(attractors, manifest[, intersect(c("seed_id", "trajectory_regime", "mode_label", "objective", "delta_objective", "mean_terminal_ploidy", "mean_late_drop_amplitude"), names(manifest)), drop = FALSE], by = "seed_id", all.x = TRUE, sort = FALSE)
  fo2_write_tsv(attractors, file.path(out_dir, "tables", "fixed_o2_attractors_by_seed.tsv"))
  fo2_write_tsv(fo2_mode_seed_stack_table(attractors), file.path(out_dir, "tables", "fixed_o2_dominant_ploidy_all_seed_stack_mode1_mode2.tsv"))
  fo2_write_mode_comparison_tables(attractors, out_dir)
  fo2_write_spectral_gap_outputs(attractors, out_dir, generate_figures = generate_figures)

  regime_summary <- do.call(rbind, lapply(split(attractors, list(attractors$O2_pct, attractors$trajectory_regime), drop = TRUE), function(d) {
    data.frame(
      O2_pct = d$O2_pct[[1]],
      trajectory_regime = d$trajectory_regime[[1]],
      n_seed = length(unique(d$seed_id)),
      dominant_mean_ploidy_median = stats::median(d$dominant_mean_ploidy, na.rm = TRUE),
      dominant_mean_ploidy_q25 = as.numeric(stats::quantile(d$dominant_mean_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
      dominant_mean_ploidy_q75 = as.numeric(stats::quantile(d$dominant_mean_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
      fraction_N_le_25_median = stats::median(d$dominant_fraction_N_le_25, na.rm = TRUE),
      fraction_N_below_44_median = stats::median(d$dominant_fraction_N_below_44, na.rm = TRUE),
      selection_22_vs_44_median = stats::median(d$selection_22_vs_44, na.rm = TRUE),
      selection_44_vs_88_median = stats::median(d$selection_44_vs_88, na.rm = TRUE),
      spectral_gap_median = stats::median(d$spectral_gap, na.rm = TRUE),
      eigenvector_nonnegative_fraction = mean(d$eigenvector_nonnegative %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  regime_summary <- regime_summary[order(regime_summary$O2_pct, regime_summary$trajectory_regime), , drop = FALSE]
  fo2_write_tsv(regime_summary, file.path(out_dir, "tables", "fixed_o2_attractor_regime_summary.tsv"))

  tests <- fo2_attractor_regime_tests(attractors)
  fo2_write_tsv(tests, file.path(out_dir, "tables", "fixed_o2_attractor_regime_tests.tsv"))

  corr <- fo2_parameter_correlations(attractors, inputs$params_long)
  fo2_write_tsv(corr[order(corr$O2_pct, corr$metric, -corr$abs_spearman_rho), , drop = FALSE], file.path(out_dir, "tables", "parameter_attractor_correlations.tsv"))

  if (generate_figures) fo2_plot_outputs(attractors, out_dir)
  fo2_write_report(out_dir, run_dir, label_file, attractors, regime_summary, tests, corr)
  message("Completed fixed-O2 attractor analysis: ", out_dir)
}

if (sys.nframe() == 0L) {
  main()
}
