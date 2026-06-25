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

options(error = function() {
  traceback(2)
  quit(status = 1)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

mo2_mkdirs <- function(out_dir) {
  invisible(vapply(
    file.path(out_dir, c("tables", "figures", "logs", "report")),
    dir.create,
    logical(1),
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

mo2_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

mo2_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

mo2_parse_windows <- function(x) {
  txt <- o2ipa_as_chr(x, "0.5:2,1:5")
  parts <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  rows <- lapply(parts, function(p) {
    vals <- suppressWarnings(as.numeric(trimws(strsplit(p, ":", fixed = TRUE)[[1]])))
    if (length(vals) != 2L || any(!is.finite(vals))) {
      stop("Bad window specification: ", p, ". Use low:high, e.g. 0.5:2,1:5.")
    }
    low <- min(vals)
    high <- max(vals)
    low_label <- gsub("[^0-9]+", "p", format(low, scientific = FALSE, trim = TRUE))
    high_label <- gsub("[^0-9]+", "p", format(high, scientific = FALSE, trim = TRUE))
    data.frame(
      low = low,
      high = high,
      label = paste0("O2_", low_label, "_to_", high_label, "pct"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mo2_read_labels <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(data.frame(seed_id = character()))
  tab <- mo2_read_tsv(path)
  if (!"seed_id" %in% names(tab)) stop("label_file must contain seed_id: ", path)
  keep <- intersect(c(
    "seed_id", "trajectory_regime", "mode_label", "objective", "delta_objective",
    "mean_o2_day0_pct", "mean_time_o2_near_floor", "mean_terminal_ploidy",
    "mean_late_drop_amplitude"
  ), names(tab))
  tab[, keep, drop = FALSE]
}

mo2_read_o2_curve <- function(seed_id, seed_dir, o2_source = "eff") {
  if (is.na(seed_dir) || !nzchar(seed_dir) || !dir.exists(seed_dir)) {
    return(list(curve = data.frame(), status = "missing_seed_dir", source_file = NA_character_))
  }
  lag_path <- file.path(seed_dir, "viz", "o2_lag_timecourse.tsv")
  burden_path <- file.path(seed_dir, "viz", "predict_burden_vs_o2.tsv")
  path <- if (file.exists(lag_path)) lag_path else if (file.exists(burden_path)) burden_path else NA_character_
  if (is.na(path)) {
    return(list(curve = data.frame(), status = "missing_o2_table", source_file = NA_character_))
  }
  tab <- tryCatch(mo2_read_tsv(path), error = function(e) data.frame())
  if (!nrow(tab) || !all(c("cohort", "day") %in% names(tab))) {
    return(list(curve = data.frame(), status = "bad_o2_table", source_file = path))
  }
  source_col <- NA_character_
  if (identical(basename(path), "o2_lag_timecourse.tsv")) {
    preferred <- if (identical(o2_source, "target")) "o2_target_pct" else "o2_eff_pct"
    fallback <- if (identical(preferred, "o2_eff_pct")) "o2_target_pct" else "o2_eff_pct"
    source_col <- if (preferred %in% names(tab)) preferred else if (fallback %in% names(tab)) fallback else NA_character_
  } else if ("o2_pct" %in% names(tab)) {
    source_col <- "o2_pct"
  }
  if (is.na(source_col)) {
    return(list(curve = data.frame(), status = "missing_o2_column", source_file = path))
  }
  out <- data.frame(
    seed_id = seed_id,
    cohort = as.character(tab$cohort),
    day = suppressWarnings(as.numeric(tab$day)),
    o2_pct = suppressWarnings(as.numeric(tab[[source_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[out$cohort %in% c("2N", "4N") & is.finite(out$day) & is.finite(out$o2_pct), , drop = FALSE]
  if (!nrow(out)) {
    return(list(curve = data.frame(), status = "no_valid_o2_rows", source_file = path))
  }
  agg <- aggregate(out$o2_pct, by = list(seed_id = out$seed_id, cohort = out$cohort, day = out$day), FUN = mean, na.rm = TRUE)
  names(agg)[4] <- "o2_pct"
  agg$source_file <- normalizePath(path, mustWork = FALSE)
  agg$o2_source_column <- source_col
  list(curve = agg[order(agg$cohort, agg$day), , drop = FALSE], status = "ok", source_file = path)
}

mo2_interp_linear <- function(x0, x1, y0, y1, x) {
  if (!is.finite(x1 - x0) || abs(x1 - x0) < 1e-12) return(NA_real_)
  y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

mo2_segment_window_features <- function(day, value, low, high) {
  ok <- is.finite(day) & is.finite(value)
  day <- day[ok]
  value <- value[ok]
  if (length(day) < 2L) {
    return(c(time_in_window = NA_real_, first_enter = NA_real_, last_exit = NA_real_, auc_in_window = NA_real_, mean_o2_in_window = NA_real_))
  }
  ord <- order(day)
  day <- day[ord]
  value <- value[ord]
  total_time <- 0
  total_auc <- 0
  first_enter <- NA_real_
  last_exit <- NA_real_
  for (i in seq_len(length(day) - 1L)) {
    x0 <- day[[i]]
    x1 <- day[[i + 1L]]
    y0 <- value[[i]]
    y1 <- value[[i + 1L]]
    if (!is.finite(x0) || !is.finite(x1) || !is.finite(y0) || !is.finite(y1) || x1 <= x0) next
    pts <- c(x0, x1)
    if (abs(y1 - y0) > 1e-12) {
      for (thr in c(low, high)) {
        if ((thr - y0) * (thr - y1) <= 0) {
          xc <- x0 + (thr - y0) * (x1 - x0) / (y1 - y0)
          if (is.finite(xc) && xc >= x0 && xc <= x1) pts <- c(pts, xc)
        }
      }
    }
    pts <- sort(unique(pts))
    if (length(pts) < 2L) next
    for (j in seq_len(length(pts) - 1L)) {
      a <- pts[[j]]
      b <- pts[[j + 1L]]
      if (b <= a) next
      mid <- (a + b) / 2
      ymid <- mo2_interp_linear(x0, x1, y0, y1, mid)
      if (is.finite(ymid) && ymid >= low && ymid <= high) {
        ya <- mo2_interp_linear(x0, x1, y0, y1, a)
        yb <- mo2_interp_linear(x0, x1, y0, y1, b)
        total_time <- total_time + (b - a)
        total_auc <- total_auc + (b - a) * (ya + yb) / 2
        if (!is.finite(first_enter)) first_enter <- a
        last_exit <- b
      }
    }
  }
  c(
    time_in_window = if (total_time > 0) total_time else 0,
    first_enter = first_enter,
    last_exit = last_exit,
    auc_in_window = if (total_time > 0) total_auc else 0,
    mean_o2_in_window = if (total_time > 0) total_auc / total_time else NA_real_
  )
}

mo2_summarize_cohort <- function(seed_id, cohort, d, windows) {
  d <- d[order(d$day), , drop = FALSE]
  out <- data.frame(
    seed_id = seed_id,
    cohort = cohort,
    n_timepoints = nrow(d),
    day_min = if (nrow(d)) min(d$day, na.rm = TRUE) else NA_real_,
    day_max = if (nrow(d)) max(d$day, na.rm = TRUE) else NA_real_,
    o2_day0_pct = if (nrow(d)) d$o2_pct[which.min(abs(d$day - 0))] else NA_real_,
    o2_auc_0_to_end = o2ipa_auc(d$day, d$o2_pct),
    source_file = if (nrow(d)) d$source_file[[1]] else NA_character_,
    o2_source_column = if (nrow(d)) d$o2_source_column[[1]] else NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(windows))) {
    w <- windows[i, , drop = FALSE]
    feat <- mo2_segment_window_features(d$day, d$o2_pct, w$low, w$high)
    lab <- w$label[[1]]
    out[[paste0("time_in_", lab)]] <- feat[["time_in_window"]]
    out[[paste0("first_enter_", lab)]] <- feat[["first_enter"]]
    out[[paste0("last_exit_", lab)]] <- feat[["last_exit"]]
    out[[paste0("AUC_", lab)]] <- feat[["auc_in_window"]]
    out[[paste0("mean_", lab)]] <- feat[["mean_o2_in_window"]]
  }
  out
}

mo2_wide_metric <- function(tab, seed_ids, cohort, metric) {
  d <- tab[tab$cohort == cohort, , drop = FALSE]
  vals <- d[[metric]]
  names(vals) <- d$seed_id
  as.numeric(vals[seed_ids])
}

mo2_make_seed_features <- function(manifest, cohort_features, labels) {
  seed_features <- manifest
  if (nrow(labels)) {
    seed_features <- merge(seed_features, labels, by = "seed_id", all.x = TRUE, sort = FALSE)
  }
  metrics <- setdiff(names(cohort_features), c("seed_id", "cohort", "source_file", "o2_source_column"))
  for (cohort in c("2N", "4N")) {
    for (metric in metrics) {
      seed_features[[paste0(cohort, "__", metric)]] <- mo2_wide_metric(cohort_features, seed_features$seed_id, cohort, metric)
    }
  }
  for (metric in metrics) {
    x <- seed_features[[paste0("2N__", metric)]]
    y <- seed_features[[paste0("4N__", metric)]]
    if (is.numeric(x) && is.numeric(y)) {
      seed_features[[paste0("mean__", metric)]] <- rowMeans(cbind(x, y), na.rm = TRUE)
    }
  }
  seed_features
}

mo2_wilcox_row <- function(tab, value_col, group_col = "trajectory_regime",
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

mo2_feature_tests <- function(seed_features) {
  feature_cols <- grep("(^2N__|^4N__|^mean__)(time_in_|first_enter_|last_exit_|AUC_|mean_O2_)", names(seed_features), value = TRUE)
  if (!length(feature_cols)) return(data.frame())
  out <- do.call(rbind, lapply(feature_cols, function(v) mo2_wilcox_row(seed_features, v, feature_name = v)))
  out$BH_adjusted_p_value <- stats::p.adjust(out$wilcox_p_value, method = "BH")
  out[order(out$BH_adjusted_p_value, out$feature), , drop = FALSE]
}

mo2_parameter_correlations <- function(seed_features, params_long) {
  if (!nrow(seed_features) || !nrow(params_long)) return(data.frame())
  params <- params_long
  params$transformed_value <- mapply(o2ipa_transform_parameter_value, params$parameter, params$value)
  params_wide <- o2ipa_params_wide(params, value_col = "transformed_value")
  params_wide$seed_id <- rownames(params_wide)
  feature_cols <- grep("(^2N__|^4N__|^mean__)(time_in_|first_enter_|last_exit_|AUC_|mean_O2_)", names(seed_features), value = TRUE)
  merged <- merge(seed_features[, c("seed_id", "trajectory_regime", feature_cols), drop = FALSE], params_wide, by = "seed_id", all.x = TRUE)
  scopes <- list(
    all_seeds = merged,
    mode1_mode2_only = merged[merged$trajectory_regime %in% c("mode1_ploidy_stable", "mode2_second_ploidy_collapse"), , drop = FALSE]
  )
  rows <- list()
  counter <- 0L
  for (scope_name in names(scopes)) {
    dscope <- scopes[[scope_name]]
    for (feature in feature_cols) {
      y <- dscope[[feature]]
      for (param in o2ipa_target_params()) {
        x <- dscope[[param]]
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < 5L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) next
        ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          correlation_scope = scope_name,
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
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$BH_adjusted_p_value <- ave(out$p_value, out$correlation_scope, FUN = function(p) stats::p.adjust(p, method = "BH"))
  out[order(out$correlation_scope, -out$abs_spearman_rho, out$BH_adjusted_p_value), , drop = FALSE]
}

mo2_duration_summary_by_regime <- function(seed_features, windows) {
  if (!nrow(seed_features) || !"trajectory_regime" %in% names(seed_features)) return(data.frame())
  summary_cols <- paste0("mean__time_in_", windows$label)
  summary_cols <- summary_cols[summary_cols %in% names(seed_features)]
  if (!length(summary_cols)) return(data.frame())
  rows <- list()
  counter <- 0L
  for (col in summary_cols) {
    for (reg in c("mode1_ploidy_stable", "mode2_second_ploidy_collapse", "ambiguous")) {
      x <- seed_features[[col]][seed_features$trajectory_regime == reg]
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        feature = col,
        trajectory_regime = reg,
        n_seed = sum(is.finite(x)),
        median = stats::median(x, na.rm = TRUE),
        mean = mean(x, na.rm = TRUE),
        q25 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)),
        q75 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

mo2_plot_duration_boxplots <- function(seed_features, fig_dir, windows) {
  if (!nrow(seed_features) || !"trajectory_regime" %in% names(seed_features)) return(invisible(FALSE))
  cols <- paste0("mean__time_in_", windows$label)
  cols <- cols[cols %in% names(seed_features)]
  if (!length(cols)) return(invisible(FALSE))
  grDevices::pdf(file.path(fig_dir, "medium_o2_duration_by_regime.pdf"), width = 8, height = 4 + 2 * length(cols))
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(length(cols), 1), mar = c(8, 5, 3, 1))
  for (col in cols) {
    d <- seed_features[seed_features$trajectory_regime %in% c("mode1_ploidy_stable", "mode2_second_ploidy_collapse") & is.finite(seed_features[[col]]), , drop = FALSE]
    if (!nrow(d)) next
    boxplot(d[[col]] ~ d$trajectory_regime, las = 2, ylab = "days", main = col)
    stripchart(d[[col]] ~ d$trajectory_regime, vertical = TRUE, method = "jitter", pch = 16, cex = 0.4, col = grDevices::adjustcolor("black", alpha.f = 0.35), add = TRUE)
  }
  invisible(TRUE)
}

mo2_report <- function(out_dir, args, windows, manifest, cohort_features, seed_features, tests, correlations) {
  regime_counts <- if ("trajectory_regime" %in% names(seed_features)) {
    as.data.frame(table(seed_features$trajectory_regime), stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  names(regime_counts) <- c("trajectory_regime", "n_seed")
  top_tests <- tests[is.finite(tests$BH_adjusted_p_value), , drop = FALSE]
  top_tests <- head(top_tests, 20L)
  top_corr <- correlations[is.finite(correlations$abs_spearman_rho), , drop = FALSE]
  top_corr <- head(top_corr, 30L)
  summary_cols <- paste0("mean__time_in_", windows$label)
  summary_cols <- summary_cols[summary_cols %in% names(seed_features)]
  window_summary <- do.call(rbind, lapply(summary_cols, function(col) {
    do.call(rbind, lapply(c("mode1_ploidy_stable", "mode2_second_ploidy_collapse"), function(reg) {
      x <- seed_features[[col]][seed_features$trajectory_regime == reg]
      data.frame(
        feature = col,
        trajectory_regime = reg,
        n_seed = sum(is.finite(x)),
        median = stats::median(x, na.rm = TRUE),
        mean = mean(x, na.rm = TRUE),
        q25 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)),
        q75 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)),
        stringsAsFactors = FALSE
      )
    }))
  }))
  lines <- c(
    "# In vivo medium-O2 window analysis",
    "",
    paste0("- run_dir: `", o2ipa_as_chr(args$run_dir, ""), "`"),
    paste0("- label_file: `", o2ipa_as_chr(args$label_file, ""), "`"),
    paste0("- O2 source: `", o2ipa_as_chr(args$o2_source, "eff"), "`"),
    paste0("- windows: `", paste(paste0(windows$low, "-", windows$high), collapse = ", "), "`"),
    paste0("- seeds discovered: ", nrow(manifest)),
    paste0("- cohort feature rows: ", nrow(cohort_features)),
    "",
    "## Regime counts",
    "",
    paste(capture.output(print(regime_counts, row.names = FALSE)), collapse = "\n"),
    "",
    "## Medium-O2 duration summary",
    "",
    paste(capture.output(print(window_summary, row.names = FALSE)), collapse = "\n"),
    "",
    "## Top mode1-vs-mode2 medium-O2 tests",
    "",
    paste(capture.output(print(top_tests, row.names = FALSE)), collapse = "\n"),
    "",
    "## Top parameter correlations",
    "",
    paste(capture.output(print(top_corr, row.names = FALSE)), collapse = "\n"),
    "",
    "## Notes",
    "",
    "- Duration and AUC are computed by linear interpolation along each observed O2 timecourse segment.",
    "- `time_in_*` is the exact interpolated time in days spent inside the requested O2 range.",
    "- `AUC_*` is the O2 percent-days accumulated only while the trajectory is inside that range.",
    "- Parameter correlations use the same transformed parameter scale as prior process-fingerprint analyses."
  )
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}

main <- function() {
  args <- o2ipa_parse_args()
  repo_root <- o2ipa_repo_root(SCRIPT_DIR)
  run_dir <- normalizePath(o2ipa_as_chr(
    args$run_dir,
    file.path(repo_root, "oxygen", "results", "fit_invivo_O2_buffering_500seed")
  ), mustWork = FALSE)
  label_file <- normalizePath(o2ipa_as_chr(
    args$label_file,
    file.path(repo_root, "oxygen", "results", "analysis", "invivo_o2_ploidy_event_coupling_500seed", "tables", "seed_event_summary.tsv")
  ), mustWork = FALSE)
  out_dir <- normalizePath(o2ipa_as_chr(
    args$out_dir,
    file.path(repo_root, "oxygen", "results", "analysis", "invivo_medium_o2_windows_500seed")
  ), mustWork = FALSE)
  windows <- mo2_parse_windows(args$windows)
  max_seeds <- o2ipa_as_int(args$max_seeds, 0L)
  o2_source <- match.arg(o2ipa_as_chr(args$o2_source, "eff"), c("eff", "target"))
  generate_figures <- o2ipa_as_bool(args$generate_figures, TRUE)

  mo2_mkdirs(out_dir)
  log_file <- file.path(out_dir, "logs", "analyze_invivo_medium_o2_windows.log")
  sink(log_file, split = TRUE)
  on.exit(sink(), add = TRUE)
  message("run_dir: ", run_dir)
  message("label_file: ", label_file)
  message("out_dir: ", out_dir)

  seed_inputs <- o2ipa_collect_seed_inputs(run_dir, objective_source = "auto")
  manifest <- seed_inputs$manifest
  params_long <- seed_inputs$params_long
  if (max_seeds > 0L && nrow(manifest) > max_seeds) {
    manifest <- manifest[seq_len(max_seeds), , drop = FALSE]
    params_long <- params_long[params_long$seed_id %in% manifest$seed_id, , drop = FALSE]
  }
  labels <- mo2_read_labels(label_file)
  if (nrow(labels)) {
    manifest <- merge(manifest, labels, by = "seed_id", all.x = TRUE, sort = FALSE)
  }

  curve_rows <- list()
  feature_rows <- list()
  status_rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    seed_id <- manifest$seed_id[[i]]
    seed_dir <- manifest$seed_dir[[i]]
    if (i == 1L || i %% 25L == 0L) {
      message("Processing ", seed_id, " (", i, "/", nrow(manifest), ")")
    }
    read <- mo2_read_o2_curve(seed_id, seed_dir, o2_source = o2_source)
    status_rows[[i]] <- data.frame(
      seed_id = seed_id,
      o2_read_status = read$status,
      o2_source_file = read$source_file,
      stringsAsFactors = FALSE
    )
    if (!nrow(read$curve)) next
    curve_rows[[seed_id]] <- read$curve
    for (cohort in c("2N", "4N")) {
      d <- read$curve[read$curve$cohort == cohort, , drop = FALSE]
      if (nrow(d)) {
        feature_rows[[paste(seed_id, cohort, sep = "||")]] <- mo2_summarize_cohort(seed_id, cohort, d, windows)
      }
    }
  }
  status <- do.call(rbind, status_rows)
  o2_curves <- if (length(curve_rows)) do.call(rbind, curve_rows) else data.frame()
  cohort_features <- if (length(feature_rows)) do.call(rbind, feature_rows) else data.frame()

  seed_features <- mo2_make_seed_features(
    manifest[, setdiff(names(manifest), intersect(names(labels), c("trajectory_regime", "mode_label", "objective", "delta_objective", "mean_o2_day0_pct", "mean_time_o2_near_floor", "mean_terminal_ploidy", "mean_late_drop_amplitude"))), drop = FALSE],
    cohort_features,
    labels
  )
  seed_features <- merge(seed_features, status, by = "seed_id", all.x = TRUE, sort = FALSE)

  tests <- mo2_feature_tests(seed_features)
  correlations <- mo2_parameter_correlations(seed_features, params_long)
  duration_summary <- mo2_duration_summary_by_regime(seed_features, windows)

  mo2_write_tsv(data.frame(
    argument = c("run_dir", "label_file", "out_dir", "windows", "o2_source", "max_seeds", "generate_figures"),
    value = c(run_dir, label_file, out_dir, o2ipa_as_chr(args$windows, "0.5:2,1:5"), o2_source, as.character(max_seeds), as.character(generate_figures)),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  mo2_write_tsv(windows, file.path(out_dir, "tables", "medium_o2_windows.tsv"))
  mo2_write_tsv(status, file.path(out_dir, "tables", "o2_read_status.tsv"))
  mo2_write_tsv(cohort_features, file.path(out_dir, "tables", "medium_o2_features_by_cohort.tsv"))
  mo2_write_tsv(seed_features, file.path(out_dir, "tables", "medium_o2_features_by_seed.tsv"))
  mo2_write_tsv(duration_summary, file.path(out_dir, "tables", "medium_o2_duration_summary_by_regime.tsv"))
  mo2_write_tsv(tests, file.path(out_dir, "tables", "medium_o2_feature_regime_tests.tsv"))
  mo2_write_tsv(correlations, file.path(out_dir, "tables", "medium_o2_parameter_correlations.tsv"))
  mo2_write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))
  if (nrow(o2_curves)) {
    mo2_write_tsv(o2_curves, file.path(out_dir, "tables", "o2_timecourses_long.tsv"))
  }
  if (isTRUE(generate_figures)) {
    mo2_plot_duration_boxplots(seed_features, file.path(out_dir, "figures"), windows)
  }
  mo2_report(out_dir, args, windows, manifest, cohort_features, seed_features, tests, correlations)
  message("Completed medium-O2 window analysis: ", out_dir)
}

if (identical(environment(), globalenv())) {
  main()
}
