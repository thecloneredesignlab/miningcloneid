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

`%||%` <- o2ipa_null_coalesce

mo2_mkdirs <- function(out_dir) {
  invisible(vapply(
    file.path(out_dir, c("tables", "logs")),
    dir.create,
    logical(1),
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

mo2_read_tsv <- o2ipa_read_tsv
mo2_write_tsv <- o2ipa_write_tsv

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

mo2_interp_linear <- function(x0, x1, y0, y1, x) {
  if (!is.finite(x1 - x0) || abs(x1 - x0) < 1e-12) return(NA_real_)
  y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

mo2_auc <- o2ipa_auc

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
    o2_auc_0_to_end = mo2_auc(d$day, d$o2_pct),
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

mo2_safe_read <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 1L) return(data.frame())
  tryCatch(mo2_read_tsv(path), error = function(e) data.frame())
}

mo2_write_analysis_manifest <- function(out_dir) {
  files <- list.files(file.path(out_dir, "tables"), pattern = "[.]tsv$", full.names = FALSE)
  rows <- lapply(files, function(file) {
    tab <- mo2_safe_read(file.path(out_dir, "tables", file))
    data.frame(
      artifact = tools::file_path_sans_ext(file), relative_path = file.path("tables", file),
      role = "analysis_table", rows = nrow(tab), columns = ncol(tab),
      path = normalizePath(file.path(out_dir, "tables", file), mustWork = TRUE),
      exists = TRUE, stringsAsFactors = FALSE
    )
  })
  mo2_write_tsv(do.call(rbind, rows), file.path(out_dir, "analysis_manifest.tsv"))
}

run_medium_o2_window_analysis <- function(argv = o2ipa_parse_args()) {
  simulation_dir <- o2ipa_as_chr(argv$simulation_dir)
  if (!nzchar(simulation_dir) || !dir.exists(simulation_dir)) stop("Missing or invalid --simulation_dir.")
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  out_dir <- o2ipa_as_chr(argv$analysis_dir %||% argv$out_dir)
  if (!nzchar(out_dir)) stop("Missing --analysis_dir (or --out_dir).")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  required <- c(
    "event_seed_manifest.tsv", "event_parameter_values_long.tsv",
    "event_o2_timecourses.tsv", "event_input_status.tsv", "event_input_labels.tsv"
  )
  if (!file.exists(file.path(simulation_dir, "simulation_manifest.tsv"))) stop("Missing O2/ploidy event simulation manifest.")
  missing <- required[!file.exists(file.path(simulation_dir, "tables", required))]
  if (length(missing)) stop("Missing medium-O2 simulation inputs: ", paste(missing, collapse = ", "))
  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

  windows <- mo2_parse_windows(argv$windows)
  o2_source <- match.arg(o2ipa_as_chr(argv$o2_source, "eff"), c("eff", "target"))
  manifest <- mo2_safe_read(file.path(simulation_dir, "tables", "event_seed_manifest.tsv"))
  params_long <- mo2_safe_read(file.path(simulation_dir, "tables", "event_parameter_values_long.tsv"))
  labels <- mo2_safe_read(file.path(simulation_dir, "tables", "event_input_labels.tsv"))
  status <- mo2_safe_read(file.path(simulation_dir, "tables", "event_input_status.tsv"))
  curves <- mo2_safe_read(file.path(simulation_dir, "tables", "event_o2_timecourses.tsv"))
  if (identical(o2_source, "target") && "o2_target_pct" %in% names(curves) && any(is.finite(curves$o2_target_pct))) {
    curves$o2_pct <- ifelse(is.finite(curves$o2_target_pct), curves$o2_target_pct, curves$o2_pct)
    curves$o2_source_column <- "o2_target_pct"
  } else {
    curves$o2_source_column <- "o2_pct"
  }

  feature_rows <- list()
  for (idx in split(seq_len(nrow(curves)), paste(curves$seed_id, curves$cohort, sep = "||"))) {
    d <- curves[idx, , drop = FALSE]
    feature_rows[[length(feature_rows) + 1L]] <- mo2_summarize_cohort(d$seed_id[[1]], d$cohort[[1]], d, windows)
  }
  cohort_features <- if (length(feature_rows)) do.call(rbind, feature_rows) else data.frame()
  label_cols <- intersect(names(labels), c("trajectory_regime", "mode_label", "objective", "delta_objective", "mean_o2_day0_pct", "mean_time_o2_near_floor", "mean_terminal_ploidy", "mean_late_drop_amplitude"))
  manifest_base <- manifest[, setdiff(names(manifest), label_cols), drop = FALSE]
  seed_features <- mo2_make_seed_features(manifest_base, cohort_features, labels)
  status_keep <- intersect(c("seed_id", "o2_available", "o2_source_file"), names(status))
  if (length(status_keep) > 1L) seed_features <- merge(seed_features, status[, status_keep, drop = FALSE], by = "seed_id", all.x = TRUE, sort = FALSE)
  tests <- mo2_feature_tests(seed_features)
  correlations <- mo2_parameter_correlations(seed_features, params_long)
  duration_summary <- mo2_duration_summary_by_regime(seed_features, windows)

  mo2_write_tsv(data.frame(
    argument = c("simulation_dir", "windows", "o2_source"),
    value = c(simulation_dir, o2ipa_as_chr(argv$windows, "0.5:2,1:5"), o2_source),
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
  mo2_write_tsv(curves, file.path(out_dir, "tables", "o2_timecourses_long.tsv"))
  mo2_write_analysis_manifest(out_dir)
  message("Completed medium-O2 window analysis: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_medium_o2_window_analysis()
