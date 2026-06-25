#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!grepl("^--", arg)) {
      i <- i + 1L
      next
    }
    kv <- sub("^--", "", arg)
    eq <- regexpr("=", kv, fixed = TRUE)
    if (eq > 0L) {
      out[[substr(kv, 1L, eq - 1L)]] <- substr(kv, eq + 1L, nchar(kv))
      i <- i + 1L
    } else {
      key <- kv
      if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      } else {
        out[[key]] <- TRUE
        i <- i + 1L
      }
    }
  }
  out
}

as_chr <- function(x, default = "") {
  val <- x %||% default
  val <- as.character(val[[1]])
  if (!nzchar(val)) default else val
}

as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric((x %||% default)[[1]]))
  if (is.finite(val)) val else default
}

as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer((x %||% default)[[1]]))
  if (!is.na(val)) val else default
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  if (is.logical(x[[1]])) return(isTRUE(x[[1]]))
  tolower(trimws(as.character(x[[1]]))) %in% c("1", "true", "t", "yes", "y", "on")
}

as_num_pair <- function(x, default) {
  txt <- as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2L) default else vals[seq_len(2L)]
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

seed_number <- function(seed_id) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

auc_window <- function(day, value, from, to) {
  ok <- is.finite(day) & is.finite(value) & day >= from & day <= to
  x <- day[ok]
  y <- value[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

interp_value <- function(day, value, at_day) {
  ok <- is.finite(day) & is.finite(value)
  x <- day[ok]
  y <- value[ok]
  if (!is.finite(at_day) || length(x) < 2L || at_day < min(x) || at_day > max(x)) return(NA_real_)
  ord <- order(x)
  as.numeric(stats::approx(x[ord], y[ord], xout = at_day, ties = "ordered", rule = 1)$y)
}

crossing_time <- function(day, value, threshold, direction = "down", after_day = -Inf) {
  ok <- is.finite(day) & is.finite(value) & day >= after_day
  x <- day[ok]
  y <- value[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  hit <- if (identical(direction, "down")) which(y <= threshold) else which(y >= threshold)
  if (!length(hit)) return(NA_real_)
  i <- hit[[1]]
  if (i == 1L) return(x[[1]])
  x0 <- x[[i - 1L]]
  x1 <- x[[i]]
  y0 <- y[[i - 1L]]
  y1 <- y[[i]]
  if (!is.finite(y1 - y0) || abs(y1 - y0) < 1e-12) return(x1)
  x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)
}

grid_curve <- function(day, value, by = 1) {
  ok <- is.finite(day) & is.finite(value)
  x <- day[ok]
  y <- value[ok]
  if (length(x) < 2L) return(data.frame(day = numeric(), value = numeric()))
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  x2 <- seq(ceiling(min(x)), floor(max(x)), by = by)
  data.frame(day = x2, value = as.numeric(stats::approx(x, y, xout = x2, ties = "ordered", rule = 1)$y))
}

slope_event <- function(day, value, from, to, lag_days = 20, direction = "down", threshold = NA_real_) {
  g <- grid_curve(day, value, by = 1)
  if (nrow(g) <= lag_days) {
    return(list(time = NA_real_, slope = NA_real_, value = NA_real_, status = "insufficient_grid"))
  }
  slopes <- (g$value[(lag_days + 1L):nrow(g)] - g$value[seq_len(nrow(g) - lag_days)]) / lag_days
  slope_day <- g$day[(lag_days + 1L):nrow(g)]
  ok <- is.finite(slopes) & slope_day >= from & slope_day <= to
  if (!any(ok)) return(list(time = NA_real_, slope = NA_real_, value = NA_real_, status = "no_window_values"))
  slopes2 <- slopes[ok]
  days2 <- slope_day[ok]
  if (identical(direction, "down")) {
    idxs <- if (is.finite(threshold)) which(slopes2 <= threshold) else integer()
    if (length(idxs)) idx <- idxs[[1]] else idx <- which.min(slopes2)
  } else {
    idxs <- if (is.finite(threshold)) which(slopes2 >= threshold) else integer()
    if (length(idxs)) idx <- idxs[[1]] else idx <- which.max(slopes2)
  }
  tt <- days2[[idx]]
  list(time = tt, slope = slopes2[[idx]], value = interp_value(day, value, tt), status = if (length(idxs)) "threshold_event" else "largest_slope")
}

metric_map <- function(path) {
  if (!file.exists(path)) return(list())
  tab <- tryCatch(read_tsv(path), error = function(e) data.frame())
  if (!all(c("metric", "value") %in% names(tab))) return(list())
  vals <- as.list(tab$value)
  names(vals) <- tab$metric
  vals
}

map_num <- function(map, keys) {
  for (key in keys) {
    val <- safe_num(map[[key]])
    if (length(val) && is.finite(val[[1]])) return(val[[1]])
  }
  NA_real_
}

map_chr <- function(map, keys) {
  for (key in keys) {
    val <- map[[key]]
    if (!is.null(val) && length(val) && !is.na(val[[1]]) && nzchar(as.character(val[[1]]))) return(as.character(val[[1]]))
  }
  NA_character_
}

read_params <- function(path) {
  if (!file.exists(path)) return(data.frame(parameter = character(), value = numeric()))
  tab <- tryCatch(read_tsv(path), error = function(e) data.frame())
  if (!all(c("parameter", "value") %in% names(tab))) return(data.frame(parameter = character(), value = numeric()))
  data.frame(parameter = trimws(as.character(tab$parameter)), value = safe_num(tab$value), stringsAsFactors = FALSE)
}

param_value <- function(params, name, default = NA_real_) {
  hit <- params$value[match(name, params$parameter)]
  if (length(hit) && is.finite(hit[[1]])) hit[[1]] else default
}

read_seed_manifest_row <- function(seed_id, seed_dir) {
  m <- metric_map(file.path(seed_dir, "fit_summary.tsv"))
  burden_raw <- map_num(m, c("objective_burden_neg2loglik_raw"))
  ploidy_raw <- map_num(m, c("objective_ploidy_neg2loglik_raw"))
  raw_sum <- if (is.finite(burden_raw) && is.finite(ploidy_raw)) burden_raw + ploidy_raw else NA_real_
  objective <- raw_sum
  objective_source <- "raw_likelihood_sum"
  if (!is.finite(objective)) {
    objective <- map_num(m, c("objective_data", "objective_total", "objective", "optimizer_local_objective", "optimizer_deoptim_objective"))
    objective_source <- "available_objective"
  }
  data.frame(
    seed_id = seed_id,
    seed_dir = normalizePath(seed_dir, mustWork = FALSE),
    objective = objective,
    objective_source = objective_source,
    optimizer_local_objective = map_num(m, c("optimizer_local_objective")),
    optimizer_deoptim_objective = map_num(m, c("optimizer_deoptim_objective")),
    objective_burden_raw = burden_raw,
    objective_ploidy_raw = ploidy_raw,
    optimizer_iter_completed = map_num(m, c("optimizer_iter_completed", "deoptim_iter_completed", "itermax")),
    optimizer_stop_reason = map_chr(m, c("deoptim_stop_reason", "optimizer_stop_reason", "optimizer_local_convergence")),
    stringsAsFactors = FALSE
  )
}

read_ploidy_curve <- function(seed_id, seed_dir, n_unit = 22) {
  path <- file.path(seed_dir, "viz", "predict_ploidy_weighted_mean_0_1000day.tsv")
  if (!file.exists(path)) path <- file.path(seed_dir, "viz", "ploidy_weighted_mean_timecourse.tsv")
  if (!file.exists(path)) return(data.frame())
  tab <- tryCatch(read_tsv(path), error = function(e) data.frame())
  vcol <- intersect(c("weighted_mean_N", "weighted_mean_endpoint", "ploidy_value", "weighted_mean_ploidy"), names(tab))
  if (!length(vcol) || !all(c("cohort", "day") %in% names(tab))) return(data.frame())
  val <- safe_num(tab[[vcol[[1]]]])
  ploidy <- if (stats::median(val, na.rm = TRUE) > 10) val / n_unit else val
  d <- data.frame(
    seed_id = seed_id,
    cohort = as.character(tab$cohort),
    day = safe_num(tab$day),
    ploidy = ploidy,
    source_file = normalizePath(path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  d <- d[d$cohort %in% c("2N", "4N") & is.finite(d$day) & is.finite(d$ploidy), , drop = FALSE]
  aggregate(d$ploidy, by = list(seed_id = d$seed_id, cohort = d$cohort, day = d$day, source_file = d$source_file), FUN = mean)
}

read_o2_curve <- function(seed_id, seed_dir) {
  candidates <- c(
    file.path(seed_dir, "viz", "predict_burden_0_1000day.tsv"),
    file.path(seed_dir, "viz", "o2_lag_timecourse.tsv"),
    file.path(seed_dir, "viz", "predict_burden_vs_o2.tsv")
  )
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) return(data.frame())
  path <- hits[[1]]
  tab <- tryCatch(read_tsv(path), error = function(e) data.frame())
  if (!all(c("cohort", "day") %in% names(tab))) return(data.frame())
  o2_col <- intersect(c("pred_o2_pct", "o2_eff_pct", "o2_pct"), names(tab))
  target_col <- intersect(c("pred_o2_target_pct", "o2_target_pct"), names(tab))
  burden_col <- intersect(c("pred_burden", "pred_burden_volume_mm3", "burden_mm3", "burden_total"), names(tab))
  if (!length(o2_col)) return(data.frame())
  d <- data.frame(
    seed_id = seed_id,
    cohort = as.character(tab$cohort),
    day = safe_num(tab$day),
    o2_pct = safe_num(tab[[o2_col[[1]]]]),
    o2_target_pct = if (length(target_col)) safe_num(tab[[target_col[[1]]]]) else NA_real_,
    burden = if (length(burden_col)) safe_num(tab[[burden_col[[1]]]]) else NA_real_,
    source_file = normalizePath(path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  d <- d[d$cohort %in% c("2N", "4N") & is.finite(d$day) & is.finite(d$o2_pct), , drop = FALSE]
  if (!nrow(d)) return(d)
  rows <- lapply(split(seq_len(nrow(d)), paste(d$cohort, d$day, sep = "\r")), function(idx) {
    x <- d[idx, , drop = FALSE]
    data.frame(
      seed_id = seed_id,
      cohort = x$cohort[[1]],
      day = x$day[[1]],
      o2_pct = mean(x$o2_pct, na.rm = TRUE),
      o2_target_pct = mean(x$o2_target_pct, na.rm = TRUE),
      burden = mean(x$burden, na.rm = TRUE),
      source_file = x$source_file[[1]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_ploidy_one <- function(d, early_window, plateau_window, terminal_window, late_after_day,
                                 slope_lag_days, early_slope_threshold, late_slope_threshold) {
  day <- d$day
  y <- d$ploidy
  early <- y[day >= early_window[[1]] & day <= early_window[[2]]]
  plateau <- y[day >= plateau_window[[1]] & day <= plateau_window[[2]]]
  terminal <- y[day >= terminal_window[[1]] & day <= terminal_window[[2]]]
  e_down <- slope_event(day, y, early_window[[1]], early_window[[2]], slope_lag_days, "down", early_slope_threshold)
  e_up <- slope_event(day, y, early_window[[1]], early_window[[2]], slope_lag_days, "up", abs(early_slope_threshold))
  l_down <- slope_event(day, y, late_after_day, max(day, na.rm = TRUE), slope_lag_days, "down", late_slope_threshold)
  day0 <- interp_value(day, y, 0)
  data.frame(
    seed_id = d$seed_id[[1]],
    cohort = d$cohort[[1]],
    ploidy_day0 = day0,
    ploidy_day100 = interp_value(day, y, 100),
    ploidy_day200 = interp_value(day, y, 200),
    ploidy_day500 = interp_value(day, y, 500),
    ploidy_day1000 = interp_value(day, y, 1000),
    early_min_ploidy = if (length(early)) min(early, na.rm = TRUE) else NA_real_,
    early_max_ploidy = if (length(early)) max(early, na.rm = TRUE) else NA_real_,
    early_drop_amplitude = if (length(early) && is.finite(day0)) max(0, day0 - min(early, na.rm = TRUE)) else NA_real_,
    early_rise_amplitude = if (length(early) && is.finite(day0)) max(0, max(early, na.rm = TRUE) - day0) else NA_real_,
    time_early_largest_negative_slope = e_down$time,
    early_largest_negative_slope = e_down$slope,
    ploidy_at_early_largest_negative_slope = e_down$value,
    time_early_largest_positive_slope = e_up$time,
    early_largest_positive_slope = e_up$slope,
    ploidy_at_early_largest_positive_slope = e_up$value,
    plateau_reference_ploidy = if (length(plateau)) as.numeric(stats::quantile(plateau, 0.9, na.rm = TRUE, names = FALSE)) else NA_real_,
    terminal_median_ploidy = if (length(terminal)) stats::median(terminal, na.rm = TRUE) else NA_real_,
    terminal_min_ploidy = if (length(terminal)) min(terminal, na.rm = TRUE) else NA_real_,
    late_drop_amplitude = if (length(plateau) && length(terminal)) as.numeric(stats::quantile(plateau, 0.9, na.rm = TRUE, names = FALSE)) - stats::median(terminal, na.rm = TRUE) else NA_real_,
    time_late_largest_negative_slope = l_down$time,
    late_largest_negative_slope = l_down$slope,
    ploidy_at_late_largest_negative_slope = l_down$value,
    time_crossing_ploidy_2p0_down_after_late_start = crossing_time(day, y, 2.0, "down", late_after_day),
    time_crossing_ploidy_1p8_down_after_late_start = crossing_time(day, y, 1.8, "down", late_after_day),
    time_crossing_ploidy_1p5_down_after_late_start = crossing_time(day, y, 1.5, "down", late_after_day),
    time_crossing_ploidy_1p3_down_after_late_start = crossing_time(day, y, 1.3, "down", late_after_day),
    trajectory_auc_0_1000 = auc_window(day, y, 0, 1000),
    stringsAsFactors = FALSE
  )
}

summarize_o2_one <- function(d, o2_min, near_floor_abs, late_after_day) {
  day <- d$day
  y <- pmax(d$o2_pct, 0)
  floor_threshold <- max(o2_min + near_floor_abs, near_floor_abs, na.rm = TRUE)
  data.frame(
    seed_id = d$seed_id[[1]],
    cohort = d$cohort[[1]],
    o2_min_parameter = o2_min,
    o2_near_floor_threshold_pct = floor_threshold,
    o2_day0_pct = interp_value(day, y, 0),
    o2_day100_pct = interp_value(day, y, 100),
    o2_day200_pct = interp_value(day, y, 200),
    o2_day500_pct = interp_value(day, y, 500),
    o2_day1000_pct = interp_value(day, y, 1000),
    o2_minimum_0_1000_pct = min(y, na.rm = TRUE),
    o2_auc_0_100 = auc_window(day, y, 0, 100),
    o2_auc_0_200 = auc_window(day, y, 0, 200),
    o2_auc_200_700 = auc_window(day, y, 200, 700),
    o2_auc_700_1000 = auc_window(day, y, 700, 1000),
    time_o2_below_1pct = crossing_time(day, y, 1.0, "down", 0),
    time_o2_below_0p5pct = crossing_time(day, y, 0.5, "down", 0),
    time_o2_below_0p1pct = crossing_time(day, y, 0.1, "down", 0),
    time_o2_near_floor = crossing_time(day, y, floor_threshold, "down", 0),
    time_o2_near_floor_after_late_start = crossing_time(day, y, floor_threshold, "down", late_after_day),
    burden_day0 = interp_value(day, d$burden, 0),
    burden_day500 = interp_value(day, d$burden, 500),
    burden_day1000 = interp_value(day, d$burden, 1000),
    o2_source_file = d$source_file[[1]],
    stringsAsFactors = FALSE
  )
}

classify_seed <- function(ploidy_seed, stable_min_ploidy, stable_max_drop, collapse_terminal_max_ploidy,
                          collapse_min_late_drop) {
  p2 <- ploidy_seed[ploidy_seed$cohort == "2N", , drop = FALSE]
  p4 <- ploidy_seed[ploidy_seed$cohort == "4N", , drop = FALSE]
  if (nrow(p2) != 1L || nrow(p4) != 1L) {
    return(c(trajectory_regime = "ambiguous", mode_label = "ambiguous", mode_reason = "missing_2N_or_4N"))
  }
  terminal <- c(p2$terminal_median_ploidy, p4$terminal_median_ploidy)
  late_drop <- c(p2$late_drop_amplitude, p4$late_drop_amplitude)
  stable <- all(is.finite(terminal)) && all(is.finite(late_drop)) &&
    all(terminal >= stable_min_ploidy) && all(late_drop <= stable_max_drop)
  collapse <- all(is.finite(terminal)) && all(is.finite(late_drop)) &&
    all(terminal <= collapse_terminal_max_ploidy) && all(late_drop >= collapse_min_late_drop)
  if (stable) {
    c(trajectory_regime = "mode1_ploidy_stable", mode_label = "mode1", mode_reason = "both_cohorts_terminal_high_and_no_large_late_drop")
  } else if (collapse) {
    c(trajectory_regime = "mode2_second_ploidy_collapse", mode_label = "mode2", mode_reason = "both_cohorts_terminal_low_with_large_late_drop")
  } else {
    c(trajectory_regime = "ambiguous", mode_label = "ambiguous", mode_reason = "thresholds_not_jointly_met")
  }
}

wide_metric <- function(tab, seed_ids, cohort, metric) {
  sub <- tab[tab$cohort == cohort, , drop = FALSE]
  sub[[metric]][match(seed_ids, sub$seed_id)]
}

wilcox_row <- function(tab, value_col, group_col = "trajectory_regime",
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

plot_outputs <- function(curves, o2_curves, seed_summary, cohort_events, out_dir) {
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(curves)) {
    d <- merge(curves, seed_summary[, c("seed_id", "trajectory_regime")], by = "seed_id", all.x = TRUE)
    cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02", ambiguous = "grey60")
    grDevices::pdf(file.path(fig_dir, "ploidy_trajectories_by_regime.pdf"), width = 9, height = 6)
    plot(NA, xlim = range(d$day, na.rm = TRUE), ylim = range(d$ploidy, na.rm = TRUE), xlab = "Day", ylab = "Ploidy", main = "Ploidy trajectories by inferred regime")
    for (key in split(seq_len(nrow(d)), paste(d$seed_id, d$cohort, sep = "||"))) {
      x <- d[key, , drop = FALSE]
      col <- grDevices::adjustcolor(cols[x$trajectory_regime[[1]]] %||% "grey60", alpha.f = 0.25)
      lines(x$day, x$ploidy, col = col, lwd = 1)
    }
    legend("topright", legend = names(cols), col = cols, lwd = 2, bty = "n")
    grDevices::dev.off()
  }

  if (nrow(o2_curves)) {
    d <- merge(o2_curves, seed_summary[, c("seed_id", "trajectory_regime")], by = "seed_id", all.x = TRUE)
    cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02", ambiguous = "grey60")
    grDevices::pdf(file.path(fig_dir, "o2_timecourses_by_regime.pdf"), width = 9, height = 6)
    plot(NA, xlim = range(d$day, na.rm = TRUE), ylim = range(d$o2_pct, na.rm = TRUE), xlab = "Day", ylab = "O2 pct", main = "O2 timecourses by inferred regime")
    for (key in split(seq_len(nrow(d)), paste(d$seed_id, d$cohort, sep = "||"))) {
      x <- d[key, , drop = FALSE]
      col <- grDevices::adjustcolor(cols[x$trajectory_regime[[1]]] %||% "grey60", alpha.f = 0.25)
      lines(x$day, x$o2_pct, col = col, lwd = 1)
    }
    legend("topright", legend = names(cols), col = cols, lwd = 2, bty = "n")
    grDevices::dev.off()
  }

  if (nrow(cohort_events)) {
    grDevices::pdf(file.path(fig_dir, "late_ploidy_drop_vs_o2_floor_time.pdf"), width = 8, height = 6)
    cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02", ambiguous = "grey60")
    plot(cohort_events$time_o2_near_floor, cohort_events$time_crossing_ploidy_1p5_down_after_late_start,
      xlab = "Time O2 reaches near-floor", ylab = "Time ploidy crosses 1.5 after late start",
      main = "Late ploidy collapse timing vs O2 floor", pch = ifelse(cohort_events$cohort == "2N", 16, 17),
      col = cols[cohort_events$trajectory_regime] %||% "grey60"
    )
    abline(0, 1, lty = 2, col = "grey40")
    legend("topleft", legend = c("2N", "4N"), pch = c(16, 17), bty = "n")
    grDevices::dev.off()
  }
}

write_report <- function(out_dir, run_dir, seed_summary, coupling, objective_tests, parameter_tests, o2_tests) {
  mode_counts <- as.data.frame(table(seed_summary$trajectory_regime), stringsAsFactors = FALSE)
  names(mode_counts) <- c("trajectory_regime", "n")
  mode2 <- seed_summary[seed_summary$trajectory_regime == "mode2_second_ploidy_collapse", , drop = FALSE]
  mode1 <- seed_summary[seed_summary$trajectory_regime == "mode1_ploidy_stable", , drop = FALSE]
  top_params <- parameter_tests[is.finite(parameter_tests$BH_adjusted_p_value), , drop = FALSE]
  top_params <- top_params[order(top_params$BH_adjusted_p_value, top_params$wilcox_p_value), , drop = FALSE]
  top_o2 <- o2_tests[is.finite(o2_tests$BH_adjusted_p_value), , drop = FALSE]
  top_o2 <- top_o2[order(top_o2$BH_adjusted_p_value, top_o2$wilcox_p_value), , drop = FALSE]
  lines <- c(
    "# In vivo O2-ploidy event coupling analysis",
    "",
    paste0("- run_dir: ", run_dir),
    paste0("- analyzed seeds: ", nrow(seed_summary)),
    paste0("- mode1_ploidy_stable: ", nrow(mode1)),
    paste0("- mode2_second_ploidy_collapse: ", nrow(mode2)),
    paste0("- ambiguous: ", sum(seed_summary$trajectory_regime == "ambiguous", na.rm = TRUE)),
    "",
    "## Interpretation boundary",
    "",
    "This analysis uses existing fitted trajectories and predicted O2 timecourses. It can support or weaken mechanistic consistency claims, but it does not prove biological causality without counterfactual simulation or independent validation.",
    "",
    "## Event-coupling summary",
    "",
    paste(utils::capture.output(print(coupling, row.names = FALSE)), collapse = "\n"),
    "",
    "## Objective comparison",
    "",
    paste(utils::capture.output(print(objective_tests, row.names = FALSE)), collapse = "\n"),
    "",
    "## Top parameter associations",
    "",
    paste(utils::capture.output(print(utils::head(top_params, 15), row.names = FALSE)), collapse = "\n"),
    "",
    "## Top O2 feature associations",
    "",
    paste(utils::capture.output(print(utils::head(top_o2, 15), row.names = FALSE)), collapse = "\n")
  )
  dir.create(file.path(out_dir, "report"), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
  write_tsv(mode_counts, file.path(out_dir, "tables", "regime_counts.tsv"))
}

main <- function(argv = parse_args()) {
  run_dir <- normalizePath(as_chr(argv$run_dir), mustWork = FALSE)
  out_dir <- normalizePath(as_chr(argv$out_dir, file.path(run_dir, "..", "analysis", "invivo_o2_ploidy_event_coupling_500seed")), mustWork = FALSE)
  if (!dir.exists(run_dir)) stop("run_dir does not exist: ", run_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(vapply(file.path(out_dir, c("tables", "figures", "logs", "report")), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
  log_file <- file.path(out_dir, "logs", "analyze_invivo_o2_ploidy_event_coupling.log")
  log_con <- file(log_file, open = "wt")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  n_unit <- as_num(argv$n_unit, 22)
  max_seeds <- as_int(argv$max_seeds, NA_integer_)
  early_window <- as_num_pair(argv$early_window, c(0, 200))
  plateau_window <- as_num_pair(argv$plateau_window, c(200, 500))
  terminal_window <- as_num_pair(argv$terminal_window, c(900, 1000))
  late_after_day <- as_num(argv$late_after_day, 200)
  slope_lag_days <- as_int(argv$slope_lag_days, 20L)
  early_slope_threshold <- as_num(argv$early_slope_threshold, -0.005)
  late_slope_threshold <- as_num(argv$late_slope_threshold, -0.0025)
  near_floor_abs <- as_num(argv$o2_near_floor_abs, 0.05)
  stable_min_ploidy <- as_num(argv$stable_min_ploidy, 2.0)
  stable_max_drop <- as_num(argv$stable_max_drop, 0.4)
  collapse_terminal_max_ploidy <- as_num(argv$collapse_terminal_max_ploidy, 1.4)
  collapse_min_late_drop <- as_num(argv$collapse_min_late_drop, 0.35)
  generate_figures <- as_bool(argv$generate_figures, TRUE)

  run_args <- data.frame(
    argument = c("run_dir", "out_dir", "n_unit", "max_seeds", "early_window", "plateau_window", "terminal_window", "late_after_day", "slope_lag_days", "early_slope_threshold", "late_slope_threshold", "o2_near_floor_abs", "stable_min_ploidy", "stable_max_drop", "collapse_terminal_max_ploidy", "collapse_min_late_drop", "generate_figures"),
    value = c(run_dir, out_dir, n_unit, max_seeds, paste(early_window, collapse = ","), paste(plateau_window, collapse = ","), paste(terminal_window, collapse = ","), late_after_day, slope_lag_days, early_slope_threshold, late_slope_threshold, near_floor_abs, stable_min_ploidy, stable_max_drop, collapse_terminal_max_ploidy, collapse_min_late_drop, generate_figures),
    stringsAsFactors = FALSE
  )
  write_tsv(run_args, file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  message("run_dir: ", run_dir)
  message("out_dir: ", out_dir)

  seed_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
  seed_ids <- basename(seed_dirs)
  ord <- order(seed_number(seed_ids))
  seed_dirs <- seed_dirs[ord]
  seed_ids <- seed_ids[ord]
  if (!is.na(max_seeds) && max_seeds > 0L) {
    keep <- seq_len(min(max_seeds, length(seed_ids)))
    seed_dirs <- seed_dirs[keep]
    seed_ids <- seed_ids[keep]
  }
  message("Discovered seeds: ", length(seed_ids))

  manifest_rows <- list()
  param_rows <- list()
  ploidy_curves <- list()
  o2_curves <- list()
  ploidy_features <- list()
  o2_features <- list()

  for (i in seq_along(seed_ids)) {
    seed_id <- seed_ids[[i]]
    seed_dir <- seed_dirs[[i]]
    if (i %% 25L == 0L || i == 1L) message("Processing ", seed_id, " (", i, "/", length(seed_ids), ")")
    manifest_rows[[i]] <- read_seed_manifest_row(seed_id, seed_dir)
    params <- read_params(file.path(seed_dir, "best_params.tsv"))
    if (nrow(params)) {
      params$seed_id <- seed_id
      param_rows[[i]] <- params[, c("seed_id", "parameter", "value")]
    }
    o2_min <- param_value(params, "o2_min", 0)

    pc <- read_ploidy_curve(seed_id, seed_dir, n_unit)
    if (nrow(pc)) {
      names(pc)[names(pc) == "x"] <- "ploidy"
      ploidy_curves[[seed_id]] <- pc
      for (co in c("2N", "4N")) {
        d <- pc[pc$cohort == co, , drop = FALSE]
        if (nrow(d)) {
          ploidy_features[[paste(seed_id, co, sep = "||")]] <- summarize_ploidy_one(
            d, early_window, plateau_window, terminal_window, late_after_day,
            slope_lag_days, early_slope_threshold, late_slope_threshold
          )
        }
      }
    }

    oc <- read_o2_curve(seed_id, seed_dir)
    if (nrow(oc)) {
      o2_curves[[seed_id]] <- oc
      for (co in c("2N", "4N")) {
        d <- oc[oc$cohort == co, , drop = FALSE]
        if (nrow(d)) {
          o2_features[[paste(seed_id, co, sep = "||")]] <- summarize_o2_one(d, o2_min, near_floor_abs, late_after_day)
        }
      }
    }
  }

  manifest <- do.call(rbind, manifest_rows)
  params_long <- if (length(param_rows)) do.call(rbind, param_rows) else data.frame()
  ploidy_curve_all <- if (length(ploidy_curves)) do.call(rbind, ploidy_curves) else data.frame()
  o2_curve_all <- if (length(o2_curves)) do.call(rbind, o2_curves) else data.frame()
  ploidy_feature_all <- if (length(ploidy_features)) do.call(rbind, ploidy_features) else data.frame()
  o2_feature_all <- if (length(o2_features)) do.call(rbind, o2_features) else data.frame()

  write_tsv(manifest, file.path(out_dir, "tables", "seed_manifest.tsv"))
  write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))
  write_tsv(ploidy_feature_all, file.path(out_dir, "tables", "ploidy_event_features_by_cohort.tsv"))
  write_tsv(o2_feature_all, file.path(out_dir, "tables", "o2_event_features_by_cohort.tsv"))

  seed_summary <- manifest
  labels <- t(vapply(seed_ids, function(seed) classify_seed(
    ploidy_feature_all[ploidy_feature_all$seed_id == seed, , drop = FALSE],
    stable_min_ploidy, stable_max_drop, collapse_terminal_max_ploidy, collapse_min_late_drop
  ), character(3)))
  labels <- data.frame(seed_id = rownames(labels), labels, row.names = NULL, stringsAsFactors = FALSE)
  seed_summary <- merge(seed_summary, labels, by = "seed_id", all.x = TRUE, sort = FALSE)
  seed_summary$delta_objective <- seed_summary$objective - min(seed_summary$objective, na.rm = TRUE)

  for (co in c("2N", "4N")) {
    prefix <- paste0(co, "__")
    for (m in setdiff(names(ploidy_feature_all), c("seed_id", "cohort"))) {
      seed_summary[[paste0(prefix, m)]] <- wide_metric(ploidy_feature_all, seed_summary$seed_id, co, m)
    }
    for (m in setdiff(names(o2_feature_all), c("seed_id", "cohort", "o2_source_file"))) {
      seed_summary[[paste0(prefix, m)]] <- wide_metric(o2_feature_all, seed_summary$seed_id, co, m)
    }
  }

  seed_summary$mean_o2_day0_pct <- rowMeans(cbind(seed_summary[["2N__o2_day0_pct"]], seed_summary[["4N__o2_day0_pct"]]), na.rm = TRUE)
  seed_summary$mean_time_o2_near_floor <- rowMeans(cbind(seed_summary[["2N__time_o2_near_floor"]], seed_summary[["4N__time_o2_near_floor"]]), na.rm = TRUE)
  seed_summary$mean_terminal_ploidy <- rowMeans(cbind(seed_summary[["2N__terminal_median_ploidy"]], seed_summary[["4N__terminal_median_ploidy"]]), na.rm = TRUE)
  seed_summary$mean_late_drop_amplitude <- rowMeans(cbind(seed_summary[["2N__late_drop_amplitude"]], seed_summary[["4N__late_drop_amplitude"]]), na.rm = TRUE)
  seed_summary$late_2N_lag_p1p5_vs_o2_floor <- seed_summary[["2N__time_crossing_ploidy_1p5_down_after_late_start"]] - seed_summary[["2N__time_o2_near_floor"]]
  seed_summary$late_4N_lag_p1p5_vs_o2_floor <- seed_summary[["4N__time_crossing_ploidy_1p5_down_after_late_start"]] - seed_summary[["4N__time_o2_near_floor"]]
  seed_summary$early_4N_drop_o2_pct <- mapply(function(seed, tt) {
    d <- o2_curve_all[o2_curve_all$seed_id == seed & o2_curve_all$cohort == "4N", , drop = FALSE]
    if (!nrow(d)) return(NA_real_)
    interp_value(d$day, d$o2_pct, tt)
  }, seed_summary$seed_id, seed_summary[["4N__time_early_largest_negative_slope"]])
  seed_summary$early_4N_drop_time_minus_o2_floor <- seed_summary[["4N__time_early_largest_negative_slope"]] - seed_summary[["4N__time_o2_near_floor"]]

  write_tsv(seed_summary, file.path(out_dir, "tables", "seed_event_summary.tsv"))

  cohort_events <- merge(ploidy_feature_all, o2_feature_all, by = c("seed_id", "cohort"), all = TRUE, sort = FALSE)
  cohort_events <- merge(cohort_events, seed_summary[, c("seed_id", "trajectory_regime", "mode_label", "objective", "delta_objective")], by = "seed_id", all.x = TRUE, sort = FALSE)
  cohort_events$late_lag_p1p5_vs_o2_floor <- cohort_events$time_crossing_ploidy_1p5_down_after_late_start - cohort_events$time_o2_near_floor
  cohort_events$late_lag_largest_drop_vs_o2_floor <- cohort_events$time_late_largest_negative_slope - cohort_events$time_o2_near_floor
  cohort_events$early_lag_largest_drop_vs_o2_floor <- cohort_events$time_early_largest_negative_slope - cohort_events$time_o2_near_floor
  cohort_events$o2_at_early_largest_negative_slope <- mapply(function(seed, co, tt) {
    d <- o2_curve_all[o2_curve_all$seed_id == seed & o2_curve_all$cohort == co, , drop = FALSE]
    if (!nrow(d)) return(NA_real_)
    interp_value(d$day, d$o2_pct, tt)
  }, cohort_events$seed_id, cohort_events$cohort, cohort_events$time_early_largest_negative_slope)
  cohort_events$o2_at_late_largest_negative_slope <- mapply(function(seed, co, tt) {
    d <- o2_curve_all[o2_curve_all$seed_id == seed & o2_curve_all$cohort == co, , drop = FALSE]
    if (!nrow(d)) return(NA_real_)
    interp_value(d$day, d$o2_pct, tt)
  }, cohort_events$seed_id, cohort_events$cohort, cohort_events$time_late_largest_negative_slope)
  write_tsv(cohort_events, file.path(out_dir, "tables", "cohort_event_coupling.tsv"))

  objective_tests <- do.call(rbind, lapply(c("objective", "delta_objective", "optimizer_iter_completed", "mean_o2_day0_pct", "mean_time_o2_near_floor", "mean_terminal_ploidy", "mean_late_drop_amplitude"), function(v) {
    wilcox_row(seed_summary, v, feature_name = v)
  }))
  objective_tests$BH_adjusted_p_value <- stats::p.adjust(objective_tests$wilcox_p_value, method = "BH")
  write_tsv(objective_tests, file.path(out_dir, "tables", "objective_and_seed_level_regime_tests.tsv"))

  param_tests <- data.frame()
  if (nrow(params_long)) {
    params_wide <- reshape(params_long, idvar = "seed_id", timevar = "parameter", direction = "wide")
    names(params_wide) <- sub("^value\\.", "", names(params_wide))
    pd <- merge(seed_summary[, c("seed_id", "trajectory_regime")], params_wide, by = "seed_id", all.x = TRUE)
    pcols <- setdiff(names(pd), c("seed_id", "trajectory_regime"))
    param_tests <- do.call(rbind, lapply(pcols, function(v) {
      vals <- pd[[v]]
      if (all(vals > 0, na.rm = TRUE)) pd[[paste0(v, "__log10")]] <- log10(vals)
      use_col <- if (all(vals > 0, na.rm = TRUE)) paste0(v, "__log10") else v
      wilcox_row(pd, use_col, feature_name = v)
    }))
    param_tests$BH_adjusted_p_value <- stats::p.adjust(param_tests$wilcox_p_value, method = "BH")
  }
  write_tsv(param_tests, file.path(out_dir, "tables", "parameter_regime_association.tsv"))

  o2_cols <- grep("__(o2_|time_o2_|burden_)", names(seed_summary), value = TRUE)
  o2_tests <- if (length(o2_cols)) do.call(rbind, lapply(o2_cols, function(v) wilcox_row(seed_summary, v, feature_name = v))) else data.frame()
  if (nrow(o2_tests)) o2_tests$BH_adjusted_p_value <- stats::p.adjust(o2_tests$wilcox_p_value, method = "BH")
  write_tsv(o2_tests, file.path(out_dir, "tables", "o2_feature_regime_association.tsv"))

  mode2_cohort <- cohort_events[cohort_events$trajectory_regime == "mode2_second_ploidy_collapse", , drop = FALSE]
  coupling_rows <- list()
  coupling_row <- function(d, diagnostic, event_col, reference_col, event_o2_col = NA_character_) {
    event_time <- d[[event_col]]
    reference_time <- d[[reference_col]]
    lag <- event_time - reference_time
    ok <- is.finite(lag)
    event_o2 <- if (!is.na(event_o2_col) && event_o2_col %in% names(d)) d[[event_o2_col]] else rep(NA_real_, nrow(d))
    data.frame(
      diagnostic = diagnostic,
      event_time = event_col,
      o2_reference_time = reference_col,
      n = sum(ok),
      median_lag_days = stats::median(lag[ok], na.rm = TRUE),
      fraction_abs_lag_le_50d = mean(abs(lag[ok]) <= 50, na.rm = TRUE),
      fraction_ploidy_drop_after_or_near_o2_reference = mean(lag[ok] >= -50, na.rm = TRUE),
      spearman_time_correlation = if (sum(ok) >= 3L) suppressWarnings(stats::cor(event_time[ok], reference_time[ok], method = "spearman")) else NA_real_,
      median_event_o2_pct = stats::median(event_o2[is.finite(event_o2)], na.rm = TRUE),
      fraction_event_o2_gt_0p5pct = mean(event_o2[is.finite(event_o2)] > 0.5, na.rm = TRUE),
      fraction_event_o2_gt_1pct = mean(event_o2[is.finite(event_o2)] > 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  late_events <- c(
    late_largest_negative_slope = "time_late_largest_negative_slope",
    late_crossing_ploidy_1p5 = "time_crossing_ploidy_1p5_down_after_late_start"
  )
  o2_refs <- c(
    o2_below_1pct = "time_o2_below_1pct",
    o2_below_0p5pct = "time_o2_below_0p5pct",
    o2_below_0p1pct = "time_o2_below_0p1pct",
    o2_near_floor = "time_o2_near_floor"
  )
  for (co in c("2N", "4N")) {
    d <- mode2_cohort[mode2_cohort$cohort == co, , drop = FALSE]
    for (event_name in names(late_events)) {
      for (ref_name in names(o2_refs)) {
        coupling_rows[[paste("mode2", "late", co, event_name, ref_name, sep = "_")]] <- coupling_row(
          d,
          diagnostic = paste("mode2", "late", co, event_name, "vs", ref_name, sep = "_"),
          event_col = late_events[[event_name]],
          reference_col = o2_refs[[ref_name]],
          event_o2_col = if (identical(event_name, "late_largest_negative_slope")) "o2_at_late_largest_negative_slope" else NA_character_
        )
      }
    }
  }
  d4 <- mode2_cohort[mode2_cohort$cohort == "4N", , drop = FALSE]
  for (ref_name in names(o2_refs)) {
    coupling_rows[[paste("mode2", "early", "4N", "largest_negative_slope", ref_name, sep = "_")]] <- coupling_row(
      d4,
      diagnostic = paste("mode2", "early", "4N", "largest_negative_slope", "vs", ref_name, sep = "_"),
      event_col = "time_early_largest_negative_slope",
      reference_col = o2_refs[[ref_name]],
      event_o2_col = "o2_at_early_largest_negative_slope"
    )
  }
  coupling <- do.call(rbind, coupling_rows)
  write_tsv(coupling, file.path(out_dir, "tables", "event_coupling_diagnostics.tsv"))

  if (generate_figures) plot_outputs(ploidy_curve_all, o2_curve_all, seed_summary, cohort_events, out_dir)
  write_report(out_dir, run_dir, seed_summary, coupling, objective_tests, param_tests, o2_tests)

  message("Completed event coupling analysis: ", out_dir)
}

if (sys.nframe() == 0L) {
  main()
}
