#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  own_arg <- file_arg[
    basename(sub("^--file=", "", file_arg)) ==
      "analyze_invivo_o2_ploidy_event_coupling_tables.R"
  ]
  if (length(own_arg)) {
    dirname(normalizePath(sub("^--file=", "", own_arg[[1L]]), mustWork = FALSE))
  } else {
    frame_files <- Filter(
      nzchar,
      vapply(sys.frames(), function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      }, character(1))
    )
    own <- frame_files[
      basename(frame_files) ==
        "analyze_invivo_o2_ploidy_event_coupling_tables.R"
    ]
    if (length(own)) dirname(own[[length(own)]]) else getwd()
  }
})

source(file.path(SCRIPT_DIR, "process_fingerprint_utils.R"), local = TRUE)
source(file.path(SCRIPT_DIR, "ploidy_regime_utils.R"), local = TRUE)

`%||%` <- o2ipa_null_coalesce
parse_args <- o2ipa_parse_args
as_chr <- o2ipa_as_chr
as_num <- o2ipa_as_num
as_int <- o2ipa_as_int
as_bool <- o2ipa_as_bool

as_num_pair <- function(x, default) {
  txt <- as_chr(x, paste(default, collapse = ","))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2L) default else vals[seq_len(2L)]
}

read_tsv <- o2ipa_read_tsv
write_tsv <- o2ipa_write_tsv

safe_num <- o2ipa_numeric

seed_number <- o2ipa_seed_number

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

crossing_time <- o2pr_crossing_time

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

event_safe_read <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 1L) return(data.frame())
  tryCatch(read_tsv(path), error = function(e) data.frame())
}

event_require_simulation <- function(simulation_dir, files) {
  manifest <- file.path(simulation_dir, "simulation_manifest.tsv")
  if (!file.exists(manifest)) stop("Missing O2/ploidy event simulation manifest: ", manifest)
  missing <- files[!file.exists(file.path(simulation_dir, "tables", files))]
  if (length(missing)) stop("Missing O2/ploidy event simulation tables: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

event_write_analysis_manifest <- function(out_dir) {
  files <- list.files(file.path(out_dir, "tables"), pattern = "[.]tsv$", full.names = FALSE)
  rows <- lapply(files, function(file) {
    tab <- event_safe_read(file.path(out_dir, "tables", file))
    data.frame(
      artifact = tools::file_path_sans_ext(file), relative_path = file.path("tables", file),
      role = "analysis_table", rows = nrow(tab), columns = ncol(tab),
      path = normalizePath(file.path(out_dir, "tables", file), mustWork = TRUE),
      exists = TRUE, stringsAsFactors = FALSE
    )
  })
  write_tsv(do.call(rbind, rows), file.path(out_dir, "analysis_manifest.tsv"))
}

run_o2_ploidy_event_coupling_analysis <- function(argv = parse_args()) {
  simulation_dir <- as_chr(argv$simulation_dir)
  if (!nzchar(simulation_dir) || !dir.exists(simulation_dir)) stop("Missing or invalid --simulation_dir.")
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  out_dir <- as_chr(argv$analysis_dir %||% argv$out_dir)
  if (!nzchar(out_dir)) stop("Missing --analysis_dir (or --out_dir).")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  required <- c(
    "event_seed_manifest.tsv", "event_parameter_values_long.tsv",
    "event_ploidy_timecourses.tsv", "event_o2_timecourses.tsv"
  )
  event_require_simulation(simulation_dir, required)
  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

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

  manifest <- event_safe_read(file.path(simulation_dir, "tables", "event_seed_manifest.tsv"))
  params_long <- event_safe_read(file.path(simulation_dir, "tables", "event_parameter_values_long.tsv"))
  ploidy_curve_all <- event_safe_read(file.path(simulation_dir, "tables", "event_ploidy_timecourses.tsv"))
  o2_curve_all <- event_safe_read(file.path(simulation_dir, "tables", "event_o2_timecourses.tsv"))
  if ("x" %in% names(ploidy_curve_all) && !"ploidy" %in% names(ploidy_curve_all)) names(ploidy_curve_all)[names(ploidy_curve_all) == "x"] <- "ploidy"
  seed_ids <- as.character(manifest$seed_id)

  ploidy_features <- list()
  o2_features <- list()
  for (seed in seed_ids) {
    params <- params_long[params_long$seed_id == seed, , drop = FALSE]
    o2_min <- params$value[match("o2_min", params$parameter)]
    if (!length(o2_min) || !is.finite(o2_min)) o2_min <- 0
    for (cohort in c("2N", "4N")) {
      pd <- ploidy_curve_all[ploidy_curve_all$seed_id == seed & ploidy_curve_all$cohort == cohort, , drop = FALSE]
      if (nrow(pd)) {
        ploidy_features[[paste(seed, cohort, sep = "||")]] <- summarize_ploidy_one(
          pd, early_window, plateau_window, terminal_window, late_after_day,
          slope_lag_days, early_slope_threshold, late_slope_threshold
        )
      }
      od <- o2_curve_all[o2_curve_all$seed_id == seed & o2_curve_all$cohort == cohort, , drop = FALSE]
      if (nrow(od)) o2_features[[paste(seed, cohort, sep = "||")]] <- summarize_o2_one(od, o2_min, near_floor_abs, late_after_day)
    }
  }
  ploidy_feature_all <- if (length(ploidy_features)) do.call(rbind, ploidy_features) else data.frame()
  o2_feature_all <- if (length(o2_features)) do.call(rbind, o2_features) else data.frame()

  labels <- t(vapply(seed_ids, function(seed) classify_seed(
    ploidy_feature_all[ploidy_feature_all$seed_id == seed, , drop = FALSE],
    stable_min_ploidy, stable_max_drop, collapse_terminal_max_ploidy, collapse_min_late_drop
  ), character(3)))
  labels <- data.frame(seed_id = rownames(labels), labels, row.names = NULL, stringsAsFactors = FALSE)
  seed_summary <- merge(manifest, labels, by = "seed_id", all.x = TRUE, sort = FALSE)
  if (!"objective" %in% names(seed_summary)) seed_summary$objective <- NA_real_
  seed_summary$delta_objective <- seed_summary$objective - min(seed_summary$objective, na.rm = TRUE)
  if (!any(is.finite(seed_summary$objective))) seed_summary$delta_objective <- NA_real_
  for (cohort in c("2N", "4N")) {
    prefix <- paste0(cohort, "__")
    for (metric in setdiff(names(ploidy_feature_all), c("seed_id", "cohort"))) {
      seed_summary[[paste0(prefix, metric)]] <- wide_metric(ploidy_feature_all, seed_summary$seed_id, cohort, metric)
    }
    for (metric in setdiff(names(o2_feature_all), c("seed_id", "cohort", "o2_source_file"))) {
      seed_summary[[paste0(prefix, metric)]] <- wide_metric(o2_feature_all, seed_summary$seed_id, cohort, metric)
    }
  }
  mean_pair <- function(a, b) rowMeans(cbind(seed_summary[[a]], seed_summary[[b]]), na.rm = TRUE)
  seed_summary$mean_o2_day0_pct <- mean_pair("2N__o2_day0_pct", "4N__o2_day0_pct")
  seed_summary$mean_time_o2_near_floor <- mean_pair("2N__time_o2_near_floor", "4N__time_o2_near_floor")
  seed_summary$mean_terminal_ploidy <- mean_pair("2N__terminal_median_ploidy", "4N__terminal_median_ploidy")
  seed_summary$mean_late_drop_amplitude <- mean_pair("2N__late_drop_amplitude", "4N__late_drop_amplitude")
  seed_summary$late_2N_lag_p1p5_vs_o2_floor <- seed_summary[["2N__time_crossing_ploidy_1p5_down_after_late_start"]] - seed_summary[["2N__time_o2_near_floor"]]
  seed_summary$late_4N_lag_p1p5_vs_o2_floor <- seed_summary[["4N__time_crossing_ploidy_1p5_down_after_late_start"]] - seed_summary[["4N__time_o2_near_floor"]]
  seed_summary$early_4N_drop_o2_pct <- mapply(function(seed, tt) {
    d <- o2_curve_all[o2_curve_all$seed_id == seed & o2_curve_all$cohort == "4N", , drop = FALSE]
    if (!nrow(d)) NA_real_ else interp_value(d$day, d$o2_pct, tt)
  }, seed_summary$seed_id, seed_summary[["4N__time_early_largest_negative_slope"]])
  seed_summary$early_4N_drop_time_minus_o2_floor <- seed_summary[["4N__time_early_largest_negative_slope"]] - seed_summary[["4N__time_o2_near_floor"]]

  cohort_events <- merge(ploidy_feature_all, o2_feature_all, by = c("seed_id", "cohort"), all = TRUE, sort = FALSE)
  cohort_events <- merge(cohort_events, seed_summary[, c("seed_id", "trajectory_regime", "mode_label", "objective", "delta_objective")], by = "seed_id", all.x = TRUE, sort = FALSE)
  cohort_events$late_lag_p1p5_vs_o2_floor <- cohort_events$time_crossing_ploidy_1p5_down_after_late_start - cohort_events$time_o2_near_floor
  cohort_events$late_lag_largest_drop_vs_o2_floor <- cohort_events$time_late_largest_negative_slope - cohort_events$time_o2_near_floor
  cohort_events$early_lag_largest_drop_vs_o2_floor <- cohort_events$time_early_largest_negative_slope - cohort_events$time_o2_near_floor
  cohort_events$o2_at_early_largest_negative_slope <- mapply(function(seed, cohort, tt) {
    d <- o2_curve_all[o2_curve_all$seed_id == seed & o2_curve_all$cohort == cohort, , drop = FALSE]
    if (!nrow(d)) NA_real_ else interp_value(d$day, d$o2_pct, tt)
  }, cohort_events$seed_id, cohort_events$cohort, cohort_events$time_early_largest_negative_slope)
  cohort_events$o2_at_late_largest_negative_slope <- mapply(function(seed, cohort, tt) {
    d <- o2_curve_all[o2_curve_all$seed_id == seed & o2_curve_all$cohort == cohort, , drop = FALSE]
    if (!nrow(d)) NA_real_ else interp_value(d$day, d$o2_pct, tt)
  }, cohort_events$seed_id, cohort_events$cohort, cohort_events$time_late_largest_negative_slope)

  seed_test_vars <- c("objective", "delta_objective", "optimizer_iter_completed", "mean_o2_day0_pct", "mean_time_o2_near_floor", "mean_terminal_ploidy", "mean_late_drop_amplitude")
  for (name in setdiff(seed_test_vars, names(seed_summary))) seed_summary[[name]] <- NA_real_
  objective_tests <- do.call(rbind, lapply(seed_test_vars, function(v) wilcox_row(seed_summary, v, feature_name = v)))
  objective_tests$BH_adjusted_p_value <- stats::p.adjust(objective_tests$wilcox_p_value, method = "BH")

  param_tests <- data.frame()
  if (nrow(params_long)) {
    params_wide <- reshape(params_long[, c("seed_id", "parameter", "value")], idvar = "seed_id", timevar = "parameter", direction = "wide")
    names(params_wide) <- sub("^value\\.", "", names(params_wide))
    pd <- merge(seed_summary[, c("seed_id", "trajectory_regime")], params_wide, by = "seed_id", all.x = TRUE)
    param_tests <- do.call(rbind, lapply(setdiff(names(pd), c("seed_id", "trajectory_regime")), function(v) {
      vals <- pd[[v]]
      use_col <- v
      if (all(vals > 0, na.rm = TRUE)) {
        use_col <- paste0(v, "__log10")
        pd[[use_col]] <- log10(vals)
      }
      wilcox_row(pd, use_col, feature_name = v)
    }))
    param_tests$BH_adjusted_p_value <- stats::p.adjust(param_tests$wilcox_p_value, method = "BH")
  }
  o2_cols <- grep("__(o2_|time_o2_|burden_)", names(seed_summary), value = TRUE)
  o2_tests <- if (length(o2_cols)) do.call(rbind, lapply(o2_cols, function(v) wilcox_row(seed_summary, v, feature_name = v))) else data.frame()
  if (nrow(o2_tests)) o2_tests$BH_adjusted_p_value <- stats::p.adjust(o2_tests$wilcox_p_value, method = "BH")

  mode2 <- cohort_events[cohort_events$trajectory_regime == "mode2_second_ploidy_collapse", , drop = FALSE]
  coupling_row <- function(d, diagnostic, event_col, reference_col, event_o2_col = NA_character_) {
    event_time <- d[[event_col]]
    reference_time <- d[[reference_col]]
    lag <- event_time - reference_time
    ok <- is.finite(lag)
    event_o2 <- if (!is.na(event_o2_col) && event_o2_col %in% names(d)) d[[event_o2_col]] else rep(NA_real_, nrow(d))
    data.frame(
      diagnostic = diagnostic, event_time = event_col, o2_reference_time = reference_col,
      n = sum(ok), median_lag_days = if (any(ok)) stats::median(lag[ok]) else NA_real_,
      fraction_abs_lag_le_50d = if (any(ok)) mean(abs(lag[ok]) <= 50) else NA_real_,
      fraction_ploidy_drop_after_or_near_o2_reference = if (any(ok)) mean(lag[ok] >= -50) else NA_real_,
      spearman_time_correlation = if (sum(ok) >= 3L) suppressWarnings(stats::cor(event_time[ok], reference_time[ok], method = "spearman")) else NA_real_,
      median_event_o2_pct = if (any(is.finite(event_o2))) stats::median(event_o2[is.finite(event_o2)]) else NA_real_,
      fraction_event_o2_gt_0p5pct = if (any(is.finite(event_o2))) mean(event_o2[is.finite(event_o2)] > 0.5) else NA_real_,
      fraction_event_o2_gt_1pct = if (any(is.finite(event_o2))) mean(event_o2[is.finite(event_o2)] > 1) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  coupling_rows <- list()
  late_events <- c(late_largest_negative_slope = "time_late_largest_negative_slope", late_crossing_ploidy_1p5 = "time_crossing_ploidy_1p5_down_after_late_start")
  o2_refs <- c(o2_below_1pct = "time_o2_below_1pct", o2_below_0p5pct = "time_o2_below_0p5pct", o2_below_0p1pct = "time_o2_below_0p1pct", o2_near_floor = "time_o2_near_floor")
  for (cohort in c("2N", "4N")) for (event_name in names(late_events)) for (ref_name in names(o2_refs)) {
    d <- mode2[mode2$cohort == cohort, , drop = FALSE]
    key <- paste("mode2", "late", cohort, event_name, ref_name, sep = "_")
    coupling_rows[[key]] <- coupling_row(d, paste("mode2", "late", cohort, event_name, "vs", ref_name, sep = "_"), late_events[[event_name]], o2_refs[[ref_name]], if (event_name == "late_largest_negative_slope") "o2_at_late_largest_negative_slope" else NA_character_)
  }
  d4 <- mode2[mode2$cohort == "4N", , drop = FALSE]
  for (ref_name in names(o2_refs)) {
    key <- paste("mode2", "early", "4N", "largest_negative_slope", ref_name, sep = "_")
    coupling_rows[[key]] <- coupling_row(d4, paste("mode2", "early", "4N", "largest_negative_slope", "vs", ref_name, sep = "_"), "time_early_largest_negative_slope", o2_refs[[ref_name]], "o2_at_early_largest_negative_slope")
  }
  coupling <- do.call(rbind, coupling_rows)

  write_tsv(data.frame(
    argument = c("simulation_dir", "early_window", "plateau_window", "terminal_window", "late_after_day", "slope_lag_days", "early_slope_threshold", "late_slope_threshold", "o2_near_floor_abs", "stable_min_ploidy", "stable_max_drop", "collapse_terminal_max_ploidy", "collapse_min_late_drop"),
    value = c(simulation_dir, paste(early_window, collapse = ","), paste(plateau_window, collapse = ","), paste(terminal_window, collapse = ","), late_after_day, slope_lag_days, early_slope_threshold, late_slope_threshold, near_floor_abs, stable_min_ploidy, stable_max_drop, collapse_terminal_max_ploidy, collapse_min_late_drop),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "tables", "analysis_run_arguments.tsv"))
  write_tsv(manifest, file.path(out_dir, "tables", "seed_manifest.tsv"))
  write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))
  write_tsv(ploidy_feature_all, file.path(out_dir, "tables", "ploidy_event_features_by_cohort.tsv"))
  write_tsv(o2_feature_all, file.path(out_dir, "tables", "o2_event_features_by_cohort.tsv"))
  write_tsv(seed_summary, file.path(out_dir, "tables", "seed_event_summary.tsv"))
  write_tsv(cohort_events, file.path(out_dir, "tables", "cohort_event_coupling.tsv"))
  write_tsv(objective_tests, file.path(out_dir, "tables", "objective_and_seed_level_regime_tests.tsv"))
  write_tsv(param_tests, file.path(out_dir, "tables", "parameter_regime_association.tsv"))
  write_tsv(o2_tests, file.path(out_dir, "tables", "o2_feature_regime_association.tsv"))
  write_tsv(coupling, file.path(out_dir, "tables", "event_coupling_diagnostics.tsv"))
  write_tsv(as.data.frame(table(seed_summary$trajectory_regime), stringsAsFactors = FALSE), file.path(out_dir, "tables", "regime_counts.tsv"))
  event_write_analysis_manifest(out_dir)
  message("Completed O2-ploidy event coupling analysis: ", normalizePath(out_dir, mustWork = TRUE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_o2_ploidy_event_coupling_analysis()
