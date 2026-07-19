# Pure fixed-O2 analysis helpers.
# Functions in this module consume materialized tables and perform classification,
# statistical tests, correlations, reliability summaries, and replicate aggregation.


fixed_o2_mode_utils_loader_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  own <- frame_files[basename(frame_files) == "fixed_o2_analysis_utils.R"]
  if (length(own)) {
    dirname(own[[length(own)]])
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})
fixed_o2_mode_workflow_root <- normalizePath(
  file.path(fixed_o2_mode_utils_loader_dir, "..", ".."),
  mustWork = FALSE
)
source(
  file.path(
    fixed_o2_mode_workflow_root,
    "util",
    "o2_supply_demand_map_fixed_o2_mode_utils.R"
  ),
  local = TRUE
)
rm(
  fixed_o2_mode_utils_loader_dir,
  fixed_o2_mode_workflow_root
)


fixo2_mode_panel_names <- function() {
  c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1: dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2: dominant ploidy < 2"
  )
}


fixo2_mode_colors <- function(alpha = 1) {
  cols <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "#0072B2",
    mode2_attractor_dominant_ploidy_lt_2 = "#D55E00"
  )
  if (is.finite(alpha) && alpha < 1) cols <- grDevices::adjustcolor(cols, alpha.f = alpha)
  cols
}


fixo2_mode_labels_from_regime <- function(regime) {
  regimes <- fixo2_mode_regimes()
  out <- names(regimes)[match(as.character(regime), unname(regimes))]
  out[is.na(out)] <- NA_character_
  out
}



fixo2_mode_reference_o2 <- function(args = NULL) {
  val <- if (is.null(args)) NA_real_ else o2ipa_as_num(args$mode_reference_o2, NA_real_)
  if (!is.finite(val)) val <- 2
  val
}


fixo2_mode_fields <- function(dominant_ploidy) {
  dominant_ploidy <- suppressWarnings(as.numeric(dominant_ploidy))
  threshold <- fixo2_mode_threshold()
  regime <- ifelse(
    !is.finite(dominant_ploidy),
    NA_character_,
    ifelse(dominant_ploidy >= threshold, fixo2_mode_regimes()[["mode1"]], fixo2_mode_regimes()[["mode2"]])
  )
  data.frame(
    trajectory_regime = regime,
    mode_label = fixo2_mode_labels_from_regime(regime),
    mode_source = "fixed_o2_attractor_dominant_ploidy",
    mode_rule = "dominant_mean_ploidy >= 2 => mode1; dominant_mean_ploidy < 2 => mode2",
    mode_threshold_dominant_ploidy = threshold,
    stringsAsFactors = FALSE
  )
}



fixo2_apply_reference_modes <- function(tab, reference_modes) {
  if (!nrow(tab) || !is.data.frame(reference_modes) || !nrow(reference_modes) || !"seed_id" %in% names(tab)) return(tab)
  mode_cols <- intersect(c(
    "seed_id", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
    "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
    "mode_reference_dominant_mean_ploidy", "mode_reference_status",
    "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
  ), names(reference_modes))
  if (!all(c("seed_id", "trajectory_regime") %in% mode_cols)) return(tab)
  meta <- reference_modes[, mode_cols, drop = FALSE]
  meta <- meta[!duplicated(meta$seed_id), , drop = FALSE]
  replace_cols <- setdiff(mode_cols, "seed_id")
  tab[, intersect(replace_cols, names(tab))] <- NULL
  merge(tab, meta, by = "seed_id", all.x = TRUE, sort = FALSE)
}


fo2_load_optional_seed_metadata <- function(path, manifest) {
  if (nzchar(path) && file.exists(path)) {
    tab <- fo2_read_tsv(path)
    if (!"seed_id" %in% names(tab)) stop("label_file must contain seed_id: ", path)
    keep <- intersect(c(
      "seed_id", "trajectory_regime", "mode_label", "mode_reason", "mean_terminal_ploidy",
      "mean_late_drop_amplitude", "4N__terminal_median_ploidy",
      "4N__late_drop_amplitude", "4N__time_o2_near_floor"
    ), names(tab))
    out <- tab[, keep, drop = FALSE]
    names(out) <- sub("^trajectory_regime$", "source_trajectory_regime", names(out))
    names(out) <- sub("^mode_label$", "source_mode_label", names(out))
    names(out) <- sub("^mode_reason$", "source_mode_reason", names(out))
    return(out)
  }
  data.frame(seed_id = manifest$seed_id, stringsAsFactors = FALSE)
}


fo2_wilcox_row <- function(tab, value_col, group_col = "trajectory_regime",
                           a = "mode1_attractor_dominant_ploidy_ge_2", b = "mode2_attractor_dominant_ploidy_lt_2",
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

fo2_attractor_regime_summary <- function(attractors) {
  required <- c(
    "seed_id", "O2_pct", "trajectory_regime", "dominant_mean_ploidy",
    "dominant_fraction_N_le_25", "dominant_fraction_N_below_44",
    "selection_22_vs_44", "selection_44_vs_88", "spectral_gap",
    "eigenvector_nonnegative"
  )
  missing <- setdiff(required, names(attractors))
  if (length(missing)) {
    stop(
      "Attractor table is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  groups <- split(
    attractors,
    list(attractors$O2_pct, attractors$trajectory_regime),
    drop = TRUE
  )
  if (!length(groups)) return(data.frame())
  out <- do.call(rbind, lapply(groups, function(d) {
    data.frame(
      O2_pct = d$O2_pct[[1L]],
      trajectory_regime = d$trajectory_regime[[1L]],
      n_seed = length(unique(d$seed_id)),
      dominant_mean_ploidy_median = stats::median(
        d$dominant_mean_ploidy,
        na.rm = TRUE
      ),
      dominant_mean_ploidy_q25 = as.numeric(stats::quantile(
        d$dominant_mean_ploidy,
        0.25,
        na.rm = TRUE,
        names = FALSE
      )),
      dominant_mean_ploidy_q75 = as.numeric(stats::quantile(
        d$dominant_mean_ploidy,
        0.75,
        na.rm = TRUE,
        names = FALSE
      )),
      fraction_N_le_25_median = stats::median(
        d$dominant_fraction_N_le_25,
        na.rm = TRUE
      ),
      fraction_N_below_44_median = stats::median(
        d$dominant_fraction_N_below_44,
        na.rm = TRUE
      ),
      selection_22_vs_44_median = stats::median(
        d$selection_22_vs_44,
        na.rm = TRUE
      ),
      selection_44_vs_88_median = stats::median(
        d$selection_44_vs_88,
        na.rm = TRUE
      ),
      spectral_gap_median = stats::median(d$spectral_gap, na.rm = TRUE),
      eigenvector_nonnegative_fraction = mean(
        d$eigenvector_nonnegative %in% TRUE,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(out$O2_pct, out$trajectory_regime), , drop = FALSE]
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
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
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
    d$trajectory_regime == "mode1_attractor_dominant_ploidy_ge_2",
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
  regimes <- c(mode1 = "mode1_attractor_dominant_ploidy_ge_2", mode2 = "mode2_attractor_dominant_ploidy_lt_2")
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
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
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


fo2_ploidy_gap_reliability_violin_table <- function(gap_by_seed,
                                                    o2_values = c(0, 0.1, 0.2, 0.5, 0.75, 1, 2, 3, 4, 5)) {
  regimes <- c(mode1 = "mode1_attractor_dominant_ploidy_ge_2", mode2 = "mode2_attractor_dominant_ploidy_lt_2")
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


cf2_wilcox_row <- function(tab, value_col, group_col = "trajectory_regime",
                           a = "mode1_attractor_dominant_ploidy_ge_2", b = "mode2_attractor_dominant_ploidy_lt_2",
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
  regs <- unname(fixo2_mode_regimes())
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
    mode1_mode2_only = merged[merged$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
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


fixo2_select_best_objective_seed_by_mode <- function(metadata, mode_reference_o2) {
  required_modes <- c("mode1", "mode2")
  if (!nrow(metadata)) stop("Cannot select simulation representative seeds; no seed mode/objective metadata was available.")
  rows <- lapply(required_modes, function(mode) {
    d <- metadata[metadata$mode_label == mode & is.finite(metadata$objective), , drop = FALSE]
    if (!nrow(d)) {
      stop(
        "Cannot select a simulation representative seed for ",
        mode,
        "; no finite final objective was available among seeds assigned to this mode at O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE),
        "."
      )
    }
    d <- d[order(d$objective, o2ipa_seed_number(d$seed_id)), , drop = FALSE]
    d <- d[1L, , drop = FALSE]
    data.frame(
      selection_mode = "best_final_objective_by_reference_mode",
      selection_rule = paste0(
        "Within ",
        mode,
        ", choose the seed with the smallest final objective among seeds assigned by the fixed-O2 attractor at O2=",
        format(mode_reference_o2, scientific = FALSE, trim = TRUE)
      ),
      mode_label = mode,
      trajectory_regime = d$trajectory_regime[[1]],
      seed_id = d$seed_id[[1]],
      final_objective = d$objective[[1]],
      objective_source = d$objective_source[[1]],
      mode_reference_o2_pct = if ("mode_reference_o2_pct" %in% names(d)) d$mode_reference_o2_pct[[1]] else mode_reference_o2,
      mode_reference_dominant_mean_ploidy = if ("mode_reference_dominant_mean_ploidy" %in% names(d)) d$mode_reference_dominant_mean_ploidy[[1]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$selection_rank <- seq_len(nrow(out))
  out[, c("selection_rank", setdiff(names(out), "selection_rank")), drop = FALSE]
}


fixo2_manual_simulation_seed_table <- function(seed_ids, seed_modes, metadata, mode_reference_o2) {
  seed_ids <- o2ipa_norm_seed(seed_ids)
  rows <- lapply(seq_along(seed_ids), function(i) {
    seed <- seed_ids[[i]]
    meta <- metadata[metadata$seed_id == seed, , drop = FALSE]
    mode <- seed_modes[[i]]
    if (!nzchar(mode) && nrow(meta) && "mode_label" %in% names(meta)) mode <- as.character(meta$mode_label[[1]])
    data.frame(
      selection_rank = i,
      selection_mode = "manual_seed_ids",
      selection_rule = "Use the explicit seed_ids supplied to the simulation validation step; mode and objective are attached from the FixO2 reference mode/objective table when available.",
      mode_label = mode,
      trajectory_regime = if (nrow(meta) && "trajectory_regime" %in% names(meta)) meta$trajectory_regime[[1]] else NA_character_,
      seed_id = seed,
      final_objective = if (nrow(meta) && "objective" %in% names(meta)) meta$objective[[1]] else NA_real_,
      objective_source = if (nrow(meta) && "objective_source" %in% names(meta)) meta$objective_source[[1]] else NA_character_,
      mode_reference_o2_pct = if (nrow(meta) && "mode_reference_o2_pct" %in% names(meta)) meta$mode_reference_o2_pct[[1]] else mode_reference_o2,
      mode_reference_dominant_mean_ploidy = if (nrow(meta) && "mode_reference_dominant_mean_ploidy" %in% names(meta)) meta$mode_reference_dominant_mean_ploidy[[1]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


aggregate_replicates <- function(sim_metrics) {
  if (!nrow(sim_metrics)) return(sim_metrics)
  numeric_cols <- intersect(
    c("day", "O2_pct", "initial_ploidy", "simulation_id", "simulation_mean_N",
      "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells"),
    names(sim_metrics)
  )
  for (col in numeric_cols) sim_metrics[[col]] <- suppressWarnings(as.numeric(sim_metrics[[col]]))
  keys <- c("seed_id", "O2_pct", "O2_key", "initial_condition", "initial_ploidy", "day", "day_key")
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(sim_metrics)
    out <- dt[, .(
      simulation_n = data.table::uniqueN(simulation_id),
      simulation_mean_N = mean(simulation_mean_N, na.rm = TRUE),
      simulation_sd_replicate_mean_N = stats::sd(simulation_mean_N, na.rm = TRUE),
      simulation_mean_ploidy = mean(simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_replicate_mean_ploidy = stats::sd(simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(simulation_sd_ploidy, na.rm = TRUE),
      simulation_live_cells = mean(simulation_live_cells, na.rm = TRUE)
    ), by = keys]
    return(as.data.frame(out))
  }
  split_key <- interaction(sim_metrics[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(sim_metrics, split_key), function(x) {
    data.frame(
      x[1, keys, drop = FALSE],
      simulation_n = length(unique(x$simulation_id)),
      simulation_mean_N = mean(x$simulation_mean_N, na.rm = TRUE),
      simulation_sd_replicate_mean_N = stats::sd(x$simulation_mean_N, na.rm = TRUE),
      simulation_mean_ploidy = mean(x$simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_replicate_mean_ploidy = stats::sd(x$simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(x$simulation_sd_ploidy, na.rm = TRUE),
      simulation_live_cells = mean(x$simulation_live_cells, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


merge_scatter_data <- function(analytical, sim_summary, objectives) {
  by_cols <- c("seed_id", "O2_key", "initial_condition", "day_key")
  dat <- merge(
    analytical,
    sim_summary,
    by = by_cols,
    all = FALSE,
    suffixes = c("_analytical", "_simulation")
  )
  if ("O2_pct_analytical" %in% names(dat)) dat$O2_pct <- dat$O2_pct_analytical
  if ("day_analytical" %in% names(dat)) dat$day <- dat$day_analytical
  objective_by <- "seed_id"
  if ("O2_key" %in% names(dat) && "O2_key" %in% names(objectives)) objective_by <- c("seed_id", "O2_key")
  dat <- merge(dat, objectives, by = objective_by, all.x = TRUE, suffixes = c("", "_objective"))
  if ("mode_label_objective" %in% names(dat) && "mode_label" %in% names(dat)) {
    fill <- !nzchar(as.character(dat$mode_label)) | is.na(dat$mode_label)
    dat$mode_label[fill] <- dat$mode_label_objective[fill]
  }
  dat$objective <- as.numeric(dat$objective)
  dat[is.finite(dat$analytical_mean_ploidy) & is.finite(dat$simulation_mean_ploidy), , drop = FALSE]
}

augment_comparison_data <- function(dat) {
  required <- c("analytical_mean_ploidy", "simulation_mean_ploidy")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(
      "Agreement table is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  dat$residual_ploidy <- as.numeric(dat$simulation_mean_ploidy) -
    as.numeric(dat$analytical_mean_ploidy)
  dat$agreement_mean_ploidy <- (
    as.numeric(dat$simulation_mean_ploidy) +
      as.numeric(dat$analytical_mean_ploidy)
  ) / 2
  dat[
    is.finite(dat$residual_ploidy) &
      is.finite(dat$agreement_mean_ploidy),
    ,
    drop = FALSE
  ]
}

bland_altman_stats <- function(dat) {
  vals <- as.numeric(dat$residual_ploidy)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(c(bias = 0, lower = NA_real_, upper = NA_real_))
  }
  bias <- mean(vals)
  sd_val <- stats::sd(vals)
  if (!is.finite(sd_val)) sd_val <- 0
  c(
    bias = bias,
    lower = bias - 1.96 * sd_val,
    upper = bias + 1.96 * sd_val
  )
}

fixo2_add_agreement_limits <- function(dat) {
  dat <- augment_comparison_data(dat)
  groups <- if ("analytical_method" %in% names(dat)) {
    split(seq_len(nrow(dat)), as.character(dat$analytical_method), drop = TRUE)
  } else {
    list(analytical = seq_len(nrow(dat)))
  }
  dat$bland_altman_bias <- NA_real_
  dat$bland_altman_lower <- NA_real_
  dat$bland_altman_upper <- NA_real_
  for (idx in groups) {
    stats <- bland_altman_stats(dat[idx, , drop = FALSE])
    dat$bland_altman_bias[idx] <- unname(stats[["bias"]])
    dat$bland_altman_lower[idx] <- unname(stats[["lower"]])
    dat$bland_altman_upper[idx] <- unname(stats[["upper"]])
  }
  dat
}

fixo2_agreement_bland_altman_summary <- function(dat) {
  dat <- augment_comparison_data(dat)
  grouping <- intersect(
    c("analytical_method", "O2_pct", "day", "initial_condition"),
    names(dat)
  )
  if (!length(grouping)) {
    stats <- bland_altman_stats(dat)
    return(data.frame(
      n = nrow(dat),
      bias = unname(stats[["bias"]]),
      lower_loa = unname(stats[["lower"]]),
      upper_loa = unname(stats[["upper"]]),
      stringsAsFactors = FALSE
    ))
  }
  key <- interaction(dat[grouping], drop = TRUE, lex.order = TRUE)
  groups <- split(dat, key, drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(d) {
    stats <- bland_altman_stats(d)
    row <- d[1L, grouping, drop = FALSE]
    row$n <- nrow(d)
    row$bias <- unname(stats[["bias"]])
    row$lower_loa <- unname(stats[["lower"]])
    row$upper_loa <- unname(stats[["upper"]])
    row
  }))
  rownames(out) <- NULL
  out
}
