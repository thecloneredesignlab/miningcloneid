#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
ANALYSIS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
WORKFLOW_DIR <- normalizePath(file.path(ANALYSIS_DIR, ".."), mustWork = FALSE)
REPO_ROOT <- normalizePath(file.path(WORKFLOW_DIR, "..", "..", ".."), mustWork = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- val
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x)) return(default)
  val <- tolower(as.character(x[[1]]))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer(x[[1]] %||% default))
  if (length(val) && is.finite(val)) val else default
}

as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric(x[[1]] %||% default))
  if (length(val) && is.finite(val)) val else default
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(as.character(x[[1]]))) return(default)
  vals <- suppressWarnings(as.numeric(trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

o2_key <- function(x) {
  vapply(x, function(xx) format(signif(as.numeric(xx), 12), scientific = FALSE, trim = TRUE), character(1))
}

o2_slug <- function(x) {
  key <- o2_key(x)
  key <- gsub("-", "minus", key, fixed = TRUE)
  key <- gsub("[^0-9A-Za-z]+", "p", key)
  key <- gsub("^p+|p+$", "", key)
  ifelse(nzchar(key), key, "NA")
}

seed_number <- function(seed_id) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

default_run_dir <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "fit_invivo_O2_buffering_500seed")
}

default_out_dir <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "analysis", "dense-grid_monotonicity_classification")
}

fixo2_script_path <- function() {
  file.path(WORKFLOW_DIR, "simulation", "fix_o2_simulation.R")
}

load_fixo2_env <- function(script_path = fixo2_script_path()) {
  script_path <- normalizePath(path.expand(script_path), mustWork = TRUE)
  env <- new.env(parent = globalenv())
  env$commandArgs <- function(trailingOnly = FALSE) {
    if (isTRUE(trailingOnly)) character(0) else paste0("--file=", script_path)
  }
  old_error <- getOption("error")
  on.exit(options(error = old_error), add = TRUE)
  sys.source(script_path, envir = env, chdir = TRUE)
  options(error = old_error)
  required <- c(
    "generate_fixo2_attractor_mode_table",
    "fixo2_attractor_mode_table",
    "fixo2_attractor_mode_summary_by_seed",
    "fixo2_reference_mode_table"
  )
  missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = TRUE)]
  if (length(missing)) stop("Fixed-O2 helper environment is missing: ", paste(missing, collapse = ", "))
  env
}

source(file.path(SCRIPT_DIR, "curve_classification_utils.R"))

threshold_crossings <- function(o2, y, target) {
  o2 <- as.numeric(o2)
  y <- as.numeric(y)
  ok <- is.finite(o2) & is.finite(y)
  o2 <- o2[ok]
  y <- y[ok]
  if (length(o2) < 2L) return(list(n = 0L, first = NA_real_, all = NA_character_, status = "insufficient_grid"))
  ord <- order(o2)
  o2 <- o2[ord]
  y <- y[ord]
  d <- y - target
  hits <- numeric()
  exact <- which(abs(d) < 1e-12)
  if (length(exact)) hits <- c(hits, o2[exact])
  for (i in seq_len(length(d) - 1L)) {
    if (!is.finite(d[[i]]) || !is.finite(d[[i + 1L]])) next
    if (d[[i]] == 0 || d[[i + 1L]] == 0) next
    if (sign(d[[i]]) != sign(d[[i + 1L]])) {
      frac <- (target - y[[i]]) / (y[[i + 1L]] - y[[i]])
      hits <- c(hits, o2[[i]] + frac * (o2[[i + 1L]] - o2[[i]]))
    }
  }
  hits <- sort(unique(round(hits, 10)))
  status <- if (length(hits)) {
    "crossed"
  } else if (all(y < target, na.rm = TRUE)) {
    "always_below"
  } else if (all(y > target, na.rm = TRUE)) {
    "always_above"
  } else {
    "no_crossing_detected"
  }
  list(
    n = length(hits),
    first = if (length(hits)) hits[[1L]] else NA_real_,
    all = if (length(hits)) paste(format(hits, scientific = FALSE, trim = TRUE), collapse = ";") else NA_character_,
    status = status
  )
}

classify_one_seed <- function(curve,
                              reporting_o2,
                              gap_low,
                              gap_caution,
                              unreliable_fraction,
                              caution_fraction,
                              flat_range_threshold,
                              step_epsilon_abs,
                              step_epsilon_fraction,
                              reverse_fraction_tolerance,
                              plateau_min_points) {
  curve <- curve[order(curve$O2_pct), , drop = FALSE]
  y <- suppressWarnings(as.numeric(curve$dominant_mean_ploidy))
  gap <- suppressWarnings(as.numeric(curve$spectral_gap))
  class_info <- classify_o2_ploidy_curve(
    curve,
    value_col = "dominant_mean_ploidy",
    x_col = "O2_pct",
    id_col = "seed_id",
    flat_range_threshold = flat_range_threshold,
    step_epsilon_abs = step_epsilon_abs,
    step_epsilon_fraction = step_epsilon_fraction,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    plateau_min_points = plateau_min_points
  )
  class_summary <- class_info$summary
  curve_class <- class_summary$curve_class[[1L]]
  min_idx <- which.min(y)
  max_idx <- which.max(y)
  frac_low <- mean(gap < gap_low, na.rm = TRUE)
  frac_caution <- mean(gap < gap_caution, na.rm = TRUE)
  min_gap <- suppressWarnings(min(gap, na.rm = TRUE))
  if (!is.finite(min_gap)) min_gap <- NA_real_
  reliability <- if (is.finite(frac_low) && frac_low >= unreliable_fraction) {
    "unreliable_small_gap"
  } else if ((is.finite(frac_low) && frac_low > 0) || (is.finite(frac_caution) && frac_caution >= caution_fraction)) {
    "caution_small_gap"
  } else {
    "reliable"
  }
  out <- data.frame(
    seed_id = curve$seed_id[[1L]],
    seed_number = seed_number(curve$seed_id[[1L]]),
    n_o2 = nrow(curve),
    o2_min = min(curve$O2_pct, na.rm = TRUE),
    o2_max = max(curve$O2_pct, na.rm = TRUE),
    curve_class = curve_class,
    final_interpretation_class = if (identical(reliability, "unreliable_small_gap")) "unreliable_small_spectral_gap" else curve_class,
    sign_sequence = class_summary$sign_sequence[[1L]],
    n_sign_changes = class_summary$n_sign_changes[[1L]],
    step_epsilon = class_summary$step_epsilon[[1L]],
    slope_epsilon = class_summary$slope_epsilon[[1L]],
    flat_range_threshold = class_summary$flat_range_threshold[[1L]],
    reverse_fraction_tolerance = class_summary$reverse_fraction_tolerance[[1L]],
    ploidy_range = class_summary$ploidy_range[[1L]],
    net_ploidy_change = class_summary$net_ploidy_change[[1L]],
    max_positive_step = class_summary$max_positive_step[[1L]],
    max_negative_step = class_summary$max_negative_step[[1L]],
    positive_step_total = class_summary$positive_step_total[[1L]],
    negative_step_total = class_summary$negative_step_total[[1L]],
    fraction_positive_steps = class_summary$fraction_positive_steps[[1L]],
    fraction_negative_steps = class_summary$fraction_negative_steps[[1L]],
    fraction_zero_steps = class_summary$fraction_zero_steps[[1L]],
    low_amplitude_curve = class_summary$low_amplitude_curve[[1L]],
    terminal_plateau = class_summary$terminal_plateau[[1L]],
    terminal_plateau_for_class = class_summary$terminal_plateau_for_class[[1L]],
    classification_rule_version = class_summary$classification_rule_version[[1L]],
    min_ploidy = y[[min_idx]],
    o2_at_min_ploidy = curve$O2_pct[[min_idx]],
    max_ploidy = y[[max_idx]],
    o2_at_max_ploidy = curve$O2_pct[[max_idx]],
    dominant_growth_rate_min = suppressWarnings(min(curve$dominant_growth_rate, na.rm = TRUE)),
    dominant_growth_rate_max = suppressWarnings(max(curve$dominant_growth_rate, na.rm = TRUE)),
    min_spectral_gap = min_gap,
    median_spectral_gap = stats::median(gap, na.rm = TRUE),
    fraction_o2_gap_below_0p005 = mean(gap < 0.005, na.rm = TRUE),
    fraction_o2_gap_below_0p01 = mean(gap < 0.01, na.rm = TRUE),
    monotonicity_reliability = reliability,
    stringsAsFactors = FALSE
  )
  for (target in c(1.5, 2.0)) {
    cr <- threshold_crossings(curve$O2_pct, y, target)
    key <- o2_slug(target)
    out[[paste0("n_crossings_ploidy_", key)]] <- cr$n
    out[[paste0("first_o2_crossing_ploidy_", key)]] <- cr$first
    out[[paste0("all_o2_crossings_ploidy_", key)]] <- cr$all
    out[[paste0("crossing_status_ploidy_", key)]] <- cr$status
  }
  for (o2 in reporting_o2) {
    hit <- which(abs(curve$O2_pct - o2) < 1e-9)
    key <- paste0("o2_", o2_slug(o2))
    out[[paste0("mode_label_", key)]] <- if (length(hit)) curve$mode_label[[hit[[1L]]]] else NA_character_
    out[[paste0("dominant_mean_ploidy_", key)]] <- if (length(hit)) curve$dominant_mean_ploidy[[hit[[1L]]]] else NA_real_
    out[[paste0("spectral_gap_", key)]] <- if (length(hit)) curve$spectral_gap[[hit[[1L]]]] else NA_real_
  }
  out
}

add_curve_differences <- function(curves, by_seed, reporting_o2) {
  eps_col <- if ("step_epsilon" %in% names(by_seed)) "step_epsilon" else "slope_epsilon"
  eps <- by_seed[[eps_col]]
  names(eps) <- by_seed$seed_id
  drop_cols <- intersect(c("finite_difference_next", "local_slope_sign", "step_epsilon", "is_reporting_o2"), names(curves))
  if (length(drop_cols)) curves[, drop_cols] <- NULL
  rows <- lapply(split(curves, curves$seed_id), function(d) {
    finite_diff_curve(d, eps[[d$seed_id[[1L]]]])
  })
  diffs <- do.call(rbind, rows)
  out <- merge(curves, diffs, by = c("seed_id", "O2_pct"), all.x = TRUE, sort = FALSE)
  out$is_reporting_o2 <- vapply(out$O2_pct, function(x) any(abs(x - reporting_o2) < 1e-9), logical(1))
  out[order(seed_number(out$seed_id), out$O2_pct), , drop = FALSE]
}

curve_quantile_table <- function(curves, by_seed) {
  d <- merge(curves, by_seed[, c("seed_id", "curve_class", "final_interpretation_class", "monotonicity_reliability")],
             by = "seed_id", all.x = TRUE, sort = FALSE)
  rows <- lapply(split(d, list(d$curve_class, d$O2_pct), drop = TRUE), function(x) {
    vals <- x$dominant_mean_ploidy
    data.frame(
      curve_class = x$curve_class[[1L]],
      O2_pct = x$O2_pct[[1L]],
      n_seed = length(unique(x$seed_id)),
      median_dominant_mean_ploidy = stats::median(vals, na.rm = TRUE),
      q25_dominant_mean_ploidy = as.numeric(stats::quantile(vals, 0.25, na.rm = TRUE, names = FALSE)),
      q75_dominant_mean_ploidy = as.numeric(stats::quantile(vals, 0.75, na.rm = TRUE, names = FALSE)),
      median_spectral_gap = stats::median(x$spectral_gap, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$curve_class, out$O2_pct), , drop = FALSE]
}

representative_seed_table <- function(curves, by_seed) {
  reps <- lapply(split(by_seed, by_seed$curve_class), function(seed_tab) {
    cls <- seed_tab$curve_class[[1L]]
    class_n_seed <- nrow(seed_tab)
    low_amplitude_excluded <- FALSE
    if ("low_amplitude_curve" %in% names(seed_tab) && any(seed_tab$low_amplitude_curve %in% FALSE, na.rm = TRUE)) {
      candidate_seed_tab <- seed_tab[seed_tab$low_amplitude_curve %in% FALSE, , drop = FALSE]
      low_amplitude_excluded <- nrow(candidate_seed_tab) < class_n_seed
      seed_tab <- candidate_seed_tab
    }
    representative_candidate_n_seed <- nrow(seed_tab)
    class_curves <- curves[curves$seed_id %in% seed_tab$seed_id, , drop = FALSE]
    med <- aggregate(class_curves$dominant_mean_ploidy, by = list(O2_pct = class_curves$O2_pct),
                     FUN = stats::median, na.rm = TRUE)
    names(med)[2] <- "class_median"
    scores <- lapply(split(class_curves, class_curves$seed_id), function(curve) {
      d <- merge(curve[, c("seed_id", "O2_pct", "dominant_mean_ploidy")], med, by = "O2_pct", all.x = TRUE)
      data.frame(seed_id = curve$seed_id[[1L]], median_curve_rmse = sqrt(mean((d$dominant_mean_ploidy - d$class_median)^2, na.rm = TRUE)))
    })
    score <- do.call(rbind, scores)
    hit <- score[which.min(score$median_curve_rmse), , drop = FALSE]
    meta <- seed_tab[match(hit$seed_id, seed_tab$seed_id), , drop = FALSE]
    data.frame(
      curve_class = cls,
      representative_seed_id = hit$seed_id,
      representative_seed_number = seed_number(hit$seed_id),
      median_curve_rmse = hit$median_curve_rmse,
      class_n_seed = class_n_seed,
      representative_candidate_n_seed = representative_candidate_n_seed,
      representative_low_amplitude_excluded = low_amplitude_excluded,
      representative_reliability = meta$monotonicity_reliability,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, reps)
  out[order(out$curve_class), , drop = FALSE]
}

class_count_table <- function(by_seed) {
  tab <- as.data.frame(table(by_seed$curve_class), stringsAsFactors = FALSE)
  names(tab) <- c("curve_class", "n_seed")
  tab$fraction_seed <- tab$n_seed / sum(tab$n_seed)
  tab[order(-tab$n_seed, tab$curve_class), , drop = FALSE]
}

class_by_mode_table <- function(by_seed, reporting_o2) {
  rows <- list()
  for (o2 in reporting_o2) {
    col <- paste0("mode_label_o2_", o2_slug(o2))
    if (!col %in% names(by_seed)) next
    tab <- as.data.frame(table(by_seed$curve_class, by_seed[[col]]), stringsAsFactors = FALSE)
    names(tab) <- c("curve_class", "mode_label", "n_seed")
    tab <- tab[nzchar(tab$mode_label) & tab$n_seed > 0, , drop = FALSE]
    tab$reference_o2 <- o2
    denom <- aggregate(tab$n_seed, by = list(reference_o2 = tab$reference_o2, curve_class = tab$curve_class), FUN = sum)
    names(denom)[3] <- "class_n_at_o2"
    tab <- merge(tab, denom, by = c("reference_o2", "curve_class"), all.x = TRUE, sort = FALSE)
    tab$fraction_within_curve_class <- tab$n_seed / tab$class_n_at_o2
    rows[[as.character(o2)]] <- tab
  }
  out <- do.call(rbind, rows)
  out[order(out$reference_o2, out$curve_class, out$mode_label), , drop = FALSE]
}

objective_rank_class_table <- function(by_seed) {
  d <- by_seed
  if (!"objective" %in% names(d) || all(!is.finite(d$objective))) return(data.frame())
  d <- d[order(d$objective), , drop = FALSE]
  d$objective_rank <- seq_len(nrow(d))
  d$objective_rank_bin <- cut(
    d$objective_rank,
    breaks = unique(c(0, 25, 50, 100, 250, nrow(d))),
    include.lowest = TRUE,
    right = TRUE
  )
  tab <- as.data.frame(table(d$objective_rank_bin, d$curve_class), stringsAsFactors = FALSE)
  names(tab) <- c("objective_rank_bin", "curve_class", "n_seed")
  tab <- tab[tab$n_seed > 0, , drop = FALSE]
  denom <- aggregate(tab$n_seed, by = list(objective_rank_bin = tab$objective_rank_bin), FUN = sum)
  names(denom)[2] <- "bin_n_seed"
  tab <- merge(tab, denom, by = "objective_rank_bin", all.x = TRUE, sort = FALSE)
  tab$fraction_within_rank_bin <- tab$n_seed / tab$bin_n_seed
  tab
}

parameter_difference_table <- function(by_seed, params_long) {
  if (!nrow(params_long) || !"curve_class" %in% names(by_seed)) return(data.frame())
  d <- merge(params_long, by_seed[, c("seed_id", "curve_class")], by = "seed_id", all.x = FALSE, sort = FALSE)
  rows <- lapply(split(d, list(d$parameter, d$curve_class), drop = TRUE), function(x) {
    data.frame(
      parameter = x$parameter[[1L]],
      curve_class = x$curve_class[[1L]],
      n_seed = length(unique(x$seed_id)),
      median_value = stats::median(x$value, na.rm = TRUE),
      q25_value = as.numeric(stats::quantile(x$value, 0.25, na.rm = TRUE, names = FALSE)),
      q75_value = as.numeric(stats::quantile(x$value, 0.75, na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$parameter, out$curve_class), , drop = FALSE]
}

class_change_audit_table <- function(previous_by_seed, current_by_seed) {
  if (is.null(previous_by_seed) || !nrow(previous_by_seed)) return(data.frame())
  keep_prev <- intersect(
    c("seed_id", "curve_class", "final_interpretation_class", "monotonicity_reliability",
      "ploidy_range", "slope_epsilon", "step_epsilon", "sign_sequence", "n_sign_changes"),
    names(previous_by_seed)
  )
  keep_curr <- intersect(
    c("seed_id", "curve_class", "final_interpretation_class", "monotonicity_reliability",
      "ploidy_range", "net_ploidy_change", "slope_epsilon", "step_epsilon",
      "flat_range_threshold", "max_positive_step", "max_negative_step",
      "fraction_positive_steps", "fraction_negative_steps", "fraction_zero_steps",
      "low_amplitude_curve", "terminal_plateau_for_class", "sign_sequence", "n_sign_changes", "classification_rule_version"),
    names(current_by_seed)
  )
  prev <- previous_by_seed[, keep_prev, drop = FALSE]
  curr <- current_by_seed[, keep_curr, drop = FALSE]
  names(prev)[names(prev) != "seed_id"] <- paste0("old_", names(prev)[names(prev) != "seed_id"])
  names(curr)[names(curr) != "seed_id"] <- paste0("new_", names(curr)[names(curr) != "seed_id"])
  out <- merge(prev, curr, by = "seed_id", all = TRUE, sort = FALSE)
  if ("old_curve_class" %in% names(out) && "new_curve_class" %in% names(out)) {
    out$class_changed <- !identical(out$old_curve_class, out$new_curve_class) &
      (is.na(out$old_curve_class) | is.na(out$new_curve_class) | out$old_curve_class != out$new_curve_class)
  }
  if ("old_final_interpretation_class" %in% names(out) && "new_final_interpretation_class" %in% names(out)) {
    out$final_interpretation_changed <- is.na(out$old_final_interpretation_class) |
      is.na(out$new_final_interpretation_class) |
      out$old_final_interpretation_class != out$new_final_interpretation_class
  }
  out[order(seed_number(out$seed_id)), , drop = FALSE]
}

curve_palette <- function(classes) {
  base <- c(
    approximately_flat = "#999999",
    monotone_increasing = "#0072B2",
    monotone_decreasing = "#D55E00",
    single_transition_increase_then_plateau = "#56B4E9",
    single_transition_decrease_then_plateau = "#E69F00",
    u_shaped = "#009E73",
    inverted_u_shaped = "#CC79A7",
    complex_nonmonotone = "#000000"
  )
  missing <- setdiff(classes, names(base))
  if (length(missing)) {
    extra <- grDevices::rainbow(length(missing), s = 0.65, v = 0.75)
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[classes]
}

save_plot_pair_in_fig_dir <- function(name, fig_dir, expr, width = 8, height = 5) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(fig_dir, paste0(name, ".pdf"))
  png_path <- file.path(fig_dir, paste0(name, ".png"))
  expr <- substitute(expr)
  grDevices::pdf(pdf_path, width = width, height = height, onefile = FALSE)
  eval(expr, envir = parent.frame())
  grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = 160)
  eval(expr, envir = parent.frame())
  grDevices::dev.off()
  c(pdf = pdf_path, png = png_path)
}

save_plot_pair <- function(name, out_dir, expr, width = 8, height = 5) {
  save_plot_pair_in_fig_dir(name, file.path(out_dir, "figures"), expr, width = width, height = height)
}

plot_all_curves <- function(curves, by_seed, out_dir) {
  d <- merge(curves, by_seed[, c("seed_id", "curve_class")], by = "seed_id", all.x = TRUE, sort = FALSE)
  classes <- sort(unique(d$curve_class))
  cols <- curve_palette(classes)
  save_plot_pair("fixed_o2_dominant_ploidy_all_seed_curves_by_class", out_dir, {
    op <- par(mar = c(4.5, 4.5, 1, 1))
    on.exit(par(op), add = TRUE)
    plot(range(d$O2_pct), range(d$dominant_mean_ploidy, na.rm = TRUE), type = "n",
         xlab = "Fixed O2 (%)", ylab = "Dominant mean ploidy")
    for (seed in unique(d$seed_id[order(seed_number(d$seed_id))])) {
      z <- d[d$seed_id == seed, , drop = FALSE]
      z <- z[order(z$O2_pct), , drop = FALSE]
      lines(z$O2_pct, z$dominant_mean_ploidy, col = grDevices::adjustcolor(cols[[z$curve_class[[1L]]]], alpha.f = 0.18), lwd = 0.7)
    }
    legend("topright", legend = classes, col = cols, lwd = 2, bty = "n", cex = 0.75)
  }, width = 9, height = 5.5)
}

plot_class_seed_overlays <- function(curves,
                                     by_seed,
                                     out_dir,
                                     stable_only = FALSE,
                                     stable_reliability = "reliable") {
  classes <- sort(unique(by_seed$curve_class))
  seed_tab <- by_seed
  version <- "all"
  version_label <- "all seeds"
  if (stable_only) {
    seed_tab <- seed_tab[seed_tab$monotonicity_reliability %in% stable_reliability, , drop = FALSE]
    version <- "stable_reliable"
    version_label <- "stable/reliable seeds"
  }

  keep_cols <- intersect(c("seed_id", "curve_class", "monotonicity_reliability"), names(seed_tab))
  d <- merge(curves, seed_tab[, keep_cols, drop = FALSE], by = "seed_id", all.x = FALSE, sort = FALSE)
  cols <- curve_palette(classes)
  fig_dir <- file.path(out_dir, "figures", "all_seed_curves_by_class")
  x_range <- range(curves$O2_pct, na.rm = TRUE)

  empty_panel <- function(cls, ylab) {
    plot(x_range, c(0, 1), type = "n", axes = FALSE,
         xlab = "O2 (%)", ylab = ylab, main = paste(cls, version_label, "n=0"))
    box()
    text(mean(x_range), 0.5, "No stable/reliable seeds", cex = 0.9, col = "grey35")
  }

  base_line_style <- function(n_seed) {
    if (n_seed >= 200) return(list(alpha = 0.10, lwd = 0.30))
    if (n_seed >= 100) return(list(alpha = 0.12, lwd = 0.32))
    if (n_seed >= 50) return(list(alpha = 0.18, lwd = 0.38))
    if (n_seed >= 25) return(list(alpha = 0.25, lwd = 0.50))
    if (n_seed >= 10) return(list(alpha = 0.35, lwd = 0.65))
    if (n_seed >= 5) return(list(alpha = 0.48, lwd = 0.85))
    list(alpha = 0.65, lwd = 1.05)
  }

  curve_overlap_density <- function(z,
                                    y_col,
                                    y_range,
                                    transform = identity,
                                    n_y_bins = 24L) {
    y <- transform(z[[y_col]])
    y_range <- transform(y_range)
    ok <- is.finite(z$O2_pct) & is.finite(y)
    if (!any(ok) || !all(is.finite(y_range))) return(NA_real_)

    y_min <- min(y_range, na.rm = TRUE)
    y_max <- max(y_range, na.rm = TRUE)
    if (!is.finite(y_min) || !is.finite(y_max) || y_min == y_max) return(NA_real_)

    breaks <- seq(y_min, y_max, length.out = n_y_bins + 1L)
    x_key <- sprintf("%.8f", z$O2_pct[ok])
    y_split <- split(y[ok], x_key)
    density <- vapply(y_split, function(vals) {
      bins <- cut(vals, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      bins <- bins[is.finite(bins)]
      n_bins <- length(unique(bins))
      if (!n_bins) return(NA_real_)
      length(vals) / n_bins
    }, numeric(1L))

    stats::median(density[is.finite(density)], na.rm = TRUE)
  }

  line_style <- function(z, y_col, y_range, transform = identity) {
    n_seed <- length(unique(z$seed_id))
    style <- base_line_style(n_seed)
    density <- curve_overlap_density(z, y_col, y_range, transform = transform)
    if (!is.finite(density) || density <= 0) return(style)

    density_target <- 12
    visibility_boost <- max(1, sqrt(density_target / density))
    style$alpha <- min(0.75, style$alpha * visibility_boost)
    style$lwd <- min(1.25, style$lwd * sqrt(visibility_boost))
    style
  }

  save_plot_pair_in_fig_dir(
    paste0("fixed_o2_all_seed_curves_by_class_", version),
    fig_dir,
    {
      op <- par(mfrow = c(length(classes), 2), mar = c(3, 4, 2, 1), oma = c(1, 0, 0, 0))
      on.exit(par(op), add = TRUE)
      for (cls in classes) {
        z <- d[d$curve_class == cls, , drop = FALSE]
        n_seed <- length(unique(z$seed_id))
        if (!n_seed) {
          empty_panel(cls, "Dominant ploidy")
          empty_panel(cls, "Spectral gap")
          next
        }

        z <- z[order(seed_number(z$seed_id), z$O2_pct), , drop = FALSE]
        seed_ids <- unique(z$seed_id[order(seed_number(z$seed_id))])

        y <- z$dominant_mean_ploidy
        y_range <- range(y[is.finite(y)], c(1.5, 2.0), na.rm = TRUE)
        style <- line_style(z, "dominant_mean_ploidy", y_range)
        line_col <- grDevices::adjustcolor(cols[[cls]], alpha.f = style$alpha)
        plot(range(z$O2_pct, na.rm = TRUE), y_range, type = "n",
             xlab = "O2 (%)", ylab = "Dominant ploidy",
             main = paste(cls, version_label, paste0("n=", n_seed)))
        abline(h = c(1.5, 2.0), lty = 2, col = "grey70")
        for (seed in seed_ids) {
          zz <- z[z$seed_id == seed, , drop = FALSE]
          zz <- zz[order(zz$O2_pct), , drop = FALSE]
          lines(zz$O2_pct, zz$dominant_mean_ploidy, col = line_col, lwd = style$lwd)
        }

        gap <- z$spectral_gap
        gap_range <- range(gap[is.finite(gap) & gap > 0], c(0.005, 0.01), na.rm = TRUE)
        style <- line_style(z, "spectral_gap", gap_range,
                            transform = function(x) log10(pmax(x, .Machine$double.eps)))
        line_col <- grDevices::adjustcolor(cols[[cls]], alpha.f = style$alpha)
        plot(range(z$O2_pct, na.rm = TRUE), gap_range, type = "n", log = "y",
             xlab = "O2 (%)", ylab = "Spectral gap",
             main = paste(cls, version_label, paste0("n=", n_seed)))
        abline(h = c(0.005, 0.01), lty = 2, col = "grey70")
        for (seed in seed_ids) {
          zz <- z[z$seed_id == seed, , drop = FALSE]
          zz <- zz[order(zz$O2_pct), , drop = FALSE]
          lines(zz$O2_pct, pmax(zz$spectral_gap, .Machine$double.eps), col = line_col, lwd = style$lwd)
        }
      }
    },
    width = 10,
    height = max(6, 2.1 * length(classes))
  )
}

plot_class_summary <- function(summary, out_dir) {
  classes <- sort(unique(summary$curve_class))
  cols <- curve_palette(classes)
  save_plot_pair("fixed_o2_dominant_ploidy_median_iqr_by_class", out_dir, {
    op <- par(mar = c(4.5, 4.5, 1, 1))
    on.exit(par(op), add = TRUE)
    plot(range(summary$O2_pct), range(c(summary$q25_dominant_mean_ploidy, summary$q75_dominant_mean_ploidy), na.rm = TRUE),
         type = "n", xlab = "Fixed O2 (%)", ylab = "Dominant mean ploidy")
    for (cls in classes) {
      z <- summary[summary$curve_class == cls, , drop = FALSE]
      z <- z[order(z$O2_pct), , drop = FALSE]
      polygon(c(z$O2_pct, rev(z$O2_pct)), c(z$q25_dominant_mean_ploidy, rev(z$q75_dominant_mean_ploidy)),
              col = grDevices::adjustcolor(cols[[cls]], alpha.f = 0.20), border = NA)
      lines(z$O2_pct, z$median_dominant_mean_ploidy, col = cols[[cls]], lwd = 2)
    }
    legend("topright", legend = classes, col = cols, lwd = 2, bty = "n", cex = 0.75)
  }, width = 9, height = 5.5)
}

plot_heatmap <- function(curves, by_seed, out_dir, metric, file_name, transform = identity, zlab = metric) {
  d <- merge(curves, by_seed[, c("seed_id", "curve_class", "objective")], by = "seed_id", all.x = TRUE, sort = FALSE)
  order_seed <- by_seed[order(by_seed$curve_class, by_seed$objective, by_seed$seed_number), "seed_id"]
  o2 <- sort(unique(d$O2_pct))
  mat <- matrix(NA_real_, nrow = length(order_seed), ncol = length(o2), dimnames = list(order_seed, o2_key(o2)))
  for (i in seq_len(nrow(d))) {
    mat[d$seed_id[[i]], o2_key(d$O2_pct[[i]])] <- suppressWarnings(as.numeric(d[[metric]][[i]]))
  }
  mat <- transform(mat)
  save_plot_pair(file_name, out_dir, {
    op <- par(mar = c(4.5, 4.5, 1, 5))
    on.exit(par(op), add = TRUE)
    z <- t(mat[nrow(mat):1L, , drop = FALSE])
    image(x = o2, y = seq_len(nrow(mat)), z = z, xlab = "Fixed O2 (%)",
          ylab = "Seeds ordered by class and objective", col = grDevices::hcl.colors(80, "Viridis"))
    box()
    mtext(zlab, side = 4, line = 3)
  }, width = 8, height = 7)
}

plot_mode_heatmap <- function(curves, by_seed, out_dir) {
  d <- merge(curves, by_seed[, c("seed_id", "curve_class", "objective")], by = "seed_id", all.x = TRUE, sort = FALSE)
  order_seed <- by_seed[order(by_seed$curve_class, by_seed$objective, by_seed$seed_number), "seed_id"]
  o2 <- sort(unique(d$O2_pct))
  mat <- matrix(NA_real_, nrow = length(order_seed), ncol = length(o2), dimnames = list(order_seed, o2_key(o2)))
  vals <- c(mode2 = 1, mode1 = 2)
  for (i in seq_len(nrow(d))) {
    mat[d$seed_id[[i]], o2_key(d$O2_pct[[i]])] <- vals[[as.character(d$mode_label[[i]])]]
  }
  save_plot_pair("fixed_o2_mode_label_heatmap_ordered_by_class", out_dir, {
    op <- par(mar = c(4.5, 4.5, 1, 4))
    on.exit(par(op), add = TRUE)
    z <- t(mat[nrow(mat):1L, , drop = FALSE])
    image(x = o2, y = seq_len(nrow(mat)), z = z, xlab = "Fixed O2 (%)",
          ylab = "Seeds ordered by class and objective", col = c("#D55E00", "#0072B2"), breaks = c(0.5, 1.5, 2.5))
    box()
    legend("right", inset = -0.18, xpd = TRUE, legend = c("mode2", "mode1"), fill = c("#D55E00", "#0072B2"), bty = "n")
  }, width = 8, height = 7)
}

plot_gap_scatter <- function(by_seed, out_dir) {
  classes <- sort(unique(by_seed$curve_class))
  cols <- curve_palette(classes)
  save_plot_pair("fixed_o2_min_spectral_gap_vs_ploidy_range", out_dir, {
    op <- par(mar = c(4.5, 4.5, 1, 1))
    on.exit(par(op), add = TRUE)
    plot(by_seed$ploidy_range, by_seed$min_spectral_gap, log = "y",
         pch = 19, col = grDevices::adjustcolor(cols[by_seed$curve_class], alpha.f = 0.75),
         xlab = "Dominant ploidy range across O2", ylab = "Minimum spectral gap")
    abline(h = c(0.005, 0.01), lty = 2, col = "grey40")
    legend("topright", legend = classes, col = cols, pch = 19, bty = "n", cex = 0.75)
  }, width = 8, height = 5.5)
}

plot_representatives <- function(curves, reps, out_dir) {
  reps <- reps[order(reps$curve_class), , drop = FALSE]
  if (!nrow(reps)) return(character())
  save_plot_pair("fixed_o2_representative_curves_by_class", out_dir, {
    op <- par(mfrow = c(nrow(reps), 2), mar = c(3, 4, 2, 1), oma = c(1, 0, 0, 0))
    on.exit(par(op), add = TRUE)
    for (i in seq_len(nrow(reps))) {
      seed <- reps$representative_seed_id[[i]]
      cls <- reps$curve_class[[i]]
      z <- curves[curves$seed_id == seed, , drop = FALSE]
      z <- z[order(z$O2_pct), , drop = FALSE]
      plot(z$O2_pct, z$dominant_mean_ploidy, type = "l", lwd = 2,
           xlab = "O2 (%)", ylab = "Dominant ploidy", main = paste(cls, seed))
      abline(h = c(1.5, 2.0), lty = 2, col = "grey60")
      plot(z$O2_pct, z$spectral_gap, type = "l", lwd = 2, log = "y",
           xlab = "O2 (%)", ylab = "Spectral gap", main = paste(cls, seed))
      abline(h = c(0.005, 0.01), lty = 2, col = "grey60")
    }
  }, width = 10, height = max(6, 2.1 * nrow(reps)))
}

run_validation <- function(out_dir,
                           flat_range_threshold,
                           step_epsilon_abs,
                           step_epsilon_fraction,
                           reverse_fraction_tolerance,
                           plateau_min_points) {
  helper_rows <- run_curve_classification_validation()
  mk <- function(seed, y, gap = 0.02) {
    o2 <- seq_along(y)
    data.frame(seed_id = seed, O2_pct = o2, dominant_mean_ploidy = y, spectral_gap = rep(gap, length(y)), stringsAsFactors = FALSE)
  }
  gap_case <- mk("seed_gap", seq(1, 3, length.out = 201), gap = c(rep(0.001, 101), rep(0.02, 100)))
  gap_case$mode_label <- ifelse(gap_case$dominant_mean_ploidy >= 2, "mode1", "mode2")
  gap_case$dominant_growth_rate <- 0
  gap_res <- classify_one_seed(
    gap_case,
    reporting_o2 = numeric(),
    gap_low = 0.005,
    gap_caution = 0.01,
    unreliable_fraction = 0.25,
    caution_fraction = 0.10,
    flat_range_threshold = flat_range_threshold,
    step_epsilon_abs = step_epsilon_abs,
    step_epsilon_fraction = step_epsilon_fraction,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    plateau_min_points = plateau_min_points
  )
  gap_row <- data.frame(
    test_case = "small_gap_flag",
    expected_class = "unreliable_small_gap",
    observed_class = gap_res$monotonicity_reliability,
    passed = identical("unreliable_small_gap", gap_res$monotonicity_reliability),
    stringsAsFactors = FALSE
  )
  out <- rbind(helper_rows, gap_row)
  write_tsv(out, file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_validation.tsv"))
  if (!all(out$passed)) stop("Validation failed: ", paste(out$test_case[!out$passed], collapse = ", "))
  out
}

write_summary_report <- function(out_dir, args, run_args, class_counts, reliability_counts, by_seed, paths) {
  report_dir <- file.path(out_dir, "report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(report_dir, "fixed_o2_ploidy_monotonicity_summary.md")
  fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
  class_lines <- paste0(
    "- `", class_counts$curve_class, "`: ",
    class_counts$n_seed, " seeds (", fmt_pct(class_counts$fraction_seed), ")"
  )
  rel_tab <- as.data.frame(table(by_seed$monotonicity_reliability), stringsAsFactors = FALSE)
  names(rel_tab) <- c("reliability", "n_seed")
  rel_tab$fraction_seed <- rel_tab$n_seed / sum(rel_tab$n_seed)
  reliability_lines <- paste0(
    "- `", rel_tab$reliability, "`: ",
    rel_tab$n_seed, " seeds (", fmt_pct(rel_tab$fraction_seed), ")"
  )
  txt <- c(
    "# Fixed-O2 Ploidy Monotonicity Summary",
    "",
    "Each O2 point in this analysis was analytically evaluated from the fixed-O2 generator and dominant eigenvector. The monotonicity class is inferred numerically from the dense O2 grid; this is not stochastic simulation.",
    "",
    "## Run Arguments",
    "",
    paste0("- `", run_args$argument, "`: `", run_args$value, "`"),
    "",
    "## Output Checks",
    "",
    paste0("- Curve rows: ", nrow(read_tsv(paths$curves))),
    paste0("- Seed classification rows: ", nrow(by_seed)),
    paste0("- Expected dense rows: ", as.integer(run_args$value[run_args$argument == "expected_curve_rows"])),
    "",
    "## Curve Classes",
    "",
    class_lines,
    "",
    "## Spectral-Gap Reliability",
    "",
    reliability_lines,
    "",
    "Small spectral gaps indicate that the dominant eigenmode is weakly separated from the next mode. Non-monotone curves localized to small-gap regions should be interpreted cautiously.",
    "",
    "## Main Tables",
    "",
    paste0("- Curves: `", paths$curves, "`"),
    paste0("- Seed classification: `", paths$by_seed, "`"),
    paste0("- Mode crosswalk: `", paths$crosswalk, "`"),
    paste0("- Class counts: `", paths$class_counts, "`"),
    paste0("- Class-change audit: `", paths$class_change_audit, "`"),
    "",
    "## Figures",
    "",
    paste0("- `", list.files(file.path(out_dir, "figures"), pattern = "\\.(pdf|png)$", full.names = FALSE, recursive = TRUE), "`")
  )
  writeLines(txt, path)
  path
}

generate_outputs <- function(args = parse_args()) {
  run_dir <- normalizePath(path.expand(args$run_dir %||% default_run_dir()), mustWork = FALSE)
  out_dir <- normalizePath(path.expand(args$out_dir %||% default_out_dir()), mustWork = FALSE)
  o2_grid <- sort(unique(as_num_vec(args$o2_grid, seq(0, 5, by = 0.025))))
  reporting_o2 <- sort(unique(as_num_vec(args$reporting_o2, c(0, 0.1, 0.5, 1, 2, 5))))
  n_workers <- as_int(args$n_workers, 8L)
  max_seeds <- as_int(args$max_seeds, NA_integer_)
  gap_low <- as_num(args$gap_low, 0.005)
  gap_caution <- as_num(args$gap_caution, 0.01)
  unreliable_fraction <- as_num(args$unreliable_fraction, 0.25)
  caution_fraction <- as_num(args$caution_fraction, 0.10)
  flat_range_threshold <- as_num(args$flat_range_threshold, 0.05)
  step_epsilon_abs <- as_num(args$step_epsilon_abs, 1e-6)
  step_epsilon_fraction <- as_num(args$step_epsilon_fraction, 1e-4)
  reverse_fraction_tolerance <- as_num(args$reverse_fraction_tolerance, 0.05)
  plateau_min_points <- as_int(args$plateau_min_points, 3L)
  overwrite <- as_bool(args$overwrite, TRUE)
  generate_figures <- as_bool(args$generate_figures, TRUE)
  validate <- as_bool(args$run_validation, TRUE)
  if (!dir.exists(run_dir)) stop("Input run_dir does not exist: ", run_dir)
  if (!length(o2_grid)) stop("o2_grid must contain at least one finite value.")
  if (!length(reporting_o2)) stop("reporting_o2 must contain at least one finite value.")
  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "report"), recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    curves = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_curves.tsv"),
    by_seed = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv"),
    crosswalk = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_mode_crosswalk.tsv"),
    class_counts = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_counts.tsv"),
    class_by_mode = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_by_reporting_o2_mode.tsv"),
    class_curves = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_curve_summary.tsv"),
    representatives = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_representative_seeds.tsv"),
    objective_rank = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_by_objective_rank.tsv"),
    parameter_differences = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_parameter_differences.tsv"),
    class_change_audit = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_change_audit.tsv"),
    run_arguments = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_run_arguments.tsv")
  )
  previous_by_seed <- if (file.exists(paths$by_seed)) {
    tryCatch(read_tsv(paths$by_seed), error = function(e) data.frame())
  } else {
    data.frame()
  }
  if (!overwrite && all(file.exists(unlist(paths[c("curves", "by_seed", "crosswalk")])))) {
    message("Reusing existing monotonicity outputs under: ", out_dir)
    return(invisible(paths))
  }

  if (validate) {
    run_validation(
      out_dir = out_dir,
      flat_range_threshold = flat_range_threshold,
      step_epsilon_abs = step_epsilon_abs,
      step_epsilon_fraction = step_epsilon_fraction,
      reverse_fraction_tolerance = reverse_fraction_tolerance,
      plateau_min_points = plateau_min_points
    )
  }

  fixo2_env <- load_fixo2_env()
  inputs <- get("o2ipa_collect_seed_inputs", envir = fixo2_env, inherits = TRUE)(run_dir, objective_source = "auto")
  manifest <- inputs$manifest
  seeds <- manifest$seed_id[order(seed_number(manifest$seed_id))]
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seeds <- seeds[seq_len(min(max_seeds, length(seeds)))]
  }
  expected_rows <- length(seeds) * length(o2_grid)
  message("Computing fixed-O2 dense analytical grid: ", length(seeds), " seeds x ", length(o2_grid), " O2 values = ", expected_rows, " rows")
  curves_raw <- get("generate_fixo2_attractor_mode_table", envir = fixo2_env, inherits = TRUE)(
    run_dir = run_dir,
    o2_values = o2_grid,
    seed_ids = seeds,
    n_workers = n_workers
  )
  if (nrow(curves_raw) != expected_rows) {
    stop("Unexpected curve row count: observed ", nrow(curves_raw), ", expected ", expected_rows)
  }
  bad_status <- curves_raw$status[is.na(curves_raw$status) | curves_raw$status != "ok"]
  if (length(bad_status)) warning("Non-ok attractor rows: ", length(bad_status))

  by_seed <- do.call(rbind, lapply(split(curves_raw, curves_raw$seed_id), classify_one_seed,
                                   reporting_o2 = reporting_o2, gap_low = gap_low,
                                   gap_caution = gap_caution,
                                   unreliable_fraction = unreliable_fraction,
                                   caution_fraction = caution_fraction,
                                   flat_range_threshold = flat_range_threshold,
                                   step_epsilon_abs = step_epsilon_abs,
                                   step_epsilon_fraction = step_epsilon_fraction,
                                   reverse_fraction_tolerance = reverse_fraction_tolerance,
                                   plateau_min_points = plateau_min_points))
  by_seed <- merge(by_seed, manifest[, intersect(c("seed_id", "objective", "objective_source", "objective_total", "objective_data", "objective_burden", "objective_ploidy", "runtime", "convergence_status"), names(manifest)), drop = FALSE],
                   by = "seed_id", all.x = TRUE, sort = FALSE)
  by_seed <- by_seed[order(by_seed$seed_number), , drop = FALSE]

  curves <- add_curve_differences(curves_raw, by_seed, reporting_o2)
  class_counts <- class_count_table(by_seed)
  class_by_mode <- class_by_mode_table(by_seed, reporting_o2)
  class_curves <- curve_quantile_table(curves, by_seed)
  reps <- representative_seed_table(curves, by_seed)
  objective_rank <- objective_rank_class_table(by_seed)
  param_diff <- parameter_difference_table(by_seed, inputs$params_long)
  class_change_audit <- class_change_audit_table(previous_by_seed, by_seed)
  crosswalk_cols <- c(
    "seed_id", "seed_number", "curve_class", "final_interpretation_class",
    "monotonicity_reliability", "n_sign_changes", "ploidy_range", "net_ploidy_change",
    "step_epsilon", "flat_range_threshold", "fraction_positive_steps", "fraction_negative_steps",
    "low_amplitude_curve", "terminal_plateau_for_class", "min_spectral_gap",
    "fraction_o2_gap_below_0p005", "fraction_o2_gap_below_0p01",
    grep("^(mode_label|dominant_mean_ploidy|spectral_gap)_o2_", names(by_seed), value = TRUE)
  )
  crosswalk <- by_seed[, intersect(crosswalk_cols, names(by_seed)), drop = FALSE]

  run_args <- data.frame(
    argument = c(
      "run_dir", "out_dir", "script", "o2_grid", "n_o2", "reporting_o2",
      "n_seed", "expected_curve_rows", "n_workers", "max_seeds",
      "classification_rule_version", "flat_range_threshold", "step_epsilon_rule",
      "step_epsilon_abs", "step_epsilon_fraction", "reverse_fraction_tolerance",
      "plateau_min_points", "gap_low", "gap_caution", "unreliable_fraction",
      "caution_fraction", "analytical_method"
    ),
    value = c(
      run_dir, out_dir, normalizePath(file.path(SCRIPT_DIR, "fixed_o2_ploidy_monotonicity.R"), mustWork = FALSE),
      paste(format(o2_grid, scientific = FALSE, trim = TRUE), collapse = ","),
      as.character(length(o2_grid)),
      paste(format(reporting_o2, scientific = FALSE, trim = TRUE), collapse = ","),
      as.character(length(seeds)),
      as.character(expected_rows),
      as.character(n_workers),
      as.character(max_seeds),
      curve_classification_rule_version(),
      as.character(flat_range_threshold),
      "max(step_epsilon_abs, step_epsilon_fraction * ploidy_range)",
      as.character(step_epsilon_abs),
      as.character(step_epsilon_fraction),
      as.character(reverse_fraction_tolerance),
      as.character(plateau_min_points),
      as.character(gap_low),
      as.character(gap_caution),
      as.character(unreliable_fraction),
      as.character(caution_fraction),
      "fixed-O2 generator dominant eigenvector; analytical grid evaluation, not stochastic simulation"
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(curves, paths$curves)
  write_tsv(by_seed, paths$by_seed)
  write_tsv(crosswalk, paths$crosswalk)
  write_tsv(class_counts, paths$class_counts)
  write_tsv(class_by_mode, paths$class_by_mode)
  write_tsv(class_curves, paths$class_curves)
  write_tsv(reps, paths$representatives)
  write_tsv(objective_rank, paths$objective_rank)
  write_tsv(param_diff, paths$parameter_differences)
  write_tsv(class_change_audit, paths$class_change_audit)
  write_tsv(run_args, paths$run_arguments)

  if (generate_figures) {
    plot_all_curves(curves, by_seed, out_dir)
    plot_class_summary(class_curves, out_dir)
    plot_heatmap(curves, by_seed, out_dir, "dominant_mean_ploidy",
                 "fixed_o2_dominant_ploidy_heatmap_ordered_by_class", zlab = "Dominant mean ploidy")
    plot_mode_heatmap(curves, by_seed, out_dir)
    plot_heatmap(curves, by_seed, out_dir, "spectral_gap",
                 "fixed_o2_spectral_gap_heatmap_ordered_by_class",
                 transform = function(x) log10(pmax(x, .Machine$double.eps)), zlab = "log10 spectral gap")
    plot_gap_scatter(by_seed, out_dir)
    plot_representatives(curves, reps, out_dir)
    plot_class_seed_overlays(curves, by_seed, out_dir, stable_only = FALSE)
    plot_class_seed_overlays(curves, by_seed, out_dir, stable_only = TRUE)
  }

  write_summary_report(out_dir, args, run_args, class_counts, table(by_seed$monotonicity_reliability), by_seed, paths)
  message("Completed fixed-O2 ploidy monotonicity analysis: ", out_dir)
  invisible(paths)
}

if (identical(environment(), globalenv())) {
  generate_outputs(parse_args())
}
