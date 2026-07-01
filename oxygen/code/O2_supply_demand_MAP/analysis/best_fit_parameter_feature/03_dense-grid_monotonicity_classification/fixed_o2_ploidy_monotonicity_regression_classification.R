#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
INTEGRATED_WORKFLOW_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
ANALYSIS_DIR <- normalizePath(file.path(INTEGRATED_WORKFLOW_DIR, ".."), mustWork = FALSE)
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

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

seed_number <- function(seed_id) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

default_input_dir <- function() {
  file.path(
    REPO_ROOT, "oxygen", "results", "analysis", "best_fit_parameter_feature",
    "03_dense-grid_monotonicity_classification", "monotonicity_classification",
    "dense-grid_monotonicity_classification"
  )
}

default_out_dir <- function() {
  file.path(
    REPO_ROOT, "oxygen", "results", "analysis", "best_fit_parameter_feature",
    "03_dense-grid_monotonicity_classification", "monotonicity_classification",
    "dense-grid_monotonicity_regression_classification"
  )
}

source(file.path(INTEGRATED_WORKFLOW_DIR, "util", "curve_classification_utils.R"))

rename_if_present <- function(x, old, new) {
  hit <- match(old, names(x), nomatch = 0L)
  if (hit > 0L) names(x)[[hit]] <- new
  x
}

safe_fraction <- function(n, d) {
  if (!is.finite(d) || d <= 0) return(NA_real_)
  n / d
}

class_count_table <- function(by_seed) {
  tab <- as.data.frame(table(by_seed$smooth_curve_class), stringsAsFactors = FALSE)
  names(tab) <- c("smooth_curve_class", "n_seed")
  tab$fraction_seed <- tab$n_seed / sum(tab$n_seed)
  tab[order(-tab$n_seed, tab$smooth_curve_class), , drop = FALSE]
}

class_change_audit_table <- function(by_seed) {
  tab <- as.data.frame(table(by_seed$pointwise_curve_class, by_seed$smooth_curve_class), stringsAsFactors = FALSE)
  names(tab) <- c("pointwise_curve_class", "smooth_curve_class", "n_seed")
  tab <- tab[tab$n_seed > 0, , drop = FALSE]
  pointwise_totals <- aggregate(n_seed ~ pointwise_curve_class, tab, sum)
  names(pointwise_totals)[[2]] <- "pointwise_n_seed"
  smooth_totals <- aggregate(n_seed ~ smooth_curve_class, tab, sum)
  names(smooth_totals)[[2]] <- "smooth_n_seed"
  tab <- merge(tab, pointwise_totals, by = "pointwise_curve_class", all.x = TRUE, sort = FALSE)
  tab <- merge(tab, smooth_totals, by = "smooth_curve_class", all.x = TRUE, sort = FALSE)
  tab$fraction_of_pointwise_class <- mapply(safe_fraction, tab$n_seed, tab$pointwise_n_seed)
  tab$fraction_of_smooth_class <- mapply(safe_fraction, tab$n_seed, tab$smooth_n_seed)
  tab <- tab[, c(
    "pointwise_curve_class", "smooth_curve_class", "n_seed",
    "pointwise_n_seed", "smooth_n_seed",
    "fraction_of_pointwise_class", "fraction_of_smooth_class"
  ), drop = FALSE]
  tab[order(-tab$n_seed, tab$pointwise_curve_class, tab$smooth_curve_class), , drop = FALSE]
}

classify_all_curves <- function(curves,
                                flat_range_threshold,
                                step_epsilon_abs,
                                step_epsilon_fraction,
                                reverse_fraction_tolerance,
                                smooth_span,
                                smooth_degree,
                                smooth_family,
                                min_segment_span_fraction,
                                min_segment_amplitude_abs,
                                min_segment_amplitude_fraction,
                                min_segment_points,
                                terminal_plateau_span_fraction,
                                terminal_plateau_amplitude_fraction) {
  groups <- split(curves, curves$seed_id)
  summaries <- vector("list", length(groups))
  diffs <- vector("list", length(groups))
  segments <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    curve <- groups[[i]]
    curve <- curve[order(curve$O2_pct), , drop = FALSE]
    res <- classify_o2_ploidy_curve_smooth(
      curve,
      value_col = "dominant_mean_ploidy",
      x_col = "O2_pct",
      id_col = "seed_id",
      flat_range_threshold = flat_range_threshold,
      step_epsilon_abs = step_epsilon_abs,
      step_epsilon_fraction = step_epsilon_fraction,
      reverse_fraction_tolerance = reverse_fraction_tolerance,
      smooth_span = smooth_span,
      smooth_degree = smooth_degree,
      smooth_family = smooth_family,
      min_segment_span_fraction = min_segment_span_fraction,
      min_segment_amplitude_abs = min_segment_amplitude_abs,
      min_segment_amplitude_fraction = min_segment_amplitude_fraction,
      min_segment_points = min_segment_points,
      terminal_plateau_span_fraction = terminal_plateau_span_fraction,
      terminal_plateau_amplitude_fraction = terminal_plateau_amplitude_fraction
    )
    sm <- res$summary
    sm$seed_id <- curve$seed_id[[1L]]
    sm$seed_number <- seed_number(curve$seed_id[[1L]])
    summaries[[i]] <- sm[, c("seed_id", "seed_number", setdiff(names(sm), c("seed_id", "seed_number"))), drop = FALSE]
    diffs[[i]] <- res$differences
    if (nrow(res$segments)) {
      sg <- res$segments
      sg$seed_id <- curve$seed_id[[1L]]
      sg$seed_number <- seed_number(curve$seed_id[[1L]])
      segments[[i]] <- sg[, c("seed_id", "seed_number", setdiff(names(sg), c("seed_id", "seed_number"))), drop = FALSE]
    }
  }
  list(
    summary = do.call(rbind, summaries),
    differences = do.call(rbind, diffs),
    segments = if (length(Filter(NROW, segments))) do.call(rbind, segments) else data.frame()
  )
}

build_by_seed_table <- function(pointwise_by_seed, smooth_summary) {
  raw_cols <- intersect(c(
    "seed_id", "seed_number", "curve_class", "final_interpretation_class",
    "sign_sequence", "n_sign_changes", "ploidy_range", "net_ploidy_change",
    "step_epsilon", "classification_rule_version", "monotonicity_reliability",
    "min_spectral_gap", "median_spectral_gap", "fraction_o2_gap_below_0p005",
    "fraction_o2_gap_below_0p01", "objective", "objective_source",
    "objective_total", "objective_data", "objective_burden", "objective_ploidy",
    "runtime", "convergence_status"
  ), names(pointwise_by_seed))
  raw <- pointwise_by_seed[, raw_cols, drop = FALSE]
  raw <- rename_if_present(raw, "curve_class", "pointwise_curve_class")
  raw <- rename_if_present(raw, "final_interpretation_class", "pointwise_final_interpretation_class")
  raw <- rename_if_present(raw, "sign_sequence", "pointwise_sign_sequence")
  raw <- rename_if_present(raw, "n_sign_changes", "pointwise_n_sign_changes")
  raw <- rename_if_present(raw, "ploidy_range", "pointwise_ploidy_range")
  raw <- rename_if_present(raw, "net_ploidy_change", "pointwise_net_ploidy_change")
  raw <- rename_if_present(raw, "step_epsilon", "pointwise_step_epsilon")
  raw <- rename_if_present(raw, "classification_rule_version", "pointwise_classification_rule_version")

  sm <- smooth_summary
  sm <- rename_if_present(sm, "curve_class", "smooth_curve_class")
  sm <- rename_if_present(sm, "sign_sequence", "smooth_sign_sequence")
  sm <- rename_if_present(sm, "n_sign_changes", "smooth_n_sign_changes")
  sm <- rename_if_present(sm, "classification_rule_version", "smooth_classification_rule_version")
  sm <- rename_if_present(sm, "step_epsilon", "smooth_step_epsilon")

  out <- merge(raw, sm, by = c("seed_id", "seed_number"), all.x = TRUE, sort = FALSE)
  out$class_changed <- out$pointwise_curve_class != out$smooth_curve_class
  out$class_change <- paste(out$pointwise_curve_class, out$smooth_curve_class, sep = " -> ")
  out$curve_class <- out$smooth_curve_class
  out$classification_rule_version <- out$smooth_classification_rule_version
  out$final_interpretation_class <- ifelse(
    "monotonicity_reliability" %in% names(out) & out$monotonicity_reliability == "unreliable_small_gap",
    "unreliable_small_spectral_gap",
    out$smooth_curve_class
  )
  front <- intersect(c(
    "seed_id", "seed_number", "curve_class", "final_interpretation_class",
    "smooth_curve_class", "pointwise_curve_class", "class_changed", "class_change",
    "monotonicity_reliability", "smooth_sign_sequence", "pointwise_sign_sequence",
    "smooth_n_sign_changes", "pointwise_n_sign_changes", "fitted_ploidy_range",
    "pointwise_ploidy_range", "net_ploidy_change", "pointwise_net_ploidy_change",
    "smooth_step_epsilon", "pointwise_step_epsilon"
  ), names(out))
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]
  out[order(out$seed_number), , drop = FALSE]
}

build_curve_table <- function(curves, differences) {
  curves$.row_order <- seq_len(nrow(curves))
  keep <- differences[, intersect(c(
    "seed_id", "O2_pct", "fitted_value", "fitted_difference_next",
    "fitted_local_slope_sign", "step_epsilon"
  ), names(differences)), drop = FALSE]
  names(keep) <- sub("^fitted_value$", "smoothed_dominant_mean_ploidy", names(keep))
  names(keep) <- sub("^fitted_difference_next$", "smoothed_difference_next", names(keep))
  names(keep) <- sub("^fitted_local_slope_sign$", "smoothed_local_slope_sign", names(keep))
  names(keep) <- sub("^step_epsilon$", "smooth_step_epsilon", names(keep))
  out <- merge(curves, keep, by = c("seed_id", "O2_pct"), all.x = TRUE, sort = FALSE)
  out <- out[order(out$.row_order), , drop = FALSE]
  out$.row_order <- NULL
  out
}

summarize_by_class_curve <- function(curves, by_seed) {
  d <- merge(
    curves,
    by_seed[, c("seed_id", "smooth_curve_class"), drop = FALSE],
    by = "seed_id",
    all.x = TRUE,
    sort = FALSE
  )
  classes <- sort(unique(d$smooth_curve_class))
  rows <- list()
  for (cls in classes) {
    zc <- d[d$smooth_curve_class == cls, , drop = FALSE]
    for (o2 in sort(unique(zc$O2_pct))) {
      z <- zc$smoothed_dominant_mean_ploidy[zc$O2_pct == o2]
      rows[[length(rows) + 1L]] <- data.frame(
        smooth_curve_class = cls,
        O2_pct = o2,
        n_seed = length(z),
        q25 = stats::quantile(z, 0.25, na.rm = TRUE, names = FALSE),
        median = stats::median(z, na.rm = TRUE),
        q75 = stats::quantile(z, 0.75, na.rm = TRUE, names = FALSE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

plot_colors <- function(classes) {
  pal <- c(
    complex_nonmonotone = "#6A5ACD",
    inverted_u_shaped = "#D55E00",
    monotone_increasing = "#009E73",
    monotone_decreasing = "#0072B2",
    single_transition_increase_then_plateau = "#CC79A7",
    single_transition_decrease_then_plateau = "#56B4E9",
    u_shaped = "#E69F00",
    approximately_flat = "#666666",
    insufficient_data = "#999999"
  )
  missing <- setdiff(classes, names(pal))
  if (length(missing)) pal <- c(pal, stats::setNames(grDevices::rainbow(length(missing)), missing))
  pal[classes]
}

with_plot_pair <- function(out_dir, stub, code, width = 11, height = 8) {
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, "figures", paste0(stub, ".pdf"))
  png_path <- file.path(out_dir, "figures", paste0(stub, ".png"))
  expr <- substitute(code)
  env <- parent.frame()
  grDevices::pdf(pdf_path, width = width, height = height, onefile = TRUE)
  eval(expr, envir = env)
  grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = 180)
  eval(expr, envir = env)
  grDevices::dev.off()
  invisible(c(pdf = pdf_path, png = png_path))
}

adaptive_lwd <- function(n) {
  max(0.18, min(1.05, 2.8 / sqrt(max(n, 1))))
}

plot_all_smoothed_curves <- function(curves, by_seed, out_dir) {
  d <- merge(curves, by_seed[, c("seed_id", "smooth_curve_class"), drop = FALSE], by = "seed_id", all.x = TRUE, sort = FALSE)
  classes <- names(sort(table(by_seed$smooth_curve_class), decreasing = TRUE))
  cols <- plot_colors(classes)
  y_range <- range(d$smoothed_dominant_mean_ploidy, na.rm = TRUE)
  with_plot_pair(out_dir, "fixed_o2_regression_smoothed_all_seed_curves_by_class", {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    nr <- ceiling(length(classes) / 2)
    par(mfrow = c(nr, 2), mar = c(4, 4, 3, 1))
    for (cls in classes) {
      z <- d[d$smooth_curve_class == cls, , drop = FALSE]
      seeds <- unique(z$seed_id)
      plot(range(z$O2_pct, na.rm = TRUE), y_range, type = "n", xlab = "O2 (%)", ylab = "Smoothed dominant mean ploidy",
           main = paste0(cls, " (n=", length(seeds), ")"))
      lwd <- adaptive_lwd(length(seeds))
      col <- grDevices::adjustcolor(cols[[cls]], alpha.f = max(0.18, min(0.55, 10 / sqrt(max(length(seeds), 1)))))
      for (seed in seeds) {
        zz <- z[z$seed_id == seed, , drop = FALSE]
        lines(zz$O2_pct, zz$smoothed_dominant_mean_ploidy, col = col, lwd = lwd)
      }
      grid(col = "#DDDDDD")
    }
  })
}

plot_class_summary <- function(class_curves, out_dir) {
  classes <- unique(class_curves$smooth_curve_class)
  cols <- plot_colors(classes)
  y_range <- range(class_curves$q25, class_curves$q75, na.rm = TRUE)
  with_plot_pair(out_dir, "fixed_o2_regression_smoothed_median_iqr_by_class", {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    plot(range(class_curves$O2_pct, na.rm = TRUE), y_range, type = "n", xlab = "O2 (%)",
         ylab = "Smoothed dominant mean ploidy", main = "Regression classification: median and IQR by class")
    for (cls in classes) {
      z <- class_curves[class_curves$smooth_curve_class == cls, , drop = FALSE]
      ord <- order(z$O2_pct)
      z <- z[ord, , drop = FALSE]
      polygon(c(z$O2_pct, rev(z$O2_pct)), c(z$q25, rev(z$q75)),
              border = NA, col = grDevices::adjustcolor(cols[[cls]], alpha.f = 0.18))
      lines(z$O2_pct, z$median, col = cols[[cls]], lwd = 2)
    }
    legend("topleft", legend = classes, col = cols, lwd = 2, bty = "n", cex = 0.75)
    grid(col = "#DDDDDD")
  })
}

plot_class_transition_heatmap <- function(by_seed, out_dir) {
  raw <- sort(unique(by_seed$pointwise_curve_class))
  smooth <- sort(unique(by_seed$smooth_curve_class))
  mat <- table(
    factor(by_seed$pointwise_curve_class, levels = raw),
    factor(by_seed$smooth_curve_class, levels = smooth)
  )
  with_plot_pair(out_dir, "fixed_o2_regression_pointwise_vs_smooth_class_transition", width = 10, height = 8, {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mar = c(8, 9, 3, 1))
    image(seq_along(smooth), seq_along(raw), t(mat), axes = FALSE, xlab = "", ylab = "",
          col = grDevices::colorRampPalette(c("#FFFFFF", "#B7D7E8", "#2B6C8A"))(100),
          main = "Pointwise vs regression-smoothed class")
    axis(1, at = seq_along(smooth), labels = smooth, las = 2, cex.axis = 0.75)
    axis(2, at = seq_along(raw), labels = raw, las = 2, cex.axis = 0.75)
    mtext("Regression-smoothed class", side = 1, line = 6)
    mtext("Pointwise class", side = 2, line = 7)
    for (i in seq_along(raw)) {
      for (j in seq_along(smooth)) {
        val <- mat[i, j]
        if (val > 0) text(j, i, labels = val, cex = 0.8)
      }
    }
  })
}

plot_representative_raw_vs_smooth <- function(curves, by_seed, out_dir) {
  classes <- names(sort(table(by_seed$smooth_curve_class), decreasing = TRUE))
  with_plot_pair(out_dir, "fixed_o2_regression_raw_vs_smoothed_representative_curves", {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    nr <- ceiling(length(classes) / 2)
    par(mfrow = c(nr, 2), mar = c(4, 4, 3, 1))
    for (cls in classes) {
      z <- by_seed[by_seed$smooth_curve_class == cls, , drop = FALSE]
      target <- stats::median(z$fitted_ploidy_range, na.rm = TRUE)
      z$dist <- abs(z$fitted_ploidy_range - target)
      rep_seed <- z$seed_id[order(z$dist, z$seed_number)][[1L]]
      zz <- curves[curves$seed_id == rep_seed, , drop = FALSE]
      yr <- range(zz$dominant_mean_ploidy, zz$smoothed_dominant_mean_ploidy, na.rm = TRUE)
      plot(zz$O2_pct, zz$dominant_mean_ploidy, type = "l", col = "#888888", lwd = 1,
           ylim = yr, xlab = "O2 (%)", ylab = "Dominant mean ploidy",
           main = paste0(cls, ": ", rep_seed))
      lines(zz$O2_pct, zz$smoothed_dominant_mean_ploidy, col = "#0072B2", lwd = 2)
      legend("topleft", legend = c("raw", "smoothed"), col = c("#888888", "#0072B2"), lwd = c(1, 2), bty = "n", cex = 0.8)
      grid(col = "#DDDDDD")
    }
  })
}

write_summary_report <- function(out_dir, run_args, class_counts, class_change_audit, by_seed, paths) {
  report_dir <- file.path(out_dir, "report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(report_dir, "fixed_o2_ploidy_monotonicity_regression_summary.md")
  fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
  class_lines <- paste0(
    "- `", class_counts$smooth_curve_class, "`: ",
    class_counts$n_seed, " seeds (", fmt_pct(class_counts$fraction_seed), ")"
  )
  changed_n <- sum(by_seed$class_changed, na.rm = TRUE)
  changed_frac <- changed_n / nrow(by_seed)
  top_changes <- class_change_audit[class_change_audit$pointwise_curve_class != class_change_audit$smooth_curve_class, , drop = FALSE]
  top_changes <- head(top_changes[order(-top_changes$n_seed), , drop = FALSE], 12)
  change_lines <- if (nrow(top_changes)) {
    paste0(
      "- `", top_changes$pointwise_curve_class, "` -> `", top_changes$smooth_curve_class,
      "`: ", top_changes$n_seed, " seeds"
    )
  } else {
    "- No pointwise-to-smoothed class changes."
  }
  txt <- c(
    "# Fixed-O2 Ploidy Regression-Smoothed Classification Summary",
    "",
    "This analysis reuses the existing fixed-O2 dominant-eigenvector curve table and changes only the curve-shape classification step.",
    "Each seed is classified from a robust loess-smoothed ploidy-vs-O2 curve, and only sign segments that persist across a minimum O2 span, point count, and ploidy amplitude are allowed to define the final shape class.",
    "",
    "## Run Arguments",
    "",
    paste0("- `", run_args$argument, "`: `", run_args$value, "`"),
    "",
    "## Smoothed Curve Classes",
    "",
    class_lines,
    "",
    "## Pointwise vs Smoothed Changes",
    "",
    paste0("- Changed seeds: ", changed_n, " / ", nrow(by_seed), " (", fmt_pct(changed_frac), ")"),
    "",
    change_lines,
    "",
    "## Main Tables",
    "",
    paste0("- Seed classification: `", paths$by_seed, "`"),
    paste0("- Smoothed curves: `", paths$curves, "`"),
    paste0("- Persistent segments: `", paths$segments, "`"),
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
  input_dir <- normalizePath(path.expand(args$input_dir %||% default_input_dir()), mustWork = FALSE)
  out_dir <- normalizePath(path.expand(args$out_dir %||% default_out_dir()), mustWork = FALSE)
  flat_range_threshold <- as_num(args$flat_range_threshold, 0.05)
  step_epsilon_abs <- as_num(args$step_epsilon_abs, 1e-6)
  step_epsilon_fraction <- as_num(args$step_epsilon_fraction, 1e-4)
  reverse_fraction_tolerance <- as_num(args$reverse_fraction_tolerance, 0.05)
  smooth_span <- as_num(args$smooth_span, 0.20)
  smooth_degree <- as_int(args$smooth_degree, 2L)
  smooth_family <- args$smooth_family %||% "symmetric"
  min_segment_span_fraction <- as_num(args$min_segment_span_fraction, 0.02)
  min_segment_amplitude_abs <- as_num(args$min_segment_amplitude_abs, 0.01)
  min_segment_amplitude_fraction <- as_num(args$min_segment_amplitude_fraction, 0.03)
  min_segment_points <- as_int(args$min_segment_points, 3L)
  terminal_plateau_span_fraction <- as_num(args$terminal_plateau_span_fraction, 0.10)
  terminal_plateau_amplitude_fraction <- as_num(args$terminal_plateau_amplitude_fraction, 0.03)
  run_validation <- as_bool(args$run_validation, TRUE)
  generate_figures <- as_bool(args$generate_figures, TRUE)

  input_paths <- list(
    curves = file.path(input_dir, "tables", "fixed_o2_ploidy_monotonicity_curves.tsv"),
    by_seed = file.path(input_dir, "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv")
  )
  missing <- unlist(input_paths)[!file.exists(unlist(input_paths))]
  if (length(missing)) stop("Missing input files: ", paste(missing, collapse = ", "))

  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "report"), recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    curves = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_curves.tsv"),
    by_seed = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv"),
    segments = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_segments.tsv"),
    class_counts = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_class_counts.tsv"),
    class_change_audit = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_class_change_audit.tsv"),
    class_curves = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_class_curve_summary.tsv"),
    validation = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_validation.tsv"),
    run_arguments = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_run_arguments.tsv")
  )

  if (run_validation) {
    validation <- run_smooth_curve_classification_validation()
    write_tsv(validation, paths$validation)
    if (!all(validation$passed)) {
      stop("Smooth classification validation failed: ", paste(validation$test_case[!validation$passed], collapse = ", "))
    }
  }

  message("Reading pointwise dense-grid outputs from: ", input_dir)
  curves <- read_tsv(input_paths$curves)
  by_seed_pointwise <- read_tsv(input_paths$by_seed)
  if (!all(c("seed_id", "O2_pct", "dominant_mean_ploidy") %in% names(curves))) {
    stop("Curve table must contain seed_id, O2_pct, and dominant_mean_ploidy.")
  }
  message("Classifying ", length(unique(curves$seed_id)), " seeds with ", smooth_curve_classification_rule_version())
  class_info <- classify_all_curves(
    curves = curves,
    flat_range_threshold = flat_range_threshold,
    step_epsilon_abs = step_epsilon_abs,
    step_epsilon_fraction = step_epsilon_fraction,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    smooth_span = smooth_span,
    smooth_degree = smooth_degree,
    smooth_family = smooth_family,
    min_segment_span_fraction = min_segment_span_fraction,
    min_segment_amplitude_abs = min_segment_amplitude_abs,
    min_segment_amplitude_fraction = min_segment_amplitude_fraction,
    min_segment_points = min_segment_points,
    terminal_plateau_span_fraction = terminal_plateau_span_fraction,
    terminal_plateau_amplitude_fraction = terminal_plateau_amplitude_fraction
  )
  by_seed <- build_by_seed_table(by_seed_pointwise, class_info$summary)
  curve_table <- build_curve_table(curves, class_info$differences)
  class_counts <- class_count_table(by_seed)
  class_change_audit <- class_change_audit_table(by_seed)
  class_curves <- summarize_by_class_curve(curve_table, by_seed)

  run_args <- data.frame(
    argument = c(
      "input_dir", "out_dir", "script", "source_curve_table", "source_by_seed_table",
      "n_seed", "n_curve_row", "classification_rule_version", "smoothing_method",
      "smooth_span", "smooth_degree", "smooth_family", "flat_range_threshold",
      "step_epsilon_rule", "step_epsilon_abs", "step_epsilon_fraction",
      "reverse_fraction_tolerance", "min_segment_span_fraction",
      "min_segment_amplitude_abs", "min_segment_amplitude_fraction",
      "min_segment_points", "terminal_plateau_span_fraction",
      "terminal_plateau_amplitude_fraction"
    ),
    value = c(
      input_dir, out_dir,
      normalizePath(file.path(SCRIPT_DIR, "fixed_o2_ploidy_monotonicity_regression_classification.R"), mustWork = FALSE),
      input_paths$curves, input_paths$by_seed,
      as.character(length(unique(curves$seed_id))),
      as.character(nrow(curves)),
      smooth_curve_classification_rule_version(),
      "stats::loess robust family=symmetric, with smooth.spline fallback",
      as.character(smooth_span),
      as.character(smooth_degree),
      as.character(smooth_family),
      as.character(flat_range_threshold),
      "max(step_epsilon_abs, step_epsilon_fraction * fitted_ploidy_range)",
      as.character(step_epsilon_abs),
      as.character(step_epsilon_fraction),
      as.character(reverse_fraction_tolerance),
      as.character(min_segment_span_fraction),
      as.character(min_segment_amplitude_abs),
      as.character(min_segment_amplitude_fraction),
      as.character(min_segment_points),
      as.character(terminal_plateau_span_fraction),
      as.character(terminal_plateau_amplitude_fraction)
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(curve_table, paths$curves)
  write_tsv(by_seed, paths$by_seed)
  write_tsv(class_info$segments, paths$segments)
  write_tsv(class_counts, paths$class_counts)
  write_tsv(class_change_audit, paths$class_change_audit)
  write_tsv(class_curves, paths$class_curves)
  write_tsv(run_args, paths$run_arguments)

  if (generate_figures) {
    plot_all_smoothed_curves(curve_table, by_seed, out_dir)
    plot_class_summary(class_curves, out_dir)
    plot_class_transition_heatmap(by_seed, out_dir)
    plot_representative_raw_vs_smooth(curve_table, by_seed, out_dir)
  }

  write_summary_report(out_dir, run_args, class_counts, class_change_audit, by_seed, paths)
  message("Completed regression-smoothed fixed-O2 monotonicity classification: ", out_dir)
  invisible(paths)
}

if (identical(environment(), globalenv())) {
  generate_outputs(parse_args())
}
