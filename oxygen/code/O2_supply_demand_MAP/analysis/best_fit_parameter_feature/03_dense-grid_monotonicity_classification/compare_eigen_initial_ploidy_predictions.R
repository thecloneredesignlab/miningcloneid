#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  normalizePath(sys.frames()[[1L]]$ofile %||% "", mustWork = FALSE)
}

script_dir <- function() {
  p <- script_path()
  if (nzchar(p)) dirname(p) else getwd()
}

repo_root <- function() {
  cur <- normalizePath(script_dir(), mustWork = FALSE)
  for (i in seq_len(12L)) {
    if (dir.exists(file.path(cur, "oxygen", "code", "O2_supply_demand_MAP"))) {
      return(normalizePath(cur, mustWork = FALSE))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not locate repository root from: ", script_dir(), call. = FALSE)
}

default_result_root <- function() {
  file.path(
    repo_root(), "oxygen", "results", "analysis", "best_fit_parameter_feature",
    "03_dense-grid_monotonicity_classification", "monotonicity_classification"
  )
}

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    if (!grepl("=", arg, fixed = TRUE)) {
      out[[arg]] <- TRUE
    } else {
      key <- sub("=.*$", "", arg)
      val <- sub("^[^=]*=", "", arg)
      out[[key]] <- val
    }
  }
  out
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) return(default)
  vals <- unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE)
  vals <- suppressWarnings(as.numeric(trimws(vals)))
  vals[is.finite(vals)]
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(as.character(x[[1L]]))) return(default)
  val <- suppressWarnings(as.integer(x[[1L]]))
  if (length(val) && is.finite(val)) val else default
}

format_num <- function(x, digits = 6L) {
  formatC(as.numeric(x), digits = digits, format = "fg", flag = "#")
}

require_data_table <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("This comparison script requires the R package 'data.table' for memory-efficient table joins.")
  }
  invisible(TRUE)
}

read_dt <- function(path, select = NULL) {
  if (!file.exists(path)) stop("Missing required table: ", path)
  data.table::fread(path, sep = "\t", header = TRUE, quote = "", data.table = TRUE, select = select)
}

write_dt <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

analysis_paths <- function(root, analysis_max_day = 1000L) {
  mono_dir <- file.path(root, "dense-grid_monotonicity_classification")
  init_dir <- file.path(root, "dense-grid_initial-ploidy_trajectory")
  compare_parent <- file.path(root, "compare")
  analysis_label <- paste0(as.integer(analysis_max_day), "D")
  out_dir <- file.path(compare_parent, analysis_label)
  list(
    root = root,
    mono_dir = mono_dir,
    init_dir = init_dir,
    compare_parent = compare_parent,
    analysis_max_day = as.integer(analysis_max_day),
    analysis_label = analysis_label,
    out_dir = out_dir,
    out_tables = file.path(out_dir, "tables"),
    out_figures = file.path(out_dir, "figures"),
    eigen_curves = file.path(mono_dir, "tables", "fixed_o2_ploidy_monotonicity_curves.tsv"),
    eigen_by_seed = file.path(mono_dir, "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv"),
    eigen_validation = file.path(mono_dir, "tables", "fixed_o2_ploidy_monotonicity_validation.tsv"),
    initial_selected = file.path(init_dir, "tables", "fixed_o2_initial_ploidy_selected_time_curves.tsv"),
    initial_class = file.path(init_dir, "tables", "fixed_o2_initial_ploidy_curve_class_by_seed_time.tsv"),
    initial_validation = file.path(init_dir, "tables", "validation.tsv"),
    initial_args = file.path(init_dir, "tables", "analysis_run_arguments.tsv")
  )
}

prefix_cols <- function(dt, prefix, except) {
  cols <- setdiff(names(dt), except)
  data.table::setnames(dt, cols, paste0(prefix, cols))
  dt
}

safe_range <- function(x, fallback = c(0, 1)) {
  x <- x[is.finite(x)]
  if (!length(x)) return(fallback)
  r <- range(x)
  if (!all(is.finite(r)) || r[[1L]] == r[[2L]]) {
    mid <- if (all(is.finite(r))) r[[1L]] else mean(fallback)
    return(mid + c(-0.5, 0.5))
  }
  r
}

save_plot_pair <- function(name, figure_dir, expr, width = 8, height = 5, res = 160) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(figure_dir, paste0(name, ".pdf"))
  png_path <- file.path(figure_dir, paste0(name, ".png"))
  expr <- substitute(expr)
  grDevices::pdf(pdf_path, width = width, height = height, onefile = FALSE)
  eval(expr, envir = parent.frame())
  grDevices::dev.off()
  png_args <- list(filename = png_path, width = width, height = height, units = "in", res = res)
  if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo"
  do.call(grDevices::png, png_args)
  eval(expr, envir = parent.frame())
  grDevices::dev.off()
  c(pdf = pdf_path, png = png_path)
}

metric_summary <- function(dt, value_col) {
  x <- dt[[value_col]]
  data.table::data.table(
    median = stats::median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q90 = as.numeric(stats::quantile(x, 0.90, na.rm = TRUE, names = FALSE)),
    q95 = as.numeric(stats::quantile(x, 0.95, na.rm = TRUE, names = FALSE)),
    max = suppressWarnings(max(x, na.rm = TRUE)),
    fraction_le_0p01 = mean(x <= 0.01, na.rm = TRUE),
    fraction_le_0p05 = mean(x <= 0.05, na.rm = TRUE)
  )
}

wide_metric_summary <- function(dt, by_cols) {
  pieces <- list()
  metrics <- c(
    abs_delta_2N_minus_eigen = "abs_delta_2N_minus_eigen",
    abs_delta_4N_minus_eigen = "abs_delta_4N_minus_eigen",
    abs_delta_4N_minus_2N = "abs_delta_4N_minus_2N"
  )
  for (metric_name in names(metrics)) {
    col <- metrics[[metric_name]]
    z <- dt[, metric_summary(.SD, col), by = by_cols]
    data.table::setnames(
      z,
      c("median", "mean", "q90", "q95", "max", "fraction_le_0p01", "fraction_le_0p05"),
      paste(metric_name, c("median", "mean", "q90", "q95", "max", "fraction_le_0p01", "fraction_le_0p05"), sep = "_")
    )
    pieces[[metric_name]] <- z
  }
  out <- Reduce(function(a, b) merge(a, b, by = by_cols, all = TRUE), pieces)
  n_tab <- dt[, .(
    n_seed_o2 = .N,
    n_seed = data.table::uniqueN(seed_id),
    n_o2 = data.table::uniqueN(O2_key),
    eigen_spectral_gap_median = stats::median(eigen_spectral_gap, na.rm = TRUE),
    eigen_spectral_gap_min = suppressWarnings(min(eigen_spectral_gap, na.rm = TRUE)),
    eigen_spectral_gap_q10 = as.numeric(stats::quantile(eigen_spectral_gap, 0.10, na.rm = TRUE, names = FALSE)),
    spearman_abs_2N_eigen_vs_gap = suppressWarnings(stats::cor(abs_delta_2N_minus_eigen, eigen_spectral_gap, method = "spearman", use = "complete.obs")),
    spearman_abs_4N_eigen_vs_gap = suppressWarnings(stats::cor(abs_delta_4N_minus_eigen, eigen_spectral_gap, method = "spearman", use = "complete.obs")),
    spearman_abs_4N_2N_vs_gap = suppressWarnings(stats::cor(abs_delta_4N_minus_2N, eigen_spectral_gap, method = "spearman", use = "complete.obs"))
  ), by = by_cols]
  merge(n_tab, out, by = by_cols, all = TRUE)
}

plot_convergence <- function(summary_by_time, figure_dir) {
  d <- summary_by_time[order(day)]
  save_plot_pair("convergence_to_eigen_by_time", figure_dir, {
    op <- par(mar = c(4.5, 4.7, 1.2, 1), las = 1)
    on.exit(par(op), add = TRUE)
    y <- c(
      d$abs_delta_2N_minus_eigen_median, d$abs_delta_4N_minus_eigen_median, d$abs_delta_4N_minus_2N_median,
      d$abs_delta_2N_minus_eigen_q95, d$abs_delta_4N_minus_eigen_q95, d$abs_delta_4N_minus_2N_q95
    )
    yr <- safe_range(y[is.finite(y)], c(0, 1))
    plot(d$day, d$abs_delta_2N_minus_eigen_median, type = "n", ylim = yr,
         xlab = "Prediction day", ylab = "Absolute mean-ploidy difference")
    lines(d$day, d$abs_delta_2N_minus_eigen_median, col = "#0072B2", lwd = 2)
    lines(d$day, d$abs_delta_4N_minus_eigen_median, col = "#D55E00", lwd = 2)
    lines(d$day, d$abs_delta_4N_minus_2N_median, col = "#009E73", lwd = 2)
    lines(d$day, d$abs_delta_2N_minus_eigen_q95, col = "#0072B2", lwd = 1.5, lty = 2)
    lines(d$day, d$abs_delta_4N_minus_eigen_q95, col = "#D55E00", lwd = 1.5, lty = 2)
    lines(d$day, d$abs_delta_4N_minus_2N_q95, col = "#009E73", lwd = 1.5, lty = 2)
    abline(h = c(0.01, 0.05), col = "grey70", lty = 3)
    legend("topright",
           legend = c("|2N - eigen| median", "|4N - eigen| median", "|4N - 2N| median", "q95"),
           col = c("#0072B2", "#D55E00", "#009E73", "grey30"),
           lty = c(1, 1, 1, 2), lwd = c(2, 2, 2, 1.5), bty = "n", cex = 0.85)
  }, width = 8.5, height = 5.2)
}

plot_seed_convergence_panels <- function(comp, summary_by_time, figure_dir) {
  seed_summary <- comp[, .(
    abs_delta_2N_minus_eigen_seed_median = stats::median(abs_delta_2N_minus_eigen, na.rm = TRUE),
    abs_delta_4N_minus_eigen_seed_median = stats::median(abs_delta_4N_minus_eigen, na.rm = TRUE),
    abs_delta_4N_minus_2N_seed_median = stats::median(abs_delta_4N_minus_2N, na.rm = TRUE)
  ), by = .(seed_id, seed_number, day)]
  seed_summary <- seed_summary[order(seed_number, day)]
  d <- summary_by_time[order(day)]

  panels <- list(
    list(
      title = "|2N - eigen|",
      seed_col = "abs_delta_2N_minus_eigen_seed_median",
      median_col = "abs_delta_2N_minus_eigen_median",
      mean_col = "abs_delta_2N_minus_eigen_mean",
      color = "#0072B2",
      ylab = "Median over O2"
    ),
    list(
      title = "|4N - eigen|",
      seed_col = "abs_delta_4N_minus_eigen_seed_median",
      median_col = "abs_delta_4N_minus_eigen_median",
      mean_col = "abs_delta_4N_minus_eigen_mean",
      color = "#D55E00",
      ylab = "Median over O2"
    ),
    list(
      title = "|4N - 2N|",
      seed_col = "abs_delta_4N_minus_2N_seed_median",
      median_col = "abs_delta_4N_minus_2N_median",
      mean_col = "abs_delta_4N_minus_2N_mean",
      color = "#D7191C",
      ylab = "Median over O2"
    )
  )

  save_plot_pair("convergence_to_eigen_by_time_seed_curves_panels", figure_dir, {
    op <- par(mfrow = c(3, 1), mar = c(1.2, 6.8, 0.8, 1), oma = c(3.4, 0.2, 0.4, 0), las = 1)
    on.exit(par(op), add = TRUE)
    for (panel_i in seq_along(panels)) {
      panel <- panels[[panel_i]]
      y <- c(seed_summary[[panel$seed_col]], d[[panel$median_col]])
      yr <- safe_range(y[is.finite(y)], c(0, 1))
      plot(range(d$day), yr, type = "n",
           xlab = "", ylab = panel$ylab, main = "",
           xaxt = if (panel_i == length(panels)) "s" else "n")
      mtext(panel$title, side = 2, line = 5.2, las = 3, font = 2, col = panel$color, cex = 0.9)
      abline(h = c(0.01, 0.05), col = "grey82", lty = 3)
      for (seed in unique(seed_summary$seed_id[order(seed_summary$seed_number)])) {
        z <- seed_summary[seed_id == seed][order(day)]
        lines(z$day, z[[panel$seed_col]],
              col = grDevices::adjustcolor("grey35", alpha.f = 0.18), lwd = 0.45)
      }
      lines(d$day, d[[panel$mean_col]], col = panel$color, lwd = 2.2, lty = 2)
      lines(d$day, d[[panel$median_col]], col = panel$color, lwd = 3.4)
      points(d$day, d[[panel$median_col]], col = panel$color, pch = 16, cex = 0.7)
      legend("topright",
             legend = c("per-seed median over O2", "overall median", "overall mean"),
             col = c(grDevices::adjustcolor("grey35", alpha.f = 0.45), panel$color, panel$color),
             lwd = c(0.8, 3.4, 2.2), lty = c(1, 1, 2), pch = c(NA, 16, NA),
             bty = "n", cex = 0.82)
    }
    mtext("Prediction day", side = 1, outer = TRUE, line = 2.0)
  }, width = 8.5, height = 10.5)
}

slugify <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", as.character(x)))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "na")
}

adaptive_seed_base_line_style <- function(n_seed) {
  if (n_seed >= 200) return(list(alpha = 0.10, lwd = 0.30))
  if (n_seed >= 100) return(list(alpha = 0.12, lwd = 0.32))
  if (n_seed >= 50) return(list(alpha = 0.28, lwd = 0.58))
  if (n_seed >= 25) return(list(alpha = 0.38, lwd = 0.72))
  if (n_seed >= 10) return(list(alpha = 0.50, lwd = 0.90))
  if (n_seed >= 5) return(list(alpha = 0.62, lwd = 1.10))
  list(alpha = 0.78, lwd = 1.30)
}

adaptive_seed_overlap_density <- function(z, x_col, y_col, y_range, n_y_bins = 24L) {
  x <- z[[x_col]]
  y <- z[[y_col]]
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok) || !all(is.finite(y_range))) return(NA_real_)
  y_min <- min(y_range, na.rm = TRUE)
  y_max <- max(y_range, na.rm = TRUE)
  if (!is.finite(y_min) || !is.finite(y_max) || y_min == y_max) return(NA_real_)

  breaks <- seq(y_min, y_max, length.out = n_y_bins + 1L)
  x_key <- sprintf("%.8f", as.numeric(x[ok]))
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

adaptive_seed_line_style <- function(z, x_col, y_col, y_range) {
  n_seed <- data.table::uniqueN(z$seed_id)
  style <- adaptive_seed_base_line_style(n_seed)
  density <- adaptive_seed_overlap_density(z, x_col, y_col, y_range)
  if (!is.finite(density) || density <= 0) return(style)

  density_target <- 12
  visibility_boost <- max(1, sqrt(density_target / density))
  style$alpha <- min(0.85, style$alpha * visibility_boost)
  style$lwd <- min(1.60, style$lwd * sqrt(visibility_boost))
  style
}

plot_o2_seed_convergence_panels <- function(comp, figure_dir, o2_values = c(0, 0.1, 0.5, 1, 2, 5),
                                            plot_name = "convergence_to_eigen_by_time_o2_seed_curves_panels",
                                            plot_title = NULL,
                                            adaptive_seed_lines = FALSE,
                                            robust_y_axis = FALSE,
                                            robust_y_quantile = 0.99) {
  available_o2 <- sort(unique(comp$O2_pct[is.finite(comp$O2_pct)]))
  if (!length(available_o2) || !nrow(comp)) return(character())
  matched_o2 <- vapply(o2_values, function(o2) {
    available_o2[which.min(abs(available_o2 - o2))]
  }, numeric(1L))
  if (any(abs(matched_o2 - o2_values) > 1e-9)) {
    warning("Some requested O2 values were matched to nearest available grid value: ",
            paste(o2_values, "->", matched_o2, collapse = ", "))
  }

  d0 <- comp[O2_pct %in% matched_o2]
  d0[, o2_panel := factor(format_num(O2_pct, 6L), levels = format_num(matched_o2, 6L))]

  panels <- list(
    list(
      title = "|2N - eigen|",
      value_col = "abs_delta_2N_minus_eigen",
      color = "#0072B2"
    ),
    list(
      title = "|4N - eigen|",
      value_col = "abs_delta_4N_minus_eigen",
      color = "#D55E00"
    ),
    list(
      title = "|4N - 2N|",
      value_col = "abs_delta_4N_minus_2N",
      color = "#D7191C"
    )
  )

  summary_o2 <- data.table::rbindlist(lapply(panels, function(panel) {
    z <- d0[, .(
      median = stats::median(get(panel$value_col), na.rm = TRUE),
      mean = mean(get(panel$value_col), na.rm = TRUE)
    ), by = .(O2_pct, o2_panel, day)]
    z[, `:=`(metric = panel$title, color = panel$color, value_col = panel$value_col)]
    z
  }), use.names = TRUE)

  y_all <- unlist(lapply(panels, function(panel) d0[[panel$value_col]]), use.names = FALSE)
  finite_y <- y_all[is.finite(y_all)]
  yr <- safe_range(finite_y, c(0, 1))
  if (robust_y_axis && length(finite_y)) {
    summary_y <- c(summary_o2$median, summary_o2$mean, 0.01, 0.05)
    summary_y <- summary_y[is.finite(summary_y)]
    upper <- max(
      as.numeric(stats::quantile(finite_y, robust_y_quantile, na.rm = TRUE, names = FALSE)),
      summary_y,
      na.rm = TRUE
    )
    if (is.finite(upper) && upper > 0) yr <- c(0, upper * 1.06)
  }
  x_range <- range(d0$day, na.rm = TRUE)
  o2_labels <- paste0("O2 = ", format(matched_o2, scientific = FALSE, trim = TRUE), "%")

  save_plot_pair(plot_name, figure_dir, {
    op <- par(mfrow = c(3, length(matched_o2)),
              mar = c(1.0, 1.0, 1.2, 0.4),
              oma = c(3.4, 6.4, if (is.null(plot_title)) 1.2 else 2.2, 0.2),
              las = 1)
    on.exit(par(op), add = TRUE)

    for (row_i in seq_along(panels)) {
      panel <- panels[[row_i]]
      for (col_i in seq_along(matched_o2)) {
        o2 <- matched_o2[[col_i]]
        z <- d0[O2_pct == o2][order(seed_number, day)]
        plot(x_range, yr, type = "n",
             xlab = "", ylab = "",
             xaxt = if (row_i == length(panels)) "s" else "n",
             yaxt = if (col_i == 1L) "s" else "n")
        if (row_i == 1L) mtext(o2_labels[[col_i]], side = 3, line = 0.2, cex = 0.85)
        if (col_i == 1L) {
          mtext(panel$title, side = 2, line = 4.6, las = 3, font = 2, col = panel$color, cex = 0.85)
        }
        abline(h = c(0.01, 0.05), col = "grey84", lty = 3)
        seed_style <- if (adaptive_seed_lines) {
          adaptive_seed_line_style(z, "day", panel$value_col, yr)
        } else {
          list(alpha = 0.24, lwd = 0.42)
        }
        for (seed in unique(z$seed_id[order(z$seed_number)])) {
          zz <- z[seed_id == seed][order(day)]
          gap_col <- if (isTRUE(any(zz$eigen_spectral_gap < 0.01, na.rm = TRUE))) {
            grDevices::adjustcolor("#CC79A7", alpha.f = seed_style$alpha)
          } else {
            grDevices::adjustcolor("#009E73", alpha.f = seed_style$alpha)
          }
          lines(zz$day, zz[[panel$value_col]],
                col = gap_col, lwd = seed_style$lwd)
        }
        s <- summary_o2[O2_pct == o2 & value_col == panel$value_col][order(day)]
        lines(s$day, s$mean, col = panel$color, lwd = 1.9, lty = 2)
        lines(s$day, s$median, col = panel$color, lwd = 3.0)
        if (row_i == 1L && col_i == length(matched_o2)) {
          legend("topright",
                 legend = c("seed value, gap < 0.01", "seed value, gap >= 0.01", "overall median", "overall mean"),
                 col = c(grDevices::adjustcolor("#CC79A7", alpha.f = 0.65),
                         grDevices::adjustcolor("#009E73", alpha.f = 0.65),
                         panel$color, panel$color),
                 lwd = c(0.9, 0.7, 3.0, 1.9), lty = c(1, 1, 1, 2), bty = "n", cex = 0.70)
        }
      }
    }
    mtext("Prediction day", side = 1, outer = TRUE, line = 2.0)
    mtext("Absolute difference", side = 2, outer = TRUE, line = 2.0, las = 3)
    if (!is.null(plot_title)) mtext(plot_title, side = 3, outer = TRUE, line = 0.4, font = 2)
  }, width = 14.5, height = 9.5)
}

plot_o2_seed_convergence_panels_by_class <- function(comp, figure_dir, class_col = "eigen_curve_class",
                                                     o2_values = c(0, 0.1, 0.5, 1, 2, 5)) {
  out_dir <- file.path(figure_dir, "convergence_to_eigen_by_time_o2_seed_curves_by_class")
  classes <- sort(unique(comp[[class_col]][!is.na(comp[[class_col]])]))
  if (!length(classes)) return(character())
  out <- list()
  for (class_label in classes) {
    z <- comp[comp[[class_col]] == class_label]
    n_seed <- data.table::uniqueN(z$seed_id)
    nm <- paste0("convergence_to_eigen_by_time_o2_seed_curves_", slugify(class_label))
    out[[class_label]] <- plot_o2_seed_convergence_panels(
      z,
      out_dir,
      o2_values = o2_values,
      plot_name = nm,
      plot_title = paste0(class_label, " (n = ", n_seed, " seeds)"),
      adaptive_seed_lines = TRUE,
      robust_y_axis = TRUE,
      robust_y_quantile = 0.99
    )
  }
  unlist(out, use.names = FALSE)
}

plot_fraction_within <- function(summary_by_time, figure_dir) {
  d <- summary_by_time[order(day)]
  save_plot_pair("fraction_close_to_eigen_by_time", figure_dir, {
    op <- par(mar = c(4.5, 4.7, 1.2, 1), las = 1)
    on.exit(par(op), add = TRUE)
    plot(d$day, d$abs_delta_2N_minus_eigen_fraction_le_0p05, type = "n",
         ylim = c(0, 1), xlab = "Prediction day", ylab = "Fraction of seed-O2 pairs")
    lines(d$day, d$abs_delta_2N_minus_eigen_fraction_le_0p05, col = "#0072B2", lwd = 2)
    lines(d$day, d$abs_delta_4N_minus_eigen_fraction_le_0p05, col = "#D55E00", lwd = 2)
    lines(d$day, d$abs_delta_4N_minus_2N_fraction_le_0p05, col = "#009E73", lwd = 2)
    lines(d$day, d$abs_delta_2N_minus_eigen_fraction_le_0p01, col = "#0072B2", lwd = 1.5, lty = 2)
    lines(d$day, d$abs_delta_4N_minus_eigen_fraction_le_0p01, col = "#D55E00", lwd = 1.5, lty = 2)
    lines(d$day, d$abs_delta_4N_minus_2N_fraction_le_0p01, col = "#009E73", lwd = 1.5, lty = 2)
    legend("bottomright",
           legend = c("|2N - eigen| <=0.05", "|4N - eigen| <=0.05", "|4N - 2N| <=0.05", "<=0.01 shown dashed"),
           col = c("#0072B2", "#D55E00", "#009E73", "grey30"),
           lty = c(1, 1, 1, 2), lwd = c(2, 2, 2, 1.5), bty = "n", cex = 0.85)
  }, width = 8.5, height = 5.2)
}

plot_class_match <- function(class_summary, figure_dir) {
  d <- class_summary[order(day)]
  save_plot_pair("curve_class_match_fraction_by_time", figure_dir, {
    op <- par(mar = c(4.5, 4.7, 1.2, 1), las = 1)
    on.exit(par(op), add = TRUE)
    plot(d$day, d$fraction_init_2N_matches_eigen, type = "n",
         ylim = c(0, 1), xlab = "Prediction day", ylab = "Fraction of seeds")
    lines(d$day, d$fraction_init_2N_matches_eigen, col = "#0072B2", lwd = 2)
    lines(d$day, d$fraction_init_4N_matches_eigen, col = "#D55E00", lwd = 2)
    lines(d$day, d$fraction_init_2N_matches_init_4N, col = "#009E73", lwd = 2)
    legend("bottomright",
           legend = c("2N class = eigen class", "4N class = eigen class", "2N class = 4N class"),
           col = c("#0072B2", "#D55E00", "#009E73"), lty = 1, lwd = 2, bty = "n", cex = 0.85)
  }, width = 8.5, height = 5.2)
}

plot_terminal_class_consistency_summary <- function(class_summary, figure_dir, terminal_day) {
  d <- class_summary[order(match_group)]
  cols <- c(
    "init_2N_vs_eigen" = "#0072B2",
    "init_4N_vs_eigen" = "#D55E00",
    "init_2N_vs_init_4N" = "#009E73",
    "all_three_same" = "#7B3294"
  )
  d[, plot_label := c("2N = eigen", "4N = eigen", "2N = 4N", "all three")[match(match_group, names(cols))]]
  save_plot_pair(paste0("day", terminal_day, "_curve_class_consistency_summary"), figure_dir, {
    op <- par(mar = c(5.8, 4.6, 1.4, 0.8), las = 1)
    on.exit(par(op), add = TRUE)
    x <- barplot(d$fraction_match,
                 names.arg = d$plot_label,
                 col = cols[d$match_group],
                 border = NA,
                 ylim = c(0, 1),
                 ylab = "Fraction of seeds",
                 las = 2)
    text(x, pmin(d$fraction_match + 0.04, 0.98),
         labels = paste0(d$n_match, "/", d$n_seed),
         cex = 0.8)
    abline(h = seq(0, 1, 0.25), col = "grey88", lty = 3)
  }, width = 6.6, height = 5.2)
}

plot_terminal_class_pair_heatmaps <- function(pair_counts, figure_dir, terminal_day) {
  comparisons <- c("init_2N_vs_eigen", "init_4N_vs_eigen", "init_2N_vs_init_4N")
  d <- pair_counts[comparison %in% comparisons]
  if (!nrow(d)) return(character())
  all_labels <- sort(unique(c(d$label_a, d$label_b)))
  max_count <- max(d$n_seed, na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) max_count <- 1
  save_plot_pair(paste0("day", terminal_day, "_curve_class_pair_heatmaps"), figure_dir, {
    op <- par(mfrow = c(1, length(comparisons)), mar = c(8.8, 7.8, 2.0, 1.2), oma = c(0, 0, 0, 3.2), las = 2)
    on.exit(par(op), add = TRUE)
    cols <- grDevices::hcl.colors(80, "YlGnBu", rev = FALSE)
    for (comparison_label in comparisons) {
      z <- d[d$comparison == comparison_label]
      mat <- matrix(0, nrow = length(all_labels), ncol = length(all_labels), dimnames = list(all_labels, all_labels))
      if (nrow(z)) {
        mat[cbind(match(z$label_a, all_labels), match(z$label_b, all_labels))] <- z$n_seed
      }
      image(seq_along(all_labels), seq_along(all_labels), mat,
            col = cols,
            breaks = seq(0, max_count, length.out = length(cols) + 1L),
            xlab = "", ylab = "", xaxt = "n", yaxt = "n",
            main = switch(comparison_label,
                          init_2N_vs_eigen = "2N vs eigen",
                          init_4N_vs_eigen = "4N vs eigen",
                          init_2N_vs_init_4N = "2N vs 4N"))
      axis(1, at = seq_along(all_labels), labels = all_labels, cex.axis = 0.65)
      axis(2, at = seq_along(all_labels), labels = all_labels, cex.axis = 0.65)
      mtext("Label A", side = 2, line = 6.2, las = 0, cex = 0.85)
      mtext("Label B", side = 1, line = 7.0, cex = 0.85)
      box()
      for (i in seq_len(nrow(mat))) {
        for (j in seq_len(ncol(mat))) {
          if (mat[i, j] > 0) text(i, j, labels = mat[i, j], cex = 0.58)
        }
      }
    }
    mtext("Seed count", side = 4, outer = TRUE, line = 1.0, las = 3)
  }, width = 15, height = 6.6)
}

plot_gap_decile_heatmap <- function(summary_gap_decile, figure_dir, metric, name, zlab) {
  d <- summary_gap_decile[is.finite(gap_decile_number)]
  days <- sort(unique(d$day))
  deciles <- sort(unique(d$gap_decile_number))
  mat <- matrix(NA_real_, nrow = length(days), ncol = length(deciles), dimnames = list(days, deciles))
  for (i in seq_len(nrow(d))) {
    mat[as.character(d$day[[i]]), as.character(d$gap_decile_number[[i]])] <- d[[metric]][[i]]
  }
  save_plot_pair(name, figure_dir, {
    op <- par(mar = c(4.5, 4.7, 1.5, 4.5), las = 1)
    on.exit(par(op), add = TRUE)
    image(x = deciles, y = days, z = t(mat), col = grDevices::hcl.colors(80, "YlOrRd", rev = FALSE),
          xlab = "Eigen spectral-gap decile (1 = smallest gap)", ylab = "Prediction day")
    axis(1, at = deciles)
    axis(2, at = days)
    box()
    mtext(zlab, side = 4, line = 3)
  }, width = 8.5, height = 5.8)
}

plot_terminal_day_heatmaps <- function(comp, figure_dir, terminal_day) {
  d <- comp[day == terminal_day]
  if (!nrow(d)) return(character())
  order_seed <- unique(d[order(eigen_curve_class, eigen_seed_min_spectral_gap, seed_number)]$seed_id)
  o2 <- sort(unique(d$O2_pct))
  make_mat <- function(col) {
    mat <- matrix(NA_real_, nrow = length(order_seed), ncol = length(o2),
                  dimnames = list(order_seed, format_num(o2, 6L)))
    idx_seed <- match(d$seed_id, order_seed)
    idx_o2 <- match(format_num(d$O2_pct, 6L), colnames(mat))
    mat[cbind(idx_seed, idx_o2)] <- d[[col]]
    mat
  }
  m2 <- make_mat("delta_2N_minus_eigen")
  m4 <- make_mat("delta_4N_minus_eigen")
  mi <- make_mat("delta_4N_minus_2N")
  lim <- max(abs(c(m2, m4, mi)), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  breaks <- seq(-lim, lim, length.out = 101L)
  cols <- grDevices::colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100L)
  save_plot_pair(paste0("day", terminal_day, "_signed_delta_heatmaps"), figure_dir, {
    op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 3))
    on.exit(par(op), add = TRUE)
    panels <- list("2N - eigen" = m2, "4N - eigen" = m4, "4N - 2N" = mi)
    for (nm in names(panels)) {
      mat <- panels[[nm]]
      image(x = o2, y = seq_len(nrow(mat)), z = t(mat[nrow(mat):1L, , drop = FALSE]),
            breaks = breaks, col = cols, xlab = "O2 (%)", ylab = "Seeds ordered by eigen class/gap",
            main = nm)
      box()
    }
    mtext("Signed mean-ploidy difference", side = 4, outer = TRUE, line = 1)
  }, width = 13, height = 5.2)
}

plot_terminal_day_delta_vs_gap <- function(comp, figure_dir, terminal_day) {
  d <- comp[day == terminal_day & is.finite(eigen_spectral_gap)]
  if (!nrow(d)) return(character())
  if (nrow(d) > 60000L) {
    set.seed(1L)
    d <- d[sample.int(nrow(d), 60000L)]
  }
  save_plot_pair(paste0("day", terminal_day, "_abs_delta_vs_eigen_spectral_gap"), figure_dir, {
    op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2, 1))
    on.exit(par(op), add = TRUE)
    panels <- list(
      "|2N - eigen|" = "abs_delta_2N_minus_eigen",
      "|4N - eigen|" = "abs_delta_4N_minus_eigen",
      "|4N - 2N|" = "abs_delta_4N_minus_2N"
    )
    for (nm in names(panels)) {
      y <- pmax(d[[panels[[nm]]]], .Machine$double.eps)
      x <- pmax(d$eigen_spectral_gap, .Machine$double.eps)
      plot(x, y, log = "xy", pch = 16, cex = 0.25,
           col = grDevices::adjustcolor("#333333", alpha.f = 0.18),
           xlab = "Eigen spectral gap", ylab = "Absolute difference", main = nm)
      abline(v = c(0.005, 0.01), lty = 2, col = "grey50")
      abline(h = c(0.01, 0.05), lty = 3, col = "grey60")
    }
  }, width = 13, height = 4.8)
}

plot_representative_curves <- function(comp, day1000_seed, figure_dir, plot_days) {
  candidates <- unique(c(
    head(day1000_seed[order(-max_abs_delta_2N_minus_eigen)]$seed_id, 2L),
    head(day1000_seed[order(max_abs_delta_2N_minus_eigen)]$seed_id, 1L),
    head(day1000_seed[order(eigen_seed_min_spectral_gap)]$seed_id, 1L)
  ))
  candidates <- candidates[!is.na(candidates)]
  if (!length(candidates)) return(character())
  d <- comp[seed_id %in% candidates & day %in% plot_days]
  save_plot_pair("representative_seed_o2_curves", figure_dir, {
    op <- par(mfrow = c(length(candidates), 2), mar = c(3.5, 4, 2, 1), oma = c(1, 0, 0, 0))
    on.exit(par(op), add = TRUE)
    day_cols <- grDevices::hcl.colors(length(plot_days), "Viridis")
    names(day_cols) <- as.character(sort(plot_days))
    for (seed in candidates) {
      z <- d[seed_id == seed]
      z <- z[order(O2_pct, day)]
      seed_meta <- day1000_seed[seed_id == seed][1L]
      y_range <- safe_range(c(z$eigen_dominant_mean_ploidy, z$mean_ploidy_init_2N, z$mean_ploidy_init_4N), c(1.5, 2.5))
      plot(range(z$O2_pct), y_range, type = "n", xlab = "O2 (%)", ylab = "Mean ploidy",
           main = paste(seed, "2N start", seed_meta$eigen_curve_class))
      eigen <- unique(z[, .(O2_pct, eigen_dominant_mean_ploidy)])
      eigen <- eigen[order(O2_pct)]
      lines(eigen$O2_pct, eigen$eigen_dominant_mean_ploidy, col = "black", lwd = 2.5)
      for (day_i in sort(plot_days)) {
        zz <- z[day == day_i][order(O2_pct)]
        lines(zz$O2_pct, zz$mean_ploidy_init_2N, col = grDevices::adjustcolor(day_cols[[as.character(day_i)]], alpha.f = 0.8), lwd = 1)
      }
      abline(h = c(1.5, 2.0), col = "grey70", lty = 3)
      plot(range(z$O2_pct), y_range, type = "n", xlab = "O2 (%)", ylab = "Mean ploidy",
           main = paste(seed, "4N start", "min gap", format_num(seed_meta$eigen_seed_min_spectral_gap, 3L)))
      lines(eigen$O2_pct, eigen$eigen_dominant_mean_ploidy, col = "black", lwd = 2.5)
      for (day_i in sort(plot_days)) {
        zz <- z[day == day_i][order(O2_pct)]
        lines(zz$O2_pct, zz$mean_ploidy_init_4N, col = grDevices::adjustcolor(day_cols[[as.character(day_i)]], alpha.f = 0.8), lwd = 1)
      }
      abline(h = c(1.5, 2.0), col = "grey70", lty = 3)
    }
    legend("bottom", inset = -0.02, xpd = NA, horiz = TRUE, bty = "n", cex = 0.75,
           legend = c("eigen", paste0("day ", sort(plot_days))),
           col = c("black", day_cols[as.character(sort(plot_days))]),
           lwd = c(2.5, rep(1, length(plot_days))))
  }, width = 12, height = max(5, 3.1 * length(candidates)))
}

build_comparison <- function(paths, plot_days = c(100, 200, 300, 500, 700, 1000), analysis_max_day = 1000L) {
  require_data_table()
  analysis_max_day <- as.integer(analysis_max_day)
  if (!is.finite(analysis_max_day) || analysis_max_day < 0L) {
    stop("analysis_max_day must be a non-negative integer.")
  }
  dir.create(paths$out_tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$out_figures, recursive = TRUE, showWarnings = FALSE)

  eigen_curves <- read_dt(paths$eigen_curves, select = c(
    "seed_id", "O2_pct", "O2_key", "dominant_mean_ploidy", "trajectory_regime",
    "mode_label", "status", "dominant_growth_rate", "spectral_gap", "objective"
  ))
  data.table::setnames(
    eigen_curves,
    c("dominant_mean_ploidy", "trajectory_regime", "mode_label", "status", "dominant_growth_rate", "spectral_gap", "objective"),
    c("eigen_dominant_mean_ploidy", "eigen_trajectory_regime", "eigen_mode_label", "eigen_status", "eigen_dominant_growth_rate", "eigen_spectral_gap", "eigen_objective")
  )

  eigen_by_seed <- read_dt(paths$eigen_by_seed, select = c(
    "seed_id", "seed_number", "curve_class", "final_interpretation_class", "sign_sequence",
    "ploidy_range", "net_ploidy_change", "min_spectral_gap", "median_spectral_gap",
    "fraction_o2_gap_below_0p005", "fraction_o2_gap_below_0p01",
    "monotonicity_reliability", "objective", "objective_source"
  ))
  data.table::setnames(
    eigen_by_seed,
    c("curve_class", "final_interpretation_class", "sign_sequence", "ploidy_range",
      "net_ploidy_change", "min_spectral_gap", "median_spectral_gap",
      "fraction_o2_gap_below_0p005", "fraction_o2_gap_below_0p01",
      "monotonicity_reliability", "objective", "objective_source"),
    c("eigen_curve_class", "eigen_final_interpretation_class", "eigen_sign_sequence", "eigen_ploidy_range",
      "eigen_net_ploidy_change", "eigen_seed_min_spectral_gap", "eigen_seed_median_spectral_gap",
      "eigen_seed_fraction_o2_gap_below_0p005", "eigen_seed_fraction_o2_gap_below_0p01",
      "eigen_monotonicity_reliability", "eigen_seed_objective", "eigen_objective_source")
  )
  eigen <- merge(eigen_curves, eigen_by_seed, by = "seed_id", all.x = TRUE, sort = FALSE)
  data.table::setkey(eigen, seed_id, O2_key)

  initial <- read_dt(paths$initial_selected, select = c(
    "seed_id", "seed_number", "day", "initial_condition", "initial_ploidy",
    "O2_pct", "O2_key", "used_initial_N", "status", "trajectory_method",
    "mean_ploidy", "dominant_mean_ploidy", "spectral_gap", "relax_time_days",
    "time_to_10x_days", "time_to_100x_days", "dominance_class"
  ))
  initial <- initial[day <= analysis_max_day]
  if (!nrow(initial)) stop("No initial-ploidy rows found with day <= analysis_max_day: ", analysis_max_day)
  available_days <- sort(unique(initial$day))
  terminal_day <- max(available_days, na.rm = TRUE)
  plot_days <- sort(unique(plot_days[plot_days <= terminal_day]))
  if (!length(plot_days)) plot_days <- terminal_day

  init_2n <- initial[initial_condition == "init_2N"]
  init_4n <- initial[initial_condition == "init_4N"]
  init_2n <- prefix_cols(init_2n, "init_2N_", except = c("seed_id", "seed_number", "day", "O2_pct", "O2_key"))
  init_4n <- prefix_cols(init_4n, "init_4N_", except = c("seed_id", "seed_number", "day", "O2_pct", "O2_key"))
  init_wide <- merge(
    init_2n,
    init_4n,
    by = c("seed_id", "seed_number", "day", "O2_pct", "O2_key"),
    all = TRUE,
    sort = FALSE
  )
  rm(initial, init_2n, init_4n)
  gc()

  class_table <- read_dt(paths$initial_class, select = c(
    "seed_id", "seed_number", "day", "initial_condition", "initial_ploidy", "curve_class",
    "sign_sequence", "ploidy_range", "net_ploidy_change", "min_spectral_gap",
    "median_spectral_gap", "fraction_o2_gap_below_0p005", "fraction_o2_gap_below_0p01",
    "curve_class_2N", "curve_class_4N", "curve_class_consistent",
    "sign_sequence_2N", "sign_sequence_4N", "sign_sequence_consistent"
  ))
  class_table <- class_table[day <= analysis_max_day]
  class_2n <- class_table[initial_condition == "init_2N"]
  class_4n <- class_table[initial_condition == "init_4N"]
  class_2n <- prefix_cols(class_2n, "init_2N_", except = c("seed_id", "seed_number", "day"))
  class_4n <- prefix_cols(class_4n, "init_4N_", except = c("seed_id", "seed_number", "day"))
  class_wide <- merge(class_2n, class_4n, by = c("seed_id", "seed_number", "day"), all = TRUE, sort = FALSE)
  class_wide[, `:=`(
    curve_class_init_2N = init_2N_curve_class,
    curve_class_init_4N = init_4N_curve_class,
    sign_sequence_init_2N = init_2N_sign_sequence,
    sign_sequence_init_4N = init_4N_sign_sequence,
    curve_class_2N_vs_4N_consistent = init_2N_curve_class == init_4N_curve_class
  )]
  class_keep <- c(
    "seed_id", "seed_number", "day", "curve_class_init_2N", "curve_class_init_4N",
    "sign_sequence_init_2N", "sign_sequence_init_4N", "curve_class_2N_vs_4N_consistent",
    "init_2N_ploidy_range", "init_4N_ploidy_range", "init_2N_net_ploidy_change",
    "init_4N_net_ploidy_change", "init_2N_min_spectral_gap", "init_4N_min_spectral_gap",
    "init_2N_median_spectral_gap", "init_4N_median_spectral_gap",
    "init_2N_fraction_o2_gap_below_0p005", "init_4N_fraction_o2_gap_below_0p005",
    "init_2N_fraction_o2_gap_below_0p01", "init_4N_fraction_o2_gap_below_0p01"
  )
  class_wide <- class_wide[, ..class_keep]
  rm(class_table, class_2n, class_4n)
  gc()

  comp <- merge(init_wide, eigen, by = c("seed_id", "O2_key"), all.x = TRUE, sort = FALSE)
  if (all(c("O2_pct.x", "O2_pct.y") %in% names(comp))) {
    comp[, O2_pct := O2_pct.x]
    comp[, c("O2_pct.x", "O2_pct.y") := NULL]
  }
  if (all(c("seed_number.x", "seed_number.y") %in% names(comp))) {
    comp[, seed_number := data.table::fifelse(!is.na(seed_number.x), seed_number.x, seed_number.y)]
    comp[, c("seed_number.x", "seed_number.y") := NULL]
  }
  comp <- merge(comp, class_wide, by = c("seed_id", "seed_number", "day"), all.x = TRUE, sort = FALSE)

  comp[, `:=`(
    delta_2N_minus_eigen = init_2N_mean_ploidy - eigen_dominant_mean_ploidy,
    delta_4N_minus_eigen = init_4N_mean_ploidy - eigen_dominant_mean_ploidy,
    delta_4N_minus_2N = init_4N_mean_ploidy - init_2N_mean_ploidy
  )]
  comp[, `:=`(
    abs_delta_2N_minus_eigen = abs(delta_2N_minus_eigen),
    abs_delta_4N_minus_eigen = abs(delta_4N_minus_eigen),
    abs_delta_4N_minus_2N = abs(delta_4N_minus_2N),
    same_class_2N_vs_eigen = curve_class_init_2N == eigen_curve_class,
    same_class_4N_vs_eigen = curve_class_init_4N == eigen_curve_class,
    same_class_2N_vs_4N = curve_class_init_2N == curve_class_init_4N,
    eigen_gap_region = data.table::fifelse(
      eigen_spectral_gap < 0.005,
      "gap_lt_0p005",
      data.table::fifelse(eigen_spectral_gap < 0.01, "gap_0p005_to_0p01", "gap_ge_0p01")
    ),
    eigen_small_gap_0p005 = eigen_spectral_gap < 0.005,
    eigen_small_gap_0p01 = eigen_spectral_gap < 0.01
  )]
  gap_q <- unique(stats::quantile(comp$eigen_spectral_gap, probs = seq(0, 1, 0.1), na.rm = TRUE, names = FALSE))
  if (length(gap_q) > 2L) {
    comp[, gap_decile_number := as.integer(cut(eigen_spectral_gap, breaks = gap_q, include.lowest = TRUE, labels = FALSE))]
  } else {
    comp[, gap_decile_number := NA_integer_]
  }
  comp[, gap_decile_label := ifelse(is.na(gap_decile_number), NA_character_, paste0("D", gap_decile_number))]

  data.table::setcolorder(comp, c(
    "seed_id", "seed_number", "O2_pct", "O2_key", "day",
    "eigen_dominant_mean_ploidy", "init_2N_mean_ploidy", "init_4N_mean_ploidy",
    "delta_2N_minus_eigen", "delta_4N_minus_eigen", "delta_4N_minus_2N",
    "abs_delta_2N_minus_eigen", "abs_delta_4N_minus_eigen", "abs_delta_4N_minus_2N",
    "eigen_spectral_gap", "eigen_gap_region", "gap_decile_number", "gap_decile_label",
    "eigen_curve_class", "curve_class_init_2N", "curve_class_init_4N",
    "same_class_2N_vs_eigen", "same_class_4N_vs_eigen", "same_class_2N_vs_4N"
  ))

  summary_by_time <- wide_metric_summary(comp, "day")[order(day)]
  summary_by_gap_region <- wide_metric_summary(comp, c("day", "eigen_gap_region"))[order(day, eigen_gap_region)]
  summary_by_gap_decile <- wide_metric_summary(comp, c("day", "gap_decile_number", "gap_decile_label"))[order(day, gap_decile_number)]

  class_comp <- unique(comp[, .(
    seed_id, seed_number, day,
    eigen_curve_class, eigen_final_interpretation_class, eigen_sign_sequence,
    eigen_seed_min_spectral_gap, eigen_seed_median_spectral_gap,
    eigen_monotonicity_reliability,
    curve_class_init_2N, curve_class_init_4N,
    sign_sequence_init_2N, sign_sequence_init_4N,
    same_class_2N_vs_eigen, same_class_4N_vs_eigen, same_class_2N_vs_4N
  )])
  class_summary <- class_comp[, .(
    n_seed = .N,
    fraction_init_2N_matches_eigen = mean(same_class_2N_vs_eigen, na.rm = TRUE),
    fraction_init_4N_matches_eigen = mean(same_class_4N_vs_eigen, na.rm = TRUE),
    fraction_init_2N_matches_init_4N = mean(same_class_2N_vs_4N, na.rm = TRUE),
    n_monotone_increasing_eigen = sum(eigen_curve_class == "monotone_increasing", na.rm = TRUE),
    n_monotone_increasing_init_2N = sum(curve_class_init_2N == "monotone_increasing", na.rm = TRUE),
    n_monotone_increasing_init_4N = sum(curve_class_init_4N == "monotone_increasing", na.rm = TRUE)
  ), by = "day"][order(day)]
  class_summary[, `:=`(
    fraction_monotone_increasing_eigen = n_monotone_increasing_eigen / n_seed,
    fraction_monotone_increasing_init_2N = n_monotone_increasing_init_2N / n_seed,
    fraction_monotone_increasing_init_4N = n_monotone_increasing_init_4N / n_seed
  )]

  terminal_prefix <- paste0("day", terminal_day)
  terminal_seed <- comp[day == terminal_day, .(
    n_o2 = .N,
    eigen_curve_class = eigen_curve_class[1L],
    curve_class_init_2N = curve_class_init_2N[1L],
    curve_class_init_4N = curve_class_init_4N[1L],
    same_class_2N_vs_eigen = same_class_2N_vs_eigen[1L],
    same_class_4N_vs_eigen = same_class_4N_vs_eigen[1L],
    same_class_2N_vs_4N = same_class_2N_vs_4N[1L],
    eigen_seed_min_spectral_gap = eigen_seed_min_spectral_gap[1L],
    eigen_seed_median_spectral_gap = eigen_seed_median_spectral_gap[1L],
    mean_abs_delta_2N_minus_eigen = mean(abs_delta_2N_minus_eigen, na.rm = TRUE),
    mean_abs_delta_4N_minus_eigen = mean(abs_delta_4N_minus_eigen, na.rm = TRUE),
    mean_abs_delta_4N_minus_2N = mean(abs_delta_4N_minus_2N, na.rm = TRUE),
    max_abs_delta_2N_minus_eigen = max(abs_delta_2N_minus_eigen, na.rm = TRUE),
    max_abs_delta_4N_minus_eigen = max(abs_delta_4N_minus_eigen, na.rm = TRUE),
    max_abs_delta_4N_minus_2N = max(abs_delta_4N_minus_2N, na.rm = TRUE),
    fraction_o2_abs_2N_eigen_le_0p05 = mean(abs_delta_2N_minus_eigen <= 0.05, na.rm = TRUE),
    fraction_o2_abs_4N_eigen_le_0p05 = mean(abs_delta_4N_minus_eigen <= 0.05, na.rm = TRUE),
    fraction_o2_abs_4N_2N_le_0p05 = mean(abs_delta_4N_minus_2N <= 0.05, na.rm = TRUE)
  ), by = .(seed_id, seed_number)][order(-max_abs_delta_2N_minus_eigen)]

  top_discrepant <- comp[day == terminal_day][order(-pmax(abs_delta_2N_minus_eigen, abs_delta_4N_minus_eigen, abs_delta_4N_minus_2N))]
  top_discrepant <- top_discrepant[seq_len(min(200L, nrow(top_discrepant)))]

  terminal_class_by_seed <- terminal_seed[, .(
    seed_id,
    seed_number,
    terminal_day = terminal_day,
    eigen_curve_class,
    curve_class_init_2N,
    curve_class_init_4N,
    same_class_2N_vs_eigen,
    same_class_4N_vs_eigen,
    same_class_2N_vs_4N,
    all_three_curve_class_consistent = same_class_2N_vs_eigen & same_class_4N_vs_eigen,
    eigen_seed_min_spectral_gap,
    eigen_seed_median_spectral_gap,
    mean_abs_delta_2N_minus_eigen,
    mean_abs_delta_4N_minus_eigen,
    mean_abs_delta_4N_minus_2N,
    max_abs_delta_2N_minus_eigen,
    max_abs_delta_4N_minus_eigen,
    max_abs_delta_4N_minus_2N
  )][order(eigen_curve_class, seed_number)]

  class_n_seed <- nrow(terminal_class_by_seed)
  terminal_class_summary <- data.table::data.table(
    match_group = c("init_2N_vs_eigen", "init_4N_vs_eigen", "init_2N_vs_init_4N", "all_three_same"),
    terminal_day = terminal_day,
    n_seed = class_n_seed,
    n_match = c(
      sum(terminal_class_by_seed$same_class_2N_vs_eigen, na.rm = TRUE),
      sum(terminal_class_by_seed$same_class_4N_vs_eigen, na.rm = TRUE),
      sum(terminal_class_by_seed$same_class_2N_vs_4N, na.rm = TRUE),
      sum(terminal_class_by_seed$all_three_curve_class_consistent, na.rm = TRUE)
    ),
    definition = c(
      "curve_class_init_2N equals eigen_curve_class at terminal day",
      "curve_class_init_4N equals eigen_curve_class at terminal day",
      "curve_class_init_2N equals curve_class_init_4N at terminal day",
      "eigen_curve_class, curve_class_init_2N, and curve_class_init_4N are all equal at terminal day"
    )
  )
  terminal_class_summary[, fraction_match := n_match / n_seed]

  pair_one <- function(x, label_a_col, label_b_col, comparison_label) {
    z <- x[, .N, by = c(label_a_col, label_b_col)]
    data.table::setnames(z, c(label_a_col, label_b_col, "N"), c("label_a", "label_b", "n_seed"))
    z[, `:=`(
      comparison = comparison_label,
      terminal_day = terminal_day,
      fraction_of_seed = n_seed / class_n_seed,
      same_label = label_a == label_b
    )]
    z[, .(comparison, terminal_day, label_a, label_b, same_label, n_seed, fraction_of_seed)]
  }
  terminal_class_pair_counts <- data.table::rbindlist(list(
    pair_one(terminal_class_by_seed, "eigen_curve_class", "curve_class_init_2N", "init_2N_vs_eigen"),
    pair_one(terminal_class_by_seed, "eigen_curve_class", "curve_class_init_4N", "init_4N_vs_eigen"),
    pair_one(terminal_class_by_seed, "curve_class_init_2N", "curve_class_init_4N", "init_2N_vs_init_4N")
  ), use.names = TRUE)
  terminal_class_pair_counts <- terminal_class_pair_counts[order(comparison, label_a, label_b)]

  terminal_class_combination_counts <- terminal_class_by_seed[, .(
    n_seed = .N,
    fraction_of_seed = .N / class_n_seed,
    all_three_curve_class_consistent = all_three_curve_class_consistent[1L]
  ), by = .(terminal_day, eigen_curve_class, curve_class_init_2N, curve_class_init_4N)]
  terminal_class_combination_counts <- terminal_class_combination_counts[
    order(-all_three_curve_class_consistent, -n_seed, eigen_curve_class, curve_class_init_2N, curve_class_init_4N)
  ]

  terminal_global <- summary_by_time[day == terminal_day]
  terminal_large_gap <- summary_by_gap_region[day == terminal_day & eigen_gap_region == "gap_ge_0p01"]
  terminal_small_gap <- summary_by_gap_region[day == terminal_day & eigen_gap_region == "gap_lt_0p005"]
  class_terminal <- class_summary[day == terminal_day]
  terminal_label <- paste0("day ", terminal_day)
  evidence_summary <- data.table::data.table(
    statement = c(
      paste0("Median finite-time 2N-start ploidy reaches the eigen steady-state by ", terminal_label, " across seed-O2 pairs."),
      paste0("Median finite-time 4N-start ploidy reaches the eigen steady-state by ", terminal_label, " across seed-O2 pairs."),
      paste0("2N-start and 4N-start finite-time trajectories largely agree by ", terminal_label, " across seed-O2 pairs."),
      paste0("Large-gap seed-O2 pairs are essentially converged to the eigen steady-state by ", terminal_label, "."),
      paste0("Small-gap seed-O2 pairs retain most of the ", terminal_label, " finite-time/eigen discrepancy."),
      "Finite-time/eigen discrepancy is negatively associated with eigen spectral gap.",
      paste0("Curve-class labels for 2N and 4N finite-time curves are mostly consistent by ", terminal_label, "."),
      "Curve-class labels do not fully match the eigen classification even when ploidy values mostly converge, because finite-time curve shape can still differ in small-gap/slow-relaxation regions."
    ),
    evidence_metric = c(
      paste0(terminal_prefix, " median |2N - eigen|"),
      paste0(terminal_prefix, " median |4N - eigen|"),
      paste0(terminal_prefix, " median |4N - 2N|"),
      paste0(terminal_prefix, " large-gap median |2N - eigen|; median |4N - eigen|"),
      paste0(terminal_prefix, " small-gap median |2N - eigen|; median |4N - eigen|"),
      paste0(terminal_prefix, " Spearman(|delta|, eigen spectral gap)"),
      paste0(terminal_prefix, " fraction 2N class = 4N class"),
      paste0(terminal_prefix, " fraction 2N class = eigen class; fraction 4N class = eigen class")
    ),
    value = c(
      format_num(terminal_global$abs_delta_2N_minus_eigen_median, 8L),
      format_num(terminal_global$abs_delta_4N_minus_eigen_median, 8L),
      format_num(terminal_global$abs_delta_4N_minus_2N_median, 8L),
      paste(format_num(terminal_large_gap$abs_delta_2N_minus_eigen_median, 8L),
            format_num(terminal_large_gap$abs_delta_4N_minus_eigen_median, 8L), sep = "; "),
      paste(format_num(terminal_small_gap$abs_delta_2N_minus_eigen_median, 8L),
            format_num(terminal_small_gap$abs_delta_4N_minus_eigen_median, 8L), sep = "; "),
      paste(format_num(terminal_global$spearman_abs_2N_eigen_vs_gap, 4L),
            format_num(terminal_global$spearman_abs_4N_eigen_vs_gap, 4L), sep = "; "),
      format_num(class_terminal$fraction_init_2N_matches_init_4N, 4L),
      paste(format_num(class_terminal$fraction_init_2N_matches_eigen, 4L),
            format_num(class_terminal$fraction_init_4N_matches_eigen, 4L), sep = "; ")
    ),
    source_table = c(
      rep("summary_by_time.tsv", 3L),
      "summary_by_time_and_gap_region.tsv",
      "summary_by_time_and_gap_region.tsv",
      "summary_by_time.tsv",
      "curve_class_match_summary.tsv",
      "curve_class_match_summary.tsv"
    )
  )

  run_args <- data.table::data.table(
    argument = c(
      "script", "result_root", "monotonicity_dir", "initial_ploidy_dir", "compare_parent",
      "analysis_max_day_requested", "terminal_day", "analysis_label", "output_dir",
      "n_seed", "n_o2", "selected_days", "plot_days", "n_comparison_rows",
      "classification_label_eigen", "classification_label_initial_2N", "classification_label_initial_4N",
      "steady_state_interpretation"
    ),
    value = c(
      script_path(), paths$root, paths$mono_dir, paths$init_dir, paths$compare_parent,
      as.character(analysis_max_day), as.character(terminal_day), paths$analysis_label, paths$out_dir,
      as.character(data.table::uniqueN(comp$seed_id)), as.character(data.table::uniqueN(comp$O2_key)),
      paste(sort(unique(comp$day)), collapse = ","), paste(sort(plot_days), collapse = ","),
      as.character(nrow(comp)),
      "eigen_curve_class from dense-grid_monotonicity_classification",
      "curve_class_init_2N from dense-grid_initial-ploidy_trajectory",
      "curve_class_init_4N from dense-grid_initial-ploidy_trajectory",
      "Eigenvector result is treated as the fixed-O2 dominant-attractor steady-state ploidy; finite-time trajectories are compared against it by day and initial ploidy."
    )
  )

  expected_rows <- data.table::uniqueN(comp$seed_id) * data.table::uniqueN(comp$O2_key) * data.table::uniqueN(comp$day)
  validation <- data.table::data.table(
    test_case = c(
      "comparison_rows",
      "seed_count",
      "o2_count",
      "day_count",
      "init_2N_present",
      "init_4N_present",
      "no_missing_eigen_ploidy",
      "no_missing_2N_ploidy",
      "no_missing_4N_ploidy",
      "eigen_class_one_per_seed",
      "initial_class_one_per_seed_day_condition",
      "terminal_day_present",
      "spectral_gap_present",
      "validation_tables_passed"
    ),
    expected = c(
      as.character(expected_rows),
      "500",
      "201",
      as.character(data.table::uniqueN(comp$day)),
      "TRUE", "TRUE", "0", "0", "0", "TRUE", "TRUE", as.character(terminal_day), "TRUE", "TRUE"
    ),
    observed = c(
      as.character(nrow(comp)),
      as.character(data.table::uniqueN(comp$seed_id)),
      as.character(data.table::uniqueN(comp$O2_key)),
      as.character(data.table::uniqueN(comp$day)),
      as.character(any(!is.na(comp$init_2N_mean_ploidy))),
      as.character(any(!is.na(comp$init_4N_mean_ploidy))),
      as.character(sum(is.na(comp$eigen_dominant_mean_ploidy))),
      as.character(sum(is.na(comp$init_2N_mean_ploidy))),
      as.character(sum(is.na(comp$init_4N_mean_ploidy))),
      as.character(all(eigen_by_seed[, .N, by = seed_id]$N == 1L)),
      as.character(all(class_wide[, .N, by = .(seed_id, day)]$N == 1L)),
      as.character(max(comp$day, na.rm = TRUE)),
      as.character(any(is.finite(comp$eigen_spectral_gap))),
      as.character({
        ev <- read_dt(paths$eigen_validation)
        iv <- read_dt(paths$initial_validation)
        all(ev$passed %in% TRUE | ev$passed == "TRUE") && all(iv$passed %in% TRUE | iv$passed == "TRUE")
      })
    )
  )
  validation <- data.table::rbindlist(list(
    validation,
    data.table::data.table(
      test_case = c(
        "terminal_curve_class_by_seed_rows",
        "terminal_curve_class_summary_rows",
        "terminal_curve_class_pair_count_comparisons"
      ),
      expected = c(
        as.character(data.table::uniqueN(comp$seed_id)),
        "4",
        "3"
      ),
      observed = c(
        as.character(nrow(terminal_class_by_seed)),
        as.character(nrow(terminal_class_summary)),
        as.character(data.table::uniqueN(terminal_class_pair_counts$comparison))
      )
    )
  ), use.names = TRUE, fill = TRUE)
  validation[, passed := expected == observed | test_case %in% c("day_count", "spectral_gap_present")]
  validation[test_case == "day_count", passed := as.integer(observed) > 0L]
  validation[test_case == "spectral_gap_present", passed := observed == "TRUE"]
  validation[, status := ifelse(passed, "PASS", "FAIL")]
  data.table::setcolorder(validation, c("test_case", "status", "expected", "observed", "passed"))

  write_dt(comp, file.path(paths$out_tables, "eigen_vs_initial_ploidy_by_seed_o2_time.tsv"))
  write_dt(summary_by_time, file.path(paths$out_tables, "summary_by_time.tsv"))
  write_dt(summary_by_gap_region, file.path(paths$out_tables, "summary_by_time_and_gap_region.tsv"))
  write_dt(summary_by_gap_decile, file.path(paths$out_tables, "summary_by_time_and_gap_decile.tsv"))
  write_dt(class_comp, file.path(paths$out_tables, "curve_class_comparison_by_seed_time.tsv"))
  write_dt(class_summary, file.path(paths$out_tables, "curve_class_match_summary.tsv"))
  write_dt(terminal_seed, file.path(paths$out_tables, paste0(terminal_prefix, "_seed_level_comparison.tsv")))
  write_dt(top_discrepant, file.path(paths$out_tables, paste0(terminal_prefix, "_top_discrepant_seed_o2.tsv")))
  write_dt(terminal_class_by_seed, file.path(paths$out_tables, paste0(terminal_prefix, "_curve_class_consistency_by_seed.tsv")))
  write_dt(terminal_class_summary, file.path(paths$out_tables, paste0(terminal_prefix, "_curve_class_consistency_summary.tsv")))
  write_dt(terminal_class_pair_counts, file.path(paths$out_tables, paste0(terminal_prefix, "_curve_class_pair_counts.tsv")))
  write_dt(terminal_class_combination_counts, file.path(paths$out_tables, paste0(terminal_prefix, "_curve_class_combination_counts.tsv")))
  write_dt(evidence_summary, file.path(paths$out_tables, "steady_state_support_summary.tsv"))
  write_dt(run_args, file.path(paths$out_tables, "analysis_run_arguments.tsv"))
  write_dt(validation, file.path(paths$out_tables, "validation.tsv"))

  plot_convergence(summary_by_time, paths$out_figures)
  plot_seed_convergence_panels(comp, summary_by_time, paths$out_figures)
  plot_o2_seed_convergence_panels(comp, paths$out_figures)
  plot_o2_seed_convergence_panels_by_class(comp, paths$out_figures, class_col = "eigen_curve_class")
  plot_fraction_within(summary_by_time, paths$out_figures)
  plot_class_match(class_summary, paths$out_figures)
  plot_terminal_class_consistency_summary(terminal_class_summary, paths$out_figures, terminal_day)
  plot_terminal_class_pair_heatmaps(terminal_class_pair_counts, paths$out_figures, terminal_day)
  plot_gap_decile_heatmap(
    summary_by_gap_decile,
    paths$out_figures,
    "abs_delta_2N_minus_eigen_median",
    "median_abs_2N_minus_eigen_by_time_gap_decile",
    "Median |2N - eigen|"
  )
  plot_gap_decile_heatmap(
    summary_by_gap_decile,
    paths$out_figures,
    "abs_delta_4N_minus_eigen_median",
    "median_abs_4N_minus_eigen_by_time_gap_decile",
    "Median |4N - eigen|"
  )
  plot_gap_decile_heatmap(
    summary_by_gap_decile,
    paths$out_figures,
    "abs_delta_4N_minus_2N_median",
    "median_abs_4N_minus_2N_by_time_gap_decile",
    "Median |4N - 2N|"
  )
  plot_terminal_day_heatmaps(comp, paths$out_figures, terminal_day)
  plot_terminal_day_delta_vs_gap(comp, paths$out_figures, terminal_day)

  if (!all(validation$passed)) {
    stop("Validation failed: ", paste(validation$test_case[!validation$passed], collapse = ", "))
  }

  invisible(list(paths = paths, validation = validation))
}

main <- function(args = parse_args()) {
  root <- normalizePath(
    path.expand(args$result_root %||% args$root %||%
                  default_result_root()),
    mustWork = FALSE
  )
  analysis_max_day <- as_int(args$analysis_max_day %||% args$max_day %||% args$time_end, 1000L)
  paths <- analysis_paths(root, analysis_max_day = analysis_max_day)
  paths$out_dir <- normalizePath(path.expand(args$out_dir %||% paths$out_dir), mustWork = FALSE)
  paths$out_tables <- file.path(paths$out_dir, "tables")
  paths$out_figures <- file.path(paths$out_dir, "figures")
  plot_days <- sort(unique(as_num_vec(args$plot_days, c(100, 200, 300, 500, 700, 1000))))
  plot_days <- plot_days[plot_days <= analysis_max_day]
  message("Building eigen vs initial-ploidy comparison under: ", paths$out_dir)
  res <- build_comparison(paths, plot_days = plot_days, analysis_max_day = analysis_max_day)
  message("Completed comparison outputs: ", paths$out_dir)
  invisible(res)
}

if (identical(environment(), globalenv())) {
  main()
}
