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
  normalizePath(file.path(script_dir(), "..", "..", "..", "..", ".."), mustWork = FALSE)
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

analysis_paths <- function(root) {
  mono_dir <- file.path(root, "dense-grid_monotonicity_classification")
  init_dir <- file.path(root, "dense-grid_initial-ploidy_trajectory")
  out_dir <- file.path(root, "compare")
  list(
    root = root,
    mono_dir = mono_dir,
    init_dir = init_dir,
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

plot_o2_seed_convergence_panels <- function(comp, figure_dir, o2_values = c(0, 0.1, 0.5, 1, 2, 5)) {
  available_o2 <- sort(unique(comp$O2_pct[is.finite(comp$O2_pct)]))
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
  yr <- safe_range(y_all[is.finite(y_all)], c(0, 1))
  x_range <- range(d0$day, na.rm = TRUE)
  o2_labels <- paste0("O2 = ", format(matched_o2, scientific = FALSE, trim = TRUE), "%")

  save_plot_pair("convergence_to_eigen_by_time_o2_seed_curves_panels", figure_dir, {
    op <- par(mfrow = c(3, length(matched_o2)),
              mar = c(1.0, 1.0, 1.2, 0.4),
              oma = c(3.4, 6.4, 1.2, 0.2),
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
        for (seed in unique(z$seed_id[order(z$seed_number)])) {
          zz <- z[seed_id == seed][order(day)]
          gap_col <- if (isTRUE(any(zz$eigen_spectral_gap < 0.01, na.rm = TRUE))) {
            grDevices::adjustcolor("#CC79A7", alpha.f = 0.24)
          } else {
            grDevices::adjustcolor("#009E73", alpha.f = 0.24)
          }
          lines(zz$day, zz[[panel$value_col]],
                col = gap_col, lwd = 0.42)
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
  }, width = 14.5, height = 9.5)
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

plot_day1000_heatmaps <- function(comp, figure_dir) {
  d <- comp[day == 1000]
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
  save_plot_pair("day1000_signed_delta_heatmaps", figure_dir, {
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

plot_day1000_delta_vs_gap <- function(comp, figure_dir) {
  d <- comp[day == 1000 & is.finite(eigen_spectral_gap)]
  if (!nrow(d)) return(character())
  if (nrow(d) > 60000L) {
    set.seed(1L)
    d <- d[sample.int(nrow(d), 60000L)]
  }
  save_plot_pair("day1000_abs_delta_vs_eigen_spectral_gap", figure_dir, {
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

build_comparison <- function(paths, plot_days = c(100, 200, 300, 500, 700, 1000)) {
  require_data_table()
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

  day1000_seed <- comp[day == 1000, .(
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

  top_discrepant <- comp[day == 1000][order(-pmax(abs_delta_2N_minus_eigen, abs_delta_4N_minus_eigen, abs_delta_4N_minus_2N))]
  top_discrepant <- top_discrepant[seq_len(min(200L, nrow(top_discrepant)))]

  day1000_global <- summary_by_time[day == 1000]
  day1000_large_gap <- summary_by_gap_region[day == 1000 & eigen_gap_region == "gap_ge_0p01"]
  day1000_small_gap <- summary_by_gap_region[day == 1000 & eigen_gap_region == "gap_lt_0p005"]
  class_day1000 <- class_summary[day == 1000]
  evidence_summary <- data.table::data.table(
    statement = c(
      "Median finite-time 2N-start ploidy reaches the eigen steady-state by day 1000 across seed-O2 pairs.",
      "Median finite-time 4N-start ploidy reaches the eigen steady-state by day 1000 across seed-O2 pairs.",
      "2N-start and 4N-start finite-time trajectories largely agree by day 1000 across seed-O2 pairs.",
      "Large-gap seed-O2 pairs are essentially converged to the eigen steady-state by day 1000.",
      "Small-gap seed-O2 pairs retain most of the day-1000 finite-time/eigen discrepancy.",
      "Finite-time/eigen discrepancy is negatively associated with eigen spectral gap.",
      "Curve-class labels for 2N and 4N finite-time curves are mostly consistent by day 1000.",
      "Curve-class labels do not fully match the eigen classification even when ploidy values mostly converge, because finite-time curve shape can still differ in small-gap/slow-relaxation regions."
    ),
    evidence_metric = c(
      "day1000 median |2N - eigen|",
      "day1000 median |4N - eigen|",
      "day1000 median |4N - 2N|",
      "day1000 large-gap median |2N - eigen|; median |4N - eigen|",
      "day1000 small-gap median |2N - eigen|; median |4N - eigen|",
      "day1000 Spearman(|delta|, eigen spectral gap)",
      "day1000 fraction 2N class = 4N class",
      "day1000 fraction 2N class = eigen class; fraction 4N class = eigen class"
    ),
    value = c(
      format_num(day1000_global$abs_delta_2N_minus_eigen_median, 8L),
      format_num(day1000_global$abs_delta_4N_minus_eigen_median, 8L),
      format_num(day1000_global$abs_delta_4N_minus_2N_median, 8L),
      paste(format_num(day1000_large_gap$abs_delta_2N_minus_eigen_median, 8L),
            format_num(day1000_large_gap$abs_delta_4N_minus_eigen_median, 8L), sep = "; "),
      paste(format_num(day1000_small_gap$abs_delta_2N_minus_eigen_median, 8L),
            format_num(day1000_small_gap$abs_delta_4N_minus_eigen_median, 8L), sep = "; "),
      paste(format_num(day1000_global$spearman_abs_2N_eigen_vs_gap, 4L),
            format_num(day1000_global$spearman_abs_4N_eigen_vs_gap, 4L), sep = "; "),
      format_num(class_day1000$fraction_init_2N_matches_init_4N, 4L),
      paste(format_num(class_day1000$fraction_init_2N_matches_eigen, 4L),
            format_num(class_day1000$fraction_init_4N_matches_eigen, 4L), sep = "; ")
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
      "script", "result_root", "monotonicity_dir", "initial_ploidy_dir", "output_dir",
      "n_seed", "n_o2", "selected_days", "plot_days", "n_comparison_rows",
      "classification_label_eigen", "classification_label_initial_2N", "classification_label_initial_4N",
      "steady_state_interpretation"
    ),
    value = c(
      script_path(), paths$root, paths$mono_dir, paths$init_dir, paths$out_dir,
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
      "day1000_present",
      "spectral_gap_present",
      "validation_tables_passed"
    ),
    expected = c(
      as.character(expected_rows),
      "500",
      "201",
      as.character(data.table::uniqueN(comp$day)),
      "TRUE", "TRUE", "0", "0", "0", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE"
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
      as.character(any(comp$day == 1000)),
      as.character(any(is.finite(comp$eigen_spectral_gap))),
      as.character({
        ev <- read_dt(paths$eigen_validation)
        iv <- read_dt(paths$initial_validation)
        all(ev$passed %in% TRUE | ev$passed == "TRUE") && all(iv$passed %in% TRUE | iv$passed == "TRUE")
      })
    )
  )
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
  write_dt(day1000_seed, file.path(paths$out_tables, "day1000_seed_level_comparison.tsv"))
  write_dt(top_discrepant, file.path(paths$out_tables, "day1000_top_discrepant_seed_o2.tsv"))
  write_dt(evidence_summary, file.path(paths$out_tables, "steady_state_support_summary.tsv"))
  write_dt(run_args, file.path(paths$out_tables, "analysis_run_arguments.tsv"))
  write_dt(validation, file.path(paths$out_tables, "validation.tsv"))

  plot_convergence(summary_by_time, paths$out_figures)
  plot_seed_convergence_panels(comp, summary_by_time, paths$out_figures)
  plot_o2_seed_convergence_panels(comp, paths$out_figures)
  plot_fraction_within(summary_by_time, paths$out_figures)
  plot_class_match(class_summary, paths$out_figures)
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
  plot_day1000_heatmaps(comp, paths$out_figures)
  plot_day1000_delta_vs_gap(comp, paths$out_figures)
  plot_representative_curves(comp, day1000_seed, paths$out_figures, plot_days)

  if (!all(validation$passed)) {
    stop("Validation failed: ", paste(validation$test_case[!validation$passed], collapse = ", "))
  }

  invisible(list(paths = paths, validation = validation))
}

main <- function(args = parse_args()) {
  root <- normalizePath(
    path.expand(args$result_root %||% args$root %||%
                  file.path(repo_root(), "oxygen", "results", "analysis", "monotonicity_classification")),
    mustWork = FALSE
  )
  paths <- analysis_paths(root)
  paths$out_dir <- normalizePath(path.expand(args$out_dir %||% paths$out_dir), mustWork = FALSE)
  paths$out_tables <- file.path(paths$out_dir, "tables")
  paths$out_figures <- file.path(paths$out_dir, "figures")
  plot_days <- sort(unique(as_num_vec(args$plot_days, c(100, 200, 300, 500, 700, 1000))))
  message("Building eigen vs initial-ploidy comparison under: ", paths$out_dir)
  res <- build_comparison(paths, plot_days = plot_days)
  message("Completed comparison outputs: ", paths$out_dir)
  invisible(res)
}

if (identical(environment(), globalenv())) {
  main()
}
