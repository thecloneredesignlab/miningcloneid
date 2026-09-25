# Fixed-O2 visualization helpers.
# These functions consume materialized simulation/analysis data frames and write figures only.

.fixed_o2_plot_utils_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[basename(frame_files) == "fixed_o2_plot_utils.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(
  file.path(
    .fixed_o2_plot_utils_dir,
    "..",
    "..",
    "util",
    "o2_supply_demand_map_fixed_o2_format_utils.R"
  ),
  local = TRUE,
  chdir = TRUE
)
rm(.fixed_o2_plot_utils_dir)

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
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
  d <- gap_by_seed
  d <- d[d$trajectory_regime %in% regimes, , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  d_gap <- d[d$spectral_gap > 0, , drop = FALSE]
  d_time <- d[is.finite(d$time_to_10x_days) & d$time_to_10x_days > 0, , drop = FALSE]
  if (!nrow(d_gap) || !nrow(d_time)) return(invisible(FALSE))

  panel_names <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1: dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2: dominant ploidy < 2"
  )
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00")
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
  d <- gap_by_seed
  if (!nrow(d)) return(invisible(FALSE))
  o2_labels <- c("0", "0.1", "0.2", "0.5", "0.75", "1", "2", "3", "4", "5")
  x_centers <- seq_along(o2_labels)
  offsets <- c(mode1 = -0.18, mode2 = 0.18)
  cols <- c(mode1 = "#0072B2", mode2 = "#D55E00")
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
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
  panel_names <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1: dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2: dominant ploidy < 2"
  )
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00")
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
           col = c(cols[["mode1_attractor_dominant_ploidy_ge_2"]], cols[["mode1_attractor_dominant_ploidy_ge_2"]], cols[["mode2_attractor_dominant_ploidy_lt_2"]], cols[["mode2_attractor_dominant_ploidy_lt_2"]]),
           lty = c(2, 1, 2, 1), lwd = 2, bty = "n")
    par(oldpar)
    grDevices::dev.off()
  }
  invisible(TRUE)
}


fo2_plot_mode_seed_stack <- function(attractors, fig_dir) {
  d <- attractors
  if (!nrow(d)) return(invisible(FALSE))
  regimes <- c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2")
  panel_names <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1: dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2: dominant ploidy < 2"
  )
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00")
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
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00", unknown = "grey60")
  d <- attractors[is.finite(attractors$dominant_mean_ploidy), , drop = FALSE]
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


cf2_plot <- function(trajectory, summary_by_seed, fig_dir) {
  d <- trajectory[trajectory$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  grDevices::pdf(file.path(fig_dir, "fixed_o2_counterfactual_median_trajectories.pdf"), width = 14, height = 12, onefile = FALSE, bg = "white")
  median_pdf_open <- TRUE
  on.exit(if (isTRUE(median_pdf_open)) grDevices::dev.off(), add = TRUE)

  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  n_o2 <- length(o2_use)
  o2_cols <- if (n_o2 >= 3L) 3L else max(1L, n_o2)
  o2_rows <- ceiling(n_o2 / o2_cols)
  if (n_o2 == 6L) {
    o2_rows <- 2L
    o2_cols <- 3L
  }
  cols <- c(mode1_attractor_dominant_ploidy_ge_2 = "#0072B2", mode2_attractor_dominant_ploidy_lt_2 = "#D55E00")
  line_cols <- setNames(grDevices::adjustcolor(unname(cols), alpha.f = 0.20), names(cols))
  reg_labels <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "Mode 1 dominant ploidy >= 2",
    mode2_attractor_dominant_ploidy_lt_2 = "Mode 2 dominant ploidy < 2"
  )
  x_rng <- range(d$day, na.rm = TRUE)
  y_rng <- range(d$mean_ploidy, na.rm = TRUE)

  par(mfrow = c(o2_rows * length(init_use), o2_cols), mar = c(2.0, 4.2, 2.7, 1.0), oma = c(3.0, 3.2, 3.6, 0.5))
  for (o2_row in seq_len(o2_rows)) {
    row_o2 <- o2_use[((o2_row - 1L) * o2_cols + 1L):min(o2_row * o2_cols, n_o2)]
    for (init_idx in seq_along(init_use)) {
      init <- init_use[[init_idx]]
      for (col_idx in seq_len(o2_cols)) {
        if (col_idx > length(row_o2)) {
          plot.new()
          next
        }
        O2 <- row_o2[[col_idx]]
        d_panel <- d[d$O2_pct == O2 & d$initial_condition == init, , drop = FALSE]
        show_x <- init_idx == length(init_use)
        show_y <- col_idx == 1L
        plot(
          NA,
          xlim = x_rng,
          ylim = y_rng,
          xlab = "",
          ylab = "",
          xaxt = if (show_x) "s" else "n",
          yaxt = if (show_y) "s" else "n",
          main = if (init_idx == 1L) paste0("O2 = ", O2, "%; initial condition: ", init) else paste("Initial condition:", init),
          cex.main = 0.9
        )
        abline(h = c(1.5, 2), col = "grey80", lty = 2)
        for (reg in names(cols)) {
          sub <- d_panel[d_panel$trajectory_regime == reg, , drop = FALSE]
          if (!nrow(sub)) next
          seed_traces <- split(sub, sub$seed_id)
          for (trace in seed_traces) {
            trace <- trace[order(trace$day), , drop = FALSE]
            lines(trace$day, trace$mean_ploidy, col = line_cols[[reg]], lwd = 0.25)
          }
        }
        if (o2_row == 1L && init_idx == 1L && col_idx == 1L) {
          legend(
            "topright",
            legend = unname(reg_labels[names(cols)]),
            col = unname(cols),
            lty = 1,
            lwd = 2,
            bty = "n",
            cex = 0.72
          )
        }
      }
    }
  }
  mtext("Day", side = 1, outer = TRUE, line = 1.5)
  mtext("Mean ploidy", side = 2, outer = TRUE, line = 1.8)
  mtext("Fixed-O2 seed trajectories", side = 3, outer = TRUE, cex = 1.2, font = 2)
  grDevices::dev.off()
  median_pdf_open <- FALSE

  grDevices::pdf(file.path(fig_dir, "fixed_o2_counterfactual_terminal_boxplots.pdf"), width = 10, height = 7, bg = "white")
  d2 <- summary_by_seed[summary_by_seed$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  boxplot(terminal_mean_ploidy ~ interaction(O2_pct, initial_condition, trajectory_regime, drop = TRUE),
          data = d2, las = 2, ylab = "Terminal mean ploidy", main = "Fixed-O2 counterfactual terminal ploidy")
  grDevices::dev.off()
  invisible(TRUE)
}


vf2_plot <- function(traj, out_path, plot_dt) {
  d <- traj[abs(traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  seed_use <- unique(d$seed_id)
  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  pal <- c("#0072B2", "#D55E00", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  eigen_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(d$seed_label[d$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init <- d[d$initial_condition == init, , drop = FALSE]
    yr <- range(c(sub_init$expm_mean_ploidy, sub_init$euler_mean_ploidy), na.rm = TRUE)
    xr <- range(sub_init$day, na.rm = TRUE)
    for (seed in seed_use) {
      sub_seed <- sub_init[sub_init$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Day",
        ylab = "Mean ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed[sub_seed$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub$expm_mean_ploidy, col = eigen_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub <- sub_seed[sub_seed$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub$euler_mean_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.6, lty = 2)
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), "Expm analytical", paste0("Euler dt=", plot_dt)),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    "Fixed-O2 expm analytical vs Euler integration",
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}


vf2_plot_solution_vs_simulation <- function(solution_traj, simulation_traj, out_path, plot_dt, solution_name, solution_label, title) {
  d <- solution_traj[abs(solution_traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d) || is.null(simulation_traj) || !nrow(simulation_traj)) return(invisible(FALSE))
  solution_col <- paste0(solution_name, "_mean_ploidy")
  if (!solution_col %in% names(d)) stop("solution_traj missing column: ", solution_col)
  seed_use <- unique(d$seed_id)
  o2_use <- sort(unique(d$O2_pct))
  init_use <- unique(d$initial_condition)
  pal <- c("#0072B2", "#D55E00", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  analytical_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(d$seed_label[d$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init_sol <- d[d$initial_condition == init, , drop = FALSE]
    sub_init_sim <- simulation_traj[simulation_traj$initial_condition == init, , drop = FALSE]
    yr <- range(c(sub_init_sol[[solution_col]], sub_init_sim$simulation_mean_ploidy), na.rm = TRUE)
    xr <- range(c(sub_init_sol$day, sub_init_sim$day), na.rm = TRUE)
    for (seed in seed_use) {
      sub_seed_sol <- sub_init_sol[sub_init_sol$seed_id == seed, , drop = FALSE]
      sub_seed_sim <- sub_init_sim[sub_init_sim$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Day",
        ylab = "Mean ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed_sol[sub_seed_sol$O2_pct == O2, , drop = FALSE]
        lines(sub$day, sub[[solution_col]], col = analytical_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub_o2_sim <- sub_seed_sim[sub_seed_sim$O2_pct == O2, , drop = FALSE]
        reps <- split(sub_o2_sim, sub_o2_sim$simulation_id)
        for (rep_dat in reps) {
          rep_dat <- rep_dat[order(rep_dat$day), , drop = FALSE]
          lines(rep_dat$day, rep_dat$simulation_mean_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.2, lty = 2)
        }
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), solution_label, "Simulation reps"),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    title,
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}


vf2_range_pad <- function(x, frac = 0.04) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 1))
  xr <- range(x)
  span <- diff(xr)
  if (!is.finite(span) || span <= 0) {
    pad <- max(abs(xr[[1]]) * frac, 0.05)
  } else {
    pad <- span * frac
  }
  xr + c(-pad, pad)
}


vf2_solution_phase_plane <- function(solution_traj, plot_dt, solution_name, seed_mode_map) {
  d <- solution_traj[abs(solution_traj$dt - plot_dt) < 1e-12, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  mean_col <- paste0(solution_name, "_mean_ploidy")
  sd_col <- paste0(solution_name, "_sd_ploidy")
  if (!all(c(mean_col, sd_col) %in% names(d))) stop("solution_traj missing phase-plane columns for method: ", solution_name)
  out <- data.frame(
    seed_id = d$seed_id,
    seed_mode = unname(seed_mode_map[d$seed_id]),
    seed_label = d$seed_label,
    initial_ploidy = d$initial_ploidy,
    initial_condition = d$initial_condition,
    o2_pct = d$O2_pct,
    simulation_id = NA_integer_,
    solution_method = solution_name,
    time_days = d$day,
    mean_ploidy = d[[mean_col]],
    sd_ploidy = d[[sd_col]],
    dt = d$dt,
    stringsAsFactors = FALSE
  )
  seed_rank <- match(out$seed_id, names(seed_mode_map))
  out[order(seed_rank, out$initial_ploidy, out$o2_pct, out$time_days), , drop = FALSE]
}


vf2_simulation_phase_plane <- function(simulation_traj, seed_mode_map) {
  if (is.null(simulation_traj) || !nrow(simulation_traj)) return(data.frame())
  out <- data.frame(
    seed_id = simulation_traj$seed_id,
    seed_mode = if ("seed_mode" %in% names(simulation_traj)) simulation_traj$seed_mode else unname(seed_mode_map[simulation_traj$seed_id]),
    seed_label = simulation_traj$seed_label,
    initial_ploidy = simulation_traj$initial_ploidy,
    initial_condition = simulation_traj$initial_condition,
    o2_pct = simulation_traj$O2_pct,
    simulation_id = simulation_traj$simulation_id,
    solution_method = "simulation",
    time_days = simulation_traj$day,
    mean_ploidy = simulation_traj$simulation_mean_ploidy,
    sd_ploidy = simulation_traj$simulation_sd_ploidy,
    dt = NA_real_,
    stringsAsFactors = FALSE
  )
  seed_rank <- match(out$seed_id, names(seed_mode_map))
  out[order(seed_rank, out$initial_ploidy, out$o2_pct, out$simulation_id, out$time_days), , drop = FALSE]
}


vf2_plot_phase_plane_solution_vs_simulation <- function(solution_phase, simulation_phase, out_path, solution_label, title) {
  if (is.null(solution_phase) || !nrow(solution_phase) || is.null(simulation_phase) || !nrow(simulation_phase)) {
    return(invisible(FALSE))
  }
  seed_use <- unique(solution_phase$seed_id)
  o2_use <- sort(unique(solution_phase$o2_pct))
  init_use <- unique(solution_phase$initial_condition)
  pal <- c("#0072B2", "#D55E00", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  o2_cols <- setNames(rep(pal, length.out = length(o2_use)), as.character(o2_use))
  analytical_cols <- setNames(grDevices::adjustcolor(unname(o2_cols), alpha.f = 0.45), names(o2_cols))
  seed_label_lookup <- stats::setNames(vapply(seed_use, function(seed) {
    labs <- unique(solution_phase$seed_label[solution_phase$seed_id == seed])
    labs <- labs[nzchar(labs)]
    if (length(labs)) labs[[1]] else seed
  }, character(1)), seed_use)

  grDevices::pdf(out_path, width = if (length(seed_use) > 1L) 14 else 10, height = 8, onefile = TRUE, bg = "white")
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(length(init_use), length(seed_use)), mar = c(4.4, 5, 3, 1), oma = c(0, 0, 3, 0))
  for (init in init_use) {
    sub_init_sol <- solution_phase[solution_phase$initial_condition == init, , drop = FALSE]
    sub_init_sim <- simulation_phase[simulation_phase$initial_condition == init, , drop = FALSE]
    xr <- vf2_range_pad(c(sub_init_sol$mean_ploidy, sub_init_sim$mean_ploidy))
    yr <- vf2_range_pad(c(sub_init_sol$sd_ploidy, sub_init_sim$sd_ploidy))
    for (seed in seed_use) {
      sub_seed_sol <- sub_init_sol[sub_init_sol$seed_id == seed, , drop = FALSE]
      sub_seed_sim <- sub_init_sim[sub_init_sim$seed_id == seed, , drop = FALSE]
      plot(
        NA,
        xlim = xr,
        ylim = yr,
        xlab = "Mean ploidy",
        ylab = "SD ploidy",
        main = paste(seed_label_lookup[[seed]], init, sep = "\n")
      )
      for (O2 in o2_use) {
        sub <- sub_seed_sol[sub_seed_sol$o2_pct == O2, , drop = FALSE]
        sub <- sub[order(sub$time_days), , drop = FALSE]
        lines(sub$mean_ploidy, sub$sd_ploidy, col = analytical_cols[[as.character(O2)]], lwd = 2, lty = 1)
      }
      for (O2 in o2_use) {
        sub_o2_sim <- sub_seed_sim[sub_seed_sim$o2_pct == O2, , drop = FALSE]
        reps <- split(sub_o2_sim, sub_o2_sim$simulation_id)
        for (rep_dat in reps) {
          rep_dat <- rep_dat[order(rep_dat$time_days), , drop = FALSE]
          lines(rep_dat$mean_ploidy, rep_dat$sd_ploidy, col = o2_cols[[as.character(O2)]], lwd = 1.1, lty = 2)
        }
      }
      if (identical(init, init_use[[1]]) && identical(seed, seed_use[[length(seed_use)]])) {
        legend(
          "topright",
          legend = c(paste0("O2=", o2_use, "%"), solution_label, "Simulation reps"),
          col = c(unname(o2_cols[as.character(o2_use)]), "black", "black"),
          lwd = c(rep(2, length(o2_use)), 2, 1.6),
          lty = c(rep(1, length(o2_use)), 1, 2),
          bty = "n",
          cex = 0.78
        )
      }
    }
  }
  mtext(
    title,
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(TRUE)
}


plot_limits <- function(dat) {
  vals <- c(dat$analytical_mean_ploidy, dat$simulation_mean_ploidy)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(0, 1))
  rng <- range(vals)
  pad <- diff(rng) * 0.04
  if (!is.finite(pad) || pad <= 0) pad <- 0.05
  rng + c(-pad, pad)
}


vis_format_o2_label <- format_o2_label


prepare_agreement_plot_data <- function(dat) {
  required <- c(
    "O2_pct", "day", "initial_condition", "analytical_mean_ploidy",
    "simulation_mean_ploidy", "residual_ploidy", "agreement_mean_ploidy",
    "bland_altman_bias", "bland_altman_lower", "bland_altman_upper"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(
      "Agreement visualization table is missing analysis column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  o2_values <- sort(unique(as.numeric(dat$O2_pct)))
  day_values <- sort(unique(as.numeric(dat$day)))
  dat$O2_factor <- factor(
    vis_format_o2_label(dat$O2_pct),
    levels = vis_format_o2_label(o2_values)
  )
  dat$day_factor <- factor(
    paste0("Day ", format(as.numeric(dat$day), scientific = FALSE, trim = TRUE)),
    levels = paste0(
      "Day ",
      format(day_values, scientific = FALSE, trim = TRUE)
    )
  )
  dat$initial_condition <- factor(
    dat$initial_condition,
    levels = sort(unique(as.character(dat$initial_condition)))
  )
  mode <- mode_values(dat)
  dat$mode_factor <- factor(mode, levels = mode_levels(mode))
  if ("objective" %in% names(dat)) dat$objective <- as.numeric(dat$objective)
  dat
}


objective_aesthetic <- function(dat, transform = "identity") {
  objective <- as.numeric(dat$objective)
  label <- "Final objective"
  if (identical(transform, "log10")) {
    objective <- ifelse(is.finite(objective) & objective > 0, log10(objective), NA_real_)
    label <- "log10(final objective)"
  }
  list(value = objective, label = label)
}


mode_values <- function(dat) {
  mode <- rep("unknown", nrow(dat))
  if ("mode_label" %in% names(dat)) {
    mode <- as.character(dat$mode_label)
  } else if ("trajectory_regime" %in% names(dat)) {
    mode <- as.character(dat$trajectory_regime)
  }
  mode[is.na(mode) | !nzchar(mode)] <- "unknown"
  mode
}


mode_levels <- function(mode) {
  preferred <- c("mode1", "mode2", "unknown")
  c(intersect(preferred, unique(mode)), sort(setdiff(unique(mode), preferred)))
}


mode_palette <- function(levels) {
  base <- c(
    mode1 = "#0072B2",
    mode2 = "#D55E00",
    unknown = "#C9C9C9"
  )
  missing <- setdiff(levels, names(base))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), "Dark 3")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[levels]
}


base_scatter <- function(dat, limits, point_size = 0.9, alpha = 0.55) {
  ggplot2::ggplot(dat, ggplot2::aes(x = analytical_mean_ploidy, y = simulation_mean_ploidy)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey55", alpha = 0.28, linewidth = 0.25) +
    ggplot2::coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
    ggplot2::labs(
      x = "Analytical solution mean ploidy",
      y = "Simulation-inferred mean ploidy"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
}


discrete_guides <- function() {
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(shape = 21, color = "grey35", stroke = 0.25, alpha = 1, size = 2.6)),
    shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
  )
}


save_plot <- function(plot, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", device = "pdf")
  invisible(path)
}


plot_time_facets_color_o2 <- function(dat, path, limits) {
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = O2_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::labs(fill = "Fixed O2", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}


plot_time_facets_color_mode <- function(dat, path, limits) {
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mode_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_manual(values = mode_palette(levels(dat$mode_factor)), drop = FALSE, na.value = "grey80") +
    ggplot2::labs(fill = "Mode", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}


plot_time_facets_color_objective <- function(dat, path, limits, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
    )
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}


plot_o2_facets_color_mode <- function(dat, path, limits, title = NULL) {
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mode_factor, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_manual(values = mode_palette(levels(dat$mode_factor)), drop = FALSE, na.value = "grey80") +
    ggplot2::labs(fill = "Mode", shape = "Initial condition", title = title) +
    discrete_guides()
  save_plot(p, path, width = 11 * 2 / 3, height = 7 * 2 / 3)
}


plot_o2_facets_color_objective <- function(dat, path, limits, title = NULL, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition", title = title) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
    )
  save_plot(p, path, width = 11 * 2 / 3, height = 7 * 2 / 3)
}


residual_limits <- function(dat) {
  vals <- as.numeric(dat$residual_ploidy)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(-1, 1))
  max_abs <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 0.05
  pad <- max_abs * 0.06
  c(-(max_abs + pad), max_abs + pad)
}


base_comparison_plot <- function(dat, comparison = c("residual", "bland_altman"), title = NULL) {
  comparison <- match.arg(comparison)
  if (identical(comparison, "residual")) {
    dat$comparison_x <- as.numeric(dat$analytical_mean_ploidy)
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = comparison_x, y = residual_ploidy)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey50", alpha = 0.35, linewidth = 0.3) +
      ggplot2::labs(
        x = "Analytical solution mean ploidy",
        y = "Simulation - analytical mean ploidy",
        title = title
      )
  } else {
    dat$comparison_x <- as.numeric(dat$agreement_mean_ploidy)
    ba <- c(
      bias = dat$bland_altman_bias[[1L]],
      lower = dat$bland_altman_lower[[1L]],
      upper = dat$bland_altman_upper[[1L]]
    )
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = comparison_x, y = residual_ploidy)) +
      ggplot2::geom_hline(yintercept = ba[["bias"]], color = "grey35", alpha = 0.45, linewidth = 0.35) +
      ggplot2::geom_hline(yintercept = ba[c("lower", "upper")], color = "grey45", alpha = 0.35, linetype = 2, linewidth = 0.3) +
      ggplot2::labs(
        x = "Mean of analytical and simulation mean ploidy",
        y = "Simulation - analytical mean ploidy",
        title = title
      )
  }
  p +
    ggplot2::coord_cartesian(ylim = residual_limits(dat), expand = FALSE) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
}


plot_time_facets_color_o2_comparison <- function(dat, path, comparison) {
  p <- base_comparison_plot(dat, comparison = comparison) +
    ggplot2::geom_point(
      ggplot2::aes(fill = O2_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::labs(fill = "Fixed O2", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}


plot_time_facets_color_objective_comparison <- function(dat, path, comparison, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_comparison_plot(dat, comparison = comparison) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
    )
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}


plot_time_facets_color_mode_comparison <- function(dat, path, comparison) {
  p <- base_comparison_plot(dat, comparison = comparison) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mode_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_manual(values = mode_palette(levels(dat$mode_factor)), drop = FALSE, na.value = "grey80") +
    ggplot2::labs(fill = "Mode", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}


plot_o2_facets_color_objective_comparison <- function(dat, path, comparison, title = NULL, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_comparison_plot(dat, comparison = comparison, title = title) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey70", color = "grey35", stroke = 0.25, alpha = 1, size = 2.6))
    )
  save_plot(p, path, width = 11 * 2 / 3, height = 7 * 2 / 3)
}


plot_o2_facets_color_mode_comparison <- function(dat, path, comparison, title = NULL) {
  p <- base_comparison_plot(dat, comparison = comparison, title = title) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mode_factor, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0, color = "transparent"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_manual(values = mode_palette(levels(dat$mode_factor)), drop = FALSE, na.value = "grey80") +
    ggplot2::labs(fill = "Mode", shape = "Initial condition") +
    discrete_guides()
  save_plot(p, path, width = 11 * 2 / 3, height = 7 * 2 / 3)
}


make_comparison_outputs <- function(dat, fig_dir, comparison, prefix, objective_transform = "identity") {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  plot_time_facets_color_o2_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_time_color_o2.pdf")),
    comparison = comparison
  )
  plot_time_facets_color_objective_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_time_color_objective.pdf")),
    comparison = comparison,
    objective_transform = objective_transform
  )
  plot_time_facets_color_mode_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_time_color_mode.pdf")),
    comparison = comparison
  )
  plot_o2_facets_color_objective_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_o2_color_objective_all_times.pdf")),
    comparison = comparison,
    title = "All selected time points",
    objective_transform = objective_transform
  )
  plot_o2_facets_color_mode_comparison(
    dat,
    file.path(fig_dir, paste0(prefix, "_by_o2_color_mode_all_times.pdf")),
    comparison = comparison,
    title = "All selected time points"
  )
  for (day in sort(unique(as.numeric(dat$day)))) {
    day_dat <- dat[abs(as.numeric(dat$day) - day) < 1e-9, , drop = FALSE]
    day_label <- format(day, scientific = FALSE, trim = TRUE)
    plot_o2_facets_color_objective_comparison(
      day_dat,
      file.path(fig_dir, paste0(prefix, "_by_o2_color_objective_day", day_label, ".pdf")),
      comparison = comparison,
      title = paste0("Day ", day_label),
      objective_transform = objective_transform
    )
    plot_o2_facets_color_mode_comparison(
      day_dat,
      file.path(fig_dir, paste0(prefix, "_by_o2_color_mode_day", day_label, ".pdf")),
      comparison = comparison,
      title = paste0("Day ", day_label)
    )
  }
  invisible(fig_dir)
}


plot_analytical_solution_vs_fixed_o2 <- function(analytical, path, analytical_method = "analytical") {
  dat <- analytical
  if (!nrow(dat)) return(invisible(NULL))
  numeric_cols <- intersect(c("O2_pct", "day", "initial_ploidy", "analytical_mean_ploidy"), names(dat))
  for (col in numeric_cols) dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  dat <- dat[is.finite(dat$O2_pct) & is.finite(dat$day) & is.finite(dat$analytical_mean_ploidy), , drop = FALSE]
  if (!nrow(dat)) return(invisible(NULL))
  o2_levels <- sort(unique(dat$O2_pct))
  day_levels <- sort(unique(dat$day))
  dat$O2_factor <- factor(vis_format_o2_label(dat$O2_pct), levels = vis_format_o2_label(o2_levels))
  dat$day_factor <- factor(
    paste0("Day ", format(dat$day, scientific = FALSE, trim = TRUE)),
    levels = paste0("Day ", format(day_levels, scientific = FALSE, trim = TRUE))
  )
  init_levels <- analytical_initial_condition_levels(dat$initial_condition)
  dat$initial_condition <- factor(dat$initial_condition, levels = init_levels)
  dat$seed_initial_group <- interaction(dat$seed_id, dat$initial_condition, drop = TRUE)
  dat <- order_initial_condition_for_drawing(dat)

  summary <- stats::aggregate(
    analytical_mean_ploidy ~ O2_factor + day_factor + initial_condition,
    data = dat,
    FUN = stats::median,
    na.rm = TRUE
  )
  names(summary)[names(summary) == "analytical_mean_ploidy"] <- "median_analytical_mean_ploidy"
  summary <- order_initial_condition_for_drawing(summary)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = O2_factor, y = analytical_mean_ploidy)) +
    ggplot2::geom_line(
      ggplot2::aes(group = seed_initial_group, color = initial_condition),
      alpha = 0.08,
      linewidth = 0.22
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = initial_condition),
      alpha = 0.22,
      size = 0.35,
      stroke = 0
    ) +
    ggplot2::geom_line(
      data = summary,
      ggplot2::aes(y = median_analytical_mean_ploidy, group = initial_condition, color = initial_condition),
      linewidth = 0.75,
      alpha = 0.95
    ) +
    ggplot2::geom_point(
      data = summary,
      ggplot2::aes(y = median_analytical_mean_ploidy, color = initial_condition),
      size = 1.15,
      alpha = 0.95
    ) +
    ggplot2::scale_color_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::labs(
      x = "Fixed O2",
      y = "Analytical solution mean ploidy",
      color = "Initial condition",
      title = paste0(method_label(analytical_method), ": analytical solution across seeds")
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}


analytical_initial_condition_levels <- function(x) {
  vals <- unique(as.character(x))
  preferred <- c("init_2N", "init_4N")
  c(preferred[preferred %in% vals], sort(setdiff(vals, preferred)))
}


analytical_initial_condition_palette <- function(levels) {
  base <- c(
    init_2N = "#1B9E77",
    init_4N = "#7B3294"
  )
  missing <- setdiff(levels, names(base))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[levels]
}


order_initial_condition_for_drawing <- function(dat) {
  if (!"initial_condition" %in% names(dat)) return(dat)
  draw_order <- match(as.character(dat$initial_condition), c("init_4N", "init_2N"))
  draw_order[is.na(draw_order)] <- 0L
  dat[order(draw_order), , drop = FALSE]
}


add_analytical_mode_labels <- function(analytical, seed_metadata = NULL) {
  dat <- analytical
  if (!nrow(dat) || !"seed_id" %in% names(dat)) return(dat)
  if (!"trajectory_regime" %in% names(dat) && is.data.frame(seed_metadata) && nrow(seed_metadata) && "seed_id" %in% names(seed_metadata)) {
    cols <- intersect(c(
      "seed_id", "O2_key", "trajectory_regime", "mode_label", "mode_source", "mode_rule",
      "mode_threshold_dominant_ploidy", "mode_reference_o2_pct", "mode_reference_o2_key",
      "mode_reference_dominant_mean_ploidy", "mode_reference_status",
      "mode_reference_dominant_growth_rate", "mode_reference_spectral_gap"
    ), names(seed_metadata))
    if (all(c("seed_id", "trajectory_regime") %in% cols)) {
      meta <- seed_metadata[, cols, drop = FALSE]
      if ("O2_key" %in% names(meta) && "O2_key" %in% names(dat)) {
        meta <- meta[!duplicated(meta[, c("seed_id", "O2_key"), drop = FALSE]), , drop = FALSE]
        dat <- merge(dat, meta, by = c("seed_id", "O2_key"), all.x = TRUE, sort = FALSE)
      } else {
        meta <- meta[!duplicated(meta$seed_id), , drop = FALSE]
        dat <- merge(dat, meta, by = "seed_id", all.x = TRUE, sort = FALSE)
      }
    }
  }
  dat
}


plot_analytical_solution_vs_fixed_o2_by_mode <- function(analytical, path, analytical_method = "analytical", seed_metadata = NULL) {
  dat <- add_analytical_mode_labels(analytical, seed_metadata = seed_metadata)
  if (!nrow(dat) || !"trajectory_regime" %in% names(dat)) return(invisible(NULL))
  numeric_cols <- intersect(c("O2_pct", "day", "initial_ploidy", "analytical_mean_ploidy"), names(dat))
  for (col in numeric_cols) dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  dat <- dat[is.finite(dat$O2_pct) & is.finite(dat$day) & is.finite(dat$analytical_mean_ploidy), , drop = FALSE]
  dat <- dat[dat$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  if (!nrow(dat)) return(invisible(NULL))

  o2_levels <- sort(unique(dat$O2_pct))
  day_levels <- sort(unique(dat$day))
  day_levels_label <- paste0("Day ", format(day_levels, scientific = FALSE, trim = TRUE))
  mode_labels <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "mode1",
    mode2_attractor_dominant_ploidy_lt_2 = "mode2"
  )
  dat$O2_factor <- factor(vis_format_o2_label(dat$O2_pct), levels = vis_format_o2_label(o2_levels))
  init_levels <- analytical_initial_condition_levels(dat$initial_condition)
  dat$initial_condition <- factor(dat$initial_condition, levels = init_levels)
  dat$mode_group <- factor(mode_labels[as.character(dat$trajectory_regime)], levels = c("mode1", "mode2"))
  dat$day_label <- factor(
    paste0("Day ", format(dat$day, scientific = FALSE, trim = TRUE)),
    levels = day_levels_label
  )
  dat$facet_label <- factor(
    paste(as.character(dat$day_label), as.character(dat$mode_group), sep = "\n"),
    levels = as.vector(rbind(
      paste(day_levels_label, "mode1", sep = "\n"),
      paste(day_levels_label, "mode2", sep = "\n")
    ))
  )
  dat$seed_initial_group <- interaction(dat$seed_id, dat$initial_condition, dat$mode_group, drop = TRUE)
  dat <- order_initial_condition_for_drawing(dat)

  summary <- stats::aggregate(
    analytical_mean_ploidy ~ O2_factor + facet_label + initial_condition,
    data = dat,
    FUN = stats::median,
    na.rm = TRUE
  )
  names(summary)[names(summary) == "analytical_mean_ploidy"] <- "median_analytical_mean_ploidy"
  summary <- order_initial_condition_for_drawing(summary)
  panel_background <- unique(dat[, c("facet_label", "mode_group"), drop = FALSE])

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = O2_factor, y = analytical_mean_ploidy)) +
    ggplot2::geom_rect(
      data = panel_background,
      ggplot2::aes(fill = mode_group),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE,
      alpha = 0.35,
      color = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(group = seed_initial_group, color = initial_condition),
      alpha = 0.08,
      linewidth = 0.2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = initial_condition),
      alpha = 0.22,
      size = 0.32,
      stroke = 0
    ) +
    ggplot2::geom_line(
      data = summary,
      ggplot2::aes(y = median_analytical_mean_ploidy, group = initial_condition, color = initial_condition),
      linewidth = 0.7,
      alpha = 0.95
    ) +
    ggplot2::geom_point(
      data = summary,
      ggplot2::aes(y = median_analytical_mean_ploidy, color = initial_condition),
      size = 1.05,
      alpha = 0.95
    ) +
    ggplot2::scale_color_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::scale_fill_manual(
      values = c(mode1 = "#DCEEFF", mode2 = "#FFE8CC"),
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::facet_wrap(~facet_label, nrow = 2, ncol = 8) +
    ggplot2::labs(
      x = "Fixed O2",
      y = "Analytical solution mean ploidy",
      color = "Initial condition",
      title = paste0(method_label(analytical_method), ": analytical solution by mode")
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
  save_plot(p, path, width = 16, height = 7)
}


half_violin_polygon_data <- function(dat, x_col, y_col, side_col, group_cols, width = 0.36, n = 128L) {
  if (!nrow(dat)) return(data.frame())
  split_cols <- unique(c(group_cols, side_col))
  split_key <- interaction(dat[, split_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(dat, split_key), function(d) {
    y <- suppressWarnings(as.numeric(d[[y_col]]))
    y <- y[is.finite(y)]
    if (length(y) < 2L || length(unique(y)) < 2L) return(NULL)
    dens <- stats::density(y, n = n, from = min(y), to = max(y), na.rm = TRUE)
    if (!length(dens$y) || !is.finite(max(dens$y)) || max(dens$y) <= 0) return(NULL)
    side <- as.character(d[[side_col]][[1]])
    direction <- if (identical(side, "init_2N")) 1 else -1
    center <- suppressWarnings(as.numeric(d[[x_col]][[1]]))
    if (!is.finite(center)) return(NULL)
    outer_x <- center + direction * width * dens$y / max(dens$y)
    out <- data.frame(
      x = c(rep(center, length(dens$x)), rev(outer_x)),
      y = c(dens$x, rev(dens$x)),
      stringsAsFactors = FALSE
    )
    for (col in split_cols) out[[col]] <- d[[col]][[1]]
    out$violin_group <- paste(vapply(split_cols, function(col) as.character(d[[col]][[1]]), character(1)), collapse = "\r")
    out
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out
}


half_violin_median_data <- function(dat, x_col, y_col, side_col, group_cols, offset = 0.16) {
  if (!nrow(dat)) return(data.frame())
  split_cols <- unique(c(group_cols, side_col))
  rows <- lapply(split(dat, interaction(dat[, split_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)), function(d) {
    y <- suppressWarnings(as.numeric(d[[y_col]]))
    y <- y[is.finite(y)]
    if (!length(y)) return(NULL)
    side <- as.character(d[[side_col]][[1]])
    direction <- if (identical(side, "init_2N")) 1 else -1
    center <- suppressWarnings(as.numeric(d[[x_col]][[1]]))
    if (!is.finite(center)) return(NULL)
    out <- data.frame(
      x = center + direction * offset,
      x0 = center + direction * offset - 0.06,
      x1 = center + direction * offset + 0.06,
      median_y = stats::median(y, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for (col in split_cols) out[[col]] <- d[[col]][[1]]
    out
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out
}


plot_analytical_solution_vs_fixed_o2_init_half_violin <- function(analytical, path, analytical_method = "analytical") {
  dat <- analytical
  if (!nrow(dat)) return(invisible(NULL))
  numeric_cols <- intersect(c("O2_pct", "day", "initial_ploidy", "analytical_mean_ploidy"), names(dat))
  for (col in numeric_cols) dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  dat <- dat[is.finite(dat$O2_pct) & is.finite(dat$day) & is.finite(dat$analytical_mean_ploidy), , drop = FALSE]
  if (!nrow(dat)) return(invisible(NULL))

  o2_levels <- sort(unique(dat$O2_pct))
  day_levels <- sort(unique(dat$day))
  init_levels <- analytical_initial_condition_levels(dat$initial_condition)
  dat$O2_factor <- factor(vis_format_o2_label(dat$O2_pct), levels = vis_format_o2_label(o2_levels))
  dat$initial_condition <- factor(dat$initial_condition, levels = init_levels)
  dat$day_factor <- factor(
    paste0("Day ", format(dat$day, scientific = FALSE, trim = TRUE)),
    levels = paste0("Day ", format(day_levels, scientific = FALSE, trim = TRUE))
  )
  x_levels <- data.frame(
    O2_factor = levels(dat$O2_factor),
    x_pos = seq_along(levels(dat$O2_factor)),
    stringsAsFactors = FALSE
  )
  dat <- merge(dat, x_levels, by = "O2_factor", all.x = TRUE, sort = FALSE)
  dat <- order_initial_condition_for_drawing(dat)

  violin <- half_violin_polygon_data(
    dat = dat,
    x_col = "x_pos",
    y_col = "analytical_mean_ploidy",
    side_col = "initial_condition",
    group_cols = c("day_factor", "O2_factor")
  )
  if (!nrow(violin)) return(invisible(NULL))
  violin$initial_condition <- factor(violin$initial_condition, levels = init_levels)
  med <- half_violin_median_data(
    dat = dat,
    x_col = "x_pos",
    y_col = "analytical_mean_ploidy",
    side_col = "initial_condition",
    group_cols = c("day_factor", "O2_factor")
  )
  med$initial_condition <- factor(med$initial_condition, levels = init_levels)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = violin,
      ggplot2::aes(x = x, y = y, group = violin_group, fill = initial_condition),
      color = NA,
      alpha = 0.62
    ) +
    ggplot2::geom_segment(
      data = med,
      ggplot2::aes(x = x0, xend = x1, y = median_y, yend = median_y, color = initial_condition),
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      data = med,
      ggplot2::aes(x = x, y = median_y, color = initial_condition),
      fill = "white",
      shape = 21,
      size = 1.15,
      stroke = 0.45
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_levels$x_pos,
      labels = x_levels$O2_factor,
      expand = ggplot2::expansion(mult = c(0.03, 0.03))
    ) +
    ggplot2::scale_fill_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::scale_color_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::labs(
      x = "Fixed O2",
      y = "Analytical solution mean ploidy",
      fill = "Initial condition",
      color = "Initial condition",
      title = paste0(method_label(analytical_method), ": analytical solution distribution across fixed O2")
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
  save_plot(p, path, width = 13 * 2 / 3, height = 7 * 2 / 3)
}


plot_analytical_solution_vs_fixed_o2_mode_init_half_violin <- function(analytical, path, analytical_method = "analytical", seed_metadata = NULL) {
  dat <- add_analytical_mode_labels(analytical, seed_metadata = seed_metadata)
  if (!nrow(dat) || !"trajectory_regime" %in% names(dat)) return(invisible(NULL))
  numeric_cols <- intersect(c("O2_pct", "day", "initial_ploidy", "analytical_mean_ploidy"), names(dat))
  for (col in numeric_cols) dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
  dat <- dat[is.finite(dat$O2_pct) & is.finite(dat$day) & is.finite(dat$analytical_mean_ploidy), , drop = FALSE]
  dat <- dat[dat$trajectory_regime %in% c("mode1_attractor_dominant_ploidy_ge_2", "mode2_attractor_dominant_ploidy_lt_2"), , drop = FALSE]
  if (!nrow(dat)) return(invisible(NULL))

  o2_levels <- sort(unique(dat$O2_pct))
  day_levels <- sort(unique(dat$day))
  day_levels_label <- paste0("Day ", format(day_levels, scientific = FALSE, trim = TRUE))
  mode_labels <- c(
    mode1_attractor_dominant_ploidy_ge_2 = "mode1",
    mode2_attractor_dominant_ploidy_lt_2 = "mode2"
  )
  init_levels <- analytical_initial_condition_levels(dat$initial_condition)

  dat$O2_factor <- factor(vis_format_o2_label(dat$O2_pct), levels = vis_format_o2_label(o2_levels))
  dat$initial_condition <- factor(dat$initial_condition, levels = init_levels)
  dat$mode_group <- factor(mode_labels[as.character(dat$trajectory_regime)], levels = c("mode1", "mode2"))
  dat$day_factor <- factor(
    paste0("Day ", format(dat$day, scientific = FALSE, trim = TRUE)),
    levels = day_levels_label
  )
  x_levels <- expand.grid(
    mode_group = levels(dat$mode_group),
    O2_factor = levels(dat$O2_factor),
    stringsAsFactors = FALSE
  )
  x_levels <- x_levels[, c("O2_factor", "mode_group"), drop = FALSE]
  x_levels$x_pos <- seq_len(nrow(x_levels))
  x_levels$x_label <- paste(x_levels$O2_factor, x_levels$mode_group, sep = "\n")
  dat <- merge(dat, x_levels[, c("O2_factor", "mode_group", "x_pos"), drop = FALSE],
               by = c("O2_factor", "mode_group"), all.x = TRUE, sort = FALSE)
  dat <- order_initial_condition_for_drawing(dat)

  violin <- half_violin_polygon_data(
    dat = dat,
    x_col = "x_pos",
    y_col = "analytical_mean_ploidy",
    side_col = "initial_condition",
    group_cols = c("day_factor", "O2_factor", "mode_group")
  )
  if (!nrow(violin)) return(invisible(NULL))
  violin$initial_condition <- factor(violin$initial_condition, levels = init_levels)
  med <- half_violin_median_data(
    dat = dat,
    x_col = "x_pos",
    y_col = "analytical_mean_ploidy",
    side_col = "initial_condition",
    group_cols = c("day_factor", "O2_factor", "mode_group")
  )
  med$initial_condition <- factor(med$initial_condition, levels = init_levels)

  bg <- unique(dat[, c("day_factor", "x_pos", "mode_group"), drop = FALSE])
  bg$xmin <- bg$x_pos - 0.5
  bg$xmax <- bg$x_pos + 0.5
  bg_mode1 <- bg[bg$mode_group == "mode1", , drop = FALSE]
  bg_mode2 <- bg[bg$mode_group == "mode2", , drop = FALSE]

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = bg_mode1,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#DCEEFF",
      alpha = 0.42,
      color = NA
    ) +
    ggplot2::geom_rect(
      data = bg_mode2,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#FFE8CC",
      alpha = 0.42,
      color = NA
    ) +
    ggplot2::geom_polygon(
      data = violin,
      ggplot2::aes(x = x, y = y, group = violin_group, fill = initial_condition),
      color = NA,
      alpha = 0.58
    ) +
    ggplot2::geom_segment(
      data = med,
      ggplot2::aes(x = x0, xend = x1, y = median_y, yend = median_y, color = initial_condition),
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      data = med,
      ggplot2::aes(x = x, y = median_y, color = initial_condition),
      fill = "white",
      shape = 21,
      size = 1.05,
      stroke = 0.45
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_levels$x_pos,
      labels = x_levels$x_label,
      expand = ggplot2::expansion(mult = c(0.015, 0.015))
    ) +
    ggplot2::scale_fill_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::scale_color_manual(values = analytical_initial_condition_palette(init_levels), breaks = init_levels, drop = FALSE) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::labs(
      x = "Fixed O2 and mode",
      y = "Analytical solution mean ploidy",
      fill = "Initial condition",
      color = "Initial condition",
      title = paste0(method_label(analytical_method), ": analytical solution distribution by fixed O2, mode, and initial condition")
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6.5),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
  save_plot(p, path, width = 16, height = 8.5)
}


make_scatter_outputs <- function(dat, out_dir, objective_transform = "identity", analytical_method = "analytical") {
  dat <- prepare_agreement_plot_data(dat)
  method_dir <- file.path(out_dir, method_slug(analytical_method))
  scatter_dir <- file.path(method_dir, "scatter")
  residual_dir <- file.path(method_dir, "residual")
  bland_altman_dir <- file.path(method_dir, "bland_altman")
  dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)
  limits <- plot_limits(dat)

  plot_time_facets_color_o2(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_time_color_o2.pdf"),
    limits = limits
  )
  plot_time_facets_color_objective(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_time_color_objective.pdf"),
    limits = limits,
    objective_transform = objective_transform
  )
  plot_time_facets_color_mode(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_time_color_mode.pdf"),
    limits = limits
  )
  plot_o2_facets_color_objective(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_o2_color_objective_all_times.pdf"),
    limits = limits,
    title = "All selected time points",
    objective_transform = objective_transform
  )
  plot_o2_facets_color_mode(
    dat,
    file.path(scatter_dir, "scatter_analytical_vs_simulation_by_o2_color_mode_all_times.pdf"),
    limits = limits,
    title = "All selected time points"
  )

  for (day in sort(unique(as.numeric(dat$day)))) {
    day_dat <- dat[abs(as.numeric(dat$day) - day) < 1e-9, , drop = FALSE]
    day_label <- format(day, scientific = FALSE, trim = TRUE)
    plot_o2_facets_color_objective(
      day_dat,
      file.path(scatter_dir, paste0("scatter_analytical_vs_simulation_by_o2_color_objective_day", day_label, ".pdf")),
      limits = limits,
      title = paste0("Day ", day_label),
      objective_transform = objective_transform
    )
    plot_o2_facets_color_mode(
      day_dat,
      file.path(scatter_dir, paste0("scatter_analytical_vs_simulation_by_o2_color_mode_day", day_label, ".pdf")),
      limits = limits,
      title = paste0("Day ", day_label)
    )
  }
  make_comparison_outputs(
    dat = dat,
    fig_dir = residual_dir,
    comparison = "residual",
    prefix = "residual_simulation_minus_analytical",
    objective_transform = objective_transform
  )
  make_comparison_outputs(
    dat = dat,
    fig_dir = bland_altman_dir,
    comparison = "bland_altman",
    prefix = "bland_altman_simulation_vs_analytical",
    objective_transform = objective_transform
  )
  invisible(method_dir)
}


make_analytical_solution_outputs <- function(analytical, out_dir, analytical_method = "analytical", seed_metadata = NULL) {
  method_dir <- file.path(out_dir, method_slug(analytical_method))
  plot_dir <- file.path(method_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_analytical_solution_vs_fixed_o2(
    analytical = analytical,
    path = file.path(plot_dir, "analytical_solution_vs_fixed_o2_by_time.pdf"),
    analytical_method = analytical_method
  )
  plot_analytical_solution_vs_fixed_o2_init_half_violin(
    analytical = analytical,
    path = file.path(plot_dir, "analytical_solution_vs_fixed_o2_by_time_init_half_violin_median.pdf"),
    analytical_method = analytical_method
  )
  plot_analytical_solution_vs_fixed_o2_by_mode(
    analytical = analytical,
    path = file.path(plot_dir, "analytical_solution_vs_fixed_o2_by_time_mode1_mode2.pdf"),
    analytical_method = analytical_method,
    seed_metadata = seed_metadata
  )
  plot_analytical_solution_vs_fixed_o2_mode_init_half_violin(
    analytical = analytical,
    path = file.path(plot_dir, "analytical_solution_vs_fixed_o2_by_time_mode_init_half_violin_median.pdf"),
    analytical_method = analytical_method,
    seed_metadata = seed_metadata
  )
  invisible(plot_dir)
}
