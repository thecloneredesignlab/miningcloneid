#!/usr/bin/env Rscript

# Day-1000 trajectory-level positive-growth conditional mean, with the
# unfiltered Figure 6B mean retained as the gray reference.

f6gp_filename <- function(family) {
  number <- if (identical(family, "C01")) "20" else if (
    identical(family, "C02")) "21" else stop("Unknown Figure 6 family.")
  sprintf("supp_fig6-%s_growth_permissive_ploidy_vs_oxygen_%s", number,
          tolower(family))
}

f6gp_draw_page <- function(curves, family) {
  kept_color <- "#176A86"
  all_color <- "#A5AAB0"
  grid_color <- "#E6E9EA"
  ink <- "#222B30"
  p_values <- f6gp_p_values()
  contexts <- f6ft_context_levels()
  n_columns <- length(p_values)
  upper <- seq.int(1L, by = 2L, length.out = n_columns)
  lower <- upper + 2L * n_columns
  layout(rbind(
    upper, upper + 1L,
    lower, lower + 1L,
    rep(4L * n_columns + 1L, n_columns)
  ),
  heights = c(3.7, 0.82, 3.7, 0.82, 0.8),
  widths = rep(1, n_columns))
  par(oma = c(2.0, 2.3, 2.2, 0.45), family = "sans", fg = ink)
  for (context_index in seq_along(contexts)) {
    context <- contexts[[context_index]]
    context_rows <- curves[curves$model_context == context &
                             curves$pair_label == family, , drop = FALSE]
    for (p_index in seq_along(p_values)) {
      p <- p_values[[p_index]]
      panel <- context_rows[abs(context_rows$p_misseg - p) < 1e-12,
                            , drop = FALSE]
      panel <- panel[order(panel$initial_ploidy, panel$O2_pct), , drop = FALSE]
      if (nrow(panel) != 2L * length(f6gp_o2_values())) {
        stop("Incomplete growth-permissive panel: ", context, " ", family)
      }
      par(mar = c(0.4, if (p_index == 1L) 3.9 else 1.0, 1.9, 0.8),
          mgp = c(2.2, 0.55, 0), tcl = -0.23, xaxs = "i", yaxs = "i")
      plot.new()
      plot.window(xlim = c(0, 5), ylim = c(1, 7.1))
      abline(h = 2:7, col = grid_color, lwd = 0.65)
      if (p_index == 1L) {
        axis(2, at = 1:7, las = 1, cex.axis = 0.86, col.axis = ink)
        mtext("Mean ploidy at day 1000 (N)", side = 2, line = 2.8,
              cex = 0.83, col = ink)
      }
      axis(1, at = 0:5, labels = FALSE, tck = -0.015)
      box(col = "#7D858A", lwd = 0.75)
      for (initial in f6gp_initial_values()) {
        trace <- panel[panel$initial_ploidy == initial, , drop = FALSE]
        trace <- trace[order(trace$O2_pct), , drop = FALSE]
        style <- if (initial == 2) 1 else 2
        lines(trace$O2_pct, trace$unfiltered_mean_ploidy,
              col = all_color, lty = style, lwd = 1.55)
        lines(trace$O2_pct, trace$growth_positive_mean_ploidy,
              col = kept_color, lty = style, lwd = 2.5)
      }
      if (context_index == 1L) {
        mtext(paste0("p_misseg = ", format(p, scientific = FALSE,
                                          trim = TRUE)), side = 3, line = 0.52,
              cex = 0.92, font = 2, col = ink)
      }
      if (p_index == 1L) {
        text(0.12, 6.87, context, adj = c(0, 1), font = 2,
             cex = 0.92, col = ink)
      }
      par(mar = c(if (context_index == 2L) 2.0 else 1.5,
                  if (p_index == 1L) 3.9 else 1.0,
                  0.25, 0.8),
          mgp = c(2.1, 0.55, 0), tcl = -0.23, xaxs = "i", yaxs = "i")
      plot.new()
      plot.window(xlim = c(0, 5), ylim = c(0, 1.06))
      abline(h = c(0, 0.5, 1), col = grid_color, lwd = 0.65)
      for (initial in f6gp_initial_values()) {
        fraction <- panel[panel$initial_ploidy == initial,
                          c("O2_pct", "growth_positive_fraction"),
                          drop = FALSE]
        fraction <- fraction[order(fraction$O2_pct), , drop = FALSE]
        lines(fraction$O2_pct, fraction$growth_positive_fraction,
              type = "s", lty = if (initial == 2) 1 else 2,
              lwd = 2.1, col = kept_color)
      }
      axis(1, at = 0:5, cex.axis = 0.76, col.axis = ink)
      if (p_index == 1L) {
        axis(2, at = c(0, 1), las = 1, cex.axis = 0.76, col.axis = ink)
        mtext("Growth > 0 fraction", side = 2, line = 2.2, cex = 0.76,
              col = ink)
      }
      box(col = "#7D858A", lwd = 0.75)
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", inset = 0, horiz = TRUE, xpd = NA,
         legend = c("Growth > 0: 2N", "Growth > 0: 4N",
                    "All: 2N", "All: 4N"),
         col = c(kept_color, kept_color, all_color, all_color),
         lty = c(1, 2, 1, 2), lwd = c(2.5, 2.5, 1.55, 1.55),
         cex = 0.88, bty = "n", seg.len = 2.8)
  mtext(sprintf(
    "Day-1000 mean ploidy versus oxygen, conditional on positive net growth (%s)",
    family
  ), side = 3, outer = TRUE, line = 0.85,
        font = 2, cex = 1.14, col = ink)
  mtext("Fixed O2 (%)", side = 1, outer = TRUE, line = 0.65,
        cex = 0.94, col = ink)
}

f6gp_draw <- function(workspace_root, family) {
  paths <- f6r_paths(workspace_root)
  output <- f6gp_paths(paths, create = TRUE)
  f6r_require_files(c(output$curve, output$source_validation),
                    "day-1000 positive-growth plotting data")
  validation <- f6r_read_tsv(output$source_validation)
  if (nrow(validation) != 4L || !all(validation$passed)) {
    stop("Day-1000 Figure 6B replay validation did not pass.")
  }
  curves <- f6r_read_tsv(output$curve)
  filename <- f6gp_filename(family)
  result <- list(
    pdf = file.path(paths$figures, paste0(filename, ".pdf")),
    png = file.path(paths$figures, paste0(filename, ".png"))
  )
  dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)
  page_width <- 14.6 * length(f6gp_p_values()) / 3
  grDevices::cairo_pdf(result$pdf, width = page_width, height = 9.5,
                       family = "sans", onefile = TRUE)
  tryCatch(f6gp_draw_page(curves, family),
           finally = grDevices::dev.off())
  grDevices::png(result$png, width = page_width, height = 9.5,
                 units = "in", res = 220, type = "cairo")
  tryCatch(f6gp_draw_page(curves, family),
           finally = grDevices::dev.off())
  if (!all(file.exists(unlist(result)))) {
    stop("Growth-permissive PDF/PNG export is incomplete.")
  }
  published <- c(
    pdf = f6r_publish(result$pdf, file.path(
      paths$manuscript_figures, basename(result$pdf))),
    png = f6r_publish(result$png, file.path(
      paths$manuscript_figures, basename(result$png)))
  )
  if (!all(vapply(names(result), function(format) {
    identical(f6r_md5(result[[format]]), f6r_md5(published[[format]]))
  }, logical(1L)))) {
    stop("Growth-permissive publication copies do not match.")
  }
  invisible(result)
}
