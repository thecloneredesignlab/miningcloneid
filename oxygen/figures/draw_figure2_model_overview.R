#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(grid))

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})

panel_fill <- "#F7F8FA"
panel_border <- "#C9CED6"
ink <- "#20242A"
muted <- "#616873"
blue <- "#2C7BB6"
blue_light <- "#DCEEF7"
red <- "#C43C39"
red_light <- "#F8E0DE"
amber <- "#C98A00"
amber_light <- "#FFF0C7"
green <- "#2E7D5B"
green_light <- "#DCEFE7"

gp_text <- function(size = 9, col = ink, fontface = "plain") {
  gpar(fontsize = size, col = col, fontface = fontface, fontfamily = "Helvetica")
}

draw_box <- function(x, y, w, h, label, fill = "white", border = panel_border,
                     size = 8.5, fontface = "plain", r = unit(2, "mm")) {
  grid.roundrect(x, y, w, h, r = r, gp = gpar(fill = fill, col = border, lwd = 1))
  grid.text(label, x, y, gp = gp_text(size, ink, fontface), just = "centre")
}

draw_arrow <- function(x0, y0, x1, y1, col = muted, lwd = 1.4, type = "closed") {
  grid.lines(
    x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
    arrow = arrow(length = unit(2.2, "mm"), type = type),
    gp = gpar(col = col, lwd = lwd, lineend = "round")
  )
}

draw_panel <- function(x, y, w, h, label, title) {
  grid.roundrect(x, y, w, h, r = unit(2.2, "mm"), gp = gpar(fill = panel_fill, col = panel_border, lwd = 1.1))
  grid.text(label, x - w / 2 + 0.015, y + h / 2 - 0.018, just = c("left", "top"), gp = gp_text(12, ink, "bold"))
  grid.text(title, x - w / 2 + 0.045, y + h / 2 - 0.019, just = c("left", "top"), gp = gp_text(9.5, ink, "bold"))
}

draw_mini_cost_plot <- function(x0, y0, w, h) {
  grid.lines(unit(c(x0, x0, x0 + w), "npc"), unit(c(y0 + h, y0, y0), "npc"), gp = gpar(col = muted, lwd = 0.8))
  grid.text("Stress", x0 + w / 2, y0 - 0.018, gp = gp_text(7, muted))
  grid.text("Net expansion", x0 - 0.018, y0 + h / 2, rot = 90, gp = gp_text(7, muted))
  grid.lines(unit(c(x0 + 0.01, x0 + w - 0.01), "npc"), unit(c(y0 + h * 0.78, y0 + h * 0.48), "npc"), gp = gpar(col = blue, lwd = 2))
  grid.lines(unit(c(x0 + 0.01, x0 + w - 0.01), "npc"), unit(c(y0 + h * 0.72, y0 + h * 0.12), "npc"), gp = gpar(col = red, lwd = 2))
  grid.text("2N", x0 + w - 0.005, y0 + h * 0.51, just = "right", gp = gp_text(7, blue, "bold"))
  grid.text("4N", x0 + w - 0.005, y0 + h * 0.15, just = "right", gp = gp_text(7, red, "bold"))
}

draw_chromosomes <- function(x, y, n = 4L, col = amber, scale = 1) {
  offsets <- seq(-0.018, 0.018, length.out = n)
  for (i in seq_along(offsets)) {
    grid.lines(unit(c(x + offsets[[i]], x + offsets[[i]]), "npc"), unit(c(y - 0.016 * scale, y + 0.016 * scale), "npc"), gp = gpar(col = col, lwd = 2.3, lineend = "round"))
  }
}

draw_figure2 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  draw_panel(0.5, 0.865, 0.94, 0.22, "A", "Resource state links competing ploidy pressures")
  draw_box(0.13, 0.865, 0.16, 0.075, "Low O2\nresource limitation", red_light, red, 8.2, "bold")
  draw_box(0.50, 0.865, 0.18, 0.075, "Resource state\nO2 (%)", blue_light, blue, 8.5, "bold")
  draw_box(0.87, 0.865, 0.16, 0.075, "Higher O2\nresource permissive", blue_light, blue, 8.2, "bold")
  draw_arrow(0.22, 0.865, 0.405, 0.865, red, 1.5)
  draw_arrow(0.595, 0.865, 0.78, 0.865, blue, 1.5)
  grid.text("stress", 0.31, 0.892, gp = gp_text(7.5, red, "bold"))
  grid.text("growth capacity", 0.69, 0.892, gp = gp_text(7.5, blue, "bold"))
  grid.text("The same resource axis changes growth/death, CIN generation, and post-event survival.", 0.5, 0.795, gp = gp_text(8.2, muted))

  draw_panel(0.18, 0.535, 0.30, 0.36, "B", "Cost under stress")
  draw_mini_cost_plot(0.085, 0.405, 0.19, 0.14)
  draw_box(0.18, 0.640, 0.20, 0.050, "growth suppression", red_light, red, 7.8, "bold")
  draw_box(0.18, 0.585, 0.20, 0.050, "hypoxic death", red_light, red, 7.8, "bold")

  draw_panel(0.50, 0.535, 0.30, 0.36, "C", "Instability under stress")
  draw_box(0.50, 0.640, 0.19, 0.050, "Low O2 / death stress", red_light, red, 7.8, "bold")
  draw_arrow(0.50, 0.612, 0.50, 0.573, red, 1.4)
  draw_box(0.43, 0.535, 0.115, 0.060, "p_mis(O2,N)\nincreases", amber_light, amber, 7.5, "bold")
  draw_box(0.57, 0.535, 0.115, 0.060, "WGD\nbranch", amber_light, amber, 7.5, "bold")
  draw_arrow(0.43, 0.502, 0.43, 0.458, amber, 1.3)
  draw_arrow(0.57, 0.502, 0.57, 0.458, amber, 1.3)
  draw_chromosomes(0.43, 0.425, 3, amber)
  draw_chromosomes(0.50, 0.425, 5, amber)
  draw_chromosomes(0.58, 0.425, 7, amber)
  grid.text("chromosome-number variation", 0.50, 0.397, gp = gp_text(7.8, muted))

  draw_panel(0.82, 0.535, 0.30, 0.36, "D", "Ploidy-dependent buffering")
  draw_box(0.75, 0.630, 0.105, 0.055, "2N mother", blue_light, blue, 7.8, "bold")
  draw_box(0.89, 0.630, 0.105, 0.055, "4N mother", green_light, green, 7.8, "bold")
  draw_arrow(0.75, 0.600, 0.75, 0.550, blue, 1.3)
  draw_arrow(0.89, 0.600, 0.89, 0.550, green, 1.3)
  draw_box(0.75, 0.515, 0.115, 0.060, "N +/- m\ndaughters", "white", blue, 7.5)
  draw_box(0.89, 0.515, 0.115, 0.060, "N +/- m\ndaughters", "white", green, 7.5)
  draw_arrow(0.75, 0.482, 0.75, 0.435, red, 1.3)
  draw_arrow(0.89, 0.482, 0.89, 0.435, green, 1.3)
  grid.text("lower survival", 0.75, 0.410, gp = gp_text(7.5, red, "bold"))
  grid.text("buffered survival", 0.89, 0.410, gp = gp_text(7.5, green, "bold"))
  grid.text("More copies can buffer otherwise lethal losses", 0.82, 0.375, gp = gp_text(7.5, muted))

  draw_arrow(0.50, 0.755, 0.18, 0.725, muted, 1.0)
  draw_arrow(0.50, 0.755, 0.50, 0.725, muted, 1.0)
  draw_arrow(0.50, 0.755, 0.82, 0.725, muted, 1.0)

  draw_panel(0.50, 0.175, 0.94, 0.25, "E", "Integrated long-term outcome")
  draw_box(0.13, 0.175, 0.16, 0.075, "growth and\ndeath rates", red_light, red, 7.8, "bold")
  draw_box(0.33, 0.175, 0.16, 0.075, "CIN / WGD\ntransition kernel", amber_light, amber, 7.8, "bold")
  draw_box(0.53, 0.175, 0.16, 0.075, "post-event\nsurvival filter", green_light, green, 7.8, "bold")
  draw_arrow(0.215, 0.175, 0.245, 0.175, muted, 1.2)
  draw_arrow(0.415, 0.175, 0.445, 0.175, muted, 1.2)
  draw_arrow(0.615, 0.175, 0.655, 0.175, muted, 1.2)
  draw_box(0.72, 0.175, 0.11, 0.075, "M(O2)", blue_light, blue, 9, "bold")
  draw_arrow(0.78, 0.175, 0.815, 0.175, blue, 1.3)
  draw_box(0.87, 0.175, 0.13, 0.075, "dominant\neigenvector v1", "white", blue, 7.8, "bold")
  grid.text(
    "Attractor: M(O2) v1 = lambda1 v1;  long-term mean ploidy = sum_N (N / 22) v1,N",
    0.50, 0.105, gp = gp_text(7.4, muted)
  )
}

png_path <- file.path(script_dir, "assembled_fig2.png")
pdf_path <- file.path(script_dir, "assembled_fig2.pdf")
grDevices::png(png_path, width = 2400, height = 2100, res = 300, bg = "white", type = "cairo")
draw_figure2()
grDevices::dev.off()
grDevices::cairo_pdf(pdf_path, width = 8, height = 7, family = "Helvetica")
draw_figure2()
grDevices::dev.off()
message("Wrote ", png_path)
message("Wrote ", pdf_path)
