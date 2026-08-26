#!/usr/bin/env Rscript

if (identical(Sys.getenv("FIGURE2_DRAW_WORKER"), "1")) {

suppressPackageStartupMessages(library(grid))

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    getwd()
  }
})

workspace_root <- normalizePath(
  Sys.getenv(
    "FIGURE_WORKSPACE_ROOT",
    unset = file.path(script_dir, "..", "..")
  ),
  mustWork = TRUE
)
default_output_dir <- file.path(workspace_root, "Figures")
output_dir <- Sys.getenv("FIGURE2_OUTPUT_DIR", unset = default_output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Visual system
# -----------------------------------------------------------------------------

ink <- "#19232D"
navy <- "#173B5C"
muted <- "#5F6B78"
hairline <- "#D7DEE7"
panel_fill <- "#FBFCFE"
white <- "#FFFFFF"

blue <- "#0072B2"
blue_dark <- "#005A8D"
blue_light <- "#E7F1F8"

red <- "#D55E00"
red_dark <- "#A84300"
red_light <- "#FCEADF"

amber <- "#E69F00"
amber_dark <- "#A96E00"
amber_light <- "#FFF3D2"

green <- "#009E73"
green_dark <- "#007A59"
green_light <- "#E5F3EC"

purple <- "#6E5AA7"
purple_dark <- "#51417E"
purple_light <- "#EFECF8"

magenta <- "#B14A86"
magenta_dark <- "#873766"
magenta_light <- "#F6E8F1"

dead <- "#8B4A4A"
dead_light <- "#F3E8E7"

gp_text <- function(size = 9, col = ink, fontface = "plain") {
  gpar(
    fontsize = size,
    col = col,
    fontface = fontface,
    fontfamily = "Helvetica",
    lineheight = 1.05
  )
}

draw_panel <- function(x, y, w, h, label, title, fill = panel_fill,
                       title_size = 10.1) {
  grid.roundrect(
    x, y, w, h,
    r = unit(2.1, "mm"),
    gp = gpar(fill = fill, col = hairline, lwd = 0.9)
  )
  grid.text(
    label,
    x - w / 2 + 0.016,
    y + h / 2 - 0.017,
    just = c("left", "top"),
    gp = gp_text(13.5, ink, "bold")
  )
  grid.text(
    title,
    x - w / 2 + 0.052,
    y + h / 2 - 0.019,
    just = c("left", "top"),
    gp = gp_text(title_size, ink, "bold")
  )
}

draw_card <- function(x, y, w, h, label, fill = white, border = hairline,
                      size = 8.2, col = ink, fontface = "plain",
                      lwd = 1.0, radius_mm = 1.7) {
  grid.roundrect(
    x, y, w, h,
    r = unit(radius_mm, "mm"),
    gp = gpar(fill = fill, col = border, lwd = lwd)
  )
  grid.text(
    label, x, y,
    gp = gp_text(size, col, fontface),
    just = "centre"
  )
}

draw_labeled_card <- function(x, y, w, h, title, subtitle,
                              fill = white, border = hairline,
                              title_size = 8.2, subtitle_size = 6.7,
                              title_col = ink, subtitle_col = muted,
                              lwd = 1.15) {
  grid.roundrect(
    x, y, w, h,
    r = unit(1.7, "mm"),
    gp = gpar(fill = fill, col = border, lwd = lwd)
  )
  grid.text(
    title, x, y + h * 0.15,
    gp = gp_text(title_size, title_col, "bold")
  )
  grid.text(
    subtitle, x, y - h * 0.23,
    gp = gp_text(subtitle_size, subtitle_col)
  )
}

draw_pill <- function(x, y, w, h, label, fill, border, col = ink,
                      size = 8.0, fontface = "bold") {
  draw_card(
    x, y, w, h, label,
    fill = fill,
    border = border,
    size = size,
    col = col,
    fontface = fontface,
    lwd = 1.0,
    radius_mm = 4
  )
}

draw_arrow <- function(x0, y0, x1, y1, col = muted, lwd = 1.35,
                       lty = 1, ends = "last", length_mm = 2.2) {
  arrow_spec <- if (identical(ends, "none")) {
    NULL
  } else {
    arrow(
      length = unit(length_mm, "mm"),
      type = "closed",
      ends = ends
    )
  }
  grid.lines(
    x = unit(c(x0, x1), "npc"),
    y = unit(c(y0, y1), "npc"),
    arrow = arrow_spec,
    gp = gpar(
      col = col,
      fill = col,
      lwd = lwd,
      lty = lty,
      lineend = "round",
      linejoin = "round"
    )
  )
}

draw_poly_arrow <- function(x, y, col = muted, lwd = 1.2, lty = 1,
                            length_mm = 2.1) {
  grid.lines(
    x = unit(x, "npc"),
    y = unit(y, "npc"),
    arrow = arrow(length = unit(length_mm, "mm"), type = "closed"),
    gp = gpar(
      col = col,
      fill = col,
      lwd = lwd,
      lty = lty,
      lineend = "round",
      linejoin = "round"
    )
  )
}

draw_curve_arrow <- function(x0, y0, x1, y1, curvature = 0.3,
                             col = muted, lwd = 1.2, lty = 1,
                             length_mm = 2.1) {
  grid.curve(
    x1 = unit(x0, "npc"),
    y1 = unit(y0, "npc"),
    x2 = unit(x1, "npc"),
    y2 = unit(y1, "npc"),
    curvature = curvature,
    square = FALSE,
    ncp = 8,
    arrow = arrow(length = unit(length_mm, "mm"), type = "closed"),
    gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round")
  )
}

draw_cell <- function(x, y, r = 0.026, copies = 4L,
                      fill = white, border = blue,
                      chromosome_col = border, alpha = 1,
                      lwd = 1.3) {
  fill_use <- adjustcolor(fill, alpha.f = alpha)
  border_use <- adjustcolor(border, alpha.f = alpha)
  chromosome_use <- adjustcolor(chromosome_col, alpha.f = alpha)
  grid.circle(
    x, y, r,
    gp = gpar(fill = fill_use, col = border_use, lwd = lwd)
  )

  copies <- max(2L, as.integer(copies))
  n_cols <- ceiling(copies / 2)
  xs <- seq(x - r * 0.50, x + r * 0.50, length.out = n_cols)
  bars <- list()
  k <- 1L
  for (cx in xs) {
    for (row in seq_len(2L)) {
      if (k > copies) break
      cy <- y + if (row == 1L) r * 0.22 else -r * 0.22
      bars[[k]] <- c(cx, cy)
      k <- k + 1L
    }
  }
  for (bar in bars) {
    grid.lines(
      x = unit(c(bar[[1]], bar[[1]]), "npc"),
      y = unit(c(bar[[2]] - r * 0.16, bar[[2]] + r * 0.16), "npc"),
      gp = gpar(col = chromosome_use, lwd = 1.6, lineend = "round")
    )
  }
}

draw_rate_plot <- function(x0, y0, w, h, title, mode = c("division", "death")) {
  mode <- match.arg(mode)

  grid.text(
    title,
    x0,
    y0 + h + 0.018,
    just = c("left", "bottom"),
    gp = gp_text(7.7, ink, "bold")
  )
  grid.lines(
    x = unit(c(x0, x0, x0 + w), "npc"),
    y = unit(c(y0 + h, y0, y0), "npc"),
    gp = gpar(col = muted, lwd = 0.7)
  )

  x_frac <- c(0.04, 0.28, 0.56, 0.96)
  x_vals <- x0 + w * x_frac
  if (mode == "division") {
    lower_y <- y0 + h * c(0.46, 0.62, 0.79, 0.91)
    higher_y <- y0 + h * c(0.12, 0.34, 0.67, 0.88)
  } else {
    lower_y <- y0 + h * c(0.47, 0.31, 0.15, 0.05)
    higher_y <- y0 + h * c(0.90, 0.64, 0.26, 0.06)
  }

  grid.lines(
    x = unit(x_vals, "npc"),
    y = unit(lower_y, "npc"),
    gp = gpar(col = blue, lwd = 2.0, lineend = "round")
  )
  grid.lines(
    x = unit(x_vals, "npc"),
    y = unit(higher_y, "npc"),
    gp = gpar(col = red, lwd = 2.0, lty = 2, lineend = "round")
  )
  grid.circle(
    x_vals[[1]], lower_y[[1]], unit(0.9, "mm"),
    gp = gpar(fill = blue, col = blue)
  )
  grid.circle(
    x_vals[[1]], higher_y[[1]], unit(0.9, "mm"),
    gp = gpar(fill = red, col = red)
  )
}

draw_distribution <- function(x0, y0, w, h) {
  grid.lines(
    x = unit(c(x0, x0, x0 + w), "npc"),
    y = unit(c(y0 + h, y0, y0), "npc"),
    gp = gpar(col = muted, lwd = 0.65)
  )
  xx <- seq(0, 1, length.out = 60)
  yy <- 0.18 * exp(-0.5 * ((xx - 0.26) / 0.13)^2) +
    0.92 * exp(-0.5 * ((xx - 0.62) / 0.17)^2)
  yy <- yy / max(yy)
  grid.lines(
    x = unit(x0 + w * xx, "npc"),
    y = unit(y0 + h * yy, "npc"),
    gp = gpar(col = purple, lwd = 2.0, lineend = "round")
  )
  grid.text(
    "chromosome state N",
    x0 + w / 2,
    y0 - 0.014,
    gp = gp_text(6.7, muted)
  )
}

draw_branch <- function(x_mother, y, x_left, x_right, daughter_alpha,
                        daughter_fill, daughter_border, label,
                        daughter_copies = c(3L, 5L),
                        label_x = (x_left + x_right) / 2,
                        label_y = y - 0.070) {
  branch_x <- x_mother + 0.055
  draw_arrow(x_mother + 0.025, y, branch_x, y, col = amber, lwd = 1.2, ends = "none")
  grid.lines(
    x = unit(c(branch_x, branch_x, branch_x), "npc"),
    y = unit(c(y, y + 0.030, y - 0.030), "npc"),
    gp = gpar(col = amber, lwd = 1.2, lineend = "round")
  )
  draw_arrow(branch_x, y + 0.030, x_left - 0.022, y + 0.030,
             col = amber, lwd = 1.2)
  draw_arrow(branch_x, y - 0.030, x_right - 0.022, y - 0.030,
             col = amber, lwd = 1.2)
  draw_cell(
    x_left, y + 0.030,
    r = 0.020,
    copies = daughter_copies[[1L]],
    fill = daughter_fill,
    border = daughter_border,
    chromosome_col = daughter_border,
    alpha = daughter_alpha,
    lwd = 1.1
  )
  draw_cell(
    x_right, y - 0.030,
    r = 0.020,
    copies = daughter_copies[[2L]],
    fill = daughter_fill,
    border = daughter_border,
    chromosome_col = daughter_border,
    alpha = daughter_alpha,
    lwd = 1.1
  )
  grid.text(
    label,
    label_x,
    label_y,
    gp = gp_text(7.0, if (daughter_alpha < 0.8) dead else green_dark, "bold")
  )
}

# -----------------------------------------------------------------------------
# Panel A: death-linked CIN and adaptive feedback
# -----------------------------------------------------------------------------

draw_panel_a <- function() {
  draw_panel(
    0.50, 0.735, 0.96, 0.50,
    "A",
    "Death-linked CIN creates an adaptive feedback",
    fill = "#FCFDFE",
    title_size = 9.8
  )

  grid.text("Environmental input", 0.50, 0.925, gp = gp_text(6.8, muted))
  draw_labeled_card(
    0.50, 0.875, 0.25, 0.075,
    "Resource limitation", "modeled through O₂",
    fill = blue_light, border = blue,
    title_size = 8.6, subtitle_size = 7.0,
    title_col = navy, subtitle_col = blue_dark
  )

  draw_card(
    0.155, 0.745, 0.23, 0.065,
    "Proliferation rate decreases",
    fill = "#F4F5F6", border = muted,
    size = 7.8, col = ink, fontface = "bold"
  )
  draw_labeled_card(
    0.465, 0.745, 0.285, 0.090,
    "State-specific death hazard increases", "experienced cellular stress",
    fill = amber_light, border = amber,
    title_size = 7.7, subtitle_size = 7.0,
    title_col = amber_dark
  )
  draw_labeled_card(
    0.805, 0.745, 0.265, 0.090,
    "Chromosome missegregation increases", "effective per-chromosome probability",
    fill = magenta_light, border = magenta,
    title_size = 7.4, subtitle_size = 7.0,
    title_col = magenta_dark
  )

  draw_arrow(0.43, 0.842, 0.25, 0.780, col = blue, lwd = 1.45)
  draw_arrow(0.50, 0.837, 0.475, 0.794, col = amber, lwd = 1.45)
  draw_arrow(0.610, 0.745, 0.670, 0.745, col = magenta, lwd = 1.6)
  grid.text(
    "death-linked\nmutagenesis",
    0.640, 0.812,
    gp = gp_text(7.0, magenta_dark, "bold")
  )

  draw_card(
    0.805, 0.617, 0.245, 0.063,
    "Chromosome-number variation",
    fill = green_light, border = green,
    size = 7.7, col = green_dark, fontface = "bold"
  )
  draw_labeled_card(
    0.525, 0.617, 0.220, 0.075,
    "Post-MS survival", "+ state-dependent fitness selection",
    fill = "#F4F5F6", border = muted,
    title_size = 7.7, subtitle_size = 7.0
  )
  draw_labeled_card(
    0.285, 0.526, 0.280, 0.065,
    "Adapted karyotype /", "ploidy composition",
    fill = green_light, border = green,
    title_size = 7.8, subtitle_size = 7.0,
    title_col = green_dark, subtitle_col = green_dark
  )
  draw_labeled_card(
    0.805, 0.526, 0.245, 0.065,
    "WGD generation", "constant probability per division",
    fill = "#F4F5F6", border = muted,
    title_size = 7.7, subtitle_size = 7.0
  )

  draw_arrow(0.805, 0.699, 0.805, 0.652, col = magenta, lwd = 1.5)
  draw_arrow(0.680, 0.617, 0.640, 0.617, col = green, lwd = 1.45)
  draw_arrow(0.465, 0.699, 0.500, 0.657, col = amber, lwd = 1.35)
  draw_curve_arrow(
    0.270, 0.725, 0.410, 0.640,
    curvature = -0.18, col = muted, lwd = 1.25, length_mm = 2.1
  )
  draw_arrow(0.475, 0.580, 0.420, 0.559, col = green, lwd = 1.45)
  draw_arrow(0.805, 0.560, 0.805, 0.582, col = muted, lwd = 1.25)

  draw_poly_arrow(
    x = c(0.285, 0.305, 0.305, 0.320),
    y = c(0.560, 0.595, 0.745, 0.745),
    col = green, lwd = 1.65, length_mm = 2.4
  )
  grid.text(
    "better-adapted states lower\npopulation-average death\n(even at fixed O₂)",
    0.205, 0.630,
    gp = gp_text(7.0, green_dark, "bold")
  )
}

# -----------------------------------------------------------------------------
# Panel B: state-specific direct fitness effects
# -----------------------------------------------------------------------------

draw_panel_b <- function() {
  draw_panel(
    0.205, 0.245, 0.36, 0.45,
    "B",
    "Direct fitness effects of\nchromosome burden",
    title_size = 8.5
  )

  grid.lines(
    x = unit(c(0.070, 0.090), "npc"),
    y = unit(c(0.397, 0.397), "npc"),
    gp = gpar(col = blue, lwd = 2.0)
  )
  grid.text("lower N", 0.096, 0.397, just = "left", gp = gp_text(7.0, blue_dark, "bold"))
  grid.lines(
    x = unit(c(0.202, 0.222), "npc"),
    y = unit(c(0.397, 0.397), "npc"),
    gp = gpar(col = red, lwd = 2.0, lty = 2)
  )
  grid.text("higher N", 0.228, 0.397, just = "left", gp = gp_text(7.0, red_dark, "bold"))

  draw_rate_plot(0.070, 0.285, 0.270, 0.070, "Effective proliferation", "division")
  draw_rate_plot(0.070, 0.160, 0.270, 0.070, "Stress-associated death", "death")

  draw_arrow(0.080, 0.137, 0.330, 0.137, col = blue, lwd = 1.0)
  grid.text("O₂ increasing", 0.205, 0.117, gp = gp_text(7.0, muted, "bold"))

  draw_card(
    0.205, 0.064, 0.270, 0.052,
    "Encoded direct cost\nis stronger at higher N",
    fill = red_light, border = red,
    size = 7.0, col = red_dark, fontface = "bold", lwd = 1.0
  )
}

# -----------------------------------------------------------------------------
# Panel C: correlated daughters and ploidy-dependent survival
# -----------------------------------------------------------------------------

draw_panel_c <- function() {
  draw_panel(
    0.690, 0.245, 0.585, 0.45,
    "C",
    "Ploidy-dependent survival after\nchromosome missegregation",
    title_size = 8.5
  )

  grid.text(
    "One CIN event; complementary daughter shifts with the same |m|",
    0.690, 0.382,
    gp = gp_text(7.0, muted, "bold")
  )
  grid.text("N - m", 0.740, 0.363, gp = gp_text(7.0, magenta_dark, "bold"))
  grid.text("N + m", 0.890, 0.363, gp = gp_text(7.0, magenta_dark, "bold"))

  grid.text("lower-N\nmother", 0.420, 0.294, just = "left", gp = gp_text(7.1, blue_dark, "bold"))
  draw_cell(0.540, 0.294, r = 0.026, copies = 4L,
            fill = blue_light, border = blue, chromosome_col = blue)
  draw_branch(
    0.540, 0.294, 0.740, 0.890,
    daughter_alpha = 0.50,
    daughter_fill = dead_light,
    daughter_border = dead,
    daughter_copies = c(3L, 5L),
    label = "lower survival after the same |m|",
    label_x = 0.815,
    label_y = 0.228
  )

  grid.text("higher-N\nmother", 0.420, 0.145, just = "left", gp = gp_text(7.1, green_dark, "bold"))
  draw_cell(0.540, 0.145, r = 0.029, copies = 7L,
            fill = green_light, border = green, chromosome_col = green)
  draw_branch(
    0.540, 0.145, 0.740, 0.890,
    daughter_alpha = 1.0,
    daughter_fill = green_light,
    daughter_border = green,
    daughter_copies = c(6L, 8L),
    label = "higher survival after the same |m|",
    label_x = 0.815,
    label_y = 0.067
  )
}

draw_figure2 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = white, col = NA))
  draw_panel_a()
  draw_panel_b()
  draw_panel_c()
}

png_path <- file.path(output_dir, "assembled_fig2.png")
pdf_path <- file.path(output_dir, "assembled_fig2.pdf")

grDevices::png(
  png_path,
  width = 2130,
  height = 1905,
  res = 300,
  bg = white,
  type = "cairo"
)
draw_figure2()
grDevices::dev.off()

grDevices::cairo_pdf(
  pdf_path,
  width = 7.1,
  height = 6.35,
  family = "Helvetica",
  bg = white
)
draw_figure2()
grDevices::dev.off()

message("Wrote ", png_path)
message("Wrote ", pdf_path)

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Figure2 <- function() {
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_Figure2.R"), mustWork = TRUE),
    env = c(
      "FIGURE2_DRAW_WORKER=1",
      paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
      paste0("FIGURE2_OUTPUT_DIR=", OUTPUT_ROOT)
    )
  )
  require_files(
    file.path(OUTPUT_ROOT, c("assembled_fig2.png", "assembled_fig2.pdf")),
    "Figure 2 output"
  )
  invisible(file.path(OUTPUT_ROOT, "assembled_fig2.png"))
}

if (sys.nframe() == 0L) draw_Figure2()

}
