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

blue <- "#2B6FA3"
blue_dark <- "#1D557E"
blue_light <- "#E7F1F8"

red <- "#D24B45"
red_dark <- "#A93632"
red_light <- "#FBEAE8"

amber <- "#D99000"
amber_dark <- "#A96E00"
amber_light <- "#FFF3D2"

green <- "#278261"
green_dark <- "#1D654A"
green_light <- "#E5F3EC"

purple <- "#6E5AA7"
purple_dark <- "#51417E"
purple_light <- "#EFECF8"

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
    gp = gpar(col = red, lwd = 2.0, lineend = "round")
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
    copies = 3L,
    fill = daughter_fill,
    border = daughter_border,
    chromosome_col = daughter_border,
    alpha = daughter_alpha,
    lwd = 1.1
  )
  draw_cell(
    x_right, y - 0.030,
    r = 0.020,
    copies = 5L,
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
    gp = gp_text(6.6, if (daughter_alpha < 0.8) dead else green_dark, "bold")
  )
}

# -----------------------------------------------------------------------------
# Panel A: opposing pressures and closed feedback
# -----------------------------------------------------------------------------

draw_panel_a <- function() {
  draw_panel(
    0.50, 0.735, 0.96, 0.49,
    "A",
    "Resource limitation couples opposing pressures on ploidy evolution",
    fill = "#FCFDFE"
  )

  draw_pill(
    0.50, 0.887, 0.30, 0.058,
    "LOW O2: resource stress h(O2) rises",
    fill = red_light,
    border = red,
    col = red_dark,
    size = 8.7
  )

  draw_card(
    0.27, 0.787, 0.33, 0.105,
    "FITNESS SELECTION\nhigher N: division falls more; death rises",
    fill = blue_light,
    border = blue,
    size = 8.2,
    col = navy,
    fontface = "bold",
    lwd = 1.2
  )
  draw_card(
    0.73, 0.787, 0.33, 0.105,
    "STATE GENERATION + RETENTION\np_mis per chromosome rises via mu_eff\nN +/- m states; higher-N buffering",
    fill = amber_light,
    border = amber,
    size = 7.5,
    col = amber_dark,
    fontface = "bold",
    lwd = 1.2
  )

  draw_curve_arrow(0.445, 0.858, 0.31, 0.842, curvature = 0.12, col = red, lwd = 1.5)
  draw_curve_arrow(0.555, 0.858, 0.69, 0.842, curvature = -0.12, col = amber, lwd = 1.5)

  grid.text(
    "selects against chromosome-rich states",
    0.27, 0.716,
    gp = gp_text(7.2, red_dark, "bold")
  )
  draw_arrow(0.49, 0.684, 0.20, 0.684, col = red, lwd = 2.7, length_mm = 2.7)

  grid.text(
    "CIN explores neighboring states",
    0.70, 0.716,
    gp = gp_text(7.2, amber_dark, "bold")
  )
  draw_arrow(0.52, 0.684, 0.79, 0.684, col = amber, lwd = 1.9,
             ends = "both", length_mm = 2.4)

  grid.lines(
    x = unit(c(0.16, 0.84), "npc"),
    y = unit(c(0.635, 0.635), "npc"),
    gp = gpar(col = navy, lwd = 1.0)
  )
  for (tick in seq(0.22, 0.78, length.out = 5)) {
    grid.lines(
      x = unit(c(tick, tick), "npc"),
      y = unit(c(0.629, 0.641), "npc"),
      gp = gpar(col = navy, lwd = 0.8)
    )
  }

  draw_cell(0.15, 0.635, r = 0.031, copies = 3L,
            fill = blue_light, border = blue, chromosome_col = blue)
  draw_cell(0.85, 0.635, r = 0.036, copies = 7L,
            fill = green_light, border = green, chromosome_col = green)
  grid.text("lower N", 0.15, 0.587, gp = gp_text(7.2, blue_dark, "bold"))
  grid.text("higher N", 0.85, 0.587, gp = gp_text(7.2, green_dark, "bold"))

  grid.polygon(
    x = unit(c(0.50, 0.487, 0.50, 0.513), "npc"),
    y = unit(c(0.654, 0.635, 0.616, 0.635), "npc"),
    gp = gpar(fill = purple, col = purple_dark, lwd = 1.0)
  )
  draw_card(
    0.50, 0.555, 0.31, 0.066,
    "conditional dominant composition\nv_N(O2; parameters)",
    fill = purple_light,
    border = purple,
    size = 8.0,
    col = purple_dark,
    fontface = "bold",
    lwd = 1.2
  )
  draw_arrow(0.50, 0.615, 0.50, 0.591, col = purple, lwd = 1.35)

  draw_pill(
    0.73, 0.590, 0.20, 0.043,
    "per-division WGD: N -> 2N",
    fill = amber_light,
    border = amber,
    col = amber_dark,
    size = 7.1
  )
  draw_arrow(0.62, 0.610, 0.78, 0.610, col = amber, lwd = 1.5)

  draw_poly_arrow(
    x = c(0.655, 0.935, 0.935, 0.655),
    y = c(0.545, 0.545, 0.887, 0.887),
    col = blue,
    lwd = 1.25,
    lty = 2,
    length_mm = 2.2
  )
  grid.text(
    "ploidy-weighted viable demand\nfeeds back to effective O2",
    0.958, 0.716,
    just = "centre",
    rot = 90,
    gp = gp_text(6.4, blue_dark, "bold")
  )
}

# -----------------------------------------------------------------------------
# Panel B: state-specific fitness
# -----------------------------------------------------------------------------

draw_panel_b <- function() {
  draw_panel(
    0.17, 0.255, 0.285, 0.43,
    "B",
    "Low-O2 fitness selection",
    title_size = 9.1
  )

  grid.lines(
    x = unit(c(0.078, 0.096), "npc"),
    y = unit(c(0.402, 0.402), "npc"),
    gp = gpar(col = blue, lwd = 2.0)
  )
  grid.text("lower N", 0.101, 0.402, just = "left", gp = gp_text(6.8, blue_dark, "bold"))
  grid.lines(
    x = unit(c(0.183, 0.201), "npc"),
    y = unit(c(0.402, 0.402), "npc"),
    gp = gpar(col = red, lwd = 2.0)
  )
  grid.text("higher N", 0.206, 0.402, just = "left", gp = gp_text(6.8, red_dark, "bold"))

  draw_rate_plot(0.078, 0.285, 0.184, 0.080, "Effective division rate", "division")
  draw_rate_plot(0.078, 0.155, 0.184, 0.080, "Hypoxia-associated death", "death")

  draw_arrow(0.086, 0.128, 0.246, 0.128, col = blue, lwd = 1.0)
  grid.text("O2 increasing", 0.166, 0.108, gp = gp_text(6.8, muted, "bold"))

  draw_card(
    0.17, 0.073, 0.225, 0.046,
    "Low O2 favors lower-N expansion",
    fill = red_light,
    border = red,
    size = 7.2,
    col = red_dark,
    fontface = "bold",
    lwd = 1.0
  )
}

# -----------------------------------------------------------------------------
# Panel C: chromosome-state generation and survival filtering
# -----------------------------------------------------------------------------

draw_panel_c <- function() {
  draw_panel(
    0.515, 0.255, 0.375, 0.43,
    "C",
    "State generation + survival filtering",
    title_size = 8.8
  )

  draw_card(
    0.515, 0.392, 0.205, 0.052,
    "p_mis per chromosome rises\nvia mu_eff",
    fill = amber_light,
    border = amber,
    col = amber_dark,
    size = 6.9,
    fontface = "bold",
    lwd = 1.0,
    radius_mm = 1.8
  )

  grid.text("N - m", 0.558, 0.350, gp = gp_text(6.9, amber_dark, "bold"))
  grid.text("N + m", 0.658, 0.350, gp = gp_text(6.9, amber_dark, "bold"))

  grid.text("lower-N\nmother", 0.350, 0.300, just = "left", gp = gp_text(7.0, blue_dark, "bold"))
  draw_cell(0.430, 0.300, r = 0.026, copies = 4L,
            fill = blue_light, border = blue, chromosome_col = blue)
  draw_branch(
    0.430, 0.300, 0.558, 0.658,
    daughter_alpha = 0.55,
    daughter_fill = dead_light,
    daughter_border = dead,
    label = "both shifts: lower survival",
    label_x = 0.455,
    label_y = 0.238
  )

  grid.text("higher-N\nmother", 0.350, 0.190, just = "left", gp = gp_text(7.0, green_dark, "bold"))
  draw_cell(0.430, 0.190, r = 0.028, copies = 7L,
            fill = green_light, border = green, chromosome_col = green)
  draw_branch(
    0.430, 0.190, 0.558, 0.658,
    daughter_alpha = 1.0,
    daughter_fill = green_light,
    daughter_border = green,
    label = "both shifts: buffered survival",
    label_x = 0.455,
    label_y = 0.118
  )

  draw_card(
    0.414, 0.069, 0.155, 0.050,
    "WGD: N -> 2N\nconstant per division",
    fill = amber_light,
    border = amber,
    col = amber_dark,
    size = 6.2,
    fontface = "bold",
    lwd = 1.0,
    radius_mm = 1.5
  )
  draw_card(
    0.604, 0.069, 0.172, 0.050,
    "out-of-grid offspring\nenter nonviable pool",
    fill = dead_light,
    border = dead,
    col = dead,
    size = 5.9,
    fontface = "bold",
    lwd = 1.0,
    radius_mm = 1.5
  )
}

# -----------------------------------------------------------------------------
# Panel D: fixed-O2 population operator and conditional dominant mode
# -----------------------------------------------------------------------------

draw_panel_d <- function() {
  draw_panel(
    0.85, 0.255, 0.28, 0.43,
    "D",
    "Population integration"
  )

  draw_card(
    0.85, 0.382, 0.225, 0.055,
    "fitness + CIN/WGD + survival",
    fill = white,
    border = hairline,
    size = 7.2,
    col = ink,
    fontface = "bold"
  )
  draw_arrow(0.85, 0.352, 0.85, 0.324, col = muted, lwd = 1.2)

  draw_card(
    0.85, 0.292, 0.180, 0.058,
    "fixed-O2 operator\nM(O2)",
    fill = blue_light,
    border = blue,
    size = 7.7,
    col = navy,
    fontface = "bold",
    lwd = 1.2
  )
  draw_arrow(0.85, 0.260, 0.85, 0.226, col = purple, lwd = 1.35)

  grid.text(
    "conditional dominant composition v_N",
    0.85, 0.214,
    gp = gp_text(7.2, purple_dark, "bold")
  )
  draw_distribution(0.765, 0.125, 0.170, 0.067)

  grid.text(
    "mean ploidy = sum(N v_N) / 22",
    0.85, 0.087,
    gp = gp_text(6.7, purple_dark, "bold")
  )
  draw_card(
    0.85, 0.061, 0.230, 0.032,
    "interpret when the spectral gap resolves v_N",
    fill = purple_light,
    border = purple,
    col = purple_dark,
    size = 6.1,
    fontface = "bold",
    lwd = 0.9,
    radius_mm = 1.2
  )
}

draw_figure2 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = white, col = NA))
  draw_panel_a()
  draw_panel_b()
  draw_panel_c()
  draw_panel_d()
}

png_path <- file.path(output_dir, "assembled_fig2.png")
pdf_path <- file.path(output_dir, "assembled_fig2.pdf")

grDevices::png(
  png_path,
  width = 2400,
  height = 2100,
  res = 300,
  bg = white,
  type = "cairo"
)
draw_figure2()
grDevices::dev.off()

grDevices::cairo_pdf(
  pdf_path,
  width = 8,
  height = 7,
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
