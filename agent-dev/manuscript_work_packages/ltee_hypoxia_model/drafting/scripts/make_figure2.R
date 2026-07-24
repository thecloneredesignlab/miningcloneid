#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = getwd())
draft_root <- file.path(
  repo_root, "agent-dev", "manuscript_work_packages",
  "ltee_hypoxia_model", "drafting"
)

dark <- "#333333"
blue <- "#0072B2"
verm <- "#D55E00"
vitro_col <- "#56B4E9"
vivo_col <- "#009E73"
o2_dark <- "#084594"
o2_light <- "#8FC5E3"

type_cols <- c(
  observed = "#F3F3F3",
  state = "#DCEEF8",
  "function" = "#EFE3F2",
  transition = "#FBE9C8",
  output = "#DDF0E8"
)

box <- function(xmin, xmax, ymin, ymax, fill, colour = dark, radius = 0.08) {
  annotate(
    "rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    fill = fill, colour = colour, linewidth = 0.55
  )
}

label_text <- function(x, y, text, size = 3, fontface = "plain",
                       hjust = 0.5, vjust = 0.5, colour = dark,
                       lineheight = 1.05) {
  annotate(
    "text", x = x, y = y, label = text, size = size,
    family = "sans", fontface = fontface, hjust = hjust, vjust = vjust,
    colour = colour, lineheight = lineheight
  )
}

arrow_segment <- function(x, xend, y, yend, colour = dark,
                          linetype = "solid", linewidth = 0.65) {
  annotate(
    "segment", x = x, xend = xend, y = y, yend = yend,
    colour = colour, linetype = linetype, linewidth = linewidth,
    arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed")
  )
}

p <- ggplot() +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 10.2), clip = "off") +
  theme_void(base_family = "sans") +
  theme(
    plot.margin = margin(8, 9, 7, 9),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# A: context-specific resource input and modeled resource state
p <- p +
  box(0.30, 4.50, 4.65, 9.70, type_cols[["state"]], colour = o2_dark) +
  label_text(0.52, 9.42, "A", 3.8, "bold", hjust = 0) +
  label_text(0.92, 9.42, "Context & resource state", 3.05, "bold", hjust = 0) +
  box(0.65, 2.37, 8.02, 8.90, "white", colour = vitro_col) +
  label_text(1.51, 8.67, "external input", 2.50, "bold", colour = "#666666") +
  label_text(1.51, 8.30, "In vitro\nO\u2082 schedule", 2.50, "bold", colour = vitro_col) +
  box(2.57, 4.15, 8.02, 8.90, "white", colour = dark) +
  label_text(3.36, 8.67, "initial state", 2.50, "bold", colour = "#666666") +
  label_text(3.36, 8.30, "starting N\ndistribution", 2.50) +
  box(0.65, 4.15, 6.30, 7.72, "white", colour = o2_dark) +
  label_text(0.83, 7.45, "modeled state", 2.50, "bold", hjust = 0, colour = "#666666") +
  label_text(2.40, 7.10, "O\u2082,eff(t)", 3.55, "bold", colour = o2_dark) +
  label_text(
    2.40, 6.48,
    "culture: prescribed\ntumor: latent",
    2.50
  ) +
  box(0.65, 4.15, 5.05, 6.03, "white", colour = vivo_col) +
  label_text(2.40, 5.77, "In vivo resource path", 2.55, "bold", colour = vivo_col) +
  label_text(
    2.40, 5.36,
    "fitted supply + simulated\nlive-cell demand",
    2.50
  ) +
  label_text(
    2.40, 4.83,
    "lower O\u2082  \u2192  higher stress",
    2.55, "bold", colour = o2_dark
  ) +
  arrow_segment(1.50, 1.50, 8.02, 7.72, colour = vitro_col, linewidth = 0.55) +
  arrow_segment(3.35, 3.35, 8.02, 7.72, colour = dark, linewidth = 0.55) +
  arrow_segment(2.40, 2.40, 6.30, 6.03, colour = vivo_col, linetype = "dashed", linewidth = 0.55)

# B: chromosome-number-dependent growth and death
p <- p +
  box(5.30, 8.45, 7.08, 9.70, type_cols[["function"]]) +
  label_text(5.52, 9.42, "B", 3.8, "bold", hjust = 0) +
  label_text(5.92, 9.42, "Growth & death", 3.00, "bold", hjust = 0) +
  label_text(5.62, 8.93, "fitted functions", 2.50, "bold", hjust = 0, colour = "#666666") +
  label_text(6.88, 8.37, "\u03BB(N, O\u2082)     \u03BC(N, O\u2082)", 2.95, "bold") +
  label_text(
    6.88, 7.65,
    "N- and O\u2082-dependent\nproliferation + death",
    2.50
  ) +
  arrow_segment(4.50, 5.30, 7.32, 8.35, colour = o2_dark)

# C: stress-linked chromosome transitions
p <- p +
  box(5.30, 8.45, 3.88, 6.52, type_cols[["transition"]]) +
  label_text(5.52, 6.24, "C", 3.8, "bold", hjust = 0) +
  label_text(5.92, 6.24, "CIN & WGD", 3.00, "bold", hjust = 0) +
  label_text(5.62, 5.76, "state transitions", 2.50, "bold", hjust = 0, colour = "#666666") +
  label_text(
    6.10, 4.92,
    "Missegregation\nN \u2192 N \u00B1 \u0394N",
    2.50, "bold"
  ) +
  label_text(
    7.67, 4.92,
    "WGD\nN \u2192 2N",
    2.50, "bold"
  ) +
  label_text(
    6.88, 4.18,
    "\u0394N: multi-copy shift",
    2.50
  ) +
  arrow_segment(4.50, 5.30, 6.55, 5.25, colour = o2_dark)

# D: ploidy-dependent post-missegregation survival
p <- p +
  box(5.30, 8.45, 0.78, 3.32, type_cols[["function"]]) +
  label_text(5.52, 3.04, "D", 3.8, "bold", hjust = 0) +
  label_text(5.92, 3.04, "Post-MS survival", 2.90, "bold", hjust = 0) +
  label_text(5.62, 2.57, "fitted function", 2.50, "bold", hjust = 0, colour = "#666666") +
  label_text(6.88, 2.05, "s(N, \u0394N)", 3.00, "bold") +
  label_text(
    6.88, 1.35,
    "ploidy-dependent filtering\nof altered daughters",
    2.50
  ) +
  arrow_segment(6.88, 6.88, 3.88, 3.32, colour = dark)

# E: predicted outcome and observed fit targets
p <- p +
  box(9.55, 15.65, 1.05, 9.35, type_cols[["output"]]) +
  label_text(9.78, 9.07, "E", 3.8, "bold", hjust = 0) +
  label_text(10.18, 9.07, "Outcome & fit comparison", 3.00, "bold", hjust = 0) +
  label_text(9.88, 8.58, "predicted output", 2.50, "bold", hjust = 0, colour = "#666666") +
  label_text(
    12.60, 8.02,
    "Chromosome-state distribution",
    2.90, "bold"
  ) +
  annotate(
    "rect", xmin = 10.52, xmax = 11.00, ymin = 5.40, ymax = 6.62,
    fill = blue, colour = blue
  ) +
  annotate(
    "rect", xmin = 11.17, xmax = 11.65, ymin = 5.40, ymax = 7.58,
    fill = blue, colour = blue
  ) +
  annotate(
    "rect", xmin = 11.82, xmax = 12.30, ymin = 5.40, ymax = 6.96,
    fill = "#56B4E9", colour = "#56B4E9"
  ) +
  annotate(
    "rect", xmin = 12.47, xmax = 12.95, ymin = 5.40, ymax = 6.05,
    fill = "#E69F00", colour = "#E69F00"
  ) +
  annotate(
    "rect", xmin = 13.12, xmax = 13.60, ymin = 5.40, ymax = 6.78,
    fill = verm, colour = verm
  ) +
  annotate(
    "segment", x = 10.35, xend = 13.80, y = 5.40, yend = 5.40,
    colour = dark, linewidth = 0.42
  ) +
  label_text(12.08, 4.98, "chromosome number, N", 2.50) +
  label_text(
    12.60, 4.45,
    "predicted growth and live burden\n+ ploidy summaries",
    2.50, "bold"
  ) +
  annotate(
    "segment", x = 9.90, xend = 15.30, y = 3.95, yend = 3.95,
    colour = "#A0A0A0", linewidth = 0.40, linetype = "dashed"
  ) +
  label_text(9.88, 3.60, "observed fit targets", 2.50, "bold", hjust = 0, colour = "#666666") +
  label_text(
    10.05, 3.03,
    "In vitro",
    2.58, "bold", hjust = 0, colour = vitro_col
  ) +
  label_text(
    11.15, 3.03,
    "growth + karyotype/flow",
    2.50, hjust = 0
  ) +
  label_text(
    10.05, 2.42,
    "In vivo",
    2.58, "bold", hjust = 0, colour = vivo_col
  ) +
  label_text(
    11.15, 2.42,
    "tumor burden + terminal karyotype",
    2.50, hjust = 0
  ) +
  label_text(
    12.60, 1.65,
    "measured burden is a fit target,\nnot a resource-state input",
    2.50, colour = "#555555"
  ) +
  arrow_segment(8.45, 9.55, 8.35, 7.80, colour = dark) +
  arrow_segment(8.45, 9.55, 5.20, 6.30, colour = dark) +
  arrow_segment(8.45, 9.55, 2.05, 4.80, colour = dark)

# In-vivo resource feedback uses the simulated live population, not measured burden.
p <- p +
  annotate(
    "segment", x = 10.15, xend = 9.00, y = 1.05, yend = 0.36,
    colour = vivo_col, linewidth = 0.60
  ) +
  annotate(
    "segment", x = 9.00, xend = 4.30, y = 0.36, yend = 0.36,
    colour = vivo_col, linewidth = 0.60
  ) +
  arrow_segment(4.30, 4.30, 0.36, 5.48, colour = vivo_col, linetype = "dashed", linewidth = 0.60) +
  label_text(
    7.05, 0.18,
    "simulated live-state demand feeds latent tumor O\u2082",
    2.50, "bold", colour = vivo_col
  )

dir.create(file.path(draft_root, "initial_subpanels", "F2"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(draft_root, "refined_subpanels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(draft_root, "final_figures", "recommended"), recursive = TRUE, showWarnings = FALSE)

save_plot <- function(stem) {
  ggsave(
    paste0(stem, ".png"), p, width = 7.1, height = 4.8,
    units = "in", dpi = 300, bg = "white"
  )
  ggsave(
    paste0(stem, ".pdf"), p, width = 7.1, height = 4.8,
    units = "in", device = grDevices::cairo_pdf, bg = "white"
  )
}

save_plot(file.path(draft_root, "initial_subpanels", "F2", "model_overview"))
save_plot(file.path(draft_root, "refined_subpanels", "figure2_recommended"))
save_plot(file.path(draft_root, "final_figures", "recommended", "figure2"))

message("Figure 2 draft written under: ", draft_root)
