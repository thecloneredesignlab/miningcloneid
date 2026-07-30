#!/usr/bin/env Rscript

if (identical(Sys.getenv("FIGURE5_DRAW_WORKER"), "1")) {

# Standalone iteration3 Figure 5 generator.
#
# Scope:
#   A. Joint-fitting workflow.
#   B. Pooled separate-fit t-SNE landscape used to select joint warm starts.
#   C. Post-missegregation survival functions across the six selected winners.
#   D. Proliferation and stress-linked missegregation functions.
#   E. Direct in-sample observed-predicted summaries for both fitted contexts.
#   F. Full parameter-ratio view across the six approved July joint-fit winners.
#
# This script is self-contained: it does not source functions from another
# checkout. It reads frozen selected-winner tables stored in iteration3 and
# canonical pooled t-SNE coordinate tables from soft_coupling. It does not
# refit the model or treat across-winner variation as a confidence interval.

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c(
  "ggplot2", "dplyr", "patchwork", "scales", "ggnewscale", "shadowtext"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(ggnewscale)
  library(shadowtext)
})

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
})
workspace_root <- normalizePath(
  Sys.getenv(
    "FIGURE_WORKSPACE_ROOT",
    unset = file.path(script_dir, "..", "..")
  ),
  mustWork = TRUE
)
repo_root <- normalizePath(
  Sys.getenv(
    "HYPOXIA_REPO_ROOT",
    unset = file.path(script_dir, "..", "..", "..", "..", "..", "..")
  ),
  mustWork = TRUE
)
figure5_data_dir <- normalizePath(
  Sys.getenv(
    "FIGURE5_DATA_DIR",
    unset = file.path(
      workspace_root, "data", "Figures", "Figure5"
    )
  ),
  mustWork = TRUE
)
frozen_root <- file.path(figure5_data_dir, "figure5_frozen_inputs")
selection_path <- file.path(frozen_root, "selected_results.tsv")
parameter_meta_path <- file.path(figure5_data_dir, "parameter_function_groups.tsv")
parameter_palette_path <- file.path(
  figure5_data_dir, "parameter_function_group_palette.tsv"
)

missing_configuration <- c(
  selection_path, parameter_meta_path, parameter_palette_path
)[!file.exists(c(selection_path, parameter_meta_path, parameter_palette_path))]
if (length(missing_configuration)) {
  stop(
    "Missing Figure 5 selection/group configuration:\n",
    paste(missing_configuration, collapse = "\n")
  )
}

figure_root <- normalizePath(
  Sys.getenv(
    "FIGURE5_FIGURE_DIR",
    unset = file.path(workspace_root, "Figures")
  ),
  mustWork = FALSE
)
panel_root <- normalizePath(
  Sys.getenv(
    "FIGURE5_PANEL_DIR",
    unset = file.path(figure_root, "figure5_subpanels")
  ),
  mustWork = FALSE
)
revised_root <- figure_root
allowed_output_dirs <- c(figure_root, panel_root, revised_root)
invisible(lapply(
  allowed_output_dirs,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

selected <- read.delim(
  selection_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
selected <- selected[selected$record_type == "joint_pair_best", , drop = FALSE]
parameter_meta_all <- read.delim(
  parameter_meta_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
parameter_palette_meta <- read.delim(
  parameter_palette_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
parameter_palette_meta <- parameter_palette_meta[
  order(parameter_palette_meta$group_order),
  ,
  drop = FALSE
]
if (nrow(parameter_meta_all) != 18L ||
    nrow(parameter_palette_meta) != 5L ||
    !identical(parameter_palette_meta$group_order, seq_len(5L)) ||
    !setequal(
      unique(parameter_meta_all$parameter_group),
      parameter_palette_meta$parameter_group
    )) {
  stop("Figure 5 requires the shared ordered five-group parameter taxonomy")
}
parameter_group_colors <- setNames(
  parameter_palette_meta$color,
  parameter_palette_meta$parameter_group
)

required_selection_fields <- c(
  "warmup_label",
  "invivo_seed",
  "invitro_seed",
  "selected_seed",
  "objective",
  "bundle_dir"
)
if (!all(required_selection_fields %in% names(selected))) {
  stop(
    "Selection table lacks required fields: ",
    paste(setdiff(required_selection_fields, names(selected)), collapse = ", ")
  )
}
if (nrow(selected) != 6L) {
  stop("Figure 5 requires exactly six approved joint-pair winners; found ", nrow(selected))
}
if (length(unique(selected$invitro_seed)) != 1L ||
    unique(selected$invitro_seed) != "seed10") {
  stop("The approved Figure 5 universe must use in-vitro anchor seed10")
}

extract_region <- function(warmup_label) {
  hit <- regmatches(warmup_label, regexpr("C[0-9]{2}", warmup_label))
  ifelse(nzchar(hit), hit, "unclassified")
}

selected$region <- extract_region(selected$warmup_label)
selected$winner_label <- paste0(
  selected$region,
  " / ",
  selected$invivo_seed,
  " / ",
  selected$selected_seed
)
selected$fit_dir <- file.path(
  frozen_root, "winners", selected$warmup_label
)

missing_fit_dirs <- selected$fit_dir[!dir.exists(selected$fit_dir)]
if (length(missing_fit_dirs)) {
  stop("Missing selected fit directories:\n", paste(missing_fit_dirs, collapse = "\n"))
}

read_winner_table <- function(fit_dir, relative_path) {
  path <- file.path(fit_dir, relative_path)
  if (!file.exists(path)) {
    stop("Missing selected-winner source table: ", path)
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

context_colors <- c(
  "In vivo" = "#0072B2",
  "In vitro" = "#CC79A7"
)
region_colors <- c(
  "C01" = "#FFD700",
  "C02" = "#6A3D9A",
  "C03" = "#39FF14"
)
cluster_hull_linetype <- "solid"
cluster_label_text_color <- "#000000"
cluster_label_fontface <- "plain"
cluster_label_border_linewidth <- 0.33
tsne_coordinate_label_scale <- 0.5
in_vivo_color <- unname(context_colors[["In vivo"]])
in_vitro_color <- unname(context_colors[["In vitro"]])
best_seed_border_color <- "#FFFFFF"
ratio_marker_border_color <- "#FFFFFF"
d_reference_linewidth <- 0.275
d_range_linewidth <- 0.625
efg_individual_linewidth <- 0.19
efg_median_linewidth <- 0.56
figure_font_family <- "Helvetica"
panel_title_size <- 7.5
axis_title_size <- 6.1
axis_text_size <- 5.8
strip_text_size <- 5.8
panel_title_geom_size <- panel_title_size / ggplot2::.pt
a_height_scale <- 0.80
f_height_scale <- 1.20
assembled_height_inches <- 9.44
aligned_row_heights <- c(
  header = 0.08,
  first_plot = 0.41,
  second_plot = 0.41,
  footer = 0.10
)
c_grid_heights <- c(
  mechanism_strip = 0.055,
  first_plot = 0.4225,
  second_plot = 0.4225,
  shared_x_axis = 0.10
)
d_shared_x_positions <- c(0, 0.46, 0.54, 1)
cde_content_heights <- c(
  panel_c = 0.56,
  panels_d_e = 0.44
)
bcd_visual_scale <- c(
  width = 0.99,
  height = 0.93
)
e_layout_heights <- c(
  header = 0.08,
  content = 0.87,
  note = 0.05
)
row2_block_widths <- c(
  panels_bcd = 1.00,
  panel_e = 1.00
)
f_ratio_layout_heights <- c(
  header = 0.010,
  body = 0.910,
  note = 0.080
)
ratio_axis_inset_bounds <- c(
  left = 0.208,
  right = 0.994,
  bottom = 0.000,
  top = 0.055
)
ratio_axis_side_offset <- 2.10
f_legend_top_margin <- -5
de_legend_vertical_shift <- 0.42
context_linetypes <- c(
  "In vivo" = "solid",
  "In vitro" = "solid"
)

theme_manuscript <- function(base_size = 9) {
  theme_bw(
    base_size = base_size,
    base_family = figure_font_family
  ) +
    theme(
      text = element_text(family = figure_font_family),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E6E8EB", linewidth = 0.28),
      panel.border = element_rect(color = "#4B5563", linewidth = 0.45),
      axis.text = element_text(
        family = figure_font_family,
        size = axis_text_size,
        color = "#27313A"
      ),
      axis.title = element_text(
        family = figure_font_family,
        size = axis_title_size,
        color = "#1F2933"
      ),
      strip.background = element_rect(fill = "#F3F4F6", color = "#9CA3AF"),
      strip.text = element_text(
        family = figure_font_family,
        size = strip_text_size,
        face = "bold",
        color = "#27313A"
      ),
      legend.key.height = grid::unit(0.28, "cm"),
      legend.key.width = grid::unit(0.46, "cm"),
      legend.title = element_text(
        family = figure_font_family,
        face = "bold"
      ),
      plot.title = element_text(
        family = figure_font_family,
        size = panel_title_size,
        face = "bold",
        color = "#111827"
      ),
      plot.subtitle = element_text(color = "#4B5563"),
      plot.caption = element_text(color = "#5F6B76", hjust = 0),
      plot.margin = margin(3, 4, 3, 4)
    )
}

save_plot_pair <- function(plot, stub, width, height, dpi = 220) {
  pdf_path <- paste0(stub, ".pdf")
  png_path <- paste0(stub, ".png")
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  invisible(c(pdf = pdf_path, png = png_path))
}

make_text_band <- function(
  label,
  size = 6.2,
  face = "bold",
  hjust = 0.5,
  x = if (hjust == 0) 0 else 0.5,
  y = 0.5
) {
  ggplot() +
    annotate(
      "text",
      x = x,
      y = y,
      label = label,
      hjust = hjust,
      size = size / ggplot2::.pt,
      family = figure_font_family,
      fontface = face,
      lineheight = 0.90,
      color = "#27313A"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void(base_family = figure_font_family) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_strip_band <- function(label, size = strip_text_size) {
  ggplot() +
    annotate(
      "rect",
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      fill = "#F3F4F6",
      color = "#9CA3AF",
      linewidth = 0.30
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      size = size / ggplot2::.pt,
      family = figure_font_family,
      fontface = "bold",
      lineheight = 0.90,
      color = "#27313A"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void(base_family = figure_font_family) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_vertical_strip_band <- function(label, size = 6.2) {
  ggplot() +
    annotate(
      "rect",
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      fill = "#F3F4F6",
      color = "#9CA3AF",
      linewidth = 0.30
    ) +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      angle = 270,
      size = size / ggplot2::.pt,
      family = figure_font_family,
      fontface = "bold",
      color = "#27313A"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = figure_font_family) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_vertical_axis_title <- function(label, size = axis_title_size) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      parse = TRUE,
      angle = 90,
      size = size / ggplot2::.pt,
      family = figure_font_family,
      color = "#1F2933"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = figure_font_family) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# -------------------------------------------------------------------------
# Panel A: joint-fitting workflow.
# -------------------------------------------------------------------------

workflow_context_boxes <- data.frame(
  xmin = c(1, 1, 21, 21),
  xmax = c(16, 16, 36, 36),
  ymin = c(14.2, 4.2, 14.2, 4.2),
  ymax = c(22.2, 12.2, 22.2, 12.2),
  context = c("In vivo", "In vitro", "In vivo", "In vitro"),
  title = c(
    "Separate\nin vivo fit",
    "Separate\nin vitro fit",
    "Six landscape\nrepresentatives",
    "Common seed-10\nculture anchor"
  ),
  detail = c(
    "500-seed search",
    "500-seed search",
    "3 clusters x 2 regions",
    "shared by all six pairs"
  ),
  stringsAsFactors = FALSE
)
workflow_context_boxes$color <- context_colors[workflow_context_boxes$context]

workflow_core_boxes <- data.frame(
  xmin = c(40.5, 57.5, 77.0, 92.0),
  xmax = c(52.5, 72.0, 87.0, 99.5),
  ymin = c(8.1, 8.1, 8.1, 8.1),
  ymax = c(18.3, 18.3, 18.3, 18.3),
  title = c(
    "Six warm-start\npairs",
    "14 paired\nparameters",
    "3,000 joint\nsearches",
    "Six pair-level\nwinners"
  ),
  detail = c(
    "landscape-informed",
    "",
    "500 starts / pair",
    "pair minimum"
  ),
  stringsAsFactors = FALSE
)

workflow_output_boxes <- data.frame(
  xmin = c(58, 72.2, 84.5),
  xmax = c(70, 82.5, 99),
  ymin = c(0.2, 0.2, 0.2),
  ymax = c(4.2, 4.2, 4.2),
  title = c("E  Fit adequacy", "F  Parameter ratios", "C-D  Fitted functions"),
  stringsAsFactors = FALSE
)

workflow_lane_arrows <- data.frame(
  x = c(16.5, 16.5),
  xend = c(20.5, 20.5),
  y = c(18.2, 8.2),
  yend = c(18.2, 8.2)
)
workflow_merge_arrows <- data.frame(
  x = c(36.6, 36.6),
  xend = c(39.9, 39.9),
  y = c(18.2, 8.2),
  yend = c(15.7, 10.7)
)
workflow_core_arrows <- data.frame(
  x = c(53.0, 72.5, 87.5),
  xend = c(57.0, 76.5, 91.5),
  y = 13.2,
  yend = 13.2
)
workflow_output_branches <- data.frame(
  x = c(64, 77.35, 91.75),
  xend = c(64, 77.35, 91.75),
  y = 6.70,
  yend = 4.45
)
workflow_main_arrowhead_length <- 0.072
workflow_branch_arrowhead_length <- 0.060
workflow_welsch_layout <- c(
  left_bar_x = 61.7,
  right_bar_x = 67.8,
  bar_ymin = 12.0,
  bar_ymax = 13.5,
  connector_y = 12.75,
  label_x = 64.75,
  label_y = 10.95,
  detail_y = 9.65
)
workflow_welsch_label_size <- 1.75
workflow_welsch_detail_size <- 1.65

p_a <- ggplot() +
  geom_rect(
    data = workflow_context_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = alpha(workflow_context_boxes$color, 0.10),
    color = workflow_context_boxes$color,
    linewidth = 0.55
  ) +
  geom_rect(
    data = workflow_core_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#F8FAFC",
    color = "#94A3B8",
    linewidth = 0.5
  ) +
  geom_rect(
    data = workflow_output_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#F3F4F6",
    color = "#CBD1D8",
    linewidth = 0.4
  ) +
  geom_segment(
    data = workflow_lane_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "#475569",
    linewidth = 0.6,
    arrow = grid::arrow(
      length = grid::unit(workflow_main_arrowhead_length, "in"),
      type = "closed"
    )
  ) +
  geom_segment(
    data = workflow_merge_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "#475569",
    linewidth = 0.6,
    arrow = grid::arrow(
      length = grid::unit(workflow_main_arrowhead_length, "in"),
      type = "closed"
    )
  ) +
  geom_segment(
    data = workflow_core_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "#475569",
    linewidth = 0.6,
    arrow = grid::arrow(
      length = grid::unit(workflow_main_arrowhead_length, "in"),
      type = "closed"
    )
  ) +
  annotate(
    "segment",
    x = 95.75, xend = 95.75, y = 8.0, yend = 6.70,
    color = "#64748B",
    linewidth = 0.42
  ) +
  annotate(
    "segment",
    x = 64, xend = 95.75, y = 6.70, yend = 6.70,
    color = "#64748B",
    linewidth = 0.42
  ) +
  geom_segment(
    data = workflow_output_branches,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "#64748B",
    linewidth = 0.42,
    arrow = grid::arrow(
      length = grid::unit(workflow_branch_arrowhead_length, "in"),
      type = "closed"
    )
  ) +
  geom_text(
    data = workflow_context_boxes,
    aes(x = (xmin + xmax) / 2, y = ymax - 2.35, label = title),
    color = workflow_context_boxes$color,
    fontface = "bold",
    size = 1.88,
    lineheight = 0.92
  ) +
  geom_text(
    data = workflow_context_boxes,
    aes(x = (xmin + xmax) / 2, y = ymin + 2.15, label = detail),
    color = "#475569",
    size = 1.62
  ) +
  geom_text(
    data = workflow_core_boxes,
    aes(x = (xmin + xmax) / 2, y = ymax - 2.45, label = title),
    color = "#0F172A",
    fontface = "bold",
    size = 1.82,
    lineheight = 0.95
  ) +
  geom_text(
    data = workflow_core_boxes,
    aes(x = (xmin + xmax) / 2, y = ymin + 2.1, label = detail),
    color = "#475569",
    size = 1.52,
    lineheight = 0.95
  ) +
  annotate(
    "segment",
    x = workflow_welsch_layout[["left_bar_x"]],
    xend = workflow_welsch_layout[["left_bar_x"]],
    y = workflow_welsch_layout[["bar_ymin"]],
    yend = workflow_welsch_layout[["bar_ymax"]],
    color = in_vivo_color, linewidth = 1.45, lineend = "round"
  ) +
  annotate(
    "segment",
    x = workflow_welsch_layout[["right_bar_x"]],
    xend = workflow_welsch_layout[["right_bar_x"]],
    y = workflow_welsch_layout[["bar_ymin"]],
    yend = workflow_welsch_layout[["bar_ymax"]],
    color = in_vitro_color, linewidth = 1.45, lineend = "round"
  ) +
  annotate(
    "segment",
    x = workflow_welsch_layout[["left_bar_x"]],
    xend = workflow_welsch_layout[["right_bar_x"]],
    y = workflow_welsch_layout[["connector_y"]],
    yend = workflow_welsch_layout[["connector_y"]],
    color = "#64748B", linewidth = 0.42, linetype = "22"
  ) +
  annotate(
    "text",
    x = workflow_welsch_layout[["label_x"]],
    y = workflow_welsch_layout[["label_y"]],
    label = "Welsch penalty",
    color = "#334155",
    size = workflow_welsch_label_size,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = workflow_welsch_layout[["label_x"]],
    y = workflow_welsch_layout[["detail_y"]],
    label = "soft-coupled",
    color = "#475569",
    size = workflow_welsch_detail_size
  ) +
  geom_text(
    data = workflow_output_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = title),
    color = "#334155",
    fontface = "bold",
    size = 1.75
  ) +
  annotate(
    "text",
    x = 27.5, y = 1.8,
    label = "Same chromosome-state model; context-specific observation models",
    color = "#64748B", size = 1.85, fontface = "italic"
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 23), clip = "off") +
  labs(title = "A  Joint-fitting workflow") +
  theme_void(base_size = 8, base_family = figure_font_family) +
  theme(
    plot.title = element_text(
      family = figure_font_family,
      size = panel_title_size,
      face = "bold",
      color = "#111827"
    ),
    plot.margin = margin(2, 3, 1, 12.7)
  )

# -------------------------------------------------------------------------
# Panel E: compact, direct in-sample fit adequacy across all six winners.
# -------------------------------------------------------------------------

burden_rows <- list()
terminal_rows <- list()
growth_rows <- list()
kary_rows <- list()

for (i in seq_len(nrow(selected))) {
  fit_dir <- selected$fit_dir[[i]]
  winner <- selected$winner_label[[i]]
  region <- selected$region[[i]]

  burden <- read_winner_table(fit_dir, "invivo_burden_fit.tsv")
  burden <- burden[
    is.finite(burden$obs_norm) &
      is.finite(burden$pred_norm) &
      burden$day > 0,
    ,
    drop = FALSE
  ]
  burden$winner <- winner
  burden$region <- region
  burden$obs_id <- paste(burden$harvest, burden$day, sep = "|")
  burden_rows[[i]] <- burden[, c(
    "winner", "region", "obs_id", "obs_norm", "pred_norm"
  )]

  terminal <- read_winner_table(fit_dir, "invivo_terminal_ploidy_fit.tsv")
  terminal_split <- split(terminal, terminal$harvest)
  terminal_summary <- do.call(rbind, lapply(terminal_split, function(x) {
    observed_total <- sum(x$obs_count, na.rm = TRUE)
    predicted_total <- sum(x$pred_fraction, na.rm = TRUE)
    data.frame(
      obs_id = x$harvest[[1]],
      cohort = x$cohort[[1]],
      observed_mean_N = if (observed_total > 0) {
        sum(x$N * x$obs_count, na.rm = TRUE) / observed_total
      } else {
        NA_real_
      },
      predicted_mean_N = if (predicted_total > 0) {
        sum(x$N * x$pred_fraction, na.rm = TRUE) / predicted_total
      } else {
        NA_real_
      }
    )
  }))
  terminal_summary$winner <- winner
  terminal_summary$region <- region
  terminal_rows[[i]] <- terminal_summary

  growth <- read_winner_table(fit_dir, "invitro_growth_loglik.tsv")
  growth <- growth[
    is.finite(growth$observed_growth) &
      is.finite(growth$predicted_growth_rate),
    ,
    drop = FALSE
  ]
  growth$winner <- winner
  growth$region <- region
  growth$obs_id <- paste(
    growth$segment_id,
    growth$passage_id,
    growth$lineage_terminal_key,
    growth$lineage_passage_index,
    sep = "|"
  )
  growth_rows[[i]] <- growth[, c(
    "winner",
    "region",
    "obs_id",
    "observed_growth",
    "predicted_growth_rate"
  )]

  kary <- read_winner_table(fit_dir, "invitro_lineage_summary.tsv")
  kary <- kary[
    is.finite(kary$observed_mean_kary_N) &
      is.finite(kary$predicted_mean_kary_N),
    ,
    drop = FALSE
  ]
  kary$winner <- winner
  kary$region <- region
  kary$obs_id <- paste(
    kary$segment_id,
    kary$passage_id,
    kary$lineage_terminal_key,
    kary$lineage_passage_index,
    sep = "|"
  )
  kary_rows[[i]] <- kary[, c(
    "winner",
    "region",
    "obs_id",
    "observed_mean_kary_N",
    "predicted_mean_kary_N"
  )]
}

burden_fit <- bind_rows(burden_rows)
terminal_fit <- bind_rows(terminal_rows)
growth_fit <- bind_rows(growth_rows)
kary_fit <- bind_rows(kary_rows)

assert_six_per_observation <- function(df, label) {
  counts <- table(df$obs_id)
  if (length(counts) == 0L || any(counts != 6L)) {
    stop(label, " does not contain exactly six winner predictions per observation")
  }
}
assert_six_per_observation(burden_fit, "In-vivo burden fit")
assert_six_per_observation(terminal_fit, "In-vivo terminal-ploidy fit")
assert_six_per_observation(growth_fit, "In-vitro growth fit")
assert_six_per_observation(kary_fit, "In-vitro karyotype fit")

make_fit_scatter <- function(
  df,
  observed_col,
  predicted_col,
  title,
  x_label,
  y_label,
  color,
  limits,
  show_range = FALSE
) {
  summary_df <- df %>%
    group_by(obs_id) %>%
    summarise(
      observed = first(.data[[observed_col]]),
      predicted = median(.data[[predicted_col]], na.rm = TRUE),
      predicted_min = min(.data[[predicted_col]], na.rm = TRUE),
      predicted_max = max(.data[[predicted_col]], na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot() +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "#6B7280",
      linetype = "22",
      linewidth = 0.5
    )
  if (show_range) {
    p <- p + geom_segment(
      data = summary_df,
      aes(
        x = observed,
        xend = observed,
        y = predicted_min,
        yend = predicted_max
      ),
      color = alpha(color, 0.48),
      linewidth = 0.55
    )
  }
  p +
    geom_point(
      data = df,
      aes(x = .data[[observed_col]], y = .data[[predicted_col]]),
      color = alpha(color, 0.18),
      size = 0.9
    ) +
    geom_point(
      data = summary_df,
      aes(x = observed, y = predicted),
      shape = 21,
      fill = color,
      color = "white",
      stroke = 0.25,
      size = 1.6
    ) +
    coord_fixed(
      xlim = limits,
      ylim = limits,
      ratio = 1,
      expand = FALSE
    ) +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme_manuscript(base_size = 7) +
    theme(
      plot.title = element_text(
        family = figure_font_family,
        size = panel_title_size
      ),
      axis.title = element_text(
        family = figure_font_family,
        size = axis_title_size
      ),
      axis.text = element_text(
        family = figure_font_family,
        size = axis_text_size
      ),
      plot.margin = margin(2, 3, 2, 3)
    )
}

p_a_burden <- make_fit_scatter(
  burden_fit,
  "obs_norm",
  "pred_norm",
  "In vivo | normalized burden",
  "Observed",
  "Predicted",
  in_vivo_color,
  c(0, 1),
  show_range = FALSE
)
p_a_terminal <- make_fit_scatter(
  terminal_fit,
  "observed_mean_N",
  "predicted_mean_N",
  "In vivo | terminal mean N",
  "Observed mean N",
  "Predicted mean N",
  in_vivo_color,
  c(38, 57),
  show_range = TRUE
)
p_a_growth <- make_fit_scatter(
  growth_fit,
  "observed_growth",
  "predicted_growth_rate",
  "In vitro | growth rate",
  "Observed rate",
  "Predicted rate",
  in_vitro_color,
  c(0.2, 2.2),
  show_range = FALSE
)
p_a_kary <- make_fit_scatter(
  kary_fit,
  "observed_mean_kary_N",
  "predicted_mean_kary_N",
  "In vitro | mean N",
  "Observed mean N",
  "Predicted mean N",
  in_vitro_color,
  c(44, 92),
  show_range = TRUE
)

p_b_header <- make_text_band(
  "E  In-sample fit adequacy",
  size = panel_title_size,
  hjust = 0,
  y = 0.70
)
p_b_content <- wrap_plots(
  p_a_burden,
  p_a_terminal,
  p_a_growth,
  p_a_kary,
  ncol = 2,
  nrow = 2,
  widths = c(1, 1),
  heights = c(1, 1)
)
p_b_note <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 0.90,
    label = paste0(
      "Pale: six selected fits; filled: median. ",
      "Bars: selected-fit range where legible."
    ),
    hjust = 0,
    size = 2.0,
    color = "#4B5563"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 3, 0, 3))
p_b <- wrap_plots(
  p_b_header,
  p_b_content,
  p_b_note,
  ncol = 1,
  heights = unname(e_layout_heights)
)

# -------------------------------------------------------------------------
# Panel B: pooled separate-fit t-SNE landscape used for warm-start selection.
# -------------------------------------------------------------------------

tsne_full_path <- Sys.getenv(
  "JOINT_TSNE_FULL_COORDINATES",
  unset = file.path(
    dirname(repo_root),
    "soft_coupling",
    "oxygen",
    "results",
    "analysis",
    "best_fit_parameter_feature",
    "02_parameter_landscape_clustering",
    "pooled_invivo_invitro",
    "full_data_in_vivo_clustring",
    "Tables",
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_full_coordinates.csv"
  )
)
tsne_best_path <- Sys.getenv(
  "JOINT_TSNE_BEST_COORDINATES",
  unset = file.path(
    dirname(repo_root),
    "soft_coupling",
    "oxygen",
    "results",
    "analysis",
    "best_fit_parameter_feature",
    "02_parameter_landscape_clustering",
    "pooled_invivo_invitro",
    "full_data_in_vivo_clustring",
    "Tables",
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
  )
)
if (!file.exists(tsne_full_path) || !file.exists(tsne_best_path)) {
  stop(
    "Missing canonical pooled t-SNE coordinates. Checked:\n",
    tsne_full_path, "\n", tsne_best_path
  )
}

tsne_full <- read.csv(tsne_full_path, check.names = FALSE, stringsAsFactors = FALSE)
tsne_best <- read.csv(tsne_best_path, check.names = FALSE, stringsAsFactors = FALSE)
required_tsne_fields <- c("tSNE1", "tSNE2", "dataset", "point_type", "seed", "objective")
if (!all(required_tsne_fields %in% names(tsne_full)) ||
    !all(c(required_tsne_fields, "cluster_id") %in% names(tsne_best))) {
  stop("Canonical t-SNE coordinate tables lack required columns")
}
if (nrow(tsne_full) != 228000L ||
    sum(tsne_full$point_type == "best" & tsne_full$dataset == "invivo") != 500L ||
    sum(tsne_full$point_type == "best" & tsne_full$dataset == "invitro") != 500L) {
  stop("Unexpected pooled t-SNE universe; expected 228,000 rows and 500 best fits/context")
}

tsne_initial <- tsne_full[tsne_full$point_type == "initial", , drop = FALSE]
tsne_initial$context <- ifelse(tsne_initial$dataset == "invivo", "In vivo", "In vitro")
tsne_best_vivo <- tsne_best[tsne_best$dataset == "invivo", , drop = FALSE]
tsne_best_vitro <- tsne_best[tsne_best$dataset == "invitro", , drop = FALSE]

cluster_hulls <- do.call(rbind, lapply(
  split(tsne_best_vivo, tsne_best_vivo$cluster_id),
  function(d) {
    xy <- unique(d[, c("tSNE1", "tSNE2"), drop = FALSE])
    if (nrow(xy) < 3L) return(NULL)
    index <- grDevices::chull(xy$tSNE1, xy$tSNE2)
    hull <- xy[c(index, index[[1L]]), , drop = FALSE]
    center <- colMeans(hull)
    hull$tSNE1 <- center[[1L]] + 1.035 * (hull$tSNE1 - center[[1L]])
    hull$tSNE2 <- center[[2L]] + 1.035 * (hull$tSNE2 - center[[2L]])
    hull$cluster_id <- as.character(d$cluster_id[[1L]])
    hull
  }
))
cluster_ids <- sort(unique(as.character(cluster_hulls$cluster_id)))
cluster_colors <- setNames(
  unname(region_colors[sub("^vi_", "", cluster_ids)]),
  cluster_ids
)
if (any(is.na(cluster_colors))) {
  stop("Figure 5E contains a warm-start region without a declared color")
}
cluster_labels <- data.frame(
  cluster_id = cluster_ids,
  label = sub("^vi_", "", cluster_ids),
  point_overlap_count = 0L,
  stringsAsFactors = FALSE
)

x_range <- range(tsne_full$tSNE1, finite = TRUE)
y_range <- range(tsne_full$tSNE2, finite = TRUE)
span <- max(diff(x_range), diff(y_range)) * 1.08
x_mid <- mean(x_range)
y_mid <- mean(y_range)
tsne_xlim <- x_mid + c(-0.5, 0.5) * span
tsne_ylim <- y_mid + c(-0.5, 0.5) * span
tsne_axis_origin <- c(
  x = tsne_xlim[[1L]] + 0.025 * span,
  y = tsne_ylim[[1L]] + 0.025 * span
)
tsne_axis_length <- 0.115 * span

p_c_scatter <- ggplot() +
  geom_point(
    data = tsne_initial,
    aes(x = tSNE1, y = tSNE2, color = context, shape = context),
    size = 0.12,
    alpha = 0.38,
    stroke = 0,
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = context_colors,
    breaks = c("In vivo", "In vitro"),
    name = "Initial samples",
    guide = guide_legend(
      order = 1,
      override.aes = list(shape = c(16, 17), size = 2.5, alpha = 0.85)
    )
  ) +
  scale_shape_manual(
    values = c("In vivo" = 16, "In vitro" = 17),
    guide = "none"
  ) +
  geom_point(
    data = tsne_best_vivo,
    aes(x = tSNE1, y = tSNE2, fill = objective),
    shape = 21,
    size = 1.15,
    color = best_seed_border_color,
    stroke = 0.32,
    alpha = 0.92
  ) +
  scale_fill_gradient(
    low = "#2C7BB6",
    high = "#FDE725",
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = tsne_best_vitro,
    aes(x = tSNE1, y = tSNE2, fill = objective),
    shape = 24,
    size = 1.28,
    color = best_seed_border_color,
    stroke = 0.32,
    alpha = 0.94
  ) +
  scale_fill_gradient(
    low = "#1A9850",
    high = "#D73027",
    guide = "none"
  )
for (cluster_id in cluster_ids) {
  hull <- cluster_hulls[cluster_hulls$cluster_id == cluster_id, , drop = FALSE]
  p_c_scatter <- p_c_scatter + geom_path(
    data = hull,
    aes(x = tSNE1, y = tSNE2),
    inherit.aes = FALSE,
    color = cluster_colors[[cluster_id]],
    linewidth = 0.48,
    linetype = cluster_hull_linetype,
    show.legend = FALSE
  )
}
p_c_scatter <- p_c_scatter +
  annotate(
    "segment",
    x = tsne_axis_origin[["x"]],
    xend = tsne_axis_origin[["x"]] + tsne_axis_length,
    y = tsne_axis_origin[["y"]],
    yend = tsne_axis_origin[["y"]],
    color = "#4B5563",
    linewidth = 0.42,
    lineend = "round",
    arrow = grid::arrow(
      length = grid::unit(0.025, "inches"),
      type = "closed"
    )
  ) +
  annotate(
    "segment",
    x = tsne_axis_origin[["x"]],
    xend = tsne_axis_origin[["x"]],
    y = tsne_axis_origin[["y"]],
    yend = tsne_axis_origin[["y"]] + tsne_axis_length,
    color = "#4B5563",
    linewidth = 0.42,
    lineend = "round",
    arrow = grid::arrow(
      length = grid::unit(0.025, "inches"),
      type = "closed"
    )
  ) +
  annotate(
    "text",
    x = tsne_axis_origin[["x"]] + tsne_axis_length + 0.012 * span,
    y = tsne_axis_origin[["y"]],
    label = "t-SNE 1",
    hjust = 0,
    size = axis_title_size * tsne_coordinate_label_scale / ggplot2::.pt,
    family = figure_font_family,
    color = "#4B5563"
  ) +
  annotate(
    "text",
    x = tsne_axis_origin[["x"]],
    y = tsne_axis_origin[["y"]] + tsne_axis_length + 0.012 * span,
    label = "t-SNE 2",
    hjust = 0,
    vjust = 0,
    size = axis_title_size * tsne_coordinate_label_scale / ggplot2::.pt,
    family = figure_font_family,
    color = "#4B5563"
  ) +
  coord_equal(xlim = tsne_xlim, ylim = tsne_ylim, expand = FALSE, clip = "on") +
  labs(x = NULL, y = NULL) +
  theme_manuscript(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5)
  )

make_objective_colorbar <- function(values, low, high, title) {
  value_range <- range(values, finite = TRUE)
  bar_data <- data.frame(
    objective = seq(value_range[[1L]], value_range[[2L]], length.out = 256L),
    y = 1
  )
  ggplot(bar_data, aes(x = objective, y = y, fill = objective)) +
    geom_tile(width = diff(value_range) / 255, height = 0.18) +
    scale_fill_gradient(low = low, high = high, guide = "none") +
    scale_x_continuous(
      breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0.89, 1.11), expand = c(0, 0)) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(
        title,
        " ",
        sprintf("%.1f", value_range[[1L]]),
        "-",
        sprintf("%.1f", value_range[[2L]])
      )
    ) +
    theme_void(base_size = 6, base_family = figure_font_family) +
    theme(
      axis.text.x = element_text(size = 4.7, color = "#4B5563"),
      plot.title = element_text(
        size = 4.5,
        face = "bold",
        color = "#27313A",
        hjust = 0.5,
        margin = margin(b = 0)
      ),
      plot.margin = margin(0, 3, 0, 3)
    )
}

p_c_vivo_bar <- make_objective_colorbar(
  tsne_best_vivo$objective,
  "#2C7BB6",
  "#FDE725",
  "Vivo"
)
p_c_vitro_bar <- make_objective_colorbar(
  tsne_best_vitro$objective,
  "#1A9850",
  "#D73027",
  "Vitro"
)
p_c_objective_bars <- p_c_vivo_bar + p_c_vitro_bar +
  plot_layout(ncol = 2, widths = c(1, 1))

p_c_initial_key <- ggplot() +
  annotate(
    "point", x = 0.07, y = 0.5, shape = 16,
    size = 2.0, color = in_vivo_color
  ) +
  annotate(
    "text", x = 0.16, y = 0.5, label = "In vivo",
    hjust = 0, size = 1.85, color = "#27313A"
  ) +
  annotate(
    "point", x = 0.57, y = 0.5, shape = 17,
    size = 2.25, color = in_vitro_color
  ) +
  annotate(
    "text", x = 0.68, y = 0.5, label = "In vitro",
    hjust = 0, size = 1.85, color = "#27313A"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 2, 0, 2))

p_c_cluster_key <- ggplot(
  transform(
    cluster_labels,
    x = c(0.17, 0.50, 0.83),
    y = 0.5
  ),
  aes(x = x, y = y, label = label, color = cluster_id)
) +
  geom_label(
    fill = "white",
    text.colour = cluster_label_text_color,
    linewidth = cluster_label_border_linewidth,
    size = 1.75,
    fontface = cluster_label_fontface,
    label.padding = grid::unit(0.07, "lines"),
    label.r = grid::unit(0.05, "lines")
  ) +
  scale_color_manual(values = cluster_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 2, 0, 2))

p_c_shape_cluster_key <- wrap_plots(
  p_c_initial_key,
  p_c_cluster_key,
  ncol = 2,
  widths = c(0.52, 0.48)
)
p_landscape_key <- wrap_plots(
  p_c_objective_bars,
  p_c_shape_cluster_key,
  ncol = 1,
  heights = c(0.34, 0.66)
)
p_landscape <- p_c_scatter +
  labs(title = "B  Warm-start map") +
  theme(
    plot.title = element_text(
      family = figure_font_family,
      size = panel_title_size,
      face = "bold",
      color = "#111827",
      hjust = -0.090,
      margin = margin(b = 4)
    ),
    plot.title.position = "plot",
    plot.margin = margin(0.5, 0.5, 0.5, -5.1)
  )

# -------------------------------------------------------------------------
# Panel F: all-six soft-coupled parameter ratios on one log2 scale.
# -------------------------------------------------------------------------

ratio_rows <- list()
for (i in seq_len(nrow(selected))) {
  ratio <- read_winner_table(selected$fit_dir[[i]], "joint_soft_coupling.tsv")
  required_ratio_fields <- c(
    "parameter",
    "vivo_natural",
    "vitro_natural",
    "ratio_vivo_to_vitro",
    "feasible_at_solution",
    "projection_applied"
  )
  if (!all(required_ratio_fields %in% names(ratio))) {
    stop("joint_soft_coupling.tsv lacks required fields for ", selected$warmup_label[[i]])
  }
  ratio <- ratio[
    is.finite(ratio$ratio_vivo_to_vitro) &
      ratio$ratio_vivo_to_vitro > 0,
    ,
    drop = FALSE
  ]
  ratio$winner <- selected$winner_label[[i]]
  ratio$region <- selected$region[[i]]
  ratio$log2_ratio <- log2(ratio$ratio_vivo_to_vitro)
  ratio_rows[[i]] <- ratio
}
ratio_df <- bind_rows(ratio_rows)

if (nrow(ratio_df) != 84L ||
    length(unique(ratio_df$parameter)) != 14L ||
    any(table(ratio_df$parameter) != 6L)) {
  stop("Panel B requires 14 coupled parameters with six values each")
}
if (any(!ratio_df$feasible_at_solution) || any(ratio_df$projection_applied)) {
  stop("Unexpected infeasible or projected soft-coupling solution in Figure 5 universe")
}

coupled_parameters <- c(
  "lam_max", "p_mis_base", "p_wgd",
  "p_misseg", "k_o_mis",
  "O2_crit", "n_O",
  "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
parameter_groups <- parameter_meta_all[
  match(coupled_parameters, parameter_meta_all$parameter),
  c("parameter", "parameter_group", "parameter_order"),
  drop = FALSE
]
if (any(is.na(parameter_groups$parameter)) ||
    !identical(parameter_groups$parameter, coupled_parameters)) {
  stop("The shared parameter taxonomy does not preserve Figure 5F order")
}
mechanism_label_map <- c(
  "Division and genome change" = "Division +\ngenome change",
  "Stress-linked missegregation" = "Stress-linked\nmissegregation",
  "Oxygen response" = "Oxygen\nresponse",
  "Ploidy-linked growth and death" = "Ploidy-linked\ngrowth + death",
  "Post-missegregation survival" = "Post-MS\nsurvival"
)
parameter_groups$mechanism <- mechanism_label_map[parameter_groups$parameter_group]
mechanism_levels <- unname(
  mechanism_label_map[parameter_palette_meta$parameter_group]
)

ratio_df <- ratio_df %>%
  left_join(parameter_groups, by = "parameter")
if (any(is.na(ratio_df$mechanism)) || any(is.na(ratio_df$parameter_group))) {
  stop("Unmapped soft-coupling parameters in Panel B")
}
ratio_df$mechanism <- factor(ratio_df$mechanism, levels = mechanism_levels)
ratio_df$parameter <- factor(
  ratio_df$parameter,
  levels = rev(coupled_parameters)
)

ratio_summary <- ratio_df %>%
  group_by(mechanism, parameter) %>%
  summarise(
    ratio_min = min(log2_ratio),
    ratio_median = median(log2_ratio),
    ratio_max = max(log2_ratio),
    .groups = "drop"
  )

p_d_background <- parameter_groups %>%
  distinct(parameter_group, mechanism) %>%
  mutate(
    mechanism = factor(mechanism, levels = mechanism_levels),
    background_color = alpha(
      parameter_group_colors[parameter_group],
      0.11
    )
  )

ratio_breaks <- seq(-6, 12, by = 3)
ratio_break_labels <- vapply(ratio_breaks, function(x) {
  if (x == 0) {
    "1x"
  } else if (x < 0) {
    paste0("1/", format(2^abs(x), scientific = FALSE, trim = TRUE), "x")
  } else {
    paste0(format(2^x, scientific = FALSE, trim = TRUE), "x")
  }
}, character(1))

p_d_body_base <- ggplot(ratio_df, aes(y = parameter)) +
  geom_rect(
    data = p_d_background,
    aes(
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = parameter_group
    ),
    inherit.aes = FALSE,
    alpha = 0.11,
    color = NA
  ) +
  scale_fill_manual(
    values = parameter_group_colors,
    breaks = parameter_palette_meta$parameter_group,
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  annotate(
    "rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
    fill = alpha(in_vitro_color, 0.025)
  ) +
  annotate(
    "rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = alpha(in_vivo_color, 0.025)
  ) +
  geom_vline(
    xintercept = 0,
    color = "#111827",
    linetype = "22",
    linewidth = d_reference_linewidth
  ) +
  geom_segment(
    data = ratio_summary,
    aes(
      x = ratio_min,
      xend = ratio_max,
      y = parameter,
      yend = parameter
    ),
    inherit.aes = FALSE,
    color = "#7A8490",
    linewidth = d_range_linewidth,
    lineend = "round"
  ) +
  geom_point(
    aes(x = log2_ratio, fill = region),
    position = position_jitter(width = 0, height = 0.105, seed = 20260724),
    shape = 21,
    size = 1.95,
    color = ratio_marker_border_color,
    stroke = 0.38,
    alpha = 0.94
  ) +
  geom_point(
    data = ratio_summary,
    aes(x = ratio_median, y = parameter),
    inherit.aes = FALSE,
    shape = 23,
    fill = "#111827",
    color = "white",
    stroke = 0.35,
    size = 2.65
  ) +
  facet_grid(
    rows = vars(mechanism),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = region_colors,
    breaks = c("C01", "C02", "C03"),
    name = "Warm-start region"
  ) +
  scale_x_continuous(
    breaks = ratio_breaks,
    labels = ratio_break_labels,
    limits = c(-7, 13),
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "<- higher in vitro  in vivo / in vitro ratio  higher in vivo ->",
    y = NULL
  ) +
  theme_manuscript(base_size = 7) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "#F8FAFC", color = "#AEB6C0"),
    strip.text.y.left = element_text(
      family = figure_font_family,
      angle = 0,
      face = "bold",
      size = strip_text_size,
      color = "#374151"
    ),
    axis.text.y = element_text(
      family = figure_font_family,
      size = axis_text_size,
      face = "bold"
    ),
    axis.text.x = element_text(
      family = figure_font_family,
      size = axis_text_size
    ),
    axis.title.x = element_text(
      family = figure_font_family,
      size = axis_title_size,
      color = "transparent"
    ),
    panel.spacing.y = grid::unit(0.045, "cm"),
    legend.position = "top",
    legend.justification = "right",
    legend.box.margin = margin(f_legend_top_margin, 0, 0, 0),
    plot.margin = margin(2, 3, 2, 3)
  )

p_d_axis_band <- ggplot() +
  annotate(
    "text",
    x = -ratio_axis_side_offset,
    y = 0.5,
    label = "<- higher in vitro",
    hjust = 1,
    family = figure_font_family,
    size = axis_title_size / ggplot2::.pt,
    color = "#1F2933"
  ) +
  annotate(
    "text",
    x = 0,
    y = 0.5,
    label = "in vivo / in vitro ratio",
    hjust = 0.5,
    family = figure_font_family,
    size = axis_title_size / ggplot2::.pt,
    color = "#1F2933"
  ) +
  annotate(
    "text",
    x = ratio_axis_side_offset,
    y = 0.5,
    label = "higher in vivo ->",
    hjust = 0,
    family = figure_font_family,
    size = axis_title_size / ggplot2::.pt,
    color = "#1F2933"
  ) +
  coord_cartesian(
    xlim = c(-7, 13),
    ylim = c(0, 1),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 0, 0, 0))

p_d_body <- p_d_body_base +
  inset_element(
    p_d_axis_band,
    left = ratio_axis_inset_bounds[["left"]],
    right = ratio_axis_inset_bounds[["right"]],
    bottom = ratio_axis_inset_bounds[["bottom"]],
    top = ratio_axis_inset_bounds[["top"]],
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  )

p_d_header <- ggplot() +
  annotate(
    "text", x = -0.04, y = 0.80,
    label = "F  In-vivo / in-vitro parameter ratios",
    hjust = 0,
    size = panel_title_geom_size,
    family = figure_font_family,
    fontface = "bold",
    color = "#111827"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(1, 3, 0, 3))

p_d_note <- ggplot() +
  annotate(
    "text",
    x = -0.04,
    y = 0.68,
    label = "Six values/parameter; gray = range; black diamond = median.",
    hjust = 0,
    size = 2.0,
    color = "#4B5563"
  ) +
  annotate(
    "text",
    x = -0.04,
    y = 0.22,
    label = "Common log2 ratio scale for all 14 coupled parameters.",
    hjust = 0,
    size = 2.0,
    color = "#4B5563"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 3, 0, 3))

p_ratios <- wrap_plots(
  wrap_elements(full = p_d_header, clip = FALSE),
  wrap_elements(full = p_d_body),
  wrap_elements(full = p_d_note),
  ncol = 1,
  heights = unname(f_ratio_layout_heights)
)

# -------------------------------------------------------------------------
# Panel D: common-grid oxygen response functions.
# -------------------------------------------------------------------------

interpolate_winner_function <- function(
  fit_dir,
  context_dir,
  winner,
  region,
  value_col,
  common_grid
) {
  table_path <- file.path(
    "viz",
    context_dir,
    "functional_curve_oxygen_multi_ploidy.tsv"
  )
  curve <- read_winner_table(fit_dir, table_path)
  required_fields <- c("oxygen_pct", "cohort", value_col)
  if (!all(required_fields %in% names(curve))) {
    stop(table_path, " lacks required fields: ", paste(required_fields, collapse = ", "))
  }
  curve <- curve[
    curve$cohort %in% c("2N", "4N") &
      is.finite(curve$oxygen_pct) &
      is.finite(curve[[value_col]]),
    ,
    drop = FALSE
  ]
  context_label <- if (context_dir == "invivo") "In vivo" else "In vitro"
  out <- lapply(c("2N", "4N"), function(cohort_value) {
    cohort_curve <- curve[curve$cohort == cohort_value, , drop = FALSE]
    if (!nrow(cohort_curve)) {
      stop("Missing ", cohort_value, " curve for ", winner, " / ", context_label)
    }
    interpolated <- approx(
      x = cohort_curve$oxygen_pct,
      y = cohort_curve[[value_col]],
      xout = common_grid,
      ties = mean,
      rule = 2
    )
    data.frame(
      winner = winner,
      region = region,
      context = context_label,
      cohort = cohort_value,
      oxygen_pct = interpolated$x,
      value = interpolated$y
    )
  })
  bind_rows(out)
}

build_common_grid_function <- function(value_col, common_grid = seq(0, 5, by = 0.025)) {
  rows <- list()
  index <- 0L
  for (i in seq_len(nrow(selected))) {
    for (context_dir in c("invivo", "invitro")) {
      index <- index + 1L
      rows[[index]] <- interpolate_winner_function(
        fit_dir = selected$fit_dir[[i]],
        context_dir = context_dir,
        winner = selected$winner_label[[i]],
        region = selected$region[[i]],
        value_col = value_col,
        common_grid = common_grid
      )
    }
  }
  bind_rows(rows)
}

summarise_function <- function(df, x_col) {
  df %>%
    group_by(context, cohort, .data[[x_col]]) %>%
    summarise(
      value = median(value, na.rm = TRUE),
      value_min = min(value, na.rm = TRUE),
      value_max = max(value, na.rm = TRUE),
      n_solutions = n(),
      .groups = "drop"
    )
}

validate_six_function_solutions <- function(summary_df, label) {
  if (any(summary_df$n_solutions != 6L)) {
    stop(label, " does not contain six selected solutions at every grid point")
  }
}

proliferation_df <- build_common_grid_function("proliferation_rate")
missegregation_df <- build_common_grid_function("ms_rate")
proliferation_summary <- summarise_function(proliferation_df, "oxygen_pct")
missegregation_summary <- summarise_function(missegregation_df, "oxygen_pct")
validate_six_function_solutions(proliferation_summary, "Proliferation panel")
validate_six_function_solutions(missegregation_summary, "Missegregation panel")

function_limits <- function(raw_df, summary_df) {
  value_range <- range(
    c(raw_df$value, summary_df$value_min, summary_df$value_max),
    finite = TRUE
  )
  padding <- diff(value_range) * 0.025
  c(value_range[[1L]] - padding, value_range[[2L]] + padding)
}

make_function_cell <- function(
  raw_df,
  summary_df,
  cohort_value,
  y_limits,
  show_x_labels
) {
  raw_cell <- raw_df[raw_df$cohort == cohort_value, , drop = FALSE]
  summary_cell <- summary_df[
    summary_df$cohort == cohort_value,
    ,
    drop = FALSE
  ]
  ggplot(
    raw_cell,
    aes(
      x = oxygen_pct,
      y = value,
      color = context,
      linetype = context,
      group = interaction(context, winner)
    )
  ) +
    geom_line(
      linewidth = efg_individual_linewidth,
      alpha = 0.22,
      show.legend = FALSE
    ) +
    geom_line(
      data = summary_cell,
      aes(
        x = oxygen_pct,
        y = value,
        color = context,
        linetype = context,
        group = context
      ),
      inherit.aes = FALSE,
      linewidth = efg_median_linewidth,
      show.legend = FALSE
    ) +
    scale_color_manual(values = context_colors, guide = "none") +
    scale_linetype_manual(values = context_linetypes, guide = "none") +
    scale_x_continuous(
      breaks = c(0, 5),
      labels = c("0", "5"),
      limits = c(0, 5),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = pretty(y_limits, n = 3),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = 6) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(
        family = figure_font_family,
        size = axis_text_size
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0.8, 1.2, 0.8, 1.2)
    )
}

proliferation_limits <- function_limits(
  proliferation_df,
  proliferation_summary
)
missegregation_limits <- function_limits(
  missegregation_df,
  missegregation_summary
)

p_c_prolif_2n <- make_function_cell(
  proliferation_df,
  proliferation_summary,
  "2N",
  proliferation_limits,
  FALSE
)
p_c_ms_2n <- make_function_cell(
  missegregation_df,
  missegregation_summary,
  "2N",
  missegregation_limits,
  FALSE
)
p_c_prolif_4n <- make_function_cell(
  proliferation_df,
  proliferation_summary,
  "4N",
  proliferation_limits,
  TRUE
)
p_c_ms_4n <- make_function_cell(
  missegregation_df,
  missegregation_summary,
  "4N",
  missegregation_limits,
  TRUE
)

p_c_title <- make_text_band(
  "D  Oxygen-dependent fitted functions",
  size = panel_title_size,
  hjust = 0,
  x = -0.015
)
p_c_header_block <- p_c_title
p_c_shared_x_title <- ggplot() +
  annotate(
    "segment",
    x = d_shared_x_positions,
    xend = d_shared_x_positions,
    y = 1.00,
    yend = 0.84,
    linewidth = 0.22,
    color = "#27313A"
  ) +
  annotate(
    "text",
    x = d_shared_x_positions,
    y = 0.69,
    label = c("0", "5", "0", "5"),
    size = axis_text_size / ggplot2::.pt,
    family = figure_font_family,
    color = "#27313A"
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 0.16,
    label = "Oxygen (%)",
    size = axis_title_size / ggplot2::.pt,
    family = figure_font_family,
    color = "#1F2933"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 0, 0, 0))
p_c_grid_design <- c(
  area(t = 1, l = 2, b = 1, r = 2),
  area(t = 1, l = 4, b = 1, r = 4),
  area(t = 2, l = 1, b = 3, r = 1),
  area(t = 2, l = 2, b = 2, r = 2),
  area(t = 2, l = 3, b = 3, r = 3),
  area(t = 2, l = 4, b = 2, r = 4),
  area(t = 2, l = 5, b = 2, r = 5),
  area(t = 3, l = 2, b = 3, r = 2),
  area(t = 3, l = 4, b = 3, r = 4),
  area(t = 3, l = 5, b = 3, r = 5),
  area(t = 4, l = 2, b = 4, r = 4)
)
p_c_plot_matrix <- wrap_plots(
  make_strip_band("Proliferation", size = strip_text_size) +
    theme(plot.margin = margin(0, 0, 0, 0)),
  make_strip_band("MS", size = strip_text_size) +
    theme(plot.margin = margin(0, 0, 0, 0)),
  make_vertical_axis_title("Proliferation~rate~(day^{-1})"),
  p_c_prolif_2n,
  make_vertical_axis_title("MS~rate~(chr^{-1})"),
  p_c_ms_2n,
  make_vertical_strip_band("2N", size = 6.2),
  p_c_prolif_4n,
  p_c_ms_4n,
  make_vertical_strip_band("4N", size = 6.2),
  p_c_shared_x_title,
  design = p_c_grid_design,
  widths = c(0.015, 0.465, 0.015, 0.465, 0.04),
  heights = unname(c_grid_heights)
)
p_c <- wrap_plots(
  wrap_elements(full = p_c_header_block),
  wrap_elements(full = p_c_plot_matrix),
  ncol = 1,
  heights = c(0.10, 0.90)
)

# -------------------------------------------------------------------------
# Panel C: post-missegregation survival across the identical six winners.
# -------------------------------------------------------------------------

survival_rows <- list()
index <- 0L
for (i in seq_len(nrow(selected))) {
  for (context_dir in c("invivo", "invitro")) {
    index <- index + 1L
    survival <- read_winner_table(
      selected$fit_dir[[i]],
      file.path("viz", context_dir, "functional_curve_ploidy.tsv")
    )
    required_survival_fields <- c("N", "endpoint_value", "viability_after_ms")
    if (!all(required_survival_fields %in% names(survival))) {
      stop("functional_curve_ploidy.tsv lacks required survival fields")
    }
    survival <- survival[
      is.finite(survival$endpoint_value) &
        is.finite(survival$viability_after_ms),
      ,
      drop = FALSE
    ]
    survival$winner <- selected$winner_label[[i]]
    survival$region <- selected$region[[i]]
    survival$context <- if (context_dir == "invivo") "In vivo" else "In vitro"
    survival$cohort <- "All chromosome states"
    survival$value <- survival$viability_after_ms
    survival_rows[[index]] <- survival[, c(
      "winner",
      "region",
      "context",
      "cohort",
      "endpoint_value",
      "value"
    )]
  }
}
survival_df <- bind_rows(survival_rows)
survival_summary <- summarise_function(survival_df, "endpoint_value")
validate_six_function_solutions(survival_summary, "Post-MS survival panel")

p_d_body <- ggplot(
  survival_df,
  aes(
    x = endpoint_value,
    y = value,
    color = context,
    linetype = context,
    group = interaction(context, winner)
  )
) +
  geom_line(
    linewidth = efg_individual_linewidth,
    alpha = 0.22,
    show.legend = FALSE
  ) +
  geom_line(
    data = survival_summary,
    aes(
      x = endpoint_value,
      y = value,
      color = context,
      linetype = context,
      group = context
    ),
    inherit.aes = FALSE,
    linewidth = efg_median_linewidth
  ) +
  scale_color_manual(
    values = context_colors,
    breaks = c("In vivo", "In vitro"),
    name = "Context"
  ) +
  scale_linetype_manual(
    values = context_linetypes,
    breaks = c("In vivo", "In vitro"),
    name = "Context"
  ) +
  scale_x_continuous(
    breaks = c(22, 88),
    limits = c(22, 154),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "C  Survival",
    subtitle = NULL,
    x = "Chromosome N",
    y = "Survival"
  ) +
  theme_manuscript(base_size = 7) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      family = figure_font_family,
      size = panel_title_size,
      face = "bold",
      color = "#111827",
      margin = margin(b = 4)
    ),
    aspect.ratio = 1,
    axis.title = element_text(
      family = figure_font_family,
      size = axis_title_size
    ),
    axis.text = element_text(
      family = figure_font_family,
      size = axis_text_size
    ),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5)
  )
p_d <- p_d_body

context_key <- ggplot() +
  annotate(
    "segment", x = 0.18, xend = 0.30, y = 0.5, yend = 0.5,
    color = in_vivo_color, linewidth = efg_median_linewidth
  ) +
  annotate(
    "text", x = 0.32, y = 0.5, label = "In vivo",
    hjust = 0, size = 2.0, color = "#27313A"
  ) +
  annotate(
    "segment", x = 0.57, xend = 0.69, y = 0.5, yend = 0.5,
    color = in_vitro_color, linewidth = efg_median_linewidth
  ) +
  annotate(
    "text", x = 0.71, y = 0.5, label = "In vitro",
    hjust = 0, size = 2.0, color = "#27313A"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 0, 0, 0))

# -------------------------------------------------------------------------
# Save individual panels and the assembled Figure 5.
# -------------------------------------------------------------------------

stale_panel_stubs <- file.path(
  panel_root,
  c(
    "figure5b_fit_adequacy",
    "figure5c_o2_functions",
    "figure5d_survival",
    "figure5e_warm_start_map"
  )
)
invisible(lapply(
  c(paste0(stale_panel_stubs, ".png"), paste0(stale_panel_stubs, ".pdf")),
  function(path) {
    if (file.exists(path)) {
      unlink(path)
    }
  }
))

save_plot_pair(
  p_a,
  file.path(panel_root, "figure5a_joint_fitting_workflow"),
  width = 7.1,
  height = 1.6,
  dpi = 300
)
save_plot_pair(
  p_b,
  file.path(panel_root, "figure5e_fit_adequacy"),
  width = 4.0,
  height = 4.0,
  dpi = 300
)
save_plot_pair(
  p_c,
  file.path(panel_root, "figure5d_o2_functions"),
  width = 3.4,
  height = 2.35,
  dpi = 300
)
save_plot_pair(
  p_d,
  file.path(panel_root, "figure5c_survival"),
  width = 2.0,
  height = 2.0,
  dpi = 300
)
save_plot_pair(
  p_landscape,
  file.path(panel_root, "figure5b_warm_start_map"),
  width = 2.0,
  height = 2.0,
  dpi = 300
)
save_plot_pair(
  p_ratios,
  file.path(panel_root, "figure5f_parameter_ratios"),
  width = 7.1,
  height = 3.45 * f_height_scale,
  dpi = 300
)

compact_function_row <- wrap_plots(
  p_landscape,
  p_d,
  ncol = 2,
  widths = c(0.50, 0.50)
)
compact_right_legends <- wrap_plots(
  wrap_elements(full = p_landscape_key),
  wrap_elements(full = context_key),
  ncol = 2,
  widths = c(0.50, 0.50)
) &
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )
compact_right_legends_shifted <- plot_spacer() +
  inset_element(
    compact_right_legends,
    left = 0,
    right = 1,
    bottom = de_legend_vertical_shift,
    top = 1 + de_legend_vertical_shift,
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  )
compact_right_block <- wrap_plots(
  compact_function_row,
  compact_right_legends_shifted,
  p_c_header_block,
  p_c_plot_matrix,
  ncol = 1,
  heights = c(
    (
      aligned_row_heights[["first_plot"]] +
        aligned_row_heights[["second_plot"]]
    ) * cde_content_heights[["panels_d_e"]],
    aligned_row_heights[["footer"]],
    aligned_row_heights[["header"]],
    (
      aligned_row_heights[["first_plot"]] +
        aligned_row_heights[["second_plot"]]
    ) * cde_content_heights[["panel_c"]]
  )
)
compact_bcd_scaled <- plot_spacer() +
  inset_element(
    compact_right_block,
    left = 0,
    right = bcd_visual_scale[["width"]],
    bottom = 1 - bcd_visual_scale[["height"]],
    top = 1,
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  )

fit_and_functions_row <- wrap_plots(
  wrap_elements(full = compact_bcd_scaled, clip = FALSE),
  wrap_elements(full = p_b, clip = FALSE),
  ncol = 2,
  widths = unname(row2_block_widths)
)

figure5 <- (
    wrap_elements(full = p_a) /
    wrap_elements(full = fit_and_functions_row) /
    wrap_elements(full = p_ratios) +
    plot_layout(
      heights = c(
        0.67 * a_height_scale,
        1.55,
        1.20 * f_height_scale
      )
    )
) +
  plot_annotation(
    title = "Joint-fit design, adequacy, and context-specific fitted functions",
    subtitle = paste0(
      "Six landscape-informed joint-pair winners. Direct fits are in-sample; ",
      "across-winner ranges and traces show solution sensitivity, not confidence intervals."
    ),
    caption = paste0(
      "C-D: thin curves = six fits; thick = pointwise median. ",
      "Necrosis omitted: predictions unavailable.\n",
      "B-F use the same six winners; D is display-interpolated ",
      "to a declared 0-5% O2 grid."
    ),
    theme = theme(
      text = element_text(family = figure_font_family),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        family = figure_font_family,
        size = 11.5,
        face = "bold",
        color = "#111827",
        margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        family = figure_font_family,
        size = 7.3,
        color = "#4B5563",
        margin = margin(b = 2)
      ),
      plot.caption = element_text(
        family = figure_font_family,
        size = 7,
        color = "#5F6B76",
        hjust = 0,
        margin = margin(t = 2)
      ),
      plot.margin = margin(4, 4, 3, 4)
    )
  )

output_stub <- file.path(figure_root, "assembled_fig5")
save_plot_pair(
  figure5,
  output_stub,
  width = 7.1,
  height = assembled_height_inches,
  dpi = 300
)
for (extension in c("png", "pdf")) {
  source_path <- normalizePath(
    paste0(output_stub, ".", extension), mustWork = TRUE
  )
  destination_path <- normalizePath(
    file.path(revised_root, paste0("assembled_fig5.", extension)),
    mustWork = FALSE
  )
  if (!identical(source_path, destination_path)) {
    file.copy(source_path, destination_path, overwrite = TRUE)
  }
}

provenance <- data.frame(
  role = c(
    "selected joint winners",
    "pooled t-SNE full coordinates",
    "pooled t-SNE best-fit cluster coordinates",
    "shared parameter-function taxonomy",
    "shared parameter-function palette"
  ),
  path = normalizePath(
    c(
      selection_path,
      tsne_full_path,
      tsne_best_path,
      parameter_meta_path,
      parameter_palette_path
    ),
    mustWork = TRUE
  ),
  md5 = unname(tools::md5sum(
    c(
      selection_path,
      tsne_full_path,
      tsne_best_path,
      parameter_meta_path,
      parameter_palette_path
    )
  )),
  stringsAsFactors = FALSE
)
write.table(
  provenance,
  file.path(figure5_data_dir, "figure5_source_file_provenance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

validation <- data.frame(
  check = c(
    "selected_winners",
    "ratio_rows",
    "ratio_parameters",
    "pooled_tsne_rows",
    "in_vivo_best_tsne_rows",
    "in_vitro_best_tsne_rows",
    "displayed_in_vivo_cluster_hulls",
    "displayed_in_vitro_cluster_hulls",
    "parameter_function_groups",
    "function_titles_contain_newline",
    "figure_font_family",
    "panel_title_size",
    "axis_title_size",
    "axis_text_size",
    "strip_text_size",
    "panel_order",
    "row2_outer_panel_order",
    "compact_top_row_order",
    "compact_second_row_panel",
    "assembled_height_inches",
    "panel_a_height_scale",
    "panel_a_min_horizontal_arrow_shaft",
    "panel_a_output_branch_length",
    "panel_a_main_arrowhead_length",
    "panel_a_branch_arrowhead_length",
    "panel_a_welsch_label_below_icon",
    "panel_a_welsch_label_center",
    "panel_a_paired_parameter_box_width",
    "panel_a_welsch_label_size",
    "panel_a_welsch_detail_size",
    "panel_f_height_scale",
    "panel_b_coordinate_aspect",
    "panel_b_content_grid",
    "panel_c_equal_fill_target_aspect",
    "panel_c_ploidy_pagination",
    "panel_c_ploidy_strip_side",
    "panel_c_header_labels",
    "panel_c_headers_bound_to_plot_columns",
    "panel_c_header_strips_fill_column_width",
    "panel_c_shared_y_titles_centered",
    "panel_c_shared_x_title",
    "panel_c_ploidy_strip_angle",
    "panel_c_shared_bottom_axis",
    "panel_d_e_coordinate_aspect",
    "panel_d_e_common_header",
    "panel_d_e_top_aligned",
    "panel_bcd_e_equal_outer_size",
    "panel_bcd_visible_scale",
    "panel_b_c_e_titles_top_aligned",
    "panel_a_b_d_f_titles_left_aligned",
    "panel_d_shared_x_positions",
    "panel_d_strip_height",
    "panel_b_c_title_bottom_margin",
    "panel_e_layout_heights",
    "panel_e_note_y",
    "panel_f_layout_heights",
    "panel_f_ratio_axis_label_order",
    "panel_f_ratio_axis_center",
    "panel_f_ratio_axis_inset_bounds",
    "panel_f_ratio_axis_side_offset",
    "panel_f_legend_top_margin",
    "row2_left_right_widths",
    "panel_b_c_widths",
    "cluster_label_point_overlaps",
    "cluster_labels_outside_scatter",
    "panel_e_grid_removed",
    "panel_e_border_removed",
    "panel_e_cell_fill_fraction",
    "panel_e_legend_order",
    "objective_bar_height",
    "compact_left_panels",
    "aligned_row_heights",
    "panel_c_grid_heights",
    "panel_cde_content_heights",
    "panel_e_coordinate_example",
    "panel_de_legend_vertical_shift",
    "in_vivo_color",
    "in_vitro_color",
    "cluster_color_C01",
    "cluster_color_C02",
    "cluster_color_C03",
    "panel_b_cluster_hull_linetype",
    "panel_b_cluster_label_border_colors",
    "panel_b_cluster_label_text_color",
    "panel_b_cluster_label_fontface",
    "panel_b_cluster_label_text_outline",
    "panel_b_cluster_label_border_linewidth",
    "panel_b_tsne_coordinate_label_scale",
    "panel_b_tsne_coordinate_label_size_pt",
    "best_seed_border_color",
    "ratio_marker_border_color",
    "d_reference_linewidth",
    "d_range_linewidth",
    "efg_individual_linewidth",
    "efg_median_linewidth",
    "in_vivo_linetype",
    "in_vitro_linetype",
    "assembled_png",
    "assembled_pdf"
  ),
  observed = c(
    nrow(selected),
    nrow(ratio_df),
    length(unique(ratio_df$parameter)),
    nrow(tsne_full),
    nrow(tsne_best_vivo),
    nrow(tsne_best_vitro),
    length(unique(cluster_hulls$cluster_id)),
    sum(grepl("^vt_", unique(cluster_hulls$cluster_id))),
    length(unique(ratio_df$parameter_group)),
    any(grepl(
      "\n",
      c(
        "D  Oxygen-dependent fitted functions",
        "C  Post-MS survival"
      ),
      fixed = TRUE
    )),
    figure_font_family,
    panel_title_size,
    axis_title_size,
    axis_text_size,
    strip_text_size,
    "A,B,C,D,E,F",
    "B-C-D,E",
    "B,C",
    "D",
    assembled_height_inches,
    a_height_scale,
    min(c(
      workflow_lane_arrows$xend - workflow_lane_arrows$x,
      workflow_core_arrows$xend - workflow_core_arrows$x
    )),
    unique(workflow_output_branches$y - workflow_output_branches$yend),
    workflow_main_arrowhead_length,
    workflow_branch_arrowhead_length,
    workflow_welsch_layout[["label_y"]] <
      workflow_welsch_layout[["bar_ymin"]],
    paste(
      workflow_welsch_layout[["label_x"]],
      workflow_welsch_layout[["label_y"]],
      sep = ","
    ),
    workflow_core_boxes$xmax[[2L]] - workflow_core_boxes$xmin[[2L]],
    workflow_welsch_label_size,
    workflow_welsch_detail_size,
    f_height_scale,
    1,
    "2x2",
    "2:1",
    "2N-over-4N",
    "right",
    "Proliferation,MS",
    TRUE,
    TRUE,
    TRUE,
    "Oxygen (%)",
    270,
    TRUE,
    "1,1",
    TRUE,
    TRUE,
    TRUE,
    paste(sprintf("%.2f", bcd_visual_scale), collapse = ","),
    TRUE,
    TRUE,
    paste(d_shared_x_positions, collapse = ","),
    c_grid_heights[["mechanism_strip"]],
    4,
    paste(sprintf("%.2f", e_layout_heights), collapse = ","),
    0.90,
    paste(sprintf("%.3f", f_ratio_layout_heights), collapse = ","),
    "<- higher in vitro | in vivo / in vitro ratio | higher in vivo ->",
    0,
    paste(sprintf("%.3f", ratio_axis_inset_bounds), collapse = ","),
    ratio_axis_side_offset,
    f_legend_top_margin,
    paste(sprintf("%.2f", row2_block_widths), collapse = ","),
    "0.50,0.50",
    paste(cluster_labels$point_overlap_count, collapse = ","),
    TRUE,
    TRUE,
    TRUE,
    1.00,
    "objective-above-shape-and-cluster",
    0.18,
    3,
    paste(sprintf("%.2f", aligned_row_heights), collapse = ","),
    paste(sprintf("%.2f", c_grid_heights), collapse = ","),
    paste(sprintf("%.2f", cde_content_heights), collapse = ","),
    TRUE,
    de_legend_vertical_shift,
    context_colors[["In vivo"]],
    context_colors[["In vitro"]],
    region_colors[["C01"]],
    region_colors[["C02"]],
    region_colors[["C03"]],
    cluster_hull_linetype,
    paste(unname(cluster_colors), collapse = ","),
    cluster_label_text_color,
    cluster_label_fontface,
    FALSE,
    cluster_label_border_linewidth,
    tsne_coordinate_label_scale,
    axis_title_size * tsne_coordinate_label_scale,
    best_seed_border_color,
    ratio_marker_border_color,
    d_reference_linewidth,
    d_range_linewidth,
    efg_individual_linewidth,
    efg_median_linewidth,
    context_linetypes[["In vivo"]],
    context_linetypes[["In vitro"]],
    file.exists(paste0(output_stub, ".png")),
    file.exists(paste0(output_stub, ".pdf"))
  ),
  expected = c(
    "6", "84", "14", "228000", "500", "500",
    "3", "0", "5", "FALSE",
    "Helvetica", "7.5", "6.1", "5.8", "5.8",
    "A,B,C,D,E,F", "B-C-D,E", "B,C", "D",
    "9.44", "0.8",
    "4", "2.25", "0.072", "0.06", "TRUE", "64.75,10.95",
    "14.5", "1.75", "1.65",
    "1.2", "1",
    "2x2", "2:1", "2N-over-4N", "right",
    "Proliferation,MS", "TRUE", "TRUE", "TRUE", "Oxygen (%)", "270", "TRUE",
    "1,1", "TRUE", "TRUE", "TRUE",
    "0.99,0.93", "TRUE", "TRUE", "0,0.46,0.54,1",
    "0.055", "4", "0.08,0.87,0.05", "0.9",
    "0.010,0.910,0.080",
    "<- higher in vitro | in vivo / in vitro ratio | higher in vivo ->",
    "0", "0.208,0.994,0.000,0.055", "2.1",
    "-5",
    "1.00,1.00", "0.50,0.50", "0,0,0", "TRUE", "TRUE",
    "TRUE", "1", "objective-above-shape-and-cluster", "0.18",
    "3", "0.08,0.41,0.41,0.10",
    "0.06,0.42,0.42,0.10", "0.56,0.44", "TRUE",
    "0.42",
    "#0072B2", "#CC79A7",
    "#FFD700", "#6A3D9A", "#39FF14",
    "solid", "#FFD700,#6A3D9A,#39FF14",
    "#000000", "plain", "FALSE",
    "0.33", "0.5", "3.05",
    "#FFFFFF", "#FFFFFF",
    "0.275", "0.625", "0.19", "0.56",
    "solid", "solid",
    "TRUE", "TRUE"
  ),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(figure5_data_dir, "figure5_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Figure 5 generation complete.\n")
cat("Selected winners:", nrow(selected), "\n")
cat(
  "Total objective range:",
  sprintf("%.4f", min(selected$objective)),
  "to",
  sprintf("%.4f", max(selected$objective)),
  "\n"
)
cat("Panel F ratio rows:", nrow(ratio_df), "\n")
cat("Assembled PNG:", paste0(output_stub, ".png"), "\n")
cat("Assembled PDF:", paste0(output_stub, ".pdf"), "\n")

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Figure5 <- function() {
  data_dir <- file.path(DATA_ROOT, "Figure5")
  tsne_full <- file.path(
    data_dir,
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_full_coordinates.csv"
  )
  tsne_best <- file.path(
    data_dir,
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
  )
  require_files(c(
    file.path(data_dir, "figure5_frozen_inputs", "selected_results.tsv"),
    tsne_full, tsne_best,
    file.path(data_dir, "parameter_function_groups.tsv"),
    file.path(data_dir, "parameter_function_group_palette.tsv")
  ), "Figure 5 intermediate")
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_Figure5.R"), mustWork = TRUE),
    env = c(
      "FIGURE5_DRAW_WORKER=1",
      paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("FIGURE5_DATA_DIR=", data_dir),
      paste0("FIGURE5_FIGURE_DIR=", OUTPUT_ROOT),
      paste0("FIGURE5_PANEL_DIR=", file.path(data_dir, "panels")),
      paste0("JOINT_TSNE_FULL_COORDINATES=", tsne_full),
      paste0("JOINT_TSNE_BEST_COORDINATES=", tsne_best)
    )
  )
  require_files(
    file.path(OUTPUT_ROOT, c("assembled_fig5.png", "assembled_fig5.pdf")),
    "Figure 5 output"
  )
  invisible(file.path(OUTPUT_ROOT, "assembled_fig5.png"))
}

if (sys.nframe() == 0L) draw_Figure5()

}
