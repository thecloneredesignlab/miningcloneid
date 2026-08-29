#!/usr/bin/env Rscript

if (identical(Sys.getenv("FIGURE5_DRAW_WORKER"), "1")) {

# Standalone iteration4 Figure 5 generator.
#
# Scope:
#   A. Joint-fitting workflow.
#   B. Direct in-sample observed-predicted model-performance summaries.
#   C. Pooled separate-fit t-SNE landscape used to select joint warm starts.
#   D. Reconstructed in-vivo and in-vitro differential-evolution initial
#      populations and optimizer-endpoint distributions, separated by
#      the manifest-declared primary C families.
#   E. Proliferation and stress-linked missegregation functions.
#   F. Post-missegregation survival functions across selected winners.
#
# This script is self-contained: it does not source functions from another
# checkout. It reads only frozen inputs and derived optimizer-ensemble tables
# stored in iteration4. It does not refit the model or treat numerical starts as
# biological replication, posterior draws, confidence intervals, or proof of
# structural parameter identifiability.

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
source(file.path(script_dir, "util", "runtime", "workspace_paths.R"))
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
optimizer_solution_path <- file.path(
  figure5_data_dir, "figure5f_optimizer_solutions.tsv"
)
pair_summary_path <- file.path(
  figure5_data_dir, "figure5f_pair_summary.tsv"
)
family_summary_path <- file.path(
  figure5_data_dir, "figure5f_family_summary.tsv"
)
cross_family_summary_path <- file.path(
  figure5_data_dir, "figure5f_cross_family_summary.tsv"
)
family_density_path <- file.path(
  figure5_data_dir, "figure5f_family_density.tsv"
)
sensitivity_validation_path <- file.path(
  figure5_data_dir, "figure5f_sensitivity_validation.tsv"
)
ratio_density_path <- file.path(
  figure5_data_dir, "figure5f_prior_optimizer_density.tsv"
)
ratio_summary_path <- file.path(
  figure5_data_dir, "figure5f_prior_optimizer_summary.tsv"
)
ratio_cross_family_path <- file.path(
  figure5_data_dir, "figure5f_prior_optimizer_cross_family.tsv"
)
ratio_readiness_path <- file.path(
  figure5_data_dir, "figure5f_prior_optimizer_readiness.tsv"
)
context_density_path <- file.path(
  figure5_data_dir, "figure5f_context_initial_optimizer_density.tsv"
)
context_summary_path <- file.path(
  figure5_data_dir, "figure5f_context_initial_optimizer_summary.tsv"
)
chart_contract_path <- file.path(
  figure5_data_dir, "figure5f_chart_contract.md"
)

missing_configuration <- c(
  selection_path,
  parameter_meta_path,
  parameter_palette_path,
  optimizer_solution_path,
  pair_summary_path,
  family_summary_path,
  cross_family_summary_path,
  family_density_path,
  sensitivity_validation_path,
  ratio_density_path,
  ratio_summary_path,
  ratio_cross_family_path,
  ratio_readiness_path,
  context_density_path,
  context_summary_path,
  chart_contract_path
)[!file.exists(c(
  selection_path,
  parameter_meta_path,
  parameter_palette_path,
  optimizer_solution_path,
  pair_summary_path,
  family_summary_path,
  cross_family_summary_path,
  family_density_path,
  sensitivity_validation_path,
  ratio_density_path,
  ratio_summary_path,
  ratio_cross_family_path,
  ratio_readiness_path,
  context_density_path,
  context_summary_path,
  chart_contract_path
))]
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
    unset = file.path(figure5_data_dir, "panels")
  ),
  mustWork = FALSE
)
revised_root <- file.path(workspace_root, "manuscript", "Figures")
allowed_output_dirs <- c(figure_root, panel_root, revised_root)
invisible(lapply(
  allowed_output_dirs,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

selected_all <- read.delim(
  selection_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
selected_all <- selected_all[
  selected_all$record_type == "joint_pair_best",
  ,
  drop = FALSE
]
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
optimizer_solutions <- read.delim(
  optimizer_solution_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pair_solution_summary <- read.delim(
  pair_summary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
family_solution_summary <- read.delim(
  family_summary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cross_family_summary <- read.delim(
  cross_family_summary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
family_density <- read.delim(
  family_density_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sensitivity_validation <- read.delim(
  sensitivity_validation_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ratio_density <- read.delim(
  ratio_density_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ratio_summary <- read.delim(
  ratio_summary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ratio_cross_family <- read.delim(
  ratio_cross_family_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ratio_readiness <- read.delim(
  ratio_readiness_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
context_density <- read.delim(
  context_density_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
context_summary <- read.delim(
  context_summary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
primary_family_order <- JOINT_FAMILY_LEVELS
primary_family_count <- length(primary_family_order)
selected_pair_map <- unique(
  context_summary[, c("family", "warmup_label"), drop = FALSE]
)
selected_pair_map <- selected_pair_map[
  match(primary_family_order, selected_pair_map$family),
  ,
  drop = FALSE
]
if (nrow(selected_pair_map) != primary_family_count ||
    anyNA(selected_pair_map$family) ||
    !identical(selected_pair_map$family, primary_family_order)) {
  stop("Figure 5 requires one retained joint pair for each primary C family")
}
selected <- selected_all[
  match(selected_pair_map$warmup_label, selected_all$warmup_label),
  ,
  drop = FALSE
]
if (anyNA(selected$warmup_label)) {
  stop("A retained primary-family pair is absent from selected_results.tsv")
}
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
if (nrow(optimizer_solutions) != 500L * 14L * primary_family_count ||
    nrow(pair_solution_summary) != 14L * primary_family_count ||
    nrow(family_solution_summary) != 14L * primary_family_count ||
    nrow(cross_family_summary) != 14L ||
    nrow(family_density) != 14L * primary_family_count * 401L ||
    nrow(sensitivity_validation) != 42L ||
    length(unique(optimizer_solutions$family)) != primary_family_count ||
    length(unique(optimizer_solutions$pair_id)) != primary_family_count ||
    any(!optimizer_solutions$feasible_at_solution) ||
    any(optimizer_solutions$projection_applied)) {
  stop("Figure 5 optimizer-ensemble tables fail the declared data contract")
}
if (nrow(ratio_density) != 14L * primary_family_count * 2L * 401L ||
    nrow(ratio_summary) != 14L * primary_family_count ||
    nrow(ratio_cross_family) != 14L ||
    !nrow(ratio_readiness) || any(!ratio_readiness$passed)) {
  stop("Figure 5 DE-initial/optimizer tables fail the declared readiness contract")
}
if (nrow(context_density) != 14L * primary_family_count * 2L * 2L * 401L ||
    nrow(context_summary) != 14L * primary_family_count * 2L ||
    !setequal(unique(context_density$context), c("in vivo", "in vitro")) ||
    !setequal(
      unique(context_density$distribution_role),
      c("de_initial_population", "optimizer_endpoints")
    )) {
  stop("Figure 5 direct-context distributions fail the declared data contract")
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
if (nrow(selected) != primary_family_count) {
  stop("Figure 5 requires exactly one retained pair per primary family; found ", nrow(selected))
}
if (length(unique(selected$invitro_seed)) != 1L ||
    unique(selected$invitro_seed) != INVITRO_VISUALIZATION_SEED) {
  stop(
    "The approved Figure 5 universe must use in-vitro anchor ",
    INVITRO_VISUALIZATION_SEED
  )
}

extract_region <- function(warmup_label) {
  hit <- regmatches(warmup_label, regexpr("C[0-9]{2}", warmup_label))
  ifelse(nzchar(hit), hit, "unclassified")
}

selected$region <- extract_region(selected$warmup_label)
if (!identical(selected$region, primary_family_order)) {
  stop("Retained Figure 5 pairs are not ordered one per declared C family")
}
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

# The red edge bands in panel D are defined from the complete joint-union
# optimizer bounds, not from quantiles of either the reconstructed DE initial
# populations or the optimizer endpoints. The frozen winner tables repeat the
# same joint-union bounds for all retained C families; reading all six here
# makes that equality an explicit release check.
joint_bound_rows <- bind_rows(lapply(seq_len(nrow(selected)), function(i) {
  soft <- read_winner_table(selected$fit_dir[[i]], "joint_soft_coupling.tsv")
  required_bound_fields <- c(
    "parameter",
    "joint_union_lower_transformed", "joint_union_upper_transformed",
    "joint_union_lower_bound", "joint_union_upper_bound"
  )
  if (!all(required_bound_fields %in% names(soft))) {
    stop(
      "joint_soft_coupling.tsv lacks joint-union bounds for ",
      selected$warmup_label[[i]]
    )
  }
  out <- soft[, required_bound_fields, drop = FALSE]
  out$family <- selected$region[[i]]
  out$warmup_label <- selected$warmup_label[[i]]
  out
}))
joint_bound_numeric_fields <- c(
  "joint_union_lower_transformed", "joint_union_upper_transformed",
  "joint_union_lower_bound", "joint_union_upper_bound"
)
joint_bound_consistent <- vapply(
  split(joint_bound_rows, joint_bound_rows$parameter),
  function(group) {
    all(vapply(
      joint_bound_numeric_fields,
      function(field) length(unique(signif(group[[field]], 14L))) == 1L,
      logical(1)
    ))
  },
  logical(1)
)
if (nrow(joint_bound_rows) != 14L * primary_family_count ||
    length(joint_bound_consistent) != 14L ||
    any(!joint_bound_consistent)) {
  stop("Panel D joint-union bounds are incomplete or differ across C families")
}
joint_parameter_bounds <- joint_bound_rows[
  !duplicated(joint_bound_rows$parameter),
  c("parameter", joint_bound_numeric_fields),
  drop = FALSE
]
joint_parameter_bounds$edge_fraction <- 0.05
joint_parameter_bounds$lower_edge_end_transformed <- with(
  joint_parameter_bounds,
  joint_union_lower_transformed + edge_fraction *
    (joint_union_upper_transformed - joint_union_lower_transformed)
)
joint_parameter_bounds$upper_edge_start_transformed <- with(
  joint_parameter_bounds,
  joint_union_lower_transformed + (1 - edge_fraction) *
    (joint_union_upper_transformed - joint_union_lower_transformed)
)
joint_parameter_bounds$transformation <- parameter_meta_all$transformation[
  match(joint_parameter_bounds$parameter, parameter_meta_all$parameter)
]
if (anyNA(joint_parameter_bounds$transformation)) {
  stop("Panel D joint bounds could not be matched to parameter transformations")
}
inverse_panel_scale <- function(values, transformation) {
  if (identical(transformation, "log10")) 10^values else values
}
joint_parameter_bounds$lower_edge_end_natural <- mapply(
  inverse_panel_scale,
  joint_parameter_bounds$lower_edge_end_transformed,
  joint_parameter_bounds$transformation
)
joint_parameter_bounds$upper_edge_start_natural <- mapply(
  inverse_panel_scale,
  joint_parameter_bounds$upper_edge_start_transformed,
  joint_parameter_bounds$transformation
)
joint_parameter_bounds$band_definition <-
  "outer 5% of complete joint-union interval on optimizer/display scale"
joint_bound_band_path <- file.path(
  figure5_data_dir, "figure5d_joint_bound_edge_bands.tsv"
)
write.table(
  joint_parameter_bounds,
  joint_bound_band_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)

context_colors <- c(
  "In vivo" = "#0072B2",
  "In vitro" = "#CC79A7"
)
region_colors <- c(
  "C01" = "#C99700",
  "C02" = "#6A3D9A",
  "C03" = "#006D2C",
  "C04" = "#1F78B4",
  "C05" = "#D95F02",
  "C06" = "#009E73"
)
cluster_hull_linetype <- "solid"
cluster_label_text_color <- "#000000"
cluster_label_fontface <- "plain"
cluster_label_border_linewidth <- 0.33
tsne_coordinate_label_scale <- 0.65
in_vivo_color <- unname(context_colors[["In vivo"]])
in_vitro_color <- unname(context_colors[["In vitro"]])
best_seed_border_color <- "#FFFFFF"
c_vitro_best_border_color <- "#FFFFFF00"
c_vitro_best_border_stroke <- 0
ratio_marker_border_color <- "#FFFFFF"
d_reference_linewidth <- 0.275
d_range_linewidth <- 0.625
efg_median_linewidth <- 0.56 * 1.50
efg_individual_linewidth <- 0.80 * efg_median_linewidth
distribution_initial_linewidth <- 0.22
distribution_endpoint_underlay_linewidth <- 0.38
distribution_endpoint_linewidth <- 0.24
distribution_initial_fill <- "#B8BEC7"
joint_bound_edge_color <- "#D73027"
joint_bound_edge_alpha <- 0.13
d_endpoint_median_shape <- 8
d_endpoint_marker_y <- 0
d_endpoint_median_size <- 1.60
d_internal_font_scale <- 2.70
d_family_panel_spacing_mm <- 1.10
d_family_header_gap <- 0.025
d_parameter_annotation_ymin <- 0.00
d_parameter_background_xmax <- 1.20
d_parameter_name_size_pt <- 11.30
d_parameter_description_size_pt <- 8.60
d_parameter_name_y <- 0.73
d_parameter_description_y <- 0.27
d_optimizer_legend_label <- "* = optimizer median\non axis"
d_bound_legend_label <- "Outer 5% of\njoint bound"
d_initial_distribution_legend_label <- "DE initial distribution"
d_family_outline_legend_label <- "Family outline\ncolor"
d_natural_axis_legend_label <- "Natural scale;\nnatural tick labels"
d_family_header_size_pt <- 12.00
d_context_marker_linewidth <- 1.86
d_parameter_group_indicator_linewidth <- 1.86
figure_font_family <- "Helvetica"
subplot_background_fill <- "transparent"
bef_axis_color <- "#27313A"
bef_axis_linewidth <- 0.65
bef_axis_tick_linewidth <- 0.50
bef_axis_tick_length_mm <- 0.90
panel_title_size <- 15.5
fit_subplot_title_size <- 10.5
b_subplot_title_hjust <- 0.5
b_starting_ploidy_legend_title <- "Starting\nploidy"
c_landscape_body_widths <- c(scatter = 0.88, legend = 0.12)
c_initial_point_size <- 0.12 * 1.50
c_initial_legend_point_size <- 2.50 * 1.50
c_d_overlap_height <- 0.235
axis_title_size <- 10.2
axis_text_size <- 9.2
strip_text_size <- 9.6
panel_title_geom_size <- panel_title_size / ggplot2::.pt
a_height_scale <- 1.15
f_height_scale <- 1.20
assembled_width_scale <- 0.80
# Preserve the released per-family width from the six-family layout. Fixed
# annotation/legend space remains constant while the faceted distribution
# region scales with the manifest-declared family count.
d_distribution_fraction_at_six_families <- 0.91808 * 0.816
assembled_width_inches <- 24.00 * assembled_width_scale * (
  1 - d_distribution_fraction_at_six_families +
    d_distribution_fraction_at_six_families * primary_family_count / 6
)
assembled_height_inches <- 18.00
assembled_row_heights <- c(
  A = 1.15,
  BC = 2.35,
  D = 5.25,
  EF = 1.85
)
survival_legend_position <- c(0.72, 0.14)
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
  header = 0.045,
  body = 0.825,
  note = 0.130
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
      plot.background = element_rect(
        fill = subplot_background_fill,
        color = NA
      ),
      legend.background = element_rect(
        fill = subplot_background_fill,
        color = NA
      ),
      legend.box.background = element_rect(
        fill = subplot_background_fill,
        color = NA
      ),
      legend.key = element_rect(
        fill = subplot_background_fill,
        color = NA
      ),
      plot.margin = margin(3, 4, 3, 4)
    )
}

bef_emphasized_axis_theme <- theme(
  axis.line.x.bottom = element_line(
    color = bef_axis_color,
    linewidth = bef_axis_linewidth,
    lineend = "square"
  ),
  axis.line.y.left = element_line(
    color = bef_axis_color,
    linewidth = bef_axis_linewidth,
    lineend = "square"
  ),
  axis.ticks = element_line(
    color = bef_axis_color,
    linewidth = bef_axis_tick_linewidth,
    lineend = "square"
  ),
  axis.ticks.length = grid::unit(bef_axis_tick_length_mm, "mm")
)

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
      xmax = d_parameter_background_xmax,
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
  context = c("In vitro", "In vivo", "In vitro", "In vivo"),
  title = c(
    "Separate\nin vitro fit",
    "Separate\nin vivo fit",
    paste0("Common ", INVITRO_VISUALIZATION_SEED, "\nculture anchor"),
    paste0(primary_family_count, " landscape\nrepresentatives")
  ),
  detail = c(
    "500-seed search",
    "500-seed search",
    paste0("shared by all ", primary_family_count, " pairs"),
    "one per C family"
  ),
  stringsAsFactors = FALSE
)
workflow_context_boxes$color <- context_colors[workflow_context_boxes$context]

workflow_core_boxes <- data.frame(
  xmin = c(40.5, 57.5, 77.0, 90.0),
  xmax = c(52.5, 72.0, 87.0, 99.5),
  ymin = c(8.1, 8.1, 8.1, 8.1),
  ymax = c(18.3, 18.3, 18.3, 18.3),
  title = c(
    paste0(primary_family_count, " warm-start\npairs"),
    "14 paired\nparameters",
    paste0(format(500L * primary_family_count, big.mark = ","), " joint\nsearches"),
    paste0(primary_family_count, " retained\nwinners")
  ),
  detail = c(
    "landscape-informed",
    "",
    "500 starts / pair",
    "one selected / C family;\n500 endpoints each"
  ),
  stringsAsFactors = FALSE
)

workflow_output_boxes <- data.frame(
  xmin = c(58, 72.2, 84.5),
  xmax = c(70, 82.5, 99),
  ymin = c(0.2, 0.2, 0.2),
  ymax = c(4.2, 4.2, 4.2),
  title = c(
    "B.  Model performance",
    "E. and F.  Fitted functions",
    "D.  Context-specific\nstart vs endpoint"
  ),
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
  xend = c(57.0, 76.5, 89.5),
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
workflow_welsch_label_size <- 2.70
workflow_welsch_detail_size <- 2.45

p_a_body <- ggplot() +
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
    x = 94.75, xend = 94.75, y = 8.0, yend = 6.70,
    color = "#64748B",
    linewidth = 0.42
  ) +
  annotate(
    "segment",
    x = 64, xend = 94.75, y = 6.70, yend = 6.70,
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
    size = 2.90,
    lineheight = 0.92
  ) +
  geom_text(
    data = workflow_context_boxes,
    aes(x = (xmin + xmax) / 2, y = ymin + 2.15, label = detail),
    color = "#475569",
    size = 2.50
  ) +
  geom_text(
    data = workflow_core_boxes,
    aes(x = (xmin + xmax) / 2, y = ymax - 2.45, label = title),
    color = "#0F172A",
    fontface = "bold",
    size = 2.75,
    lineheight = 0.95
  ) +
  geom_text(
    data = workflow_core_boxes,
    aes(x = (xmin + xmax) / 2, y = ymin + 2.1, label = detail),
    color = "#475569",
    size = 2.35,
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
    size = 2.65
  ) +
  annotate(
    "text",
    x = 27.5, y = 1.8,
    label = "Same chromosome-state model; context-specific observation models",
    color = "#64748B", size = 2.80, fontface = "italic"
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 23), clip = "off") +
  theme_void(base_size = 8, base_family = figure_font_family) +
  theme(
    plot.margin = margin(2, 3, 1, 12.7)
  )
p_a_header <- make_text_band(
  "A.  Joint-fitting workflow",
  size = panel_title_size,
  hjust = 0,
  x = 0,
  y = 0.58
)
p_a <- wrap_plots(
  wrap_elements(full = p_a_header, clip = FALSE),
  wrap_elements(full = p_a_body, clip = FALSE),
  ncol = 1,
  heights = c(0.14, 0.86)
)

# -------------------------------------------------------------------------
# Panel B: direct in-sample model performance across retained winners.
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
    "winner", "region", "cohort", "obs_id", "obs_norm", "pred_norm"
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
    "cohort",
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
    "cohort",
    "obs_id",
    "observed_mean_kary_N",
    "predicted_mean_kary_N"
  )]
}

burden_fit <- bind_rows(burden_rows)
terminal_fit <- bind_rows(terminal_rows)
growth_fit <- bind_rows(growth_rows)
kary_fit <- bind_rows(kary_rows)

assert_selected_per_observation <- function(df, label) {
  counts <- table(df$obs_id)
  if (length(counts) == 0L || any(counts != nrow(selected))) {
    stop(
      label,
      " does not contain exactly ",
      nrow(selected),
      " retained winner predictions per observation"
    )
  }
}
assert_selected_per_observation(burden_fit, "In-vivo burden fit")
assert_selected_per_observation(terminal_fit, "In-vivo terminal-ploidy fit")
assert_selected_per_observation(growth_fit, "In-vitro growth fit")
assert_selected_per_observation(kary_fit, "In-vitro karyotype fit")

make_fit_scatter <- function(
  df,
  observed_col,
  predicted_col,
  title,
  x_label,
  y_label,
  color,
  limits,
  show_range = FALSE,
  show_starting_ploidy_legend = FALSE
) {
  summary_df <- df %>%
    group_by(obs_id, cohort) %>%
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
      aes(
        x = .data[[observed_col]],
        y = .data[[predicted_col]],
        shape = cohort
      ),
      fill = alpha(color, 0.14),
      color = alpha(color, 0.24),
      stroke = 0.20,
      size = 0.92,
      show.legend = FALSE
    ) +
    geom_point(
      data = summary_df,
      aes(x = observed, y = predicted, shape = cohort),
      fill = color,
      color = "white",
      stroke = 0.32,
      size = 1.75,
      show.legend = show_starting_ploidy_legend
    ) +
    scale_shape_manual(
      values = c("2N" = 21, "4N" = 24),
      breaks = c("2N", "4N"),
      name = b_starting_ploidy_legend_title,
      drop = FALSE,
      guide = guide_legend(ncol = 1, byrow = TRUE)
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
        size = fit_subplot_title_size,
        face = "bold",
        hjust = b_subplot_title_hjust,
        lineheight = 0.90,
        margin = margin(b = 0.55)
      ),
      axis.title.x = element_text(
        family = figure_font_family,
        size = axis_title_size,
        margin = margin(t = 0.55)
      ),
      axis.title.y = element_text(
        family = figure_font_family,
        size = axis_title_size,
        margin = margin(r = 0.55)
      ),
      axis.text.x = element_text(
        family = figure_font_family,
        size = axis_text_size,
        margin = margin(t = 0.25)
      ),
      axis.text.y = element_text(
        family = figure_font_family,
        size = axis_text_size,
        margin = margin(r = 0.25)
      ),
      plot.margin = margin(0.15, 0.15, 0.15, 0.15)
    ) +
    bef_emphasized_axis_theme
}

p_a_burden <- make_fit_scatter(
  burden_fit,
  "obs_norm",
  "pred_norm",
  "In vivo\nNormalized burden",
  "Observed",
  "Predicted",
  in_vivo_color,
  c(0, 1),
  show_range = FALSE,
  show_starting_ploidy_legend = TRUE
)
p_a_terminal <- make_fit_scatter(
  terminal_fit,
  "observed_mean_N",
  "predicted_mean_N",
  "In vivo\nTerminal mean N",
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
  "In vitro\nGrowth rate",
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
  "In vitro\nMean N",
  "Observed mean N",
  "Predicted mean N",
  in_vitro_color,
  c(44, 92),
  show_range = TRUE
)

p_b_header <- make_text_band(
  "B.  In-sample model performance",
  size = panel_title_size,
  hjust = 0,
  x = 0,
  y = 0.70
)
p_b_panel_order <- c(
  "in_vitro_growth", "in_vitro_mean_N",
  "in_vivo_burden", "in_vivo_terminal_mean_N"
)
p_b_content <- wrap_plots(
  p_a_growth,
  p_a_kary,
  p_a_burden,
  p_a_terminal,
  ncol = 4,
  nrow = 1,
  byrow = TRUE,
  widths = rep(1, 4),
  heights = 1,
  guides = "collect"
) &
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(
      family = figure_font_family,
      size = axis_text_size,
      face = "bold"
    ),
    legend.text = element_text(
      family = figure_font_family,
      size = axis_text_size
    ),
    legend.key.height = grid::unit(0.28, "cm"),
    legend.key.width = grid::unit(0.28, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )
p_b_note <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 0.90,
    label = paste0(
      "Circle: 2N start; triangle: 4N start. ",
      "Pale: six selected fits; filled: median; bars: selected-fit range."
    ),
    hjust = 0,
    size = 3.10,
    color = "#4B5563"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 3, 0, 3))
p_b <- wrap_plots(
  wrap_elements(full = p_b_header, clip = FALSE),
  wrap_elements(full = p_b_content, clip = FALSE),
  wrap_elements(full = p_b_note, clip = FALSE),
  ncol = 1,
  heights = c(0.060, 0.900, 0.040)
)

# -------------------------------------------------------------------------
# Panel C: pooled separate-fit t-SNE landscape used for warm-start selection.
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
tsne_cell_counts <- c(
  invivo_best = sum(
    tsne_full$dataset == "invivo" & tsne_full$point_type == "best"
  ),
  invitro_best = sum(
    tsne_full$dataset == "invitro" & tsne_full$point_type == "best"
  ),
  invivo_initial = sum(
    tsne_full$dataset == "invivo" & tsne_full$point_type == "initial"
  ),
  invitro_initial = sum(
    tsne_full$dataset == "invitro" & tsne_full$point_type == "initial"
  )
)
if (nrow(tsne_full) != 256000L ||
    !identical(unname(tsne_cell_counts), c(500L, 500L, 127500L, 127500L))) {
  stop(
    "Unexpected pooled t-SNE universe for the NP=256 separate fits; ",
    "expected 127,500 initial and 500 best coordinates per context."
  )
}

tsne_initial <- tsne_full[tsne_full$point_type == "initial", , drop = FALSE]
tsne_initial$context <- ifelse(tsne_initial$dataset == "invivo", "In vivo", "In vitro")
tsne_best_vivo <- tsne_best[tsne_best$dataset == "invivo", , drop = FALSE]
tsne_best_vitro <- tsne_full[
  tsne_full$dataset == "invitro" & tsne_full$point_type == "best",
  , drop = FALSE
]
if (nrow(tsne_best_vivo) != 500L || nrow(tsne_best_vitro) != 500L ||
    any(!is.finite(tsne_best_vivo$objective)) ||
    any(!is.finite(tsne_best_vitro$objective))) {
  stop("Figure 5 requires 500 objective-bearing best fits per context.")
}

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
    size = c_initial_point_size,
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
      override.aes = list(
        shape = c(16, 17),
        size = c_initial_legend_point_size,
        alpha = 0.85
      )
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
    color = c_vitro_best_border_color,
    stroke = c_vitro_best_border_stroke,
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
    x = rep(c(0.17, 0.50, 0.83), times = 2L),
    y = rep(c(0.68, 0.32), each = 3L)
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

make_vertical_objective_steps <- function(
    values, low, high, xmin, xmax, ymin, ymax) {
  value_range <- range(values, finite = TRUE)
  step_fraction <- seq(0, 1, length.out = 160L)
  data.frame(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin + step_fraction * (ymax - ymin),
    ymax = ymin + pmin(1, step_fraction + 1 / 159) * (ymax - ymin),
    color = scales::gradient_n_pal(c(low, high))(step_fraction),
    value = value_range[[1L]] + step_fraction * diff(value_range),
    stringsAsFactors = FALSE
  )
}
c_vivo_objective_range <- range(tsne_best_vivo$objective, finite = TRUE)
c_vitro_objective_range <- range(tsne_best_vitro$objective, finite = TRUE)
c_objective_vertical_steps <- rbind(
  make_vertical_objective_steps(
    tsne_best_vivo$objective, "#2C7BB6", "#FDE725",
    0.08, 0.18, 0.76, 0.87
  ),
  make_vertical_objective_steps(
    tsne_best_vitro$objective, "#1A9850", "#D73027",
    0.08, 0.18, 0.57, 0.68
  )
)
c_family_key_ncol <- 2L
c_vertical_cluster_key_data <- transform(
  cluster_labels,
  x = rep(c(0.14, 0.48), length.out = length(cluster_id)),
  y = rep(
    seq(0.245, 0.105, length.out = ceiling(length(cluster_id) / 2L)),
    each = 2L,
    length.out = length(cluster_id)
  )
)
p_c_vertical_key <- ggplot() +
  annotate(
    "text", x = 0.03, y = 0.995, label = "Objective",
    hjust = 0, vjust = 1, family = figure_font_family,
    fontface = "bold", size = 2.75, color = "#27313A"
  ) +
  geom_rect(
    data = c_objective_vertical_steps,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = color),
    inherit.aes = FALSE,
    color = NA
  ) +
  annotate(
    "rect", xmin = 0.08, xmax = 0.18, ymin = 0.76, ymax = 0.87,
    fill = NA, color = "#6B7280", linewidth = 0.24
  ) +
  annotate(
    "rect", xmin = 0.08, xmax = 0.18, ymin = 0.57, ymax = 0.68,
    fill = NA, color = "#6B7280", linewidth = 0.24
  ) +
  annotate(
    "text", x = 0.08, y = 0.915, label = "In vivo",
    hjust = 0, family = figure_font_family, fontface = "bold", size = 2.20,
    color = "#27313A"
  ) +
  annotate(
    "text", x = 0.08, y = 0.720, label = "In vitro",
    hjust = 0, family = figure_font_family, fontface = "bold", size = 2.20,
    color = "#27313A"
  ) +
  annotate(
    "text", x = 0.22, y = c(0.87, 0.76),
    label = sprintf("%.1f", rev(c_vivo_objective_range)),
    hjust = 0, family = figure_font_family, size = 2.05,
    color = "#4B5563"
  ) +
  annotate(
    "text", x = 0.22, y = c(0.68, 0.57),
    label = sprintf("%.1f", rev(c_vitro_objective_range)),
    hjust = 0, family = figure_font_family, size = 2.05,
    color = "#4B5563"
  ) +
  annotate(
    "text", x = 0.03, y = 0.505, label = "Initial samples",
    hjust = 0, family = figure_font_family, fontface = "bold",
    size = 2.55, color = "#27313A"
  ) +
  annotate(
    "point", x = 0.14, y = 0.435, shape = 16,
    size = 2.4, color = in_vivo_color
  ) +
  annotate(
    "text", x = 0.27, y = 0.435, label = "In vivo",
    hjust = 0, family = figure_font_family, size = 2.30,
    color = "#27313A"
  ) +
  annotate(
    "point", x = 0.14, y = 0.375, shape = 17,
    size = 2.6, color = in_vitro_color
  ) +
  annotate(
    "text", x = 0.27, y = 0.375, label = "In vitro",
    hjust = 0, family = figure_font_family, size = 2.30,
    color = "#27313A"
  ) +
  annotate(
    "text", x = 0.03, y = 0.295, label = "Warm-start family",
    hjust = 0, family = figure_font_family, fontface = "bold",
    size = 2.55, color = "#27313A"
  ) +
  geom_label(
    data = c_vertical_cluster_key_data,
    aes(x = x, y = y, label = label, color = cluster_id),
    inherit.aes = FALSE,
    fill = "white",
    text.colour = cluster_label_text_color,
    linewidth = cluster_label_border_linewidth,
    size = 2.35,
    fontface = cluster_label_fontface,
    label.padding = grid::unit(0.08, "lines"),
    label.r = grid::unit(0.05, "lines")
  ) +
  scale_fill_identity() +
  scale_color_manual(values = cluster_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.08, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(1, 1, 1, 1))

c_legend_position <- "right"
p_landscape <- p_c_scatter +
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5)
  )
p_landscape_body <- wrap_plots(
  wrap_elements(full = p_landscape, clip = FALSE),
  wrap_elements(full = p_c_vertical_key, clip = FALSE),
  ncol = 2,
  widths = unname(c_landscape_body_widths)
)
p_landscape_header <- make_text_band(
  "C.  Warm-start map",
  size = panel_title_size,
  hjust = 0,
  x = 0,
  y = 0.70
)

# -------------------------------------------------------------------------
# Archived Panel F implementation with an induced prior/reference. Retained
# only as code provenance; the active Figure 5D implementation follows it.
# The contrast is log2(in vivo / in vitro). The prior/reference distribution
# is the pushforward of the actual bounded transformed-scale base measure,
# Welsch coupling penalty, and (for five parameters) the configured in-vivo
# Gaussian soft prior. Optimizer endpoints remain numerical-search outputs,
# not Bayesian posterior samples or biological replicates.
# -------------------------------------------------------------------------
if (FALSE) {

coupled_parameters <- c(
  "lam_max", "p_mis_base", "p_wgd",
  "p_misseg", "k_o_mis",
  "O2_crit", "n_O",
  "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
parameter_lookup <- parameter_meta_all[
  match(coupled_parameters, parameter_meta_all$parameter),
  c(
    "parameter", "parameter_group", "parameter_order", "transformation"
  ),
  drop = FALSE
]
if (anyNA(parameter_lookup$parameter) ||
    !identical(parameter_lookup$parameter, coupled_parameters)) {
  stop("Figure 5D parameter order is not aligned with the shared taxonomy")
}
family_levels <- c("C01", "C02", "C03")
mechanism_short <- c(
  "Division and genome change" = "Division/genome",
  "Stress-linked missegregation" = "Stress-MS",
  "Oxygen response" = "Oxygen",
  "Ploidy-linked growth and death" = "Ploidy growth/death",
  "Post-missegregation survival" = "Post-MS survival"
)
parameter_lookup <- merge(
  parameter_lookup,
  prior_sampling_config[, c("parameter", "soft_prior_included")],
  by = "parameter",
  all.x = TRUE,
  sort = FALSE
)
parameter_lookup <- parameter_lookup[
  match(coupled_parameters, parameter_lookup$parameter),
  ,
  drop = FALSE
]
parameter_lookup$facet_label <- paste0(
  parameter_lookup$parameter,
  ifelse(parameter_lookup$soft_prior_included, "*", ""),
  "  |  ",
  mechanism_short[parameter_lookup$parameter_group]
)
parameter_lookup$facet_label <- factor(
  parameter_lookup$facet_label,
  levels = parameter_lookup$facet_label
)

add_facet_label <- function(data) {
  data$facet_label <- parameter_lookup$facet_label[
    match(data$parameter, parameter_lookup$parameter)
  ]
  data
}
prior_density <- add_facet_label(prior_density)
prior_solution_summary <- add_facet_label(prior_solution_summary)
pair_solution_summary <- add_facet_label(pair_solution_summary)
cross_family_summary <- add_facet_label(cross_family_summary)

prior_density_plot <- prior_density[
  prior_density$distribution_role == "prior",
  ,
  drop = FALSE
]
pooled_density_plot <- prior_density[
  prior_density$distribution_role == "solution_pooled",
  ,
  drop = FALSE
]
family_density_plot <- prior_density[
  prior_density$distribution_role == "solution_family",
  ,
  drop = FALSE
]
prior_interval <- prior_solution_summary[
  prior_solution_summary$distribution_role == "prior",
  ,
  drop = FALSE
]
pooled_interval <- prior_solution_summary[
  prior_solution_summary$distribution_role == "solution_pooled",
  ,
  drop = FALSE
]
panel_f_summary <- merge(
  cross_family_summary,
  pooled_interval[, c(
    "parameter", "width90_relative_to_prior"
  )],
  by = "parameter",
  all.x = TRUE,
  sort = FALSE
)
panel_f_summary <- merge(
  panel_f_summary,
  parameter_lookup[, c(
    "parameter", "soft_prior_included"
  )],
  by = "parameter",
  all.x = TRUE,
  sort = FALSE
)
direction_sensitivity <- unique(
  sensitivity_validation[, c(
    "parameter", "direction_stable_across_sensitivities"
  )]
)
panel_f_summary <- merge(
  panel_f_summary,
  direction_sensitivity,
  by = "parameter",
  all.x = TRUE,
  sort = FALSE
)
panel_f_summary <- panel_f_summary[
  match(coupled_parameters, panel_f_summary$parameter),
  ,
  drop = FALSE
]
panel_f_summary$facet_label <- factor(
  panel_f_summary$facet_label,
  levels = levels(parameter_lookup$facet_label)
)
panel_f_summary$direction_display <- ifelse(
  panel_f_summary$direction == "higher in vivo",
  "Higher",
  ifelse(
    panel_f_summary$direction == "lower in vivo",
    "Lower",
    "Mixed"
  )
)
panel_f_summary$frequent_boundary <-
  panel_f_summary$max_active_bound_fraction >= 0.10
panel_f_summary$annotation <- paste0(
  panel_f_summary$direction_display,
  " | W90 ",
  sprintf("%.2fx", panel_f_summary$width90_relative_to_prior),
  ifelse(panel_f_summary$frequent_boundary, " | bound", "")
)
direction_colors <- c(
  "Higher" = in_vivo_color,
  "Lower" = in_vitro_color,
  "Mixed" = "#66717E"
)
panel_f_summary$annotation_color <- direction_colors[
  panel_f_summary$direction_display
]

facet_background <- parameter_lookup[, c(
  "parameter", "parameter_group", "facet_label"
)]
distribution_fill_colors <- c(
  "Induced prior/reference" = "#B8BEC7",
  "Optimizer solutions" = "#416B86"
)

p_f_main <- ggplot() +
  geom_rect(
    data = facet_background,
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf,
      fill = parameter_group
    ),
    inherit.aes = FALSE,
    alpha = 0.10,
    color = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = parameter_group_colors, guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_vline(
    xintercept = 0,
    color = "#111827",
    linetype = "22",
    linewidth = 0.34
  ) +
  geom_area(
    data = prior_density_plot,
    aes(
      x = x,
      y = density_scaled,
      group = facet_label,
      fill = "Induced prior/reference"
    ),
    color = "#7D8590",
    linewidth = 0.26,
    alpha = 0.56
  ) +
  geom_area(
    data = pooled_density_plot,
    aes(
      x = x,
      y = density_scaled,
      group = facet_label,
      fill = "Optimizer solutions"
    ),
    color = "#27465A",
    linewidth = 0.30,
    alpha = 0.42
  ) +
  scale_fill_manual(
    values = distribution_fill_colors,
    breaks = names(distribution_fill_colors),
    name = "Distribution",
    guide = guide_legend(
      order = 1,
      nrow = 1,
      override.aes = list(alpha = c(0.72, 0.55))
    )
  ) +
  geom_line(
    data = family_density_plot,
    aes(
      x = x,
      y = density_scaled,
      group = interaction(facet_label, family)
    ),
    color = "#1F2933",
    linewidth = 0.52,
    alpha = 0.92
  ) +
  geom_line(
    data = family_density_plot,
    aes(
      x = x,
      y = density_scaled,
      color = family,
      group = interaction(facet_label, family)
    ),
    linewidth = 0.32,
    alpha = 1
  ) +
  geom_segment(
    data = prior_interval,
    aes(x = q05, xend = q95, y = 0.035, yend = 0.035),
    color = "#6B7280",
    linewidth = 0.55,
    lineend = "round"
  ) +
  geom_segment(
    data = pooled_interval,
    aes(x = q05, xend = q95, y = 0.080, yend = 0.080),
    color = "#27465A",
    linewidth = 0.78,
    lineend = "round"
  ) +
  geom_segment(
    data = pooled_interval,
    aes(x = median, xend = median, y = 0.045, yend = 0.155),
    color = "#111827",
    linewidth = 0.34
  ) +
  geom_point(
    data = pair_solution_summary,
    aes(
      x = median,
      y = 0.12,
      color = family
    ),
    shape = 16,
    size = 0.82,
    stroke = 0.30,
    alpha = 0.98
  ) +
  geom_text(
    data = panel_f_summary,
    aes(
      x = Inf,
      y = 1.08,
      label = annotation,
      color = direction_display
    ),
    inherit.aes = FALSE,
    hjust = 1.04,
    vjust = 1,
    family = figure_font_family,
    fontface = "bold",
    size = 1.55,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(region_colors, direction_colors),
    breaks = family_levels,
    name = "Warm-start family",
    guide = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(linewidth = 0.8, shape = NA)
    )
  ) +
  facet_wrap(
    vars(facet_label),
    ncol = 2,
    scales = "free_x",
    axes = "all_x",
    axis.labels = "all_x"
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 3),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  scale_y_continuous(
    limits = c(0, 1.10),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  labs(
    x = "log2(in vivo / in vitro ratio); 0 = equal",
    y = NULL
  ) +
  theme_manuscript(base_size = 6.2) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      color = alpha("#9AA3AD", 0.22),
      linewidth = 0.22
    ),
    panel.border = element_rect(
      fill = NA,
      color = "#AEB6C0",
      linewidth = 0.25
    ),
    strip.background = element_rect(
      fill = "#F7F9FB",
      color = "#AEB6C0",
      linewidth = 0.25
    ),
    strip.text = element_text(
      family = figure_font_family,
      face = "bold",
      size = 5.15,
      color = "#27313A",
      margin = margin(0.6, 1, 0.6, 1)
    ),
    axis.title.x = element_text(
      family = figure_font_family,
      size = 5.3,
      margin = margin(t = 1.3)
    ),
    axis.text.x = element_text(
      family = figure_font_family,
      size = 4.65,
      color = "#4B5563",
      margin = margin(t = 0.8)
    ),
    axis.ticks.x = element_line(
      color = "#6B7280",
      linewidth = 0.22
    ),
    axis.ticks.length = grid::unit(0.55, "mm"),
    panel.spacing.x = grid::unit(1.7, "mm"),
    panel.spacing.y = grid::unit(0.7, "mm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.title = element_text(
      family = figure_font_family,
      size = 5.1,
      face = "bold"
    ),
    legend.text = element_text(
      family = figure_font_family,
      size = 5.0
    ),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.width = grid::unit(4.8, "mm"),
    legend.spacing.x = grid::unit(1.4, "mm"),
    legend.margin = margin(-1.0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0.5, 2, 0, 2)
  ) +
  guides(
    shape = guide_legend(order = 3, nrow = 1),
    color = guide_legend(order = 2, nrow = 1)
  )

p_f_header <- make_text_band(
  "F  Prior-reference and fitted-solution distributions",
  size = panel_title_size,
  hjust = 0,
  x = 0,
  y = 0.60
)
p_f_note <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 0.78,
    label = paste0(
      "Gray: induced prior/reference; blue: pooled ",
      format(500L * primary_family_count, big.mark = ","),
      " optimizer endpoints; colored outlines: ",
      paste(primary_family_order, collapse = ", "),
      "; circle/triangle: pair medians."
    ),
    hjust = 0,
    size = 1.78,
    family = figure_font_family,
    color = "#4B5563"
  ) +
  annotate(
    "text",
    x = 0,
    y = 0.48,
    label = paste0(
      "* includes the configured in-vivo Gaussian soft prior. ",
      "W90 = optimizer/prior 5-95% width; 'bound' = >=10% active-bound occupancy."
    ),
    hjust = 0,
    size = 1.78,
    family = figure_font_family,
    color = "#4B5563"
  ) +
  annotate(
    "text",
    x = 0,
    y = 0.18,
    label = paste0(
      "Free x scales. Direction labels were unchanged in the lowest 10% per pair ",
      "and 947 unique endpoints; optimizer distributions are not posterior or confidence intervals."
    ),
    hjust = 0,
    size = 1.78,
    family = figure_font_family,
    color = "#4B5563"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 3, 0, 3))
p_ratios <- wrap_plots(
  wrap_elements(full = p_f_header, clip = FALSE),
  wrap_elements(full = p_f_main, clip = FALSE),
  wrap_elements(full = p_f_note, clip = FALSE),
  ncol = 1,
  heights = unname(f_ratio_layout_heights)
)
}

# -------------------------------------------------------------------------
# Panel D: direct in-vivo and in-vitro DE initial populations and optimizer
# endpoints. Each parameter block contains one subpanel per declared family. Context is
# encoded by mirrored blue/pink fill; the actual pair-specific DE initial
# population has a gray dashed outline and the 500 feasible endpoints have the
# corresponding C-family solid outline. Translucent red background bands mark
# the outer 5% at each end of the complete joint-union optimizer interval; they
# are bound-proximity guides, not distribution quantiles or uncertainty bands.
# -------------------------------------------------------------------------

coupled_parameter_universe <- c(
  "lam_max", "p_mis_base", "p_wgd",
  "p_misseg", "k_o_mis",
  "O2_crit", "n_O",
  "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
family_levels <- primary_family_order

# Order rows by the retained joint winners. For each family, use
# the in-vivo/in-vitro ratio from that family's single selected winner, take
# the arithmetic mean across C01--C06, and place parameters closest to parity
# (mean ratio = 1) first. The saved table makes this display-only ordering
# independently auditable.
best_parameter_ratio_long <- bind_rows(lapply(seq_len(nrow(selected)), function(i) {
  ratio_table <- read_winner_table(
    selected$fit_dir[[i]], "joint_soft_coupling.tsv"
  )
  required_ratio_fields <- c(
    "parameter", "vivo_natural", "vitro_natural", "ratio_vivo_to_vitro"
  )
  if (!all(required_ratio_fields %in% names(ratio_table))) {
    stop(
      "joint_soft_coupling.tsv lacks best-winner ratio fields for ",
      selected$warmup_label[[i]]
    )
  }
  ratio_table <- ratio_table[
    match(coupled_parameter_universe, ratio_table$parameter),
    required_ratio_fields,
    drop = FALSE
  ]
  if (anyNA(ratio_table$parameter) ||
      !identical(ratio_table$parameter, coupled_parameter_universe) ||
      any(!is.finite(ratio_table$vivo_natural)) ||
      any(!is.finite(ratio_table$vitro_natural)) ||
      any(ratio_table$vivo_natural <= 0) ||
      any(ratio_table$vitro_natural <= 0)) {
    stop(
      "Best-winner natural-scale values are incomplete or invalid for ",
      selected$warmup_label[[i]]
    )
  }
  ratio_table$ratio_recomputed <- with(
    ratio_table, vivo_natural / vitro_natural
  )
  if (any(!is.finite(ratio_table$ratio_vivo_to_vitro)) ||
      any(abs(log(
        ratio_table$ratio_recomputed / ratio_table$ratio_vivo_to_vitro
      )) > 1e-10)) {
    stop(
      "Stored and recomputed best-winner ratios disagree for ",
      selected$warmup_label[[i]]
    )
  }
  ratio_table$family <- selected$region[[i]]
  ratio_table$warmup_label <- selected$warmup_label[[i]]
  ratio_table$joint_seed <- selected$selected_seed[[i]]
  ratio_table[, c(
    "parameter", "family", "warmup_label", "joint_seed",
    "vivo_natural", "vitro_natural", "ratio_recomputed"
  )]
}))
if (nrow(best_parameter_ratio_long) != 14L * primary_family_count ||
    !setequal(best_parameter_ratio_long$family, family_levels) ||
    any(table(best_parameter_ratio_long$parameter) != primary_family_count)) {
  stop("Figure 5D best-winner ratio table is incomplete")
}

best_parameter_ratio_order <- aggregate(
  ratio_recomputed ~ parameter,
  best_parameter_ratio_long,
  mean
)
names(best_parameter_ratio_order)[[2L]] <- "mean_ratio_vivo_to_vitro"
family_ratio_columns <- paste0(family_levels, "_ratio_vivo_to_vitro")
for (family_index in seq_along(family_levels)) {
  family <- family_levels[[family_index]]
  family_rows <- best_parameter_ratio_long[
    best_parameter_ratio_long$family == family,
    c("parameter", "ratio_recomputed"),
    drop = FALSE
  ]
  best_parameter_ratio_order[[family_ratio_columns[[family_index]]]] <-
    family_rows$ratio_recomputed[
      match(best_parameter_ratio_order$parameter, family_rows$parameter)
    ]
}
best_parameter_ratio_order$distance_of_mean_ratio_from_1 <- abs(
  best_parameter_ratio_order$mean_ratio_vivo_to_vitro - 1
)
best_parameter_ratio_order <- best_parameter_ratio_order[
  order(
    best_parameter_ratio_order$distance_of_mean_ratio_from_1,
    best_parameter_ratio_order$parameter
  ),
  ,
  drop = FALSE
]
best_parameter_ratio_order$row_order <- seq_len(nrow(best_parameter_ratio_order))
best_parameter_ratio_order <- best_parameter_ratio_order[, c(
  "row_order", "parameter", family_ratio_columns,
  "mean_ratio_vivo_to_vitro", "distance_of_mean_ratio_from_1"
), drop = FALSE]
coupled_parameters <- best_parameter_ratio_order$parameter
parameter_row_order_path <- file.path(
  figure5_data_dir, "figure5d_parameter_row_order.tsv"
)
write.table(
  best_parameter_ratio_order,
  parameter_row_order_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
mechanism_short <- c(
  "Division and genome change" = "Division/genome",
  "Stress-linked missegregation" = "Stress-MS",
  "Oxygen response" = "Oxygen",
  "Ploidy-linked growth and death" = "Ploidy growth/death",
  "Post-missegregation survival" = "Post-MS survival"
)
parameter_description_short <- c(
  "lam_max" = "Maximum division rate",
  "p_mis_base" = "Baseline per-chromosome missegregation",
  "p_wgd" = "Whole-genome-doubling probability",
  "p_misseg" = "Maximum stress-induced missegregation",
  "k_o_mis" = "Stress-linked MS half-saturation",
  "O2_crit" = "Critical oxygen level",
  "n_O" = "Oxygen-response steepness",
  "alpha_o2" = "Stress-linked high-ploidy growth damping",
  "gamma_growth" = "Ploidy scaling of growth penalty",
  "mu_hp" = "Stress-associated death scale",
  "gamma_mu" = "Ploidy scaling of stress death",
  "buffer_smax" = "Maximum post-MS survival",
  "buffer_beta" = "Ploidy-dependent viability loss",
  "buffer_n_exp" = "Ploidy scaling of post-MS survival"
)
parameter_description_short <- parameter_description_short[coupled_parameters]
if (!identical(names(parameter_description_short), coupled_parameters)) {
  stop("Figure 5D compact parameter descriptions are incomplete or misordered")
}
parameter_lookup <- parameter_meta_all[
  match(coupled_parameters, parameter_meta_all$parameter),
  c(
    "parameter", "parameter_group", "parameter_description",
    "transformation", "parameter_order"
  ),
  drop = FALSE
]
if (anyNA(parameter_lookup$parameter) ||
    !identical(parameter_lookup$parameter, coupled_parameters)) {
  stop("Figure 5D parameter order is not aligned with the shared taxonomy")
}

context_levels <- c("in vitro", "in vivo")
context_fill_colors <- c(
  "in vivo" = "#0072B2",
  "in vitro" = "#CC79A7"
)
context_density$family <- factor(context_density$family, levels = family_levels)
context_density$context <- factor(
  context_density$context, levels = context_levels
)
context_summary$family <- factor(context_summary$family, levels = family_levels)
context_summary$context <- factor(
  context_summary$context, levels = context_levels
)

format_context_ticks <- function(x, transformation) {
  values <- if (identical(transformation, "log10")) 10^x else x
  vapply(values, function(value) {
    if (!is.finite(value)) return(NA_character_)
    if ((abs(value) > 0 && abs(value) < 0.01) || abs(value) >= 1000) {
      return(formatC(value, format = "e", digits = 1))
    }
    format(signif(value, 2), trim = TRUE, scientific = FALSE)
  }, character(1L))
}

to_panel_scale <- function(values, transformation) {
  if (identical(transformation, "log10")) log10(values) else values
}

context_axis_breaks <- function(limits) {
  span <- diff(limits)
  if (!is.finite(span) || span <= 0) return(mean(limits))
  candidates <- pretty(limits, n = 2)
  candidates <- candidates[
    candidates > limits[[1L]] + 0.10 * span &
      candidates < limits[[2L]] - 0.10 * span
  ]
  if (!length(candidates)) return(mean(limits))
  if (length(candidates) <= 2L) return(candidates)
  candidates[unique(round(seq(1, length(candidates), length.out = 2L)))]
}

make_family_distribution_panel <- function(parameter) {
  density_data <- context_density[
    context_density$parameter == parameter,
    ,
    drop = FALSE
  ]
  initial_data <- density_data[
    density_data$distribution_role == "de_initial_population",
    , drop = FALSE
  ]
  optimizer_data <- density_data[
    density_data$distribution_role == "optimizer_endpoints",
    , drop = FALSE
  ]
  interval_data <- context_summary[
    context_summary$parameter == parameter,
    ,
    drop = FALSE
  ]
  meta <- parameter_lookup[
    parameter_lookup$parameter == parameter,
    ,
    drop = FALSE
  ]
  if (nrow(density_data) != primary_family_count * 2L * 2L * 401L ||
      nrow(initial_data) != primary_family_count * 2L * 401L ||
      nrow(optimizer_data) != primary_family_count * 2L * 401L ||
      nrow(interval_data) != primary_family_count * 2L ||
      length(unique(density_data$display_from)) != 1L ||
      length(unique(density_data$display_to)) != 1L ||
      max(density_data$displayed_tail_fraction) > 0.002) {
    stop("Figure 5D nested-panel contract failed for ", parameter)
  }
  initial_data$plot_density <- ifelse(
    initial_data$context == "in vitro",
    initial_data$density_scaled_distribution,
    -initial_data$density_scaled_distribution
  )
  optimizer_data$plot_density <- ifelse(
    optimizer_data$context == "in vitro",
    optimizer_data$density_scaled_distribution,
    -optimizer_data$density_scaled_distribution
  )
  transformation <- density_data$transformation[[1L]]
  interval_data$optimizer_median_display <- to_panel_scale(
    interval_data$optimizer_median, transformation
  )
  bound_data <- joint_parameter_bounds[
    joint_parameter_bounds$parameter == parameter,
    ,
    drop = FALSE
  ]
  if (nrow(bound_data) != 1L) {
    stop("Figure 5D requires one complete joint bound for ", parameter)
  }
  display_limits <- c(
    bound_data$joint_union_lower_transformed,
    bound_data$joint_union_upper_transformed
  )
  ggplot() +
    annotate(
      "rect",
      xmin = c(
        bound_data$joint_union_lower_transformed,
        bound_data$upper_edge_start_transformed
      ),
      xmax = c(
        bound_data$lower_edge_end_transformed,
        bound_data$joint_union_upper_transformed
      ),
      ymin = -Inf,
      ymax = Inf,
      fill = joint_bound_edge_color,
      alpha = joint_bound_edge_alpha,
      color = NA
    ) +
    geom_hline(
      yintercept = 0,
      color = "#87909A",
      linewidth = 0.20
    ) +
    geom_ribbon(
      data = initial_data,
      aes(
        x = x,
        ymin = pmin(0, plot_density),
        ymax = pmax(0, plot_density),
        group = interaction(family, context)
      ),
      fill = distribution_initial_fill,
      color = NA,
      alpha = 0.30
    ) +
    geom_line(
      data = initial_data,
      aes(
        x = x,
        y = plot_density,
        group = interaction(family, context)
      ),
      color = "#6F7780",
      linewidth = distribution_initial_linewidth,
      linetype = "22"
    ) +
    geom_ribbon(
      data = optimizer_data,
      aes(
        x = x,
        ymin = pmin(0, plot_density),
        ymax = pmax(0, plot_density),
        fill = context,
        group = interaction(family, context)
      ),
      color = NA,
      alpha = 0.28
    ) +
    geom_line(
      data = optimizer_data,
      aes(
        x = x,
        y = plot_density,
        color = family,
        group = interaction(family, context)
      ),
      linewidth = distribution_endpoint_linewidth,
      linetype = "solid",
      alpha = 0.98
    ) +
    geom_point(
      data = interval_data,
      aes(
        x = optimizer_median_display,
        y = d_endpoint_marker_y
      ),
      color = unname(context_fill_colors[as.character(interval_data$context)]),
      shape = d_endpoint_median_shape,
      size = d_endpoint_median_size,
      stroke = 0.42
    ) +
    scale_fill_manual(
      values = context_fill_colors,
      limits = context_levels,
      guide = "none"
    ) +
    scale_color_manual(
      values = region_colors,
      limits = family_levels,
      guide = "none"
    ) +
    facet_grid(cols = vars(family), scales = "fixed") +
    scale_x_continuous(
      breaks = context_axis_breaks(display_limits),
      labels = function(x) format_context_ticks(x, transformation),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(-1.03, 1.03),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = display_limits,
      ylim = c(-1.03, 1.03),
      expand = FALSE,
      clip = "on"
    ) +
    labs(x = NULL, y = NULL) +
    theme_manuscript(base_size = 5.2) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(
        color = alpha("#9AA3AD", 0.20),
        linewidth = 0.18
      ),
      panel.border = element_rect(
        fill = NA,
        color = "#AEB6C0",
        linewidth = 0.22
      ),
      panel.background = element_rect(fill = "white", color = NA),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(
        family = figure_font_family,
        size = 7.00,
        color = "#4B5563",
        margin = margin(t = 0.35)
      ),
      axis.ticks.x = element_line(
        color = "#6B7280",
        linewidth = 0.18
      ),
      axis.ticks.length = grid::unit(0.35, "mm"),
      panel.spacing.x = grid::unit(d_family_panel_spacing_mm, "mm"),
      legend.position = "bottom",
      legend.title = element_text(
        family = figure_font_family,
        face = "bold",
        size = 4.2
      ),
      legend.text = element_text(
        family = figure_font_family,
        size = 4.1
      ),
      legend.key.height = grid::unit(2.0, "mm"),
      legend.key.width = grid::unit(2.4, "mm"),
      legend.margin = margin(-1, 0, 0, 0),
      plot.margin = margin(0.15, 0.35, 0.15, 0.35)
    )
}

make_context_direction_panel <- function() {
  ggplot() +
    annotate(
      "segment",
      x = 0.50,
      xend = 0.50,
      y = 0.61,
      yend = 0.94,
      color = context_fill_colors[["in vitro"]],
      linewidth = d_context_marker_linewidth,
      lineend = "round"
    ) +
    annotate(
      "segment",
      x = 0.50,
      xend = 0.50,
      y = d_parameter_annotation_ymin + 0.08,
      yend = 0.51,
      color = context_fill_colors[["in vivo"]],
      linewidth = d_context_marker_linewidth,
      lineend = "round"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void(base_family = figure_font_family) +
    theme(plot.margin = margin(0, 0.20, 0, 0.10))
}

make_parameter_description_panel <- function(parameter) {
  meta <- parameter_lookup[
    parameter_lookup$parameter == parameter,
    ,
    drop = FALSE
  ]
  description_text <- paste0(
    parameter_description_short[[parameter]],
    if (identical(meta$transformation[[1L]], "log10")) {
      "  [log x]"
    } else {
      "  [linear x]"
    }
  )
  group_color <- parameter_group_colors[[meta$parameter_group]]
  group_fill <- scales::alpha(group_color, 0.065)

  ggplot() +
    annotate(
      "rect",
      xmin = 0,
      xmax = d_parameter_background_xmax,
      ymin = d_parameter_annotation_ymin,
      ymax = 1,
      fill = group_fill,
      color = NA
    ) +
    annotate(
      "segment",
      x = 0.025,
      xend = 0.025,
      y = d_parameter_annotation_ymin + 0.08,
      yend = 0.92,
      color = group_color,
      linewidth = d_parameter_group_indicator_linewidth,
      lineend = "round"
    ) +
    annotate(
      "text",
      x = 0.075,
      y = d_parameter_name_y,
      label = parameter,
      hjust = 0,
      size = d_parameter_name_size_pt / ggplot2::.pt,
      family = figure_font_family,
      fontface = "bold",
      color = "#27313A"
    ) +
    annotate(
      "text",
      x = 0.075,
      y = d_parameter_description_y,
      label = description_text,
      hjust = 0,
      size = d_parameter_description_size_pt / ggplot2::.pt,
      family = figure_font_family,
      color = "#4B5563"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void(base_family = figure_font_family) +
    theme(plot.margin = margin(0, 0, 0, 0.35))
}

panel_d_row_widths <- c(
  context = 0.004,
  distributions = 0.816,
  label = 0.150,
  spacer = 0.030
)
family_header_data <- data.frame(
  family = family_levels,
  xmin = seq_along(family_levels) - 1L + d_family_header_gap,
  xmax = seq_along(family_levels) - d_family_header_gap,
  stringsAsFactors = FALSE
)
p_d_family_header_center <- ggplot(family_header_data) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, color = family),
    fill = "#F7F9FB",
    linewidth = 0.42
  ) +
  geom_text(
    aes(x = (xmin + xmax) / 2, y = 0.5, label = family),
    family = figure_font_family,
    fontface = "bold",
    size = d_family_header_size_pt / ggplot2::.pt,
    color = "#27313A"
  ) +
  scale_color_manual(values = region_colors, guide = "none") +
  coord_cartesian(
    xlim = c(0, primary_family_count),
    ylim = c(0, 1),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(0, 0.35, 0, 0.35))
p_d_family_header <- wrap_plots(
  plot_spacer(),
  wrap_elements(full = p_d_family_header_center, clip = FALSE),
  plot_spacer(),
  plot_spacer(),
  ncol = 4,
  widths = unname(panel_d_row_widths)
)

parameter_panels <- lapply(coupled_parameters, function(parameter) {
  wrap_plots(
    wrap_elements(
      full = make_context_direction_panel(),
      clip = FALSE
    ),
    wrap_elements(
      full = make_family_distribution_panel(parameter),
      clip = FALSE
    ),
    wrap_elements(
      full = make_parameter_description_panel(parameter),
      clip = FALSE
    ),
    plot_spacer(),
    ncol = 4,
    widths = unname(panel_d_row_widths)
  )
})
parameter_distribution_row_count <- length(parameter_panels)
distribution_standalone_height_inches <-
  0.40 * parameter_distribution_row_count + 0.90
bc_row_widths <- c(B = 0.70, C = 0.30)
ef_panel_widths <- c(E = 1, F = 1)
# Keep panel D at the full assembled-figure width while giving the compact
# legend another 20% reduction; the released width is transferred directly to
# the distributions so the total remains exactly one Figure-5 width.
d_content_widths <- c(distributions = 0.91808, legend = 0.08192)
d_total_width_scale <- 1.00
d_family_header_height <- 0.12
d_assembled_header_fraction <- 0.045
d_optimizer_black_underlay <- FALSE
p_d_context_main <- wrap_plots(
  c(list(p_d_family_header), parameter_panels),
  ncol = 1,
  heights = c(d_family_header_height, rep(1, parameter_distribution_row_count))
)

p_d_context_header <- make_text_band(
  "D.  Context-specific DE initialization and optimizer endpoints",
  size = panel_title_size,
  hjust = 0,
  x = 0,
  y = 0.60
)
d_endpoint_legend_data <- data.frame(
  family = family_levels,
  y = seq(0.295, 0.160, length.out = length(family_levels)),
  stringsAsFactors = FALSE
)
d_initial_legend_shape <- rbind(
  data.frame(
    context = "in vitro",
    x = c(0.05, 0.09, 0.14, 0.19, 0.24, 0.29, 0.34, 0.34, 0.05),
    y = c(0.665, 0.678, 0.712, 0.740, 0.712, 0.680, 0.665, 0.665, 0.665)
  ),
  data.frame(
    context = "in vivo",
    x = c(0.05, 0.09, 0.14, 0.19, 0.24, 0.29, 0.34, 0.34, 0.05),
    y = c(0.665, 0.652, 0.618, 0.590, 0.618, 0.650, 0.665, 0.665, 0.665)
  )
)
d_endpoint_legend_shape <- rbind(
  data.frame(
    context = "in vitro",
    x = c(0.05, 0.09, 0.14, 0.20, 0.25, 0.30, 0.34, 0.34, 0.05),
    y = c(0.500, 0.512, 0.540, 0.565, 0.545, 0.515, 0.500, 0.500, 0.500)
  ),
  data.frame(
    context = "in vivo",
    x = c(0.05, 0.09, 0.14, 0.20, 0.25, 0.30, 0.34, 0.34, 0.05),
    y = c(0.500, 0.488, 0.460, 0.435, 0.455, 0.485, 0.500, 0.500, 0.500)
  )
)
p_d_compact_legend <- ggplot() +
  annotate(
    "text",
    x = 0.04,
    y = 0.97,
    label = "Context",
    hjust = 0,
    size = 1.72 * d_internal_font_scale,
    family = figure_font_family,
    fontface = "bold",
    color = "#27313A"
  ) +
  annotate(
    "rect",
    xmin = 0.05,
    xmax = 0.10,
    ymin = 0.885,
    ymax = 0.925,
    fill = context_fill_colors[["in vitro"]],
    color = NA
  ) +
  annotate(
    "text",
    x = 0.14,
    y = 0.905,
    label = "in vitro (upper)",
    hjust = 0,
    size = 1.55 * d_internal_font_scale,
    family = figure_font_family,
    color = "#27313A"
  ) +
  annotate(
    "rect",
    xmin = 0.05,
    xmax = 0.10,
    ymin = 0.815,
    ymax = 0.855,
    fill = context_fill_colors[["in vivo"]],
    color = NA
  ) +
  annotate(
    "text",
    x = 0.14,
    y = 0.835,
    label = "in vivo (lower)",
    hjust = 0,
    size = 1.55 * d_internal_font_scale,
    family = figure_font_family,
    color = "#27313A"
  ) +
  annotate(
    "text",
    x = 0.04,
    y = 0.775,
    label = "Distribution",
    hjust = 0,
    size = 1.72 * d_internal_font_scale,
    family = figure_font_family,
    fontface = "bold",
    color = "#27313A"
  ) +
  geom_polygon(
    data = d_initial_legend_shape,
    aes(x = x, y = y, group = context),
    inherit.aes = FALSE,
    fill = distribution_initial_fill,
    color = "#6F7780",
    alpha = 0.48,
    linewidth = distribution_initial_linewidth,
    linetype = "22"
  ) +
  annotate(
    "text",
    x = 0.39,
    y = c(0.710, 0.620),
    label = c("in vitro", "in vivo"),
    hjust = 0,
    size = 1.25 * d_internal_font_scale,
    family = figure_font_family,
    color = "#27313A"
  ) +
  annotate(
    "text",
    x = 0.05,
    y = 0.575,
    label = d_initial_distribution_legend_label,
    hjust = 0,
    lineheight = 0.90,
    size = 1.48 * d_internal_font_scale,
    family = figure_font_family,
    color = "#27313A"
  ) +
  geom_polygon(
    data = d_endpoint_legend_shape,
    aes(x = x, y = y, group = context, fill = context),
    inherit.aes = FALSE,
    color = "#374151",
    alpha = 0.27,
    linewidth = distribution_endpoint_linewidth
  ) +
  annotate(
    "point",
    x = c(0.17, 0.23),
    y = c(0.500, 0.500),
    shape = d_endpoint_median_shape,
    color = unname(context_fill_colors[c("in vitro", "in vivo")]),
    size = d_endpoint_median_size,
    stroke = 0.42
  ) +
  annotate(
    "text",
    x = 0.39,
    y = c(0.540, 0.460),
    label = c("in vitro", "in vivo"),
    hjust = 0,
    size = 1.25 * d_internal_font_scale,
    family = figure_font_family,
    color = "#27313A"
  ) +
  annotate(
    "text",
    x = 0.05,
    y = 0.405,
    label = d_optimizer_legend_label,
    hjust = 0,
    lineheight = 0.90,
    size = 1.48 * d_internal_font_scale,
    family = figure_font_family,
    color = "#27313A"
  ) +
  annotate(
    "text",
    x = 0.04,
    y = 0.350,
    label = d_family_outline_legend_label,
    hjust = 0,
    lineheight = 0.90,
    size = 1.72 * d_internal_font_scale,
    family = figure_font_family,
    fontface = "bold",
    color = "#27313A"
  ) +
  geom_segment(
    data = d_endpoint_legend_data,
    aes(x = 0.05, xend = 0.25, y = y, yend = y),
    inherit.aes = FALSE,
    color = unname(region_colors[d_endpoint_legend_data$family]),
    linewidth = distribution_endpoint_linewidth
  ) +
  geom_text(
    data = d_endpoint_legend_data,
    aes(x = 0.30, y = y, label = family),
    inherit.aes = FALSE,
    hjust = 0,
    size = 1.48 * d_internal_font_scale,
    family = figure_font_family,
    color = "#27313A"
  ) +
  annotate(
    "rect",
    xmin = 0.05,
    xmax = 0.25,
    ymin = 0.105,
    ymax = 0.140,
    fill = joint_bound_edge_color,
    alpha = joint_bound_edge_alpha,
    color = joint_bound_edge_color,
    linewidth = 0.20
  ) +
  annotate(
    "text",
    x = 0.30,
    y = 0.123,
    label = d_bound_legend_label,
    hjust = 0,
    lineheight = 0.90,
    size = 1.48 * d_internal_font_scale,
    family = figure_font_family,
    color = "#27313A"
  ) +
  annotate(
    "text",
    x = 0.04,
    y = 0.035,
    label = d_natural_axis_legend_label,
    hjust = 0,
    lineheight = 0.90,
    size = 1.28 * d_internal_font_scale,
    family = figure_font_family,
    color = "#4B5563"
  ) +
  scale_fill_manual(values = context_fill_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = figure_font_family) +
  theme(plot.margin = margin(1, 1, 1, 1))

p_d_content <- wrap_plots(
  wrap_elements(full = p_d_context_main, clip = FALSE),
  wrap_elements(full = p_d_compact_legend, clip = FALSE),
  ncol = 2,
  widths = unname(d_content_widths)
)
p_d_content_scaled <- plot_spacer() +
  inset_element(
    p_d_content,
    left = 0,
    bottom = 0,
    right = d_total_width_scale,
    top = 1,
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  )

p_parameter_distributions <- wrap_plots(
  wrap_elements(full = p_d_context_header, clip = FALSE),
  wrap_elements(full = p_d_content_scaled, clip = FALSE),
  ncol = 1,
  heights = c(
    d_assembled_header_fraction,
    1 - d_assembled_header_fraction
  )
)

# Archived point-only six-winner implementation. It is retained in the source
# for provenance but is not executed or used in the assembled figure.
if (FALSE) {
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

if (nrow(ratio_df) != 14L * primary_family_count ||
    length(unique(ratio_df$parameter)) != 14L ||
    any(table(ratio_df$parameter) != primary_family_count)) {
  stop("Panel B requires 14 coupled parameters with one value per family")
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
    breaks = family_levels,
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
}

# -------------------------------------------------------------------------
# Panel E: common-grid oxygen response functions.
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

validate_selected_function_solutions <- function(summary_df, label) {
  if (any(summary_df$n_solutions != nrow(selected))) {
    stop(
      label,
      " does not contain ",
      nrow(selected),
      " retained solutions at every grid point"
    )
  }
}

proliferation_df <- build_common_grid_function("proliferation_rate")
missegregation_df <- build_common_grid_function("ms_rate")
proliferation_summary <- summarise_function(proliferation_df, "oxygen_pct")
missegregation_summary <- summarise_function(missegregation_df, "oxygen_pct")
validate_selected_function_solutions(proliferation_summary, "Proliferation panel")
validate_selected_function_solutions(missegregation_summary, "Missegregation panel")

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
      breaks = pretty(y_limits, n = 3),
      expand = expansion(mult = c(0, 0.01))
    ) +
    coord_cartesian(ylim = y_limits) +
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
    ) +
    bef_emphasized_axis_theme
}

proliferation_limits <- function_limits(
  proliferation_df,
  proliferation_summary
)
missegregation_limits <- function_limits(
  missegregation_df,
  missegregation_summary
)
# Keep the Figure 5E MS column on the fixed manuscript scale requested for
# direct comparison between the 2N and 4N rows.
missegregation_limits <- c(0, 0.25)

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
  "E.  Oxygen-dependent fitted functions",
  size = panel_title_size,
  hjust = 0,
  x = 0
)
p_c_header_block <- p_c_title
p_c_shared_x_title <- ggplot() +
  annotate(
    "segment",
    x = d_shared_x_positions,
    xend = d_shared_x_positions,
    y = 1.00,
    yend = 0.84,
    linewidth = bef_axis_tick_linewidth,
    color = bef_axis_color
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
  make_vertical_strip_band("2N", size = strip_text_size),
  p_c_prolif_4n,
  p_c_ms_4n,
  make_vertical_strip_band("4N", size = strip_text_size),
  p_c_shared_x_title,
  design = p_c_grid_design,
  widths = c(0.015, 0.465, 0.015, 0.465, 0.04),
  heights = unname(c_grid_heights)
)
p_c <- p_c_plot_matrix

# -------------------------------------------------------------------------
# Panel F: post-missegregation survival across the same six winners.
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
validate_selected_function_solutions(survival_summary, "Post-MS survival panel")

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
    x = "Chromosome N",
    y = "Survival"
  ) +
  theme_manuscript(base_size = 7) +
  theme(
    legend.position = survival_legend_position,
    legend.justification = c(0.5, 0.5),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(
      family = figure_font_family,
      size = axis_text_size
    ),
    legend.key.width = grid::unit(0.42, "cm"),
    legend.key.height = grid::unit(0.20, "cm"),
    legend.spacing.x = grid::unit(0.10, "cm"),
    legend.background = element_rect(
      fill = subplot_background_fill,
      color = NA
    ),
    legend.margin = margin(0.8, 1.2, 0.8, 1.2),
    axis.title = element_text(
      family = figure_font_family,
      size = axis_title_size
    ),
    axis.text = element_text(
      family = figure_font_family,
      size = axis_text_size
    ),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5)
  ) +
  bef_emphasized_axis_theme +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = efg_median_linewidth,
        alpha = 1
      )
    ),
    linetype = "none"
  )
p_d <- p_d_body
p_survival_header <- make_text_band(
  "F.  Survival",
  size = panel_title_size,
  hjust = 0,
  x = 0,
  y = 0.58
)

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
    "figure5e_fit_adequacy",
    "figure5b_warm_start_map",
    "figure5c_survival",
    "figure5c_o2_functions",
    "figure5d_survival",
    "figure5d_o2_functions",
    "figure5e_survival",
    "figure5e_warm_start_map",
    "figure5f_parameter_ratios",
    "figure5f_solution_distributions",
    "figure5f_context_initial_optimizer_endpoints"
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
  width = 24.0,
  height = 1.85,
  dpi = 300
)
save_plot_pair(
  p_b,
  file.path(panel_root, "figure5b_model_performance"),
  width = 10.40,
  height = 2.55,
  dpi = 300
)
save_plot_pair(
  wrap_plots(
    wrap_elements(full = p_landscape_header, clip = FALSE),
    wrap_elements(full = p_landscape_body, clip = FALSE),
    ncol = 1,
    heights = c(0.10, 0.90)
  ),
  file.path(panel_root, "figure5c_warm_start_map"),
  width = 3.65,
  height = 2.55,
  dpi = 300
)
save_plot_pair(
  p_parameter_distributions,
  file.path(panel_root, "figure5d_context_initial_optimizer_endpoints"),
  width = 12.0,
  height = distribution_standalone_height_inches,
  dpi = 300
)
save_plot_pair(
  wrap_plots(
    wrap_elements(full = p_c_header_block, clip = FALSE),
    wrap_elements(full = p_c, clip = FALSE),
    ncol = 1,
    heights = c(0.10, 0.90)
  ),
  file.path(panel_root, "figure5e_o2_functions"),
  width = 5.75,
  height = 2.35,
  dpi = 300
)
save_plot_pair(
  wrap_plots(
    wrap_elements(full = p_survival_header, clip = FALSE),
    wrap_elements(full = p_d, clip = FALSE),
    ncol = 1,
    heights = c(0.08, 0.92)
  ),
  file.path(panel_root, "figure5f_survival"),
  width = 5.75,
  height = 2.35,
  dpi = 300
)

# Archived iteration2 compact layout retained for provenance only.
if (FALSE) {
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
      primary_family_count, " landscape-informed joint-pair winners. Direct fits are in-sample; ",
      "across-winner ranges and traces show solution sensitivity, not confidence intervals."
    ),
    caption = paste0(
      "C-D: thin curves = ", primary_family_count, " fits; thick = pointwise median. ",
      "Necrosis omitted: predictions unavailable.\n",
      "B-F use the same six winners; D is display-interpolated ",
      "to a declared 0-5% O2 grid."
    ),
    theme = theme(
      text = element_text(family = figure_font_family),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        family = figure_font_family,
        size = 15.5,
        face = "bold",
        color = "#111827",
        margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        family = figure_font_family,
        size = 9.8,
        color = "#4B5563",
        margin = margin(b = 2)
      ),
      plot.caption = element_text(
        family = figure_font_family,
        size = 9.4,
        color = "#5F6B76",
        hjust = 0,
        margin = margin(t = 2)
      ),
      plot.margin = margin(4, 4, 3, 4)
    )
  )
}

p_landscape_complete <- wrap_plots(
  wrap_elements(full = p_landscape_header, clip = FALSE),
  wrap_elements(full = p_landscape_body, clip = FALSE),
  ncol = 1,
  heights = c(0.060, 0.940)
) &
  theme(
    plot.background = element_rect(
      fill = subplot_background_fill,
      color = NA
    ),
    legend.background = element_rect(
      fill = subplot_background_fill,
      color = NA
    ),
    legend.box.background = element_rect(
      fill = subplot_background_fill,
      color = NA
    )
  )
p_functions_complete <- wrap_plots(
  wrap_elements(full = p_c_header_block, clip = FALSE),
  wrap_elements(full = p_c, clip = FALSE),
  ncol = 1,
  heights = c(0.10, 0.90)
)
p_survival_complete <- wrap_plots(
  wrap_elements(full = p_survival_header, clip = FALSE),
  wrap_elements(full = p_d, clip = FALSE),
  ncol = 1,
  heights = c(0.08, 0.92)
)
p_bc_row <- wrap_plots(
  wrap_elements(full = p_b, clip = FALSE),
  wrap_elements(full = p_landscape_complete, clip = FALSE),
  ncol = 2,
  widths = unname(bc_row_widths)
)
bc_d_total_height <- sum(assembled_row_heights[c("BC", "D")])
bc_d_row_boundary <- assembled_row_heights[["D"]] / bc_d_total_height
d_content_top <- bc_d_row_boundary -
  assembled_row_heights[["D"]] * d_assembled_header_fraction /
    bc_d_total_height
c_overlap_bottom <- (
  assembled_row_heights[["D"]] - c_d_overlap_height
) / bc_d_total_height
p_bc_d_overlap_block <- plot_spacer() +
  inset_element(
    p_landscape_complete,
    left = bc_row_widths[["B"]],
    bottom = c_overlap_bottom,
    right = 1,
    top = 1,
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  ) +
  inset_element(
    p_b,
    left = 0,
    bottom = bc_d_row_boundary,
    right = bc_row_widths[["B"]],
    top = 1,
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  ) +
  inset_element(
    p_d_context_header,
    left = 0,
    bottom = d_content_top,
    right = bc_row_widths[["B"]],
    top = bc_d_row_boundary,
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  ) +
  inset_element(
    p_d_content_scaled,
    left = 0,
    bottom = 0,
    right = 1,
    top = d_content_top,
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  )
p_ef_row <- wrap_plots(
  wrap_elements(full = p_functions_complete, clip = FALSE),
  wrap_elements(full = p_survival_complete, clip = FALSE),
  ncol = 2,
  widths = unname(ef_panel_widths)
)

figure5 <- (
    wrap_elements(full = p_a) /
    wrap_elements(full = p_bc_d_overlap_block) /
    wrap_elements(full = p_ef_row) +
    plot_layout(
      heights = c(
        assembled_row_heights[["A"]],
        sum(assembled_row_heights[c("BC", "D")]),
        assembled_row_heights[["EF"]]
      )
    )
) +
  plot_annotation(
    title = paste0(
      "Joint-fit workflow, performance, mechanisms, and solution stability"
    ),
    subtitle = paste0(
      primary_family_count,
      " primary-family pairs; Panel D compares pair-specific DE initial populations with 500 optimizer endpoints."
    ),
    caption = paste0(
      "B: circle = 2N start; triangle = 4N start. ",
      "E-F: thin curves = ", primary_family_count,
      " fits; thick curves = pointwise medians.\n",
      "D: upper = in vitro; lower = in vivo; gray fill/dashes = DE initialization; ",
      "pink/blue fill with C-family solid outlines = optimizer endpoints.\n",
      "Red background = outer 5% of the complete joint bound.\n",
      paste(family_levels, collapse = ", "),
      " retain the same primary-family identities across panels; optimizer distributions are descriptive ",
      "numerical-search summaries, not posterior or confidence distributions."
    ),
    theme = theme(
      text = element_text(family = figure_font_family),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        family = figure_font_family,
        size = 23.0,
        face = "bold",
        color = "#111827",
        margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        family = figure_font_family,
        size = 13.0,
        color = "#4B5563",
        margin = margin(b = 2)
      ),
      plot.caption = element_text(
        family = figure_font_family,
        size = 11.0,
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
  width = assembled_width_inches,
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
    "shared parameter-function palette",
    "all optimizer solution contrasts",
    "pair-aware optimizer summary",
    "warm-start-family optimizer summary",
    "cross-family direction summary",
    "warm-start-family common-grid density table",
    "objective-filter and deduplication sensitivity",
    "Figure 5D DE-initial/optimizer common-grid densities",
    "Figure 5D DE-initial and optimizer family summaries",
    "Figure 5D optimizer cross-family summary",
    "Figure 5D DE-initial/optimizer readiness",
    "Figure 5D direct-context common-grid densities",
    "Figure 5D direct-context family summaries",
    "Figure 5D complete joint-union bounds and edge bands",
    "Figure 5D best-winner ratio-based parameter row order",
    "Figure 5D chart contract"
  ),
  path = normalizePath(
    c(
      selection_path,
      tsne_full_path,
      tsne_best_path,
      parameter_meta_path,
      parameter_palette_path,
      optimizer_solution_path,
      pair_summary_path,
      family_summary_path,
      cross_family_summary_path,
      family_density_path,
      sensitivity_validation_path,
      ratio_density_path,
      ratio_summary_path,
      ratio_cross_family_path,
      ratio_readiness_path,
      context_density_path,
      context_summary_path,
      joint_bound_band_path,
      parameter_row_order_path,
      chart_contract_path
    ),
    mustWork = TRUE
  ),
  md5 = unname(tools::md5sum(
    c(
      selection_path,
      tsne_full_path,
      tsne_best_path,
      parameter_meta_path,
      parameter_palette_path,
      optimizer_solution_path,
      pair_summary_path,
      family_summary_path,
      cross_family_summary_path,
      family_density_path,
      sensitivity_validation_path,
      ratio_density_path,
      ratio_summary_path,
      ratio_cross_family_path,
      ratio_readiness_path,
      context_density_path,
      context_summary_path,
      joint_bound_band_path,
      parameter_row_order_path,
      chart_contract_path
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

# Archived layout-specific validation retained for provenance only.
if (FALSE) {
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
    "#C99700", "#6A3D9A", "#006D2C",
    "solid", "#C99700,#6A3D9A,#006D2C",
    "#000000", "plain", "FALSE",
    "0.33", "0.5", "3.05",
    "#FFFFFF", "#FFFFFF",
    "0.275", "0.625", "0.448", "0.56",
    "solid", "solid",
    "TRUE", "TRUE"
  ),
  stringsAsFactors = FALSE
)
}

family_dependent_parameters <- paste(
  cross_family_summary$parameter[
    !cross_family_summary$directionally_stable
  ],
  collapse = ","
)
if (FALSE) {
pooled_contraction <- prior_solution_summary[
  prior_solution_summary$distribution_role == "solution_pooled",
  "width90_relative_to_prior"
]
prior_replicate_scaled_max <- max(
  prior_sampling_config$scaled_replicate_quantile_difference
)
sensitivity_direction_stable <- unique(
  sensitivity_validation[, c(
    "parameter", "direction_stable_across_sensitivities"
  )]
)
validation <- data.frame(
  check = c(
    "selected_pair_winners",
    "optimizer_solution_rows",
    "optimizer_pair_seed_combinations",
    "optimizer_pairs",
    "optimizer_families",
    "coupled_parameters",
    "induced_prior_draw_rows",
    "prior_density_rows",
    "prior_solution_summary_rows",
    "prior_sampling_parameter_rows",
    "parameters_with_gaussian_soft_prior",
    "prior_acceptance_above_half",
    "prior_replicate_quantiles_stable",
    "all_solution_width90_below_prior",
    "sensitivity_direction_stable_parameters",
    "unique_optimizer_endpoints",
    "solutions_per_pair_parameter",
    "solutions_per_family_parameter",
    "infeasible_solution_rows",
    "projected_solution_rows",
    "nonpositive_ratio_rows",
    "directionally_stable_parameters",
    "family_dependent_parameters",
    "frequent_boundary_parameters",
    "panel_order",
    "panel_b_layout",
    "panel_b_starting_ploidy_shapes",
    "panel_d_e_context_colors",
    "cluster_color_C01",
    "cluster_color_C02",
    "cluster_color_C03",
    "cluster_identity_order",
    "panel_f_interval_semantics",
    "panel_f_pair_aware",
    "panel_f_free_x_scales",
    "panel_f_posterior_claim",
    "pooled_tsne_rows",
    "displayed_in_vivo_cluster_hulls",
    "displayed_in_vitro_cluster_hulls",
    "assembled_height_inches",
    "assembled_png",
    "assembled_pdf"
  ),
  observed = c(
    nrow(selected),
    nrow(optimizer_solutions),
    length(unique(paste(
      optimizer_solutions$pair_id,
      optimizer_solutions$seed_number,
      sep = "|"
    ))),
    length(unique(optimizer_solutions$pair_id)),
    length(unique(optimizer_solutions$family)),
    length(unique(optimizer_solutions$parameter)),
    sum(prior_sampling_config$retained_draws),
    nrow(prior_density),
    nrow(prior_solution_summary),
    nrow(prior_sampling_config),
    sum(prior_sampling_config$soft_prior_included),
    all(prior_sampling_config$primary_acceptance_rate > 0.50),
    prior_replicate_scaled_max <= 0.05,
    all(pooled_contraction > 0 & pooled_contraction < 1),
    sum(sensitivity_direction_stable$direction_stable_across_sensitivities),
    unique(sensitivity_validation$n_endpoints[
      sensitivity_validation$subset == "unique_endpoints"
    ]),
    paste(sort(unique(table(
      optimizer_solutions$pair_id,
      optimizer_solutions$parameter
    ))), collapse = ","),
    paste(sort(unique(table(
      optimizer_solutions$family,
      optimizer_solutions$parameter
    ))), collapse = ","),
    sum(!optimizer_solutions$feasible_at_solution),
    sum(optimizer_solutions$projection_applied),
    sum(optimizer_solutions$ratio_vivo_to_vitro <= 0),
    sum(cross_family_summary$directionally_stable),
    family_dependent_parameters,
    sum(cross_family_summary$max_active_bound_fraction >= 0.10),
    "A,B,C,D,E,F",
    "1x4",
    "2N=circle;4N=triangle",
    paste(context_colors, collapse = ","),
    region_colors[["C01"]],
    region_colors[["C02"]],
    region_colors[["C03"]],
    paste(family_levels, collapse = ","),
    "prior_density,pooled_solution_density,family_density,pair_medians",
    all(pair_solution_summary$n_solutions == 500L) &&
      all(family_solution_summary$n_pairs == 2L),
    TRUE,
    FALSE,
    nrow(tsne_full),
    length(unique(cluster_hulls$cluster_id)),
    sum(grepl("^vt_", unique(cluster_hulls$cluster_id))),
    assembled_height_inches,
    file.exists(paste0(output_stub, ".png")),
    file.exists(paste0(output_stub, ".pdf"))
  ),
  expected = c(
    "6",
    "42000",
    "3000",
    "6",
    "3",
    "14",
    "1400000",
    "28070",
    "70",
    "14",
    "5",
    "TRUE",
    "TRUE",
    "TRUE",
    "14",
    "947",
    "500",
    "1000",
    "0",
    "0",
    "0",
    "11",
    "p_wgd,n_O,buffer_n_exp",
    "7",
    "A,B,C,D,E,F",
    "1x4",
    "2N=circle;4N=triangle",
    "#0072B2,#CC79A7",
    "#C99700",
    "#6A3D9A",
    "#006D2C",
    "C01,C02,C03",
    "prior_density,pooled_solution_density,family_density,pair_medians",
    "TRUE",
    "TRUE",
    "FALSE",
    "228000",
    "3",
    "0",
    "10.3",
    "TRUE",
    "TRUE"
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
}

if (FALSE) {
density_parameter_groups <- split(
  posterior_density,
  posterior_density$parameter
)
common_x_range_by_parameter <- vapply(
  density_parameter_groups,
  function(group) {
    length(unique(group$display_from)) == 1L &&
      length(unique(group$display_to)) == 1L
  },
  logical(1)
)
common_bandwidth_by_parameter <- vapply(
  density_parameter_groups,
  function(group) {
    cells <- split(
      group$bandwidth,
      interaction(group$family, group$distribution_role, drop = TRUE)
    )
    all(vapply(cells, function(x) length(unique(x)) == 1L, logical(1)))
  },
  logical(1)
)
common_grid_by_parameter <- vapply(
  density_parameter_groups,
  function(group) {
    family_grids <- split(
      group$x,
      interaction(group$family, group$distribution_role, drop = TRUE)
    )
    all(vapply(
      family_grids,
      function(x) isTRUE(all.equal(x, family_grids[[1L]])),
      logical(1)
    ))
  },
  logical(1)
)
density_peak_normalization <- tapply(
  posterior_density$density_scaled_distribution,
  interaction(
    posterior_density$parameter,
    posterior_density$family,
    posterior_density$distribution_role,
    drop = TRUE
  ),
  max
)
sensitivity_direction_stable <- unique(
  sensitivity_validation[, c(
    "parameter", "direction_stable_across_sensitivities"
  )]
)
density_cell_counts <- table(
  posterior_density$parameter,
  posterior_density$family
)
density_role_cell_counts <- table(
  posterior_density$parameter,
  posterior_density$family,
  posterior_density$distribution_role
)

validation <- data.frame(
  check = c(
    "selected_pair_winners",
    "optimizer_solution_rows",
    "optimizer_pair_seed_combinations",
    "optimizer_pairs",
    "optimizer_families",
    "coupled_parameters",
    "posterior_density_rows",
    "posterior_density_parameter_family_cells",
    "density_grid_points_per_cell",
    "posterior_summary_rows",
    "posterior_selected_warmup_inputs",
    "posterior_cross_initialization_rows",
    "solutions_per_pair_parameter",
    "solutions_per_family_parameter",
    "panel_f_parameter_blocks",
    "panel_f_family_subpanels",
    "panel_f_outer_layout",
    "panel_f_inner_layout",
    "panel_f_common_x_range_within_parameter",
    "panel_f_common_bandwidth_within_parameter",
    "panel_f_common_x_grid_within_parameter",
    "panel_f_each_distribution_peak_normalized",
    "panel_f_fixed_prior_displayed",
    "sensitivity_direction_stable_parameters",
    "unique_optimizer_endpoints",
    "infeasible_solution_rows",
    "projected_solution_rows",
    "nonpositive_ratio_rows",
    "directionally_stable_parameters",
    "family_dependent_parameters",
    "frequent_boundary_parameters",
    "panel_order",
    "panel_b_layout",
    "panel_b_starting_ploidy_shapes",
    "panel_d_e_context_colors",
    "cluster_color_C01",
    "cluster_color_C02",
    "cluster_color_C03",
    "cluster_identity_order",
    "panel_f_interval_semantics",
    "panel_f_pair_aware",
    "panel_f_parameter_specific_x_scales",
    "panel_f_posterior_claim",
    "pooled_tsne_rows",
    "displayed_in_vivo_cluster_hulls",
    "displayed_in_vitro_cluster_hulls",
    "assembled_height_inches",
    "assembled_png",
    "assembled_pdf"
  ),
  observed = c(
    nrow(selected),
    nrow(optimizer_solutions),
    length(unique(paste(
      optimizer_solutions$pair_id,
      optimizer_solutions$seed_number,
      sep = "|"
    ))),
    length(unique(optimizer_solutions$pair_id)),
    length(unique(optimizer_solutions$family)),
    length(unique(optimizer_solutions$parameter)),
    nrow(posterior_density),
    sum(density_cell_counts > 0),
    paste(sort(unique(as.vector(density_role_cell_counts))), collapse = ","),
    nrow(posterior_summary),
    length(unique(posterior_summary$warmup_label)),
    nrow(posterior_cross_initialization),
    paste(sort(unique(table(
      optimizer_solutions$pair_id,
      optimizer_solutions$parameter
    ))), collapse = ","),
    paste(sort(unique(table(
      optimizer_solutions$family,
      optimizer_solutions$parameter
    ))), collapse = ","),
    length(parameter_panels),
    sum(density_cell_counts > 0),
    "2x7",
    "1x3",
    all(common_x_range_by_parameter),
    all(common_bandwidth_by_parameter),
    all(common_grid_by_parameter),
    all(abs(density_peak_normalization - 1) <= 1e-12),
    TRUE,
    sum(sensitivity_direction_stable$direction_stable_across_sensitivities),
    unique(sensitivity_validation$n_endpoints[
      sensitivity_validation$subset == "unique_endpoints"
    ]),
    sum(!optimizer_solutions$feasible_at_solution),
    sum(optimizer_solutions$projection_applied),
    sum(optimizer_solutions$ratio_vivo_to_vitro <= 0),
    sum(cross_family_summary$directionally_stable),
    family_dependent_parameters,
    sum(cross_family_summary$max_active_bound_fraction >= 0.10),
    "A,B,C,D,E,F",
    "1x4",
    "2N=circle;4N=triangle",
    paste(context_colors, collapse = ","),
    region_colors[["C01"]],
    region_colors[["C02"]],
    region_colors[["C03"]],
    paste(family_levels, collapse = ","),
    "fixed_prior,T1_generalized_posterior,optimizer_endpoints,90pct_credible_interval,optimizer_5to95pct_span,medians",
    all(tapply(
      posterior_summary$warmup_label,
      posterior_summary$family,
      function(x) length(unique(x))
    ) == 1L) && all(posterior_summary$n_chains == 2L),
    TRUE,
    TRUE,
    nrow(tsne_full),
    length(unique(cluster_hulls$cluster_id)),
    sum(grepl("^vt_", unique(cluster_hulls$cluster_id))),
    assembled_height_inches,
    file.exists(paste0(output_stub, ".png")),
    file.exists(paste0(output_stub, ".pdf"))
  ),
  expected = c(
    "6",
    "42000",
    "3000",
    "6",
    "3",
    "14",
    "50526",
    "42",
    "401",
    "42",
    "3",
    "14",
    "500",
    "1000",
    "14",
    "42",
    "2x7",
    "1x3",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "14",
    "947",
    "0",
    "0",
    "0",
    "11",
    "p_wgd,n_O,buffer_n_exp",
    "7",
    "A,B,C,D,E,F",
    "1x4",
    "2N=circle;4N=triangle",
    "#0072B2,#CC79A7",
    "#C99700",
    "#6A3D9A",
    "#006D2C",
    "C01,C02,C03",
    "fixed_prior,T1_generalized_posterior,optimizer_endpoints,90pct_credible_interval,optimizer_5to95pct_span,medians",
    "TRUE",
    "TRUE",
    "TRUE",
    "228000",
    "3",
    "0",
    "10.3",
    "TRUE",
    "TRUE"
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
}

context_density_groups <- split(context_density, context_density$parameter)
common_x_range_by_parameter <- vapply(
  context_density_groups,
  function(group) {
    length(unique(group$display_from)) == 1L &&
      length(unique(group$display_to)) == 1L
  },
  logical(1)
)
common_grid_by_parameter <- vapply(
  context_density_groups,
  function(group) {
    cell_grids <- split(
      group$x,
      interaction(
        group$family, group$context, group$distribution_role, drop = TRUE
      )
    )
    all(vapply(
      cell_grids,
      function(x) isTRUE(all.equal(x, cell_grids[[1L]])),
      logical(1)
    ))
  },
  logical(1)
)
context_density_peaks <- tapply(
  context_density$density_scaled_distribution,
  interaction(
    context_density$parameter,
    context_density$family,
    context_density$context,
    context_density$distribution_role,
    drop = TRUE
  ),
  max
)
selected_d_context <- context_summary[, c("family", "warmup_label")]
selected_d_context <- unique(selected_d_context)
panel_d_bound_band_contract <-
  nrow(joint_parameter_bounds) == 14L &&
  all(joint_bound_consistent) &&
  all(joint_parameter_bounds$joint_union_upper_transformed >
        joint_parameter_bounds$joint_union_lower_transformed) &&
  all(abs(joint_parameter_bounds$edge_fraction - 0.05) <= 1e-12) &&
  all(abs(
    (joint_parameter_bounds$lower_edge_end_transformed -
       joint_parameter_bounds$joint_union_lower_transformed) /
      (joint_parameter_bounds$joint_union_upper_transformed -
         joint_parameter_bounds$joint_union_lower_transformed) - 0.05
  ) <= 1e-12) &&
  all(abs(
    (joint_parameter_bounds$joint_union_upper_transformed -
       joint_parameter_bounds$upper_edge_start_transformed) /
      (joint_parameter_bounds$joint_union_upper_transformed -
         joint_parameter_bounds$joint_union_lower_transformed) - 0.05
  ) <= 1e-12) &&
  file.exists(joint_bound_band_path)

validation <- data.frame(
  check = c(
    "selected_joint_winners_for_panels_A_B_E_F",
    "selected_pairs_for_panel_D",
    "panel_D_selected_warmup_labels",
    "panel_D_parameters",
    "panel_D_density_rows",
    "panel_D_summary_rows",
    "panel_D_cross_family_rows",
    "panel_D_grid_points_per_distribution",
    "panel_D_common_x_range_within_parameter",
    "panel_D_common_x_grid_within_parameter",
    "panel_D_each_distribution_peak_normalized",
    "panel_D_displayed_tail_fraction_at_most_0p002",
    "panel_D_bound_band_semantics",
    "panel_D_optimizer_median_marker_shape",
    "panel_D_optimizer_median_marker_color_role",
    "panel_D_optimizer_median_marker_on_axis",
    "panel_D_initial_median_marker_drawn",
    "panel_D_parameter_row_order_rule",
    "panel_D_parameter_row_order_monotone",
    "panel_D_parameter_row_order_table",
    "panel_D_parameter_row_order",
    "panel_D_optimizer_median_marker_size",
    "panel_D_non_axis_font_scale",
    "panel_D_parameter_name_size_pt",
    "panel_D_parameter_name_description_gap",
    "panel_D_distribution_role_legend",
    "panel_D_initial_distribution_fill",
    "panel_D_posterior_claim",
    "panel_D_family_identity_order",
    "panel_D_outer_layout",
    "panel_D_shared_family_header",
    "panel_D_context_direction",
    "panel_D_context_annotation_external_strip",
    "panel_D_context_strip_width_fraction",
    "panel_D_context_and_parameter_indicator_linewidth",
    "panel_D_distribution_background",
    "panel_D_family_panel_spacing_mm",
    "panel_D_parameter_annotation_ymin",
    "panel_D_parameter_background_xmax",
    "panel_D_distribution_role_legend_geometry",
    "panel_D_distribution_legend_contexts",
    "panel_D_two_line_parameter_labels",
    "panel_A_lane_order",
    "panel_A_pair_count",
    "panel_A_search_count",
    "panel_A_height_increased",
    "panel_B_layout",
    "panel_B_legend_position",
    "panel_B_starting_ploidy_legend_two_lines",
    "panel_B_context_row_order",
    "panel_B_subtitle_smaller_than_panel_title",
    "panel_B_subtitle_centered",
    "panel_C_legend_position",
    "panel_C_legend_layout",
    "panel_C_family_key_columns",
    "panel_C_best_invitro_border_removed",
    "panel_C_initial_point_size",
    "panel_C_initial_legend_point_size",
    "panel_C_maximized_in_D_noncontent_upper_right",
    "panel_D_family_header_half_height",
    "panel_D_optimizer_black_underlay",
    "panel_D_legend_right",
    "panel_D_total_content_width_scale",
    "panel_D_parameter_text_width_fraction",
    "panel_D_legend_effective_width_fraction",
    "panel_D_legend_long_labels_wrapped",
    "panel_D_family_effective_width_scale",
    "panel_E_F_side_by_side_equal_width",
    "panel_E_F_linewidth_scale",
    "panel_E_F_individual_to_median_linewidth_ratio",
    "panel_E_MS_y_upper_limit",
    "panel_order",
    "panel_D_has_largest_row_height",
    "panel_titles_numbered_with_period",
    "panel_F_legend_inside_lower_right",
    "global_panel_title_size_pt",
    "global_axis_title_size_pt",
    "global_axis_text_size_pt",
    "assembled_width_inches",
    "assembled_height_inches",
    "assembled_width_to_height_ratio",
    "cluster_color_C01",
    "cluster_color_C02",
    "cluster_color_C03",
    "panels_B_E_F_emphasized_axes",
    "subplot_background_fill",
    "assembled_png",
    "assembled_pdf"
  ),
  observed = c(
    nrow(selected),
    nrow(selected_d_context),
    paste(selected_d_context$warmup_label[match(family_levels, selected_d_context$family)], collapse = ";"),
    length(unique(context_density$parameter)),
    nrow(context_density),
    nrow(context_summary),
    nrow(ratio_cross_family),
    paste(sort(unique(as.vector(table(
      context_density$parameter,
      context_density$family,
      context_density$context,
      context_density$distribution_role
    )))), collapse = ","),
    all(common_x_range_by_parameter),
    all(common_grid_by_parameter),
    all(abs(context_density_peaks - 1) <= 1e-12),
    max(context_density$displayed_tail_fraction) <= 0.002,
    if (panel_d_bound_band_contract) {
      "outer_5pct_of_complete_joint_union_bounds_no_distribution_quantile_spans"
    } else {
      "invalid_joint_bound_band_contract"
    },
    d_endpoint_median_shape,
    paste(
      names(context_fill_colors),
      unname(context_fill_colors),
      sep = "=", collapse = ";"
    ),
    d_endpoint_marker_y,
    FALSE,
    paste0(
      "ascending_abs_arithmetic_mean_",
      paste(family_levels, collapse = "_"),
      "_best_winner_vivo_to_vitro_ratio_minus_1"
    ),
    all(diff(best_parameter_ratio_order$distance_of_mean_ratio_from_1) >= -1e-12),
    file.exists(parameter_row_order_path),
    paste(coupled_parameters, collapse = ","),
    d_endpoint_median_size,
    d_internal_font_scale,
    d_parameter_name_size_pt,
    d_parameter_name_y - d_parameter_description_y,
    "DE initial distribution;Optimizer endpoint distribution",
    distribution_initial_fill,
    FALSE,
    paste(family_levels, collapse = ","),
    paste0(parameter_distribution_row_count, "x1"),
    paste(family_header_data$family, collapse = ","),
    "upper=in vitro;lower=in vivo",
    identical(names(panel_d_row_widths)[[1L]], "context") &&
      panel_d_row_widths[["context"]] > 0,
    panel_d_row_widths[["context"]],
    paste(d_context_marker_linewidth, d_parameter_group_indicator_linewidth, sep = ","),
    "white_except_joint_bound_edge_bands",
    d_family_panel_spacing_mm,
    d_parameter_annotation_ymin,
    d_parameter_background_xmax,
    "mirrored_context_density_silhouette_and_context_colored_optimizer_axis_markers",
    paste(context_levels, collapse = ","),
    length(parameter_description_short) == parameter_distribution_row_count,
    paste(workflow_context_boxes$context, collapse = ","),
    workflow_core_boxes$title[[1L]] == paste0(primary_family_count, " warm-start\npairs"),
    workflow_core_boxes$title[[3L]] == paste0(
      format(500L * primary_family_count, big.mark = ","),
      " joint\nsearches"
    ),
    a_height_scale > 1,
    "1x4",
    "right",
    identical(b_starting_ploidy_legend_title, "Starting\nploidy"),
    paste(p_b_panel_order, collapse = ";"),
    fit_subplot_title_size < panel_title_size,
    b_subplot_title_hjust,
    c_legend_position,
    paste0(
      "scatter_left_", c_landscape_body_widths[["scatter"]],
      ";single_column_vertical_key_right_",
      c_landscape_body_widths[["legend"]]
    ),
    c_family_key_ncol,
    identical(c_vitro_best_border_color, "#FFFFFF00") &&
      c_vitro_best_border_stroke == 0,
    c_initial_point_size,
    c_initial_legend_point_size,
    c_d_overlap_height > 0 &&
      c_overlap_bottom < assembled_row_heights[["D"]] /
        sum(assembled_row_heights[c("BC", "D")]) &&
      c_overlap_bottom >= d_content_top,
    d_family_header_height,
    d_optimizer_black_underlay,
    paste(names(d_content_widths), collapse = ","),
    d_total_width_scale,
    panel_d_row_widths[["label"]],
    d_total_width_scale * d_content_widths[["legend"]],
    all(grepl("\\n", c(
      d_optimizer_legend_label,
      d_bound_legend_label,
      d_family_outline_legend_label,
      d_natural_axis_legend_label
    ))),
    assembled_width_scale * d_total_width_scale,
    length(unique(unname(ef_panel_widths))) == 1L,
    round(efg_median_linewidth / 0.56, 2),
    round(efg_individual_linewidth / efg_median_linewidth, 2),
    missegregation_limits[[2L]],
    "A,BC,D,EF",
    assembled_row_heights[["D"]] == max(assembled_row_heights),
    all(grepl(
      "^[A-F]\\.",
      c(
        "A.  Joint-fitting workflow",
        "B.  In-sample model performance",
        "C.  Warm-start map",
        "D.  Context-specific DE initialization and optimizer endpoints",
        "E.  Oxygen-dependent fitted functions",
        "F.  Survival"
      )
    )),
    survival_legend_position[[1L]] > 0.5 &&
      survival_legend_position[[1L]] < 1 &&
      survival_legend_position[[2L]] > 0 &&
      survival_legend_position[[2L]] < 0.25,
    panel_title_size,
    axis_title_size,
    axis_text_size,
    assembled_width_inches,
    assembled_height_inches,
    round(assembled_width_inches / assembled_height_inches, 3),
    region_colors[["C01"]],
    region_colors[["C02"]],
    region_colors[["C03"]],
    paste(
      bef_axis_linewidth,
      bef_axis_tick_linewidth,
      bef_axis_tick_length_mm,
      sep = ","
    ),
    subplot_background_fill,
    file.exists(paste0(output_stub, ".png")),
    file.exists(paste0(output_stub, ".pdf"))
  ),
  expected = c(
    as.character(primary_family_count),
    as.character(primary_family_count),
    paste(
      selected_d_context$warmup_label[
        match(family_levels, selected_d_context$family)
      ],
      collapse = ";"
    ),
    "14",
    as.character(14L * primary_family_count * 2L * 2L * 401L),
    as.character(14L * primary_family_count * 2L),
    "14",
    "401",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "outer_5pct_of_complete_joint_union_bounds_no_distribution_quantile_spans",
    "8",
    "in vivo=#0072B2;in vitro=#CC79A7",
    "0",
    "FALSE",
    paste0(
      "ascending_abs_arithmetic_mean_",
      paste(family_levels, collapse = "_"),
      "_best_winner_vivo_to_vitro_ratio_minus_1"
    ),
    "TRUE",
    "TRUE",
    paste(coupled_parameters, collapse = ","),
    "1.6",
    "2.7",
    "11.3",
    "0.46",
    "DE initial distribution;Optimizer endpoint distribution",
    "#B8BEC7",
    "FALSE",
    paste(family_levels, collapse = ","),
    "14x1",
    paste(family_levels, collapse = ","),
    "upper=in vitro;lower=in vivo",
    "TRUE",
    "0.004",
    "1.86,1.86",
    "white_except_joint_bound_edge_bands",
    "1.1",
    "0",
    "1.2",
    "mirrored_context_density_silhouette_and_context_colored_optimizer_axis_markers",
    "in vitro,in vivo",
    "TRUE",
    "In vitro,In vivo,In vitro,In vivo",
    "TRUE",
    "TRUE",
    "TRUE",
    "1x4",
    "right",
    "TRUE",
    "in_vitro_growth;in_vitro_mean_N;in_vivo_burden;in_vivo_terminal_mean_N",
    "TRUE",
    "0.5",
    "right",
    "scatter_left_0.88;single_column_vertical_key_right_0.12",
    "2",
    "TRUE",
    "0.18",
    "3.75",
    "TRUE",
    "0.12",
    "FALSE",
    "distributions,legend",
    "1",
    "0.15",
    "0.08192",
    "TRUE",
    "0.8",
    "TRUE",
    "1.5",
    "0.8",
    "0.25",
    "A,BC,D,EF",
    "TRUE",
    "TRUE",
    "TRUE",
    "15.5",
    "10.2",
    "9.2",
    as.character(assembled_width_inches),
    "18",
    as.character(round(assembled_width_inches / assembled_height_inches, 3)),
    "#C99700",
    "#6A3D9A",
    "#006D2C",
    "0.65,0.5,0.9",
    "transparent",
    "TRUE",
    "TRUE"
  ),
  stringsAsFactors = FALSE
)
validation$passed <- as.character(validation$observed) == validation$expected
write.table(
  validation,
  file.path(figure5_data_dir, "figure5_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
if (any(!validation$passed)) {
  stop(
    "Figure 5 validation failed: ",
    paste(validation$check[!validation$passed], collapse = ", ")
  )
}

cat("Figure 5 generation complete.\n")
cat("Selected winners:", nrow(selected), "\n")
cat(
  "Total objective range:",
  sprintf("%.4f", min(selected$objective)),
  "to",
  sprintf("%.4f", max(selected$objective)),
  "\n"
)
cat("Panel D direct-context summary rows:", nrow(context_summary), "\n")
cat("Panel D parameter rows:", length(parameter_panels), "\n")
cat("Panel D family subpanels:", 14L * primary_family_count, "\n")
cat(
  "Cross-family optimizer-direction agreement:",
  sum(ratio_cross_family$cross_family_direction_agreement),
  "of",
  nrow(ratio_cross_family),
  "\n"
)
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
    file.path(data_dir, "parameter_function_group_palette.tsv"),
    file.path(data_dir, "figure5f_optimizer_solutions.tsv"),
    file.path(data_dir, "figure5f_pair_summary.tsv"),
    file.path(data_dir, "figure5f_family_summary.tsv"),
    file.path(data_dir, "figure5f_cross_family_summary.tsv"),
    file.path(data_dir, "figure5f_family_density.tsv"),
    file.path(data_dir, "figure5f_sensitivity_validation.tsv"),
    file.path(data_dir, "figure5f_prior_optimizer_density.tsv"),
    file.path(data_dir, "figure5f_prior_optimizer_summary.tsv"),
    file.path(data_dir, "figure5f_prior_optimizer_cross_family.tsv"),
    file.path(data_dir, "figure5f_prior_optimizer_readiness.tsv"),
    file.path(data_dir, "figure5f_context_initial_optimizer_density.tsv"),
    file.path(data_dir, "figure5f_context_initial_optimizer_summary.tsv"),
    file.path(data_dir, "figure5f_chart_contract.md")
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
