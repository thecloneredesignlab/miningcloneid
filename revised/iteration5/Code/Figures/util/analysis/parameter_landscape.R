#!/usr/bin/env Rscript

# Render Figure 4B and Supplementary Figure 4-1 from the continuous fixed-O2 analysis.
# The 500 fitted rows are optimizer-derived endpoints, not posterior samples or
# biological replicates. Figure 4B groups parameters by the O2 concentration
# at their maximum absolute Spearman association (High, Medium, then Low O2).
# Within each O2 group, positive peaks precede negative peaks, and each sign
# subgroup is ordered by decreasing maximum absolute correlation.

suppressPackageStartupMessages({
  library(data.table)
  library(ggrepel)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_path <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  } else {
    normalizePath("parameter_landscape.R", mustWork = FALSE)
  }
})
script_dir <- dirname(script_path)
data_dir <- normalizePath(
  Sys.getenv("ANALYSIS_DATA_DIR", unset = file.path(script_dir, "data")),
  mustWork = TRUE
)
figure_dir <- normalizePath(
  Sys.getenv("PLOT_OUTPUT_DIR", unset = file.path(script_dir, "figures")),
  mustWork = FALSE
)
deliverable_dir <- normalizePath(
  Sys.getenv("DELIVERABLE_OUTPUT_DIR", unset = figure_dir),
  mustWork = FALSE
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deliverable_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  association = file.path(data_dir, "continuous_ploidy_spearman_by_o2.tsv"),
  ranking = file.path(data_dir, "continuous_ploidy_parameter_ranking.tsv"),
  endpoints = file.path(data_dir, "all_parameter_fitted_endpoint_values.tsv"),
  pooled_summary = file.path(data_dir, "all_parameter_pooled_distribution_summary.tsv"),
  endpoint_ranges = file.path(data_dir, "all_parameter_log10_range_summary.tsv"),
  endpoint_density = file.path(data_dir, "all_parameter_log10_density.tsv"),
  clusters = file.path(data_dir, "invivo_best_tsne_cluster_coordinates.tsv"),
  cluster_summary = file.path(data_dir, "invivo_tsne_cluster_summary.tsv"),
  tsne_metadata = file.path(data_dir, "invivo_tsne_run_metadata.tsv")
)
missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing)) {
  stop("Missing continuous Figure 4 landscape input(s): ", paste(missing, collapse = ", "))
}

association <- fread(paths$association)
ranking <- fread(paths$ranking)
parameter_long <- fread(paths$endpoints)
pooled_summary <- fread(paths$pooled_summary)
endpoint_ranges <- fread(paths$endpoint_ranges)
endpoint_density <- fread(paths$endpoint_density)
clusters_raw <- fread(paths$clusters)[dataset == "invivo"]
cluster_summary <- fread(paths$cluster_summary)[dataset == "invivo"]
tsne_metadata <- fread(paths$tsne_metadata)

if (nrow(association) != 18L * 201L ||
    uniqueN(association$parameter) != 18L ||
    uniqueN(association$O2_pct) != 201L ||
    nrow(ranking) != 18L ||
    nrow(parameter_long) != 18L * 500L ||
    uniqueN(parameter_long$seed_number) != 500L ||
    nrow(pooled_summary) != 18L ||
    uniqueN(pooled_summary$parameter) != 18L ||
    nrow(endpoint_ranges) != 18L ||
    uniqueN(endpoint_ranges$parameter) != 18L ||
    nrow(endpoint_density) < 18L * 3L ||
    uniqueN(endpoint_density$parameter) != 18L ||
    sum(parameter_long$is_lowest_objective_fit) != 18L) {
  stop("Continuous Figure 4 inputs do not satisfy the 18 x 201 / 500-endpoint contract.")
}
if (any(abs(association$spearman_rho) > 1 + 1e-12, na.rm = TRUE)) {
  stop("Spearman rho lies outside [-1, 1].")
}
if (!identical(sort(ranking$display_order), seq_len(18L)) ||
    !identical(sort(ranking$importance_rank), seq_len(18L))) {
  stop("Parameter display and importance ranks must each be 1 through 18.")
}
peak_o2_group_levels <- c("High O2", "Medium O2", "Low O2")
peak_o2_group_palette <- c(
  "Low O2" = "#2166AC",
  "Medium O2" = "#E6A400",
  "High O2" = "#B2182B"
)
ranking_order_check <- ranking[order(display_order)]
expected_peak_o2_group <- fcase(
  ranking_order_check$O2_at_max_abs >= 0 &
    ranking_order_check$O2_at_max_abs < 1.5, "Low O2",
  ranking_order_check$O2_at_max_abs >= 1.5 &
    ranking_order_check$O2_at_max_abs < 3.5, "Medium O2",
  ranking_order_check$O2_at_max_abs >= 3.5 &
    ranking_order_check$O2_at_max_abs <= 5, "High O2",
  default = NA_character_
)
if (anyNA(expected_peak_o2_group) ||
    !identical(ranking_order_check$peak_o2_group, expected_peak_o2_group) ||
    any(diff(ranking_order_check$peak_o2_group_order) < 0) ||
    any(ranking_order_check[, diff(peak_direction_order) < 0,
      by = peak_o2_group_order]$V1) ||
    any(ranking_order_check[, diff(max_abs_rho) > 1e-12,
      by = .(peak_o2_group_order, peak_direction_order)]$V1)) {
  stop(paste0(
    "Figure 4B rows are not ordered by High-to-Low peak O2, then ",
    "positive-to-negative peak rho, then descending max-|rho|."
  ))
}

tsne_metric <- function(key_name) {
  value <- tsne_metadata[tsne_metadata$key == key_name, value]
  if (length(value) != 1L) stop("Missing t-SNE metadata key: ", key_name)
  as.character(value[[1L]])
}
if (tsne_metric("fit_scope") != "separate_in_vivo_only" ||
    as.integer(tsne_metric("retained_initial_rows")) != 127500L ||
    as.integer(tsne_metric("best_endpoint_rows")) != 500L ||
    as.integer(tsne_metric("total_tsne_rows")) != 128000L ||
    as.integer(tsne_metric("tsne_parameters")) != 18L ||
    tsne_metric("excluded_variables") != "rho_2N,sigma_burden") {
  stop("The t-SNE metadata does not match the approved in-vivo-only contract.")
}

clusters <- clusters_raw[, .(
  seed_number = as.integer(seed),
  tSNE1 = as.numeric(tSNE1),
  tSNE2 = as.numeric(tSNE2),
  cluster_id,
  cluster_num = as.integer(cluster_num),
  cluster_k = as.integer(cluster_k),
  cluster_silhouette_avg = as.numeric(cluster_silhouette_avg)
)]
if (nrow(clusters) != 500L || uniqueN(clusters$seed_number) != 500L) {
  stop("Expected 500 unique exploratory t-SNE endpoint assignments.")
}

cluster_levels <- sprintf("vi_C%02d", seq_len(nrow(cluster_summary)))
cluster_sizes <- setNames(
  cluster_summary$n_seeds[match(cluster_levels, cluster_summary$cluster_id)],
  cluster_levels
)
if (anyNA(cluster_sizes) || sum(cluster_sizes) != 500L) {
  stop("Exploratory cluster sizes do not sum to 500 endpoints.")
}

available_cluster_palette <- c(
  "#2D2D2D", "#C59B2A", "#486A8C", "#8A5A44",
  "#5E7D63", "#7B668C", "#A46776", "#858585"
)
available_cluster_shapes <- c(21, 22, 24, 23, 25, 21, 22, 24)
if (length(cluster_levels) > length(available_cluster_palette)) {
  stop("No explicit palette is defined for the selected cluster count.")
}
cluster_palette <- setNames(
  available_cluster_palette[seq_along(cluster_levels)], cluster_levels
)
cluster_shapes <- setNames(
  available_cluster_shapes[seq_along(cluster_levels)], cluster_levels
)
parameter_axis <- copy(ranking)
setorder(parameter_axis, display_order)
parameter_axis[, parameter_y := 19 - display_order]
escape_axis_html <- function(value) {
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  gsub("_", "&#95;", value, fixed = TRUE)
}
parameter_axis[, c("parameter_name_axis", "parameter_description_axis") :=
                 tstrsplit(parameter_plot_label, " | ", fixed = TRUE)]
if (anyNA(parameter_axis$parameter_name_axis) ||
    anyNA(parameter_axis$parameter_description_axis)) {
  stop("Every Figure 4B parameter label must contain a name and description.")
}
parameter_axis[, parameter_axis_label := sprintf(
  paste0(
    "<span style='font-size:15.2pt;color:#000000'><b>%s</b></span><br>",
    "<span style='font-size:11.5pt;color:#202428'><b>%s</b></span>"
  ),
  escape_axis_html(parameter_name_axis),
  escape_axis_html(parameter_description_axis)
)]

association_plot <- merge(
  association,
  parameter_axis[, .(
    parameter, parameter_y, parameter_axis_label, display_order
  )],
  by = "parameter",
  all.x = TRUE,
  sort = FALSE
)
ranking_plot <- copy(parameter_axis)

endpoint_ranges <- merge(
  endpoint_ranges,
  parameter_axis[, .(parameter, expected_parameter_y = parameter_y)],
  by = "parameter",
  all.x = TRUE,
  sort = FALSE
)
setorder(endpoint_ranges, display_order)
log_floor_raw <- unique(endpoint_ranges$log_floor_raw)
if (length(log_floor_raw) != 1L || !is.finite(log_floor_raw) ||
    log_floor_raw <= 0 ||
    any(endpoint_ranges$best_seed_number != 25L) ||
    any(abs(endpoint_ranges$parameter_y - endpoint_ranges$expected_parameter_y) > 1e-12) ||
    any(!is.finite(endpoint_ranges$initial_plot)) ||
    any(!is.finite(endpoint_ranges$best_plot)) ||
    any(!is.finite(endpoint_density$x_plot)) ||
    any(!is.finite(endpoint_density$y_density))) {
  stop("Figure 4B log10 range/density intermediates fail their plotting contract.")
}

row_separator_y <- seq(1.5, 17.5, by = 1)
peak_o2_group_counts <- setNames(
  vapply(
    peak_o2_group_levels,
    function(group) sum(ranking$peak_o2_group == group),
    integer(1L)
  ),
  peak_o2_group_levels
)
nonempty_peak_o2_group_counts <- peak_o2_group_counts[
  peak_o2_group_counts > 0
]
peak_o2_group_separator_y <- if (length(nonempty_peak_o2_group_counts) > 1L) {
  18.5 - head(cumsum(nonempty_peak_o2_group_counts), -1L)
} else {
  numeric()
}
peak_sign_transition_after <- which(
  diff(ranking_order_check$peak_o2_group_order) == 0 &
    diff(ranking_order_check$peak_direction_order) != 0
)
peak_sign_separator_y <- 18.5 - peak_sign_transition_after

peak_o2_strip_layers <- lapply(seq_len(nrow(parameter_axis)), function(i) {
  annotate(
    "rect",
    xmin = -0.22, xmax = -0.155,
    ymin = parameter_axis$parameter_y[[i]] - 0.44,
    ymax = parameter_axis$parameter_y[[i]] + 0.44,
    fill = unname(peak_o2_group_palette[[
      parameter_axis$peak_o2_group[[i]]
    ]]),
    color = NA
  )
})

p_heat <- ggplot(
  association_plot,
  aes(x = O2_pct, y = parameter_y, fill = spearman_rho)
) +
  geom_tile(width = 0.025, height = 0.88) +
  peak_o2_strip_layers +
  geom_hline(
    yintercept = row_separator_y,
    linewidth = 0.24, color = "#D0D3D6"
  ) +
  geom_hline(
    yintercept = peak_o2_group_separator_y,
    linewidth = 0.80, color = "#555B61"
  ) +
  geom_hline(
    yintercept = peak_sign_separator_y,
    linewidth = 0.48, color = "#7A8086", linetype = "22"
  ) +
  scale_x_continuous(
    name = expression(paste("Fixed ", O[2], " concentration (%)")),
    limits = c(-0.25, 5.0125), breaks = 0:5, expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = parameter_axis$parameter_y,
    labels = parameter_axis$parameter_axis_label,
    limits = c(0.5, 18.5), expand = c(0, 0)
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), oob = squish,
    na.value = "#D9D9D9", guide = "none"
  ) +
  labs(title = "Parameter-ploidy association") +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    text = element_text(face = "bold", color = "#111111"),
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 1)),
    plot.subtitle = element_text(face = "bold", size = 12, color = "#34393E", margin = margin(b = 5)),
    axis.title.x = element_text(size = 15.5, margin = margin(t = 6)),
    axis.text.x = element_text(size = 13, color = "#202020"),
    axis.text.y = ggtext::element_markdown(
      size = 11.5, color = "#202020", lineheight = 0.92,
      margin = margin(r = 3)
    ),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#707070", linewidth = 0.35),
    plot.margin = margin(4, 10, 4, 16)
  )

effect_label_offset <- 0.075
ranking_plot[, effect_label_x := max_abs_rho + effect_label_offset]
effect_x_max <- 1
if (max(ranking_plot$effect_label_x) >= effect_x_max) {
  stop("Figure 4B effect labels exceed the natural 0--1 correlation range.")
}
p_effect <- ggplot(
  ranking_plot,
  aes(y = parameter_y, x = max_abs_rho)
) +
  geom_hline(
    yintercept = row_separator_y,
    linewidth = 0.24, color = "#D0D3D6"
  ) +
  geom_hline(
    yintercept = peak_o2_group_separator_y,
    linewidth = 0.80, color = "#555B61"
  ) +
  geom_hline(
    yintercept = peak_sign_separator_y,
    linewidth = 0.48, color = "#7A8086", linetype = "22"
  ) +
  geom_segment(
    aes(x = 0, xend = max_abs_rho, yend = parameter_y),
    linewidth = 0.58, color = "#B5B8BC"
  ) +
  geom_point(
    aes(fill = peak_direction), shape = 21,
    size = 2.7, stroke = 0.45, color = "#202020"
  ) +
  geom_text(
    aes(x = effect_label_x, label = sprintf("%.2f", max_abs_rho)),
    hjust = 0,
    size = 12.2 / ggplot2::.pt, fontface = "bold", color = "#202428"
  ) +
  scale_fill_manual(
    values = c(
      "Positive peak rho" = "#B2182B",
      "Negative peak rho" = "#2166AC"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    name = expression(paste("Max |", rho, "|")),
    limits = c(0, effect_x_max), breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = NULL, breaks = parameter_axis$parameter_y,
    labels = NULL, limits = c(0.5, 18.5), expand = c(0, 0)
  ) +
  labs(title = expression(paste("Peak |", rho, "|"))) +
  coord_cartesian(clip = "on") +
  theme_bw(base_size = 11.5, base_family = "Arial") +
  theme(
    text = element_text(face = "bold", color = "#111111"),
    plot.title = element_text(face = "bold", size = 15.5, margin = margin(b = 1)),
    plot.subtitle = element_text(face = "bold", size = 12, color = "#34393E", margin = margin(b = 5)),
    axis.title.x = element_text(size = 15.5, margin = margin(t = 6)),
    axis.text.x = element_text(size = 12.5),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E3E5E8", linewidth = 0.3),
    panel.border = element_rect(color = "#707070", linewidth = 0.35),
    plot.margin = margin(4, 12, 4, 10)
  )

log_floor_plot <- log10(log_floor_raw)
max_endpoint_power <- ceiling(max(endpoint_ranges$upper_plot))
positive_endpoint_powers <- seq(
  from = log_floor_plot + 1,
  to = max_endpoint_power,
  by = 2
)
if (tail(positive_endpoint_powers, 1L) != max_endpoint_power) {
  positive_endpoint_powers <- c(positive_endpoint_powers, max_endpoint_power)
}
format_endpoint_tick <- function(power) {
  if (power == 1) return("10")
  if (power == 0) return("1")
  if (power == -1) return("0.1")
  if (power == -2) return("0.01")
  sprintf("1e%d", as.integer(power))
}
endpoint_tick_powers <- if (length(positive_endpoint_powers) > 4L) {
  positive_endpoint_powers[-1L]
} else {
  positive_endpoint_powers
}
endpoint_x_breaks <- c(log_floor_plot, endpoint_tick_powers)
endpoint_x_labels <- c(
  "0",
  vapply(endpoint_tick_powers, format_endpoint_tick, character(1L))
)
range_fill <- "#D5D8DC"
range_outline <- "#5F646A"
violin_fill <- "#56B4E9"
violin_outline <- "#0072B2"
best_seed_fill <- "#009E73"

p_prior <- ggplot() +
  geom_hline(
    yintercept = row_separator_y,
    linewidth = 0.24, color = "#D0D3D6"
  ) +
  geom_hline(
    yintercept = peak_o2_group_separator_y,
    linewidth = 0.80, color = "#555B61"
  ) +
  geom_hline(
    yintercept = peak_sign_separator_y,
    linewidth = 0.48, color = "#7A8086", linetype = "22"
  ) +
  geom_ribbon(
    data = endpoint_density,
    aes(
      x = x_plot,
      ymin = y_baseline,
      ymax = y_density,
      group = parameter
    ),
    fill = violin_fill,
    color = violin_outline,
    alpha = 0.78,
    linewidth = 0.30,
    inherit.aes = FALSE
  ) +
  geom_rect(
    data = endpoint_ranges,
    aes(
      xmin = lower_plot,
      xmax = upper_plot,
      ymin = range_ymin,
      ymax = range_ymax
    ),
    fill = range_fill,
    color = range_outline,
    alpha = 0.70,
    linewidth = 0.32,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = endpoint_ranges,
    aes(
      x = lower_plot,
      xend = upper_plot,
      y = range_ycenter,
      yend = range_ycenter
    ),
    color = range_outline,
    linewidth = 0.38,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = endpoint_ranges,
    aes(
      x = lower_plot,
      xend = lower_plot,
      y = range_ymin - 0.025,
      yend = range_ymax + 0.025
    ),
    color = range_outline,
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = endpoint_ranges,
    aes(
      x = upper_plot,
      xend = upper_plot,
      y = range_ymin - 0.025,
      yend = range_ymax + 0.025
    ),
    color = range_outline,
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = endpoint_ranges,
    aes(initial_plot, range_ycenter),
    shape = 16,
    size = 1.8,
    color = "#111111",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = endpoint_ranges,
    aes(best_plot, best_point_y),
    shape = 21,
    size = 2.2,
    stroke = 0.42,
    color = "#075B46",
    fill = best_seed_fill,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    name = "Parameter value\n(log10 scale)",
    breaks = endpoint_x_breaks,
    labels = endpoint_x_labels,
    limits = c(log_floor_plot - 0.10, max_endpoint_power + 0.15),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = parameter_axis$parameter_y,
    labels = NULL,
    limits = c(0.5, 18.5),
    expand = c(0, 0)
  ) +
  labs(title = "Endpoints and fit range") +
  theme_bw(base_size = 11.5, base_family = "Arial") +
  theme(
    text = element_text(face = "bold", color = "#111111"),
    plot.title = element_text(face = "bold", size = 15.5, margin = margin(b = 1)),
    plot.subtitle = element_text(face = "bold", size = 12, color = "#34393E", margin = margin(b = 5)),
    axis.title.x = element_text(size = 14.5, margin = margin(t = 6)),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E3E5E8", linewidth = 0.3),
    panel.border = element_rect(color = "#707070", linewidth = 0.35),
    plot.margin = margin(4, 4, 4, 10)
  )

# Solid convex hulls contain every endpoint assigned to each exploratory group;
# they are descriptive membership outlines, not confidence regions.
tsne_boundaries <- clusters[, {
  hull_index <- chull(tSNE1, tSNE2)
  closed <- c(hull_index, hull_index[[1L]])
  data.table(
    region_vertex_order = seq_along(closed),
    tSNE1 = tSNE1[closed], tSNE2 = tSNE2[closed],
    n_cluster_total = .N,
    boundary_method = "Convex hull of all assigned fitted endpoints"
  )
}, by = cluster_id]
setorder(tsne_boundaries, cluster_id, region_vertex_order)
fwrite(tsne_boundaries, file.path(data_dir, "tsne_assignment_boundaries.tsv"), sep = "\t")

tsne_centers <- clusters[, .(
  tSNE1 = median(tSNE1), tSNE2 = median(tSNE2), n = .N
), by = cluster_id]
tsne_x_range <- range(clusters$tSNE1, finite = TRUE)
tsne_y_range <- range(clusters$tSNE2, finite = TRUE)
tsne_span <- max(diff(tsne_x_range), diff(tsne_y_range)) * 1.06
tsne_x_mid <- mean(tsne_x_range)
tsne_y_mid <- mean(tsne_y_range)
tsne_x_limits <- tsne_x_mid + c(-0.5, 0.5) * tsne_span
tsne_y_limits <- tsne_y_mid + c(-0.5, 0.5) * tsne_span
tsne_x_breaks <- pretty(tsne_x_limits, n = 3)
tsne_y_breaks <- pretty(tsne_y_limits, n = 3)
tsne_global <- c(tSNE1 = median(clusters$tSNE1), tSNE2 = median(clusters$tSNE2))
tsne_labels <- copy(tsne_centers)
tsne_labels[, `:=`(
  dx = tSNE1 - tsne_global[["tSNE1"]],
  dy = tSNE2 - tsne_global[["tSNE2"]
])]
tsne_labels[, norm := sqrt(dx^2 + dy^2)]
if (any(tsne_labels$norm < 1e-12)) {
  zero <- which(tsne_labels$norm < 1e-12)
  angles <- seq(
    0, 2 * pi, length.out = length(zero) + 1L
  )[seq_len(length(zero))]
  tsne_labels[zero, `:=`(dx = cos(angles), dy = sin(angles), norm = 1)]
}
tsne_labels[, `:=`(
  nudge_x = dx / norm * tsne_span * 0.10,
  nudge_y = dy / norm * tsne_span * 0.10,
  label = sprintf("ivv%s\nn=%d", sub("^vi_", "", cluster_id), n)
)]

p_tsne <- ggplot(clusters, aes(tSNE1, tSNE2)) +
  geom_path(
    data = tsne_boundaries,
    aes(tSNE1, tSNE2, group = cluster_id, color = cluster_id),
    inherit.aes = FALSE, linewidth = 0.65, alpha = 0.82
  ) +
  geom_point(
    aes(fill = cluster_id, shape = cluster_id, color = cluster_id),
    size = 1.2, alpha = 0.58, stroke = 0.2
  ) +
  geom_label_repel(
    data = tsne_labels,
    aes(tSNE1, tSNE2, label = label),
    inherit.aes = FALSE,
    nudge_x = tsne_labels$nudge_x, nudge_y = tsne_labels$nudge_y,
    seed = 44202, max.time = 2, max.iter = 50000,
    max.overlaps = Inf, force = 10, force_pull = 0.08,
    box.padding = 0.55, point.padding = 0.15,
    min.segment.length = 0,
    size = 12.2 / ggplot2::.pt, fontface = "bold",
    label.size = 0.22, fill = scales::alpha("white", 0.90),
    color = "#202020", segment.color = "#5E6368", segment.size = 0.35
  ) +
  scale_color_manual(values = cluster_palette, guide = "none") +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  scale_shape_manual(values = cluster_shapes, guide = "none") +
  scale_x_continuous(
    name = "t-SNE 1", breaks = tsne_x_breaks,
    limits = tsne_x_limits, expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "t-SNE 2", breaks = tsne_y_breaks,
    limits = tsne_y_limits, expand = c(0, 0)
  ) +
  labs(title = NULL) +
  coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  theme_bw(base_size = 12.5, base_family = "Arial") +
  theme(
    text = element_text(face = "bold", color = "#111111"),
    plot.title = element_blank(),
    plot.subtitle = element_text(face = "bold", size = 13, color = "#34393E", margin = margin(b = 5)),
    axis.title = element_text(size = 16, color = "#202020"),
    axis.text = element_text(size = 13.5, color = "#202020"),
    axis.line = element_line(
      color = "#202020", linewidth = 0.60, lineend = "square"
    ),
    axis.ticks = element_line(color = "#202020", linewidth = 0.60),
    panel.grid = element_blank(),
    panel.border = element_rect(
      fill = NA, color = "#AEB6C0", linewidth = 0.35
    ),
    plot.margin = margin(10, 10, 10, 10)
  )

# Compact vertical sidebar for Figure 4B. Cluster identities are labeled
# directly in the standalone Figure 4C and are therefore not repeated here.
legend_title_size <- 13
legend_text_size <- 10.6
peak_o2_group_labels <- c(
  "Low O2" = "Low O2\n[0, 1.5%)",
  "Medium O2" = "Medium O2\n[1.5, 3.5%)",
  "High O2" = "High O2\n[3.5, 5%]"
)
peak_o2_legend <- data.table(
  peak_o2_group = peak_o2_group_levels,
  label = unname(peak_o2_group_labels[peak_o2_group_levels]),
  color = unname(peak_o2_group_palette[peak_o2_group_levels]),
  y = seq(0.90, 0.68, length.out = length(peak_o2_group_levels))
)
rho_steps <- data.table(
  rho = seq(-1, 1, length.out = 201),
  y = seq(0.33, 0.47, length.out = 201)
)
rho_pal <- scales::gradient_n_pal(c("#2166AC", "#F7F7F7", "#B2182B"))
rho_steps[, color := rho_pal(rescale(rho, to = c(0, 1)))]
rho_ticks <- data.table(
  y = seq(0.33, 0.47, length.out = 3),
  label = c("-1", "0", "+1")
)

p_sidebar <- ggplot() +
  annotate(
    "text", x = 0.02, y = 0.985, label = "Peak O2 group",
    hjust = 0, vjust = 1, family = "Arial", fontface = "bold",
    size = legend_title_size / ggplot2::.pt, color = "#202020"
  ) +
  geom_rect(
    data = peak_o2_legend,
    aes(
      xmin = 0.05, xmax = 0.12,
      ymin = y - 0.034, ymax = y + 0.034,
      fill = color
    ),
    color = NA, inherit.aes = FALSE
  ) +
  geom_text(
    data = peak_o2_legend,
    aes(x = 0.16, y = y, label = label),
    hjust = 0, vjust = 0.5, family = "Arial",
    size = legend_text_size / ggplot2::.pt, fontface = "bold",
    lineheight = 0.92, color = "#25292D"
  ) +
  annotate(
    "text", x = 0.02, y = 0.52, label = "Spearman rho",
    hjust = 0, vjust = 1, family = "Arial", fontface = "bold",
    size = legend_title_size / ggplot2::.pt, color = "#202020"
  ) +
  geom_rect(
    data = rho_steps,
    aes(xmin = 0.07, xmax = 0.14, ymin = y - 0.00039, ymax = y + 0.00039, fill = color),
    color = NA, inherit.aes = FALSE
  ) +
  annotate(
    "rect", xmin = 0.07, xmax = 0.14, ymin = 0.33, ymax = 0.47,
    fill = NA, color = "#656A70", linewidth = 0.3
  ) +
  geom_segment(
    data = rho_ticks,
    aes(x = 0.14, xend = 0.18, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 0.35, color = "#555A60"
  ) +
  geom_text(
    data = rho_ticks,
    aes(x = 0.21, y = y, label = label),
    hjust = 0, vjust = 0.5, family = "Arial",
    size = legend_text_size / ggplot2::.pt, fontface = "bold", color = "#25292D"
  ) +
  annotate(
    "text", x = 0.02, y = 0.315, label = "Endpoint /\nrange key",
    hjust = 0, vjust = 1, family = "Arial", fontface = "bold",
    size = legend_title_size / ggplot2::.pt, color = "#202020"
  ) +
  annotate(
    "polygon",
    x = c(0.07, 0.17, 0.29, 0.41, 0.53, 0.65, 0.72, 0.72, 0.07),
    y = c(0.205, 0.225, 0.250, 0.260, 0.245, 0.220, 0.205, 0.205, 0.205),
    fill = violin_fill, color = violin_outline, linewidth = 0.30, alpha = 0.78
  ) +
  annotate(
    "point", x = 0.51, y = 0.232, shape = 21, size = 2.8,
    fill = best_seed_fill, color = "#075B46", stroke = 0.42
  ) +
  annotate(
    "text", x = 0.07, y = 0.185,
    label = "500 endpoints\n+ seed25",
    hjust = 0, vjust = 1, family = "Arial",
    size = 10.6 / ggplot2::.pt, fontface = "bold",
    lineheight = 0.96, color = "#30353A"
  ) +
  annotate(
    "rect", xmin = 0.07, xmax = 0.72, ymin = 0.105, ymax = 0.135,
    fill = range_fill, color = range_outline, linewidth = 0.32, alpha = 0.88
  ) +
  annotate(
    "segment", x = 0.07, xend = 0.72, y = 0.120, yend = 0.120,
    linewidth = 0.38, color = range_outline
  ) +
  annotate(
    "point", x = 0.48, y = 0.120, shape = 16, size = 2.8,
    color = "#111111"
  ) +
  annotate(
    "text", x = 0.07, y = 0.085,
    label = "fit range + initial",
    hjust = 0, vjust = 1, family = "Arial",
    size = 10.6 / ggplot2::.pt, fontface = "bold", color = "#30353A"
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(
    plot.margin = margin(2, 3, 4, 5),
    plot.background = element_rect(fill = "white", color = NA)
  )

# The lite sidebar retains the peak-O2 grouping and Spearman scale while
# removing the endpoint/range key together with the endpoint distribution.
p_sidebar_lite <- p_sidebar +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.32,
    fill = "white", color = NA
  )

combined_design <- c(
  area(t = 1, l = 1, b = 1, r = 11),
  area(t = 1, l = 12, b = 1, r = 15),
  area(t = 1, l = 16, b = 1, r = 21),
  area(t = 1, l = 22, b = 1, r = 27)
)
combined_core <- p_heat + p_effect + p_prior + p_sidebar +
  plot_layout(design = combined_design, widths = rep(1, 27))

combined <- combined_core

lite_design <- c(
  area(t = 1, l = 1, b = 1, r = 11),
  area(t = 1, l = 12, b = 1, r = 15),
  area(t = 1, l = 16, b = 1, r = 21)
)
combined_lite <- p_heat + p_effect + p_sidebar_lite +
  plot_layout(design = lite_design, widths = rep(1, 21))

main_base <- file.path(figure_dir, "parameter_continuous_ploidy_landscape")
main_scale_boost <- 1.10
main_width <- 15 / main_scale_boost
main_height <- 12.90 / main_scale_boost
main_raster_dpi <- 300 * main_scale_boost
ggsave(paste0(main_base, ".png"), combined, width = main_width, height = main_height,
       units = "in", dpi = main_raster_dpi, bg = "white")
ggsave(paste0(main_base, ".pdf"), combined, width = main_width, height = main_height,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(paste0(main_base, ".svg"), combined, width = main_width, height = main_height,
       units = "in", device = svglite::svglite, bg = "white")

lite_base <- file.path(figure_dir, "Figure4B_lite")
lite_width <- 12
lite_height <- 9
ggsave(paste0(lite_base, ".png"), combined_lite, width = lite_width, height = lite_height,
       units = "in", dpi = main_raster_dpi, bg = "white")
ggsave(paste0(lite_base, ".pdf"), combined_lite, width = lite_width, height = lite_height,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(paste0(lite_base, ".svg"), combined_lite, width = lite_width, height = lite_height,
       units = "in", device = svglite::svglite, bg = "white")

for (extension in c("png", "pdf", "svg")) {
  source_paths <- c(
    Figure4B = paste0(main_base, ".", extension),
    Figure4B_lite = paste0(lite_base, ".", extension)
  )
  destination_paths <- file.path(
    deliverable_dir,
    paste0(names(source_paths), ".", extension)
  )
  copied <- file.copy(
    unname(source_paths), destination_paths,
    overwrite = TRUE, copy.mode = TRUE
  )
  if (!all(copied)) {
    stop("Failed to publish Figure 4B deliverable(s): ", extension)
  }
}

# Figure 4C: the in-vivo-only fitted-endpoint t-SNE is separated from the
# association/search-range panel so its square coordinate system remains
# legible. Cluster identity and size are labeled directly inside the panel.
tsne_base <- file.path(figure_dir, "parameter_tsne_groups")
tsne_scale_boost <- 1.15
tsne_width <- 6.8 / tsne_scale_boost
tsne_height <- 6.8 / tsne_scale_boost
tsne_raster_dpi <- 300 * tsne_scale_boost
ggsave(paste0(tsne_base, ".png"), p_tsne, width = tsne_width, height = tsne_height,
       units = "in", dpi = tsne_raster_dpi, bg = "white")
ggsave(paste0(tsne_base, ".pdf"), p_tsne, width = tsne_width, height = tsne_height,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(paste0(tsne_base, ".svg"), p_tsne, width = tsne_width, height = tsne_height,
       units = "in", device = svglite::svglite, bg = "white")

# Supplementary Figure 4-1: all 18 endpoint distributions together, in the same display
# order as Figure 4B.  This is deliberately labeled as an optimizer-endpoint
# distribution rather than a posterior distribution.
si_axis <- ranking[order(display_order), .(
  parameter, parameter_plot_label, display_order
)]
wrap_si_parameter_label <- function(label) {
  parts <- strsplit(label, " | ", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) return(label)
  paste0(
    parts[[1L]], "\n",
    paste(strwrap(parts[[2L]], width = 29), collapse = "\n")
  )
}
si_axis[, parameter_facet_label := vapply(
  parameter_plot_label, wrap_si_parameter_label, character(1L)
)]
parameter_long <- merge(
  parameter_long,
  si_axis,
  by = "parameter",
  all.x = TRUE,
  suffixes = c("", ".axis"),
  sort = FALSE
)
parameter_long[, parameter_facet := factor(
  parameter_facet_label,
  levels = si_axis$parameter_facet_label
)]
cluster_factor_labels <- sprintf(
  "%s\n(n=%d)", sub("^vi_", "", cluster_levels), cluster_sizes[cluster_levels]
)
parameter_long[, cluster_factor := factor(
  cluster_id, levels = cluster_levels, labels = cluster_factor_labels
)]
figure4d_cluster_factor_labels <- sprintf(
  "ivv%s\n(n=%d)", sub("^vi_", "", cluster_levels), cluster_sizes[cluster_levels]
)

# Figure 4D uses the same 500 fitted endpoints and prior-referenced parameter
# coordinate as Supplementary Figure 4-1. The omnibus test is deliberately treated as an
# exploratory fitted-landscape separation diagnostic: the same 18 parameters
# contributed to the t-SNE and cluster assignment, so these p values are not
# independent validation of biological subtypes.
if (nrow(parameter_long) != 18L * 500L ||
    parameter_long[, uniqueN(seed_number)] != 500L ||
    parameter_long[, uniqueN(parameter)] != 18L ||
    parameter_long[, anyDuplicated(.SD), .SDcols = c("seed_number", "parameter")]) {
  stop("Figure 4D requires one complete record per seed and parameter.")
}

cluster_group_summary <- parameter_long[, .(
  n_endpoints = uniqueN(seed_number),
  median_prior_referenced_z = median(prior_referenced_z),
  mean_prior_referenced_z = mean(prior_referenced_z),
  q25_prior_referenced_z = quantile(prior_referenced_z, 0.25, names = FALSE),
  q75_prior_referenced_z = quantile(prior_referenced_z, 0.75, names = FALSE)
), by = .(parameter, cluster_id)]

cluster_parameter_tests <- parameter_long[, {
  test <- kruskal.test(prior_referenced_z ~ cluster_id)
  n_total <- .N
  n_groups <- uniqueN(cluster_id)
  statistic <- unname(test$statistic)
  .(
    n_total = n_total,
    n_groups = n_groups,
    kruskal_wallis_H = statistic,
    degrees_of_freedom = unname(test$parameter),
    p_value = test$p.value,
    epsilon_squared = max(0, (statistic - n_groups + 1) / (n_total - n_groups))
  )
}, by = .(
  parameter, parameter_plot_label, parameter_display_label, display_order
)]
cluster_parameter_tests[, bh_q_value := p.adjust(p_value, method = "BH")]
cluster_parameter_tests[, significant_bh_0p05 := bh_q_value < 0.05]

cluster_group_wide <- dcast(
  cluster_group_summary,
  parameter ~ cluster_id,
  value.var = c(
    "n_endpoints", "median_prior_referenced_z", "mean_prior_referenced_z",
    "q25_prior_referenced_z", "q75_prior_referenced_z"
  )
)
cluster_parameter_tests <- merge(
  cluster_parameter_tests,
  cluster_group_wide,
  by = "parameter",
  all.x = TRUE,
  sort = FALSE
)
setorder(
  cluster_parameter_tests,
  -epsilon_squared, bh_q_value, display_order
)
cluster_parameter_tests[, overall_effect_rank := seq_len(.N)]
cluster_parameter_tests[, significant_effect_rank := NA_integer_]
cluster_parameter_tests[significant_bh_0p05 == TRUE,
                        significant_effect_rank := seq_len(.N)]

significant_tests <- cluster_parameter_tests[
  significant_bh_0p05 == TRUE
][order(-epsilon_squared, bh_q_value, display_order)]
if (!nrow(significant_tests)) {
  stop("No parameter passes BH q < 0.05 for Figure 4D.")
}
max_significant_epsilon_squared <- max(significant_tests$epsilon_squared)
strongest_candidates <- significant_tests[
  epsilon_squared == max_significant_epsilon_squared
]
setorder(strongest_candidates, bh_q_value, display_order)
figure4d_strongest <- copy(strongest_candidates[1L])
figure4d_strongest[, `:=`(
  selected_effect_rank = 1L,
  selection_rule = paste0(
    "BH q < 0.05 eligibility gate; select the maximum Kruskal-Wallis ",
    "epsilon-squared; BH q and Figure 4B display order only break exact ",
    "epsilon-squared ties"
  )
)]
cluster_parameter_tests[, selected_for_figure4d :=
                          parameter == figure4d_strongest$parameter[[1L]]]

fwrite(
  cluster_parameter_tests,
  file.path(data_dir, "exploratory_cluster_parameter_omnibus_tests.tsv"),
  sep = "\t"
)
fwrite(
  figure4d_strongest,
  file.path(data_dir, "figure4d_strongest_parameter.tsv"),
  sep = "\t"
)

wrap_figure4d_parameter_label <- function(label) {
  parts <- strsplit(label, " | ", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) return(label)
  paste0(
    parts[[1L]], "\n",
    paste(strwrap(parts[[2L]], width = 28), collapse = "\n")
  )
}
figure4d_parameter_title <- wrap_figure4d_parameter_label(
  figure4d_strongest$parameter_plot_label[[1L]]
)
figure4d_data <- parameter_long[
  parameter == figure4d_strongest$parameter[[1L]]
]
setorder(figure4d_data, cluster_id, seed_number)
figure4d_data[, figure4d_cluster_factor := factor(
  cluster_id,
  levels = cluster_levels,
  labels = figure4d_cluster_factor_labels
)]

figure4d_plot <- ggplot(
  figure4d_data,
  aes(
    x = figure4d_cluster_factor,
    y = prior_referenced_z,
    fill = figure4d_cluster_factor
  )
) +
  geom_hline(yintercept = 0, color = "#5F646A", linewidth = 0.38) +
  geom_violin(
    trim = TRUE, scale = "width", adjust = 0.65,
    color = "#3F4348", linewidth = 0.34, alpha = 0.88
  ) +
  geom_boxplot(
    width = 0.24, outlier.shape = NA,
    color = "#272A2E", fill = "white",
    alpha = 0.62, linewidth = 0.36
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 23, size = 2.2,
    stroke = 0.45, fill = "white", color = "#111111"
  ) +
  scale_fill_manual(
    values = setNames(
      unname(cluster_palette), figure4d_cluster_factor_labels
    ),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(-1.78, 1.78), breaks = c(-1.5, 0, 1.5)) +
  labs(
    title = figure4d_parameter_title,
    x = "Exploratory t-SNE cluster",
    y = "Prior-referenced fitted parameter value"
  ) +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(
    text = element_text(face = "bold", color = "#111111"),
    plot.title = element_text(
      size = 16.5, face = "bold", color = "#202428",
      hjust = 0.5, lineheight = 0.96, margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 12.5, face = "bold", color = "#30353A", margin = margin(b = 7)
    ),
    plot.caption = element_text(
      size = 11.2, face = "bold", color = "#3E4348", hjust = 0,
      lineheight = 1.05, margin = margin(t = 7)
    ),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12.5, color = "#202020"),
    axis.line = element_line(
      color = "#202020", linewidth = 0.60, lineend = "square"
    ),
    axis.ticks = element_line(color = "#202020", linewidth = 0.60),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E3E5E8", linewidth = 0.3),
    panel.border = element_rect(
      fill = NA, color = "#AEB6C0", linewidth = 0.35
    ),
    plot.title.position = "plot",
    plot.margin = margin(10, 10, 10, 10)
  )

figure4d_base <- file.path(figure_dir, "strongest_cluster_parameter_distribution")
figure4d_scale_boost <- 1.15
figure4d_width <- 6.8 / figure4d_scale_boost
figure4d_height <- 9.95 / figure4d_scale_boost
figure4d_raster_dpi <- 300 * figure4d_scale_boost
ggsave(
  paste0(figure4d_base, ".png"), figure4d_plot,
  width = figure4d_width, height = figure4d_height,
  units = "in", dpi = figure4d_raster_dpi, bg = "white"
)
ggsave(
  paste0(figure4d_base, ".pdf"), figure4d_plot,
  width = figure4d_width, height = figure4d_height,
  units = "in", device = cairo_pdf, bg = "white"
)
ggsave(
  paste0(figure4d_base, ".svg"), figure4d_plot,
  width = figure4d_width, height = figure4d_height,
  units = "in", device = svglite::svglite, bg = "white"
)

supp_violin <- ggplot(
  parameter_long,
  aes(x = cluster_factor, y = prior_referenced_z, fill = cluster_factor)
) +
  geom_hline(yintercept = 0, color = "#5F646A", linewidth = 0.35) +
  geom_violin(
    trim = TRUE, scale = "width", adjust = 0.65,
    color = "#3F4348", linewidth = 0.3, alpha = 0.88
  ) +
  geom_boxplot(
    width = 0.24, outlier.shape = NA,
    color = "#272A2E", fill = "white",
    alpha = 0.62, linewidth = 0.32
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 23, size = 2.0,
    stroke = 0.4, fill = "white", color = "#111111"
  ) +
  facet_wrap(vars(parameter_facet), ncol = 6, nrow = 3, scales = "fixed") +
  scale_fill_manual(
    values = setNames(unname(cluster_palette), cluster_factor_labels),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(-1.78, 1.78), breaks = c(-1.5, 0, 1.5)) +
  labs(
    title = "All 18 fitted-parameter endpoint distributions across exploratory in vivo t-SNE clusters",
    subtitle = paste0(
      "Parameters follow the Figure 4B High-, Medium-, then Low-O2 peak grouping; positive precedes negative peak rho within each group, with descending max-|rho| within each sign subgroup. ",
      "All 18 parameters were inputs to the in-vivo-only t-SNE."
    ),
    x = "Exploratory t-SNE cluster",
    y = "Prior-referenced fitted parameter value",
    caption = paste0(
      "Zero is the midpoint of the configured transformed parameter interval; one unit is the SD of a uniform distribution over that interval. ",
      "Diamonds are means. These are optimizer-derived fitted endpoints, not posterior samples, credible intervals, or biological replicates."
    )
  ) +
  theme_bw(base_size = 9.2, base_family = "Arial") +
  theme(
    text = element_text(face = "bold", color = "#111111"),
    plot.title = element_text(face = "bold", size = 14, color = "#202020"),
    plot.subtitle = element_text(face = "bold", size = 9, color = "#30353A", margin = margin(b = 7)),
    plot.caption = element_text(face = "bold", size = 7.8, color = "#3E4348", hjust = 0, margin = margin(t = 7)),
    axis.title = element_text(size = 9.2),
    axis.text = element_text(size = 7.8, color = "#202020"),
    strip.text = element_text(size = 7.8, face = "bold", lineheight = 0.95),
    strip.background = element_rect(fill = "#F2F3F5", color = "#BFC3C8", linewidth = 0.35),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E3E5E8", linewidth = 0.3),
    panel.border = element_rect(color = "#70757A", linewidth = 0.35),
    plot.title.position = "plot", plot.margin = margin(12, 12, 12, 12)
  )

supp_base <- file.path(figure_dir, "all_parameter_fitted_endpoint_distributions")
supp_width <- 15
supp_height <- 7.875
supp_raster_dpi <- 400
ggsave(paste0(supp_base, ".png"), supp_violin, width = supp_width, height = supp_height,
       units = "in", dpi = supp_raster_dpi, bg = "white")
ggsave(paste0(supp_base, ".pdf"), supp_violin, width = supp_width, height = supp_height,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(paste0(supp_base, ".svg"), supp_violin, width = supp_width, height = supp_height,
       units = "in", device = svglite::svglite, bg = "white")

validation <- data.table(
  metric = c(
    "n_parameters", "n_fixed_o2_values", "n_fitted_endpoints",
    "association_metric", "continuous_ploidy_outcome",
    "binary_ploidy_class_used", "parameter_sort_primary",
    "parameter_sort_secondary", "parameter_sort_tertiary",
    "parameter_sort_tie_break", "high_o2_peak_parameter_count",
    "medium_o2_peak_parameter_count", "low_o2_peak_parameter_count",
    "row_annotation_field",
    "effect_position_field",
    "effect_fill_field", "all_parameters_in_main_distribution_plot",
    "all_parameters_in_si_distribution_plot", "main_distribution_grouping",
    "pooled_endpoints_per_parameter", "distribution_upper_glyph",
    "distribution_lower_glyph", "distribution_point",
    "distribution_point_seed", "endpoint_display_scale",
    "endpoint_zero_floor_raw", "endpoint_distributions_are_posterior",
    "endpoint_distributions_are_measurement_bounds",
    "fitting_range_fill", "tsne_boundary_method",
    "tsne_boundaries_are_confidence_regions", "figure4c_rendered",
    "figure4c_point_count", "figure4c_direct_cluster_label_count",
    "cluster_parameter_omnibus_test", "multiple_testing_adjustment",
    "cluster_parameter_effect_size", "significant_parameter_count_bh_0p05",
    "figure4d_selected_n", "figure4d_selection_rule",
    "figure4d_selected_parameter", "figure4d_rendered",
    "main_output_width_in", "main_output_height_in",
    "figure4b_lite_endpoint_distribution_rendered",
    "figure4b_lite_output_width_in", "figure4b_lite_output_height_in",
    "figure4b_lite_output_aspect_ratio",
    "figure4c_output_width_in", "figure4c_output_height_in",
    "figure4d_output_width_in", "figure4d_output_height_in"
  ),
  value = c(
    18, 201, 500, "Spearman rho", "TRUE", "FALSE",
    "peak O2 group: High [3.5,5], Medium [1.5,3.5), Low [0,1.5)",
    "peak rho sign within peak O2 group: positive first, negative second",
    "descending maximum absolute Spearman rho within peak O2/sign subgroup",
    "configured parameter order",
    unname(peak_o2_group_counts[["High O2"]]),
    unname(peak_o2_group_counts[["Medium O2"]]),
    unname(peak_o2_group_counts[["Low O2"]]),
    "peak_o2_group",
    "max_abs_rho", "peak_direction",
    "TRUE", "TRUE", "all 500 endpoints pooled", 500,
    "blue half-violin density of 500 optimizer endpoints plus green seed25 point",
    "gray configured lower-to-upper range plus black initial-value point",
    "green solid circle", 25,
    "original parameter value on shared log10 axis",
    format(log_floor_raw, scientific = TRUE),
    "FALSE", "FALSE", "gray", "solid convex hull of all assigned endpoints",
    "FALSE", "TRUE", nrow(clusters), nrow(tsne_labels),
    paste0(
      "Kruskal-Wallis rank-sum test across ",
      length(cluster_levels), " exploratory clusters"
    ),
    "Benjamini-Hochberg across 18 parameters",
    "epsilon-squared = (H - k + 1) / (n - k)",
    sum(cluster_parameter_tests$significant_bh_0p05), 1,
    figure4d_strongest$selection_rule[[1L]],
    figure4d_strongest$parameter[[1L]], "TRUE",
    main_width, main_height, "FALSE", lite_width, lite_height,
    lite_width / lite_height,
    tsne_width, tsne_height,
    figure4d_width, figure4d_height
  )
)
fwrite(validation, file.path(data_dir, "parameter_landscape_layout_validation.tsv"), sep = "\t")

source_paths <- c(
  unlist(paths),
  tsne_boundaries = file.path(data_dir, "tsne_assignment_boundaries.tsv"),
  script = script_path
)
provenance <- data.table(
  source = names(source_paths),
  path = unname(source_paths),
  md5 = unname(tools::md5sum(source_paths)),
  role = c(rep("input", length(paths)), "derived output", "script")
)
fwrite(provenance, file.path(data_dir, "parameter_landscape_source_provenance.tsv"), sep = "\t")

message("Completed continuous Figure 4 parameter landscape workflow.")
message("  Figure 4B: ", paste0(main_base, ".png"))
message("  Figure 4B lite: ", paste0(lite_base, ".png"))
message("  Figure 4C: ", paste0(tsne_base, ".png"))
message("  Figure 4D: ", paste0(figure4d_base, ".png"))
message("  Supplementary Figure 4-1: ", paste0(supp_base, ".png"))
