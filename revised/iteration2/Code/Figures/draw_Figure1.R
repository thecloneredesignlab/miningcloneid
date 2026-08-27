#!/usr/bin/env Rscript

if (identical(Sys.getenv("FIGURE1_DRAW_WORKER"), "1")) {

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
script_dir <- if (length(script_arg) == 1L) {
  normalizePath(dirname(sub("^--file=", "", script_arg)))
} else {
  normalizePath(getwd())
}
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
data_root <- normalizePath(
  Sys.getenv(
    "FIGURE1_DATA_DIR",
    unset = file.path(script_dir, "data", "figure1_frozen_inputs")
  ),
  mustWork = TRUE
)
figures_root <- normalizePath(
  Sys.getenv(
    "FIGURE1_BUILD_DIR",
    unset = file.path(script_dir, "figures")
  ),
  mustWork = FALSE
)
subpanel_root <- file.path(figures_root, "figure1_subpanels")
revised_figures_root <- normalizePath(
  Sys.getenv(
    "FIGURE_OUTPUT_DIR",
    unset = file.path(workspace_root, "Figures")
  ),
  mustWork = FALSE
)

paths <- list(
  timeline = file.path(data_root, "invitro_lineage_timeline.tsv"),
  passage = file.path(data_root, "invitro_passage_observations.tsv"),
  population = file.path(data_root, "invitro_population_measurement_events.tsv"),
  kary = file.path(data_root, "invitro_kary_cells.tsv"),
  flow_density = file.path(data_root, "invitro_flow_observed_density.tsv"),
  flow_peak = file.path(data_root, "invitro_flow_peak_ploidy.tsv"),
  burden = file.path(data_root, "invivo_burden_long.tsv"),
  invivo_ploidy = file.path(data_root, "invivo_ploidy_cells.tsv"),
  harvest = file.path(data_root, "invivo_harvest_catalog.tsv")
)
stopifnot(all(file.exists(unlist(paths))))

read_tsv <- function(path) {
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

timeline <- read_tsv(paths$timeline)
passage <- read_tsv(paths$passage)
population <- read_tsv(paths$population)
kary <- read_tsv(paths$kary)
flow_density <- read_tsv(paths$flow_density)
flow_peak <- read_tsv(paths$flow_peak)
burden <- read_tsv(paths$burden)
invivo_ploidy <- read_tsv(paths$invivo_ploidy)
harvest <- read_tsv(paths$harvest)

ploidy_cols <- c("2N" = "#009E73", "4N" = "#D55E00")
# O2 concentration is encoded from light blue (low) to dark blue (high).
# Anchors follow the concentrations actually used in the in-vitro experiment.
# Most of the perceptual range is assigned to 0--2%, whereas 2--20.5% occupies
# the shorter dark-blue end of the scale.
oxygen_colour_breaks <- c(0, 0.1, 0.2, 0.3, 0.5, 1, 2, 20.5)
oxygen_cols <- c(
  "#8FC5E3", "#80BBDD", "#72B1D7", "#64A7D1",
  "#5599C8", "#3E86BE", "#2171B5", "#084594"
)
dark <- "#333333"
light_grid <- "#E5E5E5"
common_time_breaks <- seq(0, 150, by = 30)
common_time_limits <- c(0, 153)
condition_linetypes <- c(control = "11", `oxygen-deprived` = "solid")

# Figure 1 chart contract
# Question: how do the culture and tumor experiments align in design, sampling,
# and chromosome-number readouts?
# Takeaway: the systems share starting-ploidy cohorts and comparable observation
# classes, but use distinct experimental clocks and outcome units.
# Encoding: lineage = green/orange; O2 = blue gradient; treatment = line type;
# population size / tumor volume / flow sampling / flow peak / direct chromosome
# count / scRNA-seq-derived chromosome profile = filled circle / open circle /
# open triangle / filled triangle / diamond / double circle.

condition_from_segment <- function(x) {
  residue <- gsub("_", "", gsub("20.5", "", as.character(x), fixed = TRUE), fixed = TRUE)
  ifelse(residue == "", "control", "oxygen-deprived")
}

theme_manuscript <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = dark),
      axis.text = element_text(size = 7, colour = dark),
      axis.title = element_text(size = 8),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 7.8, face = "bold", hjust = 0),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 7, 5, 7)
    )
}

save_both <- function(plot, stem, width, height) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    paste0(stem, ".png"), plot = plot, width = width, height = height,
    units = "in", dpi = 300, bg = "white"
  )
  ggsave(
    paste0(stem, ".pdf"), plot = plot, width = width, height = height,
    units = "in", device = grDevices::cairo_pdf, bg = "white"
  )
}

# -----------------------------------------------------------------------------
# A: matched experimental design
# -----------------------------------------------------------------------------

representative_keys <- timeline |>
  group_by(cohort, lineage_label, lineage_terminal_key) |>
  summarise(max_day = max(cumulative_end_day), .groups = "drop") |>
  group_by(cohort, lineage_label) |>
  slice_max(max_day, n = 1, with_ties = FALSE) |>
  ungroup()

invitro_track_positions <- c(
  "4N oxygen-deprived" = 4.00,
  "4N control" = 5.15,
  "2N oxygen-deprived" = 6.80,
  "2N control" = 7.95
)

design_ivt <- timeline |>
  inner_join(
    representative_keys |>
      select(cohort, lineage_label, lineage_terminal_key),
    by = c("cohort", "lineage_label", "lineage_terminal_key")
  ) |>
  mutate(
    condition = recode(lineage_label, deprived = "oxygen-deprived"),
    track = paste(cohort, condition),
    track = factor(
      track,
      levels = c(
        "4N oxygen-deprived", "4N control",
        "2N oxygen-deprived", "2N control"
      )
    ),
    track_y = unname(invitro_track_positions[as.character(track)])
  )

kary_events <- kary |>
  mutate(condition = condition_from_segment(segment_id)) |>
  distinct(cohort, condition, cumulative_end_day) |>
  mutate(
    track = factor(
      paste(cohort, condition),
      levels = levels(design_ivt$track)
    ),
    track_y = unname(invitro_track_positions[as.character(track)])
  )

flow_events <- passage |>
  filter(n_flow > 0) |>
  mutate(condition = condition_from_segment(segment_id)) |>
  distinct(cohort, condition, cumulative_end_day) |>
  mutate(
    track = factor(
      paste(cohort, condition),
      levels = levels(design_ivt$track)
    ),
    track_y = unname(invitro_track_positions[as.character(track)])
  )

# The raw audit confirms a positive initial and linked terminal cell count for
# every passage ID. Mark every reconstructed passage endpoint on the displayed
# representative tracks; shared pre-branch endpoints are repeated on the two
# aligned eventual branches rather than treated as independent measurements.
population_events <- design_ivt |>
  distinct(cohort, condition, track, track_y, cumulative_end_day)

design_ivt_outlines <- design_ivt |>
  group_by(cohort, condition, track, track_y) |>
  summarise(
    cumulative_start_day = min(cumulative_start_day),
    cumulative_end_day = max(cumulative_end_day),
    .groups = "drop"
  )

active_harvests <- harvest |>
  filter(included_under_config, has_ploidy) |>
  distinct(cohort, harvest_day) |>
  group_by(cohort) |>
  filter(cohort != "4N" | harvest_day == max(harvest_day)) |>
  ungroup()

burden_days <- burden |>
  filter(included_under_config) |>
  distinct(cohort, day)

invivo_tracks <- active_harvests |>
  group_by(cohort) |>
  summarise(max_day = max(harvest_day), .groups = "drop") |>
  mutate(
    base_y = ifelse(cohort == "2N", 1.55, 0.60),
    label_y = base_y + 0.30
  )

schematic_tumor_height <- function(day, max_day, cohort) {
  progress <- pmin(pmax(day / max_day, 0), 1)
  shape_exponent <- ifelse(cohort == "2N", 2.00, 3.00)
  0.06 + 0.54 * progress^shape_exponent
}

invivo_funnel <- invivo_tracks |>
  tidyr::crossing(progress = seq(0, 1, length.out = 201)) |>
  mutate(
    day = progress * max_day,
    upper_y = base_y + schematic_tumor_height(day, max_day, cohort)
  )

burden_days <- burden_days |>
  left_join(
    invivo_tracks |> select(cohort, max_day, base_y),
    by = "cohort"
  ) |>
  mutate(
    upper_y = base_y + schematic_tumor_height(day, max_day, cohort)
  )

active_harvests <- active_harvests |>
  left_join(
    invivo_tracks |> select(cohort, max_day, base_y),
    by = "cohort"
  ) |>
  mutate(
    upper_y = base_y + schematic_tumor_height(
      harvest_day, max_day, cohort
    )
  )

p_a <- ggplot() +
  geom_vline(
    xintercept = common_time_breaks,
    colour = light_grid, linewidth = 0.25
  ) +
  geom_hline(
    yintercept = 3.25,
    colour = "#CFCFCF", linewidth = 0.42
  ) +
  geom_rect(
    data = design_ivt,
    aes(
      xmin = cumulative_start_day, xmax = cumulative_end_day,
      ymin = track_y - 0.17, ymax = track_y + 0.17,
      fill = oxygen_pct
    ),
    colour = "white", linewidth = 0.22
  ) +
  geom_rect(
    data = design_ivt_outlines,
    aes(
      xmin = cumulative_start_day, xmax = cumulative_end_day,
      ymin = track_y - 0.19, ymax = track_y + 0.19,
      colour = cohort, linetype = condition
    ),
    fill = NA, linewidth = 0.62
  ) +
  geom_point(
    data = population_events,
    aes(
      x = cumulative_end_day, y = track_y + 0.32,
      colour = cohort
    ),
    shape = 16, size = 0.82, alpha = 0.88
  ) +
  geom_point(
    data = flow_events,
    aes(x = cumulative_end_day, y = track_y + 0.56, colour = cohort),
    shape = 24, size = 1.65, stroke = 0.46,
    fill = "white"
  ) +
  geom_point(
    data = kary_events,
    aes(x = cumulative_end_day, y = track_y + 0.78, colour = cohort),
    shape = 23, size = 1.58, stroke = 0.48,
    fill = "white"
  ) +
  geom_ribbon(
    data = filter(invivo_funnel, cohort == "2N"),
    aes(x = day, ymin = base_y, ymax = upper_y),
    fill = scales::alpha(unname(ploidy_cols[["2N"]]), 0.18),
    colour = NA
  ) +
  geom_ribbon(
    data = filter(invivo_funnel, cohort == "4N"),
    aes(x = day, ymin = base_y, ymax = upper_y),
    fill = scales::alpha(unname(ploidy_cols[["4N"]]), 0.18),
    colour = NA
  ) +
  geom_segment(
    data = invivo_tracks,
    aes(
      x = 0, xend = max_day,
      y = base_y, yend = base_y, colour = cohort
    ),
    linewidth = 0.30, alpha = 0.48, lineend = "round"
  ) +
  geom_point(
    data = burden_days,
    aes(x = day, y = upper_y, colour = cohort),
    shape = 21, fill = "white", size = 1.18, stroke = 0.48,
    alpha = 0.96
  ) +
  geom_point(
    data = active_harvests,
    aes(x = harvest_day, y = upper_y + 0.16, colour = cohort),
    shape = 21, fill = "white", size = 3.0, stroke = 0.65
  ) +
  geom_point(
    data = active_harvests,
    aes(x = harvest_day, y = upper_y + 0.16, colour = cohort),
    shape = 21, fill = "white", size = 1.45, stroke = 0.52
  ) +
  annotate(
    "text", x = 0, y = 9.02, label = "In vitro evolution",
    hjust = 0, size = 2.70, fontface = "bold", colour = dark
  ) +
  annotate(
    "text", x = 0, y = 2.83, label = "In vivo tumors",
    hjust = 0, size = 2.70, fontface = "bold", colour = dark
  ) +
  annotate(
    "text", x = 0, y = 2.48, label = "Inoculation",
    hjust = 0, size = 2.35, colour = dark
  ) +
  scale_colour_manual(values = ploidy_cols, guide = "none") +
  scale_linetype_manual(values = condition_linetypes, guide = "none") +
  scale_fill_gradientn(
    colours = oxygen_cols,
    values = rescale(oxygen_colour_breaks),
    limits = c(0, 20.5),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = c(0.90, 1.85, 4.00, 5.15, 6.80, 7.95),
    labels = c(
      "4N tumors", "2N tumors",
      "4N oxygen-deprived", "4N control",
      "2N oxygen-deprived", "2N control"
    ),
    limits = c(0.45, 9.20),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    breaks = common_time_breaks, limits = common_time_limits,
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "A. Matched experimental design",
    subtitle = "Shared elapsed-day scale; each experimental clock starts from its own origin",
    x = "Elapsed day from each experiment's start", y = NULL
  ) +
  theme_manuscript() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 7.8, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(
      size = 5.7, colour = "#666666", margin = margin(b = 3)
    ),
    plot.margin = margin(2, 7, 4, 13)
  )

oxygen_key_left <- 0.855
oxygen_key_right <- 0.975
oxygen_key_low_fraction <- 0.78
oxygen_key_position <- function(oxygen_pct) {
  ifelse(
    oxygen_pct <= 2,
    oxygen_key_low_fraction * oxygen_pct / 2,
    oxygen_key_low_fraction +
      (1 - oxygen_key_low_fraction) * (oxygen_pct - 2) / (20.5 - 2)
  )
}
oxygen_key_data <- data.frame(
  key_position = seq(0, 1, length.out = 256)
) |>
  mutate(
    x = oxygen_key_left +
      (oxygen_key_right - oxygen_key_left) * key_position,
    y = 0.72,
    oxygen_pct = ifelse(
      key_position <= oxygen_key_low_fraction,
      2 * key_position / oxygen_key_low_fraction,
      2 + (20.5 - 2) *
        (key_position - oxygen_key_low_fraction) /
        (1 - oxygen_key_low_fraction)
    )
  )
oxygen_key_ticks <- data.frame(
  oxygen_pct = c(0, 2, 20.5),
  label = c("0", "2", "20.5")
) |>
  mutate(
    x = oxygen_key_left +
      (oxygen_key_right - oxygen_key_left) *
        oxygen_key_position(oxygen_pct)
  )

# -----------------------------------------------------------------------------
# B: observed in-vitro ploidy trajectories and distributions
# -----------------------------------------------------------------------------

kary_plot <- kary |>
  mutate(
    condition = factor(
      condition_from_segment(segment_id),
      levels = c("control", "oxygen-deprived")
    ),
    cohort = factor(cohort, levels = c("2N", "4N")),
    observed_ploidy = observed_kary_N / 22
  )

start_rows <- kary_plot |>
  group_by(cohort) |>
  filter(cumulative_end_day == min(cumulative_end_day)) |>
  ungroup()

kary_summary <- kary_plot |>
  group_by(cohort, condition, cumulative_end_day) |>
  summarise(
    median_ploidy = median(observed_ploidy),
    q1 = quantile(observed_ploidy, 0.25),
    q3 = quantile(observed_ploidy, 0.75),
    .groups = "drop"
  )

start_summary <- start_rows |>
  group_by(cohort, cumulative_end_day) |>
  summarise(
    median_ploidy = median(observed_ploidy),
    q1 = quantile(observed_ploidy, 0.25),
    q3 = quantile(observed_ploidy, 0.75),
    .groups = "drop"
  )

# The first chromosome-count sample is the shared starting reference for both
# branches. Repeat those directly counted cells in each condition row so that
# the two separated trajectories begin from the same observed distribution.
kary_plot_display <- bind_rows(
  kary_plot,
  start_rows |>
    mutate(
      condition = factor(
        "oxygen-deprived",
        levels = c("control", "oxygen-deprived")
      )
    )
) |>
  distinct(
    cohort, condition, cumulative_end_day,
    cell_index, observed_kary_N, observed_ploidy,
    .keep_all = TRUE
  )

# Match panels C and D by displaying each directly counted cell distribution as
# a left half-violin, a central boxplot, and individual observations to the
# right. The horizontal offsets are display-only and do not change sample day.
kary_violin_width_days <- 6.0
kary_half_violin <- kary_plot_display |>
  group_by(cohort, condition, cumulative_end_day) |>
  group_modify(function(.x, .y) {
    values <- .x$observed_ploidy[is.finite(.x$observed_ploidy)]
    density_curve <- stats::density(
      values,
      n = 256,
      adjust = 0.85,
      from = max(1, min(values) - 0.10),
      to = min(5, max(values) + 0.10)
    )
    scaled_density <- density_curve$y / max(density_curve$y)
    day <- .y$cumulative_end_day[[1L]]
    tibble(
      display_x = c(
        day,
        day - kary_violin_width_days * scaled_density,
        day
      ),
      observed_ploidy = c(
        density_curve$x[[1L]],
        density_curve$x,
        density_curve$x[[length(density_curve$x)]]
      )
    )
  }) |>
  ungroup()

kary_point_display <- kary_plot_display |>
  mutate(point_display_x = cumulative_end_day + 1.8)

trajectory_summary <- bind_rows(
  kary_summary,
  start_summary |> mutate(condition = "oxygen-deprived")
) |>
  distinct(cohort, condition, cumulative_end_day, .keep_all = TRUE) |>
  mutate(condition = factor(condition, levels = c("control", "oxygen-deprived")))

flow_density_plot <- flow_density |>
  mutate(
    cohort = factor(cohort, levels = c("2N", "4N")),
    condition = factor(condition, levels = c("control", "oxygen-deprived"))
  ) |>
  group_by(cohort, condition, cumulative_end_day, ploidy) |>
  summarise(
    pooled_density_mass = mean(normalized_density_mass),
    .groups = "drop"
  )

flow_violin_width_days <- 7.5
flow_half_violin <- flow_density_plot |>
  group_by(cohort, condition, cumulative_end_day) |>
  arrange(ploidy, .by_group = TRUE) |>
  group_modify(function(.x, .y) {
    cumulative_mass <- cumsum(.x$pooled_density_mass) /
      sum(.x$pooled_density_mass)
    keep <- cumulative_mass >= 0.001 & cumulative_mass <= 0.999
    keep[[which.max(.x$pooled_density_mass)]] <- TRUE
    keep <- keep & .x$ploidy >= 1 & .x$ploidy <= 5
    .x <- .x[keep, , drop = FALSE]
    scaled <- .x$pooled_density_mass / max(.x$pooled_density_mass)
    day <- .y$cumulative_end_day[[1L]]
    tibble(
      display_x = c(day, day - flow_violin_width_days * scaled, day),
      ploidy = c(.x$ploidy[[1L]], .x$ploidy, .x$ploidy[[nrow(.x)]])
    )
  }) |>
  ungroup()

weighted_quantile <- function(x, weight, probabilities) {
  keep <- is.finite(x) & is.finite(weight) & weight > 0
  x <- x[keep]
  weight <- weight[keep]
  order_index <- order(x)
  x <- x[order_index]
  weight <- weight[order_index]
  cumulative <- cumsum(weight) / sum(weight)
  vapply(
    probabilities,
    function(probability) x[[which(cumulative >= probability)[[1L]]]],
    numeric(1L)
  )
}

flow_box_stats <- flow_density_plot |>
  group_by(cohort, condition, cumulative_end_day) |>
  group_modify(function(.x, .y) {
    quantiles <- weighted_quantile(
      .x$ploidy, .x$pooled_density_mass,
      c(0.25, 0.50, 0.75)
    )
    iqr <- quantiles[[3L]] - quantiles[[1L]]
    support <- .x$ploidy[
      .x$pooled_density_mass >= max(.x$pooled_density_mass) * 1e-6 &
        .x$ploidy >= quantiles[[1L]] - 1.5 * iqr &
        .x$ploidy <= quantiles[[3L]] + 1.5 * iqr
    ]
    tibble(
      lower_whisker = min(support),
      q1 = quantiles[[1L]],
      median = quantiles[[2L]],
      q3 = quantiles[[3L]],
      upper_whisker = max(support)
    )
  }) |>
  ungroup()

flow_peak_plot <- flow_peak |>
  mutate(
    cohort = factor(cohort, levels = c("2N", "4N")),
    condition = factor(condition, levels = c("control", "oxygen-deprived"))
  ) |>
  group_by(cohort, condition, cumulative_end_day) |>
  arrange(sample_name, .by_group = TRUE) |>
  mutate(
    peak_display_x = cumulative_end_day +
      ifelse(n() == 1L, 1.4, seq(0.7, 1.9, length.out = n()))
  ) |>
  ungroup()

measurement_facet <- facet_grid(
  condition ~ cohort,
  switch = "y",
  labeller = labeller(
    condition = c(
      control = "Control",
      `oxygen-deprived` = "O\u2082-depr."
    )
  )
)

measurement_panel_theme <- theme(
  legend.position = "none",
  axis.text = element_text(size = 5.0, colour = dark),
  axis.title = element_text(size = 5.5),
  plot.title = element_text(size = 7.8, face = "bold", hjust = 0),
  plot.tag = element_text(size = 8.0, face = "bold"),
  panel.spacing.x = grid::unit(2.0, "mm"),
  panel.spacing.y = grid::unit(1.5, "mm"),
  panel.border = element_rect(
    fill = NA, colour = "#AEB6C0", linewidth = 0.32
  ),
  strip.background = element_rect(
    fill = "#F4F6F8", colour = "#AEB6C0", linewidth = 0.32
  ),
  strip.text.x = element_text(
    size = 4.9, face = "bold",
    margin = margin(t = 0.15, b = 0.15)
  ),
  strip.placement = "outside",
  strip.text.y.left = element_text(
    angle = 90, size = 4.5, face = "bold",
    margin = margin(l = 0.4, r = 0.4)
  ),
  plot.margin = margin(1.5, 3, 1.5, 3)
)

p_b <- ggplot() +
  geom_hline(
    data = data.frame(
      cohort = factor(c("2N", "4N"), levels = c("2N", "4N")),
      reference = c(2, 4)
    ),
    aes(yintercept = reference, colour = cohort),
    linewidth = 0.38, linetype = 3, alpha = 0.65
  ) +
  geom_polygon(
    data = kary_half_violin,
    aes(
      x = display_x, y = observed_ploidy,
      fill = cohort,
      group = interaction(cohort, condition, cumulative_end_day)
    ),
    alpha = 0.60, colour = dark, linewidth = 0.10
  ) +
  geom_boxplot(
    data = kary_plot_display,
    aes(
      x = cumulative_end_day, y = observed_ploidy,
      group = interaction(cohort, condition, cumulative_end_day)
    ),
    width = 2.20, outlier.shape = NA,
    fill = "white", colour = dark, alpha = 0.86,
    linewidth = 0.10
  ) +
  geom_point(
    data = kary_point_display,
    aes(
      x = point_display_x, y = observed_ploidy,
      colour = cohort
    ),
    position = position_jitter(width = 1.0, height = 0, seed = 23),
    shape = 23, fill = "white", size = 0.41, stroke = 0.125,
    alpha = 0.62
  ) +
  measurement_facet +
  scale_colour_manual(values = ploidy_cols, guide = "none") +
  scale_fill_manual(values = ploidy_cols, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 75, 150),
    limits = c(0, 153), expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    breaks = 1:5,
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(ylim = c(1, 5), expand = FALSE) +
  labs(
    title = "B. Direct chromosome-count distributions",
    x = NULL, y = "Ploidy"
  ) +
  theme_manuscript() +
  measurement_panel_theme

# -----------------------------------------------------------------------------
# C: observed flow-cytometry ploidy distributions
# -----------------------------------------------------------------------------

p_c <- ggplot() +
  geom_hline(
    data = data.frame(
      cohort = factor(c("2N", "4N"), levels = c("2N", "4N")),
      reference = c(2, 4)
    ),
    aes(yintercept = reference, colour = cohort),
    linewidth = 0.38, linetype = 3, alpha = 0.65
  ) +
  geom_polygon(
    data = flow_half_violin,
    aes(
      x = display_x, y = ploidy,
      group = interaction(cohort, condition, cumulative_end_day),
      fill = cohort, colour = cohort
    ),
    alpha = 0.34, linewidth = 0.12
  ) +
  geom_segment(
    data = flow_box_stats,
    aes(
      x = cumulative_end_day, xend = cumulative_end_day,
      y = lower_whisker, yend = upper_whisker
    ),
    colour = dark, linewidth = 0.20
  ) +
  geom_rect(
    data = flow_box_stats,
    aes(
      xmin = cumulative_end_day - 0.65,
      xmax = cumulative_end_day + 0.65,
      ymin = q1, ymax = q3
    ),
    fill = "white", colour = dark, linewidth = 0.18
  ) +
  geom_segment(
    data = flow_box_stats,
    aes(
      x = cumulative_end_day - 0.65,
      xend = cumulative_end_day + 0.65,
      y = median, yend = median
    ),
    colour = dark, linewidth = 0.24
  ) +
  geom_point(
    data = flow_peak_plot,
    aes(
      x = peak_display_x, y = peak_ploidy,
      colour = cohort
    ),
    shape = 24, fill = "white", size = 0.825, stroke = 0.23
  ) +
  measurement_facet +
  scale_colour_manual(values = ploidy_cols, guide = "none") +
  scale_fill_manual(values = ploidy_cols, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 75, 150),
    limits = c(0, 153), expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    breaks = 1:5,
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(ylim = c(1, 5), expand = FALSE) +
  labs(
    title = "C. Flow-cytometry ploidy distributions",
    x = "Reconstructed culture day", y = "Ploidy"
  ) +
  theme_manuscript() +
  measurement_panel_theme

# -----------------------------------------------------------------------------
# D: starting cell-line reference versus terminal tumors
# -----------------------------------------------------------------------------

baseline <- start_rows |>
  transmute(
    cohort,
    stage = "Starting\nreference",
    ploidy = observed_kary_N / 22,
    harvest = NA_character_
  ) |>
  mutate(stage = factor(stage, levels = c("Starting\nreference", "Terminal\ntumors")))

terminal <- invivo_ploidy |>
  filter(used_in_endpoint_loss, is.finite(endpoint_value)) |>
  transmute(
    cohort = factor(cohort, levels = c("2N", "4N")),
    stage = "Terminal\ntumors",
    ploidy = endpoint_value / 22,
    harvest
  ) |>
  mutate(stage = factor(stage, levels = c("Starting\nreference", "Terminal\ntumors")))

all_distributions <- bind_rows(baseline, terminal) |>
  mutate(
    cohort = factor(cohort, levels = c("2N", "4N")),
    stage = factor(
      stage,
      levels = c("Starting\nreference", "Terminal\ntumors")
    ),
    stage_x = as.numeric(stage)
  )

harvest_medians <- terminal |>
  group_by(cohort, harvest) |>
  summarise(
    stage = "Terminal\ntumors",
    ploidy = median(ploidy),
    .groups = "drop"
  ) |>
  mutate(
    stage = factor(
      stage,
      levels = c("Starting\nreference", "Terminal\ntumors")
    ),
    stage_x = as.numeric(stage)
  )

make_left_half_violin <- function(data, width = 0.29, adjust = 0.85) {
  data |>
    group_by(cohort, stage) |>
    group_modify(function(.x, .y) {
      values <- .x$ploidy[
        is.finite(.x$ploidy)
      ]
      stage_x <- match(
        as.character(.y$stage[[1]]),
        levels(data$stage)
      )
      density_curve <- stats::density(
        values,
        n = 256,
        adjust = adjust,
        from = max(1, min(values) - 0.10),
        to = min(5, max(values) + 0.10)
      )
      scaled_density <- density_curve$y / max(density_curve$y)
      tibble(
        stage_x = c(
          stage_x,
          stage_x - width * scaled_density,
          stage_x
        ),
        ploidy = c(
          density_curve$x[[1]],
          density_curve$x,
          density_curve$x[[length(density_curve$x)]]
        )
      )
    }) |>
    ungroup()
}

half_violin <- make_left_half_violin(all_distributions)

p_d <- ggplot() +
  geom_hline(
    data = data.frame(
      cohort = factor(c("2N", "4N"), levels = c("2N", "4N")),
      reference = c(2, 4)
    ),
    aes(yintercept = reference, colour = cohort),
    linewidth = 0.50, linetype = 3, alpha = 0.65
  ) +
  geom_polygon(
    data = half_violin,
    aes(
      x = stage_x, y = ploidy,
      fill = cohort, group = interaction(cohort, stage)
    ),
    alpha = 0.60, colour = dark, linewidth = 0.10
  ) +
  geom_boxplot(
    data = all_distributions,
    aes(
      x = stage_x, y = ploidy,
      group = interaction(cohort, stage)
    ),
    width = 0.105, outlier.shape = NA,
    fill = "white", colour = dark, alpha = 0.86,
    linewidth = 0.10
  ) +
  geom_point(
    data = baseline,
    aes(
      x = as.numeric(stage) + 0.18,
      y = ploidy, colour = cohort
    ),
    position = position_jitter(width = 0.038, height = 0, seed = 17),
    shape = 23, fill = "white", size = 0.46, stroke = 0.19,
    alpha = 0.55
  ) +
  geom_point(
    data = harvest_medians,
    aes(
      x = stage_x + 0.18,
      y = ploidy, colour = cohort
    ),
    position = position_jitter(width = 0.038, height = 0, seed = 17),
    shape = 21, fill = "white", size = 0.95, stroke = 0.29
  ) +
  geom_point(
    data = harvest_medians,
    aes(
      x = stage_x + 0.18,
      y = ploidy, colour = cohort
    ),
    position = position_jitter(width = 0.038, height = 0, seed = 17),
    shape = 21, fill = "white", size = 0.51, stroke = 0.24
  ) +
  facet_wrap(~cohort, nrow = 1) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Starting\nreference", "Terminal\ntumors"),
    limits = c(0.55, 2.45),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_colour_manual(values = ploidy_cols, guide = "none") +
  scale_fill_manual(values = ploidy_cols, guide = "none") +
  scale_y_continuous(
    breaks = 1:5,
    limits = c(1, 5), expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "D. Starting state and terminal tumors",
    x = NULL, y = "Ploidy"
  ) +
  theme_manuscript() +
  theme(
    axis.text = element_text(size = 5.0, colour = dark),
    axis.text.x = element_text(size = 5.0, lineheight = 0.90),
    axis.title = element_text(size = 5.5),
    strip.text.x = element_text(size = 4.9, face = "bold"),
    plot.title = element_text(size = 7.8, face = "bold", hjust = 0),
    plot.tag = element_text(size = 8.0, face = "bold"),
    panel.spacing = grid::unit(2.5, "mm"),
    plot.margin = margin(1.5, 3, 1.5, 3)
  )

# Shared encoding key for the assembled figure. Keeping the semantic channels
# in one place avoids asking readers to reinterpret the same color or symbol
# across the experimental-design and outcome panels.
p_encoding_key <- ggplot() +
  geom_tile(
    data = oxygen_key_data,
    aes(x = x, y = y, fill = oxygen_pct),
    width = (oxygen_key_right - oxygen_key_left) / 255 * 1.04,
    height = 0.20
  ) +
  geom_segment(
    data = oxygen_key_ticks,
    aes(x = x, xend = x, y = 0.62, yend = 0.82),
    inherit.aes = FALSE, linewidth = 0.28, colour = "#555555"
  ) +
  geom_text(
    data = oxygen_key_ticks,
    aes(x = x, y = 0.48, label = label),
    inherit.aes = FALSE, size = 1.95, colour = dark
  ) +
  annotate(
    "text", x = 0.015, y = 0.72, label = "Starting lineage",
    hjust = 0, vjust = 0.5, size = 2.15, fontface = "bold", colour = dark
  ) +
  annotate(
    "segment", x = 0.142, xend = 0.172, y = 0.72, yend = 0.72,
    linewidth = 0.85, colour = unname(ploidy_cols[["2N"]])
  ) +
  annotate(
    "text", x = 0.178, y = 0.72, label = "2N",
    hjust = 0, vjust = 0.5, size = 2.15, colour = dark
  ) +
  annotate(
    "segment", x = 0.220, xend = 0.250, y = 0.72, yend = 0.72,
    linewidth = 0.85, colour = unname(ploidy_cols[["4N"]])
  ) +
  annotate(
    "text", x = 0.256, y = 0.72, label = "4N",
    hjust = 0, vjust = 0.5, size = 2.15, colour = dark
  ) +
  annotate(
    "text", x = 0.330, y = 0.72, label = "Condition",
    hjust = 0, vjust = 0.5, size = 2.15, fontface = "bold", colour = dark
  ) +
  annotate(
    "segment", x = 0.412, xend = 0.452, y = 0.72, yend = 0.72,
    linewidth = 0.55, linetype = "11", colour = dark
  ) +
  annotate(
    "text", x = 0.460, y = 0.72, label = "control",
    hjust = 0, vjust = 0.5, size = 2.15, colour = dark
  ) +
  annotate(
    "segment", x = 0.550, xend = 0.590, y = 0.72, yend = 0.72,
    linewidth = 0.55, linetype = "solid", colour = dark
  ) +
  annotate(
    "text", x = 0.598, y = 0.72, label = "oxygen-deprived",
    hjust = 0, vjust = 0.5, size = 2.15, colour = dark
  ) +
  annotate(
    "text", x = 0.755, y = 0.72, label = "Target O\u2082 (%)",
    hjust = 0, vjust = 0.5, size = 2.10, fontface = "bold", colour = dark
  ) +
  annotate(
    "text", x = 0.015, y = 0.34, label = "Measurements",
    hjust = 0, vjust = 0.5, size = 2.15, fontface = "bold", colour = dark
  ) +
  annotate(
    "point", x = 0.126, y = 0.34,
    shape = 16, size = 0.75, colour = dark
  ) +
  annotate(
    "text", x = 0.143, y = 0.34, label = "population size",
    hjust = 0, vjust = 0.5, size = 1.90, colour = dark
  ) +
  annotate(
    "point", x = 0.270, y = 0.34,
    shape = 21, fill = "white", size = 0.85, stroke = 0.24,
    colour = dark
  ) +
  annotate(
    "text", x = 0.287, y = 0.34, label = "volume measurement",
    hjust = 0, vjust = 0.5, size = 1.90, colour = dark
  ) +
  annotate(
    "point", x = 0.465, y = 0.34,
    shape = 23, size = 0.95, stroke = 0.24,
    fill = "white", colour = dark
  ) +
  annotate(
    "text", x = 0.482, y = 0.34, label = "direct chromosome count",
    hjust = 0, vjust = 0.5, size = 1.88, colour = dark
  ) +
  annotate(
    "point", x = 0.126, y = 0.08,
    shape = 24, size = 0.975, stroke = 0.24,
    fill = "white", colour = dark
  ) +
  annotate(
    "text", x = 0.143, y = 0.08, label = "flow cytometry sampling",
    hjust = 0, vjust = 0.5, size = 1.90, colour = dark
  ) +
  annotate(
    "point", x = 0.365, y = 0.08,
    shape = 21, size = 1.05, stroke = 0.26,
    fill = "white", colour = dark
  ) +
  annotate(
    "point", x = 0.365, y = 0.08,
    shape = 21, size = 0.55, stroke = 0.225,
    fill = "white", colour = dark
  ) +
  annotate(
    "text", x = 0.383, y = 0.08, label = "scRNA-seq profile",
    hjust = 0, vjust = 0.5, size = 1.90, colour = dark
  ) +
  scale_fill_gradientn(
    colours = oxygen_cols,
    values = rescale(oxygen_colour_breaks),
    limits = c(0, 20.5),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(1, 8, 1, 8))

# -----------------------------------------------------------------------------
# Save subpanels and assembled draft
# -----------------------------------------------------------------------------

save_both(
  p_a,
  file.path(subpanel_root, "fig1a_matched_design"),
  width = 4.37, height = 5.325
)
save_both(
  p_b,
  file.path(subpanel_root, "fig1b_direct_chromosome_ploidy"),
  width = 2.73, height = 1.72
)
save_both(
  p_c,
  file.path(subpanel_root, "fig1c_flow_cytometry_ploidy"),
  width = 2.73, height = 1.72
)
save_both(
  p_d,
  file.path(subpanel_root, "fig1d_invivo_start_terminal_distributions"),
  width = 2.73, height = 1.20
)

right_column <- p_b / plot_spacer() / p_c / plot_spacer() / p_d +
  plot_layout(heights = c(1, 0.005, 1, 0.005, 0.62))

figure1_body <- p_a | right_column
figure1_body <- figure1_body +
  plot_layout(widths = c(1.6, 1))

figure1 <- figure1_body / p_encoding_key +
  plot_layout(heights = c(1, 0.145))

save_both(
  figure1,
  file.path(figures_root, "assembled_fig1"),
  width = 7.1, height = 5.325
)

dir.create(revised_figures_root, recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(
  file.path(figures_root, "assembled_fig1.png"),
  file.path(revised_figures_root, "assembled_fig1.png"),
  overwrite = TRUE, copy.date = TRUE
))
invisible(file.copy(
  file.path(figures_root, "assembled_fig1.pdf"),
  file.path(revised_figures_root, "assembled_fig1.pdf"),
  overwrite = TRUE, copy.date = TRUE
))

summary_table <- bind_rows(
  trajectory_summary |>
    transmute(
      panel = "F1B", cohort = as.character(cohort),
      group = as.character(condition), time_day = cumulative_end_day,
      estimate = median_ploidy, lower = q1, upper = q3
    ),
  flow_peak_plot |>
    transmute(
      panel = "F1C", cohort = as.character(cohort),
      group = paste(as.character(condition), "flow cytometry modal ploidy"),
      time_day = cumulative_end_day,
      estimate = peak_ploidy, lower = NA_real_, upper = NA_real_
    ),
  harvest_medians |>
    transmute(
      panel = "F1D", cohort = as.character(cohort),
      group = harvest, time_day = NA_real_,
      estimate = ploidy, lower = NA_real_, upper = NA_real_
    )
)
write.table(
  summary_table,
  file.path(data_root, "figure1_displayed_summaries.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

read_png_dimensions <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  header <- readBin(con, what = "raw", n = 24L)
  stopifnot(length(header) == 24L)
  c(
    width = readBin(header[17:20], what = "integer", n = 1L, size = 4L, endian = "big"),
    height = readBin(header[21:24], what = "integer", n = 1L, size = 4L, endian = "big")
  )
}

figure1_png <- file.path(figures_root, "assembled_fig1.png")
figure1_dimensions <- read_png_dimensions(figure1_png)
output_metadata <- data.frame(
  output = figure1_png,
  width = unname(figure1_dimensions[["width"]]),
  height = unname(figure1_dimensions[["height"]]),
  dpi = 300,
  plotting_contract = "iteration2_o2_low_range_expanded_large_inline_panel_titles",
  stringsAsFactors = FALSE
)
write.table(
  output_metadata,
  file.path(data_root, "figure1_output_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

upstream_input_root <- file.path(
  repo_root,
  "data",
  "InVitroData_LTEE"
)
data_contract <- read_tsv(file.path(data_root, "data_contract.tsv"))
population_upstream <- data_contract$source[
  data_contract$role == "raw population-size audit source"
][[1L]]
flow_upstream <- data_contract$source[
  data_contract$role == "observed flow-cytometry density source"
][[1L]]
provenance <- data.frame(
  input_role = names(paths),
  local_path = unname(unlist(paths)),
  upstream_path = c(
    file.path(upstream_input_root, basename(paths$timeline)),
    file.path(upstream_input_root, basename(paths$passage)),
    population_upstream,
    file.path(upstream_input_root, basename(paths$kary)),
    flow_upstream,
    flow_upstream,
    file.path(upstream_input_root, basename(paths$burden)),
    file.path(upstream_input_root, basename(paths$invivo_ploidy)),
    file.path(upstream_input_root, basename(paths$harvest))
  ),
  stringsAsFactors = FALSE
)
write.table(
  provenance,
  file.path(data_root, "figure1_source_file_provenance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

measurement_validation <- data.frame(
  metric = c(
    "raw_invitro_passage_records",
    "raw_passages_with_complete_population_measurement",
    "displayed_invitro_population_markers",
    "observed_flow_samples",
    "displayed_flow_group_timepoints",
    "minimum_flow_density_mass_retained_in_ploidy_1_to_5",
    "displayed_open_triangle_flow_markers",
    "figure1B_flow_layers",
    "figure1B_half_violin_groups",
    "figure1B_direct_cell_points",
    "figure1D_sample_count_labels",
    "displayed_4N_scRNA_seq_markers",
    "chromosome_count_per_ploidy_unit",
    "figure1A_measurement_symbol_size_scale",
    "figure1BCD_measurement_symbol_size_scale",
    "right_panel_spacer_weight"
  ),
  value = c(
    nrow(population),
    sum(population$population_measurement_complete),
    nrow(population_events),
    dplyr::n_distinct(flow_density$sample_name),
    nrow(flow_box_stats),
    min(
      flow_density |>
        group_by(sample_name) |>
        summarise(
          retained = sum(
            normalized_density_mass[ploidy >= 1 & ploidy <= 5]
          ),
          .groups = "drop"
        ) |>
        pull(retained)
    ),
    nrow(flow_peak_plot),
    0,
    dplyr::n_distinct(interaction(
      kary_plot_display$cohort,
      kary_plot_display$condition,
      kary_plot_display$cumulative_end_day,
      drop = TRUE
    )),
    nrow(kary_point_display),
    0,
    sum(active_harvests$cohort == "4N"),
    22,
    1.0,
    0.5,
    0.005
  ),
  stringsAsFactors = FALSE
)
write.table(
  measurement_validation,
  file.path(data_root, "figure1_measurement_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Figure 1 written to: ", figure1_png)

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Figure1 <- function() {
  data_dir <- file.path(DATA_ROOT, "Figure1")
  require_files(
    file.path(data_dir, c(
      "invitro_kary_cells.tsv",
      "invitro_lineage_timeline.tsv",
      "invitro_passage_observations.tsv",
      "invitro_population_measurement_events.tsv",
      "invitro_flow_observed_density.tsv",
      "invitro_flow_peak_ploidy.tsv",
      "invivo_burden_long.tsv",
      "invivo_harvest_catalog.tsv",
      "invivo_ploidy_cells.tsv"
    )),
    "Figure 1 frozen input"
  )
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_Figure1.R"), mustWork = TRUE),
    env = c(
      "FIGURE1_DRAW_WORKER=1",
      paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("FIGURE1_DATA_DIR=", data_dir),
      paste0("FIGURE1_BUILD_DIR=", file.path(data_dir, "plot_build")),
      paste0("FIGURE_OUTPUT_DIR=", OUTPUT_ROOT)
    )
  )
  require_files(
    file.path(OUTPUT_ROOT, c("assembled_fig1.png", "assembled_fig1.pdf")),
    "Figure 1 output"
  )
  invisible(file.path(OUTPUT_ROOT, "assembled_fig1.png"))
}

if (sys.nframe() == 0L) draw_Figure1()

}
