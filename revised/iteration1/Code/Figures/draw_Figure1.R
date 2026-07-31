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
  kary = file.path(data_root, "invitro_kary_cells.tsv"),
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
kary <- read_tsv(paths$kary)
burden <- read_tsv(paths$burden)
invivo_ploidy <- read_tsv(paths$invivo_ploidy)
harvest <- read_tsv(paths$harvest)

ploidy_cols <- c("2N" = "#009E73", "4N" = "#D55E00")
oxygen_cols <- c("#084594", "#2171B5", "#4292C6", "#8FC5E3")
dark <- "#333333"
light_grid <- "#E5E5E5"

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
      plot.title = element_text(size = 6.3, face = "bold", hjust = 0),
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
    track_y = as.numeric(track)
  )

kary_events <- kary |>
  mutate(condition = condition_from_segment(segment_id)) |>
  distinct(cohort, condition, cumulative_end_day) |>
  mutate(
    track = factor(
      paste(cohort, condition),
      levels = levels(design_ivt$track)
    ),
    track_y = as.numeric(track)
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
    track_y = as.numeric(track)
  )

p_design_ivt <- ggplot(design_ivt) +
  geom_rect(
    aes(
      xmin = cumulative_start_day, xmax = cumulative_end_day,
      ymin = track_y - 0.17, ymax = track_y + 0.17,
      fill = oxygen_pct
    ),
    colour = "white", linewidth = 0.22
  ) +
  geom_point(
    data = flow_events,
    aes(x = cumulative_end_day, y = track_y + 0.29),
    shape = 24, size = 1.85, stroke = 0.42,
    fill = "#D9EAF5", colour = "#356B96"
  ) +
  geom_point(
    data = kary_events,
    aes(x = cumulative_end_day, y = track_y + 0.29),
    shape = 21, size = 1.50, stroke = 0.42,
    fill = "#F6C978", colour = "#8A5A13"
  ) +
  scale_fill_gradientn(
    colours = oxygen_cols,
    values = rescale(c(0, 1, 5, 20.5)),
    limits = c(0, 20.5),
    breaks = c(0, 5, 20.5),
    labels = c("0", "5", "20.5"),
    name = expression("Target O"[2]*" (%)"),
    guide = guide_colourbar(
      direction = "horizontal",
      title.position = "left",
      label.position = "bottom",
      barwidth = grid::unit(35, "mm"),
      barheight = grid::unit(1.45, "mm"),
      ticks.colour = dark
    )
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(design_ivt$track)),
    labels = levels(design_ivt$track),
    expand = expansion(add = c(0.35, 0.55))
  ) +
  scale_x_continuous(
    breaks = c(0, 30, 60, 90, 120, 150),
    limits = c(0, 153), expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    tag = "A",
    title = "Matched experimental design",
    subtitle = "In vitro evolution",
    x = "In-vitro experimental day", y = NULL
  ) +
  theme_manuscript() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 6.3, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(
      size = 6.0, face = "bold", margin = margin(b = 3)
    ),
    plot.margin = margin(2, 6, 3, 13)
  )

active_harvests <- harvest |>
  filter(included_under_config, has_ploidy) |>
  distinct(cohort, harvest_day)

burden_days <- burden |>
  filter(included_under_config) |>
  distinct(cohort, day)

invivo_tracks <- active_harvests |>
  group_by(cohort) |>
  summarise(max_day = max(harvest_day), .groups = "drop") |>
  mutate(y = ifelse(cohort == "2N", 2, 1))

burden_days <- burden_days |>
  left_join(invivo_tracks |> select(cohort, y), by = "cohort")

active_harvests <- active_harvests |>
  left_join(invivo_tracks |> select(cohort, y), by = "cohort")

p_design_vivo <- ggplot(invivo_tracks, aes(colour = cohort)) +
  geom_rect(
    aes(
      xmin = 0, xmax = max_day,
      ymin = y - 0.105, ymax = y + 0.105,
      fill = cohort
    ),
    colour = NA, alpha = 0.14
  ) +
  geom_segment(
    aes(x = 0, xend = max_day, y = y, yend = y),
    linewidth = 0.58, lineend = "round"
  ) +
  geom_point(
    data = burden_days,
    aes(x = day, y = y),
    colour = dark, size = 0.92, alpha = 0.72
  ) +
  geom_point(
    data = active_harvests,
    aes(x = harvest_day, y = y),
    shape = 23, fill = "white", size = 2.8, stroke = 0.65
  ) +
  annotate(
    "text", x = 0, y = 2.45, label = "Inoculation",
    hjust = 0, size = 2.50, colour = dark
  ) +
  scale_colour_manual(values = ploidy_cols, guide = "none") +
  scale_fill_manual(values = ploidy_cols, guide = "none") +
  scale_y_continuous(
    breaks = c(1, 2), labels = c("4N tumors", "2N tumors"),
    limits = c(0.05, 2.55)
  ) +
  scale_x_continuous(
    breaks = c(0, 30, 60, 90), limits = c(0, 94),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "In vivo tumors",
    x = "In-vivo experimental day", y = NULL
  ) +
  theme_manuscript() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 6.0, face = "bold", margin = margin(b = 1)),
    plot.margin = margin(2, 7, 3, 13)
  )

oxygen_key_left <- 0.31
oxygen_key_right <- 0.86
oxygen_key_data <- data.frame(
  x = seq(oxygen_key_left, oxygen_key_right, length.out = 256),
  y = 0.31,
  oxygen_pct = seq(0, 20.5, length.out = 256)
)
oxygen_key_ticks <- data.frame(
  oxygen_pct = c(0, 5, 20.5),
  label = c("0", "5", "20.5")
) |>
  mutate(
    x = oxygen_key_left +
      (oxygen_key_right - oxygen_key_left) * oxygen_pct / 20.5
  )

p_o2_key <- ggplot(oxygen_key_data) +
  geom_tile(
    aes(x = x, y = y, fill = oxygen_pct),
    width = (oxygen_key_right - oxygen_key_left) / 255 * 1.04,
    height = 0.32
  ) +
  geom_segment(
    data = oxygen_key_ticks,
    aes(x = x, xend = x, y = 0.15, yend = 0.47),
    inherit.aes = FALSE, linewidth = 0.30, colour = "#555555"
  ) +
  geom_text(
    data = oxygen_key_ticks,
    aes(x = x, y = 0.02, label = label),
    inherit.aes = FALSE, size = 2.25, colour = dark
  ) +
  annotate(
    "text", x = 0.31, y = 0.78,
    label = "\u25B2 flow    \u25CF chromosome count",
    hjust = 0, vjust = 0.5, size = 2.35, colour = dark
  ) +
  annotate(
    "text", x = 0.02, y = 0.31,
    label = "Target O\u2082 (%)",
    hjust = 0, vjust = 0.5, size = 2.45, colour = dark
  ) +
  scale_fill_gradientn(
    colours = oxygen_cols,
    values = rescale(c(0, 1, 5, 20.5)),
    limits = c(0, 20.5),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 6, 2, 13))

p_vivo_key <- ggplot() +
  annotate(
    "text", x = 0.02, y = 0.72,
    label = "\u25CF burden measurements",
    hjust = 0, vjust = 0.5, size = 2.00, colour = dark
  ) +
  annotate(
    "text", x = 0.52, y = 0.72,
    label = "\u25C7 terminal chromosome count",
    hjust = 0, vjust = 0.5, size = 2.00, colour = dark
  ) +
  annotate(
    "text", x = 0.02, y = 0.20,
    label = "Effective resource state is latent",
    hjust = 0, vjust = 0.5, size = 1.95, colour = "#666666"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 7, 0, 13))

p_a <- p_design_ivt + p_o2_key + plot_spacer() +
  p_design_vivo + p_vivo_key +
  plot_layout(
    ncol = 1,
    heights = c(1, 0.12, 0.06, 0.62, 0.10)
  )

# -----------------------------------------------------------------------------
# B: observed in-vitro chromosome-number trajectories
# -----------------------------------------------------------------------------

kary_plot <- kary |>
  mutate(
    condition = condition_from_segment(segment_id),
    cohort = factor(cohort, levels = c("2N", "4N"))
  )

start_rows <- kary_plot |>
  group_by(cohort) |>
  filter(cumulative_end_day == min(cumulative_end_day)) |>
  ungroup()

kary_summary <- kary_plot |>
  group_by(cohort, condition, cumulative_end_day) |>
  summarise(
    median_chr = median(observed_kary_N),
    q1 = quantile(observed_kary_N, 0.25),
    q3 = quantile(observed_kary_N, 0.75),
    .groups = "drop"
  )

start_summary <- start_rows |>
  group_by(cohort, cumulative_end_day) |>
  summarise(
    median_chr = median(observed_kary_N),
    q1 = quantile(observed_kary_N, 0.25),
    q3 = quantile(observed_kary_N, 0.75),
    .groups = "drop"
  )

trajectory_summary <- bind_rows(
  kary_summary,
  start_summary |> mutate(condition = "oxygen-deprived")
) |>
  distinct(cohort, condition, cumulative_end_day, .keep_all = TRUE) |>
  mutate(condition = factor(condition, levels = c("control", "oxygen-deprived")))

p_b <- ggplot() +
  geom_hline(
    data = data.frame(
      cohort = factor(c("2N", "4N"), levels = c("2N", "4N")),
      reference = c(44, 88)
    ),
    aes(yintercept = reference, colour = cohort),
    linewidth = 0.38, linetype = 3, alpha = 0.65
  ) +
  geom_point(
    data = kary_plot,
    aes(
      x = cumulative_end_day, y = observed_kary_N,
      colour = cohort
    ),
    position = position_jitter(width = 0.75, height = 0, seed = 23),
    size = 0.86, alpha = 0.20
  ) +
  geom_errorbar(
    data = trajectory_summary,
    aes(
      x = cumulative_end_day, ymin = q1, ymax = q3,
      colour = cohort, group = condition
    ),
    width = 2.2, linewidth = 0.40, alpha = 0.66
  ) +
  geom_line(
    data = trajectory_summary,
    aes(
      x = cumulative_end_day, y = median_chr,
      colour = cohort, linetype = condition, group = condition
    ),
    linewidth = 0.54
  ) +
  geom_point(
    data = trajectory_summary,
    aes(
      x = cumulative_end_day, y = median_chr,
      colour = cohort, shape = condition
    ),
    size = 1.65, stroke = 0.48, fill = "white"
  ) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = ploidy_cols, guide = "none") +
  scale_fill_manual(values = ploidy_cols, guide = "none") +
  scale_linetype_manual(
    values = c(control = "dashed", `oxygen-deprived` = "solid"),
    name = "Lineage"
  ) +
  scale_shape_manual(
    values = c(control = 22, `oxygen-deprived` = 21),
    name = "Lineage"
  ) +
  scale_x_continuous(
    breaks = c(0, 75, 150),
    limits = c(0, 153), expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    breaks = c(44, 66, 88, 110),
    limits = c(32, 112), expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    tag = "B",
    title = "Observed in-vitro chromosome trajectories",
    x = "Experimental day", y = "Chromosome number"
  ) +
  theme_manuscript() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.title = element_text(size = 6.3, face = "bold", hjust = 0),
    panel.spacing = grid::unit(5, "mm")
  ) +
  guides(
    linetype = guide_legend(order = 1, ncol = 1, byrow = TRUE),
    shape = guide_legend(order = 1, ncol = 1, byrow = TRUE)
  )

# -----------------------------------------------------------------------------
# C: starting cell-line reference versus terminal tumors
# -----------------------------------------------------------------------------

baseline <- start_rows |>
  transmute(
    cohort,
    stage = "Starting\nreference",
    chromosome_number = observed_kary_N,
    harvest = NA_character_
  ) |>
  mutate(stage = factor(stage, levels = c("Starting\nreference", "Terminal\ntumors")))

terminal <- invivo_ploidy |>
  filter(used_in_endpoint_loss, is.finite(endpoint_value)) |>
  transmute(
    cohort = factor(cohort, levels = c("2N", "4N")),
    stage = "Terminal\ntumors",
    chromosome_number = endpoint_value,
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
    chromosome_number = median(chromosome_number),
    .groups = "drop"
  ) |>
  mutate(
    stage = factor(
      stage,
      levels = c("Starting\nreference", "Terminal\ntumors")
    ),
    stage_x = as.numeric(stage)
  )

cohort_counts <- harvest |>
  filter(included_under_config, has_ploidy) |>
  distinct(cohort, harvest) |>
  count(cohort, name = "n_tumors") |>
  mutate(
    cohort = factor(cohort, levels = c("2N", "4N")),
    label = paste0("n = ", n_tumors, " tumors"),
    x = 2.38, y = 109
  )

make_left_half_violin <- function(data, width = 0.29, adjust = 0.85) {
  data |>
    group_by(cohort, stage) |>
    group_modify(function(.x, .y) {
      values <- .x$chromosome_number[
        is.finite(.x$chromosome_number)
      ]
      stage_x <- match(
        as.character(.y$stage[[1]]),
        levels(data$stage)
      )
      density_curve <- stats::density(
        values,
        n = 256,
        adjust = adjust,
        from = max(32, min(values) - 2),
        to = min(112, max(values) + 2)
      )
      scaled_density <- density_curve$y / max(density_curve$y)
      tibble(
        stage_x = c(
          stage_x,
          stage_x - width * scaled_density,
          stage_x
        ),
        chromosome_number = c(
          density_curve$x[[1]],
          density_curve$x,
          density_curve$x[[length(density_curve$x)]]
        )
      )
    }) |>
    ungroup()
}

half_violin <- make_left_half_violin(all_distributions)

p_c <- ggplot() +
  geom_hline(
    data = data.frame(
      cohort = factor(c("2N", "4N"), levels = c("2N", "4N")),
      reference = c(44, 88)
    ),
    aes(yintercept = reference, colour = cohort),
    linewidth = 0.38, linetype = 3, alpha = 0.65
  ) +
  geom_polygon(
    data = half_violin,
    aes(
      x = stage_x, y = chromosome_number,
      fill = cohort, group = interaction(cohort, stage)
    ),
    alpha = 0.20, colour = dark, linewidth = 0.30
  ) +
  geom_boxplot(
    data = all_distributions,
    aes(
      x = stage_x, y = chromosome_number,
      group = interaction(cohort, stage)
    ),
    width = 0.105, outlier.shape = NA,
    fill = "white", colour = dark, alpha = 0.86,
    linewidth = 0.36
  ) +
  geom_point(
    data = baseline,
    aes(
      x = as.numeric(stage) + 0.18,
      y = chromosome_number, colour = cohort
    ),
    position = position_jitter(width = 0.038, height = 0, seed = 17),
    size = 0.72, alpha = 0.34
  ) +
  geom_point(
    data = harvest_medians,
    aes(
      x = stage_x + 0.18,
      y = chromosome_number, colour = cohort
    ),
    position = position_jitter(width = 0.038, height = 0, seed = 17),
    shape = 23, fill = "white", size = 1.72, stroke = 0.58
  ) +
  geom_text(
    data = cohort_counts,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, size = 2.50, colour = dark,
    hjust = 1, vjust = 1
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
    breaks = c(44, 66, 88, 110),
    limits = c(32, 112), expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    tag = "C",
    title = "Starting state and terminal tumors",
    x = NULL, y = "Chromosome number"
  ) +
  theme_manuscript() +
  theme(
    axis.text.x = element_text(size = 7, lineheight = 0.95),
    panel.spacing = grid::unit(4, "mm")
  )

# -----------------------------------------------------------------------------
# Save subpanels and assembled draft
# -----------------------------------------------------------------------------

save_both(
  p_a,
  file.path(subpanel_root, "fig1a_matched_design"),
  width = 4.733333, height = 5.325
)
save_both(
  p_b,
  file.path(subpanel_root, "fig1b_invitro_chromosome_trajectories"),
  width = 2.366667, height = 2.6625
)
save_both(
  p_c,
  file.path(subpanel_root, "fig1c_invivo_start_terminal_distributions"),
  width = 2.366667, height = 2.6625
)

right_column <- p_b / plot_spacer() / p_c +
  plot_layout(heights = c(1, 0.08, 1))

figure1 <- p_a | right_column
figure1 <- figure1 +
  plot_layout(widths = c(2, 1))

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
      estimate = median_chr, lower = q1, upper = q3
    ),
  harvest_medians |>
    transmute(
      panel = "F1C", cohort = as.character(cohort),
      group = harvest, time_day = NA_real_,
      estimate = chromosome_number, lower = NA_real_, upper = NA_real_
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
  plotting_contract = "iteration3_refined_layout_and_raincloud_distributions",
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
provenance <- data.frame(
  input_role = names(paths),
  local_path = unname(unlist(paths)),
  upstream_path = file.path(upstream_input_root, basename(unname(unlist(paths)))),
  stringsAsFactors = FALSE
)
write.table(
  provenance,
  file.path(data_root, "figure1_source_file_provenance.tsv"),
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
