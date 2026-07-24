#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = getwd())
draft_root <- file.path(
  repo_root, "agent-dev", "manuscript_work_packages",
  "ltee_hypoxia_model", "drafting"
)

paths <- list(
  timeline = file.path(draft_root, "source_tables/frozen_inputs/F1/invitro_lineage_timeline.tsv"),
  passage = file.path(draft_root, "source_tables/frozen_inputs/F1/invitro_passage_observations.tsv"),
  kary = file.path(draft_root, "source_tables/frozen_inputs/F1/invitro_kary_cells.tsv"),
  burden = file.path(draft_root, "source_tables/frozen_inputs/F1/invivo_burden_long.tsv"),
  invivo_ploidy = file.path(draft_root, "source_tables/frozen_inputs/F1/invivo_ploidy_cells.tsv"),
  harvest = file.path(draft_root, "source_tables/frozen_inputs/F1/invivo_harvest_catalog.tsv")
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

ploidy_cols <- c("2N" = "#0072B2", "4N" = "#D55E00")
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
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
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
      ymin = track_y - 0.24, ymax = track_y + 0.24,
      fill = oxygen_pct
    ),
    colour = "white", linewidth = 0.22
  ) +
  geom_point(
    data = flow_events,
    aes(x = cumulative_end_day, y = track_y + 0.35),
    shape = 24, size = 2.2, stroke = 0.45,
    fill = "#D9EAF5", colour = "#356B96"
  ) +
  geom_point(
    data = kary_events,
    aes(x = cumulative_end_day, y = track_y + 0.35),
    shape = 21, size = 1.8, stroke = 0.45,
    fill = "#F6C978", colour = "#8A5A13"
  ) +
  scale_fill_gradientn(
    colours = oxygen_cols,
    values = rescale(c(0, 1, 5, 20.5)),
    limits = c(0, 20.5),
    breaks = c(0, 5, 20.5),
    labels = c("0", "5", "20.5"),
    name = expression("Target O"[2]*" (%)")
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
    x = "In-vitro experimental day", y = NULL
  ) +
  annotate(
    "text", x = 151, y = 0.55,
    label = "\u25B2 flow   \u25CF chromosome count",
    hjust = 1, vjust = 0, size = 2.50, colour = dark
  ) +
  theme_manuscript() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.key.width = grid::unit(10, "mm"),
    plot.margin = margin(6, 6, 5, 13)
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
  geom_segment(
    aes(x = 0, xend = max_day, y = y, yend = y),
    linewidth = 1.05, lineend = "round"
  ) +
  geom_point(
    data = burden_days,
    aes(x = day, y = y),
    size = 1.15, alpha = 0.75
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
  annotate(
    "text", x = 92, y = 0.56,
    label = "\u25CF burden measurements\n\u25C7 terminal chromosome count",
    hjust = 1, vjust = 0, size = 2.50, colour = dark, lineheight = 1.1
  ) +
  annotate(
    "text", x = 46, y = 0.22,
    label = "Effective resource state is latent",
    size = 2.50, colour = "#666666"
  ) +
  scale_colour_manual(values = ploidy_cols, guide = "none") +
  scale_y_continuous(
    breaks = c(1, 2), labels = c("4N tumors", "2N tumors"),
    limits = c(0.05, 2.55)
  ) +
  scale_x_continuous(
    breaks = c(0, 30, 60, 90), limits = c(0, 94),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "In-vivo experimental day", y = NULL) +
  theme_manuscript() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(24, 7, 5, 6)
  )

p_a <- p_design_ivt + p_design_vivo +
  plot_layout(widths = c(1.72, 1))

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
    position = position_jitter(width = 1.25, height = 0, seed = 23),
    size = 1.15, alpha = 0.28
  ) +
  geom_ribbon(
    data = trajectory_summary,
    aes(
      x = cumulative_end_day, ymin = q1, ymax = q3,
      group = condition, fill = cohort
    ),
    alpha = 0.11, colour = NA
  ) +
  geom_line(
    data = trajectory_summary,
    aes(
      x = cumulative_end_day, y = median_chr,
      colour = cohort, linetype = condition, group = condition
    ),
    linewidth = 0.72
  ) +
  geom_point(
    data = trajectory_summary,
    aes(
      x = cumulative_end_day, y = median_chr,
      colour = cohort, shape = condition
    ),
    size = 2.05, stroke = 0.55, fill = "white"
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
    breaks = c(0, 30, 60, 90, 120, 150),
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
    panel.spacing = grid::unit(5, "mm")
  ) +
  guides(
    linetype = guide_legend(order = 1),
    shape = guide_legend(order = 1)
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
    )
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
    )
  )

cohort_counts <- harvest |>
  filter(included_under_config, has_ploidy) |>
  distinct(cohort, harvest) |>
  count(cohort, name = "n_tumors") |>
  mutate(
    cohort = factor(cohort, levels = c("2N", "4N")),
    label = paste0("n = ", n_tumors, " tumors"),
    x = 2, y = 109
  )

p_c <- ggplot(all_distributions, aes(stage, chromosome_number)) +
  geom_hline(
    data = data.frame(
      cohort = factor(c("2N", "4N"), levels = c("2N", "4N")),
      reference = c(44, 88)
    ),
    aes(yintercept = reference, colour = cohort),
    linewidth = 0.38, linetype = 3, alpha = 0.65
  ) +
  geom_violin(
    data = terminal,
    aes(fill = cohort, group = interaction(cohort, stage)),
    width = 0.78, trim = TRUE, alpha = 0.18,
    colour = dark, linewidth = 0.35
  ) +
  geom_boxplot(
    aes(fill = cohort),
    width = 0.28, outlier.shape = NA, alpha = 0.38,
    linewidth = 0.4
  ) +
  geom_point(
    data = baseline,
    aes(colour = cohort),
    position = position_jitter(width = 0.085, height = 0, seed = 17),
    size = 1.2, alpha = 0.62
  ) +
  geom_point(
    data = harvest_medians,
    aes(colour = cohort),
    position = position_jitter(width = 0.075, height = 0, seed = 17),
    shape = 23, fill = "white", size = 2.15, stroke = 0.65
  ) +
  geom_text(
    data = cohort_counts,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, size = 2.50, colour = dark, vjust = 1
  ) +
  facet_wrap(~cohort, nrow = 1) +
  scale_x_discrete(limits = c("Starting\nreference", "Terminal\ntumors")) +
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
  file.path(draft_root, "initial_subpanels", "F1A", "matched_design"),
  width = 7.1, height = 2.45
)
save_both(
  p_b,
  file.path(draft_root, "initial_subpanels", "F1B", "invitro_chromosome_trajectories"),
  width = 4.2, height = 3.25
)
save_both(
  p_c,
  file.path(draft_root, "initial_subpanels", "F1C", "invivo_start_terminal_distributions"),
  width = 3.0, height = 3.25
)

figure1 <- p_a / (p_b + p_c + plot_layout(widths = c(1.42, 1))) +
  plot_layout(heights = c(0.78, 1))

save_both(
  figure1,
  file.path(draft_root, "refined_subpanels", "figure1_recommended"),
  width = 7.1, height = 6.05
)
save_both(
  figure1,
  file.path(draft_root, "final_figures", "recommended", "figure1"),
  width = 7.1, height = 6.05
)

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
  file.path(draft_root, "source_tables", "figure1_displayed_summaries.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Figure 1 draft written under: ", draft_root)
