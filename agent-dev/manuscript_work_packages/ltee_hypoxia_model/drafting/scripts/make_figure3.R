#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) {
  stop("MININGCLONEID_REPO_ROOT is not set. Run this script with scripts/agentRrunner.sh.")
}
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)

draft_root <- file.path(
  repo_root,
  "agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting"
)
seed_root <- file.path(
  draft_root,
  "source_tables/frozen_inputs/F3"
)

input_paths <- c(
  lineage = file.path(seed_root, "invitro_lineage_summary.tsv"),
  distribution = file.path(seed_root, "invitro_distribution_summary.tsv"),
  observed_kary = file.path(seed_root, "invitro_observed_kary.tsv"),
  survival_curve = file.path(seed_root, "viz/functional_curve_ploidy.tsv"),
  nonviable_curve = file.path(seed_root, "viz/functional_curve_oxygen_multi_ploidy.tsv")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing saved seed-10 input table(s): ", paste(missing_inputs, collapse = ", "))
}

output_dirs <- c(
  file.path(draft_root, "initial_subpanels/F3A"),
  file.path(draft_root, "initial_subpanels/F3B"),
  file.path(draft_root, "initial_subpanels/F3C"),
  file.path(draft_root, "initial_subpanels/F3D"),
  file.path(draft_root, "refined_subpanels"),
  file.path(draft_root, "final_figures/recommended")
)
invisible(vapply(output_dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))

read_saved_tsv <- function(path, required_columns) {
  x <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns)) {
    stop(
      basename(path), " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  x
}

lineage_df <- read_saved_tsv(
  input_paths[["lineage"]],
  c(
    "segment_id", "parent_segment_id", "cohort", "passage_index",
    "predicted_growth_rate", "observed_growth"
  )
)
distribution_df <- read_saved_tsv(
  input_paths[["distribution"]],
  c("segment_id", "cohort", "passage_index", "N", "fraction")
)
observed_kary_df <- read_saved_tsv(
  input_paths[["observed_kary"]],
  c("segment_id", "cohort", "passage_index", "passage_id", "cell_index", "observed_kary_N")
)
survival_df <- read_saved_tsv(
  input_paths[["survival_curve"]],
  c("N", "viability_after_ms")
)
nonviable_df <- read_saved_tsv(
  input_paths[["nonviable_curve"]],
  c("cohort", "ploidy_multiple", "ms_rate", "misseg_nonviable_daughter_fraction")
)

cohort_levels <- c("2N", "4N")
lineage_levels <- c("Control", "Deprived")
cohort_colors <- c("2N" = "#0072B2", "4N" = "#D55E00")
ploidy_colors <- c(
  "1.5N" = "#CC79A7",
  "2N" = "#0072B2",
  "2.5N" = "#E69F00",
  "3N" = "#009E73",
  "3.5N" = "#56B4E9",
  "4N" = "#D55E00",
  "4.5N" = "#F0E442",
  "5N" = "#000000"
)

is_control_path <- function(segment_id) {
  tokens <- strsplit(as.character(segment_id), "_", fixed = TRUE)
  vapply(
    tokens,
    function(x) {
      oxygen <- suppressWarnings(as.numeric(x))
      length(oxygen) > 0L && all(is.finite(oxygen)) && all(abs(oxygen - 20.5) < 1e-8)
    },
    logical(1)
  )
}

add_lineage_coordinates <- function(x) {
  segment_id <- as.character(x$segment_id)
  passage_from_path <- lengths(strsplit(segment_id, "_", fixed = TRUE))
  passage_index <- suppressWarnings(as.integer(x$passage_index))
  if (any(!is.finite(passage_index)) || any(passage_from_path < 1L)) {
    stop("Saved rows contain an invalid passage coordinate or segment lineage.")
  }
  x |>
    mutate(
      cohort = factor(.data$cohort, levels = cohort_levels),
      lineage = factor(
        ifelse(is_control_path(.data$segment_id), "Control", "Deprived"),
        levels = lineage_levels
      ),
      lineage_passage = passage_from_path
    )
}

assert_four_strata <- function(x, table_name) {
  observed <- unique(data.frame(
    cohort = as.character(x$cohort),
    lineage = as.character(x$lineage),
    stringsAsFactors = FALSE
  ))
  expected <- expand.grid(
    cohort = cohort_levels,
    lineage = lineage_levels,
    stringsAsFactors = FALSE
  )
  key <- function(z) paste(z$cohort, z$lineage, sep = "::")
  missing <- setdiff(key(expected), key(observed))
  if (length(missing)) {
    stop(table_name, " lacks required control/deprived strata: ", paste(missing, collapse = ", "))
  }
}

theme_manuscript <- function() {
  theme_minimal(base_size = 8.5, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E4E4E4", linewidth = 0.25),
      plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.3, colour = "#333333", margin = margin(b = 3)),
      plot.title.position = "plot",
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = "topleft",
      axis.title = element_text(size = 8.4),
      axis.text = element_text(size = 7.2, colour = "#222222"),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "#F1F1F1", colour = "#D0D0D0", linewidth = 0.3),
      legend.title = element_text(size = 7.6, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(4.5, "mm"),
      plot.margin = margin(5, 5, 4, 6)
    )
}

save_plot_pair <- function(plot, stem, width, height) {
  png_path <- paste0(stem, ".png")
  pdf_path <- paste0(stem, ".pdf")
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  c(png_path, pdf_path)
}

# Panel A: saved observed and fitted growth rates over lineage passage.
lineage_plot_df <- add_lineage_coordinates(lineage_df)
assert_four_strata(lineage_plot_df, basename(input_paths[["lineage"]]))

growth_nodes <- lineage_plot_df |>
  filter(is.finite(.data$predicted_growth_rate)) |>
  group_by(
    .data$segment_id, .data$parent_segment_id, .data$cohort,
    .data$lineage, .data$lineage_passage
  ) |>
  summarise(
    predicted_growth_rate = mean(.data$predicted_growth_rate),
    prediction_spread = diff(range(.data$predicted_growth_rate)),
    .groups = "drop"
  )
if (any(growth_nodes$prediction_spread > 1e-10, na.rm = TRUE)) {
  stop("Replicated rows disagree on a segment's saved predicted growth rate.")
}

growth_parents <- growth_nodes |>
  transmute(
    cohort = .data$cohort,
    lineage = .data$lineage,
    parent_segment_id = .data$segment_id,
    parent_passage = .data$lineage_passage,
    parent_growth = .data$predicted_growth_rate
  )
growth_edges <- growth_nodes |>
  left_join(
    growth_parents,
    by = c("cohort", "lineage", "parent_segment_id"),
    relationship = "many-to-one"
  ) |>
  filter(is.finite(.data$parent_passage), is.finite(.data$parent_growth))

growth_observed <- lineage_plot_df |>
  filter(is.finite(.data$observed_growth))

panel_a <- ggplot() +
  geom_segment(
    data = growth_edges,
    aes(
      x = .data$parent_passage,
      y = .data$parent_growth,
      xend = .data$lineage_passage,
      yend = .data$predicted_growth_rate,
      colour = .data$cohort
    ),
    linewidth = 0.72,
    lineend = "round"
  ) +
  geom_point(
    data = growth_observed,
    aes(
      x = .data$lineage_passage,
      y = .data$observed_growth,
      fill = .data$cohort
    ),
    shape = 21,
    colour = "#202020",
    stroke = 0.28,
    size = 2,
    alpha = 0.78,
    position = position_jitter(width = 0.10, height = 0, seed = 10)
  ) +
  facet_grid(
    rows = vars(cohort),
    cols = vars(lineage),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_colour_manual(values = cohort_colors, drop = FALSE) +
  scale_fill_manual(values = cohort_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = function(x) pretty(x, n = 6),
    expand = expansion(mult = c(0.025, 0.035))
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0.06, 0.10))
  ) +
  labs(
    tag = "A",
    title = "Growth-rate fit",
    subtitle = "Points: observed \u2022 lines: fitted model \u2022 selected seed 10",
    x = "Lineage passage",
    y = expression(paste("Growth rate (day"^{-1}, ")"))
  ) +
  guides(colour = "none", fill = "none") +
  theme_manuscript()

# Panel B: saved predicted chromosome-state fractions plus every observed karyotype.
distribution_plot_df <- add_lineage_coordinates(distribution_df) |>
  filter(
    is.finite(.data$N),
    is.finite(.data$fraction),
    .data$fraction >= 0
  ) |>
  mutate(N = as.integer(round(.data$N))) |>
  group_by(.data$cohort, .data$lineage, .data$lineage_passage, .data$N) |>
  summarise(fraction = mean(.data$fraction), .groups = "drop")
assert_four_strata(distribution_plot_df, basename(input_paths[["distribution"]]))

observed_kary_plot_df <- add_lineage_coordinates(observed_kary_df) |>
  filter(is.finite(.data$observed_kary_N)) |>
  arrange(
    .data$cohort, .data$lineage, .data$lineage_passage,
    .data$passage_id, .data$cell_index
  ) |>
  group_by(.data$cohort, .data$lineage, .data$lineage_passage) |>
  mutate(
    observation_index = row_number(),
    observations_at_passage = n(),
    x_observed = .data$lineage_passage +
      0.68 * (.data$observation_index - (.data$observations_at_passage + 1) / 2) /
      max(.data$observations_at_passage, 1)
  ) |>
  ungroup()
assert_four_strata(observed_kary_plot_df, basename(input_paths[["observed_kary"]]))

fraction_fill_max <- max(distribution_plot_df$fraction, na.rm = TRUE)
if (!is.finite(fraction_fill_max) || fraction_fill_max <= 0) {
  stop("Predicted chromosome-state fractions have no positive finite values.")
}

panel_b <- ggplot(
  distribution_plot_df,
  aes(x = .data$lineage_passage, y = .data$N, fill = .data$fraction)
) +
  geom_tile(width = 0.96, height = 1) +
  geom_point(
    data = observed_kary_plot_df,
    aes(x = .data$x_observed, y = .data$observed_kary_N),
    inherit.aes = FALSE,
    shape = 21,
    size = 0.82,
    stroke = 0.22,
    colour = "#111111",
    fill = "white",
    alpha = 0.92
  ) +
  facet_grid(
    rows = vars(cohort),
    cols = vars(lineage),
    scales = "free",
    space = "free_x"
  ) +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#2C7FB8", "#FFFF33"),
    values = rescale(c(0, 0.05 * fraction_fill_max, fraction_fill_max)),
    limits = c(0, fraction_fill_max),
    oob = squish,
    na.value = "white",
    name = "Predicted\nfraction"
  ) +
  scale_x_continuous(
    breaks = function(x) pretty(x, n = 6),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    barheight = grid::unit(24, "mm"),
    barwidth = grid::unit(3.2, "mm")
  )) +
  labs(
    tag = "B",
    title = "Chromosome-state distributions",
    subtitle = "Heatmap: fitted state fraction \u2022 open points: observed karyotypes",
    x = "Lineage passage",
    y = "Chromosome count"
  ) +
  theme_manuscript() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Panel C: saved fitted post-missegregation survival curve.
survival_plot_df <- survival_df |>
  transmute(
    N = as.numeric(.data$N),
    survival = as.numeric(.data$viability_after_ms)
  ) |>
  filter(is.finite(.data$N), is.finite(.data$survival)) |>
  arrange(.data$N)
if (!nrow(survival_plot_df)) {
  stop("No finite rows in the saved post-missegregation survival curve.")
}

reference_states <- data.frame(
  cohort = factor(cohort_levels, levels = cohort_levels),
  N = c(44, 88)
) |>
  rowwise() |>
  mutate(
    survival = survival_plot_df$survival[
      which.min(abs(survival_plot_df$N - .data$N))
    ]
  ) |>
  ungroup()

panel_c <- ggplot(survival_plot_df, aes(x = .data$N, y = .data$survival)) +
  geom_vline(
    data = reference_states,
    aes(xintercept = .data$N, colour = .data$cohort),
    linewidth = 0.48,
    linetype = "22",
    alpha = 0.85
  ) +
  geom_line(linewidth = 0.85, colour = "#222222", lineend = "round") +
  geom_point(
    data = reference_states,
    aes(x = .data$N, y = .data$survival, fill = .data$cohort),
    inherit.aes = FALSE,
    shape = 21,
    colour = "white",
    stroke = 0.35,
    size = 2.2
  ) +
  geom_text(
    data = reference_states,
    aes(
      x = .data$N,
      y = .data$survival,
      label = .data$cohort,
      colour = .data$cohort
    ),
    inherit.aes = FALSE,
    nudge_x = 5,
    nudge_y = 0.035,
    size = 2.55,
    fontface = "bold",
    hjust = 0
  ) +
  scale_colour_manual(values = cohort_colors, drop = FALSE) +
  scale_fill_manual(values = cohort_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    tag = "C",
    title = "Post-missegregation survival",
    x = "Chromosome number",
    y = "Survival probability"
  ) +
  guides(colour = "none", fill = "none") +
  theme_manuscript()

# Panel D: saved model-implied fraction of nonviable daughters.
nonviable_plot_df <- nonviable_df |>
  transmute(
    cohort = factor(.data$cohort, levels = names(ploidy_colors)),
    ploidy_multiple = as.numeric(.data$ploidy_multiple),
    ms_rate = as.numeric(.data$ms_rate),
    nonviable_fraction = as.numeric(.data$misseg_nonviable_daughter_fraction)
  ) |>
  filter(
    !is.na(.data$cohort),
    is.finite(.data$ms_rate),
    is.finite(.data$nonviable_fraction)
  ) |>
  arrange(.data$cohort, .data$ms_rate)
if (!all(names(ploidy_colors) %in% as.character(unique(nonviable_plot_df$cohort)))) {
  stop("Saved nonviable-daughter table does not contain all eight reference states.")
}

panel_d <- ggplot(
  nonviable_plot_df,
  aes(
    x = .data$ms_rate,
    y = .data$nonviable_fraction,
    colour = .data$cohort,
    group = .data$cohort
  )
) +
  geom_line(linewidth = 0.82, lineend = "round") +
  scale_colour_manual(values = ploidy_colors, drop = FALSE, name = "Reference state") +
  scale_x_continuous(
    labels = label_percent(accuracy = 0.01),
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0.025, 0.03))
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  guides(colour = guide_legend(
    title.position = "top",
    nrow = 2,
    byrow = TRUE,
    override.aes = list(linewidth = 1.1)
  )) +
  labs(
    tag = "D",
    title = "Nonviable daughters after missegregation",
    x = "Missegregation probability per chromosome",
    y = "Nonviable daughter fraction"
  ) +
  theme_manuscript() +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(t = -1, r = 0, b = 0, l = 0)
  )

initial_outputs <- c(
  save_plot_pair(
    panel_a,
    file.path(draft_root, "initial_subpanels/F3A/figure3A_initial"),
    width = 7.1,
    height = 2.6
  ),
  save_plot_pair(
    panel_b,
    file.path(draft_root, "initial_subpanels/F3B/figure3B_initial"),
    width = 7.1,
    height = 3.6
  ),
  save_plot_pair(
    panel_c,
    file.path(draft_root, "initial_subpanels/F3C/figure3C_initial"),
    width = 3.45,
    height = 2.7
  ),
  save_plot_pair(
    panel_d,
    file.path(draft_root, "initial_subpanels/F3D/figure3D_initial"),
    width = 3.45,
    height = 3
  )
)

bottom_row <- panel_c + panel_d + plot_layout(widths = c(1, 1))
figure_3 <- panel_a / panel_b / bottom_row +
  plot_layout(heights = c(2.45, 3.35, 2.55))

refined_outputs <- save_plot_pair(
  figure_3,
  file.path(draft_root, "refined_subpanels/figure3_refined"),
  width = 7.1,
  height = 8.8
)
final_outputs <- save_plot_pair(
  figure_3,
  file.path(draft_root, "final_figures/recommended/figure3"),
  width = 7.1,
  height = 8.8
)

message("Saved Figure 3 from immutable seed-10 tables.")
message("Input row counts:")
message("  lineage: ", nrow(lineage_df))
message("  distribution: ", nrow(distribution_df))
message("  observed karyotypes: ", nrow(observed_kary_df))
message("  survival curve: ", nrow(survival_df))
message("  nonviable-daughter curve: ", nrow(nonviable_df))
message("Outputs:")
message(paste(c(initial_outputs, refined_outputs, final_outputs), collapse = "\n"))
