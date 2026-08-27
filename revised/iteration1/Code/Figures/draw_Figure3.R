#!/usr/bin/env Rscript

if (identical(Sys.getenv("FIGURE3_DRAW_WORKER"), "1")) {

# Regenerate the data-backed panels used by iteration3 Figure 3.
# No parameter fitting is performed. All panels use the retained seed-10
# separate in-vitro fit and its materialized post-fit outputs.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
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
seed_root <- normalizePath(
  Sys.getenv(
    "FIGURE3_SEED_ROOT",
    unset = file.path(
      workspace_root, "data", "Figures",
      "Figure3", "source_seed10"
    )
  ),
  mustWork = TRUE
)
data_dir <- normalizePath(
  Sys.getenv(
    "FIGURE3_DATA_DIR",
    unset = file.path(
      workspace_root, "data", "Figures", "Figure3"
    )
  ),
  mustWork = FALSE
)
figure_dir <- normalizePath(
  Sys.getenv(
    "FIGURE3_PANEL_DIR",
    unset = file.path(
      workspace_root, "Figures", "Figure3_panels"
    )
  ),
  mustWork = FALSE
)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

flow_path <- file.path(seed_root, "invitro_flow_overlay.tsv")
distribution_path <- file.path(seed_root, "invitro_distribution_summary.tsv")
daily_path <- file.path(seed_root, "invitro_daily_counts.tsv")
lineage_path <- file.path(seed_root, "invitro_lineage_summary.tsv")
growth_path <- file.path(seed_root, "invitro_growth_loglik.tsv")
observed_kary_path <- file.path(seed_root, "invitro_observed_kary.tsv")
parameter_path <- file.path(seed_root, "best_params.tsv")
model_path <- file.path(
  repo_root, "oxygen", "code", "O2_supply_demand_MAP",
  "model", "model_O2_supply_demand_MAP.R"
)
util_root <- normalizePath(file.path(script_dir, "util"), mustWork = TRUE)
o2_grid_path <- file.path(util_root, "oxygen", "adaptive_grids.R")
functional_path <- file.path(util_root, "oxygen", "functional_response.R")
required_paths <- c(
  flow_path,
  distribution_path,
  daily_path,
  lineage_path,
  growth_path,
  observed_kary_path,
  parameter_path,
  model_path,
  o2_grid_path,
  functional_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing Figure 3 input(s): ", paste(missing_paths, collapse = ", "))
}

flow <- read.delim(flow_path, check.names = FALSE, stringsAsFactors = FALSE)
distribution <- read.delim(
  distribution_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
daily <- read.delim(daily_path, check.names = FALSE, stringsAsFactors = FALSE)
lineage_summary <- read.delim(
  lineage_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
growth <- read.delim(
  growth_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
observed_kary <- read.delim(
  observed_kary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
parameter_table <- read.delim(
  parameter_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cohort_levels <- c("2N", "4N")
lineage_levels <- c("control", "deprived")
observed_color <- "#0072B2"
predicted_color <- "#D55E00"
ploidy_color <- "#2C7BB6"
direct_death_color <- "#C43C39"
transition_death_color <- "#2E7D5B"

path_tokens <- function(x) strsplit(as.character(x), "_", fixed = TRUE)
path_length <- function(x) lengths(path_tokens(x))
lineage_class <- function(x) {
  vapply(path_tokens(x), function(tokens) {
    oxygen <- suppressWarnings(as.numeric(tokens))
    if (all(is.finite(oxygen)) && all(oxygen >= 20.5)) "control" else "deprived"
  }, character(1L))
}
path_has_low_o2 <- function(x) {
  vapply(path_tokens(x), function(tokens) {
    oxygen <- suppressWarnings(as.numeric(tokens))
    any(is.finite(oxygen) & oxygen < 20.5)
  }, logical(1L))
}

theme_figure3 <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.8),
      plot.subtitle = element_text(size = base_size - 1, color = "#50555A"),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 1.2, color = "#202020"),
      strip.text = element_text(face = "bold", size = base_size - 1),
      strip.background = element_rect(
        fill = "#F2F3F5",
        color = "#BFC3C8",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E6E8EB", linewidth = 0.3),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 1.5),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_plot_set <- function(plot, stem, width, height) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".svg")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

# Figure 3A: passage growth only. The facet geometry matches Figure 3C so
# control/deprived columns and 2N/4N rows can be compared directly.
required_lineage <- c(
  "segment_id", "cohort", "predicted_growth_rate"
)
required_growth <- c(
  "segment_id", "cohort", "observed_growth"
)
if (!all(required_lineage %in% names(lineage_summary))) {
  stop("Lineage summary is missing Figure 3A prediction columns.")
}
if (!all(required_growth %in% names(growth))) {
  stop("Growth likelihood table is missing Figure 3A observation columns.")
}

growth_prediction_segment <- unique(
  lineage_summary[
    is.finite(lineage_summary$predicted_growth_rate) &
      lineage_summary$cohort %in% cohort_levels,
    c("segment_id", "cohort", "predicted_growth_rate"),
    drop = FALSE
  ]
)
growth_prediction_segment$lineage <- lineage_class(
  growth_prediction_segment$segment_id
)
growth_prediction_segment$lineage_passage_index <- path_length(
  growth_prediction_segment$segment_id
)
growth_prediction <- aggregate(
  predicted_growth_rate ~ cohort + lineage + lineage_passage_index,
  data = growth_prediction_segment,
  FUN = mean,
  na.rm = TRUE
)
growth_observed <- growth[
  is.finite(growth$observed_growth) &
    growth$cohort %in% cohort_levels,
  c("segment_id", "cohort", "observed_growth"),
  drop = FALSE
]
growth_observed$lineage <- lineage_class(growth_observed$segment_id)
growth_observed$lineage_passage_index <- path_length(
  growth_observed$segment_id
)
growth_prediction$cohort <- factor(
  growth_prediction$cohort,
  levels = cohort_levels
)
growth_prediction$lineage <- factor(
  growth_prediction$lineage,
  levels = lineage_levels
)
growth_observed$cohort <- factor(
  growth_observed$cohort,
  levels = cohort_levels
)
growth_observed$lineage <- factor(
  growth_observed$lineage,
  levels = lineage_levels
)
write.table(
  growth_prediction,
  file.path(data_dir, "figure3a_seed10_growth_prediction.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  growth_observed,
  file.path(data_dir, "figure3a_seed10_growth_observations.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

panel_a <- ggplot() +
  geom_line(
    data = growth_prediction,
    aes(
      lineage_passage_index,
      predicted_growth_rate,
      color = "Predicted",
      group = 1
    ),
    linewidth = 1.0,
    lineend = "round"
  ) +
  geom_point(
    data = growth_observed,
    aes(
      lineage_passage_index,
      observed_growth,
      color = "Observed"
    ),
    position = position_jitter(width = 0.08, height = 0, seed = 3101),
    shape = 16,
    size = 2.1,
    alpha = 0.86
  ) +
  facet_grid(
    cohort ~ lineage,
    scales = "free_x",
    space = "free_x",
    drop = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Observed" = observed_color,
      "Predicted" = predicted_color
    ),
    breaks = c("Observed", "Predicted"),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = function(limits) {
      candidate <- pretty(limits, n = 7)
      candidate[candidate >= limits[[1L]] & candidate <= limits[[2L]]]
    },
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Observed versus predicted passage growth rates",
    x = "Lineage passage",
    y = "Growth rate",
    color = "Series"
  ) +
  theme_figure3(11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major.x = element_line(color = "#ECEEF0", linewidth = 0.25)
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linetype = c("blank", "solid"),
        shape = c(16, NA_integer_),
        linewidth = c(0, 1)
      )
    )
  )

save_plot_set(
  panel_a,
  "fig3a_invitro_growth_rate_fit",
  width = 10,
  height = 5.25
)

# Figure 3B: flow-density overlays with the same observed/predicted palette as
# the Figure 4A terminal-distribution panel.
required_flow <- c(
  "segment_id", "cohort", "passage_index", "passage_id", "sample_name",
  "grid_index", "ploidy", "density", "series"
)
if (!all(required_flow %in% names(flow))) {
  stop("Flow overlay table is missing required columns.")
}
flow_plot <- flow[
  is.finite(flow$ploidy) &
    is.finite(flow$density) &
    flow$density >= 0 &
    flow$series %in% c("Observed", "Predicted"),
  ,
  drop = FALSE
]
flow_plot$cohort <- factor(flow_plot$cohort, levels = cohort_levels)
flow_plot$lineage <- lineage_class(flow_plot$segment_id)
flow_plot$lineage_passage_index <- path_length(flow_plot$segment_id)
flow_plot$parallel_line <- ifelse(
  grepl("_O1_", flow_plot$sample_name),
  "O1",
  ifelse(grepl("_O2_", flow_plot$sample_name), "O2", "")
)
flow_plot$sample_facet_label <- paste0(
  ifelse(
    flow_plot$lineage == "control",
    paste0("Ctrl p", flow_plot$lineage_passage_index),
    paste0(
      "p",
      flow_plot$lineage_passage_index,
      " ",
      flow_plot$parallel_line
    )
  )
)
flow_plot$series <- factor(
  flow_plot$series,
  levels = c("Observed", "Predicted")
)
flow_sample_panels <- unique(
  flow_plot[
    flow_plot$series == "Observed",
    c("cohort", "sample_facet_label", "sample_name"),
    drop = FALSE
  ]
)
flow_subpanel_sample_counts <- aggregate(
  sample_name ~ cohort + sample_facet_label,
  data = flow_sample_panels,
  FUN = function(x) length(unique(x))
)

make_flow_cohort_panel <- function(cohort_value, show_x_title) {
  cohort_data <- flow_plot[
    flow_plot$cohort == cohort_value,
    ,
    drop = FALSE
  ]
  panel_order <- unique(
    cohort_data[
      cohort_data$series == "Observed",
      c(
        "lineage", "lineage_passage_index", "parallel_line",
        "sample_facet_label"
      ),
      drop = FALSE
    ]
  )
  panel_order$lineage_order <- match(
    panel_order$lineage,
    lineage_levels
  )
  panel_order$parallel_order <- match(
    panel_order$parallel_line,
    c("", "O1", "O2")
  )
  panel_order <- panel_order[
    order(
      panel_order$lineage_order,
      panel_order$lineage_passage_index,
      panel_order$parallel_order
    ),
    ,
    drop = FALSE
  ]
  cohort_data$sample_facet_label <- factor(
    cohort_data$sample_facet_label,
    levels = panel_order$sample_facet_label
  )
  observed_data <- cohort_data[
    cohort_data$series == "Observed",
    ,
    drop = FALSE
  ]
  predicted_data <- cohort_data[
    cohort_data$series == "Predicted",
    ,
    drop = FALSE
  ]

  ggplot() +
    geom_line(
      data = observed_data,
      aes(
        ploidy,
        density,
        color = series,
        group = interaction(series, sample_name)
      ),
      linewidth = 0.55,
      linetype = "solid",
      lineend = "round"
    ) +
    geom_line(
      data = predicted_data,
      aes(
        ploidy,
        density,
        color = series,
        group = interaction(series, sample_name)
      ),
      linewidth = 0.55,
      linetype = "solid",
      lineend = "round"
    ) +
    facet_wrap(
      vars(sample_facet_label),
      nrow = 1,
      scales = "free_y",
      drop = TRUE
    ) +
    scale_color_manual(
      values = c(
        "Observed" = observed_color,
        "Predicted" = predicted_color
      ),
      drop = FALSE
    ) +
    scale_x_continuous(breaks = seq(2, 8, by = 2)) +
    labs(
      title = cohort_value,
      x = if (isTRUE(show_x_title)) "Ploidy" else NULL,
      y = "Density",
      color = "Series"
    ) +
    theme_figure3(10.5) +
    theme(
      plot.title = element_text(size = 11.5),
      strip.text.x = element_text(size = 5.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      panel.grid.major.x = element_line(
        color = "#ECEEF0",
        linewidth = 0.25
      ),
      plot.margin = margin(2, 7, 2, 7)
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 0.8, linetype = "solid")
      )
    )
}

panel_b_2n <- make_flow_cohort_panel("2N", FALSE) +
  theme(text = element_text(family = "sans"))
panel_b_4n <- make_flow_cohort_panel("4N", TRUE) +
  theme(text = element_text(family = "sans"))
panel_b_legend <- cowplot::get_legend(
  panel_b_2n +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)
panel_b_title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Observed versus predicted G0/G1 ploidy-density curves",
    x = 0.01,
    hjust = 0,
    fontfamily = "sans",
    fontface = "bold",
    size = 12.8
  )
panel_b_core <- cowplot::plot_grid(
  panel_b_2n + theme(legend.position = "none"),
  panel_b_4n + theme(legend.position = "none"),
  ncol = 1,
  rel_heights = c(1, 1)
)
panel_b <- cowplot::plot_grid(
  panel_b_title,
  panel_b_core,
  panel_b_legend,
  ncol = 1,
  rel_heights = c(0.08, 1, 0.10)
)

save_plot_set(
  panel_b,
  "fig3b_invitro_flow_density_fit",
  width = 10,
  height = 5.25
)

# Figure 3C: predicted chromosome-state heatmap with a bottom color bar.
required_distribution <- c(
  "segment_id", "cohort", "passage_index", "N", "fraction"
)
if (!all(required_distribution %in% names(distribution))) {
  stop("Distribution summary is missing required columns.")
}
distribution$lineage <- lineage_class(distribution$segment_id)
distribution$lineage_passage_index <- path_length(distribution$segment_id)
distribution_plot <- aggregate(
  fraction ~ cohort + lineage + lineage_passage_index + N,
  data = distribution,
  FUN = mean,
  na.rm = TRUE
)
distribution_plot <- distribution_plot[
  is.finite(distribution_plot$N) &
    distribution_plot$N <= 110 &
    is.finite(distribution_plot$fraction) &
    distribution_plot$fraction >= 0,
  ,
  drop = FALSE
]
distribution_plot$cohort <- factor(
  distribution_plot$cohort,
  levels = cohort_levels
)
distribution_plot$lineage <- factor(
  distribution_plot$lineage,
  levels = lineage_levels
)
distribution_plot$passage_display <- factor(
  distribution_plot$lineage_passage_index,
  levels = sort(unique(distribution_plot$lineage_passage_index))
)
fraction_fill_max <- max(distribution_plot$fraction, na.rm = TRUE)

distribution_mean <- aggregate(
  cbind(
    weighted_chromosome =
      distribution_plot$N * distribution_plot$fraction,
    distribution_mass = distribution_plot$fraction
  ) ~ cohort + lineage + lineage_passage_index,
  data = distribution_plot,
  FUN = sum,
  na.rm = TRUE
)
distribution_mean$predicted_mean_kary_N <-
  distribution_mean$weighted_chromosome /
  distribution_mean$distribution_mass
distribution_mean$passage_display <- factor(
  distribution_mean$lineage_passage_index,
  levels = levels(distribution_plot$passage_display)
)

required_observed_kary <- c(
  "segment_id", "cohort", "observed_kary_N"
)
if (!all(required_observed_kary %in% names(observed_kary))) {
  stop("Observed-karyotype table is missing Figure 3C overlay columns.")
}
observed_kary_plot <- observed_kary[
  observed_kary$cohort %in% cohort_levels &
    is.finite(observed_kary$observed_kary_N) &
    observed_kary$observed_kary_N <= 110,
  required_observed_kary,
  drop = FALSE
]
observed_kary_plot$lineage <- lineage_class(
  observed_kary_plot$segment_id
)
observed_kary_plot$lineage_passage_index <- path_length(
  observed_kary_plot$segment_id
)
observed_kary_plot$cohort <- factor(
  observed_kary_plot$cohort,
  levels = cohort_levels
)
observed_kary_plot$lineage <- factor(
  observed_kary_plot$lineage,
  levels = lineage_levels
)
observed_kary_plot$passage_display <- factor(
  observed_kary_plot$lineage_passage_index,
  levels = levels(distribution_plot$passage_display)
)
observed_kary_plot <- observed_kary_plot[
  !is.na(observed_kary_plot$passage_display),
  ,
  drop = FALSE
]
write.table(
  distribution_mean,
  file.path(data_dir, "figure3c_seed10_predicted_mean_chromosome.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  observed_kary_plot,
  file.path(data_dir, "figure3c_seed10_observed_karyotype_overlay.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

panel_c <- ggplot(
  distribution_plot,
  aes(passage_display, N, fill = fraction)
) +
  geom_tile(width = 0.98, height = 1) +
  geom_line(
    data = distribution_mean,
    aes(
      passage_display,
      predicted_mean_kary_N,
      group = interaction(cohort, lineage)
    ),
    inherit.aes = FALSE,
    color = predicted_color,
    linewidth = 0.95,
    lineend = "round"
  ) +
  geom_point(
    data = observed_kary_plot,
    aes(passage_display, observed_kary_N),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.10, height = 0, seed = 3102),
    shape = 21,
    size = 1.25,
    stroke = 0.25,
    color = "#202020",
    fill = predicted_color,
    alpha = 0.88
  ) +
  facet_grid(
    cohort ~ lineage,
    scales = "free_x",
    space = "free_x",
    drop = TRUE
  ) +
  scale_x_discrete(drop = TRUE, expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(min(distribution_plot$N) - 0.5, 110.5),
    expand = c(0, 0)
  ) +
  scale_fill_gradientn(
    colors = c("#F7F7F7", "#2C7FB8", "#FFFF33"),
    values = scales::rescale(
      c(0, 0.05 * fraction_fill_max, fraction_fill_max),
      from = c(0, fraction_fill_max)
    ),
    limits = c(0, fraction_fill_max),
    oob = scales::squish,
    na.value = "white",
    name = "Fraction"
  ) +
  labs(
    title = "Predicted ploidy distribution across passages",
    x = "Passage index",
    y = "Chromosome-count state N",
    fill = "Fraction"
  ) +
  theme_figure3(11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7.6)
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(10, "cm"),
      barheight = unit(0.38, "cm")
    )
  )

save_plot_set(
  panel_c,
  "fig3c_invitro_predicted_ploidy_distribution",
  width = 10,
  height = 5.25
)

# Figures 3D and 3E: reproduce the same seed-10 functional-response curves in
# narrow, manuscript-legible formats.
source(model_path, local = .GlobalEnv)
.first_non_null_local <- function(...) {
  values <- list(...)
  for (value in values) if (!is.null(value)) return(value)
  NULL
}
source(o2_grid_path, local = .GlobalEnv)
source(functional_path, local = .GlobalEnv)

run_params <- as.list(
  setNames(parameter_table$value, parameter_table$parameter)
)
run_params$boundary <- "drop"
run_params$ploidy_O2_death <- "diploid_NULL"
cfg <- list(
  N_MIN = 22L,
  N_MAX = 154L,
  N_UNIT = 22L,
  start_with = "chr_number",
  o2_crit_init = 1,
  n_O_init = 1,
  O2_growth = TRUE,
  ploidy_O2_death = "diploid_NULL",
  p_mis_base = 1e-5,
  p_mis_base_init = 1e-5,
  buffer_smax_init = 0.8,
  buffer_beta_init = 1,
  buffer_n_exp_init = 1,
  alpha_o2_init = 0.5,
  gamma_growth_init = 2,
  mu_hp_init = 0.001,
  gamma_mu_init = 2
)
functional_tables <- generate_invivo_functional_response_tables(
  run_params,
  cfg
)
viability_curve <- functional_tables$functional_curve_ploidy
nonviable_curve <- functional_tables$functional_curve_oxygen_multi_ploidy
nonviable_curve$cohort <- factor(
  nonviable_curve$cohort,
  levels = unique(nonviable_curve$cohort)
)
write.table(
  viability_curve,
  file.path(data_dir, "figure3d_seed10_post_missegregation_survival.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  nonviable_curve,
  file.path(data_dir, "figure3e_seed10_nonviable_fraction_vs_ms_rate.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ploidy_reference <- data.frame(
  endpoint_value = c(44, 88),
  ploidy_label = c("2N", "4N"),
  reference_color = c(observed_color, predicted_color),
  stringsAsFactors = FALSE
)
ploidy_reference$viability_after_ms <- approx(
  x = viability_curve$endpoint_value,
  y = viability_curve$viability_after_ms,
  xout = ploidy_reference$endpoint_value,
  ties = "ordered"
)$y
if (any(!is.finite(ploidy_reference$viability_after_ms))) {
  stop("Could not evaluate Figure 3D survival at the 2N/4N references.")
}

panel_d <- ggplot(
  viability_curve,
  aes(endpoint_value, viability_after_ms)
) +
  geom_vline(
    data = ploidy_reference,
    aes(xintercept = endpoint_value, color = reference_color),
    linetype = "22",
    linewidth = 0.8
  ) +
  geom_line(color = "#111111", linewidth = 1.05) +
  geom_point(
    data = ploidy_reference,
    aes(
      endpoint_value,
      viability_after_ms,
      color = reference_color
    ),
    inherit.aes = FALSE,
    shape = 16,
    size = 3.1
  ) +
  geom_text(
    data = ploidy_reference,
    aes(
      endpoint_value,
      viability_after_ms,
      label = ploidy_label,
      color = reference_color
    ),
    inherit.aes = FALSE,
    nudge_x = 3.0,
    nudge_y = 0.045,
    hjust = 0,
    vjust = 0,
    family = "Arial",
    fontface = "bold",
    size = 3.6
  ) +
  scale_color_identity(guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  labs(
    title = "Post-missegregation survival",
    subtitle = "One-copy chromosome-loss event",
    x = "Chromosome number (N)",
    y = "Viability after missegregation"
  ) +
  theme_figure3(11) +
  theme(legend.position = "none")

state_colors <- setNames(
  grDevices::hcl.colors(
    length(levels(nonviable_curve$cohort)),
    palette = "Dark 3"
  ),
  levels(nonviable_curve$cohort)
)
panel_e_labels <- do.call(
  rbind,
  lapply(levels(nonviable_curve$cohort), function(reference_state) {
    state_curve <- nonviable_curve[
      nonviable_curve$cohort == reference_state &
        is.finite(nonviable_curve$ms_rate) &
        is.finite(nonviable_curve$misseg_nonviable_daughter_fraction),
      ,
      drop = FALSE
    ]
    midpoint <- mean(range(state_curve$ms_rate))
    label_row <- state_curve[
      which.min(abs(state_curve$ms_rate - midpoint)),
      ,
      drop = FALSE
    ]
    label_row$curve_label <- reference_state
    label_row
  })
)
panel_e_labels$cohort <- factor(
  panel_e_labels$cohort,
  levels = levels(nonviable_curve$cohort)
)
panel_e <- ggplot(
  nonviable_curve,
  aes(
    ms_rate,
    misseg_nonviable_daughter_fraction,
    color = cohort,
    group = cohort
  )
) +
  geom_line(linewidth = 0.9) +
  geom_label(
    data = panel_e_labels,
    aes(
      ms_rate,
      misseg_nonviable_daughter_fraction,
      label = curve_label,
      color = cohort
    ),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 3.2,
    fill = alpha("white", 0.88),
    linewidth = 0,
    label.padding = unit(0.09, "lines"),
    label.r = unit(0.05, "lines"),
    show.legend = FALSE
  ) +
  scale_color_manual(values = state_colors, drop = FALSE) +
  labs(
    title = "Nonviable daughter fraction",
    subtitle = "Reference chromosome-number states",
    x = "Missegregation rate",
    y = "Nonviable daughters / all daughters"
  ) +
  theme_figure3(11) +
  theme(
    legend.position = "none"
  )

save_plot_set(
  panel_d,
  "fig3d_invitro_viability_after_missegregation",
  width = 5,
  height = 5.25
)
save_plot_set(
  panel_e,
  "fig3e_invitro_nonviable_fraction_vs_ms_rate",
  width = 5,
  height = 5.25
)

# Figure 3F: paired terminal severe-deprivation 2N and 4N lineages.
extract_terminal_branch <- function(cohort_value) {
  daily_endpoint <- daily[
    daily$cohort == cohort_value &
      daily$day == daily$selected_day,
    ,
    drop = FALSE
  ]
  candidate_paths <- unique(
    daily_endpoint$segment_id[path_has_low_o2(daily_endpoint$segment_id)]
  )
  if (!length(candidate_paths)) {
    stop("No oxygen-deprived lineage for cohort ", cohort_value)
  }
  lengths <- path_length(candidate_paths)
  terminal_candidates <- candidate_paths[lengths == max(lengths)]
  terminal_path <- sort(terminal_candidates)[[1L]]
  terminal_tokens <- path_tokens(terminal_path)[[1L]]
  terminal_prefixes <- vapply(
    seq_along(terminal_tokens),
    function(index) {
      paste(terminal_tokens[seq_len(index)], collapse = "_")
    },
    character(1L)
  )

  daily_branch <- daily_endpoint[
    daily_endpoint$segment_id %in% terminal_prefixes,
    ,
    drop = FALSE
  ]
  distribution_branch <- distribution[
    distribution$cohort == cohort_value &
      distribution$segment_id %in% terminal_prefixes,
    ,
    drop = FALSE
  ]
  weighted <- aggregate(
    cbind(
      weighted_chromosome =
        distribution_branch$N * distribution_branch$fraction,
      distribution_mass = distribution_branch$fraction
    ),
    by = list(segment_id = distribution_branch$segment_id),
    FUN = sum
  )
  weighted$predicted_mean_chromosome <-
    weighted$weighted_chromosome / weighted$distribution_mass
  branch <- merge(
    daily_branch,
    weighted[, c("segment_id", "predicted_mean_chromosome")],
    by = "segment_id",
    all = FALSE,
    sort = FALSE
  )
  branch$lineage_passage <- path_length(branch$segment_id)
  branch$live_burden_fraction <- branch$burden_live / branch$burden_total
  branch$direct_hypoxia_death_burden_fraction <-
    branch$burden_dead_hypoxia / branch$burden_total
  branch$nonviable_daughter_burden_fraction <-
    branch$burden_dead_buffer / branch$burden_total
  branch <- branch[
    order(branch$lineage_passage),
    c(
      "segment_id",
      "cohort",
      "lineage_passage",
      "passage_index",
      "oxygen_pct",
      "selected_day",
      "predicted_mean_chromosome",
      "live_burden_fraction",
      "direct_hypoxia_death_burden_fraction",
      "nonviable_daughter_burden_fraction"
    )
  ]
  if (
    nrow(branch) != length(terminal_prefixes) ||
      anyDuplicated(branch$lineage_passage) ||
      any(!is.finite(as.matrix(branch[, 7:10])))
  ) {
    stop("Terminal branch validation failed for cohort ", cohort_value)
  }
  branch
}

branch_2n <- extract_terminal_branch("2N")
branch_4n <- extract_terminal_branch("4N")
branch_summary <- rbind(branch_2n, branch_4n)
write.table(
  branch_summary,
  file.path(data_dir, "figure3f_seed10_2n_4n_deprived_passage_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

death_long <- rbind(
  data.frame(
    cohort = branch_summary$cohort,
    lineage_passage = branch_summary$lineage_passage,
    component = "Direct hypoxia-origin death",
    burden_fraction =
      branch_summary$direct_hypoxia_death_burden_fraction
  ),
  data.frame(
    cohort = branch_summary$cohort,
    lineage_passage = branch_summary$lineage_passage,
    component = "Nonviable chromosome-transition products",
    burden_fraction =
      branch_summary$nonviable_daughter_burden_fraction
  )
)
death_long$component <- factor(
  death_long$component,
  levels = c(
    "Direct hypoxia-origin death",
    "Nonviable chromosome-transition products"
  )
)
death_colors <- c(
  "Direct hypoxia-origin death" = direct_death_color,
  "Nonviable chromosome-transition products" = transition_death_color
)

make_ploidy_panel <- function(
  branch,
  cohort_value,
  show_y_title = TRUE
) {
  shaded <- branch[branch$oxygen_pct < 20.5, , drop = FALSE]
  initial_mean <- branch$predicted_mean_chromosome[[1L]]
  final_mean <- tail(branch$predicted_mean_chromosome, 1L)
  ggplot(
    branch,
    aes(lineage_passage, predicted_mean_chromosome)
  ) +
    geom_rect(
      data = shaded,
      aes(
        xmin = lineage_passage - 0.48,
        xmax = lineage_passage + 0.48,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "#DCEEF7",
      color = NA,
      alpha = 0.58
    ) +
    geom_line(color = ploidy_color, linewidth = 1.05) +
    geom_point(
      aes(fill = oxygen_pct < 20.5),
      shape = 21,
      size = 2.2,
      stroke = 0.5,
      color = ploidy_color
    ) +
    scale_fill_manual(
      values = c("FALSE" = "white", "TRUE" = ploidy_color),
      guide = "none"
    ) +
    annotate(
      "text",
      x = branch$lineage_passage[[1L]],
      y = initial_mean,
      label = sprintf("%.1f", initial_mean),
      hjust = -0.2,
      vjust = -0.8,
      color = ploidy_color,
      size = 3.1,
      family = "sans"
    ) +
    annotate(
      "text",
      x = tail(branch$lineage_passage, 1L),
      y = final_mean,
      label = sprintf("%.1f", final_mean),
      hjust = 1.2,
      vjust = if (final_mean >= initial_mean) -0.8 else 1.2,
      color = ploidy_color,
      size = 3.1,
      family = "sans"
    ) +
    scale_x_continuous(
      limits = c(1, 25),
      breaks = c(1, 5, 9, 13, 17, 21, 25),
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    scale_y_continuous(
      limits = c(40, 90),
      breaks = seq(40, 90, by = 10),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = paste0(cohort_value, " severe-deprivation lineage"),
      subtitle = "Predicted passage-end mean chromosome number; common scale",
      x = NULL,
      y = if (show_y_title) "Predicted mean chromosome\nnumber" else NULL
    ) +
    theme_figure3(10.5) +
    theme(
      legend.position = "none",
      text = element_text(family = "sans")
    )
}

make_death_panel <- function(cohort_value, show_y_title = TRUE) {
  cohort_death <- death_long[death_long$cohort == cohort_value, , drop = FALSE]
  direct_max <- max(
    cohort_death$burden_fraction[
      cohort_death$component == "Direct hypoxia-origin death"
    ],
    na.rm = TRUE
  )
  ggplot(
    cohort_death,
    aes(lineage_passage, burden_fraction, color = component)
  ) +
    geom_line(linewidth = 0.95) +
    geom_point(size = 1.75) +
    annotate(
      "label",
      x = 24.5,
      y = max(direct_max, 2e-4),
      label = sprintf("Direct-death max = %.2f%%", 100 * direct_max),
      hjust = 1,
      vjust = -0.55,
      size = 2.8,
      family = "sans",
      color = direct_death_color,
      fill = "white",
      linewidth = 0.2
    ) +
    scale_color_manual(values = death_colors, drop = FALSE) +
    scale_x_continuous(
      limits = c(1, 25),
      breaks = c(1, 5, 9, 13, 17, 21, 25),
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    scale_y_log10(
      limits = c(1e-5, 0.7),
      breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 0.3),
      labels = c("0.01%", "0.1%", "1%", "10%", "30%")
    ) +
    labs(
      x = "Lineage passage",
      y = if (show_y_title) "Fraction of total burden\n(log scale)" else NULL,
      color = NULL
    ) +
    theme_figure3(10.5) +
    theme(
      legend.position = "bottom",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      text = element_text(family = "sans")
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}

panel_f_ploidy_2n <- make_ploidy_panel(branch_2n, "2N", TRUE)
panel_f_ploidy_4n <- make_ploidy_panel(branch_4n, "4N", FALSE)
panel_f_death_2n <- make_death_panel("2N", TRUE)
panel_f_death_4n <- make_death_panel("4N", FALSE)
panel_f_legend <- cowplot::get_legend(
  panel_f_death_2n +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      text = element_text(family = "sans")
    )
)
panel_f_core <- cowplot::plot_grid(
  panel_f_ploidy_2n,
  panel_f_ploidy_4n,
  panel_f_death_2n + theme(legend.position = "none"),
  panel_f_death_4n + theme(legend.position = "none"),
  ncol = 2,
  rel_heights = c(1.02, 1)
)
panel_f <- cowplot::plot_grid(
  panel_f_core,
  panel_f_legend,
  ncol = 1,
  rel_heights = c(1, 0.11)
)

save_plot_set(
  panel_f,
  "fig3f_invitro_2n_4n_severe_deprivation",
  width = 20,
  height = 6.65
)

# Record exact external sources and machine-checkable panel contracts.
provenance <- data.frame(
  source = c(
    "flow_overlay",
    "distribution_summary",
    "daily_counts",
    "lineage_summary",
    "growth_likelihood",
    "observed_karyotype",
    "best_parameters",
    "model_core",
    "o2_grid_helpers",
    "functional_response_generator",
    "plotting_script"
  ),
  path = c(
    flow_path,
    distribution_path,
    daily_path,
    lineage_path,
    growth_path,
    observed_kary_path,
    parameter_path,
    model_path,
    o2_grid_path,
    functional_path,
    normalizePath(
      normalizePath(script_path <- sub(
        "^--file=", "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
      ), mustWork = TRUE),
      mustWork = TRUE
    )
  ),
  stringsAsFactors = FALSE
)
provenance$md5 <- unname(tools::md5sum(provenance$path))
write.table(
  provenance,
  file.path(data_dir, "figure3_source_file_provenance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

validation <- data.frame(
  metric = c(
    "growth_prediction_rows",
    "growth_observed_rows",
    "growth_observed_color",
    "growth_predicted_color",
    "panels_a_to_e_export_height_in",
    "panels_a_to_e_height_fraction_of_previous",
    "flow_rows",
    "flow_observed_color",
    "flow_predicted_color",
    "flow_observed_linetype",
    "flow_predicted_linetype",
    "flow_linewidth",
    "flow_layer_order",
    "flow_sample_subpanels",
    "flow_2n_sample_subpanels",
    "flow_4n_sample_subpanels",
    "flow_max_observed_samples_per_subpanel",
    "flow_parallel_lines_split",
    "distribution_rows",
    "distribution_predicted_mean_rows",
    "distribution_observed_karyotype_rows",
    "distribution_chromosome_overlay_color",
    "distribution_chromosome_overlay_z_order",
    "distribution_legend_position",
    "panel_d_curve_color",
    "panel_d_2n_reference_N",
    "panel_d_2n_reference_color",
    "panel_d_2n_intersection_viability",
    "panel_d_4n_reference_N",
    "panel_d_4n_reference_color",
    "panel_d_4n_intersection_viability",
    "panel_d_reference_legend_position",
    "panel_d_width_relative_to_panel_c",
    "panel_e_width_relative_to_panel_c",
    "panel_e_direct_label_count",
    "panel_e_legend_position",
    "f_2n_passages",
    "f_4n_passages",
    "f_shared_x_min_passage",
    "f_shared_x_max_passage",
    "f_shared_ploidy_y_min_chromosome",
    "f_shared_ploidy_y_max_chromosome",
    "f_4n_ploidy_y_title",
    "f_shared_legend_position",
    "f_export_height_in",
    "f_height_fraction_of_previous",
    "f_2n_initial_mean_chromosome",
    "f_2n_final_mean_chromosome",
    "f_4n_initial_mean_chromosome",
    "f_4n_final_mean_chromosome",
    "f_2n_max_direct_hypoxia_fraction",
    "f_4n_max_direct_hypoxia_fraction"
  ),
  value = c(
    nrow(growth_prediction),
    nrow(growth_observed),
    observed_color,
    predicted_color,
    5.25,
    0.7,
    nrow(flow_plot),
    observed_color,
    predicted_color,
    "solid",
    "solid",
    0.55,
    "observed_below_predicted",
    nrow(flow_sample_panels),
    sum(flow_sample_panels$cohort == "2N"),
    sum(flow_sample_panels$cohort == "4N"),
    max(flow_subpanel_sample_counts$sample_name),
    "O1_and_O2_separate",
    nrow(distribution_plot),
    nrow(distribution_mean),
    nrow(observed_kary_plot),
    predicted_color,
    "above_heatmap_tiles",
    "bottom",
    "#111111",
    ploidy_reference$endpoint_value[[1L]],
    ploidy_reference$reference_color[[1L]],
    ploidy_reference$viability_after_ms[[1L]],
    ploidy_reference$endpoint_value[[2L]],
    ploidy_reference$reference_color[[2L]],
    ploidy_reference$viability_after_ms[[2L]],
    "none",
    0.5,
    0.5,
    nrow(panel_e_labels),
    "none",
    nrow(branch_2n),
    nrow(branch_4n),
    1,
    25,
    40,
    90,
    "omitted",
    "bottom",
    6.65,
    0.7,
    branch_2n$predicted_mean_chromosome[[1L]],
    tail(branch_2n$predicted_mean_chromosome, 1L),
    branch_4n$predicted_mean_chromosome[[1L]],
    tail(branch_4n$predicted_mean_chromosome, 1L),
    max(branch_2n$direct_hypoxia_death_burden_fraction),
    max(branch_4n$direct_hypoxia_death_burden_fraction)
  ),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(data_dir, "figure3_panel_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Completed Figure 3 data-backed panel regeneration.")
message("  A: growth-rate-only observed/predicted fit")
message("  B: thin solid blue-observed/orange-predicted flow-density curves")
message("  C: chromosome-count overlay above the distribution heatmap")
message("  D: black survival curve with highlighted 2N/4N references")
message("  E: compact directly labeled functional-response curves")
message("  F: paired 2N/4N severe-deprivation trajectories")

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Figure3 <- function() {
  data_dir <- file.path(DATA_ROOT, "Figure3")
  panel_dir <- file.path(data_dir, "panels")
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_Figure3.R"), mustWork = TRUE),
    env = c(
      "FIGURE3_DRAW_WORKER=1",
      paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("FIGURE3_SEED_ROOT=", file.path(data_dir, "source_seed10")),
      paste0("FIGURE3_DATA_DIR=", data_dir),
      paste0("FIGURE3_PANEL_DIR=", panel_dir)
    )
  )
  compose_three_row_raster(
    data_dir = data_dir,
    panel_dir = panel_dir,
    output_dir = OUTPUT_ROOT,
    panels = c(
      A = file.path(panel_dir, "fig3a_invitro_growth_rate_fit.png"),
      B = file.path(panel_dir, "fig3b_invitro_flow_density_fit.png"),
      C = file.path(
        panel_dir, "fig3c_invitro_predicted_ploidy_distribution.png"
      ),
      D = file.path(
        panel_dir, "fig3d_invitro_viability_after_missegregation.png"
      ),
      E = file.path(
        panel_dir, "fig3e_invitro_nonviable_fraction_vs_ms_rate.png"
      ),
      F = file.path(
        panel_dir, "fig3f_invitro_2n_4n_severe_deprivation.png"
      )
    ),
    output_basename = "assembled_fig3",
    validation_basename = "figure3_layout_validation.tsv"
  )
  require_files(
    file.path(
      OUTPUT_ROOT,
      c("assembled_fig3.png", "assembled_fig3.pdf")
    ),
    "Figure 3 output"
  )
  invisible(file.path(OUTPUT_ROOT, "assembled_fig3.png"))
}

if (sys.nframe() == 0L) draw_Figure3()

}
