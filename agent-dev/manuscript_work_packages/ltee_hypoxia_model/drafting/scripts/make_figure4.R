#!/usr/bin/env Rscript

# Figure 4 drafting generator.
#
# Scope:
#   A. Compact seed25 observed-predicted tumor-burden and terminal
#      chromosome-number fit block.
#   B. Target versus latent model-implied effective-O2 trajectories.
#   C. O2 = 0/1/5 univariate separator triptych.
#   D. Preserved pooled in-vivo/in-vitro t-SNE rendering with its saved point
#      universe, geometry, context markers, objective overlays, and vi/vt
#      solution-region overlays.
#   E. p_mis_base and n_O distributions across vi_C01, vi_C02, and vi_C03.
#
# This script reads saved fit/analysis artifacts only. It does not refit the
# model or recompute an embedding.

options(stringsAsFactors = FALSE, warn = 1)

if (!nzchar(Sys.getenv("XDG_CACHE_HOME", unset = ""))) {
  font_cache <- file.path(tempdir(), "ltee_hypoxia_font_cache")
  dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = font_cache)
}

required_packages <- c(
  "cluster", "dplyr", "ggplot2", "magick", "patchwork", "png", "scales",
  "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) {
  stop("MININGCLONEID_REPO_ROOT is not set. Run with scripts/agentRrunner.sh.")
}
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)

draft_root <- file.path(
  repo_root,
  "agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting"
)
frozen_root <- file.path(draft_root, "source_tables/frozen_inputs/F4")

input_paths <- c(
  burden = file.path(frozen_root, "burden_fit.tsv"),
  terminal = file.path(frozen_root, "terminal_ploidy.tsv"),
  o2 = file.path(frozen_root, "o2_lag.tsv"),
  mode = file.path(frozen_root, "fixed_o2_modes.tsv"),
  parameter_long = file.path(frozen_root, "parameter_values.tsv"),
  pooled_best = file.path(frozen_root, "pooled_coordinates.csv"),
  seed_summary = file.path(frozen_root, "seed_summary.tsv"),
  pooled_raster = file.path(frozen_root, "pooled_embedding_preserved_source.png")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing Figure 4 source artifact(s): ", paste(missing_inputs, collapse = ", "))
}

initial_root <- file.path(draft_root, "initial_subpanels")
refined_root <- file.path(draft_root, "refined_subpanels")
recommended_root <- file.path(draft_root, "final_figures/recommended")
output_dirs <- c(
  file.path(initial_root, paste0("F4", LETTERS[1:5])),
  refined_root,
  recommended_root
)
invisible(vapply(
  output_dirs,
  dir.create,
  logical(1),
  recursive = TRUE,
  showWarnings = FALSE
))

read_tsv <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(label, " lacks required columns: ", paste(missing, collapse = ", "))
  }
  invisible(x)
}

save_plot_pair <- function(plot, stub, width, height, dpi = 300) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  png_path <- paste0(stub, ".png")
  pdf_path <- paste0(stub, ".pdf")
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
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
  invisible(c(png = png_path, pdf = pdf_path))
}

dark <- "#27313A"
blue <- "#0072B2"
vermillion <- "#D55E00"
orange <- "#E69F00"
green <- "#009E73"
grey <- "#5F6B76"
cluster_colors <- c(
  "vi_C01" = "#CC79A7",
  "vi_C02" = orange,
  "vi_C03" = green
)

theme_manuscript <- function(base_size = 8.5) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E4E7E9", linewidth = 0.25),
      plot.title = element_text(size = 9.5, face = "bold", hjust = 0, colour = dark),
      plot.subtitle = element_text(
        size = 7.25, colour = grey, margin = margin(b = 3)
      ),
      plot.title.position = "plot",
      plot.tag = element_text(size = 12, face = "bold", colour = dark),
      plot.tag.position = "topleft",
      axis.title = element_text(size = 8.2, colour = dark),
      axis.text = element_text(size = 7.1, colour = dark),
      strip.text = element_text(size = 7.8, face = "bold", colour = dark),
      strip.background = element_rect(
        fill = "#F1F3F4", colour = "#CFD4D8", linewidth = 0.3
      ),
      legend.title = element_text(size = 7.4, face = "bold"),
      legend.text = element_text(size = 6.9),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(4.8, "mm"),
      plot.margin = margin(5, 6, 5, 6)
    )
}

# -------------------------------------------------------------------------
# Panel A: compact direct in-sample fit adequacy for selected seed25.
# -------------------------------------------------------------------------

burden <- read_tsv(input_paths[["burden"]])
require_columns(
  burden,
  c("harvest", "cohort", "day", "obs_burden", "pred_burden_volume_mm3"),
  basename(input_paths[["burden"]])
)
burden_plot <- burden |>
  transmute(
    harvest = as.character(.data$harvest),
    cohort = factor(as.character(.data$cohort), levels = c("2N", "4N")),
    day = as.numeric(.data$day),
    observed = as.numeric(.data$obs_burden),
    fitted = as.numeric(.data$pred_burden_volume_mm3)
  ) |>
  filter(
    is.finite(.data$observed),
    is.finite(.data$fitted),
    !is.na(.data$cohort)
  )
if (nrow(burden_plot) != 120L || n_distinct(burden_plot$harvest) != 8L) {
  stop(
    "Expected 120 saved burden observations across eight scenarios; found ",
    nrow(burden_plot), " rows and ", n_distinct(burden_plot$harvest), " scenarios."
  )
}
burden_limit <- max(c(burden_plot$observed, burden_plot$fitted), na.rm = TRUE)
burden_breaks <- c(10, 100, 1000)
burden_breaks <- burden_breaks[burden_breaks <= burden_limit * 1.04]
burden_labels <- function(x) {
  ifelse(x >= 1000, paste0(format(x / 1000, trim = TRUE), "k"), format(x, trim = TRUE))
}

p_burden <- ggplot(
  burden_plot,
  aes(x = .data$observed, y = .data$fitted, colour = .data$cohort)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "#737B82",
    linetype = "22",
    linewidth = 0.55
  ) +
  geom_point(size = 1.35, alpha = 0.72) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_colour_manual(values = c("2N" = blue, "4N" = vermillion), guide = "none") +
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1),
    breaks = burden_breaks,
    labels = burden_labels
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1),
    breaks = burden_breaks,
    labels = burden_labels
  ) +
  coord_equal() +
  labs(
    tag = "A",
    title = "Seed25 tumor burden (in-sample)",
    x = expression("Observed burden (mm"^3*")"),
    y = expression("Fitted burden (mm"^3*")")
  ) +
  theme_manuscript()

terminal <- read_tsv(input_paths[["terminal"]])
require_columns(
  terminal,
  c("harvest", "cohort", "source", "endpoint_value", "weight"),
  basename(input_paths[["terminal"]])
)
terminal_plot <- terminal |>
  transmute(
    harvest = as.character(.data$harvest),
    cohort = factor(as.character(.data$cohort), levels = c("2N", "4N")),
    source = factor(
      as.character(.data$source),
      levels = c("Observed", "Predicted"),
      labels = c("Observed", "Fitted")
    ),
    chromosome_number = as.numeric(.data$endpoint_value),
    weight = as.numeric(.data$weight)
  ) |>
  filter(
    !is.na(.data$cohort),
    !is.na(.data$source),
    is.finite(.data$chromosome_number),
    is.finite(.data$weight),
    .data$weight > 0
  )
if (n_distinct(terminal_plot$harvest) != 8L) {
  stop("Terminal chromosome distribution does not contain all eight endpoints.")
}

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- as.numeric(x[keep])
  w <- as.numeric(w[keep])
  if (!length(x)) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cumulative <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cumulative >= p)[1L]], numeric(1))
}

terminal_box <- terminal_plot |>
  group_by(.data$cohort, .data$source) |>
  summarise(
    q1 = weighted_quantile(.data$chromosome_number, .data$weight, 0.25),
    median = weighted_quantile(.data$chromosome_number, .data$weight, 0.50),
    q3 = weighted_quantile(.data$chromosome_number, .data$weight, 0.75),
    raw_min = min(.data$chromosome_number),
    raw_max = max(.data$chromosome_number),
    .groups = "drop"
  ) |>
  mutate(
    iqr = .data$q3 - .data$q1,
    ymin = pmax(.data$raw_min, .data$q1 - 1.5 * .data$iqr),
    ymax = pmin(.data$raw_max, .data$q3 + 1.5 * .data$iqr)
  )

p_terminal <- ggplot(
  terminal_plot,
  aes(
    x = .data$source,
    y = .data$chromosome_number,
    weight = .data$weight,
    fill = .data$source
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    quantiles = NULL,
    colour = "#4B5563",
    linewidth = 0.35,
    alpha = 0.9
  ) +
  geom_boxplot(
    data = terminal_box,
    aes(
      x = .data$source,
      ymin = .data$ymin,
      lower = .data$q1,
      middle = .data$median,
      upper = .data$q3,
      ymax = .data$ymax
    ),
    stat = "identity",
    inherit.aes = FALSE,
    width = 0.15,
    fill = "white",
    colour = dark,
    linewidth = 0.38
  ) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_fill_manual(
    values = c("Observed" = "#4B5563", "Fitted" = "#9DCCE7"),
    guide = "none"
  ) +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  labs(
    title = "Terminal chromosomes (eight endpoints)",
    x = NULL,
    y = "Chromosome no."
  ) +
  theme_manuscript()

panel_a <- p_burden + p_terminal + plot_layout(widths = c(1.05, 0.95))

# -------------------------------------------------------------------------
# Panel B: target versus latent model-implied effective O2.
# -------------------------------------------------------------------------

o2 <- read_tsv(input_paths[["o2"]])
require_columns(
  o2,
  c("harvest", "cohort", "day", "o2_target_pct", "o2_eff_pct"),
  basename(input_paths[["o2"]])
)

scenario_order <- unique(as.character(o2$harvest))
scenario_label <- function(x) {
  x <- sub("^SUM159-", "", as.character(x))
  x <- sub("_harvest$", "", x)
  gsub("-", " \u00b7 ", x, fixed = TRUE)
}

o2_long <- o2 |>
  transmute(
    harvest = factor(as.character(.data$harvest), levels = scenario_order),
    cohort = factor(as.character(.data$cohort), levels = c("2N", "4N")),
    day = as.numeric(.data$day),
    `Target O2` = as.numeric(.data$o2_target_pct),
    `Latent effective O2` = as.numeric(.data$o2_eff_pct)
  ) |>
  pivot_longer(
    cols = c("Target O2", "Latent effective O2"),
    names_to = "series",
    values_to = "o2_pct"
  ) |>
  mutate(
    series = factor(
      .data$series,
      levels = c("Target O2", "Latent effective O2")
    )
  ) |>
  filter(is.finite(.data$day), is.finite(.data$o2_pct))

p_b <- ggplot(
  o2_long,
  aes(
    x = .data$day,
    y = .data$o2_pct,
    group = .data$harvest
  )
) +
  geom_line(
    data = o2_long |> filter(.data$series == "Latent effective O2"),
    colour = blue,
    linetype = "solid",
    linewidth = 0.82,
    alpha = 0.92,
    lineend = "round"
  ) +
  geom_line(
    data = o2_long |> filter(.data$series == "Target O2"),
    colour = orange,
    linetype = "22",
    linewidth = 0.62,
    alpha = 1,
    lineend = "round"
  ) +
  facet_wrap(
    ~ harvest,
    ncol = 4,
    labeller = as_labeller(stats::setNames(
      scenario_label(scenario_order),
      scenario_order
    ))
  ) +
  scale_y_continuous(
    limits = c(0, 5),
    breaks = c(0, 1, 3, 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  labs(
    tag = "B",
    title = expression("Target vs latent effective O"[2]),
    subtitle = "Orange dashed: target; blue: latent model state (not measured)",
    x = "Day",
    y = expression("O"[2]*" (%)")
  ) +
  theme_manuscript() +
  theme(
    legend.position = "none"
  )

# -------------------------------------------------------------------------
# Panel C: O2 = 0/1/5 univariate attractor-mode separator triptych.
# -------------------------------------------------------------------------

mode <- read_tsv(input_paths[["mode"]])
parameter_long <- read_tsv(input_paths[["parameter_long"]])
require_columns(
  mode,
  c("seed_id", "O2_pct", "mode_label", "status"),
  basename(input_paths[["mode"]])
)
require_columns(
  parameter_long,
  c("seed_id", "parameter", "value"),
  basename(input_paths[["parameter_long"]])
)

mode_selected <- mode |>
  transmute(
    seed_id = as.character(.data$seed_id),
    O2_pct = as.numeric(.data$O2_pct),
    mode_label = as.character(.data$mode_label),
    status = as.character(.data$status)
  ) |>
  filter(
    .data$O2_pct %in% c(0, 1, 5),
    .data$status == "ok",
    .data$mode_label %in% c("mode1", "mode2")
  )

separator_input <- parameter_long |>
  transmute(
    seed_id = as.character(.data$seed_id),
    parameter = as.character(.data$parameter),
    value = as.numeric(.data$value)
  ) |>
  inner_join(mode_selected, by = "seed_id", relationship = "many-to-many") |>
  filter(is.finite(.data$value))

rank_auc <- function(value, label) {
  positive <- label == "mode1"
  n_positive <- sum(positive)
  n_negative <- sum(!positive)
  if (!n_positive || !n_negative) return(NA_real_)
  ranks <- rank(value, ties.method = "average")
  (
    sum(ranks[positive]) - n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)
}

separator_table <- separator_input |>
  group_by(.data$O2_pct, .data$parameter) |>
  summarise(
    signed_auc = rank_auc(.data$value, .data$mode_label),
    n = n(),
    .groups = "drop"
  ) |>
  mutate(
    auc = pmax(.data$signed_auc, 1 - .data$signed_auc),
    direction = ifelse(
      .data$signed_auc >= 0.5,
      "higher-ploidy mode",
      "lower-ploidy mode"
    )
  ) |>
  group_by(.data$O2_pct) |>
  arrange(desc(.data$auc), .data$parameter, .by_group = TRUE) |>
  slice_head(n = 4) |>
  ungroup()

if (!all(c("mu_hp", "p_mis_base") %in% separator_table$parameter)) {
  stop("Expected reproduced low/high-O2 separator features are absent.")
}
low_top <- separator_table |>
  filter(.data$O2_pct == 0) |>
  slice_max(.data$auc, n = 1, with_ties = FALSE)
high_top <- separator_table |>
  filter(.data$O2_pct == 5) |>
  slice_max(.data$auc, n = 1, with_ties = FALSE)
if (
  low_top$parameter != "mu_hp" ||
    abs(low_top$auc - 0.849146) > 5e-4 ||
    high_top$parameter != "p_mis_base" ||
    abs(high_top$auc - 0.903459) > 5e-4
) {
  stop("Reproduced low/high-O2 feature ranking or AUC differs from the approved evidence.")
}

separator_table <- separator_table |>
  arrange(.data$O2_pct, .data$auc) |>
  mutate(
    oxygen_facet = factor(
      .data$O2_pct,
      levels = c(0, 1, 5),
      labels = c("O2 = 0%", "O2 = 1%", "O2 = 5%")
    ),
    parameter_key = factor(
      paste(.data$parameter, .data$O2_pct, sep = "___"),
      levels = paste(.data$parameter, .data$O2_pct, sep = "___")
    )
  )

separator_display <- separator_table |>
  group_by(.data$O2_pct) |>
  slice_max(.data$auc, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    direction_short = ifelse(
      .data$direction == "higher-ploidy mode",
      "larger:\nhigher-ploidy",
      "larger:\nlower-ploidy"
    )
  )

p_c <- ggplot(
  separator_display,
  aes(x = 0.5, y = 0.5)
) +
  geom_text(
    aes(label = .data$parameter),
    x = 0.5,
    y = 0.72,
    size = 2.9,
    fontface = "italic",
    colour = dark
  ) +
  geom_text(
    aes(label = sprintf("AUC = %.3f", .data$auc)),
    x = 0.5,
    y = 0.43,
    size = 2.55,
    fontface = "bold",
    colour = dark
  ) +
  geom_text(
    aes(
      label = .data$direction_short,
      colour = .data$direction
    ),
    x = 0.5,
    y = 0.14,
    size = 2.5,
    lineheight = 0.88,
    show.legend = FALSE
  ) +
  facet_wrap(~ oxygen_facet, nrow = 1) +
  scale_colour_manual(
    values = c("higher-ploidy mode" = blue, "lower-ploidy mode" = vermillion)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    tag = "C",
    title = expression("Univariate fixed-O"[2]*" separators"),
    subtitle = "Discrimination, not causal influence",
    x = NULL,
    y = NULL
  ) +
  theme_manuscript() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#D5DADD", fill = NA, linewidth = 0.3),
    strip.background = element_rect(
      fill = "#F1F3F4", colour = "#CFD4D8", linewidth = 0.3
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

# -------------------------------------------------------------------------
# Panel D: exact preservation of the approved pooled t-SNE render.
# -------------------------------------------------------------------------

blob_size <- as.numeric(file.info(input_paths[["pooled_raster"]])$size)
if (length(blob_size) != 1L || !is.finite(blob_size) || blob_size <= 0) {
  stop("Could not resolve the frozen preserved pooled-embedding raster.")
}
blob_connection <- file(input_paths[["pooled_raster"]], open = "rb")
pooled_raw <- readBin(
  blob_connection,
  what = "raw",
  n = as.integer(blob_size)
)
close(blob_connection)
if (length(pooled_raw) != as.integer(blob_size)) {
  stop("Preserved pooled-embedding blob was read incompletely.")
}
pooled_image <- magick::image_read(pooled_raw)
pooled_info <- magick::image_info(pooled_image)
if (
  pooled_info$width != 1628L ||
    pooled_info$height != 1430L ||
    length(pooled_raw) != 884007L
) {
  stop("Preserved pooled-embedding raster dimensions or byte count changed.")
}

pooled_preserved_png <- file.path(
  initial_root,
  "F4D/figure4D_pooled_embedding_preserved_source.png"
)
pooled_output_connection <- file(pooled_preserved_png, open = "wb")
writeBin(pooled_raw, pooled_output_connection)
close(pooled_output_connection)
if (
  unname(tools::md5sum(pooled_preserved_png)) !=
    "14cecff29ab4884823b84d83f0379119"
) {
  stop("Preserved pooled-embedding PNG checksum differs from the approved source.")
}

pooled_plot_crop <- magick::image_crop(
  pooled_image,
  geometry = "1260x1260+10+70",
  repage = TRUE
)
pooled_legend_crop <- magick::image_crop(
  pooled_image,
  geometry = "328x1120+1300+145",
  repage = TRUE
)
pooled_legend_crop <- magick::image_extent(
  pooled_legend_crop,
  geometry = "328x1260",
  gravity = "center",
  color = "white"
)
pooled_layout_image <- magick::image_append(
  c(pooled_plot_crop, pooled_legend_crop),
  stack = FALSE
)
pooled_layout_info <- magick::image_info(pooled_layout_image)
if (
  pooled_layout_info$width != 1588L ||
    pooled_layout_info$height != 1260L
) {
  stop("Unexpected dimensions after lossless Figure 4D whitespace/legend reflow.")
}
pooled_layout_png <- file.path(
  initial_root,
  "F4D/figure4D_pooled_embedding_layout_source.png"
)
magick::image_write(
  pooled_layout_image,
  path = pooled_layout_png,
  format = "png"
)
pooled_layout_array <- png::readPNG(pooled_layout_png)
p_d_image <- ggplot() +
  annotation_raster(
    pooled_layout_array,
    xmin = 0,
    xmax = pooled_layout_info$width,
    ymin = 0,
    ymax = pooled_layout_info$height
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(0, pooled_layout_info$width),
    ylim = c(0, pooled_layout_info$height),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
p_d_header <- ggplot() +
  annotate(
    "text", x = 0, y = 0.72, label = "D", hjust = 0,
    size = 4.2, fontface = "bold", colour = dark
  ) +
  annotate(
    "text", x = 0.075, y = 0.72, label = "Pooled parameter embedding",
    hjust = 0, size = 3.2, fontface = "bold", colour = dark
  ) +
  annotate(
    "text", x = 0.075, y = 0.20,
    label = "Saved in-vivo/in-vitro coordinates, point universe, and overlays",
    hjust = 0, size = 2.3, colour = grey
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(0, 4, 0, 4))
p_d <- p_d_header / p_d_image +
  plot_layout(heights = c(0.075, 0.925))

# -------------------------------------------------------------------------
# Panel E: exact formal-region membership on saved coordinates, then
#          regenerated fitted-parameter distributions.
# -------------------------------------------------------------------------

pooled_best <- utils::read.csv(
  input_paths[["pooled_best"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
require_columns(
  pooled_best,
  c("tSNE1", "tSNE2", "dataset", "point_type", "seed", "objective"),
  basename(input_paths[["pooled_best"]])
)
in_vivo_best <- pooled_best |>
  filter(.data$dataset == "invivo", .data$point_type == "best") |>
  transmute(
    seed_number = as.integer(.data$seed),
    tSNE1 = as.numeric(.data$tSNE1),
    tSNE2 = as.numeric(.data$tSNE2),
    objective = as.numeric(.data$objective)
  ) |>
  arrange(.data$seed_number)
if (
  nrow(in_vivo_best) != 500L ||
    anyDuplicated(in_vivo_best$seed_number) ||
    any(!is.finite(in_vivo_best$tSNE1)) ||
    any(!is.finite(in_vivo_best$tSNE2))
) {
  stop("Expected 500 unique finite saved in-vivo best-fit t-SNE coordinates.")
}

relabel_clusters_by_median_x <- function(cluster_number, coordinates) {
  centers <- aggregate(
    cbind(tSNE1, tSNE2) ~ cluster,
    data = data.frame(cluster = cluster_number, coordinates),
    FUN = stats::median
  )
  centers <- centers[order(centers$tSNE1, centers$tSNE2), , drop = FALSE]
  relabel <- stats::setNames(seq_len(nrow(centers)), centers$cluster)
  as.integer(relabel[as.character(cluster_number)])
}

reproduce_formal_clusters <- function(coordinates, seed = 123L) {
  matrix_xy <- as.matrix(coordinates[, c("tSNE1", "tSNE2")])
  storage.mode(matrix_xy) <- "double"
  distance_xy <- stats::dist(matrix_xy)
  k_values <- 2:8
  scores <- rep(NA_real_, length(k_values))
  centers_by_k <- vector("list", length(k_values))
  for (i in seq_along(k_values)) {
    k <- k_values[[i]]
    set.seed(seed + k)
    fit <- stats::kmeans(
      matrix_xy,
      centers = k,
      nstart = 10L,
      iter.max = 100L
    )
    scores[[i]] <- mean(cluster::silhouette(fit$cluster, distance_xy)[, "sil_width"])
    centers_by_k[[i]] <- fit$centers
  }
  selected_i <- which.max(scores)
  selected_k <- k_values[[selected_i]]
  final_fit <- stats::kmeans(
    matrix_xy,
    centers = centers_by_k[[selected_i]],
    iter.max = 100L,
    algorithm = "Lloyd"
  )
  cluster_number <- relabel_clusters_by_median_x(
    final_fit$cluster,
    coordinates[, c("tSNE1", "tSNE2")]
  )
  list(
    selected_k = selected_k,
    scores = stats::setNames(scores, k_values),
    region = sprintf("vi_C%02d", cluster_number)
  )
}

formal_cluster <- reproduce_formal_clusters(in_vivo_best)
in_vivo_best$region <- formal_cluster$region
region_counts <- table(in_vivo_best$region)
expected_region_counts <- c("vi_C01" = 99L, "vi_C02" = 385L, "vi_C03" = 16L)
if (
  formal_cluster$selected_k != 3L ||
    !identical(
      as.integer(region_counts[names(expected_region_counts)]),
      as.integer(expected_region_counts)
    )
) {
  stop(
    "Saved-coordinate clustering did not reproduce formal regions ",
    "vi_C01/02/03 = 99/385/16."
  )
}

seed_summary <- read_tsv(input_paths[["seed_summary"]])
require_columns(
  seed_summary,
  c("seed", "value__p_mis_base", "value__n_O"),
  basename(input_paths[["seed_summary"]])
)
parameter_by_region <- in_vivo_best |>
  mutate(seed = paste0("seed", .data$seed_number)) |>
  inner_join(
    seed_summary |>
      transmute(
        seed = as.character(.data$seed),
        p_mis_base = as.numeric(.data$value__p_mis_base),
        n_O = as.numeric(.data$value__n_O)
      ),
    by = "seed"
  )
if (
  nrow(parameter_by_region) != 500L ||
    any(!is.finite(parameter_by_region$p_mis_base)) ||
    any(parameter_by_region$p_mis_base <= 0) ||
    any(!is.finite(parameter_by_region$n_O))
) {
  stop("Could not match all formal in-vivo regions to finite fitted parameters.")
}

region_labels <- paste0(
  names(expected_region_counts),
  "\n(n = ",
  expected_region_counts,
  ")"
)
names(region_labels) <- names(expected_region_counts)
parameter_distribution <- bind_rows(
  parameter_by_region |>
    transmute(
      region = .data$region,
      parameter = "log10(p_mis_base)",
      value = log10(.data$p_mis_base)
    ),
  parameter_by_region |>
    transmute(
      region = .data$region,
      parameter = "n_O",
      value = .data$n_O
    )
) |>
  mutate(
    region = factor(
      .data$region,
      levels = names(expected_region_counts),
      labels = unname(region_labels)
    ),
    parameter = factor(
      .data$parameter,
      levels = c("log10(p_mis_base)", "n_O"),
      labels = c(
        "log10(p_mis_base)\nbaseline missegregation",
        "n_O\nO2-response steepness"
      )
    )
  )

distribution_colors <- stats::setNames(
  unname(cluster_colors),
  unname(region_labels)
)
p_e <- ggplot(
  parameter_distribution,
  aes(y = .data$region, x = .data$value, fill = .data$region)
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    quantiles = NULL,
    orientation = "y",
    colour = "#4B5563",
    linewidth = 0.35,
    alpha = 0.78
  ) +
  geom_boxplot(
    width = 0.13,
    orientation = "y",
    outlier.shape = NA,
    fill = "white",
    colour = dark,
    linewidth = 0.38
  ) +
  facet_wrap(~ parameter, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = distribution_colors, guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  labs(
    tag = "E",
    title = "Fitted parameters across\nin-vivo regions",
    subtitle = "Solution families, not tumor subtypes",
    x = "Fitted parameter value",
    y = NULL
  ) +
  theme_manuscript() +
  theme(
    axis.text.y = element_text(size = 7.0, lineheight = 0.95),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

# -------------------------------------------------------------------------
# Save panel drafts and assemble the recommended Figure 4.
# -------------------------------------------------------------------------

initial_outputs <- c(
  save_plot_pair(
    panel_a,
    file.path(initial_root, "F4A/figure4A_fit_adequacy"),
    width = 7.1,
    height = 2.10
  ),
  save_plot_pair(
    p_b,
    file.path(initial_root, "F4B/figure4B_latent_effective_o2"),
    width = 4.15,
    height = 2.05
  ),
  save_plot_pair(
    p_c,
    file.path(initial_root, "F4C/figure4C_univariate_separators"),
    width = 2.95,
    height = 2.05
  ),
  save_plot_pair(
    p_d,
    file.path(initial_root, "F4D/figure4D_pooled_embedding"),
    width = 4.85,
    height = 4.40
  ),
  save_plot_pair(
    p_e,
    file.path(initial_root, "F4E/figure4E_region_parameters"),
    width = 2.25,
    height = 4.40
  )
)

panel_png_paths <- c(
  A = file.path(initial_root, "F4A/figure4A_fit_adequacy.png"),
  B = file.path(initial_root, "F4B/figure4B_latent_effective_o2.png"),
  C = file.path(initial_root, "F4C/figure4C_univariate_separators.png"),
  D = file.path(initial_root, "F4D/figure4D_pooled_embedding.png"),
  E = file.path(initial_root, "F4E/figure4E_region_parameters.png")
)
if (any(!file.exists(panel_png_paths))) {
  stop("One or more rendered Figure 4 panel PNGs are missing.")
}

panel_raster <- function(path, x, y, width, height) {
  grid::rasterGrob(
    png::readPNG(path),
    interpolate = TRUE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    vp = grid::viewport(
      x = grid::unit(x, "in"),
      y = grid::unit(y, "in"),
      width = grid::unit(width, "in"),
      height = grid::unit(height, "in"),
      just = c("left", "bottom")
    )
  )
}

footer_text <- paste0(
  "A: seed25 is the best weighted-MAP fit among 500 starts; necrosis excluded ",
  "(no exported predictions).\n",
  "C: univariate discrimination, not causal influence.\n",
  "D: approved pooled embedding preserved; whitespace/legend reflow only—",
  "no re-embedding or refit."
)
figure4 <- grid::grobTree(
  panel_raster(panel_png_paths[["A"]], 0, 6.90, 7.10, 2.10),
  panel_raster(panel_png_paths[["B"]], 0, 4.85, 4.15, 2.05),
  panel_raster(panel_png_paths[["C"]], 4.15, 4.85, 2.95, 2.05),
  panel_raster(panel_png_paths[["D"]], 0, 0.45, 4.85, 4.40),
  panel_raster(panel_png_paths[["E"]], 4.85, 0.45, 2.25, 4.40),
  grid::textGrob(
    footer_text,
    x = grid::unit(0.15, "in"),
    y = grid::unit(0.225, "in"),
    hjust = 0,
    vjust = 0.5,
    gp = grid::gpar(
      fontsize = 7,
      lineheight = 0.95,
      col = grey,
      fontfamily = "sans"
    )
  )
)

refined_outputs <- save_plot_pair(
  figure4,
  file.path(refined_root, "figure4_refined"),
  width = 7.1,
  height = 9.0,
  dpi = 300
)
recommended_outputs <- save_plot_pair(
  figure4,
  file.path(recommended_root, "figure4"),
  width = 7.1,
  height = 9.0,
  dpi = 300
)

cat("Figure 4 generation complete.\n")
cat("Burden observations:", nrow(burden_plot), "across", n_distinct(burden_plot$harvest), "scenarios\n")
cat("Terminal comparison rows:", nrow(terminal_plot), "\n")
cat(
  "Low-O2 top separator:",
  low_top$parameter,
  sprintf("(AUC %.6f)", low_top$auc),
  "\n"
)
cat(
  "High-O2 top separator:",
  high_top$parameter,
  sprintf("(AUC %.6f)", high_top$auc),
  "\n"
)
cat(
  "Formal region counts:",
  paste(names(region_counts), as.integer(region_counts), sep = "=", collapse = ", "),
  "\n"
)
cat("Preserved pooled PNG:", pooled_preserved_png, "\n")
cat("Panel outputs:\n", paste(initial_outputs, collapse = "\n"), "\n")
cat("Recommended outputs:\n", paste(recommended_outputs, collapse = "\n"), "\n")
