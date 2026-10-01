#!/usr/bin/env Rscript

# Visualize the O2-window comparisons used to classify Figure 4B rows.
# The resampling unit is a complete optimizer-derived fitted endpoint with its
# full 201-point O2 curve. These endpoints are not biological replicates or
# posterior samples.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_path <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
  } else {
    normalizePath("draw_Supp_Figure4_3.R", mustWork = FALSE)
  }
})
data_dir <- normalizePath(
  Sys.getenv("ANALYSIS_DATA_DIR", unset = file.path(dirname(script_path), "data")),
  mustWork = TRUE
)
output_dir <- normalizePath(
  Sys.getenv("DELIVERABLE_OUTPUT_DIR", unset = file.path(dirname(script_path), "figures")),
  mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  scores = file.path(
    data_dir, "continuous_ploidy_o2_window_absrho_scores.tsv"
  ),
  tests = file.path(
    data_dir, "continuous_ploidy_o2_window_pairwise_tests.tsv"
  ),
  classification = file.path(
    data_dir, "continuous_ploidy_o2_window_classification.tsv"
  ),
  ranking = file.path(data_dir, "continuous_ploidy_parameter_ranking.tsv")
)
missing <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing)) {
  stop("Missing Supplementary Figure 4-3 input(s): ", paste(missing, collapse = ", "))
}

scores <- fread(paths$scores)
tests <- fread(paths$tests)
classification <- fread(paths$classification)
ranking <- fread(paths$ranking)
if (nrow(scores) != 18L * 2L ||
    nrow(tests) != 18L ||
    nrow(classification) != 18L ||
    nrow(ranking) != 18L ||
    uniqueN(scores$parameter) != 18L ||
    uniqueN(tests$parameter) != 18L ||
    uniqueN(classification$parameter) != 18L) {
  stop("Supplementary Figure 4-3 requires 18 parameters and three windows/contrasts.")
}
if (any(scores$observed_mean_abs_rho < -1e-12 |
        scores$observed_mean_abs_rho > 1 + 1e-12) ||
    any(scores$bootstrap_ci_lower < -1e-12 |
        scores$bootstrap_ci_upper > 1 + 1e-12) ||
    any(tests$bh_adjusted_p_value < 0 |
        tests$bh_adjusted_p_value > 1)) {
  stop("Supplementary Figure 4-3 statistics lie outside their valid ranges.")
}

group_levels <- c("High O2", "Low O2", "O2-independent")
o2_association_group_palette <- c(
  "High O2" = "#B2182B",
  "Low O2" = "#2166AC",
  "O2-independent" = "#8A8A8A"
)
o2_association_group_text_palette <- setNames(
  rep("#FFFFFF", length(o2_association_group_palette)),
  names(o2_association_group_palette)
)
window_levels <- c("High O2", "Low O2")
contrast_levels <- "Low - High"
parameter_levels <- ranking[order(display_order), parameter]

prepare_plot_table <- function(table) {
  table[, parameter_factor := factor(
    parameter,
    levels = rev(parameter_levels)
  )]
  table[, o2_association_group := factor(
    o2_association_group,
    levels = group_levels
  )]
  table
}
scores <- prepare_plot_table(scores)
tests <- prepare_plot_table(tests)
scores[, o2_window := factor(o2_window, levels = window_levels)]
tests[, contrast := factor(contrast, levels = contrast_levels)]
tests[, q_label := fifelse(
  bh_adjusted_p_value < 0.001,
  "q<0.001",
  sprintf("q=%.3f", bh_adjusted_p_value)
)]

window_palette <- c(
  "Low O2" = "#2166AC",
  "High O2" = "#B2182B"
)
window_shapes <- c("Low O2" = 21, "High O2" = 24)
score_position <- position_dodge(width = 0.62)

theme_supp4_3 <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(face = "bold", color = "#202428"),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9.2, color = "#4B5259"),
      plot.tag = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 10.5),
      axis.text = element_text(size = 8.3, color = "#252A2F"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#ECEDEF", linewidth = 0.28),
      panel.grid.major.x = element_line(color = "#E0E2E5", linewidth = 0.28),
      panel.border = element_rect(color = "#858A90", linewidth = 0.35),
      strip.background = element_rect(fill = "#F1F2F4", color = "#A7ABB0"),
      strip.text = element_text(size = 8.3, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 8.6),
      plot.margin = margin(7, 8, 7, 8)
    )
}

panel_a <- ggplot(
  scores,
  aes(x = observed_mean_abs_rho, y = parameter_factor, color = o2_window)
) +
  geom_errorbar(
    aes(xmin = bootstrap_ci_lower, xmax = bootstrap_ci_upper),
    orientation = "y", width = 0.18, linewidth = 0.65,
    position = score_position
  ) +
  geom_point(
    aes(fill = o2_window, shape = o2_window),
    size = 2.8, stroke = 0.55, color = "#202428",
    position = score_position
  ) +
  facet_grid(
    rows = vars(o2_association_group),
    scales = "free_y", space = "free_y", switch = "y", drop = TRUE
  ) +
  scale_color_manual(values = window_palette, drop = FALSE) +
  scale_fill_manual(values = window_palette, drop = FALSE) +
  scale_shape_manual(values = window_shapes, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    tag = "A",
    title = "Window-averaged absolute parameter-ploidy association",
    subtitle = "Points are observed normalized AUC; bars are 95% endpoint-bootstrap intervals",
    x = expression(paste("Window mean |Spearman ", rho, "|")),
    y = NULL
  ) +
  theme_supp4_3() +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "top"
  )

contrast_limit <- max(abs(c(
  tests$bootstrap_ci_lower,
  tests$bootstrap_ci_upper,
  tests$observed_delta_mean_abs_rho
)), na.rm = TRUE)
contrast_limit <- max(0.05, ceiling(contrast_limit * 20) / 20)
label_offset <- contrast_limit * 0.04
tests[, q_label_x := pmin(
  contrast_limit * 1.27,
  pmax(bootstrap_ci_upper, observed_delta_mean_abs_rho) + label_offset
)]

panel_b <- ggplot(
  tests,
  aes(x = observed_delta_mean_abs_rho, y = parameter_factor)
) +
  geom_vline(xintercept = 0, color = "#4C5157", linewidth = 0.5) +
  geom_errorbar(
    aes(xmin = bootstrap_ci_lower, xmax = bootstrap_ci_upper),
    orientation = "y", width = 0.18, linewidth = 0.65,
    color = "#5B6066"
  ) +
  geom_point(
    aes(fill = significant_bh_0p05),
    shape = 21, size = 2.6, stroke = 0.55, color = "#202428"
  ) +
  geom_text(
    aes(x = q_label_x, label = q_label),
    hjust = 0, size = 2.35, fontface = "bold", color = "#30353A"
  ) +
  facet_grid(
    rows = vars(o2_association_group),
    cols = vars(contrast),
    scales = "free_y", space = "free_y", drop = TRUE
  ) +
  scale_fill_manual(
    values = c("FALSE" = "white", "TRUE" = "#6A51A3"),
    labels = c("FALSE" = "BH q >= 0.05", "TRUE" = "BH q < 0.05"),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(-contrast_limit, contrast_limit * 1.35),
    breaks = pretty(c(-contrast_limit, contrast_limit), n = 5),
    expand = expansion(mult = c(0.02, 0.01))
  ) +
  labs(
    tag = "B",
    title = "Low-High O2-window contrast",
    subtitle = "Positive values favor Low O2; labels report BH-adjusted q values",
    x = expression(paste(Delta, " window mean |Spearman ", rho, "|")),
    y = NULL,
    fill = NULL
  ) +
  theme_supp4_3() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    strip.background.y = element_blank(),
    legend.position = "top"
  )

figure <- panel_a + panel_b +
  plot_layout(widths = c(0.55, 0.45), guides = "collect") +
  plot_annotation(
    title = "O2-window association tests underlying Figure 4B parameter groups",
    subtitle = paste(
      "Low [0,1.5] and High [3,5] scores use normalized",
      paste0(
        "trapezoid AUC of |rho|; ",
        comma(unique(scores$bootstrap_reps)),
        " complete-endpoint bootstrap replicates; (1.5,3) excluded"
      )
    ),
    caption = paste0(
      "All 18 Low-High comparisons are adjusted together by Benjamini-Hochberg. ",
      "Detected contrasts assign parameters to the higher-scoring window.\n",
      "O2-independent means no detected Low-High difference, not absence of association. ",
      "Optimizer-derived endpoints are not biological replicates or posterior samples."
    ),
    theme = theme(
      plot.title = element_text(
        family = "Arial", size = 16, face = "bold", color = "#15191D"
      ),
      plot.subtitle = element_text(
        family = "Arial", size = 10.2, face = "bold", color = "#3E454C"
      ),
      plot.caption = element_text(
        family = "Arial", size = 8.2, face = "bold", color = "#4B5259",
        hjust = 0
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

set_grob_text_color <- function(grob, color) {
  if (inherits(grob, "text")) {
    grob$gp$col <- color
  }
  if (!is.null(grob$children)) {
    for (index in seq_along(grob$children)) {
      grob$children[[index]] <- set_grob_text_color(
        grob$children[[index]], color
      )
    }
  }
  if (!is.null(grob$grobs)) {
    for (index in seq_along(grob$grobs)) {
      grob$grobs[[index]] <- set_grob_text_color(grob$grobs[[index]], color)
    }
  }
  grob
}

color_panel_a_group_strips <- function(figure, group_palette, text_palette) {
  figure_grob <- patchwork::patchworkGrob(figure)
  strip_index <- grep(
    "^strip-l.*-1$", figure_grob$layout$name
  )
  if (length(strip_index) != 1L) {
    stop(
      "Expected exactly one Panel A left-strip collection; found ",
      length(strip_index), "."
    )
  }
  strip_collection <- figure_grob$grobs[[strip_index]]
  present_groups <- group_levels[
    group_levels %in% as.character(unique(scores$o2_association_group))
  ]
  if (length(strip_collection$grobs) != length(present_groups)) {
    stop(
      "Panel A strip count does not match the nonempty O2 association groups."
    )
  }
  for (index in seq_along(present_groups)) {
    group <- present_groups[[index]]
    strip <- strip_collection$grobs[[index]]
    if (length(strip$grobs) != 1L || is.null(strip$grobs[[1L]]$children)) {
      stop("Unexpected Panel A facet-strip grob structure.")
    }
    strip_tree <- strip$grobs[[1L]]
    background_index <- grep(
      "^strip\\.background", strip_tree$childrenOrder
    )
    text_index <- grep("^strip\\.text", strip_tree$childrenOrder)
    if (length(background_index) != 1L || length(text_index) != 1L) {
      stop("Panel A facet strip lacks one background and one text grob.")
    }
    strip_tree$children[[background_index]]$gp$fill <-
      unname(group_palette[[group]])
    strip_tree$children[[background_index]]$gp$col <- "#5B6066"
    strip_tree$children[[text_index]] <- set_grob_text_color(
      strip_tree$children[[text_index]], unname(text_palette[[group]])
    )
    strip$grobs[[1L]] <- strip_tree
    strip_collection$grobs[[index]] <- strip
  }
  figure_grob$grobs[[strip_index]] <- strip_collection
  figure_grob
}

stem <- file.path(output_dir, "Supp_Figure4_3")
figure_width <- 17
figure_height <- 12
build_figure_grob <- function() {
  metric_device <- file.path(
    output_dir, ".Supp_Figure4_3_grob_metrics.png"
  )
  grDevices::png(
    metric_device,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = 300,
    type = "cairo",
    bg = "white"
  )
  on.exit({
    grDevices::dev.off()
    unlink(metric_device)
  }, add = TRUE)
  color_panel_a_group_strips(
    figure,
    group_palette = o2_association_group_palette,
    text_palette = o2_association_group_text_palette
  )
}
figure_grob <- build_figure_grob()
ggsave(
  paste0(stem, ".png"), figure_grob,
  width = figure_width, height = figure_height,
  units = "in", dpi = 300, bg = "white"
)
ggsave(
  paste0(stem, ".pdf"), figure_grob,
  width = figure_width, height = figure_height,
  units = "in", device = cairo_pdf, bg = "white"
)
ggsave(
  paste0(stem, ".svg"), figure_grob,
  width = figure_width, height = figure_height,
  units = "in", device = svglite::svglite, bg = "white"
)

validation <- data.table(
  metric = c(
    "n_parameters", "n_windows", "n_pairwise_tests",
    "bootstrap_reps", "bootstrap_seed", "multiple_testing",
    "window_score", "figure_width_in", "figure_height_in",
    "left_annotation_field", "left_annotation_palette",
    "left_annotation_colored",
    "png_rendered", "pdf_rendered", "svg_rendered",
    "endpoint_rows_are_biological_replicates",
    "endpoint_rows_are_posterior_samples"
  ),
  value = c(
    uniqueN(scores$parameter), uniqueN(scores$o2_window), nrow(tests),
    paste(unique(scores$bootstrap_reps), collapse = ","),
    paste(unique(scores$bootstrap_seed), collapse = ","),
    "Benjamini-Hochberg across 18 Low-High contrasts",
    "normalized trapezoid AUC of absolute Spearman rho",
    figure_width, figure_height,
    "o2_association_group",
    paste(
      names(o2_association_group_palette),
      unname(o2_association_group_palette),
      sep = "=", collapse = ";"
    ),
    "TRUE",
    file.exists(paste0(stem, ".png")),
    file.exists(paste0(stem, ".pdf")),
    file.exists(paste0(stem, ".svg")),
    "FALSE", "FALSE"
  )
)
fwrite(
  validation,
  file.path(data_dir, "supp_figure4_3_validation.tsv"),
  sep = "\t"
)
provenance <- data.table(
  source = c(names(paths), "script"),
  path = c(unname(unlist(paths)), script_path),
  md5 = unname(tools::md5sum(c(unname(unlist(paths)), script_path))),
  role = c(rep("input", length(paths)), "script")
)
fwrite(
  provenance,
  file.path(data_dir, "supp_figure4_3_source_provenance.tsv"),
  sep = "\t"
)

message("Supplementary Figure 4-3 written to: ", paste0(stem, ".png"))
