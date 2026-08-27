#!/usr/bin/env Rscript

if (identical(Sys.getenv("SUPP_FIGURE5_1_DRAW_WORKER"), "1")) {

# Standalone Supplementary Figure 5-1 generator.
#
# Integrated display of endpoint-level directional composition and horizontal
# log2(in vivo / in vitro) distributions for C01/C02/C03. The categorical
# component uses the declared 0.8/1.2 outer-inclusive rule; the distributional
# component prevents thresholding or medians from obscuring solution spread.
#
# All plotting helpers are defined here. The script reads only the regenerated
# iteration2 tables for one prespecified primary pair per C family and performs no
# refitting.

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("ggplot2", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
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
input_root <- normalizePath(
  Sys.getenv(
    "SUPP_FIGURE5_1_DATA_DIR",
    unset = file.path(
      workspace_root, "data", "Figures", "Supp_Figure5_1"
    )
  ),
  mustWork = TRUE
)
figure_root <- normalizePath(
  Sys.getenv(
    "SUPP_FIGURE5_1_FIGURE_DIR",
    unset = file.path(workspace_root, "Figures")
  ),
  mustWork = FALSE
)
panel_root <- normalizePath(
  Sys.getenv(
    "SUPP_FIGURE5_1_PANEL_DIR",
    unset = file.path(figure_root, "supp_figure5-1_subpanels")
  ),
  mustWork = FALSE
)
revised_root <- figure_root
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_root, recursive = TRUE, showWarnings = FALSE)
dir.create(revised_root, recursive = TRUE, showWarnings = FALSE)

obsolete_panel_stubs <- file.path(
  panel_root,
  c(
    "supp_figure5-1a_pair_parameter_ratio_map",
    "supp_figure5-1b_between_pair_directional_consensus",
    "supp_figure5-1_continuous_magnitude"
  )
)
obsolete_panel_files <- as.vector(outer(
  obsolete_panel_stubs,
  c(".png", ".pdf"),
  paste0
))
invisible(file.remove(obsolete_panel_files[file.exists(obsolete_panel_files)]))

within_path <- file.path(input_root, "within_pair_parameter_stability.tsv")
between_path <- file.path(input_root, "between_pair_parameter_stability.tsv")
config_path <- file.path(input_root, "analysis_config.tsv")
pair_path <- file.path(input_root, "selected_primary_family_pairs.tsv")
master_path <- file.path(input_root, "soft_coupling_master_long.tsv")
input_paths <- c(
  within_path, between_path, config_path, pair_path, master_path
)
if (any(!file.exists(input_paths))) {
  stop("Missing Supplementary Figure 5-1 input(s):\n", paste(input_paths[!file.exists(input_paths)], collapse = "\n"))
}

read_tsv <- function(path) {
  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

parameter_levels <- function() {
  c(
    "O2_crit", "n_O", "alpha_o2", "gamma_growth", "lam_max", "mu_hp",
    "gamma_mu", "p_mis_base", "p_misseg", "k_o_mis", "p_wgd",
    "buffer_smax", "buffer_beta", "buffer_n_exp"
  )
}

theme_si <- function(base_size = 8) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#E3E7EB", linewidth = 0.25),
      axis.text = element_text(color = "#3B4650"),
      axis.title = element_text(color = "#26313A"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", color = "#1F2933"),
      plot.subtitle = element_text(color = "#56616B"),
      plot.caption = element_text(color = "#6B747C", hjust = 0),
      plot.margin = margin(5, 7, 5, 5)
    )
}

save_plot_pair <- function(plot, stub, width, height, dpi = 300) {
  ggsave(
    paste0(stub, ".pdf"),
    plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    paste0(stub, ".png"),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
}

within <- read_tsv(within_path)
between <- read_tsv(between_path)
config <- read_tsv(config_path)
pairs <- read_tsv(pair_path)
master <- read_tsv(master_path)

required_within <- c(
  "family", "pair_id", "parameter", "n_valid",
  "n_ClassA", "n_ClassB", "n_ClassC",
  "prop_ClassA", "prop_ClassB", "prop_ClassC",
  "log2_ratio_median"
)
required_between <- c(
  "parameter", "n_pairs"
)
required_master <- c(
  "family", "pair_id", "parameter", "seed_number",
  "log2_ratio_vivo_to_vitro", "ratio_class"
)
if (!all(required_within %in% names(within)) ||
    !all(required_between %in% names(between)) ||
    !all(required_master %in% names(master))) {
  stop("Supplementary Figure 5-1 tables lack required fields")
}
if (nrow(within) != 42L || nrow(between) != 14L) {
  stop("Expected 42 selected-family parameter rows and 14 between-pair rows")
}

config_value <- function(key) {
  hit <- config$value[config$key == key]
  if (!length(hit)) NA_character_ else as.character(hit[[1L]])
}
lower_bound <- as.numeric(config_value("class_lower_bound"))
upper_bound <- as.numeric(config_value("class_upper_bound"))
boundary_rule <- config_value("class_boundary_rule")
if (!isTRUE(all.equal(lower_bound, 0.8)) ||
    !isTRUE(all.equal(upper_bound, 1.2)) ||
    !identical(boundary_rule, "outer_inclusive")) {
  stop("Supplementary Figure 5-1 requires the canonical 0.8/1.2 outer-inclusive classification")
}

# Integrated parameter display.
family_order <- c("C01", "C02", "C03")
if (nrow(pairs) != 3L ||
    !identical(as.character(pairs$family), family_order) ||
    !setequal(unique(within$family), family_order) ||
    any(pairs$pair_id[match(within$pair_id, pairs$pair_id)] != within$pair_id)) {
  stop("Supplementary Figure 5-1 family assignment is incomplete or unordered")
}
if (any(within$n_valid != 500L)) {
  stop("Each selected family-parameter cell must contain 500 valid endpoints")
}
if (nrow(master) != 500L * 3L * 14L ||
    any(!is.finite(master$log2_ratio_vivo_to_vitro)) ||
    any(!master$family %in% family_order) ||
    any(!master$parameter %in% parameter_levels())) {
  stop("The endpoint-level ratio table is incomplete or contains invalid rows")
}
master_counts <- aggregate(
  seed_number ~ family + parameter,
  master,
  length
)
if (nrow(master_counts) != 42L || any(master_counts$seed_number != 500L)) {
  stop("Each family-parameter distribution must contain 500 endpoints")
}

parameter_order <- parameter_levels()
within$family <- factor(within$family, levels = family_order)
within$parameter <- factor(
  within$parameter,
  levels = rev(parameter_order)
)
master$family <- factor(master$family, levels = family_order)
master$parameter <- factor(
  master$parameter,
  levels = levels(within$parameter)
)

master_summary <- aggregate(
  log2_ratio_vivo_to_vitro ~ family + parameter,
  master,
  function(values) {
    c(
      n_unique = length(unique(values)),
      median = stats::median(values),
      sd = stats::sd(values)
    )
  }
)
master_summary <- data.frame(
  family = master_summary$family,
  parameter = master_summary$parameter,
  n_unique = master_summary$log2_ratio_vivo_to_vitro[, "n_unique"],
  median = master_summary$log2_ratio_vivo_to_vitro[, "median"],
  sd = master_summary$log2_ratio_vivo_to_vitro[, "sd"],
  stringsAsFactors = FALSE
)
master_class_summary <- aggregate(
  ratio_class ~ family + parameter,
  master,
  function(values) {
    as.integer(table(factor(
      values,
      levels = c("ClassA", "ClassB", "ClassC")
    )))
  }
)
master_class_summary <- data.frame(
  family = master_class_summary$family,
  parameter = master_class_summary$parameter,
  n_ClassA = master_class_summary$ratio_class[, 1L],
  n_ClassB = master_class_summary$ratio_class[, 2L],
  n_ClassC = master_class_summary$ratio_class[, 3L],
  stringsAsFactors = FALSE
)
within_key <- paste(within$family, within$parameter, sep = "::")
summary_key <- paste(
  master_summary$family, master_summary$parameter, sep = "::"
)
summary_match <- match(within_key, summary_key)
class_key <- paste(
  master_class_summary$family,
  master_class_summary$parameter,
  sep = "::"
)
class_match <- match(within_key, class_key)
if (anyNA(summary_match) ||
    any(abs(
      within$log2_ratio_median - master_summary$median[summary_match]
    ) > 1e-12)) {
  stop("Endpoint-level and summarized family medians do not agree")
}
if (anyNA(class_match) ||
    any(within$n_ClassA != master_class_summary$n_ClassA[class_match]) ||
    any(within$n_ClassB != master_class_summary$n_ClassB[class_match]) ||
    any(within$n_ClassC != master_class_summary$n_ClassC[class_match])) {
  stop("Endpoint-level and summarized directional classes do not agree")
}
spread_keys <- summary_key[master_summary$n_unique > 1L]
master_spread <- master[
  paste(master$family, master$parameter, sep = "::") %in% spread_keys,
  ,
  drop = FALSE
]
degenerate_cells <- sum(master_summary$n_unique == 1L)

context_colors <- c("In vivo" = "#0072B2", "In vitro" = "#CC79A7")
class_semantics <- c(
  "ClassA" = "higher in vitro",
  "ClassB" = "approximately equal",
  "ClassC" = "higher in vivo"
)
class_colors <- c(
  "higher in vitro" = context_colors[["In vitro"]],
  "approximately equal" = "#8A9299",
  "higher in vivo" = context_colors[["In vivo"]]
)

composition <- rbind(
  transform(
    within,
    direction = class_semantics[["ClassA"]],
    proportion = prop_ClassA
  ),
  transform(
    within,
    direction = class_semantics[["ClassB"]],
    proportion = prop_ClassB
  ),
  transform(
    within,
    direction = class_semantics[["ClassC"]],
    proportion = prop_ClassC
  )
)
composition$direction <- factor(
  composition$direction,
  levels = unname(class_semantics)
)
composition <- composition[
  order(composition$family, composition$parameter, composition$direction),
  ,
  drop = FALSE
]
composition_sum <- aggregate(
  proportion ~ family + parameter,
  composition,
  sum
)
if (any(abs(composition_sum$proportion - 1) > 1e-10)) {
  stop("Endpoint directional proportions do not sum to one")
}

mixed <- within[
  apply(
    within[, c("prop_ClassA", "prop_ClassB", "prop_ClassC")],
    1,
    max
  ) < 1 - 1e-12,
  ,
  drop = FALSE
]
mixed$composition_label <- apply(
  mixed[, c("prop_ClassA", "prop_ClassB", "prop_ClassC")],
  1,
  function(values) {
    short_names <- c("in vitro", "equal", "in vivo")
    active <- which(values > 0)
    paste0(
      short_names[active], " ", sprintf("%.1f%%", 100 * values[active]),
      collapse = " / "
    )
  }
)

family_colors <- c(
  "C01" = "#C99700",
  "C02" = "#6A3D9A",
  "C03" = "#006D2C"
)
family_shapes <- c("C01" = 16, "C02" = 17, "C03" = 15)

p_composition <- ggplot(
  composition,
  aes(x = proportion, y = parameter, fill = direction)
) +
  geom_col(
    width = 0.72,
    position = position_stack(reverse = TRUE),
    color = "white",
    linewidth = 0.18
  ) +
  geom_text(
    data = mixed,
    aes(x = 0.5, y = parameter, label = composition_label),
    inherit.aes = FALSE,
    size = 2.15,
    fontface = "bold",
    color = "#111827"
  ) +
  facet_grid(. ~ family) +
  scale_fill_manual(
    values = class_colors,
    breaks = unname(class_semantics),
    name = "Optimizer-endpoint class",
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = 0.5,
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Fraction of 500 optimizer endpoints (each bar = 100%)",
    y = NULL
  ) +
  theme_si(8.7) +
  theme(
    axis.text.y = element_text(size = 7.7, face = "bold"),
    axis.text.x = element_text(size = 7.2),
    axis.title.x = element_text(size = 8.2, margin = margin(t = 5)),
    strip.text = element_text(
      size = 8.5,
      face = "bold",
      color = "#26313A"
    ),
    strip.background = element_rect(
      fill = "#F2F4F6",
      color = "#D8DDE2",
      linewidth = 0.35
    ),
    panel.spacing.x = unit(5, "pt"),
    legend.position = "bottom",
    plot.margin = margin(4, 3, 5, 5)
  )

axis_left_bound <- -5.5
axis_break_lower <- 6.2
axis_break_upper <- 12.2
axis_right_bound <- 12.7
tail_display_start <- 6.55
tail_display_scale <- 1.6

broken_axis_transform <- function(values) {
  values <- as.numeric(values)
  transformed <- ifelse(
    values <= axis_break_lower,
    values,
    ifelse(
      values >= axis_break_upper,
      tail_display_start +
        (values - axis_break_upper) * tail_display_scale,
      NA_real_
    )
  )
  transformed
}

if (any(master$log2_ratio_vivo_to_vitro < axis_left_bound) ||
    any(
      master$log2_ratio_vivo_to_vitro > axis_break_lower &
        master$log2_ratio_vivo_to_vitro < axis_break_upper
    ) ||
    any(master$log2_ratio_vivo_to_vitro > axis_right_bound)) {
  stop("The declared broken-axis intervals would omit endpoint data")
}

master$display_log2_ratio <- broken_axis_transform(
  master$log2_ratio_vivo_to_vitro
)
master_spread$display_log2_ratio <- broken_axis_transform(
  master_spread$log2_ratio_vivo_to_vitro
)
within$display_log2_ratio_median <- broken_axis_transform(
  within$log2_ratio_median
)
if (anyNA(master$display_log2_ratio) ||
    anyNA(master_spread$display_log2_ratio) ||
    anyNA(within$display_log2_ratio_median)) {
  stop("Broken-axis transformation produced missing displayed values")
}

tail_display_end <- broken_axis_transform(axis_right_bound)
tail_tick_natural <- max(master$log2_ratio_vivo_to_vitro)
tail_tick_display <- broken_axis_transform(tail_tick_natural)
axis_display_breaks <- c(-5, -2.5, 0, 2.5, 5, tail_tick_display)
axis_display_labels <- c("-5", "-2.5", "0", "+2.5", "+5", "+12.43")

ratio_component_label <-
  "Optimizer-endpoint log-ratio distributions (6.2--12.2 omitted)"
within$ratio_component <- ratio_component_label
master$ratio_component <- ratio_component_label
master_spread$ratio_component <- ratio_component_label
dodge_width <- 0.64
alternating_rows <- data.frame(
  parameter = factor(
    levels(within$parameter)[seq(1L, length(levels(within$parameter)), by = 2L)],
    levels = levels(within$parameter)
  )
)

p_distribution <- ggplot(master, aes(y = parameter)) +
  annotate(
    "rect",
    xmin = axis_left_bound,
    xmax = log2(lower_bound),
    ymin = -Inf,
    ymax = Inf,
    fill = class_colors[["higher in vitro"]],
    alpha = 0.055
  ) +
  annotate(
    "rect",
    xmin = log2(lower_bound),
    xmax = log2(upper_bound),
    ymin = -Inf,
    ymax = Inf,
    fill = class_colors[["approximately equal"]],
    alpha = 0.055
  ) +
  annotate(
    "rect",
    xmin = log2(upper_bound),
    xmax = axis_break_lower,
    ymin = -Inf,
    ymax = Inf,
    fill = class_colors[["higher in vivo"]],
    alpha = 0.055
  ) +
  annotate(
    "rect",
    xmin = tail_display_start,
    xmax = tail_display_end,
    ymin = -Inf,
    ymax = Inf,
    fill = class_colors[["higher in vivo"]],
    alpha = 0.055
  ) +
  geom_tile(
    data = alternating_rows,
    aes(
      x = (axis_left_bound + tail_display_end) / 2,
      y = parameter
    ),
    inherit.aes = FALSE,
    width = tail_display_end - axis_left_bound,
    height = 1,
    fill = "#AEB5BC",
    alpha = 0.14
  ) +
  annotate(
    "rect",
    xmin = axis_break_lower,
    xmax = tail_display_start,
    ymin = -Inf,
    ymax = Inf,
    fill = "white",
    alpha = 1
  ) +
  annotate(
    "text",
    x = (axis_break_lower + tail_display_start) / 2,
    y = 1,
    label = "//",
    size = 3.0,
    fontface = "bold",
    color = "#6B747C"
  ) +
  geom_violin(
    data = master_spread,
    aes(
      x = display_log2_ratio,
      fill = family,
      color = family,
      group = interaction(parameter, family)
    ),
    orientation = "y",
    position = position_dodge(width = dodge_width),
    width = 0.56,
    scale = "width",
    trim = TRUE,
    alpha = 0.26,
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = c(log2(lower_bound), log2(upper_bound)),
    color = "#6B747C",
    linetype = 3,
    linewidth = 0.45
  ) +
  geom_vline(xintercept = 0, color = "#26313A", linewidth = 0.55) +
  geom_point(
    data = within,
    aes(
      x = display_log2_ratio_median,
      color = family,
      shape = family,
      group = family
    ),
    position = position_dodge(width = dodge_width),
    size = 1.75,
    stroke = 0.35
  ) +
  facet_grid(. ~ ratio_component) +
  scale_fill_manual(
    values = family_colors,
    breaks = family_order,
    guide = "none",
    drop = FALSE
  ) +
  scale_color_manual(
    values = family_colors,
    breaks = family_order,
    name = "C family",
    drop = FALSE
  ) +
  scale_shape_manual(
    values = family_shapes,
    breaks = family_order,
    name = "C family",
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = levels(within$parameter),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(axis_left_bound, tail_display_end),
    breaks = axis_display_breaks,
    labels = axis_display_labels,
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  labs(
    x = expression(log[2] * "(in vivo / in vitro)" * "; broken axis"),
    y = NULL
  ) +
  theme_si(8.7) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 7.2),
    axis.title.x = element_text(size = 8.2, margin = margin(t = 5)),
    strip.text = element_text(
      size = 8.5,
      face = "bold",
      color = "#26313A"
    ),
    strip.background = element_rect(
      fill = "#F2F4F6",
      color = "#D8DDE2",
      linewidth = 0.35
    ),
    legend.position = "bottom",
    plot.margin = margin(4, 5, 5, 3)
  )

supp_figure5_1 <- (p_composition | p_distribution) +
  plot_layout(widths = c(1.0, 1.5), guides = "collect") +
  plot_annotation(
    title = "Joint-fit parameter differences across three selected primary families",
    subtitle = paste0(
      "Left: endpoint-level directional composition under the 0.8/1.2 natural-scale thresholds. ",
      "Right: horizontal distributions retain within-family spread and contrast magnitude."
    ),
    caption = paste0(
      "Each cell contains 500 endpoints. Pink/gray/blue denote higher in vitro/equal/higher in vivo; ",
      "mixed-cell labels give percentages. Violin width shows endpoint density and symbols mark medians.\n",
      "A symbol without a violin denotes 500 identical ratios. The broken x axis omits the empty ",
      "6.2--12.2 interval but retains the 12.43 endpoint group. These are numerical-search summaries, not posterior draws."
    ),
    theme = theme(
      plot.title = element_text(size = 11.5, face = "bold", color = "#111827"),
      plot.subtitle = element_text(size = 7.3, color = "#4B5563"),
      plot.caption = element_text(
        size = 6.8,
        color = "#6B747C",
        hjust = 0
      ),
      plot.margin = margin(4, 4, 4, 4)
    )
  ) &
  theme(legend.position = "bottom")

embed_only <- identical(Sys.getenv("SUPP_FIGURE5_1_EMBED_ONLY"), "1")
if (!embed_only) {
save_plot_pair(
  p_composition,
  file.path(panel_root, "supp_figure5-1_directional_composition"),
  width = 7.6,
  height = 6.8
)
save_plot_pair(
  p_distribution,
  file.path(panel_root, "supp_figure5-1_endpoint_log_ratio_distributions"),
  width = 6.1,
  height = 6.8
)

output_stub <- file.path(figure_root, "supp_fig5-1_joint_parameter_stability")
save_plot_pair(supp_figure5_1, output_stub, width = 13.6, height = 7.2)
for (extension in c("png", "pdf")) {
  source_path <- normalizePath(
    paste0(output_stub, ".", extension), mustWork = TRUE
  )
  destination_path <- normalizePath(
    file.path(
      revised_root,
      paste0("supp_fig5-1_joint_parameter_stability.", extension)
    ),
    mustWork = FALSE
  )
  if (!identical(source_path, destination_path)) {
    file.copy(source_path, destination_path, overwrite = TRUE)
  }
}

provenance <- data.frame(
  role = c(
    "within-pair parameter stability",
    "between-pair parameter stability",
    "classification configuration",
    "selected primary-family pairs",
    "endpoint-level log-ratio distributions"
  ),
  path = normalizePath(input_paths, mustWork = TRUE),
  md5 = unname(tools::md5sum(input_paths)),
  stringsAsFactors = FALSE
)
write.table(
  provenance,
  file.path(input_root, "supp_figure5-1_source_file_provenance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ordered_rows <- data.frame(
  display_order_top_to_bottom = seq_along(parameter_order),
  parameter = parameter_order,
  stringsAsFactors = FALSE
)
write.table(
  ordered_rows,
  file.path(input_root, "supp_figure5-1_display_order.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

validation <- data.frame(
  check = c(
    "selected_family_parameter_rows",
    "between_parameter_rows",
    "warm_start_pairs",
    "class_lower_bound",
    "class_upper_bound",
    "boundary_rule",
    "optimizer_endpoints_per_family_parameter",
    "composition_rows",
    "composition_sums_to_one",
    "class_labels_absent",
    "direction_labels",
    "integrated_layout",
    "three_family_composition_columns",
    "endpoint_distribution_rows",
    "endpoints_per_distribution",
    "distribution_medians_match_summary",
    "distribution_classes_match_summary",
    "parameter_row_alignment",
    "broken_axis_left_omitted_rows",
    "broken_axis_gap_rows",
    "broken_axis_tail_rows",
    "nondegenerate_violin_cells",
    "degenerate_point_cells",
    "family_identity_colors",
    "assembled_png",
    "assembled_pdf"
  ),
  observed = c(
    nrow(within),
    nrow(between),
    length(unique(within$pair_id)),
    lower_bound,
    upper_bound,
    boundary_rule,
    paste(sort(unique(within$n_valid)), collapse = ","),
    nrow(composition),
    all(abs(composition_sum$proportion - 1) <= 1e-10),
    !any(grepl("^Class[A-C]$", levels(composition$direction))),
    paste(levels(composition$direction), collapse = " | "),
    "composition_left_endpoint_distributions_right",
    identical(levels(within$family), family_order),
    nrow(master),
    paste(sort(unique(master_counts$seed_number)), collapse = ","),
    all(abs(
      within$log2_ratio_median - master_summary$median[summary_match]
    ) <= 1e-12),
    all(within$n_ClassA == master_class_summary$n_ClassA[class_match]) &&
      all(within$n_ClassB == master_class_summary$n_ClassB[class_match]) &&
      all(within$n_ClassC == master_class_summary$n_ClassC[class_match]),
    identical(levels(within$parameter), levels(master$parameter)),
    sum(master$log2_ratio_vivo_to_vitro < axis_left_bound),
    sum(
      master$log2_ratio_vivo_to_vitro > axis_break_lower &
        master$log2_ratio_vivo_to_vitro < axis_break_upper
    ),
    sum(master$log2_ratio_vivo_to_vitro >= axis_break_upper),
    sum(master_summary$n_unique > 1L),
    degenerate_cells,
    paste(names(family_colors), family_colors, collapse = ";"),
    file.exists(paste0(output_stub, ".png")),
    file.exists(paste0(output_stub, ".pdf"))
  ),
  expected = c(
    "42", "14", "3", "0.8", "1.2", "outer_inclusive",
    "500", "126", "TRUE", "TRUE",
    "higher in vitro | approximately equal | higher in vivo",
    "composition_left_endpoint_distributions_right", "TRUE",
    "21000", "500", "TRUE", "TRUE", "TRUE",
    "0", "0", "500", "23", "19",
    "C01 #C99700;C02 #6A3D9A;C03 #006D2C", "TRUE", "TRUE"
  ),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(input_root, "supp_figure5-1_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Supplementary Figure 5-1 generation complete.\n")
cat("Integrated display: directional composition plus endpoint distributions.\n")
cat("Output:", paste0(output_stub, ".png"), "\n")
}

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Supp_Figure5_1 <- function() {
  data_dir <- file.path(DATA_ROOT, "Supp_Figure5_1")
  require_files(file.path(data_dir, c(
    "within_pair_parameter_stability.tsv",
    "between_pair_parameter_stability.tsv",
    "analysis_config.tsv",
    "selected_primary_family_pairs.tsv",
    "soft_coupling_master_long.tsv"
  )), "Supplementary Figure 5-1 intermediate")
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_Supp_Figure5_1.R"), mustWork = TRUE),
    env = c(
      "SUPP_FIGURE5_1_DRAW_WORKER=1",
      paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("SUPP_FIGURE5_1_DATA_DIR=", data_dir),
      paste0("SUPP_FIGURE5_1_FIGURE_DIR=", OUTPUT_ROOT),
      paste0("SUPP_FIGURE5_1_PANEL_DIR=", file.path(data_dir, "panels"))
    )
  )
  require_files(
    file.path(OUTPUT_ROOT, c(
      "supp_fig5-1_joint_parameter_stability.png",
      "supp_fig5-1_joint_parameter_stability.pdf"
    )),
    "Supplementary Figure 5-1 output"
  )
  invisible(file.path(
    OUTPUT_ROOT, "supp_fig5-1_joint_parameter_stability.png"
  ))
}

if (sys.nframe() == 0L) draw_Supp_Figure5_1()

}
