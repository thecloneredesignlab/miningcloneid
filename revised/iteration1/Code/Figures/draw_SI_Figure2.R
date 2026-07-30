#!/usr/bin/env Rscript

if (identical(Sys.getenv("SI_FIGURE2_DRAW_WORKER"), "1")) {

# Standalone iteration3 Supplementary Figure 2 generator.
#
# A. Pair-by-parameter median log2(in vivo / in vitro) ratio map.
# B. Between-pair directional consensus under the declared 0.8/1.2
#    outer-inclusive classification rule.
#
# All plotting helpers are defined here. The script reads only the frozen
# canonical analysis tables stored under iteration3 and performs no refitting.

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
    "SI_FIGURE2_DATA_DIR",
    unset = file.path(
      workspace_root, "data", "Figures", "SI_Figure2"
    )
  ),
  mustWork = TRUE
)
figure_root <- normalizePath(
  Sys.getenv(
    "SI_FIGURE2_FIGURE_DIR",
    unset = file.path(workspace_root, "Figures")
  ),
  mustWork = FALSE
)
panel_root <- normalizePath(
  Sys.getenv(
    "SI_FIGURE2_PANEL_DIR",
    unset = file.path(figure_root, "supp_figure2_subpanels")
  ),
  mustWork = FALSE
)
revised_root <- figure_root
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_root, recursive = TRUE, showWarnings = FALSE)
dir.create(revised_root, recursive = TRUE, showWarnings = FALSE)

within_path <- file.path(input_root, "within_pair_parameter_stability.tsv")
between_path <- file.path(input_root, "between_pair_parameter_stability.tsv")
config_path <- file.path(input_root, "analysis_config.tsv")
pair_path <- file.path(input_root, "ploidy_pair_category_assignment.tsv")
input_paths <- c(within_path, between_path, config_path, pair_path)
if (any(!file.exists(input_paths))) {
  stop("Missing Supplementary Figure 2 input(s):\n", paste(input_paths[!file.exists(input_paths)], collapse = "\n"))
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

pair_label <- function(pair_id) {
  pair_id <- as.character(pair_id)
  pattern <- "^fit_joint_.+_vi_seed([0-9]+)_C([0-9]+)Sc([0-9]+)_vt_seed([0-9]+)$"
  match_object <- regexec(pattern, pair_id)
  values <- regmatches(pair_id, match_object)
  vapply(seq_along(pair_id), function(i) {
    hit <- values[[i]]
    if (length(hit) != 5L) {
      pair_id[[i]]
    } else {
      sprintf("C%02dSc%02d / vi%s", as.integer(hit[[3L]]), as.integer(hit[[4L]]), hit[[2L]])
    }
  }, character(1L))
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

required_within <- c("pair_id", "parameter", "log2_ratio_median")
required_between <- c(
  "parameter", "n_pairs", "cross_pair_dominant_class",
  "cross_pair_dominant_fraction"
)
if (!all(required_within %in% names(within)) ||
    !all(required_between %in% names(between))) {
  stop("Supplementary Figure 2 tables lack required fields")
}
if (nrow(within) != 84L || nrow(between) != 14L) {
  stop("Expected 84 pair-parameter rows and 14 between-pair rows")
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
  stop("SI Figure 2 requires the canonical 0.8/1.2 outer-inclusive classification")
}

# Panel A: pair-by-parameter median ratios.
within$pair_label <- pair_label(within$pair_id)
within$ploidy_category <- pairs$pair_ploidy_category[
  match(within$pair_id, pairs$pair_id)
]
if (any(is.na(within$ploidy_category))) {
  stop("Could not map every warm-start pair to its final ploidy category")
}
within$pair_display <- paste0(
  within$pair_label,
  "\npattern ",
  sub("^Cat", "", within$ploidy_category),
  ""
)
pair_order <- unique(within[
  order(
    match(within$ploidy_category, c("CatA", "CatB", "CatC")),
    within$pair_label
  ),
  "pair_display"
])
within$pair_display <- factor(within$pair_display, levels = pair_order)
within$parameter <- factor(within$parameter, levels = rev(parameter_levels()))

ratio_bound <- max(abs(within$log2_ratio_median), na.rm = TRUE)
within$ratio_label_color <- ifelse(
  abs(within$log2_ratio_median) >= 0.55 * ratio_bound,
  "white",
  "#26313A"
)
context_colors <- c("In vivo" = "#0072B2", "In vitro" = "#CC79A7")

p_a <- ggplot(
  within,
  aes(x = pair_display, y = parameter, fill = log2_ratio_median)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%+.2f", log2_ratio_median), color = ratio_label_color),
    size = 1.85
  ) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = context_colors[["In vitro"]],
    mid = "#FFFDF7",
    high = context_colors[["In vivo"]],
    midpoint = 0,
    limits = c(-ratio_bound, ratio_bound),
    name = "Median log2\n(in vivo / in vitro)"
  ) +
  labs(
    title = "A  Pair-level parameter-ratio map",
    subtitle = "Each cell is the median across 500 numerical seeds within one warm-start pair",
    x = NULL,
    y = NULL,
    caption = NULL
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  theme_si(8) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 6.0,
      lineheight = 0.9
    ),
    axis.text.y = element_text(size = 7.2, face = "bold"),
    legend.position = "right",
    plot.title = element_text(size = 9.2),
    plot.subtitle = element_text(size = 7)
  )

# Panel B: semantic directional classes and class-block ordering.
class_semantics <- c(
  "ClassA" = "lower in vivo",
  "ClassB" = "approximately equal",
  "ClassC" = "higher in vivo"
)
class_colors <- c(
  "lower in vivo" = context_colors[["In vitro"]],
  "approximately equal" = "#8A9299",
  "higher in vivo" = context_colors[["In vivo"]]
)
between$direction <- unname(class_semantics[between$cross_pair_dominant_class])
if (any(is.na(between$direction))) {
  stop("Unexpected cross-pair dominant class in canonical stability table")
}
between$direction <- factor(
  between$direction,
  levels = c("lower in vivo", "approximately equal", "higher in vivo")
)
parameter_rank <- match(between$parameter, parameter_levels())
ordering <- between[
  order(as.integer(between$direction), parameter_rank),
  c("parameter", "direction"),
  drop = FALSE
]
between$parameter <- factor(
  between$parameter,
  levels = rev(as.character(ordering$parameter))
)
between$pair_count_label <- sprintf(
  "%.0f/%d pairs",
  between$cross_pair_dominant_fraction * between$n_pairs,
  between$n_pairs
)

direction_counts <- table(factor(
  as.character(ordering$direction),
  levels = levels(between$direction)
))
direction_separators <- nrow(ordering) -
  cumsum(as.integer(direction_counts))[seq_len(length(direction_counts) - 1L)] +
  0.5

p_b <- ggplot(
  between,
  aes(
    x = cross_pair_dominant_fraction,
    y = parameter,
    color = direction
  )
) +
  geom_segment(
    aes(x = 0, xend = cross_pair_dominant_fraction, yend = parameter),
    color = "#D8DDE2",
    linewidth = 1.2
  ) +
  geom_hline(
    yintercept = direction_separators,
    color = "#69747D",
    linewidth = 0.55
  ) +
  geom_point(size = 3.1) +
  geom_text(
    aes(label = pair_count_label),
    hjust = -0.14,
    size = 2.7,
    color = "#26313A"
  ) +
  geom_vline(
    xintercept = c(0.8, 0.9, 0.95),
    linetype = c(3, 2, 3),
    color = "#8A9299",
    linewidth = 0.3
  ) +
  scale_color_manual(
    values = class_colors,
    breaks = levels(between$direction),
    name = "Directional class",
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 1.18),
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "B  Between-pair directional consensus",
    subtitle = paste0(
      "Category blocks: lower in vivo (0 < ratio <= 0.8), approximately equal ",
      "(0.8 < ratio < 1.2), higher in vivo (ratio >= 1.2)"
    ),
    x = "Warm-start pairs sharing the parameter's dominant directional class",
    y = NULL,
    caption = NULL
  ) +
  theme_si(8.5) +
  theme(
    aspect.ratio = 1,
    axis.text.y = element_text(size = 7.3, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(size = 9.2),
    plot.subtitle = element_text(size = 6.9),
    legend.position = "right"
  )

save_plot_pair(
  p_a,
  file.path(panel_root, "supp_figure2a_pair_parameter_ratio_map"),
  width = 6.0,
  height = 7.1
)
save_plot_pair(
  p_b,
  file.path(panel_root, "supp_figure2b_between_pair_directional_consensus"),
  width = 6.5,
  height = 7.1
)

supp_figure2 <- (p_a | p_b) +
  plot_layout(widths = c(0.92, 1.08)) +
  plot_annotation(
    title = "Joint-fit scalar parameter contrasts across six warm-start pairs",
    subtitle = paste0(
      "Pair-level medians and cross-pair directional agreement are optimizer-derived summaries, ",
      "not biological confidence intervals."
    ),
    caption = paste0(
      "A: Pink indicates lower in vivo and blue indicates higher in vivo; trajectory-pattern labels ",
      "are descriptive model-output categories. B: The warm-start pair is the summarization unit; ",
      "one pair changes the fraction by 16.7 percentage points."
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
  )

output_stub <- file.path(figure_root, "supp_fig2_joint_parameter_stability")
save_plot_pair(supp_figure2, output_stub, width = 13.6, height = 6.8)
for (extension in c("png", "pdf")) {
  source_path <- normalizePath(
    paste0(output_stub, ".", extension), mustWork = TRUE
  )
  destination_path <- normalizePath(
    file.path(
      revised_root,
      paste0("supp_fig2_joint_parameter_stability.", extension)
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
    "pair trajectory categories"
  ),
  path = normalizePath(input_paths, mustWork = TRUE),
  md5 = unname(tools::md5sum(input_paths)),
  stringsAsFactors = FALSE
)
write.table(
  provenance,
  file.path(input_root, "si_figure2_source_file_provenance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ordered_rows <- data.frame(
  display_order_top_to_bottom = seq_len(nrow(ordering)),
  parameter = ordering$parameter,
  direction = as.character(ordering$direction),
  stringsAsFactors = FALSE
)
write.table(
  ordered_rows,
  file.path(input_root, "si_figure2b_display_order.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

validation <- data.frame(
  check = c(
    "pair_parameter_rows",
    "between_parameter_rows",
    "warm_start_pairs",
    "class_lower_bound",
    "class_upper_bound",
    "boundary_rule",
    "direction_block_order",
    "class_labels_absent",
    "panel_arrangement",
    "panel_a_square_cells",
    "panel_b_square_plot_area",
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
    paste(unique(as.character(ordering$direction)), collapse = " | "),
    !any(grepl("^Class[A-C]$", as.character(between$direction))),
    "left_to_right",
    TRUE,
    TRUE,
    file.exists(paste0(output_stub, ".png")),
    file.exists(paste0(output_stub, ".pdf"))
  ),
  expected = c(
    "84", "14", "6", "0.8", "1.2", "outer_inclusive",
    "lower in vivo | approximately equal | higher in vivo",
    "TRUE", "left_to_right", "TRUE", "TRUE", "TRUE", "TRUE"
  ),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(input_root, "si_figure2_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Supplementary Figure 2 generation complete.\n")
cat("Panel B directional order:", paste(unique(as.character(ordering$direction)), collapse = " -> "), "\n")
cat("Output:", paste0(output_stub, ".png"), "\n")

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_SI_Figure2 <- function() {
  data_dir <- file.path(DATA_ROOT, "SI_Figure2")
  require_files(file.path(data_dir, c(
    "within_pair_parameter_stability.tsv",
    "between_pair_parameter_stability.tsv",
    "analysis_config.tsv",
    "ploidy_pair_category_assignment.tsv"
  )), "SI Figure 2 intermediate")
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_SI_Figure2.R"), mustWork = TRUE),
    env = c(
      "SI_FIGURE2_DRAW_WORKER=1",
      paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("SI_FIGURE2_DATA_DIR=", data_dir),
      paste0("SI_FIGURE2_FIGURE_DIR=", OUTPUT_ROOT),
      paste0("SI_FIGURE2_PANEL_DIR=", file.path(data_dir, "panels"))
    )
  )
  require_files(
    file.path(OUTPUT_ROOT, c(
      "supp_fig2_joint_parameter_stability.png",
      "supp_fig2_joint_parameter_stability.pdf"
    )),
    "SI Figure 2 output"
  )
  invisible(file.path(
    OUTPUT_ROOT, "supp_fig2_joint_parameter_stability.png"
  ))
}

if (sys.nframe() == 0L) draw_SI_Figure2()

}
