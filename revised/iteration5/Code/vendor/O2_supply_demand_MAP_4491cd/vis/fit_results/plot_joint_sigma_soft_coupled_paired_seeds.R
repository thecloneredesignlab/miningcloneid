#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) dirname(frame_files[[length(frame_files)]]) else getwd()
})
SCRIPT_DIR <- normalizePath(.script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
rm(.script_dir)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)
`%||%` <- o2fr_null_coalesce

TARGET_PARAMETERS <- c("p_misseg", "O2_crit", "alpha_o2", "mu_hp")
CONTEXT_LEVELS <- c("in vivo", "in vitro")
OBJECTIVE_TYPES <- c("objective", "objective_invivo", "objective_invitro")
OUTPUT_PREFIX <- "joint_sigma_soft_coupled_paired_seed_comparison"

parse_args <- o2fr_parse_args

safe_file_token <- function(x) {
  x <- gsub("\\+", "plus", as.character(x))
  x <- gsub("-", "minus", x)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  gsub("_+", "_", x)
}

value_group_labels <- function(levels) {
  vapply(strsplit(as.character(levels), "__", fixed = TRUE), function(x) {
    paste0(x[[1]], "\n", x[[2]])
  }, character(1))
}

plot_objective_overview <- function(objective_df) {
  plot_df <- objective_df[is.finite(objective_df$objective_value), , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)
  set.seed(1)
  ggplot(plot_df, aes(x = objective_group, y = objective_value)) +
    geom_boxplot(
      fill = "white",
      color = "grey45",
      linewidth = 0.45,
      outlier.shape = NA,
      width = 0.58
    ) +
    geom_line(
      aes(group = interaction(seed, objective_type, drop = TRUE)),
      color = "grey72",
      linewidth = 0.28,
      alpha = 0.55
    ) +
    geom_point(
      position = position_jitter(width = 0.055, height = 0),
      color = "grey20",
      size = 1.55,
      alpha = 0.82
    ) +
    scale_x_discrete(labels = value_group_labels(levels(plot_df$objective_group))) +
    labs(
      title = "Joint objective components across paired seeds",
      subtitle = "Objective type is the primary group; within each type, sigma values are subgroups. Lines connect the same seed across sigma values.",
      x = NULL,
      y = "Objective value"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8.7),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

plot_parameter_values <- function(value_df, parameter) {
  plot_df <- value_df[as.character(value_df$parameter) == parameter, , drop = FALSE]
  plot_df <- plot_df[is.finite(plot_df$value), , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)
  set.seed(1)
  p <- ggplot(plot_df, aes(x = value_group, y = value)) +
    geom_boxplot(
      fill = "white",
      color = "grey45",
      linewidth = 0.45,
      outlier.shape = NA,
      width = 0.58
    ) +
    geom_line(
      aes(group = interaction(seed, context, drop = TRUE)),
      color = "grey72",
      linewidth = 0.28,
      alpha = 0.55
    )
  p <- p +
    geom_point(
      position = position_jitter(width = 0.055, height = 0),
      color = "grey20",
      size = 1.55,
      alpha = 0.82
    )
  p +
    scale_x_discrete(labels = value_group_labels(levels(plot_df$value_group))) +
    labs(
      title = paste0(parameter, ": in vivo and in vitro fitted values"),
      subtitle = "Within each context, lines connect the same seed across sigma values.",
      x = NULL,
      y = "Natural-scale parameter value"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8.7),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

plot_parameter_ratio <- function(df, parameter) {
  plot_df <- df[as.character(df$parameter) == parameter, , drop = FALSE]
  plot_df <- plot_df[is.finite(plot_df$ratio_vivo_to_vitro), , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)
  set.seed(1)
  p <- ggplot(plot_df, aes(x = sigma_label, y = ratio_vivo_to_vitro)) +
    geom_hline(yintercept = 1, color = "grey70", linetype = "dashed", linewidth = 0.35) +
    geom_boxplot(
      fill = "white",
      color = "grey45",
      linewidth = 0.45,
      outlier.shape = NA,
      width = 0.52
    ) +
    geom_line(
      aes(group = seed),
      color = "grey72",
      linewidth = 0.28,
      alpha = 0.55
    )
  p <- p +
    geom_point(
      position = position_jitter(width = 0.045, height = 0),
      color = "grey20",
      size = 1.55,
      alpha = 0.82
    )
  p +
    labs(
      title = paste0(parameter, ": vivo/vitro ratio"),
      subtitle = "Lines connect the same seed across sigma values.",
      x = NULL,
      y = "ratio_vivo_to_vitro"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

save_combined_plot <- function(left_plot, right_plot, out_png, out_pdf, width = 14, height = 6.2) {
  if (is.null(left_plot) || is.null(right_plot)) return(invisible(NULL))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(left_plot, right_plot, ncol = 2)
    ggplot2::ggsave(out_png, combined, width = width, height = height, dpi = 180)
    ggplot2::ggsave(out_pdf, combined, width = width, height = height)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined <- gridExtra::arrangeGrob(left_plot, right_plot, ncol = 2)
    ggplot2::ggsave(out_png, combined, width = width, height = height, dpi = 180)
    ggplot2::ggsave(out_pdf, combined, width = width, height = height)
  } else {
    stop("Need either patchwork or gridExtra to save side-by-side combined plots.", call. = FALSE)
  }
  invisible(out_png)
}

save_single_plot <- function(plot, out_png, out_pdf, width = 12, height = 6.2) {
  if (is.null(plot)) return(invisible(NULL))
  ggplot2::ggsave(out_png, plot, width = width, height = height, dpi = 180)
  ggplot2::ggsave(out_pdf, plot, width = width, height = height)
  invisible(out_png)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  analysis_dir <- normalizePath(argv$analysis_dir %||% argv$out_dir, mustWork = TRUE)
  if (!file.exists(file.path(analysis_dir, "analysis_manifest.tsv"))) {
    stop("Missing joint-sigma analysis manifest: ", analysis_dir, call. = FALSE)
  }
  viz_dir <- normalizePath(argv$viz_dir %||% argv$out_dir %||% analysis_dir, mustWork = FALSE)
  prefix <- OUTPUT_PREFIX
  paths <- file.path(analysis_dir, paste0(prefix, c("_summary.tsv", "_value_long.tsv", "_objective_long.tsv")))
  if (any(!file.exists(paths))) stop("Missing joint-sigma analysis tables: ", paste(paths[!file.exists(paths)], collapse = ", "), call. = FALSE)
  paired_df <- utils::read.delim(paths[[1]], check.names = FALSE, stringsAsFactors = FALSE)
  value_long <- utils::read.delim(paths[[2]], check.names = FALSE, stringsAsFactors = FALSE)
  objective_long <- utils::read.delim(paths[[3]], check.names = FALSE, stringsAsFactors = FALSE)
  sigma_levels <- unique(as.character(paired_df$sigma_label[order(suppressWarnings(as.numeric(paired_df$sigma)))]))
  value_long$sigma_label <- factor(value_long$sigma_label, levels = sigma_levels)
  value_groups <- as.vector(t(outer(CONTEXT_LEVELS, sigma_levels, paste, sep = "__")))
  value_long$value_group <- factor(value_long$value_group, levels = value_groups)
  objective_long$sigma_label <- factor(objective_long$sigma_label, levels = sigma_levels)
  objective_groups <- as.vector(t(outer(OBJECTIVE_TYPES, sigma_levels, paste, sep = "__")))
  objective_long$objective_group <- factor(objective_long$objective_group, levels = objective_groups)
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

  objective_plot <- plot_objective_overview(objective_long)
  objective_png <- file.path(viz_dir, paste0(prefix, "_objectives.png"))
  objective_pdf <- file.path(viz_dir, paste0(prefix, "_objectives.pdf"))
  save_single_plot(objective_plot, objective_png, objective_pdf)
  outputs <- c(objective_png, objective_pdf)
  for (param in TARGET_PARAMETERS) {
    token <- safe_file_token(param)
    out_png <- file.path(viz_dir, paste0(prefix, "_", token, "_values_and_ratio.png"))
    out_pdf <- file.path(viz_dir, paste0(prefix, "_", token, "_values_and_ratio.pdf"))
    save_combined_plot(plot_parameter_values(value_long, param), plot_parameter_ratio(paired_df, param), out_png, out_pdf)
    outputs <- c(outputs, out_png, out_pdf)
  }
  manifest <- data.frame(stage = "visualization", file = basename(outputs[file.exists(outputs)]), stringsAsFactors = FALSE)
  utils::write.table(manifest, file.path(viz_dir, "visualization_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote joint-sigma visualizations to: ", viz_dir)
  invisible(viz_dir)
}

if (sys.nframe() == 0L) main()
