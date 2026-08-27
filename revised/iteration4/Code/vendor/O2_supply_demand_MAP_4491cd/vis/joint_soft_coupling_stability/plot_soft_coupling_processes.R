#!/usr/bin/env Rscript

# Render biological-process x ratio-class summaries from analysis tables.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(SCRIPT_DIR, "joint_coupling_visualization_utils.R"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(name, default = NULL) { hit <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default }
  analysis_dir <- get_arg("analysis_dir") %||% stop("--analysis_dir is required", call. = FALSE)
  out_dir <- get_arg("out_dir") %||% file.path(analysis_dir, "figures")
  input_path <- file.path(analysis_dir, "parameter_class_process_summary.tsv")
  process_path <- file.path(analysis_dir, "process_class_summary.tsv")
  data <- o2jcv_read_tsv(input_path)
  process <- o2jcv_read_tsv(process_path)
  data$parameter <- factor(data$parameter, levels = rev(o2jcv_parameter_levels()))
  data$ratio_class <- factor(data$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  data$label_colour <- ifelse(data$pair_balanced_mean_proportion >= 0.65, "white", "#26313A")
  p <- ggplot2::ggplot(data, ggplot2::aes(x = ratio_class, y = parameter, fill = pair_balanced_mean_proportion)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * pair_balanced_mean_proportion), colour = label_colour), size = 2.7) +
    ggplot2::scale_colour_identity() +
    ggplot2::facet_grid(
      primary_process ~ ., scales = "free_y", space = "free_y",
      labeller = ggplot2::labeller(primary_process = ggplot2::label_wrap_gen(width = 18))
    ) +
    ggplot2::scale_fill_gradient(low = "#FFF9E8", high = "#A97616", labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Pair-balanced\nseed proportion", title = "Soft-coupling classes by biological process", subtitle = "Process labels are audited from the model parameter table") +
    o2jcv_theme(9) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5, lineheight = 0.95, size = 8),
      plot.margin = ggplot2::margin(10, 20, 10, 10)
    )
  o2jcv_save_figure(p, "process_parameter_class_map", out_dir, c(input_path), 11.5, 9, family = "Matrix & Cohort", question = "Which biological processes contain consistently lower, coupled, or higher in-vivo parameters?")

  process$ratio_class <- factor(process$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  process$primary_process <- factor(process$primary_process, levels = rev(unique(process$primary_process)))
  p2 <- ggplot2::ggplot(process, ggplot2::aes(x = mean_pair_balanced_proportion, y = primary_process, colour = ratio_class)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = mean_pair_balanced_proportion, yend = primary_process), colour = "#D8DDE2", linewidth = 1) +
    ggplot2::geom_point(size = 3) + ggplot2::facet_wrap(~ratio_class, nrow = 1) +
    ggplot2::scale_colour_manual(values = o2jcv_class_colors(), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Mean pair-balanced seed proportion across parameters", y = NULL, title = "Biological-process class profile", subtitle = "Each point averages parameters assigned to the same audited primary process") + o2jcv_theme(8.5, rotate_x = FALSE)
  o2jcv_save_figure(p2, "process_class_profile", out_dir, c(process_path), 13, 5.8, family = "Comparison & Ranking", question = "Which biological processes are most associated with ClassA, ClassB, or ClassC behavior?")
}

if (identical(environment(), globalenv())) main()
