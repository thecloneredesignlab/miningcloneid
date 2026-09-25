#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)
parse_args <- o2fr_parse_args
`%||%` <- o2fr_null_coalesce

main <- function(argv = parse_args()) {
  analysis_dir <- normalizePath(argv$analysis_dir %||% argv$extra_results_dir, mustWork = TRUE)
  if (!file.exists(file.path(analysis_dir, "analysis_manifest.tsv"))) stop("Missing extra-results analysis manifest.", call. = FALSE)
  data_path <- file.path(analysis_dir, "objective_components_long.tsv")
  if (!file.exists(data_path)) stop("Missing objective component table: ", data_path, call. = FALSE)
  out_path <- normalizePath(argv$out_path %||% file.path(analysis_dir, "objective_components_violin.pdf"), mustWork = FALSE)
  plot_df <- utils::read.delim(data_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(plot_df)) stop("No finite objective rows in: ", data_path, call. = FALSE)
  plot_df$metric_label <- factor(plot_df$metric_label, levels = unique(as.character(plot_df$metric_label)))
  p <- ggplot(plot_df, aes(metric_label, value, fill = metric_label)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, linewidth = 0.35) +
    labs(title = "Objective Distributions Across Seeds", x = NULL, y = "Objective value") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "none")
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_path, p, width = 9, height = 7)
  message("Wrote objective-components violin: ", out_path)
  invisible(out_path)
}

if (sys.nframe() == 0L) main()
