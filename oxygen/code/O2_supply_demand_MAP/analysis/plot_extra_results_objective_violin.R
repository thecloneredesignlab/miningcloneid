#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

.o2sd_bootstrap_script_dir <- local({
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
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2sd_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

parse_args <- o2sd_parse_args

default_extra_results_dir <- "/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/fit_invivo_o2_supply_demand_MAP/extra_results"

build_long_objective_table <- function(seed_summary) {
  metric_cols <- c("objective", "objective_burden", "objective_ploidy")
  missing_cols <- setdiff(metric_cols, names(seed_summary))
  if (length(missing_cols) > 0L) {
    stop("seed_summary.tsv is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  rows <- lapply(metric_cols, function(metric_name) {
    values <- suppressWarnings(as.numeric(seed_summary[[metric_name]]))
    data.frame(
      seed = as.character(seed_summary$seed),
      metric = metric_name,
      value = values,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[is.finite(out$value), , drop = FALSE]
  out$metric <- factor(out$metric, levels = metric_cols)
  out$metric_label <- factor(
    out$metric,
    levels = metric_cols,
    labels = c("objective", "objective_burden", "objective_ploidy")
  )
  out
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  extra_results_dir <- if (!is.null(argv$extra_results_dir)) {
    normalizePath(argv$extra_results_dir, mustWork = TRUE)
  } else {
    normalizePath(default_extra_results_dir, mustWork = TRUE)
  }

  seed_summary_path <- file.path(extra_results_dir, "seed_summary.tsv")
  if (!file.exists(seed_summary_path)) {
    stop("Missing seed_summary.tsv in: ", extra_results_dir)
  }

  out_path <- if (!is.null(argv$out_path) && nzchar(trimws(argv$out_path))) {
    normalizePath(argv$out_path, mustWork = FALSE)
  } else {
    file.path(extra_results_dir, "objective_components_violin.pdf")
  }
  out_table <- if (!is.null(argv$out_table) && nzchar(trimws(argv$out_table))) {
    normalizePath(argv$out_table, mustWork = FALSE)
  } else {
    file.path(extra_results_dir, "objective_components_long.tsv")
  }

  seed_summary <- utils::read.delim(seed_summary_path, check.names = FALSE, stringsAsFactors = FALSE)
  plot_df <- build_long_objective_table(seed_summary)
  if (!nrow(plot_df)) {
    stop("No finite objective values were found in: ", seed_summary_path)
  }

  utils::write.table(
    plot_df,
    file = out_table,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p <- ggplot(plot_df, aes(x = metric_label, y = value, fill = metric_label)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, linewidth = 0.35) +
    scale_fill_manual(
      values = c(
        "objective" = "#4c78a8",
        "objective_burden" = "#f58518",
        "objective_ploidy" = "#54a24b"
      ),
      guide = "none"
    ) +
    labs(
      title = "Objective Distributions Across Seeds",
      subtitle = paste0("Source: ", basename(extra_results_dir), "/seed_summary.tsv"),
      x = NULL,
      y = "Objective value"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggplot2::ggsave(out_path, p, width = 8, height = 8)

  message("Wrote long table: ", out_table)
  message("Wrote plot: ", out_path)
}

if (sys.nframe() == 0) {
  main()
}
