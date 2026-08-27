#!/usr/bin/env Rscript

# Visualization-only consumer for the live-effective-p_ms comparison.
# All plotted values must already exist in the analysis plot table.
#
# Usage:
#   Rscript oxygen/code/O2_supply_demand_MAP/vis/profile_likelihood/plot_sigma_burden_live_effective_pms.R \
#     --analysis_dir=/path/to/live_effective_pms_analysis \
#     --out_dir=/path/to/figures

suppressPackageStartupMessages(library(ggplot2))

.o2sd_live_plot_bootstrap_script_dir <- local({
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
    self <- frame_files[
      basename(frame_files) == "plot_sigma_burden_live_effective_pms.R"
    ]
    if (length(self) > 0L) return(dirname(self[[1L]]))
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2sd_live_plot_bootstrap_script_dir, mustWork = FALSE)
locate_workflow_root <- function(starts) {
  for (start in unique(starts)) {
    current <- normalizePath(start, mustWork = FALSE)
    for (depth in 0:10) {
      candidates <- c(
        current,
        file.path(current, "oxygen", "code", "O2_supply_demand_MAP"),
        file.path(current, "code", "O2_supply_demand_MAP")
      )
      hits <- candidates[file.exists(file.path(
        candidates,
        "util",
        "o2_supply_demand_map_shared.R"
      ))]
      if (length(hits)) return(normalizePath(hits[[1L]], mustWork = TRUE))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate the O2_supply_demand_MAP workflow root.")
}
WORKFLOW_ROOT <- locate_workflow_root(c(SCRIPT_DIR, getwd()))
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_live_plot_bootstrap_script_dir, locate_workflow_root)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
resolve_path_value <- o2sd_resolve_path
read_required_tsv <- o2sd_read_tsv
write_tsv <- o2sd_write_tsv

validate_analysis_contract <- function(analysis_dir) {
  manifest_path <- file.path(
    analysis_dir,
    "live_effective_pms_comparison_manifest.tsv"
  )
  manifest <- read_required_tsv(manifest_path)
  if (!all(c("key", "value") %in% names(manifest))) {
    stop("Invalid analysis manifest schema: ", manifest_path)
  }
  status <- manifest$value[manifest$key == "status"]
  if (!length(status) || !identical(as.character(status[[1]]), "complete")) {
    stop("Analysis manifest is not marked complete: ", manifest_path)
  }
  manifest_path
}

build_live_effective_pms_plot <- function(plot_df, sigma_levels = NULL) {
  required_cols <- c("sigma_cap", "seed", "estimate_type", "value")
  if (!all(required_cols %in% names(plot_df))) {
    stop(
      "Plot table is missing required columns: ",
      paste(setdiff(required_cols, names(plot_df)), collapse = ", ")
    )
  }
  plot_df$value <- suppressWarnings(as.numeric(plot_df$value))
  plot_df <- plot_df[is.finite(plot_df$value), , drop = FALSE]
  if (!nrow(plot_df)) stop("Plot table has no finite values.")
  if (is.null(sigma_levels) || !length(sigma_levels)) {
    sigma_num <- suppressWarnings(as.numeric(unique(plot_df$sigma_cap)))
    sigma_levels <- unique(as.character(plot_df$sigma_cap))
    sigma_levels <- sigma_levels[order(sigma_num, sigma_levels, na.last = TRUE)]
  }
  plot_df$sigma_cap <- factor(plot_df$sigma_cap, levels = sigma_levels)
  plot_df$estimate_type <- factor(
    plot_df$estimate_type,
    levels = c("p_misseg parameter", "live-cell effective p_misseg")
  )

  dodge <- position_dodge(width = 0.8)
  group_sizes <- table(plot_df$sigma_cap, plot_df$estimate_type)
  enough_for_violin <- all(group_sizes[group_sizes > 0L] >= 2L)
  plot_object <- ggplot(
    plot_df,
    aes(x = sigma_cap, y = value, fill = estimate_type)
  )
  if (enough_for_violin) {
    plot_object <- plot_object +
      geom_violin(position = dodge, trim = FALSE, alpha = 0.55, color = NA)
  }
  plot_object <- plot_object +
    geom_boxplot(
      position = dodge,
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.9,
      linewidth = 0.35
    )
  if (!enough_for_violin) {
    plot_object <- plot_object +
      geom_point(
        aes(color = estimate_type),
        position = dodge,
        size = 2,
        show.legend = FALSE
      )
  }
  plot_object +
    scale_fill_manual(
      values = c(
        "p_misseg parameter" = "#4c78a8",
        "live-cell effective p_misseg" = "#f58518"
      ),
      drop = FALSE
    ) +
    labs(
      title = "p_misseg vs Live-Cell Effective p_misseg by sigma_burden Upper Bound",
      subtitle = "Each sigma_burden group includes all seeds. The live-cell value is the all-sample-days live-weighted effective p_ms mean from estimate_live_effective_pms.R.",
      x = "sigma_burden upper bound",
      y = "p_misseg estimate",
      fill = "Estimate type"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  analysis_dir <- resolve_path_value(argv$analysis_dir, getwd())
  if (is.null(analysis_dir)) stop("--analysis_dir is required.")
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  analysis_manifest <- validate_analysis_contract(analysis_dir)

  plot_table_path <- resolve_path_value(argv$plot_table, analysis_dir) %||%
    file.path(analysis_dir, "sigma_burden_p_misseg_vs_live_cell_plot.tsv")
  plot_df <- read_required_tsv(plot_table_path)
  sigma_levels <- NULL
  if (!is.null(argv$sigma_caps) && nzchar(trimws(as.character(argv$sigma_caps)))) {
    sigma_levels <- trimws(
      unlist(strsplit(as.character(argv$sigma_caps), ",", fixed = TRUE))
    )
    sigma_levels <- sigma_levels[nzchar(sigma_levels)]
  }

  out_dir <- resolve_path_value(argv$out_dir, getwd()) %||%
    file.path(analysis_dir, "figures")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(
    out_dir,
    "sigma_burden_p_misseg_vs_live_cell_violin.pdf"
  )
  plot_object <- build_live_effective_pms_plot(plot_df, sigma_levels)
  ggplot2::ggsave(out_path, plot_object, width = 8, height = 8)

  manifest <- data.frame(
    key = c(
      "schema_version",
      "status",
      "analysis_manifest",
      "plot_table",
      "figure",
      "generated_at"
    ),
    value = c(
      "o2sd-live-effective-pms-visualization-v1",
      "complete",
      normalizePath(analysis_manifest, mustWork = TRUE),
      normalizePath(plot_table_path, mustWork = TRUE),
      normalizePath(out_path, mustWork = TRUE),
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(
    out_dir,
    "live_effective_pms_visualization_manifest.tsv"
  )
  write_tsv(manifest, manifest_path)

  message("Wrote figure: ", out_path)
  message("Wrote visualization manifest: ", manifest_path)
}

if (sys.nframe() == 0) {
  main()
}
