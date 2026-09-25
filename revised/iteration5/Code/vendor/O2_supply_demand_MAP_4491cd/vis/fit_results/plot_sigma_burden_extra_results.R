#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggalluvial))

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
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

 save_grouped_count_plot <- function(counts_df, out_path, sigma_levels) {
  plot_df <- counts_df
  plot_df$sigma_cap <- factor(plot_df$sigma_cap, levels = sigma_levels)
  p <- ggplot(plot_df, aes(x = sigma_cap, y = n_seeds, fill = pred_group)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "grey20") +
    geom_text(
      aes(label = n_seeds),
      position = position_dodge(width = 0.75),
      vjust = -0.25,
      size = 3.4
    ) +
    scale_fill_manual(values = c(">44" = "#1b9e77", "<=44" = "#d95f02"), drop = FALSE) +
    labs(
      title = "1000-Day Chromosome-Number Gate Counts by sigma_burden Upper Bound",
      subtitle = "Grouped by sigma_burden upper bound; gate uses both 2N and 4N 1000d predicted chromosome-number means.",
      x = "sigma_burden upper bound",
      y = "Number of seeds",
      fill = "1000d group"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 8, height = 8)
}

save_violin_box_by_subgroup <- function(seed_summary, metric_col, metric_label, out_path, sigma_levels) {
  plot_df <- seed_summary[is.finite(seed_summary[[metric_col]]), , drop = FALSE]
  plot_df$sigma_cap <- factor(plot_df$sigma_cap, levels = sigma_levels)
  dodge <- position_dodge(width = 0.8)
  p <- ggplot(plot_df, aes(x = sigma_cap, y = .data[[metric_col]], fill = pred_group)) +
    geom_violin(
      position = dodge,
      trim = FALSE,
      alpha = 0.55,
      color = NA
    ) +
    geom_boxplot(
      position = dodge,
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.85,
      linewidth = 0.35
    ) +
    scale_fill_manual(values = c(">44" = "#1b9e77", "<=44" = "#d95f02"), drop = FALSE) +
    labs(
      title = paste0(metric_label, " by sigma_burden Upper Bound and 1000-Day Gate"),
      subtitle = "Within each sigma_burden group, seeds are split by whether both 2N and 4N 1000d predicted chromosome-number means exceed 44.",
      x = "sigma_burden upper bound",
      y = metric_label,
      fill = "1000d group"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 8, height = 8)
}

save_violin_box_global <- function(seed_summary, metric_col, metric_label, out_path, sigma_levels) {
  plot_df <- seed_summary[is.finite(seed_summary[[metric_col]]), , drop = FALSE]
  plot_df$sigma_cap <- factor(plot_df$sigma_cap, levels = sigma_levels)
  p <- ggplot(plot_df, aes(x = sigma_cap, y = .data[[metric_col]], fill = sigma_cap)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title = paste0(metric_label, " by sigma_burden Upper Bound"),
      subtitle = "All seeds pooled within each sigma_burden upper-bound group.",
      x = "sigma_burden upper bound",
      y = metric_label
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 8, height = 8)
}

save_rank_sankey_plot <- function(flow_df, out_path, sigma_levels) {
  axis_levels <- c("Overall forest rank", "Filtered >44 forest rank")
  overall_max <- max(flow_df$rank_left, na.rm = TRUE)
  filtered_max <- suppressWarnings(max(flow_df$rank_right[is.finite(flow_df$rank_right)], na.rm = TRUE))
  if (!is.finite(filtered_max)) filtered_max <- 0
  left_levels <- sprintf("overall_%03d", seq_len(overall_max))
  right_levels <- c(sprintf("filtered_%03d", seq_len(filtered_max)), "filtered_excluded")
  lodes_left <- data.frame(
    sigma_cap = flow_df$sigma_cap,
    pred_group = flow_df$pred_group,
    flow_id = flow_df$flow_id,
    n_seeds = flow_df$n_seeds,
    rank_left = flow_df$rank_left,
    axis = factor(axis_levels[[1]], levels = axis_levels),
    stratum = factor(flow_df$left_stratum, levels = left_levels),
    stringsAsFactors = FALSE
  )
  lodes_right <- data.frame(
    sigma_cap = flow_df$sigma_cap,
    pred_group = flow_df$pred_group,
    flow_id = flow_df$flow_id,
    n_seeds = flow_df$n_seeds,
    rank_left = flow_df$rank_left,
    axis = factor(axis_levels[[2]], levels = axis_levels),
    stratum = factor(flow_df$right_stratum, levels = right_levels),
    stringsAsFactors = FALSE
  )
  lodes_df <- rbind(lodes_left, lodes_right)
  lodes_df$sigma_cap <- factor(lodes_df$sigma_cap, levels = sigma_levels)
  lodes_df$pred_group <- factor(lodes_df$pred_group, levels = c(">44", "<=44"))
  p <- ggplot(
    lodes_df,
    aes(x = axis, stratum = stratum, alluvium = flow_id, y = n_seeds, order = rank_left)
  ) +
    ggalluvial::geom_flow(
      aes(fill = pred_group),
      width = 0.2,
      alpha = 0.75,
      knot.pos = 0.42
    ) +
    ggalluvial::geom_stratum(
      width = 0.2,
      fill = "grey95",
      color = "grey45",
      linewidth = 0.35
    ) +
    scale_x_discrete(limits = axis_levels, expand = c(0.08, 0.08)) +
    scale_fill_manual(values = c(">44" = "#1b9e77", "<=44" = "#d95f02"), drop = FALSE) +
    labs(
      title = "Rank Flow from Overall Forest Ordering to Filtered >44 Ordering",
      subtitle = "Each seed is one flow. Left uses the natural overall forest rank ordering; right uses the filtered >44 forest rank, with <=44 seeds flowing to the exclusion node.",
      x = NULL,
      y = "Number of seeds",
      fill = "1000d group"
    ) +
    facet_wrap(~ sigma_cap, nrow = 1, ncol = 4, scales = "free_y") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey60"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom"
    )
  ggplot2::ggsave(out_path, p, width = 22, height = 11)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  analysis_dir <- normalizePath(argv$analysis_dir %||% argv$out_dir, mustWork = TRUE)
  if (!file.exists(file.path(analysis_dir, "analysis_manifest.tsv"))) stop("Missing sigma-burden analysis manifest.", call. = FALSE)
  viz_dir <- normalizePath(argv$viz_dir %||% argv$out_dir %||% analysis_dir, mustWork = FALSE)
  required <- file.path(analysis_dir, c(
    "sigma_burden_seed_summary_merged.tsv",
    "sigma_burden_pred1000_counts.tsv",
    "sigma_burden_rank_sankey_flows.tsv"
  ))
  if (any(!file.exists(required))) stop("Missing sigma-burden analysis tables: ", paste(required[!file.exists(required)], collapse = ", "), call. = FALSE)
  seed_summary <- utils::read.delim(required[[1]], check.names = FALSE, stringsAsFactors = FALSE)
  counts_df <- utils::read.delim(required[[2]], check.names = FALSE, stringsAsFactors = FALSE)
  flow_df <- utils::read.delim(required[[3]], check.names = FALSE, stringsAsFactors = FALSE)
  sigma_levels <- unique(as.character(seed_summary$sigma_cap[order(seed_summary$sigma_cap_num)]))
  seed_summary$sigma_cap <- factor(seed_summary$sigma_cap, levels = sigma_levels)
  seed_summary$pred_group <- factor(seed_summary$pred_group, levels = c(">44", "<=44"))
  counts_df$sigma_cap <- factor(counts_df$sigma_cap, levels = sigma_levels)
  counts_df$pred_group <- factor(counts_df$pred_group, levels = c(">44", "<=44"))
  flow_df$sigma_cap <- factor(flow_df$sigma_cap, levels = sigma_levels)
  flow_df$pred_group <- factor(flow_df$pred_group, levels = c(">44", "<=44"))
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  specs <- list(
    c("fig1_pred1000_gt44_counts_by_sigma_burden_cap.pdf", "counts"),
    c("fig2_total_objective_by_sigma_burden_cap_and_pred1000_gate.pdf", "subgroup", "objective", "Total objective"),
    c("fig3_burden_objective_by_sigma_burden_cap_and_pred1000_gate.pdf", "subgroup", "objective_burden", "Burden objective"),
    c("fig4_ploidy_objective_by_sigma_burden_cap_and_pred1000_gate.pdf", "subgroup", "objective_ploidy", "Ploidy objective"),
    c("fig9_p_misseg_by_sigma_burden_cap_and_pred1000_gate.pdf", "subgroup", "value__p_misseg", "p_misseg"),
    c("fig5_total_objective_by_sigma_burden_cap.pdf", "global", "objective", "Total objective"),
    c("fig6_burden_objective_by_sigma_burden_cap.pdf", "global", "objective_burden", "Burden objective"),
    c("fig7_ploidy_objective_by_sigma_burden_cap.pdf", "global", "objective_ploidy", "Ploidy objective"),
    c("fig8_rank_sankey_by_sigma_burden_cap.pdf", "sankey")
  )
  outputs <- character()
  for (spec in specs) {
    path <- file.path(viz_dir, spec[[1]])
    if (spec[[2]] == "counts") save_grouped_count_plot(counts_df, path, sigma_levels)
    if (spec[[2]] == "subgroup") save_violin_box_by_subgroup(seed_summary, spec[[3]], spec[[4]], path, sigma_levels)
    if (spec[[2]] == "global") save_violin_box_global(seed_summary, spec[[3]], spec[[4]], path, sigma_levels)
    if (spec[[2]] == "sankey") save_rank_sankey_plot(flow_df, path, sigma_levels)
    if (file.exists(path)) outputs <- c(outputs, path)
  }
  utils::write.table(data.frame(stage = "visualization", file = basename(outputs), stringsAsFactors = FALSE), file.path(viz_dir, "visualization_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote sigma-burden visualizations to: ", viz_dir)
  invisible(viz_dir)
}

if (sys.nframe() == 0L) main()
