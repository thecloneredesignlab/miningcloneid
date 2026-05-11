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
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

default_template <- "/Users/4482173/Documents/GitHub/Constant_WGD/oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma}"
default_out_dir <- "/Users/4482173/Documents/GitHub/Constant_WGD/oxygen/results/comp"
default_sigma_caps <- c("0.05", "0.15", "0.3", "0.6")

parse_sigma_caps <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(default_sigma_caps)
  vals <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  vals <- vals[nzchar(vals)]
  if (!length(vals)) stop("sigma_caps must contain at least one comma-separated value.")
  vals
}

build_run_dir <- function(template, sigma_cap) {
  if (!grepl("\\{sigma\\}", template)) {
    stop("run_dir_template must contain the placeholder {sigma}.")
  }
  gsub("\\{sigma\\}", sigma_cap, template)
}

read_seed_summary_one <- function(run_dir, sigma_cap, subgroup_threshold = 44) {
  path <- file.path(run_dir, "extra_results", "seed_summary.tsv")
  if (!file.exists(path)) {
    stop("Missing extra_results/seed_summary.tsv for sigma_burden cap ", sigma_cap, ": ", run_dir)
  }
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required_cols <- c("seed", "objective", "objective_burden", "objective_ploidy", "pred1000_both_gt44")
  missing_cols <- setdiff(required_cols, names(tab))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }
  tab$sigma_cap <- sigma_cap
  tab$sigma_cap_num <- suppressWarnings(as.numeric(sigma_cap))
  tab$objective <- suppressWarnings(as.numeric(tab$objective))
  tab$objective_burden <- suppressWarnings(as.numeric(tab$objective_burden))
  tab$objective_ploidy <- suppressWarnings(as.numeric(tab$objective_ploidy))
  pmiss_col <- "value__p_misseg"
  if (!pmiss_col %in% names(tab)) {
    stop("Missing required p_misseg column in ", path, ": ", pmiss_col)
  }
  tab[[pmiss_col]] <- suppressWarnings(as.numeric(tab[[pmiss_col]]))
  tab$pred1000_both_gt44 <- as.logical(tab$pred1000_both_gt44)
  tab$pred_group <- ifelse(!is.na(tab$pred1000_both_gt44) & tab$pred1000_both_gt44, ">44", "<=44")
  tab$pred_group <- factor(tab$pred_group, levels = c(">44", "<=44"))
  tab$pred_group_note <- paste0("Criterion uses both 2N and 4N 1000d predicted chromosome-number means > ", subgroup_threshold, ".")
  tab
}

make_counts_table <- function(seed_summary) {
  key <- paste(seed_summary$sigma_cap, seed_summary$pred_group, sep = "\r")
  tab <- as.data.frame(table(key), stringsAsFactors = FALSE)
  parts <- strsplit(tab$key, "\r", fixed = TRUE)
  tab$sigma_cap <- vapply(parts, `[[`, character(1), 1L)
  tab$pred_group <- vapply(parts, `[[`, character(1), 2L)
  tab$key <- NULL
  names(tab)[names(tab) == "Freq"] <- "n_seeds"
  tab$pred_group <- factor(tab$pred_group, levels = c(">44", "<=44"))
  tab
}

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

cut_rank_bin_left <- function(x) {
  breaks <- c(0, 10, 25, 50, 100, 200, 300, 400, Inf)
  labels <- c("1-10", "11-25", "26-50", "51-100", "101-200", "201-300", "301-400", "401+")
  cut(x, breaks = breaks, labels = labels, include.lowest = TRUE, right = TRUE)
}

cut_rank_bin_right <- function(x) {
  breaks <- c(0, 5, 10, 25, 50, 100, Inf)
  labels <- c(">44 rank 1-5", ">44 rank 6-10", ">44 rank 11-25", ">44 rank 26-50", ">44 rank 51-100", ">44 rank 101+")
  cut(x, breaks = breaks, labels = labels, include.lowest = TRUE, right = TRUE)
}

build_sankey_flow_table <- function(seed_summary, sigma_levels) {
  rank_col_left <- if ("forest_plot_rank_simple" %in% names(seed_summary)) {
    "forest_plot_rank_simple"
  } else if ("recommend_rank_burden_ploidy_boundary_first" %in% names(seed_summary)) {
    "recommend_rank_burden_ploidy_boundary_first"
  } else {
    stop("Missing forest rank column needed for Sankey plot.")
  }
  rank_col_right <- if ("forest_plot_rank_plus_ploidy_simple" %in% names(seed_summary)) {
    "forest_plot_rank_plus_ploidy_simple"
  } else {
    stop("Missing filtered forest rank column needed for Sankey plot.")
  }

  plot_df <- seed_summary
  plot_df$rank_left <- suppressWarnings(as.numeric(plot_df[[rank_col_left]]))
  plot_df$rank_right <- suppressWarnings(as.numeric(plot_df[[rank_col_right]]))
  plot_df <- plot_df[is.finite(plot_df$rank_left), , drop = FALSE]
  plot_df <- plot_df[order(plot_df$sigma_cap, plot_df$rank_left, plot_df$seed), , drop = FALSE]
  plot_df$left_stratum <- sprintf("overall_%03d", as.integer(plot_df$rank_left))
  plot_df$right_stratum <- ifelse(
    !is.na(plot_df$pred1000_both_gt44) & plot_df$pred1000_both_gt44 & is.finite(plot_df$rank_right),
    sprintf("filtered_%03d", as.integer(plot_df$rank_right)),
    "filtered_excluded"
  )
  plot_df$sigma_cap <- factor(plot_df$sigma_cap, levels = sigma_levels)
  flow_df <- plot_df[, c("seed", "sigma_cap", "pred_group", "rank_left", "rank_right", "left_stratum", "right_stratum"), drop = FALSE]
  flow_df$pred_group <- factor(flow_df$pred_group, levels = c(">44", "<=44"))
  flow_df$n_seeds <- 1
  flow_df$flow_id <- seq_len(nrow(flow_df))
  flow_df
}

save_rank_sankey_plot <- function(seed_summary, out_path, sigma_levels) {
  flow_df <- build_sankey_flow_table(seed_summary = seed_summary, sigma_levels = sigma_levels)
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

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  sigma_caps <- parse_sigma_caps(args$sigma_caps)
  run_dir_template <- args$run_dir_template %||% default_template
  out_dir <- normalizePath(args$out_dir %||% default_out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  seed_tabs <- vector("list", length(sigma_caps))
  for (i in seq_along(sigma_caps)) {
    sigma_cap <- sigma_caps[[i]]
    run_dir <- normalizePath(build_run_dir(run_dir_template, sigma_cap), mustWork = TRUE)
    seed_tabs[[i]] <- read_seed_summary_one(run_dir = run_dir, sigma_cap = sigma_cap)
  }

  seed_summary <- do.call(rbind, seed_tabs)
  sigma_levels <- sigma_caps
  seed_summary$sigma_cap <- factor(seed_summary$sigma_cap, levels = sigma_levels)
  counts_df <- make_counts_table(seed_summary)
  counts_df$sigma_cap <- factor(counts_df$sigma_cap, levels = sigma_levels)

  utils::write.table(
    seed_summary,
    file = file.path(out_dir, "sigma_burden_seed_summary_merged.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    counts_df,
    file = file.path(out_dir, "sigma_burden_pred1000_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    build_sankey_flow_table(seed_summary = seed_summary, sigma_levels = sigma_levels),
    file = file.path(out_dir, "sigma_burden_rank_sankey_flows.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  save_grouped_count_plot(
    counts_df = counts_df,
    out_path = file.path(out_dir, "fig1_pred1000_gt44_counts_by_sigma_burden_cap.pdf"),
    sigma_levels = sigma_levels
  )

  save_violin_box_by_subgroup(
    seed_summary = seed_summary,
    metric_col = "objective",
    metric_label = "Total objective",
    out_path = file.path(out_dir, "fig2_total_objective_by_sigma_burden_cap_and_pred1000_gate.pdf"),
    sigma_levels = sigma_levels
  )
  save_violin_box_by_subgroup(
    seed_summary = seed_summary,
    metric_col = "objective_burden",
    metric_label = "Burden objective",
    out_path = file.path(out_dir, "fig3_burden_objective_by_sigma_burden_cap_and_pred1000_gate.pdf"),
    sigma_levels = sigma_levels
  )
  save_violin_box_by_subgroup(
    seed_summary = seed_summary,
    metric_col = "objective_ploidy",
    metric_label = "Ploidy objective",
    out_path = file.path(out_dir, "fig4_ploidy_objective_by_sigma_burden_cap_and_pred1000_gate.pdf"),
    sigma_levels = sigma_levels
  )
  save_violin_box_by_subgroup(
    seed_summary = seed_summary,
    metric_col = "value__p_misseg",
    metric_label = "p_misseg",
    out_path = file.path(out_dir, "fig9_p_misseg_by_sigma_burden_cap_and_pred1000_gate.pdf"),
    sigma_levels = sigma_levels
  )

  save_violin_box_global(
    seed_summary = seed_summary,
    metric_col = "objective",
    metric_label = "Total objective",
    out_path = file.path(out_dir, "fig5_total_objective_by_sigma_burden_cap.pdf"),
    sigma_levels = sigma_levels
  )
  save_violin_box_global(
    seed_summary = seed_summary,
    metric_col = "objective_burden",
    metric_label = "Burden objective",
    out_path = file.path(out_dir, "fig6_burden_objective_by_sigma_burden_cap.pdf"),
    sigma_levels = sigma_levels
  )
  save_violin_box_global(
    seed_summary = seed_summary,
    metric_col = "objective_ploidy",
    metric_label = "Ploidy objective",
    out_path = file.path(out_dir, "fig7_ploidy_objective_by_sigma_burden_cap.pdf"),
    sigma_levels = sigma_levels
  )
  save_rank_sankey_plot(
    seed_summary = seed_summary,
    out_path = file.path(out_dir, "fig8_rank_sankey_by_sigma_burden_cap.pdf"),
    sigma_levels = sigma_levels
  )

  message("Wrote merged table: ", file.path(out_dir, "sigma_burden_seed_summary_merged.tsv"))
  message("Wrote counts table: ", file.path(out_dir, "sigma_burden_pred1000_counts.tsv"))
  message("Wrote Sankey flow table: ", file.path(out_dir, "sigma_burden_rank_sankey_flows.tsv"))
  message("Wrote figures to: ", out_dir)
}

if (sys.nframe() == 0) {
  main()
}
