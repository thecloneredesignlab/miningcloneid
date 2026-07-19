#!/usr/bin/env Rscript

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

RESULTS_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", "..", "results"), mustWork = FALSE)
default_template <- file.path(
  RESULTS_ROOT,
  "fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma}"
)
default_out_dir <- file.path(RESULTS_ROOT, "comp")
default_sigma_caps <- c("0.05", "0.15", "0.3", "0.6")

parse_sigma_caps <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(default_sigma_caps)
  vals <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  vals <- vals[nzchar(vals)]
  if (!length(vals)) stop("sigma_caps must contain at least one comma-separated value.")
  vals
}

build_run_dir <- o2sd_build_sigma_run_dir

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

  utils::write.table(data.frame(stage = "analysis", file = c("sigma_burden_seed_summary_merged.tsv", "sigma_burden_pred1000_counts.tsv", "sigma_burden_rank_sankey_flows.tsv"), stringsAsFactors = FALSE), file.path(out_dir, "analysis_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  message("Wrote merged table: ", file.path(out_dir, "sigma_burden_seed_summary_merged.tsv"))
  message("Wrote counts table: ", file.path(out_dir, "sigma_burden_pred1000_counts.tsv"))
  message("Wrote Sankey flow table: ", file.path(out_dir, "sigma_burden_rank_sankey_flows.tsv"))
}

if (sys.nframe() == 0) {
  main()
}
