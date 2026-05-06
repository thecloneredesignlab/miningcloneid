#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

.o2g_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2g_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
OXYGEN_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", ".."), mustWork = FALSE)
HELPER_DIR <- normalizePath(file.path(OXYGEN_ROOT, "code", "in-vitro-utils"), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
sys.source(file.path(HELPER_DIR, "plotting.R"), envir = environment(), chdir = TRUE)
rm(.o2g_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

read_tsv_optional <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

num <- function(x) suppressWarnings(as.numeric(x))

finite_rows <- function(df, cols) {
  ok <- rep(TRUE, nrow(df))
  for (col in cols) {
    ok <- ok & is.finite(df[[col]])
  }
  df[ok, , drop = FALSE]
}

summary_value <- function(summary_df, metric) {
  if (is.null(summary_df) || !all(c("metric", "value") %in% names(summary_df))) return(NA_character_)
  idx <- match(metric, summary_df$metric)
  if (is.na(idx)) NA_character_ else summary_df$value[[idx]]
}

save_plot_pair <- function(plot, out_dir, basename, width = 9, height = 5.5) {
  pdf_path <- file.path(out_dir, paste0(basename, ".pdf"))
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in")
  ggplot2::ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 180)
  unlink(file.path(out_dir, paste0(basename, ".svg")), force = TRUE)
  invisible(c(pdf = pdf_path, png = png_path))
}

save_existing_plot_png <- function(plot, out_dir, basename, width = 10, height = 7) {
  if (is.null(plot)) return(FALSE)
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 180)
  unlink(file.path(out_dir, paste0(basename, ".svg")), force = TRUE)
  TRUE
}

ensure_invitro_plot_columns <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  n <- nrow(df)
  if (!"lineage_passage_index" %in% names(df)) {
    df$lineage_passage_index <- if ("passage_index" %in% names(df)) df$passage_index else seq_len(n)
  }
  if (!"lineage_terminal_key" %in% names(df)) {
    df$lineage_terminal_key <- if ("segment_id" %in% names(df)) as.character(df$segment_id) else as.character(seq_len(n))
  }
  if (!"lineage_label" %in% names(df)) {
    df$lineage_label <- if ("cohort" %in% names(df)) as.character(df$cohort) else "lineage"
  }
  if (!"sample_name" %in% names(df)) {
    df$sample_name <- if ("passage_id" %in% names(df)) as.character(df$passage_id) else df$lineage_terminal_key
  }
  df
}

plot_remote_lineage_growth <- function(lineage_df, out_dir) {
  required <- c("predicted_growth_rate", "observed_growth", "passage_index", "cohort")
  if (is.null(lineage_df) || !all(required %in% names(lineage_df))) return(invisible(FALSE))
  df <- ensure_invitro_plot_columns(lineage_df)
  p <- ivt_plot_lineage_growth(df, primary_label = "Best fit")
  save_plot_pair(p, out_dir, "invitro_lineage_growth", width = 10, height = 6.8)
  invisible(TRUE)
}

plot_remote_lineage_ploidy <- function(lineage_df, quantile_df, observed_kary_df, out_dir) {
  if (is.null(lineage_df) || is.null(quantile_df) || is.null(observed_kary_df)) return(invisible(FALSE))
  if (!all(c("predicted_quantile_kary_N", "quantile_prob") %in% names(quantile_df))) return(invisible(FALSE))
  if (!"observed_kary_N" %in% names(observed_kary_df)) return(invisible(FALSE))
  p <- ivt_plot_lineage_ploidy(
    ensure_invitro_plot_columns(lineage_df),
    quantile_df = ensure_invitro_plot_columns(quantile_df),
    observed_kary_df = ensure_invitro_plot_columns(observed_kary_df),
    primary_label = "Best fit",
    quantile_alpha = 0.5
  )
  save_plot_pair(p, out_dir, "invitro_lineage_ploidy", width = 10, height = 6.8)
  invisible(TRUE)
}

plot_remote_flow_density <- function(flow_df, out_dir) {
  if (is.null(flow_df) || !all(c("ploidy", "density", "series", "cohort") %in% names(flow_df))) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No matched flow-density observations available")
    save_plot_pair(p, out_dir, "invitro_flow_density", width = 10, height = 6.8)
    return(invisible(TRUE))
  }
  p <- ivt_plot_lineage_flow_density(
    ensure_invitro_plot_columns(flow_df),
    max_facets = 20L
  )
  save_plot_pair(p, out_dir, "invitro_flow_density", width = 10, height = 6.8)
  invisible(TRUE)
}

plot_remote_distribution_heatmap <- function(dist_df, out_dir) {
  if (is.null(dist_df) || !all(c("N", "fraction", "passage_index", "cohort") %in% names(dist_df))) {
    return(invisible(FALSE))
  }
  p <- ivt_plot_distribution_heatmap(ensure_invitro_plot_columns(dist_df), max_N = 110)
  save_plot_pair(p, out_dir, "invitro_distribution_heatmap", width = 10, height = 6.5)
  invisible(TRUE)
}

plot_remote_loglik_by_passage <- function(df,
                                          out_dir,
                                          basename,
                                          value_col,
                                          title,
                                          y_label) {
  if (is.null(df) || !all(c("passage_index", value_col, "cohort", "segment_id") %in% names(df))) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste0("No ", tolower(gsub("In vitro ", "", title)), " rows"))
    save_plot_pair(p, out_dir, basename, width = 9.5, height = 5)
    return(invisible(TRUE))
  }
  plot_df <- ensure_invitro_plot_columns(df)
  plot_df$passage_index <- num(plot_df$passage_index)
  plot_df[[value_col]] <- num(plot_df[[value_col]])
  plot_df <- finite_rows(plot_df, c("passage_index", value_col))
  if (nrow(plot_df) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste0("No ", tolower(gsub("In vitro ", "", title)), " rows"))
    save_plot_pair(p, out_dir, basename, width = 9.5, height = 5)
    return(invisible(TRUE))
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(passage_index, .data[[value_col]], color = cohort)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(cohort, segment_id)), alpha = 0.35) +
    ggplot2::labs(
      title = title,
      x = "Passage index",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 12)
  save_plot_pair(p, out_dir, basename, width = 9.5, height = 5)
  invisible(TRUE)
}

theme_invitro <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey75"),
      legend.position = "bottom"
    )
}

plot_objective_components <- function(summary_df, out_dir) {
  if (is.null(summary_df)) return(invisible(FALSE))
  value_for <- function(metric) suppressWarnings(as.numeric(summary_value(summary_df, metric)))
  metrics <- data.frame(
    component = c(
      "Total objective",
      "- Growth logLik sum",
      "- Karyotype logLik sum",
      "- Flow logLik sum"
    ),
    value = c(
      value_for("objective_total"),
      -value_for("growth_loglik_sum"),
      -value_for("ploidy_loglik_sum"),
      -value_for("flow_loglik_sum")
    ),
    stringsAsFactors = FALSE
  )
  metrics <- metrics[is.finite(metrics$value), , drop = FALSE]
  if (nrow(metrics) == 0L) return(invisible(FALSE))
  metrics$component <- factor(metrics$component, levels = rev(metrics$component))
  p <- ggplot2::ggplot(metrics, ggplot2::aes(x = component, y = value)) +
    ggplot2::geom_col(fill = "#4B6F8A", width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "In Vitro Objective Components",
      x = NULL,
      y = "Reported objective-scale value"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_objective_components", width = 8.5, height = 4.2)
  invisible(TRUE)
}

plot_growth_rate_fit <- function(growth_df, out_dir) {
  if (is.null(growth_df) || !all(c("observed_growth", "predicted_growth_rate") %in% names(growth_df))) {
    return(invisible(FALSE))
  }
  df <- growth_df
  df$observed <- num(df$observed_growth)
  df$predicted <- num(df$predicted_growth_rate)
  df <- finite_rows(df, c("observed", "predicted"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  axis_range <- range(c(df$observed, df$predicted), finite = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = predicted, colour = cohort)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_point(alpha = 0.8, size = 2.2) +
    ggplot2::coord_equal(xlim = axis_range, ylim = axis_range) +
    ggplot2::labs(
      title = "In Vitro Growth-Rate Fit",
      x = "Observed growth rate",
      y = "Predicted growth rate",
      colour = "Cohort"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_growth_rate_fit", width = 7.2, height = 6)
  invisible(TRUE)
}

plot_growth_count_fit <- function(growth_df, out_dir) {
  if (is.null(growth_df) || !all(c("target_live_cells", "predicted_live_cells") %in% names(growth_df))) {
    return(invisible(FALSE))
  }
  df <- growth_df
  df$observed_cells <- num(df$target_live_cells)
  df$predicted_cells <- num(df$predicted_live_cells)
  df <- finite_rows(df, c("observed_cells", "predicted_cells"))
  df <- df[df$observed_cells > 0 & df$predicted_cells > 0, , drop = FALSE]
  if (nrow(df) == 0L) return(invisible(FALSE))
  df$log10_observed_cells <- log10(df$observed_cells)
  df$log10_predicted_cells <- log10(df$predicted_cells)
  axis_range <- range(c(df$log10_observed_cells, df$log10_predicted_cells), finite = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log10_observed_cells, y = log10_predicted_cells, colour = cohort)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_point(alpha = 0.75, size = 2) +
    ggplot2::coord_equal(xlim = axis_range, ylim = axis_range) +
    ggplot2::labs(
      title = "In Vitro Live-Cell Count Fit",
      x = "Observed live cells (log10)",
      y = "Predicted live cells (log10)",
      colour = "Cohort"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_live_count_fit", width = 7.2, height = 6)
  invisible(TRUE)
}

plot_ploidy_mean_fit <- function(lineage_df, out_dir) {
  if (is.null(lineage_df) || !all(c("observed_mean_kary_N", "predicted_mean_kary_N") %in% names(lineage_df))) {
    return(invisible(FALSE))
  }
  df <- lineage_df
  df$observed <- num(df$observed_mean_kary_N)
  df$predicted <- num(df$predicted_mean_kary_N)
  df <- finite_rows(df, c("observed", "predicted"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  axis_range <- range(c(df$observed, df$predicted), finite = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = predicted, colour = cohort)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_point(alpha = 0.85, size = 2.4) +
    ggplot2::coord_equal(xlim = axis_range, ylim = axis_range) +
    ggplot2::labs(
      title = "In Vitro Mean Karyotype Fit",
      x = "Observed mean chromosome number (N)",
      y = "Predicted mean chromosome number (N)",
      colour = "Cohort"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_mean_karyotype_fit", width = 7.2, height = 6)
  invisible(TRUE)
}

plot_ploidy_loglik <- function(ploidy_df, out_dir) {
  if (is.null(ploidy_df) || !all(c("passage_index", "mean_loglik") %in% names(ploidy_df))) {
    return(invisible(FALSE))
  }
  df <- ploidy_df
  df$passage <- num(df$passage_index)
  df$mean_loglik_num <- num(df$mean_loglik)
  df <- finite_rows(df, c("passage", "mean_loglik_num"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = passage, y = mean_loglik_num, colour = cohort)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::labs(
      title = "In Vitro Karyotype Log-Likelihood by Passage",
      x = "Passage index",
      y = "Mean log-likelihood per observed cell",
      colour = "Cohort"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_karyotype_loglik", width = 8.5, height = 5.2)
  invisible(TRUE)
}

plot_distribution_quantiles <- function(quantile_df, out_dir) {
  required <- c("cohort", "passage_index", "quantile_prob", "predicted_quantile_kary_N")
  if (is.null(quantile_df) || !all(required %in% names(quantile_df))) return(invisible(FALSE))
  df <- quantile_df
  df$passage <- num(df$passage_index)
  df$quantile <- num(df$quantile_prob)
  df$predicted <- num(df$predicted_quantile_kary_N)
  df <- finite_rows(df, c("passage", "quantile", "predicted"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  agg <- stats::aggregate(
    predicted ~ cohort + passage + quantile,
    data = df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  agg$quantile <- factor(agg$quantile)
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = passage, y = predicted, colour = quantile)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~cohort, scales = "free_x") +
    ggplot2::labs(
      title = "Predicted In Vitro Karyotype Quantiles",
      x = "Passage index",
      y = "Predicted chromosome number (N)",
      colour = "Quantile"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_karyotype_quantiles", width = 9.5, height = 5.4)
  invisible(TRUE)
}

plot_distribution_heatmap <- function(dist_df, out_dir) {
  required <- c("cohort", "passage_index", "N", "fraction")
  if (is.null(dist_df) || !all(required %in% names(dist_df))) return(invisible(FALSE))
  df <- dist_df
  df$passage <- num(df$passage_index)
  df$kary_N <- num(df$N)
  df$fraction_num <- num(df$fraction)
  df <- finite_rows(df, c("passage", "kary_N", "fraction_num"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  agg <- stats::aggregate(
    fraction_num ~ cohort + passage + kary_N,
    data = df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = passage, y = kary_N, fill = fraction_num)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~cohort, scales = "free_x") +
    ggplot2::scale_fill_viridis_c(option = "C", trans = "sqrt") +
    ggplot2::labs(
      title = "Predicted In Vitro Karyotype Distribution",
      x = "Passage index",
      y = "Chromosome number (N)",
      fill = "Mean fraction"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_karyotype_distribution", width = 9.5, height = 5.4)
  invisible(TRUE)
}

plot_daily_live_cells <- function(daily_df, out_dir) {
  required <- c("cohort", "segment_id", "day", "live_cells")
  if (is.null(daily_df) || !all(required %in% names(daily_df))) return(invisible(FALSE))
  df <- daily_df
  df$day_num <- num(df$day)
  df$live_cells_num <- num(df$live_cells)
  df$oxygen_pct_num <- num(df$oxygen_pct)
  df <- finite_rows(df, c("day_num", "live_cells_num"))
  df <- df[df$live_cells_num > 0, , drop = FALSE]
  if (nrow(df) == 0L) return(invisible(FALSE))
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = day_num,
      y = live_cells_num,
      group = segment_id,
      colour = oxygen_pct_num
    )
  ) +
    ggplot2::geom_line(alpha = 0.35, linewidth = 0.45) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_viridis_c(option = "B", na.value = "grey45") +
    ggplot2::facet_wrap(~cohort, scales = "free_x") +
    ggplot2::labs(
      title = "Predicted In Vitro Live-Cell Time Courses",
      x = "Day within passage",
      y = "Live cells (log10)",
      colour = "O2 (%)"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_daily_live_cells", width = 10, height = 5.8)
  invisible(TRUE)
}

plot_flow_overlay <- function(flow_df, out_dir) {
  if (is.null(flow_df) || nrow(flow_df) == 0L) return(invisible(FALSE))
  names_lower <- tolower(names(flow_df))
  x_col <- names(flow_df)[match(TRUE, names_lower %in% c("ploidy", "ploidy_n", "kary_n", "n", "x"))]
  pred_col <- names(flow_df)[match(TRUE, names_lower %in% c("predicted_density", "pred_density", "model_density"))]
  obs_col <- names(flow_df)[match(TRUE, names_lower %in% c("observed_density", "obs_density", "density"))]
  if (is.na(x_col) || is.na(pred_col) || is.na(obs_col)) return(invisible(FALSE))
  df <- flow_df
  df$x <- num(df[[x_col]])
  df$predicted <- num(df[[pred_col]])
  df$observed <- num(df[[obs_col]])
  df <- finite_rows(df, c("x", "predicted", "observed"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  long <- rbind(
    data.frame(df[, setdiff(names(df), c("predicted", "observed")), drop = FALSE], density = df$predicted, series = "Predicted"),
    data.frame(df[, setdiff(names(df), c("predicted", "observed")), drop = FALSE], density = df$observed, series = "Observed")
  )
  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = density, colour = series)) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
    ggplot2::facet_wrap(~cohort, scales = "free_y") +
    ggplot2::labs(
      title = "In Vitro Flow-Density Overlay",
      x = x_col,
      y = "Density",
      colour = NULL
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_flow_overlay", width = 9.5, height = 5.4)
  invisible(TRUE)
}

plot_functional_response_curves_if_available <- function(fit_dir, out_dir) {
  fit_result_path <- file.path(fit_dir, "fit_result.rds")
  if (!file.exists(fit_result_path)) return(FALSE)
  fit_result <- tryCatch(readRDS(fit_result_path), error = function(e) NULL)
  if (is.null(fit_result) || is.null(fit_result$cfg) || is.null(fit_result$best_params)) {
    return(FALSE)
  }
  invivo_viz_script <- file.path(WORKFLOW_ROOT, "vis", "viz_invivo_model_O2G_supply_demand_MAP_results.R")
  if (!file.exists(invivo_viz_script)) return(FALSE)
  invivo_env <- new.env(parent = globalenv())
  tryCatch(
    {
      sys.source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), envir = invivo_env, chdir = TRUE)
      Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = file.path(WORKFLOW_ROOT, "model"))
      sys.source(file.path(WORKFLOW_ROOT, "model", "model_O2G_supply_demand_MAP.R"), envir = invivo_env, chdir = TRUE)
      sys.source(invivo_viz_script, envir = invivo_env, chdir = TRUE)
    },
    error = function(e) {
      warning("Could not source in vivo viz functions for in vitro functional curves: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (!exists("plot_functional_response_curves", envir = invivo_env, inherits = FALSE)) {
    return(FALSE)
  }
  cfg <- fit_result$cfg
  run_params <- fit_result$best_params
  cfg$glucose <- FALSE
  cfg$glucose_dynamic <- FALSE
  cfg$glucose_stress_mode <- "off"
  run_params$glucose <- FALSE
  run_params$glucose_dynamic <- FALSE
  run_params$glucose_stress_mode <- "off"
  functional_plots <- tryCatch(
    invivo_env$plot_functional_response_curves(
      run_params = run_params,
      cfg = cfg,
      out_dir = out_dir
    ),
    error = function(e) {
      warning("Could not generate in vitro functional response curves: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(functional_plots)) return(FALSE)
  save_existing_plot_png(functional_plots$p_msr_o2, out_dir, "oxygen_vs_missegregation_rate")
  save_existing_plot_png(functional_plots$p_msr_o2_multi, out_dir, "oxygen_vs_missegregation_rate_multi_ploidy")
  save_existing_plot_png(functional_plots$p_msr_death, out_dir, "ms_rate_vs_death_rate")
  save_existing_plot_png(functional_plots$p_msr_buffer_death, out_dir, "ms_rate_vs_buffer_death_rate")
  save_existing_plot_png(functional_plots$p_msr_buffer_death_per_division, out_dir, "ms_rate_vs_buffer_death_per_division")
  save_existing_plot_png(functional_plots$p_msr_buffer_death_per_division, out_dir, "ms_rate_vs_nonviable_daughter_fraction")
  save_existing_plot_png(functional_plots$p_msr_nonviable_division_prob, out_dir, "ms_rate_vs_nonviable_division_probability")
  save_existing_plot_png(functional_plots$p_prolif, out_dir, "oxygen_vs_proliferation_rate")
  save_existing_plot_png(functional_plots$p_death, out_dir, "oxygen_vs_death_rate")
  save_existing_plot_png(functional_plots$p_net, out_dir, "oxygen_vs_net_growth_rate")
  save_existing_plot_png(functional_plots$p_viability, out_dir, "ploidy_vs_viability_after_ms")
  save_existing_plot_png(functional_plots$p_ploidy_prolif_o2, out_dir, "ploidy_vs_proliferation_rate_by_o2")
  save_existing_plot_png(functional_plots$p_ploidy_death_o2, out_dir, "ploidy_vs_death_rate_by_o2")
  TRUE
}

write_manifest <- function(out_dir, fit_dir, generated) {
  manifest <- data.frame(
    key = c("fit_dir", "generated_at", names(generated)),
    value = c(
      normalizePath(fit_dir, mustWork = FALSE),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      as.character(unlist(generated, use.names = FALSE))
    ),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    manifest,
    file = file.path(out_dir, "invitro_viz_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir %||% stop(
    "Usage: viz_invitro_model_O2G_supply_demand_MAP_results.R --fit_dir=/abs/path/to/seed_dir",
    call. = FALSE
  )
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  out_dir <- argv$out_dir %||% file.path(fit_dir, "viz")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(out_dir, pattern = "^invitro_.*\\.(pdf|png|svg)$", full.names = TRUE), force = TRUE)

  summary_df <- read_tsv_optional(file.path(fit_dir, "fit_summary.tsv"))
  growth_df <- read_tsv_optional(file.path(fit_dir, "invitro_growth_loglik.tsv"))
  lineage_df <- read_tsv_optional(file.path(fit_dir, "invitro_lineage_summary.tsv"))
  ploidy_df <- read_tsv_optional(file.path(fit_dir, "invitro_ploidy_loglik.tsv"))
  dist_df <- read_tsv_optional(file.path(fit_dir, "invitro_distribution_summary.tsv"))
  quantile_df <- read_tsv_optional(file.path(fit_dir, "invitro_distribution_quantiles.tsv"))
  daily_df <- read_tsv_optional(file.path(fit_dir, "invitro_daily_counts.tsv"))
  flow_df <- read_tsv_optional(file.path(fit_dir, "invitro_flow_overlay.tsv"))
  observed_kary_df <- read_tsv_optional(file.path(fit_dir, "invitro_observed_kary.tsv"))

  generated <- list(
    lineage_growth = plot_remote_lineage_growth(lineage_df, out_dir),
    lineage_ploidy = plot_remote_lineage_ploidy(lineage_df, quantile_df, observed_kary_df, out_dir),
    flow_density = plot_remote_flow_density(flow_df, out_dir),
    distribution_heatmap = plot_remote_distribution_heatmap(dist_df, out_dir),
    growth_loglik_by_passage = plot_remote_loglik_by_passage(
      growth_df,
      out_dir,
      "invitro_growth_loglik_by_passage",
      "loglik",
      "In vitro growth log likelihood by passage",
      "Log likelihood"
    ),
    ploidy_loglik_by_passage = plot_remote_loglik_by_passage(
      ploidy_df,
      out_dir,
      "invitro_ploidy_loglik_by_passage",
      "mean_loglik",
      "In vitro chromosome-count log likelihood by passage",
      "Mean log likelihood"
    ),
    flow_loglik_by_passage = plot_remote_loglik_by_passage(
      flow_df,
      out_dir,
      "invitro_flow_loglik_by_passage",
      "mean_loglik",
      "In vitro flow-density log likelihood by passage",
      "Mean log likelihood"
    )
  )
  write_manifest(out_dir, fit_dir, generated)
  message("In vitro viz written to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(normalizePath(out_dir, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  main()
}
