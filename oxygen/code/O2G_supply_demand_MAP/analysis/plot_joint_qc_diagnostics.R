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
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

format_spearman_label <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return("Spearman rho = NA")
  }
  rho <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  if (!is.finite(rho)) return("Spearman rho = NA")
  sprintf("Spearman rho = %.2f", rho)
}

add_linear_trend <- function(p, plot_df, x_col = "qc", group_col = NULL) {
  if (is.null(group_col) || !(group_col %in% names(plot_df))) {
    trend_df <- plot_df
  } else {
    trend_parts <- lapply(split(plot_df, plot_df[[group_col]], drop = TRUE), function(part) {
      x <- suppressWarnings(as.numeric(part[[x_col]]))
      y <- suppressWarnings(as.numeric(part$value))
      keep <- is.finite(x) & is.finite(y)
      if (sum(keep) < 3L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
        return(NULL)
      }
      part[keep, , drop = FALSE]
    })
    trend_df <- do.call(rbind, Filter(Negate(is.null), trend_parts))
  }
  if (is.null(trend_df) || !nrow(trend_df)) return(p)
  p + geom_smooth(
    data = trend_df,
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.7,
    alpha = 0.9
  )
}

plot_joint_qc_objective_diagnostics <- function(summary_df, out_path, run_label) {
  required <- c("seed", "value__qc", "objective", "objective_invivo", "objective_invitro")
  missing <- setdiff(required, names(summary_df))
  if (length(missing)) {
    stop("seed_summary is missing required columns for QC-objective diagnostics: ", paste(missing, collapse = ", "))
  }

  base_df <- data.frame(
    seed = as.character(summary_df$seed),
    qc = suppressWarnings(as.numeric(summary_df$value__qc)),
    total_objective = suppressWarnings(as.numeric(summary_df$objective)),
    invivo_objective = suppressWarnings(as.numeric(summary_df$objective_invivo)),
    invitro_objective = suppressWarnings(as.numeric(summary_df$objective_invitro)),
    stringsAsFactors = FALSE
  )
  comp_long <- rbind(
    data.frame(seed = base_df$seed, qc = base_df$qc, component = "Total objective", value = base_df$total_objective, stringsAsFactors = FALSE),
    data.frame(seed = base_df$seed, qc = base_df$qc, component = "In vivo objective", value = base_df$invivo_objective, stringsAsFactors = FALSE),
    data.frame(seed = base_df$seed, qc = base_df$qc, component = "In vitro objective", value = base_df$invitro_objective, stringsAsFactors = FALSE)
  )
  comp_long <- comp_long[is.finite(comp_long$qc) & is.finite(comp_long$value), , drop = FALSE]
  if (!nrow(comp_long)) {
    stop("No finite qc/objective pairs were available for QC-objective diagnostics.")
  }

  component_levels <- c("Total objective", "In vivo objective", "In vitro objective")
  rho_labels <- vapply(component_levels, function(component) {
    part <- comp_long[comp_long$component == component, , drop = FALSE]
    paste(component, format_spearman_label(part$qc, part$value), sep = "\n")
  }, character(1), USE.NAMES = TRUE)
  comp_long$component_label <- factor(rho_labels[comp_long$component], levels = unname(rho_labels))

  component_colors <- c(
    "Total objective" = "#4D4D4D",
    "In vivo objective" = "#1f77b4",
    "In vitro objective" = "#d95f02"
  )
  p <- ggplot(comp_long, aes(x = qc, y = value, color = component)) +
    geom_point(size = 2.2, alpha = 0.82) +
    facet_wrap(~ component_label, scales = "free_y", nrow = 1) +
    scale_color_manual(values = component_colors) +
    labs(
      title = paste0("QC vs Joint Objective Components: ", run_label),
      subtitle = "Each point is one seed; trend lines are linear fits. Lower objective is better.",
      x = "qc",
      y = "Objective",
      color = "Component"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 9)
    )
  p <- add_linear_trend(p, comp_long, x_col = "qc", group_col = "component")

  ggplot2::ggsave(out_path, p, width = 10.5, height = 4.2)
  invisible(out_path)
}

plot_joint_qc_parameter_diagnostics <- function(summary_df, out_path, run_label) {
  if (!("value__qc" %in% names(summary_df))) {
    stop("seed_summary is missing required column for QC-parameter diagnostics: value__qc")
  }
  parameter_map <- data.frame(
    column = c("value__p_misseg", "value__mu_hp", "value__gamma_growth", "value__buffer_n_exp"),
    parameter = c("p_misseg", "mu_hp", "gamma_growth", "buffer_n_exp"),
    stringsAsFactors = FALSE
  )
  missing <- setdiff(parameter_map$column, names(summary_df))
  if (length(missing)) {
    stop("seed_summary is missing required columns for QC-parameter diagnostics: ", paste(missing, collapse = ", "))
  }

  qc <- suppressWarnings(as.numeric(summary_df$value__qc))
  param_long <- do.call(rbind, lapply(seq_len(nrow(parameter_map)), function(i) {
    col <- parameter_map$column[[i]]
    data.frame(
      seed = as.character(summary_df$seed),
      qc = qc,
      parameter = parameter_map$parameter[[i]],
      value = suppressWarnings(as.numeric(summary_df[[col]])),
      stringsAsFactors = FALSE
    )
  }))
  param_long <- param_long[is.finite(param_long$qc) & is.finite(param_long$value), , drop = FALSE]
  if (!nrow(param_long)) {
    stop("No finite qc/parameter pairs were available for QC-parameter diagnostics.")
  }

  parameter_levels <- parameter_map$parameter
  rho_labels <- vapply(parameter_levels, function(parameter) {
    part <- param_long[param_long$parameter == parameter, , drop = FALSE]
    paste(parameter, format_spearman_label(part$qc, part$value), sep = "\n")
  }, character(1), USE.NAMES = TRUE)
  param_long$parameter_label <- factor(rho_labels[param_long$parameter], levels = unname(rho_labels))

  parameter_colors <- c(
    "p_misseg" = "#1f77b4",
    "mu_hp" = "#d95f02",
    "gamma_growth" = "#2ca02c",
    "buffer_n_exp" = "#756bb1"
  )
  p <- ggplot(param_long, aes(x = qc, y = value, color = parameter)) +
    geom_point(size = 2.2, alpha = 0.82) +
    facet_wrap(~ parameter_label, scales = "free_y", nrow = 2) +
    scale_color_manual(values = parameter_colors) +
    labs(
      title = paste0("QC vs Key Joint Parameters: ", run_label),
      subtitle = "Each point is one seed; trend lines are linear fits on natural parameter values.",
      x = "qc",
      y = "Parameter value",
      color = "Parameter"
    ) +
    theme_bw(base_size = 11) +
    theme(
      aspect.ratio = 1,
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.text = element_text(size = 9)
    )
  p <- add_linear_trend(p, param_long, x_col = "qc", group_col = "parameter")

  ggplot2::ggsave(out_path, p, width = 8.6, height = 8.8)
  invisible(out_path)
}

infer_run_label_from_paths <- function(seed_summary_path = NULL, extra_results_dir = NULL) {
  if (!is.null(extra_results_dir) && nzchar(extra_results_dir)) {
    base <- basename(normalizePath(extra_results_dir, mustWork = FALSE))
    if (identical(base, "extra_results")) return(basename(dirname(normalizePath(extra_results_dir, mustWork = FALSE))))
    return(base)
  }
  if (!is.null(seed_summary_path) && nzchar(seed_summary_path)) {
    return(basename(dirname(dirname(normalizePath(seed_summary_path, mustWork = FALSE)))))
  }
  "joint_fit"
}

resolve_seed_summary_path <- function(seed_summary_path = NULL, extra_results_dir = NULL) {
  if (!is.null(seed_summary_path) && nzchar(trimws(seed_summary_path))) {
    return(normalizePath(seed_summary_path, mustWork = TRUE))
  }
  if (is.null(extra_results_dir) || !nzchar(trimws(extra_results_dir))) {
    stop("Provide either --seed_summary=/path/to/seed_summary.tsv or --extra_results_dir=/path/to/extra_results.")
  }
  normalizePath(file.path(extra_results_dir, "seed_summary.tsv"), mustWork = TRUE)
}

plot_joint_qc_diagnostics <- function(seed_summary_path = NULL,
                                      extra_results_dir = NULL,
                                      out_dir = NULL,
                                      run_label = NULL) {
  seed_summary_path <- resolve_seed_summary_path(
    seed_summary_path = seed_summary_path,
    extra_results_dir = extra_results_dir
  )
  if (is.null(out_dir) || !nzchar(trimws(out_dir))) {
    out_dir <- dirname(seed_summary_path)
  }
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(run_label) || !nzchar(trimws(run_label))) {
    run_label <- infer_run_label_from_paths(
      seed_summary_path = seed_summary_path,
      extra_results_dir = extra_results_dir
    )
  }

  seed_summary <- utils::read.delim(seed_summary_path, check.names = FALSE, stringsAsFactors = FALSE)
  objective_path <- file.path(out_dir, "joint_qc_objective_diagnostics.pdf")
  parameter_path <- file.path(out_dir, "joint_qc_parameter_diagnostics.pdf")
  plot_joint_qc_objective_diagnostics(seed_summary, objective_path, run_label)
  plot_joint_qc_parameter_diagnostics(seed_summary, parameter_path, run_label)
  c(objective = normalizePath(objective_path, mustWork = TRUE), parameter = normalizePath(parameter_path, mustWork = TRUE))
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  paths <- plot_joint_qc_diagnostics(
    seed_summary_path = argv$seed_summary,
    extra_results_dir = argv$extra_results_dir,
    out_dir = argv$out_dir,
    run_label = argv$run_label
  )
  message("Wrote QC-objective diagnostics plot: ", paths[["objective"]])
  message("Wrote QC-parameter diagnostics plot: ", paths[["parameter"]])
}

if (sys.nframe() == 0) {
  main()
}
