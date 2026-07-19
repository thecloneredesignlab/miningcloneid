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
safe_read <- function(path) o2fr_read_tsv(path, optional = TRUE)
save_plot <- function(plot, path, width = 10, height = 7) {
  ggplot2::ggsave(path, plot, width = width, height = height, bg = "white")
  if (file.exists(path)) path else character()
}

plot_objective_distributions <- function(long, path, title = "Objective Distributions Across Seeds") {
  if (!nrow(long)) return(character())
  long$metric_label <- factor(long$metric_label, levels = unique(as.character(long$metric_label)))
  p <- ggplot(long, aes(metric_label, value, fill = metric_label)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, linewidth = 0.35) +
    labs(title = title, x = NULL, y = "Objective value") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "none")
  save_plot(p, path, 9, 7)
}

plot_objective_risk <- function(summary, path) {
  keep <- is.finite(summary$objective) & is.finite(summary$min_rel_dist_active) & summary$min_rel_dist_active > 0
  d <- summary[keep, , drop = FALSE]
  if (!nrow(d)) return(character())
  p <- ggplot(d, aes(objective, min_rel_dist_active, label = seed)) +
    geom_point(size = 2.4, color = "#2b6cb0") + geom_text(size = 2.5, nudge_y = 0.01, check_overlap = TRUE) +
    labs(title = "Objective vs Boundary Risk", x = "Objective", y = "Minimum relative distance to fitted bound") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  save_plot(p, path, 10, 7)
}

plot_boundary_forest <- function(parameters, summary, path, gated = FALSE, log_x = FALSE) {
  d <- parameters[parameters$active_in_fit %in% TRUE & is.finite(parameters$rel_pos_plot), , drop = FALSE]
  if (!nrow(d)) return(character())
  top <- if (gated) {
    s <- summary[summary$pred1000_both_gt44 %in% TRUE, , drop = FALSE]
    utils::head(as.character(s$seed[order(s$objective, s$seed)]), 3L)
  } else utils::head(as.character(summary$seed[order(summary$objective, summary$seed)]), 3L)
  d$highlight <- as.character(d$seed) %in% top
  d$parameter <- factor(d$param_prototype, levels = rev(unique(as.character(d$param_prototype))))
  if (log_x) {
    d$x <- suppressWarnings(as.numeric(d$prototype_value))
    d <- d[is.finite(d$x) & d$x > 0, , drop = FALSE]
    if (!nrow(d)) return(character())
    p <- ggplot(d, aes(x, parameter, color = highlight)) + geom_point(alpha = 0.72, size = 1.8) + scale_x_log10()
    x_label <- "Natural-scale fitted value (log10 axis)"
  } else {
    d$x <- d$rel_pos_plot
    p <- ggplot(d, aes(x, parameter, color = highlight)) + geom_vline(xintercept = c(0, 1), color = "grey75") + geom_point(alpha = 0.72, size = 1.8) + coord_cartesian(xlim = c(0, 1))
    x_label <- "Relative position within transformed bounds"
  }
  p <- p + scale_color_manual(values = c(`FALSE` = "grey58", `TRUE` = "#d7301f"), guide = "none") +
    labs(title = if (gated) "Parameter Boundary Forest: top prediction-gated seeds" else "Parameter Boundary Forest", x = x_label, y = NULL) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  save_plot(p, path, 12, 8)
}

plot_ploidy_summary <- function(summary, path) {
  if (!nrow(summary)) return(character())
  summary$cohort <- factor(summary$cohort, levels = c("2N", "4N"))
  p <- ggplot(summary, aes(day, mean_value, color = cohort, fill = cohort)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.18, color = NA) +
    geom_line(aes(y = min_value), linetype = "dashed", linewidth = 0.45) +
    geom_line(aes(y = max_value), linetype = "dashed", linewidth = 0.45) + geom_line(linewidth = 0.9) +
    scale_color_manual(values = c(`2N` = "#1f77b4", `4N` = "#d62728")) + scale_fill_manual(values = c(`2N` = "#1f77b4", `4N` = "#d62728")) +
    scale_y_continuous("Weighted mean chromosome number", sec.axis = sec_axis(~ . / 22, name = "Ploidy")) +
    labs(title = "1000-day Ploidy Prediction Mean with 95% CI", x = "Day", color = "Cohort", fill = "Cohort") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  save_plot(p, path, 12, 6)
}

plot_seed_curves <- function(curves, summary, cohort, value_col, y_label, path, log10_y = FALSE) {
  d <- curves[curves$cohort == cohort & is.finite(curves[[value_col]]), , drop = FALSE]
  s <- summary[summary$cohort == cohort & is.finite(summary$mean_value), , drop = FALSE]
  if (!nrow(d) || !nrow(s)) return(character())
  if (log10_y) {
    d <- d[d[[value_col]] > 0, , drop = FALSE]; s <- s[s$mean_value > 0, , drop = FALSE]
    d$plot_value <- log10(d[[value_col]]); s$plot_value <- log10(s$mean_value)
  } else { d$plot_value <- d[[value_col]]; s$plot_value <- s$mean_value }
  p <- ggplot() + geom_line(data = d, aes(day, plot_value, group = seed), color = "grey72", linewidth = 0.18, alpha = 0.55) +
    geom_line(data = s, aes(day, plot_value), color = if (cohort == "2N") "#1f77b4" else "#d62728", linewidth = 1) +
    labs(title = paste0("1000-day prediction seed trajectories: ", cohort), x = "Day", y = y_label) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  save_plot(p, path, 12, 6)
}

plot_joint_outputs <- function(summary, soft, out_dir) {
  outputs <- character()
  if (all(c("objective_invivo", "objective_invitro") %in% names(summary)) && any(is.finite(summary$objective_invivo) & is.finite(summary$objective_invitro))) {
    d <- summary[is.finite(summary$objective_invivo) & is.finite(summary$objective_invitro), , drop = FALSE]
    p <- ggplot(d, aes(objective_invivo, objective_invitro, color = objective, label = seed)) + geom_point(size = 2.5) + geom_text(size = 2.4, nudge_y = 0.01, check_overlap = TRUE) + theme_bw() + labs(title = "Joint Objective Tradeoff")
    outputs <- c(outputs, save_plot(p, file.path(out_dir, "joint_objective_tradeoff.pdf"), 9, 7))
  }
  if (nrow(soft)) {
    if ("delta_transformed" %in% names(soft)) {
      soft$abs_delta <- abs(suppressWarnings(as.numeric(soft$delta_transformed)))
      p <- ggplot(soft, aes(parameter, abs_delta)) + geom_boxplot() + geom_jitter(width = 0.12, alpha = 0.5) + theme_bw() + labs(title = "Joint Soft-Coupling Delta Magnitudes", x = NULL, y = "Absolute transformed delta")
      outputs <- c(outputs, save_plot(p, file.path(out_dir, "joint_soft_coupling_delta_magnitude.pdf"), 9.5, 6.2))
    }
    if ("penalty_paid" %in% names(soft)) {
      soft$penalty_paid <- suppressWarnings(as.numeric(soft$penalty_paid))
      p <- ggplot(soft, aes(parameter, penalty_paid)) + geom_boxplot() + geom_jitter(width = 0.12, alpha = 0.5) + theme_bw() + labs(title = "Joint Soft-Coupling Penalty Ranking", x = NULL, y = "Penalty paid")
      outputs <- c(outputs, save_plot(p, file.path(out_dir, "joint_soft_coupling_penalty_ranking.pdf"), 9.5, 6.2))
    }
  }
  outputs
}

main <- function(argv = parse_args()) {
  simulation_dir <- normalizePath(argv$simulation_dir, mustWork = TRUE)
  analysis_dir <- normalizePath(argv$analysis_dir, mustWork = TRUE)
  if (!file.exists(file.path(simulation_dir, "simulation_manifest.tsv")) || !file.exists(file.path(analysis_dir, "analysis_manifest.tsv"))) stop("Extra-results visualization requires completed simulation and analysis manifests.", call. = FALSE)
  out_dir <- normalizePath(argv$viz_dir %||% argv$out_dir %||% analysis_dir, mustWork = FALSE)
  summary <- safe_read(file.path(analysis_dir, "seed_summary.tsv"))
  parameters <- safe_read(file.path(analysis_dir, "parameter_boundary_long.tsv"))
  objectives <- safe_read(file.path(analysis_dir, "objective_components_long.tsv"))
  ploidy <- safe_read(file.path(simulation_dir, "predict1000_ploidy_seed_day_mean.tsv"))
  ploidy_ci <- safe_read(file.path(simulation_dir, "predict1000_ploidy_mean_ci.tsv"))
  burden <- safe_read(file.path(simulation_dir, "predict1000_burden_total_seed_day_mean.tsv"))
  burden_ci <- safe_read(file.path(simulation_dir, "predict1000_burden_total_mean_ci.tsv"))
  soft <- safe_read(file.path(analysis_dir, "joint_soft_coupling_all.tsv"))
  if (!nrow(summary) || !nrow(parameters) || !nrow(objectives)) stop("Missing core extra-results analysis tables.", call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- c(
    plot_objective_distributions(objectives, file.path(out_dir, "objective_components_violin.pdf")),
    plot_objective_risk(summary, file.path(out_dir, "objective_vs_boundary_risk.pdf")),
    plot_boundary_forest(parameters, summary, file.path(out_dir, "parameter_boundary_forest.pdf")),
    plot_boundary_forest(parameters, summary, file.path(out_dir, "parameter_boundary_forest_log_x.pdf"), log_x = TRUE),
    plot_boundary_forest(parameters, summary, file.path(out_dir, "parameter_boundary_forest_pred1000_gt44_top3.pdf"), gated = TRUE),
    plot_boundary_forest(parameters, summary, file.path(out_dir, "parameter_boundary_forest_pred1000_gt44_top3_log_x.pdf"), gated = TRUE, log_x = TRUE),
    plot_ploidy_summary(ploidy_ci, file.path(out_dir, "predict1000_ploidy_mean_ci_2N_4N.pdf")),
    plot_seed_curves(ploidy, ploidy_ci, "2N", "ploidy_value", "Weighted mean chromosome number", file.path(out_dir, "predict1000_ploidy_seed_curves_2N.pdf")),
    plot_seed_curves(ploidy, ploidy_ci, "4N", "ploidy_value", "Weighted mean chromosome number", file.path(out_dir, "predict1000_ploidy_seed_curves_4N.pdf")),
    plot_seed_curves(burden, burden_ci, "2N", "burden_value", "log10 total tumor burden (mm^3)", file.path(out_dir, "predict1000_burden_total_log_seed_mean_2N.pdf"), TRUE),
    plot_seed_curves(burden, burden_ci, "4N", "burden_value", "log10 total tumor burden (mm^3)", file.path(out_dir, "predict1000_burden_total_log_seed_mean_4N.pdf"), TRUE),
    plot_joint_outputs(summary, soft, out_dir)
  )
  outputs <- unique(outputs[file.exists(outputs)])
  utils::write.table(data.frame(stage = "visualization", file = basename(outputs), stringsAsFactors = FALSE), file.path(out_dir, "visualization_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote extra-results visualizations: ", out_dir)
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
