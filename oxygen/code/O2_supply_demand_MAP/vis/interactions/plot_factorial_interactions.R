#!/usr/bin/env Rscript

# Pure visualization consumer for factorial in vivo interaction outputs.
# It never reads fitted parameters or runs the simulation producer.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

.interaction_script_dir <- local({
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
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
})

WORKFLOW_ROOT <- normalizePath(file.path(.interaction_script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_perturbation_utils.R"), local = environment())
rm(.interaction_script_dir)

format_compact_value_label <- function(x, digits = 3) {
  vapply(as.numeric(x), function(z) {
    if (!is.finite(z)) return(as.character(z))
    use_scientific <- z != 0 && abs(z) < 1e-3
    format(signif(z, digits), trim = TRUE, scientific = use_scientific)
  }, character(1))
}

factor_labels <- function(x, digits = 4) {
  factor(format_value_label(x, digits = digits), levels = format_value_label(sort(unique(x)), digits = digits))
}

add_plot_labels <- function(df, run_params, condition_pmiss_prefix = "post p_miss=") {
  if ("p_wgd_post" %in% names(df)) {
    df$.plot_p_wgd <- as.numeric(df$p_wgd_post)
  } else {
    df$.plot_p_wgd <- as.numeric(df$p_wgd)
  }
  if ("p_mis_base_post" %in% names(df)) {
    df$.plot_p_mis_base <- as.numeric(df$p_mis_base_post)
  } else {
    df$.plot_p_mis_base <- as.numeric(df$p_mis_base)
  }

  baseline_pwgd <- as.numeric(run_params$p_wgd)
  pwgd_values <- sort(unique(df$.plot_p_wgd))
  pwgd_ratio <- pwgd_values / baseline_pwgd
  pwgd_map <- setNames(
    paste0(
      format_compact_value_label(pwgd_ratio),
      "x fitted p_wgd (",
      format_compact_value_label(pwgd_values),
      ")"
    ),
    format_value_label(pwgd_values)
  )
  out <- df %>%
    mutate(
      o2_label = factor(
        paste0("O2=", format_compact_value_label(o2_S0)),
        levels = paste0("O2=", format_compact_value_label(sort(unique(o2_S0))))
      ),
      trigger_label = factor(
        paste0("trigger=", format(signif(trigger_burden_cells, 3), scientific = TRUE)),
        levels = paste0("trigger=", format(signif(sort(unique(trigger_burden_cells)), 3), scientific = TRUE))
      ),
      pmiss_label = factor_labels(.plot_p_mis_base),
      pwgd_label_raw = format_value_label(.plot_p_wgd),
      pwgd_label = factor(unname(pwgd_map[pwgd_label_raw]), levels = unname(pwgd_map[format_value_label(pwgd_values)]))
    )
  condition_grid <- expand.grid(
    pmiss_label = levels(out$pmiss_label),
    pwgd_label = levels(out$pwgd_label),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  condition_levels <- paste(condition_grid$pwgd_label, paste0(condition_pmiss_prefix, condition_grid$pmiss_label), sep = " / ")
  out %>%
    mutate(
      condition_label = factor(paste(pwgd_label, paste0(condition_pmiss_prefix, pmiss_label), sep = " / "), levels = condition_levels),
      condition_index = as.integer(condition_label)
    )
}

scale_fill_log10_burden <- function(name, limits = NULL) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#FFF7BC", "#FEC44F", "#FB6A4A", "#CC2C7A", "#4A1486"),
    limits = limits,
    oob = scales::squish,
    na.value = "grey85",
    name = name
  )
}

scale_fill_mean_ploidy <- function(name) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#F3E8FF", "#C084FC", "#7C3AED", "#1D4ED8"),
    values = scales::rescale(c(0, 1.5, 3, 4.5, 6)),
    limits = c(0, 6),
    oob = scales::squish,
    na.value = "grey85",
    name = name
  )
}

reset_generated_dir <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

complete_timecourse_for_plot <- function(plot_df, value_col, fill_value = 0, max_day = NULL) {
  if (is.null(max_day)) {
    max_day <- ceiling(max(plot_df$day, na.rm = TRUE))
  } else {
    max_day <- ceiling(as.numeric(max_day))
  }
  if (!is.finite(max_day) || max_day < 0) max_day <- 0
  scenario_ids <- unique(as.character(plot_df$scenario_id))
  meta_cols <- c(
    "scenario_id", "experiment", "trigger_burden_cells", "trigger_day",
    "actual_trigger_burden_cells", "o2_S0", "p_mis_base_pre",
    "p_mis_base_post", "p_wgd_pre", "p_wgd_post", "p_wgd", "status", "o2_label", "trigger_label",
    "pmiss_label", "pwgd_label", "condition_label", "condition_index", "vary_label", "vary_index",
    "unit_o2_label", "unit_o2_index", "unit_vary_label", "unit_vary_index"
  )
  meta_cols <- intersect(meta_cols, names(plot_df))
  meta <- plot_df[!duplicated(as.character(plot_df$scenario_id)), meta_cols, drop = FALSE]
  meta$scenario_id <- as.character(meta$scenario_id)
  values <- plot_df[, c("scenario_id", "day", value_col), drop = FALSE]
  values$scenario_id <- as.character(values$scenario_id)
  values$.observed_value <- TRUE
  grid <- expand.grid(
    scenario_id = scenario_ids,
    day = seq(0, max_day, by = 1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  out <- merge(grid, meta, by = "scenario_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, values, by = c("scenario_id", "day"), all.x = TRUE, sort = FALSE)
  missing_from_grid <- is.na(out$.observed_value)
  missing_value <- is.na(out[[value_col]])
  if (is.numeric(out[[value_col]])) {
    missing_value <- missing_value | !is.finite(out[[value_col]])
  }
  out[[value_col]][missing_from_grid | missing_value] <- fill_value
  out$.observed_value <- NULL
  out
}

make_treatment_line_segments <- function(endpoint_plot,
                                         row_index_col = "condition_index",
                                         show_treatment_lines = TRUE) {
  if (!isTRUE(show_treatment_lines) || !nrow(endpoint_plot)) return(endpoint_plot[FALSE, , drop = FALSE])
  if (!row_index_col %in% names(endpoint_plot)) return(endpoint_plot[FALSE, , drop = FALSE])
  endpoint_plot <- endpoint_plot[is.finite(endpoint_plot$trigger_day), , drop = FALSE]
  if (!nrow(endpoint_plot)) return(endpoint_plot)

  rows <- vector("list", nrow(endpoint_plot))
  for (i in seq_len(nrow(endpoint_plot))) {
    dr <- endpoint_plot[i, , drop = FALSE]
    treatment_days <- as.numeric(dr$trigger_day)
    is_intermittent <- "intermittent_on_day" %in% names(dr) &&
      "intermittent_off_day" %in% names(dr) &&
      "post_treatment_duration_day" %in% names(dr) &&
      is.finite(as.numeric(dr$intermittent_on_day)) &&
      is.finite(as.numeric(dr$intermittent_off_day)) &&
      is.finite(as.numeric(dr$post_treatment_duration_day)) &&
      as.numeric(dr$intermittent_on_day) > 0 &&
      as.numeric(dr$intermittent_off_day) >= 0

    if (isTRUE(is_intermittent)) {
      period <- as.numeric(dr$intermittent_on_day) + as.numeric(dr$intermittent_off_day)
      if (period > 0) {
        n_starts <- max(1L, as.integer(ceiling((as.numeric(dr$post_treatment_duration_day) - 1e-9) / period)))
        treatment_days <- as.numeric(dr$trigger_day) + seq(0, by = period, length.out = n_starts)
      }
    }

    line_rows <- dr[rep(1, length(treatment_days)), , drop = FALSE]
    line_rows$treatment_day <- treatment_days
    line_rows$y0 <- as.numeric(line_rows[[row_index_col]]) - 0.48
    line_rows$y1 <- as.numeric(line_rows[[row_index_col]]) + 0.48
    rows[[i]] <- line_rows
  }

  bind_rows(rows)
}

plot_one_split_timecourse <- function(plot_df,
                                      endpoint_plot,
                                      fixed_col,
                                      fixed_value,
                                      vary_col,
                                      value_col,
                                      out_pdf,
                                      include_trigger = TRUE,
                                      show_treatment_lines = TRUE,
                                      fill_kind = c("burden", "ploidy"),
                                      y_axis_label = NULL) {
  fill_kind <- match.arg(fill_kind)
  fixed_chr <- as.character(fixed_value)
  sub <- plot_df[as.character(plot_df[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(invisible(FALSE))

  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  sub$vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)
  sub$vary_index <- as.integer(sub$vary_label)
  sub_full <- complete_timecourse_for_plot(sub, value_col, fill_value = 0)

  endpoint_sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  endpoint_sub$vary_label <- factor(as.character(endpoint_sub[[vary_col]]), levels = vary_levels)
  endpoint_sub$vary_index <- as.integer(endpoint_sub$vary_label)
  treatment_lines <- make_treatment_line_segments(
    endpoint_sub,
    row_index_col = "vary_index",
    show_treatment_lines = show_treatment_lines
  )

  facet_layer <- if (include_trigger) facet_grid(trigger_label ~ o2_label) else facet_grid(. ~ o2_label)
  y_axis_label <- y_axis_label %||% vary_col
  p <- ggplot(sub_full, aes(day, vary_index, fill = .data[[value_col]])) +
    geom_raster() +
    facet_layer +
    scale_y_continuous(breaks = seq_along(vary_levels), labels = vary_levels, expand = c(0, 0)) +
    labs(x = "Day", y = y_axis_label) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.ticks.y = element_line(linewidth = 0.2),
      strip.text = element_text(size = 7)
    )
  if (identical(fill_kind, "burden")) {
    p <- p + scale_fill_log10_burden("log10 burden")
  } else {
    p <- p + scale_fill_mean_ploidy("mean ploidy")
  }
  if (nrow(treatment_lines) > 0L) {
    p <- p + geom_segment(
      data = treatment_lines,
      aes(x = treatment_day, xend = treatment_day, y = y0, yend = y1),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.25
    )
  }

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  plot_height <- if (include_trigger) 10 else 4.5
  ggsave(out_pdf, p, width = 20, height = plot_height)
  invisible(TRUE)
}

plot_split_timecourse_heatmaps <- function(endpoint_summary,
                                           burden_all,
                                           ploidy_summary,
                                           run_params,
                                           out_dir,
                                           include_trigger = TRUE,
                                           show_treatment_lines = TRUE,
                                           condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  burden_plot <- add_plot_labels(burden_all, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  ploidy_plot <- add_plot_labels(ploidy_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)

  split_root <- file.path(out_dir, "timecourse_split")
  reset_generated_dir(split_root)

  for (pmiss in levels(burden_plot$pmiss_label)) {
    sub_dir <- file.path(split_root, "fixed_p_mis_base", safe_id("p_miss", pmiss))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(data.frame(fixed_parameter = "p_mis_base", fixed_value = pmiss), file.path(sub_dir, "fixed_value.tsv"))
    plot_one_split_timecourse(
      burden_plot, endpoint_plot,
      fixed_col = "pmiss_label", fixed_value = pmiss,
      vary_col = "pwgd_label", value_col = "pred_log10_burden_cells",
      out_pdf = file.path(sub_dir, "timecourse_log10_burden_by_p_wgd.pdf"),
      include_trigger = include_trigger,
      show_treatment_lines = show_treatment_lines,
      fill_kind = "burden",
      y_axis_label = "p_wgd relative to fitted (actual p_wgd)"
    )
    plot_one_split_timecourse(
      ploidy_plot, endpoint_plot,
      fixed_col = "pmiss_label", fixed_value = pmiss,
      vary_col = "pwgd_label", value_col = "mean_ploidy",
      out_pdf = file.path(sub_dir, "timecourse_mean_ploidy_by_p_wgd.pdf"),
      include_trigger = include_trigger,
      show_treatment_lines = show_treatment_lines,
      fill_kind = "ploidy",
      y_axis_label = "p_wgd relative to fitted (actual p_wgd)"
    )
  }

  for (pwgd in levels(burden_plot$pwgd_label)) {
    sub_dir <- file.path(split_root, "fixed_p_wgd", safe_id("p_wgd", pwgd))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(data.frame(fixed_parameter = "p_wgd", fixed_value = pwgd), file.path(sub_dir, "fixed_value.tsv"))
    plot_one_split_timecourse(
      burden_plot, endpoint_plot,
      fixed_col = "pwgd_label", fixed_value = pwgd,
      vary_col = "pmiss_label", value_col = "pred_log10_burden_cells",
      out_pdf = file.path(sub_dir, "timecourse_log10_burden_by_p_mis_base.pdf"),
      include_trigger = include_trigger,
      show_treatment_lines = show_treatment_lines,
      fill_kind = "burden",
      y_axis_label = "p_mis_base"
    )
    plot_one_split_timecourse(
      ploidy_plot, endpoint_plot,
      fixed_col = "pwgd_label", fixed_value = pwgd,
      vary_col = "pmiss_label", value_col = "mean_ploidy",
      out_pdf = file.path(sub_dir, "timecourse_mean_ploidy_by_p_mis_base.pdf"),
      include_trigger = include_trigger,
      show_treatment_lines = show_treatment_lines,
      fill_kind = "ploidy",
      y_axis_label = "p_mis_base"
    )
  }
  invisible(split_root)
}

plot_one_split_endpoint <- function(endpoint_plot,
                                    fixed_col,
                                    fixed_value,
                                    vary_col,
                                    value_col,
                                    out_pdf,
                                    include_trigger = TRUE,
                                    fill_kind = c("burden", "ploidy"),
                                    y_axis_label = NULL) {
  fill_kind <- match.arg(fill_kind)
  fixed_chr <- as.character(fixed_value)
  sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(invisible(FALSE))

  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  sub$vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)

  y_axis_label <- y_axis_label %||% vary_col
  p <- ggplot(sub, aes(o2_label, vary_label, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.15) +
    labs(x = "O2_S0", y = y_axis_label) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 7)
    )
  if (include_trigger) {
    p <- p + facet_grid(trigger_label ~ .)
  }
  if (identical(fill_kind, "burden")) {
    p <- p + scale_fill_log10_burden("Endpoint\nlog10 burden")
  } else {
    p <- p + scale_fill_mean_ploidy("Endpoint\nmean ploidy")
  }

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  plot_height <- if (include_trigger) 8 else 4.2
  ggsave(out_pdf, p, width = 8, height = plot_height)
  invisible(TRUE)
}

plot_split_endpoint_heatmaps <- function(endpoint_summary,
                                         run_params,
                                         out_dir,
                                         include_trigger = TRUE,
                                         condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  split_root <- file.path(out_dir, "endpoint_split")
  reset_generated_dir(split_root)

  for (pmiss in levels(endpoint_plot$pmiss_label)) {
    sub_dir <- file.path(split_root, "fixed_p_mis_base", safe_id("p_miss", pmiss))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(data.frame(fixed_parameter = "p_mis_base", fixed_value = pmiss), file.path(sub_dir, "fixed_value.tsv"))
    plot_one_split_endpoint(
      endpoint_plot,
      fixed_col = "pmiss_label", fixed_value = pmiss,
      vary_col = "pwgd_label", value_col = "pred_log10_burden_cells",
      out_pdf = file.path(sub_dir, "endpoint_log10_burden_by_o2.pdf"),
      include_trigger = include_trigger,
      fill_kind = "burden",
      y_axis_label = "p_wgd relative to fitted (actual p_wgd)"
    )
    plot_one_split_endpoint(
      endpoint_plot,
      fixed_col = "pmiss_label", fixed_value = pmiss,
      vary_col = "pwgd_label", value_col = "mean_ploidy",
      out_pdf = file.path(sub_dir, "endpoint_mean_ploidy_by_o2.pdf"),
      include_trigger = include_trigger,
      fill_kind = "ploidy",
      y_axis_label = "p_wgd relative to fitted (actual p_wgd)"
    )
  }

  for (pwgd in levels(endpoint_plot$pwgd_label)) {
    sub_dir <- file.path(split_root, "fixed_p_wgd", safe_id("p_wgd", pwgd))
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(data.frame(fixed_parameter = "p_wgd", fixed_value = pwgd), file.path(sub_dir, "fixed_value.tsv"))
    plot_one_split_endpoint(
      endpoint_plot,
      fixed_col = "pwgd_label", fixed_value = pwgd,
      vary_col = "pmiss_label", value_col = "pred_log10_burden_cells",
      out_pdf = file.path(sub_dir, "endpoint_log10_burden_by_o2.pdf"),
      include_trigger = include_trigger,
      fill_kind = "burden",
      y_axis_label = "p_mis_base"
    )
    plot_one_split_endpoint(
      endpoint_plot,
      fixed_col = "pwgd_label", fixed_value = pwgd,
      vary_col = "pmiss_label", value_col = "mean_ploidy",
      out_pdf = file.path(sub_dir, "endpoint_mean_ploidy_by_o2.pdf"),
      include_trigger = include_trigger,
      fill_kind = "ploidy",
      y_axis_label = "p_mis_base"
    )
  }
  invisible(split_root)
}

plot_endpoint_unit_metric <- function(sub,
                                      vary_col,
                                      value_col,
                                      fill_kind = c("burden", "ploidy"),
                                      y_axis_label = NULL,
                                      show_y_axis = TRUE,
                                      burden_limits = NULL) {
  fill_kind <- match.arg(fill_kind)
  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  sub$vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)
  y_axis_label <- y_axis_label %||% vary_col

  p <- ggplot(sub, aes(o2_label, vary_label, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.12) +
    coord_fixed(ratio = 1) +
    labs(x = "O2_S0", y = if (show_y_axis) y_axis_label else NULL) +
    theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 7, hjust = 0.5)
    )
  if (!show_y_axis) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  if (identical(fill_kind, "burden")) {
    p <- p + scale_fill_log10_burden("Endpoint\nlog10 burden", limits = burden_limits)
  } else {
    p <- p + scale_fill_mean_ploidy("Endpoint\nmean ploidy")
  }
  p
}

make_endpoint_unit_plot <- function(endpoint_plot,
                                    fixed_col,
                                    fixed_value,
                                    vary_col,
                                    fixed_axis_label,
                                    vary_axis_label,
                                    burden_limits = NULL) {
  fixed_chr <- as.character(fixed_value)
  sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  unit_title <- paste(strwrap(paste0(fixed_axis_label, " = ", fixed_chr), width = 34), collapse = "\n")

  p_burden <- plot_endpoint_unit_metric(
    sub,
    vary_col = vary_col,
    value_col = "pred_log10_burden_cells",
    fill_kind = "burden",
    y_axis_label = vary_axis_label,
    show_y_axis = TRUE,
    burden_limits = burden_limits
  ) +
    ggtitle(paste("burden", unit_title, sep = "\n"))
  p_ploidy <- plot_endpoint_unit_metric(
    sub,
    vary_col = vary_col,
    value_col = "mean_ploidy",
    fill_kind = "ploidy",
    y_axis_label = vary_axis_label,
    show_y_axis = FALSE,
    burden_limits = burden_limits
  ) +
    ggtitle("ploidy")

  wrap_plots(p_burden, p_ploidy, ncol = 2)
}

make_endpoint_unit_grid_plot <- function(endpoint_plot,
                                         fixed_col,
                                         vary_col,
                                         fixed_axis_label,
                                         vary_axis_label,
                                         fixed_levels,
                                         burden_limits = NULL,
                                         page_title = NULL) {
  unit_plots <- lapply(fixed_levels, function(fixed_value) {
    make_endpoint_unit_plot(
      endpoint_plot = endpoint_plot,
      fixed_col = fixed_col,
      fixed_value = fixed_value,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      burden_limits = burden_limits
    )
  })
  unit_plots <- Filter(Negate(is.null), unit_plots)
  if (!length(unit_plots)) return(NULL)

  wrap_plots(unit_plots, ncol = 3, guides = "collect") +
    plot_layout(guides = "collect") +
    plot_annotation(title = page_title)
}

save_endpoint_unit_grid_pages <- function(endpoint_plot,
                                          fixed_col,
                                          vary_col,
                                          fixed_axis_label,
                                          vary_axis_label,
                                          fixed_levels,
                                          out_pdf,
                                          include_trigger = TRUE) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  burden_limits <- range(endpoint_plot$pred_log10_burden_cells, na.rm = TRUE)
  if (any(!is.finite(burden_limits))) burden_limits <- NULL

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  grDevices::pdf(out_pdf, width = 22, height = 16, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_group in page_groups) {
    page_df <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_title <- if (identical(page_group, "all")) {
      paste0(fixed_axis_label, " endpoint unit grid")
    } else {
      paste0(fixed_axis_label, " endpoint unit grid / ", page_group)
    }
    page_plot <- make_endpoint_unit_grid_plot(
      endpoint_plot = page_df,
      fixed_col = fixed_col,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      fixed_levels = fixed_levels,
      burden_limits = burden_limits,
      page_title = page_title
    )
    if (!is.null(page_plot)) print(page_plot)
  }
  invisible(out_pdf)
}

save_endpoint_unit_individual_plots <- function(endpoint_plot,
                                                fixed_col,
                                                vary_col,
                                                fixed_axis_label,
                                                vary_axis_label,
                                                fixed_levels,
                                                out_dir,
                                                include_trigger = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  burden_limits <- range(endpoint_plot$pred_log10_burden_cells, na.rm = TRUE)
  if (any(!is.finite(burden_limits))) burden_limits <- NULL

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  manifest_rows <- list()
  row_i <- 1L
  for (page_group in page_groups) {
    page_df <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_dir <- file.path(out_dir, safe_id(page_group))
    dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)

    for (fixed_value in fixed_levels) {
      unit_plot <- make_endpoint_unit_plot(
        endpoint_plot = page_df,
        fixed_col = fixed_col,
        fixed_value = fixed_value,
        vary_col = vary_col,
        fixed_axis_label = fixed_axis_label,
        vary_axis_label = vary_axis_label,
        burden_limits = burden_limits
      )
      if (is.null(unit_plot)) next

      if (!identical(page_group, "all")) {
        unit_plot <- unit_plot + plot_annotation(title = page_group)
      }

      unit_pdf <- file.path(page_dir, paste0(safe_id("unit", fixed_axis_label, fixed_value), ".pdf"))
      ggsave(unit_pdf, unit_plot, width = 7.4, height = if (identical(page_group, "all")) 4.4 else 4.7)
      manifest_rows[[row_i]] <- data.frame(
        trigger = page_group,
        fixed_parameter = fixed_axis_label,
        fixed_value = fixed_value,
        file = unit_pdf,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }

  if (length(manifest_rows)) {
    write_tsv(bind_rows(manifest_rows), file.path(out_dir, "unit_files.tsv"))
  }
  invisible(out_dir)
}

map_gradient_colors <- function(x,
                                colors,
                                limits = NULL,
                                values = NULL,
                                na_color = "grey85") {
  x <- as.numeric(x)
  out <- rep(na_color, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)

  if (is.null(limits) || length(limits) != 2L || any(!is.finite(limits))) {
    limits <- range(x[ok], na.rm = TRUE)
  }
  limits <- as.numeric(limits)
  if (!is.finite(diff(limits)) || diff(limits) <= 0) {
    limits <- limits[[1]] + c(-0.5, 0.5)
  }

  scaled <- (x - limits[[1]]) / diff(limits)
  scaled <- pmin(1, pmax(0, scaled))
  pal <- scales::gradient_n_pal(colors, values = values)
  out[ok] <- pal(scaled[ok])
  out
}

make_triangle_polygon_rows <- function(tile_df,
                                       x_col,
                                       y_col,
                                       fill_col,
                                       metric_label,
                                       upper = TRUE) {
  pieces <- vector("list", nrow(tile_df))
  for (i in seq_len(nrow(tile_df))) {
    x <- as.numeric(tile_df[[x_col]][[i]])
    y <- as.numeric(tile_df[[y_col]][[i]])
    if (isTRUE(upper)) {
      px <- c(x - 0.5, x + 0.5, x - 0.5)
      py <- c(y + 0.5, y + 0.5, y - 0.5)
    } else {
      px <- c(x + 0.5, x + 0.5, x - 0.5)
      py <- c(y + 0.5, y - 0.5, y - 0.5)
    }
    pieces[[i]] <- data.frame(
      triangle_id = paste(metric_label, i, sep = "_"),
      x = px,
      y = py,
      fill_color = tile_df[[fill_col]][[i]],
      stringsAsFactors = FALSE
    )
  }
  bind_rows(pieces)
}

make_endpoint_triangle_data <- function(sub,
                                        vary_col,
                                        burden_limits = NULL) {
  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  o2_levels <- levels(droplevels(sub$o2_label))
  if (is.null(o2_levels) || !length(o2_levels)) o2_levels <- sort(unique(as.character(sub$o2_label)))

  sub$vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)
  sub$o2_label <- factor(as.character(sub$o2_label), levels = o2_levels)
  sub$.tile_x <- as.integer(sub$o2_label)
  sub$.tile_y <- as.integer(sub$vary_label)
  sub$.burden_fill <- map_gradient_colors(
    sub$pred_log10_burden_cells,
    colors = c("#FFFFFF", "#FFF7BC", "#FEC44F", "#FB6A4A", "#CC2C7A", "#4A1486"),
    limits = burden_limits
  )
  sub$.ploidy_fill <- map_gradient_colors(
    sub$mean_ploidy,
    colors = c("#FFFFFF", "#F3E8FF", "#C084FC", "#7C3AED", "#1D4ED8"),
    limits = c(0, 6),
    values = scales::rescale(c(0, 1.5, 3, 4.5, 6))
  )

  upper <- make_triangle_polygon_rows(sub, ".tile_x", ".tile_y", ".burden_fill", "burden", upper = TRUE)
  lower <- make_triangle_polygon_rows(sub, ".tile_x", ".tile_y", ".ploidy_fill", "ploidy", upper = FALSE)
  list(
    polygons = bind_rows(upper, lower),
    o2_levels = o2_levels,
    vary_levels = vary_levels
  )
}

plot_gradient_legend_identity <- function(colors,
                                          limits,
                                          title,
                                          values = NULL,
                                          breaks = NULL) {
  if (is.null(limits) || length(limits) != 2L || any(!is.finite(limits)) || diff(limits) <= 0) {
    limits <- c(0, 1)
  }
  vals <- seq(limits[[1]], limits[[2]], length.out = 160)
  df <- data.frame(
    x = 1,
    y = vals,
    fill_color = map_gradient_colors(vals, colors = colors, limits = limits, values = values),
    stringsAsFactors = FALSE
  )
  breaks <- breaks %||% pretty(limits, n = 4)
  ggplot(df, aes(x, y, fill = fill_color)) +
    geom_raster() +
    scale_fill_identity() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = breaks, position = "right", expand = c(0, 0)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 7) +
    theme(
      axis.text.y = element_text(size = 6),
      axis.ticks.y = element_line(linewidth = 0.2),
      plot.title = element_text(size = 7, hjust = 0)
    )
}

make_triangle_legend_plot <- function(burden_limits = NULL) {
  burden_legend <- plot_gradient_legend_identity(
    colors = c("#FFFFFF", "#FFF7BC", "#FEC44F", "#FB6A4A", "#CC2C7A", "#4A1486"),
    limits = burden_limits,
    title = "upper\nlog10 burden"
  )
  ploidy_legend <- plot_gradient_legend_identity(
    colors = c("#FFFFFF", "#F3E8FF", "#C084FC", "#7C3AED", "#1D4ED8"),
    limits = c(0, 6),
    values = scales::rescale(c(0, 1.5, 3, 4.5, 6)),
    breaks = 0:6,
    title = "lower\nmean ploidy"
  )
  wrap_plots(burden_legend, ploidy_legend, ncol = 1)
}

make_endpoint_triangle_unit_plot <- function(endpoint_plot,
                                             fixed_col,
                                             fixed_value,
                                             vary_col,
                                             fixed_axis_label,
                                             vary_axis_label,
                                             burden_limits = NULL) {
  fixed_chr <- as.character(fixed_value)
  sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  tri <- make_endpoint_triangle_data(sub, vary_col = vary_col, burden_limits = burden_limits)
  unit_title <- paste(strwrap(paste0(fixed_axis_label, " = ", fixed_chr), width = 34), collapse = "\n")

  ggplot(tri$polygons, aes(x, y, group = triangle_id, fill = fill_color)) +
    geom_polygon(color = "white", linewidth = 0.08) +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = seq_along(tri$o2_levels),
      labels = tri$o2_levels,
      limits = c(0.5, length(tri$o2_levels) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(tri$vary_levels),
      labels = tri$vary_levels,
      limits = c(0.5, length(tri$vary_levels) + 0.5),
      expand = c(0, 0)
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(x = "O2_S0", y = vary_axis_label, title = unit_title) +
    theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      plot.title = element_text(size = 7, hjust = 0.5)
    )
}

make_endpoint_triangle_unit_grid_plot <- function(endpoint_plot,
                                                  fixed_col,
                                                  vary_col,
                                                  fixed_axis_label,
                                                  vary_axis_label,
                                                  fixed_levels,
                                                  burden_limits = NULL,
                                                  page_title = NULL) {
  unit_plots <- lapply(fixed_levels, function(fixed_value) {
    make_endpoint_triangle_unit_plot(
      endpoint_plot = endpoint_plot,
      fixed_col = fixed_col,
      fixed_value = fixed_value,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      burden_limits = burden_limits
    )
  })
  unit_plots <- Filter(Negate(is.null), unit_plots)
  if (!length(unit_plots)) return(NULL)

  grid_plot <- wrap_plots(unit_plots, ncol = 3) +
    plot_annotation(title = page_title)
  (grid_plot | make_triangle_legend_plot(burden_limits)) +
    plot_layout(widths = c(1, 0.08))
}

save_endpoint_triangle_unit_grid_pages <- function(endpoint_plot,
                                                   fixed_col,
                                                   vary_col,
                                                   fixed_axis_label,
                                                   vary_axis_label,
                                                   fixed_levels,
                                                   out_pdf,
                                                   include_trigger = TRUE) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  burden_limits <- range(endpoint_plot$pred_log10_burden_cells, na.rm = TRUE)
  if (any(!is.finite(burden_limits))) burden_limits <- NULL

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  grDevices::pdf(out_pdf, width = 16, height = 22, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_group in page_groups) {
    page_df <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_title <- if (identical(page_group, "all")) {
      paste0(fixed_axis_label, " endpoint triangle unit grid")
    } else {
      paste0(fixed_axis_label, " endpoint triangle unit grid / ", page_group)
    }
    page_plot <- make_endpoint_triangle_unit_grid_plot(
      endpoint_plot = page_df,
      fixed_col = fixed_col,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      fixed_levels = fixed_levels,
      burden_limits = burden_limits,
      page_title = page_title
    )
    if (!is.null(page_plot)) print(page_plot)
  }
  invisible(out_pdf)
}

save_endpoint_triangle_unit_individual_plots <- function(endpoint_plot,
                                                         fixed_col,
                                                         vary_col,
                                                         fixed_axis_label,
                                                         vary_axis_label,
                                                         fixed_levels,
                                                         out_dir,
                                                         include_trigger = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  burden_limits <- range(endpoint_plot$pred_log10_burden_cells, na.rm = TRUE)
  if (any(!is.finite(burden_limits))) burden_limits <- NULL

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  manifest_rows <- list()
  row_i <- 1L
  for (page_group in page_groups) {
    page_df <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_dir <- file.path(out_dir, safe_id(page_group))
    dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)

    for (fixed_value in fixed_levels) {
      unit_plot <- make_endpoint_triangle_unit_plot(
        endpoint_plot = page_df,
        fixed_col = fixed_col,
        fixed_value = fixed_value,
        vary_col = vary_col,
        fixed_axis_label = fixed_axis_label,
        vary_axis_label = vary_axis_label,
        burden_limits = burden_limits
      )
      if (is.null(unit_plot)) next

      unit_plot <- (unit_plot | make_triangle_legend_plot(burden_limits)) +
        plot_layout(widths = c(1, 0.14))
      if (!identical(page_group, "all")) {
        unit_plot <- unit_plot + plot_annotation(title = page_group)
      }

      unit_pdf <- file.path(page_dir, paste0(safe_id("unit_triangle", fixed_axis_label, fixed_value), ".pdf"))
      ggsave(unit_pdf, unit_plot, width = 6.2, height = if (identical(page_group, "all")) 6.2 else 6.5)
      manifest_rows[[row_i]] <- data.frame(
        trigger = page_group,
        fixed_parameter = fixed_axis_label,
        fixed_value = fixed_value,
        file = unit_pdf,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }

  if (length(manifest_rows)) {
    write_tsv(bind_rows(manifest_rows), file.path(out_dir, "unit_files.tsv"))
  }
  invisible(out_dir)
}

make_timecourse_unit_treatment_lines <- function(endpoint_sub,
                                                 max_day,
                                                 cell_width = 2,
                                                 cell_height = 1,
                                                 show_treatment_lines = TRUE) {
  if (!isTRUE(show_treatment_lines) || !nrow(endpoint_sub)) return(endpoint_sub[FALSE, , drop = FALSE])
  required_cols <- c("trigger_day", "unit_o2_index", "unit_vary_index")
  if (!all(required_cols %in% names(endpoint_sub))) return(endpoint_sub[FALSE, , drop = FALSE])

  max_day <- as.numeric(max_day)
  if (!is.finite(max_day) || max_day < 0) max_day <- 0
  denom_day <- max(1, max_day)
  endpoint_sub <- endpoint_sub[is.finite(as.numeric(endpoint_sub$trigger_day)), , drop = FALSE]
  if (!nrow(endpoint_sub)) return(endpoint_sub)

  rows <- list()
  row_i <- 1L
  for (i in seq_len(nrow(endpoint_sub))) {
    dr <- endpoint_sub[i, , drop = FALSE]
    treatment_days <- as.numeric(dr$trigger_day)
    is_intermittent <- "intermittent_on_day" %in% names(dr) &&
      "intermittent_off_day" %in% names(dr) &&
      "post_treatment_duration_day" %in% names(dr) &&
      is.finite(as.numeric(dr$intermittent_on_day)) &&
      is.finite(as.numeric(dr$intermittent_off_day)) &&
      is.finite(as.numeric(dr$post_treatment_duration_day)) &&
      as.numeric(dr$intermittent_on_day) > 0 &&
      as.numeric(dr$intermittent_off_day) >= 0

    if (isTRUE(is_intermittent)) {
      period <- as.numeric(dr$intermittent_on_day) + as.numeric(dr$intermittent_off_day)
      if (period > 0) {
        n_starts <- max(1L, as.integer(ceiling((as.numeric(dr$post_treatment_duration_day) - 1e-9) / period)))
        treatment_days <- as.numeric(dr$trigger_day) + seq(0, by = period, length.out = n_starts)
      }
    }

    treatment_days <- treatment_days[is.finite(treatment_days) & treatment_days >= 0 & treatment_days <= max_day]
    if (!length(treatment_days)) next

    line_rows <- dr[rep(1, length(treatment_days)), , drop = FALSE]
    line_rows$treatment_day <- treatment_days
    line_rows$x0 <- (as.numeric(line_rows$unit_o2_index) - 1) * cell_width +
      (pmin(max_day, pmax(0, treatment_days)) / denom_day) * cell_width
    line_rows$x1 <- line_rows$x0
    line_rows$y0 <- (as.numeric(line_rows$unit_vary_index) - 1) * cell_height
    line_rows$y1 <- as.numeric(line_rows$unit_vary_index) * cell_height
    rows[[row_i]] <- line_rows
    row_i <- row_i + 1L
  }

  if (!length(rows)) return(endpoint_sub[FALSE, , drop = FALSE])
  bind_rows(rows)
}

plot_timecourse_unit_metric <- function(plot_df,
                                        endpoint_plot,
                                        fixed_col,
                                        fixed_value,
                                        vary_col,
                                        value_col,
                                        fill_kind = c("burden", "ploidy"),
                                        y_axis_label = NULL,
                                        show_y_axis = TRUE,
                                        max_day,
                                        show_treatment_lines = TRUE,
                                        burden_limits = NULL,
                                        cell_width = 2,
                                        cell_height = 1) {
  fill_kind <- match.arg(fill_kind)
  fixed_chr <- as.character(fixed_value)
  sub <- plot_df[as.character(plot_df[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (!nrow(sub)) return(NULL)

  vary_levels <- levels(droplevels(sub[[vary_col]]))
  if (is.null(vary_levels) || !length(vary_levels)) vary_levels <- sort(unique(as.character(sub[[vary_col]])))
  o2_levels <- levels(droplevels(sub$o2_label))
  if (is.null(o2_levels) || !length(o2_levels)) o2_levels <- sort(unique(as.character(sub$o2_label)))

  sub$unit_vary_label <- factor(as.character(sub[[vary_col]]), levels = vary_levels)
  sub$unit_vary_index <- as.integer(sub$unit_vary_label)
  sub$unit_o2_label <- factor(as.character(sub$o2_label), levels = o2_levels)
  sub$unit_o2_index <- as.integer(sub$unit_o2_label)

  max_day_i <- ceiling(as.numeric(max_day))
  if (!is.finite(max_day_i) || max_day_i < 0) max_day_i <- 0
  day_count <- max_day_i + 1L
  sub_full <- complete_timecourse_for_plot(sub, value_col, fill_value = 0, max_day = max_day_i)
  sub_full$.x <- (as.numeric(sub_full$unit_o2_index) - 1) * cell_width +
    (as.numeric(sub_full$day) + 0.5) * (cell_width / day_count)
  sub_full$.y <- (as.numeric(sub_full$unit_vary_index) - 0.5) * cell_height

  endpoint_sub <- endpoint_plot[as.character(endpoint_plot[[fixed_col]]) == fixed_chr, , drop = FALSE]
  if (nrow(endpoint_sub)) {
    endpoint_sub$unit_vary_label <- factor(as.character(endpoint_sub[[vary_col]]), levels = vary_levels)
    endpoint_sub$unit_vary_index <- as.integer(endpoint_sub$unit_vary_label)
    endpoint_sub$unit_o2_label <- factor(as.character(endpoint_sub$o2_label), levels = o2_levels)
    endpoint_sub$unit_o2_index <- as.integer(endpoint_sub$unit_o2_label)
  }
  treatment_lines <- make_timecourse_unit_treatment_lines(
    endpoint_sub,
    max_day = max_day_i,
    cell_width = cell_width,
    cell_height = cell_height,
    show_treatment_lines = show_treatment_lines
  )

  border <- expand.grid(
    unit_o2_index = seq_along(o2_levels),
    unit_vary_index = seq_along(vary_levels),
    KEEP.OUT.ATTRS = FALSE
  )
  border$xmin <- (border$unit_o2_index - 1) * cell_width
  border$xmax <- border$unit_o2_index * cell_width
  border$ymin <- (border$unit_vary_index - 1) * cell_height
  border$ymax <- border$unit_vary_index * cell_height

  y_axis_label <- y_axis_label %||% vary_col
  p <- ggplot(sub_full, aes(.x, .y, fill = .data[[value_col]])) +
    geom_tile(width = cell_width / day_count, height = cell_height) +
    geom_rect(
      data = border,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      color = "grey80",
      linewidth = 0.08
    ) +
    scale_x_continuous(
      breaks = (seq_along(o2_levels) - 0.5) * cell_width,
      labels = o2_levels,
      limits = c(0, length(o2_levels) * cell_width),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = (seq_along(vary_levels) - 0.5) * cell_height,
      labels = vary_levels,
      limits = c(0, length(vary_levels) * cell_height),
      expand = c(0, 0)
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(x = "O2_S0", y = if (show_y_axis) y_axis_label else NULL) +
    theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 7, hjust = 0.5)
    )
  if (!show_y_axis) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  if (identical(fill_kind, "burden")) {
    p <- p + scale_fill_log10_burden("Timecourse\nlog10 burden", limits = burden_limits)
  } else {
    p <- p + scale_fill_mean_ploidy("Timecourse\nmean ploidy")
  }
  if (nrow(treatment_lines) > 0L) {
    p <- p + geom_segment(
      data = treatment_lines,
      aes(x = x0, xend = x1, y = y0, yend = y1),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.18
    )
  }
  p
}

make_timecourse_unit_plot <- function(burden_plot,
                                      ploidy_plot,
                                      endpoint_plot,
                                      fixed_col,
                                      fixed_value,
                                      vary_col,
                                      fixed_axis_label,
                                      vary_axis_label,
                                      max_day,
                                      show_treatment_lines = TRUE,
                                      burden_limits = NULL) {
  fixed_chr <- as.character(fixed_value)
  if (!any(as.character(burden_plot[[fixed_col]]) == fixed_chr)) return(NULL)
  unit_title <- paste(strwrap(paste0(fixed_axis_label, " = ", fixed_chr), width = 34), collapse = "\n")

  p_burden <- plot_timecourse_unit_metric(
    plot_df = burden_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = fixed_col,
    fixed_value = fixed_value,
    vary_col = vary_col,
    value_col = "pred_log10_burden_cells",
    fill_kind = "burden",
    y_axis_label = vary_axis_label,
    show_y_axis = TRUE,
    max_day = max_day,
    show_treatment_lines = show_treatment_lines,
    burden_limits = burden_limits
  ) +
    ggtitle(paste("burden", unit_title, sep = "\n"))
  p_ploidy <- plot_timecourse_unit_metric(
    plot_df = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = fixed_col,
    fixed_value = fixed_value,
    vary_col = vary_col,
    value_col = "mean_ploidy",
    fill_kind = "ploidy",
    y_axis_label = vary_axis_label,
    show_y_axis = FALSE,
    max_day = max_day,
    show_treatment_lines = show_treatment_lines,
    burden_limits = burden_limits
  ) +
    ggtitle("ploidy")

  wrap_plots(p_burden, p_ploidy, ncol = 2)
}

make_timecourse_unit_grid_plot <- function(burden_plot,
                                           ploidy_plot,
                                           endpoint_plot,
                                           fixed_col,
                                           vary_col,
                                           fixed_axis_label,
                                           vary_axis_label,
                                           fixed_levels,
                                           max_day,
                                           show_treatment_lines = TRUE,
                                           burden_limits = NULL,
                                           page_title = NULL) {
  unit_plots <- lapply(fixed_levels, function(fixed_value) {
    make_timecourse_unit_plot(
      burden_plot = burden_plot,
      ploidy_plot = ploidy_plot,
      endpoint_plot = endpoint_plot,
      fixed_col = fixed_col,
      fixed_value = fixed_value,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      max_day = max_day,
      show_treatment_lines = show_treatment_lines,
      burden_limits = burden_limits
    )
  })
  unit_plots <- Filter(Negate(is.null), unit_plots)
  if (!length(unit_plots)) return(NULL)

  wrap_plots(unit_plots, ncol = 3, guides = "collect") +
    plot_layout(guides = "collect") +
    plot_annotation(title = page_title)
}

save_timecourse_unit_grid_pages <- function(burden_plot,
                                            ploidy_plot,
                                            endpoint_plot,
                                            fixed_col,
                                            vary_col,
                                            fixed_axis_label,
                                            vary_axis_label,
                                            fixed_levels,
                                            out_pdf,
                                            max_day,
                                            include_trigger = TRUE,
                                            show_treatment_lines = TRUE) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  burden_max <- max(burden_plot$pred_log10_burden_cells, na.rm = TRUE)
  burden_limits <- if (is.finite(burden_max) && burden_max > 0) c(0, burden_max) else c(0, 1)

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  grDevices::pdf(out_pdf, width = 26, height = 18, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_group in page_groups) {
    page_endpoint <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_burden <- if (identical(page_group, "all")) {
      burden_plot
    } else {
      burden_plot[as.character(burden_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_ploidy <- if (identical(page_group, "all")) {
      ploidy_plot
    } else {
      ploidy_plot[as.character(ploidy_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_title <- if (identical(page_group, "all")) {
      paste0(fixed_axis_label, " timecourse unit grid / day 0-", max_day)
    } else {
      paste0(fixed_axis_label, " timecourse unit grid / ", page_group, " / day 0-", max_day)
    }
    page_plot <- make_timecourse_unit_grid_plot(
      burden_plot = page_burden,
      ploidy_plot = page_ploidy,
      endpoint_plot = page_endpoint,
      fixed_col = fixed_col,
      vary_col = vary_col,
      fixed_axis_label = fixed_axis_label,
      vary_axis_label = vary_axis_label,
      fixed_levels = fixed_levels,
      max_day = max_day,
      show_treatment_lines = show_treatment_lines,
      burden_limits = burden_limits,
      page_title = page_title
    )
    if (!is.null(page_plot)) print(page_plot)
  }
  invisible(out_pdf)
}

save_timecourse_unit_individual_plots <- function(burden_plot,
                                                  ploidy_plot,
                                                  endpoint_plot,
                                                  fixed_col,
                                                  vary_col,
                                                  fixed_axis_label,
                                                  vary_axis_label,
                                                  fixed_levels,
                                                  out_dir,
                                                  max_day,
                                                  include_trigger = TRUE,
                                                  show_treatment_lines = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  burden_max <- max(burden_plot$pred_log10_burden_cells, na.rm = TRUE)
  burden_limits <- if (is.finite(burden_max) && burden_max > 0) c(0, burden_max) else c(0, 1)

  page_groups <- if (include_trigger && "trigger_label" %in% names(endpoint_plot)) {
    levels(droplevels(endpoint_plot$trigger_label))
  } else {
    "all"
  }

  manifest_rows <- list()
  row_i <- 1L
  for (page_group in page_groups) {
    page_endpoint <- if (identical(page_group, "all")) {
      endpoint_plot
    } else {
      endpoint_plot[as.character(endpoint_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_burden <- if (identical(page_group, "all")) {
      burden_plot
    } else {
      burden_plot[as.character(burden_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_ploidy <- if (identical(page_group, "all")) {
      ploidy_plot
    } else {
      ploidy_plot[as.character(ploidy_plot$trigger_label) == page_group, , drop = FALSE]
    }
    page_dir <- file.path(out_dir, safe_id(page_group))
    dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)

    for (fixed_value in fixed_levels) {
      unit_plot <- make_timecourse_unit_plot(
        burden_plot = page_burden,
        ploidy_plot = page_ploidy,
        endpoint_plot = page_endpoint,
        fixed_col = fixed_col,
        fixed_value = fixed_value,
        vary_col = vary_col,
        fixed_axis_label = fixed_axis_label,
        vary_axis_label = vary_axis_label,
        max_day = max_day,
        show_treatment_lines = show_treatment_lines,
        burden_limits = burden_limits
      )
      if (is.null(unit_plot)) next

      if (!identical(page_group, "all")) {
        unit_plot <- unit_plot + plot_annotation(title = page_group)
      }

      unit_pdf <- file.path(page_dir, paste0(safe_id("unit_timecourse", fixed_axis_label, fixed_value), ".pdf"))
      ggsave(unit_pdf, unit_plot, width = 12, height = if (identical(page_group, "all")) 5.2 else 5.5)
      manifest_rows[[row_i]] <- data.frame(
        trigger = page_group,
        fixed_parameter = fixed_axis_label,
        fixed_value = fixed_value,
        max_day = max_day,
        file = unit_pdf,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }

  if (length(manifest_rows)) {
    write_tsv(bind_rows(manifest_rows), file.path(out_dir, "unit_files.tsv"))
  }
  invisible(out_dir)
}

plot_timecourse_unit_grid_heatmaps <- function(endpoint_summary,
                                               burden_all,
                                               ploidy_summary,
                                               run_params,
                                               out_dir,
                                               include_trigger = TRUE,
                                               show_treatment_lines = TRUE,
                                               condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  burden_plot <- add_plot_labels(burden_all, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  ploidy_plot <- add_plot_labels(ploidy_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  max_day <- ceiling(max(c(burden_plot$day, ploidy_plot$day), na.rm = TRUE))
  if (!is.finite(max_day) || max_day < 0) max_day <- 0

  unit_root <- file.path(out_dir, "timecourse_unit_grid")
  reset_generated_dir(unit_root)
  write_tsv(data.frame(time_start_day = 0, time_end_day = max_day), file.path(unit_root, "time_window.tsv"))

  fixed_pmiss_dir <- file.path(unit_root, "fixed_p_mis_base")
  dir.create(fixed_pmiss_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    data.frame(fixed_parameter = "p_mis_base", fixed_value = levels(endpoint_plot$pmiss_label)),
    file.path(fixed_pmiss_dir, "fixed_values.tsv")
  )
  save_timecourse_unit_grid_pages(
    burden_plot = burden_plot,
    ploidy_plot = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_pdf = file.path(fixed_pmiss_dir, "timecourse_burden_ploidy_by_o2_grid.pdf"),
    max_day = max_day,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines
  )
  save_timecourse_unit_individual_plots(
    burden_plot = burden_plot,
    ploidy_plot = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_dir = file.path(fixed_pmiss_dir, "units"),
    max_day = max_day,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines
  )

  fixed_pwgd_dir <- file.path(unit_root, "fixed_p_wgd")
  dir.create(fixed_pwgd_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    data.frame(fixed_parameter = "p_wgd", fixed_value = levels(endpoint_plot$pwgd_label)),
    file.path(fixed_pwgd_dir, "fixed_values.tsv")
  )
  save_timecourse_unit_grid_pages(
    burden_plot = burden_plot,
    ploidy_plot = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_pdf = file.path(fixed_pwgd_dir, "timecourse_burden_ploidy_by_o2_grid.pdf"),
    max_day = max_day,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines
  )
  save_timecourse_unit_individual_plots(
    burden_plot = burden_plot,
    ploidy_plot = ploidy_plot,
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_dir = file.path(fixed_pwgd_dir, "units"),
    max_day = max_day,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines
  )
  invisible(unit_root)
}

plot_endpoint_unit_grid_heatmaps <- function(endpoint_summary,
                                             run_params,
                                             out_dir,
                                             include_trigger = TRUE,
                                             condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  unit_root <- file.path(out_dir, "endpoint_unit_grid")
  reset_generated_dir(unit_root)

  fixed_pmiss_dir <- file.path(unit_root, "fixed_p_mis_base")
  dir.create(fixed_pmiss_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    data.frame(fixed_parameter = "p_mis_base", fixed_value = levels(endpoint_plot$pmiss_label)),
    file.path(fixed_pmiss_dir, "fixed_values.tsv")
  )
  save_endpoint_unit_grid_pages(
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_pdf = file.path(fixed_pmiss_dir, "endpoint_burden_ploidy_by_o2_grid.pdf"),
    include_trigger = include_trigger
  )
  save_endpoint_triangle_unit_grid_pages(
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_pdf = file.path(fixed_pmiss_dir, "endpoint_burden_ploidy_triangle_by_o2_grid.pdf"),
    include_trigger = include_trigger
  )
  save_endpoint_unit_individual_plots(
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_dir = file.path(fixed_pmiss_dir, "units"),
    include_trigger = include_trigger
  )
  save_endpoint_triangle_unit_individual_plots(
    endpoint_plot = endpoint_plot,
    fixed_col = "pmiss_label",
    vary_col = "pwgd_label",
    fixed_axis_label = "p_mis_base",
    vary_axis_label = "p_wgd relative to fitted (actual p_wgd)",
    fixed_levels = levels(endpoint_plot$pmiss_label),
    out_dir = file.path(fixed_pmiss_dir, "triangle_units"),
    include_trigger = include_trigger
  )

  fixed_pwgd_dir <- file.path(unit_root, "fixed_p_wgd")
  dir.create(fixed_pwgd_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    data.frame(fixed_parameter = "p_wgd", fixed_value = levels(endpoint_plot$pwgd_label)),
    file.path(fixed_pwgd_dir, "fixed_values.tsv")
  )
  save_endpoint_unit_grid_pages(
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_pdf = file.path(fixed_pwgd_dir, "endpoint_burden_ploidy_by_o2_grid.pdf"),
    include_trigger = include_trigger
  )
  save_endpoint_triangle_unit_grid_pages(
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_pdf = file.path(fixed_pwgd_dir, "endpoint_burden_ploidy_triangle_by_o2_grid.pdf"),
    include_trigger = include_trigger
  )
  save_endpoint_unit_individual_plots(
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_dir = file.path(fixed_pwgd_dir, "units"),
    include_trigger = include_trigger
  )
  save_endpoint_triangle_unit_individual_plots(
    endpoint_plot = endpoint_plot,
    fixed_col = "pwgd_label",
    vary_col = "pmiss_label",
    fixed_axis_label = "p_wgd",
    vary_axis_label = "p_mis_base",
    fixed_levels = levels(endpoint_plot$pwgd_label),
    out_dir = file.path(fixed_pwgd_dir, "triangle_units"),
    include_trigger = include_trigger
  )
  invisible(unit_root)
}

plot_interaction_heatmaps <- function(endpoint_summary,
                                      burden_all,
                                      ploidy_summary,
                                      run_params,
                                      out_dir,
                                      include_trigger = TRUE,
                                      show_treatment_lines = TRUE,
                                      pmiss_axis_label = "post p_mis_base",
                                      condition_pmiss_prefix = "post p_miss=") {
  endpoint_plot <- add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  ref_group_cols <- if (include_trigger) c("o2_S0", "trigger_burden_cells") else "o2_S0"
  ref <- endpoint_plot %>%
    group_by(across(all_of(ref_group_cols))) %>%
    filter(
      .plot_p_mis_base == min(.plot_p_mis_base, na.rm = TRUE),
      abs(.plot_p_wgd - as.numeric(run_params$p_wgd)) == min(abs(.plot_p_wgd - as.numeric(run_params$p_wgd)), na.rm = TRUE)
    ) %>%
    slice(1) %>%
    ungroup() %>%
    select(
      all_of(ref_group_cols),
      ref_log10_burden = pred_log10_burden_cells,
      ref_mean_ploidy = mean_ploidy
    )
  endpoint_plot <- endpoint_plot %>%
    left_join(ref, by = ref_group_cols) %>%
    mutate(
      delta_log10_burden = pred_log10_burden_cells - ref_log10_burden,
      delta_mean_ploidy = mean_ploidy - ref_mean_ploidy
    )

  heat_theme <- theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 7)
    )
  facet_layer <- if (include_trigger) facet_grid(trigger_label ~ o2_label) else facet_grid(. ~ o2_label)
  endpoint_height <- if (include_trigger) 12 else 3.6
  timecourse_height <- if (include_trigger) 13 else 5.2

  p <- ggplot(endpoint_plot, aes(pmiss_label, pwgd_label, fill = pred_log10_burden_cells)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_layer +
    scale_fill_log10_burden("Endpoint\nlog10 burden") +
    labs(x = pmiss_axis_label, y = "p_wgd relative to fitted (actual p_wgd)") +
    heat_theme
  ggsave(file.path(out_dir, "endpoint_log10_burden_heatmap.pdf"), p, width = 18, height = endpoint_height)

  p <- ggplot(endpoint_plot, aes(pmiss_label, pwgd_label, fill = mean_ploidy)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_layer +
    scale_fill_mean_ploidy("Endpoint\nmean ploidy") +
    labs(x = pmiss_axis_label, y = "p_wgd relative to fitted (actual p_wgd)", fill = "Endpoint\nmean ploidy") +
    heat_theme
  ggsave(file.path(out_dir, "endpoint_mean_ploidy_heatmap.pdf"), p, width = 18, height = endpoint_height)

  p <- ggplot(endpoint_plot, aes(pmiss_label, pwgd_label, fill = delta_log10_burden)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_layer +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") +
    labs(x = pmiss_axis_label, y = "p_wgd relative to fitted (actual p_wgd)", fill = "Delta log10\nburden") +
    heat_theme
  ggsave(file.path(out_dir, "endpoint_delta_log10_burden_heatmap.pdf"), p, width = 18, height = endpoint_height)

  p <- ggplot(endpoint_plot, aes(pmiss_label, pwgd_label, fill = delta_mean_ploidy)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_layer +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") +
    labs(x = pmiss_axis_label, y = "p_wgd relative to fitted (actual p_wgd)", fill = "Delta\nmean ploidy") +
    heat_theme
  ggsave(file.path(out_dir, "endpoint_delta_mean_ploidy_heatmap.pdf"), p, width = 18, height = endpoint_height)

  burden_plot <- add_plot_labels(burden_all, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  burden_plot_full <- complete_timecourse_for_plot(burden_plot, "pred_log10_burden_cells", fill_value = 0)
  treatment_lines <- make_treatment_line_segments(
    add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix),
    row_index_col = "condition_index",
    show_treatment_lines = show_treatment_lines
  )
  condition_breaks <- seq_along(levels(burden_plot$condition_label))
  condition_labels <- levels(burden_plot$condition_label)
  p <- ggplot(burden_plot_full, aes(day, condition_index, fill = pred_log10_burden_cells)) +
    geom_raster() +
    facet_layer +
    scale_fill_log10_burden("log10 burden") +
    scale_y_continuous(breaks = condition_breaks, labels = condition_labels, expand = c(0, 0)) +
    labs(x = "Day", y = "p_wgd relative to fitted (actual p_wgd)") +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 5),
      axis.ticks.y = element_line(linewidth = 0.2),
      strip.text = element_text(size = 7)
    )
  if (nrow(treatment_lines) > 0L) {
    p <- p + geom_segment(
      data = treatment_lines,
      aes(x = treatment_day, xend = treatment_day, y = y0, yend = y1),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.25
    )
  }
  ggsave(file.path(out_dir, "timecourse_log10_burden_heatmap.pdf"), p, width = 20, height = timecourse_height)

  ploidy_plot <- add_plot_labels(ploidy_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix)
  ploidy_plot_full <- complete_timecourse_for_plot(ploidy_plot, "mean_ploidy", fill_value = 0)
  treatment_lines <- make_treatment_line_segments(
    add_plot_labels(endpoint_summary, run_params, condition_pmiss_prefix = condition_pmiss_prefix),
    row_index_col = "condition_index",
    show_treatment_lines = show_treatment_lines
  )
  condition_breaks <- seq_along(levels(ploidy_plot$condition_label))
  condition_labels <- levels(ploidy_plot$condition_label)
  p <- ggplot(ploidy_plot_full, aes(day, condition_index, fill = mean_ploidy)) +
    geom_raster() +
    facet_layer +
    scale_fill_mean_ploidy("mean ploidy") +
    scale_y_continuous(breaks = condition_breaks, labels = condition_labels, expand = c(0, 0)) +
    labs(x = "Day", y = "p_wgd relative to fitted (actual p_wgd)", fill = "mean ploidy") +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 5),
      axis.ticks.y = element_line(linewidth = 0.2),
      strip.text = element_text(size = 7)
    )
  if (nrow(treatment_lines) > 0L) {
    p <- p + geom_segment(
      data = treatment_lines,
      aes(x = treatment_day, xend = treatment_day, y = y0, yend = y1),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed",
      linewidth = 0.25
    )
  }
  ggsave(file.path(out_dir, "timecourse_mean_ploidy_heatmap.pdf"), p, width = 20, height = timecourse_height)

  plot_split_timecourse_heatmaps(
    endpoint_summary = endpoint_summary,
    burden_all = burden_all,
    ploidy_summary = ploidy_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines,
    condition_pmiss_prefix = condition_pmiss_prefix
  )
  plot_split_endpoint_heatmaps(
    endpoint_summary = endpoint_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    condition_pmiss_prefix = condition_pmiss_prefix
  )
  plot_endpoint_unit_grid_heatmaps(
    endpoint_summary = endpoint_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    condition_pmiss_prefix = condition_pmiss_prefix
  )
  plot_timecourse_unit_grid_heatmaps(
    endpoint_summary = endpoint_summary,
    burden_all = burden_all,
    ploidy_summary = ploidy_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    show_treatment_lines = show_treatment_lines,
    condition_pmiss_prefix = condition_pmiss_prefix
  )
}

first_finite_value <- function(x, label) {
  values <- suppressWarnings(as.numeric(x))
  values <- values[is.finite(values)]
  if (!length(values)) stop("Simulation design does not contain a finite ", label, ".")
  values[[1]]
}

run_factorial_interaction_visualization <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  mode <- normalize_interaction_mode(argv$mode %||% "triggered_drug")
  fit_dir <- resolve_path_value(argv$fit_dir %||% argv$run_dir, getwd())
  if (!is.null(fit_dir)) fit_dir <- normalizePath(fit_dir, mustWork = TRUE)

  simulation_dir <- resolve_path_value(argv$simulation_dir %||% argv$input_dir, getwd())
  analysis_dir <- resolve_path_value(argv$analysis_dir, getwd())
  out_dir <- resolve_path_value(argv$out_dir %||% argv$viz_dir, getwd())
  stem <- interaction_output_stem(mode)
  if (is.null(simulation_dir)) {
    if (is.null(fit_dir)) stop("Provide --simulation_dir or --fit_dir.")
    simulation_dir <- file.path(fit_dir, "simulation", "interactions", stem)
  }
  if (is.null(analysis_dir)) {
    if (is.null(fit_dir)) stop("Provide --analysis_dir when --fit_dir is omitted.")
    analysis_dir <- file.path(fit_dir, "analysis", "interactions", stem)
  }
  if (is.null(out_dir)) {
    if (is.null(fit_dir)) stop("Provide --out_dir when --fit_dir is omitted.")
    out_dir <- file.path(fit_dir, "viz", "interactions", stem)
  }
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)

  validate_artifact_manifest(
    simulation_dir,
    "simulation_manifest.tsv",
    c("interaction_design.tsv", "burden_timecourse.tsv", "ploidy_summary_timecourse.tsv"),
    "Factorial interaction visualization"
  )
  validate_artifact_manifest(
    analysis_dir,
    "analysis_manifest.tsv",
    "endpoint_summary.tsv",
    "Factorial interaction visualization"
  )

  design <- read_required_tsv(file.path(simulation_dir, "interaction_design.tsv"))
  burden_all <- read_required_tsv(file.path(simulation_dir, "burden_timecourse.tsv"))
  ploidy_summary <- read_required_tsv(file.path(simulation_dir, "ploidy_summary_timecourse.tsv"))
  endpoint_summary <- read_required_tsv(file.path(analysis_dir, "endpoint_summary.tsv"))
  run_params <- list(
    p_mis_base = first_finite_value(design$p_mis_base_pre, "p_mis_base_pre"),
    p_wgd = first_finite_value(design$p_wgd_pre, "p_wgd_pre")
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  include_trigger <- !identical(mode, "untreated_factorial")
  plot_interaction_heatmaps(
    endpoint_summary = endpoint_summary,
    burden_all = burden_all,
    ploidy_summary = ploidy_summary,
    run_params = run_params,
    out_dir = out_dir,
    include_trigger = include_trigger,
    show_treatment_lines = include_trigger,
    pmiss_axis_label = if (include_trigger) "post p_mis_base" else "p_mis_base",
    condition_pmiss_prefix = if (include_trigger) "post p_miss=" else "p_miss="
  )

  generated <- list.files(out_dir, recursive = TRUE, full.names = FALSE, all.files = FALSE)
  generated <- generated[generated != "visualization_manifest.tsv"]
  rows <- lapply(generated, function(rel) {
    data.frame(
      artifact = tools::file_path_sans_ext(gsub("[/\\\\]", "_", rel)),
      relative_path = rel,
      role = if (grepl("\\.pdf$", rel, ignore.case = TRUE)) "figure" else "visualization_metadata",
      rows = NA_integer_,
      columns = NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) stop("Visualization did not create any output artifacts.")
  write_artifact_manifest(out_dir, rows, "visualization_manifest.tsv")
  message("Done. Wrote visualization outputs to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  run_factorial_interaction_visualization()
}
