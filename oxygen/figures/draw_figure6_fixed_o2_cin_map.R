#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(patchwork)
})

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in argv) {
    if (!grepl("^--[^=]+=", arg)) next
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    out[[key]] <- sub("^[^=]+=", "", arg)
  }
  out
}

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})

read_wide <- function(path) {
  if (!file.exists(path)) stop("Missing required FixO2 table: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

o2_from_column <- function(x, prefix) {
  key <- sub(prefix, "", x, fixed = TRUE)
  key <- sub("^m", "-", key)
  suppressWarnings(as.numeric(chartr("p", ".", key)))
}

wide_to_points <- function(tab, dataset) {
  prefixes <- c(
    ploidy = "fixo2_eigen_ploidy_o2_",
    p_misseg = "fixo2_population_average_p_misseg_o2_",
    spectral_gap = "fixo2_spectral_gap_o2_"
  )
  cols <- lapply(prefixes, function(prefix) names(tab)[startsWith(names(tab), prefix)])
  counts <- vapply(cols, length, integer(1L))
  if (any(counts != 201L) || length(unique(counts)) != 1L) {
    stop("Expected 201 ploidy, p_misseg, and spectral-gap columns; found ", paste(names(counts), counts, sep = "=", collapse = ", "), call. = FALSE)
  }
  o2_values <- o2_from_column(cols$ploidy, prefixes[["ploidy"]])
  p_o2 <- o2_from_column(cols$p_misseg, prefixes[["p_misseg"]])
  gap_o2 <- o2_from_column(cols$spectral_gap, prefixes[["spectral_gap"]])
  if (any(!is.finite(o2_values)) || !isTRUE(all.equal(o2_values, p_o2)) || !isTRUE(all.equal(o2_values, gap_o2))) {
    stop("The three FixO2 metric groups do not share an identical O2 grid.", call. = FALSE)
  }
  seed <- if ("seed" %in% names(tab)) as.integer(tab$seed) else as.integer(sub("^seed", "", tab$seed_id))
  point_id <- if ("point_id" %in% names(tab)) as.character(tab$point_id) else paste(dataset, "best", paste0("seed", seed), sep = "__")
  matrix_long <- function(col_names) as.numeric(t(as.matrix(tab[, col_names, drop = FALSE])))
  out <- data.frame(
    dataset = dataset,
    seed = rep(seed, each = length(o2_values)),
    point_id = rep(point_id, each = length(o2_values)),
    O2_pct = rep(o2_values, times = nrow(tab)),
    population_average_p_misseg = matrix_long(cols$p_misseg),
    dominant_mean_ploidy = matrix_long(cols$ploidy),
    spectral_gap = matrix_long(cols$spectral_gap),
    stringsAsFactors = FALSE
  )
  out$log10_population_average_p_misseg <- log10(out$population_average_p_misseg)
  out
}

finite_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, names = FALSE, na.rm = TRUE))
}

bin_surface <- function(points, n_y_bins = 60L, min_bin_n = 3L, low_gap_threshold = 0.005) {
  d <- points[
    is.finite(points$O2_pct) &
      is.finite(points$log10_population_average_p_misseg) &
      is.finite(points$dominant_mean_ploidy) &
      is.finite(points$spectral_gap),
    , drop = FALSE
  ]
  if (!nrow(d)) stop("No finite fixed-O2 points are available for binning.", call. = FALSE)
  y_min <- floor(min(d$log10_population_average_p_misseg) * 10) / 10
  y_max <- ceiling(max(d$log10_population_average_p_misseg) * 10) / 10
  if (y_max <= y_min) y_max <- y_min + 0.2
  breaks <- seq(y_min, y_max, length.out = as.integer(n_y_bins) + 1L)
  centers <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  d$y_bin <- findInterval(d$log10_population_average_p_misseg, breaks, all.inside = TRUE)
  split_key <- interaction(d$O2_pct, d$y_bin, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(d, split_key), function(z) {
    data.frame(
      O2_pct = z$O2_pct[[1L]],
      y_bin = z$y_bin[[1L]],
      log10_p_center = centers[z$y_bin[[1L]]],
      n = nrow(z),
      median_dominant_mean_ploidy = stats::median(z$dominant_mean_ploidy),
      q25_dominant_mean_ploidy = finite_quantile(z$dominant_mean_ploidy, 0.25),
      q75_dominant_mean_ploidy = finite_quantile(z$dominant_mean_ploidy, 0.75),
      median_spectral_gap = stats::median(z$spectral_gap),
      fraction_spectral_gap_below_0p005 = mean(z$spectral_gap < low_gap_threshold),
      stringsAsFactors = FALSE
    )
  })
  occupied <- do.call(rbind, rows)
  full <- merge(
    expand.grid(O2_pct = sort(unique(d$O2_pct)), y_bin = seq_along(centers)),
    occupied,
    by = c("O2_pct", "y_bin"),
    all.x = TRUE,
    sort = TRUE
  )
  full$log10_p_center <- centers[full$y_bin]
  full$supported <- is.finite(full$n) & full$n >= min_bin_n
  full$plot_median_ploidy <- ifelse(full$supported, full$median_dominant_mean_ploidy, NA_real_)
  full$plot_fraction_low_gap <- ifelse(full$supported, full$fraction_spectral_gap_below_0p005, NA_real_)
  list(points = d, occupied = occupied, full = full, min_bin_n = min_bin_n, low_gap_threshold = low_gap_threshold)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

write_dataset_outputs <- function(dataset, points, binned, results_dir) {
  point_path <- file.path(results_dir, paste0(dataset, "_fixed_o2_cin_points.tsv"))
  bin_path <- file.path(results_dir, paste0(dataset, "_fixed_o2_cin_binned.tsv"))
  qc_path <- file.path(results_dir, paste0(dataset, "_fixed_o2_cin_qc_summary.csv"))
  write_tsv(points, point_path)
  write_tsv(binned$occupied, bin_path)
  qc <- data.frame(
    metric = c(
      "n_seed", "n_o2", "n_points_expected", "n_points_finite", "n_failed_points",
      "population_average_p_misseg_min", "population_average_p_misseg_max",
      "dominant_mean_ploidy_min", "dominant_mean_ploidy_max",
      "spectral_gap_min", "spectral_gap_median", "fraction_spectral_gap_below_0p005",
      "n_occupied_bins", "n_supported_bins", "minimum_supported_bin_n"
    ),
    value = c(
      length(unique(points$seed)), length(unique(points$O2_pct)),
      length(unique(points$seed)) * length(unique(points$O2_pct)), nrow(binned$points), nrow(points) - nrow(binned$points),
      min(binned$points$population_average_p_misseg), max(binned$points$population_average_p_misseg),
      min(binned$points$dominant_mean_ploidy), max(binned$points$dominant_mean_ploidy),
      min(binned$points$spectral_gap), stats::median(binned$points$spectral_gap),
      mean(binned$points$spectral_gap < binned$low_gap_threshold),
      nrow(binned$occupied), sum(binned$full$supported), binned$min_bin_n
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(qc, qc_path, row.names = FALSE, quote = TRUE)
  c(point_path, bin_path, qc_path)
}

scientific_labels <- function(log10_values) {
  formatC(10^log10_values, format = "e", digits = 0)
}

curve_class_order <- c(
  "complex_nonmonotone", "inverted_u_shaped", "monotone_increasing", "u_shaped",
  "approximately_flat", "single_transition_increase_then_plateau",
  "single_transition_decrease_then_plateau", "monotone_decreasing"
)

curve_class_colors <- c(
  complex_nonmonotone = "#6A5ACD",
  inverted_u_shaped = "#D55E00",
  monotone_increasing = "#009E73",
  u_shaped = "#E69F00",
  approximately_flat = "#666666",
  single_transition_increase_then_plateau = "#CC79A7",
  single_transition_decrease_then_plateau = "#56B4E9",
  monotone_decreasing = "#0072B2"
)

curve_class_labels <- c(
  complex_nonmonotone = "Complex nonmonotone",
  inverted_u_shaped = "Inverted U",
  monotone_increasing = "Monotone increasing",
  u_shaped = "U-shaped",
  approximately_flat = "Approximately flat",
  single_transition_increase_then_plateau = "Increase then plateau",
  single_transition_decrease_then_plateau = "Decrease then plateau",
  monotone_decreasing = "Monotone decreasing"
)

selection_diagnostic_plot <- function(path) {
  if (!file.exists(path)) stop("Missing Figure 6D landscape table: ", path, call. = FALSE)
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("tSNE1", "tSNE2", "seed", "objective", "curve_class")
  if (!all(required %in% names(d))) stop("Figure 6D table is missing: ", paste(setdiff(required, names(d)), collapse = ", "), call. = FALSE)
  if (nrow(d) != 500L || anyDuplicated(d$seed) || any(!is.finite(d$tSNE1)) ||
      any(!is.finite(d$tSNE2)) || any(!is.finite(d$objective)) || anyNA(d$curve_class)) {
    stop("Figure 6D requires 500 unique, finite, classified in vivo best-fit points.", call. = FALSE)
  }
  unknown <- setdiff(unique(d$curve_class), curve_class_order)
  if (length(unknown)) stop("Unknown Figure 6D curve classes: ", paste(unknown, collapse = ", "), call. = FALSE)
  d$curve_class <- factor(d$curve_class, levels = curve_class_order)
  class_n <- table(d$curve_class)
  class_axis_labels <- setNames(
    paste0(curve_class_labels[names(class_n)], " (n=", as.integer(class_n), ")"),
    names(class_n)
  )
  best_cut <- as.numeric(stats::quantile(d$objective, 0.10, names = FALSE, type = 7))
  d$lowest_objective_decile <- d$objective <= best_cut

  p_landscape <- ggplot(d, aes(x = tSNE1, y = tSNE2, color = curve_class)) +
    geom_point(size = 1.55, alpha = 0.82) +
    geom_point(
      data = d[d$lowest_objective_decile, , drop = FALSE],
      shape = 21, fill = NA, color = "#111111", stroke = 0.7, size = 2.6,
      inherit.aes = TRUE, show.legend = FALSE
    ) +
    scale_color_manual(values = curve_class_colors, labels = curve_class_labels, drop = FALSE) +
    labs(
      x = "Parameter t-SNE 1", y = "Parameter t-SNE 2",
      title = "Response classes in parameter space",
      subtitle = "Black rings: lowest objective decile",
      color = "O2-ploidy response class", tag = "D"
    ) +
    theme_classic(base_size = 9.5, base_family = "Helvetica") +
    theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 8.5, color = "#555555"),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7.2),
      legend.key.height = unit(0.3, "cm")
    ) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 2.5, alpha = 1)))

  p_objective <- ggplot(d, aes(x = objective, y = curve_class, color = curve_class)) +
    geom_vline(xintercept = best_cut, color = "#333333", linetype = "dotted", linewidth = 0.45) +
    geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.45, alpha = 0) +
    geom_jitter(height = 0.16, width = 0, size = 0.8, alpha = 0.28) +
    scale_color_manual(values = curve_class_colors, guide = "none", drop = FALSE) +
    scale_y_discrete(labels = class_axis_labels, drop = FALSE) +
    labs(
      x = "In vivo objective (lower is better)", y = NULL,
      title = "Fit quality overlaps across classes",
      subtitle = "Dotted line: lowest-decile cutoff"
    ) +
    theme_classic(base_size = 9.5, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 8.5, color = "#555555"),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7.5)
    )

  p_landscape + p_objective + plot_layout(widths = c(1.12, 1))
}

surface_plot <- function(binned, title, tag = NULL) {
  d <- binned$full
  tiles <- d[d$supported, , drop = FALSE]
  fill_range <- range(d$plot_median_ploidy, c(1, 4), na.rm = TRUE)
  y_breaks <- pretty(range(d$log10_p_center), n = 5)
  tile_width <- min(diff(sort(unique(d$O2_pct))))
  tile_height <- min(diff(sort(unique(d$log10_p_center))))
  ggplot() +
    geom_tile(
      data = tiles,
      aes(x = O2_pct, y = log10_p_center, fill = plot_median_ploidy),
      width = tile_width, height = tile_height,
      na.rm = TRUE
    ) +
    geom_contour(
      data = d,
      aes(x = O2_pct, y = log10_p_center, z = plot_median_ploidy, linetype = after_stat(factor(level))),
      breaks = c(2, 4), color = "#171717", linewidth = 0.45, na.rm = TRUE
    ) +
    geom_contour(
      data = d,
      aes(x = O2_pct, y = log10_p_center, z = plot_fraction_low_gap), breaks = 0.5,
      color = "white", linetype = "dotted", linewidth = 0.45, na.rm = TRUE,
      show.legend = FALSE
    ) +
    scale_fill_gradientn(
      colours = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
      limits = fill_range, oob = scales::squish, na.value = NA, name = "Median dominant\nmean ploidy"
    ) +
    scale_linetype_manual(
      values = c("2" = "solid", "4" = "longdash"),
      labels = c("2" = "2N", "4" = "4N"),
      name = "Ploidy contour",
      drop = FALSE
    ) +
    scale_x_continuous(breaks = 0:5, expand = c(0, 0)) +
    scale_y_continuous(breaks = y_breaks, labels = scientific_labels, expand = c(0, 0)) +
    labs(
      x = "Fixed oxygen (%)",
      y = "Population-average per-chromosome\nmissegregation probability",
      title = title,
      subtitle = "Projection of fitted solutions into the fixed-O2/CIN plane",
      caption = "Tiles require n >= 3 fitted solutions. White dotted contour marks bins where >50% have spectral gap <0.005.",
      tag = tag
    ) +
    coord_cartesian(expand = FALSE) +
    theme_classic(base_size = 10, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "#555555"),
      plot.caption = element_text(size = 7.5, color = "#555555", hjust = 0),
      plot.tag = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 8),
      axis.title = element_text(size = 9.5),
      axis.text = element_text(size = 8.5)
    )
}

save_plot_pair <- function(plot, base_path, width, height) {
  ggsave(paste0(base_path, ".png"), plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggsave(paste0(base_path, ".pdf"), plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white")
  c(paste0(base_path, ".png"), paste0(base_path, ".pdf"))
}

main <- function(argv = parse_args()) {
  results_dir <- normalizePath(
    argv$results_dir %||% file.path(script_dir, "..", "results", "figure6_fixed_o2_cin_map"),
    mustWork = TRUE
  )
  invivo_path <- file.path(results_dir, "invivo_best_fixo2_eigen_attractor_wide.csv")
  invitro_path <- file.path(results_dir, "invitro_best_fixo2_eigen_attractor_wide.csv")
  joint_path <- file.path(results_dir, "joint_best_fixo2_eigen_attractor_long.tsv")
  selection_path <- file.path(results_dir, "invivo_parameter_landscape_tsne_curve_class.csv")
  panel_a_path <- file.path(script_dir, "iteration1", "fig6a_pooled_fixed_o2_curve_class_examples.png")
  if (!file.exists(panel_a_path)) stop("Missing Figure 6A panel: ", panel_a_path, call. = FALSE)

  invivo_points <- wide_to_points(read_wide(invivo_path), "invivo")
  invitro_points <- wide_to_points(read_wide(invitro_path), "invitro")
  if (nrow(invivo_points) != 500L * 201L || nrow(invitro_points) != 500L * 201L) {
    stop("Expected 100,500 points for each separate-fit dataset.", call. = FALSE)
  }
  invivo_binned <- bin_surface(invivo_points)
  invitro_binned <- bin_surface(invitro_points)
  result_outputs <- c(
    write_dataset_outputs("invivo", invivo_points, invivo_binned, results_dir),
    write_dataset_outputs("invitro", invitro_points, invitro_binned, results_dir)
  )

  p_main <- surface_plot(invivo_binned, "In vivo fitted solutions", tag = "B")
  p_invitro <- surface_plot(invitro_binned, "In vitro sensitivity analysis")
  p_selection <- selection_diagnostic_plot(selection_path)
  panel_b_base <- file.path(script_dir, "iteration1", "fig6b_fixed_o2_cin_ploidy_map")
  supporting_dir <- file.path(script_dir, "supporting")
  dir.create(supporting_dir, recursive = TRUE, showWarnings = FALSE)
  figure_outputs <- c(
    save_plot_pair(p_main, panel_b_base, 7.4, 5.2),
    save_plot_pair(p_invitro, file.path(supporting_dir, "fig6b_invitro_sensitivity"), 7.4, 5.2),
    save_plot_pair(p_selection, file.path(script_dir, "iteration1", "fig6d_selection_diagnostic"), 8.2, 4.4)
  )

  if (file.exists(joint_path)) {
    joint <- utils::read.delim(joint_path, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("pair_id", "O2_pct", "population_average_p_misseg")
    if (!all(required %in% names(joint))) stop("Joint overlay table is missing: ", paste(setdiff(required, names(joint)), collapse = ", "), call. = FALSE)
    joint <- joint[is.finite(joint$O2_pct) & is.finite(joint$population_average_p_misseg) & joint$population_average_p_misseg > 0, , drop = FALSE]
    joint$log10_population_average_p_misseg <- log10(joint$population_average_p_misseg)
    joint$pair_label <- sub("^fit_joint_tsne_", "", joint$pair_id)
    joint$pair_label <- sub("_vt_seed10$", "", joint$pair_label)
    p_joint <- p_main +
      geom_path(
        data = joint,
        aes(x = O2_pct, y = log10_population_average_p_misseg, color = pair_label, group = pair_id),
        inherit.aes = FALSE, linewidth = 0.65
      ) +
      scale_color_discrete(name = "Joint best pair") +
      labs(title = "In vivo surface with six joint best-candidate paths", tag = NULL)
    figure_outputs <- c(
      figure_outputs,
      save_plot_pair(p_joint, file.path(supporting_dir, "fig6b_joint_best_overlay"), 8.2, 5.4)
    )
  }

  panel_a_grob <- grobTree(
    rasterGrob(png::readPNG(panel_a_path), interpolate = TRUE),
    textGrob("A", x = unit(0.008, "npc"), y = unit(0.992, "npc"), just = c("left", "top"), gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "Helvetica"))
  )
  panel_a <- patchwork::wrap_elements(full = panel_a_grob)
  assembled <- panel_a / p_main + plot_layout(heights = c(1.0, 1.05))
  figure_outputs <- c(
    figure_outputs,
    save_plot_pair(assembled, file.path(script_dir, "assembled_fig6"), 7.8, 10.0),
    save_plot_pair(p_main, file.path(script_dir, "FixedO2_vs_ploidy"), 7.4, 5.2)
  )

  manifest <- data.frame(
    output_file = normalizePath(c(result_outputs, figure_outputs), mustWork = FALSE),
    role = c(
      rep(c("analysis_points", "analysis_bins", "quality_control"), 2L),
      rep("figure", length(figure_outputs))
    ),
    generating_script = normalizePath(file.path(script_dir, "draw_figure6_fixed_o2_cin_map.R"), mustWork = TRUE),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(results_dir, "figure6_output_manifest.csv")
  utils::write.csv(manifest, manifest_path, row.names = FALSE, quote = TRUE)
  message("Wrote Figure 6 outputs and manifest: ", manifest_path)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(as.character(x[[1L]]))) y else x

if (sys.nframe() == 0L) main()
