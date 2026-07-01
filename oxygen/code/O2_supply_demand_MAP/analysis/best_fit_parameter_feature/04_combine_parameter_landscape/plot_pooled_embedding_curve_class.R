#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

UTIL_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "util"), mustWork = FALSE)
source(file.path(UTIL_DIR, "path_utils.R"))
source(file.path(UTIL_DIR, "cli_utils.R"))

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript plot_pooled_embedding_curve_class.R [options]",
      "",
      "Purpose:",
      "  Overlay dense-grid regression curve classes onto pooled in vivo/in vitro PCA, UMAP, and t-SNE plots.",
      "",
      "Inputs:",
      "  --pooled_root=DIR       Directory containing PCAs/UMAPs/TSNEs Tables and Figures.",
      "  --parameter_root=DIR    Parameter-landscape result root. Used when pooled_root is omitted.",
      "  --dense_grid_dir=DIR    Dense-grid monotonicity result root searched for the latest regression by-seed table.",
      "  --class_table=FILE      Optional explicit fixed_o2_ploidy_monotonicity_regression_by_seed.tsv.",
      "  --class_col=NAME        Classification column to map onto in vivo best points. Default: curve_class.",
      "",
      "Outputs:",
      "  --out_dir=DIR           Output root for curve-class figures and tables.",
      "  --reductions=LIST       PCAs,UMAPs,TSNEs or aliases pca,umap,tsne. Default: all three.",
      "  --variants=LIST         Full,Sampled. Default: both.",
      "  --dry_run=TRUE|FALSE    Print planned inputs/outputs without writing files.",
      sep = "\n"
    ),
    "\n"
  )
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}

read_csv_plain <- function(path) {
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_tsv_plain <- function(path) {
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = path, row.names = FALSE, na = "")
}

normalize_reductions <- function(x) {
  vals <- tolower(gsub("-", "_", bpf_split_csv(x, c("PCAs", "UMAPs", "TSNEs"))))
  out <- character()
  for (val in vals) {
    if (val %in% c("pca", "pcas")) {
      out <- c(out, "PCAs")
    } else if (val %in% c("umap", "umaps")) {
      out <- c(out, "UMAPs")
    } else if (val %in% c("tsne", "tsnes", "t_sne", "t_snes")) {
      out <- c(out, "TSNEs")
    } else if (nzchar(val)) {
      stop("Unknown reduction: ", val, call. = FALSE)
    }
  }
  unique(out)
}

normalize_variants <- function(x) {
  vals <- tolower(bpf_split_csv(x, c("Full", "Sampled")))
  out <- character()
  for (val in vals) {
    if (val %in% c("full", "all_points")) {
      out <- c(out, "Full")
    } else if (val %in% c("sampled", "sample", "sampled500")) {
      out <- c(out, "Sampled")
    } else if (nzchar(val)) {
      stop("Unknown variant: ", val, call. = FALSE)
    }
  }
  unique(out)
}

curve_class_colors <- function(classes) {
  pal <- c(
    complex_nonmonotone = "#6A5ACD",
    inverted_u_shaped = "#D55E00",
    monotone_increasing = "#009E73",
    monotone_decreasing = "#0072B2",
    single_transition_increase_then_plateau = "#CC79A7",
    single_transition_decrease_then_plateau = "#56B4E9",
    u_shaped = "#E69F00",
    approximately_flat = "#666666",
    insufficient_data = "#999999",
    missing_curve_class = "#999999"
  )
  missing <- setdiff(classes, names(pal))
  if (length(missing)) {
    pal <- c(pal, stats::setNames(grDevices::rainbow(length(missing)), missing))
  }
  pal[classes]
}

find_latest_class_table <- function(dense_grid_dir, explicit_path = NULL) {
  if (!is.null(explicit_path) && length(explicit_path) && !is.na(explicit_path[[1L]]) && nzchar(explicit_path[[1L]])) {
    path <- normalizePath(explicit_path[[1L]], mustWork = TRUE)
    return(path)
  }
  candidates <- list.files(
    dense_grid_dir,
    pattern = "^fixed_o2_ploidy_monotonicity_regression_by_seed\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(candidates)) {
    stop("No regression by-seed class table found under: ", dense_grid_dir, call. = FALSE)
  }
  info <- file.info(candidates)
  candidates[[which.max(info$mtime)]]
}

read_curve_class_map <- function(path, class_col = "curve_class") {
  tab <- read_tsv_plain(path)
  if (!"seed_number" %in% names(tab)) {
    if (!"seed_id" %in% names(tab)) {
      stop("Class table must contain seed_number or seed_id: ", path, call. = FALSE)
    }
    tab$seed_number <- suppressWarnings(as.integer(sub("^seed", "", as.character(tab$seed_id))))
  }
  if (!class_col %in% names(tab)) {
    stop("Class table is missing --class_col=", class_col, ": ", path, call. = FALSE)
  }
  tab$seed_number <- suppressWarnings(as.integer(tab$seed_number))
  tab[[class_col]] <- trimws(as.character(tab[[class_col]]))
  tab <- tab[is.finite(tab$seed_number) & nzchar(tab[[class_col]]), , drop = FALSE]
  if (!nrow(tab)) stop("Class table has no usable seed/class rows: ", path, call. = FALSE)
  dup <- duplicated(tab$seed_number)
  if (any(dup)) {
    dup_seeds <- unique(tab$seed_number[dup])
    inconsistent <- vapply(dup_seeds, function(seed) {
      length(unique(tab[[class_col]][tab$seed_number == seed])) > 1L
    }, logical(1))
    if (any(inconsistent)) {
      stop(
        "Class table has duplicate seed_number rows with conflicting classes: ",
        paste(dup_seeds[inconsistent], collapse = ", "),
        call. = FALSE
      )
    }
    tab <- tab[!duplicated(tab$seed_number), , drop = FALSE]
  }
  stats::setNames(tab[[class_col]], as.character(tab$seed_number))
}

discover_coordinate_tables <- function(pooled_root, reductions, variants) {
  rows <- list()
  for (reduction in reductions) {
    for (variant in variants) {
      table_dir <- file.path(pooled_root, reduction, "Tables", variant)
      if (!dir.exists(table_dir)) next
      files <- list.files(
        table_dir,
        pattern = "_coordinates\\.csv$",
        full.names = TRUE
      )
      files <- files[order(files)]
      for (path in files) {
        stub <- sub("_coordinates\\.csv$", "", basename(path))
        rows[[length(rows) + 1L]] <- data.frame(
          reduction = reduction,
          variant = variant,
          coordinate_table = normalizePath(path, mustWork = FALSE),
          original_png = original_png_path(path),
          stub = stub,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      reduction = character(),
      variant = character(),
      coordinate_table = character(),
      original_png = character(),
      stub = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

original_png_path <- function(coordinate_table) {
  path <- gsub("/Tables/", "/Figures/", coordinate_table, fixed = TRUE)
  sub("_coordinates\\.csv$", ".png", path)
}

coord_spec <- function(reduction, data_names) {
  candidates <- switch(
    reduction,
    PCAs = list(c("PCA1", "PCA2")),
    UMAPs = list(c("UMAP1", "UMAP2")),
    TSNEs = list(c("tSNE1", "tSNE2"), c("TSNE1", "TSNE2")),
    stop("Unknown reduction: ", reduction, call. = FALSE)
  )
  for (candidate in candidates) {
    if (all(candidate %in% data_names)) {
      labels <- switch(
        reduction,
        PCAs = c("PCA 1", "PCA 2"),
        UMAPs = c("UMAP 1", "UMAP 2"),
        TSNEs = c("t-SNE 1", "t-SNE 2")
      )
      return(list(columns = candidate, labels = labels))
    }
  }
  stop(
    "Could not find coordinate columns for ",
    reduction,
    ". Available columns: ",
    paste(data_names, collapse = ", "),
    call. = FALSE
  )
}

square_limits <- function(plot_data, coord_names) {
  x <- suppressWarnings(as.numeric(plot_data[[coord_names[[1L]]]]))
  y <- suppressWarnings(as.numeric(plot_data[[coord_names[[2L]]]]))
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) stop("No finite embedding coordinates.", call. = FALSE)
  x_range <- range(x)
  y_range <- range(y)
  x_center <- mean(x_range)
  y_center <- mean(y_range)
  span <- max(diff(x_range), diff(y_range))
  if (!is.finite(span) || span <= 0) span <- 1
  pad <- span * 0.03
  span <- span + 2 * pad
  list(
    xlim = c(x_center - span / 2, x_center + span / 2),
    ylim = c(y_center - span / 2, y_center + span / 2)
  )
}

add_curve_class <- function(plot_data, class_map) {
  required <- c("dataset", "point_type", "seed")
  missing <- setdiff(required, names(plot_data))
  if (length(missing)) {
    stop("Coordinate table is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  plot_data$seed_number <- suppressWarnings(as.integer(plot_data$seed))
  plot_data$curve_class <- NA_character_
  idx <- plot_data$dataset == "invivo" & plot_data$point_type == "best" & is.finite(plot_data$seed_number)
  plot_data$curve_class[idx] <- unname(class_map[as.character(plot_data$seed_number[idx])])
  plot_data
}

class_labels <- function(classes) {
  labels <- gsub("_", " ", classes, fixed = TRUE)
  labels[classes == "missing_curve_class"] <- "missing curve class"
  stats::setNames(labels, classes)
}

prepare_plot_data <- function(plot_data, coord_names) {
  plot_data$.embedding_x <- suppressWarnings(as.numeric(plot_data[[coord_names[[1L]]]]))
  plot_data$.embedding_y <- suppressWarnings(as.numeric(plot_data[[coord_names[[2L]]]]))
  plot_data$dataset <- as.character(plot_data$dataset)
  plot_data$point_type <- as.character(plot_data$point_type)
  plot_data
}

initial_points <- function(plot_data) {
  initial <- plot_data[plot_data$point_type == "initial", , drop = FALSE]
  initial$initial_group <- ifelse(initial$dataset == "invivo", "invivo_initial", "invitro_initial")
  initial$initial_group <- factor(initial$initial_group, levels = c("invivo_initial", "invitro_initial"))
  initial
}

initial_legend_points <- function(lims) {
  data.frame(
    .embedding_x = rep(mean(lims$xlim), 2L),
    .embedding_y = rep(mean(lims$ylim), 2L),
    initial_group = factor(c("invivo_initial", "invitro_initial"), levels = c("invivo_initial", "invitro_initial"))
  )
}

common_embedding_theme <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      legend.box = "vertical",
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

build_objective_plot <- function(plot_data,
                                 coord_names,
                                 axis_labels,
                                 lims,
                                 initial_size = 0.22,
                                 best_size = 1.25,
                                 initial_invivo_color = "#7FA2FF",
                                 initial_invitro_color = "#FF4FBF",
                                 initial_invivo_alpha = 0.48,
                                 initial_invitro_alpha = 0.45) {
  plot_data <- prepare_plot_data(plot_data, coord_names)
  initial <- initial_points(plot_data)
  initial_invivo <- initial[initial$dataset == "invivo", , drop = FALSE]
  initial_invitro <- initial[initial$dataset == "invitro", , drop = FALSE]
  initial_legend <- initial_legend_points(lims)
  best_invivo <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invivo", , drop = FALSE]
  best_invitro <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invitro", , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = initial_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      color = initial_invivo_color,
      shape = 16,
      size = initial_size,
      alpha = initial_invivo_alpha,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      color = initial_invitro_color,
      shape = 17,
      size = initial_size,
      alpha = initial_invitro_alpha,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_legend,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = initial_group, shape = initial_group),
      color = "transparent",
      alpha = 0,
      size = 0.01,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Initial",
      values = c(invivo_initial = initial_invivo_color, invitro_initial = initial_invitro_color),
      labels = c(invivo_initial = "in vivo", invitro_initial = "in vitro"),
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          alpha = 1,
          size = 3,
          shape = c(21, 24),
          color = "transparent"
        )
      )
    ) +
    ggplot2::scale_shape_manual(
      name = "Initial",
      values = c(invivo_initial = 21, invitro_initial = 24),
      labels = c(invivo_initial = "in vivo", invitro_initial = "in vitro"),
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::geom_point(
      data = best_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = objective),
      shape = 16,
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_gradient(
      name = "in vivo\nobjective",
      low = "#2C7BB6",
      high = "#FDE725",
      guide = ggplot2::guide_colorbar(order = 2, barheight = ggplot2::unit(26, "mm"), display = "rectangles")
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      data = best_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = objective),
      shape = 17,
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_gradient(
      name = "in vitro\nobjective",
      low = "#1A9850",
      high = "#D73027",
      guide = ggplot2::guide_colorbar(order = 3, barheight = ggplot2::unit(26, "mm"), display = "rectangles")
    ) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(x = axis_labels[[1L]], y = axis_labels[[2L]]) +
    common_embedding_theme()
}

build_curve_class_plot <- function(plot_data,
                                   reduction,
                                   coord_names,
                                   axis_labels,
                                   lims,
                                   initial_size = 0.22,
                                   best_size = 1.25,
                                   initial_invivo_color = "#7FA2FF",
                                   initial_invitro_color = "#FF4FBF",
                                   initial_invivo_alpha = 0.48,
                                   initial_invitro_alpha = 0.45,
                                   invitro_best_color = "#D73027") {
  plot_data <- prepare_plot_data(plot_data, coord_names)
  initial <- initial_points(plot_data)
  initial_invivo <- initial[initial$dataset == "invivo", , drop = FALSE]
  initial_invitro <- initial[initial$dataset == "invitro", , drop = FALSE]
  initial_legend <- initial_legend_points(lims)
  best_invivo <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invivo", , drop = FALSE]
  best_invitro <- plot_data[plot_data$point_type == "best" & plot_data$dataset == "invitro", , drop = FALSE]

  best_invivo$curve_class_for_plot <- ifelse(
    is.na(best_invivo$curve_class) | !nzchar(best_invivo$curve_class),
    "missing_curve_class",
    best_invivo$curve_class
  )
  class_order <- sort(setdiff(unique(best_invivo$curve_class_for_plot), "missing_curve_class"))
  if ("missing_curve_class" %in% best_invivo$curve_class_for_plot) {
    class_order <- c(class_order, "missing_curve_class")
  }
  best_invivo$curve_class_for_plot <- factor(best_invivo$curve_class_for_plot, levels = class_order)
  best_fit_legend <- data.frame(
    .embedding_x = rep(mean(lims$xlim), 2L),
    .embedding_y = rep(mean(lims$ylim), 2L),
    best_fit_group = factor(c("invivo_best", "invitro_best"), levels = c("invivo_best", "invitro_best"))
  )
  class_values <- curve_class_colors(class_order)
  labels <- class_labels(class_order)

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = initial_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      color = initial_invivo_color,
      shape = 16,
      size = initial_size,
      alpha = initial_invivo_alpha,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      color = initial_invitro_color,
      shape = 17,
      size = initial_size,
      alpha = initial_invitro_alpha,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_legend,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = initial_group, shape = initial_group),
      color = "transparent",
      alpha = 0,
      size = 0.01,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Initial",
      values = c(invivo_initial = initial_invivo_color, invitro_initial = initial_invitro_color),
      labels = c(invivo_initial = "in vivo", invitro_initial = "in vitro"),
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          alpha = 1,
          size = 3,
          shape = c(21, 24),
          color = "transparent"
        )
      )
    ) +
    ggplot2::scale_shape_manual(
      name = "Initial",
      values = c(invivo_initial = 21, invitro_initial = 24),
      labels = c(invivo_initial = "in vivo", invitro_initial = "in vitro"),
      drop = FALSE,
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
      data = best_fit_legend,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, fill = best_fit_group),
      shape = 21,
      color = "transparent",
      alpha = 0,
      size = 0.01,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Best fits",
      values = c(invivo_best = "#555555", invitro_best = invitro_best_color),
      labels = c(invivo_best = "in vivo", invitro_best = "in vitro"),
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          shape = c(16, 17),
          color = c("#555555", invitro_best_color),
          fill = c("#555555", invitro_best_color),
          alpha = 1,
          size = 3
        )
      )
    ) +
    ggplot2::geom_point(
      data = best_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 17,
      color = invitro_best_color,
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = best_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = curve_class_for_plot),
      shape = 16,
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = expression(O[2] * "-Ploidy Curve Class"),
      values = class_values,
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 3,
        override.aes = list(shape = 16, alpha = 1, size = 3)
      )
    ) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(
      x = axis_labels[[1L]],
      y = axis_labels[[2L]]
    ) +
    common_embedding_theme()
}

save_plot_pair <- function(plot, out_stem, width = 7.4, height = 6.5, dpi = 300) {
  png_path <- paste0(out_stem, ".png")
  pdf_path <- paste0(out_stem, ".pdf")
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
  c(png = png_path, pdf = pdf_path)
}

process_embedding <- function(row, class_map, out_dir) {
  plot_data <- read_csv_plain(row$coordinate_table)
  spec <- coord_spec(row$reduction, names(plot_data))
  plot_data <- add_curve_class(plot_data, class_map)
  lims <- square_limits(plot_data, spec$columns)
  best_size <- 1.25
  initial_size <- if (identical(row$variant, "Sampled")) 0.8 * best_size else 0.22
  n_invivo_best_missing <- sum(
    plot_data$dataset == "invivo" &
      plot_data$point_type == "best" &
      (is.na(plot_data$curve_class) | !nzchar(plot_data$curve_class))
  )
  original_plot <- build_objective_plot(
    plot_data = plot_data,
    coord_names = spec$columns,
    axis_labels = spec$labels,
    lims = lims,
    initial_size = initial_size,
    best_size = best_size
  )
  curve_plot <- build_curve_class_plot(
    plot_data = plot_data,
    reduction = row$reduction,
    coord_names = spec$columns,
    axis_labels = spec$labels,
    lims = lims,
    initial_size = initial_size,
    best_size = best_size
  )
  combined_plot <- original_plot + curve_plot + patchwork::plot_layout(ncol = 2, widths = c(1, 1))

  figure_dir <- file.path(out_dir, "figures", row$reduction, row$variant)
  table_dir <- file.path(out_dir, "tables", row$reduction, row$variant)
  curve_stem <- file.path(figure_dir, paste0(row$stub, "_curve_class"))
  combined_stem <- file.path(figure_dir, paste0(row$stub, "_original_vs_curve_class"))
  plot_paths <- save_plot_pair(curve_plot, curve_stem)
  combined_paths <- save_plot_pair(combined_plot, combined_stem, width = 14.8, height = 6.5)

  best_points <- plot_data[plot_data$point_type == "best", , drop = FALSE]
  write_csv(best_points, file.path(table_dir, paste0(row$stub, "_best_points_curve_class.csv")))

  counts <- as.data.frame(table(
    dataset = best_points$dataset,
    curve_class = ifelse(is.na(best_points$curve_class) | !nzchar(best_points$curve_class), "not_applicable_or_missing", best_points$curve_class),
    useNA = "no"
  ), stringsAsFactors = FALSE)
  counts <- counts[counts$Freq > 0, , drop = FALSE]
  write_tsv(counts, file.path(table_dir, paste0(row$stub, "_curve_class_counts.tsv")))

  data.frame(
    reduction = row$reduction,
    variant = row$variant,
    stub = row$stub,
    coordinate_table = row$coordinate_table,
    original_png = row$original_png,
    curve_png = unname(plot_paths[["png"]]),
    curve_pdf = unname(plot_paths[["pdf"]]),
    combined_png = unname(combined_paths[["png"]]),
    combined_pdf = unname(combined_paths[["pdf"]]),
    n_rows = nrow(plot_data),
    n_initial = sum(plot_data$point_type == "initial"),
    n_best_invivo = sum(plot_data$dataset == "invivo" & plot_data$point_type == "best"),
    n_best_invitro = sum(plot_data$dataset == "invitro" & plot_data$point_type == "best"),
    n_invivo_best_with_class = sum(plot_data$dataset == "invivo" & plot_data$point_type == "best" & !is.na(plot_data$curve_class) & nzchar(plot_data$curve_class)),
    n_invivo_best_missing_class = n_invivo_best_missing,
    stringsAsFactors = FALSE
  )
}

main <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  argv <- bpf_parse_args(raw_args)
  if (bpf_as_bool(argv$help %||% argv$h, FALSE)) {
    usage()
    return(invisible(NULL))
  }

  require_package("ggplot2")
  require_package("ggnewscale")
  require_package("patchwork")

  repo_root <- bpf_repo_root(SCRIPT_DIR)
  parameter_root <- bpf_resolve_repo_path(
    argv$parameter_root %||% bpf_parameter_landscape_result_dir(repo_root),
    repo_root = repo_root,
    mustWork = TRUE
  )
  pooled_root <- bpf_resolve_repo_path(
    argv$pooled_root %||% file.path(parameter_root, "pooled_invivo_invitro"),
    repo_root = repo_root,
    mustWork = TRUE
  )
  dense_grid_dir <- bpf_resolve_repo_path(
    argv$dense_grid_dir %||% bpf_dense_grid_result_root(repo_root),
    repo_root = repo_root,
    mustWork = TRUE
  )
  out_dir <- bpf_resolve_repo_path(
    argv$out_dir %||% file.path(bpf_combine_result_dir(repo_root), "pooled_embedding_curve_class"),
    repo_root = repo_root,
    mustWork = FALSE
  )
  class_col <- argv$class_col %||% "curve_class"
  dry_run <- bpf_as_bool(argv$dry_run, FALSE)
  reductions <- normalize_reductions(argv$reductions)
  variants <- normalize_variants(argv$variants)

  class_table <- find_latest_class_table(dense_grid_dir, argv$class_table)
  embeddings <- discover_coordinate_tables(pooled_root, reductions, variants)
  if (!nrow(embeddings)) {
    stop("No pooled coordinate tables found under: ", pooled_root, call. = FALSE)
  }

  message("Pooled root: ", pooled_root)
  message("Class table: ", class_table)
  message("Output root: ", out_dir)
  message("Embedding coordinate tables: ", nrow(embeddings))

  if (isTRUE(dry_run)) {
    print(embeddings[, c("reduction", "variant", "coordinate_table", "original_png")], row.names = FALSE)
    return(invisible(embeddings))
  }

  class_map <- read_curve_class_map(class_table, class_col = class_col)
  manifest_rows <- vector("list", nrow(embeddings))
  for (i in seq_len(nrow(embeddings))) {
    row <- embeddings[i, , drop = FALSE]
    message("[", i, "/", nrow(embeddings), "] ", row$reduction, " ", row$variant, " ", row$stub)
    manifest_rows[[i]] <- process_embedding(row, class_map, out_dir)
  }
  manifest <- do.call(rbind, manifest_rows)
  manifest$class_table <- class_table
  manifest$class_col <- class_col
  manifest_path <- file.path(out_dir, "tables", "pooled_embedding_curve_class_manifest.tsv")
  write_tsv(manifest, manifest_path)
  message("Wrote manifest: ", manifest_path)
  invisible(manifest)
}

if (identical(environment(), globalenv())) {
  main()
}
