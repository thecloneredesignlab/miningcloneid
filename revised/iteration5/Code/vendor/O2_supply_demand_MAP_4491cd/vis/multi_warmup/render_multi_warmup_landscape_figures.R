#!/usr/bin/env Rscript

# Pure visualization of materialized landscape embedding/cluster tables.

.o2mw_landscape_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_multi_warmup_landscape_figures.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_landscape_workflow_root <- normalizePath(file.path(.o2mw_landscape_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_landscape_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

normalize_reduction <- function(reduction) {
  value <- tolower(gsub("[^a-z]", "", as.character(reduction[[1L]])))
  if (value == "pca") return("pca")
  if (value == "umap") return("umap")
  if (value %in% c("tsne", "tstochasticneighborembedding")) return("tsne")
  stop("Unknown reduction: ", reduction, call. = FALSE)
}
reduction_axis_labels <- function(reduction) {
  switch(normalize_reduction(reduction), pca = c("PC1", "PC2"), umap = c("UMAP1", "UMAP2"), tsne = c("t-SNE1", "t-SNE2"))
}
reduction_coordinate_names <- o2sd_reduction_coordinate_names

require_plotting <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Required R package is not installed for cluster figures: ggplot2", call. = FALSE)
  }
}

require_pooled_objective_plotting <- function() {
  require_plotting()
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Required R package is not installed for pooled two-scale objective coloring: ggnewscale", call. = FALSE)
  }
}

plot_frame <- function(df, coord_names) {
  out <- df
  out$.embedding_x <- suppressWarnings(as.numeric(out[[coord_names[[1L]]]]))
  out$.embedding_y <- suppressWarnings(as.numeric(out[[coord_names[[2L]]]]))
  out
}

square_embedding_limits <- function(plot_data, pad_frac = 0.05, coord_names = c("UMAP1", "UMAP2")) {
  coord_names <- as.character(coord_names)
  if (!all(coord_names %in% names(plot_data))) {
    stop("Plot data is missing coordinate columns: ", paste(setdiff(coord_names, names(plot_data)), collapse = ", "), call. = FALSE)
  }
  x <- suppressWarnings(as.numeric(plot_data[[coord_names[[1L]]]]))
  y <- suppressWarnings(as.numeric(plot_data[[coord_names[[2L]]]]))
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) stop("Plot data must contain finite embedding coordinates.", call. = FALSE)
  x_range <- range(x)
  y_range <- range(y)
  span <- max(diff(x_range), diff(y_range))
  if (!is.finite(span) || span <= 0) span <- 1
  pad_frac <- suppressWarnings(as.numeric(pad_frac))
  if (!is.finite(pad_frac) || pad_frac < 0) pad_frac <- 0.05
  span <- span * (1 + 2 * pad_frac)
  x_center <- mean(x_range)
  y_center <- mean(y_range)
  list(
    xlim = c(x_center - span / 2, x_center + span / 2),
    ylim = c(y_center - span / 2, y_center + span / 2)
  )
}

add_plot_label_margin_rows <- function(plot_data, coord_names = c("UMAP1", "UMAP2"), pad_frac = 0.14) {
  coord_names <- as.character(coord_names)
  pad_frac <- suppressWarnings(as.numeric(pad_frac))
  if (!is.finite(pad_frac) || pad_frac < 0) pad_frac <- 0.14
  lims <- square_embedding_limits(plot_data, pad_frac = pad_frac, coord_names = coord_names)
  margin_rows <- plot_data[rep(1L, 4L), , drop = FALSE]
  for (nm in names(margin_rows)) margin_rows[[nm]] <- NA
  margin_rows[[coord_names[[1L]]]] <- c(lims$xlim[[1L]], lims$xlim[[2L]], lims$xlim[[1L]], lims$xlim[[2L]])
  margin_rows[[coord_names[[2L]]]] <- c(lims$ylim[[1L]], lims$ylim[[1L]], lims$ylim[[2L]], lims$ylim[[2L]])
  if ("point_type" %in% names(margin_rows)) margin_rows$point_type <- "plot_margin"
  if ("dataset" %in% names(margin_rows)) margin_rows$dataset <- "plot_margin"
  rbind(plot_data, margin_rows)
}

cluster_palette <- function(ids) {
  ids <- sort(unique(as.character(ids)))
  ids <- ids[nzchar(ids) & !is.na(ids)]
  if (!length(ids)) return(stats::setNames(character(), character()))
  values <- grDevices::hcl.colors(length(ids), palette = "Dark 3")
  stats::setNames(values, ids)
}

cluster_hull_data <- function(df, cluster_col, expand = 1.035) {
  df <- df[!is.na(df[[cluster_col]]) & nzchar(as.character(df[[cluster_col]])), , drop = FALSE]
  if (!nrow(df)) return(data.frame())
  all_span <- max(
    diff(range(df$.embedding_x, finite = TRUE)),
    diff(range(df$.embedding_y, finite = TRUE))
  )
  point_radius <- if (is.finite(all_span) && all_span > 0) all_span * 0.012 else 0.1
  pieces <- lapply(split(df, as.character(df[[cluster_col]])), function(d) {
    xy <- unique(d[, c(".embedding_x", ".embedding_y"), drop = FALSE])
    if (nrow(xy) >= 3L) {
      idx <- grDevices::chull(xy$.embedding_x, xy$.embedding_y)
      h <- xy[c(idx, idx[[1L]]), , drop = FALSE]
      center <- colMeans(h[, c(".embedding_x", ".embedding_y"), drop = FALSE])
      h$.embedding_x <- center[[1L]] + (h$.embedding_x - center[[1L]]) * expand
      h$.embedding_y <- center[[2L]] + (h$.embedding_y - center[[2L]]) * expand
    } else if (nrow(xy) >= 1L) {
      theta <- seq(0, 2 * pi, length.out = 33L)
      h <- data.frame(
        .embedding_x = xy$.embedding_x[[1L]] + point_radius * cos(theta),
        .embedding_y = xy$.embedding_y[[1L]] + point_radius * sin(theta),
        stringsAsFactors = FALSE
      )
    } else {
      return(NULL)
    }
    h$cluster_id <- as.character(d[[cluster_col]][[1L]])
    h
  })
  rbind_fill_plain(Filter(Negate(is.null), pieces))
}

external_cluster_label_data <- function(clustered_plot_data,
                                        bounds_data,
                                        coord_names = c("UMAP1", "UMAP2"),
                                        cluster_col = "cluster_id",
                                        offset_frac = 0.07) {
  coord_names <- as.character(coord_names)
  lims <- square_embedding_limits(bounds_data, coord_names = coord_names)
  all_span <- max(diff(lims$xlim), diff(lims$ylim))
  if (!is.finite(all_span) || all_span <= 0) all_span <- 1
  margin <- all_span * 0.035
  offset <- all_span * offset_frac
  bx <- suppressWarnings(as.numeric(bounds_data[[coord_names[[1L]]]]))
  by <- suppressWarnings(as.numeric(bounds_data[[coord_names[[2L]]]]))
  global_center <- c(stats::median(bx, na.rm = TRUE), stats::median(by, na.rm = TRUE))
  if (any(!is.finite(global_center))) global_center <- c(mean(lims$xlim), mean(lims$ylim))
  cluster_ids <- sort(unique(as.character(clustered_plot_data[[cluster_col]])))
  cluster_ids <- cluster_ids[nzchar(cluster_ids) & !is.na(cluster_ids)]
  labels <- lapply(cluster_ids, function(cluster_id) {
    d <- clustered_plot_data[as.character(clustered_plot_data[[cluster_col]]) == cluster_id, coord_names, drop = FALSE]
    names(d) <- c(".embedding_x", ".embedding_y")
    d$.embedding_x <- suppressWarnings(as.numeric(d$.embedding_x))
    d$.embedding_y <- suppressWarnings(as.numeric(d$.embedding_y))
    d <- d[is.finite(d$.embedding_x) & is.finite(d$.embedding_y), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    center <- c(stats::median(d$.embedding_x), stats::median(d$.embedding_y))
    direction <- center - global_center
    norm <- sqrt(sum(direction^2))
    if (!is.finite(norm) || norm <= 1e-9) {
      far_i <- which.max((d$.embedding_x - center[[1L]])^2 + (d$.embedding_y - center[[2L]])^2)
      direction <- c(d$.embedding_x[[far_i]], d$.embedding_y[[far_i]]) - center
      norm <- sqrt(sum(direction^2))
    }
    if (!is.finite(norm) || norm <= 1e-9) {
      direction <- c(1, 0)
      norm <- 1
    }
    direction <- direction / norm
    boundary_i <- which.max((d$.embedding_x - center[[1L]]) * direction[[1L]] + (d$.embedding_y - center[[2L]]) * direction[[2L]])
    target_pos <- c(d$.embedding_x[[boundary_i]], d$.embedding_y[[boundary_i]])
    label_pos <- target_pos + direction * offset
    label_pos[[1L]] <- min(max(label_pos[[1L]], lims$xlim[[1L]] + margin), lims$xlim[[2L]] - margin)
    label_pos[[2L]] <- min(max(label_pos[[2L]], lims$ylim[[1L]] + margin), lims$ylim[[2L]] - margin)
    data.frame(
      cluster_id = cluster_id,
      .embedding_x = label_pos[[1L]],
      .embedding_y = label_pos[[2L]],
      .target_x = target_pos[[1L]],
      .target_y = target_pos[[2L]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, Filter(Negate(is.null), labels))
  if (is.null(out)) out <- data.frame()
  row.names(out) <- NULL
  out
}

add_cluster_outline_layers_external_labels <- function(plot,
                                                       clustered_plot_data,
                                                       bounds_data,
                                                       coord_names = c("UMAP1", "UMAP2"),
                                                       cluster_col = "cluster_id",
                                                       linewidth = 0.72) {
  clustered_plot_data <- plot_frame(clustered_plot_data, coord_names = coord_names)
  hulls <- cluster_hull_data(clustered_plot_data, cluster_col = cluster_col)
  if (!nrow(hulls)) return(plot)
  pal <- cluster_palette(hulls$cluster_id)
  for (cluster_id in names(pal)) {
    h <- hulls[hulls$cluster_id == cluster_id, , drop = FALSE]
    plot <- plot +
      ggplot2::geom_path(
        data = h,
        ggplot2::aes(x = .embedding_x, y = .embedding_y),
        inherit.aes = FALSE,
        color = pal[[cluster_id]],
        linewidth = linewidth,
        linetype = "dashed",
        lineend = "round",
        show.legend = FALSE
      )
  }
  labels <- external_cluster_label_data(
    clustered_plot_data = clustered_plot_data,
    bounds_data = bounds_data,
    coord_names = coord_names,
    cluster_col = cluster_col
  )
  if (!nrow(labels)) return(plot)
  labels$.nudge_x <- labels$.embedding_x - labels$.target_x
  labels$.nudge_y <- labels$.embedding_y - labels$.target_y
  lims <- square_embedding_limits(bounds_data, coord_names = coord_names)
  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    plot <- plot + ggnewscale::new_scale_color()
  }
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(
      plot +
        ggrepel::geom_label_repel(
          data = labels,
          ggplot2::aes(x = .target_x, y = .target_y, label = cluster_id, color = cluster_id),
          inherit.aes = FALSE,
          fill = "white",
          label.size = 0.18,
          size = 2.4,
          fontface = "bold",
          box.padding = 0.4,
          point.padding = 0.15,
          min.segment.length = 0,
          segment.size = 0.25,
          nudge_x = labels$.nudge_x,
          nudge_y = labels$.nudge_y,
          seed = 123,
          max.overlaps = Inf,
          xlim = lims$xlim,
          ylim = lims$ylim,
          show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = pal, guide = "none")
    )
  }
  plot +
    ggplot2::geom_segment(
      data = labels,
      ggplot2::aes(x = .target_x, y = .target_y, xend = .embedding_x, yend = .embedding_y, color = cluster_id),
      inherit.aes = FALSE,
      linewidth = 0.25,
      alpha = 0.85,
      lineend = "round",
      show.legend = FALSE
    ) +
    ggplot2::geom_label(
      data = labels,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, label = cluster_id, color = cluster_id),
      inherit.aes = FALSE,
      fill = "white",
      linewidth = 0.18,
      size = 2.4,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = pal, guide = "none")
}

base_embedding_plot <- function(plot_data,
                                coord_names,
                                axis_labels,
                                initial_size = 0.22,
                                best_size = 1.25) {
  require_pooled_objective_plotting()
  lims <- square_embedding_limits(plot_data, coord_names = coord_names)
  dat <- plot_frame(plot_data, coord_names = coord_names)
  dat$dataset <- factor(as.character(dat$dataset), levels = c("invivo", "invitro", "plot_margin"))
  dat$point_type <- as.character(dat$point_type)
  initial <- dat[dat$point_type == "initial", , drop = FALSE]
  initial_invivo <- initial[initial$dataset == "invivo", , drop = FALSE]
  initial_invitro <- initial[initial$dataset == "invitro", , drop = FALSE]
  best_invivo <- dat[dat$point_type == "best" & dat$dataset == "invivo", , drop = FALSE]
  best_invitro <- dat[dat$point_type == "best" & dat$dataset == "invitro", , drop = FALSE]
  initial_invivo_color <- "#7FA2FF"
  initial_invitro_color <- "#FF4FBF"
  initial_invivo_alpha <- 0.48
  initial_invitro_alpha <- 0.45
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = initial_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 16,
      color = initial_invivo_color,
      alpha = initial_invivo_alpha,
      size = initial_size,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 17,
      color = initial_invitro_color,
      alpha = initial_invitro_alpha,
      size = initial_size,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = "Initial samples",
      values = c(invivo = 16, invitro = 17),
      labels = c(invivo = "in vivo", invitro = "in vitro")
    ) +
    ggplot2::geom_point(
      data = best_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = objective, shape = dataset),
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_gradient(
      name = "in vivo\nobjective",
      low = "#2C7BB6",
      high = "#FDE725",
      guide = ggplot2::guide_colorbar(order = 2, barheight = ggplot2::unit(26, "mm"))
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      data = best_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = objective, shape = dataset),
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_gradient(
      name = "in vitro\nobjective",
      low = "#1A9850",
      high = "#D73027",
      guide = ggplot2::guide_colorbar(order = 3, barheight = ggplot2::unit(26, "mm"))
    )
  p +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(x = axis_labels[[1L]], y = axis_labels[[2L]]) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          color = c(initial_invivo_color, initial_invitro_color),
          alpha = c(initial_invivo_alpha, initial_invitro_alpha),
          size = 3
        )
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      legend.box = "vertical",
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

subcluster_embedding_plot <- function(subclustered_best, coord_names, axis_labels, title = NULL) {
  require_plotting()
  dat <- plot_frame(subclustered_best, coord_names = coord_names)
  bounds_data <- add_plot_label_margin_rows(dat, coord_names = coord_names, pad_frac = 0.2)
  lims <- square_embedding_limits(bounds_data, coord_names = coord_names)
  dat$subcluster_id <- as.character(dat$subcluster_id)
  dat$dataset <- factor(as.character(dat$dataset), levels = c("invivo", "invitro"))
  pal <- cluster_palette(dat$subcluster_id)
  ggplot2::ggplot(dat, ggplot2::aes(x = .embedding_x, y = .embedding_y)) +
    ggplot2::geom_point(
      ggplot2::aes(color = subcluster_id, shape = dataset),
      size = 1.7,
      alpha = 0.95,
      stroke = 0.15
    ) +
    ggplot2::scale_color_manual(name = "subcluster", values = pal, guide = "none") +
    ggplot2::scale_shape_manual(
      name = NULL,
      values = c(invivo = 16, invitro = 17),
      labels = c(invivo = "in vivo", invitro = "in vitro")
    ) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(x = axis_labels[[1L]], y = axis_labels[[2L]], title = title) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

save_plot_pair <- function(plot, stem, width = 7.4, height = 6.5) {
  require_plotting()
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 220, bg = "white")
  invisible(stem)
}

write_cluster_figures <- function(full_marked,
                                  clustered_best,
                                  subclustered_best,
                                  figure_dir,
                                  output_stem,
                                  reduction,
                                  coord_names) {
  require_plotting()
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  axis_labels <- reduction_axis_labels(reduction)
  plot_bounds <- add_plot_label_margin_rows(full_marked, coord_names = coord_names, pad_frac = 0.14)
  main_plot <- base_embedding_plot(
    plot_bounds,
    coord_names = coord_names,
    axis_labels = axis_labels
  )
  main_plot <- add_cluster_outline_layers_external_labels(
    main_plot,
    clustered_best,
    bounds_data = plot_bounds,
    coord_names = coord_names,
    cluster_col = "cluster_id"
  )
  save_plot_pair(main_plot, file.path(figure_dir, output_stem), width = 7.4, height = 6.5)

  sub_dir <- file.path(figure_dir, "Subclusters")
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  split_vals <- sort(unique(as.character(subclustered_best$primary_cluster_id)))
  for (primary_id in split_vals) {
    d <- subclustered_best[as.character(subclustered_best$primary_cluster_id) == primary_id, , drop = FALSE]
    if (!nrow(d)) next
    sub_bounds <- add_plot_label_margin_rows(d, coord_names = coord_names, pad_frac = 0.2)
    p <- subcluster_embedding_plot(
      d,
      coord_names = coord_names,
      axis_labels = axis_labels,
      title = paste0(unique(d$method)[[1L]], " ", primary_id, " subclusters")
    )
    p <- add_cluster_outline_layers_external_labels(
      p,
      d,
      bounds_data = sub_bounds,
      coord_names = coord_names,
      cluster_col = "subcluster_id",
      linewidth = 0.72
    )
    save_plot_pair(
      p,
      file.path(sub_dir, paste0(reduction_file_suffix(reduction), "_", primary_id, "_subclusters")),
      width = 5.4,
      height = 4.8
    )
  }
  invisible(TRUE)
}


render_multi_warmup_landscape_figures <- function(analysis_root) {
  analysis_root <- normalizePath(path.expand(analysis_root), mustWork = TRUE)
  full_tables <- list.files(
    analysis_root, recursive = TRUE, full.names = TRUE,
    pattern = "_best_clusters_full_coordinates\\.csv$"
  )
  if (!length(full_tables)) stop("No materialized multi-warmup cluster coordinate tables under: ", analysis_root, call. = FALSE)
  for (full_path in full_tables) {
    table_dir <- dirname(full_path)
    stem <- sub("_full_coordinates\\.csv$", "", basename(full_path))
    best_path <- file.path(table_dir, paste0(stem, "_best_coordinates.csv"))
    subcluster_path <- file.path(table_dir, paste0(stem, "_best_subclusters.csv"))
    if (!file.exists(best_path) || !file.exists(subcluster_path)) next
    reduction <- if (grepl("_pca_", stem)) "pca" else if (grepl("_umap_", stem)) "umap" else "tsne"
    figure_dir <- sub("/Tables$", "/Figures", gsub("\\\\", "/", table_dir))
    write_cluster_figures(
      full_marked = utils::read.csv(full_path, check.names = FALSE, stringsAsFactors = FALSE),
      clustered_best = utils::read.csv(best_path, check.names = FALSE, stringsAsFactors = FALSE),
      subclustered_best = utils::read.csv(subcluster_path, check.names = FALSE, stringsAsFactors = FALSE),
      figure_dir = figure_dir,
      output_stem = stem,
      reduction = reduction,
      coord_names = reduction_coordinate_names(reduction)
    )
  }
  invisible(TRUE)
}
