#!/usr/bin/env Rscript

# Pure visualization consumer for materialized FixO2 eigen coordinate tables.

.fixo2ea_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "render_fixo2_eigen_attractor_figures.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  workflow_frames <- frames[
    grepl("/O2_supply_demand_MAP/", frames, fixed = TRUE)
  ]
  if (length(workflow_frames)) {
    root <- dirname(workflow_frames[[length(workflow_frames)]])
    while (!identical(basename(root), "O2_supply_demand_MAP")) {
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    if (identical(basename(root), "O2_supply_demand_MAP")) {
      return(file.path(root, "vis", "fixed_o2_eigen"))
    }
  }
  normalizePath(getwd(), mustWork = FALSE)
})
.fixo2ea_vis_workflow_root <- normalizePath(file.path(.fixo2ea_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.fixo2ea_vis_workflow_root, "util", "o2_supply_demand_map_eigen_attractor_utils.R"), local = environment())

fixo2ea_coordinate_columns <- function(df) {
  for (columns in list(c("PC1", "PC2"), c("UMAP1", "UMAP2"), c("TSNE1", "TSNE2"))) {
    if (all(columns %in% names(df))) return(columns)
  }
  stop("Coordinate table lacks a recognized two-dimensional coordinate pair.", call. = FALSE)
}

fixo2ea_figure_path_from_table <- function(path, extension) {
  normalized <- gsub("\\\\", "/", normalizePath(path, mustWork = TRUE))
  normalized <- sub("/TablesWclusters/", "/FiguresWclusters/", normalized, fixed = TRUE)
  normalized <- sub("/Tables/", "/Figures/", normalized, fixed = TRUE)
  sub("_coordinates\\.csv$", paste0(".", extension), normalized)
}

fixo2ea_embedding_plot <- function(df, coordinate_columns, title) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
  x <- coordinate_columns[[1L]]
  y <- coordinate_columns[[2L]]
  if ("cluster_id" %in% names(df)) {
    df$cluster_id <- factor(df$cluster_id, levels = sort(unique(as.character(df$cluster_id))))
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data$cluster_id)) +
        ggplot2::geom_point(alpha = 0.68, size = 0.75) +
        ggplot2::coord_equal() +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::labs(x = x, y = y, color = "Cluster", title = title)
    )
  }
  if (!"point_type" %in% names(df)) df$point_type <- "point"
  if (!"dataset" %in% names(df)) df$dataset <- "dataset"
  df$point_type <- factor(df$point_type, levels = unique(c("initial", "best", as.character(df$point_type))))
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data$dataset, shape = .data$point_type)) +
    ggplot2::geom_point(alpha = 0.58, size = 0.8) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = x, y = y, color = "Dataset", shape = "Point type", title = title)
}

fixo2ea_render_coordinate_table <- function(path, width = 6.5, height = 6.2, dpi = 220) {
  df <- read_csv_plain(path)
  coordinates <- fixo2ea_coordinate_columns(df)
  title <- gsub("_", " ", sub("_coordinates\\.csv$", "", basename(path)))
  figure <- fixo2ea_embedding_plot(df, coordinates, title)
  pdf_path <- fixo2ea_figure_path_from_table(path, "pdf")
  png_path <- fixo2ea_figure_path_from_table(path, "png")
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(pdf_path, figure, width = width, height = height, useDingbats = FALSE)
  ggplot2::ggsave(png_path, figure, width = width, height = height, dpi = dpi)
  invisible(c(pdf = pdf_path, png = png_path))
}

fixo2ea_render_all_figures <- function(result_root = fixo2ea_default_result_root(),
                                       width = 6.5,
                                       height = 6.2,
                                       dpi = 220) {
  roots <- file.path(result_root, c("invivo", "invitro", "pooled_invivo_invitro"))
  tables <- unlist(lapply(roots[dir.exists(roots)], function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE, pattern = "_coordinates\\.csv$")
  }), use.names = FALSE)
  tables <- tables[grepl("/(Tables|TablesWclusters)/", gsub("\\\\", "/", tables))]
  if (!length(tables)) stop("No materialized FixO2 eigen coordinate tables found under: ", result_root, call. = FALSE)
  outputs <- lapply(sort(unique(tables)), fixo2ea_render_coordinate_table,
                    width = width, height = height, dpi = dpi)
  manifest <- data.frame(
    coordinate_table = sort(unique(tables)),
    pdf = vapply(outputs, `[[`, character(1L), "pdf"),
    png = vapply(outputs, `[[`, character(1L), "png"),
    stringsAsFactors = FALSE
  )
  write_csv(manifest, file.path(result_root, "FixO2EigenAttractors", "Tables", "fixo2_eigen_figure_manifest.csv"))
  invisible(manifest)
}
