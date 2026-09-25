#!/usr/bin/env Rscript

.o2pl_vis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "parameter_landscape_visualization_utils.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2pl_vis_workflow_root <- normalizePath(file.path(.o2pl_vis_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_vis_workflow_root, "util", "o2_supply_demand_map_parameter_landscape_io_utils.R"), local = environment(), chdir = TRUE)

o2pl_require_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required for parameter-landscape visualization.", call. = FALSE)
}

o2pl_save_figure <- function(figure, prefix, width = 6.5, height = 6.0) {
  o2pl_require_ggplot()
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(prefix, ".pdf"), figure, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(prefix, ".png"), figure, width = width, height = height, units = "in", dpi = 300)
  invisible(c(pdf = paste0(prefix, ".pdf"), png = paste0(prefix, ".png")))
}

o2pl_embedding_table_path <- function(root_dir, scope, reduction, preprocessing) {
  file.path(o2pl_analysis_tables_dir(root_dir, scope), o2pl_normalize_reduction(reduction), preprocessing, "embedding_coordinates.tsv")
}

o2pl_embedding_figure <- function(root_dir,
                                  scope,
                                  reduction = "umap",
                                  preprocessing = "zscore",
                                  initial_size = 0.22,
                                  fitted_size = 1.2) {
  o2pl_require_ggplot()
  path <- o2pl_embedding_table_path(root_dir, scope, reduction, preprocessing)
  if (!file.exists(path)) stop("Missing materialized analysis coordinates: ", path, call. = FALSE)
  data <- read_tsv(path)
  coordinates <- reduction_coordinate_names(reduction)
  o2pl_require_columns(data, c("point_type", coordinates), path)
  data$.x <- suppressWarnings(as.numeric(data[[coordinates[[1L]]]]))
  data$.y <- suppressWarnings(as.numeric(data[[coordinates[[2L]]]]))
  initial <- data[data$point_type == "initial", , drop = FALSE]
  fitted <- data[data$point_type == "fitted", , drop = FALSE]
  figure <- ggplot2::ggplot() +
    ggplot2::geom_point(data = initial, ggplot2::aes(x = .x, y = .y), color = "grey78", alpha = 0.25, size = initial_size) +
    ggplot2::geom_point(data = fitted, ggplot2::aes(x = .x, y = .y, color = objective), alpha = 0.95, size = fitted_size) +
    ggplot2::scale_color_gradient(low = "#2C7BB6", high = "#FDE725", na.value = "#333333", name = "Objective") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = reduction_axis_labels(reduction)[[1L]], y = reduction_axis_labels(reduction)[[2L]]) +
    ggplot2::theme_classic(base_size = 12)
  if ("cluster_id" %in% names(fitted) && any(!is.na(fitted$cluster_id))) {
    hulls <- do.call(rbind, lapply(split(fitted, fitted$cluster_id), function(group) {
      if (nrow(group) < 3L) return(NULL)
      group[chull(group$.x, group$.y), , drop = FALSE]
    }))
    if (!is.null(hulls) && nrow(hulls)) {
      figure <- figure + ggplot2::geom_polygon(data = hulls, ggplot2::aes(x = .x, y = .y, group = cluster_id), fill = NA, color = "black", linewidth = 0.45)
    }
  }
  canonical_prefix <- file.path(o2pl_vis_figures_dir(root_dir, scope), reduction, preprocessing, "embedding")
  outputs <- o2pl_save_figure(figure, canonical_prefix)
  dataset <- if (scope %in% c("invivo", "invitro")) scope else NULL
  if (!is.null(dataset)) {
    legacy_prefix_name <- paste(dataset, "deoptim", "initial", "vs", "best", reduction_file_suffix(reduction), sep = "_")
    if (preprocessing != "zscore") legacy_prefix_name <- paste(legacy_prefix_name, preprocessing, sep = "_")
    legacy_prefix <- file.path(paper_reduction_figures_dir(dataset, reduction, root_dir), "Full", legacy_prefix_name)
    dir.create(dirname(legacy_prefix), recursive = TRUE, showWarnings = FALSE)
    file.copy(outputs[["pdf"]], paste0(legacy_prefix, ".pdf"), overwrite = TRUE)
    file.copy(outputs[["png"]], paste0(legacy_prefix, ".png"), overwrite = TRUE)
    clustered_prefix <- file.path(paper_reduction_figures_wclusters_dir(dataset, reduction, root_dir), reduction_coordinate_cluster_dir(reduction), "Full", paste0(legacy_prefix_name, "_clustered"))
    dir.create(dirname(clustered_prefix), recursive = TRUE, showWarnings = FALSE)
    file.copy(outputs[["pdf"]], paste0(clustered_prefix, ".pdf"), overwrite = TRUE)
    file.copy(outputs[["png"]], paste0(clustered_prefix, ".png"), overwrite = TRUE)
  } else {
    legacy_prefix_name <- paste("pooled", "invivo", "invitro", "initial", "vs", "best", reduction_file_suffix(reduction), sep = "_")
    if (preprocessing != "zscore") legacy_prefix_name <- paste(legacy_prefix_name, preprocessing, sep = "_")
    legacy_prefix <- file.path(paper_pooled_reduction_figures_dir(reduction, root_dir), "Full", legacy_prefix_name)
    dir.create(dirname(legacy_prefix), recursive = TRUE, showWarnings = FALSE)
    file.copy(outputs[["pdf"]], paste0(legacy_prefix, ".pdf"), overwrite = TRUE)
    file.copy(outputs[["png"]], paste0(legacy_prefix, ".png"), overwrite = TRUE)
  }
  o2pl_record_artifact(root_dir, "visualization", paste(scope, reduction, preprocessing, "embedding_figure", sep = "_"), outputs[["pdf"]], NULL, NA_character_, "visualize_parameter_landscape.R", scope, path)
  invisible(outputs)
}

o2pl_contribution_figure <- function(root_dir, target = "mode", reference_o2 = 2, top_n = 20L) {
  o2pl_require_ggplot()
  target <- if (tolower(target) %in% c("mode", "discrete")) "mode" else "dominant_ploidy"
  base_dir <- file.path(root_dir, if (target == "mode") "mode_parameter_contribution" else "dominant_ploidy_parameter_contribution")
  reference_dir <- file.path(base_dir, paste0("reference_o2_", fixed_o2_o2_slug(reference_o2)))
  filename <- if (target == "mode") "mode_parameter_feature_importance.csv" else "dominant_ploidy_parameter_feature_importance.csv"
  path <- file.path(reference_dir, filename)
  if (!file.exists(path)) stop("Missing materialized contribution table: ", path, call. = FALSE)
  data <- head(read_csv_plain(path), as.integer(top_n))
  o2pl_require_columns(data, c("feature_name", "contribution_score"), path)
  data$feature_name <- factor(data$feature_name, levels = rev(data$feature_name))
  figure <- ggplot2::ggplot(data, ggplot2::aes(x = contribution_score, y = feature_name, fill = feature_type)) +
    ggplot2::geom_col(width = 0.76) +
    ggplot2::scale_fill_manual(values = c(main = "#2C7BB6", interaction = "#D7191C"), na.value = "grey55") +
    ggplot2::labs(x = if (target == "mode") "Discriminative AUC" else "Explained variance contribution", y = NULL, fill = "Feature type") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  prefix <- file.path(reference_dir, "figures", paste0(target, "_parameter_feature_importance"))
  outputs <- o2pl_save_figure(figure, prefix, width = 7.2, height = 6.2)
  o2pl_record_artifact(root_dir, "visualization", paste(target, fixed_o2_o2_slug(reference_o2), "contribution_figure", sep = "_"), outputs[["pdf"]], NULL, NA_character_, "visualize_parameter_landscape.R", "invivo", path)
  invisible(outputs)
}

o2pl_fixed_o2_distribution_figure <- function(root_dir) {
  o2pl_require_ggplot()
  path <- file.path(paper_fixo2_mode_tables_dir(root_dir), "fixed_o2_attractor_mode_by_seed_o2.tsv")
  if (!file.exists(path)) stop("Missing fixed-O2 materialized table: ", path, call. = FALSE)
  data <- read_tsv(path)
  o2pl_require_columns(data, c("O2_pct", "dominant_mean_ploidy"), path)
  figure <- ggplot2::ggplot(data, ggplot2::aes(x = factor(O2_pct), y = dominant_mean_ploidy)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.15, fill = "#56B4E9") +
    ggplot2::labs(x = expression(O[2] * " (%)"), y = "Dominant mean ploidy") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  prefix <- file.path(o2pl_vis_figures_dir(root_dir, "fixed_o2"), "dominant_ploidy_distribution")
  outputs <- o2pl_save_figure(figure, prefix, width = 8.2, height = 5.2)
  o2pl_record_artifact(root_dir, "visualization", "fixed_o2_dominant_ploidy_distribution", outputs[["pdf"]], NULL, NA_character_, "visualize_parameter_landscape.R", "invivo", path)
  invisible(outputs)
}
