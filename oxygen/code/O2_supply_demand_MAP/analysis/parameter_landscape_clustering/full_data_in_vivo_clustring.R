#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
source(file.path(SCRIPT_DIR, "parameter_landscape_utils.R"))

default_input_pdf <- function(root_dir, reduction = "tsne") {
  reduction <- normalize_reduction(reduction)
  file.path(
    normalizePath(path.expand(root_dir), mustWork = FALSE),
    "pooled_invivo_invitro", reduction_dir_name(reduction), "Figures", "Full",
    paste0("pooled_invivo_invitro_initial_vs_best_", reduction_file_suffix(reduction), ".pdf")
  )
}

infer_coordinate_csv <- function(pdf_path) {
  coord_path <- gsub("/Figures/", "/Tables/", pdf_path, fixed = TRUE)
  coord_path <- sub("\\.pdf$", "_coordinates.csv", coord_path)
  if (identical(coord_path, pdf_path)) {
    stop("Could not infer coordinate CSV from PDF path: ", pdf_path, call. = FALSE)
  }
  coord_path
}

cluster_dataset_specs <- function() {
  data.frame(
    dataset = c("invivo", "invitro"),
    dataset_label = c("in vivo", "in vitro"),
    output_token = c("invivo", "invitro"),
    cluster_prefix = c("vi", "vt"),
    stringsAsFactors = FALSE
  )
}

summarize_best_clusters <- function(clustered_best, coord_names) {
  rows <- lapply(sort(unique(clustered_best$cluster_id)), function(cluster_id) {
    d <- clustered_best[clustered_best$cluster_id == cluster_id, , drop = FALSE]
    base <- data.frame(
      dataset = unique(d$dataset)[[1L]],
      dataset_label = unique(d$dataset_label)[[1L]],
      cluster_prefix = unique(d$cluster_prefix)[[1L]],
      cluster_id = cluster_id,
      cluster_base_id = unique(d$cluster_base_id)[[1L]],
      cluster_num = unique(d$cluster_num)[[1L]],
      n_seeds = nrow(d),
      seed_min = min(d$seed, na.rm = TRUE),
      seed_max = max(d$seed, na.rm = TRUE),
      objective_mean = mean(d$objective, na.rm = TRUE),
      objective_median = stats::median(d$objective, na.rm = TRUE),
      objective_min = min(d$objective, na.rm = TRUE),
      objective_max = max(d$objective, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    coord_values <- c(
      stats::median(d[[coord_names[[1L]]]], na.rm = TRUE),
      stats::median(d[[coord_names[[2L]]]], na.rm = TRUE),
      min(d[[coord_names[[1L]]]], na.rm = TRUE),
      max(d[[coord_names[[1L]]]], na.rm = TRUE),
      min(d[[coord_names[[2L]]]], na.rm = TRUE),
      max(d[[coord_names[[2L]]]], na.rm = TRUE)
    )
    coord_out <- as.data.frame(
      as.list(stats::setNames(
        coord_values,
        c(
          paste0(coord_names, "_median"),
          paste0(coord_names[[1L]], c("_min", "_max")),
          paste0(coord_names[[2L]], c("_min", "_max"))
        )
      )),
      check.names = FALSE
    )
    cbind(base, coord_out)
  })
  out <- do.call(rbind, rows)
  out[order(out$cluster_num), , drop = FALSE]
}

best_seed_group_table <- function(clustered_best, reduction, coord_names) {
  out <- data.frame(
    method = reduction,
    reduction = reduction,
    dataset = as.character(clustered_best$dataset),
    dataset_label = as.character(clustered_best$dataset_label),
    cluster_prefix = as.character(clustered_best$cluster_prefix),
    cluster_id = as.character(clustered_best$cluster_id),
    cluster_base_id = as.character(clustered_best$cluster_base_id),
    cluster_num = as.integer(clustered_best$cluster_num),
    seed = as.integer(clustered_best$seed),
    objective = suppressWarnings(as.numeric(clustered_best$objective)),
    coordinate_1_name = coord_names[[1L]],
    coordinate_2_name = coord_names[[2L]],
    stringsAsFactors = FALSE
  )
  out$coordinate_1 <- suppressWarnings(as.numeric(clustered_best[[coord_names[[1L]]]]))
  out$coordinate_2 <- suppressWarnings(as.numeric(clustered_best[[coord_names[[2L]]]]))
  out <- out[order(out$method, out$dataset, out$cluster_num, -out$objective, out$seed), , drop = FALSE]
  row.names(out) <- NULL
  out
}

seed_key <- function(dataset, seed) {
  paste(as.character(dataset), as.integer(seed), sep = "::")
}

load_pooled_best_feature_data <- function(root_dir) {
  params <- pooled_umap_parameter_set()
  log10_params <- pooled_umap_log10_parameter_set()
  read_best <- function(dataset) {
    path <- paper_best_params_csv(dataset, root_dir = root_dir)
    if (!file.exists(path)) stop("Missing best-parameter table for ", dataset, ": ", path, call. = FALSE)
    df <- read_csv_plain(path)
    required <- c("seed", params)
    missing <- setdiff(required, names(df))
    if (length(missing)) {
      stop("Best-parameter table for ", dataset, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    df$dataset <- dataset
    df$seed <- suppressWarnings(as.integer(df$seed))
    if (!"objective" %in% names(df)) df$objective <- NA_real_
    df
  }

  best_df <- rbind_fill_plain(list(read_best("invivo"), read_best("invitro")))
  feature_df <- transform_umap_features(best_df, params, log10_params)
  feature_df <- cbind(
    data.frame(
      dataset = as.character(best_df$dataset),
      seed = as.integer(best_df$seed),
      stringsAsFactors = FALSE
    ),
    feature_df
  )
  feature_df$.seed_key <- seed_key(feature_df$dataset, feature_df$seed)
  list(
    params = params,
    log10_params = log10_params,
    feature_df = feature_df
  )
}

align_best_feature_rows <- function(clustered_best, feature_data) {
  keys <- seed_key(clustered_best$dataset, clustered_best$seed)
  idx <- match(keys, feature_data$feature_df$.seed_key)
  if (any(is.na(idx))) {
    missing <- keys[is.na(idx)]
    stop("Missing raw feature rows for clustered best seeds: ", paste(head(missing, 20), collapse = ", "), call. = FALSE)
  }
  feature_data$feature_df[idx, c(feature_data$params), drop = FALSE]
}

zscore_primary_features <- function(raw_features,
                                    params,
                                    log10_params,
                                    method,
                                    dataset,
                                    primary_cluster_id) {
  mat <- as.matrix(raw_features[, params, drop = FALSE])
  storage.mode(mat) <- "double"
  centers <- colMeans(mat, na.rm = TRUE)
  scales <- apply(mat, 2L, stats::sd, na.rm = TRUE)
  used <- is.finite(centers) & is.finite(scales) & scales > 0
  metadata <- data.frame(
    method = method,
    reduction = method,
    dataset = dataset,
    primary_cluster_id = primary_cluster_id,
    feature = params,
    log10_transformed = params %in% log10_params,
    center = as.numeric(centers),
    scale = as.numeric(scales),
    used = as.logical(used),
    drop_reason = ifelse(used, "", "zero_or_nonfinite_sd"),
    stringsAsFactors = FALSE
  )
  if (!any(used)) {
    return(list(mat = matrix(nrow = nrow(mat), ncol = 0L), metadata = metadata))
  }
  scaled <- scale(mat[, used, drop = FALSE], center = centers[used], scale = scales[used])
  colnames(scaled) <- params[used]
  list(mat = scaled, metadata = metadata)
}

summarize_best_subclusters <- function(subclustered_best, coord_names) {
  if (!nrow(subclustered_best)) return(data.frame())
  rows <- lapply(sort(unique(subclustered_best$subcluster_id)), function(subcluster_id) {
    d <- subclustered_best[subclustered_best$subcluster_id == subcluster_id, , drop = FALSE]
    min_i <- which.min(d$objective)
    base <- data.frame(
      method = unique(d$method)[[1L]],
      reduction = unique(d$reduction)[[1L]],
      dataset = unique(d$dataset)[[1L]],
      dataset_label = unique(d$dataset_label)[[1L]],
      primary_cluster_id = unique(d$primary_cluster_id)[[1L]],
      primary_cluster_base_id = unique(d$primary_cluster_base_id)[[1L]],
      primary_cluster_num = unique(d$primary_cluster_num)[[1L]],
      subcluster_id = subcluster_id,
      subcluster_base_id = unique(d$subcluster_base_id)[[1L]],
      subcluster_num = unique(d$subcluster_num)[[1L]],
      n_seeds = nrow(d),
      objective_min_seed = d$seed[[min_i]],
      objective_mean = mean(d$objective, na.rm = TRUE),
      objective_median = stats::median(d$objective, na.rm = TRUE),
      objective_min = min(d$objective, na.rm = TRUE),
      objective_max = max(d$objective, na.rm = TRUE),
      subcluster_k = unique(d$subcluster_k)[[1L]],
      subcluster_silhouette_avg = unique(d$subcluster_silhouette_avg)[[1L]],
      subcluster_silhouette_sample_n = unique(d$subcluster_silhouette_sample_n)[[1L]],
      stringsAsFactors = FALSE
    )
    coord_values <- c(
      stats::median(d[[coord_names[[1L]]]], na.rm = TRUE),
      stats::median(d[[coord_names[[2L]]]], na.rm = TRUE),
      min(d[[coord_names[[1L]]]], na.rm = TRUE),
      max(d[[coord_names[[1L]]]], na.rm = TRUE),
      min(d[[coord_names[[2L]]]], na.rm = TRUE),
      max(d[[coord_names[[2L]]]], na.rm = TRUE)
    )
    coord_out <- as.data.frame(
      as.list(stats::setNames(
        coord_values,
        c(
          paste0(coord_names, "_median"),
          paste0(coord_names[[1L]], c("_min", "_max")),
          paste0(coord_names[[2L]], c("_min", "_max"))
        )
      )),
      check.names = FALSE
    )
    cbind(base, coord_out)
  })
  out <- do.call(rbind, rows)
  out[order(out$method, out$dataset, out$primary_cluster_num, out$subcluster_num), , drop = FALSE]
}

subcluster_seed_group_table <- function(subclustered_best, coord_names) {
  out <- data.frame(
    method = as.character(subclustered_best$method),
    reduction = as.character(subclustered_best$reduction),
    dataset = as.character(subclustered_best$dataset),
    dataset_label = as.character(subclustered_best$dataset_label),
    primary_cluster_id = as.character(subclustered_best$primary_cluster_id),
    primary_cluster_base_id = as.character(subclustered_best$primary_cluster_base_id),
    primary_cluster_num = as.integer(subclustered_best$primary_cluster_num),
    subcluster_id = as.character(subclustered_best$subcluster_id),
    subcluster_base_id = as.character(subclustered_best$subcluster_base_id),
    subcluster_num = as.integer(subclustered_best$subcluster_num),
    seed = as.integer(subclustered_best$seed),
    objective = suppressWarnings(as.numeric(subclustered_best$objective)),
    coordinate_1_name = coord_names[[1L]],
    coordinate_2_name = coord_names[[2L]],
    stringsAsFactors = FALSE
  )
  out$coordinate_1 <- suppressWarnings(as.numeric(subclustered_best[[coord_names[[1L]]]]))
  out$coordinate_2 <- suppressWarnings(as.numeric(subclustered_best[[coord_names[[2L]]]]))
  out <- out[order(out$method, out$dataset, out$primary_cluster_num, out$subcluster_num, -out$objective, out$seed), , drop = FALSE]
  row.names(out) <- NULL
  out
}

subcluster_zscore_feature_table <- function(subclustered_best, zscore_features_by_primary, coord_names) {
  id_cols <- subcluster_seed_group_table(subclustered_best, coord_names = coord_names)
  feature_cols <- unique(unlist(lapply(zscore_features_by_primary, colnames), use.names = FALSE))
  if (!length(feature_cols)) return(id_cols)
  feature_out <- as.data.frame(matrix(NA_real_, nrow = nrow(subclustered_best), ncol = length(feature_cols)))
  names(feature_out) <- paste0("z_", feature_cols)
  cursor <- 1L
  for (feature_mat in zscore_features_by_primary) {
    n <- nrow(feature_mat)
    if (n > 0L && ncol(feature_mat) > 0L) {
      feature_out[cursor:(cursor + n - 1L), paste0("z_", colnames(feature_mat))] <- as.data.frame(feature_mat)
    }
    cursor <- cursor + n
  }
  cbind(id_cols, feature_out)
}

build_subcluster_plot <- function(primary_data,
                                  coord_names = c("UMAP1", "UMAP2"),
                                  axis_labels = c("UMAP 1", "UMAP 2"),
                                  label_margin_frac = 0.2) {
  primary_data$.embedding_x <- suppressWarnings(as.numeric(primary_data[[coord_names[[1L]]]]))
  primary_data$.embedding_y <- suppressWarnings(as.numeric(primary_data[[coord_names[[2L]]]]))
  bounds_data <- add_plot_label_margin_rows(primary_data, coord_names = coord_names, pad_frac = label_margin_frac)
  lims <- square_umap_limits(bounds_data, coord_names = coord_names)
  pal <- cluster_palette(primary_data$subcluster_id)
  p <- ggplot2::ggplot(primary_data, ggplot2::aes(x = .embedding_x, y = .embedding_y)) +
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
    ggplot2::labs(
      x = axis_labels[[1L]],
      y = axis_labels[[2L]],
      title = paste0(unique(primary_data$method)[[1L]], " ", unique(primary_data$primary_cluster_id)[[1L]], " subclusters")
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
  add_cluster_outline_layers_external_labels(
    p,
    primary_data,
    bounds_data = bounds_data,
    coord_names = coord_names,
    cluster_col = "subcluster_id"
  )
}

write_subcluster_plots <- function(subclustered_best,
                                   figure_dir,
                                   reduction_suffix,
                                   coord_names,
                                   axis_labels,
                                   label_margin_frac = 0.2) {
  subcluster_dir <- file.path(figure_dir, "Subclusters")
  dir.create(subcluster_dir, recursive = TRUE, showWarnings = FALSE)
  primary_ids <- unique(subclustered_best$primary_cluster_id[order(subclustered_best$dataset, subclustered_best$primary_cluster_num)])
  for (primary_cluster_id in primary_ids) {
    primary_data <- subclustered_best[subclustered_best$primary_cluster_id == primary_cluster_id, , drop = FALSE]
    output_prefix <- file.path(subcluster_dir, paste0(reduction_suffix, "_", primary_cluster_id, "_subclusters"))
    save_plot_pair(
      build_subcluster_plot(
        primary_data,
        coord_names = coord_names,
        axis_labels = axis_labels,
        label_margin_frac = label_margin_frac
      ),
      output_prefix,
      width = 5.4,
      height = 4.8
    )
  }
  invisible(subcluster_dir)
}

analyze_best_subclusters <- function(clustered_best,
                                     feature_data,
                                     reduction,
                                     coord_names,
                                     subcluster_seed = 1123L,
                                     subcluster_k_min = 2L,
                                     subcluster_k_max = 6L,
                                     subcluster_min_n = 6L,
                                     silhouette_sample_n = 5000L) {
  feature_rows <- align_best_feature_rows(clustered_best, feature_data)
  clustered_best$.row_id <- seq_len(nrow(clustered_best))
  split_ids <- split(
    clustered_best$.row_id,
    interaction(clustered_best$dataset, clustered_best$cluster_id, drop = TRUE, lex.order = TRUE)
  )

  subclustered <- list()
  silhouettes <- list()
  feature_meta <- list()
  zscore_features <- list()
  cursor <- 0L
  for (i in seq_along(split_ids)) {
    idx <- split_ids[[i]]
    d <- clustered_best[idx, , drop = FALSE]
    dataset <- unique(as.character(d$dataset))[[1L]]
    primary_cluster_id <- unique(as.character(d$cluster_id))[[1L]]
    raw <- feature_rows[idx, , drop = FALSE]
    z <- zscore_primary_features(
      raw_features = raw,
      params = feature_data$params,
      log10_params = feature_data$log10_params,
      method = reduction,
      dataset = dataset,
      primary_cluster_id = primary_cluster_id
    )
    source <- paste0(reduction, "_", primary_cluster_id, "_best_raw_features_zscore")
    if (nrow(d) < as.integer(subcluster_min_n) || ncol(z$mat) == 0L) {
      assignment <- single_cluster_assignment(nrow(d), source, nrow(d))
    } else {
      assignment <- auto_silhouette_kmeans(
        basis_mat = z$mat,
        plot_data = d,
        cluster_source = source,
        coord_names = coord_names,
        seed = as.integer(subcluster_seed) + i - 1L,
        k_min = subcluster_k_min,
        k_max = min(as.integer(subcluster_k_max), nrow(d) - 1L),
        silhouette_sample_n = silhouette_sample_n
      )
    }
    selected_summary <- assignment$summary[assignment$summary$selected, , drop = FALSE]
    d$method <- reduction
    d$reduction <- reduction
    d$primary_cluster_id <- as.character(d$cluster_id)
    d$primary_cluster_base_id <- as.character(d$cluster_base_id)
    d$primary_cluster_num <- as.integer(d$cluster_num)
    d$subcluster_source <- source
    d$subcluster_base_id <- sprintf("Sc%02d", assignment$cluster_num)
    d$subcluster_id <- paste0(d$primary_cluster_id, "_", d$subcluster_base_id)
    d$subcluster_num <- as.integer(assignment$cluster_num)
    d$subcluster_k <- length(unique(assignment$cluster_num))
    d$subcluster_silhouette_avg <- selected_summary$average_silhouette[[1L]]
    d$subcluster_silhouette_sample_n <- selected_summary$sample_n[[1L]]
    d$.row_id <- NULL

    sil <- assignment$summary
    sil$method <- reduction
    sil$reduction <- reduction
    sil$dataset <- dataset
    sil$primary_cluster_id <- primary_cluster_id
    sil$primary_cluster_n <- nrow(d)
    sil <- sil[, c("method", "reduction", "dataset", "primary_cluster_id", "primary_cluster_n", setdiff(names(sil), c("method", "reduction", "dataset", "primary_cluster_id", "primary_cluster_n"))), drop = FALSE]

    cursor <- cursor + 1L
    subclustered[[cursor]] <- d
    silhouettes[[cursor]] <- sil
    feature_meta[[cursor]] <- z$metadata
    zscore_features[[cursor]] <- z$mat
  }

  subclustered_best <- do.call(rbind, subclustered)
  row.names(subclustered_best) <- NULL
  zscore_table <- subcluster_zscore_feature_table(subclustered_best, zscore_features, coord_names = coord_names)
  subclustered_best <- subclustered_best[order(subclustered_best$dataset, subclustered_best$primary_cluster_num, subclustered_best$subcluster_num, subclustered_best$seed), , drop = FALSE]
  row.names(subclustered_best) <- NULL
  zscore_table <- zscore_table[order(zscore_table$dataset, zscore_table$primary_cluster_num, zscore_table$subcluster_num, zscore_table$seed), , drop = FALSE]
  row.names(zscore_table) <- NULL
  list(
    best = subclustered_best,
    seed_groups = subcluster_seed_group_table(subclustered_best, coord_names = coord_names),
    summary = summarize_best_subclusters(subclustered_best, coord_names = coord_names),
    silhouette = do.call(rbind, silhouettes),
    feature_metadata = do.call(rbind, feature_meta),
    zscore_features = zscore_table
  )
}

external_cluster_label_data <- function(clustered_plot_data,
                                        bounds_data,
                                        coord_names = c("UMAP1", "UMAP2"),
                                        cluster_col = "cluster_id",
                                        offset_frac = 0.07) {
  coord_names <- as.character(coord_names)
  lims <- square_umap_limits(bounds_data, coord_names = coord_names)
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
  row.names(out) <- NULL
  out
}

add_plot_label_margin_rows <- function(plot_data,
                                       coord_names = c("UMAP1", "UMAP2"),
                                       pad_frac = 0.14) {
  coord_names <- as.character(coord_names)
  pad_frac <- suppressWarnings(as.numeric(pad_frac))
  if (!is.finite(pad_frac) || pad_frac < 0) pad_frac <- 0.14
  lims <- square_umap_limits(plot_data, pad_frac = pad_frac, coord_names = coord_names)
  margin_rows <- plot_data[rep(1L, 4L), , drop = FALSE]
  for (nm in names(margin_rows)) margin_rows[[nm]] <- NA
  margin_rows[[coord_names[[1L]]]] <- c(lims$xlim[[1L]], lims$xlim[[2L]], lims$xlim[[1L]], lims$xlim[[2L]])
  margin_rows[[coord_names[[2L]]]] <- c(lims$ylim[[1L]], lims$ylim[[1L]], lims$ylim[[2L]], lims$ylim[[2L]])
  if ("point_type" %in% names(margin_rows)) margin_rows$point_type <- "plot_margin"
  if ("dataset" %in% names(margin_rows)) margin_rows$dataset <- "plot_margin"
  rbind(plot_data, margin_rows)
}

add_cluster_outline_layers_external_labels <- function(plot,
                                                       clustered_plot_data,
                                                       bounds_data,
                                                       coord_names = c("UMAP1", "UMAP2"),
                                                       cluster_col = "cluster_id") {
  hulls <- cluster_hull_data(clustered_plot_data, cluster_col = cluster_col, coord_names = coord_names)
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
        linewidth = 0.72,
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
  lims <- square_umap_limits(bounds_data, coord_names = coord_names)

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

cluster_best_dataset <- function(plot_data,
                                 dataset,
                                 dataset_label,
                                 cluster_prefix,
                                 coord_names,
                                 cluster_seed,
                                 cluster_k_min,
                                 cluster_k_max,
                                 silhouette_sample_n) {
  best_idx <- which(plot_data$dataset == dataset & plot_data$point_type == "best")
  if (!length(best_idx)) {
    stop("No rows matched dataset == '", dataset, "' and point_type == 'best'.", call. = FALSE)
  }
  best <- plot_data[best_idx, , drop = FALSE]
  if (anyDuplicated(best$seed)) {
    warning("Duplicated ", dataset_label, " best seeds were found in the coordinate table.")
  }

  cluster_source <- paste0(dataset, "_best_", paste(coord_names, collapse = "_"))
  cluster_assignment <- auto_silhouette_kmeans(
    basis_mat = as.matrix(best[, coord_names, drop = FALSE]),
    plot_data = best,
    cluster_source = cluster_source,
    coord_names = coord_names,
    seed = cluster_seed,
    k_min = cluster_k_min,
    k_max = cluster_k_max,
    silhouette_sample_n = silhouette_sample_n
  )
  selected_summary <- cluster_assignment$summary[cluster_assignment$summary$selected, , drop = FALSE]
  if (!nrow(selected_summary)) stop("No selected clustering summary row was produced for ", dataset_label, ".", call. = FALSE)

  best$dataset_label <- dataset_label
  best$cluster_scope <- paste0(dataset, "_best")
  best$cluster_source <- cluster_source
  best$cluster_prefix <- cluster_prefix
  best$cluster_base_id <- cluster_assignment$cluster_id
  best$cluster_id <- paste0(cluster_prefix, "_", cluster_assignment$cluster_id)
  best$cluster_num <- cluster_assignment$cluster_num
  best$cluster_k <- length(unique(cluster_assignment$cluster_num))
  best$cluster_silhouette_avg <- selected_summary$average_silhouette[[1L]]
  best$cluster_silhouette_sample_n <- selected_summary$sample_n[[1L]]

  silhouette <- cluster_assignment$summary
  silhouette$dataset <- dataset
  silhouette$dataset_label <- dataset_label
  silhouette$cluster_prefix <- cluster_prefix
  silhouette$cluster_source <- cluster_source

  summary <- summarize_best_clusters(best, coord_names = coord_names)
  summary$cluster_source <- cluster_source
  summary$selected_k <- selected_summary$k[[1L]]
  summary$selected_average_silhouette <- selected_summary$average_silhouette[[1L]]

  list(
    dataset = dataset,
    dataset_label = dataset_label,
    cluster_prefix = cluster_prefix,
    best_idx = best_idx,
    best = best,
    silhouette = silhouette,
    summary = summary,
    selected_summary = selected_summary
  )
}

argv_value <- function(argv, key, default = NULL) {
  value <- argv[[key]]
  if (is.null(value) || length(value) == 0L || all(is.na(value))) return(default)
  value
}

cleanup_obsolete_outputs <- function(output_dir, reduction_suffix) {
  figures_dir <- file.path(output_dir, "Figures")
  tables_dir <- file.path(output_dir, "Tables")
  stems <- c(
    paste0("pooled_invivo_invitro_initial_vs_best_", reduction_suffix, "_best_clusters_best_only"),
    paste0("pooled_invivo_invitro_initial_vs_best_", reduction_suffix, "_invivo_best_clusters"),
    paste0("pooled_invivo_invitro_initial_vs_best_", reduction_suffix, "_invitro_best_clusters")
  )
  stale <- character()
  if (dir.exists(figures_dir)) {
    stale <- c(stale, unlist(lapply(stems, function(stem) {
      list.files(figures_dir, pattern = paste0("^", stem, ".*\\.(pdf|png)$"), full.names = TRUE)
    }), use.names = FALSE))
  }
  if (dir.exists(tables_dir)) {
    stale <- c(stale, unlist(lapply(stems[-1L], function(stem) {
      list.files(tables_dir, pattern = paste0("^", stem, ".*\\.csv$"), full.names = TRUE)
    }), use.names = FALSE))
  }
  stale <- unique(stale[file.exists(stale)])
  if (length(stale)) {
    unlink(stale)
    message("Removed obsolete separate-output files: ", length(stale))
  }
  invisible(stale)
}

analyze_embedding <- function(reduction, argv, root_dir, output_dir) {
  reduction <- normalize_reduction(reduction)
  coord_names <- reduction_coordinate_names(reduction)
  axis_labels <- reduction_axis_labels(reduction)
  reduction_suffix <- reduction_file_suffix(reduction)
  single_reduction_args <- length(as_char_vec(argv$reductions %||% argv$reduction, c("tsne", "umap"))) == 1L
  pdf_default <- default_input_pdf(root_dir, reduction)
  if (isTRUE(single_reduction_args) && !is.null(argv$pdf_path) && length(argv$pdf_path) > 0L && !is.na(argv$pdf_path)) {
    pdf_default <- argv$pdf_path
  }
  pdf_path <- normalizePath(
    path.expand(argv_value(argv, paste0(reduction_suffix, "_pdf_path"), pdf_default)),
    mustWork = FALSE
  )
  coordinate_default <- infer_coordinate_csv(pdf_path)
  if (isTRUE(single_reduction_args) && !is.null(argv$coordinate_csv) &&
      length(argv$coordinate_csv) > 0L && !is.na(argv$coordinate_csv)) {
    coordinate_default <- argv$coordinate_csv
  }
  coordinate_csv <- normalizePath(
    path.expand(argv_value(argv, paste0(reduction_suffix, "_coordinate_csv"), coordinate_default)),
    mustWork = FALSE
  )
  cluster_seed <- as_int(argv$cluster_seed, 123L)
  cluster_k_min <- as_int(argv$cluster_k_min, 2L)
  cluster_k_max <- as_int(argv$cluster_k_max, 8L)
  silhouette_sample_n <- as_int(argv$cluster_silhouette_sample_n, 5000L)
  label_margin_frac <- as_num(argv$label_margin_frac, 0.14)
  subcluster_enabled <- as_bool(argv$subcluster, TRUE)
  subcluster_seed <- as_int(argv$subcluster_seed, cluster_seed + 1000L)
  subcluster_k_min <- as_int(argv$subcluster_k_min, 2L)
  subcluster_k_max <- as_int(argv$subcluster_k_max, 6L)
  subcluster_min_n <- as_int(argv$subcluster_min_n, 6L)
  subcluster_label_margin_frac <- as_num(argv$subcluster_label_margin_frac, 0.2)
  specs <- cluster_dataset_specs()

  if (!file.exists(pdf_path)) stop("Missing input PDF: ", pdf_path, call. = FALSE)
  if (!file.exists(coordinate_csv)) stop("Missing coordinate CSV: ", coordinate_csv, call. = FALSE)

  message("")
  message("Reduction: ", reduction)
  message("Input PDF: ", pdf_path)
  message("Coordinate CSV: ", coordinate_csv)
  message("Output directory: ", output_dir)

  plot_data <- utils::read.csv(coordinate_csv, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(coord_names, "dataset", "point_type", "seed", "objective")
  missing <- setdiff(required, names(plot_data))
  if (length(missing)) {
    stop("Coordinate CSV is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  plot_data$dataset <- as.character(plot_data$dataset)
  plot_data$point_type <- as.character(plot_data$point_type)
  plot_data$seed <- suppressWarnings(as.integer(plot_data$seed))
  plot_data$objective <- suppressWarnings(as.numeric(plot_data$objective))
  for (nm in coord_names) {
    plot_data[[nm]] <- suppressWarnings(as.numeric(plot_data[[nm]]))
  }
  finite_coords <- stats::complete.cases(plot_data[, coord_names, drop = FALSE])
  if (!all(finite_coords)) {
    stop("Coordinate CSV contains non-finite ", reduction, " coordinates.", call. = FALSE)
  }

  cluster_results <- lapply(seq_len(nrow(specs)), function(i) {
    cluster_best_dataset(
      plot_data = plot_data,
      dataset = specs$dataset[[i]],
      dataset_label = specs$dataset_label[[i]],
      cluster_prefix = specs$cluster_prefix[[i]],
      coord_names = coord_names,
      cluster_seed = cluster_seed + i - 1L,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      silhouette_sample_n = silhouette_sample_n
    )
  })
  names(cluster_results) <- specs$dataset
  clustered_best <- do.call(rbind, lapply(cluster_results, `[[`, "best"))
  silhouette_all <- do.call(rbind, lapply(cluster_results, `[[`, "silhouette"))
  cluster_summary <- do.call(rbind, lapply(cluster_results, `[[`, "summary"))

  full_marked <- plot_data
  full_marked$dataset_label <- NA_character_
  full_marked$cluster_scope <- NA_character_
  full_marked$cluster_source <- NA_character_
  full_marked$cluster_prefix <- NA_character_
  full_marked$cluster_base_id <- NA_character_
  full_marked$cluster_id <- NA_character_
  full_marked$cluster_num <- NA_integer_
  full_marked$cluster_k <- NA_integer_
  full_marked$cluster_silhouette_avg <- NA_real_
  full_marked$cluster_silhouette_sample_n <- NA_integer_
  cluster_cols <- c(
    "dataset_label", "cluster_scope", "cluster_source", "cluster_prefix", "cluster_base_id", "cluster_id", "cluster_num",
    "cluster_k", "cluster_silhouette_avg", "cluster_silhouette_sample_n"
  )
  for (cluster_result in cluster_results) {
    full_marked[cluster_result$best_idx, cluster_cols] <- cluster_result$best[, cluster_cols, drop = FALSE]
  }

  table_dir <- file.path(output_dir, "Tables")
  figure_dir <- file.path(output_dir, "Figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  cleanup_obsolete_outputs(output_dir, reduction_suffix)

  output_stem <- paste0("pooled_invivo_invitro_initial_vs_best_", reduction_suffix, "_best_clusters")
  full_table_path <- file.path(table_dir, paste0(output_stem, "_full_coordinates.csv"))
  best_table_path <- file.path(table_dir, paste0(output_stem, "_best_coordinates.csv"))
  silhouette_path <- file.path(table_dir, paste0(output_stem, "_silhouette.csv"))
  summary_path <- file.path(table_dir, paste0(output_stem, "_cluster_summary.csv"))
  seed_groups_path <- file.path(table_dir, paste0(output_stem, "_best_seed_groups.csv"))
  metadata_path <- file.path(table_dir, paste0(output_stem, "_metadata.csv"))
  subcluster_best_path <- file.path(table_dir, paste0(output_stem, "_best_subclusters.csv"))
  subcluster_seed_groups_path <- file.path(table_dir, paste0(output_stem, "_best_subcluster_seed_groups.csv"))
  subcluster_summary_path <- file.path(table_dir, paste0(output_stem, "_best_subcluster_summary.csv"))
  subcluster_silhouette_path <- file.path(table_dir, paste0(output_stem, "_best_subcluster_silhouette.csv"))
  subcluster_feature_metadata_path <- file.path(table_dir, paste0(output_stem, "_best_subcluster_feature_metadata.csv"))
  subcluster_zscore_features_path <- file.path(table_dir, paste0(output_stem, "_best_subcluster_zscore_features.csv"))
  seed_groups <- best_seed_group_table(clustered_best, reduction = reduction, coord_names = coord_names)

  write_csv(full_marked, full_table_path)
  write_csv(clustered_best, best_table_path)
  write_csv(silhouette_all, silhouette_path)
  write_csv(cluster_summary, summary_path)
  write_csv(seed_groups, seed_groups_path)
  write_csv(
    data.frame(
      key = c(
        "reduction", "input_pdf", "coordinate_csv", "output_dir", "n_full_rows",
        "n_invivo_best_rows", "n_invitro_best_rows", "coordinate_columns",
        "cluster_seed", "cluster_k_min", "cluster_k_max", "invivo_selected_k",
        "invitro_selected_k", "invivo_selected_average_silhouette",
        "invitro_selected_average_silhouette", "label_margin_frac",
        "subcluster_enabled", "subcluster_seed", "subcluster_k_min",
        "subcluster_k_max", "subcluster_min_n", "subcluster_label_margin_frac",
        "created_at"
      ),
      value = c(
        reduction, pdf_path, coordinate_csv, output_dir, nrow(plot_data),
        nrow(cluster_results$invivo$best), nrow(cluster_results$invitro$best),
        paste(coord_names, collapse = ","),
        cluster_seed, cluster_k_min, cluster_k_max,
        cluster_results$invivo$selected_summary$k[[1L]],
        cluster_results$invitro$selected_summary$k[[1L]],
        cluster_results$invivo$selected_summary$average_silhouette[[1L]],
        cluster_results$invitro$selected_summary$average_silhouette[[1L]],
        label_margin_frac,
        subcluster_enabled,
        subcluster_seed,
        subcluster_k_min,
        subcluster_k_max,
        subcluster_min_n,
        subcluster_label_margin_frac,
        format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      ),
      stringsAsFactors = FALSE
    ),
    metadata_path
  )

  plot_bounds <- add_plot_label_margin_rows(
    full_marked,
    coord_names = coord_names,
    pad_frac = label_margin_frac
  )
  overlay_plot <- add_cluster_outline_layers_external_labels(
    build_pooled_umap_plot(
      plot_bounds,
      initial_size = as_num(argv$initial_size, 0.22),
      best_size = as_num(argv$best_size, 1.25),
      initial_alpha = as_num(argv$initial_alpha, 0.28),
      coord_names = coord_names,
      axis_labels = axis_labels
    ),
    clustered_best,
    bounds_data = plot_bounds,
    coord_names = coord_names
  )
  save_plot_pair(
    overlay_plot,
    file.path(figure_dir, output_stem),
    width = 7.4,
    height = 6.5
  )

  subcluster_result <- NULL
  if (isTRUE(subcluster_enabled)) {
    feature_data <- load_pooled_best_feature_data(root_dir)
    subcluster_result <- analyze_best_subclusters(
      clustered_best = clustered_best,
      feature_data = feature_data,
      reduction = reduction,
      coord_names = coord_names,
      subcluster_seed = subcluster_seed,
      subcluster_k_min = subcluster_k_min,
      subcluster_k_max = subcluster_k_max,
      subcluster_min_n = subcluster_min_n,
      silhouette_sample_n = silhouette_sample_n
    )
    write_csv(subcluster_result$best, subcluster_best_path)
    write_csv(subcluster_result$seed_groups, subcluster_seed_groups_path)
    write_csv(subcluster_result$summary, subcluster_summary_path)
    write_csv(subcluster_result$silhouette, subcluster_silhouette_path)
    write_csv(subcluster_result$feature_metadata, subcluster_feature_metadata_path)
    write_csv(subcluster_result$zscore_features, subcluster_zscore_features_path)
    write_subcluster_plots(
      subclustered_best = subcluster_result$best,
      figure_dir = figure_dir,
      reduction_suffix = reduction_suffix,
      coord_names = coord_names,
      axis_labels = axis_labels,
      label_margin_frac = subcluster_label_margin_frac
    )
  }

  message("Selected k by dataset:")
  for (cluster_result in cluster_results) {
    message(
      "  ",
      cluster_result$dataset_label,
      " (",
      cluster_result$cluster_prefix,
      "): k=",
      cluster_result$selected_summary$k[[1L]],
      ", average silhouette=",
      signif(cluster_result$selected_summary$average_silhouette[[1L]], 4)
    )
  }
  message("Cluster summary:")
  print(cluster_summary, row.names = FALSE)
  if (!is.null(subcluster_result)) {
    message("Subcluster summary:")
    print(subcluster_result$summary, row.names = FALSE)
  }
  invisible(
    list(
      full_coordinates = full_table_path,
      best_coordinates = best_table_path,
      silhouette = silhouette_path,
      summary = summary_path,
      seed_groups = seed_groups_path,
      seed_group_data = seed_groups,
      subcluster_best = subcluster_best_path,
      subcluster_seed_groups = subcluster_seed_groups_path,
      subcluster_summary = subcluster_summary_path,
      subcluster_silhouette = subcluster_silhouette_path,
      subcluster_feature_metadata = subcluster_feature_metadata_path,
      subcluster_zscore_features = subcluster_zscore_features_path,
      subcluster_seed_group_data = if (is.null(subcluster_result)) NULL else subcluster_result$seed_groups,
      subcluster_summary_data = if (is.null(subcluster_result)) NULL else subcluster_result$summary,
      metadata = metadata_path
    )
  )
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  root_dir <- normalizePath(
    path.expand(argv$result_root %||% default_parameter_landscape_clustering_dir()),
    mustWork = FALSE
  )
  output_dir <- normalizePath(
    path.expand(
      argv$output_dir %||% file.path(root_dir, "pooled_invivo_invitro", "full_data_in_vivo_clustring")
    ),
    mustWork = FALSE
  )
  requested_reductions <- as_char_vec(argv$reductions %||% argv$reduction, c("tsne", "umap"))
  reductions <- unique(vapply(requested_reductions, normalize_reduction, character(1)))
  results <- lapply(reductions, analyze_embedding, argv = argv, root_dir = root_dir, output_dir = output_dir)
  names(results) <- reductions
  all_seed_groups <- do.call(rbind, lapply(results, `[[`, "seed_group_data"))
  if (!is.null(all_seed_groups) && nrow(all_seed_groups)) {
    write_csv(
      all_seed_groups,
      file.path(output_dir, "Tables", "pooled_invivo_invitro_best_seed_groups_by_method.csv")
    )
  }
  collect_result_data <- function(name) {
    pieces <- lapply(results, `[[`, name)
    pieces <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0L, pieces)
    if (!length(pieces)) return(NULL)
    rbind_fill_plain(pieces)
  }
  all_subcluster_seed_groups <- collect_result_data("subcluster_seed_group_data")
  if (!is.null(all_subcluster_seed_groups) && nrow(all_subcluster_seed_groups)) {
    write_csv(
      all_subcluster_seed_groups,
      file.path(output_dir, "Tables", "pooled_invivo_invitro_best_subclusters_by_method.csv")
    )
  }
  all_subcluster_summary <- collect_result_data("subcluster_summary_data")
  if (!is.null(all_subcluster_summary) && nrow(all_subcluster_summary)) {
    write_csv(
      all_subcluster_summary,
      file.path(output_dir, "Tables", "pooled_invivo_invitro_best_subcluster_summary_by_method.csv")
    )
  }
  invisible(results)
}

main()
