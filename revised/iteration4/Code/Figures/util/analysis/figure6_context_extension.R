#!/usr/bin/env Rscript

# Context-paired Figure 6 extensions. Model calculations use the explicit
# read-only O2_supply_demand_MAP dependency resolved by f6r_paths() and loaded
# by f6r_load_response_engine(); they never use a vendored iteration4 copy.

options(stringsAsFactors = FALSE, warn = 1)

f6x_model_source_files <- function(paths) {
  c(
    file.path(paths$oxygen_code, "simulation", "fix_o2_simulation.R"),
    file.path(paths$oxygen_code, "simulation", "fix_o2_simulation_shared_utils.R"),
    file.path(paths$oxygen_code, "model", "model_O2_supply_demand_MAP.R"),
    file.path(paths$oxygen_code, "model", "model_O2_supply_demand_MAP.cpp"),
    file.path(paths$oxygen_code, "util", "o2_supply_demand_map_shared.R"),
    file.path(paths$oxygen_code, "util", "o2_supply_demand_map_common_semantics.R")
  )
}

f6x_write_model_provenance <- function(paths) {
  source_files <- f6x_model_source_files(paths)
  f6r_require_files(source_files, "current-branch oxygen model source")
  provenance <- data.frame(
    model_context = "in vitro",
    repository = paths$repository_root,
    source_file = normalizePath(source_files, mustWork = TRUE),
    md5 = vapply(source_files, f6r_md5, character(1L)),
    user_reported_sync_reference_commit =
      "83953a874401e42cd176432786f889a896adc959",
    reproduction_identity = "absolute current-branch source path plus recorded MD5",
    access_mode = "read-only current-branch model implementation",
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(
    provenance,
    file.path(paths$figure6, "invitro_model_code_provenance.tsv")
  )
}

f6x_atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = "gzip")
  if (!file.copy(temporary, path, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Failed to publish cache: ", path)
  }
  unlink(temporary)
  normalizePath(path, mustWork = TRUE)
}

f6x_separate_invitro_endpoint_table <- function(paths) {
  source <- file.path(
    paths$root, "data", "Figures", "Figure3",
    "figure3g_parameter_endpoints_500seeds.tsv"
  )
  tab <- f6r_read_tsv(source)
  required <- c("seed", "seed_number", "objective", "parameter", "value")
  if (!all(required %in% names(tab))) {
    stop("Separate in-vitro endpoint table is missing required fields.")
  }
  tab <- tab[tab$parameter %in% f6r_shared_parameters(), required, drop = FALSE]
  if (nrow(tab) != 500L * length(f6r_shared_parameters()) ||
      length(unique(tab$seed_number)) != 500L || any(!is.finite(tab$value))) {
    stop("Expected 500 complete finite separate in-vitro parameter endpoints.")
  }
  list(data = tab, path = normalizePath(source, mustWork = TRUE))
}

f6x_separate_invitro_seed_cache <- function(
    seed_number, endpoints, cache_path, model_md5, rebuild = FALSE
) {
  profile <- "separate_invitro_current_model_201_o2_v1"
  if (file.exists(cache_path) && !isTRUE(rebuild)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$metadata$profile, profile) &&
        identical(cached$metadata$model_md5, model_md5) &&
        nrow(cached$curve) == 201L &&
        isTRUE(cached$qc$operator_qc_pass[[1L]])) {
      return(cached$qc)
    }
  }
  z <- endpoints[endpoints$seed_number == seed_number, , drop = FALSE]
  if (nrow(z) != length(f6r_shared_parameters()) ||
      !setequal(z$parameter, f6r_shared_parameters())) {
    stop("Incomplete separate in-vitro parameter endpoint for seed ", seed_number)
  }
  values <- stats::setNames(as.numeric(z$value), z$parameter)
  cfg_base <- list(o2_S0_upper_bound = 5)
  run_params <- prepare_run_params(
    values, simulation = "invitro", cfg = cfg_base, fixed_o2 = 5
  )
  config <- prepare_sim_cfg(
    cfg_base, argv = list(), fixed_o2 = 5, run_params = run_params
  )
  run_params$O2_growth <- isTRUE(config$O2_growth)
  run_params$ploidy_O2_death <- config$ploidy_O2_death
  oxygen <- seq(0, 5, length.out = 201L)
  rows <- lapply(oxygen, function(o2) {
    result <- fixo2_dominant_attractor_one(
      seed_id = paste0("seed", seed_number), run_params = run_params,
      model_env = globalenv(), cfg = config, O2 = o2
    )
    data.frame(
      seed_id = paste0("seed", seed_number),
      seed_number = seed_number,
      O2_pct = as.numeric(result$O2_pct[[1L]]),
      population_average_p_misseg =
        as.numeric(result$population_average_p_misseg[[1L]]),
      dominant_mean_ploidy =
        as.numeric(result$dominant_mean_ploidy[[1L]]),
      spectral_gap = as.numeric(result$spectral_gap[[1L]]),
      dominant_growth_rate = as.numeric(result$dominant_growth_rate[[1L]]),
      status = as.character(result$status[[1L]]),
      stringsAsFactors = FALSE
    )
  })
  curve <- do.call(rbind, rows)
  valid <- nrow(curve) == 201L && all(curve$status == "ok") &&
    all(is.finite(curve$population_average_p_misseg)) &&
    all(is.finite(curve$dominant_mean_ploidy)) &&
    all(is.finite(curve$spectral_gap)) &&
    all(is.finite(curve$dominant_growth_rate))
  qc <- data.frame(
    seed_number = seed_number, n_o2 = nrow(curve),
    minimum_o2 = min(curve$O2_pct), maximum_o2 = max(curve$O2_pct),
    all_status_ok = all(curve$status == "ok"),
    all_numeric_outputs_finite = valid,
    operator_qc_pass = valid, cache_path = cache_path,
    stringsAsFactors = FALSE
  )
  f6x_atomic_save_rds(
    list(
      metadata = list(
        profile = profile, seed_number = seed_number,
        objective = unique(z$objective), model_md5 = model_md5
      ),
      curve = curve, qc = qc
    ),
    cache_path
  )
  qc$cache_path <- normalizePath(cache_path, mustWork = TRUE)
  qc
}

f6x_compute_separate_invitro <- function(
    paths, n_core = 8L, rebuild = FALSE
) {
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
  )
  f6r_require_packages(c("Matrix", "Rcpp", "data.table"))
  f6r_load_response_engine(paths)
  provenance_path <- f6x_write_model_provenance(paths)
  model_md5 <- f6r_model_source_fingerprint(paths)
  endpoint_bundle <- f6x_separate_invitro_endpoint_table(paths)
  endpoints <- endpoint_bundle$data
  cache_root <- file.path(paths$figure6, "separate_invitro_seed_cache")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  seeds <- sort(unique(endpoints$seed_number))
  compute_one <- function(seed_number) {
    # future::multisession starts an independent R process, so compiled model
    # symbols must be initialized inside every worker rather than inherited
    # from the coordinator.
    f6r_load_response_engine(paths)
    tryCatch(
      f6x_separate_invitro_seed_cache(
        seed_number, endpoints,
        file.path(cache_root, paste0("seed", seed_number, ".rds")),
        model_md5 = model_md5, rebuild = rebuild
      ),
      error = function(e) structure(
        list(seed_number = seed_number, message = conditionMessage(e)),
        class = "f6x_error"
      )
    )
  }
  results <- f6r_resilient_lapply(seeds, compute_one, n_core = n_core)
  failed <- vapply(results, inherits, logical(1L), "f6x_error")
  if (any(failed)) {
    stop(
      "Separate in-vitro fixed-O2 failures: ",
      paste(vapply(results[failed], `[[`, character(1L), "message"), collapse = "; ")
    )
  }
  qc <- do.call(rbind, results)
  if (nrow(qc) != 500L || !all(qc$operator_qc_pass)) {
    stop("Separate in-vitro fixed-O2 endpoint QC failed.")
  }
  objects <- lapply(qc$cache_path, readRDS)
  curves <- do.call(rbind, lapply(objects, `[[`, "curve"))
  curves <- curves[order(curves$seed_number, curves$O2_pct), , drop = FALSE]
  if (nrow(curves) != 100500L || anyDuplicated(
    curves[, c("seed_number", "O2_pct")]
  )) {
    stop("Separate in-vitro fixed-O2 curve grid is incomplete or duplicated.")
  }
  classified <- response_classify_all_curves(curves)
  by_seed <- classified$summary
  smooth <- classified$smooth
  segments <- classified$segments
  objective_map <- unique(endpoints[, c("seed_number", "objective")])
  best_indices <- which(objective_map$objective == min(objective_map$objective))
  if (length(best_indices) != 1L) {
    stop("Separate in-vitro objective table lacks one unique global best seed.")
  }
  expected_best_seed <- objective_map$seed_number[[best_indices]]
  by_seed <- merge(by_seed, objective_map, by = "seed_number", all.x = TRUE)
  gap_rows <- lapply(split(curves$spectral_gap, curves$seed_number), function(gap) {
    data.frame(
      fraction_gap_below_0p005 = mean(gap < 0.005),
      fraction_gap_below_0p01 = mean(gap < 0.01),
      any_gap_below_0p005 = any(gap < 0.005),
      stringsAsFactors = FALSE
    )
  })
  gap_summary <- do.call(rbind, gap_rows)
  gap_summary$seed_number <- as.integer(rownames(gap_summary))
  rownames(gap_summary) <- NULL
  gap_summary$spectral_gap_class <- ifelse(
    gap_summary$fraction_gap_below_0p005 >= 0.25, "unreliable",
    ifelse(
      gap_summary$any_gap_below_0p005 |
        gap_summary$fraction_gap_below_0p01 >= 0.10,
      "caution", "reliable"
    )
  )
  by_seed <- merge(by_seed, gap_summary, by = "seed_number", all.x = TRUE)
  by_seed$context <- "in vitro"
  smooth$context <- "in vitro"
  if (nrow(segments)) segments$context <- "in vitro"
  class_counts <- data.frame(
    smooth_curve_class = response_curve_class_order,
    n_seed = as.integer(table(factor(
      by_seed$smooth_curve_class, levels = response_curve_class_order
    ))),
    stringsAsFactors = FALSE
  )
  class_counts$fraction_seed <- class_counts$n_seed / 500
  reliability_counts <- data.frame(
    spectral_gap_class = c("reliable", "caution", "unreliable"),
    n_seed = as.integer(table(factor(
      by_seed$spectral_gap_class,
      levels = c("reliable", "caution", "unreliable")
    ))),
    stringsAsFactors = FALSE
  )
  reliability_counts$fraction_seed <- reliability_counts$n_seed / 500
  outputs <- c(
    curves = f6r_write_tsv(
      curves, file.path(paths$figure6, "response_class_invitro_raw_curves.tsv")
    ),
    smooth = f6r_write_tsv(
      smooth, file.path(paths$figure6, "response_class_invitro_smoothed_curves.tsv")
    ),
    by_seed = f6r_write_tsv(
      by_seed, file.path(paths$figure6, "response_class_invitro_curve_class_by_seed.tsv")
    ),
    segments = f6r_write_tsv(
      segments, file.path(paths$figure6, "response_class_invitro_persistent_segments.tsv")
    ),
    counts = f6r_write_tsv(
      class_counts, file.path(paths$figure6, "response_class_invitro_class_counts.tsv")
    ),
    reliability = f6r_write_tsv(
      reliability_counts,
      file.path(paths$figure6, "response_class_invitro_reliability_counts.tsv")
    ),
    qc = f6r_write_tsv(
      qc, file.path(paths$figure6, "response_class_invitro_operator_qc.tsv")
    ),
    provenance = provenance_path
  )
  validation <- data.frame(
    check = c(
      "seed_count", "oxygen_count", "curve_row_count", "class_count_total",
      "reliability_count_total", "global_best_seed", "operator_qc_pass"
    ),
    observed = c(
      length(unique(curves$seed_number)), length(unique(curves$O2_pct)),
      nrow(curves), sum(class_counts$n_seed), sum(reliability_counts$n_seed),
      by_seed$seed_number[[which.min(by_seed$objective)]], all(qc$operator_qc_pass)
    ),
    expected = c(500, 201, 100500, 500, 500, expected_best_seed, TRUE),
    stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) ==
    as.character(validation$expected)
  outputs <- c(outputs, validation = f6r_write_tsv(
    validation,
    file.path(paths$figure6, "response_class_invitro_validation.tsv")
  ))
  if (!all(validation$passed)) {
    stop("Separate in-vitro response validation failed.")
  }
  invisible(list(
    curves = curves, smooth = smooth, by_seed = by_seed,
    segments = segments, class_counts = class_counts,
    reliability_counts = reliability_counts, qc = qc,
    validation = validation, paths = outputs
  ))
}

f6x_response_column <- function(
    smooth, by_seed, context_label, panel_tag, y_ranges, png_path, pdf_path
) {
  merged <- merge(
    smooth,
    by_seed[, c("seed_id", "smooth_curve_class"), drop = FALSE],
    by = "seed_id", all.x = TRUE, sort = FALSE
  )
  draw <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::par(
      mfrow = c(8, 1), mar = c(0.45, 3.1, 0.22, 0.35),
      oma = c(2.2, 2.3, 1.45, 0.1), family = response_font_family
    )
    for (i in seq_along(response_curve_class_order)) {
      class_name <- response_curve_class_order[[i]]
      z <- merged[merged$smooth_curve_class == class_name, , drop = FALSE]
      graphics::plot(
        c(0, 5), y_ranges[[class_name]], type = "n", xlab = "", ylab = "",
        xaxt = if (i == length(response_curve_class_order)) "s" else "n",
        cex.axis = 0.72, mgp = c(1.6, 0.38, 0), tcl = -0.2
      )
      graphics::grid(col = "#E4E4E4", lty = 3)
      if (nrow(z)) {
        seeds <- unique(z$seed_id)
        color <- grDevices::adjustcolor(
          response_curve_class_colors[[class_name]],
          alpha.f = max(0.18, min(0.55, 10 / sqrt(length(seeds))))
        )
        for (seed in seeds) {
          zz <- z[z$seed_id == seed, , drop = FALSE]
          zz <- zz[order(zz$O2_pct), , drop = FALSE]
          graphics::lines(
            zz$O2_pct, zz$smoothed_dominant_mean_ploidy,
            col = color, lwd = response_adaptive_lwd(length(seeds))
          )
        }
        median_curve <- stats::aggregate(
          smoothed_dominant_mean_ploidy ~ O2_pct, z,
          stats::median, na.rm = TRUE
        )
        graphics::lines(
          median_curve$O2_pct, median_curve$smoothed_dominant_mean_ploidy,
          col = "#111111", lwd = 1.0
        )
      }
      n_seed <- length(unique(z$seed_id))
      graphics::legend(
        "topright",
        legend = paste0(response_curve_class_labels[[class_name]], " (n=", n_seed, ")"),
        bty = "n", text.font = 2, cex = 0.76,
        inset = c(0.003, 0.005), x.intersp = 0.2, y.intersp = 0.7
      )
    }
    graphics::mtext("Fixed oxygen (%)", 1, outer = TRUE, line = 1.05, cex = 0.83)
    graphics::mtext(
      "Smoothed dominant mean ploidy", 2, outer = TRUE, line = 1.05, cex = 0.83
    )
    graphics::mtext(
      paste0(panel_tag, ". ", context_label, " oxygen-ploidy response classes"),
      3, outer = TRUE, at = 0.02, adj = 0, line = 0.12,
      font = 2, cex = 0.92
    )
  }
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(
    png_path, width = 4.5, height = 10.4, units = "in", res = 300,
    pointsize = 10, bg = "white"
  )
  draw(); grDevices::dev.off()
  grDevices::cairo_pdf(
    pdf_path, width = 4.5, height = 10.4, pointsize = 10,
    family = response_font_family, bg = "white"
  )
  draw(); grDevices::dev.off()
  c(png = normalizePath(png_path), pdf = normalizePath(pdf_path))
}

f6x_density_polygon <- function(values, class_index, side, context) {
  values <- values[is.finite(values)]
  if (length(values) < 2L || diff(range(values)) <= 0) return(NULL)
  den <- stats::density(values, n = 256L, adjust = 0.85)
  width <- 0.36 * den$y / max(den$y)
  sign <- if (identical(side, "left")) -1 else 1
  data.frame(
    x = c(class_index, class_index + sign * width, class_index),
    delta_objective = c(den$x[[1L]], den$x, den$x[[length(den$x)]]),
    context = context,
    polygon_group = paste(context, class_index, sep = "::"),
    stringsAsFactors = FALSE
  )
}

f6x_draw_supplement_6_1 <- function(
    workspace_root = f6r_find_workspace_root()
) {
  f6r_require_packages(c("ggplot2", "patchwork", "magick", "scales"))
  paths <- f6r_paths(workspace_root)
  f6r_load_response_engine(paths)
  panel_dir <- file.path(paths$figure6, "panels")
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  vivo_smooth <- f6r_read_tsv(file.path(
    paths$figure6, "response_class_smoothed_curves.tsv"
  ))
  vivo_class <- f6r_read_tsv(file.path(
    paths$figure6, "response_class_curve_class_by_seed.tsv"
  ))
  vitro_smooth <- f6r_read_tsv(file.path(
    paths$figure6, "response_class_invitro_smoothed_curves.tsv"
  ))
  vitro_class <- f6r_read_tsv(file.path(
    paths$figure6, "response_class_invitro_curve_class_by_seed.tsv"
  ))
  y_ranges <- stats::setNames(lapply(response_curve_class_order, function(cls) {
    values <- c(
      vivo_smooth$smoothed_dominant_mean_ploidy[
        vivo_smooth$seed_id %in% vivo_class$seed_id[vivo_class$smooth_curve_class == cls]
      ],
      vitro_smooth$smoothed_dominant_mean_ploidy[
        vitro_smooth$seed_id %in% vitro_class$seed_id[vitro_class$smooth_curve_class == cls]
      ]
    )
    response_panel_y_range(values, min_span = 0.25)
  }), response_curve_class_order)
  panel_a <- f6x_response_column(
    vivo_smooth, vivo_class, "In vivo", "A", y_ranges,
    file.path(panel_dir, "supp_fig6-1a_invivo_response_classes_8x1.png"),
    file.path(panel_dir, "supp_fig6-1a_invivo_response_classes_8x1.pdf")
  )
  panel_b <- f6x_response_column(
    vitro_smooth, vitro_class, "In vitro", "B", y_ranges,
    file.path(panel_dir, "supp_fig6-1b_invitro_response_classes_8x1.png"),
    file.path(panel_dir, "supp_fig6-1b_invitro_response_classes_8x1.pdf")
  )

  tsne_path <- file.path(
    paths$figure5,
    "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_full_coordinates.csv"
  )
  tsne <- utils::read.csv(tsne_path, check.names = FALSE)
  tsne <- tsne[tsne$point_type == "best" & tsne$dataset %in% c("invivo", "invitro"), ]
  tsne$seed_number <- as.integer(tsne$seed)
  tsne$context <- ifelse(tsne$dataset == "invivo", "in vivo", "in vitro")
  class_lookup <- rbind(
    data.frame(
      context = "in vivo", seed_number = vivo_class$seed_number,
      curve_class = vivo_class$smooth_curve_class, stringsAsFactors = FALSE
    ),
    data.frame(
      context = "in vitro", seed_number = vitro_class$seed_number,
      curve_class = vitro_class$smooth_curve_class, stringsAsFactors = FALSE
    )
  )
  tsne <- merge(tsne, class_lookup, by = c("context", "seed_number"), all.x = TRUE)
  if (nrow(tsne) != 1000L || anyNA(tsne$curve_class) ||
      any(table(tsne$context) != 500L)) {
    stop("Supplementary Figure 6-1C requires 500 classified endpoints per context.")
  }
  tsne$curve_class <- factor(tsne$curve_class, levels = response_curve_class_order)
  tsne$class_best <- FALSE
  for (ctx in c("in vivo", "in vitro")) {
    for (cls in response_curve_class_order) {
      idx <- which(tsne$context == ctx & as.character(tsne$curve_class) == cls)
      if (length(idx)) tsne$class_best[idx[[which.min(tsne$objective[idx])]]] <- TRUE
    }
  }
  vivo_regions <- tsne[tsne$context == "in vivo", ]
  vivo_regions$warm_start_region <- sub("^vi_", "", vivo_regions$cluster_id)
  hulls <- do.call(rbind, lapply(split(vivo_regions, vivo_regions$warm_start_region), function(d) {
    xy <- unique(d[, c("tSNE1", "tSNE2")])
    h <- grDevices::chull(xy$tSNE1, xy$tSNE2)
    out <- xy[c(h, h[[1L]]), , drop = FALSE]
    center <- colMeans(out)
    out$tSNE1 <- center[[1L]] + 1.035 * (out$tSNE1 - center[[1L]])
    out$tSNE2 <- center[[2L]] + 1.035 * (out$tSNE2 - center[[2L]])
    out$warm_start_region <- d$warm_start_region[[1L]]
    out
  }))
  p_c <- ggplot2::ggplot(
    tsne, ggplot2::aes(tSNE1, tSNE2, colour = curve_class, shape = context)
  ) +
    ggplot2::geom_point(size = 1.35, alpha = 0.72) +
    ggplot2::geom_point(
      data = tsne[tsne$class_best & tsne$context == "in vivo", ],
      shape = 1, size = 3.0, stroke = 0.75,
      colour = "#111111", show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = tsne[tsne$class_best & tsne$context == "in vitro", ],
      shape = 2, size = 3.2, stroke = 0.75,
      colour = "#111111", show.legend = FALSE
    ) +
    ggplot2::geom_path(
      data = hulls,
      ggplot2::aes(tSNE1, tSNE2, group = warm_start_region,
                   colour = warm_start_region),
      inherit.aes = FALSE, linewidth = 0.55, linetype = "22",
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(response_curve_class_colors, response_region_colors),
      breaks = response_curve_class_order,
      labels = response_curve_class_labels, drop = FALSE
    ) +
    ggplot2::scale_shape_manual(values = c("in vivo" = 16, "in vitro" = 17)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "C. Response classes in pooled parameter space",
      subtitle = "Circles: in vivo; triangles: in vitro; black outlines: class-best fits",
      x = "Pooled 14-parameter t-SNE coordinate 1",
      y = "Pooled 14-parameter t-SNE coordinate 2",
      colour = "O2-ploidy response class", shape = "Context"
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        order = 1, nrow = 1, byrow = TRUE,
        title.position = "top", title.hjust = 0.5
      ),
      colour = ggplot2::guide_legend(
        order = 2, ncol = 2, byrow = TRUE,
        title.position = "top", title.hjust = 0.5
      )
    ) +
    ggplot2::theme_classic(base_size = 9, base_family = response_font_family) +
    ggplot2::theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.text = ggplot2::element_text(size = 6.4),
      legend.title = ggplot2::element_text(size = 7.1),
      legend.key.width = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(1.2, "mm"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  objective_data <- rbind(
    data.frame(
      context = "in vivo", seed_number = vivo_class$seed_number,
      curve_class = vivo_class$smooth_curve_class,
      objective = tsne$objective[match(
        paste("in vivo", vivo_class$seed_number),
        paste(tsne$context, tsne$seed_number)
      )], stringsAsFactors = FALSE
    ),
    data.frame(
      context = "in vitro", seed_number = vitro_class$seed_number,
      curve_class = vitro_class$smooth_curve_class,
      objective = vitro_class$objective, stringsAsFactors = FALSE
    )
  )
  objective_data$delta_objective <- ave(
    objective_data$objective, objective_data$context,
    FUN = function(x) x - min(x, na.rm = TRUE)
  )
  objective_data$curve_class <- factor(
    objective_data$curve_class, levels = response_curve_class_order
  )
  polygon_rows <- list()
  for (i in seq_along(response_curve_class_order)) {
    cls <- response_curve_class_order[[i]]
    for (ctx in c("in vivo", "in vitro")) {
      z <- objective_data$delta_objective[
        objective_data$context == ctx & objective_data$curve_class == cls
      ]
      polygon_rows[[length(polygon_rows) + 1L]] <- f6x_density_polygon(
        z, i, if (ctx == "in vivo") "left" else "right", ctx
      )
    }
  }
  polygons <- do.call(rbind, polygon_rows[!vapply(polygon_rows, is.null, logical(1L))])
  p_d <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = polygons,
      ggplot2::aes(x, delta_objective, group = polygon_group, fill = context),
      alpha = 0.48, colour = "#333333", linewidth = 0.25
    ) +
    ggplot2::geom_boxplot(
      data = objective_data,
      ggplot2::aes(
        x = as.integer(curve_class) + ifelse(context == "in vivo", -0.10, 0.10),
        y = delta_objective, group = interaction(curve_class, context),
        colour = context
      ),
      width = 0.12, outlier.shape = NA, linewidth = 0.35
    ) +
    ggplot2::scale_fill_manual(values = c("in vivo" = "#0072B2", "in vitro" = "#CC79A7")) +
    ggplot2::scale_colour_manual(values = c("in vivo" = "#0072B2", "in vitro" = "#CC79A7")) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(response_curve_class_order),
      labels = c("Complex\nnonmono.", "Inverted\nU", "Monotone\ninc.", "U-shaped",
                 "Approx.\nflat", "Increase\nplateau", "Decrease\nplateau",
                 "Monotone\ndec."), limits = c(0.5, 8.5)
    ) +
    ggplot2::labs(
      title = "D. Full-MAP fit quality across response classes",
      subtitle = "Context-specific objective minus the minimum within that context",
      x = NULL, y = expression(Delta*" full-MAP objective"), fill = "Context",
      colour = "Context"
    ) +
    ggplot2::theme_classic(base_size = 9, base_family = response_font_family) +
    ggplot2::theme(
      legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(size = 6.8)
    )
  panel_c <- f6r_save_plot(
    p_c, file.path(panel_dir, "supp_fig6-1c_parameter_space_contexts"),
    width = 6.4, height = 5.2
  )
  panel_d <- f6r_save_plot(
    p_d, file.path(panel_dir, "supp_fig6-1d_objective_split_violin"),
    width = 6.4, height = 5.2
  )
  column_a <- magick::image_read(panel_a[["png"]])
  column_b <- magick::image_read(panel_b[["png"]])
  column_c <- magick::image_read(panel_c[["png"]])
  column_d <- magick::image_read(panel_d[["png"]])
  right_column <- magick::image_append(c(column_c, column_d), stack = TRUE)
  target_height <- max(
    magick::image_info(column_a)$height,
    magick::image_info(column_b)$height,
    magick::image_info(right_column)$height
  )
  pad_to_height <- function(image, height) {
    info <- magick::image_info(image)
    magick::image_extent(
      image, paste0(info$width, "x", height),
      gravity = "center", color = "white"
    )
  }
  column_a <- pad_to_height(column_a, target_height)
  column_b <- pad_to_height(column_b, target_height)
  right_column <- pad_to_height(right_column, target_height)
  assembled <- magick::image_append(
    c(column_a, column_b, right_column), stack = FALSE
  )
  rendered <- file.path(paths$figure6, "rendered")
  dir.create(rendered, recursive = TRUE, showWarnings = FALSE)
  output_png <- file.path(rendered, "supp_fig6-1_response_class_diagnostics.png")
  output_pdf <- file.path(rendered, "supp_fig6-1_response_class_diagnostics.pdf")
  magick::image_write(assembled, output_png, format = "png", density = 300)
  assembled_info <- magick::image_info(assembled)
  grDevices::cairo_pdf(
    output_pdf,
    width = assembled_info$width[[1L]] / 300,
    height = assembled_info$height[[1L]] / 300,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.raster(
    as.raster(assembled), x = 0.5, y = 0.5,
    width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
    interpolate = FALSE
  )
  grDevices::dev.off()
  published <- c(
    figure_png = f6r_publish(output_png, file.path(paths$figures, basename(output_png))),
    figure_pdf = f6r_publish(output_pdf, file.path(paths$figures, basename(output_pdf))),
    manuscript_png = f6r_publish(output_png, file.path(paths$manuscript_figures, basename(output_png))),
    manuscript_pdf = f6r_publish(output_pdf, file.path(paths$manuscript_figures, basename(output_pdf)))
  )
  validation <- data.frame(
    check = c(
      "invivo_seed_count", "invitro_seed_count", "tsne_context_counts",
      "panel_a_rows", "panel_b_rows", "figure_width_px", "figure_height_px",
      "pdf_nonempty", "publication_png_md5_match", "publication_pdf_md5_match"
    ),
    observed = c(
      nrow(vivo_class), nrow(vitro_class), paste(as.integer(table(tsne$context)), collapse = ","),
      8, 8, assembled_info$width[[1L]], assembled_info$height[[1L]],
      file.info(output_pdf)$size > 10000,
      f6r_md5(output_png) == f6r_md5(published[["manuscript_png"]]),
      f6r_md5(output_pdf) == f6r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(500, 500, "500,500", 8, 8, 4620, 3120, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) == as.character(validation$expected)
  validation_path <- f6r_write_tsv(
    validation, file.path(paths$figure6, "supp_fig6-1_validation.tsv")
  )
  if (!all(validation$passed)) stop("Supplementary Figure 6-1 validation failed.")
  f6x_refresh_output_manifest(file.path(
    paths$figure6, "supp_fig6-1_output_manifest.tsv"
  ))
  invisible(list(
    output = c(png = output_png, pdf = output_pdf), published = published,
    validation = validation, validation_path = validation_path,
    panels = c(panel_a, panel_b, panel_c, panel_d)
  ))
}

f6x_joint_context_endpoint_manifest <- function(
    paths, objective_bundle, cutoff = "q20", displayed_only = FALSE,
    output_name = "joint_invitro_endpoint_manifest.tsv"
) {
  acceptance <- objective_bundle$objectives
  eligible_field <- paste0("eligible_", cutoff)
  if (!eligible_field %in% names(acceptance)) stop("Unknown cutoff: ", cutoff)
  if (!"parameter_endpoint_group_invitro" %in% names(acceptance)) {
    stop("Joint acceptance table lacks in-vitro endpoint signatures.")
  }
  if (displayed_only) {
    display <- f6r_display_pair_manifest(acceptance$pair_id, "A")
    acceptance <- acceptance[acceptance$pair_id %in% display$pair_id, ]
  } else {
    display <- data.frame(
      display_label = sub("Sc.*", "", unique(acceptance$pair_label)),
      pair_label = unique(acceptance$pair_label),
      pair_id = unique(acceptance$pair_id), stringsAsFactors = FALSE
    )
  }
  selected <- acceptance[
    acceptance[[eligible_field]] & acceptance$hard_qc_pass, , drop = FALSE
  ]
  groups <- split(
    selected,
    interaction(
      selected$pair_id, selected$parameter_endpoint_group_invitro,
      drop = TRUE, lex.order = TRUE
    )
  )
  rows <- lapply(groups, function(z) {
    z <- z[order(z$objective_rank, z$seed_number), , drop = FALSE]
    data.frame(
      display_label = display$display_label[match(z$pair_id[[1L]], display$pair_id)],
      pair_label = z$pair_label[[1L]], pair_id = z$pair_id[[1L]],
      parameter_endpoint_group = z$parameter_endpoint_group_invitro[[1L]],
      representative_seed_number = z$seed_number[[1L]],
      representative_objective_rank = z$objective_rank[[1L]],
      representative_objective = z$objective[[1L]],
      endpoint_multiplicity = nrow(z),
      endpoint_multiplicity_q05 = sum(z$eligible_q05),
      endpoint_multiplicity_q10 = sum(z$eligible_q10),
      endpoint_multiplicity_q20 = sum(z$eligible_q20),
      represented_seed_numbers = paste(z$seed_number, collapse = ","),
      model_context = "in vitro", stringsAsFactors = FALSE
    )
  })
  endpoints <- do.call(rbind, rows)
  endpoints <- endpoints[order(
    match(endpoints$pair_id, display$pair_id),
    endpoints$representative_objective_rank,
    endpoints$representative_seed_number
  ), , drop = FALSE]
  rownames(endpoints) <- NULL
  expected_per_pair <- if (cutoff == "q20") 100L else if (cutoff == "q10") 50L else 25L
  represented <- tapply(endpoints$endpoint_multiplicity, endpoints$pair_id, sum)
  if (any(represented != expected_per_pair)) {
    stop("Unexpected represented in-vitro endpoint count for ", cutoff)
  }
  path <- f6r_write_tsv(endpoints, file.path(paths$figure6, output_name))
  invisible(list(endpoints = endpoints, display_manifest = display, path = path))
}

f6x_compute_joint_invitro_cache <- function(
    paths, objective_bundle, n_core = 8L, rebuild = FALSE
) {
  f6r_require_packages(c("Matrix", "Rcpp", "matrixStats", "data.table"))
  f6r_load_response_engine(paths)
  manifest <- f6x_joint_context_endpoint_manifest(
    paths, objective_bundle, cutoff = "q20", displayed_only = FALSE,
    output_name = "joint_invitro_q20_endpoint_manifest.tsv"
  )
  endpoints <- manifest$endpoints
  cache_root <- file.path(paths$figure6, "multiseed_endpoint_cache_invitro")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  cache_paths <- stats::setNames(
    file.path(
      cache_root, endpoints$pair_label,
      paste0("endpoint_", endpoints$parameter_endpoint_group, ".rds")
    ), endpoints$parameter_endpoint_group
  )
  contexts <- lapply(
    unique(endpoints$pair_id), f6r_pair_model_context,
    selected = objective_bundle$selected, paths = paths
  )
  names(contexts) <- unique(endpoints$pair_id)
  compute_one <- function(i) {
    f6r_load_response_engine(paths)
    z <- endpoints[i, , drop = FALSE]
    tryCatch(
      f6r_compute_seed_cache(
        pair_id = z$pair_id[[1L]],
        seed_number = z$representative_seed_number[[1L]],
        objective = z$representative_objective[[1L]],
        master = objective_bundle$master,
        context = contexts[[z$pair_id[[1L]]]],
        cache_path = cache_paths[[z$parameter_endpoint_group[[1L]]]],
        parameter_source = objective_bundle$master_path,
        full_surface = z$endpoint_multiplicity_q10[[1L]] > 0L,
        force_rebuild = rebuild,
        model_context = "in vitro",
        parameter_value_column = "vitro_natural",
        simulation_mode = "invitro",
        model_source_fingerprint = f6r_model_source_fingerprint(paths)
      ),
      error = function(e) structure(
        list(index = i, message = conditionMessage(e)), class = "f6x_error"
      )
    )
  }
  index <- seq_len(nrow(endpoints))
  result <- f6r_resilient_lapply(index, compute_one, n_core = n_core)
  failed <- vapply(result, inherits, logical(1L), "f6x_error")
  if (any(failed)) {
    stop(
      "Joint in-vitro endpoint failures: ",
      paste(vapply(result[failed], `[[`, character(1L), "message"), collapse = "; ")
    )
  }
  qc <- do.call(rbind, result)
  if (!all(qc$operator_qc_pass)) stop("Joint in-vitro endpoint operator QC failed.")
  qc_path <- f6r_write_tsv(
    qc, file.path(paths$figure6, "joint_multiseed_operator_qc_invitro.tsv")
  )
  invisible(list(
    manifest = manifest, endpoints = endpoints, cache_paths = cache_paths,
    qc = qc, qc_path = qc_path
  ))
}

f6x_seed_weighted_matrix <- function(objects, values, weights) {
  matrix <- do.call(cbind, lapply(objects, function(x) x$surface[[values]]))
  matrix[, rep(seq_along(objects), times = weights), drop = FALSE]
}

f6x_summarize_joint_invitro <- function(paths, objective_bundle, cache_bundle) {
  f6r_require_packages(c("matrixStats", "data.table"))
  acceptance <- objective_bundle$objectives
  endpoints <- cache_bundle$endpoints
  surface_rows <- list()
  trajectory_rows <- list()
  seed_claim_rows <- list()
  for (pair in unique(endpoints$pair_id)) {
    metadata <- endpoints[
      endpoints$pair_id == pair & endpoints$endpoint_multiplicity_q10 > 0,
      , drop = FALSE
    ]
    objects <- lapply(
      cache_bundle$cache_paths[metadata$parameter_endpoint_group], readRDS
    )
    reference_surface <- objects[[1L]]$surface
    reference_trajectory <- objects[[1L]]$trajectory
    weights <- metadata$endpoint_multiplicity_q10
    if (sum(weights) != 50L) stop("In-vitro q10 weights do not sum to 50: ", pair)
    expand <- rep(seq_along(objects), times = weights)
    ploidy <- do.call(cbind, lapply(objects, function(x) x$surface$dominant_mean_ploidy))
    gap <- do.call(cbind, lapply(objects, function(x) x$surface$spectral_gap))
    growth <- do.call(cbind, lapply(objects, function(x) x$surface$dominant_growth_rate))
    ploidy <- ploidy[, expand, drop = FALSE]
    gap <- gap[, expand, drop = FALSE]
    growth <- growth[, expand, drop = FALSE]
    q <- matrixStats::rowQuantiles(ploidy, probs = c(0.10, 0.25, 0.75, 0.90))
    fraction_high <- rowMeans(ploidy >= 2)
    surface_rows[[pair]] <- data.frame(
      pair_id = pair, pair_label = metadata$pair_label[[1L]],
      model_context = "in vitro", O2_pct = reference_surface$O2_pct,
      effective_p_misseg = reference_surface$effective_p_misseg,
      n_seed = 50L, n_unique_parameter_endpoint = nrow(metadata),
      dominant_mean_ploidy_median = matrixStats::rowMedians(ploidy),
      dominant_mean_ploidy_q10 = q[, 1L], dominant_mean_ploidy_q25 = q[, 2L],
      dominant_mean_ploidy_q75 = q[, 3L], dominant_mean_ploidy_q90 = q[, 4L],
      proportion_dominant_mean_ploidy_ge_2 = fraction_high,
      ploidy_regime_consensus = pmax(fraction_high, 1 - fraction_high),
      spectral_gap_median = matrixStats::rowMedians(gap),
      proportion_spectral_gap_below_0p005 = rowMeans(gap < 0.005),
      dominant_growth_rate_median = matrixStats::rowMedians(growth),
      max_abs_actual_minus_requested_p_misseg = max(vapply(
        objects, function(x) x$qc$max_abs_actual_minus_requested_p_misseg[[1L]],
        numeric(1L)
      )), cutoff = "q10", stringsAsFactors = FALSE
    )
    tr_ploidy <- do.call(cbind, lapply(objects, function(x) x$trajectory$dominant_mean_ploidy))[, expand, drop = FALSE]
    tr_p <- do.call(cbind, lapply(objects, function(x) x$trajectory$population_average_p_misseg))[, expand, drop = FALSE]
    tr_gap <- do.call(cbind, lapply(objects, function(x) x$trajectory$spectral_gap))[, expand, drop = FALSE]
    tq <- matrixStats::rowQuantiles(tr_ploidy, probs = c(0.10, 0.90))
    pq <- matrixStats::rowQuantiles(tr_p, probs = c(0.10, 0.90))
    trajectory_rows[[pair]] <- data.frame(
      pair_id = pair, pair_label = metadata$pair_label[[1L]],
      model_context = "in vitro", O2_pct = reference_trajectory$O2_pct,
      n_seed = 50L,
      population_average_p_misseg_median = matrixStats::rowMedians(tr_p),
      population_average_p_misseg_q10 = pq[, 1L],
      population_average_p_misseg_q90 = pq[, 2L],
      dominant_mean_ploidy_median = matrixStats::rowMedians(tr_ploidy),
      dominant_mean_ploidy_q10 = tq[, 1L], dominant_mean_ploidy_q90 = tq[, 2L],
      spectral_gap_median = matrixStats::rowMedians(tr_gap),
      proportion_spectral_gap_below_0p005 = rowMeans(tr_gap < 0.005),
      cutoff = "q10", stringsAsFactors = FALSE
    )
  }
  surface <- do.call(rbind, surface_rows)
  trajectory <- do.call(rbind, trajectory_rows)

  # Expand one claim profile per exact endpoint back to optimizer seeds.
  for (i in seq_len(nrow(endpoints))) {
    meta <- endpoints[i, , drop = FALSE]
    object <- readRDS(cache_bundle$cache_paths[[meta$parameter_endpoint_group]])
    claim <- f6r_seed_claim_row(object)
    seed_ids <- as.integer(strsplit(meta$represented_seed_numbers, ",", fixed = TRUE)[[1L]])
    z <- acceptance[
      acceptance$pair_id == meta$pair_id & acceptance$seed_number %in% seed_ids,
      , drop = FALSE
    ]
    expanded <- claim[rep(1L, nrow(z)), , drop = FALSE]
    expanded$seed_number <- z$seed_number
    expanded$objective <- z$objective
    expanded$parameter_endpoint_group <- z$parameter_endpoint_group_invitro
    expanded$model_context <- "in vitro"
    seed_claim_rows[[i]] <- expanded
  }
  seed_claims <- do.call(rbind, seed_claim_rows)
  seed_claims <- seed_claims[order(
    seed_claims$pair_label, seed_claims$objective, seed_claims$seed_number
  ), , drop = FALSE]
  expected_q20_claims <- 100L * f6r_family_count()
  if (nrow(seed_claims) != expected_q20_claims) {
    stop(
      "Expected ", expected_q20_claims,
      " q20 in-vitro optimizer-seed claim profiles."
    )
  }

  claim_specs <- list(
    both_o2_0 = c("Both ploidy regimes present at 0% O2", "surface_both_regimes_o2_0", "boolean"),
    both_o2_1 = c("Both ploidy regimes present at 1% O2", "surface_both_regimes_o2_1", "boolean"),
    both_o2_5 = c("Both ploidy regimes present at 5% O2", "surface_both_regimes_o2_5", "boolean"),
    trajectory_direction = c("Unmodified trajectory: O2 response direction", "trajectory_direction_o2_0_to_5", "direction"),
    surface_direction_low_cin = c("Surface O2 response at low effective missegregation", "surface_direction_o2_0_to_5_cin_low", "direction"),
    surface_direction_mid_cin = c("Surface O2 response at intermediate effective missegregation", "surface_direction_o2_0_to_5_cin_mid", "direction"),
    surface_direction_high_cin = c("Surface O2 response at high effective missegregation", "surface_direction_o2_0_to_5_cin_high", "direction"),
    spectral_gap_majority_reliable = c("At least half of the prespecified diagnostic union has spectral gap >= 0.005", "surface_fraction_spectral_gap_ge_0p005", "boolean_numeric")
  )
  summarize <- function(z, cutoff, weighting) {
    do.call(rbind, lapply(names(claim_specs), function(id) {
      spec <- claim_specs[[id]]; value <- z[[spec[[2L]]]]
      if (spec[[3L]] == "direction") {
        counts <- table(value); modal <- names(counts)[which.max(counts)]
        support <- mean(value == modal); estimand <- NA_real_
      } else {
        logical_value <- if (spec[[3L]] == "boolean_numeric") value >= 0.5 else as.logical(value)
        fraction <- mean(logical_value); modal <- if (fraction >= 0.5) "TRUE" else "FALSE"
        support <- max(fraction, 1 - fraction)
        estimand <- if (spec[[3L]] == "boolean_numeric") stats::median(value) else fraction
      }
      data.frame(
        cutoff = cutoff, pair_id = z$pair_id[[1L]], pair_label = z$pair_label[[1L]],
        model_context = "in vitro", claim_id = id, claim_label = spec[[1L]],
        analysis_weighting = weighting,
        n_seed = length(unique(z$seed_number)),
        n_unique_parameter_endpoint = length(unique(z$parameter_endpoint_group)),
        n_analysis_unit = nrow(z), modal_result = modal,
        modal_support_fraction = support, median_estimand = estimand,
        support_definition = paste0("modal state/direction: ", modal),
        stringsAsFactors = FALSE
      )
    }))
  }
  robust <- list(); unique_robust <- list()
  for (cutoff in c("q05", "q10", "q20")) {
    field <- paste0("eligible_", cutoff)
    eligible <- acceptance[acceptance[[field]], c("pair_id", "seed_number")]
    z_all <- merge(seed_claims, eligible, by = c("pair_id", "seed_number"))
    for (pair in unique(z_all$pair_id)) {
      z <- z_all[z_all$pair_id == pair, , drop = FALSE]
      z <- z[order(z$objective, z$seed_number), , drop = FALSE]
      robust[[length(robust) + 1L]] <- summarize(z, cutoff, "optimizer-seed endpoints")
      zu <- z[!duplicated(z$parameter_endpoint_group), , drop = FALSE]
      unique_robust[[length(unique_robust) + 1L]] <- summarize(zu, cutoff, "unique 14-parameter endpoints")
    }
  }
  robustness <- do.call(rbind, robust)
  robustness_unique <- do.call(rbind, unique_robust)
  dedup <- merge(
    robustness, robustness_unique,
    by = c("cutoff", "pair_id", "pair_label", "model_context", "claim_id", "claim_label"),
    suffixes = c("_seed_weighted", "_unique_endpoint"), all = TRUE
  )
  dedup$same_modal_result <- dedup$modal_result_seed_weighted == dedup$modal_result_unique_endpoint
  cutoff_consistency <- do.call(rbind, lapply(
    split(robustness, interaction(robustness$pair_id, robustness$claim_id)),
    function(z) data.frame(
      pair_id = z$pair_id[[1L]], pair_label = z$pair_label[[1L]],
      model_context = "in vitro", claim_id = z$claim_id[[1L]],
      claim_label = z$claim_label[[1L]], n_cutoff = nrow(z),
      same_modal_result_all_cutoffs = length(unique(z$modal_result)) == 1L,
      minimum_modal_support_across_cutoffs = min(z$modal_support_fraction),
      stringsAsFactors = FALSE
    )
  ))
  outputs <- c(
    surface = f6r_write_tsv(surface, file.path(paths$figure6, "joint_multiseed_surface_summary_invitro.tsv")),
    trajectory = f6r_write_tsv(trajectory, file.path(paths$figure6, "joint_multiseed_trajectory_summary_invitro.tsv")),
    seed_claims = f6r_write_tsv(seed_claims, file.path(paths$figure6, "joint_seed_biological_robustness_invitro.tsv")),
    robustness = f6r_write_tsv(robustness, file.path(paths$figure6, "joint_seed_claim_robustness_invitro.tsv")),
    robustness_unique = f6r_write_tsv(robustness_unique, file.path(paths$figure6, "joint_unique_parameter_claim_robustness_invitro.tsv")),
    dedup = f6r_write_tsv(dedup, file.path(paths$figure6, "joint_seed_vs_unique_parameter_robustness_invitro.tsv")),
    cutoff_consistency = f6r_write_tsv(cutoff_consistency, file.path(paths$figure6, "joint_seed_cutoff_consistency_invitro.tsv"))
  )
  invisible(list(
    surface = surface, trajectory = trajectory, seed_claims = seed_claims,
    robustness = robustness, robustness_unique = robustness_unique,
    dedup = dedup, cutoff_consistency = cutoff_consistency, paths = outputs
  ))
}

f6x_compute_dense_invitro <- function(
    paths, objective_bundle, n_core = 8L, rebuild = FALSE
) {
  f6r_require_packages(c("Matrix", "Rcpp", "matrixStats"))
  f6r_load_response_engine(paths)
  manifest <- f6x_joint_context_endpoint_manifest(
    paths, objective_bundle, cutoff = "q10", displayed_only = TRUE,
    output_name = "figure6_invitro_dense_endpoint_manifest.tsv"
  )
  endpoints <- manifest$endpoints
  cache_root <- file.path(paths$figure6, "figure6_invitro_dense_endpoint_cache")
  cache_paths <- stats::setNames(
    file.path(
      cache_root, endpoints$pair_label,
      paste0("endpoint_", endpoints$parameter_endpoint_group, ".rds")
    ), endpoints$parameter_endpoint_group
  )
  contexts <- lapply(
    manifest$display_manifest$pair_id, f6r_pair_model_context,
    selected = objective_bundle$selected, paths = paths
  )
  names(contexts) <- manifest$display_manifest$pair_id
  compute_one <- function(i) {
    f6r_load_response_engine(paths)
    z <- endpoints[i, , drop = FALSE]
    tryCatch(
      f6r_figure6d_compute_endpoint_cache(
        metadata = z, parameters = objective_bundle$parameters_invitro,
        context = contexts[[z$pair_id[[1L]]]],
        cache_path = cache_paths[[z$parameter_endpoint_group[[1L]]]],
        parameter_source = objective_bundle$paths[["parameters_invitro"]],
        force_rebuild = rebuild, model_context = "in vitro",
        simulation_mode = "invitro",
        model_source_fingerprint = f6r_model_source_fingerprint(paths)
      ),
      error = function(e) structure(
        list(index = i, message = conditionMessage(e)), class = "f6x_error"
      )
    )
  }
  result <- f6r_resilient_lapply(
    seq_len(nrow(endpoints)), compute_one, n_core = n_core
  )
  failed <- vapply(result, inherits, logical(1L), "f6x_error")
  if (any(failed)) stop(
    "Dense in-vitro endpoint failures: ",
    paste(vapply(result[failed], `[[`, character(1L), "message"), collapse = "; ")
  )
  qc <- do.call(rbind, result)
  if (!all(qc$operator_qc_pass)) stop("Dense in-vitro endpoint QC failed.")
  summary <- f6r_figure6d_summarize_dense_caches(paths, manifest, cache_paths)
  summary$model_context <- "in vitro"
  outputs <- c(
    summary = f6r_write_tsv(summary, file.path(paths$figure6, "figure6_invitro_fixed_p_curve_family.tsv")),
    qc = f6r_write_tsv(qc, file.path(paths$figure6, "figure6_invitro_dense_endpoint_qc.tsv")),
    manifest = manifest$path
  )
  validation <- data.frame(
    check = c("unique_endpoint_count", "represented_seed_count", "summary_rows", "operator_qc"),
    observed = c(nrow(qc), sum(qc$endpoint_multiplicity_q10), nrow(summary), all(qc$operator_qc_pass)),
    expected = c(
      nrow(qc), 50L * f6r_family_count(),
      f6r_family_count() * 496L * 201L, TRUE
    ), stringsAsFactors = FALSE
  )
  validation$passed <- c(
    nrow(qc) >= f6r_family_count(),
    sum(qc$endpoint_multiplicity_q10) == 50L * f6r_family_count(),
    nrow(summary) == f6r_family_count() * 496L * 201L,
    all(qc$operator_qc_pass)
  )
  outputs <- c(outputs, validation = f6r_write_tsv(
    validation, file.path(paths$figure6, "figure6_invitro_dense_validation.tsv")
  ))
  if (!all(validation$passed)) stop("Dense in-vitro validation failed.")
  inverse <- f6r_inverse_panel_data(
    paths, rebuild = rebuild, n_core = n_core,
    dense_qc_path = outputs[["qc"]], output_prefix = "figure6_invitro",
    model_context = "in vitro"
  )
  invisible(list(
    manifest = manifest, endpoints = endpoints, cache_paths = cache_paths,
    qc = qc, summary = summary, inverse = inverse, paths = outputs
  ))
}

f6x_data <- function(
    workspace_root = f6r_find_workspace_root(), n_core = 8L,
    rebuild = FALSE
) {
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
  )
  paths <- f6r_paths(workspace_root)
  f6r_load_response_engine(paths)
  f6x_write_model_provenance(paths)
  objective_bundle <- f6r_objective_selection(paths)
  multiseed_cache <- f6x_compute_joint_invitro_cache(
    paths, objective_bundle, n_core = n_core, rebuild = rebuild
  )
  multiseed <- f6x_summarize_joint_invitro(
    paths, objective_bundle, multiseed_cache
  )
  dense <- f6x_compute_dense_invitro(
    paths, objective_bundle, n_core = n_core, rebuild = rebuild
  )
  invisible(list(
    paths = paths, objective_bundle = objective_bundle,
    multiseed_cache = multiseed_cache, multiseed = multiseed, dense = dense
  ))
}

# Supplementary Figure 6-3 needs only the rank-1 population distribution and
# its local grid behavior. Reuse the exact q10 endpoint surfaces calculated for
# the main Figure 6 rather than diagonalizing the operator a second time.
f6x_si6_rank1_object <- function(object, parameter_endpoint_group) {
  required <- c("O2_pct", "effective_p_misseg", "dominant_mean_ploidy")
  if (!all(required %in% names(object$surface))) {
    stop("Joint in-vitro endpoint cache lacks rank-1 surface fields.")
  }
  list(
    metadata = list(
      parameter_endpoint_group = as.character(parameter_endpoint_group)
    ),
    grid = object$surface[, c("O2_pct", "effective_p_misseg"), drop = FALSE],
    localization_l1 = matrix(object$surface$dominant_mean_ploidy, ncol = 1L)
  )
}

f6x_si6_pair_summary <- function(grid, model_context) {
  pair_rows <- lapply(split(grid, grid$display_label), function(z) {
    weak <- z[z$weak_gap_region, , drop = FALSE]
    data.frame(
      display_label = z$display_label[[1L]], pair_label = z$pair_label[[1L]],
      n_grid_cell = nrow(z), n_weak_gap_cell = nrow(weak),
      fraction_grid_weak_gap = nrow(weak) / nrow(z),
      fraction_weak_gap_stable_high = mean(weak$regime_class == "Stable high"),
      fraction_weak_gap_stable_low = mean(weak$regime_class == "Stable low"),
      fraction_weak_gap_stable_intermediate = mean(
        weak$regime_class == "Stable intermediate"
      ),
      fraction_weak_gap_mixed = mean(weak$regime_class == "Mixed"),
      fraction_weak_gap_consensus_ge_0p9 = mean(
        weak$ploidy_regime_consensus >= 0.90
      ),
      fraction_weak_gap_any_endpoint_local_switch = mean(
        weak$local_regime_switch_proportion > 0
      ),
      fraction_weak_gap_majority_endpoint_local_switch = mean(
        weak$local_regime_switch_proportion >= 0.50
      ),
      weak_gap_ploidy_spread_median = stats::median(
        weak$dominant_ploidy_spread_q90_q10
      ),
      weak_gap_ploidy_spread_q90 = stats::quantile(
        weak$dominant_ploidy_spread_q90_q10, 0.90, names = FALSE
      ),
      weak_gap_ploidy_spread_max = max(weak$dominant_ploidy_spread_q90_q10),
      weak_gap_local_jump_median = stats::median(
        weak$local_adjacent_ploidy_jump_median
      ),
      weak_gap_local_jump_q90 = stats::quantile(
        weak$local_adjacent_ploidy_jump_median, 0.90, names = FALSE
      ),
      weak_gap_local_jump_max = max(weak$local_adjacent_ploidy_jump_median),
      fraction_weak_gap_local_jump_median_ge_1 = mean(
        weak$local_adjacent_ploidy_jump_median >= 1
      ),
      model_context = model_context, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pair_rows)
  out <- out[
    match(f6r_family_levels(), out$display_label), , drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

f6x_si6_context_from_q20 <- function(
    manifest, cache_paths, surface, model_context
) {
  f6r_require_files(unname(cache_paths), paste0(
    model_context, " q10 joint endpoint surfaces for Supplementary Figure 6-3"
  ))
  raw <- lapply(cache_paths, readRDS)
  qc <- do.call(rbind, lapply(raw, `[[`, "qc"))
  if (!all(qc$operator_qc_pass)) {
    stop("Supplementary Figure 6-3 encountered a non-passing q20 cache.")
  }
  names(raw) <- names(cache_paths)
  rows <- lapply(manifest$display_manifest$pair_id, function(pair) {
    metadata <- manifest$endpoints[
      manifest$endpoints$pair_id == pair, , drop = FALSE
    ]
    raw_objects <- raw[metadata$parameter_endpoint_group]
    objects <- Map(
      f6x_si6_rank1_object, raw_objects, metadata$parameter_endpoint_group
    )
    si6_summarize_weak_gap_pair(metadata, objects, surface)
  })
  grid <- do.call(rbind, rows)
  rownames(grid) <- NULL
  grid$model_context <- model_context
  list(
    grid = grid,
    pair_summary = f6x_si6_pair_summary(grid, model_context),
    qc = qc
  )
}

f6x_supplement_6_3_context_data <- function(
    workspace_root = f6r_find_workspace_root()
) {
  f6r_require_packages("matrixStats")
  paths <- f6r_paths(workspace_root)
  si_paths <- si6_paths(workspace_root)
  dir.create(si_paths$data, recursive = TRUE, showWarnings = FALSE)
  objective_bundle <- f6r_objective_selection(paths)

  vivo_manifest <- f6r_figure6d_endpoint_manifest(paths, objective_bundle)
  vivo_cache_paths <- stats::setNames(
    file.path(
      paths$figure6, "multiseed_seed_cache",
      vivo_manifest$endpoints$pair_label,
      paste0("seed", vivo_manifest$endpoints$representative_seed_number, ".rds")
    ),
    vivo_manifest$endpoints$parameter_endpoint_group
  )
  vivo_surface <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_surface_summary.tsv"
  ))
  vivo <- f6x_si6_context_from_q20(
    vivo_manifest, vivo_cache_paths, vivo_surface, "in vivo"
  )

  vitro_manifest <- f6x_joint_context_endpoint_manifest(
    paths, objective_bundle, cutoff = "q10", displayed_only = TRUE,
    output_name = "supp_figure6-3_invitro_endpoint_manifest.tsv"
  )
  vitro_cache_paths <- stats::setNames(
    file.path(
      paths$figure6, "multiseed_endpoint_cache_invitro",
      vitro_manifest$endpoints$pair_label,
      paste0(
        "endpoint_", vitro_manifest$endpoints$parameter_endpoint_group, ".rds"
      )
    ),
    vitro_manifest$endpoints$parameter_endpoint_group
  )
  vitro_surface <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_surface_summary_invitro.tsv"
  ))
  vitro <- f6x_si6_context_from_q20(
    vitro_manifest, vitro_cache_paths, vitro_surface, "in vitro"
  )

  combined_grid <- rbind(vivo$grid, vitro$grid[, names(vivo$grid), drop = FALSE])
  combined_pair <- rbind(
    vivo$pair_summary,
    vitro$pair_summary[, names(vivo$pair_summary), drop = FALSE]
  )
  combined_grid$model_context <- factor(
    combined_grid$model_context, levels = c("in vivo", "in vitro")
  )
  combined_pair$model_context <- factor(
    combined_pair$model_context, levels = c("in vivo", "in vitro")
  )
  grid_path <- f6r_write_tsv(combined_grid, file.path(
    si_paths$data, "supp_figure6-3_weak_gap_regime_robustness.tsv"
  ))
  pair_path <- f6r_write_tsv(combined_pair, file.path(
    si_paths$data, "supp_figure6-3_weak_gap_pair_summary.tsv"
  ))

  endpoint_columns <- c(
    "display_label", "pair_label", "pair_id", "parameter_endpoint_group",
    "representative_seed_number", "representative_objective_rank",
    "representative_objective", "endpoint_multiplicity_q10",
    "represented_seed_numbers"
  )
  endpoint_manifest <- rbind(
    transform(
      vivo_manifest$endpoints[, endpoint_columns, drop = FALSE],
      model_context = "in vivo"
    ),
    transform(
      vitro_manifest$endpoints[, endpoint_columns, drop = FALSE],
      model_context = "in vitro"
    )
  )
  endpoint_manifest_path <- f6r_write_tsv(endpoint_manifest, file.path(
    si_paths$data, "supp_figure6-3_endpoint_manifest.tsv"
  ))

  key <- function(x) paste(
    x$model_context, x$pair_id, sprintf("%.12f", x$O2_pct),
    sprintf("%.12g", x$effective_p_misseg), sep = "|"
  )
  reference_columns <- c(
    "pair_id", "O2_pct", "effective_p_misseg",
    "dominant_mean_ploidy_median", "dominant_mean_ploidy_q10",
    "dominant_mean_ploidy_q90", "proportion_dominant_mean_ploidy_ge_2"
  )
  reference <- rbind(
    transform(
      vivo_surface[
        vivo_surface$cutoff == "q10", reference_columns, drop = FALSE
      ],
      model_context = "in vivo"
    ),
    transform(
      vitro_surface[
        vitro_surface$cutoff == "q10", reference_columns, drop = FALSE
      ],
      model_context = "in vitro"
    )
  )
  matched <- match(key(combined_grid), key(reference))
  max_difference <- function(field) {
    if (anyNA(matched)) return(Inf)
    max(abs(combined_grid[[field]] - reference[[field]][matched]))
  }
  surface_differences <- c(
    median = max_difference("dominant_mean_ploidy_median"),
    q10 = max_difference("dominant_mean_ploidy_q10"),
    q90 = max_difference("dominant_mean_ploidy_q90"),
    ge2 = max_difference("proportion_dominant_mean_ploidy_ge_2")
  )
  expected_rows <- 2L * f6r_family_count() * 201L * 60L
  validation <- data.frame(
    check = c(
      "model_context_count", "displayed_context_pair_count", "grid_row_count",
      "pair_summary_row_count", "fifty_seed_weight_per_context_pair",
      "unique_endpoint_cache_count_per_context", "all_q20_cache_qc_pass",
      "regime_proportions_sum_to_one", "weak_gap_cells_present",
      "surface_unmatched_rows", "surface_median_max_abs_difference",
      "surface_q10_max_abs_difference", "surface_q90_max_abs_difference",
      "surface_ge2_share_max_abs_difference"
    ),
    observed = c(
      length(unique(combined_grid$model_context)),
      nrow(unique(combined_grid[, c("model_context", "display_label")])),
      nrow(combined_grid), nrow(combined_pair),
      paste(sort(unique(combined_grid$n_seed)), collapse = ","),
      all(c(
        nrow(vivo_manifest$endpoints), nrow(vitro_manifest$endpoints)
      ) > 0L),
      all(vivo$qc$operator_qc_pass) && all(vitro$qc$operator_qc_pass),
      max(abs(rowSums(combined_grid[, c(
        "proportion_low_ploidy_le_2",
        "proportion_intermediate_ploidy_gt_2_lt_4",
        "proportion_high_ploidy_ge_4"
      )]) - 1)),
      all(tapply(
        combined_grid$weak_gap_region,
        interaction(combined_grid$model_context, combined_grid$display_label),
        sum
      ) > 0),
      sum(is.na(matched)), unname(surface_differences)
    ),
    expected = c(
      2, 2L * f6r_family_count(), expected_rows,
      2L * f6r_family_count(), "50", TRUE, TRUE, "<=1e-12", TRUE,
      0, "<=1e-10", "<=1e-10", "<=1e-10", "<=1e-12"
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- c(
    validation$observed[[1L]] == 2L,
    validation$observed[[2L]] == 2L * f6r_family_count(),
    validation$observed[[3L]] == expected_rows,
    validation$observed[[4L]] == 2L * f6r_family_count(),
    validation$observed[[5L]] == "50",
    identical(as.character(validation$observed[[6L]]), "TRUE"),
    identical(as.character(validation$observed[[7L]]), "TRUE"),
    as.numeric(validation$observed[[8L]]) <= 1e-12,
    identical(as.character(validation$observed[[9L]]), "TRUE"),
    as.numeric(validation$observed[[10L]]) == 0,
    as.numeric(validation$observed[[11L]]) <= 1e-10,
    as.numeric(validation$observed[[12L]]) <= 1e-10,
    as.numeric(validation$observed[[13L]]) <= 1e-10,
    as.numeric(validation$observed[[14L]]) <= 1e-12
  )
  context_validation_path <- f6r_write_tsv(
    validation, file.path(si_paths$data, "supp_figure6-3_context_validation.tsv")
  )
  data_validation_path <- f6r_write_tsv(
    validation, file.path(si_paths$data, "supp_figure6-3_data_validation.tsv")
  )
  if (!all(validation$passed)) {
    stop(
      "Context-paired Supplementary Figure 6-3 validation failed: ",
      paste(validation$check[!validation$passed], collapse = ", ")
    )
  }
  invisible(list(
    grid = grid_path, pair_summary = pair_path,
    validation = context_validation_path,
    data_validation = data_validation_path,
    endpoint_manifest = endpoint_manifest_path
  ))
}

f6x_add_display_fields <- function(data, display_manifest, context) {
  data <- data[data$pair_id %in% display_manifest$pair_id, , drop = FALSE]
  data$display_label <- display_manifest$display_label[
    match(data$pair_id, display_manifest$pair_id)
  ]
  data$model_context <- context
  data$display_label <- factor(
    data$display_label, levels = f6r_family_levels()
  )
  data$model_context <- factor(
    data$model_context, levels = c("in vivo", "in vitro")
  )
  data
}

f6x_rbind_fill <- function(...) {
  pieces <- list(...)
  columns <- unique(unlist(lapply(pieces, names), use.names = FALSE))
  pieces <- lapply(pieces, function(x) {
    missing <- setdiff(columns, names(x))
    for (column in missing) x[[column]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, pieces)
}

f6x_refresh_output_manifest <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  manifest <- f6r_read_tsv(path)
  path_column <- intersect(c("path", "output_file"), names(manifest))
  if (!length(path_column) || !"md5" %in% names(manifest)) {
    stop("Unsupported output-manifest schema: ", path)
  }
  path_column <- path_column[[1L]]
  manifest$exists <- file.exists(manifest[[path_column]])
  manifest$size_bytes <- file.info(manifest[[path_column]])$size
  manifest$md5 <- vapply(manifest[[path_column]], f6r_md5, character(1L))
  f6r_write_tsv(manifest, path)
}

f6x_main_surface_plot <- function(paths) {
  f6r_require_packages(c("ggplot2", "isoband", "scales"))
  vivo <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_surface_summary.tsv"
  ))
  vitro <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_surface_summary_invitro.tsv"
  ))
  vivo_tr <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_trajectory_summary.tsv"
  ))
  vitro_tr <- f6r_read_tsv(file.path(
    paths$figure6, "joint_multiseed_trajectory_summary_invitro.tsv"
  ))
  vivo <- vivo[vivo$cutoff == "q10", ]; vitro <- vitro[vitro$cutoff == "q10", ]
  vivo_tr <- vivo_tr[vivo_tr$cutoff == "q10", ]; vitro_tr <- vitro_tr[vitro_tr$cutoff == "q10", ]
  display <- f6r_display_pair_manifest(c(vivo$pair_id, vitro$pair_id), "A")
  surface <- f6x_rbind_fill(
    f6x_add_display_fields(vivo, display, "in vivo"),
    f6x_add_display_fields(vitro, display, "in vitro")
  )
  trajectory <- f6x_rbind_fill(
    f6x_add_display_fields(vivo_tr, display, "in vivo"),
    f6x_add_display_fields(vitro_tr, display, "in vitro")
  )
  surface$log10_effective_p_misseg <- log10(surface$effective_p_misseg)
  trajectory$log10_p_median <- log10(trajectory$population_average_p_misseg_median)
  trajectory$log10_p_q10 <- log10(trajectory$population_average_p_misseg_q10)
  trajectory$log10_p_q90 <- log10(trajectory$population_average_p_misseg_q90)
  hatch_rows <- list()
  for (ctx in levels(surface$model_context)) for (label in levels(surface$display_label)) {
    z <- surface[surface$model_context == ctx & surface$display_label == label, ]
    if (!nrow(z)) next
    h <- f6r_weak_gap_hatch_data(z)
    if (!nrow(h)) next
    h$model_context <- ctx; h$display_label <- label
    hatch_rows[[length(hatch_rows) + 1L]] <- h
  }
  hatch <- if (length(hatch_rows)) do.call(rbind, hatch_rows) else data.frame(
    O2_pct = numeric(), log10_effective_p_misseg = numeric(),
    hatch_group = integer(), model_context = character(),
    display_label = character(), stringsAsFactors = FALSE
  )
  hatch$model_context <- factor(hatch$model_context, levels = c("in vivo", "in vitro"))
  hatch$display_label <- factor(
    hatch$display_label, levels = f6r_family_levels()
  )
  low <- surface[
    surface$ploidy_regime_consensus < 0.80 &
      match(surface$O2_pct, sort(unique(surface$O2_pct))) %% 8L == 1L &
      match(surface$effective_p_misseg, sort(unique(surface$effective_p_misseg))) %% 3L == 1L,
  ]
  fill_limits <- range(
    c(surface$dominant_mean_ploidy_median, 1, 7), na.rm = TRUE
  )
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = surface,
      ggplot2::aes(O2_pct, log10_effective_p_misseg,
                   fill = dominant_mean_ploidy_median)
    ) +
    ggplot2::geom_path(
      data = hatch,
      ggplot2::aes(O2_pct, log10_effective_p_misseg,
                   group = interaction(model_context, display_label, hatch_group)),
      colour = "#9B59B6", linewidth = 0.17, alpha = 0.70
    ) +
    ggplot2::geom_contour(
      data = surface,
      ggplot2::aes(O2_pct, log10_effective_p_misseg,
                   z = proportion_spectral_gap_below_0p005),
      breaks = 0.5, colour = "#9B59B6", linetype = "dotted",
      linewidth = 0.30
    ) +
    ggplot2::geom_point(
      data = low, ggplot2::aes(O2_pct, log10_effective_p_misseg),
      shape = 4, size = 0.32, stroke = 0.22, colour = "white"
    ) +
    ggplot2::geom_ribbon(
      data = trajectory,
      ggplot2::aes(O2_pct, ymin = log10_p_q10, ymax = log10_p_q90),
      inherit.aes = FALSE, fill = "#E0E0E0", alpha = 0.62
    ) +
    ggplot2::geom_path(
      data = trajectory,
      ggplot2::aes(O2_pct, log10_p_median),
      inherit.aes = FALSE, colour = "#111111", linewidth = 0.48
    ) +
    ggplot2::facet_grid(model_context ~ display_label, switch = "y") +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "#FFFFBF", "#B2182B"),
      trans = "log10", limits = fill_limits,
      breaks = c(1, 1.5, 2, 3, 4, 6),
      name = "Median dominant\nmean ploidy\n(log colors)"
    ) +
    ggplot2::scale_y_continuous(
      breaks = log10(c(0.005, 0.01, 0.05, 0.1, 0.5)),
      labels = c("0.005", "0.01", "0.05", "0.10", "0.50")
    ) +
    ggplot2::coord_cartesian(ylim = log10(c(0.005, 0.5)), expand = FALSE) +
    ggplot2::labs(
      title = "A. Oxygen-CIN-ploidy response surfaces",
      subtitle = paste0(
        "Tiles: seed-weighted q10 median; black line and gray band: fitted median and 10-90% trajectory; ",
        "purple hatching: weak spectral gap"
      ),
      x = "Fixed oxygen (%)", y = expression("Effective "*p[miss]*" probability")
    ) +
    ggplot2::theme_classic(base_size = 8.3, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 7.2, colour = "#555555"),
      strip.text = ggplot2::element_text(face = "bold", size = 8.5),
      strip.placement = "outside", panel.spacing = grid::unit(2.2, "mm"),
      legend.position = "right", legend.title = ggplot2::element_text(size = 7.2),
      legend.text = ggplot2::element_text(size = 6.8)
    )
}

f6x_main_inverse_plot <- function(paths) {
  f6r_require_packages(c("ggplot2", "isoband", "scales"))
  vivo <- f6r_read_tsv(file.path(paths$figure6, "figure6_inverse_response_summary.tsv"))
  vitro <- f6r_read_tsv(file.path(paths$figure6, "figure6_invitro_inverse_response_summary.tsv"))
  vivo_dense <- f6r_read_tsv(file.path(paths$figure6, "figure6d_fixed_p_curve_family.tsv"))
  vitro_dense <- f6r_read_tsv(file.path(paths$figure6, "figure6_invitro_fixed_p_curve_family.tsv"))
  display <- f6r_display_pair_manifest(c(vivo$pair_id, vitro$pair_id), "B")
  inverse <- f6x_rbind_fill(
    f6x_add_display_fields(vivo, display, "in vivo"),
    f6x_add_display_fields(vitro, display, "in vitro")
  )
  dense <- f6x_rbind_fill(
    f6x_add_display_fields(vivo_dense, display, "in vivo"),
    f6x_add_display_fields(vitro_dense, display, "in vitro")
  )
  highlighted_values <- c(0.01, 0.10, 0.20, 0.30)
  highlighted <- dense[vapply(
    dense$effective_p_misseg,
    function(x) any(abs(x - highlighted_values) < 1e-12), logical(1L)
  ), ]
  highlighted$reference_label <- factor(
    sprintf("%.2f", highlighted$effective_p_misseg),
    levels = c("0.01", "0.10", "0.20", "0.30")
  )
  mean_curve <- stats::aggregate(
    dominant_mean_ploidy_median ~ model_context + display_label + pair_id + O2_pct,
    dense, mean
  )
  mean_curve$reference_label <- factor(
    "Mean across 496 fixed p_miss,eff values",
    levels = c("0.01", "0.10", "0.20", "0.30",
               "Mean across 496 fixed p_miss,eff values")
  )
  hatch_rows <- list()
  for (ctx in levels(inverse$model_context)) for (label in levels(inverse$display_label)) {
    z <- inverse[inverse$model_context == ctx & inverse$display_label == label, ]
    if (!nrow(z)) next
    h <- f6r_inverse_multivalue_hatch_data(z)
    if (!nrow(h)) next
    h$model_context <- ctx; h$display_label <- label
    hatch_rows[[length(hatch_rows) + 1L]] <- h
  }
  hatch <- if (length(hatch_rows)) do.call(rbind, hatch_rows) else data.frame(
    O2_pct = numeric(), target_ploidy = numeric(), hatch_group = integer(),
    model_context = character(), display_label = character(),
    stringsAsFactors = FALSE
  )
  hatch$model_context <- factor(hatch$model_context, levels = c("in vivo", "in vitro"))
  hatch$display_label <- factor(
    hatch$display_label, levels = f6r_family_levels()
  )
  reference_levels <- c(
    "0.01", "0.10", "0.20", "0.30",
    "Mean across 496 fixed p_miss,eff values"
  )
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = inverse, ggplot2::aes(O2_pct, target_ploidy, fill = p_display)
    ) +
    ggplot2::geom_hline(
      yintercept = c(2, 4), colour = "#666666", linewidth = 0.22,
      linetype = "longdash"
    ) +
    ggplot2::geom_path(
      data = hatch,
      ggplot2::aes(O2_pct, target_ploidy,
                   group = interaction(model_context, display_label, hatch_group)),
      colour = "#7B3294", linewidth = 0.17, alpha = 0.72
    ) +
    ggplot2::geom_contour(
      data = inverse,
      ggplot2::aes(O2_pct, target_ploidy, z = fraction_multiple_solutions),
      breaks = 0.20, colour = "#7B3294", linetype = "dotted",
      linewidth = 0.30
    ) +
    ggplot2::geom_path(
      data = highlighted,
      ggplot2::aes(O2_pct, dominant_mean_ploidy_median,
                   group = reference_label, linetype = reference_label),
      colour = "#111111", linewidth = 0.50
    ) +
    ggplot2::geom_path(
      data = mean_curve,
      ggplot2::aes(O2_pct, dominant_mean_ploidy_median,
                   group = reference_label, linetype = reference_label),
      colour = "#D62728", linewidth = 0.68
    ) +
    ggplot2::facet_grid(model_context ~ display_label, switch = "y") +
    ggplot2::scale_fill_viridis_c(
      option = "D", trans = "log10", limits = c(0.005, 0.5),
      breaks = c(0.005, 0.01, 0.05, 0.10, 0.50), na.value = "#EFEFEF",
      name = "Median required\np_miss,eff\n(log colors)"
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "0.01" = "solid", "0.10" = "F28282",
        "0.20" = "dotdash", "0.30" = "dotted",
        "Mean across 496 fixed p_miss,eff values" = "solid"
      ), breaks = reference_levels,
      labels = c(
        "p_miss,eff = 0.01", "p_miss,eff = 0.10",
        "p_miss,eff = 0.20", "p_miss,eff = 0.30",
        "Mean across 496 fixed p_miss,eff values"
      ), name = "Reference curves"
    ) +
    ggplot2::labs(
      title = "B. Effective missegregation required for target ploidy",
      subtitle = "Gray: no stable unique inverse; purple hatching: multiple inverse solutions",
      x = "Fixed oxygen (%)", y = "Target dominant mean ploidy"
    ) +
    ggplot2::coord_cartesian(ylim = c(1, 7)) +
    ggplot2::theme_classic(base_size = 8.3, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 7.2, colour = "#555555"),
      strip.text = ggplot2::element_text(face = "bold", size = 8.5),
      strip.placement = "outside", panel.spacing = grid::unit(2.2, "mm"),
      legend.position = "right", legend.title = ggplot2::element_text(size = 7.2),
      legend.text = ggplot2::element_text(size = 6.6)
    )
}

f6x_draw_main <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c("ggplot2", "patchwork", "magick", "isoband"))
  paths <- f6r_paths(workspace_root)
  panel_dir <- file.path(paths$figure6, "panels")
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  p_a <- f6x_main_surface_plot(paths)
  p_b <- f6x_main_inverse_plot(paths)
  panel_a <- f6r_save_plot(
    p_a, file.path(panel_dir, "figure6a_invivo_invitro_response_surfaces"),
    width = 18.4, height = 4.6
  )
  panel_b <- f6r_save_plot(
    p_b, file.path(panel_dir, "figure6b_invivo_invitro_inverse_response"),
    width = 18.4, height = 4.6
  )
  combined <- p_a / p_b + patchwork::plot_layout(heights = c(1, 1))
  output <- f6r_save_plot(
    combined, file.path(paths$figure6, "rendered", "assembled_fig6"),
    width = 18.8, height = 9.1
  )
  published <- c(
    figures_png = f6r_publish(output[["png"]], file.path(paths$figures, "assembled_fig6.png")),
    figures_pdf = f6r_publish(output[["pdf"]], file.path(paths$figures, "assembled_fig6.pdf")),
    manuscript_png = f6r_publish(output[["png"]], file.path(paths$manuscript_figures, "assembled_fig6.png")),
    manuscript_pdf = f6r_publish(output[["pdf"]], file.path(paths$manuscript_figures, "assembled_fig6.pdf"))
  )
  info <- magick::image_info(magick::image_read(output[["png"]]))
  validation <- data.frame(
    check = c(
      "panel_count", "context_count", "displayed_pair_count", "figure_width_px",
      "publication_png_md5_match", "publication_pdf_md5_match"
    ),
    observed = c(
      2, 2, 6, info$width[[1L]],
      f6r_md5(output[["png"]]) == f6r_md5(published[["manuscript_png"]]),
      f6r_md5(output[["pdf"]]) == f6r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(2, 2, 6, 5640, TRUE, TRUE), stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) == as.character(validation$expected)
  validation_path <- f6r_write_tsv(
    validation, file.path(paths$figure6, "figure6_context_validation.tsv")
  )
  if (!all(validation$passed)) stop("Context-paired Figure 6 validation failed.")
  f6x_refresh_output_manifest(file.path(
    paths$figure6, "figure6_output_manifest.tsv"
  ))
  invisible(list(
    output = output, published = published, validation = validation,
    validation_path = validation_path, panels = c(panel_a, panel_b)
  ))
}

f6x_draw_supplement_6_2 <- function(
    workspace_root = f6r_find_workspace_root()
) {
  f6r_require_packages(c("ggplot2", "patchwork", "scales", "magick"))
  selected_pair_labels <- f6r_family_levels()
  selected_display_labels <- stats::setNames(
    selected_pair_labels, selected_pair_labels
  )
  paths <- f6r_paths(workspace_root)
  ksel <- f6r_read_tsv(file.path(paths$figure6, "cluster_k_selection.tsv"))
  bootstrap <- f6r_read_tsv(file.path(
    paths$figure6, "cluster_k_resample_selection_frequency.tsv"
  ))
  acceptance <- f6r_read_tsv(file.path(paths$figure6, "joint_seed_acceptance.tsv"))
  vivo <- f6r_read_tsv(file.path(paths$figure6, "joint_seed_claim_robustness.tsv"))
  vivo_dedup <- f6r_read_tsv(file.path(paths$figure6, "joint_seed_vs_unique_parameter_robustness.tsv"))
  vitro <- f6r_read_tsv(file.path(paths$figure6, "joint_seed_claim_robustness_invitro.tsv"))
  vitro_dedup <- f6r_read_tsv(file.path(paths$figure6, "joint_seed_vs_unique_parameter_robustness_invitro.tsv"))
  acceptance <- acceptance[acceptance$pair_label %in% selected_pair_labels, , drop = FALSE]
  vivo <- vivo[vivo$pair_label %in% selected_pair_labels, , drop = FALSE]
  vitro <- vitro[vitro$pair_label %in% selected_pair_labels, , drop = FALSE]
  vivo_dedup <- vivo_dedup[
    vivo_dedup$pair_id %in% unique(vivo$pair_id), , drop = FALSE
  ]
  vitro_dedup <- vitro_dedup[
    vitro_dedup$pair_id %in% unique(vitro$pair_id), , drop = FALSE
  ]
  primary <- ksel[ksel$analysis_level == "primary pooled t-SNE regions", ]
  selected_k_frequency <- paste0(
    "k=", bootstrap$selected_k, ": ", bootstrap$n_subsample, "/",
    sum(bootstrap$n_subsample), collapse = "; "
  )
  p_a <- ggplot2::ggplot(primary, ggplot2::aes(k, average_silhouette)) +
    ggplot2::geom_line(colour = "#555555", linewidth = 0.55) +
    ggplot2::geom_point(
      ggplot2::aes(fill = selected_for_warm_starts), shape = 21,
      size = 2.4, colour = "#222222", stroke = 0.45
    ) +
    ggplot2::annotate(
      "text", x = f6r_family_count(),
      y = primary$average_silhouette[
        primary$k == f6r_family_count()
      ] - 0.015,
      label = paste0("saved k=", f6r_family_count()), size = 2.45
    ) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#0072B2"), guide = "none") +
    ggplot2::scale_x_continuous(breaks = 2:8) +
    ggplot2::labs(
      tag = "A", title = "Shared primary warm-start-region selection",
      subtitle = paste0(
        "The joint workflow uses the saved k=", f6r_family_count(),
        " primary partition.\n",
        "80% subsample silhouette maxima: ", selected_k_frequency, "."
      ), x = "Number of regions (k)", y = "Average silhouette"
    ) + f6r_theme()
  acceptance$delta_display <- pmax(acceptance$delta_objective, 1e-4)
  counts <- aggregate(
    seed_number ~ pair_label, acceptance[acceptance$eligible_q10, ], length
  )
  acceptance$display_label <- unname(
    selected_display_labels[acceptance$pair_label]
  )
  acceptance$facet_label <- paste0(
    acceptance$display_label,
    "\n", counts$seed_number[match(acceptance$pair_label, counts$pair_label)],
    " seeds"
  )
  p_b <- ggplot2::ggplot(acceptance, ggplot2::aes(objective_rank, delta_display)) +
    ggplot2::annotate("rect", xmin = 0.5, xmax = 50.5, ymin = 1e-4, ymax = Inf,
                      fill = "#0072B2", alpha = 0.10) +
    ggplot2::geom_line(colour = "#666666", linewidth = 0.42) +
    ggplot2::geom_vline(
      data = data.frame(
        xintercept = c(25.5, 50.5, 100.5),
        cutoff_linetype = c("dotted", "solid", "dashed")
      ),
      ggplot2::aes(xintercept = xintercept, linetype = cutoff_linetype),
      linewidth = 0.30, show.legend = FALSE
    ) +
    ggplot2::scale_linetype_identity() +
    ggplot2::facet_wrap(~facet_label, nrow = 2L, scales = "free_y") +
    ggplot2::scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
    ggplot2::scale_x_continuous(breaks = c(1, 100, 250, 500)) +
    ggplot2::labs(
      tag = "B", title = "Selected-family joint-objective eligibility",
      subtitle = "Blue region: lowest 10%; vertical lines mark 5%, 10%, and 20% sets.",
      x = "Objective rank within warm-start pair",
      y = expression(Delta*" joint full-MAP objective (log10)")
    ) + f6r_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6.2))
  prepare_robust <- function(robustness, dedup, context) {
    z <- robustness[robustness$cutoff == "q10", , drop = FALSE]
    dz <- dedup[dedup$cutoff == "q10", c("pair_id", "claim_id", "same_modal_result")]
    z <- merge(z, dz, by = c("pair_id", "claim_id"), all.x = TRUE)
    z$model_context <- context
    z
  }
  robust <- rbind(
    prepare_robust(vivo, vivo_dedup, "in vivo"),
    prepare_robust(vitro, vitro_dedup, "in vitro")
  )
  robust$model_context <- factor(robust$model_context, levels = c("in vivo", "in vitro"))
  robust$display_label <- factor(
    unname(selected_display_labels[as.character(robust$pair_label)]),
    levels = selected_pair_labels
  )
  claim_order <- unique(vivo$claim_label)
  robust$claim_label <- factor(robust$claim_label, levels = rev(claim_order))
  robust$display_result <- ifelse(
    robust$modal_result == "TRUE", "yes",
    ifelse(robust$modal_result == "FALSE", "no", robust$modal_result)
  )
  robust$display <- paste0(
    robust$display_result, "\n",
    sprintf("%.0f%%", 100 * robust$modal_support_fraction)
  )
  expected_robustness_rows_per_context <-
    8L * f6r_family_count()
  p_c <- ggplot2::ggplot(
    robust,
    ggplot2::aes(display_label, claim_label, fill = modal_support_fraction)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.30) +
    ggplot2::geom_text(ggplot2::aes(label = display), size = 1.70, lineheight = 0.82) +
    ggplot2::geom_point(
      data = robust[!robust$same_modal_result, ], shape = 4,
      size = 1.7, stroke = 0.5, colour = "#B2182B"
    ) +
    ggplot2::facet_grid(model_context ~ ., switch = "y") +
    ggplot2::scale_fill_gradient(
      low = "#F2F2F2", high = "#0072B2", limits = c(0.5, 1),
      name = "Modal support within\nthe top-50 ensemble"
    ) +
    ggplot2::scale_y_discrete(labels = scales::label_wrap(34)) +
    ggplot2::labs(
      tag = "C", title = "Context-specific top-50 endpoint robustness",
      subtitle = paste0(
        "Text gives the modal result and exact support; a red cross marks a change after ",
        "collapsing duplicate parameter endpoints."
      ),
      x = "Selected warm-start family", y = NULL
    ) + f6r_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 5.8),
      strip.text.y = ggplot2::element_text(face = "bold")
    )
  combined <- ((p_a | p_b) + patchwork::plot_layout(widths = c(0.68, 1.32))) /
    p_c + patchwork::plot_layout(heights = c(0.78, 1.42)) +
    patchwork::plot_annotation(
      caption = paste0(
        "A records the ", f6r_family_count(),
        " primary warm-start regions. B-C include ",
        paste(f6r_family_levels(), collapse = ", "),
        ", the pairs displayed in Figure 6; optimizer seeds are not biological replicates."
      )
    )
  output <- f6r_save_plot(
    combined,
    file.path(paths$supp6_2, "supp_fig6-2_joint_ensemble_robustness"),
    width = 23.0, height = 8.8
  )
  published <- c(
    figures_png = f6r_publish(output[["png"]], file.path(paths$figures, "supp_fig6-2_joint_ensemble_robustness.png")),
    figures_pdf = f6r_publish(output[["pdf"]], file.path(paths$figures, "supp_fig6-2_joint_ensemble_robustness.pdf")),
    manuscript_png = f6r_publish(output[["png"]], file.path(paths$manuscript_figures, "supp_fig6-2_joint_ensemble_robustness.png")),
    manuscript_pdf = f6r_publish(output[["pdf"]], file.path(paths$manuscript_figures, "supp_fig6-2_joint_ensemble_robustness.pdf"))
  )
  validation <- data.frame(
    check = c(
      "panel_count", "selected_pair_count", "display_labels", "context_count_panel_c",
      "robustness_rows_per_context", "visible_subcluster_labels_absent",
      "publication_png_md5", "publication_pdf_md5"
    ),
    observed = c(
      3, length(unique(acceptance$pair_label)),
      paste(sort(unique(acceptance$display_label)), collapse = ","),
      length(unique(robust$model_context)),
      paste(as.integer(table(robust$model_context)), collapse = ","),
      !any(grepl("Sc[0-9]+", levels(robust$display_label))),
      f6r_md5(output[["png"]]) == f6r_md5(published[["manuscript_png"]]),
      f6r_md5(output[["pdf"]]) == f6r_md5(published[["manuscript_pdf"]])
    ),
    expected = c(
      3, f6r_family_count(), paste(selected_pair_labels, collapse = ","), 2,
      paste(
        rep(expected_robustness_rows_per_context, 2L),
        collapse = ","
      ),
      TRUE, TRUE, TRUE
    ),
    stringsAsFactors = FALSE
  )
  validation$passed <- as.character(validation$observed) == as.character(validation$expected)
  f6r_write_tsv(validation, file.path(paths$supp6_2, "figure_validation.tsv"))
  if (!all(validation$passed)) stop("Supplementary Figure 6-2 validation failed.")
  f6x_refresh_output_manifest(file.path(
    paths$supp6_2, "supp_figure6-2_output_manifest.tsv"
  ))
  invisible(list(output = output, published = published, validation = validation))
}
