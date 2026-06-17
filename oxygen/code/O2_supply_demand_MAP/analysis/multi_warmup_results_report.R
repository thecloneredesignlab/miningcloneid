#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos > 0L) out[[substr(kv, 1L, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv)) else out[[kv]] <- TRUE
  }
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) y else x
as_chr <- function(x, default = "") as.character((x %||% default)[[1]])

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
}

SCRIPT_DIR <- normalizePath(get_script_dir(), mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "report", "run_provenance_report.R"), local = environment())

escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

read_table_optional <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}

num_col <- function(x) suppressWarnings(as.numeric(x))

format_unique_values <- function(x, digits = 6L, scientific = FALSE) {
  if (is.null(x) || !length(x)) return("NA")
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return("NA")
  num <- suppressWarnings(as.numeric(x))
  if (all(is.finite(num))) {
    vals <- unique(signif(num, digits))
    vals <- vals[order(vals)]
    return(paste(format(vals, scientific = scientific, trim = TRUE), collapse = ", "))
  }
  paste(unique(as.character(x)), collapse = ", ")
}

format_count_summary <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  ux <- sort(unique(x))
  if (length(ux) == 1L) return(as.character(ux[[1]]))
  paste0(min(ux), "-", max(ux), " (mixed)")
}

strip_project_path_prefix <- function(x) {
  prefix <- "/share/lab_crd/lab_crd/taoli/Project"
  if (!length(x)) return(x)
  gsub(prefix, "", as.character(x), fixed = TRUE)
}

strip_project_path_prefix_table <- function(tab) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return(tab)
  out <- tab
  chr_cols <- vapply(out, is.character, logical(1))
  for (col in names(out)[chr_cols]) {
    out[[col]] <- strip_project_path_prefix(out[[col]])
  }
  out
}

drop_report_columns <- function(tab, cols) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return(tab)
  tab[, setdiff(names(tab), cols), drop = FALSE]
}

table_html <- function(tab, max_rows = 50L) {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return('<p class="empty">No data available.</p>')
  tab <- utils::head(tab, max_rows)
  header <- paste0("<tr>", paste(sprintf("<th>%s</th>", escape_html(names(tab))), collapse = ""), "</tr>")
  rows <- apply(tab, 1L, function(row) paste0("<tr>", paste(sprintf("<td>%s</td>", escape_html(row)), collapse = ""), "</tr>"))
  paste0('<table class="report-table">', header, paste(rows, collapse = ""), "</table>")
}

file_to_data_uri <- function(path, mime = "application/pdf") {
  if (requireNamespace("base64enc", quietly = TRUE)) {
    encoded <- base64enc::base64encode(path)
  } else if (requireNamespace("openssl", quietly = TRUE)) {
    bytes <- readBin(path, what = "raw", n = file.info(path)$size)
    encoded <- openssl::base64_encode(bytes)
  } else {
    return(NULL)
  }
  paste0("data:", mime, ";base64,", encoded)
}

magick_available <- function() {
  requireNamespace("magick", quietly = TRUE)
}

ghostscript_bin <- function() {
  candidates <- c(Sys.which("gs"), Sys.which("gswin64c"), Sys.which("gswin32c"))
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) candidates[[1]] else ""
}

ghostscript_available <- function() {
  nzchar(ghostscript_bin())
}

render_pdf_preview_png_gs <- function(src_pdf, dest_png, density = 180) {
  gs_bin <- ghostscript_bin()
  if (!nzchar(gs_bin)) {
    stop("Ghostscript is not available.", call. = FALSE)
  }
  src_pdf_use <- normalizePath(src_pdf, mustWork = TRUE)
  dest_png_use <- normalizePath(dest_png, mustWork = FALSE)
  density_use <- suppressWarnings(as.integer(density))
  if (!is.finite(density_use) || density_use <= 0L) density_use <- 180L
  args <- c(
    "-dSAFER",
    "-dBATCH",
    "-dNOPAUSE",
    "-sDEVICE=pngalpha",
    sprintf("-r%d", density_use),
    sprintf("-sOutputFile=%s", shQuote(dest_png_use)),
    shQuote(src_pdf_use)
  )
  out <- suppressWarnings(system2(gs_bin, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("Ghostscript failed while rendering ", src_pdf, ": ", paste(out, collapse = "\n"), call. = FALSE)
  }
  if (!file.exists(dest_png_use)) {
    stop("Ghostscript did not create expected PNG preview: ", dest_png_use, call. = FALSE)
  }
  normalizePath(dest_png_use, mustWork = TRUE)
}

pdf_to_image_data_uri <- function(pdf_path, density = 180) {
  if (!magick_available() && !ghostscript_available()) return(NULL)
  png_path <- tempfile("multi_warmup_pdf_preview_", fileext = ".png")
  on.exit(unlink(png_path, force = TRUE), add = TRUE)
  tryCatch({
    if (magick_available()) {
      img <- magick::image_read(pdf_path, density = density)
      if (length(img) > 1L) img <- img[1]
      magick::image_write(img, path = png_path, format = "png")
    } else {
      render_pdf_preview_png_gs(pdf_path, png_path, density = density)
    }
    file_to_data_uri(png_path, mime = "image/png")
  }, error = function(e) NULL)
}

figure_html <- function(root, filename, title, caption, section_id = NULL) {
  path <- file.path(root, filename)
  if (!file.exists(path)) return("")
  section_attr <- if (!is.null(section_id) && nzchar(as.character(section_id))) {
    paste0(' id="', escape_html(section_id), '"')
  } else {
    ""
  }
  figure_data <- pdf_to_image_data_uri(path)
  figure_body <- if (!is.null(figure_data) && nzchar(figure_data)) {
    paste0(
      '<img class="figure-image" src="', escape_html(figure_data), '" alt="', escape_html(title), '"/>'
    )
  } else {
    paste0('<p class="empty">PDF preview rendering is unavailable. Open <code>', escape_html(filename), '</code>.</p>')
  }
  paste0(
    '<section class="card"', section_attr, '>',
    '<h2>', escape_html(title), '</h2>',
    '<p>', escape_html(caption), '</p>',
    figure_body,
    '</section>'
  )
}

insert_after_first_regex <- function(x, pattern, insert) {
  if (!nzchar(insert)) return(x)
  hit <- regexpr(pattern, x, perl = TRUE)
  if (hit[[1]] < 0L) return(x)
  end <- hit[[1]] + attr(hit, "match.length") - 1L
  paste0(substr(x, 1L, end), insert, substr(x, end + 1L, nchar(x)))
}

plot_integrated_tradeoff_by_invivo_warmup <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  summary_path <- file.path(root, "multi_warmup_integrated_joint_run", "extra_results", "seed_summary.tsv")
  summary_df <- read_table_optional(summary_path)
  required <- c("seed", "objective_invivo", "objective_invitro")
  if (is.null(summary_df) || !is.data.frame(summary_df) || !all(required %in% names(summary_df))) {
    return(invisible(NULL))
  }
  plot_df <- summary_df[, required, drop = FALSE]
  plot_df$objective_invivo <- num_col(plot_df$objective_invivo)
  plot_df$objective_invitro <- num_col(plot_df$objective_invitro)
  plot_df$warmup_label <- sub("__.*$", "", as.character(plot_df$seed))
  plot_df$invivo_warmup <- toupper(sub("_.*$", "", plot_df$warmup_label))
  plot_df <- plot_df[
    is.finite(plot_df$objective_invivo) &
      is.finite(plot_df$objective_invitro) &
      nzchar(plot_df$invivo_warmup),
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))

  manifest <- read_table_optional(file.path(root, "multi_warmup_manifest.tsv"))
  if (is.data.frame(manifest) && nrow(manifest) && "warmup_label" %in% names(manifest)) {
    warmup_levels <- unique(toupper(sub("_.*$", "", as.character(manifest$warmup_label))))
    warmup_levels <- warmup_levels[nzchar(warmup_levels)]
  } else {
    warmup_levels <- unique(plot_df$invivo_warmup)
  }
  plot_df$invivo_warmup <- factor(plot_df$invivo_warmup, levels = warmup_levels)

  plot_data_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.tsv")
  utils::write.table(plot_df, plot_data_path, sep = "\t", quote = FALSE, row.names = FALSE)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = objective_invivo, y = objective_invitro, color = invivo_warmup)
  ) +
    ggplot2::geom_point(size = 1, alpha = 0.9) +
    ggplot2::labs(
      title = "Joint Objective Tradeoff by In Vivo Warm-Up Family",
      subtitle = "Each point is one integrated joint seed; color shows Vxx warm-up source.",
      x = "In vivo objective",
      y = "In vitro objective",
      color = "In vivo warm-up"
    ) +
    scale_color_invivo_warmup(warmup_levels) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  pdf_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.pdf")
  png_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.png")
  ggplot2::ggsave(pdf_path, p, width = 9, height = 7, bg = "white")
  ggplot2::ggsave(png_path, p, width = 9, height = 7, dpi = 220, bg = "white")
  invisible(pdf_path)
}

integrated_tradeoff_by_invivo_warmup_html <- function(root) {
  png_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.png")
  pdf_path <- file.path(root, "multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.pdf")
  if (!file.exists(png_path) && !file.exists(pdf_path)) return("")
  image_data <- if (file.exists(png_path)) {
    file_to_data_uri(png_path, mime = "image/png")
  } else {
    pdf_to_image_data_uri(pdf_path)
  }
  figure_body <- if (!is.null(image_data) && nzchar(image_data)) {
    paste0(
      '<div class="report-figure"><img src="',
      escape_html(image_data),
      '" alt="Joint In Vivo vs In Vitro Objective Tradeoff by In Vivo Warm-Up Family" class="report-figure-image"/></div>'
    )
  } else {
    '<p class="report-empty">Plot preview rendering is unavailable.</p>'
  }
  paste0(
    '<section class="report-section" id="integrated-joint-figure-4b">',
    figure_body,
    '<h2 class="report-figure-title">Figure 2.4b. Joint In Vivo vs In Vitro Objective Tradeoff by In Vivo Warm-Up Family</h2>',
    '<p class="report-figure-legend">Seed-level in vivo objective against in vitro objective; color shows Vxx in vivo warm-up source. Points are one third of the Figure 2.4 point size and are not text-labeled.</p>',
    '<p class="report-figure-file"><code>multi_warmup_integrated_objective_tradeoff_by_invivo_warmup.pdf</code></p>',
    '</section>'
  )
}

integrated_extra_results_dir <- function(root) {
  file.path(root, "multi_warmup_integrated_joint_run", "extra_results")
}

warmup_levels_from_manifest <- function(root, fallback = character(0)) {
  manifest <- read_table_optional(file.path(root, "multi_warmup_manifest.tsv"))
  if (is.data.frame(manifest) && nrow(manifest) && "warmup_label" %in% names(manifest)) {
    out <- unique(toupper(sub("_.*$", "", as.character(manifest$warmup_label))))
    out <- out[nzchar(out)]
    if (length(out)) return(out)
  }
  unique(as.character(fallback[nzchar(fallback)]))
}

add_invivo_warmup_from_seed <- function(tab, root, seed_col = "seed") {
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab) || !(seed_col %in% names(tab))) return(tab)
  tab$warmup_label <- sub("__.*$", "", as.character(tab[[seed_col]]))
  tab$invivo_warmup <- toupper(sub("_.*$", "", tab$warmup_label))
  levels <- warmup_levels_from_manifest(root, fallback = tab$invivo_warmup)
  tab$invivo_warmup <- factor(tab$invivo_warmup, levels = levels)
  tab
}

invivo_warmup_palette <- function(levels) {
  levels <- as.character(levels)
  levels <- levels[nzchar(levels)]
  palette <- c(
    "#009E73", # green
    "#6A3D9A", # purple
    "#E69F00", # orange
    "#000000", # black
    "#8C6D31", # dark brown
    "#F0E442", # yellow
    "#66A61E", # olive green
    "#7F7F7F", # grey
    "#B07AA1", # muted purple
    "#A6761D"  # brown
  )
  if (length(levels) > length(palette)) {
    palette <- grDevices::colorRampPalette(palette)(length(levels))
  } else {
    palette <- palette[seq_along(levels)]
  }
  stats::setNames(palette, levels)
}

scale_color_invivo_warmup <- function(levels, name = "In vivo warm-up", breaks = levels, drop = FALSE) {
  ggplot2::scale_color_manual(
    values = invivo_warmup_palette(levels),
    limits = levels,
    breaks = breaks,
    drop = drop,
    name = name
  )
}

save_plot_pdf_png <- function(plot, root, stem, width, height, dpi = 220) {
  pdf_path <- file.path(root, paste0(stem, ".pdf"))
  png_path <- file.path(root, paste0(stem, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(c(pdf = pdf_path, png = png_path))
}

point_shape_values <- function(n) {
  base <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 5, 6, 9, 10, 12, 13, 14)
  rep(base, length.out = n)
}

square_limits <- function(x, y, pad_frac = 0.08) {
  x <- num_col(x)
  y <- num_col(y)
  xr <- range(x[is.finite(x)], na.rm = TRUE)
  yr <- range(y[is.finite(y)], na.rm = TRUE)
  cx <- mean(xr)
  cy <- mean(yr)
  span <- max(diff(xr), diff(yr))
  if (!is.finite(span) || span <= 0) span <- 1
  half <- span * (0.5 + pad_frac)
  list(x = c(cx - half, cx + half), y = c(cy - half, cy + half))
}

wrap_label <- function(x, width = 68L) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

plot_part1_objective_figures <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  summary_df <- read_table_optional(file.path(root, "multi_warmup_best_seed_summary.tsv"))
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df)) return(invisible(NULL))
  if (!all(c("warmup_label", "objective", "invivo_family") %in% names(summary_df))) return(invisible(NULL))
  summary_df$objective <- num_col(summary_df$objective)
  family_levels <- unique(as.character(summary_df$invivo_family))
  family_levels <- family_levels[nzchar(family_levels)]
  summary_df$invivo_family <- factor(as.character(summary_df$invivo_family), levels = family_levels)

  p1 <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = stats::reorder(warmup_label, objective), y = objective, color = invivo_family)
  ) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::coord_flip() +
    scale_color_invivo_warmup(family_levels, name = "In vivo family") +
    ggplot2::labs(
      title = "Multi-Warm-Up Best Joint Objective",
      x = "Warm-up pair",
      y = "Best total objective"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  save_plot_pdf_png(p1, root, "multi_warmup_objective_summary", width = 8, height = 5)

  if (all(c("objective_invivo", "objective_invitro") %in% names(summary_df))) {
    summary_df$objective_invivo <- num_col(summary_df$objective_invivo)
    summary_df$objective_invitro <- num_col(summary_df$objective_invitro)
    plot_df <- summary_df[is.finite(summary_df$objective_invivo) & is.finite(summary_df$objective_invitro), , drop = FALSE]
    if (nrow(plot_df)) {
      p2 <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(objective_invivo, objective_invitro, color = invivo_family, label = warmup_label)
      ) +
        ggplot2::geom_point(size = 2.8) +
        scale_color_invivo_warmup(family_levels, name = "In vivo family") +
        ggplot2::labs(
          title = "In Vivo vs In Vitro Objective by Warm-Up Pair",
          x = "Best seed in vivo objective",
          y = "Best seed in vitro objective"
        ) +
        ggplot2::theme_minimal(base_size = 11)
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p2 <- p2 + ggrepel::geom_text_repel(size = 3, max.overlaps = Inf)
      } else {
        p2 <- p2 + ggplot2::geom_text(size = 3, vjust = -0.6)
      }
      save_plot_pdf_png(p2, root, "multi_warmup_invivo_invitro_objective_scatter", width = 7, height = 5)
    }
  }
  invisible(NULL)
}

plot_part1_cluster_umap <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  coords <- read_table_optional(file.path(root, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_coords.tsv"))
  required <- c("pair_id", "invivo_cluster", "invitro_rank", "UMAP1", "UMAP2")
  if (is.null(coords) || !is.data.frame(coords) || !all(required %in% names(coords))) return(invisible(NULL))
  plot_df <- coords
  cluster_levels <- unique(as.character(plot_df$invivo_cluster))
  cluster_levels <- cluster_levels[nzchar(cluster_levels)]
  plot_df$invivo_cluster <- factor(as.character(plot_df$invivo_cluster), levels = cluster_levels)
  top_n <- max(suppressWarnings(as.integer(plot_df$invitro_rank)), na.rm = TRUE)
  if (!is.finite(top_n) || top_n <= 0L) top_n <- length(unique(plot_df$invitro_rank))
  plot_df$invitro_rank_factor <- factor(plot_df$invitro_rank, levels = seq_len(top_n))
  if (!("is_invivo_representative" %in% names(plot_df))) plot_df$is_invivo_representative <- FALSE
  plot_df$is_invivo_representative <- tolower(as.character(plot_df$is_invivo_representative)) %in% c("true", "t", "1", "yes")
  limits <- square_limits(plot_df$UMAP1, plot_df$UMAP2)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(UMAP1, UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = invivo_cluster, shape = invitro_rank_factor), size = 3.0, stroke = 0.9) +
    scale_color_invivo_warmup(cluster_levels, name = "In Vivo Cluster") +
    ggplot2::scale_shape_manual(values = point_shape_values(top_n), name = "In Vitro Rank") +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3.2)),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::labs(
      title = "18D Warm-Start Profile UMAP by In Vivo Cluster",
      subtitle = wrap_label(
        "Same 18D paired warm-start profile UMAP as the standalone view, colored by profile-space in vivo cluster. Red labels mark selected representative in vivo ranks.",
        width = 68L
      ),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    ggplot2::coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE, clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical"
    )
  label_nonrep <- plot_df[!plot_df$is_invivo_representative, , drop = FALSE]
  label_rep <- plot_df[plot_df$is_invivo_representative, , drop = FALSE]
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    if (nrow(label_nonrep)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_nonrep,
        ggplot2::aes(label = pair_id),
        size = 2.1,
        color = "grey15",
        max.overlaps = Inf,
        seed = 1L
      )
    }
    if (nrow(label_rep)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_rep,
        ggplot2::aes(label = pair_id),
        size = 2.35,
        color = "#d62728",
        fontface = "bold",
        max.overlaps = Inf,
        seed = 1L
      )
    }
  } else {
    if (nrow(label_nonrep)) {
      p <- p + ggplot2::geom_text(data = label_nonrep, ggplot2::aes(label = pair_id), size = 2.1, color = "grey15", vjust = -0.6)
    }
    if (nrow(label_rep)) {
      p <- p + ggplot2::geom_text(data = label_rep, ggplot2::aes(label = pair_id), size = 2.35, color = "#d62728", fontface = "bold", vjust = -0.6)
    }
  }
  save_plot_pdf_png(p, root, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_500seed", width = 7, height = 7)
  invisible(NULL)
}

plot_part1_color_aligned_figures <- function(root) {
  plot_part1_cluster_umap(root)
  plot_part1_objective_figures(root)
  invisible(NULL)
}

mw_boundary_axis_config <- function(x,
                                    near_thresh,
                                    x_scale = c("relative", "log10_original"),
                                    raw_value = NULL,
                                    raw_lower = NULL,
                                    raw_upper = NULL) {
  x_scale <- match.arg(x_scale)
  if (identical(x_scale, "relative")) {
    return(list(
      axis_type = "relative",
      x = x,
      lower_limit = 0,
      upper_limit = 1,
      vlines = c(0, near_thresh, 0.5, 1 - near_thresh, 1),
      vline_colors = c("grey50", "grey70", "grey82", "grey70", "grey50"),
      vline_linetypes = c("solid", "dashed", "dotted", "dashed", "solid"),
      lower_rect = c(0, near_thresh),
      upper_rect = c(1 - near_thresh, 1),
      scale = ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, near_thresh, 0.5, 1 - near_thresh, 1)),
      x_label = "Relative position in transformed fit range",
      subtitle_note = ""
    ))
  }

  raw_value <- suppressWarnings(as.numeric(raw_value))
  raw_lower <- suppressWarnings(as.numeric(raw_lower))
  raw_upper <- suppressWarnings(as.numeric(raw_upper))
  positive_x <- c(raw_value, raw_lower, raw_upper)
  positive_x <- positive_x[is.finite(positive_x) & positive_x > 0]
  log_floor <- min(positive_x, na.rm = TRUE) / 10
  if (!is.finite(log_floor) || log_floor <= 0) log_floor <- 1e-6
  lower_plot <- ifelse(is.finite(raw_lower), pmax(raw_lower, log_floor), NA_real_)
  upper_plot <- ifelse(is.finite(raw_upper), pmax(raw_upper, log_floor), NA_real_)
  value_plot <- ifelse(is.finite(raw_value), pmax(raw_value, log_floor), NA_real_)
  axis_upper <- max(c(upper_plot, value_plot), na.rm = TRUE)
  if (!is.finite(axis_upper) || axis_upper <= log_floor) axis_upper <- log_floor * 10
  breaks <- 10^seq(floor(log10(log_floor)), ceiling(log10(axis_upper)), by = 1)
  breaks <- breaks[breaks >= log_floor & breaks <= axis_upper]
  needs_floor_label <- any(
    is.finite(c(raw_value, raw_lower, raw_upper)) & c(raw_value, raw_lower, raw_upper) <= 0,
    na.rm = TRUE
  )
  labels <- function(vals) {
    vapply(
      vals,
      function(v) {
        if (isTRUE(needs_floor_label) && isTRUE(all.equal(v, log_floor, tolerance = 1e-12))) return("0")
        formatC(v, format = "fg", digits = 4)
      },
      character(1)
    )
  }
  list(
    axis_type = "log10_original",
    x = value_plot,
    lower_limit = log_floor,
    upper_limit = axis_upper,
    lower_plot = lower_plot,
    upper_plot = upper_plot,
    scale = ggplot2::scale_x_log10(limits = c(log_floor, axis_upper), breaks = breaks, labels = labels),
    x_label = "Original parameter value (log10 scale)",
    subtitle_note = if (isTRUE(needs_floor_label)) " Values <=0 shown at log floor." else ""
  )
}

mw_top_ranked_seeds <- function(summary_df, n = 3L, pred1000_gate = FALSE) {
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df) || !all(c("seed", "objective") %in% names(summary_df))) {
    return(character(0))
  }
  keep <- is.finite(num_col(summary_df$objective))
  if (isTRUE(pred1000_gate) && "pred1000_both_gt44" %in% names(summary_df)) {
    gate <- summary_df$pred1000_both_gt44
    if (!is.logical(gate)) gate <- tolower(as.character(gate)) %in% c("true", "t", "1", "yes")
    keep <- keep & isTRUE(TRUE) & gate
  }
  ranked <- summary_df[keep, , drop = FALSE]
  if (!nrow(ranked)) return(character(0))
  ranked$objective_num <- num_col(ranked$objective)
  ranked <- ranked[order(ranked$objective_num), , drop = FALSE]
  head(as.character(ranked$seed), n)
}

plot_multi_warmup_boundary_forest <- function(root,
                                              stem,
                                              title_suffix = NULL,
                                              pred1000_top3_only = FALSE,
                                              x_scale = c("relative", "log10_original"),
                                              near_thresh = 0.05) {
  x_scale <- match.arg(x_scale)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  extra_dir <- integrated_extra_results_dir(root)
  long_df <- read_table_optional(file.path(extra_dir, "parameter_boundary_long.tsv"))
  summary_df <- read_table_optional(file.path(extra_dir, "seed_summary.tsv"))
  required <- c(
    "seed",
    "active_in_fit",
    "rel_pos_plot",
    "rel_dist_to_nearest",
    "param_prototype",
    "prototype_value",
    "prototype_lower_bound",
    "prototype_upper_bound"
  )
  if (is.null(long_df) || !is.data.frame(long_df) || !all(required %in% names(long_df))) return(invisible(NULL))
  plot_df <- long_df[as.logical(long_df$active_in_fit) & is.finite(num_col(long_df$rel_pos_plot)), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  top3_seeds <- mw_top_ranked_seeds(summary_df, n = 3L, pred1000_gate = pred1000_top3_only)
  plot_df <- add_invivo_warmup_from_seed(plot_df, root)
  plot_df$rel_pos_plot <- num_col(plot_df$rel_pos_plot)
  plot_df$rel_dist_to_nearest <- num_col(plot_df$rel_dist_to_nearest)
  plot_df$prototype_value <- num_col(plot_df$prototype_value)
  plot_df$prototype_lower_bound <- num_col(plot_df$prototype_lower_bound)
  plot_df$prototype_upper_bound <- num_col(plot_df$prototype_upper_bound)
  axis_cfg <- mw_boundary_axis_config(
    plot_df$rel_pos_plot,
    near_thresh = near_thresh,
    x_scale = x_scale,
    raw_value = plot_df$prototype_value,
    raw_lower = plot_df$prototype_lower_bound,
    raw_upper = plot_df$prototype_upper_bound
  )
  plot_df$boundary_x_plot <- axis_cfg$x
  if (identical(axis_cfg$axis_type, "log10_original")) {
    plot_df$boundary_x_lower <- axis_cfg$lower_plot
    plot_df$boundary_x_upper <- axis_cfg$upper_plot
    plot_df <- plot_df[
      is.finite(plot_df$boundary_x_plot) &
        is.finite(plot_df$boundary_x_lower) &
        is.finite(plot_df$boundary_x_upper) &
        plot_df$boundary_x_upper > plot_df$boundary_x_lower,
      ,
      drop = FALSE
    ]
    if (!nrow(plot_df)) return(invisible(NULL))
  }
  param_rank <- tapply(plot_df$rel_dist_to_nearest, plot_df$param_prototype, min, na.rm = TRUE)
  param_levels <- names(sort(param_rank, decreasing = FALSE))
  plot_df$param_prototype <- factor(plot_df$param_prototype, levels = rev(param_levels))
  point_pos <- ggplot2::position_jitter(height = 0.14, width = 0, seed = 1L)
  title_detail <- if (!is.null(title_suffix) && nzchar(title_suffix)) paste0(" (", title_suffix, ")") else ""
  top3_label <- if (length(top3_seeds)) {
    if (isTRUE(pred1000_top3_only)) {
      " Black symbols mark eligible objective top 3."
    } else {
      " Black symbols mark objective top 3."
    }
  } else if (isTRUE(pred1000_top3_only)) {
    "No seeds met the 2N/4N 1000-day prediction gate."
  } else {
    ""
  }
  color_levels <- levels(plot_df$invivo_warmup)
  color_breaks <- color_levels
  point_size <- 2.0
  point_alpha <- 0.62
  top_df <- plot_df[as.character(plot_df$seed) %in% top3_seeds, , drop = FALSE]
  if (nrow(top_df)) {
    marker_levels <- c(
      if (length(top3_seeds) >= 1L) paste0("Top 1: ", top3_seeds[[1]], " (*)"),
      if (length(top3_seeds) >= 2L) paste0("Top 2: ", top3_seeds[[2]], " (triangle)"),
      if (length(top3_seeds) >= 3L) paste0("Top 3: ", top3_seeds[[3]], " (black dot)")
    )
    names(marker_levels) <- top3_seeds
    top_df$top3_marker <- factor(marker_levels[as.character(top_df$seed)], levels = unname(marker_levels))
    top_shape_values <- c(8, 17, 16)[seq_along(marker_levels)]
    names(top_shape_values) <- unname(marker_levels)
  } else {
    top_shape_values <- numeric(0)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = boundary_x_plot, y = param_prototype))
  if (identical(axis_cfg$axis_type, "relative")) {
    ref_df <- unique(plot_df["param_prototype"])
    ref_df$boundary_x_start <- axis_cfg$lower_limit
    ref_df$boundary_x_end <- axis_cfg$upper_limit
    p <- p +
      ggplot2::annotate("rect", xmin = axis_cfg$lower_rect[[1]], xmax = axis_cfg$lower_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#fddbc7", alpha = 0.28) +
      ggplot2::annotate("rect", xmin = axis_cfg$upper_rect[[1]], xmax = axis_cfg$upper_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#d1e5f0", alpha = 0.28) +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = boundary_x_start, xend = boundary_x_end, y = param_prototype, yend = param_prototype),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      ) +
      ggplot2::geom_vline(xintercept = axis_cfg$vlines, color = axis_cfg$vline_colors, linetype = axis_cfg$vline_linetypes, linewidth = 0.35)
  } else {
    ref_df <- unique(plot_df[c("param_prototype", "boundary_x_lower", "boundary_x_upper")])
    p <- p +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = boundary_x_lower, xend = boundary_x_upper, y = param_prototype, yend = param_prototype),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      )
  }

  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(color = invivo_warmup),
      shape = 16,
      size = point_size,
      alpha = point_alpha,
      position = point_pos
    ) +
    axis_cfg$scale +
    ggplot2::labs(
      title = paste0("Parameter Positions Within Fitted Bounds", title_detail, ": integrated multi-warm-up"),
      subtitle = paste0(
        if (identical(axis_cfg$axis_type, "relative")) {
          paste0("0=lower, 1=upper; shaded zones are within ", sprintf("%.0f", 100 * near_thresh), "% of a bound. ")
        } else {
          "Horizontal lines span original bounds. "
        },
        "Points are colored by Vxx warm-up family.",
        top3_label,
        axis_cfg$subtitle_note
      ),
      x = axis_cfg$x_label,
      y = NULL,
      color = "In vivo warm-up"
    ) +
    scale_color_invivo_warmup(color_levels, breaks = color_breaks) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )
  if (nrow(top_df)) {
    p <- p +
      ggplot2::geom_point(
        data = top_df,
        ggplot2::aes(shape = top3_marker),
        color = "black",
        size = 3.0,
        position = point_pos
      ) +
      ggplot2::scale_shape_manual(values = top_shape_values, drop = FALSE) +
      ggplot2::labs(shape = "Top seed marker") +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1), order = 1),
        shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE, order = 2)
      )
  }
  save_plot_pdf_png(p, root, stem, width = 12, height = 5.2)
}

plot_multi_warmup_paired_boundary_scatter <- function(root,
                                                      stem,
                                                      title_suffix = NULL,
                                                      pred1000_top3_only = FALSE,
                                                      x_scale = c("relative", "log10_original"),
                                                      near_thresh = 0.05) {
  x_scale <- match.arg(x_scale)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  extra_dir <- integrated_extra_results_dir(root)
  soft_df <- read_table_optional(file.path(extra_dir, "joint_soft_coupling_all.tsv"))
  summary_df <- read_table_optional(file.path(extra_dir, "seed_summary.tsv"))
  required <- c(
    "seed",
    "parameter",
    "vivo_transformed",
    "vitro_transformed",
    "vivo_natural",
    "vitro_natural"
  )
  if (is.null(soft_df) || !is.data.frame(soft_df) || !nrow(soft_df) || !all(required %in% names(soft_df))) {
    return(invisible(NULL))
  }
  lower_t_col <- if ("joint_union_lower_transformed" %in% names(soft_df)) "joint_union_lower_transformed" else "center_lower_transformed"
  upper_t_col <- if ("joint_union_upper_transformed" %in% names(soft_df)) "joint_union_upper_transformed" else "center_upper_transformed"
  lower_nat_col <- if ("joint_union_lower_bound" %in% names(soft_df)) "joint_union_lower_bound" else "center_lower_bound"
  upper_nat_col <- if ("joint_union_upper_bound" %in% names(soft_df)) "joint_union_upper_bound" else "center_upper_bound"
  bound_cols <- c(lower_t_col, upper_t_col, lower_nat_col, upper_nat_col)
  if (!all(bound_cols %in% names(soft_df))) return(invisible(NULL))

  make_context_df <- function(context, value_t_col, value_nat_col) {
    data.frame(
      seed = as.character(soft_df$seed),
      parameter = as.character(soft_df$parameter),
      context = context,
      value_transformed = num_col(soft_df[[value_t_col]]),
      value_natural = num_col(soft_df[[value_nat_col]]),
      bound_lower_transformed = num_col(soft_df[[lower_t_col]]),
      bound_upper_transformed = num_col(soft_df[[upper_t_col]]),
      bound_lower_natural = num_col(soft_df[[lower_nat_col]]),
      bound_upper_natural = num_col(soft_df[[upper_nat_col]]),
      stringsAsFactors = FALSE
    )
  }
  plot_df <- rbind(
    make_context_df("in vivo", "vivo_transformed", "vivo_natural"),
    make_context_df("in vitro", "vitro_transformed", "vitro_natural")
  )
  width_t <- plot_df$bound_upper_transformed - plot_df$bound_lower_transformed
  plot_df$rel_pos_in_range <- (plot_df$value_transformed - plot_df$bound_lower_transformed) / width_t
  plot_df$rel_pos_plot <- pmin(pmax(plot_df$rel_pos_in_range, 0), 1)
  plot_df$rel_dist_to_lower <- (plot_df$value_transformed - plot_df$bound_lower_transformed) / width_t
  plot_df$rel_dist_to_upper <- (plot_df$bound_upper_transformed - plot_df$value_transformed) / width_t
  plot_df$rel_dist_to_nearest <- pmin(plot_df$rel_dist_to_lower, plot_df$rel_dist_to_upper)
  plot_df <- plot_df[
    is.finite(plot_df$rel_pos_plot) &
      is.finite(plot_df$rel_dist_to_nearest) &
      is.finite(width_t) &
      width_t > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))

  plot_df <- add_invivo_warmup_from_seed(plot_df, root)
  family_levels <- levels(plot_df$invivo_warmup)
  top3_seeds <- mw_top_ranked_seeds(summary_df, n = 3L, pred1000_gate = pred1000_top3_only)

  axis_cfg <- mw_boundary_axis_config(
    plot_df$rel_pos_plot,
    near_thresh = near_thresh,
    x_scale = x_scale,
    raw_value = plot_df$value_natural,
    raw_lower = plot_df$bound_lower_natural,
    raw_upper = plot_df$bound_upper_natural
  )
  plot_df$boundary_x_plot <- axis_cfg$x
  if (identical(axis_cfg$axis_type, "log10_original")) {
    plot_df$boundary_x_lower <- axis_cfg$lower_plot
    plot_df$boundary_x_upper <- axis_cfg$upper_plot
    plot_df <- plot_df[
      is.finite(plot_df$boundary_x_plot) &
        is.finite(plot_df$boundary_x_lower) &
        is.finite(plot_df$boundary_x_upper) &
        plot_df$boundary_x_upper > plot_df$boundary_x_lower,
      ,
      drop = FALSE
    ]
    if (!nrow(plot_df)) return(invisible(NULL))
  }

  plot_df$parameter_label <- as.character(plot_df$parameter)
  param_rank <- tapply(plot_df$rel_dist_to_nearest, plot_df$parameter_label, min, na.rm = TRUE)
  param_levels <- names(sort(param_rank, decreasing = FALSE))
  plot_df$parameter <- factor(plot_df$parameter_label, levels = rev(param_levels))
  y_breaks <- seq_along(levels(plot_df$parameter))
  plot_df$y_base <- as.numeric(plot_df$parameter)
  plot_df$context <- factor(plot_df$context, levels = c("in vivo", "in vitro"))
  pair_key <- paste(plot_df$seed, plot_df$parameter_label, sep = "\r")
  pair_hash <- vapply(
    pair_key,
    function(key) {
      ints <- utf8ToInt(key)
      if (!length(ints)) return(0)
      sum((seq_along(ints) * ints) %% 997)
    },
    numeric(1)
  )
  plot_df$pair_jitter <- ((pair_hash %% 101) / 100 - 0.5) * 0.08
  plot_df$context_offset <- ifelse(as.character(plot_df$context) == "in vivo", 0.18, -0.18)
  plot_df$y_plot <- plot_df$y_base + plot_df$context_offset + plot_df$pair_jitter
  ref_df <- if (identical(axis_cfg$axis_type, "log10_original")) {
    unique(plot_df[c("parameter", "y_base", "boundary_x_lower", "boundary_x_upper")])
  } else {
    ref <- unique(plot_df[c("parameter", "y_base")])
    ref$boundary_x_start <- axis_cfg$lower_limit
    ref$boundary_x_end <- axis_cfg$upper_limit
    ref
  }

  plot_df$seed_marker <- "Other seeds"
  if (length(top3_seeds) >= 1L) plot_df$seed_marker[plot_df$seed == top3_seeds[[1]]] <- paste0("Top 1: ", top3_seeds[[1]], " (*)")
  if (length(top3_seeds) >= 2L) plot_df$seed_marker[plot_df$seed == top3_seeds[[2]]] <- paste0("Top 2: ", top3_seeds[[2]], " (triangle)")
  if (length(top3_seeds) >= 3L) plot_df$seed_marker[plot_df$seed == top3_seeds[[3]]] <- paste0("Top 3: ", top3_seeds[[3]], " (black dot)")
  other_df <- plot_df[plot_df$seed_marker == "Other seeds", , drop = FALSE]
  top_df <- plot_df[plot_df$seed_marker != "Other seeds", , drop = FALSE]
  other_vivo_df <- other_df[as.character(other_df$context) == "in vivo", , drop = FALSE]
  other_vitro_df <- other_df[as.character(other_df$context) == "in vitro", , drop = FALSE]
  top_breaks <- c(
    if (length(top3_seeds) >= 1L) paste0("Top 1: ", top3_seeds[[1]], " (*)"),
    if (length(top3_seeds) >= 2L) paste0("Top 2: ", top3_seeds[[2]], " (triangle)"),
    if (length(top3_seeds) >= 3L) paste0("Top 3: ", top3_seeds[[3]], " (black dot)")
  )
  shape_values <- c(
    if (length(top3_seeds) >= 1L) setNames(8, paste0("Top 1: ", top3_seeds[[1]], " (*)")),
    if (length(top3_seeds) >= 2L) setNames(17, paste0("Top 2: ", top3_seeds[[2]], " (triangle)")),
    if (length(top3_seeds) >= 3L) setNames(16, paste0("Top 3: ", top3_seeds[[3]], " (black dot)"))
  )
  top3_label <- if (length(top3_seeds)) {
    if (isTRUE(pred1000_top3_only)) {
      paste0("Black symbols mark eligible objective top 3: ", paste(top3_seeds, collapse = "; "), ".")
    } else {
      paste0("Black symbols mark objective top 3: ", paste(top3_seeds, collapse = "; "), ".")
    }
  } else if (isTRUE(pred1000_top3_only)) {
    "No seeds met the 2N/4N 1000-day prediction gate."
  } else {
    ""
  }

  vivo_pair_df <- plot_df[
    as.character(plot_df$context) == "in vivo",
    c("seed", "parameter_label", "boundary_x_plot", "y_plot", "seed_marker"),
    drop = FALSE
  ]
  vitro_pair_df <- plot_df[
    as.character(plot_df$context) == "in vitro",
    c("seed", "parameter_label", "boundary_x_plot", "y_plot", "seed_marker"),
    drop = FALSE
  ]
  pair_df <- merge(
    vivo_pair_df,
    vitro_pair_df,
    by = c("seed", "parameter_label"),
    suffixes = c("_vivo", "_vitro"),
    all = FALSE,
    sort = FALSE
  )
  pair_df$is_top_seed <- pair_df$seed %in% top3_seeds
  other_pair_df <- pair_df[!pair_df$is_top_seed, , drop = FALSE]
  top_pair_df <- pair_df[pair_df$is_top_seed, , drop = FALSE]

  title_detail <- if (!is.null(title_suffix) && nzchar(title_suffix)) paste0(" (", title_suffix, ")") else ""
  subtitle_line1 <- if (identical(axis_cfg$axis_type, "relative")) {
    paste0(
      "0 = joint union lower bound, 1 = joint union upper bound; shaded zones are within ",
      sprintf("%.0f", 100 * near_thresh),
      "% of a bound."
    )
  } else {
    "Horizontal lines span natural joint union lower-to-upper bounds."
  }
  subtitle_text <- paste(
    c(
      subtitle_line1,
      "In vivo points are colored by Vxx warm-up family; in vitro points are blue; lines connect paired seed-parameter values.",
      top3_label,
      trimws(axis_cfg$subtitle_note)
    )[nzchar(c(
      subtitle_line1,
      "In vivo points are colored by Vxx warm-up family; in vitro points are blue; lines connect paired seed-parameter values.",
      top3_label,
      trimws(axis_cfg$subtitle_note)
    ))],
    collapse = " "
  )
  subtitle_text <- paste(strwrap(subtitle_text, width = 125L), collapse = "\n")

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = boundary_x_plot, y = y_plot))
  if (identical(axis_cfg$axis_type, "relative")) {
    p <- p +
      ggplot2::annotate("rect", xmin = axis_cfg$lower_rect[[1]], xmax = axis_cfg$lower_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#fddbc7", alpha = 0.28) +
      ggplot2::annotate("rect", xmin = axis_cfg$upper_rect[[1]], xmax = axis_cfg$upper_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#d1e5f0", alpha = 0.28) +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = boundary_x_start, xend = boundary_x_end, y = y_base, yend = y_base),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      ) +
      ggplot2::geom_vline(xintercept = axis_cfg$vlines, color = axis_cfg$vline_colors, linetype = axis_cfg$vline_linetypes, linewidth = 0.35)
  } else {
    p <- p +
      ggplot2::geom_segment(
        data = ref_df,
        ggplot2::aes(x = boundary_x_lower, xend = boundary_x_upper, y = y_base, yend = y_base),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      )
  }
  if (nrow(other_pair_df)) {
    p <- p +
      ggplot2::geom_segment(
        data = other_pair_df,
        ggplot2::aes(
          x = boundary_x_plot_vivo,
          xend = boundary_x_plot_vitro,
          y = y_plot_vivo,
          yend = y_plot_vitro
        ),
        inherit.aes = FALSE,
        color = "grey55",
        alpha = 0.12,
        linewidth = 0.18
      )
  }
  if (nrow(top_pair_df)) {
    p <- p +
      ggplot2::geom_segment(
        data = top_pair_df,
        ggplot2::aes(
          x = boundary_x_plot_vivo,
          xend = boundary_x_plot_vitro,
          y = y_plot_vivo,
          yend = y_plot_vitro
        ),
        inherit.aes = FALSE,
        color = "grey25",
        alpha = 0.70,
        linewidth = 0.45
      )
  }

  p <- p +
    ggplot2::geom_point(
      data = other_vitro_df,
      shape = 16,
      size = 2.1,
      color = "#0072B2",
      alpha = 0.42
    ) +
    ggplot2::geom_point(
      data = other_vivo_df,
      ggplot2::aes(color = invivo_warmup),
      shape = 16,
      size = 2.1,
      alpha = 0.62
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = levels(plot_df$parameter),
      expand = ggplot2::expansion(add = 0.45)
    ) +
    axis_cfg$scale +
    ggplot2::labs(
      title = paste0("Joint Soft-Coupled In Vivo/In Vitro Paired Parameter Positions", title_detail),
      subtitle = subtitle_text,
      x = axis_cfg$x_label,
      y = NULL,
      color = "In vivo warm-up"
    ) +
    scale_color_invivo_warmup(family_levels) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1), order = 1)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )

  if (nrow(top_df)) {
    p <- p +
      ggplot2::geom_point(
        data = top_df,
        ggplot2::aes(shape = seed_marker),
        size = 3.0,
        color = "black"
      ) +
      ggplot2::scale_shape_manual(values = shape_values, breaks = top_breaks, drop = FALSE) +
      ggplot2::labs(shape = if (isTRUE(pred1000_top3_only)) "Top 3 eligible seeds" else "Objective top 3 seeds") +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1), order = 1),
        shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE, order = 2)
      )
  }
  save_plot_pdf_png(p, root, stem, width = 12, height = 6.2)
}

plot_multi_warmup_seed_trajectory <- function(root,
                                              stem,
                                              data_file,
                                              summary_file,
                                              cohort_value,
                                              y_col,
                                              summary_transform = identity,
                                              y_transform = identity,
                                              title,
                                              subtitle,
                                              y_label,
                                              add_ploidy_axis = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  extra_dir <- integrated_extra_results_dir(root)
  seed_df <- read_table_optional(file.path(extra_dir, data_file))
  mean_df <- read_table_optional(file.path(extra_dir, summary_file))
  if (is.null(seed_df) || !is.data.frame(seed_df) || !all(c("seed", "cohort", "day", y_col) %in% names(seed_df))) return(invisible(NULL))
  if (is.null(mean_df) || !is.data.frame(mean_df) || !all(c("cohort", "day", "mean_value", "median_value") %in% names(mean_df))) return(invisible(NULL))
  seed_df <- seed_df[as.character(seed_df$cohort) == cohort_value, , drop = FALSE]
  mean_df <- mean_df[as.character(mean_df$cohort) == cohort_value, , drop = FALSE]
  if (!nrow(seed_df) || !nrow(mean_df)) return(invisible(NULL))
  seed_df <- add_invivo_warmup_from_seed(seed_df, root)
  seed_df$day <- num_col(seed_df$day)
  seed_df$plot_value <- y_transform(num_col(seed_df[[y_col]]))
  mean_df$day <- num_col(mean_df$day)
  mean_df$plot_mean <- summary_transform(num_col(mean_df$mean_value))
  mean_df$plot_median <- summary_transform(num_col(mean_df$median_value))
  seed_df <- seed_df[is.finite(seed_df$day) & is.finite(seed_df$plot_value), , drop = FALSE]
  mean_df <- mean_df[is.finite(mean_df$day) & is.finite(mean_df$plot_mean), , drop = FALSE]
  if (!nrow(seed_df) || !nrow(mean_df)) return(invisible(NULL))
  summary_line_color <- if (identical(cohort_value, "2N")) "#2166AC" else "#B2182B"
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = seed_df,
      ggplot2::aes(x = day, y = plot_value, group = seed, color = invivo_warmup),
      linewidth = 0.28,
      alpha = 0.58
    ) +
    ggplot2::geom_line(
      data = mean_df,
      ggplot2::aes(x = day, y = plot_mean),
      color = summary_line_color,
      linewidth = 1.1
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Day",
      y = y_label,
      color = "In vivo warm-up"
    ) +
    scale_color_invivo_warmup(levels(seed_df$invivo_warmup)) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.84, alpha = 1))) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
  if (isTRUE(add_ploidy_axis)) {
    p <- p +
      ggplot2::geom_line(
        data = mean_df[is.finite(mean_df$plot_median), , drop = FALSE],
        ggplot2::aes(x = day, y = plot_median),
        color = summary_line_color,
        linewidth = 0.75,
        linetype = "dashed"
      ) +
      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ . / 22, name = "Ploidy"))
  }
  save_plot_pdf_png(p, root, stem, width = 12, height = 6)
}

plot_integrated_colored_seed_figures <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_multi_warmup_boundary_forest(
    root,
    stem = "multi_warmup_integrated_parameter_boundary_forest_by_invivo_warmup",
    x_scale = "relative"
  )
  plot_multi_warmup_paired_boundary_scatter(
    root,
    stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_by_invivo_warmup",
    x_scale = "relative"
  )
  plot_multi_warmup_boundary_forest(
    root,
    stem = "multi_warmup_integrated_parameter_boundary_forest_pred1000_gt44_top3_by_invivo_warmup",
    title_suffix = "Pred1000 > 44 Top 3",
    pred1000_top3_only = TRUE,
    x_scale = "relative"
  )
  plot_multi_warmup_paired_boundary_scatter(
    root,
    stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_pred1000_gt44_top3_by_invivo_warmup",
    title_suffix = "Pred1000 > 44 Top 3",
    pred1000_top3_only = TRUE,
    x_scale = "relative"
  )
  plot_multi_warmup_boundary_forest(
    root,
    stem = "multi_warmup_integrated_parameter_boundary_forest_log_x_by_invivo_warmup",
    title_suffix = "Log X Scale",
    x_scale = "log10_original"
  )
  plot_multi_warmup_paired_boundary_scatter(
    root,
    stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_log_x_by_invivo_warmup",
    title_suffix = "Log X Scale",
    x_scale = "log10_original"
  )
  plot_multi_warmup_boundary_forest(
    root,
    stem = "multi_warmup_integrated_parameter_boundary_forest_pred1000_gt44_top3_log_x_by_invivo_warmup",
    title_suffix = "Pred1000 > 44 Top 3, Log X Scale",
    pred1000_top3_only = TRUE,
    x_scale = "log10_original"
  )
  plot_multi_warmup_paired_boundary_scatter(
    root,
    stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_pred1000_gt44_top3_log_x_by_invivo_warmup",
    title_suffix = "Pred1000 > 44 Top 3, Log X Scale",
    pred1000_top3_only = TRUE,
    x_scale = "log10_original"
  )
  plot_multi_warmup_seed_trajectory(
    root,
    stem = "multi_warmup_integrated_burden_total_log_seed_mean_2N_by_invivo_warmup",
    data_file = "predict1000_burden_total_seed_day_mean.tsv",
    summary_file = "predict1000_burden_total_mean_ci.tsv",
    cohort_value = "2N",
    y_col = "burden_value",
    y_transform = function(x) log10(x),
    summary_transform = function(x) log10(x),
    title = "1000-day Total Burden Prediction by In Vivo Warm-Up Family: 2N",
    subtitle = "Colored hairlines = individual seed trajectories by Vxx warm-up source; solid line = cross-seed mean.",
    y_label = "log10 total tumor burden (mm^3)"
  )
  plot_multi_warmup_seed_trajectory(
    root,
    stem = "multi_warmup_integrated_burden_total_log_seed_mean_4N_by_invivo_warmup",
    data_file = "predict1000_burden_total_seed_day_mean.tsv",
    summary_file = "predict1000_burden_total_mean_ci.tsv",
    cohort_value = "4N",
    y_col = "burden_value",
    y_transform = function(x) log10(x),
    summary_transform = function(x) log10(x),
    title = "1000-day Total Burden Prediction by In Vivo Warm-Up Family: 4N",
    subtitle = "Colored hairlines = individual seed trajectories by Vxx warm-up source; solid line = cross-seed mean.",
    y_label = "log10 total tumor burden (mm^3)"
  )
  plot_multi_warmup_seed_trajectory(
    root,
    stem = "multi_warmup_integrated_ploidy_seed_curves_2N_by_invivo_warmup",
    data_file = "predict1000_ploidy_seed_day_mean.tsv",
    summary_file = "predict1000_ploidy_mean_ci.tsv",
    cohort_value = "2N",
    y_col = "ploidy_value",
    title = "1000-day Ploidy Seed Trajectories by In Vivo Warm-Up Family: 2N",
    subtitle = "Colored hairlines = individual seed trajectories by Vxx warm-up source; solid line = cross-seed mean; dashed line = median.",
    y_label = "Mean chromosome count",
    add_ploidy_axis = TRUE
  )
  plot_multi_warmup_seed_trajectory(
    root,
    stem = "multi_warmup_integrated_ploidy_seed_curves_4N_by_invivo_warmup",
    data_file = "predict1000_ploidy_seed_day_mean.tsv",
    summary_file = "predict1000_ploidy_mean_ci.tsv",
    cohort_value = "4N",
    y_col = "ploidy_value",
    title = "1000-day Ploidy Seed Trajectories by In Vivo Warm-Up Family: 4N",
    subtitle = "Colored hairlines = individual seed trajectories by Vxx warm-up source; solid line = cross-seed mean; dashed line = median.",
    y_label = "Mean chromosome count",
    add_ploidy_axis = TRUE
  )
  invisible(NULL)
}

report_figure_section_html <- function(root, section_id, title, legend, filename_stem) {
  png_path <- file.path(root, paste0(filename_stem, ".png"))
  pdf_path <- file.path(root, paste0(filename_stem, ".pdf"))
  if (!file.exists(png_path) && !file.exists(pdf_path)) return("")
  image_data <- if (file.exists(png_path)) {
    file_to_data_uri(png_path, mime = "image/png")
  } else {
    pdf_to_image_data_uri(pdf_path)
  }
  figure_body <- if (!is.null(image_data) && nzchar(image_data)) {
    paste0(
      '<div class="report-figure"><img src="',
      escape_html(image_data),
      '" alt="', escape_html(title), '" class="report-figure-image"/></div>'
    )
  } else {
    '<p class="report-empty">Plot preview rendering is unavailable.</p>'
  }
  paste0(
    '<section class="report-section" id="', escape_html(section_id), '">',
    figure_body,
    '<h2 class="report-figure-title">', title, '</h2>',
    '<p class="report-figure-legend">', escape_html(legend), '</p>',
    '<p class="report-figure-file"><code>', escape_html(paste0(filename_stem, ".pdf")), '</code></p>',
    '</section>'
  )
}

replace_report_section_by_id <- function(fragment, section_id, replacement) {
  if (!nzchar(replacement)) return(fragment)
  pattern <- paste0('<section class="report-section" id="', section_id, '">.*?</section>')
  sub(pattern, replacement, fragment, perl = TRUE)
}

replace_integrated_colored_seed_sections <- function(fragment, root) {
  replacements <- list(
    "integrated-joint-figure-5" = list(
      title = "Figure 2.5. Parameter Boundary Forest",
      legend = "Fitted active-parameter positions within bounds. Seed-level points are colored by Vxx in vivo warm-up family.",
      stem = "multi_warmup_integrated_parameter_boundary_forest_by_invivo_warmup"
    ),
    "integrated-joint-figure-6" = list(
      title = "Figure 2.5b. Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter",
      legend = "Soft-coupled paired in vivo and in vitro context-specific parameter positions. In vivo points are colored by Vxx warm-up family; in vitro points are blue; lines connect paired seed-parameter values.",
      stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_by_invivo_warmup"
    ),
    "integrated-joint-figure-7" = list(
      title = "Figure 2.6. Parameter Boundary Forest (Pred1000 &gt; 44 Top 3)",
      legend = "Top three integrated seeds passing the 2N/4N 1000-day prediction gate, with seed-level points colored by Vxx in vivo warm-up family.",
      stem = "multi_warmup_integrated_parameter_boundary_forest_pred1000_gt44_top3_by_invivo_warmup"
    ),
    "integrated-joint-figure-8" = list(
      title = "Figure 2.6b. Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Pred1000 &gt; 44 Top 3)",
      legend = "Soft-coupled paired context-specific parameter positions for all integrated seeds, with eligible prediction-gated objective top 3 seeds highlighted. In vivo points are colored by Vxx warm-up family; in vitro points are blue.",
      stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_pred1000_gt44_top3_by_invivo_warmup"
    ),
    "integrated-joint-figure-9" = list(
      title = "Figure 2.7. Parameter Boundary Forest (Log X Scale)",
      legend = "Fitted active-parameter positions on the original parameter scale. Seed-level points are colored by Vxx in vivo warm-up family.",
      stem = "multi_warmup_integrated_parameter_boundary_forest_log_x_by_invivo_warmup"
    ),
    "integrated-joint-figure-10" = list(
      title = "Figure 2.7b. Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Log X Scale)",
      legend = "Natural-scale soft-coupled paired in vivo and in vitro context-specific parameter positions on a log10 x-axis. In vivo points are colored by Vxx warm-up family; in vitro points are blue.",
      stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_log_x_by_invivo_warmup"
    ),
    "integrated-joint-figure-11" = list(
      title = "Figure 2.8. Parameter Boundary Forest (Pred1000 &gt; 44 Top 3, Log X Scale)",
      legend = "Top three integrated seeds passing the 2N/4N 1000-day prediction gate on the original parameter scale, with seed-level points colored by Vxx in vivo warm-up family.",
      stem = "multi_warmup_integrated_parameter_boundary_forest_pred1000_gt44_top3_log_x_by_invivo_warmup"
    ),
    "integrated-joint-figure-12" = list(
      title = "Figure 2.8b. Joint Soft-Coupling In Vivo/In Vitro Paired Boundary Scatter (Pred1000 &gt; 44 Top 3, Log X Scale)",
      legend = "Natural-scale soft-coupled paired context-specific parameter positions for all integrated seeds, with eligible prediction-gated objective top 3 seeds highlighted. In vivo points are colored by Vxx warm-up family; in vitro points are blue.",
      stem = "multi_warmup_integrated_joint_soft_context_boundary_forest_pred1000_gt44_top3_log_x_by_invivo_warmup"
    ),
    "integrated-joint-figure-14" = list(
      title = "Figure 2.10. Cross-seed 1000-day Total Burden Prediction: 2N",
      legend = "Seed-level 2N burden trajectories colored by Vxx in vivo warm-up family; the solid line shows the cross-seed mean.",
      stem = "multi_warmup_integrated_burden_total_log_seed_mean_2N_by_invivo_warmup"
    ),
    "integrated-joint-figure-15" = list(
      title = "Figure 2.11. Cross-seed 1000-day Total Burden Prediction: 4N",
      legend = "Seed-level 4N burden trajectories colored by Vxx in vivo warm-up family; the solid line shows the cross-seed mean.",
      stem = "multi_warmup_integrated_burden_total_log_seed_mean_4N_by_invivo_warmup"
    ),
    "integrated-joint-figure-16" = list(
      title = "Figure 2.12. Cross-seed 1000-day Ploidy Seed Trajectories: 2N",
      legend = "Seed-level 2N ploidy trajectories colored by Vxx in vivo warm-up family; solid and dashed lines show cross-seed mean and median.",
      stem = "multi_warmup_integrated_ploidy_seed_curves_2N_by_invivo_warmup"
    ),
    "integrated-joint-figure-17" = list(
      title = "Figure 2.13. Cross-seed 1000-day Ploidy Seed Trajectories: 4N",
      legend = "Seed-level 4N ploidy trajectories colored by Vxx in vivo warm-up family; solid and dashed lines show cross-seed mean and median.",
      stem = "multi_warmup_integrated_ploidy_seed_curves_4N_by_invivo_warmup"
    )
  )
  for (section_id in names(replacements)) {
    item <- replacements[[section_id]]
    fragment <- replace_report_section_by_id(
      fragment,
      section_id,
      report_figure_section_html(root, section_id, item$title, item$legend, item$stem)
    )
  }
  fragment
}

nav_link <- function(href, label, class = "mw-nav-link") {
  paste0('<li><a class="', class, '" href="#', escape_html(href), '">', escape_html(label), '</a></li>')
}

sidebar_html <- function() {
  part1_items <- c(
    nav_link("part-1", "Part 1. Warm-Up Sweep Summary", "mw-nav-link mw-nav-part"),
    nav_link("run-status", "1.1 Run Status"),
    nav_link("umap-profile", "1.2 18D Warm-Start Profile UMAP"),
    nav_link("umap-cluster", "1.3 18D UMAP by In Vivo Cluster"),
    nav_link("objective-summary", "1.4 Objective Summary"),
    nav_link("objective-tradeoff", "1.5 In Vivo vs In Vitro Objective"),
    nav_link("final-ratios", "1.6 Final Parameter Ratios"),
    nav_link("invivo-only-params", "1.7 In Vivo-Only Warm-Start Parameters"),
    nav_link("deoptim-iterations", "1.8 DEoptim Iteration Frequency"),
    nav_link("warmup-manifest", "1.9 Warm-Up Manifest"),
    nav_link("invivo-cluster-summary", "1.10 In Vivo Cluster Summary"),
    nav_link("best-seed-summary", "1.11 Best Seed Summary"),
    nav_link("final-basin-assignments", "1.12 Final Basin Assignments")
  )
  part2_items <- c(
    nav_link("part-2", "Part 2. Integrated Joint Fitting Results", "mw-nav-link mw-nav-part"),
    nav_link("integrated-joint-chapter-1", "2.1 Joint Summary"),
    nav_link("integrated-joint-convergence-summary", "2.1.1 Convergence Summary"),
    nav_link("integrated-joint-parameter-summary", "2.1.2 Parameter Summary"),
    nav_link("integrated-joint-chapter-2", "2.2 In Vivo"),
    nav_link("integrated-joint-chapter-3", "2.3 In Vitro")
  )
  part3_items <- c(
    nav_link("run-provenance", "Part 3. Run Provenance", "mw-nav-link mw-nav-part")
  )
  paste0(
    '<aside class="mw-sidebar" aria-label="Multi-warm-up report navigation">',
    '<div class="mw-sidebar-header">',
    '<div class="mw-kicker">Navigation</div>',
    '<div class="mw-sidebar-title">Multi-Warm-Up Results</div>',
    '</div>',
    '<nav class="mw-nav"><ul>',
    paste(c(part1_items, part2_items, part3_items), collapse = ""),
    '</ul></nav>',
    '</aside>'
  )
}

regex_extract <- function(x, pattern) {
  m <- regexec(pattern, x, perl = TRUE)
  hit <- regmatches(x, m)[[1]]
  if (length(hit) >= 2L) hit[[2]] else ""
}

rewrite_integrated_joint_fragment <- function(fragment, root) {
  if (!nzchar(fragment)) return(fragment)
  fragment <- sub(
    '(?s)<section class="report-hero">.*?</section>',
    paste0(
      '<section class="report-hero" id="part-2">',
      '<h1>Part 2. Integrated Joint Fitting Results</h1>',
      '<p class="report-meta">',
      'This section reuses the standard joint <code>extra_results.R</code> report logic after combining all warm-up pairs into one prefixed-seed run. ',
      'Seed IDs are written as <code>warmup_label__seed</code> so each row and figure can still be traced back to its warm-up pair.<br/>',
      '<strong>Integrated source:</strong> <code>',
      escape_html(file.path(root, "multi_warmup_integrated_joint_run", "extra_results")),
      '</code>',
      '</p></section>'
    ),
    fragment,
    perl = TRUE
  )
  replacements <- c(
    "1. Joint Summary" = "2.1 Joint Summary",
    "2. In Vivo" = "2.2 In Vivo",
    "3. In Vitro" = "2.3 In Vitro"
  )
  for (old in names(replacements)) {
    fragment <- gsub(old, replacements[[old]], fragment, fixed = TRUE)
  }
  fragment <- gsub('id="chapter-([0-9]+)"', 'id="integrated-joint-chapter-\\1"', fragment, perl = TRUE)
  fragment <- gsub('id="figure-([0-9]+)"', 'id="integrated-joint-figure-\\1"', fragment, perl = TRUE)
  fragment <- gsub('id="convergence-summary"', 'id="integrated-joint-convergence-summary"', fragment, fixed = TRUE)
  fragment <- gsub('id="parameter-summary"', 'id="integrated-joint-parameter-summary"', fragment, fixed = TRUE)
  fragment <- gsub(
    '<h2 class="report-figure-title">Convergence Summary</h2>',
    '<h2 class="report-figure-title">2.1.1 Convergence Summary</h2>',
    fragment,
    fixed = TRUE
  )
  fragment <- gsub(
    '<h2 class="report-figure-title">Parameter Summary</h2>',
    '<h2 class="report-figure-title">2.1.2 Parameter Summary</h2>',
    fragment,
    fixed = TRUE
  )
  fragment <- gsub("Figure ([0-9A-Za-z.]+)\\.", "Figure 2.\\1.", fragment, perl = TRUE)
  fragment <- insert_after_first_regex(
    fragment,
    '<section class="report-section" id="integrated-joint-figure-4">.*?</section>',
    integrated_tradeoff_by_invivo_warmup_html(root)
  )
  replace_integrated_colored_seed_sections(fragment, root)
}

integrated_joint_results_fragment <- function(root) {
  report_path <- file.path(root, "multi_warmup_integrated_joint_run", "extra_results", "extra_results_report.html")
  if (!file.exists(report_path)) {
    return(list(
      style = "",
      html = paste0(
        '<section class="card">',
        '<h2>Part 2. Integrated Joint Fitting Results</h2>',
        '<p class="empty">Integrated joint extra_results report was not found. Run ',
        '<code>collect_multi_warmup_results.R</code> first to build ',
        '<code>multi_warmup_integrated_joint_run/extra_results</code>.</p>',
        '</section>'
      )
    ))
  }
  txt <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
  style <- regex_extract(txt, "(?s)<style>(.*?)</style>")
  main <- regex_extract(txt, '(?s)<main class="report-main">(.*)</main>\\s*</div>\\s*</body>\\s*</html>')
  if (!nzchar(main)) {
    return(list(
      style = style,
      html = paste0(
        '<section class="card">',
        '<h2>Part 2. Integrated Joint Fitting Results</h2>',
        '<p class="empty">Integrated joint report exists but its main content could not be extracted: <code>',
        escape_html(report_path),
        '</code>.</p></section>'
      )
    ))
  }
  list(
    style = style,
    html = paste0('<div class="integrated-joint-report report-main">', rewrite_integrated_joint_fragment(main, root), '</div>')
  )
}

read_pair_seed_summaries <- function(root, manifest) {
  if (is.null(manifest) || !is.data.frame(manifest) || !nrow(manifest) || !("joint_run_prefix" %in% names(manifest))) {
    return(data.frame())
  }
  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    warmup_label <- if ("warmup_label" %in% names(manifest)) as.character(manifest$warmup_label[[i]]) else as.character(i)
    run_dir <- file.path(root, as.character(manifest$joint_run_prefix[[i]]))
    seed_summary <- read_table_optional(file.path(run_dir, "extra_results", "seed_summary.tsv"))
    if (is.null(seed_summary) || !is.data.frame(seed_summary) || !nrow(seed_summary)) next
    seed_summary$warmup_label <- warmup_label
    seed_summary$joint_run_prefix <- as.character(manifest$joint_run_prefix[[i]])
    rows[[length(rows) + 1L]] <- seed_summary
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

run_status_table <- function(root, summary_df, manifest) {
  tasks <- read_table_optional(file.path(root, "multi_warmup_tasks.tsv"))
  seed_summary_all <- read_pair_seed_summaries(root, manifest)
  total <- if (is.data.frame(manifest)) nrow(manifest) else 0L
  completed <- if (is.data.frame(summary_df)) nrow(summary_df) else 0L
  planned_seed_tasks <- if (is.data.frame(tasks) && nrow(tasks)) nrow(tasks) else NA_integer_
  planned_per_pair <- if (is.data.frame(tasks) && nrow(tasks) && all(c("warmup_label", "joint_seed") %in% names(tasks))) {
    stats::aggregate(joint_seed ~ warmup_label, data = tasks, FUN = function(x) length(unique(x)))
  } else {
    data.frame()
  }
  completed_per_pair <- if (is.data.frame(seed_summary_all) && nrow(seed_summary_all) && all(c("warmup_label", "seed") %in% names(seed_summary_all))) {
    stats::aggregate(seed ~ warmup_label, data = seed_summary_all, FUN = function(x) length(unique(x)))
  } else {
    data.frame()
  }
  sigma_values <- c(
    if (is.data.frame(tasks) && "joint_soft_coupling_sigma_default" %in% names(tasks)) tasks$joint_soft_coupling_sigma_default else NULL,
    if (is.data.frame(seed_summary_all) && "joint_soft_coupling_sigma_default" %in% names(seed_summary_all)) seed_summary_all$joint_soft_coupling_sigma_default else NULL
  )
  itermax_values <- c(
    if (is.data.frame(tasks) && "itermax" %in% names(tasks)) tasks$itermax else NULL,
    if (is.data.frame(seed_summary_all) && "itermax" %in% names(seed_summary_all)) seed_summary_all$itermax else NULL
  )
  n_soft_params <- if (is.data.frame(seed_summary_all) && "joint_soft_coupling_n_params" %in% names(seed_summary_all)) {
    seed_summary_all$joint_soft_coupling_n_params
  } else {
    NULL
  }
  status <- data.frame(
    metric = c(
      "planned_pairs",
      "completed_pairs",
      "missing_or_failed_pairs",
      "planned_seed_tasks",
      "planned_seeds_per_pair",
      "completed_seed_results",
      "completed_seeds_per_pair",
      "soft_coupling_sigma",
      "soft_coupling_parameter_count",
      "DEoptim_itermax",
      "qos",
      "walltime"
    ),
    value = c(
      total,
      completed,
      max(total - completed, 0L),
      if (is.na(planned_seed_tasks)) "NA" else planned_seed_tasks,
      format_count_summary(if (nrow(planned_per_pair)) planned_per_pair$joint_seed else integer(0)),
      if (is.data.frame(seed_summary_all) && nrow(seed_summary_all)) nrow(seed_summary_all) else 0L,
      format_count_summary(if (nrow(completed_per_pair)) completed_per_pair$seed else integer(0)),
      format_unique_values(sigma_values),
      format_unique_values(n_soft_params),
      format_unique_values(itermax_values),
      if (is.data.frame(tasks) && "qos" %in% names(tasks)) format_unique_values(tasks$qos) else "NA",
      if (is.data.frame(tasks) && "walltime" %in% names(tasks)) format_unique_values(tasks$walltime) else "NA"
    ),
    source = c(
      "multi_warmup_manifest.tsv",
      "multi_warmup_best_seed_summary.tsv",
      "manifest vs completed summary",
      "multi_warmup_tasks.tsv",
      "multi_warmup_tasks.tsv",
      "per-pair extra_results/seed_summary.tsv",
      "per-pair extra_results/seed_summary.tsv",
      "multi_warmup_tasks.tsv / seed_summary.tsv",
      "per-pair extra_results/seed_summary.tsv",
      "multi_warmup_tasks.tsv / seed_summary.tsv",
      "multi_warmup_tasks.tsv",
      "multi_warmup_tasks.tsv"
    ),
    stringsAsFactors = FALSE
  )
  status
}

status_counts <- function(summary_df, manifest) {
  total <- if (is.data.frame(manifest)) nrow(manifest) else 0L
  completed <- if (is.data.frame(summary_df)) nrow(summary_df) else 0L
  data.frame(
    metric = c("planned_pairs", "completed_pairs", "missing_or_failed_pairs"),
    value = c(total, completed, max(total - completed, 0L)),
    stringsAsFactors = FALSE
  )
}

decision_text <- function(summary_df) {
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df)) {
    return("No completed warm-up pair summaries were available at report time.")
  }
  if ("final_basin_id" %in% names(summary_df)) {
    basins <- unique(na.omit(summary_df$final_basin_id))
    if (length(basins) <= 1L) {
      return("All completed warm-up pairs assigned to one final parameter basin; this supports warm-start stability for the completed runs.")
    }
    return(paste0("Completed warm-up pairs assigned to ", length(basins), " final parameter basins; compare objective and biological diagnostics before selecting a final result."))
  }
  "Final basin assignments were not available; use objective and diagnostic summaries for the completed runs."
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root <- normalizePath(as_chr(argv$multi_warmup_root, as_chr(argv$out_dir)), mustWork = TRUE)
  out_path <- normalizePath(as_chr(argv$out, file.path(root, "multi-warm-up_results.html")), mustWork = FALSE)
  manifest <- read_table_optional(file.path(root, "multi_warmup_manifest.tsv"))
  invivo_clusters <- read_table_optional(file.path(root, "multi_warmup_invivo_cluster_summary.tsv"))
  summary_df <- read_table_optional(file.path(root, "multi_warmup_best_seed_summary.tsv"))
  basins <- read_table_optional(file.path(root, "multi_warmup_final_basin_assignments.tsv"))

  manifest_display <- strip_project_path_prefix_table(manifest)
  hidden_summary_cols <- c(
    "seed_dir",
    "cluster_feature_set",
    "invivo_seed_dir",
    "invitro_seed_dir",
    "joint_soft_coupling_parameters_table",
    "joint_run_dir"
  )
  invivo_clusters_display <- drop_report_columns(invivo_clusters, hidden_summary_cols)
  summary_display <- drop_report_columns(summary_df, hidden_summary_cols)
  plot_part1_color_aligned_figures(root)
  plot_integrated_tradeoff_by_invivo_warmup(root)
  plot_integrated_colored_seed_figures(root)
  integrated_joint <- integrated_joint_results_fragment(root)
  provenance_block <- o2sd_run_provenance_html(root, title = "Part 3. Run Provenance", section_id = "run-provenance")

  html <- paste0(
    '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
    '<title>Multi-Warm-Up Results</title>',
    '<style>',
    'body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;background:#f4f7fa;color:#1b2a38;}',
    '.mw-layout{display:flex;gap:24px;max-width:1780px;margin:0 auto;padding:24px;}',
    '.shell{flex:1;min-width:0;max-width:none;padding:0;}',
    '.mw-sidebar{position:sticky;top:24px;align-self:flex-start;width:305px;max-height:calc(100vh - 48px);overflow:auto;border:1px solid #d6dde6;border-radius:10px;background:#fff;box-shadow:0 8px 22px rgba(0,0,0,.05);}',
    '.mw-sidebar-header{padding:14px 16px;background:#1f3348;color:#fff;}',
    '.mw-kicker{font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.08em;opacity:.8;}',
    '.mw-sidebar-title{margin-top:4px;font-size:17px;font-weight:700;line-height:1.2;}',
    '.mw-nav{padding:10px;}.mw-nav ul{margin:0;padding:0;list-style:none;}.mw-nav li{margin:3px 0;}',
    '.mw-nav-link{display:block;border-radius:7px;padding:6px 9px;color:#17324c;text-decoration:none;font-size:12px;font-weight:600;line-height:1.3;}',
    '.mw-nav-part{margin-top:6px;background:#eef3f8;font-size:13px;font-weight:800;}',
    '.mw-nav-link:hover{background:rgba(47,110,164,.08);}',
    '.hero,.card{background:#fff;border:1px solid #d6dde6;border-radius:10px;box-shadow:0 8px 22px rgba(0,0,0,.05);padding:16px;margin-bottom:18px;}',
    '.report-section{background:#fff;border:1px solid #d6dde6;border-radius:10px;box-shadow:0 8px 22px rgba(0,0,0,.05);padding:16px;margin-bottom:18px;}',
    'h1{margin:0 0 8px;font-size:28px;}h2{margin:0 0 8px;font-size:18px;}p{color:#425365;line-height:1.45;}',
    '.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px;}',
    '.part-card h2{font-size:20px;}.part-card p{margin-bottom:0;}',
    '.integrated-joint-report{margin-top:18px;}',
    '.integrated-joint-report .report-hero h1{font-size:25px;}',
    '.figure-image{display:block;width:100%;height:auto;border:1px solid #d7dee7;border-radius:8px;background:#fff;}',
    '.report-table{width:100%;border-collapse:collapse;font-size:12px;}.report-table th,.report-table td{padding:7px 9px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}',
    '.report-table th{background:#f7f9fb;font-weight:700;}.empty{font-style:italic;color:#657789;}code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;}',
    '@media(max-width:1100px){.mw-layout{display:block;padding:14px}.mw-sidebar{position:relative;top:auto;width:auto;max-height:none;margin-bottom:16px}.grid{grid-template-columns:1fr}}',
    integrated_joint$style,
    '</style></head><body><div class="mw-layout">',
    sidebar_html(),
    '<main class="shell">',
    '<section class="hero"><h1>Multi-Warm-Up Results</h1>',
    '<p><strong>Root:</strong> <code>', escape_html(root), '</code><br/>',
    '<strong>Generated:</strong> ', escape_html(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), '<br/>',
    '<strong>Decision summary:</strong> ', escape_html(decision_text(summary_df)), '</p></section>',
    '<section class="card part-card" id="part-1"><h2>Part 1. Warm-Up Sweep Summary</h2>',
    '<p>Cross-pair warm-start selection, run status, objective summaries, paired-parameter ratios, and final basin diagnostics.</p></section>',
    '<section class="card" id="run-status"><h2>1.1 Run Status</h2>', table_html(run_status_table(root, summary_df, manifest), max_rows = 30), '</section>',
    '<div class="grid">',
    figure_html(root, "joint_soft_coupling_18d_profile_umap_500seed.pdf", "1.2 18D Warm-Start Profile UMAP", "Standalone paired 18D warm-start profile UMAP. In vivo rank is shown by color and in vitro rank by point shape.", section_id = "umap-profile"),
    figure_html(root, "joint_soft_coupling_18d_profile_umap_by_invivo_cluster_500seed.pdf", "1.3 18D Warm-Start Profile UMAP by In Vivo Cluster", "Same 18D UMAP with in vivo cluster overlay. Red labels mark selected representative in vivo ranks.", section_id = "umap-cluster"),
    '</div>',
    '<div class="grid">',
    figure_html(root, "multi_warmup_objective_summary.pdf", "1.4 Objective Summary", "Best total objective for each completed warm-up pair.", section_id = "objective-summary"),
    figure_html(root, "multi_warmup_invivo_invitro_objective_scatter.pdf", "1.5 In Vivo vs In Vitro Objective", "Best-seed in vivo and in vitro objective components by warm-up pair.", section_id = "objective-tradeoff"),
    figure_html(root, "multi_warmup_final_parameter_ratio_heatmap.pdf", "1.6 Final Parameter Ratios", "Best-seed final in vivo/in vitro parameter ratios.", section_id = "final-ratios"),
    figure_html(root, "multi_warmup_invivo_only_parameter_heatmap.pdf", "1.7 In Vivo-Only Warm-Start Parameters", "Source in vivo seed values for the four parameters not represented in paired soft-coupling ratios.", section_id = "invivo-only-params"),
    '</div>',
    figure_html(root, "multi_warmup_deoptim_iteration_histograms.pdf", "1.8 DEoptim Iteration Frequency by Warm-Up Pair", "Facet histogram of seed-level DEoptim iterations completed for each warm-up pair.", section_id = "deoptim-iterations"),
    '<section class="card" id="warmup-manifest"><h2>1.9 Warm-Up Manifest</h2>', table_html(manifest_display, max_rows = 200), '</section>',
    '<section class="card" id="invivo-cluster-summary"><h2>1.10 In Vivo Cluster Summary</h2>', table_html(invivo_clusters_display, max_rows = 200), '</section>',
    '<section class="card" id="best-seed-summary"><h2>1.11 Best Seed Summary</h2>', table_html(summary_display, max_rows = 200), '</section>',
    '<section class="card" id="final-basin-assignments"><h2>1.12 Final Basin Assignments</h2>', table_html(basins, max_rows = 200), '</section>',
    integrated_joint$html,
    provenance_block,
    '</main></div></body></html>'
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_path, useBytes = TRUE)
  message("Wrote multi-warm-up report: ", out_path)
}

if (identical(environment(), globalenv())) {
  main()
}
