#!/usr/bin/env Rscript

# Visualization-only helpers shared by the joint coupling figure scripts.

o2jcv_read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing visualization input: ", path, call. = FALSE)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

o2jcv_theme <- function(base_size = 10, rotate_x = TRUE) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "#E3E7EB", linewidth = 0.25),
      strip.background = ggplot2::element_rect(fill = "#F6F7F8", colour = "#D5DADF", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold", colour = "#25313B"),
      axis.text.x = ggplot2::element_text(angle = if (rotate_x) 35 else 0, hjust = if (rotate_x) 1 else 0.5, colour = "#3B4650"),
      axis.text.y = ggplot2::element_text(colour = "#3B4650"),
      plot.title = ggplot2::element_text(face = "bold", colour = "#1F2933", size = ggplot2::rel(1.15)),
      plot.subtitle = ggplot2::element_text(colour = "#56616B", margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(colour = "#6B747C", hjust = 0),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    )
}

o2jcv_class_colors <- function() {
  c(ClassA = "#2F6B9A", ClassB = "#D4A72C", ClassC = "#D56A27")
}

o2jcv_category_colors <- function() {
  c(CatA = "#2A9D8F", CatB = "#E07A3F", CatC = "#7A5AA6", CatU = "#8A9299")
}

o2jcv_parameter_levels <- function() {
  c(
    "O2_crit", "n_O", "alpha_o2", "gamma_growth", "lam_max", "mu_hp",
    "gamma_mu", "p_mis_base", "p_misseg", "k_o_mis", "p_wgd",
    "buffer_smax", "buffer_beta", "buffer_n_exp"
  )
}

o2jcv_pair_label <- function(pair_id) {
  pair_id <- as.character(pair_id)
  pattern <- "^fit_joint_.+_vi_seed([0-9]+)_C([0-9]+)Sc([0-9]+)_vt_seed([0-9]+)$"
  match <- regexec(pattern, pair_id)
  values <- regmatches(pair_id, match)
  vapply(seq_along(pair_id), function(i) {
    hit <- values[[i]]
    if (length(hit) != 5L) pair_id[[i]] else sprintf("C%02dSc%02d / vi%s", as.integer(hit[[3L]]), as.integer(hit[[4L]]), hit[[2L]])
  }, character(1L))
}

o2jcv_add_pair_label <- function(data) {
  if ("pair_id" %in% names(data)) data$pair_label <- o2jcv_pair_label(data$pair_id)
  data
}

o2jcv_classification_spec <- function(config) {
  value <- function(key, default = NA_character_) {
    hit <- config$value[config$key == key]
    if (length(hit) && nzchar(hit[[1L]])) hit[[1L]] else default
  }
  threshold <- suppressWarnings(as.numeric(value("class_threshold", "1.1")))
  lower <- suppressWarnings(as.numeric(value("class_lower_bound", as.character(1 / threshold))))
  upper <- suppressWarnings(as.numeric(value("class_upper_bound", as.character(threshold))))
  boundary <- value("class_boundary_rule", "classb_inclusive")
  if (!is.finite(lower) || !is.finite(upper) || lower <= 0 || lower >= upper) {
    stop("Invalid ClassB bounds in analysis_config.tsv", call. = FALSE)
  }
  list(lower = lower, upper = upper, boundary_rule = boundary)
}

o2jcv_classb_label <- function(spec, digits = 6L) {
  lower <- format(signif(spec$lower, digits), trim = TRUE, scientific = FALSE)
  upper <- format(signif(spec$upper, digits), trim = TRUE, scientific = FALSE)
  if (identical(spec$boundary_rule, "outer_inclusive")) {
    paste0(lower, " < in vivo / in vitro < ", upper)
  } else {
    paste0(lower, " ≤ in vivo / in vitro ≤ ", upper)
  }
}

o2jcv_class_band <- function(lower_bound, upper_bound = NULL) {
  if (is.null(upper_bound)) {
    upper_bound <- lower_bound
    lower_bound <- 1 / lower_bound
  }
  ggplot2::annotate(
    "rect", xmin = log2(lower_bound), xmax = log2(upper_bound), ymin = -Inf, ymax = Inf,
    fill = o2jcv_class_colors()[["ClassB"]], alpha = 0.12
  )
}

o2jcv_save_figure <- function(figure, stem, out_dir, input_tables, width = 10, height = 6, dpi = 300, family = "unspecified", question = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, figure, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_path, figure, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  manifest_path <- file.path(out_dir, "visualization_manifest.tsv")
  rows <- data.frame(
    figure_id = stem,
    format = c("png", "pdf"),
    path = normalizePath(c(png_path, pdf_path), mustWork = TRUE),
    input_tables = paste(normalizePath(input_tables, mustWork = TRUE), collapse = ";"),
    chart_family = family,
    analytical_question = question,
    width_inches = width,
    height_inches = height,
    stringsAsFactors = FALSE
  )
  if (file.exists(manifest_path)) {
    old <- utils::read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
    old <- old[old$figure_id != stem, , drop = FALSE]
    rows <- rbind(old, rows)
  }
  utils::write.table(rows, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(c(png_path, pdf_path))
}
