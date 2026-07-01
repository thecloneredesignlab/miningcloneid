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
source(file.path(UTIL_DIR, "table_io_utils.R"))

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript plot_objective_by_curve_class.R [options]",
      "",
      "Purpose:",
      "  Plot objective distributions by dense-grid curve class with significant pairwise differences.",
      "",
      "Inputs:",
      "  --result_dir=DIR        Dense-grid monotonicity result directory.",
      "  --class_table=FILE      Optional fixed_o2_ploidy_monotonicity_regression_by_seed.tsv.",
      "  --fit_root=DIR          Seed-level fit directory used to read final objective from fit_summary.tsv.",
      "  --class_col=NAME        Class column. Default: curve_class.",
      "  --objective_col=NAME    Optional explicit class-table objective column. Default: final fit_summary objective.",
      "",
      "Outputs:",
      "  --figure_stub=NAME      Output figure/table stem. Default: fixed_o2_regression_objective_by_curve_class_boxplot.",
      "  --alpha=NUM            Adjusted p-value significance threshold. Default: 0.05.",
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

default_result_dir <- function(repo_root = bpf_repo_root()) {
  file.path(
    bpf_dense_grid_result_root(repo_root),
    "dense-grid_monotonicity_classification"
  )
}

default_fit_root <- function(repo_root = bpf_repo_root()) {
  file.path(repo_root, "oxygen", "results", "fit_invivo_O2_buffering_500seed")
}

resolve_class_table <- function(result_dir, class_table = NULL) {
  if (!is.null(class_table) && length(class_table) && !is.na(class_table[[1L]]) && nzchar(as.character(class_table[[1L]]))) {
    return(normalizePath(path.expand(as.character(class_table[[1L]])), mustWork = TRUE))
  }
  candidates <- c(
    file.path(result_dir, "tables", "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv"),
    file.path(result_dir, "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("No supported by-seed class table found under: ", result_dir, call. = FALSE)
  }
  hit[[1L]]
}

read_fit_summary_metric <- function(path, metrics) {
  if (!file.exists(path)) return(list(value = NA_real_, source = NA_character_))
  tab <- tryCatch(bpf_read_tsv(path), error = function(e) data.frame())
  if (!all(c("metric", "value") %in% names(tab))) return(list(value = NA_real_, source = NA_character_))
  for (metric in metrics) {
    idx <- which(tab$metric == metric)
    if (!length(idx)) next
    val <- suppressWarnings(as.numeric(tab$value[[idx[[1L]]]]))
    if (is.finite(val)) return(list(value = val, source = metric))
  }
  list(value = NA_real_, source = NA_character_)
}

read_final_objective_by_seed <- function(fit_root, seed_numbers) {
  seed_numbers <- sort(unique(seed_numbers[is.finite(seed_numbers)]))
  rows <- lapply(seed_numbers, function(seed) {
    seed_dir <- file.path(fit_root, paste0("seed", as.integer(seed)))
    hit <- read_fit_summary_metric(
      file.path(seed_dir, "fit_summary.tsv"),
      # `objective` is the post-local-search final objective in the fitting summary.
      metrics = c("objective", "optimizer_local_objective", "objective_total")
    )
    data.frame(
      seed_number = as.integer(seed),
      final_objective = hit$value,
      final_objective_source = hit$source,
      fit_summary = file.path(seed_dir, "fit_summary.tsv"),
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

select_objective_column <- function(tab, requested = NULL) {
  if (!is.null(requested) && length(requested) && !is.na(requested[[1L]]) && nzchar(as.character(requested[[1L]]))) {
    requested <- as.character(requested[[1L]])
    if (!requested %in% names(tab)) stop("Requested objective column is missing: ", requested, call. = FALSE)
    val <- suppressWarnings(as.numeric(tab[[requested]]))
    if (!any(is.finite(val))) stop("Requested objective column has no finite values: ", requested, call. = FALSE)
    return(requested)
  }
  candidates <- c("final_objective", "objective_total", "objective", "objective_data", "objective_burden", "objective_ploidy")
  for (nm in candidates) {
    if (!nm %in% names(tab)) next
    val <- suppressWarnings(as.numeric(tab[[nm]]))
    if (any(is.finite(val))) return(nm)
  }
  stop("No finite objective column found. Checked: ", paste(candidates, collapse = ", "), call. = FALSE)
}

summarize_objective <- function(plot_df) {
  do.call(rbind, lapply(split(plot_df, plot_df$curve_class), function(d) {
    q <- stats::quantile(d$objective_value, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    data.frame(
      curve_class = d$curve_class[[1L]],
      n_seed = nrow(d),
      objective_min = min(d$objective_value, na.rm = TRUE),
      objective_q25 = q[[1L]],
      objective_median = q[[2L]],
      objective_q75 = q[[3L]],
      objective_max = max(d$objective_value, na.rm = TRUE),
      objective_mean = mean(d$objective_value, na.rm = TRUE),
      objective_sd = stats::sd(d$objective_value, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

p_label <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 5e-2) return("*")
  "ns"
}

pairwise_tests <- function(plot_df, alpha = 0.05) {
  classes <- levels(plot_df$curve_class)
  pairs <- utils::combn(classes, 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    x <- plot_df$objective_value[plot_df$curve_class == pair[[1L]]]
    y <- plot_df$objective_value[plot_df$curve_class == pair[[2L]]]
    p <- tryCatch(
      suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value),
      error = function(e) NA_real_
    )
    data.frame(
      class1 = pair[[1L]],
      class2 = pair[[2L]],
      n1 = length(x),
      n2 = length(y),
      median1 = stats::median(x, na.rm = TRUE),
      median2 = stats::median(y, na.rm = TRUE),
      median_difference_class2_minus_class1 = stats::median(y, na.rm = TRUE) - stats::median(x, na.rm = TRUE),
      p_value = p,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out$significance <- vapply(out$p_adj, p_label, character(1))
  out$significant <- is.finite(out$p_adj) & out$p_adj < alpha
  out[order(out$p_adj, out$class1, out$class2), , drop = FALSE]
}

global_test <- function(plot_df) {
  kt <- tryCatch(
    stats::kruskal.test(objective_value ~ curve_class, data = plot_df),
    error = function(e) NULL
  )
  if (is.null(kt)) {
    return(data.frame(test = "Kruskal-Wallis", statistic = NA_real_, df = NA_real_, p_value = NA_real_))
  }
  data.frame(
    test = "Kruskal-Wallis",
    statistic = unname(kt$statistic),
    df = unname(kt$parameter),
    p_value = kt$p.value,
    stringsAsFactors = FALSE
  )
}

build_bracket_data <- function(pair_tab, plot_df) {
  sig <- pair_tab[pair_tab$significant, , drop = FALSE]
  if (!nrow(sig)) return(data.frame())
  class_levels <- levels(plot_df$curve_class)
  y_min <- min(plot_df$objective_value, na.rm = TRUE)
  y_max <- max(plot_df$objective_value, na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range <= 0) y_range <- max(abs(y_max), 1)
  sig$x1 <- match(sig$class1, class_levels)
  sig$x2 <- match(sig$class2, class_levels)
  sig <- sig[order(abs(sig$x2 - sig$x1), sig$p_adj, sig$x1, sig$x2), , drop = FALSE]
  step <- 0.075 * y_range
  sig$y <- y_max + step * seq_len(nrow(sig))
  sig$tip_y <- sig$y - 0.025 * y_range
  sig$label <- sig$significance
  sig
}

plot_objective_by_class <- function(plot_df, pair_tab, objective_col, out_stem) {
  require_package("ggplot2")
  require_package("viridisLite")
  require_package("scales")
  bracket_df <- build_bracket_data(pair_tab, plot_df)
  y_min <- min(plot_df$objective_value, na.rm = TRUE)
  y_max <- max(plot_df$objective_value, na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range <= 0) y_range <- max(abs(y_max), 1)
  y_upper <- if (nrow(bracket_df)) max(bracket_df$y, na.rm = TRUE) + 0.08 * y_range else y_max + 0.08 * y_range
  class_counts <- as.data.frame(table(plot_df$curve_class), stringsAsFactors = FALSE)
  names(class_counts) <- c("curve_class", "n_seed")
  class_labels <- stats::setNames(
    paste0(gsub("_", "\n", class_counts$curve_class), "\n(n=", class_counts$n_seed, ")"),
    class_counts$curve_class
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = curve_class, y = objective_value)) +
    ggplot2::geom_boxplot(width = 0.58, outlier.shape = NA, fill = "white", color = "#333333", linewidth = 0.45) +
    ggplot2::geom_jitter(
      ggplot2::aes(fill = objective_value),
      shape = 21,
      color = "#333333",
      stroke = 0.22,
      width = 0.18,
      size = 2.1,
      alpha = 0.88
    ) +
    ggplot2::scale_fill_gradientn(
      colors = viridisLite::viridis(256, option = "viridis"),
      name = objective_col,
      labels = scales::comma
    ) +
    ggplot2::scale_x_discrete(labels = class_labels) +
    ggplot2::scale_y_continuous(labels = scales::comma, limits = c(y_min - 0.03 * y_range, y_upper)) +
    ggplot2::labs(
      x = "Regression-smoothed curve class",
      y = objective_col,
      title = "Objective distribution by dense-grid curve class"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8.5, angle = 0, vjust = 0.5),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold")
    )

  if (nrow(bracket_df)) {
    p <- p +
      ggplot2::geom_segment(
        data = bracket_df,
        ggplot2::aes(x = x1, xend = x2, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = bracket_df,
        ggplot2::aes(x = x1, xend = x1, y = tip_y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = bracket_df,
        ggplot2::aes(x = x2, xend = x2, y = tip_y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      ggplot2::geom_text(
        data = bracket_df,
        ggplot2::aes(x = (x1 + x2) / 2, y = y + 0.012 * y_range, label = label),
        inherit.aes = FALSE,
        size = 3.0,
        fontface = "bold"
      )
  }

  dir.create(dirname(out_stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(out_stem, ".pdf"), p, width = 11.5, height = 8.2, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(out_stem, ".png"), p, width = 11.5, height = 8.2, units = "in", dpi = 220)
  invisible(p)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  argv <- bpf_parse_args(args)
  if (bpf_as_bool(argv$help, FALSE)) {
    usage()
    return(invisible(0L))
  }
  repo_root <- bpf_repo_root()
  result_dir <- bpf_resolve_repo_path(argv$result_dir %||% default_result_dir(repo_root), repo_root = repo_root, mustWork = TRUE)
  fit_root <- bpf_resolve_repo_path(argv$fit_root %||% default_fit_root(repo_root), repo_root = repo_root, mustWork = TRUE)
  class_table <- resolve_class_table(result_dir, argv$class_table)
  class_col <- argv$class_col %||% "curve_class"
  alpha <- bpf_as_num(argv$alpha, 0.05)
  figure_stub <- argv$figure_stub %||% "fixed_o2_regression_objective_by_curve_class_boxplot"

  tab <- bpf_read_tsv(class_table)
  if (!class_col %in% names(tab)) stop("Class table is missing --class_col=", class_col, call. = FALSE)
  if (!"seed_number" %in% names(tab)) {
    if (!"seed_id" %in% names(tab)) stop("Class table must contain seed_number or seed_id.", call. = FALSE)
    tab$seed_number <- suppressWarnings(as.integer(sub("^seed", "", as.character(tab$seed_id))))
  }
  final_objectives <- read_final_objective_by_seed(fit_root, suppressWarnings(as.integer(tab$seed_number)))
  tab <- merge(tab, final_objectives, by = "seed_number", all.x = TRUE, sort = FALSE)
  objective_col <- select_objective_column(tab, argv$objective_col)
  objective_value <- suppressWarnings(as.numeric(tab[[objective_col]]))
  plot_df <- data.frame(
    seed_id = if ("seed_id" %in% names(tab)) tab$seed_id else NA_character_,
    seed_number = if ("seed_number" %in% names(tab)) suppressWarnings(as.integer(tab$seed_number)) else NA_integer_,
    curve_class = trimws(as.character(tab[[class_col]])),
    objective_value = objective_value,
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[is.finite(plot_df$objective_value) & nzchar(plot_df$curve_class), , drop = FALSE]
  if (!nrow(plot_df)) stop("No finite objective/class rows available for plotting.", call. = FALSE)
  if (length(unique(plot_df$curve_class)) < 2L) stop("Need at least two curve classes for comparison.", call. = FALSE)

  summary_tab <- summarize_objective(plot_df)
  class_order <- summary_tab$curve_class[order(summary_tab$objective_median, summary_tab$curve_class)]
  plot_df$curve_class <- factor(plot_df$curve_class, levels = class_order)
  summary_tab <- summary_tab[match(class_order, summary_tab$curve_class), , drop = FALSE]
  summary_tab$curve_class <- as.character(summary_tab$curve_class)

  pair_tab <- pairwise_tests(plot_df, alpha = alpha)
  global_tab <- global_test(plot_df)
  run_args <- data.frame(
    argument = c("result_dir", "class_table", "fit_root", "class_col", "objective_col", "n_seed", "n_class", "alpha", "test", "p_adjust_method"),
    value = c(result_dir, class_table, fit_root, class_col, objective_col, as.character(nrow(plot_df)), as.character(length(class_order)), as.character(alpha), "Wilcoxon rank-sum pairwise; Kruskal-Wallis global", "BH"),
    stringsAsFactors = FALSE
  )

  table_dir <- file.path(result_dir, "tables")
  figure_dir <- file.path(result_dir, "figures")
  stem <- file.path(figure_dir, figure_stub)
  bpf_write_tsv(summary_tab, file.path(table_dir, paste0(figure_stub, "_summary.tsv")))
  bpf_write_tsv(pair_tab, file.path(table_dir, paste0(figure_stub, "_pairwise_tests.tsv")))
  bpf_write_tsv(global_tab, file.path(table_dir, paste0(figure_stub, "_global_test.tsv")))
  bpf_write_tsv(
    plot_df[, c("seed_id", "seed_number", "curve_class", "objective_value"), drop = FALSE],
    file.path(table_dir, paste0(figure_stub, "_plot_data.tsv"))
  )
  bpf_write_tsv(run_args, file.path(table_dir, paste0(figure_stub, "_run_arguments.tsv")))
  plot_objective_by_class(plot_df, pair_tab, objective_col, stem)

  message("Objective column: ", objective_col)
  message("Rows plotted: ", nrow(plot_df))
  message("Significant pairwise comparisons: ", sum(pair_tab$significant, na.rm = TRUE))
  message("Wrote figure: ", paste0(stem, ".pdf"))
  invisible(list(summary = summary_tab, pairwise = pair_tab, global = global_tab))
}

if (identical(environment(), globalenv())) {
  main()
}
