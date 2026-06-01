#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) dirname(frame_files[[length(frame_files)]]) else getwd()
})
SCRIPT_DIR <- normalizePath(.script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
rm(.script_dir)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

TARGET_PARAMETERS <- c("p_misseg", "O2_crit", "alpha_o2", "mu_hp")
CONTEXT_LEVELS <- c("in vivo", "in vitro")
OBJECTIVE_TYPES <- c("objective", "objective_invivo", "objective_invitro")
OUTPUT_PREFIX <- "joint_sigma_soft_coupled_paired_seed_comparison"

parse_args <- function(args) {
  out <- list()
  if (!length(args)) return(out)
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  out <- suppressWarnings(as.numeric(x[[1]]))
  if (length(out) && is.finite(out)) out[[1]] else default
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  tolower(trimws(as.character(x[[1]]))) %in% c("true", "t", "1", "yes", "y", "on")
}

html_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

seed_order_key <- function(x) {
  x <- basename(as.character(x))
  m <- regexec("^seed([0-9]+)$", x)
  hit <- regmatches(x, m)
  out <- rep(Inf, length(x))
  for (i in seq_along(hit)) {
    if (length(hit[[i]]) == 2L) out[[i]] <- suppressWarnings(as.numeric(hit[[i]][[2]]))
  }
  out
}

format_sigma <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || !is.finite(x[[1]])) return("unknown")
  x <- x[[1]]
  if (abs(x) >= 1000 || (abs(x) > 0 && abs(x) < 0.001)) {
    format(x, scientific = TRUE, digits = 4, trim = TRUE)
  } else {
    format(x, scientific = FALSE, digits = 4, trim = TRUE)
  }
}

safe_file_token <- function(x) {
  x <- gsub("\\+", "plus", as.character(x))
  x <- gsub("-", "minus", x)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  gsub("_+", "_", x)
}

read_metric_map_optional <- function(path) {
  if (!file.exists(path)) return(character(0))
  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || !all(c("metric", "value") %in% names(tab))) return(character(0))
  setNames(tab$value, tab$metric)
}

first_existing_col <- function(tab, candidates) {
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit)) hit[[1]] else NA_character_
}

resolve_default_results_dir <- function() {
  normalizePath(file.path(WORKFLOW_ROOT, "..", "..", "results"), mustWork = FALSE)
}

resolve_run_dirs <- function(results_dir, run_dirs_arg = NULL) {
  if (!is.null(run_dirs_arg) && nzchar(trimws(run_dirs_arg))) {
    raw <- trimws(unlist(strsplit(run_dirs_arg, ",", fixed = TRUE)))
    raw <- raw[nzchar(raw)]
    dirs <- vapply(raw, function(x) {
      candidate <- x
      if (!dir.exists(candidate)) candidate <- file.path(results_dir, x)
      normalizePath(candidate, mustWork = TRUE)
    }, character(1))
    return(unname(dirs))
  }

  dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
  keep <- vapply(dirs, function(d) {
    seed_dirs <- list.dirs(d, recursive = FALSE, full.names = TRUE)
    seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
    length(seed_dirs) > 0L && any(file.exists(file.path(seed_dirs, "joint_soft_coupling.tsv")))
  }, logical(1))
  dirs <- dirs[keep]
  if (!length(dirs)) {
    stop("No run directories with seed*/joint_soft_coupling.tsv found under: ", results_dir, call. = FALSE)
  }
  normalizePath(dirs[order(basename(dirs))], mustWork = TRUE)
}

read_seed_soft_table <- function(run_dir, seed_dir, target_parameters = TARGET_PARAMETERS) {
  soft_path <- file.path(seed_dir, "joint_soft_coupling.tsv")
  if (!file.exists(soft_path)) return(NULL)
  soft_tab <- tryCatch(
    utils::read.delim(soft_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(soft_tab) || !nrow(soft_tab) || !"parameter" %in% names(soft_tab)) return(NULL)
  soft_tab <- soft_tab[as.character(soft_tab$parameter) %in% target_parameters, , drop = FALSE]
  if (!nrow(soft_tab)) return(NULL)

  seed <- basename(seed_dir)
  fit_summary <- read_metric_map_optional(file.path(seed_dir, "fit_summary.tsv"))
  run_sigma <- as_num(fit_summary[["joint_soft_coupling_sigma_default"]], NA_real_)
  if (!is.finite(run_sigma) && "regularization_sigma" %in% names(soft_tab)) {
    run_sigma <- as_num(soft_tab$regularization_sigma[is.finite(suppressWarnings(as.numeric(soft_tab$regularization_sigma)))][[1]], NA_real_)
  }

  get_num <- function(col) {
    if (col %in% names(soft_tab)) suppressWarnings(as.numeric(soft_tab[[col]])) else rep(NA_real_, nrow(soft_tab))
  }
  get_chr <- function(col) {
    if (col %in% names(soft_tab)) as.character(soft_tab[[col]]) else rep(NA_character_, nrow(soft_tab))
  }
  get_bool <- function(col) {
    if (!(col %in% names(soft_tab))) return(rep(FALSE, nrow(soft_tab)))
    vapply(soft_tab[[col]], as_bool, logical(1), default = FALSE)
  }

  vivo_natural <- get_num("vivo_natural")
  vitro_natural <- get_num("vitro_natural")
  ratio <- get_num("ratio_vivo_to_vitro")
  missing_ratio <- !is.finite(ratio) & is.finite(vivo_natural) & is.finite(vitro_natural) & vitro_natural != 0
  ratio[missing_ratio] <- vivo_natural[missing_ratio] / vitro_natural[missing_ratio]

  data.frame(
    run_label = basename(run_dir),
    run_dir = normalizePath(run_dir, mustWork = FALSE),
    seed = seed,
    seed_id = seed_order_key(seed),
    parameter = as.character(soft_tab$parameter),
    parameter_order = match(as.character(soft_tab$parameter), target_parameters),
    sigma = run_sigma,
    sigma_label = paste0("sigma=", format_sigma(run_sigma)),
    regularization_sigma = get_num("regularization_sigma"),
    split_enabled = get_bool("split_enabled"),
    center_name = get_chr("center_name"),
    delta_name = get_chr("delta_name"),
    center_transformed = get_num("center_transformed"),
    delta_transformed = get_num("delta_transformed"),
    vivo_transformed = get_num("vivo_transformed"),
    vitro_transformed = get_num("vitro_transformed"),
    center_natural = get_num("center_natural"),
    vivo_natural = vivo_natural,
    vitro_natural = vitro_natural,
    delta_interpretable = get_num("delta_interpretable"),
    natural_difference_vivo_to_vitro = get_num("natural_difference_vivo_to_vitro"),
    transformed_difference_vivo_to_vitro = get_num("transformed_difference_vivo_to_vitro"),
    log10_ratio_vivo_to_vitro = get_num("log10_ratio_vivo_to_vitro"),
    fold_change_vivo_to_vitro = get_num("fold_change_vivo_to_vitro"),
    ratio_vivo_to_vitro = ratio,
    odds_ratio_vivo_to_vitro = get_num("odds_ratio_vivo_to_vitro"),
    penalty_paid = get_num("penalty_paid"),
    center_lower_bound = get_num("center_lower_bound"),
    center_upper_bound = get_num("center_upper_bound"),
    center_lower_transformed = get_num("center_lower_transformed"),
    center_upper_transformed = get_num("center_upper_transformed"),
    vivo_clipped = get_bool("vivo_clipped"),
    vitro_clipped = get_bool("vitro_clipped"),
    any_clipped = get_bool("vivo_clipped") | get_bool("vitro_clipped"),
    boundary_status_vivo = get_chr("boundary_status_vivo"),
    boundary_status_vitro = get_chr("boundary_status_vitro"),
    objective = as_num(fit_summary[["objective"]], NA_real_),
    objective_invivo = as_num(fit_summary[["objective_invivo"]], NA_real_),
    objective_invitro = as_num(fit_summary[["objective_invitro"]], NA_real_),
    objective_soft_coupling = as_num(fit_summary[["objective_soft_coupling"]], NA_real_),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_all_soft_tables <- function(run_dirs, target_parameters = TARGET_PARAMETERS) {
  rows <- list()
  k <- 0L
  for (run_dir in run_dirs) {
    seed_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
    seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
    seed_dirs <- seed_dirs[order(seed_order_key(basename(seed_dirs)), basename(seed_dirs))]
    for (seed_dir in seed_dirs) {
      row <- read_seed_soft_table(run_dir, seed_dir, target_parameters = target_parameters)
      if (is.null(row) || !nrow(row)) next
      k <- k + 1L
      rows[[k]] <- row
    }
  }
  if (!length(rows)) {
    stop("No soft-coupling rows were read from run directories.", call. = FALSE)
  }
  out <- do.call(rbind, rows)
  out$parameter <- factor(as.character(out$parameter), levels = target_parameters)
  out <- out[order(out$parameter_order, out$seed_id, out$seed, out$sigma, out$run_label), , drop = FALSE]
  row.names(out) <- NULL
  out
}

filter_complete_paired_seeds <- function(df, target_parameters = TARGET_PARAMETERS) {
  run_labels <- unique(as.character(df$run_label))
  needed <- length(run_labels) * length(target_parameters)
  keys <- unique(df[, c("seed", "run_label", "parameter"), drop = FALSE])
  counts <- stats::aggregate(
    run_label ~ seed,
    data = keys,
    FUN = length
  )
  names(counts)[names(counts) == "run_label"] <- "complete_cells"
  complete_seeds <- as.character(counts$seed[counts$complete_cells == needed])
  df <- df[as.character(df$seed) %in% complete_seeds, , drop = FALSE]
  df <- df[order(df$parameter_order, df$seed_id, df$seed, df$sigma, df$run_label), , drop = FALSE]
  row.names(df) <- NULL
  attr(df, "complete_seed_count") <- length(unique(as.character(df$seed)))
  attr(df, "all_seed_count") <- length(unique(as.character(keys$seed)))
  df
}

sigma_levels_from_df <- function(df) {
  sigma_map <- unique(df[, c("sigma", "sigma_label"), drop = FALSE])
  sigma_map <- sigma_map[order(sigma_map$sigma, sigma_map$sigma_label, na.last = TRUE), , drop = FALSE]
  as.character(sigma_map$sigma_label)
}

make_value_long <- function(df, sigma_levels) {
  vivo <- data.frame(
    df[, c("run_label", "run_dir", "seed", "seed_id", "parameter", "parameter_order", "sigma", "sigma_label"), drop = FALSE],
    context = "in vivo",
    value = df$vivo_natural,
    clipped = df$vivo_clipped,
    boundary_status = df$boundary_status_vivo,
    stringsAsFactors = FALSE
  )
  vitro <- data.frame(
    df[, c("run_label", "run_dir", "seed", "seed_id", "parameter", "parameter_order", "sigma", "sigma_label"), drop = FALSE],
    context = "in vitro",
    value = df$vitro_natural,
    clipped = df$vitro_clipped,
    boundary_status = df$boundary_status_vitro,
    stringsAsFactors = FALSE
  )
  out <- rbind(vivo, vitro)
  out$context <- factor(out$context, levels = CONTEXT_LEVELS)
  out$sigma_label <- factor(as.character(out$sigma_label), levels = sigma_levels)
  group_levels <- as.vector(t(outer(CONTEXT_LEVELS, sigma_levels, paste, sep = "__")))
  out$value_group <- factor(
    paste(as.character(out$context), as.character(out$sigma_label), sep = "__"),
    levels = group_levels
  )
  out
}

value_group_labels <- function(levels) {
  vapply(strsplit(as.character(levels), "__", fixed = TRUE), function(x) {
    paste0(x[[1]], "\n", x[[2]])
  }, character(1))
}

make_objective_long <- function(df, sigma_levels) {
  seed_df <- unique(df[, c(
    "run_label",
    "run_dir",
    "seed",
    "seed_id",
    "sigma",
    "sigma_label",
    OBJECTIVE_TYPES
  ), drop = FALSE])
  seed_df <- seed_df[order(seed_df$seed_id, seed_df$seed, seed_df$sigma), , drop = FALSE]
  rows <- lapply(OBJECTIVE_TYPES, function(objective_type) {
    data.frame(
      seed_df[, c("run_label", "run_dir", "seed", "seed_id", "sigma", "sigma_label"), drop = FALSE],
      objective_type = objective_type,
      objective_value = suppressWarnings(as.numeric(seed_df[[objective_type]])),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$objective_type <- factor(out$objective_type, levels = OBJECTIVE_TYPES)
  out$sigma_label <- factor(as.character(out$sigma_label), levels = sigma_levels)
  group_levels <- as.vector(t(outer(OBJECTIVE_TYPES, sigma_levels, paste, sep = "__")))
  out$objective_group <- factor(
    paste(as.character(out$objective_type), as.character(out$sigma_label), sep = "__"),
    levels = group_levels
  )
  out <- out[order(out$objective_type, out$seed_id, out$seed, out$sigma), , drop = FALSE]
  row.names(out) <- NULL
  out
}

plot_objective_overview <- function(objective_df) {
  plot_df <- objective_df[is.finite(objective_df$objective_value), , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)
  set.seed(1)
  ggplot(plot_df, aes(x = objective_group, y = objective_value)) +
    geom_boxplot(
      fill = "white",
      color = "grey45",
      linewidth = 0.45,
      outlier.shape = NA,
      width = 0.58
    ) +
    geom_line(
      aes(group = interaction(seed, objective_type, drop = TRUE)),
      color = "grey72",
      linewidth = 0.28,
      alpha = 0.55
    ) +
    geom_point(
      position = position_jitter(width = 0.055, height = 0),
      color = "grey20",
      size = 1.55,
      alpha = 0.82
    ) +
    scale_x_discrete(labels = value_group_labels(levels(plot_df$objective_group))) +
    labs(
      title = "Joint objective components across paired seeds",
      subtitle = "Objective type is the primary group; within each type, sigma values are subgroups. Lines connect the same seed across sigma values.",
      x = NULL,
      y = "Objective value"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8.7),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

plot_parameter_values <- function(value_df, parameter, clipped = FALSE) {
  plot_df <- value_df[as.character(value_df$parameter) == parameter, , drop = FALSE]
  plot_df <- plot_df[is.finite(plot_df$value), , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)
  set.seed(1)
  p <- ggplot(plot_df, aes(x = value_group, y = value)) +
    geom_boxplot(
      fill = "white",
      color = "grey45",
      linewidth = 0.45,
      outlier.shape = NA,
      width = 0.58
    ) +
    geom_line(
      aes(group = interaction(seed, context, drop = TRUE)),
      color = "grey72",
      linewidth = 0.28,
      alpha = 0.55
    )
  if (isTRUE(clipped)) {
    p <- p +
      geom_point(
        aes(color = clipped),
        position = position_jitter(width = 0.055, height = 0),
        size = 1.8,
        alpha = 0.9
      ) +
      scale_color_manual(
        name = "Status",
        values = c("FALSE" = "grey78", "TRUE" = "black"),
        breaks = c(FALSE, TRUE),
        labels = c("not clipped", "clipped"),
        drop = FALSE
      )
  } else {
    p <- p +
      geom_point(
        position = position_jitter(width = 0.055, height = 0),
        color = "grey20",
        size = 1.55,
        alpha = 0.82
      )
  }
  p +
    scale_x_discrete(labels = value_group_labels(levels(plot_df$value_group))) +
    labs(
      title = if (isTRUE(clipped)) {
        paste0(parameter, ": fitted values with clipped status")
      } else {
        paste0(parameter, ": in vivo and in vitro fitted values")
      },
      subtitle = "Within each context, lines connect the same seed across sigma values.",
      x = NULL,
      y = "Natural-scale parameter value"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8.7),
      legend.position = if (isTRUE(clipped)) "bottom" else "none",
      plot.title = element_text(face = "bold")
    )
}

plot_parameter_ratio <- function(df, parameter, clipped = FALSE) {
  plot_df <- df[as.character(df$parameter) == parameter, , drop = FALSE]
  plot_df <- plot_df[is.finite(plot_df$ratio_vivo_to_vitro), , drop = FALSE]
  if (!nrow(plot_df)) return(NULL)
  set.seed(1)
  p <- ggplot(plot_df, aes(x = sigma_label, y = ratio_vivo_to_vitro)) +
    geom_hline(yintercept = 1, color = "grey70", linetype = "dashed", linewidth = 0.35) +
    geom_boxplot(
      fill = "white",
      color = "grey45",
      linewidth = 0.45,
      outlier.shape = NA,
      width = 0.52
    ) +
    geom_line(
      aes(group = seed),
      color = "grey72",
      linewidth = 0.28,
      alpha = 0.55
    )
  if (isTRUE(clipped)) {
    p <- p +
      geom_point(
        aes(color = any_clipped),
        position = position_jitter(width = 0.045, height = 0),
        size = 1.8,
        alpha = 0.9
      ) +
      scale_color_manual(
        name = "Status",
        values = c("FALSE" = "grey78", "TRUE" = "black"),
        breaks = c(FALSE, TRUE),
        labels = c("not clipped", "vivo or vitro clipped"),
        drop = FALSE
      )
  } else {
    p <- p +
      geom_point(
        position = position_jitter(width = 0.045, height = 0),
        color = "grey20",
        size = 1.55,
        alpha = 0.82
      )
  }
  p +
    labs(
      title = if (isTRUE(clipped)) {
        paste0(parameter, ": ratio with clipped status")
      } else {
        paste0(parameter, ": vivo/vitro ratio")
      },
      subtitle = "Lines connect the same seed across sigma values.",
      x = NULL,
      y = "ratio_vivo_to_vitro"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = if (isTRUE(clipped)) "bottom" else "none",
      plot.title = element_text(face = "bold")
    )
}

save_combined_plot <- function(left_plot, right_plot, out_png, out_pdf, width = 14, height = 6.2) {
  if (is.null(left_plot) || is.null(right_plot)) return(invisible(NULL))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(left_plot, right_plot, ncol = 2)
    ggplot2::ggsave(out_png, combined, width = width, height = height, dpi = 180)
    ggplot2::ggsave(out_pdf, combined, width = width, height = height)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined <- gridExtra::arrangeGrob(left_plot, right_plot, ncol = 2)
    ggplot2::ggsave(out_png, combined, width = width, height = height, dpi = 180)
    ggplot2::ggsave(out_pdf, combined, width = width, height = height)
  } else {
    stop("Need either patchwork or gridExtra to save side-by-side combined plots.", call. = FALSE)
  }
  invisible(out_png)
}

save_single_plot <- function(plot, out_png, out_pdf, width = 12, height = 6.2) {
  if (is.null(plot)) return(invisible(NULL))
  ggplot2::ggsave(out_png, plot, width = width, height = height, dpi = 180)
  ggplot2::ggsave(out_pdf, plot, width = width, height = height)
  invisible(out_png)
}

image_data_uri <- function(path, mime = "image/png") {
  if (!file.exists(path)) stop("Missing image to embed in HTML: ", path, call. = FALSE)
  paste0("data:", mime, ";base64,", base64enc::base64encode(path))
}

build_two_sigma_contrast <- function(df) {
  sigma_levels <- sigma_levels_from_df(df)
  if (length(sigma_levels) != 2L) return(data.frame())
  low_label <- sigma_levels[[1]]
  high_label <- sigma_levels[[2]]
  low <- df[as.character(df$sigma_label) == low_label, , drop = FALSE]
  high <- df[as.character(df$sigma_label) == high_label, , drop = FALSE]
  key <- c("seed", "parameter")
  merged <- merge(
    low,
    high,
    by = key,
    suffixes = c("_low_sigma", "_high_sigma"),
    all = FALSE,
    sort = FALSE
  )
  if (!nrow(merged)) return(data.frame())
  ratio_safe <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b != 0, a / b, NA_real_)
  out <- data.frame(
    seed = merged$seed,
    seed_id = merged$seed_id_low_sigma,
    parameter = merged$parameter,
    parameter_order = merged$parameter_order_low_sigma,
    sigma_low = merged$sigma_low_sigma,
    sigma_high = merged$sigma_high_sigma,
    sigma_label_low = merged$sigma_label_low_sigma,
    sigma_label_high = merged$sigma_label_high_sigma,
    run_low = merged$run_label_low_sigma,
    run_high = merged$run_label_high_sigma,
    vivo_natural_low_sigma = merged$vivo_natural_low_sigma,
    vivo_natural_high_sigma = merged$vivo_natural_high_sigma,
    vivo_natural_delta_high_minus_low = merged$vivo_natural_high_sigma - merged$vivo_natural_low_sigma,
    vivo_natural_ratio_high_to_low = ratio_safe(merged$vivo_natural_high_sigma, merged$vivo_natural_low_sigma),
    vitro_natural_low_sigma = merged$vitro_natural_low_sigma,
    vitro_natural_high_sigma = merged$vitro_natural_high_sigma,
    vitro_natural_delta_high_minus_low = merged$vitro_natural_high_sigma - merged$vitro_natural_low_sigma,
    vitro_natural_ratio_high_to_low = ratio_safe(merged$vitro_natural_high_sigma, merged$vitro_natural_low_sigma),
    ratio_vivo_to_vitro_low_sigma = merged$ratio_vivo_to_vitro_low_sigma,
    ratio_vivo_to_vitro_high_sigma = merged$ratio_vivo_to_vitro_high_sigma,
    ratio_vivo_to_vitro_delta_high_minus_low = merged$ratio_vivo_to_vitro_high_sigma - merged$ratio_vivo_to_vitro_low_sigma,
    ratio_vivo_to_vitro_ratio_high_to_low = ratio_safe(merged$ratio_vivo_to_vitro_high_sigma, merged$ratio_vivo_to_vitro_low_sigma),
    vivo_clipped_low_sigma = merged$vivo_clipped_low_sigma,
    vivo_clipped_high_sigma = merged$vivo_clipped_high_sigma,
    vitro_clipped_low_sigma = merged$vitro_clipped_low_sigma,
    vitro_clipped_high_sigma = merged$vitro_clipped_high_sigma,
    any_clipped_low_sigma = merged$any_clipped_low_sigma,
    any_clipped_high_sigma = merged$any_clipped_high_sigma,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[order(out$parameter_order, out$seed_id, out$seed), , drop = FALSE]
  row.names(out) <- NULL
  out
}

write_html_report <- function(out_path, out_dir, summary_path, contrast_path, objective_long_path, objective_spec, plot_specs, df, contrast_df) {
  rel <- function(path) html_escape(basename(path))
  sigma_rows <- unique(df[, c("run_label", "sigma", "sigma_label"), drop = FALSE])
  sigma_rows <- sigma_rows[order(sigma_rows$sigma, sigma_rows$run_label, na.last = TRUE), , drop = FALSE]
  sigma_items <- paste0(
    "<li><code>", html_escape(sigma_rows$run_label), "</code>: ",
    html_escape(sigma_rows$sigma_label), "</li>",
    collapse = "\n"
  )
  section_html <- vapply(TARGET_PARAMETERS, function(param) {
    spec <- plot_specs[[param]]
    if (is.null(spec)) return("")
    main_src <- image_data_uri(spec$main_png)
    clipped_src <- image_data_uri(spec$clipped_png)
    paste0(
      '<section class="parameter-section" id="', html_escape(param), '">\n',
      "<h2>", html_escape(param), "</h2>\n",
      '<figure><img src="', main_src, '" alt="', html_escape(param), ' paired values and ratios">',
      "<figcaption>Natural-scale in vivo/in vitro parameter values and ratio_vivo_to_vitro. Lines connect the same seed across sigma values.</figcaption></figure>\n",
      '<figure><img src="', clipped_src, '" alt="', html_escape(param), ' clipped status">',
      "<figcaption>Black points mark clipped estimates. For ratio panels, black means either vivo_clipped or vitro_clipped is TRUE.</figcaption></figure>\n",
      "</section>\n"
    )
  }, character(1))
  clipped_count <- sum(df$any_clipped, na.rm = TRUE)
  summary_rows <- nrow(df)
  objective_section <- ""
  if (!is.null(objective_spec) && file.exists(objective_spec$png)) {
    objective_src <- image_data_uri(objective_spec$png)
    objective_section <- paste0(
      '<section class="objective-section" id="objectives">\n',
      "<h2>Objectives</h2>\n",
      '<figure><img src="', objective_src, '" alt="Objective components across paired seeds">',
      "<figcaption>Total objective, objective_invivo, and objective_invitro grouped by objective type and split by sigma. Lines connect the same seed across sigma values.</figcaption></figure>\n",
      "</section>\n"
    )
  }
  contrast_link <- if (is.data.frame(contrast_df) && nrow(contrast_df)) {
    paste0('<li>Two-sigma paired contrast table: <a href="', rel(contrast_path), '">', html_escape(basename(contrast_path)), "</a></li>")
  } else {
    ""
  }
  html <- paste0(
    "<!doctype html>\n<html><head><meta charset=\"utf-8\">\n",
    "<title>Joint Sigma Soft-Coupled Paired Seed Comparison</title>\n",
    "<style>\n",
    "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:28px;color:#222;line-height:1.45;background:#fff;}\n",
    "main{max-width:1280px;margin:0 auto;}\n",
    "h1{font-size:28px;margin-bottom:8px;} h2{font-size:22px;margin-top:34px;border-top:1px solid #ddd;padding-top:22px;}\n",
    "p,li{font-size:14px;} code{background:#f4f4f4;padding:2px 4px;border-radius:3px;}\n",
    "figure{margin:18px 0 26px 0;} img{max-width:100%;height:auto;border:1px solid #ddd;}\n",
    "figcaption{font-size:13px;color:#555;margin-top:7px;}\n",
    ".meta{background:#f7f7f7;border:1px solid #ddd;padding:12px 16px;margin:18px 0;}\n",
    "</style></head><body><main>\n",
    "<h1>Joint Sigma Soft-Coupled Paired Seed Comparison</h1>\n",
    "<p>Paired-seed comparison across joint soft-coupling sigma settings. Parameter sections are ordered as p_misseg, O2_crit, alpha_o2, then mu_hp.</p>\n",
    "<div class=\"meta\"><ul>\n",
    "<li>Long summary table: <a href=\"", rel(summary_path), "\">", html_escape(basename(summary_path)), "</a></li>\n",
    "<li>Objective long table: <a href=\"", rel(objective_long_path), "\">", html_escape(basename(objective_long_path)), "</a></li>\n",
    contrast_link, "\n",
    "<li>Rows in long summary: ", summary_rows, "</li>\n",
    "<li>Paired seed count: ", length(unique(as.character(df$seed))), "</li>\n",
    "<li>Rows with vivo_clipped or vitro_clipped TRUE: ", clipped_count, "</li>\n",
    "</ul><p>Runs:</p><ul>", sigma_items, "</ul></div>\n",
    objective_section,
    paste(section_html, collapse = "\n"),
    "</main></body></html>\n"
  )
  writeLines(html, out_path, useBytes = TRUE)
  invisible(out_path)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  results_dir <- normalizePath(args$results_dir %||% resolve_default_results_dir(), mustWork = TRUE)
  out_dir <- normalizePath(args$out_dir %||% results_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_dirs <- resolve_run_dirs(results_dir, args$run_dirs %||% NULL)
  if (length(run_dirs) < 2L) {
    stop("Need at least two run directories for paired sigma comparison.", call. = FALSE)
  }

  message("Reading run directories: ", paste(basename(run_dirs), collapse = ", "))
  raw_df <- read_all_soft_tables(run_dirs, target_parameters = TARGET_PARAMETERS)
  paired_df <- filter_complete_paired_seeds(raw_df, target_parameters = TARGET_PARAMETERS)
  if (!nrow(paired_df)) {
    stop("No complete paired seeds remained after requiring all target parameters in all runs.", call. = FALSE)
  }
  all_seed_count <- attr(paired_df, "all_seed_count")
  complete_seed_count <- attr(paired_df, "complete_seed_count")
  if (!is.null(all_seed_count) && !is.null(complete_seed_count) && all_seed_count != complete_seed_count) {
    message("Keeping ", complete_seed_count, " complete paired seeds out of ", all_seed_count, " seeds observed in any run.")
  } else {
    message("Keeping ", length(unique(as.character(paired_df$seed))), " complete paired seeds.")
  }

  sigma_levels <- sigma_levels_from_df(paired_df)
  paired_df$sigma_label <- factor(as.character(paired_df$sigma_label), levels = sigma_levels)
  paired_df$parameter <- factor(as.character(paired_df$parameter), levels = TARGET_PARAMETERS)
  paired_df <- paired_df[order(paired_df$parameter_order, paired_df$seed_id, paired_df$seed, paired_df$sigma), , drop = FALSE]

  summary_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_summary.tsv"))
  utils::write.table(paired_df, file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

  contrast_df <- build_two_sigma_contrast(paired_df)
  contrast_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_two_sigma_contrast.tsv"))
  if (nrow(contrast_df)) {
    utils::write.table(contrast_df, file = contrast_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  value_long <- make_value_long(paired_df, sigma_levels = sigma_levels)
  value_long_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_value_long.tsv"))
  utils::write.table(value_long, file = value_long_path, sep = "\t", quote = FALSE, row.names = FALSE)

  objective_long <- make_objective_long(paired_df, sigma_levels = sigma_levels)
  objective_long_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_objective_long.tsv"))
  utils::write.table(objective_long, file = objective_long_path, sep = "\t", quote = FALSE, row.names = FALSE)
  objective_plot <- plot_objective_overview(objective_long)
  objective_png <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_objectives.png"))
  objective_pdf <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_objectives.pdf"))
  save_single_plot(objective_plot, objective_png, objective_pdf)
  objective_spec <- list(png = objective_png, pdf = objective_pdf)

  plot_specs <- list()
  for (param in TARGET_PARAMETERS) {
    token <- safe_file_token(param)
    value_plot <- plot_parameter_values(value_long, param, clipped = FALSE)
    ratio_plot <- plot_parameter_ratio(paired_df, param, clipped = FALSE)
    clipped_value_plot <- plot_parameter_values(value_long, param, clipped = TRUE)
    clipped_ratio_plot <- plot_parameter_ratio(paired_df, param, clipped = TRUE)
    main_png <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_", token, "_values_and_ratio.png"))
    main_pdf <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_", token, "_values_and_ratio.pdf"))
    clipped_png <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_", token, "_clipped_values_and_ratio.png"))
    clipped_pdf <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_", token, "_clipped_values_and_ratio.pdf"))
    save_combined_plot(value_plot, ratio_plot, main_png, main_pdf)
    save_combined_plot(clipped_value_plot, clipped_ratio_plot, clipped_png, clipped_pdf)
    plot_specs[[param]] <- list(
      main_png = main_png,
      main_pdf = main_pdf,
      clipped_png = clipped_png,
      clipped_pdf = clipped_pdf
    )
  }

  html_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, ".html"))
  write_html_report(
    out_path = html_path,
    out_dir = out_dir,
    summary_path = summary_path,
    contrast_path = contrast_path,
    objective_long_path = objective_long_path,
    objective_spec = objective_spec,
    plot_specs = plot_specs,
    df = paired_df,
    contrast_df = contrast_df
  )

  message("Wrote summary table: ", summary_path)
  if (nrow(contrast_df)) message("Wrote paired contrast table: ", contrast_path)
  message("Wrote value-long table: ", value_long_path)
  message("Wrote objective-long table: ", objective_long_path)
  if (!is.null(objective_plot) && file.exists(objective_png)) message("Wrote objective plot: ", objective_png)
  message("Wrote HTML report: ", html_path)
}

if (sys.nframe() == 0L) {
  main()
}
