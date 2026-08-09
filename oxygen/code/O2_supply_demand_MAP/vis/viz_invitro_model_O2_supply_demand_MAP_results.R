#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

.o2_bootstrap_script_dir <- local({
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
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
OXYGEN_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
sys.source(
  file.path(WORKFLOW_ROOT, "vis", "invitro", "o2_supply_demand_map_invitro_plot_utils.R"),
  envir = environment(),
  chdir = TRUE
)
sys.source(
  file.path(WORKFLOW_ROOT, "vis", "invitro", "o2_supply_demand_map_invitro_diagnostic_plot_utils.R"),
  envir = environment(),
  chdir = TRUE
)
rm(.o2_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

read_tsv_optional <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
}

num <- o2sd_numeric

finite_rows <- function(df, cols) {
  ok <- rep(TRUE, nrow(df))
  for (col in cols) {
    ok <- ok & is.finite(df[[col]])
  }
  df[ok, , drop = FALSE]
}

summary_value <- function(summary_df, metric) {
  if (is.null(summary_df) || !all(c("metric", "value") %in% names(summary_df))) return(NA_character_)
  idx <- match(metric, summary_df$metric)
  if (is.na(idx)) NA_character_ else summary_df$value[[idx]]
}

headless_png_device <- function(filename, width, height, units, res, bg, ...) {
  args <- list(
    filename = filename,
    width = width,
    height = height,
    units = units,
    res = res,
    bg = bg,
    ...
  )
  if (isTRUE(capabilities("cairo"))) args$type <- "cairo"
  do.call(grDevices::png, args)
}

save_plot_pair <- function(plot, out_dir, basename, width = 9, height = 5.5) {
  pdf_path <- file.path(out_dir, paste0(basename, ".pdf"))
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 180,
    bg = "white",
    device = headless_png_device
  )
  unlink(file.path(out_dir, paste0(basename, ".svg")), force = TRUE)
  invisible(c(pdf = pdf_path, png = png_path))
}

save_existing_plot_png <- function(plot, out_dir, basename, width = 10, height = 7) {
  if (is.null(plot)) return(FALSE)
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  ggplot2::ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 180,
    bg = "white",
    device = headless_png_device
  )
  unlink(file.path(out_dir, paste0(basename, ".svg")), force = TRUE)
  TRUE
}

wrap_ggplot_title <- function(plot, width = 34) {
  title <- plot$labels$title %||% ""
  if (!nzchar(title)) return(plot)
  plot + ggplot2::labs(title = paste(strwrap(title, width = width), collapse = "\n"))
}

derive_invitro_lineage_label <- function(key) {
  key <- as.character(key)
  out <- rep(NA_character_, length(key))
  ok <- !is.na(key) & nzchar(key)
  if (!any(ok)) return(out)
  parts <- strsplit(key[ok], "_", fixed = TRUE)
  is_control <- vapply(
    parts,
    function(x) length(x) > 0L && all(trimws(x) == "20.5"),
    logical(1)
  )
  out[ok] <- ifelse(is_control, "control", "deprived")
  out
}

derive_invitro_lineage_passage_index <- function(key, fallback = NA_real_) {
  key <- as.character(key)
  out <- rep(NA_real_, length(key))
  ok <- !is.na(key) & nzchar(key)
  if (any(ok)) {
    out[ok] <- vapply(strsplit(key[ok], "_", fixed = TRUE), length, integer(1))
  }
  fallback <- suppressWarnings(as.numeric(fallback))
  use_fallback <- !is.finite(out) & is.finite(fallback)
  out[use_fallback] <- fallback[use_fallback]
  out
}

ensure_invitro_plot_columns <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  n <- nrow(df)
  if (!"lineage_terminal_key" %in% names(df)) {
    df$lineage_terminal_key <- if ("segment_id" %in% names(df)) as.character(df$segment_id) else as.character(seq_len(n))
  }
  if (!"lineage_passage_index" %in% names(df)) {
    fallback <- if ("passage_index" %in% names(df)) df$passage_index else seq_len(n)
    df$lineage_passage_index <- derive_invitro_lineage_passage_index(df$lineage_terminal_key, fallback)
  }
  if (!"lineage_label" %in% names(df) || all(is.na(df$lineage_label))) {
    key <- if ("lineage_terminal_key" %in% names(df)) {
      df$lineage_terminal_key
    } else if ("segment_id" %in% names(df)) {
      df$segment_id
    } else {
      rep(NA_character_, n)
    }
    df$lineage_label <- derive_invitro_lineage_label(key)
  }
  missing_label <- is.na(df$lineage_label) | !nzchar(as.character(df$lineage_label))
  if (any(missing_label)) {
    df$lineage_label[missing_label] <- if ("cohort" %in% names(df)) as.character(df$cohort[missing_label]) else "lineage"
  }
  if (!"sample_name" %in% names(df)) {
    df$sample_name <- if ("passage_id" %in% names(df)) as.character(df$passage_id) else df$lineage_terminal_key
  }
  df
}

plot_remote_lineage_counts <- function(lineage_df, out_dir) {
  required <- c("predicted_live_cells", "passage_index", "cohort")
  if (is.null(lineage_df) || !all(required %in% names(lineage_df))) return(invisible(FALSE))
  p <- ivt_plot_lineage_counts(ensure_invitro_plot_columns(lineage_df))
  save_plot_pair(p, out_dir, "invitro_lineage_counts", width = 10, height = 5.6)
  invisible(TRUE)
}

plot_remote_daily_counts <- function(daily_df, out_dir) {
  required <- c("day", "live_cells", "selected_day", "passage_index", "cohort")
  if (is.null(daily_df) || !all(required %in% names(daily_df))) return(invisible(FALSE))
  p <- ivt_plot_daily_counts(ensure_invitro_plot_columns(daily_df))
  save_plot_pair(p, out_dir, "invitro_daily_counts", width = 15, height = 12)
  invisible(TRUE)
}

plot_remote_burden_decomposition <- function(daily_df, out_dir) {
  required <- c(
    "cohort", "segment_id", "passage_index", "day",
    "live_cells", "dead_hypoxia_cells", "dead_buffer_cells"
  )
  if (is.null(daily_df) || !all(required %in% names(daily_df))) return(invisible(FALSE))
  df <- ensure_invitro_plot_columns(daily_df)
  df$passage_index_num <- num(df$lineage_passage_index)
  df$day_num <- num(df$day)
  df$live_cells <- num(df$live_cells)
  df$dead_hypoxia_cells <- num(df$dead_hypoxia_cells)
  df$dead_buffer_cells <- num(df$dead_buffer_cells)
  df <- finite_rows(df, c("passage_index_num", "day_num", "live_cells", "dead_hypoxia_cells", "dead_buffer_cells"))
  df$component_total_cells <- pmax(df$live_cells, 0) + pmax(df$dead_hypoxia_cells, 0) + pmax(df$dead_buffer_cells, 0)
  df <- df[df$component_total_cells > 0, , drop = FALSE]
  if (!nrow(df)) return(invisible(FALSE))
  df$live_fraction <- pmax(df$live_cells, 0) / df$component_total_cells
  df$dead_hypoxia_fraction <- pmax(df$dead_hypoxia_cells, 0) / df$component_total_cells
  df$dead_buffer_fraction <- pmax(df$dead_buffer_cells, 0) / df$component_total_cells
  df$total_fraction <- 1

  seg_duration <- stats::aggregate(
    day_num ~ cohort + lineage_label + segment_id + passage_index_num,
    data = df,
    FUN = function(x) max(x, na.rm = TRUE)
  )
  names(seg_duration)[names(seg_duration) == "day_num"] <- "duration_days"
  seg_duration$duration_days[!is.finite(seg_duration$duration_days) | seg_duration$duration_days <= 0] <- 1
  df <- merge(df, seg_duration, by = c("cohort", "lineage_label", "segment_id", "passage_index_num"), all.x = TRUE, sort = FALSE)
  df$x_passage <- df$passage_index_num - 1 + pmin(pmax(df$day_num / pmax(df$duration_days, 1e-12), 0), 1)
  df$lineage_label <- factor(as.character(df$lineage_label), levels = unique(as.character(df$lineage_label)))
  facet_levels <- unique(paste(as.character(df$cohort), as.character(df$lineage_label), sep = " / "))
  df$facet_label <- factor(paste(as.character(df$cohort), as.character(df$lineage_label), sep = " / "), levels = facet_levels)
  n_facet_cols <- max(1L, length(unique(as.character(df$lineage_label))))

  long <- dplyr::bind_rows(
    data.frame(df[, setdiff(names(df), c("live_fraction", "dead_hypoxia_fraction", "dead_buffer_fraction")), drop = FALSE], component = "Viable", value = df$live_fraction),
    data.frame(df[, setdiff(names(df), c("live_fraction", "dead_hypoxia_fraction", "dead_buffer_fraction")), drop = FALSE], component = "Hypoxia-Origin Dead", value = df$dead_hypoxia_fraction),
    data.frame(df[, setdiff(names(df), c("live_fraction", "dead_hypoxia_fraction", "dead_buffer_fraction")), drop = FALSE], component = "CIN-Associated Dead", value = df$dead_buffer_fraction)
  )
  long$value <- pmin(pmax(num(long$value), 0), 1)
  long$component <- factor(long$component, levels = c("Viable", "Hypoxia-Origin Dead", "CIN-Associated Dead"))
  long <- long |>
    dplyr::filter(is.finite(.data$x_passage), is.finite(.data$value)) |>
    dplyr::group_by(.data$facet_label, .data$segment_id, .data$x_passage, .data$component) |>
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$facet_label, .data$segment_id, .data$x_passage, .data$component) |>
    dplyr::group_by(.data$facet_label, .data$segment_id, .data$x_passage) |>
    dplyr::mutate(
      ymin = pmin(cumsum(dplyr::lag(.data$value, default = 0)), 1),
      ymax = pmin(.data$ymin + .data$value, 1)
    ) |>
    dplyr::ungroup()
  total_line <- df |>
    dplyr::filter(is.finite(.data$x_passage)) |>
    dplyr::group_by(.data$facet_label, .data$segment_id, .data$x_passage) |>
    dplyr::summarise(total_fraction = 1, .groups = "drop")

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$x_passage)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$component,
        group = interaction(.data$segment_id, .data$component)
      ),
      alpha = 0.55,
      color = NA
    ) +
    ggplot2::geom_line(
      data = total_line,
      ggplot2::aes(x = .data$x_passage, y = .data$total_fraction, group = .data$segment_id),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.55
    ) +
    ggplot2::scale_fill_manual(
      values = c("Viable" = "#1f77b4", "Hypoxia-Origin Dead" = "#d62728", "CIN-Associated Dead" = "#2ca02c")
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::facet_wrap(~facet_label, scales = "free_x", ncol = n_facet_cols) +
    ggplot2::labs(
      title = "In Vitro Predicted Viable/Dead Burden Fractions",
      subtitle = "Component fractions are normalized by the displayed viable/dead predicted cells; total fraction (black) is 1",
      x = "Lineage passage",
      y = "Fraction of total cells",
      fill = "Cell component"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_burden_live_dead_decomposition", width = 12, height = 7)
  invisible(TRUE)
}

order_invitro_cohort <- function(x) {
  x_chr <- as.character(x)
  preferred <- c("2N", "4N")
  levels <- c(preferred[preferred %in% x_chr], sort(setdiff(unique(x_chr), preferred)))
  factor(x_chr, levels = unique(levels))
}

order_invitro_lineage <- function(x) {
  x_chr <- as.character(x)
  preferred <- c("C", "O1", "O2", "control", "deprived")
  levels <- c(preferred[preferred %in% x_chr], sort(setdiff(unique(x_chr), preferred)))
  factor(x_chr, levels = unique(levels))
}

invitro_lineage_palette <- function(levels) {
  base <- c(
    C = "#4e79a7",
    O1 = "#d95f02",
    O2 = "#59a14f",
    control = "#4e79a7",
    deprived = "#d95f02",
    `2N` = "#4e79a7",
    `4N` = "#d95f02"
  )
  missing <- setdiff(levels, names(base))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[levels]
}

make_invitro_facet_label <- function(cohort, panel) {
  paste(as.character(cohort), as.character(panel), sep = "\n")
}

format_invitro_axis_oxygen <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  vapply(x, function(v) {
    if (!is.finite(v)) return("fixed O2=NA")
    paste0("fixed O2=", trimws(formatC(v, format = "fg", digits = 4)), "%")
  }, character(1))
}

first_finite_num <- o2sd_first_finite_numeric

min_finite_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

first_nonempty_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else NA_character_
}

invitro_branch_axis_source <- function(df) {
  empty <- data.frame(
    cohort = character(),
    lineage_label = character(),
    segment_id = character(),
    parent_segment_id = character(),
    lineage_passage_index = numeric(),
    passage_index = numeric(),
    oxygen_pct = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(empty)
  df <- ensure_invitro_plot_columns(df)
  if (!"segment_id" %in% names(df)) {
    df$segment_id <- df$lineage_terminal_key %||% as.character(seq_len(nrow(df)))
  }
  if (!"oxygen_pct" %in% names(df)) df$oxygen_pct <- NA_real_
  if (!"passage_index" %in% names(df)) df$passage_index <- df$lineage_passage_index
  out <- data.frame(
    cohort = as.character(df$cohort),
    lineage_label = as.character(df$lineage_label),
    segment_id = as.character(df$segment_id),
    parent_segment_id = if ("parent_segment_id" %in% names(df)) as.character(df$parent_segment_id) else NA_character_,
    lineage_passage_index = num(df$lineage_passage_index),
    passage_index = num(df$passage_index),
    oxygen_pct = num(df$oxygen_pct),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$cohort) & nzchar(out$lineage_label) & nzchar(out$segment_id), , drop = FALSE]
  out
}

build_invitro_branch_axis_map <- function(...) {
  src <- dplyr::bind_rows(lapply(list(...), invitro_branch_axis_source))
  if (is.null(src) || !nrow(src)) {
    return(data.frame())
  }
  src <- src |>
    dplyr::group_by(.data$cohort, .data$lineage_label, .data$segment_id) |>
    dplyr::summarise(
      parent_segment_id = first_nonempty_chr(.data$parent_segment_id),
      lineage_passage_index = min_finite_num(.data$lineage_passage_index),
      passage_index = min_finite_num(.data$passage_index),
      oxygen_pct = first_finite_num(.data$oxygen_pct),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      lineage_label_order = as.integer(order_invitro_lineage(.data$lineage_label)),
      cohort_order = as.integer(order_invitro_cohort(.data$cohort))
    ) |>
    dplyr::arrange(
      .data$cohort_order,
      .data$lineage_label_order,
      .data$lineage_passage_index,
      .data$passage_index,
      .data$segment_id
    ) |>
    dplyr::group_by(.data$cohort, .data$lineage_label) |>
    dplyr::mutate(x_passage_axis = dplyr::row_number()) |>
    dplyr::ungroup()

  src <- src |>
    dplyr::group_by(.data$cohort, .data$lineage_label, .data$lineage_passage_index) |>
    dplyr::mutate(branch_duplicate = dplyr::n() > 1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      passage_label = paste0("p", format(.data$lineage_passage_index, trim = TRUE, scientific = FALSE)),
      x_label_axis = ifelse(
        .data$branch_duplicate,
        paste0(
          .data$passage_label,
          "\nP",
          format(.data$passage_index, trim = TRUE, scientific = FALSE),
          " ",
          format_invitro_axis_oxygen(.data$oxygen_pct)
        ),
        .data$passage_label
      )
    )

  src[, c(
    "cohort", "lineage_label", "segment_id", "parent_segment_id", "lineage_passage_index",
    "passage_index", "oxygen_pct", "x_passage_axis", "x_label_axis", "branch_duplicate"
  ), drop = FALSE]
}

attach_invitro_branch_axis <- function(df, axis_map) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df) || is.null(axis_map) || !nrow(axis_map)) {
    if (!is.null(df) && is.data.frame(df) && !"x_passage" %in% names(df)) {
      df$x_passage <- if ("lineage_passage_index" %in% names(df)) num(df$lineage_passage_index) else num(df$passage_index)
      df$x_label_axis <- paste0("p", format(df$x_passage, trim = TRUE, scientific = FALSE))
    }
    return(df)
  }
  if (!"segment_id" %in% names(df)) {
    df$segment_id <- if ("lineage_terminal_key" %in% names(df)) as.character(df$lineage_terminal_key) else as.character(seq_len(nrow(df)))
  }
  df$.axis_row <- seq_len(nrow(df))
  df$.axis_cohort <- as.character(df$cohort)
  df$.axis_lineage_label <- as.character(df$lineage_label)
  df$.axis_segment_id <- as.character(df$segment_id)

  map <- axis_map[, c("cohort", "lineage_label", "segment_id", "parent_segment_id", "x_passage_axis", "x_label_axis", "branch_duplicate"), drop = FALSE]
  map$.axis_cohort <- as.character(map$cohort)
  map$.axis_lineage_label <- as.character(map$lineage_label)
  map$.axis_segment_id <- as.character(map$segment_id)
  map$axis_parent_segment_id <- as.character(map$parent_segment_id)
  map <- map[, c(".axis_cohort", ".axis_lineage_label", ".axis_segment_id", "axis_parent_segment_id", "x_passage_axis", "x_label_axis", "branch_duplicate"), drop = FALSE]

  out <- merge(df, map, by = c(".axis_cohort", ".axis_lineage_label", ".axis_segment_id"), all.x = TRUE, sort = FALSE)
  out <- out[order(out$.axis_row), , drop = FALSE]
  if (!"parent_segment_id" %in% names(out)) out$parent_segment_id <- NA_character_
  missing_parent <- is.na(out$parent_segment_id) | !nzchar(as.character(out$parent_segment_id))
  out$parent_segment_id[missing_parent] <- out$axis_parent_segment_id[missing_parent]
  out$x_passage <- num(out$x_passage_axis)
  fallback <- if ("lineage_passage_index" %in% names(out)) num(out$lineage_passage_index) else num(out$passage_index)
  missing_x <- !is.finite(out$x_passage)
  out$x_passage[missing_x] <- fallback[missing_x]
  missing_label <- is.na(out$x_label_axis) | !nzchar(as.character(out$x_label_axis))
  out$x_label_axis[missing_label] <- paste0("p", format(out$x_passage[missing_label], trim = TRUE, scientific = FALSE))
  out[, setdiff(names(out), c(".axis_row", ".axis_cohort", ".axis_lineage_label", ".axis_segment_id", "axis_parent_segment_id")), drop = FALSE]
}

nearest_quantile_summary <- function(quantile_df, target, value_name) {
  if (is.null(quantile_df) || !nrow(quantile_df)) {
    return(data.frame())
  }
  quantile_df$distance_to_target <- abs(num(quantile_df$quantile_prob) - as.numeric(target))
  out <- quantile_df |>
    dplyr::group_by(.data$cohort, .data$lineage_label, .data$x_passage) |>
    dplyr::filter(.data$distance_to_target == min(.data$distance_to_target, na.rm = TRUE)) |>
    dplyr::summarise(value = mean(.data$predicted_quantile_kary_N, na.rm = TRUE), .groups = "drop")
  names(out)[names(out) == "value"] <- value_name
  out
}

.invitro_scenario_kary_rows <- function(observed_kary_df) {
  if (is.null(observed_kary_df) || !is.data.frame(observed_kary_df) || !nrow(observed_kary_df)) {
    return(observed_kary_df)
  }
  is_initial <- rep(FALSE, nrow(observed_kary_df))
  for (column in intersect(c("lineage_id", "lineage_label"), names(observed_kary_df))) {
    is_initial <- is_initial |
      toupper(trimws(as.character(observed_kary_df[[column]]))) == "INITIAL"
  }
  if ("lineage_group" %in% names(observed_kary_df)) {
    is_initial <- is_initial |
      tolower(trimws(as.character(observed_kary_df$lineage_group))) == "initial"
  }
  if ("scenario_id" %in% names(observed_kary_df)) {
    is_initial <- is_initial |
      grepl("(^|-)INITIAL$", as.character(observed_kary_df$scenario_id), ignore.case = TRUE)
  }
  is_initial[is.na(is_initial)] <- FALSE
  observed_kary_df[!is_initial, , drop = FALSE]
}

plot_remote_growth_ploidy_burden_composite <- function(lineage_df,
                                                       quantile_df,
                                                       observed_kary_df,
                                                       daily_df,
                                                       out_dir) {
  if (is.null(lineage_df) || is.null(quantile_df) || is.null(daily_df)) {
    return(invisible(FALSE))
  }
  if (!all(c("cohort", "predicted_growth_rate", "passage_index") %in% names(lineage_df)) ||
      !all(c("cohort", "quantile_prob", "predicted_quantile_kary_N", "passage_index") %in% names(quantile_df)) ||
      !all(c("cohort", "segment_id", "passage_index", "day", "live_cells", "dead_hypoxia_cells", "dead_buffer_cells") %in% names(daily_df))) {
    return(invisible(FALSE))
  }

  scenario_observed_kary_df <- .invitro_scenario_kary_rows(observed_kary_df)
  axis_map <- build_invitro_branch_axis_map(lineage_df, quantile_df, scenario_observed_kary_df, daily_df)

  lin <- ensure_invitro_plot_columns(lineage_df)
  lin$cohort <- order_invitro_cohort(lin$cohort)
  lin$lineage_label <- order_invitro_lineage(lin$lineage_label)
  lin <- attach_invitro_branch_axis(lin, axis_map)
  lin$x_passage <- num(lin$x_passage)
  lin$predicted_growth_rate <- num(lin$predicted_growth_rate)
  lin$observed_growth <- if ("observed_growth" %in% names(lin)) num(lin$observed_growth) else NA_real_

  summarise_segment_nodes <- function(df, value_col, panel_name, extra_group_cols = character()) {
    if (is.null(df) || !nrow(df) || !value_col %in% names(df)) return(data.frame())
    if (!"segment_id" %in% names(df)) df$segment_id <- as.character(seq_len(nrow(df)))
    if (!"parent_segment_id" %in% names(df)) df$parent_segment_id <- NA_character_
    df$segment_id <- as.character(df$segment_id)
    df$parent_segment_id <- as.character(df$parent_segment_id)
    df[[value_col]] <- num(df[[value_col]])
    group_cols <- unique(c("cohort", "lineage_label", "segment_id", "parent_segment_id", "x_passage", extra_group_cols))
    if (length(setdiff(group_cols, names(df)))) return(data.frame())
    df |>
      dplyr::filter(is.finite(.data$x_passage), is.finite(.data[[value_col]])) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(is.finite(.data$value)) |>
      dplyr::mutate(panel = panel_name)
  }

  make_parent_child_edges <- function(nodes, extra_group_cols = character()) {
    required <- unique(c("cohort", "lineage_label", "segment_id", "parent_segment_id", "x_passage", "value", extra_group_cols))
    if (is.null(nodes) || !nrow(nodes) || length(setdiff(required, names(nodes)))) return(data.frame())
    nodes <- nodes[is.finite(num(nodes$x_passage)) & is.finite(num(nodes$value)), , drop = FALSE]
    nodes$segment_id <- as.character(nodes$segment_id)
    nodes$parent_segment_id <- as.character(nodes$parent_segment_id)
    child <- nodes[!is.na(nodes$parent_segment_id) & nzchar(nodes$parent_segment_id), , drop = FALSE]
    if (!nrow(child)) return(data.frame())
    parent <- nodes[, unique(c("cohort", "lineage_label", "segment_id", extra_group_cols, "x_passage", "value")), drop = FALSE]
    names(parent)[names(parent) == "segment_id"] <- "parent_segment_id"
    names(parent)[names(parent) == "x_passage"] <- "x_parent"
    names(parent)[names(parent) == "value"] <- "y_parent"
    by_cols <- unique(c("cohort", "lineage_label", "parent_segment_id", extra_group_cols))
    edges <- merge(child, parent, by = by_cols, all = FALSE, sort = FALSE)
    if (!nrow(edges)) return(data.frame())
    edges$x_parent <- num(edges$x_parent)
    edges$y_parent <- num(edges$y_parent)
    edges$x_passage <- num(edges$x_passage)
    edges$value <- num(edges$value)
    edges <- edges[
      is.finite(edges$x_parent) & is.finite(edges$y_parent) &
        is.finite(edges$x_passage) & is.finite(edges$value),
      ,
      drop = FALSE
    ]
    edges$edge_id <- seq_len(nrow(edges))
    edges
  }

  growth_pred <- summarise_segment_nodes(lin, "predicted_growth_rate", "Growth Rate Fit")
  growth_pred_edges <- make_parent_child_edges(growth_pred)
  growth_obs <- lin |>
    dplyr::filter(is.finite(.data$x_passage), is.finite(.data$observed_growth)) |>
    dplyr::transmute(
      cohort = .data$cohort,
      lineage_label = .data$lineage_label,
      x_passage = .data$x_passage,
      value = .data$observed_growth
    ) |>
    dplyr::mutate(panel = "Growth Rate Fit")

  qdf <- ensure_invitro_plot_columns(quantile_df)
  qdf$cohort <- order_invitro_cohort(qdf$cohort)
  qdf$lineage_label <- order_invitro_lineage(qdf$lineage_label)
  qdf <- attach_invitro_branch_axis(qdf, axis_map)
  qdf$x_passage <- num(qdf$x_passage)
  qdf$quantile_prob <- num(qdf$quantile_prob)
  qdf$predicted_quantile_kary_N <- num(qdf$predicted_quantile_kary_N)
  qdf <- qdf |>
    dplyr::filter(is.finite(.data$x_passage), is.finite(.data$quantile_prob), is.finite(.data$predicted_quantile_kary_N))
  if (!nrow(qdf)) return(invisible(FALSE))
  ploidy_lines <- summarise_segment_nodes(
    qdf,
    "predicted_quantile_kary_N",
    "Chromosome-Number Quantile Fit",
    extra_group_cols = "quantile_prob"
  )
  if (!nrow(ploidy_lines)) return(invisible(FALSE))
  ploidy_line_edges <- make_parent_child_edges(ploidy_lines, extra_group_cols = "quantile_prob")

  ploidy_obs <- data.frame()
  if (!is.null(scenario_observed_kary_df) && "observed_kary_N" %in% names(scenario_observed_kary_df)) {
    obs <- ensure_invitro_plot_columns(scenario_observed_kary_df)
    obs$cohort <- order_invitro_cohort(obs$cohort)
    obs$lineage_label <- order_invitro_lineage(obs$lineage_label)
    obs <- attach_invitro_branch_axis(obs, axis_map)
    obs$x_passage <- num(obs$x_passage)
    obs$observed_kary_N <- num(obs$observed_kary_N)
    ploidy_obs <- obs |>
      dplyr::filter(is.finite(.data$x_passage), is.finite(.data$observed_kary_N)) |>
      dplyr::distinct(cohort, lineage_label, x_passage, passage_id, cell_index, observed_kary_N) |>
      dplyr::transmute(
        cohort = .data$cohort,
        lineage_label = .data$lineage_label,
        x_passage = .data$x_passage,
        value = .data$observed_kary_N
      ) |>
      dplyr::mutate(panel = "Chromosome-Number Quantile Fit")
  } else if ("observed_mean_kary_N" %in% names(lin)) {
    ploidy_obs <- lin |>
      dplyr::filter(is.finite(.data$x_passage), is.finite(.data$observed_mean_kary_N)) |>
      dplyr::transmute(
        cohort = .data$cohort,
        lineage_label = .data$lineage_label,
        x_passage = .data$x_passage,
        value = .data$observed_mean_kary_N
      ) |>
      dplyr::mutate(panel = "Chromosome-Number Quantile Fit")
  }

  burden_df <- ensure_invitro_plot_columns(daily_df)
  burden_df$cohort <- order_invitro_cohort(burden_df$cohort)
  burden_df$lineage_label <- order_invitro_lineage(burden_df$lineage_label)
  burden_df <- attach_invitro_branch_axis(burden_df, axis_map)
  burden_df$branch_x_passage <- num(burden_df$x_passage)
  burden_df$passage_index_num <- num(burden_df$lineage_passage_index)
  burden_df$day_num <- num(burden_df$day)
  burden_df$live_cells <- num(burden_df$live_cells)
  burden_df$dead_hypoxia_cells <- num(burden_df$dead_hypoxia_cells)
  burden_df$dead_buffer_cells <- num(burden_df$dead_buffer_cells)
  burden_df <- finite_rows(burden_df, c("branch_x_passage", "passage_index_num", "day_num", "live_cells", "dead_hypoxia_cells", "dead_buffer_cells"))
  burden_df$component_total_cells <- pmax(burden_df$live_cells, 0) +
    pmax(burden_df$dead_hypoxia_cells, 0) +
    pmax(burden_df$dead_buffer_cells, 0)
  burden_df <- burden_df[burden_df$component_total_cells > 0, , drop = FALSE]
  if (!nrow(burden_df)) return(invisible(FALSE))
  burden_df$live_fraction <- pmax(burden_df$live_cells, 0) / burden_df$component_total_cells
  burden_df$dead_hypoxia_fraction <- pmax(burden_df$dead_hypoxia_cells, 0) / burden_df$component_total_cells
  burden_df$dead_buffer_fraction <- pmax(burden_df$dead_buffer_cells, 0) / burden_df$component_total_cells
  seg_duration <- stats::aggregate(
    day_num ~ cohort + lineage_label + segment_id + passage_index_num,
    data = burden_df,
    FUN = function(x) max(x, na.rm = TRUE)
  )
  names(seg_duration)[names(seg_duration) == "day_num"] <- "duration_days"
  seg_duration$duration_days[!is.finite(seg_duration$duration_days) | seg_duration$duration_days <= 0] <- 1
  burden_df <- merge(burden_df, seg_duration, by = c("cohort", "lineage_label", "segment_id", "passage_index_num"), all.x = TRUE, sort = FALSE)
  burden_df$x_passage <- burden_df$branch_x_passage - 1 + pmin(pmax(burden_df$day_num / pmax(burden_df$duration_days, 1e-12), 0), 1)
  burden_long <- dplyr::bind_rows(
    data.frame(burden_df[, c("cohort", "lineage_label", "segment_id", "x_passage"), drop = FALSE], component = "Viable", value = burden_df$live_fraction),
    data.frame(burden_df[, c("cohort", "lineage_label", "segment_id", "x_passage"), drop = FALSE], component = "Hypoxia-Origin Dead", value = burden_df$dead_hypoxia_fraction),
    data.frame(burden_df[, c("cohort", "lineage_label", "segment_id", "x_passage"), drop = FALSE], component = "CIN-Associated Dead", value = burden_df$dead_buffer_fraction)
  )
  burden_long$value <- pmin(pmax(num(burden_long$value), 0), 1)
  burden_long$component <- factor(burden_long$component, levels = c("Viable", "Hypoxia-Origin Dead", "CIN-Associated Dead"))
  burden_agg <- burden_long |>
    dplyr::filter(is.finite(.data$x_passage), is.finite(.data$value)) |>
    dplyr::group_by(.data$cohort, .data$lineage_label, .data$segment_id, .data$x_passage, .data$component) |>
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$cohort, .data$lineage_label, .data$segment_id, .data$x_passage, .data$component) |>
    dplyr::group_by(.data$cohort, .data$lineage_label, .data$segment_id, .data$x_passage) |>
    dplyr::mutate(ymin = cumsum(dplyr::lag(.data$value, default = 0)), ymax = .data$ymin + .data$value) |>
    dplyr::ungroup() |>
    dplyr::mutate(panel = "Predicted Viable/Dead Burden Fractions")
  burden_total <- burden_agg |>
    dplyr::group_by(.data$cohort, .data$lineage_label, .data$segment_id, .data$x_passage, .data$panel) |>
    dplyr::summarise(value = pmin(max(.data$ymax, na.rm = TRUE), 1), .groups = "drop")

  panel_labels <- c(
    "Growth Rate Fit" = "Growth rate",
    "Chromosome-Number Quantile Fit" = "Chr number",
    "Predicted Viable/Dead Burden Fractions" = "Burden fraction"
  )
  branch_markers <- data.frame(
    cohort = character(),
    lineage_label = character(),
    xintercept = numeric(),
    panel_label = factor(character(), levels = unname(panel_labels)),
    stringsAsFactors = FALSE
  )
  if (!is.null(axis_map) && nrow(axis_map)) {
    duplicate_axis <- axis_map[isTRUE(axis_map$branch_duplicate) | as.logical(axis_map$branch_duplicate), , drop = FALSE]
    duplicate_axis <- duplicate_axis[!is.na(duplicate_axis$branch_duplicate) & duplicate_axis$branch_duplicate, , drop = FALSE]
    if (nrow(duplicate_axis)) {
      branch_markers <- do.call(rbind, lapply(unname(panel_labels), function(panel_value) {
        data.frame(
          cohort = duplicate_axis$cohort,
          lineage_label = duplicate_axis$lineage_label,
          xintercept = duplicate_axis$x_passage_axis,
          panel_label = factor(panel_value, levels = unname(panel_labels)),
          stringsAsFactors = FALSE
        )
      }))
    }
  }
  cohort_levels <- levels(order_invitro_cohort(c(as.character(lin$cohort), as.character(qdf$cohort), as.character(burden_df$cohort))))
  normalise_plot_df <- function(df) {
    if (is.null(df) || !nrow(df)) return(df)
    df$cohort <- order_invitro_cohort(df$cohort)
    df$lineage_label <- order_invitro_lineage(df$lineage_label)
    if ("panel" %in% names(df)) {
      labels <- unname(panel_labels[as.character(df$panel)])
      labels[is.na(labels)] <- as.character(df$panel)[is.na(labels)]
      df$panel_label <- factor(labels, levels = unname(panel_labels))
    }
    df
  }
  growth_pred <- normalise_plot_df(growth_pred)
  growth_pred_edges <- normalise_plot_df(growth_pred_edges)
  growth_obs <- normalise_plot_df(growth_obs)
  ploidy_lines <- normalise_plot_df(ploidy_lines)
  ploidy_line_edges <- normalise_plot_df(ploidy_line_edges)
  ploidy_obs <- normalise_plot_df(ploidy_obs)
  burden_agg <- normalise_plot_df(burden_agg)
  burden_total <- normalise_plot_df(burden_total)

  padded_range <- function(x, include_zero = FALSE, pad_fraction = 0.05) {
    x <- num(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(c(0, 1))
    rng <- range(x, na.rm = TRUE)
    if (isTRUE(include_zero)) rng[1] <- min(rng[1], 0)
    span <- diff(rng)
    if (!is.finite(span) || span <= 0) {
      delta <- max(abs(rng[1]) * 0.10, 0.5)
      rng <- rng + c(-delta, delta)
    } else {
      rng <- rng + c(-span, span) * pad_fraction
    }
    rng
  }
  growth_y_limits <- padded_range(c(growth_pred$value, growth_obs$value), include_zero = TRUE)
  ploidy_y_limits <- padded_range(c(ploidy_lines$value, ploidy_obs$value), include_zero = FALSE)
  growth_y_breaks <- pretty(growth_y_limits, n = 5)
  ploidy_y_breaks <- pretty(ploidy_y_limits, n = 5)
  subset_cohort <- function(df, cohort_value) {
    if (is.null(df) || !nrow(df) || !"cohort" %in% names(df)) return(df)
    df[as.character(df$cohort) == as.character(cohort_value), , drop = FALSE]
  }
  subset_cohort_lineage <- function(df, cohort_value, lineage_value) {
    if (is.null(df) || !nrow(df) || !all(c("cohort", "lineage_label") %in% names(df))) return(df)
    df[
      as.character(df$cohort) == as.character(cohort_value) &
        as.character(df$lineage_label) == as.character(lineage_value),
      ,
      drop = FALSE
    ]
  }
  make_facet_scale_formula <- function(condition, scale, parent = parent.frame()) {
    env <- new.env(parent = parent)
    env$scale <- scale
    f <- stats::as.formula(paste(condition, "~ scale"))
    environment(f) <- env
    f
  }
  x_rows_for <- function(cohort_value = NULL) {
    x_col_subset <- function(df) {
      cols <- c("cohort", "lineage_label", "x_passage")
      if (is.null(df) || !nrow(df) || !all(cols %in% names(df))) {
        return(data.frame(cohort = character(), lineage_label = character(), x_passage = numeric()))
      }
      df[, cols, drop = FALSE]
    }
    x_rows <- dplyr::bind_rows(
      x_col_subset(growth_pred),
      x_col_subset(growth_obs),
      x_col_subset(ploidy_lines),
      x_col_subset(ploidy_obs),
      x_col_subset(burden_agg)
    )
    if (!is.null(cohort_value) && nrow(x_rows)) {
      x_rows <- x_rows[as.character(x_rows$cohort) == as.character(cohort_value), , drop = FALSE]
    }
    x_rows
  }
  lineage_x_meta <- function(cohort_value = NULL) {
    x_rows <- x_rows_for(cohort_value)
    if (!nrow(x_rows)) {
      return(data.frame())
    }
    lineage_levels <- levels(order_invitro_lineage(x_rows$lineage_label))
    lineage_levels <- lineage_levels[lineage_levels %in% as.character(x_rows$lineage_label)]
    rows <- lapply(lineage_levels, function(lineage_value) {
      lineage_x <- num(x_rows$x_passage[as.character(x_rows$lineage_label) == lineage_value])
      lineage_x <- lineage_x[is.finite(lineage_x)]
      if (!length(lineage_x)) return(NULL)
      x_lower <- min(0, floor(min(lineage_x, na.rm = TRUE)))
      x_upper <- max(1, ceiling(max(lineage_x, na.rm = TRUE)))
      x_break_by <- if (x_upper > 25) 5 else if (x_upper > 12) 2 else 1
      data.frame(
        lineage_label = lineage_value,
        x_lower = x_lower,
        x_upper = x_upper,
        x_break_by = x_break_by,
        span = max(1, x_upper - x_lower),
        stringsAsFactors = FALSE
      )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(data.frame())
    do.call(rbind, rows)
  }
  global_lineage_x_meta <- lineage_x_meta()
  global_lineage_width <- function(lineage_value, fallback) {
    idx <- which(as.character(global_lineage_x_meta$lineage_label) == as.character(lineage_value))
    if (length(idx) && is.finite(global_lineage_x_meta$span[idx[[1]]])) {
      return(global_lineage_x_meta$span[idx[[1]]])
    }
    fallback
  }
  axis_ticks_for <- function(cohort_value, lineage_value, x_break_by = 1) {
    if (is.null(axis_map) || !nrow(axis_map)) {
      return(data.frame(x_passage = numeric(), x_label = character(), stringsAsFactors = FALSE))
    }
    ticks <- axis_map[
      as.character(axis_map$cohort) == as.character(cohort_value) &
        as.character(axis_map$lineage_label) == as.character(lineage_value),
      ,
      drop = FALSE
    ]
    if (!nrow(ticks)) {
      return(data.frame(x_passage = numeric(), x_label = character(), stringsAsFactors = FALSE))
    }
    ticks <- ticks[order(ticks$x_passage_axis), , drop = FALSE]
    x_break_by <- suppressWarnings(as.numeric(x_break_by))
    if (!is.finite(x_break_by) || x_break_by < 1) x_break_by <- 1
    show_tick <- rep(TRUE, nrow(ticks))
    if (x_break_by > 1) {
      branch_tick <- as.logical(ticks$branch_duplicate)
      branch_tick[is.na(branch_tick)] <- FALSE
      show_tick <- abs(ticks$x_passage_axis %% x_break_by) < 1e-8 | branch_tick
      show_tick[is.na(show_tick)] <- FALSE
      if (!any(show_tick)) show_tick <- rep(TRUE, nrow(ticks))
    }
    data.frame(
      x_passage = ticks$x_passage_axis[show_tick],
      x_label = ticks$x_label_axis[show_tick],
      stringsAsFactors = FALSE
    )
  }
  y_scale_formulas <- function() {
    list(
      make_facet_scale_formula(
        "panel_label == 'Growth rate'",
        ggplot2::scale_y_continuous(limits = growth_y_limits, breaks = growth_y_breaks)
      ),
      make_facet_scale_formula(
        "panel_label == 'Chr number'",
        ggplot2::scale_y_continuous(limits = ploidy_y_limits, breaks = ploidy_y_breaks)
      ),
      make_facet_scale_formula(
        "panel_label == 'Burden fraction'",
        ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25))
      )
    )
  }
  make_lineage_plot <- function(cohort_value, lineage_value, x_lower, x_upper, x_break_by, show_legend = FALSE) {
    lineage_panel_label <- as.character(lineage_value)
    with_lineage_panel <- function(df) {
      out <- subset_cohort_lineage(df, cohort_value, lineage_value)
      if (!is.null(out)) {
        out$lineage_panel_label <- factor(rep(lineage_panel_label, nrow(out)), levels = lineage_panel_label)
      }
      out
    }
    burden_agg_lineage <- with_lineage_panel(burden_agg)
    burden_total_lineage <- with_lineage_panel(burden_total)
    ploidy_lines_lineage <- with_lineage_panel(ploidy_lines)
    ploidy_line_edges_lineage <- with_lineage_panel(ploidy_line_edges)
    ploidy_obs_lineage <- with_lineage_panel(ploidy_obs)
    growth_pred_lineage <- with_lineage_panel(growth_pred)
    growth_pred_edges_lineage <- with_lineage_panel(growth_pred_edges)
    growth_obs_lineage <- with_lineage_panel(growth_obs)
    branch_markers_lineage <- with_lineage_panel(branch_markers)
    axis_ticks <- axis_ticks_for(cohort_value, lineage_value, x_break_by)
    x_breaks <- if (nrow(axis_ticks)) axis_ticks$x_passage else seq(x_lower, x_upper, by = x_break_by)
    x_labels <- if (nrow(axis_ticks)) axis_ticks$x_label else x_breaks
    zero_line <- data.frame(
      panel_label = factor("Growth rate", levels = unname(panel_labels)),
      lineage_panel_label = factor(lineage_panel_label, levels = lineage_panel_label),
      yintercept = 0,
      stringsAsFactors = FALSE
    )
    chr_ref_line <- data.frame(
      panel_label = factor("Chr count", levels = unname(panel_labels)),
      lineage_panel_label = factor(lineage_panel_label, levels = lineage_panel_label),
      yintercept = c(44, 88),
      stringsAsFactors = FALSE
    )
    ggplot2::ggplot() +
      {
        if (!is.null(branch_markers_lineage) && nrow(branch_markers_lineage)) {
          ggplot2::geom_vline(
            data = branch_markers_lineage,
            ggplot2::aes(xintercept = .data$xintercept),
            color = "black",
            linetype = "22",
            linewidth = 0.28,
            alpha = 0.45
          )
        } else {
          NULL
        }
      } +
      ggplot2::geom_ribbon(
        data = burden_agg_lineage,
        ggplot2::aes(
          x = .data$x_passage,
          ymin = .data$ymin,
          ymax = .data$ymax,
          fill = .data$component,
          group = interaction(.data$segment_id, .data$component)
        ),
        alpha = 0.30,
        color = NA
      ) +
      ggplot2::geom_line(
        data = burden_total_lineage,
        ggplot2::aes(x = .data$x_passage, y = .data$value, group = .data$segment_id),
        color = "black",
        linewidth = 0.55,
        alpha = 0.9
      ) +
      ggplot2::geom_hline(
        data = zero_line,
        ggplot2::aes(yintercept = .data$yintercept),
        color = "black",
        linewidth = 0.25
      ) +
      ggplot2::geom_segment(
        data = ploidy_line_edges_lineage,
        ggplot2::aes(
          x = .data$x_parent,
          y = .data$y_parent,
          xend = .data$x_passage,
          yend = .data$value
        ),
        color = "#1b9e77",
        linewidth = 0.38,
        alpha = 0.42
      ) +
      ggplot2::geom_point(
        data = ploidy_obs_lineage,
        ggplot2::aes(x = .data$x_passage, y = .data$value),
        size = 0.85,
        alpha = 0.62,
        color = "#d95f02",
        position = ggplot2::position_jitter(width = 0.10, height = 0)
      ) +
      ggplot2::geom_segment(
        data = growth_pred_edges_lineage,
        ggplot2::aes(
          x = .data$x_parent,
          y = .data$y_parent,
          xend = .data$x_passage,
          yend = .data$value
        ),
        color = "#1b9e77",
        linewidth = 0.8
      ) +
      ggplot2::geom_point(
        data = growth_pred_lineage,
        ggplot2::aes(x = .data$x_passage, y = .data$value),
        color = "#1b9e77",
        size = 1.3
      ) +
      ggplot2::geom_point(
        data = growth_obs_lineage,
        ggplot2::aes(x = .data$x_passage, y = .data$value),
        color = "#7570b3",
        size = 1.7,
        alpha = 0.82,
        shape = 17,
        position = ggplot2::position_jitter(width = 0.08, height = 0)
      ) +
      ggplot2::geom_hline(
        data = chr_ref_line,
        ggplot2::aes(yintercept = .data$yintercept),
        color = "black",
        linewidth = 0.28,
        alpha = 0.85
      ) +
      ggplot2::facet_grid(panel_label ~ lineage_panel_label, scales = "free_y") +
      ggplot2::scale_x_continuous(
        limits = c(x_lower - 0.25, x_upper + 0.25),
        breaks = x_breaks,
        labels = x_labels,
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      {
        if (requireNamespace("ggh4x", quietly = TRUE)) {
          ggh4x::facetted_pos_scales(y = y_scale_formulas())
        } else {
          NULL
        }
      } +
      ggplot2::scale_fill_manual(
        values = c("Viable" = "#1f77b4", "Hypoxia-Origin Dead" = "#d62728", "CIN-Associated Dead" = "#2ca02c"),
        drop = FALSE
      ) +
      ggplot2::labs(
        x = "Lineage passage / branch",
        y = NULL,
        fill = "Cell component"
      ) +
      theme_invitro() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
        strip.text = ggplot2::element_text(size = 9.2, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 32, hjust = 1, vjust = 1, size = 6.7, lineheight = 0.86),
        legend.position = if (isTRUE(show_legend)) "bottom" else "none",
        panel.spacing.x = grid::unit(0.65, "lines"),
        panel.spacing.y = grid::unit(0.55, "lines"),
        plot.margin = grid::unit(c(0.02, 0.06, 0.02, 0.06), "in")
      )
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) return(invisible(FALSE))
  make_cohort_plot <- function(cohort_value, show_legend = FALSE) {
    x_meta <- lineage_x_meta(cohort_value)
    if (!nrow(x_meta)) return(NULL)
    x_meta$width <- vapply(
      x_meta$lineage_label,
      function(lineage_value) global_lineage_width(lineage_value, x_meta$span[x_meta$lineage_label == lineage_value][[1]]),
      numeric(1)
    )
    lineage_plots <- Map(
      function(lineage_value, x_lower, x_upper, x_break_by, lineage_idx) {
        make_lineage_plot(
          cohort_value = cohort_value,
          lineage_value = lineage_value,
          x_lower = x_lower,
          x_upper = x_upper,
          x_break_by = x_break_by,
          show_legend = show_legend && lineage_idx == length(x_meta$lineage_label)
        )
      },
      x_meta$lineage_label,
      x_meta$x_lower,
      x_meta$x_upper,
      x_meta$x_break_by,
      seq_along(x_meta$lineage_label)
    )
    body_plot <- patchwork::wrap_plots(lineage_plots, nrow = 1, widths = x_meta$width)
    cohort_label_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0.5,
        label = as.character(cohort_value),
        hjust = 0,
        vjust = 0.5,
        fontface = "bold",
        size = 4.2
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0.06, -0.01, 0.06), "in"))
    patchwork::wrap_plots(cohort_label_plot, body_plot, ncol = 1, heights = c(0.055, 1))
  }
  cohort_plots <- Map(
    function(cohort_value, idx) make_cohort_plot(cohort_value, show_legend = idx == length(cohort_levels)),
    cohort_levels,
    seq_along(cohort_levels)
  )
  cohort_plots <- cohort_plots[!vapply(cohort_plots, is.null, logical(1))]
  if (length(cohort_plots) > 1L) {
    p <- patchwork::wrap_plots(cohort_plots, ncol = 1) +
      patchwork::plot_annotation(
        title = "Aligned In Vitro Growth, Chromosome-Number, and Burden Fits",
        subtitle = "Rows share each lineage x-axis; repeated lineage passages are split into branch-specific fixed-oxygen labels; growth and chromosome-number lines follow parent-child lineage links."
      )
  } else if (length(cohort_plots) == 1L) {
    p <- cohort_plots[[1]] +
      ggplot2::labs(
        title = "Aligned In Vitro Growth, Chromosome-Number, and Burden Fits",
        subtitle = "Rows share each lineage x-axis; repeated lineage passages are split into branch-specific fixed-oxygen labels; growth and chromosome-number lines follow parent-child lineage links."
      )
  } else {
    return(invisible(FALSE))
  }
  save_plot_pair(p, out_dir, "invitro_growth_ploidy_burden_composite", width = 12, height = 9)
  invisible(TRUE)
}

oxygen_label_values <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- num(vals)
  vals <- vals[is.finite(vals)]
  sort(unique(format(signif(vals, 3), trim = TRUE)))
}

invitro_branch_axis_ticks <- function(axis_map, cohort_value, lineage_value, x_break_by = 1) {
  if (is.null(axis_map) || !nrow(axis_map)) {
    return(data.frame(x_passage = numeric(), x_label = character(), stringsAsFactors = FALSE))
  }
  ticks <- axis_map[
    as.character(axis_map$cohort) == as.character(cohort_value) &
      as.character(axis_map$lineage_label) == as.character(lineage_value),
    ,
    drop = FALSE
  ]
  if (!nrow(ticks)) {
    return(data.frame(x_passage = numeric(), x_label = character(), stringsAsFactors = FALSE))
  }
  ticks <- ticks[order(ticks$x_passage_axis), , drop = FALSE]
  x_break_by <- suppressWarnings(as.numeric(x_break_by))
  if (!is.finite(x_break_by) || x_break_by < 1) x_break_by <- 1
  show_tick <- rep(TRUE, nrow(ticks))
  if (x_break_by > 1) {
    branch_tick <- as.logical(ticks$branch_duplicate)
    branch_tick[is.na(branch_tick)] <- FALSE
    show_tick <- abs(ticks$x_passage_axis %% x_break_by) < 1e-8 | branch_tick
    show_tick[is.na(show_tick)] <- FALSE
    if (!any(show_tick)) show_tick <- rep(TRUE, nrow(ticks))
  }
  data.frame(
    x_passage = ticks$x_passage_axis[show_tick],
    x_label = ticks$x_label_axis[show_tick],
    stringsAsFactors = FALSE
  )
}

invitro_branch_lineage_x_meta <- function(axis_map, cohort_value = NULL) {
  if (is.null(axis_map) || !nrow(axis_map)) return(data.frame())
  rows <- axis_map
  if (!is.null(cohort_value)) {
    rows <- rows[as.character(rows$cohort) == as.character(cohort_value), , drop = FALSE]
  }
  if (!nrow(rows)) return(data.frame())
  lineage_levels <- levels(order_invitro_lineage(rows$lineage_label))
  lineage_levels <- lineage_levels[lineage_levels %in% as.character(rows$lineage_label)]
  out <- lapply(lineage_levels, function(lineage_value) {
    lineage_x <- num(rows$x_passage_axis[as.character(rows$lineage_label) == lineage_value])
    lineage_x <- lineage_x[is.finite(lineage_x)]
    if (!length(lineage_x)) return(NULL)
    x_lower <- min(0, floor(min(lineage_x, na.rm = TRUE)))
    x_upper <- max(1, ceiling(max(lineage_x, na.rm = TRUE)))
    x_break_by <- if (x_upper > 25) 5 else if (x_upper > 12) 2 else 1
    data.frame(
      lineage_label = lineage_value,
      x_lower = x_lower,
      x_upper = x_upper,
      x_break_by = x_break_by,
      span = max(1, x_upper - x_lower),
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

make_invitro_branch_edges <- function(nodes) {
  required <- c("cohort", "lineage_label", "segment_id", "parent_segment_id", "x_passage", "value")
  if (is.null(nodes) || !nrow(nodes) || length(setdiff(required, names(nodes)))) return(data.frame())
  nodes <- nodes[is.finite(num(nodes$x_passage)) & is.finite(num(nodes$value)), , drop = FALSE]
  nodes$segment_id <- as.character(nodes$segment_id)
  nodes$parent_segment_id <- as.character(nodes$parent_segment_id)
  child <- nodes[!is.na(nodes$parent_segment_id) & nzchar(nodes$parent_segment_id), , drop = FALSE]
  if (!nrow(child)) return(data.frame())
  parent <- nodes[, c("cohort", "lineage_label", "segment_id", "x_passage", "value"), drop = FALSE]
  names(parent)[names(parent) == "segment_id"] <- "parent_segment_id"
  names(parent)[names(parent) == "x_passage"] <- "x_parent"
  names(parent)[names(parent) == "value"] <- "y_parent"
  edges <- merge(child, parent, by = c("cohort", "lineage_label", "parent_segment_id"), all = FALSE, sort = FALSE)
  if (!nrow(edges)) return(data.frame())
  edges$x_parent <- num(edges$x_parent)
  edges$y_parent <- num(edges$y_parent)
  edges$x_passage <- num(edges$x_passage)
  edges$value <- num(edges$value)
  edges <- edges[
    is.finite(edges$x_parent) & is.finite(edges$y_parent) &
      is.finite(edges$x_passage) & is.finite(edges$value),
    ,
    drop = FALSE
  ]
  edges$edge_id <- seq_len(nrow(edges))
  edges
}

build_branch_aware_o2_selected_live_plot <- function(daily_df, lineage_df, oxygen_levels = NULL) {
  if (is.null(lineage_df) || !all(c("cohort", "segment_id", "oxygen_pct", "predicted_live_cells") %in% names(lineage_df))) {
    return(NULL)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)

  axis_map <- build_invitro_branch_axis_map(lineage_df, daily_df)
  if (is.null(axis_map) || !nrow(axis_map)) return(NULL)
  if (is.null(oxygen_levels) || !length(oxygen_levels)) {
    oxygen_levels <- oxygen_label_values(axis_map$oxygen_pct, lineage_df$oxygen_pct)
  }
  if (!length(oxygen_levels)) return(NULL)

  axis_nodes <- axis_map
  axis_nodes$cohort <- order_invitro_cohort(axis_nodes$cohort)
  axis_nodes$lineage_label <- order_invitro_lineage(axis_nodes$lineage_label)
  axis_nodes$x_passage <- num(axis_nodes$x_passage_axis)
  axis_nodes$value <- num(axis_nodes$oxygen_pct)
  axis_nodes$oxygen_label <- format(signif(axis_nodes$value, 3), trim = TRUE)
  axis_nodes$oxygen_factor <- factor(axis_nodes$oxygen_label, levels = oxygen_levels)
  axis_nodes <- axis_nodes |>
    dplyr::filter(is.finite(.data$x_passage), is.finite(.data$value))
  oxygen_edges <- make_invitro_branch_edges(axis_nodes)

  selected_nodes <- ensure_invitro_plot_columns(lineage_df)
  selected_nodes$cohort <- order_invitro_cohort(selected_nodes$cohort)
  selected_nodes$lineage_label <- order_invitro_lineage(selected_nodes$lineage_label)
  selected_nodes <- attach_invitro_branch_axis(selected_nodes, axis_map)
  selected_nodes$x_passage <- num(selected_nodes$x_passage)
  selected_nodes$value <- num(selected_nodes$predicted_live_cells)
  selected_nodes$oxygen_pct_num <- num(selected_nodes$oxygen_pct)
  selected_nodes$oxygen_label <- format(signif(selected_nodes$oxygen_pct_num, 3), trim = TRUE)
  selected_nodes$oxygen_factor <- factor(selected_nodes$oxygen_label, levels = oxygen_levels)
  if (!"parent_segment_id" %in% names(selected_nodes)) selected_nodes$parent_segment_id <- NA_character_
  selected_nodes <- selected_nodes |>
    dplyr::filter(is.finite(.data$x_passage), is.finite(.data$value), .data$value > 0) |>
    dplyr::group_by(
      .data$cohort,
      .data$lineage_label,
      .data$segment_id,
      .data$parent_segment_id,
      .data$x_passage,
      .data$oxygen_factor
    ) |>
    dplyr::summarise(value = mean(.data$value, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(is.finite(.data$value), .data$value > 0)
  if (!nrow(selected_nodes)) return(NULL)
  selected_edges <- make_invitro_branch_edges(selected_nodes)

  observed_nodes <- data.frame()
  if ("observed_final_cells" %in% names(lineage_df)) {
    observed_nodes <- ensure_invitro_plot_columns(lineage_df)
    observed_nodes$cohort <- order_invitro_cohort(observed_nodes$cohort)
    observed_nodes$lineage_label <- order_invitro_lineage(observed_nodes$lineage_label)
    observed_nodes <- attach_invitro_branch_axis(observed_nodes, axis_map)
    observed_nodes$x_passage <- num(observed_nodes$x_passage)
    observed_nodes$value <- num(observed_nodes$observed_final_cells)
    observed_nodes <- observed_nodes |>
      dplyr::filter(is.finite(.data$x_passage), is.finite(.data$value), .data$value > 0) |>
      dplyr::group_by(
        .data$cohort,
        .data$lineage_label,
        .data$segment_id,
        .data$x_passage
      ) |>
      dplyr::summarise(value = mean(.data$value, na.rm = TRUE), .groups = "drop")
  }

  protocol_lines <- data.frame()
  if ("protocol_threshold_cells" %in% names(lineage_df)) {
    protocol_lines <- ensure_invitro_plot_columns(lineage_df)
    protocol_lines$cohort <- order_invitro_cohort(protocol_lines$cohort)
    protocol_lines$lineage_label <- order_invitro_lineage(protocol_lines$lineage_label)
    protocol_lines$value <- num(protocol_lines$protocol_threshold_cells)
    protocol_lines <- protocol_lines |>
      dplyr::filter(is.finite(.data$value), .data$value > 0) |>
      dplyr::group_by(.data$cohort, .data$lineage_label) |>
      dplyr::summarise(value = mean(.data$value, na.rm = TRUE), .groups = "drop")
  }

  cohort_levels <- levels(order_invitro_cohort(c(as.character(axis_nodes$cohort), as.character(selected_nodes$cohort))))
  cohort_levels <- cohort_levels[nzchar(cohort_levels)]
  if (!length(cohort_levels)) return(NULL)
  global_lineage_x_meta <- invitro_branch_lineage_x_meta(axis_map)
  global_lineage_width <- function(lineage_value, fallback) {
    idx <- which(as.character(global_lineage_x_meta$lineage_label) == as.character(lineage_value))
    if (length(idx) && is.finite(global_lineage_x_meta$span[idx[[1]]])) {
      return(global_lineage_x_meta$span[idx[[1]]])
    }
    fallback
  }
  subset_cohort_lineage <- function(df, cohort_value, lineage_value) {
    if (is.null(df) || !nrow(df) || !all(c("cohort", "lineage_label") %in% names(df))) return(df)
    df[
      as.character(df$cohort) == as.character(cohort_value) &
        as.character(df$lineage_label) == as.character(lineage_value),
      ,
      drop = FALSE
    ]
  }
  branch_markers_for <- function(cohort_value, lineage_value) {
    out <- axis_map[
      as.character(axis_map$cohort) == as.character(cohort_value) &
        as.character(axis_map$lineage_label) == as.character(lineage_value) &
        !is.na(axis_map$branch_duplicate) &
        as.logical(axis_map$branch_duplicate),
      ,
      drop = FALSE
    ]
    if (!nrow(out)) return(data.frame())
    data.frame(xintercept = num(out$x_passage_axis), stringsAsFactors = FALSE)
  }
  make_lineage_plot <- function(cohort_value, lineage_value, x_lower, x_upper, x_break_by, show_legend = FALSE) {
    o2_nodes <- subset_cohort_lineage(axis_nodes, cohort_value, lineage_value)
    o2_edges <- subset_cohort_lineage(oxygen_edges, cohort_value, lineage_value)
    live_nodes <- subset_cohort_lineage(selected_nodes, cohort_value, lineage_value)
    live_edges <- subset_cohort_lineage(selected_edges, cohort_value, lineage_value)
    observed_live_nodes <- subset_cohort_lineage(
      observed_nodes,
      cohort_value,
      lineage_value
    )
    protocol_line <- subset_cohort_lineage(
      protocol_lines,
      cohort_value,
      lineage_value
    )
    branch_markers <- branch_markers_for(cohort_value, lineage_value)
    axis_ticks <- invitro_branch_axis_ticks(axis_map, cohort_value, lineage_value, x_break_by)
    x_breaks <- if (nrow(axis_ticks)) axis_ticks$x_passage else seq(x_lower, x_upper, by = x_break_by)
    x_labels <- if (nrow(axis_ticks)) axis_ticks$x_label else x_breaks

    common_x <- ggplot2::scale_x_continuous(
      limits = c(x_lower - 0.25, x_upper + 0.25),
      breaks = x_breaks,
      labels = x_labels,
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    )
    branch_vline <- if (nrow(branch_markers)) {
      ggplot2::geom_vline(
        data = branch_markers,
        ggplot2::aes(xintercept = .data$xintercept),
        color = "black",
        linetype = "22",
        linewidth = 0.28,
        alpha = 0.45
      )
    } else {
      NULL
    }

    p_o2 <- ggplot2::ggplot() +
      branch_vline +
      ggplot2::geom_segment(
        data = o2_edges,
        ggplot2::aes(
          x = .data$x_parent,
          y = .data$y_parent,
          xend = .data$x_passage,
          yend = .data$value,
          color = .data$oxygen_factor
        ),
        linewidth = 0.65,
        alpha = 0.58
      ) +
      ggplot2::geom_point(
        data = o2_nodes,
        ggplot2::aes(x = .data$x_passage, y = .data$value, color = .data$oxygen_factor),
        size = 1.7
      ) +
      common_x +
      ggplot2::scale_color_viridis_d(option = "B", limits = oxygen_levels, drop = FALSE) +
      ggplot2::labs(title = as.character(lineage_value), x = NULL, y = "Fixed oxygen (%)", color = "Fixed oxygen (%)") +
      theme_invitro() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 10, hjust = 0.5),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none",
        plot.margin = grid::unit(c(0.02, 0.06, 0.01, 0.06), "in")
      )

    p_live <- ggplot2::ggplot() +
      branch_vline +
      if (nrow(protocol_line)) {
        ggplot2::geom_hline(
          data = protocol_line,
          ggplot2::aes(yintercept = .data$value),
          inherit.aes = FALSE,
          color = "grey35",
          linetype = "dashed",
          linewidth = 0.5
        )
      } else {
        NULL
      } +
      ggplot2::geom_segment(
        data = live_edges,
        ggplot2::aes(
          x = .data$x_parent,
          y = .data$y_parent,
          xend = .data$x_passage,
          yend = .data$value,
          color = .data$oxygen_factor
        ),
        linewidth = 0.75,
        alpha = 0.72
      ) +
      ggplot2::geom_point(
        data = live_nodes,
        ggplot2::aes(x = .data$x_passage, y = .data$value, color = .data$oxygen_factor),
        size = 1.75
      ) +
      if (nrow(observed_live_nodes)) {
        ggplot2::geom_point(
          data = observed_live_nodes,
          ggplot2::aes(x = .data$x_passage, y = .data$value),
          inherit.aes = FALSE,
          shape = 1,
          color = "black",
          stroke = 0.55,
          size = 2.2
        )
      } else {
        NULL
      } +
      common_x +
      ggplot2::scale_y_log10() +
      ggplot2::scale_color_viridis_d(option = "B", limits = oxygen_levels, drop = FALSE) +
      ggplot2::labs(x = "Lineage passage / branch", y = "Last-observation viable cells", color = "Fixed oxygen (%)") +
      theme_invitro() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 32, hjust = 1, vjust = 1, size = 6.7, lineheight = 0.86),
        legend.position = if (isTRUE(show_legend)) "bottom" else "none",
        plot.margin = grid::unit(c(0.01, 0.06, 0.02, 0.06), "in")
      )

    patchwork::wrap_plots(p_o2, p_live, ncol = 1, heights = c(0.46, 0.64))
  }
  make_cohort_plot <- function(cohort_value, show_legend = FALSE) {
    x_meta <- invitro_branch_lineage_x_meta(axis_map, cohort_value)
    if (!nrow(x_meta)) return(NULL)
    x_meta$width <- vapply(
      x_meta$lineage_label,
      function(lineage_value) global_lineage_width(lineage_value, x_meta$span[x_meta$lineage_label == lineage_value][[1]]),
      numeric(1)
    )
    lineage_plots <- Map(
      function(lineage_value, x_lower, x_upper, x_break_by, lineage_idx) {
        make_lineage_plot(
          cohort_value = cohort_value,
          lineage_value = lineage_value,
          x_lower = x_lower,
          x_upper = x_upper,
          x_break_by = x_break_by,
          show_legend = show_legend && lineage_idx == length(x_meta$lineage_label)
        )
      },
      x_meta$lineage_label,
      x_meta$x_lower,
      x_meta$x_upper,
      x_meta$x_break_by,
      seq_along(x_meta$lineage_label)
    )
    body_plot <- patchwork::wrap_plots(lineage_plots, nrow = 1, widths = x_meta$width)
    cohort_label_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0.5,
        label = as.character(cohort_value),
        hjust = 0,
        vjust = 0.5,
        fontface = "bold",
        size = 4.2
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0.06, -0.01, 0.06), "in"))
    patchwork::wrap_plots(cohort_label_plot, body_plot, ncol = 1, heights = c(0.065, 1))
  }
  cohort_plots <- Map(
    function(cohort_value, idx) make_cohort_plot(cohort_value, show_legend = idx == length(cohort_levels)),
    cohort_levels,
    seq_along(cohort_levels)
  )
  cohort_plots <- cohort_plots[!vapply(cohort_plots, is.null, logical(1))]
  if (!length(cohort_plots)) return(NULL)
  patchwork::wrap_plots(cohort_plots, ncol = 1) +
    patchwork::plot_annotation(
      title = "Independent-Lineage Fixed Oxygen and Passage-Endpoint Viable Cells",
      subtitle = paste(
        "Filled points are predicted counts at the last observation; open points are observed final counts;",
        "dashed lines are cohort-level T75 80% protocol thresholds."
      )
    )
}

build_constant_external_oxygen_plot <- function(daily_df, lineage_df, oxygen_levels = NULL) {
  build_branch_aware_o2_selected_live_plot(daily_df, lineage_df, oxygen_levels = oxygen_levels)
}

plot_constant_external_oxygen <- function(daily_df, lineage_df, out_dir) {
  p <- build_constant_external_oxygen_plot(daily_df, lineage_df)
  if (is.null(p)) return(invisible(FALSE))
  save_plot_pair(p, out_dir, "invitro_constant_external_oxygen", width = 14, height = 7.6)
  invisible(TRUE)
}

build_selected_day_live_cells_plot <- function(lineage_df, oxygen_levels = NULL) {
  build_branch_aware_o2_selected_live_plot(NULL, lineage_df, oxygen_levels = oxygen_levels)
}

plot_remote_o2_selected_live_panels <- function(daily_df, lineage_df, out_dir) {
  oxygen_levels <- oxygen_label_values(
    if (!is.null(daily_df) && "oxygen_pct" %in% names(daily_df)) daily_df$oxygen_pct else NULL,
    if (!is.null(lineage_df) && "oxygen_pct" %in% names(lineage_df)) lineage_df$oxygen_pct else NULL
  )
  p <- build_branch_aware_o2_selected_live_plot(daily_df, lineage_df, oxygen_levels = oxygen_levels)
  if (is.null(p)) return(invisible(FALSE))
  save_plot_pair(p, out_dir, "invitro_o2_selected_live_panels", width = 14, height = 7.6)
  invisible(TRUE)
}

plot_invitro_identifiability_from_tables <- function(variance_df, loadings_df, out_dir) {
  plot <- ivt_plot_identifiability_diagnostics(variance_df, loadings_df)
  if (is.null(plot)) return(invisible(FALSE))
  save_plot_pair(plot, out_dir, "invitro_identifiability_diagnostics", width = 11, height = 9)
  invisible(TRUE)
}

plot_remote_lineage_growth <- function(lineage_df, out_dir) {
  required <- c("predicted_growth_rate", "observed_growth", "passage_index", "cohort")
  if (is.null(lineage_df) || !all(required %in% names(lineage_df))) return(invisible(FALSE))
  df <- ensure_invitro_plot_columns(lineage_df)
  p <- ivt_plot_lineage_growth(df, primary_label = "Best fit")
  save_plot_pair(p, out_dir, "invitro_lineage_growth", width = 10, height = 6.8)
  invisible(TRUE)
}

plot_remote_lineage_ploidy <- function(lineage_df, quantile_df, observed_kary_df, out_dir) {
  if (is.null(lineage_df) || is.null(quantile_df)) return(invisible(FALSE))
  if (!all(c("predicted_quantile_kary_N", "quantile_prob") %in% names(quantile_df))) return(invisible(FALSE))
  scenario_observed_kary_df <- .invitro_scenario_kary_rows(observed_kary_df)
  if (!is.null(scenario_observed_kary_df) &&
      !"observed_kary_N" %in% names(scenario_observed_kary_df)) {
    return(invisible(FALSE))
  }
  p <- ivt_plot_lineage_ploidy(
    ensure_invitro_plot_columns(lineage_df),
    quantile_df = ensure_invitro_plot_columns(quantile_df),
    observed_kary_df = ensure_invitro_plot_columns(scenario_observed_kary_df),
    primary_label = "Best fit",
    quantile_alpha = 0.5
  )
  save_plot_pair(p, out_dir, "invitro_lineage_ploidy", width = 10, height = 6.8)
  invisible(TRUE)
}

plot_remote_flow_density <- function(flow_df, out_dir) {
  if (is.null(flow_df) || !all(c("ploidy", "density", "series", "cohort") %in% names(flow_df))) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No matched flow-density observations available")
    save_plot_pair(p, out_dir, "invitro_flow_density", width = 10, height = 6.8)
    return(invisible(TRUE))
  }
  p <- ivt_plot_lineage_flow_density(
    ensure_invitro_plot_columns(flow_df),
    max_facets = 20L
  )
  save_plot_pair(p, out_dir, "invitro_flow_density", width = 10, height = 6.8)
  invisible(TRUE)
}

plot_remote_distribution_heatmap <- function(dist_df, out_dir) {
  if (is.null(dist_df) || !all(c("N", "fraction", "passage_index", "cohort") %in% names(dist_df))) {
    return(invisible(FALSE))
  }
  p <- ivt_plot_distribution_heatmap(ensure_invitro_plot_columns(dist_df), max_N = 110)
  save_plot_pair(p, out_dir, "invitro_distribution_heatmap", width = 10, height = 6.5)
  invisible(TRUE)
}

plot_remote_loglik_by_passage <- function(df,
                                          out_dir,
                                          basename,
                                          value_col,
                                          title,
                                          y_label) {
  if (is.null(df) || !all(c("passage_index", value_col, "cohort", "segment_id") %in% names(df))) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste0("No ", tolower(gsub("In vitro ", "", title)), " rows"))
    save_plot_pair(p, out_dir, basename, width = 9.5, height = 5)
    return(invisible(TRUE))
  }
  plot_df <- ensure_invitro_plot_columns(df)
  plot_df$passage_index <- num(plot_df$passage_index)
  plot_df[[value_col]] <- num(plot_df[[value_col]])
  plot_df <- finite_rows(plot_df, c("passage_index", value_col))
  if (nrow(plot_df) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste0("No ", tolower(gsub("In vitro ", "", title)), " rows"))
    save_plot_pair(p, out_dir, basename, width = 9.5, height = 5)
    return(invisible(TRUE))
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(passage_index, .data[[value_col]], color = cohort)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(cohort, segment_id)), alpha = 0.35) +
    ggplot2::labs(
      title = title,
      x = "Passage index",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 12)
  save_plot_pair(p, out_dir, basename, width = 9.5, height = 5)
  invisible(TRUE)
}

theme_invitro <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey75"),
      legend.position = "bottom"
    )
}

.invitro_objective_component_metrics <- function(summary_df) {
  value_for <- function(metric) {
    direct <- suppressWarnings(as.numeric(summary_value(summary_df, metric)))
    if (length(direct) == 1L && is.finite(direct)) return(direct)
    suppressWarnings(as.numeric(summary_value(summary_df, paste0("invitro_", metric))))
  }
  first_finite <- function(...) {
    values <- suppressWarnings(as.numeric(unlist(list(...), use.names = FALSE)))
    values <- values[is.finite(values)]
    if (length(values)) values[[1L]] else NA_real_
  }
  modalities <- c("growth", "ploidy", "flow")
  modality_labels <- c("Growth", "Karyotype", "Flow")
  if (any(c("death_loglik", "invitro_death_loglik") %in% as.character(summary_df$metric))) {
    modalities <- c(modalities, "death")
    modality_labels <- c(modality_labels, "Death fraction")
  }
  if (any(c("passage_time_loglik", "invitro_passage_time_loglik") %in% as.character(summary_df$metric))) {
    modalities <- c(modalities, "passage_time")
    modality_labels <- c(modality_labels, "Passage time")
  }
  hierarchical_loglik <- vapply(
    paste0(modalities, "_loglik"),
    value_for,
    numeric(1)
  )
  weights <- vapply(modalities, function(modality) {
    first_finite(
      value_for(paste0(modality, "_weight")),
      value_for(paste0("joint_invitro_", modality, "_weight")),
      1
    )
  }, numeric(1))

  if (all(is.finite(hierarchical_loglik)) && all(is.finite(weights))) {
    contributions <- -weights * hierarchical_loglik
    component_labels <- paste0(
      "- ", modality_labels, " hierarchical logLik (weight=",
      format(weights, trim = TRUE), ")"
    )
    prior_penalty <- value_for("buffer_prior_penalty")
    if (length(prior_penalty) == 1L && is.finite(prior_penalty)) {
      contributions <- c(contributions, prior_penalty)
      component_labels <- c(component_labels, "Buffer soft-prior penalty")
    }
    total_objective <- first_finite(
      value_for("objective_total"),
      value_for("objective_invitro"),
      sum(contributions)
    )
    out <- data.frame(
      component = c(
        "Total objective",
        component_labels
      ),
      value = c(total_objective, contributions),
      stringsAsFactors = FALSE
    )
    attr(out, "title") <- "In Vitro Hierarchical Objective Components"
    attr(out, "y_label") <- "Weighted objective contribution"
    attr(out, "scale_mode") <- "hierarchical_objective"
    return(out)
  }

  raw_sums <- vapply(
    paste0(modalities, "_loglik_sum"),
    value_for,
    numeric(1)
  )
  out <- data.frame(
    component = paste0("- ", modality_labels, " raw logLik sum"),
    value = -raw_sums,
    stringsAsFactors = FALSE
  )
  attr(out, "title") <- "In Vitro Raw Log-Likelihood Sums (Legacy Diagnostic)"
  attr(out, "y_label") <- "Raw negative log-likelihood sum (not objective scale)"
  attr(out, "scale_mode") <- "legacy_raw_sum"
  out
}

plot_objective_components <- function(summary_df, out_dir) {
  if (is.null(summary_df)) return(invisible(FALSE))
  metrics <- .invitro_objective_component_metrics(summary_df)
  plot_title <- attr(metrics, "title")
  y_label <- attr(metrics, "y_label")
  metrics <- metrics[is.finite(metrics$value), , drop = FALSE]
  if (nrow(metrics) == 0L) return(invisible(FALSE))
  metrics$component <- factor(metrics$component, levels = rev(metrics$component))
  p <- ggplot2::ggplot(metrics, ggplot2::aes(x = component, y = value)) +
    ggplot2::geom_col(fill = "#4B6F8A", width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = plot_title,
      x = NULL,
      y = y_label
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_objective_components", width = 8.5, height = 4.2)
  invisible(TRUE)
}

plot_growth_rate_fit <- function(growth_df, out_dir) {
  if (is.null(growth_df) || !all(c("observed_growth", "predicted_growth_rate") %in% names(growth_df))) {
    return(invisible(FALSE))
  }
  df <- growth_df
  df$observed <- num(df$observed_growth)
  df$predicted <- num(df$predicted_growth_rate)
  df <- finite_rows(df, c("observed", "predicted"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  axis_range <- range(c(df$observed, df$predicted), finite = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = predicted, colour = cohort)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_point(alpha = 0.8, size = 2.2) +
    ggplot2::coord_equal(xlim = axis_range, ylim = axis_range) +
    ggplot2::labs(
      title = "In Vitro Growth-Rate Fit",
      x = "Observed growth rate",
      y = "Predicted growth rate",
      colour = "Cohort"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_growth_rate_fit", width = 7.2, height = 6)
  invisible(TRUE)
}

plot_death_fraction_fit <- function(death_df, out_dir) {
  required <- c("observed_dead_fraction", "predicted_dead_fraction", "cohort", "lineage_id")
  if (is.null(death_df) || !all(required %in% names(death_df))) return(invisible(FALSE))
  df <- death_df
  df$observed <- num(df$observed_dead_fraction)
  df$predicted <- num(df$predicted_dead_fraction)
  df <- finite_rows(df, c("observed", "predicted"))
  df <- df[df$observed > 0 & df$predicted > 0, , drop = FALSE]
  if (!nrow(df)) return(invisible(FALSE))
  limits <- range(c(df$observed, df$predicted), finite = TRUE)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = observed, y = predicted, colour = cohort, shape = lineage_id)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_point(alpha = 0.82, size = 2.3) +
    ggplot2::scale_x_log10(limits = limits) +
    ggplot2::scale_y_log10(limits = limits) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "In Vitro Death-Fraction Fit",
      x = "Observed dead-cell fraction (log10 scale)",
      y = "Predicted dead-cell fraction (log10 scale)",
      colour = "Cohort",
      shape = "Lineage"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_death_fraction_fit", width = 7.4, height = 6.2)
  invisible(TRUE)
}

plot_death_logit_residual <- function(death_df, out_dir) {
  required <- c("lineage_passage_index", "logit_residual", "cohort", "lineage_id")
  if (is.null(death_df) || !all(required %in% names(death_df))) return(invisible(FALSE))
  df <- death_df
  df$passage <- num(df$lineage_passage_index)
  df$residual <- num(df$logit_residual)
  df <- finite_rows(df, c("passage", "residual"))
  if (!nrow(df)) return(invisible(FALSE))
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = passage, y = residual, colour = cohort, shape = lineage_id)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_line(ggplot2::aes(group = interaction(cohort, lineage_id)), alpha = 0.45) +
    ggplot2::geom_point(alpha = 0.82, size = 2.1) +
    ggplot2::facet_wrap(~cohort, scales = "free_x") +
    ggplot2::labs(
      title = "In Vitro Death-Fraction Logit Residuals",
      x = "Passage index",
      y = "logit(observed) - logit(predicted)",
      colour = "Cohort",
      shape = "Lineage"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_death_logit_residual", width = 9.2, height = 5.4)
  invisible(TRUE)
}

plot_growth_count_fit <- function(growth_df, out_dir) {
  required <- c(
    "observed_live_cells_at_observation",
    "predicted_live_cells_at_observation"
  )
  if (is.null(growth_df) || !all(required %in% names(growth_df))) {
    return(invisible(FALSE))
  }
  df <- growth_df
  df$observed_cells <- num(df$observed_live_cells_at_observation)
  df$predicted_cells <- num(df$predicted_live_cells_at_observation)
  df <- finite_rows(df, c("observed_cells", "predicted_cells"))
  df <- df[df$observed_cells > 0 & df$predicted_cells > 0, , drop = FALSE]
  if (nrow(df) == 0L) return(invisible(FALSE))
  df$log10_observed_cells <- log10(df$observed_cells)
  df$log10_predicted_cells <- log10(df$predicted_cells)
  axis_range <- range(c(df$log10_observed_cells, df$log10_predicted_cells), finite = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log10_observed_cells, y = log10_predicted_cells, colour = cohort)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_point(alpha = 0.75, size = 2) +
    ggplot2::coord_equal(xlim = axis_range, ylim = axis_range) +
    ggplot2::labs(
      title = "In Vitro Measured-Timepoint Viable-Cell Count Fit",
      x = "Observed viable cells at measured timepoints (log10)",
      y = "Predicted viable cells at measured timepoints (log10)",
      colour = "Cohort"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_live_count_fit", width = 7.2, height = 6)
  invisible(TRUE)
}

plot_ploidy_mean_fit <- function(lineage_df, out_dir) {
  if (is.null(lineage_df) || !all(c("observed_mean_kary_N", "predicted_mean_kary_N") %in% names(lineage_df))) {
    return(invisible(FALSE))
  }
  df <- lineage_df
  df$observed <- num(df$observed_mean_kary_N)
  df$predicted <- num(df$predicted_mean_kary_N)
  df <- finite_rows(df, c("observed", "predicted"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  axis_range <- range(c(df$observed, df$predicted), finite = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = observed, y = predicted, colour = cohort)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_point(alpha = 0.85, size = 2.4) +
    ggplot2::coord_equal(xlim = axis_range, ylim = axis_range) +
    ggplot2::labs(
      title = "In Vitro Mean Karyotype Fit",
      x = "Observed mean chromosome number (N)",
      y = "Predicted mean chromosome number (N)",
      colour = "Cohort"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_mean_karyotype_fit", width = 7.2, height = 6)
  invisible(TRUE)
}

plot_ploidy_loglik <- function(ploidy_df, out_dir) {
  if (is.null(ploidy_df) || !all(c("passage_index", "mean_loglik") %in% names(ploidy_df))) {
    return(invisible(FALSE))
  }
  df <- ploidy_df
  df$passage <- num(df$passage_index)
  df$mean_loglik_num <- num(df$mean_loglik)
  df <- finite_rows(df, c("passage", "mean_loglik_num"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = passage, y = mean_loglik_num, colour = cohort)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::labs(
      title = "In Vitro Karyotype Log-Likelihood by Passage",
      x = "Passage index",
      y = "Mean log-likelihood per observed cell",
      colour = "Cohort"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_karyotype_loglik", width = 8.5, height = 5.2)
  invisible(TRUE)
}

plot_distribution_quantiles <- function(quantile_df, out_dir) {
  required <- c("cohort", "passage_index", "quantile_prob", "predicted_quantile_kary_N")
  if (is.null(quantile_df) || !all(required %in% names(quantile_df))) return(invisible(FALSE))
  df <- quantile_df
  df$passage <- num(df$passage_index)
  df$quantile <- num(df$quantile_prob)
  df$predicted <- num(df$predicted_quantile_kary_N)
  df <- finite_rows(df, c("passage", "quantile", "predicted"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  agg <- stats::aggregate(
    predicted ~ cohort + passage + quantile,
    data = df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  agg$quantile <- factor(agg$quantile)
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = passage, y = predicted, colour = quantile)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~cohort, scales = "free_x") +
    ggplot2::labs(
      title = "Predicted In Vitro Karyotype Quantiles",
      x = "Passage index",
      y = "Predicted chromosome number (N)",
      colour = "Quantile"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_karyotype_quantiles", width = 9.5, height = 5.4)
  invisible(TRUE)
}

plot_distribution_heatmap <- function(dist_df, out_dir) {
  required <- c("cohort", "passage_index", "N", "fraction")
  if (is.null(dist_df) || !all(required %in% names(dist_df))) return(invisible(FALSE))
  df <- dist_df
  df$passage <- num(df$passage_index)
  df$kary_N <- num(df$N)
  df$fraction_num <- num(df$fraction)
  df <- finite_rows(df, c("passage", "kary_N", "fraction_num"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  agg <- stats::aggregate(
    fraction_num ~ cohort + passage + kary_N,
    data = df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  fraction_fill_max <- max(agg$fraction_num, na.rm = TRUE)
  if (!is.finite(fraction_fill_max) || fraction_fill_max <= 0) {
    fraction_fill_max <- 1
  }
  p <- ggplot2::ggplot(agg, ggplot2::aes(x = passage, y = kary_N, fill = fraction_num)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~cohort, scales = "free_x") +
    ivt_ploidy_fraction_fill_scale(fraction_fill_max, name = "Mean fraction") +
    ggplot2::labs(
      title = "Predicted In Vitro Karyotype Distribution",
      x = "Passage index",
      y = "Chromosome number (N)",
      fill = "Mean fraction"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_karyotype_distribution", width = 9.5, height = 5.4)
  invisible(TRUE)
}

plot_invitro_missegregation_from_table <- function(timecourse_df, out_dir) {
  plot <- ivt_plot_missegregation_probability_timecourse(timecourse_df)
  if (is.null(plot)) return(invisible(FALSE))
  save_plot_pair(
    plot,
    out_dir,
    "invitro_missegregation_probability_over_passage",
    width = 11,
    height = 6.6
  )
  invisible(TRUE)
}

plot_daily_live_cells <- function(daily_df, out_dir) {
  required <- c("cohort", "segment_id", "day", "live_cells")
  if (is.null(daily_df) || !all(required %in% names(daily_df))) return(invisible(FALSE))
  df <- daily_df
  df$day_num <- num(df$day)
  df$live_cells_num <- num(df$live_cells)
  df$oxygen_pct_num <- num(df$oxygen_pct)
  df <- finite_rows(df, c("day_num", "live_cells_num"))
  df <- df[df$live_cells_num > 0, , drop = FALSE]
  if (nrow(df) == 0L) return(invisible(FALSE))
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = day_num,
      y = live_cells_num,
      group = segment_id,
      colour = oxygen_pct_num
    )
  ) +
    ggplot2::geom_line(alpha = 0.35, linewidth = 0.45) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_viridis_c(option = "B", na.value = "grey45") +
    ggplot2::facet_wrap(~cohort, scales = "free_x") +
    ggplot2::labs(
      title = "Predicted In Vitro Viable-Cell Time Courses",
      x = "Day within passage",
      y = "Viable cells (log10)",
      colour = "Fixed oxygen (%)"
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_daily_live_cells", width = 10, height = 5.8)
  invisible(TRUE)
}

plot_flow_overlay <- function(flow_df, out_dir) {
  if (is.null(flow_df) || nrow(flow_df) == 0L) return(invisible(FALSE))
  names_lower <- tolower(names(flow_df))
  x_col <- names(flow_df)[match(TRUE, names_lower %in% c("ploidy", "ploidy_n", "kary_n", "n", "x"))]
  pred_col <- names(flow_df)[match(TRUE, names_lower %in% c("predicted_density", "pred_density", "model_density"))]
  obs_col <- names(flow_df)[match(TRUE, names_lower %in% c("observed_density", "obs_density", "density"))]
  if (is.na(x_col) || is.na(pred_col) || is.na(obs_col)) return(invisible(FALSE))
  df <- flow_df
  df$x <- num(df[[x_col]])
  df$predicted <- num(df[[pred_col]])
  df$observed <- num(df[[obs_col]])
  df <- finite_rows(df, c("x", "predicted", "observed"))
  if (nrow(df) == 0L) return(invisible(FALSE))
  long <- rbind(
    data.frame(df[, setdiff(names(df), c("predicted", "observed")), drop = FALSE], density = df$predicted, series = "Predicted"),
    data.frame(df[, setdiff(names(df), c("predicted", "observed")), drop = FALSE], density = df$observed, series = "Observed")
  )
  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = density, colour = series)) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
    ggplot2::facet_wrap(~cohort, scales = "free_y") +
    ggplot2::labs(
      title = "In Vitro Flow-Density Overlay",
      x = x_col,
      y = "Density",
      colour = NULL
    ) +
    theme_invitro()
  save_plot_pair(p, out_dir, "invitro_flow_overlay", width = 9.5, height = 5.4)
  invisible(TRUE)
}

plot_invitro_functional_response_from_tables <- function(o2_curve_df,
                                                           viability_df,
                                                           ploidy_o2_df,
                                                           out_dir) {
  plot <- ivt_plot_functional_response_diagnostics(
    o2_curve_df,
    viability_df,
    ploidy_o2_df
  )
  if (is.null(plot)) return(invisible(FALSE))
  save_plot_pair(plot, out_dir, "invitro_rate_function_diagnostics", width = 14, height = 8.8)
  invisible(TRUE)
}

read_preferred_tsv <- function(primary_dir, fallback_dir, basename) {
  primary <- file.path(primary_dir, basename)
  if (file.exists(primary)) return(read_tsv_optional(primary))
  read_tsv_optional(file.path(fallback_dir, basename))
}

write_manifest <- function(out_dir, fit_dir, simulation_dir, analysis_dir, generated) {
  plot_files <- sort(basename(list.files(
    out_dir,
    pattern = "\\.(pdf|png|svg)$",
    full.names = FALSE
  )))
  manifest <- data.frame(
    key = c(
      "fit_dir", "simulation_dir", "analysis_dir", "out_dir", "generated_at",
      "plot_file_count", names(generated), paste0("plot.", seq_along(plot_files))
    ),
    value = c(
      normalizePath(fit_dir, mustWork = FALSE),
      normalizePath(simulation_dir, mustWork = FALSE),
      normalizePath(analysis_dir, mustWork = FALSE),
      normalizePath(out_dir, mustWork = FALSE),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      as.character(length(plot_files)),
      as.character(unlist(generated, use.names = FALSE)),
      plot_files
    ),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    manifest,
    file = file.path(out_dir, "invitro_viz_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir %||% stop(
    "Usage: viz_invitro_model_O2_supply_demand_MAP_results.R --fit_dir=/abs/path/to/seed_dir",
    call. = FALSE
  )
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  simulation_dir <- argv$simulation_dir %||% file.path(fit_dir, "simulation", "invitro")
  analysis_dir <- argv$analysis_dir %||% file.path(fit_dir, "analysis", "invitro")
  out_dir <- argv$out_dir %||% file.path(fit_dir, "viz", "invitro")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  unlink(list.files(out_dir, pattern = "\\.(pdf|png|svg)$", full.names = TRUE), force = TRUE)

  lineage_df <- read_preferred_tsv(simulation_dir, fit_dir, "invitro_lineage_summary.tsv")
  dist_df <- read_preferred_tsv(simulation_dir, fit_dir, "invitro_distribution_summary.tsv")
  quantile_df <- read_preferred_tsv(simulation_dir, fit_dir, "invitro_distribution_quantiles.tsv")
  daily_df <- read_preferred_tsv(simulation_dir, fit_dir, "invitro_daily_counts.tsv")
  flow_df <- read_tsv_optional(file.path(fit_dir, "invitro_flow_overlay.tsv"))
  observed_kary_df <- read_tsv_optional(file.path(fit_dir, "invitro_observed_kary.tsv"))
  growth_df <- read_tsv_optional(file.path(fit_dir, "invitro_growth_loglik.tsv"))
  death_df <- read_tsv_optional(file.path(fit_dir, "invitro_death_loglik.tsv"))
  summary_df <- read_tsv_optional(file.path(fit_dir, "fit_summary.tsv"))
  misseg_df <- read_tsv_optional(file.path(
    simulation_dir,
    "invitro_missegregation_probability_timecourse.tsv"
  ))
  functional_o2_df <- read_tsv_optional(file.path(
    simulation_dir,
    "functional_curve_oxygen_multi_ploidy.tsv"
  ))
  functional_viability_df <- read_tsv_optional(file.path(
    simulation_dir,
    "functional_curve_ploidy.tsv"
  ))
  functional_ploidy_o2_df <- read_tsv_optional(file.path(
    simulation_dir,
    "functional_curve_ploidy_by_o2.tsv"
  ))
  identifiability_variance_df <- read_tsv_optional(file.path(
    analysis_dir,
    "invitro_identifiability_variance.tsv"
  ))
  identifiability_loadings_df <- read_tsv_optional(file.path(
    analysis_dir,
    "invitro_identifiability_loadings.tsv"
  ))

  generated <- list(
    objective_components = plot_objective_components(summary_df, out_dir),
    growth_count_fit = plot_growth_count_fit(growth_df, out_dir),
    death_fraction_fit = plot_death_fraction_fit(death_df, out_dir),
    death_logit_residual = plot_death_logit_residual(death_df, out_dir),
    identifiability_diagnostics = plot_invitro_identifiability_from_tables(
      identifiability_variance_df,
      identifiability_loadings_df,
      out_dir
    ),
    o2_selected_live_panels = plot_remote_o2_selected_live_panels(daily_df, lineage_df, out_dir),
    missegregation_probability = plot_invitro_missegregation_from_table(misseg_df, out_dir),
    rate_function_diagnostics = plot_invitro_functional_response_from_tables(
      functional_o2_df,
      functional_viability_df,
      functional_ploidy_o2_df,
      out_dir
    ),
    daily_counts = plot_remote_daily_counts(daily_df, out_dir),
    lineage_growth = plot_remote_lineage_growth(lineage_df, out_dir),
    lineage_ploidy = plot_remote_lineage_ploidy(lineage_df, quantile_df, observed_kary_df, out_dir),
    burden_decomposition = plot_remote_burden_decomposition(daily_df, out_dir),
    growth_ploidy_burden_composite = plot_remote_growth_ploidy_burden_composite(lineage_df, quantile_df, observed_kary_df, daily_df, out_dir),
    flow_density = plot_remote_flow_density(flow_df, out_dir),
    distribution_heatmap = plot_remote_distribution_heatmap(dist_df, out_dir)
  )
  write_manifest(out_dir, fit_dir, simulation_dir, analysis_dir, generated)
  message("In vitro viz written to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(normalizePath(out_dir, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  main()
}
