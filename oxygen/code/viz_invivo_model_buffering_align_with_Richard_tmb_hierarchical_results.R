#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(argv) {
  out <- list()
  if (length(argv) == 0) return(out)
  for (a in argv) {
    if (!startsWith(a, "--")) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.numeric(x))
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.integer(x))
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

get_script_dir_self <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grepl("^--file=", args)]
  if (length(farg) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", farg[[1]])))
}

find_latest_fit_dir <- function(results_root) {
  dirs <- list.dirs(results_root, recursive = FALSE, full.names = TRUE)
  if (length(dirs) == 0) stop("No result directories under: ", results_root)
  dirs[[which.max(file.info(dirs)$mtime)]]
}

read_tsv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

step_dirs_under <- function(run_dir) {
  dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[basename(dirs) != "viz_tmb"]
  keep <- grepl("^step[0-9]+_wb", basename(dirs))
  sort(dirs[keep])
}

is_tmb_hier_run_dir <- function(d) {
  file.exists(file.path(d, "fit_config.rds")) &&
    length(step_dirs_under(d)) > 0
}

find_tmb_run_dirs_under <- function(root_dir) {
  all_dirs <- unique(c(root_dir, list.dirs(root_dir, recursive = TRUE, full.names = TRUE)))
  runs <- all_dirs[vapply(all_dirs, is_tmb_hier_run_dir, logical(1))]
  sort(unique(normalizePath(runs, mustWork = TRUE)))
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

count_na_cells <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(0L)
  total <- 0L
  for (nm in names(df)) {
    x <- df[[nm]]
    total <- total + sum(is.na(x))
    if (is.character(x)) total <- total + sum(x == "NA", na.rm = TRUE)
  }
  as.integer(total)
}

collect_step_na_summary <- function(step_dir) {
  tsvs <- sort(list.files(step_dir, pattern = "\\.tsv$", full.names = TRUE))
  if (length(tsvs) == 0) return(data.frame())
  rows <- lapply(tsvs, function(f) {
    df <- read_tsv_if_exists(f)
    data.frame(
      file = basename(f),
      rows = if (is.null(df)) 0L else nrow(df),
      cols = if (is.null(df)) 0L else ncol(df),
      na_cells = count_na_cells(df),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

pick_selected_step <- function(run_dir, argv) {
  step_dirs <- step_dirs_under(run_dir)
  if (length(step_dirs) == 0) stop("No step directories in: ", run_dir)

  if (!is.null(argv$step_dir)) {
    s <- argv$step_dir
    cand <- if (grepl("^/", s)) s else file.path(run_dir, s)
    if (!dir.exists(cand)) stop("Requested --step_dir not found: ", cand)
    return(normalizePath(cand))
  }

  if (!is.null(argv$step_name)) {
    cand <- file.path(run_dir, argv$step_name)
    if (!dir.exists(cand)) stop("Requested --step_name not found: ", cand)
    return(normalizePath(cand))
  }

  sel_file <- file.path(run_dir, "selected_best_step.tsv")
  if (file.exists(sel_file)) {
    sel <- read.delim(sel_file, check.names = FALSE, stringsAsFactors = FALSE)
    if ("step_dir" %in% names(sel) && nrow(sel) >= 1) {
      s <- as.character(sel$step_dir[[1]])
      if (nzchar(s) && dir.exists(s)) return(normalizePath(s))
    }
    if (all(c("pass", "w_burden", "w_ploidy") %in% names(sel)) && nrow(sel) >= 1) {
      pass_i <- as.integer(sel$pass[[1]])
      wb <- as.character(sel$w_burden[[1]])
      wp <- as.character(sel$w_ploidy[[1]])
      cand_pat <- paste0("^step", sprintf("%02d", pass_i), "_wb")
      cand <- step_dirs[grepl(cand_pat, basename(step_dirs))]
      if (length(cand) >= 1) return(normalizePath(cand[[1]]))
    }
  }

  chain_file <- file.path(run_dir, "chain_global_summary.tsv")
  if (file.exists(chain_file)) {
    ch <- read.delim(chain_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(ch) > 0 && "step_dir" %in% names(ch) && "objective_data" %in% names(ch)) {
      ord <- order(safe_numeric(ch$objective_data), na.last = TRUE)
      for (i in ord) {
        s <- as.character(ch$step_dir[[i]])
        if (nzchar(s) && dir.exists(s)) return(normalizePath(s))
      }
    }
  }

  step_dirs[[which.max(file.info(step_dirs)$mtime)]]
}

make_chain_plots <- function(run_dir, out_dir) {
  chain_file <- file.path(run_dir, "chain_global_summary.tsv")
  if (!file.exists(chain_file)) return(invisible(NULL))
  ch <- read.delim(chain_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(ch) == 0) return(invisible(NULL))

  num_cols <- intersect(
    c("pass", "w_burden", "w_ploidy", "objective_data", "objective_penalized",
      "objective_burden_scaled", "objective_ploidy_scaled", "global_active_workers",
      "global_requested_workers", "global_started_workers"),
    names(ch)
  )
  for (nm in num_cols) ch[[nm]] <- safe_numeric(ch[[nm]])
  if ("global_fallback" %in% names(ch)) ch$global_fallback <- as.logical(ch$global_fallback)

  if (all(c("w_burden", "w_ploidy", "pass") %in% names(ch))) {
    p_weight <- ggplot(ch, aes(x = w_burden, y = w_ploidy)) +
      geom_path(color = "#377eb8", linewidth = 0.7, alpha = 0.8) +
      geom_point(aes(color = objective_data), size = 2) +
      geom_text(aes(label = pass), nudge_y = 0.02, size = 3, check_overlap = TRUE) +
      scale_color_viridis_c(option = "C", na.value = "grey70") +
      labs(
        title = "TMB Hierarchical Chain Weight Path",
        subtitle = basename(run_dir),
        x = "w_burden",
        y = "w_ploidy",
        color = "objective_data"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "chain_weight_path.pdf"), p_weight, width = 8, height = 6)
  }

  metric_cols <- intersect(c("objective_data", "objective_burden_scaled", "objective_ploidy_scaled"), names(ch))
  if (length(metric_cols) > 0 && "pass" %in% names(ch)) {
    long <- ch %>%
      select(any_of(c("pass", metric_cols))) %>%
      pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value")
    p_obj <- ggplot(long, aes(x = pass, y = value, color = metric)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.6) +
      labs(
        title = "Chain Objectives Across Weight Steps",
        subtitle = basename(run_dir),
        x = "Pass",
        y = "Value",
        color = "Metric"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "chain_objectives.pdf"), p_obj, width = 10, height = 6)
  }

  if (all(c("objective_burden_scaled", "objective_ploidy_scaled", "pass") %in% names(ch))) {
    p_tradeoff <- ggplot(ch, aes(x = objective_burden_scaled, y = objective_ploidy_scaled)) +
      geom_path(color = "grey60", linewidth = 0.7) +
      geom_point(aes(color = pass), size = 2.2) +
      geom_text(aes(label = pass), nudge_y = 0.01, size = 3, check_overlap = TRUE) +
      scale_color_viridis_c(option = "D") +
      labs(
        title = "Chain Tradeoff (Scaled Burden vs Scaled Ploidy)",
        subtitle = basename(run_dir),
        x = "objective_burden_scaled",
        y = "objective_ploidy_scaled",
        color = "Pass"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "chain_tradeoff.pdf"), p_tradeoff, width = 8, height = 6)
  }

  if (all(c("pass", "global_active_workers") %in% names(ch))) {
    if (!("global_fallback" %in% names(ch))) ch$global_fallback <- FALSE
    p_workers <- ggplot(ch, aes(x = pass, y = global_active_workers)) +
      geom_line(color = "#1b9e77", linewidth = 0.8) +
      geom_point(aes(shape = global_fallback), size = 2) +
      labs(
        title = "Global DEoptim Worker Usage Across Chain",
        subtitle = basename(run_dir),
        x = "Pass",
        y = "global_active_workers",
        shape = "global_fallback"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "chain_global_workers.pdf"), p_workers, width = 9, height = 5)
  }
}

plot_burden_fit <- function(step_dir, out_dir) {
  f <- file.path(step_dir, "global_burden_fit.tsv")
  if (!file.exists(f)) return(invisible(NULL))
  df <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible(NULL))
  for (nm in intersect(c("day", "pred_pop", "pred_burden_volume_mm3", "obs_burden", "obs_norm", "pred_norm", "dose"), names(df))) {
    df[[nm]] <- safe_numeric(df[[nm]])
  }

  if (!("scenario_id" %in% names(df))) {
    df$scenario_id <- paste(df$harvest %||% "", df$cohort %||% "", df$dose %||% "", sep = "__")
  }

  if (all(c("pred_norm", "obs_norm") %in% names(df))) {
    p_norm <- ggplot(df, aes(x = day)) +
      geom_line(aes(y = pred_norm), color = "#1f77b4", linewidth = 0.8) +
      geom_line(aes(y = obs_norm), color = "black", linetype = "dashed", linewidth = 0.5, na.rm = TRUE) +
      geom_point(aes(y = obs_norm), color = "black", size = 0.9, na.rm = TRUE) +
      facet_wrap(~ harvest, scales = "free_x") +
      coord_cartesian(ylim = c(-1, 1)) +
      labs(
        title = "Selected Step Burden Fit (Normalized)",
        subtitle = basename(step_dir),
        x = "Day",
        y = "Normalized burden"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "selected_burden_fit_normalized.pdf"), p_norm, width = 12, height = 8)
  }

  pred_abs_col <- if ("pred_burden_volume_mm3" %in% names(df)) "pred_burden_volume_mm3" else "pred_pop"
  abs_y_label <- if (identical(pred_abs_col, "pred_burden_volume_mm3")) "Tumor burden (mm^3)" else "Tumor burden / pop"

  p_abs <- ggplot(df, aes(x = day)) +
    geom_line(aes_string(y = pred_abs_col), color = "#1f77b4", linewidth = 0.8) +
    geom_line(aes(y = obs_burden), color = "black", linetype = "dashed", linewidth = 0.5, na.rm = TRUE) +
    geom_point(aes(y = obs_burden), color = "black", size = 0.9, na.rm = TRUE) +
    facet_wrap(~ harvest, scales = "free_y") +
    labs(
      title = "Selected Step Burden Fit (Absolute)",
      subtitle = basename(step_dir),
      x = "Day",
      y = abs_y_label
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "selected_burden_fit_absolute.pdf"), p_abs, width = 12, height = 8)
}

plot_ploidy_fit <- function(step_dir, out_dir, top_n = 8L) {
  f <- file.path(step_dir, "global_terminal_ploidy_fit.tsv")
  if (!file.exists(f)) return(invisible(NULL))
  df <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible(NULL))
  for (nm in intersect(c("N", "pred_fraction", "obs_count", "dose"), names(df))) {
    df[[nm]] <- safe_numeric(df[[nm]])
  }
  df$harvest_f <- factor(df$harvest, levels = sort(unique(df$harvest)))
  df$obs_fraction <- df %>%
    group_by(scenario_id) %>%
    mutate(obs_fraction = ifelse(sum(obs_count, na.rm = TRUE) > 0, obs_count / sum(obs_count, na.rm = TRUE), NA_real_)) %>%
    pull(obs_fraction)

  p_heat <- ggplot(df, aes(x = N, y = harvest_f, fill = pred_fraction)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", na.value = "grey90") +
    labs(
      title = "Selected Step Terminal Ploidy Prediction (Heatmap)",
      subtitle = basename(step_dir),
      x = "Chromosome number (N)",
      y = "Harvest",
      fill = "Pred frac"
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "selected_terminal_ploidy_heatmap.pdf"), p_heat, width = 12, height = 7)

  top_states <- df %>%
    group_by(N) %>%
    summarise(mean_pred = mean(safe_numeric(pred_fraction), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_pred)) %>%
    slice_head(n = top_n) %>%
    pull(N)
  dft <- df %>% filter(N %in% top_states) %>% mutate(N = factor(N, levels = sort(unique(top_states))))

  p_lines <- ggplot(dft, aes(x = N, group = 1)) +
    geom_col(aes(y = pred_fraction), fill = "#4C78A8", alpha = 0.7) +
    geom_point(aes(y = obs_fraction), color = "black", size = 0.8, na.rm = TRUE) +
    geom_line(aes(y = obs_fraction, group = scenario_id), color = "black", linewidth = 0.3, na.rm = TRUE) +
    facet_wrap(~ harvest, scales = "free_y") +
    labs(
      title = paste0("Selected Step Terminal Ploidy (Top ", length(unique(dft$N)), " states)"),
      subtitle = "Bars=predicted fraction, black points/lines=observed fraction",
      x = "N",
      y = "Fraction"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(file.path(out_dir, "selected_terminal_ploidy_top_states.pdf"), p_lines, width = 13, height = 8)
}

plot_per_sample_loss <- function(step_dir, out_dir) {
  f <- file.path(step_dir, "per_sample_loss.tsv")
  if (!file.exists(f)) return(invisible(NULL))
  df <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible(NULL))
  num_cols <- intersect(c("objective_data", "objective_penalty", "objective_total", "objective_burden_scaled", "objective_ploidy_scaled", "data_weight"), names(df))
  for (nm in num_cols) df[[nm]] <- safe_numeric(df[[nm]])

  metric_cols <- intersect(c("objective_data", "objective_burden_scaled", "objective_ploidy_scaled", "objective_penalty", "objective_total"), names(df))
  if (length(metric_cols) == 0) return(invisible(NULL))
  long <- df %>%
    select(any_of(c("scenario_id", "harvest", metric_cols))) %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "value")

  p <- ggplot(long, aes(x = reorder(harvest, value, FUN = median, na.rm = TRUE), y = value, fill = metric)) +
    geom_col(position = "dodge", na.rm = TRUE) +
    coord_flip() +
    labs(
      title = "Selected Step Per-sample Loss Components",
      subtitle = basename(step_dir),
      x = "Harvest",
      y = "Value",
      fill = "Metric"
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "selected_per_sample_loss.pdf"), p, width = 12, height = 8)
}

plot_theta_distributions <- function(step_dir, out_dir) {
  plot_one <- function(f_in, stem) {
    if (!file.exists(f_in)) return(invisible(NULL))
    df <- read.delim(f_in, check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(df) == 0) return(invisible(NULL))
    meta_cols <- intersect(c("scenario_id", "harvest", "cohort", "dose", "n_obs_burden", "n_obs_ploidy", "data_weight"), names(df))
    par_cols <- setdiff(names(df), meta_cols)
    if (length(par_cols) == 0) return(invisible(NULL))

    for (nm in par_cols) df[[nm]] <- safe_numeric(df[[nm]])
    long <- df %>%
      select(all_of(c(meta_cols, par_cols))) %>%
      pivot_longer(cols = all_of(par_cols), names_to = "parameter", values_to = "value")

    finite_long <- long %>% filter(is.finite(value))
    if (nrow(finite_long) == 0) return(invisible(NULL))

    p <- ggplot(finite_long, aes(x = parameter, y = value)) +
      geom_boxplot(outlier.shape = NA, fill = "#9ecae1", alpha = 0.7) +
      geom_jitter(width = 0.15, size = 1.2, alpha = 0.75, color = "#08519c") +
      coord_flip() +
      labs(
        title = paste0("Selected Step Per-sample Parameters (", stem, ")"),
        subtitle = basename(step_dir),
        x = NULL,
        y = "Value"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, paste0("selected_", stem, "_distribution.pdf")), p, width = 10, height = 7)
  }

  plot_one(file.path(step_dir, "per_sample_theta_i.tsv"), "theta_natural")
  plot_one(file.path(step_dir, "per_sample_theta_i_transformed.tsv"), "theta_transformed")
}

plot_theta0 <- function(step_dir, out_dir) {
  f <- file.path(step_dir, "theta0_robust.tsv")
  if (!file.exists(f)) return(invisible(NULL))
  df <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible(NULL))
  if (!all(c("space", "parameter", "value") %in% names(df))) return(invisible(NULL))
  df$value <- safe_numeric(df$value)
  df <- df %>% filter(is.finite(value))
  if (nrow(df) == 0) return(invisible(NULL))

  p <- ggplot(df, aes(x = reorder(parameter, value), y = value, fill = space)) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(
      title = "Selected Step theta0_robust",
      subtitle = basename(step_dir),
      x = NULL,
      y = "Value",
      fill = "Space"
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "selected_theta0_robust.pdf"), p, width = 10, height = 7)
}

write_selected_step_manifest <- function(run_dir, step_dir, out_dir, na_sum) {
  chain <- read_tsv_if_exists(file.path(run_dir, "chain_global_summary.tsv"))
  sel <- read_tsv_if_exists(file.path(run_dir, "selected_best_step.tsv"))

  manifest <- data.frame(
    key = c(
      "run_dir", "selected_step_dir", "selected_step_basename",
      "has_chain_global_summary", "has_selected_best_step", "step_na_files_with_na", "step_na_total_cells"
    ),
    value = c(
      normalizePath(run_dir, mustWork = FALSE),
      normalizePath(step_dir, mustWork = FALSE),
      basename(step_dir),
      as.character(!is.null(chain)),
      as.character(!is.null(sel)),
      as.character(sum(na_sum$na_cells > 0)),
      as.character(sum(na_sum$na_cells))
    ),
    stringsAsFactors = FALSE
  )
  write.table(manifest, file = file.path(out_dir, "viz_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

run_viz_for_run_dir <- function(run_dir, argv, top_n = 8L) {
  out_dir <- file.path(run_dir, "viz_tmb")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  step_dir <- pick_selected_step(run_dir, argv)

  na_sum <- collect_step_na_summary(step_dir)
  if (nrow(na_sum) > 0) {
    write.table(na_sum, file = file.path(out_dir, "selected_step_na_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    p_na <- ggplot(na_sum, aes(x = reorder(file, na_cells), y = na_cells, fill = na_cells > 0)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c(`TRUE` = "#d7301f", `FALSE` = "#74c476"), guide = "none") +
      labs(
        title = "Selected Step NA Cell Counts by File",
        subtitle = paste0(basename(run_dir), " | ", basename(step_dir)),
        x = NULL,
        y = "NA cells"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "selected_step_na_summary.pdf"), p_na, width = 10, height = 6)
  }

  make_chain_plots(run_dir, out_dir)
  plot_burden_fit(step_dir, out_dir)
  plot_ploidy_fit(step_dir, out_dir, top_n = top_n)
  plot_per_sample_loss(step_dir, out_dir)
  plot_theta_distributions(step_dir, out_dir)
  plot_theta0(step_dir, out_dir)
  write_selected_step_manifest(run_dir, step_dir, out_dir, na_sum)

  normalizePath(out_dir)
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  script_dir <- get_script_dir_self()
  results_root <- normalizePath(file.path(script_dir, "..", "results"), mustWork = FALSE)

  fit_root <- if (!is.null(argv$fit_dir)) {
    normalizePath(argv$fit_dir, mustWork = TRUE)
  } else {
    normalizePath(find_latest_fit_dir(results_root), mustWork = TRUE)
  }

  top_n <- as_int(argv$top_n, 8L)
  if (!is.finite(top_n) || top_n < 1L) stop("top_n must be >= 1")

  run_dirs <- find_tmb_run_dirs_under(fit_root)
  if (length(run_dirs) == 0) {
    stop("No TMB hierarchical run directories found under: ", fit_root)
  }

  message("Found ", length(run_dirs), " TMB run directories under: ", fit_root)
  ok <- character(0)
  failed <- character(0)
  for (i in seq_along(run_dirs)) {
    run_dir <- run_dirs[[i]]
    message("[", i, "/", length(run_dirs), "] Processing: ", run_dir)
    tryCatch(
      {
        out_dir <- run_viz_for_run_dir(run_dir, argv = argv, top_n = top_n)
        ok <- c(ok, run_dir)
        message("  Done: ", out_dir)
      },
      error = function(e) {
        failed <<- c(failed, paste0(run_dir, " :: ", conditionMessage(e)))
        message("  Failed: ", conditionMessage(e))
      }
    )
  }

  if (length(ok) == 0) stop("All runs failed.")
  message("Finished. Success: ", length(ok), " | Failed: ", length(failed))
  if (length(failed) > 0) {
    message("Failed runs:")
    for (x in failed) message("  - ", x)
  }
}

if (sys.nframe() == 0) {
  main()
}
