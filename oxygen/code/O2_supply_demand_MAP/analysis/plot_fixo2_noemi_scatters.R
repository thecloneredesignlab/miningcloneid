#!/usr/bin/env Rscript

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    frame_files <- Filter(
      nzchar,
      vapply(sys.frames(), function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      }, character(1))
    )
    if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
  }
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    if (!grepl("=", arg, fixed = TRUE)) {
      out[[arg]] <- TRUE
    } else {
      key <- sub("=.*$", "", arg)
      val <- sub("^[^=]*=", "", arg)
      out[[key]] <- val
    }
  }
  out
}

repo_root <- function() {
  normalizePath(file.path(SCRIPT_DIR, "..", "..", "..", ".."), mustWork = FALSE)
}

resolve_repo_path <- function(path, root = repo_root(), mustWork = FALSE) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path)) return(path)
  path <- as.character(path[[1]])
  if (startsWith(path, "~/")) {
    path <- file.path(root, substring(path, 3L))
  } else if (identical(path, "~")) {
    path <- root
  } else if (!grepl("^/", path)) {
    path <- file.path(root, path)
  }
  normalizePath(path, mustWork = mustWork)
}

split_csv <- function(x, default = character()) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(trimws(as.character(x[[1]])))) {
    return(default)
  }
  vals <- trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

as_num_vec <- function(x, default) {
  vals <- suppressWarnings(as.numeric(split_csv(x, as.character(default))))
  vals <- vals[is.finite(vals)]
  if (length(vals)) unique(vals) else default
}

as_int_vec <- function(x, default) {
  vals <- suppressWarnings(as.integer(split_csv(x, as.character(default))))
  vals <- vals[!is.na(vals)]
  if (length(vals)) unique(vals) else default
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  tolower(as.character(x[[1]])) %in% c("true", "1", "yes", "y", "on")
}

num_key <- function(x) {
  format(signif(as.numeric(x), 12), scientific = FALSE, trim = TRUE)
}

num_path_tag <- function(x) {
  val <- suppressWarnings(as.numeric(x))
  if (length(val) != 1L || !is.finite(val)) return("NA")
  txt <- format(val, scientific = FALSE, trim = TRUE)
  txt <- sub("^-", "m", txt)
  txt <- gsub("\\.", "p", txt)
  txt <- gsub("[^A-Za-z0-9]+", "", txt)
  if (!nzchar(txt)) "NA" else txt
}

read_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing file: ", path)
    return(data.frame())
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE)))
  }
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

format_o2_label <- function(x) {
  paste0("O2 = ", format(as.numeric(x), scientific = FALSE, trim = TRUE))
}

initial_condition_from_ploidy <- function(x) {
  paste0("init_", format(as.numeric(x), scientific = FALSE, trim = TRUE), "N")
}

filter_by_values <- function(df, col, values) {
  if (!nrow(df) || !(col %in% names(df))) return(df)
  key <- num_key(df[[col]])
  keep <- key %in% num_key(values)
  df[keep, , drop = FALSE]
}

read_analytical_trajectories <- function(analysis_dir, time_points, o2_values, analytical_table = NULL) {
  if (is.null(analytical_table) || !nzchar(analytical_table)) {
    analytical_table <- file.path(
      analysis_dir,
      "counterfactual_trajectories",
      "tables",
      "fixed_o2_counterfactual_trajectories.tsv"
    )
  }
  tab <- read_tsv(analytical_table)
  required <- c("seed_id", "O2_pct", "initial_condition", "day", "mean_ploidy")
  missing <- setdiff(required, names(tab))
  if (length(missing)) stop("Analytical trajectory table is missing column(s): ", paste(missing, collapse = ", "))
  if ("status" %in% names(tab)) tab <- tab[tab$status %in% "ok", , drop = FALSE]
  tab <- filter_by_values(tab, "day", time_points)
  tab <- filter_by_values(tab, "O2_pct", o2_values)
  tab$day_key <- num_key(tab$day)
  tab$O2_key <- num_key(tab$O2_pct)
  tab$analytical_mean_ploidy <- as.numeric(tab$mean_ploidy)
  keep <- intersect(
    c("seed_id", "trajectory_regime", "mode_label", "O2_pct", "O2_key", "initial_condition", "day", "day_key", "analytical_mean_ploidy"),
    names(tab)
  )
  out <- tab[, keep, drop = FALSE]
  out[is.finite(out$analytical_mean_ploidy), , drop = FALSE]
}

read_seed_objectives <- function(analysis_dir, fit_dir = NULL) {
  attractor_path <- file.path(
    analysis_dir,
    "attractors",
    "tables",
    "fixed_o2_attractor_spectral_gap_by_seed.tsv"
  )
  if (file.exists(attractor_path)) {
    tab <- read_tsv(attractor_path)
    cols <- intersect(c("seed_id", "trajectory_regime", "mode_label", "objective", "delta_objective"), names(tab))
    if (all(c("seed_id", "objective") %in% cols)) {
      tab <- tab[, cols, drop = FALSE]
      tab <- tab[!duplicated(tab$seed_id), , drop = FALSE]
      tab$objective <- as.numeric(tab$objective)
      return(tab)
    }
  }

  if (is.null(fit_dir) || !nzchar(fit_dir) || !dir.exists(fit_dir)) {
    warning("Seed objective values were not found in the analysis directory, and fit_dir is unavailable.")
    return(data.frame(seed_id = character(), objective = numeric(), stringsAsFactors = FALSE))
  }

  seed_dirs <- list.dirs(fit_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("seed[0-9]+$", basename(seed_dirs))]
  rows <- lapply(seed_dirs, function(seed_dir) {
    candidates <- c(
      file.path(seed_dir, "best_params.tsv"),
      file.path(seed_dir, "best_parameters.tsv"),
      file.path(seed_dir, "parameter_table_input.csv")
    )
    hits <- candidates[file.exists(candidates)]
    path <- if (length(hits)) hits[[1]] else NA_character_
    if (is.na(path)) return(NULL)
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    tab <- tryCatch(
      utils::read.table(path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) data.frame()
    )
    objective_cols <- intersect(
      c("objective", "optimizer_local_objective", "optimizer_deoptim_objective", "loss", "best_objective"),
      names(tab)
    )
    if (!length(objective_cols) || !nrow(tab)) return(NULL)
    data.frame(
      seed_id = basename(seed_dir),
      objective = suppressWarnings(as.numeric(tab[[objective_cols[[1]]]][[1]])),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) out <- data.frame(seed_id = character(), objective = numeric(), stringsAsFactors = FALSE)
  out
}

task_table_path <- function(simulation_dir, simulation_mode) {
  file.path(simulation_dir, simulation_mode, "task_list.tsv")
}

build_task_table_from_paths <- function(simulation_dir, simulation_mode, seed_ids, o2_values, initial_ploidy_values, simulation_ids) {
  rows <- list()
  k <- 0L
  for (seed_id in seed_ids) {
    for (O2 in o2_values) {
      for (initial_ploidy in initial_ploidy_values) {
        for (simulation_id in simulation_ids) {
          k <- k + 1L
          output_dir <- file.path(
            simulation_dir,
            simulation_mode,
            paste0("O2_", num_path_tag(O2)),
            paste0("ploidy", num_path_tag(initial_ploidy)),
            seed_id,
            paste0("sim", simulation_id)
          )
          rows[[k]] <- data.frame(
            task_id = k,
            seed_id = seed_id,
            fixed_o2_pct = O2,
            initial_ploidy = initial_ploidy,
            initial_condition = initial_condition_from_ploidy(initial_ploidy),
            simulation_id = simulation_id,
            output_dir = output_dir,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

read_simulation_tasks <- function(simulation_dir, simulation_mode, analytical, o2_values, initial_ploidy_values, simulation_ids) {
  path <- task_table_path(simulation_dir, simulation_mode)
  if (file.exists(path)) {
    tasks <- read_tsv(path)
    if (!"fixed_o2_pct" %in% names(tasks) && "O2_pct" %in% names(tasks)) {
      names(tasks)[names(tasks) == "O2_pct"] <- "fixed_o2_pct"
    }
    required <- c("seed_id", "fixed_o2_pct", "initial_ploidy", "simulation_id", "output_dir")
    missing <- setdiff(required, names(tasks))
    if (length(missing)) stop("Simulation task table is missing column(s): ", paste(missing, collapse = ", "))
    if (!"initial_condition" %in% names(tasks)) {
      tasks$initial_condition <- initial_condition_from_ploidy(tasks$initial_ploidy)
    }
  } else {
    seed_ids <- sort(unique(analytical$seed_id))
    tasks <- build_task_table_from_paths(
      simulation_dir = simulation_dir,
      simulation_mode = simulation_mode,
      seed_ids = seed_ids,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      simulation_ids = simulation_ids
    )
  }
  tasks$fixed_o2_pct <- as.numeric(tasks$fixed_o2_pct)
  tasks$initial_ploidy <- as.numeric(tasks$initial_ploidy)
  tasks$simulation_id <- as.integer(tasks$simulation_id)
  expected_output_dir <- file.path(
    simulation_dir,
    simulation_mode,
    paste0("O2_", vapply(tasks$fixed_o2_pct, num_path_tag, character(1))),
    paste0("ploidy", vapply(tasks$initial_ploidy, num_path_tag, character(1))),
    tasks$seed_id,
    paste0("sim", tasks$simulation_id)
  )
  use_expected <- !dir.exists(tasks$output_dir) & dir.exists(expected_output_dir)
  tasks$output_dir[use_expected] <- expected_output_dir[use_expected]
  tasks <- filter_by_values(tasks, "fixed_o2_pct", o2_values)
  tasks <- filter_by_values(tasks, "initial_ploidy", initial_ploidy_values)
  tasks <- tasks[tasks$simulation_id %in% simulation_ids, , drop = FALSE]
  tasks <- tasks[tasks$seed_id %in% unique(analytical$seed_id), , drop = FALSE]
  tasks$state_file <- file.path(tasks$output_dir, "state_trajectory.tsv.gz")
  tasks
}

read_state_metrics_awk <- function(path, time_points) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    return(data.frame())
  }
  gzip <- Sys.which("gzip")
  awk <- Sys.which("awk")
  if (!nzchar(gzip) || !nzchar(awk)) {
    stop("Reading state trajectories requires gzip and awk on PATH.")
  }
  days_arg <- paste(num_key(time_points), collapse = ",")
  awk_script <- paste(
    'BEGIN {',
    '  split(days, d, ",");',
    '  for (i in d) keep[sprintf("%.10g", d[i] + 0)] = 1;',
    '  OFS = "\t";',
    '}',
    'NR == 1 {',
    '  for (i = 1; i <= NF; i++) idx[$i] = i;',
    '  next;',
    '}',
    '{',
    '  day = sprintf("%.10g", $(idx["day"]) + 0);',
    '  if ((day in keep) && $(idx["status"]) == "live") {',
    '    c = $(idx["cell_count"]) + 0;',
    '    n = $(idx["N"]) + 0;',
    '    p = $(idx["ploidy"]) + 0;',
    '    sumc[day] += c;',
    '    sumn[day] += n * c;',
    '    sump[day] += p * c;',
    '    sump2[day] += p * p * c;',
    '  }',
    '}',
    'END {',
    '  for (day in keep) {',
    '    if (sumc[day] > 0) {',
    '      meanp = sump[day] / sumc[day];',
    '      varp = sump2[day] / sumc[day] - meanp * meanp;',
    '      if (varp < 0 && varp > -1e-9) varp = 0;',
    '      print day, sumn[day] / sumc[day], meanp, sqrt(varp), sumc[day];',
    '    }',
    '  }',
    '}',
    sep = "\n"
  )
  cmd <- paste(
    shQuote(gzip), "-cd", shQuote(path), "|",
    shQuote(awk), "-v", shQuote(paste0("days=", days_arg)), shQuote(awk_script)
  )
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())
  if (!length(out)) return(data.frame())
  con <- textConnection(out)
  on.exit(close(con), add = TRUE)
  tab <- utils::read.table(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(tab) <- c("day", "simulation_mean_N", "simulation_mean_ploidy", "simulation_sd_ploidy", "simulation_live_cells")
  tab$day <- as.numeric(tab$day)
  tab
}

read_simulation_metrics <- function(tasks, time_points, progress_every = 100L) {
  rows <- vector("list", nrow(tasks))
  missing <- 0L
  for (i in seq_len(nrow(tasks))) {
    if (progress_every > 0L && (i == 1L || i %% progress_every == 0L || i == nrow(tasks))) {
      message("Reading simulation state metrics: ", i, "/", nrow(tasks))
    }
    task <- tasks[i, , drop = FALSE]
    metric <- read_state_metrics_awk(task$state_file[[1]], time_points)
    if (!nrow(metric)) {
      missing <- missing + 1L
      next
    }
    metric$seed_id <- task$seed_id[[1]]
    metric$O2_pct <- task$fixed_o2_pct[[1]]
    metric$O2_key <- num_key(metric$O2_pct)
    metric$initial_ploidy <- task$initial_ploidy[[1]]
    metric$initial_condition <- task$initial_condition[[1]]
    metric$simulation_id <- task$simulation_id[[1]]
    metric$day_key <- num_key(metric$day)
    rows[[i]] <- metric
  }
  if (missing) warning("Missing or unreadable simulation state files: ", missing)
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) out <- data.frame()
  out
}

aggregate_replicates <- function(sim_metrics) {
  if (!nrow(sim_metrics)) return(sim_metrics)
  keys <- c("seed_id", "O2_pct", "O2_key", "initial_condition", "initial_ploidy", "day", "day_key")
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(sim_metrics)
    out <- dt[, .(
      simulation_n = data.table::uniqueN(simulation_id),
      simulation_mean_N = mean(simulation_mean_N, na.rm = TRUE),
      simulation_sd_replicate_mean_N = stats::sd(simulation_mean_N, na.rm = TRUE),
      simulation_mean_ploidy = mean(simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_replicate_mean_ploidy = stats::sd(simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(simulation_sd_ploidy, na.rm = TRUE),
      simulation_live_cells = mean(simulation_live_cells, na.rm = TRUE)
    ), by = keys]
    return(as.data.frame(out))
  }
  split_key <- interaction(sim_metrics[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(sim_metrics, split_key), function(x) {
    data.frame(
      x[1, keys, drop = FALSE],
      simulation_n = length(unique(x$simulation_id)),
      simulation_mean_N = mean(x$simulation_mean_N, na.rm = TRUE),
      simulation_sd_replicate_mean_N = stats::sd(x$simulation_mean_N, na.rm = TRUE),
      simulation_mean_ploidy = mean(x$simulation_mean_ploidy, na.rm = TRUE),
      simulation_sd_replicate_mean_ploidy = stats::sd(x$simulation_mean_ploidy, na.rm = TRUE),
      simulation_mean_sd_ploidy = mean(x$simulation_sd_ploidy, na.rm = TRUE),
      simulation_live_cells = mean(x$simulation_live_cells, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

merge_scatter_data <- function(analytical, sim_summary, objectives) {
  by_cols <- c("seed_id", "O2_key", "initial_condition", "day_key")
  dat <- merge(
    analytical,
    sim_summary,
    by = by_cols,
    all = FALSE,
    suffixes = c("_analytical", "_simulation")
  )
  if ("O2_pct_analytical" %in% names(dat)) dat$O2_pct <- dat$O2_pct_analytical
  if ("day_analytical" %in% names(dat)) dat$day <- dat$day_analytical
  dat <- merge(dat, objectives, by = "seed_id", all.x = TRUE, suffixes = c("", "_objective"))
  if ("mode_label_objective" %in% names(dat) && "mode_label" %in% names(dat)) {
    fill <- !nzchar(as.character(dat$mode_label)) | is.na(dat$mode_label)
    dat$mode_label[fill] <- dat$mode_label_objective[fill]
  }
  dat$O2_factor <- factor(format_o2_label(dat$O2_pct), levels = format_o2_label(sort(unique(dat$O2_pct))))
  dat$day_factor <- factor(paste0("Day ", format(as.numeric(dat$day), scientific = FALSE, trim = TRUE)),
                           levels = paste0("Day ", format(sort(unique(as.numeric(dat$day))), scientific = FALSE, trim = TRUE)))
  dat$initial_condition <- factor(dat$initial_condition, levels = sort(unique(dat$initial_condition)))
  dat$objective <- as.numeric(dat$objective)
  dat[is.finite(dat$analytical_mean_ploidy) & is.finite(dat$simulation_mean_ploidy), , drop = FALSE]
}

plot_limits <- function(dat) {
  vals <- c(dat$analytical_mean_ploidy, dat$simulation_mean_ploidy)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(c(0, 1))
  rng <- range(vals)
  pad <- diff(rng) * 0.04
  if (!is.finite(pad) || pad <= 0) pad <- 0.05
  rng + c(-pad, pad)
}

objective_aesthetic <- function(dat, transform = "identity") {
  objective <- as.numeric(dat$objective)
  label <- "Objective"
  if (identical(transform, "log10")) {
    objective <- ifelse(is.finite(objective) & objective > 0, log10(objective), NA_real_)
    label <- "log10(objective)"
  }
  list(value = objective, label = label)
}

base_scatter <- function(dat, limits, point_size = 0.9, alpha = 0.55) {
  ggplot2::ggplot(dat, ggplot2::aes(x = analytical_mean_ploidy, y = simulation_mean_ploidy)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey45", linetype = 2, linewidth = 0.35) +
    ggplot2::coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
    ggplot2::labs(
      x = "Analytical solution mean ploidy",
      y = "Simulation-inferred mean ploidy"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65"),
      legend.position = "right"
    )
}

save_plot <- function(plot, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", device = "pdf")
  invisible(path)
}

plot_time_facets_color_o2 <- function(dat, path, limits) {
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = O2_factor, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0.15, color = "grey25"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::labs(fill = "Fixed O2", shape = "Initial condition")
  save_plot(p, path, width = 13, height = 7)
}

plot_time_facets_color_objective <- function(dat, path, limits, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 1.0, alpha = 0.55, stroke = 0.15, color = "grey25"
    ) +
    ggplot2::facet_wrap(~day_factor, nrow = 2, ncol = 4) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition")
  save_plot(p, path, width = 13, height = 7)
}

plot_o2_facets_color_objective <- function(dat, path, limits, title = NULL, objective_transform = "identity") {
  obj <- objective_aesthetic(dat, objective_transform)
  dat$objective_color_value <- obj$value
  p <- base_scatter(dat, limits) +
    ggplot2::geom_point(
      ggplot2::aes(fill = objective_color_value, shape = initial_condition),
      size = 0.75, alpha = 0.42, stroke = 0.12, color = "grey25"
    ) +
    ggplot2::facet_wrap(~O2_factor, nrow = 2, ncol = 3) +
    ggplot2::scale_shape_manual(values = c(21, 24, 22, 23)) +
    ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(9, "viridis"), na.value = "grey80") +
    ggplot2::labs(fill = obj$label, shape = "Initial condition", title = title)
  save_plot(p, path, width = 11, height = 7)
}

make_scatter_outputs <- function(dat, out_dir, objective_transform = "identity") {
  fig_dir <- file.path(out_dir, "simulation", "scatters")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  limits <- plot_limits(dat)

  plot_time_facets_color_o2(
    dat,
    file.path(fig_dir, "scatter_analytical_vs_simulation_by_time_color_o2.pdf"),
    limits = limits
  )
  plot_time_facets_color_objective(
    dat,
    file.path(fig_dir, "scatter_analytical_vs_simulation_by_time_color_objective.pdf"),
    limits = limits,
    objective_transform = objective_transform
  )
  plot_o2_facets_color_objective(
    dat,
    file.path(fig_dir, "scatter_analytical_vs_simulation_by_o2_color_objective_all_times.pdf"),
    limits = limits,
    title = "All selected time points",
    objective_transform = objective_transform
  )

  for (day in sort(unique(as.numeric(dat$day)))) {
    day_dat <- dat[abs(as.numeric(dat$day) - day) < 1e-9, , drop = FALSE]
    day_label <- format(day, scientific = FALSE, trim = TRUE)
    plot_o2_facets_color_objective(
      day_dat,
      file.path(fig_dir, paste0("scatter_analytical_vs_simulation_by_o2_color_objective_day", day_label, ".pdf")),
      limits = limits,
      title = paste0("Day ", day_label),
      objective_transform = objective_transform
    )
  }
  invisible(fig_dir)
}

main <- function(argv = parse_args()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package is required for plotting.")

  root <- resolve_repo_path(argv$repo_root %||% "~", root = repo_root(), mustWork = TRUE)
  simulation_dir <- resolve_repo_path(argv$simulation_dir %||% "~/oxygen/results/O2_fixed_simulation", root = root, mustWork = TRUE)
  analysis_dir <- resolve_repo_path(argv$analysis_dir %||% "~/oxygen/results/analysis/FixO2_invivo_500seed", root = root, mustWork = TRUE)
  fit_dir <- resolve_repo_path(argv$fit_dir %||% "~/oxygen/results/fit_invitro_O2_buffering_500seed", root = root, mustWork = FALSE)
  out_dir <- resolve_repo_path(argv$out_dir %||% "~/oxygen/results/analysis/FixO2_invivo_500seed", root = root, mustWork = FALSE)
  simulation_mode <- argv$simulation_mode %||% "invivo"
  time_points <- sort(as_num_vec(argv$time_points, c(25, 50, 100, 200, 300, 500, 700, 1000)))
  o2_values <- sort(as_num_vec(argv$o2_values, c(0, 0.1, 0.5, 1, 2, 5)))
  initial_ploidy_values <- sort(as_num_vec(argv$initial_ploidy_values, c(2, 4)))
  simulation_ids <- sort(as_int_vec(argv$simulation_ids, c(1L, 2L, 3L)))
  objective_transform <- argv$objective_transform %||% "identity"
  if (!objective_transform %in% c("identity", "log10")) stop("--objective_transform must be identity or log10.")
  recompute <- as_bool(argv$recompute, FALSE)

  table_dir <- file.path(out_dir, "simulation", "scatters", "tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  sim_metric_path <- argv$simulation_metric_table %||% file.path(table_dir, "scatter_simulation_selected_time_metrics_by_replicate.tsv")
  sim_summary_path <- argv$simulation_summary_table %||% file.path(table_dir, "scatter_simulation_selected_time_metrics.tsv")
  scatter_data_path <- argv$scatter_data_table %||% file.path(table_dir, "scatter_analytical_vs_simulation_data.tsv")

  message("Reading analytical trajectories.")
  analytical <- read_analytical_trajectories(
    analysis_dir = analysis_dir,
    time_points = time_points,
    o2_values = o2_values,
    analytical_table = argv$analytical_table %||% NULL
  )
  if (!nrow(analytical)) stop("No analytical trajectory rows were found for the requested O2/time grid.")

  message("Reading seed objective values.")
  objectives <- read_seed_objectives(analysis_dir = analysis_dir, fit_dir = fit_dir)

  if (!recompute && file.exists(sim_summary_path)) {
    message("Reading cached simulation summary: ", sim_summary_path)
    sim_summary <- read_tsv(sim_summary_path)
  } else {
    message("Reading simulation task table.")
    tasks <- read_simulation_tasks(
      simulation_dir = simulation_dir,
      simulation_mode = simulation_mode,
      analytical = analytical,
      o2_values = o2_values,
      initial_ploidy_values = initial_ploidy_values,
      simulation_ids = simulation_ids
    )
    if (!nrow(tasks)) stop("No simulation tasks matched the requested seed/O2/initial ploidy/simulation id filters.")
    message("Matched simulation tasks: ", nrow(tasks))

    sim_metrics <- read_simulation_metrics(
      tasks = tasks,
      time_points = time_points,
      progress_every = as.integer(argv$progress_every %||% 100L)
    )
    if (!nrow(sim_metrics)) stop("No simulation metrics were read from state trajectories.")
    write_tsv(sim_metrics, sim_metric_path)

    sim_summary <- aggregate_replicates(sim_metrics)
    write_tsv(sim_summary, sim_summary_path)
  }

  message("Merging analytical and simulation summaries.")
  scatter_data <- merge_scatter_data(analytical, sim_summary, objectives)
  if (!nrow(scatter_data)) stop("No merged analytical-vs-simulation rows were produced.")
  write_tsv(scatter_data, scatter_data_path)

  message("Drawing scatter plots.")
  fig_dir <- make_scatter_outputs(
    dat = scatter_data,
    out_dir = out_dir,
    objective_transform = objective_transform
  )

  manifest <- data.frame(
    field = c(
      "simulation_dir", "analysis_dir", "fit_dir", "out_dir", "simulation_mode",
      "time_points", "o2_values", "initial_ploidy_values", "simulation_ids",
      "objective_transform", "simulation_metric_table", "simulation_summary_table",
      "scatter_data_table", "figure_dir"
    ),
    value = c(
      simulation_dir, analysis_dir, fit_dir, out_dir, simulation_mode,
      paste(time_points, collapse = ","), paste(o2_values, collapse = ","),
      paste(initial_ploidy_values, collapse = ","), paste(simulation_ids, collapse = ","),
      objective_transform, sim_metric_path, sim_summary_path, scatter_data_path, fig_dir
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(table_dir, "scatter_run_manifest.tsv"))
  message("Done. Scatter figures written to: ", fig_dir)
  invisible(fig_dir)
}

if (identical(environment(), globalenv())) {
  main()
}
