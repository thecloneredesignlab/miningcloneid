#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

.o2sd_live_compare_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2sd_live_compare_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_live_compare_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_bool <- o2sd_as_bool

default_template <- "/Users/4482173/Documents/GitHub/Constant_WGD/oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma}"
default_out_dir <- "/Users/4482173/Documents/GitHub/Constant_WGD/oxygen/results/comp_live_effective_pms"
default_sigma_caps <- c("0.05", "0.15", "0.3", "0.6")

build_task_manifest <- function(sigma_caps, run_dir_template, max_seeds = NULL) {
  task_rows <- list()
  idx <- 0L
  for (sigma_cap in sigma_caps) {
    run_dir <- normalizePath(build_run_dir(run_dir_template, sigma_cap), mustWork = TRUE)
    seed_dirs <- list_seed_dirs(run_dir = run_dir, max_seeds = max_seeds)
    for (seed_dir in seed_dirs) {
      idx <- idx + 1L
      task_rows[[idx]] <- data.frame(
        task_id = idx,
        sigma_cap = sigma_cap,
        run_dir = run_dir,
        seed = basename(seed_dir),
        seed_dir = normalizePath(seed_dir, mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, task_rows)
}

resolve_path_value <- function(path_value, base_dir = getwd()) {
  txt <- path_value
  if (is.null(txt) || !length(txt)) return(NULL)
  txt <- as.character(txt[[1]])
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  if (startsWith(txt, "~")) return(normalizePath(path.expand(txt), mustWork = FALSE))
  if (grepl("^(/|[A-Za-z]:[/\\\\])", txt)) return(normalizePath(txt, mustWork = FALSE))
  normalizePath(file.path(base_dir, txt), mustWork = FALSE)
}

parse_sigma_caps <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(default_sigma_caps)
  vals <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  vals <- vals[nzchar(vals)]
  if (!length(vals)) stop("sigma_caps must contain at least one comma-separated value.")
  vals
}

build_run_dir <- function(template, sigma_cap) {
  if (!grepl("\\{sigma\\}", template)) {
    stop("run_dir_template must contain the placeholder {sigma}.")
  }
  gsub("\\{sigma\\}", sigma_cap, template)
}

list_seed_dirs <- function(run_dir, max_seeds = NULL) {
  seed_dirs <- list.dirs(run_dir, full.names = TRUE, recursive = FALSE)
  seed_dirs <- seed_dirs[grepl("/seed[0-9]+$", seed_dirs) | grepl("\\\\seed[0-9]+$", seed_dirs)]
  if (!length(seed_dirs)) {
    stop("No seed directories were found in ", run_dir)
  }
  seed_ids <- basename(seed_dirs)
  seed_num <- suppressWarnings(as.numeric(sub("^seed", "", seed_ids)))
  ord <- order(seed_num, seed_ids, na.last = TRUE)
  seed_dirs <- seed_dirs[ord]
  if (!is.null(max_seeds) && is.finite(max_seeds) && max_seeds > 0) {
    seed_dirs <- head(seed_dirs, as.integer(max_seeds))
  }
  seed_dirs
}

run_estimate_if_needed <- function(seed_dir, estimate_script, rerun_estimate = FALSE) {
  out_dir <- file.path(seed_dir, "viz", "live_effective_pms")
  overall_path <- file.path(out_dir, "live_effective_pms_overall.tsv")
  if (!rerun_estimate && file.exists(overall_path)) {
    return(out_dir)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  args <- c(
    estimate_script,
    paste0("--seed_dir=", seed_dir),
    paste0("--out_dir=", out_dir)
  )
  status <- system2(rscript_bin, args = args)
  if (!identical(status, 0L)) {
    stop("estimate_live_effective_pms.R failed for ", seed_dir, " with exit status ", status)
  }
  if (!file.exists(overall_path)) {
    stop("estimate_live_effective_pms.R completed but expected output was not found: ", overall_path)
  }
  out_dir
}

read_required_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file was not found: ", path)
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

read_one_seed_summary <- function(seed_dir, sigma_cap, estimate_script, rerun_estimate = FALSE) {
  out_dir <- run_estimate_if_needed(
    seed_dir = seed_dir,
    estimate_script = estimate_script,
    rerun_estimate = rerun_estimate
  )
  overall_path <- file.path(out_dir, "live_effective_pms_overall.tsv")
  harvest_path <- file.path(out_dir, "live_effective_pms_harvest_only.tsv")
  overall_tab <- read_required_tsv(overall_path)
  harvest_tab <- read_required_tsv(harvest_path)
  required_cols <- c("summary_scope", "p_misseg_parameter", "live_weighted_effective_p_ms_mean")
  if (!all(required_cols %in% names(overall_tab))) {
    stop("Missing required columns in ", overall_path, ": ", paste(setdiff(required_cols, names(overall_tab)), collapse = ", "))
  }
  if (!all(required_cols %in% names(harvest_tab))) {
    stop("Missing required columns in ", harvest_path, ": ", paste(setdiff(required_cols, names(harvest_tab)), collapse = ", "))
  }
  overall_row <- overall_tab[overall_tab$summary_scope == "all_sample_days", , drop = FALSE]
  if (!nrow(overall_row)) overall_row <- overall_tab[1, , drop = FALSE]
  harvest_row <- harvest_tab[harvest_tab$summary_scope == "harvest_only", , drop = FALSE]
  if (!nrow(harvest_row)) harvest_row <- harvest_tab[1, , drop = FALSE]

  p_misseg_parameter <- suppressWarnings(as.numeric(overall_row$p_misseg_parameter[[1]]))
  live_cell_p_misseg <- suppressWarnings(as.numeric(overall_row$live_weighted_effective_p_ms_mean[[1]]))
  harvest_live_cell_p_misseg <- suppressWarnings(as.numeric(harvest_row$live_weighted_effective_p_ms_mean[[1]]))

  data.frame(
    sigma_cap = sigma_cap,
    sigma_cap_num = suppressWarnings(as.numeric(sigma_cap)),
    seed = basename(seed_dir),
    seed_dir = normalizePath(seed_dir, mustWork = TRUE),
    live_effective_pms_dir = normalizePath(out_dir, mustWork = TRUE),
    p_misseg_parameter = p_misseg_parameter,
    live_cell_p_misseg = live_cell_p_misseg,
    harvest_live_cell_p_misseg = harvest_live_cell_p_misseg,
    abs_diff_live_minus_parameter = live_cell_p_misseg - p_misseg_parameter,
    ratio_live_over_parameter = live_cell_p_misseg / pmax(p_misseg_parameter, 1e-12),
    stringsAsFactors = FALSE
  )
}

process_one_task <- function(task_row, estimate_script, rerun_estimate = FALSE) {
  task_seed_dir <- as.character(task_row$seed_dir[[1]])
  task_sigma_cap <- as.character(task_row$sigma_cap[[1]])
  task_seed <- as.character(task_row$seed[[1]])
  tryCatch(
    {
      message("Processing sigma_burden=", task_sigma_cap, " seed=", task_seed)
      out <- read_one_seed_summary(
        seed_dir = task_seed_dir,
        sigma_cap = task_sigma_cap,
        estimate_script = estimate_script,
        rerun_estimate = rerun_estimate
      )
      out$task_status <- "ok"
      out$task_error <- ""
      out
    },
    error = function(e) {
      data.frame(
        sigma_cap = task_sigma_cap,
        sigma_cap_num = suppressWarnings(as.numeric(task_sigma_cap)),
        seed = task_seed,
        seed_dir = task_seed_dir,
        live_effective_pms_dir = NA_character_,
        p_misseg_parameter = NA_real_,
        live_cell_p_misseg = NA_real_,
        harvest_live_cell_p_misseg = NA_real_,
        abs_diff_live_minus_parameter = NA_real_,
        ratio_live_over_parameter = NA_real_,
        task_status = "error",
        task_error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
}

save_violin_plot <- function(seed_summary, out_path, sigma_levels) {
  plot_df <- rbind(
    data.frame(
      sigma_cap = seed_summary$sigma_cap,
      estimate_type = "p_misseg parameter",
      value = seed_summary$p_misseg_parameter,
      stringsAsFactors = FALSE
    ),
    data.frame(
      sigma_cap = seed_summary$sigma_cap,
      estimate_type = "live-cell effective p_misseg",
      value = seed_summary$live_cell_p_misseg,
      stringsAsFactors = FALSE
    )
  )
  plot_df <- plot_df[is.finite(plot_df$value), , drop = FALSE]
  plot_df$sigma_cap <- factor(plot_df$sigma_cap, levels = sigma_levels)
  plot_df$estimate_type <- factor(
    plot_df$estimate_type,
    levels = c("p_misseg parameter", "live-cell effective p_misseg")
  )
  dodge <- position_dodge(width = 0.8)
  p <- ggplot(plot_df, aes(x = sigma_cap, y = value, fill = estimate_type)) +
    geom_violin(position = dodge, trim = FALSE, alpha = 0.55, color = NA) +
    geom_boxplot(
      position = dodge,
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.9,
      linewidth = 0.35
    ) +
    scale_fill_manual(
      values = c(
        "p_misseg parameter" = "#4c78a8",
        "live-cell effective p_misseg" = "#f58518"
      ),
      drop = FALSE
    ) +
    labs(
      title = "p_misseg vs Live-Cell Effective p_misseg by sigma_burden Upper Bound",
      subtitle = "Each sigma_burden group includes all seeds. The live-cell value is the all-sample-days live-weighted effective p_ms mean from estimate_live_effective_pms.R.",
      x = "sigma_burden upper bound",
      y = "p_misseg estimate",
      fill = "Estimate type"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 8, height = 8)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  sigma_caps <- parse_sigma_caps(argv$sigma_caps)
  run_dir_template <- argv$run_dir_template %||% default_template
  out_dir <- resolve_path_value(argv$out_dir, getwd()) %||% default_out_dir
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  estimate_script <- resolve_path_value(argv$estimate_script, SCRIPT_DIR) %||%
    normalizePath(file.path(SCRIPT_DIR, "estimate_live_effective_pms.R"), mustWork = TRUE)
  max_seeds <- NULL
  if (!is.null(argv$max_seeds) && nzchar(trimws(as.character(argv$max_seeds)))) {
    max_seeds <- suppressWarnings(as.integer(argv$max_seeds))
    if (!is.finite(max_seeds) || max_seeds <= 0) {
      stop("--max_seeds must be a positive integer when supplied.")
    }
  }
  rerun_estimate <- as_bool(argv$rerun_estimate, FALSE)
  n_cores <- suppressWarnings(as.integer(argv$n_cores %||% 1L))
  if (!is.finite(n_cores) || n_cores <= 0L) {
    stop("--n_cores must be a positive integer.")
  }

  task_manifest <- build_task_manifest(
    sigma_caps = sigma_caps,
    run_dir_template = run_dir_template,
    max_seeds = max_seeds
  )
  utils::write.table(
    task_manifest,
    file = file.path(out_dir, "sigma_burden_live_effective_pms_task_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (n_cores == 1L || nrow(task_manifest) <= 1L) {
    seed_rows <- lapply(
      seq_len(nrow(task_manifest)),
      function(i) process_one_task(
        task_row = task_manifest[i, , drop = FALSE],
        estimate_script = estimate_script,
        rerun_estimate = rerun_estimate
      )
    )
  } else if (.Platform$OS.type == "unix") {
    worker_count <- min(n_cores, nrow(task_manifest))
    seed_rows <- parallel::mclapply(
      X = seq_len(nrow(task_manifest)),
      FUN = function(i) process_one_task(
        task_row = task_manifest[i, , drop = FALSE],
        estimate_script = estimate_script,
        rerun_estimate = rerun_estimate
      ),
      mc.cores = worker_count,
      mc.preschedule = FALSE
    )
  } else {
    worker_count <- min(n_cores, nrow(task_manifest))
    cl <- parallel::makeCluster(worker_count)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c(
        "read_required_tsv",
        "run_estimate_if_needed",
        "read_one_seed_summary",
        "process_one_task"
      ),
      envir = environment()
    )
    seed_rows <- parallel::parLapplyLB(
      cl,
      X = seq_len(nrow(task_manifest)),
      fun = function(i, task_tab, estimate_script_arg, rerun_estimate_arg) {
        process_one_task(
          task_row = task_tab[i, , drop = FALSE],
          estimate_script = estimate_script_arg,
          rerun_estimate = rerun_estimate_arg
        )
      },
      task_tab = task_manifest,
      estimate_script_arg = estimate_script,
      rerun_estimate_arg = rerun_estimate
    )
  }

  seed_summary <- do.call(rbind, seed_rows)
  utils::write.table(
    seed_summary,
    file = file.path(out_dir, "sigma_burden_live_effective_pms_task_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  failed_rows <- seed_summary[seed_summary$task_status != "ok", , drop = FALSE]
  if (nrow(failed_rows) > 0L) {
    stop(
      "One or more tasks failed. See ",
      file.path(out_dir, "sigma_burden_live_effective_pms_task_results.tsv"),
      ". First error: sigma_burden=", failed_rows$sigma_cap[[1]],
      " seed=", failed_rows$seed[[1]],
      " message=", failed_rows$task_error[[1]]
    )
  }
  seed_summary$task_status <- NULL
  seed_summary$task_error <- NULL
  sigma_levels <- sigma_caps
  seed_summary$sigma_cap <- factor(seed_summary$sigma_cap, levels = sigma_levels)
  seed_summary <- seed_summary[order(seed_summary$sigma_cap, seed_summary$seed), , drop = FALSE]

  summary_by_sigma <- do.call(
    rbind,
    lapply(
      split(seed_summary, seed_summary$sigma_cap),
      function(df) data.frame(
        sigma_cap = unique(as.character(df$sigma_cap)),
        n_seeds = nrow(df),
        p_misseg_parameter_mean = mean(df$p_misseg_parameter, na.rm = TRUE),
        p_misseg_parameter_median = stats::median(df$p_misseg_parameter, na.rm = TRUE),
        live_cell_p_misseg_mean = mean(df$live_cell_p_misseg, na.rm = TRUE),
        live_cell_p_misseg_median = stats::median(df$live_cell_p_misseg, na.rm = TRUE),
        harvest_live_cell_p_misseg_mean = mean(df$harvest_live_cell_p_misseg, na.rm = TRUE),
        abs_diff_live_minus_parameter_mean = mean(df$abs_diff_live_minus_parameter, na.rm = TRUE),
        ratio_live_over_parameter_mean = mean(df$ratio_live_over_parameter, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  )

  utils::write.table(
    seed_summary,
    file = file.path(out_dir, "sigma_burden_live_effective_pms_by_seed.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    summary_by_sigma,
    file = file.path(out_dir, "sigma_burden_live_effective_pms_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  save_violin_plot(
    seed_summary = seed_summary,
    out_path = file.path(out_dir, "sigma_burden_p_misseg_vs_live_cell_violin.pdf"),
    sigma_levels = sigma_levels
  )

  message("Wrote by-seed table: ", file.path(out_dir, "sigma_burden_live_effective_pms_by_seed.tsv"))
  message("Wrote summary table: ", file.path(out_dir, "sigma_burden_live_effective_pms_summary.tsv"))
  message("Wrote figure: ", file.path(out_dir, "sigma_burden_p_misseg_vs_live_cell_violin.pdf"))
  message("Wrote task manifest: ", file.path(out_dir, "sigma_burden_live_effective_pms_task_manifest.tsv"))
  message("Wrote task results: ", file.path(out_dir, "sigma_burden_live_effective_pms_task_results.tsv"))
}

if (sys.nframe() == 0) {
  main()
}
