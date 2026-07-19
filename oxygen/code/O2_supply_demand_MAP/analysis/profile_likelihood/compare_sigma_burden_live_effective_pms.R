#!/usr/bin/env Rscript

# Analysis-only comparison of p_misseg and live-cell effective p_ms across
# sigma_burden fits. This script consumes materialized simulation tables and
# never launches a simulation or creates a figure.
#
# Upstream contract for every seed:
#   simulation/invivo/cin/live_effective_pms/
#     live_effective_pms_manifest.tsv
#     live_effective_pms_overall.tsv
#     live_effective_pms_harvest_only.tsv
#
# Usage:
#   Rscript oxygen/code/O2_supply_demand_MAP/analysis/profile_likelihood/compare_sigma_burden_live_effective_pms.R \
#     --run_dir_template=/path/to/run_sigma_{sigma} \
#     --sigma_caps=0.05,0.15,0.3,0.6 \
#     --out_dir=/path/to/analysis/profile_likelihood/live_effective_pms

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
    self <- frame_files[
      basename(frame_files) == "compare_sigma_burden_live_effective_pms.R"
    ]
    if (length(self) > 0L) return(dirname(self[[1L]]))
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2sd_live_compare_bootstrap_script_dir, mustWork = FALSE)
locate_workflow_root <- function(starts) {
  for (start in unique(starts)) {
    current <- normalizePath(start, mustWork = FALSE)
    for (depth in 0:10) {
      candidates <- c(
        current,
        file.path(current, "oxygen", "code", "O2_supply_demand_MAP"),
        file.path(current, "code", "O2_supply_demand_MAP")
      )
      hits <- candidates[file.exists(file.path(
        candidates,
        "util",
        "o2_supply_demand_map_shared.R"
      ))]
      if (length(hits)) return(normalizePath(hits[[1L]], mustWork = TRUE))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate the O2_supply_demand_MAP workflow root.")
}
WORKFLOW_ROOT <- locate_workflow_root(c(SCRIPT_DIR, getwd()))
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_live_compare_bootstrap_script_dir, locate_workflow_root)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
resolve_path_value <- o2sd_resolve_path
read_required_tsv <- o2sd_read_tsv
write_tsv <- o2sd_write_tsv

REPO_ROOT <- normalizePath(
  file.path(WORKFLOW_ROOT, "..", "..", ".."),
  mustWork = TRUE
)
default_template <- file.path(
  REPO_ROOT,
  "oxygen",
  "results",
  "fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma}"
)
default_out_dir <- file.path(
  REPO_ROOT,
  "oxygen",
  "results",
  "comp_live_effective_pms"
)
default_sigma_caps <- c("0.05", "0.15", "0.3", "0.6")
default_live_effective_subdir <- file.path(
  "simulation",
  "invivo",
  "cin",
  "live_effective_pms"
)

parse_sigma_caps <- function(x) {
  if (is.null(x) || !nzchar(trimws(as.character(x)))) return(default_sigma_caps)
  vals <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  vals <- vals[nzchar(vals)]
  if (!length(vals)) stop("sigma_caps must contain at least one comma-separated value.")
  vals
}

build_run_dir <- o2sd_build_sigma_run_dir

list_seed_dirs <- function(run_dir, max_seeds = NULL) {
  seed_dirs <- list.dirs(run_dir, full.names = TRUE, recursive = FALSE)
  seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
  if (!length(seed_dirs)) stop("No seed directories were found in ", run_dir)
  seed_ids <- basename(seed_dirs)
  seed_num <- suppressWarnings(as.numeric(sub("^seed", "", seed_ids)))
  seed_dirs <- seed_dirs[order(seed_num, seed_ids, na.last = TRUE)]
  if (!is.null(max_seeds)) seed_dirs <- head(seed_dirs, max_seeds)
  seed_dirs
}

build_task_manifest <- function(
  sigma_caps,
  run_dir_template,
  live_effective_subdir = default_live_effective_subdir,
  max_seeds = NULL
) {
  task_rows <- list()
  idx <- 0L
  for (sigma_cap in sigma_caps) {
    run_dir <- normalizePath(build_run_dir(run_dir_template, sigma_cap), mustWork = TRUE)
    seed_dirs <- list_seed_dirs(run_dir = run_dir, max_seeds = max_seeds)
    for (seed_dir in seed_dirs) {
      idx <- idx + 1L
      live_dir <- file.path(seed_dir, live_effective_subdir)
      task_rows[[idx]] <- data.frame(
        task_id = idx,
        sigma_cap = sigma_cap,
        run_dir = run_dir,
        seed = basename(seed_dir),
        seed_dir = normalizePath(seed_dir, mustWork = TRUE),
        live_effective_pms_dir = normalizePath(live_dir, mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, task_rows)
}

validate_live_effective_contract <- function(live_dir) {
  manifest_path <- file.path(live_dir, "live_effective_pms_manifest.tsv")
  if (!file.exists(manifest_path)) {
    stop(
      "Missing live-effective-p_ms simulation manifest: ",
      manifest_path,
      ". Run simulation/invivo/cin/generate_live_effective_pms_outputs.R ",
      "for this seed before running the comparison analysis."
    )
  }
  manifest <- read_required_tsv(manifest_path)
  if (!all(c("key", "value") %in% names(manifest))) {
    stop("Invalid live-effective-p_ms manifest schema: ", manifest_path)
  }
  status <- manifest$value[manifest$key == "status"]
  if (!length(status) || !identical(as.character(status[[1]]), "complete")) {
    stop("Live-effective-p_ms simulation is not marked complete: ", manifest_path)
  }
  invisible(manifest_path)
}

read_one_seed_summary <- function(task_row) {
  seed_dir <- as.character(task_row$seed_dir[[1]])
  sigma_cap <- as.character(task_row$sigma_cap[[1]])
  live_dir <- as.character(task_row$live_effective_pms_dir[[1]])
  validate_live_effective_contract(live_dir)

  overall_path <- file.path(live_dir, "live_effective_pms_overall.tsv")
  harvest_path <- file.path(live_dir, "live_effective_pms_harvest_only.tsv")
  overall_tab <- read_required_tsv(overall_path)
  harvest_tab <- read_required_tsv(harvest_path)
  required_cols <- c(
    "summary_scope",
    "p_misseg_parameter",
    "live_weighted_effective_p_ms_mean"
  )
  if (!all(required_cols %in% names(overall_tab))) {
    stop(
      "Missing required columns in ",
      overall_path,
      ": ",
      paste(setdiff(required_cols, names(overall_tab)), collapse = ", ")
    )
  }
  if (!all(required_cols %in% names(harvest_tab))) {
    stop(
      "Missing required columns in ",
      harvest_path,
      ": ",
      paste(setdiff(required_cols, names(harvest_tab)), collapse = ", ")
    )
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
    live_effective_pms_dir = normalizePath(live_dir, mustWork = TRUE),
    p_misseg_parameter = p_misseg_parameter,
    live_cell_p_misseg = live_cell_p_misseg,
    harvest_live_cell_p_misseg = harvest_live_cell_p_misseg,
    abs_diff_live_minus_parameter = live_cell_p_misseg - p_misseg_parameter,
    ratio_live_over_parameter = live_cell_p_misseg / pmax(p_misseg_parameter, 1e-12),
    stringsAsFactors = FALSE
  )
}

process_one_task <- function(task_row) {
  tryCatch(
    {
      message(
        "Analyzing sigma_burden=",
        task_row$sigma_cap[[1]],
        " seed=",
        task_row$seed[[1]]
      )
      out <- read_one_seed_summary(task_row)
      out$task_status <- "ok"
      out$task_error <- ""
      out
    },
    error = function(e) {
      data.frame(
        sigma_cap = as.character(task_row$sigma_cap[[1]]),
        sigma_cap_num = suppressWarnings(as.numeric(task_row$sigma_cap[[1]])),
        seed = as.character(task_row$seed[[1]]),
        seed_dir = as.character(task_row$seed_dir[[1]]),
        live_effective_pms_dir = as.character(task_row$live_effective_pms_dir[[1]]),
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

summarise_by_sigma <- function(seed_summary, sigma_levels) {
  rows <- lapply(
    sigma_levels,
    function(sigma_cap) {
      df <- seed_summary[as.character(seed_summary$sigma_cap) == sigma_cap, , drop = FALSE]
      data.frame(
        sigma_cap = sigma_cap,
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
    }
  )
  do.call(rbind, rows)
}

make_plot_table <- function(seed_summary) {
  rbind(
    data.frame(
      sigma_cap = as.character(seed_summary$sigma_cap),
      seed = seed_summary$seed,
      estimate_type = "p_misseg parameter",
      value = seed_summary$p_misseg_parameter,
      stringsAsFactors = FALSE
    ),
    data.frame(
      sigma_cap = as.character(seed_summary$sigma_cap),
      seed = seed_summary$seed,
      estimate_type = "live-cell effective p_misseg",
      value = seed_summary$live_cell_p_misseg,
      stringsAsFactors = FALSE
    )
  )
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  if (!is.null(argv$estimate_script) || !is.null(argv$rerun_estimate)) {
    stop(
      "--estimate_script and --rerun_estimate are no longer supported by the ",
      "analysis layer. Materialize every seed first with ",
      "simulation/invivo/cin/generate_live_effective_pms_outputs.R."
    )
  }

  sigma_caps <- parse_sigma_caps(argv$sigma_caps)
  run_dir_template <- argv$run_dir_template %||% default_template
  out_dir <- resolve_path_value(argv$out_dir, getwd()) %||% default_out_dir
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  live_effective_subdir <- as.character(
    argv$live_effective_subdir %||% default_live_effective_subdir
  )
  if (!nzchar(trimws(live_effective_subdir))) {
    stop("--live_effective_subdir must not be empty.")
  }
  if (grepl("^(/|[A-Za-z]:[/\\\\])", live_effective_subdir)) {
    stop("--live_effective_subdir must be relative to each seed directory.")
  }

  max_seeds <- NULL
  if (!is.null(argv$max_seeds) && nzchar(trimws(as.character(argv$max_seeds)))) {
    max_seeds <- suppressWarnings(as.integer(argv$max_seeds))
    if (!is.finite(max_seeds) || max_seeds <= 0L) {
      stop("--max_seeds must be a positive integer when supplied.")
    }
  }
  n_cores <- suppressWarnings(as.integer(argv$n_cores %||% 1L))
  if (!is.finite(n_cores) || n_cores <= 0L) {
    stop("--n_cores must be a positive integer.")
  }

  task_manifest <- build_task_manifest(
    sigma_caps = sigma_caps,
    run_dir_template = run_dir_template,
    live_effective_subdir = live_effective_subdir,
    max_seeds = max_seeds
  )
  task_manifest_path <- file.path(
    out_dir,
    "sigma_burden_live_effective_pms_task_manifest.tsv"
  )
  write_tsv(task_manifest, task_manifest_path)

  if (n_cores == 1L || nrow(task_manifest) <= 1L) {
    seed_rows <- lapply(
      seq_len(nrow(task_manifest)),
      function(i) process_one_task(task_manifest[i, , drop = FALSE])
    )
  } else if (.Platform$OS.type == "unix") {
    seed_rows <- parallel::mclapply(
      seq_len(nrow(task_manifest)),
      function(i) process_one_task(task_manifest[i, , drop = FALSE]),
      mc.cores = min(n_cores, nrow(task_manifest)),
      mc.preschedule = FALSE
    )
  } else {
    cl <- parallel::makeCluster(min(n_cores, nrow(task_manifest)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c(
        "read_required_tsv",
        "validate_live_effective_contract",
        "read_one_seed_summary",
        "process_one_task"
      ),
      envir = environment()
    )
    seed_rows <- parallel::parLapplyLB(
      cl,
      seq_len(nrow(task_manifest)),
      function(i, task_tab) process_one_task(task_tab[i, , drop = FALSE]),
      task_tab = task_manifest
    )
  }

  task_results <- do.call(rbind, seed_rows)
  task_results_path <- file.path(
    out_dir,
    "sigma_burden_live_effective_pms_task_results.tsv"
  )
  write_tsv(task_results, task_results_path)
  failed_rows <- task_results[task_results$task_status != "ok", , drop = FALSE]
  if (nrow(failed_rows) > 0L) {
    stop(
      "One or more analysis inputs failed validation. See ",
      task_results_path,
      ". First error: sigma_burden=",
      failed_rows$sigma_cap[[1]],
      " seed=",
      failed_rows$seed[[1]],
      " message=",
      failed_rows$task_error[[1]]
    )
  }

  seed_summary <- task_results
  seed_summary$task_status <- NULL
  seed_summary$task_error <- NULL
  sigma_factor <- factor(seed_summary$sigma_cap, levels = sigma_caps)
  seed_num <- suppressWarnings(as.numeric(sub("^seed", "", seed_summary$seed)))
  seed_summary <- seed_summary[order(sigma_factor, seed_num, seed_summary$seed), , drop = FALSE]

  summary_by_sigma <- summarise_by_sigma(seed_summary, sigma_caps)
  plot_tab <- make_plot_table(seed_summary)
  plot_sigma <- factor(plot_tab$sigma_cap, levels = sigma_caps)
  plot_seed_num <- suppressWarnings(as.numeric(sub("^seed", "", plot_tab$seed)))
  plot_tab <- plot_tab[
    order(plot_sigma, plot_tab$estimate_type, plot_seed_num, plot_tab$seed),
    ,
    drop = FALSE
  ]

  by_seed_path <- file.path(
    out_dir,
    "sigma_burden_live_effective_pms_by_seed.tsv"
  )
  summary_path <- file.path(
    out_dir,
    "sigma_burden_live_effective_pms_summary.tsv"
  )
  plot_path <- file.path(
    out_dir,
    "sigma_burden_p_misseg_vs_live_cell_plot.tsv"
  )
  write_tsv(seed_summary, by_seed_path)
  write_tsv(summary_by_sigma, summary_path)
  write_tsv(plot_tab, plot_path)

  output_tables <- list(
    sigma_burden_live_effective_pms_task_manifest = task_manifest,
    sigma_burden_live_effective_pms_task_results = task_results,
    sigma_burden_live_effective_pms_by_seed = seed_summary,
    sigma_burden_live_effective_pms_summary = summary_by_sigma,
    sigma_burden_p_misseg_vs_live_cell_plot = plot_tab
  )
  schema <- do.call(
    rbind,
    lapply(names(output_tables), function(name) {
      tab <- output_tables[[name]]
      data.frame(
        table = paste0(name, ".tsv"),
        rows = nrow(tab),
        columns = paste(names(tab), collapse = ","),
        stringsAsFactors = FALSE
      )
    })
  )
  schema_path <- file.path(out_dir, "live_effective_pms_comparison_schema.tsv")
  write_tsv(schema, schema_path)
  manifest <- data.frame(
    key = c(
      "schema_version",
      "status",
      "run_dir_template",
      "sigma_caps",
      "live_effective_subdir",
      "output_dir",
      "generated_at"
    ),
    value = c(
      "o2sd-live-effective-pms-analysis-v1",
      "complete",
      run_dir_template,
      paste(sigma_caps, collapse = ","),
      live_effective_subdir,
      out_dir,
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(out_dir, "live_effective_pms_comparison_manifest.tsv")
  write_tsv(manifest, manifest_path)

  message("Wrote by-seed analysis: ", by_seed_path)
  message("Wrote sigma summary: ", summary_path)
  message("Wrote plot-ready data: ", plot_path)
  message("Wrote analysis manifest: ", manifest_path)
}

if (sys.nframe() == 0) {
  main()
}
