#!/usr/bin/env Rscript

.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)) else getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)
parse_args <- o2fr_parse_args
first_non_null <- o2fr_null_coalesce
as_num <- o2fr_as_num
as_chr <- o2fr_as_chr
read_tsv <- function(path) o2fr_read_tsv(path, optional = FALSE)

metric_value <- function(summary_df, metric_names) {
  if (!is.data.frame(summary_df) || !all(c("metric", "value") %in% names(summary_df))) {
    return(NA_real_)
  }
  idx <- match(metric_names, as.character(summary_df$metric), nomatch = 0L)
  idx <- idx[idx > 0L]
  if (!length(idx)) return(NA_real_)
  suppressWarnings(as.numeric(summary_df$value[[idx[[1]]]]))
}

seed_number_from_dir <- function(seed_dir) {
  suppressWarnings(as.integer(sub("^seed", "", basename(seed_dir))))
}

select_prediction_path <- function(seed_dir, horizon) {
  hz <- as.integer(round(horizon))
  candidates <- c(
    file.path(seed_dir, "simulation", "invivo", sprintf("predict_ploidy_weighted_mean_0_%sday.tsv", hz)),
    file.path(seed_dir, "simulation", "invivo", "ploidy_weighted_mean_timecourse.tsv"),
    file.path(seed_dir, "viz", sprintf("predict_ploidy_weighted_mean_0_%sday.tsv", hz)),
    file.path(seed_dir, "viz", "ploidy_weighted_mean_timecourse.tsv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1]] else NA_character_
}

long_term_value <- function(pred_path, horizon, cohort = "2N", dose = "ALL") {
  if (!is.finite(horizon) || horizon <= 0) {
    return(list(value = NA_real_, day = NA_real_, rows = 0L, status = "invalid_horizon"))
  }
  if (!is.character(pred_path) || !nzchar(pred_path) || !file.exists(pred_path)) {
    return(list(value = NA_real_, day = NA_real_, rows = 0L, status = "missing_prediction"))
  }

  tab <- tryCatch(read_tsv(pred_path), error = function(e) NULL)
  if (is.null(tab) || !nrow(tab)) {
    return(list(value = NA_real_, day = NA_real_, rows = 0L, status = "empty_prediction"))
  }
  if (!"day" %in% names(tab)) {
    return(list(value = NA_real_, day = NA_real_, rows = 0L, status = "missing_day"))
  }
  tab$day <- suppressWarnings(as.numeric(tab$day))
  tab <- tab[is.finite(tab$day), , drop = FALSE]
  if (!nrow(tab)) {
    return(list(value = NA_real_, day = NA_real_, rows = 0L, status = "no_finite_day"))
  }

  if (!identical(toupper(cohort), "ALL") && "cohort" %in% names(tab)) {
    tab <- tab[as.character(tab$cohort) == cohort, , drop = FALSE]
  }
  if (!identical(toupper(dose), "ALL") && "dose" %in% names(tab)) {
    tab <- tab[as.character(tab$dose) == dose, , drop = FALSE]
  }
  if (!nrow(tab)) {
    return(list(value = NA_real_, day = NA_real_, rows = 0L, status = "no_rows_after_filter"))
  }

  day_use <- max(tab$day[tab$day <= horizon], na.rm = TRUE)
  if (!is.finite(day_use)) day_use <- max(tab$day, na.rm = TRUE)
  tab <- tab[abs(tab$day - day_use) <= 1e-8, , drop = FALSE]
  if (!nrow(tab)) {
    return(list(value = NA_real_, day = day_use, rows = 0L, status = "no_rows_at_day"))
  }

  value_col <- intersect(c("weighted_mean_N", "weighted_mean_endpoint", "weighted_mean_ploidy"), names(tab))
  if (!length(value_col)) {
    return(list(value = NA_real_, day = day_use, rows = nrow(tab), status = "missing_weighted_mean"))
  }
  vals <- suppressWarnings(as.numeric(tab[[value_col[[1]]]]))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(list(value = NA_real_, day = day_use, rows = nrow(tab), status = "no_finite_weighted_mean"))
  }
  list(value = max(vals, na.rm = TRUE), day = day_use, rows = nrow(tab), status = "ok")
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  invivo_dir <- normalizePath(as_chr(argv$invivo_dir), mustWork = FALSE)
  if (!dir.exists(invivo_dir)) stop("invivo_dir does not exist: ", invivo_dir, call. = FALSE)
  out_dir <- normalizePath(
    as_chr(first_non_null(argv$simulation_dir, argv$out_dir), file.path(invivo_dir, "simulation", "fit_results", "long_ploidy")),
    mustWork = FALSE
  )
  horizon <- as_num(argv$horizon, 1000)
  cohort <- as_chr(argv$cohort, "2N")
  dose <- as_chr(argv$dose, "ALL")
  seed_dirs <- list.dirs(invivo_dir, recursive = FALSE, full.names = TRUE)
  seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
  seed_nums <- vapply(seed_dirs, seed_number_from_dir, integer(1))
  keep <- is.finite(seed_nums)
  seed_dirs <- seed_dirs[keep]
  seed_nums <- seed_nums[keep]
  ord <- order(seed_nums)
  seed_dirs <- seed_dirs[ord]
  seed_nums <- seed_nums[ord]
  if (!length(seed_dirs)) stop("No seed directories found under: ", invivo_dir, call. = FALSE)

  rows <- lapply(seq_along(seed_dirs), function(i) {
    seed_dir <- seed_dirs[[i]]
    summary_path <- file.path(seed_dir, "fit_summary.tsv")
    summary_df <- if (file.exists(summary_path)) tryCatch(read_tsv(summary_path), error = function(e) NULL) else NULL
    objective <- metric_value(summary_df, c("objective", "optimizer_local_objective", "optimizer_deoptim_objective"))
    pred_path <- select_prediction_path(seed_dir, horizon)
    lt <- long_term_value(pred_path, horizon = horizon, cohort = cohort, dose = dose)
    data.frame(
      seed = seed_nums[[i]], seed_dir = normalizePath(seed_dir, mustWork = FALSE),
      objective = objective, long_term_day = as.numeric(lt$day),
      long_term_weighted_mean_N = as.numeric(lt$value), long_term_rows = as.integer(lt$rows),
      long_term_status = as.character(lt$status), prediction_path = as.character(pred_path),
      selection_horizon = horizon, selection_cohort = cohort, selection_dose = dose,
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, rows)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, "invivo_long_ploidy_metrics.tsv")
  utils::write.table(metrics, path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(stage = "simulation", file = basename(path), stringsAsFactors = FALSE), file.path(out_dir, "simulation_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote long-ploidy simulation metrics: ", path)
  invisible(out_dir)
}

if (sys.nframe() == 0L) main()
