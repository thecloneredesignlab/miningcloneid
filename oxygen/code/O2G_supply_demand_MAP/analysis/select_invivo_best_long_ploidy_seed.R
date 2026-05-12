#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos < 0) {
      out[[kv]] <- TRUE
    } else {
      key <- substr(kv, 1, pos - 1)
      val <- substr(kv, pos + 1, nchar(kv))
      out[[key]] <- val
    }
  }
  out
}

first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric(first_non_null(x, default)))
  if (!is.finite(val)) default else val
}

as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer(first_non_null(x, default)))
  if (!is.finite(val)) default else val
}

as_chr <- function(x, default = "") {
  val <- as.character(first_non_null(x, default))
  if (!length(val) || !nzchar(val[[1]])) default else val[[1]]
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

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
  if (!dir.exists(invivo_dir)) {
    stop("invivo_dir does not exist: ", invivo_dir, call. = FALSE)
  }
  out_tsv <- normalizePath(as_chr(argv$out_tsv, file.path(invivo_dir, "best_long_ploidy_gt2_seed.tsv")), mustWork = FALSE)
  best_dir_file <- normalizePath(as_chr(argv$best_dir_file, file.path(invivo_dir, "best_long_ploidy_gt2_seed.dir")), mustWork = FALSE)
  horizon <- as_num(argv$horizon, 1000)
  threshold_n <- as_num(argv$threshold_N, 44)
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
  if (!length(seed_dirs)) {
    stop("No seed directories found under: ", invivo_dir, call. = FALSE)
  }

  rows <- lapply(seq_along(seed_dirs), function(i) {
    seed_dir <- seed_dirs[[i]]
    summary_path <- file.path(seed_dir, "fit_summary.tsv")
    summary_df <- if (file.exists(summary_path)) {
      tryCatch(read_tsv(summary_path), error = function(e) NULL)
    } else {
      NULL
    }
    objective <- metric_value(
      summary_df,
      c("objective", "optimizer_local_objective", "optimizer_deoptim_objective")
    )
    pred_path <- select_prediction_path(seed_dir, horizon)
    lt <- long_term_value(pred_path, horizon = horizon, cohort = cohort, dose = dose)
    data.frame(
      seed = seed_nums[[i]],
      seed_dir = normalizePath(seed_dir, mustWork = FALSE),
      objective = objective,
      long_term_day = as.numeric(lt$day),
      long_term_weighted_mean_N = as.numeric(lt$value),
      long_term_rows = as.integer(lt$rows),
      long_term_status = as.character(lt$status),
      eligible = is.finite(objective) &&
        identical(lt$status, "ok") &&
        is.finite(as.numeric(lt$value)) &&
        as.numeric(lt$value) > threshold_n,
      prediction_path = as.character(pred_path),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, rows)
  res$selected <- FALSE
  res$selection_horizon <- horizon
  res$selection_threshold_N <- threshold_n
  res$selection_cohort <- cohort
  res$selection_dose <- dose

  eligible <- res[res$eligible %in% TRUE, , drop = FALSE]
  if (!nrow(eligible)) {
    dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(res, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    stop(
      "No eligible seed found with long-term weighted_mean_N > ", threshold_n,
      " at horizon <= ", horizon, " for cohort=", cohort, ", dose=", dose,
      ". Summary written to: ", out_tsv,
      call. = FALSE
    )
  }
  eligible <- eligible[order(eligible$objective, -eligible$long_term_weighted_mean_N, eligible$seed), , drop = FALSE]
  best_seed <- eligible$seed[[1]]
  res$selected[res$seed == best_seed] <- TRUE

  dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(res, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(res$seed_dir[res$selected][[1]], con = best_dir_file)

  message("Selected in vivo seed: seed", best_seed)
  message("  seed_dir: ", res$seed_dir[res$selected][[1]])
  message("  objective: ", signif(res$objective[res$selected][[1]], 8))
  message("  long_term_weighted_mean_N: ", signif(res$long_term_weighted_mean_N[res$selected][[1]], 8))
  message("  summary: ", out_tsv)
  message("  best_dir_file: ", best_dir_file)
}

main()
