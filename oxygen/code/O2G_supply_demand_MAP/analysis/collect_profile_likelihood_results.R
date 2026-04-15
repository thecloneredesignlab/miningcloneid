#!/usr/bin/env Rscript

# Usage:
#   Rscript oxygen/code/O2G_supply_demand_MAP/analysis/collect_profile_likelihood_results.R \
#     --output_root=/abs/path/to/profile_output_root

.o2sd_profile_collect_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2sd_profile_collect_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_profile_collect_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args

resolve_path_value <- function(path_value, base_dir) {
  txt <- path_value
  if (is.null(txt) || !length(txt)) return(NULL)
  txt <- as.character(txt[[1]])
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  if (startsWith(txt, "~")) {
    return(normalizePath(path.expand(txt), mustWork = FALSE))
  }
  if (grepl("^(/|[A-Za-z]:[/\\\\])", txt)) {
    return(normalizePath(txt, mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, txt), mustWork = FALSE)
}

read_required_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file was not found: ", path)
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_table_tsv <- function(tab, path) {
  utils::write.table(
    tab,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(path)
}

escape_markdown_cell <- function(x) {
  txt <- ifelse(is.na(x), "", as.character(x))
  txt <- gsub("\\|", "\\\\|", txt)
  txt <- gsub("\n", " ", txt, fixed = TRUE)
  txt
}

data_frame_to_markdown <- function(tab) {
  if (!ncol(tab)) return("No columns.\n")
  cols <- names(tab)
  header <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
  if (!nrow(tab)) return(paste(c(header, sep), collapse = "\n"))
  rows <- apply(
    tab,
    1,
    function(row) {
      paste0("| ", paste(vapply(row, escape_markdown_cell, character(1)), collapse = " | "), " |")
    }
  )
  paste(c(header, sep, rows), collapse = "\n")
}

first_numeric_or_na <- function(tab, col_name) {
  if (!(col_name %in% names(tab)) || !nrow(tab)) return(NA_real_)
  suppressWarnings(as.numeric(tab[[col_name]][[1]]))
}

with_delta_from_metric <- function(tab, metric_col, baseline_metric, delta_col) {
  out <- tab
  if (!(metric_col %in% names(out)) || !is.finite(baseline_metric)) {
    out[[delta_col]] <- NA_real_
    return(out)
  }
  out[[delta_col]] <- suppressWarnings(as.numeric(out[[metric_col]])) - baseline_metric
  out
}

direction_profile_order <- function(direction_name, direction_df) {
  vals <- as.numeric(direction_df$fixed_value)
  if (identical(direction_name, "decreasing")) {
    ord <- order(-vals, direction_df$path_index %||% seq_along(vals), na.last = TRUE)
  } else if (identical(direction_name, "increasing")) {
    ord <- order(vals, direction_df$path_index %||% seq_along(vals), na.last = TRUE)
  } else {
    ord <- order(vals, direction_df$path_index %||% seq_along(vals), na.last = TRUE)
  }
  as.integer(ord)
}

estimate_direction_ci_from_delta <- function(direction_df, baseline_value, delta_col, ci_delta_threshold, direction_name, direction_summary_row) {
  if (nrow(direction_summary_row) != 1L) {
    stop("Expected exactly one direction_summary row for ", direction_name)
  }
  if (isTRUE(direction_summary_row$start_blocked[[1]])) {
    return(list(
      ci_value = NA_real_,
      ci_status = "blocked_at_start_boundary"
    ))
  }
  if (!(delta_col %in% names(direction_df))) {
    return(list(
      ci_value = NA_real_,
      ci_status = paste0("missing_metric:", delta_col)
    ))
  }

  complete <- direction_df[is.finite(direction_df[[delta_col]]), , drop = FALSE]
  if (!nrow(complete)) {
    return(list(
      ci_value = NA_real_,
      ci_status = as.character(direction_summary_row$termination_reason[[1]] %||% "no_complete_steps")
    ))
  }
  complete <- complete[direction_profile_order(direction_name, complete), , drop = FALSE]

  deltas <- as.numeric(complete[[delta_col]])
  crossing_idx <- which(is.finite(deltas) & deltas >= ci_delta_threshold)
  if (length(crossing_idx) == 0L) {
    return(list(
      ci_value = NA_real_,
      ci_status = paste0("threshold_not_reached:", as.character(direction_summary_row$termination_reason[[1]]))
    ))
  }

  k <- crossing_idx[[1]]
  x1 <- as.numeric(complete$fixed_value[[k]])
  y1 <- as.numeric(deltas[[k]])
  if (k == 1L) {
    x0 <- as.numeric(baseline_value)
    y0 <- 0
  } else {
    x0 <- as.numeric(complete$fixed_value[[k - 1L]])
    y0 <- as.numeric(deltas[[k - 1L]])
  }

  if (!is.finite(y1) || !is.finite(x1)) {
    return(list(ci_value = NA_real_, ci_status = "crossing_non_finite"))
  }
  if (!is.finite(y0) || !is.finite(x0) || abs(y1 - y0) <= .Machine$double.eps) {
    return(list(ci_value = x1, ci_status = "crossing_no_interpolation"))
  }

  frac <- (ci_delta_threshold - y0) / (y1 - y0)
  frac <- min(max(frac, 0), 1)
  ci_value <- x0 + frac * (x1 - x0)
  list(ci_value = ci_value, ci_status = "interpolated_crossing")
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  output_root <- resolve_path_value(argv$output_root, getwd())
  if (is.null(output_root)) {
    stop("--output_root is required.")
  }
  output_root <- normalizePath(output_root, mustWork = TRUE)

  target_path <- file.path(output_root, "parameter_targets.tsv")
  submission_manifest_path <- file.path(output_root, "submission_manifest.tsv")
  parameter_dirs <- list.dirs(output_root, recursive = FALSE, full.names = TRUE)
  parameter_dirs <- parameter_dirs[file.exists(file.path(parameter_dirs, "profile_path_combined.tsv"))]
  if (length(parameter_dirs) == 0L) {
    stop("No parameter directories with profile_path_combined.tsv were found under ", output_root)
  }

  combined_list <- lapply(parameter_dirs, function(d) read_required_tsv(file.path(d, "profile_path_combined.tsv")))
  manifest_list <- lapply(parameter_dirs, function(d) read_required_tsv(file.path(d, "parameter_manifest.tsv")))
  direction_summary_list <- lapply(parameter_dirs, function(d) read_required_tsv(file.path(d, "direction_summary.tsv")))
  point_summary_list <- lapply(
    parameter_dirs,
    function(d) {
      p <- file.path(d, "profile_point_summary.tsv")
      if (file.exists(p)) read_required_tsv(p) else NULL
    }
  )
  relation_list <- lapply(
    parameter_dirs,
    function(d) {
      p <- file.path(d, "parameter_relation_long.tsv")
      if (file.exists(p)) read_required_tsv(p) else NULL
    }
  )

  profile_all <- do.call(rbind, combined_list)
  manifest_all <- do.call(rbind, manifest_list)
  direction_summary_all <- do.call(rbind, direction_summary_list)
  point_summary_all <- do.call(rbind, point_summary_list[!vapply(point_summary_list, is.null, logical(1))])
  relation_all <- do.call(rbind, relation_list[!vapply(relation_list, is.null, logical(1))])

  split_profiles <- split(profile_all, profile_all$param_symbol)
  split_point_summaries <- if (!is.null(point_summary_all) && nrow(point_summary_all) > 0L) split(point_summary_all, point_summary_all$param_symbol) else list()
  summary_rows <- lapply(
    names(split_profiles),
    function(param_symbol) {
      tab <- split_profiles[[param_symbol]]
      point_tab <- if (param_symbol %in% names(split_point_summaries)) split_point_summaries[[param_symbol]] else NULL
      manifest_row <- manifest_all[manifest_all$param_symbol == param_symbol, , drop = FALSE]
      if (nrow(manifest_row) != 1L) {
        stop("Expected exactly one parameter_manifest row for ", param_symbol)
      }
      dir_summary <- direction_summary_all[direction_summary_all$param_symbol == param_symbol, , drop = FALSE]
      if (nrow(dir_summary) != 2L) {
        stop("Expected two direction_summary rows for ", param_symbol)
      }

      baseline_row <- tab[tab$direction == "baseline", , drop = FALSE]
      if (nrow(baseline_row) != 1L) {
        stop("Expected exactly one baseline row for ", param_symbol)
      }
      complete <- tab[is.finite(tab$objective), , drop = FALSE]
      best_idx <- which.min(complete$objective)
      best_row <- complete[best_idx, , drop = FALSE]

      decreasing_tab <- tab[tab$direction == "decreasing", , drop = FALSE]
      increasing_tab <- tab[tab$direction == "increasing", , drop = FALSE]
      decreasing_summary <- dir_summary[dir_summary$direction == "decreasing", , drop = FALSE]
      increasing_summary <- dir_summary[dir_summary$direction == "increasing", , drop = FALSE]
      point_decreasing <- if (!is.null(point_tab)) point_tab[point_tab$direction == "decreasing", , drop = FALSE] else data.frame()
      point_increasing <- if (!is.null(point_tab)) point_tab[point_tab$direction == "increasing", , drop = FALSE] else data.frame()
      raw_neg2loglik_ci_threshold <- if ("raw_neg2loglik_ci_threshold" %in% names(manifest_row)) {
        as.numeric(manifest_row$raw_neg2loglik_ci_threshold[[1]])
      } else {
        3.84
      }
      decreasing_tab_burden_raw <- with_delta_from_metric(decreasing_tab, "objective_burden_neg2loglik_raw", first_numeric_or_na(baseline_row, "objective_burden_neg2loglik_raw"), "delta_objective_burden_neg2loglik_raw")
      increasing_tab_burden_raw <- with_delta_from_metric(increasing_tab, "objective_burden_neg2loglik_raw", first_numeric_or_na(baseline_row, "objective_burden_neg2loglik_raw"), "delta_objective_burden_neg2loglik_raw")
      decreasing_tab_ploidy_raw <- with_delta_from_metric(decreasing_tab, "objective_ploidy_neg2loglik_raw", first_numeric_or_na(baseline_row, "objective_ploidy_neg2loglik_raw"), "delta_objective_ploidy_neg2loglik_raw")
      increasing_tab_ploidy_raw <- with_delta_from_metric(increasing_tab, "objective_ploidy_neg2loglik_raw", first_numeric_or_na(baseline_row, "objective_ploidy_neg2loglik_raw"), "delta_objective_ploidy_neg2loglik_raw")

      lower_ci <- estimate_direction_ci_from_delta(
        direction_df = decreasing_tab,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_vs_baseline",
        ci_delta_threshold = manifest_row$ci_delta_threshold[[1]],
        direction_name = "decreasing",
        direction_summary_row = decreasing_summary
      )
      upper_ci <- estimate_direction_ci_from_delta(
        direction_df = increasing_tab,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_vs_baseline",
        ci_delta_threshold = manifest_row$ci_delta_threshold[[1]],
        direction_name = "increasing",
        direction_summary_row = increasing_summary
      )
      burden_best_lower_ci <- estimate_direction_ci_from_delta(
        direction_df = decreasing_tab_burden_raw,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_burden_neg2loglik_raw",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "decreasing",
        direction_summary_row = decreasing_summary
      )
      burden_best_upper_ci <- estimate_direction_ci_from_delta(
        direction_df = increasing_tab_burden_raw,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_burden_neg2loglik_raw",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "increasing",
        direction_summary_row = increasing_summary
      )
      ploidy_best_lower_ci <- estimate_direction_ci_from_delta(
        direction_df = decreasing_tab_ploidy_raw,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_ploidy_neg2loglik_raw",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "decreasing",
        direction_summary_row = decreasing_summary
      )
      ploidy_best_upper_ci <- estimate_direction_ci_from_delta(
        direction_df = increasing_tab_ploidy_raw,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_ploidy_neg2loglik_raw",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "increasing",
        direction_summary_row = increasing_summary
      )
      burden_all_lower_ci_near <- estimate_direction_ci_from_delta(
        direction_df = point_decreasing,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_burden_neg2loglik_raw_max",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "decreasing",
        direction_summary_row = decreasing_summary
      )
      burden_all_lower_ci_far <- estimate_direction_ci_from_delta(
        direction_df = point_decreasing,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_burden_neg2loglik_raw_min",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "decreasing",
        direction_summary_row = decreasing_summary
      )
      burden_all_upper_ci_near <- estimate_direction_ci_from_delta(
        direction_df = point_increasing,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_burden_neg2loglik_raw_max",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "increasing",
        direction_summary_row = increasing_summary
      )
      burden_all_upper_ci_far <- estimate_direction_ci_from_delta(
        direction_df = point_increasing,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_burden_neg2loglik_raw_min",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "increasing",
        direction_summary_row = increasing_summary
      )
      ploidy_all_lower_ci_near <- estimate_direction_ci_from_delta(
        direction_df = point_decreasing,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_ploidy_neg2loglik_raw_max",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "decreasing",
        direction_summary_row = decreasing_summary
      )
      ploidy_all_lower_ci_far <- estimate_direction_ci_from_delta(
        direction_df = point_decreasing,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_ploidy_neg2loglik_raw_min",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "decreasing",
        direction_summary_row = decreasing_summary
      )
      ploidy_all_upper_ci_near <- estimate_direction_ci_from_delta(
        direction_df = point_increasing,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_ploidy_neg2loglik_raw_max",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "increasing",
        direction_summary_row = increasing_summary
      )
      ploidy_all_upper_ci_far <- estimate_direction_ci_from_delta(
        direction_df = point_increasing,
        baseline_value = manifest_row$baseline_value[[1]],
        delta_col = "delta_objective_ploidy_neg2loglik_raw_min",
        ci_delta_threshold = raw_neg2loglik_ci_threshold,
        direction_name = "increasing",
        direction_summary_row = increasing_summary
      )

      data.frame(
        param_index = as.integer(manifest_row$param_index[[1]]),
        param_symbol = as.character(param_symbol),
        lower_bound = as.numeric(manifest_row$lower_bound[[1]]),
        upper_bound = as.numeric(manifest_row$upper_bound[[1]]),
        baseline_value = as.numeric(manifest_row$baseline_value[[1]]),
        baseline_objective = as.numeric(manifest_row$baseline_objective[[1]]),
        baseline_objective_burden = first_numeric_or_na(baseline_row, "objective_burden"),
        baseline_objective_ploidy = first_numeric_or_na(baseline_row, "objective_ploidy"),
        baseline_objective_burden_neg2loglik_raw = first_numeric_or_na(baseline_row, "objective_burden_neg2loglik_raw"),
        baseline_objective_ploidy_neg2loglik_raw = first_numeric_or_na(baseline_row, "objective_ploidy_neg2loglik_raw"),
        use_soft_prior_for_profile = isTRUE(manifest_row$use_soft_prior_for_profile[[1]]),
        lambda_prior_for_profile = as.numeric(manifest_row$lambda_prior_for_profile[[1]]),
        max_steps_per_direction = as.integer(manifest_row$max_steps_per_direction[[1]]),
        seeds_per_step = as.integer(manifest_row$seeds_per_step[[1]]),
        n_cores = as.integer(manifest_row$n_cores[[1]]),
        ci_delta_threshold = as.numeric(manifest_row$ci_delta_threshold[[1]]),
        raw_neg2loglik_ci_threshold = as.numeric(raw_neg2loglik_ci_threshold),
        steps_completed_decreasing = as.integer(decreasing_summary$steps_completed[[1]]),
        steps_completed_increasing = as.integer(increasing_summary$steps_completed[[1]]),
        decreasing_termination_reason = as.character(decreasing_summary$termination_reason[[1]]),
        increasing_termination_reason = as.character(increasing_summary$termination_reason[[1]]),
        lower_profile_ci = as.numeric(lower_ci$ci_value),
        lower_profile_ci_status = as.character(lower_ci$ci_status),
        upper_profile_ci = as.numeric(upper_ci$ci_value),
        upper_profile_ci_status = as.character(upper_ci$ci_status),
        best_profile_value = as.numeric(best_row$fixed_value[[1]]),
        best_profile_objective = as.numeric(best_row$objective[[1]]),
        best_profile_objective_burden = first_numeric_or_na(best_row, "objective_burden"),
        best_profile_objective_ploidy = first_numeric_or_na(best_row, "objective_ploidy"),
        best_profile_objective_burden_neg2loglik_raw = first_numeric_or_na(best_row, "objective_burden_neg2loglik_raw"),
        best_profile_objective_ploidy_neg2loglik_raw = first_numeric_or_na(best_row, "objective_ploidy_neg2loglik_raw"),
        lower_profile_ci_burden_neg2loglik_raw_best = as.numeric(burden_best_lower_ci$ci_value),
        lower_profile_ci_burden_neg2loglik_raw_best_status = as.character(burden_best_lower_ci$ci_status),
        upper_profile_ci_burden_neg2loglik_raw_best = as.numeric(burden_best_upper_ci$ci_value),
        upper_profile_ci_burden_neg2loglik_raw_best_status = as.character(burden_best_upper_ci$ci_status),
        lower_profile_ci_ploidy_neg2loglik_raw_best = as.numeric(ploidy_best_lower_ci$ci_value),
        lower_profile_ci_ploidy_neg2loglik_raw_best_status = as.character(ploidy_best_lower_ci$ci_status),
        upper_profile_ci_ploidy_neg2loglik_raw_best = as.numeric(ploidy_best_upper_ci$ci_value),
        upper_profile_ci_ploidy_neg2loglik_raw_best_status = as.character(ploidy_best_upper_ci$ci_status),
        lower_profile_ci_burden_neg2loglik_raw_allseeds_near = as.numeric(burden_all_lower_ci_near$ci_value),
        lower_profile_ci_burden_neg2loglik_raw_allseeds_near_status = as.character(burden_all_lower_ci_near$ci_status),
        lower_profile_ci_burden_neg2loglik_raw_allseeds_far = as.numeric(burden_all_lower_ci_far$ci_value),
        lower_profile_ci_burden_neg2loglik_raw_allseeds_far_status = as.character(burden_all_lower_ci_far$ci_status),
        upper_profile_ci_burden_neg2loglik_raw_allseeds_near = as.numeric(burden_all_upper_ci_near$ci_value),
        upper_profile_ci_burden_neg2loglik_raw_allseeds_near_status = as.character(burden_all_upper_ci_near$ci_status),
        upper_profile_ci_burden_neg2loglik_raw_allseeds_far = as.numeric(burden_all_upper_ci_far$ci_value),
        upper_profile_ci_burden_neg2loglik_raw_allseeds_far_status = as.character(burden_all_upper_ci_far$ci_status),
        lower_profile_ci_ploidy_neg2loglik_raw_allseeds_near = as.numeric(ploidy_all_lower_ci_near$ci_value),
        lower_profile_ci_ploidy_neg2loglik_raw_allseeds_near_status = as.character(ploidy_all_lower_ci_near$ci_status),
        lower_profile_ci_ploidy_neg2loglik_raw_allseeds_far = as.numeric(ploidy_all_lower_ci_far$ci_value),
        lower_profile_ci_ploidy_neg2loglik_raw_allseeds_far_status = as.character(ploidy_all_lower_ci_far$ci_status),
        upper_profile_ci_ploidy_neg2loglik_raw_allseeds_near = as.numeric(ploidy_all_upper_ci_near$ci_value),
        upper_profile_ci_ploidy_neg2loglik_raw_allseeds_near_status = as.character(ploidy_all_upper_ci_near$ci_status),
        upper_profile_ci_ploidy_neg2loglik_raw_allseeds_far = as.numeric(ploidy_all_upper_ci_far$ci_value),
        upper_profile_ci_ploidy_neg2loglik_raw_allseeds_far_status = as.character(ploidy_all_upper_ci_far$ci_status),
        best_profile_direction = as.character(best_row$direction[[1]]),
        best_profile_step_index = as.integer(best_row$step_index[[1]]),
        best_profile_seed = as.integer(best_row$best_seed[[1]]),
        stringsAsFactors = FALSE
      )
    }
  )
  parameter_summary <- do.call(rbind, summary_rows)
  parameter_summary <- parameter_summary[order(parameter_summary$param_index), , drop = FALSE]

  profile_all_path <- file.path(output_root, "profile_likelihood_all.tsv")
  direction_summary_path <- file.path(output_root, "profile_direction_summary_all.tsv")
  point_summary_path <- file.path(output_root, "profile_point_summary_all.tsv")
  relation_all_path <- file.path(output_root, "profile_parameter_relations_all.tsv")
  parameter_summary_path <- file.path(output_root, "profile_parameter_summary.tsv")
  report_path <- file.path(output_root, "profile_likelihood_report.md")

  write_table_tsv(profile_all, profile_all_path)
  write_table_tsv(direction_summary_all, direction_summary_path)
  if (!is.null(point_summary_all) && nrow(point_summary_all) > 0L) write_table_tsv(point_summary_all, point_summary_path)
  write_table_tsv(parameter_summary, parameter_summary_path)
  if (!is.null(relation_all) && nrow(relation_all) > 0L) {
    write_table_tsv(relation_all, relation_all_path)
  }

  target_table <- if (file.exists(target_path)) read_required_tsv(target_path) else NULL
  submission_manifest <- if (file.exists(submission_manifest_path)) read_required_tsv(submission_manifest_path) else NULL

  missing_targets <- NULL
  if (!is.null(target_table) && all(c("param_index", "param_symbol") %in% names(target_table))) {
    missing_targets <- target_table[!target_table$param_symbol %in% parameter_summary$param_symbol, , drop = FALSE]
  }

  report_lines <- c(
    "# Profile Likelihood Report",
    "",
    paste0("- Output root: `", normalizePath(output_root, mustWork = FALSE), "`"),
    paste0("- Generated on: `", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "`"),
    ""
  )

  if (!is.null(submission_manifest) && nrow(submission_manifest) > 0L) {
    report_lines <- c(
      report_lines,
      "## Submission Manifest",
      "",
      data_frame_to_markdown(submission_manifest),
      ""
    )
  }

  if (!is.null(target_table) && nrow(target_table) > 0L) {
    report_lines <- c(
      report_lines,
      "## Target Parameters",
      "",
      data_frame_to_markdown(target_table),
      ""
    )
  }

  report_lines <- c(
    report_lines,
    "## Parameter Summary",
    "",
    data_frame_to_markdown(parameter_summary),
    "",
    "## Direction Summary",
    "",
    data_frame_to_markdown(direction_summary_all),
    ""
  )

  if (!is.null(missing_targets) && nrow(missing_targets) > 0L) {
    report_lines <- c(
      report_lines,
      "## Missing Parameter Jobs",
      "",
      data_frame_to_markdown(missing_targets),
      ""
    )
  }

  report_lines <- c(
    report_lines,
    "## Combined Profile Table",
    "",
    data_frame_to_markdown(profile_all),
    ""
  )
  writeLines(report_lines, con = report_path, sep = "\n")

  message("Wrote combined profile table: ", profile_all_path)
  message("Wrote direction summary: ", direction_summary_path)
  message("Wrote parameter summary: ", parameter_summary_path)
  if (!is.null(relation_all) && nrow(relation_all) > 0L) {
    message("Wrote parameter relations: ", relation_all_path)
  }
  message("Wrote report: ", report_path)
}

if (sys.nframe() == 0) {
  main()
}
