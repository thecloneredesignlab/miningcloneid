#!/usr/bin/env Rscript

# Canonical analysis consumer for the combined FixO2 eigen-attractor workflow.

SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

WORKFLOW_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
UTIL_DIR <- file.path(WORKFLOW_DIR, "util")
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_path_utils.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_bpf_cli_utils.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_shared.R"))
source(file.path(UTIL_DIR, "o2_supply_demand_map_combined_fixo2_eigen_utils.R"))

read_tsv <- function(path, required = TRUE) {
  path <- normalizePath(path.expand(path), mustWork = FALSE)
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Input table does not exist: ", path, call. = FALSE)
    return(data.frame())
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- o2cfe_write_tsv

require_columns <- function(x, cols, table_name) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) {
    stop(table_name, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

as_numeric_column <- o2sd_numeric

seed_number_from_id <- function(seed_id) {
  seed_id <- as.character(seed_id)
  out <- sub("^.*?(\\d+)$", "\\1", seed_id)
  out[identical(out, seed_id) & !grepl("\\d+$", seed_id)] <- NA_character_
  suppressWarnings(as.integer(out))
}

first_finite <- function(x, default = NA_real_) {
  x <- x[is.finite(x)]
  if (length(x)) x[[1L]] else default
}

safe_fraction <- function(numerator, denominator) {
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  numerator / denominator
}

default_curve_table <- function(repo_root = bpf_repo_root(SCRIPT_DIR)) {
  file.path(
    bpf_dense_grid_result_root(repo_root),
    "dense-grid_monotonicity_classification",
    "tables",
    "fixed_o2_ploidy_monotonicity_regression_curves.tsv"
  )
}

default_by_seed_table <- function(repo_root = bpf_repo_root(SCRIPT_DIR)) {
  file.path(
    bpf_dense_grid_result_root(repo_root),
    "dense-grid_monotonicity_classification",
    "tables",
    "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv"
  )
}

default_out_dir <- function(repo_root = bpf_repo_root(SCRIPT_DIR)) {
  file.path(bpf_combine_fixo2_eigen_attractor_result_dir(repo_root), "pooled_embedding_curve_class", "tables")
}

default_output_file <- function(out_dir) {
  file.path(out_dir, "fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv")
}

summarize_seed_curve <- function(seed_data) {
  seed_data$O2_pct <- as_numeric_column(seed_data$O2_pct)
  seed_data$smoothed_dominant_mean_ploidy <- as_numeric_column(seed_data$smoothed_dominant_mean_ploidy)
  if ("smooth_step_epsilon" %in% names(seed_data)) {
    seed_data$smooth_step_epsilon <- as_numeric_column(seed_data$smooth_step_epsilon)
  } else {
    seed_data$smooth_step_epsilon <- NA_real_
  }
  seed_data <- seed_data[is.finite(seed_data$O2_pct) & is.finite(seed_data$smoothed_dominant_mean_ploidy), , drop = FALSE]
  seed_data <- seed_data[order(seed_data$O2_pct), , drop = FALSE]

  seed_id <- as.character(seed_data$seed_id[[1L]])
  n_o2 <- nrow(seed_data)
  if (n_o2 < 2L) {
    return(data.frame(
      seed_id = seed_id,
      seed_number = seed_number_from_id(seed_id),
      n_o2_points = n_o2,
      o2_min = NA_real_,
      o2_max = NA_real_,
      smoothed_ploidy_at_o2_min = NA_real_,
      smoothed_ploidy_at_o2_max = NA_real_,
      smoothed_ploidy_range = NA_real_,
      net_smoothed_ploidy_change = NA_real_,
      average_slope = NA_real_,
      mean_interval_slope = NA_real_,
      weighted_mean_interval_slope = NA_real_,
      median_interval_slope = NA_real_,
      min_interval_slope = NA_real_,
      max_interval_slope = NA_real_,
      positive_interval_fraction = NA_real_,
      negative_interval_fraction = NA_real_,
      near_flat_interval_fraction = NA_real_,
      n_intervals = 0L,
      smooth_step_epsilon = first_finite(seed_data$smooth_step_epsilon),
      stringsAsFactors = FALSE
    ))
  }

  x <- seed_data$O2_pct
  y <- seed_data$smoothed_dominant_mean_ploidy
  dx <- diff(x)
  dy <- diff(y)
  ok <- is.finite(dx) & is.finite(dy) & dx > 0
  dx <- dx[ok]
  dy <- dy[ok]
  slopes <- dy / dx
  step_epsilon <- first_finite(seed_data$smooth_step_epsilon, default = 0)
  n_intervals <- length(slopes)
  o2_range <- x[[length(x)]] - x[[1L]]
  net_change <- y[[length(y)]] - y[[1L]]

  data.frame(
    seed_id = seed_id,
    seed_number = seed_number_from_id(seed_id),
    n_o2_points = n_o2,
    o2_min = x[[1L]],
    o2_max = x[[length(x)]],
    smoothed_ploidy_at_o2_min = y[[1L]],
    smoothed_ploidy_at_o2_max = y[[length(y)]],
    smoothed_ploidy_range = max(y, na.rm = TRUE) - min(y, na.rm = TRUE),
    net_smoothed_ploidy_change = net_change,
    average_slope = if (is.finite(o2_range) && o2_range > 0) net_change / o2_range else NA_real_,
    mean_interval_slope = if (n_intervals) mean(slopes, na.rm = TRUE) else NA_real_,
    weighted_mean_interval_slope = if (n_intervals) stats::weighted.mean(slopes, dx, na.rm = TRUE) else NA_real_,
    median_interval_slope = if (n_intervals) stats::median(slopes, na.rm = TRUE) else NA_real_,
    min_interval_slope = if (n_intervals) min(slopes, na.rm = TRUE) else NA_real_,
    max_interval_slope = if (n_intervals) max(slopes, na.rm = TRUE) else NA_real_,
    positive_interval_fraction = safe_fraction(sum(dy > step_epsilon, na.rm = TRUE), n_intervals),
    negative_interval_fraction = safe_fraction(sum(dy < -step_epsilon, na.rm = TRUE), n_intervals),
    near_flat_interval_fraction = safe_fraction(sum(abs(dy) <= step_epsilon, na.rm = TRUE), n_intervals),
    n_intervals = n_intervals,
    smooth_step_epsilon = step_epsilon,
    stringsAsFactors = FALSE
  )
}

compute_average_slope_table <- function(curves, by_seed = data.frame()) {
  require_columns(curves, c("seed_id", "O2_pct", "smoothed_dominant_mean_ploidy"), "curve table")
  groups <- split(curves, curves$seed_id)
  out <- do.call(rbind, lapply(groups, summarize_seed_curve))
  row.names(out) <- NULL

  if (nrow(by_seed)) {
    keep <- intersect(
      c(
        "seed_id", "seed_number", "curve_class", "final_interpretation_class",
        "smooth_curve_class", "pointwise_curve_class", "class_changed",
        "classification_rule_version", "smooth_classification_rule_version",
        "objective", "objective_source", "objective_total", "objective_data",
        "objective_burden", "objective_ploidy", "convergence_status",
        "monotonicity_reliability"
      ),
      names(by_seed)
    )
    if ("seed_id" %in% keep) {
      annot <- by_seed[, keep, drop = FALSE]
      annot <- annot[!duplicated(annot$seed_id), , drop = FALSE]
      if ("seed_number" %in% names(annot)) {
        names(annot)[names(annot) == "seed_number"] <- "seed_number_annotation"
      }
      out <- merge(out, annot, by = "seed_id", all.x = TRUE, sort = FALSE)
      if ("seed_number_annotation" %in% names(out)) {
        out$seed_number <- ifelse(is.na(out$seed_number), out$seed_number_annotation, out$seed_number)
        out$seed_number_annotation <- NULL
      }
    }
  }

  out <- out[order(out$seed_number, out$seed_id), , drop = FALSE]
  row.names(out) <- NULL
  out
}

write_average_slope_analysis <- function(curve_table = default_curve_table(),
                                         by_seed_table = default_by_seed_table(),
                                         out_dir = default_out_dir(),
                                         output_file = default_output_file(out_dir),
                                         dry_run = FALSE) {
  curve_table <- normalizePath(path.expand(curve_table), mustWork = FALSE)
  by_seed_table <- normalizePath(path.expand(by_seed_table), mustWork = FALSE)
  out_dir <- normalizePath(path.expand(out_dir), mustWork = FALSE)
  output_file <- normalizePath(path.expand(output_file), mustWork = FALSE)

  message("Curve table: ", curve_table)
  message("By-seed annotation table: ", by_seed_table)
  message("Output table: ", output_file)
  if (isTRUE(dry_run)) return(invisible(output_file))

  curves <- read_tsv(curve_table)
  by_seed <- read_tsv(by_seed_table, required = FALSE)
  out <- compute_average_slope_table(curves, by_seed)
  write_tsv(out, output_file)
  message("Wrote ", nrow(out), " seed-level average-slope rows: ", output_file)
  invisible(output_file)
}

main <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  args <- bpf_parse_args(raw_args)
  repo_root <- bpf_repo_root(SCRIPT_DIR)
  out_dir <- bpf_resolve_repo_path(args$out_dir %||% default_out_dir(repo_root), repo_root)
  curve_table <- bpf_resolve_repo_path(args$curve_table %||% default_curve_table(repo_root), repo_root)
  by_seed_table <- bpf_resolve_repo_path(args$by_seed_table %||% default_by_seed_table(repo_root), repo_root)
  output_file <- bpf_resolve_repo_path(args$output_file %||% default_output_file(out_dir), repo_root)
  dry_run <- bpf_as_bool(args$dry_run, FALSE)

  write_average_slope_analysis(
    curve_table = curve_table,
    by_seed_table = by_seed_table,
    out_dir = out_dir,
    output_file = output_file,
    dry_run = dry_run
  )
}

if (identical(environment(), globalenv())) {
  main()
}
