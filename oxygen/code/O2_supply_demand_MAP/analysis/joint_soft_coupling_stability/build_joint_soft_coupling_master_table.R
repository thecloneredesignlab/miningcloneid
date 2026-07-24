#!/usr/bin/env Rscript

# Materialize and validate the pair x seed x parameter table used by all joint
# soft-coupling analyses. This script reads completed fit summaries only.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

build_joint_soft_coupling_master_table <- function(
    result_root,
    out_dir,
    class_threshold = 1.1,
    class_lower_bound = NA_real_,
    class_upper_bound = NA_real_,
    class_boundary_rule = "classb_inclusive",
    pair_pattern = "^fit_joint_.*_vt_seed[0-9]+$",
    max_pairs = NA_integer_,
    max_seeds = NA_integer_) {
  class_spec <- o2jca_classification_spec(
    class_threshold, class_lower_bound, class_upper_bound, class_boundary_rule
  )
  pair_dirs <- o2jca_discover_pair_dirs(result_root, pair_pattern, max_pairs)
  all_rows <- list()
  pair_quality <- list()
  pair_metadata <- list()

  soft_keep <- c(
    "parameter", "vivo_natural", "vitro_natural", "ratio_vivo_to_vitro",
    "penalty_paid", "penalty_region", "abs_delta_over_sigma",
    "welsch_saturation_fraction", "feasible_at_solution",
    "feasible_before_projection", "projection_applied", "projection_action", "seed"
  )
  summary_keep <- c(
    "seed", "objective", "objective_invivo", "objective_invitro",
    "objective_soft_coupling", "objective_constraints", "joint_constraints_pass",
    "optimizer_interrupted", "n_at_bound_active", "n_near_bound_only_active",
    "boundary_penalty_active", "min_rel_dist_active", "objective_rank"
  )

  for (pair_dir in pair_dirs) {
    pair_id <- basename(pair_dir)
    message("Reading ", pair_id)
    soft_path <- file.path(pair_dir, "extra_results", "joint_soft_coupling_all.tsv")
    summary_path <- file.path(pair_dir, "extra_results", "seed_summary.tsv")
    soft <- o2jca_read_tsv(soft_path)
    summary <- o2jca_read_tsv(summary_path)
    o2jca_assert_columns(soft, c("parameter", "vivo_natural", "vitro_natural", "ratio_vivo_to_vitro", "seed"), pair_id)
    o2jca_assert_columns(summary, c("seed", "objective"), paste0(pair_id, " seed summary"))
    soft <- soft[, intersect(soft_keep, names(soft)), drop = FALSE]
    summary <- summary[, intersect(summary_keep, names(summary)), drop = FALSE]
    soft$seed <- o2jca_normalize_seed(soft$seed)
    summary$seed <- o2jca_normalize_seed(summary$seed)
    if (is.finite(max_seeds)) {
      keep_seeds <- paste0("seed", seq_len(max_seeds))
      soft <- soft[soft$seed %in% keep_seeds, , drop = FALSE]
      summary <- summary[summary$seed %in% keep_seeds, , drop = FALSE]
    }
    summary <- summary[!duplicated(summary$seed), , drop = FALSE]
    idx <- match(soft$seed, summary$seed)
    for (column in setdiff(names(summary), "seed")) soft[[column]] <- summary[[column]][idx]
    soft$pair_id <- pair_id
    soft$seed_number <- o2jca_seed_number(soft$seed)
    soft$vivo_natural <- suppressWarnings(as.numeric(soft$vivo_natural))
    soft$vitro_natural <- suppressWarnings(as.numeric(soft$vitro_natural))
    soft$ratio_reported <- suppressWarnings(as.numeric(soft$ratio_vivo_to_vitro))
    soft$ratio_vivo_to_vitro <- soft$vivo_natural / soft$vitro_natural
    soft$ratio_relative_error <- abs(soft$ratio_reported - soft$ratio_vivo_to_vitro) /
      pmax(abs(soft$ratio_vivo_to_vitro), .Machine$double.eps)
    soft$log2_ratio_vivo_to_vitro <- ifelse(
      is.finite(soft$ratio_vivo_to_vitro) & soft$ratio_vivo_to_vitro > 0,
      log2(soft$ratio_vivo_to_vitro),
      NA_real_
    )
    soft$ratio_class <- as.character(o2jca_ratio_class(
      soft$ratio_vivo_to_vitro,
      threshold = class_spec$threshold,
      lower_bound = class_spec$lower_bound,
      upper_bound = class_spec$upper_bound,
      boundary_rule = class_spec$boundary_rule
    ))
    soft$class_threshold <- class_spec$threshold
    soft$class_lower_bound <- class_spec$lower_bound
    soft$class_upper_bound <- class_spec$upper_bound
    soft$class_boundary_rule <- class_spec$boundary_rule
    soft$class_scheme <- class_spec$scheme
    soft$class_rule <- class_spec$class_rule
    meta <- o2jca_pair_metadata(pair_id)
    for (column in setdiff(names(meta), "pair_id")) soft[[column]] <- meta[[column]][[1L]]
    all_rows[[length(all_rows) + 1L]] <- soft
    pair_metadata[[length(pair_metadata) + 1L]] <- meta

    valid_ratio <- is.finite(soft$ratio_vivo_to_vitro) & soft$ratio_vivo_to_vitro > 0
    expected <- length(unique(soft$seed)) * length(o2jca_parameter_levels())
    pair_quality[[length(pair_quality) + 1L]] <- data.frame(
      pair_id = pair_id,
      n_seed = length(unique(soft$seed)),
      n_parameter = length(unique(soft$parameter)),
      n_rows = nrow(soft),
      expected_rows = expected,
      row_count_pass = nrow(soft) == expected,
      parameter_set_pass = setequal(unique(soft$parameter), o2jca_parameter_levels()),
      n_invalid_ratio = sum(!valid_ratio),
      max_ratio_relative_error = if (any(valid_ratio)) max(soft$ratio_relative_error[valid_ratio], na.rm = TRUE) else NA_real_,
      ratio_recalculation_pass = all(!valid_ratio | soft$ratio_relative_error <= 1e-10, na.rm = TRUE),
      n_missing_objective = sum(!is.finite(suppressWarnings(as.numeric(soft$objective)))),
      stringsAsFactors = FALSE
    )
  }

  master <- do.call(rbind, all_rows)
  master <- master[order(master$pair_id, master$seed_number, match(master$parameter, o2jca_parameter_levels())), , drop = FALSE]
  o2jca_assert_unique_key(master, c("pair_id", "seed", "parameter"), "joint soft-coupling master table")
  quality <- do.call(rbind, pair_quality)
  if (!all(quality$row_count_pass) || !all(quality$parameter_set_pass) || !all(quality$ratio_recalculation_pass)) {
    stop("Input quality contract failed; inspect input_quality_summary.tsv after resolving the source data", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(
    o2jca_write_tsv(master, file.path(out_dir, "soft_coupling_master_long.tsv")),
    o2jca_write_tsv(master[master$ratio_class == "Invalid", , drop = FALSE], file.path(out_dir, "excluded_ratio_rows.tsv")),
    o2jca_write_tsv(quality, file.path(out_dir, "input_quality_summary.tsv")),
    o2jca_write_tsv(do.call(rbind, pair_metadata), file.path(out_dir, "pair_metadata.tsv")),
    o2jca_write_tsv(o2jca_process_map(), file.path(out_dir, "parameter_process_map.tsv"))
  )
  config <- data.frame(
    key = c(
      "result_root", "analysis_label", "class_threshold", "class_lower_bound",
      "class_upper_bound", "class_boundary_rule", "class_scheme", "class_rule",
      "max_pairs", "max_seeds", "n_pairs", "n_rows"
    ),
    value = c(
      normalizePath(result_root, mustWork = TRUE), basename(dirname(dirname(normalizePath(out_dir, mustWork = TRUE)))),
      class_spec$threshold, class_spec$lower_bound, class_spec$upper_bound,
      class_spec$boundary_rule, class_spec$scheme, class_spec$class_rule,
      max_pairs, max_seeds, length(pair_dirs), nrow(master)
    ),
    stringsAsFactors = FALSE
  )
  files <- c(files, o2jca_write_tsv(config, file.path(out_dir, "analysis_config.tsv")))
  o2jca_write_manifest(files, file.path(out_dir, "master_table_manifest.tsv"))
  invisible(master)
}

main <- function() {
  args <- o2jca_parse_args()
  result_root <- args$result_root %||% stop("--result_root is required", call. = FALSE)
  output_root <- args$output_root %||% o2jca_default_output_root(result_root, WORKFLOW_ROOT)
  out_dir <- args$out_dir %||% file.path(output_root, "tables", "soft_coupling")
  o2jca_assert_separate_output_root(result_root, out_dir)
  build_joint_soft_coupling_master_table(
    result_root = result_root,
    out_dir = out_dir,
    class_threshold = o2jca_as_num(args$class_threshold, 1.1),
    class_lower_bound = o2jca_as_num(args$class_lower_bound, NA_real_),
    class_upper_bound = o2jca_as_num(args$class_upper_bound, NA_real_),
    class_boundary_rule = args$class_boundary_rule %||% "classb_inclusive",
    pair_pattern = args$pair_pattern %||% "^fit_joint_.*_vt_seed[0-9]+$",
    max_pairs = o2jca_as_int(args$max_pairs, NA_integer_),
    max_seeds = o2jca_as_int(args$max_seeds, NA_integer_)
  )
}

if (identical(environment(), globalenv())) main()
