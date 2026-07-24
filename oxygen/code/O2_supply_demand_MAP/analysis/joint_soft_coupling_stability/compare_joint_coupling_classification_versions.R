#!/usr/bin/env Rscript

# Compare a threshold-sensitivity result with its materialized source version.
# The comparison is descriptive and keeps pair-level and cross-pair units
# separate from the underlying seed-by-parameter transition counts.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

classification_text <- function(config) {
  value <- function(key, fallback = NA_character_) {
    hit <- config$value[config$key == key]
    if (length(hit)) as.character(hit[[1L]]) else fallback
  }
  threshold <- suppressWarnings(as.numeric(value("class_threshold", "1.1")))
  lower <- suppressWarnings(as.numeric(value("class_lower_bound", as.character(1 / threshold))))
  upper <- suppressWarnings(as.numeric(value("class_upper_bound", as.character(threshold))))
  boundary <- value("class_boundary_rule", "classb_inclusive")
  paste0("lower=", signif(lower, 8), ";upper=", signif(upper, 8), ";boundary=", boundary)
}

compare_joint_coupling_classification_versions <- function(source_analysis_dir, current_analysis_dir, out_dir) {
  source_analysis_dir <- normalizePath(source_analysis_dir, mustWork = TRUE)
  current_analysis_dir <- normalizePath(current_analysis_dir, mustWork = TRUE)
  if (identical(source_analysis_dir, current_analysis_dir)) stop("Source and current analysis directories must differ", call. = FALSE)
  old <- o2jca_read_tsv(file.path(source_analysis_dir, "soft_coupling_master_long.tsv"))
  new <- o2jca_read_tsv(file.path(current_analysis_dir, "soft_coupling_master_long.tsv"))
  required <- c("pair_id", "seed", "parameter", "ratio_vivo_to_vitro", "ratio_class")
  o2jca_assert_columns(old, required, "source master")
  o2jca_assert_columns(new, required, "current master")
  key_columns <- c("pair_id", "seed", "parameter")
  o2jca_assert_unique_key(old, key_columns, "source master")
  o2jca_assert_unique_key(new, key_columns, "current master")
  old$key <- do.call(paste, c(old[key_columns], sep = "\r"))
  new$key <- do.call(paste, c(new[key_columns], sep = "\r"))
  if (!setequal(old$key, new$key)) stop("Source and current master keys differ", call. = FALSE)
  old <- old[match(new$key, old$key), , drop = FALSE]
  ratio_delta <- abs(suppressWarnings(as.numeric(old$ratio_vivo_to_vitro)) - suppressWarnings(as.numeric(new$ratio_vivo_to_vitro)))
  if (any(is.finite(ratio_delta) & ratio_delta > 1e-12)) stop("Ratios changed between classification versions", call. = FALSE)

  classes <- c("ClassA", "ClassB", "ClassC", "Invalid")
  transition_table <- table(
    source_class = factor(old$ratio_class, levels = classes),
    current_class = factor(new$ratio_class, levels = classes)
  )
  transitions <- as.data.frame(transition_table, stringsAsFactors = FALSE)
  names(transitions)[[3L]] <- "n_rows"
  transitions$proportion_all_rows <- transitions$n_rows / nrow(new)
  transitions$class_changed <- transitions$source_class != transitions$current_class

  old_within <- o2jca_read_tsv(file.path(source_analysis_dir, "within_pair_parameter_stability.tsv"))
  new_within <- o2jca_read_tsv(file.path(current_analysis_dir, "within_pair_parameter_stability.tsv"))
  within_key <- paste(old_within$pair_id, old_within$parameter, sep = "\r")
  current_key <- paste(new_within$pair_id, new_within$parameter, sep = "\r")
  if (!setequal(within_key, current_key)) stop("Within-pair summary keys differ", call. = FALSE)
  old_within <- old_within[match(current_key, within_key), , drop = FALSE]
  pair_parameter <- data.frame(
    pair_id = new_within$pair_id,
    parameter = new_within$parameter,
    source_dominant_class = old_within$dominant_class,
    current_dominant_class = new_within$dominant_class,
    source_dominant_fraction = old_within$dominant_fraction,
    current_dominant_fraction = new_within$dominant_fraction,
    delta_prop_ClassA = new_within$prop_ClassA - old_within$prop_ClassA,
    delta_prop_ClassB = new_within$prop_ClassB - old_within$prop_ClassB,
    delta_prop_ClassC = new_within$prop_ClassC - old_within$prop_ClassC,
    dominant_class_changed = old_within$dominant_class != new_within$dominant_class,
    stringsAsFactors = FALSE
  )

  old_between <- o2jca_read_tsv(file.path(source_analysis_dir, "between_pair_parameter_stability.tsv"))
  new_between <- o2jca_read_tsv(file.path(current_analysis_dir, "between_pair_parameter_stability.tsv"))
  old_between <- old_between[match(new_between$parameter, old_between$parameter), , drop = FALSE]
  cross_pair <- data.frame(
    parameter = new_between$parameter,
    source_cross_pair_dominant_class = old_between$cross_pair_dominant_class,
    current_cross_pair_dominant_class = new_between$cross_pair_dominant_class,
    source_cross_pair_dominant_fraction = old_between$cross_pair_dominant_fraction,
    current_cross_pair_dominant_fraction = new_between$cross_pair_dominant_fraction,
    delta_cross_pair_dominant_fraction = new_between$cross_pair_dominant_fraction - old_between$cross_pair_dominant_fraction,
    cross_pair_dominant_class_changed = old_between$cross_pair_dominant_class != new_between$cross_pair_dominant_class,
    stringsAsFactors = FALSE
  )
  process_map <- o2jca_read_tsv(file.path(current_analysis_dir, "parameter_process_map.tsv"))
  cross_pair$primary_process <- process_map$primary_process[match(cross_pair$parameter, process_map$parameter)]

  old_class <- o2jca_read_tsv(file.path(source_analysis_dir, "between_pair_class_summary.tsv"))
  new_class <- o2jca_read_tsv(file.path(current_analysis_dir, "between_pair_class_summary.tsv"))
  class_key <- paste(old_class$parameter, old_class$ratio_class, sep = "\r")
  new_class_key <- paste(new_class$parameter, new_class$ratio_class, sep = "\r")
  old_class <- old_class[match(new_class_key, class_key), , drop = FALSE]
  class_share <- data.frame(
    parameter = new_class$parameter,
    ratio_class = new_class$ratio_class,
    source_pair_balanced_mean_proportion = old_class$pair_balanced_mean_proportion,
    current_pair_balanced_mean_proportion = new_class$pair_balanced_mean_proportion,
    delta_pair_balanced_mean_proportion = new_class$pair_balanced_mean_proportion - old_class$pair_balanced_mean_proportion,
    source_n_pairs_stable90 = old_class$n_pairs_stable90,
    current_n_pairs_stable90 = new_class$n_pairs_stable90,
    delta_n_pairs_stable90 = new_class$n_pairs_stable90 - old_class$n_pairs_stable90,
    stringsAsFactors = FALSE
  )

  source_config <- o2jca_read_tsv(file.path(source_analysis_dir, "analysis_config.tsv"))
  current_config <- o2jca_read_tsv(file.path(current_analysis_dir, "analysis_config.tsv"))
  changed_rows <- sum(old$ratio_class != new$ratio_class)
  summary <- data.frame(
    metric = c(
      "source_classification", "current_classification", "n_seed_parameter_rows",
      "n_seed_parameter_rows_changed", "proportion_seed_parameter_rows_changed",
      "n_pair_parameter_units", "n_pair_parameter_dominant_class_changed",
      "n_cross_pair_parameter_units", "n_cross_pair_dominant_class_changed"
    ),
    value = c(
      classification_text(source_config), classification_text(current_config), nrow(new),
      changed_rows, changed_rows / nrow(new), nrow(pair_parameter),
      sum(pair_parameter$dominant_class_changed), nrow(cross_pair),
      sum(cross_pair$cross_pair_dominant_class_changed)
    ),
    stringsAsFactors = FALSE
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(
    o2jca_write_tsv(summary, file.path(out_dir, "classification_version_comparison_summary.tsv")),
    o2jca_write_tsv(transitions, file.path(out_dir, "seed_parameter_class_transition_counts.tsv")),
    o2jca_write_tsv(pair_parameter, file.path(out_dir, "pair_parameter_dominant_class_changes.tsv")),
    o2jca_write_tsv(cross_pair, file.path(out_dir, "cross_pair_dominant_class_changes.tsv")),
    o2jca_write_tsv(class_share, file.path(out_dir, "parameter_class_share_changes.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "classification_version_comparison_manifest.tsv"))
  invisible(summary)
}

main <- function() {
  args <- o2jca_parse_args()
  compare_joint_coupling_classification_versions(
    source_analysis_dir = args$source_analysis_dir %||% stop("--source_analysis_dir is required", call. = FALSE),
    current_analysis_dir = args$current_analysis_dir %||% stop("--current_analysis_dir is required", call. = FALSE),
    out_dir = args$out_dir %||% stop("--out_dir is required", call. = FALSE)
  )
}

if (identical(environment(), globalenv())) main()
