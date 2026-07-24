#!/usr/bin/env Rscript

# Within-pair ClassA/B/C stability analysis over independently fitted seeds.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

summarize_groups <- function(data, group_columns, lower_bound, upper_bound) {
  groups <- o2jca_group_split(data, group_columns)
  rows <- lapply(groups, function(group) {
    data.frame(
      group[1L, group_columns, drop = FALSE],
      o2jca_summarize_parameter_group(
        group, threshold = upper_bound, lower_bound = lower_bound, upper_bound = upper_bound
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

seed_similarity_summary <- function(data) {
  rows <- list()
  for (pair_id in unique(data$pair_id)) {
    pair <- data[data$pair_id == pair_id & data$ratio_class %in% c("ClassA", "ClassB", "ClassC"), , drop = FALSE]
    seeds <- sort(unique(pair$seed_number))
    parameters <- o2jca_parameter_levels()
    class_matrix <- matrix(NA_character_, length(seeds), length(parameters), dimnames = list(seeds, parameters))
    class_matrix[cbind(match(pair$seed_number, seeds), match(pair$parameter, parameters))] <- pair$ratio_class
    exact <- matrix(0, length(seeds), length(seeds))
    for (j in seq_len(ncol(class_matrix))) exact <- exact + outer(class_matrix[, j], class_matrix[, j], "==")
    values <- exact[upper.tri(exact)] / ncol(class_matrix)
    rows[[length(rows) + 1L]] <- data.frame(
      pair_id = pair_id, similarity_metric = "fraction_parameters_same_class",
      n_seed_pairs = length(values), mean = mean(values), sd = stats::sd(values),
      q05 = o2jca_quantile(values, 0.05), median = stats::median(values), q95 = o2jca_quantile(values, 0.95),
      stringsAsFactors = FALSE
    )
    for (class_name in c("ClassA", "ClassB", "ClassC")) {
      binary <- class_matrix == class_name
      intersection <- tcrossprod(binary)
      sizes <- rowSums(binary)
      union <- outer(sizes, sizes, "+") - intersection
      jaccard <- ifelse(union > 0, intersection / union, 1)
      values <- jaccard[upper.tri(jaccard)]
      rows[[length(rows) + 1L]] <- data.frame(
        pair_id = pair_id, similarity_metric = paste0(class_name, "_parameter_set_jaccard"),
        n_seed_pairs = length(values), mean = mean(values), sd = stats::sd(values),
        q05 = o2jca_quantile(values, 0.05), median = stats::median(values), q95 = o2jca_quantile(values, 0.95),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

flag_rate <- function(x) {
  if (is.null(x)) return(NA_real_)
  x <- tolower(trimws(as.character(x)))
  value <- ifelse(x %in% c("true", "t", "1", "yes"), 1, ifelse(x %in% c("false", "f", "0", "no"), 0, NA_real_))
  if (all(!is.finite(value))) NA_real_ else mean(value, na.rm = TRUE)
}

classb_direction_summary <- function(data) {
  rows <- lapply(o2jca_group_split(data, c("pair_id", "parameter")), function(group) {
    valid <- group[group$ratio_class %in% c("ClassA", "ClassB", "ClassC"), , drop = FALSE]
    classb <- valid[valid$ratio_class == "ClassB", , drop = FALSE]
    ratio <- suppressWarnings(as.numeric(classb$ratio_vivo_to_vitro))
    data.frame(
      pair_id = group$pair_id[[1L]], parameter = group$parameter[[1L]],
      n_valid = nrow(valid), n_ClassB = nrow(classb),
      classB_share_all_valid = if (nrow(valid)) nrow(classb) / nrow(valid) else NA_real_,
      ClassB_prop_ratio_lt1 = if (length(ratio)) mean(ratio < 1) else NA_real_,
      ClassB_prop_ratio_eq1 = if (length(ratio)) mean(abs(log2(ratio)) <= 1e-12) else NA_real_,
      ClassB_prop_ratio_gt1 = if (length(ratio)) mean(ratio > 1) else NA_real_,
      ClassB_ratio_mean = o2jca_safe_stat(ratio, mean),
      ClassB_ratio_variance = o2jca_safe_stat(ratio, function(z) if (length(z) > 1L) stats::var(z) else 0),
      ClassB_log2_ratio_mean = o2jca_safe_stat(log2(ratio), mean),
      ClassB_log2_ratio_sd = o2jca_safe_stat(log2(ratio), function(z) if (length(z) > 1L) stats::sd(z) else 0),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

boundary_projection_summary <- function(data) {
  rows <- lapply(o2jca_group_split(data, c("pair_id", "parameter")), function(group) {
    data.frame(
      pair_id = group$pair_id[[1L]], parameter = group$parameter[[1L]], n_seed = nrow(group),
      projection_applied_rate = flag_rate(group$projection_applied),
      feasible_before_projection_rate = flag_rate(group$feasible_before_projection),
      feasible_at_solution_rate = flag_rate(group$feasible_at_solution),
      any_bound_active_rate = if ("n_at_bound_active" %in% names(group)) mean(suppressWarnings(as.numeric(group$n_at_bound_active)) > 0, na.rm = TRUE) else NA_real_,
      mean_penalty_paid = if ("penalty_paid" %in% names(group)) o2jca_safe_stat(group$penalty_paid, mean) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

analyze_within_pair_soft_coupling <- function(master_file, out_dir, thresholds = c(1.05, 1.075, 1.1, 1.15, 1.2)) {
  data <- o2jca_read_joint_master(master_file)
  legacy_threshold <- unique(suppressWarnings(as.numeric(data$class_threshold)))
  legacy_threshold <- legacy_threshold[is.finite(legacy_threshold)][[1L]]
  primary_lower <- if ("class_lower_bound" %in% names(data)) unique(suppressWarnings(as.numeric(data$class_lower_bound))) else 1 / legacy_threshold
  primary_upper <- if ("class_upper_bound" %in% names(data)) unique(suppressWarnings(as.numeric(data$class_upper_bound))) else legacy_threshold
  primary_lower <- primary_lower[is.finite(primary_lower)][[1L]]
  primary_upper <- primary_upper[is.finite(primary_upper)][[1L]]
  primary_boundary_rule <- if ("class_boundary_rule" %in% names(data)) unique(as.character(data$class_boundary_rule))[[1L]] else "classb_inclusive"
  primary_spec <- o2jca_classification_spec(
    primary_upper, primary_lower, primary_upper, primary_boundary_rule
  )
  primary_lower <- primary_spec$lower_bound
  primary_upper <- primary_spec$upper_bound
  primary_boundary_rule <- primary_spec$boundary_rule
  summary <- summarize_groups(data, c("pair_id", "parameter"), primary_lower, primary_upper)
  class_long <- o2jca_class_long(summary, c("pair_id", "parameter"))

  threshold_rows <- list()
  for (threshold in sort(unique(thresholds))) {
    spec <- o2jca_classification_spec(threshold)
    temp <- data
    temp$ratio_class <- as.character(o2jca_ratio_class(temp$ratio_vivo_to_vitro, threshold))
    part <- summarize_groups(temp, c("pair_id", "parameter"), spec$lower_bound, spec$upper_bound)
    part <- o2jca_class_long(part, c("pair_id", "parameter"))
    part$class_threshold <- threshold
    part$class_lower_bound <- spec$lower_bound
    part$class_upper_bound <- spec$upper_bound
    part$class_boundary_rule <- spec$boundary_rule
    part$classification_label <- sprintf("symmetric %.3g", threshold)
    part$is_primary <- isTRUE(all.equal(spec$lower_bound, primary_lower, tolerance = 1e-12)) &&
      isTRUE(all.equal(spec$upper_bound, primary_upper, tolerance = 1e-12)) &&
      identical(spec$boundary_rule, primary_boundary_rule)
    threshold_rows[[length(threshold_rows) + 1L]] <- part
  }
  primary_already_present <- any(vapply(threshold_rows, function(x) any(x$is_primary), logical(1L)))
  if (!primary_already_present) {
    primary <- o2jca_class_long(summary, c("pair_id", "parameter"))
    primary$class_threshold <- primary_upper
    primary$class_lower_bound <- primary_lower
    primary$class_upper_bound <- primary_upper
    primary$class_boundary_rule <- primary_boundary_rule
    primary$classification_label <- sprintf("primary %.3g–%.3g", primary_lower, primary_upper)
    primary$is_primary <- TRUE
    threshold_rows[[length(threshold_rows) + 1L]] <- primary
  }
  threshold_sensitivity <- do.call(rbind, threshold_rows)

  stratum_rows <- list()
  for (pair_id in unique(data$pair_id)) {
    pair <- data[data$pair_id == pair_id, , drop = FALSE]
    seed_summary <- unique(pair[, c("seed", "objective"), drop = FALSE])
    strata <- o2jca_objective_strata(seed_summary)
    for (stratum in unique(strata$objective_stratum)) {
      chosen <- strata$seed[strata$objective_stratum == stratum]
      subset <- pair[pair$seed %in% chosen, , drop = FALSE]
      part <- summarize_groups(subset, c("pair_id", "parameter"), primary_lower, primary_upper)
      part <- o2jca_class_long(part, c("pair_id", "parameter"))
      part$objective_stratum <- stratum
      part$n_selected_seeds <- length(unique(chosen))
      stratum_rows[[length(stratum_rows) + 1L]] <- part
    }
  }
  objective_sensitivity <- do.call(rbind, stratum_rows)
  similarity <- seed_similarity_summary(data)
  classb_direction <- classb_direction_summary(data)
  boundary_projection <- boundary_projection_summary(data)

  class_status <- class_long
  class_status$stability_level <- ifelse(
    class_status$strict_intersection_all_seeds, "Strict 100%",
    ifelse(class_status$stable95, "Stable >=95%",
      ifelse(class_status$stable90, "Stable >=90%",
        ifelse(class_status$stable80, "Stable >=80%", ifelse(class_status$union_any_seed, "Union only", "Absent"))
      )
    )
  )

  threshold_pair_balanced <- do.call(rbind, lapply(o2jca_group_split(
    threshold_sensitivity,
    c("parameter", "ratio_class", "classification_label")
  ), function(group) {
    data.frame(
      parameter = group$parameter[[1L]], ratio_class = group$ratio_class[[1L]],
      classification_label = group$classification_label[[1L]],
      class_threshold = group$class_threshold[[1L]],
      class_lower_bound = group$class_lower_bound[[1L]],
      class_upper_bound = group$class_upper_bound[[1L]],
      class_boundary_rule = group$class_boundary_rule[[1L]],
      is_primary = group$is_primary[[1L]],
      n_pairs = nrow(group), pair_balanced_mean_proportion = mean(group$proportion, na.rm = TRUE),
      pair_min = min(group$proportion, na.rm = TRUE), pair_max = max(group$proportion, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  objective_pair_balanced <- do.call(rbind, lapply(o2jca_group_split(objective_sensitivity, c("parameter", "ratio_class", "objective_stratum")), function(group) {
    data.frame(
      parameter = group$parameter[[1L]], ratio_class = group$ratio_class[[1L]], objective_stratum = group$objective_stratum[[1L]],
      n_pairs = nrow(group), pair_balanced_mean_proportion = mean(group$proportion, na.rm = TRUE),
      pair_min = min(group$proportion, na.rm = TRUE), pair_max = max(group$proportion, na.rm = TRUE),
      total_selected_seeds = sum(group$n_selected_seeds, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  files <- c(
    o2jca_write_tsv(summary, file.path(out_dir, "within_pair_parameter_stability.tsv")),
    o2jca_write_tsv(class_long, file.path(out_dir, "within_pair_class_summary.tsv")),
    o2jca_write_tsv(threshold_sensitivity, file.path(out_dir, "within_pair_threshold_sensitivity.tsv")),
    o2jca_write_tsv(objective_sensitivity, file.path(out_dir, "within_pair_objective_sensitivity.tsv")),
    o2jca_write_tsv(similarity, file.path(out_dir, "within_pair_seed_similarity_summary.tsv")),
    o2jca_write_tsv(classb_direction, file.path(out_dir, "within_pair_classB_direction_summary.tsv")),
    o2jca_write_tsv(class_status, file.path(out_dir, "within_pair_class_membership.tsv")),
    o2jca_write_tsv(boundary_projection, file.path(out_dir, "within_pair_boundary_projection_summary.tsv")),
    o2jca_write_tsv(threshold_pair_balanced, file.path(out_dir, "threshold_sensitivity_pair_balanced.tsv")),
    o2jca_write_tsv(objective_pair_balanced, file.path(out_dir, "objective_sensitivity_pair_balanced.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "within_pair_analysis_manifest.tsv"))
  invisible(summary)
}

main <- function() {
  args <- o2jca_parse_args()
  out_dir <- args$out_dir %||% stop("--out_dir is required", call. = FALSE)
  master_file <- args$master_file %||% file.path(out_dir, "soft_coupling_master_long.tsv")
  thresholds <- suppressWarnings(as.numeric(o2jca_split_csv(args$thresholds, c("1.05", "1.075", "1.1", "1.15", "1.2"))))
  analyze_within_pair_soft_coupling(master_file, out_dir, thresholds)
}

if (identical(environment(), globalenv())) main()
