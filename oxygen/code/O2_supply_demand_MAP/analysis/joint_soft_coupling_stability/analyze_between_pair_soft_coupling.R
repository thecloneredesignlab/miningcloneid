#!/usr/bin/env Rscript

# Pair-balanced cross-pair stability analysis. A pair is the inferential unit;
# pooled-seed summaries are retained only as secondary descriptive outputs.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

estimate_icc <- function(data) {
  data <- data[is.finite(data$log2_ratio_vivo_to_vitro), , drop = FALSE]
  group_sizes <- table(data$pair_id)
  k <- length(group_sizes); n <- nrow(data)
  if (k < 2L || n <= k) return(c(between_variance = NA, within_variance = NA, icc = NA))
  grand <- mean(data$log2_ratio_vivo_to_vitro)
  means <- tapply(data$log2_ratio_vivo_to_vitro, data$pair_id, mean)
  ss_between <- sum(group_sizes[names(means)] * (means - grand)^2)
  ss_within <- sum(unlist(lapply(split(data$log2_ratio_vivo_to_vitro, data$pair_id), function(x) sum((x - mean(x))^2))))
  ms_between <- ss_between / (k - 1L)
  ms_within <- ss_within / (n - k)
  n_bar <- (n - sum(group_sizes^2) / n) / (k - 1L)
  between <- max((ms_between - ms_within) / n_bar, 0)
  within <- ms_within
  c(between_variance = between, within_variance = within, icc = between / (between + within))
}

analyze_between_pair_soft_coupling <- function(master_file, within_file, out_dir, n_boot = 2000L) {
  master <- o2jca_read_joint_master(master_file)
  within <- o2jca_read_tsv(within_file)
  o2jca_assert_columns(within, c("pair_id", "parameter", "dominant_class", "prop_ClassA", "prop_ClassB", "prop_ClassC"), "within-pair stability")
  class_rows <- list(); stability_rows <- list(); loo_rows <- list()
  for (parameter in o2jca_parameter_levels()) {
    param <- within[within$parameter == parameter, , drop = FALSE]
    dominant_counts <- table(factor(param$dominant_class, levels = c("ClassA", "ClassB", "ClassC")))
    dominant <- names(dominant_counts)[which.max(dominant_counts)]
    dominant_fraction <- max(dominant_counts) / sum(dominant_counts)
    raw <- master[master$parameter == parameter, , drop = FALSE]
    icc <- estimate_icc(raw)
    pair_medians <- tapply(raw$log2_ratio_vivo_to_vitro, raw$pair_id, stats::median, na.rm = TRUE)
    stability_rows[[length(stability_rows) + 1L]] <- data.frame(
      parameter = parameter,
      n_pairs = nrow(param),
      cross_pair_dominant_class = dominant,
      cross_pair_dominant_fraction = dominant_fraction,
      all_pairs_same_dominant_class = length(unique(param$dominant_class)) == 1L,
      n_pairs_dominant_ClassA = unname(dominant_counts[["ClassA"]]),
      n_pairs_dominant_ClassB = unname(dominant_counts[["ClassB"]]),
      n_pairs_dominant_ClassC = unname(dominant_counts[["ClassC"]]),
      n_pairs_stable90_ClassA = sum(param$stable90_ClassA, na.rm = TRUE),
      n_pairs_stable90_ClassB = sum(param$stable90_ClassB, na.rm = TRUE),
      n_pairs_stable90_ClassC = sum(param$stable90_ClassC, na.rm = TRUE),
      n_pairs_strict_ClassA = sum(param$intersection_ClassA, na.rm = TRUE),
      n_pairs_strict_ClassB = sum(param$intersection_ClassB, na.rm = TRUE),
      n_pairs_strict_ClassC = sum(param$intersection_ClassC, na.rm = TRUE),
      pair_median_direction_consistency = max(mean(pair_medians < 0), mean(pair_medians > 0)),
      between_pair_variance_log2_ratio = icc[["between_variance"]],
      within_pair_variance_log2_ratio = icc[["within_variance"]],
      intraclass_correlation = icc[["icc"]],
      shared_invitro_anchor = length(unique(raw$invitro_seed)) == 1L,
      invitro_anchor_seed = paste(sort(unique(raw$invitro_seed)), collapse = ","),
      stringsAsFactors = FALSE
    )
    for (class_name in c("ClassA", "ClassB", "ClassC")) {
      values <- param[[paste0("prop_", class_name)]]
      ci <- o2jca_bootstrap_pair_mean(values, n_boot = n_boot, seed = 5826L + match(parameter, o2jca_parameter_levels()) * 10L + match(class_name, c("ClassA", "ClassB", "ClassC")))
      class_rows[[length(class_rows) + 1L]] <- data.frame(
        parameter = parameter, ratio_class = class_name, n_pairs = sum(is.finite(values)),
        pair_balanced_mean_proportion = ci[["mean"]], bootstrap_ci_lower = ci[["lower"]],
        bootstrap_ci_upper = ci[["upper"]], pair_sd = stats::sd(values, na.rm = TRUE),
        pair_min = min(values, na.rm = TRUE), pair_max = max(values, na.rm = TRUE),
        n_pairs_union = sum(param[[paste0("union_", class_name)]], na.rm = TRUE),
        n_pairs_stable80 = sum(param[[paste0("stable80_", class_name)]], na.rm = TRUE),
        n_pairs_stable90 = sum(param[[paste0("stable90_", class_name)]], na.rm = TRUE),
        n_pairs_stable95 = sum(param[[paste0("stable95_", class_name)]], na.rm = TRUE),
        n_pairs_strict_intersection = sum(param[[paste0("intersection_", class_name)]], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    for (excluded_pair in param$pair_id) {
      kept <- param[param$pair_id != excluded_pair, , drop = FALSE]
      counts <- table(factor(kept$dominant_class, levels = c("ClassA", "ClassB", "ClassC")))
      loo_rows[[length(loo_rows) + 1L]] <- data.frame(
        parameter = parameter, excluded_pair = excluded_pair,
        leave_one_out_dominant_class = names(counts)[which.max(counts)],
        leave_one_out_dominant_fraction = max(counts) / sum(counts),
        agrees_with_full_result = names(counts)[which.max(counts)] == dominant,
        stringsAsFactors = FALSE
      )
    }
  }
  class_summary <- do.call(rbind, class_rows)
  stability <- do.call(rbind, stability_rows)
  leave_one_out <- do.call(rbind, loo_rows)

  pair_ids <- sort(unique(within$pair_id))
  similarity_rows <- list()
  for (i in seq_len(length(pair_ids) - 1L)) {
    for (j in (i + 1L):length(pair_ids)) {
      x <- within[within$pair_id == pair_ids[[i]], , drop = FALSE]
      y <- within[within$pair_id == pair_ids[[j]], , drop = FALSE]
      y <- y[match(x$parameter, y$parameter), , drop = FALSE]
      js <- vapply(seq_len(nrow(x)), function(k) {
        o2jca_js_divergence(
          unlist(x[k, c("prop_ClassA", "prop_ClassB", "prop_ClassC")]),
          unlist(y[k, c("prop_ClassA", "prop_ClassB", "prop_ClassC")])
        )
      }, numeric(1L))
      similarity_rows[[length(similarity_rows) + 1L]] <- data.frame(
        pair_1 = pair_ids[[i]], pair_2 = pair_ids[[j]],
        dominant_class_agreement = mean(x$dominant_class == y$dominant_class),
        mean_jensen_shannon_divergence = mean(js, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  pair_similarity <- do.call(rbind, similarity_rows)

  pooled_rows <- lapply(o2jca_group_split(master, c("parameter")), function(group) {
    threshold <- unique(suppressWarnings(as.numeric(group$class_threshold)))[[1L]]
    lower_bound <- if ("class_lower_bound" %in% names(group)) unique(suppressWarnings(as.numeric(group$class_lower_bound)))[[1L]] else 1 / threshold
    upper_bound <- if ("class_upper_bound" %in% names(group)) unique(suppressWarnings(as.numeric(group$class_upper_bound)))[[1L]] else threshold
    data.frame(
      parameter = group$parameter[[1L]],
      o2jca_summarize_parameter_group(group, threshold, lower_bound, upper_bound)
    )
  })
  pooled <- do.call(rbind, pooled_rows)
  pooled$inference_role <- "secondary_pooled_seed_description"

  membership_rows <- list(); pattern_rows <- list()
  pair_levels <- sort(unique(within$pair_id))
  criteria <- c(union = "union", stable80 = "stable80", stable90 = "stable90", stable95 = "stable95", strict = "intersection")
  for (class_name in c("ClassA", "ClassB", "ClassC")) {
    for (criterion in names(criteria)) {
      column <- paste0(criteria[[criterion]], "_", class_name)
      for (parameter in o2jca_parameter_levels()) {
        param <- within[within$parameter == parameter, , drop = FALSE]
        included <- as.logical(param[[column]])
        included[is.na(included)] <- FALSE
        selected_pairs <- sort(param$pair_id[included])
        pattern_id <- if (length(selected_pairs)) paste(selected_pairs, collapse = "|") else "none"
        pattern_rows[[length(pattern_rows) + 1L]] <- data.frame(
          ratio_class = class_name, criterion = criterion, parameter = parameter,
          n_pairs = length(selected_pairs), pair_pattern = pattern_id,
          pairs = paste(vapply(selected_pairs, o2jca_pair_short_label, character(1L)), collapse = "; "),
          stringsAsFactors = FALSE
        )
        membership_rows[[length(membership_rows) + 1L]] <- data.frame(
          ratio_class = class_name, criterion = criterion, parameter = parameter,
          pair_id = pair_levels, included = pair_levels %in% selected_pairs,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  intersection_patterns <- do.call(rbind, pattern_rows)
  set_membership <- do.call(rbind, membership_rows)

  files <- c(
    o2jca_write_tsv(class_summary, file.path(out_dir, "between_pair_class_summary.tsv")),
    o2jca_write_tsv(stability, file.path(out_dir, "between_pair_parameter_stability.tsv")),
    o2jca_write_tsv(leave_one_out, file.path(out_dir, "between_pair_leave_one_pair_out.tsv")),
    o2jca_write_tsv(pair_similarity, file.path(out_dir, "between_pair_pair_similarity.tsv")),
    o2jca_write_tsv(pooled, file.path(out_dir, "between_pair_pooled_seed_secondary.tsv")),
    o2jca_write_tsv(intersection_patterns, file.path(out_dir, "class_intersection_pattern_summary.tsv")),
    o2jca_write_tsv(set_membership, file.path(out_dir, "class_pair_set_membership.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "between_pair_analysis_manifest.tsv"))
  invisible(stability)
}

main <- function() {
  args <- o2jca_parse_args(); out_dir <- args$out_dir %||% stop("--out_dir is required", call. = FALSE)
  analyze_between_pair_soft_coupling(
    master_file = args$master_file %||% file.path(out_dir, "soft_coupling_master_long.tsv"),
    within_file = args$within_file %||% file.path(out_dir, "within_pair_parameter_stability.tsv"),
    out_dir = out_dir,
    n_boot = o2jca_as_int(args$n_boot, 2000L)
  )
}

if (identical(environment(), globalenv())) main()
