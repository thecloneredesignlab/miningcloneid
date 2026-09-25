#!/usr/bin/env Rscript

# Pair-balanced synthesis of CatA/B/C parameter associations across joint-fit
# warm-up pairs.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

analyze_between_pair_ploidy_coupling <- function(out_dir, n_boot = 2000L) {
  summary <- o2jca_read_tsv(file.path(out_dir, "within_pair_cat_parameter_summary.tsv"))
  effects <- o2jca_read_tsv(file.path(out_dir, "within_pair_cat_parameter_effects.tsv"))
  groups <- o2jca_group_split(summary, c("ploidy_category", "parameter", "metric"))
  between_rows <- lapply(groups, function(group) {
    ci <- o2jca_bootstrap_pair_mean(group$mean, n_boot, seed = 5826L + match(group$parameter[[1L]], o2jca_parameter_levels()))
    signs <- sign(group$median)
    data.frame(
      ploidy_category = group$ploidy_category[[1L]], parameter = group$parameter[[1L]], metric = group$metric[[1L]],
      n_pairs = sum(is.finite(group$mean)), pair_balanced_mean = ci[["mean"]],
      bootstrap_ci_lower = ci[["lower"]], bootstrap_ci_upper = ci[["upper"]],
      between_pair_sd = stats::sd(group$mean, na.rm = TRUE),
      pair_median_direction_consistency = if (any(signs != 0, na.rm = TRUE)) max(mean(signs < 0, na.rm = TRUE), mean(signs > 0, na.rm = TRUE)) else 1,
      stringsAsFactors = FALSE
    )
  })
  between_summary <- do.call(rbind, between_rows)

  effect_groups <- o2jca_group_split(effects, c("parameter", "metric"))
  effect_summary <- do.call(rbind, lapply(effect_groups, function(group) {
    data.frame(
      parameter = group$parameter[[1L]], metric = group$metric[[1L]], n_pairs = sum(is.finite(group$epsilon_squared)),
      mean_epsilon_squared = mean(group$epsilon_squared, na.rm = TRUE),
      median_epsilon_squared = stats::median(group$epsilon_squared, na.rm = TRUE),
      n_pairs_fdr_lt_0p05 = sum(group$fdr_within_pair_metric < 0.05, na.rm = TRUE),
      n_pairs_same_effect_direction = NA_integer_,
      stringsAsFactors = FALSE
    )
  }))

  loo_rows <- list()
  for (group in groups) {
    for (excluded_pair in group$pair_id) {
      kept <- group[group$pair_id != excluded_pair, , drop = FALSE]
      loo_rows[[length(loo_rows) + 1L]] <- data.frame(
        ploidy_category = group$ploidy_category[[1L]], parameter = group$parameter[[1L]], metric = group$metric[[1L]],
        excluded_pair = excluded_pair, leave_one_out_pair_mean = mean(kept$mean, na.rm = TRUE),
        leave_one_out_pair_sd = stats::sd(kept$mean, na.rm = TRUE), stringsAsFactors = FALSE
      )
    }
  }
  leave_one_out <- do.call(rbind, loo_rows)
  files <- c(
    o2jca_write_tsv(between_summary, file.path(out_dir, "between_pair_cat_parameter_summary.tsv")),
    o2jca_write_tsv(effect_summary, file.path(out_dir, "between_pair_cat_parameter_effects.tsv")),
    o2jca_write_tsv(leave_one_out, file.path(out_dir, "between_pair_cat_leave_one_pair_out.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "between_pair_ploidy_coupling_manifest.tsv"))
  invisible(between_summary)
}

main <- function() {
  args <- o2jca_parse_args(); out_dir <- args$out_dir %||% stop("--out_dir is required", call. = FALSE)
  analyze_between_pair_ploidy_coupling(out_dir, o2jca_as_int(args$n_boot, 2000L))
}

if (identical(environment(), globalenv())) main()
