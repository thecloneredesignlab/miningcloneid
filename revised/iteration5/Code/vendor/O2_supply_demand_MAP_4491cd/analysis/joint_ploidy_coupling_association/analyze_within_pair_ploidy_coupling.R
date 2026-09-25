#!/usr/bin/env Rscript

# Within-pair comparison of fitted in-vivo values, in-vitro values, and their
# coupling ratios across CatA/B/C ploidy trajectories.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

numeric_summary <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  c(
    n = length(x), mean = if (length(x)) mean(x) else NA_real_,
    sd = if (length(x) > 1L) stats::sd(x) else if (length(x)) 0 else NA_real_,
    median = if (length(x)) stats::median(x) else NA_real_,
    q25 = if (length(x)) o2jca_quantile(x, 0.25) else NA_real_,
    q75 = if (length(x)) o2jca_quantile(x, 0.75) else NA_real_
  )
}

kruskal_effect <- function(value, category) {
  keep <- is.finite(value) & category %in% c("CatA", "CatB", "CatC")
  value <- value[keep]; category <- droplevels(factor(category[keep]))
  if (length(value) < 3L || nlevels(category) < 2L) return(c(statistic = NA, df = NA, p_value = NA, epsilon_squared = NA))
  test <- stats::kruskal.test(value ~ category)
  h <- unname(test$statistic); k <- nlevels(category); n <- length(value)
  c(statistic = h, df = unname(test$parameter), p_value = test$p.value, epsilon_squared = max((h - k + 1) / (n - k), 0))
}

analyze_within_pair_ploidy_coupling <- function(master_file, category_file, out_dir) {
  master <- o2jca_read_joint_master(master_file)
  category <- o2jca_read_tsv(category_file)
  o2jca_assert_columns(category, c("pair_id", "seed", "ploidy_category"), "ploidy seed categories")
  data <- merge(master, category[, c("pair_id", "seed", "ploidy_category", "category_reason")], by = c("pair_id", "seed"), all.x = TRUE, sort = FALSE)
  data$natural_delta_vivo_minus_vitro <- data$vivo_natural - data$vitro_natural
  data$analysis_included <- data$ploidy_category %in% c("CatA", "CatB", "CatC")
  included <- data[data$analysis_included, , drop = FALSE]

  metrics <- c("vivo_natural", "vitro_natural", "ratio_vivo_to_vitro", "log2_ratio_vivo_to_vitro", "natural_delta_vivo_minus_vitro", "penalty_paid", "objective")
  groups <- o2jca_group_split(included, c("pair_id", "ploidy_category", "parameter"))
  summary_rows <- list()
  for (group in groups) {
    for (metric in metrics) {
      if (!(metric %in% names(group))) next
      stats <- numeric_summary(group[[metric]])
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        pair_id = group$pair_id[[1L]], ploidy_category = group$ploidy_category[[1L]],
        parameter = group$parameter[[1L]], metric = metric,
        n = stats[["n"]], mean = stats[["mean"]], sd = stats[["sd"]], median = stats[["median"]],
        q25 = stats[["q25"]], q75 = stats[["q75"]], stringsAsFactors = FALSE
      )
    }
  }
  parameter_summary <- do.call(rbind, summary_rows)

  class_counts <- as.data.frame(table(
    included$pair_id, included$ploidy_category, included$parameter,
    factor(included$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
  ), stringsAsFactors = FALSE)
  names(class_counts) <- c("pair_id", "ploidy_category", "parameter", "ratio_class", "n")
  class_counts$proportion_within_category_parameter <- ave(class_counts$n, class_counts$pair_id, class_counts$ploidy_category, class_counts$parameter, FUN = function(x) if (sum(x)) x / sum(x) else 0)

  effect_rows <- list()
  for (group in o2jca_group_split(included, c("pair_id", "parameter"))) {
    n_observed_categories <- length(unique(group$ploidy_category))
    for (metric in metrics) {
      if (!(metric %in% names(group))) next
      effect <- kruskal_effect(suppressWarnings(as.numeric(group[[metric]])), group$ploidy_category)
      effect_rows[[length(effect_rows) + 1L]] <- data.frame(
        pair_id = group$pair_id[[1L]], parameter = group$parameter[[1L]], metric = metric,
        statistic = effect[["statistic"]], df = effect[["df"]], p_value = effect[["p_value"]],
        epsilon_squared = effect[["epsilon_squared"]],
        n_observed_categories = n_observed_categories,
        within_pair_category_comparison_estimable = n_observed_categories >= 2L,
        stringsAsFactors = FALSE
      )
    }
  }
  effects <- do.call(rbind, effect_rows)
  effects$fdr_within_pair_metric <- ave(effects$p_value, effects$pair_id, effects$metric, FUN = function(x) stats::p.adjust(x, method = "BH"))

  files <- c(
    o2jca_write_tsv(data, file.path(out_dir, "ploidy_coupling_master_long.tsv")),
    o2jca_write_tsv(parameter_summary, file.path(out_dir, "within_pair_cat_parameter_summary.tsv")),
    o2jca_write_tsv(class_counts, file.path(out_dir, "within_pair_cat_class_summary.tsv")),
    o2jca_write_tsv(effects, file.path(out_dir, "within_pair_cat_parameter_effects.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "within_pair_ploidy_coupling_manifest.tsv"))
  invisible(data)
}

main <- function() {
  args <- o2jca_parse_args(); out_dir <- args$out_dir %||% stop("--out_dir is required", call. = FALSE)
  ratio_dir <- args$ratio_analysis_dir %||% stop("--ratio_analysis_dir is required", call. = FALSE)
  analyze_within_pair_ploidy_coupling(
    master_file = args$master_file %||% file.path(ratio_dir, "soft_coupling_master_long.tsv"),
    category_file = args$category_file %||% file.path(out_dir, "ploidy_seed_categories.tsv"),
    out_dir = out_dir
  )
}

if (identical(environment(), globalenv())) main()
