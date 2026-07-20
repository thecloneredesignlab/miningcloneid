#!/usr/bin/env Rscript

# CatA/B/C x ClassA/B/C association analysis with pair-stratified permutation,
# effect sizes, multiple-testing correction, and metric comparison.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))

contingency_details <- function(data) {
  cats <- c("CatA", "CatB", "CatC"); classes <- c("ClassA", "ClassB", "ClassC")
  tab <- table(factor(data$ploidy_category, levels = cats), factor(data$ratio_class, levels = classes))
  n <- sum(tab); row_total <- rowSums(tab); col_total <- colSums(tab)
  expected <- if (n) outer(row_total, col_total) / n else matrix(0, 3L, 3L)
  denominator <- sqrt(expected * outer(1 - row_total / max(n, 1), 1 - col_total / max(n, 1)))
  residual <- ifelse(denominator > 0, (tab - expected) / denominator, NA_real_)
  out <- expand.grid(ploidy_category = cats, ratio_class = classes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out$observed <- as.numeric(tab[cbind(match(out$ploidy_category, cats), match(out$ratio_class, classes))])
  out$expected <- as.numeric(expected[cbind(match(out$ploidy_category, cats), match(out$ratio_class, classes))])
  out$standardized_residual <- as.numeric(residual[cbind(match(out$ploidy_category, cats), match(out$ratio_class, classes))])
  out$row_percentage <- ifelse(row_total[match(out$ploidy_category, cats)] > 0, out$observed / row_total[match(out$ploidy_category, cats)], NA_real_)
  out$column_percentage <- ifelse(col_total[match(out$ratio_class, classes)] > 0, out$observed / col_total[match(out$ratio_class, classes)], NA_real_)
  list(table = tab, details = out)
}

stratified_permutation_p <- function(data, observed, permutations, seed) {
  if (!is.finite(observed) || permutations <= 0L) return(NA_real_)
  set.seed(seed)
  pair_groups <- split(seq_len(nrow(data)), data$pair_id, drop = TRUE)
  if (!any(vapply(pair_groups, function(index) length(unique(data$ploidy_category[index])) > 1L, logical(1L)))) {
    return(NA_real_)
  }
  permuted <- numeric(permutations)
  for (b in seq_len(permutations)) {
    shuffled <- data$ploidy_category
    for (index in pair_groups) shuffled[index] <- sample(shuffled[index], length(index), replace = FALSE)
    tab <- table(
      factor(shuffled, levels = c("CatA", "CatB", "CatC")),
      factor(data$ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
    )
    permuted[[b]] <- o2jca_cramers_v(tab)
  }
  (1 + sum(permuted >= observed, na.rm = TRUE)) / (permutations + 1)
}

metric_effect <- function(value, category) {
  keep <- is.finite(value) & category %in% c("CatA", "CatB", "CatC")
  value <- value[keep]; category <- droplevels(factor(category[keep]))
  if (length(value) < 3L || nlevels(category) < 2L) return(c(p_value = NA, epsilon_squared = NA))
  test <- stats::kruskal.test(value ~ category)
  h <- unname(test$statistic); k <- nlevels(category); n <- length(value)
  c(p_value = test$p.value, epsilon_squared = max((h - k + 1) / (n - k), 0))
}

analyze_ploidy_ratio_class_association <- function(out_dir, permutations = 999L) {
  data <- o2jca_read_tsv(file.path(out_dir, "ploidy_coupling_master_long.tsv"))
  data <- data[data$ploidy_category %in% c("CatA", "CatB", "CatC") & data$ratio_class %in% c("ClassA", "ClassB", "ClassC"), , drop = FALSE]
  seed_categories <- unique(data[, c("pair_id", "seed", "ploidy_category"), drop = FALSE])
  category_estimability <- do.call(rbind, lapply(split(seed_categories, seed_categories$pair_id, drop = TRUE), function(pair) {
    counts <- table(factor(pair$ploidy_category, levels = c("CatA", "CatB", "CatC")))
    data.frame(
      pair_id = pair$pair_id[[1L]], n_seed = nrow(pair), n_observed_categories = sum(counts > 0),
      n_CatA = unname(counts[["CatA"]]), n_CatB = unname(counts[["CatB"]]), n_CatC = unname(counts[["CatC"]]),
      within_pair_category_comparison_estimable = sum(counts > 0) >= 2L,
      stringsAsFactors = FALSE
    )
  }))

  pair_parameter <- do.call(rbind, lapply(o2jca_group_split(data, c("pair_id", "parameter")), function(group) {
    cat_counts <- table(factor(group$ploidy_category, levels = c("CatA", "CatB", "CatC")))
    class_counts <- table(factor(group$ratio_class, levels = c("ClassA", "ClassB", "ClassC")))
    n <- sum(class_counts)
    data.frame(
      pair_id = group$pair_id[[1L]], pair_label = o2jca_pair_short_label(group$pair_id[[1L]]),
      ploidy_category = names(cat_counts)[which.max(cat_counts)], parameter = group$parameter[[1L]],
      class_threshold = unique(suppressWarnings(as.numeric(group$class_threshold)))[[1L]],
      n_seed = length(unique(group$seed)),
      mean_vivo_natural = mean(suppressWarnings(as.numeric(group$vivo_natural)), na.rm = TRUE),
      mean_vitro_natural = mean(suppressWarnings(as.numeric(group$vitro_natural)), na.rm = TRUE),
      median_log2_ratio = stats::median(suppressWarnings(as.numeric(group$log2_ratio_vivo_to_vitro)), na.rm = TRUE),
      mean_log2_ratio = mean(suppressWarnings(as.numeric(group$log2_ratio_vivo_to_vitro)), na.rm = TRUE),
      q05_log2_ratio = o2jca_quantile(group$log2_ratio_vivo_to_vitro, 0.05),
      q95_log2_ratio = o2jca_quantile(group$log2_ratio_vivo_to_vitro, 0.95),
      prop_ClassA = unname(class_counts[["ClassA"]]) / n,
      prop_ClassB = unname(class_counts[["ClassB"]]) / n,
      prop_ClassC = unname(class_counts[["ClassC"]]) / n,
      dominant_ratio_class = names(class_counts)[which.max(class_counts)],
      dominant_ratio_class_fraction = max(class_counts) / n,
      stringsAsFactors = FALSE
    )
  }))
  pair_class_cells <- do.call(rbind, lapply(seq_len(nrow(pair_parameter)), function(i) {
    row <- pair_parameter[i, , drop = FALSE]
    row_core <- row[rep(1L, 3L), c("pair_id", "pair_label", "ploidy_category", "parameter"), drop = FALSE]
    rownames(row_core) <- NULL
    data.frame(
      row_core,
      ratio_class = c("ClassA", "ClassB", "ClassC"),
      proportion = c(row$prop_ClassA, row$prop_ClassB, row$prop_ClassC),
      stringsAsFactors = FALSE
    )
  }))
  pair_balanced_cells <- do.call(rbind, lapply(o2jca_group_split(pair_class_cells, c("ploidy_category", "parameter", "ratio_class")), function(group) {
    data.frame(
      ploidy_category = group$ploidy_category[[1L]], parameter = group$parameter[[1L]], ratio_class = group$ratio_class[[1L]],
      n_pairs = nrow(group), pair_balanced_mean_proportion = mean(group$proportion),
      pair_min = min(group$proportion), pair_max = max(group$proportion),
      stringsAsFactors = FALSE
    )
  }))
  pair_category_values <- do.call(rbind, lapply(o2jca_group_split(pair_parameter, c("ploidy_category", "parameter")), function(group) {
    data.frame(
      ploidy_category = group$ploidy_category[[1L]], parameter = group$parameter[[1L]], n_pairs = nrow(group),
      pair_balanced_mean_vivo = mean(group$mean_vivo_natural),
      pair_balanced_mean_vitro = mean(group$mean_vitro_natural),
      pair_balanced_mean_log2_ratio = mean(group$mean_log2_ratio),
      pair_min_log2_ratio = min(group$mean_log2_ratio), pair_max_log2_ratio = max(group$mean_log2_ratio),
      stringsAsFactors = FALSE
    )
  }))
  pair_level_association <- do.call(rbind, lapply(split(pair_parameter, pair_parameter$parameter, drop = TRUE), function(group) {
    tab <- table(
      factor(group$ploidy_category, levels = c("CatA", "CatB", "CatC")),
      factor(group$dominant_ratio_class, levels = c("ClassA", "ClassB", "ClassC"))
    )
    data.frame(
      parameter = group$parameter[[1L]], n_pairs = nrow(group),
      pair_level_cramers_v = o2jca_cramers_v(tab),
      pair_level_normalized_mutual_information = o2jca_normalized_mutual_information(tab),
      association_scope = "descriptive pair-level association; only six pairs and ploidy category is pair-confounded",
      stringsAsFactors = FALSE
    )
  }))
  global_rows <- list(); cell_rows <- list(); pair_rows <- list(); metric_rows <- list()
  metrics <- c("vivo_natural", "vitro_natural", "log2_ratio_vivo_to_vitro")
  for (parameter in o2jca_parameter_levels()) {
    group <- data[data$parameter == parameter, , drop = FALSE]
    details <- contingency_details(group)
    v <- o2jca_cramers_v(details$table)
    nmi <- o2jca_normalized_mutual_information(details$table)
    chi_table <- details$table[rowSums(details$table) > 0, colSums(details$table) > 0, drop = FALSE]
    chi <- if (min(dim(chi_table)) >= 2L) tryCatch(suppressWarnings(stats::chisq.test(chi_table, correct = FALSE)), error = function(e) NULL) else NULL
    pairs_with_variation <- sum(vapply(split(group$ploidy_category, group$pair_id, drop = TRUE), function(x) length(unique(x)) > 1L, logical(1L)))
    stratified_estimable <- pairs_with_variation > 0L
    array <- array(0, dim = c(3L, 3L, length(unique(group$pair_id))))
    for (i in seq_along(unique(group$pair_id))) {
      pair <- group[group$pair_id == unique(group$pair_id)[[i]], , drop = FALSE]
      array[, , i] <- contingency_details(pair)$table
    }
    cmh <- tryCatch(suppressWarnings(stats::mantelhaen.test(array)), error = function(e) NULL)
    perm_p <- stratified_permutation_p(group, v, permutations, 5826L + match(parameter, o2jca_parameter_levels()))
    global_rows[[length(global_rows) + 1L]] <- data.frame(
      parameter = parameter, n = nrow(group), n_pairs = length(unique(group$pair_id)),
      cramers_v = v, normalized_mutual_information = nmi,
      chi_square = if (is.null(chi)) NA_real_ else unname(chi$statistic),
      chi_square_p = if (is.null(chi)) NA_real_ else chi$p.value,
      pair_stratified_test_p = if (!stratified_estimable || is.null(cmh)) NA_real_ else cmh$p.value,
      pair_stratified_permutation_p = perm_p,
      n_pairs_with_within_pair_category_variation = pairs_with_variation,
      pair_stratified_association_estimable = stratified_estimable,
      association_scope = if (stratified_estimable) "pair-stratified model-derived association" else "descriptive cross-pair pattern; ploidy category is confounded with pair",
      permutations = permutations,
      stringsAsFactors = FALSE
    )
    detail <- details$details
    detail$parameter <- parameter
    cell_rows[[length(cell_rows) + 1L]] <- detail
    for (pair_id in unique(group$pair_id)) {
      pair <- group[group$pair_id == pair_id, , drop = FALSE]
      pair_details <- contingency_details(pair)
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        pair_id = pair_id, parameter = parameter, n = nrow(pair),
        cramers_v = o2jca_cramers_v(pair_details$table),
        normalized_mutual_information = o2jca_normalized_mutual_information(pair_details$table),
        stringsAsFactors = FALSE
      )
    }
    for (metric in metrics) {
      global_effect <- metric_effect(suppressWarnings(as.numeric(group[[metric]])), group$ploidy_category)
      pair_effects <- vapply(split(group, group$pair_id, drop = TRUE), function(pair) {
        metric_effect(suppressWarnings(as.numeric(pair[[metric]])), pair$ploidy_category)[["epsilon_squared"]]
      }, numeric(1L))
      metric_rows[[length(metric_rows) + 1L]] <- data.frame(
        parameter = parameter, metric = metric,
        global_kruskal_p = global_effect[["p_value"]], global_epsilon_squared = global_effect[["epsilon_squared"]],
        pair_balanced_mean_epsilon_squared = mean(pair_effects, na.rm = TRUE),
        pair_balanced_median_epsilon_squared = stats::median(pair_effects, na.rm = TRUE),
        n_pairs_estimable = sum(is.finite(pair_effects)),
        interpretation_scope = "model-derived association; not causal or external biological validation",
        stringsAsFactors = FALSE
      )
    }
  }
  global <- do.call(rbind, global_rows)
  global$fdr_permutation <- stats::p.adjust(global$pair_stratified_permutation_p, method = "BH")
  global$fdr_chi_square <- stats::p.adjust(global$chi_square_p, method = "BH")
  global$significance_status <- ifelse(
    !global$pair_stratified_association_estimable,
    "not_estimable_no_within_pair_category_variation",
    ifelse(global$fdr_permutation < 0.05, "pair_stratified_FDR_lt_0.05", "pair_stratified_FDR_ge_0.05")
  )
  cells <- do.call(rbind, cell_rows)
  pair_effects <- do.call(rbind, pair_rows)
  metric_comparison <- do.call(rbind, metric_rows)
  metric_comparison$fdr_global_kruskal <- ave(metric_comparison$global_kruskal_p, metric_comparison$metric, FUN = function(x) stats::p.adjust(x, method = "BH"))

  report_summary <- merge(global, o2jca_process_map(), by = "parameter", all.x = TRUE, sort = FALSE)
  report_summary <- report_summary[order(-report_summary$cramers_v), , drop = FALSE]
  files <- c(
    o2jca_write_tsv(global, file.path(out_dir, "cat_ratio_class_global_association.tsv")),
    o2jca_write_tsv(category_estimability, file.path(out_dir, "ploidy_category_estimability.tsv")),
    o2jca_write_tsv(cells, file.path(out_dir, "cat_ratio_class_contingency_cells.tsv")),
    o2jca_write_tsv(pair_effects, file.path(out_dir, "cat_ratio_class_within_pair_effects.tsv")),
    o2jca_write_tsv(metric_comparison, file.path(out_dir, "cat_parameter_metric_information_comparison.tsv")),
    o2jca_write_tsv(pair_parameter, file.path(out_dir, "pair_level_cat_parameter_summary.tsv")),
    o2jca_write_tsv(pair_class_cells, file.path(out_dir, "pair_level_cat_class_cells.tsv")),
    o2jca_write_tsv(pair_balanced_cells, file.path(out_dir, "cat_ratio_class_pair_balanced_cells.tsv")),
    o2jca_write_tsv(pair_category_values, file.path(out_dir, "cat_parameter_pair_balanced_values.tsv")),
    o2jca_write_tsv(pair_level_association, file.path(out_dir, "cat_ratio_class_pair_level_association.tsv")),
    o2jca_write_tsv(report_summary, paste0(out_dir, "/ploidy_coupling_report_summary.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "ploidy_ratio_class_association_manifest.tsv"))
  invisible(global)
}

main <- function() {
  args <- o2jca_parse_args(); out_dir <- args$out_dir %||% stop("--out_dir is required", call. = FALSE)
  analyze_ploidy_ratio_class_association(out_dir, o2jca_as_int(args$permutations, 999L))
}

if (identical(environment(), globalenv())) main()
