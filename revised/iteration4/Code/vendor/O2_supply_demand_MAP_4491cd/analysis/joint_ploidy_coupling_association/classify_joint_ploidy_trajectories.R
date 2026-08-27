#!/usr/bin/env Rscript

# Classify materialized joint-fit mean chromosome trajectories into CatA/B/C/U.
# No model or simulation code is invoked here.

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_coupling_analysis_utils.R"))
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_ploidy_category_utils.R"))

reclassify_from_features <- function(features, high_tolerance, low_endpoint, plateau_min_days, bic_cutoff) {
  groups <- o2jca_group_split(features, c("pair_id", "seed"))
  rows <- lapply(groups, function(group) {
    cohort_category <- rep("CatU", nrow(group))
    high <- group$minimum_value >= 44 - high_tolerance
    low <- group$terminal_median <= low_endpoint
    catc <- low & group$n_drop_episodes >= 2L & group$plateau_duration >= plateau_min_days &
      group$delta_bic_two_minus_one <= bic_cutoff
    catb <- low & group$n_drop_episodes >= 1L & !catc
    cohort_category[high] <- "CatA"
    cohort_category[catb] <- "CatB"
    cohort_category[catc] <- "CatC"
    consensus <- o2jpc_seed_consensus(stats::setNames(cohort_category, group$cohort))
    category <- consensus[["category"]]
    data.frame(pair_id = group$pair_id[[1L]], seed = group$seed[[1L]], category = category, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

classify_joint_ploidy_trajectories <- function(
    result_root,
    out_dir,
    pair_pattern = "^fit_joint_.*_vt_seed[0-9]+$",
    max_pairs = NA_integer_,
    max_seeds = NA_integer_,
    trajectory_step = 10L) {
  pair_dirs <- o2jca_discover_pair_dirs(result_root, pair_pattern, max_pairs)
  all_features <- list(); all_categories <- list(); plot_rows <- list(); quality_rows <- list()
  for (pair_dir in pair_dirs) {
    pair_id <- basename(pair_dir)
    path <- file.path(pair_dir, "extra_results", "predict1000_ploidy_seed_day_mean.tsv")
    message("Classifying ploidy trajectories for ", pair_id)
    data <- o2jca_read_tsv(path)
    o2jca_assert_columns(data, c("seed", "cohort", "day", "ploidy_value"), paste0(pair_id, " ploidy trajectories"))
    data$seed <- o2jca_normalize_seed(data$seed)
    data$seed_number <- o2jca_seed_number(data$seed)
    data$day <- suppressWarnings(as.numeric(data$day))
    data$ploidy_value <- suppressWarnings(as.numeric(data$ploidy_value))
    if (is.finite(max_seeds)) data <- data[data$seed_number <= max_seeds, , drop = FALSE]
    data <- data[data$day >= 0 & data$day <= 1000 & data$cohort %in% c("2N", "4N"), , drop = FALSE]
    data <- data[order(data$seed_number, data$cohort, data$day), , drop = FALSE]
    group_key <- paste(data$seed, data$cohort, sep = "\r")
    groups <- split(seq_len(nrow(data)), group_key, drop = TRUE)
    feature_rows <- lapply(groups, function(index) {
      group <- data[index, , drop = FALSE]
      data.frame(
        pair_id = pair_id, seed = group$seed[[1L]], seed_number = group$seed_number[[1L]], cohort = group$cohort[[1L]],
        o2jpc_curve_features(group$day, group$ploidy_value),
        stringsAsFactors = FALSE
      )
    })
    features <- do.call(rbind, feature_rows)
    seed_rows <- lapply(split(features, features$seed, drop = TRUE), function(group) {
      consensus <- o2jpc_seed_consensus(stats::setNames(group$category, group$cohort))
      data.frame(
        pair_id = pair_id, seed = group$seed[[1L]], seed_number = group$seed_number[[1L]],
        ploidy_category = consensus[["category"]], category_reason = consensus[["reason"]],
        category_2N = group$category[match("2N", group$cohort)],
        category_4N = group$category[match("4N", group$cohort)],
        terminal_2N = group$terminal_median[match("2N", group$cohort)],
        terminal_4N = group$terminal_median[match("4N", group$cohort)],
        stringsAsFactors = FALSE
      )
    })
    categories <- do.call(rbind, seed_rows)
    category_map <- categories[, c("seed", "ploidy_category"), drop = FALSE]
    selected <- data[data$day %% trajectory_step == 0, c("seed", "seed_number", "cohort", "day", "ploidy_value"), drop = FALSE]
    selected$pair_id <- pair_id
    selected$ploidy_category <- category_map$ploidy_category[match(selected$seed, category_map$seed)]
    all_features[[length(all_features) + 1L]] <- features
    all_categories[[length(all_categories) + 1L]] <- categories
    plot_rows[[length(plot_rows) + 1L]] <- selected
    expected_rows <- length(unique(data$seed)) * 2L * 1001L
    quality_rows[[length(quality_rows) + 1L]] <- data.frame(
      pair_id = pair_id, n_seed = length(unique(data$seed)), n_cohort = length(unique(data$cohort)),
      n_rows = nrow(data), expected_rows = expected_rows, row_count_pass = nrow(data) == expected_rows,
      n_nonfinite = sum(!is.finite(data$day) | !is.finite(data$ploidy_value)),
      n_uncertain = sum(categories$ploidy_category == "CatU"),
      stringsAsFactors = FALSE
    )
  }
  features <- do.call(rbind, all_features)
  categories <- do.call(rbind, all_categories)
  plot_data <- do.call(rbind, plot_rows)
  quality <- do.call(rbind, quality_rows)
  o2jca_assert_unique_key(features, c("pair_id", "seed", "cohort"), "ploidy cohort features")
  o2jca_assert_unique_key(categories, c("pair_id", "seed"), "ploidy seed categories")
  if (!all(quality$row_count_pass) || any(quality$n_nonfinite > 0)) stop("Ploidy trajectory input contract failed", call. = FALSE)

  category_summary <- as.data.frame(table(categories$pair_id, factor(categories$ploidy_category, levels = c("CatA", "CatB", "CatC", "CatU"))), stringsAsFactors = FALSE)
  names(category_summary) <- c("pair_id", "ploidy_category", "n_seed")
  category_summary$proportion <- ave(category_summary$n_seed, category_summary$pair_id, FUN = function(x) x / sum(x))

  trajectory_groups <- o2jca_group_split(plot_data, c("pair_id", "ploidy_category", "cohort", "day"))
  trajectory_summary <- do.call(rbind, lapply(trajectory_groups, function(group) {
    data.frame(
      pair_id = group$pair_id[[1L]], ploidy_category = group$ploidy_category[[1L]],
      cohort = group$cohort[[1L]], day = group$day[[1L]], n_seed = nrow(group),
      median_ploidy = stats::median(group$ploidy_value, na.rm = TRUE),
      q25_ploidy = o2jca_quantile(group$ploidy_value, 0.25),
      q75_ploidy = o2jca_quantile(group$ploidy_value, 0.75),
      stringsAsFactors = FALSE
    )
  }))

  sensitivity_rows <- list()
  grid <- expand.grid(
    high_tolerance = c(0.25, 0.5, 1), low_endpoint = c(28, 30, 32),
    plateau_min_days = c(40, 60, 80), bic_delta_cutoff = c(-6, -10, -14),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(grid))) {
    classified <- reclassify_from_features(features, grid$high_tolerance[[i]], grid$low_endpoint[[i]], grid$plateau_min_days[[i]], grid$bic_delta_cutoff[[i]])
    counts <- as.data.frame(table(classified$pair_id, factor(classified$category, levels = c("CatA", "CatB", "CatC", "CatU"))), stringsAsFactors = FALSE)
    names(counts) <- c("pair_id", "ploidy_category", "n_seed")
    counts$proportion <- ave(counts$n_seed, counts$pair_id, FUN = function(x) x / sum(x))
    counts$high_tolerance <- grid$high_tolerance[[i]]
    counts$low_endpoint <- grid$low_endpoint[[i]]
    counts$plateau_min_days <- grid$plateau_min_days[[i]]
    counts$bic_delta_cutoff <- grid$bic_delta_cutoff[[i]]
    counts$is_primary <- grid$high_tolerance[[i]] == 0.5 && grid$low_endpoint[[i]] == 30 && grid$plateau_min_days[[i]] == 60 && grid$bic_delta_cutoff[[i]] == -10
    sensitivity_rows[[length(sensitivity_rows) + 1L]] <- counts
  }
  sensitivity <- do.call(rbind, sensitivity_rows)

  pair_assignment <- do.call(rbind, lapply(split(categories, categories$pair_id, drop = TRUE), function(group) {
    counts <- table(factor(group$ploidy_category, levels = c("CatA", "CatB", "CatC", "CatU")))
    winner <- names(counts)[which.max(counts)]
    data.frame(
      pair_id = group$pair_id[[1L]], pair_label = o2jca_pair_short_label(group$pair_id[[1L]]),
      pair_ploidy_category = winner, n_seed = nrow(group),
      dominant_fraction = max(counts) / sum(counts),
      n_observed_categories = sum(counts > 0),
      within_pair_category_comparison_estimable = sum(counts > 0) >= 2L,
      stringsAsFactors = FALSE
    )
  }))

  archetype_groups <- o2jca_group_split(trajectory_summary, c("ploidy_category", "cohort", "day"))
  archetypes <- do.call(rbind, lapply(archetype_groups, function(group) {
    data.frame(
      ploidy_category = group$ploidy_category[[1L]], cohort = group$cohort[[1L]], day = group$day[[1L]],
      n_pairs = nrow(group), pair_balanced_median = mean(group$median_ploidy, na.rm = TRUE),
      pair_q25 = o2jca_quantile(group$median_ploidy, 0.25), pair_q75 = o2jca_quantile(group$median_ploidy, 0.75),
      stringsAsFactors = FALSE
    )
  }))

  representative_seeds <- do.call(rbind, lapply(split(categories, categories$pair_id, drop = TRUE), function(group) {
    z2 <- abs(group$terminal_2N - stats::median(group$terminal_2N, na.rm = TRUE))
    z4 <- abs(group$terminal_4N - stats::median(group$terminal_4N, na.rm = TRUE))
    scale2 <- stats::mad(group$terminal_2N, constant = 1, na.rm = TRUE); if (!is.finite(scale2) || scale2 == 0) scale2 <- 1
    scale4 <- stats::mad(group$terminal_4N, constant = 1, na.rm = TRUE); if (!is.finite(scale4) || scale4 == 0) scale4 <- 1
    score <- z2 / scale2 + z4 / scale4
    chosen <- group[order(score, group$seed_number), , drop = FALSE][1L, ]
    data.frame(
      pair_id = chosen$pair_id, seed = chosen$seed, seed_number = chosen$seed_number,
      ploidy_category = chosen$ploidy_category, representative_rule = "closest_to_pair_median_terminal_2N_and_4N",
      stringsAsFactors = FALSE
    )
  }))
  representative <- merge(plot_data, representative_seeds, by = c("pair_id", "seed", "seed_number", "ploidy_category"), sort = FALSE)
  feature_annotations <- features[, c(
    "pair_id", "seed", "cohort", "first_drop_start_day", "first_drop_end_day",
    "second_drop_start_day", "second_drop_end_day", "plateau_duration",
    "delta_bic_two_minus_one", "category"
  )]
  representative <- merge(representative, feature_annotations, by = c("pair_id", "seed", "cohort"), all.x = TRUE, sort = FALSE)

  confidence <- features
  confidence$margin_high_floor_chr <- confidence$minimum_value - (44 - 0.5)
  confidence$margin_low_endpoint_chr <- 30 - confidence$terminal_median
  confidence$margin_plateau_days <- confidence$plateau_duration - 60
  confidence$margin_two_transition_bic <- -10 - confidence$delta_bic_two_minus_one
  confidence <- confidence[, c(
    "pair_id", "seed", "seed_number", "cohort", "category",
    "margin_high_floor_chr", "margin_low_endpoint_chr", "margin_plateau_days",
    "margin_two_transition_bic", "minimum_value", "terminal_median",
    "plateau_duration", "delta_bic_two_minus_one"
  )]

  sensitivity$key <- do.call(paste, c(sensitivity[, c("pair_id", "high_tolerance", "low_endpoint", "plateau_min_days", "bic_delta_cutoff")], sep = "\r"))
  sensitivity_agreement <- do.call(rbind, lapply(split(sensitivity, sensitivity$key, drop = TRUE), function(group) {
    group <- group[order(match(group$ploidy_category, c("CatA", "CatB", "CatC", "CatU"))), , drop = FALSE]
    winner <- group$ploidy_category[[which.max(group$n_seed)]]
    primary <- pair_assignment$pair_ploidy_category[match(group$pair_id[[1L]], pair_assignment$pair_id)]
    data.frame(
      pair_id = group$pair_id[[1L]], pair_label = o2jca_pair_short_label(group$pair_id[[1L]]),
      high_tolerance = group$high_tolerance[[1L]], low_endpoint = group$low_endpoint[[1L]],
      plateau_min_days = group$plateau_min_days[[1L]], bic_delta_cutoff = group$bic_delta_cutoff[[1L]],
      scenario_category = winner, primary_category = primary,
      dominant_fraction = max(group$n_seed) / sum(group$n_seed), agrees_with_primary = winner == primary,
      is_primary = group$is_primary[[1L]], stringsAsFactors = FALSE
    )
  }))
  sensitivity$key <- NULL
  sensitivity_pair_summary <- do.call(rbind, lapply(split(sensitivity_agreement, sensitivity_agreement$pair_id, drop = TRUE), function(group) {
    counts <- table(factor(group$scenario_category, levels = c("CatA", "CatB", "CatC", "CatU")))
    do.call(rbind, lapply(names(counts), function(category) {
      data.frame(
        pair_id = group$pair_id[[1L]], pair_label = group$pair_label[[1L]],
        primary_category = group$primary_category[[1L]], scenario_category = category,
        n_scenarios = nrow(group), n_scenarios_in_category = unname(counts[[category]]),
        proportion_scenarios = unname(counts[[category]]) / nrow(group),
        overall_primary_agreement_rate = mean(group$agrees_with_primary),
        stringsAsFactors = FALSE
      )
    }))
  }))
  files <- c(
    o2jca_write_tsv(features, file.path(out_dir, "ploidy_cohort_features.tsv")),
    o2jca_write_tsv(categories, file.path(out_dir, "ploidy_seed_categories.tsv")),
    o2jca_write_tsv(category_summary, file.path(out_dir, "ploidy_category_summary.tsv")),
    o2jca_write_tsv(plot_data, file.path(out_dir, "ploidy_trajectory_plot_data.tsv")),
    o2jca_write_tsv(trajectory_summary, file.path(out_dir, "ploidy_trajectory_category_summary.tsv")),
    o2jca_write_tsv(sensitivity, file.path(out_dir, "ploidy_category_sensitivity.tsv")),
    o2jca_write_tsv(sensitivity_agreement, file.path(out_dir, "ploidy_category_sensitivity_agreement.tsv")),
    o2jca_write_tsv(sensitivity_pair_summary, file.path(out_dir, "ploidy_category_sensitivity_pair_summary.tsv")),
    o2jca_write_tsv(o2jpc_category_definition(), file.path(out_dir, "ploidy_category_definition.tsv")),
    o2jca_write_tsv(quality, file.path(out_dir, "ploidy_input_quality_summary.tsv")),
    o2jca_write_tsv(pair_assignment, file.path(out_dir, "ploidy_pair_category_assignment.tsv")),
    o2jca_write_tsv(archetypes, file.path(out_dir, "ploidy_category_archetype_summary.tsv")),
    o2jca_write_tsv(representative, file.path(out_dir, "ploidy_representative_trajectory_plot_data.tsv")),
    o2jca_write_tsv(confidence, file.path(out_dir, "ploidy_classification_confidence.tsv"))
  )
  o2jca_write_manifest(files, file.path(out_dir, "ploidy_classification_manifest.tsv"))
  invisible(categories)
}

main <- function() {
  args <- o2jca_parse_args(); result_root <- args$result_root %||% stop("--result_root is required", call. = FALSE)
  output_root <- args$output_root %||% o2jca_default_output_root(result_root, WORKFLOW_ROOT)
  out_dir <- args$out_dir %||% file.path(output_root, "tables", "ploidy_coupling")
  o2jca_assert_separate_output_root(result_root, out_dir)
  classify_joint_ploidy_trajectories(
    result_root = result_root, out_dir = out_dir,
    pair_pattern = args$pair_pattern %||% "^fit_joint_.*_vt_seed[0-9]+$",
    max_pairs = o2jca_as_int(args$max_pairs, NA_integer_),
    max_seeds = o2jca_as_int(args$max_seeds, NA_integer_),
    trajectory_step = o2jca_as_int(args$trajectory_step, 10L)
  )
}

if (identical(environment(), globalenv())) main()
