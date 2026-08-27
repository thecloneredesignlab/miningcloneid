#!/usr/bin/env Rscript

.o2pl_contribution_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "parameter_contribution_analysis.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(.o2pl_contribution_dir, "parameter_landscape_analysis_utils.R"), local = environment(), chdir = TRUE)

o2pl_contribution_target <- function(target) {
  target <- tolower(trimws(as.character(target %||% "mode")))
  if (target %in% c("mode", "discrete", "classification", "logistic")) return("mode")
  if (target %in% c("dominant_ploidy", "dominant_mean_ploidy", "ploidy", "continuous", "regression")) return("dominant_ploidy")
  stop("Contribution target must be mode or dominant_ploidy.", call. = FALSE)
}

o2pl_contribution_output_dir <- function(root_dir, target) {
  file.path(root_dir, if (o2pl_contribution_target(target) == "mode") "mode_parameter_contribution" else "dominant_ploidy_parameter_contribution")
}

o2pl_auc <- function(y, score) {
  keep <- is.finite(score) & !is.na(y)
  y <- as.integer(y[keep])
  score <- as.numeric(score[keep])
  if (length(unique(y)) != 2L) return(NA_real_)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  (sum(rank(score, ties.method = "average")[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

o2pl_reference_response <- function(root_dir, reference_o2) {
  path <- file.path(paper_fixo2_mode_tables_dir(root_dir), "fixed_o2_attractor_mode_by_seed_o2.tsv")
  if (!file.exists(path)) stop("Missing materialized fixed-O2 response table: ", path, call. = FALSE)
  table <- read_tsv(path)
  o2pl_require_columns(table, c("seed_id", "O2_pct", "dominant_mean_ploidy"), path)
  table$O2_pct <- suppressWarnings(as.numeric(table$O2_pct))
  table <- table[is.finite(table$O2_pct) & abs(table$O2_pct - reference_o2) < 1e-9, , drop = FALSE]
  if (!nrow(table)) stop("No response rows at O2=", reference_o2, call. = FALSE)
  table$seed <- suppressWarnings(as.integer(sub("^seed", "", table$seed_id)))
  table$mode_label <- if ("mode_label" %in% names(table)) as.character(table$mode_label) else ifelse(as.numeric(table$dominant_mean_ploidy) >= 2, "mode1", "mode2")
  table[!duplicated(table$seed), , drop = FALSE]
}

o2pl_mode_main_effects <- function(features, y) {
  rows <- lapply(names(features), function(parameter) {
    x <- as.numeric(features[[parameter]])
    fit <- tryCatch(stats::glm(y ~ x, family = stats::binomial()), error = function(e) NULL)
    coefficient <- if (!is.null(fit) && nrow(summary(fit)$coefficients) >= 2L) summary(fit)$coefficients[2L, 1L] else NA_real_
    p_value <- if (!is.null(fit) && nrow(summary(fit)$coefficients) >= 2L) summary(fit)$coefficients[2L, 4L] else NA_real_
    auc <- o2pl_auc(y, x)
    data.frame(
      feature_name = parameter,
      feature_type = "main",
      parameter = parameter,
      parameter_a = parameter,
      parameter_b = NA_character_,
      coefficient = coefficient,
      p_value = p_value,
      auc = auc,
      discriminative_auc = if (is.finite(auc)) max(auc, 1 - auc) else NA_real_,
      direction = if (is.finite(coefficient) && coefficient >= 0) "mode1_higher" else "mode2_higher",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

o2pl_mode_interactions <- function(features, y) {
  pairs <- utils::combn(names(features), 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    a <- as.numeric(features[[pair[[1L]]]])
    b <- as.numeric(features[[pair[[2L]]]])
    fit <- tryCatch(stats::glm(y ~ a + b + a:b, family = stats::binomial()), error = function(e) NULL)
    coefficients <- if (!is.null(fit)) summary(fit)$coefficients else matrix(numeric(), 0, 0)
    term <- grep("a:b", rownames(coefficients), fixed = TRUE)
    coefficient <- if (length(term)) coefficients[term[[1L]], 1L] else NA_real_
    p_value <- if (length(term)) coefficients[term[[1L]], 4L] else NA_real_
    score <- a * b
    auc <- o2pl_auc(y, score)
    data.frame(
      feature_name = paste(pair, collapse = " x "),
      feature_type = "interaction",
      parameter = NA_character_,
      parameter_a = pair[[1L]],
      parameter_b = pair[[2L]],
      coefficient = coefficient,
      p_value = p_value,
      auc = auc,
      discriminative_auc = if (is.finite(auc)) max(auc, 1 - auc) else NA_real_,
      direction = if (is.finite(coefficient) && coefficient >= 0) "mode1_higher" else "mode2_higher",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

o2pl_continuous_main_effects <- function(features, y) {
  rows <- lapply(names(features), function(parameter) {
    x <- as.numeric(features[[parameter]])
    fit <- stats::lm(y ~ x)
    coefficients <- summary(fit)$coefficients
    data.frame(
      feature_name = parameter,
      feature_type = "main",
      parameter = parameter,
      parameter_a = parameter,
      parameter_b = NA_character_,
      coefficient = coefficients[2L, 1L],
      p_value = coefficients[2L, 4L],
      r2 = summary(fit)$r.squared,
      spearman_rho = suppressWarnings(stats::cor(x, y, method = "spearman")),
      direction = if (coefficients[2L, 1L] >= 0) "positive" else "negative",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

o2pl_continuous_interactions <- function(features, y) {
  pairs <- utils::combn(names(features), 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    a <- as.numeric(features[[pair[[1L]]]])
    b <- as.numeric(features[[pair[[2L]]]])
    main <- stats::lm(y ~ a + b)
    full <- stats::lm(y ~ a + b + a:b)
    coefficients <- summary(full)$coefficients
    term <- grep("a:b", rownames(coefficients), fixed = TRUE)
    data.frame(
      feature_name = paste(pair, collapse = " x "),
      feature_type = "interaction",
      parameter = NA_character_,
      parameter_a = pair[[1L]],
      parameter_b = pair[[2L]],
      coefficient = coefficients[term[[1L]], 1L],
      p_value = coefficients[term[[1L]], 4L],
      r2 = summary(full)$r.squared,
      delta_r2 = summary(full)$r.squared - summary(main)$r.squared,
      direction = if (coefficients[term[[1L]], 1L] >= 0) "positive" else "negative",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out
}

o2pl_reference_contribution <- function(root_dir, target, reference_o2, max_seeds = NA_integer_) {
  target <- o2pl_contribution_target(target)
  fitted <- o2pl_read_materialized_features(root_dir, "invivo", "fitted")
  response <- o2pl_reference_response(root_dir, reference_o2)
  joined <- merge(fitted, response[, c("seed", "seed_id", "mode_label", "dominant_mean_ploidy"), drop = FALSE], by = "seed", all = FALSE, sort = FALSE)
  if (is.finite(max_seeds) && max_seeds > 0L) joined <- head(joined[order(joined$seed), , drop = FALSE], max_seeds)
  transformed <- o2pl_transform_features(joined, "invivo")
  standardized <- as.data.frame(o2pl_scale_features(transformed, "zscore")$matrix, check.names = FALSE)
  output_dir <- o2pl_contribution_output_dir(root_dir, target)
  reference_dir <- file.path(output_dir, paste0("reference_o2_", fixed_o2_o2_slug(reference_o2)))
  dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
  if (target == "mode") {
    y <- as.integer(joined$mode_label == "mode1")
    if (length(unique(y)) < 2L) stop("Mode response has only one class at O2=", reference_o2, call. = FALSE)
    main <- o2pl_mode_main_effects(standardized, y)
    interactions <- o2pl_mode_interactions(standardized, y)
    combined <- rbind_fill_plain(list(main, interactions))
    combined$contribution_score <- combined$discriminative_auc
    combined <- combined[order(-combined$contribution_score, combined$fdr, combined$feature_name), , drop = FALSE]
    combined$rank <- seq_len(nrow(combined))
    seed_table <- data.frame(seed = joined$seed, seed_id = joined$seed_id, mode_label = joined$mode_label, mode_reference_o2 = reference_o2, stringsAsFactors = FALSE)
    write_csv(seed_table, file.path(reference_dir, "mode_parameter_seed_labels.csv"))
    write_csv(main, file.path(reference_dir, "mode_parameter_main_effects.csv"))
    write_csv(interactions, file.path(reference_dir, "mode_parameter_pairwise_interactions.csv"))
    write_csv(combined, file.path(reference_dir, "mode_parameter_feature_importance.csv"))
    summary <- data.frame(mode_reference_o2 = reference_o2, n_seed = nrow(joined), n_mode1 = sum(y == 1L), n_mode2 = sum(y == 0L), top_feature = combined$feature_name[[1L]], stringsAsFactors = FALSE)
    summary_name <- "mode_parameter_contribution_summary.tsv"
  } else {
    y <- suppressWarnings(as.numeric(joined$dominant_mean_ploidy))
    main <- o2pl_continuous_main_effects(standardized, y)
    interactions <- o2pl_continuous_interactions(standardized, y)
    combined <- rbind_fill_plain(list(main, interactions))
    combined$contribution_score <- ifelse(combined$feature_type == "main", combined$r2, combined$delta_r2)
    combined <- combined[order(-combined$contribution_score, combined$fdr, combined$feature_name), , drop = FALSE]
    combined$rank <- seq_len(nrow(combined))
    seed_table <- data.frame(seed = joined$seed, seed_id = joined$seed_id, mode_reference_dominant_mean_ploidy = y, mode_reference_o2 = reference_o2, stringsAsFactors = FALSE)
    write_csv(seed_table, file.path(reference_dir, "dominant_ploidy_parameter_seed_values.csv"))
    write_csv(main, file.path(reference_dir, "dominant_ploidy_parameter_main_effects.csv"))
    write_csv(interactions, file.path(reference_dir, "dominant_ploidy_parameter_pairwise_interactions.csv"))
    write_csv(combined, file.path(reference_dir, "dominant_ploidy_parameter_feature_importance.csv"))
    summary <- data.frame(mode_reference_o2 = reference_o2, n_seed = nrow(joined), response_mean = mean(y), response_sd = stats::sd(y), top_feature = combined$feature_name[[1L]], stringsAsFactors = FALSE)
    summary_name <- "dominant_ploidy_parameter_contribution_summary.tsv"
  }
  summary_path <- file.path(reference_dir, summary_name)
  write_tsv_plain(summary, summary_path)
  schema <- file.path(reference_dir, sub("[.]tsv$", ".schema.tsv", summary_name))
  o2pl_write_schema(summary, paste(target, fixed_o2_o2_slug(reference_o2), "contribution_summary", sep = "_"), schema)
  o2pl_record_artifact(root_dir, "analysis", paste(target, fixed_o2_o2_slug(reference_o2), "contribution", sep = "_"), summary_path, summary, schema, "parameter_contribution_analysis.R", "invivo", paste(o2pl_materialized_feature_path(root_dir, "invivo", "fitted"), file.path(paper_fixo2_mode_tables_dir(root_dir), "fixed_o2_attractor_mode_by_seed_o2.tsv"), sep = ";"))
  list(summary = summary, importance = head(combined, 30L), output_dir = reference_dir)
}

o2pl_merge_contributions <- function(root_dir, target, reference_o2) {
  target <- o2pl_contribution_target(target)
  output_dir <- o2pl_contribution_output_dir(root_dir, target)
  runs <- lapply(reference_o2, function(o2) {
    reference_dir <- file.path(output_dir, paste0("reference_o2_", fixed_o2_o2_slug(o2)))
    summary_name <- if (target == "mode") "mode_parameter_contribution_summary.tsv" else "dominant_ploidy_parameter_contribution_summary.tsv"
    importance_name <- if (target == "mode") "mode_parameter_feature_importance.csv" else "dominant_ploidy_parameter_feature_importance.csv"
    list(summary = read_tsv(file.path(reference_dir, summary_name)), importance = head(read_csv_plain(file.path(reference_dir, importance_name)), 30L), output_dir = reference_dir)
  })
  index <- do.call(rbind, lapply(runs, function(run) cbind(run$summary, output_dir = run$output_dir, stringsAsFactors = FALSE)))
  importance <- do.call(rbind, lapply(seq_along(runs), function(i) cbind(mode_reference_o2 = reference_o2[[i]], runs[[i]]$importance, stringsAsFactors = FALSE)))
  index_name <- if (target == "mode") "mode_parameter_contribution_index.csv" else "dominant_ploidy_parameter_contribution_index.csv"
  importance_name <- if (target == "mode") "mode_parameter_top_features_across_reference_o2.csv" else "dominant_ploidy_parameter_top_features_across_reference_o2.csv"
  write_csv(index, file.path(output_dir, index_name))
  write_csv(importance, file.path(output_dir, importance_name))
  invisible(file.path(output_dir, index_name))
}

o2pl_contribution_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE)), forced_target = NULL) {
  root_dir <- normalizePath(path.expand(argv$result_root %||% o2pl_default_result_root()), mustWork = FALSE)
  target <- o2pl_contribution_target(forced_target %||% argv$mode_contribution_target %||% argv$contribution_target %||% "mode")
  reference_o2 <- sort(unique(as_num_vec(argv$mode_reference_o2_values, as_num(argv$mode_reference_o2, 2))))
  merge_only <- as_bool(argv$merge_only, FALSE)
  if (!merge_only) {
    invisible(lapply(reference_o2, function(o2) o2pl_reference_contribution(root_dir, target, o2, as_int(argv$max_seeds, NA_integer_))))
  }
  o2pl_merge_contributions(root_dir, target, reference_o2)
  message("Parameter contribution analysis complete: ", o2pl_contribution_output_dir(root_dir, target))
}

if (sys.nframe() == 0L) o2pl_contribution_main()
