#!/usr/bin/env Rscript

# Canonical shared ploidy-regime helpers used by simulation and analysis.

suppressPackageStartupMessages({
  if (requireNamespace("Matrix", quietly = TRUE)) library(Matrix)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

o2pr_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

o2pr_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

o2pr_mkdirs <- function(out_dir) {
  invisible(vapply(
    file.path(out_dir, c("tables", "cache", "logs")),
    dir.create,
    logical(1),
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

o2pr_crossing_time <- function(day, value, threshold, direction = "down", after_day = -Inf) {
  ok <- is.finite(day) & is.finite(value) & day >= after_day
  x <- day[ok]
  y <- value[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  hit <- if (identical(direction, "down")) which(y <= threshold) else which(y >= threshold)
  if (!length(hit)) return(NA_real_)
  i <- hit[[1]]
  if (i == 1L) return(x[[1]])
  x0 <- x[[i - 1L]]
  x1 <- x[[i]]
  y0 <- y[[i - 1L]]
  y1 <- y[[i]]
  if (!is.finite(y1 - y0) || abs(y1 - y0) < 1e-12) return(x1)
  x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)
}

o2pr_auc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

o2pr_cohort_features <- function(seed_id, cohort, d, mid_window, late_window, nmin) {
  d <- d[order(d$day), , drop = FALSE]
  day <- as.numeric(d$day)
  y <- as.numeric(d$trajectory_value)
  mid <- y[day >= mid_window[[1]] & day <= mid_window[[2]]]
  late <- y[day >= late_window[[1]] & day <= late_window[[2]]]
  late_day <- day[day >= late_window[[1]] & day <= late_window[[2]]]
  late_level <- if (length(late)) stats::median(late, na.rm = TRUE) else NA_real_
  mid_reference <- if (length(mid)) as.numeric(stats::quantile(mid, 0.9, na.rm = TRUE, names = FALSE)) else NA_real_
  late_slope <- if (length(late) >= 2L) coef(stats::lm(late ~ late_day))[[2]] else NA_real_
  post <- day >= mid_window[[2]]
  dy <- diff(y[post])
  dx <- diff(day[post])
  slopes <- dy / dx
  largest_negative <- if (length(slopes)) min(slopes, na.rm = TRUE) else NA_real_
  t_largest_negative <- if (length(slopes) && is.finite(largest_negative)) {
    day[post][-1L][which.min(slopes)]
  } else {
    NA_real_
  }
  data.frame(
    seed_id = seed_id,
    cohort = cohort,
    late_level = late_level,
    mid_reference = mid_reference,
    late_drop_amplitude = mid_reference - late_level,
    late_slope = late_slope,
    largest_negative_slope_after_midpoint = largest_negative,
    time_of_largest_negative_slope = t_largest_negative,
    time_crossing_N35 = o2pr_crossing_time(day, y, 35, "down", mid_window[[1]]),
    time_crossing_N30 = o2pr_crossing_time(day, y, 30, "down", mid_window[[1]]),
    time_crossing_N25 = o2pr_crossing_time(day, y, 25, "down", mid_window[[1]]),
    minimum_late_value = if (length(late)) min(late, na.rm = TRUE) else NA_real_,
    late_fraction_at_or_near_Nmin = if (length(late)) mean(late <= nmin + 2, na.rm = TRUE) else NA_real_,
    terminal_distance_from_Nmin = late_level - nmin,
    trajectory_AUC = o2pr_auc(day, y),
    initial_to_terminal_change = if (length(y)) tail(y, 1) - y[[1]] else NA_real_,
    stringsAsFactors = FALSE
  )
}

o2pr_features_wide <- function(feat_long, labels) {
  if (!nrow(feat_long)) return(data.frame())
  ids <- unique(feat_long$seed_id)
  out <- data.frame(seed_id = ids, stringsAsFactors = FALSE)
  metric_cols <- setdiff(names(feat_long), c("seed_id", "cohort"))
  for (co in c("2N", "4N")) {
    sub <- feat_long[feat_long$cohort == co, , drop = FALSE]
    for (m in metric_cols) {
      out[[paste0(co, "__", m)]] <- sub[[m]][match(out$seed_id, sub$seed_id)]
    }
  }
  if (all(c("2N", "4N") %in% feat_long$cohort)) {
    l2 <- feat_long[feat_long$cohort == "2N", , drop = FALSE]
    l4 <- feat_long[feat_long$cohort == "4N", , drop = FALSE]
    out$abs_terminal_difference_2N_4N <- abs(l2$late_level[match(out$seed_id, l2$seed_id)] - l4$late_level[match(out$seed_id, l4$seed_id)])
    out$late_drop_mean <- rowMeans(cbind(l2$late_drop_amplitude[match(out$seed_id, l2$seed_id)], l4$late_drop_amplitude[match(out$seed_id, l4$seed_id)]), na.rm = TRUE)
    out$late_drop_min <- pmin(l2$late_drop_amplitude[match(out$seed_id, l2$seed_id)], l4$late_drop_amplitude[match(out$seed_id, l4$seed_id)], na.rm = TRUE)
    out$both_cohorts_cross_N35 <- is.finite(l2$time_crossing_N35[match(out$seed_id, l2$seed_id)]) & is.finite(l4$time_crossing_N35[match(out$seed_id, l4$seed_id)])
    out$both_cohorts_near_Nmin <- l2$late_fraction_at_or_near_Nmin[match(out$seed_id, l2$seed_id)] > 0.5 & l4$late_fraction_at_or_near_Nmin[match(out$seed_id, l4$seed_id)] > 0.5
    out$cohort_convergence_score <- 1 / (1 + out$abs_terminal_difference_2N_4N)
  }
  if (nrow(labels)) out <- merge(out, labels[, c("seed_id", "trajectory_regime", "trajectory_drop_score")], by = "seed_id", all.x = TRUE)
  out
}

o2pr_common_grid_matrix <- function(curve_all, trajectory_end_day) {
  if (!nrow(curve_all)) return(list(matrix = matrix(numeric(0), 0, 0), grid = numeric(0)))
  seeds <- sort(unique(curve_all$seed_id))
  grid <- seq(0, trajectory_end_day, by = 1)
  mat <- matrix(NA_real_, nrow = length(seeds), ncol = length(grid) * 2L)
  rownames(mat) <- seeds
  colnames(mat) <- c(paste0("2N_day", grid), paste0("4N_day", grid))
  for (seed in seeds) {
    vals <- numeric(0)
    ok_seed <- TRUE
    for (co in c("2N", "4N")) {
      d <- curve_all[curve_all$seed_id == seed & curve_all$cohort == co, , drop = FALSE]
      d <- d[order(d$day), , drop = FALSE]
      if (!nrow(d) || min(d$day) > min(grid) || max(d$day) < max(grid)) {
        ok_seed <- FALSE
        vals <- c(vals, rep(NA_real_, length(grid)))
      } else {
        vals <- c(vals, as.numeric(stats::approx(d$day, d$trajectory_value, xout = grid, ties = "ordered", rule = 1)$y))
      }
    }
    mat[seed, ] <- vals
    if (!ok_seed) mat[seed, ] <- NA_real_
  }
  mat <- mat[rowSums(!is.finite(mat)) == 0, , drop = FALSE]
  list(matrix = mat, grid = grid)
}

o2pr_bootstrap_trajectory_stability <- function(mat, k = 2L, n_bootstrap = 100L, random_seed = 1L) {
  if (nrow(mat) > 0L && ncol(mat) > 0L) {
    sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
    mat <- mat[, is.finite(sds) & sds > 0, drop = FALSE]
  }
  if (nrow(mat) < 3L || ncol(mat) < 2L || n_bootstrap <= 0L) {
    ids <- rownames(mat)
    if (is.null(ids)) ids <- character(0)
    return(data.frame(seed_id = ids, trajectory_cluster_stability = rep(NA_real_, length(ids)), stringsAsFactors = FALSE))
  }
  set.seed(random_seed)
  pc <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
  base <- stats::kmeans(pc$x[, seq_len(min(5L, ncol(pc$x))), drop = FALSE], centers = k, nstart = 20)$cluster
  co <- matrix(0, nrow(mat), nrow(mat), dimnames = list(rownames(mat), rownames(mat)))
  used <- 0L
  for (b in seq_len(n_bootstrap)) {
    cols <- sample(seq_len(ncol(mat)), ncol(mat), replace = TRUE)
    mm <- mat[, cols, drop = FALSE]
    res <- tryCatch({
      p <- stats::prcomp(mm, center = TRUE, scale. = TRUE)
      stats::kmeans(p$x[, seq_len(min(5L, ncol(p$x))), drop = FALSE], centers = k, nstart = 10)$cluster
    }, error = function(e) NULL)
    if (is.null(res)) next
    co <- co + outer(res, res, "==")
    used <- used + 1L
  }
  if (!used) return(data.frame(seed_id = rownames(mat), trajectory_cluster_stability = NA_real_, stringsAsFactors = FALSE))
  co <- co / used
  data.frame(
    seed_id = rownames(mat),
    trajectory_cluster = as.integer(base),
    trajectory_cluster_stability = vapply(seq_along(base), function(i) {
      same <- base == base[[i]]
      mean(co[i, same], na.rm = TRUE)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

o2pr_make_trajectory_labels <- function(curve_all, feat_long, status_all, cfg, trajectory_end_day, mid_window, late_window,
                                        collapse_terminal_offset = 8, min_late_drop = 8, stable_min_chr = 44) {
  if (!nrow(status_all)) return(data.frame())
  nmin <- as.numeric(cfg$N_MIN %||% 22)
  ids <- status_all$seed_id
  out <- data.frame(
    seed_id = ids,
    rule_regime = "ambiguous",
    trajectory_regime = "ambiguous",
    trajectory_drop_score = NA_real_,
    terminal_mean_N = NA_real_,
    rule_reason = "insufficient_2N_4N_support",
    stringsAsFactors = FALSE
  )
  for (seed in ids) {
    f2 <- feat_long[feat_long$seed_id == seed & feat_long$cohort == "2N", , drop = FALSE]
    f4 <- feat_long[feat_long$seed_id == seed & feat_long$cohort == "4N", , drop = FALSE]
    if (nrow(f2) != 1L || nrow(f4) != 1L) next
    idx <- match(seed, out$seed_id)
    terminal <- mean(c(f2$late_level, f4$late_level), na.rm = TRUE)
    drop <- mean(c(f2$late_drop_amplitude, f4$late_drop_amplitude), na.rm = TRUE)
    out$terminal_mean_N[idx] <- terminal
    out$trajectory_drop_score[idx] <- drop
    stable <- all(c(f2$late_level, f4$late_level) >= stable_min_chr, na.rm = FALSE) &&
      all(c(f2$late_drop_amplitude, f4$late_drop_amplitude) < min_late_drop, na.rm = FALSE)
    collapse <- all(c(f2$late_level, f4$late_level) <= nmin + collapse_terminal_offset, na.rm = FALSE) &&
      all(c(f2$late_drop_amplitude, f4$late_drop_amplitude) >= min_late_drop, na.rm = FALSE)
    if (isTRUE(stable)) {
      out$rule_regime[idx] <- "stable_high_chr"
      out$trajectory_regime[idx] <- "stable_high_chr"
      out$rule_reason[idx] <- "both_cohorts_late_high_and_no_large_late_drop"
    } else if (isTRUE(collapse)) {
      out$rule_regime[idx] <- "late_collapse_to_low_chr"
      out$trajectory_regime[idx] <- "late_collapse_to_low_chr"
      out$rule_reason[idx] <- "both_cohorts_late_near_Nmin_with_large_late_drop"
    } else {
      out$rule_reason[idx] <- "rule_thresholds_not_jointly_met"
    }
  }
  mat <- o2pr_common_grid_matrix(curve_all, trajectory_end_day)$matrix
  if (nrow(mat) > 0L && ncol(mat) > 0L) {
    sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
    mat <- mat[, is.finite(sds) & sds > 0, drop = FALSE]
  }
  if (nrow(mat) >= 3L) {
    pc <- tryCatch(stats::prcomp(mat, center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (!is.null(pc) && nrow(pc$x) >= 2L) {
      km <- stats::kmeans(pc$x[, seq_len(min(5L, ncol(pc$x))), drop = FALSE], centers = 2L, nstart = 50)
      cluster_terminal <- tapply(out$terminal_mean_N[match(rownames(mat), out$seed_id)], km$cluster, mean, na.rm = TRUE)
      collapse_cluster <- as.integer(names(cluster_terminal)[which.min(cluster_terminal)])
      cluster_label <- ifelse(km$cluster == collapse_cluster, "late_collapse_to_low_chr", "stable_high_chr")
      out$trajectory_cluster <- NA_integer_
      out$trajectory_cluster_label <- NA_character_
      out$trajectory_cluster[match(rownames(mat), out$seed_id)] <- km$cluster
      out$trajectory_cluster_label[match(rownames(mat), out$seed_id)] <- cluster_label
      amb <- out$trajectory_regime == "ambiguous" & !is.na(out$trajectory_cluster_label)
      out$trajectory_regime[amb] <- out$trajectory_cluster_label[amb]
      out$rule_reason[amb] <- "assigned_by_trajectory_functional_clustering"
    }
  }
  out
}

o2pr_parameter_outputs <- function(inputs, out_dir) {
  params_long <- inputs$params_long
  params_long$module <- vapply(params_long$parameter, o2ipa_param_module, character(1))
  params_long$transformed_value <- mapply(o2ipa_transform_parameter_value, params_long$parameter, params_long$value)
  o2pr_write_tsv(params_long, file.path(out_dir, "tables", "parameter_values_long.tsv"))

  raw <- o2ipa_params_wide(params_long, "value")
  raw_out <- data.frame(seed_id = rownames(raw), raw, check.names = FALSE)
  rownames(raw_out) <- NULL
  o2pr_write_tsv(raw_out, file.path(out_dir, "tables", "parameter_matrix_raw.tsv"))

  trans <- o2ipa_params_wide(params_long, "transformed_value")
  trans_out <- data.frame(seed_id = rownames(trans), trans, check.names = FALSE)
  rownames(trans_out) <- NULL
  o2pr_write_tsv(trans_out, file.path(out_dir, "tables", "parameter_matrix_transformed.tsv"))

  meta <- o2ipa_transform_metadata(list())
  o2pr_write_tsv(meta, file.path(out_dir, "tables", "parameter_transform_metadata.tsv"))
  list(raw = raw_out, transformed = trans_out, metadata = meta, long = params_long)
}

o2pr_scale_long_features_safe <- function(long_df, missing_feature_max = 0.5, missing_seed_max = 0.25) {
  out <- tryCatch(
    o2ipa_scale_long_features(long_df, missing_feature_max = missing_feature_max, missing_seed_max = missing_seed_max),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  if (!is.data.frame(long_df) || !nrow(long_df)) {
    return(list(long = data.frame(), wide = data.frame(), metadata = data.frame(), missingness = data.frame(), excluded_seeds = character()))
  }
  long_df$transformed_value <- o2ipa_feature_transform(long_df$raw_value, long_df$feature_type)
  long_df$feature_key <- paste(long_df$fingerprint_scope, long_df$module, long_df$feature, sep = "||")
  meta <- unique(long_df[, c("feature_key", "fingerprint_scope", "module", "feature", "feature_type")])
  meta$missing_fraction <- vapply(meta$feature_key, function(k) mean(is.na(long_df$transformed_value[long_df$feature_key == k])), numeric(1))
  meta$center <- vapply(meta$feature_key, function(k) stats::median(long_df$transformed_value[long_df$feature_key == k], na.rm = TRUE), numeric(1))
  meta$scale <- NA_real_
  meta$zero_variance <- TRUE
  meta$retained_for_clustering <- FALSE
  long_df$center <- meta$center[match(long_df$feature_key, meta$feature_key)]
  long_df$scale <- NA_real_
  long_df$retained_for_clustering <- FALSE
  long_df$scaled_value <- NA_real_
  long_df$module_weight <- NA_real_
  long_df$clustering_value <- NA_real_
  missing <- data.frame(seed_id = unique(long_df$seed_id), missing_fraction_retained_features = NA_real_, stringsAsFactors = FALSE)
  list(long = long_df, wide = data.frame(), metadata = meta, missingness = missing, excluded_seeds = character())
}

o2pr_scale_and_write_process <- function(proc, out_dir) {
  scaled <- list()
  for (nm in names(proc)[names(proc) %in% c("intrinsic_full", "intrinsic_18", "exposure_full", "exposure_18")]) {
    sc <- o2pr_scale_long_features_safe(proc[[nm]], missing_feature_max = 0.5, missing_seed_max = 0.25)
    scaled[[nm]] <- sc
    outbase <- switch(nm,
      intrinsic_full = "process_intrinsic_full",
      intrinsic_18 = "process_intrinsic_18only",
      exposure_full = "process_exposure_full",
      exposure_18 = "process_exposure_18only"
    )
    o2pr_write_tsv(sc$wide, file.path(out_dir, "tables", paste0(outbase, "_scaled.tsv")))
    o2pr_write_tsv(sc$long, file.path(out_dir, "tables", paste0(outbase, "_scaled_long.tsv")))
    if (identical(nm, "intrinsic_full")) {
      o2pr_write_tsv(sc$wide, file.path(out_dir, "tables", "process_intrinsic_scaled.tsv"))
      o2pr_write_tsv(o2ipa_long_to_wide(transform(proc[[nm]], feature_key = paste(module, feature, sep = "||")), "feature_key", "raw_value"), file.path(out_dir, "tables", "process_intrinsic_wide_raw.tsv"))
    }
    if (identical(nm, "exposure_full")) {
      o2pr_write_tsv(sc$wide, file.path(out_dir, "tables", "process_exposure_scaled.tsv"))
      o2pr_write_tsv(o2ipa_long_to_wide(transform(proc[[nm]], feature_key = paste(module, feature, sep = "||")), "feature_key", "raw_value"), file.path(out_dir, "tables", "process_exposure_wide_raw.tsv"))
    }
  }
  for (scope in c("full", "18only")) {
    il <- proc[[paste0("intrinsic_", if (scope == "full") "full" else "18")]]
    el <- proc[[paste0("exposure_", if (scope == "full") "full" else "18")]]
    comb <- rbind(il, el)
    comb$fingerprint_scope <- paste0("combined_", scope)
    sc <- o2pr_scale_long_features_safe(comb, missing_feature_max = 0.5, missing_seed_max = 0.25)
    scaled[[paste0("combined_", scope)]] <- sc
    o2pr_write_tsv(sc$wide, file.path(out_dir, "tables", paste0("process_combined_", scope, "_scaled.tsv")))
    o2pr_write_tsv(sc$long, file.path(out_dir, "tables", paste0("process_combined_", scope, "_scaled_long.tsv")))
  }
  meta <- do.call(rbind, lapply(names(scaled), function(nm) {
    x <- scaled[[nm]]$metadata
    x$matrix <- nm
    x
  }))
  missing <- do.call(rbind, lapply(names(scaled), function(nm) {
    x <- scaled[[nm]]$missingness
    x$matrix <- nm
    x
  }))
  zero <- meta[meta$zero_variance %in% TRUE, , drop = FALSE]
  module_counts <- aggregate(meta$retained_for_clustering %in% TRUE, by = list(matrix = meta$matrix, module = meta$module), FUN = sum)
  names(module_counts)[3] <- "n_retained_features"
  module_counts$module_weight <- ifelse(module_counts$n_retained_features > 0, 1 / sqrt(module_counts$n_retained_features), NA_real_)
  o2pr_write_tsv(meta, file.path(out_dir, "tables", "process_feature_metadata.tsv"))
  o2pr_write_tsv(missing, file.path(out_dir, "tables", "process_feature_missingness.tsv"))
  o2pr_write_tsv(zero, file.path(out_dir, "tables", "process_feature_zero_variance.tsv"))
  o2pr_write_tsv(module_counts, file.path(out_dir, "tables", "module_weights.tsv"))
  scaled
}

o2pr_rbind_fill <- function(xs) {
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(cols, names(x))
    for (m in missing) x[[m]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, xs)
}

o2pr_cluster_process <- function(scaled, out_dir, n_bootstrap, random_seed) {
  sources <- c("intrinsic_full", "intrinsic_18", "exposure_full", "exposure_18", "combined_full", "combined_18only")
  diag <- list()
  memb <- list()
  med <- list()
  stab <- list()
  consensus_written <- FALSE
  clusters <- list()
  for (src in sources) {
    if (is.null(scaled[[src]])) next
    d <- o2ipa_distance(scaled[[src]]$wide)
    o2ipa_write_distance(d, paste0("process_", src), out_dir)
    cl <- o2ipa_cluster_distance(d, n_bootstrap = n_bootstrap, random_seed = random_seed, feature_wide = scaled[[src]]$wide, feature_meta = scaled[[src]]$metadata)
    clusters[[src]] <- cl
    dd <- cl$diagnostics
    if (nrow(dd)) dd$cluster_source <- src
    mm <- cl$membership
    if (nrow(mm)) mm$cluster_source <- src
    md <- cl$medoids
    if (nrow(md)) md$cluster_source <- src
    ss <- cl$stability
    if (nrow(ss)) ss$cluster_source <- src
    diag[[src]] <- dd
    memb[[src]] <- mm
    med[[src]] <- md
    stab[[src]] <- ss
    if (!consensus_written && !is.null(cl$consensus)) {
      o2pr_write_tsv(data.frame(seed_id = rownames(cl$consensus), cl$consensus, check.names = FALSE), file.path(out_dir, "tables", "process_cluster_consensus_matrix.tsv"))
      consensus_written <- TRUE
    }
  }
  o2pr_write_tsv(o2pr_rbind_fill(diag), file.path(out_dir, "tables", "process_cluster_k_diagnostics.tsv"))
  o2pr_write_tsv(o2pr_rbind_fill(memb), file.path(out_dir, "tables", "process_cluster_membership.tsv"))
  o2pr_write_tsv(o2pr_rbind_fill(med), file.path(out_dir, "tables", "process_cluster_medoids.tsv"))
  o2pr_write_tsv(o2pr_rbind_fill(stab), file.path(out_dir, "tables", "process_cluster_stability.tsv"))
  if (!consensus_written) o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "process_cluster_consensus_matrix.tsv"))
  clusters
}

o2pr_adjusted_rand <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- as.factor(a[ok])
  b <- as.factor(b[ok])
  n <- length(a)
  if (n < 2L) return(NA_real_)
  tab <- table(a, b)
  comb2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(comb2(tab))
  sum_a <- sum(comb2(rowSums(tab)))
  sum_b <- sum(comb2(colSums(tab)))
  total <- comb2(n)
  expected <- sum_a * sum_b / total
  maxv <- 0.5 * (sum_a + sum_b)
  if (abs(maxv - expected) < 1e-12) return(NA_real_)
  (sum_ij - expected) / (maxv - expected)
}

o2pr_nmi <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- as.factor(a[ok])
  b <- as.factor(b[ok])
  if (length(a) < 2L) return(NA_real_)
  tab <- table(a, b)
  p <- tab / sum(tab)
  pa <- rowSums(p)
  pb <- colSums(p)
  mi <- 0
  for (i in seq_len(nrow(p))) for (j in seq_len(ncol(p))) {
    if (p[i, j] > 0) mi <- mi + p[i, j] * log(p[i, j] / (pa[i] * pb[j]))
  }
  ha <- -sum(pa[pa > 0] * log(pa[pa > 0]))
  hb <- -sum(pb[pb > 0] * log(pb[pb > 0]))
  if (ha <= 0 || hb <= 0) return(NA_real_)
  mi / sqrt(ha * hb)
}

o2pr_concordance_outputs <- function(labels, clusters, scaled, out_dir, n_permutations, random_seed) {
  labels2 <- labels[labels$trajectory_regime %in% c("stable_high_chr", "late_collapse_to_low_chr"), , drop = FALSE]
  rows <- list()
  cont_rows <- list()
  set.seed(random_seed)
  for (src in names(clusters)) {
    cl <- clusters[[src]]
    if (is.null(cl) || is.na(cl$recommended_k) || !nrow(cl$membership)) next
    memb <- cl$membership[cl$membership$k == cl$recommended_k, , drop = FALSE]
    d <- merge(labels2, memb, by = "seed_id", all = FALSE)
    if (nrow(d) < 3L) next
    ari <- o2pr_adjusted_rand(d$trajectory_regime, d$cluster)
    nmi <- o2pr_nmi(d$trajectory_regime, d$cluster)
    p_perm <- NA_real_
    if (n_permutations > 0L && is.finite(ari)) {
      null <- replicate(n_permutations, o2pr_adjusted_rand(sample(d$trajectory_regime), d$cluster))
      p_perm <- mean(abs(null) >= abs(ari), na.rm = TRUE)
    }
    tab <- table(d$trajectory_regime, d$cluster)
    p_test <- tryCatch({
      if (all(dim(tab) == c(2, 2))) stats::fisher.test(tab)$p.value else suppressWarnings(stats::chisq.test(tab)$p.value)
    }, error = function(e) NA_real_)
    rows[[src]] <- data.frame(cluster_source = src, n_seed = nrow(d), k = cl$recommended_k, adjusted_rand_index = ari, normalized_mutual_information = nmi, contingency_p_value = p_test, permutation_p_value = p_perm, stringsAsFactors = FALSE)
    ct <- as.data.frame(tab)
    names(ct) <- c("trajectory_regime", "process_cluster", "n")
    ct$cluster_source <- src
    cont_rows[[src]] <- ct
  }
  o2pr_write_tsv(do.call(rbind, rows), file.path(out_dir, "tables", "trajectory_process_concordance.tsv"))
  o2pr_write_tsv(do.call(rbind, cont_rows), file.path(out_dir, "tables", "trajectory_process_contingency.tsv"))
  o2pr_feature_differences(labels2, scaled$combined_full$long, out_dir)
  invisible(rows)
}

o2pr_feature_differences <- function(labels, scaled_long, out_dir) {
  if (is.null(scaled_long) || !nrow(scaled_long) || !nrow(labels)) {
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "process_regime_feature_differences.tsv"))
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "process_regime_module_differences.tsv"))
    return(invisible(NULL))
  }
  d <- merge(labels[, c("seed_id", "trajectory_regime")], scaled_long, by = "seed_id")
  feats <- split(d, d$feature_key)
  rows <- lapply(feats, function(x) {
    g <- split(x$scaled_value, x$trajectory_regime)
    a <- g[["stable_high_chr"]]
    b <- g[["late_collapse_to_low_chr"]]
    p <- if (length(a) >= 2L && length(b) >= 2L) tryCatch(stats::wilcox.test(a, b)$p.value, error = function(e) NA_real_) else NA_real_
    data.frame(
      feature_key = x$feature_key[[1]],
      module = x$module[[1]],
      feature = x$feature[[1]],
      median_stable_high_chr = stats::median(a, na.rm = TRUE),
      median_late_collapse_to_low_chr = stats::median(b, na.rm = TRUE),
      standardized_effect_size = stats::median(b, na.rm = TRUE) - stats::median(a, na.rm = TRUE),
      rank_p_value = p,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$BH_adjusted_p_value <- stats::p.adjust(out$rank_p_value, method = "BH")
  o2pr_write_tsv(out[order(out$BH_adjusted_p_value), , drop = FALSE], file.path(out_dir, "tables", "process_regime_feature_differences.tsv"))
  mod <- aggregate(out$standardized_effect_size, by = list(module = out$module), FUN = function(z) stats::median(abs(z), na.rm = TRUE))
  names(mod)[2] <- "median_abs_standardized_effect_size"
  o2pr_write_tsv(mod, file.path(out_dir, "tables", "process_regime_module_differences.tsv"))
}

o2pr_boundary_diagnostics <- function(traj, cfg, out_dir) {
  density <- traj$density
  nmin <- as.numeric(cfg$N_MIN %||% 22)
  if (!is.data.frame(density) || !nrow(density)) {
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "boundary_diagnostics_by_seed.tsv"))
    o2pr_write_tsv(data.frame(test = c("expanded_lower_state_grid", "alternative_boundary_handling", "smaller_simulation_step", "stricter_integration_tolerance"), reason = "not_available_without_safe_refit_or_simulation_api"), file.path(out_dir, "tables", "unavailable_boundary_tests.tsv"))
    return(invisible(NULL))
  }
  density$day <- as.numeric(density$day)
  density$N <- as.numeric(density$N)
  density$fraction <- as.numeric(density$fraction)
  rows <- lapply(split(density, density$seed_id), function(d) {
    dd <- d[d$N == nmin, , drop = FALSE]
    by_day <- if (nrow(dd)) {
      aggregate(dd$fraction, by = list(day = dd$day), FUN = sum, na.rm = TRUE)
    } else {
      data.frame(day = numeric(), x = numeric())
    }
    max_frac <- if (nrow(by_day)) max(by_day$x, na.rm = TRUE) else NA_real_
    term_day <- if (nrow(by_day)) max(by_day$day, na.rm = TRUE) else NA_real_
    term_frac <- if (nrow(by_day)) by_day$x[which.max(by_day$day)] else NA_real_
    time_first <- if (nrow(by_day) && any(by_day$x > 0.05, na.rm = TRUE)) min(by_day$day[by_day$x > 0.05], na.rm = TRUE) else NA_real_
    data.frame(
      seed_id = d$seed_id[[1]],
      terminal_fraction_at_Nmin = term_frac,
      maximum_fraction_at_Nmin = max_frac,
      time_first_accumulating_at_Nmin = time_first,
      terminal_distance_to_Nmin = NA_real_,
      transition_mass_dropped_at_Nmin = NA_real_,
      whether_terminal_mean_is_within_2_of_Nmin = NA,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (nrow(traj$labels)) {
    out <- merge(out, traj$labels[, c("seed_id", "terminal_mean_N", "trajectory_regime")], by = "seed_id", all.x = TRUE)
    out$terminal_distance_to_Nmin <- out$terminal_mean_N - nmin
    out$whether_terminal_mean_is_within_2_of_Nmin <- out$terminal_distance_to_Nmin <= 2
  }
  o2pr_write_tsv(out, file.path(out_dir, "tables", "boundary_diagnostics_by_seed.tsv"))
  enrich <- aggregate(out$maximum_fraction_at_Nmin, by = list(trajectory_regime = out$trajectory_regime), FUN = mean, na.rm = TRUE)
  names(enrich)[2] <- "mean_maximum_fraction_at_Nmin"
  o2pr_write_tsv(enrich, file.path(out_dir, "tables", "boundary_regime_enrichment.tsv"))
  o2pr_write_tsv(data.frame(test = c("expanded_lower_state_grid", "alternative_boundary_handling", "smaller_simulation_step", "stricter_integration_tolerance"), reason = "not_run_by_core_pipeline; requires explicit diagnostic simulation API and independent out_dir"), file.path(out_dir, "tables", "unavailable_boundary_tests.tsv"))
  o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "boundary_sensitivity_results.tsv"))
}

o2pr_oxygen_switch_diagnostics <- function(proc, out_dir) {
  exp <- proc$exposure_full
  keep_features <- c(
    "minimum_realized_O2",
    "terminal_O2",
    "O2_AUC",
    "switch_available",
    "switch_margin",
    "switch_crossed",
    "ever_below_O2_2pct",
    "ever_below_O2_1pct",
    "ever_below_O2_0p5pct",
    "time_below_O2_2pct",
    "time_below_O2_1pct",
    "time_below_O2_0p5pct"
  )
  if (!is.data.frame(exp) || !nrow(exp)) {
    o2pr_write_tsv(data.frame(), file.path(out_dir, "tables", "oxygen_switch_diagnostics.tsv"))
    return(invisible(NULL))
  }
  out <- exp[exp$module == "oxygen_feedback" & exp$feature %in% keep_features, , drop = FALSE]
  o2pr_write_tsv(out, file.path(out_dir, "tables", "oxygen_switch_diagnostics.tsv"))
}

o2pr_placeholder_advanced <- function(out_dir,
                                      analysis_level = "core",
                                      module_factorial = FALSE,
                                      n_local_perturbations = 0L,
                                      local_perturbation_sd = NA_real_) {
  context <- paste0(
    "analysis_level=", analysis_level,
    "; module_factorial=", module_factorial,
    "; n_local_perturbations=", n_local_perturbations,
    "; local_perturbation_sd=", local_perturbation_sd
  )
  o2pr_write_tsv(data.frame(status = "not_run_in_current_pipeline", reason = "objective-compatible interpolation requires safe objective evaluator and actual trajectory simulator wiring", context = context), file.path(out_dir, "tables", "parameter_interpolation_path.tsv"))
  o2pr_write_tsv(data.frame(status = "unavailable", reason = "objective evaluator not called by current ploidy regime pipeline", context = context), file.path(out_dir, "tables", "interpolation_objective.tsv"))
  o2pr_write_tsv(data.frame(status = "unavailable", reason = "interpolation process path not run by current pipeline", context = context), file.path(out_dir, "tables", "interpolation_process_features.tsv"))
  o2pr_write_tsv(data.frame(status = "unavailable", reason = "actual 1000-day trajectory interpolation not run by current pipeline", context = context), file.path(out_dir, "tables", "interpolation_trajectories.tsv"))
  o2pr_write_tsv(data.frame(status = "not_run_in_current_pipeline", reason = "module swaps require actual trajectory re-simulation to claim causal mode switch", context = context), file.path(out_dir, "tables", "counterfactual_module_swaps.tsv"))
  o2pr_write_tsv(data.frame(status = "unavailable", reason = "counterfactual trajectory objective not recalculated", context = context), file.path(out_dir, "tables", "counterfactual_module_effects.tsv"))
  o2pr_write_tsv(data.frame(status = if (module_factorial) "requested_but_unavailable" else "disabled", reason = "factorial module swaps require safe hybrid trajectory simulation before producing effects", context = context), file.path(out_dir, "tables", "counterfactual_factorial_results.tsv"))
  o2pr_write_tsv(data.frame(status = if (n_local_perturbations > 0L) "requested_but_unavailable" else "not_run_in_current_pipeline", reason = "local perturbations require objective evaluator for near-optimal claims", context = context), file.path(out_dir, "tables", "local_perturbation_results.tsv"))
  o2pr_write_tsv(data.frame(status = "unavailable", reason = "near-optimal perturbation subset not defined without objective recomputation", context = context), file.path(out_dir, "tables", "local_mode_stability.tsv"))
  o2pr_write_tsv(data.frame(status = "unavailable", reason = "repeated stratified process-PC prediction CV is not implemented in the current pipeline", context = context), file.path(out_dir, "tables", "process_mode_prediction_cv.tsv"))
  o2pr_write_tsv(data.frame(status = "unavailable", reason = "leave-one-module-out prediction requires the CV predictor, which is not implemented in the current pipeline", context = context), file.path(out_dir, "tables", "leave_one_module_out_results.tsv"))
}

o2pr_evidence_summary <- function(out_dir, manifest, labels, clusters, objective_deltas = c(2, 5, 10)) {
  cl <- clusters$combined_full
  stable_cluster <- !is.null(cl) && !is.na(cl$recommended_k)
  rows <- lapply(objective_deltas, function(cutoff) {
    keep <- manifest$seed_id[manifest$fit_success %in% TRUE & manifest$delta_objective <= cutoff]
    lab <- labels[labels$seed_id %in% keep, , drop = FALSE]
    n_stable <- sum(lab$trajectory_regime == "stable_high_chr", na.rm = TRUE)
    n_collapse <- sum(lab$trajectory_regime == "late_collapse_to_low_chr", na.rm = TRUE)
    classification <- if (n_stable > 0 && n_collapse > 0 && stable_cluster) {
      "unresolved"
    } else if (n_stable > 0 && n_collapse > 0) {
      "long_horizon_prediction_nonidentifiability"
    } else {
      "unresolved"
    }
    data.frame(
      objective_cutoff = cutoff,
      n_valid_seed = length(keep),
      n_stable_high_chr = n_stable,
      n_late_collapse_to_low_chr = n_collapse,
      n_ambiguous = sum(lab$trajectory_regime == "ambiguous", na.rm = TRUE),
      evidence_classification = classification,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
