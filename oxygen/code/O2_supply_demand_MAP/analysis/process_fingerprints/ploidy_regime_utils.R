#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("Matrix", quietly = TRUE)) library(Matrix)
})

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

o2pr_source_process_utils <- function(script_dir) {
  source(file.path(script_dir, "process_fingerprint_utils.R"), local = TRUE)
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
    file.path(out_dir, c("tables", "figures", "cache", "logs", "report")),
    dir.create,
    logical(1),
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

o2pr_git_value <- function(args, repo_root) {
  out <- tryCatch(
    system2("git", c("-C", repo_root, args), stdout = TRUE, stderr = TRUE),
    error = function(e) NA_character_
  )
  if (!length(out)) NA_character_ else paste(out, collapse = "\n")
}

o2pr_first_seed_cfg <- function(manifest) {
  for (p in manifest$config_file) {
    if (!is.na(p) && file.exists(p)) {
      cfg <- tryCatch(readRDS(p), error = function(e) NULL)
      if (!is.null(cfg)) return(cfg)
    }
  }
  list(N_MIN = 22L, N_MAX = 154L, N_UNIT = 22L, DT = 0.05, start_with = "chr_number")
}

o2pr_cfg_metadata <- function(cfg) {
  data.frame(
    metric = c(
      "Nmin", "Nmax", "N_UNIT", "DT", "start_with", "boundary_default",
      "o2_burden_feedback", "O2_growth", "ploidy_O2_death",
      "o2_S0_upper_bound", "o2_Nref", "trajectory_value_semantics"
    ),
    value = c(
      cfg$N_MIN %||% NA, cfg$N_MAX %||% NA, cfg$N_UNIT %||% NA,
      cfg$DT %||% NA, cfg$start_with %||% NA, cfg$boundary %||% "drop",
      cfg$o2_burden_feedback %||% NA, cfg$O2_growth %||% NA,
      cfg$ploidy_O2_death %||% NA, cfg$o2_S0_upper_bound %||% NA,
      cfg$o2_Nref %||% NA,
      if (identical(as.character(cfg$start_with %||% "ploidy"), "chr_number")) {
        "weighted mean chromosome number N"
      } else {
        "weighted mean ploidy converted to N when needed"
      }
    ),
    stringsAsFactors = FALSE
  )
}

o2pr_generate_viz_if_requested <- function(seed_dir, script_dir, generate_viz = FALSE) {
  if (!isTRUE(generate_viz)) return(invisible(FALSE))
  viz_script <- file.path(dirname(script_dir), "vis", "viz_invivo_model_O2_supply_demand_MAP_results.R")
  if (!file.exists(viz_script)) return(invisible(FALSE))
  status <- tryCatch(
    system2("Rscript", c(viz_script, "--fit_dir", seed_dir, "--report_dt", "1"), stdout = TRUE, stderr = TRUE),
    error = function(e) conditionMessage(e)
  )
  invisible(status)
}

o2pr_find_trajectory_path <- function(seed_dir) {
  candidates <- file.path(seed_dir, "viz", c(
    "predict_ploidy_weighted_mean_0_1000day.tsv",
    "ploidy_weighted_mean_timecourse.tsv",
    "predict_ploidy_weighted_mean_0_300day.tsv",
    "predict_ploidy_0_1000day.tsv",
    "ploidy_timecourse.tsv",
    "predict_ploidy_0_300day.tsv"
  ))
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1]] else NA_character_
}

o2pr_value_col_for_trajectory <- function(tab) {
  for (nm in c("weighted_mean_N", "weighted_mean_endpoint", "ploidy_value", "weighted_mean_ploidy")) {
    if (nm %in% names(tab)) return(nm)
  }
  NA_character_
}

o2pr_weighted_from_density <- function(tab) {
  needed <- c("cohort", "day", "N", "fraction")
  if (!all(needed %in% names(tab))) return(NULL)
  tab$day <- suppressWarnings(as.numeric(tab$day))
  tab$N <- suppressWarnings(as.numeric(tab$N))
  tab$fraction <- suppressWarnings(as.numeric(tab$fraction))
  tab <- tab[is.finite(tab$day) & is.finite(tab$N) & is.finite(tab$fraction), , drop = FALSE]
  if (!nrow(tab)) return(NULL)
  key <- paste(tab$cohort, tab$day, sep = "\r")
  rows <- lapply(split(seq_len(nrow(tab)), key), function(idx) {
    d <- tab[idx, , drop = FALSE]
    denom <- sum(d$fraction, na.rm = TRUE)
    data.frame(
      cohort = as.character(d$cohort[[1]]),
      day = as.numeric(d$day[[1]]),
      value = if (is.finite(denom) && denom > 0) sum(d$N * d$fraction, na.rm = TRUE) / denom else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

o2pr_read_trajectory_seed <- function(seed_id, seed_dir, cfg, script_dir, generate_viz = FALSE) {
  if (is.na(seed_dir) || !nzchar(seed_dir) || !dir.exists(seed_dir)) {
    return(list(curve = data.frame(), status = "missing_seed_dir", source_file = NA_character_))
  }
  path <- o2pr_find_trajectory_path(seed_dir)
  if (is.na(path)) {
    o2pr_generate_viz_if_requested(seed_dir, script_dir, generate_viz)
    path <- o2pr_find_trajectory_path(seed_dir)
  }
  if (is.na(path)) {
    return(list(curve = data.frame(), status = "missing_trajectory_table", source_file = NA_character_))
  }
  tab <- tryCatch(o2pr_read_tsv(path), error = function(e) NULL)
  if (is.null(tab) || !nrow(tab) || !"cohort" %in% names(tab) || !"day" %in% names(tab)) {
    return(list(curve = data.frame(), status = "bad_trajectory_table", source_file = path))
  }
  vcol <- o2pr_value_col_for_trajectory(tab)
  if (is.na(vcol)) {
    dens <- o2pr_weighted_from_density(tab)
    if (is.null(dens)) return(list(curve = data.frame(), status = "no_weighted_value_column", source_file = path))
    tab2 <- dens
  } else {
    tab2 <- data.frame(
      cohort = as.character(tab$cohort),
      day = suppressWarnings(as.numeric(tab$day)),
      value = suppressWarnings(as.numeric(tab[[vcol]])),
      stringsAsFactors = FALSE
    )
    if (identical(vcol, "weighted_mean_ploidy") && !identical(as.character(cfg$start_with %||% "ploidy"), "chr_number")) {
      tab2$value <- tab2$value * as.numeric(cfg$N_UNIT %||% 22)
    }
  }
  tab2 <- tab2[tab2$cohort %in% c("2N", "4N") & is.finite(tab2$day) & is.finite(tab2$value), , drop = FALSE]
  if (!nrow(tab2)) return(list(curve = data.frame(), status = "no_valid_2N_4N_rows", source_file = path))
  agg <- aggregate(tab2$value, by = list(cohort = tab2$cohort, day = tab2$day), FUN = mean, na.rm = TRUE)
  names(agg)[3] <- "trajectory_value"
  agg$seed_id <- seed_id
  agg$source_file <- normalizePath(path, mustWork = FALSE)
  agg <- agg[, c("seed_id", "cohort", "day", "trajectory_value", "source_file")]
  list(curve = agg[order(agg$cohort, agg$day), , drop = FALSE], status = "ok", source_file = path)
}

o2pr_read_density_seed <- function(seed_id, seed_dir) {
  if (is.na(seed_dir) || !nzchar(seed_dir)) return(data.frame())
  path <- file.path(seed_dir, "viz", "ploidy_timecourse.tsv")
  if (!file.exists(path)) path <- file.path(seed_dir, "viz", "predict_ploidy_0_1000day.tsv")
  if (!file.exists(path)) return(data.frame())
  tab <- tryCatch(o2pr_read_tsv(path), error = function(e) data.frame())
  needed <- c("cohort", "day", "N", "fraction")
  if (!all(needed %in% names(tab))) return(data.frame())
  tab$seed_id <- seed_id
  tab[, c("seed_id", needed), drop = FALSE]
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

o2pr_collect_trajectories <- function(inputs, script_dir, trajectory_end_day, mid_window, late_window, generate_viz = FALSE, cfg = NULL) {
  cfg <- cfg %||% o2pr_first_seed_cfg(inputs$manifest)
  curves <- list()
  status <- list()
  density <- list()
  manifest <- inputs$manifest
  for (i in seq_len(nrow(manifest))) {
    seed <- manifest$seed_id[[i]]
    seed_dir <- manifest$seed_dir[[i]]
    tr <- o2pr_read_trajectory_seed(seed, seed_dir, cfg, script_dir, generate_viz)
    if (nrow(tr$curve)) curves[[seed]] <- tr$curve
    density[[seed]] <- o2pr_read_density_seed(seed, seed_dir)
    status[[seed]] <- data.frame(
      seed_id = seed,
      trajectory_available = identical(tr$status, "ok"),
      trajectory_status = tr$status,
      trajectory_source_file = tr$source_file,
      stringsAsFactors = FALSE
    )
  }
  curve_all <- if (length(curves)) do.call(rbind, curves) else data.frame()
  status_all <- do.call(rbind, status)
  density_all <- do.call(rbind, density)
  nmin <- as.numeric(cfg$N_MIN %||% 22)
  long_features <- list()
  if (nrow(curve_all)) {
    for (key in split(seq_len(nrow(curve_all)), paste(curve_all$seed_id, curve_all$cohort, sep = "||"))) {
      d <- curve_all[key, , drop = FALSE]
      if (!nrow(d)) next
      long_features[[length(long_features) + 1L]] <- o2pr_cohort_features(
        d$seed_id[[1]], d$cohort[[1]], d, mid_window, late_window, nmin
      )
    }
  }
  feat_long <- if (length(long_features)) do.call(rbind, long_features) else data.frame()
  labels <- o2pr_make_trajectory_labels(curve_all, feat_long, status_all, cfg, trajectory_end_day, mid_window, late_window)
  wide <- o2pr_features_wide(feat_long, labels)
  list(curves = curve_all, density = density_all, status = status_all, features_long = feat_long, features_wide = wide, labels = labels)
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

o2pr_build_G <- function(model_env, cfg, run_params, O2) {
  fn <- get("cpp_o2simps_build_G_for_o2_triplet", envir = model_env, inherits = TRUE)
  tri <- fn(
    O2 = as.numeric(O2),
    O2_crit = as.numeric(run_params$O2_crit %||% cfg$o2_crit_init %||% 1),
    N0min = as.integer(cfg$N_MIN %||% 22L),
    N0max = as.integer(cfg$N_MAX %||% 154L),
    N1min = as.integer(cfg$N_MIN %||% 22L),
    N1max = as.integer(cfg$N_MAX %||% 154L),
    lam_max = as.numeric(run_params$lam_max),
    p_mis_base = as.numeric(run_params$p_mis_base %||% cfg$p_mis_base %||% 1e-5),
    p_misseg = as.numeric(run_params$p_misseg %||% 0),
    k_o_mis = as.numeric(run_params$k_o_mis %||% 50),
    p_wgd = as.numeric(run_params$p_wgd %||% 0),
    boundary = as.character(run_params$boundary %||% "drop"),
    eps_tail = 1e-8,
    buffer_smax = as.numeric(run_params$buffer_smax %||% 1),
    buffer_beta = as.numeric(run_params$buffer_beta %||% 0),
    buffer_n_exp = as.numeric(run_params$buffer_n_exp %||% 1),
    N_unit = as.integer(cfg$N_UNIT %||% 22L),
    beta_size = 0,
    O2_growth = isTRUE(run_params$O2_growth %||% cfg$O2_growth %||% TRUE),
    alpha_o2 = as.numeric(run_params$alpha_o2 %||% 0),
    gamma_growth = as.numeric(run_params$gamma_growth %||% 1),
    mu_hp = as.numeric(run_params$mu_hp %||% 0),
    gamma_mu = as.numeric(run_params$gamma_mu %||% 1),
    n_O = as.numeric(run_params$n_O %||% 1),
    ploidy_O2_death = as.character(run_params$ploidy_O2_death %||% cfg$ploidy_O2_death %||% "diploid_NULL")
  )
  G <- Matrix::sparseMatrix(
    i = as.integer(tri$i),
    j = as.integer(tri$j),
    x = as.numeric(tri$x),
    dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
    repr = "C"
  )
  attr(G, "triplet") <- tri
  G
}

o2pr_run_params_from_vec <- function(vec, cfg) {
  rp <- as.list(vec)
  rp$o2_min <- if ("o2_min" %in% names(vec)) vec[["o2_min"]] else (cfg$o2_min %||% 0)
  rp$o2_S0_upper_bound <- cfg$o2_S0_upper_bound %||% 5
  rp$o2_Nref <- cfg$o2_Nref %||% 1e6
  rp$O2_growth <- cfg$O2_growth %||% TRUE
  rp$ploidy_O2_death <- cfg$ploidy_O2_death %||% "ploidy_related"
  rp$boundary <- cfg$boundary %||% "drop"
  rp
}

o2pr_build_process_one <- function(seed_id, run_params, model_env, cfg, o2_grid, scope = "full") {
  rows <- list()
  add <- function(module, feature, value, type = "positive", source = "model_probe") {
    rows[[length(rows) + 1L]] <<- o2ipa_feature_row(seed_id, scope, module, feature, value, type, source)
  }
  states <- as.integer(pmin(pmax(c(22, 44, 66, 88), cfg$N_MIN %||% 22), cfg$N_MAX %||% 154))
  states <- sort(unique(states))
  state_note <- data.frame(requested_state = c(22, 44, 66, 88), used_state = as.integer(pmin(pmax(c(22, 44, 66, 88), cfg$N_MIN %||% 22), cfg$N_MAX %||% 154)))
  ngrid <- seq.int(as.integer(cfg$N_MIN %||% 22), as.integer(cfg$N_MAX %||% 154))
  for (O2 in o2_grid) {
    lam <- o2ipa_call_model(model_env, ".lambda_eff_of_O2", O2 = rep(O2, length(states)), run_params = run_params, N = states, O2_growth = isTRUE(run_params$O2_growth %||% TRUE))
    mu <- o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(states)), run_params = run_params, N = states)
    pm <- o2ipa_call_model(model_env, ".pmisseg_of_O2", O2 = rep(O2, length(states)), run_params = run_params, N = states)
    stress <- o2ipa_call_model(model_env, ".resource_stress_of_O2", O2 = O2, run_params = run_params)
    for (i in seq_along(states)) {
      suffix <- paste0("_O2_", gsub("[^0-9]+", "p", O2), "_N_", states[[i]])
      add("hypoxia_sensing", paste0("hypoxia_stress", suffix), stress, "fraction")
      add("growth_selection", paste0("lambda", suffix), lam[[i]], "positive")
      add("death_selection", paste0("mu", suffix), mu[[i]], "positive")
      add("CIN", paste0("p_misseg", suffix), pm[[i]], "probability")
      add("CIN", paste0("missegregation_flux", suffix), lam[[i]] * pm[[i]], "positive")
      add("WGD", paste0("WGD_flux", suffix), lam[[i]] * as.numeric(run_params$p_wgd %||% 0), "positive")
    }
    G <- o2pr_build_G(model_env, cfg, run_params, O2)
    mu_all <- as.numeric(o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = rep(O2, length(ngrid)), run_params = run_params, N = ngrid))
    M <- G - Matrix::Diagonal(x = mu_all)
    eff <- vapply(states, function(N) {
      col <- as.integer(N - (cfg$N_MIN %||% 22) + 1L)
      if (col < 1L || col > ncol(G)) return(NA_real_)
      sum(G[, col]) - mu_all[[col]]
    }, numeric(1))
    names(eff) <- states
    if (all(c("22", "44") %in% names(eff))) add("attractor", paste0("selection_22_vs_44_O2_", gsub("[^0-9]+", "p", O2)), eff[["22"]] - eff[["44"]], "signed")
    if (all(c("44", "88") %in% names(eff))) add("attractor", paste0("selection_44_vs_88_O2_", gsub("[^0-9]+", "p", O2)), eff[["44"]] - eff[["88"]], "signed")
    eig <- tryCatch(eigen(as.matrix(M), only.values = FALSE), error = function(e) NULL)
    if (!is.null(eig)) {
      idx <- which.max(Re(eig$values))
      v <- Re(eig$vectors[, idx])
      if (sum(v, na.rm = TRUE) < 0) v <- -v
      nonneg <- all(v >= -1e-8, na.rm = TRUE)
      v <- pmax(v, 0)
      if (sum(v) > 0) v <- v / sum(v)
      lambda1 <- Re(eig$values[[idx]])
      lambda2 <- sort(Re(eig$values), decreasing = TRUE)[min(2L, length(eig$values))]
      add("attractor", paste0("dominant_mean_N_O2_", gsub("[^0-9]+", "p", O2)), sum(ngrid * v), "positive", "eigen")
      add("attractor", paste0("dominant_fraction_N_le_25_O2_", gsub("[^0-9]+", "p", O2)), sum(v[ngrid <= 25]), "fraction", "eigen")
      add("attractor", paste0("dominant_fraction_N_below_44_O2_", gsub("[^0-9]+", "p", O2)), sum(v[ngrid < 44]), "fraction", "eigen")
      add("attractor", paste0("dominant_growth_rate_O2_", gsub("[^0-9]+", "p", O2)), lambda1, "signed", "eigen")
      add("attractor", paste0("spectral_gap_O2_", gsub("[^0-9]+", "p", O2)), lambda1 - lambda2, "signed", "eigen")
      add("attractor", paste0("eigenvector_nonnegative_O2_", gsub("[^0-9]+", "p", O2)), as.numeric(nonneg), "binary", "eigen")
    }
    tri <- attr(G, "triplet")
    if (!is.null(tri$boundary_dropped_rate)) {
      for (N in states) {
        col <- as.integer(N - (cfg$N_MIN %||% 22) + 1L)
        add("buffering", paste0("transition_mass_dropped_at_N_", N, "_O2_", gsub("[^0-9]+", "p", O2)), tri$boundary_dropped_rate[[col]], "positive")
      }
    }
  }
  for (N in states) {
    pN <- as.numeric(o2ipa_call_model(model_env, ".pmisseg_of_O2", O2 = min(o2_grid), run_params = run_params, N = N))
    delta <- o2ipa_call_model(model_env, ".pr_delta_vec",
      N = N, p = pN, eps_tail = 1e-8, N_unit = as.integer(cfg$N_UNIT %||% 22),
      buffer_smax = as.numeric(run_params$buffer_smax %||% 1),
      buffer_beta = as.numeric(run_params$buffer_beta %||% 0),
      buffer_n_exp = as.numeric(run_params$buffer_n_exp %||% 1)
    )
    shifts <- suppressWarnings(as.numeric(names(delta)))
    prob <- as.numeric(delta)
    add("buffering", paste0("surviving_transition_mass_N_", N), sum(prob, na.rm = TRUE), "positive")
    add("buffering", paste0("mass_dropped_N_", N), as.numeric(attr(delta, "mass_dropped") %||% NA), "fraction")
    add("buffering", paste0("expected_delta_N_N_", N), sum(shifts * prob, na.rm = TRUE), "signed")
    add("buffering", paste0("expected_absolute_delta_N_N_", N), sum(abs(shifts) * prob, na.rm = TRUE), "positive")
    add("buffering", paste0("probability_surviving_below_source_N_", N), sum(prob[shifts < 0], na.rm = TRUE), "fraction")
    add("buffering", paste0("probability_surviving_above_source_N_", N), sum(prob[shifts > 0], na.rm = TRUE), "fraction")
    add("buffering", paste0("probability_entering_N_below_44_from_N_", N), sum(prob[N + shifts < 44], na.rm = TRUE), "fraction")
    add("buffering", paste0("probability_entering_N_le_25_from_N_", N), sum(prob[N + shifts <= 25], na.rm = TRUE), "fraction")
  }
  out <- do.call(rbind, rows)
  attr(out, "state_note") <- state_note
  out
}

o2pr_o2_exposure_one <- function(seed_id, run_params, model_env, cfg, burden_grid, realized_o2 = NULL, switch_o2 = NA_real_, scope = "full") {
  rows <- list()
  add <- function(module, feature, value, type = "positive", source = "model_probe") {
    rows[[length(rows) + 1L]] <<- o2ipa_feature_row(seed_id, scope, module, feature, value, type, source)
  }
  demand_grid <- burden_grid * as.numeric(run_params$rho_2N %||% 1) * ((44 / 44)^as.numeric(run_params$eta_o2 %||% 1))
  o2 <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = demand_grid, run_params = run_params, o2_Nref = as.numeric(cfg$o2_Nref %||% 1e6))
  for (i in seq_along(burden_grid)) {
    suffix <- paste0("_burden_", gsub("[^0-9]+", "p", burden_grid[[i]]))
    add("oxygen_feedback", paste0("O2_at_common", suffix), o2[[i]], "oxygen")
    add("oxygen_feedback", paste0("effective_demand_at_common", suffix), demand_grid[[i]], "positive")
  }
  for (target in c(2, 1, 0.5)) {
    dense_b <- 10^seq(log10(min(burden_grid)), log10(max(burden_grid)), length.out = 300L)
    dense_d <- dense_b * as.numeric(run_params$rho_2N %||% 1) * ((44 / 44)^as.numeric(run_params$eta_o2 %||% 1))
    dense_o <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = dense_d, run_params = run_params, o2_Nref = as.numeric(cfg$o2_Nref %||% 1e6))
    cr <- o2ipa_crossing(log10(dense_b), dense_o, target, decreasing = TRUE)
    target_key <- gsub("[^0-9]+", "p", target)
    add("oxygen_feedback", paste0("burden_reaches_O2_", target_key, "pct"), as.numeric(is.finite(cr$value)), "binary")
    add("oxygen_feedback", paste0("log10_burden_at_O2_", target_key, "pct"), cr$value, "signed")
  }
  add("oxygen_feedback", "dO2_dlogBurden_median", stats::median(diff(o2) / diff(log10(demand_grid + 1)), na.rm = TRUE), "signed")
  add("oxygen_feedback", "minimum_O2_floor", min(o2, na.rm = TRUE), "oxygen")
  if (!is.null(realized_o2) && nrow(realized_o2)) {
    r <- realized_o2[realized_o2$seed_id == seed_id, , drop = FALSE]
    if (nrow(r)) {
      add("oxygen_feedback", "minimum_realized_O2", min(r$pred_o2_pct, na.rm = TRUE), "oxygen", "realized_prediction")
      add("oxygen_feedback", "terminal_O2", tail(r$pred_o2_pct[order(r$day)], 1), "oxygen", "realized_prediction")
      add("oxygen_feedback", "O2_AUC", o2pr_auc(r$day, r$pred_o2_pct), "positive", "realized_prediction")
      for (target in c(2, 1, 0.5)) {
        rr <- r[order(r$day), , drop = FALSE]
        hit <- rr$day[rr$pred_o2_pct < target]
        first_hit <- if (length(hit)) hit[[1]] else NA_real_
        target_key <- gsub("[^0-9]+", "p", target)
        add("oxygen_feedback", paste0("ever_below_O2_", target_key, "pct"), as.numeric(is.finite(first_hit)), "binary", "realized_prediction")
        add("oxygen_feedback", paste0("time_below_O2_", target_key, "pct"), first_hit, "time", "realized_prediction")
      }
      min_realized_o2 <- min(r$pred_o2_pct, na.rm = TRUE)
      switch_available <- is.finite(switch_o2)
      add("oxygen_feedback", "switch_available", as.numeric(switch_available), "binary", "switch_diagnostic")
      add("oxygen_feedback", "switch_margin", if (switch_available && is.finite(min_realized_o2)) min_realized_o2 - switch_o2 else NA_real_, "signed", "switch_diagnostic")
      add("oxygen_feedback", "switch_crossed", as.numeric(switch_available && is.finite(min_realized_o2) && min_realized_o2 <= switch_o2), "binary", "switch_diagnostic")
    }
  }
  do.call(rbind, rows)
}

o2pr_read_realized_o2 <- function(inputs) {
  rows <- list()
  for (i in seq_len(nrow(inputs$manifest))) {
    seed <- inputs$manifest$seed_id[[i]]
    d <- inputs$manifest$seed_dir[[i]]
    p <- file.path(d, "viz", "predict_burden_0_1000day.tsv")
    if (!file.exists(p)) p <- file.path(d, "viz", "predict_burden_vs_o2.tsv")
    if (!file.exists(p)) next
    tab <- tryCatch(o2pr_read_tsv(p), error = function(e) data.frame())
    if (!all(c("day", "pred_o2_pct") %in% names(tab))) next
    tab$seed_id <- seed
    rows[[seed]] <- tab[, intersect(c("seed_id", "cohort", "day", "pred_o2_pct", "pred_o2_target_pct"), names(tab)), drop = FALSE]
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

o2pr_build_process_outputs <- function(inputs, out_dir, model_env, cfg, o2_grid, force = FALSE) {
  manifest <- inputs$manifest[inputs$manifest$fit_success %in% TRUE, , drop = FALSE]
  param_mat <- o2ipa_params_wide(inputs$params_long, "value")
  best_seed <- manifest$seed_id[which.min(manifest$objective)]
  ref <- as.numeric(param_mat[best_seed, , drop = TRUE])
  names(ref) <- colnames(param_mat)
  target <- o2ipa_target_params()
  real_o2 <- o2pr_read_realized_o2(inputs)
  burden_grid <- c(0.01, 0.1, 1, 10, 100)
  intr_full <- list()
  intr_18 <- list()
  exp_full <- list()
  exp_18 <- list()
  state_note <- NULL
  for (seed in manifest$seed_id) {
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    full_rp <- o2pr_run_params_from_vec(pvec, cfg)
    only_rp <- ref
    only_rp[target] <- pvec[target]
    only_rp <- o2pr_run_params_from_vec(only_rp, cfg)
    full_intr <- o2pr_build_process_one(seed, full_rp, model_env, cfg, o2_grid, "intrinsic_full")
    only_intr <- o2pr_build_process_one(seed, only_rp, model_env, cfg, o2_grid, "intrinsic_18only")
    if (is.null(state_note)) state_note <- attr(full_intr, "state_note")
    intr_full[[seed]] <- full_intr
    intr_18[[seed]] <- only_intr
    switch_o2 <- o2pr_switch_from_intrinsic(full_intr)
    exp_full[[seed]] <- o2pr_o2_exposure_one(seed, full_rp, model_env, cfg, burden_grid, real_o2, switch_o2, "exposure_full")
    exp_18[[seed]] <- o2pr_o2_exposure_one(seed, only_rp, model_env, cfg, burden_grid, real_o2, switch_o2, "exposure_18only")
  }
  intrinsic_full <- do.call(rbind, intr_full)
  intrinsic_18 <- do.call(rbind, intr_18)
  exposure_full <- do.call(rbind, exp_full)
  exposure_18 <- do.call(rbind, exp_18)
  o2pr_write_tsv(intrinsic_full, file.path(out_dir, "tables", "process_intrinsic_full_long.tsv"))
  o2pr_write_tsv(intrinsic_18, file.path(out_dir, "tables", "process_intrinsic_18only_long.tsv"))
  o2pr_write_tsv(exposure_full, file.path(out_dir, "tables", "process_exposure_full_long.tsv"))
  o2pr_write_tsv(exposure_18, file.path(out_dir, "tables", "process_exposure_18only_long.tsv"))
  o2pr_write_tsv(intrinsic_full, file.path(out_dir, "tables", "process_intrinsic_long.tsv"))
  o2pr_write_tsv(exposure_full, file.path(out_dir, "tables", "process_exposure_long.tsv"))
  o2pr_write_tsv(state_note, file.path(out_dir, "tables", "state_anchor_mapping.tsv"))
  ref_df <- data.frame(parameter = names(ref), value = as.numeric(ref), reference_seed = best_seed, stringsAsFactors = FALSE)
  o2pr_write_tsv(ref_df, file.path(out_dir, "tables", "reference_parameter_set.tsv"))
  list(
    intrinsic_full = intrinsic_full,
    intrinsic_18 = intrinsic_18,
    exposure_full = exposure_full,
    exposure_18 = exposure_18,
    realized_o2 = real_o2
  )
}

o2pr_switch_from_intrinsic <- function(long_df) {
  d <- long_df[long_df$feature_type == "signed" & grepl("^selection_22_vs_44_O2_", long_df$feature), , drop = FALSE]
  if (nrow(d) < 2L) return(NA_real_)
  o2 <- suppressWarnings(as.numeric(sub("^.*_O2_", "", gsub("p", ".", d$feature))))
  y <- d$raw_value
  ok <- is.finite(o2) & is.finite(y)
  o2 <- o2[ok]
  y <- y[ok]
  if (length(o2) < 2L || min(y) > 0 || max(y) < 0) return(NA_real_)
  ord <- order(o2)
  o2 <- o2[ord]
  y <- y[ord]
  as.numeric(stats::approx(y, o2, xout = 0, ties = "ordered")$y)
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
    by_day <- aggregate(dd$fraction, by = list(day = dd$day), FUN = sum, na.rm = TRUE)
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

o2pr_basic_figures <- function(traj, scaled, clusters, labels, out_dir) {
  fig <- function(name) file.path(out_dir, "figures", paste0(name, ".pdf"))
  if (is.data.frame(traj$curves) && nrow(traj$curves)) {
    d <- merge(traj$curves, labels[, c("seed_id", "trajectory_regime")], by = "seed_id", all.x = TRUE)
    grDevices::pdf(fig("all_trajectories_by_regime"), width = 9, height = 6)
    cols <- c(stable_high_chr = "#1b9e77", late_collapse_to_low_chr = "#d95f02", ambiguous = "grey60")
    plot(NA, xlim = range(d$day, na.rm = TRUE), ylim = range(d$trajectory_value, na.rm = TRUE), xlab = "Day", ylab = "Weighted mean chromosome number", main = "Trajectories by regime")
    for (key in split(seq_len(nrow(d)), paste(d$seed_id, d$cohort))) {
      x <- d[key, ]
      lines(x$day, x$trajectory_value, col = adjustcolor(cols[x$trajectory_regime[[1]]] %||% "grey60", alpha.f = 0.35))
    }
    legend("topright", legend = names(cols), col = cols, lty = 1, bty = "n")
    grDevices::dev.off()
  }
  if (!is.null(scaled$combined_full) && nrow(scaled$combined_full$wide) >= 2L && ncol(scaled$combined_full$wide) >= 3L) {
    mat <- as.matrix(scaled$combined_full$wide[, setdiff(names(scaled$combined_full$wide), "seed_id"), drop = FALSE])
    rownames(mat) <- scaled$combined_full$wide$seed_id
    mat <- mat[rowSums(!is.finite(mat)) == 0, colSums(!is.finite(mat)) == 0, drop = FALSE]
    if (nrow(mat) >= 2L && ncol(mat) >= 2L) {
      grDevices::pdf(fig("process_PCA"), width = 6, height = 5)
      pc <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
      lab <- labels$trajectory_regime[match(rownames(mat), labels$seed_id)]
      col <- c(stable_high_chr = "#1b9e77", late_collapse_to_low_chr = "#d95f02", ambiguous = "grey60")[lab]
      plot(pc$x[, 1], pc$x[, 2], pch = 16, col = col, xlab = "PC1", ylab = "PC2", main = "Process PCA")
      grDevices::dev.off()
      grDevices::pdf(fig("intrinsic_process_fingerprint_heatmap"), width = 9, height = 7)
      heatmap(mat, scale = "none", margins = c(5, 5), main = "Combined process fingerprint")
      grDevices::dev.off()
            if (nrow(mat) >= 3L) {
              grDevices::pdf(fig("process_dendrogram"), width = 7, height = 5)
              plot(stats::hclust(stats::dist(mat), method = "ward.D2"), main = "Process dendrogram", xlab = "")
              grDevices::dev.off()
            }
        }
    }
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

o2pr_report <- function(out_dir, run_dir, repo_root, manifest, cfg_meta, labels, unavailable = character()) {
  branch <- o2pr_git_value(c("branch", "--show-current"), repo_root)
  commit <- o2pr_git_value(c("rev-parse", "HEAD"), repo_root)
  lines <- c(
    "# In vivo ploidy regime analysis summary",
    "",
    paste0("- git branch: ", branch),
    paste0("- git commit: ", commit),
    paste0("- input run_dir: ", run_dir),
    paste0("- output out_dir: ", out_dir),
    paste0("- discovered seeds: ", nrow(manifest)),
    paste0("- valid seeds: ", sum(manifest$fit_success %in% TRUE, na.rm = TRUE)),
    paste0("- trajectory stable_high_chr: ", sum(labels$trajectory_regime == "stable_high_chr", na.rm = TRUE)),
    paste0("- trajectory late_collapse_to_low_chr: ", sum(labels$trajectory_regime == "late_collapse_to_low_chr", na.rm = TRUE)),
    paste0("- trajectory ambiguous: ", sum(labels$trajectory_regime == "ambiguous", na.rm = TRUE)),
    "",
    "## Model/schema notes",
    "",
    paste0("- ", cfg_meta$metric, ": ", cfg_meta$value),
    "",
    "## Scientific interpretation boundaries",
    "",
    "- Different optimizer seeds are not posterior samples; seed proportions cannot be interpreted as biological occurrence probabilities.",
    "- If two modes separate only beyond the observed time window, the result indicates long-horizon prediction nonidentifiability rather than proof of two true biological subtypes.",
    "- Process clustering can only identify mechanistic regimes supported by this model; real biological heterogeneity still requires independent data validation.",
    "",
    "## Unavailable or deferred diagnostics",
    "",
    if (length(unavailable)) paste0("- ", unavailable) else "- See unavailable_boundary_tests.tsv and placeholder advanced-analysis tables."
  )
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}
