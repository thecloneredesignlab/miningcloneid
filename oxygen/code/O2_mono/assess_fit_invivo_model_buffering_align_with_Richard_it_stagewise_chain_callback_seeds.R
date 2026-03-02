#!/usr/bin/env Rscript

`%||%` <- function(a, b) if (is.null(a)) b else a

parse_args <- function(argv) {
  out <- list()
  if (length(argv) == 0) return(out)
  for (a in argv) {
    if (!startsWith(a, "--")) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.numeric(x))
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.integer(x))
}

safe_read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.delim(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}

read_named_metric <- function(path, key_col, val_col) {
  tab <- safe_read_tsv(path)
  if (is.null(tab)) return(NULL)
  if (!all(c(key_col, val_col) %in% names(tab))) return(NULL)
  vals <- suppressWarnings(as.numeric(tab[[val_col]]))
  setNames(vals, as.character(tab[[key_col]]))
}

first_scalar <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0) return(default)
  x[[1]]
}

rank_num <- function(x) {
  out <- rep(NA_real_, length(x))
  idx <- which(is.finite(x))
  if (length(idx) > 0) out[idx] <- rank(x[idx], ties.method = "min")
  out
}

join_flags <- function(flags) {
  flags <- unique(flags[nzchar(flags)])
  if (length(flags) == 0) "" else paste(flags, collapse = ";")
}

derive_prefix_root <- function(metrics_tsv = NULL, run_prefix_root = NULL) {
  if (!is.null(run_prefix_root)) return(normalizePath(run_prefix_root, mustWork = FALSE))
  if (is.null(metrics_tsv)) stop("Need --metrics_tsv or --run_prefix_root")
  m <- normalizePath(metrics_tsv, mustWork = TRUE)
  sub("_callback_metrics\\.tsv$", "", m)
}

default_out_dir <- function(prefix_root) {
  paste0(prefix_root, "_seed_review")
}

extract_seed_from_dirname <- function(seed_dir) {
  bn <- basename(seed_dir)
  sub("^.*_seed", "", bn)
}

compute_callback_metrics_row <- function(callback_dir, seed) {
  sum_path <- file.path(callback_dir, "fit_summary.tsv")
  if (!file.exists(sum_path)) {
    return(data.frame(
      seed = as.character(seed),
      objective = NA_real_,
      objective_burden = NA_real_,
      objective_ploidy = NA_real_,
      rmse_4N_burden = NA_real_,
      mean_nll_4N_ploidy = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  vals <- read_named_metric(sum_path, "metric", "value")
  obj <- if (!is.null(vals)) as.numeric(vals[["objective"]]) else NA_real_
  obj_b <- if (!is.null(vals)) as.numeric(vals[["objective_burden"]]) else NA_real_
  obj_p <- if (!is.null(vals)) as.numeric(vals[["objective_ploidy"]]) else NA_real_

  rmse_4n <- NA_real_
  bf <- file.path(callback_dir, "burden_fit.tsv")
  b <- safe_read_tsv(bf)
  if (!is.null(b) && all(c("cohort", "pred_norm", "obs_norm") %in% names(b))) {
    idx <- which(as.character(b$cohort) == "4N" & is.finite(b$pred_norm) & is.finite(b$obs_norm))
    if (length(idx) > 0) {
      rmse_4n <- sqrt(mean((as.numeric(b$pred_norm[idx]) - as.numeric(b$obs_norm[idx]))^2))
    }
  }

  nll_4n <- NA_real_
  pf <- file.path(callback_dir, "terminal_ploidy_fit.tsv")
  p <- safe_read_tsv(pf)
  if (!is.null(p) && all(c("cohort", "pred_fraction", "obs_count") %in% names(p))) {
    p <- p[as.character(p$cohort) == "4N", , drop = FALSE]
    if (nrow(p) > 0) {
      if (all(c("harvest", "dose") %in% names(p))) {
        key <- interaction(p$harvest, p$dose, drop = TRUE)
      } else {
        key <- rep(1L, nrow(p))
      }
      split_idx <- split(seq_len(nrow(p)), key)
      nll_each <- vapply(split_idx, function(ix) {
        x <- p[ix, , drop = FALSE]
        obs <- as.numeric(x$obs_count)
        if (!any(is.finite(obs)) || sum(obs, na.rm = TRUE) <= 0) return(NA_real_)
        pred <- pmax(as.numeric(x$pred_fraction), 1e-12)
        -sum(obs * log(pred), na.rm = TRUE) / sum(obs, na.rm = TRUE)
      }, numeric(1))
      nll_4n <- mean(nll_each, na.rm = TRUE)
      if (!is.finite(nll_4n)) nll_4n <- NA_real_
    }
  }

  data.frame(
    seed = as.character(seed),
    objective = obj,
    objective_burden = obj_b,
    objective_ploidy = obj_p,
    rmse_4N_burden = rmse_4n,
    mean_nll_4N_ploidy = nll_4n,
    stringsAsFactors = FALSE
  )
}

rebuild_metrics_from_callbacks <- function(prefix_root, callback_dir_name = "callback_equal") {
  patt <- paste0(prefix_root, "_seed*")
  seed_dirs <- Sys.glob(patt)
  seed_dirs <- seed_dirs[dir.exists(seed_dirs)]
  if (length(seed_dirs) == 0) {
    stop("No seed directories found for prefix_root: ", prefix_root)
  }
  rows <- vector("list", length(seed_dirs))
  for (i in seq_along(seed_dirs)) {
    sd <- seed_dirs[[i]]
    seed <- extract_seed_from_dirname(sd)
    cb <- file.path(sd, callback_dir_name)
    rows[[i]] <- compute_callback_metrics_row(cb, seed = seed)
  }
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) stop("Failed to rebuild callback metrics from callback directories.")
  out$seed_num <- suppressWarnings(as.numeric(out$seed))
  ord <- order(ifelse(is.finite(out$seed_num), out$seed_num, Inf), out$seed)
  out <- out[ord, c("seed", "objective", "objective_burden", "objective_ploidy", "rmse_4N_burden", "mean_nll_4N_ploidy"), drop = FALSE]
  rownames(out) <- NULL
  out
}

extract_params <- function(callback_dir, rho_2N_min, rho_2N_max, legacy_c_vol_2N_mm3_fallback) {
  best_path <- file.path(callback_dir, "best_params.tsv")
  cfg_path <- file.path(callback_dir, "fit_config.rds")
  fit_sum_path <- file.path(callback_dir, "fit_summary.tsv")

  flags <- character(0)
  pvals <- read_named_metric(best_path, "parameter", "value")
  if (is.null(pvals)) {
    flags <- c(flags, "missing_best_params")
    return(list(
      rho_2N = NA_real_, rho_source = NA_character_, beta_size = NA_real_, ratio_4N_2N = NA_real_,
      lam_min = NA_real_, lam_max = NA_real_, p_wgd = NA_real_,
      rho_near_lower = NA, rho_near_upper = NA, rho_out_of_range = NA,
      beta_size_near_lower = NA, beta_size_near_upper = NA,
      flags = join_flags(flags)
    ))
  }

  cfg_rds <- NULL
  if (file.exists(cfg_path)) {
    cfg_rds <- tryCatch(readRDS(cfg_path), error = function(e) NULL)
  }
  legacy_c_vol <- first_scalar(
    c(
      if (!is.null(cfg_rds)) cfg_rds$legacy_c_vol_2N_mm3 else NULL,
      if (!is.null(cfg_rds)) cfg_rds$c_vol_2N_mm3 else NULL,
      legacy_c_vol_2N_mm3_fallback
    ),
    default = legacy_c_vol_2N_mm3_fallback
  )
  if (!is.finite(legacy_c_vol) || legacy_c_vol <= 0) legacy_c_vol <- legacy_c_vol_2N_mm3_fallback

  rho_2N <- NA_real_
  rho_source <- NA_character_
  if ("rho_2N" %in% names(pvals) && is.finite(pvals[["rho_2N"]]) && pvals[["rho_2N"]] > 0) {
    rho_2N <- as.numeric(pvals[["rho_2N"]])
    rho_source <- "best_params.rho_2N"
  } else if ("c_scale" %in% names(pvals) && is.finite(pvals[["c_scale"]]) && pvals[["c_scale"]] > 0) {
    rho_2N <- 1 / (legacy_c_vol * as.numeric(pvals[["c_scale"]]))
    rho_source <- "legacy_c_scale_converted"
    flags <- c(flags, "legacy_c_scale_param")
  } else {
    flags <- c(flags, "missing_rho_2N")
  }

  beta_size <- if ("beta_size" %in% names(pvals)) as.numeric(pvals[["beta_size"]]) else NA_real_
  ratio_4N_2N <- if (is.finite(beta_size)) 2^beta_size else NA_real_

  # Boundary heuristics
  boundary_rel_tol <- 0.02
  rho_near_lower <- if (is.finite(rho_2N)) rho_2N <= rho_2N_min * (1 + boundary_rel_tol) else NA
  rho_near_upper <- if (is.finite(rho_2N)) rho_2N >= rho_2N_max * (1 - boundary_rel_tol) else NA
  rho_out_of_range <- if (is.finite(rho_2N)) (rho_2N < rho_2N_min || rho_2N > rho_2N_max) else NA
  beta_size_near_lower <- if (is.finite(beta_size)) beta_size <= 0 + 0.02 else NA
  beta_size_near_upper <- if (is.finite(beta_size)) beta_size >= 2 - 0.02 else NA

  if (isTRUE(rho_out_of_range)) flags <- c(flags, "rho_2N_out_of_range")
  if (isTRUE(rho_near_lower)) flags <- c(flags, "rho_2N_near_lower")
  if (isTRUE(rho_near_upper)) flags <- c(flags, "rho_2N_near_upper")
  if (isTRUE(beta_size_near_lower)) flags <- c(flags, "beta_size_near_lower")
  if (isTRUE(beta_size_near_upper)) flags <- c(flags, "beta_size_near_upper")

  # Fit summary sanity
  fsum <- read_named_metric(fit_sum_path, "metric", "value")
  if (is.null(fsum)) {
    flags <- c(flags, "missing_fit_summary")
  }

  list(
    rho_2N = rho_2N,
    rho_source = rho_source,
    beta_size = beta_size,
    ratio_4N_2N = ratio_4N_2N,
    lam_min = if ("lam_min" %in% names(pvals)) as.numeric(pvals[["lam_min"]]) else NA_real_,
    lam_max = if ("lam_max" %in% names(pvals)) as.numeric(pvals[["lam_max"]]) else NA_real_,
    p_wgd = if ("p_wgd" %in% names(pvals)) as.numeric(pvals[["p_wgd"]]) else NA_real_,
    rho_near_lower = rho_near_lower,
    rho_near_upper = rho_near_upper,
    rho_out_of_range = rho_out_of_range,
    beta_size_near_lower = beta_size_near_lower,
    beta_size_near_upper = beta_size_near_upper,
    flags = join_flags(flags)
  )
}

make_recommendation <- function(df, objective_rel_tol = 0.03) {
  if (nrow(df) == 0) return(list(recommended_seed = NA, note = "No rows."))
  finite_obj <- is.finite(df$objective)
  if (!any(finite_obj)) return(list(recommended_seed = NA, note = "No finite objective values."))

  best_obj <- min(df$objective[finite_obj], na.rm = TRUE)
  candidates <- which(finite_obj & df$objective <= best_obj * (1 + objective_rel_tol))
  if (length(candidates) == 0) candidates <- which.min(df$objective)

  ord <- order(
    df$penalty_count[candidates],
    df$score[candidates],
    df$objective[candidates],
    na.last = TRUE
  )
  pick <- candidates[[ord[[1]]]]
  note <- if (length(candidates) == 1) {
    sprintf("Unique candidate within %.1f%% objective tolerance.", 100 * objective_rel_tol)
  } else {
    sprintf("%d seeds within %.1f%% objective tolerance; selected by penalty_count -> score -> objective.", length(candidates), 100 * objective_rel_tol)
  }
  list(recommended_seed = df$seed[[pick]], best_objective_seed = df$seed[[which.min(df$objective)]], note = note, candidate_rows = candidates)
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  metrics_tsv <- argv$metrics_tsv %||% NULL
  prefix_root <- derive_prefix_root(metrics_tsv = metrics_tsv, run_prefix_root = argv$run_prefix_root %||% NULL)

  callback_dir_name <- as.character(argv$callback_dir_name %||% "callback_equal")
  out_dir <- normalizePath(argv$out_dir %||% default_out_dir(prefix_root), mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rho_2N_min <- as_num(argv$rho_2N_min, 3.2e4)
  rho_2N_max <- as_num(argv$rho_2N_max, 5.6e4)
  if (!is.finite(rho_2N_min) || rho_2N_min <= 0) rho_2N_min <- 3.2e4
  if (!is.finite(rho_2N_max) || rho_2N_max <= 0) rho_2N_max <- 5.6e4
  if (rho_2N_min > rho_2N_max) {
    tmp <- rho_2N_min; rho_2N_min <- rho_2N_max; rho_2N_max <- tmp
  }
  legacy_c_vol_2N_mm3 <- as_num(argv$legacy_c_vol_2N_mm3, 4.19e-06)
  if (!is.finite(legacy_c_vol_2N_mm3) || legacy_c_vol_2N_mm3 <= 0) legacy_c_vol_2N_mm3 <- 4.19e-06
  objective_rel_tol <- as_num(argv$objective_rel_tol, 0.03)
  if (!is.finite(objective_rel_tol) || objective_rel_tol < 0) objective_rel_tol <- 0.03

  if (is.null(metrics_tsv)) metrics_tsv <- paste0(prefix_root, "_callback_metrics.tsv")
  metrics_path_exists <- file.exists(metrics_tsv)
  m <- NULL
  metrics_source <- NULL
  rebuilt_metrics_tsv <- NA_character_
  if (metrics_path_exists) {
    metrics_tsv <- normalizePath(metrics_tsv, mustWork = TRUE)
    m <- safe_read_tsv(metrics_tsv)
    metrics_source <- "file"
    if (is.null(m)) stop("Failed to read metrics_tsv: ", metrics_tsv)
  } else {
    message("callback_metrics.tsv not found; rebuilding from per-seed callback directories.")
    m <- rebuild_metrics_from_callbacks(prefix_root = prefix_root, callback_dir_name = callback_dir_name)
    rebuilt_metrics_tsv <- file.path(out_dir, "callback_metrics_rebuilt.tsv")
    write.table(m, file = rebuilt_metrics_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    metrics_source <- "rebuilt_from_callbacks"
    metrics_tsv <- normalizePath(metrics_tsv, mustWork = FALSE)
  }
  req <- c("seed", "objective", "objective_burden", "objective_ploidy", "rmse_4N_burden", "mean_nll_4N_ploidy")
  miss <- setdiff(req, names(m))
  if (length(miss) > 0) stop("metrics_tsv missing columns: ", paste(miss, collapse = ", "))

  rows <- vector("list", nrow(m))
  for (i in seq_len(nrow(m))) {
    seed_chr <- as.character(m$seed[[i]])
    seed_dir <- paste0(prefix_root, "_seed", seed_chr)
    callback_dir <- file.path(seed_dir, callback_dir_name)

    flags <- character(0)
    if (!dir.exists(seed_dir)) flags <- c(flags, "missing_seed_dir")
    if (!dir.exists(callback_dir)) flags <- c(flags, "missing_callback_dir")

    p <- extract_params(
      callback_dir = callback_dir,
      rho_2N_min = rho_2N_min,
      rho_2N_max = rho_2N_max,
      legacy_c_vol_2N_mm3_fallback = legacy_c_vol_2N_mm3
    )
    if (nzchar(p$flags)) flags <- c(flags, strsplit(p$flags, ";", fixed = TRUE)[[1]])

    rows[[i]] <- data.frame(
      seed = seed_chr,
      seed_dir = seed_dir,
      callback_dir = callback_dir,
      objective = as_num(m$objective[[i]], NA_real_),
      objective_burden = as_num(m$objective_burden[[i]], NA_real_),
      objective_ploidy = as_num(m$objective_ploidy[[i]], NA_real_),
      rmse_4N_burden = as_num(m$rmse_4N_burden[[i]], NA_real_),
      mean_nll_4N_ploidy = as_num(m$mean_nll_4N_ploidy[[i]], NA_real_),
      rho_2N = p$rho_2N,
      rho_source = p$rho_source %||% NA_character_,
      beta_size = p$beta_size,
      ratio_4N_2N = p$ratio_4N_2N,
      lam_min = p$lam_min,
      lam_max = p$lam_max,
      p_wgd = p$p_wgd,
      rho_near_lower = p$rho_near_lower,
      rho_near_upper = p$rho_near_upper,
      rho_out_of_range = p$rho_out_of_range,
      beta_size_near_lower = p$beta_size_near_lower,
      beta_size_near_upper = p$beta_size_near_upper,
      flags = join_flags(flags),
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, rows)

  df$rank_objective <- rank_num(df$objective)
  df$rank_rmse_4N_burden <- rank_num(df$rmse_4N_burden)
  df$rank_mean_nll_4N_ploidy <- rank_num(df$mean_nll_4N_ploidy)
  df$penalty_count <- vapply(strsplit(df$flags, ";", fixed = TRUE), function(x) {
    x <- x[nzchar(x)]
    # Count only quality-risk flags (exclude legacy format flag).
    risk <- setdiff(x, "legacy_c_scale_param")
    length(risk)
  }, integer(1))

  df$score <- with(df,
    rank_objective +
      0.5 * ifelse(is.finite(rank_rmse_4N_burden), rank_rmse_4N_burden, max(rank_rmse_4N_burden, na.rm = TRUE) + 1) +
      0.5 * ifelse(is.finite(rank_mean_nll_4N_ploidy), rank_mean_nll_4N_ploidy, max(rank_mean_nll_4N_ploidy, na.rm = TRUE) + 1) +
      2.0 * penalty_count
  )
  df$score_rank <- rank_num(df$score)

  rec <- make_recommendation(df, objective_rel_tol = objective_rel_tol)
  df$within_objective_tol <- FALSE
  if (!is.null(rec$candidate_rows) && length(rec$candidate_rows) > 0) {
    df$within_objective_tol[rec$candidate_rows] <- TRUE
  }
  df$recommended <- df$seed == rec$recommended_seed
  df$best_objective_seed <- df$seed == rec$best_objective_seed

  ord <- order(df$recommended * -1, df$score, df$objective, na.last = TRUE)
  ranked <- df[ord, , drop = FALSE]

  summary_tsv <- file.path(out_dir, "seed_quality_summary.tsv")
  ranked_tsv <- file.path(out_dir, "seed_quality_ranked.tsv")
  rec_txt <- file.path(out_dir, "seed_quality_recommendation.txt")

  write.table(df, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(ranked, file = ranked_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  lines <- c(
    paste0("metrics_tsv\t", metrics_tsv),
    paste0("metrics_source\t", metrics_source),
    paste0("rebuilt_metrics_tsv\t", ifelse(is.na(rebuilt_metrics_tsv), "", rebuilt_metrics_tsv)),
    paste0("prefix_root\t", prefix_root),
    paste0("callback_dir_name\t", callback_dir_name),
    paste0("rho_2N_range_cells_per_mm3\t[", signif(rho_2N_min, 6), ", ", signif(rho_2N_max, 6), "]"),
    paste0("objective_rel_tol\t", objective_rel_tol),
    paste0("recommended_seed\t", rec$recommended_seed %||% "NA"),
    paste0("best_objective_seed\t", rec$best_objective_seed %||% "NA"),
    paste0("note\t", rec$note %||% ""),
    "",
    "Top ranked seeds (score, objective, penalties, flags):"
  )
  top_n <- min(5L, nrow(ranked))
  for (i in seq_len(top_n)) {
    r <- ranked[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "seed=%s\tscore=%.3f\tobjective=%s\tpenalty_count=%d\twithin_tol=%s\tflags=%s\trho_2N=%s\tbeta_size=%s\tratio_4N_2N=%s",
      r$seed[[1]],
      as.numeric(r$score[[1]]),
      ifelse(is.finite(r$objective[[1]]), format(r$objective[[1]], digits = 8), "NA"),
      as.integer(r$penalty_count[[1]]),
      as.character(r$within_objective_tol[[1]]),
      r$flags[[1]],
      ifelse(is.finite(r$rho_2N[[1]]), format(r$rho_2N[[1]], digits = 8), "NA"),
      ifelse(is.finite(r$beta_size[[1]]), format(r$beta_size[[1]], digits = 6), "NA"),
      ifelse(is.finite(r$ratio_4N_2N[[1]]), format(r$ratio_4N_2N[[1]], digits = 6), "NA")
    ))
  }
  writeLines(lines, con = rec_txt)

  message("Wrote: ", normalizePath(summary_tsv, mustWork = FALSE))
  message("Wrote: ", normalizePath(ranked_tsv, mustWork = FALSE))
  message("Wrote: ", normalizePath(rec_txt, mustWork = FALSE))
  if (!is.na(rebuilt_metrics_tsv)) message("Wrote: ", normalizePath(rebuilt_metrics_tsv, mustWork = FALSE))
  message("Recommended seed: ", rec$recommended_seed %||% "NA")
}

if (sys.nframe() == 0) {
  main()
}
