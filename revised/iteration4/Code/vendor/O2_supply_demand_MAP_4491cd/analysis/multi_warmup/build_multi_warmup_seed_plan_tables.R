#!/usr/bin/env Rscript

# Table-only seed-plan analysis.  Coordinate tables are materialized here;
# figure rendering is owned by vis/multi_warmup/.

.o2mw_seed_analysis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "build_multi_warmup_seed_plan_tables.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_seed_workflow_root <- normalizePath(
  file.path(.o2mw_seed_analysis_dir, "..", ".."),
  mustWork = TRUE
)
source(
  file.path(.o2mw_seed_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"),
  local = environment()
)

validate_top_n <- function(name, value) {
  if (!is.finite(value) || value < 0L) {
    stop(name, " must be a non-negative integer, got: ", value, call. = FALSE)
  }
  as.integer(value)
}

seed_integer <- function(x) {
  out <- seed_order_key(x)
  out[!is.finite(out)] <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", x[!is.finite(out)])))
  as.integer(out)
}

ratio_parameters <- c(
  "O2_crit", "alpha_o2", "mu_hp", "p_misseg", "k_o_mis",
  "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O",
  "gamma_growth", "lam_max", "p_mis_base", "p_wgd", "gamma_mu"
)
ratio_parameter_scales <- c(
  O2_crit = "log10",
  alpha_o2 = "log10",
  mu_hp = "log10",
  p_misseg = "log10",
  k_o_mis = "log10",
  buffer_smax = "identity",
  buffer_beta = "log10",
  buffer_n_exp = "log10",
  n_O = "identity",
  gamma_growth = "identity",
  lam_max = "log10",
  p_mis_base = "log10",
  p_wgd = "log10",
  gamma_mu = "identity"
)
ratio_feature_cols <- paste0("log10_ratio_", ratio_parameters)
invivo_only_oxygen_parameters <- c("o2_S0", "kappa_O", "eta_o2", "rho_2N")
invivo_only_oxygen_feature_cols <- paste0("log10_invivo_", invivo_only_oxygen_parameters)
invivo_profile_feature_cols <- c(ratio_feature_cols, invivo_only_oxygen_feature_cols)

single_side_ratio_feature_cols <- function(side) {
  side <- match.arg(side, c("invivo", "invitro"))
  scales <- ratio_parameter_scales[ratio_parameters]
  paste0(ifelse(scales == "log10", "log10_", ""), side, "_", ratio_parameters)
}

single_side_feature_cols <- function(side) {
  side <- match.arg(side, c("invivo", "invitro"))
  out <- single_side_ratio_feature_cols(side)
  if (identical(side, "invivo")) out <- c(out, invivo_only_oxygen_feature_cols)
  out
}

single_side_feature_sources <- function(side) {
  side <- match.arg(side, c("invivo", "invitro"))
  scales <- ratio_parameter_scales[ratio_parameters]
  out <- paste0(side, "_", ifelse(scales == "log10", "log10_best_param", "best_param"))
  if (identical(side, "invivo")) {
    out <- c(out, rep("invivo_only_log10_best_param", length(invivo_only_oxygen_feature_cols)))
  }
  out
}

read_metric_map <- function(path) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("metric", "value") %in% names(tab))) stop("Missing metric/value columns: ", path, call. = FALSE)
  setNames(tab$value, tab$metric)
}

metric_value <- function(vals, key, default = NA_character_) {
  if (is.null(vals) || !(key %in% names(vals))) return(default)
  vals[[key]]
}

find_seed_dirs <- function(run_dir) {
  dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "fit_summary.tsv")) & file.exists(file.path(dirs, "best_params.tsv"))]
  dirs[order(seed_order_key(basename(dirs)), basename(dirs))]
}

read_seed_summary <- function(seed_dir) {
  vals <- read_metric_map(file.path(seed_dir, "fit_summary.tsv"))
  objective <- suppressWarnings(as.numeric(metric_value(vals, "objective", metric_value(vals, "objective_total", NA))))
  data.frame(
    seed = basename(seed_dir),
    seed_int = seed_integer(basename(seed_dir)),
    seed_dir = normalizePath(seed_dir, mustWork = TRUE),
    objective = objective,
    stringsAsFactors = FALSE
  )
}

read_best_params <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("parameter", "value") %in% names(tab))) stop("Missing parameter/value columns: ", path, call. = FALSE)
  vals <- suppressWarnings(as.numeric(tab$value))
  names(vals) <- as.character(tab$parameter)
  vals
}

read_top_seed_table <- function(run_dir, top_n, label) {
  seed_dirs <- find_seed_dirs(run_dir)
  if (!length(seed_dirs)) stop("No valid seed dirs found for ", label, ": ", run_dir, call. = FALSE)
  tab <- do.call(rbind, lapply(seed_dirs, read_seed_summary))
  tab <- tab[is.finite(tab$objective), , drop = FALSE]
  if (nrow(tab) < top_n) {
    stop(label, " source run has only ", nrow(tab), " finite-objective seeds; need at least ", top_n, ".", call. = FALSE)
  }
  tab <- tab[order(tab$objective, tab$seed_int, tab$seed), , drop = FALSE]
  tab <- utils::head(tab, top_n)
  tab$rank <- seq_len(nrow(tab))
  row.names(tab) <- NULL
  tab
}

params_for_seed <- function(seed_dir) {
  vals <- read_best_params(seed_dir)
  out <- suppressWarnings(as.numeric(vals[ratio_parameters]))
  if (any(!is.finite(out))) {
    missing <- ratio_parameters[!is.finite(out)]
    stop("Missing/non-finite ratio parameters in ", seed_dir, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out
}

single_side_ratio_features_for_seed <- function(seed_dir, side) {
  side <- match.arg(side, c("invivo", "invitro"))
  vals <- read_best_params(seed_dir)
  raw <- suppressWarnings(as.numeric(vals[ratio_parameters]))
  if (any(!is.finite(raw))) {
    missing <- ratio_parameters[!is.finite(raw)]
    stop("Missing/non-finite ", side, " parameters in ", seed_dir, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  scales <- ratio_parameter_scales[ratio_parameters]
  log_idx <- scales == "log10"
  if (any(log_idx & raw <= 0)) {
    missing <- ratio_parameters[log_idx & raw <= 0]
    stop("Non-positive ", side, " parameters for log10 transform in ", seed_dir, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- raw
  out[log_idx] <- log10(out[log_idx])
  names(out) <- single_side_ratio_feature_cols(side)
  out
}

invivo_only_oxygen_features_for_seed <- function(seed_dir) {
  vals <- read_best_params(seed_dir)
  raw <- suppressWarnings(as.numeric(vals[invivo_only_oxygen_parameters]))
  if (any(!is.finite(raw) | raw <= 0)) {
    missing <- invivo_only_oxygen_parameters[!is.finite(raw) | raw <= 0]
    stop(
      "Missing/non-positive in vivo-only oxygen parameters in ",
      seed_dir, ": ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  out <- log10(raw)
  names(out) <- invivo_only_oxygen_feature_cols
  out
}

single_side_features_for_seed <- function(seed_dir, side) {
  side <- match.arg(side, c("invivo", "invitro"))
  out <- single_side_ratio_features_for_seed(seed_dir, side)
  if (identical(side, "invivo")) {
    out <- c(out, invivo_only_oxygen_features_for_seed(seed_dir))
  }
  out
}

build_pair_tables <- function(invivo_top, invitro_top) {
  invivo_params <- lapply(invivo_top$seed_dir, params_for_seed)
  invivo_only_features <- lapply(invivo_top$seed_dir, invivo_only_oxygen_features_for_seed)
  invitro_params <- lapply(invitro_top$seed_dir, params_for_seed)
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(invivo_top))) {
    for (j in seq_len(nrow(invitro_top))) {
      ratio_vals <- invivo_params[[i]] / invitro_params[[j]]
      log_vals <- log10(ratio_vals)
      if (any(!is.finite(log_vals))) next
      row <- data.frame(
        pair_id = sprintf("V%02d-I%02d", invivo_top$rank[[i]], invitro_top$rank[[j]]),
        invivo_rank = invivo_top$rank[[i]],
        invitro_rank = invitro_top$rank[[j]],
        invivo_seed = invivo_top$seed_int[[i]],
        invitro_seed = invitro_top$seed_int[[j]],
        invivo_seed_dir = invivo_top$seed_dir[[i]],
        invitro_seed_dir = invitro_top$seed_dir[[j]],
        invivo_objective = invivo_top$objective[[i]],
        invitro_objective = invitro_top$objective[[j]],
        stringsAsFactors = FALSE
      )
      for (p in ratio_parameters) row[[paste0("ratio_", p)]] <- ratio_vals[[match(p, ratio_parameters)]]
      for (p in ratio_parameters) row[[paste0("log10_ratio_", p)]] <- log_vals[[match(p, ratio_parameters)]]
      for (col in invivo_only_oxygen_feature_cols) row[[col]] <- invivo_only_features[[i]][[col]]
      rows[[idx]] <- row
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

build_single_side_feature_table <- function(top, side) {
  side <- match.arg(side, c("invivo", "invitro"))
  keys <- data.frame(
    rank = top$rank,
    seed = top$seed_int,
    objective = top$objective,
    seed_dir = top$seed_dir,
    stringsAsFactors = FALSE
  )
  features <- do.call(rbind, lapply(top$seed_dir, single_side_features_for_seed, side = side))
  features <- as.data.frame(features, check.names = FALSE, stringsAsFactors = FALSE)
  cbind(keys, features)
}

choose_k <- function(x, k_spec = "auto", k_max = 8L) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  n <- nrow(x)
  if (n <= 2L) return(1L)
  if (!identical(k_spec, "auto")) {
    k <- suppressWarnings(as.integer(k_spec))
    return(max(1L, min(k, n)))
  }
  max_k <- max(2L, min(k_max, n - 1L))
  d <- stats::dist(x)
  hc <- stats::hclust(d, method = "ward.D2")
  if (!requireNamespace("cluster", quietly = TRUE)) return(min(5L, max_k))
  scores <- sapply(2L:max_k, function(k) {
    cl <- stats::cutree(hc, k = k)
    mean(cluster::silhouette(cl, d)[, "sil_width"], na.rm = TRUE)
  })
  as.integer((2L:max_k)[which.max(scores)])
}

umap_embedding <- function(x, n_neighbors, umap_seed) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  if (nrow(x) <= 1L) {
    emb <- matrix(c(0, 0), nrow = 1L)
    colnames(emb) <- c("UMAP1", "UMAP2")
    return(emb)
  }
  if (nrow(x) < 3L) {
    emb <- stats::cmdscale(stats::dist(x), k = min(2L, nrow(x) - 1L))
    emb <- as.matrix(emb)
    if (ncol(emb) < 2L) emb <- cbind(emb[, 1L], 0)
    colnames(emb) <- c("UMAP1", "UMAP2")
    return(emb)
  }
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Package 'uwot' is required for UMAP.", call. = FALSE)
  set.seed(umap_seed)
  umap_args <- list(
    X = x,
    n_neighbors = min(n_neighbors, nrow(x) - 1L),
    min_dist = 0.1,
    metric = "euclidean",
    n_components = 2L,
    n_threads = 1L,
    ret_model = FALSE,
    verbose = FALSE
  )
  if ("n_sgd_threads" %in% names(formals(uwot::umap))) {
    umap_args$n_sgd_threads <- 1L
  }
  emb <- do.call(uwot::umap, umap_args)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb
}

cluster_feature_table <- function(means, feature_cols, k_spec, k_max) {
  missing <- setdiff(feature_cols, names(means))
  if (length(missing)) stop("Cluster feature columns are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x <- as.matrix(means[, feature_cols, drop = FALSE])
  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  k <- choose_k(x_scaled, k_spec = k_spec, k_max = k_max)
  cl <- if (k <= 1L) rep(1L, nrow(means)) else stats::cutree(stats::hclust(stats::dist(x_scaled), method = "ward.D2"), k = k)
  means$cluster <- as.integer(cl)
  means$cluster_space <- "profile"
  means$cluster_feature_count <- length(feature_cols)
  means$cluster_feature_set <- paste(feature_cols, collapse = ",")
  means
}

cluster_means <- function(pairs,
                          rank_col,
                          seed_col,
                          objective_col,
                          seed_dir_col,
                          k_spec,
                          k_max,
                          extra_feature_fn = NULL,
                          extra_feature_cols = character()) {
  log_cols <- ratio_feature_cols
  keys <- unique(pairs[, c(rank_col, seed_col, objective_col, seed_dir_col), drop = FALSE])
  names(keys) <- c("rank", "seed", "objective", "seed_dir")
  means <- aggregate(pairs[, log_cols], list(rank = pairs[[rank_col]]), mean)
  means <- merge(keys, means, by = "rank", all.x = TRUE, sort = FALSE)
  means <- means[order(means$rank), , drop = FALSE]
  if (!is.null(extra_feature_fn)) {
    extra <- do.call(rbind, lapply(means$seed_dir, extra_feature_fn))
    extra <- as.data.frame(extra, check.names = FALSE, stringsAsFactors = FALSE)
    missing_extra <- setdiff(extra_feature_cols, names(extra))
    if (length(missing_extra)) {
      stop("Extra feature function did not return columns: ", paste(missing_extra, collapse = ", "), call. = FALSE)
    }
    for (col in extra_feature_cols) means[[col]] <- suppressWarnings(as.numeric(extra[[col]]))
  }
  feature_cols <- c(log_cols, extra_feature_cols)
  cluster_feature_table(means, feature_cols, k_spec = k_spec, k_max = k_max)
}

representatives_from_clusters <- function(means, feature_cols = ratio_feature_cols) {
  missing <- setdiff(feature_cols, names(means))
  if (length(missing)) stop("Representative feature columns are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x_all <- scale(as.matrix(means[, feature_cols, drop = FALSE]))
  x_all[!is.finite(x_all)] <- 0
  means$.representative_row_id <- seq_len(nrow(means))
  reps <- do.call(rbind, lapply(split(means, means$cluster), function(df) {
    x <- x_all[df$.representative_row_id, , drop = FALSE]
    centroid <- colMeans(x)
    df$centroid_distance <- sqrt(rowSums(sweep(x, 2L, centroid, "-")^2))
    df <- df[order(df$centroid_distance, df$objective, df$rank), , drop = FALSE]
    df[1L, , drop = FALSE]
  }))
  reps$.representative_row_id <- NULL
  reps[order(reps$cluster), , drop = FALSE]
}

parse_rank_list <- function(spec, available_ranks) {
  spec <- trimws(spec)
  if (!nzchar(spec)) return(integer(0))
  vals <- suppressWarnings(as.integer(trimws(strsplit(spec, ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  unique(vals[vals %in% available_ranks])
}

auto_phase2_invitro_ranks <- function(invitro_means) {
  log_cols <- ratio_feature_cols
  x <- as.matrix(invitro_means[, log_cols, drop = FALSE])
  global <- colMeans(x)
  dist_global <- sqrt(rowSums(sweep(x, 2L, global, "-")^2))
  outlier_rank <- invitro_means$rank[[which.max(dist_global)]]
  middle_candidates <- invitro_means[!invitro_means$rank %in% c(1L, outlier_rank), , drop = FALSE]
  if (nrow(middle_candidates)) {
    middle_x <- as.matrix(middle_candidates[, log_cols, drop = FALSE])
    middle_rank <- middle_candidates$rank[[which.min(sqrt(rowSums(sweep(middle_x, 2L, global, "-")^2)))]]
  } else {
    middle_rank <- integer(0)
  }
  unique(c(1L, middle_rank, outlier_rank))
}

manifest_columns <- c(
  "warmup_label", "phase", "invivo_family", "invivo_rank", "invivo_seed", "invivo_seed_dir",
  "invitro_family", "invitro_rank", "invitro_seed", "invitro_seed_dir",
  "selection_reason", "joint_run_prefix", "joint_soft_coupling_parameters_table"
)

empty_manifest <- function() {
  out <- as.data.frame(setNames(rep(list(character()), length(manifest_columns)), manifest_columns), stringsAsFactors = FALSE)
  out[0L, , drop = FALSE]
}

write_seed_plan_mode <- function(out_dir, mode, invivo_top_n, invitro_top_n, warmup_pairs) {
  write_tsv(
    data.frame(
      field = c("mode", "invivo_top_n", "invitro_top_n", "warmup_pairs"),
      value = as.character(c(mode, invivo_top_n, invitro_top_n, warmup_pairs)),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "multi_warmup_seed_plan_mode.tsv")
  )
}





paired_profile_umap_coords <- function(pairs, out_dir, umap_seed) {
  feature_cols <- invivo_profile_feature_cols
  required <- c(
    "pair_id", "invivo_rank", "invitro_rank", "invivo_seed", "invitro_seed",
    "invivo_seed_dir", "invitro_seed_dir", "invivo_objective", "invitro_objective",
    feature_cols
  )
  missing <- setdiff(required, names(pairs))
  if (length(missing)) stop("Paired profile table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

  x <- as.matrix(pairs[, feature_cols, drop = FALSE])
  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  emb <- umap_embedding(x_scaled, n_neighbors = min(15L, nrow(x_scaled) - 1L), umap_seed = umap_seed)
  coords <- cbind(
    pairs[, c(
      "pair_id", "invivo_rank", "invitro_rank", "invivo_seed", "invitro_seed",
      "invivo_seed_dir", "invitro_seed_dir", "invivo_objective", "invitro_objective"
    )],
    data.frame(UMAP1 = emb[, 1L], UMAP2 = emb[, 2L], stringsAsFactors = FALSE)
  )
  write_tsv(coords, file.path(out_dir, "joint_soft_coupling_18d_profile_umap_coords.tsv"))
  coords
}



single_side_umap_coords <- function(means, side, feature_cols, out_dir, umap_seed) {
  side <- match.arg(side, c("invivo", "invitro"))
  required <- c("rank", "seed", "seed_dir", "objective", "cluster", feature_cols)
  missing <- setdiff(required, names(means))
  if (length(missing)) stop(side, " cluster table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

  x <- as.matrix(means[, feature_cols, drop = FALSE])
  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  emb <- umap_embedding(x_scaled, n_neighbors = min(15L, nrow(x_scaled) - 1L), umap_seed = umap_seed)
  coords <- cbind(
    means[, c("rank", "seed", "seed_dir", "objective", "cluster"), drop = FALSE],
    data.frame(UMAP1 = emb[, 1L], UMAP2 = emb[, 2L], stringsAsFactors = FALSE)
  )
  write_tsv(coords, file.path(out_dir, paste0("multi_warmup_", side, "_only_profile_umap_coords.tsv")))
  coords
}


run_single_side_mode <- function(run_dir, out_dir, side, top_n, k_spec, umap_seed) {
  top <- read_top_seed_table(run_dir, top_n, if (identical(side, "invivo")) "in vivo" else "in vitro")
  feature_cols <- single_side_feature_cols(side)
  means <- build_single_side_feature_table(top, side)
  means <- cluster_feature_table(means, feature_cols, k_spec = k_spec, k_max = max(1L, top_n - 1L))
  reps <- representatives_from_clusters(means, feature_cols = feature_cols)
  write_tsv(means, file.path(out_dir, paste0("multi_warmup_", side, "_cluster_summary.tsv")))
  write_tsv(reps, file.path(out_dir, paste0("multi_warmup_", side, "_representatives.tsv")))
  write_tsv(
    data.frame(feature = feature_cols, source = single_side_feature_sources(side), stringsAsFactors = FALSE),
    file.path(out_dir, paste0("multi_warmup_", side, "_feature_columns.tsv"))
  )
  coords <- single_side_umap_coords(means, side = side, feature_cols = feature_cols, out_dir = out_dir, umap_seed = umap_seed)
  write_tsv(empty_manifest(), file.path(out_dir, "multi_warmup_manifest.tsv"))
  mode <- paste0(side, "_only")
  write_seed_plan_mode(out_dir, mode = mode, invivo_top_n = if (identical(side, "invivo")) top_n else 0L, invitro_top_n = if (identical(side, "invitro")) top_n else 0L, warmup_pairs = 0L)
  message("Wrote ", side, "-only cluster summary and representatives under: ", out_dir)
  message("No warm-up pairs were generated because ", if (identical(side, "invivo")) "--invitro_top_n=0" else "--invivo_top_n=0", ".")
  invisible(NULL)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  out_dir <- normalizePath(as_chr(argv$out_dir), mustWork = FALSE)
  if (!nzchar(out_dir)) stop("--out_dir is required.", call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  top_n <- as_int(argv$top_n, 10L)
  invivo_top_n <- validate_top_n("invivo_top_n", as_int(argv$invivo_top_n %||% argv$multi_warmup_invivo_top_n, top_n))
  invitro_top_n <- validate_top_n("invitro_top_n", as_int(argv$invitro_top_n %||% argv$multi_warmup_invitro_top_n, top_n))
  if (invivo_top_n == 0L && invitro_top_n == 0L) {
    stop("At least one of --invivo_top_n or --invitro_top_n must be greater than 0.", call. = FALSE)
  }
  umap_seed <- as_int(argv$umap_seed %||% argv$multi_warmup_umap_seed, 1L)
  invivo_k <- as_chr(argv$invivo_k, "auto")
  invitro_k <- as_chr(argv$invitro_k, "auto")
  invivo_cluster_space <- tolower(as_chr(argv$invivo_cluster_space, "profile"))
  if (!identical(invivo_cluster_space, "profile")) {
    stop("--invivo_cluster_space=umap has been removed; multi-warmup representative selection now always uses the 18D profile space.", call. = FALSE)
  }
  include_phase2 <- as_bool(argv$include_phase2, FALSE)

  invivo_run_arg <- as_chr(argv$invivo_run_dir)
  invitro_run_arg <- as_chr(argv$invitro_run_dir)
  invivo_run_dir <- ""
  invitro_run_dir <- ""
  if (invivo_top_n > 0L) {
    if (!nzchar(invivo_run_arg)) stop("--invivo_run_dir is required when --invivo_top_n > 0.", call. = FALSE)
    invivo_run_dir <- normalizePath(invivo_run_arg, mustWork = TRUE)
  }
  if (invitro_top_n > 0L) {
    if (!nzchar(invitro_run_arg)) stop("--invitro_run_dir is required when --invitro_top_n > 0.", call. = FALSE)
    invitro_run_dir <- normalizePath(invitro_run_arg, mustWork = TRUE)
  }

  if (invitro_top_n == 0L) {
    return(run_single_side_mode(invivo_run_dir, out_dir, side = "invivo", top_n = invivo_top_n, k_spec = invivo_k, umap_seed = umap_seed))
  }
  if (invivo_top_n == 0L) {
    return(run_single_side_mode(invitro_run_dir, out_dir, side = "invitro", top_n = invitro_top_n, k_spec = invitro_k, umap_seed = umap_seed))
  }

  anchor_ranks <- parse_rank_list(as_chr(argv$invitro_anchor_ranks, "1"), seq_len(invitro_top_n))

  invivo_top <- read_top_seed_table(invivo_run_dir, invivo_top_n, "in vivo")
  invitro_top <- read_top_seed_table(invitro_run_dir, invitro_top_n, "in vitro")
  pairs <- build_pair_tables(invivo_top, invitro_top)

  invivo_means <- cluster_means(
    pairs, "invivo_rank", "invivo_seed", "invivo_objective", "invivo_seed_dir",
    invivo_k,
    k_max = invivo_top_n - 1L,
    extra_feature_fn = invivo_only_oxygen_features_for_seed,
    extra_feature_cols = invivo_only_oxygen_feature_cols
  )
  invitro_means <- cluster_means(
    pairs, "invitro_rank", "invitro_seed", "invitro_objective", "invitro_seed_dir",
    invitro_k, k_max = min(3L, invitro_top_n - 1L)
  )
  invivo_reps <- representatives_from_clusters(invivo_means, feature_cols = invivo_profile_feature_cols)
  write_tsv(invivo_means, file.path(out_dir, "multi_warmup_invivo_cluster_summary.tsv"))
  write_tsv(invitro_means, file.path(out_dir, "multi_warmup_invitro_cluster_summary.tsv"))
  write_tsv(invivo_reps, file.path(out_dir, "multi_warmup_invivo_representatives.tsv"))
  write_tsv(
    data.frame(
      feature = invivo_profile_feature_cols,
      source = c(rep("mean_paired_log10_invivo_over_invitro_ratio", length(ratio_feature_cols)), rep("invivo_only_log10_best_param", length(invivo_only_oxygen_feature_cols))),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "multi_warmup_invivo_feature_columns.tsv")
  )
  paired_coords <- paired_profile_umap_coords(pairs, out_dir = out_dir, umap_seed = umap_seed)

  if (!length(anchor_ranks)) anchor_ranks <- 1L
  if (include_phase2) {
    phase2_spec <- as_chr(argv$phase2_invitro_anchor_ranks, "auto")
    phase2_ranks <- if (identical(phase2_spec, "auto")) auto_phase2_invitro_ranks(invitro_means) else parse_rank_list(phase2_spec, seq_len(invitro_top_n))
    anchor_ranks <- unique(c(anchor_ranks, phase2_ranks))
  }
  anchors <- invitro_means[invitro_means$rank %in% anchor_ranks, , drop = FALSE]
  anchors <- anchors[order(match(anchors$rank, anchor_ranks)), , drop = FALSE]

  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(invivo_reps))) {
    for (j in seq_len(nrow(anchors))) {
      warmup_label <- sprintf("v%02d_i%02d", invivo_reps$rank[[i]], anchors$rank[[j]])
      rows[[idx]] <- data.frame(
        warmup_label = warmup_label,
        phase = if (anchors$rank[[j]] == 1L) "phase1" else "phase2",
        invivo_family = sprintf("C%02d", invivo_reps$cluster[[i]]),
        invivo_rank = invivo_reps$rank[[i]],
        invivo_seed = invivo_reps$seed[[i]],
        invivo_seed_dir = invivo_reps$seed_dir[[i]],
        invitro_family = sprintf("I%02d", anchors$cluster[[j]]),
        invitro_rank = anchors$rank[[j]],
        invitro_seed = anchors$seed[[j]],
        invitro_seed_dir = anchors$seed_dir[[j]],
        selection_reason = "invivo_nearest_18d_profile_cluster_centroid",
        joint_run_prefix = paste0("fit_joint_", warmup_label),
        joint_soft_coupling_parameters_table = file.path(out_dir, "joint_soft_coupling_tables", paste0("joint_soft_coupling_parameters_table__", warmup_label, ".csv")),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  manifest <- do.call(rbind, rows)
  write_tsv(manifest, file.path(out_dir, "multi_warmup_manifest.tsv"))
  write_seed_plan_mode(out_dir, mode = "paired", invivo_top_n = invivo_top_n, invitro_top_n = invitro_top_n, warmup_pairs = nrow(manifest))
  message("Wrote multi-warmup manifest: ", file.path(out_dir, "multi_warmup_manifest.tsv"))
  message("Selected ", length(unique(manifest$invivo_family)), " in vivo clusters and ", nrow(manifest), " warm-up pairs.")
}
