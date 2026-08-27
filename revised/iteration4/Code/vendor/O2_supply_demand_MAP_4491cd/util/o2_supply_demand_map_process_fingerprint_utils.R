#!/usr/bin/env Rscript

# Canonical shared process-fingerprint parsing, transformation, distance, and
# clustering helpers. Simulation and analysis source this module rather than
# carrying private copies.

.o2ipa_util_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[
    basename(frame_files) == "o2_supply_demand_map_process_fingerprint_utils.R"
  ]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(
  file.path(.o2ipa_util_dir, "o2_supply_demand_map_shared.R"),
  local = TRUE,
  chdir = TRUE
)
rm(.o2ipa_util_dir)

o2ipa_null_coalesce <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

o2ipa_numeric <- o2sd_numeric

o2ipa_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!grepl("^--", arg)) {
      i <- i + 1L
      next
    }
    kv <- sub("^--", "", arg)
    eq <- regexpr("=", kv, fixed = TRUE)
    if (eq > 0L) {
      key <- substr(kv, 1L, eq - 1L)
      val <- substr(kv, eq + 1L, nchar(kv))
      out[[key]] <- val
      i <- i + 1L
    } else {
      key <- kv
      if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      } else {
        out[[key]] <- TRUE
        i <- i + 1L
      }
    }
  }
  out
}

o2ipa_as_chr <- function(x, default = "") {
  val <- o2ipa_null_coalesce(x, default)
  val <- as.character(val[[1]])
  if (!nzchar(val)) default else val
}

o2ipa_as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric(o2ipa_null_coalesce(x, default)[[1]]))
  if (is.finite(val)) val else default
}

o2ipa_as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer(o2ipa_null_coalesce(x, default)[[1]]))
  if (!is.na(val)) val else default
}

o2ipa_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  if (is.logical(x[[1]])) return(isTRUE(x[[1]]))
  tolower(trimws(as.character(x[[1]]))) %in% c("1", "true", "t", "yes", "y", "on")
}

o2ipa_split_csv <- function(x, default = character()) {
  txt <- trimws(o2ipa_as_chr(x, paste(default, collapse = ",")))
  if (!nzchar(txt)) return(default)
  vals <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

o2ipa_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
}

o2ipa_find_workflow_root <- function(path = o2ipa_script_dir()) {
  cur <- normalizePath(path, mustWork = FALSE)
  if (file.exists(cur) && !dir.exists(cur)) cur <- dirname(cur)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(cur, "util", "o2_supply_demand_map_shared.R"))) {
      return(normalizePath(cur, mustWork = FALSE))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(file.path(path, "..", ".."), mustWork = FALSE)
}

o2ipa_workflow_root <- function(script_dir = o2ipa_script_dir()) {
  o2ipa_find_workflow_root(script_dir)
}

o2ipa_repo_root <- function(script_dir = o2ipa_script_dir()) {
  normalizePath(file.path(o2ipa_workflow_root(script_dir), "..", "..", ".."), mustWork = FALSE)
}

o2ipa_default_out_dir <- function(script_dir = o2ipa_script_dir()) {
  file.path(o2ipa_repo_root(script_dir), "oxygen", "results", "analysis")
}

o2ipa_mkdirs <- function(out_dir) {
  dirs <- file.path(out_dir, c("tables", "cache", "logs"))
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

o2ipa_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

o2ipa_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

o2ipa_read_csv_or_tsv <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  if (length(first) && grepl(",", first, fixed = TRUE) && !grepl("\t", first, fixed = TRUE)) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    o2ipa_read_tsv(path)
  }
}

o2ipa_auc <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

o2ipa_seed_number <- function(seed_id) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

o2ipa_norm_seed <- function(x) {
  s <- as.character(x)
  s <- trimws(s)
  out <- ifelse(grepl("^seed[0-9]+$", s), s, ifelse(grepl("^[0-9]+$", s), paste0("seed", s), s))
  out
}

o2ipa_order_seeds <- function(seed_id) {
  n <- o2ipa_seed_number(seed_id)
  order(ifelse(is.na(n), Inf, n), seed_id)
}

o2ipa_target_params <- function() {
  c(
    "O2_crit", "alpha_o2", "gamma_growth", "mu_hp", "p_misseg", "k_o_mis",
    "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O", "lam_max",
    "p_mis_base", "p_wgd", "gamma_mu", "o2_S0", "kappa_O", "eta_o2", "rho_2N"
  )
}

o2ipa_param_module <- function(parameter) {
  map <- c(
    O2_crit = "hypoxia_sensing", alpha_o2 = "proliferation", gamma_growth = "proliferation",
    mu_hp = "death", p_misseg = "CIN_missegregation", k_o_mis = "CIN_missegregation",
    buffer_smax = "aneuploidy_buffering", buffer_beta = "aneuploidy_buffering",
    buffer_n_exp = "aneuploidy_buffering", n_O = "hypoxia_sensing", lam_max = "proliferation",
    p_mis_base = "CIN_missegregation", p_wgd = "WGD", gamma_mu = "death",
    o2_S0 = "O2_supply_demand", kappa_O = "O2_supply_demand",
    eta_o2 = "O2_supply_demand", rho_2N = "O2_supply_demand"
  )
  unname(map[parameter])
}

o2ipa_param_aliases <- function() {
  list(
    O2_crit = c("O2_crit", "o2_crit", "O2crit", "o2crit"),
    alpha_o2 = c("alpha_o2", "alpha_O2"),
    gamma_growth = c("gamma_growth"),
    mu_hp = c("mu_hp"),
    p_misseg = c("p_misseg", "p_mis", "p_missegregation"),
    k_o_mis = c("k_o_mis", "ko_mis", "k_O_mis"),
    buffer_smax = c("buffer_smax", "buffer_s_max"),
    buffer_beta = c("buffer_beta"),
    buffer_n_exp = c("buffer_n_exp", "buffer_n"),
    n_O = c("n_O", "n_o", "nO"),
    lam_max = c("lam_max", "lambda_max"),
    p_mis_base = c("p_mis_base", "p_misseg_base"),
    p_wgd = c("p_wgd", "wgd_prob"),
    gamma_mu = c("gamma_mu"),
    o2_S0 = c("o2_S0", "O2_S0", "S0_o2"),
    kappa_O = c("kappa_O", "kappa_o"),
    eta_o2 = c("eta_o2", "eta_O2"),
    rho_2N = c("rho_2N", "rho2N")
  )
}

o2ipa_optimizer_transform <- function(parameter) {
  log10_params <- c(
    "O2_crit", "alpha_o2", "mu_hp", "p_misseg", "k_o_mis",
    "buffer_beta", "buffer_n_exp", "lam_max", "p_mis_base", "p_wgd",
    "o2_S0", "kappa_O", "eta_o2", "rho_2N"
  )
  if (parameter %in% log10_params) "log10" else "identity"
}

o2ipa_transform_parameter_value <- function(parameter, value, epsilon = 1e-12) {
  v <- as.numeric(value)
  if (!is.finite(v)) return(NA_real_)
  tr <- o2ipa_optimizer_transform(parameter)
  if (identical(tr, "log10")) {
    if (v <= 0) return(NA_real_)
    return(log10(v))
  }
  if (identical(tr, "logit")) {
    vv <- min(max(v, epsilon), 1 - epsilon)
    return(stats::qlogis(vv))
  }
  v
}

o2ipa_transform_metadata <- function(parameter_tables = list()) {
  params <- o2ipa_target_params()
  out <- data.frame(
    parameter = params,
    module = vapply(params, o2ipa_param_module, character(1)),
    raw_scale = ifelse(params %in% c("p_misseg", "p_mis_base", "p_wgd", "buffer_smax"), "bounded_or_probability", "positive_or_identity"),
    transform = vapply(params, o2ipa_optimizer_transform, character(1)),
    epsilon = ifelse(params %in% c("p_misseg", "p_mis_base", "p_wgd", "buffer_smax"), 1e-12, NA_real_),
    optimizer_scale = vapply(params, o2ipa_optimizer_transform, character(1)),
    source_file = NA_character_,
    stringsAsFactors = FALSE
  )
  for (pt in parameter_tables) {
    if (!is.data.frame(pt) || !all(c("param_prototype", "param_name") %in% names(pt))) next
    idx <- match(out$parameter, pt$param_prototype)
    hit <- which(!is.na(idx))
    if (length(hit)) {
      out$source_file[hit] <- pt$source_file[idx[hit]]
      out$optimizer_scale[hit] <- ifelse(grepl("^log10_", pt$param_name[idx[hit]]), "log10", out$optimizer_scale[hit])
    }
  }
  out
}

o2ipa_extract_param <- function(seed_id, parameter, best_vals, summary_row, matrix_row) {
  aliases <- o2ipa_param_aliases()[[parameter]]
  sources <- list()
  if (length(best_vals)) {
    sources$preferred_values <- best_vals
  }
  if (!is.null(summary_row) && nrow(summary_row) == 1L) {
    vals <- numeric(0)
    for (a in aliases) {
      col <- paste0("value__", a)
      if (col %in% names(summary_row)) vals[a] <- suppressWarnings(as.numeric(summary_row[[col]][[1]]))
    }
    sources$seed_summary_value_cols <- vals
  }
  if (!is.null(matrix_row) && nrow(matrix_row) == 1L) {
    vals <- numeric(0)
    for (a in aliases) {
      for (col in c(paste0("final_", a), a)) {
        if (col %in% names(matrix_row)) vals[a] <- suppressWarnings(as.numeric(matrix_row[[col]][[1]]))
      }
    }
    sources$param_matrix <- vals
  }

  for (src in names(sources)) {
    vals <- sources[[src]]
    for (a in aliases) {
      if (a %in% names(vals)) {
        v <- suppressWarnings(as.numeric(vals[[a]]))
        if (is.finite(v)) {
          return(list(value = v, source = src, alias = a))
        }
      }
    }
  }
  list(value = NA_real_, source = NA_character_, alias = NA_character_)
}

o2ipa_params_wide <- function(params_long, value_col = "value") {
  seeds <- unique(params_long$seed_id)
  params <- o2ipa_target_params()
  mat <- matrix(NA_real_, nrow = length(seeds), ncol = length(params), dimnames = list(seeds, params))
  for (i in seq_len(nrow(params_long))) {
    mat[params_long$seed_id[[i]], params_long$parameter[[i]]] <- params_long[[value_col]][[i]]
  }
  as.data.frame(mat, check.names = FALSE)
}

o2ipa_feature_transform <- function(value, feature_type, oxygen_floor = 1e-6) {
  v <- as.numeric(value)
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v)
  ft <- as.character(feature_type)
  for (tp in unique(ft)) {
    idx <- ok & ft == tp
    if (!any(idx)) next
    if (tp %in% c("probability", "fraction", "binary")) {
      vv <- pmin(pmax(v[idx], 1e-12), 1 - 1e-12)
      out[idx] <- stats::qlogis(vv)
    } else if (tp %in% c("positive", "ratio")) {
      out[idx] <- ifelse(v[idx] > 0, log10(v[idx]), NA_real_)
    } else if (tp %in% c("time")) {
      out[idx] <- ifelse(v[idx] >= 0, log1p(v[idx]), NA_real_)
    } else if (tp %in% c("oxygen")) {
      out[idx] <- ifelse(v[idx] >= 0, log10(pmax(v[idx], oxygen_floor)), NA_real_)
    } else if (tp %in% c("threshold")) {
      out[idx] <- ifelse(v[idx] > 0, log10(v[idx]), NA_real_)
    } else {
      out[idx] <- v[idx]
    }
  }
  out
}

o2ipa_scale_long_features <- function(long_df, missing_feature_max = 0.5, missing_seed_max = 0.05, scale_floor = sqrt(.Machine$double.eps)) {
  if (!nrow(long_df)) return(list(long = long_df, wide = data.frame(), metadata = data.frame(), missingness = data.frame(), excluded_seeds = character()))
  long_df$transformed_value <- o2ipa_feature_transform(long_df$raw_value, long_df$feature_type)
  long_df$feature_key <- paste(long_df$fingerprint_scope, long_df$module, long_df$feature, sep = "||")
  meta_base <- unique(long_df[, c("feature_key", "fingerprint_scope", "module", "feature", "feature_type")])
  seeds <- unique(long_df$seed_id)
  trans_wide <- o2ipa_long_to_wide(long_df[, c("seed_id", "feature_key", "transformed_value")], "feature_key", "transformed_value")
  trans_mat <- as.matrix(trans_wide[, setdiff(names(trans_wide), "seed_id"), drop = FALSE])
  rownames(trans_mat) <- trans_wide$seed_id
  trans_mat <- trans_mat[seeds, meta_base$feature_key, drop = FALSE]
  stats_rows <- lapply(meta_base$feature_key, function(k) {
    vals <- trans_mat[, k]
    miss <- mean(!is.finite(vals))
    finite_vals <- vals[is.finite(vals)]
    med <- stats::median(finite_vals, na.rm = TRUE)
    mad <- stats::mad(finite_vals, center = med, constant = 1, na.rm = TRUE)
    if (!is.finite(med)) med <- NA_real_
    if (!is.finite(mad)) mad <- NA_real_
    data.frame(
      feature_key = k,
      n_seed = length(seeds),
      n_observed = sum(is.finite(vals)),
      n_missing = sum(!is.finite(vals)),
      missing_fraction = miss,
      center = med,
      scale = mad,
      zero_variance = !is.finite(mad) || mad <= scale_floor,
      stringsAsFactors = FALSE
    )
  })
  stat <- do.call(rbind, stats_rows)
  meta <- merge(meta_base, stat, by = "feature_key", all.x = TRUE)
  meta$scale_floor <- scale_floor
  meta$retained_for_clustering <- meta$missing_fraction <= missing_feature_max & !meta$zero_variance
  long_df <- merge(long_df, meta[, c("feature_key", "center", "scale", "retained_for_clustering")], by = "feature_key", all.x = TRUE)
  retained_keys <- meta$feature_key[meta$retained_for_clustering %in% TRUE]
  long_df$missing_before_imputation <- !is.finite(long_df$transformed_value)
  long_df$imputed_for_clustering <- long_df$retained_for_clustering %in% TRUE & long_df$missing_before_imputation
  long_df$clustering_transformed_value <- long_df$transformed_value
  imp_idx <- long_df$imputed_for_clustering %in% TRUE
  long_df$clustering_transformed_value[imp_idx] <- long_df$center[imp_idx]
  long_df$scaled_value <- with(long_df, ifelse(retained_for_clustering %in% TRUE, (clustering_transformed_value - center) / scale, NA_real_))
  retained <- long_df[long_df$feature_key %in% retained_keys, , drop = FALSE]
  if (length(retained_keys)) {
    retained_raw_mat <- trans_mat[, retained_keys, drop = FALSE]
    seed_miss <- data.frame(
      seed_id = rownames(retained_raw_mat),
      missing_fraction_retained_features = rowMeans(!is.finite(retained_raw_mat)),
      stringsAsFactors = FALSE
    )
  } else {
    seed_miss <- data.frame(seed_id = seeds, missing_fraction_retained_features = NA_real_, stringsAsFactors = FALSE)
  }
  excluded <- seed_miss$seed_id[is.finite(seed_miss$missing_fraction_retained_features) & seed_miss$missing_fraction_retained_features > missing_seed_max]
  module_counts <- table(meta$module[meta$retained_for_clustering %in% TRUE])
  weight_map <- setNames(1 / sqrt(as.numeric(module_counts)), names(module_counts))
  long_df$module_weight <- unname(weight_map[long_df$module])
  long_df$module_weight[is.na(long_df$module_weight)] <- NA_real_
  long_df$clustering_value <- long_df$scaled_value * long_df$module_weight
  if (length(retained_keys)) {
    center_map <- setNames(meta$center, meta$feature_key)
    scale_map <- setNames(meta$scale, meta$feature_key)
    module_map <- setNames(meta$module, meta$feature_key)
    weight_key_map <- setNames(as.numeric(weight_map[module_map[retained_keys]]), retained_keys)
    retained_mat <- trans_mat[, retained_keys, drop = FALSE]
    for (k in retained_keys) {
      miss <- !is.finite(retained_mat[, k])
      retained_mat[miss, k] <- center_map[[k]]
      retained_mat[, k] <- (retained_mat[, k] - center_map[[k]]) / scale_map[[k]]
      retained_mat[, k] <- retained_mat[, k] * as.numeric(weight_key_map[[k]])
    }
    keep_seed <- !(rownames(retained_mat) %in% excluded)
    retained_mat <- retained_mat[keep_seed, , drop = FALSE]
    wide <- data.frame(seed_id = rownames(retained_mat), retained_mat, check.names = FALSE)
    rownames(wide) <- NULL
  } else {
    wide <- data.frame()
  }
  list(long = long_df, wide = wide, metadata = meta, missingness = seed_miss, excluded_seeds = excluded)
}

o2ipa_long_to_wide <- function(df, feature_col, value_col) {
  if (!nrow(df)) return(data.frame())
  seeds <- unique(df$seed_id)
  feats <- unique(df[[feature_col]])
  mat <- matrix(NA_real_, nrow = length(seeds), ncol = length(feats), dimnames = list(seeds, feats))
  for (i in seq_len(nrow(df))) {
    mat[df$seed_id[[i]], df[[feature_col]][[i]]] <- as.numeric(df[[value_col]][[i]])
  }
  out <- data.frame(seed_id = rownames(mat), mat, check.names = FALSE)
  rownames(out) <- NULL
  out
}

o2ipa_distance <- function(wide) {
  if (!is.data.frame(wide) || nrow(wide) < 2L || ncol(wide) < 3L) return(NULL)
  mat <- as.matrix(wide[, setdiff(names(wide), "seed_id"), drop = FALSE])
  rownames(mat) <- wide$seed_id
  keep <- rowSums(!is.finite(mat)) == 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 1L) return(NULL)
  d <- as.matrix(stats::dist(mat, method = "euclidean"))
  diag(d) <- 0
  d
}

o2ipa_distance_long <- function(dmat) {
  if (is.null(dmat)) return(data.frame())
  idx <- which(upper.tri(dmat), arr.ind = TRUE)
  data.frame(
    seed_i = rownames(dmat)[idx[, 1]],
    seed_j = colnames(dmat)[idx[, 2]],
    distance = dmat[idx],
    stringsAsFactors = FALSE
  )
}

o2ipa_write_distance <- function(dmat, name, out_dir) {
  if (is.null(dmat)) {
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", paste0("distance_", name, "_matrix.tsv")))
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", paste0("distance_", name, "_long.tsv")))
    return(invisible(NULL))
  }
  mat_df <- data.frame(seed_id = rownames(dmat), dmat, check.names = FALSE)
  o2ipa_write_tsv(mat_df, file.path(out_dir, "tables", paste0("distance_", name, "_matrix.tsv")))
  o2ipa_write_tsv(o2ipa_distance_long(dmat), file.path(out_dir, "tables", paste0("distance_", name, "_long.tsv")))
}

o2ipa_upper_vec <- function(dmat) {
  if (is.null(dmat)) return(numeric(0))
  dmat[upper.tri(dmat)]
}

o2ipa_distance_correlations <- function(dlist) {
  pairs <- combn(names(dlist), 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    a <- dlist[[p[[1]]]]
    b <- dlist[[p[[2]]]]
    common <- intersect(rownames(a), rownames(b))
    if (length(common) < 3L) {
      rho <- NA_real_
      n <- 0L
    } else {
      aa <- a[common, common, drop = FALSE]
      bb <- b[common, common, drop = FALSE]
      av <- o2ipa_upper_vec(aa)
      bv <- o2ipa_upper_vec(bb)
      ok <- is.finite(av) & is.finite(bv)
      n <- sum(ok)
      rho <- if (n >= 3L) suppressWarnings(stats::cor(av[ok], bv[ok], method = "spearman")) else NA_real_
    }
    data.frame(distance_a = p[[1]], distance_b = p[[2]], n_pairs = n, spearman_rho = rho, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

o2ipa_silhouette_mean <- function(dmat, clusters) {
  n <- length(clusters)
  if (n < 3L || length(unique(clusters)) < 2L) return(NA_real_)
  vals <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    same <- clusters == clusters[[i]]
    other_clusters <- setdiff(unique(clusters), clusters[[i]])
    a <- if (sum(same) > 1L) mean(dmat[i, same & (seq_len(n) != i)]) else 0
    b <- min(vapply(other_clusters, function(cl) mean(dmat[i, clusters == cl]), numeric(1)))
    vals[[i]] <- if (max(a, b) > 0) (b - a) / max(a, b) else 0
  }
  mean(vals, na.rm = TRUE)
}

o2ipa_cluster_distance <- function(dmat, n_bootstrap = 200L, random_seed = 1L, feature_wide = NULL, feature_meta = NULL) {
  if (is.null(dmat) || nrow(dmat) < 3L) {
    return(list(diagnostics = data.frame(), membership = data.frame(), medoids = data.frame(), consensus = NULL, stability = data.frame(), recommended_k = NA_integer_, reason = "fewer_than_3_seeds"))
  }
  seeds <- rownames(dmat)
  hc <- stats::hclust(stats::as.dist(dmat), method = "ward.D2")
  kmax <- min(10L, nrow(dmat) - 1L)
  ks <- 2L:kmax
  diag_rows <- lapply(ks, function(k) {
    cl <- stats::cutree(hc, k = k)
    sizes <- as.integer(table(cl))
    data.frame(k = k, mean_silhouette = o2ipa_silhouette_mean(dmat, cl), min_cluster_size = min(sizes), max_cluster_size = max(sizes), cluster_size_imbalance = max(sizes) / min(sizes), stringsAsFactors = FALSE)
  })
  diag <- do.call(rbind, diag_rows)
  eligible <- diag[diag$min_cluster_size >= 2L & is.finite(diag$mean_silhouette), , drop = FALSE]
  membership <- do.call(rbind, lapply(ks, function(k) data.frame(seed_id = seeds, k = k, cluster = stats::cutree(hc, k = k), stringsAsFactors = FALSE)))
  if (!nrow(eligible)) {
    return(list(
      diagnostics = diag,
      membership = membership,
      medoids = data.frame(),
      consensus = NULL,
      stability = data.frame(seed_id = seeds, assignment_stability = NA_real_, stringsAsFactors = FALSE),
      recommended_k = NA_integer_,
      reason = "no_candidate_without_singleton_clusters",
      hclust = hc
    ))
  }
  recommended_k <- eligible$k[which.max(eligible$mean_silhouette)]
  cl_rec <- stats::cutree(hc, k = recommended_k)
  medoids <- o2ipa_medoids(dmat, cl_rec, recommended_k)

  consensus <- matrix(NA_real_, nrow = length(seeds), ncol = length(seeds), dimnames = list(seeds, seeds))
  stability <- data.frame(seed_id = seeds, assignment_stability = NA_real_, stringsAsFactors = FALSE)
  if (!is.null(feature_wide) && n_bootstrap > 0L && ncol(feature_wide) > 2L) {
    set.seed(random_seed)
    feature_mat <- as.matrix(feature_wide[match(seeds, feature_wide$seed_id), setdiff(names(feature_wide), "seed_id"), drop = FALSE])
    co <- matrix(0, nrow = length(seeds), ncol = length(seeds), dimnames = list(seeds, seeds))
    used <- 0L
    for (b in seq_len(n_bootstrap)) {
      cols <- sample(seq_len(ncol(feature_mat)), size = ncol(feature_mat), replace = TRUE)
      boot_mat <- feature_mat[, cols, drop = FALSE]
      if (ncol(boot_mat) < 1L || any(!is.finite(boot_mat))) next
      bd <- as.matrix(stats::dist(boot_mat))
      rownames(bd) <- colnames(bd) <- seeds
      bhc <- stats::hclust(stats::as.dist(bd), method = "ward.D2")
      bcl <- stats::cutree(bhc, k = recommended_k)
      co <- co + outer(bcl, bcl, "==")
      used <- used + 1L
    }
    if (used > 0L) {
      consensus <- co / used
      diag(consensus) <- 1
      stability$assignment_stability <- vapply(seq_along(seeds), function(i) {
        same <- cl_rec == cl_rec[[i]]
        if (sum(same) <= 1L) return(1)
        mean(consensus[i, same], na.rm = TRUE)
      }, numeric(1))
      diag$mean_within_cluster_consensus <- vapply(diag$k, function(k) {
        cl <- stats::cutree(hc, k = k)
        idx <- which(upper.tri(consensus), arr.ind = TRUE)
        same <- cl[idx[, 1]] == cl[idx[, 2]]
        if (!any(same)) NA_real_ else mean(consensus[cbind(idx[same, 1], idx[same, 2])], na.rm = TRUE)
      }, numeric(1))
    }
  }
  list(diagnostics = diag, membership = membership, medoids = medoids, consensus = consensus, stability = stability, recommended_k = recommended_k, reason = "ok", hclust = hc)
}

o2ipa_medoids <- function(dmat, clusters, k) {
  rows <- lapply(sort(unique(clusters)), function(cl) {
    idx <- which(clusters == cl)
    sub <- dmat[idx, idx, drop = FALSE]
    means <- if (length(idx) == 1L) 0 else rowMeans(sub)
    seed <- rownames(sub)[which.min(means)]
    data.frame(k = k, cluster = cl, medoid_seed = seed, mean_within_cluster_distance = min(means), cluster_size = length(idx), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

o2ipa_pairwise_diagnostics <- function(d_parameter, d_static, d_static18, d_dynamic, manifest, static_membership = NULL, dynamic_membership = NULL) {
  if (is.null(d_parameter) || is.null(d_static)) return(data.frame())
  common <- Reduce(intersect, Filter(Negate(is.null), list(rownames(d_parameter), rownames(d_static), if (!is.null(d_static18)) rownames(d_static18), if (!is.null(d_dynamic)) rownames(d_dynamic))))
  if (length(common) < 2L) return(data.frame())
  dp <- d_parameter[common, common, drop = FALSE]
  ds <- d_static[common, common, drop = FALSE]
  ds18 <- if (!is.null(d_static18)) d_static18[common, common, drop = FALSE] else matrix(NA_real_, length(common), length(common), dimnames = list(common, common))
  dd <- if (!is.null(d_dynamic)) d_dynamic[common, common, drop = FALSE] else matrix(NA_real_, length(common), length(common), dimnames = list(common, common))
  pv <- o2ipa_upper_vec(dp)
  sv <- o2ipa_upper_vec(ds)
  dv <- o2ipa_upper_vec(dd)
  th <- data.frame(
    distance = c("parameter", "static_process", "dynamic_phenotype"),
    close_threshold = c(stats::quantile(pv, 0.25, na.rm = TRUE), stats::quantile(sv, 0.25, na.rm = TRUE), stats::quantile(dv, 0.25, na.rm = TRUE)),
    far_threshold = c(stats::quantile(pv, 0.75, na.rm = TRUE), stats::quantile(sv, 0.75, na.rm = TRUE), stats::quantile(dv, 0.75, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  idx <- which(upper.tri(dp), arr.ind = TRUE)
  rows <- lapply(seq_len(nrow(idx)), function(ii) {
    i <- idx[ii, 1]
    j <- idx[ii, 2]
    si <- rownames(dp)[i]
    sj <- colnames(dp)[j]
    cls <- "intermediate"
    if (dp[i, j] >= th$far_threshold[th$distance == "parameter"] && ds[i, j] <= th$close_threshold[th$distance == "static_process"]) cls <- "parameter_far_process_close"
    if (ds[i, j] >= th$far_threshold[th$distance == "static_process"] && is.finite(dd[i, j]) && dd[i, j] <= th$close_threshold[th$distance == "dynamic_phenotype"]) cls <- "process_far_dynamic_close"
    if (ds[i, j] >= th$far_threshold[th$distance == "static_process"] && is.finite(dd[i, j]) && dd[i, j] >= th$far_threshold[th$distance == "dynamic_phenotype"]) cls <- "process_far_dynamic_far"
    if (dp[i, j] <= th$close_threshold[th$distance == "parameter"] && ds[i, j] >= th$far_threshold[th$distance == "static_process"]) cls <- "parameter_close_process_far"
    oi <- manifest$objective[match(si, manifest$seed_id)]
    oj <- manifest$objective[match(sj, manifest$seed_id)]
    data.frame(seed_i = si, seed_j = sj, parameter_distance = dp[i, j], static_process_distance = ds[i, j], static_18only_distance = ds18[i, j], dynamic_distance = dd[i, j], objective_i = oi, objective_j = oj, max_delta_objective = max(oi, oj, na.rm = TRUE) - min(manifest$objective, na.rm = TRUE), same_static_cluster = NA, same_dynamic_cluster = NA, diagnostic_class = cls, stringsAsFactors = FALSE)
  })
  list(pairs = do.call(rbind, rows), thresholds = th)
}
