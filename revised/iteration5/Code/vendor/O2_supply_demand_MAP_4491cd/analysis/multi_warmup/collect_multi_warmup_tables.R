#!/usr/bin/env Rscript

# Table-only collection and statistical diagnostics for completed multi-warmup
# fits. Plotting is owned by vis/multi_warmup/ and integrated cross-stage
# orchestration by runner/multi_warmup/.

.o2mw_collect_analysis_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "collect_multi_warmup_tables.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_collect_workflow_root <- normalizePath(
  file.path(.o2mw_collect_analysis_dir, "..", ".."),
  mustWork = TRUE
)
source(
  file.path(.o2mw_collect_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"),
  local = environment()
)

invivo_only_oxygen_parameters <- c("o2_S0", "kappa_O", "eta_o2", "rho_2N")
invivo_only_oxygen_feature_cols <- paste0("log10_invivo_", invivo_only_oxygen_parameters)


read_param_value_map <- function(seed_dir) {
  tab <- read_table_optional(file.path(seed_dir, "best_params.tsv"))
  if (is.null(tab) || !all(c("parameter", "value") %in% names(tab))) return(NULL)
  vals <- num_col(tab$value)
  names(vals) <- as.character(tab$parameter)
  vals
}

choose_k <- function(x, k_max = 6L) {
  n <- nrow(x)
  if (n <= 2L) return(1L)
  d <- stats::dist(x)
  hc <- stats::hclust(d, method = "ward.D2")
  max_k <- max(2L, min(k_max, n - 1L))
  if (!requireNamespace("cluster", quietly = TRUE)) return(min(2L, max_k))
  scores <- sapply(2L:max_k, function(k) {
    cl <- stats::cutree(hc, k = k)
    mean(cluster::silhouette(cl, d)[, "sil_width"], na.rm = TRUE)
  })
  as.integer((2L:max_k)[which.max(scores)])
}

best_seed_row <- function(seed_summary) {
  if (is.null(seed_summary) || !nrow(seed_summary) || !("objective" %in% names(seed_summary))) return(NULL)
  seed_summary$objective <- num_col(seed_summary$objective)
  seed_summary <- seed_summary[is.finite(seed_summary$objective), , drop = FALSE]
  if (!nrow(seed_summary)) return(NULL)
  seed_summary[order(seed_summary$objective, seed_summary$seed), , drop = FALSE][1L, , drop = FALSE]
}

collect_parameter_ratios <- function(soft_tab, best_seed, warmup_label) {
  if (is.null(soft_tab) || !nrow(soft_tab)) return(NULL)
  if (!all(c("seed", "parameter") %in% names(soft_tab))) return(NULL)
  soft_tab <- soft_tab[as.character(soft_tab$seed) == as.character(best_seed), , drop = FALSE]
  if (!nrow(soft_tab)) return(NULL)
  keep <- intersect(
    c(
      "parameter", "ratio_vivo_to_vitro", "log10_ratio_vivo_to_vitro",
      "fold_change_vivo_to_vitro", "penalty_paid", "vivo_natural", "vitro_natural",
      "vivo_transformed", "vitro_transformed", "regularization_sigma",
      "penalty_type", "welsch_c", "welsch_transition_delta", "abs_delta_over_sigma",
      "welsch_saturation_fraction", "penalty_region"
    ),
    names(soft_tab)
  )
  out <- soft_tab[, keep, drop = FALSE]
  out$warmup_label <- warmup_label
  out$best_joint_seed <- best_seed
  out
}

add_final_log10_ratio <- function(ratio_long) {
  if (is.null(ratio_long) || !is.data.frame(ratio_long) || !nrow(ratio_long)) return(ratio_long)
  log_ratio <- if ("log10_ratio_vivo_to_vitro" %in% names(ratio_long)) {
    num_col(ratio_long$log10_ratio_vivo_to_vitro)
  } else {
    rep(NA_real_, nrow(ratio_long))
  }
  if ("ratio_vivo_to_vitro" %in% names(ratio_long)) {
    natural_ratio <- num_col(ratio_long$ratio_vivo_to_vitro)
    fill <- !is.finite(log_ratio) & is.finite(natural_ratio) & natural_ratio > 0
    log_ratio[fill] <- log10(natural_ratio[fill])
  }
  if ("fold_change_vivo_to_vitro" %in% names(ratio_long)) {
    fold_change <- num_col(ratio_long$fold_change_vivo_to_vitro)
    fill <- !is.finite(log_ratio) & is.finite(fold_change) & fold_change > 0
    log_ratio[fill] <- log10(fold_change[fill])
  }
  ratio_long$final_log10_ratio_vivo_to_vitro <- log_ratio
  ratio_long
}


collect_deoptim_iterations <- function(seed_summary, manifest_row, run_dir) {
  if (is.null(seed_summary) || !is.data.frame(seed_summary) || !nrow(seed_summary)) return(NULL)
  if (!("optimizer_iter_completed" %in% names(seed_summary))) return(NULL)
  iter_completed <- num_col(seed_summary$optimizer_iter_completed)
  keep <- is.finite(iter_completed)
  if (!any(keep)) return(NULL)
  seed <- if ("seed" %in% names(seed_summary)) as.character(seed_summary$seed) else as.character(seq_len(nrow(seed_summary)))
  out <- data.frame(
    warmup_label = as.character(manifest_row$warmup_label[[1]]),
    invivo_family = if ("invivo_family" %in% names(manifest_row)) as.character(manifest_row$invivo_family[[1]]) else NA_character_,
    invitro_family = if ("invitro_family" %in% names(manifest_row)) as.character(manifest_row$invitro_family[[1]]) else NA_character_,
    joint_run_prefix = if ("joint_run_prefix" %in% names(manifest_row)) as.character(manifest_row$joint_run_prefix[[1]]) else basename(run_dir),
    joint_run_dir = run_dir,
    seed = seed,
    optimizer_iter_completed = iter_completed,
    stringsAsFactors = FALSE
  )
  optional_cols <- intersect(
    c(
      "optimizer_iter_target", "deoptim_stop_reason", "optimizer_deoptim_objective",
      "optimizer_local_accepted", "optimizer_local_convergence", "objective"
    ),
    names(seed_summary)
  )
  for (col in optional_cols) {
    out[[col]] <- seed_summary[[col]]
  }
  out[keep, , drop = FALSE]
}


collect_invivo_only_warmstart_parameters <- function(manifest, root) {
  if (is.null(manifest) || !is.data.frame(manifest) || !nrow(manifest)) return(data.frame())
  if (!("invivo_seed_dir" %in% names(manifest))) return(data.frame())
  key_cols <- intersect(c("invivo_rank", "invivo_seed", "invivo_seed_dir"), names(manifest))
  key <- if (length(key_cols)) do.call(paste, c(manifest[, key_cols, drop = FALSE], sep = "\r")) else seq_len(nrow(manifest))
  used <- manifest[!duplicated(key), , drop = FALSE]
  profile <- read_table_optional(file.path(root, "multi_warmup_invivo_representatives.tsv"))
  if (is.null(profile) || !nrow(profile)) {
    profile <- read_table_optional(file.path(root, "multi_warmup_invivo_cluster_summary.tsv"))
  }

  rows <- list()
  for (i in seq_len(nrow(used))) {
    invivo_rank <- if ("invivo_rank" %in% names(used)) suppressWarnings(as.integer(used$invivo_rank[[i]])) else NA_integer_
    invivo_seed <- if ("invivo_seed" %in% names(used)) suppressWarnings(as.integer(used$invivo_seed[[i]])) else NA_integer_
    invivo_family <- if ("invivo_family" %in% names(used)) as.character(used$invivo_family[[i]]) else NA_character_
    seed_dir <- as.character(used$invivo_seed_dir[[i]])
    seed_label <- if (is.finite(invivo_rank) && is.finite(invivo_seed)) {
      sprintf("V%02d seed%s", invivo_rank, invivo_seed)
    } else if (is.finite(invivo_rank)) {
      sprintf("V%02d", invivo_rank)
    } else if (is.finite(invivo_seed)) {
      sprintf("seed%s", invivo_seed)
    } else {
      basename(seed_dir)
    }

    log_values <- rep(NA_real_, length(invivo_only_oxygen_parameters))
    names(log_values) <- invivo_only_oxygen_parameters
    if (!is.null(profile) && nrow(profile)) {
      match_idx <- rep(FALSE, nrow(profile))
      if (is.finite(invivo_rank) && "rank" %in% names(profile)) {
        match_idx <- match_idx | (suppressWarnings(as.integer(profile$rank)) == invivo_rank)
      }
      if (is.finite(invivo_seed) && "seed" %in% names(profile)) {
        match_idx <- match_idx | (suppressWarnings(as.integer(profile$seed)) == invivo_seed)
      }
      if (any(match_idx)) {
        prof_row <- profile[which(match_idx)[1L], , drop = FALSE]
        for (j in seq_along(invivo_only_oxygen_parameters)) {
          col <- invivo_only_oxygen_feature_cols[[j]]
          if (col %in% names(prof_row)) log_values[[j]] <- num_col(prof_row[[col]])[[1]]
        }
      }
    }
    if (any(!is.finite(log_values))) {
      vals <- read_param_value_map(seed_dir)
      if (!is.null(vals)) {
        for (pname in invivo_only_oxygen_parameters[!is.finite(log_values)]) {
          value <- suppressWarnings(as.numeric(vals[[pname]]))
          if (is.finite(value) && value > 0) log_values[[pname]] <- log10(value)
        }
      }
    }

    for (pname in invivo_only_oxygen_parameters) {
      log_value <- log_values[[pname]]
      value <- if (is.finite(log_value)) 10^log_value else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        invivo_label = seed_label,
        invivo_rank = invivo_rank,
        invivo_seed = invivo_seed,
        invivo_family = invivo_family,
        invivo_seed_dir = seed_dir,
        parameter = pname,
        value = value,
        log10_value = log_value,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(out)) return(out)
  out$parameter <- factor(out$parameter, levels = invivo_only_oxygen_parameters)
  out$invivo_label <- factor(out$invivo_label, levels = rev(unique(out$invivo_label)))
  out$scaled_log10_value <- ave(out$log10_value, out$parameter, FUN = function(z) {
    z <- suppressWarnings(as.numeric(z))
    if (sum(is.finite(z)) <= 1L) return(rep(0, length(z)))
    s <- stats::sd(z, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(z)))
    (z - mean(z, na.rm = TRUE)) / s
  })
  out$display_value <- format_heatmap_value(out$value)
  out$label_color <- ifelse(abs(out$scaled_log10_value) >= 1.25, "white", "#1b2a38")
  out
}



assign_final_basins <- function(ratio_long, out_dir) {
  if (is.null(ratio_long) || !nrow(ratio_long)) return(NULL)
  ratio_long <- add_final_log10_ratio(ratio_long)
  ratio_long <- ratio_long[is.finite(ratio_long$final_log10_ratio_vivo_to_vitro), , drop = FALSE]
  if (!nrow(ratio_long)) return(NULL)
  wide <- reshape(
    ratio_long[, c("warmup_label", "parameter", "final_log10_ratio_vivo_to_vitro")],
    idvar = "warmup_label",
    timevar = "parameter",
    direction = "wide"
  )
  labels <- wide$warmup_label
  x <- as.matrix(wide[, setdiff(names(wide), "warmup_label"), drop = FALSE])
  complete <- stats::complete.cases(x)
  x <- x[complete, , drop = FALSE]
  labels <- labels[complete]
  if (!nrow(x)) return(NULL)
  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  k <- choose_k(x_scaled)
  basin <- if (k <= 1L || nrow(x_scaled) <= 2L) rep(1L, nrow(x_scaled)) else stats::cutree(stats::hclust(stats::dist(x_scaled), method = "ward.D2"), k = k)
  out <- data.frame(warmup_label = labels, final_basin_id = paste0("B", sprintf("%02d", basin)), stringsAsFactors = FALSE)
  out
}

profile_feature_cols <- c(paste0("log10_ratio_", c(
  "O2_crit", "alpha_o2", "mu_hp", "p_misseg", "k_o_mis",
  "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O",
  "gamma_growth", "lam_max", "p_mis_base", "p_wgd", "gamma_mu"
)), invivo_only_oxygen_feature_cols)
profile_ratio_parameters <- sub("^log10_ratio_", "", profile_feature_cols[grepl("^log10_ratio_", profile_feature_cols)])

logical_col <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

single_finite_value <- o2sd_first_finite_numeric

read_warmup_profile_from_seed <- function(seed_dir) {
  tab <- read_table_optional(file.path(seed_dir, "joint_soft_coupling_initial_values.tsv"))
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) {
    tab <- read_table_optional(file.path(seed_dir, "joint_warmup_initial_values.tsv"))
  }
  if (is.null(tab) || !is.data.frame(tab) || !nrow(tab)) return(NULL)
  required <- c("param_name", "parameter", "is_delta", "scale", "optimizer_value")
  if (!all(required %in% names(tab)) && all(c("optimizer_name", "warmup_init") %in% names(tab))) {
    optimizer_name <- as.character(tab$optimizer_name)
    base_name <- sub("^delta__", "", optimizer_name)
    parameter <- sub("^log10_", "", base_name)
    tab <- data.frame(
      param_name = optimizer_name,
      optimizer_name = optimizer_name,
      parameter = parameter,
      is_delta = startsWith(optimizer_name, "delta__"),
      scale = ifelse(startsWith(base_name, "log10_"), "log10", "identity"),
      optimizer_value = num_col(tab$warmup_init),
      stringsAsFactors = FALSE
    )
  }
  if (!all(required %in% names(tab))) return(NULL)
  is_delta <- logical_col(tab$is_delta)
  out <- rep(NA_real_, length(profile_feature_cols))
  names(out) <- profile_feature_cols

  for (pname in profile_ratio_parameters) {
    delta_row <- tab[is_delta & as.character(tab$parameter) == pname, , drop = FALSE]
    center_row <- tab[!is_delta & as.character(tab$parameter) == pname, , drop = FALSE]
    if (!nrow(delta_row) || !nrow(center_row)) next
    delta <- single_finite_value(delta_row$optimizer_value)
    center <- single_finite_value(center_row$optimizer_value)
    scale <- as.character(delta_row$scale[[1]])
    feature_name <- paste0("log10_ratio_", pname)
    if (!is.finite(delta)) next
    if (identical(scale, "log10")) {
      out[[feature_name]] <- delta
    } else if (is.finite(center)) {
      vivo <- center + delta / 2
      vitro <- center - delta / 2
      if (is.finite(vivo) && is.finite(vitro) && vivo > 0 && vitro > 0) {
        out[[feature_name]] <- log10(vivo / vitro)
      }
    }
  }

  for (pname in invivo_only_oxygen_parameters) {
    transformed_name <- paste0("log10_", pname)
    feature_name <- paste0("log10_invivo_", pname)
    row <- tab[!is_delta & as.character(tab$param_name) == transformed_name, , drop = FALSE]
    if (!nrow(row) && "optimizer_name" %in% names(tab)) {
      row <- tab[!is_delta & as.character(tab$optimizer_name) == transformed_name, , drop = FALSE]
    }
    if (nrow(row)) out[[feature_name]] <- single_finite_value(row$optimizer_value)
  }

  if (any(!is.finite(out[invivo_only_oxygen_feature_cols]))) {
    warmup_tab <- read_table_optional(file.path(seed_dir, "joint_warmup_initial_values.tsv"))
    if (!is.null(warmup_tab) && all(c("optimizer_name", "warmup_init") %in% names(warmup_tab))) {
      for (pname in invivo_only_oxygen_parameters) {
        feature_name <- paste0("log10_invivo_", pname)
        if (is.finite(out[[feature_name]])) next
        transformed_name <- paste0("log10_", pname)
        row <- warmup_tab[as.character(warmup_tab$optimizer_name) == transformed_name, , drop = FALSE]
        if (nrow(row)) out[[feature_name]] <- single_finite_value(row$warmup_init)
      }
    }
  }

  out
}

read_final_profile_from_seed <- function(seed_dir) {
  out <- rep(NA_real_, length(profile_feature_cols))
  names(out) <- profile_feature_cols

  soft <- read_table_optional(file.path(seed_dir, "joint_soft_coupling.tsv"))
  if (!is.null(soft) && is.data.frame(soft) && nrow(soft) && "parameter" %in% names(soft)) {
    for (pname in profile_ratio_parameters) {
      row <- soft[as.character(soft$parameter) == pname, , drop = FALSE]
      if (!nrow(row)) next
      val <- if ("log10_ratio_vivo_to_vitro" %in% names(row)) single_finite_value(row$log10_ratio_vivo_to_vitro) else NA_real_
      if (!is.finite(val) && "ratio_vivo_to_vitro" %in% names(row)) {
        ratio <- single_finite_value(row$ratio_vivo_to_vitro)
        if (is.finite(ratio) && ratio > 0) val <- log10(ratio)
      }
      if (!is.finite(val) && "fold_change_vivo_to_vitro" %in% names(row)) {
        ratio <- single_finite_value(row$fold_change_vivo_to_vitro)
        if (is.finite(ratio) && ratio > 0) val <- log10(ratio)
      }
      out[[paste0("log10_ratio_", pname)]] <- val
    }
  }

  transformed <- read_table_optional(file.path(seed_dir, "best_params_transformed.tsv"))
  transformed_vals <- NULL
  if (!is.null(transformed) && all(c("transformed_parameter", "transformed_value") %in% names(transformed))) {
    transformed_vals <- num_col(transformed$transformed_value)
    names(transformed_vals) <- as.character(transformed$transformed_parameter)
  }
  raw_vals <- NULL
  if (is.null(transformed_vals)) raw_vals <- read_param_value_map(seed_dir)
  for (pname in invivo_only_oxygen_parameters) {
    transformed_name <- paste0("log10_", pname)
    feature_name <- paste0("log10_invivo_", pname)
    val <- if (!is.null(transformed_vals)) suppressWarnings(as.numeric(transformed_vals[[transformed_name]])) else NA_real_
    if (!is.finite(val) && !is.null(raw_vals)) {
      raw <- suppressWarnings(as.numeric(raw_vals[[pname]]))
      if (is.finite(raw) && raw > 0) val <- log10(raw)
    }
    out[[feature_name]] <- val
  }

  out
}

scale_profile_matrices <- function(warmup_mat, final_mat) {
  combined <- rbind(warmup_mat, final_mat)
  scaled <- combined
  for (j in seq_len(ncol(combined))) {
    z <- suppressWarnings(as.numeric(combined[, j]))
    finite <- is.finite(z)
    if (!any(finite)) {
      scaled[, j] <- 0
      next
    }
    mu <- mean(z[finite])
    sdv <- stats::sd(z[finite])
    if (!is.finite(sdv) || sdv == 0) {
      scaled[, j] <- 0
    } else {
      scaled[, j] <- (z - mu) / sdv
      scaled[!finite, j] <- 0
    }
  }
  list(
    warmup = scaled[seq_len(nrow(warmup_mat)), , drop = FALSE],
    final = scaled[seq_len(nrow(final_mat)) + nrow(warmup_mat), , drop = FALSE]
  )
}

profile_distance_matrix <- function(a, b) {
  out <- matrix(NA_real_, nrow = nrow(a), ncol = nrow(b), dimnames = list(rownames(a), rownames(b)))
  for (j in seq_len(nrow(b))) {
    diff <- sweep(a, 2L, b[j, ], "-")
    out[, j] <- sqrt(rowMeans(diff^2))
  }
  out
}

assign_seed_final_basins <- function(final_scaled, objective) {
  if (!nrow(final_scaled)) return(character(0))
  if (nrow(final_scaled) <= 2L) {
    return(rep("B01", nrow(final_scaled)))
  }
  k <- choose_k(final_scaled, k_max = min(8L, nrow(final_scaled) - 1L))
  if (k <= 1L) {
    basin_num <- rep(1L, nrow(final_scaled))
  } else {
    basin_num <- stats::cutree(stats::hclust(stats::dist(final_scaled), method = "ward.D2"), k = k)
  }
  obj <- suppressWarnings(as.numeric(objective))
  cluster_order <- tapply(seq_along(basin_num), basin_num, function(idx) {
    vals <- obj[idx]
    if (any(is.finite(vals))) min(vals[is.finite(vals)]) else Inf
  })
  cluster_order <- names(sort(cluster_order, na.last = TRUE))
  label_map <- setNames(paste0("B", sprintf("%02d", seq_along(cluster_order))), cluster_order)
  unname(label_map[as.character(basin_num)])
}


collect_initial_final_distance_diagnostics <- function(manifest, root, summary_df, out_dir) {
  if (is.null(manifest) || !is.data.frame(manifest) || !nrow(manifest) || !("joint_run_prefix" %in% names(manifest))) {
    return(invisible(NULL))
  }
  warmup_profiles <- list()
  final_profiles <- list()
  meta_rows <- list()

  for (i in seq_len(nrow(manifest))) {
    warmup_label <- as.character(manifest$warmup_label[[i]])
    run_dir <- file.path(root, manifest$joint_run_prefix[[i]])
    seed_dirs <- find_valid_seed_dirs(run_dir)
    if (!length(seed_dirs)) next
    seed_dirs <- seed_dirs[order(seed_order_key(basename(seed_dirs)), basename(seed_dirs))]
    warmup_profile <- read_warmup_profile_from_seed(seed_dirs[[1]])
    if (!is.null(warmup_profile)) warmup_profiles[[warmup_label]] <- warmup_profile

    seed_summary <- read_table_optional(file.path(run_dir, "extra_results", "seed_summary.tsv"))
    for (seed_dir in seed_dirs) {
      seed <- basename(seed_dir)
      final_profile <- read_final_profile_from_seed(seed_dir)
      if (is.null(final_profile) || !any(is.finite(final_profile))) next
      seed_key <- paste0(warmup_label, "__", seed)
      final_profiles[[seed_key]] <- final_profile
      ss <- NULL
      if (!is.null(seed_summary) && is.data.frame(seed_summary) && nrow(seed_summary) && "seed" %in% names(seed_summary)) {
        hit <- as.character(seed_summary$seed) == seed
        if (any(hit)) ss <- seed_summary[which(hit)[1L], , drop = FALSE]
      }
      objective <- if (!is.null(ss) && "objective" %in% names(ss)) single_finite_value(ss$objective) else NA_real_
      meta_rows[[length(meta_rows) + 1L]] <- data.frame(
        seed_key = seed_key,
        warmup_label = warmup_label,
        invivo_family = if ("invivo_family" %in% names(manifest)) as.character(manifest$invivo_family[[i]]) else NA_character_,
        invitro_family = if ("invitro_family" %in% names(manifest)) as.character(manifest$invitro_family[[i]]) else NA_character_,
        joint_seed = seed,
        joint_run_prefix = as.character(manifest$joint_run_prefix[[i]]),
        objective = objective,
        objective_invivo = if (!is.null(ss) && "objective_invivo" %in% names(ss)) single_finite_value(ss$objective_invivo) else NA_real_,
        objective_invitro = if (!is.null(ss) && "objective_invitro" %in% names(ss)) single_finite_value(ss$objective_invitro) else NA_real_,
        objective_soft_coupling = if (!is.null(ss) && "objective_soft_coupling" %in% names(ss)) single_finite_value(ss$objective_soft_coupling) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(warmup_profiles) || !length(final_profiles) || !length(meta_rows)) return(invisible(NULL))
  warmup_mat <- do.call(rbind, warmup_profiles)
  final_mat <- do.call(rbind, final_profiles)
  warmup_mat <- warmup_mat[, profile_feature_cols, drop = FALSE]
  final_mat <- final_mat[, profile_feature_cols, drop = FALSE]
  scaled <- scale_profile_matrices(warmup_mat, final_mat)
  distances <- profile_distance_matrix(scaled$final, scaled$warmup)
  meta <- do.call(rbind, meta_rows)
  meta <- meta[match(rownames(distances), meta$seed_key), , drop = FALSE]
  meta$final_basin_id <- assign_seed_final_basins(scaled$final, meta$objective)

  own_idx <- match(meta$warmup_label, colnames(distances))
  distance_to_own <- rep(NA_real_, nrow(meta))
  has_own <- is.finite(own_idx)
  distance_to_own[has_own] <- distances[cbind(which(has_own), own_idx[has_own])]
  ordered_idx <- t(apply(distances, 1L, order))
  nearest_idx <- ordered_idx[, 1L]
  second_idx <- if (ncol(ordered_idx) >= 2L) ordered_idx[, 2L] else rep(NA_integer_, nrow(ordered_idx))
  meta$distance_to_own_warmup <- distance_to_own
  meta$nearest_warmup_label <- colnames(distances)[nearest_idx]
  meta$distance_to_nearest_warmup <- distances[cbind(seq_len(nrow(distances)), nearest_idx)]
  meta$distance_to_second_nearest_warmup <- ifelse(
    is.finite(second_idx),
    distances[cbind(seq_len(nrow(distances)), second_idx)],
    NA_real_
  )
  meta$is_nearest_own_warmup <- as.character(meta$nearest_warmup_label) == as.character(meta$warmup_label)
  meta$warmup_order <- match(meta$warmup_label, rownames(warmup_mat))
  meta <- meta[order(meta$objective, meta$warmup_label, seed_order_key(meta$joint_seed), meta$joint_seed), , drop = FALSE]

  warmup_pair_dist <- if (nrow(scaled$warmup) > 1L) {
    wd <- as.matrix(stats::dist(scaled$warmup))
    wd[upper.tri(wd)]
  } else {
    numeric(0)
  }
  warmup_dist_q25 <- if (length(warmup_pair_dist)) as.numeric(stats::quantile(warmup_pair_dist, 0.25, na.rm = TRUE)) else NA_real_
  warmup_dist_median <- if (length(warmup_pair_dist)) stats::median(warmup_pair_dist, na.rm = TRUE) else NA_real_

  basin_rows <- list()
  for (basin in unique(meta$final_basin_id)) {
    sub <- meta[as.character(meta$final_basin_id) == basin, , drop = FALSE]
    sub <- sub[order(sub$objective, sub$warmup_label, seed_order_key(sub$joint_seed), sub$joint_seed), , drop = FALSE]
    rep <- sub[1L, , drop = FALSE]
    basin_rows[[length(basin_rows) + 1L]] <- data.frame(
      final_basin_id = basin,
      n_seeds = nrow(sub),
      best_warmup_label = rep$warmup_label,
      best_joint_seed = rep$joint_seed,
      best_objective = rep$objective,
      distance_to_own_warmup = rep$distance_to_own_warmup,
      nearest_warmup_label = rep$nearest_warmup_label,
      distance_to_nearest_warmup = rep$distance_to_nearest_warmup,
      distance_to_second_nearest_warmup = rep$distance_to_second_nearest_warmup,
      warmup_distance_reference_q25 = warmup_dist_q25,
      warmup_distance_reference_median = warmup_dist_median,
      is_new_space_candidate = is.finite(rep$distance_to_nearest_warmup) &&
        is.finite(warmup_dist_q25) &&
        rep$distance_to_nearest_warmup > warmup_dist_q25,
      stringsAsFactors = FALSE
    )
  }
  basin_summary <- if (length(basin_rows)) do.call(rbind, basin_rows) else data.frame()
  if (nrow(basin_summary)) basin_summary <- basin_summary[order(basin_summary$best_objective), , drop = FALSE]

  heat_rows <- list()
  if (nrow(basin_summary)) {
    for (i in seq_len(nrow(basin_summary))) {
      seed_key <- paste0(basin_summary$best_warmup_label[[i]], "__", basin_summary$best_joint_seed[[i]])
      if (!(seed_key %in% rownames(distances))) next
      row_label <- paste0(
        basin_summary$final_basin_id[[i]],
        " / ", basin_summary$best_warmup_label[[i]],
        " / ", basin_summary$best_joint_seed[[i]]
      )
      for (warmup_label in colnames(distances)) {
        heat_rows[[length(heat_rows) + 1L]] <- data.frame(
          final_basin_id = basin_summary$final_basin_id[[i]],
          final_representative = row_label,
          representative_seed_key = seed_key,
          warmup_label = warmup_label,
          distance = distances[seed_key, warmup_label],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  heatmap_long <- if (length(heat_rows)) do.call(rbind, heat_rows) else data.frame()

  full_distance <- as.data.frame(distances, check.names = FALSE)
  full_distance$seed_key <- rownames(distances)
  full_distance <- full_distance[, c("seed_key", colnames(distances)), drop = FALSE]

  write_tsv(meta, file.path(out_dir, "multi_warmup_initial_final_distance_per_seed.tsv"))
  write_tsv(basin_summary, file.path(out_dir, "multi_warmup_final_basin_distance_summary.tsv"))
  write_tsv(heatmap_long, file.path(out_dir, "multi_warmup_final_to_warmup_distance_heatmap.tsv"))
  write_tsv(full_distance, file.path(out_dir, "multi_warmup_final_to_warmup_distance_matrix.tsv"))
  invisible(list(per_seed = meta, basin_summary = basin_summary, heatmap = heatmap_long))
}



find_valid_seed_dirs <- function(run_dir) {
  if (!dir.exists(run_dir)) return(character(0))
  dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  keep <- vapply(
    dirs,
    function(d) file.exists(file.path(d, "fit_summary.tsv")) && file.exists(file.path(d, "best_params.tsv")),
    logical(1)
  )
  dirs[keep]
}






main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root <- normalizePath(as_chr(argv$multi_warmup_root, as_chr(argv$out_dir)), mustWork = TRUE)
  manifest_path <- normalizePath(as_chr(argv$manifest, file.path(root, "multi_warmup_manifest.tsv")), mustWork = TRUE)
  out_dir <- normalizePath(as_chr(argv$out_dir, root), mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- utils::read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
  rows <- list()
  ratio_rows <- list()
  iteration_rows <- list()
  for (i in seq_len(nrow(manifest))) {
    run_dir <- file.path(root, manifest$joint_run_prefix[[i]])
    seed_summary <- read_table_optional(file.path(run_dir, "extra_results", "seed_summary.tsv"))
    iter_tab <- collect_deoptim_iterations(seed_summary, manifest[i, , drop = FALSE], run_dir)
    if (!is.null(iter_tab)) iteration_rows[[length(iteration_rows) + 1L]] <- iter_tab
    best <- best_seed_row(seed_summary)
    if (is.null(best)) next
    base_cols <- intersect(
      c(
        "seed", "objective", "objective_invivo", "objective_invitro", "objective_soft_coupling",
        "objective_burden", "objective_ploidy", "objective_necrosis",
        "boundary_penalty_active", "min_rel_dist_active", "pred1000_2N", "pred1000_4N"
      ),
      names(best)
    )
    row <- cbind(
      manifest[i, , drop = FALSE],
      data.frame(joint_run_dir = run_dir, stringsAsFactors = FALSE),
      best[, base_cols, drop = FALSE]
    )
    names(row)[names(row) == "seed"] <- "best_joint_seed"
    rows[[length(rows) + 1L]] <- row
    soft <- read_table_optional(file.path(run_dir, "extra_results", "joint_soft_coupling_all.tsv"))
    ratio_tab <- collect_parameter_ratios(soft, best$seed[[1]], manifest$warmup_label[[i]])
    if (!is.null(ratio_tab)) ratio_rows[[length(ratio_rows) + 1L]] <- ratio_tab
  }
  summary_df <- if (length(rows)) do.call(rbind, rows) else data.frame()
  ratio_long <- if (length(ratio_rows)) do.call(rbind, ratio_rows) else data.frame()
  iteration_long <- if (length(iteration_rows)) do.call(rbind, iteration_rows) else data.frame()
  invivo_only_long <- collect_invivo_only_warmstart_parameters(manifest, root)
  ratio_long <- add_final_log10_ratio(ratio_long)
  basin_df <- assign_final_basins(ratio_long, out_dir)
  if (!is.null(basin_df) && nrow(summary_df)) {
    summary_df <- merge(summary_df, basin_df, by = "warmup_label", all.x = TRUE, sort = FALSE)
  }

  if (nrow(summary_df)) {
    for (col in intersect(c("objective", "objective_invivo", "objective_invitro", "objective_soft_coupling", "objective_burden", "objective_ploidy", "objective_necrosis", "boundary_penalty_active", "min_rel_dist_active"), names(summary_df))) {
      summary_df[[col]] <- num_col(summary_df[[col]])
    }
  }
  write_tsv(summary_df, file.path(out_dir, "multi_warmup_best_seed_summary.tsv"))
  write_tsv(ratio_long, file.path(out_dir, "multi_warmup_parameter_ratio_long.tsv"))
  write_tsv(iteration_long, file.path(out_dir, "multi_warmup_deoptim_iteration_summary.tsv"))
  if (nrow(iteration_long)) {
    iteration_long$optimizer_iter_completed <- num_col(iteration_long$optimizer_iter_completed)
    finite_iterations <- iteration_long[is.finite(iteration_long$optimizer_iter_completed), , drop = FALSE]
    iteration_counts <- if (nrow(finite_iterations)) {
      out <- stats::aggregate(
        seed ~ warmup_label + optimizer_iter_completed,
        data = finite_iterations,
        FUN = length
      )
      names(out)[names(out) == "seed"] <- "seed_frequency"
      out
    } else data.frame()
    write_tsv(iteration_counts, file.path(out_dir, "multi_warmup_deoptim_iteration_counts.tsv"))
  }
  write_tsv(invivo_only_long, file.path(out_dir, "multi_warmup_invivo_only_parameter_long.tsv"))
  if (!is.null(basin_df)) write_tsv(basin_df, file.path(out_dir, "multi_warmup_final_basin_assignments.tsv"))
  collect_initial_final_distance_diagnostics(manifest, root, summary_df, out_dir)
  message("Wrote multi-warmup summary: ", file.path(out_dir, "multi_warmup_best_seed_summary.tsv"))
}
