#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

invivo_only_oxygen_parameters <- c("o2_S0", "kappa_O", "eta_o2", "rho_2N")
invivo_only_oxygen_feature_cols <- paste0("log10_invivo_", invivo_only_oxygen_parameters)

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos > 0L) out[[substr(kv, 1L, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv)) else out[[kv]] <- TRUE
  }
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) y else x
as_chr <- function(x, default = "") as.character((x %||% default)[[1]])

truthy <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]])) return(default)
  val <- tolower(trimws(as.character(x[[1]])))
  if (!nzchar(val)) return(default)
  val %in% c("1", "true", "t", "yes", "y", "on")
}

read_table_optional <- function(path, sep = "\t") {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}

num_col <- function(x) suppressWarnings(as.numeric(x))

format_heatmap_value <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.finite(x),
    ifelse(abs(x) >= 1000 | abs(x) < 0.01, formatC(x, format = "e", digits = 2), format(signif(x, 3), trim = TRUE)),
    ""
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

family_palette <- function(levels) {
  levels <- as.character(levels)
  levels <- levels[nzchar(levels)]
  palette <- c(
    "#009E73",
    "#6A3D9A",
    "#E69F00",
    "#000000",
    "#8C6D31",
    "#F0E442",
    "#66A61E",
    "#7F7F7F",
    "#B07AA1",
    "#A6761D"
  )
  if (length(levels) > length(palette)) {
    palette <- grDevices::colorRampPalette(palette)(length(levels))
  } else {
    palette <- palette[seq_along(levels)]
  }
  stats::setNames(palette, levels)
}

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
      "huber_k", "huber_threshold", "abs_delta_over_sigma", "penalty_region"
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

plot_objectives <- function(summary_df, out_dir) {
  if (!nrow(summary_df)) return(invisible(NULL))
  family_levels <- unique(as.character(summary_df$invivo_family))
  family_levels <- family_levels[nzchar(family_levels)]
  summary_df$invivo_family <- factor(as.character(summary_df$invivo_family), levels = family_levels)
  p1 <- ggplot(summary_df, aes(x = reorder(warmup_label, objective), y = objective, color = invivo_family)) +
    geom_point(size = 2.8) +
    coord_flip() +
    scale_color_manual(values = family_palette(family_levels), drop = FALSE) +
    labs(title = "Multi-Warm-Up Best Joint Objective", x = "Warm-up pair", y = "Best total objective", color = "In vivo family") +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "multi_warmup_objective_summary.pdf"), p1, width = 8, height = 5, bg = "white")

  if (all(c("objective_invivo", "objective_invitro") %in% names(summary_df))) {
    p2 <- ggplot(summary_df, aes(objective_invivo, objective_invitro, color = invivo_family, label = warmup_label)) +
      geom_point(size = 2.8) +
      scale_color_manual(values = family_palette(family_levels), drop = FALSE) +
      labs(title = "In Vivo vs In Vitro Objective by Warm-Up Pair", x = "Best seed in vivo objective", y = "Best seed in vitro objective", color = "In vivo family") +
      theme_minimal(base_size = 11)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p2 <- p2 + ggrepel::geom_text_repel(size = 3, max.overlaps = Inf)
    } else {
      p2 <- p2 + geom_text(size = 3, vjust = -0.6)
    }
    ggsave(file.path(out_dir, "multi_warmup_invivo_invitro_objective_scatter.pdf"), p2, width = 7, height = 5, bg = "white")
  }
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

plot_deoptim_iteration_histograms <- function(iter_df, out_dir) {
  if (is.null(iter_df) || !is.data.frame(iter_df) || !nrow(iter_df)) return(invisible(NULL))
  iter_df$optimizer_iter_completed <- num_col(iter_df$optimizer_iter_completed)
  iter_df <- iter_df[is.finite(iter_df$optimizer_iter_completed), , drop = FALSE]
  if (!nrow(iter_df)) return(invisible(NULL))
  pair_levels <- unique(iter_df$warmup_label)
  iter_df$warmup_label <- factor(iter_df$warmup_label, levels = pair_levels)
  counts <- stats::aggregate(
    seed ~ warmup_label + optimizer_iter_completed,
    data = iter_df,
    FUN = length
  )
  names(counts)[names(counts) == "seed"] <- "seed_frequency"
  write_tsv(counts, file.path(out_dir, "multi_warmup_deoptim_iteration_counts.tsv"))
  n_pairs <- length(pair_levels)
  ncol_use <- min(3L, max(1L, n_pairs))
  height_use <- max(4.8, 2.4 * ceiling(n_pairs / ncol_use))
  width_use <- max(8, 3.8 * ncol_use)
  iter_breaks <- sort(unique(counts$optimizer_iter_completed))
  bar_width <- if (length(iter_breaks) <= 1L) 0.35 else 0.85
  x_scale <- if (length(iter_breaks) <= 12L) {
    scale_x_continuous(breaks = iter_breaks)
  } else {
    scale_x_continuous(n.breaks = 5)
  }
  p <- ggplot(counts, aes(x = optimizer_iter_completed, y = seed_frequency)) +
    geom_col(width = bar_width, fill = "#2b6cb0", alpha = 0.82) +
    facet_wrap(~warmup_label, scales = "free_y", ncol = ncol_use) +
    x_scale +
    labs(
      title = "DEoptim Iteration Frequency by Warm-Up Pair",
      x = "Completed DEoptim iterations",
      y = "Seed frequency"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
  ggsave(file.path(out_dir, "multi_warmup_deoptim_iteration_histograms.pdf"), p, width = width_use, height = height_use, bg = "white")
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

plot_invivo_only_warmstart_parameters <- function(invivo_only_long, out_dir) {
  if (is.null(invivo_only_long) || !is.data.frame(invivo_only_long) || !nrow(invivo_only_long)) return(invisible(NULL))
  plot_df <- invivo_only_long[is.finite(invivo_only_long$scaled_log10_value), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  n_seed <- length(unique(plot_df$invivo_label))
  height_use <- max(4.5, 0.42 * n_seed + 2.2)
  p <- ggplot(plot_df, aes(parameter, invivo_label, fill = scaled_log10_value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = display_value, color = label_color), size = 2.7) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#d33f3f", midpoint = 0, name = "scaled\nlog10 value") +
    scale_color_identity() +
    labs(
      title = "In Vivo-Only Warm-Start Parameters",
      subtitle = "Source in vivo seed values for parameters not represented in paired soft-coupling ratios.",
      x = "Parameter",
      y = "In vivo warm-start seed"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  ggsave(file.path(out_dir, "multi_warmup_invivo_only_parameter_heatmap.pdf"), p, width = 7.5, height = height_use, bg = "white")
}

plot_parameter_ratios <- function(ratio_long, out_dir) {
  if (is.null(ratio_long) || !nrow(ratio_long)) return(invisible(NULL))
  ratio_long <- add_final_log10_ratio(ratio_long)
  ratio_long <- ratio_long[is.finite(ratio_long$final_log10_ratio_vivo_to_vitro), , drop = FALSE]
  if (!nrow(ratio_long)) return(invisible(NULL))
  p <- ggplot(ratio_long, aes(parameter, warmup_label, fill = final_log10_ratio_vivo_to_vitro)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#d33f3f", midpoint = 0, name = "log10 ratio") +
    labs(title = "Best-Seed Final In Vivo / In Vitro Parameter Ratios", x = "Parameter", y = "Warm-up pair") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "multi_warmup_final_parameter_ratio_heatmap.pdf"), p, width = 9, height = 5.5, bg = "white")
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

seed_order_key <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(sub("^seed", "", x)))
  out[!is.finite(out)] <- Inf
  out
}

single_finite_value <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) x[[1]] else NA_real_
}

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

plot_initial_final_distance_diagnostics <- function(per_seed, heatmap_long, out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (is.data.frame(per_seed) && nrow(per_seed)) {
    pair_levels <- if ("warmup_order" %in% names(per_seed)) {
      unique(as.character(per_seed$warmup_label[order(per_seed$warmup_order)]))
    } else {
      unique(as.character(per_seed$warmup_label))
    }
    per_seed$warmup_label <- factor(as.character(per_seed$warmup_label), levels = pair_levels)
    per_seed$final_basin_id <- factor(as.character(per_seed$final_basin_id), levels = unique(as.character(per_seed$final_basin_id)))
    p <- ggplot(per_seed, aes(x = warmup_label, y = distance_to_own_warmup, color = final_basin_id)) +
      geom_jitter(width = 0.18, height = 0, size = 2.2, alpha = 0.82) +
      labs(
        title = "Initial-to-Final 18D Profile Distance by Joint Seed",
        subtitle = "Distance is computed in scaled 18D warm-start profile space.",
        x = "Warm-up pair",
        y = "Distance to own warm-up profile",
        color = "Final basin"
      ) +
      theme_minimal(base_size = 11) +
      guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3, alpha = 1))) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom"
      )
    square_size <- 7.2
    ggsave(file.path(out_dir, "multi_warmup_initial_to_final_distance.pdf"), p, width = square_size, height = square_size, bg = "white")
    ggsave(file.path(out_dir, "multi_warmup_initial_to_final_distance.png"), p, width = square_size, height = square_size, dpi = 220, bg = "white")
  }

  if (is.data.frame(heatmap_long) && nrow(heatmap_long)) {
    heatmap_long$final_representative <- factor(
      as.character(heatmap_long$final_representative),
      levels = rev(unique(as.character(heatmap_long$final_representative)))
    )
    heatmap_long$warmup_label <- factor(as.character(heatmap_long$warmup_label), levels = unique(as.character(heatmap_long$warmup_label)))
    heatmap_long$distance_label <- format_heatmap_value(heatmap_long$distance)
    p <- ggplot(heatmap_long, aes(x = warmup_label, y = final_representative, fill = distance)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = distance_label), size = 2.8) +
      scale_fill_gradient(low = "#f7fbff", high = "#4a1486", guide = "none") +
      labs(
        title = "Final Basin Representative vs Warm-Up 18D Profile Distance",
        subtitle = "Rows are best-objective representatives from each final basin;\ncolumns are initial warm-up profiles.",
        x = "Warm-up profile",
        y = "Final basin representative"
      ) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
    square_size <- 7.2
    ggsave(file.path(out_dir, "multi_warmup_final_to_warmup_distance_heatmap.pdf"), p, width = square_size, height = square_size, bg = "white")
    ggsave(file.path(out_dir, "multi_warmup_final_to_warmup_distance_heatmap.png"), p, width = square_size, height = square_size, dpi = 220, bg = "white")
  }
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
  plot_initial_final_distance_diagnostics(meta, heatmap_long, out_dir)
  invisible(list(per_seed = meta, basin_summary = basin_summary, heatmap = heatmap_long))
}

cleanup_removed_basin_pca_files <- function(out_dir) {
  unlink(
    file.path(
      out_dir,
      c("multi_warmup_final_parameter_basin_pca.pdf", "multi_warmup_final_parameter_pca_coords.tsv")
    ),
    force = TRUE
  )
}

safe_component <- function(x) {
  x <- trimws(as.character(x %||% ""))
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
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

remove_generated_integrated_run <- function(integrated_run_dir, root) {
  integrated_run_dir <- normalizePath(integrated_run_dir, mustWork = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  if (!startsWith(integrated_run_dir, paste0(root, .Platform$file.sep))) {
    stop("Refusing to remove integrated run outside multi-warmup root: ", integrated_run_dir)
  }
  if (!identical(basename(integrated_run_dir), "multi_warmup_integrated_joint_run")) {
    stop("Refusing to remove unexpected integrated run path: ", integrated_run_dir)
  }
  if (dir.exists(integrated_run_dir)) unlink(integrated_run_dir, recursive = TRUE, force = TRUE)
}

current_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }
  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) as.character(frame$ofile %||% ""), character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) return(dirname(normalizePath(ofiles[[length(ofiles)]], mustWork = TRUE)))
  normalizePath(".", mustWork = TRUE)
}

create_integrated_seed_links <- function(manifest, root, out_dir) {
  integrated_run_dir <- file.path(out_dir, "multi_warmup_integrated_joint_run")
  remove_generated_integrated_run(integrated_run_dir, out_dir)
  dir.create(integrated_run_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    warmup_label <- safe_component(if ("warmup_label" %in% names(manifest)) manifest$warmup_label[[i]] else i)
    joint_prefix <- if ("joint_run_prefix" %in% names(manifest)) as.character(manifest$joint_run_prefix[[i]]) else paste0("fit_joint_", warmup_label)
    source_run_dir <- normalizePath(file.path(root, joint_prefix), mustWork = FALSE)
    seed_dirs <- find_valid_seed_dirs(source_run_dir)
    if (!length(seed_dirs)) next
    seed_dirs <- seed_dirs[order(seed_dirs)]
    for (seed_dir in seed_dirs) {
      raw_seed <- basename(seed_dir)
      integrated_seed <- paste0(warmup_label, "__", raw_seed)
      target_dir <- file.path(integrated_run_dir, integrated_seed)
      if (!file.symlink(normalizePath(seed_dir, mustWork = TRUE), target_dir)) {
        stop("Failed to create integrated seed symlink: ", target_dir, " -> ", seed_dir)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        warmup_label = warmup_label,
        joint_run_prefix = joint_prefix,
        raw_seed = raw_seed,
        integrated_seed = integrated_seed,
        source_seed_dir = normalizePath(seed_dir, mustWork = TRUE),
        integrated_seed_dir = target_dir,
        stringsAsFactors = FALSE
      )
    }
  }

  manifest_out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write_tsv(manifest_out, file.path(out_dir, "multi_warmup_integrated_seed_manifest.tsv"))
  list(run_dir = integrated_run_dir, seed_manifest = manifest_out)
}

run_integrated_extra_results <- function(integrated_run_dir) {
  if (!dir.exists(integrated_run_dir)) return(invisible(NULL))
  seed_dirs <- find_valid_seed_dirs(integrated_run_dir)
  if (!length(seed_dirs)) return(invisible(NULL))
  script_dir <- current_script_dir()
  extra_results_script <- normalizePath(file.path(script_dir, "extra_results.R"), mustWork = FALSE)
  if (!file.exists(extra_results_script)) {
    extra_results_script <- normalizePath(file.path(getwd(), "oxygen/code/O2_supply_demand_MAP/analysis/extra_results.R"), mustWork = FALSE)
  }
  if (!file.exists(extra_results_script)) {
    stop("extra_results.R was not found for integrated multi-warmup report generation.")
  }
  out_dir <- file.path(integrated_run_dir, "extra_results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  args <- c(
    paste0("--run_dir=", integrated_run_dir),
    paste0("--out_dir=", out_dir),
    "--allow_partial_seed_dirs=TRUE"
  )
  log_path <- file.path(integrated_run_dir, "integrated_extra_results_run.log")
  rscript <- file.path(R.home("bin"), "Rscript")
  message("Running integrated joint extra_results.R for ", length(seed_dirs), " prefixed seeds.")
  out <- suppressWarnings(system2(rscript, args = c(extra_results_script, args), stdout = TRUE, stderr = TRUE))
  writeLines(out, log_path, useBytes = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("Integrated extra_results.R failed with status ", status, ". See: ", log_path)
  }
  invisible(out_dir)
}

build_integrated_joint_extra_results <- function(manifest, root, out_dir, enabled = TRUE) {
  if (!isTRUE(enabled)) return(invisible(NULL))
  if (is.null(manifest) || !is.data.frame(manifest) || !nrow(manifest)) return(invisible(NULL))
  if (!("joint_run_prefix" %in% names(manifest))) return(invisible(NULL))
  integrated <- create_integrated_seed_links(manifest = manifest, root = root, out_dir = out_dir)
  if (!is.data.frame(integrated$seed_manifest) || !nrow(integrated$seed_manifest)) {
    message("No completed seed directories were available for integrated joint extra_results.")
    return(invisible(NULL))
  }
  run_integrated_extra_results(integrated$run_dir)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root <- normalizePath(as_chr(argv$multi_warmup_root, as_chr(argv$out_dir)), mustWork = TRUE)
  manifest_path <- normalizePath(as_chr(argv$manifest, file.path(root, "multi_warmup_manifest.tsv")), mustWork = TRUE)
  out_dir <- normalizePath(as_chr(argv$out_dir, root), mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  build_integrated_joint_results <- truthy(argv$build_integrated_joint_results, default = TRUE)

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
  write_tsv(invivo_only_long, file.path(out_dir, "multi_warmup_invivo_only_parameter_long.tsv"))
  if (!is.null(basin_df)) write_tsv(basin_df, file.path(out_dir, "multi_warmup_final_basin_assignments.tsv"))
  cleanup_removed_basin_pca_files(out_dir)
  plot_objectives(summary_df, out_dir)
  plot_deoptim_iteration_histograms(iteration_long, out_dir)
  plot_parameter_ratios(ratio_long, out_dir)
  plot_invivo_only_warmstart_parameters(invivo_only_long, out_dir)
  collect_initial_final_distance_diagnostics(manifest, root, summary_df, out_dir)
  build_integrated_joint_extra_results(
    manifest = manifest,
    root = root,
    out_dir = out_dir,
    enabled = build_integrated_joint_results
  )
  message("Wrote multi-warmup summary: ", file.path(out_dir, "multi_warmup_best_seed_summary.tsv"))
}

if (identical(environment(), globalenv())) {
  main()
}
