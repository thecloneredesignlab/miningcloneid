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
      "vivo_transformed", "vitro_transformed"
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
  p1 <- ggplot(summary_df, aes(x = reorder(warmup_label, objective), y = objective, color = invivo_family)) +
    geom_point(size = 2.8) +
    coord_flip() +
    labs(title = "Multi-Warm-Up Best Joint Objective", x = "Warm-up pair", y = "Best total objective", color = "In vivo family") +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "multi_warmup_objective_summary.pdf"), p1, width = 8, height = 5, bg = "white")

  if (all(c("objective_invivo", "objective_invitro") %in% names(summary_df))) {
    p2 <- ggplot(summary_df, aes(objective_invivo, objective_invitro, color = invivo_family, label = warmup_label)) +
      geom_point(size = 2.8) +
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
