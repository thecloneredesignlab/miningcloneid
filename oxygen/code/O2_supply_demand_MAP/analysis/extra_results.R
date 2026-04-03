#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

.o2sd_bootstrap_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2sd_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_bool <- o2sd_as_bool

seed_order_key <- function(x) {
  m <- regexec("^seed([0-9]+)$", x)
  hit <- regmatches(x, m)
  out <- rep(Inf, length(x))
  for (i in seq_along(hit)) {
    if (length(hit[[i]]) == 2L) out[[i]] <- as.numeric(hit[[i]][[2]])
  }
  out
}

find_seed_dirs <- function(run_dir) {
  dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  keep <- vapply(
    dirs,
    function(d) {
      file.exists(file.path(d, "fit_summary.tsv")) &&
        file.exists(file.path(d, "best_params.tsv"))
    },
    logical(1)
  )
  dirs <- dirs[keep]
  if (length(dirs) == 0L &&
      file.exists(file.path(run_dir, "fit_summary.tsv")) &&
      file.exists(file.path(run_dir, "best_params.tsv"))) {
    dirs <- run_dir
  }
  dirs[order(seed_order_key(basename(dirs)), basename(dirs))]
}

read_metric_map <- function(path, key_col, value_col) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c(key_col, value_col) %in% names(tab))) {
    stop("Missing required columns in ", path, ": ", key_col, ", ", value_col)
  }
  setNames(tab[[value_col]], tab[[key_col]])
}

lookup_named_num <- function(x, key) {
  if (!length(x) || is.null(key) || !nzchar(key) || !(key %in% names(x))) return(NA_real_)
  val <- suppressWarnings(as.numeric(x[[key]]))
  if (!is.finite(val)) NA_real_ else val
}

derive_delta_lam <- function(best_vals) {
  lam_min <- lookup_named_num(best_vals, "lam_min")
  lam_max <- lookup_named_num(best_vals, "lam_max")
  gap <- lam_max - lam_min
  if (!is.finite(gap) || gap <= 0) return(NA_real_)
  log(gap)
}

derive_prototype_value <- function(param_name, param_prototype, best_vals) {
  val <- lookup_named_num(best_vals, param_prototype)
  if (is.finite(val)) return(val)
  if (identical(param_name, "delta_lam") || identical(param_prototype, "delta_lam")) {
    return(derive_delta_lam(best_vals))
  }
  NA_real_
}

derive_transformed_value <- function(param_name, param_prototype, best_vals) {
  if (identical(param_name, "delta_lam")) {
    return(derive_delta_lam(best_vals))
  }
  proto_val <- derive_prototype_value(param_name, param_prototype, best_vals)
  if (!is.finite(proto_val)) return(NA_real_)
  if (startsWith(param_name, "log10_")) {
    if (proto_val <= 0) return(NA_real_)
    return(log10(proto_val))
  }
  proto_val
}

is_active_parameter <- function(param_prototype, estimate_flag, fit_summary_vals) {
  active <- isTRUE(estimate_flag)
  if (!active) return(FALSE)
  if (identical(param_prototype, "tau_O2")) {
    return(as_bool(fit_summary_vals[["fit_tau_O2"]], FALSE))
  }
  if (identical(param_prototype, "alpha") || identical(param_prototype, "gamma")) {
    return(as_bool(fit_summary_vals[["fit_treatment"]], FALSE))
  }
  if (identical(param_prototype, "alpha_o2") || identical(param_prototype, "gamma_growth")) {
    return(as_bool(fit_summary_vals[["O2_growth"]], TRUE))
  }
  TRUE
}

compute_boundary_metrics <- function(x, lower, upper, near_thresh = 0.05, tol = 1e-10) {
  if (!is.finite(x) || !is.finite(lower) || !is.finite(upper) || upper <= lower) {
    return(list(
      rel_pos_in_range = NA_real_,
      rel_pos_plot = NA_real_,
      rel_dist_to_lower = NA_real_,
      rel_dist_to_upper = NA_real_,
      rel_dist_to_nearest = NA_real_,
      nearest_bound = NA_character_,
      at_lower = FALSE,
      at_upper = FALSE,
      near_lower = FALSE,
      near_upper = FALSE,
      bound_status = "unavailable"
    ))
  }
  width <- upper - lower
  dist_lower <- x - lower
  dist_upper <- upper - x
  rel_pos <- dist_lower / width
  rel_lower <- dist_lower / width
  rel_upper <- dist_upper / width
  at_lower <- abs(dist_lower) <= max(tol, tol * width)
  at_upper <- abs(dist_upper) <= max(tol, tol * width)
  near_lower <- !at_lower && rel_lower <= near_thresh
  near_upper <- !at_upper && rel_upper <= near_thresh
  nearest <- if (dist_lower <= dist_upper) "lower" else "upper"
  status <- if (at_lower) {
    "at_lower"
  } else if (at_upper) {
    "at_upper"
  } else if (near_lower) {
    "near_lower"
  } else if (near_upper) {
    "near_upper"
  } else {
    "interior"
  }
  list(
    rel_pos_in_range = rel_pos,
    rel_pos_plot = min(max(rel_pos, 0), 1),
    rel_dist_to_lower = rel_lower,
    rel_dist_to_upper = rel_upper,
    rel_dist_to_nearest = min(rel_lower, rel_upper),
    nearest_bound = nearest,
    at_lower = at_lower,
    at_upper = at_upper,
    near_lower = at_lower || near_lower,
    near_upper = at_upper || near_upper,
    bound_status = status
  )
}

build_parameter_long_table <- function(seed, objective, fit_summary_vals, best_vals, param_table, near_thresh) {
  rows <- vector("list", nrow(param_table))
  for (i in seq_len(nrow(param_table))) {
    param_name <- trimws(as.character(param_table$param_name[[i]]))
    param_prototype <- trimws(as.character(param_table$param_prototype[[i]]))
    estimate_flag <- isTRUE(param_table$estimate[[i]])
    prototype_value <- derive_prototype_value(param_name, param_prototype, best_vals)
    transformed_value <- derive_transformed_value(param_name, param_prototype, best_vals)
    lower_bound <- suppressWarnings(as.numeric(param_table$lower_bound[[i]]))
    upper_bound <- suppressWarnings(as.numeric(param_table$upper_bound[[i]]))
    prototype_lower_bound <- suppressWarnings(as.numeric(param_table$prototype_lower_bound[[i]]))
    prototype_upper_bound <- suppressWarnings(as.numeric(param_table$prototype_upper_bound[[i]]))
    active_in_fit <- is_active_parameter(param_prototype, estimate_flag, fit_summary_vals)
    met <- compute_boundary_metrics(transformed_value, lower_bound, upper_bound, near_thresh = near_thresh)
    rows[[i]] <- data.frame(
      seed = seed,
      objective = objective,
      param_name = param_name,
      param_prototype = param_prototype,
      estimate = estimate_flag,
      active_in_fit = active_in_fit,
      prototype_value = prototype_value,
      transformed_value = transformed_value,
      lower_bound_transformed = lower_bound,
      upper_bound_transformed = upper_bound,
      prototype_lower_bound = prototype_lower_bound,
      prototype_upper_bound = prototype_upper_bound,
      rel_pos_in_range = met$rel_pos_in_range,
      rel_pos_plot = met$rel_pos_plot,
      rel_dist_to_lower = met$rel_dist_to_lower,
      rel_dist_to_upper = met$rel_dist_to_upper,
      rel_dist_to_nearest = met$rel_dist_to_nearest,
      nearest_bound = met$nearest_bound,
      at_lower = met$at_lower,
      at_upper = met$at_upper,
      near_lower = met$near_lower,
      near_upper = met$near_upper,
      bound_status = met$bound_status,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

build_seed_summary_record <- function(seed, fit_summary_vals, best_vals, parameter_long) {
  active_rows <- parameter_long[parameter_long$active_in_fit & is.finite(parameter_long$rel_dist_to_nearest), , drop = FALSE]
  active_no_sigma <- active_rows[active_rows$param_prototype != "sigma_burden", , drop = FALSE]
  n_at_bound_active <- sum(active_rows$at_lower | active_rows$at_upper, na.rm = TRUE)
  n_near_bound_only_active <- sum(active_rows$bound_status %in% c("near_lower", "near_upper"), na.rm = TRUE)
  n_at_bound_active_excl_sigma <- sum(active_no_sigma$at_lower | active_no_sigma$at_upper, na.rm = TRUE)
  n_near_bound_only_active_excl_sigma <- sum(active_no_sigma$bound_status %in% c("near_lower", "near_upper"), na.rm = TRUE)
  record <- list(
    seed = seed,
    objective = as_num(fit_summary_vals[["objective"]]),
    objective_data = as_num(fit_summary_vals[["objective_data"]]),
    objective_prior = as_num(fit_summary_vals[["objective_prior"]]),
    objective_burden = as_num(fit_summary_vals[["objective_burden"]]),
    objective_ploidy = as_num(fit_summary_vals[["objective_ploidy"]]),
    objective_burden_neg2loglik_raw = as_num(fit_summary_vals[["objective_burden_neg2loglik_raw"]]),
    objective_ploidy_neg2loglik_raw = as_num(fit_summary_vals[["objective_ploidy_neg2loglik_raw"]]),
    optimizer_interrupted = as.character(fit_summary_vals[["optimizer_interrupted"]] %||% NA_character_),
    optimizer_iter_completed = as_num(fit_summary_vals[["optimizer_iter_completed"]]),
    n_active_params = nrow(active_rows),
    n_at_bound_active = n_at_bound_active,
    n_at_lower_active = sum(active_rows$at_lower, na.rm = TRUE),
    n_at_upper_active = sum(active_rows$at_upper, na.rm = TRUE),
    n_near_bound_only_active = n_near_bound_only_active,
    n_near_lower_active = sum(active_rows$near_lower, na.rm = TRUE),
    n_near_upper_active = sum(active_rows$near_upper, na.rm = TRUE),
    boundary_penalty_active = 2 * n_at_bound_active + n_near_bound_only_active,
    min_rel_dist_active = if (nrow(active_rows)) min(active_rows$rel_dist_to_nearest, na.rm = TRUE) else NA_real_,
    median_rel_dist_active = if (nrow(active_rows)) median(active_rows$rel_dist_to_nearest, na.rm = TRUE) else NA_real_,
    worst_param_active = if (nrow(active_rows)) active_rows$param_prototype[[which.min(active_rows$rel_dist_to_nearest)]] else NA_character_,
    n_active_params_excl_sigma_burden = nrow(active_no_sigma),
    n_at_bound_active_excl_sigma_burden = n_at_bound_active_excl_sigma,
    n_near_bound_only_active_excl_sigma_burden = n_near_bound_only_active_excl_sigma,
    boundary_penalty_active_excl_sigma_burden = 2 * n_at_bound_active_excl_sigma + n_near_bound_only_active_excl_sigma,
    min_rel_dist_active_excl_sigma_burden = if (nrow(active_no_sigma)) min(active_no_sigma$rel_dist_to_nearest, na.rm = TRUE) else NA_real_,
    median_rel_dist_active_excl_sigma_burden = if (nrow(active_no_sigma)) median(active_no_sigma$rel_dist_to_nearest, na.rm = TRUE) else NA_real_,
    worst_param_active_excl_sigma_burden = if (nrow(active_no_sigma)) active_no_sigma$param_prototype[[which.min(active_no_sigma$rel_dist_to_nearest)]] else NA_character_
  )

  bp_names <- sort(names(best_vals))
  for (nm in bp_names) {
    record[[paste0("value__", nm)]] <- suppressWarnings(as.numeric(best_vals[[nm]]))
  }

  for (i in seq_len(nrow(parameter_long))) {
    proto <- parameter_long$param_prototype[[i]]
    record[[paste0("active__", proto)]] <- parameter_long$active_in_fit[[i]]
    record[[paste0("value__", proto, "__prototype")]] <- parameter_long$prototype_value[[i]]
    record[[paste0("value__", proto, "__transformed")]] <- parameter_long$transformed_value[[i]]
    record[[paste0("status__", proto)]] <- parameter_long$bound_status[[i]]
    record[[paste0("nearest_bound__", proto)]] <- parameter_long$nearest_bound[[i]]
    record[[paste0("rel_pos__", proto)]] <- parameter_long$rel_pos_in_range[[i]]
    record[[paste0("rel_dist__", proto)]] <- parameter_long$rel_dist_to_nearest[[i]]
    record[[paste0("near_lower__", proto)]] <- parameter_long$near_lower[[i]]
    record[[paste0("near_upper__", proto)]] <- parameter_long$near_upper[[i]]
  }

  record
}

bind_records <- function(records) {
  all_names <- unique(unlist(lapply(records, names), use.names = FALSE))
  out <- vector("list", length(records))
  for (i in seq_along(records)) {
    x <- records[[i]]
    row <- as.list(rep(NA, length(all_names)))
    names(row) <- all_names
    row[names(x)] <- x
    out[[i]] <- as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  }
  do.call(rbind, out)
}

plot_parameter_boundary_forest <- function(long_df, summary_df, out_path, run_label, near_thresh = 0.05) {
  plot_df <- long_df[long_df$active_in_fit & is.finite(long_df$rel_pos_plot), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))

  rank_col <- if ("recommend_rank_burden_ploidy_boundary_first" %in% names(summary_df)) {
    "recommend_rank_burden_ploidy_boundary_first"
  } else if ("recommend_rank_ploidy_burden_boundary_first" %in% names(summary_df)) {
    "recommend_rank_ploidy_burden_boundary_first"
  } else if ("recommend_rank_ploidy_boundary_first" %in% names(summary_df)) {
    "recommend_rank_ploidy_boundary_first"
  } else if ("recommend_rank_ploidy_first" %in% names(summary_df)) {
    "recommend_rank_ploidy_first"
  } else if ("objective_burden_rank" %in% names(summary_df)) {
    "objective_burden_rank"
  } else if ("objective_ploidy_rank" %in% names(summary_df)) {
    "objective_ploidy_rank"
  } else {
    "objective"
  }
  rank_ord <- order(summary_df[[rank_col]], summary_df$objective, summary_df$seed, na.last = TRUE)
  ranked_seeds <- summary_df$seed[rank_ord]
  top3_seeds <- head(ranked_seeds, 3L)

  param_rank <- tapply(plot_df$rel_dist_to_nearest, plot_df$param_prototype, min, na.rm = TRUE)
  param_levels <- names(sort(param_rank, decreasing = FALSE))
  plot_df$param_prototype <- factor(plot_df$param_prototype, levels = rev(param_levels))
  ref_df <- unique(plot_df["param_prototype"])
  plot_df$seed_marker <- "Other seeds"
  if (length(top3_seeds) >= 1L) plot_df$seed_marker[plot_df$seed == top3_seeds[[1]]] <- paste0("Top 1: ", top3_seeds[[1]], " (*)")
  if (length(top3_seeds) >= 2L) plot_df$seed_marker[plot_df$seed == top3_seeds[[2]]] <- paste0("Top 2: ", top3_seeds[[2]], " (triangle)")
  if (length(top3_seeds) >= 3L) plot_df$seed_marker[plot_df$seed == top3_seeds[[3]]] <- paste0("Top 3: ", top3_seeds[[3]], " (black dot)")
  other_df <- plot_df[plot_df$seed_marker == "Other seeds", , drop = FALSE]
  top_df <- plot_df[plot_df$seed_marker != "Other seeds", , drop = FALSE]
  top_breaks <- c(
    if (length(top3_seeds) >= 1L) paste0("Top 1: ", top3_seeds[[1]], " (*)"),
    if (length(top3_seeds) >= 2L) paste0("Top 2: ", top3_seeds[[2]], " (triangle)"),
    if (length(top3_seeds) >= 3L) paste0("Top 3: ", top3_seeds[[3]], " (black dot)")
  )
  shape_values <- c(
    if (length(top3_seeds) >= 1L) setNames(8, paste0("Top 1: ", top3_seeds[[1]], " (*)")),
    if (length(top3_seeds) >= 2L) setNames(17, paste0("Top 2: ", top3_seeds[[2]], " (triangle)")),
    if (length(top3_seeds) >= 3L) setNames(16, paste0("Top 3: ", top3_seeds[[3]], " (black dot)"))
  )
  point_pos <- position_jitter(height = 0.14, width = 0)

  p <- ggplot(plot_df, aes(x = rel_pos_plot, y = param_prototype)) +
    annotate("rect", xmin = 0, xmax = near_thresh, ymin = -Inf, ymax = Inf, fill = "#fddbc7", alpha = 0.28) +
    annotate("rect", xmin = 1 - near_thresh, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#d1e5f0", alpha = 0.28) +
    geom_segment(
      data = ref_df,
      aes(x = 0, xend = 1, y = param_prototype, yend = param_prototype),
      inherit.aes = FALSE,
      color = "grey78",
      linewidth = 0.5
    ) +
    geom_vline(xintercept = c(0, near_thresh, 0.5, 1 - near_thresh, 1), color = c("grey50", "grey70", "grey82", "grey70", "grey50"), linetype = c("solid", "dashed", "dotted", "dashed", "solid"), linewidth = 0.35) +
    geom_point(
      data = other_df,
      shape = 16,
      size = 2.1,
      alpha = 0.5,
      color = "grey65",
      position = point_pos
    ) +
    geom_point(
      data = top_df,
      aes(shape = seed_marker),
      size = 3.0,
      color = "black",
      position = point_pos
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, near_thresh, 0.5, 1 - near_thresh, 1)) +
    scale_shape_manual(values = shape_values, breaks = top_breaks, drop = FALSE) +
    labs(
      title = paste0("Parameter Positions Within Fitted Bounds: ", run_label),
      subtitle = paste0("0 = lower bound, 1 = upper bound; shaded zones are within ", sprintf("%.0f", 100 * near_thresh), "% of a bound"),
      x = "Relative position in transformed fit range",
      y = NULL,
      shape = "Recommended Top 3 Seeds"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  ggplot2::ggsave(out_path, p, width = 13, height = max(8, 0.42 * length(param_levels) + 3))
  invisible(out_path)
}

plot_objective_vs_boundary_risk <- function(summary_df, out_path, run_label) {
  plot_df <- summary_df[
    is.finite(summary_df$objective) &
      is.finite(summary_df$min_rel_dist_active_excl_sigma_burden) &
      summary_df$min_rel_dist_active_excl_sigma_burden > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))
  plot_df <- plot_df[order(plot_df$objective), , drop = FALSE]
  plot_df$seed <- factor(plot_df$seed, levels = plot_df$seed)

  p <- ggplot(plot_df, aes(x = objective, y = min_rel_dist_active_excl_sigma_burden, color = seed, label = seed)) +
    geom_point(size = 2.8) +
    geom_text(nudge_y = 0.015, show.legend = FALSE, size = 3) +
    scale_y_log10() +
    labs(
      title = paste0("Objective vs Boundary Risk: ", run_label),
      subtitle = "Boundary risk shown as minimum relative distance to a fitted bound, excluding sigma_burden",
      x = "Objective",
      y = "Min relative distance to nearest bound (log10 scale)",
      color = "Seed"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggplot2::ggsave(out_path, p, width = 10, height = 7)
  invisible(out_path)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  run_dir <- args$run_dir
  if (is.null(run_dir) || !nzchar(trimws(run_dir))) {
    stop("Usage: Rscript extra_results.R --run_dir=/path/to/run [--out_dir=/path/to/out] [--near_thresh=0.05]")
  }
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  out_dir <- normalizePath(args$out_dir %||% file.path(run_dir, "extra_results"), mustWork = FALSE)
  near_thresh <- as_num(args$near_thresh, 0.05)
  if (!is.finite(near_thresh) || near_thresh <= 0 || near_thresh >= 0.5) {
    stop("near_thresh must be in (0, 0.5).")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  parameter_table_path <- file.path(run_dir, "parameter_table.csv")
  if (!file.exists(parameter_table_path)) {
    stop("Missing parameter_table.csv in run directory: ", run_dir)
  }
  param_table <- utils::read.csv(parameter_table_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("param_name", "estimate", "lower_bound", "upper_bound", "param_prototype", "prototype_lower_bound", "prototype_upper_bound")
  if (!all(required_cols %in% names(param_table))) {
    stop("parameter_table.csv is missing required columns: ", paste(setdiff(required_cols, names(param_table)), collapse = ", "))
  }
  param_table$estimate <- vapply(param_table$estimate, function(x) as_bool(x, FALSE), logical(1))

  seed_dirs <- find_seed_dirs(run_dir)
  if (!length(seed_dirs)) {
    stop("No valid seed directories found under: ", run_dir)
  }

  long_rows <- vector("list", length(seed_dirs))
  summary_records <- vector("list", length(seed_dirs))

  for (i in seq_along(seed_dirs)) {
    seed_dir <- seed_dirs[[i]]
    seed <- basename(seed_dir)
    fit_summary_vals <- read_metric_map(file.path(seed_dir, "fit_summary.tsv"), "metric", "value")
    best_vals <- read_metric_map(file.path(seed_dir, "best_params.tsv"), "parameter", "value")
    objective <- as_num(fit_summary_vals[["objective"]])
    long_df <- build_parameter_long_table(
      seed = seed,
      objective = objective,
      fit_summary_vals = fit_summary_vals,
      best_vals = best_vals,
      param_table = param_table,
      near_thresh = near_thresh
    )
    long_rows[[i]] <- long_df
    summary_records[[i]] <- build_seed_summary_record(
      seed = seed,
      fit_summary_vals = fit_summary_vals,
      best_vals = best_vals,
      parameter_long = long_df
    )
  }

  parameter_long <- do.call(rbind, long_rows)
  seed_summary <- bind_records(summary_records)
  seed_summary$objective <- suppressWarnings(as.numeric(seed_summary$objective))
  seed_summary$objective_ploidy <- suppressWarnings(as.numeric(seed_summary$objective_ploidy))
  seed_summary$objective_burden <- suppressWarnings(as.numeric(seed_summary$objective_burden))
  seed_summary$objective_ploidy_neg2loglik_raw <- suppressWarnings(as.numeric(seed_summary$objective_ploidy_neg2loglik_raw))
  seed_summary$objective_burden_neg2loglik_raw <- suppressWarnings(as.numeric(seed_summary$objective_burden_neg2loglik_raw))
  seed_summary$boundary_penalty_active <- suppressWarnings(as.numeric(seed_summary$boundary_penalty_active))
  seed_summary$min_rel_dist_active <- suppressWarnings(as.numeric(seed_summary$min_rel_dist_active))
  seed_summary$boundary_penalty_active_excl_sigma_burden <- suppressWarnings(as.numeric(seed_summary$boundary_penalty_active_excl_sigma_burden))
  seed_summary$min_rel_dist_active_excl_sigma_burden <- suppressWarnings(as.numeric(seed_summary$min_rel_dist_active_excl_sigma_burden))
  seed_summary$objective_rank <- rank(seed_summary$objective, ties.method = "first", na.last = "keep")
  seed_summary$objective_ploidy_rank <- rank(seed_summary$objective_ploidy, ties.method = "first", na.last = "keep")
  seed_summary$objective_burden_rank <- rank(seed_summary$objective_burden, ties.method = "first", na.last = "keep")
  boundary_order <- order(
    seed_summary$boundary_penalty_active,
    -seed_summary$min_rel_dist_active,
    seed_summary$objective_burden,
    seed_summary$objective,
    seed_summary$seed,
    na.last = TRUE
  )
  seed_summary$boundary_rank_active_support <- NA_integer_
  seed_summary$boundary_rank_active_support[boundary_order] <- seq_len(nrow(seed_summary))
  seed_summary$recommend_score_burden_ploidy_boundary <- with(
    seed_summary,
    objective_burden_rank + 0.2 * objective_ploidy_rank + 0.1 * boundary_rank_active_support
  )
  recommend_order <- order(
    seed_summary$objective_burden,
    seed_summary$objective_ploidy,
    seed_summary$boundary_penalty_active,
    -seed_summary$min_rel_dist_active,
    seed_summary$objective,
    seed_summary$seed,
    na.last = TRUE
  )
  seed_summary$recommend_score_ploidy_burden_boundary <- seed_summary$recommend_score_burden_ploidy_boundary
  seed_summary$recommend_score_ploidy_boundary <- seed_summary$recommend_score_burden_ploidy_boundary
  seed_summary$recommend_rank_burden_ploidy_boundary_first <- NA_integer_
  seed_summary$recommend_rank_burden_ploidy_boundary_first[recommend_order] <- seq_len(nrow(seed_summary))
  seed_summary$recommend_rank_ploidy_burden_boundary_first <- NA_integer_
  seed_summary$recommend_rank_ploidy_burden_boundary_first[recommend_order] <- seq_len(nrow(seed_summary))
  seed_summary$recommend_rank_ploidy_boundary_first <- NA_integer_
  seed_summary$recommend_rank_ploidy_boundary_first[recommend_order] <- seq_len(nrow(seed_summary))
  seed_summary$recommend_rank_ploidy_first <- seed_summary$recommend_rank_burden_ploidy_boundary_first
  seed_summary <- seed_summary[order(seed_summary$objective, seed_summary$seed), , drop = FALSE]
  row.names(seed_summary) <- NULL

  utils::write.table(
    seed_summary,
    file = file.path(out_dir, "seed_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    parameter_long,
    file = file.path(out_dir, "parameter_boundary_long.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  plot_parameter_boundary_forest(
    long_df = parameter_long,
    summary_df = seed_summary,
    out_path = file.path(out_dir, "parameter_boundary_forest.pdf"),
    run_label = basename(run_dir),
    near_thresh = near_thresh
  )
  plot_objective_vs_boundary_risk(
    summary_df = seed_summary,
    out_path = file.path(out_dir, "objective_vs_boundary_risk.pdf"),
    run_label = basename(run_dir)
  )

  message("Wrote summary table: ", file.path(out_dir, "seed_summary.tsv"))
  message("Wrote parameter long table: ", file.path(out_dir, "parameter_boundary_long.tsv"))
  message("Wrote forest plot: ", file.path(out_dir, "parameter_boundary_forest.pdf"))
  message("Wrote objective-risk plot: ", file.path(out_dir, "objective_vs_boundary_risk.pdf"))
}

if (sys.nframe() == 0) {
  main()
}
