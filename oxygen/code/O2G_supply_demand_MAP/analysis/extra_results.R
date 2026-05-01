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
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_bool <- o2sd_as_bool

summary_flag_true <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  tolower(trimws(as.character(x[[1]]))) %in% c("true", "t", "1", "yes", "y", "on")
}

summary_metric_value <- function(fit_summary_vals, key, default = NULL) {
  if (is.null(fit_summary_vals) || !length(fit_summary_vals) || is.null(key) || !nzchar(key)) {
    return(default)
  }
  if (!(key %in% names(fit_summary_vals))) {
    return(default)
  }
  fit_summary_vals[[key]]
}

filter_best_vals_for_output <- function(best_vals, fit_summary_vals) {
  glucose_use <- canonical_glucose_enabled(summary_metric_value(fit_summary_vals, "glucose", TRUE), TRUE)
  glucose_dynamic_use <- canonical_glucose_dynamic(summary_metric_value(fit_summary_vals, "glucose_dynamic", FALSE), FALSE)
  loss_mode <- canonical_misseg_loss_survival_mode(
    summary_metric_value(fit_summary_vals, "misseg_loss_survival", "nullisomy"),
    "nullisomy"
  )
  harvest_use <- summary_flag_true(summary_metric_value(fit_summary_vals, "harvest_init_multiplier", FALSE), default = FALSE)
  out <- filter_family_specific_run_params_for_output_common(
    run_params = best_vals,
    glucose = glucose_use,
    glucose_dynamic = glucose_dynamic_use,
    misseg_loss_survival = loss_mode
  )
  if (!isTRUE(harvest_use)) {
    out <- out[setdiff(names(out), grep("^init_mult_", names(out), value = TRUE))]
  }
  out
}

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
  if (startsWith(param_name, "logit_")) {
    if (proto_val <= 0 || proto_val >= 1) return(NA_real_)
    return(qlogis(proto_val))
  }
  proto_val
}

is_active_parameter <- function(param_name, param_prototype, estimate_flag, fit_summary_vals) {
  active <- isTRUE(estimate_flag)
  if (!active) return(FALSE)
  glucose_use <- canonical_glucose_enabled(summary_metric_value(fit_summary_vals, "glucose", TRUE), TRUE)
  glucose_dynamic_use <- canonical_glucose_dynamic(summary_metric_value(fit_summary_vals, "glucose_dynamic", FALSE), FALSE)
  loss_mode <- canonical_misseg_loss_survival_mode(
    summary_metric_value(fit_summary_vals, "misseg_loss_survival", "nullisomy"),
    "nullisomy"
  )
  harvest_use <- summary_flag_true(summary_metric_value(fit_summary_vals, "harvest_init_multiplier", FALSE), FALSE)
  if (identical(param_prototype, "tau_O2")) {
    return(as_bool(summary_metric_value(fit_summary_vals, "fit_tau_O2", FALSE), FALSE))
  }
  if (identical(param_prototype, "alpha") || identical(param_prototype, "gamma")) {
    return(as_bool(summary_metric_value(fit_summary_vals, "fit_treatment", FALSE), FALSE))
  }
  if (identical(param_prototype, "alpha_o2") || identical(param_prototype, "gamma_growth")) {
    return(as_bool(summary_metric_value(fit_summary_vals, "O2_growth", TRUE), TRUE))
  }
  if (identical(param_prototype, "k_o")) {
    return(!isTRUE(glucose_use) || !isTRUE(glucose_dynamic_use))
  }
  if (identical(param_prototype, "p_wgd")) {
    return(!isTRUE(glucose_use))
  }
  if (identical(param_prototype, "p_wgd_max") || identical(param_prototype, "O2_wgd")) {
    return(isTRUE(glucose_use))
  }
  if (param_prototype %in% c("G_S0", "kappa_G", "eta_G", "G_c", "tau_G")) {
    return(isTRUE(glucose_use) && isTRUE(glucose_dynamic_use))
  }
  if (identical(param_prototype, "gamma_loss")) {
    return(identical(loss_mode, "nullisomy"))
  }
  if (param_prototype %in% c("buffer_smax", "buffer_beta", "buffer_n_exp")) {
    return(identical(loss_mode, "buffering"))
  }
  if (identical(param_prototype, "harvest_init_multiplier") || startsWith(param_name, "log_init_mult_")) {
    return(isTRUE(harvest_use))
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
    active_in_fit <- is_active_parameter(param_name, param_prototype, estimate_flag, fit_summary_vals)
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

build_seed_summary_record <- function(seed, fit_summary_vals, best_vals, parameter_long, pred_gate_metrics = NULL) {
  best_vals <- filter_best_vals_for_output(best_vals, fit_summary_vals)
  active_rows <- parameter_long[parameter_long$active_in_fit & is.finite(parameter_long$rel_dist_to_nearest), , drop = FALSE]
  active_no_sigma <- active_rows[active_rows$param_prototype != "sigma_burden", , drop = FALSE]
  n_at_bound_active <- sum(active_rows$at_lower | active_rows$at_upper, na.rm = TRUE)
  n_near_bound_only_active <- sum(active_rows$bound_status %in% c("near_lower", "near_upper"), na.rm = TRUE)
  n_at_bound_active_excl_sigma <- sum(active_no_sigma$at_lower | active_no_sigma$at_upper, na.rm = TRUE)
  n_near_bound_only_active_excl_sigma <- sum(active_no_sigma$bound_status %in% c("near_lower", "near_upper"), na.rm = TRUE)
  pred_gate_metrics <- pred_gate_metrics %||% list()
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
    worst_param_active_excl_sigma_burden = if (nrow(active_no_sigma)) active_no_sigma$param_prototype[[which.min(active_no_sigma$rel_dist_to_nearest)]] else NA_character_,
    pred1000_2N = suppressWarnings(as.numeric(pred_gate_metrics$pred1000_2N %||% NA_real_)),
    pred1000_4N = suppressWarnings(as.numeric(pred_gate_metrics$pred1000_4N %||% NA_real_)),
    pred1000_both_gt44 = isTRUE(pred_gate_metrics$pred1000_both_gt44 %||% FALSE)
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

read_1000day_ploidy_gate_metrics <- function(seed_dir, target_day = 1000, threshold = 44) {
  path <- file.path(seed_dir, "viz", "predict_ploidy_weighted_mean_0_1000day.tsv")
  out <- list(
    pred1000_2N = NA_real_,
    pred1000_4N = NA_real_,
    pred1000_both_gt44 = FALSE
  )
  if (!file.exists(path)) return(out)

  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || !nrow(tab) || !all(c("cohort", "day") %in% names(tab))) return(out)
  value_col <- if ("weighted_mean_endpoint" %in% names(tab)) {
    "weighted_mean_endpoint"
  } else if ("weighted_mean_ploidy" %in% names(tab)) {
    "weighted_mean_ploidy"
  } else if ("weighted_mean_N" %in% names(tab)) {
    "weighted_mean_N"
  } else {
    NULL
  }
  if (is.null(value_col)) return(out)

  target_rows <- tab[as.numeric(tab$day) == as.numeric(target_day) & tab$cohort %in% c("2N", "4N"), , drop = FALSE]
  if (!nrow(target_rows)) return(out)
  target_rows[[value_col]] <- suppressWarnings(as.numeric(target_rows[[value_col]]))
  cohort_means <- tapply(target_rows[[value_col]], target_rows$cohort, function(x) mean(x[is.finite(x)], na.rm = TRUE))

  v2 <- suppressWarnings(as.numeric(cohort_means[["2N"]]))
  v4 <- suppressWarnings(as.numeric(cohort_means[["4N"]]))
  if (!is.finite(v2)) v2 <- NA_real_
  if (!is.finite(v4)) v4 <- NA_real_

  out$pred1000_2N <- v2
  out$pred1000_4N <- v4
  out$pred1000_both_gt44 <- isTRUE(is.finite(v2) && is.finite(v4) && v2 > threshold && v4 > threshold)
  out
}

get_recommend_rank_col <- function(summary_df) {
  if ("recommend_rank_burden_ploidy_boundary_first" %in% names(summary_df)) {
    return("recommend_rank_burden_ploidy_boundary_first")
  }
  if ("recommend_rank_ploidy_burden_boundary_first" %in% names(summary_df)) {
    return("recommend_rank_ploidy_burden_boundary_first")
  }
  if ("recommend_rank_ploidy_boundary_first" %in% names(summary_df)) {
    return("recommend_rank_ploidy_boundary_first")
  }
  if ("recommend_rank_ploidy_first" %in% names(summary_df)) {
    return("recommend_rank_ploidy_first")
  }
  if ("objective_burden_rank" %in% names(summary_df)) {
    return("objective_burden_rank")
  }
  if ("objective_ploidy_rank" %in% names(summary_df)) {
    return("objective_ploidy_rank")
  }
  "objective"
}

get_top_ranked_seeds <- function(summary_df, n = 3L, rank_col = NULL, eligible_mask = NULL) {
  if (is.null(rank_col) || !nzchar(rank_col)) rank_col <- get_recommend_rank_col(summary_df)
  keep <- rep(TRUE, nrow(summary_df))
  if (!is.null(eligible_mask)) {
    keep <- keep & !is.na(eligible_mask) & as.logical(eligible_mask)
  }
  candidate_df <- summary_df[keep, , drop = FALSE]
  if (!nrow(candidate_df)) return(character(0))
  ord <- order(candidate_df[[rank_col]], candidate_df$objective, candidate_df$seed, na.last = TRUE)
  as.character(utils::head(candidate_df$seed[ord], as.integer(n)))
}

plot_parameter_boundary_forest <- function(long_df, summary_df, out_path, run_label, near_thresh = 0.05, top3_seeds = NULL, title_suffix = NULL, legend_title = "Recommended Top 3 Seeds") {
  plot_df <- long_df[long_df$active_in_fit & is.finite(long_df$rel_pos_plot), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))

  if (is.null(top3_seeds)) top3_seeds <- get_top_ranked_seeds(summary_df, n = 3L)
  top3_seeds <- as.character(top3_seeds)

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
    scale_x_continuous(limits = c(0, 1), breaks = c(0, near_thresh, 0.5, 1 - near_thresh, 1)) +
    labs(
      title = paste0("Parameter Positions Within Fitted Bounds", if (!is.null(title_suffix) && nzchar(title_suffix)) paste0(" (", title_suffix, ")") else "", ": ", run_label),
      subtitle = paste0("0 = lower bound, 1 = upper bound; shaded zones are within ", sprintf("%.0f", 100 * near_thresh), "% of a bound"),
      x = "Relative position in transformed fit range",
      y = NULL,
      shape = legend_title
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  if (nrow(top_df)) {
    p <- p +
      geom_point(
        data = top_df,
        aes(shape = seed_marker),
        size = 3.0,
        color = "black",
        position = point_pos
      ) +
      scale_shape_manual(values = shape_values, breaks = top_breaks, drop = FALSE)
  }

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

  p <- ggplot(plot_df, aes(x = objective, y = min_rel_dist_active_excl_sigma_burden, label = seed)) +
    geom_point(size = 2.8, color = "#2c7fb8") +
    geom_text(nudge_y = 0.015, show.legend = FALSE, size = 3) +
    scale_y_log10() +
    labs(
      title = paste0("Objective vs Boundary Risk: ", run_label),
      subtitle = "Boundary risk shown as minimum relative distance to a fitted bound, excluding sigma_burden",
      x = "Objective",
      y = "Min relative distance to nearest bound (log10 scale)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  ggplot2::ggsave(out_path, p, width = 10, height = 7)
  invisible(out_path)
}

read_prediction_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
}

first_existing_col <- function(tab, candidates) {
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit)) hit[[1]] else NA_character_
}

prediction_seed_day_mean <- function(tab, seed, value_col, value_name = "value") {
  if (is.null(tab) || !nrow(tab) || is.na(value_col) || !nzchar(value_col)) {
    return(data.frame())
  }
  if (!all(c("cohort", "day", value_col) %in% names(tab))) {
    return(data.frame())
  }
  value <- suppressWarnings(as.numeric(tab[[value_col]]))
  day <- suppressWarnings(as.numeric(tab$day))
  cohort <- as.character(tab$cohort)
  keep <- cohort %in% c("2N", "4N") & is.finite(day) & is.finite(value)
  if (!any(keep)) return(data.frame())
  raw_df <- data.frame(
    seed = as.character(seed),
    cohort = cohort[keep],
    day = day[keep],
    value = value[keep],
    stringsAsFactors = FALSE
  )
  agg <- stats::aggregate(
    value ~ seed + cohort + day,
    data = raw_df,
    FUN = function(x) mean(x[is.finite(x)], na.rm = TRUE)
  )
  names(agg)[names(agg) == "value"] <- value_name
  agg[order(seed_order_key(agg$seed), agg$seed, agg$cohort, agg$day), , drop = FALSE]
}

summarise_seed_prediction_ci <- function(seed_day_df, value_col = "value", group_cols = c("cohort", "day")) {
  if (is.null(seed_day_df) || !nrow(seed_day_df) || !(value_col %in% names(seed_day_df))) {
    return(data.frame())
  }
  seed_day_df[[value_col]] <- suppressWarnings(as.numeric(seed_day_df[[value_col]]))
  keep <- is.finite(seed_day_df[[value_col]])
  for (col in group_cols) {
    keep <- keep & !is.na(seed_day_df[[col]])
  }
  seed_day_df <- seed_day_df[keep, , drop = FALSE]
  if (!nrow(seed_day_df)) return(data.frame())

  group_key <- do.call(
    interaction,
    c(seed_day_df[group_cols], list(drop = TRUE, sep = "\r"))
  )
  parts <- split(seq_len(nrow(seed_day_df)), group_key)
  rows <- lapply(parts, function(idx) {
    sub <- seed_day_df[idx, , drop = FALSE]
    vals <- sub[[value_col]]
    n_seed <- length(vals)
    mean_value <- mean(vals, na.rm = TRUE)
    sd_value <- if (n_seed > 1L) stats::sd(vals, na.rm = TRUE) else 0
    se_value <- if (n_seed > 1L) sd_value / sqrt(n_seed) else 0
    row <- sub[1L, group_cols, drop = FALSE]
    row$n_seed <- n_seed
    row$mean_value <- mean_value
    row$median_value <- stats::median(vals, na.rm = TRUE)
    row$sd_value <- sd_value
    row$se_value <- se_value
    row$ci_low <- mean_value - 1.96 * se_value
    row$ci_high <- mean_value + 1.96 * se_value
    row$min_value <- min(vals, na.rm = TRUE)
    row$max_value <- max(vals, na.rm = TRUE)
    row
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  if ("day" %in% names(out)) out$day <- suppressWarnings(as.numeric(out$day))
  out[do.call(order, c(out[group_cols], list(na.last = TRUE))), , drop = FALSE]
}

prediction_y_limits <- function(summary_df, include_zero = TRUE) {
  value_cols <- intersect(c("ci_low", "ci_high", "min_value", "max_value", "mean_value", "median_value", "ploidy_value"), names(summary_df))
  vals <- suppressWarnings(as.numeric(unlist(summary_df[value_cols], use.names = FALSE)))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(NULL)
  lower <- min(vals, na.rm = TRUE)
  upper <- max(vals, na.rm = TRUE)
  if (isTRUE(include_zero)) lower <- min(0, lower)
  if (!is.finite(lower) || !is.finite(upper)) return(NULL)
  if (identical(lower, upper)) {
    pad <- max(abs(upper) * 0.05, 1)
  } else {
    pad <- 0.04 * (upper - lower)
  }
  c(lower - pad, upper + pad)
}

plot_prediction_mean_ci <- function(summary_df, out_path, cohort, title, subtitle, y_label, color = "#2c7fb8", y_limits = NULL) {
  if (is.null(summary_df) || !nrow(summary_df)) return(invisible(NULL))
  plot_df <- summary_df[as.character(summary_df$cohort) == as.character(cohort), , drop = FALSE]
  plot_df <- plot_df[
    is.finite(plot_df$day) &
      is.finite(plot_df$mean_value) &
      is.finite(plot_df$ci_low) &
      is.finite(plot_df$ci_high),
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))
  plot_df <- plot_df[order(plot_df$day), , drop = FALSE]
  p <- ggplot(plot_df, aes(x = day, y = mean_value)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = color, alpha = 0.22) +
    geom_line(aes(y = min_value), color = color, linewidth = 0.55, linetype = "dashed", alpha = 0.78) +
    geom_line(aes(y = max_value), color = color, linewidth = 0.55, linetype = "dashed", alpha = 0.78) +
    geom_line(color = color, linewidth = 0.9) +
    coord_cartesian(xlim = range(plot_df$day, na.rm = TRUE), ylim = y_limits) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Day",
      y = y_label
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 10, height = 7)
  invisible(out_path)
}

plot_ploidy_prediction_mean_ci_combined <- function(summary_df, out_path, run_label, n_unit = 22) {
  if (is.null(summary_df) || !nrow(summary_df)) return(invisible(NULL))
  plot_df <- summary_df[as.character(summary_df$cohort) %in% c("2N", "4N"), , drop = FALSE]
  plot_df <- plot_df[
    is.finite(plot_df$day) &
      is.finite(plot_df$mean_value) &
      is.finite(plot_df$ci_low) &
      is.finite(plot_df$ci_high) &
      is.finite(plot_df$min_value) &
      is.finite(plot_df$max_value),
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))
  plot_df$cohort <- factor(as.character(plot_df$cohort), levels = c("2N", "4N"))
  plot_df <- plot_df[order(plot_df$cohort, plot_df$day), , drop = FALSE]
  y_limits <- prediction_y_limits(plot_df, include_zero = FALSE)
  colors <- c("2N" = "#1f77b4", "4N" = "#d62728")

  p <- ggplot(plot_df, aes(x = day, y = mean_value, color = cohort, fill = cohort)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.18, color = NA) +
    geom_line(aes(y = min_value), linewidth = 0.55, linetype = "dashed", alpha = 0.75) +
    geom_line(aes(y = max_value), linewidth = 0.55, linetype = "dashed", alpha = 0.75) +
    geom_line(linewidth = 0.95) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_y_continuous(
      name = "Weighted mean chromosome number",
      sec.axis = sec_axis(~ . / n_unit, name = "Ploidy")
    ) +
    coord_cartesian(xlim = range(plot_df$day, na.rm = TRUE), ylim = y_limits) +
    labs(
      title = "1000-day Ploidy Prediction Mean with 95% CI: 2N and 4N",
      subtitle = paste0("Solid lines = cross-seed mean; ribbons = 95% CI; dashed lines = cross-seed min/max envelope | run=", run_label),
      x = "Day",
      color = "Cohort",
      fill = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 11, height = 7)
  invisible(out_path)
}

plot_ploidy_seed_curves_by_cohort <- function(seed_day_df, summary_df, out_path, cohort, run_label, n_unit = 22) {
  if (is.null(seed_day_df) || !nrow(seed_day_df) || is.null(summary_df) || !nrow(summary_df)) {
    return(invisible(NULL))
  }
  curve_df <- seed_day_df[as.character(seed_day_df$cohort) == as.character(cohort), , drop = FALSE]
  curve_df <- curve_df[is.finite(curve_df$day) & is.finite(curve_df$ploidy_value), , drop = FALSE]
  summary_plot_df <- summary_df[as.character(summary_df$cohort) == as.character(cohort), , drop = FALSE]
  summary_plot_df <- summary_plot_df[
    is.finite(summary_plot_df$day) &
      is.finite(summary_plot_df$mean_value) &
      is.finite(summary_plot_df$median_value),
    ,
    drop = FALSE
  ]
  if (!nrow(curve_df) || !nrow(summary_plot_df)) return(invisible(NULL))
  curve_df <- curve_df[order(seed_order_key(curve_df$seed), curve_df$seed, curve_df$day), , drop = FALSE]
  summary_plot_df <- summary_plot_df[order(summary_plot_df$day), , drop = FALSE]
  color <- c("2N" = "#1f77b4", "4N" = "#d62728")[[as.character(cohort)]]
  if (is.null(color) || is.na(color)) color <- "#2c7fb8"
  y_limits <- prediction_y_limits(
    rbind(
      data.frame(ploidy_value = curve_df$ploidy_value, stringsAsFactors = FALSE),
      data.frame(ploidy_value = summary_plot_df$mean_value, stringsAsFactors = FALSE),
      data.frame(ploidy_value = summary_plot_df$median_value, stringsAsFactors = FALSE)
    ),
    include_zero = FALSE
  )

  p <- ggplot() +
    geom_line(
      data = curve_df,
      aes(x = day, y = ploidy_value, group = seed),
      color = "grey72",
      linewidth = 0.18,
      alpha = 0.55
    ) +
    geom_line(
      data = summary_plot_df,
      aes(x = day, y = mean_value),
      color = color,
      linewidth = 1.0,
      linetype = "solid"
    ) +
    geom_line(
      data = summary_plot_df,
      aes(x = day, y = median_value),
      color = color,
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    scale_y_continuous(
      name = "Weighted mean chromosome number",
      sec.axis = sec_axis(~ . / n_unit, name = "Ploidy")
    ) +
    coord_cartesian(xlim = range(curve_df$day, na.rm = TRUE), ylim = y_limits) +
    labs(
      title = paste0("1000-day Ploidy Seed Trajectories: ", cohort),
      subtitle = paste0("Grey hairlines = individual seed trajectories; solid = cross-seed mean; dashed = cross-seed median | run=", run_label),
      x = "Day"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 11, height = 7)
  invisible(out_path)
}

cohort_file_tag <- function(cohort) {
  gsub("[^A-Za-z0-9]+", "", as.character(cohort))
}

extract_total_burden <- function(tab) {
  if (is.null(tab) || !nrow(tab)) return(data.frame())
  if (!all(c("cohort", "day") %in% names(tab))) return(data.frame())
  value_col <- first_existing_col(tab, c("pred_burden_volume_mm3", "pred_burden", "burden_total"))
  if (is.na(value_col)) return(data.frame())
  data.frame(
    cohort = as.character(tab$cohort),
    day = suppressWarnings(as.numeric(tab$day)),
    burden_value = suppressWarnings(as.numeric(tab[[value_col]])),
    stringsAsFactors = FALSE
  )
}

read_ploidy_seed_day_predictions <- function(seed_dirs, horizon_tag = "0_1000day") {
  rows <- vector("list", length(seed_dirs))
  for (i in seq_along(seed_dirs)) {
    seed_dir <- seed_dirs[[i]]
    seed <- basename(seed_dir)
    path <- file.path(seed_dir, "viz", paste0("predict_ploidy_weighted_mean_", horizon_tag, ".tsv"))
    tab <- read_prediction_tsv(path)
    if (is.null(tab)) next
    value_col <- first_existing_col(tab, c("weighted_mean_endpoint", "weighted_mean_N", "weighted_mean_ploidy"))
    rows[[i]] <- prediction_seed_day_mean(tab, seed = seed, value_col = value_col, value_name = "ploidy_value")
  }
  rows <- rows[vapply(rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

read_burden_seed_day_predictions <- function(seed_dirs, horizon_tag = "0_1000day") {
  rows <- list()
  k <- 0L
  for (seed_dir in seed_dirs) {
    seed <- basename(seed_dir)
    path <- file.path(seed_dir, "viz", paste0("predict_burden_", horizon_tag, ".tsv"))
    tab <- read_prediction_tsv(path)
    burden_df <- extract_total_burden(tab)
    if (!nrow(burden_df)) next
    k <- k + 1L
    seed_burden <- prediction_seed_day_mean(
      tab = burden_df,
      seed = seed,
      value_col = "burden_value",
      value_name = "burden_value"
    )
    if (!nrow(seed_burden)) next
    rows[[k]] <- seed_burden[, c("seed", "cohort", "day", "burden_value"), drop = FALSE]
  }
  rows <- rows[vapply(rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

plot_burden_seed_curves_by_cohort <- function(seed_day_df, summary_df, out_path, cohort, run_label, y_limits = NULL) {
  if (is.null(seed_day_df) || !nrow(seed_day_df) || is.null(summary_df) || !nrow(summary_df)) {
    return(invisible(NULL))
  }
  curve_df <- seed_day_df[as.character(seed_day_df$cohort) == as.character(cohort), , drop = FALSE]
  curve_df <- curve_df[is.finite(curve_df$day) & is.finite(curve_df$burden_value), , drop = FALSE]
  summary_plot_df <- summary_df[as.character(summary_df$cohort) == as.character(cohort), , drop = FALSE]
  summary_plot_df <- summary_plot_df[
    is.finite(summary_plot_df$day) &
      is.finite(summary_plot_df$mean_value) &
      is.finite(summary_plot_df$median_value),
    ,
    drop = FALSE
  ]
  if (!nrow(curve_df) || !nrow(summary_plot_df)) return(invisible(NULL))
  curve_df <- curve_df[order(seed_order_key(curve_df$seed), curve_df$seed, curve_df$day), , drop = FALSE]
  summary_plot_df <- summary_plot_df[order(summary_plot_df$day), , drop = FALSE]
  color <- c("2N" = "#1f77b4", "4N" = "#d62728")[[as.character(cohort)]]
  if (is.null(color) || is.na(color)) color <- "#2c7fb8"

  p <- ggplot() +
    geom_line(
      data = curve_df,
      aes(x = day, y = burden_value, group = seed),
      color = "grey72",
      linewidth = 0.18,
      alpha = 0.55
    ) +
    geom_line(
      data = summary_plot_df,
      aes(x = day, y = mean_value),
      color = color,
      linewidth = 1.0,
      linetype = "solid"
    ) +
    geom_line(
      data = summary_plot_df,
      aes(x = day, y = median_value),
      color = color,
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    coord_cartesian(xlim = range(curve_df$day, na.rm = TRUE), ylim = y_limits) +
    labs(
      title = paste0("1000-day Total Burden Seed Trajectories: ", cohort),
      subtitle = paste0("Grey hairlines = individual seed trajectories; solid = cross-seed mean; dashed = cross-seed median | run=", run_label),
      x = "Day",
      y = "Total tumor burden (mm^3)"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 11, height = 7)
  invisible(out_path)
}

cleanup_stale_burden_component_outputs <- function(out_dir) {
  stale <- c(
    "predict1000_burden_components_seed_day_mean.tsv",
    "predict1000_burden_components_mean_ci.tsv",
    "predict1000_burden_total_mean_ci_2N_4N.pdf"
  )
  stale <- c(
    stale,
    list.files(
      out_dir,
      pattern = "^predict1000_burden_(live|dead_hypoxia|dead_buffer|dead_total)_mean_ci_(2N|4N)\\.pdf$",
      full.names = FALSE
    )
  )
  unlink(file.path(out_dir, stale), force = TRUE)
  invisible(TRUE)
}

plot_prediction_summaries <- function(seed_dirs, out_dir, run_label, horizon_tag = "0_1000day") {
  written <- character(0)
  colors <- c("2N" = "#1f77b4", "4N" = "#d62728")

  ploidy_seed_day <- read_ploidy_seed_day_predictions(seed_dirs, horizon_tag = horizon_tag)
  if (nrow(ploidy_seed_day)) {
    ploidy_summary <- summarise_seed_prediction_ci(
      seed_day_df = ploidy_seed_day,
      value_col = "ploidy_value",
      group_cols = c("cohort", "day")
    )
    utils::write.table(
      ploidy_seed_day,
      file = file.path(out_dir, "predict1000_ploidy_seed_day_mean.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    utils::write.table(
      ploidy_summary,
      file = file.path(out_dir, "predict1000_ploidy_mean_ci.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    written <- c(
      written,
      file.path(out_dir, "predict1000_ploidy_seed_day_mean.tsv"),
      file.path(out_dir, "predict1000_ploidy_mean_ci.tsv")
    )
    unlink(file.path(out_dir, paste0("predict1000_ploidy_mean_ci_", c("2N", "4N"), ".pdf")), force = TRUE)
    out_path <- file.path(out_dir, "predict1000_ploidy_mean_ci_2N_4N.pdf")
    res <- plot_ploidy_prediction_mean_ci_combined(
      summary_df = ploidy_summary,
      out_path = out_path,
      run_label = run_label
    )
    if (!is.null(res)) written <- c(written, out_path)
    for (cohort in c("2N", "4N")) {
      out_path <- file.path(out_dir, paste0("predict1000_ploidy_seed_curves_", cohort_file_tag(cohort), ".pdf"))
      res <- plot_ploidy_seed_curves_by_cohort(
        seed_day_df = ploidy_seed_day,
        summary_df = ploidy_summary,
        out_path = out_path,
        cohort = cohort,
        run_label = run_label
      )
      if (!is.null(res)) written <- c(written, out_path)
    }
  }

  burden_seed_day <- read_burden_seed_day_predictions(seed_dirs, horizon_tag = horizon_tag)
  cleanup_stale_burden_component_outputs(out_dir)
  if (nrow(burden_seed_day)) {
    burden_summary <- summarise_seed_prediction_ci(
      seed_day_df = burden_seed_day,
      value_col = "burden_value",
      group_cols = c("cohort", "day")
    )
    utils::write.table(
      burden_seed_day,
      file = file.path(out_dir, "predict1000_burden_total_seed_day_mean.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    utils::write.table(
      burden_summary,
      file = file.path(out_dir, "predict1000_burden_total_mean_ci.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    written <- c(
      written,
      file.path(out_dir, "predict1000_burden_total_seed_day_mean.tsv"),
      file.path(out_dir, "predict1000_burden_total_mean_ci.tsv")
    )
    for (cohort in c("2N", "4N")) {
      cohort_summary <- burden_summary[as.character(burden_summary$cohort) == cohort, , drop = FALSE]
      burden_y_limits <- prediction_y_limits(cohort_summary, include_zero = TRUE)

      out_path <- file.path(out_dir, paste0("predict1000_burden_total_mean_ci_", cohort_file_tag(cohort), ".pdf"))
      res <- plot_prediction_mean_ci(
        summary_df = burden_summary,
        out_path = out_path,
        cohort = cohort,
        title = paste0("1000-day Total Burden Prediction Mean with 95% CI: ", cohort),
        subtitle = paste0("Seed-level scenario means aggregated across seeds; dashed lines = min/max envelope | run=", run_label),
        y_label = "Total tumor burden (mm^3)",
        color = colors[[cohort]],
        y_limits = burden_y_limits
      )
      if (!is.null(res)) written <- c(written, out_path)

      out_path <- file.path(out_dir, paste0("predict1000_burden_total_seed_curves_", cohort_file_tag(cohort), ".pdf"))
      res <- plot_burden_seed_curves_by_cohort(
        seed_day_df = burden_seed_day,
        summary_df = burden_summary,
        out_path = out_path,
        cohort = cohort,
        run_label = run_label,
        y_limits = burden_y_limits
      )
      if (!is.null(res)) written <- c(written, out_path)
    }
  }

  unique(written)
}

run_rscript_helper <- function(script_path, args, label) {
  if (!file.exists(script_path)) {
    stop("Missing helper script for ", label, ": ", script_path)
  }
  output <- system2(
    "Rscript",
    args = c(normalizePath(script_path, mustWork = TRUE), args),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) {
    detail <- if (length(output)) paste(utils::tail(output, 20L), collapse = "\n") else "(no helper output)"
    stop(label, " failed with exit status ", status, ". Last output:\n", detail)
  }
  if (length(output)) {
    message(paste(output, collapse = "\n"))
  }
  invisible(TRUE)
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
    pred_gate_metrics <- read_1000day_ploidy_gate_metrics(seed_dir)
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
      parameter_long = long_df,
      pred_gate_metrics = pred_gate_metrics
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
  seed_summary$pred1000_2N <- suppressWarnings(as.numeric(seed_summary$pred1000_2N))
  seed_summary$pred1000_4N <- suppressWarnings(as.numeric(seed_summary$pred1000_4N))
  seed_summary$pred1000_both_gt44 <- as.logical(seed_summary$pred1000_both_gt44)
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
  forest_rank_col <- get_recommend_rank_col(seed_summary)
  forest_rank_simple <- suppressWarnings(as.integer(seed_summary[[forest_rank_col]]))
  forest_rank_plus_ploidy_simple <- rep(NA_integer_, nrow(seed_summary))
  eligible_plot_idx <- which(!is.na(seed_summary$pred1000_both_gt44) & seed_summary$pred1000_both_gt44)
  if (length(eligible_plot_idx) > 0L) {
    eligible_plot_ord <- eligible_plot_idx[order(
      seed_summary[[forest_rank_col]][eligible_plot_idx],
      seed_summary$objective[eligible_plot_idx],
      seed_summary$seed[eligible_plot_idx],
      na.last = TRUE
    )]
    forest_rank_plus_ploidy_simple[eligible_plot_ord] <- seq_along(eligible_plot_ord)
  }
  pred_gate_top3_seeds <- get_top_ranked_seeds(
    summary_df = seed_summary,
    n = 3L,
    eligible_mask = seed_summary$pred1000_both_gt44
  )
  seed_summary$forest_plot_rank_simple <- forest_rank_simple
  seed_summary$forest_plot_rank_plus_ploidy_simple <- forest_rank_plus_ploidy_simple
  seed_summary <- seed_summary[order(seed_summary$objective, seed_summary$seed), , drop = FALSE]
  row.names(seed_summary) <- NULL

  objective_simple <- seed_summary[, c("seed", "objective", "objective_burden", "objective_ploidy"), drop = FALSE]
  objective_simple$objective_rank <- suppressWarnings(as.integer(seed_summary$forest_plot_rank_simple))
  objective_simple$objective_rank_plus_ploidy <- suppressWarnings(as.integer(seed_summary$forest_plot_rank_plus_ploidy_simple))
  objective_simple <- objective_simple[, c("seed", "objective_rank", "objective_rank_plus_ploidy", "objective", "objective_burden", "objective_ploidy"), drop = FALSE]
  objective_simple <- objective_simple[order(objective_simple$objective_rank, objective_simple$objective, objective_simple$seed, na.last = TRUE), , drop = FALSE]
  row.names(objective_simple) <- NULL

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
  utils::write.table(
    objective_simple,
    file = file.path(out_dir, "seed_objective_simple.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  forest_out <- plot_parameter_boundary_forest(
    long_df = parameter_long,
    summary_df = seed_summary,
    out_path = file.path(out_dir, "parameter_boundary_forest.pdf"),
    run_label = basename(run_dir),
    near_thresh = near_thresh
  )
  objective_risk_out <- plot_objective_vs_boundary_risk(
    summary_df = seed_summary,
    out_path = file.path(out_dir, "objective_vs_boundary_risk.pdf"),
    run_label = basename(run_dir)
  )
  forest_filtered_out <- plot_parameter_boundary_forest(
    long_df = parameter_long,
    summary_df = seed_summary,
    out_path = file.path(out_dir, "parameter_boundary_forest_pred1000_gt44_top3.pdf"),
    run_label = basename(run_dir),
    near_thresh = near_thresh,
    top3_seeds = pred_gate_top3_seeds,
    title_suffix = "Top 3 among seeds with 2N/4N 1000d predictions > 44",
    legend_title = "Top 3 Seeds with 2N/4N 1000d > 44"
  )
  prediction_outputs <- plot_prediction_summaries(
    seed_dirs = seed_dirs,
    out_dir = out_dir,
    run_label = basename(run_dir),
    horizon_tag = "0_1000day"
  )

  objective_violin_script <- normalizePath(file.path(SCRIPT_DIR, "plot_extra_results_objective_violin.R"), mustWork = FALSE)
  extra_results_report_script <- normalizePath(file.path(SCRIPT_DIR, "extra_results_report.R"), mustWork = FALSE)
  run_rscript_helper(
    script_path = objective_violin_script,
    args = c(paste0("--extra_results_dir=", out_dir)),
    label = "plot_extra_results_objective_violin"
  )
  run_rscript_helper(
    script_path = extra_results_report_script,
    args = c(paste0("--extra_results_dir=", out_dir)),
    label = "extra_results_report"
  )

  message("Wrote summary table: ", file.path(out_dir, "seed_summary.tsv"))
  message("Wrote parameter long table: ", file.path(out_dir, "parameter_boundary_long.tsv"))
  message("Wrote objective simple table: ", file.path(out_dir, "seed_objective_simple.tsv"))
  if (!is.null(forest_out) && file.exists(forest_out)) {
    message("Wrote forest plot: ", forest_out)
  } else {
    message("Skipped forest plot because no active fitted parameters were available.")
  }
  if (!is.null(forest_filtered_out) && file.exists(forest_filtered_out)) {
    message("Wrote filtered forest plot: ", forest_filtered_out)
  } else {
    message("Skipped filtered forest plot because no eligible plotted parameters were available.")
  }
  if (!is.null(objective_risk_out) && file.exists(objective_risk_out)) {
    message("Wrote objective-risk plot: ", objective_risk_out)
  } else {
    message("Skipped objective-risk plot because no seeds had a finite positive boundary-distance metric.")
  }
  message("Wrote objective-components violin: ", file.path(out_dir, "objective_components_violin.pdf"))
  if (length(prediction_outputs)) {
    message("Wrote cross-seed prediction outputs: ", paste(basename(prediction_outputs), collapse = ", "))
  }
  message("Wrote extra results report: ", file.path(out_dir, "extra_results_report.html"))
}

if (sys.nframe() == 0) {
  main()
}
