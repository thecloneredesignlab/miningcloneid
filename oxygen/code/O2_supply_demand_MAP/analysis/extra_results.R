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
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_bool <- o2sd_as_bool

summary_flag_true <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  tolower(trimws(as.character(x[[1]]))) %in% c("true", "t", "1", "yes", "y", "on")
}

summary_flag_na <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(NA)
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

fit_mode_value <- function(fit_summary_vals, default = "fit_invivo") {
  if (is.null(fit_summary_vals) || !length(fit_summary_vals)) return(default)
  inferred_default <- if ("objective_total" %in% names(fit_summary_vals)) {
    "fit_invitro"
  } else if ("objective_invivo" %in% names(fit_summary_vals) || "objective_invitro" %in% names(fit_summary_vals)) {
    "fit_joint"
  } else {
    default
  }
  mode <- summary_metric_value(fit_summary_vals, "fit_mode", inferred_default)
  mode <- trimws(as.character(mode %||% default))
  if (!nzchar(mode)) inferred_default else mode
}

is_joint_fit_summary <- function(fit_summary_vals) {
  identical(fit_mode_value(fit_summary_vals), "fit_joint")
}

is_invitro_fit_summary <- function(fit_summary_vals) {
  identical(fit_mode_value(fit_summary_vals), "fit_invitro")
}

filter_best_vals_for_output <- function(best_vals, fit_summary_vals) {
  harvest_use <- summary_flag_true(summary_metric_value(fit_summary_vals, "harvest_init_multiplier", FALSE), default = FALSE)
  out <- filter_family_specific_run_params_for_output_common(run_params = best_vals)
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

read_metric_map_optional <- function(path, key_col, value_col) {
  if (!file.exists(path)) return(character(0))
  tryCatch(
    read_metric_map(path, key_col = key_col, value_col = value_col),
    error = function(e) character(0)
  )
}

supplement_joint_invitro_metrics <- function(fit_summary_vals, seed_dir) {
  if (!is_joint_fit_summary(fit_summary_vals)) return(fit_summary_vals)
  joint_components <- read_metric_map_optional(file.path(seed_dir, "joint_components.tsv"), "component", "value")
  if (!length(joint_components)) return(fit_summary_vals)
  metric_map <- c(
    growth_loglik = "invitro_growth_loglik",
    ploidy_loglik = "invitro_ploidy_loglik",
    flow_loglik = "invitro_flow_loglik"
  )
  for (metric in names(metric_map)) {
    src <- metric_map[[metric]]
    if (src %in% names(joint_components) && (!metric %in% names(fit_summary_vals) || !is.finite(as_num(fit_summary_vals[[metric]], NA_real_)))) {
      fit_summary_vals[[metric]] <- joint_components[[src]]
    }
  }
  passthrough_metrics <- c(
    "joint_constraint_penalty",
    "joint_constraint_failures",
    "joint_constraint_penalty_total",
    "joint_constraints_pass",
    "objective_soft_coupling",
    "objective_constraints",
    "joint_soft_coupling_enabled",
    "joint_soft_coupling_params",
    "joint_soft_coupling_n_params",
    "joint_soft_coupling_sigma_default",
    "joint_restriction",
    "joint_require_invivo_pred1000_ploidy_gt2",
    "joint_require_invitro_growth_nonnegative",
    "joint_require_invitro_ploidy_phenotype",
    "invitro_growth_nonnegative_pass",
    "invitro_n_growth_negative_pred",
    "invivo_pred1000_ploidy_pass",
    "invivo_pred1000_min_ploidy_fold",
    "invivo_pred1000_threshold_ploidy_fold",
    "invivo_pred1000_n_rows",
    "invivo_pred1000_status",
    "invitro_ploidy_phenotype_pass",
    "invitro_2N_deprived_wgd_pass",
    "invitro_2N_deprived_wgd_max_fraction",
    "invitro_2N_deprived_wgd_min_N",
    "invitro_2N_deprived_wgd_min_fraction",
    "invitro_2N_deprived_max_mean_N",
    "invitro_4N_deprived_chr_drop_pass",
    "invitro_4N_deprived_chr_drop",
    "invitro_4N_deprived_min_chr_drop_required",
    "invitro_4N_deprived_initial_mean_N",
    "invitro_4N_deprived_min_mean_N"
  )
  for (metric in passthrough_metrics) {
    if (metric %in% names(joint_components) && (!metric %in% names(fit_summary_vals) || !nzchar(as.character(fit_summary_vals[[metric]])))) {
      fit_summary_vals[[metric]] <- joint_components[[metric]]
    }
  }
  if ("objective_invitro" %in% names(joint_components)) {
    objective_invitro <- as_num(joint_components[["objective_invitro"]], NA_real_)
    if (is.finite(objective_invitro) && (!"total_loglik" %in% names(fit_summary_vals) || !is.finite(as_num(fit_summary_vals[["total_loglik"]], NA_real_)))) {
      fit_summary_vals[["total_loglik"]] <- as.character(-objective_invitro)
    }
  }
  fit_summary_vals
}

bool_from_table_value <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  tolower(trimws(as.character(x[[1]]))) %in% c("true", "t", "1", "yes", "y", "on")
}

invitro_parameter_transform_map <- function() {
  data.frame(
    param_symbol = c(
      "lam_max", "p_misseg", "k_o_mis",
      "buffer_smax", "buffer_beta", "buffer_n_exp",
      "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
      "O2_crit", "n_O", "p_mis_base", "sigma_growth", "sigma_kary",
      "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N"
    ),
    param_name = c(
      "log10_lam_max", "log10_p_misseg", "log10_k_o_mis",
      "buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp",
      "log10_p_wgd", "log10_alpha_o2", "gamma_growth", "log10_mu_hp", "gamma_mu",
      "log10_O2_crit", "n_O", "log10_p_mis_base", "log10_sigma_growth", "log10_sigma_kary",
      "init_mean_2N", "log10_init_sd_2N", "init_mean_4N", "log10_init_sd_4N"
    ),
    transform = c(
      "log10", "log10", "log10",
      "identity", "log10", "log10",
      "log10", "log10", "identity", "log10", "identity",
      "log10", "identity", "log10", "log10", "log10",
      "identity", "log10", "identity", "log10"
    ),
    stringsAsFactors = FALSE
  )
}

transform_invitro_bound <- function(x, transform, param_symbol) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return(NA_real_)
  if (identical(transform, "log10")) {
    if (x <= 0) stop("In vitro parameter '", param_symbol, "' must have positive bounds for log10 transform.")
    return(log10(x))
  }
  if (identical(transform, "logit")) {
    if (x <= 0 || x >= 1) stop("In vitro parameter '", param_symbol, "' must have bounds inside (0, 1) for logit transform.")
    return(qlogis(x))
  }
  x
}

convert_invitro_parameter_table <- function(tab, path) {
  required_cols <- c("param_symbol", "init_value", "lower_bound", "upper_bound")
  if (!all(required_cols %in% names(tab))) {
    stop(basename(path), " is missing required in vitro columns: ", paste(setdiff(required_cols, names(tab)), collapse = ", "))
  }
  if (!("use_invitro_fit" %in% names(tab)) && !("estimate" %in% names(tab))) {
    stop(basename(path), " must contain either 'use_invitro_fit' or 'estimate' for in vitro boundary reporting.")
  }
  tab$param_symbol <- trimws(as.character(tab$param_symbol))
  map <- invitro_parameter_transform_map()
  tab <- merge(map, tab, by = "param_symbol", all.x = FALSE, all.y = FALSE, sort = FALSE)
  tab <- tab[match(map$param_symbol, tab$param_symbol, nomatch = 0L), , drop = FALSE]
  if (!nrow(tab)) {
    stop("No supported in vitro optimizer parameters found in: ", path)
  }
  estimate_col <- if ("use_invitro_fit" %in% names(tab)) "use_invitro_fit" else "estimate"
  estimate_flag <- vapply(tab[[estimate_col]], bool_from_table_value, logical(1))
  lower_t <- mapply(transform_invitro_bound, tab$lower_bound, tab$transform, tab$param_symbol)
  upper_t <- mapply(transform_invitro_bound, tab$upper_bound, tab$transform, tab$param_symbol)
  data.frame(
    param_name = tab$param_name,
    estimate = estimate_flag,
    init_value = mapply(transform_invitro_bound, tab$init_value, tab$transform, tab$param_symbol),
    lower_bound = as.numeric(lower_t),
    upper_bound = as.numeric(upper_t),
    param_prototype = tab$param_symbol,
    prototype_init_value = suppressWarnings(as.numeric(tab$init_value)),
    prototype_lower_bound = suppressWarnings(as.numeric(tab$lower_bound)),
    prototype_upper_bound = suppressWarnings(as.numeric(tab$upper_bound)),
    source = if ("source" %in% names(tab)) as.character(tab$source) else "parameter_table_input",
    note = if ("description" %in% names(tab)) as.character(tab$description) else if ("note" %in% names(tab)) as.character(tab$note) else "",
    stringsAsFactors = FALSE
  )
}

read_parameter_table_checked <- function(path) {
  param_table <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("param_name", "estimate", "lower_bound", "upper_bound", "param_prototype", "prototype_lower_bound", "prototype_upper_bound")
  if (!all(required_cols %in% names(param_table)) && "param_symbol" %in% names(param_table)) {
    param_table <- convert_invitro_parameter_table(param_table, path)
  }
  if (!all(required_cols %in% names(param_table))) {
    stop(basename(path), " is missing required columns: ", paste(setdiff(required_cols, names(param_table)), collapse = ", "))
  }
  param_table$estimate <- vapply(param_table$estimate, function(x) as_bool(x, FALSE), logical(1))
  param_table
}

find_parameter_table_for_seed <- function(run_dir, seed_dir, fit_summary_vals) {
  if (is_joint_fit_summary(fit_summary_vals)) {
    candidates <- c(
      file.path(seed_dir, "parameter_table_invivo_transformed.csv"),
      file.path(run_dir, "parameter_table.csv"),
      file.path(seed_dir, "parameter_table.csv")
    )
  } else if (is_invitro_fit_summary(fit_summary_vals)) {
    candidates <- c(
      file.path(seed_dir, "parameter_table_input.csv"),
      file.path(run_dir, "parameter_table_input.csv"),
      file.path(seed_dir, "parameter_table.csv"),
      file.path(run_dir, "parameter_table.csv")
    )
  } else {
    candidates <- c(
      file.path(run_dir, "parameter_table.csv"),
      file.path(seed_dir, "parameter_table.csv"),
      file.path(seed_dir, "parameter_table_invivo_transformed.csv")
    )
  }
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop(
      "Missing parameter table for ", basename(seed_dir), ". Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[[1]]
}

find_joint_invitro_parameter_table_for_seed <- function(run_dir, seed_dir) {
  candidates <- c(
    file.path(seed_dir, "parameter_table_input_invitro.csv"),
    file.path(run_dir, "parameter_table_input_invitro.csv"),
    file.path(seed_dir, "parameter_table_input.csv"),
    file.path(run_dir, "parameter_table_input.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1]] else NA_character_
}

lookup_named_num <- function(x, key) {
  if (!length(x) || is.null(key) || !nzchar(key) || !(key %in% names(x))) return(NA_real_)
  val <- suppressWarnings(as.numeric(x[[key]]))
  if (!is.finite(val)) NA_real_ else val
}

derive_prototype_value <- function(param_name, param_prototype, best_vals) {
  val <- lookup_named_num(best_vals, param_prototype)
  if (is.finite(val)) return(val)
  NA_real_
}

derive_transformed_value <- function(param_name, param_prototype, best_vals) {
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
  if (identical(param_prototype, "p_wgd")) {
    return(TRUE)
  }
  if (param_prototype %in% c("buffer_smax", "buffer_beta", "buffer_n_exp")) {
    return(TRUE)
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
    parameter_description <- ""
    for (desc_col in intersect(c("parameter_description", "note", "description"), names(param_table))) {
      candidate <- trimws(as.character(param_table[[desc_col]][[i]]))
      if (nzchar(candidate)) {
        parameter_description <- candidate
        break
      }
    }
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
      parameter_description = parameter_description,
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
    fit_mode = fit_mode_value(fit_summary_vals),
    objective = as_num(summary_metric_value(fit_summary_vals, "objective", summary_metric_value(fit_summary_vals, "objective_total", NA_real_))),
    objective_total = as_num(summary_metric_value(fit_summary_vals, "objective_total", NA_real_)),
    objective_data = as_num(summary_metric_value(fit_summary_vals, "objective_data", NA_real_)),
    objective_prior = as_num(summary_metric_value(fit_summary_vals, "objective_prior", NA_real_)),
    objective_burden = as_num(summary_metric_value(fit_summary_vals, "objective_burden", NA_real_)),
    objective_ploidy = as_num(summary_metric_value(fit_summary_vals, "objective_ploidy", NA_real_)),
    objective_invivo = as_num(summary_metric_value(fit_summary_vals, "objective_invivo", NA_real_)),
    objective_invitro = as_num(summary_metric_value(fit_summary_vals, "objective_invitro", NA_real_)),
    objective_soft_coupling = as_num(summary_metric_value(fit_summary_vals, "objective_soft_coupling", NA_real_)),
    objective_constraints = as_num(summary_metric_value(fit_summary_vals, "objective_constraints", NA_real_)),
    total_loglik = as_num(summary_metric_value(fit_summary_vals, "total_loglik", NA_real_)),
    growth_loglik = as_num(summary_metric_value(fit_summary_vals, "growth_loglik", NA_real_)),
    ploidy_loglik = as_num(summary_metric_value(fit_summary_vals, "ploidy_loglik", NA_real_)),
    flow_loglik = as_num(summary_metric_value(fit_summary_vals, "flow_loglik", NA_real_)),
    growth_loglik_sum = as_num(summary_metric_value(fit_summary_vals, "growth_loglik_sum", NA_real_)),
    ploidy_loglik_sum = as_num(summary_metric_value(fit_summary_vals, "ploidy_loglik_sum", NA_real_)),
    flow_loglik_sum = as_num(summary_metric_value(fit_summary_vals, "flow_loglik_sum", NA_real_)),
    sigma_growth = as_num(summary_metric_value(fit_summary_vals, "sigma_growth", NA_real_)),
    sigma_kary = as_num(summary_metric_value(fit_summary_vals, "sigma_kary", NA_real_)),
    sigma_flow_ploidy = as_num(summary_metric_value(fit_summary_vals, "sigma_flow_ploidy", NA_real_)),
    n_growth = as_num(summary_metric_value(fit_summary_vals, "n_growth", NA_real_)),
    n_ploidy_passages = as_num(summary_metric_value(fit_summary_vals, "n_ploidy_passages", NA_real_)),
    n_kary_cells = as_num(summary_metric_value(fit_summary_vals, "n_kary_cells", NA_real_)),
    n_flow_passages = as_num(summary_metric_value(fit_summary_vals, "n_flow_passages", NA_real_)),
    n_flow_samples = as_num(summary_metric_value(fit_summary_vals, "n_flow_samples", NA_real_)),
    joint_weight_invivo = as_num(summary_metric_value(fit_summary_vals, "joint_weight_invivo", NA_real_)),
    joint_weight_invitro = as_num(summary_metric_value(fit_summary_vals, "joint_weight_invitro", NA_real_)),
    joint_invitro_growth_weight = as_num(summary_metric_value(fit_summary_vals, "joint_invitro_growth_weight", NA_real_)),
    joint_invitro_ploidy_weight = as_num(summary_metric_value(fit_summary_vals, "joint_invitro_ploidy_weight", NA_real_)),
    joint_invitro_flow_weight = as_num(summary_metric_value(fit_summary_vals, "joint_invitro_flow_weight", NA_real_)),
    n_cores_requested = as_num(summary_metric_value(fit_summary_vals, "n_cores_requested", NA_real_)),
    n_cores_used = as_num(summary_metric_value(fit_summary_vals, "n_cores_used", NA_real_)),
    n_parameters = as_num(summary_metric_value(fit_summary_vals, "n_parameters", NA_real_)),
    n_invivo_scenarios = as_num(summary_metric_value(fit_summary_vals, "n_invivo_scenarios", NA_real_)),
    objective_burden_neg2loglik_raw = as_num(summary_metric_value(fit_summary_vals, "objective_burden_neg2loglik_raw", NA_real_)),
    objective_ploidy_neg2loglik_raw = as_num(summary_metric_value(fit_summary_vals, "objective_ploidy_neg2loglik_raw", NA_real_)),
    itermax = as_num(summary_metric_value(fit_summary_vals, "itermax", NA_real_)),
    optimizer_deoptim_objective = as_num(summary_metric_value(fit_summary_vals, "optimizer_deoptim_objective", NA_real_)),
    optimizer_local_objective = as_num(summary_metric_value(fit_summary_vals, "optimizer_local_objective", NA_real_)),
    optimizer_local_attempted = summary_flag_na(summary_metric_value(fit_summary_vals, "optimizer_local_attempted", NA)),
    optimizer_local_accepted = summary_flag_na(summary_metric_value(fit_summary_vals, "optimizer_local_accepted", NA)),
    optimizer_local_convergence = as_num(summary_metric_value(fit_summary_vals, "optimizer_local_convergence", NA_real_)),
    optimizer_local_maxit = as_num(summary_metric_value(fit_summary_vals, "optimizer_local_maxit", NA_real_)),
    optimizer_interrupted = as.character(summary_metric_value(fit_summary_vals, "optimizer_interrupted", NA_character_)),
    optimizer_iter_completed = as_num(summary_metric_value(fit_summary_vals, "optimizer_iter_completed", NA_real_)),
    optimizer_iter_target = as_num(summary_metric_value(fit_summary_vals, "optimizer_iter_target", NA_real_)),
    deoptim_stop_reason = as.character(summary_metric_value(fit_summary_vals, "deoptim_stop_reason", NA_character_)),
    joint_constraint_penalty = as_num(summary_metric_value(fit_summary_vals, "joint_constraint_penalty", NA_real_)),
    joint_constraint_failures = as_num(summary_metric_value(fit_summary_vals, "joint_constraint_failures", NA_real_)),
    joint_constraint_penalty_total = as_num(summary_metric_value(fit_summary_vals, "joint_constraint_penalty_total", NA_real_)),
    joint_constraints_pass = summary_flag_na(summary_metric_value(fit_summary_vals, "joint_constraints_pass", NA)),
    joint_soft_coupling_enabled = summary_flag_na(summary_metric_value(fit_summary_vals, "joint_soft_coupling_enabled", NA)),
    joint_soft_coupling_params = as.character(summary_metric_value(fit_summary_vals, "joint_soft_coupling_params", NA_character_)),
    joint_soft_coupling_n_params = as_num(summary_metric_value(fit_summary_vals, "joint_soft_coupling_n_params", NA_real_)),
    joint_soft_coupling_sigma_default = as_num(summary_metric_value(fit_summary_vals, "joint_soft_coupling_sigma_default", NA_real_)),
    joint_restriction = summary_flag_na(summary_metric_value(fit_summary_vals, "joint_restriction", NA)),
    joint_require_invivo_pred1000_ploidy_gt2 = summary_flag_na(summary_metric_value(fit_summary_vals, "joint_require_invivo_pred1000_ploidy_gt2", NA)),
    joint_require_invitro_growth_nonnegative = summary_flag_na(summary_metric_value(fit_summary_vals, "joint_require_invitro_growth_nonnegative", NA)),
    joint_require_invitro_ploidy_phenotype = summary_flag_na(summary_metric_value(fit_summary_vals, "joint_require_invitro_ploidy_phenotype", NA)),
    invitro_growth_nonnegative_pass = summary_flag_na(summary_metric_value(fit_summary_vals, "invitro_growth_nonnegative_pass", NA)),
    invitro_n_growth_negative_pred = as_num(summary_metric_value(fit_summary_vals, "invitro_n_growth_negative_pred", NA_real_)),
    invivo_pred1000_ploidy_pass = summary_flag_na(summary_metric_value(fit_summary_vals, "invivo_pred1000_ploidy_pass", NA)),
    invivo_pred1000_min_ploidy_fold = as_num(summary_metric_value(fit_summary_vals, "invivo_pred1000_min_ploidy_fold", NA_real_)),
    invivo_pred1000_threshold_ploidy_fold = as_num(summary_metric_value(fit_summary_vals, "invivo_pred1000_threshold_ploidy_fold", NA_real_)),
    invivo_pred1000_n_rows = as_num(summary_metric_value(fit_summary_vals, "invivo_pred1000_n_rows", NA_real_)),
    invivo_pred1000_status = as.character(summary_metric_value(fit_summary_vals, "invivo_pred1000_status", NA_character_)),
    invitro_ploidy_phenotype_pass = summary_flag_na(summary_metric_value(fit_summary_vals, "invitro_ploidy_phenotype_pass", NA)),
    invitro_2N_deprived_wgd_pass = summary_flag_na(summary_metric_value(fit_summary_vals, "invitro_2N_deprived_wgd_pass", NA)),
    invitro_2N_deprived_wgd_max_fraction = as_num(summary_metric_value(fit_summary_vals, "invitro_2N_deprived_wgd_max_fraction", NA_real_)),
    invitro_2N_deprived_wgd_min_N = as_num(summary_metric_value(fit_summary_vals, "invitro_2N_deprived_wgd_min_N", NA_real_)),
    invitro_2N_deprived_wgd_min_fraction = as_num(summary_metric_value(fit_summary_vals, "invitro_2N_deprived_wgd_min_fraction", NA_real_)),
    invitro_2N_deprived_max_mean_N = as_num(summary_metric_value(fit_summary_vals, "invitro_2N_deprived_max_mean_N", NA_real_)),
    invitro_4N_deprived_chr_drop_pass = summary_flag_na(summary_metric_value(fit_summary_vals, "invitro_4N_deprived_chr_drop_pass", NA)),
    invitro_4N_deprived_chr_drop = as_num(summary_metric_value(fit_summary_vals, "invitro_4N_deprived_chr_drop", NA_real_)),
    invitro_4N_deprived_min_chr_drop_required = as_num(summary_metric_value(fit_summary_vals, "invitro_4N_deprived_min_chr_drop_required", NA_real_)),
    invitro_4N_deprived_initial_mean_N = as_num(summary_metric_value(fit_summary_vals, "invitro_4N_deprived_initial_mean_N", NA_real_)),
    invitro_4N_deprived_min_mean_N = as_num(summary_metric_value(fit_summary_vals, "invitro_4N_deprived_min_mean_N", NA_real_)),
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

seed_prediction_path <- function(seed_dir, filename) {
  candidates <- c(
    file.path(seed_dir, "viz", filename),
    file.path(seed_dir, "viz", "invivo", filename)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1]] else candidates[[1]]
}

read_1000day_ploidy_gate_metrics <- function(seed_dir, target_day = 1000, threshold = 44) {
  path <- seed_prediction_path(seed_dir, "predict_ploidy_weighted_mean_0_1000day.tsv")
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

  tab$day <- suppressWarnings(as.numeric(tab$day))
  tab <- tab[is.finite(tab$day) & tab$cohort %in% c("2N", "4N"), , drop = FALSE]
  if (!nrow(tab)) return(out)

  target_day <- suppressWarnings(as.numeric(target_day))
  target_rows <- tab[abs(tab$day - target_day) <= 1e-8, , drop = FALSE]
  if (!nrow(target_rows)) {
    day_use <- suppressWarnings(max(tab$day[tab$day <= target_day], na.rm = TRUE))
    if (!is.finite(day_use)) day_use <- suppressWarnings(max(tab$day, na.rm = TRUE))
    if (!is.finite(day_use)) return(out)
    target_rows <- tab[abs(tab$day - day_use) <= 1e-8, , drop = FALSE]
  }
  if (!nrow(target_rows)) return(out)
  target_rows[[value_col]] <- suppressWarnings(as.numeric(target_rows[[value_col]]))
  cohort_values <- lapply(c("2N", "4N"), function(cohort) {
    vals <- target_rows[[value_col]][as.character(target_rows$cohort) == cohort]
    vals[is.finite(vals)]
  })
  names(cohort_values) <- c("2N", "4N")

  v2 <- if (length(cohort_values[["2N"]])) min(cohort_values[["2N"]], na.rm = TRUE) else NA_real_
  v4 <- if (length(cohort_values[["4N"]])) min(cohort_values[["4N"]], na.rm = TRUE) else NA_real_
  if (!is.finite(v2)) v2 <- NA_real_
  if (!is.finite(v4)) v4 <- NA_real_

  out$pred1000_2N <- v2
  out$pred1000_4N <- v4
  out$pred1000_both_gt44 <- isTRUE(
    length(cohort_values[["2N"]]) > 0L &&
      length(cohort_values[["4N"]]) > 0L &&
      all(cohort_values[["2N"]] > threshold) &&
      all(cohort_values[["4N"]] > threshold)
  )
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

get_unfiltered_forest_top_seeds <- function(summary_df, is_joint_run = FALSE, is_invitro_run = FALSE, n = 3L) {
  rank_col <- if ("objective_rank" %in% names(summary_df)) "objective_rank" else "objective"
  get_top_ranked_seeds(summary_df, n = n, rank_col = rank_col)
}

WARM_START_RATIO_PARAMETERS <- c(
  "O2_crit",
  "alpha_o2",
  "mu_hp",
  "p_misseg",
  "k_o_mis",
  "buffer_smax",
  "buffer_beta",
  "buffer_n_exp",
  "n_O",
  "gamma_growth",
  "lam_max",
  "p_mis_base",
  "p_wgd",
  "gamma_mu"
)

truthy_vector <- function(x) {
  if (is.null(x)) return(logical(0))
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y", "on")
}

has_nonmissing_text <- function(x) {
  if (is.null(x)) return(logical(0))
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & !tolower(x) %in% c("na", "nan", "null")
}

deoptim_converged_flags <- function(seed_summary) {
  n <- nrow(seed_summary)
  if (!n) return(logical(0))
  de_obj <- if ("optimizer_deoptim_objective" %in% names(seed_summary)) {
    suppressWarnings(as.numeric(seed_summary$optimizer_deoptim_objective))
  } else {
    rep(NA_real_, n)
  }
  if ("objective" %in% names(seed_summary)) {
    objective <- suppressWarnings(as.numeric(seed_summary$objective))
    fallback_idx <- !is.finite(de_obj) & is.finite(objective)
    de_obj[fallback_idx] <- objective[fallback_idx]
  }
  has_objective <- is.finite(de_obj)
  converged <- rep(FALSE, n)
  has_convergence_evidence <- rep(FALSE, n)

  if ("optimizer_interrupted" %in% names(seed_summary)) {
    interrupted <- truthy_vector(seed_summary$optimizer_interrupted)
  } else {
    interrupted <- rep(FALSE, n)
  }
  if (length(interrupted) != n) interrupted <- rep(FALSE, n)

  if ("deoptim_stop_reason" %in% names(seed_summary)) {
    reason <- trimws(as.character(seed_summary$deoptim_stop_reason))
    has_reason <- has_nonmissing_text(reason)
    reason_lower <- tolower(reason)
    converged[has_reason] <- reason_lower[has_reason] %in% c("early_stop_reltol_or_steptol")
    has_convergence_evidence[has_reason] <- TRUE
  }

  iter_completed <- if ("optimizer_iter_completed" %in% names(seed_summary)) {
    suppressWarnings(as.numeric(seed_summary$optimizer_iter_completed))
  } else {
    rep(NA_real_, n)
  }
  iter_target <- if ("optimizer_iter_target" %in% names(seed_summary)) {
    suppressWarnings(as.numeric(seed_summary$optimizer_iter_target))
  } else {
    rep(NA_real_, n)
  }
  if ("itermax" %in% names(seed_summary)) {
    itermax <- suppressWarnings(as.numeric(seed_summary$itermax))
    fill_itermax <- !is.finite(iter_target) & is.finite(itermax)
    iter_target[fill_itermax] <- itermax[fill_itermax]
  }
  legacy_target_idx <- !is.finite(iter_target) & is.finite(iter_completed)
  iter_target[legacy_target_idx] <- 500

  has_iter_evidence <- is.finite(iter_completed) & is.finite(iter_target)
  fill_from_iter <- has_iter_evidence & !has_convergence_evidence
  converged[fill_from_iter] <- iter_completed[fill_from_iter] < iter_target[fill_from_iter]
  has_convergence_evidence[fill_from_iter] <- TRUE
  max_iter_reached <- has_iter_evidence & iter_completed >= iter_target
  converged[max_iter_reached] <- FALSE
  has_convergence_evidence[max_iter_reached] <- TRUE

  converged <- converged & has_convergence_evidence & has_objective & !interrupted
  converged[is.na(converged)] <- FALSE
  converged
}

local_refinement_accepted_flags <- function(seed_summary) {
  n <- nrow(seed_summary)
  if (!n) return(logical(0))
  if (!("optimizer_local_accepted" %in% names(seed_summary))) return(rep(FALSE, n))
  accepted <- truthy_vector(seed_summary$optimizer_local_accepted)
  if (length(accepted) != n) rep(FALSE, n) else accepted
}

is_missing_summary_value <- function(x) {
  if (!length(x)) return(logical(0))
  text <- trimws(as.character(x))
  is.na(x) | !nzchar(text) | tolower(text) %in% c("na", "nan", "null")
}

supplement_optimizer_fields_from_refinement_csv <- function(seed_summary, run_dir) {
  if (is.null(seed_summary) || !is.data.frame(seed_summary) || !nrow(seed_summary)) return(seed_summary)
  if (!("seed" %in% names(seed_summary))) return(seed_summary)
  refinement_path <- file.path(run_dir, "lbfgsb_refinement_accepted_seeds.csv")
  if (!file.exists(refinement_path)) return(seed_summary)
  refinement <- tryCatch(
    utils::read.csv(refinement_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(refinement) || !is.data.frame(refinement) || !nrow(refinement) || !("seed" %in% names(refinement))) {
    return(seed_summary)
  }
  optimizer_cols <- intersect(
    c(
      "optimizer_deoptim_objective",
      "optimizer_local_objective",
      "optimizer_local_attempted",
      "optimizer_local_accepted",
      "optimizer_local_convergence",
      "optimizer_local_maxit"
    ),
    names(refinement)
  )
  if (!length(optimizer_cols)) return(seed_summary)
  match_idx <- match(as.character(seed_summary$seed), as.character(refinement$seed))
  has_match <- !is.na(match_idx)
  if (!any(has_match)) return(seed_summary)
  for (col in optimizer_cols) {
    if (!(col %in% names(seed_summary))) seed_summary[[col]] <- NA
    replacement <- rep(NA, nrow(seed_summary))
    replacement[has_match] <- refinement[[col]][match_idx[has_match]]
    fill_idx <- has_match & is_missing_summary_value(seed_summary[[col]]) & !is_missing_summary_value(replacement)
    if (any(fill_idx)) seed_summary[[col]][fill_idx] <- replacement[fill_idx]
  }
  seed_summary
}

fit_label_for_summary <- function(seed_summary) {
  fit_mode <- if ("fit_mode" %in% names(seed_summary)) unique(as.character(seed_summary$fit_mode)) else character(0)
  fit_mode <- fit_mode[has_nonmissing_text(fit_mode)]
  if (any(fit_mode == "fit_joint")) return("joint fitting")
  if (any(fit_mode == "fit_invitro")) return("in vitro")
  "in vivo"
}

build_convergence_summary <- function(seed_summary) {
  de_converged <- deoptim_converged_flags(seed_summary)
  fit_label <- fit_label_for_summary(seed_summary)
  if (!identical(fit_label, "joint fitting")) {
    return(data.frame(
      Fit = fit_label,
      `Total seeds` = nrow(seed_summary),
      `DEoptim converged` = sum(de_converged, na.rm = TRUE),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  local_accepted <- local_refinement_accepted_flags(seed_summary)
  both <- de_converged & local_accepted
  data.frame(
    Fit = fit_label,
    `Total seeds` = nrow(seed_summary),
    `DEoptim converged` = sum(de_converged, na.rm = TRUE),
    `L-BFGS-B accepted` = sum(local_accepted, na.rm = TRUE),
    `Converged and accepted` = sum(both, na.rm = TRUE),
    `Converged only` = sum(de_converged & !local_accepted, na.rm = TRUE),
    `Accepted only` = sum(!de_converged & local_accepted, na.rm = TRUE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

build_convergence_venn_counts <- function(seed_summary) {
  de_converged <- deoptim_converged_flags(seed_summary)
  local_accepted <- local_refinement_accepted_flags(seed_summary)
  data.frame(
    Fit = fit_label_for_summary(seed_summary),
    `DEoptim only` = sum(de_converged & !local_accepted, na.rm = TRUE),
    `L-BFGS-B only` = sum(!de_converged & local_accepted, na.rm = TRUE),
    Both = sum(de_converged & local_accepted, na.rm = TRUE),
    Neither = sum(!de_converged & !local_accepted, na.rm = TRUE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

plot_convergence_venn <- function(seed_summary,
                                  out_dir,
                                  run_label,
                                  filename_prefix = "convergence",
                                  width = 7,
                                  height = 7) {
  if (is.null(seed_summary) || !is.data.frame(seed_summary) || !nrow(seed_summary)) return(invisible(NULL))
  counts <- build_convergence_venn_counts(seed_summary)
  counts_path <- file.path(out_dir, paste0(filename_prefix, "_venn_counts.tsv"))
  utils::write.table(counts, file = counts_path, sep = "\t", quote = FALSE, row.names = FALSE)

  theta <- seq(0, 2 * pi, length.out = 361L)
  circle_df <- function(cx, cy, r, set_name) {
    data.frame(
      x = cx + r * cos(theta),
      y = cy + r * sin(theta),
      set = set_name,
      stringsAsFactors = FALSE
    )
  }
  circles <- rbind(
    circle_df(-0.55, 0, 1, "DEoptim converged"),
    circle_df(0.55, 0, 1, "LBFGSB accepted")
  )
  de_only <- counts[["DEoptim only"]][[1]]
  accepted_only <- counts[["L-BFGS-B only"]][[1]]
  both <- counts[["Both"]][[1]]
  neither <- counts[["Neither"]][[1]]
  total <- nrow(seed_summary)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = circles,
      ggplot2::aes(x = x, y = y, group = set, fill = set),
      alpha = 0.26,
      color = NA
    ) +
    ggplot2::geom_path(
      data = circles,
      ggplot2::aes(x = x, y = y, group = set, color = set),
      linewidth = 1.1
    ) +
    ggplot2::annotate("text", x = -1.08, y = 0.03, label = de_only, size = 7, fontface = "bold", color = "#1f4e79") +
    ggplot2::annotate("text", x = 0, y = 0.03, label = both, size = 7, fontface = "bold", color = "#384860") +
    ggplot2::annotate("text", x = 1.08, y = 0.03, label = accepted_only, size = 7, fontface = "bold", color = "#8f3f4d") +
    ggplot2::annotate("text", x = -0.55, y = 1.17, label = "DEoptim converged", size = 3.7, color = "#1f4e79") +
    ggplot2::annotate("text", x = 0.55, y = 1.17, label = "LBFGSB accepted", size = 3.7, color = "#8f3f4d") +
    ggplot2::annotate("text", x = 0, y = -1.32, label = paste0("Neither: ", neither, "   Total seeds: ", total), size = 3.5, color = "grey25") +
    ggplot2::scale_fill_manual(values = c("DEoptim converged" = "#2b6cb0", "LBFGSB accepted" = "#d33f3f")) +
    ggplot2::scale_color_manual(values = c("DEoptim converged" = "#2b6cb0", "LBFGSB accepted" = "#d33f3f")) +
    ggplot2::coord_fixed(xlim = c(-1.75, 1.75), ylim = c(-1.48, 1.36), clip = "off") +
    ggplot2::labs(
      title = paste0(fit_label_for_summary(seed_summary), " Convergence Venn"),
      subtitle = "",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "plain", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5),
      legend.position = "none",
      plot.margin = ggplot2::margin(20, 20, 20, 20)
    )

  pdf_path <- file.path(out_dir, paste0(filename_prefix, "_venn.pdf"))
  png_path <- file.path(out_dir, paste0(filename_prefix, "_venn.png"))
  ggplot2::ggsave(pdf_path, p, width = width, height = height, bg = "white")
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 220, bg = "white")
  invisible(pdf_path)
}

square_umap_limits <- function(x, y, pad_fraction = 0.12) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) {
    return(list(x = c(-1, 1), y = c(-1, 1)))
  }
  x_mid <- mean(range(x))
  y_mid <- mean(range(y))
  span <- max(diff(range(x)), diff(range(y)), 1e-6)
  half_span <- span * (0.5 + pad_fraction)
  list(
    x = c(x_mid - half_span, x_mid + half_span),
    y = c(y_mid - half_span, y_mid + half_span)
  )
}

plot_top_seed_parameter_umap <- function(summary_df,
                                         parameter_long,
                                         out_dir,
                                         run_label,
                                         fit_label,
                                         filename_prefix,
                                         top_n = 20L,
                                         umap_seed = 1L) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    warning("Skipping ", fit_label, " top", top_n, " parameter UMAP because package 'uwot' is not available.", call. = FALSE)
    return(invisible(NULL))
  }
  if (is.null(summary_df) || !nrow(summary_df) || is.null(parameter_long) || !nrow(parameter_long)) {
    return(invisible(NULL))
  }
  if (!all(c("seed", "objective") %in% names(summary_df))) return(invisible(NULL))
  if (!all(c("seed", "transformed_value") %in% names(parameter_long))) return(invisible(NULL))

  summary_df$objective <- suppressWarnings(as.numeric(summary_df$objective))
  candidates <- summary_df[is.finite(summary_df$objective), , drop = FALSE]
  if (!nrow(candidates)) return(invisible(NULL))
  candidates <- candidates[order(candidates$objective, seed_order_key(candidates$seed), candidates$seed, na.last = TRUE), , drop = FALSE]
  candidates <- utils::head(candidates, as.integer(top_n))
  candidates$umap_rank <- seq_len(nrow(candidates))
  seed_levels <- as.character(candidates$seed)

  parameter_long$seed <- as.character(parameter_long$seed)
  param_id_col <- if ("param_name" %in% names(parameter_long)) "param_name" else if ("param_prototype" %in% names(parameter_long)) "param_prototype" else NA_character_
  if (is.na(param_id_col)) return(invisible(NULL))

  umap_long <- parameter_long[
    parameter_long$seed %in% seed_levels &
      has_nonmissing_text(parameter_long[[param_id_col]]),
    ,
    drop = FALSE
  ]
  if (!nrow(umap_long)) return(invisible(NULL))
  umap_long$transformed_value <- suppressWarnings(as.numeric(umap_long$transformed_value))
  umap_long <- umap_long[is.finite(umap_long$transformed_value), , drop = FALSE]
  if (!nrow(umap_long)) return(invisible(NULL))

  params <- unique(as.character(umap_long[[param_id_col]]))
  x <- matrix(
    NA_real_,
    nrow = length(seed_levels),
    ncol = length(params),
    dimnames = list(seed_levels, params)
  )
  for (i in seq_len(nrow(umap_long))) {
    x[as.character(umap_long$seed[[i]]), as.character(umap_long[[param_id_col]][[i]])] <- umap_long$transformed_value[[i]]
  }
  complete_cols <- colSums(is.finite(x)) == nrow(x)
  x <- x[, complete_cols, drop = FALSE]
  if (ncol(x) < 2L || nrow(x) < 3L) return(invisible(NULL))
  varying_cols <- apply(x, 2L, function(v) stats::sd(v, na.rm = TRUE) > 0)
  x <- x[, varying_cols, drop = FALSE]
  if (ncol(x) < 2L) return(invisible(NULL))

  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  set.seed(as.integer(umap_seed))
  embedding <- uwot::umap(
    x_scaled,
    n_neighbors = min(15L, nrow(x_scaled) - 1L),
    min_dist = 0.1,
    metric = "euclidean",
    n_components = 2L,
    n_threads = 1L,
    ret_model = FALSE,
    verbose = FALSE
  )

  coords <- data.frame(
    seed = seed_levels,
    objective = candidates$objective,
    objective_rank = candidates$umap_rank,
    UMAP1 = embedding[, 1],
    UMAP2 = embedding[, 2],
    stringsAsFactors = FALSE
  )
  umap_limits <- square_umap_limits(coords$UMAP1, coords$UMAP2)
  coords_path <- file.path(out_dir, paste0(filename_prefix, "_top20_seed_parameter_umap_coords.tsv"))
  matrix_path <- file.path(out_dir, paste0(filename_prefix, "_top20_seed_parameter_umap_matrix.tsv"))
  utils::write.table(coords, file = coords_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(
    data.frame(seed = rownames(x), x, check.names = FALSE),
    file = matrix_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p <- ggplot2::ggplot(coords, ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = objective_rank), size = 4.2, alpha = 0.95) +
    ggplot2::scale_color_gradient(
      low = "#2b6cb0",
      high = "#d33f3f",
      breaks = unique(round(seq(1, nrow(coords), length.out = min(5L, nrow(coords))), 0)),
      name = "Objective Rank"
    ) +
    ggplot2::coord_fixed(xlim = umap_limits$x, ylim = umap_limits$y, expand = FALSE, clip = "off") +
    ggplot2::labs(
      title = paste0(fit_label, " Top 20 Seed Parameter UMAP"),
      subtitle = paste0(nrow(coords), " seeds; ", ncol(x), " transformed parameters; ", run_label),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "plain"),
      plot.subtitle = ggplot2::element_text(size = 9),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = seed),
      size = 3.0,
      color = "grey15",
      min.segment.length = 0,
      segment.alpha = 0.35,
      box.padding = 0.25,
      point.padding = 0.15,
      seed = as.integer(umap_seed),
      max.overlaps = Inf
    )
  } else {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = seed), size = 3.0, vjust = -0.8, color = "grey15")
  }

  pdf_path <- file.path(out_dir, paste0(filename_prefix, "_top20_seed_parameter_umap.pdf"))
  png_path <- file.path(out_dir, paste0(filename_prefix, "_top20_seed_parameter_umap.png"))
  ggplot2::ggsave(pdf_path, p, width = 7, height = 7, bg = "white")
  ggplot2::ggsave(png_path, p, width = 7, height = 7, dpi = 220, bg = "white")
  invisible(pdf_path)
}

seed_label_to_integer <- function(x) {
  x <- as.character(x)
  out <- seed_order_key(x)
  missing <- !is.finite(out)
  out[missing] <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", x[missing])))
  as.integer(out)
}

read_fit_summary_values_optional <- function(seed_dir) {
  read_metric_map_optional(file.path(seed_dir, "fit_summary.tsv"), "metric", "value")
}

localize_project_path <- function(path, project_root) {
  if (is.null(path) || !length(path)) return(NA_character_)
  path <- trimws(as.character(path[[1]]))
  if (!nzchar(path) || is.na(path)) return(NA_character_)
  if (file.exists(path)) return(normalizePath(path, mustWork = TRUE))
  if (!startsWith(path, "/") && file.exists(file.path(project_root, path))) {
    return(normalizePath(file.path(project_root, path), mustWork = TRUE))
  }
  marker <- "oxygen/results/"
  pos <- regexpr(marker, path, fixed = TRUE)
  if (pos > 0L) {
    suffix <- substr(path, pos, nchar(path))
    candidate <- file.path(project_root, suffix)
    if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
    return(normalizePath(candidate, mustWork = FALSE))
  }
  normalizePath(path, mustWork = FALSE)
}

joint_warmup_source_paths <- function(run_dir) {
  seed_dirs <- find_seed_dirs(run_dir)
  if (!length(seed_dirs)) seed_dirs <- run_dir
  for (seed_dir in seed_dirs) {
    vals <- read_fit_summary_values_optional(seed_dir)
    invivo <- summary_metric_value(vals, "joint_warmup_invivo_seed_dir", NA_character_)
    invitro <- summary_metric_value(vals, "joint_warmup_invitro_seed_dir", NA_character_)
    if (has_nonmissing_text(invivo) && has_nonmissing_text(invitro)) {
      return(list(invivo = as.character(invivo), invitro = as.character(invitro)))
    }
  }
  list(invivo = NA_character_, invitro = NA_character_)
}

resolve_warmup_source_run_dir <- function(warmup_seed_dir, project_root) {
  if (!has_nonmissing_text(warmup_seed_dir)) return(NA_character_)
  local_seed_dir <- localize_project_path(warmup_seed_dir, project_root)
  source_run_dir <- dirname(local_seed_dir)
  if (dir.exists(source_run_dir)) {
    normalizePath(source_run_dir, mustWork = TRUE)
  } else {
    normalizePath(source_run_dir, mustWork = FALSE)
  }
}

read_seed_params_for_ratio <- function(seed_dir, params) {
  param_path <- file.path(seed_dir, "best_params.tsv")
  if (file.exists(param_path)) {
    tab <- tryCatch(
      utils::read.delim(param_path, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(tab) && all(c("parameter", "value") %in% names(tab))) {
      vals <- setNames(suppressWarnings(as.numeric(tab$value)), as.character(tab$parameter))
      out <- as.numeric(vals[params])
      if (all(is.finite(out))) return(out)
    }
  }
  rep(NA_real_, length(params))
}

read_warmup_source_seed_summary <- function(source_run_dir, kind) {
  if (!has_nonmissing_text(source_run_dir) || !dir.exists(source_run_dir)) {
    warning("Skipping joint ratio UMAP source ", kind, " because warmup source run directory was not found: ", source_run_dir, call. = FALSE)
    return(NULL)
  }
  seed_dirs <- find_seed_dirs(source_run_dir)
  if (!length(seed_dirs)) {
    warning("Skipping joint ratio UMAP source ", kind, " because no valid seed directories were found under: ", source_run_dir, call. = FALSE)
    return(NULL)
  }
  rows <- vector("list", length(seed_dirs))
  for (i in seq_along(seed_dirs)) {
    seed_dir <- seed_dirs[[i]]
    fit_summary_vals <- read_fit_summary_values_optional(seed_dir)
    objective <- as_num(
      summary_metric_value(fit_summary_vals, "objective", summary_metric_value(fit_summary_vals, "objective_total", NA_real_)),
      NA_real_
    )
    rows[[i]] <- data.frame(
      seed = basename(seed_dir),
      seed_dir = normalizePath(seed_dir, mustWork = TRUE),
      objective = objective,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out <- out[is.finite(out$objective), , drop = FALSE]
  if (!nrow(out)) {
    warning("Skipping joint ratio UMAP source ", kind, " because no finite objectives were found under: ", source_run_dir, call. = FALSE)
    return(NULL)
  }
  out[order(out$objective, seed_order_key(out$seed), out$seed, na.last = TRUE), , drop = FALSE]
}

read_ratio_seed_ranking_from_source_run <- function(source_run_dir, kind, top_n = 10L) {
  summary_df <- read_warmup_source_seed_summary(source_run_dir, kind)
  if (is.null(summary_df) || !is.data.frame(summary_df) || !nrow(summary_df)) return(NULL)
  top_n <- as.integer(top_n)
  if (!is.finite(top_n) || top_n < 1L) top_n <- 10L
  if (nrow(summary_df) < top_n) {
    warning(
      "Skipping joint ratio UMAP source ", kind,
      " because warmup source run has only ", nrow(summary_df),
      " valid seeds; need ", top_n, " to build source top", top_n, ".",
      call. = FALSE
    )
    return(NULL)
  }
  top <- utils::head(summary_df, as.integer(top_n))
  rows <- vector("list", nrow(top))
  for (i in seq_len(nrow(top))) {
    seed_label <- as.character(top$seed[[i]])
    param_vals <- read_seed_params_for_ratio(top$seed_dir[[i]], WARM_START_RATIO_PARAMETERS)
    if (any(!is.finite(param_vals))) {
      missing <- WARM_START_RATIO_PARAMETERS[!is.finite(param_vals)]
      warning(
        "Skipping joint ratio UMAP source ", kind, " seed ", seed_label,
        " because parameter values are missing or non-finite: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
      return(NULL)
    }
    row <- data.frame(
      rank = i,
      seed = seed_label_to_integer(seed_label),
      objective = as.numeric(top$objective[[i]]),
      stringsAsFactors = FALSE
    )
    for (p_name in WARM_START_RATIO_PARAMETERS) {
      row[[p_name]] <- param_vals[[match(p_name, WARM_START_RATIO_PARAMETERS)]]
    }
    rows[[i]] <- row
  }
  do.call(rbind, rows)
}

read_joint_warmup_ratio_rankings <- function(run_dir, top_n = 10L) {
  project_root <- normalizePath(file.path(WORKFLOW_ROOT, "../../.."), mustWork = FALSE)
  warmup_paths <- joint_warmup_source_paths(run_dir)
  invivo_source <- resolve_warmup_source_run_dir(warmup_paths$invivo, project_root)
  invitro_source <- resolve_warmup_source_run_dir(warmup_paths$invitro, project_root)
  invivo <- read_ratio_seed_ranking_from_source_run(invivo_source, "invivo", top_n = top_n)
  invitro <- read_ratio_seed_ranking_from_source_run(invitro_source, "invitro", top_n = top_n)
  if (is.data.frame(invivo) && nrow(invivo) && is.data.frame(invitro) && nrow(invitro)) {
    return(list(
      invivo = invivo,
      invitro = invitro,
      sources = data.frame(
        kind = c("invivo", "invitro"),
        source_dir = c(invivo_source, invitro_source),
        warmup_seed_dir = c(warmup_paths$invivo, warmup_paths$invitro),
        source_type = "joint_warmup_source_run",
        stringsAsFactors = FALSE
      )
    ))
  }
  NULL
}

joint_ratio_umap_output_paths <- function(out_dir) {
  file.path(
    out_dir,
    c(
      "cross_paired_top10_ratio_matrix.tsv",
      "cross_paired_top10_umap_coords.tsv",
      "cross_paired_top10_ratio_sources.tsv",
      "joint_soft_coupling_ratio_umap_500seed.png",
      "joint_soft_coupling_ratio_umap_500seed.pdf"
    )
  )
}

clear_joint_ratio_umap_outputs <- function(out_dir) {
  paths <- joint_ratio_umap_output_paths(out_dir)
  unlink(paths[file.exists(paths)], force = TRUE)
}

plot_joint_ratio_umap <- function(run_dir,
                                  out_dir,
                                  run_label,
                                  top_n = 10L,
                                  umap_seed = 1L) {
  clear_joint_ratio_umap_outputs(out_dir)
  if (!requireNamespace("uwot", quietly = TRUE)) {
    warning("Skipping joint ratio UMAP because package 'uwot' is not available.", call. = FALSE)
    return(invisible(NULL))
  }
  rankings <- read_joint_warmup_ratio_rankings(run_dir, top_n = top_n)
  if (is.null(rankings) || !is.data.frame(rankings$invivo) || !is.data.frame(rankings$invitro)) {
    warning("Skipping joint ratio UMAP because joint warm-start source top-seed rankings were not found.", call. = FALSE)
    return(invisible(NULL))
  }
  invivo <- rankings$invivo
  invitro <- rankings$invitro
  missing_vivo <- setdiff(WARM_START_RATIO_PARAMETERS, names(invivo))
  missing_vitro <- setdiff(WARM_START_RATIO_PARAMETERS, names(invitro))
  if (length(missing_vivo) || length(missing_vitro)) {
    warning(
      "Skipping joint ratio UMAP because warmup source rankings are missing parameter columns: ",
      paste(c(missing_vivo, missing_vitro), collapse = ", "),
      call. = FALSE
    )
    return(invisible(NULL))
  }
  if (!all(c("rank", "seed", "objective") %in% names(invivo)) || !all(c("rank", "seed", "objective") %in% names(invitro))) {
    warning("Skipping joint ratio UMAP because warmup source rankings are missing rank, seed, or objective columns.", call. = FALSE)
    return(invisible(NULL))
  }
  invivo$rank <- suppressWarnings(as.integer(invivo$rank))
  invitro$rank <- suppressWarnings(as.integer(invitro$rank))
  invivo_top <- invivo[order(invivo$rank), , drop = FALSE]
  invitro_top <- invitro[order(invitro$rank), , drop = FALSE]
  invivo_top <- utils::head(invivo_top[is.finite(invivo_top$rank), , drop = FALSE], as.integer(top_n))
  invitro_top <- utils::head(invitro_top[is.finite(invitro_top$rank), , drop = FALSE], as.integer(top_n))
  if (!nrow(invivo_top) || !nrow(invitro_top)) return(invisible(NULL))

  pair_rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(invivo_top))) {
    for (j in seq_len(nrow(invitro_top))) {
      vivo_vals <- suppressWarnings(as.numeric(invivo_top[i, WARM_START_RATIO_PARAMETERS]))
      vitro_vals <- suppressWarnings(as.numeric(invitro_top[j, WARM_START_RATIO_PARAMETERS]))
      ratio_vals <- vivo_vals / vitro_vals
      log_ratio_vals <- log10(ratio_vals)
      if (any(!is.finite(log_ratio_vals))) next
      row <- data.frame(
        pair_id = sprintf("V%02d-I%02d", invivo_top$rank[[i]], invitro_top$rank[[j]]),
        invivo_rank = as.integer(invivo_top$rank[[i]]),
        invitro_rank = as.integer(invitro_top$rank[[j]]),
        invivo_seed = as.integer(invivo_top$seed[[i]]),
        invitro_seed = as.integer(invitro_top$seed[[j]]),
        invivo_objective = as.numeric(invivo_top$objective[[i]]),
        invitro_objective = as.numeric(invitro_top$objective[[j]]),
        stringsAsFactors = FALSE
      )
      for (p_name in WARM_START_RATIO_PARAMETERS) row[[paste0("ratio_", p_name)]] <- ratio_vals[[match(p_name, WARM_START_RATIO_PARAMETERS)]]
      for (p_name in WARM_START_RATIO_PARAMETERS) row[[paste0("log10_ratio_", p_name)]] <- log_ratio_vals[[match(p_name, WARM_START_RATIO_PARAMETERS)]]
      pair_rows[[idx]] <- row
      idx <- idx + 1L
    }
  }
  if (!length(pair_rows)) return(invisible(NULL))
  pairs <- do.call(rbind, pair_rows)
  matrix_cols <- paste0("log10_ratio_", WARM_START_RATIO_PARAMETERS)
  x <- as.matrix(pairs[, matrix_cols, drop = FALSE])
  colnames(x) <- WARM_START_RATIO_PARAMETERS
  if (nrow(x) < 3L) {
    warning("Skipping joint ratio UMAP because fewer than 3 warmup source cross-pairs were available.", call. = FALSE)
    return(invisible(NULL))
  }

  set.seed(as.integer(umap_seed))
  embedding <- uwot::umap(
    x,
    n_neighbors = min(15L, nrow(x) - 1L),
    min_dist = 0.1,
    metric = "euclidean",
    n_components = 2L,
    n_threads = 1L,
    ret_model = FALSE,
    verbose = FALSE
  )
  coords <- cbind(
    pairs[, c("pair_id", "invivo_rank", "invitro_rank", "invivo_seed", "invitro_seed", "invivo_objective", "invitro_objective")],
    data.frame(UMAP1 = embedding[, 1], UMAP2 = embedding[, 2], stringsAsFactors = FALSE)
  )
  umap_limits <- square_umap_limits(coords$UMAP1, coords$UMAP2)
  plot_df <- coords
  plot_df$invitro_rank_factor <- factor(plot_df$invitro_rank, levels = seq_len(top_n))
  shape_values <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4)
  names(shape_values) <- as.character(seq_along(shape_values))

  subtitle <- paste0(nrow(pairs), " joint warm-start source top10 x top10 pairings; 14 dimensions from log10(in vivo / in vitro parameter ratio)")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = invivo_rank, shape = invitro_rank_factor),
      size = 2.8,
      stroke = 0.9
    ) +
    ggplot2::scale_color_gradient(
      low = "#2b6cb0",
      high = "#d33f3f",
      limits = c(1, top_n),
      breaks = unique(round(seq(1, top_n, length.out = 5), 1)),
      name = "In Vivo Rank"
    ) +
    ggplot2::scale_shape_manual(values = shape_values[as.character(seq_len(top_n))], name = "In Vitro Rank") +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(order = 1),
      shape = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::coord_fixed(xlim = umap_limits$x, ylim = umap_limits$y, expand = FALSE, clip = "off") +
    ggplot2::labs(
      title = "Joint Soft Coupling Ratio UMAP",
      subtitle = subtitle,
      x = "UMAP1",
      y = "UMAP2"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "plain"),
      plot.subtitle = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = pair_id),
      color = "grey15",
      size = 2.1,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.alpha = 0.35,
      box.padding = 0.2,
      point.padding = 0.1,
      seed = as.integer(umap_seed)
    )
  } else {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = pair_id), size = 2.1, vjust = -0.6, color = "grey15")
  }

  ratio_path <- file.path(out_dir, "cross_paired_top10_ratio_matrix.tsv")
  coords_path <- file.path(out_dir, "cross_paired_top10_umap_coords.tsv")
  sources_path <- file.path(out_dir, "cross_paired_top10_ratio_sources.tsv")
  png_path <- file.path(out_dir, "joint_soft_coupling_ratio_umap_500seed.png")
  pdf_path <- file.path(out_dir, "joint_soft_coupling_ratio_umap_500seed.pdf")
  utils::write.table(pairs, ratio_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(coords, coords_path, sep = "\t", quote = FALSE, row.names = FALSE)
  if (is.data.frame(rankings$sources)) {
    utils::write.table(rankings$sources, sources_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  ggplot2::ggsave(pdf_path, p, width = 7, height = 7, bg = "white")
  ggplot2::ggsave(png_path, p, width = 7, height = 7, dpi = 220, bg = "white")
  invisible(pdf_path)
}

boundary_axis_config <- function(x,
                                 near_thresh,
                                 x_scale = c("relative", "log10_original"),
                                 raw_value = NULL,
                                 raw_lower = NULL,
                                 raw_upper = NULL) {
  x_scale <- match.arg(x_scale)
  if (identical(x_scale, "relative")) {
    return(list(
      axis_type = "relative",
      x = x,
      lower_limit = 0,
      upper_limit = 1,
      vlines = c(0, near_thresh, 0.5, 1 - near_thresh, 1),
      vline_colors = c("grey50", "grey70", "grey82", "grey70", "grey50"),
      vline_linetypes = c("solid", "dashed", "dotted", "dashed", "solid"),
      lower_rect = c(0, near_thresh),
      upper_rect = c(1 - near_thresh, 1),
      scale = ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, near_thresh, 0.5, 1 - near_thresh, 1)),
      x_label = "Relative position in transformed fit range",
      subtitle_note = ""
    ))
  }

  raw_value <- suppressWarnings(as.numeric(raw_value))
  raw_lower <- suppressWarnings(as.numeric(raw_lower))
  raw_upper <- suppressWarnings(as.numeric(raw_upper))
  positive_x <- c(raw_value, raw_lower, raw_upper)
  positive_x <- positive_x[is.finite(positive_x) & positive_x > 0]
  log_floor <- min(positive_x, na.rm = TRUE) / 10
  if (!is.finite(log_floor) || log_floor <= 0) log_floor <- 1e-6
  lower_plot <- ifelse(is.finite(raw_lower), pmax(raw_lower, log_floor), NA_real_)
  upper_plot <- ifelse(is.finite(raw_upper), pmax(raw_upper, log_floor), NA_real_)
  value_plot <- ifelse(is.finite(raw_value), pmax(raw_value, log_floor), NA_real_)
  near_lower_raw <- raw_lower + near_thresh * (raw_upper - raw_lower)
  near_upper_raw <- raw_upper - near_thresh * (raw_upper - raw_lower)
  lower_near_plot <- ifelse(is.finite(near_lower_raw), pmax(near_lower_raw, log_floor), NA_real_)
  upper_near_plot <- ifelse(is.finite(near_upper_raw), pmax(near_upper_raw, log_floor), NA_real_)
  axis_upper <- max(c(upper_plot, value_plot), na.rm = TRUE)
  if (!is.finite(axis_upper) || axis_upper <= log_floor) axis_upper <- log_floor * 10
  breaks <- 10^seq(floor(log10(log_floor)), ceiling(log10(axis_upper)), by = 1)
  breaks <- breaks[breaks >= log_floor & breaks <= axis_upper]
  needs_floor_label <- any(
    is.finite(c(raw_value, raw_lower, raw_upper)) & c(raw_value, raw_lower, raw_upper) <= 0,
    na.rm = TRUE
  )
  labels <- function(vals) {
    vapply(
      vals,
      function(v) {
        if (isTRUE(needs_floor_label) && isTRUE(all.equal(v, log_floor, tolerance = 1e-12))) return("0")
        formatC(v, format = "fg", digits = 4)
      },
      character(1)
    )
  }
  list(
    axis_type = "log10_original",
    x = value_plot,
    lower_limit = log_floor,
    upper_limit = axis_upper,
    lower_plot = lower_plot,
    upper_plot = upper_plot,
    lower_near_plot = lower_near_plot,
    upper_near_plot = upper_near_plot,
    scale = ggplot2::scale_x_log10(limits = c(log_floor, axis_upper), breaks = breaks, labels = labels),
    x_label = "Original parameter value (log10 scale)",
    subtitle_note = if (isTRUE(needs_floor_label)) " Non-positive raw values or bounds are shown at the log floor labeled 0." else ""
  )
}

plot_parameter_boundary_forest <- function(long_df,
                                           summary_df,
                                           out_path,
                                           run_label,
                                           near_thresh = 0.05,
                                           top3_seeds = NULL,
                                           title_suffix = NULL,
                                           legend_title = "Objective Top 3 Seeds",
                                           x_scale = c("relative", "log10_original")) {
  x_scale <- match.arg(x_scale)
  plot_df <- long_df[long_df$active_in_fit & is.finite(long_df$rel_pos_plot), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))

  if (is.null(top3_seeds)) top3_seeds <- get_top_ranked_seeds(summary_df, n = 3L)
  top3_seeds <- as.character(top3_seeds)

  axis_cfg <- boundary_axis_config(
    plot_df$rel_pos_plot,
    near_thresh = near_thresh,
    x_scale = x_scale,
    raw_value = plot_df$prototype_value,
    raw_lower = plot_df$prototype_lower_bound,
    raw_upper = plot_df$prototype_upper_bound
  )
  plot_df$boundary_x_plot <- axis_cfg$x
  if (identical(axis_cfg$axis_type, "log10_original")) {
    plot_df$boundary_x_lower <- axis_cfg$lower_plot
    plot_df$boundary_x_upper <- axis_cfg$upper_plot
    plot_df$boundary_x_lower_near <- axis_cfg$lower_near_plot
    plot_df$boundary_x_upper_near <- axis_cfg$upper_near_plot
    plot_df <- plot_df[
      is.finite(plot_df$boundary_x_plot) &
        is.finite(plot_df$boundary_x_lower) &
        is.finite(plot_df$boundary_x_upper) &
        plot_df$boundary_x_upper > plot_df$boundary_x_lower,
      ,
      drop = FALSE
    ]
    if (!nrow(plot_df)) return(invisible(NULL))
  }
  param_rank <- tapply(plot_df$rel_dist_to_nearest, plot_df$param_prototype, min, na.rm = TRUE)
  param_levels <- names(sort(param_rank, decreasing = FALSE))
  plot_df$param_prototype <- factor(plot_df$param_prototype, levels = rev(param_levels))
  ref_df <- if (identical(axis_cfg$axis_type, "log10_original")) {
    unique(plot_df[c("param_prototype", "boundary_x_lower", "boundary_x_upper")])
  } else {
    ref <- unique(plot_df["param_prototype"])
    ref$boundary_x_start <- axis_cfg$lower_limit
    ref$boundary_x_end <- axis_cfg$upper_limit
    ref
  }
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
  top3_label <- if (length(top3_seeds)) {
    paste(paste0("Top ", seq_along(top3_seeds), ": ", top3_seeds), collapse = "; ")
  } else {
    "No seeds met the 2N/4N 1000-day prediction gate."
  }

  p <- ggplot(plot_df, aes(x = boundary_x_plot, y = param_prototype))
  if (identical(axis_cfg$axis_type, "relative")) {
    p <- p +
      annotate("rect", xmin = axis_cfg$lower_rect[[1]], xmax = axis_cfg$lower_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#fddbc7", alpha = 0.28) +
      annotate("rect", xmin = axis_cfg$upper_rect[[1]], xmax = axis_cfg$upper_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#d1e5f0", alpha = 0.28) +
      geom_segment(
        data = ref_df,
        aes(x = boundary_x_start, xend = boundary_x_end, y = param_prototype, yend = param_prototype),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      ) +
      geom_vline(xintercept = axis_cfg$vlines, color = axis_cfg$vline_colors, linetype = axis_cfg$vline_linetypes, linewidth = 0.35)
  } else {
    p <- p +
      geom_segment(
        data = ref_df,
        aes(x = boundary_x_lower, xend = boundary_x_upper, y = param_prototype, yend = param_prototype),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      )
  }

  p <- p +
    geom_point(
      data = other_df,
      shape = 16,
      size = 2.1,
      alpha = 0.5,
      color = "grey65",
      position = point_pos
    ) +
    axis_cfg$scale +
    labs(
      title = paste0("Parameter Positions Within Fitted Bounds", if (!is.null(title_suffix) && nzchar(title_suffix)) paste0(" (", title_suffix, ")") else "", ": ", run_label),
      subtitle = paste0(
        if (identical(axis_cfg$axis_type, "relative")) {
          paste0(
            "0 = lower bound, 1 = upper bound; shaded zones are within ",
            sprintf("%.0f", 100 * near_thresh),
            "% of a bound | "
          )
        } else {
          "Horizontal lines span original lower-to-upper parameter bounds | "
        },
        top3_label,
        axis_cfg$subtitle_note
      ),
      x = axis_cfg$x_label,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
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
      scale_shape_manual(values = shape_values, breaks = top_breaks, drop = FALSE) +
      labs(shape = legend_title) +
      guides(shape = guide_legend(nrow = 1, byrow = TRUE))
  }

  ggplot2::ggsave(out_path, p, width = 12, height = 4.8)
  invisible(out_path)
}

plot_joint_soft_context_boundary_forest <- function(soft_df,
                                                    summary_df,
                                                    out_path,
                                                    run_label,
                                                    near_thresh = 0.05,
                                                    top3_seeds = NULL,
                                                    title_suffix = NULL,
                                                    legend_title = "Objective Top 3 Seeds",
                                                    x_scale = c("relative", "log10_original")) {
  x_scale <- match.arg(x_scale)
  required <- c(
    "seed",
    "parameter",
    "vivo_transformed",
    "vitro_transformed",
    "vivo_natural",
    "vitro_natural"
  )
  if (!is.data.frame(soft_df) || !nrow(soft_df) || !all(required %in% names(soft_df))) {
    return(invisible(NULL))
  }
  lower_t_col <- if ("joint_union_lower_transformed" %in% names(soft_df)) "joint_union_lower_transformed" else "center_lower_transformed"
  upper_t_col <- if ("joint_union_upper_transformed" %in% names(soft_df)) "joint_union_upper_transformed" else "center_upper_transformed"
  lower_nat_col <- if ("joint_union_lower_bound" %in% names(soft_df)) "joint_union_lower_bound" else "center_lower_bound"
  upper_nat_col <- if ("joint_union_upper_bound" %in% names(soft_df)) "joint_union_upper_bound" else "center_upper_bound"
  bound_cols <- c(lower_t_col, upper_t_col, lower_nat_col, upper_nat_col)
  if (!all(bound_cols %in% names(soft_df))) {
    return(invisible(NULL))
  }

  make_context_df <- function(context, value_t_col, value_nat_col) {
    data.frame(
      seed = as.character(soft_df$seed),
      parameter = as.character(soft_df$parameter),
      context = context,
      value_transformed = suppressWarnings(as.numeric(soft_df[[value_t_col]])),
      value_natural = suppressWarnings(as.numeric(soft_df[[value_nat_col]])),
      bound_lower_transformed = suppressWarnings(as.numeric(soft_df[[lower_t_col]])),
      bound_upper_transformed = suppressWarnings(as.numeric(soft_df[[upper_t_col]])),
      bound_lower_natural = suppressWarnings(as.numeric(soft_df[[lower_nat_col]])),
      bound_upper_natural = suppressWarnings(as.numeric(soft_df[[upper_nat_col]])),
      stringsAsFactors = FALSE
    )
  }
  plot_df <- dplyr::bind_rows(
    make_context_df("in vivo", "vivo_transformed", "vivo_natural"),
    make_context_df("in vitro", "vitro_transformed", "vitro_natural")
  )
  width_t <- plot_df$bound_upper_transformed - plot_df$bound_lower_transformed
  plot_df$rel_pos_in_range <- (plot_df$value_transformed - plot_df$bound_lower_transformed) / width_t
  plot_df$rel_pos_plot <- pmin(pmax(plot_df$rel_pos_in_range, 0), 1)
  plot_df$rel_dist_to_lower <- (plot_df$value_transformed - plot_df$bound_lower_transformed) / width_t
  plot_df$rel_dist_to_upper <- (plot_df$bound_upper_transformed - plot_df$value_transformed) / width_t
  plot_df$rel_dist_to_nearest <- pmin(plot_df$rel_dist_to_lower, plot_df$rel_dist_to_upper)
  plot_df <- plot_df[
    is.finite(plot_df$rel_pos_plot) &
      is.finite(plot_df$rel_dist_to_nearest) &
      is.finite(width_t) &
      width_t > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))

  if (is.null(top3_seeds)) top3_seeds <- get_top_ranked_seeds(summary_df, n = 3L)
  top3_seeds <- as.character(top3_seeds)

  axis_cfg <- boundary_axis_config(
    plot_df$rel_pos_plot,
    near_thresh = near_thresh,
    x_scale = x_scale,
    raw_value = plot_df$value_natural,
    raw_lower = plot_df$bound_lower_natural,
    raw_upper = plot_df$bound_upper_natural
  )
  plot_df$boundary_x_plot <- axis_cfg$x
  if (identical(axis_cfg$axis_type, "log10_original")) {
    plot_df$boundary_x_lower <- axis_cfg$lower_plot
    plot_df$boundary_x_upper <- axis_cfg$upper_plot
    plot_df <- plot_df[
      is.finite(plot_df$boundary_x_plot) &
        is.finite(plot_df$boundary_x_lower) &
        is.finite(plot_df$boundary_x_upper) &
        plot_df$boundary_x_upper > plot_df$boundary_x_lower,
      ,
      drop = FALSE
    ]
    if (!nrow(plot_df)) return(invisible(NULL))
  }

  plot_df$parameter_label <- as.character(plot_df$parameter)
  param_rank <- tapply(plot_df$rel_dist_to_nearest, plot_df$parameter_label, min, na.rm = TRUE)
  param_levels <- names(sort(param_rank, decreasing = FALSE))
  plot_df$parameter <- factor(plot_df$parameter_label, levels = rev(param_levels))
  y_breaks <- seq_along(levels(plot_df$parameter))
  plot_df$y_base <- as.numeric(plot_df$parameter)
  plot_df$context <- factor(plot_df$context, levels = c("in vivo", "in vitro"))
  pair_key <- paste(plot_df$seed, plot_df$parameter_label, sep = "\r")
  pair_hash <- vapply(
    pair_key,
    function(key) {
      ints <- utf8ToInt(key)
      if (!length(ints)) return(0)
      sum((seq_along(ints) * ints) %% 997)
    },
    numeric(1)
  )
  plot_df$pair_jitter <- ((pair_hash %% 101) / 100 - 0.5) * 0.08
  plot_df$context_offset <- ifelse(as.character(plot_df$context) == "in vivo", 0.18, -0.18)
  plot_df$y_plot <- plot_df$y_base + plot_df$context_offset + plot_df$pair_jitter
  ref_df <- if (identical(axis_cfg$axis_type, "log10_original")) {
    unique(plot_df[c("parameter", "y_base", "boundary_x_lower", "boundary_x_upper")])
  } else {
    ref <- unique(plot_df[c("parameter", "y_base")])
    ref$boundary_x_start <- axis_cfg$lower_limit
    ref$boundary_x_end <- axis_cfg$upper_limit
    ref
  }

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
  top3_label <- if (length(top3_seeds)) {
    paste(paste0("Top ", seq_along(top3_seeds), ": ", top3_seeds), collapse = "; ")
  } else {
    "No seeds met the 2N/4N 1000-day prediction gate."
  }
  vivo_pair_df <- plot_df[
    as.character(plot_df$context) == "in vivo",
    c("seed", "parameter_label", "boundary_x_plot", "y_plot", "seed_marker"),
    drop = FALSE
  ]
  vitro_pair_df <- plot_df[
    as.character(plot_df$context) == "in vitro",
    c("seed", "parameter_label", "boundary_x_plot", "y_plot", "seed_marker"),
    drop = FALSE
  ]
  pair_df <- merge(
    vivo_pair_df,
    vitro_pair_df,
    by = c("seed", "parameter_label"),
    suffixes = c("_vivo", "_vitro"),
    all = FALSE,
    sort = FALSE
  )
  pair_df$is_top_seed <- pair_df$seed %in% top3_seeds
  other_pair_df <- pair_df[!pair_df$is_top_seed, , drop = FALSE]
  top_pair_df <- pair_df[pair_df$is_top_seed, , drop = FALSE]
  title_detail <- if (!is.null(title_suffix) && nzchar(title_suffix)) title_suffix else ""
  title_text <- paste0(
    "Joint Soft-Coupled In Vivo/In Vitro Paired Parameter Positions",
    ": ",
    run_label
  )
  subtitle_line1 <- if (identical(axis_cfg$axis_type, "relative")) {
    paste0(
      "0 = joint union lower bound, 1 = joint union upper bound; shaded zones are within ",
      sprintf("%.0f", 100 * near_thresh),
      "% of a bound"
    )
  } else {
    "Horizontal lines span natural joint union lower-to-upper bounds"
  }
  subtitle_line2 <- paste(
    c(
      title_detail,
      "Context values are center +/- delta / 2; lines connect paired seed-parameter values",
      top3_label
    )[nzchar(c(
      title_detail,
      "Context values are center +/- delta / 2; lines connect paired seed-parameter values",
      top3_label
    ))],
    collapse = " | "
  )
  subtitle_lines <- c(
    subtitle_line1,
    subtitle_line2,
    trimws(axis_cfg$subtitle_note)
  )
  subtitle_lines <- subtitle_lines[nzchar(subtitle_lines)]
  subtitle_text <- paste(
    vapply(
      subtitle_lines,
      function(line) paste(strwrap(line, width = 135), collapse = "\n"),
      character(1)
    ),
    collapse = "\n"
  )

  p <- ggplot(plot_df, aes(x = boundary_x_plot, y = y_plot))
  if (identical(axis_cfg$axis_type, "relative")) {
    p <- p +
      annotate("rect", xmin = axis_cfg$lower_rect[[1]], xmax = axis_cfg$lower_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#fddbc7", alpha = 0.28) +
      annotate("rect", xmin = axis_cfg$upper_rect[[1]], xmax = axis_cfg$upper_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#d1e5f0", alpha = 0.28) +
      geom_segment(
        data = ref_df,
        aes(x = boundary_x_start, xend = boundary_x_end, y = y_base, yend = y_base),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      ) +
      geom_vline(xintercept = axis_cfg$vlines, color = axis_cfg$vline_colors, linetype = axis_cfg$vline_linetypes, linewidth = 0.35)
  } else {
    p <- p +
      geom_segment(
        data = ref_df,
        aes(x = boundary_x_lower, xend = boundary_x_upper, y = y_base, yend = y_base),
        inherit.aes = FALSE,
        color = "grey78",
        linewidth = 0.5
      )
  }

  if (nrow(other_pair_df)) {
    p <- p +
      geom_segment(
        data = other_pair_df,
        aes(
          x = boundary_x_plot_vivo,
          xend = boundary_x_plot_vitro,
          y = y_plot_vivo,
          yend = y_plot_vitro
        ),
        inherit.aes = FALSE,
        color = "grey55",
        alpha = 0.12,
        linewidth = 0.18
      )
  }
  if (nrow(top_pair_df)) {
    p <- p +
      geom_segment(
        data = top_pair_df,
        aes(
          x = boundary_x_plot_vivo,
          xend = boundary_x_plot_vitro,
          y = y_plot_vivo,
          yend = y_plot_vitro
        ),
        inherit.aes = FALSE,
        color = "grey25",
        alpha = 0.70,
        linewidth = 0.45
      )
  }

  p <- p +
    geom_point(
      data = other_df,
      aes(color = context),
      shape = 16,
      size = 2.1,
      alpha = 0.45
    ) +
    scale_color_manual(
      values = c("in vivo" = "#1b9e77", "in vitro" = "#d95f02"),
      breaks = c("in vivo", "in vitro"),
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = levels(plot_df$parameter),
      expand = ggplot2::expansion(add = 0.45)
    ) +
    axis_cfg$scale +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = axis_cfg$x_label,
      y = NULL,
      color = "Context"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  if (nrow(top_df)) {
    p <- p +
      geom_point(
        data = top_df,
        aes(shape = seed_marker),
        size = 3.0,
        color = "black"
      ) +
      scale_shape_manual(values = shape_values, breaks = top_breaks, drop = FALSE) +
      labs(shape = legend_title) +
      guides(
        color = guide_legend(nrow = 1, byrow = TRUE),
        shape = guide_legend(nrow = 1, byrow = TRUE)
      )
  } else {
    p <- p + guides(color = guide_legend(nrow = 1, byrow = TRUE))
  }

  ggplot2::ggsave(out_path, p, width = 12, height = 4.8 * 1.2)
  invisible(out_path)
}

plot_objective_vs_boundary_risk <- function(summary_df, out_path, run_label) {
  distance_col <- "min_rel_dist_active_excl_sigma_burden"
  if (!(distance_col %in% names(summary_df)) ||
      !any(is.finite(suppressWarnings(as.numeric(summary_df[[distance_col]]))) &
        suppressWarnings(as.numeric(summary_df[[distance_col]])) > 0, na.rm = TRUE)) {
    distance_col <- "min_rel_dist_active"
  }
  use_distance <- (distance_col %in% names(summary_df)) &&
    any(is.finite(suppressWarnings(as.numeric(summary_df[[distance_col]]))) &
      suppressWarnings(as.numeric(summary_df[[distance_col]])) > 0, na.rm = TRUE)
  if (isTRUE(use_distance)) {
    summary_df$boundary_risk_plot <- suppressWarnings(as.numeric(summary_df[[distance_col]]))
    keep_risk <- is.finite(summary_df$boundary_risk_plot) & summary_df$boundary_risk_plot > 0
    distance_label <- if (identical(distance_col, "min_rel_dist_active_excl_sigma_burden")) {
      "Boundary risk shown as minimum relative distance to a fitted bound, excluding sigma_burden"
    } else {
      "Boundary risk shown as minimum relative distance to a fitted active parameter bound"
    }
    y_label <- "Min relative distance to nearest bound (log10 scale)"
  } else {
    penalty_col <- "boundary_penalty_active"
    if (!(penalty_col %in% names(summary_df))) return(invisible(NULL))
    summary_df$boundary_risk_plot <- suppressWarnings(as.numeric(summary_df[[penalty_col]]))
    keep_risk <- is.finite(summary_df$boundary_risk_plot)
    distance_label <- "Boundary risk shown as active-parameter boundary penalty because all seeds have zero minimum relative distance."
    y_label <- "Active-parameter boundary penalty"
  }
  plot_df <- summary_df[
    is.finite(summary_df$objective) & keep_risk,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))
  plot_df <- plot_df[order(plot_df$objective), , drop = FALSE]
  plot_df$seed <- factor(plot_df$seed, levels = plot_df$seed)

  p <- ggplot(plot_df, aes(x = objective, y = boundary_risk_plot, label = seed)) +
    geom_point(size = 2.8, color = "#2c7fb8") +
    geom_text(nudge_y = 0.015, show.legend = FALSE, size = 3) +
    labs(
      title = paste0("Objective vs Boundary Risk: ", run_label),
      subtitle = distance_label,
      x = "Objective",
      y = y_label
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  if (isTRUE(use_distance)) {
    p <- p + scale_y_log10()
  }

  ggplot2::ggsave(out_path, p, width = 10, height = 7)
  invisible(out_path)
}

plot_joint_objective_components <- function(summary_df, out_path, run_label) {
  required <- c("seed", "objective", "objective_invivo", "objective_invitro")
  if (!all(required %in% names(summary_df))) return(invisible(NULL))
  plot_df <- summary_df[
    is.finite(summary_df$objective) |
      is.finite(summary_df$objective_invivo) |
      is.finite(summary_df$objective_invitro),
    required,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))
  plot_df <- plot_df[order(plot_df$objective, seed_order_key(plot_df$seed), plot_df$seed, na.last = TRUE), , drop = FALSE]
  plot_df$seed <- factor(as.character(plot_df$seed), levels = as.character(plot_df$seed))

  comp_long <- rbind(
    data.frame(seed = plot_df$seed, component = "Joint total", objective_value = plot_df$objective, stringsAsFactors = FALSE),
    data.frame(seed = plot_df$seed, component = "In vivo", objective_value = plot_df$objective_invivo, stringsAsFactors = FALSE),
    data.frame(seed = plot_df$seed, component = "In vitro", objective_value = plot_df$objective_invitro, stringsAsFactors = FALSE)
  )
  comp_long <- comp_long[is.finite(comp_long$objective_value), , drop = FALSE]
  if (!nrow(comp_long)) return(invisible(NULL))
  comp_long$component <- factor(comp_long$component, levels = c("Joint total", "In vivo", "In vitro"))
  component_colors <- c("Joint total" = "#4D4D4D", "In vivo" = "#1f77b4", "In vitro" = "#d95f02")

  p <- ggplot(comp_long, aes(x = component, y = objective_value, fill = component)) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_violin(width = 0.86, alpha = 0.68, trim = FALSE, color = "#53606b", linewidth = 0.35) +
    geom_boxplot(width = 0.2, outlier.size = 0.8, outlier.alpha = 0.5, fill = "white", color = "#202830", linewidth = 0.35) +
    scale_fill_manual(values = component_colors) +
    labs(
      title = paste0("Joint Objective Components: ", run_label),
      subtitle = "Across-seed distributions of joint total, in vivo, and in vitro objectives. Lower is better.",
      x = NULL,
      y = "Objective",
      fill = "Component"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 10)
    )
  ggplot2::ggsave(out_path, p, width = 9, height = 7)
  invisible(out_path)
}

plot_joint_objective_tradeoff <- function(summary_df, out_path, run_label) {
  required <- c("seed", "objective", "objective_invivo", "objective_invitro")
  if (!all(required %in% names(summary_df))) return(invisible(NULL))
  plot_df <- summary_df[
    is.finite(summary_df$objective_invivo) & is.finite(summary_df$objective_invitro),
    required,
    drop = FALSE
  ]
  if (!nrow(plot_df)) return(invisible(NULL))
  plot_df <- plot_df[order(plot_df$objective, seed_order_key(plot_df$seed), plot_df$seed, na.last = TRUE), , drop = FALSE]

  p <- ggplot(plot_df, aes(x = objective_invivo, y = objective_invitro, color = objective, label = seed)) +
    geom_point(size = 3) +
    geom_text(nudge_y = 0.02 * diff(range(plot_df$objective_invitro, na.rm = TRUE)), size = 3, show.legend = FALSE) +
    scale_color_gradient(low = "#2c7fb8", high = "#d7191c") +
    labs(
      title = paste0("Joint Objective Tradeoff: ", run_label),
      subtitle = "Each point is one seed; color shows total joint objective.",
      x = "In vivo objective",
      y = "In vitro objective",
      color = "Joint objective"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 9, height = 7)
  invisible(out_path)
}

plot_joint_soft_delta_magnitude <- function(soft_df, out_path, run_label) {
  required <- c("seed", "parameter", "delta_transformed")
  if (!is.data.frame(soft_df) || !all(required %in% names(soft_df))) return(invisible(NULL))
  plot_df <- soft_df[, required, drop = FALSE]
  plot_df$delta_transformed <- suppressWarnings(as.numeric(plot_df$delta_transformed))
  plot_df$abs_delta <- abs(plot_df$delta_transformed)
  plot_df <- plot_df[is.finite(plot_df$abs_delta), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  med <- stats::aggregate(abs_delta ~ parameter, data = plot_df, FUN = median, na.rm = TRUE)
  levels_use <- med$parameter[order(med$abs_delta, decreasing = TRUE)]
  plot_df$parameter <- factor(as.character(plot_df$parameter), levels = levels_use)
  p <- ggplot(plot_df, aes(x = parameter, y = abs_delta)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, fill = "#d9ead3", color = "#334155") +
    geom_point(aes(color = seed), size = 2.2, alpha = 0.8, position = position_jitter(width = 0.08, height = 0)) +
    coord_flip() +
    labs(
      title = paste0("Joint Soft-Coupling Delta Magnitudes: ", run_label),
      subtitle = "Absolute context split on the transformed optimizer scale.",
      x = NULL,
      y = "|delta| on transformed scale",
      color = "Seed"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 9.5, height = 6.2)
  invisible(out_path)
}

plot_joint_soft_penalty_ranking <- function(soft_df, out_path, run_label) {
  required <- c("seed", "parameter", "penalty_paid")
  if (!is.data.frame(soft_df) || !all(required %in% names(soft_df))) return(invisible(NULL))
  plot_df <- soft_df[, required, drop = FALSE]
  plot_df$penalty_paid <- suppressWarnings(as.numeric(plot_df$penalty_paid))
  plot_df <- plot_df[is.finite(plot_df$penalty_paid), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  med <- stats::aggregate(penalty_paid ~ parameter, data = plot_df, FUN = median, na.rm = TRUE)
  levels_use <- med$parameter[order(med$penalty_paid, decreasing = TRUE)]
  plot_df$parameter <- factor(as.character(plot_df$parameter), levels = levels_use)
  p <- ggplot(plot_df, aes(x = parameter, y = penalty_paid)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, fill = "#fee8c8", color = "#334155") +
    geom_point(aes(color = seed), size = 2.2, alpha = 0.8, position = position_jitter(width = 0.08, height = 0)) +
    coord_flip() +
    labs(
      title = paste0("Joint Soft-Coupling Penalty Ranking: ", run_label),
      subtitle = "Penalty contribution delta^2 / (2 sigma^2) by parameter.",
      x = NULL,
      y = "Penalty paid",
      color = "Seed"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 9.5, height = 6.2)
  invisible(out_path)
}

plot_joint_soft_vivo_vitro_pairs <- function(soft_df, out_path, run_label) {
  required <- c("seed", "parameter", "vivo_natural", "vitro_natural")
  if (!is.data.frame(soft_df) || !all(required %in% names(soft_df))) return(invisible(NULL))
  plot_df <- soft_df[, required, drop = FALSE]
  plot_df$vivo_natural <- suppressWarnings(as.numeric(plot_df$vivo_natural))
  plot_df$vitro_natural <- suppressWarnings(as.numeric(plot_df$vitro_natural))
  plot_df <- plot_df[is.finite(plot_df$vivo_natural) & is.finite(plot_df$vitro_natural), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  long_df <- rbind(
    data.frame(seed = plot_df$seed, parameter = plot_df$parameter, context = "in vivo", value = plot_df$vivo_natural, stringsAsFactors = FALSE),
    data.frame(seed = plot_df$seed, parameter = plot_df$parameter, context = "in vitro", value = plot_df$vitro_natural, stringsAsFactors = FALSE)
  )
  long_df$context <- factor(long_df$context, levels = c("in vivo", "in vitro"))
  p <- ggplot(long_df, aes(x = context, y = value, group = seed, color = seed)) +
    geom_line(alpha = 0.65, linewidth = 0.35) +
    geom_point(size = 2.2, alpha = 0.85) +
    facet_wrap(~ parameter, scales = "free_y") +
    labs(
      title = paste0("Joint Soft-Coupled In Vivo vs In Vitro Parameters: ", run_label),
      subtitle = "Natural-scale paired values derived from center +/- delta/2.",
      x = NULL,
      y = "Natural-scale value",
      color = "Seed"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 10.5, height = 7)
  invisible(out_path)
}

plot_invitro_objective_components <- function(summary_df, out_path, run_label) {
  required <- c("seed", "objective", "growth_loglik", "ploidy_loglik", "flow_loglik")
  if (!all(required %in% names(summary_df))) return(invisible(NULL))
  plot_df <- summary_df[is.finite(summary_df$objective), required, drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  plot_df <- plot_df[order(plot_df$objective, seed_order_key(plot_df$seed), plot_df$seed, na.last = TRUE), , drop = FALSE]
  plot_df$seed <- factor(as.character(plot_df$seed), levels = as.character(plot_df$seed))
  comp_long <- rbind(
    data.frame(seed = plot_df$seed, component = "Total objective", objective_value = plot_df$objective, stringsAsFactors = FALSE),
    data.frame(seed = plot_df$seed, component = "Growth -logLik", objective_value = -plot_df$growth_loglik, stringsAsFactors = FALSE),
    data.frame(seed = plot_df$seed, component = "Ploidy -logLik", objective_value = -plot_df$ploidy_loglik, stringsAsFactors = FALSE),
    data.frame(seed = plot_df$seed, component = "Flow -logLik", objective_value = -plot_df$flow_loglik, stringsAsFactors = FALSE)
  )
  comp_long <- comp_long[is.finite(comp_long$objective_value), , drop = FALSE]
  if (!nrow(comp_long)) return(invisible(NULL))
  comp_long$component <- factor(comp_long$component, levels = c("Total objective", "Growth -logLik", "Ploidy -logLik", "Flow -logLik"))
  n_seed <- length(levels(comp_long$seed))

  p <- ggplot(comp_long, aes(x = seed, y = objective_value, color = component, group = component)) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_line(linewidth = 0.45, alpha = 0.75) +
    geom_point(size = if (n_seed > 150L) 0.65 else 1.5, alpha = if (n_seed > 150L) 0.65 else 0.9) +
    scale_color_manual(values = c("Total objective" = "#4D4D4D", "Growth -logLik" = "#1f77b4", "Ploidy -logLik" = "#2ca02c", "Flow -logLik" = "#d95f02")) +
    labs(
      title = paste0("In Vitro Objective Components: ", run_label),
      subtitle = "Seeds are ordered by total objective. Lower total objective is better; component traces are -logLik contributions.",
      x = "Seed rank by total objective",
      y = "Objective contribution",
      color = "Component"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = if (n_seed > 80L) element_blank() else element_text(angle = 60, hjust = 1),
      axis.ticks.x = if (n_seed > 80L) element_blank() else element_line()
    )
  ggplot2::ggsave(out_path, p, width = 12, height = 4.2)
  invisible(out_path)
}

prepare_invitro_summary_for_plot <- function(seed_summary, is_joint_run = FALSE) {
  plot_summary <- seed_summary
  if (isTRUE(is_joint_run) && "objective_invitro" %in% names(plot_summary)) {
    invitro_objective <- suppressWarnings(as.numeric(plot_summary$objective_invitro))
    plot_summary$objective <- invitro_objective
    plot_summary$objective_total <- invitro_objective
    plot_summary$total_loglik <- -invitro_objective
    plot_summary$fit_mode <- "fit_invitro"
  }
  plot_summary$objective <- suppressWarnings(as.numeric(plot_summary$objective))
  plot_summary$objective_total <- suppressWarnings(as.numeric(plot_summary$objective_total))
  plot_summary$objective_rank <- rank(plot_summary$objective, ties.method = "first", na.last = "keep")
  plot_summary[order(plot_summary$objective, seed_order_key(plot_summary$seed), plot_summary$seed, na.last = TRUE), , drop = FALSE]
}

prepare_joint_invitro_fit_summary_values <- function(fit_summary_vals) {
  out <- fit_summary_vals
  objective_invitro <- as_num(summary_metric_value(out, "objective_invitro", NA_real_), NA_real_)
  out[["fit_mode"]] <- "fit_invitro"
  if (is.finite(objective_invitro)) {
    out[["objective"]] <- as.character(objective_invitro)
    out[["objective_total"]] <- as.character(objective_invitro)
    out[["total_loglik"]] <- as.character(-objective_invitro)
  }
  out
}

plot_invitro_objective_component_distributions <- function(summary_df, out_path, out_table, run_label) {
  required <- c("seed", "objective", "growth_loglik", "ploidy_loglik", "flow_loglik")
  if (!all(required %in% names(summary_df))) return(invisible(NULL))
  plot_df <- summary_df
  plot_df$objective <- suppressWarnings(as.numeric(plot_df$objective))
  plot_df$growth_loglik <- suppressWarnings(as.numeric(plot_df$growth_loglik))
  plot_df$ploidy_loglik <- suppressWarnings(as.numeric(plot_df$ploidy_loglik))
  plot_df$flow_loglik <- suppressWarnings(as.numeric(plot_df$flow_loglik))
  rows <- list(
    data.frame(seed = as.character(plot_df$seed), metric = "objective_total_plot", metric_label = "objective_total", value = plot_df$objective, stringsAsFactors = FALSE),
    data.frame(seed = as.character(plot_df$seed), metric = "growth_negloglik", metric_label = "growth_-logLik", value = -plot_df$growth_loglik, stringsAsFactors = FALSE),
    data.frame(seed = as.character(plot_df$seed), metric = "ploidy_negloglik", metric_label = "ploidy_-logLik", value = -plot_df$ploidy_loglik, stringsAsFactors = FALSE),
    data.frame(seed = as.character(plot_df$seed), metric = "flow_negloglik", metric_label = "flow_-logLik", value = -plot_df$flow_loglik, stringsAsFactors = FALSE)
  )
  comp_long <- do.call(rbind, rows)
  comp_long <- comp_long[is.finite(comp_long$value), , drop = FALSE]
  if (!nrow(comp_long)) return(invisible(NULL))
  metric_levels <- c("objective_total_plot", "growth_negloglik", "ploidy_negloglik", "flow_negloglik")
  label_levels <- c("objective_total", "growth_-logLik", "ploidy_-logLik", "flow_-logLik")
  comp_long$metric <- factor(comp_long$metric, levels = metric_levels)
  comp_long$metric_label <- factor(comp_long$metric_label, levels = label_levels)
  utils::write.table(
    comp_long,
    file = out_table,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p <- ggplot(comp_long, aes(x = metric_label, y = value, fill = metric_label)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, linewidth = 0.35) +
    scale_fill_manual(
      values = stats::setNames(c("#4c78a8", "#f58518", "#54a24b", "#d95f02"), label_levels),
      guide = "none"
    ) +
    labs(
      title = "Objective Distributions Across Seeds",
      subtitle = paste0("Source: ", run_label, " joint-fit in vitro components"),
      x = NULL,
      y = "Objective value"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 10, height = 7)
  invisible(out_path)
}

blank_extra_diagnostic_plot <- function(title, message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message, size = 4, color = "#4b5563") +
    labs(title = title) +
    theme_void(base_size = 11)
}

plot_invitro_optimization_diagnostics <- function(summary_df,
                                                  parameter_long,
                                                  out_path,
                                                  run_label,
                                                  near_thresh = 0.05,
                                                  parameter_x_scale = c("relative", "log10_original"),
                                                  parameter_panel_only = FALSE) {
  parameter_x_scale <- match.arg(parameter_x_scale)
  if (is.null(summary_df) || !nrow(summary_df) || !"seed" %in% names(summary_df) || !"objective" %in% names(summary_df)) {
    return(invisible(NULL))
  }
  plot_df <- summary_df[is.finite(summary_df$objective), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  plot_df <- plot_df[order(plot_df$objective, seed_order_key(plot_df$seed), plot_df$seed, na.last = TRUE), , drop = FALSE]
  plot_df$objective_rank <- seq_len(nrow(plot_df))
  plot_df$seed_marker <- "Other seeds"
  plot_df$seed_marker[plot_df$objective_rank == 1L] <- paste0("Best: ", plot_df$seed[[1]])
  if (nrow(plot_df) >= 2L) plot_df$seed_marker[plot_df$objective_rank == 2L] <- paste0("Rank 2: ", plot_df$seed[[2]])
  if (nrow(plot_df) >= 3L) plot_df$seed_marker[plot_df$objective_rank == 3L] <- paste0("Rank 3: ", plot_df$seed[[3]])
  marker_levels <- c(paste0("Best: ", plot_df$seed[[1]]))
  if (nrow(plot_df) >= 2L) marker_levels <- c(marker_levels, paste0("Rank 2: ", plot_df$seed[[2]]))
  if (nrow(plot_df) >= 3L) marker_levels <- c(marker_levels, paste0("Rank 3: ", plot_df$seed[[3]]))
  marker_levels <- c(marker_levels, "Other seeds")
  plot_df$seed_marker <- factor(plot_df$seed_marker, levels = marker_levels)

  p_rank <- ggplot(plot_df, aes(x = objective_rank, y = objective)) +
    geom_line(color = "grey55", linewidth = 0.5) +
    geom_point(aes(color = seed_marker), size = 2.3, alpha = 0.9) +
    scale_color_manual(
      values = c(
        setNames("#1b9e77", marker_levels[[1]]),
        if (length(marker_levels) >= 3L) setNames("#377eb8", marker_levels[[2]]) else NULL,
        if (length(marker_levels) >= 4L) setNames("#4e79a7", marker_levels[[3]]) else NULL,
        "Other seeds" = "grey70"
      ),
      drop = FALSE
    ) +
    labs(
      title = "Objective by Seed Rank",
      subtitle = paste0("Seeds ordered by total in vitro objective | run=", run_label),
      x = "Seed rank",
      y = "Total objective",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  best <- plot_df[1L, , drop = FALSE]
  component_names <- c(
    "Total objective" = "objective",
    "Objective total" = "objective_total",
    "Growth -logLik" = "growth_loglik",
    "Ploidy -logLik" = "ploidy_loglik",
    "Flow -logLik" = "flow_loglik"
  )
  comp_rows <- lapply(names(component_names), function(label) {
    col <- component_names[[label]]
    if (!(col %in% names(best))) return(NULL)
    value <- suppressWarnings(as.numeric(best[[col]][[1]]))
    if (!is.finite(value)) return(NULL)
    if (grepl("-logLik", label, fixed = TRUE)) value <- -value
    data.frame(component = label, value = value, stringsAsFactors = FALSE)
  })
  comp_df <- do.call(rbind, comp_rows[!vapply(comp_rows, is.null, logical(1))])
  if (!is.null(comp_df) && nrow(comp_df)) {
    comp_df$component <- factor(comp_df$component, levels = rev(comp_df$component))
    p_components <- ggplot(comp_df, aes(x = component, y = value)) +
      geom_col(fill = "#4B6F8A", width = 0.72) +
      coord_flip() +
      labs(
        title = paste0("Best Seed Objective Components: ", best$seed[[1]]),
        x = NULL,
        y = "Objective-scale value"
      ) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  } else {
    p_components <- blank_extra_diagnostic_plot(
      "Best Seed Objective Components",
      "Objective component columns were unavailable."
    )
  }

  param_plot <- parameter_long[
    parameter_long$active_in_fit &
      is.finite(parameter_long$rel_pos_plot) &
      !is.na(parameter_long$param_prototype),
    ,
    drop = FALSE
  ]
  if (nrow(param_plot)) {
    param_plot$objective_rank <- plot_df$objective_rank[match(param_plot$seed, plot_df$seed)]
    param_plot <- param_plot[is.finite(param_plot$objective_rank), , drop = FALSE]
  }
  if (nrow(param_plot)) {
    param_rank <- tapply(param_plot$rel_dist_to_nearest, param_plot$param_prototype, min, na.rm = TRUE)
    param_levels <- names(sort(param_rank, decreasing = FALSE))
    param_plot$param_prototype <- factor(param_plot$param_prototype, levels = rev(param_levels))
    axis_cfg <- boundary_axis_config(
      param_plot$rel_pos_plot,
      near_thresh = near_thresh,
      x_scale = parameter_x_scale,
      raw_value = param_plot$prototype_value,
      raw_lower = param_plot$prototype_lower_bound,
      raw_upper = param_plot$prototype_upper_bound
    )
    param_plot$boundary_x_plot <- axis_cfg$x
    if (identical(axis_cfg$axis_type, "log10_original")) {
      param_plot$boundary_x_lower <- axis_cfg$lower_plot
      param_plot$boundary_x_upper <- axis_cfg$upper_plot
      param_plot$boundary_x_lower_near <- axis_cfg$lower_near_plot
      param_plot$boundary_x_upper_near <- axis_cfg$upper_near_plot
      param_plot <- param_plot[
        is.finite(param_plot$boundary_x_plot) &
          is.finite(param_plot$boundary_x_lower) &
          is.finite(param_plot$boundary_x_upper) &
          param_plot$boundary_x_upper > param_plot$boundary_x_lower,
        ,
        drop = FALSE
      ]
      if (!nrow(param_plot)) {
        p_params <- blank_extra_diagnostic_plot(
          "Fitted Parameter Positions Across Seeds",
          "No active fitted parameter boundary positions were available."
        )
      }
    }
  }
  if (nrow(param_plot)) {
    param_plot$rank_group <- ifelse(param_plot$objective_rank == 1L, "Best seed", "Other seeds")
    best_param_plot <- param_plot[param_plot$objective_rank == 1L, , drop = FALSE]
    best_param_plot <- best_param_plot[!duplicated(as.character(best_param_plot$param_prototype)), , drop = FALSE]
    best_param_plot <- best_param_plot[order(as.character(best_param_plot$param_prototype)), , drop = FALSE]
    best_position_subtitle <- ""
    if (nrow(best_param_plot)) {
      best_position_values <- if (identical(axis_cfg$axis_type, "log10_original")) {
        formatC(best_param_plot$prototype_value, format = "fg", digits = 4)
      } else {
        formatC(best_param_plot$rel_pos_plot, format = "f", digits = 3)
      }
      best_position_text <- paste0(
        as.character(best_param_plot$param_prototype),
        "=",
        best_position_values
      )
      best_position_subtitle <- paste(
        strwrap(
          paste0(
            "Best seed ",
            best$seed[[1]],
            if (identical(axis_cfg$axis_type, "log10_original")) " fitted raw values: " else " fitted positions: ",
            paste(best_position_text, collapse = "; ")
          ),
          width = 120
        ),
        collapse = "\n"
      )
    }
    param_subtitle <- paste0(
      "Green points are the best seed; ",
      if (identical(axis_cfg$axis_type, "relative")) "shaded zones" else "colored end segments",
      " are within ",
      sprintf("%.0f", 100 * near_thresh),
      "% of a fitted bound."
    )
    if (nzchar(best_position_subtitle)) {
      param_subtitle <- paste(param_subtitle, best_position_subtitle, sep = "\n")
    }
    if (nzchar(axis_cfg$subtitle_note)) {
      param_subtitle <- paste0(param_subtitle, axis_cfg$subtitle_note)
    }
    ref_df <- if (identical(axis_cfg$axis_type, "log10_original")) {
      unique(param_plot[c("param_prototype", "boundary_x_lower", "boundary_x_upper", "boundary_x_lower_near", "boundary_x_upper_near")])
    } else {
      ref <- unique(param_plot["param_prototype"])
      ref$boundary_x_start <- axis_cfg$lower_limit
      ref$boundary_x_end <- axis_cfg$upper_limit
      ref
    }
    p_params <- ggplot(param_plot, aes(x = boundary_x_plot, y = param_prototype))
    if (identical(axis_cfg$axis_type, "relative")) {
      p_params <- p_params +
        annotate("rect", xmin = axis_cfg$lower_rect[[1]], xmax = axis_cfg$lower_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#fddbc7", alpha = 0.28) +
        annotate("rect", xmin = axis_cfg$upper_rect[[1]], xmax = axis_cfg$upper_rect[[2]], ymin = -Inf, ymax = Inf, fill = "#d1e5f0", alpha = 0.28) +
        geom_segment(
          data = ref_df,
          aes(x = boundary_x_start, xend = boundary_x_end, y = param_prototype, yend = param_prototype),
          inherit.aes = FALSE,
          color = "grey78",
          linewidth = 0.5
        ) +
        geom_vline(xintercept = axis_cfg$vlines, color = axis_cfg$vline_colors, linetype = axis_cfg$vline_linetypes, linewidth = 0.35)
    } else {
      p_params <- p_params +
        geom_segment(
          data = ref_df,
          aes(x = boundary_x_lower, xend = boundary_x_upper, y = param_prototype, yend = param_prototype),
          inherit.aes = FALSE,
          color = "grey78",
          linewidth = 0.5
        ) +
        geom_segment(
          data = ref_df,
          aes(x = boundary_x_lower, xend = boundary_x_lower_near, y = param_prototype, yend = param_prototype),
          inherit.aes = FALSE,
          color = "#fddbc7",
          linewidth = 2.4,
          alpha = 0.68
        ) +
        geom_segment(
          data = ref_df,
          aes(x = boundary_x_upper_near, xend = boundary_x_upper, y = param_prototype, yend = param_prototype),
          inherit.aes = FALSE,
          color = "#d1e5f0",
          linewidth = 2.4,
          alpha = 0.68
        )
    }
    p_params <- p_params +
      geom_point(
        data = param_plot[param_plot$rank_group == "Other seeds", , drop = FALSE],
        color = "grey70",
        size = 1.7,
        alpha = 0.5,
        position = position_jitter(height = 0.13, width = 0)
      ) +
      geom_point(
        data = param_plot[param_plot$rank_group == "Best seed", , drop = FALSE],
        color = "#1b9e77",
        size = 2.8,
        alpha = 0.95,
        position = position_jitter(height = 0.13, width = 0)
      ) +
      axis_cfg$scale +
      labs(
        title = "Fitted Parameter Positions Across Seeds",
        subtitle = param_subtitle,
        x = axis_cfg$x_label,
        y = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  } else if (!exists("p_params", inherits = FALSE)) {
    p_params <- blank_extra_diagnostic_plot(
      "Fitted Parameter Positions Across Seeds",
      "No active fitted parameter boundary positions were available."
    )
  }

  if (isTRUE(parameter_panel_only)) {
    ggplot2::ggsave(out_path, p_params, width = 12, height = 6.8)
    return(invisible(out_path))
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- (p_rank + p_components) / p_params +
      patchwork::plot_layout(heights = c(1, 1.55)) +
      patchwork::plot_annotation(
        title = paste0(
          "Optimization diagnostics",
          if (identical(parameter_x_scale, "log10_original")) " (parameter x-axis raw log10 scale)" else ""
        )
      )
    ggplot2::ggsave(out_path, p, width = 13, height = 10.5)
  } else {
    ggplot2::ggsave(out_path, p_rank, width = 11, height = 7)
  }
  invisible(out_path)
}

derive_invitro_lineage_label <- function(key) {
  key <- as.character(key)
  out <- rep(NA_character_, length(key))
  ok <- !is.na(key) & nzchar(key)
  if (!any(ok)) return(out)
  parts <- strsplit(key[ok], "_", fixed = TRUE)
  is_control <- vapply(
    parts,
    function(x) length(x) > 0L && all(trimws(x) == "20.5"),
    logical(1)
  )
  out[ok] <- ifelse(is_control, "control", "deprived")
  out
}

derive_invitro_lineage_passage_index <- function(key, fallback = NA_real_) {
  key <- as.character(key)
  out <- rep(NA_real_, length(key))
  ok <- !is.na(key) & nzchar(key)
  if (any(ok)) {
    out[ok] <- vapply(strsplit(key[ok], "_", fixed = TRUE), length, integer(1))
  }
  fallback <- suppressWarnings(as.numeric(fallback))
  use_fallback <- !is.finite(out) & is.finite(fallback)
  out[use_fallback] <- fallback[use_fallback]
  out
}

ensure_invitro_likelihood_columns <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  n <- nrow(df)
  if (!"lineage_terminal_key" %in% names(df)) {
    df$lineage_terminal_key <- if ("segment_id" %in% names(df)) as.character(df$segment_id) else as.character(seq_len(n))
  }
  if (!"lineage_passage_index" %in% names(df)) {
    df$lineage_passage_index <- if ("passage_index" %in% names(df)) suppressWarnings(as.numeric(df$passage_index)) else seq_len(n)
  }
  if (!"lineage_label" %in% names(df) || all(is.na(df$lineage_label))) {
    df$lineage_label <- derive_invitro_lineage_label(df$lineage_terminal_key)
  }
  missing_label <- is.na(df$lineage_label) | !nzchar(as.character(df$lineage_label))
  if (any(missing_label)) {
    df$lineage_label[missing_label] <- if ("cohort" %in% names(df)) as.character(df$cohort[missing_label]) else "lineage"
  }
  if (!"cohort" %in% names(df)) df$cohort <- "fit"
  if (!"segment_id" %in% names(df)) df$segment_id <- df$lineage_terminal_key
  df
}

read_seed_likelihood_table <- function(seed_dir, filename, fit_label, seed) {
  path <- file.path(seed_dir, filename)
  if (!file.exists(path)) return(NULL)
  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || !nrow(tab)) return(NULL)
  tab$fit_label <- fit_label
  tab$seed <- as.character(seed)
  ensure_invitro_likelihood_columns(tab)
}

plot_invitro_likelihood_comparison <- function(plot_df,
                                               out_path,
                                               title,
                                               y_label,
                                               size_col,
                                               size_label) {
  if (is.null(plot_df) || !nrow(plot_df) || !"mean_loglik" %in% names(plot_df)) return(invisible(NULL))
  plot_df$mean_loglik <- suppressWarnings(as.numeric(plot_df$mean_loglik))
  plot_df$lineage_passage_index <- suppressWarnings(as.numeric(plot_df$lineage_passage_index))
  if (size_col %in% names(plot_df)) {
    plot_df$size_value <- suppressWarnings(as.numeric(plot_df[[size_col]]))
  } else {
    plot_df$size_value <- NA_real_
  }
  plot_df <- plot_df[is.finite(plot_df$mean_loglik) & is.finite(plot_df$lineage_passage_index), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))
  fit_levels <- unique(as.character(plot_df$fit_label))
  plot_df$fit_label <- factor(as.character(plot_df$fit_label), levels = fit_levels)
  palette <- stats::setNames(c("#1b9e77", "#e15759", "#4e79a7", "#f28e2b", "#76b7b2")[seq_along(fit_levels)], fit_levels)

  p <- ggplot(plot_df, aes(x = lineage_passage_index, y = mean_loglik, color = fit_label)) +
    geom_line(aes(group = interaction(fit_label, cohort, lineage_label)), linewidth = 0.75, alpha = 0.82) +
    geom_point(aes(size = size_value), alpha = 0.78) +
    facet_grid(cohort ~ lineage_label, scales = "free_x", space = "free_x") +
    scale_color_manual(values = palette, drop = FALSE) +
    labs(
      title = title,
      x = "Lineage passage index",
      y = y_label,
      color = "Parameter set",
      size = size_label
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  ggplot2::ggsave(out_path, p, width = 12, height = 7.2)
  invisible(out_path)
}

plot_invitro_best_fit_likelihood_comparison <- function(seed_dirs,
                                                        seed_summary,
                                                        out_dir,
                                                        run_label,
                                                        max_fits = 2L) {
  if (is.null(seed_summary) || !nrow(seed_summary) || !"seed" %in% names(seed_summary)) {
    return(character(0))
  }
  seed_path <- stats::setNames(seed_dirs, basename(seed_dirs))
  available <- names(seed_path)[
    file.exists(file.path(seed_path, "invitro_ploidy_loglik.tsv")) |
      file.exists(file.path(seed_path, "invitro_flow_loglik.tsv"))
  ]
  if (!length(available)) return(character(0))
  candidates <- seed_summary[as.character(seed_summary$seed) %in% available, , drop = FALSE]
  if (!nrow(candidates)) return(character(0))
  candidates$objective <- suppressWarnings(as.numeric(candidates$objective))
  candidates <- candidates[order(candidates$objective, seed_order_key(candidates$seed), candidates$seed, na.last = TRUE), , drop = FALSE]
  selected <- utils::head(as.character(candidates$seed), max(1L, as.integer(max_fits)))
  labels <- vapply(seq_along(selected), function(i) {
    seed <- selected[[i]]
    if (i == 1L) return(paste0("Best: ", seed))
    rank <- candidates$objective_rank[match(seed, candidates$seed)]
    if (length(rank) && is.finite(suppressWarnings(as.numeric(rank)))) {
      paste0("Rank ", as.integer(rank), ": ", seed)
    } else {
      paste0("Comparison: ", seed)
    }
  }, character(1))
  names(labels) <- selected

  summary_cols <- intersect(
    c(
      "seed", "objective_rank", "objective", "objective_total", "total_loglik",
      "growth_loglik", "growth_loglik_sum", "ploidy_loglik", "ploidy_loglik_sum",
      "flow_loglik", "flow_loglik_sum", "n_growth", "n_ploidy_passages",
      "n_kary_cells", "n_flow_passages"
    ),
    names(candidates)
  )
  comparison_summary <- candidates[match(selected, candidates$seed), summary_cols, drop = FALSE]
  comparison_summary$fit_label <- labels[as.character(comparison_summary$seed)]
  comparison_summary <- comparison_summary[, c("fit_label", setdiff(names(comparison_summary), "fit_label")), drop = FALSE]
  utils::write.table(
    comparison_summary,
    file = file.path(out_dir, "best_fit_likelihood_comparison.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ploidy_rows <- lapply(selected, function(seed) {
    read_seed_likelihood_table(seed_path[[seed]], "invitro_ploidy_loglik.tsv", labels[[seed]], seed)
  })
  flow_rows <- lapply(selected, function(seed) {
    read_seed_likelihood_table(seed_path[[seed]], "invitro_flow_loglik.tsv", labels[[seed]], seed)
  })
  ploidy_df <- do.call(rbind, ploidy_rows[!vapply(ploidy_rows, is.null, logical(1))])
  flow_df <- do.call(rbind, flow_rows[!vapply(flow_rows, is.null, logical(1))])
  ploidy_long_path <- file.path(out_dir, "best_fit_ploidy_likelihood_comparison_long.tsv")
  flow_long_path <- file.path(out_dir, "best_fit_flow_likelihood_comparison_long.tsv")
  if (!is.null(ploidy_df) && nrow(ploidy_df)) {
    utils::write.table(
      ploidy_df,
      file = ploidy_long_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    ploidy_df <- utils::read.delim(ploidy_long_path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.null(flow_df) && nrow(flow_df)) {
    utils::write.table(
      flow_df,
      file = flow_long_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    flow_df <- utils::read.delim(flow_long_path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  written <- character(0)
  ploidy_out <- plot_invitro_likelihood_comparison(
    plot_df = ploidy_df,
    out_path = file.path(out_dir, "best_fit_ploidy_likelihood_comparison.pdf"),
    title = "Passage-level ploidy likelihood comparison",
    y_label = "Mean ploidy log-likelihood",
    size_col = "n_cells",
    size_label = "Observed cells"
  )
  if (!is.null(ploidy_out) && file.exists(ploidy_out)) written <- c(written, ploidy_out)
  flow_out <- plot_invitro_likelihood_comparison(
    plot_df = flow_df,
    out_path = file.path(out_dir, "best_fit_flow_likelihood_comparison.pdf"),
    title = "Passage-level flow-density likelihood comparison",
    y_label = "Mean flow-density log-likelihood",
    size_col = "n_grid",
    size_label = "Grid points"
  )
  if (!is.null(flow_out) && file.exists(flow_out)) written <- c(written, flow_out)
  unique(written)
}

read_seed_distribution_quantiles_table <- function(seed_dir, seed) {
  path <- file.path(seed_dir, "invitro_distribution_quantiles.tsv")
  if (!file.exists(path)) return(NULL)
  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  required <- c("cohort", "passage_index", "quantile_prob", "predicted_quantile_kary_N")
  if (is.null(tab) || !nrow(tab) || !all(required %in% names(tab))) return(NULL)

  tab$seed <- as.character(seed)
  tab$cohort <- as.character(tab$cohort)
  tab$passage_index <- suppressWarnings(as.numeric(tab$passage_index))
  tab$quantile_prob <- suppressWarnings(as.numeric(tab$quantile_prob))
  tab$predicted_quantile_kary_N <- suppressWarnings(as.numeric(tab$predicted_quantile_kary_N))
  if (!"segment_id" %in% names(tab)) tab$segment_id <- NA_character_
  if (!"oxygen_pct" %in% names(tab)) tab$oxygen_pct <- NA_real_
  if (!"selected_day" %in% names(tab)) tab$selected_day <- NA_real_
  if (!"lineage_label" %in% names(tab)) tab$lineage_label <- derive_invitro_lineage_label(tab$segment_id)
  if (!"lineage_passage_index" %in% names(tab)) {
    tab$lineage_passage_index <- derive_invitro_lineage_passage_index(tab$segment_id, tab$passage_index)
  }
  tab$oxygen_pct <- suppressWarnings(as.numeric(tab$oxygen_pct))
  tab$selected_day <- suppressWarnings(as.numeric(tab$selected_day))
  tab$lineage_passage_index <- suppressWarnings(as.numeric(tab$lineage_passage_index))
  missing_label <- is.na(tab$lineage_label) | !nzchar(as.character(tab$lineage_label))
  if (any(missing_label)) tab$lineage_label[missing_label] <- as.character(tab$cohort[missing_label])

  keep <- is.finite(tab$passage_index) &
    is.finite(tab$lineage_passage_index) &
    is.finite(tab$quantile_prob) &
    is.finite(tab$predicted_quantile_kary_N)
  tab <- tab[keep, , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  keep_cols <- intersect(
    c(
      "seed",
      "segment_id",
      "cohort",
      "lineage_label",
      "lineage_passage_index",
      "passage_index",
      "oxygen_pct",
      "selected_day",
      "quantile_prob",
      "predicted_quantile_kary_N"
    ),
    names(tab)
  )
  tab[, keep_cols, drop = FALSE]
}

read_seed_observed_kary_table <- function(seed_dir, seed) {
  path <- file.path(seed_dir, "invitro_observed_kary.tsv")
  if (!file.exists(path)) return(NULL)
  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  required <- c("cohort", "passage_index", "observed_kary_N")
  if (is.null(tab) || !nrow(tab) || !all(required %in% names(tab))) return(NULL)

  tab$seed <- as.character(seed)
  tab$cohort <- as.character(tab$cohort)
  tab$passage_index <- suppressWarnings(as.numeric(tab$passage_index))
  tab$observed_kary_N <- suppressWarnings(as.numeric(tab$observed_kary_N))
  if (!"segment_id" %in% names(tab)) tab$segment_id <- NA_character_
  if (!"lineage_label" %in% names(tab)) tab$lineage_label <- derive_invitro_lineage_label(tab$segment_id)
  if (!"lineage_passage_index" %in% names(tab)) {
    tab$lineage_passage_index <- derive_invitro_lineage_passage_index(tab$segment_id, tab$passage_index)
  }
  tab$lineage_passage_index <- suppressWarnings(as.numeric(tab$lineage_passage_index))
  missing_label <- is.na(tab$lineage_label) | !nzchar(as.character(tab$lineage_label))
  if (any(missing_label)) tab$lineage_label[missing_label] <- as.character(tab$cohort[missing_label])
  if (!"passage_id" %in% names(tab)) tab$passage_id <- NA_character_
  if (!"cell_index" %in% names(tab)) tab$cell_index <- NA_integer_

  keep <- is.finite(tab$passage_index) & is.finite(tab$lineage_passage_index) & is.finite(tab$observed_kary_N)
  tab <- tab[keep, , drop = FALSE]
  if (!nrow(tab)) return(NULL)

  keep_cols <- intersect(
    c(
      "segment_id",
      "cohort",
      "lineage_label",
      "lineage_passage_index",
      "passage_index",
      "passage_id",
      "cell_index",
      "observed_kary_N"
    ),
    names(tab)
  )
  unique(tab[, keep_cols, drop = FALSE])
}

plot_invitro_distribution_quantiles_multiseed <- function(seed_dirs,
                                                          out_path,
                                                          out_table,
                                                          run_label) {
  seed_path <- stats::setNames(seed_dirs, basename(seed_dirs))
  rows <- lapply(names(seed_path), function(seed) {
    read_seed_distribution_quantiles_table(seed_path[[seed]], seed)
  })
  quantile_df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(quantile_df) || !nrow(quantile_df)) return(invisible(NULL))

  obs_rows <- lapply(names(seed_path), function(seed) {
    read_seed_observed_kary_table(seed_path[[seed]], seed)
  })
  observed_df <- do.call(rbind, obs_rows[!vapply(obs_rows, is.null, logical(1))])
  if (!is.null(observed_df) && nrow(observed_df)) {
    observed_df <- unique(observed_df)
  }

  utils::write.table(
    quantile_df,
    file = out_table,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  mean_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else mean(x)
  }

  seed_quantiles <- stats::aggregate(
    predicted_quantile_kary_N ~ seed + cohort + lineage_label + lineage_passage_index + passage_index + quantile_prob,
    data = quantile_df,
    FUN = mean_finite
  )
  seed_quantiles <- seed_quantiles[is.finite(seed_quantiles$predicted_quantile_kary_N), , drop = FALSE]
  if (!nrow(seed_quantiles)) return(invisible(NULL))
  seed_quantiles <- seed_quantiles[
    order(seed_order_key(seed_quantiles$seed), seed_quantiles$seed, seed_quantiles$cohort, seed_quantiles$lineage_label, seed_quantiles$lineage_passage_index, seed_quantiles$quantile_prob),
    ,
    drop = FALSE
  ]

  plot_df <- stats::aggregate(
    predicted_quantile_kary_N ~ cohort + lineage_label + lineage_passage_index + quantile_prob,
    data = seed_quantiles,
    FUN = mean_finite
  )
  plot_df <- plot_df[is.finite(plot_df$predicted_quantile_kary_N), , drop = FALSE]
  if (!nrow(plot_df)) return(invisible(NULL))

  cohort_levels <- unique(c("2N", "4N", as.character(plot_df$cohort)))
  cohort_levels <- cohort_levels[cohort_levels %in% as.character(plot_df$cohort)]
  if (length(cohort_levels)) {
    plot_df$cohort <- factor(as.character(plot_df$cohort), levels = cohort_levels)
  }
  lineage_levels <- unique(c("control", "deprived", as.character(plot_df$lineage_label)))
  lineage_levels <- lineage_levels[lineage_levels %in% as.character(plot_df$lineage_label)]
  if (length(lineage_levels)) {
    plot_df$lineage_label <- factor(as.character(plot_df$lineage_label), levels = lineage_levels)
  }
  if (!is.null(observed_df) && nrow(observed_df)) {
    if (length(cohort_levels)) observed_df$cohort <- factor(as.character(observed_df$cohort), levels = cohort_levels)
    if (length(lineage_levels)) observed_df$lineage_label <- factor(as.character(observed_df$lineage_label), levels = lineage_levels)
  }
  plot_df$fit_label <- "Across-seed mean"

  p <- ggplot() +
    geom_line(
      data = plot_df,
      aes(
        x = lineage_passage_index,
        y = predicted_quantile_kary_N,
        color = fit_label,
        group = interaction(cohort, lineage_label, quantile_prob, drop = TRUE)
      ),
      linewidth = 0.7,
      alpha = 0.55
    )
  if (!is.null(observed_df) && nrow(observed_df)) {
    p <- p +
      geom_point(
        data = observed_df,
        aes(x = lineage_passage_index, y = observed_kary_N),
        color = "#d95f02",
        size = 1.25,
        alpha = 0.7,
        position = position_jitter(width = 0.12, height = 0)
      )
  }
  p <- p +
    facet_grid(cohort ~ lineage_label, scales = "free_x", space = "free_x") +
    scale_color_manual(values = c("Across-seed mean" = "#1b9e77"), drop = FALSE) +
    labs(
      title = "Multi-seed predicted chromosome-count quantiles versus observed cells by passage",
      subtitle = "Green curves show across-seed mean predicted quantiles; orange points are observed single-cell chromosome counts.",
      x = "Passage index",
      y = "Chromosome count (N)",
      color = "Predicted fit"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  ggplot2::ggsave(out_path, p, width = 10, height = 6.8)
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
  ggplot2::ggsave(out_path, p, width = 12, height = 6)
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
  ggplot2::ggsave(out_path, p, width = 12, height = 6)
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
    path <- seed_prediction_path(seed_dir, paste0("predict_ploidy_weighted_mean_", horizon_tag, ".tsv"))
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
    path <- seed_prediction_path(seed_dir, paste0("predict_burden_", horizon_tag, ".tsv"))
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

plot_burden_log_seed_mean_by_cohort <- function(seed_day_df, summary_df, out_path, cohort, run_label) {
  if (is.null(seed_day_df) || !nrow(seed_day_df) || is.null(summary_df) || !nrow(summary_df)) {
    return(invisible(NULL))
  }
  curve_df <- seed_day_df[as.character(seed_day_df$cohort) == as.character(cohort), , drop = FALSE]
  curve_df <- curve_df[is.finite(curve_df$day) & is.finite(curve_df$burden_value) & curve_df$burden_value > 0, , drop = FALSE]
  summary_plot_df <- summary_df[as.character(summary_df$cohort) == as.character(cohort), , drop = FALSE]
  summary_plot_df <- summary_plot_df[
    is.finite(summary_plot_df$day) &
      is.finite(summary_plot_df$mean_value) &
      summary_plot_df$mean_value > 0,
    ,
    drop = FALSE
  ]
  if (!nrow(curve_df) || !nrow(summary_plot_df)) return(invisible(NULL))
  curve_df$log10_burden_value <- log10(curve_df$burden_value)
  summary_plot_df$log10_mean_value <- log10(summary_plot_df$mean_value)
  curve_df <- curve_df[order(seed_order_key(curve_df$seed), curve_df$seed, curve_df$day), , drop = FALSE]
  summary_plot_df <- summary_plot_df[order(summary_plot_df$day), , drop = FALSE]
  color <- c("2N" = "#1f77b4", "4N" = "#d62728")[[as.character(cohort)]]
  if (is.null(color) || is.na(color)) color <- "#2c7fb8"

  y_vals <- c(curve_df$log10_burden_value, summary_plot_df$log10_mean_value)
  y_vals <- y_vals[is.finite(y_vals)]
  y_limits <- if (length(y_vals)) {
    rng <- range(y_vals, na.rm = TRUE)
    pad <- if (identical(rng[[1]], rng[[2]])) 0.5 else 0.04 * diff(rng)
    c(rng[[1]] - pad, rng[[2]] + pad)
  } else {
    NULL
  }

  p <- ggplot() +
    geom_line(
      data = curve_df,
      aes(x = day, y = log10_burden_value, group = seed),
      color = "grey72",
      linewidth = 0.18,
      alpha = 0.55
    ) +
    geom_line(
      data = summary_plot_df,
      aes(x = day, y = log10_mean_value),
      color = color,
      linewidth = 1.05,
      linetype = "solid"
    ) +
    coord_cartesian(xlim = range(curve_df$day, na.rm = TRUE), ylim = y_limits) +
    labs(
      title = paste0("1000-day Total Burden Prediction: ", cohort),
      subtitle = paste0("Grey hairlines = individual seed trajectories; colored line = cross-seed mean | run=", run_label),
      x = "Day",
      y = "log10 total tumor burden (mm^3)"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggplot2::ggsave(out_path, p, width = 12, height = 6)
  invisible(out_path)
}

cleanup_stale_burden_component_outputs <- function(out_dir) {
  stale <- c(
    "predict1000_burden_components_seed_day_mean.tsv",
    "predict1000_burden_components_mean_ci.tsv",
    "predict1000_burden_total_mean_ci_2N_4N.pdf",
    "predict1000_burden_total_mean_ci_2N.pdf",
    "predict1000_burden_total_seed_curves_2N.pdf",
    "predict1000_burden_total_mean_ci_4N.pdf",
    "predict1000_burden_total_seed_curves_4N.pdf"
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
      out_path <- file.path(out_dir, paste0("predict1000_burden_total_log_seed_mean_", cohort_file_tag(cohort), ".pdf"))
      res <- plot_burden_log_seed_mean_by_cohort(
        seed_day_df = burden_seed_day,
        summary_df = burden_summary,
        out_path = out_path,
        cohort = cohort,
        run_label = run_label
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
    stop("Usage: Rscript extra_results.R --run_dir=/path/to/run [--out_dir=/path/to/out] [--near_thresh=0.05] [--allow_partial_seed_dirs=FALSE]")
  }
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  out_dir <- normalizePath(args$out_dir %||% file.path(run_dir, "extra_results"), mustWork = FALSE)
  allow_partial_seed_dirs <- as_bool(args$allow_partial_seed_dirs, FALSE)
  near_thresh <- as_num(args$near_thresh, 0.05)
  if (!is.finite(near_thresh) || near_thresh <= 0 || near_thresh >= 0.5) {
    stop("near_thresh must be in (0, 0.5).")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  seed_dirs <- find_seed_dirs(run_dir)
  if (!length(seed_dirs)) {
    stop("No valid seed directories found under: ", run_dir)
  }
  existing_seed_summary_path <- file.path(out_dir, "seed_summary.tsv")
  if (!isTRUE(allow_partial_seed_dirs) && file.exists(existing_seed_summary_path)) {
    existing_seed_summary <- tryCatch(
      utils::read.delim(existing_seed_summary_path, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(existing_seed_summary) &&
        "seed" %in% names(existing_seed_summary) &&
        length(unique(as.character(existing_seed_summary$seed))) > length(seed_dirs)) {
      stop(
        "Refusing to overwrite an existing all-seed extra_results summary with fewer seed directories. ",
        "Existing seed_summary.tsv has ",
        length(unique(as.character(existing_seed_summary$seed))),
        " seeds but only ",
        length(seed_dirs),
        " seed directories were found under run_dir. Re-run with --allow_partial_seed_dirs=TRUE only if this truncation is intentional."
      )
    }
  }

  long_rows <- vector("list", length(seed_dirs))
  summary_records <- vector("list", length(seed_dirs))
  joint_invitro_long_rows <- vector("list", length(seed_dirs))
  joint_soft_rows <- vector("list", length(seed_dirs))
  first_invitro_param_table <- NULL

  for (i in seq_along(seed_dirs)) {
    seed_dir <- seed_dirs[[i]]
    seed <- basename(seed_dir)
    fit_summary_vals <- read_metric_map(file.path(seed_dir, "fit_summary.tsv"), "metric", "value")
    fit_summary_vals <- supplement_joint_invitro_metrics(fit_summary_vals, seed_dir)
    best_vals <- read_metric_map(file.path(seed_dir, "best_params.tsv"), "parameter", "value")
    pred_gate_metrics <- read_1000day_ploidy_gate_metrics(seed_dir)
    param_table_path <- find_parameter_table_for_seed(run_dir, seed_dir, fit_summary_vals)
    param_table <- read_parameter_table_checked(param_table_path)
    if (is.null(first_invitro_param_table) && is_invitro_fit_summary(fit_summary_vals)) {
      first_invitro_param_table <- param_table
    }
    objective <- as_num(summary_metric_value(fit_summary_vals, "objective", summary_metric_value(fit_summary_vals, "objective_total", NA_real_)))
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
    if (is_joint_fit_summary(fit_summary_vals)) {
      soft_path <- file.path(seed_dir, "joint_soft_coupling.tsv")
      if (file.exists(soft_path)) {
        soft_tab <- tryCatch(
          utils::read.delim(soft_path, check.names = FALSE, stringsAsFactors = FALSE),
          error = function(e) NULL
        )
        if (is.data.frame(soft_tab) && nrow(soft_tab) > 0L) {
          soft_tab$seed <- seed
          joint_soft_rows[[i]] <- soft_tab
        }
      }
      invitro_param_table_path <- find_joint_invitro_parameter_table_for_seed(run_dir, seed_dir)
      if (!is.na(invitro_param_table_path) && nzchar(invitro_param_table_path)) {
        invitro_param_table <- read_parameter_table_checked(invitro_param_table_path)
        invitro_fit_summary_vals <- prepare_joint_invitro_fit_summary_values(fit_summary_vals)
        invitro_objective <- as_num(summary_metric_value(invitro_fit_summary_vals, "objective", NA_real_), NA_real_)
        joint_invitro_long_rows[[i]] <- build_parameter_long_table(
          seed = seed,
          objective = invitro_objective,
          fit_summary_vals = invitro_fit_summary_vals,
          best_vals = best_vals,
          param_table = invitro_param_table,
          near_thresh = near_thresh
        )
      }
    }
  }

  parameter_long <- do.call(rbind, long_rows)
  joint_invitro_long_rows <- joint_invitro_long_rows[vapply(joint_invitro_long_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
  joint_invitro_parameter_long <- if (length(joint_invitro_long_rows)) do.call(rbind, joint_invitro_long_rows) else NULL
  joint_soft_rows <- joint_soft_rows[vapply(joint_soft_rows, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
  joint_soft_coupling_all <- if (length(joint_soft_rows)) do.call(rbind, joint_soft_rows) else NULL
  seed_summary <- bind_records(summary_records)
  for (col in c(
    "fit_mode",
    "objective_total",
    "total_loglik",
    "growth_loglik",
    "ploidy_loglik",
    "flow_loglik",
    "growth_loglik_sum",
    "ploidy_loglik_sum",
    "flow_loglik_sum",
    "sigma_growth",
    "sigma_kary",
    "sigma_flow_ploidy",
    "n_growth",
    "n_ploidy_passages",
    "n_kary_cells",
    "n_flow_passages",
    "n_flow_samples",
    "objective_invivo",
    "objective_invitro",
    "objective_soft_coupling",
    "objective_constraints",
    "joint_weight_invivo",
    "joint_weight_invitro",
    "joint_invitro_growth_weight",
    "joint_invitro_ploidy_weight",
    "joint_invitro_flow_weight",
    "joint_soft_coupling_enabled",
    "joint_soft_coupling_params",
    "joint_soft_coupling_n_params",
    "joint_soft_coupling_sigma_default",
    "n_cores_requested",
    "n_cores_used",
    "n_parameters",
    "n_invivo_scenarios",
    "itermax",
    "optimizer_deoptim_objective",
    "optimizer_local_objective",
    "optimizer_local_attempted",
    "optimizer_local_accepted",
    "optimizer_local_convergence",
    "optimizer_local_maxit",
    "optimizer_interrupted",
    "optimizer_iter_completed",
    "optimizer_iter_target",
    "deoptim_stop_reason"
  )) {
    if (!(col %in% names(seed_summary))) seed_summary[[col]] <- NA
  }
  seed_summary <- supplement_optimizer_fields_from_refinement_csv(seed_summary, run_dir)
  seed_summary$fit_mode <- as.character(seed_summary$fit_mode)
  seed_summary$objective <- suppressWarnings(as.numeric(seed_summary$objective))
  seed_summary$objective_total <- suppressWarnings(as.numeric(seed_summary$objective_total))
  seed_summary$total_loglik <- suppressWarnings(as.numeric(seed_summary$total_loglik))
  seed_summary$growth_loglik <- suppressWarnings(as.numeric(seed_summary$growth_loglik))
  seed_summary$ploidy_loglik <- suppressWarnings(as.numeric(seed_summary$ploidy_loglik))
  seed_summary$flow_loglik <- suppressWarnings(as.numeric(seed_summary$flow_loglik))
  seed_summary$growth_loglik_sum <- suppressWarnings(as.numeric(seed_summary$growth_loglik_sum))
  seed_summary$ploidy_loglik_sum <- suppressWarnings(as.numeric(seed_summary$ploidy_loglik_sum))
  seed_summary$flow_loglik_sum <- suppressWarnings(as.numeric(seed_summary$flow_loglik_sum))
  seed_summary$sigma_growth <- suppressWarnings(as.numeric(seed_summary$sigma_growth))
  seed_summary$sigma_kary <- suppressWarnings(as.numeric(seed_summary$sigma_kary))
  seed_summary$sigma_flow_ploidy <- suppressWarnings(as.numeric(seed_summary$sigma_flow_ploidy))
  seed_summary$n_growth <- suppressWarnings(as.numeric(seed_summary$n_growth))
  seed_summary$n_ploidy_passages <- suppressWarnings(as.numeric(seed_summary$n_ploidy_passages))
  seed_summary$n_kary_cells <- suppressWarnings(as.numeric(seed_summary$n_kary_cells))
  seed_summary$n_flow_passages <- suppressWarnings(as.numeric(seed_summary$n_flow_passages))
  seed_summary$n_flow_samples <- suppressWarnings(as.numeric(seed_summary$n_flow_samples))
  seed_summary$objective_invivo <- suppressWarnings(as.numeric(seed_summary$objective_invivo))
  seed_summary$objective_invitro <- suppressWarnings(as.numeric(seed_summary$objective_invitro))
  seed_summary$objective_soft_coupling <- suppressWarnings(as.numeric(seed_summary$objective_soft_coupling))
  seed_summary$objective_constraints <- suppressWarnings(as.numeric(seed_summary$objective_constraints))
  seed_summary$joint_weight_invivo <- suppressWarnings(as.numeric(seed_summary$joint_weight_invivo))
  seed_summary$joint_weight_invitro <- suppressWarnings(as.numeric(seed_summary$joint_weight_invitro))
  seed_summary$joint_invitro_growth_weight <- suppressWarnings(as.numeric(seed_summary$joint_invitro_growth_weight))
  seed_summary$joint_invitro_ploidy_weight <- suppressWarnings(as.numeric(seed_summary$joint_invitro_ploidy_weight))
  seed_summary$joint_invitro_flow_weight <- suppressWarnings(as.numeric(seed_summary$joint_invitro_flow_weight))
  seed_summary$joint_soft_coupling_enabled <- as.logical(seed_summary$joint_soft_coupling_enabled)
  seed_summary$joint_soft_coupling_n_params <- suppressWarnings(as.numeric(seed_summary$joint_soft_coupling_n_params))
  seed_summary$joint_soft_coupling_sigma_default <- suppressWarnings(as.numeric(seed_summary$joint_soft_coupling_sigma_default))
  seed_summary$n_cores_requested <- suppressWarnings(as.numeric(seed_summary$n_cores_requested))
  seed_summary$n_cores_used <- suppressWarnings(as.numeric(seed_summary$n_cores_used))
  seed_summary$n_parameters <- suppressWarnings(as.numeric(seed_summary$n_parameters))
  seed_summary$n_invivo_scenarios <- suppressWarnings(as.numeric(seed_summary$n_invivo_scenarios))
  seed_summary$itermax <- suppressWarnings(as.numeric(seed_summary$itermax))
  seed_summary$optimizer_deoptim_objective <- suppressWarnings(as.numeric(seed_summary$optimizer_deoptim_objective))
  seed_summary$optimizer_local_objective <- suppressWarnings(as.numeric(seed_summary$optimizer_local_objective))
  seed_summary$optimizer_local_attempted <- as.logical(seed_summary$optimizer_local_attempted)
  seed_summary$optimizer_local_accepted <- as.logical(seed_summary$optimizer_local_accepted)
  seed_summary$optimizer_local_convergence <- suppressWarnings(as.numeric(seed_summary$optimizer_local_convergence))
  seed_summary$optimizer_local_maxit <- suppressWarnings(as.numeric(seed_summary$optimizer_local_maxit))
  seed_summary$optimizer_interrupted <- as.character(seed_summary$optimizer_interrupted)
  seed_summary$optimizer_iter_completed <- suppressWarnings(as.numeric(seed_summary$optimizer_iter_completed))
  seed_summary$optimizer_iter_target <- suppressWarnings(as.numeric(seed_summary$optimizer_iter_target))
  seed_summary$deoptim_stop_reason <- as.character(seed_summary$deoptim_stop_reason)
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
  is_joint_run <- any(seed_summary$fit_mode == "fit_joint", na.rm = TRUE) ||
    any(is.finite(seed_summary$objective_invivo) | is.finite(seed_summary$objective_invitro))
  is_invitro_run <- any(seed_summary$fit_mode == "fit_invitro", na.rm = TRUE) ||
    any(is.finite(seed_summary$objective_total) | is.finite(seed_summary$growth_loglik) | is.finite(seed_summary$ploidy_loglik))
  is_invitro_only_run <- isTRUE(is_invitro_run) && !isTRUE(is_joint_run)
  boundary_order <- if (isTRUE(is_joint_run) || isTRUE(is_invitro_run)) {
    order(
      seed_summary$boundary_penalty_active,
      -seed_summary$min_rel_dist_active,
      seed_summary$objective,
      seed_summary$seed,
      na.last = TRUE
    )
  } else {
    order(
      seed_summary$boundary_penalty_active,
      -seed_summary$min_rel_dist_active,
      seed_summary$objective_burden,
      seed_summary$objective,
      seed_summary$seed,
      na.last = TRUE
    )
  }
  seed_summary$boundary_rank_active_support <- NA_integer_
  seed_summary$boundary_rank_active_support[boundary_order] <- seq_len(nrow(seed_summary))
  if (isTRUE(is_joint_run)) {
    joint_objective_order <- function(idx) {
      if (!length(idx)) return(integer(0))
      idx[order(
        seed_summary$objective[idx],
        seed_order_key(seed_summary$seed[idx]),
        seed_summary$seed[idx],
        na.last = TRUE
      )]
    }
    joint_pred_gate <- !is.na(seed_summary$pred1000_both_gt44) & seed_summary$pred1000_both_gt44
    joint_eligible_idx <- which(joint_pred_gate & is.finite(seed_summary$objective))
    if (length(joint_eligible_idx)) {
      joint_ineligible_idx <- setdiff(seq_len(nrow(seed_summary)), joint_eligible_idx)
      recommend_order <- c(
        joint_objective_order(joint_eligible_idx),
        joint_objective_order(joint_ineligible_idx)
      )
    } else {
      recommend_order <- joint_objective_order(seq_len(nrow(seed_summary)))
    }
    seed_summary$recommend_score_burden_ploidy_boundary <- rep(NA_real_, nrow(seed_summary))
    seed_summary$recommend_score_burden_ploidy_boundary[recommend_order] <- seq_along(recommend_order)
  } else if (isTRUE(is_invitro_run)) {
    seed_summary$recommend_score_burden_ploidy_boundary <- with(
      seed_summary,
      objective_rank + 0.1 * boundary_rank_active_support
    )
    recommend_order <- order(
      seed_summary$objective,
      seed_summary$boundary_penalty_active,
      -seed_summary$min_rel_dist_active,
      seed_summary$seed,
      na.last = TRUE
    )
  } else {
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
  }
  seed_summary$recommend_score_ploidy_burden_boundary <- seed_summary$recommend_score_burden_ploidy_boundary
  seed_summary$recommend_score_ploidy_boundary <- seed_summary$recommend_score_burden_ploidy_boundary
  seed_summary$recommend_rank_burden_ploidy_boundary_first <- NA_integer_
  seed_summary$recommend_rank_burden_ploidy_boundary_first[recommend_order] <- seq_len(nrow(seed_summary))
  seed_summary$recommend_rank_ploidy_burden_boundary_first <- NA_integer_
  seed_summary$recommend_rank_ploidy_burden_boundary_first[recommend_order] <- seq_len(nrow(seed_summary))
  seed_summary$recommend_rank_ploidy_boundary_first <- NA_integer_
  seed_summary$recommend_rank_ploidy_boundary_first[recommend_order] <- seq_len(nrow(seed_summary))
  seed_summary$recommend_rank_ploidy_first <- seed_summary$recommend_rank_burden_ploidy_boundary_first
  forest_rank_col <- if ("objective_rank" %in% names(seed_summary)) "objective_rank" else "objective"
  forest_rank_simple <- suppressWarnings(as.integer(seed_summary[[forest_rank_col]]))
  forest_rank_plus_ploidy_simple <- rep(NA_integer_, nrow(seed_summary))
  eligible_plot_idx <- which(!is.na(seed_summary$pred1000_both_gt44) & seed_summary$pred1000_both_gt44)
  if (length(eligible_plot_idx) > 0L) {
    eligible_plot_ord <- eligible_plot_idx[order(
      seed_summary$objective[eligible_plot_idx],
      seed_summary$seed[eligible_plot_idx],
      na.last = TRUE
    )]
    forest_rank_plus_ploidy_simple[eligible_plot_ord] <- seq_along(eligible_plot_ord)
  }
  pred_gate_top3_seeds <- get_top_ranked_seeds(
    summary_df = seed_summary,
    n = 3L,
    rank_col = "objective",
    eligible_mask = seed_summary$pred1000_both_gt44
  )
  forest_top3_seeds <- get_unfiltered_forest_top_seeds(
    summary_df = seed_summary,
    is_joint_run = is_joint_run,
    is_invitro_run = is_invitro_run,
    n = 3L
  )
  seed_summary$forest_plot_rank_simple <- forest_rank_simple
  seed_summary$forest_plot_rank_plus_ploidy_simple <- forest_rank_plus_ploidy_simple
  seed_summary <- seed_summary[order(seed_summary$objective, seed_summary$seed), , drop = FALSE]
  row.names(seed_summary) <- NULL
  convergence_summary <- build_convergence_summary(seed_summary)
  convergence_venn_out <- NULL
  if (isTRUE(is_joint_run)) {
    convergence_venn_out <- plot_convergence_venn(
      seed_summary = seed_summary,
      out_dir = out_dir,
      run_label = basename(run_dir)
    )
  }
  top20_parameter_umap_out <- NULL
  joint_ratio_umap_out <- NULL
  if (isTRUE(is_joint_run)) {
    joint_ratio_umap_out <- plot_joint_ratio_umap(
      run_dir = run_dir,
      out_dir = out_dir,
      run_label = basename(run_dir),
      top_n = 10L
    )
  } else {
    if (isTRUE(is_invitro_only_run)) {
      top20_parameter_umap_out <- plot_top_seed_parameter_umap(
        summary_df = seed_summary,
        parameter_long = parameter_long,
        out_dir = out_dir,
        run_label = basename(run_dir),
        fit_label = "In Vitro",
        filename_prefix = "invitro"
      )
    } else {
      top20_parameter_umap_out <- plot_top_seed_parameter_umap(
        summary_df = seed_summary,
        parameter_long = parameter_long,
        out_dir = out_dir,
        run_label = basename(run_dir),
        fit_label = "In Vivo",
        filename_prefix = "invivo"
      )
    }
  }

  objective_cols <- c("seed", "objective")
  if (isTRUE(is_joint_run)) {
    objective_cols <- c(objective_cols, "objective_invivo", "objective_invitro")
  }
  if (isTRUE(is_invitro_run)) {
    objective_cols <- c(objective_cols, "objective_total", "total_loglik", "growth_loglik", "ploidy_loglik", "flow_loglik")
  }
  objective_cols <- c(objective_cols, intersect(c("objective_burden", "objective_ploidy"), names(seed_summary)))
  objective_simple <- seed_summary[, objective_cols, drop = FALSE]
  objective_simple$objective_rank <- suppressWarnings(as.integer(seed_summary$forest_plot_rank_simple))
  objective_simple$objective_rank_plus_ploidy <- suppressWarnings(as.integer(seed_summary$forest_plot_rank_plus_ploidy_simple))
  objective_simple <- objective_simple[, c("seed", "objective_rank", "objective_rank_plus_ploidy", setdiff(names(objective_simple), c("seed", "objective_rank", "objective_rank_plus_ploidy"))), drop = FALSE]
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
    convergence_summary,
    file = file.path(out_dir, "convergence_summary.tsv"),
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
  if (!is.null(joint_invitro_parameter_long) && nrow(joint_invitro_parameter_long)) {
    utils::write.table(
      joint_invitro_parameter_long,
      file = file.path(out_dir, "invitro_parameter_boundary_long.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  if (!is.null(joint_soft_coupling_all) && nrow(joint_soft_coupling_all)) {
    utils::write.table(
      joint_soft_coupling_all,
      file = file.path(out_dir, "joint_soft_coupling_all.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
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
    near_thresh = near_thresh,
    top3_seeds = forest_top3_seeds
  )
  forest_log_out <- NULL
  if (!isTRUE(is_invitro_only_run)) {
    forest_log_out <- plot_parameter_boundary_forest(
      long_df = parameter_long,
      summary_df = seed_summary,
      out_path = file.path(out_dir, "parameter_boundary_forest_log_x.pdf"),
      run_label = basename(run_dir),
      near_thresh = near_thresh,
      top3_seeds = forest_top3_seeds,
      title_suffix = "Original values on log10 x-axis",
      x_scale = "log10_original"
    )
  }
  objective_risk_out <- plot_objective_vs_boundary_risk(
    summary_df = seed_summary,
    out_path = file.path(out_dir, "objective_vs_boundary_risk.pdf"),
    run_label = basename(run_dir)
  )
  joint_objective_components_out <- NULL
  joint_objective_tradeoff_out <- NULL
  joint_soft_delta_out <- NULL
  joint_soft_penalty_out <- NULL
  joint_soft_pairs_out <- NULL
  joint_soft_context_forest_out <- NULL
  joint_soft_context_forest_log_out <- NULL
  joint_soft_context_forest_filtered_out <- NULL
  joint_soft_context_forest_filtered_log_out <- NULL
  invitro_objective_components_out <- NULL
  invitro_objective_component_distributions_out <- NULL
  invitro_objective_risk_out <- NULL
  invitro_optimization_diagnostics_out <- NULL
  invitro_parameter_positions_log_out <- NULL
  invitro_best_fit_likelihood_out <- character(0)
  invitro_distribution_quantiles_out <- NULL
  unlink(
    file.path(out_dir, c(
      "optimization_diagnostics_objective_draws.tsv",
      "optimization_diagnostics_parameter_long.tsv",
      "optimization_diagnostics_log_x.pdf"
    )),
    force = TRUE
  )
  if (isTRUE(is_joint_run)) {
    joint_objective_components_out <- plot_joint_objective_components(
      summary_df = seed_summary,
      out_path = file.path(out_dir, "joint_objective_components.pdf"),
      run_label = basename(run_dir)
    )
    joint_objective_tradeoff_out <- plot_joint_objective_tradeoff(
      summary_df = seed_summary,
      out_path = file.path(out_dir, "joint_objective_tradeoff.pdf"),
      run_label = basename(run_dir)
    )
    if (!is.null(joint_soft_coupling_all) && nrow(joint_soft_coupling_all)) {
      joint_soft_delta_out <- plot_joint_soft_delta_magnitude(
        soft_df = joint_soft_coupling_all,
        out_path = file.path(out_dir, "joint_soft_coupling_delta_magnitude.pdf"),
        run_label = basename(run_dir)
      )
      joint_soft_penalty_out <- plot_joint_soft_penalty_ranking(
        soft_df = joint_soft_coupling_all,
        out_path = file.path(out_dir, "joint_soft_coupling_penalty_ranking.pdf"),
        run_label = basename(run_dir)
      )
      joint_soft_pairs_out <- plot_joint_soft_vivo_vitro_pairs(
        soft_df = joint_soft_coupling_all,
        out_path = file.path(out_dir, "joint_soft_coupling_vivo_vitro_pairs.pdf"),
        run_label = basename(run_dir)
      )
      joint_soft_context_forest_out <- plot_joint_soft_context_boundary_forest(
        soft_df = joint_soft_coupling_all,
        summary_df = seed_summary,
        out_path = file.path(out_dir, "joint_soft_coupling_context_boundary_forest.pdf"),
        run_label = basename(run_dir),
        near_thresh = near_thresh,
        top3_seeds = forest_top3_seeds
      )
      joint_soft_context_forest_log_out <- plot_joint_soft_context_boundary_forest(
        soft_df = joint_soft_coupling_all,
        summary_df = seed_summary,
        out_path = file.path(out_dir, "joint_soft_coupling_context_boundary_forest_log_x.pdf"),
        run_label = basename(run_dir),
        near_thresh = near_thresh,
        top3_seeds = forest_top3_seeds,
        title_suffix = "Original values on log10 x-axis",
        x_scale = "log10_original"
      )
    }
    joint_cols <- intersect(
      c(
        "seed",
        "objective_rank",
        "objective",
        "objective_invivo",
        "objective_invitro",
        "objective_soft_coupling",
        "objective_constraints",
        "joint_soft_coupling_enabled",
        "joint_soft_coupling_params",
        "joint_weight_invivo",
        "joint_weight_invitro",
        "joint_invitro_growth_weight",
        "joint_invitro_ploidy_weight",
        "joint_invitro_flow_weight",
        "boundary_penalty_active",
        "min_rel_dist_active"
      ),
      names(seed_summary)
    )
    joint_simple <- seed_summary[, joint_cols, drop = FALSE]
    if ("objective_rank" %in% names(joint_simple)) {
      joint_simple$objective_rank <- suppressWarnings(as.integer(joint_simple$objective_rank))
    }
    joint_simple <- joint_simple[order(joint_simple$objective, joint_simple$seed, na.last = TRUE), , drop = FALSE]
    utils::write.table(
      joint_simple,
      file = file.path(out_dir, "joint_objective_simple.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  if (isTRUE(is_invitro_run)) {
    invitro_summary_for_plot <- if (isTRUE(is_joint_run)) {
      prepare_invitro_summary_for_plot(seed_summary, is_joint_run = TRUE)
    } else {
      seed_summary
    }
    invitro_parameter_long_for_plot <- if (isTRUE(is_joint_run) && !is.null(joint_invitro_parameter_long) && nrow(joint_invitro_parameter_long)) {
      joint_invitro_parameter_long
    } else {
      parameter_long
    }
    invitro_optimization_diagnostics_out <- plot_invitro_optimization_diagnostics(
      summary_df = invitro_summary_for_plot,
      parameter_long = invitro_parameter_long_for_plot,
      out_path = file.path(out_dir, "optimization_diagnostics.pdf"),
      run_label = basename(run_dir),
      near_thresh = near_thresh
    )
    invitro_parameter_positions_log_out <- plot_invitro_optimization_diagnostics(
      summary_df = invitro_summary_for_plot,
      parameter_long = invitro_parameter_long_for_plot,
      out_path = file.path(out_dir, "invitro_parameter_positions_log_x.pdf"),
      run_label = basename(run_dir),
      near_thresh = near_thresh,
      parameter_x_scale = "log10_original",
      parameter_panel_only = TRUE
    )
    invitro_objective_components_out <- plot_invitro_objective_components(
      summary_df = invitro_summary_for_plot,
      out_path = file.path(out_dir, "invitro_objective_components.pdf"),
      run_label = basename(run_dir)
    )
    invitro_best_fit_likelihood_out <- plot_invitro_best_fit_likelihood_comparison(
      seed_dirs = seed_dirs,
      seed_summary = invitro_summary_for_plot,
      out_dir = out_dir,
      run_label = basename(run_dir)
    )
    invitro_distribution_quantiles_out <- plot_invitro_distribution_quantiles_multiseed(
      seed_dirs = seed_dirs,
      out_path = file.path(out_dir, "invitro_karyotype_quantiles_multiseed.pdf"),
      out_table = file.path(out_dir, "invitro_karyotype_quantiles_multiseed.tsv"),
      run_label = basename(run_dir)
    )
    if (isTRUE(is_joint_run)) {
      invitro_objective_component_distributions_out <- plot_invitro_objective_component_distributions(
        summary_df = invitro_summary_for_plot,
        out_path = file.path(out_dir, "invitro_objective_components_violin.pdf"),
        out_table = file.path(out_dir, "invitro_objective_components_long.tsv"),
        run_label = basename(run_dir)
      )
      invitro_objective_risk_out <- plot_objective_vs_boundary_risk(
        summary_df = invitro_summary_for_plot,
        out_path = file.path(out_dir, "invitro_objective_vs_boundary_risk.pdf"),
        run_label = paste0(basename(run_dir), " in vitro")
      )
    }
    invitro_cols <- intersect(
      c(
        "seed",
        "objective_rank",
        "objective",
        "objective_total",
        "total_loglik",
        "growth_loglik",
        "ploidy_loglik",
        "flow_loglik",
        "sigma_growth",
        "sigma_kary",
        "sigma_flow_ploidy",
        "n_growth",
        "n_ploidy_passages",
        "n_kary_cells",
        "n_flow_passages",
        "n_flow_samples",
        "boundary_penalty_active",
        "min_rel_dist_active"
      ),
      names(seed_summary)
    )
    invitro_simple <- invitro_summary_for_plot[, invitro_cols, drop = FALSE]
    if ("objective_rank" %in% names(invitro_simple)) {
      invitro_simple$objective_rank <- suppressWarnings(as.integer(invitro_simple$objective_rank))
    }
    invitro_simple <- invitro_simple[order(invitro_simple$objective, invitro_simple$seed, na.last = TRUE), , drop = FALSE]
    utils::write.table(
      invitro_simple,
      file = file.path(out_dir, "invitro_objective_simple.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  if (isTRUE(is_invitro_only_run)) {
    forest_filtered_out <- NULL
    forest_filtered_log_out <- NULL
    prediction_outputs <- character(0)
  } else {
    forest_filtered_out <- plot_parameter_boundary_forest(
      long_df = parameter_long,
      summary_df = seed_summary,
      out_path = file.path(out_dir, "parameter_boundary_forest_pred1000_gt44_top3.pdf"),
      run_label = basename(run_dir),
      near_thresh = near_thresh,
      top3_seeds = pred_gate_top3_seeds,
      title_suffix = "All seeds shown; top 3 with 2N/4N 1000d predictions > 44 highlighted",
      legend_title = "Top 3 Seeds with 2N/4N 1000d > 44"
    )
    forest_filtered_log_out <- plot_parameter_boundary_forest(
      long_df = parameter_long,
      summary_df = seed_summary,
      out_path = file.path(out_dir, "parameter_boundary_forest_pred1000_gt44_top3_log_x.pdf"),
      run_label = basename(run_dir),
      near_thresh = near_thresh,
      top3_seeds = pred_gate_top3_seeds,
      title_suffix = "All seeds shown; top 3 with 2N/4N 1000d predictions > 44 highlighted; original values on log10 x-axis",
      legend_title = "Top 3 Seeds with 2N/4N 1000d > 44",
      x_scale = "log10_original"
    )
    if (isTRUE(is_joint_run) && !is.null(joint_soft_coupling_all) && nrow(joint_soft_coupling_all)) {
      joint_soft_context_forest_filtered_out <- plot_joint_soft_context_boundary_forest(
        soft_df = joint_soft_coupling_all,
        summary_df = seed_summary,
        out_path = file.path(out_dir, "joint_soft_coupling_context_boundary_forest_pred1000_gt44_top3.pdf"),
        run_label = basename(run_dir),
        near_thresh = near_thresh,
        top3_seeds = pred_gate_top3_seeds,
        title_suffix = "All seeds shown; top 3 with 2N/4N 1000d predictions > 44 highlighted",
        legend_title = "Top 3 Seeds with 2N/4N 1000d > 44"
      )
      joint_soft_context_forest_filtered_log_out <- plot_joint_soft_context_boundary_forest(
        soft_df = joint_soft_coupling_all,
        summary_df = seed_summary,
        out_path = file.path(out_dir, "joint_soft_coupling_context_boundary_forest_pred1000_gt44_top3_log_x.pdf"),
        run_label = basename(run_dir),
        near_thresh = near_thresh,
        top3_seeds = pred_gate_top3_seeds,
        title_suffix = "All seeds shown; top 3 with 2N/4N 1000d predictions > 44 highlighted; original values on log10 x-axis",
        legend_title = "Top 3 Seeds with 2N/4N 1000d > 44",
        x_scale = "log10_original"
      )
    }
    prediction_outputs <- plot_prediction_summaries(
      seed_dirs = seed_dirs,
      out_dir = out_dir,
      run_label = basename(run_dir),
      horizon_tag = "0_1000day"
    )
  }

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
  message("Wrote convergence summary table: ", file.path(out_dir, "convergence_summary.tsv"))
  message("Wrote parameter long table: ", file.path(out_dir, "parameter_boundary_long.tsv"))
  message("Wrote objective simple table: ", file.path(out_dir, "seed_objective_simple.tsv"))
  if (!is.null(convergence_venn_out) && file.exists(convergence_venn_out)) {
    message("Wrote convergence Venn diagram: ", convergence_venn_out)
  }
  if (!is.null(top20_parameter_umap_out) && file.exists(top20_parameter_umap_out)) {
    message("Wrote top20 seed parameter UMAP: ", top20_parameter_umap_out)
  }
  if (!is.null(joint_ratio_umap_out) && file.exists(joint_ratio_umap_out)) {
    message("Wrote joint ratio UMAP: ", joint_ratio_umap_out)
  }
  if (!is.null(forest_out) && file.exists(forest_out)) {
    message("Wrote forest plot: ", forest_out)
  } else {
    message("Skipped forest plot because no active fitted parameters were available.")
  }
  if (!is.null(forest_log_out) && file.exists(forest_log_out)) {
    message("Wrote log-x forest plot: ", forest_log_out)
  }
  if (!is.null(forest_filtered_out) && file.exists(forest_filtered_out)) {
    message("Wrote pred-gated top3-highlight forest plot: ", forest_filtered_out)
  } else {
    message("Skipped pred-gated top3-highlight forest plot because no plotted parameters were available.")
  }
  if (!is.null(forest_filtered_log_out) && file.exists(forest_filtered_log_out)) {
    message("Wrote pred-gated top3-highlight log-x forest plot: ", forest_filtered_log_out)
  }
  if (!is.null(objective_risk_out) && file.exists(objective_risk_out)) {
    message("Wrote objective-risk plot: ", objective_risk_out)
  } else {
    message("Skipped objective-risk plot because no seeds had a finite positive boundary-distance metric.")
  }
  if (isTRUE(is_joint_run)) {
    if (!is.null(joint_objective_components_out) && file.exists(joint_objective_components_out)) {
      message("Wrote joint objective-components plot: ", joint_objective_components_out)
    } else {
      message("Skipped joint objective-components plot because no finite joint objective fields were available.")
    }
    if (!is.null(joint_objective_tradeoff_out) && file.exists(joint_objective_tradeoff_out)) {
      message("Wrote joint objective tradeoff plot: ", joint_objective_tradeoff_out)
    } else {
      message("Skipped joint objective tradeoff plot because no finite in vivo/in vitro objective pair was available.")
    }
    if (!is.null(joint_soft_coupling_all) && nrow(joint_soft_coupling_all)) {
      message("Wrote joint soft-coupling table: ", file.path(out_dir, "joint_soft_coupling_all.tsv"))
      if (!is.null(joint_soft_delta_out) && file.exists(joint_soft_delta_out)) {
        message("Wrote joint soft-coupling delta plot: ", joint_soft_delta_out)
      }
      if (!is.null(joint_soft_penalty_out) && file.exists(joint_soft_penalty_out)) {
        message("Wrote joint soft-coupling penalty plot: ", joint_soft_penalty_out)
      }
      if (!is.null(joint_soft_pairs_out) && file.exists(joint_soft_pairs_out)) {
        message("Wrote joint soft-coupling vivo/vitro pair plot: ", joint_soft_pairs_out)
      }
      if (!is.null(joint_soft_context_forest_out) && file.exists(joint_soft_context_forest_out)) {
        message("Wrote joint soft-coupling context boundary forest plot: ", joint_soft_context_forest_out)
      }
      if (!is.null(joint_soft_context_forest_filtered_out) && file.exists(joint_soft_context_forest_filtered_out)) {
        message("Wrote pred-gated joint soft-coupling context boundary forest plot: ", joint_soft_context_forest_filtered_out)
      }
      if (!is.null(joint_soft_context_forest_log_out) && file.exists(joint_soft_context_forest_log_out)) {
        message("Wrote joint soft-coupling context boundary forest log-x plot: ", joint_soft_context_forest_log_out)
      }
      if (!is.null(joint_soft_context_forest_filtered_log_out) && file.exists(joint_soft_context_forest_filtered_log_out)) {
        message("Wrote pred-gated joint soft-coupling context boundary forest log-x plot: ", joint_soft_context_forest_filtered_log_out)
      }
    }
    message("Wrote joint objective simple table: ", file.path(out_dir, "joint_objective_simple.tsv"))
  }
  if (isTRUE(is_invitro_run)) {
    if (!is.null(invitro_optimization_diagnostics_out) && file.exists(invitro_optimization_diagnostics_out)) {
      message("Wrote in vitro optimization-diagnostics plot: ", invitro_optimization_diagnostics_out)
    } else {
      message("Skipped in vitro optimization-diagnostics plot because no finite in vitro objective fields were available.")
    }
    if (!is.null(invitro_parameter_positions_log_out) && file.exists(invitro_parameter_positions_log_out)) {
      message("Wrote in vitro log-x parameter-position plot: ", invitro_parameter_positions_log_out)
    }
    if (!is.null(invitro_objective_components_out) && file.exists(invitro_objective_components_out)) {
      message("Wrote in vitro objective-components plot: ", invitro_objective_components_out)
    } else {
      message("Skipped in vitro objective-components plot because no finite in vitro objective fields were available.")
    }
    if (!is.null(invitro_objective_component_distributions_out) && file.exists(invitro_objective_component_distributions_out)) {
      message("Wrote in vitro objective-component distributions plot: ", invitro_objective_component_distributions_out)
    }
    if (!is.null(invitro_objective_risk_out) && file.exists(invitro_objective_risk_out)) {
      message("Wrote in vitro objective-risk plot: ", invitro_objective_risk_out)
    }
    if (length(invitro_best_fit_likelihood_out)) {
      message("Wrote in vitro best-fit likelihood comparison plots: ", paste(basename(invitro_best_fit_likelihood_out), collapse = ", "))
    } else {
      message("Skipped in vitro best-fit likelihood comparison plots because no selected seeds had passage-level likelihood tables.")
    }
    if (!is.null(invitro_distribution_quantiles_out) && file.exists(invitro_distribution_quantiles_out)) {
      message("Wrote in vitro multi-seed chromosome-count quantiles: ", invitro_distribution_quantiles_out)
    } else {
      message("Skipped in vitro multi-seed chromosome-count quantiles because no distribution quantile tables were available.")
    }
    message("Wrote in vitro objective simple table: ", file.path(out_dir, "invitro_objective_simple.tsv"))
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
