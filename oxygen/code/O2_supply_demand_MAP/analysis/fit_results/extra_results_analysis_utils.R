#!/usr/bin/env Rscript

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
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_invitro_parameter_utils.R"), local = environment())
rm(.o2sd_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_bool <- o2sd_as_bool

summary_flag_true <- o2sd_flag_true
summary_flag_na <- o2sd_flag_na

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
    flow_loglik = "invitro_flow_loglik",
    death_loglik = "invitro_death_loglik"
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
    "joint_soft_coupling_welsch_c",
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

bool_from_table_value <- o2sd_flag_true

invitro_parameter_transform_map <- o2sd_invitro_parameter_transform_map

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
    death_loglik = as_num(summary_metric_value(fit_summary_vals, "death_loglik", summary_metric_value(fit_summary_vals, "invitro_death_loglik", NA_real_))),
    growth_loglik_sum = as_num(summary_metric_value(fit_summary_vals, "growth_loglik_sum", NA_real_)),
    ploidy_loglik_sum = as_num(summary_metric_value(fit_summary_vals, "ploidy_loglik_sum", NA_real_)),
    flow_loglik_sum = as_num(summary_metric_value(fit_summary_vals, "flow_loglik_sum", NA_real_)),
    death_loglik_sum = as_num(summary_metric_value(fit_summary_vals, "death_loglik_sum", summary_metric_value(fit_summary_vals, "invitro_death_loglik_sum", NA_real_))),
    sigma_growth = as_num(summary_metric_value(fit_summary_vals, "sigma_growth", NA_real_)),
    sigma_kary = as_num(summary_metric_value(fit_summary_vals, "sigma_kary", NA_real_)),
    sigma_flow_ploidy = as_num(summary_metric_value(fit_summary_vals, "sigma_flow_ploidy", NA_real_)),
    sigma_death_logit = as_num(summary_metric_value(fit_summary_vals, "sigma_death_logit", summary_metric_value(fit_summary_vals, "invitro_sigma_death_logit", NA_real_))),
    death_fraction_eps = as_num(summary_metric_value(fit_summary_vals, "death_fraction_eps", summary_metric_value(fit_summary_vals, "invitro_death_fraction_eps", NA_real_))),
    n_growth = as_num(summary_metric_value(fit_summary_vals, "n_growth", NA_real_)),
    n_ploidy_passages = as_num(summary_metric_value(fit_summary_vals, "n_ploidy_passages", NA_real_)),
    n_kary_cells = as_num(summary_metric_value(fit_summary_vals, "n_kary_cells", NA_real_)),
    n_flow_passages = as_num(summary_metric_value(fit_summary_vals, "n_flow_passages", NA_real_)),
    n_flow_samples = as_num(summary_metric_value(fit_summary_vals, "n_flow_samples", NA_real_)),
    n_death_passages = as_num(summary_metric_value(fit_summary_vals, "n_death_passages", summary_metric_value(fit_summary_vals, "invitro_n_death_passages", NA_real_))),
    joint_weight_invivo = as_num(summary_metric_value(fit_summary_vals, "joint_weight_invivo", NA_real_)),
    joint_weight_invitro = as_num(summary_metric_value(fit_summary_vals, "joint_weight_invitro", NA_real_)),
    joint_invitro_growth_weight = as_num(summary_metric_value(fit_summary_vals, "joint_invitro_growth_weight", NA_real_)),
    joint_invitro_ploidy_weight = as_num(summary_metric_value(fit_summary_vals, "joint_invitro_ploidy_weight", NA_real_)),
    joint_invitro_flow_weight = as_num(summary_metric_value(fit_summary_vals, "joint_invitro_flow_weight", NA_real_)),
    joint_invitro_death_weight = as_num(summary_metric_value(fit_summary_vals, "joint_invitro_death_weight", summary_metric_value(fit_summary_vals, "death_weight", NA_real_))),
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
    joint_soft_coupling_welsch_c = as_num(summary_metric_value(fit_summary_vals, "joint_soft_coupling_welsch_c", NA_real_)),
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

is_missing_summary_value <- o2sd_missing_summary_value

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

build_long_objective_table <- function(seed_summary) {
  fit_mode <- if ("fit_mode" %in% names(seed_summary)) unique(as.character(seed_summary$fit_mode)) else character(0)
  is_invitro <- any(fit_mode == "fit_invitro", na.rm = TRUE) ||
    ("objective_total" %in% names(seed_summary) && any(is.finite(suppressWarnings(as.numeric(seed_summary$objective_total)))))
  is_joint <- any(fit_mode == "fit_joint", na.rm = TRUE) ||
    ("objective_invivo" %in% names(seed_summary) && any(is.finite(suppressWarnings(as.numeric(seed_summary$objective_invivo)))))
  if (is_invitro) {
    seed_summary$objective_total_plot <- suppressWarnings(as.numeric(seed_summary$objective_total))
    if (!any(is.finite(seed_summary$objective_total_plot))) seed_summary$objective_total_plot <- suppressWarnings(as.numeric(seed_summary$objective))
    seed_summary$growth_negloglik <- -suppressWarnings(as.numeric(seed_summary$growth_loglik))
    seed_summary$ploidy_negloglik <- -suppressWarnings(as.numeric(seed_summary$ploidy_loglik))
    seed_summary$flow_negloglik <- -suppressWarnings(as.numeric(seed_summary$flow_loglik))
    seed_summary$death_negloglik <- -suppressWarnings(as.numeric(seed_summary$death_loglik))
    metric_cols <- c("objective_total_plot", "growth_negloglik", "ploidy_negloglik", "flow_negloglik", "death_negloglik")
    metric_labels <- c("objective_total", "growth_-logLik", "ploidy_-logLik", "flow_-logLik", "death_-logLik")
  } else if (is_joint) {
    metric_cols <- c("objective", "objective_invivo", "objective_invitro")
    metric_labels <- metric_cols
  } else {
    metric_cols <- c("objective", "objective_burden", "objective_ploidy")
    metric_labels <- metric_cols
  }
  missing <- setdiff(metric_cols, names(seed_summary))
  if (length(missing)) stop("seed summary is missing required objective columns: ", paste(missing, collapse = ", "))
  rows <- lapply(metric_cols, function(metric) data.frame(seed = as.character(seed_summary$seed), metric = metric, value = suppressWarnings(as.numeric(seed_summary[[metric]])), stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  out <- out[is.finite(out$value), , drop = FALSE]
  out$metric <- factor(out$metric, levels = metric_cols)
  out$metric_label <- factor(out$metric, levels = metric_cols, labels = metric_labels)
  out
}
