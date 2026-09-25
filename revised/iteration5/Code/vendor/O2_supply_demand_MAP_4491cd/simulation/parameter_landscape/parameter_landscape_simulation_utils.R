#!/usr/bin/env Rscript

.o2pl_simulation_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "parameter_landscape_simulation_utils.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.o2pl_workflow_root <- normalizePath(file.path(.o2pl_simulation_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2pl_workflow_root, "util", "o2_supply_demand_map_parameter_landscape_io_utils.R"), local = environment(), chdir = TRUE)
source(file.path(.o2pl_workflow_root, "util", "o2_supply_demand_map_invitro_parameter_utils.R"), local = environment(), chdir = TRUE)
source(file.path(.o2pl_workflow_root, "util", "o2_supply_demand_map_shared.R"), local = environment(), chdir = TRUE)

metric_map <- function(path) {
  tab <- read_tsv(path)
  o2pl_require_columns(tab, c("metric", "value"), path)
  stats::setNames(as.character(tab$value), as.character(tab$metric))
}

metric_value <- function(metrics, keys, default = NA_character_) {
  hit <- keys[keys %in% names(metrics)]
  if (!length(hit)) return(default)
  value <- metrics[[hit[[1L]]]]
  if (!length(value) || is.na(value) || !nzchar(trimws(as.character(value)))) default else value
}

metric_numeric <- function(metrics, keys, seed, label) {
  value <- suppressWarnings(as.numeric(metric_value(metrics, keys, NA_character_)))
  if (!is.finite(value)) {
    warning("Missing/non-finite ", label, " in seed", seed, call. = FALSE)
    return(NA_real_)
  }
  value
}

read_objective_from_summary <- function(seed_dir) {
  path <- file.path(seed_dir, "fit_summary.tsv")
  if (!file.exists(path)) return(NA_real_)
  metric_numeric(
    metric_map(path),
    c("objective", "objective_total", "optimizer_local_objective", "optimizer_deoptim_objective"),
    seed_from_dir(seed_dir),
    "objective"
  )
}

attach_objective <- function(fitted_df, objective_seed_dir) {
  if ("objective" %in% names(fitted_df)) {
    fitted_df$objective <- suppressWarnings(as.numeric(fitted_df$objective))
    return(fitted_df)
  }
  if (is.null(objective_seed_dir) || !dir.exists(objective_seed_dir)) {
    stop("Materialized fitted-parameter table has no objective column and objective_seed_dir is unavailable: ", objective_seed_dir, call. = FALSE)
  }
  o2pl_require_columns(fitted_df, "seed", "materialized fitted-parameter table")
  fitted_df$objective <- vapply(
    fitted_df$seed,
    function(seed) read_objective_from_summary(file.path(objective_seed_dir, paste0("seed", as.integer(seed)))),
    numeric(1)
  )
  if (any(!is.finite(fitted_df$objective))) {
    stop("Missing/non-finite objective values for seeds: ", paste(head(fitted_df$seed[!is.finite(fitted_df$objective)], 20L), collapse = ", "), call. = FALSE)
  }
  fitted_df
}

seed_from_dir <- o2sd_seed_from_dir
list_seed_dirs <- function(input_dir) {
  dirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl("^seed[0-9]+$", basename(dirs))]
  dirs[order(vapply(dirs, seed_from_dir, integer(1)))]
}

read_best_params <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  if (!file.exists(path)) stop("Missing best_params.tsv: ", path, call. = FALSE)
  tab <- read_tsv(path)
  o2pl_require_columns(tab, c("parameter", "value"), path)
  stats::setNames(suppressWarnings(as.numeric(tab$value)), as.character(tab$parameter))
}

invitro_parameter_transform_map <- o2sd_invitro_parameter_transform_map

transform_param_slot <- function(value, transform, param_symbol, slot_label) {
  value <- suppressWarnings(as.numeric(value))
  if (!is.finite(value)) stop("Non-finite ", slot_label, " for parameter ", param_symbol, call. = FALSE)
  if (identical(transform, "log10")) {
    if (value <= 0) stop(param_symbol, " ", slot_label, " must be > 0 for log10 transform.", call. = FALSE)
    return(log10(value))
  }
  value
}

convert_invitro_parameter_table <- function(tab, path) {
  o2pl_require_columns(tab, c("param_symbol", "init_value", "lower_bound", "upper_bound"), path)
  estimate_col <- if ("use_invitro_fit" %in% names(tab)) "use_invitro_fit" else "estimate"
  if (!estimate_col %in% names(tab)) stop("Missing use_invitro_fit/estimate in ", path, call. = FALSE)
  tab$param_symbol <- trimws(as.character(tab$param_symbol))
  map <- invitro_parameter_transform_map()
  tab <- merge(map, tab, by = "param_symbol", all = FALSE, sort = FALSE)
  tab <- tab[match(map$param_symbol, tab$param_symbol, nomatch = 0L), , drop = FALSE]
  data.frame(
    param_name = tab$param_name,
    estimate = vapply(tab[[estimate_col]], as_bool, logical(1), default = FALSE),
    init_value = mapply(transform_param_slot, tab$init_value, tab$transform, tab$param_symbol, "init_value"),
    lower_bound = mapply(transform_param_slot, tab$lower_bound, tab$transform, tab$param_symbol, "lower_bound"),
    upper_bound = mapply(transform_param_slot, tab$upper_bound, tab$transform, tab$param_symbol, "upper_bound"),
    param_prototype = tab$param_symbol,
    stringsAsFactors = FALSE
  )
}

read_param_table <- function(seed_dir, dataset = "invivo") {
  candidates <- c(
    file.path(seed_dir, "parameter_table.csv"), file.path(dirname(seed_dir), "parameter_table.csv"),
    file.path(seed_dir, "parameter_table_input.csv"), file.path(dirname(seed_dir), "parameter_table_input.csv"),
    file.path(seed_dir, "parameter_table_input_invitro.csv"), file.path(dirname(seed_dir), "parameter_table_input_invitro.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) stop("Missing parameter table for ", seed_dir, call. = FALSE)
  path <- existing[[1L]]
  tab <- read_csv_plain(path)
  required <- c("param_name", "estimate", "init_value", "lower_bound", "upper_bound", "param_prototype")
  if (!all(required %in% names(tab))) {
    if (identical(dataset, "invitro") && "param_symbol" %in% names(tab)) return(convert_invitro_parameter_table(tab, path))
    o2pl_require_columns(tab, required, path)
  }
  tab$estimate <- vapply(tab$estimate, as_bool, logical(1), default = FALSE)
  for (name in c("init_value", "lower_bound", "upper_bound")) tab[[name]] <- suppressWarnings(as.numeric(tab[[name]]))
  tab
}

read_np <- function(metrics, param_table) {
  np <- as_int(metric_value(metrics, c("NP", "NP_used", "NP_requested")), NA_integer_)
  if (!is.finite(np) || np < 1L) np <- max(10L * sum(param_table$estimate), 1L)
  np
}

sample_uniform_box <- function(n, lower, upper) {
  n <- as.integer(n)
  if (n <= 0L || !length(lower)) return(matrix(numeric(), nrow = 0L, ncol = length(lower)))
  out <- sweep(matrix(stats::runif(n * length(lower)), nrow = n), 2L, upper - lower, `*`)
  out <- sweep(out, 2L, lower, `+`)
  colnames(out) <- names(lower)
  out
}

sample_truncnorm_box <- function(n, center, lower, upper, sigma_frac = 0.1) {
  n <- as.integer(n)
  if (n <= 0L || !length(lower)) return(matrix(numeric(), nrow = 0L, ncol = length(lower)))
  sd_vec <- pmax((upper - lower) * as.numeric(sigma_frac), 1e-12)
  out <- sweep(matrix(stats::rnorm(n * length(lower)), nrow = n), 2L, sd_vec, `*`)
  out <- sweep(out, 2L, center, `+`)
  out <- sweep(out, 2L, lower, pmax)
  out <- sweep(out, 2L, upper, pmin)
  colnames(out) <- names(lower)
  out
}

build_de_initialpop <- function(np, lower, upper, init_use, mode = "uniform", uniform_frac = 0.3, sigma_frac = 0.1) {
  np <- as.integer(np)
  mode <- tolower(trimws(mode))
  if (np < 1L || !mode %in% c("uniform", "hybrid")) stop("Invalid DE initial-population settings.", call. = FALSE)
  pop <- sample_uniform_box(np, lower, upper)
  init_vec <- pmin(pmax(as.numeric(init_use[names(lower)]), lower), upper)
  pop[1L, ] <- init_vec
  remaining <- np - 1L
  if (remaining <= 0L || mode == "uniform") {
    if (remaining > 0L) pop[2L:np, ] <- sample_uniform_box(remaining, lower, upper)
    return(pop)
  }
  n_uniform <- max(0L, min(remaining, as.integer(round(remaining * uniform_frac))))
  n_local <- remaining - n_uniform
  cursor <- 2L
  if (n_local > 0L) {
    pop[cursor:(cursor + n_local - 1L), ] <- sample_truncnorm_box(n_local, init_vec, lower, upper, sigma_frac)
    cursor <- cursor + n_local
  }
  if (n_uniform > 0L) pop[cursor:(cursor + n_uniform - 1L), ] <- sample_uniform_box(n_uniform, lower, upper)
  if (remaining > 1L) pop[2L:np, ] <- pop[1L + sample.int(remaining), , drop = FALSE]
  pop
}

inverse_transform_column <- function(param_name, x) if (startsWith(param_name, "log10_")) 10^as.numeric(x) else as.numeric(x)

initial_population_natural <- function(seed_dir, seed, metrics, best_names, dataset = "invivo") {
  param_table <- read_param_table(seed_dir, dataset)
  opt <- param_table[param_table$estimate, , drop = FALSE]
  if (!nrow(opt)) stop("No estimated parameters for ", seed_dir, call. = FALSE)
  param_names <- as.character(opt$param_name)
  init <- stats::setNames(as.numeric(opt$init_value), param_names)
  lower <- stats::setNames(as.numeric(opt$lower_bound), param_names)
  upper <- stats::setNames(as.numeric(opt$upper_bound), param_names)
  set.seed(as.integer(seed))
  transformed <- if (dataset == "invitro") sample_uniform_box(read_np(metrics, param_table), lower, upper) else build_de_initialpop(read_np(metrics, param_table), lower, upper, init, "uniform")
  out <- data.frame(matrix(nrow = nrow(transformed), ncol = 0L))
  for (j in seq_along(param_names)) out[[as.character(opt$param_prototype[[j]])]] <- inverse_transform_column(param_names[[j]], transformed[, j])
  if ("c_vol_2N_eff_mm3" %in% best_names && "rho_2N" %in% names(out)) out$c_vol_2N_eff_mm3 <- 1 / out$rho_2N
  best <- read_best_params(seed_dir)
  for (name in setdiff(best_names, names(out))) out[[name]] <- if (name %in% names(best)) as.numeric(best[[name]]) else NA_real_
  out <- out[, best_names, drop = FALSE]
  out$seed <- as.integer(seed)
  out
}

seed_prediction_path <- function(seed_dir, filename) {
  candidates <- c(file.path(seed_dir, "viz", filename), file.path(seed_dir, "viz", "invivo", filename))
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1L]] else candidates[[1L]]
}
prediction_value_column <- function(tab) {
  hit <- c("weighted_mean_N", "weighted_mean_endpoint", "weighted_mean_ploidy")
  hit <- hit[hit %in% names(tab)]
  if (length(hit)) hit[[1L]] else NA_character_
}
prediction_initial_ploidy <- function(tab, value_col, cohort) {
  values <- suppressWarnings(as.numeric(tab[[value_col]]))
  start_with <- if ("start_with" %in% names(tab)) tolower(as.character(tab$start_with)) else ""
  chr_scale <- value_col == "weighted_mean_N" || any(grepl("chr", start_with), na.rm = TRUE) || any(values > 10, na.rm = TRUE)
  if (chr_scale) if (cohort == "4N") 88 else 44 else if (cohort == "4N") 4 else 2
}

read_pred1000_ploidy_ratios <- function(seed_dir, seed, target_day = 1000) {
  out <- data.frame(pred1000_ploidy_ratio_2N = NA_real_, pred1000_ploidy_ratio_4N = NA_real_)
  path <- seed_prediction_path(seed_dir, "predict_ploidy_weighted_mean_0_1000day.tsv")
  if (!file.exists(path)) return(out)
  tab <- tryCatch(read_tsv(path), error = function(e) NULL)
  if (is.null(tab) || !nrow(tab) || !all(c("cohort", "day") %in% names(tab))) return(out)
  value_col <- prediction_value_column(tab)
  if (is.na(value_col)) return(out)
  for (cohort in c("2N", "4N")) {
    rows <- tab[as.character(tab$cohort) == cohort & is.finite(suppressWarnings(as.numeric(tab$day))), , drop = FALSE]
    if (!nrow(rows)) next
    idx <- which.min(abs(as.numeric(rows$day) - target_day))
    value <- suppressWarnings(as.numeric(rows[[value_col]][[idx]]))
    initial <- prediction_initial_ploidy(tab, value_col, cohort)
    out[[paste0("pred1000_ploidy_ratio_", cohort)]] <- value / initial
  }
  out
}

read_pred1000_both_gt44 <- function(seed_dir, seed, target_day = 1000) {
  ratios <- read_pred1000_ploidy_ratios(seed_dir, seed, target_day)
  values <- suppressWarnings(as.numeric(unlist(ratios[1, ], use.names = FALSE)))
  if (any(!is.finite(values))) return(NA)
  ratios$pred1000_ploidy_ratio_2N > 1 && ratios$pred1000_ploidy_ratio_4N > 0.5
}

best_params_row <- function(seed_dir, seed, best_names, dataset = "invivo") {
  values <- read_best_params(seed_dir)
  out <- as.data.frame(as.list(stats::setNames(rep(NA_real_, length(best_names)), best_names)), check.names = FALSE)
  for (name in intersect(best_names, names(values))) out[[name]] <- as.numeric(values[[name]])
  metrics <- metric_map(file.path(seed_dir, "fit_summary.tsv"))
  out$objective <- metric_numeric(metrics, c("objective", "objective_total", "optimizer_local_objective", "optimizer_deoptim_objective"), seed, "objective")
  if (dataset == "invivo") {
    out$objective_ploidy <- metric_numeric(metrics, "objective_ploidy", seed, "objective_ploidy")
    out$objective_burden <- metric_numeric(metrics, "objective_burden", seed, "objective_burden")
  } else {
    out$objective_growth <- -metric_numeric(metrics, "growth_loglik", seed, "growth_loglik")
    out$objective_ploidy <- -metric_numeric(metrics, "ploidy_loglik", seed, "ploidy_loglik")
    out$objective_flow <- -metric_numeric(metrics, "flow_loglik", seed, "flow_loglik")
  }
  out$deoptim_iter_completed <- as_int(metric_value(metrics, c("optimizer_iter_completed", "deoptim_iter_completed", "deoptim_stop_iteration", "iter_completed", "itermax")), NA_integer_)
  if (dataset == "invivo") {
    out$pred1000_both_gt44 <- read_pred1000_both_gt44(seed_dir, seed)
    out <- cbind(out, read_pred1000_ploidy_ratios(seed_dir, seed))
  }
  out$seed <- as.integer(seed)
  out
}

o2pl_materialize_seed_parameter_tables <- function(dataset = "invivo",
                                                    input_dir = default_dataset_input_dir(dataset),
                                                    root_dir = o2pl_default_result_root(),
                                                    max_seeds = NA_integer_,
                                                    write_best = TRUE,
                                                    write_initial = TRUE) {
  dataset <- o2pl_normalize_dataset(dataset)
  input_dir <- normalizePath(path.expand(input_dir), mustWork = TRUE)
  root_dir <- normalizePath(path.expand(root_dir), mustWork = FALSE)
  seed_dirs <- list_seed_dirs(input_dir)
  if (is.finite(max_seeds) && max_seeds > 0L) seed_dirs <- head(seed_dirs, max_seeds)
  if (!length(seed_dirs)) stop("No seed directories under ", input_dir, call. = FALSE)
  best_vectors <- lapply(seed_dirs, read_best_params)
  best_names <- unique(unlist(lapply(best_vectors, names), use.names = FALSE))
  best_rows <- if (write_best) vector("list", length(seed_dirs)) else NULL
  initial_rows <- if (write_initial) vector("list", length(seed_dirs)) else NULL
  for (i in seq_along(seed_dirs)) {
    seed <- seed_from_dir(seed_dirs[[i]])
    metrics <- metric_map(file.path(seed_dirs[[i]], "fit_summary.tsv"))
    if (write_best) best_rows[[i]] <- best_params_row(seed_dirs[[i]], seed, best_names, dataset)
    if (write_initial) initial_rows[[i]] <- initial_population_natural(seed_dirs[[i]], seed, metrics, best_names, dataset)
  }
  canonical_dir <- o2pl_simulation_tables_dir(root_dir, dataset)
  dir.create(canonical_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- list()
  if (write_initial) {
    data <- do.call(rbind, initial_rows)
    data <- data[, c(best_names, "seed"), drop = FALSE]
    legacy <- paper_initial_population_csv(dataset, root_dir)
    canonical <- file.path(canonical_dir, "deoptim_initial_population.tsv")
    schema <- file.path(canonical_dir, "deoptim_initial_population.schema.tsv")
    write_csv(data, legacy)
    write_tsv_plain(data, canonical)
    o2pl_write_schema(data, paste0(dataset, "_deoptim_initial_population"), schema)
    o2pl_record_artifact(root_dir, "simulation", paste0(dataset, "_deoptim_initial_population"), canonical, data, schema, "generate_parameter_landscape_simulation_tables.R", dataset, input_dir)
    outputs$initial_population <- legacy
    outputs$initial_population_materialized <- canonical
  }
  if (write_best) {
    data <- do.call(rbind, best_rows)
    extra <- setdiff(names(data), best_names)
    data <- data[, c(best_names, extra), drop = FALSE]
    legacy <- paper_best_params_csv(dataset, root_dir)
    canonical <- file.path(canonical_dir, "fitted_parameter_features.tsv")
    schema <- file.path(canonical_dir, "fitted_parameter_features.schema.tsv")
    write_csv(data, legacy)
    write_tsv_plain(data, canonical)
    o2pl_write_schema(data, paste0(dataset, "_fitted_parameter_features"), schema)
    o2pl_record_artifact(root_dir, "simulation", paste0(dataset, "_fitted_parameter_features"), canonical, data, schema, "generate_parameter_landscape_simulation_tables.R", dataset, input_dir)
    outputs$best_params <- legacy
    outputs$best_params_materialized <- canonical
  }
  invisible(outputs)
}

o2pl_fixed_o2_environment <- local({
  cache <- NULL
  workflow_root <- .o2pl_workflow_root
  function() {
    if (!is.null(cache)) return(cache)
    script <- file.path(workflow_root, "simulation", "o2", "fixed_o2", "run_fixed_o2_simulation.R")
    env <- new.env(parent = globalenv())
    env$commandArgs <- function(trailingOnly = FALSE) character()
    sys.source(script, envir = env, chdir = TRUE)
    required <- c("generate_fixo2_attractor_mode_table", "fixo2_attractor_mode_summary_by_seed", "fixo2_reference_mode_table")
    missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = TRUE)]
    if (length(missing)) stop("Fixed-O2 simulator is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    cache <<- env
    cache
  }
})

o2pl_materialize_fixed_o2_modes <- function(input_dir,
                                            root_dir,
                                            best_csv = paper_best_params_csv("invivo", root_dir),
                                            attractor_o2_grid = paper_default_attractor_o2_grid(),
                                            summary_o2 = paper_default_mode_summary_o2(),
                                            reference_o2 = 2,
                                            max_seeds = NA_integer_,
                                            n_workers = 1L,
                                            overwrite = FALSE) {
  mode_dir <- paper_fixo2_mode_tables_dir(root_dir)
  common_path <- file.path(mode_dir, "fixed_o2_attractor_mode_by_seed_o2.tsv")
  summary_path <- file.path(mode_dir, "fixed_o2_attractor_mode_summary_by_seed.tsv")
  reference_o2 <- sort(unique(as_num_vec(reference_o2, 2)))
  attractor_o2_grid <- sort(unique(c(as_num_vec(attractor_o2_grid, paper_default_attractor_o2_grid()), reference_o2)))
  if (!overwrite && file.exists(common_path) && file.exists(summary_path)) return(invisible(list(mode_by_seed_o2 = common_path, mode_summary = summary_path)))
  seed_dirs <- list_seed_dirs(input_dir)
  if (is.finite(max_seeds) && max_seeds > 0L) seed_dirs <- head(seed_dirs, max_seeds)
  seed_ids <- basename(seed_dirs)
  env <- o2pl_fixed_o2_environment()
  mode <- get("generate_fixo2_attractor_mode_table", env)(input_dir, attractor_o2_grid, seed_ids, as_int(n_workers, 1L))
  mode$is_mode_reference_o2 <- vapply(as.numeric(mode$O2_pct), function(x) any(abs(x - reference_o2) < 1e-9), logical(1))
  summary <- get("fixo2_attractor_mode_summary_by_seed", env)(mode, as_num_vec(summary_o2, paper_default_mode_summary_o2()))
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv_plain(mode, common_path)
  write_tsv_plain(summary, summary_path)
  best <- if (file.exists(best_csv)) read_csv_plain(best_csv) else data.frame()
  for (o2 in reference_o2) {
    ref <- get("fixo2_reference_mode_table", env)(mode, o2)
    ref_dir <- fixed_o2_reference_dir(mode_dir, o2)
    dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv_plain(ref, file.path(ref_dir, "fixed_o2_mode_reference_by_seed.tsv"))
    write_tsv_plain(ref, file.path(ref_dir, "fixed_o2_mode_by_seed.tsv"))
    write_tsv_plain(summary, file.path(ref_dir, "fixed_o2_mode_summary_by_seed.tsv"))
    if (nrow(best) && "seed" %in% names(best)) {
      ref$seed <- suppressWarnings(as.integer(sub("^seed", "", ref$seed_id)))
      joined <- merge(best, ref, by = "seed", all.x = TRUE, sort = FALSE)
      write_csv(joined, file.path(ref_dir, "invivo_best_params_by_seed_with_fixed_o2_mode.csv"))
    }
  }
  schema <- file.path(mode_dir, "fixed_o2_attractor_mode_by_seed_o2.schema.tsv")
  o2pl_write_schema(mode, "fixed_o2_attractor_mode_by_seed_o2", schema)
  o2pl_record_artifact(root_dir, "simulation", "fixed_o2_attractor_mode_by_seed_o2", common_path, mode, schema, "generate_parameter_landscape_simulation_tables.R", "invivo", input_dir)
  invisible(list(mode_by_seed_o2 = common_path, mode_summary = summary_path))
}

rm(.o2pl_simulation_dir, .o2pl_workflow_root)
