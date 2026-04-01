#!/usr/bin/env Rscript

# Usage:
#   Rscript oxygen/code/O2_supply_demand_MAP/profile_likelihood_O2_supply_demand_MAP.R \
#     --config=/abs/path/to/O2_supply_demand.yaml \
#     --baseline_seed_dir=/abs/path/to/baseline/seed_dir \
#     --profile_bounds_table=/abs/path/to/parameter_table_input.csv \
#     --output_root=/abs/path/to/output_root \
#     --param_index=1 \
#     --max_steps_per_direction=20 \
#     --seeds_per_step=5 \
#     --n_cores=62
#
# Optional:
#   --param_name=<parameter name>
#   --target_delta_objective=0.2
#   --ci_delta_threshold=1.92
#   --step_fraction_initial=0.10
#   --step_fraction_min=1e-6
#   --step_fraction_max=0.30
#   --boundary_start_tolerance=1e-8
#   --max_attempts_per_step=5
#   --use_soft_prior_for_profile=TRUE|FALSE
#   --lambda_prior_for_profile=<numeric>

.o2sd_profile_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2sd_profile_bootstrap_script_dir, mustWork = FALSE)
source(file.path(.o2sd_profile_bootstrap_script_dir, "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_profile_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool

DEFAULT_MAX_STEPS_PER_DIRECTION <- 20L
DEFAULT_SEEDS_PER_STEP <- 5L
DEFAULT_N_CORES <- 62L
DEFAULT_TARGET_DELTA_OBJECTIVE <- 0.2
DEFAULT_CI_DELTA_THRESHOLD <- 1.92
DEFAULT_STEP_FRACTION_INITIAL <- 0.10
DEFAULT_STEP_FRACTION_MIN <- 1e-6
DEFAULT_STEP_FRACTION_MAX <- 0.30
DEFAULT_BOUNDARY_START_TOLERANCE <- 1e-8
DEFAULT_MAX_ATTEMPTS_PER_STEP <- 5L
DEFAULT_SEED_BASE <- 300000L

default_config_path <- function(script_dir = SCRIPT_DIR) {
  normalizePath(file.path(script_dir, "..", "..", "config", "O2_supply_demand.yaml"), mustWork = FALSE)
}

default_results_root <- function(script_dir = SCRIPT_DIR) {
  normalizePath(file.path(script_dir, "..", "..", "results"), mustWork = FALSE)
}

default_baseline_seed_dir <- function(script_dir = SCRIPT_DIR) {
  normalizePath(
    file.path(
      script_dir,
      "..", "..", "results",
      "fit_invivo_o2_supply_demand_eq21_20260331_011709",
      "seed2"
    ),
    mustWork = FALSE
  )
}

resolve_path_value <- function(path_value, base_dir) {
  txt <- path_value
  if (is.null(txt) || !length(txt)) return(NULL)
  txt <- as.character(txt[[1]])
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  if (startsWith(txt, "~")) {
    return(normalizePath(path.expand(txt), mustWork = FALSE))
  }
  if (grepl("^(/|[A-Za-z]:[/\\\\])", txt)) {
    return(normalizePath(txt, mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, txt), mustWork = FALSE)
}

normalize_estimate_column <- function(x) {
  vapply(
    x,
    function(one) isTRUE(as_bool(one, FALSE)),
    logical(1)
  )
}

read_natural_parameter_table <- function(path) {
  if (!file.exists(path)) {
    stop("Natural-scale parameter table was not found: ", path)
  }
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("param_symbol", "estimate", "init_value", "lower_bound", "upper_bound")
  missing_cols <- setdiff(required, names(tab))
  if (length(missing_cols) > 0L) {
    stop("Missing required parameter-table columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }
  if (!("description" %in% names(tab))) tab$description <- ""
  if (!("source" %in% names(tab))) tab$source <- ""
  tab$param_symbol <- trimws(as.character(tab$param_symbol))
  if (any(!nzchar(tab$param_symbol))) {
    stop("Blank param_symbol values were found in ", path)
  }
  if (any(duplicated(tab$param_symbol))) {
    dup <- unique(tab$param_symbol[duplicated(tab$param_symbol)])
    stop("Duplicated param_symbol values were found in ", path, ": ", paste(dup, collapse = ", "))
  }
  tab$estimate <- normalize_estimate_column(tab$estimate)
  for (col in c("init_value", "lower_bound", "upper_bound")) {
    tab[[col]] <- suppressWarnings(as.numeric(tab[[col]]))
    if (any(!is.finite(tab[[col]]))) {
      bad_rows <- which(!is.finite(tab[[col]]))
      stop("Non-finite values were found in ", col, " at rows: ", paste(bad_rows, collapse = ", "))
    }
  }
  tab$description <- as.character(tab$description)
  tab$source <- as.character(tab$source)
  tab
}

write_natural_parameter_table <- function(tab, path) {
  out <- tab
  out$estimate <- ifelse(out$estimate, "TRUE", "FALSE")
  utils::write.table(
    out,
    file = path,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(path)
}

read_metric_map <- function(path, key_col, value_col) {
  if (!file.exists(path)) {
    stop("Required file was not found: ", path)
  }
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c(key_col, value_col) %in% names(tab))) {
    stop("Missing required columns in ", path, ": ", key_col, ", ", value_col)
  }
  setNames(tab[[value_col]], tab[[key_col]])
}

sync_table_init_values_from_best <- function(parameter_table, best_params_map) {
  tab <- parameter_table
  for (i in seq_len(nrow(tab))) {
    if (!isTRUE(tab$estimate[[i]])) next
    sym <- as.character(tab$param_symbol[[i]])
    if (!sym %in% names(best_params_map)) next
    val <- suppressWarnings(as.numeric(best_params_map[[sym]]))
    if (!is.finite(val)) {
      stop("Non-finite best-fit value was found for parameter ", sym)
    }
    tab$init_value[[i]] <- val
  }
  tab
}

target_parameters_from_table <- function(parameter_table) {
  target_rows <- which(vapply(parameter_table$estimate, isTRUE, logical(1)))
  if (length(target_rows) == 0L) {
    stop("No estimate=TRUE parameters were found in the profile bounds table.")
  }
  data.frame(
    param_index = seq_along(target_rows),
    row_index = target_rows,
    param_symbol = as.character(parameter_table$param_symbol[target_rows]),
    stringsAsFactors = FALSE
  )
}

resolve_profile_parameter <- function(target_table, argv) {
  if (!is.null(argv$param_name) && nzchar(trimws(as.character(argv$param_name)))) {
    param_name <- trimws(as.character(argv$param_name))
    hit <- which(as.character(target_table$param_symbol) == param_name)
    if (length(hit) != 1L) {
      stop("Parameter ", param_name, " was not found in the estimated-parameter target list.")
    }
    return(target_table[hit, , drop = FALSE])
  }
  param_index <- as_int(argv$param_index, NA_integer_)
  if (!is.finite(param_index) || param_index < 1L || param_index > nrow(target_table)) {
    stop("param_index must be between 1 and ", nrow(target_table), ".")
  }
  target_table[param_index, , drop = FALSE]
}

boundary_tolerance_abs <- function(lower_bound, upper_bound, tolerance_fraction) {
  width <- upper_bound - lower_bound
  max(as.numeric(tolerance_fraction) * width, sqrt(.Machine$double.eps))
}

classify_boundary_position <- function(value, lower_bound, upper_bound, tolerance_fraction) {
  tol_abs <- boundary_tolerance_abs(lower_bound, upper_bound, tolerance_fraction)
  at_lower <- abs(value - lower_bound) <= tol_abs
  at_upper <- abs(value - upper_bound) <= tol_abs
  rel <- (value - lower_bound) / (upper_bound - lower_bound)
  list(
    rel = rel,
    at_lower = at_lower,
    at_upper = at_upper,
    tol_abs = tol_abs
  )
}

direction_name_to_sign <- function(direction_name) {
  if (identical(direction_name, "decreasing")) return(-1L)
  if (identical(direction_name, "increasing")) return(1L)
  stop("Unsupported direction: ", direction_name)
}

direction_code <- function(direction_name) {
  if (identical(direction_name, "decreasing")) return(0L)
  if (identical(direction_name, "increasing")) return(1L)
  stop("Unsupported direction: ", direction_name)
}

direction_blocked_at_start <- function(direction_name, boundary_info) {
  if (identical(direction_name, "decreasing") && isTRUE(boundary_info$at_lower)) return(TRUE)
  if (identical(direction_name, "increasing") && isTRUE(boundary_info$at_upper)) return(TRUE)
  FALSE
}

start_boundary_side <- function(direction_name, boundary_info) {
  if (direction_blocked_at_start(direction_name, boundary_info)) {
    if (identical(direction_name, "decreasing")) return("lower")
    return("upper")
  }
  "none"
}

generate_profile_seeds <- function(param_index, direction_name, step_index, max_steps_per_direction, seeds_per_step) {
  direction_width <- as.integer(max(1000L, as.integer(max_steps_per_direction) * as.integer(seeds_per_step) + 100L))
  param_width <- as.integer(2L * direction_width + 1000L)
  base <- as.integer(
    DEFAULT_SEED_BASE +
      (as.integer(param_index) - 1L) * param_width +
      direction_code(direction_name) * direction_width +
      (as.integer(step_index) - 1L) * as.integer(seeds_per_step)
  )
  as.integer(base + seq_len(as.integer(seeds_per_step)))
}

write_seed_manifest <- function(seeds, path) {
  utils::write.table(
    data.frame(seed = as.integer(seeds), stringsAsFactors = FALSE),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  invisible(path)
}

make_fixed_parameter_table <- function(base_table, target_symbol, fixed_value) {
  tab <- base_table
  hit <- which(as.character(tab$param_symbol) == target_symbol)
  if (length(hit) != 1L) {
    stop("Target parameter ", target_symbol, " was not found exactly once in the parameter table.")
  }
  idx <- hit[[1]]
  tab$init_value[[idx]] <- fixed_value
  tab$lower_bound[[idx]] <- fixed_value
  tab$upper_bound[[idx]] <- fixed_value
  tab$estimate[[idx]] <- FALSE
  tab
}

propose_step_candidate <- function(current_value, direction_name, step_fraction, lower_bound, upper_bound, tolerance_fraction) {
  sign <- direction_name_to_sign(direction_name)
  width <- upper_bound - lower_bound
  tol_abs <- boundary_tolerance_abs(lower_bound, upper_bound, tolerance_fraction)
  proposed <- current_value + sign * as.numeric(step_fraction) * width
  clipped <- min(max(proposed, lower_bound), upper_bound)
  boundary_hit <- FALSE
  boundary_side <- "none"
  if (clipped <= lower_bound + tol_abs) {
    clipped <- lower_bound
    boundary_hit <- TRUE
    boundary_side <- "lower"
  } else if (clipped >= upper_bound - tol_abs) {
    clipped <- upper_bound
    boundary_hit <- TRUE
    boundary_side <- "upper"
  }
  if (abs(clipped - current_value) <= tol_abs) {
    return(list(
      can_move = FALSE,
      candidate_value = current_value,
      boundary_hit = boundary_hit,
      boundary_side = boundary_side,
      stop_reason = "step_below_tolerance_or_no_room"
    ))
  }
  list(
    can_move = TRUE,
    candidate_value = clipped,
    boundary_hit = boundary_hit,
    boundary_side = boundary_side,
    stop_reason = ""
  )
}

adapt_step_fraction <- function(current_fraction, delta_from_previous, target_delta, step_fraction_min, step_fraction_max) {
  if (!is.finite(delta_from_previous) || delta_from_previous <= 1e-12) {
    return(min(step_fraction_max, max(step_fraction_min, current_fraction * 2)))
  }
  ratio <- target_delta / delta_from_previous
  if (!is.finite(ratio) || ratio <= 0) {
    return(min(step_fraction_max, max(step_fraction_min, current_fraction)))
  }
  out <- current_fraction * sqrt(ratio)
  min(step_fraction_max, max(step_fraction_min, out))
}

run_single_seed_fit <- function(
  fit_script,
  config_path,
  parameter_table_path,
  runs_root,
  logs_root,
  run_prefix,
  seed,
  n_cores,
  warm_start_path = NULL,
  use_soft_prior_for_profile = NULL,
  lambda_prior_for_profile = NULL
) {
  dir.create(runs_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_root, recursive = TRUE, showWarnings = FALSE)
  run_dir <- file.path(runs_root, run_prefix)
  seed_dir <- file.path(run_dir, paste0("seed", as.integer(seed)))
  call_log <- file.path(logs_root, paste0(run_prefix, ".log"))
  args <- c(
    normalizePath(fit_script, mustWork = FALSE),
    "--mode=run",
    paste0("--config=", normalizePath(config_path, mustWork = FALSE)),
    paste0("--parameter_table=", normalizePath(parameter_table_path, mustWork = FALSE)),
    paste0("--out_root=", normalizePath(runs_root, mustWork = FALSE)),
    paste0("--run_prefix=", run_prefix),
    paste0("--seeds_csv=", as.integer(seed)),
    "--append_run_prefix_timestamp=FALSE",
    "--auto_viz=FALSE",
    "--predict_n_cores=1",
    paste0("--n_cores=", as.integer(n_cores))
  )
  if (!is.null(warm_start_path) && nzchar(trimws(as.character(warm_start_path)))) {
    args <- c(args, paste0("--init_params_tsv=", normalizePath(warm_start_path, mustWork = FALSE)))
  }
  if (!is.null(use_soft_prior_for_profile)) {
    args <- c(args, paste0("--use_soft_prior=", ifelse(isTRUE(use_soft_prior_for_profile), "TRUE", "FALSE")))
  }
  if (!is.null(lambda_prior_for_profile) && is.finite(lambda_prior_for_profile)) {
    args <- c(args, paste0("--lambda_prior=", as.numeric(lambda_prior_for_profile)))
  }
  if (as.integer(n_cores) > 1L) {
    args <- c(args, "--deoptim_parallel=TRUE")
  }
  status <- tryCatch(
    system2("Rscript", args = args, stdout = call_log, stderr = call_log, wait = TRUE),
    error = function(e) {
      writeLines(paste0("Rscript launch failed: ", conditionMessage(e)), con = call_log, sep = "\n")
      999L
    }
  )
  list(
    status = as.integer(status %||% 0L),
    run_prefix = run_prefix,
    run_dir = run_dir,
    seed_dir = seed_dir,
    call_log = call_log
  )
}

read_single_seed_result <- function(
  seed_run_info,
  target_symbol,
  param_index,
  param_symbol,
  direction_name,
  step_index,
  attempt_index,
  fixed_value,
  step_fraction_requested
) {
  fit_summary_path <- file.path(seed_run_info$seed_dir, "fit_summary.tsv")
  best_params_path <- file.path(seed_run_info$seed_dir, "best_params.tsv")
  fit_summary_present <- file.exists(fit_summary_path)
  best_params_present <- file.exists(best_params_path)

  objective <- NA_real_
  objective_data <- NA_real_
  objective_prior <- NA_real_
  objective_burden <- NA_real_
  objective_ploidy <- NA_real_
  selected_value_reported <- NA_real_
  status_note <- character(0)
  best_params_map <- NULL

  if (fit_summary_present) {
    fit_vals <- tryCatch(
      read_metric_map(fit_summary_path, "metric", "value"),
      error = function(e) NULL
    )
    if (!is.null(fit_vals)) {
      objective <- suppressWarnings(as.numeric(fit_vals[["objective"]]))
      objective_data <- suppressWarnings(as.numeric(fit_vals[["objective_data"]]))
      objective_prior <- suppressWarnings(as.numeric(fit_vals[["objective_prior"]]))
      objective_burden <- suppressWarnings(as.numeric(fit_vals[["objective_burden"]]))
      objective_ploidy <- suppressWarnings(as.numeric(fit_vals[["objective_ploidy"]]))
    } else {
      status_note <- c(status_note, "fit_summary_parse_failed")
    }
  } else {
    status_note <- c(status_note, "missing_fit_summary")
  }

  if (best_params_present) {
    best_params_map <- tryCatch(
      read_metric_map(best_params_path, "parameter", "value"),
      error = function(e) NULL
    )
    if (!is.null(best_params_map) && target_symbol %in% names(best_params_map)) {
      selected_value_reported <- suppressWarnings(as.numeric(best_params_map[[target_symbol]]))
    } else {
      status_note <- c(status_note, "missing_target_in_best_params")
    }
  } else {
    status_note <- c(status_note, "missing_best_params")
  }

  complete <- is.finite(objective) && is.finite(selected_value_reported)
  if (!complete && length(status_note) == 0L) {
    status_note <- "non_finite_metrics"
  }

  row <- data.frame(
    param_index = as.integer(param_index),
    param_symbol = as.character(param_symbol),
    direction = as.character(direction_name),
    step_index = as.integer(step_index),
    attempt_index = as.integer(attempt_index),
    fixed_value = as.numeric(fixed_value),
    step_fraction_requested = as.numeric(step_fraction_requested),
    seed = as.integer(sub("^seed", "", basename(seed_run_info$seed_dir))),
    run_exit_status = as.integer(seed_run_info$status),
    run_prefix = as.character(seed_run_info$run_prefix),
    run_dir = normalizePath(seed_run_info$run_dir, mustWork = FALSE),
    seed_dir = normalizePath(seed_run_info$seed_dir, mustWork = FALSE),
    call_log = normalizePath(seed_run_info$call_log, mustWork = FALSE),
    fit_summary_present = fit_summary_present,
    best_params_present = best_params_present,
    complete = complete,
    objective = objective,
    objective_data = objective_data,
    objective_prior = objective_prior,
    objective_burden = objective_burden,
    objective_ploidy = objective_ploidy,
    selected_value_reported = selected_value_reported,
    best_params_path = if (best_params_present) normalizePath(best_params_path, mustWork = FALSE) else NA_character_,
    status_note = paste(status_note, collapse = ";"),
    stringsAsFactors = FALSE
  )
  list(row = row, best_params_map = best_params_map)
}

best_params_to_long <- function(best_params_map, param_index, target_symbol, direction_name, step_index, fixed_value, point_type) {
  data.frame(
    param_index = as.integer(param_index),
    target_param_symbol = as.character(target_symbol),
    direction = as.character(direction_name),
    step_index = as.integer(step_index),
    fixed_value = as.numeric(fixed_value),
    point_type = as.character(point_type),
    parameter = names(best_params_map),
    value = suppressWarnings(as.numeric(best_params_map)),
    stringsAsFactors = FALSE
  )
}

select_best_complete_seed <- function(seed_df) {
  complete_idx <- which(vapply(seed_df$complete, isTRUE, logical(1)))
  if (length(complete_idx) == 0L) return(NULL)
  ord <- order(seed_df$objective[complete_idx], seed_df$seed[complete_idx], na.last = TRUE)
  seed_df[complete_idx[ord[[1]]], , drop = FALSE]
}

run_step_attempt <- function(
  param_index,
  param_symbol,
  direction_name,
  step_index,
  attempt_index,
  candidate_value,
  step_fraction_requested,
  seeds,
  working_table,
  step_dir,
  fit_script,
  config_path,
  runs_root,
  logs_root,
  n_cores,
  warm_start_path,
  use_soft_prior_for_profile,
  lambda_prior_for_profile
) {
  dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
  parameter_table_path <- file.path(step_dir, "parameter_table.fixed.csv")
  temp_table <- make_fixed_parameter_table(working_table, param_symbol, candidate_value)
  write_natural_parameter_table(temp_table, parameter_table_path)
  write_seed_manifest(seeds, file.path(step_dir, "step_seeds.tsv"))

  seed_rows <- vector("list", length(seeds))
  best_maps <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    seed <- as.integer(seeds[[i]])
    run_prefix <- paste0(
      direction_name,
      "_step", sprintf("%02d", as.integer(step_index)),
      "_att", sprintf("%02d", as.integer(attempt_index)),
      "_seed", as.integer(seed)
    )
    seed_run <- run_single_seed_fit(
      fit_script = fit_script,
      config_path = config_path,
      parameter_table_path = parameter_table_path,
      runs_root = runs_root,
      logs_root = logs_root,
      run_prefix = run_prefix,
      seed = seed,
      n_cores = n_cores,
      warm_start_path = warm_start_path,
      use_soft_prior_for_profile = use_soft_prior_for_profile,
      lambda_prior_for_profile = lambda_prior_for_profile
    )
    parsed <- read_single_seed_result(
      seed_run_info = seed_run,
      target_symbol = param_symbol,
      param_index = param_index,
      param_symbol = param_symbol,
      direction_name = direction_name,
      step_index = step_index,
      attempt_index = attempt_index,
      fixed_value = candidate_value,
      step_fraction_requested = step_fraction_requested
    )
    seed_rows[[i]] <- parsed$row
    best_maps[[i]] <- parsed$best_params_map
  }
  seed_df <- do.call(rbind, seed_rows)
  best_row <- select_best_complete_seed(seed_df)
  best_map <- NULL
  if (!is.null(best_row)) {
    best_hit <- which(seed_df$seed == best_row$seed[[1]])[[1]]
    best_map <- best_maps[[best_hit]]
  }
  list(
    seed_df = seed_df,
    best_row = best_row,
    best_map = best_map,
    parameter_table_path = parameter_table_path
  )
}

run_direction_profile <- function(
  direction_name,
  param_index,
  param_symbol,
  lower_bound,
  upper_bound,
  baseline_value,
  baseline_objective,
  baseline_best_path,
  baseline_best_map,
  bounds_table,
  param_dir,
  fit_script,
  config_path,
  n_cores,
  max_steps_per_direction,
  seeds_per_step,
  step_fraction_initial,
  step_fraction_min,
  step_fraction_max,
  target_delta_objective,
  ci_delta_threshold,
  boundary_start_tolerance,
  max_attempts_per_step,
  use_soft_prior_for_profile,
  lambda_prior_for_profile
) {
  boundary_info <- classify_boundary_position(
    value = baseline_value,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    tolerance_fraction = boundary_start_tolerance
  )
  blocked_at_start <- direction_blocked_at_start(direction_name, boundary_info)
  start_side <- start_boundary_side(direction_name, boundary_info)
  steps_dir <- file.path(param_dir, "steps", direction_name)
  runs_root <- file.path(param_dir, "runs")
  logs_root <- file.path(param_dir, "logs")
  dir.create(steps_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(runs_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_root, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(blocked_at_start)) {
    summary_row <- data.frame(
      param_index = as.integer(param_index),
      param_symbol = as.character(param_symbol),
      direction = as.character(direction_name),
      start_blocked = TRUE,
      start_boundary_side = as.character(start_side),
      steps_completed = 0L,
      termination_reason = "blocked_at_start_boundary",
      last_fixed_value = baseline_value,
      last_objective = baseline_objective,
      max_delta_objective_vs_baseline = 0,
      threshold_reached = FALSE,
      boundary_reached = TRUE,
      stringsAsFactors = FALSE
    )
    return(list(
      step_df = data.frame(),
      seed_df = data.frame(),
      relation_df = data.frame(),
      summary_df = summary_row
    ))
  }

  current_value <- baseline_value
  current_objective <- baseline_objective
  current_table <- sync_table_init_values_from_best(bounds_table, baseline_best_map)
  current_best_map <- baseline_best_map
  current_warm_start_path <- baseline_best_path
  step_fraction_current <- as.numeric(step_fraction_initial)
  step_rows <- list()
  seed_rows <- list()
  relation_rows <- list()
  termination_reason <- "max_steps_reached"
  threshold_reached <- FALSE
  boundary_reached <- FALSE

  for (step_index in seq_len(as.integer(max_steps_per_direction))) {
    attempt_fraction <- step_fraction_current
    step_success <- FALSE
    for (attempt_index in seq_len(as.integer(max_attempts_per_step))) {
      candidate <- propose_step_candidate(
        current_value = current_value,
        direction_name = direction_name,
        step_fraction = attempt_fraction,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        tolerance_fraction = boundary_start_tolerance
      )
      if (!isTRUE(candidate$can_move)) {
        termination_reason <- candidate$stop_reason
        step_success <- FALSE
        break
      }

      step_dir <- file.path(
        steps_dir,
        paste0("step_", sprintf("%02d", as.integer(step_index))),
        paste0("attempt_", sprintf("%02d", as.integer(attempt_index)))
      )
      seeds <- generate_profile_seeds(
        param_index = param_index,
        direction_name = direction_name,
        step_index = step_index,
        max_steps_per_direction = max_steps_per_direction,
        seeds_per_step = seeds_per_step
      )
      attempted <- run_step_attempt(
        param_index = param_index,
        param_symbol = param_symbol,
        direction_name = direction_name,
        step_index = step_index,
        attempt_index = attempt_index,
        candidate_value = candidate$candidate_value,
        step_fraction_requested = attempt_fraction,
        seeds = seeds,
        working_table = current_table,
        step_dir = step_dir,
        fit_script = fit_script,
        config_path = config_path,
        runs_root = runs_root,
        logs_root = logs_root,
        n_cores = n_cores,
        warm_start_path = current_warm_start_path,
        use_soft_prior_for_profile = use_soft_prior_for_profile,
        lambda_prior_for_profile = lambda_prior_for_profile
      )
      seed_rows[[length(seed_rows) + 1L]] <- attempted$seed_df

      if (is.null(attempted$best_row)) {
        if (attempt_fraction <= step_fraction_min * (1 + 1e-12)) {
          termination_reason <- "no_complete_seed_at_min_step"
          step_success <- FALSE
          break
        }
        attempt_fraction <- max(step_fraction_min, attempt_fraction / 2)
        next
      }

      best_row <- attempted$best_row
      best_map <- attempted$best_map
      if (is.null(best_map)) {
        termination_reason <- "best_seed_missing_best_params_map"
        step_success <- FALSE
        break
      }

      delta_vs_baseline <- as.numeric(best_row$objective[[1]]) - baseline_objective
      delta_vs_previous <- as.numeric(best_row$objective[[1]]) - current_objective
      threshold_reached_now <- is.finite(delta_vs_baseline) && delta_vs_baseline >= ci_delta_threshold
      next_step_fraction <- adapt_step_fraction(
        current_fraction = attempt_fraction,
        delta_from_previous = delta_vs_previous,
        target_delta = target_delta_objective,
        step_fraction_min = step_fraction_min,
        step_fraction_max = step_fraction_max
      )

      step_row <- data.frame(
        param_index = as.integer(param_index),
        param_symbol = as.character(param_symbol),
        direction = as.character(direction_name),
        step_index = as.integer(step_index),
        attempt_index_used = as.integer(attempt_index),
        fixed_value_previous = as.numeric(current_value),
        fixed_value = as.numeric(candidate$candidate_value),
        step_fraction_requested = as.numeric(step_fraction_current),
        step_fraction_used = as.numeric(attempt_fraction),
        step_fraction_next = as.numeric(next_step_fraction),
        requested_seeds = paste(as.integer(seeds), collapse = ","),
        requested_seed_count = as.integer(length(seeds)),
        completed_seed_count = as.integer(sum(vapply(attempted$seed_df$complete, isTRUE, logical(1)))),
        best_seed = as.integer(best_row$seed[[1]]),
        best_run_dir = as.character(best_row$run_dir[[1]]),
        best_seed_dir = as.character(best_row$seed_dir[[1]]),
        best_params_path = as.character(best_row$best_params_path[[1]]),
        warm_start_path = normalizePath(current_warm_start_path, mustWork = FALSE),
        objective = as.numeric(best_row$objective[[1]]),
        objective_data = as.numeric(best_row$objective_data[[1]]),
        objective_prior = as.numeric(best_row$objective_prior[[1]]),
        objective_burden = as.numeric(best_row$objective_burden[[1]]),
        objective_ploidy = as.numeric(best_row$objective_ploidy[[1]]),
        selected_value_reported = as.numeric(best_row$selected_value_reported[[1]]),
        delta_objective_vs_baseline = as.numeric(delta_vs_baseline),
        delta_objective_vs_previous = as.numeric(delta_vs_previous),
        boundary_hit = isTRUE(candidate$boundary_hit),
        boundary_side = as.character(candidate$boundary_side),
        threshold_reached = isTRUE(threshold_reached_now),
        stop_after_step = FALSE,
        step_status = "complete",
        stringsAsFactors = FALSE
      )

      relation_rows[[length(relation_rows) + 1L]] <- best_params_to_long(
        best_params_map = best_map,
        param_index = param_index,
        target_symbol = param_symbol,
        direction_name = direction_name,
        step_index = step_index,
        fixed_value = candidate$candidate_value,
        point_type = "profile_step"
      )
      step_rows[[length(step_rows) + 1L]] <- step_row

      current_value <- as.numeric(candidate$candidate_value)
      current_objective <- as.numeric(best_row$objective[[1]])
      current_best_map <- best_map
      current_table <- sync_table_init_values_from_best(current_table, current_best_map)
      current_warm_start_path <- as.character(best_row$best_params_path[[1]])
      step_fraction_current <- as.numeric(next_step_fraction)
      step_success <- TRUE

      if (isTRUE(candidate$boundary_hit)) {
        termination_reason <- "boundary_reached"
        boundary_reached <- TRUE
      }
      if (isTRUE(threshold_reached_now)) {
        termination_reason <- "ci_threshold_reached"
        threshold_reached <- TRUE
      }
      if (isTRUE(candidate$boundary_hit) || isTRUE(threshold_reached_now)) {
        step_rows[[length(step_rows)]][["stop_after_step"]] <- TRUE
        break
      }
      break
    }

    if (!isTRUE(step_success)) {
      break
    }
    if (identical(termination_reason, "boundary_reached") || identical(termination_reason, "ci_threshold_reached")) {
      break
    }
  }

  step_df <- if (length(step_rows) > 0L) do.call(rbind, step_rows) else data.frame()
  seed_df <- if (length(seed_rows) > 0L) do.call(rbind, seed_rows) else data.frame()
  relation_df <- if (length(relation_rows) > 0L) do.call(rbind, relation_rows) else data.frame()
  steps_completed <- if (nrow(step_df)) nrow(step_df) else 0L
  last_fixed_value <- if (nrow(step_df)) as.numeric(step_df$fixed_value[[nrow(step_df)]]) else baseline_value
  last_objective <- if (nrow(step_df)) as.numeric(step_df$objective[[nrow(step_df)]]) else baseline_objective
  max_delta <- if (nrow(step_df)) max(step_df$delta_objective_vs_baseline, na.rm = TRUE) else 0
  if (!is.finite(max_delta)) max_delta <- 0

  summary_row <- data.frame(
    param_index = as.integer(param_index),
    param_symbol = as.character(param_symbol),
    direction = as.character(direction_name),
    start_blocked = FALSE,
    start_boundary_side = "none",
    steps_completed = as.integer(steps_completed),
    termination_reason = as.character(termination_reason),
    last_fixed_value = as.numeric(last_fixed_value),
    last_objective = as.numeric(last_objective),
    max_delta_objective_vs_baseline = as.numeric(max_delta),
    threshold_reached = isTRUE(threshold_reached),
    boundary_reached = isTRUE(boundary_reached),
    stringsAsFactors = FALSE
  )

  list(
    step_df = step_df,
    seed_df = seed_df,
    relation_df = relation_df,
    summary_df = summary_row
  )
}

write_table_tsv <- function(tab, path) {
  utils::write.table(
    tab,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(path)
}

write_profile_plots <- function(combined_df, param_dir, param_symbol, lower_bound, upper_bound, baseline_value) {
  complete <- combined_df[is.finite(combined_df$objective), , drop = FALSE]
  if (!nrow(complete)) return(invisible(NULL))

  ord <- order(complete$fixed_value)
  x <- as.numeric(complete$fixed_value[ord])
  y_objective <- as.numeric(complete$objective[ord])
  y_delta <- as.numeric(complete$delta_objective_vs_baseline[ord])

  pdf(file.path(param_dir, "profile_curve.pdf"), width = 8, height = 6)
  plot(
    x,
    y_objective,
    type = "b",
    pch = 19,
    xlab = param_symbol,
    ylab = "Objective",
    main = paste0("Profile objective path: ", param_symbol)
  )
  abline(v = baseline_value, col = "red", lty = 2, lwd = 2)
  abline(v = lower_bound, col = "grey60", lty = 3)
  abline(v = upper_bound, col = "grey60", lty = 3)
  dev.off()

  pdf(file.path(param_dir, "profile_delta_curve.pdf"), width = 8, height = 6)
  plot(
    x,
    y_delta,
    type = "b",
    pch = 19,
    xlab = param_symbol,
    ylab = "Delta objective vs baseline",
    main = paste0("Profile delta path: ", param_symbol)
  )
  abline(h = 0, col = "grey40", lty = 2)
  abline(v = baseline_value, col = "red", lty = 2, lwd = 2)
  abline(v = lower_bound, col = "grey60", lty = 3)
  abline(v = upper_bound, col = "grey60", lty = 3)
  dev.off()
  invisible(NULL)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  config_path <- resolve_path_value(argv$config, getwd()) %||% default_config_path()
  config_path <- normalizePath(config_path, mustWork = TRUE)

  baseline_seed_dir <- resolve_path_value(argv$baseline_seed_dir, getwd()) %||% default_baseline_seed_dir()
  baseline_seed_dir <- normalizePath(baseline_seed_dir, mustWork = TRUE)

  profile_bounds_table <- resolve_path_value(argv$profile_bounds_table, getwd())
  if (is.null(profile_bounds_table)) {
    profile_bounds_table <- file.path(baseline_seed_dir, "parameter_table_input.csv")
  }
  profile_bounds_table <- normalizePath(profile_bounds_table, mustWork = TRUE)

  output_root <- resolve_path_value(argv$output_root, getwd())
  if (is.null(output_root)) {
    output_root <- file.path(
      default_results_root(),
      paste0("profile_likelihood_O2_supply_demand_MAP_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
  }
  output_root <- normalizePath(output_root, mustWork = FALSE)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  max_steps_per_direction <- as_int(argv$max_steps_per_direction, as_int(argv$grid_points, DEFAULT_MAX_STEPS_PER_DIRECTION))
  seeds_per_step <- as_int(argv$seeds_per_step, as_int(argv$seeds_per_point, DEFAULT_SEEDS_PER_STEP))
  n_cores <- as_int(argv$n_cores, DEFAULT_N_CORES)
  target_delta_objective <- as_num(argv$target_delta_objective, DEFAULT_TARGET_DELTA_OBJECTIVE)
  ci_delta_threshold <- as_num(argv$ci_delta_threshold, DEFAULT_CI_DELTA_THRESHOLD)
  step_fraction_initial <- as_num(argv$step_fraction_initial, DEFAULT_STEP_FRACTION_INITIAL)
  step_fraction_min <- as_num(argv$step_fraction_min, DEFAULT_STEP_FRACTION_MIN)
  step_fraction_max <- as_num(argv$step_fraction_max, DEFAULT_STEP_FRACTION_MAX)
  boundary_start_tolerance <- as_num(argv$boundary_start_tolerance, DEFAULT_BOUNDARY_START_TOLERANCE)
  max_attempts_per_step <- as_int(argv$max_attempts_per_step, DEFAULT_MAX_ATTEMPTS_PER_STEP)

  if (!is.finite(max_steps_per_direction) || max_steps_per_direction < 1L) stop("max_steps_per_direction must be >= 1.")
  if (!is.finite(seeds_per_step) || seeds_per_step < 1L) stop("seeds_per_step must be >= 1.")
  if (!is.finite(n_cores) || n_cores < 1L) stop("n_cores must be >= 1.")
  if (!is.finite(target_delta_objective) || target_delta_objective <= 0) stop("target_delta_objective must be > 0.")
  if (!is.finite(ci_delta_threshold) || ci_delta_threshold <= 0) stop("ci_delta_threshold must be > 0.")
  if (!is.finite(step_fraction_initial) || step_fraction_initial <= 0) stop("step_fraction_initial must be > 0.")
  if (!is.finite(step_fraction_min) || step_fraction_min <= 0) stop("step_fraction_min must be > 0.")
  if (!is.finite(step_fraction_max) || step_fraction_max < step_fraction_initial) stop("step_fraction_max must be >= step_fraction_initial.")
  if (!is.finite(boundary_start_tolerance) || boundary_start_tolerance < 0) stop("boundary_start_tolerance must be >= 0.")
  if (!is.finite(max_attempts_per_step) || max_attempts_per_step < 1L) stop("max_attempts_per_step must be >= 1.")

  fit_script <- normalizePath(file.path(SCRIPT_DIR, "fit_invivo_model_O2_supply_demand_MAP.R"), mustWork = TRUE)
  baseline_best_path <- file.path(baseline_seed_dir, "best_params.tsv")
  baseline_fit_summary_path <- file.path(baseline_seed_dir, "fit_summary.tsv")
  baseline_best_map <- read_metric_map(baseline_best_path, "parameter", "value")
  baseline_fit_summary <- read_metric_map(baseline_fit_summary_path, "metric", "value")

  baseline_objective <- suppressWarnings(as.numeric(baseline_fit_summary[["objective"]]))
  baseline_objective_data <- suppressWarnings(as.numeric(baseline_fit_summary[["objective_data"]]))
  baseline_objective_prior <- suppressWarnings(as.numeric(baseline_fit_summary[["objective_prior"]]))
  baseline_objective_burden <- suppressWarnings(as.numeric(baseline_fit_summary[["objective_burden"]]))
  baseline_objective_ploidy <- suppressWarnings(as.numeric(baseline_fit_summary[["objective_ploidy"]]))
  if (!is.finite(baseline_objective)) {
    stop("Baseline objective could not be read from ", baseline_fit_summary_path)
  }

  baseline_use_soft_prior <- as_bool(baseline_fit_summary[["use_soft_prior"]], TRUE)
  use_soft_prior_for_profile <- if (!is.null(argv$use_soft_prior_for_profile)) {
    as_bool(argv$use_soft_prior_for_profile, baseline_use_soft_prior)
  } else {
    baseline_use_soft_prior
  }
  baseline_lambda_prior <- as_num(baseline_fit_summary[["lambda_prior"]], NA_real_)
  lambda_prior_for_profile <- if (!is.null(argv$lambda_prior_for_profile)) {
    as_num(argv$lambda_prior_for_profile, baseline_lambda_prior)
  } else {
    baseline_lambda_prior
  }

  bounds_table <- read_natural_parameter_table(profile_bounds_table)
  bounds_table <- sync_table_init_values_from_best(bounds_table, baseline_best_map)
  target_table <- target_parameters_from_table(bounds_table)
  target_row <- resolve_profile_parameter(target_table, argv)
  param_index <- as.integer(target_row$param_index[[1]])
  row_index <- as.integer(target_row$row_index[[1]])
  param_symbol <- as.character(target_row$param_symbol[[1]])

  lower_bound <- as.numeric(bounds_table$lower_bound[[row_index]])
  upper_bound <- as.numeric(bounds_table$upper_bound[[row_index]])
  baseline_value <- as.numeric(bounds_table$init_value[[row_index]])
  if (!is.finite(lower_bound) || !is.finite(upper_bound) || upper_bound <= lower_bound) {
    stop("Invalid profile bounds were found for parameter ", param_symbol)
  }
  if (!is.finite(baseline_value)) {
    stop("Baseline value could not be resolved for parameter ", param_symbol)
  }
  if (baseline_value < lower_bound || baseline_value > upper_bound) {
    stop(
      "Baseline value ",
      signif(baseline_value, 8),
      " is outside profile bounds [",
      signif(lower_bound, 8),
      ", ",
      signif(upper_bound, 8),
      "] for parameter ",
      param_symbol,
      "."
    )
  }

  param_dir <- file.path(output_root, sprintf("%02d_%s", as.integer(param_index), param_symbol))
  dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)

  start_boundary <- classify_boundary_position(
    value = baseline_value,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    tolerance_fraction = boundary_start_tolerance
  )

  manifest_row <- data.frame(
    param_index = as.integer(param_index),
    param_symbol = param_symbol,
    max_steps_per_direction = as.integer(max_steps_per_direction),
    seeds_per_step = as.integer(seeds_per_step),
    n_cores = as.integer(n_cores),
    target_delta_objective = as.numeric(target_delta_objective),
    ci_delta_threshold = as.numeric(ci_delta_threshold),
    step_fraction_initial = as.numeric(step_fraction_initial),
    step_fraction_min = as.numeric(step_fraction_min),
    step_fraction_max = as.numeric(step_fraction_max),
    boundary_start_tolerance = as.numeric(boundary_start_tolerance),
    max_attempts_per_step = as.integer(max_attempts_per_step),
    lower_bound = as.numeric(lower_bound),
    upper_bound = as.numeric(upper_bound),
    baseline_value = as.numeric(baseline_value),
    baseline_objective = as.numeric(baseline_objective),
    baseline_objective_data = as.numeric(baseline_objective_data),
    baseline_objective_prior = as.numeric(baseline_objective_prior),
    baseline_objective_burden = as.numeric(baseline_objective_burden),
    baseline_objective_ploidy = as.numeric(baseline_objective_ploidy),
    baseline_use_soft_prior = isTRUE(baseline_use_soft_prior),
    use_soft_prior_for_profile = isTRUE(use_soft_prior_for_profile),
    lambda_prior_for_profile = as.numeric(lambda_prior_for_profile),
    start_at_lower_boundary = isTRUE(start_boundary$at_lower),
    start_at_upper_boundary = isTRUE(start_boundary$at_upper),
    config_path = normalizePath(config_path, mustWork = FALSE),
    baseline_seed_dir = normalizePath(baseline_seed_dir, mustWork = FALSE),
    profile_bounds_table = normalizePath(profile_bounds_table, mustWork = FALSE),
    fit_script = normalizePath(fit_script, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  write_table_tsv(manifest_row, file.path(param_dir, "parameter_manifest.tsv"))

  baseline_profile_row <- data.frame(
    param_index = as.integer(param_index),
    param_symbol = param_symbol,
    direction = "baseline",
    step_index = 0L,
    attempt_index_used = 0L,
    fixed_value_previous = as.numeric(baseline_value),
    fixed_value = as.numeric(baseline_value),
    step_fraction_requested = NA_real_,
    step_fraction_used = NA_real_,
    step_fraction_next = as.numeric(step_fraction_initial),
    requested_seeds = "",
    requested_seed_count = 0L,
    completed_seed_count = 1L,
    best_seed = as.integer(sub("^seed", "", basename(baseline_seed_dir))),
    best_run_dir = normalizePath(dirname(baseline_seed_dir), mustWork = FALSE),
    best_seed_dir = normalizePath(baseline_seed_dir, mustWork = FALSE),
    best_params_path = normalizePath(baseline_best_path, mustWork = FALSE),
    warm_start_path = normalizePath(baseline_best_path, mustWork = FALSE),
    objective = as.numeric(baseline_objective),
    objective_data = as.numeric(baseline_objective_data),
    objective_prior = as.numeric(baseline_objective_prior),
    objective_burden = as.numeric(baseline_objective_burden),
    objective_ploidy = as.numeric(baseline_objective_ploidy),
    selected_value_reported = as.numeric(baseline_value),
    delta_objective_vs_baseline = 0,
    delta_objective_vs_previous = 0,
    boundary_hit = isTRUE(start_boundary$at_lower) || isTRUE(start_boundary$at_upper),
    boundary_side = if (isTRUE(start_boundary$at_lower)) "lower" else if (isTRUE(start_boundary$at_upper)) "upper" else "none",
    threshold_reached = FALSE,
    stop_after_step = FALSE,
    step_status = "baseline",
    stringsAsFactors = FALSE
  )

  baseline_relation_df <- best_params_to_long(
    best_params_map = baseline_best_map,
    param_index = param_index,
    target_symbol = param_symbol,
    direction_name = "baseline",
    step_index = 0L,
    fixed_value = baseline_value,
    point_type = "baseline"
  )

  decreasing <- run_direction_profile(
    direction_name = "decreasing",
    param_index = param_index,
    param_symbol = param_symbol,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    baseline_value = baseline_value,
    baseline_objective = baseline_objective,
    baseline_best_path = baseline_best_path,
    baseline_best_map = baseline_best_map,
    bounds_table = bounds_table,
    param_dir = param_dir,
    fit_script = fit_script,
    config_path = config_path,
    n_cores = as.integer(n_cores),
    max_steps_per_direction = as.integer(max_steps_per_direction),
    seeds_per_step = as.integer(seeds_per_step),
    step_fraction_initial = as.numeric(step_fraction_initial),
    step_fraction_min = as.numeric(step_fraction_min),
    step_fraction_max = as.numeric(step_fraction_max),
    target_delta_objective = as.numeric(target_delta_objective),
    ci_delta_threshold = as.numeric(ci_delta_threshold),
    boundary_start_tolerance = as.numeric(boundary_start_tolerance),
    max_attempts_per_step = as.integer(max_attempts_per_step),
    use_soft_prior_for_profile = use_soft_prior_for_profile,
    lambda_prior_for_profile = lambda_prior_for_profile
  )
  increasing <- run_direction_profile(
    direction_name = "increasing",
    param_index = param_index,
    param_symbol = param_symbol,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    baseline_value = baseline_value,
    baseline_objective = baseline_objective,
    baseline_best_path = baseline_best_path,
    baseline_best_map = baseline_best_map,
    bounds_table = bounds_table,
    param_dir = param_dir,
    fit_script = fit_script,
    config_path = config_path,
    n_cores = as.integer(n_cores),
    max_steps_per_direction = as.integer(max_steps_per_direction),
    seeds_per_step = as.integer(seeds_per_step),
    step_fraction_initial = as.numeric(step_fraction_initial),
    step_fraction_min = as.numeric(step_fraction_min),
    step_fraction_max = as.numeric(step_fraction_max),
    target_delta_objective = as.numeric(target_delta_objective),
    ci_delta_threshold = as.numeric(ci_delta_threshold),
    boundary_start_tolerance = as.numeric(boundary_start_tolerance),
    max_attempts_per_step = as.integer(max_attempts_per_step),
    use_soft_prior_for_profile = use_soft_prior_for_profile,
    lambda_prior_for_profile = lambda_prior_for_profile
  )

  decreasing_path <- decreasing$step_df
  increasing_path <- increasing$step_df
  combined_path <- baseline_profile_row
  if (nrow(decreasing_path)) combined_path <- rbind(combined_path, decreasing_path)
  if (nrow(increasing_path)) combined_path <- rbind(combined_path, increasing_path)

  seed_df <- data.frame()
  if (nrow(decreasing$seed_df)) seed_df <- rbind(seed_df, decreasing$seed_df)
  if (nrow(increasing$seed_df)) seed_df <- rbind(seed_df, increasing$seed_df)

  relation_df <- baseline_relation_df
  if (nrow(decreasing$relation_df)) relation_df <- rbind(relation_df, decreasing$relation_df)
  if (nrow(increasing$relation_df)) relation_df <- rbind(relation_df, increasing$relation_df)

  direction_summary_df <- rbind(decreasing$summary_df, increasing$summary_df)

  write_table_tsv(decreasing_path, file.path(param_dir, "profile_path_decreasing.tsv"))
  write_table_tsv(increasing_path, file.path(param_dir, "profile_path_increasing.tsv"))
  write_table_tsv(combined_path, file.path(param_dir, "profile_path_combined.tsv"))
  write_table_tsv(combined_path, file.path(param_dir, "profile_likelihood.tsv"))
  write_table_tsv(seed_df, file.path(param_dir, "profile_seed_results.tsv"))
  write_table_tsv(relation_df, file.path(param_dir, "parameter_relation_long.tsv"))
  write_table_tsv(direction_summary_df, file.path(param_dir, "direction_summary.tsv"))

  write_profile_plots(
    combined_df = combined_path,
    param_dir = param_dir,
    param_symbol = param_symbol,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    baseline_value = baseline_value
  )

  message("Completed supplement-style profile path for parameter ", param_symbol)
  message("Parameter directory: ", normalizePath(param_dir, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  main()
}
