#!/usr/bin/env Rscript

# Usage:
#   Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/auto_calibrate_boundary_params.R
#   Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/auto_calibrate_boundary_params.R --config=/abs/path/to/O2G_supply_demand.yaml
#   Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/auto_calibrate_boundary_params.R --config=/abs/path/to/O2G_supply_demand.yaml --output_root=/abs/path/to/output_root
#   Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/auto_calibrate_boundary_params.R --config=/abs/path/to/O2G_supply_demand.yaml --parameter_table=/abs/path/to/parameter_table.csv --max_rounds_per_parameter=3 --seeds_per_round=10 --boundary_expand_fraction=0.10 --boundary_stick_lower_threshold=0.05 --boundary_stick_upper_threshold=0.95 --objective_threshold_total=9 --objective_threshold_burden=1.5 --objective_threshold_ploidy=7.5 --day1000_min_burden_threshold=2

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
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool
get_script_dir <- o2sd_get_script_dir

TARGET_PARAMETER_ORDER <- c(
  "rho_2N",
  "sigma_burden",
  "eta_o2",
  "mu_hp",
  "gamma_mu",
  "k_o_mis",
  "p_misseg",
  "p_wgd_max",
  "O2_wgd"
)

DEFAULT_MAX_ROUNDS_PER_PARAMETER <- 3L
DEFAULT_SEEDS_PER_ROUND <- 10L
DEFAULT_BOUNDARY_EXPAND_FRACTION <- 0.10
DEFAULT_BOUNDARY_STICK_LOWER_THRESHOLD <- 0.05
DEFAULT_BOUNDARY_STICK_UPPER_THRESHOLD <- 0.95
DEFAULT_OBJECTIVE_THRESHOLD_TOTAL <- 9.0
DEFAULT_OBJECTIVE_THRESHOLD_BURDEN <- 1.5
DEFAULT_OBJECTIVE_THRESHOLD_PLOIDY <- 7.5
DEFAULT_DAY1000_MIN_BURDEN_THRESHOLD <- 2.0
DEFAULT_PREDICTION_DAY <- 1000L
DEFAULT_SEED_BASE <- 100000L

default_config_path <- function(script_dir = SCRIPT_DIR) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  normalizePath(file.path(workflow_root, "..", "..", "config", "O2G_supply_demand.yaml"), mustWork = FALSE)
}

default_results_root <- function(script_dir = SCRIPT_DIR) {
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  normalizePath(file.path(workflow_root, "..", "..", "results"), mustWork = FALSE)
}

default_parameter_table_path <- function(script_dir = SCRIPT_DIR, glucose = TRUE) {
  default_o2g_parameter_table_path_common(
    script_dir = script_dir,
    glucose = glucose,
    must_exist = TRUE
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

read_yaml_config_file <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but is not installed.")
  }
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg)) cfg <- list()
  if (!is.list(cfg) || is.null(names(cfg))) {
    stop("Config file must contain a named YAML mapping: ", path)
  }
  cfg
}

resolve_parameter_table_path_from_config <- function(config_path, script_dir = SCRIPT_DIR) {
  cfg <- read_yaml_config_file(config_path)
  cfg_dir <- dirname(config_path)
  pt <- resolve_path_value(cfg$parameter_table, cfg_dir)
  if (is.null(pt)) {
    glucose_use <- canonical_glucose_enabled(
      .first_non_null_local(cfg$glucose, TRUE),
      default = TRUE
    )
    pt <- default_parameter_table_path(script_dir, glucose = glucose_use)
  }
  pt
}

resolve_runtime_settings <- function(argv) {
  settings <- list(
    max_rounds_per_parameter = as_int(argv$max_rounds_per_parameter, DEFAULT_MAX_ROUNDS_PER_PARAMETER),
    seeds_per_round = as_int(argv$seeds_per_round, DEFAULT_SEEDS_PER_ROUND),
    boundary_expand_fraction = as_num(argv$boundary_expand_fraction, DEFAULT_BOUNDARY_EXPAND_FRACTION),
    boundary_stick_lower_threshold = as_num(argv$boundary_stick_lower_threshold, DEFAULT_BOUNDARY_STICK_LOWER_THRESHOLD),
    boundary_stick_upper_threshold = as_num(argv$boundary_stick_upper_threshold, DEFAULT_BOUNDARY_STICK_UPPER_THRESHOLD),
    objective_threshold_total = as_num(argv$objective_threshold_total, DEFAULT_OBJECTIVE_THRESHOLD_TOTAL),
    objective_threshold_burden = as_num(argv$objective_threshold_burden, DEFAULT_OBJECTIVE_THRESHOLD_BURDEN),
    objective_threshold_ploidy = as_num(argv$objective_threshold_ploidy, DEFAULT_OBJECTIVE_THRESHOLD_PLOIDY),
    day1000_min_burden_threshold = as_num(argv$day1000_min_burden_threshold, DEFAULT_DAY1000_MIN_BURDEN_THRESHOLD),
    prediction_day = as_int(argv$prediction_day, DEFAULT_PREDICTION_DAY)
  )

  if (!is.finite(settings$max_rounds_per_parameter) || settings$max_rounds_per_parameter < 1) {
    stop("max_rounds_per_parameter must be a positive integer.")
  }
  if (!is.finite(settings$seeds_per_round) || settings$seeds_per_round < 1) {
    stop("seeds_per_round must be a positive integer.")
  }
  if (!is.finite(settings$boundary_expand_fraction) || settings$boundary_expand_fraction <= 0) {
    stop("boundary_expand_fraction must be > 0.")
  }
  if (!is.finite(settings$boundary_stick_lower_threshold) ||
      !is.finite(settings$boundary_stick_upper_threshold) ||
      settings$boundary_stick_lower_threshold < 0 ||
      settings$boundary_stick_upper_threshold > 1 ||
      settings$boundary_stick_lower_threshold >= settings$boundary_stick_upper_threshold) {
    stop("boundary stick thresholds must satisfy 0 <= lower < upper <= 1.")
  }
  if (!is.finite(settings$objective_threshold_total) ||
      !is.finite(settings$objective_threshold_burden) ||
      !is.finite(settings$objective_threshold_ploidy)) {
    stop("Objective thresholds must be finite.")
  }
  if (!is.finite(settings$day1000_min_burden_threshold)) {
    stop("day1000_min_burden_threshold must be finite.")
  }
  if (!is.finite(settings$prediction_day) || settings$prediction_day < 1) {
    stop("prediction_day must be a positive integer.")
  }

  settings$max_rounds_per_parameter <- as.integer(settings$max_rounds_per_parameter)
  settings$seeds_per_round <- as.integer(settings$seeds_per_round)
  settings$prediction_day <- as.integer(settings$prediction_day)
  settings$seed_block_width <- as.integer(
    max(
      1000L,
      settings$seeds_per_round * settings$max_rounds_per_parameter + 100L
    )
  )
  settings
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
  required <- c("param_symbol", "estimate", "init_value", "lower_bound", "upper_bound", "description")
  missing_cols <- setdiff(required, names(tab))
  if (length(missing_cols) > 0L) {
    stop("Missing required parameter-table columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }
  if (!("source" %in% names(tab))) {
    tab$source <- ""
  }
  tab$param_symbol <- trimws(as.character(tab$param_symbol))
  if (any(!nzchar(tab$param_symbol))) {
    stop("Blank param_symbol values were found in ", path)
  }
  if (any(duplicated(tab$param_symbol))) {
    dup <- unique(tab$param_symbol[duplicated(tab$param_symbol)])
    stop("Duplicated param_symbol values were found in ", path, ": ", paste(dup, collapse = ", "))
  }
  tab$estimate <- normalize_estimate_column(tab$estimate)
  numeric_cols <- c("init_value", "lower_bound", "upper_bound")
  for (col in numeric_cols) {
    tab[[col]] <- suppressWarnings(as.numeric(tab[[col]]))
    if (any(!is.finite(tab[[col]]))) {
      bad_rows <- which(!is.finite(tab[[col]]))
      stop("Non-finite values were found in column ", col, " at rows: ", paste(bad_rows, collapse = ", "))
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

find_parameter_index <- function(tab, param_name) {
  idx <- which(as.character(tab$param_symbol) == param_name)
  if (length(idx) != 1L) {
    stop("Parameter ", param_name, " was not found exactly once in the parameter table.")
  }
  idx[[1]]
}

validate_target_parameter_row <- function(tab, param_name) {
  idx <- find_parameter_index(tab, param_name)
  lower <- suppressWarnings(as.numeric(tab$lower_bound[[idx]]))
  upper <- suppressWarnings(as.numeric(tab$upper_bound[[idx]]))
  init <- suppressWarnings(as.numeric(tab$init_value[[idx]]))
  if (!is.finite(lower) || !is.finite(upper) || !is.finite(init)) {
    stop("Non-finite init/bounds were found for parameter ", param_name)
  }
  if (lower < 0 || upper < 0) {
    stop("Negative bounds were found for nonnegative calibration parameter ", param_name)
  }
  if (upper <= lower) {
    stop("Invalid bounds were found for parameter ", param_name, ": upper_bound must be > lower_bound.")
  }
  invisible(idx)
}

build_round_parameter_table <- function(working_table, target_param) {
  idx <- find_parameter_index(working_table, target_param)
  validate_target_parameter_row(working_table, target_param)
  round_table <- working_table
  round_table$estimate <- rep(FALSE, nrow(round_table))
  fixed_idx <- setdiff(seq_len(nrow(round_table)), idx)
  if (length(fixed_idx) > 0L) {
    round_table$lower_bound[fixed_idx] <- round_table$init_value[fixed_idx]
    round_table$upper_bound[fixed_idx] <- round_table$init_value[fixed_idx]
  }
  round_table$estimate[[idx]] <- TRUE
  round_table
}

generate_deterministic_seeds <- function(param_index, round_index, settings) {
  base <- as.integer(
    DEFAULT_SEED_BASE +
      as.integer(param_index) * as.integer(settings$seed_block_width) +
      (as.integer(round_index) - 1L) * as.integer(settings$seeds_per_round)
  )
  as.integer(base + seq_len(as.integer(settings$seeds_per_round)))
}

write_seed_file <- function(seeds, path) {
  writeLines(as.character(as.integer(seeds)), con = path, sep = "\n")
  invisible(path)
}

build_round_paths <- function(output_root, param_index, param_name, round_index) {
  param_dir <- file.path(output_root, sprintf("%02d_%s", as.integer(param_index), param_name))
  round_dir <- file.path(param_dir, sprintf("round_%02d", as.integer(round_index)))
  list(
    param_dir = param_dir,
    round_dir = round_dir,
    temp_parameter_table = file.path(round_dir, "parameter_table.round_input.csv"),
    temp_config = file.path(round_dir, "config.round.yaml"),
    temp_seed_file = file.path(round_dir, "seeds.round.csv")
  )
}

write_round_config <- function(base_config_path, out_path, round_dir, temp_parameter_table, temp_seed_file) {
  cfg <- read_yaml_config_file(base_config_path)
  cfg_dir <- dirname(base_config_path)
  data_dir <- resolve_path_value(cfg$data_dir, cfg_dir)
  if (is.null(data_dir)) {
    stop("data_dir was not found in config: ", base_config_path)
  }
  cfg$run_prefix <- basename(round_dir)
  cfg$out_root <- normalizePath(dirname(round_dir), mustWork = FALSE)
  cfg$data_dir <- data_dir
  cfg$seeds_file <- normalizePath(temp_seed_file, mustWork = FALSE)
  cfg$seeds_csv <- ""
  cfg$parameter_table <- normalizePath(temp_parameter_table, mustWork = FALSE)
  cfg$append_run_prefix_timestamp <- FALSE
  cfg$auto_viz <- TRUE
  yaml::write_yaml(cfg, out_path)
  invisible(out_path)
}

launch_fitting_batch <- function(fit_script, config_path) {
  if (!file.exists(fit_script)) {
    stop("Fitter script was not found: ", fit_script)
  }
  status <- system2(
    "Rscript",
    args = c(normalizePath(fit_script, mustWork = FALSE), "--fit_invivo", "--mode=run", paste0("--config=", normalizePath(config_path, mustWork = FALSE))),
    stdout = "",
    stderr = "",
    wait = TRUE
  )
  if (!is.null(status) && status != 0L) {
    stop("Fitting batch failed with exit status ", status, " for config ", config_path)
  }
  invisible(TRUE)
}

launch_extra_results <- function(extra_results_script, run_dir) {
  if (!file.exists(extra_results_script)) {
    stop("extra_results.R was not found: ", extra_results_script)
  }
  status <- system2(
    "Rscript",
    args = c(normalizePath(extra_results_script, mustWork = FALSE), paste0("--run_dir=", normalizePath(run_dir, mustWork = FALSE))),
    stdout = "",
    stderr = "",
    wait = TRUE
  )
  if (!is.null(status) && status != 0L) {
    stop("extra_results.R failed with exit status ", status, " for run_dir ", run_dir)
  }
  out_path <- file.path(run_dir, "extra_results", "seed_summary.tsv")
  if (!file.exists(out_path)) {
    stop("seed_summary.tsv was not produced by extra_results.R: ", out_path)
  }
  out_path
}

seed_numeric_key <- function(seed_ids) {
  out <- rep(Inf, length(seed_ids))
  m <- regexec("^seed([0-9]+)$", seed_ids)
  hit <- regmatches(seed_ids, m)
  for (i in seq_along(hit)) {
    if (length(hit[[i]]) == 2L) {
      out[[i]] <- suppressWarnings(as.numeric(hit[[i]][[2]]))
    }
  }
  out
}

read_seed_summary <- function(path) {
  if (!file.exists(path)) {
    stop("seed_summary.tsv was not found: ", path)
  }
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("seed", "objective", "objective_burden", "objective_ploidy")
  missing_cols <- setdiff(required, names(tab))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in seed_summary.tsv: ", paste(missing_cols, collapse = ", "))
  }
  if (!nrow(tab)) {
    stop("seed_summary.tsv contains no complete seeds: ", path)
  }
  tab$seed <- as.character(tab$seed)
  numeric_cols <- intersect(
    c(
      "objective",
      "objective_burden",
      "objective_ploidy",
      "recommend_rank_burden_ploidy_boundary_first",
      "recommend_rank_ploidy_burden_boundary_first",
      "recommend_rank_ploidy_boundary_first",
      "recommend_rank_ploidy_first",
      "boundary_rank_active_support",
      "objective_rank"
    ),
    names(tab)
  )
  for (col in numeric_cols) {
    tab[[col]] <- suppressWarnings(as.numeric(tab[[col]]))
  }
  rank_candidates <- c(
    "recommend_rank_burden_ploidy_boundary_first",
    "recommend_rank_ploidy_burden_boundary_first",
    "recommend_rank_ploidy_boundary_first",
    "recommend_rank_ploidy_first",
    "boundary_rank_active_support",
    "objective_rank"
  )
  rank_matches <- rank_candidates[rank_candidates %in% names(tab)]
  if (length(rank_matches) == 0L) {
    stop("No usable rank column was found in seed_summary.tsv.")
  }
  rank_col <- rank_matches[[1]]
  list(summary = tab, rank_col = rank_col)
}

read_named_value_table <- function(path, key_col, value_col) {
  if (!file.exists(path)) {
    stop("Required table was not found: ", path)
  }
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c(key_col, value_col) %in% names(tab))) {
    stop("Missing required columns in ", path, ": ", key_col, ", ", value_col)
  }
  setNames(tab[[value_col]], tab[[key_col]])
}

read_best_param_value <- function(seed_dir, param_name) {
  path <- file.path(seed_dir, "best_params.tsv")
  vals <- read_named_value_table(path, "parameter", "value")
  if (!(param_name %in% names(vals))) {
    stop("Parameter ", param_name, " was not found in best_params.tsv: ", path)
  }
  value <- suppressWarnings(as.numeric(vals[[param_name]]))
  if (!is.finite(value)) {
    stop("Non-finite fitted value was found for parameter ", param_name, " in ", path)
  }
  value
}

safe_read_best_param_value <- function(seed_dir, param_name) {
  tryCatch(
    {
      list(ok = TRUE, value = read_best_param_value(seed_dir, param_name), message = "")
    },
    error = function(e) {
      list(ok = FALSE, value = NA_real_, message = conditionMessage(e))
    }
  )
}

read_seed_prediction_output <- function(seed_dir, settings) {
  seed_id <- basename(seed_dir)
  path <- file.path(seed_dir, "viz", paste0("predict_burden_0_", as.integer(settings$prediction_day), "day.tsv"))
  out <- data.frame(
    seed = seed_id,
    prediction_file_present = FALSE,
    min_pred_burden_2N_day1000 = NA_real_,
    min_pred_burden_4N_day1000 = NA_real_,
    prediction_pass = FALSE,
    stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    return(out)
  }
  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab)) {
    return(out)
  }
  required <- c("day", "cohort", "pred_burden")
  if (!all(required %in% names(tab))) {
    return(out)
  }
  tab$day <- suppressWarnings(as.numeric(tab$day))
  tab$pred_burden <- suppressWarnings(as.numeric(tab$pred_burden))
  day_keep <- is.finite(tab$day) & abs(tab$day - as.numeric(settings$prediction_day)) <= 1e-9
  tab <- tab[day_keep, , drop = FALSE]
  if (nrow(tab) == 0L) {
    out$prediction_file_present <- TRUE
    return(out)
  }
  vals_2n <- tab$pred_burden[tab$cohort == "2N" & is.finite(tab$pred_burden)]
  vals_4n <- tab$pred_burden[tab$cohort == "4N" & is.finite(tab$pred_burden)]
  min_2n <- if (length(vals_2n) > 0L) min(vals_2n) else NA_real_
  min_4n <- if (length(vals_4n) > 0L) min(vals_4n) else NA_real_
  out$prediction_file_present <- TRUE
  out$min_pred_burden_2N_day1000 <- min_2n
  out$min_pred_burden_4N_day1000 <- min_4n
  out$prediction_pass <- is.finite(min_2n) && is.finite(min_4n) &&
    (min_2n > settings$day1000_min_burden_threshold) &&
    (min_4n > settings$day1000_min_burden_threshold)
  out
}

read_prediction_batch <- function(run_dir, seeds, settings) {
  rows <- lapply(
    as.integer(seeds),
    function(seed) {
      read_seed_prediction_output(file.path(run_dir, paste0("seed", seed)), settings = settings)
    }
  )
  do.call(rbind, rows)
}

select_reference_seed <- function(seed_summary, rank_col, prediction_info, settings) {
  if (!nrow(seed_summary)) {
    stop("No complete seeds were found in seed_summary.tsv.")
  }
  merged <- merge(seed_summary, prediction_info, by = "seed", all.x = TRUE, sort = FALSE)
  merged$rank_value <- suppressWarnings(as.numeric(merged[[rank_col]]))
  merged$objective <- suppressWarnings(as.numeric(merged$objective))
  merged$objective_burden <- suppressWarnings(as.numeric(merged$objective_burden))
  merged$objective_ploidy <- suppressWarnings(as.numeric(merged$objective_ploidy))
  merged$prediction_file_present[is.na(merged$prediction_file_present)] <- FALSE
  merged$prediction_pass[is.na(merged$prediction_pass)] <- FALSE
  merged$passes_objective_thresholds <- (
    is.finite(merged$objective) &
      is.finite(merged$objective_burden) &
      is.finite(merged$objective_ploidy) &
      merged$objective <= settings$objective_threshold_total &
      merged$objective_burden <= settings$objective_threshold_burden &
      merged$objective_ploidy <= settings$objective_threshold_ploidy
  )
  merged$eligible <- is.finite(merged$rank_value) &
    merged$passes_objective_thresholds &
    merged$prediction_pass

  ord <- order(merged$rank_value, merged$objective, seed_numeric_key(merged$seed), na.last = TRUE)
  merged <- merged[ord, , drop = FALSE]
  row.names(merged) <- NULL

  if (!any(is.finite(merged$rank_value))) {
    stop("No finite values were found in rank column ", rank_col, ".")
  }

  diagnostic_idx <- which(is.finite(merged$rank_value))[[1]]
  diagnostic_row <- merged[diagnostic_idx, , drop = FALSE]

  eligible_rows <- which(merged$eligible)
  reference_row <- if (length(eligible_rows) > 0L) merged[eligible_rows[[1]], , drop = FALSE] else NULL

  list(
    merged = merged,
    diagnostic = diagnostic_row,
    reference = reference_row,
    valid_reference_found = !is.null(reference_row),
    n_complete_seeds = nrow(seed_summary),
    n_eligible_seeds = sum(merged$eligible, na.rm = TRUE),
    n_prediction_pass_seeds = sum(merged$prediction_pass, na.rm = TRUE)
  )
}

check_boundary_stickiness <- function(value, lower_bound, upper_bound, settings) {
  if (!is.finite(value)) {
    stop("Selected parameter value is not finite.")
  }
  if (!is.finite(lower_bound) || !is.finite(upper_bound) || upper_bound <= lower_bound) {
    stop("Invalid bounds were supplied for boundary-stickiness checking.")
  }
  rel <- (value - lower_bound) / (upper_bound - lower_bound)
  near_lower <- is.finite(rel) && (rel <= settings$boundary_stick_lower_threshold)
  near_upper <- is.finite(rel) && (rel >= settings$boundary_stick_upper_threshold)
  status <- if (near_lower && !near_upper) {
    "near_lower"
  } else if (near_upper && !near_lower) {
    "near_upper"
  } else if (near_lower && near_upper) {
    "near_both"
  } else {
    "interior"
  }
  list(
    rel = rel,
    near_lower = near_lower,
    near_upper = near_upper,
    sticky = near_lower || near_upper,
    status = status
  )
}

expand_bounds <- function(param_name, lower_bound, upper_bound, boundary_info, settings) {
  width <- upper_bound - lower_bound
  if (!is.finite(width) || width <= 0) {
    stop("A non-positive bound width was found for parameter ", param_name)
  }
  lower_new <- lower_bound
  upper_new <- upper_bound
  lower_blocked <- FALSE
  expanded <- FALSE
  expansion_side <- "none"

  if (isTRUE(boundary_info$near_upper)) {
    upper_new <- upper_bound + settings$boundary_expand_fraction * width
    expanded <- TRUE
    expansion_side <- "upper"
  } else if (isTRUE(boundary_info$near_lower)) {
    candidate <- lower_bound - settings$boundary_expand_fraction * width
    expansion_side <- "lower"
    if (candidate < 0) {
      lower_blocked <- TRUE
      lower_new <- lower_bound
    } else {
      lower_new <- candidate
      expanded <- TRUE
    }
  }

  list(
    lower_new = lower_new,
    upper_new = upper_new,
    lower_blocked_by_zero_floor = lower_blocked,
    expanded = expanded,
    expansion_side = expansion_side
  )
}

append_round_logs <- function(round_log_df, out_path) {
  utils::write.table(
    round_log_df,
    file = out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(out_path)
}

write_parameter_summary <- function(parameter_summary_df, out_path) {
  utils::write.table(
    parameter_summary_df,
    file = out_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(out_path)
}

write_table_snapshot <- function(tab, path) {
  write_natural_parameter_table(tab, path)
  invisible(path)
}

write_rollback_record <- function(
  out_path,
  param_name,
  param_index,
  rounds_attempted,
  failure_reason,
  pre_parameter_snapshot_path,
  rollback_snapshot_path,
  last_round_dir = NA_character_
) {
  lines <- c(
    paste0("parameter\t", param_name),
    paste0("parameter_index\t", as.integer(param_index)),
    paste0("rounds_attempted\t", as.integer(rounds_attempted)),
    paste0("failure_reason\t", failure_reason),
    paste0("pre_parameter_snapshot_path\t", pre_parameter_snapshot_path),
    paste0("rollback_snapshot_path\t", rollback_snapshot_path),
    paste0("last_round_dir\t", ifelse(is.na(last_round_dir), "", last_round_dir))
  )
  writeLines(lines, con = out_path, sep = "\n")
  invisible(out_path)
}

escape_markdown_cell <- function(x) {
  txt <- ifelse(is.na(x), "", as.character(x))
  txt <- gsub("\\|", "\\\\|", txt)
  txt <- gsub("\n", " ", txt, fixed = TRUE)
  txt
}

data_frame_to_markdown <- function(tab) {
  if (!ncol(tab)) {
    return("No rows were produced.\n")
  }
  cols <- names(tab)
  header <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
  if (!nrow(tab)) {
    return(paste(c(header, sep), collapse = "\n"))
  }
  rows <- apply(
    tab,
    1,
    function(row) {
      paste0("| ", paste(vapply(row, escape_markdown_cell, character(1)), collapse = " | "), " |")
    }
  )
  paste(c(header, sep, rows), collapse = "\n")
}

settings_to_data_frame <- function(settings) {
  data.frame(
    setting = c(
      "max_rounds_per_parameter",
      "seeds_per_round",
      "boundary_expand_fraction",
      "boundary_stick_lower_threshold",
      "boundary_stick_upper_threshold",
      "objective_threshold_total",
      "objective_threshold_burden",
      "objective_threshold_ploidy",
      "day1000_min_burden_threshold",
      "prediction_day",
      "seed_block_width"
    ),
    value = c(
      as.character(settings$max_rounds_per_parameter),
      as.character(settings$seeds_per_round),
      as.character(settings$boundary_expand_fraction),
      as.character(settings$boundary_stick_lower_threshold),
      as.character(settings$boundary_stick_upper_threshold),
      as.character(settings$objective_threshold_total),
      as.character(settings$objective_threshold_burden),
      as.character(settings$objective_threshold_ploidy),
      as.character(settings$day1000_min_burden_threshold),
      as.character(settings$prediction_day),
      as.character(settings$seed_block_width)
    ),
    stringsAsFactors = FALSE
  )
}

write_final_report <- function(
  out_path,
  config_path,
  parameter_table_path,
  output_root,
  settings,
  round_log_df,
  parameter_summary_df,
  final_parameter_table
) {
  calibration_order_df <- data.frame(
    order = seq_along(TARGET_PARAMETER_ORDER),
    parameter = TARGET_PARAMETER_ORDER,
    stringsAsFactors = FALSE
  )
  final_table_md <- final_parameter_table[, c("param_symbol", "init_value", "lower_bound", "upper_bound", "estimate"), drop = FALSE]
  final_table_md$estimate <- ifelse(final_table_md$estimate, "TRUE", "FALSE")

  report_lines <- c(
    "# Auto Calibration Report",
    "",
    paste0("- Config: `", normalizePath(config_path, mustWork = FALSE), "`"),
    paste0("- Parameter table input: `", normalizePath(parameter_table_path, mustWork = FALSE), "`"),
    paste0("- Output root: `", normalizePath(output_root, mustWork = FALSE), "`"),
    "- Failed parameters were reverted to the pre-parameter snapshot and were not propagated downstream.",
    "",
    "## Runtime Settings",
    "",
    data_frame_to_markdown(settings_to_data_frame(settings)),
    "",
    "## Calibration Order",
    "",
    data_frame_to_markdown(calibration_order_df),
    "",
    "## Parameter Summary",
    "",
    data_frame_to_markdown(parameter_summary_df),
    "",
    "## Round Log",
    "",
    data_frame_to_markdown(round_log_df),
    "",
    "## Final Parameter Table",
    "",
    data_frame_to_markdown(final_table_md),
    ""
  )
  writeLines(report_lines, con = out_path, sep = "\n")
  invisible(out_path)
}

make_empty_round_log <- function() {
  data.frame(
    parameter = character(0),
    parameter_index = integer(0),
    round_index = integer(0),
    configured_max_rounds_per_parameter = integer(0),
    configured_seeds_per_round = integer(0),
    seeds_used = character(0),
    run_dir = character(0),
    rank_column = character(0),
    n_complete_seeds = integer(0),
    n_prediction_pass_seeds = integer(0),
    n_eligible_seeds = integer(0),
    valid_eligible_reference_seed_found = logical(0),
    eligible_reference_seed = character(0),
    eligible_reference_rank = numeric(0),
    eligible_reference_value = numeric(0),
    eligible_reference_objective = numeric(0),
    eligible_reference_objective_burden = numeric(0),
    eligible_reference_objective_ploidy = numeric(0),
    eligible_reference_min_pred_burden_2N_day1000 = numeric(0),
    eligible_reference_min_pred_burden_4N_day1000 = numeric(0),
    eligible_reference_prediction_file_present = logical(0),
    eligible_reference_boundary_rel = numeric(0),
    eligible_reference_boundary_status = character(0),
    eligible_reference_boundary_sticky = logical(0),
    best_ranked_diagnostic_seed = character(0),
    best_ranked_diagnostic_rank = numeric(0),
    best_ranked_diagnostic_value = numeric(0),
    best_ranked_diagnostic_objective = numeric(0),
    best_ranked_diagnostic_objective_burden = numeric(0),
    best_ranked_diagnostic_objective_ploidy = numeric(0),
    best_ranked_diagnostic_min_pred_burden_2N_day1000 = numeric(0),
    best_ranked_diagnostic_min_pred_burden_4N_day1000 = numeric(0),
    best_ranked_diagnostic_prediction_file_present = logical(0),
    old_lower_bound = numeric(0),
    old_upper_bound = numeric(0),
    new_lower_bound = numeric(0),
    new_upper_bound = numeric(0),
    lower_blocked_by_zero_floor = logical(0),
    expansion_side = character(0),
    round_outcome = character(0),
    parameter_completed_successfully = logical(0),
    parameter_failed = logical(0),
    rollback_applied_after_round = logical(0),
    failure_reason_if_any = character(0),
    stringsAsFactors = FALSE
  )
}

make_empty_parameter_summary <- function() {
  data.frame(
    parameter = character(0),
    parameter_index = integer(0),
    configured_max_rounds_per_parameter = integer(0),
    configured_seeds_per_round = integer(0),
    rounds_attempted = integer(0),
    status = character(0),
    completed_successfully = logical(0),
    failed = logical(0),
    rollback_applied = logical(0),
    failure_reason = character(0),
    valid_reference_seed_found_in_any_round = logical(0),
    successful_reference_seed = character(0),
    successful_reference_rank = numeric(0),
    successful_selected_value = numeric(0),
    successful_selected_objective = numeric(0),
    successful_selected_objective_burden = numeric(0),
    successful_selected_objective_ploidy = numeric(0),
    successful_min_pred_burden_2N_day1000 = numeric(0),
    successful_min_pred_burden_4N_day1000 = numeric(0),
    successful_boundary_rel = numeric(0),
    successful_boundary_status = character(0),
    final_working_init_value = numeric(0),
    final_working_lower_bound = numeric(0),
    final_working_upper_bound = numeric(0),
    pre_parameter_snapshot_path = character(0),
    final_snapshot_path = character(0),
    rollback_record_path = character(0),
    final_run_dir = character(0),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  script_dir <- SCRIPT_DIR
  argv <- parse_args(commandArgs(trailingOnly = TRUE))
  settings <- resolve_runtime_settings(argv)

  fit_script <- normalizePath(file.path(script_dir, "fit_model_O2G_supply_demand_MAP.R"), mustWork = TRUE)
  extra_results_script <- normalizePath(file.path(WORKFLOW_ROOT, "analysis", "extra_results.R"), mustWork = TRUE)

  config_path <- resolve_path_value(argv$config, getwd()) %||% default_config_path(script_dir)
  config_path <- normalizePath(config_path, mustWork = TRUE)

  parameter_table_path <- resolve_path_value(argv$parameter_table, getwd())
  if (is.null(parameter_table_path)) {
    parameter_table_path <- resolve_parameter_table_path_from_config(config_path, script_dir = script_dir)
  }
  parameter_table_path <- normalizePath(parameter_table_path, mustWork = TRUE)

  output_root <- resolve_path_value(argv$output_root, getwd())
  if (is.null(output_root)) {
    output_root <- file.path(
      default_results_root(script_dir),
      paste0("auto_calibrate_boundary_params_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
  }
  output_root <- normalizePath(output_root, mustWork = FALSE)
  if (dir.exists(output_root) && length(list.files(output_root, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("Output root already exists and is not empty: ", output_root)
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  snapshot_dir <- file.path(output_root, "parameter_table_snapshots")
  rollback_dir <- file.path(output_root, "rollback_records")
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rollback_dir, recursive = TRUE, showWarnings = FALSE)

  backup_path <- file.path(output_root, "parameter_table.original.backup.csv")
  ok_backup <- file.copy(parameter_table_path, backup_path, overwrite = TRUE)
  if (!isTRUE(ok_backup)) {
    stop("Original parameter table could not be backed up to ", backup_path)
  }

  working_table <- read_natural_parameter_table(parameter_table_path)
  original_estimate <- working_table$estimate

  round_log_path <- file.path(output_root, "auto_calibration_round_log.tsv")
  parameter_summary_path <- file.path(output_root, "auto_calibration_parameter_summary.tsv")
  report_path <- file.path(output_root, "auto_calibration_report.md")
  updated_parameter_table_path <- file.path(output_root, "parameter_table.updated.csv")

  round_log_df <- make_empty_round_log()
  parameter_summary_df <- make_empty_parameter_summary()

  write_table_snapshot(working_table, file.path(snapshot_dir, "parameter_table_initial.csv"))
  append_round_logs(round_log_df, round_log_path)
  write_parameter_summary(parameter_summary_df, parameter_summary_path)

  for (param_index in seq_along(TARGET_PARAMETER_ORDER)) {
    param_name <- TARGET_PARAMETER_ORDER[[param_index]]
    validate_target_parameter_row(working_table, param_name)

    pre_parameter_table <- working_table
    pre_parameter_snapshot_path <- file.path(
      snapshot_dir,
      sprintf("parameter_table_before_parameter_%02d_%s.csv", as.integer(param_index), param_name)
    )
    write_table_snapshot(pre_parameter_table, pre_parameter_snapshot_path)

    parameter_table_state <- pre_parameter_table
    parameter_table_state$estimate <- original_estimate

    rounds_attempted <- 0L
    completed_successfully <- FALSE
    failed <- FALSE
    rollback_applied <- FALSE
    failure_reason <- ""
    valid_reference_seed_found_in_any_round <- FALSE
    successful_round_record <- NULL
    last_round_record <- NULL
    final_snapshot_path <- NA_character_
    rollback_record_path <- NA_character_

    for (round_index in seq_len(settings$max_rounds_per_parameter)) {
      rounds_attempted <- as.integer(round_index)
      round_paths <- build_round_paths(output_root, param_index, param_name, round_index)
      dir.create(round_paths$param_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(round_paths$round_dir, recursive = TRUE, showWarnings = FALSE)

      seeds <- generate_deterministic_seeds(param_index = param_index, round_index = round_index, settings = settings)
      round_table <- build_round_parameter_table(parameter_table_state, param_name)
      write_natural_parameter_table(round_table, round_paths$temp_parameter_table)
      write_table_snapshot(
        round_table,
        file.path(
          snapshot_dir,
          sprintf("parameter_table_round_input_%02d_%s_round_%02d.csv", as.integer(param_index), param_name, as.integer(round_index))
        )
      )
      write_seed_file(seeds, round_paths$temp_seed_file)
      write_round_config(
        base_config_path = config_path,
        out_path = round_paths$temp_config,
        round_dir = round_paths$round_dir,
        temp_parameter_table = round_paths$temp_parameter_table,
        temp_seed_file = round_paths$temp_seed_file
      )

      launch_fitting_batch(fit_script = fit_script, config_path = round_paths$temp_config)
      seed_summary_path <- launch_extra_results(extra_results_script = extra_results_script, run_dir = round_paths$round_dir)
      summary_info <- read_seed_summary(seed_summary_path)
      prediction_info <- read_prediction_batch(round_paths$round_dir, seeds, settings = settings)
      selection <- select_reference_seed(
        seed_summary = summary_info$summary,
        rank_col = summary_info$rank_col,
        prediction_info = prediction_info,
        settings = settings
      )

      row_idx <- find_parameter_index(parameter_table_state, param_name)
      old_lower <- as.numeric(parameter_table_state$lower_bound[[row_idx]])
      old_upper <- as.numeric(parameter_table_state$upper_bound[[row_idx]])
      new_lower <- old_lower
      new_upper <- old_upper
      lower_blocked_by_zero_floor <- FALSE
      expansion_side <- "none"

      diagnostic <- selection$diagnostic
      diagnostic_seed <- as.character(diagnostic$seed[[1]])
      diagnostic_seed_dir <- file.path(round_paths$round_dir, diagnostic_seed)
      diagnostic_value_res <- safe_read_best_param_value(diagnostic_seed_dir, param_name)
      diagnostic_value <- if (isTRUE(diagnostic_value_res$ok)) diagnostic_value_res$value else NA_real_

      reference_row <- selection$reference
      valid_reference_seed_found <- isTRUE(selection$valid_reference_found)
      eligible_reference_seed <- if (valid_reference_seed_found) as.character(reference_row$seed[[1]]) else NA_character_
      eligible_reference_rank <- if (valid_reference_seed_found) suppressWarnings(as.numeric(reference_row$rank_value[[1]])) else NA_real_
      eligible_reference_value <- NA_real_
      eligible_reference_boundary_rel <- NA_real_
      eligible_reference_boundary_status <- NA_character_
      eligible_reference_boundary_sticky <- NA
      eligible_reference_prediction_file_present <- if (valid_reference_seed_found) isTRUE(reference_row$prediction_file_present[[1]]) else NA

      if (valid_reference_seed_found) {
        valid_reference_seed_found_in_any_round <- TRUE
      }

      round_outcome <- "no_valid_eligible_reference_seed"

      if (valid_reference_seed_found) {
        eligible_seed_dir <- file.path(round_paths$round_dir, eligible_reference_seed)
        eligible_value_res <- safe_read_best_param_value(eligible_seed_dir, param_name)
        if (!isTRUE(eligible_value_res$ok)) {
          failed <- TRUE
          failure_reason <- paste0("eligible_reference_value_unavailable: ", eligible_value_res$message)
          round_outcome <- "eligible_reference_value_unavailable"
        } else {
          eligible_reference_value <- eligible_value_res$value
          eligible_boundary <- check_boundary_stickiness(
            value = eligible_reference_value,
            lower_bound = old_lower,
            upper_bound = old_upper,
            settings = settings
          )
          eligible_reference_boundary_rel <- as.numeric(eligible_boundary$rel)
          eligible_reference_boundary_status <- as.character(eligible_boundary$status)
          eligible_reference_boundary_sticky <- isTRUE(eligible_boundary$sticky)

          if (!isTRUE(eligible_boundary$sticky)) {
            parameter_table_state$init_value[[row_idx]] <- eligible_reference_value
            parameter_table_state$estimate <- original_estimate
            completed_successfully <- TRUE
            round_outcome <- "successful_eligible_interior_reference"
          } else {
            expansion <- expand_bounds(
              param_name = param_name,
              lower_bound = old_lower,
              upper_bound = old_upper,
              boundary_info = eligible_boundary,
              settings = settings
            )
            parameter_table_state$init_value[[row_idx]] <- eligible_reference_value
            parameter_table_state$lower_bound[[row_idx]] <- expansion$lower_new
            parameter_table_state$upper_bound[[row_idx]] <- expansion$upper_new
            parameter_table_state$estimate <- original_estimate

            new_lower <- as.numeric(expansion$lower_new)
            new_upper <- as.numeric(expansion$upper_new)
            lower_blocked_by_zero_floor <- isTRUE(expansion$lower_blocked_by_zero_floor)
            expansion_side <- as.character(expansion$expansion_side)
            round_outcome <- if (isTRUE(expansion$expanded)) {
              "eligible_reference_boundary_expanded"
            } else {
              "eligible_reference_boundary_zero_floor_blocked"
            }
          }
        }
      } else {
        parameter_table_state$estimate <- original_estimate
      }

      if (!completed_successfully && !failed && round_index >= settings$max_rounds_per_parameter) {
        failed <- TRUE
        if (!isTRUE(valid_reference_seed_found_in_any_round)) {
          failure_reason <- paste0(
            "no_valid_eligible_reference_seed_within_",
            as.integer(settings$max_rounds_per_parameter),
            "_rounds"
          )
        } else {
          failure_reason <- paste0(
            "no_valid_eligible_non_boundary_sticking_result_within_",
            as.integer(settings$max_rounds_per_parameter),
            "_rounds"
          )
        }
      }

      round_record <- data.frame(
        parameter = param_name,
        parameter_index = as.integer(param_index),
        round_index = as.integer(round_index),
        configured_max_rounds_per_parameter = as.integer(settings$max_rounds_per_parameter),
        configured_seeds_per_round = as.integer(settings$seeds_per_round),
        seeds_used = paste(as.integer(seeds), collapse = ","),
        run_dir = normalizePath(round_paths$round_dir, mustWork = FALSE),
        rank_column = summary_info$rank_col,
        n_complete_seeds = as.integer(selection$n_complete_seeds),
        n_prediction_pass_seeds = as.integer(selection$n_prediction_pass_seeds),
        n_eligible_seeds = as.integer(selection$n_eligible_seeds),
        valid_eligible_reference_seed_found = isTRUE(valid_reference_seed_found),
        eligible_reference_seed = eligible_reference_seed,
        eligible_reference_rank = eligible_reference_rank,
        eligible_reference_value = eligible_reference_value,
        eligible_reference_objective = if (valid_reference_seed_found) suppressWarnings(as.numeric(reference_row$objective[[1]])) else NA_real_,
        eligible_reference_objective_burden = if (valid_reference_seed_found) suppressWarnings(as.numeric(reference_row$objective_burden[[1]])) else NA_real_,
        eligible_reference_objective_ploidy = if (valid_reference_seed_found) suppressWarnings(as.numeric(reference_row$objective_ploidy[[1]])) else NA_real_,
        eligible_reference_min_pred_burden_2N_day1000 = if (valid_reference_seed_found) suppressWarnings(as.numeric(reference_row$min_pred_burden_2N_day1000[[1]])) else NA_real_,
        eligible_reference_min_pred_burden_4N_day1000 = if (valid_reference_seed_found) suppressWarnings(as.numeric(reference_row$min_pred_burden_4N_day1000[[1]])) else NA_real_,
        eligible_reference_prediction_file_present = eligible_reference_prediction_file_present,
        eligible_reference_boundary_rel = eligible_reference_boundary_rel,
        eligible_reference_boundary_status = eligible_reference_boundary_status,
        eligible_reference_boundary_sticky = eligible_reference_boundary_sticky,
        best_ranked_diagnostic_seed = diagnostic_seed,
        best_ranked_diagnostic_rank = suppressWarnings(as.numeric(diagnostic$rank_value[[1]])),
        best_ranked_diagnostic_value = diagnostic_value,
        best_ranked_diagnostic_objective = suppressWarnings(as.numeric(diagnostic$objective[[1]])),
        best_ranked_diagnostic_objective_burden = suppressWarnings(as.numeric(diagnostic$objective_burden[[1]])),
        best_ranked_diagnostic_objective_ploidy = suppressWarnings(as.numeric(diagnostic$objective_ploidy[[1]])),
        best_ranked_diagnostic_min_pred_burden_2N_day1000 = suppressWarnings(as.numeric(diagnostic$min_pred_burden_2N_day1000[[1]])),
        best_ranked_diagnostic_min_pred_burden_4N_day1000 = suppressWarnings(as.numeric(diagnostic$min_pred_burden_4N_day1000[[1]])),
        best_ranked_diagnostic_prediction_file_present = isTRUE(diagnostic$prediction_file_present[[1]]),
        old_lower_bound = old_lower,
        old_upper_bound = old_upper,
        new_lower_bound = new_lower,
        new_upper_bound = new_upper,
        lower_blocked_by_zero_floor = isTRUE(lower_blocked_by_zero_floor),
        expansion_side = expansion_side,
        round_outcome = round_outcome,
        parameter_completed_successfully = isTRUE(completed_successfully),
        parameter_failed = isTRUE(failed),
        rollback_applied_after_round = isTRUE(failed),
        failure_reason_if_any = if (isTRUE(failed)) as.character(failure_reason) else "",
        stringsAsFactors = FALSE
      )

      if (isTRUE(completed_successfully)) {
        successful_round_record <- round_record
      }

      round_log_df <- rbind(round_log_df, round_record)
      append_round_logs(round_log_df, round_log_path)
      last_round_record <- round_record

      write_table_snapshot(
        parameter_table_state,
        file.path(
          snapshot_dir,
          sprintf("parameter_table_after_round_%02d_%s_round_%02d.csv", as.integer(param_index), param_name, as.integer(round_index))
        )
      )

      if (isTRUE(completed_successfully) || isTRUE(failed)) {
        break
      }
    }

    if (isTRUE(completed_successfully)) {
      working_table <- parameter_table_state
      working_table$estimate <- original_estimate
      final_snapshot_path <- file.path(
        snapshot_dir,
        sprintf("parameter_table_after_successful_parameter_%02d_%s.csv", as.integer(param_index), param_name)
      )
      write_table_snapshot(working_table, final_snapshot_path)
    } else {
      rollback_applied <- TRUE
      working_table <- pre_parameter_table
      working_table$estimate <- original_estimate
      final_snapshot_path <- file.path(
        snapshot_dir,
        sprintf("parameter_table_after_failed_parameter_rollback_%02d_%s.csv", as.integer(param_index), param_name)
      )
      write_table_snapshot(working_table, final_snapshot_path)
      rollback_record_path <- file.path(
        rollback_dir,
        sprintf("rollback_record_%02d_%s.txt", as.integer(param_index), param_name)
      )
      write_rollback_record(
        out_path = rollback_record_path,
        param_name = param_name,
        param_index = param_index,
        rounds_attempted = rounds_attempted,
        failure_reason = failure_reason,
        pre_parameter_snapshot_path = pre_parameter_snapshot_path,
        rollback_snapshot_path = final_snapshot_path,
        last_round_dir = if (is.null(last_round_record)) NA_character_ else as.character(last_round_record$run_dir[[1]])
      )
    }

    row_idx <- find_parameter_index(working_table, param_name)
    summary_record <- data.frame(
      parameter = param_name,
      parameter_index = as.integer(param_index),
      configured_max_rounds_per_parameter = as.integer(settings$max_rounds_per_parameter),
      configured_seeds_per_round = as.integer(settings$seeds_per_round),
      rounds_attempted = as.integer(rounds_attempted),
      status = if (isTRUE(completed_successfully)) "completed_successfully" else "failed_rolled_back",
      completed_successfully = isTRUE(completed_successfully),
      failed = !isTRUE(completed_successfully),
      rollback_applied = isTRUE(rollback_applied),
      failure_reason = if (isTRUE(completed_successfully)) "" else as.character(failure_reason),
      valid_reference_seed_found_in_any_round = isTRUE(valid_reference_seed_found_in_any_round),
      successful_reference_seed = if (is.null(successful_round_record)) NA_character_ else as.character(successful_round_record$eligible_reference_seed[[1]]),
      successful_reference_rank = if (is.null(successful_round_record)) NA_real_ else as.numeric(successful_round_record$eligible_reference_rank[[1]]),
      successful_selected_value = if (is.null(successful_round_record)) NA_real_ else as.numeric(successful_round_record$eligible_reference_value[[1]]),
      successful_selected_objective = if (is.null(successful_round_record)) NA_real_ else as.numeric(successful_round_record$eligible_reference_objective[[1]]),
      successful_selected_objective_burden = if (is.null(successful_round_record)) NA_real_ else as.numeric(successful_round_record$eligible_reference_objective_burden[[1]]),
      successful_selected_objective_ploidy = if (is.null(successful_round_record)) NA_real_ else as.numeric(successful_round_record$eligible_reference_objective_ploidy[[1]]),
      successful_min_pred_burden_2N_day1000 = if (is.null(successful_round_record)) NA_real_ else as.numeric(successful_round_record$eligible_reference_min_pred_burden_2N_day1000[[1]]),
      successful_min_pred_burden_4N_day1000 = if (is.null(successful_round_record)) NA_real_ else as.numeric(successful_round_record$eligible_reference_min_pred_burden_4N_day1000[[1]]),
      successful_boundary_rel = if (is.null(successful_round_record)) NA_real_ else as.numeric(successful_round_record$eligible_reference_boundary_rel[[1]]),
      successful_boundary_status = if (is.null(successful_round_record)) NA_character_ else as.character(successful_round_record$eligible_reference_boundary_status[[1]]),
      final_working_init_value = as.numeric(working_table$init_value[[row_idx]]),
      final_working_lower_bound = as.numeric(working_table$lower_bound[[row_idx]]),
      final_working_upper_bound = as.numeric(working_table$upper_bound[[row_idx]]),
      pre_parameter_snapshot_path = pre_parameter_snapshot_path,
      final_snapshot_path = final_snapshot_path,
      rollback_record_path = if (!is.na(rollback_record_path) && nzchar(rollback_record_path)) rollback_record_path else NA_character_,
      final_run_dir = if (is.null(last_round_record)) NA_character_ else as.character(last_round_record$run_dir[[1]]),
      stringsAsFactors = FALSE
    )
    parameter_summary_df <- rbind(parameter_summary_df, summary_record)
    write_parameter_summary(parameter_summary_df, parameter_summary_path)
  }

  working_table$estimate <- original_estimate
  write_natural_parameter_table(working_table, updated_parameter_table_path)
  write_natural_parameter_table(working_table, parameter_table_path)

  write_final_report(
    out_path = report_path,
    config_path = config_path,
    parameter_table_path = parameter_table_path,
    output_root = output_root,
    settings = settings,
    round_log_df = round_log_df,
    parameter_summary_df = parameter_summary_df,
    final_parameter_table = working_table
  )

  message("Wrote round log: ", round_log_path)
  message("Wrote parameter summary: ", parameter_summary_path)
  message("Wrote report: ", report_path)
  message("Wrote updated parameter table: ", updated_parameter_table_path)
  message("Overwrote original parameter table: ", parameter_table_path)
}

if (sys.nframe() == 0) {
  main()
}
