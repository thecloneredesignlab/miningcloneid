#!/usr/bin/env Rscript

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}

SCRIPT_DIR <- local_script_dir()
ANALYSIS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
WORKFLOW_DIR <- normalizePath(file.path(ANALYSIS_DIR, ".."), mustWork = FALSE)
REPO_ROOT <- normalizePath(file.path(WORKFLOW_DIR, "..", "..", ".."), mustWork = FALSE)
FIXO2_SCRIPT_PATH <- file.path(WORKFLOW_DIR, "simulation", "fix_o2_simulation.R")
MONOTONICITY_SCRIPT_PATH <- file.path(SCRIPT_DIR, "fixed_o2_ploidy_monotonicity.R")
INITIAL_PLOIDY_SCRIPT_PATH <- file.path(SCRIPT_DIR, "fixed_o2_initial_ploidy_trajectory.R")
CURVE_UTILS_PATH <- file.path(SCRIPT_DIR, "curve_classification_utils.R")

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[[1]], fixed = TRUE)
    val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x)) return(default)
  val <- tolower(as.character(x[[1]]))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer(x[[1]] %||% default))
  if (length(val) && is.finite(val)) val else default
}

as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric(x[[1]] %||% default))
  if (length(val) && is.finite(val)) val else default
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(as.character(x[[1]]))) return(default)
  vals <- suppressWarnings(as.numeric(trimws(strsplit(as.character(x[[1]]), ",", fixed = TRUE)[[1]])))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals else default
}

format_num <- function(x) {
  format(as.numeric(x), scientific = FALSE, trim = TRUE)
}

default_simulation_times <- function() {
  c(0, 1, 2, 5, 10, 20, 50, 100, 200, 300, 500, 700, 1000)
}

default_plot_times <- function() {
  c(100, 200, 300, 500, 700, 1000)
}

seed_number <- function(seed_id) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

normalize_seed_ids <- function(seed_ids) {
  seed_ids <- as.character(seed_ids)
  ifelse(grepl("^seed", seed_ids), seed_ids, paste0("seed", seed_ids))
}

num_key <- function(x) {
  vapply(x, function(xx) format(signif(as.numeric(xx), 12), scientific = FALSE, trim = TRUE), character(1))
}

seed_file_stem <- function(seed_id) {
  sn <- seed_number(seed_id)
  if (is.finite(sn)) sprintf("seed%03d", sn) else as.character(seed_id)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

read_tsv_auto <- function(path) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  } else {
    read_tsv(path)
  }
}

rbind_nonempty <- function(rows) {
  rows <- rows[vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

source_script_env <- function(script_path) {
  script_path <- normalizePath(path.expand(script_path), mustWork = TRUE)
  env <- new.env(parent = globalenv())
  env$commandArgs <- function(trailingOnly = FALSE) {
    character(0)
  }
  old_error <- getOption("error")
  on.exit(options(error = old_error), add = TRUE)
  source(script_path, local = env, chdir = TRUE)
  options(error = old_error)
  env
}

load_fixo2_env <- function(script_path = FIXO2_SCRIPT_PATH) {
  workflow_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
  shared_path <- file.path(workflow_root, "util", "o2_supply_demand_map_shared.R")
  common_path <- file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R")
  if (file.exists(shared_path)) source(shared_path, local = globalenv(), chdir = FALSE)
  if (file.exists(common_path)) source(common_path, local = globalenv(), chdir = FALSE)
  env <- source_script_env(script_path)
  required <- c(
    "fixo2_dominant_attractor_one", "fixo2_assign_attractor_modes", "fixo2_attractor_mode_table",
    "o2ipa_collect_seed_inputs", "o2ipa_params_wide", "o2pr_first_seed_cfg",
    "o2pr_run_params_from_vec", "o2ipa_source_model", "cf2_fixed_matrix", "cf2_init_vector",
    "fixo2_trajectory_with_fallback", "fixo2_dominant_from_eig"
  )
  missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = TRUE)]
  if (length(missing)) stop("FixO2 helper environment is missing: ", paste(missing, collapse = ", "))
  env
}

default_run_dir <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "fit_invivo_O2_buffering_500seed")
}

default_result_root <- function() {
  file.path(REPO_ROOT, "oxygen", "results", "analysis", "monotonicity_classification")
}

default_out_dir <- function(part) {
  if (identical(part, "initial_ploidy")) {
    file.path(default_result_root(), "dense-grid_initial-ploidy_trajectory")
  } else {
    file.path(default_result_root(), "dense-grid_monotonicity_classification")
  }
}

normalize_part <- function(part) {
  part <- tolower(gsub("-", "_", part %||% "monotonicity"))
  if (part %in% c("monotonicity", "dense_grid", "classification", "ploidy_monotonicity")) return("monotonicity")
  if (part %in% c("initial_ploidy", "initial", "trajectory", "initial_ploidy_trajectory")) return("initial_ploidy")
  stop("Unknown part: ", part)
}

task_list_paths <- function(out_dir, part) {
  task_dir <- file.path(out_dir, "hpc", "task_lists")
  list(
    task_dir = task_dir,
    task_file = file.path(task_dir, paste0(part, "_seed_o2_tasks.tsv")),
    metadata = file.path(task_dir, paste0(part, "_task_metadata.tsv"))
  )
}

chunk_paths <- function(out_dir, part) {
  if (identical(part, "initial_ploidy")) {
    list(
      root = file.path(out_dir, "tables", "task_chunks", part),
      daily = file.path(out_dir, "tables", "task_chunks", part, "daily"),
      manifest = file.path(out_dir, "tables", "task_chunks", part, "manifest")
    )
  } else {
    list(
      root = file.path(out_dir, "tables", "task_chunks", part),
      curves = file.path(out_dir, "tables", "task_chunks", part, "curves")
    )
  }
}

selected_seed_ids <- function(run_dir, max_seeds = NA_integer_, fixo2_env = load_fixo2_env()) {
  inputs <- get("o2ipa_collect_seed_inputs", envir = fixo2_env, inherits = TRUE)(run_dir, objective_source = "auto")
  param_mat <- get("o2ipa_params_wide", envir = fixo2_env, inherits = TRUE)(inputs$params_long, "value")
  seeds <- rownames(param_mat)
  seeds <- seeds[order(seed_number(seeds))]
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seeds <- seeds[seq_len(min(max_seeds, length(seeds)))]
  }
  if (!length(seeds)) stop("No fitted seed parameters found under: ", run_dir)
  seeds
}

build_tasks <- function(args) {
  part <- normalize_part(args$part %||% args$workflow_part)
  run_dir <- normalizePath(path.expand(args$run_dir %||% args$fit_root %||% default_run_dir()), mustWork = FALSE)
  out_dir <- normalizePath(path.expand(args$out_dir %||% args$output_root %||% default_out_dir(part)), mustWork = FALSE)
  if (!dir.exists(run_dir)) stop("run_dir does not exist: ", run_dir)

  o2_min <- as_num(args$o2_min, 0)
  o2_max <- as_num(args$o2_max, 5)
  o2_by <- as_num(args$o2_by, 0.025)
  o2_grid <- sort(unique(as_num_vec(args$o2_grid, seq(o2_min, o2_max, by = o2_by))))
  if (!length(o2_grid)) stop("O2 grid is empty.")
  max_seeds <- as_int(args$max_seeds, NA_integer_)
  tasks_per_array_task <- max(1L, as_int(args$tasks_per_array_task, 10L))
  force <- as_bool(args$force %||% args$overwrite, TRUE)

  fixo2_env <- load_fixo2_env()
  seeds <- selected_seed_ids(run_dir, max_seeds, fixo2_env)
  task_grid <- expand.grid(seed_id = seeds, O2_pct = o2_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  task_grid$seed_number <- seed_number(task_grid$seed_id)
  task_grid <- task_grid[order(task_grid$seed_number, task_grid$O2_pct), , drop = FALSE]
  task_grid$task_id <- seq_len(nrow(task_grid))
  task_grid$O2_key <- num_key(task_grid$O2_pct)
  task_grid$array_task_id <- ((task_grid$task_id - 1L) %/% tasks_per_array_task) + 1L
  task_grid <- task_grid[, c("task_id", "array_task_id", "seed_id", "seed_number", "O2_pct", "O2_key")]

  paths <- task_list_paths(out_dir, part)
  dir.create(paths$task_dir, recursive = TRUE, showWarnings = FALSE)
  chunks <- chunk_paths(out_dir, part)
  if (force && dir.exists(chunks$root)) unlink(chunks$root, recursive = TRUE, force = TRUE)
  if (force && identical(part, "initial_ploidy")) {
    daily_dir <- file.path(out_dir, "tables", "daily_trajectories")
    if (dir.exists(daily_dir)) unlink(daily_dir, recursive = TRUE, force = TRUE)
  }
  write_tsv(task_grid, paths$task_file)

  metadata <- data.frame(
    key = c(
      "part", "run_dir", "out_dir", "task_file", "n_seed", "n_o2", "n_tasks",
      "tasks_per_array_task", "n_array_tasks", "o2_grid", "max_seeds"
    ),
    value = c(
      part, run_dir, out_dir, paths$task_file, as.character(length(seeds)), as.character(length(o2_grid)),
      as.character(nrow(task_grid)), as.character(tasks_per_array_task),
      as.character(max(task_grid$array_task_id)), paste(format_num(o2_grid), collapse = ","),
      as.character(max_seeds)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(metadata, paths$metadata)
  message("Built ", part, " task list: ", nrow(task_grid), " seed x O2 tasks across ", max(task_grid$array_task_id), " array tasks")
  invisible(paths)
}

read_metadata_value <- function(metadata_path, key) {
  metadata <- read_tsv(metadata_path)
  val <- metadata$value[metadata$key == key]
  if (length(val)) val[[1L]] else NA_character_
}

task_slice <- function(args, part) {
  out_dir <- normalizePath(path.expand(args$out_dir %||% args$output_root %||% default_out_dir(part)), mustWork = FALSE)
  paths <- task_list_paths(out_dir, part)
  task_file <- normalizePath(path.expand(args$task_file %||% paths$task_file), mustWork = TRUE)
  tasks <- read_tsv(task_file)
  tasks_per_array_task <- as_int(args$tasks_per_array_task, NA_integer_)
  if (!is.finite(tasks_per_array_task) || is.na(tasks_per_array_task)) {
    meta_path <- paths$metadata
    if (file.exists(meta_path)) tasks_per_array_task <- as_int(read_metadata_value(meta_path, "tasks_per_array_task"), 10L)
  }
  tasks_per_array_task <- max(1L, tasks_per_array_task)
  array_task_id <- as_int(args$array_task_id, as_int(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_), NA_integer_))
  if (!is.finite(array_task_id) || is.na(array_task_id) || array_task_id < 1L) {
    stop("--array_task_id or SLURM_ARRAY_TASK_ID is required for run_tasks mode.")
  }
  start <- (array_task_id - 1L) * tasks_per_array_task + 1L
  end <- min(nrow(tasks), array_task_id * tasks_per_array_task)
  if (start > nrow(tasks)) {
    return(list(out_dir = out_dir, task_file = task_file, array_task_id = array_task_id, tasks = tasks[0, , drop = FALSE]))
  }
  list(out_dir = out_dir, task_file = task_file, array_task_id = array_task_id, tasks = tasks[seq.int(start, end), , drop = FALSE])
}

run_tasks_monotonicity <- function(args) {
  part <- "monotonicity"
  slice <- task_slice(args, part)
  tasks <- slice$tasks
  chunks <- chunk_paths(slice$out_dir, part)
  dir.create(chunks$curves, recursive = TRUE, showWarnings = FALSE)
  chunk_file <- file.path(chunks$curves, sprintf("array_%06d.tsv", slice$array_task_id))
  force <- as_bool(args$force %||% args$overwrite, TRUE)
  if (file.exists(chunk_file) && !force) stop("Chunk already exists; rerun with --force=TRUE: ", chunk_file)
  if (!nrow(tasks)) {
    write_tsv(data.frame(), chunk_file)
    return(invisible(chunk_file))
  }

  run_dir <- normalizePath(path.expand(args$run_dir %||% args$fit_root %||% default_run_dir()), mustWork = TRUE)
  fixo2_env <- load_fixo2_env()
  inputs <- get("o2ipa_collect_seed_inputs", envir = fixo2_env, inherits = TRUE)(run_dir, objective_source = "auto")
  param_mat <- get("o2ipa_params_wide", envir = fixo2_env, inherits = TRUE)(inputs$params_long, "value")
  cfg <- get("o2pr_first_seed_cfg", envir = fixo2_env, inherits = TRUE)(inputs$manifest)
  model_env <- get("o2ipa_source_model", envir = fixo2_env, inherits = TRUE)(dirname(FIXO2_SCRIPT_PATH))
  run_params_fn <- get("o2pr_run_params_from_vec", envir = fixo2_env, inherits = TRUE)
  dominant_fn <- get("fixo2_dominant_attractor_one", envir = fixo2_env, inherits = TRUE)
  mode_table_fn <- get("fixo2_attractor_mode_table", envir = fixo2_env, inherits = TRUE)

  rows <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    seed <- as.character(tasks$seed_id[[i]])
    if (!seed %in% rownames(param_mat)) stop("Seed not present in parameter table: ", seed)
    pvec <- as.numeric(param_mat[seed, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- run_params_fn(pvec, cfg)
    rows[[i]] <- dominant_fn(seed, run_params, model_env, cfg, as.numeric(tasks$O2_pct[[i]]))
  }
  attractors <- rbind_nonempty(rows)
  manifest <- inputs$manifest
  if ("objective" %in% names(manifest)) {
    manifest$delta_objective <- manifest$objective - min(manifest$objective, na.rm = TRUE)
  }
  attractors <- merge(
    attractors,
    manifest[, intersect(c("seed_id", "objective", "delta_objective"), names(manifest)), drop = FALSE],
    by = "seed_id",
    all.x = TRUE,
    sort = FALSE
  )
  curves <- mode_table_fn(attractors)
  curves$array_task_id <- slice$array_task_id
  curves$task_id <- tasks$task_id[match(paste(curves$seed_id, curves$O2_key), paste(tasks$seed_id, tasks$O2_key))]
  curves <- curves[order(curves$task_id), , drop = FALSE]
  write_tsv(curves, chunk_file)
  message("Wrote monotonicity chunk: ", chunk_file, " rows=", nrow(curves))
  invisible(chunk_file)
}

initial_task_rows <- function(tasks, param_mat, cfg, model_env, fixed_matrix_fn, init_vector_fn,
                              run_params_fn, trajectory_fn, dominant_fn, initial_specs_fn,
                              time_grid, selected_times) {
  rows <- list()
  manifest_rows <- list()
  row_i <- 0L
  man_i <- 0L
  n_unit <- as.numeric(cfg$N_UNIT %||% 22)
  init_specs <- initial_specs_fn(n_unit)
  for (i in seq_len(nrow(tasks))) {
    seed_id <- as.character(tasks$seed_id[[i]])
    seed_num <- seed_number(seed_id)
    O2 <- as.numeric(tasks$O2_pct[[i]])
    task_status <- character()
    task_rows <- list()
    task_i <- 0L
    init_used <- data.frame()
    pvec <- as.numeric(param_mat[seed_id, , drop = TRUE])
    names(pvec) <- colnames(param_mat)
    run_params <- run_params_fn(pvec, cfg)
    fm <- tryCatch(fixed_matrix_fn(model_env, cfg, run_params, O2), error = function(e) e)
    if (inherits(fm, "error")) {
      task_status <- c(task_status, paste0("matrix_error:", conditionMessage(fm)))
    } else {
      eig <- tryCatch(eigen(fm$M, only.values = FALSE), error = function(e) e)
      if (inherits(eig, "error")) {
        task_status <- c(task_status, paste0("eigen_error:", conditionMessage(eig)))
      } else {
        dom <- dominant_fn(eig, fm$ngrid, n_unit)
        for (j in seq_len(nrow(init_specs))) {
          init <- init_vector_fn(fm$ngrid, init_specs$requested_initial_N[[j]])
          init_used <- rbind(init_used, data.frame(
            initial_condition = init_specs$initial_condition[[j]],
            requested_initial_ploidy = init_specs$initial_ploidy[[j]],
            requested_initial_N = init_specs$requested_initial_N[[j]],
            used_initial_N = init$used_N,
            used_initial_ploidy = init$used_ploidy,
            stringsAsFactors = FALSE
          ))
          sim <- trajectory_fn(fm$M, eig, fm$ngrid, init$vector, time_grid, n_unit)
          task_status <- c(task_status, sim$status)
          tr <- sim$trajectory
          if (!nrow(tr)) next
          tr$seed_id <- seed_id
          tr$seed_number <- seed_num
          tr$O2_pct <- O2
          tr$O2_key <- num_key(O2)
          tr$initial_condition <- init_specs$initial_condition[[j]]
          tr$initial_ploidy <- init_specs$initial_ploidy[[j]]
          tr$requested_initial_N <- init_specs$requested_initial_N[[j]]
          tr$used_initial_N <- init$used_N
          tr$used_initial_ploidy <- init$used_ploidy
          tr$status <- sim$status
          tr$trajectory_method <- sim$trajectory_method
          tr$dominant_mean_N <- dom$dominant_mean_N[[1L]]
          tr$dominant_mean_ploidy <- dom$dominant_mean_ploidy[[1L]]
          tr$dominant_fraction_N_le_25 <- dom$dominant_fraction_N_le_25[[1L]]
          tr$dominant_fraction_N_below_44 <- dom$dominant_fraction_N_below_44[[1L]]
          tr$dominant_fraction_N_ge_44 <- dom$dominant_fraction_N_ge_44[[1L]]
          tr$dominant_growth_rate <- dom$dominant_growth_rate[[1L]]
          tr$second_growth_rate <- dom$second_growth_rate[[1L]]
          tr$spectral_gap <- dom$spectral_gap[[1L]]
          tr$relative_spectral_gap <- dom$relative_spectral_gap[[1L]]
          tr$relax_time_days <- dom$relax_time_days[[1L]]
          tr$time_to_10x_days <- dom$time_to_10x_days[[1L]]
          tr$time_to_100x_days <- dom$time_to_100x_days[[1L]]
          tr$log10_advantage_1000d <- dom$log10_advantage_1000d[[1L]]
          tr$dominance_class <- dom$dominance_class[[1L]]
          task_i <- task_i + 1L
          task_rows[[task_i]] <- tr[, c(
            "seed_id", "seed_number", "O2_pct", "O2_key", "initial_condition", "initial_ploidy",
            "requested_initial_N", "used_initial_N", "used_initial_ploidy", "status", "day",
            "trajectory_method", "mean_N", "mean_ploidy", "fraction_N_le_25", "fraction_N_below_44",
            "fraction_N_ge_44", "fraction_N_ge_66", "fraction_N_ge_88", "dominant_mean_N",
            "dominant_mean_ploidy", "dominant_fraction_N_le_25", "dominant_fraction_N_below_44",
            "dominant_fraction_N_ge_44", "dominant_growth_rate", "second_growth_rate", "spectral_gap",
            "relative_spectral_gap", "relax_time_days", "time_to_10x_days", "time_to_100x_days",
            "log10_advantage_1000d", "dominance_class"
          )]
        }
      }
    }
    task_daily <- rbind_nonempty(task_rows)
    expected_rows <- 2L * length(time_grid)
    status <- if (!nrow(task_daily)) {
      "failed_no_rows"
    } else if (nrow(task_daily) != expected_rows) {
      "partial"
    } else if (all(unique(task_status) == "ok")) {
      "ok"
    } else {
      "ok_with_warnings"
    }
    if (nrow(task_daily)) {
      row_i <- row_i + 1L
      rows[[row_i]] <- task_daily
    }
    man_i <- man_i + 1L
    manifest_rows[[man_i]] <- data.frame(
      task_id = tasks$task_id[[i]],
      array_task_id = tasks$array_task_id[[i]],
      seed_id = seed_id,
      seed_number = seed_num,
      O2_pct = O2,
      O2_key = num_key(O2),
      expected_rows = expected_rows,
      observed_rows = nrow(task_daily),
      selected_rows = if (nrow(task_daily)) sum(task_daily$day %in% selected_times) else 0L,
      status = status,
      status_values = paste(sort(unique(task_status)), collapse = ";"),
      n_expm_fallback_trajectories = if (nrow(task_daily)) {
        length(unique(paste(task_daily$O2_pct, task_daily$initial_condition, sep = "\r")[task_daily$trajectory_method == "expm_fallback"]))
      } else {
        0L
      },
      init_2N_used_initial_N = if (nrow(init_used)) init_used$used_initial_N[init_used$initial_condition == "init_2N"][[1L]] else NA_real_,
      init_4N_used_initial_N = if (nrow(init_used)) init_used$used_initial_N[init_used$initial_condition == "init_4N"][[1L]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  list(daily = rbind_nonempty(rows), manifest = rbind_nonempty(manifest_rows))
}

run_tasks_initial_ploidy <- function(args) {
  part <- "initial_ploidy"
  slice <- task_slice(args, part)
  tasks <- slice$tasks
  chunks <- chunk_paths(slice$out_dir, part)
  dir.create(chunks$daily, recursive = TRUE, showWarnings = FALSE)
  dir.create(chunks$manifest, recursive = TRUE, showWarnings = FALSE)
  daily_file <- file.path(chunks$daily, sprintf("array_%06d.tsv.gz", slice$array_task_id))
  manifest_file <- file.path(chunks$manifest, sprintf("array_%06d.tsv", slice$array_task_id))
  force <- as_bool(args$force %||% args$overwrite, TRUE)
  if ((file.exists(daily_file) || file.exists(manifest_file)) && !force) {
    stop("Initial-ploidy chunk already exists; rerun with --force=TRUE: ", daily_file)
  }
  if (!nrow(tasks)) {
    write_tsv_gz(data.frame(), daily_file)
    write_tsv(data.frame(), manifest_file)
    return(invisible(c(daily_file, manifest_file)))
  }

  run_dir <- normalizePath(path.expand(args$run_dir %||% args$fit_root %||% default_run_dir()), mustWork = TRUE)
  time_end <- as_int(args$time_end, 1000L)
  time_grid <- seq.int(0L, time_end)
  selected_times <- sort(unique(as_num_vec(args$simulation_times %||% args$selected_times, default_simulation_times())))
  selected_times <- selected_times[selected_times %in% time_grid]
  if (!length(selected_times)) stop("selected_times is empty after intersecting with time_grid.")

  fixo2_env <- load_fixo2_env()
  init_env <- source_script_env(INITIAL_PLOIDY_SCRIPT_PATH)
  assign("fixo2_helper", function(name) get(name, envir = fixo2_env, inherits = TRUE), envir = init_env)
  inputs <- get("o2ipa_collect_seed_inputs", envir = fixo2_env, inherits = TRUE)(run_dir, objective_source = "auto")
  param_mat <- get("o2ipa_params_wide", envir = fixo2_env, inherits = TRUE)(inputs$params_long, "value")
  cfg <- get("o2pr_first_seed_cfg", envir = fixo2_env, inherits = TRUE)(inputs$manifest)
  model_env <- get("o2ipa_source_model", envir = fixo2_env, inherits = TRUE)(dirname(FIXO2_SCRIPT_PATH))

  pvec_names <- colnames(param_mat)
  run_params_fn <- get("o2pr_run_params_from_vec", envir = fixo2_env, inherits = TRUE)
  wrapped_run_params <- function(pvec, cfg) {
    names(pvec) <- pvec_names
    run_params_fn(pvec, cfg)
  }
  out <- initial_task_rows(
    tasks = tasks,
    param_mat = param_mat,
    cfg = cfg,
    model_env = model_env,
    fixed_matrix_fn = get("cf2_fixed_matrix", envir = fixo2_env, inherits = TRUE),
    init_vector_fn = get("cf2_init_vector", envir = fixo2_env, inherits = TRUE),
    run_params_fn = wrapped_run_params,
    trajectory_fn = get("fixo2_trajectory_with_fallback", envir = fixo2_env, inherits = TRUE),
    dominant_fn = get("fixo2_dominant_from_eig", envir = fixo2_env, inherits = TRUE),
    initial_specs_fn = get("initial_specs", envir = init_env, inherits = TRUE),
    time_grid = time_grid,
    selected_times = selected_times
  )
  if (nrow(out$daily)) {
    out$daily <- out$daily[order(out$daily$seed_number, out$daily$O2_pct, out$daily$initial_ploidy, out$daily$day), , drop = FALSE]
  }
  write_tsv_gz(out$daily, daily_file)
  out$manifest$chunk_daily_file <- daily_file
  write_tsv(out$manifest, manifest_file)
  message("Wrote initial-ploidy chunk: ", daily_file, " rows=", nrow(out$daily))
  invisible(c(daily_file, manifest_file))
}

run_tasks <- function(args) {
  part <- normalize_part(args$part %||% args$workflow_part)
  if (identical(part, "initial_ploidy")) run_tasks_initial_ploidy(args) else run_tasks_monotonicity(args)
}

monotonicity_paths <- function(out_dir) {
  list(
    curves = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_curves.tsv"),
    by_seed = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv"),
    crosswalk = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_mode_crosswalk.tsv"),
    class_counts = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_counts.tsv"),
    class_by_mode = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_by_reporting_o2_mode.tsv"),
    class_curves = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_curve_summary.tsv"),
    representatives = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_representative_seeds.tsv"),
    objective_rank = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_by_objective_rank.tsv"),
    parameter_differences = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_parameter_differences.tsv"),
    class_change_audit = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_class_change_audit.tsv"),
    run_arguments = file.path(out_dir, "tables", "fixed_o2_ploidy_monotonicity_run_arguments.tsv")
  )
}

merge_monotonicity <- function(args) {
  part <- "monotonicity"
  run_dir <- normalizePath(path.expand(args$run_dir %||% default_run_dir()), mustWork = TRUE)
  out_dir <- normalizePath(path.expand(args$out_dir %||% default_out_dir(part)), mustWork = FALSE)
  o2_min <- as_num(args$o2_min, 0)
  o2_max <- as_num(args$o2_max, 5)
  o2_by <- as_num(args$o2_by, 0.025)
  o2_grid <- sort(unique(as_num_vec(args$o2_grid, seq(o2_min, o2_max, by = o2_by))))
  reporting_o2 <- sort(unique(as_num_vec(args$reporting_o2, c(0, 0.1, 0.5, 1, 2, 5))))
  max_seeds <- as_int(args$max_seeds, NA_integer_)
  gap_low <- as_num(args$gap_low, 0.005)
  gap_caution <- as_num(args$gap_caution, 0.01)
  unreliable_fraction <- as_num(args$unreliable_fraction, 0.25)
  caution_fraction <- as_num(args$caution_fraction, 0.10)
  flat_range_threshold <- as_num(args$flat_range_threshold, 0.05)
  step_epsilon_abs <- as_num(args$step_epsilon_abs, 1e-6)
  step_epsilon_fraction <- as_num(args$step_epsilon_fraction, 1e-4)
  reverse_fraction_tolerance <- as_num(args$reverse_fraction_tolerance, 0.05)
  plateau_min_points <- as_int(args$plateau_min_points, 3L)
  generate_figures <- as_bool(args$generate_figures, TRUE)
  validate <- as_bool(args$run_validation, TRUE)

  dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "report"), recursive = TRUE, showWarnings = FALSE)
  paths <- monotonicity_paths(out_dir)
  chunks <- chunk_paths(out_dir, part)
  chunk_files <- list.files(chunks$curves, pattern = "^array_[0-9]+\\.tsv$", full.names = TRUE)
  if (!length(chunk_files)) stop("No monotonicity chunk files found under: ", chunks$curves)
  chunk_files <- chunk_files[order(as.integer(gsub("\\D", "", basename(chunk_files))))]
  curves_raw <- rbind_nonempty(lapply(chunk_files, read_tsv))
  if (!nrow(curves_raw)) stop("No monotonicity rows were generated.")
  curves_raw <- curves_raw[order(seed_number(curves_raw$seed_id), curves_raw$O2_pct), , drop = FALSE]
  duplicate_key <- duplicated(paste(curves_raw$seed_id, curves_raw$O2_key))
  if (any(duplicate_key)) stop("Duplicate monotonicity seed/O2 rows detected: ", sum(duplicate_key))

  fixo2_env <- load_fixo2_env()
  mono_env <- source_script_env(MONOTONICITY_SCRIPT_PATH)
  if (validate) {
    get("run_validation", envir = mono_env, inherits = TRUE)(
      out_dir = out_dir,
      flat_range_threshold = flat_range_threshold,
      step_epsilon_abs = step_epsilon_abs,
      step_epsilon_fraction = step_epsilon_fraction,
      reverse_fraction_tolerance = reverse_fraction_tolerance,
      plateau_min_points = plateau_min_points
    )
  }
  inputs <- get("o2ipa_collect_seed_inputs", envir = fixo2_env, inherits = TRUE)(run_dir, objective_source = "auto")
  manifest <- inputs$manifest
  seeds <- selected_seed_ids(run_dir, max_seeds, fixo2_env)
  expected_rows <- length(seeds) * length(o2_grid)
  if (nrow(curves_raw) != expected_rows) {
    stop("Unexpected monotonicity row count: observed ", nrow(curves_raw), ", expected ", expected_rows)
  }
  bad_status <- curves_raw$status[is.na(curves_raw$status) | curves_raw$status != "ok"]
  if (length(bad_status)) warning("Non-ok attractor rows: ", length(bad_status))

  previous_by_seed <- if (file.exists(paths$by_seed)) {
    tryCatch(read_tsv(paths$by_seed), error = function(e) data.frame())
  } else {
    data.frame()
  }
  by_seed <- do.call(rbind, lapply(split(curves_raw, curves_raw$seed_id), get("classify_one_seed", envir = mono_env, inherits = TRUE),
                                   reporting_o2 = reporting_o2, gap_low = gap_low,
                                   gap_caution = gap_caution,
                                   unreliable_fraction = unreliable_fraction,
                                   caution_fraction = caution_fraction,
                                   flat_range_threshold = flat_range_threshold,
                                   step_epsilon_abs = step_epsilon_abs,
                                   step_epsilon_fraction = step_epsilon_fraction,
                                   reverse_fraction_tolerance = reverse_fraction_tolerance,
                                   plateau_min_points = plateau_min_points))
  by_seed <- merge(by_seed, manifest[, intersect(c("seed_id", "objective", "objective_source", "objective_total", "objective_data", "objective_burden", "objective_ploidy", "runtime", "convergence_status"), names(manifest)), drop = FALSE],
                   by = "seed_id", all.x = TRUE, sort = FALSE)
  by_seed <- by_seed[order(by_seed$seed_number), , drop = FALSE]

  curves <- get("add_curve_differences", envir = mono_env, inherits = TRUE)(curves_raw, by_seed, reporting_o2)
  class_counts <- get("class_count_table", envir = mono_env, inherits = TRUE)(by_seed)
  class_by_mode <- get("class_by_mode_table", envir = mono_env, inherits = TRUE)(by_seed, reporting_o2)
  class_curves <- get("curve_quantile_table", envir = mono_env, inherits = TRUE)(curves, by_seed)
  reps <- get("representative_seed_table", envir = mono_env, inherits = TRUE)(curves, by_seed)
  objective_rank <- get("objective_rank_class_table", envir = mono_env, inherits = TRUE)(by_seed)
  param_diff <- get("parameter_difference_table", envir = mono_env, inherits = TRUE)(by_seed, inputs$params_long)
  class_change_audit <- get("class_change_audit_table", envir = mono_env, inherits = TRUE)(previous_by_seed, by_seed)
  crosswalk_cols <- c(
    "seed_id", "seed_number", "curve_class", "final_interpretation_class",
    "monotonicity_reliability", "n_sign_changes", "ploidy_range", "net_ploidy_change",
    "step_epsilon", "flat_range_threshold", "fraction_positive_steps", "fraction_negative_steps",
    "low_amplitude_curve", "terminal_plateau_for_class", "min_spectral_gap",
    "fraction_o2_gap_below_0p005", "fraction_o2_gap_below_0p01",
    grep("^(mode_label|dominant_mean_ploidy|spectral_gap)_o2_", names(by_seed), value = TRUE)
  )
  crosswalk <- by_seed[, intersect(crosswalk_cols, names(by_seed)), drop = FALSE]
  run_args <- data.frame(
    argument = c(
      "run_dir", "out_dir", "script", "array_backend", "task_file", "o2_grid", "n_o2", "reporting_o2",
      "n_seed", "expected_curve_rows", "max_seeds", "classification_rule_version",
      "flat_range_threshold", "step_epsilon_rule", "step_epsilon_abs", "step_epsilon_fraction",
      "reverse_fraction_tolerance", "plateau_min_points", "gap_low", "gap_caution",
      "unreliable_fraction", "caution_fraction", "analytical_method"
    ),
    value = c(
      run_dir, out_dir, normalizePath(MONOTONICITY_SCRIPT_PATH, mustWork = FALSE),
      normalizePath(file.path(SCRIPT_DIR, "dense_grid_monotonicity_array_backend.R"), mustWork = FALSE),
      task_list_paths(out_dir, part)$task_file,
      paste(format_num(o2_grid), collapse = ","), as.character(length(o2_grid)),
      paste(format_num(reporting_o2), collapse = ","), as.character(length(seeds)),
      as.character(expected_rows), as.character(max_seeds),
      get("curve_classification_rule_version", envir = mono_env, inherits = TRUE)(),
      as.character(flat_range_threshold), "max(step_epsilon_abs, step_epsilon_fraction * ploidy_range)",
      as.character(step_epsilon_abs), as.character(step_epsilon_fraction),
      as.character(reverse_fraction_tolerance), as.character(plateau_min_points),
      as.character(gap_low), as.character(gap_caution),
      as.character(unreliable_fraction), as.character(caution_fraction),
      "fixed-O2 generator dominant eigenvector; seed x O2 Slurm-array analytical grid evaluation, not stochastic simulation"
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(curves, paths$curves)
  write_tsv(by_seed, paths$by_seed)
  write_tsv(crosswalk, paths$crosswalk)
  write_tsv(class_counts, paths$class_counts)
  write_tsv(class_by_mode, paths$class_by_mode)
  write_tsv(class_curves, paths$class_curves)
  write_tsv(reps, paths$representatives)
  write_tsv(objective_rank, paths$objective_rank)
  write_tsv(param_diff, paths$parameter_differences)
  write_tsv(class_change_audit, paths$class_change_audit)
  write_tsv(run_args, paths$run_arguments)

  if (generate_figures) {
    get("plot_all_curves", envir = mono_env, inherits = TRUE)(curves, by_seed, out_dir)
    get("plot_class_summary", envir = mono_env, inherits = TRUE)(class_curves, out_dir)
    get("plot_heatmap", envir = mono_env, inherits = TRUE)(curves, by_seed, out_dir, "dominant_mean_ploidy",
                                                           "fixed_o2_dominant_ploidy_heatmap_ordered_by_class", zlab = "Dominant mean ploidy")
    get("plot_mode_heatmap", envir = mono_env, inherits = TRUE)(curves, by_seed, out_dir)
    get("plot_heatmap", envir = mono_env, inherits = TRUE)(curves, by_seed, out_dir, "spectral_gap",
                                                           "fixed_o2_spectral_gap_heatmap_ordered_by_class",
                                                           transform = function(x) log10(pmax(x, .Machine$double.eps)), zlab = "log10 spectral gap")
    get("plot_gap_scatter", envir = mono_env, inherits = TRUE)(by_seed, out_dir)
    get("plot_representatives", envir = mono_env, inherits = TRUE)(curves, reps, out_dir)
    get("plot_class_seed_overlays", envir = mono_env, inherits = TRUE)(curves, by_seed, out_dir, stable_only = FALSE)
    get("plot_class_seed_overlays", envir = mono_env, inherits = TRUE)(curves, by_seed, out_dir, stable_only = TRUE)
  }
  get("write_summary_report", envir = mono_env, inherits = TRUE)(out_dir, args, run_args, class_counts, table(by_seed$monotonicity_reliability), by_seed, paths)
  message("Merged monotonicity array outputs: ", out_dir)
  invisible(paths)
}

initial_paths <- function(init_env, out_dir) {
  get("analysis_paths", envir = init_env, inherits = TRUE)(out_dir)
}

flush_seed_daily <- function(seed_rows, seed_id, selected_times, daily_dir, expected_rows) {
  daily <- rbind_nonempty(seed_rows)
  seed_num <- seed_number(seed_id)
  daily_path <- file.path(daily_dir, paste0(seed_file_stem(seed_id), ".tsv.gz"))
  if (nrow(daily)) {
    daily <- daily[order(daily$seed_number, daily$O2_pct, daily$initial_ploidy, daily$day), , drop = FALSE]
    write_tsv_gz(daily, daily_path)
  }
  status_values <- if (nrow(daily)) sort(unique(daily$status)) else character()
  status <- if (!nrow(daily)) {
    "failed_no_rows"
  } else if (nrow(daily) != expected_rows) {
    "partial"
  } else if (all(status_values == "ok")) {
    "ok"
  } else {
    "ok_with_warnings"
  }
  data.frame(
    seed_id = seed_id,
    seed_number = seed_num,
    daily_file = daily_path,
    expected_rows = expected_rows,
    observed_rows = nrow(daily),
    selected_rows = if (nrow(daily)) sum(daily$day %in% selected_times) else 0L,
    file_exists = file.exists(daily_path),
    status = status,
    status_values = paste(status_values, collapse = ";"),
    n_expm_fallback_trajectories = if (nrow(daily) && "trajectory_method" %in% names(daily)) {
      length(unique(paste(daily$O2_pct, daily$initial_condition, sep = "\r")[daily$trajectory_method == "expm_fallback"]))
    } else {
      0L
    },
    init_2N_used_initial_N = if (nrow(daily)) unique(daily$used_initial_N[daily$initial_condition == "init_2N"])[[1L]] else NA_real_,
    init_4N_used_initial_N = if (nrow(daily)) unique(daily$used_initial_N[daily$initial_condition == "init_4N"])[[1L]] else NA_real_,
    stringsAsFactors = FALSE
  )
}

merge_initial_daily_chunks <- function(out_dir, chunk_files, seeds, o2_grid, time_grid, selected_times) {
  daily_dir <- file.path(out_dir, "tables", "daily_trajectories")
  dir.create(daily_dir, recursive = TRUE, showWarnings = FALSE)
  expected_rows <- length(o2_grid) * 2L * length(time_grid)
  current_seed <- NULL
  current_rows <- list()
  selected_rows <- list()
  manifest_rows <- list()
  selected_i <- 0L
  manifest_i <- 0L

  flush_current <- function() {
    if (is.null(current_seed)) return(NULL)
    manifest_i <<- manifest_i + 1L
    manifest_rows[[manifest_i]] <<- flush_seed_daily(current_rows, current_seed, selected_times, daily_dir, expected_rows)
    current_rows <<- list()
    current_seed <<- NULL
    NULL
  }

  for (file in chunk_files) {
    daily <- read_tsv_auto(file)
    if (!nrow(daily)) next
    daily <- daily[order(daily$seed_number, daily$O2_pct, daily$initial_ploidy, daily$day), , drop = FALSE]
    sel <- daily[daily$day %in% selected_times, , drop = FALSE]
    if (nrow(sel)) {
      selected_i <- selected_i + 1L
      selected_rows[[selected_i]] <- sel
    }
    seed_order <- unique(daily$seed_id)
    for (seed in seed_order) {
      z <- daily[daily$seed_id == seed, , drop = FALSE]
      if (is.null(current_seed)) {
        current_seed <- seed
      } else if (!identical(current_seed, seed)) {
        flush_current()
        current_seed <- seed
      }
      current_rows[[length(current_rows) + 1L]] <- z
    }
  }
  flush_current()

  daily_manifest <- rbind_nonempty(manifest_rows)
  missing_seeds <- setdiff(seeds, daily_manifest$seed_id)
  if (length(missing_seeds)) {
    missing_rows <- lapply(missing_seeds, function(seed) data.frame(
      seed_id = seed,
      seed_number = seed_number(seed),
      daily_file = file.path(daily_dir, paste0(seed_file_stem(seed), ".tsv.gz")),
      expected_rows = expected_rows,
      observed_rows = 0L,
      selected_rows = 0L,
      file_exists = FALSE,
      status = "missing_seed_rows",
      status_values = "missing_seed_rows",
      n_expm_fallback_trajectories = 0L,
      init_2N_used_initial_N = NA_real_,
      init_4N_used_initial_N = NA_real_,
      stringsAsFactors = FALSE
    ))
    daily_manifest <- rbind_nonempty(c(list(daily_manifest), missing_rows))
  }
  daily_manifest <- daily_manifest[order(daily_manifest$seed_number), , drop = FALSE]
  selected <- rbind_nonempty(selected_rows)
  if (nrow(selected)) {
    selected <- selected[order(selected$seed_number, selected$O2_pct, selected$initial_ploidy, selected$day), , drop = FALSE]
  }
  list(daily_manifest = daily_manifest, selected = selected)
}

merge_initial_ploidy <- function(args) {
  part <- "initial_ploidy"
  fit_root <- normalizePath(path.expand(args$fit_root %||% args$run_dir %||% default_run_dir()), mustWork = TRUE)
  output_root <- normalizePath(path.expand(args$output_root %||% args$out_dir %||% default_out_dir(part)), mustWork = FALSE)
  o2_min <- as_num(args$o2_min, 0)
  o2_max <- as_num(args$o2_max, 5)
  o2_by <- as_num(args$o2_by, 0.025)
  o2_grid <- sort(unique(as_num_vec(args$o2_grid, seq(o2_min, o2_max, by = o2_by))))
  time_end <- as_int(args$time_end, 1000L)
  time_grid <- seq.int(0L, time_end)
  selected_times <- sort(unique(as_num_vec(args$simulation_times %||% args$selected_times, default_simulation_times())))
  selected_times <- selected_times[selected_times %in% time_grid]
  plot_times <- sort(unique(as_num_vec(args$plot_times, default_plot_times())))
  plot_times <- plot_times[plot_times %in% selected_times]
  max_seeds <- as_int(args$max_seeds, NA_integer_)
  force <- as_bool(args$force %||% args$overwrite, TRUE)
  flat_range_threshold <- as_num(args$flat_range_threshold, 0.05)
  step_epsilon_abs <- as_num(args$step_epsilon_abs, 1e-6)
  step_epsilon_fraction <- as_num(args$step_epsilon_fraction, 1e-4)
  reverse_fraction_tolerance <- as_num(args$reverse_fraction_tolerance, 0.05)
  plateau_min_points <- as_int(args$plateau_min_points, 3L)
  run_validation <- as_bool(args$run_validation, TRUE)
  generate_figures <- as_bool(args$generate_figures, TRUE)

  if (!length(o2_grid)) stop("O2 grid is empty.")
  if (!length(selected_times)) stop("selected_times is empty after intersecting with time_grid.")
  if (!length(plot_times)) stop("plot_times is empty after intersecting with selected_times.")
  table_dir <- file.path(output_root, "tables")
  daily_dir <- file.path(table_dir, "daily_trajectories")
  if (force && dir.exists(daily_dir)) unlink(daily_dir, recursive = TRUE, force = TRUE)
  dir.create(daily_dir, recursive = TRUE, showWarnings = FALSE)

  init_env <- source_script_env(INITIAL_PLOIDY_SCRIPT_PATH)
  fixo2_env <- load_fixo2_env()
  assign("fixo2_helper", function(name) get(name, envir = fixo2_env, inherits = TRUE), envir = init_env)
  paths <- initial_paths(init_env, output_root)
  chunks <- chunk_paths(output_root, part)
  chunk_files <- list.files(chunks$daily, pattern = "^array_[0-9]+\\.tsv\\.gz$", full.names = TRUE)
  if (!length(chunk_files)) stop("No initial-ploidy daily chunk files found under: ", chunks$daily)
  chunk_files <- chunk_files[order(as.integer(gsub("\\D", "", basename(chunk_files))))]

  seeds <- selected_seed_ids(fit_root, max_seeds, fixo2_env)
  merged <- merge_initial_daily_chunks(output_root, chunk_files, seeds, o2_grid, time_grid, selected_times)
  daily_manifest <- merged$daily_manifest
  selected_raw <- merged$selected
  if (!nrow(selected_raw)) stop("No selected trajectory rows were generated.")

  class_info <- get("classify_selected_curves", envir = init_env, inherits = TRUE)(
    selected_raw,
    flat_range_threshold = flat_range_threshold,
    step_epsilon_abs = step_epsilon_abs,
    step_epsilon_fraction = step_epsilon_fraction,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    plateau_min_points = plateau_min_points
  )
  selected <- get("attach_curve_differences", envir = init_env, inherits = TRUE)(selected_raw, class_info$differences)
  delta <- get("build_delta_table", envir = init_env, inherits = TRUE)(selected_raw)
  convergence_summary <- get("summarize_convergence", envir = init_env, inherits = TRUE)(delta, time_end = time_end)

  run_args <- data.frame(
    argument = c(
      "script", "array_backend", "fit_root", "output_root", "fixo2_script", "curve_classification_utils",
      "n_seed", "max_seeds", "o2_grid", "n_o2", "time_grid", "n_time",
      "selected_times", "n_selected_time", "plot_times", "n_plot_time", "initial_ploidy_values",
      "expected_daily_rows_per_seed", "expected_selected_rows", "expected_delta_rows",
      "expected_class_rows", "force", "classification_rule_version",
      "flat_range_threshold", "step_epsilon_abs", "step_epsilon_fraction",
      "reverse_fraction_tolerance", "plateau_min_points", "generate_figures", "analytical_method",
      "daily_file_pattern", "task_file"
    ),
    value = c(
      normalizePath(INITIAL_PLOIDY_SCRIPT_PATH, mustWork = FALSE),
      normalizePath(file.path(SCRIPT_DIR, "dense_grid_monotonicity_array_backend.R"), mustWork = FALSE),
      fit_root,
      output_root,
      normalizePath(FIXO2_SCRIPT_PATH, mustWork = FALSE),
      normalizePath(CURVE_UTILS_PATH, mustWork = FALSE),
      as.character(length(seeds)),
      as.character(max_seeds),
      paste(format_num(o2_grid), collapse = ","),
      as.character(length(o2_grid)),
      paste(range(time_grid), collapse = ":"),
      as.character(length(time_grid)),
      paste(format_num(selected_times), collapse = ","),
      as.character(length(selected_times)),
      paste(format_num(plot_times), collapse = ","),
      as.character(length(plot_times)),
      "2N,4N",
      as.character(length(o2_grid) * 2L * length(time_grid)),
      as.character(length(seeds) * length(o2_grid) * 2L * length(selected_times)),
      as.character(length(seeds) * length(o2_grid) * length(selected_times)),
      as.character(length(seeds) * length(selected_times) * 2L),
      as.character(force),
      get("curve_classification_rule_version", envir = init_env, inherits = TRUE)(),
      as.character(flat_range_threshold),
      as.character(step_epsilon_abs),
      as.character(step_epsilon_fraction),
      as.character(reverse_fraction_tolerance),
      as.character(plateau_min_points),
      as.character(generate_figures),
      "cf2_fixed_matrix + cf2_init_vector + cached eigen propagation equivalent to cf2_eigen_trajectory; seed x O2 Slurm-array finite-time analytical trajectory; no stochastic simulation; no refitting",
      file.path(daily_dir, "seed%03d.tsv.gz"),
      task_list_paths(output_root, part)$task_file
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(daily_manifest, paths$daily_manifest)
  write_tsv(selected, paths$selected)
  write_tsv(delta, paths$delta)
  write_tsv(class_info$class_table, paths$class_by_seed_time)
  write_tsv(convergence_summary, paths$convergence)
  write_tsv(run_args, paths$run_arguments)

  if (run_validation) {
    inputs <- get("o2ipa_collect_seed_inputs", envir = fixo2_env, inherits = TRUE)(fit_root, objective_source = "auto")
    param_mat <- get("o2ipa_params_wide", envir = fixo2_env, inherits = TRUE)(inputs$params_long, "value")
    cfg <- get("o2pr_first_seed_cfg", envir = fixo2_env, inherits = TRUE)(inputs$manifest)
    model_env <- get("o2ipa_source_model", envir = fixo2_env, inherits = TRUE)(dirname(FIXO2_SCRIPT_PATH))
    init_audit <- unique(selected_raw[, intersect(c(
      "seed_id", "seed_number", "initial_condition", "requested_initial_ploidy",
      "requested_initial_N", "used_initial_N", "used_initial_ploidy"
    ), names(selected_raw)), drop = FALSE])
    validation <- get("build_validation", envir = init_env, inherits = TRUE)(
      daily_manifest = daily_manifest,
      selected = selected,
      delta = delta,
      class_table = class_info$class_table,
      convergence_summary = convergence_summary,
      seeds = seeds,
      o2_grid = o2_grid,
      time_grid = time_grid,
      selected_times = selected_times,
      init_audit = init_audit,
      fixo2_env = fixo2_env,
      model_env = model_env,
      cfg = cfg,
      param_mat = param_mat,
      curve_utils_path = normalizePath(CURVE_UTILS_PATH, mustWork = FALSE)
    )
    write_tsv(validation, paths$validation)
    if (!all(validation$passed)) {
      stop("Validation failed: ", paste(validation$test_case[!validation$passed], collapse = ", "))
    }
  }
  if (generate_figures) {
    get("write_figures", envir = init_env, inherits = TRUE)(selected, delta, class_info$class_table, convergence_summary, daily_manifest, output_root, plot_times = plot_times)
  }
  message("Merged initial-ploidy array outputs: ", output_root)
  invisible(paths)
}

merge_outputs <- function(args) {
  part <- normalize_part(args$part %||% args$workflow_part)
  if (identical(part, "initial_ploidy")) merge_initial_ploidy(args) else merge_monotonicity(args)
}

main <- function(args = parse_args()) {
  mode <- tolower(gsub("-", "_", args$mode %||% ""))
  if (identical(mode, "build_tasks")) return(build_tasks(args))
  if (identical(mode, "run_tasks")) return(run_tasks(args))
  if (mode %in% c("merge", "merge_outputs", "finalize")) return(merge_outputs(args))
  stop("Unknown --mode. Use build_tasks, run_tasks, or merge.")
}

if (identical(environment(), globalenv())) {
  main(parse_args())
}
