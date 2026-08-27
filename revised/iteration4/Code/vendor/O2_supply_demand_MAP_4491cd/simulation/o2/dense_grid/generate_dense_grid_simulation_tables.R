#!/usr/bin/env Rscript

.dense_sim_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "generate_dense_grid_simulation_tables.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.dense_workflow_root <- normalizePath(file.path(.dense_sim_dir, "..", "..", ".."), mustWork = TRUE)
source(file.path(.dense_workflow_root, "util", "o2_supply_demand_map_dense_grid_utils.R"), local = environment(), chdir = TRUE)

dense_fixed_o2_environment <- local({
  cache <- NULL
  workflow_root <- .dense_workflow_root
  function(include_trajectory = FALSE) {
    if (!is.null(cache) && (!include_trajectory || exists("generate_analytical_trajectories", envir = cache, inherits = FALSE))) return(cache)
    if (is.null(cache)) {
      cache <<- new.env(parent = globalenv())
      cache$commandArgs <- function(trailingOnly = FALSE) character()
      sys.source(file.path(workflow_root, "simulation", "o2", "fixed_o2", "run_fixed_o2_simulation.R"), envir = cache, chdir = TRUE)
    }
    if (include_trajectory && !exists("generate_analytical_trajectories", envir = cache, inherits = FALSE)) {
      # Preserve the self-contained attractor producer from the canonical
      # fixed-O2 simulation entrypoint.  The trajectory module carries a
      # same-named legacy helper with additional runner-only dependencies.
      cache$dense_attractor_mode_table <- cache$generate_fixo2_attractor_mode_table
      sys.source(file.path(workflow_root, "util", "o2_supply_demand_map_fixed_o2_format_utils.R"), envir = cache, chdir = TRUE)
      # The numerical producer intentionally expects the legacy scalar-key
      # helper.  Keep that formatting dependency local to this simulation
      # environment instead of importing the legacy mixed-layer runner.
      cache$num_key <- dense_num_key
      sys.source(file.path(workflow_root, "simulation", "o2", "fixed_o2", "fixed_o2_numerical_producers.R"), envir = cache, chdir = TRUE)
    }
    cache
  }
})

dense_selected_seeds <- function(run_dir, max_seeds = NA_integer_, seed_manifest = NULL) {
  seeds <- character()
  if (!is.null(seed_manifest) && nzchar(seed_manifest) && file.exists(seed_manifest)) {
    manifest <- dense_read_tsv(seed_manifest)
    candidate <- intersect(c("seed_id", "seed", "seed_label"), names(manifest))
    if (length(candidate)) seeds <- dense_normalize_seed_ids(manifest[[candidate[[1L]]]])
  }
  if (!length(seeds)) {
    dirs <- list.dirs(run_dir, recursive = FALSE, full.names = FALSE)
    seeds <- dirs[grepl("^seed", dirs) & file.exists(file.path(run_dir, dirs, "best_params.tsv"))]
  }
  seeds <- unique(seeds[order(dense_seed_number(seeds), seeds)])
  if (is.finite(max_seeds) && max_seeds > 0L) seeds <- head(seeds, max_seeds)
  if (!length(seeds)) stop("No fitted seed directories found under ", run_dir, call. = FALSE)
  seeds
}

dense_task_paths <- function(out_dir, part) {
  task_dir <- file.path(out_dir, "hpc", "task_lists")
  list(
    task_file = file.path(task_dir, paste0(part, "_seed_o2_tasks.tsv")),
    metadata = file.path(task_dir, paste0(part, "_task_metadata.tsv")),
    chunk_dir = file.path(dense_grid_simulation_dir(out_dir), "task_chunks", part)
  )
}

dense_build_tasks <- function(argv, part, run_dir, out_dir, o2_grid) {
  seeds <- dense_selected_seeds(run_dir, dense_as_int(argv$max_seeds, NA_integer_), argv$seed_manifest %||% argv$seed_manifest_file)
  tasks_per <- max(1L, dense_as_int(argv$tasks_per_array_task, 10L))
  tasks <- expand.grid(seed_id = seeds, O2_pct = o2_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  tasks$seed_number <- dense_seed_number(tasks$seed_id)
  tasks <- tasks[order(tasks$seed_number, tasks$O2_pct), , drop = FALSE]
  tasks$task_id <- seq_len(nrow(tasks))
  tasks$array_task_id <- ((tasks$task_id - 1L) %/% tasks_per) + 1L
  tasks$O2_key <- dense_num_key(tasks$O2_pct)
  tasks <- tasks[, c("task_id", "array_task_id", "seed_id", "seed_number", "O2_pct", "O2_key")]
  paths <- dense_task_paths(out_dir, part)
  if (dense_as_bool(argv$overwrite, TRUE) && dir.exists(paths$chunk_dir)) unlink(paths$chunk_dir, recursive = TRUE, force = TRUE)
  dense_write_tsv(tasks, paths$task_file)
  metadata <- data.frame(
    key = c("part", "run_dir", "out_dir", "task_file", "n_seed", "n_o2", "n_tasks", "tasks_per_array_task", "n_array_tasks", "o2_grid"),
    value = c(part, run_dir, out_dir, paths$task_file, length(seeds), length(o2_grid), nrow(tasks), tasks_per, max(tasks$array_task_id), paste(o2_grid, collapse = ",")),
    stringsAsFactors = FALSE
  )
  dense_write_tsv(metadata, paths$metadata)
  message("Built ", part, " task list: ", nrow(tasks), " tasks across ", max(tasks$array_task_id), " array tasks")
  invisible(paths)
}

dense_generate_rows <- function(part, run_dir, seeds, o2_values, time_points, initial_ploidies, n_workers = 1L) {
  if (part == "monotonicity") {
    env <- dense_fixed_o2_environment(FALSE)
    return(get("generate_fixo2_attractor_mode_table", env)(run_dir, o2_values, seeds, n_workers))
  }
  env <- dense_fixed_o2_environment(TRUE)
  trajectories <- get("generate_analytical_trajectories", env)(
    run_dir = run_dir,
    time_points = time_points,
    o2_values = o2_values,
    initial_ploidy_values = initial_ploidies,
    analytical_methods = "eigen",
    n_workers = n_workers,
    seed_ids = seeds
  )
  attractors <- get("dense_attractor_mode_table", env)(run_dir, o2_values, seeds, n_workers)
  keep <- intersect(c(
    "seed_id", "O2_pct", "dominant_mean_N", "dominant_mean_ploidy",
    "dominant_growth_rate", "spectral_gap", "mode_label", "trajectory_regime",
    "objective", "delta_objective"
  ), names(attractors))
  merge(trajectories, attractors[, keep, drop = FALSE], by = c("seed_id", "O2_pct"), all.x = TRUE, sort = FALSE)
}

dense_run_task_slice <- function(argv, part, run_dir, out_dir, time_points, initial_ploidies) {
  paths <- dense_task_paths(out_dir, part)
  task_file <- normalizePath(path.expand(argv$task_file %||% paths$task_file), mustWork = TRUE)
  tasks <- dense_read_tsv(task_file)
  array_id <- dense_as_int(argv$array_task_id %||% Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"), 1L)
  tasks <- tasks[tasks$array_task_id == array_id, , drop = FALSE]
  if (!nrow(tasks)) stop("No dense-grid tasks for array_task_id=", array_id, call. = FALSE)
  dir.create(paths$chunk_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, , drop = FALSE]
    output <- file.path(paths$chunk_dir, sprintf("task_%08d.tsv", as.integer(task$task_id)))
    if (file.exists(output) && !dense_as_bool(argv$overwrite, TRUE)) next
    rows <- dense_generate_rows(part, run_dir, task$seed_id, as.numeric(task$O2_pct), time_points, initial_ploidies, 1L)
    rows$dense_grid_task_id <- as.integer(task$task_id)
    dense_write_tsv(rows, output)
  }
  message("Completed ", nrow(tasks), " ", part, " task(s) for array id ", array_id)
}

dense_materialized_path <- function(out_dir, part) file.path(dense_grid_simulation_dir(out_dir), if (part == "monotonicity") "fixed_o2_dense_grid_attractor_curves.tsv" else "fixed_o2_initial_ploidy_trajectories.tsv")

dense_write_materialized <- function(rows, out_dir, part, source) {
  path <- dense_materialized_path(out_dir, part)
  schema <- sub("[.]tsv$", ".schema.tsv", path)
  dense_write_tsv(rows, path)
  dense_write_schema(rows, paste0("dense_grid_", part, "_simulation"), schema)
  dense_record_artifact(out_dir, "simulation", paste0("dense_grid_", part, "_simulation"), path, rows, "generate_dense_grid_simulation_tables.R", source)
  message("Wrote materialized dense-grid simulation table: ", path)
  invisible(path)
}

dense_merge_chunks <- function(out_dir, part) {
  paths <- dense_task_paths(out_dir, part)
  chunks <- sort(list.files(paths$chunk_dir, pattern = "[.]tsv$", full.names = TRUE))
  if (!length(chunks)) stop("No task chunks under ", paths$chunk_dir, call. = FALSE)
  rows <- dense_rbind_fill(lapply(chunks, dense_read_tsv))
  if (nrow(rows)) rows <- rows[order(dense_seed_number(rows$seed_id), as.numeric(rows$O2_pct), if ("day" %in% names(rows)) as.numeric(rows$day) else seq_len(nrow(rows))), , drop = FALSE]
  dense_write_materialized(rows, out_dir, part, paste(chunks, collapse = ";"))
}

dense_simulation_main <- function(argv = dense_parse_args()) {
  part <- dense_grid_normalize_part(argv$part %||% argv$workflow_part)
  mode <- tolower(gsub("-", "_", argv$mode %||% argv$action %||% "run", fixed = TRUE))
  run_dir <- normalizePath(path.expand(argv$run_dir %||% argv$fit_root %||% dense_grid_default_fit_root()), mustWork = TRUE)
  out_dir <- normalizePath(path.expand(argv$out_dir %||% argv$output_root %||% dense_grid_default_out_dir(part)), mustWork = FALSE)
  o2_grid <- sort(unique(dense_as_num_vec(argv$o2_grid, seq(dense_as_num(argv$o2_min, 0), dense_as_num(argv$o2_max, 5), by = dense_as_num(argv$o2_by, 0.025)))))
  time_points <- sort(unique(dense_as_num_vec(argv$simulation_times %||% argv$time_points, dense_default_times())))
  initial_ploidies <- sort(unique(dense_as_num_vec(argv$initial_ploidy_values, c(2, 4))))
  if (mode == "build_tasks") return(dense_build_tasks(argv, part, run_dir, out_dir, o2_grid))
  if (mode == "run_tasks") return(dense_run_task_slice(argv, part, run_dir, out_dir, time_points, initial_ploidies))
  if (mode == "merge_daily_seed") {
    message("Per-seed daily materialization is already contained in dense-grid task chunks; no extra merge is required.")
    return(invisible(TRUE))
  }
  if (mode == "merge") return(dense_merge_chunks(out_dir, part))
  if (!mode %in% c("run", "generate", "all")) stop("Unknown simulation mode: ", mode, call. = FALSE)
  seeds <- dense_selected_seeds(run_dir, dense_as_int(argv$max_seeds, NA_integer_), argv$seed_manifest %||% argv$seed_manifest_file)
  rows <- dense_generate_rows(part, run_dir, seeds, o2_grid, time_points, initial_ploidies, max(1L, dense_as_int(argv$n_workers %||% argv$n_threads, 1L)))
  dense_write_materialized(rows, out_dir, part, run_dir)
}

if (sys.nframe() == 0L) dense_simulation_main()
