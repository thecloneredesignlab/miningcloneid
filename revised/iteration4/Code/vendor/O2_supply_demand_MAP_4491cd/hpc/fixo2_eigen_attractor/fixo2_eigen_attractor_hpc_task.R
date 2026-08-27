#!/usr/bin/env Rscript

hpc_task_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "fixo2_eigen_attractor_hpc_task.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

workflow_root <- normalizePath(
  file.path(hpc_task_dir, "..", ".."),
  mustWork = TRUE
)
sys.source(
  file.path(
    workflow_root,
    "simulation", "o2", "fixed_o2", "eigen_attractor",
    "build_fixo2_eigen_attractor_features.R"
  ),
  envir = environment(),
  chdir = TRUE
)

fixo2ea_hpc_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  mode <- tolower(trimws(as.character(argv$mode %||% "run_task")))
  result_root <- normalizePath(path.expand(argv$result_root %||% fixo2ea_default_result_root()), mustWork = FALSE)
  source_root <- normalizePath(path.expand(argv$source_root %||% fixo2ea_default_parameter_source_root()), mustWork = FALSE)
  invivo_input <- normalizePath(path.expand(argv$invivo_input %||% default_dataset_input_dir("invivo")), mustWork = FALSE)
  invitro_input <- normalizePath(path.expand(argv$invitro_input %||% default_dataset_input_dir("invitro")), mustWork = FALSE)
  task_table <- normalizePath(path.expand(argv$task_table %||% fixo2ea_task_table_path(result_root)), mustWork = FALSE)
  o2_values <- fixo2ea_o2_grid(
    o2_values = argv$o2_values,
    o2_min = argv$o2_min %||% 0,
    o2_max = argv$o2_max %||% 5,
    o2_n = argv$o2_n %||% 201
  )

  if (identical(mode, "build_tasks")) {
    fixo2ea_build_task_table(
      result_root = result_root,
      source_root = source_root,
      invivo_input = invivo_input,
      invitro_input = invitro_input,
      invivo_best_csv = argv$invivo_best_csv %||% fixo2ea_default_best_csv("invivo", source_root),
      invitro_best_csv = argv$invitro_best_csv %||% fixo2ea_default_best_csv("invitro", source_root),
      invivo_initial_csv = argv$invivo_initial_csv %||% fixo2ea_default_initial_csv("invivo", source_root),
      invitro_initial_csv = argv$invitro_initial_csv %||% fixo2ea_default_initial_csv("invitro", source_root),
      datasets = as_char_vec(argv$datasets, c("invivo", "invitro")),
      point_types = as_char_vec(argv$point_types, c("best", "initial")),
      seeds = argv$seeds,
      max_seeds = as_int(argv$max_seeds, NA_integer_),
      force = as_bool(argv$force, FALSE)
    )
    return(invisible(TRUE))
  }

  if (identical(mode, "run_task")) {
    task_ids <- fixo2ea_parse_seed_selector(argv$task_ids)
    if (!length(task_ids)) {
      array_task_id <- argv$array_task_id %||% Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_)
      max_task_id <- as_int(argv$max_task_id, NA_integer_)
      if (!is.finite(max_task_id) || is.na(max_task_id)) {
        max_task_id <- nrow(fixo2ea_read_tsv_plain(task_table))
      }
      task_ids <- fixo2ea_array_task_ids(
        array_task_id = array_task_id,
        points_per_task = as_int(argv$points_per_task, 1L),
        max_task_id = max_task_id
      )
    }
    if (!length(task_ids)) {
      message("No task ids selected for this array task; exiting.")
      return(invisible(TRUE))
    }
    fixo2ea_run_task_ids(
      task_table = task_table,
      result_root = result_root,
      task_ids = task_ids,
      invivo_input = invivo_input,
      invitro_input = invitro_input,
      o2_values = o2_values,
      force = as_bool(argv$force, FALSE)
    )
    return(invisible(TRUE))
  }

  if (identical(mode, "merge_tasks")) {
    fixo2ea_merge_task_results(
      result_root = result_root,
      task_table = task_table,
      force = as_bool(argv$force, FALSE)
    )
    return(invisible(TRUE))
  }

  stop("--mode must be one of: build_tasks, run_task, merge_tasks.", call. = FALSE)
}

if (sys.nframe() == 0L) fixo2ea_hpc_main()
