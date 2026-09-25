#!/usr/bin/env Rscript

.dense_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "run_dense_grid_monotonicity.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
.dense_runner_root <- normalizePath(file.path(.dense_runner_dir, "..", ".."), mustWork = TRUE)
source(file.path(.dense_runner_root, "util", "o2_supply_demand_map_dense_grid_utils.R"), local = environment(), chdir = TRUE)

dense_runner_arg <- function(key, value) {
  if (is.null(value) || !length(value) || all(is.na(value))) character() else paste0("--", key, "=", paste(value, collapse = ","))
}
dense_runner_child <- function(label, script, args, dry_run = FALSE) {
  command <- c("--vanilla", normalizePath(script, mustWork = FALSE), args)
  message("[", label, "] ", paste(shQuote(c(file.path(R.home("bin"), "Rscript"), command)), collapse = " "))
  if (dry_run) return(invisible(0L))
  status <- system2(file.path(R.home("bin"), "Rscript"), command)
  if (!identical(status, 0L)) stop(label, " failed with status ", status, call. = FALSE)
  invisible(status)
}
dense_runner_replace_arg <- function(raw, key, value) c(raw[!grepl(paste0("^--", key, "="), raw)], dense_runner_arg(key, value))

dense_runner_postprocess <- function(part, argv, raw_args, dry_run = FALSE,
                                     run_analysis = TRUE, run_visualization = TRUE,
                                     run_report = TRUE) {
  out_dir <- normalizePath(path.expand(argv$out_dir %||% argv$output_root %||% dense_grid_default_out_dir(part)), mustWork = FALSE)
  args <- dense_runner_replace_arg(raw_args, "part", part)
  args <- dense_runner_replace_arg(args, "out_dir", out_dir)
  analysis <- file.path(.dense_runner_root, "analysis", "dense_grid_monotonicity", "analyze_dense_grid_tables.R")
  visualization <- file.path(.dense_runner_root, "vis", "dense_grid_monotonicity", "render_dense_grid_figures.R")
  report <- file.path(.dense_runner_root, "report", "dense_grid_monotonicity", "render_dense_grid_report.R")
  if (run_analysis) dense_runner_child(paste("analyze", part), analysis, args, dry_run)
  if (run_visualization && dense_as_bool(argv$generate_figures, TRUE)) dense_runner_child(paste("visualize", part), visualization, args, dry_run)
  if (run_report && dense_as_bool(argv$run_report, TRUE)) dense_runner_child(paste("assemble", part, "report"), report, args, dry_run)
  invisible(TRUE)
}

dense_runner_one_part <- function(part, argv, raw_args, dry_run = FALSE, simulation_mode = "run") {
  out_dir <- normalizePath(path.expand(argv$out_dir %||% argv$output_root %||% dense_grid_default_out_dir(part)), mustWork = FALSE)
  args <- dense_runner_replace_arg(raw_args, "part", part)
  args <- dense_runner_replace_arg(args, "out_dir", out_dir)
  simulation <- file.path(.dense_runner_root, "simulation", "o2", "dense_grid", "generate_dense_grid_simulation_tables.R")
  dense_runner_child(paste("materialize", part, "simulation"), simulation, dense_runner_replace_arg(args, "mode", simulation_mode), dry_run)
  if (simulation_mode %in% c("build_tasks", "run_tasks", "merge_daily_seed")) return(invisible(TRUE))
  dense_runner_postprocess(part, argv, args, dry_run)
}

dense_runner_main <- function(argv = dense_parse_args(), raw_args = commandArgs(trailingOnly = TRUE)) {
  dry_run <- dense_as_bool(argv$dry_run, FALSE)
  mode <- tolower(gsub("-", "_", argv$mode %||% argv$action %||% "run", fixed = TRUE))
  if (mode %in% c("build_tasks", "run_tasks", "merge_daily_seed", "merge")) {
    return(dense_runner_one_part(dense_grid_normalize_part(argv$part %||% argv$workflow_part), argv, raw_args, dry_run, mode))
  }
  if (mode %in% c("analyze", "analysis", "postprocess", "visualize", "visualise", "report")) {
    part <- dense_grid_normalize_part(argv$part %||% argv$workflow_part)
    return(dense_runner_postprocess(
      part, argv, raw_args, dry_run,
      run_analysis = mode %in% c("analyze", "analysis", "postprocess"),
      run_visualization = mode %in% c("postprocess", "visualize", "visualise"),
      run_report = mode %in% c("postprocess", "report")
    ))
  }
  parts <- tolower(dense_as_char_vec(argv$run_parts %||% argv$part, "all"))
  if ("all" %in% parts) parts <- c("monotonicity", "initial_ploidy")
  parts <- unique(vapply(parts, dense_grid_normalize_part, character(1)))
  for (part in parts) dense_runner_one_part(part, argv, raw_args, dry_run, "run")
  message("Dense-grid runner complete.")
}

if (sys.nframe() == 0L) dense_runner_main()
