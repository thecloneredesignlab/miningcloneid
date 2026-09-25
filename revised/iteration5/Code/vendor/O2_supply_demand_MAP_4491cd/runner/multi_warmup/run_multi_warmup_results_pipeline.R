#!/usr/bin/env Rscript

# Cross-stage collection pipeline: analysis tables, optional integrated
# extra-results execution, pure visualization, and report assembly.

.o2mw_results_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "run_multi_warmup_results_pipeline.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  workflow_frames <- frames[
    grepl("/O2_supply_demand_MAP/", frames, fixed = TRUE)
  ]
  if (length(workflow_frames)) {
    root <- dirname(workflow_frames[[length(workflow_frames)]])
    while (!identical(basename(root), "O2_supply_demand_MAP")) {
      parent <- dirname(root)
      if (identical(parent, root)) break
      root <- parent
    }
    if (identical(basename(root), "O2_supply_demand_MAP")) {
      return(file.path(root, "runner", "multi_warmup"))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
.o2mw_results_workflow_root <- normalizePath(file.path(.o2mw_results_runner_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_results_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

.o2mw_results_analysis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_results_workflow_root, "analysis", "multi_warmup", "collect_multi_warmup_tables.R"),
           envir = .o2mw_results_analysis, chdir = TRUE)
.o2mw_results_vis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_results_workflow_root, "vis", "multi_warmup", "render_multi_warmup_collected_figures.R"),
           envir = .o2mw_results_vis, chdir = TRUE)
find_valid_seed_dirs <- .o2mw_results_analysis$find_valid_seed_dirs
current_script_dir <- function() file.path(.o2mw_results_workflow_root, "analysis", "multi_warmup")

safe_component <- function(x) {
  x <- trimws(as.character(x %||% ""))
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

remove_generated_integrated_run <- function(integrated_run_dir, root) {
  integrated_run_dir <- normalizePath(integrated_run_dir, mustWork = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  if (!startsWith(integrated_run_dir, paste0(root, .Platform$file.sep))) {
    stop("Refusing to remove integrated run outside multi-warmup root: ", integrated_run_dir)
  }
  if (!identical(basename(integrated_run_dir), "multi_warmup_integrated_joint_run")) {
    stop("Refusing to remove unexpected integrated run path: ", integrated_run_dir)
  }
  if (dir.exists(integrated_run_dir)) unlink(integrated_run_dir, recursive = TRUE, force = TRUE)
}

create_integrated_seed_links <- function(manifest, root, out_dir) {
  integrated_run_dir <- file.path(out_dir, "multi_warmup_integrated_joint_run")
  remove_generated_integrated_run(integrated_run_dir, out_dir)
  dir.create(integrated_run_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    warmup_label <- safe_component(if ("warmup_label" %in% names(manifest)) manifest$warmup_label[[i]] else i)
    joint_prefix <- if ("joint_run_prefix" %in% names(manifest)) as.character(manifest$joint_run_prefix[[i]]) else paste0("fit_joint_", warmup_label)
    source_run_dir <- normalizePath(file.path(root, joint_prefix), mustWork = FALSE)
    seed_dirs <- find_valid_seed_dirs(source_run_dir)
    if (!length(seed_dirs)) next
    seed_dirs <- seed_dirs[order(seed_dirs)]
    for (seed_dir in seed_dirs) {
      raw_seed <- basename(seed_dir)
      integrated_seed <- paste0(warmup_label, "__", raw_seed)
      target_dir <- file.path(integrated_run_dir, integrated_seed)
      if (!file.symlink(normalizePath(seed_dir, mustWork = TRUE), target_dir)) {
        stop("Failed to create integrated seed symlink: ", target_dir, " -> ", seed_dir)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        warmup_label = warmup_label,
        joint_run_prefix = joint_prefix,
        raw_seed = raw_seed,
        integrated_seed = integrated_seed,
        source_seed_dir = normalizePath(seed_dir, mustWork = TRUE),
        integrated_seed_dir = target_dir,
        stringsAsFactors = FALSE
      )
    }
  }

  manifest_out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write_tsv(manifest_out, file.path(out_dir, "multi_warmup_integrated_seed_manifest.tsv"))
  list(run_dir = integrated_run_dir, seed_manifest = manifest_out)
}

run_integrated_extra_results <- function(integrated_run_dir) {
  if (!dir.exists(integrated_run_dir)) return(invisible(NULL))
  seed_dirs <- find_valid_seed_dirs(integrated_run_dir)
  if (!length(seed_dirs)) return(invisible(NULL))
  script_dir <- current_script_dir()
  extra_results_script <- normalizePath(file.path(script_dir, "..", "fit_results", "extra_results.R"), mustWork = FALSE)
  if (!file.exists(extra_results_script)) {
    extra_results_script <- normalizePath(file.path(getwd(), "oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results.R"), mustWork = FALSE)
  }
  if (!file.exists(extra_results_script)) {
    stop("extra_results.R was not found for integrated multi-warmup report generation.")
  }
  out_dir <- file.path(integrated_run_dir, "extra_results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  args <- c(
    paste0("--run_dir=", integrated_run_dir),
    paste0("--out_dir=", out_dir),
    "--allow_partial_seed_dirs=TRUE"
  )
  log_path <- file.path(integrated_run_dir, "integrated_extra_results_run.log")
  rscript <- file.path(R.home("bin"), "Rscript")
  message("Running integrated joint extra_results.R for ", length(seed_dirs), " prefixed seeds.")
  out <- suppressWarnings(system2(rscript, args = c(extra_results_script, args), stdout = TRUE, stderr = TRUE))
  writeLines(out, log_path, useBytes = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("Integrated extra_results.R failed with status ", status, ". See: ", log_path)
  }
  invisible(out_dir)
}

build_integrated_joint_extra_results <- function(manifest, root, out_dir, enabled = TRUE) {
  if (!isTRUE(enabled)) return(invisible(NULL))
  if (is.null(manifest) || !is.data.frame(manifest) || !nrow(manifest)) return(invisible(NULL))
  if (!("joint_run_prefix" %in% names(manifest))) return(invisible(NULL))
  integrated <- create_integrated_seed_links(manifest = manifest, root = root, out_dir = out_dir)
  if (!is.data.frame(integrated$seed_manifest) || !nrow(integrated$seed_manifest)) {
    message("No completed seed directories were available for integrated joint extra_results.")
    return(invisible(NULL))
  }
  run_integrated_extra_results(integrated$run_dir)
}


run_multi_warmup_results_pipeline <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  .o2mw_results_analysis$main(argv)
  root <- normalizePath(path.expand(as_chr(argv$multi_warmup_root, as_chr(argv$out_dir))), mustWork = TRUE)
  out_dir <- normalizePath(path.expand(as_chr(argv$out_dir, root)), mustWork = TRUE)
  manifest_path <- normalizePath(
    path.expand(as_chr(argv$manifest, file.path(root, "multi_warmup_manifest.tsv"))),
    mustWork = TRUE
  )
  manifest <- read_tsv(manifest_path)
  if (as_bool(argv$build_integrated_joint_results, TRUE)) {
    build_integrated_joint_extra_results(manifest, root, out_dir, enabled = TRUE)
  }
  if (as_bool(argv$render_figures, TRUE)) {
    .o2mw_results_vis$render_multi_warmup_collected_figures(out_dir)
  }
  invisible(TRUE)
}

main <- run_multi_warmup_results_pipeline
if (sys.nframe() == 0L) run_multi_warmup_results_pipeline()
