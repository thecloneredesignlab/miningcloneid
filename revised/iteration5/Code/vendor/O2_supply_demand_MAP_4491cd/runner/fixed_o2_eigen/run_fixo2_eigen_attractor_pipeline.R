#!/usr/bin/env Rscript

# Cross-stage FixO2 eigen-attractor pipeline orchestrator.

.fixo2ea_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "run_fixo2_eigen_attractor_pipeline.R"]
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
      return(file.path(root, "runner", "fixed_o2_eigen"))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else normalizePath(getwd(), mustWork = FALSE)
})
.fixo2ea_runner_workflow_root <- normalizePath(file.path(.fixo2ea_runner_dir, "..", ".."), mustWork = TRUE)

sys.source(
  file.path(.fixo2ea_runner_workflow_root, "simulation", "o2", "fixed_o2", "eigen_attractor", "build_fixo2_eigen_attractor_features.R"),
  envir = environment(), chdir = TRUE
)
sys.source(
  file.path(.fixo2ea_runner_workflow_root, "analysis", "fixed_o2_eigen", "analyze_fixo2_eigen_attractor_embeddings.R"),
  envir = environment(), chdir = TRUE
)
sys.source(
  file.path(.fixo2ea_runner_workflow_root, "vis", "fixed_o2_eigen", "render_fixo2_eigen_attractor_figures.R"),
  envir = environment(), chdir = TRUE
)
sys.source(
  file.path(.fixo2ea_runner_workflow_root, "report", "fixed_o2_eigen", "render_fixo2_eigen_attractor_report.R"),
  envir = environment(), chdir = TRUE
)

fixo2ea_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  result_root <- normalizePath(path.expand(argv$result_root %||% fixo2ea_default_result_root()), mustWork = FALSE)
  source_root <- normalizePath(path.expand(argv$source_root %||% fixo2ea_default_parameter_source_root()), mustWork = FALSE)
  run_parts <- tolower(as_char_vec(argv$run_parts, c("features", "reductions", "figures", "report")))
  datasets <- as_char_vec(argv$datasets, c("invivo", "invitro"))
  reductions <- as_char_vec(argv$reductions, c("pca", "umap", "tsne"))
  o2_values <- fixo2ea_o2_grid(
    o2_values = argv$o2_values,
    o2_min = argv$o2_min %||% 0,
    o2_max = argv$o2_max %||% 5,
    o2_n = argv$o2_n %||% 201
  )
  if ("features" %in% run_parts) {
    fixo2ea_build_feature_tables(
      result_root = result_root,
      source_root = source_root,
      invivo_input = argv$invivo_input %||% default_dataset_input_dir("invivo"),
      invitro_input = argv$invitro_input %||% default_dataset_input_dir("invitro"),
      invivo_best_csv = argv$invivo_best_csv %||% fixo2ea_default_best_csv("invivo", source_root),
      invitro_best_csv = argv$invitro_best_csv %||% fixo2ea_default_best_csv("invitro", source_root),
      invivo_initial_csv = argv$invivo_initial_csv %||% fixo2ea_default_initial_csv("invivo", source_root),
      invitro_initial_csv = argv$invitro_initial_csv %||% fixo2ea_default_initial_csv("invitro", source_root),
      o2_values = o2_values,
      datasets = datasets,
      point_types = as_char_vec(argv$point_types, c("best", "initial")),
      seeds = argv$seeds,
      max_seeds = as_int(argv$max_seeds, NA_integer_),
      n_workers = as_int(argv$n_workers, 1L),
      chunk_size = as_int(argv$chunk_size, 250L),
      write_initial_long = as_bool(argv$write_initial_long, FALSE),
      force = as_bool(argv$force, FALSE)
    )
  }
  if ("reductions" %in% run_parts || "analysis" %in% run_parts) {
    fixo2ea_run_all_reductions(
      result_root = result_root,
      datasets = datasets,
      reductions = reductions,
      run_single = as_bool(argv$run_single, TRUE),
      run_pooled = as_bool(argv$run_pooled, TRUE),
      run_initial_vs_best = as_bool(argv$run_initial_vs_best, TRUE),
      run_best_only = as_bool(argv$run_best_only, TRUE),
      run_clustered = as_bool(argv$run_clustered, TRUE),
      umap_seed = as_int(argv$umap_seed, 123L),
      n_neighbors = as_int(argv$n_neighbors, 80L),
      min_dist = as_num(argv$min_dist, 0.1),
      n_threads = as_int(argv$n_threads, 1L),
      tsne_seed = as_int(argv$tsne_seed, 123L),
      tsne_perplexity = as_num(argv$tsne_perplexity, 30),
      tsne_theta = as_num(argv$tsne_theta, 0.5),
      tsne_max_iter = as_int(argv$tsne_max_iter, 1000L),
      cluster_seed = as_int(argv$cluster_seed, 123L),
      cluster_k_min = as_int(argv$cluster_k_min, 2L),
      cluster_k_max = as_int(argv$cluster_k_max, 8L),
      cluster_silhouette_sample_n = as_int(argv$cluster_silhouette_sample_n, 5000L)
    )
  }
  if ("figures" %in% run_parts || "vis" %in% run_parts) {
    fixo2ea_render_all_figures(
      result_root = result_root,
      width = as_num(argv$figure_width, 6.5),
      height = as_num(argv$figure_height, 6.2),
      dpi = as_int(argv$figure_dpi, 220L)
    )
  }
  if ("report" %in% run_parts) fixo2ea_write_report(result_root = result_root)
  invisible(TRUE)
}

if (sys.nframe() == 0L) fixo2ea_main()
