#!/usr/bin/env Rscript

# Cross-stage orchestration for warm-up joint prepare, embedding, curve-class,
# summary, and visualization stages.

.o2mw_warm_runner_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "run_warm_up_joint_results_pipeline.R"]
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
.o2mw_warm_workflow_root <- normalizePath(file.path(.o2mw_warm_runner_dir, "..", ".."), mustWork = TRUE)
source(file.path(.o2mw_warm_workflow_root, "util", "o2_supply_demand_map_multi_warmup_utils.R"), local = environment())

.o2mw_warm_analysis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_warm_workflow_root, "analysis", "multi_warmup", "analyze_warm_up_joint_results.R"),
           envir = .o2mw_warm_analysis, chdir = TRUE)
.o2mw_warm_vis <- new.env(parent = globalenv())
sys.source(file.path(.o2mw_warm_workflow_root, "vis", "multi_warmup", "render_warm_up_joint_figures.R"),
           envir = .o2mw_warm_vis, chdir = TRUE)

run_warm_up_curve_classification <- function(args) {
  input_root <- normalizePath(path.expand(args$input_root %||% .o2mw_warm_analysis$default_input_root()), mustWork = FALSE)
  out_root <- normalizePath(path.expand(args$output_root %||% .o2mw_warm_analysis$default_output_root(input_root)), mustWork = FALSE)
  curve_manifest <- file.path(out_root, "tables", "joint_best_curve_seed_manifest.tsv")
  if (!file.exists(curve_manifest)) stop("Missing curve seed manifest. Run prepare first: ", curve_manifest, call. = FALSE)
  dense_dir <- .o2mw_warm_analysis$DENSE_DIR
  mono_env <- .o2mw_warm_analysis$source_env(file.path(dense_dir, "fixed_o2_ploidy_monotonicity.R"))
  dense_out <- file.path(out_root, "curve_classification", "dense-grid_monotonicity_classification")
  o2_grid <- .o2mw_warm_analysis$as_num_vec(args$o2_grid, seq(0, 5, by = 0.025))
  reporting_o2 <- .o2mw_warm_analysis$as_num_vec(args$reporting_o2, c(0, 0.1, 0.5, 1, 2, 5))
  mono_env$generate_outputs(list(
    run_dir = input_root,
    seed_manifest = curve_manifest,
    out_dir = dense_out,
    o2_grid = paste(format(o2_grid, scientific = FALSE, trim = TRUE), collapse = ","),
    reporting_o2 = paste(format(reporting_o2, scientific = FALSE, trim = TRUE), collapse = ","),
    n_workers = as.character(as_int(args$n_workers, 8L)),
    max_seeds = as.character(as_int(args$curve_max_seeds, NA_integer_)),
    overwrite = as.character(as_bool(args$overwrite, TRUE)),
    generate_figures = "FALSE",
    run_validation = as.character(as_bool(args$run_validation, TRUE))
  ))
  run_warm_up_curve_regression(args)
}

run_warm_up_curve_regression <- function(args) {
  input_root <- normalizePath(path.expand(args$input_root %||% .o2mw_warm_analysis$default_input_root()), mustWork = FALSE)
  out_root <- normalizePath(path.expand(args$output_root %||% .o2mw_warm_analysis$default_output_root(input_root)), mustWork = FALSE)
  dense_out <- file.path(out_root, "curve_classification", "dense-grid_monotonicity_classification")
  reg_out <- file.path(out_root, "curve_classification", "dense-grid_monotonicity_regression_classification")
  required <- file.path(dense_out, "tables", c(
    "fixed_o2_ploidy_monotonicity_by_seed.tsv",
    "fixed_o2_ploidy_monotonicity_curves.tsv"
  ))
  if (!all(file.exists(required))) stop("Missing dense-grid monotonicity tables under: ", dense_out, call. = FALSE)
  reg_env <- .o2mw_warm_analysis$source_env(file.path(
    .o2mw_warm_analysis$DENSE_DIR,
    "fixed_o2_ploidy_monotonicity_regression_classification.R"
  ))
  reg_env$generate_outputs(list(
    input_dir = dense_out,
    out_dir = reg_out,
    overwrite = as.character(as_bool(args$overwrite, TRUE)),
    generate_figures = "FALSE"
  ))
  invisible(list(pointwise = dense_out, regression = reg_out))
}

run_warm_up_joint_results_pipeline <- function(argv = .o2mw_warm_analysis$parse_args()) {
  stages <- .o2mw_warm_analysis$stage_values(argv$stage %||% "all")
  if ("prepare" %in% stages) .o2mw_warm_analysis$build_prepare_outputs(argv)
  if ("embedding" %in% stages) .o2mw_warm_analysis$build_embedding_outputs(argv)
  if ("curve" %in% stages) run_warm_up_curve_classification(argv)
  if ("curve_regression" %in% stages && !"curve" %in% stages) run_warm_up_curve_regression(argv)
  if ("summary" %in% stages) .o2mw_warm_analysis$build_summary_outputs(argv)
  if (as_bool(argv$render_figures, TRUE) && any(stages %in% c("embedding", "summary"))) {
    input_root <- normalizePath(path.expand(argv$input_root %||% .o2mw_warm_analysis$default_input_root()), mustWork = FALSE)
    out_root <- normalizePath(path.expand(argv$output_root %||% .o2mw_warm_analysis$default_output_root(input_root)), mustWork = TRUE)
    .o2mw_warm_vis$render_warm_up_joint_figures(out_root)
  }
  invisible(TRUE)
}

main <- run_warm_up_joint_results_pipeline
if (sys.nframe() == 0L) run_warm_up_joint_results_pipeline()
