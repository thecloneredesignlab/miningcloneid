#!/usr/bin/env Rscript

.o2pl_vis_entry_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "visualize_parameter_landscape.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(.o2pl_vis_entry_dir, "parameter_landscape_visualization_utils.R"), local = environment(), chdir = TRUE)

o2pl_visualization_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  root_dir <- normalizePath(path.expand(argv$result_root %||% o2pl_default_result_root()), mustWork = FALSE)
  part <- tolower(trimws(as.character(argv$visualization_part %||% argv$analysis_part %||% argv$part %||% "all")))
  part <- switch(part,
    invivo_reductions = "invivo",
    invitro_reductions = "invitro",
    pooled_reductions = "pooled",
    reductions = "all",
    part
  )
  scopes <- if (part == "all") c("invivo", "invitro", "pooled") else part
  if (all(scopes %in% c("invivo", "invitro", "pooled"))) {
    reductions <- as_char_vec(argv$reductions, "umap")
    preprocessing <- as_char_vec(argv$preprocess_modes %||% argv$preprocessing, "zscore")
    for (scope in scopes) for (method in preprocessing) for (reduction in reductions) {
      o2pl_embedding_figure(
        root_dir = root_dir,
        scope = scope,
        reduction = reduction,
        preprocessing = method,
        initial_size = as_num(argv$initial_size, 0.22),
        fitted_size = as_num(argv$best_size %||% argv$fitted_size, 1.2)
      )
    }
  }
  if (part %in% c("mode_contribution", "contribution", "contributions")) {
    for (o2 in as_num_vec(argv$mode_reference_o2_values, as_num(argv$mode_reference_o2, 2))) {
      o2pl_contribution_figure(root_dir, argv$mode_contribution_target %||% "mode", o2, as_int(argv$top_n, 20L))
    }
  }
  if (part %in% c("fixed_o2", "fixed_o2_distribution")) o2pl_fixed_o2_distribution_figure(root_dir)
  message("Parameter-landscape visualization complete: ", root_dir)
  invisible(TRUE)
}

if (sys.nframe() == 0L) o2pl_visualization_main()
