#!/usr/bin/env Rscript

.o2pl_analysis_entry_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "analyze_parameter_landscape.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(.o2pl_analysis_entry_dir, "parameter_landscape_analysis_utils.R"), local = environment(), chdir = TRUE)

o2pl_analysis_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  part <- tolower(trimws(as.character(argv$analysis_part %||% argv$part %||% "all")))
  part <- switch(part,
    invivo_reductions = "invivo",
    invitro_reductions = "invitro",
    pooled_reductions = "pooled",
    reductions = "all",
    part
  )
  if (!part %in% c("invivo", "invitro", "pooled", "all")) stop("Unknown analysis part: ", part, call. = FALSE)
  root_dir <- normalizePath(path.expand(argv$result_root %||% o2pl_default_result_root()), mustWork = FALSE)
  reductions <- unique(vapply(as_char_vec(argv$reductions, "umap"), o2pl_normalize_reduction, character(1)))
  preprocessing <- as_char_vec(argv$preprocess_modes %||% argv$preprocessing, "zscore")
  scopes <- if (part == "all") c("invivo", "invitro", "pooled") else part
  outputs <- list()
  for (scope in scopes) {
    for (method in preprocessing) {
      for (reduction in reductions) {
        key <- paste(scope, method, reduction, sep = "::")
        outputs[[key]] <- if (scope == "pooled") {
          o2pl_materialize_pooled_embedding(
            root_dir = root_dir,
            reduction = reduction,
            preprocessing = method,
            seed = as_int(argv$umap_seed %||% argv$tsne_seed, 123L),
            n_neighbors = as_int(argv$n_neighbors, 80L),
            min_dist = as_num(argv$min_dist, 0.1),
            n_threads = as_int(argv$n_threads, 1L),
            cluster_seed = as_int(argv$cluster_seed, 123L),
            cluster_k_min = as_int(argv$cluster_k_min, 2L),
            cluster_k_max = as_int(argv$cluster_k_max, 8L)
          )
        } else {
          o2pl_materialize_embedding(
            root_dir = root_dir,
            dataset = scope,
            reduction = reduction,
            preprocessing = method,
            seed = as_int(argv$umap_seed %||% argv$tsne_seed, 123L),
            n_neighbors = as_int(argv$n_neighbors, 80L),
            min_dist = as_num(argv$min_dist, 0.1),
            n_threads = as_int(argv$n_threads, 1L),
            cluster_seed = as_int(argv$cluster_seed, 123L),
            cluster_k_min = as_int(argv$cluster_k_min, 2L),
            cluster_k_max = as_int(argv$cluster_k_max, 8L)
          )
        }
      }
    }
  }
  message("Parameter-landscape analysis complete: ", root_dir)
  invisible(outputs)
}

if (sys.nframe() == 0L) o2pl_analysis_main()
