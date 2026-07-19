#!/usr/bin/env Rscript

.dense_entry_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "analyze_dense_grid_tables.R"]
  if (length(own)) dirname(own[[length(own)]]) else dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = FALSE))
})
source(file.path(.dense_entry_dir, "dense_grid_analysis_utils.R"), local = environment(), chdir = TRUE)

dense_analysis_main <- function(argv = dense_parse_args()) {
  part <- dense_grid_normalize_part(argv$part %||% argv$workflow_part)
  out_dir <- normalizePath(path.expand(argv$out_dir %||% argv$output_root %||% dense_grid_default_out_dir(part)), mustWork = FALSE)
  if (part == "monotonicity") dense_materialize_monotonicity_analysis(out_dir, argv) else dense_materialize_initial_analysis(out_dir, argv)
  message("Dense-grid analysis complete: ", out_dir)
}

if (sys.nframe() == 0L) dense_analysis_main()
