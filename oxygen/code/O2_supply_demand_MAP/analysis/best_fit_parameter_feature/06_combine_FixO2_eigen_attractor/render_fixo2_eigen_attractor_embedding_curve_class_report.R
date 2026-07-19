#!/usr/bin/env Rscript

# Deprecated compatibility wrapper. Canonical report implementation:
# report/combined_fixo2_eigen/render_fixo2_eigen_attractor_embedding_curve_class_report.R
.wrapper_file <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "render_fixo2_eigen_attractor_embedding_curve_class_report.R"]
  if (length(own)) return(own[[length(own)]])
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  }
  normalizePath(
    file.path(getwd(), "render_fixo2_eigen_attractor_embedding_curve_class_report.R"),
    mustWork = FALSE
  )
})
.canonical <- normalizePath(
  file.path(dirname(.wrapper_file), "..", "..", "..", "report", "combined_fixo2_eigen",
            "render_fixo2_eigen_attractor_embedding_curve_class_report.R"),
  mustWork = TRUE
)
source(.canonical, local = environment(), chdir = TRUE)
rm(.wrapper_file, .canonical)
