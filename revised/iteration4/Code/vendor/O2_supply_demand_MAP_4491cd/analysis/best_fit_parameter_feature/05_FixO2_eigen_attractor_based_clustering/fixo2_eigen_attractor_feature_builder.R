#!/usr/bin/env Rscript

# Deprecated compatibility wrapper.  Numerical feature production now lives
# under simulation/o2/fixed_o2/eigen_attractor/.

.fixo2ea_legacy_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "fixo2_eigen_attractor_feature_builder.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
})

.fixo2ea_workflow_root <- normalizePath(file.path(.fixo2ea_legacy_dir, "..", "..", ".."), mustWork = FALSE)
.fixo2ea_builder <- file.path(
  .fixo2ea_workflow_root,
  "simulation", "o2", "fixed_o2", "eigen_attractor",
  "build_fixo2_eigen_attractor_features.R"
)
if (!file.exists(.fixo2ea_builder)) stop("Missing canonical FixO2 eigen builder: ", .fixo2ea_builder, call. = FALSE)
sys.source(.fixo2ea_builder, envir = environment(), chdir = TRUE)
rm(.fixo2ea_legacy_dir, .fixo2ea_workflow_root, .fixo2ea_builder)
