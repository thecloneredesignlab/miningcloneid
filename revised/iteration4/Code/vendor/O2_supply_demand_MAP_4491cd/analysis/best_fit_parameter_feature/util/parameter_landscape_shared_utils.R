#!/usr/bin/env Rscript

# Deprecated compatibility loader.
#
# Temporary bridge for fixed-O2 eigen-attractor and warm-up callers.  New code
# must source its owning canonical layer directly.  This loader intentionally
# contains no model, analysis, or visualization implementation.
.o2pl_legacy_root <- local({
  starts <- unique(c(normalizePath(getwd(), mustWork = FALSE), Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else dirname(normalizePath(env$ofile, mustWork = FALSE))
  }, character(1)))))
  for (start in starts) {
    current <- start
    for (i in seq_len(8L)) {
      if (dir.exists(file.path(current, "simulation", "parameter_landscape")) && dir.exists(file.path(current, "analysis", "parameter_landscape_clustering"))) {
        return(normalizePath(current, mustWork = TRUE))
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Cannot locate O2_supply_demand_MAP for the deprecated parameter-landscape compatibility loader.", call. = FALSE)
})
source(file.path(.o2pl_legacy_root, "simulation", "parameter_landscape", "parameter_landscape_simulation_utils.R"), local = environment(), chdir = TRUE)
source(file.path(.o2pl_legacy_root, "simulation", "parameter_landscape", "parameter_landscape_invivo_feature_simulation.R"), local = environment(), chdir = TRUE)
source(file.path(.o2pl_legacy_root, "analysis", "parameter_landscape_clustering", "parameter_landscape_analysis_utils.R"), local = environment(), chdir = TRUE)
