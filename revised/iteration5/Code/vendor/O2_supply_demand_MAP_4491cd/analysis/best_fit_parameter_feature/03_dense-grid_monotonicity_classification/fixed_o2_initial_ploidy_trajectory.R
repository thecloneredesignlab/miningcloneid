#!/usr/bin/env Rscript
# Deprecated compatibility entrypoint. Trajectories are materialized by simulation/o2/dense_grid.
.legacy_file <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "fixed_o2_initial_ploidy_trajectory.R"]
  if (length(own)) return(own[[length(own)]])
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  candidates <- c(getwd(), file.path(getwd(), "oxygen", "code", "O2_supply_demand_MAP", "analysis", "best_fit_parameter_feature", "03_dense-grid-monotonicity_classification"))
  hit <- candidates[file.exists(file.path(candidates, "deprecated_forward.R"))]
  if (!length(hit)) stop("Could not locate deprecated dense-grid wrapper directory.", call. = FALSE)
  file.path(hit[[1L]], "fixed_o2_initial_ploidy_trajectory.R")
})
source(file.path(dirname(.legacy_file), "deprecated_forward.R"), local = environment())
generate_outputs <- function(args = list()) deprecated_dense_run(args, "run", "initial_ploidy")
main <- generate_outputs
if (sys.nframe() == 0L) deprecated_dense_direct("run", "initial_ploidy")
