#!/usr/bin/env Rscript

# Backward-compatible thin loader for the canonical ploidy-regime utilities.
# The historical analysis path remains source-compatible.

ploidy_regime_utils_loader_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files)) {
    dirname(frame_files[[length(frame_files)]])
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})
ploidy_regime_workflow_root <- normalizePath(
  file.path(ploidy_regime_utils_loader_dir, "..", ".."),
  mustWork = FALSE
)
source(
  file.path(
    ploidy_regime_workflow_root,
    "util",
    "o2_supply_demand_map_ploidy_regime_utils.R"
  ),
  local = TRUE
)
rm(
  ploidy_regime_utils_loader_dir,
  ploidy_regime_workflow_root
)
