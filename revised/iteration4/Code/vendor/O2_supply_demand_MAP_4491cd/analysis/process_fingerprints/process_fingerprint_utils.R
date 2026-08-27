#!/usr/bin/env Rscript

# Backward-compatible thin loader for the canonical process-fingerprint
# utilities. The historical analysis path remains source-compatible.

process_fingerprint_utils_loader_dir <- local({
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
process_fingerprint_workflow_root <- normalizePath(
  file.path(process_fingerprint_utils_loader_dir, "..", ".."),
  mustWork = FALSE
)
source(
  file.path(
    process_fingerprint_workflow_root,
    "util",
    "o2_supply_demand_map_process_fingerprint_utils.R"
  ),
  local = TRUE
)
rm(
  process_fingerprint_utils_loader_dir,
  process_fingerprint_workflow_root
)
