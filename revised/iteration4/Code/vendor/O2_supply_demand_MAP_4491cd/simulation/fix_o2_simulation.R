#!/usr/bin/env Rscript

# Backward-compatible entrypoint. The canonical implementation is organized
# under simulation/o2/fixed_o2/.

.legacy_fix_o2_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
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
  own <- frame_files[basename(frame_files) == "fix_o2_simulation.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  getwd()
})

.canonical_fix_o2 <- file.path(
  .legacy_fix_o2_dir,
  "o2",
  "fixed_o2",
  "run_fixed_o2_simulation.R"
)
if (!file.exists(.canonical_fix_o2)) {
  stop("Missing canonical fixed-O2 simulator: ", .canonical_fix_o2)
}
sys.source(.canonical_fix_o2, envir = environment(), chdir = TRUE)
rm(.legacy_fix_o2_dir, .canonical_fix_o2)

if (sys.nframe() == 0L) {
  main()
}
