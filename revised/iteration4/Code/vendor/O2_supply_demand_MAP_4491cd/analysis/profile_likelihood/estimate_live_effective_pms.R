#!/usr/bin/env Rscript

# Deprecated compatibility wrapper.
#
# The numerical producer moved to:
#   simulation/invivo/cin/generate_live_effective_pms_outputs.R
#
# This wrapper forwards all CLI arguments unchanged. The canonical producer
# requires an existing in-vivo simulation manifest and tables.

.o2sd_live_pms_wrapper_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2sd_live_pms_wrapper_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  producer <- normalizePath(
    file.path(
      WORKFLOW_ROOT,
      "simulation",
      "invivo",
      "cin",
      "generate_live_effective_pms_outputs.R"
    ),
    mustWork = TRUE
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = shQuote(c(producer, argv))
  )
  if (!identical(status, 0L)) {
    stop(
      "Canonical live-effective-p_ms simulation producer failed with status ",
      status
    )
  }
}

if (sys.nframe() == 0) {
  main()
}
