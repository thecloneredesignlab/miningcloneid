#!/usr/bin/env Rscript

# Dedicated entrypoint for the historical standalone in vitro extra-results
# artifact set. The generic staged runner remains the canonical joint/in vivo
# workflow and is intentionally not changed by this compatibility entrypoint.
.script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    getwd()
  }
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
reference_script <- file.path(
  WORKFLOW_ROOT,
  "analysis",
  "fit_results",
  "extra_results_invitro_reference.R"
)

if (sys.nframe() == 0L) {
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(normalizePath(reference_script, mustWork = TRUE), commandArgs(trailingOnly = TRUE))
  )
  if (!identical(status, 0L)) quit(save = "no", status = status)
}
