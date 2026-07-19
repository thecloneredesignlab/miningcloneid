#!/usr/bin/env Rscript

# Deprecated compatibility wrapper for the report-only canonical module.

.fixo2ea_report_legacy_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1L)))
  own <- frames[basename(frames) == "fixo2_eigen_attractor_report.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
.fixo2ea_report_workflow_root <- normalizePath(file.path(.fixo2ea_report_legacy_dir, "..", "..", ".."), mustWork = TRUE)
.fixo2ea_report_canonical <- file.path(
  .fixo2ea_report_workflow_root, "report", "fixed_o2_eigen",
  "render_fixo2_eigen_attractor_report.R"
)
sys.source(.fixo2ea_report_canonical, envir = environment(), chdir = TRUE)
rm(.fixo2ea_report_legacy_dir, .fixo2ea_report_workflow_root, .fixo2ea_report_canonical)
