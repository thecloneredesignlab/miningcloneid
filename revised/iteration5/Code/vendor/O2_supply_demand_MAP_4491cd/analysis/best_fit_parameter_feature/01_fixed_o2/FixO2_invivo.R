#!/usr/bin/env Rscript

# Deprecated compatibility orchestrator for the historical
# best_fit_parameter_feature/01_fixed_o2 output location.

.fixo2_entry_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      path <- env$ofile
      if (is.null(path)) "" else normalizePath(path, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[basename(frame_files) == "FixO2_invivo.R"]
  if (length(own)) {
    dirname(own[[length(own)]])
  } else {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) {
      dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
    } else {
      normalizePath(getwd(), mustWork = FALSE)
    }
  }
})

.fixo2_loader_path <- normalizePath(
  file.path(
    .fixo2_entry_dir,
    "..",
    "..",
    "..",
    "runner",
    "fixed_o2",
    "fixed_o2_pipeline_loader.R"
  ),
  mustWork = TRUE
)
sys.source(.fixo2_loader_path, envir = environment(), chdir = TRUE)

.fixo2_legacy_main <- fixo2_main
fixo2_main <- function(args = o2ipa_parse_args()) {
  if (is.null(args$out_dir) || !length(args$out_dir) || !nzchar(as.character(args$out_dir[[1L]]))) {
    repo_root <- o2ipa_repo_root(SCRIPT_DIR)
    args$out_dir <- file.path(
      repo_root,
      "oxygen",
      "results",
      "analysis",
      "best_fit_parameter_feature",
      "01_fixed_o2",
      "FixO2_invivo_500seed"
    )
  }
  message(
    fixo2_deprecation_message(
      "analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo.R"
    )
  )
  .fixo2_legacy_main(args)
}

if (sys.nframe() == 0L) {
  fixo2_main()
}
