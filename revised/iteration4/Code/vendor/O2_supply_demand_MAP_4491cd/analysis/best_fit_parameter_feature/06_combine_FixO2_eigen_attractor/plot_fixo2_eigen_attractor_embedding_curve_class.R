#!/usr/bin/env Rscript

# Deprecated compatibility orchestrator. It materializes analysis tables first,
# then invokes the consume-only canonical visualization implementation.
.wrapper_file <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "plot_fixo2_eigen_attractor_embedding_curve_class.R"]
  if (length(own)) return(own[[length(own)]])
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  }
  normalizePath(
    file.path(getwd(), "plot_fixo2_eigen_attractor_embedding_curve_class.R"),
    mustWork = FALSE
  )
})
.workflow_dir <- normalizePath(
  file.path(dirname(.wrapper_file), "..", "..", ".."),
  mustWork = TRUE
)
.analysis_script <- file.path(
  .workflow_dir,
  "analysis", "combined_fixo2_eigen", "prepare_fixo2_eigen_curve_class_tables.R"
)
.vis_script <- file.path(
  .workflow_dir,
  "vis", "combined_fixo2_eigen", "plot_fixo2_eigen_attractor_embedding_curve_class.R"
)

# Preserve the legacy source-level function surface without executing either CLI.
.target_env <- environment()
for (.script in c(.analysis_script, .vis_script)) {
  .api_env <- new.env(parent = .target_env)
  sys.source(.script, envir = .api_env, chdir = TRUE)
  .functions <- ls(.api_env)[vapply(ls(.api_env), function(name) {
    is.function(get(name, envir = .api_env, inherits = FALSE))
  }, logical(1))]
  list2env(mget(.functions, envir = .api_env, inherits = FALSE), envir = .target_env)
}

main <- function(raw_args = commandArgs(trailingOnly = TRUE)) {
  quoted <- vapply(c(.analysis_script, raw_args), shQuote, character(1))
  status <- system2("Rscript", args = quoted)
  if (!identical(as.integer(status), 0L)) {
    stop("Combined FixO2 eigen-attractor analysis failed with exit status ", status, call. = FALSE)
  }
  parsed <- bpf_parse_args(raw_args)
  if (bpf_as_bool(parsed$dry_run, FALSE)) return(invisible(status))
  quoted <- vapply(c(.vis_script, raw_args), shQuote, character(1))
  status <- system2("Rscript", args = quoted)
  if (!identical(as.integer(status), 0L)) {
    stop("Combined FixO2 eigen-attractor visualization failed with exit status ", status, call. = FALSE)
  }
  invisible(status)
}

rm(.api_env, .functions, .script, .target_env)
if (identical(environment(), globalenv())) main()
