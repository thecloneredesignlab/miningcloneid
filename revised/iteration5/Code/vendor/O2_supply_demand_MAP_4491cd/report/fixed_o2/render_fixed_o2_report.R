#!/usr/bin/env Rscript

# Assembly-only staged entrypoint for the existing fixed-O2 HTML renderer.

.fixo2_report_script_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      path <- env$ofile
      if (is.null(path)) "" else normalizePath(path, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[basename(frame_files) == "render_fixed_o2_report.R"]
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

fixo2_report_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      i <- i + 1L
      next
    }
    item <- substring(arg, 3L)
    pos <- regexpr("=", item, fixed = TRUE)
    if (pos > 0L) {
      out[[substring(item, 1L, pos - 1L)]] <- substring(item, pos + 1L)
      i <- i + 1L
    } else if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
      out[[item]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      out[[item]] <- TRUE
      i <- i + 1L
    }
  }
  out
}

fixo2_report_usage <- function() {
  paste(
    "Usage:",
    "  Rscript render_fixed_o2_report.R --results_dir DIR [options]",
    "",
    "Required:",
    "  --results_dir DIR      Unified directory containing materialized tables and figures.",
    "",
    "Options:",
    "  --out_dir DIR          HTML output directory (default: DIR/html_report).",
    "  --report_basename NAME Output basename (default: index).",
    "  --help                 Print this message.",
    sep = "\n"
  )
}

.fixo2_report_renderer_env <- new.env(parent = globalenv())
sys.source(
  file.path(.fixo2_report_script_dir, "..", "render_fixo2_invivo_report.R"),
  envir = .fixo2_report_renderer_env,
  chdir = TRUE
)

`%||%` <- .fixo2_report_renderer_env$o2sd_report_null_coalesce

fixo2_report_main <- function(args = fixo2_report_parse_args()) {
  if (isTRUE(args$help)) {
    cat(fixo2_report_usage(), "\n")
    return(invisible(NULL))
  }
  results_dir <- as.character(args$results_dir %||% args$analysis_dir %||% "")
  if (!nzchar(results_dir)) {
    stop("--results_dir is required.\n", fixo2_report_usage(), call. = FALSE)
  }
  results_dir <- normalizePath(results_dir, mustWork = TRUE)
  out_dir <- normalizePath(
    as.character(args$out_dir %||% file.path(results_dir, "html_report")),
    mustWork = FALSE
  )
  basename <- as.character(args$report_basename %||% "index")
  path <- .fixo2_report_renderer_env$write_html_report(
    analysis_dir = results_dir,
    out_dir = out_dir,
    report_basename = basename
  )
  message("Fixed-O2 HTML assembled from materialized artifacts: ", path)
  invisible(path)
}

if (sys.nframe() == 0L) {
  fixo2_report_main()
}
