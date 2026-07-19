#!/usr/bin/env Rscript

# Deprecated compatibility helpers for the pre-stage-3 dense-grid entrypoints.
# New code must call runner/dense_grid_monotonicity/run_dense_grid_monotonicity.R.

.deprecated_dense_file <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "deprecated_forward.R"]
  if (length(own)) own[[length(own)]] else ""
})

deprecated_dense_root <- function() {
  current <- dirname(.deprecated_dense_file)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(current, "runner", "dense_grid_monotonicity", "run_dense_grid_monotonicity.R"))) return(normalizePath(current, mustWork = TRUE))
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate O2_supply_demand_MAP from deprecated dense-grid entrypoint.", call. = FALSE)
}

deprecated_dense_flags <- function(args) {
  if (!length(args)) return(character())
  unlist(Map(function(name, value) {
    if (is.null(value) || !length(value) || all(is.na(value))) return(character())
    paste0("--", name, "=", paste(value, collapse = ","))
  }, names(args), args), use.names = FALSE)
}

deprecated_dense_parse_args <- function(raw = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in raw[startsWith(raw, "--")]) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[gsub("-", "_", pair[[1L]], fixed = TRUE)]] <- if (length(pair) > 1L) paste(pair[-1L], collapse = "=") else "TRUE"
  }
  out
}

deprecated_dense_run <- function(args = list(), mode = "run", part = "monotonicity") {
  if (is.null(args$out_dir)) args$out_dir <- args$result_dir %||% args$output_root
  if (is.null(args$mode)) args$mode <- mode
  if (is.null(args$part)) args$part <- part
  script <- file.path(deprecated_dense_root(), "runner", "dense_grid_monotonicity", "run_dense_grid_monotonicity.R")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, deprecated_dense_flags(args)))
  if (!identical(status, 0L)) stop("Canonical dense-grid runner failed with status ", status, call. = FALSE)
  invisible(TRUE)
}

deprecated_dense_direct <- function(mode, part) {
  raw <- commandArgs(trailingOnly = TRUE)
  if (!any(grepl("^--mode=", raw))) raw <- c(raw, paste0("--mode=", mode))
  if (!any(grepl("^--part=", raw))) raw <- c(raw, paste0("--part=", part))
  script <- file.path(deprecated_dense_root(), "runner", "dense_grid_monotonicity", "run_dense_grid_monotonicity.R")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, raw))
  if (!identical(status, 0L)) stop("Canonical dense-grid runner failed with status ", status, call. = FALSE)
  invisible(TRUE)
}

deprecated_dense_regression_outputs <- function(args = list()) {
  input_dir <- normalizePath(path.expand(args$input_dir %||% args$result_dir), mustWork = TRUE)
  out_dir <- normalizePath(path.expand(args$out_dir %||% input_dir), mustWork = FALSE)
  source_table <- file.path(input_dir, "tables", "fixed_o2_ploidy_monotonicity_curves.tsv")
  if (!file.exists(source_table)) stop("Missing pointwise curve table: ", source_table, call. = FALSE)
  simulation_dir <- file.path(out_dir, "simulation_tables")
  dir.create(simulation_dir, recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(source_table, file.path(simulation_dir, "fixed_o2_dense_grid_attractor_curves.tsv"), overwrite = TRUE)
  if (!isTRUE(copied)) stop("Could not stage the legacy pointwise curve table.", call. = FALSE)
  args$input_dir <- NULL
  args$out_dir <- out_dir
  deprecated_dense_run(args, mode = "postprocess", part = "monotonicity")
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || (length(x) == 1L && (is.na(x) || !nzchar(as.character(x))))) y else x
