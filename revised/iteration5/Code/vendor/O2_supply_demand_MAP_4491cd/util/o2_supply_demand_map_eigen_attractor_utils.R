#!/usr/bin/env Rscript

# Shared, side-effect-free contracts for fixed-O2 eigen-attractor stages.
# This file contains only CLI, path, seed, and table helpers.  It does not
# source a model, run a simulation, perform an embedding, or draw a figure.

o2ea_util_script_dir <- function() {
  frames <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE)
    }, character(1L))
  )
  own <- frames[basename(frames) == "o2_supply_demand_map_eigen_attractor_utils.R"]
  if (length(own)) return(dirname(own[[length(own)]]))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  normalizePath(getwd(), mustWork = FALSE)
}

source(
  file.path(o2ea_util_script_dir(), "o2_supply_demand_map_shared.R"),
  local = TRUE,
  chdir = TRUE
)

o2ea_workflow_root <- function() {
  normalizePath(file.path(o2ea_util_script_dir(), ".."), mustWork = FALSE)
}

default_oxygen_dir <- function() {
  normalizePath(file.path(o2ea_workflow_root(), "..", ".."), mustWork = FALSE)
}

fixo2ea_default_result_root <- function() {
  file.path(
    default_oxygen_dir(),
    "results",
    "analysis",
    "best_fit_parameter_feature",
    "05_FixO2_eigen_attractor_based_clustering"
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[kv[[1L]]]] <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- tolower(trimws(as.character(x[[1L]])))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- suppressWarnings(as.integer(x[[1L]]))
  if (is.na(value)) default else value
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- suppressWarnings(as.numeric(x[[1L]]))
  if (!is.finite(value)) default else value
}

as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(default)
  values <- suppressWarnings(as.numeric(trimws(strsplit(paste(x, collapse = ","), ",", fixed = TRUE)[[1L]])))
  values <- values[is.finite(values)]
  if (length(values)) values else default
}

as_char_vec <- function(x, default = character()) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(default)
  values <- trimws(strsplit(paste(x, collapse = ","), ",", fixed = TRUE)[[1L]])
  values <- values[nzchar(values)]
  if (length(values)) values else default
}

read_tsv <- function(path) {
  utils::read.table(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
}

fixo2ea_read_tsv_plain <- read_tsv

read_csv_plain <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, quote = FALSE, row.names = FALSE)
  invisible(path)
}

rbind_fill_plain <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x), rows)
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(df) {
    for (column in setdiff(columns, names(df))) df[[column]] <- NA
    df[, columns, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

normalize_dataset <- function(dataset) {
  dataset <- tolower(trimws(as.character(dataset %||% "invivo")))
  if (!dataset %in% c("invivo", "invitro")) stop("dataset must be invivo or invitro.", call. = FALSE)
  dataset
}

seed_from_dir <- o2sd_seed_from_dir

list_seed_dirs <- function(input_dir) {
  dirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl("^seed[0-9]+$", basename(dirs))]
  dirs[order(vapply(dirs, seed_from_dir, integer(1L)))]
}

default_dataset_input_dir <- function(dataset) {
  dataset <- normalize_dataset(dataset)
  file.path(
    default_oxygen_dir(),
    "results",
    if (identical(dataset, "invivo")) "fit_invivo_O2_buffering_500seed" else "fit_invitro_O2_buffering_500seed"
  )
}

normalize_reduction <- function(reduction) {
  reduction <- tolower(trimws(as.character(reduction[[1L]])))
  aliases <- c(pca = "pca", umap = "umap", tsne = "tsne", `t-sne` = "tsne", `t_sne` = "tsne")
  if (!reduction %in% names(aliases)) stop("reduction must be pca, umap, or tsne.", call. = FALSE)
  unname(aliases[[reduction]])
}

reduction_dir_name <- function(reduction) {
  switch(normalize_reduction(reduction), pca = "PCAs", umap = "UMAPs", tsne = "TSNEs")
}

reduction_coordinate_names <- o2sd_reduction_coordinate_names
