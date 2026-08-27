#!/usr/bin/env Rscript

# Shared CLI, path, and tabular I/O contracts for combined embedding workflows.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && (is.na(x) || !nzchar(trimws(as.character(x)))))) y else x
}

.combined_utils_path <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else normalizePath(x$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "o2_supply_demand_map_combined_landscape_utils.R"]
  if (length(own)) own[[length(own)]] else ""
})
combined_utils_file <- function() .combined_utils_path

combined_workflow_root <- function() normalizePath(file.path(dirname(combined_utils_file()), ".."), mustWork = TRUE)
combined_repo_root <- function() normalizePath(file.path(combined_workflow_root(), "..", "..", ".."), mustWork = TRUE)
combined_results_root <- function() file.path(combined_repo_root(), "oxygen", "results", "analysis", "best_fit_parameter_feature")
combined_default_pooled_root <- function() file.path(combined_results_root(), "02_parameter_landscape_clustering", "pooled_invivo_invitro")
combined_default_dense_root <- function() file.path(combined_results_root(), "03_dense-grid_monotonicity_classification", "monotonicity_classification", "dense-grid_monotonicity_classification")
combined_default_out_dir <- function() file.path(combined_results_root(), "04_combine_parameter_landscape", "pooled_embedding_curve_class")
combined_analysis_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "analysis_tables")
combined_figure_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "figures")
combined_report_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "report")
combined_legacy_table_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "tables")

combined_parse_args <- function(raw = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in raw[startsWith(raw, "--")]) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[gsub("-", "_", pair[[1L]], fixed = TRUE)]] <- if (length(pair) > 1L) paste(pair[-1L], collapse = "=") else "TRUE"
  }
  out
}
combined_as_bool <- function(x, default = FALSE) {
  value <- tolower(trimws(as.character(x[[1L]] %||% default)))
  if (value %in% c("true", "t", "1", "yes", "y")) TRUE else if (value %in% c("false", "f", "0", "no", "n")) FALSE else default
}
combined_split <- function(x, default) {
  if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(x[[1L]])) return(default)
  trimws(unlist(strsplit(as.character(x[[1L]]), ",", fixed = TRUE)))
}
combined_normalize_reductions <- function(x) {
  values <- tolower(gsub("-", "_", combined_split(x, c("PCAs", "UMAPs", "TSNEs")), fixed = TRUE))
  map <- c(pca = "PCAs", pcas = "PCAs", umap = "UMAPs", umaps = "UMAPs", tsne = "TSNEs", tsnes = "TSNEs", t_sne = "TSNEs", t_snes = "TSNEs")
  out <- unname(map[values])
  if (anyNA(out)) stop("Unknown reduction(s): ", paste(values[is.na(out)], collapse = ", "), call. = FALSE)
  unique(out)
}
combined_normalize_variants <- function(x) {
  values <- tolower(combined_split(x, c("Full", "Sampled")))
  map <- c(full = "Full", all_points = "Full", sampled = "Sampled", sample = "Sampled", sampled500 = "Sampled")
  out <- unname(map[values])
  if (anyNA(out)) stop("Unknown variant(s): ", paste(values[is.na(out)], collapse = ", "), call. = FALSE)
  unique(out)
}

combined_read_tsv <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
combined_read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
combined_write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}
combined_write_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data, path, row.names = FALSE, na = "NA")
  invisible(path)
}
combined_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(data)
}
combined_seed_number <- function(x) suppressWarnings(as.integer(sub("^.*?([0-9]+)$", "\\1", as.character(x))))
combined_find_latest <- function(root, basename_pattern, explicit = NULL) {
  if (!is.null(explicit) && length(explicit) && !is.na(explicit[[1L]]) && nzchar(explicit[[1L]])) return(normalizePath(path.expand(explicit[[1L]]), mustWork = TRUE))
  candidates <- list.files(root, pattern = basename_pattern, recursive = TRUE, full.names = TRUE)
  if (!length(candidates)) stop("No input matching ", basename_pattern, " under ", root, call. = FALSE)
  candidates[[which.max(file.info(candidates)$mtime)]]
}

combined_discover_coordinates <- function(pooled_root, reductions, variants) {
  rows <- list()
  for (reduction in reductions) for (variant in variants) {
    table_dir <- file.path(pooled_root, reduction, "Tables", variant)
    if (!dir.exists(table_dir)) next
    files <- sort(list.files(table_dir, pattern = "_coordinates[.]csv$", full.names = TRUE))
    for (path in files) {
      stub <- sub("_coordinates[.]csv$", "", basename(path))
      original <- sub("_coordinates[.]csv$", ".png", gsub("/Tables/", "/Figures/", normalizePath(path, mustWork = FALSE), fixed = TRUE))
      rows[[length(rows) + 1L]] <- data.frame(reduction = reduction, variant = variant, stub = stub, coordinate_table = normalizePath(path, mustWork = FALSE), original_png = original, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(reduction = character(), variant = character(), stub = character(), coordinate_table = character(), original_png = character(), stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

combined_coordinate_columns <- function(reduction, names) {
  candidates <- switch(reduction, PCAs = list(c("PCA1", "PCA2")), UMAPs = list(c("UMAP1", "UMAP2")), TSNEs = list(c("tSNE1", "tSNE2"), c("TSNE1", "TSNE2")))
  for (candidate in candidates) if (all(candidate %in% names)) return(candidate)
  stop("Missing ", reduction, " coordinate columns.", call. = FALSE)
}
