#!/usr/bin/env Rscript

# Shared CLI, path, and table-contract helpers for dense fixed-O2 workflows.
# Numerical production, statistics, graphics, and reports live in their owning
# layers and are intentionally absent here.
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && (is.na(x) || !nzchar(trimws(as.character(x)))))) y else x
}

dense_grid_source_file <- function() {
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE), character(1)))
  own <- frames[basename(frames) == "o2_supply_demand_map_dense_grid_utils.R"]
  if (length(own)) own[[length(own)]] else ""
}

dense_grid_workflow_root <- function() {
  starts <- unique(c(dirname(dense_grid_source_file()), normalizePath(getwd(), mustWork = FALSE)))
  for (start in starts[nzchar(starts)]) {
    current <- start
    for (i in seq_len(10L)) {
      if (dir.exists(file.path(current, "simulation")) && dir.exists(file.path(current, "analysis")) && dir.exists(file.path(current, "util"))) return(normalizePath(current, mustWork = TRUE))
      candidate <- file.path(current, "oxygen", "code", "O2_supply_demand_MAP")
      if (dir.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate O2_supply_demand_MAP.", call. = FALSE)
}

dense_grid_repo_root <- function() normalizePath(file.path(dense_grid_workflow_root(), "..", "..", ".."), mustWork = FALSE)
dense_grid_default_fit_root <- function() file.path(dense_grid_repo_root(), "oxygen", "results", "fit_invivo_O2_buffering_500seed")
dense_grid_default_result_root <- function() file.path(dense_grid_repo_root(), "oxygen", "results", "analysis", "best_fit_parameter_feature", "03_dense-grid_monotonicity_classification", "monotonicity_classification")
dense_grid_normalize_part <- function(part) {
  value <- tolower(gsub("-", "_", as.character(part %||% "monotonicity"), fixed = TRUE))
  if (value %in% c("monotonicity", "dense_grid", "classification", "ploidy_monotonicity")) return("monotonicity")
  if (value %in% c("initial_ploidy", "initial", "trajectory", "initial_ploidy_trajectory")) return("initial_ploidy")
  stop("Unknown dense-grid part: ", value, call. = FALSE)
}
dense_grid_default_out_dir <- function(part = "monotonicity", result_root = dense_grid_default_result_root()) file.path(result_root, if (dense_grid_normalize_part(part) == "monotonicity") "dense-grid_monotonicity_classification" else "dense-grid_initial-ploidy_trajectory")
dense_grid_simulation_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "simulation_tables")
dense_grid_analysis_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "analysis_tables")
dense_grid_figure_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "figures")
dense_grid_report_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "report")
dense_grid_legacy_table_dir <- function(out_dir) file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), "tables")

dense_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[gsub("-", "_", kv[[1L]], fixed = TRUE)]] <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
  }
  out
}
dense_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- tolower(trimws(as.character(x[[1L]])))
  if (value %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)
  default
}
dense_as_int <- function(x, default = NA_integer_) {
  value <- suppressWarnings(as.integer(x[[1L]] %||% default))
  if (length(value) && is.finite(value)) value else default
}
dense_as_num <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(x[[1L]] %||% default))
  if (length(value) && is.finite(value)) value else default
}
dense_as_num_vec <- function(x, default = numeric()) {
  if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(trimws(as.character(x[[1L]])))) return(default)
  value <- suppressWarnings(as.numeric(trimws(unlist(strsplit(as.character(x[[1L]]), ",", fixed = TRUE)))))
  value <- value[is.finite(value)]
  if (length(value)) value else default
}
dense_as_char_vec <- function(x, default = character()) {
  if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(trimws(as.character(x[[1L]])))) return(default)
  value <- trimws(unlist(strsplit(as.character(x[[1L]]), ",", fixed = TRUE)))
  value <- value[nzchar(value)]
  if (length(value)) value else default
}
dense_seed_number <- function(seed_id) suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
dense_normalize_seed_ids <- function(seed_ids) {
  seed_ids <- as.character(seed_ids)
  ifelse(grepl("^seed", seed_ids), seed_ids, paste0("seed", seed_ids))
}
dense_num_key <- function(x) vapply(x, function(value) format(signif(as.numeric(value), 12), scientific = FALSE, trim = TRUE), character(1))
dense_default_o2_grid <- function() seq(0, 5, by = 0.05)
dense_default_times <- function() c(0, 1, 2, 5, 10, 20, 50, 100, 200, 300, 500, 700, 1000)

dense_read_tsv <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
dense_write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}
dense_rbind_fill <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x), rows)
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(data) {
    for (column in setdiff(columns, names(data))) data[[column]] <- rep(NA, nrow(data))
    data[, columns, drop = FALSE]
  })
  do.call(rbind, rows)
}
dense_require_columns <- function(data, columns, artifact = "table") {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(artifact, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(data)
}
dense_write_schema <- function(data, artifact_id, path) {
  schema <- data.frame(artifact_id = artifact_id, column_index = seq_along(data), column_name = names(data), r_type = vapply(data, function(x) paste(class(x), collapse = "/"), character(1)), required = TRUE, stringsAsFactors = FALSE)
  dense_write_tsv(schema, path)
}
dense_record_artifact <- function(out_dir, layer, artifact_id, path, data = NULL, producer, source = NA_character_) {
  manifest <- file.path(normalizePath(path.expand(out_dir), mustWork = FALSE), paste0(layer, "_artifact_manifest.tsv"))
  row <- data.frame(artifact_id = artifact_id, layer = layer, path = normalizePath(path, mustWork = FALSE), rows = if (is.data.frame(data)) nrow(data) else NA_integer_, columns = if (is.data.frame(data)) ncol(data) else NA_integer_, producer = producer, source = source, created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), stringsAsFactors = FALSE)
  old <- if (file.exists(manifest)) dense_read_tsv(manifest) else data.frame()
  if (nrow(old)) old <- old[old$artifact_id != artifact_id, , drop = FALSE]
  dense_write_tsv(dense_rbind_fill(list(old, row)), manifest)
  invisible(path)
}
