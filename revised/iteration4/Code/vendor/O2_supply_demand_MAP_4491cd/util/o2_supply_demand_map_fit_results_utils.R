# Shared, layer-neutral helpers for fit-results post-fit workflows.

o2fr_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args[startsWith(args, "--")]) {
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    out[[parts[[1]]]] <- if (length(parts) > 1L) paste(parts[-1], collapse = "=") else "TRUE"
  }
  out
}

o2fr_null_coalesce <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) y else x
}

o2fr_as_num <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(o2fr_null_coalesce(x, default)))
  if (length(value) && is.finite(value[[1]])) value[[1]] else default
}

o2fr_as_int <- function(x, default = NA_integer_) {
  value <- suppressWarnings(as.integer(o2fr_null_coalesce(x, default)))
  if (length(value) && is.finite(value[[1]])) value[[1]] else default
}

o2fr_as_chr <- function(x, default = "") {
  value <- as.character(o2fr_null_coalesce(x, default))
  if (length(value) && nzchar(value[[1]])) value[[1]] else default
}

o2fr_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  tolower(trimws(as.character(x[[1]]))) %in% c("true", "t", "1", "yes", "y", "on")
}

o2fr_first_existing_col <- function(tab, candidates) {
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit)) hit[[1]] else NA_character_
}

o2fr_read_tsv <- function(path, optional = FALSE) {
  if (!file.exists(path) || file.info(path)$size <= 1L) {
    if (optional) return(data.frame())
    stop("Missing or empty TSV: ", path, call. = FALSE)
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

o2fr_write_tsv <- function(tab, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

o2fr_require_manifests <- function(paths, stage) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop(stage, " requires completed upstream manifests: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(paths)
}

o2fr_write_manifest <- function(out_dir, stage, files, extra = list()) {
  manifest <- data.frame(stage = stage, file = basename(files), stringsAsFactors = FALSE)
  if (length(extra)) for (name in names(extra)) manifest[[name]] <- extra[[name]]
  o2fr_write_tsv(manifest, file.path(out_dir, paste0(stage, "_manifest.tsv")))
}

o2fr_run_rscript_stage <- function(script, args, label) {
  out <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(args)), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (length(out)) message(paste(out, collapse = "\n"))
  if (as.integer(status) != 0L) stop(label, " failed with status ", status, call. = FALSE)
  invisible(out)
}
