#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

bpf_parse_args <- function(args = commandArgs(trailingOnly = TRUE), normalize_keys = TRUE) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- kv[[1L]]
    if (isTRUE(normalize_keys)) key <- gsub("-", "_", key, fixed = TRUE)
    val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

bpf_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  val <- tolower(trimws(as.character(x[[1L]])))
  if (!nzchar(val)) return(default)
  if (val %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n", "off")) return(FALSE)
  default
}

bpf_as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  val <- suppressWarnings(as.integer(x[[1L]]))
  if (length(val) && is.finite(val)) val else default
}

bpf_as_num <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  val <- suppressWarnings(as.numeric(x[[1L]]))
  if (length(val) && is.finite(val)) val else default
}

bpf_split_csv <- function(x, default = character()) {
  if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(as.character(x[[1L]]))) return(default)
  vals <- trimws(strsplit(as.character(x[[1L]]), ",", fixed = TRUE)[[1L]])
  vals[nzchar(vals)]
}
