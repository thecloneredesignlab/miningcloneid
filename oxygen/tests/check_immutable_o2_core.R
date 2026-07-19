#!/usr/bin/env Rscript

resolve_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

sha256_file <- function(path) {
  shasum <- Sys.which("shasum")
  if (!nzchar(shasum)) stop("Required command 'shasum' was not found.")
  out <- system2(shasum, c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("shasum failed for ", path, ": ", paste(out, collapse = "\n"))
  }
  sub("[[:space:]].*$", "", out[[1]])
}

script_dir <- resolve_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
workflow_root <- file.path(repo_root, "oxygen", "code", "O2_supply_demand_MAP")
full_mode <- "--full" %in% commandArgs(trailingOnly = TRUE)

if (full_mode) {
  manifest_path <- file.path(
    workflow_root,
    "docs",
    "immutable_core_full_manifest.tsv"
  )
  manifest <- utils::read.delim(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c("path", "size", "mtime", "sha256")
  if (!all(required_columns %in% names(manifest))) {
    stop("Full immutable-core manifest has invalid columns: ", manifest_path)
  }
  protected_roots <- file.path(workflow_root, c("model", "optimizer"))
  current_files <- sort(unlist(lapply(
    protected_roots,
    list.files,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  ), use.names = FALSE))
  current_files <- current_files[file.exists(current_files) & !dir.exists(current_files)]
  expected_files <- normalizePath(
    file.path(repo_root, manifest$path),
    mustWork = FALSE
  )
  if (!identical(current_files, sort(expected_files))) {
    stop(
      "Protected-tree file set changed. Missing: ",
      paste(setdiff(expected_files, current_files), collapse = ", "),
      "; added: ",
      paste(setdiff(current_files, expected_files), collapse = ", ")
    )
  }
  paths <- expected_files
} else {
  manifest_path <- file.path(workflow_root, "docs", "immutable_core_sha256.tsv")
  manifest <- utils::read.delim(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c("path", "sha256")
  if (!all(required_columns %in% names(manifest))) {
    stop("Tracked immutable-core manifest has invalid columns: ", manifest_path)
  }
  paths <- normalizePath(
    file.path(workflow_root, manifest$path),
    mustWork = FALSE
  )
}

missing <- paths[!file.exists(paths)]
if (length(missing)) {
  stop("Missing protected file(s): ", paste(missing, collapse = ", "))
}

actual_sha <- vapply(paths, sha256_file, character(1))
sha_ok <- all(tolower(unname(actual_sha)) == tolower(manifest$sha256))
if (!sha_ok) {
  bad <- which(tolower(actual_sha) != tolower(manifest$sha256))
  stop(
    "Protected SHA-256 mismatch: ",
    paste(manifest$path[bad], collapse = ", ")
  )
}

if (full_mode) {
  info <- file.info(paths)
  actual_size <- as.numeric(info$size)
  actual_mtime <- format(
    info$mtime,
    format = "%Y-%m-%dT%H:%M:%OS6%z",
    tz = "America/New_York"
  )
  size_ok <- all(actual_size == as.numeric(manifest$size))
  mtime_ok <- all(actual_mtime == manifest$mtime)
  if (!size_ok || !mtime_ok) {
    bad <- which(
      actual_size != as.numeric(manifest$size) |
        actual_mtime != manifest$mtime
    )
    stop(
      "Protected size/mtime mismatch: ",
      paste(manifest$path[bad], collapse = ", ")
    )
  }
}

cat(
  "IMMUTABLE_CORE_CHECK=PASS\n",
  "mode=", if (full_mode) "full" else "tracked", "\n",
  "files=", length(paths), "\n",
  "manifest=", normalizePath(manifest_path, mustWork = TRUE), "\n",
  sep = ""
)
