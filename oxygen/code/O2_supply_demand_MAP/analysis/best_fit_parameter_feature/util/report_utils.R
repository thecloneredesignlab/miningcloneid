#!/usr/bin/env Rscript

bpf_rel_path <- function(path, from_dir) {
  path <- normalizePath(path, mustWork = FALSE)
  from_dir <- normalizePath(from_dir, mustWork = FALSE)
  if (requireNamespace("tools", quietly = TRUE)) {
    return(tools::file_path_as_absolute(path))
  }
  path
}

bpf_html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}
