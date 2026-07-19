#!/usr/bin/env Rscript

# Shared, side-effect-free HTML text escaping for report assemblers.

o2sd_html_escape_standard <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

o2sd_html_escape_full <- function(x) {
  if (is.null(x) || !length(x)) x <- ""
  x <- o2sd_html_escape_standard(x)
  gsub("'", "&#39;", x, fixed = TRUE)
}

o2sd_html_escape_minimal <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
