#!/usr/bin/env Rscript

# Shared, source-only table and CLI normalization helpers for the combined
# FixO2/eigen-attractor analysis and visualization stages.

o2cfe_read_csv_plain <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

o2cfe_read_tsv_plain <- function(path) {
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

o2cfe_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

o2cfe_split_csv <- function(x, default = character()) {
  if (is.null(x) || !length(x) || is.na(x[[1L]]) ||
      !nzchar(as.character(x[[1L]]))) {
    return(default)
  }
  values <- trimws(strsplit(as.character(x[[1L]]), ",", fixed = TRUE)[[1L]])
  values[nzchar(values)]
}

o2cfe_normalize_reductions <- function(x) {
  values <- tolower(gsub(
    "-",
    "_",
    o2cfe_split_csv(x, c("PCAs", "UMAPs", "TSNEs"))
  ))
  out <- character()
  for (value in values) {
    if (value %in% c("pca", "pcas")) {
      out <- c(out, "PCAs")
    } else if (value %in% c("umap", "umaps")) {
      out <- c(out, "UMAPs")
    } else if (value %in% c("tsne", "tsnes", "t_sne", "t_snes")) {
      out <- c(out, "TSNEs")
    } else if (nzchar(value)) {
      stop("Unknown reduction: ", value, call. = FALSE)
    }
  }
  unique(out)
}

o2cfe_normalize_variants <- function(x) {
  values <- tolower(o2cfe_split_csv(x, c("Full", "BestOnly")))
  out <- character()
  for (value in values) {
    if (value %in% c("full", "all_points")) {
      out <- c(out, "Full")
    } else if (value %in% c("bestonly", "best_only", "best")) {
      out <- c(out, "BestOnly")
    } else if (value %in% c("sampled", "sample", "sampled500")) {
      out <- c(out, "Sampled")
    } else if (nzchar(value)) {
      stop("Unknown variant: ", value, call. = FALSE)
    }
  }
  unique(out)
}
