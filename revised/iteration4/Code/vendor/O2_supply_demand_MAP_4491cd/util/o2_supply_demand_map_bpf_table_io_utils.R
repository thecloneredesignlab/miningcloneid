#!/usr/bin/env Rscript

bpf_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.null(x)) x <- data.frame()
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

bpf_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

bpf_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = path, quote = TRUE, row.names = FALSE, na = "NA")
  invisible(path)
}

bpf_write_run_arguments <- function(args, path) {
  if (is.null(args)) args <- list()
  dat <- data.frame(
    argument = names(args),
    value = vapply(args, function(x) paste(as.character(x), collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
  bpf_write_tsv(dat, path)
}
