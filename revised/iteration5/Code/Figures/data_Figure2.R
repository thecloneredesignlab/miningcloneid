#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

data_Figure2 <- function() {
  contract <- data.frame(
    role = "native vector drawing",
    source = NA_character_,
    local_file = normalizePath(
      file.path(CODE_ROOT, "draw_Figure2.R"), mustWork = TRUE
    ),
    source_md5 = NA_character_,
    local_md5 = unname(tools::md5sum(
      normalizePath(
        file.path(CODE_ROOT, "draw_Figure2.R"), mustWork = TRUE
      )
    )),
    stringsAsFactors = FALSE
  )
  write_data_contract("Figure2", contract)
  invisible(contract)
}

if (sys.nframe() == 0L) data_Figure2()
