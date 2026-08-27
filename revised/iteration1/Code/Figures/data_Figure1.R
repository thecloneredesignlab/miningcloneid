#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

data_Figure1 <- function() {
  destination_root <- file.path(DATA_ROOT, "Figure1")
  frozen_source <- LTEE_DATA_ROOT
  source_files <- sort(list.files(
    frozen_source, full.names = TRUE, recursive = FALSE
  ))
  source_files <- source_files[file.info(source_files)$isdir %in% FALSE]
  expected <- c(
    "invitro_kary_cells.tsv",
    "invitro_lineage_timeline.tsv",
    "invitro_passage_observations.tsv",
    "invivo_burden_long.tsv",
    "invivo_harvest_catalog.tsv",
    "invivo_ploidy_cells.tsv"
  )
  if (!identical(basename(source_files), sort(expected))) {
    stop("Figure 1 frozen-input set differs from the approved six-file contract.")
  }
  destinations <- file.path(destination_root, basename(source_files))
  copied <- Map(copy_input, source_files, destinations)
  contract <- data.frame(
    role = "approved Figure 1 frozen input",
    source = normalizePath(source_files, mustWork = TRUE),
    local_file = unlist(copied, use.names = FALSE),
    source_md5 = unname(tools::md5sum(source_files)),
    local_md5 = unname(tools::md5sum(unlist(copied, use.names = FALSE))),
    stringsAsFactors = FALSE
  )
  write_data_contract("Figure1", contract)
  invisible(contract)
}

if (sys.nframe() == 0L) data_Figure1()
