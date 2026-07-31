#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

data_Figure3 <- function() {
  source_root <- file.path(INVITRO_RESULT_ROOT, "seed10")
  local_root <- file.path(DATA_ROOT, "Figure3", "source_seed10")
  files <- c(
    "invitro_flow_overlay.tsv",
    "invitro_distribution_summary.tsv",
    "invitro_daily_counts.tsv",
    "invitro_lineage_summary.tsv",
    "invitro_growth_loglik.tsv",
    "invitro_observed_kary.tsv",
    "best_params.tsv"
  )
  sources <- file.path(source_root, files)
  require_files(sources, "Figure 3 approved seed10 input")
  destinations <- file.path(local_root, files)
  copied <- Map(copy_input, sources, destinations)
  contract <- data.frame(
    role = "approved separate in-vitro seed10 post-fit input",
    source = normalizePath(sources, mustWork = TRUE),
    local_file = unlist(copied, use.names = FALSE),
    source_md5 = unname(tools::md5sum(sources)),
    local_md5 = unname(tools::md5sum(unlist(copied, use.names = FALSE))),
    stringsAsFactors = FALSE
  )
  write_data_contract("Figure3", contract)
  invisible(contract)
}

if (sys.nframe() == 0L) data_Figure3()
