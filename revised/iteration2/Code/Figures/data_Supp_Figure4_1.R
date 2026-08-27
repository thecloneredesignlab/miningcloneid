#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

data_Supp_Figure4_1 <- function() {
  figure4_dir <- file.path(DATA_ROOT, "Figure4")
  si_dir <- file.path(DATA_ROOT, "Supp_Figure4_1")
  required <- file.path(figure4_dir, c(
    "fixed_o2_dominant_ploidy_201grid.tsv",
    "invivo_best_parameters_500seeds.tsv",
    "parameter_function_groups.tsv",
    "parameter_function_group_palette.tsv",
    "invivo_parameter_table_seed25.csv",
    "invivo_best_tsne_cluster_coordinates.tsv",
    "invivo_tsne_cluster_summary.tsv",
    "invivo_tsne_run_metadata.tsv",
    "figure4a_seed25_burden_timecourse.tsv",
    "invivo_fit_objective_ranking_500seeds.tsv"
  ))
  require_files(required, "Supplementary Figure 4-1 shared Figure 4 input")
  local_files <- file.path(si_dir, basename(required))
  for (i in seq_along(required)) {
    ok <- file.copy(required[[i]], local_files[[i]], overwrite = TRUE)
    if (!ok) stop("Failed to stage Supplementary Figure 4-1 intermediate: ", required[[i]])
  }
  contract <- data.frame(
    role = "shared Figure 4B derived input",
    source = normalizePath(required, mustWork = TRUE),
    local_file = normalizePath(local_files, mustWork = TRUE),
    source_md5 = unname(tools::md5sum(required)),
    local_md5 = unname(tools::md5sum(local_files)),
    stringsAsFactors = FALSE
  )
  write_data_contract("Supp_Figure4_1", contract)
  invisible(contract)
}

if (sys.nframe() == 0L) data_Supp_Figure4_1()
