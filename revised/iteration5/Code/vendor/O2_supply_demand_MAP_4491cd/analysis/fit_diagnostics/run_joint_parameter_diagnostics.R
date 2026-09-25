#!/usr/bin/env Rscript

.joint_analysis_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  getwd()
})

SCRIPT_DIR <- normalizePath(.joint_analysis_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"),
  local = environment()
)
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_joint_parameter_utils.R"),
  local = environment()
)
rm(.joint_analysis_script_dir)

`%||%` <- o2sd_null_coalesce

write_joint_analysis_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  invisible(path)
}

main <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir %||% stop(
    paste(
      "Usage: run_joint_parameter_diagnostics.R",
      "--fit_dir=/abs/seed",
      "[--analysis_dir=/abs/output]"
    ),
    call. = FALSE
  )
  fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  analysis_dir <- argv$analysis_dir %||% argv$out_dir %||%
    file.path(fit_dir, "analysis", "joint_parameters")
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)

  ratio_table <- joint_parameter_ratio_build_data(fit_dir)
  if (!is.data.frame(ratio_table) || !nrow(ratio_table)) {
    stop(
      "No joint in-vivo/in-vitro parameter comparison could be materialized for ",
      fit_dir,
      call. = FALSE
    )
  }

  summary_map <- joint_parameter_ratio_read_summary(fit_dir)
  metadata <- data.frame(
    key = c(
      "schema_version",
      "fit_dir",
      "analysis_dir",
      "fit_label",
      "output_basename",
      "row_count",
      "generated_at"
    ),
    value = c(
      "o2sd-joint-parameter-analysis-v1",
      fit_dir,
      analysis_dir,
      joint_parameter_ratio_fit_report_label(fit_dir, summary_map),
      joint_parameter_ratio_output_basename(fit_dir),
      as.character(nrow(ratio_table)),
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    stringsAsFactors = FALSE
  )

  table_path <- file.path(analysis_dir, "joint_parameter_ratio.tsv")
  metadata_path <- file.path(
    analysis_dir,
    "joint_parameter_ratio_metadata.tsv"
  )
  write_joint_analysis_tsv(ratio_table, table_path)
  write_joint_analysis_tsv(metadata, metadata_path)

  manifest <- data.frame(
    artifact = c(basename(table_path), basename(metadata_path)),
    role = c("parameter_ratio_analysis", "analysis_metadata"),
    rows = c(nrow(ratio_table), nrow(metadata)),
    columns = c(
      paste(names(ratio_table), collapse = ","),
      paste(names(metadata), collapse = ",")
    ),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(
    analysis_dir,
    "joint_parameter_analysis_manifest.tsv"
  )
  write_joint_analysis_tsv(manifest, manifest_path)
  message("Joint parameter diagnostics written to: ", analysis_dir)
  invisible(c(table_path, metadata_path, manifest_path))
}

if (sys.nframe() == 0L) {
  main()
}
