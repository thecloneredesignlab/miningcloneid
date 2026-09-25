#!/usr/bin/env Rscript

.joint_parameter_vis_dir <- local({
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

SCRIPT_DIR <- normalizePath(.joint_parameter_vis_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"),
  local = environment()
)
source(
  file.path(SCRIPT_DIR, "o2_supply_demand_map_joint_parameter_plot_utils.R"),
  local = environment()
)
rm(.joint_parameter_vis_dir)

`%||%` <- o2sd_null_coalesce

read_joint_parameter_metadata <- function(path) {
  tab <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("key", "value") %in% names(tab))) {
    stop("Invalid joint parameter metadata table: ", path, call. = FALSE)
  }
  values <- as.character(tab$value)
  names(values) <- as.character(tab$key)
  values
}

main <- function(argv = o2sd_parse_args(commandArgs(trailingOnly = TRUE))) {
  fit_dir <- argv$fit_dir %||% argv$run_dir
  if (!is.null(fit_dir)) {
    fit_dir <- normalizePath(fit_dir, mustWork = TRUE)
  }
  analysis_dir <- argv$analysis_dir %||% (
    if (!is.null(fit_dir)) {
      file.path(fit_dir, "analysis", "joint_parameters")
    } else {
      NULL
    }
  )
  out_dir <- argv$out_dir %||% argv$viz_dir %||% (
    if (!is.null(fit_dir)) file.path(fit_dir, "viz", "joint_parameters") else NULL
  )
  if (is.null(analysis_dir) || is.null(out_dir)) {
    stop(
      paste(
        "Usage: plot_joint_parameter_ratios.R",
        "--fit_dir=/abs/seed",
        "[--analysis_dir=/abs/analysis/joint_parameters]",
        "[--out_dir=/abs/viz/joint_parameters]"
      ),
      call. = FALSE
    )
  }
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)

  table_path <- file.path(analysis_dir, "joint_parameter_ratio.tsv")
  metadata_path <- file.path(
    analysis_dir,
    "joint_parameter_ratio_metadata.tsv"
  )
  if (!file.exists(table_path) || !file.exists(metadata_path)) {
    stop(
      "Missing joint parameter analysis outputs under ",
      analysis_dir,
      ". Run run_joint_parameter_diagnostics.R first.",
      call. = FALSE
    )
  }
  ratio_table <- utils::read.delim(
    table_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  metadata <- read_joint_parameter_metadata(metadata_path)
  basename_use <- metadata[["output_basename"]] %||%
    "fit_report_ratio_vivo_to_vitro_all_soft"
  fit_label <- metadata[["fit_label"]] %||% NULL

  plot <- joint_parameter_ratio_plot(ratio_table, fit_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(basename_use, ".pdf"))
  png_path <- file.path(out_dir, paste0(basename_use, ".png"))
  ggplot2::ggsave(
    pdf_path,
    plot,
    width = 18,
    height = 6.5,
    units = "in",
    bg = "white"
  )
  ggplot2::ggsave(
    png_path,
    plot,
    width = 18,
    height = 6.5,
    units = "in",
    dpi = 180,
    bg = "white"
  )
  manifest <- data.frame(
    artifact = basename(c(pdf_path, png_path)),
    format = c("pdf", "png"),
    source_analysis_table = normalizePath(table_path, mustWork = TRUE),
    bytes = as.numeric(file.info(c(pdf_path, png_path))$size),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(
    out_dir,
    "joint_parameter_visualization_manifest.tsv"
  )
  utils::write.table(
    manifest,
    manifest_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  message("Joint parameter figures written to: ", out_dir)
  invisible(c(pdf_path, png_path, manifest_path))
}

if (sys.nframe() == 0L) {
  main()
}
