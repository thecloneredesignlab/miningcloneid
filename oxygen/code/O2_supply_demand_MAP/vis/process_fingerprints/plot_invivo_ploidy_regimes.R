#!/usr/bin/env Rscript

# Pure visualization consumer for ploidy-regime tables.

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE), character(1)))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(.script_dir, "process_fingerprint_visualization_utils.R"), local = TRUE)

run_invivo_ploidy_regime_visualization <- function(argv = pfv_parse_args()) {
  simulation_dir <- pfv_as_chr(argv$simulation_dir)
  analysis_dir <- pfv_as_chr(argv$analysis_dir)
  out_dir <- pfv_as_chr(argv$viz_dir %||% argv$out_dir)
  if (!file.exists(file.path(simulation_dir, "simulation_manifest.tsv"))) stop("Missing process simulation manifest.")
  if (!file.exists(file.path(analysis_dir, "analysis_manifest.tsv"))) stop("Missing ploidy-regime analysis manifest.")
  if (!nzchar(out_dir)) stop("Missing --viz_dir (or --out_dir).")
  required <- c(
    "trajectory_curves_common_input.tsv", "trajectory_regime_labels.tsv",
    "process_combined_full_scaled.tsv"
  )
  missing <- required[!file.exists(file.path(analysis_dir, "tables", required))]
  if (length(missing)) stop("Ploidy-regime visualization missing analysis tables: ", paste(missing, collapse = ", "))
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  read_a <- function(name) pfv_read_tsv(file.path(analysis_dir, "tables", name))
  traj <- list(curves = read_a("trajectory_curves_common_input.tsv"))
  scaled <- list(combined_full = list(wide = read_a("process_combined_full_scaled.tsv")))
  labels <- read_a("trajectory_regime_labels.tsv")
  pfv_plot_ploidy_regime_outputs(traj = traj, scaled = scaled, labels = labels, out_dir = out_dir)
  figures <- list.files(file.path(out_dir, "figures"), pattern = "[.]pdf$", full.names = FALSE)
  if (!length(figures)) stop("Ploidy-regime visualization produced no figures.")
  manifest <- data.frame(
    artifact = tools::file_path_sans_ext(figures),
    relative_path = file.path("figures", figures), role = "figure",
    path = normalizePath(file.path(out_dir, "figures", figures), mustWork = TRUE),
    exists = TRUE, stringsAsFactors = FALSE
  )
  pfv_write_tsv(manifest, file.path(out_dir, "visualization_manifest.tsv"))
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_invivo_ploidy_regime_visualization()
