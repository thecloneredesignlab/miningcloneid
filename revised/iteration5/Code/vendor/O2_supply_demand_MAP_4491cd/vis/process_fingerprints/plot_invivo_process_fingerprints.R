#!/usr/bin/env Rscript

# Pure visualization consumer for process-fingerprint simulation and analysis
# tables.  It never reads fit artifacts or starts simulation/analysis.

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE), character(1)))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  getwd()
})
WORKFLOW_ROOT <- normalizePath(file.path(.script_dir, "..", ".."), mustWork = FALSE)
source(file.path(.script_dir, "process_fingerprint_visualization_utils.R"), local = TRUE)

pfv_require <- function(base, manifest_name, files, stage) {
  mp <- file.path(base, manifest_name)
  if (!file.exists(mp)) stop(stage, " missing manifest: ", mp)
  manifest <- pfv_read_tsv(mp)
  missing <- files[!file.exists(file.path(base, files))]
  if (length(missing)) stop(stage, " missing required upstream tables: ", paste(missing, collapse = ", "))
  invisible(manifest)
}

pfv_read_matrix <- function(path) {
  tab <- pfv_read_tsv(path)
  if (!nrow(tab) || ncol(tab) < 2L) return(NULL)
  ids <- as.character(tab[[1]])
  mat <- as.matrix(data.frame(lapply(tab[-1], as.numeric), check.names = FALSE))
  rownames(mat) <- ids
  colnames(mat) <- names(tab)[-1]
  mat
}

pfv_recommended_k <- function(diagnostics) {
  eligible <- diagnostics[diagnostics$min_cluster_size >= 2L & is.finite(diagnostics$mean_silhouette), , drop = FALSE]
  if (!nrow(eligible)) return(NA_integer_)
  as.integer(eligible$k[[which.max(eligible$mean_silhouette)]])
}

run_invivo_process_fingerprint_visualization <- function(argv = pfv_parse_args()) {
  simulation_dir <- pfv_as_chr(argv$simulation_dir)
  analysis_dir <- pfv_as_chr(argv$analysis_dir)
  out_dir <- pfv_as_chr(argv$viz_dir %||% argv$out_dir)
  if (!dir.exists(simulation_dir)) stop("Missing --simulation_dir.")
  if (!dir.exists(analysis_dir)) stop("Missing --analysis_dir.")
  if (!nzchar(out_dir)) stop("Missing --viz_dir (or --out_dir).")
  simulation_dir <- normalizePath(simulation_dir, mustWork = TRUE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)

  sim_files <- file.path("tables", c("process_fingerprint_static_full_long.tsv"))
  ana_files <- file.path("tables", c(
    "seed_manifest.tsv", "process_fingerprint_static_scaled.tsv",
    "process_fingerprint_dynamic_scaled.tsv", "distance_parameter_matrix.tsv",
    "distance_static_full_matrix.tsv", "distance_static_18only_matrix.tsv",
    "distance_dynamic_matrix.tsv", "cluster_k_diagnostics.tsv",
    "cluster_membership_static_full.tsv", "cluster_medoids.tsv",
    "cluster_consensus_matrix.tsv"
  ))
  pfv_require(simulation_dir, "simulation_manifest.tsv", sim_files, "Process visualization")
  pfv_require(analysis_dir, "analysis_manifest.tsv", ana_files, "Process visualization")
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  read_a <- function(name) pfv_read_tsv(file.path(analysis_dir, "tables", name))
  diagnostics <- read_a("cluster_k_diagnostics.tsv")
  static_diag <- diagnostics[diagnostics$cluster_source == "static_full", , drop = FALSE]
  membership <- read_a("cluster_membership_static_full.tsv")
  medoids <- read_a("cluster_medoids.tsv")
  medoids <- medoids[medoids$cluster_source == "static_full", , drop = FALSE]
  consensus <- pfv_read_matrix(file.path(analysis_dir, "tables", "cluster_consensus_matrix.tsv"))
  cluster <- list(
    membership = membership,
    medoids = medoids,
    consensus = consensus,
    recommended_k = pfv_recommended_k(static_diag)
  )
  pfv_plot_process_outputs(
    out_dir,
    pfv_read_matrix(file.path(analysis_dir, "tables", "distance_parameter_matrix.tsv")),
    pfv_read_matrix(file.path(analysis_dir, "tables", "distance_static_full_matrix.tsv")),
    pfv_read_matrix(file.path(analysis_dir, "tables", "distance_static_18only_matrix.tsv")),
    pfv_read_matrix(file.path(analysis_dir, "tables", "distance_dynamic_matrix.tsv")),
    read_a("process_fingerprint_static_scaled.tsv"),
    read_a("process_fingerprint_dynamic_scaled.tsv"),
    cluster,
    read_a("seed_manifest.tsv"),
    pfv_read_tsv(file.path(simulation_dir, "tables", "process_fingerprint_static_full_long.tsv"))
  )
  figures <- list.files(file.path(out_dir, "figures"), pattern = "[.]pdf$", full.names = FALSE)
  if (!length(figures)) stop("Process visualization produced no figures.")
  manifest <- data.frame(
    artifact = tools::file_path_sans_ext(figures),
    relative_path = file.path("figures", figures), role = "figure",
    path = normalizePath(file.path(out_dir, "figures", figures), mustWork = TRUE),
    exists = TRUE, stringsAsFactors = FALSE
  )
  pfv_write_tsv(manifest, file.path(out_dir, "visualization_manifest.tsv"))
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_invivo_process_fingerprint_visualization()
