#!/usr/bin/env Rscript

bpf_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  }
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files)) frame_files[[length(frame_files)]] else ""
}

bpf_script_dir <- function(default = getwd()) {
  path <- bpf_script_path()
  if (nzchar(path)) dirname(path) else normalizePath(default, mustWork = FALSE)
}

bpf_repo_root <- function(start = bpf_script_dir()) {
  cur <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(12L)) {
    if (dir.exists(file.path(cur, "oxygen", "code", "O2_supply_demand_MAP"))) {
      return(normalizePath(cur, mustWork = FALSE))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not locate repository root from: ", start, call. = FALSE)
}

bpf_map_root <- function(repo_root = bpf_repo_root()) {
  file.path(repo_root, "oxygen", "code", "O2_supply_demand_MAP")
}

bpf_analysis_root <- function(repo_root = bpf_repo_root()) {
  file.path(bpf_map_root(repo_root), "analysis")
}

bpf_workflow_root <- function(repo_root = bpf_repo_root()) {
  file.path(bpf_analysis_root(repo_root), "best_fit_parameter_feature")
}

bpf_result_root <- function(repo_root = bpf_repo_root()) {
  file.path(repo_root, "oxygen", "results", "analysis", "best_fit_parameter_feature")
}

bpf_fixed_o2_code_dir <- function(repo_root = bpf_repo_root()) {
  file.path(bpf_workflow_root(repo_root), "01_fixed_o2")
}

bpf_parameter_landscape_code_dir <- function(repo_root = bpf_repo_root()) {
  file.path(bpf_workflow_root(repo_root), "02_parameter_landscape_clustering")
}

bpf_dense_grid_code_dir <- function(repo_root = bpf_repo_root()) {
  file.path(bpf_workflow_root(repo_root), "03_dense-grid_monotonicity_classification")
}

bpf_fixed_o2_result_dir <- function(repo_root = bpf_repo_root()) {
  file.path(bpf_result_root(repo_root), "01_fixed_o2", "FixO2_invivo_500seed")
}

bpf_parameter_landscape_result_dir <- function(repo_root = bpf_repo_root()) {
  file.path(bpf_result_root(repo_root), "02_parameter_landscape_clustering", "parameter_landscape")
}

bpf_dense_grid_result_root <- function(repo_root = bpf_repo_root()) {
  file.path(bpf_result_root(repo_root), "03_dense-grid_monotonicity_classification", "monotonicity_classification")
}

bpf_resolve_repo_path <- function(path, repo_root = bpf_repo_root(), mustWork = FALSE) {
  if (is.null(path) || !length(path) || is.na(path[[1]]) || !nzchar(as.character(path[[1]]))) {
    return(path)
  }
  path <- as.character(path[[1]])
  if (identical(path, "~")) {
    path <- repo_root
  } else if (startsWith(path, "~/")) {
    path <- file.path(repo_root, substring(path, 3L))
  } else if (!grepl("^/", path)) {
    path <- file.path(repo_root, path)
  }
  normalizePath(path, mustWork = mustWork)
}
