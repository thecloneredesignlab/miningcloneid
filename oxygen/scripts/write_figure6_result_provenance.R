#!/usr/bin/env Rscript

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
result_dir <- file.path(repo_root, "oxygen", "results", "figure6_fixed_o2_cin_map")

hpc_project <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling"
hpc_result <- file.path(hpc_project, "oxygen/results/analysis/figure6_fixed_o2_cin_map_20260721")
hpc_tables <- file.path(hpc_result, "FixO2EigenAttractors/Tables")
hpc_tasks <- file.path(hpc_result, "FixO2EigenAttractors/HPC")
hpc_joint <- file.path(hpc_result, "Joint")
local_analysis <- "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification/dense-grid_monotonicity_classification/tables"
local_parameter_landscape <- paste0(
  "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/analysis/",
  "best_fit_parameter_feature/04_combine_parameter_landscape/pooled_embedding_curve_class/",
  "tables/TSNEs/Full/",
  "pooled_invivo_invitro_initial_vs_best_context_prior_unit_tsne_best_points_curve_class.csv"
)
figure_script <- file.path(repo_root, "oxygen/figures/draw_figure6_fixed_o2_cin_map.R")
landscape_import_script <- file.path(repo_root, "oxygen/scripts/import_figure6d_parameter_landscape.R")

row <- function(destination_file, source_host, origin_path, transfer_type, purpose,
                source_git_commit, slurm_job_id = NA_character_, generating_script = NA_character_) {
  data.frame(
    destination_file = destination_file,
    source_host = source_host,
    origin_directory = dirname(origin_path),
    origin_file = basename(origin_path),
    origin_path = origin_path,
    transfer_type = transfer_type,
    purpose = purpose,
    source_git_commit = source_git_commit,
    slurm_job_id = slurm_job_id,
    generating_script = generating_script,
    stringsAsFactors = FALSE
  )
}

copied <- do.call(rbind, list(
  row("invivo_best_fixo2_eigen_attractor_wide.csv", "red.moffitt.org", file.path(hpc_tables, "invivo_best_fixo2_eigen_attractor_wide.csv"), "copied_hpc", "Figure 6B main surface", "951e90b+session_patch", "19188339;19188340"),
  row("invitro_best_fixo2_eigen_attractor_wide.csv", "red.moffitt.org", file.path(hpc_tables, "invitro_best_fixo2_eigen_attractor_wide.csv"), "copied_hpc", "Figure 6B in-vitro sensitivity", "951e90b+session_patch", "19188339;19188340"),
  row("invivo_best_fixo2_eigen_attractor_status_summary.csv", "red.moffitt.org", file.path(hpc_tables, "invivo_best_fixo2_eigen_attractor_status_summary.csv"), "copied_hpc", "Figure 6B completion QC", "951e90b+session_patch", "19188339;19188340"),
  row("invitro_best_fixo2_eigen_attractor_status_summary.csv", "red.moffitt.org", file.path(hpc_tables, "invitro_best_fixo2_eigen_attractor_status_summary.csv"), "copied_hpc", "Figure 6B sensitivity completion QC", "951e90b+session_patch", "19188339;19188340"),
  row("fixo2_eigen_feature_manifest.csv", "red.moffitt.org", file.path(hpc_tables, "fixo2_eigen_feature_manifest.csv"), "copied_hpc", "Merged-output audit", "951e90b+session_patch", "19188340"),
  row("fixo2_eigen_parameter_tasks.tsv", "red.moffitt.org", file.path(hpc_tasks, "fixo2_eigen_parameter_tasks.tsv"), "copied_hpc", "Exact 1000-point parameter input contract", "951e90b+session_patch", "19188339"),
  row("fixo2_eigen_parameter_task_manifest.tsv", "red.moffitt.org", file.path(hpc_tasks, "fixo2_eigen_parameter_task_manifest.tsv"), "copied_hpc", "Dataset and task-count audit", "951e90b+session_patch", "19188339"),
  row("joint_best_fixo2_eigen_attractor_long.tsv", "red.moffitt.org", file.path(hpc_joint, "joint_best_fixo2_eigen_attractor_long.tsv"), "copied_hpc", "Six joint best-candidate QC overlay", "951e90b+session_patch", "19189386"),
  row("joint_best_seed_summary.tsv", "red.moffitt.org", file.path(hpc_joint, "joint_best_seed_summary.tsv"), "copied_hpc", "Joint best-seed selection audit", "951e90b+session_patch", "19189386"),
  row("invivo_fixed_o2_ploidy_monotonicity_regression_class_counts.tsv", "local", file.path(local_analysis, "fixed_o2_ploidy_monotonicity_regression_class_counts.tsv"), "copied_local", "Figure 6A class counts", "f8abc65"),
  row("invivo_fixed_o2_ploidy_monotonicity_regression_curves.tsv", "local", file.path(local_analysis, "fixed_o2_ploidy_monotonicity_regression_curves.tsv"), "copied_local", "Figure 6A class curves", "f8abc65")
))

derived <- do.call(rbind, list(
  row("invivo_parameter_landscape_tsne_curve_class.csv", "local", local_parameter_landscape, "filtered_local", "Figure 6D prior-scaled parameter landscape, response class, and fit objective", "f8abc65", generating_script = landscape_import_script),
  row("invivo_fixed_o2_cin_points.tsv", "local", file.path(result_dir, "invivo_best_fixo2_eigen_attractor_wide.csv"), "derived_local", "Long-form Figure 6B point table", "7e72dca+working_tree", generating_script = figure_script),
  row("invivo_fixed_o2_cin_binned.tsv", "local", file.path(result_dir, "invivo_fixed_o2_cin_points.tsv"), "derived_local", "Figure 6B occupied-bin statistics", "7e72dca+working_tree", generating_script = figure_script),
  row("invivo_fixed_o2_cin_qc_summary.csv", "local", file.path(result_dir, "invivo_fixed_o2_cin_points.tsv"), "derived_local", "Figure 6B row, range, gap, and coverage QC", "7e72dca+working_tree", generating_script = figure_script),
  row("invitro_fixed_o2_cin_points.tsv", "local", file.path(result_dir, "invitro_best_fixo2_eigen_attractor_wide.csv"), "derived_local", "Long-form in-vitro sensitivity point table", "7e72dca+working_tree", generating_script = figure_script),
  row("invitro_fixed_o2_cin_binned.tsv", "local", file.path(result_dir, "invitro_fixed_o2_cin_points.tsv"), "derived_local", "In-vitro sensitivity occupied-bin statistics", "7e72dca+working_tree", generating_script = figure_script),
  row("invitro_fixed_o2_cin_qc_summary.csv", "local", file.path(result_dir, "invitro_fixed_o2_cin_points.tsv"), "derived_local", "In-vitro sensitivity row, range, gap, and coverage QC", "7e72dca+working_tree", generating_script = figure_script),
  row("figure6_output_manifest.csv", "local", figure_script, "generated_manifest", "Generated result and figure inventory", "7e72dca+working_tree", generating_script = figure_script)
))

provenance <- rbind(copied, derived)
actual <- basename(list.files(result_dir, full.names = TRUE))
actual <- setdiff(actual, "source_file_provenance.csv")
missing_registration <- setdiff(actual, provenance$destination_file)
missing_file <- setdiff(provenance$destination_file, actual)
if (length(missing_registration) || length(missing_file)) {
  stop(
    "Provenance registration mismatch. Unregistered: ", paste(missing_registration, collapse = ", "),
    "; missing files: ", paste(missing_file, collapse = ", "), call. = FALSE
  )
}

sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  if (!length(output)) return(NA_character_)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

destination_paths <- file.path(result_dir, provenance$destination_file)
provenance$destination_sha256 <- vapply(destination_paths, sha256, character(1L))
provenance$destination_bytes <- as.numeric(file.info(destination_paths)$size)
provenance$recorded_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
provenance <- provenance[order(provenance$destination_file), , drop = FALSE]

out_path <- file.path(result_dir, "source_file_provenance.csv")
utils::write.csv(provenance, out_path, row.names = FALSE, quote = TRUE, na = "")
message("Wrote ", out_path, " [", nrow(provenance), " files]")
