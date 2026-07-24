#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) {
  stop("MININGCLONEID_REPO_ROOT is not set. Run with scripts/agentRrunner.sh.")
}
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
draft_rel <- "agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting"
draft_root <- file.path(repo_root, draft_rel)
frozen_rel_root <- file.path(draft_rel, "source_tables", "frozen_inputs")
frozen_root <- file.path(repo_root, frozen_rel_root)
bundle_rel <- "oxygen/results/fitting_output_bundle_20260722"
bundle_root <- file.path(repo_root, bundle_rel)

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, ": ", paste(output, collapse = "\n"))
  }
  sub("[[:space:]].*$", "", output[[1]])
}

records <- list()
record_index <- 0L

copy_input <- function(figure, original_rel, frozen_rel, role = "frozen_generator_input") {
  original_abs <- file.path(repo_root, original_rel)
  frozen_abs <- file.path(repo_root, frozen_rel)
  if (!file.exists(original_abs)) stop("Missing upstream input: ", original_rel)
  dir.create(dirname(frozen_abs), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(original_abs, frozen_abs, overwrite = TRUE, copy.mode = FALSE)) {
    stop("Could not copy ", original_rel, " to ", frozen_rel)
  }
  original_sha <- sha256_file(original_abs)
  frozen_sha <- sha256_file(frozen_abs)
  if (!identical(original_sha, frozen_sha)) {
    stop("Frozen copy checksum mismatch: ", frozen_rel)
  }
  record_index <<- record_index + 1L
  records[[record_index]] <<- data.frame(
    figure = figure,
    role = role,
    original_path = gsub("\\\\", "/", original_rel),
    frozen_path = gsub("\\\\", "/", frozen_rel),
    sha256 = frozen_sha,
    bytes = as.numeric(file.info(frozen_abs)$size),
    stringsAsFactors = FALSE
  )
}

copy_many <- function(figure, originals, frozen_dir, role = "frozen_generator_input") {
  for (original in originals) {
    copy_input(
      figure,
      original,
      file.path(frozen_dir, basename(original)),
      role
    )
  }
}

# Figure 1 observed tables.
figure1_sources <- file.path(
  "figures",
  c(
    "invitro_lineage_timeline.tsv",
    "invitro_passage_observations.tsv",
    "invitro_kary_cells.tsv",
    "invivo_burden_long.tsv",
    "invivo_ploidy_cells.tsv",
    "invivo_harvest_catalog.tsv"
  )
)
copy_many("Figure 1", figure1_sources, file.path(frozen_rel_root, "F1"))

# Figure 3 selected saved tables, preserving the viz/ subdirectory.
figure3_base <- file.path(
  bundle_rel, "separate_fits/invitro/selected_seeds/seed10"
)
figure3_relpaths <- c(
  "invitro_lineage_summary.tsv",
  "invitro_distribution_summary.tsv",
  "invitro_observed_kary.tsv",
  "viz/functional_curve_ploidy.tsv",
  "viz/functional_curve_oxygen_multi_ploidy.tsv"
)
for (relative_path in figure3_relpaths) {
  copy_input(
    "Figure 3",
    file.path(figure3_base, relative_path),
    file.path(frozen_rel_root, "F3", relative_path)
  )
}

# Figure 4 selected fit, fixed-O2 analysis, pooled coordinates, and pinned raster.
figure4_sources <- c(
  burden_fit = file.path(
    bundle_rel, "separate_fits/invivo/selected_seeds/seed25/burden_fit.tsv"
  ),
  terminal_ploidy = file.path(
    bundle_rel,
    "separate_fits/invivo/selected_seeds/seed25/viz/terminal_ploidy_observed_vs_predicted.tsv"
  ),
  o2_lag = file.path(
    bundle_rel,
    "separate_fits/invivo/selected_seeds/seed25/viz/o2_lag_timecourse.tsv"
  ),
  fixed_o2_modes = file.path(
    bundle_rel,
    "supporting_analysis/fixed_o2_analysis_tables/analytical_attractors/fixed_o2_attractor_mode_by_seed_o2.tsv"
  ),
  parameter_values = file.path(
    bundle_rel,
    "supporting_analysis/fixed_o2_analysis_tables/analytical_attractors/parameter_values_long.tsv"
  ),
  pooled_coordinates = file.path(
    bundle_rel,
    paste0(
      "supporting_analysis/pooled_embedding_curve_class/tables/TSNEs/Full/",
      "pooled_invivo_invitro_initial_vs_best_tsne_best_points_curve_class.csv"
    )
  ),
  seed_summary = file.path(
    bundle_rel, "separate_fits/invivo/run/extra_results/seed_summary.tsv"
  )
)
for (source_name in names(figure4_sources)) {
  extension <- tools::file_ext(figure4_sources[[source_name]])
  copy_input(
    "Figure 4",
    figure4_sources[[source_name]],
    file.path(frozen_rel_root, "F4", paste0(source_name, ".", extension))
  )
}

pinned_blob <- paste0(
  "7e72dca88caf784fc61271d87a1c0dfb564b8303:",
  "oxygen/figures/iteration1/fig4f_landscape_tsne_clusters.png"
)
pinned_output_rel <- file.path(
  frozen_rel_root, "F4", "pooled_embedding_preserved_source.png"
)
pinned_output_abs <- file.path(repo_root, pinned_output_rel)
dir.create(dirname(pinned_output_abs), recursive = TRUE, showWarnings = FALSE)
blob_size <- as.numeric(system2(
  "git", c("cat-file", "-s", pinned_blob),
  stdout = TRUE, stderr = TRUE
))
if (length(blob_size) != 1L || !is.finite(blob_size) || blob_size <= 0) {
  stop("Could not resolve pinned Figure 4 raster.")
}
blob_connection <- pipe(
  paste("git cat-file blob", shQuote(pinned_blob)),
  open = "rb"
)
blob_raw <- readBin(blob_connection, what = "raw", n = as.integer(blob_size))
close(blob_connection)
if (length(blob_raw) != as.integer(blob_size)) stop("Pinned raster read was incomplete.")
output_connection <- file(pinned_output_abs, open = "wb")
writeBin(blob_raw, output_connection)
close(output_connection)
if (unname(tools::md5sum(pinned_output_abs)) != "14cecff29ab4884823b84d83f0379119") {
  stop("Pinned Figure 4 raster MD5 mismatch.")
}
record_index <- record_index + 1L
records[[record_index]] <- data.frame(
  figure = "Figure 4",
  role = "frozen_pinned_git_blob",
  original_path = paste0("git:", pinned_blob),
  frozen_path = gsub("\\\\", "/", pinned_output_rel),
  sha256 = sha256_file(pinned_output_abs),
  bytes = as.numeric(file.info(pinned_output_abs)$size),
  stringsAsFactors = FALSE
)

# Figure 5 exact six-winner selection and all tables read by the generator.
selection_original_rel <- file.path(bundle_rel, "selected_results.tsv")
selection_frozen_rel <- file.path(frozen_rel_root, "F5", "selected_results.tsv")
copy_input("Figure 5", selection_original_rel, selection_frozen_rel, "frozen_selection_table")
selection <- read.delim(
  file.path(repo_root, selection_original_rel),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
selection <- selection[selection$record_type == "joint_pair_best", , drop = FALSE]
if (nrow(selection) != 6L) stop("Expected exactly six joint_pair_best records.")
figure5_relpaths <- c(
  "invivo_burden_fit.tsv",
  "invivo_terminal_ploidy_fit.tsv",
  "invitro_growth_loglik.tsv",
  "invitro_lineage_summary.tsv",
  "joint_soft_coupling.tsv",
  "viz/invivo/functional_curve_oxygen_multi_ploidy.tsv",
  "viz/invitro/functional_curve_oxygen_multi_ploidy.tsv",
  "viz/invivo/functional_curve_ploidy.tsv",
  "viz/invitro/functional_curve_ploidy.tsv"
)
for (i in seq_len(nrow(selection))) {
  for (relative_path in figure5_relpaths) {
    copy_input(
      "Figure 5",
      file.path(bundle_rel, selection$bundle_dir[[i]], relative_path),
      file.path(
        frozen_rel_root, "F5", "winners",
        selection$warmup_label[[i]], relative_path
      ),
      paste0("frozen winner input: ", selection$warmup_label[[i]])
    )
  }
}

# Supplementary diagnostics: exact saved objective and run-status tables.
diagnostic_files <- c(
  "seed_objective_simple.tsv", "seed_summary.tsv", "convergence_summary.tsv"
)
for (context in c("invitro", "invivo")) {
  for (filename in diagnostic_files) {
    copy_input(
      "Supplementary diagnostics",
      file.path(
        bundle_rel, "separate_fits", context, "run/extra_results", filename
      ),
      file.path(
        frozen_rel_root, "diagnostics", "separate", context, filename
      )
    )
  }
}
for (i in seq_len(nrow(selection))) {
  for (filename in diagnostic_files) {
    copy_input(
      "Supplementary diagnostics",
      file.path(
        bundle_rel, "joint_fit/pairs", selection$warmup_label[[i]],
        "run/extra_results", filename
      ),
      file.path(
        frozen_rel_root, "diagnostics", "joint",
        selection$warmup_label[[i]], filename
      ),
      paste0("frozen optimizer input: ", selection$warmup_label[[i]])
    )
  }
}

manifest <- do.call(rbind, records)
manifest <- manifest[order(manifest$figure, manifest$frozen_path), , drop = FALSE]
manifest_path <- file.path(draft_root, "source_tables", "frozen_input_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE, quote = TRUE)
message(
  "Materialized ", nrow(manifest), " frozen inputs (",
  sprintf("%.1f", sum(manifest$bytes) / 1024^2), " MB)."
)
