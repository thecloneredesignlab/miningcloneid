#!/usr/bin/env Rscript

# Visualization-only consumer for medium-O2 analysis tables.

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE), character(1)))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  getwd()
})
source(file.path(.script_dir, "process_fingerprint_visualization_utils.R"), local = TRUE)

run_medium_o2_window_visualization <- function(argv = pfv_parse_args()) {
  simulation_dir <- pfv_as_chr(argv$simulation_dir)
  analysis_dir <- pfv_as_chr(argv$analysis_dir)
  out_dir <- pfv_as_chr(argv$viz_dir %||% argv$out_dir)
  if (!file.exists(file.path(simulation_dir, "simulation_manifest.tsv"))) stop("Missing event simulation manifest.")
  if (!file.exists(file.path(analysis_dir, "analysis_manifest.tsv"))) stop("Missing medium-O2 analysis manifest.")
  if (!nzchar(out_dir)) stop("Missing --viz_dir (or --out_dir).")
  feature_path <- file.path(analysis_dir, "tables", "medium_o2_features_by_seed.tsv")
  window_path <- file.path(analysis_dir, "tables", "medium_o2_windows.tsv")
  if (!file.exists(feature_path) || !file.exists(window_path)) stop("Missing medium-O2 visualization tables.")
  seed_features <- pfv_read_tsv(feature_path)
  windows <- pfv_read_tsv(window_path)
  cols <- paste0("mean__time_in_", windows$label)
  cols <- cols[cols %in% names(seed_features)]
  if (!length(cols) || !"trajectory_regime" %in% names(seed_features)) stop("No medium-O2 regime durations available for plotting.")
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(fig_dir, "medium_o2_duration_by_regime.pdf")
  grDevices::pdf(path, width = 8, height = 4 + 2 * length(cols))
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(oldpar); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(length(cols), 1), mar = c(8, 5, 3, 1))
  for (col in cols) {
    d <- seed_features[seed_features$trajectory_regime %in% c("mode1_ploidy_stable", "mode2_second_ploidy_collapse") & is.finite(seed_features[[col]]), , drop = FALSE]
    if (!nrow(d)) next
    graphics::boxplot(d[[col]] ~ d$trajectory_regime, las = 2, ylab = "days", main = col)
    graphics::stripchart(d[[col]] ~ d$trajectory_regime, vertical = TRUE, method = "jitter", pch = 16, cex = 0.4, col = grDevices::adjustcolor("black", alpha.f = 0.35), add = TRUE)
  }
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  manifest <- data.frame(
    artifact = "medium_o2_duration_by_regime", relative_path = file.path("figures", basename(path)),
    role = "figure", path = normalizePath(path, mustWork = TRUE), exists = TRUE,
    stringsAsFactors = FALSE
  )
  pfv_write_tsv(manifest, file.path(out_dir, "visualization_manifest.tsv"))
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_medium_o2_window_visualization()
