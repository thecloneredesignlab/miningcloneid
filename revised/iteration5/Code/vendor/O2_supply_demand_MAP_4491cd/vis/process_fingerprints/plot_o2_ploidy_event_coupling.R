#!/usr/bin/env Rscript

# Visualization-only consumer for materialized O2/ploidy event tables.

.script_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) if (is.null(env$ofile)) "" else normalizePath(env$ofile, mustWork = FALSE), character(1)))
  if (length(frame_files)) return(dirname(frame_files[[length(frame_files)]]))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  getwd()
})
source(file.path(.script_dir, "process_fingerprint_visualization_utils.R"), local = TRUE)

run_o2_ploidy_event_visualization <- function(argv = pfv_parse_args()) {
  simulation_dir <- pfv_as_chr(argv$simulation_dir)
  analysis_dir <- pfv_as_chr(argv$analysis_dir)
  out_dir <- pfv_as_chr(argv$viz_dir %||% argv$out_dir)
  if (!file.exists(file.path(simulation_dir, "simulation_manifest.tsv"))) stop("Missing event simulation manifest.")
  if (!file.exists(file.path(analysis_dir, "analysis_manifest.tsv"))) stop("Missing event analysis manifest.")
  if (!nzchar(out_dir)) stop("Missing --viz_dir (or --out_dir).")
  required_sim <- c("event_ploidy_timecourses.tsv", "event_o2_timecourses.tsv")
  required_analysis <- c("seed_event_summary.tsv", "cohort_event_coupling.tsv")
  missing <- c(
    file.path(simulation_dir, "tables", required_sim)[!file.exists(file.path(simulation_dir, "tables", required_sim))],
    file.path(analysis_dir, "tables", required_analysis)[!file.exists(file.path(analysis_dir, "tables", required_analysis))]
  )
  if (length(missing)) stop("Missing event visualization inputs: ", paste(missing, collapse = ", "))
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  ploidy <- pfv_read_tsv(file.path(simulation_dir, "tables", "event_ploidy_timecourses.tsv"))
  if ("x" %in% names(ploidy) && !"ploidy" %in% names(ploidy)) names(ploidy)[names(ploidy) == "x"] <- "ploidy"
  oxygen <- pfv_read_tsv(file.path(simulation_dir, "tables", "event_o2_timecourses.tsv"))
  summary <- pfv_read_tsv(file.path(analysis_dir, "tables", "seed_event_summary.tsv"))
  events <- pfv_read_tsv(file.path(analysis_dir, "tables", "cohort_event_coupling.tsv"))
  cols <- c(mode1_ploidy_stable = "#1b9e77", mode2_second_ploidy_collapse = "#d95f02", ambiguous = "grey60")
  color_for <- function(label, alpha = 1) {
    col <- unname(cols[as.character(label)])
    if (!length(col) || is.na(col)) col <- "grey60"
    grDevices::adjustcolor(col, alpha.f = alpha)
  }

  if (nrow(ploidy)) {
    d <- merge(ploidy, summary[, c("seed_id", "trajectory_regime")], by = "seed_id", all.x = TRUE)
    grDevices::pdf(file.path(fig_dir, "ploidy_trajectories_by_regime.pdf"), width = 9, height = 6)
    graphics::plot(NA, xlim = range(d$day, na.rm = TRUE), ylim = range(d$ploidy, na.rm = TRUE), xlab = "Day", ylab = "Ploidy", main = "Ploidy trajectories by inferred regime")
    for (key in split(seq_len(nrow(d)), paste(d$seed_id, d$cohort, sep = "||"))) {
      x <- d[key, , drop = FALSE]
      graphics::lines(x$day, x$ploidy, col = color_for(x$trajectory_regime[[1]], 0.25), lwd = 1)
    }
    graphics::legend("topright", legend = names(cols), col = cols, lwd = 2, bty = "n")
    grDevices::dev.off()
  }
  if (nrow(oxygen)) {
    d <- merge(oxygen, summary[, c("seed_id", "trajectory_regime")], by = "seed_id", all.x = TRUE)
    grDevices::pdf(file.path(fig_dir, "o2_timecourses_by_regime.pdf"), width = 9, height = 6)
    graphics::plot(NA, xlim = range(d$day, na.rm = TRUE), ylim = range(d$o2_pct, na.rm = TRUE), xlab = "Day", ylab = "O2 pct", main = "O2 timecourses by inferred regime")
    for (key in split(seq_len(nrow(d)), paste(d$seed_id, d$cohort, sep = "||"))) {
      x <- d[key, , drop = FALSE]
      graphics::lines(x$day, x$o2_pct, col = color_for(x$trajectory_regime[[1]], 0.25), lwd = 1)
    }
    graphics::legend("topright", legend = names(cols), col = cols, lwd = 2, bty = "n")
    grDevices::dev.off()
  }
  event_ok <- nrow(events) && any(
    is.finite(events$time_o2_near_floor) &
      is.finite(events$time_crossing_ploidy_1p5_down_after_late_start)
  )
  if (event_ok) {
    grDevices::pdf(file.path(fig_dir, "late_ploidy_drop_vs_o2_floor_time.pdf"), width = 8, height = 6)
    graphics::plot(
      events$time_o2_near_floor, events$time_crossing_ploidy_1p5_down_after_late_start,
      xlab = "Time O2 reaches near-floor", ylab = "Time ploidy crosses 1.5 after late start",
      main = "Late ploidy collapse timing vs O2 floor", pch = ifelse(events$cohort == "2N", 16, 17),
      col = vapply(events$trajectory_regime, color_for, character(1))
    )
    graphics::abline(0, 1, lty = 2, col = "grey40")
    graphics::legend("topleft", legend = c("2N", "4N"), pch = c(16, 17), bty = "n")
    grDevices::dev.off()
  }
  figures <- list.files(fig_dir, pattern = "[.]pdf$", full.names = FALSE)
  if (!length(figures)) stop("Event visualization produced no figures.")
  manifest <- data.frame(
    artifact = tools::file_path_sans_ext(figures), relative_path = file.path("figures", figures),
    role = "figure", path = normalizePath(file.path(fig_dir, figures), mustWork = TRUE),
    exists = TRUE, stringsAsFactors = FALSE
  )
  pfv_write_tsv(manifest, file.path(out_dir, "visualization_manifest.tsv"))
  invisible(out_dir)
}

if (sys.nframe() == 0L) run_o2_ploidy_event_visualization()
