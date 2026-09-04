#!/usr/bin/env Rscript
# Layout-only preview; never publishes production outputs or changes pointers.
source("Code/Figures/draw_Figure7.R")
JOINT_FAMILY_LEVELS <- c("C01", "C02")
paths <- f7r_paths()
run <- f7g_paths(paths)
objects <- list()
for (key in c("C", "D", "E", "F")) {
  context <- if (key %in% c("C", "D")) "in vivo" else "in vitro"
  mode <- if (context == "in vivo") "continuous" else "passage"
  family <- if (key %in% c("C", "E")) "C01" else "C02"
  objects[[key]] <- readRDS(file.path(run$run_root, f7g_panel_filename(context, mode, family)))
}
layout <- f7ab_build(paths, objects)
output <- file.path(paths$root, "audit", "figure7_segment_revision_20260904", "layout_preview_NOT_FINAL")
f7r_save_plot(layout$plot, output, layout$width, layout$height)
f7ft_atomic_write_tsv(layout$geometry, paste0(output, "_geometry.tsv"))
