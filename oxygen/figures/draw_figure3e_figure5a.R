#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(magick)
  library(patchwork)
})

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)) else getwd()
})
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
result_dir <- file.path(repo_root, "oxygen", "results", "figure3e_figure5a")
iteration_dir <- file.path(script_dir, "iteration1")
dir.create(iteration_dir, recursive = TRUE, showWarnings = FALSE)

figure3e_data_path <- file.path(result_dir, "figure3e_4n_deprived_passage_summary.csv")
figure5a_pairs_path <- file.path(result_dir, "figure5a_joint_warmup_pairs.csv")
figure5a_soft_path <- file.path(result_dir, "figure5a_soft_coupled_parameters.csv")
required <- c(figure3e_data_path, figure5a_pairs_path, figure5a_soft_path)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop(
    "Missing selected result tables: ", paste(missing, collapse = ", "),
    ". Run oxygen/scripts/import_figure3e_figure5a_results.R first.",
    call. = FALSE
  )
}

ink <- "#20242A"
muted <- "#606873"
border <- "#C9CED6"
panel_fill <- "#F7F8FA"
blue <- "#2C7BB6"
blue_light <- "#DCEEF7"
red <- "#C43C39"
red_light <- "#F8E0DE"
amber <- "#C98A00"
amber_light <- "#FFF0C7"
green <- "#2E7D5B"
green_light <- "#DCEFE7"

theme_manuscript <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = base_size + 1.5, face = "bold", color = ink),
      plot.subtitle = element_text(size = base_size - 0.4, color = muted),
      axis.title = element_text(size = base_size, color = ink),
      axis.text = element_text(size = base_size - 0.7, color = ink),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E8EBEF", linewidth = 0.3),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.8),
      plot.margin = margin(7, 10, 7, 7)
    )
}

draw_figure3e <- function(data_path) {
  d <- utils::read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_columns <- c(
    "lineage_passage", "oxygen_pct", "predicted_mean_chromosome",
    "direct_hypoxia_death_burden_fraction", "nonviable_daughter_burden_fraction"
  )
  if (!all(required_columns %in% names(d)) || nrow(d) < 20L || anyDuplicated(d$lineage_passage)) {
    stop("Figure 3E input does not contain one complete terminal 4N lineage.", call. = FALSE)
  }
  d <- d[order(d$lineage_passage), , drop = FALSE]
  shaded <- d[d$oxygen_pct < 20.5, , drop = FALSE]
  initial_mean <- d$predicted_mean_chromosome[[1L]]
  final_mean <- tail(d$predicted_mean_chromosome, 1L)
  max_direct <- max(d$direct_hypoxia_death_burden_fraction)

  p_ploidy <- ggplot(d, aes(lineage_passage, predicted_mean_chromosome)) +
    geom_rect(
      data = shaded,
      aes(xmin = lineage_passage - 0.48, xmax = lineage_passage + 0.48, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = blue_light, color = NA, alpha = 0.58
    ) +
    geom_line(color = blue, linewidth = 0.9) +
    geom_point(aes(fill = oxygen_pct < 20.5), shape = 21, size = 2.1, stroke = 0.5, color = blue) +
    scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = blue), guide = "none") +
    annotate(
      "text", x = d$lineage_passage[[1L]], y = initial_mean,
      label = sprintf("%.1f", initial_mean), hjust = -0.2, vjust = -0.8,
      color = blue, size = 2.8, family = "Helvetica"
    ) +
    annotate(
      "text", x = tail(d$lineage_passage, 1L), y = final_mean,
      label = sprintf("%.1f", final_mean), hjust = 1.2, vjust = 1.2,
      color = blue, size = 2.8, family = "Helvetica"
    ) +
    labs(
      x = NULL,
      y = "Predicted mean\nchromosome number\n(focused scale)",
      title = "4N severe-deprivation lineage: predicted ploidy and death components",
      subtitle = "Best separate in vitro fit (seed 10); blue shading marks passages below 20.5% O2"
    ) +
    scale_x_continuous(breaks = seq(1, max(d$lineage_passage), by = 4), expand = expansion(mult = c(0.015, 0.015))) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.13))) +
    theme_manuscript(9)

  fractions <- rbind(
    data.frame(
      lineage_passage = d$lineage_passage,
      component = "Direct hypoxia-origin death",
      burden_fraction = d$direct_hypoxia_death_burden_fraction
    ),
    data.frame(
      lineage_passage = d$lineage_passage,
      component = "Nonviable chromosome-transition products",
      burden_fraction = d$nonviable_daughter_burden_fraction
    )
  )
  component_colors <- c(
    "Direct hypoxia-origin death" = red,
    "Nonviable chromosome-transition products" = green
  )
  shade_y_min <- min(fractions$burden_fraction) * 0.8
  shade_y_max <- max(fractions$burden_fraction) * 1.2
  p_death <- ggplot(fractions, aes(lineage_passage, burden_fraction, color = component)) +
    geom_rect(
      data = shaded,
      aes(xmin = lineage_passage - 0.48, xmax = lineage_passage + 0.48),
      ymin = shade_y_min, ymax = shade_y_max,
      inherit.aes = FALSE, fill = blue_light, color = NA, alpha = 0.58
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    annotate(
      "label", x = max(d$lineage_passage) - 0.3, y = max_direct,
      label = sprintf("Direct death <= %.1f%%", 100 * max_direct),
      hjust = 1, vjust = -0.5, size = 2.7, family = "Helvetica",
      color = red, fill = "white", linewidth = 0.2
    ) +
    scale_color_manual(values = component_colors) +
    scale_y_log10(
      breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 3e-1),
      labels = c("0.01%", "0.1%", "1%", "10%", "30%")
    ) +
    scale_x_continuous(breaks = seq(1, max(d$lineage_passage), by = 4), expand = expansion(mult = c(0.015, 0.015))) +
    labs(
      x = "Lineage passage",
      y = "Fraction of total burden\n(log scale)"
    ) +
    theme_manuscript(9) +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())

  p_ploidy / p_death + plot_layout(heights = c(1.05, 1))
}

gp_text <- function(size = 9, col = ink, fontface = "plain") {
  gpar(fontsize = size, col = col, fontface = fontface, fontfamily = "Helvetica")
}

draw_box <- function(x, y, w, h, label, fill = "white", line = border,
                     size = 8.2, fontface = "plain") {
  grid.roundrect(
    x, y, w, h, r = unit(2, "mm"),
    gp = gpar(fill = fill, col = line, lwd = 1.05)
  )
  grid.text(label, x, y, gp = gp_text(size, ink, fontface))
}

draw_arrow <- function(x0, y0, x1, y1, color = muted, lwd = 1.4) {
  grid.lines(
    unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
    arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
    gp = gpar(col = color, lwd = lwd, lineend = "round")
  )
}

draw_figure5a <- function(pair_path, soft_path) {
  pairs <- utils::read.csv(pair_path, check.names = FALSE, stringsAsFactors = FALSE)
  soft <- utils::read.csv(soft_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(pairs) != 6L || nrow(soft) != 14L || !all(pairs$optimizer_seed_count == 500L)) {
    stop("Figure 5A setup tables do not match the validated six-pair, 14-parameter design.", call. = FALSE)
  }
  pair_labels <- vapply(sort(unique(pairs$invivo_cluster)), function(cluster) {
    z <- pairs[pairs$invivo_cluster == cluster, , drop = FALSE]
    paste0("C", cluster, ": seeds ", paste(z$invivo_seed[order(z$invivo_subcluster)], collapse = ", "))
  }, character(1L))
  process_counts <- table(soft$primary_process)
  process_order <- c(
    "Hypoxia sensing/resource stress", "Growth", "Hypoxia-associated death",
    "CIN/missegregation", "Whole-genome doubling", "Post-missegregation survival/buffering"
  )
  short_process <- c(
    "Resource sensing", "Growth", "Stress-associated death",
    "CIN / missegregation", "WGD", "Post-MS survival"
  )
  process_labels <- paste0(short_process, " (", as.integer(process_counts[process_order]), ")")

  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text("Joint soft-coupled fitting design", 0.04, 0.94, just = c("left", "top"), gp = gp_text(12, ink, "bold"))
  grid.text(
    "Six landscape-informed warm starts; each joint run uses 500 optimizer seeds",
    0.04, 0.875, just = c("left", "top"), gp = gp_text(8.2, muted)
  )

  grid.text("Warm-start selection", 0.17, 0.78, gp = gp_text(9.2, ink, "bold"))
  draw_box(
    0.17, 0.61, 0.25, 0.22,
    paste(c("In vivo separate-fit landscape", pair_labels), collapse = "\n"),
    blue_light, blue, 7.8, "bold"
  )
  draw_box(
    0.17, 0.31, 0.25, 0.16,
    "In vitro separate fit\nbest fit: seed 10\nshared across all pairs",
    green_light, green, 8.0, "bold"
  )

  draw_arrow(0.30, 0.61, 0.355, 0.61, blue)
  draw_arrow(0.30, 0.31, 0.355, 0.40, green)

  grid.roundrect(
    0.56, 0.48, 0.38, 0.66, r = unit(2.5, "mm"),
    gp = gpar(fill = panel_fill, col = border, lwd = 1.2)
  )
  grid.text("Matched joint objective", 0.56, 0.755, gp = gp_text(9.5, ink, "bold"))
  draw_box(0.46, 0.645, 0.16, 0.12, "Tumor data\nburden | chromosome\nnecrosis", blue_light, blue, 7.4, "bold")
  draw_box(0.66, 0.645, 0.16, 0.12, "Culture data\ngrowth | chromosome\nflow density", green_light, green, 7.4, "bold")
  grid.text("One mechanistic model; context-resolved parameter values", 0.56, 0.545, gp = gp_text(7.8, muted, "bold"))
  grid.text("theta_vivo = center + delta / 2", 0.48, 0.485, gp = gp_text(7.7, blue, "bold"))
  grid.text("theta_vitro = center - delta / 2", 0.65, 0.485, gp = gp_text(7.7, green, "bold"))
  grid.text(
    sprintf("14 soft-coupled parameters | Welsch penalty: sigma = %.2f, c = %.1f", pairs$soft_coupling_sigma_default[[1L]], pairs$soft_coupling_welsch_c[[1L]]),
    0.56, 0.420, gp = gp_text(7.7, ink, "bold")
  )

  chip_x <- c(0.445, 0.56, 0.675)
  chip_y <- c(0.335, 0.245)
  chip_fill <- c(blue_light, amber_light, red_light, amber_light, blue_light, green_light)
  chip_line <- c(blue, amber, red, amber, blue, green)
  for (i in seq_along(process_labels)) {
    row <- if (i <= 3L) 1L else 2L
    col <- if (i <= 3L) i else i - 3L
    draw_box(chip_x[[col]], chip_y[[row]], 0.105, 0.072, process_labels[[i]], chip_fill[[i]], chip_line[[i]], 5.9, "bold")
  }

  draw_arrow(0.755, 0.48, 0.80, 0.48, muted)
  grid.text("Context contrasts", 0.88, 0.78, gp = gp_text(9.2, ink, "bold"))
  draw_box(0.88, 0.61, 0.18, 0.16, "Six joint warm-start pairs\nx 500 optimizer seeds\n(best fits retained)", panel_fill, border, 7.8, "bold")
  draw_arrow(0.88, 0.52, 0.88, 0.45, muted)
  draw_box(0.88, 0.36, 0.18, 0.12, "Context-specific\nresponse functions\n(Fig. 5B-D)", blue_light, blue, 7.6, "bold")
  draw_box(0.88, 0.19, 0.18, 0.12, "In vivo / in vitro\nparameter contrasts\n(Fig. 5E)", green_light, green, 7.6, "bold")
}

save_plot_pair <- function(plot, stem, width, height) {
  png_path <- file.path(iteration_dir, paste0(stem, ".png"))
  pdf_path <- file.path(iteration_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white")
  c(png_path, pdf_path)
}

save_grid_pair <- function(draw_function, stem, width, height) {
  png_path <- file.path(iteration_dir, paste0(stem, ".png"))
  pdf_path <- file.path(iteration_dir, paste0(stem, ".pdf"))
  grDevices::png(png_path, width = width * 300, height = height * 300, res = 300, bg = "white", type = "cairo")
  draw_function()
  grDevices::dev.off()
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "Helvetica")
  draw_function()
  grDevices::dev.off()
  c(png_path, pdf_path)
}

label_panel <- function(path, label, width_px) {
  img <- image_read(path)
  img <- image_scale(img, paste0(width_px, "x"))
  image_annotate(
    img, label, gravity = "northwest", location = "+22+12",
    size = max(50, round(width_px * 0.055)), color = "#111111",
    boxcolor = "#FFFFFFDD", font = "Helvetica"
  )
}

spacer_image <- function(width, height) image_blank(width, height, color = "white")

make_row <- function(images, gap = 36L) {
  heights <- vapply(images, function(img) image_info(img)$height[[1L]], integer(1L))
  max_height <- max(heights)
  fitted <- lapply(images, function(img) image_extent(img, paste0(image_info(img)$width[[1L]], "x", max_height), gravity = "center", color = "white"))
  with_gaps <- list()
  for (i in seq_along(fitted)) {
    if (i > 1L) with_gaps[[length(with_gaps) + 1L]] <- spacer_image(gap, max_height)
    with_gaps[[length(with_gaps) + 1L]] <- fitted[[i]]
  }
  image_append(do.call(c, with_gaps), stack = FALSE)
}

stack_rows <- function(rows, gap = 44L) {
  widths <- vapply(rows, function(img) image_info(img)$width[[1L]], integer(1L))
  max_width <- max(widths)
  fitted <- lapply(rows, function(img) image_extent(img, paste0(max_width, "x", image_info(img)$height[[1L]]), gravity = "center", color = "white"))
  with_gaps <- list()
  for (i in seq_along(fitted)) {
    if (i > 1L) with_gaps[[length(with_gaps) + 1L]] <- spacer_image(max_width, gap)
    with_gaps[[length(with_gaps) + 1L]] <- fitted[[i]]
  }
  image_append(do.call(c, with_gaps), stack = TRUE)
}

write_assembled <- function(image, stem) {
  png_path <- file.path(script_dir, paste0(stem, ".png"))
  pdf_path <- file.path(script_dir, paste0(stem, ".pdf"))
  image_write(image, png_path, format = "png", density = "300x300")

  info <- image_info(image)
  pdf_width <- 7
  pdf_height <- pdf_width * info$height[[1L]] / info$width[[1L]]
  grDevices::cairo_pdf(pdf_path, width = pdf_width, height = pdf_height, bg = "white")
  grid.newpage()
  grid.raster(
    as.raster(image_read(png_path)),
    x = unit(0.5, "npc"), y = unit(0.5, "npc"),
    width = unit(1, "npc"), height = unit(1, "npc"),
    interpolate = TRUE
  )
  grDevices::dev.off()

  pdftoppm <- Sys.which("pdftoppm")
  if (nzchar(pdftoppm)) {
    render_prefix <- tempfile(paste0(stem, "_pdf_check_"))
    render_path <- paste0(render_prefix, ".png")
    render_log <- system2(
      pdftoppm,
      c("-f", "1", "-singlefile", "-png", "-r", "72", pdf_path, render_prefix),
      stdout = TRUE, stderr = TRUE
    )
    if (!file.exists(render_path)) {
      stop("Could not render assembled PDF for validation: ", paste(render_log, collapse = "\n"), call. = FALSE)
    }
    rendered <- image_data(image_read(render_path), channels = "gray")
    nonwhite_fraction <- mean(as.integer(rendered) < 250L)
    unlink(render_path)
    if (!is.finite(nonwhite_fraction) || nonwhite_fraction < 0.005) {
      stop(sprintf("Assembled PDF %s rendered blank (nonwhite fraction %.5f).", pdf_path, nonwhite_fraction), call. = FALSE)
    }
  }
  c(png_path, pdf_path)
}

figure3e <- draw_figure3e(figure3e_data_path)
new_figure3_paths <- save_plot_pair(figure3e, "fig3e_rejected_strong_elimination", 7.0, 5.4)
new_figure5_paths <- save_grid_pair(
  function() draw_figure5a(figure5a_pairs_path, figure5a_soft_path),
  "fig5a_joint_fit_setup", 10.0, 4.0
)

fig3_top <- make_row(list(
  label_panel(file.path(iteration_dir, "fig3a_invitro_growth_ploidy_burden_fit.png"), "A", 1800L),
  label_panel(file.path(iteration_dir, "fig3b_invitro_predicted_ploidy_distribution.png"), "B", 1800L)
))
fig3_bottom <- make_row(list(
  label_panel(file.path(iteration_dir, "fig3c_invitro_viability_after_missegregation.png"), "C", 1188L),
  label_panel(file.path(iteration_dir, "fig3d_invitro_nonviable_fraction_vs_ms_rate.png"), "D", 1188L),
  label_panel(new_figure3_paths[[1L]], "E", 1188L)
))
assembled3 <- stack_rows(list(fig3_top, fig3_bottom))
assembled3_paths <- write_assembled(assembled3, "assembled_fig3")

fig5_top <- label_panel(new_figure5_paths[[1L]], "A", 3600L)
fig5_middle <- make_row(list(
  label_panel(file.path(iteration_dir, "fig5b_joint_o2_vs_proliferation_rate.png"), "B", 1782L),
  label_panel(file.path(iteration_dir, "fig5c_joint_o2_vs_missegregation_rate.png"), "C", 1782L)
))
fig5_bottom <- make_row(list(
  label_panel(file.path(iteration_dir, "fig5d_joint_post_ms_survival_in_vivo_vs_invitro.png"), "D", 1782L),
  label_panel(file.path(iteration_dir, "fig5e_joint_parameter_ratio_all_soft.png"), "E", 1782L)
))
assembled5 <- stack_rows(list(fig5_top, fig5_middle, fig5_bottom))
assembled5_paths <- write_assembled(assembled5, "assembled_fig5")

message("Wrote Figure 3E: ", paste(new_figure3_paths, collapse = ", "))
message("Wrote Figure 5A: ", paste(new_figure5_paths, collapse = ", "))
message("Rebuilt Figure 3: ", paste(assembled3_paths, collapse = ", "))
message("Rebuilt Figure 5: ", paste(assembled5_paths, collapse = ", "))
