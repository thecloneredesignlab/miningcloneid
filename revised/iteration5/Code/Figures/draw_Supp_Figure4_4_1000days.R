#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", arg[1L])))
root <- normalizePath(file.path(here, "../.."))
source(file.path(here, "util/analysis/supp_figure4_4.R"))
supp4_4_assert(basename(root) == "iteration5", "Only iteration5 is supported.")
data_dir <- file.path(root, "data/Figures/Supp_Figure4_4_1000days")
quantile_columns <- c("q25", "median", "q75")
dt <- fread(file.path(data_dir, "timecourse_median_iqr.tsv"), colClasses = list(character = quantile_columns))
for (column in quantile_columns) {
  value <- suppressWarnings(as.numeric(dt[[column]]))
  supp4_4_assert(all(is.finite(value)), paste("Invalid numeric quantile", column))
  set(dt, j = column, value = value)
}
supp4_4_assert(all(dt$q25 <= dt$median & dt$median <= dt$q75), "Invalid quantile ordering.")
observed <- fread(file.path(data_dir, "observed_tumor_size.tsv"))
supp4_4_assert(nrow(dt) == 6006L && all(dt$n_seeds == 500L) && all(dt$n_tumors == 4L), "Incomplete forecast summaries.")
dt[, cohort := factor(cohort, levels = c("2N", "4N"))]
observed[, cohort := factor(cohort, levels = c("2N", "4N"))]
ends <- data.table(cohort = factor(c("2N", "4N"), levels = c("2N", "4N")), day = c(56, 91))
colors <- c(ploidy = "#7B3294", effective_o2_pct = "#E69F00", tumor_size_mm3 = "#303840")
titles <- c(ploidy = "A  Mean viable ploidy", effective_o2_pct = "B  Effective O2", tumor_size_mm3 = "C  Tumor size")
y_titles <- c(ploidy = "Mean viable ploidy (N/22)", effective_o2_pct = "Effective O2 (%)", tumor_size_mm3 = "Total tumor volume (mm³; pseudo-log10)")
volume_max <- max(c(dt[metric == "tumor_size_mm3", q75], observed$observed_mean_mm3))
top_power <- max(5, ceiling(log10(volume_max)))
volume_breaks <- c(0, 10^unique(c(5 * seq_len(floor(top_power / 5)), top_power)))
volume_labels <- function(x) ifelse(x == 0, "0", ifelse(x == 1, "1", paste0("1e", round(log10(pmax(x, 1))))))
panels <- lapply(supp4_4_metrics, function(metric_name) {
  dat <- dt[metric == metric_name]
  p <- ggplot(dat, aes(day, median, group = cohort)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = colors[[metric_name]], alpha = .22, color = NA) +
    geom_line(color = colors[[metric_name]], linewidth = .85) +
    geom_vline(data = ends, aes(xintercept = day), linetype = "dashed", color = "#777777", linewidth = .38) +
    facet_grid(rows = vars(cohort), drop = FALSE) +
    scale_x_continuous(limits = c(0, 1000), breaks = c(0, 250, 500, 750, 1000), expand = expansion(mult = c(.005, .005))) +
    labs(title = titles[[metric_name]], x = "Day", y = y_titles[[metric_name]]) +
    theme_bw(base_size = 11, base_family = "Arial") +
    theme(text = element_text(color = "#202428"), plot.title = element_text(face = "bold", size = 13, margin = margin(b = 9)),
          axis.title = element_text(face = "bold", size = 11), axis.text = element_text(size = 10, color = "#252A2F"),
          panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#E5E6E8", linewidth = .28),
          panel.border = element_rect(color = "#92979D", linewidth = .4),
          strip.background = element_rect(fill = "#ECECEC", color = "#B8B8B8", linewidth = .35),
          strip.text.y = element_text(face = "bold", size = 11), panel.spacing.y = grid::unit(12, "pt"),
          plot.margin = margin(6, 8, 6, 6))
  if (metric_name == "tumor_size_mm3") {
    p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                                limits = c(0, 10^top_power), breaks = volume_breaks, labels = volume_labels,
                                expand = expansion(mult = c(.015, .025))) +
      geom_point(data = observed, aes(day, observed_mean_mm3), inherit.aes = FALSE,
                  shape = 21, size = 1.65, stroke = .55, fill = "white", color = "#111111")
  } else {
    ymax <- if (metric_name == "ploidy") max(5, ceiling(max(dat$q75))) else 5
    p <- p + scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(.01, .025)))
  }
  built <- ggplot_build(p)
  supp4_4_assert(nrow(built$data[[1]]) == 2002L && all(is.finite(built$data[[1]]$ymin)) && all(is.finite(built$data[[1]]$ymax)), "Missing/clipped forecast IQR bands.")
  p
})
figure <- wrap_plots(panels, nrow = 1) + plot_annotation(
  title = "In-vivo dynamics across 500 fitted seeds: 1,000-day model extrapolation",
  subtitle = "Solid line: median; shaded band: interquartile range (25th-75th percentiles) across all 500 optimizer endpoints",
  caption = paste(
    "Within each seed, the four original untreated tumors per cohort are continued to day 1,000 and weighted equally. Open points: observed cohort-mean tumor volume.",
    "Dashed lines mark the last observed follow-up (2N: day 56; 4N: day 91). Later values are unvalidated model extrapolations, not observed long-term tumor growth.",
    "Bands describe fitted-endpoint variation, not biological confidence intervals or posterior uncertainty. Tumor volume uses pseudo-log10 (sigma = 1 mm³) to retain zero.", sep = "\n"),
  theme = theme(text = element_text(family = "Arial", color = "#202428"),
                plot.title = element_text(face = "bold", size = 17, hjust = 0, margin = margin(b = 6)),
                plot.subtitle = element_text(size = 11.2, hjust = 0, margin = margin(b = 12)),
                plot.caption = element_text(size = 8.3, hjust = 0, lineheight = 1.15, margin = margin(t = 12)),
                plot.margin = margin(12, 15, 12, 12)))
stem <- file.path(root, "Figures/Supp_Figure4_4_1000days")
ggsave(paste0(stem, ".pdf"), figure, width = 15, height = 8, device = cairo_pdf, bg = "white")
ggsave(paste0(stem, ".png"), figure, width = 15, height = 8, dpi = 300, bg = "white")
ggsave(paste0(stem, ".svg"), figure, width = 15, height = 8, device = svglite::svglite, bg = "white")
fwrite(data.table(metric = c("horizon_days", "n_seeds", "width_in", "height_in", "volume_axis", "volume_sigma_mm3", "ploidy_o2_axes", "original_short_window_modified"),
                  value = c(1000, 500, 15, 8, "pseudo-log10", 1, "linear", FALSE)), file.path(data_dir, "plot_validation.tsv"), sep = "\t")
cat("supp_figure4_4_1000days_render_ok\n")
