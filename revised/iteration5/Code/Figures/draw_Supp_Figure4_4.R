#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", arg[1L])))
root <- normalizePath(file.path(here, "../.."))
source(file.path(here, "util/analysis/supp_figure4_4.R"))
supp4_4_assert(basename(root) == "iteration5", "Only iteration5 is supported.")
data_dir <- file.path(root, "data/Figures/Supp_Figure4_4")
output <- file.path(root, "Figures")
dt <- fread(file.path(data_dir, "timecourse_median_iqr.tsv"))
observed <- fread(file.path(data_dir, "observed_tumor_size.tsv"))
supp4_4_assert(nrow(dt) == 447L && all(dt$n_seeds == 500L), "Plot requires complete 500-seed summaries.")
dt[, cohort := factor(cohort, levels = c("2N", "4N"))]
observed[, cohort := factor(cohort, levels = c("2N", "4N"))]
colors <- c(ploidy = "#7B3294", effective_o2_pct = "#E69F00", tumor_size_mm3 = "#303840")
titles <- c(ploidy = "A  Mean viable ploidy", effective_o2_pct = "B  Effective O2", tumor_size_mm3 = "C  Tumor size")
y_titles <- c(ploidy = "Mean viable ploidy (N/22)", effective_o2_pct = "Effective O2 (%)", tumor_size_mm3 = "Total tumor volume (mm³)")
y_max <- c(ploidy = max(5, ceiling(max(dt[metric == "ploidy", q75]))),
           effective_o2_pct = 5,
           tumor_size_mm3 = max(5000, ceiling(max(c(dt[metric == "tumor_size_mm3", q75], observed$observed_mean_mm3)) / 1000) * 1000))
panels <- lapply(supp4_4_metrics, function(metric_name) {
  dat <- dt[metric == metric_name]
  ymax <- unname(y_max[[metric_name]])
  membership <- data.table(cohort = factor("4N", levels = c("2N", "4N")), day = 77.5)
  label <- data.table(cohort = factor("4N", levels = c("2N", "4N")), x = 98, y = ymax * 0.96,
                      text = "4 to 2 tumors")
  p <- ggplot(dat, aes(day, median, group = interaction(cohort, segment))) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = colors[[metric_name]], alpha = 0.22, color = NA) +
    geom_line(color = colors[[metric_name]], linewidth = 0.85) +
    geom_vline(data = membership, aes(xintercept = day), inherit.aes = FALSE,
               linetype = "dashed", linewidth = 0.38, color = "#777777") +
    geom_text(data = label, aes(x, y, label = text), inherit.aes = FALSE,
              hjust = 1, vjust = 1, size = 3.05, color = "#555555", family = "Arial") +
    facet_grid(rows = vars(cohort), drop = FALSE) +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.025))) +
    labs(title = titles[[metric_name]], x = "Day", y = y_titles[[metric_name]]) +
    theme_bw(base_size = 11, base_family = "Arial") +
    theme(text = element_text(color = "#202428"),
          plot.title = element_text(face = "bold", size = 13, margin = margin(b = 9)),
          axis.title = element_text(face = "bold", size = 11),
          axis.text = element_text(size = 10, color = "#252A2F"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#E5E6E8", linewidth = 0.28),
          panel.border = element_rect(color = "#92979D", linewidth = 0.4),
          strip.background = element_rect(fill = "#ECECEC", color = "#B8B8B8", linewidth = 0.35),
          strip.text.y = element_text(face = "bold", size = 11),
          panel.spacing.y = grid::unit(12, "pt"),
          plot.margin = margin(6, 8, 6, 6))
  if (metric_name == "tumor_size_mm3") {
    p <- p + geom_point(data = observed, aes(day, observed_mean_mm3), inherit.aes = FALSE,
                        shape = 21, size = 1.85, stroke = 0.6, fill = "white", color = "#111111")
  }
  built <- ggplot_build(p)
  supp4_4_assert(nrow(built$data[[1L]]) == 149L && all(is.finite(built$data[[1L]]$ymin)) &&
                  all(is.finite(built$data[[1L]]$ymax)), paste("Missing/clipped IQR band for", metric_name))
  p
})
figure <- wrap_plots(panels, nrow = 1) +
  plot_annotation(
    title = "In-vivo dynamics across 500 fitted seeds",
    subtitle = "Solid line: median; shaded band: interquartile range (25th-75th percentiles) across all 500 optimizer endpoints",
    caption = paste(
      "Within each seed, tumors are weighted equally within the 2N or 4N cohort. Open points: observed cohort-mean tumor volume.",
      "Original follow-up is retained: 2N to day 56; 4N to day 91. In 4N, model membership changes from 4 tumors to 2 after day 77 (dashed line; curves split).",
      "Bands describe variation across fitted optimizer endpoints, not biological confidence intervals or posterior uncertainty. No fit-quality filtering or extrapolation.",
      sep = "\n"),
    theme = theme(text = element_text(family = "Arial", color = "#202428"),
                  plot.title = element_text(face = "bold", size = 18, hjust = 0, margin = margin(b = 6)),
                  plot.subtitle = element_text(size = 11.2, hjust = 0, margin = margin(b = 12)),
                  plot.caption = element_text(size = 8.6, hjust = 0, lineheight = 1.15, margin = margin(t = 12)),
                  plot.margin = margin(12, 15, 12, 12)))
stem <- file.path(output, "Supp_Figure4_4")
ggsave(paste0(stem, ".pdf"), figure, width = 15, height = 8, units = "in", device = cairo_pdf, bg = "white")
ggsave(paste0(stem, ".png"), figure, width = 15, height = 8, units = "in", dpi = 300, bg = "white")
ggsave(paste0(stem, ".svg"), figure, width = 15, height = 8, units = "in", device = svglite::svglite, bg = "white")
fwrite(data.table(metric = c("figure_width_in", "figure_height_in", "n_metric_columns", "n_cohort_rows", "n_seeds_per_point",
                             "ribbon_lower", "center", "ribbon_upper", "ploidy_color", "effective_o2_color", "tumor_size_color",
                             "model_membership_break_day", "observations_used_in_seed_quantiles", "fixed_o2_scan_used", "existing_figure4_modified"),
                  value = c(15, 8, 3, 2, 500, "Q25", "Q50", "Q75", colors, 77.5, FALSE, FALSE, FALSE)),
       file.path(data_dir, "plot_validation.tsv"), sep = "\t")
cat("supp_figure4_4_render_ok\n")
