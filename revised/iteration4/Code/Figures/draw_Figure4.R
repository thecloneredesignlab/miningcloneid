#!/usr/bin/env Rscript

if (identical(Sys.getenv("FIGURE4_DRAW_WORKER"), "1")) {

# Build the active Figure 4A as a 2 x 2 composite:
#   left column (70%): cohort-level burden/O2/ploidy dynamics;
#   right column (30%): terminal observed-vs-predicted ploidy distributions.
# The two rows are the 2N and 4N cohorts. All chromosome-number values are
# converted to ploidy as N / N_UNIT, with N_UNIT = 22.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyr))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
script_dir <- dirname(script_path)
data_dir <- normalizePath(
  Sys.getenv("FIGURE4_DATA_DIR", unset = file.path(script_dir, "data")),
  mustWork = TRUE
)
figure_dir <- normalizePath(
  Sys.getenv(
    "FIGURE4_PANEL_DIR",
    unset = file.path(script_dir, "figures")
  ),
  mustWork = FALSE
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

burden_path <- file.path(data_dir, "figure4a_seed25_burden_timecourse.tsv")
chromosome_path <- file.path(
  data_dir, "figure4a_seed25_weighted_mean_chromosome_timecourse.tsv"
)
terminal_path <- file.path(
  data_dir, "figure4a_seed25_terminal_chromosome_observed_vs_predicted.tsv"
)
burden_audit_path <- file.path(
  data_dir, "figure4a_burden_measurement_audit.tsv"
)

required_paths <- c(
  burden_path, chromosome_path, terminal_path, burden_audit_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing Figure 4A input(s): ", paste(missing_paths, collapse = ", "))
}

burden <- read.delim(burden_path, check.names = FALSE, stringsAsFactors = FALSE)
chromosome <- read.delim(
  chromosome_path, check.names = FALSE, stringsAsFactors = FALSE
)
terminal <- read.delim(
  terminal_path, check.names = FALSE, stringsAsFactors = FALSE
)
burden_audit <- read.delim(
  burden_audit_path, check.names = FALSE, stringsAsFactors = FALSE
)

required_burden <- c(
  "harvest", "cohort", "dose", "day", "pred_burden_live_volume_mm3",
  "pred_burden_dead_hypoxia_volume_mm3",
  "pred_burden_dead_buffer_volume_mm3", "pred_burden_volume_mm3",
  "pred_o2_pct", "obs_burden"
)
required_chromosome <- c("harvest", "cohort", "dose", "day", "weighted_mean_N")
required_terminal <- c(
  "harvest", "cohort", "source", "endpoint_value", "weight"
)
required_audit <- c(
  "harvest", "cohort", "day", "figure4_obs_burden_mm3",
  "figure1_obs_burden_mm3", "exact_package_internal_match",
  "cohort_day_n", "observed_line_segment"
)
if (!all(required_burden %in% names(burden))) {
  stop("Burden time course is missing required columns.")
}
if (!all(required_chromosome %in% names(chromosome))) {
  stop("Chromosome time course is missing required columns.")
}
if (!all(required_terminal %in% names(terminal))) {
  stop("Terminal comparison table is missing required columns.")
}
if (!all(required_audit %in% names(burden_audit))) {
  stop("Figure 4A burden audit is missing required columns.")
}
if (any(!as.logical(burden_audit$exact_package_internal_match)) ||
    unique(burden_audit$cohort_day_n[
      burden_audit$cohort == "4N" & burden_audit$day == 77
    ]) != 4L ||
    unique(burden_audit$cohort_day_n[
      burden_audit$cohort == "4N" & burden_audit$day == 81
    ]) != 2L) {
  stop("Figure 4A burden audit does not confirm the day-77/day-81 measurements.")
}

n_unit <- 22
shared_axis_max <- 5
terminal_axis_display_max <- 5
burden_axis_max <- 5000
burden_to_shared <- burden_axis_max / shared_axis_max
cohort_levels <- c("2N", "4N")

key <- c("harvest", "cohort", "dose", "day")
trajectory <- inner_join(
  burden[, required_burden],
  chromosome[, required_chromosome],
  by = key
) %>%
  filter(day >= 0, day <= 100) %>%
  mutate(mean_viable_ploidy = weighted_mean_N / n_unit)

if (nrow(trajectory) != nrow(burden)) {
  stop("Burden and chromosome trajectory tables do not match one-to-one.")
}
if (length(unique(trajectory$harvest)) != 8L) {
  stop("Figure 4A requires exactly eight paired in-vivo trajectories.")
}
if (!setequal(unique(trajectory$cohort), cohort_levels)) {
  stop("Figure 4A requires exactly the 2N and 4N cohorts.")
}

predicted <- trajectory %>%
  group_by(cohort, day) %>%
  summarise(
    n_model_trajectories = n_distinct(harvest),
    burden_live_mm3 = mean(pred_burden_live_volume_mm3),
    burden_dead_hypoxia_mm3 = mean(pred_burden_dead_hypoxia_volume_mm3),
    burden_dead_cin_mm3 = mean(pred_burden_dead_buffer_volume_mm3),
    burden_total_mm3 = mean(pred_burden_volume_mm3),
    effective_o2_pct = mean(pred_o2_pct),
    mean_viable_ploidy = mean(mean_viable_ploidy),
    .groups = "drop"
  ) %>%
  mutate(cohort = factor(cohort, levels = cohort_levels))

observed <- trajectory %>%
  filter(is.finite(obs_burden)) %>%
  group_by(cohort, day) %>%
  summarise(
    n_observed = n_distinct(harvest),
    observed_burden_mean_mm3 = mean(obs_burden),
    observed_burden_sd_mm3 = ifelse(n() > 1L, sd(obs_burden), NA_real_),
    observed_burden_min_mm3 = min(obs_burden),
    observed_burden_max_mm3 = max(obs_burden),
    .groups = "drop"
  ) %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    observed_line_segment = ifelse(
      cohort == "4N" & day >= 81,
      "4N late n=2",
      ifelse(cohort == "4N", "4N early n=4", "2N n=4")
    ),
    observed_burden_lower_mm3 = pmax(
      observed_burden_mean_mm3 - observed_burden_sd_mm3, 0
    ),
    observed_burden_upper_mm3 =
      observed_burden_mean_mm3 + observed_burden_sd_mm3
  )

if (max(
  c(predicted$burden_total_mm3, observed$observed_burden_upper_mm3),
  na.rm = TRUE
) > burden_axis_max) {
  stop("Burden data exceed the prespecified 0--5000 mm3 display range.")
}
if (any(predicted$effective_o2_pct < 0 |
        predicted$effective_o2_pct > shared_axis_max) ||
    any(predicted$mean_viable_ploidy < 0 |
        predicted$mean_viable_ploidy > shared_axis_max)) {
  stop("O2 or mean ploidy lies outside the shared 0--5 display range.")
}

predicted_long <- predicted %>%
  pivot_longer(
    cols = c(
      burden_live_mm3,
      burden_dead_hypoxia_mm3,
      burden_dead_cin_mm3
    ),
    names_to = "component",
    values_to = "burden_mm3"
  ) %>%
  mutate(
    component = factor(
      component,
      levels = c(
        "burden_live_mm3",
        "burden_dead_hypoxia_mm3",
        "burden_dead_cin_mm3"
      ),
      labels = c("Viable", "Stress dead (model)", "CIN-associated dead")
    )
  )

trajectory_lines <- bind_rows(
  predicted %>%
    transmute(
      cohort, day, series = "Effective O2 (%)",
      shared_axis_value = effective_o2_pct
    ),
  predicted %>%
    transmute(
      cohort, day, series = "Mean viable ploidy",
      shared_axis_value = mean_viable_ploidy
    )
) %>%
  mutate(
    series = factor(
      series,
      levels = c("Effective O2 (%)", "Mean viable ploidy")
    )
  )

# Figure 2.2 stored chromosome number under endpoint_value because the active
# fit used start_with=chr_number. Each tumor is assigned equal total mass
# within its cohort. This avoids letting tumors with more observed cells
# dominate the cohort-level terminal comparison.
terminal_plot <- terminal %>%
  filter(
    cohort %in% cohort_levels,
    source %in% c("Observed", "Predicted"),
    is.finite(endpoint_value),
    is.finite(weight),
    weight > 0
  ) %>%
  group_by(cohort, source, harvest) %>%
  mutate(
    tumor_weight_total = sum(weight),
    equal_tumor_weight = weight / tumor_weight_total
  ) %>%
  ungroup() %>%
  group_by(cohort, source) %>%
  mutate(plot_weight = equal_tumor_weight / n_distinct(harvest)) %>%
  ungroup() %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    source = factor(source, levels = c("Observed", "Predicted")),
    source_x = as.numeric(source),
    terminal_ploidy = endpoint_value / n_unit
  )

terminal_weight_check <- terminal_plot %>%
  group_by(cohort, source) %>%
  summarise(total_weight = sum(plot_weight), .groups = "drop")
if (any(abs(terminal_weight_check$total_weight - 1) > 1e-10)) {
  stop("Terminal cohort/source distributions are not normalized to unit mass.")
}
if (any(terminal_plot$terminal_ploidy < 1 |
        terminal_plot$terminal_ploidy > 7)) {
  stop("Terminal ploidy lies outside the modeled 1--7 range.")
}

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- as.numeric(x[keep])
  w <- as.numeric(w[keep])
  if (!length(x)) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cumulative <- cumsum(w) / sum(w)
  vapply(
    probs,
    function(prob) {
      index <- which(cumulative >= prob)[1L]
      if (is.na(index)) x[[length(x)]] else x[[index]]
    },
    numeric(1)
  )
}

terminal_box <- terminal_plot %>%
  group_by(cohort, source) %>%
  summarise(
    q1 = weighted_quantile(terminal_ploidy, plot_weight, 0.25),
    median = weighted_quantile(terminal_ploidy, plot_weight, 0.50),
    q3 = weighted_quantile(terminal_ploidy, plot_weight, 0.75),
    raw_min = min(terminal_ploidy),
    raw_max = max(terminal_ploidy),
    .groups = "drop"
  ) %>%
  mutate(
    iqr = q3 - q1,
    ymin = pmax(raw_min, q1 - 1.5 * iqr),
    ymax = pmin(raw_max, q3 + 1.5 * iqr),
    source_x = as.numeric(source),
    box_x = source_x + 0.13,
    box_half_width = 0.055,
    whisker_half_width = 0.035
  )

half_violin_width <- 0.42
terminal_half_violin <- terminal_plot %>%
  group_by(cohort, source) %>%
  group_split() %>%
  lapply(function(dat) {
    normalized_weight <- dat$plot_weight / sum(dat$plot_weight)
    weighted_center <- sum(normalized_weight * dat$terminal_ploidy)
    variance_denominator <- 1 - sum(normalized_weight^2)
    weighted_sd <- sqrt(
      sum(
        normalized_weight *
          (dat$terminal_ploidy - weighted_center)^2
      ) / variance_denominator
    )
    weighted_iqr <- diff(weighted_quantile(
      dat$terminal_ploidy,
      normalized_weight,
      c(0.25, 0.75)
    ))
    bandwidth_scale <- min(weighted_sd, weighted_iqr / 1.34)
    if (!is.finite(bandwidth_scale) || bandwidth_scale <= 0) {
      bandwidth_scale <- weighted_sd
    }
    if (!is.finite(bandwidth_scale) || bandwidth_scale <= 0) {
      bandwidth_scale <- diff(range(dat$terminal_ploidy)) / 4
    }
    effective_n <- 1 / sum(normalized_weight^2)
    weighted_nrd0_bandwidth <-
      0.9 * bandwidth_scale * effective_n^(-1 / 5)
    density_fit <- stats::density(
      dat$terminal_ploidy,
      weights = normalized_weight,
      bw = weighted_nrd0_bandwidth,
      adjust = 1,
      from = min(dat$terminal_ploidy),
      to = max(dat$terminal_ploidy),
      n = 512,
      na.rm = TRUE
    )
    scaled_density <- density_fit$y / max(density_fit$y)
    center_x <- unique(dat$source_x)
    bind_rows(
      tibble(
        cohort = unique(dat$cohort),
        source = unique(dat$source),
        polygon_order = 1L,
        terminal_ploidy = density_fit$x[[1L]],
        violin_x = center_x
      ),
      tibble(
        cohort = unique(dat$cohort),
        source = unique(dat$source),
        polygon_order = seq_along(density_fit$x) + 1L,
        terminal_ploidy = density_fit$x,
        violin_x = center_x - half_violin_width * scaled_density
      ),
      tibble(
        cohort = unique(dat$cohort),
        source = unique(dat$source),
        polygon_order = length(density_fit$x) + 2L,
        terminal_ploidy = density_fit$x[[length(density_fit$x)]],
        violin_x = center_x
      )
    )
  }) %>%
  bind_rows() %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    source = factor(source, levels = c("Observed", "Predicted"))
  ) %>%
  arrange(cohort, source, polygon_order)

write.table(
  predicted,
  file.path(data_dir, "figure4a_cohort_model_dynamics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  observed,
  file.path(data_dir, "figure4a_observed_burden_mean_sd.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  terminal_plot,
  file.path(data_dir, "figure4a_terminal_ploidy_equal_tumor_weight.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

component_colors <- c(
  "Viable" = "#56B4E9",
  "Stress dead (model)" = "#D55E00",
  "CIN-associated dead" = "#009E73"
)
trajectory_colors <- c(
  "Effective O2 (%)" = "#E69F00",
  "Mean viable ploidy" = "#7B3294"
)
terminal_colors <- c(
  "Observed" = "#0072B2",
  "Predicted" = "#D55E00"
)

base_theme <- theme_bw(base_size = 14) +
  theme(
    plot.subtitle = element_text(size = 12.6, color = "#444444"),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12.6, color = "#222222"),
    axis.line = element_line(
      color = "#202020", linewidth = 0.60, lineend = "square"
    ),
    axis.ticks = element_line(color = "#202020", linewidth = 0.60),
    panel.border = element_rect(
      fill = NA, color = "#AEB6C0", linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(face = "bold", size = 12.8),
    legend.text = element_text(size = 11.7),
    legend.key.height = unit(0.75, "lines"),
    legend.key.width = unit(1.45, "lines"),
    legend.position = "none",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.justification = "left",
    plot.margin = margin(5, 9, 5, 8)
  )

make_column_strip <- function(label) {
  ggplot() +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = "#E6E6E6", color = "#B8B8B8", linewidth = 0.45
    ) +
    annotate(
      "text",
      x = 0.5, y = 0.5, label = label,
      fontface = "bold", size = 6.0, color = "#222222"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 3, 0))
}

make_row_strip <- function(label) {
  ggplot() +
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = "#E6E6E6", color = "#B8B8B8", linewidth = 0.45
    ) +
    annotate(
      "text",
      x = 0.5, y = 0.5, label = label,
      angle = -90, fontface = "bold", size = 6.0, color = "#222222"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 3))
}

make_axis_title <- function(label, angle = 0) {
  ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.5, label = label,
      angle = angle, fontface = "bold", size = 5.4, color = "#222222"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_dynamics_plot <- function(cohort_value, show_x_title) {
  cohort_factor <- factor(cohort_value, levels = cohort_levels)
  area_data <- predicted_long %>% filter(cohort == cohort_factor)
  line_data <- trajectory_lines %>% filter(cohort == cohort_factor)
  obs_data <- observed %>% filter(cohort == cohort_factor)
  sample_note <- if (identical(cohort_value, "2N")) {
    "Observed burden: n=4 at every displayed measurement day"
  } else {
    "Observed burden: n=4 through day 77; n=2 at days 81–91"
  }

  ggplot() +
    geom_area(
      data = area_data,
      aes(
        x = day,
        y = burden_mm3 / burden_to_shared,
        fill = component,
        group = component
      ),
      position = position_stack(reverse = TRUE),
      alpha = 0.75,
      linewidth = 0
    ) +
    geom_ribbon(
      data = obs_data,
      aes(
        x = day,
        ymin = observed_burden_lower_mm3 / burden_to_shared,
        ymax = observed_burden_upper_mm3 / burden_to_shared,
        group = observed_line_segment
      ),
      fill = "#666666",
      alpha = 0.18,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = obs_data,
      aes(
        x = day,
        y = observed_burden_mean_mm3 / burden_to_shared,
        group = observed_line_segment
      ),
      color = "black",
      linewidth = 0.9,
      lineend = "round"
    ) +
    geom_point(
      data = obs_data,
      aes(x = day, y = observed_burden_mean_mm3 / burden_to_shared),
      color = "black",
      fill = "white",
      shape = 21,
      stroke = 0.55,
      size = 1.65
    ) +
    {
      if (identical(cohort_value, "4N")) {
        list(
          geom_vline(
            xintercept = 79,
            color = "#8C510A",
            linewidth = 0.48,
            linetype = "22"
          ),
          annotate(
            "text", x = 79, y = 4.78,
            label = "n=4 to n=2\nobserved line break",
            hjust = 0.5, vjust = 1,
            size = 3.55, fontface = "bold", color = "#6B3D0C"
          )
        )
      }
    } +
    geom_line(
      data = line_data,
      aes(
        x = day,
        y = shared_axis_value,
        color = series,
        group = series
      ),
      linewidth = 1.0,
      lineend = "round"
    ) +
    scale_fill_manual(
      values = component_colors,
      name = "Model-predicted burden component",
      drop = FALSE
    ) +
    scale_color_manual(
      values = trajectory_colors,
      name = "Model-predicted trajectory",
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = c(0, 50, 100),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, shared_axis_max),
      breaks = c(0, 2, 4, 5),
      expand = expansion(mult = c(0, 0.012)),
      sec.axis = sec_axis(
        ~ . * burden_to_shared,
        breaks = c(0, 2000, 4000, 5000),
        labels = scales::label_number(big.mark = ",", accuracy = 1),
        name = NULL
      )
    ) +
    coord_cartesian(xlim = c(0, 100), expand = FALSE) +
    labs(
      subtitle = paste0(
        sample_note,
        "; gray ribbon = between-tumor sample SD (not error, CI, or model uncertainty)"
      ),
      x = NULL,
      y = NULL
    ) +
    guides(
      fill = guide_legend(order = 1, nrow = 1, byrow = TRUE),
      color = guide_legend(
        order = 2,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 1.45)
      )
    ) +
    base_theme +
    theme(
      axis.text.x = if (isTRUE(show_x_title)) {
        element_text(size = 12.6, color = "#222222")
      } else {
        element_blank()
      },
      axis.ticks.x = if (isTRUE(show_x_title)) {
        element_line(color = "#202020", linewidth = 0.60)
      } else {
        element_blank()
      },
      axis.line.x = if (isTRUE(show_x_title)) {
        element_line(
          color = "#202020", linewidth = 0.60, lineend = "square"
        )
      } else {
        element_blank()
      }
    )
}

make_terminal_plot <- function(cohort_value, show_x_title) {
  cohort_factor <- factor(cohort_value, levels = cohort_levels)
  dat <- terminal_plot %>% filter(cohort == cohort_factor)
  boxes <- terminal_box %>% filter(cohort == cohort_factor)

  ggplot(
    dat,
    aes(
      x = source_x,
      y = terminal_ploidy,
      weight = plot_weight,
      fill = source,
      group = source
    )
  ) +
    geom_polygon(
      data = terminal_half_violin %>% filter(cohort == cohort_factor),
      aes(
        x = violin_x,
        y = terminal_ploidy,
        fill = source,
        group = source
      ),
      inherit.aes = FALSE,
      color = "#4A4A4A",
      linewidth = 0.25,
      alpha = 0.92
    ) +
    geom_segment(
      data = boxes,
      aes(
        x = box_x, xend = box_x,
        y = ymin, yend = q1
      ),
      inherit.aes = FALSE,
      color = "#333333",
      linewidth = 0.28
    ) +
    geom_segment(
      data = boxes,
      aes(
        x = box_x, xend = box_x,
        y = q3, yend = ymax
      ),
      inherit.aes = FALSE,
      color = "#333333",
      linewidth = 0.28
    ) +
    geom_segment(
      data = boxes,
      aes(
        x = box_x - whisker_half_width,
        xend = box_x + whisker_half_width,
        y = ymin, yend = ymin
      ),
      inherit.aes = FALSE,
      color = "#333333",
      linewidth = 0.28
    ) +
    geom_segment(
      data = boxes,
      aes(
        x = box_x - whisker_half_width,
        xend = box_x + whisker_half_width,
        y = ymax, yend = ymax
      ),
      inherit.aes = FALSE,
      color = "#333333",
      linewidth = 0.28
    ) +
    geom_rect(
      data = boxes,
      aes(
        xmin = box_x - box_half_width,
        xmax = box_x + box_half_width,
        ymin = q1,
        ymax = q3,
        fill = source
      ),
      inherit.aes = FALSE,
      color = "#333333",
      linewidth = 0.28,
      alpha = 0.26,
      show.legend = FALSE
    ) +
    geom_segment(
      data = boxes,
      aes(
        x = box_x - box_half_width,
        xend = box_x + box_half_width,
        y = median, yend = median
      ),
      inherit.aes = FALSE,
      color = "#111111",
      linewidth = 0.35
    ) +
    scale_fill_manual(
      values = terminal_colors,
      name = "Terminal distribution",
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = c(1, 3, terminal_axis_display_max),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_x_continuous(
      limits = c(0.45, 2.45),
      breaks = c(1, 2),
      labels = c("Observed", "Predicted"),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      ylim = c(1, terminal_axis_display_max),
      expand = FALSE,
      clip = "on"
    ) +
    labs(
      subtitle = "Observed vs predicted; tumors weighted equally",
      x = NULL,
      y = NULL
    ) +
    guides(fill = guide_legend(order = 3, nrow = 1, byrow = TRUE)) +
    base_theme +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      axis.text.x = if (isTRUE(show_x_title)) {
        element_text(size = 12.6, color = "#222222")
      } else {
        element_blank()
      },
      axis.ticks.x = if (isTRUE(show_x_title)) {
        element_line(color = "#202020", linewidth = 0.60)
      } else {
        element_blank()
      },
      axis.line.x = if (isTRUE(show_x_title)) {
        element_line(
          color = "#202020", linewidth = 0.60, lineend = "square"
        )
      } else {
        element_blank()
      }
    )
}

p_2n_dynamics <- make_dynamics_plot("2N", show_x_title = FALSE)
p_4n_dynamics <- make_dynamics_plot("4N", show_x_title = TRUE)
p_2n_terminal <- make_terminal_plot("2N", show_x_title = FALSE)
p_4n_terminal <- make_terminal_plot("4N", show_x_title = TRUE)

column_strip_dynamics <- make_column_strip("Longitudinal dynamics")
column_strip_terminal <- make_column_strip("Terminal ploidy distribution")
row_strip_2n <- make_row_strip("2N")
row_strip_4n <- make_row_strip("4N")
shared_left_axis_title <- make_axis_title(
  "Effective O2 (%) / mean viable ploidy",
  angle = 90
)
shared_burden_axis_title <- make_axis_title(
  "Tumor burden (mm³)",
  angle = -90
)
shared_terminal_axis_title <- make_axis_title(
  "Terminal ploidy (N/22)",
  angle = 90
)
shared_day_axis_title <- make_axis_title("Day")

# Draw the four semantic legend groups as one horizontal row. A custom legend
# avoids patchwork stacking separately collected fill/color guides even when
# enough horizontal space is available.
legend_text_size <- 10.8 / ggplot2::.pt
legend_title_size <- 11.6 / ggplot2::.pt
legend_common <- list(
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off"),
  theme_void(),
  theme(plot.margin = margin(0, 4, 0, 4))
)

legend_observed <- ggplot() +
  annotate(
    "text", x = 0.00, y = 0.5,
    label = "Observed burden",
    hjust = 0, vjust = 0.5, fontface = "bold",
    size = legend_title_size, color = "#202020"
  ) +
  annotate(
    "segment", x = 0.335, xend = 0.425, y = 0.5, yend = 0.5,
    color = "black", linewidth = 0.9, lineend = "round"
  ) +
  annotate(
    "point", x = 0.38, y = 0.5,
    shape = 21, size = 2.3, stroke = 0.55,
    fill = "white", color = "black"
  ) +
  annotate(
    "text", x = 0.445, y = 0.5, label = "cohort mean",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  annotate(
    "rect", xmin = 0.690, xmax = 0.760, ymin = 0.30, ymax = 0.70,
    fill = "#666666", alpha = 0.18, color = "#666666", linewidth = 0.2
  ) +
  annotate(
    "text", x = 0.780, y = 0.5, label = "±1 between-tumor SD",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  legend_common

legend_burden <- ggplot() +
  annotate(
    "text", x = 0.00, y = 0.5,
    label = "Burden components",
    hjust = 0, vjust = 0.5, fontface = "bold",
    size = legend_title_size, color = "#202020"
  ) +
  annotate(
    "rect", xmin = 0.285, xmax = 0.320, ymin = 0.30, ymax = 0.70,
    fill = component_colors[["Viable"]], color = NA
  ) +
  annotate(
    "text", x = 0.330, y = 0.5, label = "Viable",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  annotate(
    "rect", xmin = 0.430, xmax = 0.465, ymin = 0.30, ymax = 0.70,
    fill = component_colors[["Stress dead (model)"]], color = NA
  ) +
  annotate(
    "text", x = 0.475, y = 0.5, label = "Stress dead (model)",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  annotate(
    "rect", xmin = 0.720, xmax = 0.755, ymin = 0.30, ymax = 0.70,
    fill = component_colors[["CIN-associated dead"]], color = NA
  ) +
  annotate(
    "text", x = 0.765, y = 0.5, label = "CIN-associated dead",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  legend_common

legend_trajectory <- ggplot() +
  annotate(
    "text", x = 0.00, y = 0.5,
    label = "Trajectories",
    hjust = 0, vjust = 0.5, fontface = "bold",
    size = legend_title_size, color = "#202020"
  ) +
  annotate(
    "segment", x = 0.260, xend = 0.320, y = 0.5, yend = 0.5,
    color = trajectory_colors[["Effective O2 (%)"]],
    linewidth = 1.35, lineend = "round"
  ) +
  annotate(
    "text", x = 0.335, y = 0.5, label = "Effective O2 (%)",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  annotate(
    "segment", x = 0.600, xend = 0.660, y = 0.5, yend = 0.5,
    color = trajectory_colors[["Mean viable ploidy"]],
    linewidth = 1.35, lineend = "round"
  ) +
  annotate(
    "text", x = 0.675, y = 0.5, label = "Mean viable ploidy",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  legend_common

legend_terminal <- ggplot() +
  annotate(
    "text", x = 0.00, y = 0.5, label = "Terminal",
    hjust = 0, vjust = 0.5, fontface = "bold",
    size = legend_title_size, color = "#202020"
  ) +
  annotate(
    "rect", xmin = 0.405, xmax = 0.465, ymin = 0.30, ymax = 0.70,
    fill = terminal_colors[["Observed"]], color = "#333333", linewidth = 0.25
  ) +
  annotate(
    "text", x = 0.485, y = 0.5, label = "Observed",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  annotate(
    "rect", xmin = 0.690, xmax = 0.750, ymin = 0.30, ymax = 0.70,
    fill = terminal_colors[["Predicted"]], color = "#333333", linewidth = 0.25
  ) +
  annotate(
    "text", x = 0.770, y = 0.5, label = "Predicted",
    hjust = 0, vjust = 0.5, size = legend_text_size, color = "#202020"
  ) +
  legend_common

single_row_legend <-
  legend_observed + legend_burden + legend_trajectory + legend_terminal +
  plot_layout(widths = c(1.15, 1.50, 1.25, 0.65))

# The seven dynamics columns and three terminal columns retain a strict 7:3
# ratio. Shared axis-title columns and the row-strip column are independent
# gutters and are excluded from that analytical-panel ratio. This lets every
# title span the relevant pair of cohort panels rather than attaching it to
# only the lower plot.
layout_design <- c(
  area(t = 1, l = 2, b = 1, r = 8),
  area(t = 1, l = 11, b = 1, r = 13),
  area(t = 2, l = 1, b = 3, r = 1),
  area(t = 2, l = 2, b = 2, r = 8),
  area(t = 2, l = 9, b = 3, r = 9),
  area(t = 2, l = 10, b = 3, r = 10),
  area(t = 2, l = 11, b = 2, r = 13),
  area(t = 2, l = 14, b = 2, r = 14),
  area(t = 3, l = 2, b = 3, r = 8),
  area(t = 3, l = 11, b = 3, r = 13),
  area(t = 3, l = 14, b = 3, r = 14),
  area(t = 4, l = 2, b = 4, r = 8),
  area(t = 5, l = 1, b = 5, r = 14),
  area(t = 6, l = 1, b = 6, r = 14)
)

column_strip_height <- 0.18
dynamics_row_height <- 1
axis_title_row_height <- 0.14
legend_gap_row_height <- 0.03
legend_row_height <- 0.28
axis_title_width <- 0.34
row_strip_width <- 0.35
output_width_in <- 20
# Match the non-label content box of the square final Figure 4 assembly:
# 5376 px wide by 1665 px high within the 5400 x 1780 px panel-A slot.
output_height_in <- output_width_in * 1665 / 5376

composite <-
  column_strip_dynamics +
  column_strip_terminal +
  shared_left_axis_title +
  p_2n_dynamics +
  shared_burden_axis_title +
  shared_terminal_axis_title +
  p_2n_terminal +
  row_strip_2n +
  p_4n_dynamics +
  p_4n_terminal +
  row_strip_4n +
  shared_day_axis_title +
  plot_spacer() +
  single_row_legend +
  plot_layout(
    design = layout_design,
    widths = c(
      axis_title_width,
      rep(1, 7),
      axis_title_width,
      axis_title_width,
      rep(1, 3),
      row_strip_width
    ),
    heights = c(
      column_strip_height,
      dynamics_row_height,
      dynamics_row_height,
      axis_title_row_height,
      legend_gap_row_height,
      legend_row_height
    )
  )

output_base <- file.path(figure_dir, "fig4a_combined_invivo_dynamics")
ggsave(
  paste0(output_base, ".png"),
  composite,
  width = output_width_in,
  height = output_height_in,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  paste0(output_base, ".pdf"),
  composite,
  width = output_width_in,
  height = output_height_in,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  paste0(output_base, ".svg"),
  composite,
  width = output_width_in,
  height = output_height_in,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

validation <- data.frame(
  metric = c(
    "n_tumors",
    "n_cohorts",
    "n_model_daily_rows",
    "n_observed_cohort_day_rows",
    "N_UNIT",
    "dynamics_x_min_day",
    "dynamics_x_max_day",
    "shared_axis_min",
    "shared_axis_max",
    "terminal_axis_display_max",
    "burden_axis_max_mm3",
    "max_observed_mean_plus_sd_mm3",
    "terminal_observed_rows",
    "terminal_predicted_state_rows",
    "terminal_equal_tumor_weighting",
    "terminal_distribution_geometry",
    "terminal_half_violin_side",
    "terminal_half_violin_linewidth",
    "terminal_box_whisker_linewidth",
    "terminal_box_median_linewidth",
    "terminal_raw_points_drawn",
    "legend_rows",
    "legend_semantic_groups",
    "legend_is_custom_fixed_layout",
    "terminal_observed_color",
    "terminal_predicted_color",
    "left_to_right_width_ratio",
    "shared_axis_titles_use_independent_layout_grobs",
    "shared_day_title_scope",
    "shared_left_title_scope",
    "shared_burden_title_scope",
    "shared_terminal_title_scope",
    "plot_native_axis_titles_removed",
    "row_strip_width_relative_to_main_column",
    "column_strip_longitudinal_dynamics",
    "column_strip_terminal_ploidy_distribution",
    "row_strip_top",
    "row_strip_bottom",
    "individual_titles_include_cohort_labels",
    "output_width_in",
    "output_height_in",
    "output_width_to_height_ratio",
    "black_model_total_burden_line_drawn",
    "black_observed_burden_line_drawn",
    "gray_observed_sd_ribbon_drawn"
    ,"observed_mean_line_broken_at_4N_composition_change"
    ,"figure4a_package_internal_measurements_verified"
    ,"fourN_day77_observed_n"
    ,"fourN_day81_observed_n"
    ,"gray_ribbon_is_between_tumor_sample_sd"
    ,"gray_ribbon_is_measurement_error_or_model_uncertainty"
  ),
  value = c(
    length(unique(trajectory$harvest)),
    length(unique(trajectory$cohort)),
    nrow(predicted),
    nrow(observed),
    n_unit,
    0,
    100,
    0,
    shared_axis_max,
    terminal_axis_display_max,
    burden_axis_max,
    max(observed$observed_burden_upper_mm3, na.rm = TRUE),
    sum(terminal_plot$source == "Observed"),
    sum(terminal_plot$source == "Predicted"),
    "TRUE",
    "left half violin plus right weighted box-and-whisker",
    "left",
    0.25,
    0.28,
    0.35,
    "FALSE",
    1,
    4,
    "TRUE",
    unname(terminal_colors[["Observed"]]),
    unname(terminal_colors[["Predicted"]]),
    "7:3",
    "TRUE",
    "longitudinal dynamics column only",
    "both cohort rows",
    "both longitudinal cohort rows",
    "both terminal cohort rows",
    "TRUE",
    row_strip_width,
    "TRUE",
    "TRUE",
    "2N",
    "4N",
    "FALSE",
    output_width_in,
    output_height_in,
    output_width_in / output_height_in,
    "FALSE",
    "TRUE",
    "TRUE",
    "TRUE",
    all(as.logical(burden_audit$exact_package_internal_match)),
    4,
    2,
    "TRUE",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(data_dir, "figure4a_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Figure 4A written: ", paste0(output_base, ".{png,pdf,svg}"))

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Figure4 <- function() {
  data_dir <- file.path(DATA_ROOT, "Figure4")
  panel_dir <- file.path(data_dir, "panels")
  required <- file.path(data_dir, c(
    "figure4a_seed25_burden_timecourse.tsv",
    "figure4a_seed25_weighted_mean_chromosome_timecourse.tsv",
    "figure4a_seed25_terminal_chromosome_observed_vs_predicted.tsv",
    "invivo_best_parameters_500seeds.tsv",
    "fixed_o2_dominant_ploidy_201grid.tsv",
    "parameter_function_groups.tsv",
    "parameter_function_group_palette.tsv",
    "invivo_best_tsne_cluster_coordinates.tsv",
    "invivo_tsne_cluster_summary.tsv",
    "invivo_tsne_run_metadata.tsv",
    "invivo_parameter_table_seed25.csv",
    "invivo_fit_objective_ranking_500seeds.tsv"
  ))
  require_files(required, "Figure 4 intermediate")
  run_figure4_continuous_association(data_dir = data_dir)
  continuous_required <- file.path(data_dir, c(
    "figure4a_burden_measurement_audit.tsv",
    "continuous_ploidy_spearman_by_o2.tsv",
    "continuous_ploidy_parameter_ranking.tsv",
    "all_parameter_fitted_endpoint_values.tsv",
    "all_parameter_pooled_distribution_summary.tsv",
    "all_parameter_log10_range_summary.tsv",
    "all_parameter_log10_density.tsv",
    "parameter_display_labels.tsv"
  ))
  require_files(continuous_required, "Figure 4 continuous-ploidy intermediate")
  env <- c(
    paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
    paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
    paste0("FIGURE4_DATA_DIR=", data_dir),
    paste0("FIGURE4_PANEL_DIR=", panel_dir)
  )
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_Figure4.R"), mustWork = TRUE),
    env = c("FIGURE4_DRAW_WORKER=1", env)
  )
  render_parameter_landscape(
    data_dir = data_dir,
    panel_dir = panel_dir,
    output_dir = OUTPUT_ROOT
  )
  require_files(
    c(
      file.path(data_dir, c(
        "exploratory_cluster_parameter_omnibus_tests.tsv",
        "figure4d_top6_parameters.tsv"
      )),
      file.path(panel_dir, c(
        "parameter_continuous_ploidy_landscape.png",
        "parameter_tsne_groups.png",
        "top6_cluster_parameter_distributions.png"
      ))
    ),
    "Figure 4B-D parameter-landscape output"
  )
  compose_top_with_left_and_stacked_right_panels(
    data_dir = data_dir,
    panel_dir = panel_dir,
    output_dir = OUTPUT_ROOT,
    top_panel = file.path(
      panel_dir, "fig4a_combined_invivo_dynamics.png"
    ),
    bottom_left_panel = file.path(
      panel_dir, "parameter_continuous_ploidy_landscape.png"
    ),
    bottom_right_top_panel = file.path(
      panel_dir, "parameter_tsne_groups.png"
    ),
    bottom_right_bottom_panel = file.path(
      panel_dir, "top6_cluster_parameter_distributions.png"
    ),
    supplementary_source = file.path(
      panel_dir, "all_parameter_fitted_endpoint_distributions"
    ),
    supplementary_destination = file.path(
      OUTPUT_ROOT, "supp_fig4-1_all18_cluster_prior_violins"
    ),
    output_basename = "assembled_fig4",
    validation_basename = "figure4_layout_validation.tsv"
  )
  require_files(
    file.path(OUTPUT_ROOT, c(
      "assembled_fig4.png", "assembled_fig4.pdf",
      "supp_fig4-1_all18_cluster_prior_violins.png",
      "supp_fig4-1_all18_cluster_prior_violins.pdf"
    )),
    "Figure 4/Supplementary Figure 4-1 output"
  )
  invisible(file.path(OUTPUT_ROOT, "assembled_fig4.png"))
}

if (sys.nframe() == 0L) draw_Figure4()

}
