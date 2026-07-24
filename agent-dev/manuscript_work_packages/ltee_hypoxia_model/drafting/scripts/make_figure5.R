#!/usr/bin/env Rscript

# Figure 5 drafting generator.
#
# Scope:
#   A. Direct in-sample observed-predicted summaries for both fitted contexts.
#   B. Full parameter-ratio view across the six approved July joint-fit winners.
#   C. Proliferation functions across the same six winners.
#   D. Stress-linked missegregation functions across the same six winners.
#   E. Post-missegregation survival functions across the same six winners.
#
# This script reads existing fit outputs only. It does not refit the model,
# modify production analysis code, or treat across-winner variation as a
# confidence interval.

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("ggplot2", "dplyr", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

repo_root <- Sys.getenv("MININGCLONEID_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) {
  repo_root <- normalizePath(".", mustWork = TRUE)
}

draft_root <- file.path(
  repo_root,
  "agent-dev",
  "manuscript_work_packages",
  "ltee_hypoxia_model",
  "drafting"
)
frozen_root <- file.path(draft_root, "source_tables", "frozen_inputs", "F5")
selection_path <- file.path(frozen_root, "selected_results.tsv")

if (!file.exists(selection_path)) {
  stop("Missing Figure 5 selection table: ", selection_path)
}

initial_root <- file.path(draft_root, "initial_subpanels")
refined_root <- file.path(draft_root, "refined_subpanels")
recommended_root <- file.path(draft_root, "final_figures", "recommended")

allowed_output_dirs <- c(
  file.path(initial_root, paste0("F5", LETTERS[1:5])),
  refined_root,
  recommended_root
)
invisible(lapply(
  allowed_output_dirs,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

selected <- read.delim(
  selection_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
selected <- selected[selected$record_type == "joint_pair_best", , drop = FALSE]

required_selection_fields <- c(
  "warmup_label",
  "invivo_seed",
  "invitro_seed",
  "selected_seed",
  "objective",
  "bundle_dir"
)
if (!all(required_selection_fields %in% names(selected))) {
  stop(
    "Selection table lacks required fields: ",
    paste(setdiff(required_selection_fields, names(selected)), collapse = ", ")
  )
}
if (nrow(selected) != 6L) {
  stop("Figure 5 requires exactly six approved joint-pair winners; found ", nrow(selected))
}
if (length(unique(selected$invitro_seed)) != 1L ||
    unique(selected$invitro_seed) != "seed10") {
  stop("The approved Figure 5 universe must use in-vitro anchor seed10")
}

extract_region <- function(warmup_label) {
  hit <- regmatches(warmup_label, regexpr("C[0-9]{2}", warmup_label))
  ifelse(nzchar(hit), hit, "unclassified")
}

selected$region <- extract_region(selected$warmup_label)
selected$winner_label <- paste0(
  selected$region,
  " / ",
  selected$invivo_seed,
  " / ",
  selected$selected_seed
)
selected$fit_dir <- file.path(
  frozen_root, "winners", selected$warmup_label
)

missing_fit_dirs <- selected$fit_dir[!dir.exists(selected$fit_dir)]
if (length(missing_fit_dirs)) {
  stop("Missing selected fit directories:\n", paste(missing_fit_dirs, collapse = "\n"))
}

read_winner_table <- function(fit_dir, relative_path) {
  path <- file.path(fit_dir, relative_path)
  if (!file.exists(path)) {
    stop("Missing selected-winner source table: ", path)
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

context_colors <- c(
  "In vivo" = "#D55E00",
  "In vitro" = "#0072B2"
)
region_colors <- c(
  "C01" = "#009E73",
  "C02" = "#0072B2",
  "C03" = "#CC79A7"
)
in_vivo_color <- unname(context_colors[["In vivo"]])
in_vitro_color <- unname(context_colors[["In vitro"]])

theme_manuscript <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E6E8EB", linewidth = 0.28),
      panel.border = element_rect(color = "#4B5563", linewidth = 0.45),
      axis.text = element_text(color = "#27313A"),
      axis.title = element_text(color = "#1F2933"),
      strip.background = element_rect(fill = "#F3F4F6", color = "#9CA3AF"),
      strip.text = element_text(face = "bold", color = "#27313A"),
      legend.key.height = grid::unit(0.28, "cm"),
      legend.key.width = grid::unit(0.46, "cm"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", color = "#111827"),
      plot.subtitle = element_text(color = "#4B5563"),
      plot.caption = element_text(color = "#5F6B76", hjust = 0),
      plot.margin = margin(3, 4, 3, 4)
    )
}

save_plot_pair <- function(plot, stub, width, height, dpi = 220) {
  pdf_path <- paste0(stub, ".pdf")
  png_path <- paste0(stub, ".png")
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  invisible(c(pdf = pdf_path, png = png_path))
}

# -------------------------------------------------------------------------
# Panel A: compact, direct in-sample fit adequacy across all six winners.
# -------------------------------------------------------------------------

burden_rows <- list()
terminal_rows <- list()
growth_rows <- list()
kary_rows <- list()

for (i in seq_len(nrow(selected))) {
  fit_dir <- selected$fit_dir[[i]]
  winner <- selected$winner_label[[i]]
  region <- selected$region[[i]]

  burden <- read_winner_table(fit_dir, "invivo_burden_fit.tsv")
  burden <- burden[
    is.finite(burden$obs_norm) &
      is.finite(burden$pred_norm) &
      burden$day > 0,
    ,
    drop = FALSE
  ]
  burden$winner <- winner
  burden$region <- region
  burden$obs_id <- paste(burden$harvest, burden$day, sep = "|")
  burden_rows[[i]] <- burden[, c(
    "winner", "region", "obs_id", "obs_norm", "pred_norm"
  )]

  terminal <- read_winner_table(fit_dir, "invivo_terminal_ploidy_fit.tsv")
  terminal_split <- split(terminal, terminal$harvest)
  terminal_summary <- do.call(rbind, lapply(terminal_split, function(x) {
    observed_total <- sum(x$obs_count, na.rm = TRUE)
    predicted_total <- sum(x$pred_fraction, na.rm = TRUE)
    data.frame(
      obs_id = x$harvest[[1]],
      cohort = x$cohort[[1]],
      observed_mean_N = if (observed_total > 0) {
        sum(x$N * x$obs_count, na.rm = TRUE) / observed_total
      } else {
        NA_real_
      },
      predicted_mean_N = if (predicted_total > 0) {
        sum(x$N * x$pred_fraction, na.rm = TRUE) / predicted_total
      } else {
        NA_real_
      }
    )
  }))
  terminal_summary$winner <- winner
  terminal_summary$region <- region
  terminal_rows[[i]] <- terminal_summary

  growth <- read_winner_table(fit_dir, "invitro_growth_loglik.tsv")
  growth <- growth[
    is.finite(growth$observed_growth) &
      is.finite(growth$predicted_growth_rate),
    ,
    drop = FALSE
  ]
  growth$winner <- winner
  growth$region <- region
  growth$obs_id <- paste(
    growth$segment_id,
    growth$passage_id,
    growth$lineage_terminal_key,
    growth$lineage_passage_index,
    sep = "|"
  )
  growth_rows[[i]] <- growth[, c(
    "winner",
    "region",
    "obs_id",
    "observed_growth",
    "predicted_growth_rate"
  )]

  kary <- read_winner_table(fit_dir, "invitro_lineage_summary.tsv")
  kary <- kary[
    is.finite(kary$observed_mean_kary_N) &
      is.finite(kary$predicted_mean_kary_N),
    ,
    drop = FALSE
  ]
  kary$winner <- winner
  kary$region <- region
  kary$obs_id <- paste(
    kary$segment_id,
    kary$passage_id,
    kary$lineage_terminal_key,
    kary$lineage_passage_index,
    sep = "|"
  )
  kary_rows[[i]] <- kary[, c(
    "winner",
    "region",
    "obs_id",
    "observed_mean_kary_N",
    "predicted_mean_kary_N"
  )]
}

burden_fit <- bind_rows(burden_rows)
terminal_fit <- bind_rows(terminal_rows)
growth_fit <- bind_rows(growth_rows)
kary_fit <- bind_rows(kary_rows)

assert_six_per_observation <- function(df, label) {
  counts <- table(df$obs_id)
  if (length(counts) == 0L || any(counts != 6L)) {
    stop(label, " does not contain exactly six winner predictions per observation")
  }
}
assert_six_per_observation(burden_fit, "In-vivo burden fit")
assert_six_per_observation(terminal_fit, "In-vivo terminal-ploidy fit")
assert_six_per_observation(growth_fit, "In-vitro growth fit")
assert_six_per_observation(kary_fit, "In-vitro karyotype fit")

make_fit_scatter <- function(
  df,
  observed_col,
  predicted_col,
  title,
  x_label,
  y_label,
  color,
  limits,
  show_range = FALSE
) {
  summary_df <- df %>%
    group_by(obs_id) %>%
    summarise(
      observed = first(.data[[observed_col]]),
      predicted = median(.data[[predicted_col]], na.rm = TRUE),
      predicted_min = min(.data[[predicted_col]], na.rm = TRUE),
      predicted_max = max(.data[[predicted_col]], na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot() +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "#6B7280",
      linetype = "22",
      linewidth = 0.5
    )
  if (show_range) {
    p <- p + geom_segment(
      data = summary_df,
      aes(
        x = observed,
        xend = observed,
        y = predicted_min,
        yend = predicted_max
      ),
      color = alpha(color, 0.48),
      linewidth = 0.55
    )
  }
  p +
    geom_point(
      data = df,
      aes(x = .data[[observed_col]], y = .data[[predicted_col]]),
      color = alpha(color, 0.18),
      size = 0.9
    ) +
    geom_point(
      data = summary_df,
      aes(x = observed, y = predicted),
      shape = 21,
      fill = color,
      color = "white",
      stroke = 0.25,
      size = 1.6
    ) +
    coord_cartesian(xlim = limits, ylim = limits, expand = FALSE) +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme_manuscript(base_size = 7) +
    theme(
      plot.title = element_text(size = 7.5),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      plot.margin = margin(2, 3, 2, 3)
    )
}

p_a_burden <- make_fit_scatter(
  burden_fit,
  "obs_norm",
  "pred_norm",
  "In vivo | normalized burden",
  "Observed",
  "Predicted",
  in_vivo_color,
  c(0, 1),
  show_range = FALSE
)
p_a_terminal <- make_fit_scatter(
  terminal_fit,
  "observed_mean_N",
  "predicted_mean_N",
  "In vivo | terminal mean N",
  "Observed mean N",
  "Predicted mean N",
  in_vivo_color,
  c(38, 57),
  show_range = TRUE
)
p_a_growth <- make_fit_scatter(
  growth_fit,
  "observed_growth",
  "predicted_growth_rate",
  "In vitro | growth rate",
  "Observed rate",
  "Predicted rate",
  in_vitro_color,
  c(0.2, 2.2),
  show_range = FALSE
)
p_a_kary <- make_fit_scatter(
  kary_fit,
  "observed_mean_kary_N",
  "predicted_mean_kary_N",
  "In vitro | mean N",
  "Observed mean N",
  "Predicted mean N",
  in_vitro_color,
  c(44, 92),
  show_range = TRUE
)

p_a <- wrap_plots(
  p_a_burden,
  p_a_terminal,
  p_a_growth,
  p_a_kary,
  ncol = 2,
  nrow = 2
) +
  plot_annotation(
    title = "A  In-sample fit adequacy",
    subtitle = paste0(
      "Pale: six selected fits; filled: median. ",
      "Bars: selected-fit range where legible."
    ),
    theme = theme(
      plot.title = element_text(size = 9, face = "bold", color = "#111827"),
      plot.subtitle = element_text(size = 7, color = "#4B5563"),
      plot.margin = margin(1, 3, 1, 3)
    )
  )

# -------------------------------------------------------------------------
# Panel B: all-six soft-coupled parameter ratios on one log2 scale.
# -------------------------------------------------------------------------

ratio_rows <- list()
for (i in seq_len(nrow(selected))) {
  ratio <- read_winner_table(selected$fit_dir[[i]], "joint_soft_coupling.tsv")
  required_ratio_fields <- c(
    "parameter",
    "vivo_natural",
    "vitro_natural",
    "ratio_vivo_to_vitro",
    "feasible_at_solution",
    "projection_applied"
  )
  if (!all(required_ratio_fields %in% names(ratio))) {
    stop("joint_soft_coupling.tsv lacks required fields for ", selected$warmup_label[[i]])
  }
  ratio <- ratio[
    is.finite(ratio$ratio_vivo_to_vitro) &
      ratio$ratio_vivo_to_vitro > 0,
    ,
    drop = FALSE
  ]
  ratio$winner <- selected$winner_label[[i]]
  ratio$region <- selected$region[[i]]
  ratio$log2_ratio <- log2(ratio$ratio_vivo_to_vitro)
  ratio_rows[[i]] <- ratio
}
ratio_df <- bind_rows(ratio_rows)

if (nrow(ratio_df) != 84L ||
    length(unique(ratio_df$parameter)) != 14L ||
    any(table(ratio_df$parameter) != 6L)) {
  stop("Panel B requires 14 coupled parameters with six values each")
}
if (any(!ratio_df$feasible_at_solution) || any(ratio_df$projection_applied)) {
  stop("Unexpected infeasible or projected soft-coupling solution in Figure 5 universe")
}

parameter_groups <- data.frame(
  parameter = c(
    "lam_max", "p_mis_base", "p_wgd",
    "p_misseg", "k_o_mis",
    "O2_crit", "n_O",
    "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
    "buffer_smax", "buffer_beta", "buffer_n_exp"
  ),
  mechanism = c(
    rep("Division +\ngenome change", 3),
    rep("Stress-linked\nmissegregation", 2),
    rep("Oxygen\nresponse", 2),
    rep("Ploidy-linked\ngrowth + death", 4),
    rep("Post-MS\nsurvival", 3)
  ),
  parameter_order = seq_len(14),
  stringsAsFactors = FALSE
)
mechanism_levels <- unique(parameter_groups$mechanism)

ratio_df <- ratio_df %>%
  left_join(parameter_groups, by = "parameter")
if (any(is.na(ratio_df$mechanism))) {
  stop("Unmapped soft-coupling parameters in Panel B")
}
ratio_df$mechanism <- factor(ratio_df$mechanism, levels = mechanism_levels)
ratio_df$parameter <- factor(
  ratio_df$parameter,
  levels = rev(parameter_groups$parameter)
)

ratio_summary <- ratio_df %>%
  group_by(mechanism, parameter) %>%
  summarise(
    ratio_min = min(log2_ratio),
    ratio_median = median(log2_ratio),
    ratio_max = max(log2_ratio),
    .groups = "drop"
  )

ratio_breaks <- seq(-6, 12, by = 3)
ratio_break_labels <- vapply(ratio_breaks, function(x) {
  if (x == 0) {
    "1x"
  } else if (x < 0) {
    paste0("1/", format(2^abs(x), scientific = FALSE, trim = TRUE), "x")
  } else {
    paste0(format(2^x, scientific = FALSE, trim = TRUE), "x")
  }
}, character(1))

p_b <- ggplot(ratio_df, aes(y = parameter)) +
  geom_vline(
    xintercept = 0,
    color = "#111827",
    linetype = "22",
    linewidth = 0.55
  ) +
  geom_segment(
    data = ratio_summary,
    aes(
      x = ratio_min,
      xend = ratio_max,
      y = parameter,
      yend = parameter
    ),
    inherit.aes = FALSE,
    color = "#7A8490",
    linewidth = 1.25,
    lineend = "round"
  ) +
  geom_point(
    aes(x = log2_ratio, color = region),
    position = position_jitter(width = 0, height = 0.105, seed = 20260724),
    size = 1.8,
    alpha = 0.84
  ) +
  geom_point(
    data = ratio_summary,
    aes(x = ratio_median, y = parameter),
    inherit.aes = FALSE,
    shape = 23,
    fill = "#111827",
    color = "white",
    stroke = 0.35,
    size = 2.65
  ) +
  facet_grid(
    rows = vars(mechanism),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_color_manual(
    values = region_colors,
    breaks = c("C01", "C02", "C03"),
    name = "Warm-start region"
  ) +
  scale_x_continuous(
    breaks = ratio_breaks,
    labels = ratio_break_labels,
    limits = c(-7, 13),
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  labs(
    title = "B  In-vivo / in-vitro parameter ratios",
    subtitle = paste0(
      "Six values/parameter; gray = range; black diamond = median.\n",
      "Common log2 ratio scale for all 14 coupled parameters."
    ),
    x = "in vivo / in vitro ratio  <- higher in vitro      higher in vivo ->",
    y = NULL
  ) +
  theme_manuscript(base_size = 7) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "#F4F5F7", color = "#C1C7CD"),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 7,
      color = "#374151"
    ),
    axis.text.y = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    panel.spacing.y = grid::unit(0.045, "cm"),
    legend.position = "top",
    legend.justification = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 9),
    plot.subtitle = element_text(size = 7),
    plot.margin = margin(2, 3, 2, 3)
  )

# -------------------------------------------------------------------------
# Panels C-D: common-grid oxygen response functions.
# -------------------------------------------------------------------------

interpolate_winner_function <- function(
  fit_dir,
  context_dir,
  winner,
  region,
  value_col,
  common_grid
) {
  table_path <- file.path(
    "viz",
    context_dir,
    "functional_curve_oxygen_multi_ploidy.tsv"
  )
  curve <- read_winner_table(fit_dir, table_path)
  required_fields <- c("oxygen_pct", "cohort", value_col)
  if (!all(required_fields %in% names(curve))) {
    stop(table_path, " lacks required fields: ", paste(required_fields, collapse = ", "))
  }
  curve <- curve[
    curve$cohort %in% c("2N", "4N") &
      is.finite(curve$oxygen_pct) &
      is.finite(curve[[value_col]]),
    ,
    drop = FALSE
  ]
  context_label <- if (context_dir == "invivo") "In vivo" else "In vitro"
  out <- lapply(c("2N", "4N"), function(cohort_value) {
    cohort_curve <- curve[curve$cohort == cohort_value, , drop = FALSE]
    if (!nrow(cohort_curve)) {
      stop("Missing ", cohort_value, " curve for ", winner, " / ", context_label)
    }
    interpolated <- approx(
      x = cohort_curve$oxygen_pct,
      y = cohort_curve[[value_col]],
      xout = common_grid,
      ties = mean,
      rule = 2
    )
    data.frame(
      winner = winner,
      region = region,
      context = context_label,
      cohort = cohort_value,
      oxygen_pct = interpolated$x,
      value = interpolated$y
    )
  })
  bind_rows(out)
}

build_common_grid_function <- function(value_col, common_grid = seq(0, 5, by = 0.025)) {
  rows <- list()
  index <- 0L
  for (i in seq_len(nrow(selected))) {
    for (context_dir in c("invivo", "invitro")) {
      index <- index + 1L
      rows[[index]] <- interpolate_winner_function(
        fit_dir = selected$fit_dir[[i]],
        context_dir = context_dir,
        winner = selected$winner_label[[i]],
        region = selected$region[[i]],
        value_col = value_col,
        common_grid = common_grid
      )
    }
  }
  bind_rows(rows)
}

summarise_function <- function(df, x_col) {
  df %>%
    group_by(context, cohort, .data[[x_col]]) %>%
    summarise(
      value = median(value, na.rm = TRUE),
      value_min = min(value, na.rm = TRUE),
      value_max = max(value, na.rm = TRUE),
      n_solutions = n(),
      .groups = "drop"
    )
}

validate_six_function_solutions <- function(summary_df, label) {
  if (any(summary_df$n_solutions != 6L)) {
    stop(label, " does not contain six selected solutions at every grid point")
  }
}

proliferation_df <- build_common_grid_function("proliferation_rate")
missegregation_df <- build_common_grid_function("ms_rate")
proliferation_summary <- summarise_function(proliferation_df, "oxygen_pct")
missegregation_summary <- summarise_function(missegregation_df, "oxygen_pct")
validate_six_function_solutions(proliferation_summary, "Proliferation panel")
validate_six_function_solutions(missegregation_summary, "Missegregation panel")

make_o2_function_plot <- function(
  raw_df,
  summary_df,
  title,
  y_label
) {
  ggplot(
    raw_df,
    aes(
      x = oxygen_pct,
      y = value,
      color = context,
      linetype = context,
      group = interaction(context, winner)
    )
  ) +
    geom_line(
      linewidth = 0.38,
      alpha = 0.22,
      show.legend = FALSE
    ) +
    geom_line(
      data = summary_df,
      aes(
        x = oxygen_pct,
        y = value,
        color = context,
        linetype = context,
        group = context
      ),
      inherit.aes = FALSE,
      linewidth = 1.12
    ) +
    facet_wrap(vars(cohort), nrow = 1) +
    scale_color_manual(
      values = context_colors,
      breaks = c("In vivo", "In vitro"),
      name = "Context"
    ) +
    scale_linetype_manual(
      values = c("In vivo" = "solid", "In vitro" = "31"),
      breaks = c("In vivo", "In vitro"),
      name = "Context"
    ) +
    scale_x_continuous(
      breaks = c(0, 5),
      labels = c("0", "5"),
      limits = c(0, 5),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title = title,
      subtitle = NULL,
      x = "Oxygen (%)",
      y = y_label
    ) +
    theme_manuscript(base_size = 7) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 8.2),
      strip.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      panel.spacing.x = grid::unit(0.24, "cm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

p_c <- make_o2_function_plot(
  proliferation_df,
  proliferation_summary,
  "C  Proliferation",
  "Rate (day^-1)"
)
p_d <- make_o2_function_plot(
  missegregation_df,
  missegregation_summary,
  "D  Stress-linked\nmissegregation",
  "MS rate / chromosome"
)

# -------------------------------------------------------------------------
# Panel E: post-missegregation survival across the identical six winners.
# -------------------------------------------------------------------------

survival_rows <- list()
index <- 0L
for (i in seq_len(nrow(selected))) {
  for (context_dir in c("invivo", "invitro")) {
    index <- index + 1L
    survival <- read_winner_table(
      selected$fit_dir[[i]],
      file.path("viz", context_dir, "functional_curve_ploidy.tsv")
    )
    required_survival_fields <- c("N", "endpoint_value", "viability_after_ms")
    if (!all(required_survival_fields %in% names(survival))) {
      stop("functional_curve_ploidy.tsv lacks required survival fields")
    }
    survival <- survival[
      is.finite(survival$endpoint_value) &
        is.finite(survival$viability_after_ms),
      ,
      drop = FALSE
    ]
    survival$winner <- selected$winner_label[[i]]
    survival$region <- selected$region[[i]]
    survival$context <- if (context_dir == "invivo") "In vivo" else "In vitro"
    survival$cohort <- "All N"
    survival$value <- survival$viability_after_ms
    survival_rows[[index]] <- survival[, c(
      "winner",
      "region",
      "context",
      "cohort",
      "endpoint_value",
      "value"
    )]
  }
}
survival_df <- bind_rows(survival_rows)
survival_summary <- summarise_function(survival_df, "endpoint_value")
validate_six_function_solutions(survival_summary, "Post-MS survival panel")

p_e <- ggplot(
  survival_df,
  aes(
    x = endpoint_value,
    y = value,
    color = context,
    linetype = context,
    group = interaction(context, winner)
  )
) +
  geom_line(
    linewidth = 0.38,
    alpha = 0.22,
    show.legend = FALSE
  ) +
  geom_line(
    data = survival_summary,
    aes(
      x = endpoint_value,
      y = value,
      color = context,
      linetype = context,
      group = context
    ),
    inherit.aes = FALSE,
    linewidth = 1.12
  ) +
  scale_color_manual(
    values = context_colors,
    breaks = c("In vivo", "In vitro"),
    name = "Context"
  ) +
  scale_linetype_manual(
    values = c("In vivo" = "solid", "In vitro" = "31"),
    breaks = c("In vivo", "In vitro"),
    name = "Context"
  ) +
  scale_x_continuous(
    breaks = c(22, 88, 154),
    limits = c(22, 154),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "E  Post-missegregation\nsurvival",
    subtitle = NULL,
    x = "Chromosome number (N)",
    y = "Survival"
  ) +
  theme_manuscript(base_size = 7) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 8.2),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    plot.margin = margin(2, 2, 2, 2)
  )

# -------------------------------------------------------------------------
# Save individual draft panels and the assembled Figure 5.
# -------------------------------------------------------------------------

save_plot_pair(
  p_a,
  file.path(initial_root, "F5A", "figure5A_fit_adequacy"),
  width = 7.1,
  height = 4.4
)
save_plot_pair(
  p_b,
  file.path(initial_root, "F5B", "figure5B_parameter_ratios"),
  width = 7.1,
  height = 3.7
)
save_plot_pair(
  p_c,
  file.path(initial_root, "F5C", "figure5C_proliferation"),
  width = 3.25,
  height = 2.5
)
save_plot_pair(
  p_d,
  file.path(initial_root, "F5D", "figure5D_missegregation"),
  width = 3.25,
  height = 2.5
)
save_plot_pair(
  p_e,
  file.path(initial_root, "F5E", "figure5E_survival"),
  width = 3.25,
  height = 2.5
)

bottom_row <- (
  p_c +
    p_d +
    p_e +
    plot_layout(ncol = 3, widths = c(1.04, 1.04, 0.92), guides = "collect")
) &
  theme(legend.position = "bottom")

panel_a_for_assembly <- wrap_elements(plot = p_a)

figure5 <- (
  panel_a_for_assembly /
    p_b /
    bottom_row +
    plot_layout(heights = c(1.18, 1.25, 0.70))
) +
  plot_annotation(
    title = "Joint-fit adequacy and context-specific fitted functions",
    subtitle = paste0(
      "Six approved July joint-pair winners. Direct fits are in-sample; ",
      "across-winner ranges and traces show solution sensitivity, not confidence intervals."
    ),
    caption = paste0(
      "C-E: thin curves = six fits; thick = pointwise median. ",
      "Necrosis omitted: predictions unavailable.\n",
      "B-E use identical six winners; C-D are display-interpolated ",
      "to a declared 0-5% O2 grid."
    ),
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        size = 11.5,
        face = "bold",
        color = "#111827",
        margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        size = 7.3,
        color = "#4B5563",
        margin = margin(b = 2)
      ),
      plot.caption = element_text(
        size = 7,
        color = "#5F6B76",
        hjust = 0,
        margin = margin(t = 2)
      ),
      plot.margin = margin(5, 5, 4, 5)
    )
  )

refined_stub <- file.path(
  refined_root,
  "figure5_joint_fit_adequacy_and_context_functions"
)
recommended_stub <- file.path(
  recommended_root,
  "figure5_joint_fit_adequacy_and_context_functions"
)

save_plot_pair(
  figure5,
  refined_stub,
  width = 7.1,
  height = 9,
  dpi = 220
)
save_plot_pair(
  figure5,
  recommended_stub,
  width = 7.1,
  height = 9,
  dpi = 300
)

cat("Figure 5 generation complete.\n")
cat("Selected winners:", nrow(selected), "\n")
cat(
  "Total objective range:",
  sprintf("%.4f", min(selected$objective)),
  "to",
  sprintf("%.4f", max(selected$objective)),
  "\n"
)
cat("Panel B ratio rows:", nrow(ratio_df), "\n")
cat("Recommended PNG:", paste0(recommended_stub, ".png"), "\n")
cat("Recommended PDF:", paste0(recommended_stub, ".pdf"), "\n")
