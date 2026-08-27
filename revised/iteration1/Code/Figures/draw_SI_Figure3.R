#!/usr/bin/env Rscript

if (identical(Sys.getenv("SI_FIGURE3_DRAW_WORKER"), "1")) {

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("dplyr", "ggplot2", "patchwork", "scales", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
})
workspace_root <- normalizePath(
  Sys.getenv(
    "FIGURE_WORKSPACE_ROOT",
    unset = file.path(script_dir, "..", "..")
  ),
  mustWork = TRUE
)
repo_root <- normalizePath(
  Sys.getenv(
    "HYPOXIA_REPO_ROOT",
    unset = file.path(script_dir, "..", "..", "..", "..", "..", "..")
  ),
  mustWork = TRUE
)
frozen_root <- normalizePath(
  Sys.getenv(
    "SI_FIGURE3_DATA_DIR",
    unset = file.path(
      workspace_root, "data", "Figures", "SI_Figure3"
    )
  ),
  mustWork = TRUE
)
diagnostic_root <- frozen_root
out_dir <- normalizePath(
  Sys.getenv(
    "SI_FIGURE3_FIGURE_DIR",
    unset = file.path(workspace_root, "Figures")
  ),
  mustWork = FALSE
)
revised_root <- out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(revised_root, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  stopifnot(file.exists(path))
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

objective_display_floor <- 1e-4
context_colors <- c(
  "Separate in vivo" = "#0072B2",
  "Separate in vitro" = "#CC79A7"
)

theme_ms <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = "#333333"),
      axis.text = element_text(size = 7, colour = "#333333"),
      axis.title = element_text(size = 8),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      legend.title = element_text(size = 7.5, face = "bold"),
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7.5, colour = "#555555"),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 7, 5, 7)
    )
}

ivt_path <- file.path(
  diagnostic_root, "separate", "invitro",
  "seed_objective_simple.tsv"
)
viv_path <- file.path(
  diagnostic_root, "separate", "invivo",
  "seed_objective_simple.tsv"
)
ivt_seed_path <- file.path(
  diagnostic_root, "separate", "invitro",
  "seed_summary.tsv"
)
viv_seed_path <- file.path(
  diagnostic_root, "separate", "invivo",
  "seed_summary.tsv"
)
ivt_convergence_path <- file.path(
  diagnostic_root, "separate", "invitro",
  "convergence_summary.tsv"
)
viv_convergence_path <- file.path(
  diagnostic_root, "separate", "invivo",
  "convergence_summary.tsv"
)

ivt <- read_tsv(ivt_path) |>
  transmute(
    context = "Separate in vitro",
    seed, rank = objective_rank, objective,
    selected = seed == "seed10"
  )
viv <- read_tsv(viv_path) |>
  transmute(
    context = "Separate in vivo",
    seed, rank = objective_rank, objective,
    selected = seed == "seed25"
  )

separate <- bind_rows(ivt, viv) |>
  group_by(context) |>
  mutate(
    delta = objective - min(objective, na.rm = TRUE),
    delta_display = pmax(delta, objective_display_floor)
  ) |>
  ungroup()

p_a <- ggplot(separate, aes(rank, delta_display)) +
  geom_line(colour = "#777777", linewidth = 0.5) +
  geom_point(
    data = separate |> filter(selected),
    aes(fill = context),
    shape = 21, size = 2.7, stroke = 0.65,
    colour = "#333333"
  ) +
  scale_fill_manual(values = context_colors, guide = "none") +
  facet_wrap(~context, nrow = 1, scales = "free_y") +
  scale_y_log10(
    labels = label_number(accuracy = 0.0001),
    breaks = trans_breaks("log10", function(x) 10^x)
  ) +
  scale_x_continuous(breaks = c(1, 100, 250, 500)) +
  labs(
    tag = "A",
    title = "Separate-fit objective landscapes",
    subtitle = "Log10 y-axis; exact zero is displayed at 10^-4. Selected fits are highlighted.",
    x = "Objective rank among 500 starts",
    y = expression(Delta*" objective (log10; zero at 10"^{-4}*")")
  ) +
  theme_ms() +
  theme(panel.spacing = grid::unit(6, "mm"))

selection_path <- file.path(
  normalizePath(
    Sys.getenv(
      "FIGURE5_DATA_DIR",
      unset = file.path(
        workspace_root, "data", "Figures", "Figure5"
      )
    ),
    mustWork = TRUE
  ),
  "figure5_frozen_inputs",
  "selected_results.tsv"
)
selected <- read_tsv(selection_path) |>
  filter(record_type == "joint_pair_best") |>
  mutate(
    pair_label = sub(
      ".*_(C[0-9]+Sc[0-9]+)_vt.*",
      "\\1", warmup_label
    )
  )

joint_rows <- lapply(seq_len(nrow(selected)), function(i) {
  pair_dir <- file.path(
    diagnostic_root, "joint", selected$warmup_label[[i]],
    "seed_objective_simple.tsv"
  )
  read_tsv(pair_dir) |>
    transmute(
      pair_label = selected$pair_label[[i]],
      seed, rank = objective_rank, objective,
      selected = seed == selected$selected_seed[[i]]
    )
}) |>
  bind_rows() |>
  group_by(pair_label) |>
  mutate(
    delta = objective - min(objective, na.rm = TRUE),
    delta_display = pmax(delta, objective_display_floor)
  ) |>
  ungroup()

selected_invitro <- read_tsv(ivt_seed_path) |>
  filter(seed == "seed10")
selected_invivo <- read_tsv(viv_seed_path) |>
  filter(seed == "seed25")
if (nrow(selected_invitro) != 1L || nrow(selected_invivo) != 1L) {
  stop("Expected one saved summary row each for seed10 and seed25.")
}

invitro_convergence <- read_tsv(ivt_convergence_path)
invivo_convergence <- read_tsv(viv_convergence_path)
if (nrow(invitro_convergence) != 1L || nrow(invivo_convergence) != 1L) {
  stop("Expected one separate-fit convergence summary row per context.")
}

joint_selected_summaries <- lapply(seq_len(nrow(selected)), function(i) {
  summary_path <- file.path(
    diagnostic_root, "joint", selected$warmup_label[[i]],
    "seed_summary.tsv"
  )
  summary_row <- read_tsv(summary_path) |>
    filter(seed == selected$selected_seed[[i]])
  if (nrow(summary_row) != 1L) {
    stop("Expected one selected joint summary row for ", selected$warmup_label[[i]])
  }
  summary_row$pair_label <- selected$pair_label[[i]]
  summary_row
}) |>
  bind_rows()

joint_convergence <- lapply(seq_len(nrow(selected)), function(i) {
  convergence_path <- file.path(
    diagnostic_root, "joint", selected$warmup_label[[i]],
    "convergence_summary.tsv"
  )
  convergence_row <- read_tsv(convergence_path)
  if (nrow(convergence_row) != 1L) {
    stop("Expected one joint convergence row for ", selected$warmup_label[[i]])
  }
  convergence_row$pair_label <- selected$pair_label[[i]]
  convergence_row
}) |>
  bind_rows()

joint_cols <- c(
  C01Sc01 = "#009E73",
  C01Sc02 = "#CC79A7",
  C02Sc01 = "#E69F00",
  C02Sc02 = "#56B4E9",
  C03Sc01 = "#7A6F00",
  C03Sc02 = "#000000"
)

p_b <- ggplot(joint_rows, aes(rank, delta_display, colour = pair_label)) +
  geom_line(linewidth = 0.55, alpha = 0.82) +
  geom_point(
    data = joint_rows |> filter(selected),
    shape = 21, size = 2.2, stroke = 0.55,
    fill = "white"
  ) +
  facet_wrap(~pair_label, nrow = 2) +
  scale_colour_manual(values = joint_cols, guide = "none") +
  scale_y_log10(
    labels = label_number(accuracy = 0.0001),
    breaks = trans_breaks("log10", function(x) 10^x)
  ) +
  scale_x_continuous(breaks = c(1, 250, 500)) +
  labs(
    tag = "B",
    title = "Joint-fit objective landscapes",
    subtitle = paste(
      "Log10 y-axis; exact zero is displayed at 10^-4.",
      "Facets keep the six rank-1 selected winners distinct."
    ),
    x = "Objective rank within pair",
    y = expression(Delta*" objective (log10; zero at 10"^{-4}*")")
  ) +
  theme_ms() +
  theme(
    legend.position = "none",
    panel.spacing = grid::unit(3.5, "mm"),
    strip.text = element_text(size = 7, face = "bold")
  )

format_de <- function(x) {
  sprintf("%s / %s", x[["DEoptim converged"]], x[["Total seeds"]])
}

diag <- data.frame(
  fit = factor(
    c("Separate in vitro", "Separate in vivo", "Joint: six pair runs"),
    levels = rev(c("Separate in vitro", "Separate in vivo", "Joint: six pair runs"))
  ),
  starts = c(
    as.character(invitro_convergence$`Total seeds`[[1]]),
    as.character(invivo_convergence$`Total seeds`[[1]]),
    paste0(nrow(joint_convergence), " x ", unique(joint_convergence$`Total seeds`))
  ),
  de_flag = c(
    format_de(invitro_convergence[1, , drop = FALSE]),
    format_de(invivo_convergence[1, , drop = FALSE]),
    paste0(
      unique(joint_convergence$`DEoptim converged`), " / ",
      unique(joint_convergence$`Total seeds`), "\neach"
    )
  ),
  local = c(
    paste0(
      ifelse(selected_invitro$optimizer_local_accepted[[1]], "accepted", "not accepted"),
      "; code ", selected_invitro$optimizer_local_convergence[[1]]
    ),
    paste0(
      ifelse(selected_invivo$optimizer_local_accepted[[1]], "accepted", "not accepted"),
      "; code ", selected_invivo$optimizer_local_convergence[[1]]
    ),
    paste0(
      sum(joint_selected_summaries$optimizer_local_accepted, na.rm = TRUE),
      " / ", nrow(joint_selected_summaries), " accepted"
    )
  ),
  competitive = c(
    paste0(sum(ivt$objective - min(ivt$objective) <= 0.01), " within 0.01"),
    paste0(
      sum(viv$objective <= min(viv$objective) * 1.01), " within 1%\n",
      sum(viv$objective <= min(viv$objective) * 1.05), " within 5%"
    ),
    sprintf(
      "winner totals\n%.2f–%.2f",
      min(selected$objective), max(selected$objective)
    )
  ),
  bounds = c(
    paste0(
      selected_invitro$n_at_bound_active[[1]], " / ",
      selected_invitro$n_active_params[[1]]
    ),
    paste0(
      selected_invivo$n_at_bound_active[[1]], " / ",
      selected_invivo$n_active_params[[1]]
    ),
    paste0(
      min(joint_selected_summaries$n_at_bound_active), "–",
      max(joint_selected_summaries$n_at_bound_active), " / ",
      unique(joint_selected_summaries$n_active_params)
    )
  ),
  stringsAsFactors = FALSE
)

diag_long <- diag |>
  tidyr::pivot_longer(
    cols = -fit,
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("starts", "de_flag", "local", "competitive", "bounds"),
      labels = c(
        "Starts", "DEoptim\nconverged", "Local step",
        "Competitive\nsolutions", "At bounds"
      )
    )
  )

p_c <- ggplot(diag_long, aes(metric, fit)) +
  geom_tile(aes(fill = metric), colour = "white", linewidth = 0.6, alpha = 0.18) +
  geom_text(aes(label = value), size = 2.45, family = "sans", lineheight = 1.0) +
  scale_fill_manual(
    values = c(
      "#56B4E9", "#009E73", "#CC79A7", "#E69F00", "#999999"
    ),
    guide = "none"
  ) +
  labs(
    tag = "C",
    title = "Optimizer diagnostics and solution multiplicity",
    subtitle = paste(
      "Flags and optimizer codes are run diagnostics, not confidence intervals;\n",
      "none of these fits has bootstrap, profile-likelihood, or held-out uncertainty."
    ),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7, face = "bold", colour = "#333333"),
    axis.text.y = element_text(size = 7.5, colour = "#333333"),
    plot.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "#555555"),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(7, 7, 7, 7)
  )

figure <- (p_a / p_b / p_c) +
  plot_layout(heights = c(1, 1.25, 0.78)) +
  plot_annotation(
    title = "Fit-quality landscapes and optimizer diagnostics",
    theme = theme(
      plot.title = element_text(size = 11.5, face = "bold", color = "#111827"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

stem <- file.path(out_dir, "supp_fig3_fit_quality_optimizer_diagnostics")
ggsave(
  paste0(stem, ".png"), figure, width = 7.1, height = 7.7,
  units = "in", dpi = 300, bg = "white"
)
ggsave(
  paste0(stem, ".pdf"), figure, width = 7.1, height = 7.7,
  units = "in", device = grDevices::cairo_pdf, bg = "white"
)

write.table(
  diag,
  file.path(frozen_root, "si_figure3_optimizer_diagnostics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

for (extension in c("png", "pdf")) {
  source_path <- normalizePath(
    paste0(stem, ".", extension), mustWork = TRUE
  )
  destination_path <- normalizePath(
    file.path(
      revised_root,
      paste0(
        "supp_fig3_fit_quality_optimizer_diagnostics.", extension
      )
    ),
    mustWork = FALSE
  )
  if (!identical(source_path, destination_path)) {
    file.copy(source_path, destination_path, overwrite = TRUE)
  }
}

joint_input_paths <- unlist(lapply(seq_len(nrow(selected)), function(i) {
  file.path(
    diagnostic_root,
    "joint",
    selected$warmup_label[[i]],
    c("seed_objective_simple.tsv", "seed_summary.tsv", "convergence_summary.tsv")
  )
}), use.names = FALSE)
diagnostic_inputs <- c(
  selection_path,
  ivt_path,
  viv_path,
  ivt_seed_path,
  viv_seed_path,
  ivt_convergence_path,
  viv_convergence_path,
  joint_input_paths
)
if (any(!file.exists(diagnostic_inputs))) {
  stop("A diagnostic source disappeared during plotting")
}
provenance <- data.frame(
  role = c(
    "selected joint winner table",
    rep("optimizer diagnostic source", length(diagnostic_inputs) - 1L)
  ),
  path = normalizePath(diagnostic_inputs, mustWork = TRUE),
  md5 = unname(tools::md5sum(diagnostic_inputs)),
  stringsAsFactors = FALSE
)
write.table(
  provenance,
  file.path(frozen_root, "si_figure3_source_file_provenance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

validation <- data.frame(
  check = c(
    "separate_invitro_starts",
    "separate_invivo_starts",
    "joint_pairs",
    "joint_starts_each",
    "joint_selected_winners",
    "assembled_png",
    "assembled_pdf"
  ),
  observed = c(
    nrow(ivt),
    nrow(viv),
    nrow(joint_convergence),
    paste(unique(joint_convergence$`Total seeds`), collapse = ","),
    nrow(joint_selected_summaries),
    file.exists(paste0(stem, ".png")),
    file.exists(paste0(stem, ".pdf"))
  ),
  expected = c("500", "500", "6", "500", "6", "TRUE", "TRUE"),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(frozen_root, "si_figure3_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Supplementary optimizer diagnostics written to: ", paste0(stem, ".png"))

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_SI_Figure3 <- function() {
  data_dir <- file.path(DATA_ROOT, "SI_Figure3")
  require_files(file.path(
    data_dir, "separate", "invitro", "seed_objective_simple.tsv"
  ), "SI Figure 3 intermediate")
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_SI_Figure3.R"), mustWork = TRUE),
    env = c(
      "SI_FIGURE3_DRAW_WORKER=1",
      paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
      paste0("HYPOXIA_REPO_ROOT=", REPO_ROOT),
      paste0("SI_FIGURE3_DATA_DIR=", data_dir),
      paste0("SI_FIGURE3_FIGURE_DIR=", OUTPUT_ROOT),
      paste0("FIGURE5_DATA_DIR=", file.path(DATA_ROOT, "Figure5"))
    )
  )
  require_files(
    file.path(OUTPUT_ROOT, c(
      "supp_fig3_fit_quality_optimizer_diagnostics.png",
      "supp_fig3_fit_quality_optimizer_diagnostics.pdf"
    )),
    "SI Figure 3 output"
  )
  invisible(file.path(
    OUTPUT_ROOT,
    "supp_fig3_fit_quality_optimizer_diagnostics.png"
  ))
}

if (sys.nframe() == 0L) draw_SI_Figure3()

}
