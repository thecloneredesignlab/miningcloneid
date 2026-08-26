#!/usr/bin/env Rscript

if (identical(Sys.getenv("SUPP_FIGURE5_2_DRAW_WORKER"), "1")) {

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
workspace_root <- normalizePath(
  Sys.getenv("FIGURE_WORKSPACE_ROOT"), mustWork = TRUE
)
data_dir <- normalizePath(Sys.getenv("SUPP_FIGURE5_2_DATA_DIR"), mustWork = TRUE)
figure5_dir <- normalizePath(Sys.getenv("FIGURE5_DATA_DIR"), mustWork = TRUE)
out_dir <- normalizePath(
  Sys.getenv("SUPP_FIGURE5_2_FIGURE_DIR"), mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing Supplementary Figure 5-2 input: ", path)
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

selection_path <- file.path(
  figure5_dir, "figure5f_selected_pair_inputs.tsv"
)
selected <- read_tsv(selection_path)
required_selection <- c(
  "family", "warmup_label", "selected_seed", "invitro_seed",
  "selected_for_figure5f"
)
if (!all(required_selection %in% names(selected))) {
  stop("Figure 5 selected-pair table lacks required fields.")
}
selected <- selected |>
  filter(as.logical(selected_for_figure5f)) |>
  mutate(pair_label = family)
family_order <- c("C01", "C02", "C03")
selected <- selected[order(match(selected$family, family_order)), , drop = FALSE]
expected_selection <- c(
  C01 = "tsne_vi_seed366_C01Sc01_vt_seed10",
  C02 = "tsne_vi_seed25_C02Sc01_vt_seed10",
  C03 = "tsne_vi_seed311_C03Sc02_vt_seed10"
)
observed_selection <- setNames(selected$warmup_label, selected$family)
if (nrow(selected) != 3L ||
    !identical(selected$family, family_order) ||
    !identical(observed_selection[family_order], expected_selection) ||
    any(selected$invitro_seed != "seed10")) {
  stop(
    "Supplementary Figure 5-2 selected pairs do not match the approved ",
    "C01/C02/C03 primary-family inputs."
  )
}

objective_floor <- 1e-4
input_paths <- selection_path
joint_rows <- list()
diagnostic_rows <- list()
for (i in seq_len(nrow(selected))) {
  pair_root <- file.path(data_dir, "joint", selected$warmup_label[[i]])
  objective_path <- file.path(pair_root, "seed_objective_simple.tsv")
  optimizer_path <- file.path(pair_root, "seed_optimizer_diagnostics.tsv")
  seed_summary_path <- file.path(pair_root, "seed_summary.tsv")
  convergence_path <- file.path(pair_root, "convergence_summary.tsv")
  input_paths <- c(
    input_paths,
    objective_path,
    optimizer_path,
    seed_summary_path,
    convergence_path
  )
  objective <- read_tsv(objective_path) |>
    mutate(
      pair_label = selected$pair_label[[i]],
      selected = seed == selected$selected_seed[[i]],
      delta = objective - min(objective, na.rm = TRUE),
      delta_display = pmax(delta, objective_floor)
    )
  optimizer <- read_tsv(optimizer_path)
  seed_summary <- read_tsv(seed_summary_path) |>
    filter(seed == selected$selected_seed[[i]])
  convergence <- read_tsv(convergence_path)
  selected_optimizer <- optimizer |>
    filter(seed == selected$selected_seed[[i]])
  if (nrow(objective) != 500L || nrow(seed_summary) != 1L ||
      nrow(convergence) != 1L || nrow(selected_optimizer) != 1L) {
    stop("Incomplete joint diagnostic bundle for ", selected$warmup_label[[i]])
  }
  best <- min(objective$objective)
  relative_excess <- objective$objective / best - 1
  target_iterations <- unique(optimizer$optimizer_iter_target[
    is.finite(optimizer$optimizer_iter_target)
  ])
  if (length(target_iterations) != 1L) target_iterations <- NA_real_
  diagnostic_rows[[i]] <- data.frame(
    fit = selected$pair_label[[i]],
    starts = as.character(nrow(objective)),
    de_flag = paste0(
      convergence$`DEoptim convergence flag`[[1L]], "/",
      convergence$`Total seeds`[[1L]]
    ),
    iterations = paste0(
      min(optimizer$optimizer_iter_completed, na.rm = TRUE), "--",
      max(optimizer$optimizer_iter_completed, na.rm = TRUE), "/",
      target_iterations
    ),
    local = paste0(
      ifelse(
        seed_summary$optimizer_local_accepted[[1L]],
        "accepted", "not accepted"
      ),
      "; code ", seed_summary$optimizer_local_convergence[[1L]]
    ),
    competitive = paste0(
      sum(relative_excess <= 0.01), " <=1%\n",
      sum(relative_excess <= 0.05), " <=5%"
    ),
    bounds = paste0(
      seed_summary$n_at_bound_active[[1L]], "/",
      seed_summary$n_active_params[[1L]]
    ),
    stringsAsFactors = FALSE
  )
  joint_rows[[i]] <- objective
}
joint <- bind_rows(joint_rows)
diagnostic <- bind_rows(diagnostic_rows)

pair_colors <- c(
  C01 = "#C99700",
  C02 = "#6A3D9A",
  C03 = "#006D2C"
)

theme_si3 <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(size = base_size - 1, color = "#333333"),
      axis.title = element_text(size = base_size),
      strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB"),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 1.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.5, color = "#555555"),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 7, 6, 7)
    )
}

panel_a <- ggplot(joint, aes(objective_rank, delta_display, color = pair_label)) +
  geom_line(linewidth = 0.58, alpha = 0.86) +
  geom_point(
    data = joint |> filter(selected),
    shape = 21,
    size = 2.5,
    stroke = 0.65,
    fill = "white"
  ) +
  facet_wrap(~pair_label, nrow = 1) +
  scale_color_manual(values = pair_colors, guide = "none") +
  scale_x_continuous(breaks = c(1, 250, 500)) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_number(accuracy = 0.0001)
  ) +
  labs(
    tag = "A",
    title = "Joint-fit objective landscapes",
    subtitle = paste(
      "Three selected primary-family pairs, each with 500 numerical starts;",
      "the retained lowest objective is highlighted"
    ),
    x = "Objective rank within warm-start pair",
    y = expression(Delta*" objective (log"[10]*"; zero at 10"^{-4}*")")
  ) +
  theme_si3() +
  theme(panel.spacing = grid::unit(3.5, "mm"))

diagnostic_long <- diagnostic |>
  pivot_longer(
    cols = -fit,
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    fit = factor(fit, levels = rev(family_order)),
    metric = factor(
      metric,
      levels = c(
        "starts", "de_flag", "iterations", "local", "competitive", "bounds"
      ),
      labels = c(
        "Starts", "DEoptim\nflag", "Iterations", "Selected\nlocal step",
        "Near-optimal\nstarts", "Selected\nat bounds"
      )
    )
  )

panel_b <- ggplot(diagnostic_long, aes(metric, fit)) +
  geom_tile(aes(fill = metric), color = "white", linewidth = 0.6, alpha = 0.20) +
  geom_text(aes(label = value), size = 2.45, lineheight = 1.0) +
  scale_fill_manual(
    values = c("#56B4E9", "#009E73", "#7A6F00", "#CC79A7", "#E69F00", "#999999"),
    guide = "none"
  ) +
  labs(
    tag = "B",
    title = "Joint-search diagnostics by warm-start pair",
    subtitle = paste0(
      "Flags, codes, numerical-start spread, and bounds describe search behavior.\n",
      "They are not posterior uncertainty or biological replication."
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 8.5, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7.4, face = "bold", color = "#333333"),
    axis.text.y = element_text(size = 7.6, color = "#333333"),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 7.7, color = "#555555"),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(7, 7, 7, 7)
  )

figure <- (panel_a / panel_b) +
  plot_layout(heights = c(1.45, 0.85)) +
  plot_annotation(
    title = "Joint-fit numerical-search performance",
    theme = ggplot2::theme(
      plot.title = element_text(size = 12, face = "bold", color = "#111827"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

stem <- file.path(out_dir, "supp_fig5-2_joint_fit_optimizer_diagnostics")
ggsave(
  paste0(stem, ".png"), figure,
  width = 8.2, height = 5.8, units = "in", dpi = 300, bg = "white"
)
ggsave(
  paste0(stem, ".pdf"), figure,
  width = 8.2, height = 5.8, units = "in",
  device = grDevices::cairo_pdf, bg = "white"
)

write.table(
  diagnostic,
  file.path(data_dir, "supp_figure5-2_joint_optimizer_diagnostics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
provenance <- data.frame(
  role = c(
    "selected joint winner table",
    rep("joint optimizer diagnostic intermediate", length(input_paths) - 1L)
  ),
  path = normalizePath(input_paths, mustWork = TRUE),
  md5 = unname(tools::md5sum(input_paths)),
  stringsAsFactors = FALSE
)
write.table(
  provenance,
  file.path(data_dir, "supp_figure5-2_source_file_provenance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
validation <- data.frame(
  check = c(
    "content_scope", "joint_pairs", "joint_starts_each",
    "selected_winners", "assembled_png", "assembled_pdf"
  ),
  observed = c(
    "joint_only", nrow(diagnostic),
    paste(unique(diagnostic$starts), collapse = ","),
    sum(joint$selected),
    file.exists(paste0(stem, ".png")),
    file.exists(paste0(stem, ".pdf"))
  ),
  expected = c("joint_only", "3", "500", "3", "TRUE", "TRUE"),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(data_dir, "supp_figure5-2_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("Joint-only Supplementary Figure 5-2 written to: ", paste0(stem, ".png"))

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Supp_Figure5_2 <- function() {
  data_dir <- file.path(DATA_ROOT, "Supp_Figure5_2")
  require_files(
    file.path(
      DATA_ROOT, "Figure5", "figure5f_selected_pair_inputs.tsv"
    ),
    "Supplementary Figure 5-2 selected primary-family input"
  )
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_Supp_Figure5_2.R"), mustWork = TRUE),
    env = c(
      "SUPP_FIGURE5_2_DRAW_WORKER=1",
      paste0("FIGURE_WORKSPACE_ROOT=", WORKSPACE_ROOT),
      paste0("SUPP_FIGURE5_2_DATA_DIR=", data_dir),
      paste0("SUPP_FIGURE5_2_FIGURE_DIR=", OUTPUT_ROOT),
      paste0("FIGURE5_DATA_DIR=", file.path(DATA_ROOT, "Figure5"))
    )
  )
  outputs <- file.path(OUTPUT_ROOT, paste0(
    "supp_fig5-2_joint_fit_optimizer_diagnostics.", c("png", "pdf")
  ))
  require_files(outputs, "Supplementary Figure 5-2 output")
  invisible(outputs[[1L]])
}

if (sys.nframe() == 0L) draw_Supp_Figure5_2()

}
