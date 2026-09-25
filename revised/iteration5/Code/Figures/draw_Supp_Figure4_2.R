#!/usr/bin/env Rscript

if (identical(Sys.getenv("SUPP_FIGURE4_2_DRAW_WORKER"), "1")) {

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

data_dir <- normalizePath(Sys.getenv("SUPP_FIGURE4_2_DATA_DIR"), mustWork = TRUE)
out_dir <- normalizePath(
  Sys.getenv("SUPP_FIGURE4_2_FIGURE_DIR"), mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(name) {
  path <- file.path(data_dir, name)
  if (!file.exists(path)) stop("Missing Supplementary Figure 4-2 input: ", path)
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

objective_path <- file.path(data_dir, "seed_objective_simple.tsv")
optimizer_path <- file.path(data_dir, "seed_optimizer_diagnostics.tsv")
selected_path <- file.path(data_dir, "seed_summary.tsv")
convergence_path <- file.path(data_dir, "convergence_summary.tsv")
boundary_path <- file.path(data_dir, "seed_boundary_summary.tsv")
near_path <- file.path(data_dir, "objective_near_optimal_summary.tsv")
termination_path <- file.path(data_dir, "termination_summary.tsv")
input_paths <- c(
  objective_path, optimizer_path, selected_path, convergence_path,
  boundary_path, near_path, termination_path
)
if (any(!file.exists(input_paths))) {
  stop("Incomplete Supplementary Figure 4-2 diagnostic bundle.")
}

objective <- read_tsv("seed_objective_simple.tsv")
optimizer <- read_tsv("seed_optimizer_diagnostics.tsv")
selected <- read_tsv("seed_summary.tsv")
convergence <- read_tsv("convergence_summary.tsv")
boundary <- read_tsv("seed_boundary_summary.tsv")
near <- read_tsv("objective_near_optimal_summary.tsv")
termination <- read_tsv("termination_summary.tsv")
if (nrow(objective) != 500L || nrow(optimizer) != 500L ||
    nrow(boundary) != 500L || nrow(selected) != 1L ||
    selected$seed[[1L]] != "seed25") {
  stop("Supplementary Figure 4-2 requires 500 in-vivo starts and selected seed25.")
}

objective_floor <- 1e-4
best_objective <- min(objective$objective)
objective <- objective |>
  mutate(
    delta = objective - best_objective,
    delta_display = pmax(delta, objective_floor),
    relative_excess_pct = 100 * (objective / best_objective - 1),
    selected = seed == "seed25"
  )
relative_q99 <- unname(quantile(
  objective$relative_excess_pct, 0.99, names = FALSE, type = 8
))
near_1 <- near$n_seeds[
  near$threshold_type == "relative" & abs(near$threshold - 0.01) < 1e-12
][[1L]]
near_5 <- near$n_seeds[
  near$threshold_type == "relative" & abs(near$threshold - 0.05) < 1e-12
][[1L]]

theme_si4 <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(size = base_size - 1, color = "#333333"),
      axis.title = element_text(size = base_size),
      plot.title = element_text(size = base_size + 1.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.5, color = "#555555"),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.caption = element_text(size = base_size - 1, color = "#555555", hjust = 0),
      plot.margin = margin(7, 8, 6, 7)
    )
}

panel_a <- ggplot(objective, aes(objective_rank, delta_display)) +
  geom_line(color = "#666666", linewidth = 0.65) +
  geom_point(
    data = objective |> filter(selected),
    shape = 21, size = 3, stroke = 0.7,
    color = "#111111", fill = "#0072B2"
  ) +
  scale_x_continuous(breaks = c(1, 100, 250, 500)) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_number(accuracy = 0.0001)
  ) +
  labs(
    tag = "A",
    title = "In-vivo objective landscape",
    subtitle = "seed25 is the lowest retained objective among 500 starts",
    x = "Objective rank",
    y = expression(Delta*" objective (log"[10]*"; zero at 10"^{-4}*")")
  ) +
  theme_si4()

panel_b <- ggplot(objective, aes(relative_excess_pct)) +
  stat_ecdf(geom = "step", color = "#0072B2", linewidth = 0.85) +
  geom_vline(
    xintercept = c(1, 5),
    color = c("#D55E00", "#009E73"),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  annotate(
    "label",
    x = c(1, 5),
    y = c(0.19, 0.57),
    label = c(
      paste0(near_1, "/500 within 1%"),
      paste0(near_5, "/500 within 5%")
    ),
    hjust = c(-0.05, -0.05),
    size = 2.7,
    color = c("#D55E00", "#009E73"),
    fill = "white",
    linewidth = 0.2
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, relative_q99), clip = "on") +
  labs(
    tag = "B",
    title = "Near-optimal solution multiplicity",
    subtitle = "ECDF with a 99th-percentile coordinate zoom",
    x = "Objective excess above seed25 (%)",
    y = "Fraction of numerical starts",
    caption = "Coordinate zoom retains all 500 starts in the analysis."
  ) +
  theme_si4()

selected_optimizer <- optimizer |>
  filter(seed == "seed25")
stop_reason_counts <- termination |>
  filter(diagnostic == "DEoptim stop reason")
stop_reason_label <- paste0(
  stop_reason_counts$category, " (", stop_reason_counts$n_seeds, ")",
  collapse = "; "
)
iteration_range <- range(
  optimizer$optimizer_iter_completed,
  na.rm = TRUE
)
iteration_target <- unique(optimizer$optimizer_iter_target[
  is.finite(optimizer$optimizer_iter_target)
])
local_codes <- termination |>
  filter(diagnostic == "Local convergence code")
local_code_label <- paste0(
  local_codes$category, " (", local_codes$n_seeds, ")",
  collapse = "; "
)
diagnostic <- data.frame(
  metric = factor(
    c(
      "Starts", "DEoptim flag", "Stop reason", "Iterations",
      "Local accepted", "Local codes", "Selected seed", "Selected bounds"
    ),
    levels = rev(c(
      "Starts", "DEoptim flag", "Stop reason", "Iterations",
      "Local accepted", "Local codes", "Selected seed", "Selected bounds"
    ))
  ),
  value = c(
    as.character(convergence$`Total seeds`[[1L]]),
    paste0(
      convergence$`DEoptim convergence flag`[[1L]], "/",
      convergence$`Total seeds`[[1L]]
    ),
    stop_reason_label,
    paste0(iteration_range[[1L]], "--", iteration_range[[2L]], "/", iteration_target),
    paste0(
      convergence$`Local accepted`[[1L]], "/",
      convergence$`Total seeds`[[1L]]
    ),
    local_code_label,
    paste0(
      "seed25; local ",
      ifelse(selected$optimizer_local_accepted[[1L]], "accepted", "not accepted"),
      "; code ", selected$optimizer_local_convergence[[1L]]
    ),
    paste0(
      selected$n_at_bound_active[[1L]], "/",
      selected$n_active_params[[1L]], " active parameters"
    )
  ),
  stringsAsFactors = FALSE
)

panel_c <- ggplot(diagnostic, aes(1, metric)) +
  geom_tile(fill = "#E8F1F8", color = "white", linewidth = 0.7) +
  geom_text(aes(label = value), size = 2.65, lineheight = 1.0) +
  scale_x_continuous(breaks = NULL) +
  labs(
    tag = "C",
    title = "Optimizer diagnostics",
    subtitle = "Flags and codes describe numerical stopping,\nnot global convergence",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 8.5, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7.6, face = "bold", color = "#333333"),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 7.7, color = "#555555"),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(7, 8, 6, 7)
  )

boundary$local_status <- ifelse(
  boundary$optimizer_local_accepted, "Local accepted", "Local not accepted"
)
panel_d <- ggplot(
  boundary,
  aes(n_at_bound_active, fill = local_status)
) +
  geom_bar(color = "white", linewidth = 0.3) +
  geom_vline(
    xintercept = selected$n_at_bound_active[[1L]],
    color = "#D55E00", linetype = "dashed", linewidth = 0.75
  ) +
  annotate(
    "label",
    x = selected$n_at_bound_active[[1L]],
    y = Inf,
    label = paste0(
      "seed25: ", selected$n_at_bound_active[[1L]], "/",
      selected$n_active_params[[1L]]
    ),
    vjust = 1.3,
    size = 2.7,
    color = "#D55E00",
    fill = "white",
    linewidth = 0.2
  ) +
  scale_x_continuous(breaks = sort(unique(boundary$n_at_bound_active))) +
  scale_fill_manual(
    values = c(
      "Local accepted" = "#0072B2",
      "Local not accepted" = "#B8B8B8"
    ),
    drop = FALSE
  ) +
  labs(
    tag = "D",
    title = "Boundary dependence",
    subtitle = "Exact active-bound counts across 500 starts",
    x = "Active parameters exactly at a configured bound",
    y = "Number of starts",
    fill = NULL
  ) +
  theme_si4() +
  theme(legend.position = "bottom")

figure <- ((panel_a | panel_b) / (panel_c | panel_d)) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Separate in-vivo numerical-search performance",
    subtitle = paste(
      "Optimizer starts are numerical endpoints, not posterior draws,",
      "confidence intervals, or biological replicates"
    ),
    theme = ggplot2::theme(
      plot.title = element_text(size = 12, face = "bold", color = "#111827"),
      plot.subtitle = element_text(size = 8.5, color = "#4B5563"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

stem <- file.path(out_dir, "supp_fig4-2_invivo_optimizer_diagnostics")
ggsave(
  paste0(stem, ".png"), figure,
  width = 7.1, height = 7.1, units = "in", dpi = 300, bg = "white"
)
ggsave(
  paste0(stem, ".pdf"), figure,
  width = 7.1, height = 7.1, units = "in",
  device = grDevices::cairo_pdf, bg = "white"
)

write.table(
  diagnostic,
  file.path(data_dir, "supp_figure4-2_optimizer_diagnostics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
display_range <- data.frame(
  variable = "relative_objective_excess_percent",
  display_lower = 0,
  display_upper = relative_q99,
  percentile = 0.99,
  rule = "coordinate zoom only; all starts retained",
  stringsAsFactors = FALSE
)
write.table(
  display_range,
  file.path(data_dir, "supp_figure4-2_display_range.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
provenance <- data.frame(
  role = "in-vivo optimizer diagnostic intermediate",
  path = normalizePath(input_paths, mustWork = TRUE),
  md5 = unname(tools::md5sum(input_paths)),
  stringsAsFactors = FALSE
)
write.table(
  provenance,
  file.path(data_dir, "supp_figure4-2_source_file_provenance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
validation <- data.frame(
  check = c(
    "in_vivo_starts", "selected_seed", "near_1pct", "near_5pct",
    "boundary_rows", "assembled_png", "assembled_pdf", "evidence_type"
  ),
  observed = c(
    nrow(objective), selected$seed[[1L]], near_1, near_5,
    nrow(boundary), file.exists(paste0(stem, ".png")),
    file.exists(paste0(stem, ".pdf")), "numerical_search_diagnostics"
  ),
  expected = c(
    "500", "seed25", "29", "241", "500", "TRUE", "TRUE",
    "numerical_search_diagnostics"
  ),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(data_dir, "supp_figure4-2_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("Supplementary Figure 4-2 written to: ", paste0(stem, ".png"))

} else {

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

draw_Supp_Figure4_2 <- function() {
  data_dir <- file.path(DATA_ROOT, "Supp_Figure4_2")
  require_files(
    file.path(data_dir, "seed_boundary_summary.tsv"),
    "Supplementary Figure 4-2 intermediate"
  )
  run_process(
    "Rscript",
    normalizePath(file.path(script_dir, "draw_Supp_Figure4_2.R"), mustWork = TRUE),
    env = c(
      "SUPP_FIGURE4_2_DRAW_WORKER=1",
      paste0("SUPP_FIGURE4_2_DATA_DIR=", data_dir),
      paste0("SUPP_FIGURE4_2_FIGURE_DIR=", OUTPUT_ROOT)
    )
  )
  outputs <- file.path(OUTPUT_ROOT, paste0(
    "supp_fig4-2_invivo_optimizer_diagnostics.", c("png", "pdf")
  ))
  require_files(outputs, "Supplementary Figure 4-2 output")
  invisible(outputs[[1L]])
}

if (sys.nframe() == 0L) draw_Supp_Figure4_2()

}
