#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)
out_root <- file.path(repo_root, "revised", "add", "round3_mandatory")
result_dir <- file.path(out_root, "results", "05_boundary_necrosis_audit")
figure_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
provenance_dir <- file.path(out_root, "provenance")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, name) write.table(
  x, file.path(result_dir, name), sep = "\t", row.names = FALSE,
  quote = FALSE, na = "NA"
)
get_metric <- function(tab, key) {
  x <- suppressWarnings(as.numeric(tab$value[tab$metric == key]))
  if (length(x)) x[[1]] else NA_real_
}

frozen_root <- file.path(
  repo_root, "revised", "iteration1", "data", "Figures", "Figure5",
  "figure5_frozen_inputs"
)
selection <- read_tsv(file.path(frozen_root, "selected_results.tsv"))
selection <- selection[selection$record_type == "joint_pair_best", ]
if (nrow(selection) != 6L) stop("Expected six joint winners")

# Boundary-loss decomposition from the already exported functional response grids.
boundary_rows <- list()
for (i in seq_len(nrow(selection))) {
  label <- selection$warmup_label[[i]]
  d <- file.path(frozen_root, "winners", label, "viz")
  for (context in c("invivo", "invitro")) {
    x <- read_tsv(file.path(d, context, "functional_curve_oxygen_multi_ploidy.tsv"))
    x$solution <- label
    x$context <- context
    x$objective <- selection$objective[[i]]
    x$boundary_share_of_cin_loss <- x$boundary_dropped_rate /
      pmax(x$boundary_dropped_rate + x$misseg_nonviable_rate, .Machine$double.xmin)
    boundary_rows[[length(boundary_rows) + 1L]] <- x
  }
}
boundary <- do.call(rbind, boundary_rows)
write_tsv(
  boundary[, c(
    "solution", "objective", "context", "oxygen_pct", "cohort", "N_ref",
    "misseg_nonviable_rate", "boundary_dropped_rate",
    "boundary_share_of_cin_loss"
  )],
  "boundary_loss_full_function_grid.tsv"
)

boundary_key <- interaction(
  boundary$solution, boundary$context, boundary$cohort, drop = TRUE
)
boundary_summary <- do.call(rbind, lapply(split(boundary, boundary_key), function(d) {
  data.frame(
    solution = d$solution[[1]],
    context = d$context[[1]],
    cohort = d$cohort[[1]],
    N_ref = d$N_ref[[1]],
    n_oxygen_grid = nrow(d),
    max_misseg_nonviable_rate = max(d$misseg_nonviable_rate),
    max_boundary_dropped_rate = max(d$boundary_dropped_rate),
    max_boundary_share_of_cin_loss = max(d$boundary_share_of_cin_loss),
    median_boundary_share_of_cin_loss = median(d$boundary_share_of_cin_loss)
  )
}))
write_tsv(boundary_summary, "boundary_loss_summary_by_solution_context_ploidy.tsv")

boundary_context_summary <- do.call(
  rbind,
  lapply(split(boundary, boundary$context), function(d) {
    data.frame(
      context = d$context[[1]],
      n_solutions = length(unique(d$solution)),
      n_grid_rows = nrow(d),
      maximum_boundary_dropped_rate = max(d$boundary_dropped_rate),
      maximum_misseg_nonviable_rate = max(d$misseg_nonviable_rate),
      maximum_boundary_share_of_cin_loss = max(d$boundary_share_of_cin_loss),
      fraction_grid_rows_boundary_share_ge_0p01 =
        mean(d$boundary_share_of_cin_loss >= 0.01),
      fraction_grid_rows_boundary_share_ge_0p001 =
        mean(d$boundary_share_of_cin_loss >= 0.001)
    )
  })
)
write_tsv(boundary_context_summary, "boundary_loss_context_summary.tsv")

# Reconstruct necrosis predictions from the exported terminal burden decomposition.
# This is algebraically identical to the C++ objective: dead volume / total volume.
invivo_root <- file.path(
  "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results",
  "fit_invivo_O2_buffering_500seed", "seed25"
)
necrosis_obs <- read_tsv(file.path(invivo_root, "necrosis_fit.tsv"))[
  , c("harvest", "cohort", "dose", "day", "obs_necrosis_fraction",
      "n_necrosis_obs")
]
eps <- 1e-4
sigma_necrosis <- 0.75

reconstruct_one <- function(
    solution, burden_path, fit_summary_path, objective_metric,
    source_type, objective_total
) {
  b <- read_tsv(burden_path)
  m <- merge(
    necrosis_obs,
    b[, c(
      "harvest", "day", "pred_burden_dead_total_volume_mm3",
      "pred_burden_dead_hypoxia_volume_mm3",
      "pred_burden_dead_buffer_volume_mm3",
      "pred_burden_volume_mm3"
    )],
    by = c("harvest", "day"), all.x = TRUE
  )
  m$pred_necrosis_fraction_reconstructed <-
    m$pred_burden_dead_total_volume_mm3 / m$pred_burden_volume_mm3
  m$residual_logit_reconstructed <-
    qlogis(pmin(pmax(m$pred_necrosis_fraction_reconstructed, eps), 1 - eps)) -
    qlogis(pmin(pmax(m$obs_necrosis_fraction, eps), 1 - eps))
  m$standardized_residual <-
    m$residual_logit_reconstructed / sigma_necrosis
  m$squared_standardized_residual <- m$standardized_residual^2
  m$solution <- solution
  m$source_type <- source_type
  m$objective_total <- objective_total
  fit_summary <- read_tsv(fit_summary_path)
  reported <- get_metric(fit_summary, objective_metric)
  reconstructed <- mean(m$squared_standardized_residual)
  objective_check <- data.frame(
    solution = solution,
    source_type = source_type,
    n_necrosis = nrow(m),
    sigma_necrosis_logit = sigma_necrosis,
    objective_necrosis_reported = reported,
    objective_necrosis_reconstructed = reconstructed,
    absolute_difference = abs(reported - reconstructed),
    neg2loglik_reported = if (source_type == "standalone") {
      get_metric(fit_summary, "objective_necrosis_neg2loglik_raw")
    } else {
      get_metric(fit_summary, "objective_invivo_necrosis_neg2loglik_raw")
    },
    neg2loglik_reconstructed = 2 * sum(m$squared_standardized_residual)
  )
  list(rows = m, check = objective_check)
}

necrosis_sets <- list()
necrosis_sets[[1]] <- reconstruct_one(
  "seed25_standalone",
  file.path(invivo_root, "burden_fit.tsv"),
  file.path(invivo_root, "fit_summary.tsv"),
  "objective_necrosis",
  "standalone",
  14.1193277940156
)
for (i in seq_len(nrow(selection))) {
  label <- selection$warmup_label[[i]]
  d <- file.path(frozen_root, "winners", label)
  necrosis_sets[[length(necrosis_sets) + 1L]] <- reconstruct_one(
    label,
    file.path(d, "invivo_burden_fit.tsv"),
    file.path(selection$source_dir[[i]], "fit_summary.tsv"),
    "objective_invivo_necrosis",
    "joint_winner",
    selection$objective[[i]]
  )
}
necrosis_rows <- do.call(rbind, lapply(necrosis_sets, `[[`, "rows"))
necrosis_checks <- do.call(rbind, lapply(necrosis_sets, `[[`, "check"))
write_tsv(necrosis_rows, "necrosis_predictions_reconstructed.tsv")
write_tsv(necrosis_checks, "necrosis_objective_reconstruction_check.tsv")

if (any(!is.finite(necrosis_rows$pred_necrosis_fraction_reconstructed))) {
  stop("Reconstructed necrosis predictions contain non-finite values")
}
if (max(necrosis_checks$absolute_difference) > 1e-8) {
  stop("Necrosis objective reconstruction did not reproduce reported objective")
}

# Explicit model implementation record.
soft_root <- "/Users/4482173/Documents/GitHub/soft_coupling"
soft_sha <- tryCatch(
  system2("git", c("-C", soft_root, "rev-parse", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_
)
implementation <- data.frame(
  component = c(
    "WGD probability",
    "WGD state mapping",
    "WGD branch offspring weight",
    "WGD generator combination",
    "Boundary handling",
    "Necrosis prediction",
    "Necrosis objective"
  ),
  implementation = c(
    "p_wgd is a constant per-division probability in the main model",
    "a source chromosome state N maps to doubled state 2N",
    "the WGD branch contributes weight +1 rather than +2",
    "non-WGD and WGD branch operators are mixed by 1-p_wgd and p_wgd before subtracting the source division operator",
    "default mode is drop; out-of-grid transition mass is not absorbed",
    "predicted necrosis is terminal dead volume divided by terminal total volume",
    "mean squared standardized logit residual with sigma 0.75; reported -2 log-likelihood is twice the summed squared standardized residual"
  ),
  source_file = c(
    "oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R",
    "oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp",
    "oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R",
    "oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R",
    "oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R",
    "oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invivo_backend.R",
    "oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp"
  ),
  source_lines = c(
    "837-849", "838-866", "15; 958-965", "998-1024",
    "16; 950-965", "2937-2971", "2352-2366; 2426-2429"
  ),
  stringsAsFactors = FALSE
)
write_tsv(implementation, "model_implementation_audit.tsv")
writeLines(
  c(
    paste0("soft_coupling_commit\t", soft_sha[[1]]),
    paste0("audit_script\t", script_path)
  ),
  file.path(provenance_dir, "model_implementation_audit_provenance.tsv")
)

plot_audit <- function(device, filename, width, height, res = NULL) {
  if (device == "pdf") pdf(filename, width = width, height = height, useDingbats = FALSE)
  else png(filename, width = width, height = height, units = "in", res = res)
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.2, 0.8), las = 1)
  plot_context_order <- c("invivo", "invitro")
  vals <- 100 * boundary_context_summary$maximum_boundary_share_of_cin_loss[
    match(plot_context_order, boundary_context_summary$context)
  ]
  barplot(
    vals, names.arg = c("in vivo", "in vitro"),
    col = c("#4E79A7", "#F28E2B"), log = "y",
    ylab = "Maximum boundary share of CIN loss (%)",
    main = "Grid-boundary loss audit"
  )
  abline(h = 1, lty = 2, col = "#555555")
  mtext("Dashed line: 1% materiality threshold", side = 3, line = 0.15, cex = 0.72)

  joint_nec <- necrosis_rows[necrosis_rows$source_type == "joint_winner", ]
  nec_env <- do.call(rbind, lapply(split(joint_nec, joint_nec$harvest), function(d) {
    data.frame(
      harvest = d$harvest[[1]],
      observed = d$obs_necrosis_fraction[[1]],
      pred_min = min(d$pred_necrosis_fraction_reconstructed),
      pred_median = median(d$pred_necrosis_fraction_reconstructed),
      pred_max = max(d$pred_necrosis_fraction_reconstructed)
    )
  }))
  plot(
    nec_env$observed, nec_env$pred_median,
    pch = 21, bg = "#59A14F", col = "white", cex = 1.2,
    xlab = "Observed necrosis fraction",
    ylab = "Median predicted necrosis fraction",
    main = "Reconstructed necrosis audit", xlim = c(0.5, 1), ylim = c(0.5, 1)
  )
  segments(
    nec_env$observed, nec_env$pred_min,
    nec_env$observed, nec_env$pred_max,
    col = adjustcolor("#333333", 0.45)
  )
  points(
    nec_env$observed, nec_env$pred_median,
    pch = 21, bg = "#59A14F", col = "white", cex = 1.2
  )
  abline(0, 1, lty = 2, col = "#555555")
  mtext("Bars: min-max across six joint winners", side = 3, line = 0.15, cex = 0.72)
}
plot_audit(
  "pdf", file.path(figure_dir, "round3_boundary_necrosis_audit.pdf"),
  7.5, 3.6
)
plot_audit(
  "png", file.path(figure_dir, "round3_boundary_necrosis_audit.png"),
  7.5, 3.6, 300
)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
standalone_check <- necrosis_checks[necrosis_checks$source_type == "standalone", ]
tex <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Boundary-loss, necrosis-export, and WGD implementation audit. Boundary fractions are evaluated over every O$_2$ and reference-ploidy row exported for the six frozen joint winners. Necrosis predictions were reconstructed from exported terminal dead and total volumes and reproduce the fitted objective.}",
  "\\label{tab:round3_boundary_necrosis_wgd}",
  "\\small",
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Diagnostic & In vivo & In vitro \\\\",
  "\\midrule",
  sprintf(
    "Maximum boundary-dropped rate & %.3g & %.3g \\\\",
    boundary_context_summary$maximum_boundary_dropped_rate[
      boundary_context_summary$context == "invivo"
    ],
    boundary_context_summary$maximum_boundary_dropped_rate[
      boundary_context_summary$context == "invitro"
    ]
  ),
  sprintf(
    "Maximum boundary share of CIN-associated loss & %.3g & %.3g \\\\",
    boundary_context_summary$maximum_boundary_share_of_cin_loss[
      boundary_context_summary$context == "invivo"
    ],
    boundary_context_summary$maximum_boundary_share_of_cin_loss[
      boundary_context_summary$context == "invitro"
    ]
  ),
  sprintf(
    "Standalone necrosis objective, reported/reconstructed & %s / %s & -- \\\\",
    fmt(standalone_check$objective_necrosis_reported),
    fmt(standalone_check$objective_necrosis_reconstructed)
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(tex, file.path(table_dir, "round3_boundary_necrosis_wgd_audit.tex"))

max_boundary_share <- max(boundary_context_summary$maximum_boundary_share_of_cin_loss)
max_boundary_fraction_material <- max(
  boundary_context_summary$fraction_grid_rows_boundary_share_ge_0p01
)
summary <- c(
  "# Boundary loss, necrosis, and WGD audit",
  "",
  sprintf(
    "Across the complete functional-response grids for all six frozen joint winners, the maximum boundary contribution is %.2f%% of the combined boundary plus within-grid CIN-associated nonviable loss, and as many as %.2f%% of rows in a context exceed the 1%% materiality threshold.",
    100 * max_boundary_share,
    100 * max_boundary_fraction_material
  ),
  sprintf(
    "The published seed-25 `necrosis_fit.tsv` contains six observed values but six missing predictions. Reconstructing predictions as terminal dead volume / terminal total volume yields six finite predictions and reproduces the reported normalized necrosis objective to maximum absolute error %.3g across the standalone fit and six joint winners.",
    max(necrosis_checks$absolute_difference)
  ),
  "",
  "WGD is implemented as a per-division competing branch: source state N maps to 2N with branch weight +1, mixed with the non-WGD branch by p_wgd, and out-of-grid transitions are dropped under the default boundary mode. The exact files and line ranges are recorded in `model_implementation_audit.tsv`.",
  "",
  "Interpretation: boundary truncation is material for a subset of in-vivo functional-grid evaluations and therefore requires the separate expanded-grid sensitivity audit. The missing necrosis export was also a real audit defect; the reconstructed table repairs auditability without refitting and does not change the objective."
)
writeLines(summary, file.path(result_dir, "boundary_necrosis_wgd_audit_summary.md"))
message("Wrote boundary, necrosis, and WGD audit")
