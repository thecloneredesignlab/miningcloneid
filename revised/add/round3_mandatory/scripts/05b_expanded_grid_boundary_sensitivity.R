#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)
out_root <- file.path(repo_root, "revised", "add", "round3_mandatory")
result_dir <- file.path(out_root, "results", "05_boundary_necrosis_audit")
figure_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

soft_root <- "/Users/4482173/Documents/GitHub/soft_coupling"
model_file <- file.path(
  soft_root, "oxygen", "code", "O2_supply_demand_MAP", "model",
  "model_O2_supply_demand_MAP.R"
)
source(model_file)
if (!requireNamespace("RSpectra", quietly = TRUE)) {
  stop("RSpectra is required for sparse dominant-eigenvalue comparison")
}

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, name) write.table(
  x, file.path(result_dir, name), sep = "\t", row.names = FALSE,
  quote = FALSE, na = "NA"
)
as_params <- function(path) {
  x <- read_tsv(path)
  setNames(as.numeric(x$value), x$parameter)
}
param <- function(p, name, fallback = NA_real_) {
  x <- p[[name]]
  if (is.null(x) || !length(x) || !is.finite(x[[1]])) fallback else x[[1]]
}

call_core <- function(p, cfg, O2, n_min, n_max) {
  cpp_o2simps_build_G_for_o2_triplet(
    O2 = as.numeric(O2),
    O2_crit = param(p, "O2_crit"),
    N0min = as.integer(n_min),
    N0max = as.integer(n_max),
    N1min = as.integer(n_min),
    N1max = as.integer(n_max),
    lam_max = param(p, "lam_max"),
    p_mis_base = param(p, "p_mis_base", 1e-5),
    p_misseg = param(p, "p_misseg", 0),
    k_o_mis = param(p, "k_o_mis", 50),
    p_wgd = param(p, "p_wgd", 0),
    boundary = "drop",
    eps_tail = 1e-8,
    buffer_smax = param(p, "buffer_smax", 1),
    buffer_beta = param(p, "buffer_beta", 0),
    buffer_n_exp = param(p, "buffer_n_exp", 1),
    N_unit = as.integer(cfg$N_UNIT),
    beta_size = 0,
    O2_growth = isTRUE(cfg$O2_growth),
    alpha_o2 = param(p, "alpha_o2", 0.5),
    gamma_growth = param(p, "gamma_growth", 2),
    mu_hp = param(p, "mu_hp", 0),
    gamma_mu = param(p, "gamma_mu", 0),
    n_O = param(p, "n_O", 1.5),
    ploidy_O2_death = as.character(cfg$ploidy_O2_death)
  )
}

frozen_root <- file.path(
  repo_root, "revised", "iteration1", "data", "Figures", "Figure5",
  "figure5_frozen_inputs"
)
selection <- read_tsv(file.path(frozen_root, "selected_results.tsv"))
selection <- selection[selection$record_type == "joint_pair_best", ]

expanded_rows <- list()
validation_rows <- list()
for (i in seq_len(nrow(selection))) {
  label <- selection$warmup_label[[i]]
  d <- file.path(frozen_root, "winners", label)
  cfg <- readRDS(file.path(d, "fit_config.rds"))
  for (context in c("invivo", "invitro")) {
    p <- if (context == "invivo") {
      as_params(file.path(d, "best_params.tsv"))
    } else {
      as_params(file.path(d, "invitro_effective_params.tsv"))
    }
    curve <- read_tsv(file.path(
      d, "viz", context, "functional_curve_oxygen_multi_ploidy.tsv"
    ))
    row_groups <- split(seq_len(nrow(curve)), sprintf("%.12f", curve$oxygen_pct))
    current_calc_boundary <- rep(NA_real_, nrow(curve))
    current_calc_nonviable <- rep(NA_real_, nrow(curve))
    expanded_boundary <- rep(NA_real_, nrow(curve))
    expanded_nonviable <- rep(NA_real_, nrow(curve))
    for (idx in row_groups) {
      O2 <- curve$oxygen_pct[idx[[1]]]
      tri_current <- call_core(p, cfg, O2, cfg$N_MIN, cfg$N_MAX)
      tri_expanded <- call_core(p, cfg, O2, cfg$N_MIN, 308L)
      current_state <- as.integer(round(curve$N_ref[idx])) - cfg$N_MIN + 1L
      expanded_state <- as.integer(round(curve$N_ref[idx])) - cfg$N_MIN + 1L
      current_calc_boundary[idx] <- tri_current$boundary_dropped_rate[current_state]
      current_calc_nonviable[idx] <- tri_current$misseg_nonviable_rate[current_state]
      expanded_boundary[idx] <- tri_expanded$boundary_dropped_rate[expanded_state]
      expanded_nonviable[idx] <- tri_expanded$misseg_nonviable_rate[expanded_state]
    }
    validation_rows[[length(validation_rows) + 1L]] <- data.frame(
      solution = label,
      context = context,
      n_rows = nrow(curve),
      max_abs_boundary_export_vs_recomputed = max(abs(
        curve$boundary_dropped_rate - current_calc_boundary
      )),
      max_abs_nonviable_export_vs_recomputed = max(abs(
        curve$misseg_nonviable_rate - current_calc_nonviable
      ))
    )
    expanded_rows[[length(expanded_rows) + 1L]] <- data.frame(
      solution = label,
      objective = selection$objective[[i]],
      context = context,
      oxygen_pct = curve$oxygen_pct,
      cohort = curve$cohort,
      N_ref = curve$N_ref,
      current_N_min = cfg$N_MIN,
      current_N_max = cfg$N_MAX,
      expanded_N_min = cfg$N_MIN,
      expanded_N_max = 308L,
      current_misseg_nonviable_rate = current_calc_nonviable,
      expanded_misseg_nonviable_rate = expanded_nonviable,
      current_boundary_dropped_rate = current_calc_boundary,
      expanded_boundary_dropped_rate = expanded_boundary,
      upper_truncation_rate_removed = pmax(
        current_calc_boundary - expanded_boundary, 0
      ),
      current_boundary_share = current_calc_boundary / pmax(
        current_calc_boundary + current_calc_nonviable,
        .Machine$double.xmin
      ),
      expanded_boundary_share = expanded_boundary / pmax(
        expanded_boundary + expanded_nonviable,
        .Machine$double.xmin
      ),
      upper_truncation_share_of_current_cin_loss = pmax(
        current_calc_boundary - expanded_boundary, 0
      ) / pmax(
        current_calc_boundary + current_calc_nonviable,
        .Machine$double.xmin
      )
    )
  }
}
expanded <- do.call(rbind, expanded_rows)
validation <- do.call(rbind, validation_rows)
write_tsv(expanded, "expanded_grid_functional_boundary_comparison.tsv")
write_tsv(validation, "expanded_grid_recomputation_validation.tsv")
if (max(validation$max_abs_boundary_export_vs_recomputed) > 2e-10 ||
    max(validation$max_abs_nonviable_export_vs_recomputed) > 2e-10) {
  stop("Current-grid recomputation does not reproduce exported functional curves")
}

expanded_summary <- do.call(rbind, lapply(split(expanded, expanded$context), function(d) {
  data.frame(
    context = d$context[[1]],
    n_solutions = length(unique(d$solution)),
    n_grid_rows = nrow(d),
    current_max_boundary_share = max(d$current_boundary_share),
    expanded_max_boundary_share = max(d$expanded_boundary_share),
    max_upper_truncation_share_removed =
      max(d$upper_truncation_share_of_current_cin_loss),
    current_fraction_boundary_share_ge_0p01 = mean(d$current_boundary_share >= 0.01),
    expanded_fraction_boundary_share_ge_0p01 = mean(d$expanded_boundary_share >= 0.01),
    fraction_upper_truncation_share_ge_0p01 =
      mean(d$upper_truncation_share_of_current_cin_loss >= 0.01),
    max_abs_change_misseg_nonviable_rate = max(abs(
      d$expanded_misseg_nonviable_rate - d$current_misseg_nonviable_rate
    ))
  )
}))
write_tsv(expanded_summary, "expanded_grid_functional_boundary_summary.tsv")

dominant_attractor <- function(p, cfg, O2, n_min, n_max) {
  tri <- call_core(p, cfg, O2, n_min, n_max)
  G <- Matrix::sparseMatrix(
    i = as.integer(tri$i), j = as.integer(tri$j), x = as.numeric(tri$x),
    dims = c(as.integer(tri$nrow), as.integer(tri$ncol))
  )
  eg <- tryCatch(
    suppressWarnings(RSpectra::eigs(
      G, k = 2L, which = "LR",
      opts = list(
        ncv = min(60L, nrow(G) - 1L),
        maxitr = 100000L,
        tol = 1e-10
      )
    )),
    error = function(e) NULL
  )
  if (is.null(eg) || length(eg$values) < 2L) {
    dense <- eigen(as.matrix(G))
    eg <- list(values = dense$values, vectors = dense$vectors)
  }
  ord <- order(Re(eg$values), decreasing = TRUE)
  values <- eg$values[ord]
  v <- Re(eg$vectors[, ord[[1]]])
  if (sum(v) < 0) v <- -v
  if (any(v < -1e-7)) {
    return(data.frame(
      status = "dominant_vector_negative_components",
      dominant_mean_ploidy = NA_real_,
      dominant_growth_rate = Re(values[[1]]),
      spectral_gap = Re(values[[1]]) - Re(values[[2]])
    ))
  }
  v <- pmax(v, 0)
  v <- v / sum(v)
  N <- seq.int(n_min, n_max)
  data.frame(
    status = "ok",
    dominant_mean_ploidy = sum((N / cfg$N_UNIT) * v),
    dominant_growth_rate = Re(values[[1]]),
    spectral_gap = Re(values[[1]]) - Re(values[[2]])
  )
}

reporting_o2 <- c(0, 0.1, 0.5, 1, 2, 5)
eigen_rows <- list()
for (i in seq_len(nrow(selection))) {
  label <- selection$warmup_label[[i]]
  d <- file.path(frozen_root, "winners", label)
  cfg <- readRDS(file.path(d, "fit_config.rds"))
  p <- as_params(file.path(d, "best_params.tsv"))
  for (O2 in reporting_o2) {
    current <- dominant_attractor(p, cfg, O2, cfg$N_MIN, cfg$N_MAX)
    expanded_one <- dominant_attractor(p, cfg, O2, cfg$N_MIN, 308L)
    eigen_rows[[length(eigen_rows) + 1L]] <- data.frame(
      solution = label,
      objective = selection$objective[[i]],
      oxygen_pct = O2,
      current_N_min = cfg$N_MIN,
      current_N_max = cfg$N_MAX,
      expanded_N_min = cfg$N_MIN,
      expanded_N_max = 308L,
      current_status = current$status,
      expanded_status = expanded_one$status,
      current_dominant_mean_ploidy = current$dominant_mean_ploidy,
      expanded_dominant_mean_ploidy = expanded_one$dominant_mean_ploidy,
      expanded_minus_current_ploidy =
        expanded_one$dominant_mean_ploidy - current$dominant_mean_ploidy,
      current_dominant_growth_rate = current$dominant_growth_rate,
      expanded_dominant_growth_rate = expanded_one$dominant_growth_rate,
      current_spectral_gap = current$spectral_gap,
      expanded_spectral_gap = expanded_one$spectral_gap
    )
  }
}
eigen_compare <- do.call(rbind, eigen_rows)
write_tsv(eigen_compare, "expanded_grid_dominant_attractor_comparison.tsv")

eigen_summary <- do.call(
  rbind,
  lapply(split(eigen_compare, eigen_compare$oxygen_pct), function(d) {
    data.frame(
      oxygen_pct = d$oxygen_pct[[1]],
      n_solutions = nrow(d),
      median_current_ploidy = median(d$current_dominant_mean_ploidy),
      median_expanded_ploidy = median(d$expanded_dominant_mean_ploidy),
      median_expanded_minus_current = median(d$expanded_minus_current_ploidy),
      max_abs_expanded_minus_current = max(abs(d$expanded_minus_current_ploidy)),
      direction_fraction_expanded_gt_current =
        mean(d$expanded_minus_current_ploidy > 0)
    )
  })
)
write_tsv(eigen_summary, "expanded_grid_dominant_attractor_summary.tsv")

plot_expanded <- function(device, filename, width, height, res = NULL) {
  if (device == "pdf") pdf(filename, width = width, height = height, useDingbats = FALSE)
  else png(filename, width = width, height = height, units = "in", res = res)
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.4, 4.5, 2.2, 0.8), las = 1)
  low_boundary_vals <- 100 * expanded_summary$expanded_max_boundary_share
  upper_truncation_vals <- 100 * expanded_summary$max_upper_truncation_share_removed
  mat <- rbind(pmax(low_boundary_vals, 1e-12), pmax(upper_truncation_vals, 1e-12))
  colnames(mat) <- c("in vitro", "in vivo")[
    match(expanded_summary$context, c("invitro", "invivo"))
  ]
  barplot(
    mat, beside = TRUE, log = "y",
    col = c("#E15759", "#59A14F"),
    ylab = "Maximum boundary share of CIN loss (%)",
    main = "Boundary component audit"
  )
  abline(h = 1, lty = 2, col = "#555555")
  legend(
    "topright",
    legend = c("Retained lower-boundary loss", "Removed upper truncation"),
    fill = c("#E15759", "#59A14F"), bty = "n", cex = 0.72
  )
  plot(
    eigen_compare$current_dominant_mean_ploidy,
    eigen_compare$expanded_dominant_mean_ploidy,
    pch = 21, bg = "#4E79A7", col = "white",
    xlab = "Current-grid dominant mean ploidy",
    ylab = "Expanded-grid dominant mean ploidy",
    main = "Six joint winners at six O2 levels"
  )
  abline(0, 1, lty = 2, col = "#555555")
}
plot_expanded(
  "pdf", file.path(figure_dir, "round3_expanded_grid_boundary_sensitivity.pdf"),
  7.5, 3.6
)
plot_expanded(
  "png", file.path(figure_dir, "round3_expanded_grid_boundary_sensitivity.png"),
  7.5, 3.6, 300
)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
invivo_summary <- expanded_summary[expanded_summary$context == "invivo", ]
tex <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Upper-boundary chromosome-grid sensitivity using the six frozen joint-winner parameter sets without refitting. The original $N=22$--154 grid was recomputed exactly and compared with $N=22$--308, preserving the biological lower bound while testing upper truncation.}",
  "\\label{tab:round3_expanded_grid}",
  "\\small",
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Diagnostic & Current grid & Expanded grid \\\\",
  "\\midrule",
  sprintf(
    "Maximum in-vivo boundary share of CIN loss & %s & %s \\\\",
    fmt(invivo_summary$current_max_boundary_share),
    fmt(invivo_summary$expanded_max_boundary_share)
  ),
  sprintf(
    "Fraction of in-vivo response rows above 1\\%% total boundary share & %s & %s \\\\",
    fmt(invivo_summary$current_fraction_boundary_share_ge_0p01),
    fmt(invivo_summary$expanded_fraction_boundary_share_ge_0p01)
  ),
  sprintf(
    "Maximum upper-truncation share removed & -- & %s \\\\",
    fmt(invivo_summary$max_upper_truncation_share_removed)
  ),
  sprintf(
    "Maximum absolute change in within-grid nonviable rate & -- & %.3g \\\\",
    invivo_summary$max_abs_change_misseg_nonviable_rate
  ),
  sprintf(
    "Maximum absolute change in dominant mean ploidy & -- & %s \\\\",
    fmt(max(abs(eigen_compare$expanded_minus_current_ploidy)))
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(tex, file.path(table_dir, "round3_expanded_grid_boundary_sensitivity.tex"))

summary <- c(
  "# Expanded-grid boundary sensitivity",
  "",
  sprintf(
    "Direct recomputation reproduced the exported current-grid boundary and within-grid nonviable rates to maximum absolute errors %.3g and %.3g.",
    max(validation$max_abs_boundary_export_vs_recomputed),
    max(validation$max_abs_nonviable_export_vs_recomputed)
  ),
  sprintf(
    "On the current in-vivo grid, the maximum total boundary share is %.1f%% and %.1f%% of evaluated rows exceed 1%%. Expanding only the upper boundary from N=154 to N=308, while retaining the biological lower bound N=22, attributes at most %.3g%% of current CIN-associated loss to removable upper truncation and changes the within-grid nonviable rate by at most %.3g.",
    100 * invivo_summary$current_max_boundary_share,
    100 * invivo_summary$current_fraction_boundary_share_ge_0p01,
    100 * invivo_summary$max_upper_truncation_share_removed,
    invivo_summary$max_abs_change_misseg_nonviable_rate
  ),
  sprintf(
    "Across the six joint winners and O2={0,0.1,0.5,1,2,5}%%, the maximum absolute change in dominant mean ploidy is %.3f ploidy units.",
    max(abs(eigen_compare$expanded_minus_current_ploidy))
  ),
  "",
  "The retained lower-boundary component represents daughter states below the model's viable chromosome range and is a biological nonviability convention, not a numerical state to restore. Interpretation of upper-grid robustness must follow the dominant-attractor comparison. No parameters were refitted and no original result was overwritten."
)
writeLines(summary, file.path(result_dir, "expanded_grid_boundary_sensitivity_summary.md"))
message("Wrote expanded-grid functional and dominant-attractor sensitivity audit")
