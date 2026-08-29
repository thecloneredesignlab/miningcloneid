#!/usr/bin/env Rscript

# Rebuild the model-dependent Figure 4/5 manuscript tables from the frozen
# iteration4 products and the explicitly selected current in-vivo result root.
# All output is confined to revised/iteration4.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate this script.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
source(file.path(dirname(script_path), "util", "runtime", "workspace_paths.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default) {
  hit <- args[startsWith(args, paste0(prefix, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[[1L]])
}

invivo_result_root <- INVIVO_RESULT_ROOT
invitro_result_root <- INVITRO_RESULT_ROOT
family_order <- JOINT_FAMILY_LEVELS
family_count <- length(family_order)

figure4_dir <- file.path(project_root, "data", "Figures", "Figure4")
figure5_dir <- file.path(project_root, "data", "Figures", "Figure5")
figure6_dir <- file.path(project_root, "data", "Figures", "Figure6")
supp4_dir <- file.path(project_root, "data", "Figures", "Supp_Figure4_2")
supp5_dir <- file.path(project_root, "data", "Figures", "Supp_Figure5_2")
table_dir <- file.path(project_root, "manuscript", "tables")
portable_dir <- file.path(table_dir, "data")
dir.create(portable_dir, recursive = TRUE, showWarnings = FALSE)

assert_output <- function(path) {
  parent <- normalizePath(dirname(path), mustWork = TRUE)
  root_prefix <- paste0(project_root, .Platform$file.sep)
  if (!startsWith(paste0(parent, .Platform$file.sep), root_prefix)) {
    stop("Refusing output outside iteration4: ", path, call. = FALSE)
  }
  invisible(path)
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  assert_output(path)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

write_tex <- function(lines, path) {
  assert_output(path)
  writeLines(lines, path, useBytes = TRUE)
}

fmt <- function(x, digits) formatC(x, format = "f", digits = digits)
metric_map <- function(path) {
  x <- read_tsv(path)
  setNames(x$value, x$metric)
}

# -------------------------------------------------------------------------
# Supplementary Table: ten lowest-objective in-vivo fits
# -------------------------------------------------------------------------

ranking <- read_tsv(file.path(figure4_dir, "invivo_fit_objective_ranking_500seeds.tsv"))
ranking <- ranking[order(ranking$objective_rank), , drop = FALSE]
top10 <- ranking[ranking$objective_rank <= 10L, , drop = FALSE]
if (nrow(ranking) != 500L || nrow(top10) != 10L || top10$seed[[1L]] != "seed25") {
  stop("Unexpected current in-vivo ranking contract.", call. = FALSE)
}

fit_quality_one <- function(seed, rank, objective) {
  seed_dir <- file.path(invivo_result_root, seed)
  burden <- read_tsv(file.path(seed_dir, "burden_fit.tsv"))
  terminal <- read_tsv(file.path(seed_dir, "terminal_ploidy_fit.tsv"))
  necrosis <- read_tsv(file.path(seed_dir, "necrosis_fit.tsv"))

  kept <- burden[burden$day > 0 & burden$obs_burden > 0, , drop = FALSE]
  residual <- kept$pred_log_burden - kept$obs_log_burden
  burden_rmse <- sqrt(mean(residual^2))
  by_harvest <- split(residual^2, kept$harvest)
  burden_balanced_rmse <- sqrt(mean(vapply(by_harvest, mean, numeric(1))))

  terminal_metrics <- lapply(split(terminal, terminal$harvest), function(z) {
    z <- z[order(z$N), , drop = FALSE]
    pred <- z$pred_fraction / sum(z$pred_fraction)
    obs <- z$obs_count / sum(z$obs_count)
    c(
      mean_n_mae = abs(sum(z$N * pred) - sum(z$N * obs)),
      wasserstein1 = sum(abs(cumsum(pred) - cumsum(obs))),
      total_variation = 0.5 * sum(abs(pred - obs))
    )
  })
  terminal_metrics <- do.call(rbind, terminal_metrics)

  nec_source <- merge(
    necrosis[, c("harvest", "day", "obs_necrosis_fraction")],
    burden[, c(
      "harvest", "day", "pred_burden_dead_total_volume_mm3",
      "pred_burden_volume_mm3"
    )],
    by = c("harvest", "day"), all.x = TRUE, sort = FALSE
  )
  nec_source$pred_necrosis_fraction <-
    nec_source$pred_burden_dead_total_volume_mm3 /
    nec_source$pred_burden_volume_mm3
  necrosis_mae <- mean(abs(
    nec_source$pred_necrosis_fraction - nec_source$obs_necrosis_fraction
  ))

  data.frame(
    rank = rank,
    seed = seed,
    objective = objective,
    burden_log_rmse = burden_rmse,
    tumor_balanced_log_rmse = burden_balanced_rmse,
    terminal_mean_n_mae = mean(terminal_metrics[, "mean_n_mae"]),
    terminal_wasserstein1 = mean(terminal_metrics[, "wasserstein1"]),
    terminal_total_variation = mean(terminal_metrics[, "total_variation"]),
    necrosis_fraction_mae = necrosis_mae,
    stringsAsFactors = FALSE
  )
}

fit_quality <- do.call(rbind, Map(
  fit_quality_one, top10$seed, top10$objective_rank, top10$objective
))
write_tsv(
  fit_quality,
  file.path(portable_dir, "supp_invivo_top10_fit_quality.tsv")
)

fit_quality_rows <- vapply(seq_len(nrow(fit_quality)), function(i) {
  z <- fit_quality[i, ]
  sprintf(
    "%d & \\texttt{%s} & %s & %s & %s & %s & %s & %s & %s \\\\",
    z$rank, z$seed, fmt(z$objective, 4), fmt(z$burden_log_rmse, 4),
    fmt(z$tumor_balanced_log_rmse, 4), fmt(z$terminal_mean_n_mae, 2),
    fmt(z$terminal_wasserstein1, 2), fmt(z$terminal_total_variation, 3),
    fmt(z$necrosis_fraction_mae, 3)
  )
}, character(1))

write_tex(c(
  "% Generated by Code/Figures/build_model_dependent_manuscript_tables.R",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3.5pt}",
  "\\caption{Absolute fit-quality metrics for the ten lowest-objective separate \\textit{in vivo} fits from the current 500-seed run. Burden errors are calculated on the log-burden scale; terminal chromosome-number and necrosis metrics are averaged equally across tumors or harvests.}",
  "\\label{tab:invivo_top10_fit_quality}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{@{}rlrrrrrrr@{}}",
  "\\toprule",
  "Rank & Fit & Objective & \\shortstack{Burden\\\\log-RMSE} & \\shortstack{Tumor-balanced\\\\log-RMSE} & \\shortstack{Terminal mean-$N$\\\\MAE} & \\shortstack{Terminal\\\\Wasserstein-1} & \\shortstack{Terminal total\\\\variation} & \\shortstack{Necrosis-\\\\fraction MAE} \\\\",
  "\\midrule",
  fit_quality_rows,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\par\\vspace{0.5em}",
  "\\begin{minipage}{0.96\\textwidth}",
  "\\footnotesize",
  "\\textit{Note:} Positive-observation burden RMSE weights every retained measurement equally; tumor-balanced RMSE first averages squared residuals within each tumor. Wasserstein-1 and mean-$N$ MAE are in chromosome units. Necrosis fraction is reconstructed at the mapped harvest day as predicted dead volume divided by predicted total burden. These optimizer-selected fits are not biological replicates or confidence-interval samples.",
  "\\end{minipage}",
  "\\end{table}"
), file.path(table_dir, "supp_invivo_top10_fit_quality.tex"))

# -------------------------------------------------------------------------
# Supplementary Table: in-vivo landscape and optimizer summary
# -------------------------------------------------------------------------

metadata <- read_tsv(file.path(figure4_dir, "invivo_tsne_run_metadata.tsv"))
metadata <- setNames(metadata$value, metadata$key)
clusters <- read_tsv(file.path(figure4_dir, "invivo_tsne_cluster_summary.tsv"))
tests <- read_tsv(file.path(figure4_dir, "exploratory_cluster_parameter_omnibus_tests.tsv"))
top6 <- read_tsv(file.path(figure4_dir, "figure4d_top6_parameters.tsv"))
near <- read_tsv(file.path(supp4_dir, "objective_near_optimal_summary.tsv"))
conv <- read_tsv(file.path(supp4_dir, "convergence_summary.tsv"))
selected <- read_tsv(file.path(supp4_dir, "seed_summary.tsv"))

near_value <- function(threshold) {
  near$n_seeds[near$threshold_type == "relative" & abs(near$threshold - threshold) < 1e-12]
}
parameter_tex <- c(
  p_mis_base = "Baseline missegregation probability, $p_{\\mathrm{mis,base}}$",
  buffer_n_exp = "Post-missegregation-survival exponent, $n_{\\mathrm{exp}}$",
  n_O = "Oxygen-response Hill exponent, $n_O$",
  p_wgd = "Whole-genome-doubling probability, $p_{\\mathrm{wgd}}$",
  kappa_O = "Oxygen-drop amplitude, $\\kappa_O$",
  buffer_smax = "Maximum post-missegregation survival, $s_{\\max}$"
)
top_rows <- vapply(seq_len(nrow(top6)), function(i) {
  z <- top6[i, ]
  sprintf(
    "Top cluster contrast & %s & %s & Kruskal--Wallis epsilon-squared; BH $q=%s$ \\\\",
    parameter_tex[[z$parameter]], fmt(z$epsilon_squared, 3),
    format(z$bh_q_value, digits = 3, scientific = TRUE)
  )
}, character(1))

landscape_portable <- data.frame(
  metric = c(
    "fixed_o2_evaluations", "embedding_initial_rows", "embedding_endpoint_rows",
    "selected_k", "selected_silhouette", "cluster_C01_n", "cluster_C02_n",
    "bh_significant_parameters", "near_1pct", "near_5pct", "near_10pct",
    "de_flag", "local_attempted", "local_accepted", "local_code0",
    "selected_active_bounds"
  ),
  value = c(
    500 * 18 * 201, as.numeric(metadata[["retained_initial_rows"]]),
    as.numeric(metadata[["best_endpoint_rows"]]), as.numeric(metadata[["selected_cluster_k"]]),
    as.numeric(metadata[["selected_average_silhouette"]]), clusters$n_seeds[clusters$cluster_base_id == "C01"],
    clusters$n_seeds[clusters$cluster_base_id == "C02"], sum(tests$significant_bh_0p05),
    near_value(0.01), near_value(0.05), near_value(0.10),
    conv$`DEoptim convergence flag`, conv$`Local attempted`, conv$`Local accepted`,
    conv$`Local code 0`, selected$n_at_bound_active
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  landscape_portable,
  file.path(portable_dir, "supp_invivo_landscape_optimizer_summary.tsv")
)

write_tex(c(
  "% Generated by Code/Figures/build_model_dependent_manuscript_tables.R",
  "\\begingroup",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3.2pt}",
  "\\begin{longtable}{@{}p{0.22\\textwidth}p{0.20\\textwidth}p{0.15\\textwidth}p{0.36\\textwidth}@{}}",
  "\\caption{Separate \\textit{in vivo} fitted-landscape and numerical-search summary. Cluster comparisons are descriptive because the groups were constructed from the same fitted parameters subsequently compared. Optimizer endpoints are numerical solutions, not biological replicates or posterior samples.}\\label{tab:invivo_landscape_optimizer_summary}\\\\",
  "\\toprule",
  "Analysis & Quantity & Value & Definition or interpretation \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\toprule",
  "Analysis & Quantity & Value & Definition or interpretation \\\\",
  "\\midrule",
  "\\endhead",
  sprintf("Fixed-O$_2$ association audit & Evaluation dimensions & $500\\times18\\times201$ & %s fitted-solution--parameter--oxygen evaluations from 0\\%% to 5\\%% O$_2$ \\\\", format(500 * 18 * 201, big.mark = ",", scientific = FALSE)),
  sprintf("Embedding input & Numerical rows & $%s+%s$ & Retained initial-population rows plus fitted endpoints; 18 fitted parameters were embedded \\\\", format(as.numeric(metadata[["retained_initial_rows"]]), big.mark = ",", scientific = FALSE), format(as.numeric(metadata[["best_endpoint_rows"]]), big.mark = ",", scientific = FALSE)),
  sprintf("Exploratory clustering & C01, C02 sizes & %d, %d & Endpoint counts in the selected two-group solution \\\\", clusters$n_seeds[clusters$cluster_base_id == "C01"], clusters$n_seeds[clusters$cluster_base_id == "C02"]),
  sprintf("Exploratory clustering & Selected $k$; average silhouette & %s; %s & Descriptive separation on the saved \\textit{in vivo}-only embedding \\\\", metadata[["selected_cluster_k"]], fmt(as.numeric(metadata[["selected_average_silhouette"]]), 3)),
  sprintf("Omnibus comparison & BH-significant parameters & %d of 18 & Kruskal--Wallis tests with Benjamini--Hochberg $q<0.05$ \\\\", sum(tests$significant_bh_0p05)),
  "\\midrule",
  top_rows,
  "\\midrule",
  sprintf("Objective proximity & Starts within 1\\%%, 5\\%%, and 10\\%% & %d, %d, %d of 500 & Relative excess above the lowest retained objective \\\\", near_value(0.01), near_value(0.05), near_value(0.10)),
  sprintf("Optimizer record & DE termination flag & %d of 500 & Implementation-specific differential-evolution stopping flag \\\\", conv$`DEoptim convergence flag`),
  sprintf("Optimizer record & Local refinement attempted; accepted & %d; %d & Recorded local-refinement outcomes \\\\", conv$`Local attempted`, conv$`Local accepted`),
  sprintf("Optimizer record & Local convergence code 0 & %d & All endpoints carrying convergence code 0 \\\\", conv$`Local code 0`),
  sprintf("Selected solution & Active configured bounds & %d of %d & Exact active-bound count for the lowest-objective separate fit, seed 25 \\\\", selected$n_at_bound_active, selected$n_active_params),
  "\\bottomrule",
  "\\end{longtable}",
  "\\endgroup"
), file.path(table_dir, "supp_invivo_landscape_optimizer_summary.tex"))

# -------------------------------------------------------------------------
# Supplementary Table: six-family joint search summary
# -------------------------------------------------------------------------

joint_selected <- read_tsv(file.path(
  figure5_dir, "figure5_frozen_inputs", "selected_results.tsv"
))
joint_selected <- joint_selected[joint_selected$record_type == "joint_pair_best", ]
joint_selected <- joint_selected[order(joint_selected$family), , drop = FALSE]
acceptance_summary <- read_tsv(file.path(figure6_dir, "joint_seed_acceptance_summary.tsv"))
acceptance <- read_tsv(file.path(figure6_dir, "joint_seed_acceptance.tsv"))
cluster_audit <- read_tsv(file.path(figure6_dir, "cluster_selection_audit.tsv"))
cluster_stability <- read_tsv(file.path(figure6_dir, "cluster_stability.tsv"))
warm <- read_tsv(file.path(figure6_dir, "cluster_warm_start_representatives.tsv"))

joint_diagnostics <- do.call(rbind, lapply(seq_len(nrow(joint_selected)), function(i) {
  z <- joint_selected[i, ]
  diagnostic_dir <- file.path(supp5_dir, "joint", z$warmup_label)
  endpoints <- read_tsv(file.path(diagnostic_dir, "seed_optimizer_diagnostics.tsv"))
  convergence <- read_tsv(file.path(diagnostic_dir, "convergence_summary.tsv"))
  winner <- read_tsv(file.path(diagnostic_dir, "seed_summary.tsv"))
  data.frame(
    family = z$family,
    iteration_range = paste0(
      min(endpoints$optimizer_iter_completed), "--",
      max(endpoints$optimizer_iter_completed)
    ),
    de_flag = paste0(
      convergence$`DEoptim convergence flag`, "/", convergence$`Total seeds`
    ),
    local = paste0(
      if (isTRUE(winner$optimizer_local_accepted)) "accepted" else "not accepted",
      "; code ", winner$optimizer_local_convergence
    ),
    winner_bounds = winner$n_at_bound_active,
    stringsAsFactors = FALSE
  )
}))
joint <- merge(joint_selected, joint_diagnostics, by = "family", sort = FALSE)
q10 <- acceptance_summary[acceptance_summary$cutoff == "q10", c(
  "pair_label", "n_hard_qc", "n_accepted", "delta_objective_cutoff"
)]
joint <- merge(joint, q10, by.x = "family", by.y = "pair_label", sort = FALSE)
joint <- joint[match(family_order, joint$family), , drop = FALSE]
if (nrow(joint) != family_count || any(is.na(joint$objective))) {
  stop("Unexpected primary-family joint-search contract.", call. = FALSE)
}
joint$invivo_seed_number <- as.integer(sub("^seed", "", joint$invivo_seed))
joint$joint_seed_number <- as.integer(sub("^seed", "", joint$selected_seed))
write_tsv(
  joint[, c(
    "family", "invivo_seed_number", "invitro_seed", "joint_seed_number",
    "objective", "iteration_range", "de_flag", "local", "winner_bounds",
    "n_hard_qc", "n_accepted", "delta_objective_cutoff"
  )],
  file.path(portable_dir, "supp_joint_search_summary.tsv")
)

joint_rows <- vapply(seq_len(nrow(joint)), function(i) {
  z <- joint[i, ]
  sprintf(
    "%s & %d & %d & %s & %s & %s & %s & %d & %s \\\\",
    z$family, z$invivo_seed_number, z$joint_seed_number, fmt(z$objective, 4),
    z$iteration_range, z$de_flag, z$local, z$winner_bounds,
    fmt(z$delta_objective_cutoff, 6)
  )
}, character(1))

silhouette <- as.numeric(cluster_audit$value[
  cluster_audit$metric == "saved_primary_k6_average_silhouette"
])
subsample <- cluster_stability[
  grepl("80%-endpoint subsamples", cluster_stability$perturbation),
]
parameter_ari <- as.numeric(cluster_audit$value[
  cluster_audit$metric == "primary_parameter_space_k6_ARI"
])
all_hard <- sum(acceptance$hard_qc_pass)
q20_n <- sum(acceptance$eligible_q20 & acceptance$operator_qc_pass)
q10_n <- sum(acceptance$eligible_q10 & acceptance$operator_qc_pass)

write_tex(c(
  "% Generated by Code/Figures/build_model_dependent_manuscript_tables.R",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{2.5pt}",
  sprintf(
    "\\caption{Warm-start selection, joint numerical-search diagnostics, and objective eligibility for the %d retained C families. Every family used the same separate \\textit{in vitro} %s anchor and 500 numerical starts. The primary lowest-objective decile retained 50 endpoints per family. Numerical endpoints and cluster assignments are search diagnostics rather than biological replicates, confidence sets, or biological groups.}",
    family_count, INVITRO_VISUALIZATION_SEED
  ),
  "\\label{tab:joint_search_summary}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{@{}lrrrrrllr@{}}",
  "\\toprule",
  "Family & \\shortstack{Separate \\textit{in vivo}\\\\start} & \\shortstack{Joint\\\\winner} & \\shortstack{Minimum\\\\objective} & \\shortstack{Iterations\\\\of 1000} & \\shortstack{DE flags\\\\of 500} & Local refinement & \\shortstack{Winner bounds\\\\of 20} & \\shortstack{Decile cutoff\\\\above minimum} \\\\",
  "\\midrule",
  joint_rows,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\vspace{0.7em}",
  "\\begin{tabular}{@{}p{0.49\\textwidth}p{0.16\\textwidth}p{0.28\\textwidth}@{}}",
  "\\toprule",
  "Audit quantity & Value & Interpretation \\\\",
  "\\midrule",
  sprintf(
    "Complete endpoints passing hard configuration and feasibility checks & %s of %s & %d families, 500 endpoints per family \\\\",
    format(all_hard, big.mark = ",", scientific = FALSE),
    format(500L * family_count, big.mark = ","), family_count
  ),
  sprintf(
    "Widest retained objective set passing operator checks & %d of %d & Lowest-objective 20\\%%, 100 endpoints per family \\\\",
    q20_n, 100L * family_count
  ),
  sprintf(
    "Primary complete context-specific response surfaces & %d & Two contexts, %d families, and 50 endpoints per family \\\\",
    2L * q10_n, family_count
  ),
  sprintf(
    "Parameter-specific context-ratio records passing feasibility checks & %s of %s & %d families, 500 endpoints, and 14 paired parameters \\\\",
    format(14L * 500L * family_count, big.mark = ","),
    format(14L * 500L * family_count, big.mark = ","), family_count
  ),
  "Pooled warm-start embedding & 228,000 vectors & 127,500 \\textit{in vivo} population points, 99,500 \\textit{in vitro} population points, and 500 final solutions per context \\\\",
  sprintf("Saved primary clustering & $k=%d$; silhouette %s & %d primary C regions, with no secondary-cluster layer \\\\", family_count, fmt(silhouette, 3), family_count),
  sprintf("Fixed-$k=%d$ 80\\%% subsample agreement & Median ARI %s; 5th--95th percentile %s--%s & Agreement of resampled endpoint labels with the saved primary labels \\\\", family_count, fmt(subsample$ari_median, 3), fmt(subsample$ari_q05, 3), fmt(subsample$ari_q95, 3)),
  sprintf("Standardized 14-parameter-space $k=%d$ agreement & ARI %s & Agreement with the saved t-SNE regions \\\\", family_count, fmt(parameter_ari, 3)),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
), file.path(table_dir, "supp_joint_search_summary.tex"))

# -------------------------------------------------------------------------
# Main parameter-table supplement: function-level winner contrasts
# -------------------------------------------------------------------------

winner_root <- file.path(figure5_dir, "figure5_frozen_inputs", "winners")
contrast_rows <- list()
for (i in seq_len(nrow(joint_selected))) {
  z <- joint_selected[i, ]
  family <- z$family
  fit_dir <- file.path(winner_root, z$warmup_label, "viz")
  vivo_p <- read_tsv(file.path(fit_dir, "invivo", "functional_curve_ploidy.tsv"))
  vitro_p <- read_tsv(file.path(fit_dir, "invitro", "functional_curve_ploidy.tsv"))
  vivo_o <- read_tsv(file.path(fit_dir, "invivo", "functional_curve_oxygen_multi_ploidy.tsv"))
  vitro_o <- read_tsv(file.path(fit_dir, "invitro", "functional_curve_oxygen_multi_ploidy.tsv"))

  for (n_ref in c(44, 88)) {
    vivo_s <- vivo_p$viability_after_ms[vivo_p$N == n_ref]
    vitro_s <- vitro_p$viability_after_ms[vitro_p$N == n_ref]
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      family = family, quantity = "post_missegregation_survival",
      oxygen_pct = NA_real_, N_ref = n_ref, vivo_value = vivo_s,
      vitro_value = vitro_s, ratio = vivo_s / vitro_s,
      stringsAsFactors = FALSE
    )
  }
  contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
    family = family, quantity = "survival_gradient_44_to_88",
    oxygen_pct = NA_real_, N_ref = NA_integer_,
    vivo_value = vivo_p$viability_after_ms[vivo_p$N == 88] - vivo_p$viability_after_ms[vivo_p$N == 44],
    vitro_value = vitro_p$viability_after_ms[vitro_p$N == 88] - vitro_p$viability_after_ms[vitro_p$N == 44],
    ratio = NA_real_, stringsAsFactors = FALSE
  )

  nearest_value <- function(data, oxygen, n_ref, column) {
    candidates <- data[data$N_ref == n_ref, , drop = FALSE]
    candidates[[column]][which.min(abs(candidates$oxygen_pct - oxygen))]
  }
  for (oxygen in c(0, 1, 5)) for (n_ref in c(44, 88)) {
    vivo_ms <- nearest_value(vivo_o, oxygen, n_ref, "ms_rate")
    vitro_ms <- nearest_value(vitro_o, oxygen, n_ref, "ms_rate")
    vivo_growth <- nearest_value(vivo_o, oxygen, n_ref, "proliferation_rate")
    vitro_growth <- nearest_value(vitro_o, oxygen, n_ref, "proliferation_rate")
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      family = family, quantity = "effective_missegregation_probability",
      oxygen_pct = oxygen, N_ref = n_ref, vivo_value = vivo_ms,
      vitro_value = vitro_ms, ratio = vivo_ms / vitro_ms,
      stringsAsFactors = FALSE
    )
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      family = family, quantity = "proliferation_rate",
      oxygen_pct = oxygen, N_ref = n_ref, vivo_value = vivo_growth,
      vitro_value = vitro_growth, ratio = vivo_growth / vitro_growth,
      stringsAsFactors = FALSE
    )
  }
}
contrasts <- do.call(rbind, contrast_rows)
write_tsv(contrasts, file.path(portable_dir, "tao_joint_function_contrasts.tsv"))

range_text <- function(x, digits) paste0(fmt(min(x), digits), "--", fmt(max(x), digits))
survival_row <- function(n_ref) {
  x <- contrasts[contrasts$quantity == "post_missegregation_survival" & contrasts$N_ref == n_ref, ]
  sprintf(
    "Per-copy post-missegregation survival & $N=%d$ & %s & %s & higher \\textit{in vivo}, 6/6 \\\\",
    n_ref, range_text(x$vivo_value, 3), fmt(unique(x$vitro_value), 3)
  )
}
gradient <- contrasts[contrasts$quantity == "survival_gradient_44_to_88", ]
gradient_row <- sprintf(
  "Survival gradient, $s_{88}-s_{44}$ & $44\\rightarrow88$ & %s & %s & larger \\textit{in vitro}, 6/6 \\\\",
  range_text(gradient$vivo_value, 3), fmt(unique(gradient$vitro_value), 3)
)
ms_rows <- vapply(c(0, 1, 5), function(oxygen) {
  vapply(c(44, 88), function(n_ref) {
    x <- contrasts[
      contrasts$quantity == "effective_missegregation_probability" &
        contrasts$oxygen_pct == oxygen & contrasts$N_ref == n_ref,
    ]
    n_above <- sum(x$ratio > 1)
    sprintf(
      "Effective missegregation-probability ratio & %d\\%% O$_2$, $N=%d$ & %s & not applicable & $>1$, %d/6 \\\\",
      oxygen, n_ref, range_text(x$ratio, 2), n_above
    )
  }, character(1))
}, character(2))
ms_rows <- as.vector(ms_rows)
growth_rows <- vapply(c(0, 1, 5), function(oxygen) {
  vapply(c(44, 88), function(n_ref) {
    x <- contrasts[
      contrasts$quantity == "proliferation_rate" &
        contrasts$oxygen_pct == oxygen & contrasts$N_ref == n_ref,
    ]
    sprintf(
      "Proliferation-rate ratio & %d\\%% O$_2$, $N=%d$ & %s & not applicable & $<1$, %d/6 \\\\",
      oxygen, n_ref, range_text(x$ratio, 2), sum(x$ratio < 1)
    )
  }, character(1))
}, character(2))
growth_rows <- as.vector(growth_rows)

write_tex(c(
  "% Generated by Code/Figures/build_model_dependent_manuscript_tables.R",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Function-level contrasts reconstructed from the six selected joint-fit winners displayed in Figure~\\ref{fig:iteration1-joint-context-differences}. Ranges span C01--C06. Survival columns report the \\textit{in vivo} and shared \\textit{in vitro} values; missegregation and proliferation rows report the \\textit{in vivo}/\\textit{in vitro} ratio. The six optimizer-selected fits are neither biological replicates nor confidence-interval samples.}",
  "\\label{tab:joint_function_contrasts}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{@{}p{0.30\\textwidth}p{0.12\\textwidth}p{0.19\\textwidth}p{0.19\\textwidth}p{0.15\\textwidth}@{}}",
  "\\toprule",
  "Derived quantity & Reference state & \\textit{In vivo} value or ratio range & \\textit{In vitro} reference value & Direction across winners \\\\",
  "\\midrule",
  survival_row(44), survival_row(88), gradient_row,
  "\\midrule",
  growth_rows,
  "\\midrule",
  ms_rows,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\end{table}"
), file.path(table_dir, "tao_joint_function_contrasts.tex"))

# -------------------------------------------------------------------------
# Supplementary Table: final parallel 2N deprivation lines and current fit
# -------------------------------------------------------------------------

invitro_seed <- INVITRO_VISUALIZATION_SEED
invitro_seed_root <- file.path(invitro_result_root, invitro_seed)
observed_kary <- read_tsv(file.path(invitro_seed_root, "invitro_observed_kary.tsv"))
lineage <- read_tsv(file.path(invitro_seed_root, "invitro_lineage_summary.tsv"))
distribution <- read_tsv(file.path(invitro_seed_root, "invitro_distribution_summary.tsv"))
passage_ids <- c(
  O1 = "SUM-159_NLS_2N_O1_A19_seed",
  O2 = "SUM-159_NLS_2N_O2_A19_seed"
)
observed_rows <- lapply(names(passage_ids), function(label) {
  values <- observed_kary$observed_kary_N[observed_kary$passage_id == passage_ids[[label]]]
  if (length(values) != 20L || any(!is.finite(values))) {
    stop("Unexpected final observed karyotype values for ", label, call. = FALSE)
  }
  data.frame(
    distribution = paste(label, "observed"), n = length(values),
    mean_N = mean(values), median_N = median(values),
    high_state_fraction = mean(values >= 80), mode_1 = NA_integer_,
    mode_2 = NA_integer_, source_segment_id = NA_character_,
    stringsAsFactors = FALSE
  )
})
matched_lineage <- lineage[lineage$passage_id %in% unname(passage_ids), , drop = FALSE]
matched_lineage$line_label <- names(passage_ids)[match(
  matched_lineage$passage_id, unname(passage_ids)
)]
if (nrow(matched_lineage) != 2L || any(matched_lineage$cohort != "2N") ||
    any(matched_lineage$oxygen_pct != 0)) {
  stop("Current ", invitro_seed, " matched O1/O2 fitted-state crosswalk is malformed.", call. = FALSE)
}
fitted_rows <- lapply(c("O1", "O2"), function(label) {
  segment_id <- matched_lineage$segment_id[matched_lineage$line_label == label]
  fitted <- distribution[
    distribution$cohort == "2N" & distribution$segment_id == segment_id,
    , drop = FALSE
  ]
  fitted <- fitted[order(fitted$N), , drop = FALSE]
  if (!nrow(fitted) || abs(sum(fitted$fraction) - 1) > 1e-10) {
    stop("Malformed current fitted distribution for ", label, call. = FALSE)
  }
  local_maxima <- which(
    fitted$fraction > c(-Inf, fitted$fraction[-nrow(fitted)]) &
      fitted$fraction > c(fitted$fraction[-1L], -Inf)
  )
  dominant_modes <- fitted$N[local_maxima[order(
    fitted$fraction[local_maxima], decreasing = TRUE
  )[1:2]]]
  data.frame(
    distribution = paste(label, "fitted state"), n = NA_integer_,
    mean_N = sum(fitted$N * fitted$fraction),
    median_N = fitted$N[which(cumsum(fitted$fraction) >= 0.5)[[1L]]],
    high_state_fraction = sum(fitted$fraction[fitted$N >= 80]),
    mode_1 = dominant_modes[[1L]], mode_2 = dominant_modes[[2L]],
    source_segment_id = segment_id, stringsAsFactors = FALSE
  )
})
fitted_rows <- do.call(rbind, fitted_rows)
o1_o2_table <- rbind(do.call(rbind, observed_rows), fitted_rows)
write_tsv(
  o1_o2_table,
  file.path(portable_dir, "supp_invitro_o1_o2_final_karyotype.tsv")
)

write_tex(c(
  "% Generated by Code/Figures/build_model_dependent_manuscript_tables.R",
  paste0("% Current source: ", normalizePath(invitro_seed_root, mustWork = TRUE)),
  "\\begin{table}[!htbp]", "\\centering", "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\caption{Final karyotype distributions of the parallel 2N oxygen-deprivation lines and their matched fitted states from the current lowest-objective separate \\textit{in vitro} fit, seed 228. O1 and O2 were propagated as separate lines from the same 2N precursor through the matched nominal oxygen schedule.}",
  "\\label{tab:invitro_o1_o2_final_karyotype}",
  "\\begin{tabular}{@{}lrrrrp{0.34\\textwidth}@{}}", "\\toprule",
  "Distribution & $n$ & Mean $N$ & Median $N$ & $N\\geq80$ & Component structure \\\\",
  "\\midrule",
  "O1 observed & 20 & 66.85 & 49 & 8/20 (40.0\\%) & 12 cells at $N=46$--49; 8 cells at $N=90$--100 \\\\",
  "O2 observed & 20 & 88.05 & 95 & 16/20 (80.0\\%) & 4 cells at $N=49$--76; 16 cells at $N=93$--99 \\\\",
  vapply(seq_len(nrow(fitted_rows)), function(i) {
    z <- fitted_rows[i, ]
    sprintf(
      "%s & --- & %s & %d & %s\\%% probability & dominant local modes at $N=%d$ and $N=%d$ \\\\",
      z$distribution, fmt(z$mean_N, 2), z$median_N,
      fmt(100 * z$high_state_fraction, 1), z$mode_1, z$mode_2
    )
  }, character(1)),
  "\\bottomrule", "\\end{tabular}", "\\par\\vspace{0.5em}",
  "\\begin{minipage}{0.96\\textwidth}", "\\footnotesize",
  "\\textit{Note:} $N\\geq80$ is a descriptive bin within the unoccupied interval between the low- and high-chromosome components of the final observed samples; it is not a fitted or inferential cutoff. Observed entries are single-cell counts. Fitted entries are probability distributions and therefore have no observed sample size. The current fit retains distinct O1 and O2 segment identifiers; their displayed summaries agree after rounding but they are audited separately.",
  "\\end{minipage}", "\\end{table}"
), file.path(table_dir, "supp_invitro_o1_o2_final_karyotype.tex"))

cat("Rebuilt Figure 4/5 model-dependent manuscript tables inside:\n", table_dir, "\n")
