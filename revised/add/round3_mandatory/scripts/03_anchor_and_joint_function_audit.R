#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)
out_root <- file.path(repo_root, "revised", "add", "round3_mandatory")
result_dir <- file.path(out_root, "results", "03_joint_sensitivity")
figure_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
hpc_dir <- file.path(out_root, "hpc")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hpc_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(
  x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
get_metric <- function(tab, key) {
  x <- suppressWarnings(as.numeric(tab$value[tab$metric == key]))
  if (length(x)) x[[1]] else NA_real_
}
rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

# Standalone in-vitro anchors: global best, same-cluster near-best, other-cluster best.
vitro_base <- file.path(
  "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results",
  "fit_invitro_O2_buffering_500seed"
)
anchor_spec <- data.frame(
  seed = c(10L, 132L, 157L),
  role = c("global_best", "same_cluster_near_best", "other_cluster_best"),
  cluster = c("vt_C02", "vt_C02", "vt_C01"),
  stringsAsFactors = FALSE
)
anchor_metrics <- list()
anchor_params <- list()
for (i in seq_len(nrow(anchor_spec))) {
  seed <- anchor_spec$seed[[i]]
  d <- file.path(vitro_base, paste0("seed", seed))
  fit_summary <- read_tsv(file.path(d, "fit_summary.tsv"))
  growth <- read_tsv(file.path(d, "invitro_growth_loglik.tsv"))
  growth <- growth[
    is.finite(growth$observed_growth) & is.finite(growth$predicted_growth_rate), ]
  ploidy <- read_tsv(file.path(d, "invitro_ploidy_loglik.tsv"))
  flow <- read_tsv(file.path(d, "invitro_flow_loglik.tsv"))
  anchor_metrics[[i]] <- data.frame(
    seed = seed,
    role = anchor_spec$role[[i]],
    cluster = anchor_spec$cluster[[i]],
    objective_total = get_metric(fit_summary, "objective_total"),
    objective_delta_from_seed10 = get_metric(fit_summary, "objective_total") -
      3.85253526260594,
    growth_component = get_metric(fit_summary, "growth_loglik"),
    karyotype_component = get_metric(fit_summary, "ploidy_loglik"),
    flow_component = get_metric(fit_summary, "flow_loglik"),
    growth_rmse = rmse(growth$observed_growth - growth$predicted_growth_rate),
    growth_correlation = safe_cor(growth$observed_growth, growth$predicted_growth_rate),
    mean_karyotype_negative_log_score = -sum(ploidy$total_loglik) / sum(ploidy$n_cells),
    mean_flow_negative_log_score = -mean(flow$mean_loglik)
  )
  p <- read_tsv(file.path(d, "best_params.tsv"))
  p$seed <- seed
  p$role <- anchor_spec$role[[i]]
  p$cluster <- anchor_spec$cluster[[i]]
  anchor_params[[i]] <- p
}
anchor_metrics <- do.call(rbind, anchor_metrics)
anchor_params <- do.call(rbind, anchor_params)
write_tsv(anchor_metrics, file.path(result_dir, "invitro_anchor_fit_metrics.tsv"))

param_wide <- reshape(
  anchor_params[, c("parameter", "seed", "value")],
  idvar = "parameter", timevar = "seed", direction = "wide"
)
names(param_wide) <- sub("^value\\.", "seed", names(param_wide))
for (s in c(132L, 157L)) {
  param_wide[[paste0("ratio_seed", s, "_to_seed10")]] <-
    param_wide[[paste0("seed", s)]] / param_wide$seed10
  param_wide[[paste0("log10_ratio_seed", s, "_to_seed10")]] <-
    log10(param_wide[[paste0("ratio_seed", s, "_to_seed10")]])
}
write_tsv(param_wide, file.path(result_dir, "invitro_anchor_parameter_comparison.tsv"))

# Existing six joint winners: parameter-direction and penalty-saturation audit.
frozen_root <- file.path(
  repo_root, "revised", "iteration1", "data", "Figures", "Figure5",
  "figure5_frozen_inputs"
)
selection <- read_tsv(file.path(frozen_root, "selected_results.tsv"))
selection <- selection[selection$record_type == "joint_pair_best", ]
joint_soft <- list()
function_pairs <- list()
target_o2 <- c(0, 1, 5, 20.5)
target_N <- c(44, 88)
function_names <- c(
  "ms_rate", "proliferation_rate", "death_rate",
  "buffer_death_rate", "misseg_nonviable_rate", "boundary_dropped_rate",
  "net_growth_rate"
)

for (i in seq_len(nrow(selection))) {
  label <- selection$warmup_label[[i]]
  d <- file.path(frozen_root, "winners", label)
  soft <- read_tsv(file.path(d, "joint_soft_coupling.tsv"))
  soft$solution <- label
  soft$objective <- selection$objective[[i]]
  soft$objective_delta <- soft$objective - min(selection$objective)
  joint_soft[[i]] <- soft

  vivo <- read_tsv(file.path(d, "viz", "invivo", "functional_curve_oxygen_multi_ploidy.tsv"))
  vitro <- read_tsv(file.path(d, "viz", "invitro", "functional_curve_oxygen_multi_ploidy.tsv"))
  for (N0 in target_N) {
    for (O0 in target_o2) {
      vv <- vivo[vivo$N_ref == N0, ]
      vt <- vitro[vitro$N_ref == N0, ]
      vv <- vv[which.min(abs(vv$oxygen_pct - O0)), ]
      vt <- vt[which.min(abs(vt$oxygen_pct - O0)), ]
      for (fn in function_names) {
        vivo_val <- vv[[fn]]
        vitro_val <- vt[[fn]]
        function_pairs[[length(function_pairs) + 1L]] <- data.frame(
          solution = label,
          objective = selection$objective[[i]],
          N_ref = N0,
          target_oxygen_pct = O0,
          evaluated_oxygen_pct_vivo = vv$oxygen_pct,
          evaluated_oxygen_pct_vitro = vt$oxygen_pct,
          function_name = fn,
          vivo_value = vivo_val,
          vitro_value = vitro_val,
          vivo_minus_vitro = vivo_val - vitro_val,
          vivo_to_vitro_ratio = if (is.finite(vitro_val) && vitro_val != 0) {
            vivo_val / vitro_val
          } else {
            NA_real_
          }
        )
      }
    }
  }
}
joint_soft <- do.call(rbind, joint_soft)
function_pairs <- do.call(rbind, function_pairs)
write_tsv(joint_soft, file.path(result_dir, "joint_winner_soft_coupling_all_parameters.tsv"))
write_tsv(function_pairs, file.path(result_dir, "joint_winner_function_contrasts.tsv"))

soft_summary <- do.call(rbind, lapply(split(joint_soft, joint_soft$parameter), function(d) {
  signs <- sign(d$natural_difference_vivo_to_vitro)
  dominant_sign <- if (sum(signs > 0, na.rm = TRUE) >= sum(signs < 0, na.rm = TRUE)) 1 else -1
  data.frame(
    parameter = d$parameter[[1]],
    n_solutions = nrow(d),
    median_fold_change_vivo_to_vitro = median(d$fold_change_vivo_to_vitro, na.rm = TRUE),
    min_fold_change_vivo_to_vitro = min(d$fold_change_vivo_to_vitro, na.rm = TRUE),
    max_fold_change_vivo_to_vitro = max(d$fold_change_vivo_to_vitro, na.rm = TRUE),
    dominant_direction = if (dominant_sign > 0) "vivo_gt_vitro" else "vivo_lt_vitro",
    dominant_direction_fraction = mean(signs == dominant_sign, na.rm = TRUE),
    saturating_penalty_fraction = mean(d$penalty_region == "saturating", na.rm = TRUE),
    median_penalty_paid = median(d$penalty_paid, na.rm = TRUE)
  )
}))
write_tsv(soft_summary, file.path(result_dir, "joint_winner_parameter_direction_stability.tsv"))

function_key <- interaction(
  function_pairs$N_ref, function_pairs$target_oxygen_pct,
  function_pairs$function_name, drop = TRUE
)
function_summary <- do.call(rbind, lapply(split(function_pairs, function_key), function(d) {
  signs <- sign(d$vivo_minus_vitro)
  pos <- sum(signs > 0, na.rm = TRUE)
  neg <- sum(signs < 0, na.rm = TRUE)
  dominant <- if (pos >= neg) 1 else -1
  data.frame(
    N_ref = d$N_ref[[1]],
    target_oxygen_pct = d$target_oxygen_pct[[1]],
    function_name = d$function_name[[1]],
    n_solutions = nrow(d),
    median_vivo_value = median(d$vivo_value),
    median_vitro_value = median(d$vitro_value),
    median_vivo_minus_vitro = median(d$vivo_minus_vitro),
    median_vivo_to_vitro_ratio = median(d$vivo_to_vitro_ratio, na.rm = TRUE),
    dominant_direction = if (dominant > 0) "vivo_gt_vitro" else "vivo_lt_vitro",
    dominant_direction_fraction = mean(signs == dominant, na.rm = TRUE)
  )
}))
write_tsv(function_summary, file.path(result_dir, "joint_winner_function_direction_stability.tsv"))

# Preflight refit matrix only: no submission commands or job IDs.
invivo_seed <- as.integer(sub("^seed", "", selection$invivo_seed))
anchor_rows <- expand.grid(
  invivo_seed = invivo_seed,
  invitro_seed = c(132L, 157L),
  stringsAsFactors = FALSE
)
anchor_rows$sensitivity_axis <- "alternative_invitro_anchor"
anchor_rows$joint_soft_coupling_sigma_default <- 0.65
anchor_rows$joint_soft_coupling_welsch_c <- 0.4
anchor_rows$priority <- "P0"
anchor_rows$rationale <- ifelse(
  anchor_rows$invitro_seed == 132L,
  "same-cluster near-optimal anchor",
  "other-cluster best anchor"
)

function_rows <- expand.grid(
  invivo_seed = invivo_seed,
  invitro_seed = 10L,
  joint_soft_coupling_welsch_c = c(1, 10),
  stringsAsFactors = FALSE
)
function_rows$sensitivity_axis <- "less_saturating_welsch"
function_rows$joint_soft_coupling_sigma_default <- 0.65
function_rows$priority <- "P1"
function_rows$rationale <- ifelse(
  function_rows$joint_soft_coupling_welsch_c == 1,
  "broader near-quadratic region",
  "approximately quadratic over current delta scale"
)
task_matrix <- rbind(
  anchor_rows[, c(
    "invivo_seed", "invitro_seed", "sensitivity_axis",
    "joint_soft_coupling_sigma_default", "joint_soft_coupling_welsch_c",
    "priority", "rationale"
  )],
  function_rows[, c(
    "invivo_seed", "invitro_seed", "sensitivity_axis",
    "joint_soft_coupling_sigma_default", "joint_soft_coupling_welsch_c",
    "priority", "rationale"
  )]
)
task_matrix$configuration_id <- sprintf(
  "vi%03d_vt%03d_sig%s_c%s",
  task_matrix$invivo_seed, task_matrix$invitro_seed,
  task_matrix$joint_soft_coupling_sigma_default,
  task_matrix$joint_soft_coupling_welsch_c
)
task_matrix$submission_status <- "NOT_SUBMITTED_REQUIRES_USER_CONFIRMATION"
write_tsv(task_matrix, file.path(hpc_dir, "joint_anchor_function_sensitivity_task_matrix.tsv"))

cv_matrix <- rbind(
  data.frame(
    validation_axis = "leave_one_tumor_out",
    fold = sort(unique(basename(
      read_tsv(file.path(
        frozen_root, "winners", selection$warmup_label[[1]],
        "invivo_burden_fit.tsv"
      ))$harvest
    ))),
    requires_code_extension = TRUE
  ),
  data.frame(
    validation_axis = "leave_one_lineage_out",
    fold = c("invitro_2N_deprived_O1", "invitro_2N_deprived_O2",
             "invitro_4N_deprived_O1", "invitro_4N_deprived_O2"),
    requires_code_extension = TRUE
  )
)
cv_matrix$status <- "DESIGN_ONLY_NOT_SUBMITTED"
write_tsv(cv_matrix, file.path(hpc_dir, "heldout_validation_design_matrix.tsv"))

plot_anchor <- function(device, filename, width, height, res = NULL) {
  if (device == "pdf") pdf(filename, width = width, height = height, useDingbats = FALSE)
  else png(filename, width = width, height = height, units = "in", res = res)
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  par(mfrow = c(1, 2), mar = c(5.4, 4.4, 2.2, 0.8), las = 1)
  cols <- c("#F28E2B", "#59A14F", "#4E79A7")
  barplot(
    anchor_metrics$objective_total,
    names.arg = paste0("seed", anchor_metrics$seed, "\n", anchor_metrics$cluster),
    col = cols, ylab = "In-vitro objective", ylim = c(0, 4.8),
    main = "Anchor fit quality"
  )
  abline(h = anchor_metrics$objective_total[[1]], lty = 2, col = "#555555")
  selected_params <- c(
    "lam_max", "p_misseg", "k_o_mis", "buffer_beta",
    "buffer_n_exp", "p_wgd"
  )
  p <- param_wide[param_wide$parameter %in% selected_params, ]
  mat <- rbind(
    p$log10_ratio_seed132_to_seed10,
    p$log10_ratio_seed157_to_seed10
  )
  colnames(mat) <- p$parameter
  bp <- barplot(
    mat, beside = TRUE, col = c("#59A14F", "#4E79A7"),
    ylab = "log10 parameter ratio to seed10",
    main = "Anchor parameter displacement", las = 2, cex.names = 0.72
  )
  abline(h = 0, lty = 2, col = "#555555")
  legend(
    "topright", legend = c("seed132, same cluster", "seed157, other cluster"),
    fill = c("#59A14F", "#4E79A7"), bty = "n", cex = 0.75
  )
}
plot_anchor(
  "pdf", file.path(figure_dir, "round3_invitro_anchor_sensitivity.pdf"),
  7.5, 3.7
)
plot_anchor(
  "png", file.path(figure_dir, "round3_invitro_anchor_sensitivity.png"),
  7.5, 3.7, 300
)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
key_functions <- function_summary[
  function_summary$N_ref %in% c(44, 88) &
    function_summary$target_oxygen_pct %in% c(0, 5) &
    function_summary$function_name %in% c("ms_rate", "death_rate", "net_growth_rate"), ]
tex <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{In-vitro anchor and existing joint-winner stability audit. Seed 132 is a same-cluster near-optimal anchor and seed 157 is the best anchor from the other in-vitro solution cluster. Parameter and functional stability across the six current joint winners does not replace the pending alternative-anchor refits.}",
  "\\label{tab:round3_anchor_joint_sensitivity}",
  "\\small",
  "\\begin{tabular}{llrrr}",
  "\\toprule",
  "Anchor & Cluster & Objective & Growth RMSE & Mean karyotype NLS \\\\",
  "\\midrule",
  sprintf(
    "seed%d & %s & %s & %s & %s \\\\",
    anchor_metrics$seed,
    gsub("_", paste0("\\", "_"), anchor_metrics$cluster, fixed = TRUE),
    fmt(anchor_metrics$objective_total), fmt(anchor_metrics$growth_rmse),
    fmt(anchor_metrics$mean_karyotype_negative_log_score)
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(tex, file.path(table_dir, "round3_anchor_joint_sensitivity.tex"))

summary <- c(
  "# Anchor and joint-function sensitivity",
  "",
  sprintf(
    "Seed 132 is a same-cluster near-optimal in-vitro anchor (objective delta %.4f from seed 10), whereas seed 157 is the best anchor in the other in-vitro cluster (objective delta %.3f).",
    anchor_metrics$objective_delta_from_seed10[anchor_metrics$seed == 132],
    anchor_metrics$objective_delta_from_seed10[anchor_metrics$seed == 157]
  ),
  sprintf(
    "Across the six frozen joint winners, %d/%d split parameters have a common in-vivo-versus-in-vitro direction in all six solutions, but a median %.1f%% of active parameter penalties are in the saturating Welsch region.",
    sum(soft_summary$dominant_direction_fraction == 1),
    nrow(soft_summary),
    100 * median(soft_summary$saturating_penalty_fraction)
  ),
  "",
  paste0(
    "The existing winner ensemble therefore supports a numerical direction-stability statement ",
    "conditional on the seed-10 anchor and current Welsch penalty. It does not establish anchor ",
    "or coupling-form invariance. The P0/P1 refit configurations are written to ",
    "`hpc/joint_anchor_function_sensitivity_task_matrix.tsv` and remain unsubmitted."
  )
)
writeLines(summary, file.path(result_dir, "anchor_and_joint_function_audit_summary.md"))
message("Wrote anchor/function audit and unsubmitted HPC matrices")
