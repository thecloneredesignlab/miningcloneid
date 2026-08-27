#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)
out_root <- file.path(repo_root, "revised", "add", "round3_mandatory")
result_dir <- file.path(out_root, "results", "02_predictive_adequacy")
figure_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, name) write.table(
  x, file.path(result_dir, name), sep = "\t", row.names = FALSE,
  quote = FALSE, na = "NA"
)
rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))
mae <- function(x) mean(abs(x), na.rm = TRUE)
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}
quant <- function(x, p) unname(quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))

# Standalone in-vitro fit: all observed growth, karyotype, and flow samples.
vitro_root <- file.path(
  "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results",
  "fit_invitro_O2_buffering_500seed", "seed10"
)
required_vitro <- c(
  "invitro_growth_loglik.tsv", "invitro_ploidy_loglik.tsv",
  "invitro_flow_loglik.tsv", "invitro_flow_overlay.tsv"
)
if (!all(file.exists(file.path(vitro_root, required_vitro)))) {
  stop("Missing standalone in-vitro fit outputs under ", vitro_root)
}
growth <- read_tsv(file.path(vitro_root, "invitro_growth_loglik.tsv"))
growth <- growth[is.finite(growth$observed_growth) & is.finite(growth$predicted_growth_rate), ]
growth$raw_residual <- growth$observed_growth - growth$predicted_growth_rate
growth$absolute_residual <- abs(growth$raw_residual)
write_tsv(growth, "standalone_invitro_growth_all_samples.tsv")

growth_summary <- data.frame(
  endpoint = "Growth rate",
  n = nrow(growth),
  rmse = rmse(growth$raw_residual),
  mae = mae(growth$raw_residual),
  mean_residual = mean(growth$raw_residual),
  observed_predicted_correlation = safe_cor(
    growth$observed_growth, growth$predicted_growth_rate
  )
)

ploidy_ll <- read_tsv(file.path(vitro_root, "invitro_ploidy_loglik.tsv"))
flow_ll <- read_tsv(file.path(vitro_root, "invitro_flow_loglik.tsv"))
write_tsv(ploidy_ll, "standalone_invitro_karyotype_log_scores_all_passages.tsv")
write_tsv(flow_ll, "standalone_invitro_flow_log_scores_all_samples.tsv")

overlay <- read_tsv(file.path(vitro_root, "invitro_flow_overlay.tsv"))
flow_split <- split(overlay, overlay$sample_name)
flow_metrics <- do.call(rbind, lapply(flow_split, function(d) {
  obs <- d[d$series == "Observed", ]
  pred <- d[d$series == "Predicted", ]
  m <- merge(
    obs[, c("segment_id", "cohort", "passage_index", "oxygen_pct", "passage_id",
             "sample_name", "grid_index", "ploidy", "density")],
    pred[, c("sample_name", "grid_index", "density")],
    by = c("sample_name", "grid_index"), suffixes = c("_observed", "_predicted")
  )
  m <- m[order(m$grid_index), ]
  dx <- if (nrow(m) > 1L) median(diff(m$ploidy)) else NA_real_
  data.frame(
    segment_id = m$segment_id[[1]],
    cohort = m$cohort[[1]],
    passage_index = m$passage_index[[1]],
    oxygen_pct = m$oxygen_pct[[1]],
    passage_id = m$passage_id[[1]],
    sample_name = m$sample_name[[1]],
    n_grid = nrow(m),
    density_rmse = rmse(m$density_observed - m$density_predicted),
    density_integrated_absolute_error = sum(
      abs(m$density_observed - m$density_predicted)
    ) * dx,
    density_correlation = safe_cor(m$density_observed, m$density_predicted)
  )
}))
flow_metrics <- merge(
  flow_metrics, flow_ll[, c("sample_name", "mean_loglik")],
  by = "sample_name", all.x = TRUE
)
write_tsv(flow_metrics, "standalone_invitro_flow_curve_metrics_all_samples.tsv")

standalone_summary <- rbind(
  growth_summary,
  data.frame(
    endpoint = "Karyotype cell log score",
    n = sum(ploidy_ll$n_cells),
    rmse = NA_real_,
    mae = NA_real_,
    mean_residual = -sum(ploidy_ll$total_loglik) / sum(ploidy_ll$n_cells),
    observed_predicted_correlation = NA_real_
  ),
  data.frame(
    endpoint = "Flow-density curve",
    n = nrow(flow_metrics),
    rmse = median(flow_metrics$density_rmse),
    mae = median(flow_metrics$density_integrated_absolute_error),
    mean_residual = -mean(flow_metrics$mean_loglik),
    observed_predicted_correlation = median(
      flow_metrics$density_correlation, na.rm = TRUE
    )
  )
)
write_tsv(standalone_summary, "standalone_invitro_predictive_summary.tsv")

# Six objective-near joint winners: all tumors and terminal distributions.
frozen_root <- file.path(
  repo_root, "revised", "iteration1", "data", "Figures", "Figure5",
  "figure5_frozen_inputs"
)
selection <- read_tsv(file.path(frozen_root, "selected_results.tsv"))
selection <- selection[selection$record_type == "joint_pair_best", ]
if (nrow(selection) != 6L) stop("Expected six frozen joint winners; found ", nrow(selection))

burden_all <- list()
terminal_all <- list()
joint_growth_all <- list()
for (i in seq_len(nrow(selection))) {
  label <- selection$warmup_label[[i]]
  d <- file.path(frozen_root, "winners", label)
  b <- read_tsv(file.path(d, "invivo_burden_fit.tsv"))
  b$solution <- label
  b$objective <- selection$objective[[i]]
  b$objective_delta <- b$objective - min(selection$objective)
  burden_all[[i]] <- b

  p <- read_tsv(file.path(d, "invivo_terminal_ploidy_fit.tsv"))
  p$solution <- label
  p$objective <- selection$objective[[i]]
  p$objective_delta <- p$objective - min(selection$objective)
  terminal_all[[i]] <- p

  g <- read_tsv(file.path(d, "invitro_growth_loglik.tsv"))
  g <- unique(g[, c(
    "segment_id", "cohort", "passage_index", "oxygen_pct", "passage_id",
    "predicted_growth_rate", "observed_growth", "loglik"
  )])
  g$solution <- label
  g$objective <- selection$objective[[i]]
  joint_growth_all[[i]] <- g
}
burden <- do.call(rbind, burden_all)
terminal <- do.call(rbind, terminal_all)
joint_growth <- do.call(rbind, joint_growth_all)

burden_nonzero <- burden[
  burden$day > 0 & burden$obs_burden > 0 &
    burden$pred_burden_volume_mm3 > 0, ]
burden_nonzero$log_residual <- log(burden_nonzero$obs_burden) -
  log(burden_nonzero$pred_burden_volume_mm3)
burden_metrics <- do.call(rbind, lapply(split(burden_nonzero, burden_nonzero$solution), function(d) {
  data.frame(
    solution = d$solution[[1]],
    objective = d$objective[[1]],
    objective_delta = d$objective_delta[[1]],
    n_tumor_timepoints = nrow(d),
    n_tumors = length(unique(d$harvest)),
    log_burden_rmse = rmse(d$log_residual),
    log_burden_mae = mae(d$log_residual),
    log_burden_bias_obs_minus_pred = mean(d$log_residual),
    observed_predicted_log_correlation = safe_cor(
      log(d$obs_burden), log(d$pred_burden_volume_mm3)
    )
  )
}))
write_tsv(burden_metrics, "joint_winner_burden_metrics.tsv")

key_burden <- interaction(
  burden_nonzero$harvest, burden_nonzero$day, drop = TRUE, lex.order = TRUE
)
burden_envelope <- do.call(rbind, lapply(split(burden_nonzero, key_burden), function(d) {
  preds <- d$pred_burden_volume_mm3
  data.frame(
    harvest = d$harvest[[1]],
    cohort = d$cohort[[1]],
    dose = d$dose[[1]],
    day = d$day[[1]],
    obs_burden = d$obs_burden[[1]],
    pred_min = min(preds),
    pred_q025 = quant(preds, 0.025),
    pred_q25 = quant(preds, 0.25),
    pred_median = median(preds),
    pred_q75 = quant(preds, 0.75),
    pred_q975 = quant(preds, 0.975),
    pred_max = max(preds),
    observed_in_minmax_envelope = d$obs_burden[[1]] >= min(preds) &
      d$obs_burden[[1]] <= max(preds),
    observed_in_q025_q975_envelope = d$obs_burden[[1]] >= quant(preds, 0.025) &
      d$obs_burden[[1]] <= quant(preds, 0.975)
  )
}))
burden_envelope$median_log_residual <- log(burden_envelope$obs_burden) -
  log(burden_envelope$pred_median)
write_tsv(burden_envelope, "joint_winner_burden_objective_near_envelope.tsv")

burden_by_tumor <- do.call(rbind, lapply(split(burden_envelope, burden_envelope$harvest), function(d) {
  data.frame(
    harvest = d$harvest[[1]],
    cohort = d$cohort[[1]],
    dose = d$dose[[1]],
    n_timepoints = nrow(d),
    median_prediction_log_rmse = rmse(d$median_log_residual),
    median_prediction_log_bias = mean(d$median_log_residual),
    minmax_coverage = mean(d$observed_in_minmax_envelope),
    q025_q975_coverage = mean(d$observed_in_q025_q975_envelope)
  )
}))
write_tsv(burden_by_tumor, "joint_winner_burden_metrics_all_tumors.tsv")

terminal_key <- interaction(terminal$solution, terminal$harvest, drop = TRUE)
terminal_metrics <- do.call(rbind, lapply(split(terminal, terminal_key), function(d) {
  p <- d$pred_fraction / sum(d$pred_fraction)
  n_obs <- sum(d$obs_count)
  q <- if (n_obs > 0) d$obs_count / n_obs else rep(NA_real_, nrow(d))
  cdf_p <- cumsum(p)
  cdf_q <- cumsum(q)
  m <- 0.5 * (p + q)
  kl <- function(a, b) sum(ifelse(a > 0, a * log(a / pmax(b, .Machine$double.xmin)), 0))
  js <- 0.5 * kl(p, m) + 0.5 * kl(q, m)
  data.frame(
    solution = d$solution[[1]],
    objective = d$objective[[1]],
    objective_delta = d$objective_delta[[1]],
    harvest = d$harvest[[1]],
    cohort = d$cohort[[1]],
    dose = d$dose[[1]],
    n_observed_cells = n_obs,
    observed_mean_N = sum(d$N * q),
    predicted_mean_N = sum(d$N * p),
    mean_N_residual_obs_minus_pred = sum(d$N * q) - sum(d$N * p),
    total_variation_distance = 0.5 * sum(abs(q - p)),
    wasserstein1_chromosomes = sum(abs(cdf_q - cdf_p)),
    jensen_shannon_distance = sqrt(js)
  )
}))
write_tsv(terminal_metrics, "joint_winner_terminal_ploidy_metrics_all_tumors.tsv")

terminal_envelope <- do.call(rbind, lapply(split(terminal_metrics, terminal_metrics$harvest), function(d) {
  data.frame(
    harvest = d$harvest[[1]],
    cohort = d$cohort[[1]],
    dose = d$dose[[1]],
    n_observed_cells = d$n_observed_cells[[1]],
    observed_mean_N = d$observed_mean_N[[1]],
    predicted_mean_min = min(d$predicted_mean_N),
    predicted_mean_median = median(d$predicted_mean_N),
    predicted_mean_max = max(d$predicted_mean_N),
    observed_mean_in_minmax_envelope = d$observed_mean_N[[1]] >= min(d$predicted_mean_N) &
      d$observed_mean_N[[1]] <= max(d$predicted_mean_N),
    median_total_variation = median(d$total_variation_distance),
    median_wasserstein1 = median(d$wasserstein1_chromosomes),
    median_jensen_shannon_distance = median(d$jensen_shannon_distance)
  )
}))
write_tsv(terminal_envelope, "joint_winner_terminal_ploidy_objective_near_envelope.tsv")

joint_growth <- joint_growth[
  is.finite(joint_growth$observed_growth) & is.finite(joint_growth$predicted_growth_rate), ]
joint_growth$residual <- joint_growth$observed_growth - joint_growth$predicted_growth_rate
joint_growth_metrics <- do.call(rbind, lapply(split(joint_growth, joint_growth$solution), function(d) {
  data.frame(
    solution = d$solution[[1]],
    objective = d$objective[[1]],
    n_growth_observations = nrow(d),
    growth_rmse = rmse(d$residual),
    growth_mae = mae(d$residual),
    growth_bias_obs_minus_pred = mean(d$residual),
    observed_predicted_correlation = safe_cor(
      d$observed_growth, d$predicted_growth_rate
    )
  )
}))
write_tsv(joint_growth_metrics, "joint_winner_invitro_growth_metrics.tsv")

coverage_summary <- data.frame(
  quantity = c(
    "Burden timepoints within six-winner min-max envelope",
    "Burden timepoints within six-winner 2.5%-97.5% envelope",
    "Terminal observed means within six-winner min-max envelope"
  ),
  numerator = c(
    sum(burden_envelope$observed_in_minmax_envelope),
    sum(burden_envelope$observed_in_q025_q975_envelope),
    sum(terminal_envelope$observed_mean_in_minmax_envelope)
  ),
  denominator = c(
    nrow(burden_envelope), nrow(burden_envelope), nrow(terminal_envelope)
  )
)
coverage_summary$fraction <- coverage_summary$numerator / coverage_summary$denominator
write_tsv(coverage_summary, "objective_near_envelope_coverage_summary.tsv")

plot_predictive <- function(device, filename, width, height, res = NULL) {
  if (device == "pdf") pdf(filename, width = width, height = height, useDingbats = FALSE)
  else png(filename, width = width, height = height, units = "in", res = res)
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.3, 4.5, 2.2, 0.8), las = 1)
  x <- burden_envelope$obs_burden
  y <- burden_envelope$pred_median
  xr <- range(c(x, y, burden_envelope$pred_min, burden_envelope$pred_max), positive = TRUE)
  plot(
    x, y, log = "xy", pch = 21, bg = "#4E79A7", col = "white",
    xlab = "Observed burden (mm3)", ylab = "Median predicted burden (mm3)",
    main = "All in-vivo timepoints", xlim = xr, ylim = xr
  )
  segments(
    x, burden_envelope$pred_min, x, burden_envelope$pred_max,
    col = adjustcolor("#4E79A7", 0.35)
  )
  points(x, y, pch = 21, bg = "#4E79A7", col = "white")
  abline(0, 1, lty = 2, col = "#555555")
  mtext("Bars: min-max across six joint winners", side = 3, line = 0.15, cex = 0.72)

  x2 <- terminal_envelope$observed_mean_N
  y2 <- terminal_envelope$predicted_mean_median
  xr2 <- range(c(
    x2, y2, terminal_envelope$predicted_mean_min,
    terminal_envelope$predicted_mean_max
  ))
  plot(
    x2, y2, pch = 21,
    bg = ifelse(terminal_envelope$cohort == "2N", "#59A14F", "#E15759"),
    col = "white", cex = 1.1,
    xlab = "Observed terminal mean N", ylab = "Median predicted terminal mean N",
    main = "All terminal tumors", xlim = xr2, ylim = xr2
  )
  segments(
    x2, terminal_envelope$predicted_mean_min,
    x2, terminal_envelope$predicted_mean_max,
    col = adjustcolor("#333333", 0.35)
  )
  points(
    x2, y2, pch = 21,
    bg = ifelse(terminal_envelope$cohort == "2N", "#59A14F", "#E15759"),
    col = "white", cex = 1.1
  )
  abline(0, 1, lty = 2, col = "#555555")
  legend(
    "topleft", legend = c("2N", "4N"), pt.bg = c("#59A14F", "#E15759"),
    pch = 21, bty = "n", cex = 0.8
  )
}
plot_predictive(
  "pdf", file.path(figure_dir, "round3_full_sample_predictive_adequacy.pdf"),
  7.5, 3.6
)
plot_predictive(
  "png", file.path(figure_dir, "round3_full_sample_predictive_adequacy.png"),
  7.5, 3.6, 300
)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
tex <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Full-sample in-sample adequacy audit. The six-winner ranges summarize numerical solution multiplicity among the frozen joint cluster winners and are not confidence or posterior intervals. This is not held-out predictive validation.}",
  "\\label{tab:round3_full_sample_adequacy}",
  "\\small",
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Diagnostic & Estimate & Evaluation units \\\\",
  "\\midrule",
  sprintf(
    "Standalone in-vitro growth RMSE & %s & %d passage records \\\\",
    fmt(growth_summary$rmse), growth_summary$n
  ),
  sprintf(
    "Standalone flow median integrated absolute error & %s & %d sample curves \\\\",
    fmt(median(flow_metrics$density_integrated_absolute_error)), nrow(flow_metrics)
  ),
  sprintf(
    "Joint-winner burden log-RMSE, median (range) & %s (%s--%s) & %d tumor-timepoints \\\\",
    fmt(median(burden_metrics$log_burden_rmse)),
    fmt(min(burden_metrics$log_burden_rmse)),
    fmt(max(burden_metrics$log_burden_rmse)),
    nrow(burden_envelope)
  ),
  sprintf(
    "Burden observed within winner min--max & %s & %d/%d timepoints \\\\",
    fmt(coverage_summary$fraction[[1]]),
    coverage_summary$numerator[[1]], coverage_summary$denominator[[1]]
  ),
  sprintf(
    "Terminal mean observed within winner min--max & %s & %d/%d tumors \\\\",
    fmt(coverage_summary$fraction[[3]]),
    coverage_summary$numerator[[3]], coverage_summary$denominator[[3]]
  ),
  sprintf(
    "Terminal ploidy W$_1$, median across winners/tumors & %s & %d solution--tumor pairs \\\\",
    fmt(median(terminal_metrics$wasserstein1_chromosomes), 2),
    nrow(terminal_metrics)
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(tex, file.path(table_dir, "round3_full_sample_predictive_adequacy.tex"))

summary <- c(
  "# Full-sample predictive adequacy",
  "",
  sprintf(
    "Standalone in-vitro diagnostics include %d growth records, %d directly karyotyped cells across %d passages, and %d flow-density samples.",
    nrow(growth), sum(ploidy_ll$n_cells), nrow(ploidy_ll), nrow(flow_metrics)
  ),
  sprintf(
    "Across the six frozen joint winners, burden diagnostics include %d nonzero tumor-timepoints from %d tumors and terminal-ploidy diagnostics include %d tumors.",
    nrow(burden_envelope), length(unique(burden_envelope$harvest)),
    nrow(terminal_envelope)
  ),
  sprintf(
    "Observed burden falls inside the six-winner min-max envelope at %d/%d timepoints; observed terminal mean chromosome state falls inside the corresponding envelope for %d/%d tumors.",
    coverage_summary$numerator[[1]], coverage_summary$denominator[[1]],
    coverage_summary$numerator[[3]], coverage_summary$denominator[[3]]
  ),
  "",
  "These are in-sample predictive checks. The six-winner ranges reflect numerical solution multiplicity, not sampling uncertainty. Leave-one-tumor-out or leave-one-lineage-out validation requires refitting and is not claimed here."
)
writeLines(summary, file.path(result_dir, "full_sample_predictive_adequacy_summary.md"))
message("Wrote full-sample predictive adequacy outputs to ", result_dir)
