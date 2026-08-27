#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)
out_root <- file.path(repo_root, "revised", "add", "round3_mandatory")
result_dir <- file.path(out_root, "results", "01_o1_o2_adequacy")
figure_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

source_dir <- file.path(
  repo_root, "revised", "iteration1", "data", "Figures", "Figure3", "source_seed10"
)
obs <- read.delim(file.path(source_dir, "invitro_observed_kary.tsv"), check.names = FALSE)
fit <- read.delim(file.path(source_dir, "invitro_distribution_summary.tsv"), check.names = FALSE)

target_ids <- c(
  O1 = "SUM-159_NLS_2N_O1_A19_seed",
  O2 = "SUM-159_NLS_2N_O2_A19_seed"
)
target_obs <- obs[obs$passage_id %in% unname(target_ids), , drop = FALSE]
if (nrow(target_obs) != 40L) stop("Expected 40 O1/O2 cells; found ", nrow(target_obs))

fit_target <- fit[
  fit$cohort == "2N" & fit$passage_index == 34 & fit$oxygen_pct == 0,
  c("N", "fraction"), drop = FALSE
]
if (nrow(fit_target) == 0L) stop("No shared fitted distribution at 2N passage 34, 0% O2")
if (anyDuplicated(fit_target$N)) {
  fit_target <- aggregate(fraction ~ N, fit_target, sum)
}
fit_target$fraction <- fit_target$fraction / sum(fit_target$fraction)

threshold <- 80L
eps <- .Machine$double.xmin
support <- seq.int(
  min(c(fit_target$N, target_obs$observed_kary_N)),
  max(c(fit_target$N, target_obs$observed_kary_N))
)
q <- numeric(length(support))
q[match(fit_target$N, support)] <- fit_target$fraction
q <- q / sum(q)
p_high <- sum(q[support >= threshold])
fit_mean <- sum(support * q)
fit_median <- support[which(cumsum(q) >= 0.5)[1]]

safe_js <- function(p, q) {
  m <- 0.5 * (p + q)
  kl <- function(a, b) sum(ifelse(a > 0, a * log(a / pmax(b, eps)), 0))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

one_lineage <- function(label, passage_id) {
  x <- target_obs$observed_kary_N[target_obs$passage_id == passage_id]
  n <- length(x)
  k <- sum(x >= threshold)
  p <- tabulate(match(x, support), nbins = length(support)) / n
  cdf_p <- cumsum(p)
  cdf_q <- cumsum(q)
  exact <- binom.test(k, n, p_high)
  fitted_prob_at_obs <- q[match(x, support)]
  dev <- if (k == 0 || k == n) {
    2 * (
      if (k > 0) k * log(k / (n * p_high)) else 0
    ) + 2 * (
      if (n - k > 0) (n - k) * log((n - k) / (n * (1 - p_high))) else 0
    )
  } else {
    2 * (
      k * log(k / (n * p_high)) +
        (n - k) * log((n - k) / (n * (1 - p_high)))
    )
  }
  data.frame(
    lineage = label,
    passage_id = passage_id,
    n_cells = n,
    observed_mean_N = mean(x),
    observed_median_N = median(x),
    high_state_threshold_N = threshold,
    observed_high_count = k,
    observed_high_fraction = k / n,
    observed_high_exact_ci_low = exact$conf.int[[1]],
    observed_high_exact_ci_high = exact$conf.int[[2]],
    fitted_mean_N = fit_mean,
    fitted_median_N = fit_median,
    fitted_high_fraction = p_high,
    high_fraction_residual = k / n - p_high,
    high_state_binomial_deviance = dev,
    high_state_two_sided_predictive_p = exact$p.value,
    mean_negative_log_score = -mean(log(pmax(fitted_prob_at_obs, eps))),
    total_variation_distance = 0.5 * sum(abs(p - q)),
    wasserstein1_chromosomes = sum(abs(cdf_p - cdf_q)),
    jensen_shannon_divergence_nats = safe_js(p, q),
    jensen_shannon_distance = sqrt(safe_js(p, q)),
    stringsAsFactors = FALSE
  )
}

metrics <- do.call(
  rbind,
  Map(one_lineage, names(target_ids), unname(target_ids))
)
write.table(
  metrics,
  file.path(result_dir, "o1_o2_endpoint_adequacy_metrics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

dist_rows <- do.call(rbind, lapply(names(target_ids), function(label) {
  x <- target_obs$observed_kary_N[target_obs$passage_id == target_ids[[label]]]
  p <- tabulate(match(x, support), nbins = length(support)) / length(x)
  data.frame(
    lineage = label,
    N = support,
    empirical_fraction = p,
    fitted_fraction = q,
    empirical_cdf = cumsum(p),
    fitted_cdf = cumsum(q)
  )
}))
write.table(
  dist_rows,
  file.path(result_dir, "o1_o2_endpoint_distribution_comparison.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

plot_endpoint <- function(device, filename, width, height, res = NULL) {
  if (device == "pdf") {
    pdf(filename, width = width, height = height, useDingbats = FALSE)
  } else {
    png(filename, width = width, height = height, units = "in", res = res)
  }
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.1, 0.8), las = 1)
  vals <- c(metrics$observed_high_fraction, p_high)
  bp <- barplot(
    vals,
    names.arg = c("O1", "O2", "Shared\nfit"),
    col = c("#2C7FB8", "#7FCDBB", "#F28E2B"),
    ylim = c(0, 1),
    ylab = sprintf("Fraction with N >= %d", threshold),
    main = "High-state occupancy"
  )
  arrows(
    bp[1:2], metrics$observed_high_exact_ci_low,
    bp[1:2], metrics$observed_high_exact_ci_high,
    angle = 90, code = 3, length = 0.05
  )
  text(bp, vals + 0.07, labels = sprintf("%.2f", vals), cex = 0.85)
  plot(
    NA, xlim = range(support), ylim = c(0, 1),
    xlab = "Chromosome state N", ylab = "Cumulative probability",
    main = "Endpoint distribution adequacy"
  )
  lines(support, cumsum(q), col = "#F28E2B", lwd = 2.4)
  for (i in seq_along(target_ids)) {
    d <- dist_rows[dist_rows$lineage == names(target_ids)[i], ]
    lines(
      d$N, d$empirical_cdf,
      col = c("#2C7FB8", "#238B45")[i], lwd = 2, lty = c(1, 2)[i],
      type = "s"
    )
  }
  legend(
    "bottomright",
    legend = c("Shared fit", "O1 empirical", "O2 empirical"),
    col = c("#F28E2B", "#2C7FB8", "#238B45"),
    lty = c(1, 1, 2), lwd = c(2.4, 2, 2), bty = "n", cex = 0.85
  )
}

plot_endpoint(
  "pdf", file.path(figure_dir, "round3_o1_o2_endpoint_adequacy.pdf"),
  width = 7.5, height = 3.5
)
plot_endpoint(
  "png", file.path(figure_dir, "round3_o1_o2_endpoint_adequacy.png"),
  width = 7.5, height = 3.5, res = 300
)

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
tex_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Quantitative adequacy of the shared fitted terminal distribution for the parallel O1 and O2 lineages. Distances compare each 20-cell empirical chromosome-count distribution with the same fitted passage-34 distribution. The predictive $p$ value is a two-sided exact binomial check for the number of cells with $N\\geq80$ under the fitted high-state probability; it is a model check, not a test of biological replication.}",
  "\\label{tab:round3_o1_o2_adequacy}",
  "\\small",
  "\\begin{tabular}{lrrrrrr}",
  "\\toprule",
  "Lineage & High count & Obs. high frac. & Fit high frac. & TV & W$_1$ & JS dist. \\\\",
  "\\midrule",
  sprintf(
    "%s & %d/%d & %s & %s & %s & %s & %s \\\\",
    metrics$lineage, metrics$observed_high_count, metrics$n_cells,
    fmt(metrics$observed_high_fraction), fmt(metrics$fitted_high_fraction),
    fmt(metrics$total_variation_distance),
    fmt(metrics$wasserstein1_chromosomes, 2),
    fmt(metrics$jensen_shannon_distance)
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(tex_lines, file.path(table_dir, "round3_o1_o2_endpoint_adequacy.tex"))

summary_lines <- c(
  "# O1/O2 endpoint adequacy",
  "",
  sprintf(
    "The shared fitted passage-34 distribution assigns %.3f probability to the high state (N >= %d).",
    p_high, threshold
  ),
  "",
  paste0(
    "- O1: observed ", metrics$observed_high_count[1], "/", metrics$n_cells[1],
    " high-state cells; residual ", fmt(metrics$high_fraction_residual[1]),
    "; TV ", fmt(metrics$total_variation_distance[1]),
    "; W1 ", fmt(metrics$wasserstein1_chromosomes[1], 2),
    " chromosomes; exact predictive p=", format(metrics$high_state_two_sided_predictive_p[1], digits = 3), "."
  ),
  paste0(
    "- O2: observed ", metrics$observed_high_count[2], "/", metrics$n_cells[2],
    " high-state cells; residual ", fmt(metrics$high_fraction_residual[2]),
    "; TV ", fmt(metrics$total_variation_distance[2]),
    "; W1 ", fmt(metrics$wasserstein1_chromosomes[2], 2),
    " chromosomes; exact predictive p=", format(metrics$high_state_two_sided_predictive_p[2], digits = 3), "."
  ),
  "",
  "Interpretation: the shared process supports access to both low- and high-chromosome states but does not quantitatively reproduce the lineage-specific endpoint mixture weights. These diagnostics are in-sample model checks and do not establish independent biological replication."
)
writeLines(summary_lines, file.path(result_dir, "o1_o2_endpoint_adequacy_summary.md"))
message("Wrote O1/O2 adequacy outputs to ", result_dir)

