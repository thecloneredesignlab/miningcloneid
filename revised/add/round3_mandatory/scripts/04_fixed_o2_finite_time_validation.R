#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."), mustWork = TRUE)
out_root <- file.path(repo_root, "revised", "add", "round3_mandatory")
result_dir <- file.path(out_root, "results", "04_finite_time_validation")
figure_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

source_dir <- file.path(
  "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/analysis",
  "best_fit_parameter_feature/03_dense-grid_monotonicity_classification",
  "monotonicity_classification/compare/1000D/tables"
)
read_tsv <- function(name) read.delim(
  file.path(source_dir, name), check.names = FALSE, stringsAsFactors = FALSE
)
write_tsv <- function(x, name) write.table(
  x, file.path(result_dir, name), sep = "\t", row.names = FALSE,
  quote = FALSE, na = "NA"
)

required <- c(
  "validation.tsv", "summary_by_time.tsv",
  "summary_by_time_and_gap_region.tsv", "curve_class_match_summary.tsv",
  "steady_state_support_summary.tsv", "analysis_run_arguments.tsv"
)
if (!all(file.exists(file.path(source_dir, required)))) {
  stop("Missing finite-time comparison source tables under ", source_dir)
}
validation <- read_tsv("validation.tsv")
if (!all(validation$passed)) stop("Finite-time source validation contains failures")
by_time <- read_tsv("summary_by_time.tsv")
by_gap <- read_tsv("summary_by_time_and_gap_region.tsv")
class_match <- read_tsv("curve_class_match_summary.tsv")
support <- read_tsv("steady_state_support_summary.tsv")
run_args <- read_tsv("analysis_run_arguments.tsv")

selected_days <- c(100, 300, 500, 1000)
audit <- merge(
  by_time[by_time$day %in% selected_days, c(
    "day", "n_seed_o2", "n_seed", "n_o2",
    "abs_delta_2N_minus_eigen_median",
    "abs_delta_2N_minus_eigen_fraction_le_0p05",
    "abs_delta_4N_minus_eigen_median",
    "abs_delta_4N_minus_eigen_fraction_le_0p05",
    "abs_delta_4N_minus_2N_median",
    "spearman_abs_2N_eigen_vs_gap",
    "spearman_abs_4N_eigen_vs_gap"
  )],
  class_match[class_match$day %in% selected_days, c(
    "day", "fraction_init_2N_matches_eigen",
    "fraction_init_4N_matches_eigen",
    "fraction_init_2N_matches_init_4N"
  )],
  by = "day", all = TRUE
)
write_tsv(audit, "finite_time_eigen_validation_by_day.tsv")

gap_audit <- by_gap[
  by_gap$day == 1000,
  c(
    "day", "eigen_gap_region", "n_seed_o2",
    "eigen_spectral_gap_median",
    "abs_delta_2N_minus_eigen_median",
    "abs_delta_2N_minus_eigen_fraction_le_0p05",
    "abs_delta_4N_minus_eigen_median",
    "abs_delta_4N_minus_eigen_fraction_le_0p05",
    "abs_delta_4N_minus_2N_median"
  )
]
write_tsv(gap_audit, "day1000_validation_by_spectral_gap_region.tsv")
write_tsv(validation, "source_validation.tsv")
write_tsv(support, "steady_state_support_summary.tsv")
write_tsv(run_args, "source_analysis_run_arguments.tsv")

plot_validation <- function(device, filename, width, height, res = NULL) {
  if (device == "pdf") pdf(filename, width = width, height = height, useDingbats = FALSE)
  else png(filename, width = width, height = height, units = "in", res = res)
  old <- par(no.readonly = TRUE)
  on.exit({par(old); dev.off()}, add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.3, 4.5, 2.2, 0.8), las = 1)
  plot(
    by_time$day, by_time$abs_delta_2N_minus_eigen_median,
    type = "o", log = "y", pch = 16, col = "#4E79A7", lwd = 2,
    xlab = "Finite-time horizon (days)",
    ylab = "Median absolute ploidy difference",
    main = "Convergence to eigen attractor",
    ylim = range(c(
      by_time$abs_delta_2N_minus_eigen_median,
      by_time$abs_delta_4N_minus_eigen_median
    )[c(
      by_time$abs_delta_2N_minus_eigen_median,
      by_time$abs_delta_4N_minus_eigen_median
    ) > 0])
  )
  lines(
    by_time$day, by_time$abs_delta_4N_minus_eigen_median,
    type = "o", pch = 17, col = "#E15759", lwd = 2
  )
  abline(h = 0.05, lty = 2, col = "#555555")
  legend(
    "topright", legend = c("2N start", "4N start", "0.05 threshold"),
    col = c("#4E79A7", "#E15759", "#555555"),
    pch = c(16, 17, NA), lty = c(1, 1, 2), bty = "n", cex = 0.8
  )
  plot(
    class_match$day, class_match$fraction_init_2N_matches_eigen,
    type = "o", pch = 16, col = "#4E79A7", lwd = 2,
    xlab = "Finite-time horizon (days)", ylab = "Curve-class match fraction",
    main = "Shape classification agreement", ylim = c(0, 1)
  )
  lines(
    class_match$day, class_match$fraction_init_4N_matches_eigen,
    type = "o", pch = 17, col = "#E15759", lwd = 2
  )
  lines(
    class_match$day, class_match$fraction_init_2N_matches_init_4N,
    type = "o", pch = 15, col = "#59A14F", lwd = 2
  )
  legend(
    "bottomright",
    legend = c("2N vs eigen", "4N vs eigen", "2N vs 4N"),
    col = c("#4E79A7", "#E15759", "#59A14F"),
    pch = c(16, 17, 15), lty = 1, bty = "n", cex = 0.78
  )
}
plot_validation(
  "pdf", file.path(figure_dir, "round3_fixed_o2_finite_time_validation.pdf"),
  7.5, 3.6
)
plot_validation(
  "png", file.path(figure_dir, "round3_fixed_o2_finite_time_validation.png"),
  7.5, 3.6, 300
)

day1000 <- audit[audit$day == 1000, ]
small_gap <- gap_audit[gap_audit$eigen_gap_region == "gap_lt_0p005", ]
large_gap <- gap_audit[gap_audit$eigen_gap_region == "gap_ge_0p01", ]
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
tex <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Finite-time validation of fixed-O$_2$ dominant-eigenvector summaries. Matrix-exponential trajectories were evaluated for 500 seeds, 201 O$_2$ values, both 2N and 4N initial states, and 13 time points through day 1000. Small-gap and large-gap strata use eigenvalue gaps below 0.005 and at least 0.01, respectively.}",
  "\\label{tab:round3_fixed_o2_finite_time}",
  "\\small",
  "\\begin{tabular}{lrrr}",
  "\\toprule",
  "Diagnostic & 2N start & 4N start & 2N--4N \\\\",
  "\\midrule",
  sprintf(
    "Day-1000 median absolute difference from eigen & %.3g & %.3g & %.3g \\\\",
    day1000$abs_delta_2N_minus_eigen_median,
    day1000$abs_delta_4N_minus_eigen_median,
    day1000$abs_delta_4N_minus_2N_median
  ),
  sprintf(
    "Day-1000 curve-class match fraction & %s & %s & %s \\\\",
    fmt(day1000$fraction_init_2N_matches_eigen),
    fmt(day1000$fraction_init_4N_matches_eigen),
    fmt(day1000$fraction_init_2N_matches_init_4N)
  ),
  sprintf(
    "Small-gap day-1000 median difference from eigen & %s & %s & -- \\\\",
    fmt(small_gap$abs_delta_2N_minus_eigen_median),
    fmt(small_gap$abs_delta_4N_minus_eigen_median)
  ),
  sprintf(
    "Large-gap day-1000 median difference from eigen & %.3g & %.3g & -- \\\\",
    large_gap$abs_delta_2N_minus_eigen_median,
    large_gap$abs_delta_4N_minus_eigen_median
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(tex, file.path(table_dir, "round3_fixed_o2_finite_time_validation.tex"))

summary <- c(
  "# Fixed-O2 finite-time validation",
  "",
  sprintf(
    "All %d source-validation checks passed. The comparison contains %d seed-O2 pairs per time point (%d seeds x %d O2 values), both 2N and 4N starts, and 13 horizons through day 1000.",
    nrow(validation), day1000$n_seed_o2, day1000$n_seed, day1000$n_o2
  ),
  sprintf(
    "At day 1000 the overall median absolute difference from the dominant-eigenvector ploidy is %.3g for 2N starts and %.3g for 4N starts.",
    day1000$abs_delta_2N_minus_eigen_median,
    day1000$abs_delta_4N_minus_eigen_median
  ),
  sprintf(
    "This convergence is not uniform: in the spectral-gap <0.005 stratum, the corresponding medians are %.3f and %.3f, whereas the >=0.01 stratum is essentially converged.",
    small_gap$abs_delta_2N_minus_eigen_median,
    small_gap$abs_delta_4N_minus_eigen_median
  ),
  sprintf(
    "At day 1000, finite-time curve classes match the eigen class for %.1f%% of 2N starts and %.1f%% of 4N starts; the two finite-time initial conditions agree with each other for %.1f%% of seeds.",
    100 * day1000$fraction_init_2N_matches_eigen,
    100 * day1000$fraction_init_4N_matches_eigen,
    100 * day1000$fraction_init_2N_matches_init_4N
  ),
  "",
  "Interpretation: dominant eigenvectors are valid steady-state attractor diagnostics for well-separated spectra, but they are not universal finite-time predictions. Small-gap results and curve-shape labels must be qualified or excluded from strong interpretation."
)
writeLines(summary, file.path(result_dir, "fixed_o2_finite_time_validation_summary.md"))
message("Wrote finite-time validation audit to ", result_dir)

