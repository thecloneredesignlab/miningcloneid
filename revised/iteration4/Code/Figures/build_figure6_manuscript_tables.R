#!/usr/bin/env Rscript

# Rebuild every model-dependent Figure 6 manuscript table from the fresh
# iteration4 products. Inputs outside iteration4 are never written.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot locate this script.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
source(file.path(dirname(script_path), "util", "runtime", "workspace_paths.R"))
figure4_dir <- file.path(project_root, "data", "Figures", "Figure4")
figure5_dir <- file.path(project_root, "data", "Figures", "Figure5")
figure6_dir <- file.path(project_root, "data", "Figures", "Figure6")
supp6_3_dir <- file.path(project_root, "data", "Figures", "Supp_Figure6_3")
table_dir <- file.path(project_root, "manuscript", "tables")
portable_dir <- file.path(table_dir, "data")
dir.create(portable_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing fresh Figure 6 input: ", path, call. = FALSE)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
assert_output <- function(path) {
  root_prefix <- paste0(project_root, .Platform$file.sep)
  normalized_parent <- paste0(normalizePath(dirname(path), mustWork = TRUE), .Platform$file.sep)
  if (!startsWith(normalized_parent, root_prefix)) {
    stop("Refusing output outside iteration4: ", path, call. = FALSE)
  }
  invisible(path)
}
write_tsv <- function(x, path) {
  assert_output(path)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_tex <- function(lines, path) {
  assert_output(path)
  writeLines(lines, path, useBytes = TRUE)
}
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_signed <- function(x, digits = 3) paste0(ifelse(x >= 0, "+", ""), fmt(x, digits))
mode_value <- function(x) names(sort(table(x), decreasing = TRUE))[[1L]]
family_order <- JOINT_FAMILY_LEVELS
family_count <- length(family_order)
context_order <- c("in vivo", "in vitro")

# -------------------------------------------------------------------------
# Separate-fit fixed-O2 response classes and numerical reliability
# -------------------------------------------------------------------------

vivo_class <- read_tsv(file.path(figure6_dir, "response_class_class_counts.tsv"))
vitro_class <- read_tsv(file.path(figure6_dir, "response_class_invitro_class_counts.tsv"))
vivo_reliability <- read_tsv(file.path(figure6_dir, "response_class_reliability_counts.tsv"))
vitro_reliability <- read_tsv(file.path(figure6_dir, "response_class_invitro_reliability_counts.tsv"))
vivo_seed <- read_tsv(file.path(figure6_dir, "response_class_curve_class_by_seed.tsv"))
vitro_seed <- read_tsv(file.path(figure6_dir, "response_class_invitro_curve_class_by_seed.tsv"))
representatives <- read_tsv(file.path(figure6_dir, "cluster_warm_start_representatives.tsv"))
selected <- read_tsv(file.path(
  figure5_dir, "figure5_frozen_inputs", "selected_results.tsv"
))
vitro_best <- as.integer(sub(
  "^seed", "", selected$invitro_seed[selected$record_type == "separate_invitro_global_best"]
))
if (length(vitro_best) != 1L || vitro_best != 228L) {
  stop("Unexpected current separate in-vitro winner.", call. = FALSE)
}

vivo_class$context <- "in vivo"
vitro_class$context <- "in vitro"
vivo_reliability$context <- "in vivo"
vitro_reliability$context <- "in vitro"
separate_portable <- rbind(
  transform(
    rbind(vivo_class, vitro_class), summary_type = "response_shape",
    category = smooth_curve_class
  )[, c("context", "summary_type", "category", "n_seed", "fraction_seed")],
  transform(
    rbind(vivo_reliability, vitro_reliability), summary_type = "spectral_gap",
    category = spectral_gap_class
  )[, c("context", "summary_type", "category", "n_seed", "fraction_seed")]
)
write_tsv(separate_portable, file.path(portable_dir, "supp_fixed_o2_separate_summary.tsv"))

class_labels <- c(
  complex_nonmonotone = "Complex nonmonotone",
  inverted_u_shaped = "Inverted U-shaped",
  monotone_increasing = "Monotone increasing",
  u_shaped = "U-shaped",
  approximately_flat = "Approximately flat",
  single_transition_increase_then_plateau = "Increase then plateau",
  single_transition_decrease_then_plateau = "Decrease then plateau",
  monotone_decreasing = "Monotone decreasing"
)
all_classes <- unique(c(vivo_class$smooth_curve_class, vitro_class$smooth_curve_class))
class_rows <- vapply(all_classes, function(class_id) {
  v <- vivo_class[vivo_class$smooth_curve_class == class_id, , drop = FALSE]
  t <- vitro_class[vitro_class$smooth_curve_class == class_id, , drop = FALSE]
  vn <- if (nrow(v)) v$n_seed else 0L
  tn <- if (nrow(t)) t$n_seed else 0L
  sprintf(
    "%s & %d & %s & %d & %s \\\\",
    class_labels[[class_id]], vn, fmt(100 * vn / 500, 1),
    tn, fmt(100 * tn / 500, 1)
  )
}, character(1))
reliability_rows <- vapply(c("reliable", "caution", "unreliable"), function(category) {
  v <- vivo_reliability[vivo_reliability$spectral_gap_class == category, ]
  t <- vitro_reliability[vitro_reliability$spectral_gap_class == category, ]
  sprintf(
    "%s & %d & %s & %d & %s \\\\",
    tools::toTitleCase(category), v$n_seed, fmt(100 * v$fraction_seed, 1),
    t$n_seed, fmt(100 * t$fraction_seed, 1)
  )
}, character(1))
representative_rows <- vapply(seq_len(nrow(representatives)), function(i) {
  z <- representatives[i, ]
  cls <- vivo_seed$spectral_gap_class[vivo_seed$seed_number == z$seed]
  sprintf("%s & %d & %s \\\\", z$primary_region, z$seed, cls)
}, character(1))
representative_rows <- c(
  representative_rows,
  sprintf(
    "Shared \\textit{in vitro} anchor & %d & %s \\\\",
    vitro_best,
    vitro_seed$spectral_gap_class[vitro_seed$seed_number == vitro_best]
  )
)

write_tex(c(
  "% Generated by Code/Figures/build_figure6_manuscript_tables.R",
  "\\begin{table}[!htbp]", "\\centering", "\\scriptsize",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Fixed-O$_2$ response-shape and spectral-gap classifications for the 500 fresh separate fitted solutions in each context. Reliability was classified from the fraction of the 201 oxygen-grid evaluations below the prespecified spectral-gap thresholds. These categories describe asymptotic numerical diagnostics and are not additional fitted-data hypothesis tests.}",
  "\\label{tab:fixed_o2_separate_summary}",
  "\\begin{tabular}{@{}p{0.52\\textwidth}p{0.34\\textwidth}@{}}", "\\toprule",
  "Evaluation quantity & Value \\\\", "\\midrule",
  "Separate fitted solutions and oxygen grid & 500 solutions per context; 201 values from 0\\% to 5\\% in 0.025-percentage-point increments \\\\",
  "Total seed--oxygen evaluations & 100,500 per context \\\\",
  "Unreliable criterion & At least 25\\% of grid values with spectral gap below 0.005 \\\\",
  "Caution criterion & Any grid value below 0.005, or at least 10\\% below 0.01 \\\\",
  "\\bottomrule", "\\end{tabular}", "\\vspace{0.8em}",
  "\\begin{tabular}{@{}lrrrr@{}}", "\\toprule",
  "Response-shape class & \\multicolumn{2}{c}{\\textit{In vivo}} & \\multicolumn{2}{c}{\\textit{In vitro}} \\\\",
  " & Solutions & Fraction (\\%) & Solutions & Fraction (\\%) \\\\", "\\midrule",
  class_rows, "\\bottomrule", "\\end{tabular}", "\\hspace{1.2em}",
  "\\begin{tabular}{@{}lrrrr@{}}", "\\toprule",
  "Spectral-gap class & \\multicolumn{2}{c}{\\textit{In vivo}} & \\multicolumn{2}{c}{\\textit{In vitro}} \\\\",
  " & Solutions & Fraction (\\%) & Solutions & Fraction (\\%) \\\\", "\\midrule",
  reliability_rows, "\\bottomrule", "\\end{tabular}", "\\vspace{0.8em}",
  "\\begin{tabular}{@{}lrl@{}}", "\\toprule",
  "Warm-start representative & Separate-fit seed & Spectral-gap class \\\\", "\\midrule",
  representative_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"
), file.path(table_dir, "supp_fixed_o2_separate_summary.tex"))

# -------------------------------------------------------------------------
# Joint fixed-input trajectories, inverse maps, and robustness audits
# -------------------------------------------------------------------------

acceptance <- read_tsv(file.path(figure6_dir, "joint_seed_acceptance.tsv"))

summarize_context <- function(path, context) {
  d <- read_tsv(path)
  if ("endpoint_multiplicity_q10" %in% names(d)) {
    d <- d[d$endpoint_multiplicity_q10 > 0, , drop = FALSE]
  } else {
    q10 <- acceptance[, c("pair_id", "seed_number", "eligible_q10")]
    d <- merge(
      d, q10, by = c("pair_id", "seed_number"), all.x = TRUE, sort = FALSE
    )
    if (anyNA(d$eligible_q10)) {
      stop("Cannot map context-specific robustness rows to q10 eligibility.", call. = FALSE)
    }
    d <- d[d$eligible_q10, , drop = FALSE]
  }
  rows <- lapply(split(d, d$pair_label), function(z) data.frame(
    context = context, family = z$pair_label[[1L]], n_seed = nrow(z),
    trajectory_delta = median(z$trajectory_delta_ploidy_o2_5_minus_0),
    fixed_low_delta = median(z$surface_delta_ploidy_o2_5_minus_0_cin_low),
    fixed_mid_delta = median(z$surface_delta_ploidy_o2_5_minus_0_cin_mid),
    fixed_high_delta = median(z$surface_delta_ploidy_o2_5_minus_0_cin_high),
    trajectory_direction = mode_value(z$trajectory_direction_o2_0_to_5),
    p_misseg_o2_0 = median(z$trajectory_p_misseg_o2_0),
    p_misseg_o2_5 = median(z$trajectory_p_misseg_o2_5),
    stringsAsFactors = FALSE
  ))
  do.call(rbind, rows)
}
forward <- rbind(
  summarize_context(file.path(figure6_dir, "joint_seed_biological_robustness.tsv"), "in vivo"),
  summarize_context(file.path(figure6_dir, "joint_seed_biological_robustness_invitro.tsv"), "in vitro")
)
forward <- forward[order(match(forward$context, context_order), match(forward$family, family_order)), ]
if (nrow(forward) != 12L || any(forward$n_seed != 50L)) {
  stop("Unexpected six-family paired-context forward summary.", call. = FALSE)
}
write_tsv(forward, file.path(portable_dir, "supp_joint_fixed_o2_summary.tsv"))

inverse_class <- rbind(
  transform(read_tsv(file.path(figure6_dir, "figure6_inverse_class_summary.tsv")), context = "in vivo"),
  transform(read_tsv(file.path(figure6_dir, "figure6_invitro_inverse_class_summary.tsv")), context = "in vitro")
)
inverse_anchor <- rbind(
  transform(read_tsv(file.path(figure6_dir, "figure6_inverse_ploidy4_anchor_summary.tsv")), context = "in vivo"),
  transform(read_tsv(file.path(figure6_dir, "figure6_invitro_inverse_ploidy4_anchor_summary.tsv")), context = "in vitro")
)
inverse_class <- inverse_class[order(match(inverse_class$context, context_order), match(inverse_class$pair_label, family_order)), ]
inverse_anchor <- inverse_anchor[order(match(inverse_anchor$context, context_order), match(inverse_anchor$pair_label, family_order), inverse_anchor$O2_pct), ]
write_tsv(inverse_class, file.path(portable_dir, "supp_joint_inverse_class_summary.tsv"))
write_tsv(inverse_anchor, file.path(portable_dir, "supp_joint_inverse_ploidy4_summary.tsv"))

vivo_dense_qc <- read_tsv(file.path(figure6_dir, "figure6d_dense_endpoint_qc.tsv"))
vitro_dense_qc <- read_tsv(file.path(figure6_dir, "figure6_invitro_dense_endpoint_qc.tsv"))
cutoff_vivo <- read_tsv(file.path(figure6_dir, "joint_seed_cutoff_consistency.tsv"))
cutoff_vitro <- read_tsv(file.path(figure6_dir, "joint_seed_cutoff_consistency_invitro.tsv"))
dedup_vivo <- read_tsv(file.path(figure6_dir, "joint_seed_vs_unique_parameter_robustness.tsv"))
dedup_vitro <- read_tsv(file.path(figure6_dir, "joint_seed_vs_unique_parameter_robustness_invitro.tsv"))
support_vivo <- read_tsv(file.path(figure6_dir, "joint_seed_claim_robustness.tsv"))
support_vitro <- read_tsv(file.path(figure6_dir, "joint_seed_claim_robustness_invitro.tsv"))

count_consistent <- function(d) {
  candidate <- intersect(c(
    "same_modal_result_all_cutoffs", "consistent_across_cutoffs",
    "modal_result_consistent"
  ), names(d))
  if (!length(candidate)) stop("Missing cutoff-consistency flag.", call. = FALSE)
  sum(d[[candidate[[1L]]]] %in% TRUE)
}
count_dedup_same <- function(d) {
  candidate <- intersect(c("same_modal_result", "modal_result_agrees"), names(d))
  if (!length(candidate)) stop("Missing deduplication-consistency flag.", call. = FALSE)
  sum(d[[candidate[[1L]]]] %in% TRUE)
}
logical_surface_evaluations <- 2L * family_count * 50L * 201L * 60L
distinct_dense_evaluations <-
  (nrow(vivo_dense_qc) + nrow(vitro_dense_qc)) * 496L * 201L
unique_by_panel <- c(
  tapply(vivo_dense_qc$endpoint_multiplicity_q10, vivo_dense_qc$pair_label, length),
  tapply(vitro_dense_qc$endpoint_multiplicity_q10, vitro_dense_qc$pair_label, length)
)
audit_rows <- data.frame(
  quantity = c(
    "Complete fixed-input surface per endpoint", "Logical complete-surface evaluations",
    "Dense fixed-missegregation line grid", "Distinct dense-grid operator evaluations",
    "Family-diagnostic modal results unchanged across objective cutoffs",
    "Seed-weighted versus unique-parameter comparisons with the same modal result",
    "Minimum within-set modal support", "Unique context-specific endpoints in a 50-seed ensemble"
  ),
  value = c(
    "201 x 60", format(logical_surface_evaluations, big.mark = ",", scientific = FALSE),
    "496 values", format(distinct_dense_evaluations, big.mark = ",", scientific = FALSE),
    paste0(count_consistent(cutoff_vivo) + count_consistent(cutoff_vitro), "/", nrow(cutoff_vivo) + nrow(cutoff_vitro)),
    paste0(count_dedup_same(dedup_vivo) + count_dedup_same(dedup_vitro), "/", nrow(dedup_vivo) + nrow(dedup_vitro)),
    paste0(fmt(100 * min(c(support_vivo$modal_support_fraction, support_vitro$modal_support_fraction)), 1), "\\%"),
    paste0(min(unique_by_panel), "--", max(unique_by_panel))
  ),
  scope = c(
    "Oxygen by fixed effective-missegregation grid",
    paste0("Two contexts, ", family_count, " families, 50 represented endpoints"),
    "0.005--0.500 in increments of 0.001", "Exact context-specific parameter duplicates evaluated once",
    paste0("Two contexts, ", family_count, " families, eight qualitative diagnostics"),
    paste0("Two contexts, ", family_count, " families, eight diagnostics, three cutoffs"),
    "All context--family--diagnostic--cutoff combinations",
    paste0("Range across the ", 2L * family_count, " context--family panels")
  ), stringsAsFactors = FALSE
)
write_tsv(audit_rows, file.path(portable_dir, "supp_joint_fixed_o2_audit.tsv"))

forward_rows <- vapply(seq_len(nrow(forward)), function(i) {
  z <- forward[i, ]
  sprintf(
    "\\textit{%s} & %s & $%s$ & $%s$ & $%s$ & $%s$ & %s \\\\",
    tools::toTitleCase(z$context), z$family,
    fmt_signed(z$trajectory_delta, 3), fmt_signed(z$fixed_low_delta, 3),
    fmt_signed(z$fixed_mid_delta, 3), fmt_signed(z$fixed_high_delta, 3),
    z$trajectory_direction
  )
}, character(1))
p_rows <- vapply(seq_len(nrow(forward)), function(i) {
  z <- forward[i, ]
  sprintf(
    "\\textit{%s} & %s & %s & %s & & & \\\\",
    tools::toTitleCase(z$context), z$family,
    format(z$p_misseg_o2_0, digits = 4), format(z$p_misseg_o2_5, digits = 4)
  )
}, character(1))
inverse_rows <- vapply(seq_len(nrow(inverse_class)), function(i) {
  z <- inverse_class[i, ]
  sprintf(
    "\\textit{%s} & %s & %s & %s & %s & & \\\\",
    tools::toTitleCase(z$context), z$pair_label,
    fmt(100 * z$fraction_stable_unique_inverse, 2),
    fmt(100 * z$fraction_multiple_solutions, 2),
    fmt(100 * z$fraction_no_stable_unique_inverse, 2)
  )
}, character(1))
anchor_cell <- function(z) {
  if (!nrow(z) || z$inverse_class != "stable unique inverse" || !is.finite(z$p_unique_median)) return("--")
  sprintf(
    "%s (%s--%s)", format(z$p_unique_median, digits = 4),
    format(z$p_unique_q10, digits = 4), format(z$p_unique_q90, digits = 4)
  )
}
anchor_rows <- unlist(lapply(context_order, function(context) {
  vapply(family_order, function(family) {
    z <- inverse_anchor[inverse_anchor$context == context & inverse_anchor$pair_label == family, ]
    sprintf(
      "\\textit{%s} & %s & %s & %s & %s & & \\\\",
      tools::toTitleCase(context), family,
      anchor_cell(z[z$O2_pct == 0, ]), anchor_cell(z[z$O2_pct == 1, ]),
      anchor_cell(z[z$O2_pct == 5, ])
    )
  }, character(1))
}), use.names = FALSE)
audit_tex_rows <- vapply(seq_len(nrow(audit_rows)), function(i) {
  z <- audit_rows[i, ]
  sprintf("\\multicolumn{2}{l}{%s} & %s & \\multicolumn{4}{p{0.47\\textwidth}}{%s} \\\\", z$quantity, z$value, z$scope)
}, character(1))

write_tex(c(
  "% Generated by Code/Figures/build_figure6_manuscript_tables.R",
  "\\begingroup", "\\scriptsize", "\\setlength{\\tabcolsep}{2.8pt}",
  "\\renewcommand{\\arraystretch}{0.92}",
  "\\begin{longtable}{@{}llrrrrp{0.18\\textwidth}@{}}",
  "\\caption{Context-specific joint-ensemble fixed-O$_2$ and inverse-response summaries for the lowest-objective decile of C01--C06. Forward entries report the median change in dominant mean ploidy from 0\\% to 5\\% O$_2$ across 50 represented endpoints. Inverse entries summarize the Figure~\\ref{fig:iteration1-o2-linked-model-selection}B grid.}\\label{tab:joint_fixed_o2_summary}\\\\",
  "\\toprule", "Context & Family & Unmodified & Fixed low & Fixed intermediate & Fixed high & Direction \\\\", "\\midrule", "\\endfirsthead",
  "\\toprule", "Context & Family & Unmodified & Fixed low & Fixed intermediate & Fixed high & Direction \\\\", "\\midrule", "\\endhead",
  "\\multicolumn{7}{l}{\\textit{Forward dominant-ploidy change, 5\\% minus 0\\% O$_2$}} \\\\",
  forward_rows, "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Median population-average fitted $p_{\\mathrm{miss,eff}}$ along the unmodified trajectory}} \\\\",
  "Context & Family & At 0\\% O$_2$ & At 5\\% O$_2$ & & & \\\\", p_rows, "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Inverse-grid class fractions (percent of oxygen--target-ploidy cells)}} \\\\",
  "Context & Family & Stable unique & Multiple & No stable unique & & \\\\", inverse_rows, "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Stable unique $p_{\\mathrm{miss,eff}}$ for target ploidy 4: median (10th--90th percentile)}} \\\\",
  "Context & Family & At 0\\% O$_2$ & At 1\\% O$_2$ & At 5\\% O$_2$ & & \\\\", anchor_rows, "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Design and robustness audit}} \\\\", audit_tex_rows,
  "\\bottomrule", "\\end{longtable}",
  "\\noindent\\footnotesize Numerical endpoints and exact duplicates retain their optimizer-seed multiplicities in the primary summaries. These operational ensembles are not posterior samples, confidence regions, or biological replicates.",
  "\\endgroup"
), file.path(table_dir, "supp_joint_fixed_o2_summary.tex"))

# -------------------------------------------------------------------------
# Surface topology
# -------------------------------------------------------------------------

surface_summary <- function(path, context) {
  d <- read_tsv(path)
  d <- d[d$cutoff == "q10", , drop = FALSE]
  d$regime <- ifelse(
    d$dominant_mean_ploidy_median <= 2, "low",
    ifelse(d$dominant_mean_ploidy_median < 4, "intermediate", "high")
  )
  cells <- split(d, interaction(d$O2_pct, d$effective_p_misseg, drop = TRUE))
  agreement <- mean(vapply(cells, function(z) length(unique(z$regime)) == 1L, logical(1L)))
  rows <- lapply(split(d, d$pair_label), function(z) data.frame(
    context = context, family = z$pair_label[[1L]], n_grid = nrow(z),
    low_fraction = mean(z$regime == "low"),
    intermediate_fraction = mean(z$regime == "intermediate"),
    high_fraction = mean(z$regime == "high"),
    across_family_agreement = agreement, stringsAsFactors = FALSE
  ))
  do.call(rbind, rows)
}
topology <- rbind(
  surface_summary(file.path(figure6_dir, "joint_multiseed_surface_summary.tsv"), "in vivo"),
  surface_summary(file.path(figure6_dir, "joint_multiseed_surface_summary_invitro.tsv"), "in vitro")
)
topology <- topology[order(match(topology$context, context_order), match(topology$family, family_order)), ]
write_tsv(topology, file.path(portable_dir, "supp_figure6_topology_summary.tsv"))
topology_rows <- vapply(seq_len(nrow(topology)), function(i) {
  z <- topology[i, ]
  sprintf(
    "\\textit{%s} & %s & %s & %s & %s & %s \\\\",
    tools::toTitleCase(z$context), z$family, fmt(100 * z$low_fraction, 2),
    fmt(100 * z$intermediate_fraction, 2), fmt(100 * z$high_fraction, 2),
    fmt(100 * z$across_family_agreement, 2)
  )
}, character(1))
write_tex(c(
  "% Generated by Code/Figures/build_figure6_manuscript_tables.R",
  "\\begin{table}[!htbp]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{4pt}",
  sprintf(
    "\\caption{Descriptive topology of the %d-family joint-ensemble oxygen--CIN--ploidy response. Percentages classify the seed-weighted median dominant ploidy at each of 12,060 matched cells as low ($\\leq2$), intermediate ($>2$ and $<4$), or high ($\\geq4$). Across-family agreement is the percentage of matched cells assigned to the same class in all %d families.}",
    family_count, family_count
  ),
  "\\label{tab:figure6_topology_summary}",
  "\\begin{tabular}{@{}llrrrr@{}}", "\\toprule",
  sprintf(
    "Context & Family & Low (\\%%) & Intermediate (\\%%) & High (\\%%) & %d-family agreement (\\%%) \\\\",
    family_count
  ), "\\midrule",
  topology_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"
), file.path(table_dir, "supp_figure6_topology_summary.tex"))

# -------------------------------------------------------------------------
# Dense fixed-input reference-curve features
# -------------------------------------------------------------------------

dense <- read_tsv(file.path(figure6_dir, "figure6d_fixed_p_curve_family.tsv"))
curve_direction <- function(z, tolerance = 1e-8) {
  delta <- diff(z$dominant_mean_ploidy_median[order(z$O2_pct)])
  if (all(delta >= -tolerance) && any(delta > tolerance)) "monotone increase" else
    if (all(delta <= tolerance) && any(delta < -tolerance)) "monotone decrease" else
      if (max(abs(delta)) <= tolerance) "approximately flat" else "nonmonotone"
}
sharpest_downturn <- function(z) {
  z <- z[order(z$O2_pct), ]
  z$O2_pct[-1L][which.min(diff(z$dominant_mean_ploidy_median))]
}
fixed_features <- do.call(rbind, lapply(family_order, function(family) {
  z <- dense[dense$pair_label == family, ]
  low <- z[abs(z$effective_p_misseg - 0.01) < 1e-12, ]
  mid <- z[abs(z$effective_p_misseg - 0.20) < 1e-12, ]
  high <- z[abs(z$effective_p_misseg - 0.30) < 1e-12, ]
  data.frame(
    family = family, low_input_direction = curve_direction(low),
    downturn_o2_p0p20 = sharpest_downturn(mid),
    downturn_o2_p0p30 = sharpest_downturn(high), stringsAsFactors = FALSE
  )
}))
write_tsv(fixed_features, file.path(portable_dir, "supp_figure6_fixed_input_summary.tsv"))
fixed_rows <- vapply(seq_len(nrow(fixed_features)), function(i) {
  z <- fixed_features[i, ]
  sprintf("%s & %s & %s & %s \\\\", z$family, z$low_input_direction, fmt(z$downturn_o2_p0p20, 3), fmt(z$downturn_o2_p0p30, 3))
}, character(1))
write_tex(c(
  "% Generated by Code/Figures/build_figure6_manuscript_tables.R",
  "\\begin{table}[!htbp]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Descriptive features of the six \\textit{in vivo} fixed-input reference curves in Figure~\\ref{fig:iteration1-o2-linked-model-selection}B. The sharpest sampled downturn is the oxygen coordinate at the end of the largest negative adjacent-grid change and is not a fitted biological threshold.}",
  "\\label{tab:figure6_fixed_input_summary}",
  "\\begin{tabular}{@{}p{0.13\\textwidth}p{0.23\\textwidth}p{0.27\\textwidth}p{0.27\\textwidth}@{}}", "\\toprule",
  "Family & $p_{\\mathrm{miss,eff}}=0.01$ over 0--5\\% O$_2$ & Sharpest downturn O$_2$ (\\%) at $p_{\\mathrm{miss,eff}}=0.20$ & Sharpest downturn O$_2$ (\\%) at $p_{\\mathrm{miss,eff}}=0.30$ \\\\", "\\midrule",
  fixed_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"
), file.path(table_dir, "supp_figure6_fixed_input_summary.tex"))

# -------------------------------------------------------------------------
# Weak-gap robustness
# -------------------------------------------------------------------------

weak <- read_tsv(file.path(supp6_3_dir, "supp_figure6-3_weak_gap_pair_summary.tsv"))
weak <- weak[order(match(weak$model_context, context_order), match(weak$display_label, family_order)), ]
write_tsv(weak, file.path(portable_dir, "supp_weak_gap_summary.tsv"))
weak_rows_a <- vapply(seq_len(nrow(weak)), function(i) {
  z <- weak[i, ]
  sprintf(
    "\\textit{%s} & %s & %s & %s & %s & %s & %s & %s \\\\",
    tools::toTitleCase(z$model_context), z$display_label,
    format(z$n_weak_gap_cell, big.mark = ",", scientific = FALSE),
    fmt(100 * z$fraction_weak_gap_stable_low, 1),
    fmt(100 * z$fraction_weak_gap_stable_intermediate, 1),
    fmt(100 * z$fraction_weak_gap_stable_high, 1),
    fmt(100 * z$fraction_weak_gap_mixed, 2),
    fmt(100 * z$fraction_weak_gap_majority_endpoint_local_switch, 2)
  )
}, character(1))
weak_rows_b <- vapply(seq_len(nrow(weak)), function(i) {
  z <- weak[i, ]
  sprintf(
    "\\textit{%s} & %s & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\",
    tools::toTitleCase(z$model_context), z$display_label,
    format(z$weak_gap_ploidy_spread_median, digits = 4),
    format(z$weak_gap_ploidy_spread_q90, digits = 4),
    fmt(100 * z$fraction_weak_gap_local_jump_median_ge_1, 2)
  )
}, character(1))
write_tex(c(
  "% Generated by Code/Figures/build_figure6_manuscript_tables.R",
  "\\begingroup", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}", "\\renewcommand{\\arraystretch}{0.90}",
  "\\begin{longtable}{@{}llrrrrrr@{}}",
  "\\caption{Context-specific dominant-ploidy robustness within weak-spectral-gap regions of the six displayed joint ensembles. Stable classes require at least 90\\% endpoint agreement. Local-switch and local-jump quantities compare immediately adjacent oxygen and fixed-missegregation cells. These are numerical diagnostics, not biological transition probabilities.}\\label{tab:weak_gap_summary}\\\\",
  "\\toprule", "Context & Family & Weak cells & Stable low (\\%) & Stable intermediate (\\%) & Stable high (\\%) & Mixed (\\%) & Majority local switch (\\%) \\\\", "\\midrule", "\\endfirsthead",
  "\\toprule", "Context & Family & Weak cells & Stable low (\\%) & Stable intermediate (\\%) & Stable high (\\%) & Mixed (\\%) & Majority local switch (\\%) \\\\", "\\midrule", "\\endhead",
  weak_rows_a, "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Across-endpoint spread and local jump summaries}} \\\\",
  "Context & Family & \\multicolumn{2}{c}{Median spread} & \\multicolumn{2}{c}{90th-percentile spread} & \\multicolumn{2}{c}{Local jump $\\geq1$ (\\%)} \\\\",
  weak_rows_b,
  "\\bottomrule", "\\end{longtable}", "\\endgroup"
), file.path(table_dir, "supp_weak_gap_summary.tex"))

cat("Rebuilt all Figure 6 model-dependent manuscript tables inside:\n", table_dir, "\n")
