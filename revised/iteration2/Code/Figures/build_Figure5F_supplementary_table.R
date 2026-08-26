#!/usr/bin/env Rscript

# Build the active Supplementary Table supporting Figure 5F. Both the DE
# initial-population spans and optimizer-endpoint spans are descriptive
# numerical summaries, not posterior or confidence intervals.

options(stringsAsFactors = FALSE, warn = 1)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
manuscript_dir <- file.path(iteration_root, "manuscript")

ratio_summary_path <- file.path(
  figure5_dir, "figure5f_prior_optimizer_summary.tsv"
)
context_summary_path <- file.path(
  figure5_dir, "figure5f_context_initial_optimizer_summary.tsv"
)
cross_path <- file.path(figure5_dir, "figure5f_prior_optimizer_cross_family.tsv")
readiness_path <- file.path(figure5_dir, "figure5f_prior_optimizer_readiness.tsv")
selection_path <- file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv")
initial_path <- file.path(
  figure5_dir, "figure5f_de_initial_population_context_values.rds"
)
initial_config_path <- file.path(
  figure5_dir, "figure5f_de_initial_population_config.tsv"
)
optimizer_path <- file.path(figure5_dir, "figure5f_optimizer_solutions.tsv")
output_tsv <- file.path(
  manuscript_dir, "tables", "data", "supp_figure5f_prior_optimizer_values.tsv"
)
output_tex <- file.path(
  manuscript_dir, "tables", "supp_figure5f_prior_optimizer.tex"
)
dir.create(dirname(output_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_tex), recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  ratio_summary_path, context_summary_path, cross_path, readiness_path,
  selection_path, initial_path, initial_config_path, optimizer_path
)
missing <- inputs[!file.exists(inputs)]
if (length(missing)) {
  stop("Missing Figure 5F SI input(s):\n", paste(missing, collapse = "\n"))
}
read_tsv <- function(path) utils::read.delim(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
ratio_summary <- read_tsv(ratio_summary_path)
context_summary <- read_tsv(context_summary_path)
cross <- read_tsv(cross_path)
readiness <- read_tsv(readiness_path)
selection <- read_tsv(selection_path)
if (any(!as.logical(readiness$passed))) {
  stop("Figure 5F DE-initial/optimizer readiness did not pass.")
}

parameters <- c(
  "lam_max", "p_mis_base", "p_wgd", "p_misseg", "k_o_mis",
  "O2_crit", "n_O", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
families <- c("C01", "C02", "C03")
required_summary <- c(
  "parameter", "family", "context", "parameter_group", "transformation",
  "de_initial_median", "de_initial_q05", "de_initial_q95",
  "optimizer_median", "optimizer_q05", "optimizer_q95",
  "optimizer_active_bound_fraction", "n_de_initial_values",
  "n_optimizer_endpoints"
)
if (nrow(context_summary) != 84L || nrow(ratio_summary) != 42L ||
    nrow(cross) != 14L ||
    !all(required_summary %in% names(context_summary)) ||
    !setequal(context_summary$parameter, parameters) ||
    !setequal(context_summary$family, families) ||
    !setequal(context_summary$context, c("in vivo", "in vitro"))) {
  stop("Figure 5F SI-table coverage is invalid.")
}

fmt <- function(x, digits = 3L) {
  ifelse(is.na(x), "NA", formatC(x, format = "fg", digits = digits, flag = "#"))
}
fmt_percent <- function(x) {
  ifelse(is.na(x), "NA", trimws(formatC(x, format = "fg", digits = 3L)))
}
escape_tex <- function(x) {
  x <- gsub("_", paste0("\\", "_"), x, fixed = TRUE)
  x <- gsub("%", paste0("\\", "%"), x, fixed = TRUE)
  gsub("&", paste0("\\", "&"), x, fixed = TRUE)
}
interval <- function(median, q05, q95) {
  paste0(fmt(median), " [", fmt(q05), ", ", fmt(q95), "]")
}

rows <- lapply(parameters, function(parameter) {
  parameter_rows <- context_summary[
    context_summary$parameter == parameter, , drop = FALSE
  ]
  cross_row <- cross[cross$parameter == parameter, , drop = FALSE]
  family_text <- vapply(seq_along(families), function(i) {
    family <- families[[i]]
    vivo <- parameter_rows[
      parameter_rows$family == family & parameter_rows$context == "in vivo",
      , drop = FALSE
    ]
    vitro <- parameter_rows[
      parameter_rows$family == family & parameter_rows$context == "in vitro",
      , drop = FALSE
    ]
    paste0(
      "in vivo: DE init ", interval(
        vivo$de_initial_median, vivo$de_initial_q05, vivo$de_initial_q95
      ),
      "; endpoint ", interval(
        vivo$optimizer_median, vivo$optimizer_q05, vivo$optimizer_q95
      ),
      "; bound=", fmt_percent(100 * vivo$optimizer_active_bound_fraction),
      "%; in vitro: DE init ", interval(
        vitro$de_initial_median, vitro$de_initial_q05, vitro$de_initial_q95
      ),
      "; endpoint ", interval(
        vitro$optimizer_median, vitro$optimizer_q05, vitro$optimizer_q95
      ),
      "; bound=", fmt_percent(100 * vitro$optimizer_active_bound_fraction), "%"
    )
  }, character(1))
  data.frame(
    parameter = parameter,
    parameter_group = parameter_rows$parameter_group[[1L]],
    transformation = parameter_rows$transformation[[1L]],
    C01_de_initial_and_optimizer = family_text[[1L]],
    C02_de_initial_and_optimizer = family_text[[2L]],
    C03_de_initial_and_optimizer = family_text[[3L]],
    cross_family_direction = cross_row$cross_family_direction,
    direction_agreement = cross_row$cross_family_direction_agreement,
    paired_endpoint_to_initial_width90_ratio_max =
      cross_row$endpoint_to_initial_width90_ratio_max,
    max_active_bound_fraction =
      cross_row$optimizer_active_bound_fraction_max,
    interpretation_limit = cross_row$inference_limit,
    stringsAsFactors = FALSE
  )
})
table_data <- do.call(rbind, rows)
rownames(table_data) <- NULL
utils::write.table(
  table_data, output_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)

selected <- selection[as.logical(selection$selected_for_figure5f), , drop = FALSE]
selected <- selected[order(match(selected$family, families)), , drop = FALSE]
source_comments <- c(
  "% AUTO-GENERATED: do not edit numerical values by hand.",
  paste0("% Direct-context summary source: ", normalizePath(context_summary_path, mustWork = TRUE)),
  paste0("% Ratio direction source: ", normalizePath(ratio_summary_path, mustWork = TRUE)),
  paste0("% Cross-family source: ", normalizePath(cross_path, mustWork = TRUE)),
  paste0("% Selection source: ", normalizePath(selection_path, mustWork = TRUE)),
  paste0("% DE-initial population source: ", normalizePath(initial_path, mustWork = TRUE)),
  paste0("% DE-initial configuration source: ", normalizePath(initial_config_path, mustWork = TRUE)),
  paste0("% Optimizer-endpoint source: ", normalizePath(optimizer_path, mustWork = TRUE)),
  paste0(
    "% Selected in-vivo warm-start seeds: C01=",
    selected$separate_invivo_seed_number[selected$family == "C01"],
    ", C02=", selected$separate_invivo_seed_number[selected$family == "C02"],
    ", C03=", selected$separate_invivo_seed_number[selected$family == "C03"], "."
  ),
  "% All three use the common in-vitro seed-10 anchor.",
  "% C-family definition: C01, C02, and C03 use in-vivo seeds 366, 25, and 311, respectively, with the common in-vitro seed-10 anchor.",
  "% DE initialization: for each C family, seed1--seed500 were replayed with NP=400 and joint_warmup_sigmaN=0.1216 using joint_deoptim_initial_population().",
  "% Each context-specific DE-initial entry reports the natural-scale median and type-8 5th/95th percentiles of 200,000 initial population members (500 seeds x 400 members).",
  "% Member 1 of every population is the exact family-specific warm start; other members use warm-start-centered truncated normals and center-dependent feasible delta bounds.",
  "% Optimizer entries report the median and type-8 5th/95th percentiles of 500 feasible, unprojected endpoints for the same C family.",
  "% bound is calculated separately for in vivo and in vitro and gives the fraction of endpoints at that context's active parameter bound.",
  "% Cross-family direction is a sign-based audit of the paired endpoint log2(in vivo/in vitro) ratios; it does not apply the 0.8/1.2 equivalence band used in the family-specific endpoint-class figure.",
  "% Neither distribution is a posterior/confidence interval, and optimizer concentration does not prove structural identifiability."
)
body <- vapply(seq_len(nrow(table_data)), function(i) {
  row <- table_data[i, , drop = FALSE]
  paste0(
    "\\texttt{", escape_tex(row$parameter), "} & ",
    escape_tex(row$C01_de_initial_and_optimizer), " & ",
    escape_tex(row$C02_de_initial_and_optimizer), " & ",
    escape_tex(row$C03_de_initial_and_optimizer), " & ",
    escape_tex(paste0(
      row$cross_family_direction,
      "; max paired W_E/W_I=",
      fmt(row$paired_endpoint_to_initial_width90_ratio_max),
      "; max bound=", fmt_percent(100 * row$max_active_bound_fraction), "%"
    )), " \\\\"
  )
}, character(1))

tex <- c(
  source_comments,
  "\\begingroup",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{2.2pt}",
  "\\begin{longtable}{p{0.08\\textwidth}p{0.22\\textwidth}p{0.22\\textwidth}p{0.22\\textwidth}p{0.15\\textwidth}}",
  "\\caption{Context-specific differential-evolution initial-population and optimizer-endpoint distributions for the 14 paired parameters. Each C-family cell reports the in vivo and in vitro natural-scale parameter distributions separately. The DE-initial distribution comprises 200,000 members reconstructed from the 500 runs (400 members per run); the endpoint distribution comprises the corresponding 500 feasible final solutions. Entries report median [5th, 95th percentile], followed by context-specific active-bound occupancy. Cross-family direction is an auxiliary sign-based audit of the paired endpoint ratios and does not apply the 0.8/1.2 equivalence band used for the family-specific classes in Supplementary Figure~\\ref{fig:supp_joint_parameter_stability}. These numerical-search distributions are not posterior or confidence intervals.}\\label{tab:figure5f_prior_optimizer}\\\\",
  "\\toprule",
  "Parameter & C01: in vivo; in vitro & C02: in vivo; in vitro & C03: in vivo; in vitro & Directional/search audit \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\toprule",
  "Parameter & C01: in vivo; in vitro & C02: in vivo; in vitro & C03: in vivo; in vitro & Directional/search audit \\\\",
  "\\midrule",
  "\\endhead",
  body,
  "\\bottomrule",
  "\\end{longtable}",
  "\\endgroup"
)
writeLines(tex, output_tex, useBytes = TRUE)

cat("Saved Figure 5F Supplementary Table:\n", output_tex, "\n", output_tsv, "\n")
