#!/usr/bin/env Rscript

# Figure 5F input audit for the fixed-prior/generalized-posterior/optimizer-
# endpoint comparison.
#
# This script deliberately does not reinterpret optimizer endpoints as
# posterior draws. It selects one joint warm-start pair per C01/C02/C03 family
# using the objective of the original separate-in-vivo seed, then checks whether
# convergence-gated draws are available for three distinct posteriors
# conditioned on frozen C01/C02/C03 family basins.

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
  } else {
    getwd()
  }
})

workspace_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  mustWork = TRUE
)
figure5_dir <- file.path(workspace_root, "data", "Figures", "Figure5")
ranking_path <- file.path(
  workspace_root,
  "data", "Figures", "Figure4",
  "invivo_fit_objective_ranking_500seeds.tsv"
)
pair_path <- file.path(
  figure5_dir,
  "figure5_frozen_inputs", "selected_results.tsv"
)
posterior_path <- file.path(figure5_dir, "figure5f_posterior_draws.tsv")
optimizer_path <- file.path(figure5_dir, "figure5f_optimizer_solutions.tsv")
sampling_readiness_path <- file.path(
  figure5_dir, "generalized_posterior_family_conditioned",
  "figure5f_generalized_posterior_readiness.tsv"
)
sampling_configuration_path <- file.path(
  figure5_dir, "generalized_posterior_family_conditioned",
  "figure5f_generalized_posterior_config.tsv"
)
frozen_hash_sources <- c(
  generator = file.path(script_dir, "generate_Figure5F_generalized_posterior.R"),
  aggregator = file.path(script_dir, "aggregate_Figure5F_family_posterior.R"),
  combiner = file.path(script_dir, "combine_Figure5F_family_posteriors.R"),
  family_contract = file.path(
    figure5_dir, "generalized_posterior_family_conditioned",
    "figure5f_family_basin_contract.rds"
  ),
  prior = file.path(
    figure5_dir, "generalized_posterior_family_conditioned",
    "figure5f_family_conditioned_prior_draws.tsv"
  )
)
selection_output <- file.path(
  figure5_dir,
  "figure5f_selected_pair_inputs.tsv"
)
readiness_output <- file.path(
  figure5_dir,
  "figure5f_prior_posterior_readiness.tsv"
)

require_files <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing ", label, ":\n", paste(missing, collapse = "\n"))
  }
}

write_tsv <- function(x, path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  if (!startsWith(normalized, paste0(figure5_dir, .Platform$file.sep))) {
    stop("Refusing to write outside iteration2 Figure5 data: ", normalized)
  }
  utils::write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

require_files(
  c(ranking_path, pair_path),
  "Figure 5F pair-selection input"
)

ranking <- utils::read.delim(
  ranking_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pairs <- utils::read.delim(
  pair_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_ranking <- c("seed_number", "seed", "objective_rank", "objective")
required_pairs <- c(
  "record_type", "warmup_label", "invivo_seed", "invitro_seed",
  "selected_seed", "objective", "objective_invivo", "objective_invitro"
)
if (!all(required_ranking %in% names(ranking))) {
  stop(
    "Malformed separate-in-vivo objective ranking; missing: ",
    paste(setdiff(required_ranking, names(ranking)), collapse = ", ")
  )
}
if (!all(required_pairs %in% names(pairs))) {
  stop(
    "Malformed frozen pair table; missing: ",
    paste(setdiff(required_pairs, names(pairs)), collapse = ", ")
  )
}

candidate <- pairs[
  pairs$record_type == "joint_pair_best",
  required_pairs[-1L],
  drop = FALSE
]
candidate$family <- regmatches(
  candidate$warmup_label,
  regexpr("C[0-9]{2}", candidate$warmup_label)
)
candidate$subcluster <- regmatches(
  candidate$warmup_label,
  regexpr("Sc[0-9]{2}", candidate$warmup_label)
)
candidate$separate_invivo_seed_number <- as.integer(
  sub("^seed", "", candidate$invivo_seed)
)

ranking_lookup <- ranking[, c(
  "seed_number", "objective_rank", "objective"
)]
names(ranking_lookup) <- c(
  "separate_invivo_seed_number",
  "separate_invivo_objective_rank",
  "separate_invivo_objective"
)
candidate <- merge(
  candidate,
  ranking_lookup,
  by = "separate_invivo_seed_number",
  all.x = TRUE,
  sort = FALSE
)
if (nrow(candidate) != 6L ||
    anyNA(candidate$separate_invivo_objective) ||
    length(table(candidate$family)) != 3L ||
    !all(table(candidate$family) == 2L)) {
  stop("Figure 5F requires exactly two objective-matched candidates per family")
}

candidate$family_selection_rank <- ave(
  candidate$separate_invivo_objective,
  candidate$family,
  FUN = function(x) rank(x, ties.method = "first")
)
candidate$selected_for_figure5f <- candidate$family_selection_rank == 1L
candidate$selection_metric <- "original separate-in-vivo objective"
candidate$selection_rule <- paste0(
  "minimum original separate-in-vivo objective within ",
  "each C01/C02/C03 family"
)
candidate$ranking_source <- normalizePath(ranking_path, mustWork = TRUE)
candidate$ranking_source_md5 <- unname(tools::md5sum(ranking_path))
candidate$pair_source <- normalizePath(pair_path, mustWork = TRUE)
candidate$pair_source_md5 <- unname(tools::md5sum(pair_path))

candidate <- candidate[order(
  match(candidate$family, c("C01", "C02", "C03")),
  candidate$family_selection_rank
), , drop = FALSE]
rownames(candidate) <- NULL

selected <- candidate[candidate$selected_for_figure5f, , drop = FALSE]
expected_selection <- c(
  C01 = "tsne_vi_seed366_C01Sc01_vt_seed10",
  C02 = "tsne_vi_seed25_C02Sc01_vt_seed10",
  C03 = "tsne_vi_seed311_C03Sc02_vt_seed10"
)
observed_selection <- setNames(selected$warmup_label, selected$family)
if (!identical(observed_selection[names(expected_selection)], expected_selection)) {
  stop(
    "Unexpected Figure 5F selection: ",
    paste(
      names(observed_selection), observed_selection,
      sep = "=", collapse = "; "
    )
  )
}

write_tsv(candidate, selection_output)

required_posterior_columns <- c(
  "family", "warmup_label", "chain", "draw", "parameter",
  "temperature", "vivo_natural", "vitro_natural", "log2_ratio",
  "target_version"
)
required_optimizer_columns <- c(
  "family", "warmup_label", "seed_number", "parameter", "log2_ratio",
  "feasible_at_solution", "projection_applied"
)
optimizer_exists <- file.exists(optimizer_path)
optimizer_schema_valid <- FALSE
optimizer_selected_pair_coverage <- FALSE
optimizer_exactly_500_per_family_parameter <- FALSE
optimizer_numeric_and_feasibility_valid <- FALSE
if (optimizer_exists) {
  optimizer <- utils::read.delim(
    optimizer_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  optimizer_schema_valid <- all(required_optimizer_columns %in% names(optimizer))
  if (optimizer_schema_valid) {
    optimizer <- optimizer[
      optimizer$warmup_label %in% selected$warmup_label,
      , drop = FALSE
    ]
    optimizer_selected_pair_coverage <- setequal(
      unique(optimizer$warmup_label), selected$warmup_label
    ) && setequal(unique(optimizer$family), c("C01", "C02", "C03"))
    optimizer_counts <- table(optimizer$family, optimizer$parameter)
    optimizer_exactly_500_per_family_parameter <-
      length(optimizer_counts) == 42L && all(optimizer_counts == 500L)
    optimizer_numeric_and_feasibility_valid <-
      all(is.finite(optimizer$log2_ratio)) &&
      all(optimizer$feasible_at_solution) &&
      !any(optimizer$projection_applied)
  }
}
posterior_exists <- file.exists(posterior_path)
sampling_readiness_exists <- file.exists(sampling_readiness_path)
sampling_configuration_exists <- file.exists(sampling_configuration_path)
sampling_frozen_hashes_match <- FALSE
sampling_frozen_hash_evidence <- "missing sampler configuration or frozen source"
if (sampling_configuration_exists && all(file.exists(frozen_hash_sources))) {
  sampling_configuration <- utils::read.delim(
    sampling_configuration_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  configuration_lookup <- setNames(
    sampling_configuration$value,
    sampling_configuration$key
  )
  frozen_hash_contract <- c(
    generator = "generator_md5",
    aggregator = "aggregator_md5",
    combiner = "combiner_md5",
    family_contract = "contract_md5",
    prior = "conditioned_prior_md5"
  )
  if (all(unname(frozen_hash_contract) %in% names(configuration_lookup))) {
    observed_hashes <- unname(tools::md5sum(
      frozen_hash_sources[names(frozen_hash_contract)]
    ))
    expected_hashes <- unname(configuration_lookup[
      unname(frozen_hash_contract)
    ])
    sampling_frozen_hashes_match <-
      !anyNA(observed_hashes) && identical(observed_hashes, expected_hashes) &&
      identical(configuration_lookup[["family_targets_distinct"]], "TRUE") &&
      identical(
        configuration_lookup[["target_version"]],
        "family_conditioned_selected_map_voronoi_v1"
      )
    sampling_frozen_hash_evidence <- if (sampling_frozen_hashes_match) {
      paste(names(frozen_hash_contract), observed_hashes, sep = "=", collapse = "; ")
    } else {
      mismatched <- names(frozen_hash_contract)[
        is.na(observed_hashes) | observed_hashes != expected_hashes
      ]
      paste("mismatch:", paste(mismatched, collapse = ", "))
    }
  } else {
    sampling_frozen_hash_evidence <- paste0(
      "configuration missing hash keys: ",
      paste(
        setdiff(unname(frozen_hash_contract), names(configuration_lookup)),
        collapse = ", "
      )
    )
  }
}
posterior_schema_valid <- FALSE
posterior_selected_pair_coverage <- FALSE
posterior_parameter_coverage <- FALSE
posterior_minimum_draws_per_chain <- FALSE
posterior_numeric_values_valid <- FALSE
posterior_temperature_is_t1 <- FALSE
posterior_family_targets_distinct <- FALSE
sampling_readiness_valid <- FALSE
posterior_message <- paste0(
  "Missing released generalized-posterior draws. Optimizer endpoints are not posterior ",
  "samples and must not be substituted."
)

if (posterior_exists) {
  posterior <- utils::read.delim(
    posterior_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  posterior_schema_valid <- all(
    required_posterior_columns %in% names(posterior)
  )
  if (posterior_schema_valid) {
    posterior_selected_pair_coverage <- setequal(
      unique(posterior$warmup_label),
      selected$warmup_label
    ) && setequal(unique(posterior$family), c("C01", "C02", "C03"))
    expected_parameters <- c(
      "lam_max", "p_mis_base", "p_wgd", "p_misseg", "k_o_mis",
      "O2_crit", "n_O", "alpha_o2", "gamma_growth", "mu_hp",
      "gamma_mu", "buffer_smax", "buffer_beta", "buffer_n_exp"
    )
    posterior_parameter_coverage <- setequal(
      unique(posterior$parameter), expected_parameters
    )
    draw_counts <- table(
      interaction(
        posterior$warmup_label,
        posterior$chain,
        posterior$parameter,
        drop = TRUE
      )
    )
    posterior_minimum_draws_per_chain <-
      length(draw_counts) > 0L && all(draw_counts >= 500L)
    posterior_numeric_values_valid <-
      all(is.finite(posterior$vivo_natural)) &&
      all(is.finite(posterior$vitro_natural)) &&
      all(is.finite(posterior$log2_ratio)) &&
      all(posterior$vivo_natural > 0) &&
      all(posterior$vitro_natural > 0)
    posterior_temperature_is_t1 <- all(posterior$temperature == 1)
    posterior_family_targets_distinct <-
      all(posterior$target_version == "family_conditioned_selected_map_voronoi_v1")
    if (sampling_readiness_exists) {
      sampling_readiness <- utils::read.delim(
        sampling_readiness_path,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      sampling_readiness_valid <- nrow(sampling_readiness) > 0L &&
        all(sampling_readiness$passed)
    }
    posterior_message <- if (
      posterior_selected_pair_coverage &&
      posterior_parameter_coverage &&
      posterior_minimum_draws_per_chain &&
      posterior_numeric_values_valid &&
      posterior_temperature_is_t1 &&
      posterior_family_targets_distinct &&
      sampling_readiness_valid &&
      sampling_frozen_hashes_match
    ) {
      "Released T=1 draws from three distinct family-conditioned generalized posteriors pass structural and sampler-readiness checks."
    } else {
      "Posterior draw table exists but fails one or more minimum structural checks."
    }
  } else {
    posterior_message <- paste0(
      "Posterior draw table exists but lacks required columns: ",
      paste(
        setdiff(required_posterior_columns, names(posterior)),
        collapse = ", "
      )
    )
  }
}

readiness <- data.frame(
  check = c(
    "six_candidate_pairs_available",
    "two_candidates_per_family",
    "selection_metric_is_original_separate_invivo_objective",
    "exactly_one_pair_selected_per_family",
    "optimizer_endpoint_table_exists",
    "optimizer_endpoint_schema_valid",
    "optimizer_selected_pair_coverage",
    "optimizer_exactly_500_per_family_parameter",
    "optimizer_values_feasible_unprojected_and_finite",
    "posterior_draw_table_exists",
    "posterior_schema_valid",
    "posterior_selected_pair_coverage",
    "posterior_parameter_coverage",
    "posterior_minimum_500_draws_per_chain_parameter",
    "posterior_natural_values_positive_and_finite",
    "posterior_temperature_is_t1",
    "posterior_targets_are_distinct_family_conditioned_posteriors",
    "sampler_readiness_all_passed",
    "sampler_configuration_exists",
    "post_sampling_frozen_hashes_match",
    "safe_to_label_figure5f_as_three_distribution_comparison"
  ),
  passed = c(
    nrow(candidate) == 6L,
    all(table(candidate$family) == 2L),
    all(candidate$selection_metric ==
          "original separate-in-vivo objective"),
    nrow(selected) == 3L && all(table(selected$family) == 1L),
    optimizer_exists,
    optimizer_schema_valid,
    optimizer_selected_pair_coverage,
    optimizer_exactly_500_per_family_parameter,
    optimizer_numeric_and_feasibility_valid,
    posterior_exists,
    posterior_schema_valid,
    posterior_selected_pair_coverage,
    posterior_parameter_coverage,
    posterior_minimum_draws_per_chain,
    posterior_numeric_values_valid,
    posterior_temperature_is_t1,
    posterior_family_targets_distinct,
    sampling_readiness_valid,
    sampling_configuration_exists,
    sampling_frozen_hashes_match,
    optimizer_exists && optimizer_schema_valid &&
      optimizer_selected_pair_coverage &&
      optimizer_exactly_500_per_family_parameter &&
      optimizer_numeric_and_feasibility_valid &&
      posterior_exists && posterior_schema_valid &&
      posterior_selected_pair_coverage &&
      posterior_parameter_coverage &&
      posterior_minimum_draws_per_chain &&
      posterior_numeric_values_valid &&
      posterior_temperature_is_t1 &&
      posterior_family_targets_distinct &&
      sampling_readiness_valid &&
      sampling_configuration_exists &&
      sampling_frozen_hashes_match
  ),
  evidence = c(
    paste(nrow(candidate), "candidate joint pairs"),
    paste(names(table(candidate$family)), table(candidate$family), collapse = "; "),
    "joined to Figure4/invivo_fit_objective_ranking_500seeds.tsv",
    paste(selected$family, selected$warmup_label, sep = "=", collapse = "; "),
    if (optimizer_exists) normalizePath(optimizer_path) else optimizer_path,
    paste(required_optimizer_columns, collapse = ", "),
    paste(selected$warmup_label, collapse = "; "),
    "500 endpoints for each of 3 selected families x 14 parameters",
    "all selected endpoints feasible, unprojected, and finite on log2 ratio scale",
    if (posterior_exists) normalizePath(posterior_path) else posterior_path,
    paste(required_posterior_columns, collapse = ", "),
    paste(selected$warmup_label, collapse = "; "),
    "14 soft-coupled parameters",
    ">=500 retained draws for every selected pair x chain x parameter",
    "vivo_natural, vitro_natural, and log2_ratio must be finite; natural values >0",
    "released primary table must contain T=1 only",
    "target_version must be family_conditioned_selected_map_voronoi_v1 for every draw",
    if (sampling_readiness_exists) sampling_readiness_path else "missing sampler readiness",
    if (sampling_configuration_exists) sampling_configuration_path else "missing sampler configuration",
    sampling_frozen_hash_evidence,
    posterior_message
  ),
  stringsAsFactors = FALSE
)
write_tsv(readiness, readiness_output)

cat("Figure 5F selected pairs:\n")
print(selected[, c(
  "family", "warmup_label", "invivo_seed",
  "separate_invivo_objective", "separate_invivo_objective_rank"
)])
cat("\nPosterior readiness:\n")
print(readiness[, c("check", "passed")], row.names = FALSE)

if (!readiness$passed[readiness$check ==
                      "safe_to_label_figure5f_as_three_distribution_comparison"]) {
  message(
    "Figure 5F three-distribution rendering is intentionally blocked: ",
    posterior_message
  )
}
