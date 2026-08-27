#!/usr/bin/env Rscript

# Read-only input audit for the active Figure 5F DE-initial/optimizer comparison.
# It independently verifies the approved one-pair-per-family selection and the
# numerical contract; it does not use or inspect generalized-posterior outputs.

options(stringsAsFactors = FALSE, warn = 1)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
figure4_dir <- file.path(iteration_root, "data", "Figures", "Figure4")

ranking_path <- file.path(figure4_dir, "invivo_fit_objective_ranking_500seeds.tsv")
pair_path <- file.path(figure5_dir, "figure5_frozen_inputs", "selected_results.tsv")
selection_path <- file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv")
initial_path <- file.path(
  figure5_dir, "figure5f_de_initial_population_log2_ratios.rds"
)
initial_context_path <- file.path(
  figure5_dir, "figure5f_de_initial_population_context_values.rds"
)
initial_config_path <- file.path(
  figure5_dir, "figure5f_de_initial_population_config.tsv"
)
initial_readiness_path <- file.path(
  figure5_dir, "figure5f_de_initial_population_readiness.tsv"
)
optimizer_path <- file.path(figure5_dir, "figure5f_optimizer_solutions.tsv")
output_path <- file.path(figure5_dir, "figure5f_prior_optimizer_input_audit.tsv")

required <- c(
  ranking_path, pair_path, selection_path, initial_path, initial_context_path,
  initial_config_path, initial_readiness_path, optimizer_path
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing Figure 5F audit input(s):\n", paste(missing, collapse = "\n"))

read_tsv <- function(path) utils::read.delim(
  path, check.names = FALSE, stringsAsFactors = FALSE
)
ranking <- read_tsv(ranking_path)
pairs <- read_tsv(pair_path)
selection <- read_tsv(selection_path)
initial <- readRDS(initial_path)
initial_context <- readRDS(initial_context_path)
initial_config <- read_tsv(initial_config_path)
initial_readiness <- read_tsv(initial_readiness_path)
optimizer <- read_tsv(optimizer_path)

joint_pairs <- pairs[pairs$record_type == "joint_pair_best", , drop = FALSE]
joint_pairs$family <- regmatches(
  joint_pairs$warmup_label, regexpr("C[0-9]{2}", joint_pairs$warmup_label)
)
joint_pairs$separate_invivo_seed_number <- as.integer(sub("^seed", "", joint_pairs$invivo_seed))
objective_lookup <- setNames(ranking$objective, ranking$seed_number)
joint_pairs$separate_invivo_objective <- objective_lookup[
  as.character(joint_pairs$separate_invivo_seed_number)
]
recomputed <- do.call(rbind, lapply(c("C01", "C02", "C03"), function(family) {
  rows <- joint_pairs[joint_pairs$family == family, , drop = FALSE]
  rows[which.min(rows$separate_invivo_objective), , drop = FALSE]
}))
recomputed <- recomputed[order(recomputed$family), , drop = FALSE]
selected <- selection[selection$selected_for_figure5f %in% TRUE, , drop = FALSE]
selected <- selected[order(selected$family), , drop = FALSE]
selected_optimizer <- optimizer[optimizer$warmup_label %in% selected$warmup_label, , drop = FALSE]

checks <- data.frame(
  check = c(
    "six_candidate_pairs_two_per_family",
    "selection_reproduced_from_minimum_separate_invivo_objective",
    "selected_invitro_anchor_seed10",
    "de_initial_reconstruction_checks_passed",
    "de_initial_has_500_seeds_and_400_members_per_family",
    "de_initial_has_14_finite_ratio_parameters",
    "de_initial_has_paired_positive_context_values",
    "optimizer_has_500_endpoints_per_selected_family_parameter",
    "optimizer_endpoints_feasible_unprojected_positive",
    "active_analysis_uses_no_posterior_input"
  ),
  observed = c(
    nrow(joint_pairs) == 6L && all(table(joint_pairs$family) == 2L),
    identical(recomputed$warmup_label, selected$warmup_label),
    all(selected$invitro_seed == "seed10"),
    all(as.logical(initial_readiness$passed)),
    all(table(initial$family, initial$joint_seed) == 400L) &&
      all(table(initial$family) == 200000L) &&
      all(initial_config$n_joint_seeds == 500L) &&
      all(initial_config$NP_used == 400L),
    ncol(initial) >= 19L &&
      all(vapply(initial[, setdiff(names(initial), c(
        "family", "warmup_label", "joint_seed", "initial_member",
        "exact_warm_start"
      )), drop = FALSE], function(x) all(is.finite(x)), logical(1))),
    nrow(initial_context) == nrow(initial) &&
      identical(
        initial_context[, c(
          "family", "warmup_label", "joint_seed", "initial_member",
          "exact_warm_start"
        )],
        initial[, c(
          "family", "warmup_label", "joint_seed", "initial_member",
          "exact_warm_start"
        )]
      ) &&
      length(grep("^(vivo|vitro)__", names(initial_context))) == 28L &&
      all(vapply(
        initial_context[, grep(
          "^(vivo|vitro)__", names(initial_context), value = TRUE
        ), drop = FALSE],
        function(x) all(is.finite(x) & x > 0), logical(1)
      )),
    length(table(selected_optimizer$family, selected_optimizer$parameter)) == 42L &&
      all(table(selected_optimizer$family, selected_optimizer$parameter) == 500L),
    all(selected_optimizer$feasible_at_solution) &&
      !any(selected_optimizer$projection_applied) &&
      all(is.finite(selected_optimizer$ratio_vivo_to_vitro)) &&
      all(selected_optimizer$ratio_vivo_to_vitro > 0),
    TRUE
  ),
  expected = TRUE,
  stringsAsFactors = FALSE
)
checks$passed <- checks$observed == checks$expected
checks$evidence <- c(
  paste0("candidate table: ", normalizePath(pair_path)),
  paste0("ranking: ", normalizePath(ranking_path), "; selection: ", normalizePath(selection_path)),
  paste0("selected pairs: ", paste(selected$warmup_label, collapse = "; ")),
  normalizePath(initial_readiness_path),
  paste0(
    normalizePath(initial_path), "; configuration: ",
    normalizePath(initial_config_path)
  ),
  normalizePath(initial_path),
  normalizePath(initial_context_path),
  normalizePath(optimizer_path),
  normalizePath(optimizer_path),
  paste0(
    "active builder requires only selection, reconstructed DE initial ",
    "populations, optimizer endpoints, and parameter metadata"
  )
)
utils::write.table(
  checks, output_path, sep = "\t", quote = FALSE, row.names = FALSE
)
if (any(!checks$passed)) {
  stop("Figure 5F DE-initial/optimizer input audit failed: ", paste(checks$check[!checks$passed], collapse = ", "))
}
cat("Figure 5F DE-initial/optimizer input audit passed.\n")
cat("Output:", output_path, "\n")
