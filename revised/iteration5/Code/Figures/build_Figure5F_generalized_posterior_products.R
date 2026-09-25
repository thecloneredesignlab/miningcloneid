#!/usr/bin/env Rscript

# Build plot-ready Figure 5F products from validated generalized-posterior
# draws. This script never substitutes optimizer endpoints for posterior draws.
# It requires the sampler-level convergence gate to pass before releasing any
# plot products.

options(stringsAsFactors = FALSE, warn = 1)

resolve_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE))
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
}

script_dir <- resolve_script_dir()
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
legacy_sampling_dir <- file.path(figure5_dir, "generalized_posterior")
sampling_dir <- file.path(figure5_dir, "generalized_posterior_family_conditioned")

paths <- list(
  generator = file.path(script_dir, "generate_Figure5F_generalized_posterior.R"),
  aggregator = file.path(script_dir, "aggregate_Figure5F_family_posterior.R"),
  combiner = file.path(script_dir, "combine_Figure5F_family_posteriors.R"),
  proposal_calibration = file.path(script_dir, "calibrate_Figure5F_pilot_proposal.R"),
  selection = file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv"),
  map_audit = file.path(figure5_dir, "figure5f_map_replay_audit.tsv"),
  prior = file.path(
    sampling_dir, "figure5f_family_conditioned_prior_draws.tsv"
  ),
  family_contract = file.path(
    sampling_dir, "figure5f_family_basin_contract.rds"
  ),
  prior_config = file.path(figure5_dir, "figure5f_prior_sampling_config.tsv"),
  optimizer = file.path(figure5_dir, "figure5f_optimizer_solutions.tsv"),
  proposal_covariance = file.path(
    legacy_sampling_dir, "figure5f_pilot2_proposal_covariance.rds"
  ),
  primary = file.path(figure5_dir, "figure5f_posterior_draws.tsv"),
  all_temperatures = file.path(
    sampling_dir, "figure5f_generalized_posterior_draws.tsv"
  ),
  diagnostics = file.path(
    sampling_dir, "figure5f_generalized_posterior_diagnostics.tsv"
  ),
  readiness = file.path(
    sampling_dir, "figure5f_generalized_posterior_readiness.tsv"
  ),
  configuration = file.path(
    sampling_dir, "figure5f_generalized_posterior_config.tsv"
  ),
  parameter_meta = file.path(figure5_dir, "parameter_function_groups.tsv")
)
missing <- unlist(paths, use.names = FALSE)
missing <- missing[!file.exists(missing)]
if (length(missing)) {
  stop("Missing validated Figure 5F product input(s):\n", paste(missing, collapse = "\n"))
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  allowed <- normalizePath(figure5_dir, mustWork = TRUE)
  if (!startsWith(normalized, paste0(allowed, .Platform$file.sep))) {
    stop("Refusing to write outside iteration2 Figure5 data: ", normalized)
  }
  utils::write.table(
    x, normalized, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
  invisible(normalized)
}

prior <- read_tsv(paths$prior)
prior_config <- read_tsv(paths$prior_config)
selection <- read_tsv(paths$selection)
optimizer <- read_tsv(paths$optimizer)
primary <- read_tsv(paths$primary)
all_temperatures <- read_tsv(paths$all_temperatures)
diagnostics <- read_tsv(paths$diagnostics)
readiness <- read_tsv(paths$readiness)
configuration <- read_tsv(paths$configuration)
parameter_meta <- read_tsv(paths$parameter_meta)

if (!nrow(readiness) || any(!readiness$passed)) {
  failed <- readiness$check[!readiness$passed]
  stop(
    "Figure 5F generalized-posterior readiness gate has not passed: ",
    paste(failed, collapse = ", ")
  )
}
configuration_lookup <- setNames(configuration$value, configuration$key)
if (!identical(configuration_lookup[["safe_for_primary_figure"]], "TRUE") ||
    !identical(configuration_lookup[["family_targets_distinct"]], "TRUE") ||
    !identical(
      configuration_lookup[["target_version"]],
      "family_conditioned_selected_map_voronoi_v1"
    ) ||
    !grepl("generalized posterior", configuration_lookup[["target_name"]], fixed = TRUE)) {
  stop("Figure 5F configuration does not certify the declared generalized posterior.")
}
frozen_hash_contract <- c(
  generator = "generator_md5",
  aggregator = "aggregator_md5",
  combiner = "combiner_md5",
  family_contract = "contract_md5",
  prior = "conditioned_prior_md5"
)
missing_hash_keys <- setdiff(
  unname(frozen_hash_contract),
  names(configuration_lookup)
)
if (length(missing_hash_keys)) {
  stop(
    "Figure 5F configuration lacks frozen-hash key(s): ",
    paste(missing_hash_keys, collapse = ", ")
  )
}
observed_frozen_hashes <- unname(tools::md5sum(
  unlist(paths[names(frozen_hash_contract)], use.names = FALSE)
))
expected_frozen_hashes <- unname(configuration_lookup[
  unname(frozen_hash_contract)
])
if (anyNA(observed_frozen_hashes) ||
    !identical(observed_frozen_hashes, expected_frozen_hashes)) {
  mismatched <- names(frozen_hash_contract)[
    is.na(observed_frozen_hashes) |
      observed_frozen_hashes != expected_frozen_hashes
  ]
  stop(
    "Figure 5F frozen input/code hash mismatch after sampling: ",
    paste(mismatched, collapse = ", ")
  )
}

parameters <- c(
  "lam_max", "p_mis_base", "p_wgd",
  "p_misseg", "k_o_mis",
  "O2_crit", "n_O",
  "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
families <- c("C01", "C02", "C03")
display_temperatures <- c(0.5, 1, 2)

required_prior <- c("family", "parameter", "parameter_order", "draw_id", "log2_ratio")
required_primary <- c(
  "family", "warmup_label", "selected_seed", "chain", "replicate", "draw",
  "temperature", "parameter", "vivo_transformed", "vitro_transformed",
  "vivo_natural", "vitro_natural", "log2_ratio", "target_version"
)
required_optimizer <- c(
  "family", "warmup_label", "seed_number", "parameter", "log2_ratio",
  "either_context_at_active_bound"
)
required_diag <- c(
  "variable", "scope", "temperature", "family", "rhat", "ess_bulk",
  "ess_tail", "mcse_mean"
)
if (!all(required_prior %in% names(prior)) ||
    !all(required_optimizer %in% names(optimizer)) ||
    !all(required_primary %in% names(primary)) ||
    !all(required_primary %in% names(all_temperatures)) ||
    !all(required_diag %in% names(diagnostics))) {
  stop("Figure 5F generalized-posterior input schema is incomplete.")
}
selected <- selection[selection$selected_for_figure5f %in% TRUE, , drop = FALSE]
selected <- selected[order(match(selected$family, families)), , drop = FALSE]
optimizer <- optimizer[
  optimizer$warmup_label %in% selected$warmup_label,
  , drop = FALSE
]
if (!setequal(unique(prior$parameter), parameters) ||
    !setequal(unique(prior$family), families) ||
    !setequal(unique(primary$parameter), parameters) ||
    !setequal(unique(primary$family), families) ||
    !all(primary$temperature == 1) ||
    !all(primary$target_version == "family_conditioned_selected_map_voronoi_v1") ||
    any(!is.finite(primary$log2_ratio))) {
  stop("Figure 5F primary posterior draw coverage is invalid.")
}
optimizer_counts <- table(optimizer$family, optimizer$parameter)
if (nrow(selected) != 3L || !identical(selected$family, families) ||
    !setequal(unique(optimizer$family), families) ||
    !setequal(unique(optimizer$parameter), parameters) ||
    any(optimizer_counts != 500L) ||
    any(!is.finite(optimizer$log2_ratio))) {
  stop("Figure 5F selected optimizer-endpoint coverage is invalid.")
}
primary_counts <- table(primary$family, primary$chain, primary$parameter)
primary_counts <- primary_counts[primary_counts > 0]
if (!length(primary_counts) || any(primary_counts < 500L)) {
  stop("Figure 5F requires at least 500 retained T=1 draws per chain and parameter.")
}
if (!all(display_temperatures %in% unique(all_temperatures$temperature))) {
  stop("Figure 5F temperature-sensitivity draws are incomplete.")
}

meta <- parameter_meta[
  match(parameters, parameter_meta$parameter),
  c("parameter", "parameter_group", "parameter_order", "transformation"),
  drop = FALSE
]
if (anyNA(meta$parameter) || !identical(meta$parameter, parameters)) {
  stop("Figure 5F parameter metadata does not match the 14 coupled parameters.")
}
prior_config <- prior_config[
  match(parameters, prior_config$parameter),
  , drop = FALSE
]
if (anyNA(prior_config$parameter) || !identical(prior_config$parameter, parameters)) {
  stop("Figure 5F prior configuration does not match the 14 coupled parameters.")
}

quantile_row <- function(values, prefix = "") {
  q <- stats::quantile(
    values, c(0.001, 0.05, 0.25, 0.50, 0.75, 0.95, 0.999),
    names = FALSE, type = 8
  )
  out <- data.frame(
    q001 = q[[1L]], q05 = q[[2L]], q25 = q[[3L]], median = q[[4L]],
    q75 = q[[5L]], q95 = q[[6L]], q999 = q[[7L]],
    width90 = q[[6L]] - q[[2L]],
    stringsAsFactors = FALSE
  )
  if (nzchar(prefix)) names(out) <- paste0(prefix, names(out))
  out
}

prior_groups <- split(prior, interaction(prior$parameter, prior$family, drop = TRUE))
prior_summary <- do.call(rbind, lapply(prior_groups, function(group) {
  cbind(
    data.frame(
      parameter = group$parameter[[1L]],
      family = group$family[[1L]],
      stringsAsFactors = FALSE
    ),
    quantile_row(group$log2_ratio, "prior_")
  )
}))
rownames(prior_summary) <- NULL

primary_groups <- split(
  primary,
  interaction(primary$parameter, primary$family, drop = TRUE)
)
posterior_summary <- do.call(rbind, lapply(primary_groups, function(group) {
  q <- quantile_row(group$log2_ratio)
  config_row <- prior_config[
    prior_config$parameter == group$parameter[[1L]],
    , drop = FALSE
  ]
  lower <- config_row$joint_lower_transformed[[1L]]
  upper <- config_row$joint_upper_transformed[[1L]]
  width <- upper - lower
  vivo_unit <- (group$vivo_transformed - lower) / width
  vitro_unit <- (group$vitro_transformed - lower) / width
  boundary_fraction_1pct <- mean(
    vivo_unit <= 0.01 | vivo_unit >= 0.99 |
      vitro_unit <= 0.01 | vitro_unit >= 0.99
  )
  data.frame(
    parameter = group$parameter[[1L]],
    family = group$family[[1L]],
    warmup_label = group$warmup_label[[1L]],
    selected_seed = group$selected_seed[[1L]],
    n_chains = length(unique(group$chain)),
    n_draws = nrow(group),
    q001 = q$q001,
    q05 = q$q05,
    q25 = q$q25,
    median = q$median,
    q75 = q$q75,
    q95 = q$q95,
    q999 = q$q999,
    width90 = q$width90,
    boundary_fraction_1pct = boundary_fraction_1pct,
    probability_higher_invivo = mean(group$log2_ratio > 0),
    direction = if (q$q05 > 0) {
      "higher in vivo"
    } else if (q$q95 < 0) {
      "lower in vivo"
    } else {
      "overlaps equality"
    },
    stringsAsFactors = FALSE
  )
}))
posterior_summary <- merge(
  posterior_summary, prior_summary,
  by = c("parameter", "family"), all.x = TRUE, sort = FALSE
)
posterior_summary <- merge(
  posterior_summary, meta,
  by = "parameter", all.x = TRUE, sort = FALSE
)
posterior_summary$contraction_ratio90 <-
  posterior_summary$width90 / posterior_summary$prior_width90

optimizer_groups <- split(
  optimizer,
  interaction(optimizer$parameter, optimizer$family, drop = TRUE)
)
optimizer_summary <- do.call(rbind, lapply(optimizer_groups, function(group) {
  q <- quantile_row(group$log2_ratio, "optimizer_")
  data.frame(
    parameter = group$parameter[[1L]],
    family = group$family[[1L]],
    optimizer_n_endpoints = nrow(group),
    optimizer_q001 = q$optimizer_q001,
    optimizer_q05 = q$optimizer_q05,
    optimizer_q25 = q$optimizer_q25,
    optimizer_median = q$optimizer_median,
    optimizer_q75 = q$optimizer_q75,
    optimizer_q95 = q$optimizer_q95,
    optimizer_q999 = q$optimizer_q999,
    optimizer_width90 = q$optimizer_width90,
    optimizer_probability_higher_invivo = mean(group$log2_ratio > 0),
    optimizer_active_bound_fraction = mean(group$either_context_at_active_bound),
    optimizer_direction = if (q$optimizer_q05 > 0) {
      "higher in vivo"
    } else if (q$optimizer_q95 < 0) {
      "lower in vivo"
    } else {
      "overlaps equality"
    },
    stringsAsFactors = FALSE
  )
}))
posterior_summary <- merge(
  posterior_summary,
  optimizer_summary,
  by = c("parameter", "family"),
  all.x = TRUE,
  sort = FALSE
)
posterior_summary$posterior_to_optimizer_width90 <-
  posterior_summary$width90 / posterior_summary$optimizer_width90
posterior_summary <- posterior_summary[order(
  posterior_summary$parameter_order,
  match(posterior_summary$family, families)
), , drop = FALSE]
rownames(posterior_summary) <- NULL

ratio_diag <- diagnostics[
  diagnostics$scope == "paired_log2_ratios" & diagnostics$temperature == 1,
  , drop = FALSE
]
ratio_diag$parameter <- ratio_diag$variable
ratio_diag <- ratio_diag[order(
  match(ratio_diag$parameter, parameters), match(ratio_diag$family, families)
), , drop = FALSE]
if (nrow(ratio_diag) != length(parameters) * length(families) ||
    anyDuplicated(ratio_diag[, c("parameter", "family")]) ||
    !setequal(unique(ratio_diag$parameter), parameters) ||
    !setequal(unique(ratio_diag$family), families) ||
    any(!is.finite(ratio_diag$rhat))) {
  stop("Figure 5F family-specific ratio diagnostics are incomplete.")
}

cross_summary <- do.call(rbind, lapply(parameters, function(parameter) {
  group <- posterior_summary[posterior_summary$parameter == parameter, , drop = FALSE]
  directions <- unique(group$direction)
  direction_agreement <- length(directions) == 1L &&
    !identical(directions[[1L]], "overlaps equality")
  optimizer_directions <- unique(group$optimizer_direction)
  optimizer_direction_agreement <- length(optimizer_directions) == 1L &&
    !identical(optimizer_directions[[1L]], "overlaps equality")
  familywise_posterior_optimizer_direction_matches <- sum(
    group$direction == group$optimizer_direction
  )
  diag_row <- ratio_diag[ratio_diag$parameter == parameter, , drop = FALSE]
  concentrated <- max(group$contraction_ratio90) <= 0.50
  not_boundary_limited <- max(group$boundary_fraction_1pct) <= 0.10
  all_family_convergence_passed <-
    all(diag_row$rhat <= 1.01) &&
    all(diag_row$ess_bulk >= 400) &&
    all(diag_row$ess_tail >= 400)
  data.frame(
    parameter = parameter,
    parameter_group = group$parameter_group[[1L]],
    parameter_order = group$parameter_order[[1L]],
    max_contraction_ratio90 = max(group$contraction_ratio90),
    min_contraction_ratio90 = min(group$contraction_ratio90),
    max_boundary_fraction_1pct = max(group$boundary_fraction_1pct),
    family_median_min = min(group$median),
    family_median_max = max(group$median),
    family_median_span = max(group$median) - min(group$median),
    max_family_rhat = max(diag_row$rhat),
    min_family_ess_bulk = min(diag_row$ess_bulk),
    min_family_ess_tail = min(diag_row$ess_tail),
    direction_agreement = direction_agreement,
    direction = if (direction_agreement) directions[[1L]] else "not directionally resolved",
    optimizer_direction_agreement = optimizer_direction_agreement,
    optimizer_direction = if (optimizer_direction_agreement) {
      optimizer_directions[[1L]]
    } else {
      "not directionally resolved"
    },
    posterior_optimizer_direction_concordance =
      direction_agreement && optimizer_direction_agreement &&
      identical(directions[[1L]], optimizer_directions[[1L]]),
    familywise_posterior_optimizer_direction_matches_of_3 =
      familywise_posterior_optimizer_direction_matches,
    posterior_concentrated_vs_prior = concentrated,
    not_boundary_limited = not_boundary_limited,
    all_family_convergence_passed = all_family_convergence_passed,
    practical_identifiability_summary = if (
      concentrated && not_boundary_limited && all_family_convergence_passed
    ) {
      "concentrated in all family-conditioned targets"
    } else {
      "weakly resolved"
    },
    stringsAsFactors = FALSE
  )
}))
rownames(cross_summary) <- NULL

support_for_parameter <- function(parameter) {
  row <- prior_config[prior_config$parameter == parameter, , drop = FALSE]
  if (identical(row$transformation[[1L]], "log10")) {
    width <- (row$joint_upper_transformed[[1L]] -
      row$joint_lower_transformed[[1L]]) * log2(10)
    c(-width, width)
  } else if (row$joint_lower_natural[[1L]] > 0) {
    width <- log2(
      row$joint_upper_natural[[1L]] / row$joint_lower_natural[[1L]]
    )
    c(-width, width)
  } else {
    c(-Inf, Inf)
  }
}

bounded_density <- function(values, grid, support) {
  values <- values[is.finite(values)]
  if (length(values) < 2L || length(unique(values)) < 2L) {
    y <- numeric(length(grid))
    y[[which.min(abs(grid - values[[1L]]))]] <- 1
    return(data.frame(x = grid, density = y, bandwidth = NA_real_))
  }
  bw <- stats::bw.nrd0(values)
  if (!is.finite(bw) || bw <= 0) {
    bw <- max(diff(range(grid)) / 100, .Machine$double.eps^0.25)
  }
  augmented <- values
  multiplier <- 1L
  if (is.finite(support[[1L]])) {
    augmented <- c(augmented, 2 * support[[1L]] - values)
    multiplier <- multiplier + 1L
  }
  if (is.finite(support[[2L]])) {
    augmented <- c(augmented, 2 * support[[2L]] - values)
    multiplier <- multiplier + 1L
  }
  dens <- stats::density(
    augmented, bw = bw, from = min(grid), to = max(grid),
    n = length(grid), cut = 0
  )
  data.frame(x = dens$x, density = dens$y * multiplier, bandwidth = bw)
}

density_rows <- list()
density_index <- 1L
for (parameter in parameters) {
  prior_values <- prior$log2_ratio[prior$parameter == parameter]
  posterior_values <- primary$log2_ratio[primary$parameter == parameter]
  optimizer_values <- optimizer$log2_ratio[optimizer$parameter == parameter]
  support <- support_for_parameter(parameter)
  prior_q <- stats::quantile(prior_values, c(0.001, 0.999), names = FALSE, type = 8)
  post_q <- stats::quantile(posterior_values, c(0.001, 0.999), names = FALSE, type = 8)
  optimizer_q <- stats::quantile(optimizer_values, c(0.001, 0.999), names = FALSE, type = 8)
  display_from <- min(prior_q[[1L]], post_q[[1L]], optimizer_q[[1L]], 0)
  display_to <- max(prior_q[[2L]], post_q[[2L]], optimizer_q[[2L]], 0)
  padding <- max(0.025 * (display_to - display_from), 0.02)
  display_from <- max(display_from - padding, support[[1L]])
  display_to <- min(display_to + padding, support[[2L]])
  grid <- seq(display_from, display_to, length.out = 401L)
  parameter_rows <- list()
  parameter_index <- 1L
  for (family in families) {
    prior_family <- prior$log2_ratio[
      prior$parameter == parameter & prior$family == family
    ]
    posterior_family <- primary$log2_ratio[
      primary$parameter == parameter & primary$family == family
    ]
    optimizer_family <- optimizer$log2_ratio[
      optimizer$parameter == parameter & optimizer$family == family
    ]
    for (role in c("prior", "generalized_posterior", "optimizer_endpoints")) {
      values <- switch(
        role,
        prior = prior_family,
        generalized_posterior = posterior_family,
        optimizer_endpoints = optimizer_family
      )
      dens <- bounded_density(values, grid, support)
      dens$density_scaled_distribution <- dens$density / max(dens$density)
      parameter_rows[[parameter_index]] <- cbind(
        data.frame(
          parameter = parameter,
          family = family,
          distribution_role = role,
          distribution = switch(
            role,
            prior = "Prior",
            generalized_posterior = "Generalized posterior (T=1)",
            optimizer_endpoints = "Optimizer endpoints"
          ),
          display_from = display_from,
          display_to = display_to,
          support_lower = support[[1L]],
          support_upper = support[[2L]],
          displayed_tail_fraction = mean(values < display_from | values > display_to),
          stringsAsFactors = FALSE
        ),
        dens
      )
      parameter_index <- parameter_index + 1L
    }
  }
  parameter_table <- do.call(rbind, parameter_rows)
  peak <- max(parameter_table$density)
  parameter_table$density_scaled_parameter <- parameter_table$density / peak
  density_rows[[density_index]] <- parameter_table
  density_index <- density_index + 1L
}
density_table <- do.call(rbind, density_rows)
density_table <- merge(
  density_table, meta,
  by = "parameter", all.x = TRUE, sort = FALSE
)
density_table <- density_table[order(
  density_table$parameter_order,
  match(density_table$family, families),
  match(
    density_table$distribution_role,
    c("prior", "generalized_posterior", "optimizer_endpoints")
  ),
  density_table$x
), , drop = FALSE]
rownames(density_table) <- NULL

temperature_groups <- split(
  all_temperatures[all_temperatures$temperature %in% display_temperatures, , drop = FALSE],
  interaction(
    all_temperatures$temperature,
    all_temperatures$parameter,
    all_temperatures$family,
    drop = TRUE
  )
)
temperature_summary <- do.call(rbind, lapply(temperature_groups, function(group) {
  q <- quantile_row(group$log2_ratio)
  data.frame(
    temperature = group$temperature[[1L]],
    parameter = group$parameter[[1L]],
    family = group$family[[1L]],
    n_draws = nrow(group),
    q05 = q$q05,
    median = q$median,
    q95 = q$q95,
    width90 = q$width90,
    probability_higher_invivo = mean(group$log2_ratio > 0),
    direction = if (q$q05 > 0) {
      "higher in vivo"
    } else if (q$q95 < 0) {
      "lower in vivo"
    } else {
      "overlaps equality"
    },
    stringsAsFactors = FALSE
  )
}))
temperature_summary <- merge(
  temperature_summary, meta,
  by = "parameter", all.x = TRUE, sort = FALSE
)

# Learning-temperature direction summaries are secondary sensitivity results,
# not part of the primary T=1 release gate.  They are only interpretable when
# the corresponding two chains within that family-conditioned target have adequate
# diagnostics at that temperature.  Keep failed cells in the output and mark
# them ineligible instead of silently classifying them as directionally stable.
temperature_diagnostics <- diagnostics[
  diagnostics$scope == "paired_log2_ratios" &
    diagnostics$temperature %in% display_temperatures,
  c("variable", "temperature", "family", "rhat", "ess_bulk", "ess_tail"),
  drop = FALSE
]
names(temperature_diagnostics)[names(temperature_diagnostics) == "variable"] <-
  "parameter"
if (nrow(temperature_diagnostics) !=
      length(parameters) * length(families) * length(display_temperatures) ||
    anyDuplicated(temperature_diagnostics[, c("parameter", "temperature", "family")]) ||
    !setequal(unique(temperature_diagnostics$parameter), parameters) ||
    !setequal(unique(temperature_diagnostics$family), families)) {
  stop("Figure 5F learning-temperature diagnostics are incomplete or duplicated.")
}
names(temperature_diagnostics)[names(temperature_diagnostics) == "rhat"] <-
  "temperature_rhat_within_family"
names(temperature_diagnostics)[names(temperature_diagnostics) == "ess_bulk"] <-
  "temperature_ess_bulk_within_family"
names(temperature_diagnostics)[names(temperature_diagnostics) == "ess_tail"] <-
  "temperature_ess_tail_within_family"
temperature_summary <- merge(
  temperature_summary,
  temperature_diagnostics,
  by = c("parameter", "temperature", "family"),
  all.x = TRUE,
  sort = FALSE
)
temperature_summary$temperature_diagnostic_evaluable <- with(
  temperature_summary,
  is.finite(temperature_rhat_within_family) &
    is.finite(temperature_ess_bulk_within_family) &
    is.finite(temperature_ess_tail_within_family) &
    temperature_rhat_within_family <= 1.05 &
    temperature_ess_bulk_within_family >= 100 &
    temperature_ess_tail_within_family >= 100
)
temperature_summary <- temperature_summary[order(
  temperature_summary$parameter_order,
  match(temperature_summary$family, families),
  temperature_summary$temperature
), , drop = FALSE]
rownames(temperature_summary) <- NULL

temperature_cells <- split(
  temperature_summary,
  interaction(temperature_summary$parameter, temperature_summary$family, drop = TRUE)
)
temperature_stability_cells <- do.call(rbind, lapply(
  temperature_cells,
  function(group) {
    complete_temperatures <- setequal(
      group$temperature,
      display_temperatures
    )
    evaluable <- complete_temperatures &&
      all(group$temperature_diagnostic_evaluable)
    data.frame(
      parameter = group$parameter[[1L]],
      family = group$family[[1L]],
      temperature_sensitivity_evaluable = evaluable,
      temperature_direction_stable =
        evaluable && length(unique(group$direction)) == 1L,
      stringsAsFactors = FALSE
    )
  }
))
rownames(temperature_stability_cells) <- NULL
temperature_stability_summary <- aggregate(
  cbind(
    temperature_sensitivity_evaluable,
    temperature_direction_stable
  ) ~ parameter,
  data = temperature_stability_cells,
  FUN = sum
)
names(temperature_stability_summary)[
  names(temperature_stability_summary) == "temperature_sensitivity_evaluable"
] <- "temperature_evaluable_family_count_of_3"
names(temperature_stability_summary)[
  names(temperature_stability_summary) == "temperature_direction_stable"
] <- "temperature_stable_family_count_of_3"
cross_summary <- merge(
  cross_summary,
  temperature_stability_summary,
  by = "parameter",
  all.x = TRUE,
  sort = FALSE
)
cross_summary$all_family_directions_stable_T0p5_T1_T2 <-
  cross_summary$temperature_evaluable_family_count_of_3 == 3L &
  cross_summary$temperature_stable_family_count_of_3 == 3L
cross_summary <- cross_summary[order(cross_summary$parameter_order), , drop = FALSE]
rownames(cross_summary) <- NULL

density_cell_counts <- table(
  density_table$parameter,
  density_table$family,
  density_table$distribution_role
)
density_cell_peaks <- tapply(
  density_table$density_scaled_distribution,
  interaction(
    density_table$parameter,
    density_table$family,
    density_table$distribution_role,
    drop = TRUE
  ),
  max
)

plot_readiness <- data.frame(
  check = c(
    "sampler_readiness_all_passed",
    "primary_temperature_is_t1",
    "three_distinct_family_conditioned_targets",
    "fourteen_paired_parameters",
    "minimum_500_draws_per_chain_parameter",
    "exactly_500_optimizer_endpoints_per_family_parameter",
    "density_rows_14x3x3x401",
    "density_roles_prior_posterior_optimizer_complete",
    "density_grid_401_per_parameter_family_role",
    "each_density_curve_peak_normalized",
    "displayed_tail_fraction_max_le_0p002",
    "family_ratio_rhat_max_le_1p01",
    "temperature_sensitivity_cells_complete",
    "temperature_sensitivity_diagnostic_rows_complete",
    "temperature_direction_stability_cells"
  ),
  observed = c(
    all(readiness$passed),
    paste(sort(unique(primary$temperature)), collapse = ","),
    paste(sort(unique(primary$family)), collapse = ","),
    length(unique(primary$parameter)),
    min(primary_counts),
    paste(sort(unique(as.vector(optimizer_counts))), collapse = ","),
    nrow(density_table),
    paste(sort(unique(density_table$distribution_role)), collapse = ","),
    paste(sort(unique(as.vector(density_cell_counts))), collapse = ","),
    max(abs(density_cell_peaks - 1)),
    max(density_table$displayed_tail_fraction),
    max(ratio_diag$rhat),
    length(temperature_cells),
    nrow(temperature_diagnostics),
    paste(
      sum(temperature_stability_cells$temperature_direction_stable),
      "stable of",
      sum(temperature_stability_cells$temperature_sensitivity_evaluable),
      "evaluable /",
      nrow(temperature_stability_cells),
      "total"
    )
  ),
  passed = c(
    all(readiness$passed),
    all(primary$temperature == 1),
    setequal(unique(primary$family), families),
    setequal(unique(primary$parameter), parameters),
    min(primary_counts) >= 500L,
    all(optimizer_counts == 500L),
    nrow(density_table) == 14L * 3L * 3L * 401L,
    setequal(
      unique(density_table$distribution_role),
      c("prior", "generalized_posterior", "optimizer_endpoints")
    ),
    length(density_cell_counts) == 14L * 3L * 3L &&
      all(density_cell_counts == 401L),
    all(is.finite(density_table$density_scaled_distribution)) &&
      all(density_table$density_scaled_distribution >= 0) &&
      all(density_table$density_scaled_distribution <= 1 + 1e-12) &&
      max(abs(density_cell_peaks - 1)) <= 1e-12,
    max(density_table$displayed_tail_fraction) <= 0.002,
    max(ratio_diag$rhat) <= 1.01,
    length(temperature_cells) == 14L * 3L,
    nrow(temperature_diagnostics) == 14L * 3L * 3L,
    nrow(temperature_stability_cells) == 14L * 3L &&
      all(!is.na(temperature_stability_cells$temperature_sensitivity_evaluable)) &&
      all(!is.na(temperature_stability_cells$temperature_direction_stable))
  ),
  stringsAsFactors = FALSE
)
if (any(!plot_readiness$passed)) {
  stop(
    "Figure 5F plot-product validation failed: ",
    paste(plot_readiness$check[!plot_readiness$passed], collapse = ", ")
  )
}

write_tsv(
  density_table,
  file.path(figure5_dir, "figure5f_three_distribution_density.tsv")
)
write_tsv(
  posterior_summary,
  file.path(figure5_dir, "figure5f_three_distribution_summary.tsv")
)
write_tsv(
  cross_summary,
  file.path(figure5_dir, "figure5f_generalized_posterior_cross_family.tsv")
)
write_tsv(
  temperature_summary,
  file.path(figure5_dir, "figure5f_generalized_posterior_temperature_sensitivity.tsv")
)
write_tsv(
  plot_readiness,
  file.path(figure5_dir, "figure5f_generalized_posterior_plot_readiness.tsv")
)

provenance_paths <- unlist(paths, use.names = FALSE)
provenance <- data.frame(
  role = names(paths),
  path = normalizePath(provenance_paths, mustWork = TRUE),
  md5 = unname(tools::md5sum(provenance_paths)),
  stringsAsFactors = FALSE
)
write_tsv(
  provenance,
  file.path(figure5_dir, "figure5f_generalized_posterior_product_provenance.tsv")
)

cat("Figure 5F generalized-posterior plot products passed all checks.\n")
print(cross_summary[, c(
  "parameter", "max_contraction_ratio90", "max_family_rhat",
  "direction", "practical_identifiability_summary"
)], row.names = FALSE)
