#!/usr/bin/env Rscript

# Freeze the C01/C02/C03 posterior-conditioning contract and generate the
# corresponding family-conditioned induced priors.  Family membership is the
# nearest selected family MAP in the exact bound-normalized 18-dimensional
# in-vivo mechanistic parameter space.  The three Voronoi cells form a fixed,
# exhaustive partition (ties use C01 < C02 < C03).  No optimizer endpoint is
# used as a posterior or prior draw.

options(stringsAsFactors = FALSE, warn = 1)

resolve_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Unable to resolve script path.")
  dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE))
}

write_tsv_atomic <- function(x, path) {
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  utils::write.table(
    x, temporary, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
  if (!file.rename(temporary, path)) stop("Failed to atomically write: ", path)
}

save_rds_atomic <- function(x, path) {
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(x, temporary, version = 3)
  if (!file.rename(temporary, path)) stop("Failed to atomically write: ", path)
}

script_dir <- resolve_script_dir()
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
legacy_sampling_dir <- file.path(figure5_dir, "generalized_posterior")
output_dir <- file.path(figure5_dir, "generalized_posterior_family_conditioned")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

selection_path <- file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv")
coordinate_path <- file.path(
  legacy_sampling_dir, "figure5f_sampling_coordinate_contract.tsv"
)
prior_config_path <- file.path(figure5_dir, "figure5f_prior_sampling_config.tsv")
contract_rds_path <- file.path(output_dir, "figure5f_family_basin_contract.rds")
contract_tsv_path <- file.path(output_dir, "figure5f_family_basin_centers.tsv")
prior_draw_path <- file.path(output_dir, "figure5f_family_conditioned_prior_draws.tsv")
prior_config_output <- file.path(
  output_dir, "figure5f_family_conditioned_prior_config.tsv"
)

result_root <- Sys.getenv(
  "FIGURE5F_RESULT_ROOT",
  unset = paste0(
    "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/",
    "fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540"
  )
)
result_root <- normalizePath(result_root, mustWork = TRUE)

required <- c(selection_path, coordinate_path, prior_config_path)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing family-conditioning input(s):\n", paste(missing, collapse = "\n"))
}

families <- c("C01", "C02", "C03")
selection <- utils::read.delim(selection_path, check.names = FALSE)
selection <- selection[selection$selected_for_figure5f %in% TRUE, , drop = FALSE]
selection <- selection[order(match(selection$family, families)), , drop = FALSE]
if (nrow(selection) != 3L || !identical(selection$family, families)) {
  stop("Expected one selected MAP for each of C01, C02, and C03.")
}

coordinate <- utils::read.delim(coordinate_path, check.names = FALSE)
required_coordinate <- c(
  "coordinate", "role", "parameter", "original_name", "lower", "upper", "width"
)
if (nrow(coordinate) != 40L ||
    !all(required_coordinate %in% names(coordinate)) ||
    any(!is.finite(coordinate$width)) || any(coordinate$width <= 0)) {
  stop("The 40-dimensional sampling-coordinate contract is invalid.")
}

ordinary_in_vivo_features <- c(
  "log10_o2_S0", "log10_kappa_O", "log10_eta_o2", "log10_k_clear"
)
paired_in_vivo_features <- coordinate$coordinate[coordinate$role == "vivo"]
family_features <- c(ordinary_in_vivo_features, paired_in_vivo_features)
if (length(family_features) != 18L ||
    any(!family_features %in% coordinate$coordinate) ||
    anyDuplicated(family_features)) {
  stop("Family conditioning requires exactly 18 in-vivo mechanistic coordinates.")
}

selected_run_dir <- function(row) {
  file.path(
    result_root,
    paste0("fit_joint_", row$warmup_label[[1L]]),
    as.character(row$selected_seed[[1L]])
  )
}

read_map_unit <- function(row) {
  path <- file.path(selected_run_dir(row), "best_params_transformed.tsv")
  if (!file.exists(path)) stop("Missing selected MAP table: ", path)
  tab <- utils::read.delim(path, check.names = FALSE)
  par <- setNames(as.numeric(tab$transformed_value), tab$transformed_parameter)
  actual <- setNames(numeric(nrow(coordinate)), coordinate$coordinate)
  ordinary <- coordinate$role == "ordinary"
  actual[ordinary] <- par[coordinate$original_name[ordinary]]
  for (i in which(coordinate$role %in% c("vivo", "vitro"))) {
    center_name <- coordinate$original_name[[i]]
    delta_name <- paste0("delta__", center_name)
    actual[[i]] <- if (coordinate$role[[i]] == "vivo") {
      par[[center_name]] + par[[delta_name]] / 2
    } else {
      par[[center_name]] - par[[delta_name]] / 2
    }
  }
  unit <- (actual[coordinate$coordinate] - coordinate$lower) / coordinate$width
  if (any(!is.finite(unit)) || any(unit < -1e-8) || any(unit > 1 + 1e-8)) {
    stop("Selected MAP is outside the sampling box for ", row$family[[1L]])
  }
  setNames(pmin(pmax(unit, 0), 1), coordinate$coordinate)
}

map_states <- t(vapply(
  seq_len(nrow(selection)),
  function(i) read_map_unit(selection[i, , drop = FALSE]),
  numeric(nrow(coordinate))
))
rownames(map_states) <- families
colnames(map_states) <- coordinate$coordinate
family_centers <- map_states[, family_features, drop = FALSE]

classify_family <- function(unit_matrix) {
  unit_matrix <- as.matrix(unit_matrix)
  if (is.null(colnames(unit_matrix))) {
    if (ncol(unit_matrix) != nrow(coordinate)) {
      stop("Unlabelled family-classification matrix has the wrong dimension.")
    }
    colnames(unit_matrix) <- coordinate$coordinate
  }
  x <- unit_matrix[, family_features, drop = FALSE]
  distances <- vapply(families, function(family) {
    rowSums((x - matrix(
      family_centers[family, ], nrow(x), length(family_features), byrow = TRUE
    ))^2)
  }, numeric(nrow(x)))
  families[max.col(-distances, ties.method = "first")]
}

center_assignments <- classify_family(map_states)
if (!identical(unname(center_assignments), families)) {
  stop("A selected family MAP is not assigned to its own frozen basin.")
}

contract <- list(
  target_version = "family_conditioned_selected_map_voronoi_v1",
  definition = paste0(
    "nearest selected family MAP by squared Euclidean distance in the exact ",
    "bound-normalized 18-dimensional in-vivo mechanistic parameter space"
  ),
  families = families,
  tie_order = families,
  feature_names = family_features,
  coordinate_names = coordinate$coordinate,
  coordinate_md5 = unname(tools::md5sum(coordinate_path)),
  selection_md5 = unname(tools::md5sum(selection_path)),
  centers = family_centers,
  selected_warmup_labels = setNames(selection$warmup_label, selection$family),
  selected_seeds = setNames(as.character(selection$selected_seed), selection$family)
)
save_rds_atomic(contract, contract_rds_path)

center_table <- do.call(rbind, lapply(families, function(family) {
  data.frame(
    target_version = contract$target_version,
    family = family,
    warmup_label = contract$selected_warmup_labels[[family]],
    selected_seed = contract$selected_seeds[[family]],
    feature_order = seq_along(family_features),
    coordinate = family_features,
    center_unit = as.numeric(family_centers[family, ]),
    distance_metric = "squared_euclidean_on_unit_bound_normalized_features",
    tie_order = paste(families, collapse = ","),
    stringsAsFactors = FALSE
  )
}))
center_table$contract_md5 <- unname(tools::md5sum(contract_rds_path))
write_tsv_atomic(center_table, contract_tsv_path)

prior_config <- utils::read.delim(prior_config_path, check.names = FALSE)
parameters <- c(
  "lam_max", "p_mis_base", "p_wgd", "p_misseg", "k_o_mis",
  "O2_crit", "n_O", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
prior_config <- prior_config[match(parameters, prior_config$parameter), , drop = FALSE]
required_prior <- c(
  "parameter", "parameter_order", "transformation",
  "joint_lower_transformed", "joint_upper_transformed",
  "regularization_sigma", "welsch_c", "soft_prior_included",
  "soft_prior_mu", "soft_prior_sd", "lambda_prior"
)
if (anyNA(prior_config$parameter) ||
    !identical(prior_config$parameter, parameters) ||
    !all(required_prior %in% names(prior_config))) {
  stop("The induced-prior configuration does not match the 14 paired parameters.")
}

sample_pair_prior <- function(row, n) {
  vivo <- numeric()
  vitro <- numeric()
  n_candidates <- 0L
  while (length(vivo) < n) {
    remaining <- n - length(vivo)
    batch_n <- max(2000L, as.integer(ceiling(remaining / 0.55)))
    vivo_z <- stats::runif(
      batch_n, row$joint_lower_transformed, row$joint_upper_transformed
    )
    vitro_z <- stats::runif(
      batch_n, row$joint_lower_transformed, row$joint_upper_transformed
    )
    delta <- vivo_z - vitro_z
    penalty <- (row$welsch_c^2 / 2) * (
      1 - exp(-((delta / (row$welsch_c * row$regularization_sigma))^2))
    )
    if (isTRUE(row$soft_prior_included)) {
      penalty <- penalty + row$lambda_prior * 0.5 * (
        (vivo_z - row$soft_prior_mu) / row$soft_prior_sd
      )^2
    }
    keep <- stats::runif(batch_n) <= exp(-penalty)
    n_candidates <- n_candidates + batch_n
    if (any(keep)) {
      vivo <- c(vivo, vivo_z[keep])
      vitro <- c(vitro, vitro_z[keep])
    }
  }
  list(
    vivo = vivo[seq_len(n)],
    vitro = vitro[seq_len(n)],
    candidates = n_candidates
  )
}

draw_conditioned_prior <- function(target_per_family, seed, batch_size = 25000L) {
  set.seed(seed)
  retained <- setNames(lapply(families, function(x) {
    list(
      vivo = matrix(numeric(), nrow = 0L, ncol = length(parameters)),
      vitro = matrix(numeric(), nrow = 0L, ncol = length(parameters))
    )
  }), families)
  for (family in families) {
    colnames(retained[[family]]$vivo) <- parameters
    colnames(retained[[family]]$vitro) <- parameters
  }
  generated <- 0L
  assignment_counts <- setNames(integer(length(families)), families)
  pair_candidate_counts <- setNames(numeric(length(parameters)), parameters)
  max_candidates <- max(2000000L, 100L * target_per_family)

  while (any(vapply(retained, function(x) nrow(x$vivo), integer(1)) < target_per_family)) {
    current <- vapply(retained, function(x) nrow(x$vivo), integer(1))
    n_batch <- min(batch_size, max(target_per_family - current) * 4L)
    n_batch <- max(2000L, as.integer(n_batch))
    vivo <- matrix(NA_real_, nrow = n_batch, ncol = length(parameters),
                   dimnames = list(NULL, parameters))
    vitro <- vivo
    for (i in seq_along(parameters)) {
      sampled <- sample_pair_prior(prior_config[i, , drop = FALSE], n_batch)
      vivo[, i] <- sampled$vivo
      vitro[, i] <- sampled$vitro
      pair_candidate_counts[[i]] <- pair_candidate_counts[[i]] + sampled$candidates
    }

    unit <- matrix(
      NA_real_, nrow = n_batch, ncol = nrow(coordinate),
      dimnames = list(NULL, coordinate$coordinate)
    )
    for (name in ordinary_in_vivo_features) {
      row <- coordinate[coordinate$coordinate == name, , drop = FALSE]
      unit[, name] <- stats::runif(n_batch)
    }
    for (i in seq_along(parameters)) {
      coord_name <- paired_in_vivo_features[
        coordinate$parameter[match(paired_in_vivo_features, coordinate$coordinate)] ==
          parameters[[i]]
      ]
      if (length(coord_name) != 1L) {
        stop("Cannot identify the in-vivo basin coordinate for ", parameters[[i]])
      }
      row <- coordinate[coordinate$coordinate == coord_name, , drop = FALSE]
      unit[, coord_name] <- (vivo[, i] - row$lower) / row$width
    }
    assignment <- classify_family(unit)
    generated <- generated + n_batch
    assignment_counts <- assignment_counts + table(factor(assignment, levels = families))

    for (family in families) {
      need <- target_per_family - nrow(retained[[family]]$vivo)
      if (need <= 0L) next
      idx <- which(assignment == family)
      if (!length(idx)) next
      idx <- head(idx, need)
      retained[[family]]$vivo <- rbind(retained[[family]]$vivo, vivo[idx, , drop = FALSE])
      retained[[family]]$vitro <- rbind(retained[[family]]$vitro, vitro[idx, , drop = FALSE])
    }
    if (generated > max_candidates) {
      stop("Family-conditioned prior rejection exceeded its candidate safety limit.")
    }
  }
  list(
    retained = retained,
    generated_candidates = generated,
    assignment_counts = assignment_counts,
    pair_candidate_counts = pair_candidate_counts
  )
}

primary_n <- as.integer(Sys.getenv("FIGURE5F_FAMILY_PRIOR_DRAWS", "100000"))
secondary_n <- as.integer(Sys.getenv("FIGURE5F_FAMILY_PRIOR_QA_DRAWS", "20000"))
if (!is.finite(primary_n) || primary_n < 1000L ||
    !is.finite(secondary_n) || secondary_n < 1000L) {
  stop("Family-conditioned prior draw counts are invalid.")
}

message("Generating primary family-conditioned induced prior...")
primary <- draw_conditioned_prior(primary_n, seed = 20260817L)
message("Generating independent QA replicate...")
secondary <- draw_conditioned_prior(secondary_n, seed = 20260917L)

draw_rows <- vector("list", length(families) * length(parameters))
config_rows <- vector("list", length(draw_rows))
row_i <- 1L
for (family in families) {
  for (parameter_i in seq_along(parameters)) {
    parameter <- parameters[[parameter_i]]
    vivo <- primary$retained[[family]]$vivo[, parameter_i]
    vitro <- primary$retained[[family]]$vitro[, parameter_i]
    vivo_qa <- secondary$retained[[family]]$vivo[, parameter_i]
    vitro_qa <- secondary$retained[[family]]$vitro[, parameter_i]
    transformation <- prior_config$transformation[[parameter_i]]
    ratio <- if (transformation == "log10") {
      (vivo - vitro) * log2(10)
    } else {
      log2(vivo / vitro)
    }
    ratio_qa <- if (transformation == "log10") {
      (vivo_qa - vitro_qa) * log2(10)
    } else {
      log2(vivo_qa / vitro_qa)
    }
    probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
    q <- stats::quantile(ratio, probs, names = FALSE, type = 8)
    q_qa <- stats::quantile(ratio_qa, probs, names = FALSE, type = 8)
    scaled_difference <- max(abs(q - q_qa)) / max(q[[5L]] - q[[1L]], 1e-12)
    draw_rows[[row_i]] <- data.frame(
      family = family,
      parameter = parameter,
      parameter_order = prior_config$parameter_order[[parameter_i]],
      draw_id = seq_len(primary_n),
      vivo_transformed = vivo,
      vitro_transformed = vitro,
      log2_ratio = ratio,
      stringsAsFactors = FALSE
    )
    config_rows[[row_i]] <- data.frame(
      target_version = contract$target_version,
      family = family,
      parameter = parameter,
      parameter_order = prior_config$parameter_order[[parameter_i]],
      transformation = transformation,
      retained_draws = primary_n,
      qa_draws = secondary_n,
      max_scaled_qa_quantile_difference = scaled_difference,
      primary_total_joint_candidates = primary$generated_candidates,
      primary_family_assignment_fraction =
        primary$assignment_counts[[family]] / primary$generated_candidates,
      contract_md5 = unname(tools::md5sum(contract_rds_path)),
      source_unconditioned_prior_config_md5 =
        unname(tools::md5sum(prior_config_path)),
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
}

draws <- do.call(rbind, draw_rows)
config <- do.call(rbind, config_rows)
if (nrow(draws) != length(families) * length(parameters) * primary_n ||
    any(!is.finite(draws$log2_ratio)) ||
    any(config$max_scaled_qa_quantile_difference > 0.05)) {
  failing <- config[
    !is.finite(config$max_scaled_qa_quantile_difference) |
      config$max_scaled_qa_quantile_difference > 0.05,
    c("family", "parameter", "max_scaled_qa_quantile_difference"),
    drop = FALSE
  ]
  stop(
    "Family-conditioned prior generation failed its numerical QA contract.\n",
    paste(capture.output(print(failing, row.names = FALSE)), collapse = "\n")
  )
}
write_tsv_atomic(draws, prior_draw_path)
config$prior_draw_md5 <- unname(tools::md5sum(prior_draw_path))
write_tsv_atomic(config, prior_config_output)

cat("Saved frozen Figure 5F family-conditioning contract and priors:\n")
cat(contract_rds_path, "\n", contract_tsv_path, "\n", prior_draw_path, "\n", sep = "")
print(unique(config[, c(
  "family", "primary_total_joint_candidates",
  "primary_family_assignment_fraction", "max_scaled_qa_quantile_difference",
  "contract_md5", "prior_draw_md5"
)]))
