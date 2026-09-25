#!/usr/bin/env Rscript

# Diagnose and release exactly one Figure 5F family-conditioned posterior.
# Convergence is assessed from that family's independent R1/R2 ladders only.

options(stringsAsFactors = FALSE, warn = 1)

resolve_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Unable to resolve script path.")
  dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE))
}

parse_cli <- function(x) {
  out <- list()
  for (arg in x) {
    if (!startsWith(arg, "--")) next
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[bits[[1L]]]] <- if (length(bits) > 1L) {
      paste(bits[-1L], collapse = "=")
    } else {
      "TRUE"
    }
  }
  out
}

as_integer <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop("--", name, " must be a positive integer.")
  }
  value
}

write_tsv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  utils::write.table(
    x, temporary, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
  if (!file.rename(temporary, path)) stop("Failed to atomically write: ", path)
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(x, temporary, version = 3)
  if (!file.rename(temporary, path)) stop("Failed to atomically write: ", path)
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
family <- toupper(trimws(if (is.null(cli$family)) "" else cli$family))
if (!family %in% c("C01", "C02", "C03")) {
  stop("--family must be C01, C02, or C03.")
}
n_iter <- as_integer(if (is.null(cli$n_iter)) "" else cli$n_iter, "n_iter")
burnin <- as_integer(if (is.null(cli$burnin)) "3000" else cli$burnin, "burnin")
thin <- as_integer(if (is.null(cli$thin)) "1" else cli$thin, "thin")
if (burnin >= n_iter) stop("burnin must be smaller than n_iter.")
n_save <- floor((n_iter - burnin) / thin)
if (n_save < 100L) stop("Fewer than 100 post-warm-up draws are available per chain.")

script_dir <- resolve_script_dir()
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
legacy_output_dir <- file.path(figure5_dir, "generalized_posterior")
output_dir <- file.path(figure5_dir, "generalized_posterior_family_conditioned")
checkpoint_dir <- file.path(
  iteration_root, "tmp", "figure5f_generalized_posterior_family_conditioned"
)
coordinate_path <- file.path(
  legacy_output_dir, "figure5f_sampling_coordinate_contract.tsv"
)
contract_path <- file.path(output_dir, "figure5f_family_basin_contract.rds")
prior_path <- file.path(output_dir, "figure5f_family_conditioned_prior_draws.tsv")
required <- c(coordinate_path, contract_path, prior_path)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing family aggregation input(s):\n", paste(missing, collapse = "\n"))

coordinate <- utils::read.delim(coordinate_path, check.names = FALSE)
contract <- readRDS(contract_path)
target_version <- "family_conditioned_selected_map_voronoi_v1"
if (nrow(coordinate) != 40L ||
    !identical(coordinate$coordinate, contract$coordinate_names) ||
    !identical(contract$target_version, target_version) ||
    !family %in% contract$families) {
  stop("Coordinate metadata and frozen family contract are inconsistent.")
}

classify_family_matrix <- function(unit) {
  unit <- as.matrix(unit)
  x <- unit[, contract$feature_names, drop = FALSE]
  distance <- vapply(contract$families, function(candidate_family) {
    rowSums((x - matrix(
      contract$centers[candidate_family, ],
      nrow(x), length(contract$feature_names), byrow = TRUE
    ))^2)
  }, numeric(nrow(x)))
  contract$tie_order[max.col(-distance, ties.method = "first")]
}

read_ladder <- function(replicate_id) {
  chain <- paste0(family, "_R", replicate_id)
  path <- file.path(checkpoint_dir, paste0("ladder_", chain, "_complete.rds"))
  if (!file.exists(path)) stop("Missing completed ladder: ", path)
  x <- readRDS(path)
  if (!is.list(x) || !identical(x$family, family) ||
      !identical(as.integer(x$replicate), as.integer(replicate_id)) ||
      !identical(x$chain, chain) ||
      !identical(x$target_version, target_version) ||
      !identical(as.integer(x$completed_iteration), n_iter) ||
      !identical(names(x$retained), c("0.5", "1", "2"))) {
    stop("Completed ladder identity or iteration contract failed for ", chain, ".")
  }
  for (temperature in names(x$retained)) {
    retained <- x$retained[[temperature]]
    if (!is.matrix(retained$unit) || !is.matrix(retained$loss) ||
        nrow(retained$unit) != n_save || nrow(retained$loss) != n_save ||
        !identical(colnames(retained$unit), coordinate$coordinate) ||
        any(!is.finite(retained$unit)) || any(!is.finite(retained$loss))) {
      stop("Completed ladder numerical contract failed for ", chain, " at T=", temperature, ".")
    }
  }
  x
}

ladders <- lapply(1:2, read_ladder)
if (!identical(ladders[[1L]]$config_signature, ladders[[2L]]$config_signature)) {
  stop("R1/R2 configuration signatures differ for ", family, ".")
}

split_chain_matrix <- function(x) {
  x <- as.matrix(x)
  half <- floor(nrow(x) / 2L)
  if (half < 3L || ncol(x) != 2L) stop("Split diagnostics require two adequate chains.")
  cbind(
    x[seq_len(half), , drop = FALSE],
    x[seq.int(nrow(x) - half + 1L, nrow(x)), , drop = FALSE]
  )
}

rank_normalize_matrix <- function(x) {
  flat <- as.vector(x)
  ranks <- rank(flat, ties.method = "average")
  matrix(
    stats::qnorm((ranks - 3 / 8) / (length(ranks) + 1 / 4)),
    nrow = nrow(x), ncol = ncol(x)
  )
}

basic_rhat <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  within <- mean(apply(x, 2L, stats::var))
  if (!is.finite(within) || within <= 0) {
    return(if (stats::sd(as.vector(x)) == 0) 1 else Inf)
  }
  between <- n * stats::var(colMeans(x))
  sqrt((((n - 1) / n) * within + between / n) / within)
}

basic_ess <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- nrow(x)
  m <- ncol(x)
  total <- n * m
  within <- mean(apply(x, 2L, stats::var))
  between <- n * stats::var(colMeans(x))
  variance_plus <- (n - 1) / n * within + between / n
  if (!is.finite(variance_plus) || variance_plus <= 0 || !is.finite(within)) {
    return(if (stats::sd(as.vector(x)) == 0) total else NA_real_)
  }
  autocovariance <- vapply(seq_len(m), function(chain_i) {
    as.numeric(stats::acf(
      x[, chain_i], lag.max = n - 1L,
      type = "covariance", plot = FALSE, demean = TRUE
    )$acf)
  }, numeric(n))
  mean_autocovariance <- rowMeans(autocovariance)
  rho <- 1 - (within - mean_autocovariance) / variance_plus
  rho[[1L]] <- 1
  pair_sums <- numeric()
  for (left in seq.int(1L, length(rho) - 1L, by = 2L)) {
    value <- rho[[left]] + rho[[left + 1L]]
    if (!is.finite(value) || value < 0) break
    pair_sums <- c(pair_sums, value)
  }
  if (!length(pair_sums)) return(min(total, m))
  pair_sums <- cummin(pair_sums)
  tau <- -1 + 2 * sum(pair_sums)
  if (!is.finite(tau) || tau <= 0) return(total)
  min(total, total / tau)
}

summarize_matrix <- function(x, scope, variables) {
  rows <- lapply(seq_along(variables), function(variable_i) {
    original <- x[, , variable_i, drop = TRUE]
    if (is.null(dim(original))) original <- matrix(original, ncol = 2L)
    split <- split_chain_matrix(original)
    rank_normalized <- rank_normalize_matrix(split)
    folded <- rank_normalize_matrix(abs(split - stats::median(as.vector(split))))
    tails <- stats::quantile(as.vector(split), c(0.05, 0.95), names = FALSE, type = 8)
    ess_bulk <- basic_ess(rank_normalized)
    ess_tail <- min(
      basic_ess(split <= tails[[1L]]),
      basic_ess(split >= tails[[2L]])
    )
    values <- as.vector(original)
    data.frame(
      family = family,
      temperature = 1,
      scope = scope,
      variable = variables[[variable_i]],
      mean = mean(values),
      sd = stats::sd(values),
      rhat = max(basic_rhat(rank_normalized), basic_rhat(folded)),
      ess_bulk = ess_bulk,
      ess_tail = ess_tail,
      mcse_mean = stats::sd(values) / sqrt(ess_bulk),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

parameters <- unique(coordinate$parameter[coordinate$role == "vivo"])
if (length(parameters) != 14L) stop("Expected 14 paired parameters.")
state_array <- array(
  NA_real_, dim = c(n_save, 2L, nrow(coordinate)),
  dimnames = list(NULL, vapply(ladders, `[[`, character(1), "chain"), coordinate$coordinate)
)
ratio_array <- array(
  NA_real_, dim = c(n_save, 2L, length(parameters)),
  dimnames = list(NULL, vapply(ladders, `[[`, character(1), "chain"), parameters)
)

ratios_from_unit <- function(unit) {
  actual <- sweep(unit, 2L, coordinate$width, `*`)
  actual <- sweep(actual, 2L, coordinate$lower, `+`)
  colnames(actual) <- coordinate$coordinate
  ratio <- matrix(
    NA_real_, nrow = nrow(unit), ncol = length(parameters),
    dimnames = list(NULL, parameters)
  )
  for (parameter in parameters) {
    vivo_row <- coordinate[coordinate$role == "vivo" & coordinate$parameter == parameter, ]
    vitro_row <- coordinate[coordinate$role == "vitro" & coordinate$parameter == parameter, ]
    if (nrow(vivo_row) != 1L || nrow(vitro_row) != 1L ||
        !identical(vivo_row$transform, vitro_row$transform)) {
      stop("Malformed coordinate pair for ", parameter, ".")
    }
    vivo <- actual[, vivo_row$coordinate]
    vitro <- actual[, vitro_row$coordinate]
    ratio[, parameter] <- if (identical(vivo_row$transform, "log10")) {
      (vivo - vitro) * log2(10)
    } else {
      log2(vivo / vitro)
    }
  }
  ratio
}

for (chain_i in seq_along(ladders)) {
  unit <- ladders[[chain_i]]$retained[["1"]]$unit
  state_array[, chain_i, ] <- unit
  ratio_array[, chain_i, ] <- ratios_from_unit(unit)
}

diagnostics <- rbind(
  summarize_matrix(ratio_array, "paired_log2_ratios", parameters),
  summarize_matrix(state_array, "all_40_sampling_coordinates", coordinate$coordinate)
)
ratio_diagnostics <- diagnostics[diagnostics$scope == "paired_log2_ratios", ]
state_diagnostics <- diagnostics[diagnostics$scope == "all_40_sampling_coordinates", ]

sampler_rows <- list()
sampler_i <- 1L
for (ladder in ladders) {
  for (temperature_i in seq_along(ladder$temperatures)) {
    for (kernel in colnames(ladder$local_attempts)) {
      attempts <- ladder$local_attempts[temperature_i, kernel]
      accepts <- ladder$local_accepts[temperature_i, kernel]
      sampler_rows[[sampler_i]] <- data.frame(
        family = family, chain = ladder$chain,
        temperature = as.character(ladder$temperatures[[temperature_i]]),
        diagnostic = "local_proposal", kernel = kernel,
        attempts = attempts, accepts = accepts,
        acceptance_rate = accepts / max(1, attempts),
        final_scale = ladder$final_scales[temperature_i, kernel],
        basin_boundary_rejections = ladder$basin_rejections[[temperature_i]],
        stringsAsFactors = FALSE
      )
      sampler_i <- sampler_i + 1L
    }
  }
  for (left in seq_len(length(ladder$temperatures) - 1L)) {
    attempts <- ladder$swap_attempts[[left]]
    accepts <- ladder$swap_accepts[[left]]
    sampler_rows[[sampler_i]] <- data.frame(
      family = family, chain = ladder$chain,
      temperature = paste0(ladder$temperatures[[left]], "<->", ladder$temperatures[[left + 1L]]),
      diagnostic = "replica_swap", kernel = "adjacent",
      attempts = attempts, accepts = accepts,
      acceptance_rate = accepts / max(1, attempts),
      final_scale = NA_real_, basin_boundary_rejections = NA_integer_,
      stringsAsFactors = FALSE
    )
    sampler_i <- sampler_i + 1L
  }
}
sampler_diagnostics <- do.call(rbind, sampler_rows)

all_assignments_correct <- all(vapply(ladders, function(ladder) {
  all(vapply(ladder$retained, function(retained) {
    all(classify_family_matrix(retained$unit) == family)
  }, logical(1)))
}, logical(1)))
local_acceptance <- sampler_diagnostics$acceptance_rate[
  sampler_diagnostics$diagnostic == "local_proposal"
]
swap_acceptance <- sampler_diagnostics$acceptance_rate[
  sampler_diagnostics$diagnostic == "replica_swap"
]

checks <- data.frame(
  check = c(
    "family_ratio_rhat_max_le_1p01",
    "family_ratio_ess_bulk_min_ge_400",
    "family_ratio_ess_tail_min_ge_400",
    "family_state_rhat_max_le_1p05",
    "family_state_ess_bulk_min_ge_100",
    "all_retained_states_in_frozen_family_basin",
    "local_acceptance_finite",
    "adjacent_swap_acceptance_min_ge_0p05"
  ),
  observed = c(
    max(ratio_diagnostics$rhat),
    min(ratio_diagnostics$ess_bulk),
    min(ratio_diagnostics$ess_tail),
    max(state_diagnostics$rhat),
    min(state_diagnostics$ess_bulk),
    as.numeric(all_assignments_correct),
    as.numeric(length(local_acceptance) > 0L && all(is.finite(local_acceptance))),
    min(swap_acceptance)
  ),
  threshold = c("<=1.01", ">=400", ">=400", "<=1.05", ">=100", "TRUE", "TRUE", ">=0.05"),
  passed = c(
    all(is.finite(ratio_diagnostics$rhat)) && all(ratio_diagnostics$rhat <= 1.01),
    all(is.finite(ratio_diagnostics$ess_bulk)) && all(ratio_diagnostics$ess_bulk >= 400),
    all(is.finite(ratio_diagnostics$ess_tail)) && all(ratio_diagnostics$ess_tail >= 400),
    all(is.finite(state_diagnostics$rhat)) && all(state_diagnostics$rhat <= 1.05),
    all(is.finite(state_diagnostics$ess_bulk)) && all(state_diagnostics$ess_bulk >= 100),
    all_assignments_correct,
    length(local_acceptance) > 0L && all(is.finite(local_acceptance)),
    length(swap_acceptance) > 0L && all(is.finite(swap_acceptance)) &&
      all(swap_acceptance >= 0.05)
  ),
  stringsAsFactors = FALSE
)
checks$family <- family
checks$n_iter <- n_iter
checks$retained_draws_per_chain <- n_save
checks$target_version <- target_version
checks$config_signature <- ladders[[1L]]$config_signature
checks$checked_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
passed <- all(checks$passed)

prefix <- file.path(output_dir, paste0("figure5f_", family))
write_tsv_atomic(diagnostics, paste0(prefix, "_diagnostics.tsv"))
write_tsv_atomic(sampler_diagnostics, paste0(prefix, "_sampler_diagnostics.tsv"))
write_tsv_atomic(checks, paste0(prefix, "_readiness.tsv"))
save_rds_atomic(
  list(
    family = family, n_iter = n_iter, burnin = burnin, thin = thin,
    target_version = target_version,
    config_signature = ladders[[1L]]$config_signature,
    passed = passed, diagnostics = diagnostics, readiness = checks
  ),
  paste0(prefix, "_diagnostic_state.rds")
)

quantile_summary <- function(x) {
  q <- stats::quantile(x, c(0.05, 0.25, 0.50, 0.75, 0.95), names = FALSE, type = 8)
  setNames(q, c("q05", "q25", "median", "q75", "q95"))
}

if (passed) {
  sensitivity_diagnostics <- list()
  sensitivity_i <- 1L
  for (temperature in c(0.5, 2)) {
    key <- as.character(temperature)
    sensitivity_ratio_array <- array(
      NA_real_, dim = c(n_save, 2L, length(parameters)),
      dimnames = list(NULL, vapply(ladders, `[[`, character(1), "chain"), parameters)
    )
    for (chain_i in seq_along(ladders)) {
      sensitivity_ratio_array[, chain_i, ] <-
        ratios_from_unit(ladders[[chain_i]]$retained[[key]]$unit)
    }
    item <- summarize_matrix(
      sensitivity_ratio_array, "paired_log2_ratios", parameters
    )
    item$temperature <- temperature
    sensitivity_diagnostics[[sensitivity_i]] <- item
    sensitivity_i <- sensitivity_i + 1L
  }
  diagnostics <- rbind(diagnostics, do.call(rbind, sensitivity_diagnostics))
  diagnostics <- diagnostics[order(
    diagnostics$temperature,
    match(diagnostics$scope, c("paired_log2_ratios", "all_40_sampling_coordinates")),
    diagnostics$variable
  ), , drop = FALSE]
  write_tsv_atomic(diagnostics, paste0(prefix, "_diagnostics.tsv"))

  draw_rows <- list()
  row_i <- 1L
  for (ladder in ladders) {
    for (temperature in c(0.5, 1, 2)) {
      key <- as.character(temperature)
      unit <- ladder$retained[[key]]$unit
      loss <- ladder$retained[[key]]$loss
      ratios <- ratios_from_unit(unit)
      actual <- sweep(unit, 2L, coordinate$width, `*`)
      actual <- sweep(actual, 2L, coordinate$lower, `+`)
      colnames(actual) <- coordinate$coordinate
      for (parameter in parameters) {
        vivo_row <- coordinate[coordinate$role == "vivo" & coordinate$parameter == parameter, ]
        vitro_row <- coordinate[coordinate$role == "vitro" & coordinate$parameter == parameter, ]
        vivo_t <- actual[, vivo_row$coordinate]
        vitro_t <- actual[, vitro_row$coordinate]
        vivo_n <- if (identical(vivo_row$transform, "log10")) 10^vivo_t else vivo_t
        vitro_n <- if (identical(vitro_row$transform, "log10")) 10^vitro_t else vitro_t
        draw_rows[[row_i]] <- data.frame(
          family = family,
          warmup_label = ladder$warmup_label,
          selected_seed = ladder$selected_seed,
          chain = ladder$chain,
          replicate = ladder$replicate,
          draw = seq_len(n_save),
          temperature = temperature,
          parameter = parameter,
          vivo_transformed = vivo_t,
          vitro_transformed = vitro_t,
          vivo_natural = vivo_n,
          vitro_natural = vitro_n,
          log2_ratio = ratios[, parameter],
          data_loss = loss[, "data_loss"],
          regularization_loss = loss[, "regularization_loss"],
          target_loss = loss[, "target_loss"],
          target_version = target_version,
          final_family_n_iter = n_iter,
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1L
      }
    }
  }
  draws <- do.call(rbind, draw_rows)
  expected_rows <- 2L * 3L * 14L * n_save
  if (nrow(draws) != expected_rows || any(!is.finite(draws$log2_ratio))) {
    stop("Released family draw table failed its structural contract.")
  }
  write_tsv_atomic(draws, paste0(prefix, "_posterior_draws.tsv"))

  prior <- utils::read.delim(prior_path, check.names = FALSE)
  prior <- prior[prior$family == family, , drop = FALSE]
  summary_groups <- split(draws, interaction(draws$temperature, draws$parameter, drop = TRUE))
  summary_table <- do.call(rbind, lapply(summary_groups, function(group) {
    q <- quantile_summary(group$log2_ratio)
    prior_values <- prior$log2_ratio[prior$parameter == group$parameter[[1L]]]
    prior_q <- quantile_summary(prior_values)
    data.frame(
      family = family,
      temperature = group$temperature[[1L]],
      parameter = group$parameter[[1L]],
      n_draws = nrow(group),
      q05 = q[["q05"]], q25 = q[["q25"]], median = q[["median"]],
      q75 = q[["q75"]], q95 = q[["q95"]],
      width90 = q[["q95"]] - q[["q05"]],
      prior_q05 = prior_q[["q05"]], prior_median = prior_q[["median"]],
      prior_q95 = prior_q[["q95"]],
      prior_width90 = prior_q[["q95"]] - prior_q[["q05"]],
      contraction_ratio90 =
        (q[["q95"]] - q[["q05"]]) / (prior_q[["q95"]] - prior_q[["q05"]]),
      final_family_n_iter = n_iter,
      target_version = target_version,
      stringsAsFactors = FALSE
    )
  }))
  summary_table <- summary_table[order(
    summary_table$temperature, match(summary_table$parameter, parameters)
  ), , drop = FALSE]
  write_tsv_atomic(summary_table, paste0(prefix, "_posterior_summary.tsv"))
}

print(checks[, c("check", "observed", "threshold", "passed")])
message(
  family, " family-conditioned posterior diagnostics at n_iter=", n_iter,
  if (passed) " passed." else " did not pass; only this family should be extended."
)
quit(save = "no", status = 0L, runLast = FALSE)
