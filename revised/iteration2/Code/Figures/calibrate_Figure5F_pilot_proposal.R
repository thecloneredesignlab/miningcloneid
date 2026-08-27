#!/usr/bin/env Rscript

# Derive a fixed proposal preconditioner from the failed Figure 5F pilot 2.
# These matrices affect proposal efficiency only. Pilot draws never enter the
# reported generalized posterior, and adaptation is frozen before production.

resolve_own_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) stop("Unable to resolve this script path.")
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

script_path <- resolve_own_file()
iteration_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
pilot_dir <- file.path(figure5_dir, "generalized_posterior", "pilot2_map_offset")
pilot_state_path <- Sys.getenv(
  "FIGURE5F_PILOT_STATE",
  unset = file.path(pilot_dir, "figure5f_generalized_posterior_state_draws.rds")
)
selection_path <- file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv")
output_rds <- file.path(
  figure5_dir, "generalized_posterior", "figure5f_pilot2_proposal_covariance.rds"
)
output_tsv <- file.path(
  figure5_dir, "generalized_posterior", "figure5f_pilot2_proposal_covariance_diagnostics.tsv"
)
result_root <- Sys.getenv(
  "FIGURE5F_RESULT_ROOT",
  unset = paste0(
    "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/",
    "fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540"
  )
)
result_root <- normalizePath(result_root, mustWork = TRUE)

required <- c(pilot_state_path, selection_path)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing proposal-calibration input(s):\n", paste(missing, collapse = "\n"))

pilot <- readRDS(pilot_state_path)
coord <- pilot$coordinate_metadata
records <- pilot$state_records
if (!is.data.frame(coord) || nrow(coord) != 40L || !length(records)) {
  stop("Pilot 2 state archive does not satisfy the 40-dimensional contract.")
}

selection <- utils::read.delim(selection_path, check.names = FALSE, stringsAsFactors = FALSE)
selection <- selection[selection$selected_for_figure5f %in% TRUE, , drop = FALSE]
selection <- selection[order(match(selection$family, c("C01", "C02", "C03"))), , drop = FALSE]
if (nrow(selection) != 3L || !identical(selection$family, c("C01", "C02", "C03"))) {
  stop("Expected one selected MAP for each of C01, C02, and C03.")
}

read_map_unit <- function(row) {
  run_dir <- file.path(
    result_root,
    paste0("fit_joint_", row$warmup_label[[1L]]),
    as.character(row$selected_seed[[1L]])
  )
  path <- file.path(run_dir, "best_params_transformed.tsv")
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  par <- setNames(as.numeric(tab$transformed_value), tab$transformed_parameter)
  actual <- setNames(numeric(nrow(coord)), coord$coordinate)
  ordinary <- coord$role == "ordinary"
  actual[ordinary] <- par[coord$original_name[ordinary]]
  for (i in which(coord$role != "ordinary")) {
    center_name <- coord$original_name[[i]]
    delta_name <- paste0("delta__", center_name)
    center <- par[[center_name]]
    delta <- par[[delta_name]]
    actual[[i]] <- if (coord$role[[i]] == "vivo") center + delta / 2 else center - delta / 2
  }
  unit <- (actual - coord$lower) / coord$width
  if (any(!is.finite(unit)) || any(unit < -1e-8) || any(unit > 1 + 1e-8)) {
    stop("Selected MAP could not be mapped into the pilot coordinate box for ", row$family[[1L]])
  }
  setNames(pmin(pmax(unit, 0), 1), coord$coordinate)
}

map_states <- t(vapply(seq_len(nrow(selection)), function(i) {
  read_map_unit(selection[i, , drop = FALSE])
}, numeric(nrow(coord))))
rownames(map_states) <- selection$family
colnames(map_states) <- coord$coordinate

nearest_map <- function(x) {
  distances <- vapply(
    seq_len(nrow(map_states)),
    function(i) rowSums((x - matrix(map_states[i, ], nrow(x), ncol(x), byrow = TRUE))^2),
    numeric(nrow(x))
  )
  colnames(distances) <- rownames(map_states)
  colnames(distances)[max.col(-distances, ties.method = "first")]
}

temperatures <- sort(unique(vapply(records, function(x) as.numeric(x$temperature), numeric(1))))
covariance <- list()
diagnostic_rows <- list()
assignment_rows <- list()
row_i <- 1L
assign_i <- 1L
shrinkage <- 0.20

for (temperature in temperatures) {
  selected <- records[vapply(
    records,
    function(x) isTRUE(all.equal(as.numeric(x$temperature), temperature)),
    logical(1)
  )]
  residual_blocks <- list()
  block_i <- 1L
  for (record in selected) {
    unit <- record$unit
    if (!is.matrix(unit) || ncol(unit) != nrow(coord)) stop("Malformed pilot state record.")
    mode <- nearest_map(unit)
    for (family in rownames(map_states)) {
      idx <- which(mode == family)
      assignment_rows[[assign_i]] <- data.frame(
        temperature = temperature,
        chain = record$chain,
        nearest_map = family,
        n_states = length(idx),
        stringsAsFactors = FALSE
      )
      assign_i <- assign_i + 1L
      if (length(idx) < 5L) next
      block <- unit[idx, , drop = FALSE]
      # Remove each visited local neighborhood's location. This retains local
      # correlation geometry without turning between-mode displacement into a
      # random-walk scale; explicit symmetric MAP-offset moves handle that.
      center <- apply(block, 2L, stats::median)
      residual_blocks[[block_i]] <- sweep(block, 2L, center, "-")
      block_i <- block_i + 1L
    }
  }
  residual <- do.call(rbind, residual_blocks)
  if (is.null(residual) || nrow(residual) < 5L * nrow(coord)) {
    stop("Insufficient pilot residuals at T = ", temperature)
  }
  raw_cov <- stats::cov(residual)
  diagonal_target <- diag(diag(raw_cov), nrow(raw_cov))
  stabilized <- (1 - shrinkage) * raw_cov + shrinkage * diagonal_target
  ridge <- max(diag(stabilized)) * 1e-6
  stabilized <- stabilized + diag(ridge, nrow(stabilized))
  dimnames(stabilized) <- list(coord$coordinate, coord$coordinate)
  eigenvalues <- eigen(stabilized, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eigenvalues)) || min(eigenvalues) <= 0) {
    stop("Proposal covariance is not positive definite at T = ", temperature)
  }
  covariance[[as.character(temperature)]] <- stabilized
  diagnostic_rows[[row_i]] <- data.frame(
    temperature = temperature,
    n_residual_states = nrow(residual),
    n_coordinates = ncol(residual),
    shrinkage_to_diagonal = shrinkage,
    ridge = ridge,
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    condition_number = max(eigenvalues) / min(eigenvalues),
    source_pilot_safe_for_primary_figure = FALSE,
    role = "fixed_proposal_preconditioner_only",
    stringsAsFactors = FALSE
  )
  row_i <- row_i + 1L
}

object <- list(
  covariance = covariance,
  coordinate_names = coord$coordinate,
  temperatures = temperatures,
  selected_map_states = map_states,
  shrinkage_to_diagonal = shrinkage,
  source_pilot = normalizePath(pilot_state_path, mustWork = TRUE),
  source_pilot_md5 = unname(tools::md5sum(pilot_state_path)),
  selection_md5 = unname(tools::md5sum(selection_path)),
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  scientific_draws_included = FALSE
)

temporary <- paste0(output_rds, ".tmp_", Sys.getpid())
saveRDS(object, temporary)
if (!file.rename(temporary, output_rds)) stop("Failed to atomically save proposal covariance.")

diagnostics <- do.call(rbind, diagnostic_rows)
assignments <- do.call(rbind, assignment_rows)
assignments <- assignments[order(assignments$temperature, assignments$chain, assignments$nearest_map), ]
diagnostics$proposal_covariance_md5 <- unname(tools::md5sum(output_rds))
diagnostics$nearest_map_assignment_counts <- vapply(
  diagnostics$temperature,
  function(temp) {
    x <- assignments[assignments$temperature == temp & assignments$n_states > 0, ]
    paste(paste(x$chain, x$nearest_map, x$n_states, sep = ":"), collapse = ";")
  },
  character(1)
)
utils::write.table(
  diagnostics, output_tsv, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)

cat("Saved fixed proposal covariance:", output_rds, "\n")
print(diagnostics)
