#!/usr/bin/env Rscript

# Read-only checkpoint monitor for the Figure 5F production sampler. It never
# modifies checkpoints or releases draws. Before warm-up ends it reports only
# progress and proposal diagnostics; afterward it also computes provisional
# T=1 ratio diagnostics from the retained prefix shared by all six ladders.

resolve_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Unable to resolve script path.")
  dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE))
}

script_dir <- resolve_script_dir()
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
checkpoint_dir <- file.path(iteration_root, "tmp", "figure5f_generalized_posterior")
coord_path <- file.path(
  iteration_root, "data", "Figures", "Figure5", "generalized_posterior",
  "figure5f_sampling_coordinate_contract.tsv"
)

paths <- sort(list.files(
  checkpoint_dir,
  pattern = "^ladder_C0[123]_R[12][.]rds$",
  full.names = TRUE
))
if (length(paths) != 6L) stop("Expected six active ladder checkpoints.")
checkpoints <- lapply(paths, readRDS)
progress <- do.call(rbind, lapply(checkpoints, function(x) {
  data.frame(
    chain = x$chain,
    iteration = x$iteration,
    target_iteration = x$target_iteration,
    retained_draws = x$save_index,
    map_jump_attempts = sum(x$local_attempts[, "map_jump"]),
    map_jump_accepts = sum(x$local_accepts[, "map_jump"]),
    minimum_swap_acceptance = min(x$swap_accepts / pmax(1, x$swap_attempts)),
    stringsAsFactors = FALSE
  )
}))
progress <- progress[order(progress$chain), ]
print(progress, row.names = FALSE)

n_common <- min(progress$retained_draws)
if (n_common < 100L) {
  cat("\nProvisional convergence diagnostics withheld: fewer than 100 common retained draws per chain.\n")
  quit(save = "no", status = 0L, runLast = FALSE)
}
if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Package 'posterior' is required for provisional diagnostics.")
}

coord <- utils::read.delim(coord_path, check.names = FALSE, stringsAsFactors = FALSE)
paired <- unique(coord[coord$role %in% c("vivo", "vitro"), c("parameter", "original_name", "transform")])
parameters <- c(
  "lam_max", "p_mis_base", "p_wgd", "p_misseg", "k_o_mis",
  "O2_crit", "n_O", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
paired <- paired[match(parameters, paired$parameter), , drop = FALSE]
if (anyNA(paired$parameter)) stop("Paired-coordinate metadata is incomplete.")

ratio_array <- array(
  NA_real_,
  dim = c(n_common, length(checkpoints), length(parameters)),
  dimnames = list(NULL, vapply(checkpoints, `[[`, character(1), "chain"), parameters)
)
for (chain_i in seq_along(checkpoints)) {
  unit <- checkpoints[[chain_i]]$retained[["1"]]$unit[seq_len(n_common), , drop = FALSE]
  actual <- sweep(unit, 2L, coord$width, "*")
  actual <- sweep(actual, 2L, coord$lower, "+")
  colnames(actual) <- coord$coordinate
  for (parameter_i in seq_along(parameters)) {
    row <- paired[parameter_i, , drop = FALSE]
    vivo <- actual[, paste0("vivo__", row$original_name)]
    vitro <- actual[, paste0("vitro__", row$original_name)]
    ratio_array[, chain_i, parameter_i] <- if (row$transform == "identity") {
      log2(vivo / vitro)
    } else {
      (vivo - vitro) * log2(10)
    }
  }
}

diagnostics <- as.data.frame(posterior::summarise_draws(
  posterior::as_draws_array(ratio_array),
  "rhat", "ess_bulk", "ess_tail"
))
diagnostics <- diagnostics[order(-diagnostics$rhat), ]
cat("\nProvisional T=1 paired-ratio diagnostics (not a release decision):\n")
print(diagnostics, row.names = FALSE)
cat(
  "\nmax R-hat =", max(diagnostics$rhat),
  "; min bulk ESS =", min(diagnostics$ess_bulk),
  "; min tail ESS =", min(diagnostics$ess_tail), "\n"
)
