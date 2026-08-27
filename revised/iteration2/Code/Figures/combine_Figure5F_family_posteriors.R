#!/usr/bin/env Rscript

# Combine three independently diagnosed family-conditioned posteriors only
# after every family has passed its frozen gates twice consecutively.

options(stringsAsFactors = FALSE, warn = 1)

resolve_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Unable to resolve script path.")
  dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE))
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

script_dir <- resolve_script_dir()
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
output_dir <- file.path(figure5_dir, "generalized_posterior_family_conditioned")
checkpoint_dir <- file.path(
  iteration_root, "tmp", "figure5f_generalized_posterior_family_conditioned"
)
families <- c("C01", "C02", "C03")
target_version <- "family_conditioned_selected_map_voronoi_v1"
contract_path <- file.path(output_dir, "figure5f_family_basin_contract.rds")
prior_path <- file.path(output_dir, "figure5f_family_conditioned_prior_draws.tsv")
required_frozen_inputs <- c(contract_path, prior_path)
missing_frozen_inputs <- required_frozen_inputs[!file.exists(required_frozen_inputs)]
if (length(missing_frozen_inputs)) {
  stop(
    "Missing frozen family-combination input(s):\n",
    paste(missing_frozen_inputs, collapse = "\n")
  )
}

family_objects <- lapply(families, function(family) {
  prefix <- file.path(output_dir, paste0("figure5f_", family))
  paths <- c(
    marker = file.path(checkpoint_dir, paste0("family_", family, "_complete.tsv")),
    readiness = paste0(prefix, "_readiness.tsv"),
    diagnostics = paste0(prefix, "_diagnostics.tsv"),
    sampler = paste0(prefix, "_sampler_diagnostics.tsv"),
    draws = paste0(prefix, "_posterior_draws.tsv"),
    summary = paste0(prefix, "_posterior_summary.tsv")
  )
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing completed-family artifact(s) for ", family, ":\n", paste(missing, collapse = "\n"))
  }
  marker <- utils::read.delim(paths[["marker"]], check.names = FALSE)
  readiness <- utils::read.delim(paths[["readiness"]], check.names = FALSE)
  if (nrow(marker) != 1L || marker$family[[1L]] != family ||
      marker$consecutive_passes[[1L]] < 2L ||
      !all(readiness$family == family) || !all(readiness$passed %in% TRUE) ||
      !all(readiness$target_version == target_version) ||
      !all(readiness$n_iter == marker$n_iter[[1L]])) {
    stop("Frozen completion/readiness contract failed for ", family, ".")
  }
  list(
    family = family,
    n_iter = marker$n_iter[[1L]],
    marker = marker,
    readiness = readiness,
    diagnostics = utils::read.delim(paths[["diagnostics"]], check.names = FALSE),
    sampler = utils::read.delim(paths[["sampler"]], check.names = FALSE),
    draws = utils::read.delim(paths[["draws"]], check.names = FALSE),
    summary = utils::read.delim(paths[["summary"]], check.names = FALSE),
    paths = paths
  )
})
names(family_objects) <- families

family_signatures <- vapply(
  family_objects,
  function(x) unique(as.character(x$readiness$config_signature)),
  character(1)
)
if (length(unique(family_signatures)) != 1L) {
  stop("C01/C02/C03 were not sampled under one frozen configuration signature.")
}

for (family in families) {
  item <- family_objects[[family]]
  if (!all(item$draws$family == family) ||
      !all(item$draws$target_version == target_version) ||
      !setequal(unique(item$draws$temperature), c(0.5, 1, 2)) ||
      !all(item$summary$family == family) ||
      !all(item$diagnostics$family == family)) {
    stop("Released family outputs fail identity checks for ", family, ".")
  }
}

draws <- do.call(rbind, lapply(family_objects, `[[`, "draws"))
diagnostics <- do.call(rbind, lapply(family_objects, `[[`, "diagnostics"))
sampler_diagnostics <- do.call(rbind, lapply(family_objects, `[[`, "sampler"))
family_summary <- do.call(rbind, lapply(family_objects, `[[`, "summary"))
readiness <- do.call(rbind, lapply(family_objects, `[[`, "readiness"))
rownames(draws) <- rownames(diagnostics) <- rownames(sampler_diagnostics) <-
  rownames(family_summary) <- rownames(readiness) <- NULL

primary <- draws[draws$temperature == 1, , drop = FALSE]
if (!nrow(primary) || !all(primary$family %in% families) ||
    any(!is.finite(primary$log2_ratio)) ||
    any(!is.finite(primary$vivo_natural)) ||
    any(!is.finite(primary$vitro_natural))) {
  stop("Combined T=1 family posterior draw contract failed.")
}

configuration <- data.frame(
  key = c(
    "target_name", "target_version", "family_targets_distinct",
    "family_support_definition", "family_tie_order",
    "C01_final_n_iter", "C02_final_n_iter", "C03_final_n_iter",
    "burnin", "thin", "independent_chains_per_family",
    "temperature_ladder", "display_temperatures",
    "dynamic_stop_unit", "required_consecutive_passes",
    "generator_md5", "aggregator_md5", "combiner_md5",
    "contract_md5", "conditioned_prior_md5", "safe_for_primary_figure"
  ),
  value = c(
    paste0(
      "three family-conditioned generalized posteriors: ",
      "exp(-L_data/T) x exp(-L_regularization) x I(theta in family basin)"
    ),
    target_version,
    TRUE,
    paste0(
      "nearest selected family MAP in exact bound-normalized ",
      "18-dimensional in-vivo mechanistic parameter space"
    ),
    paste(families, collapse = ","),
    vapply(family_objects, function(x) as.numeric(x$n_iter), numeric(1)),
    3000, 1, 2,
    "0.5,1,2,4,8", "0.5,1,2",
    "one family R1/R2 pair", 2,
    unname(tools::md5sum(file.path(script_dir, "generate_Figure5F_generalized_posterior.R"))),
    unname(tools::md5sum(file.path(script_dir, "aggregate_Figure5F_family_posterior.R"))),
    unname(tools::md5sum(file.path(script_dir, "combine_Figure5F_family_posteriors.R"))),
    unname(tools::md5sum(contract_path)),
    unname(tools::md5sum(prior_path)),
    TRUE
  ),
  stringsAsFactors = FALSE
)

write_tsv_atomic(
  draws,
  file.path(output_dir, "figure5f_generalized_posterior_draws.tsv")
)
write_tsv_atomic(
  diagnostics,
  file.path(output_dir, "figure5f_generalized_posterior_diagnostics.tsv")
)
write_tsv_atomic(
  sampler_diagnostics,
  file.path(output_dir, "figure5f_sampler_diagnostics.tsv")
)
write_tsv_atomic(
  family_summary,
  file.path(output_dir, "figure5f_generalized_posterior_family_summary.tsv")
)
write_tsv_atomic(
  readiness,
  file.path(output_dir, "figure5f_generalized_posterior_readiness.tsv")
)
write_tsv_atomic(
  configuration,
  file.path(output_dir, "figure5f_generalized_posterior_config.tsv")
)
write_tsv_atomic(primary, file.path(figure5_dir, "figure5f_posterior_draws.tsv"))

release <- data.frame(
  target_version = target_version,
  family = families,
  final_n_iter = vapply(family_objects, function(x) as.numeric(x$n_iter), numeric(1)),
  consecutive_passes = 2L,
  posterior_draw_rows = vapply(family_objects, function(x) nrow(x$draws), integer(1)),
  primary_draw_rows = vapply(family_objects, function(x) sum(x$draws$temperature == 1), integer(1)),
  released_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  stringsAsFactors = FALSE
)
write_tsv_atomic(
  release,
  file.path(output_dir, "figure5f_family_conditioned_release.tsv")
)

print(release)
message("All three distinct family-conditioned posteriors passed and were combined.")
