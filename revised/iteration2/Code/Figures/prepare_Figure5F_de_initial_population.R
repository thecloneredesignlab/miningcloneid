#!/usr/bin/env Rscript

# Reconstruct the exact differential-evolution (DE) initial populations used by
# the selected C01/C02/C03 joint warm-start pairs. The plotted reference is the
# optimizer's actual starting-population distribution, not a Bayesian prior
# and not a distribution induced by objective penalties.

options(stringsAsFactors = FALSE, warn = 1)

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    dirname(normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE))
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
})
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")

source_result_root <- normalizePath(
  "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540",
  mustWork = TRUE
)
source_repo_root <- normalizePath(
  "/Users/4482173/Documents/GitHub/soft_coupling",
  mustWork = TRUE
)
code_snapshot_root <- normalizePath(
  file.path(iteration_root, "tmp", "miningcloneid_83953a8"),
  mustWork = TRUE
)
analysis_script <- file.path(
  code_snapshot_root,
  "oxygen", "code", "O2_supply_demand_MAP", "analysis", "multi_warmup",
  "analyze_warm_up_joint_results.R"
)
backend_script <- file.path(
  code_snapshot_root,
  "oxygen", "code", "O2_supply_demand_MAP", "util",
  "o2_supply_demand_map_fit_joint_backend.R"
)
selection_path <- file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv")

required_paths <- c(analysis_script, backend_script, selection_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop(
    "Missing Figure 5F DE-initialization input(s):\n",
    paste(missing_paths, collapse = "\n")
  )
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

families <- c("C01", "C02", "C03")
parameters <- c(
  "lam_max", "p_mis_base", "p_wgd",
  "p_misseg", "k_o_mis",
  "O2_crit", "n_O",
  "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
  "buffer_smax", "buffer_beta", "buffer_n_exp"
)
expected_selection <- c(
  C01 = "tsne_vi_seed366_C01Sc01_vt_seed10",
  C02 = "tsne_vi_seed25_C02Sc01_vt_seed10",
  C03 = "tsne_vi_seed311_C03Sc02_vt_seed10"
)

selection <- read_tsv(selection_path)
if (!all(c(
  "family", "warmup_label", "selected_for_figure5f",
  "invivo_seed", "invitro_seed"
) %in% names(selection))) {
  stop("Figure 5F selected-pair table has an incomplete schema.")
}
selected <- selection[as.logical(selection$selected_for_figure5f), , drop = FALSE]
selected <- selected[order(match(selected$family, families)), , drop = FALSE]
observed_selection <- setNames(selected$warmup_label, selected$family)
if (nrow(selected) != 3L ||
    !identical(selected$family, families) ||
    !identical(observed_selection[families], expected_selection) ||
    any(selected$invitro_seed != "seed10")) {
  stop("Figure 5F selected-pair contract does not match the approved inputs.")
}

# Load the exact code snapshot identified for these fits.  The analysis helper
# provides the same deterministic reconstruction route used by the original
# landscape analysis: build one context per pair, set each joint seed, and call
# joint_deoptim_initial_population().
analysis_env <- new.env(parent = globalenv())
analysis_env$commandArgs <- function(trailingOnly = FALSE) {
  if (isTRUE(trailingOnly)) character(0) else paste0("--file=", analysis_script)
}
sys.source(analysis_script, envir = analysis_env, chdir = TRUE)
joint_env <- analysis_env$load_joint_backend_env()

hpc_repo_prefix <-
  "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling"

seed_number <- function(path) {
  suppressWarnings(as.integer(sub("^seed", "", basename(path))))
}

pair_seed_dirs <- function(pair_dir) {
  dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^seed[0-9]+$", basename(dirs))]
  dirs <- dirs[file.exists(file.path(dirs, "best_params.tsv"))]
  dirs[order(seed_number(dirs))]
}

transformation_for_center <- function(center_name) {
  if (startsWith(center_name, "log10_")) "log10" else "identity"
}

ratio_from_center_delta <- function(center, delta, transformation) {
  if (identical(transformation, "log10")) {
    return(delta * log2(10))
  }
  vivo <- center + delta / 2
  vitro <- center - delta / 2
  out <- log2(vivo / vitro)
  out[!is.finite(out) | vivo <= 0 | vitro <= 0] <- NA_real_
  out
}

population_rows <- vector("list", length(families))
context_value_rows <- vector("list", length(families))
config_rows <- list()
provenance_rows <- list()
config_index <- 1L
provenance_index <- 1L
requested_cores <- suppressWarnings(as.integer(
  Sys.getenv("FIGURE5F_DE_CORES", "4")
))
if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 1L
replay_cores <- min(requested_cores, 8L)

for (family_index in seq_along(families)) {
  family <- families[[family_index]]
  selected_row <- selected[selected$family == family, , drop = FALSE]
  warmup_label <- selected_row$warmup_label[[1L]]
  pair_dir <- file.path(source_result_root, paste0("fit_joint_", warmup_label))
  if (!dir.exists(pair_dir)) stop("Missing selected joint pair: ", pair_dir)

  seed_dirs <- pair_seed_dirs(pair_dir)
  joint_seeds <- seed_number(seed_dirs)
  if (length(seed_dirs) != 500L || !identical(joint_seeds, 1:500)) {
    stop("Selected pair does not contain the complete seed1--seed500 set: ", pair_dir)
  }

  # The run snapshots are the authoritative parameter-table inputs.  Override
  # mutable repository paths in the effective command with the per-run copies.
  seed1_dir <- seed_dirs[[1L]]
  argv <- analysis_env$read_effective_argv(
    seed1_dir,
    path_map_from = hpc_repo_prefix,
    path_map_to = source_repo_root
  )
  argv$parameter_table <- file.path(seed1_dir, "parameter_table_input_invivo.csv")
  argv$invitro_parameter_table <-
    file.path(seed1_dir, "parameter_table_input_invitro.csv")
  argv$joint_soft_coupling_parameters_table <-
    file.path(seed1_dir, "joint_soft_coupling_parameters_table_input.csv")

  context <- joint_env$build_joint_context(argv)
  NP_used <- max(
    as.integer(context$NP),
    as.integer(context$joint_np_min_factor) * length(context$init)
  )
  sigmaN <- as.numeric(context$joint_warmup$sigmaN)
  metadata <- context$joint_soft_coupling$metadata
  metadata <- metadata[match(parameters, metadata$parameter), , drop = FALSE]
  if (!isTRUE(context$joint_warmup$enabled) ||
      NP_used != 400L ||
      !isTRUE(all.equal(sigmaN, 0.1216, tolerance = 1e-12)) ||
      length(context$init) != 40L ||
      anyNA(metadata$parameter) ||
      !identical(metadata$parameter, parameters)) {
    stop("Unexpected DE initialization contract for ", family, ": ", warmup_label)
  }

  n_rows <- length(joint_seeds) * NP_used
  row_seed <- rep(joint_seeds, each = NP_used)
  row_member <- rep(seq_len(NP_used), times = length(joint_seeds))

  seed_population_blocks <- parallel::mclapply(seq_along(joint_seeds), function(seed_index) {
    current_seed <- joint_seeds[[seed_index]]
    seed_context <- context
    seed_context$seed <- current_seed
    seed_context$raw$seed <- as.character(current_seed)
    set.seed(current_seed)
    population <- joint_env$joint_deoptim_initial_population(
      seed_context, NP_used
    )
    ratio_block <- matrix(
      NA_real_, nrow = NP_used, ncol = length(parameters),
      dimnames = list(NULL, parameters)
    )
    vivo_block <- ratio_block
    vitro_block <- ratio_block
    for (parameter_index in seq_along(parameters)) {
      meta <- metadata[parameter_index, , drop = FALSE]
      center_name <- meta$center_name[[1L]]
      delta_name <- meta$delta_name[[1L]]
      transformation <- transformation_for_center(center_name)
      center <- population[, center_name]
      delta <- population[, delta_name]
      vivo_transformed <- center + delta / 2
      vitro_transformed <- center - delta / 2
      ratio_block[, parameter_index] <- ratio_from_center_delta(
        center, delta, transformation
      )
      if (identical(transformation, "log10")) {
        vivo_block[, parameter_index] <- 10^vivo_transformed
        vitro_block[, parameter_index] <- 10^vitro_transformed
      } else {
        vivo_block[, parameter_index] <- vivo_transformed
        vitro_block[, parameter_index] <- vitro_transformed
      }
    }
    list(ratio = ratio_block, vivo = vivo_block, vitro = vitro_block)
  }, mc.cores = replay_cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
  if (any(vapply(seed_population_blocks, inherits, logical(1), what = "try-error"))) {
    stop("Parallel DE initial-population replay failed for ", family)
  }
  ratios <- do.call(rbind, lapply(seed_population_blocks, `[[`, "ratio"))
  vivo_values <- do.call(rbind, lapply(seed_population_blocks, `[[`, "vivo"))
  vitro_values <- do.call(rbind, lapply(seed_population_blocks, `[[`, "vitro"))

  if (any(!is.finite(ratios)) ||
      any(!is.finite(vivo_values)) ||
      any(!is.finite(vitro_values)) ||
      any(vivo_values <= 0) || any(vitro_values <= 0)) {
    bad <- which(!is.finite(ratios), arr.ind = TRUE)
    if (!nrow(bad)) bad <- matrix(c(1L, 1L), nrow = 1L)
    stop(
      "Invalid reconstructed DE context value or ratio for ", family, ", row ",
      bad[[1L]], ", parameter ", colnames(ratios)[bad[[2L]]]
    )
  }

  population_rows[[family_index]] <- data.frame(
    family = family,
    warmup_label = warmup_label,
    joint_seed = row_seed,
    initial_member = row_member,
    exact_warm_start = row_member == 1L,
    as.data.frame(ratios, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  context_value_rows[[family_index]] <- data.frame(
    family = family,
    warmup_label = warmup_label,
    joint_seed = row_seed,
    initial_member = row_member,
    exact_warm_start = row_member == 1L,
    setNames(as.data.frame(vivo_values, check.names = FALSE), paste0("vivo__", parameters)),
    setNames(as.data.frame(vitro_values, check.names = FALSE), paste0("vitro__", parameters)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  for (parameter_index in seq_along(parameters)) {
    meta <- metadata[parameter_index, , drop = FALSE]
    center_name <- meta$center_name[[1L]]
    delta_name <- meta$delta_name[[1L]]
    transformation <- transformation_for_center(center_name)
    center_init <- as.numeric(context$init[[center_name]])
    delta_init <- as.numeric(context$init[[delta_name]])
    center_scale_ref <- if (abs(center_init) > 1e-12) {
      abs(center_init)
    } else {
      as.numeric(context$upper[[center_name]] - context$lower[[center_name]])
    }
    delta_scale_ref <- if (abs(delta_init) > 1e-12) {
      abs(delta_init)
    } else {
      as.numeric(context$upper[[delta_name]] - context$lower[[delta_name]])
    }
    config_rows[[config_index]] <- data.frame(
      family = family,
      warmup_label = warmup_label,
      invivo_seed = selected_row$invivo_seed[[1L]],
      invitro_seed = selected_row$invitro_seed[[1L]],
      parameter = parameters[[parameter_index]],
      transformation = transformation,
      center_name = center_name,
      delta_name = delta_name,
      center_init = center_init,
      delta_init = delta_init,
      center_lower = as.numeric(context$lower[[center_name]]),
      center_upper = as.numeric(context$upper[[center_name]]),
      delta_lower = as.numeric(context$lower[[delta_name]]),
      delta_upper = as.numeric(context$upper[[delta_name]]),
      center_sampling_sd = center_scale_ref * sigmaN,
      delta_sampling_sd = delta_scale_ref * sigmaN,
      joint_union_lower_transformed =
        as.numeric(meta$joint_union_lower_t[[1L]]),
      joint_union_upper_transformed =
        as.numeric(meta$joint_union_upper_t[[1L]]),
      joint_warmup_sigmaN = sigmaN,
      n_joint_seeds = length(joint_seeds),
      NP_used = NP_used,
      replay_cores = replay_cores,
      n_initial_population_values = n_rows,
      first_member_is_exact_warm_start = TRUE,
      sampling_definition = paste0(
        "center and delta drawn by joint_deoptim_initial_population from ",
        "warm-start-centered truncated normals; delta bounds are recomputed ",
        "from each sampled center; member 1 is the exact warm start"
      ),
      displayed_quantity = "log2(in vivo / in vitro)",
      context_quantity = "natural-scale in vivo and in vitro parameter values",
      ratio_formula = if (identical(transformation, "log10")) {
        "delta * log2(10)"
      } else {
        "log2((center + delta/2) / (center - delta/2))"
      },
      stringsAsFactors = FALSE
    )
    config_index <- config_index + 1L
  }

  provenance_inputs <- c(
    file.path(seed1_dir, "run_effective_args.tsv"),
    file.path(seed1_dir, "fit_summary.tsv"),
    file.path(seed1_dir, "joint_warmup_initial_values.tsv"),
    file.path(seed1_dir, "joint_soft_coupling_initial_values.tsv"),
    file.path(seed1_dir, "joint_soft_coupling_parameters_table_input.csv"),
    file.path(seed1_dir, "parameter_table_input_invivo.csv"),
    file.path(seed1_dir, "parameter_table_input_invitro.csv")
  )
  for (input_path in provenance_inputs) {
    provenance_rows[[provenance_index]] <- data.frame(
      family = family,
      warmup_label = warmup_label,
      role = basename(input_path),
      source_path = normalizePath(input_path, mustWork = TRUE),
      md5 = unname(tools::md5sum(input_path)),
      stringsAsFactors = FALSE
    )
    provenance_index <- provenance_index + 1L
  }
}

initial_population <- do.call(rbind, population_rows)
rownames(initial_population) <- NULL
initial_context_values <- do.call(rbind, context_value_rows)
rownames(initial_context_values) <- NULL
config <- do.call(rbind, config_rows)
rownames(config) <- NULL
provenance <- do.call(rbind, provenance_rows)
rownames(provenance) <- NULL

expected_rows <- 3L * 500L * 400L
if (nrow(initial_population) != expected_rows ||
    nrow(initial_context_values) != expected_rows ||
    nrow(config) != 3L * 14L ||
    any(table(initial_population$family) != 500L * 400L) ||
    any(table(initial_population$family, initial_population$joint_seed) != 400L) ||
    any(table(initial_population$family, initial_population$exact_warm_start)[, "FALSE"] !=
          199500L) ||
    any(table(initial_population$family, initial_population$exact_warm_start)[, "TRUE"] !=
          500L) ||
    any(!vapply(
      initial_context_values[, c(
        paste0("vivo__", parameters), paste0("vitro__", parameters)
      ), drop = FALSE],
      function(x) all(is.finite(x) & x > 0), logical(1)
    ))) {
  stop("Reconstructed DE initial-population coverage failed.")
}

population_path <- file.path(
  figure5_dir, "figure5f_de_initial_population_log2_ratios.rds"
)
saveRDS(initial_population, population_path, compress = "xz")
context_population_path <- file.path(
  figure5_dir, "figure5f_de_initial_population_context_values.rds"
)
saveRDS(initial_context_values, context_population_path, compress = "gzip")
write_tsv(
  config,
  file.path(figure5_dir, "figure5f_de_initial_population_config.tsv")
)
write_tsv(
  provenance,
  file.path(figure5_dir, "figure5f_de_initial_population_provenance.tsv")
)

release_checks <- data.frame(
  check = c(
    "selected_one_pair_per_family",
    "selected_pairs_match_approved_contract",
    "complete_joint_seed_sets",
    "NP_used",
    "joint_warmup_sigmaN",
    "initial_population_rows",
    "initial_context_value_rows",
    "initial_values_per_family_parameter",
    "exact_warm_start_members_per_family",
    "log2_ratios_finite",
    "context_values_positive_and_finite",
    "backend_snapshot_is_commit_83953a8"
  ),
  observed = c(
    nrow(selected),
    paste(observed_selection[families], collapse = ";"),
    "seed1--seed500",
    paste(sort(unique(config$NP_used)), collapse = ","),
    paste(sort(unique(config$joint_warmup_sigmaN)), collapse = ","),
    nrow(initial_population),
    nrow(initial_context_values),
    500L * 400L,
    paste(as.integer(table(
      initial_population$family[initial_population$exact_warm_start]
    )[families]), collapse = ","),
    all(vapply(
      initial_population[, parameters, drop = FALSE],
      function(x) all(is.finite(x)), logical(1)
    )),
    all(vapply(
      initial_context_values[, c(
        paste0("vivo__", parameters), paste0("vitro__", parameters)
      ), drop = FALSE],
      function(x) all(is.finite(x) & x > 0), logical(1)
    )),
    grepl("miningcloneid_83953a8", backend_script, fixed = TRUE)
  ),
  expected = c(
    "3",
    paste(expected_selection[families], collapse = ";"),
    "seed1--seed500",
    "400",
    "0.1216",
    as.character(expected_rows),
    as.character(expected_rows),
    as.character(500L * 400L),
    "500,500,500",
    "TRUE",
    "TRUE",
    "TRUE"
  ),
  stringsAsFactors = FALSE
)
release_checks$passed <-
  as.character(release_checks$observed) == release_checks$expected
write_tsv(
  release_checks,
  file.path(figure5_dir, "figure5f_de_initial_population_readiness.tsv")
)
if (any(!release_checks$passed)) {
  stop(
    "Figure 5F DE initial-population release checks failed: ",
    paste(release_checks$check[!release_checks$passed], collapse = ", ")
  )
}

cat("Reconstructed exact DE initial populations for Figure 5F.\n")
cat("Families:", paste(families, collapse = ", "), "\n")
cat("Joint seeds per family: 500\n")
cat("DE population members per seed: 400\n")
cat("Rows:", nrow(initial_population), "\n")
cat("Natural-scale context rows:", nrow(initial_context_values), "\n")
