#!/usr/bin/env Rscript

# Generate three family-conditioned generalized posteriors for Figure 5F from
# the exact joint-fit objective reconstructed from miningcloneid commit 83953a8.
#
# Target at learning temperature T:
#   pi_T(theta | data) proportional to
#     exp[-L_data(theta) / T] * exp[-L_regularization(theta)] * I(theta in B)
#
# where L_data is the weighted in-vivo plus in-vitro empirical loss,
# L_regularization is the configured in-vivo soft-prior penalty plus the Welsch
# soft-coupling penalty, and B is the exact bounded feasible joint parameter
# space intersected with exactly one frozen family basin V_f.
# C01, C02, and C03 therefore define distinct conditional targets. T=1
# reproduces the saved optimization objective inside the selected basin;
# T=0.5 and T=2 are sensitivity targets. Hotter replicas at T=4 and T=8
# support within-family replica exchange but are not scientific outputs.
# Two independent replica-exchange ladders are run for each family. Their
# convergence and stopping decisions are made per family, never across all
# three targets.

resolve_own_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files)) return(frame_files[[length(frame_files)]])
  stop("Unable to resolve the Figure 5F generator path.")
}

generator_path <- resolve_own_file()
generator_dir <- dirname(generator_path)
iteration_root <- normalizePath(file.path(generator_dir, "..", ".."), mustWork = TRUE)
snapshot_root <- file.path(iteration_root, "tmp", "miningcloneid_83953a8")
joint_backend_path <- file.path(
  snapshot_root,
  "oxygen", "code", "O2_supply_demand_MAP", "util",
  "o2_supply_demand_map_fit_joint_backend.R"
)

# Historical backends resolve SCRIPT_DIR from commandArgs before source frames.
# Relaunch through Rscript -e so the backend itself is the only active source
# frame when it is loaded. This keeps the tracked historical source untouched
# and allows this generator to remain a directly executable standalone script.
if (!identical(Sys.getenv("FIGURE5F_BACKEND_PRELOADED"), "1")) {
  if (!file.exists(joint_backend_path)) {
    stop("Missing exact historical joint backend: ", joint_backend_path)
  }
  expr <- paste0(
    "joint_env <- new.env(parent=globalenv()); ",
    "sys.source(", encodeString(joint_backend_path, quote = '"'),
    ", envir=joint_env, chdir=TRUE); ",
    "Sys.setenv(FIGURE5F_BACKEND_PRELOADED='1'); ",
    "source(", encodeString(generator_path, quote = '"'), ")"
  )
  trailing <- commandArgs(trailingOnly = TRUE)
  child_args <- c("-e", shQuote(expr))
  if (length(trailing)) {
    child_args <- c(child_args, "--args", vapply(trailing, shQuote, character(1)))
  }
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = child_args,
    stdout = "",
    stderr = ""
  )
  quit(save = "no", status = status, runLast = FALSE)
}

if (!exists("joint_env", envir = .GlobalEnv, inherits = FALSE) ||
    !is.environment(get("joint_env", envir = .GlobalEnv, inherits = FALSE))) {
  stop("Historical joint backend preload failed.")
}
joint_env <- get("joint_env", envir = .GlobalEnv, inherits = FALSE)

parse_cli <- function(x) {
  out <- list()
  for (arg in x) {
    if (!startsWith(arg, "--")) next
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- bits[[1L]]
    value <- if (length(bits) > 1L) paste(bits[-1L], collapse = "=") else "TRUE"
    out[[key]] <- value
  }
  out
}

as_integer <- function(x, default) {
  value <- suppressWarnings(as.integer(x))
  if (!length(value) || is.na(value) || value < 1L) as.integer(default) else value
}

as_number <- function(x, default) {
  value <- suppressWarnings(as.numeric(x))
  if (!length(value) || !is.finite(value)) as.numeric(default) else value
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(isTRUE(default))
  tolower(trimws(as.character(x[[1L]]))) %in% c("1", "true", "t", "yes", "y")
}

parse_numeric_vector <- function(x, default) {
  if (is.null(x) || !length(x)) return(as.numeric(default))
  values <- suppressWarnings(as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1L]]))
  if (!length(values) || any(!is.finite(values))) {
    stop("Invalid numeric vector: ", as.character(x))
  }
  values
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
mode <- tolower(if (is.null(cli$mode)) "sample" else as.character(cli$mode))
if (!mode %in% c("benchmark", "pilot", "sample")) {
  stop("--mode must be benchmark, pilot, or sample.")
}
ladder_family <- toupper(trimws(
  if (is.null(cli$ladder_family)) "" else as.character(cli$ladder_family)
))
ladder_replicate <- if (is.null(cli$ladder_replicate)) {
  NA_integer_
} else {
  as_integer(cli$ladder_replicate, NA_integer_)
}
aggregate_only <- as_flag(cli$aggregate_only, FALSE)
single_ladder_mode <- nzchar(ladder_family)
if (single_ladder_mode && aggregate_only) {
  stop("--ladder_family and --aggregate_only cannot be used together.")
}
if (single_ladder_mode && !ladder_family %in% c("C01", "C02", "C03")) {
  stop("--ladder_family must be C01, C02, or C03.")
}

default_iter <- if (mode == "pilot") 600L else 6000L
default_burn <- if (mode == "pilot") 300L else 3000L
# Retain every post-warm-up state. Thinning does not improve Markov-chain ESS
# per target evaluation and would discard information from this expensive model.
default_thin <- 1L
cfg <- list(
  mode = mode,
  n_iter = as_integer(cli$n_iter, default_iter),
  burnin = as_integer(cli$burnin, default_burn),
  thin = as_integer(cli$thin, default_thin),
  ladders_per_family = as_integer(cli$ladders_per_family, 2L),
  temperatures = parse_numeric_vector(cli$temperatures, c(0.5, 1, 2, 4, 8)),
  display_temperatures = parse_numeric_vector(cli$display_temperatures, c(0.5, 1, 2)),
  cores = as_integer(
    cli$cores,
    max(1L, min(6L, parallel::detectCores(logical = FALSE) - 1L))
  ),
  temperature_cores_per_ladder = as_integer(
    cli$temperature_cores_per_ladder,
    1L
  ),
  checkpoint_every = as_integer(cli$checkpoint_every, 250L),
  base_seed = as_integer(cli$base_seed, 20260813L),
  benchmark_evals = as_integer(cli$benchmark_evals, 20L),
  force = as_flag(cli$force, FALSE),
  resume = as_flag(cli$resume, TRUE)
)
if (cfg$burnin >= cfg$n_iter) stop("burnin must be smaller than n_iter.")
if (!all(cfg$display_temperatures %in% cfg$temperatures)) {
  stop("Every display temperature must occur in the replica-exchange ladder.")
}
if (any(cfg$temperatures <= 0) || is.unsorted(cfg$temperatures, strictly = TRUE)) {
  stop("temperatures must be strictly increasing and positive.")
}
if (!1 %in% cfg$temperatures) stop("The temperature ladder must include T=1.")
if (cfg$cores < 1L || cfg$temperature_cores_per_ladder < 1L) {
  stop("cores and temperature_cores_per_ladder must both be positive.")
}
if (cfg$temperature_cores_per_ladder > length(cfg$temperatures)) {
  stop("temperature_cores_per_ladder cannot exceed the temperature count.")
}
if (single_ladder_mode &&
    (is.na(ladder_replicate) || ladder_replicate > cfg$ladders_per_family)) {
  stop("--ladder_replicate must identify one configured ladder replicate.")
}

figure5_dir <- file.path(iteration_root, "data", "Figures", "Figure5")
legacy_output_dir <- file.path(figure5_dir, "generalized_posterior")
output_dir <- file.path(figure5_dir, "generalized_posterior_family_conditioned")
checkpoint_dir <- file.path(
  iteration_root, "tmp", "figure5f_generalized_posterior_family_conditioned"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

selection_path <- file.path(figure5_dir, "figure5f_selected_pair_inputs.tsv")
prior_draw_path <- file.path(
  output_dir, "figure5f_family_conditioned_prior_draws.tsv"
)
family_contract_path <- file.path(output_dir, "figure5f_family_basin_contract.rds")
map_audit_path <- file.path(figure5_dir, "figure5f_map_replay_audit.tsv")
proposal_covariance_path <- file.path(
  legacy_output_dir, "figure5f_pilot2_proposal_covariance.rds"
)
proposal_calibration_script_path <- file.path(
  iteration_root, "Code", "Figures", "calibrate_Figure5F_pilot_proposal.R"
)
result_root <- Sys.getenv(
  "FIGURE5F_RESULT_ROOT",
  unset = paste0(
    "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/",
    "fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540"
  )
)
result_root <- normalizePath(result_root, mustWork = TRUE)

required_files <- c(
  selection_path,
  prior_draw_path,
  family_contract_path,
  map_audit_path,
  proposal_covariance_path,
  proposal_calibration_script_path,
  file.path(snapshot_root, "oxygen", "config", "O2_supply_demand.yaml"),
  file.path(snapshot_root, "data", "InVivoData_Gemcitabine", "dt_Gem_VT_20260209_v5.xlsx"),
  file.path(snapshot_root, "data", "InVivoData_Gemcitabine", "all_ploidy.csv"),
  file.path(
    snapshot_root, "data", "InVivoData_Gemcitabine",
    "histology_to_dt_Gem_VT_20260209_v5_mapping.csv"
  ),
  file.path(snapshot_root, "oxygen", "ploidyOxygen", "data", "fit_objects", "fit_data.Rds"),
  file.path(snapshot_root, "oxygen", "ploidyOxygen", "data", "fit_objects", "jobs_2N.Rds"),
  file.path(snapshot_root, "oxygen", "ploidyOxygen", "data", "fit_objects", "jobs_4N.Rds"),
  file.path(snapshot_root, "oxygen", "data", "g0g1_ploidy_density_grid.csv")
)
missing <- required_files[!file.exists(required_files)]
if (length(missing)) stop("Missing Figure 5F sampling input(s):\n", paste(missing, collapse = "\n"))

map_audit <- utils::read.delim(map_audit_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!nrow(map_audit) || !all(map_audit$passed %in% TRUE) ||
    !all(map_audit$replay_code_commit == "83953a874401e42cd176432786f889a896adc959")) {
  stop("Figure 5F MAP replay audit is absent or failed; sampling is blocked.")
}

selection <- utils::read.delim(selection_path, check.names = FALSE, stringsAsFactors = FALSE)
selection <- selection[selection$selected_for_figure5f %in% TRUE, , drop = FALSE]
selection <- selection[order(match(selection$family, c("C01", "C02", "C03"))), , drop = FALSE]
if (nrow(selection) != 3L || !identical(selection$family, c("C01", "C02", "C03"))) {
  stop("Expected exactly one selected MAP input for C01, C02, and C03.")
}

read_effective_args <- function(path) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  tab <- tab[tab$source == "fit_command", c("key", "value"), drop = FALSE]
  if (!nrow(tab) || anyDuplicated(tab$key)) stop("Malformed effective-args table: ", path)
  out <- as.list(tab$value)
  names(out) <- tab$key
  out
}

selected_run_dir <- function(row) {
  file.path(
    result_root,
    paste0("fit_joint_", row$warmup_label[[1L]]),
    as.character(row$selected_seed[[1L]])
  )
}

context_source_dir <- selected_run_dir(selection[1L, , drop = FALSE])
context_inputs <- file.path(
  context_source_dir,
  c(
    "run_effective_args.tsv",
    "parameter_table_input_invivo.csv",
    "parameter_table_input_invitro.csv",
    "joint_soft_coupling_parameters_table_input.csv"
  )
)
if (any(!file.exists(context_inputs))) {
  stop("Missing exact context input(s) under: ", context_source_dir)
}

argv <- read_effective_args(file.path(context_source_dir, "run_effective_args.tsv"))
argv$config <- file.path(snapshot_root, "oxygen", "config", "O2_supply_demand.yaml")
argv$data_dir <- file.path(snapshot_root, "data", "InVivoData_Gemcitabine")
argv$parameter_table <- file.path(context_source_dir, "parameter_table_input_invivo.csv")
argv$invitro_parameter_table <- file.path(context_source_dir, "parameter_table_input_invitro.csv")
argv$fit_objects_dir <- file.path(snapshot_root, "oxygen", "ploidyOxygen", "data", "fit_objects")
argv$flow_density_path <- file.path(snapshot_root, "oxygen", "data", "g0g1_ploidy_density_grid.csv")
argv$necrosis_mapping_csv <- file.path(
  snapshot_root, "data", "InVivoData_Gemcitabine",
  "histology_to_dt_Gem_VT_20260209_v5_mapping.csv"
)
argv$joint_soft_coupling_parameters_table <- file.path(
  context_source_dir, "joint_soft_coupling_parameters_table_input.csv"
)
argv$joint_warmup_enable <- "FALSE"
argv$n_cores <- "1"
argv$joint_n_cores <- "1"
argv$predict_n_cores <- "1"
argv$deoptim_parallel <- "FALSE"
context_label <- if (single_ladder_mode) {
  paste0("context_", ladder_family, "_R", ladder_replicate)
} else if (aggregate_only) {
  "context_aggregate"
} else {
  "context_only"
}
argv$out_dir <- file.path(checkpoint_dir, context_label)

message("Building the exact T=1 joint target context...")
ctx <- joint_env$build_joint_context(argv)
soft_meta <- ctx$joint_soft_coupling$metadata
if (!is.data.frame(soft_meta) || nrow(soft_meta) != 14L) {
  stop("Expected exactly 14 active soft-coupled parameters.")
}

# Sampling coordinates use the context-specific transformed values directly.
# The center/delta -> in-vivo/in-vitro map has constant unit Jacobian, and this
# coordinate system enforces the dynamic feasibility constraints by design.
center_names <- as.character(soft_meta$center_name)
delta_names <- as.character(soft_meta$delta_name)
ordinary_names <- setdiff(names(ctx$init), c(center_names, delta_names))
coord_rows <- list()
for (name in ordinary_names) {
  coord_rows[[length(coord_rows) + 1L]] <- data.frame(
    coordinate = name,
    role = "ordinary",
    parameter = name,
    original_name = name,
    lower = as.numeric(ctx$lower[[name]]),
    upper = as.numeric(ctx$upper[[name]]),
    transform = NA_character_,
    stringsAsFactors = FALSE
  )
}
for (i in seq_len(nrow(soft_meta))) {
  for (context_name in c("vivo", "vitro")) {
    coord_rows[[length(coord_rows) + 1L]] <- data.frame(
      coordinate = paste0(context_name, "__", soft_meta$center_name[[i]]),
      role = context_name,
      parameter = as.character(soft_meta$parameter[[i]]),
      original_name = as.character(soft_meta$center_name[[i]]),
      lower = as.numeric(soft_meta$joint_union_lower_t[[i]]),
      upper = as.numeric(soft_meta$joint_union_upper_t[[i]]),
      transform = as.character(soft_meta$transform[[i]]),
      stringsAsFactors = FALSE
    )
  }
}
coord_meta <- do.call(rbind, coord_rows)
coord_meta$width <- coord_meta$upper - coord_meta$lower
if (nrow(coord_meta) != length(ctx$init) || anyDuplicated(coord_meta$coordinate) ||
    any(!is.finite(coord_meta$lower)) || any(!is.finite(coord_meta$upper)) ||
    any(coord_meta$width <= 0)) {
  stop("Invalid or non-bijective Figure 5F sampling-coordinate map.")
}

# A failed, separately archived pilot supplies fixed covariance geometry for a
# symmetric multivariate random-walk proposal. The pilot contributes no draws,
# weights, centers, or likelihood terms to the production target. Freezing this
# matrix before production is equivalent to selecting a proposal scale during
# a preliminary tuning run and does not adapt the production target.
proposal_preconditioner <- readRDS(proposal_covariance_path)
proposal_coordinates <- as.character(proposal_preconditioner$coordinate_names)
if (!isFALSE(proposal_preconditioner$scientific_draws_included) ||
    !identical(proposal_coordinates, coord_meta$coordinate) ||
    !is.list(proposal_preconditioner$covariance) ||
    !length(proposal_preconditioner$covariance)) {
  coordinate_evidence <- if (setequal(proposal_coordinates, coord_meta$coordinate)) {
    paste0(
      "same names but different order; first expected/observed mismatch: ",
      which(proposal_coordinates != coord_meta$coordinate)[[1L]], " = ",
      proposal_coordinates[which(proposal_coordinates != coord_meta$coordinate)[[1L]]],
      " / ", coord_meta$coordinate[which(proposal_coordinates != coord_meta$coordinate)[[1L]]]
    )
  } else {
    paste0(
      "missing=", paste(setdiff(coord_meta$coordinate, proposal_coordinates), collapse = ","),
      "; extra=", paste(setdiff(proposal_coordinates, coord_meta$coordinate), collapse = ",")
    )
  }
  stop(
    "Fixed pilot proposal preconditioner fails its provenance or coordinate contract: ",
    coordinate_evidence
  )
}
proposal_covariance_keys <- suppressWarnings(as.numeric(names(proposal_preconditioner$covariance)))
if (any(!is.finite(proposal_covariance_keys)) || any(proposal_covariance_keys <= 0)) {
  stop("Pilot proposal covariance temperatures are invalid.")
}
proposal_cholesky <- lapply(proposal_preconditioner$covariance, function(matrix_value) {
  matrix_value <- as.matrix(matrix_value)
  if (!identical(rownames(matrix_value), coord_meta$coordinate) ||
      !identical(colnames(matrix_value), coord_meta$coordinate) ||
      any(!is.finite(matrix_value))) {
    stop("Pilot proposal covariance matrix is malformed.")
  }
  chol(matrix_value)
})
proposal_cholesky_for_temperature <- function(temperature) {
  key_i <- which.min(abs(log(temperature) - log(proposal_covariance_keys)))
  proposal_cholesky[[key_i]]
}

sampling_to_optimizer <- function(unit_state) {
  if (length(unit_state) != nrow(coord_meta) || any(!is.finite(unit_state)) ||
      any(unit_state < -1e-12) || any(unit_state > 1 + 1e-12)) {
    stop("Invalid unit sampling state.")
  }
  actual <- coord_meta$lower + unit_state * coord_meta$width
  names(actual) <- coord_meta$coordinate
  par_t <- setNames(as.numeric(ctx$init), names(ctx$init))
  ordinary <- coord_meta$role == "ordinary"
  par_t[coord_meta$original_name[ordinary]] <- actual[coord_meta$coordinate[ordinary]]
  for (i in seq_len(nrow(soft_meta))) {
    center_name <- as.character(soft_meta$center_name[[i]])
    vivo <- actual[[paste0("vivo__", center_name)]]
    vitro <- actual[[paste0("vitro__", center_name)]]
    par_t[[center_name]] <- (vivo + vitro) / 2
    par_t[[as.character(soft_meta$delta_name[[i]])]] <- vivo - vitro
  }
  par_t[names(ctx$init)]
}

optimizer_to_sampling <- function(par_t) {
  par_t <- par_t[names(ctx$init)]
  actual <- setNames(numeric(nrow(coord_meta)), coord_meta$coordinate)
  ordinary <- coord_meta$role == "ordinary"
  actual[coord_meta$coordinate[ordinary]] <- par_t[coord_meta$original_name[ordinary]]
  for (i in seq_len(nrow(soft_meta))) {
    center_name <- as.character(soft_meta$center_name[[i]])
    delta_name <- as.character(soft_meta$delta_name[[i]])
    center <- as.numeric(par_t[[center_name]])
    delta <- as.numeric(par_t[[delta_name]])
    actual[[paste0("vivo__", center_name)]] <- center + delta / 2
    actual[[paste0("vitro__", center_name)]] <- center - delta / 2
  }
  unit <- (actual[coord_meta$coordinate] - coord_meta$lower) / coord_meta$width
  unit[abs(unit) < 1e-12] <- 0
  unit[abs(unit - 1) < 1e-12] <- 1
  if (any(unit < -1e-9) || any(unit > 1 + 1e-9)) {
    stop("Saved MAP lies outside the reconstructed sampling box.")
  }
  pmin(pmax(unit, 0), 1)
}

evaluate_loss <- function(unit_state) {
  par_t <- sampling_to_optimizer(unit_state)
  comp <- tryCatch(
    joint_env$joint_objective_components(par_t, ctx),
    error = function(e) NULL
  )
  if (is.null(comp)) {
    return(c(data_loss = 1e9, regularization_loss = 0, objective_t1 = 1e9))
  }
  data_loss <- ctx$joint_weight_invivo * as.numeric(comp$invivo$L_data) +
    ctx$joint_weight_invitro * as.numeric(comp$invitro$objective)
  constraint_loss <- as.numeric(
    comp$constraint_metrics$joint_constraint_penalty_total
  )
  if (!length(constraint_loss) || !is.finite(constraint_loss)) {
    constraint_loss <- 0
  }
  regularization_loss <- ctx$joint_weight_invivo * as.numeric(comp$invivo$L_prior) +
    as.numeric(comp$objective_soft_coupling) + constraint_loss
  objective_t1 <- data_loss + regularization_loss
  if (!all(is.finite(c(data_loss, regularization_loss, objective_t1))) ||
      objective_t1 >= 1e8) {
    return(c(data_loss = 1e9, regularization_loss = 0, objective_t1 = 1e9))
  }
  full_objective <- as.numeric(comp$objective)
  if (!length(full_objective) || !is.finite(full_objective) ||
      abs(objective_t1 - full_objective) > 1e-8) {
    stop(
      "Generalized-posterior T=1 decomposition no longer reproduces the ",
      "historical joint objective."
    )
  }
  c(
    data_loss = data_loss,
    regularization_loss = regularization_loss,
    objective_t1 = objective_t1
  )
}

map_states <- list()
map_loss_rows <- list()
for (i in seq_len(nrow(selection))) {
  row <- selection[i, , drop = FALSE]
  run_dir <- selected_run_dir(row)
  par_path <- file.path(run_dir, "best_params_transformed.tsv")
  if (!file.exists(par_path)) stop("Missing selected MAP parameter table: ", par_path)
  tab <- utils::read.delim(par_path, check.names = FALSE, stringsAsFactors = FALSE)
  par_t <- setNames(as.numeric(tab$transformed_value), tab$transformed_parameter)
  if (!setequal(names(par_t), names(ctx$init))) {
    stop("Selected MAP parameter names differ for ", row$family[[1L]])
  }
  unit <- optimizer_to_sampling(par_t)
  loss <- evaluate_loss(unit)
  saved <- map_audit$saved_value[
    map_audit$family == row$family[[1L]] & map_audit$component == "objective_joint"
  ]
  if (length(saved) != 1L || abs(loss[["objective_t1"]] - saved) > 1e-8) {
    stop("Base objective context does not replay ", row$family[[1L]], " MAP.")
  }
  map_states[[row$family[[1L]]]] <- unit
  map_loss_rows[[i]] <- data.frame(
    family = row$family[[1L]],
    warmup_label = row$warmup_label[[1L]],
    selected_seed = row$selected_seed[[1L]],
    data_loss = loss[["data_loss"]],
    regularization_loss = loss[["regularization_loss"]],
    objective_t1 = loss[["objective_t1"]],
    saved_objective = saved,
    absolute_difference = abs(loss[["objective_t1"]] - saved),
    stringsAsFactors = FALSE
  )
}
map_loss_table <- do.call(rbind, map_loss_rows)

family_contract <- readRDS(family_contract_path)
expected_target_version <- "family_conditioned_selected_map_voronoi_v1"
expected_family_features <- c(
  "log10_o2_S0", "log10_kappa_O", "log10_eta_o2", "log10_k_clear",
  coord_meta$coordinate[coord_meta$role == "vivo"]
)
if (!is.list(family_contract) ||
    !identical(family_contract$target_version, expected_target_version) ||
    !identical(family_contract$families, c("C01", "C02", "C03")) ||
    !identical(family_contract$tie_order, c("C01", "C02", "C03")) ||
    !identical(family_contract$coordinate_names, coord_meta$coordinate) ||
    !identical(family_contract$feature_names, expected_family_features) ||
    !identical(family_contract$selection_md5, unname(tools::md5sum(selection_path))) ||
    !is.matrix(family_contract$centers) ||
    !identical(rownames(family_contract$centers), c("C01", "C02", "C03")) ||
    !identical(colnames(family_contract$centers), expected_family_features)) {
  stop("Frozen Figure 5F family-basin contract is malformed or stale.")
}
expected_centers <- t(vapply(
  c("C01", "C02", "C03"),
  function(family) map_states[[family]][expected_family_features],
  numeric(length(expected_family_features))
))
if (max(abs(expected_centers - family_contract$centers)) > 1e-12) {
  stop("Frozen family-basin centers no longer match the selected MAP states.")
}

classify_family_state <- function(unit_state) {
  x <- as.numeric(unit_state[family_contract$feature_names])
  if (length(x) != length(family_contract$feature_names) || any(!is.finite(x))) {
    stop("Cannot classify a malformed sampling state.")
  }
  distances <- vapply(
    family_contract$families,
    function(family) sum((x - family_contract$centers[family, ])^2),
    numeric(1)
  )
  family_contract$tie_order[[which.min(distances)]]
}

if (!identical(
  unname(vapply(map_states, classify_family_state, character(1))),
  family_contract$families
)) {
  stop("At least one selected MAP does not belong to its own frozen basin.")
}

write_tsv <- function(x, path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  allowed <- normalizePath(figure5_dir, mustWork = TRUE)
  if (!startsWith(normalized, paste0(allowed, .Platform$file.sep))) {
    stop("Refusing to write outside iteration2 Figure5 data: ", normalized)
  }
  dir.create(dirname(normalized), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, normalized, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
  invisible(normalized)
}

if (!single_ladder_mode) {
  write_tsv(map_loss_table, file.path(output_dir, "figure5f_selected_map_target_replay.tsv"))
  write_tsv(coord_meta, file.path(output_dir, "figure5f_sampling_coordinate_contract.tsv"))
}

reflect_unit <- function(x) {
  y <- x %% 2
  y[y > 1] <- 2 - y[y > 1]
  pmin(pmax(y, 0), 1)
}

if (mode == "benchmark") {
  set.seed(cfg$base_seed)
  state <- map_states[["C01"]]
  proposals <- lapply(seq_len(cfg$benchmark_evals), function(i) {
    dims <- sample.int(length(state), 5L)
    proposal <- state
    proposal[dims] <- reflect_unit(
      proposal[dims] + stats::rnorm(5L, 0, 0.002)
    )
    proposal
  })
  serial_elapsed <- system.time({
    serial_losses <- lapply(proposals, evaluate_loss)
  })[["elapsed"]]
  parallel_elapsed <- serial_elapsed
  parallel_losses <- serial_losses
  benchmark_workers <- min(
    cfg$temperature_cores_per_ladder,
    length(cfg$temperatures),
    length(proposals)
  )
  if (.Platform$OS.type == "unix" && benchmark_workers > 1L) {
    benchmark_cluster <- parallel::makeForkCluster(
      benchmark_workers,
      port = 20000L + (as.integer(Sys.getpid()) %% 30000L)
    )
    parallel_elapsed <- system.time({
      parallel_losses <- parallel::parLapply(
        benchmark_cluster,
        proposals,
        function(proposal) evaluate_loss(proposal)
      )
    })[["elapsed"]]
    parallel::stopCluster(benchmark_cluster)
  }
  serial_matrix <- do.call(rbind, serial_losses)
  parallel_matrix <- do.call(rbind, parallel_losses)
  parity_max_abs_difference <- max(abs(serial_matrix - parallel_matrix))
  if (!is.finite(parity_max_abs_difference) ||
      parity_max_abs_difference > 1e-12) {
    stop(
      "Temperature-parallel objective parity failed: max abs difference = ",
      format(parity_max_abs_difference, scientific = TRUE)
    )
  }
  bench <- data.frame(
    evaluation = seq_len(cfg$benchmark_evals),
    data_loss = serial_matrix[, "data_loss"],
    regularization_loss = serial_matrix[, "regularization_loss"],
    objective_t1 = serial_matrix[, "objective_t1"],
    serial_total_seconds = serial_elapsed,
    parallel_total_seconds = parallel_elapsed,
    parallel_worker_count = benchmark_workers,
    parallel_speedup = serial_elapsed / parallel_elapsed,
    serial_parallel_max_abs_difference = parity_max_abs_difference,
    stringsAsFactors = FALSE
  )
  bench$mean_serial_seconds_per_evaluation <-
    serial_elapsed / cfg$benchmark_evals
  parallel_batches <- ceiling(cfg$benchmark_evals / benchmark_workers)
  bench$mean_parallel_seconds_per_temperature_batch <-
    parallel_elapsed / parallel_batches
  bench$projected_evaluations_per_ladder <- cfg$n_iter * length(cfg$temperatures)
  bench$projected_parallel_wall_minutes <-
    bench$mean_parallel_seconds_per_temperature_batch * cfg$n_iter / 60
  write_tsv(bench, file.path(output_dir, "figure5f_objective_benchmark.tsv"))
  print(bench)
  message(
    "Temperature-parallel benchmark used ", benchmark_workers,
    " workers; serial/parallel speedup = ",
    signif(serial_elapsed / parallel_elapsed, 4),
    "; parity max abs difference = ",
    format(parity_max_abs_difference, scientific = TRUE),
    "; projected wall minutes for one full ladder = ",
    signif(bench$projected_parallel_wall_minutes[[1L]], 4)
  )
  quit(save = "no", status = 0L, runLast = FALSE)
}

# Every proposal is a symmetric random walk. The hard family indicator is part
# of the target, so proposals crossing a frozen basin boundary are rejected
# before the expensive model objective is evaluated. Cross-family MAP jumps are
# invalid for these conditional targets and are deliberately absent.
kernel_names <- c("single", "block8", "pilot_cov", "full")
kernel_probabilities <- c(
  single = 0.40, block8 = 0.25, pilot_cov = 0.30, full = 0.05
)
kernel_targets <- c(
  single = 0.44, block8 = 0.30, pilot_cov = 0.234, full = 0.234
)
initial_scales <- c(
  single = 0.030, block8 = 0.012, pilot_cov = 0.40, full = 0.004
)
minimum_scales <- c(
  single = 1e-5, block8 = 1e-5, pilot_cov = 0.02, full = 1e-5
)
maximum_scales <- c(
  single = 0.35, block8 = 0.35, pilot_cov = 2.5, full = 0.35
)
n_save <- floor((cfg$n_iter - cfg$burnin) / cfg$thin)
if (n_save < 100L) stop("Sampling configuration retains fewer than 100 draws per chain.")

generator_md5_at_start <- unname(tools::md5sum(generator_path))
selection_md5_at_start <- unname(tools::md5sum(selection_path))
map_audit_md5_at_start <- unname(tools::md5sum(map_audit_path))
prior_draw_md5_at_start <- unname(tools::md5sum(prior_draw_path))
family_contract_md5_at_start <- unname(tools::md5sum(family_contract_path))
proposal_covariance_md5_at_start <- unname(tools::md5sum(proposal_covariance_path))
proposal_calibration_script_md5_at_start <- unname(
  tools::md5sum(proposal_calibration_script_path)
)
config_signature <- paste(
  generator_md5_at_start,
  selection_md5_at_start,
  map_audit_md5_at_start,
  prior_draw_md5_at_start,
  family_contract_md5_at_start,
  proposal_covariance_md5_at_start,
  proposal_calibration_script_md5_at_start,
  paste(cfg$temperatures, collapse = ","),
  # n_iter is intentionally excluded: a frozen production chain may be
  # extended from its final atomic checkpoint without changing its target,
  # warm-up, proposal, thinning, or random-number stream.
  cfg$burnin, cfg$thin, cfg$ladders_per_family,
  cfg$temperature_cores_per_ladder,
  cfg$base_seed, sep = "|"
)

expand_retained_to <- function(retained, target_rows) {
  lapply(retained, function(item) {
    current_rows <- nrow(item$unit)
    if (current_rows > target_rows || nrow(item$loss) != current_rows) {
      stop("Checkpoint retention shape is incompatible with the requested extension.")
    }
    if (current_rows < target_rows) {
      add_n <- target_rows - current_rows
      item$unit <- rbind(
        item$unit,
        matrix(
          NA_real_, nrow = add_n, ncol = ncol(item$unit),
          dimnames = list(NULL, colnames(item$unit))
        )
      )
      item$loss <- rbind(
        item$loss,
        matrix(
          NA_real_, nrow = add_n, ncol = ncol(item$loss),
          dimnames = list(NULL, colnames(item$loss))
        )
      )
    }
    item
  })
}

run_ladder <- function(family, replicate_id) {
  chain_id <- paste0(family, "_R", replicate_id)
  chain_seed <- cfg$base_seed +
    match(family, c("C01", "C02", "C03")) * 1000L + replicate_id
  checkpoint_path <- file.path(checkpoint_dir, paste0("ladder_", chain_id, ".rds"))
  result_path <- file.path(checkpoint_dir, paste0("ladder_", chain_id, "_complete.rds"))
  if (file.exists(result_path) && !cfg$force) {
    existing <- readRDS(result_path)
    if (identical(existing$config_signature, config_signature) &&
        !is.null(existing$completed_iteration) &&
        existing$completed_iteration >= cfg$n_iter) {
      return(existing)
    }
  }

  temperature_cluster <- NULL
  temperature_worker_count <- min(
    cfg$temperature_cores_per_ladder,
    length(cfg$temperatures)
  )
  if (.Platform$OS.type == "unix" && temperature_worker_count > 1L) {
    # Create the nested workers once per ladder. Re-forking five objective
    # workers at every MCMC iteration would overwhelm the actual model work.
    # Workers only evaluate deterministic proposals; all RNG, accept/reject,
    # adaptation, swaps, and checkpoint state remain in the ladder parent.
    # R_PARALLEL_PORT can force every nested cluster to request the same port.
    # Ladder parents have distinct process IDs, so derive an explicit port from
    # the PID and keep it separate from the scientific MCMC random stream.
    cluster_port <- 20000L + (as.integer(Sys.getpid()) %% 30000L)
    temperature_cluster <- parallel::makeForkCluster(
      temperature_worker_count,
      port = cluster_port
    )
    on.exit(parallel::stopCluster(temperature_cluster), add = TRUE)
  }

  set.seed(chain_seed)
  n_temp <- length(cfg$temperatures)
  n_dim <- nrow(coord_meta)
  states <- matrix(rep(map_states[[family]], each = n_temp), nrow = n_temp, byrow = FALSE)
  colnames(states) <- coord_meta$coordinate
  for (temperature_index in seq_len(n_temp)) {
    jitter_sd <- if (replicate_id == 1L && cfg$temperatures[[temperature_index]] == 1) {
      0
    } else {
      0.0015 * sqrt(cfg$temperatures[[temperature_index]]) * (1 + 0.35 * (replicate_id - 1L))
    }
    if (jitter_sd > 0) {
      candidate <- states[temperature_index, ]
      accepted_jitter <- FALSE
      for (jitter_attempt in seq_len(100L)) {
        candidate <- reflect_unit(
          map_states[[family]] + stats::rnorm(n_dim, 0, jitter_sd)
        )
        names(candidate) <- coord_meta$coordinate
        if (identical(classify_family_state(candidate), family)) {
          accepted_jitter <- TRUE
          break
        }
      }
      if (accepted_jitter) states[temperature_index, ] <- candidate
    }
  }
  if (any(vapply(
    seq_len(n_temp),
    function(i) !identical(classify_family_state(states[i, ]), family),
    logical(1)
  ))) {
    stop("Initial state escaped the frozen family basin for ", chain_id)
  }
  losses <- t(vapply(seq_len(n_temp), function(i) evaluate_loss(states[i, ]), numeric(3)))
  colnames(losses) <- c("data_loss", "regularization_loss", "objective_t1")
  scales <- matrix(
    rep(initial_scales, each = n_temp),
    nrow = n_temp,
    dimnames = list(as.character(cfg$temperatures), kernel_names)
  )
  local_attempts <- matrix(0L, nrow = n_temp, ncol = length(kernel_names), dimnames = dimnames(scales))
  local_accepts <- local_attempts
  basin_rejections <- integer(n_temp)
  swap_attempts <- integer(n_temp - 1L)
  swap_accepts <- integer(n_temp - 1L)
  retained <- lapply(cfg$display_temperatures, function(x) {
    list(
      unit = matrix(NA_real_, nrow = n_save, ncol = n_dim, dimnames = list(NULL, coord_meta$coordinate)),
      loss = matrix(
        NA_real_, nrow = n_save, ncol = 3L,
        dimnames = list(NULL, c("data_loss", "regularization_loss", "target_loss"))
      )
    )
  })
  names(retained) <- as.character(cfg$display_temperatures)
  start_iter <- 1L
  save_index <- 0L

  if (cfg$resume && file.exists(checkpoint_path) && !cfg$force) {
    checkpoint <- readRDS(checkpoint_path)
    if (identical(checkpoint$config_signature, config_signature) &&
        checkpoint$iteration < cfg$n_iter) {
      states <- checkpoint$states
      losses <- checkpoint$losses
      scales <- checkpoint$scales
      local_attempts <- checkpoint$local_attempts
      local_accepts <- checkpoint$local_accepts
      basin_rejections <- checkpoint$basin_rejections
      if (is.null(basin_rejections)) basin_rejections <- integer(n_temp)
      swap_attempts <- checkpoint$swap_attempts
      swap_accepts <- checkpoint$swap_accepts
      retained <- expand_retained_to(checkpoint$retained, n_save)
      save_index <- checkpoint$save_index
      expected_checkpoint_saves <- max(
        0L,
        floor((checkpoint$iteration - cfg$burnin) / cfg$thin)
      )
      if (save_index != expected_checkpoint_saves || save_index > n_save) {
        stop("Checkpoint retained-draw index is inconsistent for ", chain_id)
      }
      start_iter <- checkpoint$iteration + 1L
      assign(".Random.seed", checkpoint$random_seed, envir = .GlobalEnv)
      message("Resuming ", chain_id, " at iteration ", start_iter)
    }
  }

  checkpoint_state <- function(iteration) {
    object <- list(
      config_signature = config_signature,
      target_version = family_contract$target_version,
      family = family,
      replicate = replicate_id,
      chain = chain_id,
      iteration = iteration,
      target_iteration = cfg$n_iter,
      states = states,
      losses = losses,
      scales = scales,
      local_attempts = local_attempts,
      local_accepts = local_accepts,
      basin_rejections = basin_rejections,
      swap_attempts = swap_attempts,
      swap_accepts = swap_accepts,
      retained = retained,
      save_index = save_index,
      random_seed = get(".Random.seed", envir = .GlobalEnv)
    )
    temporary <- paste0(checkpoint_path, ".tmp_", Sys.getpid())
    saveRDS(object, temporary)
    if (!file.rename(temporary, checkpoint_path)) {
      stop("Failed to atomically update checkpoint: ", checkpoint_path)
    }
  }

  progress_path <- file.path(checkpoint_dir, paste0("progress_", chain_id, ".tsv"))
  write_progress <- function(iteration) {
    progress <- data.frame(
      chain = chain_id,
      family = family,
      replicate = replicate_id,
      iteration = iteration,
      n_iter = cfg$n_iter,
      warmup_complete = iteration >= cfg$burnin,
      retained_draws = save_index,
      family_target_version = family_contract$target_version,
      basin_boundary_rejections = sum(basin_rejections),
      minimum_swap_acceptance = min(
        swap_accepts / pmax(1L, swap_attempts)
      ),
      updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    )
    temporary <- paste0(progress_path, ".tmp_", Sys.getpid())
    utils::write.table(
      progress, temporary, sep = "\t", quote = FALSE,
      row.names = FALSE, na = "NA"
    )
    if (!file.rename(temporary, progress_path)) {
      stop("Failed to atomically update progress file: ", progress_path)
    }
  }

  for (iteration in seq.int(start_iter, cfg$n_iter)) {
    proposal_records <- vector("list", n_temp)
    for (temperature_index in seq_len(n_temp)) {
      kernel <- sample(kernel_names, 1L, prob = kernel_probabilities)
      if (identical(kernel, "pilot_cov")) {
        proposal <- reflect_unit(
          states[temperature_index, ] +
            as.numeric(
              stats::rnorm(n_dim) %*%
                proposal_cholesky_for_temperature(cfg$temperatures[[temperature_index]])
            ) * scales[temperature_index, kernel]
        )
      } else {
        dims <- switch(
          kernel,
          single = sample.int(n_dim, 1L),
          block8 = sample.int(n_dim, min(8L, n_dim)),
          full = seq_len(n_dim)
        )
        proposal <- states[temperature_index, ]
        proposal[dims] <- reflect_unit(
          proposal[dims] + stats::rnorm(
            length(dims), 0, scales[temperature_index, kernel]
          )
        )
      }
      proposal_records[[temperature_index]] <- list(
        kernel = kernel,
        proposal = proposal,
        in_family = identical(classify_family_state(proposal), family),
        log_acceptance_uniform = log(stats::runif(1L))
      )
    }

    if (!is.null(temperature_cluster)) {
      proposal_losses <- parallel::parLapply(
        temperature_cluster,
        proposal_records,
        function(record) {
          if (!record$in_family) {
            c(data_loss = 1e9, regularization_loss = 0, objective_t1 = 1e9)
          } else {
            evaluate_loss(record$proposal)
          }
        }
      )
    } else {
      proposal_losses <- lapply(
        proposal_records,
        function(record) {
          if (!record$in_family) {
            c(data_loss = 1e9, regularization_loss = 0, objective_t1 = 1e9)
          } else {
            evaluate_loss(record$proposal)
          }
        }
      )
    }
    if (length(proposal_losses) != n_temp ||
        any(vapply(proposal_losses, inherits, logical(1), what = "try-error"))) {
      stop("One or more temperature-parallel objective evaluations failed for ", chain_id)
    }

    for (temperature_index in seq_len(n_temp)) {
      record <- proposal_records[[temperature_index]]
      kernel <- record$kernel
      proposal <- record$proposal
      proposal_loss <- proposal_losses[[temperature_index]]
      if (!is.numeric(proposal_loss) ||
          !identical(names(proposal_loss), c(
            "data_loss", "regularization_loss", "objective_t1"
          ))) {
        stop("Malformed temperature-parallel objective result for ", chain_id)
      }
      temperature <- cfg$temperatures[[temperature_index]]
      current_target <- losses[temperature_index, "data_loss"] / temperature +
        losses[temperature_index, "regularization_loss"]
      proposal_target <- proposal_loss[["data_loss"]] / temperature +
        proposal_loss[["regularization_loss"]]
      log_alpha <- current_target - proposal_target
      accepted <- record$in_family && is.finite(log_alpha) &&
        record$log_acceptance_uniform < min(0, log_alpha)
      local_attempts[temperature_index, kernel] <- local_attempts[temperature_index, kernel] + 1L
      if (!record$in_family) {
        basin_rejections[[temperature_index]] <- basin_rejections[[temperature_index]] + 1L
      }
      if (accepted) {
        states[temperature_index, ] <- proposal
        losses[temperature_index, ] <- proposal_loss
        local_accepts[temperature_index, kernel] <- local_accepts[temperature_index, kernel] + 1L
      }
      if (iteration <= cfg$burnin &&
          is.finite(scales[temperature_index, kernel])) {
        gain <- 1 / (50 + local_attempts[temperature_index, kernel])^0.60
        scales[temperature_index, kernel] <- exp(
          log(scales[temperature_index, kernel]) +
            gain * (as.numeric(accepted) - kernel_targets[[kernel]])
        )
        scales[temperature_index, kernel] <- min(
          max(scales[temperature_index, kernel], minimum_scales[[kernel]]),
          maximum_scales[[kernel]]
        )
      }
    }

    first_pair <- if (iteration %% 2L == 1L) 1L else 2L
    for (left in seq.int(first_pair, n_temp - 1L, by = 2L)) {
      right <- left + 1L
      swap_attempts[[left]] <- swap_attempts[[left]] + 1L
      inv_temp_difference <- 1 / cfg$temperatures[[left]] - 1 / cfg$temperatures[[right]]
      log_alpha_swap <- inv_temp_difference * (
        losses[left, "data_loss"] - losses[right, "data_loss"]
      )
      if (is.finite(log_alpha_swap) && log(stats::runif(1L)) < min(0, log_alpha_swap)) {
        state_tmp <- states[left, ]
        states[left, ] <- states[right, ]
        states[right, ] <- state_tmp
        loss_tmp <- losses[left, ]
        losses[left, ] <- losses[right, ]
        losses[right, ] <- loss_tmp
        swap_accepts[[left]] <- swap_accepts[[left]] + 1L
      }
    }

    if (iteration > cfg$burnin && (iteration - cfg$burnin) %% cfg$thin == 0L) {
      save_index <- save_index + 1L
      for (temperature in cfg$display_temperatures) {
        temperature_index <- match(temperature, cfg$temperatures)
        key <- as.character(temperature)
        retained[[key]]$unit[save_index, ] <- states[temperature_index, ]
        retained[[key]]$loss[save_index, ] <- c(
          losses[temperature_index, "data_loss"],
          losses[temperature_index, "regularization_loss"],
          losses[temperature_index, "data_loss"] / temperature +
            losses[temperature_index, "regularization_loss"]
        )
      }
    }

    if (iteration %% cfg$checkpoint_every == 0L && iteration < cfg$n_iter) {
      checkpoint_state(iteration)
      write_progress(iteration)
    }
  }
  if (save_index != n_save) stop("Unexpected retained-draw count for ", chain_id)
  # Retain a final continuation checkpoint so a failed formal diagnostic stage
  # can be extended without restarting or changing the frozen production chain.
  checkpoint_state(cfg$n_iter)
  write_progress(cfg$n_iter)

  result <- list(
    config_signature = config_signature,
    target_version = family_contract$target_version,
    completed_iteration = cfg$n_iter,
    family = family,
    replicate = replicate_id,
    chain = chain_id,
    warmup_label = selection$warmup_label[selection$family == family],
    selected_seed = selection$selected_seed[selection$family == family],
    temperatures = cfg$temperatures,
    display_temperatures = cfg$display_temperatures,
    retained = retained,
    local_attempts = local_attempts,
    local_accepts = local_accepts,
    basin_rejections = basin_rejections,
    final_scales = scales,
    swap_attempts = swap_attempts,
    swap_accepts = swap_accepts
  )
  temporary <- paste0(result_path, ".tmp_", Sys.getpid())
  saveRDS(result, temporary)
  if (!file.rename(temporary, result_path)) stop("Failed to save completed ladder: ", result_path)
  result
}

ladder_plan <- expand.grid(
  family = c("C01", "C02", "C03"),
  replicate = seq_len(cfg$ladders_per_family),
  stringsAsFactors = FALSE
)
ladder_plan <- ladder_plan[order(
  match(ladder_plan$family, c("C01", "C02", "C03")),
  ladder_plan$replicate
), , drop = FALSE]

validate_ladder_result <- function(result, family, replicate_id) {
  expected_chain <- paste0(family, "_R", replicate_id)
  if (!is.list(result)) {
    stop("Completed ladder result is not a list for ", expected_chain, ".")
  }
  if (!identical(result$config_signature, config_signature)) {
    stop("Completed ladder configuration signature mismatch for ", expected_chain, ".")
  }
  if (!identical(result$target_version, family_contract$target_version)) {
    stop("Completed ladder target-version mismatch for ", expected_chain, ".")
  }
  if (!identical(result$family, family) ||
      !identical(as.integer(result$replicate), as.integer(replicate_id)) ||
      !identical(result$chain, expected_chain)) {
    stop("Completed ladder identity mismatch for ", expected_chain, ".")
  }
  if (!identical(as.integer(result$completed_iteration), as.integer(cfg$n_iter))) {
    stop("Completed ladder iteration mismatch for ", expected_chain, ".")
  }
  expected_temperatures <- as.character(cfg$display_temperatures)
  if (!identical(names(result$retained), expected_temperatures)) {
    stop("Completed ladder temperature keys mismatch for ", expected_chain, ".")
  }
  for (temperature in expected_temperatures) {
    retained <- result$retained[[temperature]]
    if (!is.list(retained) || !is.matrix(retained$unit) ||
        !is.matrix(retained$loss) || nrow(retained$unit) != n_save ||
        ncol(retained$unit) != nrow(coord_meta) ||
        nrow(retained$loss) != n_save || ncol(retained$loss) != 3L ||
        any(!is.finite(retained$unit)) || any(!is.finite(retained$loss))) {
      stop(
        "Completed ladder retained-state contract failed for ",
        expected_chain, " at T=", temperature, "."
      )
    }
    assignments <- apply(retained$unit, 1L, classify_family_state)
    if (any(assignments != family)) {
      stop(
        "Completed ladder contains a state outside the frozen basin for ",
        expected_chain, " at T=", temperature, "."
      )
    }
  }
  result
}

load_completed_ladders <- function() {
  lapply(seq_len(nrow(ladder_plan)), function(i) {
    family <- ladder_plan$family[[i]]
    replicate_id <- ladder_plan$replicate[[i]]
    chain <- paste0(family, "_R", replicate_id)
    path <- file.path(checkpoint_dir, paste0("ladder_", chain, "_complete.rds"))
    if (!file.exists(path)) stop("Missing completed ladder file: ", path)
    result <- tryCatch(readRDS(path), error = identity)
    if (inherits(result, "error")) {
      stop("Failed to read completed ladder ", chain, ": ", conditionMessage(result))
    }
    validate_ladder_result(result, family, replicate_id)
  })
}

if (single_ladder_mode) {
  message(
    "Running independent ladder ", ladder_family, "_R", ladder_replicate,
    " across ", length(cfg$temperatures), " temperatures using ",
    cfg$temperature_cores_per_ladder, " temperature workers."
  )
  ladder_result <- run_ladder(ladder_family, ladder_replicate)
  validate_ladder_result(ladder_result, ladder_family, ladder_replicate)
  message("Independent ladder completed and passed its file contract: ", ladder_result$chain)
  quit(save = "no", status = 0L, runLast = FALSE)
}

stop(
  "Cross-family aggregation is disabled for the family-conditioned targets. ",
  "Run aggregate_Figure5F_family_posterior.R separately for each family."
)

if (aggregate_only) {
  message("Loading and validating six independently completed ladder files.")
  ladders <- load_completed_ladders()
} else {
  if (cfg$cores > 1L) {
    stop(
      "Nested fork execution is disabled. Run each ladder independently with ",
      "--ladder_family and --ladder_replicate, then use --aggregate_only=TRUE."
    )
  }
  message(
    "Running ", nrow(ladder_plan),
    " ladders sequentially; production 30-CPU runs use independent Slurm tasks."
  )
  ladders <- lapply(seq_len(nrow(ladder_plan)), function(i) {
    result <- run_ladder(ladder_plan$family[[i]], ladder_plan$replicate[[i]])
    validate_ladder_result(result, ladder_plan$family[[i]], ladder_plan$replicate[[i]])
  })
}

ladder_failed <- vapply(
  ladders,
  function(x) is.null(x) || inherits(x, "try-error") || !is.list(x),
  logical(1)
)
if (any(ladder_failed)) {
  stop("One or more Figure 5F ladder results were absent or malformed.")
}

transform_natural <- function(value, transform) {
  if (identical(transform, "identity")) value else 10^value
}

draw_rows <- list()
state_records <- list()
draw_index <- 1L
state_index <- 1L
for (ladder in ladders) {
  for (temperature in cfg$display_temperatures) {
    key <- as.character(temperature)
    unit <- ladder$retained[[key]]$unit
    loss <- ladder$retained[[key]]$loss
    state_records[[state_index]] <- list(
      family = ladder$family,
      chain = ladder$chain,
      replicate = ladder$replicate,
      temperature = temperature,
      unit = unit,
      loss = loss
    )
    state_index <- state_index + 1L
    actual <- sweep(unit, 2L, coord_meta$width, "*")
    actual <- sweep(actual, 2L, coord_meta$lower, "+")
    colnames(actual) <- coord_meta$coordinate
    for (i in seq_len(nrow(soft_meta))) {
      parameter <- as.character(soft_meta$parameter[[i]])
      center_name <- as.character(soft_meta$center_name[[i]])
      transform <- as.character(soft_meta$transform[[i]])
      vivo_t <- actual[, paste0("vivo__", center_name)]
      vitro_t <- actual[, paste0("vitro__", center_name)]
      vivo_natural <- transform_natural(vivo_t, transform)
      vitro_natural <- transform_natural(vitro_t, transform)
      draw_rows[[draw_index]] <- data.frame(
        family = ladder$family,
        warmup_label = ladder$warmup_label,
        selected_seed = ladder$selected_seed,
        chain = ladder$chain,
        replicate = ladder$replicate,
        draw = seq_len(nrow(unit)),
        temperature = temperature,
        parameter = parameter,
        vivo_transformed = vivo_t,
        vitro_transformed = vitro_t,
        vivo_natural = vivo_natural,
        vitro_natural = vitro_natural,
        log2_ratio = log2(vivo_natural / vitro_natural),
        data_loss = loss[, "data_loss"],
        regularization_loss = loss[, "regularization_loss"],
        target_loss = loss[, "target_loss"],
        stringsAsFactors = FALSE
      )
      draw_index <- draw_index + 1L
    }
  }
}
draws <- do.call(rbind, draw_rows)
rownames(draws) <- NULL

expected_draw_rows <- nrow(ladder_plan) * length(cfg$display_temperatures) * n_save * 14L
if (nrow(draws) != expected_draw_rows || any(!is.finite(draws$log2_ratio)) ||
    any(draws$vivo_natural <= 0) || any(draws$vitro_natural <= 0)) {
  stop("Generated Figure 5F draw table fails its structural contract.")
}

split_chain_matrix <- function(x) {
  x <- as.matrix(x)
  half <- floor(nrow(x) / 2L)
  if (half < 3L || ncol(x) < 2L) stop("Insufficient chains for split diagnostics.")
  cbind(
    x[seq_len(half), , drop = FALSE],
    x[seq.int(nrow(x) - half + 1L, nrow(x)), , drop = FALSE]
  )
}

rank_normalize_matrix <- function(x) {
  flat <- as.vector(x)
  ranks <- rank(flat, ties.method = "average")
  z <- stats::qnorm((ranks - 3 / 8) / (length(ranks) + 1 / 4))
  matrix(z, nrow = nrow(x), ncol = ncol(x))
}

basic_rhat <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- nrow(x)
  chain_variances <- apply(x, 2L, stats::var)
  within <- mean(chain_variances)
  if (!is.finite(within) || within <= 0) {
    return(if (stats::sd(as.vector(x)) == 0) 1 else Inf)
  }
  between <- n * stats::var(colMeans(x))
  variance_plus <- (n - 1) / n * within + between / n
  sqrt(variance_plus / within)
}

basic_ess <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- nrow(x)
  m <- ncol(x)
  total <- n * m
  chain_variances <- apply(x, 2L, stats::var)
  within <- mean(chain_variances)
  between <- n * stats::var(colMeans(x))
  variance_plus <- (n - 1) / n * within + between / n
  if (!is.finite(variance_plus) || variance_plus <= 0 || !is.finite(within)) {
    return(if (stats::sd(as.vector(x)) == 0) total else NA_real_)
  }
  autocovariance <- vapply(
    seq_len(m),
    function(chain_i) {
      as.numeric(stats::acf(
        x[, chain_i], lag.max = n - 1L, type = "covariance",
        plot = FALSE, demean = TRUE
      )$acf)
    },
    numeric(n)
  )
  mean_autocovariance <- rowMeans(autocovariance)
  rho <- 1 - (within - mean_autocovariance) / variance_plus
  rho[[1L]] <- 1
  max_pair <- floor((length(rho) - 1L) / 2L)
  pair_sums <- numeric(0)
  if (max_pair >= 0L) {
    for (pair_i in 0:max_pair) {
      left <- 2L * pair_i + 1L
      right <- left + 1L
      if (right > length(rho)) break
      value <- rho[[left]] + rho[[right]]
      if (!is.finite(value) || value < 0) break
      pair_sums <- c(pair_sums, value)
    }
  }
  if (!length(pair_sums)) return(min(total, m))
  pair_sums <- cummin(pair_sums)
  tau <- -1 + 2 * sum(pair_sums)
  if (!is.finite(tau) || tau <= 0) return(total)
  min(total, total / tau)
}

base_draw_summary <- function(arr) {
  variables <- dimnames(arr)[[3L]]
  if (is.null(variables)) variables <- paste0("V", seq_len(dim(arr)[[3L]]))
  rows <- lapply(seq_along(variables), function(variable_i) {
    original <- arr[, , variable_i, drop = TRUE]
    if (is.null(dim(original))) original <- matrix(original, ncol = dim(arr)[[2L]])
    split <- split_chain_matrix(original)
    rank_normalized <- rank_normalize_matrix(split)
    folded_rank_normalized <- rank_normalize_matrix(
      abs(split - stats::median(as.vector(split)))
    )
    q_tail <- stats::quantile(as.vector(split), c(0.05, 0.95), names = FALSE, type = 8)
    ess_bulk <- basic_ess(rank_normalized)
    ess_tail <- min(
      basic_ess(split <= q_tail[[1L]]),
      basic_ess(split >= q_tail[[2L]])
    )
    values <- as.vector(original)
    data.frame(
      variable = variables[[variable_i]],
      mean = mean(values),
      sd = stats::sd(values),
      rhat = max(basic_rhat(rank_normalized), basic_rhat(folded_rank_normalized)),
      ess_bulk = ess_bulk,
      ess_tail = ess_tail,
      mcse_mean = stats::sd(values) / sqrt(ess_bulk),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_array <- function(arr, scope, temperature, family = "all") {
  out <- base_draw_summary(arr)
  out$scope <- scope
  out$temperature <- temperature
  out$family <- family
  out
}

diagnostic_rows <- list()
diagnostic_index <- 1L
for (temperature in cfg$display_temperatures) {
  state_subset <- state_records[vapply(
    state_records,
    function(x) identical(as.numeric(x$temperature), as.numeric(temperature)),
    logical(1)
  )]
  chain_names <- vapply(state_subset, `[[`, character(1), "chain")
  state_array <- array(
    NA_real_,
    dim = c(n_save, length(state_subset), nrow(coord_meta)),
    dimnames = list(NULL, chain_names, coord_meta$coordinate)
  )
  for (chain_i in seq_along(state_subset)) {
    state_array[, chain_i, ] <- state_subset[[chain_i]]$unit
  }
  diagnostic_rows[[diagnostic_index]] <- summarize_array(
    state_array, "all_40_sampling_coordinates", temperature
  )
  diagnostic_index <- diagnostic_index + 1L

  temp_draws <- draws[draws$temperature == temperature, , drop = FALSE]
  parameter_levels <- as.character(soft_meta$parameter)
  ratio_array <- array(
    NA_real_,
    dim = c(n_save, length(chain_names), length(parameter_levels)),
    dimnames = list(NULL, chain_names, parameter_levels)
  )
  for (chain_i in seq_along(chain_names)) {
    chain_data <- temp_draws[temp_draws$chain == chain_names[[chain_i]], , drop = FALSE]
    for (parameter_i in seq_along(parameter_levels)) {
      values <- chain_data$log2_ratio[chain_data$parameter == parameter_levels[[parameter_i]]]
      ratio_array[, chain_i, parameter_i] <- values[order(
        chain_data$draw[chain_data$parameter == parameter_levels[[parameter_i]]]
      )]
    }
  }
  diagnostic_rows[[diagnostic_index]] <- summarize_array(
    ratio_array, "paired_log2_ratios", temperature
  )
  diagnostic_index <- diagnostic_index + 1L

  for (family in c("C01", "C02", "C03")) {
    family_chains <- grep(paste0("^", family, "_"), chain_names)
    diagnostic_rows[[diagnostic_index]] <- summarize_array(
      ratio_array[, family_chains, , drop = FALSE],
      "paired_log2_ratios_within_family",
      temperature,
      family
    )
    diagnostic_index <- diagnostic_index + 1L
  }
}
diagnostics <- do.call(rbind, diagnostic_rows)
rownames(diagnostics) <- NULL

sampler_rows <- list()
sampler_index <- 1L
for (ladder in ladders) {
  for (temperature_index in seq_along(cfg$temperatures)) {
    for (kernel in kernel_names) {
      sampler_rows[[sampler_index]] <- data.frame(
        family = ladder$family,
        chain = ladder$chain,
        temperature = cfg$temperatures[[temperature_index]],
        diagnostic = "local_proposal",
        kernel = kernel,
        attempts = ladder$local_attempts[temperature_index, kernel],
        accepts = ladder$local_accepts[temperature_index, kernel],
        acceptance_rate = ladder$local_accepts[temperature_index, kernel] /
          max(1, ladder$local_attempts[temperature_index, kernel]),
        final_scale = ladder$final_scales[temperature_index, kernel],
        stringsAsFactors = FALSE
      )
      sampler_index <- sampler_index + 1L
    }
  }
  for (left in seq_len(length(cfg$temperatures) - 1L)) {
    sampler_rows[[sampler_index]] <- data.frame(
      family = ladder$family,
      chain = ladder$chain,
      temperature = paste0(cfg$temperatures[[left]], "<->", cfg$temperatures[[left + 1L]]),
      diagnostic = "replica_swap",
      kernel = "adjacent",
      attempts = ladder$swap_attempts[[left]],
      accepts = ladder$swap_accepts[[left]],
      acceptance_rate = ladder$swap_accepts[[left]] / max(1, ladder$swap_attempts[[left]]),
      final_scale = NA_real_,
      stringsAsFactors = FALSE
    )
    sampler_index <- sampler_index + 1L
  }
}
sampler_diagnostics <- do.call(rbind, sampler_rows)

quantile_summary <- function(x) {
  q <- stats::quantile(x, c(0.05, 0.25, 0.50, 0.75, 0.95), names = FALSE, type = 8)
  c(q05 = q[[1L]], q25 = q[[2L]], median = q[[3L]], q75 = q[[4L]], q95 = q[[5L]])
}

family_groups <- split(
  draws,
  interaction(draws$temperature, draws$parameter, draws$family, drop = TRUE)
)
family_summary <- do.call(rbind, lapply(family_groups, function(group) {
  q <- quantile_summary(group$log2_ratio)
  data.frame(
    temperature = group$temperature[[1L]],
    parameter = group$parameter[[1L]],
    family = group$family[[1L]],
    n_draws = nrow(group),
    q05 = q[["q05"]],
    q25 = q[["q25"]],
    median = q[["median"]],
    q75 = q[["q75"]],
    q95 = q[["q95"]],
    width90 = q[["q95"]] - q[["q05"]],
    direction = if (q[["q05"]] > 0) {
      "higher in vivo"
    } else if (q[["q95"]] < 0) {
      "lower in vivo"
    } else {
      "overlaps equality"
    },
    stringsAsFactors = FALSE
  )
}))
rownames(family_summary) <- NULL

prior_draws <- utils::read.delim(prior_draw_path, check.names = FALSE, stringsAsFactors = FALSE)
prior_groups <- split(prior_draws, prior_draws$parameter)
prior_summary <- do.call(rbind, lapply(prior_groups, function(group) {
  q <- quantile_summary(group$log2_ratio)
  data.frame(
    parameter = group$parameter[[1L]],
    prior_q05 = q[["q05"]],
    prior_median = q[["median"]],
    prior_q95 = q[["q95"]],
    prior_width90 = q[["q95"]] - q[["q05"]],
    stringsAsFactors = FALSE
  )
}))
family_summary <- merge(family_summary, prior_summary, by = "parameter", all.x = TRUE, sort = FALSE)
family_summary$contraction_ratio90 <- family_summary$width90 / family_summary$prior_width90
family_summary <- family_summary[order(
  family_summary$temperature,
  match(family_summary$parameter, as.character(soft_meta$parameter)),
  match(family_summary$family, c("C01", "C02", "C03"))
), , drop = FALSE]

primary_ratio_diag <- diagnostics[
  diagnostics$scope == "paired_log2_ratios" & diagnostics$temperature == 1,
  , drop = FALSE
]
primary_state_diag <- diagnostics[
  diagnostics$scope == "all_40_sampling_coordinates" & diagnostics$temperature == 1,
  , drop = FALSE
]
within_family_diag <- diagnostics[
  diagnostics$scope == "paired_log2_ratios_within_family" & diagnostics$temperature == 1,
  , drop = FALSE
]
all_finite_le <- function(x, threshold) {
  length(x) > 0L && all(is.finite(x)) && all(x <= threshold)
}
all_finite_ge <- function(x, threshold) {
  length(x) > 0L && all(is.finite(x)) && all(x >= threshold)
}
finite_max_or_na <- function(x) {
  if (!length(x) || any(!is.finite(x))) NA_real_ else max(x)
}
finite_min_or_na <- function(x) {
  if (!length(x) || any(!is.finite(x))) NA_real_ else min(x)
}
local_acceptance <- sampler_diagnostics$acceptance_rate[
  sampler_diagnostics$diagnostic == "local_proposal"
]
swap_acceptance <- sampler_diagnostics$acceptance_rate[
  sampler_diagnostics$diagnostic == "replica_swap"
]
map_jump_diagnostics <- sampler_diagnostics[
  sampler_diagnostics$diagnostic == "local_proposal" &
    sampler_diagnostics$kernel == "map_jump",
  , drop = FALSE
]
map_jump_accepts_by_chain <- tapply(
  map_jump_diagnostics$accepts,
  map_jump_diagnostics$chain,
  sum
)
map_jump_hot_acceptance <- map_jump_diagnostics$acceptance_rate[
  as.character(map_jump_diagnostics$temperature) == "8"
]
convergence_checks <- data.frame(
  check = c(
    "primary_ratio_rhat_max_le_1p01",
    "primary_ratio_ess_bulk_min_ge_400",
    "primary_ratio_ess_tail_min_ge_400",
    "primary_state_rhat_max_le_1p05",
    "primary_state_ess_bulk_min_ge_100",
    "within_family_ratio_rhat_max_le_1p05",
    "within_family_ratio_ess_bulk_min_ge_100",
    "local_acceptance_finite",
    "map_jump_accepts_min_per_chain_ge_10",
    "map_jump_T8_acceptance_min_ge_0p01",
    "adjacent_swap_acceptance_min_ge_0p05"
  ),
  observed = c(
    finite_max_or_na(primary_ratio_diag$rhat),
    finite_min_or_na(primary_ratio_diag$ess_bulk),
    finite_min_or_na(primary_ratio_diag$ess_tail),
    finite_max_or_na(primary_state_diag$rhat),
    finite_min_or_na(primary_state_diag$ess_bulk),
    finite_max_or_na(within_family_diag$rhat),
    finite_min_or_na(within_family_diag$ess_bulk),
    length(local_acceptance) > 0L && all(is.finite(local_acceptance)),
    finite_min_or_na(map_jump_accepts_by_chain),
    finite_min_or_na(map_jump_hot_acceptance),
    finite_min_or_na(swap_acceptance)
  ),
  passed = c(
    all_finite_le(primary_ratio_diag$rhat, 1.01),
    all_finite_ge(primary_ratio_diag$ess_bulk, 400),
    all_finite_ge(primary_ratio_diag$ess_tail, 400),
    all_finite_le(primary_state_diag$rhat, 1.05),
    all_finite_ge(primary_state_diag$ess_bulk, 100),
    all_finite_le(within_family_diag$rhat, 1.05),
    all_finite_ge(within_family_diag$ess_bulk, 100),
    length(local_acceptance) > 0L && all(is.finite(local_acceptance)),
    all_finite_ge(map_jump_accepts_by_chain, 10),
    all_finite_ge(map_jump_hot_acceptance, 0.01),
    all_finite_ge(swap_acceptance, 0.05)
  ),
  stringsAsFactors = FALSE
)
safe_for_figure <- all(convergence_checks$passed)

end_md5 <- c(
  generator = unname(tools::md5sum(generator_path)),
  selection = unname(tools::md5sum(selection_path)),
  map_audit = unname(tools::md5sum(map_audit_path)),
  prior_draws = unname(tools::md5sum(prior_draw_path)),
  proposal_covariance = unname(tools::md5sum(proposal_covariance_path)),
  proposal_calibration_script = unname(tools::md5sum(proposal_calibration_script_path))
)
start_md5 <- c(
  generator = generator_md5_at_start,
  selection = selection_md5_at_start,
  map_audit = map_audit_md5_at_start,
  prior_draws = prior_draw_md5_at_start,
  proposal_covariance = proposal_covariance_md5_at_start,
  proposal_calibration_script = proposal_calibration_script_md5_at_start
)
if (!identical(end_md5, start_md5)) {
  stop(
    "Figure 5F sampler code or an input changed while chains were running; ",
    "results are not being released."
  )
}

configuration <- data.frame(
  key = c(
    "target_name", "code_commit", "mode", "n_iter", "burnin", "thin",
    "retained_draws_per_chain_temperature", "ladders_per_family",
    "families", "temperature_ladder", "display_temperatures", "cores",
    "temperature_cores_per_ladder", "maximum_parallel_objective_evaluations",
    "parallel_execution_model", "base_seed", "parameter_dimension", "paired_parameter_count",
    "proposal_kernels", "proposal_probabilities",
    "map_jump_noise_scales", "map_jump_reverse_symmetric",
    "pilot_covariance_role", "pilot_covariance_md5",
    "pilot_calibration_script_md5", "production_chain_extendable",
    "generator_md5", "selection_md5", "map_replay_audit_md5", "prior_draw_md5",
    "safe_for_primary_figure"
  ),
  value = c(
    paste0(
      "generalized posterior: exp(-L_data/T) x exp(-L_regularization) ",
      "x I(bounded feasible support)"
    ),
    "83953a874401e42cd176432786f889a896adc959",
    cfg$mode, cfg$n_iter, cfg$burnin, cfg$thin, n_save,
    cfg$ladders_per_family, "C01,C02,C03",
    paste(cfg$temperatures, collapse = ","),
    paste(cfg$display_temperatures, collapse = ","),
    cfg$cores, cfg$temperature_cores_per_ladder,
    min(cfg$cores, 3L * cfg$ladders_per_family) *
      min(cfg$temperature_cores_per_ladder, length(cfg$temperatures)),
    "six independent Slurm ladder tasks x five forked temperature workers; file-backed aggregation",
    cfg$base_seed, nrow(coord_meta), nrow(soft_meta),
    paste(kernel_names, collapse = ","),
    paste(names(kernel_probabilities), kernel_probabilities, sep = "=", collapse = ","),
    paste(map_jump_noise_scales, collapse = ","), all(offset_has_reverse),
    "fixed symmetric proposal preconditioner; no pilot draws enter target or output",
    proposal_covariance_md5_at_start,
    proposal_calibration_script_md5_at_start,
    TRUE,
    generator_md5_at_start,
    selection_md5_at_start,
    map_audit_md5_at_start,
    prior_draw_md5_at_start,
    safe_for_figure
  ),
  stringsAsFactors = FALSE
)

draws_path <- file.path(output_dir, "figure5f_generalized_posterior_draws.tsv")
write_tsv(draws, draws_path)
write_tsv(diagnostics, file.path(output_dir, "figure5f_generalized_posterior_diagnostics.tsv"))
write_tsv(sampler_diagnostics, file.path(output_dir, "figure5f_sampler_diagnostics.tsv"))
write_tsv(family_summary, file.path(output_dir, "figure5f_generalized_posterior_family_summary.tsv"))
write_tsv(convergence_checks, file.path(output_dir, "figure5f_generalized_posterior_readiness.tsv"))
write_tsv(configuration, file.path(output_dir, "figure5f_generalized_posterior_config.tsv"))
saveRDS(
  list(
    config = cfg,
    config_signature = config_signature,
    coordinate_metadata = coord_meta,
    state_records = state_records,
    ladders = ladders
  ),
  file.path(output_dir, "figure5f_generalized_posterior_state_draws.rds")
)

primary_draws <- draws[draws$temperature == 1, , drop = FALSE]
if (safe_for_figure) {
  write_tsv(primary_draws, file.path(figure5_dir, "figure5f_posterior_draws.tsv"))
  message("Figure 5F generalized posterior passed all readiness checks.")
} else {
  message(
    "Figure 5F draws were saved, but one or more convergence/readiness checks failed. ",
    "The primary figure draw table was not released."
  )
}
print(convergence_checks)
quit(save = "no", status = if (safe_for_figure || mode == "pilot") 0L else 2L, runLast = FALSE)
