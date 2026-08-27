#!/usr/bin/env Rscript

# Small, deterministic parameter replay for the in-vitro structure fix.
#
# This script intentionally:
# - evaluates only seed 10 and seed 340;
# - reuses their existing parameter vectors without optimization;
# - writes only below revised/add/invitro_structure_fix;
# - defaults to results/ and accepts a safe alternate directory name through
#   INVITRO_STRUCTURE_FIX_RESULTS_DIR;
# - refuses to overwrite an existing results directory.

soft_coupling_root <- normalizePath(
  "/Users/4482173/Documents/GitHub/soft_coupling",
  mustWork = TRUE
)
hypoxia_figures_root <- normalizePath(
  "/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures",
  mustWork = TRUE
)
delivery_root <- file.path(
  hypoxia_figures_root,
  "revised",
  "add",
  "invitro_structure_fix"
)
results_dir_name <- Sys.getenv(
  "INVITRO_STRUCTURE_FIX_RESULTS_DIR",
  unset = "results"
)
if (!nzchar(results_dir_name) ||
    !identical(results_dir_name, basename(results_dir_name)) ||
    results_dir_name %in% c(".", "..")) {
  stop(
    "INVITRO_STRUCTURE_FIX_RESULTS_DIR must be one safe directory name.",
    call. = FALSE
  )
}
results_root <- file.path(delivery_root, results_dir_name)

if (dir.exists(results_root) || file.exists(results_root)) {
  stop(
    "Refusing to overwrite existing validation results: ",
    results_root,
    call. = FALSE
  )
}

temporary_results_root <- file.path(
  delivery_root,
  paste0(".", results_dir_name, ".tmp-", Sys.getpid())
)
if (dir.exists(temporary_results_root) || file.exists(temporary_results_root)) {
  stop(
    "Temporary validation path already exists: ",
    temporary_results_root,
    call. = FALSE
  )
}
dir.create(temporary_results_root, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

oxygen_root <- file.path(soft_coupling_root, "oxygen")
workflow_root <- file.path(
  oxygen_root,
  "code",
  "O2_supply_demand_MAP"
)
shared_path <- file.path(
  workflow_root,
  "util",
  "o2_supply_demand_map_shared.R"
)
common_path <- file.path(
  workflow_root,
  "util",
  "o2_supply_demand_map_common_semantics.R"
)
loader_path <- file.path(
  workflow_root,
  "util",
  "o2_supply_demand_map_invitro_utils.R"
)
model_path <- file.path(
  workflow_root,
  "model",
  "model_O2_supply_demand_MAP.R"
)

support_env <- new.env(parent = globalenv())
sys.source(shared_path, envir = support_env, chdir = TRUE)
sys.source(common_path, envir = support_env, chdir = TRUE)
previous_oxygen_code_dir <- Sys.getenv(
  "MININGCLONEID_OXYGEN_CODE_DIR",
  unset = ""
)
Sys.setenv(
  MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path)
)
sys.source(model_path, envir = support_env, chdir = TRUE)
if (nzchar(previous_oxygen_code_dir)) {
  Sys.setenv(
    MININGCLONEID_OXYGEN_CODE_DIR = previous_oxygen_code_dir
  )
} else {
  Sys.unsetenv("MININGCLONEID_OXYGEN_CODE_DIR")
}
api_env <- new.env(parent = support_env)
source(loader_path, local = api_env, chdir = TRUE)

parameter_table <- file.path(
  oxygen_root,
  "data",
  "O2_supply_demand",
  "parameter_table_invitro_buffering.csv"
)
fit_objects_dir <- file.path(
  oxygen_root,
  "ploidyOxygen",
  "data",
  "fit_objects"
)
flow_density_path <- file.path(
  oxygen_root,
  "data",
  "g0g1_ploidy_density_grid.csv"
)

cfg <- api_env$ivt_build_default_cfg(
  repo_root = oxygen_root,
  dt = 0.05,
  init_total_size = 1e6,
  o2_upper_bound = 21,
  fixed_oxygen = TRUE
)
cfg$parameter_table <- normalizePath(parameter_table, mustWork = TRUE)
cfg <- get(
  "normalize_sim_cfg_common",
  envir = api_env,
  inherits = TRUE
)(cfg, context = "fit")

fit_objects <- api_env$ivt_load_fit_objects(
  repo_root = oxygen_root,
  fit_objects_dir = fit_objects_dir,
  flow_csv_path = flow_density_path
)

seed_parameter_paths <- c(
  `10` = file.path(
    oxygen_root,
    "results",
    "fit_invitro_O2_buffering_500seed",
    "seed10",
    "best_params.tsv"
  ),
  `340` = file.path(
    oxygen_root,
    "results",
    "fit_invitro_O2_buffering_500seed",
    "seed340",
    "best_params.tsv"
  )
)
if (any(!file.exists(seed_parameter_paths))) {
  stop(
    "Missing local parameter replay input(s): ",
    paste(seed_parameter_paths[!file.exists(seed_parameter_paths)], collapse = ", "),
    call. = FALSE
  )
}

read_replay_params <- function(path) {
  parameter_df <- utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("parameter", "value")
  missing <- setdiff(required, names(parameter_df))
  if (length(missing) > 0L) {
    stop(
      "Parameter table is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(parameter_df$parameter)) {
    stop("Parameter table contains duplicate parameter names: ", path, call. = FALSE)
  }
  values <- stats::setNames(
    as.numeric(parameter_df$value),
    as.character(parameter_df$parameter)
  )
  if (any(!is.finite(values))) {
    stop("Parameter table contains non-finite values: ", path, call. = FALSE)
  }
  run_params <- api_env$ivt_load_default_run_params(cfg)
  for (parameter_name in names(values)) {
    run_params[[parameter_name]] <- unname(values[[parameter_name]])
  }
  run_params <- get(
    "normalize_run_params_common",
    envir = api_env,
    inherits = TRUE
  )(run_params, cfg = cfg)
  list(
    run_params = run_params,
    input_parameters = values
  )
}

passage_audit_columns <- c(
  "seed",
  "cohort",
  "lineage_id",
  "scenario_id",
  "passage_id",
  "passage_index",
  "passage_duration",
  "endpoint_day",
  "selected_day",
  "closest_day_diagnostic",
  "predicted_initial_cells",
  "predicted_final_cells",
  "observed_initial_cells",
  "observed_final_cells",
  "predicted_growth",
  "observed_growth",
  "passage_recorded",
  "available_cells",
  "required_cells",
  "supply_ratio",
  "reseed_mode",
  "boundary_scale",
  "cell_number_before",
  "cell_number_after",
  "cumulative_time"
)

collect_boundary_outputs <- function(comp, seed) {
  runs <- list(comp$run_2N, comp$run_4N)
  segment_results <- unlist(
    lapply(runs, function(run) run$segment_results),
    recursive = FALSE,
    use.names = FALSE
  )
  boundary_checks <- dplyr::bind_rows(lapply(segment_results, function(result) {
    segment <- result$segment
    selection <- result$selection
    endpoint_state <- as.numeric(selection$endpoint_state)
    reseeded_state <- as.numeric(selection$reseeded_state)
    expected_state <- if (identical(
      selection$reseed_mode,
      "downsample_to_observed_inoculum"
    )) {
      endpoint_state * as.numeric(selection$boundary_scale)
    } else {
      endpoint_state
    }
    max_abs_state_difference <- if (length(endpoint_state) > 0L) {
      max(abs(reseeded_state - expected_state))
    } else {
      NA_real_
    }
    data.frame(
      seed = as.integer(seed),
      cohort = as.character(segment$cohort),
      lineage_id = as.character(segment$lineage_id),
      scenario_id = as.character(segment$scenario_id),
      passage_id = as.character(segment$passage_id),
      passage_recorded = isTRUE(selection$passage_recorded),
      reseed_mode = as.character(selection$reseed_mode),
      available_cells = as.numeric(selection$available_cells),
      required_cells = as.numeric(selection$required_cells),
      supply_ratio = as.numeric(selection$supply_ratio),
      boundary_scale = as.numeric(selection$boundary_scale),
      cell_number_before = as.numeric(selection$cell_number_before),
      cell_number_after = as.numeric(selection$cell_number_after),
      n_state_components = length(endpoint_state),
      max_abs_state_difference = max_abs_state_difference,
      insufficient_state_identical = if (identical(
        selection$reseed_mode,
        "carry_forward_insufficient"
      )) {
        identical(endpoint_state, reseeded_state)
      } else {
        NA
      },
      proportional_downsample_exact = if (identical(
        selection$reseed_mode,
        "downsample_to_observed_inoculum"
      )) {
        identical(reseeded_state, expected_state)
      } else {
        NA
      },
      boundary_scale_at_most_one = isTRUE(
        as.numeric(selection$boundary_scale) <= 1
      ),
      stringsAsFactors = FALSE
    )
  }))

  state_vectors <- dplyr::bind_rows(lapply(segment_results, function(result) {
    segment <- result$segment
    selection <- result$selection
    data.frame(
      seed = as.integer(seed),
      cohort = as.character(segment$cohort),
      lineage_id = as.character(segment$lineage_id),
      scenario_id = as.character(segment$scenario_id),
      passage_id = as.character(segment$passage_id),
      chromosome_number = as.numeric(
        if (identical(segment$cohort, "2N")) {
          comp$run_2N$grid_pre
        } else {
          comp$run_4N$grid_pre
        }
      ),
      endpoint_state = as.numeric(selection$endpoint_state),
      reseeded_state = as.numeric(selection$reseeded_state),
      reseed_mode = as.character(selection$reseed_mode),
      boundary_scale = as.numeric(selection$boundary_scale),
      stringsAsFactors = FALSE
    )
  }))

  list(
    checks = boundary_checks,
    state_vectors = state_vectors
  )
}

collect_scenario_time <- function(summary_df, seed) {
  grouping <- unique(summary_df[, c(
    "cohort",
    "lineage_id",
    "scenario_id"
  )])
  durations <- stats::aggregate(
    passage_duration ~ cohort + lineage_id + scenario_id,
    data = summary_df,
    FUN = sum
  )
  names(durations)[names(durations) == "passage_duration"] <- "cumulative_passage_days"
  endpoint_days <- stats::aggregate(
    endpoint_day ~ cohort + lineage_id + scenario_id,
    data = summary_df,
    FUN = sum
  )
  names(endpoint_days)[names(endpoint_days) == "endpoint_day"] <- "cumulative_endpoint_days"
  selected_days <- stats::aggregate(
    selected_day ~ cohort + lineage_id + scenario_id,
    data = summary_df,
    FUN = sum
  )
  names(selected_days)[names(selected_days) == "selected_day"] <- "cumulative_selected_days"
  closest_days <- stats::aggregate(
    closest_day_diagnostic ~ cohort + lineage_id + scenario_id,
    data = summary_df,
    FUN = sum
  )
  names(closest_days)[names(closest_days) == "closest_day_diagnostic"] <-
    "diagnostic_closest_days_total"
  out <- Reduce(
    function(left, right) merge(
      left,
      right,
      by = c("cohort", "lineage_id", "scenario_id"),
      all = TRUE,
      sort = FALSE
    ),
    list(grouping, durations, endpoint_days, selected_days, closest_days)
  )
  out$seed <- as.integer(seed)
  out[, c("seed", setdiff(names(out), "seed")), drop = FALSE]
}

collect_likelihood_uniqueness <- function(comp, seed) {
  tables <- list(
    growth = comp$growth_df,
    karyotype = comp$ploidy_df,
    flow = comp$flow_df
  )
  dplyr::bind_rows(lapply(names(tables), function(modality) {
    table <- tables[[modality]]
    passage_ids <- as.character(table$passage_id)
    data.frame(
      seed = as.integer(seed),
      modality = modality,
      n_likelihood_rows = nrow(table),
      n_unique_passage_ids = length(unique(passage_ids)),
      n_duplicate_passage_ids = sum(duplicated(passage_ids)),
      every_passage_id_enters_once = !anyDuplicated(passage_ids),
      stringsAsFactors = FALSE
    )
  }))
}

collect_shared_parameter_check <- function(comp,
                                           seed,
                                           input_parameters,
                                           run_params) {
  scenario_ids <- sort(unique(comp$summary$scenario_id))
  run_by_cohort <- list(
    `2N` = comp$run_2N,
    `4N` = comp$run_4N
  )
  reference_values <- vapply(
    names(input_parameters),
    function(parameter_name) as.numeric(run_params[[parameter_name]]),
    numeric(1)
  )
  dplyr::bind_rows(lapply(scenario_ids, function(scenario_id) {
    cohort <- sub("-.*$", "", scenario_id)
    scenario_values <- vapply(
      names(input_parameters),
      function(parameter_name) {
        as.numeric(run_by_cohort[[cohort]]$shared_run_params[[parameter_name]])
      },
      numeric(1)
    )
    data.frame(
      seed = as.integer(seed),
      scenario_id = scenario_id,
      parameter = names(input_parameters),
      value = scenario_values,
      reference_value = reference_values,
      identical_to_shared_reference = scenario_values == reference_values,
      stringsAsFactors = FALSE
    )
  }))
}

overview_rows <- list()
scenario_time_rows <- list()
uniqueness_rows <- list()
shared_parameter_rows <- list()

for (seed_name in names(seed_parameter_paths)) {
  seed <- as.integer(seed_name)
  replay <- read_replay_params(seed_parameter_paths[[seed_name]])
  comp <- api_env$ivt_objective_components(
    run_params = replay$run_params,
    fit_objects = fit_objects,
    cfg = cfg,
    fallback_max_passage_days = 14,
    growth_weight = 1,
    ploidy_weight = 1,
    flow_weight = 1
  )

  seed_output_dir <- file.path(
    temporary_results_root,
    paste0("seed", seed)
  )
  dir.create(seed_output_dir, recursive = TRUE, showWarnings = FALSE)

  summary_df <- comp$summary
  summary_df$seed <- seed
  summary_df <- summary_df[, c(
    "seed",
    setdiff(names(summary_df), "seed")
  ), drop = FALSE]
  missing_audit_columns <- setdiff(passage_audit_columns, names(summary_df))
  if (length(missing_audit_columns) > 0L) {
    stop(
      "Passage audit is missing required columns: ",
      paste(missing_audit_columns, collapse = ", "),
      call. = FALSE
    )
  }
  passage_audit <- summary_df[, passage_audit_columns, drop = FALSE]

  boundary <- collect_boundary_outputs(comp, seed)
  scenario_time <- collect_scenario_time(summary_df, seed)
  uniqueness <- collect_likelihood_uniqueness(comp, seed)
  shared_parameter_check <- collect_shared_parameter_check(
    comp = comp,
    seed = seed,
    input_parameters = replay$input_parameters,
    run_params = replay$run_params
  )

  write_tsv(
    summary_df,
    file.path(seed_output_dir, "invitro_lineage_summary.tsv")
  )
  write_tsv(
    passage_audit,
    file.path(seed_output_dir, "invitro_passage_audit.tsv")
  )
  write_tsv(
    comp$growth_df,
    file.path(seed_output_dir, "invitro_growth_loglik.tsv")
  )
  write_tsv(
    comp$ploidy_df,
    file.path(seed_output_dir, "invitro_karyotype_loglik.tsv")
  )
  write_tsv(
    comp$flow_df,
    file.path(seed_output_dir, "invitro_flow_loglik.tsv")
  )
  write_tsv(
    comp$objective_hierarchy,
    file.path(seed_output_dir, "invitro_objective_hierarchy.tsv")
  )
  write_tsv(
    boundary$checks,
    file.path(seed_output_dir, "passage_boundary_checks.tsv")
  )
  write_tsv(
    boundary$state_vectors,
    file.path(seed_output_dir, "passage_boundary_state_vectors.tsv")
  )
  write_tsv(
    shared_parameter_check,
    file.path(seed_output_dir, "shared_parameter_check.tsv")
  )

  overview_rows[[seed_name]] <- data.frame(
    seed = seed,
    parameter_source = normalizePath(
      seed_parameter_paths[[seed_name]],
      mustWork = TRUE
    ),
    objective = as.numeric(comp$objective),
    total_loglik = as.numeric(comp$total_loglik),
    growth_loglik = as.numeric(comp$growth_loglik),
    karyotype_loglik = as.numeric(comp$ploidy_loglik),
    flow_loglik = as.numeric(comp$flow_loglik),
    n_passage_predictions = nrow(summary_df),
    n_scenarios = as.integer(comp$n_scenarios),
    n_growth_likelihood_units = nrow(comp$growth_df),
    n_karyotype_likelihood_units = nrow(comp$ploidy_df),
    n_flow_likelihood_units = nrow(comp$flow_df),
    n_insufficient_boundaries = as.integer(comp$n_insufficient_boundaries),
    n_downsample_boundaries = sum(
      summary_df$reseed_mode == "downsample_to_observed_inoculum"
    ),
    n_terminal_passages = sum(
      summary_df$reseed_mode == "terminal_no_reseed"
    ),
    max_boundary_scale = as.numeric(comp$max_boundary_scale),
    all_boundary_scales_at_most_one = all(
      boundary$checks$boundary_scale_at_most_one
    ),
    all_insufficient_states_identical = all(
      boundary$checks$insufficient_state_identical[
        !is.na(boundary$checks$insufficient_state_identical)
      ]
    ),
    all_downsamples_proportional = all(
      boundary$checks$proportional_downsample_exact[
        !is.na(boundary$checks$proportional_downsample_exact)
      ]
    ),
    all_likelihood_passage_ids_unique = all(
      uniqueness$every_passage_id_enters_once
    ),
    all_scenarios_use_shared_parameters = all(
      shared_parameter_check$identical_to_shared_reference
    ),
    stringsAsFactors = FALSE
  )
  scenario_time_rows[[seed_name]] <- scenario_time
  uniqueness_rows[[seed_name]] <- uniqueness
  shared_parameter_rows[[seed_name]] <- shared_parameter_check
}

overview <- dplyr::bind_rows(overview_rows)
scenario_time <- dplyr::bind_rows(scenario_time_rows)
likelihood_uniqueness <- dplyr::bind_rows(uniqueness_rows)
shared_parameters <- dplyr::bind_rows(shared_parameter_rows)

expected_cumulative_days <- c(
  `2N-C` = 29,
  `4N-C` = 29,
  `2N-O1` = 122,
  `2N-O2` = 122,
  `4N-O1` = 125,
  `4N-O2` = 121
)
scenario_time$expected_cumulative_days <- unname(
  expected_cumulative_days[scenario_time$scenario_id]
)
scenario_time$cumulative_time_matches_expected <-
  scenario_time$cumulative_passage_days ==
  scenario_time$expected_cumulative_days &
  scenario_time$cumulative_endpoint_days ==
  scenario_time$expected_cumulative_days &
  scenario_time$cumulative_selected_days ==
  scenario_time$expected_cumulative_days

shared_parameter_summary <- stats::aggregate(
  identical_to_shared_reference ~ seed + parameter,
  data = shared_parameters,
  FUN = all
)
names(shared_parameter_summary)[
  names(shared_parameter_summary) == "identical_to_shared_reference"
] <- "all_scenarios_identical_to_shared_reference"

write_tsv(
  overview,
  file.path(temporary_results_root, "validation_overview.tsv")
)
write_tsv(
  scenario_time,
  file.path(temporary_results_root, "scenario_cumulative_time.tsv")
)
write_tsv(
  likelihood_uniqueness,
  file.path(temporary_results_root, "likelihood_uniqueness.tsv")
)
write_tsv(
  shared_parameter_summary,
  file.path(temporary_results_root, "shared_parameter_summary.tsv")
)

hpc_baseline <- paste0(
  "/share/lab_crd/lab_crd/taoli/Project/",
  "miningcloneid_soft_coupling/oxygen/results/",
  "fit_invitro_O2_buffering_500seed"
)
source_files <- c(
  shared_path,
  common_path,
  model_path,
  loader_path,
  file.path(
    workflow_root,
    "util",
    "o2_supply_demand_map_invitro_lineage_utils.R"
  ),
  file.path(
    workflow_root,
    "util",
    "o2_supply_demand_map_invitro_lineage_simulation_utils.R"
  ),
  file.path(
    workflow_root,
    "util",
    "o2_supply_demand_map_invitro_objective_utils.R"
  ),
  file.path(
    workflow_root,
    "util",
    "o2_supply_demand_map_invitro_summary_utils.R"
  ),
  parameter_table,
  flow_density_path,
  unname(seed_parameter_paths)
)
provenance <- data.frame(
  path = normalizePath(source_files, mustWork = TRUE),
  md5 = unname(tools::md5sum(source_files)),
  stringsAsFactors = FALSE
)
write_tsv(
  provenance,
  file.path(temporary_results_root, "input_and_code_provenance.tsv")
)

run_context <- data.frame(
  field = c(
    "validation_type",
    "soft_coupling_root",
    "git_branch",
    "git_head",
    "hpc_baseline_path",
    "hpc_baseline_accessible",
    "result_policy",
    "dt",
    "init_total_size",
    "fixed_oxygen"
  ),
  value = c(
    "parameter replay only; no optimization",
    soft_coupling_root,
    system2(
      "git",
      c("-C", soft_coupling_root, "branch", "--show-current"),
      stdout = TRUE
    ),
    system2(
      "git",
      c("-C", soft_coupling_root, "rev-parse", "HEAD"),
      stdout = TRUE
    ),
    hpc_baseline,
    as.character(dir.exists(hpc_baseline)),
    "new directory; existing results never overwritten",
    as.character(cfg$DT),
    as.character(cfg$init_total_size),
    as.character(!isTRUE(cfg$o2_burden_feedback))
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  run_context,
  file.path(temporary_results_root, "run_context.tsv")
)

session_lines <- sub(
  "[[:blank:]]+$",
  "",
  capture.output(utils::sessionInfo())
)
writeLines(
  session_lines,
  con = file.path(temporary_results_root, "session_info.txt")
)

if (!all(scenario_time$cumulative_time_matches_expected)) {
  stop("Scenario cumulative-time validation failed.", call. = FALSE)
}
if (!all(overview$all_boundary_scales_at_most_one)) {
  stop("At least one boundary scale exceeds one.", call. = FALSE)
}
if (!all(overview$all_insufficient_states_identical)) {
  stop("At least one insufficient boundary changed its full state vector.", call. = FALSE)
}
if (!all(overview$all_downsamples_proportional)) {
  stop("At least one sufficient boundary was not proportionally downsampled.", call. = FALSE)
}
if (!all(overview$all_likelihood_passage_ids_unique)) {
  stop("At least one likelihood passage ID entered more than once.", call. = FALSE)
}
if (!all(overview$all_scenarios_use_shared_parameters)) {
  stop("At least one scenario did not use the shared parameter vector.", call. = FALSE)
}

if (!file.rename(temporary_results_root, results_root)) {
  stop(
    "Validation succeeded but the temporary result directory could not be moved to: ",
    results_root,
    call. = FALSE
  )
}

message("Validation results written to: ", results_root)
