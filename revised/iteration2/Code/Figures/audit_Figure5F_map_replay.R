#!/usr/bin/env Rscript

# Replay the saved joint-fit MAP endpoints used as Figure 5F chain starts.
#
# Scientific guardrail:
#   This script only checks whether commit 83953a8 reproduces the objective
#   values saved with the three selected C01/C02/C03 joint fits. It does not
#   generate posterior draws and does not treat optimizer endpoints as draws.
#
# Invocation (important because the historical backend resolves its own path
# from commandArgs before source-frame paths): first preload the backend from
# an Rscript -e top-level expression, then source this audit in the same R
# process. The exact copy-paste command is printed when it is not preloaded.

find_current_source <- function() {
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
  if (!length(frame_files)) {
    stop(
      "Source this audit with Rscript -e so its location and the historical ",
      "backend location can be resolved exactly."
    )
  }
  frame_files[[length(frame_files)]]
}

script_path <- find_current_source()
script_dir <- dirname(script_path)
iteration_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
snapshot_root <- file.path(iteration_root, "tmp", "miningcloneid_83953a8")
snapshot_sha <- "83953a874401e42cd176432786f889a896adc959"
joint_backend_path <- file.path(
  snapshot_root,
  "oxygen", "code", "O2_supply_demand_MAP", "util",
  "o2_supply_demand_map_fit_joint_backend.R"
)
result_root <- paste0(
  "/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/",
  "fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540"
)
selection_path <- file.path(
  iteration_root,
  "data", "Figures", "Figure5", "figure5f_selected_pair_inputs.tsv"
)
output_path <- file.path(
  iteration_root,
  "data", "Figures", "Figure5", "figure5f_map_replay_audit.tsv"
)

required_snapshot_inputs <- c(
  joint_backend_path,
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
  file.path(snapshot_root, "oxygen", "data", "g0g1_ploidy_density_grid.csv"),
  selection_path
)
missing <- required_snapshot_inputs[!file.exists(required_snapshot_inputs)]
if (length(missing)) {
  stop("Missing Figure 5F MAP-replay input(s):\n", paste(missing, collapse = "\n"))
}

selection <- utils::read.delim(
  selection_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
selection <- selection[selection$selected_for_figure5f %in% TRUE, , drop = FALSE]
selection <- selection[order(match(selection$family, c("C01", "C02", "C03"))), , drop = FALSE]
if (nrow(selection) != 3L || !identical(selection$family, c("C01", "C02", "C03"))) {
  stop("Expected exactly one selected Figure 5F input for C01, C02, and C03.")
}

# Source once so all three replays use the exact same compiled model backend.
# The backend must be preloaded at Rscript -e top level. If it were sourced from
# inside this file, the historical commandArgs-first bootstrap would mistake
# this audit's directory for the workflow util directory.
if (!exists("joint_env", envir = .GlobalEnv, inherits = FALSE) ||
    !is.environment(get("joint_env", envir = .GlobalEnv, inherits = FALSE)) ||
    !exists(
      "build_joint_context",
      envir = get("joint_env", envir = .GlobalEnv, inherits = FALSE),
      inherits = FALSE
    )) {
  quoted_backend <- encodeString(joint_backend_path, quote = '"')
  quoted_audit <- encodeString(script_path, quote = '"')
  stop(
    "Historical backend is not preloaded. Run exactly:\n",
    "Rscript -e 'joint_env <- new.env(parent=globalenv()); ",
    "sys.source(", quoted_backend, ", envir=joint_env, chdir=TRUE); ",
    "source(", quoted_audit, ")'"
  )
}
joint_env <- get("joint_env", envir = .GlobalEnv, inherits = FALSE)

read_effective_args <- function(path) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  tab <- tab[tab$source == "fit_command", c("key", "value"), drop = FALSE]
  if (!nrow(tab) || anyDuplicated(tab$key)) {
    stop("Malformed run_effective_args.tsv: ", path)
  }
  out <- as.list(tab$value)
  names(out) <- tab$key
  out
}

read_component <- function(path, name) {
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  hit <- suppressWarnings(as.numeric(tab$value[tab$component == name]))
  if (length(hit) != 1L || !is.finite(hit)) {
    stop("Missing/non-finite saved component '", name, "' in ", path)
  }
  hit
}

replay_one <- function(row) {
  run_name <- paste0("fit_joint_", row$warmup_label[[1]])
  seed_name <- as.character(row$selected_seed[[1]])
  run_dir <- file.path(result_root, run_name, seed_name)
  needed <- file.path(
    run_dir,
    c(
      "run_effective_args.tsv",
      "run_provenance.tsv",
      "best_params_transformed.tsv",
      "joint_components.tsv",
      "parameter_table_input_invivo.csv",
      "parameter_table_input_invitro.csv",
      "joint_soft_coupling_parameters_table_input.csv"
    )
  )
  if (any(!file.exists(needed))) {
    stop("Missing selected joint-fit file(s) under: ", run_dir)
  }

  argv <- read_effective_args(file.path(run_dir, "run_effective_args.tsv"))
  argv$config <- file.path(snapshot_root, "oxygen", "config", "O2_supply_demand.yaml")
  argv$data_dir <- file.path(snapshot_root, "data", "InVivoData_Gemcitabine")
  argv$parameter_table <- file.path(run_dir, "parameter_table_input_invivo.csv")
  argv$invitro_parameter_table <- file.path(run_dir, "parameter_table_input_invitro.csv")
  argv$fit_objects_dir <- file.path(
    snapshot_root, "oxygen", "ploidyOxygen", "data", "fit_objects"
  )
  argv$flow_density_path <- file.path(
    snapshot_root, "oxygen", "data", "g0g1_ploidy_density_grid.csv"
  )
  argv$necrosis_mapping_csv <- file.path(
    snapshot_root, "data", "InVivoData_Gemcitabine",
    "histology_to_dt_Gem_VT_20260209_v5_mapping.csv"
  )
  argv$joint_soft_coupling_parameters_table <- file.path(
    run_dir, "joint_soft_coupling_parameters_table_input.csv"
  )
  # Warm-up changes optimizer initialization only, never the objective target.
  # Disabling it avoids any dependency on separate-fit result directories.
  argv$joint_warmup_enable <- "FALSE"
  argv$n_cores <- "1"
  argv$joint_n_cores <- "1"
  argv$predict_n_cores <- "1"
  argv$deoptim_parallel <- "FALSE"
  argv$out_dir <- file.path(iteration_root, "tmp", "figure5f_map_replay", row$family[[1]])

  ctx <- joint_env$build_joint_context(argv)
  saved_par <- utils::read.delim(
    file.path(run_dir, "best_params_transformed.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expected_cols <- c("transformed_parameter", "transformed_value")
  if (!all(expected_cols %in% names(saved_par))) {
    stop("Malformed best_params_transformed.tsv in ", run_dir)
  }
  par_t <- setNames(
    as.numeric(saved_par$transformed_value),
    as.character(saved_par$transformed_parameter)
  )
  if (!setequal(names(par_t), names(ctx$init))) {
    stop(
      "Saved optimizer parameter names differ from commit 83953a8 context for ",
      row$family[[1]], ". Missing from saved={",
      paste(setdiff(names(ctx$init), names(par_t)), collapse = ","),
      "}; extra={", paste(setdiff(names(par_t), names(ctx$init)), collapse = ","), "}."
    )
  }
  par_t <- par_t[names(ctx$init)]
  if (any(!is.finite(par_t)) || any(par_t < ctx$lower - 1e-10) || any(par_t > ctx$upper + 1e-10)) {
    stop("Saved MAP endpoint is non-finite or outside commit 83953a8 bounds for ", row$family[[1]])
  }

  replay <- joint_env$joint_objective_components(par_t, ctx)
  component_path <- file.path(run_dir, "joint_components.tsv")
  values <- data.frame(
    component = c(
      "objective_joint",
      "objective_joint_unpenalized",
      "objective_soft_coupling",
      "objective_invivo",
      "objective_invitro"
    ),
    replayed_value = c(
      replay$objective,
      replay$objective_unpenalized,
      replay$objective_soft_coupling,
      replay$invivo$L,
      replay$invitro$objective
    ),
    stringsAsFactors = FALSE
  )
  values$saved_value <- vapply(
    values$component,
    function(x) read_component(component_path, x),
    numeric(1)
  )
  values$absolute_difference <- abs(values$replayed_value - values$saved_value)
  values$tolerance <- 1e-8
  values$passed <- is.finite(values$absolute_difference) &
    values$absolute_difference <= values$tolerance

  provenance <- utils::read.delim(
    file.path(run_dir, "run_provenance.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  git_commit <- provenance$value[
    provenance$section == "git" & provenance$key == "commit"
  ]
  git_status <- provenance$value[
    provenance$section == "git" & provenance$key == "dirty_status"
  ]
  if (!length(git_commit)) git_commit <- NA_character_
  if (!length(git_status)) git_status <- NA_character_

  data.frame(
    family = row$family[[1]],
    warmup_label = row$warmup_label[[1]],
    selected_seed = seed_name,
    selection_metric = row$selection_metric[[1]],
    separate_invivo_objective = row$separate_invivo_objective[[1]],
    replay_code_commit = snapshot_sha,
    saved_run_git_commit = git_commit[[1]],
    saved_run_git_status = git_status[[1]],
    component_source = normalizePath(component_path, mustWork = TRUE),
    transformed_parameter_source = normalizePath(
      file.path(run_dir, "best_params_transformed.tsv"), mustWork = TRUE
    ),
    values,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

audit <- do.call(
  rbind,
  lapply(seq_len(nrow(selection)), function(i) replay_one(selection[i, , drop = FALSE]))
)
rownames(audit) <- NULL

output_norm <- normalizePath(output_path, mustWork = FALSE)
allowed_root <- normalizePath(file.path(iteration_root, "data", "Figures", "Figure5"), mustWork = TRUE)
if (!startsWith(output_norm, paste0(allowed_root, .Platform$file.sep))) {
  stop("Refusing to write outside iteration2 Figure5 data: ", output_norm)
}
utils::write.table(
  audit,
  output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

print(audit[, c(
  "family", "component", "saved_value", "replayed_value",
  "absolute_difference", "passed"
)])
if (!all(audit$passed)) {
  stop(
    "Commit 83953a8 did not reproduce every saved Figure 5F MAP component ",
    "within 1e-8. Posterior sampling must not proceed from this snapshot."
  )
}
message("MAP replay passed for C01, C02, and C03: ", output_path)
