# Loader for the staged fixed-O2 pipeline and legacy compatibility entrypoints.

.fixo2_loader_dir <- local({
  if (exists(".fixo2_loader_path", inherits = TRUE)) {
    loader_path <- get(".fixo2_loader_path", inherits = TRUE)
    return(dirname(normalizePath(loader_path, mustWork = TRUE)))
  }
  frames <- sys.frames()
  frame_files <- Filter(
    nzchar,
    vapply(frames, function(env) {
      path <- env$ofile
      if (is.null(path)) "" else normalizePath(path, mustWork = FALSE)
    }, character(1))
  )
  own <- frame_files[
    basename(frame_files) == "fixed_o2_pipeline_loader.R"
  ]
  if (length(own)) {
    dirname(own[[length(own)]])
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) {
      dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
    } else {
      getwd()
    }
  }
})

FIXO2_WORKFLOW_ROOT <- normalizePath(
  file.path(.fixo2_loader_dir, "..", ".."),
  mustWork = TRUE
)
FIXO2_REPO_ROOT <- normalizePath(
  file.path(FIXO2_WORKFLOW_ROOT, "..", "..", ".."),
  mustWork = TRUE
)
# Preserve the historical process-helper path semantics for extracted functions.
SCRIPT_DIR <- normalizePath(
  file.path(FIXO2_WORKFLOW_ROOT, "analysis", "fixed_o2"),
  mustWork = TRUE
)

source(
  file.path(
    FIXO2_WORKFLOW_ROOT,
    "analysis",
    "process_fingerprints",
    "process_fingerprint_utils.R"
  ),
  local = environment()
)
source(
  file.path(
    FIXO2_WORKFLOW_ROOT,
    "analysis",
    "process_fingerprints",
    "ploidy_regime_utils.R"
  ),
  local = environment()
)

FIXO2_SIMULATION_SCRIPT <- normalizePath(
  file.path(
    FIXO2_WORKFLOW_ROOT,
    "simulation",
    "o2",
    "fixed_o2",
    "run_fixed_o2_simulation.R"
  ),
  mustWork = TRUE
)
sys.source(
  file.path(
    FIXO2_WORKFLOW_ROOT,
    "simulation",
    "o2",
    "fixed_o2",
    "fixed_o2_simulation_utils.R"
  ),
  envir = environment(),
  chdir = TRUE
)
sys.source(
  file.path(
    FIXO2_WORKFLOW_ROOT,
    "util",
    "o2_supply_demand_map_fixed_o2_format_utils.R"
  ),
  envir = environment(),
  chdir = TRUE
)

# The canonical simulator is a CLI and deliberately sources the model at its
# top level.  Compatibility analyses only need a small set of its pure helper
# functions.  Load function definitions from its parse tree without executing
# any CLI setup, model source, or main().
.fixo2_runner_function_env <- new.env(parent = environment())
.fixo2_runner_exprs <- parse(FIXO2_SIMULATION_SCRIPT, keep.source = FALSE)
for (.fixo2_runner_expr in .fixo2_runner_exprs) {
  if (
    is.call(.fixo2_runner_expr) &&
      identical(.fixo2_runner_expr[[1L]], as.name("<-")) &&
      is.symbol(.fixo2_runner_expr[[2L]]) &&
      is.call(.fixo2_runner_expr[[3L]]) &&
      identical(.fixo2_runner_expr[[3L]][[1L]], as.name("function"))
  ) {
    eval(.fixo2_runner_expr, envir = .fixo2_runner_function_env)
  }
}
.fixo2_required_runner_functions <- c(
  "fixo2_fixed_matrix",
  "fixo2_init_vector",
  "fixo2_normalize_state",
  "fixo2_eigen_states",
  "fixo2_eigen_trajectory",
  "fixo2_expm_states",
  "fixo2_trajectory_with_fallback",
  "fixo2_dominant_from_eig",
  "fixo2_dominant_one",
  "fixo2_simulation_output_complete",
  "fixo2_filter_missing_simulation_tasks"
)
.fixo2_missing_runner_functions <- .fixo2_required_runner_functions[
  !vapply(
    .fixo2_required_runner_functions,
    exists,
    logical(1),
    envir = .fixo2_runner_function_env,
    inherits = FALSE
  )
]
if (length(.fixo2_missing_runner_functions)) {
  stop(
    "Canonical fixed-O2 simulator is missing required helper(s): ",
    paste(.fixo2_missing_runner_functions, collapse = ", "),
    call. = FALSE
  )
}
for (.fixo2_runner_name in .fixo2_required_runner_functions) {
  assign(
    .fixo2_runner_name,
    get(
      .fixo2_runner_name,
      envir = .fixo2_runner_function_env,
      inherits = FALSE
    ),
    envir = environment()
  )
}
rm(
  .fixo2_runner_expr,
  .fixo2_runner_exprs,
  .fixo2_runner_function_env,
  .fixo2_runner_name,
  .fixo2_required_runner_functions,
  .fixo2_missing_runner_functions
)

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required")
  }
})

`%||%` <- o2ipa_null_coalesce

.fixo2_module_paths <- c(
  file.path(
    FIXO2_WORKFLOW_ROOT,
    "util",
    "o2_supply_demand_map_fixed_o2_validation_utils.R"
  ),
  file.path(
    FIXO2_WORKFLOW_ROOT,
    "util",
    "o2_supply_demand_map_fixed_o2_table_utils.R"
  ),
  file.path(
    FIXO2_WORKFLOW_ROOT,
    "simulation",
    "o2",
    "fixed_o2",
    "fixed_o2_numerical_producers.R"
  ),
  file.path(FIXO2_WORKFLOW_ROOT, "analysis", "fixed_o2", "fixed_o2_analysis_utils.R"),
  file.path(FIXO2_WORKFLOW_ROOT, "vis", "fixed_o2", "fixed_o2_plot_utils.R"),
  file.path(FIXO2_WORKFLOW_ROOT, "report", "fixed_o2", "fixed_o2_report_utils.R"),
  file.path(.fixo2_loader_dir, "fixed_o2_legacy_pipeline_functions.R")
)
.fixo2_missing_modules <- .fixo2_module_paths[!file.exists(.fixo2_module_paths)]
if (length(.fixo2_missing_modules)) {
  stop(
    "Missing fixed-O2 pipeline module(s): ",
    paste(.fixo2_missing_modules, collapse = ", "),
    call. = FALSE
  )
}
for (.fixo2_module_path in .fixo2_module_paths) {
  sys.source(.fixo2_module_path, envir = environment(), chdir = TRUE)
}
rm(.fixo2_module_path, .fixo2_module_paths, .fixo2_missing_modules)

fixo2_deprecation_message <- function(entrypoint) {
  paste0(
    entrypoint,
    " is a deprecated compatibility orchestrator. New workflows should run ",
    "simulation/o2/fixed_o2, analysis/fixed_o2, vis/fixed_o2, and report ",
    "entrypoints explicitly."
  )
}
