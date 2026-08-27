Sys.setenv(
  OMP_NUM_THREADS = "1",
  KMP_USE_SHM = "0",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

resolve_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  }
  frames <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(frame) {
        path <- frame$ofile
        if (is.null(path)) "" else normalizePath(path, mustWork = FALSE)
      },
      character(1L)
    )
  )
  if (length(frames)) return(frames[[length(frames)]])
  normalizePath(getwd(), mustWork = TRUE)
}

RUNTIME_PATH_OPTION_ENV <- c(
  "invitro-result-root" = "FIGURE_INVITRO_RESULT_ROOT",
  "invivo-result-root" = "FIGURE_INVIVO_RESULT_ROOT",
  "joint-result-root" = "FIGURE_JOINT_RESULT_ROOT",
  "gemcitabine-data-root" = "FIGURE_GEMCITABINE_DATA_ROOT",
  "ltee-data-root" = "FIGURE_LTEE_DATA_ROOT"
)

parse_runtime_path_options <- function(
    args = commandArgs(trailingOnly = TRUE)
) {
  values <- stats::setNames(
    rep(NA_character_, length(RUNTIME_PATH_OPTION_ENV)),
    names(RUNTIME_PATH_OPTION_ENV)
  )
  for (option in names(RUNTIME_PATH_OPTION_ENV)) {
    prefix <- paste0("--", option, "=")
    matched <- args[startsWith(args, prefix)]
    if (length(matched) > 1L) {
      stop("Runtime path option was provided more than once: --", option)
    }
    if (length(matched) == 1L) {
      value <- substring(matched[[1L]], nchar(prefix) + 1L)
      if (!nzchar(trimws(value))) {
        stop("Runtime path option cannot be empty: --", option)
      }
      values[[option]] <- value
    }
  }
  values
}

RUNTIME_PATH_OPTION_VALUES <- parse_runtime_path_options()

resolve_runtime_input_path <- function(option) {
  value <- RUNTIME_PATH_OPTION_VALUES[[option]]
  if (is.na(value) || !nzchar(trimws(value))) {
    value <- Sys.getenv(RUNTIME_PATH_OPTION_ENV[[option]], unset = "")
  }
  if (!nzchar(trimws(value))) {
    stop(
      "Missing required runtime path option: --", option,
      "=PATH"
    )
  }
  if (!dir.exists(value)) {
    stop(
      "Runtime path option is not an existing directory: --",
      option, "=", value
    )
  }
  normalizePath(value, mustWork = TRUE)
}

is_runtime_path_argument <- function(argument) {
  any(vapply(
    names(RUNTIME_PATH_OPTION_ENV),
    function(option) startsWith(argument, paste0("--", option, "=")),
    logical(1L)
  ))
}

find_workspace_root <- function(start = dirname(resolve_script_path())) {
  explicit <- trimws(Sys.getenv("FIGURE_WORKSPACE_ROOT", unset = ""))
  if (nzchar(explicit)) {
    explicit <- normalizePath(explicit, mustWork = TRUE)
    if (file.exists(file.path(explicit, "manager.sh")) &&
        dir.exists(file.path(explicit, "Code", "Figures"))) {
      return(explicit)
    }
    stop(
      "FIGURE_WORKSPACE_ROOT does not contain manager.sh and Code/Figures: ",
      explicit
    )
  }
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "manager.sh")) &&
        dir.exists(file.path(current, "Code", "Figures"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the figure workspace from: ", start)
    }
    current <- parent
  }
}

WORKSPACE_ROOT <- find_workspace_root()
# iteration4 owns every generated output and all figure-specific code. The
# mechanistic model is the one explicit read-only external dependency required
# by the analysis contract below.
REPO_ROOT <- WORKSPACE_ROOT
INVITRO_RESULT_ROOT <- resolve_runtime_input_path("invitro-result-root")
INVIVO_RESULT_ROOT <- resolve_runtime_input_path("invivo-result-root")
JOINT_RESULT_ROOT <- resolve_runtime_input_path("joint-result-root")
GEMCITABINE_DATA_ROOT <- resolve_runtime_input_path(
  "gemcitabine-data-root"
)
LTEE_DATA_ROOT <- resolve_runtime_input_path("ltee-data-root")
DEFAULT_MODEL_CODE_ROOT <-
  "/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP"
MODEL_CODE_ROOT <- trimws(Sys.getenv(
  "FIGURE_MODEL_CODE_ROOT", unset = DEFAULT_MODEL_CODE_ROOT
))
if (!nzchar(MODEL_CODE_ROOT)) MODEL_CODE_ROOT <- DEFAULT_MODEL_CODE_ROOT
MODEL_CODE_ROOT <- normalizePath(MODEL_CODE_ROOT, mustWork = TRUE)

do.call(
  Sys.setenv,
  as.list(c(
    HYPOXIA_REPO_ROOT = REPO_ROOT,
    FIGURE_MODEL_CODE_ROOT = MODEL_CODE_ROOT,
    FIGURE_INVITRO_RESULT_ROOT = INVITRO_RESULT_ROOT,
    FIGURE_INVIVO_RESULT_ROOT = INVIVO_RESULT_ROOT,
    FIGURE_JOINT_RESULT_ROOT = JOINT_RESULT_ROOT,
    FIGURE_GEMCITABINE_DATA_ROOT = GEMCITABINE_DATA_ROOT,
    FIGURE_LTEE_DATA_ROOT = LTEE_DATA_ROOT
  ))
)

runtime_path_arguments <- function() {
  paste0(
    "--", names(RUNTIME_PATH_OPTION_ENV), "=",
    c(
      INVITRO_RESULT_ROOT,
      INVIVO_RESULT_ROOT,
      JOINT_RESULT_ROOT,
      GEMCITABINE_DATA_ROOT,
      LTEE_DATA_ROOT
    )
  )
}

CODE_CONTAINER_ROOT <- file.path(WORKSPACE_ROOT, "Code")
CODE_ROOT <- file.path(WORKSPACE_ROOT, "Code", "Figures")
CONFIG_ROOT <- file.path(CODE_CONTAINER_ROOT, "config")
TEMPLATE_ROOT <- file.path(CODE_CONTAINER_ROOT, "templates")
DATA_ROOT <- file.path(WORKSPACE_ROOT, "data", "Figures")
OUTPUT_ROOT <- file.path(WORKSPACE_ROOT, "Figures")
MANUSCRIPT_ROOT <- file.path(WORKSPACE_ROOT, "manuscript")
PUBLISHED_OUTPUT_ROOT <- file.path(MANUSCRIPT_ROOT, "Figures")
AUDIT_ROOT <- file.path(WORKSPACE_ROOT, "audit")

if (!dir.exists(MODEL_CODE_ROOT)) {
  stop("Missing required external model code: ", MODEL_CODE_ROOT)
}
required_model_files <- file.path(
  MODEL_CODE_ROOT,
  c(
    "model/model_O2_supply_demand_MAP.R",
    "model/model_O2_supply_demand_MAP.cpp",
    "util/o2_supply_demand_map_shared.R",
    "util/o2_supply_demand_map_common_semantics.R",
    "simulation/o2/fixed_o2/run_fixed_o2_simulation.R"
  )
)
if (any(!file.exists(required_model_files))) {
  stop(
    "Required external model implementation is incomplete:\n",
    paste(required_model_files[!file.exists(required_model_files)], collapse = "\n")
  )
}

INVIVO_VISUALIZATION_SEED <- "seed25"
INVITRO_VISUALIZATION_SEED <- "seed228"
JOINT_FAMILY_LEVELS <- sprintf("C%02d", seq_len(6L))

joint_family_from_label <- function(x) {
  out <- regmatches(as.character(x), regexpr("C[0-9]{2}", as.character(x)))
  out[!nzchar(out)] <- NA_character_
  out
}

read_joint_warmup_manifest <- function() {
  path <- file.path(JOINT_RESULT_ROOT, "multi_warmup_manifest.tsv")
  if (!file.exists(path)) stop("Missing joint warm-up manifest: ", path)
  manifest <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE
  )
  required <- c(
    "warmup_label", "invivo_seed", "invitro_seed", "joint_run_prefix"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("Joint warm-up manifest is missing: ", paste(missing, collapse = ", "))
  }
  manifest$family <- joint_family_from_label(manifest$warmup_label)
  manifest$invivo_seed <- paste0("seed", as.integer(manifest$invivo_seed))
  manifest$invitro_seed <- paste0("seed", as.integer(manifest$invitro_seed))
  manifest <- manifest[
    order(match(manifest$family, JOINT_FAMILY_LEVELS)), , drop = FALSE
  ]
  if (nrow(manifest) != 6L ||
      !identical(manifest$family, JOINT_FAMILY_LEVELS) ||
      any(manifest$invitro_seed != INVITRO_VISUALIZATION_SEED) ||
      anyDuplicated(manifest$warmup_label)) {
    stop(
      "Joint result must contain exactly one C01-C06 primary-cluster pair ",
      "anchored to ", INVITRO_VISUALIZATION_SEED, "."
    )
  }
  manifest
}

copy_directory_tree <- function(source, destination) {
  source <- normalizePath(source, mustWork = TRUE)
  if (dir.exists(destination)) return(invisible(FALSE))
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  source_dirs <- list.dirs(source, recursive = TRUE, full.names = TRUE)
  for (directory in source_dirs) {
    relative <- substring(directory, nchar(source) + 2L)
    if (nzchar(relative)) {
      dir.create(
        file.path(destination, relative),
        recursive = TRUE,
        showWarnings = FALSE
      )
    }
  }
  source_files <- list.files(
    source,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE,
    include.dirs = FALSE
  )
  for (source_file in source_files) {
    relative <- substring(source_file, nchar(source) + 2L)
    destination_file <- file.path(destination, relative)
    dir.create(
      dirname(destination_file),
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!file.copy(
      source_file,
      destination_file,
      overwrite = FALSE,
      copy.mode = TRUE,
      copy.date = TRUE
    )) {
      stop("Failed to bootstrap manuscript file: ", source_file)
    }
  }
  invisible(TRUE)
}

publish_bootstrap_file <- function(source, destination) {
  if (!file.exists(source)) stop("Missing bootstrap configuration: ", source)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(
    source, destination, overwrite = TRUE, copy.mode = TRUE
  )) {
    stop("Failed to publish bootstrap configuration: ", source)
  }
  invisible(destination)
}

write_resolved_input_manifest <- function() {
  source <- file.path(
    CONFIG_ROOT, "manifests", "allowed_scientific_inputs.txt"
  )
  lines <- readLines(source, warn = FALSE)
  replacements <- c(
    "${INVITRO_RESULT_ROOT}" = INVITRO_RESULT_ROOT,
    "${INVIVO_RESULT_ROOT}" = INVIVO_RESULT_ROOT,
    "${JOINT_RESULT_ROOT}" = JOINT_RESULT_ROOT,
    "${GEMCITABINE_DATA_ROOT}" = GEMCITABINE_DATA_ROOT,
    "${LTEE_DATA_ROOT}" = LTEE_DATA_ROOT
  )
  for (token in names(replacements)) {
    lines <- gsub(token, replacements[[token]], lines, fixed = TRUE)
  }
  destination <- file.path(
    AUDIT_ROOT, "manifests", "allowed_scientific_inputs.txt"
  )
  writeLines(lines, destination, useBytes = TRUE)
  invisible(destination)
}

write_resolved_pipeline_parameters <- function() {
  parameters <- data.frame(
    parameter = c(
      "workspace_root", "repository_root",
      "baseline_logic", "invitro_result_root", "invivo_result_root",
      "joint_result_root", "gemcitabine_data_root", "ltee_data_root",
      "primary_frozen_input_source", "model_code_root",
      "invivo_visualization_seed", "invitro_visualization_seed",
      "fixed_o2_min", "fixed_o2_max",
      "fixed_o2_n", "si2_class_lower", "si2_class_upper",
      "si2_boundary_rule"
    ),
    value = c(
      WORKSPACE_ROOT, REPO_ROOT,
      "current packaged figure logic", INVITRO_RESULT_ROOT,
      INVIVO_RESULT_ROOT, JOINT_RESULT_ROOT, GEMCITABINE_DATA_ROOT,
      LTEE_DATA_ROOT, LTEE_DATA_ROOT,
      MODEL_CODE_ROOT, INVIVO_VISUALIZATION_SEED,
      INVITRO_VISUALIZATION_SEED, "0", "5", "201",
      "0.8", "1.2", "outer_inclusive"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    parameters,
    file.path(AUDIT_ROOT, "parameters", "pipeline_parameters.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  invisible(parameters)
}

bootstrap_workspace <- function() {
  required_config <- c(
    file.path(CONFIG_ROOT, "workspace_layout.tsv"),
    file.path(CONFIG_ROOT, "manifests", c(
      "allowed_scientific_inputs.txt",
      "expected_scientific_input_md5.tsv",
      "figure_entrypoints.tsv",
      "figure_logic_provenance.tsv"
    )),
    file.path(CONFIG_ROOT, "parameters", c(
      "parameter_function_groups.tsv",
      "parameter_function_group_palette.tsv"
    )),
    file.path(TEMPLATE_ROOT, "README.md")
  )
  missing_config <- required_config[!file.exists(required_config)]
  if (length(missing_config)) {
    stop(
      "Portable Code bundle is incomplete: ",
      paste(missing_config, collapse = ", ")
    )
  }

  dir.create(
    file.path(AUDIT_ROOT, "md5"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  if (!dir.exists(MANUSCRIPT_ROOT)) {
    stop("iteration4 manuscript directory is missing: ", MANUSCRIPT_ROOT)
  }

  layout <- utils::read.delim(
    file.path(CONFIG_ROOT, "workspace_layout.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!identical(names(layout), c("path", "type")) ||
      any(layout$type != "directory")) {
    stop("Invalid Code/config/workspace_layout.tsv.")
  }
  invisible(lapply(
    file.path(WORKSPACE_ROOT, layout$path),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  publish_bootstrap_file(
    file.path(CONFIG_ROOT, "manifests", "figure_entrypoints.tsv"),
    file.path(AUDIT_ROOT, "manifests", "figure_entrypoints.tsv")
  )
  publish_bootstrap_file(
    file.path(CONFIG_ROOT, "manifests", "figure_logic_provenance.tsv"),
    file.path(AUDIT_ROOT, "manifests", "figure_logic_provenance.tsv")
  )
  write_resolved_input_manifest()
  for (name in c(
    "parameter_function_groups.tsv",
    "parameter_function_group_palette.tsv"
  )) {
    publish_bootstrap_file(
      file.path(CONFIG_ROOT, "parameters", name),
      file.path(AUDIT_ROOT, "parameters", name)
    )
  }
  write_resolved_pipeline_parameters()

  entrypoints <- utils::read.delim(
    file.path(AUDIT_ROOT, "manifests", "figure_entrypoints.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  invisible(lapply(
    file.path(WORKSPACE_ROOT, entrypoints$data_directory),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  readme <- file.path(WORKSPACE_ROOT, "README.md")
  if (!file.exists(readme)) {
    publish_bootstrap_file(file.path(TEMPLATE_ROOT, "README.md"), readme)
  }
  required_manuscript <- file.path(
    MANUSCRIPT_ROOT,
    c("ltee_hypoxia_model.tex", "references_Zotero_TaoLi.bib")
  )
  if (any(!file.exists(required_manuscript))) {
    stop(
      "Workspace manuscript exists but is incomplete; missing: ",
      paste(required_manuscript[!file.exists(required_manuscript)], collapse = ", ")
    )
  }
  invisible(WORKSPACE_ROOT)
}

bootstrap_workspace()

ALLOWED_INPUT_MANIFEST <- file.path(
  AUDIT_ROOT, "manifests", "allowed_scientific_inputs.txt"
)
ALLOWED_SCIENTIFIC_INPUT_ROOTS <- normalizePath(
  Filter(
    nzchar,
    trimws(readLines(ALLOWED_INPUT_MANIFEST, warn = FALSE))
  ),
  mustWork = TRUE
)

ensure_workspace_directories <- function() {
  paths <- c(
    DATA_ROOT,
    OUTPUT_ROOT,
    PUBLISHED_OUTPUT_ROOT,
    file.path(AUDIT_ROOT, c(
      "md5", "logs", "parameters", "manifests", "reports"
    ))
  )
  entrypoint_manifest <- file.path(
    AUDIT_ROOT, "manifests", "figure_entrypoints.tsv"
  )
  if (file.exists(entrypoint_manifest)) {
    entrypoints <- utils::read.delim(
      entrypoint_manifest,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if ("data_directory" %in% names(entrypoints)) {
      paths <- c(
        paths,
        file.path(WORKSPACE_ROOT, entrypoints$data_directory)
      )
    }
  }
  invisible(lapply(
    paths, dir.create, recursive = TRUE, showWarnings = FALSE
  ))
}

assert_allowed_input <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  allowed <- vapply(
    ALLOWED_SCIENTIFIC_INPUT_ROOTS,
    function(root) {
      identical(normalized, root) ||
        startsWith(normalized, paste0(root, .Platform$file.sep))
    },
    logical(1L)
  )
  if (!any(allowed)) {
    stop("Scientific input is outside the approved roots: ", normalized)
  }
  normalized
}

assert_local_intermediate <- function(path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(DATA_ROOT, mustWork = TRUE)
  if (!identical(normalized, root) &&
      !startsWith(normalized, paste0(root, .Platform$file.sep))) {
    stop("Intermediate path is outside the workspace data root: ", normalized)
  }
  normalized
}

copy_input <- function(source, destination, overwrite = TRUE) {
  source <- assert_allowed_input(source)
  destination <- assert_local_intermediate(destination)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(
    source, destination, overwrite = overwrite, copy.mode = TRUE
  )
  if (!ok) stop("Failed to copy input: ", source, " -> ", destination)
  normalizePath(destination, mustWork = TRUE)
}

write_intermediate_tsv <- function(x, path) {
  path <- assert_local_intermediate(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  normalizePath(path, mustWork = TRUE)
}

read_metric_map <- function(path) {
  tab <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!all(c("metric", "value") %in% names(tab))) {
    stop("Expected metric/value table: ", path)
  }
  tab <- tab[!duplicated(tab$metric), , drop = FALSE]
  stats::setNames(as.character(tab$value), tab$metric)
}

as_boolean <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  tolower(trimws(as.character(x[[1L]]))) %in%
    c("true", "t", "1", "yes", "y", "on")
}

run_process <- function(command, args = character(), env = character()) {
  status <- system2(command, args = args, env = env)
  if (!identical(status, 0L)) {
    stop(
      "Command failed with status ", status, ": ",
      paste(c(command, args), collapse = " ")
    )
  }
  invisible(status)
}

ensure_workspace_directories()
