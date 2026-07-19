#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("oxygen/tests/check_o2_layer_boundaries.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
workflow_root <- file.path(repo_root, "oxygen", "code", "O2_supply_demand_MAP")

read_code <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  paste(lines, collapse = "\n")
}

is_compatibility_wrapper <- function(path) {
  header <- paste(head(readLines(path, warn = FALSE), 30L), collapse = "\n")
  grepl(
    paste(
      c(
        "Backward-compatible",
        "Compatibility pipeline wrapper",
        "compatibility orchestrator",
        "deprecated compatibility",
        "thin compatibility wrapper"
      ),
      collapse = "|"
    ),
    header,
    ignore.case = TRUE,
    perl = TRUE
  )
}

r_files <- function(dir) {
  if (!dir.exists(dir)) return(character())
  sort(list.files(dir, recursive = TRUE, full.names = TRUE, pattern = "\\.[Rr]$"))
}

assert_no_patterns <- function(paths, patterns, layer, allow_compat = FALSE) {
  failures <- character()
  checked <- 0L
  skipped <- 0L
  for (path in paths) {
    if (allow_compat && is_compatibility_wrapper(path)) {
      skipped <- skipped + 1L
      wrapper_lines <- length(readLines(path, warn = FALSE))
      if (wrapper_lines > 250L) {
        failures <- c(
          failures,
          paste0(
            layer,
            " compatibility file is not thin (",
            wrapper_lines,
            " lines): ",
            substring(path, nchar(repo_root) + 2L)
          )
        )
      }
      next
    }
    checked <- checked + 1L
    text <- read_code(path)
    for (label in names(patterns)) {
      if (grepl(patterns[[label]], text, perl = TRUE, ignore.case = TRUE)) {
        failures <- c(
          failures,
          paste0(
            layer,
            " boundary violation [",
            label,
            "]: ",
            substring(path, nchar(repo_root) + 2L)
          )
        )
      }
    }
  }
  if (length(failures)) {
    failures <- unique(failures)
    cat(paste(failures, collapse = "\n"), "\n", file = stderr(), sep = "")
    stop(
      length(failures),
      " ",
      layer,
      " boundary violation(s); full list printed above.",
      call. = FALSE
    )
  }
  c(checked = checked, compatibility_wrappers = skipped)
}

plot_patterns <- c(
  ggplot = "(^|[^A-Za-z0-9_.])(ggplot|ggsave)\\s*\\(",
  ggplot_namespace = "ggplot2::",
  base_plot = paste0(
    "(^|[^A-Za-z0-9_.])",
    "(plot|lines|points|hist|heatmap|image|boxplot|legend|axis|abline|",
    "polygon|segments|matplot)\\s*\\("
  ),
  graphics_device = "(^|[^A-Za-z0-9_.])((grDevices::)?(pdf|png|jpeg|svg))\\s*\\(",
  report_render = "(rmarkdown|quarto)::render\\s*\\("
)

simulation_paths <- r_files(file.path(workflow_root, "simulation"))
analysis_paths <- r_files(file.path(workflow_root, "analysis"))
vis_paths <- r_files(file.path(workflow_root, "vis"))
report_paths <- r_files(file.path(workflow_root, "report"))
util_paths <- r_files(file.path(workflow_root, "util"))

simulation_result <- assert_no_patterns(
  simulation_paths,
  c(
    plot_patterns,
    downstream_output_dir = "file\\.path\\s*\\([^\\n]*out_dir[^\\n]*(figures|report)",
    source_analysis = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)analysis",
    source_vis = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)vis",
    source_report = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)report"
  ),
  layer = "simulation",
  allow_compat = TRUE
)

analysis_result <- assert_no_patterns(
  analysis_paths,
  c(
    plot_patterns,
    downstream_output_dir = "file\\.path\\s*\\([^\\n]*out_dir[^\\n]*(figures|report)",
    source_model = "(source|sys\\.source)\\s*\\([^\\n]*model_O2_supply_demand_MAP",
    source_simulation = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)simulation",
    source_vis = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)vis",
    invoke_simulator = "(system2|system|source|sys\\.source)\\s*\\([^\\n]*(simulate_|generate_.*simulation|run_fixed_o2_simulation)"
  ),
  layer = "analysis",
  allow_compat = TRUE
)

vis_result <- assert_no_patterns(
  vis_paths,
  c(
    read_fit_object = "readRDS\\s*\\(",
    best_params = "best_params\\.(tsv|rds)",
    fit_result = "fit_result\\.(rds|Rds|RDS)",
    source_model = "(source|sys\\.source)\\s*\\([^\\n]*model_O2_supply_demand_MAP",
    source_optimizer = "(source|sys\\.source)\\s*\\([^\\n]*optimizer",
    source_simulation = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)simulation",
    source_analysis = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)analysis"
  ),
  layer = "vis"
)

report_result <- assert_no_patterns(
  report_paths,
  c(
    plot_patterns[names(plot_patterns) != "report_render"],
    read_fit_object = "readRDS\\s*\\(",
    source_model = "(source|sys\\.source)\\s*\\([^\\n]*model_O2_supply_demand_MAP",
    source_simulation = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)simulation",
    source_analysis = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)analysis",
    source_vis = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)vis"
  ),
  layer = "report"
)

util_result <- assert_no_patterns(
  util_paths,
  c(
    plot_patterns,
    source_simulation = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)simulation",
    source_analysis = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)analysis",
    source_vis = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)vis",
    source_report = "(source|sys\\.source)\\s*\\([^\\n]*(/|file\\.path\\([^\\n]*)report"
  ),
  layer = "util"
)

portable_r_paths <- c(
  simulation_paths,
  analysis_paths,
  vis_paths,
  report_paths,
  util_paths,
  r_files(file.path(workflow_root, "runner"))
)
hardcoded_macos_users <- portable_r_paths[vapply(
  portable_r_paths,
  function(path) {
    grepl(
      "/Users/[^/[:space:]\"']+/",
      read_code(path),
      perl = TRUE
    )
  },
  logical(1)
)]
if (length(hardcoded_macos_users)) {
  stop(
    "Hard-coded macOS user paths outside HPC/protected code:\n",
    paste(
      substring(hardcoded_macos_users, nchar(repo_root) + 2L),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

misplaced_shell <- list.files(
  workflow_root,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.(sh|sub)$"
)
allowed_shell_roots <- c(
  normalizePath(file.path(workflow_root, "hpc"), mustWork = TRUE),
  normalizePath(file.path(workflow_root, "runner"), mustWork = TRUE)
)
misplaced_shell <- misplaced_shell[
  !vapply(
    normalizePath(misplaced_shell, mustWork = TRUE),
    function(path) {
      identical(
        path,
        normalizePath(
          file.path(
            workflow_root,
            "util",
            "o2_supply_demand_map_shell_utils.sh"
          ),
          mustWork = TRUE
        )
      ) || any(startsWith(path, paste0(allowed_shell_roots, "/")))
    },
    logical(1)
  )
]
if (length(misplaced_shell)) {
  stop(
    "Shell/HPC entrypoints outside hpc/, runner/, or the canonical util shell library:\n",
    paste(misplaced_shell, collapse = "\n"),
    call. = FALSE
  )
}

misplaced_hpc_names <- list.files(
  workflow_root,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "hpc"
)
misplaced_hpc_names <- misplaced_hpc_names[file.info(misplaced_hpc_names)$isdir %in% FALSE]
misplaced_hpc_names <- misplaced_hpc_names[
  !startsWith(
    normalizePath(misplaced_hpc_names, mustWork = TRUE),
    paste0(normalizePath(file.path(workflow_root, "hpc"), mustWork = TRUE), "/")
  )
]
if (length(misplaced_hpc_names)) {
  stop(
    "HPC-named entrypoints outside hpc/:\n",
    paste(misplaced_hpc_names, collapse = "\n"),
    call. = FALSE
  )
}

tree_entries <- list.files(
  workflow_root,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  include.dirs = TRUE,
  no.. = TRUE
)
protected_artifact_roots <- normalizePath(
  file.path(workflow_root, c("model", "optimizer")),
  mustWork = TRUE
)
is_under_protected_root <- function(path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  any(vapply(
    protected_artifact_roots,
    function(root) {
      identical(normalized, root) || startsWith(normalized, paste0(root, "/"))
    },
    logical(1)
  ))
}
nonprotected_ds_store <- tree_entries[
  basename(tree_entries) == ".DS_Store" &
    !vapply(tree_entries, is_under_protected_root, logical(1))
]
if (length(nonprotected_ds_store)) {
  warning(
    "Ignoring macOS metadata files that can be recreated asynchronously:\n",
    paste(
      substring(nonprotected_ds_store, nchar(repo_root) + 2L),
      collapse = "\n"
    ),
    call. = FALSE
  )
}
forbidden_artifacts <- tree_entries[
  basename(tree_entries) %in% c("Rplots.pdf", "__pycache__") |
    grepl("\\.pyc$", basename(tree_entries), ignore.case = TRUE)
]
if (length(forbidden_artifacts)) {
  stop(
    "Generated/source-tree artifacts outside the protected core:\n",
    paste(
      substring(forbidden_artifacts, nchar(repo_root) + 2L),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

legacy_invitro_dir <- file.path(repo_root, "oxygen", "code", "in-vitro-utils")
if (dir.exists(legacy_invitro_dir)) {
  stop("Removed in-vitro-utils directory still exists: ", legacy_invitro_dir, call. = FALSE)
}

cat("O2_LAYER_BOUNDARY_CHECK=PASS\n")
cat("simulation_files=", simulation_result[["checked"]], "\n", sep = "")
cat("analysis_files=", analysis_result[["checked"]], "\n", sep = "")
cat("vis_files=", vis_result[["checked"]], "\n", sep = "")
cat("report_files=", report_result[["checked"]], "\n", sep = "")
cat("util_files=", util_result[["checked"]], "\n", sep = "")
cat(
  "compatibility_wrappers=",
  simulation_result[["compatibility_wrappers"]] +
    analysis_result[["compatibility_wrappers"]],
  "\n",
  sep = ""
)
