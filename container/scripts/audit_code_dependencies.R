#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: audit_code_dependencies.R PROJECT_ROOT OUT_DIR", call. = FALSE)
}

project_root <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_rel <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(project_root, .Platform$file.sep)
  sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", prefix)), "", path)
}

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA"
  )
}

classify_scope <- function(rel) {
  if (grepl("^oxygen/code/O2_supply_demand_MAP/analysis/", rel)) return("analysis")
  if (grepl("^oxygen/code/O2_supply_demand_MAP/simulation/", rel)) return("simulation")
  if (grepl("^oxygen/code/O2_supply_demand_MAP/vis/", rel)) return("visualization")
  if (grepl("^oxygen/code/O2_supply_demand_MAP/report/", rel)) return("report")
  if (grepl("^oxygen/code/O2_supply_demand_MAP/runner/", rel)) return("runner")
  if (grepl("^oxygen/code/O2_supply_demand_MAP/model/", rel)) return("model")
  if (grepl("^oxygen/code/O2_supply_demand_MAP/util/", rel)) return("util")
  if (grepl("^oxygen/code/O2_supply_demand_MAP/hpc/", rel)) return("hpc")
  if (grepl("^oxygen/tests/", rel)) return("test")
  if (grepl("^oxygen/", rel)) return("oxygen_auxiliary")
  if (grepl("\\.py$", rel, ignore.case = TRUE)) return("repository_python_other")
  "repository_other"
}

scan_roots <- c(
  file.path(project_root, "oxygen", "code", "O2_supply_demand_MAP"),
  file.path(project_root, "oxygen", "scripts"),
  file.path(project_root, "oxygen", "tests"),
  file.path(project_root, "code")
)
scan_roots <- scan_roots[dir.exists(scan_roots)]
all_paths <- unlist(lapply(
  scan_roots,
  list.files,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
), use.names = FALSE)
all_paths <- c(
  all_paths,
  list.files(project_root, recursive = FALSE, full.names = TRUE, all.files = TRUE, no.. = TRUE),
  list.files(
    file.path(project_root, "oxygen", "figures"),
    pattern = "\\.py$",
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    file.path(project_root, "docs"),
    pattern = "\\.py$",
    recursive = TRUE,
    full.names = TRUE
  )
)
all_paths <- sort(unique(all_paths[file.exists(all_paths) & !dir.exists(all_paths)]))
all_rel <- vapply(all_paths, normalize_rel, character(1))
keep <- !grepl("(^|/)(\\.git|container|oxygen/results|\\.rcpp_cache[^/]*)(/|$)", all_rel) &
  grepl("\\.(R|Rmd|qmd|py|sh|sub|sbatch|slurm|cpp|h|hpp|yaml|yml)$", all_rel, ignore.case = TRUE)
code_paths <- all_paths[keep]
code_rel <- all_rel[keep]

code_scope <- data.frame(
  file = code_rel,
  scope = vapply(code_rel, classify_scope, character(1)),
  extension = tolower(tools::file_ext(code_rel)),
  size_bytes = as.numeric(file.info(code_paths)$size),
  stringsAsFactors = FALSE
)
code_scope <- code_scope[order(code_scope$scope, code_scope$file), , drop = FALSE]
write_tsv(code_scope, file.path(out_dir, "code-scope.tsv"))

empty_usage <- data.frame(
  file = character(),
  line = integer(),
  package = character(),
  call_style = character(),
  detection_method = character(),
  code_scope = character(),
  dependency_scope = character(),
  stringsAsFactors = FALSE
)
usage_rows <- list()
source_rows <- list()
external_rows <- list()

add_usage <- function(file, line, package, style, method) {
  if (!nzchar(package) || package %in% c("NA", "R")) return(invisible(NULL))
  rel <- normalize_rel(file)
  scope <- classify_scope(rel)
  dependency_scope <- if (scope == "analysis") {
    "analysis_direct"
  } else if (scope %in% c("simulation", "visualization", "report", "runner", "model", "util", "hpc")) {
    "o2_runtime"
  } else if (scope == "test") {
    "o2_test"
  } else if (scope == "oxygen_auxiliary") {
    "o2_auxiliary"
  } else {
    "repository_other"
  }
  usage_rows[[length(usage_rows) + 1L]] <<- data.frame(
    file = rel,
    line = as.integer(line),
    package = package,
    call_style = style,
    detection_method = method,
    code_scope = scope,
    dependency_scope = dependency_scope,
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

add_source <- function(file, line, target, style) {
  source_rows[[length(source_rows) + 1L]] <<- data.frame(
    file = normalize_rel(file),
    line = as.integer(line),
    target_expression = target,
    call_style = style,
    stringsAsFactors = FALSE
  )
}

add_external <- function(file, line, command, style) {
  external_rows[[length(external_rows) + 1L]] <<- data.frame(
    file = normalize_rel(file),
    line = as.integer(line),
    command = command,
    call_style = style,
    stringsAsFactors = FALSE
  )
}

call_name <- function(x) {
  if (!is.call(x)) return("")
  head <- x[[1L]]
  if (is.symbol(head)) return(as.character(head))
  ""
}

scalar_text <- function(x) {
  if (is.character(x) && length(x) == 1L) return(x)
  if (is.symbol(x)) return(as.character(x))
  paste(deparse(x, width.cutoff = 500L), collapse = "")
}

walk_expr <- function(x, file, inherited_line = NA_integer_) {
  sr <- attr(x, "srcref")
  line <- inherited_line
  if (!is.null(sr) && length(sr) >= 1L) line <- as.integer(sr[[1L]])
  if (is.call(x)) {
    op <- call_name(x)
    if (op %in% c("::", ":::") && length(x) >= 3L) {
      add_usage(file, line, scalar_text(x[[2L]]), op, "R_AST")
    }
    if (op %in% c("library", "require", "requireNamespace", "loadNamespace", "getNamespace") &&
        length(x) >= 2L) {
      package_arg <- x[[2L]]
      if (is.character(package_arg) ||
          (is.symbol(package_arg) && op %in% c("library", "require"))) {
        add_usage(file, line, scalar_text(package_arg), op, "R_AST")
      }
    }
    if (op %in% c("source", "sys.source") && length(x) >= 2L) {
      add_source(file, line, scalar_text(x[[2L]]), op)
    }
    if (op %in% c("system", "system2") && length(x) >= 2L) {
      add_external(file, line, scalar_text(x[[2L]]), op)
    }
    for (i in seq_along(x)) walk_expr(x[[i]], file, line)
  } else if (is.pairlist(x) || is.expression(x)) {
    for (i in seq_along(x)) walk_expr(x[[i]], file, line)
  }
  invisible(NULL)
}

scan_lines <- function(file, lines) {
  namespace_pattern <- "([A-Za-z][A-Za-z0-9.]*)[[:space:]]*:::{0,1}"
  loader_pattern <- "(library|require|requireNamespace|loadNamespace|getNamespace)[[:space:]]*\\([[:space:]]*['\"]?([A-Za-z][A-Za-z0-9.]*)"
  for (i in seq_along(lines)) {
    line <- lines[[i]]
    ns_hits <- gregexpr(namespace_pattern, line, perl = TRUE)
    hits <- regmatches(line, ns_hits)[[1L]]
    if (length(hits) && !identical(hits, character(0)) && !identical(hits, "-1")) {
      for (hit in hits) {
        package <- sub("[[:space:]]*:::{0,1}$", "", hit, perl = TRUE)
        add_usage(file, i, package, "namespace_operator", "line_scan")
      }
    }
    loader_hit <- regexec(loader_pattern, line, perl = TRUE)
    groups <- regmatches(line, loader_hit)[[1L]]
    if (length(groups) >= 3L) add_usage(file, i, groups[[3L]], groups[[2L]], "line_scan")
  }
}

r_paths <- code_paths[grepl("\\.(R|Rmd|qmd)$", code_paths, ignore.case = TRUE)]
for (file in r_paths) {
  lines <- readLines(file, warn = FALSE)
  scan_lines(file, lines)
  ext <- tolower(tools::file_ext(file))
  if (ext == "r") {
    parsed <- tryCatch(parse(file, keep.source = TRUE), error = function(e) NULL)
    if (!is.null(parsed)) walk_expr(parsed, file)
  } else {
    in_chunk <- FALSE
    chunk <- character()
    chunk_start <- 0L
    for (i in seq_along(lines)) {
      if (!in_chunk && grepl("^```\\{[rR]([,}]|$)", lines[[i]])) {
        in_chunk <- TRUE
        chunk <- character()
        chunk_start <- i
      } else if (in_chunk && grepl("^```[[:space:]]*$", lines[[i]])) {
        parsed <- tryCatch(parse(text = chunk, keep.source = TRUE), error = function(e) NULL)
        if (!is.null(parsed)) walk_expr(parsed, file, chunk_start + 1L)
        in_chunk <- FALSE
      } else if (in_chunk) {
        chunk <- c(chunk, lines[[i]])
      }
    }
  }
}

usage <- if (length(usage_rows)) do.call(rbind, usage_rows) else empty_usage
if (nrow(usage)) {
  usage <- unique(usage)
  ast_keys <- paste(
    usage$file[usage$detection_method == "R_AST"],
    usage$package[usage$detection_method == "R_AST"],
    sep = "\r"
  )
  usage_keys <- paste(usage$file, usage$package, sep = "\r")
  usage <- usage[usage$detection_method == "R_AST" | usage_keys %in% ast_keys, , drop = FALSE]
  usage <- usage[order(usage$file, usage$line, usage$package, usage$detection_method), , drop = FALSE]
}
write_tsv(usage, file.path(out_dir, "r-package-usage-by-file.tsv"))

direct <- if (nrow(usage)) {
  aggregate(
    file ~ package + dependency_scope,
    data = usage,
    FUN = function(x) paste(sort(unique(x)), collapse = ";")
  )
} else {
  data.frame(package = character(), dependency_scope = character(), file = character())
}
names(direct)[names(direct) == "file"] <- "referencing_files"
direct <- direct[order(direct$dependency_scope, direct$package), , drop = FALSE]
write_tsv(direct, file.path(out_dir, "r-direct-packages.tsv"))

sources <- if (length(source_rows)) unique(do.call(rbind, source_rows)) else data.frame(
  file = character(), line = integer(), target_expression = character(), call_style = character()
)
sources <- sources[order(sources$file, sources$line), , drop = FALSE]
write_tsv(sources, file.path(out_dir, "r-source-calls.tsv"))

external <- if (length(external_rows)) unique(do.call(rbind, external_rows)) else data.frame(
  file = character(), line = integer(), command = character(), call_style = character()
)
external <- external[order(external$file, external$line), , drop = FALSE]
write_tsv(external, file.path(out_dir, "r-external-program-usage.tsv"))

summary <- data.frame(
  metric = c(
    "analysis_R_files",
    "analysis_Python_files",
    "all_scanned_R_like_files",
    "R_package_usage_rows",
    "analysis_direct_packages",
    "all_detected_packages",
    "source_calls",
    "external_program_calls"
  ),
  value = c(
    sum(code_scope$scope == "analysis" & code_scope$extension == "r"),
    sum(code_scope$scope == "analysis" & code_scope$extension == "py"),
    length(r_paths),
    nrow(usage),
    length(unique(usage$package[usage$dependency_scope == "analysis_direct"])),
    length(unique(usage$package)),
    nrow(sources),
    nrow(external)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, file.path(out_dir, "code-audit-summary.tsv"))
