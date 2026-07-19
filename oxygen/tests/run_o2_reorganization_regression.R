#!/usr/bin/env Rscript

# End-to-end static, unit, layer, and immutable-core regression gate for the
# O2_supply_demand_MAP reorganization.
#
# Usage:
#   Rscript oxygen/tests/run_o2_reorganization_regression.R
#   Rscript oxygen/tests/run_o2_reorganization_regression.R --skip_unit=TRUE

parse_args <- function(argv) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[bits[[1L]]]] <- if (length(bits) > 1L) {
      paste(bits[-1L], collapse = "=")
    } else {
      "TRUE"
    }
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x)) return(default)
  value <- tolower(trimws(as.character(x[[1L]])))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

script_path <- function() {
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (!length(file_arg)) {
    return(normalizePath("oxygen/tests/run_o2_reorganization_regression.R"))
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

write_log <- function(lines, path) {
  writeLines(enc2utf8(as.character(lines)), path, useBytes = TRUE)
  invisible(path)
}

run_command <- function(label,
                        command,
                        args = character(),
                        log_path,
                        env = character()) {
  cat("[RUN] ", label, "\n", sep = "")
  output <- system2(
    command,
    args = shQuote(args),
    stdout = TRUE,
    stderr = TRUE,
    env = env
  )
  write_log(output, log_path)
  status <- attr(output, "status")
  if (!is.null(status)) {
    stop(
      label,
      " failed with exit status ",
      status,
      ". Log: ",
      log_path,
      "\n",
      paste(utils::tail(output, 60L), collapse = "\n"),
      call. = FALSE
    )
  }
  cat("[PASS] ", label, "\n", sep = "")
  invisible(output)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  test_script <- script_path()
  repo_root <- normalizePath(
    file.path(dirname(test_script), "..", ".."),
    mustWork = TRUE
  )
  workflow_root <- file.path(
    repo_root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP"
  )
  if (!dir.exists(workflow_root)) {
    stop("Workflow root was not found: ", workflow_root)
  }

  log_dir <- argv$log_dir
  if (is.null(log_dir) || !nzchar(trimws(log_dir))) {
    log_dir <- tempfile("o2_reorganization_regression_")
  } else if (!grepl("^/", log_dir)) {
    log_dir <- file.path(repo_root, log_dir)
  }
  log_dir <- normalizePath(log_dir, mustWork = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  cat("REPO_ROOT=", repo_root, "\n", sep = "")
  cat("WORKFLOW_ROOT=", workflow_root, "\n", sep = "")
  cat("LOG_DIR=", log_dir, "\n", sep = "")

  r_files <- list.files(
    workflow_root,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.[Rr]$"
  )
  r_files <- r_files[!grepl("/\\.rcpp_cache", r_files)]
  parse_errors <- list()
  for (path in sort(r_files)) {
    tryCatch(
      parse(file = path),
      error = function(e) {
        parse_errors[[path]] <<- conditionMessage(e)
      }
    )
  }
  if (length(parse_errors)) {
    details <- paste(
      names(parse_errors),
      unlist(parse_errors, use.names = FALSE),
      sep = ": ",
      collapse = "\n"
    )
    stop("R parse gate failed:\n", details, call. = FALSE)
  }
  cat("[PASS] R parse ", length(r_files), "/", length(r_files), "\n", sep = "")

  shell_files <- list.files(
    workflow_root,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(sh|sub)$"
  )
  for (path in sort(shell_files)) {
    run_command(
      paste0("bash syntax: ", substring(path, nchar(repo_root) + 2L)),
      "bash",
      c("-n", path),
      file.path(log_dir, paste0("bash_", basename(path), ".log"))
    )
  }
  cat(
    "[PASS] shell/Slurm syntax ",
    length(shell_files),
    "/",
    length(shell_files),
    "\n",
    sep = ""
  )

  python_files <- list.files(
    workflow_root,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.py$"
  )
  if (length(python_files)) {
    pycache <- file.path(log_dir, "python_bytecode")
    dir.create(pycache, recursive = TRUE, showWarnings = FALSE)
    run_command(
      "Python compile",
      Sys.which("python3"),
      c("-m", "py_compile", sort(python_files)),
      file.path(log_dir, "python_compile.log"),
      env = paste0("PYTHONPYCACHEPREFIX=", pycache)
    )
  }

  rscript <- file.path(R.home("bin"), "Rscript")
  immutable_script <- file.path(
    repo_root,
    "oxygen",
    "tests",
    "check_immutable_o2_core.R"
  )
  boundary_script <- file.path(
    repo_root,
    "oxygen",
    "tests",
    "check_o2_layer_boundaries.R"
  )
  unit_script <- file.path(
    repo_root,
    "oxygen",
    "tests",
    "run_unit_tests.R"
  )

  run_command(
    "layer boundary checker",
    rscript,
    boundary_script,
    file.path(log_dir, "layer_boundaries.log")
  )
  run_command(
    "immutable core before unit tests",
    rscript,
    c(immutable_script, "--full"),
    file.path(log_dir, "immutable_pre_unit.log")
  )

  if (!as_bool(argv$skip_unit, FALSE)) {
    run_command(
      "full unit regression",
      rscript,
      unit_script,
      file.path(log_dir, "unit_tests.log"),
      env = "MININGCLONEID_TEST_REPORTER=summary"
    )
  } else {
    cat("[SKIP] full unit regression\n")
  }

  run_command(
    "immutable core after unit tests",
    rscript,
    c(immutable_script, "--full"),
    file.path(log_dir, "immutable_post_unit.log")
  )
  run_command(
    "git whitespace check",
    "git",
    c(
      "-C",
      repo_root,
      "diff",
      "--check",
      "--",
      "oxygen/code/O2_supply_demand_MAP",
      "oxygen/code/in-vitro-utils",
      "oxygen/tests"
    ),
    file.path(log_dir, "git_diff_check.log")
  )

  cat("O2_REORGANIZATION_REGRESSION=PASS\n")
  cat("R_FILES=", length(r_files), "\n", sep = "")
  cat("SHELL_FILES=", length(shell_files), "\n", sep = "")
  cat("PYTHON_FILES=", length(python_files), "\n", sep = "")
  cat("LOG_DIR=", log_dir, "\n", sep = "")
}

if (sys.nframe() == 0L) {
  main()
}
