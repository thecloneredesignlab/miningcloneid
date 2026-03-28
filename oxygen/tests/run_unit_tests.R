if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Package 'testthat' is required. Please install it before running unit tests.")
}

resolve_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

resolve_test_dir <- function(script_dir) {
  candidates <- c(
    Sys.getenv("MININGCLONEID_TEST_DIR", unset = ""),
    file.path(script_dir, "testthat"),
    file.path(script_dir, "tests", "testthat"),
    file.path(getwd(), "tests", "testthat")
  )
  candidates <- candidates[nzchar(candidates)]
  hits <- candidates[dir.exists(candidates)]
  if (length(hits) == 0L) {
    stop(
      "No testthat directory found. Checked: ",
      paste(normalizePath(candidates, mustWork = FALSE), collapse = "; ")
    )
  }
  normalizePath(hits[[1]], mustWork = TRUE)
}

script_dir <- resolve_script_dir()
test_dir <- resolve_test_dir(script_dir)
reporter <- Sys.getenv("MININGCLONEID_TEST_REPORTER", unset = "location")
message("Using test dir: ", test_dir)
message("Using reporter: ", reporter)
testthat::test_dir(test_dir, reporter = reporter)
