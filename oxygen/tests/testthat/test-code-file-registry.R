registry_workflow_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)
registry_script <- file.path(
  registry_workflow_root,
  "runner",
  "documentation",
  "generate_code_file_registry.R"
)

testthat::test_that("per-file registry generator parses and documents all source layers", {
  testthat::expect_true(file.exists(registry_script))
  testthat::expect_silent(parse(registry_script))

  out_dir <- tempfile("o2_file_registry_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  out_md <- file.path(out_dir, "CODE_FILE_REGISTRY.md")
  out_tsv <- file.path(out_dir, "code_file_registry.tsv")
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    shQuote(c(
      registry_script,
      paste0("--workflow_root=", registry_workflow_root),
      paste0("--out_md=", out_md),
      paste0("--out_tsv=", out_tsv)
    )),
    stdout = TRUE,
    stderr = TRUE
  )
  testthat::expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  testthat::expect_true(file.exists(out_md))
  testthat::expect_true(file.exists(out_tsv))

  registry <- read.delim(out_tsv, check.names = FALSE)
  testthat::expect_true(all(c(
    "file",
    "layer",
    "kind",
    "purpose",
    "inputs",
    "outputs",
    "functions",
    "tests",
    "lines"
  ) %in% names(registry)))
  testthat::expect_true(all(nzchar(registry$purpose)))
  testthat::expect_true(all(c(
    "model",
    "optimizer",
    "util",
    "simulation",
    "analysis",
    "vis",
    "report",
    "runner",
    "hpc"
  ) %in% registry$layer))
  testthat::expect_true(all(c(
    "model/model_O2_supply_demand_MAP.R",
    "optimizer/fit_model_O2_supply_demand_MAP.R",
    "runner/documentation/generate_code_file_registry.R"
  ) %in% registry$file))
})

testthat::test_that("published per-file registry matches the current source tree", {
  published_md <- file.path(
    registry_workflow_root,
    "docs",
    "CODE_FILE_REGISTRY.md"
  )
  published_tsv <- file.path(
    registry_workflow_root,
    "docs",
    "code_file_registry.tsv"
  )
  testthat::expect_true(file.exists(published_md))
  testthat::expect_true(file.exists(published_tsv))

  source_files <- list.files(
    registry_workflow_root,
    recursive = TRUE,
    full.names = TRUE
  )
  source_files <- source_files[
    file.info(source_files)$isdir %in% FALSE &
      grepl("(\\.R|\\.Rmd|\\.py|\\.sh|\\.sub|\\.cpp)$", source_files) &
      !grepl("/\\.rcpp_cache", source_files)
  ]
  source_files <- sort(sub(
    paste0("^", registry_workflow_root, "/"),
    "",
    normalizePath(source_files, mustWork = TRUE)
  ))

  registry <- utils::read.delim(
    published_tsv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  testthat::expect_false(anyDuplicated(registry$file) > 0L)
  testthat::expect_setequal(registry$file, source_files)

  markdown <- paste(readLines(published_md, warn = FALSE), collapse = "\n")
  missing_markdown <- source_files[
    !vapply(
      source_files,
      function(path) grepl(path, markdown, fixed = TRUE),
      logical(1)
    )
  ]
  testthat::expect_length(missing_markdown, 0L)
})
