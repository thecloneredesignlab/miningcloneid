postprocess_relocation_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

testthat::test_that("extra-results workers survive Slurm script relocation", {
  scripts <- c(
    file.path(
      postprocess_relocation_root,
      "hpc",
      "postprocess",
      "postprocess_extra_results.sh"
    ),
    file.path(
      postprocess_relocation_root,
      "Docker",
      "hpc",
      "postprocess",
      "postprocess_extra_results.sh"
    )
  )
  testthat::expect_true(all(file.exists(scripts)))

  run_dir <- tempfile("o2sd-extra-results-run-")
  dir.create(run_dir)
  on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)

  for (script in scripts) {
    relocated <- tempfile("slurm-spool-script-", fileext = ".sh")
    testthat::expect_true(file.copy(script, relocated, overwrite = TRUE))
    output <- system2(
      "bash",
      c(
        relocated,
        paste0("--project_root=", repo_info$root),
        paste0("--run_dir=", run_dir),
        "--dry_run=TRUE"
      ),
      stdout = TRUE,
      stderr = TRUE,
      env = c(
        "O2_POSTPROCESS_LOGIN_SHELL=1",
        paste0("PROJECT_ROOT=", repo_info$root)
      )
    )
    testthat::expect_null(
      attr(output, "status"),
      info = paste(script, paste(output, collapse = "\n"))
    )
    testthat::expect_true(any(grepl(
      "DRY_RUN=TRUE; not running extra_results.R",
      output,
      fixed = TRUE
    )))
  }
})
