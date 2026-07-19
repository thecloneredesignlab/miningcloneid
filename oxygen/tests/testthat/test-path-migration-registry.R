migration_workflow_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)
migration_registry_path <- file.path(
  migration_workflow_root,
  "docs",
  "path_migration_table.tsv"
)

testthat::test_that("path migration registry is complete and points to live canonical files", {
  testthat::expect_true(file.exists(migration_registry_path))
  registry <- utils::read.delim(
    migration_registry_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required <- c(
    "stage",
    "old_path",
    "canonical_path",
    "responsibility",
    "status"
  )
  testthat::expect_identical(names(registry), required)
  testthat::expect_true(nrow(registry) > 0L)
  testthat::expect_true(all(registry$stage %in% 1:5))
  testthat::expect_true(all(nzchar(registry$old_path)))
  testthat::expect_true(all(nzchar(registry$canonical_path)))
  testthat::expect_true(all(nzchar(registry$responsibility)))
  testthat::expect_true(all(nzchar(registry$status)))
  testthat::expect_false(anyDuplicated(
    paste(registry$stage, registry$old_path, registry$canonical_path, sep = "\t")
  ) > 0L)
  expected_order <- order(
    registry$stage,
    registry$old_path,
    registry$canonical_path,
    registry$status
  )
  testthat::expect_identical(expected_order, seq_len(nrow(registry)))

  is_literal <- !grepl("[*?\\[]", registry$canonical_path)
  canonical_absolute <- file.path(
    repo_info$root,
    registry$canonical_path[is_literal]
  )
  missing_canonical <- registry$canonical_path[is_literal][
    !file.exists(canonical_absolute)
  ]
  testthat::expect_length(missing_canonical, 0L)

  retained_status <- c(
    "compatibility-wrapper",
    "deprecated-compatibility-wrapper",
    "split-in-place"
  )
  retained <- registry$status %in% retained_status &
    !grepl("[*?\\[]", registry$old_path)
  old_absolute <- file.path(repo_info$root, registry$old_path[retained])
  missing_retained <- registry$old_path[retained][!file.exists(old_absolute)]
  testthat::expect_length(missing_retained, 0L)

  protected_targets <- grepl(
    "/(model|optimizer)/",
    paste0("/", registry$canonical_path)
  )
  testthat::expect_false(any(protected_targets))
})
