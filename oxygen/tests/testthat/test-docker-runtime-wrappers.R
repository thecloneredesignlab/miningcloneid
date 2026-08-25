docker_workflow_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)
docker_hpc_root <- file.path(docker_workflow_root, "Docker", "hpc")
docker_local_root <- file.path(docker_workflow_root, "Docker", "local")

testthat::test_that("Docker HPC scripts have complete one-to-one parity", {
  verifier <- file.path(docker_hpc_root, "verify_hpc_parity.sh")
  mapping_path <- file.path(docker_hpc_root, "hpc_script_mapping.tsv")
  testthat::expect_true(file.exists(verifier))
  testthat::expect_true(file.exists(mapping_path))

  mapping <- utils::read.delim(
    mapping_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  testthat::expect_equal(nrow(mapping), 28L)
  testthat::expect_false(anyDuplicated(mapping$module_path) > 0L)
  testthat::expect_false(anyDuplicated(mapping$container_path) > 0L)

  output <- system2("bash", verifier, stdout = TRUE, stderr = TRUE)
  testthat::expect_null(
    attr(output, "status"),
    info = paste(output, collapse = "\n")
  )
  testthat::expect_true(any(grepl("28 mapped files", output, fixed = TRUE)))
})

testthat::test_that("Docker runtimes use only the published image artifacts", {
  runtime_lock <- utils::read.delim(
    file.path(docker_workflow_root, "Docker", "image_runtime_lock.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  hpc_runtime <- paste(
    readLines(
      file.path(
        docker_hpc_root,
        "util",
        "o2_supply_demand_map_apptainer_runtime.sh"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  local_runtime <- paste(
    readLines(file.path(docker_local_root, "docker_runtime.sh"), warn = FALSE),
    collapse = "\n"
  )

  testthat::expect_match(
    hpc_runtime,
    "/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif",
    fixed = TRUE
  )
  testthat::expect_match(hpc_runtime, "apptainer exec", fixed = TRUE)
  testthat::expect_match(hpc_runtime, "--cleanenv", fixed = TRUE)
  testthat::expect_match(hpc_runtime, "O2SD_CONTAINER_BITMAP_TYPE", fixed = TRUE)
  testthat::expect_match(hpc_runtime, 'R_BITMAP_TYPE=${O2SD_CONTAINER_BITMAP_TYPE}', fixed = TRUE)
  testthat::expect_match(hpc_runtime, "O2SD_CONTAINER_R_PROFILE", fixed = TRUE)
  testthat::expect_match(hpc_runtime, "R_PROFILE_USER=${O2SD_CONTAINER_R_PROFILE}", fixed = TRUE)
  testthat::expect_match(hpc_runtime, "o2sd_container_r_sanity_check", fixed = TRUE)
  testthat::expect_match(hpc_runtime, "grDevices::png", fixed = TRUE)
  profile_text <- paste(
    readLines(
      file.path(
        docker_hpc_root,
        "util",
        "o2_supply_demand_map_container.Rprofile"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  testthat::expect_match(profile_text, 'options(bitmapType = "cairo")', fixed = TRUE)
  testthat::expect_match(
    local_runtime,
    "zafiro/o2_supply_demand_map:r44",
    fixed = TRUE
  )
  testthat::expect_match(local_runtime, "--platform", fixed = TRUE)
  testthat::expect_setequal(
    runtime_lock$runtime,
    c("docker_index", "docker_amd64_manifest", "apptainer_sif")
  )
  testthat::expect_true(all(grepl("^[0-9a-f]{64}$", runtime_lock$sha256)))
  testthat::expect_false(grepl(
    "\\.(pem|key|crt)|\\.ssh|docker\\.json",
    paste(hpc_runtime, local_runtime),
    ignore.case = TRUE
  ))
})

testthat::test_that("Docker fitting workers use the retrying R and PNG sanity check", {
  worker_paths <- file.path(
    docker_hpc_root,
    "array_workers",
    c(
      "submit_fit_seed_array_buffering.sub",
      "submit_fit_seed_array_invitro_buffering.sub",
      "submit_fit_seed_array_joint_buffering.sub"
    )
  )
  testthat::expect_true(all(file.exists(worker_paths)))
  worker_text <- vapply(
    worker_paths,
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1)
  )
  testthat::expect_true(all(grepl(
    "o2sd_container_r_sanity_check",
    worker_text,
    fixed = TRUE
  )))
  testthat::expect_false(any(grepl(
    'Rscript -e \'cat("Container R sanity check OK',
    worker_text,
    fixed = TRUE
  )))
})

testthat::test_that("local Docker full workflows retain production entrypoints", {
  full_fit <- paste(
    readLines(file.path(docker_local_root, "run_full_fit_docker.sh"), warn = FALSE),
    collapse = "\n"
  )
  full_analysis <- paste(
    readLines(
      file.path(docker_local_root, "run_full_analysis_docker.sh"),
      warn = FALSE
    ),
    collapse = "\n"
  )

  testthat::expect_match(full_fit, "--fitting_mode=all", fixed = TRUE)
  testthat::expect_match(full_fit, "complete fitting chain", fixed = TRUE)
  testthat::expect_false(grepl("--joint_fitting_mode=", full_fit, fixed = TRUE))
  testthat::expect_false(grepl("--joint_warmup_enable=", full_fit, fixed = TRUE))
  testthat::expect_match(full_fit, "--invivo_total_seeds=", fixed = TRUE)
  testthat::expect_match(full_fit, "--invitro_total_seeds=", fixed = TRUE)
  testthat::expect_match(full_fit, "--joint_total_seeds=", fixed = TRUE)
  testthat::expect_match(
    full_analysis,
    "runner/run_postfit_pipeline.R",
    fixed = TRUE
  )
  testthat::expect_match(
    full_analysis,
    "runner/best_fit_parameter_feature/runner.R",
    fixed = TRUE
  )
  testthat::expect_match(
    full_analysis,
    "runner/joint_coupling/run_joint_coupling_pipeline.R",
    fixed = TRUE
  )
})
