joint_primary_root <- file.path(
  repo_info$root,
  "oxygen",
  "code",
  "O2_supply_demand_MAP"
)

joint_primary_analysis_path <- file.path(
  joint_primary_root,
  "analysis",
  "multi_warmup",
  "build_multi_warmup_landscape_tables.R"
)

joint_primary_env <- new.env(parent = globalenv())
source(joint_primary_analysis_path, local = joint_primary_env, chdir = TRUE)

testthat::test_that("primary-cluster summaries select objective-minimum seeds with deterministic ties", {
  clustered <- data.frame(
    dataset = rep("invivo", 5L),
    dataset_label = rep("in vivo", 5L),
    cluster_prefix = rep("vi", 5L),
    cluster_id = c("vi_C01", "vi_C01", "vi_C02", "vi_C02", "vi_C02"),
    cluster_base_id = c("C01", "C01", "C02", "C02", "C02"),
    cluster_num = c(1L, 1L, 2L, 2L, 2L),
    seed = c(3L, 1L, 7L, 4L, 5L),
    objective = c(4, 2, 9, 1, 1),
    tSNE1 = c(-2, -1, 1, 2, 3),
    tSNE2 = c(0, 1, 0, 1, 2),
    stringsAsFactors = FALSE
  )

  summary <- joint_primary_env$summarize_best_clusters(
    clustered,
    coord_names = c("tSNE1", "tSNE2")
  )

  testthat::expect_equal(summary$cluster_base_id, c("C01", "C02"))
  testthat::expect_equal(summary$objective_min_seed, c(1L, 4L))
  testthat::expect_equal(summary$objective_min, c(2, 1))
})

testthat::test_that("the shared embedding clusters only in-vivo best points", {
  coordinate_csv <- tempfile(fileext = ".csv")
  output_dir <- tempfile()
  coords <- data.frame(
    tSNE1 = c(-4, -3, 3, 4, -4, -3, 3, 4),
    tSNE2 = c(-4, -3, 3, 4, 4, 3, -3, -4),
    dataset = rep(c("invivo", "invitro"), each = 4L),
    point_type = "best",
    seed = rep(1:4, 2L),
    objective = c(4, 1, 3, 2, 8, 6, 5, 7),
    stringsAsFactors = FALSE
  )
  utils::write.csv(coords, coordinate_csv, row.names = FALSE)

  result <- joint_primary_env$analyze_embedding(
    reduction = "tsne",
    coordinate_csv = coordinate_csv,
    output_dir = output_dir,
    cluster_seed = 123L,
    cluster_k_min = 2L,
    cluster_k_max = 2L,
    silhouette_sample_n = 8L
  )

  testthat::expect_identical(unique(result$seed_groups$dataset), "invivo")
  testthat::expect_identical(unique(result$cluster_summary$dataset), "invivo")
  testthat::expect_true(all(startsWith(result$seed_groups$cluster_id, "vi_")))
  full_coordinates <- utils::read.csv(
    file.path(output_dir, "Tables", "pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_full_coordinates.csv"),
    stringsAsFactors = FALSE
  )
  testthat::expect_true(all(is.na(full_coordinates$cluster_id[full_coordinates$dataset == "invitro"])))
})

testthat::test_that("global in-vitro anchor uses objective minimum with a seed tie-break", {
  best_csv <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(seed = c(228L, 91L, 10L), objective = c(4, 2, 2)),
    best_csv,
    row.names = FALSE
  )
  seed_dirs <- data.frame(
    dataset = "invitro",
    seed = c(10L, 91L, 228L),
    seed_dir = paste0("/results/invitro/seed", c(10L, 91L, 228L)),
    stringsAsFactors = FALSE
  )
  anchor <- joint_primary_env$global_invitro_best_anchor(best_csv, seed_dirs)
  testthat::expect_equal(anchor$seed, 10L)
  testthat::expect_equal(anchor$objective, 2)
  testthat::expect_equal(anchor$selection_rule, "global_objective_min_seed")
})

testthat::test_that("each in-vivo primary representative pairs with one global in-vitro anchor", {
  reps <- data.frame(
    method = rep("tsne", 3L),
    dataset = rep("invivo", 3L),
    cluster_id = c("vi_C01", "vi_C02", "vi_C03"),
    cluster_base_id = c("C01", "C02", "C03"),
    cluster_num = 1:3,
    representative_rank = 1:3,
    seed = c(366L, 25L, 311L),
    seed_dir = paste0("/results/invivo/seed", c(366L, 25L, 311L)),
    stringsAsFactors = FALSE
  )
  anchor <- data.frame(
    method = "tsne",
    dataset = "invitro",
    cluster_id = "global_best",
    cluster_base_id = "global_best",
    cluster_num = 1L,
    representative_rank = 1L,
    seed = 10L,
    seed_dir = "/results/invitro/seed10",
    stringsAsFactors = FALSE
  )

  manifest <- joint_primary_env$build_manifest_from_representatives(
    invivo_reps = reps,
    invitro_anchor = anchor,
    out_dir = "/results/joint"
  )

  testthat::expect_equal(nrow(manifest), 3L)
  testthat::expect_equal(manifest$invivo_seed, c(366L, 25L, 311L))
  testthat::expect_equal(manifest$invitro_seed, rep(10L, 3L))
  testthat::expect_equal(manifest$invivo_family, c("vi_C01", "vi_C02", "vi_C03"))
  testthat::expect_equal(manifest$invitro_family, rep("global_best", 3L))
  testthat::expect_true(all(grepl("_vt_seed", manifest$warmup_label, fixed = TRUE)))
  testthat::expect_false(any(grepl("Sc", manifest$warmup_label, fixed = TRUE)))
  testthat::expect_false(any(grepl("vt_seed10_C", manifest$warmup_label, fixed = TRUE)))
  testthat::expect_equal(unique(manifest$selection_reason), "invivo_primary_cluster_objective_min_to_global_invitro_objective_min")
})

testthat::test_that("unified local and HPC joint entries expose only the fixed workflow", {
  paths <- c(
    file.path(joint_primary_root, "runner", "run_o2_fit.sh"),
    file.path(joint_primary_root, "hpc", "submit", "submit_o2_fit.sh"),
    file.path(joint_primary_root, "Docker", "hpc", "submit", "submit_o2_fit.sh")
  )
  texts <- vapply(paths, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1L))

  testthat::expect_true(all(grepl("--joint_fitting_mode has been removed", texts, fixed = TRUE)))
  testthat::expect_true(all(grepl("--invivo_run_dir", texts, fixed = TRUE)))
  testthat::expect_true(all(grepl("--invitro_run_dir", texts, fixed = TRUE)))
  testthat::expect_true(all(grepl("global-best in-vitro seed", texts, fixed = TRUE)))
  testthat::expect_false(any(grepl("case \"${JOINT_FITTING_MODE}\"", texts, fixed = TRUE)))
})

testthat::test_that("Docker array workers resolve runtime paths after Slurm spooling", {
  workers <- c(
    file.path(joint_primary_root, "Docker", "hpc", "array_workers", "run_landscape_seed_space_task.sub"),
    file.path(joint_primary_root, "Docker", "hpc", "array_workers", "run_multi_warmup_task_table_array.sub")
  )
  for (worker in workers) {
    spooled <- tempfile(pattern = "slurm_script_")
    testthat::expect_true(file.copy(worker, spooled))
    output <- suppressWarnings(system2(
      "env",
      c("-u", "O2SD_DOCKER_HPC_ROOT", paste0("PROJECT_ROOT=", repo_info$root), "bash", spooled),
      stdout = TRUE,
      stderr = TRUE
    ))
    testthat::expect_equal(attr(output, "status"), 2L, info = worker)
    testthat::expect_true(any(grepl("TASKS_TSV is required", output, fixed = TRUE)), info = worker)
    testthat::expect_false(any(grepl("No such file or directory", output, fixed = TRUE)), info = worker)
  }

  submitters <- c(
    file.path(joint_primary_root, "Docker", "hpc", "submit", "submit_o2_fit.sh"),
    file.path(joint_primary_root, "Docker", "hpc", "submit", "submit_multi_warmup_task_table.sh")
  )
  submitter_text <- paste(unlist(lapply(submitters, readLines, warn = FALSE)), collapse = "\n")
  testthat::expect_match(submitter_text, "O2SD_DOCKER_HPC_ROOT=", fixed = TRUE)
  testthat::expect_match(submitter_text, "O2SD_CONTAINER_IMAGE=", fixed = TRUE)
})

testthat::test_that("fitting_mode=all derives source directories and preserves stage dependencies", {
  capture_bash <- function(script, args) {
    output <- system2("bash", c(script, args), stdout = TRUE, stderr = TRUE)
    testthat::expect_null(attr(output, "status"), info = paste(script, args, collapse = " "))
    paste(output, collapse = "\n")
  }
  stage_positions <- function(text, labels) {
    vapply(labels, function(label) regexpr(label, text, fixed = TRUE)[[1L]], integer(1L))
  }

  local_out <- file.path(normalizePath(tempdir()), "joint-primary-all-local")
  local_text <- capture_bash(
    file.path(joint_primary_root, "runner", "run_o2_fit.sh"),
    c(
      "--fitting_mode=all",
      paste0("--project_root=", repo_info$root),
      paste0("--out_root=", local_out),
      "--invivo_run_prefix=all_invivo",
      "--invitro_run_prefix=all_invitro",
      "--joint_run_prefix=all_joint",
      "--invivo_total_seeds=1",
      "--invitro_total_seeds=1",
      "--joint_total_seeds=1",
      "--run_extra_results=FALSE",
      "--dry_run=TRUE"
    )
  )
  local_positions <- stage_positions(
    local_text,
    c("Run in vivo fit:", "Run in vitro seed 1:", "Run joint primary-cluster sweep:")
  )
  testthat::expect_true(all(local_positions > 0L))
  testthat::expect_true(all(diff(local_positions) > 0L))
  testthat::expect_match(local_text, paste0("--invivo_run_dir=", local_out, "/all_invivo"), fixed = TRUE)
  testthat::expect_match(local_text, paste0("--invitro_run_dir=", local_out, "/all_invitro"), fixed = TRUE)

  hpc_scripts <- c(
    file.path(joint_primary_root, "hpc", "submit", "submit_o2_fit.sh"),
    file.path(joint_primary_root, "Docker", "hpc", "submit", "submit_o2_fit.sh")
  )
  for (idx in seq_along(hpc_scripts)) {
    out_root <- file.path(normalizePath(tempdir()), paste0("joint-primary-all-hpc-", idx))
    text <- capture_bash(
      hpc_scripts[[idx]],
      c(
        "--fitting_mode=all",
        paste0("--project_root=", repo_info$root),
        paste0("--out_root=", out_root),
        "--invivo_run_prefix=all_invivo",
        "--invitro_run_prefix=all_invitro",
        "--joint_run_prefix=all_joint",
        "--invivo_total_seeds=1",
        "--invivo_array_tasks=1",
        "--invivo_seeds_per_task=1",
        "--invitro_total_seeds=1",
        "--invitro_array_tasks=1",
        "--invitro_seeds_per_task=1",
        "--joint_total_seeds=1",
        "--dry_run=TRUE"
      )
    )
    positions <- stage_positions(
      text,
      c("Submit in vivo:", "Submit in vitro:", "Submit multi-warm-up controller:")
    )
    testthat::expect_true(all(positions > 0L), info = hpc_scripts[[idx]])
    testthat::expect_true(all(diff(positions) > 0L), info = hpc_scripts[[idx]])
    testthat::expect_match(text, "--dependency=afterok:DRYRUN_INVIVO_JOB", fixed = TRUE)
    testthat::expect_match(text, "--dependency=afterok:DRYRUN_INVITRO_JOB", fixed = TRUE)
    testthat::expect_match(text, paste0("--invivo_run_dir=", out_root, "/all_invivo"), fixed = TRUE)
    testthat::expect_match(text, paste0("--invitro_run_dir=", out_root, "/all_invitro"), fixed = TRUE)
  }
})
