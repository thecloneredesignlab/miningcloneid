# Targeted, independent Monte Carlo precision rescue for Figure 6.
#
# The failed full-grid run is used only as a read-only pilot and cache source.
# The one failed grid condition is recomputed from entirely new RNG streams at
# a fixed repeat count.  No production draw from the pilot is reused.

f6pr_condition <- function() list(
  model_context = "in vitro", pair_label = "C02", initial_ploidy = 3,
  O2_pct = 0.7, p_misseg = 0.01
)

f6pr_assert_base_root <- function(paths, base_run_root) {
  run_base <- normalizePath(
    file.path(paths$figure6, "finite_time_full_q10_runs"), mustWork = TRUE
  )
  base_run_root <- normalizePath(base_run_root, mustWork = TRUE)
  prefix <- paste0(run_base, .Platform$file.sep)
  if (!startsWith(paste0(base_run_root, .Platform$file.sep), prefix)) {
    stop("Base run must be inside iteration4 Figure6 full-range runs.")
  }
  base_run_root
}

f6pr_compare_manifest <- function(observed, expected, label) {
  common <- intersect(names(observed), names(expected))
  common <- setdiff(common, "cache_path")
  if (!length(common) || nrow(observed) != nrow(expected) ||
      !isTRUE(all.equal(
        observed[common], expected[common], check.attributes = FALSE
      ))) {
    stop("Rescue ", label, " does not match the immutable base run.")
  }
  invisible(TRUE)
}

f6pr_rescue_fingerprint <- function(
    base_fingerprint, master_seed, replicates, paths
) {
  source_paths <- c(
    file.path(paths$code, "util", "analysis", "figure6_stochastic_passage.R"),
    file.path(paths$code, "util", "analysis", "figure6_full_range_q10.R"),
    file.path(paths$code, "util", "analysis", "figure6_precision_rescue.R")
  )
  paste(
    f6g_profile(), "precision_rescue_v1", paste0("base=", base_fingerprint),
    paste0("master_seed=", master_seed), paste0("replicates=", replicates),
    paste0("sources=", paste(unname(tools::md5sum(source_paths)), collapse = ":")),
    sep = "|"
  )
}

f6pr_validate_stochastic_kernel <- function(
    paths, run_paths, endpoint_manifest, objective_bundle, contexts,
    passage_bundle
) {
  source(file.path(
    paths$code, "util", "analysis", "figure6_stochastic_validation.R"
  ))
  cases <- list()
  for (family in f6ft_family_levels()) {
    endpoint <- endpoint_manifest$endpoints[
      endpoint_manifest$endpoints$model_context == "in vitro" &
        endpoint_manifest$endpoints$pair_label == family, , drop = FALSE
    ][1L, ]
    prepared <- f6ft_prepare_endpoint(endpoint, objective_bundle, contexts)
    fixed <- fixo2_fixed_matrix(
      globalenv(), prepared$config,
      figure6_force_p_misseg(prepared$run_params, 0.01), O2 = 0.7
    )
    unit <- prepared$config$N_UNIT %||% 22L
    cases[[length(cases) + 1L]] <- list(
      name = paste0(family, "_precision_rescue"),
      step = as.matrix(Matrix::expm(fixed$M)),
      initial = f6ft_initial_matrix(fixed$ngrid, unit, f6ft_initial_ploidy()),
      ploidy = as.numeric(fixed$ngrid) / unit, days = 25L,
      seed = passage_bundle$rule$seed_cells[[1L]],
      target = passage_bundle$rule$target_live_cells[[1L]],
      duration = as.integer(
        passage_bundle$rule$representative_passage_duration_day[[1L]]
      )
    )
  }
  f6s_validate(paths, run_paths, cases)
}

f6pr_compute_endpoint <- function(
    endpoint, objective_bundle, contexts, passage_bundle, fingerprint,
    checkpoint_root, replicates, condition
) {
  prepared <- f6ft_prepare_endpoint(endpoint, objective_bundle, contexts)
  forced <- figure6_force_p_misseg(prepared$run_params, condition$p_misseg)
  formula <- figure6_p_misseg_formula_qc(
    prepared$run_params, condition$p_misseg
  )
  fixed <- fixo2_fixed_matrix(
    globalenv(), prepared$config, forced, O2 = condition$O2_pct
  )
  unit <- prepared$config$N_UNIT %||% 22L
  initial <- f6ft_initial_matrix(
    fixed$ngrid, unit, condition$initial_ploidy
  )
  checkpoint <- file.path(
    checkpoint_root,
    paste0("endpoint_", endpoint$endpoint_group[[1L]], ".rds")
  )
  result <- f6s_operator(
    as.matrix(Matrix::expm(fixed$M)), initial,
    as.numeric(fixed$ngrid) / unit, endpoint, condition$O2_pct,
    condition$p_misseg, max(f6g_days(FALSE)), passage_bundle,
    checkpoint_path = checkpoint, fingerprint = paste(
      fingerprint, endpoint$endpoint_group[[1L]], sep = "|"
    ), keep_trace = FALSE, replicates = replicates,
    initial_ploidy_values = condition$initial_ploidy
  )
  list(
    endpoint_group = endpoint$endpoint_group[[1L]],
    endpoint_multiplicity_q10 = endpoint$endpoint_multiplicity_q10[[1L]],
    formula_error = formula$maximum_direct_formula_error,
    result = result, checkpoint = checkpoint
  )
}

f6pr_combine_condition <- function(results, condition, replicates, target) {
  represented <- sum(vapply(
    results, `[[`, numeric(1L), "endpoint_multiplicity_q10"
  ))
  if (represented != 50L) stop("Rescue endpoints do not restore q10 weight 50.")
  reduce_field <- function(field) Reduce(
    `+`, lapply(results, function(x) x$result[[field]])
  )
  sum_mean <- reduce_field("sum")
  sum_square <- reduce_field("sum_squared_endpoint_mean")
  sum_within <- reduce_field("sum_within_variance")
  sum_mc_variance <- reduce_field("sum_mc_variance")
  mean_ploidy <- as.numeric(sum_mean[1L, ] / represented)
  within_variance <- as.numeric(sum_within[1L, ] / represented)
  mcse <- as.numeric(sqrt(sum_mc_variance[1L, ]) / represented)
  between_variance <- as.numeric(pmax(
    (sum_square[1L, ] - represented * mean_ploidy^2) /
      (represented - 1L), 0
  ))
  summary <- do.call(rbind, lapply(results, function(x) x$result$summary))
  validation <- data.frame(
    model_context = condition$model_context,
    pair_label = condition$pair_label,
    initial_ploidy = condition$initial_ploidy,
    O2_pct = condition$O2_pct,
    p_misseg = condition$p_misseg,
    optimizer_endpoint_weight = represented,
    n_unique_parameter_identity = length(results),
    n_replicate_per_optimizer_endpoint = replicates,
    production_trajectory_count = represented * replicates,
    maximum_mcse_N = max(mcse),
    maximum_mcse_day = which.max(mcse) - 1L,
    mcse_final_N = tail(mcse, 1L),
    within_endpoint_sd_final_N = sqrt(tail(within_variance, 1L)),
    mcse_target_N = target,
    maximum_direct_formula_error = max(vapply(
      results, `[[`, numeric(1L), "formula_error"
    )),
    independent_from_base_production = TRUE,
    passed = max(mcse) <= target,
    stringsAsFactors = FALSE
  )
  list(
    mean = mean_ploidy, within = within_variance, mcse = mcse,
    between = between_variance, summary = summary,
    validation = validation
  )
}

f6pr_patch_task <- function(
    base_cache, output, combined, condition, rescue_fingerprint, target
) {
  object <- readRDS(base_cache)
  initial_index <- match(
    condition$initial_ploidy, f6ft_initial_ploidy()
  )
  oxygen_index <- which(abs(object$o2_values - condition$O2_pct) < 1e-12)
  if (length(initial_index) != 1L || length(oxygen_index) != 1L ||
      abs(object$task$p_misseg[[1L]] - condition$p_misseg) > 1e-12 ||
      object$task$pair_label[[1L]] != condition$pair_label) {
    stop("Rescue condition does not map uniquely into the selected task cache.")
  }
  object$passage_weighted_mean[initial_index, , oxygen_index] <- combined$mean
  object$passage_feasible_weight[initial_index, , oxygen_index] <- 50L
  object$stochastic_within_endpoint_variance[
    initial_index, , oxygen_index
  ] <- combined$within
  object$stochastic_mcse[initial_index, , oxygen_index] <- combined$mcse
  object$between_endpoint_mean_variance[
    initial_index, , oxygen_index
  ] <- combined$between

  replace <- object$passage$pair_label == condition$pair_label &
    abs(object$passage$p_misseg - condition$p_misseg) < 1e-12 &
    abs(object$passage$O2_pct - condition$O2_pct) < 1e-12 &
    object$passage$initial_ploidy == condition$initial_ploidy
  if (sum(replace) != 50L || nrow(combined$summary) != 50L) {
    stop("Expected exactly 50 optimizer-endpoint summary rows to replace.")
  }
  old_passage_count <- sum(object$passage$passage_count[replace])
  object$passage <- rbind(object$passage[!replace, , drop = FALSE], combined$summary)
  new_passage_count <- sum(combined$summary$passage_count)
  object$qc$total_passage_events <- object$qc$total_passage_events -
    old_passage_count + new_passage_count
  object$qc$maximum_mcse_N <- max(object$stochastic_mcse)
  object$qc$mcse_target_met <- object$qc$maximum_mcse_N <= target
  object$qc$n_stochastic_replicate <- max(object$passage$n_replicate)
  object$qc$minimum_stochastic_replicate <- min(object$passage$n_replicate)
  object$qc$passed <- isTRUE(object$qc$passed) && object$qc$mcse_target_met
  object$qc$cache_path <- output
  object$fingerprint <- rescue_fingerprint
  object$precision_rescue <- list(
    condition = condition, repeat_count = unique(combined$summary$n_replicate),
    independent_random_streams = TRUE,
    rescued_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    source_cache = normalizePath(base_cache, mustWork = TRUE),
    source_cache_md5 = unname(tools::md5sum(base_cache))
  )
  f6ft_atomic_save_rds(object, output, compress = FALSE)
  invisible(object$qc)
}

f6pr_data <- function(
    workspace_root = f6r_find_workspace_root(), base_run_root, run_id,
    n_core = 16L, replicates = 200L, master_seed = 20260907L,
    publish_current = TRUE
) {
  n_core <- as.integer(n_core)
  replicates <- as.integer(replicates)
  master_seed <- as.integer(master_seed)
  stopifnot(n_core >= 1L, replicates == 200L, master_seed > 0L)
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1",
    FIGURE6_STOCHASTIC_REPLICATES = as.character(replicates),
    FIGURE6_STOCHASTIC_ALLOCATION = "fixed",
    FIGURE6_STOCHASTIC_MASTER_SEED = as.character(master_seed)
  )
  f6r_require_packages(c(
    "Matrix", "Rcpp", "data.table", "future", "future.apply"
  ))
  paths <- f6r_paths(workspace_root)
  base_run_root <- f6pr_assert_base_root(paths, base_run_root)
  run_paths <- f6g_paths(paths, run_id = run_id, create = TRUE)
  if (identical(normalizePath(run_paths$run_root, mustWork = TRUE), base_run_root)) {
    stop("Rescue run must not overwrite the failed base run.")
  }
  condition <- f6pr_condition()
  f6r_load_response_engine(paths)
  f6g_load_propagator(paths)
  objective_bundle <- f6r_objective_selection(paths)
  endpoint_manifest <- f6ft_build_endpoint_manifest(
    paths, objective_bundle, run_paths
  )
  f6pr_compare_manifest(
    endpoint_manifest$endpoints,
    f6r_read_tsv(file.path(base_run_root, "q10_unique_endpoint_manifest.tsv")),
    "unique endpoint manifest"
  )
  f6pr_compare_manifest(
    endpoint_manifest$expanded,
    f6r_read_tsv(file.path(base_run_root, "q10_optimizer_seed_manifest.tsv")),
    "expanded endpoint manifest"
  )
  f6g_chart_contract(run_paths)
  schedule <- f6p_extract_schedule(run_paths, paths)
  passage_bundle <- f6g_canonical_passage_rule(schedule, run_paths)
  passage_bundle$stochastic <- f6s_prepare(endpoint_manifest, run_paths)

  base_tasks <- f6r_read_tsv(file.path(base_run_root, "full_range_task_manifest.tsv"))
  target_o2_index <- which.min(abs(
    f6g_o2(condition$model_context) - condition$O2_pct
  ))
  if (abs(f6g_o2(condition$model_context)[target_o2_index] -
      condition$O2_pct) > 1e-12) {
    stop("Rescue oxygen value is absent from the full-range grid.")
  }
  selected <- base_tasks[
    base_tasks$model_context == condition$model_context &
      base_tasks$pair_label == condition$pair_label &
      abs(base_tasks$p_misseg - condition$p_misseg) < 1e-12 &
      base_tasks$o2_index_start <= target_o2_index &
      base_tasks$o2_index_end >= target_o2_index, , drop = FALSE
  ]
  if (nrow(selected) != 1L) stop("Failed rescue task was not uniquely identified.")
  base_cache <- file.path(
    base_run_root, "task_cache", basename(selected$cache_path[[1L]])
  )
  f6r_require_files(base_cache, "failed task cache")
  base_object <- readRDS(base_cache)
  base_fingerprint <- base_object$fingerprint
  if (isTRUE(base_object$qc$mcse_target_met) ||
      abs(base_object$qc$maximum_mcse_N - 0.0185020743134283) > 1e-12) {
    stop("Base cache is not the audited Figure 6 precision failure.")
  }
  rescue_fingerprint <- f6pr_rescue_fingerprint(
    base_fingerprint, master_seed, replicates, paths
  )
  contexts <- lapply(
    unique(endpoint_manifest$endpoints$pair_id),
    f6r_pair_model_context, selected = objective_bundle$selected, paths = paths
  )
  names(contexts) <- unique(endpoint_manifest$endpoints$pair_id)
  f6pr_validate_stochastic_kernel(
    paths, run_paths, endpoint_manifest, objective_bundle, contexts,
    passage_bundle
  )
  endpoint_indices <- as.integer(strsplit(
    selected$endpoint_indices[[1L]], ",", fixed = TRUE
  )[[1L]])
  endpoints <- endpoint_manifest$endpoints[
    match(endpoint_indices, endpoint_manifest$endpoints$endpoint_index),
    , drop = FALSE
  ]
  if (anyNA(endpoints$endpoint_index) || nrow(endpoints) != 16L ||
      sum(endpoints$endpoint_multiplicity_q10) != 50L) {
    stop("Unexpected optimizer-endpoint composition for the rescue condition.")
  }
  checkpoint_root <- file.path(run_paths$run_root, "precision_rescue_checkpoints")
  dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
  endpoint_list <- split(endpoints, seq_len(nrow(endpoints)))
  compute_one <- function(endpoint) tryCatch(
    f6pr_compute_endpoint(
      endpoint, objective_bundle, contexts, passage_bundle,
      rescue_fingerprint, checkpoint_root, replicates, condition
    ), error = function(e) structure(
      list(endpoint_group = endpoint$endpoint_group[[1L]],
           message = conditionMessage(e)), class = "f6pr_error"
    )
  )
  message(
    "Figure 6 precision rescue: condition=C02/3N/O2=0.7/p=0.01, ",
    "unique_parameter_identities=", nrow(endpoints), ", endpoints=50, R=",
    replicates, ", workers=", min(n_core, nrow(endpoints))
  )
  results <- f6ft_parallel_lapply(endpoint_list, compute_one, n_core = n_core)
  failed <- vapply(results, inherits, logical(1L), "f6pr_error")
  if (any(failed)) stop(
    "Figure 6 precision-rescue endpoint failures: ",
    paste(vapply(results[failed], `[[`, character(1L), "message"), collapse = "; ")
  )
  combined <- f6pr_combine_condition(
    results, condition, replicates,
    passage_bundle$stochastic$config$mcse_target_N
  )
  f6ft_atomic_save_rds(
    combined,
    file.path(run_paths$run_root, "precision_rescue_condition_result.rds"),
    compress = "gzip"
  )
  validation_path <- f6ft_atomic_write_tsv(
    combined$validation,
    file.path(run_paths$run_root, "precision_rescue_condition_validation.tsv")
  )
  if (!isTRUE(combined$validation$passed[[1L]])) {
    stop(
      "Independent R=200 rescue did not reach MCSE 0.01N; observed ",
      signif(combined$validation$maximum_mcse_N[[1L]], 6)
    )
  }

  patched_cache <- file.path(
    run_paths$cache, basename(selected$cache_path[[1L]])
  )
  f6pr_patch_task(
    base_cache, patched_cache, combined, condition, rescue_fingerprint,
    passage_bundle$stochastic$config$mcse_target_N
  )
  tasks <- base_tasks
  tasks$cache_path <- file.path(
    base_run_root, "task_cache", basename(tasks$cache_path)
  )
  tasks$cache_path[tasks$task_id == selected$task_id[[1L]]] <- patched_cache
  f6ft_atomic_write_tsv(
    tasks, file.path(run_paths$run_root, "full_range_task_manifest.tsv")
  )
  aggregation <- f6g_aggregate(
    tasks, run_paths, rescue_fingerprint, passage_bundle, smoke = FALSE,
    accepted_cache_fingerprints = c(base_fingerprint, rescue_fingerprint)
  )
  provenance <- data.frame(
    key = c(
      "profile", "run_id", "rescue_policy", "base_run_root",
      "base_fingerprint", "rescue_fingerprint", "model_code_root",
      "model_source_fingerprint", "condition", "replicates",
      "master_seed", "rng_independence", "n_core", "validation"
    ),
    value = c(
      f6g_profile(), run_paths$run_id,
      "immutable pilot; one failed condition; fixed-R independent production",
      base_run_root, base_fingerprint, rescue_fingerprint,
      normalizePath(paths$oxygen_code, mustWork = TRUE),
      f6r_model_source_fingerprint(paths),
      "in vitro|C02|initial=3N|O2=0.7|p_misseg=0.01",
      as.character(replicates), as.character(master_seed),
      "new L'Ecuyer-CMRG stream catalog; no pilot production draw reused",
      as.character(n_core), validation_path
    ), stringsAsFactors = FALSE
  )
  provenance_path <- f6ft_atomic_write_tsv(
    provenance,
    file.path(run_paths$run_root, "precision_rescue_run_provenance.tsv")
  )
  if (isTRUE(publish_current)) {
    f6g_publish_current(paths, run_paths, rescue_fingerprint)
  }
  unlink(vapply(results, `[[`, character(1L), "checkpoint"))
  invisible(list(
    paths = run_paths, condition = combined$validation,
    aggregation = aggregation, provenance = provenance_path
  ))
}
