#!/usr/bin/env Rscript
# Explicit, fail-closed dispatch for canonical natural-scale targets.
fixed_for <- function(task, study) {
  required <- c("case", "point", "target", "value")
  if (!all(required %in% names(task))) stop("Incomplete target dispatch")
  if (identical(task$case, "unrestricted")) {
    if (!identical(task$point, "base") || !identical(task$target, "none") || !is.na(task$value)) stop("Conflicting unrestricted target")
    return(list())
  }
  if (!identical(task$case, "fixed_target")) stop("Unknown case kind")
  registry <- read_table(file.path(study, "cases.tsv")); row <- registry[registry$point == task$point, , drop = FALSE]
  if (nrow(row) != 1L || !identical(task$target, row$target[[1]]) || length(task$value) != 1L ||
      !is.finite(task$value) || !identical(as.numeric(task$value), row$value[[1]])) stop("Conflicting canonical fixed target")
  if (!task$target %in% c("p_misseg", "buffer_beta")) stop("Unknown scientific target")
  setNames(list(task$value), task$target)
}

main <- function() {
  args <- commandArgs(TRUE); stopifnot(length(args) >= 4)
  project <- normalizePath(args[1]); study <- normalizePath(args[2]); mode <- args[3]; tid <- args[4]
  code <- file.path(study, "implementation"); modules <- file.path(code, c("helpers.R", "runtime.R"))
  for (p in modules) source(p)
  verify_snapshot(study)
  settings <- read_kv(file.path(study, "settings.tsv"))
  stopifnot(settings$authorized, settings$core_hour_cap == 4000, settings$qos == "xxlarge", settings$output_language == "English",
    identical(Sys.getenv("O2PROF_IMAGE_SHA256"), settings$image_sha256))
  if (mode == "unit") {
    source(file.path(code, "test_runtime.R"), local = environment())
    return(invisible(NULL))
  }
  if (settings$wp0_input_manifest_sha256 != o2prof_hash(file.path(project, settings$wp0_relative_path, "manifests/input_manifest.tsv"))) stop("Source manifest fingerprint mismatch")
  sources <- read_table(file.path(study, "manifests/source_manifest.tsv"))
  if (!all(vapply(file.path(project, sources$path), o2prof_hash, character(1)) == sources$sha256)) stop("Source fingerprint mismatch")
  ctx <- o2prof_load_culture(file.path(project, settings$wp0_relative_path), modules)
  if (mode == "report") { source(file.path(code, "report.R"), local = environment()); return(invisible(NULL)) }
  if (mode == "replay") {
    id <- strsplit(tid, "_rep", fixed = TRUE)[[1]][1]
    imp <- read_table(file.path(study, "imports.tsv")); row <- imp[imp$import_id == id, , drop = FALSE]
    stopifnot(nrow(row) == 1L)
    registry <- read_table(file.path(study, "cases.tsv")); point <- row$point[[1]]
    fixed <- if (point == "base") list() else setNames(list(registry$value[registry$point == point]), registry$target[registry$point == point])
    map <- original_map(ctx, fixed)
    replay_fit(ctx, map, readRDS(file.path(study, row$relative_path[[1]])), file.path(study, "tests/replays", tid))
    return(invisible(NULL))
  }
  fixture <- startsWith(tid, "fixture_")
  task <- read_kv(file.path(study, "tasks", paste0(tid, ".tsv")))
  if (mode %in% c("global", "local")) {
    expected <- if (mode == "global") 8L else 1L
    if (as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) != expected) stop("Wrong actual CPU allocation")
  }
  fixed <- fixed_for(task, study)
  if (!fixture && mode != "check") {
    gate <- read_kv(file.path(study, "tests/validation.tsv")); stopifnot(gate$pass_gate, gate$code_sha256 == o2prof_hash(file.path(study, "manifests/code_manifest.tsv")))
  }
  map <- original_map(ctx, fixed); out <- file.path(study, if (fixture) "tests/fixtures" else "fits", tid)
  task_path <- file.path(study, "tasks", paste0(tid, ".tsv"))
  if (mode == "diagnostic") {
    source(file.path(code, "diagnostics.R"))
    d <- file.path(study, "diagnostics", task$point); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    fit <- readRDS(file.path(project, task$source_fit)); fit$cfg <- ctx$cfg
    replay_fit(ctx, map, fit, file.path(d, "source_replay"))
    kernels <- new.env(parent = emptyenv()); summaries <- list(); series <- list(); rates <- list()
    scored <- ctx$manifest$observation_id[ctx$manifest$growth_scored]
    for (cohort in c("2N", "4N")) for (segment in fit$best_components[[paste0("run_", cohort)]]$segment_results) {
      key <- as.character(segment$segment$oxygen_pct)
      if (!exists(key, kernels, inherits = FALSE)) {
        kernel <- pilotdiag_rates(ctx, fit$best_params, ctx$cfg, segment$segment$oxygen_pct)
        assign(key, kernel, kernels); rates[[length(rates) + 1L]] <- kernel$rates
      }
      checked <- pilotdiag_dense_flux(ctx, fit, cohort, segment, get(key, kernels))
      checked$summary$recording_scope <- if (any(segment$segment$data_ids %in% scored)) "recorded_growth" else "other_protocol_segments"
      checked$summary$saved_time_difference <- checked$daily_replay_maximum_difference
      summaries[[length(summaries) + 1L]] <- checked$summary; series[[length(series) + 1L]] <- checked$flux
    }
    o2prof_write_tsv(do.call(rbind, summaries), file.path(d, "segment_event_ledger.tsv"))
    o2prof_write_tsv(do.call(rbind, series), file.path(d, "dense_flux.tsv"))
    o2prof_write_tsv(do.call(rbind, rates), file.path(d, "state_daughter_rates.tsv"))
    write_kv(list(pass = TRUE, source_fit = task$source_fit, source_sha256 = o2prof_hash(file.path(project, task$source_fit)),
      segments = length(summaries), utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(d, "validation.tsv"))
    assert_runtime(d); writeLines("Verified kernel, dense trajectory and event balance", file.path(d, "COMPLETE"))
    return(invisible(NULL))
  }
  if (mode == "check") {
    payload <- read_table(file.path(out, "payload_manifest.tsv"))
    stopifnot(all(vapply(file.path(out, payload$path), o2prof_hash, character(1)) == payload$sha256))
    inputs <- read_kv(file.path(out, "run_inputs.tsv")); result <- read_kv(file.path(out, "result.tsv"))
    stopifnot(inputs$task_sha256 == o2prof_hash(task_path), inputs$code_sha256 == o2prof_hash(file.path(study, "manifests/code_manifest.tsv")),
      inputs$settings_sha256 == o2prof_hash(file.path(study, "settings.tsv")))
    if (!is.na(inputs$source_fit)) stopifnot(inputs$source_sha256 == o2prof_hash(file.path(project, inputs$source_fit)))
    if (!is.na(inputs$fallback_source)) stopifnot(inputs$fallback_sha256 == o2prof_hash(file.path(project, inputs$fallback_source)))
    if (isTRUE(result$skipped)) {
      stopifnot(task$phase == "neighbor", file.exists(file.path(out, "TRANSFER_SKIPPED")))
      prepared <- prepare_local_start(ctx, map, task, project, file.path(out, "external_replay"))
      stopifnot(prepared$skipped)
      write_kv(list(pass = TRUE, result_sha256 = o2prof_hash(file.path(out, "result.tsv")), accepted_fit_sha256 = NA_character_,
        payload_manifest_sha256 = o2prof_hash(file.path(out, "payload_manifest.tsv")), fresh_process = TRUE), file.path(out, "artifact_validation.tsv"))
      writeLines("Independently verified infeasible neighbor transfer", file.path(out, "COMPLETE")); return(invisible(NULL))
    }
    if (!is.na(inputs$fallback_source)) stopifnot(inputs$fallback_sha256 == o2prof_hash(file.path(project, inputs$fallback_source)))
    audit_improvements(ctx, map, out, result$fingerprint)
    fit <- readRDS(file.path(out, "accepted_fit.rds")); replay_fit(ctx, map, fit, file.path(out, "external_replay"))
    stopifnot(abs(fit$best_components$objective - result$objective) <= 1e-8)
    if (task$mode == "global") {
      stopifnot(identical(readRDS(file.path(out, "native_de_initial_population.rds")), readRDS(file.path(study, "populations", tid, "population.rds"))),
        result$actual_generations >= task$min_generations, readRDS(file.path(out, "local_intercept_calls.rds")) == 1L)
    } else {
      prepared <- prepare_local_start(ctx, map, task, project, file.path(out, "external_start_check"))
      stopifnot(!prepared$skipped, identical(readRDS(file.path(out, "actual_local_start.rds")), prepared$fit$best_params),
        identical(as.numeric(readRDS(file.path(out, "native_sbplx_start.rds"))), as.numeric(o2review_unit_spec(map, readRDS(file.path(out, "actual_local_start.rds")))$unit)))
    }
    write_kv(list(pass = TRUE, result_sha256 = o2prof_hash(file.path(out, "result.tsv")),
      accepted_fit_sha256 = o2prof_hash(file.path(out, "accepted_fit.rds")), payload_manifest_sha256 = o2prof_hash(file.path(out, "payload_manifest.tsv")), fresh_process = TRUE), file.path(out, "artifact_validation.tsv"))
    writeLines("Complete independently validated payload", file.path(out, "COMPLETE"))
    return(invisible(NULL))
  }
  if (mode == "screen") {
    folder <- file.path(study, "populations", tid); stopifnot(!dir.exists(folder)); dir.create(folder, recursive = TRUE)
    set.seed(task$seed); saveRDS(.Random.seed, file.path(folder, "rng_initial.rds"))
    start <- as.numeric(Sys.time()); accepted <- list(); types <- character(); records <- list(); proposals <- list()
    target <- if (task$family == "B") c(continuous = 32L, boundary = 32L) else c(continuous = 64L, boundary = 0L)
    counts <- c(continuous = 0L, boundary = 0L)
    for (i in seq_len(5000L)) {
      if (as.numeric(Sys.time()) - start >= 3500 || all(counts >= target)) break
      type <- if (task$family == "B" && runif(1) < .5) "boundary" else "continuous"
      proposal <- if (i %% 2L) "full_domain" else "growth_favoring"
      p <- o2review_draw(map, proposal, i)
      if (type == "boundary") {
        boundary <- c("alpha_o2", "gamma_growth", "mu_hp", "gamma_mu", "p_mis_base", "buffer_smax")
        selected <- if (runif(1) < 1/3) boundary else sample(boundary, sample(1:3, 1))
        for (n in selected) p[[n]] <- map$table[[if (n == "buffer_smax") "upper" else "lower"]][match(n, map$table$parameter)]
      } else selected <- character()
      reason <- ""; q <- tryCatch({comp <- ctx$raw_objective(p, ctx$objects, ctx$cfg); o2prof_check_samples(comp, ctx$manifest); comp$objective}, error = function(e) {reason <<- conditionMessage(e); Inf})
      keep <- is.finite(q) && counts[type] < target[type]
      records[[i]] <- data.frame(draw = i, branch = type, proposal = proposal, boundary_set = paste(selected, collapse = ";"), objective = q, accepted = keep, error = gsub("[\t\r\n]", " ", reason))
      proposals[[i]] <- p
      if (keep) { accepted[[length(accepted) + 1L]] <- p; types <- c(types, type); counts[type] <- counts[type] + 1L }
    }
    saveRDS(proposals, file.path(folder, "proposals.rds")); o2prof_write_tsv(do.call(rbind, records), file.path(folder, "screening_draws.tsv"))
    pass <- if (task$family == "A") counts["continuous"] >= 8 else all(counts >= c(continuous = 8L, boundary = 8L))
    write_kv(list(pass = unname(pass), continuous = unname(counts[1]), boundary = unname(counts[2]), proposals = length(records), seed = task$seed), file.path(folder, "screening.tsv"))
    if (!pass) stop("Initialization unresolved; prospective origin attempt consumed")
    anchor <- accepted[[1]]; pop <- o2prof_initial_population(map, 256L, list(anchor = anchor))$population
    for (i in seq_along(accepted)) pop[i, ] <- o2prof_to_free(map, accepted[[i]])
    saveRDS(pop, file.path(folder, "population.rds")); saveRDS(anchor, file.path(folder, "anchor.rds")); saveRDS(.Random.seed, file.path(folder, "rng_final.rds"))
    o2prof_write_tsv(data.frame(member = seq_len(256), provenance = c(types, rep("full_domain_random", 256 - length(types)))), file.path(folder, "member_provenance.tsv"))
    write_kv(list(population_sha256 = o2prof_hash(file.path(folder, "population.rds")), anchor_sha256 = o2prof_hash(file.path(folder, "anchor.rds"))), file.path(folder, "population_hashes.tsv"))
    assert_runtime(folder); writeLines("Completed independent feasible-bank initialization", file.path(folder, "COMPLETE")); return(invisible(NULL))
  }
  stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
  source_fit <- if (task$mode == "local") task$source_fit else NA_character_
  write_kv(list(task_sha256 = o2prof_hash(task_path), code_sha256 = o2prof_hash(file.path(study, "manifests/code_manifest.tsv")),
    settings_sha256 = o2prof_hash(file.path(study, "settings.tsv")), allocated_cpus = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), source_fit = source_fit,
    source_sha256 = if (!is.na(source_fit)) o2prof_hash(file.path(project, source_fit)) else NA_character_,
    fallback_source = task$fallback_source, fallback_sha256 = if (!is.na(task$fallback_source)) o2prof_hash(file.path(project, task$fallback_source)) else NA_character_), file.path(out, "run_inputs.tsv"))
  fingerprint <- o2prof_hash(file.path(out, "run_inputs.tsv")); saveRDS(map, file.path(out, "parameter_map.rds")); o2prof_write_tsv(map$table, file.path(out, "parameters.tsv"))
  if (task$mode == "local") {
    prepared <- prepare_local_start(ctx, map, task, project, file.path(out, "prepared_start"))
    if (prepared$skipped) {
      writeLines("No feasible registered neighbor mixture", file.path(out, "TRANSFER_SKIPPED"))
      write_kv(list(task_id = task$task_id, point = task$point, origin = task$origin, family = task$family, mode = task$mode,
        objective = NA_real_, skipped = TRUE, local_stable = FALSE, start_objective = NA_real_, fingerprint = fingerprint), file.path(out, "result.tsv"))
      paths <- list.files(out, full.names = TRUE, recursive = TRUE)
      o2prof_write_tsv(data.frame(path = substring(paths, nchar(out) + 2L), sha256 = vapply(paths, o2prof_hash, character(1))), file.path(out, "payload_manifest.tsv"))
      return(invisible(NULL))
    }
    task$start_objective <- prepared$fit$best_components$objective
    saveRDS(prepared$fit, file.path(out, "prepared_fit.rds"))
    run_local(ctx, map, file.path(out, "prepared_fit.rds"), out, fingerprint, task, fixture)
  } else if (task$mode == "global") {
    folder <- file.path(study, "populations", tid)
    if (fixture) {
      dir.create(folder, recursive = TRUE)
      fit <- readRDS(file.path(project, task$source_fit)); anchor <- o2prof_project_start(map, fit$best_params)$params
      set.seed(task$seed); population <- o2prof_initial_population(map, 256L, list(anchor = anchor))$population
      saveRDS(population, file.path(folder, "population.rds")); saveRDS(anchor, file.path(folder, "anchor.rds"))
    } else {
      stopifnot(file.exists(file.path(folder, "COMPLETE")))
      hashes <- read_kv(file.path(folder, "population_hashes.tsv"))
      stopifnot(hashes$population_sha256 == o2prof_hash(file.path(folder, "population.rds")), hashes$anchor_sha256 == o2prof_hash(file.path(folder, "anchor.rds")))
      population <- readRDS(file.path(folder, "population.rds")); anchor <- readRDS(file.path(folder, "anchor.rds"))
    }
    ctx <- o2pilot_track_context(ctx, file.path(out, "best_candidates"), fingerprint)
    o2prof_install_culture_map(ctx, map, list(anchor = anchor), out, population, TRUE)
    intercepted <- new.env(parent = emptyenv()); intercepted$calls <- 0L; ctx$backend$optim <- o2review_local_stub(intercepted)
    cleanup <- o2prof_trace_native_population(out); on.exit(cleanup(), add = TRUE)
    control <- o2review_de_controls(task); write_kv(control, file.path(out, "execution_controls.tsv"))
    set.seed(task$seed + 17L); saveRDS(.Random.seed, file.path(out, "rng_initial.rds"))
    ctx$backend$main(list(parameter_table = ctx$cfg$parameter_table,
      fit_objects_dir = file.path(ctx$wp0_root, "source/inputs/oxygen/ploidyOxygen/data/fit_objects"),
      flow_density_path = file.path(ctx$wp0_root, "source/inputs/oxygen/data/g0g1_ploidy_density_grid.csv"),
      out_dir = out, seed = task$seed + 17L, n_cores = 8L, itermax = control$itermax,
      itermax_max = control$itermax, de_reltol = control$reltol, de_steptol = control$steptol, auto_viz = FALSE))
    saveRDS(.Random.seed, file.path(out, "rng_final.rds")); saveRDS(intercepted$calls, file.path(out, "local_intercept_calls.rds"))
    workers <- unique(read_table(file.path(out, "worker_loaded_namespaces.tsv"))$pid)
    stopifnot(length(workers) == 8L)
    live <- character()
    for (retry in seq_len(50)) {
      live <- suppressWarnings(system2("ps", c("-p", paste(workers, collapse = ","), "-o", "pid=,stat="), stdout = TRUE))
      live <- live[!grepl("Z", live)]
      if (!length(live)) break
      Sys.sleep(.1)
    }
    if (length(live)) stop("Global workers remained active after backend return")
    o2prof_write_tsv(data.frame(pid = workers, active_after_global = FALSE), file.path(out, "worker_shutdown.tsv"))
    write_kv(list(method = "Native DEoptim only", inline_local_optimization = FALSE, intercepted_evaluation_entries = intercepted$calls), file.path(out, "method_sequence.tsv"))
    native <- readRDS(file.path(out, "fit_result.rds"))
    native$backend_reported_optimizer_trace <- native$optimizer_trace
    native$optimizer_trace$method <- "DEoptim_parallel_global_only"
    native$optimizer_trace$local_attempted <- FALSE
    native$optimizer_trace$local_accepted <- FALSE
    native$optimizer_trace$local_message <- "Backend local entry intercepted for one score evaluation; no inline local optimization"
    native$local_intercept_evaluation <- native$local_optim
    native$local_optim <- NULL
    saveRDS(native, file.path(out, "fit_result.rds"))
    stopifnot(native$optimizer_trace$iter_completed >= task$min_generations, intercepted$calls == 1L)
    finalize_fit(ctx, map, out, fingerprint, task, native = native)
  } else stop("Unknown mode")
}
main()
