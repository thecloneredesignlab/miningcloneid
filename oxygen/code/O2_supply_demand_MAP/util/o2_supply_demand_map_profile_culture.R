# Isolated culture adapter. The frozen production backend owns DEoptim, worker
# preflight, local refinement, penalties, and the simulation/observation model.

o2prof_load_culture <- function(wp0_root, module_paths, verify = TRUE) {
  wp0_root <- normalizePath(wp0_root, mustWork = TRUE)
  contract <- jsonlite::read_json(file.path(wp0_root, "manifests", "production_contract.json"), simplifyVector = TRUE)
  if (!isTRUE(contract$replay$pass) || contract$production_reference$commit != "1a5d861c9c4d2af72ea0d07def0023831cc531e4") stop("Unverified production reference")
  if (verify) {
    frozen <- read.delim(file.path(wp0_root, "manifests", "input_manifest.tsv"), stringsAsFactors = FALSE)
    actual <- vapply(file.path(wp0_root, frozen$path), o2prof_hash, character(1))
    if (any(actual != frozen$sha256)) stop("Frozen production source/input fingerprint changed")
  }
  source_root <- file.path(wp0_root, "source", "code", "oxygen", "code", "O2_supply_demand_MAP")
  saved <- readRDS(file.path(wp0_root, "source", "seed144", "fit_result.rds"))
  cfg <- saved$cfg
  cfg$parameter_table <- file.path(wp0_root, "source", "seed144", "parameter_table_input.csv")
  backend <- new.env(parent = .GlobalEnv)
  # Production's bootstrap prefers --file over source frames. Mask only that
  # bootstrap lookup so it resolves its own frozen source directory correctly.
  backend$commandArgs <- function(trailingOnly = FALSE) {
    x <- base::commandArgs(trailingOnly = trailingOnly)
    if (!trailingOnly) x <- x[!startsWith(x, "--file=")]
    x
  }
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = file.path(source_root, "model"))
  source(file.path(source_root, "util", "o2_supply_demand_map_fit_invitro_backend.R"), local = backend, chdir = TRUE)
  stopifnot(identical(backend$INVITRO_PASSAGE_IMPLEMENTATION, "v2"))
  inputs <- file.path(wp0_root, "source", "inputs", "oxygen")
  objects <- backend$ivt_load_fit_objects(
    repo_root = inputs, fit_objects_dir = file.path(inputs, "ploidyOxygen", "data", "fit_objects"),
    flow_csv_path = file.path(inputs, "data", "g0g1_ploidy_density_grid.csv"))
  manifest <- read.delim(file.path(wp0_root, "manifests", "observation_manifest.tsv"), stringsAsFactors = FALSE)
  list(backend = backend, cfg = cfg, objects = objects, saved = saved, manifest = manifest,
       wp0_root = wp0_root, source_root = source_root, contract = contract,
       raw_objective = backend$ivt_objective_components, module_paths = module_paths)
}

o2prof_independent_score <- function(ctx, comp) {
  # Independent arithmetic on raw observations and predicted distributions.
  values <- list(growth = numeric(), karyotype = numeric(), flow = numeric())
  for (cohort in c("2N", "4N")) {
    run <- comp[[paste0("run_", cohort)]]
    for (s in run$segment_results) {
      p <- pmax(as.numeric(s$selection$selected_frac), 0)
      stopifnot(all(is.finite(p)), sum(p) > 0)
      p <- p / sum(p)
      K <- outer(run$grid_pre, run$grid_pre, function(a, b) stats::dnorm(a, b, comp$sigma_kary))
      colsum <- colSums(K)
      colsum[!is.finite(colsum) | colsum <= 0] <- 1
      smoothed <- as.numeric(sweep(K, 2, colsum, "/") %*% p)
      names(smoothed) <- as.character(run$grid_pre)
      for (pid in s$segment$data_ids) {
        obs <- ctx$backend$ivt_observed_passage_summary(ctx$objects$fit_data[[pid]])
        if (is.finite(obs$observed_growth)) {
          rate <- (log(s$selection$selected_live_cells) - log(s$sim$Ntot_live_obs[[1]])) / s$selection$selected_day
          sigma <- comp$sigma_growth
          values$growth <- c(values$growth, -log(sigma) - log(2 * pi) / 2 - ((obs$observed_growth - rate) / sigma)^2 / 2)
        }
        if (length(obs$observed_kary)) {
          probability <- smoothed[as.character(obs$observed_kary)]
          probability[is.na(probability)] <- 1e-12
          values$karyotype <- c(values$karyotype, mean(log(pmax(probability, 1e-12))))
        }
        f <- obs$observed_flow
        if (!is.null(f) && nrow(f)) {
          g <- as.numeric(f$ploidy)
          density <- as.numeric(f$density)
          stopifnot(all(is.finite(g)), all(is.finite(density)), all(density >= 0), sum(density) > 0)
          step <- if (length(g) > 1L) stats::median(diff(sort(unique(g)))) else 1
          mass <- density * step
          mass <- mass / sum(mass)
          predicted <- vapply(g, function(x) sum(p * stats::dnorm(x, run$grid_pre / ctx$cfg$N_UNIT, comp$sigma_flow_ploidy)), numeric(1))
          predicted <- predicted / (sum(predicted) * step)
          values$flow <- c(values$flow, sum(mass * log(pmax(predicted, 1e-12))))
        }
      }
    }
  }
  expected <- c(growth = sum(ctx$manifest$growth_scored), karyotype = sum(ctx$manifest$karyotype_scored), flow = sum(ctx$manifest$flow_scored))
  if (!identical(vapply(values, length, integer(1)), as.integer(expected) |> setNames(names(expected)))) stop("Independent scoring sample mismatch")
  if (any(!is.finite(unlist(values)))) stop("Nonfinite independently recomputed score")
  means <- vapply(values, mean, numeric(1))
  list(objective = -sum(means), means = means, counts = expected)
}

o2prof_install_culture_map <- function(ctx, map, starts, out_dir = NULL, provided_population = NULL,
                                      count_evaluations = FALSE) {
  backend <- ctx$backend
  projected <- lapply(starts, function(x) o2prof_project_start(map, x))
  first <- projected[[1]]$params
  backend$ivt_optimizer_spec <- function(cfg) {
    if (!identical(cfg, ctx$cfg)) stop("Immutable culture configuration changed")
    o2prof_spec(map, first)
  }
  backend$build_invitro_cfg <- function(parameter_table, dt = 0.05, init_total_size = 1e6, o2_upper_bound = 21, fixed_oxygen = TRUE) {
    if (normalizePath(parameter_table) != normalizePath(ctx$cfg$parameter_table) || dt != ctx$cfg$DT ||
        init_total_size != ctx$cfg$init_total_size || o2_upper_bound != ctx$cfg$o2_S0_upper_bound || !isTRUE(fixed_oxygen)) stop("Attempt to alter frozen model settings")
    ctx$cfg
  }
  backend$ivt_optim_par_to_run_params <- function(par_t, cfg) {
    params <- o2prof_from_free(map, par_t)
    normalized <- backend$normalize_run_params_common(params, cfg = cfg)
    o2prof_to_free(map, normalized)
    normalized
  }
  backend$ivt_run_params_to_optim_par <- function(run_params, cfg) o2prof_to_free(map, run_params)
  backend$ivt_load_fit_objects_compat <- function(fit_objects_dir, flow_density_path = NULL) ctx$objects
  state <- new.env(parent = emptyenv())
  state$pid <- NA_integer_
  state$count <- 0L
  raw <- ctx$raw_objective
  backend$ivt_objective_components <- function(run_params, fit_objects, cfg, ...) {
    if (!identical(Sys.getpid(), state$pid)) { state$pid <- Sys.getpid(); state$count <- 0L }
    state$count <- state$count + 1L
    start <- proc.time()[[3]]
    error_text <- ""
    value <- NA_real_
    on.exit({
      if (count_evaluations) {
        path <- file.path(out_dir, "evaluation_logs", paste0("pid_", Sys.getpid(), ".tsv"))
        row <- data.frame(pid = Sys.getpid(), evaluation = state$count, phase = getOption("o2prof.phase", "preflight"),
                          seconds = proc.time()[[3]] - start, objective = value,
                          error = gsub("[\t\r\n]", " ", error_text))
        write.table(row, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = !file.exists(path), append = file.exists(path))
      }
    }, add = TRUE)
    tryCatch({
      if (!identical(cfg, ctx$cfg)) stop("Immutable culture configuration changed during evaluation")
      o2prof_to_free(map, run_params)
      comp <- raw(run_params = run_params, fit_objects = fit_objects, cfg = cfg, ...)
      o2prof_check_samples(comp, ctx$manifest)
      value <- comp$objective
      comp
    }, error = function(e) { error_text <<- conditionMessage(e); stop(e) })
  }
  backend$build_invitro_de_initial_population <- function(NP, lower, upper, init) {
    answer <- o2prof_initial_population(map, NP, starts, provided_population)
    if (!identical(unname(answer$population[1, ]), as.numeric(init))) stop("Actual initial vector was not passed through native backend")
    if (!is.null(out_dir)) {
      saveRDS(answer$population, file.path(out_dir, "requested_initial_population.rds"))
      o2prof_write_tsv(as.data.frame(answer$population), file.path(out_dir, "actual_initial_population_free.tsv"))
      natural <- do.call(rbind, lapply(seq_len(nrow(answer$population)), function(i) {
        p <- o2prof_from_free(map, setNames(answer$population[i, ], colnames(answer$population)))
        as.data.frame(as.list(unlist(p[o2prof_symbols])))
      }))
      natural$member <- seq_len(nrow(natural))
      natural$source <- c(names(starts), rep("independent_uniform_transformed", nrow(natural) - length(starts)))
      o2prof_write_tsv(natural, file.path(out_dir, "actual_initial_population_natural.tsv"))
      for (i in seq_along(answer$starts)) o2prof_write_tsv(answer$starts[[i]]$projection, file.path(out_dir, paste0("start_projection_", i, ".tsv")))
    }
    answer$population
  }
  original_initializer <- backend$initialize_invitro_deoptim_workers
  backend$initialize_invitro_deoptim_workers <- function(cluster, cpp_info) {
    answer <- original_initializer(cluster, cpp_info)
    parallel::clusterCall(cluster, function(paths) {
      for (path in paths) sys.source(path, envir = .GlobalEnv)
      options(o2prof.phase = "preflight")
      TRUE
    }, ctx$module_paths)
    answer
  }
  backend$optim <- function(...) {
    old <- options(o2prof.phase = "local")
    on.exit(options(old), add = TRUE)
    stats::optim(...)
  }
  if (!is.null(out_dir)) dir.create(file.path(out_dir, "evaluation_logs"), recursive = TRUE, showWarnings = FALSE)
  invisible(list(map = map, projected_starts = projected, backend = backend))
}

o2prof_trace_native_population <- function(out_dir) {
  # Observe the actual DEoptim entry, not merely the population generator.
  recorder <- function(control) {
    actual <- control$initialpop
    expected <- readRDS(file.path(out_dir, "requested_initial_population.rds"))
    if (!identical(actual, expected)) stop("Native DEoptim received the wrong population")
    saveRDS(actual, file.path(out_dir, "native_de_initial_population.rds"))
    options(o2prof.phase = "de")
    if (!is.null(control$cluster)) parallel::clusterCall(control$cluster, function() { options(o2prof.phase = "de"); TRUE })
  }
  previous <- getOption("o2prof.de_recorder")
  options(o2prof.de_recorder = recorder)
  trace("DEoptim", where = asNamespace("DEoptim"), print = FALSE,
        tracer = quote(getOption("o2prof.de_recorder")(control)),
        exit = quote(options(o2prof.phase = "post_de")))
  function() {
    untrace("DEoptim", where = asNamespace("DEoptim"))
    options(o2prof.de_recorder = previous)
  }
}
