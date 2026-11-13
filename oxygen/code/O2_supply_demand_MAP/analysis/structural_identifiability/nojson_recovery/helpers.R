# Shared natural/free parameter mapping and validation for diagnostic profiles.
# This module does not source a simulator or implement an optimizer.

o2prof_symbols <- c(
  "lam_max", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta", "buffer_n_exp",
  "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu", "O2_crit", "n_O",
  "p_mis_base", "sigma_growth", "sigma_kary", "init_mean_2N", "init_sd_2N",
  "init_mean_4N", "init_sd_4N")

o2prof_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

o2prof_hash <- function(path) unname(digest::digest(file = path, algo = "sha256"))

o2prof_map <- function(parameter_table, base_params, domain = c("production", "zero_extended"),
                       fixed = list(), penalty_active = character(),
                       zero_scales = c(p_misseg = 0.015, buffer_beta = 1)) {
  domain <- match.arg(domain)
  tab <- if (is.data.frame(parameter_table)) parameter_table else read.csv(parameter_table, stringsAsFactors = FALSE)
  required <- c("param_symbol", "init_value", "lower_bound", "upper_bound", "use_invitro_fit")
  if (!all(required %in% names(tab))) stop("Incomplete production parameter table")
  flagged <- tolower(trimws(as.character(tab$use_invitro_fit))) %in% c("true", "t", "1", "yes", "y")
  if (anyDuplicated(tab$param_symbol) || !setequal(tab$param_symbol[flagged], o2prof_symbols)) stop("Unexpected production fitted-parameter specification")
  if (!all(o2prof_symbols %in% names(base_params))) stop("Missing full-precision baseline parameters")
  if (length(fixed) && (is.null(names(fixed)) || anyDuplicated(names(fixed)) || any(!names(fixed) %in% o2prof_symbols))) stop("Unknown or unnamed fixed parameter")
  if (any(!penalty_active %in% o2prof_symbols)) stop("Unknown penalty-active parameter")
  if (!setequal(names(zero_scales), c("p_misseg", "buffer_beta")) || any(!is.finite(zero_scales) | zero_scales <= 0)) stop("Zero-transform scales must be fixed positive constants")
  p <- tab[match(o2prof_symbols, tab$param_symbol), , drop = FALSE]
  rownames(p) <- NULL
  linear <- c("buffer_smax", "gamma_growth", "gamma_mu", "n_O", "init_mean_2N", "init_mean_4N")
  p$parameter <- p$param_symbol
  p$lower <- as.numeric(p$lower_bound)
  p$upper <- as.numeric(p$upper_bound)
  p$transform <- ifelse(p$parameter %in% linear, "identity", "log10")
  p$scale <- 1
  if (domain == "zero_extended") {
    hit <- match(names(zero_scales), p$parameter)
    p$lower[hit] <- 0
    p$transform[hit] <- "log1p"
    p$scale[hit] <- zero_scales
  }
  if (any(!is.finite(p$lower) | !is.finite(p$upper) | p$lower >= p$upper) || any(p$lower[p$transform == "log10"] <= 0)) stop("Invalid parameter domain")
  p$status <- "active"
  p$reason <- "fitted nuisance"
  p$data_active <- TRUE
  p$penalty_active <- p$parameter %in% penalty_active
  p$fixed_value <- NA_real_
  for (n in names(fixed)) {
    i <- match(n, p$parameter)
    value <- as.numeric(fixed[[n]])
    if (length(value) != 1L || !is.finite(value) || value < p$lower[i] || value > p$upper[i]) stop("Fixed value outside declared domain: ", n)
    p$status[i] <- "fixed"
    p$fixed_value[i] <- value
    p$reason[i] <- "exact target constraint"
  }
  inactive <- character()
  if ("p_misseg" %in% names(fixed) && fixed$p_misseg == 0) inactive <- c(inactive, k_o_mis = "zero inducible amplitude")
  if ("buffer_beta" %in% names(fixed) && fixed$buffer_beta == 0) inactive <- c(inactive, buffer_n_exp = "zero ploidy buffering strength")
  for (n in names(inactive)) {
    i <- match(n, p$parameter)
    p$data_active[i] <- FALSE
    if (p$status[i] == "active" && !p$penalty_active[i]) {
      p$status[i] <- "inactive"
      p$fixed_value[i] <- as.numeric(base_params[[n]])
      p$reason[i] <- paste(inactive[[n]], "and no complete-objective contribution")
    } else if (p$penalty_active[i]) p$reason[i] <- "data-inactive but retained for objective penalty"
  }
  prefix <- ifelse(p$transform == "identity", "", paste0(p$transform, "_"))
  p$coordinate <- paste0(prefix, p$parameter)
  map <- list(table = p[, c("parameter", "status", "reason", "data_active", "penalty_active", "lower", "upper", "transform", "scale", "coordinate", "fixed_value")],
              base_params = base_params, domain = domain, fixed = fixed, zero_scales = zero_scales)
  class(map) <- "o2prof_parameter_map"
  map
}

o2prof_transform <- function(x, table, inverse = FALSE) {
  out <- as.numeric(x)
  log <- table$transform == "log10"
  zero <- table$transform == "log1p"
  if (inverse) {
    out[log] <- 10^out[log]
    out[zero] <- table$scale[zero] * expm1(out[zero])
  } else {
    if (any(out[log] <= 0) || any(out[zero] < 0)) stop("Invalid natural value for transform")
    out[log] <- log10(out[log])
    out[zero] <- log1p(out[zero] / table$scale[zero])
  }
  out
}

o2prof_project_start <- function(map, candidate, clip = FALSE) {
  p <- map$table
  if (!all(p$parameter %in% names(candidate))) stop("Start must contain all fitted natural parameters")
  values <- vapply(p$parameter, function(n) as.numeric(candidate[[n]]), numeric(1))
  if (any(!is.finite(values))) stop("Nonfinite start")
  before <- values
  constrained <- p$status != "active"
  values[constrained] <- p$fixed_value[constrained]
  outside <- values < p$lower | values > p$upper
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(p$lower), abs(p$upper))
  materially_outside <- values < p$lower - tolerance | values > p$upper + tolerance
  if (any(materially_outside) && !clip) stop("Start outside declared bounds")
  if (any(outside)) values <- pmin(pmax(values, p$lower), p$upper)
  full <- map$base_params
  for (i in seq_along(values)) full[[p$parameter[i]]] <- unname(values[i])
  reason <- ifelse(constrained, p$status, ifelse(materially_outside, "projected_to_bound", ifelse(outside, "roundoff_to_bound", "unchanged")))
  list(params = full, projection = data.frame(parameter = p$parameter, provided_value = before, actual_value = values, reason = reason))
}

o2prof_to_free <- function(map, params) {
  p <- map$table
  v <- vapply(p$parameter, function(n) as.numeric(params[[n]]), numeric(1))
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(p$lower), abs(p$upper))
  if (length(v) != nrow(p) || any(!is.finite(v)) || any(v < p$lower - tolerance | v > p$upper + tolerance)) stop("Infeasible full parameter vector")
  fixed <- p$status != "active"
  if (any(v[fixed] != p$fixed_value[fixed])) stop("Full vector disagrees with fixed/inactive values")
  active <- !fixed
  v[active] <- pmin(pmax(v[active], p$lower[active]), p$upper[active])
  setNames(o2prof_transform(v[active], p[active, , drop = FALSE]), p$coordinate[active])
}

o2prof_from_free <- function(map, free) {
  p <- map$table
  active <- p$status == "active"
  wanted <- p$coordinate[active]
  if (length(free) != length(wanted) || (length(names(free)) && !identical(names(free), wanted))) stop("Wrong free-vector coordinates")
  v <- p$fixed_value
  v[active] <- o2prof_transform(free, p[active, , drop = FALSE], inverse = TRUE)
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(p$lower), abs(p$upper))
  if (any(!is.finite(v)) || any(v < p$lower - tolerance | v > p$upper + tolerance)) stop("Free vector outside declared bounds")
  # Roundoff at transformed endpoints is projected only onto those exact bounds.
  v[active] <- pmin(pmax(v[active], p$lower[active]), p$upper[active])
  full <- map$base_params
  for (i in seq_along(v)) full[[p$parameter[i]]] <- unname(v[i])
  full
}

o2prof_spec <- function(map, initial_params) {
  p <- map$table[map$table$status == "active", , drop = FALSE]
  list(param_name = p$coordinate, lower = o2prof_transform(p$lower, p),
       upper = o2prof_transform(p$upper, p), init = unname(o2prof_to_free(map, initial_params)))
}

o2prof_initial_population <- function(map, NP, starts, provided = NULL) {
  if (!length(starts) || is.null(names(starts)) || any(!nzchar(names(starts))) || anyDuplicated(names(starts))) stop("Provide uniquely named start vectors")
  projected <- lapply(starts, function(x) o2prof_project_start(map, x))
  vectors <- lapply(projected, function(x) o2prof_to_free(map, x$params))
  spec <- o2prof_spec(map, projected[[1]]$params)
  d <- length(spec$param_name)
  if (NP < length(starts) || NP < 5L || d < 1L) stop("Invalid initial population size")
  if (is.null(provided)) {
    population <- matrix(stats::runif(NP * d), nrow = NP, ncol = d)
    population <- sweep(sweep(population, 2, spec$upper - spec$lower, "*"), 2, spec$lower, "+")
    for (i in seq_along(vectors)) population[i, ] <- vectors[[i]]
  } else {
    population <- as.matrix(provided)
    if (!identical(dim(population), c(as.integer(NP), as.integer(d)))) stop("Wrong provided population dimensions")
    if (!identical(colnames(population), spec$param_name)) stop("Wrong provided population coordinates")
    for (i in seq_along(vectors)) if (!identical(as.numeric(population[i, ]), unname(vectors[[i]]))) stop("Provided population disagrees with declared injected starts")
  }
  colnames(population) <- spec$param_name
  if (any(!is.finite(population)) || any(sweep(population, 2, spec$lower, "<")) || any(sweep(population, 2, spec$upper, ">"))) stop("Infeasible initial population")
  list(population = population, starts = projected, labels = names(starts))
}

o2prof_check_samples <- function(comp, manifest) {
  expected <- list(growth_df = manifest$observation_id[manifest$growth_scored],
                   ploidy_df = manifest$observation_id[manifest$karyotype_scored],
                   flow_df = manifest$observation_id[manifest$flow_scored])
  for (n in names(expected)) {
    tab <- comp[[n]]
    column <- if (n == "growth_df") "loglik" else "mean_loglik"
    if (is.null(tab) || !all(c("passage_id", column) %in% names(tab)) ||
        anyDuplicated(tab$passage_id) || !setequal(as.character(tab$passage_id), expected[[n]]) ||
        any(!is.finite(tab[[column]]))) stop("Retained sample contract violated: ", n)
  }
  if (!is.finite(comp$objective) || comp$n_growth_missing_pred != 0L) stop("Nonfinite objective or missing prediction")
  invisible(TRUE)
}

o2prof_accept_candidate <- function(candidate, expected_fingerprint, map, initial_population_hash,
                                    score_tolerance = 1e-8) {
  if (!identical(candidate$fingerprint, expected_fingerprint)) stop("Stale or wrong-run fingerprint")
  if (!isTRUE(candidate$launch_completed) || !identical(candidate$launch_exit_code, 0L)) stop("Interrupted or failed launch")
  if (!identical(candidate$initial_population_hash, initial_population_hash)) stop("Wrong actual starting population")
  o2prof_to_free(map, candidate$params)
  if (!isTRUE(candidate$samples_valid) || !is.finite(candidate$objective) || !is.finite(candidate$reevaluated_objective) ||
      abs(candidate$objective - candidate$reevaluated_objective) > score_tolerance) stop("Failed independent candidate validation")
  list(accepted = TRUE, search_status = if (isTRUE(candidate$budget_limited)) "budget_limited_provisional" else "completed_search_candidate")
}


# Isolated culture adapter. The frozen production backend owns DEoptim, worker
# preflight, local refinement, penalties, and the simulation/observation model.

o2prof_load_culture <- function(wp0_root, module_paths, verify = TRUE) {
  wp0_root <- normalizePath(wp0_root, mustWork = TRUE)
  replay <- read.delim(file.path(wp0_root, "audit/baseline_replay.tsv"), stringsAsFactors = FALSE)
  stopifnot(all(is.finite(replay$difference)), max(abs(replay$difference)) < 1e-8)
  contract <- list(source_commit = "1a5d861c9c4d2af72ea0d07def0023831cc531e4", replay_verified = TRUE)
  if (verify) {
    frozen <- read.delim(file.path(wp0_root, "manifests", "input_manifest.tsv"), stringsAsFactors = FALSE)
    # Archived environment metadata is not an execution dependency.
    frozen <- frozen[!grepl("\\.(json|jsonl|ndjson)$", frozen$path, ignore.case = TRUE), ]
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
    options(o2prof.evaluation_id = state$count)
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
      natural$source <- c(names(starts), rep(if (is.null(provided_population)) "independent_uniform_transformed" else "provided_population_member", nrow(natural) - length(starts)))
      o2prof_write_tsv(natural, file.path(out_dir, "actual_initial_population_natural.tsv"))
      for (i in seq_along(answer$starts)) o2prof_write_tsv(answer$starts[[i]]$projection, file.path(out_dir, paste0("start_projection_", i, ".tsv")))
    }
    answer$population
  }
  original_penalty <- backend$make_penalty_components
  backend$make_penalty_components <- function(objective = 1e9, reason = "penalty") {
    if (grepl("checkpoint_write_failed|Stale candidate|Immutable|Wrong free-vector", reason)) stop(reason)
    original_penalty(objective = objective, reason = reason)
  }
  original_initializer <- backend$initialize_invitro_deoptim_workers
  backend$initialize_invitro_deoptim_workers <- function(cluster, cpp_info) {
    answer <- original_initializer(cluster, cpp_info)
    parallel::clusterCall(cluster, function(paths) {
      for (path in paths) sys.source(path, envir = .GlobalEnv)
      options(o2prof.phase = "preflight")
      TRUE
      forbidden <- intersect(loadedNamespaces(), c("rJava", "V8", "jsonlite", "rjson", "RJSONIO"))
      if (length(forbidden)) stop("Forbidden worker namespace")
      data.frame(pid = Sys.getpid(), namespace = loadedNamespaces())
    }, ctx$module_paths) -> loaded
    if (!is.null(out_dir)) o2prof_write_tsv(do.call(rbind, loaded), file.path(out_dir, "worker_loaded_namespaces.tsv"))
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

# Pilot bookkeeping only: the frozen production backend remains the fitter.

o2pilot_atomic_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".candidate_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary)
  if (!file.rename(temporary, path)) stop("checkpoint_write_failed: atomic rename")
  invisible(path)
}

o2pilot_checkpoint <- function(params, objective, directory, fingerprint,
    phase = getOption("o2prof.phase", "unknown"), components = NULL,
    evaluation_index = getOption("o2prof.evaluation_id", NA_integer_)) {
  if (!is.finite(objective)) stop("Nonfinite candidate cannot be checkpointed")
  path <- file.path(directory, paste0("pid_", Sys.getpid(), ".rds"))
  if (!exists(".o2pilot_best_cache", envir = .GlobalEnv, inherits = FALSE))
    assign(".o2pilot_best_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
  cache <- get(".o2pilot_best_cache", envir = .GlobalEnv)
  key <- paste(fingerprint, path, sep = "|")
  if (exists(key, envir = cache, inherits = FALSE) && get(key, envir = cache) <= objective) return(invisible(FALSE))
  previous <- if (file.exists(path)) readRDS(path) else NULL
  if (!is.null(previous)) {
    if (!identical(previous$fingerprint, fingerprint)) stop("Stale candidate checkpoint")
    assign(key, previous$objective, envir = cache)
    if (previous$objective <= objective) return(invisible(FALSE))
  }
  index <- if (is.null(previous$record_index)) 1L else previous$record_index + 1L
  selected <- NULL; scores <- NULL; component_means <- NULL
  if (!is.null(components)) {
    scores <- components[c("growth_df", "ploidy_df", "flow_df")]
    component_means <- c(growth = mean(components$growth_df$loglik),
      karyotype = mean(components$ploidy_df$mean_loglik), flow = mean(components$flow_df$mean_loglik))
    if (any(!is.finite(component_means)) || abs(-sum(component_means) - objective) > 1e-8) stop("checkpoint_write_failed: component arithmetic")
    selected <- do.call(rbind, lapply(c("2N", "4N"), function(cohort) {
      z <- components[[paste0("run_", cohort)]]$segment_results
      data.frame(cohort = cohort, segment_order = seq_along(z),
        segment_id = vapply(z, function(v) v$segment$segment_id, character(1)),
        selected_day = vapply(z, function(v) v$selection$selected_day, numeric(1)))
    }))
    if (any(!is.finite(selected$selected_day))) stop("checkpoint_write_failed: selected days")
  }
  id <- paste0("pid_", Sys.getpid(), "_", sprintf("%08d", index))
  relative <- file.path("improvements", paste0(id, ".rds"))
  archive <- file.path(directory, relative)
  if (file.exists(archive)) stop("checkpoint_write_failed: existing immutable improvement")
  record <- list(params = params, objective = objective, phase = phase, fingerprint = fingerprint,
    pid = Sys.getpid(), record_index = index, evaluation_index = evaluation_index,
    recorded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    complete_components = !is.null(components), component_means = component_means,
    score_tables = scores, selected_days = selected)
  tryCatch({
    # Publish the immutable improvement before replacing the best pointer.
    o2pilot_atomic_rds(record, archive)
    pointer <- record[c("params", "objective", "phase", "fingerprint", "pid", "record_index", "evaluation_index", "recorded_utc")]
    pointer$archive_relative_path <- relative; pointer$archive_sha256 <- o2prof_hash(archive)
    o2pilot_atomic_rds(pointer, path)
  }, error = function(e) stop("checkpoint_write_failed: ", conditionMessage(e)))
  assign(key, objective, envir = cache)
  invisible(TRUE)
}

o2pilot_track_context <- function(ctx, directory, fingerprint) {
  raw <- ctx$raw_objective
  manifest <- ctx$manifest
  ctx$raw_objective <- function(run_params, fit_objects, cfg, ...) {
    comp <- raw(run_params = run_params, fit_objects = fit_objects, cfg = cfg, ...)
    o2prof_check_samples(comp, manifest)
    o2pilot_checkpoint(run_params, comp$objective, directory, fingerprint, components = comp)
    comp
  }
  ctx
}


o2review_unit_spec <- function(map, params) {
  s <- o2prof_spec(map, params)
  s$span <- s$upper - s$lower
  s$unit <- unname((o2prof_to_free(map, params) - s$lower) / s$span)
  s
}

o2review_directions <- function(map) {
  tab <- map$table[map$table$status == "active", ]
  d <- nrow(tab)
  directions <- diag(d)
  labels <- tab$parameter
  pairs <- list(c("p_misseg", "k_o_mis"), c("p_misseg", "buffer_beta"), c("p_misseg", "mu_hp"), c("p_misseg", "gamma_mu"),
                c("p_mis_base", "buffer_smax"), c("p_mis_base", "buffer_beta"), c("p_mis_base", "buffer_n_exp"),
                c("buffer_beta", "buffer_n_exp"), c("buffer_beta", "buffer_smax"),
                c("lam_max", "mu_hp"), c("mu_hp", "gamma_mu"),
                c("O2_crit", "n_O"), c("init_mean_2N", "init_sd_2N"),
                c("init_mean_4N", "init_sd_4N"))
  for (pair in pairs) if (all(pair %in% tab$parameter)) for (sign in c(-1, 1)) {
    v <- numeric(d); v[match(pair, tab$parameter)] <- c(1, sign)
    directions <- rbind(directions, v)
    labels <- c(labels, paste(pair, collapse = if (sign == 1) "+" else "-"))
  }
  if (all(c("buffer_smax", "buffer_beta", "buffer_n_exp") %in% tab$parameter)) {
    for (sign in c(-1, 1)) {
      v <- numeric(d); v[match(c("buffer_smax", "buffer_beta", "buffer_n_exp"), tab$parameter)] <- c(1, sign, -sign)
      directions <- rbind(directions, v); labels <- c(labels, paste0("joint_survival_", sign))
    }
  }
  rownames(directions) <- labels
  directions
}

o2review_pattern <- function(fn, initial, directions, max_evaluations = 20000L,
                             deadline = Inf, mesh = c(.02, .005, .001, .0001),
                             clock = function() as.numeric(Sys.time())) {
  evaluations <- 0L; exhausted <- FALSE; history <- list(); sweeps <- list()
  evaluate <- function(x) {
    if (evaluations >= max_evaluations || clock() >= deadline) {
      exhausted <<- TRUE; return(Inf)
    }
    evaluations <<- evaluations + 1L
    value <- fn(x)
    if (length(value) != 1L || !is.finite(value)) value <- Inf
    value
  }
  x <- pmin(pmax(initial, 0), 1); best <- evaluate(x)
  if (!is.finite(best)) stop("Infeasible local start")
  # Sequential best-improvement sweeps, with an explicit final stability sweep.
  for (h in mesh) {
    repeat {
      before <- best; max_improvement <- 0; complete <- TRUE
      candidate <- x; candidate_q <- best
      for (i in seq_len(nrow(directions))) for (sign in c(-1, 1)) {
        z <- pmin(pmax(x + sign * h * directions[i, ], 0), 1)
        if (identical(z, x)) next
        q <- evaluate(z)
        if (exhausted) { complete <- FALSE; break }
        max_improvement <- max(max_improvement, best - q)
        if (q < candidate_q) { candidate <- z; candidate_q <- q }
      }
      improvement <- best - candidate_q
      if (candidate_q < best) { x <- candidate; best <- candidate_q }
      sweeps[[length(sweeps) + 1L]] <- data.frame(mesh = h, complete = complete,
          before = before, after = best, improvement = improvement,
          maximum_tested_improvement = max_improvement, evaluations = evaluations)
      history[[length(history) + 1L]] <- x
      if (exhausted || improvement < 1e-5) break
    }
    if (exhausted) break
  }
  audit <- do.call(rbind, sweeps)
  fine <- lapply(tail(mesh, 2L), function(h) {
    z <- audit[audit$mesh == h, , drop = FALSE]
    if (!nrow(z)) return(FALSE)
    z <- z[nrow(z), ]
    isTRUE(z$complete) && z$maximum_tested_improvement <= 1e-4 && z$improvement < 1e-5
  })
  list(par = x, value = best, evaluations = evaluations, sweeps = audit,
       tested_directional_stability = !exhausted && all(unlist(fine)),
       budget_limited = exhausted, history = history)
}

o2review_local_stub <- function(counter) {
  function(par, fn, ...) {
    counter$calls <- counter$calls + 1L
    list(par = par, value = fn(par), convergence = NA_integer_,
         message = "Local optimization deferred to a separately allocated job")
  }
}

o2review_de_controls <- function(task) {
  # Native DE needs steptol >= the minimum to prohibit its early-stop path.
  # This deliberately strengthens the nominal 150-generation patience.
  list(itermax = as.integer(task$max_generations), reltol = 1e-7,
       steptol = as.integer(max(task$min_generations, 150L)))
}

o2review_draw <- function(map, proposal, draw) {
  s <- o2prof_spec(map, o2prof_project_start(map, map$base_params)$params)
  z <- setNames(runif(length(s$lower), s$lower, s$upper), s$param_name)
  p <- o2prof_from_free(map, z)
  if (proposal == "growth_favoring") {
    ranges <- list(lam_max = c(1, 1.5), alpha_o2 = c(.5, .75), gamma_growth = c(1.5, 2),
                   mu_hp = c(.005, .01), gamma_mu = c(1.2, 1.6), buffer_smax = c(.8, 1))
    for (name in names(ranges)) if (map$table$status[match(name, map$table$parameter)] == "active")
      p[[name]] <- runif(1, ranges[[name]][1], ranges[[name]][2])
    if (map$table$status[match("p_misseg", map$table$parameter)] == "active") {
      # Alternate a feasible low-error proposal with full-range induction.
      if (draw %% 4L == 0L) p$p_misseg <- runif(1, 0, .003)
    }
    if (map$table$status[match("p_mis_base", map$table$parameter)] == "active" && draw %% 4L == 0L)
      p$p_mis_base <- 10^runif(1, -6, -4)
  }
  o2prof_project_start(map, p)$params
}

o2review_evaluator <- function(ctx, map, out, fingerprint, phase, deadline = Inf) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  spec <- o2prof_spec(map, o2prof_project_start(map, map$base_params)$params)
  counts <- new.env(parent = emptyenv()); counts$n <- 0L; counts$failed <- 0L
  fn <- function(unit) {
    if (as.numeric(Sys.time()) >= deadline) return(Inf)
    counts$n <- counts$n + 1L
    stage <- getOption("o2prof.local_stage", phase)
    started <- proc.time()[[3]]; reason <- ""; q <- Inf; days <- ""
    tryCatch({
      free <- setNames(spec$lower + pmin(pmax(unit, 0), 1) * (spec$upper - spec$lower), spec$param_name)
      params <- o2prof_from_free(map, free)
      params <- ctx$backend$normalize_run_params_common(params, cfg = ctx$cfg)
      o2prof_to_free(map, params)
      comp <- ctx$raw_objective(params, ctx$objects, ctx$cfg)
      o2prof_check_samples(comp, ctx$manifest); q <- comp$objective
      day_values <- unlist(lapply(c("2N", "4N"), function(cohort)
        vapply(comp[[paste0("run_", cohort)]]$segment_results, function(z) z$selection$selected_day, numeric(1))))
      days <- digest::digest(day_values, algo = "sha256")
      day_path <- file.path(out, "selected_day_vectors", paste0(days, ".rds"))
      if (!file.exists(day_path)) o2pilot_atomic_rds(day_values, day_path)
      o2pilot_checkpoint(params, q, file.path(out, "best_candidates"), fingerprint, stage, comp, counts$n)
    }, error = function(e) {
      reason <<- conditionMessage(e); counts$failed <- counts$failed + 1L
      if (grepl("checkpoint_write_failed|Stale candidate|Immutable|Wrong free-vector", reason)) stop(e)
    })
    row <- data.frame(evaluation = counts$n, phase = stage, objective = q,
        seconds = proc.time()[[3]] - started, error = gsub("[\t\r\n]", " ", reason), selected_day_vector_sha256 = days)
    path <- file.path(out, "local_evaluations.tsv")
    write.table(row, path, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = !file.exists(path), append = file.exists(path))
    q
  }
  list(fn = fn, counts = counts, spec = spec)
}

