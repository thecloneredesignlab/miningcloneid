# Typed, atomic artifacts for the isolated seven-point culture study.
o2prof_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  old <- options(digits = 17, scipen = 0); on.exit(options(old))
  x[] <- lapply(x, function(v) if (is.numeric(v) && !is.integer(v)) {
    z <- sprintf("%.17g", v); z[is.na(v)] <- NA_character_; z
  } else v)
  tmp <- tempfile(".table_", dirname(path)); on.exit(unlink(tmp), add = TRUE)
  write.table(x, tmp, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA", qmethod = "double")
  if (!file.rename(tmp, path)) stop("Atomic table publication failed")
}
read_table <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA")
write_kv <- function(values, path) {
  stopifnot(!is.null(names(values)), !anyDuplicated(names(values)))
  types <- vapply(values, function(x) if (is.logical(x)) "boolean" else if (is.integer(x)) "integer" else if (is.numeric(x)) "double" else "string", character(1))
  vals <- vapply(values, function(x) {
    stopifnot(length(x) == 1L)
    if (is.na(x)) return(NA_character_)
    if (is.logical(x)) return(if (x) "true" else "false")
    if (is.numeric(x)) return(sprintf("%.17g", x))
    as.character(x)
  }, character(1))
  o2prof_write_tsv(data.frame(key = names(values), type = types, value = vals), path)
}
read_kv <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")
  stopifnot(identical(names(x), c("key", "type", "value")), !anyDuplicated(x$key))
  ans <- lapply(seq_len(nrow(x)), function(i) switch(x$type[i],
    boolean = {stopifnot(is.na(x$value[i]) || x$value[i] %in% c("true", "false")); as.logical(x$value[i])},
    integer = as.integer(x$value[i]), double = as.numeric(x$value[i]), string = x$value[i], stop("Unknown scalar type")))
  setNames(ans, x$key)
}
verify_snapshot <- function(study) {
  controls <- read_table(file.path(study, "manifests/control_manifest.tsv"))
  for (i in seq_len(nrow(controls))) if (o2prof_hash(file.path(study, controls$path[i])) != controls$sha256[i]) stop("Stale study control")
  m <- read_table(file.path(study, "manifests/code_manifest.tsv"))
  stopifnot(nrow(m) > 0)
  for (i in seq_len(nrow(m))) if (o2prof_hash(file.path(study, "implementation", m$path[i])) != m$sha256[i]) stop("Stale implementation")
}
assert_runtime <- function(out) {
  forbidden <- intersect(loadedNamespaces(), c("rJava", "V8", "jsonlite", "rjson", "RJSONIO"))
  if (length(forbidden)) stop("Forbidden namespace loaded: ", paste(forbidden, collapse = ","))
  o2prof_write_tsv(data.frame(namespace = sort(loadedNamespaces())), file.path(out, "loaded_namespaces.tsv"))
}
audit_native_strings <- function(x) {
  if (is.character(x)) {
    z <- x[!is.na(x)]
    if (any(grepl("^[[:space:]]*(\\{[[:space:]]*\"|\\[[[:space:]]*[\\{\"0-9-])", z))) stop("Embedded serialized transport in RDS")
    if (any(grepl("[\u3400-\u4dbf\u4e00-\u9fff]", z))) stop("Non-English generated RDS text")
  } else if (is.list(x)) invisible(lapply(x, audit_native_strings))
  invisible(TRUE)
}
assert_predictions <- function(comp, original = NULL, comparison_directory = NULL) {
  telemetry <- c("cache_g_build", "cache_g_hit", "cache_g_hysteresis")
  collect <- function(z) lapply(c("2N", "4N"), function(cohort) lapply(z[[paste0("run_", cohort)]]$segment_results, function(s) {
    stopifnot(is.finite(s$selection$selected_day), all(is.finite(s$selection$selected_frac)),
              all(is.finite(s$selection$reseeded_state)), all(is.finite(s$sim$Ntot_live_obs)),
              all(is.finite(s$sim$live_state_obs)), all(is.finite(s$sim$dead_hypoxia_state_obs)),
              all(is.finite(s$sim$dead_buffer_state_obs)))
    # Preserve every saved selection field, propagated state and full simulated distributions.
    list(selection = s$selection, state = s$sim[setdiff(names(s$sim), telemetry)])
  }))
  actual <- collect(comp)
  if (!is.null(original)) {
    expected <- collect(original)
    if (!is.null(comparison_directory)) {
      rows <- list()
      for (cohort in c("2N", "4N")) {
        a <- comp[[paste0("run_", cohort)]]$segment_results
        b <- original[[paste0("run_", cohort)]]$segment_results
        stopifnot(identical(vapply(a, function(z) z$selection$selected_day, numeric(1)), vapply(b, function(z) z$selection$selected_day, numeric(1))))
        for (i in seq_along(a)) for (n in telemetry) rows[[length(rows) + 1L]] <- data.frame(cohort = cohort, segment_index = i, metric = n, actual = a[[i]]$sim[[n]], archived = b[[i]]$sim[[n]])
      }
      o2prof_write_tsv(do.call(rbind, rows), file.path(comparison_directory, "cache_telemetry_comparison.tsv"))
    }
    comparison <- all.equal(actual, expected, tolerance = 1e-12, check.attributes = TRUE)
    if (!isTRUE(comparison)) stop("Prediction/state replay mismatch: ", paste(comparison, collapse = "; "))
  }
  invisible(TRUE)
}
score_components <- function(comp) c(objective = comp$objective, growth = mean(comp$growth_df$loglik),
  karyotype = mean(comp$ploidy_df$mean_loglik), flow = mean(comp$flow_df$mean_loglik))
replay_fit <- function(ctx, map, fit, out) {
  p <- o2prof_project_start(map, fit$best_params)$params
  o2prof_to_free(map, p)
  comp <- ctx$raw_objective(p, ctx$objects, ctx$cfg)
  o2prof_check_samples(comp, ctx$manifest); assert_predictions(comp, fit$best_components, out)
  independent <- o2prof_independent_score(ctx, comp)
  stopifnot(max(abs(score_components(comp) - score_components(fit$best_components))) <= 1e-8,
            abs(independent$objective - comp$objective) <= 1e-8)
  if (length(map$fixed)) for (n in names(map$fixed)) stopifnot(p[[n]] == map$fixed[[n]])
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_kv(c(as.list(score_components(comp)), list(independent = independent$objective, pass = TRUE)), file.path(out, "replay.tsv"))
  assert_runtime(out)
  comp
}
original_map <- function(ctx, fixed = list()) o2prof_map(ctx$cfg$parameter_table, ctx$saved$best_params, "zero_extended", fixed)

prepare_local_start <- function(ctx, map, task, project, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  input <- readRDS(file.path(project, task$source_fit))
  projected <- o2prof_project_start(map, input$best_params)
  o2prof_write_tsv(projected$projection, file.path(out, "projection.tsv"))
  if (task$phase != "neighbor") {
    comp <- replay_fit(ctx, map, input, file.path(out, "source_replay"))
    return(list(skipped = FALSE, fit = list(best_params = projected$params, best_components = comp)))
  }
  stopifnot(!is.na(task$fallback_source))
  fallback <- readRDS(file.path(project, task$fallback_source))
  destination <- o2prof_project_start(map, fallback$best_params)$params
  source_unit <- o2review_unit_spec(map, projected$params)$unit
  dest_unit <- o2review_unit_spec(map, destination)$unit
  spec <- o2prof_spec(map, destination); attempts <- list()
  for (fraction in c(1, .75, .5, .25)) {
    free <- setNames(spec$lower + (fraction * source_unit + (1 - fraction) * dest_unit) * (spec$upper - spec$lower), spec$param_name)
    p <- if (fraction == 1) projected$params else o2prof_from_free(map, free)
    error <- ""
    comp <- tryCatch({z <- ctx$raw_objective(p, ctx$objects, ctx$cfg); o2prof_check_samples(z, ctx$manifest); z},
      error = function(e) {if (!inherits(e, "invitro_protocol_infeasible")) stop(e); error <<- conditionMessage(e); NULL})
    attempts[[length(attempts) + 1L]] <- data.frame(source_fraction = fraction, feasible = !is.null(comp),
      objective = if (is.null(comp)) Inf else comp$objective, error = gsub("[\t\r\n]", " ", error))
    saveRDS(p, file.path(out, paste0("proposal_", length(attempts), ".rds")))
    if (!is.null(comp)) {
      assert_predictions(comp); independently <- o2prof_independent_score(ctx, comp)
      stopifnot(abs(comp$objective - independently$objective) <= 1e-8)
      o2prof_write_tsv(do.call(rbind, attempts), file.path(out, "transfer_attempts.tsv"))
      return(list(skipped = FALSE, fit = list(best_params = p, best_components = comp)))
    }
  }
  o2prof_write_tsv(do.call(rbind, attempts), file.path(out, "transfer_attempts.tsv"))
  list(skipped = TRUE)
}

finalize_fit <- function(ctx, map, out, fingerprint, task, native = NULL, local = NULL) {
  files <- list.files(file.path(out, "best_candidates"), pattern = "\\.rds$", full.names = TRUE)
  stopifnot(length(files) > 0)
  candidates <- lapply(files, readRDS); checks <- list(); winner <- NULL
  for (i in seq_along(candidates)) {
    z <- candidates[[i]]
    if (!identical(z$fingerprint, fingerprint)) stop("Stale checkpoint")
    o2prof_to_free(map, z$params)
    comp <- ctx$raw_objective(z$params, ctx$objects, ctx$cfg); o2prof_check_samples(comp, ctx$manifest)
    independent <- o2prof_independent_score(ctx, comp); assert_predictions(comp)
    stopifnot(abs(z$objective - comp$objective) <= 1e-8, abs(comp$objective - independent$objective) <= 1e-8)
    checks[[i]] <- data.frame(phase = z$phase, saved = z$objective, replay = comp$objective, independent = independent$objective)
    if (is.null(winner) || comp$objective < winner$objective) winner <- list(best_params = z$params, best_components = comp, cfg = ctx$cfg, objective = comp$objective)
  }
  o2pilot_atomic_rds(winner, file.path(out, "accepted_fit.rds"))
  o2prof_write_tsv(do.call(rbind, checks), file.path(out, "candidate_checks.tsv"))
  tab <- map$table; tab$estimate <- vapply(tab$parameter, function(n) winner$best_params[[n]], numeric(1))
  tab$data_inactive_at_endpoint <- (tab$parameter == "k_o_mis" & winner$best_params$p_misseg == 0) |
    (tab$parameter == "buffer_n_exp" & winner$best_params$buffer_beta == 0) |
    (tab$parameter %in% c("buffer_beta", "buffer_n_exp") & winner$best_params$buffer_smax == 0)
  tab$exact_lower <- tab$estimate == tab$lower; tab$exact_upper <- tab$estimate == tab$upper
  spec <- o2prof_spec(map, winner$best_params); u <- o2review_unit_spec(map, winner$best_params)$unit
  tab$near_lower <- NA; tab$near_upper <- NA
  tab$near_lower[tab$status == "active"] <- u <= .001; tab$near_upper[tab$status == "active"] <- u >= .999
  o2prof_write_tsv(tab, file.path(out, "parameter_result.tsv"))
  result <- c(as.list(score_components(winner$best_components)), list(task_id = task$task_id, origin = task$origin,
    family = task$family, mode = task$mode, fingerprint = fingerprint, point = task$point,
    target = task$target, target_value = task$value, skipped = FALSE, alternative = task$alternative,
    start_objective = if (is.null(task$start_objective)) NA_real_ else task$start_objective,
    local_stable = !is.null(local) && isTRUE(local$tested_directional_stability),
    budget_limited = if (!is.null(local)) isTRUE(local$budget_limited) else native$optimizer_trace$iter_completed >= task$max_generations,
    native_global_objective = if (is.null(native)) NA_real_ else native$optimizer_trace$deoptim_objective,
    actual_generations = if (is.null(native)) 0L else as.integer(native$optimizer_trace$iter_completed)))
  write_kv(result, file.path(out, "result.tsv"))
  assert_runtime(out)
  audit_improvements(ctx, map, out, fingerprint, write_index = TRUE)
  paths <- list.files(out, full.names = TRUE, recursive = TRUE)
  paths <- paths[!dir.exists(paths) & !grepl("/evaluation_logs/", paths)]
  o2prof_write_tsv(data.frame(path = substring(paths, nchar(out) + 2L), sha256 = vapply(paths, o2prof_hash, character(1))), file.path(out, "payload_manifest.tsv"))
  invisible(result)
}
run_local <- function(ctx, map, source, out, fingerprint, task, fixture = FALSE) {
  fit <- readRDS(source); p <- o2prof_project_start(map, fit$best_params)$params
  replay_fit(ctx, map, fit, file.path(out, "start_replay"))
  saveRDS(p, file.path(out, "actual_local_start.rds"))
  spec <- o2review_unit_spec(map, p); start <- as.numeric(Sys.time()); deadline <- start + if (fixture) 90 else 20700
  ev <- o2review_evaluator(ctx, map, out, fingerprint, "local", deadline)
  # Every nloptr invocation is observed inside the callback; this is the actual optimizer entry.
  initial_seen <- FALSE
  fn <- function(x) {
    if (!initial_seen) {
      stopifnot(identical(as.numeric(x), as.numeric(spec$unit)))
      initial_seen <<- TRUE; saveRDS(x, file.path(out, "native_sbplx_start.rds"))
    }
    ev$fn(x)
  }
  opts <- list(algorithm = "NLOPT_LN_SBPLX", maxeval = if (fixture) 2L else 9990L,
    maxtime = if (fixture) 10 else min(10800, deadline - as.numeric(Sys.time())),
    xtol_abs = rep(1e-4, length(spec$unit)), ftol_abs = 1e-8, print_level = 0)
  old_stage <- options(o2prof.local_stage = "sbplx"); on.exit(options(old_stage), add = TRUE)
  sb <- nloptr::nloptr(spec$unit, fn, lb = rep(0, length(spec$unit)), ub = rep(1, length(spec$unit)), opts = opts)
  saveRDS(sb, file.path(out, "sbplx.rds")); stopifnot(initial_seen, ev$counts$n <= if (fixture) 10L else 10000L)
  # Restart pattern search from the best observed feasible point, including rejected native trials.
  files <- list.files(file.path(out, "best_candidates"), full.names = TRUE, pattern = "\\.rds$")
  all <- lapply(files, readRDS); best <- all[[which.min(vapply(all, function(z) z$objective, numeric(1)))]]
  unit <- o2review_unit_spec(map, best$params)$unit
  dirs <- o2review_directions(map)
  if (isTRUE(task$alternative)) dirs <- dirs[rev(seq_len(nrow(dirs))), , drop = FALSE]
  mesh <- if (isTRUE(task$alternative)) c(.01, .005, .001, .0001) else c(.02, .005, .001, .0001)
  options(o2prof.local_stage = "pattern")
  pat <- o2review_pattern(ev$fn, unit, dirs, if (fixture) 5L else min(9500L, 19500L - ev$counts$n), deadline, mesh)
  options(o2prof.local_stage = "boundary")
  boundary <- list(); improved <- FALSE
  for (i in seq_along(pat$par)) if (pat$par[i] <= .001 || pat$par[i] >= .999) for (h in c(.001, .01, .05)) {
    if (ev$counts$n >= 20000L || as.numeric(Sys.time()) >= deadline) { pat$budget_limited <- TRUE; break }
    z <- pat$par; z[i] <- if (z[i] < .5) h else 1 - h
    q <- ev$fn(z)
    boundary[[length(boundary) + 1L]] <- data.frame(coordinate = spec$param_name[i], inward_fraction = h, endpoint = pat$value, probe = q)
    if (q < pat$value - 1e-5) improved <- TRUE
  }
  if (improved || pat$budget_limited) pat$tested_directional_stability <- FALSE
  # A best checkpoint below the audited endpoint invalidates the directional claim.
  all <- lapply(files, readRDS)
  if (min(vapply(all, function(z) z$objective, numeric(1))) < pat$value - 1e-8) pat$tested_directional_stability <- FALSE
  pat$evaluations <- ev$counts$n; pat$sbplx_status <- sb$status
  o2prof_write_tsv(pat$sweeps, file.path(out, "local_sweeps.tsv"))
  if (length(boundary)) o2prof_write_tsv(do.call(rbind, boundary), file.path(out, "boundary_probes.tsv"))
  saveRDS(pat, file.path(out, "local_search.rds"))
  write_kv(list(method = "NLOPT_LN_SBPLX_then_bounded_pattern_then_boundary", sbplx_q = best$objective,
    sbplx_calls = sb$iterations, total_calls = ev$counts$n, deadline_seconds = deadline - start), file.path(out, "method_sequence.tsv"))
  finalize_fit(ctx, map, out, fingerprint, task, local = pat)
}

# Validate append-only records without simulating their vectors again.
audit_improvements <- function(ctx, map, out, fingerprint, write_index = FALSE) {
  directory <- file.path(out, "best_candidates")
  paths <- list.files(file.path(directory, "improvements"), pattern = "\\.rds$", full.names = TRUE)
  if (!length(paths)) stop("Missing immutable improvement archive")
  expected_days <- do.call(rbind, lapply(c("2N", "4N"), function(cohort) {
    z <- ctx$saved$best_components[[paste0("run_", cohort)]]$segment_results
    data.frame(cohort = cohort, segment_order = seq_along(z), segment_id = vapply(z, function(v) v$segment$segment_id, character(1)))
  }))
  rows <- lapply(paths, function(path) {
    r <- readRDS(path)
    audit_native_strings(r)
    if (!identical(r$fingerprint, fingerprint) || !isTRUE(r$complete_components)) stop("Stale or incomplete improvement archive")
    o2prof_to_free(map, r$params)
    comp <- c(list(objective = r$objective, n_growth_missing_pred = 0L), r$score_tables)
    o2prof_check_samples(comp, ctx$manifest)
    means <- c(growth = mean(comp$growth_df$loglik), karyotype = mean(comp$ploidy_df$mean_loglik), flow = mean(comp$flow_df$mean_loglik))
    stopifnot(identical(r$component_means, means), abs(-sum(means) - r$objective) <= 1e-8,
      identical(r$selected_days[, names(expected_days)], expected_days), all(is.finite(r$selected_days$selected_day)))
    data.frame(pid = r$pid, record_index = r$record_index, evaluation_index = r$evaluation_index,
      phase = r$phase, objective = r$objective, growth = means[1], karyotype = means[2], flow = means[3],
      path = substring(path, nchar(out) + 2L), sha256 = o2prof_hash(path), stringsAsFactors = FALSE)
  })
  index <- do.call(rbind, rows)
  for (pid in unique(index$pid)) {
    z <- index[index$pid == pid, ]; z <- z[order(z$record_index), ]
    stopifnot(identical(as.integer(z$record_index), seq_len(nrow(z))), all(diff(z$objective) < 0))
    pointer <- readRDS(file.path(directory, paste0("pid_", pid, ".rds")))
    stopifnot(pointer$record_index == nrow(z), pointer$objective == tail(z$objective, 1),
      pointer$archive_sha256 == o2prof_hash(file.path(directory, pointer$archive_relative_path)))
    archived <- readRDS(file.path(directory, pointer$archive_relative_path))
    stopifnot(identical(pointer$params, archived$params), identical(pointer$fingerprint, fingerprint))
  }
  if (write_index) o2prof_write_tsv(index[order(index$pid, index$record_index), ], file.path(out, "improvements.tsv"))
  invisible(index)
}
