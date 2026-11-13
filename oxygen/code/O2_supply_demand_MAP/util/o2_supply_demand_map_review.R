# Optimization review controls. The production model and objective are unchanged.
o2review_map <- function(ctx, domain = "original", fixed = list()) {
  tab <- read.csv(ctx$cfg$parameter_table, stringsAsFactors = FALSE)
  if (!domain %in% c("original", "D_base", "D_sat")) stop("Unknown review domain")
  if (domain == "D_base") tab$upper_bound[tab$param_symbol == "p_mis_base"] <- .02
  if (domain == "D_sat") tab$upper_bound[tab$param_symbol == "k_o_mis"] <- .2
  o2prof_map(tab, ctx$saved$best_params, "zero_extended", fixed = fixed)
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
  pairs <- list(c("p_misseg", "k_o_mis"), c("p_misseg", "buffer_beta"),
                c("buffer_beta", "buffer_n_exp"), c("buffer_beta", "buffer_smax"),
                c("lam_max", "mu_hp"), c("mu_hp", "gamma_mu"),
                c("O2_crit", "n_O"), c("init_mean_2N", "init_sd_2N"),
                c("init_mean_4N", "init_sd_4N"))
  for (pair in pairs) if (all(pair %in% tab$parameter)) for (sign in c(-1, 1)) {
    v <- numeric(d); v[match(pair, tab$parameter)] <- c(1, sign)
    directions <- rbind(directions, v)
    labels <- c(labels, paste(pair, collapse = if (sign == 1) "+" else "-"))
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
    started <- proc.time()[[3]]; reason <- ""; q <- Inf; days <- ""
    tryCatch({
      free <- setNames(spec$lower + pmin(pmax(unit, 0), 1) * (spec$upper - spec$lower), spec$param_name)
      params <- o2prof_from_free(map, free)
      params <- ctx$backend$normalize_run_params_common(params, cfg = ctx$cfg)
      o2prof_to_free(map, params)
      comp <- ctx$raw_objective(params, ctx$objects, ctx$cfg)
      o2prof_check_samples(comp, ctx$manifest); q <- comp$objective
      days <- paste(unlist(lapply(c("2N", "4N"), function(cohort)
        vapply(comp[[paste0("run_", cohort)]]$segment_results, function(z) z$selection$selected_day, numeric(1)))), collapse = ",")
      o2pilot_checkpoint(params, q, file.path(out, "best_candidates"), fingerprint, phase)
    }, error = function(e) {
      reason <<- conditionMessage(e); counts$failed <- counts$failed + 1L
      if (grepl("checkpoint_write_failed|Stale candidate|Immutable|Wrong free-vector", reason)) stop(e)
    })
    row <- data.frame(evaluation = counts$n, phase = phase, objective = q,
        seconds = proc.time()[[3]] - started, error = gsub("[\t\r\n]", " ", reason), selected_days = days)
    path <- file.path(out, "local_evaluations.tsv")
    write.table(row, path, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = !file.exists(path), append = file.exists(path))
    q
  }
  list(fn = fn, counts = counts, spec = spec)
}

o2review_finalize <- function(ctx, map, out, fingerprint, task, native = NULL, local = NULL) {
  files <- list.files(file.path(out, "best_candidates"), pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No complete feasible checkpoint")
  candidates <- lapply(files, readRDS)
  candidates <- candidates[order(vapply(candidates, function(x) x$objective, numeric(1)))]
  winner <- NULL; checks <- list()
  for (i in seq_along(candidates)) {
    c <- candidates[[i]]
    if (!identical(c$fingerprint, fingerprint)) stop("Stale checkpoint fingerprint")
    o2prof_to_free(map, c$params)
    comp <- ctx$raw_objective(c$params, ctx$objects, ctx$cfg)
    o2prof_check_samples(comp, ctx$manifest)
    independent <- o2prof_independent_score(ctx, comp)
    if (abs(comp$objective - c$objective) > 1e-8 || abs(independent$objective - comp$objective) > 1e-8)
      stop("Independent score replay mismatch")
    checks[[i]] <- data.frame(phase = c$phase, checkpoint = c$objective,
                              replay = comp$objective, independent = independent$objective)
    if (is.null(winner) || comp$objective < winner$objective) {
      winner <- c; winner$objective <- comp$objective; winner$components <- comp
    }
  }
  o2prof_write_tsv(do.call(rbind, checks), file.path(out, "candidate_checks.tsv"))
  saveRDS(list(best_params = winner$params, best_components = winner$components,
               cfg = ctx$cfg, origin_phase = winner$phase), file.path(out, "accepted_fit.rds"))
  tab <- map$table
  tab$estimate <- vapply(tab$parameter, function(n) winner$params[[n]], numeric(1))
  tab$data_inactive_at_endpoint <- (tab$parameter == "k_o_mis" & winner$params$p_misseg == 0) |
    (tab$parameter == "buffer_n_exp" & winner$params$buffer_beta == 0) |
    (tab$parameter %in% c("buffer_beta", "buffer_n_exp") & winner$params$buffer_smax == 0)
  tab$at_lower <- abs(tab$estimate - tab$lower) <= 1e-8 * pmax(1, abs(tab$lower))
  tab$at_upper <- abs(tab$estimate - tab$upper) <= 1e-8 * pmax(1, abs(tab$upper))
  o2prof_write_tsv(tab, file.path(out, "parameter_result.tsv"))
  numeric_params <- winner$params[vapply(winner$params, function(x) is.numeric(x) && length(x) == 1 && is.finite(x), logical(1))]
  payloads <- c("accepted_fit.rds", "parameter_map.rds", "parameter_result.tsv", "candidate_checks.tsv", "run_inputs.json")
  additional <- c("fit_result.rds", "native_de_initial_population.rds", "requested_initial_population.rds",
                  "execution_controls.json", "local_intercept_calls.rds", "local_search.rds", "local_sweeps.tsv", "boundary_probes.tsv")
  payloads <- c(payloads, additional[file.exists(file.path(out, additional))])
  result <- list(task_id = task$task_id, fingerprint = fingerprint, objective = winner$objective,
      numeric_params_exact = lapply(numeric_params, function(x) sprintf("%.17g", x)),
      component_means = list(growth = mean(winner$components$growth_df$loglik),
          karyotype = mean(winner$components$ploidy_df$mean_loglik), flow = mean(winner$components$flow_df$mean_loglik)),
      domain = task$domain, point_id = task$point_id, fixed = task$fixed, mode = task$mode,
      start_class = task$start_class, pure_independent = isTRUE(task$pure_independent),
      source_fit = task$source_fit, local_diagnostic = local,
      local_improvement = if (!is.null(local)) local$sweeps$before[1] - winner$objective else NULL,
      tested_directional_stability = !is.null(local) && isTRUE(local$tested_directional_stability),
      budget_limited = if (!is.null(local)) local$budget_limited else
        (!is.null(native) && native$optimizer_trace$iter_completed >= task$max_generations),
      native_optimizer_trace = if (!is.null(native)) native$optimizer_trace else NULL,
      actual_generations = if (!is.null(native)) native$optimizer_trace$iter_completed else 0L,
      min_generations = task$min_generations,
      optimizer_convergence_certified = FALSE, samples_valid = TRUE,
      files_sha256 = as.list(setNames(vapply(file.path(out, payloads), o2prof_hash, character(1)), payloads)),
      finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  o2pilot_json(result, file.path(out, "candidate_result.json"))
  invisible(result)
}
