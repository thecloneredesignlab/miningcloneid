args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
archive <- normalizePath(args[1L])
pilot <- normalizePath(args[2L])
task <- args[3L]
pilot_release <- normalizePath(args[4L], mustWork = TRUE)
out <- file.path(pilot, 'verification', task)
release <- file.path(archive, 'release')
for (file in c('production_bindings.R', 'interface.R', 'feasibility.R', 'search_helpers.R', 'core.R')) {
  source(file.path(release, file))
}
source(file.path(pilot_release, 'numerics.R'))

main <- function() {
  on.exit(sf_inventory(out), add = TRUE)
  env <- sf_backend(archive, out, native = TRUE)
  launch <- sf_tab(file.path(pilot_release, 'origins.tsv'))
  stopifnot(identical(launch$origin, c('D_C', 'D_E')))
  cap <- c(C2 = 800L, E2 = 798L, P2 = 16L, V3 = 4L, V4 = 4L)[[task]]
  stopifnot(cap > 0, !file.exists(file.path(out, 'attempts.tsv')))
  ledger <- file(file.path(out, 'attempts.tsv'), 'wt')
  on.exit(close(ledger), add = TRUE)
  fields <- c('call', 'origin', 'kind', 'event', 'key', 'feasible', 'D', 'seconds', 'branch',
              'same_branch', 'coordinate', 'signed_step', 'cutoff_margin', paste0('u', 1:19))
  writeLines(paste(fields, collapse = '\t'), ledger)
  fmt <- function(x) sprintf('%.17g', x)
  call <- 0L
  times <- numeric()
  remaining <- function() cap - call
  checks <- list()
  checked <- function(name, condition, detail = '') {
    checks[[length(checks) + 1L]] <<- data.frame(check = name, status = if (isTRUE(condition)) 'PASS' else 'FAIL', detail = detail)
    sf_emit(do.call(rbind, checks), file.path(out, 'checks.tsv'))
    if (!isTRUE(condition)) stop('Failed check: ', name)
  }

  # Keep the model closures in their defining environment and reset its wrappers between origins.
  base_functions <- as.list(env, all.names = TRUE)
  context_for <- function(origin) {
    env$ivt_select_segment_observation <- base_functions$ivt_select_segment_observation
    env$ivt_extract_passage_end_state <- base_functions$ivt_extract_passage_end_state
    row <- launch[launch$origin == origin, , drop = FALSE]
    original <- readRDS(file.path(archive, row$state_file))$candidate
    native <- d_read(archive, origin)
    context <- d_context(env, native$case, native$contract, native$id)
    stopifnot(identical(original$key, row$key), identical(original$origin, origin),
              identical(original$key, sf_hash(original$unit)), identical(original$identity, context$identity),
              identical(original$params, d_decode(context, original$unit)[cg_roles]))
    grid <- joint_mask_grid(context)
    evaluate <- function(unit, kind, parent_branch = NA_character_, coordinate = NA_integer_,
                         step_size = NA_real_, independent = FALSE, keep_states = FALSE) {
      stopifnot(remaining() > 0L)
      unit <- as.numeric(unit)
      stopifnot(length(unit) == 19L, all(is.finite(unit)), all(unit >= 0 & unit <= 1))
      call <<- call + 1L
      key <- sf_hash(unit)
      emit <- function(event, ok = NA, score = NA_real_, elapsed = NA_real_, branch = NA_character_,
                       same_branch = NA, cutoff_margin = NA_real_) {
        writeLines(paste(c(call, origin, kind, event, key, ok, fmt(score), fmt(elapsed), branch,
                           same_branch, coordinate, fmt(step_size), fmt(cutoff_margin), fmt(unit)),
                         collapse = '\t'), ledger)
        flush(ledger)
      }
      emit('START')
      begin <- proc.time()[3L]
      params <- d_decode(context, unit)
      native <- d_physics(context, params, reference = independent)
      if (!native$feasible) {
        elapsed <- proc.time()[3L] - begin
        times <<- c(times, elapsed)
        emit('END', FALSE, Inf, elapsed)
        failure <- if (independent) native$error else lapply(native$lineages, function(x) x$error)
        sf_rds(failure, file.path(out, 'failures', sprintf('%04d.rds', call)))
        return(list(ok = FALSE))
      }
      residual <- joint_residual(context, native$pred, params)
      stopifnot(abs(sum(residual^2) / 2 - native$loss$value) <= 1e-12 * max(1, native$loss$value))
      pruning <- joint_pruning(context, params, grid)
      branch <- if (independent) NA_character_ else joint_branch(native, pruning, params)
      candidate <- d_candidate(context, unit, native, origin, task)
      elapsed <- proc.time()[3L] - begin
      times <<- c(times, elapsed)
      emit('END', TRUE, candidate$score, elapsed, branch,
           if (is.na(parent_branch) || is.na(branch)) NA else identical(parent_branch, branch),
           pruning$cutoff_margin)
      list(ok = TRUE, candidate = candidate, residual = residual, branch = branch,
           cutoff_margin = pruning$cutoff_margin, prediction = if (keep_states) native$pred else NULL)
    }
    list(context = context, original = original, grid = grid, evaluate = evaluate)
  }

  confirm_seed <- function(model, point) {
    initial <- model$original
    checked(paste0(initial$origin, '_archived_seed'), point$ok &&
      identical(point$candidate$key, initial$key) && identical(point$candidate$params, initial$params) &&
      abs(point$candidate$score - initial$score) <= 1e-8)
  }

  if (task %in% c('C2', 'E2')) {
    origin <- paste0('D_', substr(task, 1L, 1L))
    model <- context_for(origin)
    initial <- model$original
    point <- model$evaluate(initial$unit, 'archived_seed')
    confirm_seed(model, point)
    previous <- readRDS(file.path(pilot, 'verification', substr(task, 1L, 1L), 'checkpoint.rds'))
    stopifnot(previous$status == 'NO_ACCEPTABLE_JOINT_STEP', previous$calls + cap == 850L,
              identical(initial$key, previous$best$key))
    start_score <- point$candidate$score
    gate <- NULL
    register <- function(item) {
      if (item$ok && joint_pass(item$candidate$loss) && (is.null(gate) || item$candidate$score < gate$score)) {
        gate <<- item$candidate
      }
    }
    checkpoint <- function(status) sf_rds(list(origin = origin, start = initial, best = point$candidate,
                                                gate = gate, calls = call, status = status), file.path(out, 'checkpoint.rds'))
    # The first pass already exhausted the prespecified width grid at this exact point.
    radius <- .03
    damping <- 1e-4
    iterations <- list()
    derivatives <- list()
    status <- 'EVALUATION_ALLOWANCE'
    iteration <- 0L
    checkpoint('RUNNING')
    while (remaining() > 39L) {
      iteration <- iteration + 1L
      anchor <- point
      measure <- function(...) {
        p <- model$evaluate(...)
        register(p)
        p
      }
      derivative <- joint_jacobian(anchor, measure, remaining)
      derivatives[[iteration]] <- data.frame(iteration = iteration, key = anchor$candidate$key,
        parameter = model$context$free$parameter, scheme = derivative$scheme,
        unit_step = derivative$step, column_norm = sqrt(colSums(derivative$matrix^2)))
      if (!derivative$complete || all(derivative$scheme == 'blocked')) {
        status <- 'BRANCH_OR_ALLOWANCE_LIMIT'
        break
      }
      if (iteration == 1L) {
        previous_helpers <- new.env(parent = .GlobalEnv)
        sys.source(file.path(pilot, 'release', 'numerics.R'), previous_helpers)
        original_step <- previous_helpers$joint_step(anchor$candidate$unit, anchor$residual,
                                                     derivative$matrix, radius, damping)
        revised_step <- joint_step(anchor$candidate$unit, anchor$residual, derivative$matrix, radius, damping)
        sf_emit(data.frame(origin = origin, original_step = original_step$step,
                           original_predicted_gain = original_step$predicted_gain,
                           revised_step = revised_step$step, revised_predicted_gain = revised_step$predicted_gain,
                           projected_gradient = revised_step$projected_gradient,
                           revised_solver_status = revised_step$solver_status), file.path(out, 'inner_solver_comparison.tsv'))
        checked(paste0(origin, '_joint_step_regression'), revised_step$step > 1e-9 && revised_step$predicted_gain > 1e-14)
      }
      accepted <- FALSE
      for (attempt in seq_len(8L)) {
        if (remaining() <= 0L) break
        step <- joint_step(anchor$candidate$unit, anchor$residual, derivative$matrix, radius, damping)
        if (step$step < 1e-9 || step$predicted_gain <= 1e-14) break
        candidate <- measure(step$proposal, 'joint_trial', parent_branch = anchor$branch, step_size = step$step)
        gain <- if (candidate$ok) anchor$candidate$score - candidate$candidate$score else -Inf
        ratio <- gain / step$predicted_gain
        accepted <- gain > 1e-14
        switched <- candidate$ok && !identical(candidate$branch, anchor$branch)
        iterations[[length(iterations) + 1L]] <- data.frame(iteration = iteration, trial = attempt, call = call,
          anchor = anchor$candidate$key, candidate = if (candidate$ok) candidate$candidate$key else NA_character_,
          radius = radius, damping = damping, predicted_gain = step$predicted_gain, actual_gain = gain,
          branch_changed = switched, ratio = ratio, accepted = accepted,
          solver_status = step$solver_status, projected_gradient = step$projected_gradient)
        if (accepted) {
          point <- candidate
          if (ratio > .75 && !switched) radius <- min(.15, radius * 1.7)
          if (ratio < .25 || switched) radius <- max(1e-6, radius * .6)
          damping <- if (ratio > .5) max(1e-8, damping / 3) else min(10, damping * 3)
          checkpoint('RUNNING')
          break
        }
        radius <- max(1e-6, radius * .4)
        damping <- min(10, damping * 4)
      }
      if (!accepted) {
        status <- 'NO_ACCEPTABLE_JOINT_STEP'
        break
      }
    }
    checked(paste0(origin, '_bounded_monotone'), point$candidate$score <= start_score && call <= cap)
    checkpoint(status)
    if (length(iterations)) sf_emit(do.call(rbind, iterations), file.path(out, 'joint_trials.tsv'))
    if (length(derivatives)) sf_emit(do.call(rbind, derivatives), file.path(out, 'jacobian_columns.tsv'))
  }

  if (task == 'P2') {
    to_replay <- list()
    rows <- list()
    for (origin in launch$origin) {
      model <- context_for(origin)
      status <- readRDS(file.path(pilot, 'verification', paste0(substring(origin, 3L), '2'), 'checkpoint.rds'))
      stopifnot(status$origin == origin, status$best$score <= model$original$score)
      anchor <- model$evaluate(status$best$unit, 'local_anchor')
      checked(paste0(origin, '_poll_anchor'), anchor$ok && identical(anchor$candidate$key, status$best$key) &&
              abs(anchor$candidate$score - status$best$score) <= 1e-8)
      verified <- anchor$candidate
      gate <- status$gate
      dirs <- r_directions(model$context)
      # Seven targeted directions per continuation fit inside the unspent verification allowance.
      coordinate <- match(c('lam_max', 'buffer_n_exp', 'sigma_kary'), model$context$free$parameter)
      direction_id <- c(rep(coordinate, each = 2L), 25L)
      signs <- c(rep(c(-1, 1), 3L), 1)
      seen <- character()
      h <- 2e-5
      for (t in seq_along(direction_id)) {
          j <- direction_id[t]
          sign <- signs[t]
          u <- r_admissible(anchor$candidate$unit, sign * dirs[j, ], h)
          if (is.null(u)) next
          key <- sf_hash(u)
          if (key %in% seen) next
          seen <- c(seen, key)
          result <- model$evaluate(u, 'registered_local_poll', parent_branch = anchor$branch,
                                   coordinate = j, step_size = sign * h)
          gain <- if (result$ok) anchor$candidate$score - result$candidate$score else -Inf
          rows[[length(rows) + 1L]] <- data.frame(origin = origin, anchor_key = anchor$candidate$key,
            key = key, mesh = h, direction = j, sign = sign, gain = gain,
            same_branch = result$ok && identical(result$branch, anchor$branch))
          if (!result$ok) next
          candidate <- result$candidate
          if (candidate$score < verified$score) verified <- candidate
          if (joint_pass(candidate$loss) && (is.null(gate) || candidate$score < gate$score)) gate <- candidate
      }
      to_replay[[origin]] <- list(start = model$original, anchor = anchor$candidate, lowest = verified,
                                  gate = gate, poll_calls = length(seen),
                                  poll_max_gain = anchor$candidate$score - verified$score)
    }
    sf_emit(do.call(rbind, rows), file.path(out, 'local_polls.tsv'))
    sf_rds(to_replay, file.path(out, 'replay_selection.rds'))
  }

  if (task %in% c('V3', 'V4')) {
    selected <- readRDS(file.path(pilot, 'verification', 'P2', 'replay_selection.rds'))
    results <- list()
    for (origin in names(selected)) {
      model <- context_for(origin)
      choice <- selected[[origin]]
      candidates <- c(choice['lowest'], if (is.null(choice$gate)) choice['anchor'] else choice['gate'])
      keys <- vapply(candidates, function(candidate) candidate$key, character(1))
      for (key in unique(keys)) {
        candidate <- candidates[[which(keys == key)[1L]]]
        saved <- model$evaluate(candidate$unit, 'fresh_replay', independent = task == 'V4', keep_states = TRUE)
        id <- paste(origin, key, sep = '_')
        checked(paste0(id, '_replay'), saved$ok && identical(saved$candidate$key, key) &&
                identical(saved$candidate$params, candidate$params) &&
                abs(saved$candidate$score - candidate$score) <= 1e-8)
        checked(paste0(id, '_independent_score'), max(abs(d_independent_loss(model$context, saved$prediction) -
                unlist(saved$candidate$loss[c('growth', 'kary', 'value')], use.names = FALSE))) <= 1e-8)
        point_file <- paste0(id, '.rds')
        if (task == 'V4') {
          first <- readRDS(file.path(pilot, 'verification', 'V3', point_file))
          checked(paste0(id, '_fresh_process'), !identical(first$task, task))
          f_compare(first$prediction, saved$prediction, checked, id)
        }
        sf_rds(list(task = task, prediction = saved$prediction, candidate = saved$candidate), file.path(out, point_file))
        loss <- saved$candidate$loss
        kary <- model$context$case$data$prob
        a19 <- grep('_2N_O[12]_A19_', names(kary), value = TRUE)
        stopifnot(length(a19) == 2L)
        a19_tv <- vapply(a19, function(n) sum(abs(saved$prediction$kary[[n]] - kary[[n]])) / 2, numeric(1))
        results[[length(results) + 1L]] <- data.frame(origin = origin, key = key,
          role = paste(names(candidates)[keys == key], collapse = '|'), start_D = choice$start$score,
          D = loss$value, D_G = loss$growth, D_K = loss$kary,
          growth_RMS_over_SD = loss$rms_growth, max_growth_over_SD = loss$max_growth,
          max_karyotype_TV = loss$max_kary_tv, A19_O1_TV = a19_tv[1L], A19_O2_TV = a19_tv[2L],
          sigma_kary = saved$candidate$params$sigma_kary, prediction_pass = joint_pass(loss),
          local_poll_anchor = identical(key, choice$anchor$key), polls_at_anchor = choice$poll_calls,
          maximum_gain_at_anchor = choice$poll_max_gain)
      }
    }
    sf_emit(do.call(rbind, results), file.path(out, 'quality.tsv'))
  }

  checked('charged_calls_within_allotment', call <= cap && length(times) == call)
  sf_emit(data.frame(task = task, status = 'PASS', forward_attempts = call, cap = cap,
                     N20 = 0L, profiles = 0L, tumors = 0L), file.path(out, 'application_result.tsv'))
}

main()
