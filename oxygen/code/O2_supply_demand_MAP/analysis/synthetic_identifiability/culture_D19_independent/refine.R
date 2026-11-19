# Per-origin refinements. Every model call passes through the charged evaluator.
blind_refine <- function(context, evaluator, out, stage) {
  state <- evaluator$state
  if (stage == 'restore' && !is.null(state$gate)) return(list(refinement_status = 'ALREADY_STRICT_PASS'))
  grid <- joint_mask_grid(context)
  left <- function() d_limits[[stage]] - state$counts[[stage]]
  probes <- list()
  evaluate <- function(unit, kind, parent_branch = NA_character_, coordinate = NA_integer_,
                       step_size = NA_real_, bucket = stage) {
    evaluator$fn(unit, bucket = bucket)
    z <- state$last
    if (!z$feasible) {
      probes[[length(probes) + 1L]] <<- data.frame(call = state$n, kind = kind, coordinate = coordinate,
        signed_step = step_size, parent_branch = parent_branch, branch = NA_character_, feasible = FALSE, same_branch = FALSE)
      return(list(ok = FALSE))
    }
    p <- d_decode(context, unit)
    pruning <- joint_pruning(context, p, grid)
    residual <- joint_residual(context, z$pred, p)
    stopifnot(abs(sum(residual^2) / 2 - z$loss$value) <= 1e-12)
    point <- list(ok = TRUE, candidate = d_candidate(context, unit, z, state$best$origin, stage),
                  branch = joint_branch(z, pruning, p), cutoff_margin = pruning$cutoff_margin,
                  metrics = exact_metrics(context, z), residual = residual, kary = z$pred$kary)
    probes[[length(probes) + 1L]] <<- data.frame(call = state$n, kind = kind, coordinate = coordinate,
      signed_step = step_size, parent_branch = parent_branch, branch = point$branch,
      feasible = TRUE, same_branch = !is.na(parent_branch) && identical(point$branch, parent_branch))
    point$kary <- lapply(point$kary, function(q) q / sum(q))
    point
  }
  best <- state$best
  point <- evaluate(best$unit, 'exact_handoff', bucket = 'handoff')
  stopifnot(point$ok, identical(best$key, point$candidate$key),
            abs(best$score - point$candidate$score) <= 1e-8)
  radius <- if (stage == 'joint') .03 else .006
  damping <- 1e-4
  priority <- match(fit_priority, context$free$parameter)
  stopifnot(!anyNA(priority))
  status <- 'EVALUATION_LIMIT'
  iteration <- 0L
  derivatives <- proposals <- list()
  while (left() >= 40L) {
    if (stage == 'restore' && !is.null(state$gate)) {
      status <- 'STRICT_PASS_FOUND'
      break
    }
    iteration <- iteration + 1L
    derivative <- joint_jacobian(point, evaluate, function() max(0L, left() - 2L))
    derivatives[[iteration]] <- data.frame(iteration = iteration, anchor_key = point$candidate$key,
      parameter = context$free$parameter, scheme = derivative$scheme, step = derivative$step,
      norm = sqrt(colSums(derivative$matrix^2)))
    if (!derivative$complete || all(derivative$scheme == 'blocked')) {
      status <- 'BRANCH_OR_EVALUATION_LIMIT'
      break
    }
    accepted <- FALSE
    for (attempt in seq_len(7L)) {
      if (left() <= 0L || radius < 2e-6) break
      choices <- if (stage == 'restore' && iteration == 1L)
        list(root_and_width = priority, full = seq_len(19L)) else list(full = seq_len(19L))
      winners <- list()
      for (name in names(choices)) {
        if (left() <= 0L) break
        if (stage == 'joint') {
          proposal <- joint_step(point$candidate$unit, point$residual, derivative$matrix, radius, damping)
          unit <- proposal$proposal
          predicted <- proposal$predicted_gain
          solver <- proposal$solver_status
        } else {
          proposal <- bounded_proposal(point, derivative$matrix, context$case$data$prob, choices[[name]], radius)
          unit <- proposal$unit
          predicted <- proposal$original_linear_violation - proposal$linear_violation
          solver <- proposal$status
        }
        if (proposal$step <= 1e-10) next
        next_point <- evaluate(unit, 'proposal', parent_branch = point$branch)
        switched <- next_point$ok && !identical(next_point$branch, point$branch)
        gain <- if (next_point$ok) point$candidate$score - next_point$candidate$score else -Inf
        proposals[[length(proposals) + 1L]] <- data.frame(iteration = iteration, attempt = attempt,
          call = state$n, proposal = name, radius = radius, step = proposal$step,
          predicted_improvement = predicted, actual_D_gain = gain, solver_status = solver,
          branch_changed = switched, feasible = next_point$ok,
          strict_pass = next_point$ok && next_point$metrics$pass,
          true_violation = if (next_point$ok) next_point$metrics$violation else Inf)
        if (if (stage == 'joint') gain > 0 else is_improved(next_point, point))
          winners[[length(winners) + 1L]] <- list(value = next_point, gain = gain, predicted = predicted, switched = switched)
      }
      if (length(winners)) {
        index <- 1L
        if (length(winners) > 1L && is_improved(winners[[2L]]$value, winners[[1L]]$value)) index <- 2L
        winner <- winners[[index]]
        if (stage == 'joint') {
          ratio <- winner$gain / max(winner$predicted, 1e-16)
          if (ratio > .75 && !winner$switched) radius <- min(.15, radius * 1.7)
          if (ratio < .25 || winner$switched) radius <- max(1e-6, radius * .6)
          damping <- if (ratio > .5) max(1e-8, damping / 3) else min(10, damping * 3)
        } else {
          ratio <- (point$metrics$violation - winner$value$metrics$violation) / max(point$metrics$violation, 1e-16)
          if (winner$switched) radius <- max(2e-6, radius * .75)
          else if (ratio > .15) radius <- min(.03, radius * 1.5)
        }
        point <- winner$value
        accepted <- TRUE
        break
      }
      radius <- radius * .4
      damping <- min(10, damping * 4)
    }
    if (!accepted) {
      status <- 'NO_EXACT_IMPROVING_STEP'
      break
    }
  }
  if (stage == 'restore' && !is.null(state$gate)) status <- 'STRICT_PASS_FOUND'
  sf_emit(do.call(rbind, probes), file.path(out, 'branch_probes.tsv'))
  if (length(derivatives)) sf_emit(do.call(rbind, derivatives), file.path(out, 'derivatives.tsv'))
  if (length(proposals)) sf_emit(do.call(rbind, proposals), file.path(out, 'proposals.tsv'))
  list(refinement_status = status)
}
