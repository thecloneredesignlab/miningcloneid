# Residual and bounded local-update helpers for the C/E method pilot.
joint_residual <- function(context, prediction, params) {
  stopifnot(identical(params$sigma_growth, context$case$known_sigma_growth))
  growth <- (prediction$growth - context$case$data$growth) / params$sigma_growth
  target <- context$case$data$prob
  karyotype <- lapply(names(target), function(id) {
    probability <- prediction$kary[[id]]
    sqrt(probability / sum(probability)) - sqrt(target[[id]])
  })
  residual <- c(unname(growth) / sqrt(length(growth)),
                unlist(karyotype, use.names = FALSE) / sqrt(length(target)))
  stopifnot(length(residual) == 1710L, all(is.finite(residual)))
  residual
}

joint_mask_grid <- function(context) {
  grid <- expand.grid(N = context$core$grid_pre, oxygen = sort(unique(context$contract$edges$oxygen_pct)))
  multiplicities <- lapply(grid$N, function(n) seq.int(0L, n))
  list(N = grid$N, oxygen = grid$oxygen, multiplicity = unlist(multiplicities),
       event_row = rep(seq_along(grid$N), lengths(multiplicities)))
}

joint_pruning <- function(context, params, grid) {
  probability <- context$env$.pmisseg_of_O2(grid$oxygen, params, grid$N)
  event <- dbinom(grid$multiplicity, grid$N[grid$event_row], probability[grid$event_row])
  stopifnot(length(probability) == length(grid$N), all(is.finite(event)))
  list(pruned = event < 1e-8, rate = probability,
       cutoff_margin = min(abs(event - 1e-8)) / 1e-8)
}

joint_branch <- function(propagated, pruning, params) {
  decisions <- lapply(propagated$lineages, function(a) {
    vapply(a$edges, function(edge) edge$selection$reseed_mode, character(1))
  })
  sf_hash(list(selection = propagated$pred$selection, reset_mode = decisions,
               pruned = pruning$pruned, exact_survival_zero = params$buffer_smax == 0,
               karyotype_zero = lapply(propagated$pred$kary, function(q) q == 0)))
}

joint_pass <- function(loss) {
  loss$value <= 1e-4 && loss$rms_growth <= .01 &&
    loss$max_growth <= .03 && loss$max_kary_tv <= .005
}

joint_jacobian <- function(anchor, evaluate, left, step = 2e-5) {
  unit <- anchor$candidate$unit
  jacobian <- matrix(0, length(anchor$residual), length(unit))
  scheme <- rep('blocked', length(unit))
  step_size <- numeric(length(unit))
  if (anchor$cutoff_margin < 1e-8) {
    return(list(matrix = jacobian, scheme = scheme, step = step_size, complete = FALSE))
  }
  for (i in seq_along(unit)) {
    for (h in step * c(1, .25, .0625, .015625)) {
      side <- c(-1, 1)[unit[i] + h * c(-1, 1) >= 0 & unit[i] + h * c(-1, 1) <= 1]
      if (left() < length(side)) {
        return(list(matrix = jacobian, scheme = scheme, step = step_size, complete = FALSE))
      }
      probes <- lapply(side, function(sign) {
        proposal <- unit
        proposal[i] <- proposal[i] + sign * h
        evaluate(proposal, 'jacobian', parent_branch = anchor$branch, coordinate = i,
                 step_size = sign * h)
      })
      usable <- vapply(probes, function(p) {
        p$ok && identical(p$branch, anchor$branch) && p$cutoff_margin >= 1e-8
      }, logical(1))
      side <- side[usable]
      probes <- probes[usable]
      if (!length(side)) next
      if (length(side) == 2L) {
        jacobian[, i] <- (probes[[2L]]$residual - probes[[1L]]$residual) / (2 * h)
        scheme[i] <- 'central'
      } else {
        jacobian[, i] <- (probes[[1L]]$residual - anchor$residual) / (side * h)
        scheme[i] <- if (side > 0) 'forward' else 'backward'
      }
      step_size[i] <- h
      break
    }
  }
  stopifnot(all(is.finite(jacobian)))
  list(matrix = jacobian, scheme = scheme, step = step_size, complete = TRUE)
}

joint_step <- function(unit, residual, jacobian, radius, damping) {
  gradient <- as.numeric(crossprod(jacobian, residual))
  curvature <- crossprod(jacobian)
  column_scale <- pmax(diag(curvature), max(diag(curvature), 1e-12) * 1e-8)
  regularized <- curvature + diag(damping * column_scale, nrow(curvature))
  lower <- pmax(-radius, -unit)
  upper <- pmin(radius, 1 - unit)
  active <- which(diag(curvature) > 0 & lower < upper)
  delta <- numeric(length(unit))
  solver_status <- NA_integer_
  if (length(active)) {
    H <- regularized[active, active, drop = FALSE]
    g <- gradient[active]
    # Scale the box problem by column norms and the positive residual baseline.
    # An explicit projected-gradient candidate protects against an inner solver returning zero.
    scale <- sqrt(column_scale[active])
    base <- max(sum(residual^2) / 2, 1e-16)
    model <- function(s) sum(g * s) + .5 * sum(s * (H %*% s))
    fit <- optim(rep(0, length(active)), method = 'L-BFGS-B',
                 fn = function(s) 1 + model(s / scale) / base,
                 gr = function(s) (g + as.numeric(H %*% (s / scale))) / scale / base,
                 lower = lower[active] * scale, upper = upper[active] * scale,
                 control = list(maxit = 140L, factr = 10, pgtol = 1e-10))
    stopifnot(all(is.finite(fit$par)))
    clamp <- function(s) pmin(upper[active], pmax(lower[active], s))
    direction <- clamp(-g / diag(H))
    direction_curvature <- sum(direction * (H %*% direction))
    multiplier <- if (direction_curvature > 0) max(0, min(1, -sum(g * direction) / direction_curvature)) else 0
    candidates <- list(clamp(fit$par / scale), multiplier * direction)
    values <- vapply(candidates, model, numeric(1))
    delta[active] <- candidates[[which.min(values)]]
    solver_status <- fit$convergence
  }
  proposal <- pmin(1, pmax(0, unit + delta))
  modeled <- residual + as.numeric(jacobian %*% (proposal - unit))
  list(proposal = proposal, predicted_gain = (sum(residual^2) - sum(modeled^2)) / 2,
       step = max(abs(proposal - unit)), solver_status = solver_status,
       projected_gradient = max(abs(pmin(upper, pmax(lower, -gradient)))))
}
