# The original D19 metric remains the reported objective.
growth_count <- 114L
kary_count <- 12L
bins <- 133L
limits <- c(D = 1e-4, growth_RMS = .01, growth_max = .03, kary_TV = .005)
fit_priority <- c('init_mean_4N', 'init_sd_4N', 'sigma_kary')

exact_metrics <- function(context, native) {
  target <- context$case$data$prob
  prediction <- native$pred$kary
  stopifnot(identical(names(target), names(prediction)))
  tv <- vapply(names(target), function(id) {
    q <- prediction[[id]]
    sum(abs(q / sum(q) - target[[id]])) / 2
  }, numeric(1))
  loss <- native$loss
  measure <- c(D = loss$value, growth_RMS = loss$rms_growth,
               growth_max = loss$max_growth, kary_TV = max(tv))
  stopifnot(abs(unname(measure['kary_TV']) - loss$max_kary_tv) <= 1e-12)
  list(values = measure, tv = tv, violation = max(0, measure / limits - 1),
       pass = all(measure <= limits))
}

is_improved <- function(a, b) {
  if (!a$ok) return(FALSE)
  if (a$metrics$pass != b$metrics$pass) return(a$metrics$pass)
  if (a$metrics$pass || a$metrics$violation == b$metrics$violation)
    return(a$candidate$score < b$candidate$score)
  a$metrics$violation < b$metrics$violation
}

linear_model <- function(anchor, derivative, target, columns, radius) {
  active_J <- derivative[, columns, drop = FALSE] * radius
  growth_rows <- seq_len(growth_count)
  Jg <- active_J[growth_rows, , drop = FALSE] * sqrt(growth_count)
  offset <- growth_count + seq_len(kary_count * bins)
  q <- unlist(anchor$kary, use.names = FALSE)
  q0 <- unlist(target, use.names = FALSE)
  K <- active_J[offset, , drop = FALSE] * (2 * sqrt(kary_count * q))
  groups <- split(seq_along(q), rep(seq_len(kary_count), each = bins))
  r <- anchor$residual
  smooth <- 1e-8
  function(x) {
    u <- x[seq_along(columns)]
    t <- x[length(x)]
    rp <- r + as.numeric(active_J %*% u)
    d <- as.numeric(q - q0 + K %*% u)
    smooth_abs <- sqrt(d^2 + smooth^2)
    sg <- rp[growth_rows] * sqrt(growth_count)
    growth_rms <- sqrt(mean(sg^2))
    tv <- vapply(groups, function(i) sum(smooth_abs[i]) / 2, numeric(1))
    tv_j <- do.call(rbind, lapply(groups, function(i) {
      as.numeric(crossprod(K[i, , drop = FALSE], d[i] / smooth_abs[i])) / 2
    }))
    # Proposal buffer only. Exact gate uses the original .005 threshold.
    tv_bound <- .005 * .997
    g <- c(sum(rp^2) / (2 * limits['D']), growth_rms / limits['growth_RMS'],
           sg / limits['growth_max'], -sg / limits['growth_max'], tv / tv_bound) - 1 - t
    gJ <- rbind(as.numeric(crossprod(active_J, rp)) / limits['D'],
                as.numeric(crossprod(Jg, sg)) / (growth_count * max(growth_rms, 1e-16) * limits['growth_RMS']),
                Jg / limits['growth_max'], -Jg / limits['growth_max'], tv_j / tv_bound)
    list(constraints = g, jacobian = cbind(gJ, -1))
  }
}

bounded_proposal <- function(anchor, derivative, target, columns, radius) {
  unit <- anchor$candidate$unit
  model <- linear_model(anchor, derivative, target, columns, radius)
  n <- length(columns)
  z <- rep(0, n)
  start_slack <- max(0, model(c(z, 0))$constraints) + 1e-6
  objective <- function(x) list(objective = x[n + 1L] + 1e-4 * sum(x[seq_len(n)]^2),
                                gradient = c(2e-4 * x[seq_len(n)], 1))
  constraint <- function(x) {
    vals <- model(x)
    list(constraints = vals$constraints, jacobian = vals$jacobian)
  }
  result <- nloptr::nloptr(
    x0 = c(z, start_slack), eval_f = objective, eval_g_ineq = constraint,
    lb = c(pmax(-1, -unit[columns] / radius), 0),
    ub = c(pmin(1, (1 - unit[columns]) / radius), start_slack + 1),
    opts = list(algorithm = 'NLOPT_LD_SLSQP', xtol_rel = 1e-8, maxeval = 160L,
                maxtime = 0, tol_constraints_ineq = rep(1e-9, 2L + 2L * growth_count + kary_count))
  )
  stopifnot(all(is.finite(result$solution)))
  proposal <- unit
  proposal[columns] <- pmin(1, pmax(0, unit[columns] + radius * result$solution[seq_len(n)]))
  normalized_step <- (proposal[columns] - unit[columns]) / radius
  predicted <- max(0, model(c(normalized_step, 0))$constraints)
  list(unit = proposal, step = max(abs(proposal - unit)), status = result$status,
       linear_violation = predicted, original_linear_violation = start_slack - 1e-6)
}

interface_fixture <- function() {
  reference <- rep(list(c(.5, .5, rep(0, bins - 2L))), kary_count)
  fit <- reference
  fit[[1L]][1:2] <- c(.5075, .4925)
  r <- c(rep(0, growth_count), unlist(lapply(seq_along(fit), function(i) {
    (sqrt(fit[[i]]) - sqrt(reference[[i]])) / sqrt(kary_count)
  })))
  anchor <- list(residual = r, kary = fit, candidate = list(unit = rep(.5, 19)))
  j <- matrix(0, length(r), 19L)
  j[growth_count + 1:2, 1L] <- c(1,-1) / (2 * sqrt(kary_count * fit[[1L]][1:2]))
  proposal <- bounded_proposal(anchor, j, reference, 1:19, .01)
  stopifnot(proposal$step > 0, proposal$linear_violation < 1e-6,
            proposal$original_linear_violation > .49, proposal$unit[1L] < .5)
  invisible(TRUE)
}
