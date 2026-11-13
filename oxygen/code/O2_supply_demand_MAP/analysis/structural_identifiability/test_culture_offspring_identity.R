#!/usr/bin/env Rscript
# Independent PMF-sum oracle for the diagnostic identity, before numerical loss.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: test_culture_offspring_identity.R <diagnostic-helper.R>")
source(args[1])
cases <- 0L
for (N in c(22L, 44L, 88L, 154L)) for (p in c(0, 1e-6, .01, .4, 1)) for (s in c(0, .2, .8, 1)) {
  m <- 0:N
  survival <- s^m; survival[1] <- 1
  weights <- dbinom(m, N, p) * survival
  total <- sum(weights)
  actual <- pilotdiag_offspring_identity(N, p, s)
  stopifnot(abs(actual$biological_viable_daughters_per_ordinary_division - 2 * total) <= 1e-11 * max(total, 1e-300))
  if (total > 0) {
    conditional <- weights / total
    expected <- dbinom(m, N, actual$effective_error_probability_conditional_on_viability)
    stopifnot(max(abs(conditional - expected)) < 1e-10)
  } else stopifnot(is.na(actual$effective_error_probability_conditional_on_viability))
  cases <- cases + 1L
}
cat("Binomial-survival diagnostic identity passed", cases, "PMF-sum cases; no simulations or optimization.\n")
