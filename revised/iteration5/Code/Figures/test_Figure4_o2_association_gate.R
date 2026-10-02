#!/usr/bin/env Rscript
# Fast regression checks; no bootstrap, figure generation, or file writes.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
source(normalizePath(args[[1L]], mustWork = TRUE))
check <- function(delta, q, peak, group, passes) {
  result <- figure4_o2_association_decision(delta, q, peak)
  stopifnot(identical(result$group, group), identical(result$peak_passes, passes),
            identical(result$minimum_peak, 0.3))
}
check(.2, .001, .254246771, "O2-independent", FALSE) # k_o_mis
check(.2, .001, .428290065, "Low O2", TRUE) # gamma_mu
check(-.2, .001, .4, "High O2", TRUE)
check(.2, .001, .3, "O2-independent", FALSE) # strict > boundary
check(.2, .001, .3 + 1e-10, "Low O2", TRUE)
check(.2, .05, .4, "O2-independent", TRUE) # q must be strict < .05
check(.2, .2, .8, "O2-independent", TRUE)
check(-.2, .001, .1, "O2-independent", FALSE) # gate applies to High too
check(0, .5, .4, "O2-independent", TRUE)
stopifnot(inherits(try(figure4_o2_association_decision(0, .001, .4), silent = TRUE), "try-error"),
          inherits(try(figure4_o2_association_decision(.2, NA_real_, .4), silent = TRUE), "try-error"),
          inherits(try(figure4_o2_association_decision(.2, .001, 1.1), silent = TRUE), "try-error"))
cat("figure4_o2_global_peak_gate_unit_tests_ok\n")
