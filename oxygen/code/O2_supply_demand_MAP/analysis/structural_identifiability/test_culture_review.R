# Necessary algorithm and execution-policy tests; no biological fitting.
stopifnot(requireNamespace("testthat", quietly = TRUE), requireNamespace("DEoptim", quietly = TRUE))
testthat::test_that("bounded pattern search handles coupled nonsmooth objectives", {
  directions <- rbind(diag(2), c(1, 1), c(1, -1))
  fn <- function(x) abs(x[1] - .3) + abs(x[2] - .7)
  x <- o2review_pattern(fn, c(.9, .1), directions, 5000L)
  testthat::expect_lt(x$value, 1e-4)
  testthat::expect_true(x$tested_directional_stability)
  testthat::expect_true(all(x$par >= 0 & x$par <= 1))
})
testthat::test_that("evaluation and wall limits cannot masquerade as convergence", {
  x <- o2review_pattern(function(x) sum(x^2), c(.8, .8), diag(2), 3L)
  testthat::expect_true(x$budget_limited)
  testthat::expect_false(x$tested_directional_stability)
  testthat::expect_lte(x$evaluations, 3L)
  counter <- 0
  clock <- function() { counter <<- counter + 1; counter }
  x <- o2review_pattern(function(x) sum(x), c(.5, .5), diag(2), deadline = 5, clock = clock)
  testthat::expect_true(x$budget_limited)
})
testthat::test_that("physical boundary minima pass one-sided stability only", {
  x <- o2review_pattern(function(x) sum(x), c(.1, .1), diag(2), 5000L)
  testthat::expect_equal(x$par, c(0, 0), tolerance = 1e-10)
  testthat::expect_true(x$tested_directional_stability)
})
testthat::test_that("global-only interception performs no local optimization", {
  counter <- new.env(); counter$calls <- 0L
  calls <- 0L
  stub <- o2review_local_stub(counter)
  answer <- stub(c(1, 2), function(p) { calls <<- calls + 1L; sum(p^2) }, method = "L-BFGS-B")
  testthat::expect_equal(counter$calls, 1L)
  testthat::expect_equal(calls, 1L)
  testthat::expect_equal(answer$par, c(1, 2))
  testthat::expect_true(is.na(answer$convergence))
})
testthat::test_that("the real DE early-stop path respects the enforced minimum", {
  task <- list(min_generations = 300L, max_generations = 350L)
  controls <- o2review_de_controls(task)
  # Real DE, with a constant objective that would otherwise stop immediately.
  controls$NP <- 20L; controls$trace <- FALSE; controls$parallelType <- "none"
  set.seed(41)
  fit <- DEoptim::DEoptim(function(x) 1, c(0, 0), c(1, 1), control = controls)
  testthat::expect_gte(fit$optim$iter, task$min_generations)
  testthat::expect_lte(fit$optim$iter, task$max_generations)
  testthat::expect_equal(o2review_de_controls(list(min_generations = 300L, max_generations = 1000L))$steptol, 300L)
})
testthat::test_that("zero constraints and sensitivity domains have separate maps", {
  if (exists("ctx", inherits = TRUE)) {
    ctx_test <- get("ctx", inherits = TRUE)
    p <- ctx_test$saved$best_params
    for (target in c("p_misseg", "buffer_beta")) {
      fixed <- setNames(list(0), target)
      m <- o2review_map(ctx_test, "original", fixed)
      x <- o2prof_project_start(m, p)$params
      testthat::expect_identical(x[[target]], 0)
      testthat::expect_identical(o2prof_from_free(m, o2prof_to_free(m, x))[[target]], 0)
    }
    original <- o2review_map(ctx_test)$table
    for (domain in c("D_base", "D_sat")) {
      t <- o2review_map(ctx_test, domain)$table
      changed <- which(t$upper != original$upper)
      testthat::expect_length(changed, 1L)
      testthat::expect_equal(t$upper[changed], if (domain == "D_base") .02 else .2)
      testthat::expect_identical(t$lower, original$lower)
    }
  }
})
cat("Review algorithm and execution tests passed.\n")
