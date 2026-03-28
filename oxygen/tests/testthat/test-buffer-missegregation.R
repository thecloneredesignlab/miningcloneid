testthat::test_that("nullisomy oracle matches production and required invariants", {
  gamma <- 1.0
  N_unit <- 22L

  s33_prod <- cpp_o2simps_loss_survival_nullisomy(33L, 1L, gamma_loss = gamma, N_unit = N_unit)
  s33_oracle <- oracle_Sloss(33L, 1L, gamma_loss = gamma, N_unit = N_unit)
  testthat::expect_lt(s33_prod, 1.0)
  testthat::expect_equal(s33_prod, s33_oracle, tolerance = 1e-12)

  s44_prod <- cpp_o2simps_loss_survival_nullisomy(44L, 1L, gamma_loss = gamma, N_unit = N_unit)
  testthat::expect_equal(s44_prod, 1.0, tolerance = 1e-12)

  s88_1 <- cpp_o2simps_loss_survival_nullisomy(88L, 1L, gamma_loss = gamma, N_unit = N_unit)
  s88_2 <- cpp_o2simps_loss_survival_nullisomy(88L, 2L, gamma_loss = gamma, N_unit = N_unit)
  s88_3 <- cpp_o2simps_loss_survival_nullisomy(88L, 3L, gamma_loss = gamma, N_unit = N_unit)
  testthat::expect_equal(s88_1, 1.0, tolerance = 1e-12)
  testthat::expect_equal(s88_2, 1.0, tolerance = 1e-12)
  testthat::expect_equal(s88_3, 1.0, tolerance = 1e-12)

  m_grid <- 1:5
  s33_grid <- vapply(
    m_grid,
    function(m) cpp_o2simps_loss_survival_nullisomy(33L, as.integer(m), gamma_loss = gamma, N_unit = N_unit),
    numeric(1)
  )
  testthat::expect_true(all(diff(s33_grid) <= 1e-12))

  safe_33_1 <- oracle_prob_no_null_after_loss(33L, 1L, N_unit = N_unit)
  testthat::expect_gt(safe_33_1, 0.0)
  testthat::expect_lt(safe_33_1, 1.0)
  s_gamma1 <- cpp_o2simps_loss_survival_nullisomy(33L, 1L, gamma_loss = 1.0, N_unit = N_unit)
  s_gamma2 <- cpp_o2simps_loss_survival_nullisomy(33L, 1L, gamma_loss = 2.0, N_unit = N_unit)
  testthat::expect_lt(s_gamma2, s_gamma1)
})

testthat::test_that("signed daughter asymmetry is preserved in kernel and B matrix", {
  N <- 33L
  n <- 1L
  p <- 0.03
  gamma <- 1.0
  N_unit <- 22L

  delta <- cpp_o2simps_pr_delta_vec(N, p, eps_tail = 0.0, gamma_loss = gamma, N_unit = N_unit)
  w_plus <- shift_weight(delta, n)
  w_minus <- shift_weight(delta, -n)
  s_loss <- oracle_Sloss(N, n, gamma_loss = gamma, N_unit = N_unit)

  testthat::expect_gt(w_plus, 0.0)
  testthat::expect_lt(w_minus, w_plus)
  testthat::expect_equal(w_minus / w_plus, s_loss, tolerance = 1e-10)

  Nmin <- 30L
  Nmax <- 36L
  R <- Nmax - Nmin + 1L
  p_vec <- rep(0.0, R)
  p_vec[N - Nmin + 1L] <- p
  tri <- cpp_o2simps_build_B_total_triplet(
    Nmin, Nmax, p_vec,
    boundary = "drop",
    eps_tail = 0.0,
    gamma_loss = gamma,
    N_unit = N_unit
  )
  B <- triplet_to_sparse(tri)
  col_idx <- N - Nmin + 1L

  testthat::expect_equal(as.numeric(B[N + n - Nmin + 1L, col_idx]), w_plus, tolerance = 1e-12)
  testthat::expect_equal(as.numeric(B[N - n - Nmin + 1L, col_idx]), w_minus, tolerance = 1e-12)
})

testthat::test_that("per-division offspring accounting is conserved", {
  N <- 33L
  p <- 0.07
  gamma <- 0.8
  N_unit <- 22L

  delta <- cpp_o2simps_pr_delta_vec(N, p, eps_tail = 0.0, gamma_loss = gamma, N_unit = N_unit)
  live_model <- sum(as.numeric(delta$prob))
  dropped_daughters_model <- 2.0 * as.numeric(delta$mass_dropped)

  live_oracle <- oracle_live_mass_per_division(N, p, gamma_loss = gamma, N_unit = N_unit, eps_tail = 0.0)
  dropped_daughters_oracle <- 2.0 - live_oracle

  testthat::expect_equal(live_model + dropped_daughters_model, 2.0, tolerance = 1e-12)
  testthat::expect_equal(live_model, live_oracle, tolerance = 1e-10)
  testthat::expect_equal(dropped_daughters_model, dropped_daughters_oracle, tolerance = 1e-10)
})

testthat::test_that("zero-missegregation limit keeps all daughters on self state", {
  N <- 35L
  N_unit <- 22L
  delta <- cpp_o2simps_pr_delta_vec(N, p = 0.0, eps_tail = 0.0, gamma_loss = 1.0, N_unit = N_unit)

  testthat::expect_equal(shift_weight(delta, 0L), 2.0, tolerance = 1e-12)
  off_idx <- which(as.integer(delta$ts) != 0L)
  testthat::expect_equal(sum(as.numeric(delta$prob[off_idx])), 0.0, tolerance = 1e-12)
  testthat::expect_equal(as.numeric(delta$mass_dropped), 0.0, tolerance = 1e-12)

  Nmin <- 30L
  Nmax <- 36L
  R <- Nmax - Nmin + 1L
  tri <- cpp_o2simps_build_B_total_triplet(
    Nmin, Nmax, p_vec = rep(0.0, R),
    boundary = "drop",
    eps_tail = 0.0,
    gamma_loss = 1.0,
    N_unit = N_unit
  )
  B <- as.matrix(triplet_to_sparse(tri))
  testthat::expect_equal(B, diag(2.0, R), tolerance = 1e-12)
})

testthat::test_that("boundary drop vs absorb_minmax only differs by out-of-grid routing", {
  Nmin <- 30L
  Nmax <- 36L
  N <- 35L
  p <- 0.20
  gamma <- 1.0
  N_unit <- 22L

  R <- Nmax - Nmin + 1L
  p_vec <- rep(0.0, R)
  p_vec[N - Nmin + 1L] <- p

  tri_drop <- cpp_o2simps_build_B_total_triplet(
    Nmin, Nmax, p_vec,
    boundary = "drop",
    eps_tail = 0.0,
    gamma_loss = gamma,
    N_unit = N_unit
  )
  tri_absorb <- cpp_o2simps_build_B_total_triplet(
    Nmin, Nmax, p_vec,
    boundary = "absorb_minmax",
    eps_tail = 0.0,
    gamma_loss = gamma,
    N_unit = N_unit
  )

  B_drop <- triplet_to_sparse(tri_drop)
  B_absorb <- triplet_to_sparse(tri_absorb)
  col_idx <- N - Nmin + 1L
  v_drop <- as.numeric(B_drop[, col_idx])
  v_absorb <- as.numeric(B_absorb[, col_idx])

  delta <- cpp_o2simps_pr_delta_vec(N, p, eps_tail = 0.0, gamma_loss = gamma, N_unit = N_unit)
  expected_drop <- rep(0.0, R)
  expected_absorb <- rep(0.0, R)
  for (k in seq_along(delta$ts)) {
    t <- as.integer(delta$ts[[k]])
    w <- as.numeric(delta$prob[[k]])
    Np <- N + t
    if (Np >= Nmin && Np <= Nmax) {
      idx <- Np - Nmin + 1L
      expected_drop[idx] <- expected_drop[idx] + w
      expected_absorb[idx] <- expected_absorb[idx] + w
    } else {
      idx <- min(max(Np, Nmin), Nmax) - Nmin + 1L
      expected_absorb[idx] <- expected_absorb[idx] + w
    }
  }

  testthat::expect_equal(v_drop, expected_drop, tolerance = 1e-12)
  testthat::expect_equal(v_absorb, expected_absorb, tolerance = 1e-12)

  interior_states <- (Nmin + 1L):(Nmax - 1L)
  interior_idx <- interior_states - Nmin + 1L
  testthat::expect_equal(v_drop[interior_idx], v_absorb[interior_idx], tolerance = 1e-12)
  max_idx <- Nmax - Nmin + 1L
  testthat::expect_gt(v_absorb[max_idx], v_drop[max_idx])
  testthat::expect_lt(sum(v_drop), sum(v_absorb))
})

testthat::test_that("loss-only penalty invariance panel", {
  gamma <- 1.0
  p <- 0.05
  N_unit <- 22L
  panel <- data.frame(
    N = c(44L, 33L, 88L, 88L),
    n = c(1L, 1L, 3L, 4L)
  )

  for (idx in seq_len(nrow(panel))) {
    N <- panel$N[[idx]]
    n <- panel$n[[idx]]
    delta <- cpp_o2simps_pr_delta_vec(N, p, eps_tail = 0.0, gamma_loss = gamma, N_unit = N_unit)
    w_plus <- shift_weight(delta, n)
    w_minus <- shift_weight(delta, -n)
    s_loss <- oracle_Sloss(N, n, gamma_loss = gamma, N_unit = N_unit)

    testthat::expect_gt(w_plus, 0.0)
    if (abs(s_loss - 1.0) < 1e-12) {
      testthat::expect_equal(w_minus, w_plus, tolerance = 1e-10)
    } else {
      testthat::expect_lt(w_minus, w_plus)
      testthat::expect_equal(w_minus / w_plus, s_loss, tolerance = 1e-8)
    }
  }
})

testthat::test_that("boundary=drop routes out-of-grid offspring into dead buffer rate", {
  N <- 44L
  p <- 0.10
  gamma <- 1.0
  N_unit <- 22L

  delta <- cpp_o2simps_pr_delta_vec(N, p, eps_tail = 0.0, gamma_loss = gamma, N_unit = N_unit)
  expected_dead_buffer_rate <- 2.0 * as.numeric(delta$mass_dropped) +
    sum(as.numeric(delta$prob[as.integer(delta$ts) != 0L]))

  tri_g <- cpp_o2simps_build_G_for_o2_triplet(
    O2 = 1.0,
    O2_crit = 1.0,
    N0min = N,
    N0max = N,
    N1min = N,
    N1max = N,
    lam_min = 1.0,
    lam_max = 1.0,
    k_o = 1.0,
    has_p_misseg = TRUE,
    p_mis_base = p,
    p_misseg = 0.0,
    k_o_mis = 1.0,
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0.0,
    pmis_O2_1 = 0.0,
    p_const = 0.0,
    p_wgd = 0.0,
    boundary = "drop",
    eps_tail = 0.0,
    gamma_loss = gamma,
    N_unit = N_unit,
    beta_size = 0.0,
    alpha_o2 = 0.0,
    gamma_growth = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    n_O = 1.0,
    ploidy_O2_death = "ploidy_related"
  )

  testthat::expect_equal(
    as.numeric(tri_g$dead_buffer_rate[[1]]),
    expected_dead_buffer_rate,
    tolerance = 1e-10
  )
})
