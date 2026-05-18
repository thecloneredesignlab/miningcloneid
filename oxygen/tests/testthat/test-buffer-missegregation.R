testthat::test_that("buffering kernel applies symmetric daughter survival", {
  N <- 33L
  n <- 1L
  p <- 0.03
  N_unit <- 22L
  buffer_smax <- 0.8
  buffer_beta <- 0.5
  buffer_n_exp <- 1.2

  delta <- cpp_o2simps_pr_delta_vec(
    N,
    p,
    eps_tail = 0.0,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp,
    N_unit = N_unit
  )
  w_plus <- shift_weight(delta, n)
  w_minus <- shift_weight(delta, -n)
  s_buf <- oracle_buffering_survival(
    N,
    n,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp,
    N_unit = N_unit
  )

  testthat::expect_gt(w_plus, 0.0)
  testthat::expect_equal(w_minus, w_plus, tolerance = 1e-12)
  testthat::expect_equal(w_plus / stats::dbinom(n, N, p), s_buf, tolerance = 1e-10)

  Nmin <- 30L
  Nmax <- 36L
  R <- Nmax - Nmin + 1L
  p_vec <- rep(0.0, R)
  p_vec[N - Nmin + 1L] <- p
  tri <- cpp_o2simps_build_B_total_triplet(
    Nmin, Nmax, p_vec,
    boundary = "drop",
    eps_tail = 0.0,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp,
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
  N_unit <- 22L
  buffer_smax <- 0.85
  buffer_beta <- 0.4
  buffer_n_exp <- 1.1

  delta <- cpp_o2simps_pr_delta_vec(
    N,
    p,
    eps_tail = 0.0,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp,
    N_unit = N_unit
  )
  live_model <- sum(as.numeric(delta$prob))
  dropped_daughters_model <- 2.0 * as.numeric(delta$mass_dropped)

  live_oracle <- oracle_live_mass_per_division(
    N,
    p,
    N_unit = N_unit,
    eps_tail = 0.0,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp
  )
  dropped_daughters_oracle <- 2.0 - live_oracle

  testthat::expect_equal(live_model + dropped_daughters_model, 2.0, tolerance = 1e-12)
  testthat::expect_equal(live_model, live_oracle, tolerance = 1e-10)
  testthat::expect_equal(dropped_daughters_model, dropped_daughters_oracle, tolerance = 1e-10)
})

testthat::test_that("zero-missegregation limit keeps all daughters on self state", {
  N <- 35L
  N_unit <- 22L
  delta <- cpp_o2simps_pr_delta_vec(N, p = 0.0, eps_tail = 0.0, N_unit = N_unit)

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
  N_unit <- 22L

  R <- Nmax - Nmin + 1L
  p_vec <- rep(0.0, R)
  p_vec[N - Nmin + 1L] <- p

  tri_drop <- cpp_o2simps_build_B_total_triplet(
    Nmin, Nmax, p_vec,
    boundary = "drop",
    eps_tail = 0.0,
    N_unit = N_unit
  )
  tri_absorb <- cpp_o2simps_build_B_total_triplet(
    Nmin, Nmax, p_vec,
    boundary = "absorb_minmax",
    eps_tail = 0.0,
    N_unit = N_unit
  )

  B_drop <- triplet_to_sparse(tri_drop)
  B_absorb <- triplet_to_sparse(tri_absorb)
  col_idx <- N - Nmin + 1L
  v_drop <- as.numeric(B_drop[, col_idx])
  v_absorb <- as.numeric(B_absorb[, col_idx])

  delta <- cpp_o2simps_pr_delta_vec(N, p, eps_tail = 0.0, N_unit = N_unit)
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

testthat::test_that("boundary=drop routes out-of-grid offspring into dead buffer rate", {
  N <- 44L
  p <- 0.10
  N_unit <- 22L

  delta <- cpp_o2simps_pr_delta_vec(N, p, eps_tail = 0.0, N_unit = N_unit)
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

testthat::test_that("buffering dead-buffer rate preserves expected dead daughters per division above one", {
  N <- 88L
  p <- 0.12
  N_unit <- 22L
  buffer_smax <- 0.804673298774287
  buffer_beta <- 0.304512642053469
  buffer_n_exp <- 7.80109879059526

  delta <- cpp_o2simps_pr_delta_vec(
    N,
    p,
    eps_tail = 0.0,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp,
    N_unit = N_unit
  )
  expected_dead_daughters_per_division <- 2.0 * as.numeric(delta$mass_dropped)
  testthat::expect_gt(expected_dead_daughters_per_division, 1.0)

  tri_g <- cpp_o2simps_build_G_for_o2_triplet(
    O2 = 1.0,
    O2_crit = 1.0,
    N0min = 0L,
    N0max = 220L,
    N1min = 0L,
    N1max = 220L,
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
    p_wgd_max = 0.0,
    O2_wgd = 0.1,
    boundary = "drop",
    eps_tail = 0.0,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp,
    N_unit = N_unit,
    beta_size = 0.0,
    O2_growth = TRUE,
    alpha_o2 = 0.0,
    gamma_growth = 1.0,
    mu_hp = 0.0,
    gamma_mu = 1.0,
    n_O = 1.0,
    ploidy_O2_death = "ploidy_related"
  )
  idx <- N + 1L

  testthat::expect_equal(
    as.numeric(tri_g$misseg_nonviable_daughters_per_division[[idx]]),
    expected_dead_daughters_per_division,
    tolerance = 1e-10
  )
  testthat::expect_equal(
    as.numeric(tri_g$dead_buffer_rate[[idx]]),
    expected_dead_daughters_per_division,
    tolerance = 1e-10
  )
})

testthat::test_that("glucose-enabled growth uses coupled O2 resource stress", {
  run_params_base <- list(
    lam_min = 0.04,
    lam_max = 0.48,
    k_o = 2.0,
    alpha_o2 = 1.7,
    gamma_growth = 2.2,
    O2_crit = 1.3,
    n_O = 1.4
  )
  O2 <- c(0.6, 1.8)
  N <- c(44, 88)

  static_coupled <- .lambda_eff_of_O2(
    O2 = O2,
    run_params = c(
      run_params_base,
      list(glucose = TRUE)
    ),
    N = N,
    O2_growth = TRUE
  )
  oxygen_only <- .lambda_eff_of_O2(
    O2 = O2,
    run_params = c(
      run_params_base,
      list(glucose = FALSE)
    ),
    N = N,
    O2_growth = TRUE
  )
  coupled_q3 <- .lambda_eff_of_O2(
    O2 = O2,
    run_params = c(
      run_params_base,
      list(glucose = TRUE, qc = 3)
    ),
    N = N,
    O2_growth = TRUE
  )
  oxygen_only_q5 <- .lambda_eff_of_O2(
    O2 = O2,
    run_params = c(
      run_params_base,
      list(glucose = FALSE, qc = 5)
    ),
    N = N,
    O2_growth = TRUE
  )

  testthat::expect_gt(max(abs(static_coupled - oxygen_only)), 1e-4)
  testthat::expect_gt(max(abs(coupled_q3 - static_coupled)), 1e-4)
  testthat::expect_equal(oxygen_only_q5, oxygen_only, tolerance = 1e-12)
})

testthat::test_that("C++ glucose generator uses coupled O2 resource growth", {
  N <- 66L
  O2 <- 0.9
  O2_crit <- 1.4
  n_O <- 1.6
  lam_min <- 0.03
  lam_max <- 0.52
  alpha_o2 <- 1.5
  gamma_growth <- 2.0

  generator_diag <- function(glucose, qc = 2) {
    tri <- cpp_o2simps_build_G_for_o2_triplet(
      O2 = O2,
      O2_crit = O2_crit,
      N0min = N,
      N0max = N,
      N1min = N,
      N1max = N,
      lam_min = lam_min,
      lam_max = lam_max,
      k_o = 2.0,
      has_p_misseg = TRUE,
      p_mis_base = 0.0,
      p_misseg = 0.0,
      k_o_mis = 1.0,
      has_pmis_endpoints = FALSE,
      pmis_O2_0 = 0.0,
      pmis_O2_1 = 0.0,
      p_const = 0.0,
      glucose = glucose,
      p_wgd = 0.0,
      p_wgd_max = 0.0,
      O2_wgd = 0.1,
      boundary = "drop",
      eps_tail = 0.0,
      buffer_smax = 0.8,
      buffer_beta = 1.0,
      buffer_n_exp = 1.0,
      N_unit = 22L,
      beta_size = 0.0,
      O2_growth = TRUE,
      alpha_o2 = alpha_o2,
      gamma_growth = gamma_growth,
      mu_hp = 0.0,
      gamma_mu = 1.0,
      n_O = n_O,
      qc = qc,
      ploidy_O2_death = "ploidy_related"
    )
    as.numeric(triplet_to_sparse(tri)[1, 1])
  }

  h_o2 <- (O2_crit^n_O) / (O2_crit^n_O + O2^n_O)
  R <- (1 - h_o2)^2
  h_resource <- 1 - R
  expected <- (lam_min + (lam_max - lam_min) * R) /
    (1 + alpha_o2 * h_resource * (N / 44)^gamma_growth)
  R_q1 <- 1 - h_o2
  expected_q1 <- (lam_min + (lam_max - lam_min) * R_q1) /
    (1 + alpha_o2 * h_o2 * (N / 44)^gamma_growth)

  static_coupled <- generator_diag(glucose = TRUE)
  coupled_q1 <- generator_diag(glucose = TRUE, qc = 1)
  oxygen_only <- generator_diag(glucose = FALSE)
  oxygen_only_q5 <- generator_diag(glucose = FALSE, qc = 5)

  testthat::expect_equal(static_coupled, expected, tolerance = 1e-12)
  testthat::expect_equal(coupled_q1, expected_q1, tolerance = 1e-12)
  testthat::expect_gt(abs(static_coupled - oxygen_only), 1e-4)
  testthat::expect_equal(oxygen_only_q5, oxygen_only, tolerance = 1e-12)
})
