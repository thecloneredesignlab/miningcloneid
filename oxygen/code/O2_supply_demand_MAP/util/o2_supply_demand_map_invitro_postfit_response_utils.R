# Shared fitted-response engine used by the in-vitro O2 and ploidy simulation
# modules.
# It evaluates model rates but does not choose or write domain-specific tables.

ivt_sim_build_functional_response_context <- function(run_params, cfg) {
  start_mode <- assert_canonical_start_with_mode(
    ivt_sim_first_non_null(cfg$start_with, "ploidy")
  )
  o2_crit <- as.numeric(ivt_sim_first_non_null(run_params$O2_crit, cfg$o2_crit_init, 1))
  if (!is.finite(o2_crit) || o2_crit < 0) o2_crit <- 1
  o2_growth <- isTRUE(ivt_sim_first_non_null(cfg$O2_growth, TRUE))
  ploidy_death <- assert_canonical_ploidy_o2_death_mode(
    ivt_sim_first_non_null(
      cfg$ploidy_O2_death,
      run_params$ploidy_O2_death,
      "diploid_NULL"
    )
  )
  n_o <- as.numeric(ivt_sim_first_non_null(run_params$n_O, cfg$n_O_init, 1))
  mu_hp <- pmax(
    as.numeric(ivt_sim_first_non_null(run_params$mu_hp, cfg$mu_hp_init, 1e-3)),
    0
  )
  gamma_mu <- pmax(
    as.numeric(ivt_sim_first_non_null(run_params$gamma_mu, cfg$gamma_mu_init, 2)),
    1e-12
  )
  boundary <- as.character(ivt_sim_first_non_null(run_params$boundary, "drop"))
  p_wgd <- as.numeric(ivt_sim_first_non_null(run_params$p_wgd, 0))
  buffer_smax <- as.numeric(
    ivt_sim_first_non_null(run_params$buffer_smax, cfg$buffer_smax_init, 1)
  )
  buffer_beta <- as.numeric(
    ivt_sim_first_non_null(run_params$buffer_beta, cfg$buffer_beta_init, 0)
  )
  buffer_n_exp <- as.numeric(
    ivt_sim_first_non_null(run_params$buffer_n_exp, cfg$buffer_n_exp_init, 1)
  )
  n_state <- as.integer(cfg$N_MAX) - as.integer(cfg$N_MIN) + 1L
  .require_cpp_o2simps_fn("cpp_o2simps_build_G_for_o2_triplet")
  cache <- new.env(parent = emptyenv())

  rate_curves_at_o2 <- function(o2) {
    o2 <- .assert_o2_pct(as.numeric(o2), label = "O2")
    key <- sprintf("%.6f", o2)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      tri <- cpp_o2simps_build_G_for_o2_triplet(
        O2 = o2,
        O2_crit = o2_crit,
        N0min = as.integer(cfg$N_MIN),
        N0max = as.integer(cfg$N_MAX),
        N1min = as.integer(cfg$N_MIN),
        N1max = as.integer(cfg$N_MAX),
        lam_max = as.numeric(run_params$lam_max),
        p_mis_base = as.numeric(ivt_sim_first_non_null(
          run_params$p_mis_base,
          cfg$p_mis_base,
          cfg$p_mis_base_init,
          1e-5
        )),
        p_misseg = as.numeric(ivt_sim_first_non_null(run_params$p_misseg, 0)),
        k_o_mis = as.numeric(ivt_sim_first_non_null(run_params$k_o_mis, 50)),
        p_wgd = p_wgd,
        boundary = boundary,
        eps_tail = 1e-8,
        buffer_smax = buffer_smax,
        buffer_beta = buffer_beta,
        buffer_n_exp = buffer_n_exp,
        N_unit = as.integer(cfg$N_UNIT),
        beta_size = 0,
        O2_growth = o2_growth,
        alpha_o2 = as.numeric(ivt_sim_first_non_null(
          run_params$alpha_o2,
          cfg$alpha_o2_init,
          0.5
        )),
        gamma_growth = as.numeric(ivt_sim_first_non_null(
          run_params$gamma_growth,
          cfg$gamma_growth_init,
          2
        )),
        mu_hp = mu_hp,
        gamma_mu = gamma_mu,
        n_O = n_o,
        ploidy_O2_death = ploidy_death
      )
      names_needed <- c(
        "dead_buffer_rate",
        "misseg_nonviable_rate",
        "boundary_dropped_rate",
        "misseg_nonviable_division_prob",
        "misseg_nonviable_daughters_per_division"
      )
      curves <- lapply(names_needed, function(name) {
        value <- as.numeric(tri[[name]])
        if (length(value) != n_state) {
          stop("Unexpected ", name, " length from model core.")
        }
        value
      })
      names(curves) <- names_needed
      curves$misseg_nonviable_daughter_fraction <- pmax(
        pmin(0.5 * curves$misseg_nonviable_daughters_per_division, 1),
        0
      )
      assign(key, curves, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }

  compute_rates <- function(o2, n) {
    o2 <- rep_len(
      .assert_o2_pct(as.numeric(o2), label = "O2"),
      max(length(o2), length(n))
    )
    n <- rep_len(as.numeric(n), length(o2))
    proliferation <- as.numeric(.lambda_eff_of_O2(
      O2 = o2,
      run_params = run_params,
      N = n,
      O2_crit = o2_crit,
      O2_growth = o2_growth
    ))
    death <- as.numeric(.mu_eff_of_O2(
      O2 = o2,
      run_params = run_params,
      N = n,
      O2_crit = o2_crit
    ))
    extras <- matrix(NA_real_, nrow = length(o2), ncol = 6L)
    colnames(extras) <- c(
      "dead_buffer_rate",
      "misseg_nonviable_rate",
      "boundary_dropped_rate",
      "misseg_nonviable_division_prob",
      "misseg_nonviable_daughters_per_division",
      "misseg_nonviable_daughter_fraction"
    )
    groups <- split(seq_along(o2), sprintf("%.6f", o2))
    for (index in groups) {
      curves <- rate_curves_at_o2(o2[[index[[1L]]]])
      state_index <- as.integer(round(n[index])) - as.integer(cfg$N_MIN) + 1L
      valid <- is.finite(state_index) & state_index >= 1L & state_index <= n_state
      if (any(valid)) {
        for (name in colnames(extras)) {
          extras[index[valid], name] <- curves[[name]][state_index[valid]]
        }
      }
    }
    proliferation_rate <- pmax(proliferation, 0)
    death_rate <- pmax(death, 0)
    buffer_death_rate <- pmax(extras[, "dead_buffer_rate"], 0)
    # Match the C++ live-state balance: dead-buffer flow is lost from live cells.
    data.frame(
      O2 = o2,
      N = n,
      proliferation_rate = proliferation_rate,
      death_rate = death_rate,
      buffer_death_rate = buffer_death_rate,
      buffer_death_per_division = buffer_death_rate /
        pmax(proliferation_rate, 1e-12),
      misseg_nonviable_rate = pmax(extras[, "misseg_nonviable_rate"], 0),
      boundary_dropped_rate = pmax(extras[, "boundary_dropped_rate"], 0),
      misseg_nonviable_division_prob = pmax(
        pmin(extras[, "misseg_nonviable_division_prob"], 1),
        0
      ),
      misseg_nonviable_daughters_per_division = pmax(
        pmin(extras[, "misseg_nonviable_daughters_per_division"], 2),
        0
      ),
      misseg_nonviable_daughter_fraction = pmax(
        pmin(extras[, "misseg_nonviable_daughter_fraction"], 1),
        0
      ),
      net_growth_rate = proliferation_rate - death_rate - buffer_death_rate,
      stringsAsFactors = FALSE
    )
  }

  list(
    start_mode = start_mode,
    o2_crit = o2_crit,
    n_unit = as.numeric(cfg$N_UNIT),
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp,
    compute_rates = compute_rates
  )
}
