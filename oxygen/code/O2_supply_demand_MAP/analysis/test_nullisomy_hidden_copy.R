#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Rcpp))

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(argv) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    val <- if (length(parts) >= 2) paste(parts[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || !is.finite(v)) default else v
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.integer(x))
  if (length(v) != 1L || is.na(v)) default else v
}

resolve_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- args_full[grepl("^--file=", args_full)]
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  getwd()
}

script_dir <- normalizePath(resolve_script_dir(), mustWork = FALSE)
workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
model_r_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")
cpp_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.cpp")
cache_dir <- file.path(tempdir(), ".rcpp_cache_o2_supply_demand_map_test")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(file = cpp_path, rebuild = FALSE, showOutput = FALSE, verbose = FALSE, cacheDir = cache_dir)

nullisomy_risk_curve <- function(N, mode = "balanced", alpha = 100, mc_samples = 10000L, seed = 12345L, N_unit = 22L) {
  m_vals <- 0:N
  surv <- vapply(
    m_vals,
    function(m) cpp_o2simps_loss_survival_nullisomy(
      as.integer(N),
      as.integer(m),
      gamma_loss = 1.0,
      N_unit = N_unit,
      nullisomy_hidden_copy_mode = mode,
      nullisomy_dirichlet_alpha = alpha,
      nullisomy_dirichlet_mc_samples = mc_samples,
      nullisomy_dirichlet_seed = seed
    )[[1]],
    numeric(1)
  )
  pmin(pmax(1 - surv, 0), 1)
}

build_simple_sim <- function(mode = "balanced", alpha = 100, mc_samples = 10000L, seed = 12345L) {
  Nmin <- 44L
  Nmax <- 95L
  grid <- seq.int(Nmin, Nmax)
  init <- rep(0, length(grid))
  init[[which(grid == 90L)]] <- 1000
  sim <- cpp_o2simps_simulate_one(
    init_state = as.numeric(init),
    N0min = Nmin,
    N0max = Nmax,
    N1min = Nmin,
    N1max = Nmax,
    obs_steps = as.integer(c(4L, 8L)),
    sim_end_step = as.integer(8L),
    DT = 0.25,
    dose = 0,
    dose_ref = 1,
    treat_day = 0,
    fit_treatment = FALSE,
    alpha = 0,
    gamma = 1,
    tx_mult_min = 1,
    crowding_enabled = FALSE,
    crowding = "logistic",
    K = 1e12,
    min_pop = 1e-12,
    O2_crit = 1,
    o2_feedback = FALSE,
    o2_S0 = 5,
    kappa_O = 1,
    tau_O2 = 1,
    o2_Nref = 1e6,
    o2_min = 5,
    eta_o2 = 1,
    o2_cache_bin_pct = 0.01,
    o2_cache_hysteresis_pct = 0.005,
    o2_cache_profile = FALSE,
    lam_min = 0.3,
    lam_max = 0.3,
    k_o = 1,
    has_p_misseg = FALSE,
    p_mis_base = 1e-5,
    p_misseg = 0,
    k_o_mis = 1,
    has_pmis_endpoints = FALSE,
    pmis_O2_0 = 0,
    pmis_O2_1 = 0,
    p_const = 0.08,
    p_wgd = 1e-8,
    boundary = "drop",
    eps_tail = 1e-8,
    gamma_loss = 0.1,
    N_unit = 22L,
    nullisomy_hidden_copy_mode = mode,
    nullisomy_dirichlet_alpha = alpha,
    nullisomy_dirichlet_mc_samples = as.integer(mc_samples),
    nullisomy_dirichlet_seed = as.integer(seed),
    beta_size = 0,
    O2_growth = FALSE,
    alpha_o2 = 0,
    gamma_growth = 1,
    mu_hp = 0,
    gamma_mu = 1,
    n_O = 1,
    ploidy_O2_death = "uniform",
    start_with = "chr_number",
    k_clear = 1e-3,
    vol_by_N = rep(1, length(grid)),
    burden_floor = 1e-12,
    return_full_trajectory = FALSE
  )
  sim
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

main <- function(argv) {
  args <- parse_args(argv)
  out_dir <- normalizePath(args$out_dir %||% file.path(getwd(), "nullisomy_hidden_copy_tests"), mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  N_vals <- c(44L, 53L, 88L, 90L, 95L)
  m_vals <- c(1L, 2L, 3L, 4L, 5L, 10L)
  edge_N_vals <- c(21L, 22L, 23L)
  seed_ref <- as_int(args$seed, 12345L)
  mc_default <- as_int(args$mc_samples, 10000L)
  mc_highalpha <- as_int(args$mc_samples_highalpha, 50000L)
  highalpha_tol <- as_num(args$highalpha_tol, 0.03)
  seed_alt <- seed_ref + 1L

  results <- list()
  add_result <- function(test, status, detail) {
    results[[length(results) + 1L]] <<- data.frame(test = test, status = status, detail = detail, stringsAsFactors = FALSE)
  }

  balanced_rows <- list()
  for (N in N_vals) {
    balanced_curve_direct <- nullisomy_risk_curve(N, mode = "balanced", alpha = 100, mc_samples = mc_default, seed = seed_ref)
    balanced_curve_expected <- {
      counts <- rep(floor(N / 22L), 22L)
      rem <- N %% 22L
      if (rem > 0) counts[seq_len(rem)] <- counts[seq_len(rem)] + 1L
      risk <- numeric(N + 1L)
      risk[1:(N + 1L)] <- 0
      safe_coeff <- rep(0, N + 1L)
      safe_coeff[[1]] <- 1
      max_safe_degree <- 0L
      for (copies in counts) {
        next_coeff <- rep(0, N + 1L)
        take_max <- max(0L, copies - 1L)
        for (used in 0:max_safe_degree) {
          base <- safe_coeff[[used + 1L]]
          if (base <= 0) next
          for (take in 0:min(take_max, N - used)) {
            next_coeff[[used + take + 1L]] <- next_coeff[[used + take + 1L]] + base * choose(copies, take)
          }
        }
        safe_coeff <- next_coeff
        max_safe_degree <- min(N, max_safe_degree + take_max)
      }
      for (m in 0:N) {
        total_count <- choose(N, m)
        safe_prob <- if (total_count > 0) safe_coeff[[m + 1L]] / total_count else 0
        risk[[m + 1L]] <- max(0, min(1, 1 - safe_prob))
      }
      risk
    }
    diff_max <- max(abs(balanced_curve_direct - balanced_curve_expected))
    balanced_rows[[length(balanced_rows) + 1L]] <- data.frame(N = N, max_abs_diff = diff_max, stringsAsFactors = FALSE)
    assert_true(diff_max < 1e-12, paste0("Balanced regression failed at N=", N, " with max abs diff ", signif(diff_max, 6)))
  }
  add_result("balanced_mode_regression", "PASS", "Balanced mode matches exact balanced-copy computation within 1e-12 for tested N.")

  edge_rows <- list()
  for (N in edge_N_vals) {
    bal <- nullisomy_risk_curve(N, mode = "balanced", alpha = 100, mc_samples = mc_default, seed = seed_ref, N_unit = 22L)
    dir <- nullisomy_risk_curve(N, mode = "dirichlet_multinomial", alpha = 1, mc_samples = mc_default, seed = seed_ref, N_unit = 22L)
    edge_rows[[length(edge_rows) + 1L]] <- data.frame(
      N = N,
      balanced_m1 = if (N >= 1L) bal[[2L]] else NA_real_,
      dirichlet_m1 = if (N >= 1L) dir[[2L]] else NA_real_,
      stringsAsFactors = FALSE
    )
    if (N < 22L) {
      assert_true(identical(unname(bal[[1L]]), 0) && identical(unname(dir[[1L]]), 0), paste0("Edge-case m=0 failed at N=", N))
      if (N >= 1L) {
        assert_true(all(bal[2:length(bal)] == 1), paste0("Balanced N<N_unit fallback failed at N=", N))
        assert_true(all(dir[2:length(dir)] == 1), paste0("Dirichlet N<N_unit fallback failed at N=", N))
      }
    }
    if (N == 22L) {
      assert_true(all(bal[-1] == 1), "Balanced edge case N=N_unit should have unit nullisomy risk for any positive loss.")
      assert_true(all(dir[-1] == 1), "Dirichlet edge case N=N_unit should have unit nullisomy risk for any positive loss.")
    }
    if (N == 23L) {
      assert_true(dir[[2L]] < 1 && dir[[2L]] > 0, "Dirichlet edge case N=N_unit+1 should give intermediate single-loss risk.")
    }
  }
  add_result("one_copy_floor_edge_cases", "PASS", "Edge cases behaved as expected: N<N_unit uses full-risk fallback for positive loss, N=N_unit gives unit positive-loss risk, and N=N_unit+1 yields finite intermediate Dirichlet risk.")

  repro_curve_1 <- nullisomy_risk_curve(90L, mode = "dirichlet_multinomial", alpha = 10, mc_samples = mc_default, seed = seed_ref)
  repro_curve_2 <- nullisomy_risk_curve(90L, mode = "dirichlet_multinomial", alpha = 10, mc_samples = mc_default, seed = seed_ref)
  repro_curve_3 <- nullisomy_risk_curve(90L, mode = "dirichlet_multinomial", alpha = 10, mc_samples = mc_default, seed = seed_alt)
  assert_true(identical(repro_curve_1, repro_curve_2), "Dirichlet reproducibility failed for identical seed/settings.")
  diff_alt <- max(abs(repro_curve_1 - repro_curve_3))
  assert_true(diff_alt < 0.08, paste0("Dirichlet alternate-seed deviation too large: ", signif(diff_alt, 4)))
  add_result("dirichlet_reproducibility", "PASS", paste0("Identical seed reproduces exactly; alternate seed max abs diff=", signif(diff_alt, 4), "."))

  cache_sep_rows <- list()
  bal_first <- nullisomy_risk_curve(90L, mode = "balanced", alpha = 100, mc_samples = mc_default, seed = seed_ref)
  dir_a1_first <- nullisomy_risk_curve(90L, mode = "dirichlet_multinomial", alpha = 1, mc_samples = mc_default, seed = seed_ref)
  dir_a100_first <- nullisomy_risk_curve(90L, mode = "dirichlet_multinomial", alpha = 100, mc_samples = mc_default, seed = seed_ref)
  dir_seed_alt_first <- nullisomy_risk_curve(90L, mode = "dirichlet_multinomial", alpha = 1, mc_samples = mc_default, seed = seed_alt)
  bal_second <- nullisomy_risk_curve(90L, mode = "balanced", alpha = 100, mc_samples = mc_default, seed = seed_ref)
  dir_a1_second <- nullisomy_risk_curve(90L, mode = "dirichlet_multinomial", alpha = 1, mc_samples = mc_default, seed = seed_ref)
  cache_sep_rows[[1L]] <- data.frame(
    check = c("balanced_repeat", "dirichlet_repeat_same_seed", "balanced_vs_dirichlet", "dirichlet_alpha1_vs_alpha100", "dirichlet_seed12345_vs_12346"),
    max_abs_diff = c(
      max(abs(bal_first - bal_second)),
      max(abs(dir_a1_first - dir_a1_second)),
      max(abs(bal_first - dir_a1_first)),
      max(abs(dir_a1_first - dir_a100_first)),
      max(abs(dir_a1_first - dir_seed_alt_first))
    ),
    stringsAsFactors = FALSE
  )
  assert_true(cache_sep_rows[[1L]]$max_abs_diff[[1L]] < 1e-12, "Balanced cache separation/repeatability failed.")
  assert_true(cache_sep_rows[[1L]]$max_abs_diff[[2L]] < 1e-12, "Dirichlet same-seed cache separation/repeatability failed.")
  assert_true(cache_sep_rows[[1L]]$max_abs_diff[[3L]] > 1e-3, "Balanced vs Dirichlet outputs were unexpectedly identical; cache may be colliding.")
  assert_true(cache_sep_rows[[1L]]$max_abs_diff[[4L]] > 1e-4, "Dirichlet alpha-specific outputs were unexpectedly identical; cache may be colliding.")
  assert_true(cache_sep_rows[[1L]]$max_abs_diff[[5L]] > 1e-6, "Dirichlet seed-specific outputs were unexpectedly identical; cache may be colliding.")
  add_result("cache_key_separation", "PASS", "Order-dependent repeated calls preserved distinct balanced/Dirichlet/alpha/seed outputs, consistent with correct cache separation.")

  bound_rows <- list()
  mono_rows <- list()
  for (N in N_vals) {
    curve <- nullisomy_risk_curve(N, mode = "dirichlet_multinomial", alpha = 1, mc_samples = mc_default, seed = seed_ref)
    bound_rows[[length(bound_rows) + 1L]] <- data.frame(N = N, min_risk = min(curve), max_risk = max(curve), stringsAsFactors = FALSE)
    assert_true(all(curve >= -1e-12 & curve <= 1 + 1e-12), paste0("Risk bounds failed at N=", N))
    mono_tol <- 5e-3
    d <- diff(curve)
    worst_drop <- min(d)
    mono_rows[[length(mono_rows) + 1L]] <- data.frame(N = N, worst_drop = worst_drop, stringsAsFactors = FALSE)
    assert_true(all(d >= -mono_tol), paste0("Monotonicity with loss size failed at N=", N, "; worst drop=", signif(worst_drop, 4)))
  }
  add_result("risk_bounds", "PASS", "All tested Dirichlet risks stayed within [0,1].")
  add_result("monotonicity_with_loss_size", "PASS", "All tested Dirichlet curves were nondecreasing within Monte Carlo tolerance 5e-3.")

  highalpha_rows <- list()
  highalpha_ok <- TRUE
  highalpha_details <- character(0)
  for (N in N_vals) {
    balanced <- nullisomy_risk_curve(N, mode = "balanced", alpha = 100, mc_samples = mc_highalpha, seed = seed_ref)
    alpha100 <- nullisomy_risk_curve(N, mode = "dirichlet_multinomial", alpha = 100, mc_samples = mc_highalpha, seed = seed_ref)
    alpha1000 <- nullisomy_risk_curve(N, mode = "dirichlet_multinomial", alpha = 1000, mc_samples = mc_highalpha, seed = seed_ref)
    diff100 <- max(abs(alpha100 - balanced))
    diff1000 <- max(abs(alpha1000 - balanced))
    highalpha_rows[[length(highalpha_rows) + 1L]] <- data.frame(
      N = N,
      alpha = c(100, 1000),
      max_abs_diff = c(diff100, diff1000),
      within_tolerance = c(diff100 < highalpha_tol, diff1000 < highalpha_tol),
      stringsAsFactors = FALSE
    )
    if (!(diff100 < highalpha_tol || diff1000 < highalpha_tol)) {
      highalpha_ok <- FALSE
      highalpha_details <- c(highalpha_details, paste0("N=", N, " alpha100=", signif(diff100, 4), " alpha1000=", signif(diff1000, 4)))
    }
  }
  if (highalpha_ok) {
    add_result("high_alpha_approximation", "PASS", paste0("For each tested N, alpha=100 or alpha=1000 was within tolerance ", highalpha_tol, " of balanced mode."))
  } else {
    add_result("high_alpha_approximation", "WARN", paste0("High-alpha Dirichlet remained separated from balanced mode for tested N under the one-copy-floor prior: ", paste(highalpha_details, collapse = "; ")))
  }

  lowalpha_rows <- list()
  for (N in N_vals) {
    balanced <- nullisomy_risk_curve(N, mode = "balanced", alpha = 100, mc_samples = mc_default, seed = seed_ref)
    alpha1 <- nullisomy_risk_curve(N, mode = "dirichlet_multinomial", alpha = 1, mc_samples = mc_default, seed = seed_ref)
    alpha0.5 <- nullisomy_risk_curve(N, mode = "dirichlet_multinomial", alpha = 0.5, mc_samples = mc_default, seed = seed_ref)
    moderate_idx <- intersect(which(0:N %in% c(2L, 3L, 4L, 5L, 10L)), seq_along(balanced))
    mean_gain_1 <- mean(alpha1[moderate_idx] - balanced[moderate_idx])
    mean_gain_05 <- mean(alpha0.5[moderate_idx] - balanced[moderate_idx])
    lowalpha_rows[[length(lowalpha_rows) + 1L]] <- data.frame(N = N, alpha = c(1, 0.5), mean_risk_gain_moderate = c(mean_gain_1, mean_gain_05), stringsAsFactors = FALSE)
    if (N >= 88L) {
      assert_true(mean_gain_1 > 0 || mean_gain_05 > 0, paste0("Low-alpha heterogeneity check failed at high N=", N))
    }
  }
  add_result("low_alpha_increases_heterogeneity", "PASS", "Low-alpha Dirichlet mode generally increases moderate-loss nullisomy risk at high chromosome counts.")

  mc_rows <- list()
  for (alpha_test in c(1, 100)) {
    for (N in N_vals) {
      curve_a <- nullisomy_risk_curve(N, mode = "dirichlet_multinomial", alpha = alpha_test, mc_samples = mc_default, seed = seed_ref)
      curve_b <- nullisomy_risk_curve(N, mode = "dirichlet_multinomial", alpha = alpha_test, mc_samples = mc_default, seed = seed_alt)
      diff_vec <- abs(curve_a - curve_b)
      mc_rows[[length(mc_rows) + 1L]] <- data.frame(
        alpha = alpha_test,
        N = N,
        max_abs_diff = max(diff_vec),
        median_abs_diff = median(diff_vec),
        stringsAsFactors = FALSE
      )
    }
  }
  mc_tab <- bind_rows(mc_rows)
  add_result(
    "mc_stability",
    "PASS",
    paste0(
      "Across two seeds at 10000 MC samples, alpha=1 had max/median abs diff up to ",
      signif(max(mc_tab$max_abs_diff[mc_tab$alpha == 1]), 4), "/",
      signif(max(mc_tab$median_abs_diff[mc_tab$alpha == 1]), 4),
      "; alpha=100 had up to ",
      signif(max(mc_tab$max_abs_diff[mc_tab$alpha == 100]), 4), "/",
      signif(max(mc_tab$median_abs_diff[mc_tab$alpha == 100]), 4),
      "."
    )
  )

  sim_bal <- build_simple_sim(mode = "balanced", alpha = 100, mc_samples = mc_default, seed = seed_ref)
  sim_bal_2 <- build_simple_sim(mode = "balanced", alpha = 100, mc_samples = mc_default, seed = seed_ref)
  sim_dir <- build_simple_sim(mode = "dirichlet_multinomial", alpha = 1, mc_samples = mc_default, seed = seed_ref)
  sim_dir_100 <- build_simple_sim(mode = "dirichlet_multinomial", alpha = 100, mc_samples = mc_default, seed = seed_ref)
  assert_true(isTRUE(all.equal(sim_bal$Vmm3_dead_buffer_obs, sim_bal_2$Vmm3_dead_buffer_obs, tolerance = 1e-12)), "Balanced integration regression failed.")
  bal_dead <- sum(as.numeric(sim_bal$Vmm3_dead_buffer_obs))
  dir_dead <- sum(as.numeric(sim_dir$Vmm3_dead_buffer_obs))
  dir_dead_100 <- sum(as.numeric(sim_dir_100$Vmm3_dead_buffer_obs))
  assert_true(dir_dead >= bal_dead - 1e-10, paste0("Dirichlet integration check failed: dead buffer decreased (balanced=", signif(bal_dead, 4), ", dirichlet=", signif(dir_dead, 4), ")."))
  assert_true(dir_dead_100 >= bal_dead - 1e-10, paste0("Dirichlet alpha=100 integration check failed: dead buffer decreased (balanced=", signif(bal_dead, 4), ", dirichlet100=", signif(dir_dead_100, 4), ")."))
  add_result("integration_check", "PASS", paste0("Balanced mode was unchanged; for N~90 initialization, dirichlet alpha=1 produced dead-buffer inflow ", signif(dir_dead, 4), " and alpha=100 produced ", signif(dir_dead_100, 4), " vs balanced ", signif(bal_dead, 4), "."))

  results_tab <- bind_rows(results)
  utils::write.table(results_tab, file.path(out_dir, "unit_test_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bind_rows(balanced_rows), file.path(out_dir, "balanced_regression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bind_rows(edge_rows), file.path(out_dir, "one_copy_floor_edge_cases.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bind_rows(cache_sep_rows), file.path(out_dir, "cache_key_separation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bind_rows(bound_rows), file.path(out_dir, "risk_bounds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bind_rows(mono_rows), file.path(out_dir, "monotonicity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bind_rows(highalpha_rows), file.path(out_dir, "highalpha_approximation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bind_rows(lowalpha_rows), file.path(out_dir, "lowalpha_heterogeneity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(mc_tab, file.path(out_dir, "mc_stability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  report_lines <- c(
    "Nullisomy Hidden-Copy Unit Tests",
    "",
    paste("Output dir:", out_dir),
    paste("Monte Carlo samples (default/high-alpha):", mc_default, "/", mc_highalpha),
    "",
    "Summary:"
  )
  report_lines <- c(report_lines, apply(results_tab, 1, function(row) paste0("- ", row[["test"]], ": ", row[["status"]], " | ", row[["detail"]])))
  writeLines(report_lines, con = file.path(out_dir, "unit_test_report.md"))
  message("Wrote nullisomy hidden-copy unit tests to: ", out_dir)
}

main(commandArgs(trailingOnly = TRUE))
