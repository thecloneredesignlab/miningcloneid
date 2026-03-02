`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/O2_mono/.rcpp_cache_richard/shared/sourceCpp-aarch64-apple-darwin20-1.1.0/sourcecpp_876616ffdf53/sourceCpp_2.so')

cpp_richard_pr_delta_vec <- Rcpp:::sourceCppFunction(function(N, p, eps_tail = 1e-8, beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_richard_pr_delta_vec')
cpp_richard_build_B_total_triplet <- Rcpp:::sourceCppFunction(function(Nmin, Nmax, p_vec, boundary = "drop", eps_tail = 1e-8, beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_richard_build_B_total_triplet')
cpp_richard_build_B_WGD_triplet <- Rcpp:::sourceCppFunction(function(N0min, N0max, N1min, N1max, boundary = "drop", wgd_value = 1.0) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_richard_build_B_WGD_triplet')
cpp_richard_simulate_one <- Rcpp:::sourceCppFunction(function(init_state, N0min, N0max, N1min, N1max, obs_steps, sim_end_step, DT, dose, dose_ref, treat_day, fit_treatment, alpha, gamma, tx_mult_min, crowding, K, min_pop, O2_base, o2_feedback, o2_min, h_O2, K_O2, lam_min, lam_max, k_o, p_misseg, k_o_mis, beta_buffer, n_exp, smax, p_wgd, N_unit, vol_by_N, burden_floor) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_richard_simulate_one')

rm(`.sourceCpp_1_DLLInfo`)
