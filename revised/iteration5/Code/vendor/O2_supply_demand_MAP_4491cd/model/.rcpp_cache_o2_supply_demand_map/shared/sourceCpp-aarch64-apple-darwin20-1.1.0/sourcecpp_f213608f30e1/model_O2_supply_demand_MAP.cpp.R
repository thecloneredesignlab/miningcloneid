`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map/shared/sourceCpp-aarch64-apple-darwin20-1.1.0/sourcecpp_f213608f30e1/sourceCpp_2.so')

cpp_o2simps_pr_delta_vec <- Rcpp:::sourceCppFunction(function(N, p, eps_tail = 1e-8, buffer_smax = 1.0, buffer_beta = 0.0, buffer_n_exp = 1.0, N_unit = 22L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_o2simps_pr_delta_vec')
cpp_o2simps_o2_window_supply <- Rcpp:::sourceCppFunction(function(Ntot, o2_S0 = 0.5, kappa_O = 1.0, o2_Nref = 1e6, o2_min = 0.0, o2_S0_upper_bound = 5.0) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_o2simps_o2_window_supply')
cpp_o2simps_build_B_total_triplet <- Rcpp:::sourceCppFunction(function(Nmin, Nmax, p_vec, boundary = "drop", eps_tail = 1e-8, buffer_smax = 1.0, buffer_beta = 0.0, buffer_n_exp = 1.0, N_unit = 22L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_o2simps_build_B_total_triplet')
cpp_o2simps_build_B_WGD_triplet <- Rcpp:::sourceCppFunction(function(N0min, N0max, N1min, N1max, boundary = "drop", wgd_value = 1.0) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_o2simps_build_B_WGD_triplet')
cpp_o2simps_build_G_for_o2_triplet <- Rcpp:::sourceCppFunction(function(O2, O2_crit, N0min, N0max, N1min, N1max, lam_max, p_mis_base, p_misseg, k_o_mis, p_wgd = 0.0, boundary = "drop", eps_tail = 1e-8, buffer_smax = 1.0, buffer_beta = 0.0, buffer_n_exp = 1.0, N_unit = 22L, beta_size = 0.0, O2_growth = TRUE, alpha_o2 = 0.0, gamma_growth = 1.0, mu_hp = 0.0, gamma_mu = 1.0, n_O = 1.0, ploidy_O2_death = "diploid_NULL") {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_o2simps_build_G_for_o2_triplet')
cpp_o2simps_simulate_one <- Rcpp:::sourceCppFunction(function(sim_args) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_o2simps_simulate_one')
cpp_o2simps_objective_components_map <- Rcpp:::sourceCppFunction(function(scenario_data, objective_data, state_data, sim_args) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_o2simps_objective_components_map')

rm(`.sourceCpp_1_DLLInfo`)
