`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/.rcpp_cache_richard/sourceCpp-aarch64-apple-darwin20-1.1.0/sourcecpp_75c552490c4d/sourceCpp_2.so')

cpp_richard_pr_delta_vec <- Rcpp:::sourceCppFunction(function(N, p, eps_tail = 1e-8, beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_richard_pr_delta_vec')
cpp_richard_build_B_total_triplet <- Rcpp:::sourceCppFunction(function(Nmin, Nmax, p_vec, boundary = "drop", eps_tail = 1e-8, beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_richard_build_B_total_triplet')
cpp_richard_build_B_WGD_triplet <- Rcpp:::sourceCppFunction(function(N0min, N0max, N1min, N1max, boundary = "drop", wgd_value = 1.0) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_cpp_richard_build_B_WGD_triplet')

rm(`.sourceCpp_1_DLLInfo`)
