#!/usr/bin/env Rscript

# Canonical mapping between natural-scale in-vitro parameters and optimizer
# coordinates. Fit, simulation, and post-fit analysis must share this contract.

o2sd_invitro_parameter_transform_map <- function() {
  data.frame(
    param_symbol = c(
      "lam_max", "p_misseg", "k_o_mis",
      "buffer_smax", "buffer_beta", "buffer_n_exp",
      "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
      "O2_crit", "n_O", "p_mis_base", "sigma_growth", "sigma_kary",
      "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N"
    ),
    param_name = c(
      "log10_lam_max", "log10_p_misseg", "log10_k_o_mis",
      "buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp",
      "log10_p_wgd", "log10_alpha_o2", "gamma_growth", "log10_mu_hp", "gamma_mu",
      "log10_O2_crit", "n_O", "log10_p_mis_base", "log10_sigma_growth", "log10_sigma_kary",
      "init_mean_2N", "log10_init_sd_2N", "init_mean_4N", "log10_init_sd_4N"
    ),
    transform = c(
      "log10", "log10", "log10",
      "identity", "log10", "log10",
      "log10", "log10", "identity", "log10", "identity",
      "log10", "identity", "log10", "log10", "log10",
      "identity", "log10", "identity", "log10"
    ),
    stringsAsFactors = FALSE
  )
}
