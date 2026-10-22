#!/usr/bin/env Rscript

# Canonical fitted-configuration and fixed-O2 model-probe helpers.

o2pr_first_seed_cfg <- function(manifest) {
    for (p in manifest$config_file) {
        if (!is.na(p) && file.exists(p)) {
            cfg <- tryCatch(readRDS(p), error = function(e) NULL)
            if (!is.null(cfg))
                return(cfg)
        }
    }
    list(N_MIN = 22L, N_MAX = 154L, N_UNIT = 22L, DT = 0.05, start_with = "chr_number")
}


o2pr_cfg_metadata <- function(cfg) {
    data.frame(metric = c("Nmin", "Nmax", "N_UNIT", "DT", "start_with", "boundary_default", "o2_burden_feedback", "O2_growth",
        "ploidy_O2_death", "o2_S0_upper_bound", "o2_Nref", "trajectory_value_semantics"), value = c(cfg$N_MIN %||% NA, cfg$N_MAX %||%
        NA, cfg$N_UNIT %||% NA, cfg$DT %||% NA, cfg$start_with %||% NA, cfg$boundary %||% "drop", cfg$o2_burden_feedback %||%
        NA, cfg$O2_growth %||% NA, cfg$ploidy_O2_death %||% NA, cfg$o2_S0_upper_bound %||% NA, cfg$o2_Nref %||% NA, if (identical(as.character(cfg$start_with %||%
        "ploidy"), "chr_number")) {
        "weighted mean chromosome number N"
    } else {
        "weighted mean ploidy converted to N when needed"
    }), stringsAsFactors = FALSE)
}


o2pr_build_G <- function(model_env, cfg, run_params, O2) {
    fn <- get("cpp_o2simps_build_G_for_o2_triplet", envir = model_env, inherits = TRUE)
    tri <- fn(O2 = as.numeric(O2), O2_crit = as.numeric(run_params$O2_crit %||% cfg$o2_crit_init %||% 1), N0min = as.integer(cfg$N_MIN %||%
        22L), N0max = as.integer(cfg$N_MAX %||% 154L), N1min = as.integer(cfg$N_MIN %||% 22L), N1max = as.integer(cfg$N_MAX %||%
        154L), lam_max = as.numeric(run_params$lam_max), p_mis_base = as.numeric(run_params$p_mis_base %||% cfg$p_mis_base %||%
        1e-05), p_misseg = as.numeric(run_params$p_misseg %||% 0), k_o_mis = as.numeric(run_params$k_o_mis %||% 50), p_wgd = as.numeric(run_params$p_wgd %||%
        0), boundary = as.character(run_params$boundary %||% "drop"), eps_tail = 1e-08, buffer_smax = as.numeric(run_params$buffer_smax %||%
        1), buffer_beta = as.numeric(run_params$buffer_beta %||% 0), buffer_n_exp = as.numeric(run_params$buffer_n_exp %||%
        1), N_unit = as.integer(cfg$N_UNIT %||% 22L), beta_size = 0, O2_growth = isTRUE(run_params$O2_growth %||% cfg$O2_growth %||%
        TRUE), alpha_o2 = as.numeric(run_params$alpha_o2 %||% 0), gamma_growth = as.numeric(run_params$gamma_growth %||%
        1), mu_hp = as.numeric(run_params$mu_hp %||% 0), gamma_mu = as.numeric(run_params$gamma_mu %||% 1), n_O = as.numeric(run_params$n_O %||%
        1), ploidy_O2_death = as.character(run_params$ploidy_O2_death %||% cfg$ploidy_O2_death %||% "diploid_NULL"))
    G <- Matrix::sparseMatrix(i = as.integer(tri$i), j = as.integer(tri$j), x = as.numeric(tri$x), dims = c(as.integer(tri$nrow),
        as.integer(tri$ncol)), repr = "C")
    attr(G, "triplet") <- tri
    G
}


o2pr_run_params_from_vec <- function(vec, cfg) {
    rp <- as.list(vec)
    rp$o2_min <- if ("o2_min" %in% names(vec))
        vec[["o2_min"]]
    else (cfg$o2_min %||% 0)
    rp$o2_S0_upper_bound <- cfg$o2_S0_upper_bound %||% 5
    rp$o2_Nref <- cfg$o2_Nref %||% 1e+06
    rp$O2_growth <- cfg$O2_growth %||% TRUE
    rp$ploidy_O2_death <- cfg$ploidy_O2_death %||% "ploidy_related"
    rp$boundary <- cfg$boundary %||% "drop"
    rp
}
