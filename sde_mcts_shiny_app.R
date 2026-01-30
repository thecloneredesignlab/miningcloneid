library(shiny)
library(deSolve)

eps <- 1e-6

drug_dosing_schedules <- list(
    "volasertib"    = "IV",
    "umi-77"        = "IV",
    "tegafur"       = "Oral",
    "tas"           = "Oral",
    "osi-027"       = "Oral",
    "alisertib"     = "IV",
    "5-azacytidine" = "Oral",
    "abt-199"       = "Oral",
    "abt-263"       = "Oral",
    "capecitabine"  = "Oral",
    "ceralasertib"  = "Oral",
    "cytarabine"    = "IV",
    "gemcitabine"   = "IV",
    "bay1895344"    = "IV",
    "ispinesib"     = "IV",
    "navitoclax"    = "Oral",
    "adavosertib"   = "Oral",
    "none"          = "None"
)

IV_DEFAULTS <- list(C_peak=1.0, half_life=0.5, period=7.0)

ORAL_DEFAULTS <- list(
    dose = 100.0, F = 0.6, Vd = 60.0,
    ka_day = 2.0, ke_day = 0.7,
    period = 1.0, tlag = 0.0
)

unknown_value <- 1.0
PER_DRUG <- list(
    # IV drugs
    "volasertib"   = list(C_peak = 1.0, half_life = 4.0, period = 7.0),
    "alisertib"    = list(C_peak = 42.5, half_life = 0.875, period = 0.5),
    "cytarabine"   = list(C_peak = 1.0, half_life = 0.2, period = 3.5),
    "gemcitabine"  = list(C_peak = 239, half_life = 0.05, period = 7.0),
    "ispinesib"    = list(C_peak = 2.1, half_life = 1.04, period = 7),
    "umi-77"       = list(C_peak = 1.0, half_life = 0.8, period = 7.0),
    "navitoclax"   = list(C_peak = 1.0, half_life = 0.73, period = 1),
    "bay1895344"   = list(C_peak = 6.2, half_life = 0.50, period = 0.5),
    "none"         = list(C_peak = 0, half_life = 0.50, period = 7),
    "Topotecan"    = list(C_peak = unknown_value, half_life = unknown_value, period = unknown_value),
    "Doxorubicin"  = list(C_peak = unknown_value, half_life = unknown_value, period = unknown_value),

    # Oral drugs
    "osi-027"      = list(dose = 50, F = 0.5, Vd = 80, ka_day = 1.8, ke_day = 0.5, period = 1.0),
    "abt-199"      = list(dose = 100, F = 0.6, Vd = 250, ka_day = 1.0, ke_day = 0.3, period = 1.0),
    "abt-263"      = list(dose = 100, F = 0.6, Vd = 120, ka_day = 2.0, ke_day = 0.5, period = 1.0),
    "ceralasertib" = list(dose = 80,  F = 0.5, Vd = 100,ka_day = 2.0, ke_day = 0.5, period = 1.0),
    "adavosertib"  = list(dose = 100, F = 0.6, Vd = 65, ka_day = 2.4, ke_day = 0.6, period = 1.0),
    "tas"          = list(dose = 60,  F = 0.5, Vd = 40, ka_day = 2.4, ke_day = 0.7, period = 1.0),
    "tegafur"      = list(dose = 40, F = 0.5,Vd = 45, ka_day = 1.6, ke_day = 0.5, period = 1.0),
    "capecitabine" = list(dose = 100, F = 0.8, Vd = 40, ka_day = 3.0, ke_day = 0.6, period = 1.0),
    "5-azacytidine"= list(dose = 100, F = 0.2, Vd = 40, ka_day = 3.0, ke_day = 2.0, period = 1.0),
    "Elimusertib"  = list(dose = 40, F = unknown_value, Vd = unknown_value, ka_day = unknown_value, ke_day = unknown_value, period = unknown_value)
)

# Drug cycle lengths (in days)
# Default is 14.0 days
default <- 14.0
cycle_lengths <- list(
    # IV drugs
    "volasertib"   = default,
    "alisertib"    = 21.0,
    "cytarabine"   = default,
    "gemcitabine"  = 28.0,
    "ispinesib"    = default,
    "umi-77"       = default,
    "navitoclax"   = default,
    "bay1895344"   = default,
    "Topotecan"    = 28.0,
    "Doxorubicin"  = 21.0,

    # Oral drugs
    "osi-027"      = default,
    "abt-199"      = default,
    "abt-263"      = default,
    "ceralasertib" = default,
    "adavosertib"  = default,
    "tas"          = default,
    "tegafur"      = default,
    "capecitabine" = default,
    "5-azacytidine"= default,
    "Elimusertib"  = 7.0,

    "none"         = 1.0
)


# --- PD Function: phi_Hill ---
phi_Hill <- function(C, EC50, n, Emax=1.0) {
    # Hill-type kill rate function
    Emax * (C^n) / (EC50^n + C^n)
}

# --- PD Function: f(ploidy, drug) ---
f_pd_params <- function(ploidy, drug) {
    drug <- tolower(drug)

    # Helper functions from Python
    clamp_ec50 <- function(x) { max(x, 1e-12) }
    clamp_n <- function(x) { max(x, 0.1) }

    if (drug == "bay1895344") {
        n_out <- clamp_n(3.85 * exp(-0.861 * ploidy) + 0.81)
        ec50  <- clamp_ec50(1.04 * exp(0.35 * ploidy) - 2.05)
        list(EC50=ec50, n=n_out, Emax=1.0)

    } else if (drug == "alisertib") {
        n_out <- clamp_n(1.0)
        ec50  <- clamp_ec50(51.02 * exp(-0.62 * ploidy) - 4.78)
        list(EC50=ec50, n=n_out, Emax=1.0)

    } else if (drug == "ispinesib") {
        n_out <- clamp_n(0.94 * exp(-0.303 * ploidy) - 0.73)
        ec50  <- clamp_ec50(1.185 * exp(-0.21 * ploidy) - 0.56)
        list(EC50=ec50, n=n_out, Emax=1.0)

    } else if (drug == "gemcitabine") {
        n_out <- clamp_n(28.92 * exp(-0.94 * ploidy) + 0.92)
        ec50  <- clamp_ec50(0.004 * exp(0.78 * ploidy) - 0.01) # This formula is updated
        list(EC50=ec50, n=n_out, Emax=1.0)

    } else if (drug == "none") {   # <--- ADD THIS BLOCK
        list(EC50=1.0, n=1.0, Emax=0.0)

    } else {
        # Default fallback for unknown drugs
        warning(paste("Unknown drug:", drug, "- using default parameters (EC50=1.0, n=1.0)."))
        list(EC50=1.0, n=1.0, Emax=1.0)
    }
}

# --- PK Function: pulsed_dose ---
pulsed_dose <- function(C_peak=5.0, half_life=2.0, period=7.0) {
    lam <- log(2) / max(half_life, 1e-12)
    function(t) {
        t <- as.numeric(t)
        modt <- t %% period
        C_peak * exp(-lam * modt)
    }
}

# --- PK Function: oral_pulsed_ss_days ---
oral_pulsed_ss_days <- function(dose=100.0, F=0.7, Vd=70.0, ka_day=1.2, ke_day=0.3, period=1.0, tlag=0.0) {
    if (abs(ka_day - ke_day) < 1e-12) {
        # safe limiting form if ka ≈ ke
        return(function(t) {
            t <- as.numeric(t)
            tau <- as.numeric(period)
            tstar <- (t - tlag) %% tau
            num <- exp(-ke_day * tstar)
            den <- max(1.0 - exp(-ke_day * tau), 1e-12)
            (F * dose / Vd) * (ke_day * tstar) * num / den
        })
    }

    A <- (F * dose * ka_day) / (Vd * (ka_day - ke_day))
    function(t) {
        t <- as.numeric(t)
        tau <- as.numeric(period)
        tstar <- (t - tlag) %% tau
        term_elim <- exp(-ke_day * tstar) / max(1.0 - exp(-ke_day * tau), 1e-12)
        term_abs  <- exp(-ka_day * tstar) / max(1.0 - exp(-ka_day * tau), 1e-12)
        A * (term_elim - term_abs)
    }
}

# --- PK Function: get_concentration_curve ---
get_concentration_curve <- function(drug_name, ...) {
    drug_name <- tolower(drug_name)
    route <- drug_dosing_schedules[[drug_name]]

    if (is.null(route)) {
        warning(paste("Unknown drug:", drug_name, "- using IV defaults."))
        route <- "IV"
    }

    # Gather per-drug params (if any), then apply overrides
    per_drug <- PER_DRUG[[drug_name]]
    if (is.null(per_drug)) per_drug <- list()

    overrides <- list(...)

    if (route == "IV") {
        base_params <- IV_DEFAULTS
        params <- c(per_drug, overrides) # Overrides take precedence
        final_params <- modifyList(base_params, params)
        return(do.call(pulsed_dose, final_params))

    } else if (route == "Oral") {
        base_params <- ORAL_DEFAULTS
        params <- c(per_drug, overrides) # Overrides take precedence
        final_params <- modifyList(base_params, params)
        return(do.call(oral_pulsed_ss_days, final_params))

    } else if (route == "None") {
        return(function(t) 0.0)

    } else {
        stop(paste("Unsupported route", route))
    }
}

# Simulate SDE using Exponential Euler
simulate_sde <- function(ploidy_status, drug, T_span, dt=0.1, r=0.4, K=6e10, n_sims=1000, beta_by_ploidy=NULL) {
    ploidies <- sort(as.numeric(names(ploidy_status)))
    T0 <- as.numeric(ploidy_status[as.character(ploidies)])
    M <- length(ploidies)

    # Get PD parameters for each ploidy
    phi_params <- lapply(ploidies, function(p) f_pd_params(p, drug))

    # Growth rates
    r_vec <- rep(r, M)

    # Beta (volatility) for each ploidy
    if (is.null(beta_by_ploidy)) {
        beta_vec <- rep(0.02, M)
    } else {
        beta_vec <- sapply(ploidies, function(p) beta_by_ploidy[[as.character(p)]])
    }

    # Get concentration function
    C_func <- get_concentration_curve(drug)

    # Time grid
    times <- seq(T_span[1], T_span[2], by=dt)
    N <- length(times)

    # Initialize paths: (n_sims, M, N)
    Tpaths <- array(0, dim=c(n_sims, M, N))
    Tpaths[, , 1] <- matrix(rep(T0, each=n_sims), nrow=n_sims, ncol=M)

    # Exponential Euler update
    for (k in 1:(N-1)) {
        t <- times[k]
        dtk <- times[k+1] - t
        sqrt_dtk <- sqrt(max(dtk, eps))

        C <- C_func(t)
        phi_vals <- sapply(phi_params, function(pp) phi_Hill(C, pp$EC50, pp$n, pp$Emax))

        T_curr <- Tpaths[, , k]  # (n_sims, M)
        Tsum <- rowSums(T_curr)  # (n_sims)

        # Drift term: mu = r*(1 - Tsum/K) - phi
        mu <- sweep(matrix(r_vec, nrow=n_sims, ncol=M, byrow=TRUE), 1, Tsum/(K + 1e-6), function(x,y) x*(1-y)) -
              matrix(phi_vals, nrow=n_sims, ncol=M, byrow=TRUE)

        # Diffusion term
        Z <- matrix(rnorm(n_sims * M), nrow=n_sims, ncol=M)
        dW <- matrix(beta_vec, nrow=n_sims, ncol=M, byrow=TRUE) * sqrt_dtk * Z
        beta_eff_sq <- matrix(beta_vec^2, nrow=n_sims, ncol=M, byrow=TRUE)

        # Exponential Euler step
        expo <- (mu - 0.5 * beta_eff_sq) * dtk + dW
        T_next <- T_curr * exp(expo)

        Tpaths[, , k+1] <- T_next
    }

    return(list(times=times, ploidies=ploidies, Tpaths=Tpaths))
}

# Main forecasting function (simplified R version of ploidy_forcast)
ploidy_forcast <- function(ploidy_cell_count, drug, T, R_BASE=0.4, K_CAP=6e10,
                           beta_by_ploidy=NULL, DT=0.1, N_SIMS=1000) {

    result <- simulate_sde(ploidy_cell_count, drug, T_span=c(0, T), dt=DT,
                          r=R_BASE, K=K_CAP, n_sims=N_SIMS, beta_by_ploidy=beta_by_ploidy)

    # Return ploidies, t_ode (approx as t_sde), T_mat_ode (mean trajectory), t_sde, Tpaths
    mean_trajectory <- apply(result$Tpaths, c(2, 3), mean)  # (M, N)

    return(list(
        ploidies = result$ploidies,
        t_ode = result$times,
        T_mat_ode = mean_trajectory,
        t_sde = result$times,
        Tpaths = result$Tpaths
    ))
}

# Node class for MCTS
Node <- setRefClass("Node",
  fields = list(
    ploidy_status = "list",
    cycle = "numeric",
    parent = "ANY",
    children = "list",
    N = "numeric",
    W = "numeric"
  ),
  methods = list(
    initialize = function(ploidy_status_val, cycle_val, parent_val=NULL) {
      ploidy_status <<- ploidy_status_val
      cycle <<- cycle_val
      parent <<- parent_val
      children <<- list()
      N <<- 0
      W <<- 0.0
    },
    is_terminal = function(total_cycles, min_size, max_size) {
      total <- sum(unlist(ploidy_status))
      return(cycle >= total_cycles || total < min_size || total > max_size)
    },
    is_fully_expanded = function(drugs) {
      return(length(children) == length(drugs))
    }
  )
)

# Simulate next state with uncertainty calculation (from sde_mcts_shiny_app.R)
simulate_next_state_mcts <- function(ploidy_status, drug, d_switch, N_SIMS=1000) {
    result <- ploidy_forcast(ploidy_status, drug, T=d_switch, N_SIMS=N_SIMS)

    Tpaths <- result$Tpaths
    ploidies <- result$ploidies

    # Calculate final total burdens across all simulations
    final_total_burdens <- rowSums(Tpaths[, , dim(Tpaths)[3]])

    # Freedman-Diaconis Rule for Binning
    n <- length(final_total_burdens)
    if (n > 1) {
        q75 <- quantile(final_total_burdens, 0.75)
        q25 <- quantile(final_total_burdens, 0.25)
        iqr <- q75 - q25

        if (iqr > eps) {
            bin_width <- 2 * iqr * (n ^ (-1/3))
            data_range <- max(final_total_burdens) - min(final_total_burdens)
            # Prevent integer overflow and ensure reasonable number of bins
            num_bins <- min(as.integer(ceiling(data_range / (bin_width + eps))), 100)
            num_bins <- max(num_bins, 1)
        } else {
            num_bins <- 1
        }
    } else {
        num_bins <- 1
    }

    # Calculate Shannon Entropy with safe histogram
    uncertainty_score <- 0.0
    tryCatch({
        # Create explicit breaks to avoid integer overflow
        min_val <- min(final_total_burdens)
        max_val <- max(final_total_burdens)

        # Use explicit sequence of breaks instead of just specifying number
        if (max_val > min_val) {
            breaks <- seq(min_val, max_val, length.out = min(num_bins + 1, 101))
            hist_result <- hist(final_total_burdens, breaks=breaks, plot=FALSE)
            counts <- hist_result$counts

            # Filter out zero counts
            counts <- counts[counts > 0]

            if (length(counts) > 0) {
                probs <- counts / (sum(counts) + eps)
                probs <- probs[probs > eps]

                if (length(probs) > 0) {
                    entropy <- -sum(probs * log(probs + eps))
                    max_entropy <- log(length(counts) + eps)
                    if (max_entropy > 1e-9) {
                        uncertainty_score <- entropy / (max_entropy + eps)
                    }
                }
            }
        }
    }, error = function(e) {
        # If histogram fails, set uncertainty to 0
        uncertainty_score <- 0.0
    })
    uncertainty_score <- min(max(uncertainty_score, 0.0), 1.0)  # Clamp to [0, 1]

    # Calculate mean final ploidy distribution
    final_per_ploidy <- Tpaths[, , dim(Tpaths)[3]]  # (n_sims, num_ploidies)
    mean_sde_per_ploidy <- colMeans(final_per_ploidy)

    # Mean trajectory over time
    mean_trajectory <- apply(Tpaths, c(2, 3), mean)  # (num_ploidies, time_points)
    y <- t(mean_trajectory[, -1])  # Transpose and skip first time point

    # New status
    new_status <- as.list(mean_sde_per_ploidy)
    names(new_status) <- as.character(ploidies)

    confidence <- 1 - uncertainty_score

    return(list(new_status=new_status, y=y, confidence=confidence))
}

# Select best child using UCB1
select_best_child_to_explore <- function(node, c) {
    best_score <- -Inf
    best_child <- NULL

    for (drug in names(node$children)) {
        child <- node$children[[drug]]
        Q <- child$W / (child$N + 1e-6)
        U <- c * sqrt(log(node$N + 1e-6) / (child$N + 1e-6))
        score <- Q + U

        if (score > best_score) {
            best_score <- score
            best_child <- child
        }
    }

    return(best_child)
}

# Expand node by trying an untried drug
expand_mcts <- function(node, drugs, d_switch) {
    untried <- setdiff(drugs, names(node$children))
    drug <- sample(untried, 1)

    result <- simulate_next_state_mcts(node$ploidy_status, drug, d_switch)
    new_ploidy <- result$new_status

    child <- Node$new(new_ploidy, node$cycle + 1, node)
    node$children[[drug]] <- child

    return(child)
}

# Rollout simulation from a node
rollout_mcts <- function(node, rollout_depth, drugs, d_switch, min_size, max_size,
                        alpha=0.01, p_order=3, safe_size=1e10, beta=1.0) {
    ploidy <- node$ploidy_status
    confidence <- 1.0

    # Initialize path_burdens with the starting burden
    path_burdens <- c(sum(unlist(ploidy)))
    confidences <- c(confidence)

    for (step in 1:rollout_depth) {
        total <- sum(unlist(ploidy))
        if (total < min_size || total > max_size) {
            break
        }

        drug <- sample(drugs, 1)
        result <- simulate_next_state_mcts(ploidy, drug, d_switch)
        ploidy <- result$new_status
        confidence <- result$confidence

        # Collect burden along the path
        path_burdens <- c(path_burdens, sum(unlist(ploidy)))
        confidences <- c(confidences, confidence)
    }

    # Calculate reward using ODE formula
    reward <- 0
    for (burden in path_burdens) {
        rollout_confidence = confidences[which(path_burdens == burden)[1]]
        reward <- reward - ( rollout_confidence * ( (burden / (max_size + eps)) - alpha * ((max(0, burden - safe_size) / (max_size - safe_size)) ^ p_order) ) )
    }

    reward <- reward / length(path_burdens)

    # if extinction is achieved, give a bonus
    if (tail(path_burdens, 1) < min_size) {
        t_idx <- which(path_burdens < min_size)
        steps_to_extinct <- t_idx[1] - 1
        bonus <- beta * ((rollout_depth - steps_to_extinct) / max(1, rollout_depth))
        reward <- reward + bonus
    }

    return(reward)
}

# Backpropagate reward up the tree
backpropagate <- function(node, reward) {
    current <- node
    while (!is.null(current)) {
        current$N <- current$N + 1
        current$W <- current$W + reward
        current <- current$parent
    }
}

# Main MCTS function (from sde_mcts_shiny_app.R)
run_mcts_new <- function(ploidy_status, cycle, drugs, d_switch, total_cycles, min_size, max_size,
                     depth=30, num_rollouts=100, c=sqrt(2),
                     alpha=0.01, p_order=3, safe_size=1e10, beta=1.0) {

    # Initialize root node for current state
    root <- Node$new(ploidy_status, cycle)

    # Run MCTS rollouts
    for (i in 1:num_rollouts) {
        node <- root

        # Selection
        while (node$is_fully_expanded(drugs) && !node$is_terminal(total_cycles, min_size, max_size)) {
            node <- select_best_child_to_explore(node, c)
        }

        # Expansion
        if (!node$is_terminal(total_cycles, min_size, max_size)) {
            child <- expand_mcts(node, drugs, d_switch)
            reward <- rollout_mcts(child, depth, drugs, d_switch, min_size, max_size,
                                 alpha, p_order, safe_size, beta)
            backpropagate(child, reward)
        } else {
            reward <- rollout_mcts(node, 0, drugs, d_switch, min_size, max_size,
                                 alpha, p_order, safe_size, beta)
            backpropagate(node, reward)
        }
    }

    # Pick best drug based on Q-value
    get_q <- function(child) {
        child$W / (child$N + eps)
    }

    best_drug <- NULL
    best_q <- -Inf
    for (drug in names(root$children)) {
        q <- get_q(root$children[[drug]])
        if (q > best_q) {
            best_q <- q
            best_drug <- drug
        }
    }

    return(best_drug)
}

# ==============================================================================
# ====================== Core Shiny App Code (ORIGINAL) ========================
# ==============================================================================

# ---- ODE model ----
model_ode_fn <- function(t, state, parms) {
    # state is a named vector c(B1 = ..., B2 = ...)
    B_total <- sum(state)

    # Get parameters
    r_vec <- parms[["r_vec"]]
    K <- parms[["K"]]
    C_func <- parms[["C_func"]]
    phi_params_list <- parms[["phi_params_list"]]
    k_multiplier <- parms[["k_multiplier"]]

    # 1. Calculate competition
    competition_term <- max(0, (1 - B_total / K))

    # 2. Calculate dynamic kill rate based on C(t)
    C <- C_func(t)

    # Calculate phi_val for each ploidy
    phi_vals <- sapply(phi_params_list, function(p) {
        phi_Hill(C, EC50 = p$EC50, n = p$n, Emax = p$Emax)
    })

    # Apply the fitted multiplier (if any)
    kill_rates <- phi_vals * k_multiplier

    # 3. Calculate net growth
    net_growth_rates <- (r_vec * competition_term) - kill_rates

    # dBi = (ri * (1-B_tot/K) - phi_i(C(t))) * Bi
    dB_vec <- net_growth_rates * state

    list(dB_vec)
}

# --- run_one_cycle ----
run_one_cycle <- function(ploidy_fracs, B0, drug_name = "A",
                          K = 1e9, days = 28, dt = 0.1, k_multiplier_base = 1.0) {

    n <- length(ploidy_fracs)

    # Growth rates (still hard-coded, but separate from drug)
    base_r  <- c(0.020, 0.015, 0.010, 0.008, 0.006)
    r_vec <- rep_len(base_r, n)

    # Using a crude ploidy proxy (e.g., 2, 3, 4... for n=3)
    ploidy_proxy <- seq(2, length.out = n)

    phi_params_list <- lapply(ploidy_proxy, f_pd_params, drug = drug_name)

    # 2. Get PK function C(t)
    C_func <- get_concentration_curve(drug_name)


    initial_state <- ploidy_fracs * B0
    names(initial_state) <- paste0("B", 1:n)

    parms <- list(r_vec = r_vec,
                  K = K,
                  phi_params_list = phi_params_list,
                  C_func = C_func,
                  k_multiplier = k_multiplier_base)

    times <- seq(0, days, by = dt)
    out <- ode(y = initial_state, times = times, func = model_ode_fn, parms = parms)
    out <- as.data.frame(out)

    b_cols <- grep("^B[0-9]+$", names(out), value = TRUE)
    out$B_total <- rowSums(out[, b_cols, drop = FALSE])

    # Store parameters for fitting and summary
    attr(out, "r_vec") <- r_vec
    attr(out, "K")     <- K
    attr(out, "drug")  <- drug_name
    attr(out, "phi_params_list") <- phi_params_list
    attr(out, "C_func") <- C_func
    attr(out, "k_multiplier_base") <- k_multiplier_base

    out
}

# ---- Parsing helpers ----
parse_fractions <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(NULL)
    clean <- gsub("[A-Za-z_=]", " ", txt)
    parts <- unlist(strsplit(clean, "[,\\s]+"))
    parts <- parts[nzchar(parts)]
    vals <- suppressWarnings(as.numeric(parts))
    vals <- vals[!is.na(vals)]
    if (!length(vals)) return(NULL)
    s <- sum(vals)
    if (s <= 0) return(NULL)
    vals / s
}

parse_measurements <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(NULL)
    txt <- gsub("[;]", "\n", txt)
    lines <- unlist(strsplit(txt, "[\n]+"))
    meas <- lapply(lines, function(line) {
        parts <- unlist(strsplit(line, "[,\\s]+"))
        parts <- parts[nzchar(parts)]
        if (length(parts) >= 2) {
            val <- suppressWarnings(as.numeric(parts[1:2]))
            if (all(!is.na(val))) return(data.frame(percent = val[1], time = val[2]))
        }
        NULL
    })
    res <- do.call(rbind, meas)
    if (is.null(res)) data.frame(percent=numeric(0), time=numeric(0)) else res
}

# ==============================================================================
# ====================== UI (ORIGINAL) =========================================
# ==============================================================================

ui <- fluidPage(
    tags$head(tags$style(HTML("
    .container-fluid { max-width: 850px; }
    .small-note { color:#666; font-size: 12px; }
  "))),
    titlePanel("Tumor Burden Across Treatment Cycles (PK/PD Model)"),
    sidebarLayout(
        sidebarPanel(
            width = 5,
            textInput("fractions", "Ploidy composition fractions",
                      placeholder = "e.g. 0.6,0.3,0.1", value = "0.6,0.3,0.1"),
            numericInput("B0", "Initial tumor burden (cells)", value = 1e7, min = 1, step = 1e6),
            # --- UPDATED: Drug name is now a dropdown menu ---
            selectInput("drug", "Drug name",
                        choices = c("gemcitabine", "bay1895344", "alisertib", "ispinesib", "none"),
                        selected = "gemcitabine"),
            numericInput("cycleDays", "Cycle Length (days)", value = 28, min = 1, max = 100, step = 1),

            actionButton("runCycle", "Run Next Cycle", class = "btn-primary"),
            hr(),
            actionButton("runMCTS", "Run MCTS Optimal Drug Prediction", class = "btn-primary"),
            uiOutput("mcts_prediction"),
            hr(),
            textAreaInput("measurements", "Add tumor-burden measurements (optional)",
                          placeholder = "Format: %burden time\nExample:\n90 5\n80 10\n70 20",
                          height = "120px"),
            actionButton("addMeas", "Add Measurements"),
            actionButton("correctModel", "Correct Model", class = "btn-info"),
            hr(),
            div(class = "small-note",
                tags$ul(
                    tags$li("Try drugs: gemcitabine, alisertib, bay1895344, ispinesib."),
                    tags$li("Add measurements, then click 'Correct Model' to fit."),
                    tags$li("Top-left: total burden; Top-right: composition.")
                )
            )
        ),
        mainPanel(
            width = 7,
            fluidRow(
                column(6, plotOutput("trajPlot_left", height = "300px")),
                column(6, plotOutput("trajPlot_right", height = "300px"))
            ),
            plotOutput("accumPlot", height = "320px"),
            verbatimTextOutput("summaryRates")
        )
    )
)

# ==============================================================================
# ====================== Server (WITH NEW MCTS) ================================
# ==============================================================================

server <- function(input, output, session) {
    all_cycles <- reactiveVal(list())
    last_B <- reactiveVal(NULL)
    current_k_mult <- reactiveVal(1.0)
    prediction_result <- reactiveVal(NULL)

    observeEvent(input$runMCTS, {
        cycles <- all_cycles()
        if (length(cycles) > 0) {
            last_df <- cycles[[length(cycles)]]$df
            b_cols <- grep("^B[0-9]+$", names(last_df), value = TRUE)
            initial_counts <- unlist(tail(last_df, 1)[, b_cols])
        } else {
            fr <- parse_fractions(input$fractions); if (is.null(fr)) return(NULL)
            initial_counts <- fr * input$B0
            names(initial_counts) <- paste0("B", seq_along(fr))
        }

        # Convert initial_counts to ploidy_status format for new MCTS
        # Assuming B1 = 2n, B2 = 3n, B3 = 4n
        ploidy_labels <- seq(2, length.out=length(initial_counts))
        ploidy_status <- as.list(as.numeric(initial_counts))
        names(ploidy_status) <- as.character(ploidy_labels)

        # Parameters for new MCTS
        drugs_mcts <- c("gemcitabine", "bay1895344", "alisertib", "ispinesib", "none")
        total_cycles_plan <- 5
        num_rollouts <- 100  # Reduced from 1000 for faster execution
        rollout_depth <- 10
        min_size <- 1e5
        max_size <- 2e10
        c_param <- sqrt(2)
        d_switch <- 7  # Drug switch time in days

        # Reward function parameters (matching ODE app)
        alpha <- 0.01
        p_order <- 3
        safe_size <- 1e10
        beta <- 1.0

        predicted_sequence <- c()

        withProgress(message = 'Calculating Optimal Sequence...', value = 0, {
            for (step in 1:total_cycles_plan) {
                incProgress(1/total_cycles_plan, detail = paste("Simulating Cycle", step))

                if (sum(unlist(ploidy_status)) < min_size) {
                    predicted_sequence <- c(predicted_sequence, "Extinct")
                    next
                }

                # Run new MCTS
                best_drug <- run_mcts_new(ploidy_status, step-1, drugs_mcts, d_switch,
                                         total_cycles_plan, min_size, max_size,
                                         rollout_depth, num_rollouts, c_param,
                                         alpha, p_order, safe_size, beta)

                predicted_sequence <- c(predicted_sequence, best_drug)

                # Update state with best drug
                result <- simulate_next_state_mcts(ploidy_status, best_drug, d_switch)
                ploidy_status <- result$new_status
            }

            prediction_result(predicted_sequence)
        })
    })

    output$mcts_prediction <- renderUI({
        pred <- prediction_result()
        if (is.null(pred)) return(NULL)

        div(
            style = "background-color: #e8f4f8; padding: 10px; margin-top: 10px; border-radius: 5px;",
            tags$strong("MCTS Predicted Sequence:"),
            tags$ul(
                lapply(seq_along(pred), function(i) {
                    tags$li(paste("Cycle", i, ":", pred[i]))
                })
            )
        )
    })

    observeEvent(input$runCycle, {
        fr <- parse_fractions(input$fractions)
        if (is.null(fr)) {
            showNotification("Invalid ploidy fractions.", type = "error")
            return(NULL)
        }

        B_start <- if (!is.null(last_B())) last_B() else input$B0

        drug_chosen <- input$drug
        cycle_days <- input$cycleDays
        K_val <- 1e9

        df <- run_one_cycle(ploidy_fracs = fr,
                           B0 = B_start,
                           drug_name = drug_chosen,
                           K = K_val,
                           days = cycle_days,
                           dt = 0.1,
                           k_multiplier_base = current_k_mult())

        df_fitted_out <- NULL
        k_mult_fitted <- NULL

        cum_time_start <- 0
        cycles_so_far <- all_cycles()
        if (length(cycles_so_far) > 0) {
            last_cycle <- cycles_so_far[[length(cycles_so_far)]]
            cum_time_start <- max(last_cycle$df$cum_time)
        }
        df$cum_time <- df$time + cum_time_start

        cur_cycle <- list(
            df = df,
            meas = data.frame(time=numeric(0), B=numeric(0), cum_time=numeric(0)),
            drug = drug_chosen,
            r_vec = attr(df, "r_vec"),
            K = K_val,
            phi_params_list = attr(df, "phi_params_list"),
            C_func = attr(df, "C_func"),
            k_multiplier_base = attr(df, "k_multiplier_base"),
            k_multiplier_fitted = k_mult_fitted,
            df_fitted = df_fitted_out
        )

        cycles <- c(all_cycles(), list(cur_cycle))
        all_cycles(cycles)

        last_B(tail(df$B_total, 1))
    })

    observeEvent(input$addMeas, {
        cycles <- all_cycles()
        validate(need(length(cycles) > 0, "Run a cycle first."))

        cur_idx <- length(cycles)
        cur <- cycles[[cur_idx]]

        meas_raw <- parse_measurements(input$measurements)
        if (is.null(meas_raw) || nrow(meas_raw) == 0) return(NULL)

        meas_raw$B <- (meas_raw$percent / 100) * cur$df$B_total[1]
        meas_raw$cum_time <- meas_raw$time + min(cur$df$cum_time)

        cur$meas <- rbind(cur$meas, meas_raw[, c("time","B","cum_time")])

        cycles[[cur_idx]] <- cur
        all_cycles(cycles)
    })

    observeEvent(input$correctModel, {
        cycles <- all_cycles()
        validate(need(length(cycles) > 0, "Run a cycle first."))

        cur_idx <- length(cycles)
        cur <- cycles[[cur_idx]]

        validate(need(nrow(cur$meas) > 0, "Add measurements before correcting model."))

        cost_fn <- function(k_multiplier, data, B0_vec, r_vec, K, phi_params_list, C_func) {

            parms_fit <- list(r_vec = r_vec,
                              K = K,
                              phi_params_list = phi_params_list,
                              C_func = C_func,
                              k_multiplier = k_multiplier)

            sim_times <- seq(0, max(data$time), by = 0.1)
            out <- ode(y = B0_vec, times = sim_times, func = model_ode_fn, parms = parms_fit)
            out_df <- as.data.frame(out)

            b_cols <- grep("^B[0-9]+$", names(out_df), value = TRUE)
            out_df$B_total <- rowSums(out_df[, b_cols, drop = FALSE])

            pred_B <- approx(x = out_df$time, y = out_df$B_total, xout = data$time)$y

            sse <- sum((pred_B - data$B)^2)
            return(sse)
        }

        b_cols <- grep("^B[0-9]+$", names(cur$df), value = TRUE)
        B0_vector <- unlist(cur$df[1, b_cols])

        fit <- optim(
            par = cur$k_multiplier_base,
            fn = cost_fn,
            data = cur$meas,
            B0_vec = B0_vector,
            r_vec = cur$r_vec,
            K = cur$K,
            phi_params_list = cur$phi_params_list,
            C_func = cur$C_func,
            method = "Brent",
            lower = 0.0,
            upper = 5.0
        )

        fitted_k_multiplier <- fit$par

        parms_fitted_final <- list(r_vec = cur$r_vec,
                                   K = cur$K,
                                   phi_params_list = cur$phi_params_list,
                                   C_func = cur$C_func,
                                   k_multiplier = fitted_k_multiplier)

        all_times <- cur$df$time
        df_fitted_out <- ode(y = B0_vector, times = all_times, func = model_ode_fn, parms = parms_fitted_final)
        df_fitted_out <- as.data.frame(df_fitted_out)

        df_fitted_out$B_total <- rowSums(df_fitted_out[, b_cols, drop = FALSE])

        cur$df_fitted <- df_fitted_out
        cur$k_multiplier_fitted <- fitted_k_multiplier
        cycles[[cur_idx]] <- cur
        all_cycles(cycles)

        current_k_mult(fitted_k_multiplier)

        last_B(tail(df_fitted_out$B_total, 1))
    })


    output$trajPlot_left <- renderPlot({
        cycles <- all_cycles()
        validate(need(length(cycles) > 0, "Run a cycle to see results."))
        cur <- cycles[[length(cycles)]]
        df <- cur$df

        max_y <- max(c(df$B_total, cur$meas$B, cur$df_fitted$B_total), na.rm = TRUE)

        par(mar = c(4, 4, 5, 1))
        plot(df$time, df$B_total, type = "l", lwd = 2, col = "firebrick",
             xlab = "Time (days, current cycle)", ylab = "Tumor burden (B)",
             main = paste0("Total Burden\n(Drug: ", cur$drug, ")"),
             ylim = c(0, max_y * 1.05))
        grid()

        if (nrow(cur$meas) > 0) {
            points(cur$meas$time, cur$meas$B, pch = 19, col = "darkred", cex = 1.2)
        }

        if (!is.null(cur$df_fitted)) {
            lines(cur$df_fitted$time, cur$df_fitted$B_total, lty = 2, col = "blue", lwd = 2)
            legend("topright",
                   legend = c("Original", "Fitted"),
                   col = c("firebrick", "blue"),
                   lty = c(1, 2), lwd = 2, bty = "n", cex=0.9)
        }
    })

    output$trajPlot_right <- renderPlot({
        cycles <- all_cycles()
        validate(need(length(cycles) > 0, "Run a cycle to see results."))
        cur <- cycles[[length(cycles)]]
        df <- cur$df

        b_cols <- grep("^B[0-9]+$", names(df), value = TRUE)
        n_cols <- length(b_cols)

        # Handle case with B_total = 0 to avoid NaN
        df_fracs <- df[, b_cols, drop = FALSE] / (df$B_total + 1e-12)

        par(mar = c(4, 4, 3, 1))

        matplot(df$time, df_fracs, type = "l", lty = 1, lwd = 2,
                xlab = "Time (days, current cycle)", ylab = "Fraction of Tumor",
                main = "Ploidy Composition",
                ylim = c(0, 1),
                col = 1:n_cols)
        grid()

        if (!is.null(cur$df_fitted)) {
            df_fitted_fracs <- cur$df_fitted[, b_cols, drop = FALSE] / (cur$df_fitted$B_total + 1e-12)

            matlines(cur$df_fitted$time, df_fitted_fracs,
                     lty = 2, lwd = 2, col = 1:n_cols)
        }

        legend("topright",
               legend = b_cols,
               col = 1:n_cols,
               lty = 1, lwd = 2, cex = 0.9,
               title = "Ploidy Type")

        if (!is.null(cur$df_fitted)) {
            legend("topleft",
                   legend = c("Original", "Fitted"),
                   lty = c(1, 2), lwd = 2,
                   col = "black", bty = "n", cex = 0.9)
        }
    })


    output$accumPlot <- renderPlot({
        cycles <- all_cycles()
        validate(need(length(cycles) > 0, "Run a cycle to see cumulative results."))
        par(mar = c(4, 4, 2, 1))

        all_B_vals <- sapply(cycles, function(x) c(x$df$B_total, x$meas$B, x$df_fitted$B_total))

        max_time <- max(sapply(cycles, function(x) max(x$df$cum_time)))
        max_B <- max(unlist(all_B_vals), na.rm = TRUE)

        plot(NULL, xlim = c(0, max_time), ylim = c(0, max_B * 1.05),
             xlab = "Cumulative Time (days)", ylab = "Tumor burden (B)",
             main = "Accumulated Tumor Burden Over Cycles")

        colors <- rainbow(length(cycles))

        for (i in seq_along(cycles)) {
            lines(cycles[[i]]$df$cum_time, cycles[[i]]$df$B_total, col = colors[i], lwd = 2)

            if (nrow(cycles[[i]]$meas) > 0) {
                points(cycles[[i]]$meas$cum_time, cycles[[i]]$meas$B, pch = 19, col = colors[i], cex = 1.2)
            }

            if (!is.null(cycles[[i]]$df_fitted)) {
                base_cum_t0 <- min(cycles[[i]]$df$cum_time)
                cum_time_fitted <- cycles[[i]]$df_fitted$time + base_cum_t0
                lines(cum_time_fitted, cycles[[i]]$df_fitted$B_total, col = colors[i], lty = 2, lwd = 2)
            }
        }
        legend("topright", legend = paste("Cycle", seq_along(cycles)),
               col = colors, lwd = 2, cex = 0.8)
        grid()
    })

    output$summaryRates <- renderText({
        cycles <- all_cycles()
        if (length(cycles) == 0) return("No cycles run yet.")
        cur <- cycles[[length(cycles)]]

        b_cols <- grep("^B[0-9]+$", names(cur$df), value = TRUE)

        b_initial <- unlist(cur$df[1, b_cols])
        frac_initial <- b_initial / (sum(b_initial) + 1e-12)
        frac_initial_str <- paste(sprintf("%.2f", frac_initial), collapse = ", ")

        b_final <- unlist(tail(cur$df, 1)[, b_cols])
        frac_final <- b_final / (sum(b_final) + 1e-12)
        frac_final_str <- paste(sprintf("%.2f", frac_final), collapse = ", ")

        base_summary <- paste0(
            "Cycle ", length(cycles), " summary (Original Model)\n",
            "Base k_multiplier: ", sprintf("%.3f", cur$k_multiplier_base), "\n",
            "Ploidy types: ", paste(b_cols, collapse = ", "), "\n",
            "Initial Comp: ", frac_initial_str, "\n",
            "Final Comp:   ", frac_final_str, "\n",
            "Final TB: ", sum(b_final)
        )

        if (!is.null(cur$k_multiplier_fitted)) {

            b_final_fit <- unlist(tail(cur$df_fitted, 1)[, b_cols])
            frac_final_fit <- b_final_fit / (sum(b_final_fit) + 1e-12)
            frac_final_fit_str <- paste(sprintf("%.2f", frac_final_fit), collapse = ", ")

            base_summary <- paste0(
                base_summary, "\n\n",
                "--- Fitted Model ---\n",
                "Fitted k_multiplier: ", sprintf("%.3f", cur$k_multiplier_fitted), "\n",
                "Final Comp (Fitted): ", frac_final_fit_str
            )
        }

        return(base_summary)
    })
}

shinyApp(ui, server)