suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

# ----------------------------------------------------------------------------
## ---- 1. MODEL & HELPER FUNCTIONS
# ----------------------------------------------------------------------------

#' Mechanistically-based growth rate function
#'
#' @param O2 Oxygen level (0 to 1).
#' @param N Chromosome number.
#' @param R Base growth rate (at P=1, O2=1).
#' @param beta O2-dependent cost parameter. Controls how quickly the
#'             O2-dependent penalty rises with ploidy.
#' @param N_unit Chromosome number of the baseline ploidy (e.g., 22 for P=1).
#' @param eta O2-independent cost parameter. Controls the baseline,
#'            divisive cost of ploidy.
#'
growth_lambda <- function(O2, N, R = 1.0, beta = 0.35, N_unit = 22L, eta = 0.01) {
  
  # P is relative ploidy (P=1 for N=22, P=2 for N=44, etc.)
  # We use pmax to avoid P < 1 if N < N_unit
  P <- pmax(1.0, N / N_unit)
  
  # --- O2-dependent Component (Numerator) ---
  #
  # 1. Calculate the O2-dependent *penalty*.
  # This is a saturating (hyperbolic) function of (P-1).
  # - If P=1, penalty_P = 0.
  # - As P -> Inf, penalty_P -> 1.
  # 'beta' now controls the "steepness" of this penalty.
  # A higher beta means the penalty saturates at lower ploidy.
  penalty_P <- (beta * (P - 1)) / (1 + beta * (P - 1))
  
  # 2. Modulate this penalty by oxygen availability.
  # (1 - O2) is the "anoxia factor": 0 at full O2, 1 at zero O2.
  # O2_factor is the growth multiplier (1 at full O2, decreases at low O2)
  O2_factor <- 1 - (1 - O2) * penalty_P
  
  # --- O2-independent Component (Denominator) ---
  #
  # This is the original baseline cost of ploidy (e.g., metabolic load,
  # mitotic problems) that applies regardless of oxygen.
  # - If P=1, baseline_cost_denom = 1 (no cost).
  # - If P > 1, baseline_cost_denom > 1 (applies a cost).
  baseline_cost_denom <- (1 + eta * (P - 1))
  
  # --- Combine ---
  # G = R * (O2-dependent part) / (O2-independent part)
  G <- R * O2_factor / baseline_cost_denom
  
  # Ensure non-negative growth
  pmax(G, 0)
}

.pr_delta_vec <- function(N, p, eps_tail = 1e-8, mr_lethality = 0.9){
  if (p <= 0 || N == 0) { out <- c("0"=1); attr(out,"mass_dropped") <- 0; return(out) }
  sd <- sqrt(N * p)
  if (sd == 0) { out <- c("0"=1); attr(out,"mass_dropped") <- 0; return(out) }
  z  <- stats::qnorm(1 - eps_tail/2)
  T  <- min(N, max(0L, ceiling(z * sd)))
  ts <- (-T):T
  out <- numeric(length(ts))
  for (idx in seq_along(ts)) {
    t  <- ts[idx]
    ks <- seq.int(abs(t), N, by = 2)
    if (length(ks)) {
      pk <- stats::dbinom(ks, size = N, prob = p)
      m  <- (ks + t)/2
      qm <- stats::dbinom(m, size = ks, prob = 0.5)
      surv_m <- (1 - mr_lethality)^ks
      out[idx] <- sum(pk * qm * surv_m)
    }
  }
  names(out) <- as.character(ts)
  attr(out,"mass_dropped") <- max(0, 1 - sum(out))
  out
}

.build_B_total <- function(Nmin, Nmax, p_vec, mr_lethality = 0.9,
                           boundary = c("drop","absorb_minmax"),
                           eps_tail = 1e-8, return_sparse = TRUE){
  boundary <- match.arg(boundary)
  R <- Nmax - Nmin + 1L
  if (length(p_vec) == 1L) p_vec <- rep(p_vec, R)
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (col in seq_len(R)) {
    N  <- Nmin + col - 1L
    pN <- p_vec[col]
    prDelta <- .pr_delta_vec(N, pN, eps_tail = eps_tail, mr_lethality = mr_lethality)
    ts <- as.integer(names(prDelta)); pr <- as.numeric(prDelta)
    for (k in seq_along(ts)) {
      t <- ts[k]; w <- pr[k]; if (w == 0) next
      if (t == 0L) {
        Np <- N; val <- 2*w
        if (Np < Nmin || Np > Nmax) {
          if (boundary == "absorb_minmax"){
            Np2 <- max(min(Np,Nmax),Nmin); ii <- c(ii, Np2-Nmin+1L); jj <- c(jj, col); xx <- c(xx, val)
          }
        } else {
          ii <- c(ii, N - Nmin + 1L); jj <- c(jj, col); xx <- c(xx, val)
        }
      } else {
        for (Np in c(N+t, N-t)) {
          if (Np < Nmin || Np > Nmax) {
            if (boundary == "absorb_minmax"){
              Np2 <- max(min(Np,Nmax),Nmin); ii <- c(ii, Np2-Nmin+1L); jj <- c(jj, col); xx <- c(xx, w)
            }
          } else {
            ii <- c(ii, Np - Nmin + 1L); jj <- c(jj, col); xx <- c(xx, w)
          }
        }
      }
    }
  }
  sparseMatrix(i = ii, j = jj, x = xx, dims = c(R,R), repr = "C")
}

.build_B_WGD <- function(N0min, N0max, N1min, N1max,
                         boundary = c("drop","absorb_minmax"),
                         return_sparse = TRUE){
  boundary <- match.arg(boundary)
  R0 <- N0max - N0min + 1L
  R1 <- N1max - N1min + 1L
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (col in seq_len(R0)) {
    N0 <- N0min + col - 1L
    Np <- 2L * N0
    val <- 2.0
    if (Np < N1min || Np > N1max) {
      if (boundary == "absorb_minmax"){
        Np2 <- max(min(Np, N1max), N1min)
        ii <- c(ii, Np2 - N1min + 1L)
        jj <- c(jj, col)
        xx <- c(xx, val)
      }
    } else {
      ii <- c(ii, Np - N1min + 1L)
      jj <- c(jj, col)
      xx <- c(xx, val)
    }
  }
  sparseMatrix(i = ii, j = jj, x = xx, dims = c(R1, R0), repr = "C")
}

.build_G_with_WGD <- function(
    N0min, N0max, lambda0_vec, p0_vec, wgd_prob_vec,
    N1min, N1max, lambda1_vec, p1_vec,
    mr_lethality0 = 0.9, mr_lethality1 = 0.9,
    boundary = "absorb_minmax", eps_tail = 1e-8
){
  R0 <- N0max - N0min + 1L
  R1 <- N1max - N1min + 1L
  if (length(lambda0_vec) == 1L) lambda0_vec <- rep(lambda0_vec, R0)
  if (length(p0_vec) == 1L) p0_vec <- rep(p0_vec, R0)
  if (length(wgd_prob_vec) == 1L) wgd_prob_vec <- rep(wgd_prob_vec, R0)
  if (length(lambda1_vec) == 1L) lambda1_vec <- rep(lambda1_vec, R1)
  if (length(p1_vec) == 1L) p1_vec <- rep(p1_vec, R1)
  wgd_prob_vec <- pmin(pmax(wgd_prob_vec, 0), 1)
  
  B0 <- .build_B_total(N0min, N0max, p_vec = p0_vec, mr_lethality = mr_lethality0, boundary = boundary, eps_tail = eps_tail)
  B1 <- .build_B_total(N1min, N1max, p_vec = p1_vec, mr_lethality = mr_lethality1, boundary = boundary, eps_tail = eps_tail)
  BW <- .build_B_WGD (N0min, N0max, N1min, N1max, boundary = boundary)
  
  L0 <- Diagonal(x = lambda0_vec)
  L1 <- Diagonal(x = lambda1_vec)
  S0 <- Diagonal(x = (1 - wgd_prob_vec))
  SW <- Diagonal(x = wgd_prob_vec)
  
  UL  <- (B0 %*% S0) %*% L0 - L0
  LR  <- (B1 %*% L1) - L1
  LL  <- (BW %*% SW) %*% L0
  UR  <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(R0, R1))
  
  G <- rbind(cbind(UL, UR), cbind(LL, LR))
  G
}

step_dt <- function(G, x, dt, steps = 1L, normalize = FALSE){
  I <- Diagonal(n = nrow(G))
  A <- I + dt * G
  v <- as.numeric(x)
  for (k in seq_len(steps)) {
    v <- as.numeric(A %*% v)
    if (normalize) v <- v / sum(v)
  }
  v
}

#' Helper to create initial distribution from ploidy data
create_initial_dist <- function(ploidy_values, N_grid, N_unit = 22L) {
  N_values <- round(ploidy_values * N_unit)
  N_counts <- table(N_values)
  N_fracs <- N_counts / sum(N_counts)
  x_vec <- rep(0, length(N_grid))
  names(x_vec) <- N_grid
  valid_names <- names(N_fracs)[names(N_fracs) %in% names(x_vec)]
  x_vec[valid_names] <- N_fracs[valid_names]
  x_vec
}

# ----------------------------------------------------------------------------
## ---- 3. SIMULATION RUNNER FUNCTION
# ----------------------------------------------------------------------------

#' Runs all 4 simulations for a given parameter set
#' @param run_params A list of REAL-SCALE parameters
#' @return A list containing all_dists and all_passage_times data.frames
#' @details This function is designed to be called from a parallel context
#'          and relies on several GLOBAL variables being exported:
#'          x, sim_configs, N_UNIT, N_MIN, N_MAX, DT, 
#'          POP_GROWTH_FACTOR, PASSAGES_TO_RUN, REPORT_PASSAGES, 
#'          grid_pre, grid_post, R0, R1
run_all_sims <- function(run_params) {
  
  all_results_list <- list()
  passage_times <- list()
  
  # --- Get Initial Ploidy Data ---
  init_P_2N <- x$P[x$passage == 0 & x$ploidy == "2N"]
  init_P_4N <- x$P[x$passage == 0 & x$ploidy == "4N"]
  
  for (sim in sim_configs) {
    
    O2_LEVEL <- sim$O2
    pmis_const <- ifelse(O2_LEVEL == 1, run_params$pmis_O2_1, run_params$pmis_O2_0)
    
    # --- Create initial state vectors ---
    x0_pre  <- rep(0, R0); names(x0_pre)  <- grid_pre
    x0_post <- rep(0, R1); names(x0_post) <- grid_post
    
    if (sim$init_layer == "pre") {
      init_P_values <- if (sim$init_ploidy == "2N") init_P_2N else init_P_4N
      x0_pre <- create_initial_dist(init_P_values, grid_pre, N_UNIT)
    } else {
      init_P_values <- if (sim$init_ploidy == "2N") init_P_2N else init_P_4N
      x0_post <- create_initial_dist(init_P_values, grid_post, N_UNIT)
    }
    
    x_current <- c(x0_pre, x0_post)
    x_current <- x_current / sum(x_current) # Start with total pop = 1
    
    # --- Build the Generator Matrix G ---
    lambda0_vec <- growth_lambda(O2_LEVEL, grid_pre,  R = run_params$R, beta = run_params$beta, eta = run_params$eta, N_unit = N_UNIT)
    lambda1_vec <- growth_lambda(O2_LEVEL, grid_post, R = run_params$R, beta = run_params$beta, eta = run_params$eta, N_unit = N_UNIT)
    
    G <- .build_G_with_WGD(
      N0min = N_MIN, N0max = N_MAX,
      lambda0_vec = lambda0_vec,
      p0_vec = pmis_const,
      wgd_prob_vec = run_params$pwgd,
      N1min = N_MIN, N1max = N_MAX,
      lambda1_vec = lambda1_vec,
      p1_vec = pmis_const,
      mr_lethality0 = run_params$mr_lethality0,
      mr_lethality1 = run_params$mr_lethality1
    )
    
    # --- Run Passage Loop ---
    sim_passage_times <- numeric(PASSAGES_TO_RUN)
    
    for (p in 1:PASSAGES_TO_RUN) {
      pop_start  <- sum(x_current)
      pop_target <- pop_start * POP_GROWTH_FACTOR
      time_in_passage <- 0.0
      
      while (sum(x_current) < pop_target) {
        x_current <- step_dt(G, x_current, DT, 1L)
        time_in_passage <- time_in_passage + DT
        if (sum(x_current) < pop_start * 1e-3 || time_in_passage > 1000) {
          break 
        }
      }
      sim_passage_times[p] <- time_in_passage
      
      if (p %in% REPORT_PASSAGES) {
        pop_total <- sum(x_current)
        dist_df <- data.frame(
          sim_id  = sim$id,
          passage = p,
          layer   = c(rep("pre", R0), rep("post", R1)),
          N       = c(grid_pre, grid_post),
          fraction = x_current / pop_total
        )
        all_results_list[[length(all_results_list) + 1]] <- dist_df
      }
      x_current <- x_current / sum(x_current) * pop_start
    }
    passage_times[[sim$id]] <- data.frame(
      sim_id = sim$id,
      passage = 1:PASSAGES_TO_RUN,
      duration = sim_passage_times
    )
  }
  
  all_dists <- do.call(rbind, all_results_list)
  all_passage_times <- do.call(rbind, passage_times)
  
  return(list(all_dists = all_dists, all_passage_times = all_passage_times))
}

# ----------------------------------------------------------------------------
## ---- 4. COST FUNCTION CALCULATOR
# ----------------------------------------------------------------------------

#' Internal helper for cost calculation
calculate_nll <- function(pred_dist, obs_P_values, N_unit, sim_configs) {
  if (length(obs_P_values) == 0) {
    warning(paste("No observed data for", unique(pred_dist$sim_id), "p", unique(pred_dist$passage)))
    return(0)
  }
  obs_N <- round(obs_P_values * N_unit)
  obs_counts <- as.data.frame(table(N = obs_N))
  obs_counts$N <- as.integer(as.character(obs_counts$N))
  
  # --- MODIFICATION ---
  # Aggregate pre and post layers by N to get total ploidy distribution
  pred_dist_combined <- pred_dist %>%
    group_by(N) %>%
    summarise(fraction = sum(fraction), .groups = 'drop')
  
  # Join observations with combined predictions
  comp_df <- merge(obs_counts, pred_dist_combined, by = "N", all.x = TRUE)
  # --- END MODIFICATION ---
  
  min_prob <- 1e-10
  comp_df$fraction[is.na(comp_df$fraction)] <- min_prob
  
  nll <- -sum(comp_df$Freq * log(comp_df$fraction + min_prob))
  nll
}

#' Calculates total cost (NLL + SSE) from simulation output
calculate_total_cost <- function(sim_output, x_data, m_data, sim_configs, N_unit, pop_growth_factor) {
  
  # 1. Cost from Ploidy Distributions (NLL)
  all_dists <- sim_output$all_dists
  cost_list_nll <- list()
  
  # Sim 1
  pred_s1_p17 <- subset(all_dists, sim_id == "Sim 1: 2N, O2=1" & passage == 17)
  obs_s1_p0   <- x_data$P[x_data$passage == 0 & x_data$ploidy == "2N"]
  cost_list_nll$s1_p17 <- calculate_nll(pred_s1_p17, obs_s1_p0, N_unit, sim_configs)
  
  # Sim 2
  pred_s2_p7  <- subset(all_dists, sim_id == "Sim 2: 2N, O2=0" & passage == 7)
  obs_s2_p7   <- x_data$P[x_data$passage == 7 & x_data$ploidy == "2N"]
  cost_list_nll$s2_p7 <- calculate_nll(pred_s2_p7, obs_s2_p7, N_unit, sim_configs)
  pred_s2_p17 <- subset(all_dists, sim_id == "Sim 2: 2N, O2=0" & passage == 17)
  obs_s2_p17  <- x_data$P[x_data$passage == 17 & x_data$ploidy == "2N"]
  cost_list_nll$s2_p17 <- calculate_nll(pred_s2_p17, obs_s2_p17, N_unit, sim_configs)
  
  # Sim 3
  pred_s3_p17 <- subset(all_dists, sim_id == "Sim 3: 4N, O2=1" & passage == 17)
  obs_s3_p0   <- x_data$P[x_data$passage == 0 & x_data$ploidy == "4N"]
  cost_list_nll$s3_p17 <- calculate_nll(pred_s3_p17, obs_s3_p0, N_unit, sim_configs)
  
  # Sim 4
  pred_s4_p7  <- subset(all_dists, sim_id == "Sim 4: 4N, O2=0" & passage == 7)
  obs_s4_p7   <- x_data$P[x_data$passage == 7 & x_data$ploidy == "4N"]
  cost_list_nll$s4_p7 <- calculate_nll(pred_s4_p7, obs_s4_p7, N_unit, sim_configs)
  pred_s4_p17 <- subset(all_dists, sim_id == "Sim 4: 4N, O2=0" & passage == 17)
  obs_s4_p17  <- x_data$P[x_data$passage == 17 & x_data$ploidy == "4N"]
  cost_list_nll$s4_p17 <- calculate_nll(pred_s4_p17, obs_s4_p17, N_unit, sim_configs)
  
  total_cost_nll <- sum(unlist(cost_list_nll))
  
  # 2. Cost from Growth Rates (SSE)
  sim_passage_times <- sim_output$all_passage_times
  
  # Calculate simulated growth rate
  sim_passage_times$g_sim <- log(pop_growth_factor) / sim_passage_times$duration
  
  # Create mapping from sim_id to ploidy/o2
  mapping <- data.frame(
    sim_id = sapply(sim_configs, `[[`, "id"),
    o2 = c(1, 0, 1, 0),
    ploidy = c("2N", "2N", "4N", "4N")
  )
  
  # Join simulated and observed growth data
  comp_g <- inner_join(sim_passage_times, mapping, by = "sim_id")
  m_data_renamed <- dplyr::rename(m_data, g_obs = g)
  comp_g <- inner_join(comp_g, m_data_renamed, by = c("ploidy", "o2", "passage"))
  
  # Handle potential Inf/-Inf from duration=0 or duration=Inf
  comp_g <- dplyr::filter(comp_g, is.finite(g_sim) & is.finite(g_obs))
  
  # Calculate Sum of Squared Errors
  total_cost_g <- sum((comp_g$g_sim - comp_g$g_obs)^2)
  
  # 3. Total Cost
  total_cost <- total_cost_nll + total_cost_g
  
  # Print breakdown for debugging
  cat(sprintf("  NLL Cost: %.2f | Growth Cost (SSE): %.2f | Total Cost: %.2f\n", 
              total_cost_nll, total_cost_g, total_cost))
  
  return(total_cost)
}


# ----------------------------------------------------------------------------
## ---- 5. PLOTTING FUNCTIONS
# ----------------------------------------------------------------------------

#' Generates the Model vs. Observed data comparison plot
#' Generates the Model vs. Observed data comparison plots (Normoxia and Anoxia)
#' This function now prints two plots directly as requested.
#' Generates the Model vs. Observed data comparison plots (Normoxia and Anoxia)
#' This function now prints two plots directly as requested.
#' Generates the Model vs. Observed data comparison plots (2N and 4N)
#' This function now prints two plots directly as requested.
#' Generates the Model vs. Observed data comparison plots (2N and 4N)
#' This function now prints two plots directly as requested, using boxplots
#' for both observed data and sampled model data.
#' Generates the Model vs. Observed data comparison plots
#' This function now prints one combined plot, using boxplots
#' for both observed data and sampled model data.
#' Generates the Model vs. Observed data comparison plots
#' This function now prints one combined plot, using dodged boxplots
#' to show time progression on the x-axis.
#' Generates the Model vs. Observed data comparison plots
#' This function now prints one combined plot, using dodged boxplots
#' faceted by ploidy (rows) and oxygen (cols).
generate_comparison_plot <- function(final_sim_data, x_data, sim_configs, N_unit, N_max) {
  
  # Number of samples to draw from the model distribution for the boxplot
  n_samples_for_model <- 10000
  
  # --- 1. Prepare Model Data ---
  
  # A) Get the model's probability distributions
  model_distributions <- final_sim_data$all_dists %>%
    mutate(
      O2 = ifelse(grepl("O2=0", sim_id), 0, 1),
      P = N / N_unit,
      week = passage, 
      ploidy = ifelse(grepl("2N", sim_id), "2N", "4N"),
      O2_label = ifelse(O2 == 1, "Normoxia", "Anoxia")
    ) %>%
    # Aggregate pre and post layers to get total ploidy distribution
    group_by(sim_id, O2, O2_label, P, week, ploidy, passage) %>%
    summarise(fraction = sum(fraction), .groups = 'drop') %>%
    # Ensure fractions sum to 1 within each group (handles potential float errors)
    group_by(O2_label, week, ploidy) %>%
    mutate(fraction = fraction / sum(fraction)) %>%
    ungroup()
  
  # B) Create a sampled dataset from the model distributions
  model_data_sampled <- model_distributions %>%
    group_by(O2_label, week, ploidy) %>%
    reframe(
      P = sample(
        x = P, 
        size = n_samples_for_model, 
        replace = TRUE, 
        prob = fraction
      ),
      source = "Model"
    )
  
  # --- 2. Prepare Observed Data ---
  # This is the experimental data
  
  # Separate P0 data, as it's the start for *both* conditions
  obs_P0 <- x_data %>% filter(passage == 0)
  obs_P_other <- x_data %>% filter(passage != 0)
  
  # Duplicate P0 data for both O2=1 (Normoxia) and O2=0 (Anoxia)
  obs_P0_norm <- obs_P0 %>% mutate(O2 = 1, week = 0)
  obs_P0_anox <- obs_P0 %>% mutate(O2 = 0, week = 0)
  
  # All other observed data (P7, P17) is from the Anoxia (O2=0) experiment
  obs_P_other <- obs_P_other %>% mutate(O2 = 0, week = passage)
  
  # Combine all observed data
  obs_plot_data <- bind_rows(obs_P0_norm, obs_P0_anox, obs_P_other) %>%
    mutate(source = "Observed") %>%
    mutate(O2_label = ifelse(O2 == 1, "Normoxia", "Anoxia")) %>%
    select(P, ploidy, O2, O2_label, week, source, passage)
  
  # --- 3. Combine and Factor Data ---
  all_plot_data <- bind_rows(model_data_sampled, obs_plot_data) %>%
    mutate(
      O2_label = factor(O2_label, levels = c("Normoxia", "Anoxia")),
      source = factor(source, levels = c("Observed", "Model"))
    )
  
  # --- 4. Create Combined Plot ---
  
  # *** MODIFICATION 2: Mapped 'color = source' ***
  p_combined <- ggplot(all_plot_data, aes(x = factor(week), y = P, fill = source, color = source)) +
    geom_boxplot(
      outlier.shape = NA, 
      position = position_dodge(preserve = "single")
    ) +
    # *** MODIFICATION 1: Facet by ploidy (rows) ~ O2_label (cols) ***
    facet_grid(
      ploidy ~ O2_label, 
      scales = "free_y"
    ) +
    labs(
      x = "Week",
      y = "Ploidy",
      # *** MODIFICATION 2: Updated legend labs ***
      fill = "Data Source",
      color = "Data Source"
    ) +
    scale_fill_manual(values = c("Observed" = "grey", "Model" = "skyblue")) +
    # *** MODIFICATION 2: Added color scale ***
    scale_color_manual(values = c("Observed" = "grey", "Model" = "skyblue")) +
    scale_y_continuous(limits = c(1.8, 4.2)) + 
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0), 
      legend.position = "bottom"
    )
  
  print(p_combined)
  
  invisible(p_combined)
}

#' Generates the Passage Durations plot (Model vs. Observed)
generate_passage_plot <- function(final_sim_data, m_data, sim_configs, pop_growth_factor) {
  
  # Simulated data
  sim_passage_times <- final_sim_data$all_passage_times
  sim_passage_times$g_sim <- log(pop_growth_factor) / sim_passage_times$duration
  
  # Mapping
  mapping <- data.frame(
    sim_id = sapply(sim_configs, `[[`, "id"),
    o2 = c(1, 0, 1, 0),
    ploidy = c("2N", "2N", "4N", "4N")
  )
  
  # Combine sim and obs data
  plot_data_sim <- inner_join(sim_passage_times, mapping, by = "sim_id") %>%
    dplyr::rename(g = g_sim) %>%
    mutate(source = "Model")
  
  plot_data_obs <- inner_join(dplyr::rename(m_data, g_obs = g), mapping, by = c("ploidy", "o2")) %>%
    dplyr::rename(g = g_obs) %>%
    mutate(source = "Observed")
  
  plot_data_combined <- bind_rows(
    select(plot_data_sim, sim_id, passage, g, source),
    select(plot_data_obs, sim_id, passage, g, source)
  )
  
  p_passage_times <- ggplot(plot_data_combined, aes(x = passage, y = g, color = sim_id)) +
    geom_line(data = . %>% dplyr::filter(source == "Model"), linewidth = 1) +
    geom_point(data = . %>% dplyr::filter(source == "Observed"), size = 3, shape = 1) +
    labs(
      title = "Passage Growth Rate (Model = line, Observed = point)",
      x = "Passage Number",
      y = "Growth Rate (g = log(10) / duration)",
      color = "Simulation"
    ) +
    theme_bw()
  
  return(p_passage_times)
}

# Build explicit initial state vectors (no globals, no data lookups)
# - pass either:
#   * ploidy = 2 or 4  (for a pure delta at N=44 or 88 on the chosen layer), OR
#   * init_pre/init_post numeric named vectors summing to 1 (fractions on grids)
# total_size is the starting population size (cells)
make_init_state <- function(grid_pre, grid_post,
                            ploidy = c(2,4), layer = c("pre","post"),
                            init_pre = NULL, init_post = NULL,
                            N_UNIT = 22L, total_size = 1e6) {
  layer  <- match.arg(layer)
  ploidy <- match.arg(as.character(ploidy), choices = c("2","4"))
  Pnum   <- as.numeric(ploidy)
  
  R0 <- length(grid_pre); R1 <- length(grid_post)
  x_pre  <- rep(0, R0); names(x_pre)  <- grid_pre
  x_post <- rep(0, R1); names(x_post) <- grid_post
  
  if (!is.null(init_pre) || !is.null(init_post)) {
    if (!is.null(init_pre))  x_pre[names(init_pre)]   <- as.numeric(init_pre)
    if (!is.null(init_post)) x_post[names(init_post)] <- as.numeric(init_post)
  } else {
    N_delta <- as.integer(Pnum * N_UNIT)
    if (layer == "pre") {
      stopifnot(N_delta %in% grid_pre)
      x_pre[as.character(N_delta)] <- 1
    } else {
      stopifnot(N_delta %in% grid_post)
      x_post[as.character(N_delta)] <- 1
    }
  }
  # normalize to fractions, then scale to size
  s <- sum(x_pre) + sum(x_post); if (s <= 0) stop("Init mass is zero.")
  v <- c(x_pre, x_post) / s * total_size
  return(v)
}


run_in_vivo_crowd <- function(run_params,
                              O2_schedule = list(c(t0=0, t1=Inf, O2=1.0)),
                              T_end = 28, sample_days = c(0,7,14,21,28),
                              N_UNIT = 22L, DT = 0.1,
                              K = 1e9, crowding = c("logistic","gompertz"),
                              grid_pre = get("grid_pre", inherits=TRUE),
                              grid_post = get("grid_post", inherits=TRUE),
                              init_state) {
  
  crowding <- match.arg(crowding)
  
  # --- grids & bounds (needed by the G builder closure) ---
  R0 <- length(grid_pre);  R1 <- length(grid_post)
  N0min <- min(grid_pre);  N0max <- max(grid_pre)
  N1min <- min(grid_post); N1max <- max(grid_post)
  
  # --- init state (explicit; no globals) ---
  stopifnot(length(init_state) == (R0 + R1))
  v <- as.numeric(init_state)
  
  # --- helpers ---
  get_O2 <- function(t){
    for (seg in O2_schedule) if (t >= seg["t0"] && t < seg["t1"]) return(as.numeric(seg["O2"]))
    as.numeric(O2_schedule[[length(O2_schedule)]]["O2"])
  }
  pmis_of_O2 <- function(O2){
    # log-linear between fitted endpoints
    logp <- (1 - O2) * log10(run_params$pmis_O2_0) + O2 * log10(run_params$pmis_O2_1)
    10^logp
  }
  
  # cache G by rounded O2
  G_cache <- new.env(parent = emptyenv())
  build_G_for_O2 <- function(O2){
    key <- sprintf("%.3f", O2)
    if (!exists(key, envir = G_cache)) {
      lambda0 <- growth_lambda(O2, grid_pre,  R=run_params$R, beta=run_params$beta, eta=run_params$eta, N_unit=N_UNIT)
      lambda1 <- growth_lambda(O2, grid_post, R=run_params$R, beta=run_params$beta, eta=run_params$eta, N_unit=N_UNIT)
      p_mis   <- pmis_of_O2(O2)
      G <- .build_G_with_WGD(
        N0min, N0max, lambda0, p0_vec = p_mis, wgd_prob_vec = run_params$pwgd,
        N1min, N1max, lambda1, p1_vec = p_mis,
        mr_lethality0 = run_params$mr_lethality0,
        mr_lethality1 = run_params$mr_lethality1
      )
      assign(key, G, envir = G_cache)
    }
    get(key, envir = G_cache)
  }
  
  crowd <- function(Ntot){
    if (crowding == "logistic") return(max(0, 1 - Ntot / K))
    if (crowding == "gompertz") return(exp(- Ntot / K))
  }
  
  I <- Diagonal(n = length(v))
  times <- seq(0, T_end, by = DT)
  snapshots <- list()
  size_trace <- data.frame(day = 0, Ntot = sum(v))
  
  for (t in times){
    if (t %in% sample_days){
      snapshots[[as.character(t)]] <- data.frame(
        day = t,
        layer = c(rep("pre", R0), rep("post", R1)),
        N = c(grid_pre, grid_post),
        fraction = v / sum(v),
        pop = sum(v)
      )
    }
    if (t >= T_end) break
    O2t <- get_O2(t)
    G   <- build_G_for_O2(O2t)
    cfac <- crowd(sum(v))
    v <- as.numeric((I + DT * (cfac * G)) %*% v)
    size_trace <- rbind(size_trace, data.frame(day = t + DT, Ntot = sum(v)))
    if (sum(v) <= 1e-9) break
  }
  
  list(
    all_dists  = do.call(rbind, snapshots),
    tumor_size = size_trace
  )
}

plot_misseg_interp <- function(par, o2_ref = 20.5){
  stopifnot(is.numeric(par$pmis_O2_0), is.numeric(par$pmis_O2_1),
            par$pmis_O2_0 > 0, par$pmis_O2_1 > 0)
  p0 <- as.numeric(par$pmis_O2_0); p1 <- as.numeric(par$pmis_O2_1)
  O2 <- seq(0, 1, length.out = 401)
  p  <- exp((1 - O2) * log(p0) + O2 * log(p1))  # log-linear interpolation (base-agnostic)
  df <- data.frame(O2 = O2, O2_pct = O2 * o2_ref, p = p)
  ggplot(df, aes(O2_pct, p)) +
    geom_line(linewidth = 1, color = "black") +
    geom_point(data = df[c(1, nrow(df)), ], size = 2, color = "red") +
    labs(x = "Oxygen (%)", y = "Missegregation rate") +
    theme_bw()
}


