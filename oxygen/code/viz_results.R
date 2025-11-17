suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix)) # Needed for sourced functions
suppressPackageStartupMessages(library(data.table)) # For fread

# ----------------------------------------------------------------------------
## ---- 0. SOURCE FUNCTIONS AND LOAD DATA
# ----------------------------------------------------------------------------
source("code/model_functions.R")

# Load ploidy distribution data
tryCatch({
  x <- data.table::fread("data/ploidy_distribution.csv")
  x$P <- x$ploidy
  x$hypoxia <- grepl("_O",x$id)
  x$ploidy <- rep("2N",nrow(x))
  x$ploidy[grepl("4N",x$id)] <- "4N"
  x$passage <- 0
  x$passage[grepl("_A7",x$id) & x$hypoxia] <- 7
  x$passage[grepl("_A17",x$id)|grepl("_A19",x$id)] <- 17
  
  print("--- Successfully loaded ploidy data (x) ---")
}, error = function(e) {
  stop("Could not load ploidy data: ", e$message)
})

# Load growth rate data
tryCatch({
  m_data <- readRDS("data/fit_g.Rds")
  print("--- Successfully loaded growth rate data (m_data) ---")
}, error = function(e) {
  stop("Could not load 'data/fit_g.Rds'. Error: ", e$message)
})

# ----------------------------------------------------------------------------
## ---- 1. SET UP GLOBAL CONFIGS
# ----------------------------------------------------------------------------

# --- Model & Grid Parameters (GLOBAL) ---
N_UNIT        <- 22L
N_MIN         <- 22L
N_MAX         <- 154L
DT            <- 0.1    # Time step for simulation
POP_GROWTH_FACTOR <- 10.0 # Factor by which pop grows each passage
PASSAGES_TO_RUN   <- 17
REPORT_PASSAGES   <- c(7, 17)

grid_pre  <- N_MIN:N_MAX
grid_post <- N_MIN:N_MAX
R0 <- length(grid_pre)
R1 <- length(grid_post)

# --- Simulation Configurations (GLOBAL) ---
sim_configs <- list(
  list(id = "Sim 1: 2N, O2=1", O2 = 1, init_ploidy = "2N", init_layer = "pre"),
  list(id = "Sim 2: 2N, O2=0", O2 = 0, init_ploidy = "2N", init_layer = "pre"),
  list(id = "Sim 3: 4N, O2=1", O2 = 1, init_ploidy = "4N", init_layer = "post"),
  list(id = "Sim 4: 4N, O2=0", O2 = 0, init_ploidy = "4N", init_layer = "post")
)


# ----------------------------------------------------------------------------
## ---- 2. LOAD OPTIMIZATION RESULTS
# ----------------------------------------------------------------------------

print("--- Loading optimization results from optim_results.RData ---")

tryCatch({
  load("optim_results.RData")
  
  print("--- Running final simulation with optimized parameters... ---")
  
  # 2. Run simulation with loaded parameters
  final_sim_data <- run_all_sims(optim_params_real)
  
  # 3. Print final cost
  final_cost <- calculate_total_cost(final_sim_data, x, m_data, sim_configs, N_UNIT, POP_GROWTH_FACTOR)
  print(paste("--- FINAL OPTIMIZED COST (NLL + SSE) = ", final_cost, "---"))
  print("Optimized Real-Scale Parameters:")
  print(as.data.frame(optim_params_real))
  
  
  # ----------------------------------------------------------------------------
  ## ---- 3. PLOT RESULTS
  # ----------------------------------------------------------------------------
  
  # --- Plot 1: Model vs. Observed Data Comparison ---
  # MODIFICATION: This function now prints the two requested plots directly
  # and does not return a value.
  generate_comparison_plot(final_sim_data, x, sim_configs, N_UNIT, N_MAX)
  
  
  # --- Plot 2: Passage Durations ---
  p_passage_times <- generate_passage_plot(final_sim_data, m_data, sim_configs, POP_GROWTH_FACTOR)
  print(p_passage_times)
  
  # usage (after optim_params_real is loaded):
  p_misseg_interp <- plot_misseg_interp(optim_params_real)
  print(p_misseg_interp)
  
  
}, error = function(e) {
  stop("Could not load 'optim_results.RData'. Please run 'run_optimization.R' first. Error: ", e$message)
})
