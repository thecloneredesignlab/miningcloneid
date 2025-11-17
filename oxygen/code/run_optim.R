suppressPackageStartupMessages(library(DEoptim))
suppressPackageStartupMessages(library(Matrix)) # Ensure Matrix is loaded for objective_function
suppressPackageStartupMessages(library(ggplot2)) # For saving plots
suppressPackageStartupMessages(library(tidyr))   # For saving plots
suppressPackageStartupMessages(library(dplyr))   # For data manipulation

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
  print(head(m_data))
}, error = function(e) {
  stop("Could not load 'data/fit_g.Rds'. Error: ", e$message)
})


# ----------------------------------------------------------------------------
## ---- 1. SET UP GLOBAL CONFIGS AND OPTIMIZATION
# ----------------------------------------------------------------------------

# --- Model & Grid Parameters (GLOBAL) ---
N_UNIT       <- 22L
N_MIN        <- 22L
N_MAX        <- 154L
DT           <- 0.1   # Time step for simulation
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

# --- Parameter Definitions (GLOBAL) ---
param_names <- c("R", "beta", "log10_pwgd", "mr_lethality0", "mr_lethality1", 
                 "log10_pmis_O2_1", "log10_pmis_O2_0", "log10_eta")

# --- Variables to export for parallel processing ---
vars_to_export <- c(
  "x", "m_data", "sim_configs", "param_names", "N_UNIT", "N_MIN", "N_MAX", "DT", 
  "POP_GROWTH_FACTOR", "PASSAGES_TO_RUN", "REPORT_PASSAGES", 
  "grid_pre", "grid_post", "R0", "R1",
  "run_all_sims", "calculate_total_cost", "calculate_nll", 
  "growth_lambda", ".pr_delta_vec", ".build_B_total", ".build_B_WGD", 
  ".build_G_with_WGD", "step_dt", "create_initial_dist"
)

# --- Iterative Optimization Parameters ---
TOTAL_GENERATIONS    <- 100 # Total number of generations to run
GENERATIONS_PER_LOOP <- 10  # Run this many generations, then plot

# Create plot directory
dir.create("plots/tmp", recursive = TRUE, showWarnings = FALSE)


# ----------------------------------------------------------------------------
## ---- 2. OPTIMIZATION OBJECTIVE FUNCTION
# ----------------------------------------------------------------------------

#' Objective function for optim()
#' Takes a transformed parameter vector, runs the model, returns the cost
objective_function <- function(par_transformed) {
  
  names(par_transformed) <- param_names
  
  # 1. Un-transform parameters to their real scale
  run_params <- list(
    R             = par_transformed["R"],
    beta          = par_transformed["beta"],
    pwgd          = 10^par_transformed["log10_pwgd"],
    mr_lethality0 = par_transformed["mr_lethality0"],
    mr_lethality1 = par_transformed["mr_lethality1"],
    pmis_O2_1     = 10^par_transformed["log10_pmis_O2_1"],
    pmis_O2_0     = 10^par_transformed["log10_pmis_O2_0"],
    eta           = 10^par_transformed["log10_eta"]
  )
  
  # 2. Run all simulations
  sim_data <- run_all_sims(run_params)
  
  # 3. Calculate total cost
  # 'x', 'm_data', 'sim_configs', etc. are exported to the cluster environment
  total_cost <- calculate_total_cost(sim_data, x, m_data, sim_configs, N_UNIT, POP_GROWTH_FACTOR)
  
  return(total_cost)
}

# ----------------------------------------------------------------------------
## ---- 3. RUN OPTIMIZATION
# ----------------------------------------------------------------------------

# --- Lower Bounds (Transformed Scale) ---
lower_bounds <- c(
  R               = 0.5,
  beta            = 0.0,
  log10_pwgd      = log10(1e-5),
  mr_lethality0   = 0.0,
  mr_lethality1   = 0.0,
  log10_pmis_O2_1 = log10(1e-5),
  log10_pmis_O2_0 = log10(1e-5),
  log10_eta       = log10(1e-5)
)

# --- Upper Bounds (Transformed Scale) ---
upper_bounds <- c(
  R               = 1.0,
  beta            = 1.0,
  log10_pwgd      = log10(1e-2),
  mr_lethality0   = 1.0,
  mr_lethality1   = 1.0,
  log10_pmis_O2_1 = log10(1e-2),
  log10_pmis_O2_0 = log10(1e-2),
  log10_eta       = log10(0.1)
)

print("--- Starting Iterative Optimization with DEoptim ---")

current_pop <- NULL # Start with no initial population
optim_result <- NULL      # Store the result

num_loops <- TOTAL_GENERATIONS / GENERATIONS_PER_LOOP

for (i in 1:num_loops) {
  current_generation_chunk <- i * GENERATIONS_PER_LOOP
  cat(paste("\n--- Running generations", (current_generation_chunk - GENERATIONS_PER_LOOP + 1), 
            "to", current_generation_chunk, "---\n"))
  
  optim_result <- DEoptim(
    fn = objective_function,
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(
      trace = 1, # Print progress
      parallelType = "parallel",
      packages = c("Matrix", "dplyr", "tidyr"), # Add dplyr/tidyr for cost function
      parVar = vars_to_export,
      itermax = GENERATIONS_PER_LOOP, # Run for N generations
      initialpop = current_pop       # Resume from last population
    )
  )
  
  # Update current population for the next loop
  current_pop <- optim_result$member$pop
  
  # --- Generate and Save Intermediate Plots ---
  cat("--- Generating plots for generation", current_generation_chunk, "---\n")
  
  # Get best parameters so far
  optim_par_transformed <- optim_result$optim$bestmem
  names(optim_par_transformed) <- param_names
  optim_params_real <- list(
    R             = optim_par_transformed["R"],
    beta          = optim_par_transformed["beta"],
    pwgd          = 10^optim_par_transformed["log10_pwgd"],
    mr_lethality0 = optim_par_transformed["mr_lethality0"],
    mr_lethality1 = optim_par_transformed["mr_lethality1"],
    pmis_O2_1     = 10^optim_par_transformed["log10_pmis_O2_1"],
    pmis_O2_0     = 10^optim_par_transformed["log10_pmis_O2_0"],
    eta           = 10^optim_par_transformed["log10_eta"]
  )
  
  # Run simulation with best-so-far parameters
  sim_data <- run_all_sims(optim_params_real)
  
  # Generate plot objects
  p_comparison <- generate_comparison_plot(sim_data, x, sim_configs, N_UNIT, N_MAX)
  p_passage    <- generate_passage_plot(sim_data, m_data, sim_configs, POP_GROWTH_FACTOR)
  
  # Save plots
  plot_file_comp <- paste0("plots/tmp/comparison_gen_", current_generation_chunk, ".png")
  plot_file_pass <- paste0("plots/tmp/passage_time_gen_", current_generation_chunk, ".png")
  
  ggsave(filename = plot_file_comp, plot = p_comparison, width = 10, height = 12)
  ggsave(filename = plot_file_pass, plot = p_passage, width = 10, height = 6)
  
  cat("--- Plots saved to", plot_file_comp, "and", plot_file_pass, "---\n")
}


print("--- Optimization Complete ---")

# --- Process and Save Final Results ---
print(optim_result)
print(paste("Best cost:", optim_result$optim$bestval))

# Get best parameters (already have from last loop)
print("Optimized Real-Scale Parameters:")
print(as.data.frame(optim_params_real))

# Save final results for visualization
save(optim_result, optim_params_real, file = "optim_results.RData")
print("--- Final optimization results saved to optim_results.RData ---")