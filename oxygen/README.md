**WGD–Missegregation Model and Parameter Fitting**

This repository contains the implementation of a continuous-time model of whole-genome doubling (WGD) and chromosome missegregation, along with the experimental datasets, parameter-fitting workflow, and plotting utilities. A full mathematical specification of the model is provided in the file model_description.tex.

**Repository Structure**

- code/  
  - model_functions.R – core model (growth function, missegregation kernels, WGD transitions, generator construction)  
  - main_script.R – optimization loop using DEoptim and diagnostic plotting  

- data/  
  - ploidy_distribution.csv – observed ploidy distributions  
  - fit_g.Rds – observed passage-growth rates  

- plots/  
  - tmp/ – diagnostic plots created during optimization  

- model_description.tex – detailed mathematical description  
- README.md – this file  

**Requirements**

The analysis uses R with the following packages:  
Matrix, dplyr, tidyr, ggplot2, data.table, DEoptim

Install them using the standard R function install.packages("pkgname").

**Running the Optimization**

1. Open an R session in the repository root.  
2. Load the model and main script:

   source("code/model_functions.R")  
   source("code/main_script.R")

3. The script will automatically:  
   - load the experimental datasets  
   - construct the WGD–missegregation generator  
   - run DEoptim in iterative batches  
   - save diagnostic plots to plots/tmp  
   - save the best-fit parameters into optim_results.RData

**Reproducing Model Predictions**

After optimization completes, you can regenerate the simulated ploidy trajectories via:

   load("optim_results.RData")  
   sim <- run_all_sims(optim_params_real)

This returns ploidy distributions and passage-growth curves for all experimental conditions.

**Documentation**

The full mathematical specification of the model is provided in model_description.tex, including the growth-rate function, missegregation distribution, generator construction, and details of the optimization procedure.



