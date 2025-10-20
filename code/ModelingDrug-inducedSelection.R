# Title: R script to plot and model ploidy-dependent drug response curves
# Author: Gemini
# Date: 2023-10-27
# Description: This script reads drug response data, fits a ploidy-dependent
#              Hill function model, calculates key parameters (beta, gamma),
#              visualizes the results, and saves parameters to an Excel file.

# --- 0. Load Libraries ---
# Used for string manipulation and writing to Excel files.
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(stringr)
library(openxlsx)

# --- Main Analysis Function ---
analyze_drug_response <- function(file_path, ploidy_map, output_xlsx_path = "fitted_parameters.xlsx", Rmax_optional = NULL) {
  
  # --- 1. Load Data ---
  cat(paste("Reading data from", file_path, "...\n"))
  if (!file.exists(file_path)) {
    stop("Error: Input file not found at the specified path.")
  }
  data <- read.delim(file_path, header = TRUE, sep = "\t", check.names = FALSE)
  data <- standardize_colnames(data)
  
  
  # --- 2. Pre-process Data ---
  # Group columns by ploidy level and average the replicates.
  concentration_col <- colnames(data)[1]
  response_cols <- colnames(data)[2:ncol(data)]
  
  # Use the provided ploidy_map to assign ploidy levels to columns
  if (!all(response_cols %in% names(ploidy_map))) {
    print(setdiff(response_cols, names(ploidy_map)))
    stop("Error: Not all response column names are present in the ploidy_map.")
  }
  ploidy_levels <- ploidy_map[response_cols] # Map column names to ploidy values
  unique_ploidies <- unique(ploidy_levels)
  
  processed_data <- data.frame(Dose = data[[concentration_col]])
  
  cat("Processing and averaging data for each ploidy level...\n")
  for (ploidy in unique_ploidies) {
    # Find columns corresponding to the current ploidy level
    cols_for_ploidy <- response_cols[which(ploidy_levels == ploidy)]
    if (length(cols_for_ploidy) > 1) {
      avg_response <- rowMeans(data[, cols_for_ploidy, drop = FALSE], na.rm = TRUE)
    } else {
      avg_response <- data[, cols_for_ploidy]
    }
    # Name the column in processed_data with the numeric ploidy
    processed_data[[as.character(ploidy)]] <- avg_response
  }
  
  
  # --- 3. Define the Model and Fit ---
  hill_eq <- function(D, R_max, EC50, h) {
    R_max * (1 - (D^h) / (EC50^h + D^h))
  }
  
  fitted_params <- list()
  
  cat("Fitting non-linear model for each ploidy level...\n")
  for (ploidy in unique_ploidies) {
    response_vec <- processed_data[[as.character(ploidy)]]
    dose_vec <- processed_data$Dose
    
    # Starting parameter estimation
    R_max_start <- response_vec[1]
    half_max_response <- R_max_start / 2
    ec50_start_index <- which.min(abs(response_vec - half_max_response))
    ec50_start <- dose_vec[ec50_start_index]
    
    tryCatch({
      
      if (is.null(Rmax_optional)) {
        # --- Fit R_max, EC50, and h ---
        fit <- nls(
          response_vec ~ hill_eq(dose_vec, R_max, EC50, h),
          start = list(R_max = R_max_start, EC50 = ec50_start, h = 1.5),
          control = nls.control(maxiter = 100, warnOnly = TRUE)
        )
        fitted_params[[as.character(ploidy)]] <- coef(fit)
        
      } else {
        # --- R_max is fixed, only fit EC50 and h ---
        hill_eq_fixed_rmax <- function(D, EC50, h) {
          Rmax_optional * (1 - (D^h) / (EC50^h + D^h))
        }
        
        fit <- nls(
          response_vec ~ hill_eq_fixed_rmax(dose_vec, EC50, h),
          start = list(EC50 = ec50_start, h = 1.5),
          control = nls.control(maxiter = 100, warnOnly = TRUE)
        )
        
        params <- coef(fit)
        params['R_max'] <- Rmax_optional
        fitted_params[[as.character(ploidy)]] <- params
      }
      
      cat(paste("Successfully fitted model for", ploidy, "N\n"))
    }, error = function(e) {
      cat(paste("Failed to fit model for", ploidy, "N:", e$message, "\n"))
    })
  }
  
  if (length(fitted_params) < 2) {
    stop("Could not fit the model to at least two ploidy levels. Cannot calculate beta and gamma.")
  }
  
  # --- 4. Calculate Ploidy-Dependent Parameters (beta, gamma) ---
  cat("Calculating ploidy-dependent parameters (beta, gamma)...\n")
  # Assuming the first two unique ploidies are the ones to compare
  P_1 <- unique_ploidies[1]
  P_2 <- unique_ploidies[2]
  
  params_1 <- fitted_params[[as.character(P_1)]]
  params_2 <- fitted_params[[as.character(P_2)]]
  
  # Determine which is the 2N baseline
  if(abs(P_1- 2)<0.05 ){
    EC50_2N <- params_1['EC50']
    h_2N <- params_1['h']
    Rmax_2N <- params_1['R_max']
  } else if (abs(P_2- 2)<0.05) {
    EC50_2N <- params_2['EC50']
    h_2N <- params_2['h']
    Rmax_2N <- params_2['R_max']
  } else {
    # Fallback if no 2N is found, use the lower ploidy as baseline
    cat("Warning: No 2N ploidy found. Using lower ploidy as baseline for EC50_2N and h_2N.\n")
    baseline_idx <- which.min(c(P_1, P_2))
    baseline_ploidy <- unique_ploidies[baseline_idx]
    EC50_2N <- fitted_params[[as.character(baseline_ploidy)]]['EC50']
    h_2N <- fitted_params[[as.character(baseline_ploidy)]]['h']
    Rmax_2N <- fitted_params[[as.character(baseline_ploidy)]]['R_max']
  }
  
  beta <- (params_2['EC50'] - params_1['EC50']) / (P_2 - P_1)
  gamma <- (params_2['h'] - params_1['h']) / (P_2 - P_1)
  
  
  # --- 5. Generate Combined Plot with Fitted Curves ---
  plot_dir <- "plots"
  if (!dir.exists(plot_dir)) dir.create(plot_dir)
  
  plot_filename <- paste0(plot_dir, "/", tools::file_path_sans_ext(basename(file_path)), "_fit.png")
  
  png(plot_filename, width = 1200, height = 800, res = 100)
  cat(paste("Generating plot:", plot_filename, "\n"))
  
  colors <- c("blue", "red", "darkgreen", "purple")
  y_range <- range(data[, response_cols], na.rm = TRUE)
  x_range <- range(data[[concentration_col]][data[[concentration_col]] > 0], na.rm = TRUE)
  
  plot(NULL, xlim = x_range, ylim = y_range, log = "x",
       main = "Ploidy-Dependent Drug Response with Model Fit",
       xlab = concentration_col, ylab = "Cell Viability/Response",
       cex.lab = 1.2, cex.main = 1.4, yaxt = "n")
  
  axis_ticks <- pretty(y_range)
  axis(2, at = axis_ticks, labels = format(axis_ticks, scientific = TRUE), las = 1)
  grid()
  
  for (i in 1:length(response_cols)) {
    points(data[[concentration_col]], data[[response_cols[i]]],
           pch = 19, col = adjustcolor(colors[i], alpha.f = 0.5))
  }
  
  dose_curve <- exp(seq(log(min(x_range)), log(max(x_range)), length.out = 200))
  
  legend_names_fit <- c()
  legend_cols_fit <- c()
  for (i in 1:length(unique_ploidies)) {
    ploidy <- unique_ploidies[i]
    params <- fitted_params[[as.character(ploidy)]]
    if (!is.null(params)) {
      predicted_response <- hill_eq(dose_curve, params['R_max'], params['EC50'], params['h'])
      lines(dose_curve, predicted_response, col = colors[i], lwd = 3)
      legend_names_fit <- c(legend_names_fit, paste0(ploidy, "N Fit")) # Modified legend
      legend_cols_fit <- c(legend_cols_fit, colors[i])
    }
  }
  
  legend("bottomleft", legend = response_cols, col = colors[1:length(response_cols)],
         pch = 19, bty = "n", cex = 1.0, title = "Raw Data", inset = c(0.02, 0.1))
  legend("bottomleft", legend = legend_names_fit, col = legend_cols_fit,
         lwd = 3, bty = "n", cex = 1.0, title = "Model Fit", inset = c(0.02, 0.02))
  
  dev.off()
  
  
  # --- 6. Output Final Parameters to Console and XLSX ---
  results_df <- data.frame(
    SourceFile = basename(file_path),
    Rmax_2N = Rmax_2N,
    EC50_2N = EC50_2N,
    h_2N = h_2N,
    Beta = beta,
    Gamma = gamma,
    Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  cat("\n--- Model Fitting Summary ---\n")
  print(results_df)
  cat("---------------------------\n")
  
  # --- Write to Excel ---
  cat(paste("Writing parameters to", output_xlsx_path, "...\n"))
  sheet_name <- tools::file_path_sans_ext(basename(file_path))
  
  if (file.exists(output_xlsx_path)) {
    wb <- loadWorkbook(output_xlsx_path)
  } else {
    wb <- createWorkbook()
  }
  
  if (!sheet_name %in% names(wb)) {
    addWorksheet(wb, sheet_name)
  }
  
  writeData(wb, sheet = sheet_name, x = results_df, startCol = 1, startRow = 1)
  saveWorkbook(wb, output_xlsx_path, overwrite = TRUE)
  
  cat("Analysis complete.\n\n")
}



standardize_colnames <- function(tab) {
  # Normalize SUM159 spacing before pattern matching
  colnames(tab) <- gsub("159(", "159 (", colnames(tab), fixed = TRUE)
  colnames(tab) <- gsub("SUM159", "SUM-159", colnames(tab), fixed = TRUE)
  
  # Standardize colnames to match naming in rownames(ploidy)
  colnames(tab) <- gsub("MDA[- ]?MB[- ]?231", "MDAMB231", colnames(tab), ignore.case = TRUE)
  colnames(tab) <- gsub("SUM[- ]?159", "SUM-159", colnames(tab), ignore.case = TRUE)
  colnames(tab) <- gsub("HS?578[ -]?T", "HS578T", colnames(tab), ignore.case = TRUE)
  colnames(tab) <- gsub("MCF10[ -]?DCIS", "MCFdcis", colnames(tab), ignore.case = TRUE)
  colnames(tab) <- gsub("SUM[- ]?149", "SUM-149", colnames(tab), ignore.case = TRUE)
  
  # Trim whitespace
  colnames(tab) <- trimws(colnames(tab))
  
  return(tab)
}
