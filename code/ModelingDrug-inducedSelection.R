# Title: R script for ROBUST GLOBAL fitting of ploidy-dependent drug response using nls
# Author: Gemini
# Date: 2025-10-20
# Description: This script normalizes data to have a global R_max of 1. It then
#              uses a robust hybrid nls approach to fit a global model.

# --- 0. Load Libraries ---
# (Libraries remain the same)
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(tidyr)
library(dplyr)
library(openxlsx)

analyze_drug_response_global <- function(file_path, ploidy_map, output_xlsx_path = "fitted_parameters_global.xlsx") {
  
  # --- 1. Load and Normalize Data ---
  cat(paste("Reading data from", file_path, "...\n"))
  data_raw <- read.delim(file_path, header = TRUE, sep = "\t", check.names = FALSE)
  data_raw <- standardize_colnames(data_raw)
  
  concentration_col <- colnames(data_raw)[1]
  response_cols <- colnames(data_raw)[2:ncol(data_raw)]
  
  # Keep the original concentration column
  concentrations <- data_raw[, concentration_col, drop = FALSE]
  
  # Normalize the response columns so that R_max is 1 for all
  cat("Normalizing response data to force global R_max = 1...\n")
  response_normalized <- sweep(data_raw[, response_cols], MARGIN = 2, FUN = "/", apply(data_raw[, response_cols], 2, max))
  colnames(response_normalized) = response_cols
  
  # Recombine the data
  data <- cbind(concentrations, response_normalized)
  
  # --- 2. Reshape Data to Long Format ---
  if (!all(response_cols %in% names(ploidy_map))) {
    stop("Error: Not all response column names are present in the ploidy_map.")
  }
  long_data <- data %>%
    pivot_longer(cols = all_of(response_cols), names_to = "Sample", values_to = "Response") %>%
    mutate(Ploidy = ploidy_map[Sample]) %>%
    dplyr::rename(Dose = all_of(concentration_col)) %>%
    dplyr::select(Dose, Ploidy, Response) %>%
    na.omit()
  
  # NOTE: The ploidy-specific r_max_estimates section has been removed.
  
  # --- 3. Pre-fitting Step using nls ---
  cat("Performing independent nls fits to find robust starting parameters...\n")
  independent_params <- list()
  unique_ploidies <- sort(unique(long_data$Ploidy))
  
  # This equation now assumes R_max is always 1
  hill_eq_norm <- function(D, EC50, h) {
    1 * (1 - (D^h) / (EC50^h + D^h))
  }
  
  for (ploidy in unique_ploidies) {
    subset_data <- filter(long_data, Ploidy == ploidy & Dose > 0)
    if (nrow(subset_data) == 0) next
    
    # Starting EC50 is where the response is closest to 0.5
    ec50_start_idx <- which.min(abs(subset_data$Response - 0.5))
    ec50_start <- if(length(ec50_start_idx) > 0) subset_data$Dose[ec50_start_idx] else median(subset_data$Dose)
    
    tryCatch({
      fit <- nls(
        Response ~ hill_eq_norm(Dose, EC50, h),
        data = subset_data,
        start = list(EC50 = ec50_start, h = 1.5),
        control = nls.control(maxiter = 100, warnOnly = TRUE)
      )
      independent_params[[as.character(ploidy)]] <- coef(fit)
    }, error = function(e) {
      cat(paste("Warning: Pre-fit with nls failed for ploidy", ploidy, "\n"))
    })
  }
  
  if (length(independent_params) < 2) {
    stop("Could not generate stable pre-fits for at least two ploidy levels. Cannot proceed.")
  }
  
  param_df_indep <- do.call(rbind, lapply(names(independent_params), function(p) {
    data.frame(Ploidy = as.numeric(p), EC50 = unname(independent_params[[p]]['EC50']), h = unname(independent_params[[p]]['h']))
  }))
  param_df_indep <- filter(param_df_indep, is.finite(EC50) & is.finite(h))
  
  if (nrow(param_df_indep) < 2) {
    stop("Fewer than two valid pre-fits were generated after cleaning. Cannot calculate beta/gamma slopes.")
  }
  
  fit_beta_start <- lm(EC50 ~ I(Ploidy - 2), data = param_df_indep)
  fit_gamma_start <- lm(h ~ I(Ploidy - 2), data = param_df_indep)
  
  start_params <- list(
    EC50_2N = coef(fit_beta_start)['(Intercept)'],
    beta = coef(fit_beta_start)['I(Ploidy - 2)'],
    h_2N = coef(fit_gamma_start)['(Intercept)'],
    gamma = coef(fit_gamma_start)['I(Ploidy - 2)']
  )
  start_params[is.na(start_params)] <- 0
  
  # --- 4. Define and Fit the Global Model ---
  cat("Fitting global non-linear model with R_max = 1...\n")
  fit_result <- NULL
  
  fitting_data <- filter(long_data, Dose > 0)
  
  tryCatch({
    # The model formula now has R_max hard-coded as 1
    fit_result <- nls(
      Response ~ 1 * (1 - Dose^(h_2N + gamma * (Ploidy - 2)) / 
                        ((EC50_2N + beta * (Ploidy - 2))^(h_2N + gamma * (Ploidy - 2)) + Dose^(h_2N + gamma * (Ploidy - 2)))),
      data = fitting_data,
      start = start_params,
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    )
    cat("Successfully fitted the global model.\n")
    fitted_params <- coef(fit_result)
    if(is.null(dim(fitted_params))){
      names(fitted_params) = sapply(strsplit(names(fitted_params),".", fixed = T),"[[",1)
    }
    fit_method <- "Global Fit"
  }, error = function(e) {
    cat(paste("Warning: Global model fit failed:", e$message, "\n"))
    cat(">>> FALLING BACK to parameters derived from independent fits.\n")
    fitted_params <- c(EC50_2N = start_params$EC50_2N, h_2N = start_params$h_2N, beta = start_params$beta, gamma = start_params$gamma)
    fit_method <- "Independent Fits (Fallback)"
  })
  
  # --- 5. Generate Conditional Plots ---
  # ... (Plotting code remains the same, but will now use R_max = 1 implicitly) ...
  plot_dir <- "plots"
  if (!dir.exists(plot_dir)) dir.create(plot_dir)
  plot_filename <- paste0(plot_dir, "/", tools::file_path_sans_ext(basename(file_path)), "_fit_global.png")
  png(plot_filename, width = 1600, height = 1200, res = 120)
  cat(paste("Generating plot:", plot_filename, "\n"))
  
  if (fit_method == "Global Fit") { layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)) } else { layout(1) }
  par(mar = c(5, 5, 4, 2))
  
  # -- PLOT 1: Main Dose-Response Curve --
  x_range <- range(long_data$Dose[long_data$Dose > 0], na.rm = TRUE)
  y_range <- range(long_data$Response, na.rm = TRUE)
  plot_title <- paste("Ploidy-Dependent Drug Response (", fit_method, ")", sep="")
  plot(NULL, xlim = x_range, ylim = y_range, log = "x", main = plot_title,
       xlab = concentration_col, ylab = "Normalized Response", cex.main=1.5, cex.lab=1.2)
  grid()
  
  num_ploidies <- length(unique_ploidies)
  plot_colors <- hcl.colors(num_ploidies, palette = "viridis")
  color_map <- setNames(plot_colors, unique_ploidies)
  
  for (ploidy in unique_ploidies) {
    subset_data <- filter(long_data, Ploidy == ploidy)
    points(subset_data$Dose, subset_data$Response, pch = 19, col = adjustcolor(color_map[as.character(ploidy)], alpha.f = 0.4))
  }
  
  dose_curve <- exp(seq(log(min(x_range)), log(max(x_range)), length.out = 200))
  
  for (ploidy in unique_ploidies) {
    EC50_P <- fitted_params["EC50_2N"] + fitted_params["beta"] * (ploidy - 2)
    h_P <- fitted_params["h_2N"] + fitted_params["gamma"] * (ploidy - 2)
    # R_max is now always 1
    predicted_response <- 1 * (1 - dose_curve^h_P / (EC50_P^h_P + dose_curve^h_P))
    lines(dose_curve, predicted_response, col = color_map[as.character(ploidy)], lwd = 3)
  }
  
  legend_labels <- sapply(unique_ploidies, function(p) {
    sample_name <- names(ploidy_map)[match(p, ploidy_map)]
    rounded_ploidy <- format(round(p, 1), nsmall = 1)
    paste0(sample_name, " (", rounded_ploidy, "N)")
  })
  legend("bottomleft", legend = legend_labels, col = unname(color_map), pch = 19, lwd = 3, bty = "n", cex = 1.0, title = "Condition")
  
  # -- PLOTS 2 & 3: Regression Plots --
  if (fit_method == "Global Fit") {
    # (This section is unchanged)
    param_df_global <- data.frame(Ploidy = unique_ploidies)
    param_df_global$EC50 <- fitted_params["EC50_2N"] + fitted_params["beta"] * (param_df_global$Ploidy - 2)
    param_df_global$h <- fitted_params["h_2N"] + fitted_params["gamma"] * (param_df_global$Ploidy - 2)
    plot(param_df_global$Ploidy, param_df_global$EC50, main = "EC50 vs. Ploidy (from Global Fit)", xlab = "Ploidy (N)", ylab = "EC50", pch = 19, col = "dodgerblue", cex = 1.5, cex.main=1.5, cex.lab=1.2)
    grid(); abline(a = fitted_params["EC50_2N"] - 2 * fitted_params["beta"], b = fitted_params["beta"], col = "firebrick", lwd = 2)
    legend("topleft", legend = paste0("Beta = ", format(fitted_params["beta"], digits = 3)), bty = "n", cex = 1.2)
    plot(param_df_global$Ploidy, param_df_global$h, main = "Hill Slope vs. Ploidy (from Global Fit)", xlab = "Ploidy (N)", ylab = "Hill Slope (h)", pch = 19, col = "dodgerblue", cex = 1.5, cex.main=1.5, cex.lab=1.2)
    grid(); abline(a = fitted_params["h_2N"] - 2 * fitted_params["gamma"], b = fitted_params["gamma"], col = "firebrick", lwd = 2)
    legend("topleft", legend = paste0("Gamma = ", format(fitted_params["gamma"], digits = 3)), bty = "n", cex = 1.2)
  }
  dev.off()
  
  # --- 6. Output Final Parameters ---
  # (This section is unchanged)
  results_df <- data.frame(
    SourceFile = basename(file_path),
    EC50_2N = fitted_params["EC50_2N"],
    h_2N = fitted_params["h_2N"],
    Beta = fitted_params["beta"],
    Gamma = fitted_params["gamma"],
    FitMethod = fit_method,
    Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  cat("\n--- Model Fitting Summary ---\n"); print(results_df); cat("------------------------------------\n")
  
  # (Code to save to XLSX would follow here)
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
