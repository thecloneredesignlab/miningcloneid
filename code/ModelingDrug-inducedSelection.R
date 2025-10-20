
# Title: R script for ROBUST GLOBAL fitting of ploidy-dependent drug response using nls
# Author: Gemini
# Date: 2025-10-20
# Description: This script uses a robust hybrid approach with nls. It first performs
#              independent nls fits to find intelligent starting parameters, then
#              attempts a global nls model fit. If the global fit fails, it
#              falls back to the parameters derived from the independent fits.
#              The output includes a multi-panel plot with regression plots
#              only if the global fit is successful. The plot legend now shows
#              the sample name and rounded ploidy value.

# --- 0. Load Libraries ---
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(tidyr)
library(dplyr)
library(openxlsx)

analyze_drug_response_global <- function(file_path, ploidy_map, output_xlsx_path = "fitted_parameters_global.xlsx") {
  condition=fileparts(file_path)$name
  
  # --- 1. Load Data ---
  cat(paste("Reading data from", file_path, "...\n"))
  data <- read.delim(file_path, header = TRUE, sep = "\t", check.names = FALSE)
  data <- standardize_colnames(data)
  
  concentration_col <- colnames(data)[1]
  response_cols <- colnames(data)[2:ncol(data)]
  
  # Normalize the response columns so that R_max is 1 for all
  cat("Normalizing response data to force global R_max = 1...\n")
  data[,2:ncol(data)] <- sweep(data[, response_cols], MARGIN = 2, FUN = "/", apply(data[, response_cols], 2, max))
  
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
  
  r_max_estimates <- long_data %>%
    filter(Dose == 0) %>%
    group_by(Ploidy) %>%
    summarise(R_max = mean(Response, na.rm = TRUE), .groups = 'drop')
  
  long_data <- left_join(long_data, r_max_estimates, by = "Ploidy")
  
  # --- 3. Pre-fitting Step using nls to Get Robust Starting Values ---
  cat("Performing independent nls fits to find robust starting parameters...\n")
  independent_params <- list()
  unique_ploidies <- sort(unique(long_data$Ploidy))
  
  hill_eq_simple <- function(D, R_max, EC50, h) {
    R_max * (1 - (D^h) / (EC50^h + D^h))
  }
  
  for (ploidy in unique_ploidies) {
    subset_data <- filter(long_data, Ploidy == ploidy & Dose > 0)
    if (nrow(subset_data) == 0) next
    
    r_max_start <- r_max_estimates$R_max[r_max_estimates$Ploidy == ploidy]
    half_max_resp <- r_max_start / 2
    ec50_start_idx <- which.min(abs(subset_data$Response - half_max_resp))
    ec50_start <- if(length(ec50_start_idx) > 0) subset_data$Dose[ec50_start_idx] else median(subset_data$Dose)
    
    tryCatch({
      # --- First Attempt: nls ---
      fit <- nls(
        Response ~ hill_eq_simple(Dose, R_max, EC50, h),
        data = subset_data,
        start = list(R_max = r_max_start, EC50 = ec50_start, h = 1.5),
        control = nls.control(maxiter = 100, warnOnly = TRUE)
      )
      # Use <<- for consistency, though <- often works in the main 'try' block
      independent_params[[as.character(ploidy)]] <- coef(fit)
      
    }, error = function(e) {
      # --- Fallback: drm ---
      cat(paste("Warning: nls pre-fit failed for ploidy", ploidy, ". Trying drm fallback...\n"))
      
      tryCatch({
        # Use the three-parameter log-logistic model from drc
        fit_drm <- drc::drm(Response ~ Dose, data = subset_data, fct = drc::LL.3())
        
        # Extract and rename coefficients
        drm_coefs <- coef(fit_drm)
        params <- c(
          'R_max' = as.numeric(drm_coefs['d:(Intercept)']),
          'EC50'  = as.numeric(drm_coefs['e:(Intercept)']),
          'h'     = as.numeric(-drm_coefs['b:(Intercept)'])
        )
        
        # # FIX: Use the superassignment operator to modify the parent variable
        independent_params[[as.character(ploidy)]] <<- params
        
        cat(paste("Success: drm fallback fit succeeded for ploidy", ploidy, "\n"))
        
      }, error = function(e_drm) {
        # This runs if drm also fails
        cat(paste("Error: drm fallback also failed for ploidy", ploidy, ":", e_drm$message, "\n"))
      })
    })
  }
  if (length(independent_params) < 2) {
    stop("Could not generate stable pre-fits for at least two ploidy levels. Cannot proceed.")
  }
  
  # Get the dose range from the data, excluding zero concentrations
  dose_range <- range(long_data$Dose[long_data$Dose > 0], na.rm = TRUE)
  param_df_indep <- do.call(rbind, lapply(names(independent_params), function(p) {
    ploidy_val <- as.numeric(p)
    params <- independent_params[[p]]
    # Get the parameters for this ploidy from the pre-fit
    ec50_val <- unname(params['EC50'])
    h_val <- unname(params['h'])
    # Get the corresponding R_max estimated from the raw data
    r_max_val <- r_max_estimates$R_max[r_max_estimates$Ploidy == ploidy_val]
    # Calculate AUC using the parameters from this independent fit
    auc_val <- NA # Default to NA
    if (all(is.finite(c(ec50_val, h_val, r_max_val)))) {
      auc_val <- calculate_log_auc(ec50_val, h_val, r_max_val, dose_range[1], dose_range[2])
    }
    # Create the data frame row with the new AUC column
    data.frame(
      lineage = names(ploidy_map)[ploidy_map == ploidy_val][1], # [1] handles duplicate ploidies
      Ploidy = ploidy_val,
      EC50 = ec50_val,
      h = h_val,
      AUC = auc_val
    )
  }))
  
  param_df_indep <- param_df_indep %>%
    filter(is.finite(EC50) & is.finite(h))
  
  if (nrow(param_df_indep) < 2) {
    stop("Fewer than two valid pre-fits were generated after cleaning. Cannot calculate beta/gamma slopes.")
  }
  
  # Load workbook or create a new one
  if (file.exists(output_xlsx_path)) {
    wb <- loadWorkbook(output_xlsx_path)
  } else {
    wb <- createWorkbook()
  }
  if (!condition %in% names(wb)) { addWorksheet(wb, condition) }
  writeData(wb, sheet = condition, x = param_df_indep)
  saveWorkbook(wb, output_xlsx_path, overwrite = TRUE)
  
  
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
  cat("Fitting global non-linear model with intelligent starting values...\n")
  fit_result <- NULL
  
  fitting_data <- filter(long_data, Dose > 0)
  
  tryCatch({
    fit_result <- nls(
      Response ~ R_max * (1 - Dose^(h_2N + gamma * (Ploidy - 2)) / 
                            ((EC50_2N + beta * (Ploidy - 2))^(h_2N + gamma * (Ploidy - 2)) + Dose^(h_2N + gamma * (Ploidy - 2)))),
      data = fitting_data,
      start = start_params,
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    )
    cat("Successfully fitted the global model.\n")
    fitted_params <- coef(fit_result)
    fit_method <- "Global Fit"
  }, error = function(e) {
    cat(paste("Warning: Global model fit failed:", e$message, "\n"))
    cat(">>> FALLING BACK to parameters derived from independent fits.\n")
    fitted_params <- c(EC50_2N = start_params$EC50_2N, h_2N = start_params$h_2N, beta = start_params$beta, gamma = start_params$gamma)
    fit_method <- "Independent Fits (Fallback)"
  })
  
  # --- 5. [CORRECTED] Generate Multi-Panel Plots ---
  plot_dir <- "plots"; if (!dir.exists(plot_dir)) dir.create(plot_dir)
  plot_filename <- paste0(plot_dir, "/", condition, "_fit_comparison.png")
  png(plot_filename, width = 1800, height = 1600, res = 120)
  cat(paste("Generating plot:", plot_filename, "\n"))
  
  if (fit_method == "Global Fit") { layout(matrix(1:4, 2, 2, byrow = TRUE)) } else { par(mfrow = c(1, 2)) }
  par(mar = c(5, 5, 4, 2))
  
  num_ploidies <- length(unique_ploidies); plot_colors <- hcl.colors(num_ploidies, palette = "viridis")
  color_map <- setNames(plot_colors, unique_ploidies); x_range <- range(long_data$Dose[long_data$Dose > 0], na.rm = TRUE)
  y_range <- range(long_data$Response, na.rm = TRUE)
  dose_curve <- exp(seq(log(min(x_range)), log(max(x_range)), length.out = 200))
  legend_labels <- sapply(unique_ploidies, function(p) {
    sample_name <- names(ploidy_map)[match(p, ploidy_map)]; rounded_ploidy <- format(round(p, 1), nsmall = 1)
    paste0(sample_name, " (", rounded_ploidy, "N)")
  })
  
  # -- PLOT 1 (LEFT): Primary Model Fit --
  plot_title_left <- paste("Primary Model (", fit_method, ")\n", condition, sep="")
  plot(NULL, xlim = x_range, ylim = y_range, log = "x", main = plot_title_left, xlab = concentration_col, ylab = "Response", yaxt = "n", cex.main=1.5, cex.lab=1.2)
  axis(2, at = pretty(y_range), labels = format(pretty(y_range), scientific = TRUE), las = 1); grid()
  for (ploidy in unique_ploidies) {
    # FIX: Use formula notation `Response ~ Dose`
    points(Response ~ Dose, data = filter(long_data, Ploidy == ploidy), pch = 19, col = adjustcolor(color_map[as.character(ploidy)], alpha.f = 0.4))
  }
  for (ploidy in unique_ploidies) {
    EC50_P <- fitted_params["EC50_2N"] + fitted_params["beta"] * (ploidy - 2); h_P <- fitted_params["h_2N"] + fitted_params["gamma"] * (ploidy - 2)
    R_max_P <- r_max_estimates$R_max[r_max_estimates$Ploidy == ploidy]
    if(!is.na(R_max_P)) {
      predicted_response <- R_max_P * (1 - dose_curve^h_P / (EC50_P^h_P + dose_curve^h_P))
      lines(dose_curve, predicted_response, col = color_map[as.character(ploidy)], lwd = 3, lty = 1)
    }
  }
  legend("bottomleft", legend = legend_labels, col = unname(color_map), pch = 19, lwd = 3, bty = "n", cex = 1.0, title = "Condition")
  
  # -- PLOT 2 (RIGHT): Independent Pre-Fits --
  plot(NULL, xlim = x_range, ylim = y_range, log = "x", main = "Independent Pre-Fits", xlab = concentration_col, ylab = "Response", yaxt = "n", cex.main=1.5, cex.lab=1.2)
  axis(2, at = pretty(y_range), labels = format(pretty(y_range), scientific = TRUE), las = 1); grid()
  for (ploidy in unique_ploidies) {
    # FIX: Use formula notation `Response ~ Dose`
    points(Response ~ Dose, data = filter(long_data, Ploidy == ploidy), pch = 19, col = adjustcolor(color_map[as.character(ploidy)], alpha.f = 0.4))
  }
  for (p_str in names(independent_params)) {
    params <- independent_params[[p_str]]; R_max_P <- params['R_max']; EC50_P <- params['EC50']; h_P <- params['h']
    if(all(!is.na(c(R_max_P, EC50_P, h_P)))) {
      predicted_response <- R_max_P * (1 - dose_curve^h_P / (EC50_P^h_P + dose_curve^h_P))
      lines(dose_curve, predicted_response, col = color_map[p_str], lwd = 3, lty = 2)
    }
  }
  legend("bottomleft", legend = legend_labels, col = unname(color_map), pch = 19, lwd = 3, lty=2, bty = "n", cex = 1.0, title = "Condition")
  
  # -- PLOTS 3 & 4: Regression Plots (Only if Global Fit Succeeded) --
  if (fit_method == "Global Fit") {
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
  results_df <- data.frame(
    SourceFile = basename(file_path),
    EC50_2N = fitted_params["EC50_2N"],
    h_2N = fitted_params["h_2N"],
    Beta = fitted_params["beta"],
    Gamma = fitted_params["gamma"],
    FitMethod = fit_method,
    Timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  cat("\n--- Model Fitting Summary ---\n")
  print(results_df)
  cat("------------------------------------\n")
  
  # (Code to save to XLSX would follow here)
  return(fit_result)
}

standardize_colnames <- function(tab) {
  # Normalize SUM159 spacing before pattern matching
  colnames(tab) <- gsub("159(", "159 (", colnames(tab), fixed = TRUE)
  colnames(tab) <- gsub("SUM159", "SUM-159", colnames(tab), fixed = TRUE)
  
  colnames(tab) <- gsub("MDA-MB-231(", "MDA-MB-231 (", colnames(tab), fixed = TRUE)
  
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

# Title: R function to compare an arbitrary number of ploidy-dependent regression fits
# Author: Gemini
# Date: 2025-10-20
# Description: This function takes a named list of model summary outputs and
#              plots their EC50 vs. Ploidy and Hill Slope vs. Ploidy regression
#              lines on the same graphs for direct comparison. It dynamically
#              handles any number of models, using the user's specified
#              bracket notation for parameter access.

compare_model_fits <- function(model_list, output_filename = "regression_comparison.png",legloc="topleft",ploidy_range = seq(2, 8, by = 0.1)) {
  
  cat(paste("Generating comparison plot for", length(model_list), "models:", output_filename, "\n"))
  
  # --- Setup plot aesthetics for multiple models ---
  num_models <- length(model_list)
  if (num_models == 0) {
    stop("Input 'model_list' is empty. Nothing to plot.")
  }
  plot_colors <- hcl.colors(num_models, palette = "viridis")
  plot_ltys <- 1:num_models
  
  # --- Setup the PNG device and layout ---
  png(output_filename, width = 1600, height = 800, res = 120)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), cex.axis = 1, cex.lab = 1.2, cex.main = 1.5)
  
  
  # --- 1. EC50 vs. Ploidy Comparison Plot ---
  
  # Calculate all predictions using bracket notation
  all_ec50_preds <- lapply(model_list, function(model) {
    # Using [[...]] which is safest for lists/dataframes within lapply
    if (all(is.finite(c(model[["EC50_2N"]], model[["beta"]])))) {
      return(model[["EC50_2N"]] + model[["beta"]] * (ploidy_range - 2))
    }
    return(NA)
  })
  
  y_min_ec50 <- min(unlist(all_ec50_preds), na.rm = TRUE)
  y_max_ec50 <- max(unlist(all_ec50_preds), na.rm = TRUE)
  
  if (!is.finite(y_min_ec50)) {
    plot(1, type="n", main="Comparison of EC50 vs. Ploidy", xlab="Ploidy (N)", ylab="EC50")
    text(1, 1, "No valid model data to plot.", cex=1.5)
  } else {
    plot(NULL, xlim = range(ploidy_range), ylim = c(y_min_ec50, y_max_ec50),
         main = "Comparison of EC50 vs. Ploidy", xlab = "Ploidy (N)", ylab = "EC50")
    grid()
    
    for (i in 1:num_models) {
      if (!all(is.na(all_ec50_preds[[i]]))) {
        lines(ploidy_range, all_ec50_preds[[i]], col = plot_colors[i], lty = plot_ltys[i], lwd = 3)
      }
    }
    
    valid_indices <- which(!sapply(all_ec50_preds, function(x) all(is.na(x))))
    legend_texts <- sapply(valid_indices, function(i) {
      paste(names(model_list)[i], "(Beta =", format(model_list[[i]][["beta"]], digits = 2), ")")
    })
    legend(legloc, legend = legend_texts, col = plot_colors[valid_indices],
           lwd = 3, lty = plot_ltys[valid_indices], bty = "n", cex = 1.0)
  }
  
  
  # --- 2. Hill Slope (h) vs. Ploidy Comparison Plot ---
  
  all_h_preds <- lapply(model_list, function(model) {
    if (all(is.finite(c(model[["h_2N"]], model[["gamma"]])))) {
      return(model[["h_2N"]] + model[["gamma"]] * (ploidy_range - 2))
    }
    return(NA)
  })
  
  y_min_h <- min(unlist(all_h_preds), na.rm = TRUE)
  y_max_h <- max(unlist(all_h_preds), na.rm = TRUE)
  
  if (!is.finite(y_min_h)) {
    plot(1, type="n", main="Comparison of Hill Slope vs. Ploidy", xlab="Ploidy (N)", ylab="Hill Slope (h)")
    text(1, 1, "No valid model data to plot.", cex=1.5)
  } else {
    plot(NULL, xlim = range(ploidy_range), ylim = c(y_min_h, y_max_h),
         main = "Comparison of Hill Slope vs. Ploidy", xlab = "Ploidy (N)", ylab = "Hill Slope (h)")
    grid()
    
    for (i in 1:num_models) {
      if (!all(is.na(all_h_preds[[i]]))) {
        lines(ploidy_range, all_h_preds[[i]], col = plot_colors[i], lty = plot_ltys[i], lwd = 3)
      }
    }
    
    valid_indices <- which(!sapply(all_h_preds, function(x) all(is.na(x))))
    legend_texts <- sapply(valid_indices, function(i) {
      paste(names(model_list)[i], "(gamma =", format(model_list[[i]][["gamma"]], digits = 2), ")")
    })
    legend("topleft", legend = legend_texts, col = plot_colors[valid_indices],
           lwd = 3, lty = plot_ltys[valid_indices], bty = "n", cex = 1.0)
  }
  
  dev.off()
  cat("Plot saved successfully.\n")
}

# Helper function for calculating Area Under the Curve
calculate_log_auc <- function(EC50, h, R_max, min_dose, max_dose) {
  # Define the response function based on a log-dose scale
  response_func_log <- function(log_dose) {
    dose <- exp(log_dose)
    R_max * (1 - (dose^h) / (EC50^h + dose^h))
  }
  # Integrate over the log-transformed dose range
  result <- try(integrate(response_func_log, lower = log(min_dose), upper = log(max_dose)), silent = TRUE)
  
  # Return NA if integration fails, otherwise return the value
  if (inherits(result, "try-error")) return(NA) else return(result$value)
}