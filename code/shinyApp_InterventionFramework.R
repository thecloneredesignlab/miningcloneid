library(shiny)
library(deSolve)

# ==============================================================================
# === Python Model Translations (UPDATED) ======================================
# ==============================================================================

# --- Drug Dosing Definitions (UPDATED) ---

drug_dosing_schedules <- list(
    "volasertib"    = "IV",
    "umi-77"        = "IV",
    "tegafur"       = "Oral",
    "tas"           = "Oral",
    "osi-027"       = "Oral",
    "alisertib"     = "IV", # Changed
    "5-azacytidine" = "Oral",
    "abt-199"       = "Oral",
    "abt-263"       = "Oral",
    "capecitabine"  = "Oral",
    "ceralasertib"  = "Oral",
    "cytarabine"    = "IV",
    "gemcitabine"   = "IV",
    "bay1895344"    = "IV", # Changed
    "ispinesib"     = "IV",
    "navitoclax"    = "Oral",
    "adavosertib"   = "Oral"
)

IV_DEFAULTS <- list(C_peak=1.0, half_life=0.5, period=7.0)

ORAL_DEFAULTS <- list(
    dose = 100.0, F = 0.6, Vd = 60.0,
    ka_day = 2.0, ke_day = 0.7,
    period = 1.0, tlag = 0.0
)

# --- PER_DRUG (UPDATED) ---
PER_DRUG <- list(
    # IV examples
    "volasertib"   = list(C_peak = 1.0, half_life = 4.0, period = 7.0),
    "alisertib"    = list(C_peak = 42.5, half_life = 0.875, period = 0.5), # Updated
    "cytarabine"   = list(C_peak = 1.0, half_life = 0.2, period = 3.5),
    "gemcitabine"  = list(C_peak = 239, half_life = 0.05, period = 7.0), # Updated
    "ispinesib"    = list(C_peak = 2.1, half_life = 1.04, period = 7), # New
    "umi-77"       = list(C_peak = 1.0, half_life = 0.8, period = 7.0),
    "navitoclax"   = list(C_peak = 1.0, half_life = 0.73, period = 1), # New
    "bay1895344"   = list(C_peak = 6.2, half_life = 0.50, period = 0.5), # New
    
    # Oral examples
    "osi-027"      = list(dose = 50, F = 0.5, Vd = 80, ka_day = 1.8, ke_day = 0.5, period = 1.0),
    "abt-199"      = list(dose = 100, F = 0.6, Vd = 250, ka_day = 1.0, ke_day = 0.3, period = 1.0),
    "abt-263"      = list(dose = 100, F = 0.6, Vd = 120, ka_day = 2.0, ke_day = 0.5, period = 1.0),
    "ceralasertib" = list(dose = 80,  F = 0.5, Vd = 100,ka_day = 2.0, ke_day = 0.5, period = 1.0),
    "adavosertib"  = list(dose = 100, F = 0.6, Vd = 65, ka_day = 2.4, ke_day = 0.6, period = 1.0),
    "tas"          = list(dose = 60,  F = 0.5, Vd = 40, ka_day = 2.4, ke_day = 0.7, period = 1.0),
    "tegafur"      = list(dose = 40, F = 0.5,Vd = 45, ka_day = 1.6, ke_day = 0.5, period = 1.0),
    "capecitabine" = list(dose = 100, F = 0.8, Vd = 40, ka_day = 3.0, ke_day = 0.6, period = 1.0),
    "5-azacytidine"= list(dose = 100, F = 0.2, Vd = 40, ka_day = 3.0, ke_day = 2.0, period = 1.0)
)


# --- PD Function: phi_Hill (Unchanged) ---
phi_Hill <- function(C, EC50, n, Emax=1.0) {
    # Hill-type kill rate function
    Emax * (C^n) / (EC50^n + C^n)
}

# --- PD Function: f(ploidy, drug) (UPDATED) ---
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
        
    } else {
        # Default fallback for unknown drugs
        warning(paste("Unknown drug:", drug, "- using default parameters (EC50=1.0, n=1.0)."))
        list(EC50=1.0, n=1.0, Emax=1.0)
    }
}

# --- PK Function: pulsed_dose (Unchanged) ---
pulsed_dose <- function(C_peak=5.0, half_life=2.0, period=7.0) {
    lam <- log(2) / max(half_life, 1e-12)
    function(t) {
        t <- as.numeric(t)
        modt <- t %% period
        C_peak * exp(-lam * modt)
    }
}

# --- PK Function: oral_pulsed_ss_days (Unchanged) ---
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

# --- PK Function: get_concentration_curve (Unchanged) ---
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
        
    } else {
        stop(paste("Unsupported route", route))
    }
}


# ==============================================================================
# === Core Shiny App Code (Unchanged) ==========================================
# ==============================================================================

# ---- ODE model (Unchanged) ----
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

# --- run_one_cycle (Unchanged) ----
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

# ---- Parsing helpers (Unchanged) ----
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

# ---- UI (UPDATED) ----
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
                      placeholder = "e.g. 0.6,0.3,0.1"),
            numericInput("B0", "Initial tumor burden (cells)", value = 1e7, min = 1, step = 1e6),
            textInput("drug", "Drug name", value = "gemcitabine"), # Default drug updated
            
            numericInput("cycleDays", "Cycle Length (days)", value = 28, min = 1, max = 100, step = 1),
            
            actionButton("runCycle", "Run Next Cycle", class = "btn-primary"),
            hr(),
            textAreaInput("measurements", "Add tumor-burden measurements (optional)",
                          placeholder = "Format: %burden time\nExample:\n90 5\n80 10\n70 20",
                          height = "120px"),
            actionButton("addMeas", "Add Measurements"),
            actionButton("correctModel", "Correct Model", class = "btn-info"), 
            hr(),
            div(class = "small-note",
                tags$ul(
                    # Note updated with new drugs
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

# ---- Server (Unchanged) ----
server <- function(input, output, session) {
    
    all_cycles <- reactiveVal(list())
    last_B <- reactiveVal(NULL)
    
    current_k_mult <- reactiveVal(1.0)
    
    observeEvent(input$runCycle, {
        fr <- parse_fractions(input$fractions)
        validate(
            need(!is.null(fr), "Please enter valid ploidy fractions (e.g., 0.6,0.3,0.1)."),
            need(input$B0 > 0, "Initial tumor burden must be > 0."),
            need(input$cycleDays > 0, "Cycle length must be > 0.")
        )
        
        B_start <- if (is.null(last_B())) input$B0 else last_B()
        
        k_mult_base <- current_k_mult()
        
        ode_df <- run_one_cycle(ploidy_fracs = fr, 
                                B0 = B_start, 
                                drug_name = input$drug,
                                days = input$cycleDays,
                                k_multiplier_base = k_mult_base)
        
        cycles <- all_cycles()
        idx <- length(cycles) + 1
        
        if (idx > 1) {
            last_time <- max(cycles[[idx - 1]]$df$cum_time)
            ode_df$cum_time <- ode_df$time + last_time
        } else {
            ode_df$cum_time <- ode_df$time
        }
        
        meas_df <- data.frame(percent = numeric(0), time = numeric(0),
                              B = numeric(0), cum_time = numeric(0))
        
        cycles[[idx]] <- list(df = ode_df, meas = meas_df,
                              r_vec = attr(ode_df, "r_vec"),
                              K = attr(ode_df, "K"),
                              drug  = attr(ode_df, "drug"),
                              phi_params_list = attr(ode_df, "phi_params_list"), 
                              C_func = attr(ode_df, "C_func"),
                              k_multiplier_base = k_mult_base,
                              df_fitted = NULL,
                              k_multiplier_fitted = NULL) 
        all_cycles(cycles)
        
        last_B(tail(ode_df$B_total, 1))
    })
    
    observeEvent(input$addMeas, {
        meas <- parse_measurements(input$measurements)
        validate(need(nrow(meas) > 0, "Please enter valid measurement pairs (%burden time)."))
        
        cycles <- all_cycles()
        validate(need(length(cycles) > 0, "Run a cycle before adding measurements."))
        
        cur_idx <- length(cycles)
        cur <- cycles[[cur_idx]]
        
        B_start <- cur$df$B_total[1]
        meas$B <- (meas$percent / 100) * B_start
        
        base_cum_t0 <- min(cur$df$cum_time)
        meas$cum_time <- meas$time + base_cum_t0
        meas <- subset(meas, meas$time >= 0 & meas$time <= max(cur$df$time))
        
        cur$meas <- meas 
        cur$df_fitted <- NULL
        cur$k_multiplier_fitted <- NULL
        
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
        
        par(mar = c(4, 4, 3, 1))
        plot(df$time, df$B_total, type = "l", lwd = 2, col = "firebrick",
             xlab = "Time (days, current cycle)", ylab = "Tumor burden (B)",
             main = paste("Total Burden (Drug:", cur$drug, ")"),
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
            "Final Comp:   ", frac_final_str
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
