library(shiny)
library(deSolve)

drug_dosing_schedules <- list(
  "alisertib" = "IV", "gemcitabine" = "IV", "bay1895344" = "IV",
  "ispinesib" = "IV", "none" = "IV"
)

IV_DEFAULTS <- list(C_peak=1.0, half_life=0.5, period=7.0)
PER_DRUG <- list(
  "alisertib"    = list(C_peak = 42.5, half_life = 0.875, period = 0.5),
  "gemcitabine"  = list(C_peak = 239, half_life = 0.05, period = 7.0),
  "ispinesib"    = list(C_peak = 2.1, half_life = 1.04, period = 7),
  "bay1895344"   = list(C_peak = 6.2, half_life = 0.50, period = 0.5),
  "none"         = list(C_peak = 0, half_life = 1, period = 7)
)

phi_Hill <- function(C, EC50, n, Emax=1.0) Emax * (C^n) / (EC50^n + C^n)

f_pd_params <- function(ploidy, drug) {
  drug <- tolower(drug)
  clamp_ec50 <- function(x) max(x, 1e-12)
  clamp_n <- function(x) max(x, 0.1)

  if (drug == "bay1895344") {
    list(EC50=clamp_ec50(1.04 * exp(0.35 * ploidy) - 2.05), n=clamp_n(3.85 * exp(-0.861 * ploidy) + 0.81), Emax=1.0)
  } else if (drug == "alisertib") {
    list(EC50=clamp_ec50(51.02 * exp(-0.62 * ploidy) - 4.78), n=clamp_n(1.0), Emax=1.0)
  } else if (drug == "ispinesib") {
    list(EC50=clamp_ec50(1.185 * exp(-0.21 * ploidy) - 0.56), n=clamp_n(0.94 * exp(-0.303 * ploidy) - 0.73), Emax=1.0)
  } else if (drug == "gemcitabine") {
    list(EC50=clamp_ec50(0.004 * exp(0.78 * ploidy) - 0.01), n=clamp_n(28.92 * exp(-0.94 * ploidy) + 0.92), Emax=1.0)
  } else {
    list(EC50=1e10, n=1.0, Emax=0.0) # No drug effect
  }
}

pulsed_dose <- function(C_peak=5.0, half_life=2.0, period=7.0) {
  lam <- log(2) / max(half_life, 1e-12)
  function(t) { modt <- as.numeric(t) %% period; C_peak * exp(-lam * modt) }
}

get_concentration_curve <- function(drug_name) {
  p <- PER_DRUG[[tolower(drug_name)]]
  if (is.null(p)) p <- list(C_peak=0, half_life=1, period=7)
  do.call(pulsed_dose, p)
}

model_ode_fn <- function(t, state, parms) {
  B_total <- sum(state)
  phi_vals <- sapply(parms$phi_params_list, function(p) phi_Hill(parms$C_func(t), p$EC50, p$n, p$Emax))
  dB_vec <- ((parms$r_vec * max(0, 1 - B_total / parms$K)) - (phi_vals * parms$k_multiplier)) * state
  list(dB_vec)
}


# Helper to run a fast simulation for MCTS steps
simulate_step <- function(current_state, drug, days=28, k_mult=1.0) {
  n <- length(current_state)
  r_vec <- rep_len(c(0.020, 0.015, 0.010), n)
  phi_params <- lapply(seq(2, length.out=n), f_pd_params, drug=drug)

  out <- ode(y=current_state, times=c(0, days), func=model_ode_fn,
             parms=list(r_vec=r_vec, K=2e10, phi_params_list=phi_params,
                        C_func=get_concentration_curve(drug), k_multiplier=k_mult))
  return(as.numeric(out[2, 2:(n+1)]))
}

# MCTS Node as an Environment (for mutable state in R)
create_node <- function(state, cycle, parent=NULL) {
  node <- new.env()
  node$state <- state
  node$cycle <- cycle
  node$parent <- parent
  node$children <- list()
  node$N <- 0
  node$W <- 0.0
  node
}

run_mcts <- function(root_state, drugs, num_rollouts=20, depth=3) {
  root <- create_node(root_state, 0)
  min_size <- 1e5
  max_size <- 2e10

  for (i in 1:num_rollouts) {
    node <- root

    # 1. Selection
    while (length(node$children) == length(drugs) && node$cycle < depth) {
      scores <- sapply(node$children, function(child) {
        (child$W / (child$N + 1e-6)) + 1.4 * sqrt(log(node$N + 1) / (child$N + 1e-6))
      })
      node <- node$children[[which.max(scores)]]
    }

    # 2. Expansion
    if (node$cycle < depth) {
      untried <- setdiff(drugs, names(node$children))
      drug <- sample(untried, 1)
      new_state <- simulate_step(node$state, drug)
      child <- create_node(new_state, node$cycle + 1, parent=node)
      node$children[[drug]] <- child
      node <- child
    }

    # 3. Rollout (Simplified)
    curr_s <- node$state
    for (d in 1:(depth - node$cycle)) {
      if (sum(curr_s) < min_size || sum(curr_s) > max_size) break
      curr_s <- simulate_step(curr_s, sample(drugs, 1))
    }

    # 4. Backprop
    reward <- -sum(curr_s)
    if (sum(curr_s) < min_size) reward <- 1e12 # Extinction bonus

    curr_back <- node
    while (!is.null(curr_back)) {
      curr_back$N <- curr_back$N + 1
      curr_back$W <- curr_back$W + reward
      curr_back <- curr_back$parent
    }
  }

  # Return drug with highest average reward
  avg_rewards <- sapply(root$children, function(nc) nc$W / nc$N)
  return(names(avg_rewards)[which.max(avg_rewards)])
}


ui <- fluidPage(
  titlePanel("Tumor Model with MCTS Drug Optimization"),
  sidebarLayout(
    sidebarPanel(
      textInput("fractions", "Ploidy fractions", value="0.6,0.3,0.1"),
      numericInput("B0", "Initial burden", value=1e9),
      selectInput("drug", "Manual Drug Selection", choices=names(PER_DRUG)),
      actionButton("runCycle", "Run Selected Drug", class="btn-primary"),
      hr(),
      helpText("Reccomended Drug (MCTS):"),
      actionButton("optimizeDrug", "Find & Run Best Drug", class="btn-success"),
      hr(),
      actionButton("reset", "Reset Simulation")
    ),
    mainPanel(
      plotOutput("trajPlot"),
      verbatimTextOutput("log")
    )
  )
)

server <- function(input, output, session) {
  history <- reactiveVal(list())
  current_state <- reactiveVal(NULL)

  observeEvent(input$runCycle, { run_logic(input$drug) })

  observeEvent(input$optimizeDrug, {
    showNotification("MCTS engine calculating best drug...")
    state <- if(is.null(current_state())) parse_fractions(input$fractions)*input$B0 else current_state()
    best <- run_mcts(state, names(PER_DRUG), num_rollouts=50)
    updateSelectInput(session, "drug", selected=best)
    run_logic(best)
  })

  run_logic <- function(drug_name) {
    fr <- if(is.null(current_state())) parse_fractions(input$fractions) else current_state()/sum(current_state())
    B_start <- if(is.null(current_state())) input$B0 else sum(current_state())

    n <- length(fr)
    phi_list <- lapply(seq(2, length.out=n), f_pd_params, drug=drug_name)

    parms <- list(r_vec=rep_len(c(0.02, 0.015, 0.01), n), K=2e10,
                  phi_params_list=phi_list, C_func=get_concentration_curve(drug_name), k_multiplier=1.0)

    out <- as.data.frame(ode(y=fr*B_start, times=seq(0, 28, by=0.5), func=model_ode_fn, parms=parms))

    # Update global state
    h <- history()
    out$cycle <- length(h) + 1
    out$drug <- drug_name
    h[[length(h) + 1]] <- out
    history(h)
    current_state(as.numeric(tail(out, 1)[, 2:(n+1)]))
  }

  parse_fractions <- function(txt) {
    v <- as.numeric(unlist(strsplit(txt, ","))); v/sum(v)
  }

  observeEvent(input$reset, { history(list()); current_state(NULL) })

  output$trajPlot <- renderPlot({
    h <- history()
    if(length(h) == 0) return(NULL)
    all_data <- do.call(rbind, h)
    all_data$total <- rowSums(all_data[, 2:(ncol(all_data)-2)])
    all_data$global_time <- 1:nrow(all_data) # Proxy for simplicity

    plot(all_data$total, type="l", log="y", lwd=2, col="blue",
         xlab="Time Steps", ylab="Total Cells", main="Treatment Trajectory")
    abline(h=1e5, col="red", lty=2) # Extinction line
  })

  output$log <- renderText({
    h <- history()
    if(length(h) == 0) return("No cycles run.")
    paste("Cycle History:", paste(sapply(h, function(x) unique(x$drug)), collapse=" -> "))
  })
}

shinyApp(ui, server)