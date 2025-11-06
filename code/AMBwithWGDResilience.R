Dc <- 0
Coxy <- 1
P <- c(2,3,4)
DS <- 0.6
DS1 <- 0.4
MSR <- 0.05
SP <- 0.5
SD <- 0.2
time_steps <- 20
space_size <- 50

G <- function(ploidy, oxygen) {
  return(0.5 * oxygen / ploidy)
}

#---------------------------
# Initialize Cells
#---------------------------
N <- 100
cells <- data.frame(
  id = 1:N,
  ploidy = sample(P, N, replace=TRUE),
  WGD = FALSE,
  t_wait = 0,
  alive = TRUE,
  X = runif(N, 0, space_size),
  Y = runif(N, 0, space_size)
)

#---------------------------
# Simulation Loop
#---------------------------
history <- list()

for (t in 1:time_steps) {
  new_cells <- data.frame()
  
  for (i in 1:nrow(cells)) {
    if (!cells$alive[i]) next
    
    # ---------------------------
    # Cell death
    # ---------------------------
    misseg <- rbinom(1, 1, MSR * (1 - Coxy))
    death_prob <- Dc
    if (misseg == 1) {
      if (cells$WGD[i]) death_prob <- death_prob + SD * (1 - SP) else death_prob <- death_prob + SD
    }
    if (runif(1) < death_prob) {
      cells$alive[i] <- FALSE
      next
    }
    
    # ---------------------------
    # Enter cell cycle
    # ---------------------------
    if (runif(1) < G(cells$ploidy[i], Coxy)) {
      if (cells$t_wait[i] == 0) cells$t_wait[i] <- rpois(1, 1)
    }
    
    # ---------------------------
    # Event occurs if t_wait = 0
    # ---------------------------
    if (cells$t_wait[i] <= 0) {
      if (runif(1) < DS / (DS + DS1)) {
        # ---------- Division with possible missegregation ----------
        daughter_ploidy <- cells$ploidy[i]
        
        # Missegregation occurs with probability MSR
        if (rbinom(1, 1, MSR) == 1) {
          # Reduce daughter ploidy randomly between 50% and 100% of parent
          daughter_ploidy <- round(daughter_ploidy * runif(1, 0.5, 1))
        }
        
        # Check viability thresholds
        threshold <- ifelse(cells$WGD[i], 0.15 * cells$ploidy[i], 0.3 * cells$ploidy[i])
        if (daughter_ploidy < threshold) {
          # Daughter dies immediately, skip adding
          cells$t_wait[i] <- 0
          next
        }
        
        # Create daughter cell
        new_cell <- data.frame(
          id = max(cells$id) + nrow(new_cells) + 1,
          ploidy = daughter_ploidy,
          WGD = cells$WGD[i],
          t_wait = 0,
          alive = TRUE,
          X = cells$X[i] + rnorm(1, 0, 1),
          Y = cells$Y[i] + rnorm(1, 0, 1)
        )
        new_cell$X <- pmin(pmax(new_cell$X, 0), space_size)
        new_cell$Y <- pmin(pmax(new_cell$Y, 0), space_size)
        new_cells <- rbind(new_cells, new_cell)
        
      } else {
        # ---------- WGD ----------
        cells$ploidy[i] <- 2 * cells$ploidy[i]
        cells$WGD[i] <- TRUE
      }
      cells$t_wait[i] <- 0
    } else {
      cells$t_wait[i] <- cells$t_wait[i] - 1
    }
  }
  
  # Add new cells
  if (nrow(new_cells) > 0) cells <- rbind(cells, new_cells)
  
  # Save snapshot for animation
  history[[t]] <- data.frame(cells, time = t)
}

# Combine history into one data frame
history_df <- do.call(rbind, history)

history_df$ploidy <- pmax(history_df$ploidy, 2)

# Convert WGD to a factor for coloring
history_df$WGD <- factor(history_df$WGD, levels = c(FALSE, TRUE), labels = c("No WGD", "WGD"))

# Plot: spatial positions of cells over time
ggplot(history_df, aes(x = X, y = Y, color = WGD, size = ploidy)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("blue", "red")) +
  scale_size_continuous(range = c(2, 6)) +
  facet_wrap(~time, ncol = time_steps) +
  theme_minimal() +
  labs(
    title = "Cell Population Over Time",
    subtitle = "Color: WGD status, Size: Ploidy",
    x = "X position",
    y = "Y position"
  ) +
  theme(legend.position = "bottom")