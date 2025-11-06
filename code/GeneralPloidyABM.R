#---------------------------
# Parameters
#---------------------------
Dc <- 0          # baseline death probability
Coxy <- 1        # oxygen level (0-1)
P <- c(2,3,4)    # initial ploidy distribution
DS <- 0.6        # division probability weight
DS1 <- 0.4       # WGD probability weight
MSR <- 0.05      # base missegregation probability
SP <- 0.5        # WGD protection factor
SD <- 0.2        # missegregation death factor
time_steps <- 5 # total simulation time

# Function: rate of entering cell cycle
G <- function(ploidy, oxygen) {
  return(0.5 * oxygen / ploidy) # example scaling
}

#---------------------------
# Initialize Cells
#---------------------------
N <- 100
cells <- data.frame(
  ploidy = sample(P, N, replace=TRUE),
  WGD = FALSE,
  t_wait = 0,
  alive = TRUE
)

#---------------------------
# Simulation Loop
#---------------------------
for (t in 1:time_steps) {
  for (i in 1:nrow(cells)) {
    if (!cells$alive[i]) next
    
    #-----------------------
    # Cell Death
    #-----------------------
    misseg <- rbinom(1, 1, MSR * (1 - Coxy)) # higher missegregation in low O2
    death_prob <- Dc
    if (misseg == 1) {
      if (cells$WGD[i]) {
        death_prob <- death_prob + SD * (1 - SP)
      } else {
        death_prob <- death_prob + SD
      }
    }
    if (runif(1) < death_prob) {
      cells$alive[i] <- FALSE
      next
    }
    
    #-----------------------
    # Enter cell cycle
    #-----------------------
    if (runif(1) < G(cells$ploidy[i], Coxy)) {
      # Set waiting time if not already set
      if (cells$t_wait[i] == 0) {
        lambda <- 1 # rate of division/WGD
        cells$t_wait[i] <- rpois(1, lambda)
      }
    }
    
    #-----------------------
    # Event occurs if t_wait = 0
    #-----------------------
    if (cells$t_wait[i] <= 0) {
      # Decide Division vs WGD
      if (runif(1) < DS / (DS + DS1)) {
        # Division
        cells$ploidy[i] <- cells$ploidy[i] # parent keeps ploidy
        # Optionally add daughter cell
        cells <- rbind(cells, data.frame(ploidy=cells$ploidy[i], WGD=cells$WGD[i], t_wait=0, alive=TRUE))
      } else {
        # WGD
        cells$ploidy[i] <- 2 * cells$ploidy[i]
        cells$WGD[i] <- TRUE
      }
      # Reset waiting time
      cells$t_wait[i] <- 0
    } else {
      cells$t_wait[i] <- cells$t_wait[i] - 1
    }
  }
  
  cat("Time step:", t, "Alive cells:", sum(cells$alive), "\n")
}

#---------------------------
# Quick Visualization
#---------------------------
library(ggplot2)
ggplot(cells, aes(x=ploidy, fill=WGD)) +
  geom_histogram(binwidth=1, position="dodge") +
  theme_minimal() +
  labs(title="Ploidy Distribution After Simulation", x="Ploidy", y="Cell Count")


ggplot(cells, aes(x=X, y=Y, color=P, shape=WGD)) +
  geom_point(alpha=0.7) +
  scale_color_viridis_c(option = "plasma") +
  scale_size(range=c(2,6)) +
  scale_shape_manual(values=c(16,17)) +
  labs(title="ABM: Oxygen-dependent Division, WGD, Missegregation",
       x="X", y="Y", color="Ploidy", shape="WGD") +
  theme_minimal() +
  coord_fixed()

library(ggplot2)
library(viridis)

#---------------------------
# Parameters
#---------------------------
Dc <- 0
Coxy <- 1
P <- c(2,3,4)
DS <- 0.6
DS1 <- 0.4
MSR <- 0.05
SP <- 0.5
SD <- 0.2
time_steps <- 50
space_size <- 50 # for X,Y coordinates

G <- function(ploidy, oxygen) {
  return(0.5 * oxygen / ploidy)
}

#---------------------------
# Initialize Cells with Spatial Coordinates
#---------------------------
N <- 100
cells <- data.frame(
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
for (t in 1:time_steps) {
  for (i in 1:nrow(cells)) {
    if (!cells$alive[i]) next
    
    # Cell death
    misseg <- rbinom(1, 1, MSR * (1 - Coxy))
    death_prob <- Dc
    if (misseg == 1) {
      if (cells$WGD[i]) death_prob <- death_prob + SD * (1 - SP) else death_prob <- death_prob + SD
    }
    if (runif(1) < death_prob) {
      cells$alive[i] <- FALSE
      next
    }
    
    # Enter cell cycle
    if (runif(1) < G(cells$ploidy[i], Coxy)) {
      if (cells$t_wait[i] == 0) cells$t_wait[i] <- rpois(1, 1)
    }
    
    # Event occurs if t_wait = 0
    if (cells$t_wait[i] <= 0) {
      if (runif(1) < DS / (DS + DS1)) {
        # Division: create daughter near parent
        new_cell <- data.frame(
          ploidy = cells$ploidy[i],
          WGD = cells$WGD[i],
          t_wait = 0,
          alive = TRUE,
          X = cells$X[i] + rnorm(1, 0, 1),
          Y = cells$Y[i] + rnorm(1, 0, 1)
        )
        # Keep within bounds
        new_cell$X <- pmin(pmax(new_cell$X, 0), space_size)
        new_cell$Y <- pmin(pmax(new_cell$Y, 0), space_size)
        cells <- rbind(cells, new_cell)
      } else {
        # WGD
        cells$ploidy[i] <- 2 * cells$ploidy[i]
        cells$WGD[i] <- TRUE
      }
      cells$t_wait[i] <- 0
    } else {
      cells$t_wait[i] <- cells$t_wait[i] - 1
    }
  }
}

#---------------------------
# Visualization
#---------------------------
ggplot(cells, aes(x=X, y=Y, color=ploidy, shape=WGD)) +
  geom_point(alpha=0.7, size=3) +
  scale_color_viridis_c(option = "plasma") +
  scale_shape_manual(values=c(16,17)) +
  labs(title="ABM: Oxygen-dependent Division, WGD, Missegregation",
       x="X", y="Y", color="Ploidy", shape="WGD") +
  theme_minimal() +
  coord_fixed()


library(ggplot2)
library(viridis)
library(gganimate)

#---------------------------
# Parameters
#---------------------------
Dc <- 0
Coxy <- 1
P <- c(2,3,4)
DS <- 0.6
DS1 <- 0.4
MSR <- 0.05
SP <- 0.5
SD <- 0.2
time_steps <- 5
space_size <- 50

G <- function(ploidy, oxygen) {
  return(0.5 * oxygen / ploidy)
}

#---------------------------
# Initialize Cells with Spatial Coordinates
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
    
    # Cell death
    misseg <- rbinom(1, 1, MSR * (1 - Coxy))
    death_prob <- Dc
    if (misseg == 1) {
      if (cells$WGD[i]) death_prob <- death_prob + SD * (1 - SP) else death_prob <- death_prob + SD
    }
    if (runif(1) < death_prob) {
      cells$alive[i] <- FALSE
      next
    }
    
    # Enter cell cycle
    if (runif(1) < G(cells$ploidy[i], Coxy)) {
      if (cells$t_wait[i] == 0) cells$t_wait[i] <- rpois(1, 1)
    }
    
    # Event occurs if t_wait = 0
    if (cells$t_wait[i] <= 0) {
      if (runif(1) < DS / (DS + DS1)) {
        # Division: create daughter near parent
        new_cell <- data.frame(
          id = max(cells$id) + nrow(new_cells) + 1,
          ploidy = cells$ploidy[i],
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
        # WGD
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

#---------------------------
# Animated Visualization
#---------------------------
library(ggplot2)

# Convert WGD to a factor for coloring
history_df$WGD <- factor(history_df$WGD, levels=c(FALSE, TRUE), labels=c("No WGD", "WGD"))

# Plot
ggplot(history_df, aes(x = X, y = Y, color = WGD, size = ploidy)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("blue", "red")) +
  scale_size_continuous(range = c(2, 6)) +
  facet_wrap(~time, ncol = 5) +
  theme_minimal() +
  labs(title = "Cell Population Over Time",
       subtitle = "Color: WGD status, Size: Ploidy",
       x = "X position",
       y = "Y position") +
  theme(legend.position = "bottom")
