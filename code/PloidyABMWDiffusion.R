library(ggplot2)
library(dplyr)

#---------------------------
# Parameters
#---------------------------
Dc <- 0.01           # baseline death probability
Coxy <- 0            # oxygen level (1 = high, 0 = low)
P <- c(2,3,4)        # starting ploidy
DS <- 0.06
DS1 <- 0.004
MSR <- 0.05          # missegregation probability
SP <- 0.5            # survival probability adjustment for WGD
SD <- 0.2
time_steps <- 5
space_size <- 50
max_ploidy <- 10

# Function: probability to enter cell cycle depends on oxygen and ploidy
G <- function(ploidy, oxygen) {
  0.5 * oxygen / ploidy
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
  X = sample(0:space_size, N, replace=TRUE),
  Y = sample(0:space_size, N, replace=TRUE)
)

# Initialize oxygen field (constant for simplicity)
oxygen <- matrix(Coxy, nrow = space_size+1, ncol = space_size+1)

history <- list()

#---------------------------
# Simulation Loop
#---------------------------
for (t in 1:time_steps) {
  new_cells <- data.frame()
  
  for (i in 1:nrow(cells)) {
    if (!cells$alive[i]) next
    
    # Get local oxygen
    ox <- oxygen[cells$X[i]+1, cells$Y[i]+1]
    
    # Cell death due to missegregation
    misseg <- rbinom(1, 1, MSR * (1 - ox))
    death_prob <- Dc
    if (misseg == 1) {
      death_prob <- death_prob + ifelse(cells$WGD[i], SD*(1-SP), SD)
    }
    if (runif(1) < death_prob) {
      cells$alive[i] <- FALSE
      next
    }
    
    # Enter cell cycle
    if (runif(1) < G(cells$ploidy[i], ox)) {
      if (cells$t_wait[i] == 0) cells$t_wait[i] <- rpois(1, 1)
    }
    
    # Event occurs if t_wait = 0
    if (cells$t_wait[i] <= 0) {
      if (runif(1) < DS / (DS + DS1)) {
        # Division
        daughter_ploidy <- cells$ploidy[i]
        if (rbinom(1,1,MSR)==1) {
          daughter_ploidy <- round(daughter_ploidy * runif(1,0.5,1))
        }
        threshold <- ifelse(cells$WGD[i], 0.15*cells$ploidy[i], 0.3*cells$ploidy[i])
        
        # Check viability and max ploidy
        if (!is.na(daughter_ploidy) && daughter_ploidy >= threshold && daughter_ploidy <= max_ploidy) {
          # Snap daughter to nearby lattice site
          dx <- sample(c(-1,0,1),1)
          dy <- sample(c(-1,0,1),1)
          new_cell <- data.frame(
            id = max(cells$id) + nrow(new_cells) + 1,
            ploidy = daughter_ploidy,
            WGD = cells$WGD[i],
            t_wait = 0,
            alive = TRUE,
            X = min(max(cells$X[i] + dx, 0), space_size),
            Y = min(max(cells$Y[i] + dy, 0), space_size)
          )
          new_cells <- rbind(new_cells, new_cell)
        }
      } else {
        # Whole Genome Doubling
        cells$ploidy[i] <- min(cells$ploidy[i]*2, max_ploidy)
        cells$WGD[i] <- TRUE
      }
      cells$t_wait[i] <- 0
    } else {
      cells$t_wait[i] <- cells$t_wait[i] - 1
    }
  }
  
  if (nrow(new_cells) > 0) cells <- rbind(cells, new_cells)
  history[[t]] <- data.frame(cells, time=t)
}

# Combine history
history_df <- do.call(rbind, history)
history_df$ploidy <- pmin(history_df$ploidy, max_ploidy)
history_df$WGD <- factor(history_df$WGD, levels=c(FALSE,TRUE), labels=c("No WGD","WGD"))

#---------------------------
# Plot: Cells on lattice with oxygen as background
#---------------------------
# Create oxygen data frame for plotting
oxygen_df <- expand.grid(X=0:space_size, Y=0:space_size)
oxygen_df$oxygen <- as.vector(oxygen)

ggplot() +
  geom_tile(data=oxygen_df, aes(x=X, y=Y, fill=oxygen), alpha=0.3) +
  scale_fill_gradient(low="lightyellow", high="darkgreen", name="Oxygen") +
  geom_jitter(data=history_df, aes(x=X, y=Y, color=WGD, size=ploidy), width=0.3, height=0.3, alpha=0.7) +
  scale_color_manual(values=c("blue","red")) +
  scale_size_continuous(range=c(2,6)) +
  facet_wrap(~time, ncol=5) +
  theme_minimal() +
  labs(title="Cell Population on Lattice with Oxygen",
       subtitle="Color: WGD status, Size: Ploidy",
       x="X lattice coordinate",
       y="Y lattice coordinate",
       fill="Oxygen") +
  theme(legend.position="bottom")
library(dplyr)
library(ggplot2)

# Summarize over time
summary_df <- history_df %>%
  group_by(time) %>%
  summarise(
    mean_ploidy = mean(ploidy, na.rm = TRUE),
    median_ploidy = median(ploidy, na.rm = TRUE),
    WGD_fraction = mean(WGD == "WGD")
  )

# Plot mean and median ploidy
p1 <- ggplot(summary_df, aes(x = time)) +
  geom_line(aes(y = mean_ploidy, color = "Mean ploidy"), size=1) +
  geom_line(aes(y = median_ploidy, color = "Median ploidy"), size=1, linetype="dashed") +
  scale_color_manual(values=c("Mean ploidy"="blue","Median ploidy"="red")) +
  labs(title="Ploidy over Time", y="Ploidy", color="") +
  theme_minimal()

# Plot WGD fraction
p2 <- ggplot(summary_df, aes(x = time, y = WGD_fraction)) +
  geom_line(color="darkgreen", size=1) +
  labs(title="Fraction of Cells with WGD over Time",
       y="Fraction WGD",
       x="Time") +
  theme_minimal()

# Display plots
library(gridExtra)
grid.arrange(p1, p2, ncol=1)