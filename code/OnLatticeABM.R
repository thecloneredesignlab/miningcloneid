library(ggplot2)
library(dplyr)
##Change function so that the Daughter cell dies with a random distribution
##Less likely with whole genome doubling
# Parameters
Dc <- .01
Coxy <- 0
P <- c(2,3,4)
DS <- 0.6
DS1 <- 0.4
MSR <- 0.05
SP <- 0.5
SD <- 0.2
time_steps <- 10
space_size <- 50  # lattice is 0..50 in both X and Y

G <- function(ploidy, oxygen) {
  0.5 * oxygen / ploidy
}

# Initialize cells on lattice
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

history <- list()

for (t in 1:time_steps) {
  new_cells <- data.frame()
  
  for (i in 1:nrow(cells)) {
    if (!cells$alive[i]) next
    
    # Cell death
    misseg <- rbinom(1, 1, MSR * (1 - Coxy))
    death_prob <- Dc
    if (misseg == 1) {
      death_prob <- death_prob + ifelse(cells$WGD[i], SD*(1-SP), SD)
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
        # Division with lattice placement
        daughter_ploidy <- cells$ploidy[i]
        if (rbinom(1,1,MSR)==1) daughter_ploidy <- round(daughter_ploidy * runif(1,0.5,1))
        
        threshold <- ifelse(cells$WGD[i], 0.15*cells$ploidy[i], 0.3*cells$ploidy[i])
        if (daughter_ploidy >= threshold) {
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
        # WGD
        cells$ploidy[i] <- 2*cells$ploidy[i]
        cells$WGD[i] <- TRUE
      }
      cells$t_wait[i] <- 0
    } else {
      cells$t_wait[i] <- cells$t_wait[i] - 1
    }
  }
  
  if (nrow(new_cells)>0) cells <- rbind(cells, new_cells)
  history[[t]] <- data.frame(cells, time=t)
}

history_df <- do.call(rbind, history)
history_df$ploidy <- pmax(history_df$ploidy, 2)
history_df$WGD <- factor(history_df$WGD, levels=c(FALSE,TRUE), labels=c("No WGD","WGD"))

# Plot on lattice
ggplot(history_df, aes(x=X, y=Y, color=WGD, size=ploidy)) +
  geom_jitter(width=0.3, height=0.3, alpha=0.7) +  # jitter to show overlapping lattice cells
  scale_color_manual(values=c("blue","red")) +
  scale_size_continuous(range=c(2,6)) +
  facet_wrap(~time, ncol=5) +
  theme_minimal() +
  labs(title="Cell Population on Lattice Over Time",
       subtitle="Color: WGD status, Size: Ploidy",
       x="X lattice coordinate",
       y="Y lattice coordinate") +
  theme(legend.position="bottom")
