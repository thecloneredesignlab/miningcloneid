# Required libraries
library(ggplot2)
library(gganimate)
library(viridis)

# Example ABM cells dataframe
# Replace this with your actual simulation output
set.seed(123)
cells <- data.frame(
  X = runif(200, 0, 100),
  Y = runif(200, 0, 100),
  P = sample(2:8, 200, replace = TRUE),
  WGD = sample(c(0,1), 200, replace = TRUE),
  time = rep(1:10, each = 20)  # 10 time steps
)

# Base ggplot
p <- ggplot(cells, aes(x = X, y = Y, color = P, shape = factor(WGD))) +
  geom_point(alpha = 0.7, size = 4) +
  scale_color_viridis_c(option = "plasma") +
  scale_shape_manual(values = c(16, 17), labels = c("No WGD", "WGD")) +
  labs(
    title = "ABM: Oxygen-dependent Division, WGD, Missegregation",
    x = "X",
    y = "Y",
    color = "Ploidy",
    shape = "WGD"
  ) +
  theme_minimal() +
  coord_fixed() +
  theme(plot.title = element_text(hjust = 0.5))

# Add animation
p_anim <- p +
  transition_states(time, state_length = 1) +
  ease_aes('linear') +
  labs(subtitle = 'Time: {closest_state}')

# Render animation
animate(p_anim, nframes = 100, fps = 10)
