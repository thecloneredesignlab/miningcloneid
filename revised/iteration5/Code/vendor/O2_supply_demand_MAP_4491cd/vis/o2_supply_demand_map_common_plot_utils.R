#!/usr/bin/env Rscript

# Visualization-only scales shared by in-vivo and in-vitro plot constructors.
# No fitted inputs, simulations, analyses, or reports are read or produced here.

o2sd_ploidy_fraction_fill_scale <- function(fill_max, name = "Fraction") {
  fill_max <- suppressWarnings(as.numeric(fill_max))
  if (!is.finite(fill_max) || fill_max <= 0) fill_max <- 1
  ggplot2::scale_fill_gradientn(
    colors = c("#f7f7f7", "#2c7fb8", "#ffff33"),
    values = scales::rescale(
      c(0, 0.05 * fill_max, fill_max),
      from = c(0, fill_max)
    ),
    limits = c(0, fill_max),
    oob = scales::squish,
    na.value = "white",
    name = name
  )
}
