# In-vivo o2 simulation products. This module is data-only.

make_o2_crit_adaptive_grid <-
function(o2_crit, grid_min = 0, grid_max = 5, base_step = 0.02, dense_n = 220L) {
    grid_min <- as.numeric(grid_min[[1L]])
    grid_max <- as.numeric(grid_max[[1L]])
    base_step <- as.numeric(base_step[[1L]])
    dense_n <- as.integer(dense_n[[1L]])
    if (!is.finite(grid_min) || !is.finite(grid_max))
        stop("O2 adaptive grid bounds must be finite.")
    if (grid_max < grid_min) {
        tmp <- grid_min
        grid_min <- grid_max
        grid_max <- tmp
    }
    if (!is.finite(base_step) || base_step <= 0)
        base_step <- max((grid_max - grid_min)/250, 1e-06)
    if (!is.finite(dense_n) || dense_n < 2L)
        dense_n <- 220L
    base_grid <- seq(grid_min, grid_max, by = base_step)
    base_grid <- c(base_grid, grid_min, grid_max)
    crit <- suppressWarnings(as.numeric(o2_crit[[1L]]))
    if (!is.finite(crit) || crit <= 0) {
        return(sort(unique(signif(base_grid[is.finite(base_grid)], 14))))
    }
    near_upper <- min(grid_max, max(base_step, 25 * crit))
    near_zero_grid <- if (near_upper > grid_min) {
        seq(grid_min, near_upper, length.out = dense_n)
    }
    else {
        numeric(0)
    }
    transition_lower <- max(grid_min, crit * 0.05)
    transition_upper <- min(grid_max, crit * 20)
    transition_grid <- if (transition_upper > transition_lower) {
        seq(transition_lower, transition_upper, length.out = dense_n)
    }
    else {
        numeric(0)
    }
    log_grid <- crit * 10^seq(-4, 2, length.out = dense_n)
    multiplier_grid <- crit * c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 5,
        8, 10, 15, 20, 30, 50, 100)
    grid <- c(base_grid, near_zero_grid, transition_grid, log_grid, multiplier_grid)
    grid <- grid[is.finite(grid) & grid >= grid_min & grid <= grid_max]
    sort(unique(signif(grid, 14)))
}

make_o2_crit_reference_levels <-
function(o2_crit, grid_min = 0, grid_max = 5, coarse_step = 0.5) {
    grid_min <- as.numeric(grid_min[[1L]])
    grid_max <- as.numeric(grid_max[[1L]])
    coarse_step <- as.numeric(coarse_step[[1L]])
    if (!is.finite(grid_min) || !is.finite(grid_max))
        stop("O2 reference level bounds must be finite.")
    if (grid_max < grid_min) {
        tmp <- grid_min
        grid_min <- grid_max
        grid_max <- tmp
    }
    if (!is.finite(coarse_step) || coarse_step <= 0)
        coarse_step <- max((grid_max - grid_min)/10, 1e-06)
    coarse_levels <- seq(grid_min, grid_max, by = coarse_step)
    crit <- suppressWarnings(as.numeric(o2_crit[[1L]]))
    crit_levels <- if (is.finite(crit) && crit > 0) {
        crit * c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 25)
    }
    else {
        numeric(0)
    }
    levels <- c(grid_min, grid_max, coarse_levels, crit_levels)
    levels <- levels[is.finite(levels) & levels >= grid_min & levels <= grid_max]
    sort(unique(signif(levels, 12)))
}
