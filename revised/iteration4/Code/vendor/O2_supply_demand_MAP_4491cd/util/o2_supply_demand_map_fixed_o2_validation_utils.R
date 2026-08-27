#!/usr/bin/env Rscript

# Canonical fixed-O2 CLI validation shared by simulation and legacy runners.

fixo2_validate_mode_reference_o2 <- function(mode_reference_o2, attractor_o2_grid) {
  mode_reference_o2 <- suppressWarnings(as.numeric(mode_reference_o2))
  attractor_o2_grid <- sort(unique(suppressWarnings(as.numeric(attractor_o2_grid))))
  attractor_o2_grid <- attractor_o2_grid[is.finite(attractor_o2_grid)]
  if (!is.finite(mode_reference_o2)) {
    stop("--mode_reference_o2 must be a finite numeric O2 value.", call. = FALSE)
  }
  if (!length(attractor_o2_grid)) {
    stop("--attractor_o2_grid must contain at least one finite numeric O2 value.", call. = FALSE)
  }
  if (!any(abs(attractor_o2_grid - mode_reference_o2) < 1e-9)) {
    stop(
      "--mode_reference_o2=",
      format(mode_reference_o2, scientific = FALSE, trim = TRUE),
      " is invalid. It must exactly match one value in --attractor_o2_grid. Available attractor O2 values: ",
      fixo2_format_o2_list(attractor_o2_grid),
      call. = FALSE
    )
  }
  invisible(mode_reference_o2)
}
