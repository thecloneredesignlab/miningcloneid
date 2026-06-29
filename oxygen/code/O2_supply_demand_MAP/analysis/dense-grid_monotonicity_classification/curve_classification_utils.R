curve_classification_rule_version <- function() {
  "shape_first_v3"
}

curve_step_epsilon <- function(ploidy_range,
                               step_epsilon_abs = 1e-6,
                               step_epsilon_fraction = 1e-4) {
  scale <- if (is.finite(ploidy_range) && ploidy_range > 0) ploidy_range else 1
  max(step_epsilon_abs, step_epsilon_fraction * scale)
}

collapse_signs <- function(signs) {
  s <- signs[is.finite(signs) & !is.na(signs)]
  nz <- s[s != 0L]
  if (!length(nz)) return(integer(0))
  rle(nz)$values
}

finite_diff_curve <- function(curve,
                              step_epsilon = NULL,
                              value_col = "dominant_mean_ploidy",
                              x_col = "O2_pct",
                              id_col = "seed_id",
                              step_epsilon_abs = 1e-6,
                              step_epsilon_fraction = 1e-4) {
  if (!value_col %in% names(curve)) stop("Missing value column: ", value_col, call. = FALSE)
  if (!x_col %in% names(curve)) stop("Missing x column: ", x_col, call. = FALSE)
  curve <- curve[order(curve[[x_col]]), , drop = FALSE]
  y <- suppressWarnings(as.numeric(curve[[value_col]]))
  pr <- suppressWarnings(max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
  if (!is.finite(pr)) pr <- NA_real_
  if (is.null(step_epsilon) || !is.finite(step_epsilon)) {
    step_epsilon <- curve_step_epsilon(pr, step_epsilon_abs, step_epsilon_fraction)
  }
  dy <- c(diff(y), NA_real_)
  sign <- rep(NA_integer_, length(dy))
  sign[is.finite(dy) & dy > step_epsilon] <- 1L
  sign[is.finite(dy) & dy < -step_epsilon] <- -1L
  sign[is.finite(dy) & abs(dy) <= step_epsilon] <- 0L
  out <- data.frame(
    O2_pct = curve[[x_col]],
    finite_difference_next = dy,
    local_slope_sign = sign,
    step_epsilon = rep(step_epsilon, length(dy)),
    stringsAsFactors = FALSE
  )
  if (!is.null(id_col) && id_col %in% names(curve)) {
    out <- cbind(data.frame(seed_id = curve[[id_col]], stringsAsFactors = FALSE), out)
  }
  out
}

terminal_plateau_from_signs <- function(signs, plateau_min_points = 3L) {
  signs <- signs[!is.na(signs)]
  if (!length(signs)) return(FALSE)
  last_nz <- suppressWarnings(max(which(signs != 0L)))
  if (!is.finite(last_nz)) return(FALSE)
  (length(signs) - last_nz) >= plateau_min_points
}

classify_o2_ploidy_curve <- function(curve,
                                     value_col = "dominant_mean_ploidy",
                                     x_col = "O2_pct",
                                     id_col = "seed_id",
                                     flat_range_threshold = 0.05,
                                     step_epsilon_abs = 1e-6,
                                     step_epsilon_fraction = 1e-4,
                                     reverse_fraction_tolerance = 0.05,
                                     plateau_min_points = 3L) {
  if (!value_col %in% names(curve)) stop("Missing value column: ", value_col, call. = FALSE)
  if (!x_col %in% names(curve)) stop("Missing x column: ", x_col, call. = FALSE)
  curve <- curve[order(curve[[x_col]]), , drop = FALSE]
  y <- suppressWarnings(as.numeric(curve[[value_col]]))
  finite_y <- y[is.finite(y)]
  if (length(finite_y) < 2L) {
    step_epsilon <- curve_step_epsilon(NA_real_, step_epsilon_abs, step_epsilon_fraction)
    diffs <- finite_diff_curve(curve, step_epsilon, value_col, x_col, id_col, step_epsilon_abs, step_epsilon_fraction)
    summary <- data.frame(
      curve_class = "insufficient_data",
      sign_sequence = "",
      n_sign_changes = NA_integer_,
      step_epsilon = step_epsilon,
      slope_epsilon = step_epsilon,
      flat_range_threshold = flat_range_threshold,
      reverse_fraction_tolerance = reverse_fraction_tolerance,
      ploidy_range = NA_real_,
      net_ploidy_change = NA_real_,
      max_positive_step = NA_real_,
      max_negative_step = NA_real_,
      positive_step_total = NA_real_,
      negative_step_total = NA_real_,
      fraction_positive_steps = NA_real_,
      fraction_negative_steps = NA_real_,
      fraction_zero_steps = NA_real_,
      low_amplitude_curve = NA,
      terminal_plateau = NA,
      terminal_plateau_for_class = NA,
      classification_rule_version = curve_classification_rule_version(),
      stringsAsFactors = FALSE
    )
    return(list(summary = summary, differences = diffs, collapsed_signs = integer(0)))
  }

  ploidy_range <- max(finite_y, na.rm = TRUE) - min(finite_y, na.rm = TRUE)
  finite_idx <- which(is.finite(y))
  net_change <- y[finite_idx[[length(finite_idx)]]] - y[finite_idx[[1L]]]
  step_epsilon <- curve_step_epsilon(ploidy_range, step_epsilon_abs, step_epsilon_fraction)
  diffs <- finite_diff_curve(curve, step_epsilon, value_col, x_col, id_col, step_epsilon_abs, step_epsilon_fraction)
  dy <- diffs$finite_difference_next
  signs <- diffs$local_slope_sign
  step_ok <- is.finite(dy) & !is.na(signs)
  dy <- dy[step_ok]
  signs <- signs[step_ok]
  collapsed <- collapse_signs(signs)
  n_changes <- if (length(collapsed) <= 1L) 0L else length(collapsed) - 1L
  pos_total <- sum(dy[signs > 0], na.rm = TRUE)
  neg_total <- sum(abs(dy[signs < 0]), na.rm = TRUE)
  frac_pos <- if (length(signs)) mean(signs > 0) else NA_real_
  frac_neg <- if (length(signs)) mean(signs < 0) else NA_real_
  frac_zero <- if (length(signs)) mean(signs == 0) else NA_real_
  max_pos <- suppressWarnings(max(dy[signs > 0], na.rm = TRUE))
  max_neg <- suppressWarnings(max(abs(dy[signs < 0]), na.rm = TRUE))
  if (!is.finite(max_pos)) max_pos <- 0
  if (!is.finite(max_neg)) max_neg <- 0
  terminal_plateau <- terminal_plateau_from_signs(signs, plateau_min_points)
  low_amplitude_curve <- is.finite(ploidy_range) && ploidy_range <= flat_range_threshold
  terminal_plateau_for_class <- terminal_plateau && !low_amplitude_curve

  reverse_ok_for_increase <- neg_total <= reverse_fraction_tolerance * max(pos_total, ploidy_range, .Machine$double.eps)
  reverse_ok_for_decrease <- pos_total <= reverse_fraction_tolerance * max(neg_total, ploidy_range, .Machine$double.eps)

  curve_class <- "complex_nonmonotone"
  if (!length(collapsed)) {
    curve_class <- if (is.finite(net_change) && net_change > 0) {
      "monotone_increasing"
    } else if (is.finite(net_change) && net_change < 0) {
      "monotone_decreasing"
    } else {
      "approximately_flat"
    }
  } else if (length(collapsed) == 1L && collapsed[[1L]] == 1L) {
    curve_class <- if (terminal_plateau_for_class) "single_transition_increase_then_plateau" else "monotone_increasing"
  } else if (length(collapsed) == 1L && collapsed[[1L]] == -1L) {
    curve_class <- if (terminal_plateau_for_class) "single_transition_decrease_then_plateau" else "monotone_decreasing"
  } else if (length(collapsed) == 2L && identical(as.integer(collapsed), c(-1L, 1L))) {
    curve_class <- "u_shaped"
  } else if (length(collapsed) == 2L && identical(as.integer(collapsed), c(1L, -1L))) {
    curve_class <- "inverted_u_shaped"
  } else if (is.finite(net_change) && net_change > 0 && pos_total > 0 && reverse_ok_for_increase) {
    curve_class <- "monotone_increasing"
  } else if (is.finite(net_change) && net_change < 0 && neg_total > 0 && reverse_ok_for_decrease) {
    curve_class <- "monotone_decreasing"
  }

  summary <- data.frame(
    curve_class = curve_class,
    sign_sequence = if (length(collapsed)) paste(collapsed, collapse = ",") else "",
    n_sign_changes = n_changes,
    step_epsilon = step_epsilon,
    slope_epsilon = step_epsilon,
    flat_range_threshold = flat_range_threshold,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    ploidy_range = ploidy_range,
    net_ploidy_change = net_change,
    max_positive_step = max_pos,
    max_negative_step = max_neg,
    positive_step_total = pos_total,
    negative_step_total = neg_total,
    fraction_positive_steps = frac_pos,
    fraction_negative_steps = frac_neg,
    fraction_zero_steps = frac_zero,
    low_amplitude_curve = low_amplitude_curve,
    terminal_plateau = terminal_plateau,
    terminal_plateau_for_class = terminal_plateau_for_class,
    classification_rule_version = curve_classification_rule_version(),
    stringsAsFactors = FALSE
  )
  list(summary = summary, differences = diffs, collapsed_signs = collapsed)
}

run_curve_classification_validation <- function() {
  mk <- function(seed, y, gap = 0.02) {
    data.frame(
      seed_id = seed,
      O2_pct = seq_along(y),
      dominant_mean_ploidy = y,
      spectral_gap = rep(gap, length(y)),
      stringsAsFactors = FALSE
    )
  }
  cases <- list(
    flat = list(y = rep(2, 201), expected = "approximately_flat"),
    tiny_increasing = list(y = seq(2, 2.0002, length.out = 201), expected = "monotone_increasing"),
    tiny_u = list(y = c(seq(2.0004, 2, length.out = 101), seq(2.000004, 2.0002, length.out = 100)), expected = "u_shaped"),
    small_increasing = list(y = seq(2, 2.2, length.out = 201), expected = "monotone_increasing"),
    small_decreasing = list(y = seq(2.2, 2, length.out = 201), expected = "monotone_decreasing"),
    increasing = list(y = seq(1, 3, length.out = 201), expected = "monotone_increasing"),
    decreasing = list(y = seq(3, 1, length.out = 201), expected = "monotone_decreasing"),
    increase_plateau = list(y = c(seq(1, 2.2, length.out = 120), rep(2.2, 81)), expected = "single_transition_increase_then_plateau"),
    decrease_plateau = list(y = c(seq(3, 1.8, length.out = 120), rep(1.8, 81)), expected = "single_transition_decrease_then_plateau"),
    u = list(y = c(seq(3, 1.7, length.out = 101), seq(1.72, 3, length.out = 100)), expected = "u_shaped"),
    inverted_u = list(y = c(seq(1, 2.8, length.out = 101), seq(2.78, 1, length.out = 100)), expected = "inverted_u_shaped"),
    complex = list(y = c(seq(1, 2, length.out = 67), seq(1.98, 1.4, length.out = 67), seq(1.42, 2.4, length.out = 67)), expected = "complex_nonmonotone")
  )
  rows <- lapply(names(cases), function(nm) {
    z <- mk(paste0("seed_", nm), cases[[nm]]$y)
    obs <- classify_o2_ploidy_curve(z)$summary$curve_class[[1L]]
    data.frame(
      test_case = nm,
      expected_class = cases[[nm]]$expected,
      observed_class = obs,
      passed = identical(cases[[nm]]$expected, obs),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
