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

smooth_curve_classification_rule_version <- function() {
  "loess_persistent_v1"
}

smooth_curve_values <- function(x, y,
                                span = 0.20,
                                degree = 2L,
                                family = "symmetric",
                                spline_spar = 0.65) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  out <- rep(NA_real_, length(y))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 3L) {
    out[ok] <- y[ok]
    return(out)
  }

  dat <- data.frame(x = x[ok], y = y[ok])
  span <- min(max(as.numeric(span), 0.05), 1)
  degree <- as.integer(degree)
  if (!degree %in% c(0L, 1L, 2L)) degree <- 2L
  fit <- tryCatch(
    suppressWarnings(stats::loess(
      y ~ x,
      data = dat,
      span = span,
      degree = degree,
      family = family,
      control = stats::loess.control(surface = "direct", trace.hat = "approximate")
    )),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    pred <- tryCatch(
      suppressWarnings(as.numeric(stats::predict(fit, newdata = data.frame(x = x[ok])))),
      error = function(e) rep(NA_real_, sum(ok))
    )
    out[ok] <- pred
  }

  if (any(!is.finite(out[ok]))) {
    spl <- tryCatch(
      stats::smooth.spline(dat$x, dat$y, spar = spline_spar),
      error = function(e) NULL
    )
    if (!is.null(spl)) {
      pred <- tryCatch(
        as.numeric(stats::predict(spl, x[ok])$y),
        error = function(e) rep(NA_real_, sum(ok))
      )
      out[ok] <- pred
    }
  }

  bad <- which(ok & !is.finite(out))
  if (length(bad)) out[bad] <- y[bad]
  out
}

persistent_sign_segments <- function(x, y,
                                     signs,
                                     min_segment_points,
                                     min_segment_span,
                                     min_segment_amplitude) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  signs <- as.integer(signs)
  n <- length(x)
  if (n < 2L) {
    return(data.frame())
  }
  signs <- signs[seq_len(n - 1L)]
  signs[is.na(signs)] <- 0L
  runs <- rle(signs)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  rows <- list()
  for (i in seq_along(runs$values)) {
    s <- runs$values[[i]]
    if (!s %in% c(-1L, 1L)) next
    st <- starts[[i]]
    en <- ends[[i]]
    if (en + 1L > n) next
    x_span <- x[[en + 1L]] - x[[st]]
    signed_amplitude <- y[[en + 1L]] - y[[st]]
    amplitude <- abs(signed_amplitude)
    persistent <- is.finite(x_span) && is.finite(amplitude) &&
      (en - st + 1L) >= min_segment_points &&
      x_span >= min_segment_span &&
      amplitude >= min_segment_amplitude
    rows[[length(rows) + 1L]] <- data.frame(
      segment_index = length(rows) + 1L,
      sign = s,
      start_interval = st,
      end_interval = en,
      n_steps = en - st + 1L,
      x_start = x[[st]],
      x_end = x[[en + 1L]],
      x_span = x_span,
      signed_amplitude = signed_amplitude,
      amplitude = amplitude,
      persistent = persistent,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$kept_index <- NA_integer_
  keep <- which(out$persistent)
  if (length(keep)) out$kept_index[keep] <- seq_along(keep)
  out
}

classify_o2_ploidy_curve_smooth <- function(curve,
                                            value_col = "dominant_mean_ploidy",
                                            x_col = "O2_pct",
                                            id_col = "seed_id",
                                            flat_range_threshold = 0.05,
                                            step_epsilon_abs = 1e-6,
                                            step_epsilon_fraction = 1e-4,
                                            reverse_fraction_tolerance = 0.05,
                                            smooth_span = 0.20,
                                            smooth_degree = 2L,
                                            smooth_family = "symmetric",
                                            min_segment_span_fraction = 0.02,
                                            min_segment_amplitude_abs = 0.01,
                                            min_segment_amplitude_fraction = 0.03,
                                            min_segment_points = 3L,
                                            terminal_plateau_span_fraction = 0.10,
                                            terminal_plateau_amplitude_fraction = 0.03) {
  if (!value_col %in% names(curve)) stop("Missing value column: ", value_col, call. = FALSE)
  if (!x_col %in% names(curve)) stop("Missing x column: ", x_col, call. = FALSE)
  curve <- curve[order(curve[[x_col]]), , drop = FALSE]
  x <- suppressWarnings(as.numeric(curve[[x_col]]))
  y <- suppressWarnings(as.numeric(curve[[value_col]]))
  finite_y <- y[is.finite(y)]
  if (length(finite_y) < 2L) {
    step_epsilon <- curve_step_epsilon(NA_real_, step_epsilon_abs, step_epsilon_fraction)
    diffs <- data.frame(
      O2_pct = x,
      raw_value = y,
      fitted_value = y,
      finite_difference_next = c(diff(y), NA_real_),
      fitted_difference_next = c(diff(y), NA_real_),
      local_slope_sign = rep(NA_integer_, length(y)),
      fitted_local_slope_sign = rep(NA_integer_, length(y)),
      step_epsilon = rep(step_epsilon, length(y)),
      stringsAsFactors = FALSE
    )
    if (!is.null(id_col) && id_col %in% names(curve)) {
      diffs <- cbind(data.frame(seed_id = curve[[id_col]], stringsAsFactors = FALSE), diffs)
    }
    summary <- data.frame(
      curve_class = "insufficient_data",
      sign_sequence = "",
      n_sign_changes = NA_integer_,
      classification_rule_version = smooth_curve_classification_rule_version(),
      raw_ploidy_range = NA_real_,
      fitted_ploidy_range = NA_real_,
      net_ploidy_change = NA_real_,
      step_epsilon = step_epsilon,
      flat_range_threshold = flat_range_threshold,
      reverse_fraction_tolerance = reverse_fraction_tolerance,
      smooth_span = smooth_span,
      smooth_degree = smooth_degree,
      min_segment_points = min_segment_points,
      min_segment_span = NA_real_,
      min_segment_amplitude = NA_real_,
      terminal_plateau = NA,
      terminal_plateau_for_class = NA,
      positive_persistent_total = NA_real_,
      negative_persistent_total = NA_real_,
      n_persistent_segments = NA_integer_,
      stringsAsFactors = FALSE
    )
    return(list(summary = summary, differences = diffs, segments = data.frame()))
  }

  fitted <- smooth_curve_values(
    x = x,
    y = y,
    span = smooth_span,
    degree = smooth_degree,
    family = smooth_family
  )
  raw_range <- max(finite_y, na.rm = TRUE) - min(finite_y, na.rm = TRUE)
  finite_fitted <- fitted[is.finite(fitted)]
  fitted_range <- max(finite_fitted, na.rm = TRUE) - min(finite_fitted, na.rm = TRUE)
  if (!is.finite(fitted_range)) fitted_range <- NA_real_
  finite_idx <- which(is.finite(fitted))
  net_change <- fitted[finite_idx[[length(finite_idx)]]] - fitted[finite_idx[[1L]]]
  step_epsilon <- curve_step_epsilon(fitted_range, step_epsilon_abs, step_epsilon_fraction)
  raw_dy <- c(diff(y), NA_real_)
  fitted_dy <- c(diff(fitted), NA_real_)
  raw_sign <- rep(NA_integer_, length(raw_dy))
  fitted_sign <- rep(NA_integer_, length(fitted_dy))
  raw_sign[is.finite(raw_dy) & raw_dy > step_epsilon] <- 1L
  raw_sign[is.finite(raw_dy) & raw_dy < -step_epsilon] <- -1L
  raw_sign[is.finite(raw_dy) & abs(raw_dy) <= step_epsilon] <- 0L
  fitted_sign[is.finite(fitted_dy) & fitted_dy > step_epsilon] <- 1L
  fitted_sign[is.finite(fitted_dy) & fitted_dy < -step_epsilon] <- -1L
  fitted_sign[is.finite(fitted_dy) & abs(fitted_dy) <= step_epsilon] <- 0L

  x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (!is.finite(x_range) || x_range <= 0) x_range <- 1
  min_segment_points_eff <- max(as.integer(min_segment_points), ceiling(min_segment_span_fraction * max(length(x) - 1L, 1L)))
  min_segment_span <- min_segment_span_fraction * x_range
  min_segment_amplitude <- max(min_segment_amplitude_abs, min_segment_amplitude_fraction * fitted_range)
  if (!is.finite(min_segment_amplitude)) min_segment_amplitude <- min_segment_amplitude_abs

  seg <- persistent_sign_segments(
    x = x,
    y = fitted,
    signs = fitted_sign,
    min_segment_points = min_segment_points_eff,
    min_segment_span = min_segment_span,
    min_segment_amplitude = min_segment_amplitude
  )
  kept <- if (nrow(seg)) seg[seg$persistent, , drop = FALSE] else data.frame()
  kept_signs <- if (nrow(kept)) as.integer(kept$sign) else integer(0)
  collapsed <- if (length(kept_signs)) rle(kept_signs)$values else integer(0)
  n_changes <- if (length(collapsed) <= 1L) 0L else length(collapsed) - 1L
  pos_total <- if (nrow(kept)) sum(kept$amplitude[kept$sign > 0], na.rm = TRUE) else 0
  neg_total <- if (nrow(kept)) sum(kept$amplitude[kept$sign < 0], na.rm = TRUE) else 0

  terminal_plateau <- FALSE
  if (nrow(kept)) {
    last_end <- kept$end_interval[[nrow(kept)]]
    anchor <- min(last_end + 1L, length(x))
    trailing_span <- x[[length(x)]] - x[[anchor]]
    trailing_amplitude <- abs(fitted[[length(fitted)]] - fitted[[anchor]])
    terminal_plateau <- is.finite(trailing_span) && is.finite(trailing_amplitude) &&
      trailing_span >= terminal_plateau_span_fraction * x_range &&
      trailing_amplitude <= max(min_segment_amplitude, terminal_plateau_amplitude_fraction * fitted_range)
  }
  low_amplitude_curve <- is.finite(fitted_range) && fitted_range <= flat_range_threshold
  terminal_plateau_for_class <- terminal_plateau && !low_amplitude_curve
  reverse_ok_for_increase <- neg_total <= reverse_fraction_tolerance * max(pos_total, fitted_range, .Machine$double.eps)
  reverse_ok_for_decrease <- pos_total <= reverse_fraction_tolerance * max(neg_total, fitted_range, .Machine$double.eps)

  curve_class <- "complex_nonmonotone"
  if (low_amplitude_curve) {
    curve_class <- "approximately_flat"
  } else if (!length(collapsed)) {
    curve_class <- if (is.finite(net_change) && net_change > flat_range_threshold) {
      "monotone_increasing"
    } else if (is.finite(net_change) && net_change < -flat_range_threshold) {
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

  diffs <- data.frame(
    O2_pct = x,
    raw_value = y,
    fitted_value = fitted,
    finite_difference_next = raw_dy,
    fitted_difference_next = fitted_dy,
    local_slope_sign = raw_sign,
    fitted_local_slope_sign = fitted_sign,
    step_epsilon = rep(step_epsilon, length(y)),
    stringsAsFactors = FALSE
  )
  if (!is.null(id_col) && id_col %in% names(curve)) {
    diffs <- cbind(data.frame(seed_id = curve[[id_col]], stringsAsFactors = FALSE), diffs)
  }
  summary <- data.frame(
    curve_class = curve_class,
    sign_sequence = if (length(collapsed)) paste(collapsed, collapse = ",") else "",
    n_sign_changes = n_changes,
    classification_rule_version = smooth_curve_classification_rule_version(),
    raw_ploidy_range = raw_range,
    fitted_ploidy_range = fitted_range,
    net_ploidy_change = net_change,
    step_epsilon = step_epsilon,
    flat_range_threshold = flat_range_threshold,
    reverse_fraction_tolerance = reverse_fraction_tolerance,
    smooth_span = smooth_span,
    smooth_degree = smooth_degree,
    min_segment_points = min_segment_points_eff,
    min_segment_span = min_segment_span,
    min_segment_amplitude = min_segment_amplitude,
    terminal_plateau = terminal_plateau,
    terminal_plateau_for_class = terminal_plateau_for_class,
    positive_persistent_total = pos_total,
    negative_persistent_total = neg_total,
    n_persistent_segments = nrow(kept),
    stringsAsFactors = FALSE
  )
  list(summary = summary, differences = diffs, segments = seg)
}

run_smooth_curve_classification_validation <- function() {
  mk <- function(seed, y) {
    data.frame(
      seed_id = seed,
      O2_pct = seq(0, 5, length.out = length(y)),
      dominant_mean_ploidy = y,
      stringsAsFactors = FALSE
    )
  }
  y_inc_with_spike <- seq(1, 3, length.out = 201)
  y_inc_with_spike[100] <- y_inc_with_spike[100] - 0.35
  cases <- list(
    flat_noise = list(y = 2 + 0.004 * sin(seq(0, 8 * pi, length.out = 201)), expected = "approximately_flat"),
    increasing_spike = list(y = y_inc_with_spike, expected = "monotone_increasing"),
    increasing = list(y = seq(1, 3, length.out = 201), expected = "monotone_increasing"),
    decreasing = list(y = seq(3, 1, length.out = 201), expected = "monotone_decreasing"),
    u = list(y = c(seq(3, 1.6, length.out = 101), seq(1.62, 3, length.out = 100)), expected = "u_shaped"),
    inverted_u = list(y = c(seq(1, 2.8, length.out = 101), seq(2.78, 1, length.out = 100)), expected = "inverted_u_shaped"),
    increase_plateau = list(y = c(seq(1, 2.2, length.out = 115), rep(2.2, 86)), expected = "single_transition_increase_then_plateau")
  )
  rows <- lapply(names(cases), function(nm) {
    z <- mk(paste0("seed_", nm), cases[[nm]]$y)
    obs <- classify_o2_ploidy_curve_smooth(z)$summary$curve_class[[1L]]
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
