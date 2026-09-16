#!/usr/bin/env Rscript

# Full-range extension of Supplementary Figure 6-14. This file is sourced
# after figure6_net_growth_q10.R and deliberately overrides only the range,
# storage, publication, and aggregation contracts. The short-range run and
# its current pointer remain untouched.

options(stringsAsFactors = FALSE, warn = 1)

f6ng_profile <- function() {
  "figure6_population_net_live_growth_full_range_q10_v1"
}
f6ng_initial_ploidy <- function() c(2, 3, 4, 5, 6)
f6ng_days <- function(smoke = FALSE) {
  if (isTRUE(smoke)) 0:10 else 0:10000
}
f6ng_p_values <- function(smoke = FALSE) {
  if (isTRUE(smoke)) f6ft_p_values()[1:2] else f6ft_p_values()
}
f6ng_o2 <- function(context, smoke = FALSE) {
  if (identical(context, "in vivo")) {
    if (isTRUE(smoke)) c(0, 2.5, 5) else seq(0, 5, length.out = 201L)
  } else if (identical(context, "in vitro")) {
    if (isTRUE(smoke)) c(0, 10, 20) else seq(0, 20, by = 0.1)
  } else {
    stop("Unknown Figure 6 net-growth model context: ", context)
  }
}

f6ng_paths <- function(paths, run_id = NULL, create = FALSE) {
  base_name <- "net_growth_full_range_q10_runs"
  base <- file.path(paths$figure6, base_name)
  current <- file.path(
    paths$figure6, "net_growth_full_range_q10_current.tsv"
  )
  if (is.null(run_id)) {
    f6r_require_files(current, "current full-range Figure 6 net-growth pointer")
    pointer <- f6r_read_tsv(current)
    required <- c("run_id", "relative_run_path", "profile", "fingerprint")
    if (nrow(pointer) != 1L || !all(required %in% names(pointer)) ||
        !identical(as.character(pointer$profile[[1L]]), f6ng_profile())) {
      stop("Malformed full-range Figure 6 net-growth pointer: ", current)
    }
    run_id <- f6ft_sanitize_run_id(pointer$run_id[[1L]])
    run_root <- normalizePath(
      file.path(paths$figure6, pointer$relative_run_path[[1L]]),
      mustWork = TRUE
    )
    expected_root <- normalizePath(file.path(base, run_id), mustWork = TRUE)
    if (!identical(run_root, expected_root)) {
      stop("Full-range net-growth pointer escaped its dedicated run root.")
    }
  } else {
    run_id <- f6ft_sanitize_run_id(run_id)
    run_root <- file.path(base, run_id)
  }
  out <- list(
    run_id = run_id, base_name = base_name, base = base,
    run_root = run_root, cache = file.path(run_root, "task_cache"),
    current = current, rendered = file.path(run_root, "rendered")
  )
  if (isTRUE(create)) {
    invisible(lapply(
      unname(unlist(out[c("base", "run_root", "cache", "rendered")])),
      dir.create, recursive = TRUE, showWarnings = FALSE
    ))
  }
  out
}

f6ng_fingerprint <- function(paths, source, propagator, smoke) {
  inputs <- c(
    source$files,
    file.path(paths$code, "util", "analysis", "figure6_net_growth_q10.R"),
    file.path(
      paths$code, "util", "analysis",
      "figure6_net_growth_full_range_q10.R"
    ),
    propagator
  )
  paste(
    f6ng_profile(), paste0("model=", source$model_fingerprint),
    paste0("source_run=", source$run_paths$run_id),
    paste0("inputs=", paste(unname(tools::md5sum(inputs)), collapse = ":")),
    paste0("initial=", paste(f6ng_initial_ploidy(), collapse = ",")),
    paste0("p=", paste(f6ng_p_values(smoke), collapse = ",")),
    paste0("day=", paste(range(f6ng_days(smoke)), collapse = ":"),
           ":", length(f6ng_days(smoke))),
    paste0("invivo_o2=", paste(range(f6ng_o2("in vivo", smoke)), collapse = ":"),
           ":", length(f6ng_o2("in vivo", smoke))),
    paste0("invitro_o2=", paste(range(f6ng_o2("in vitro", smoke)), collapse = ":"),
           ":", length(f6ng_o2("in vitro", smoke))),
    paste0("replicates=", source$config$replicates),
    paste0("master_seed=", source$config$master_seed), sep = "|"
  )
}

f6ng_short_run <- function(paths) {
  current <- file.path(paths$figure6, "net_growth_q10_current.tsv")
  f6r_require_files(current, "short-range Figure 6 net-growth pointer")
  pointer <- f6r_read_tsv(current)
  if (nrow(pointer) != 1L ||
      !all(c("run_id", "relative_run_path") %in% names(pointer))) {
    stop("Malformed short-range Figure 6 net-growth pointer: ", current)
  }
  list(
    run_id = as.character(pointer$run_id[[1L]]),
    run_root = normalizePath(
      file.path(paths$figure6, pointer$relative_run_path[[1L]]), mustWork = TRUE
    )
  )
}

f6ng_short_panel <- function(paths, context, family) {
  run <- f6ng_short_run(paths)
  path <- file.path(run$run_root, f6ng_panel_filename(context, family))
  f6r_require_files(path, paste(context, family, "short-range net-growth panel"))
  object <- readRDS(path)
  if (!identical(object$model_context, context) ||
      !identical(object$pair_label, family)) {
    stop("Unexpected short-range net-growth panel: ", path)
  }
  object
}

f6ng_overlap_row <- function(paths, full) {
    context <- full$model_context
    family <- full$pair_label
    short_run <- f6ng_short_run(paths)
    short <- f6ng_short_panel(paths, context, family)
    initial <- match(short$initial_ploidy, full$initial_ploidy)
    day <- match(short$day_values, full$day_values)
    oxygen <- match(
      sprintf("%.12f", short$o2_values), sprintf("%.12f", full$o2_values)
    )
    p <- match(
      sprintf("%.12f", short$p_misseg), sprintf("%.12f", full$p_misseg)
    )
    if (anyNA(c(initial, day, oxygen, p))) {
      stop("Short-range net-growth grid is not nested in the full-range grid.")
    }
    subset <- full$mean_net_growth_rate[
      initial, day, oxygen, p, drop = FALSE
    ]
    delta <- abs(subset - short$mean_net_growth_rate)
    data.frame(
      model_context = context, pair_label = family,
      short_run_id = short_run$run_id,
      n_compared = length(delta), maximum_absolute_difference = max(delta),
      tolerance = 1e-12, passed = all(is.finite(delta)) && max(delta) <= 1e-12,
      stringsAsFactors = FALSE
    )
}

# The full-range panel files retain only the metric needed for drawing and its
# Monte Carlo error. Mean-ploidy replay is checked against the established
# full-range Figure 6 source before being discarded, avoiding four redundant
# 400-MB arrays per panel.
f6ng_aggregate <- function(tasks, source, run_paths, fingerprint, smoke = FALSE) {
  task_qc <- do.call(rbind, lapply(
    tasks$cache_path, function(path) readRDS(path)$qc
  ))
  task_qc_path <- f6ft_atomic_write_tsv(
    task_qc, file.path(run_paths$run_root, "net_growth_task_validation.tsv")
  )
  if (!all(task_qc$passed)) stop("One or more net-growth tasks failed validation.")

  days <- f6ng_days(smoke)
  p_values <- f6ng_p_values(smoke)
  panel_qc <- list()
  panel_paths <- character()
  observed_range <- c(Inf, -Inf)
  overlap_rows <- list()
  panel_index <- 0L
  for (context in f6ft_context_levels()) for (family in f6ft_family_levels()) {
    oxygen <- f6ng_o2(context, smoke)
    shape <- c(
      length(f6ng_initial_ploidy()), length(days), length(oxygen),
      length(p_values)
    )
    dimension_names <- list(
      initial_ploidy = paste0(f6ng_initial_ploidy(), "N"),
      day = as.character(days), O2_pct = as.character(oxygen),
      p_misseg = format(p_values, scientific = FALSE, trim = TRUE)
    )
    rate <- ploidy <- mcse <- array(
      NA_real_, dim = shape, dimnames = dimension_names
    )
    selected <- tasks[
      tasks$model_context == context & tasks$pair_label == family,
      , drop = FALSE
    ]
    for (task_index in seq_len(nrow(selected))) {
      task <- selected[task_index, , drop = FALSE]
      cached <- readRDS(task$cache_path[[1L]])
      if (!identical(cached$fingerprint, fingerprint)) {
        stop("Full-range net-growth aggregate found an incompatible task cache.")
      }
      pi <- match(as.numeric(task$p_misseg[[1L]]), p_values)
      rate[, , cached$oxygen_index, pi] <- cached$mean_net_growth_rate
      ploidy[, , cached$oxygen_index, pi] <- cached$mean_ploidy_replay
      mcse[, , cached$oxygen_index, pi] <- cached$stochastic_mcse
    }
    if (any(!is.finite(rate)) || any(!is.finite(ploidy)) || any(!is.finite(mcse))) {
      stop("Incomplete full-range net-growth aggregate for ", context, " ", family)
    }
    represented <- if (isTRUE(smoke)) {
      unique(selected$represented_optimizer_endpoint)[[1L]]
    } else 50L
    replay_error <- 0
    source_mask_match <- TRUE
    replay_tolerance <- NA_real_
    if (!isTRUE(smoke)) {
      source_object <- f6g_read_panel(
        source$run_paths, context, f6ng_mode(context), family
      )
      source_ploidy <- f6ng_source_subset(
        source_object, f6ng_initial_ploidy(), days, oxygen, p_values
      )
      replay_error <- max(abs(ploidy - source_ploidy), na.rm = TRUE)
      source_mask_match <- identical(is.na(ploidy), is.na(source_ploidy))
      replay_tolerance <- if (identical(context, "in vivo")) 1e-10 else 0.05
      rm(source_object, source_ploidy)
    }
    maximum_mcse <- max(mcse)
    mcse_target <- if (identical(context, "in vivo")) 0 else
      as.numeric(source$config$mcse_target_N %||% 0.01)
    object <- list(
      profile = f6ng_profile(), fingerprint = fingerprint,
      source_full_range_run_id = source$run_paths$run_id,
      model_context = context, propagation_mode = f6ng_mode(context),
      pair_label = family, initial_ploidy = f6ng_initial_ploidy(),
      day_values = days, o2_values = oxygen, p_misseg = p_values,
      optimizer_endpoint_weight = rep(represented, length(p_values)),
      mean_net_growth_rate = rate, stochastic_mcse = mcse,
      metric_definition = paste0(
        "sum_N f_N(t)*colSums(M)_N = 1^T M x(t)/1^T x(t); day^-1; ",
        "passage-day state is post-sampling and dilution is excluded"
      )
    )
    panel_path <- file.path(
      run_paths$run_root, f6ng_panel_filename(context, family)
    )
    f6ft_atomic_save_rds(object, panel_path, compress = "gzip")
    observed_range <- range(
      observed_range, object$mean_net_growth_rate, finite = TRUE
    )
    if (!isTRUE(smoke)) {
      overlap_rows[[panel_index + 1L]] <- f6ng_overlap_row(
        f6r_paths(), object
      )
    }
    panel_paths <- c(panel_paths, panel_path)
    panel_index <- panel_index + 1L
    panel_qc[[panel_index]] <- data.frame(
      model_context = context, propagation_mode = f6ng_mode(context),
      pair_label = family, n_initial = length(f6ng_initial_ploidy()),
      n_day = length(days), n_o2 = length(oxygen), n_p = length(p_values),
      optimizer_endpoint_weight = represented,
      minimum_net_growth_per_day = min(rate),
      maximum_net_growth_per_day = max(rate),
      maximum_growth_mcse_per_day = maximum_mcse,
      growth_mcse_target_per_day = mcse_target,
      maximum_ploidy_replay_error = replay_error,
      ploidy_replay_tolerance = replay_tolerance,
      source_missing_mask_match = source_mask_match,
      validation_policy = if (isTRUE(smoke)) "smoke kernel identity" else
        if (identical(context, "in vivo")) "exact deterministic replay" else
          "fixed-200 stochastic replay versus established full-range source",
      passed = all(is.finite(rate)) && source_mask_match &&
        (isTRUE(smoke) || replay_error <= replay_tolerance) &&
        (identical(context, "in vivo") || maximum_mcse <= mcse_target),
      stringsAsFactors = FALSE
    )
    rm(rate, ploidy, mcse, object)
    invisible(gc(FALSE))
  }
  panel_qc <- do.call(rbind, panel_qc)
  panel_qc_path <- f6ft_atomic_write_tsv(
    panel_qc, file.path(run_paths$run_root, "net_growth_panel_validation.tsv")
  )
  if (!all(panel_qc$passed)) stop("Full-range net-growth panel validation failed.")

  observed <- observed_range
  bound <- if (observed[[1L]] < 0 && observed[[2L]] > 0) {
    max(abs(observed))
  } else max(abs(observed))
  color_contract <- data.frame(
    scale = "signed_pseudo_log10", sigma_per_day = 0.01,
    observed_minimum_per_day = observed[[1L]],
    observed_maximum_per_day = observed[[2L]],
    displayed_minimum_per_day = -bound,
    displayed_maximum_per_day = bound,
    shared_across = "in vivo C01/C02 and in vitro C01/C02",
    stringsAsFactors = FALSE
  )
  color_path <- f6ft_atomic_write_tsv(
    color_contract,
    file.path(run_paths$run_root, "net_growth_full_range_color_contract.tsv")
  )
  overlap_path <- NA_character_
  if (!isTRUE(smoke)) {
    overlap_validation <- do.call(rbind, overlap_rows)
    overlap_path <- f6ft_atomic_write_tsv(
      overlap_validation,
      file.path(
        run_paths$run_root,
        "net_growth_short_range_overlap_validation.tsv"
      )
    )
    if (!all(overlap_validation$passed)) {
      stop("Full-range net-growth values do not exactly replay the short-range run.")
    }
  }
  list(
    task_qc = task_qc_path, panel_qc = panel_qc_path,
    panels = normalizePath(panel_paths, mustWork = TRUE),
    overlap = overlap_path, color_contract = color_path
  )
}

f6ng_publish_current <- function(paths, run_paths, fingerprint) {
  pointer <- data.frame(
    run_id = run_paths$run_id,
    relative_run_path = file.path(run_paths$base_name, run_paths$run_id),
    profile = f6ng_profile(), fingerprint = fingerprint,
    published_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  f6ft_atomic_write_tsv(pointer, run_paths$current)
}
