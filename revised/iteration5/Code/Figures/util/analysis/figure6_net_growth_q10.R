#!/usr/bin/env Rscript

# Supplementary Figure 6-14: population-weighted net-live growth rate under
# the same fixed-p_misseg, finite-time, and stochastic-passage contracts used
# by Figure 6B. Existing Figure 6 panel objects are read-only validation inputs.

options(stringsAsFactors = FALSE, warn = 1)

f6ng_profile <- function() "figure6_population_net_live_growth_q10_v1"
f6ng_initial_ploidy <- function() c(2, 4)
f6ng_days <- function(smoke = FALSE) if (isTRUE(smoke)) 0:20 else 0:1000
f6ng_p_values <- function(smoke = FALSE) {
  if (isTRUE(smoke)) f6ft_p_values()[1:2] else f6ft_p_values()
}
f6ng_o2 <- function(context, smoke = FALSE) {
  if (isTRUE(smoke)) return(c(0, 0.5, 2))
  values <- f6g_o2(context)
  values[values >= -1e-12 & values <= 2 + 1e-12]
}
f6ng_mode <- function(context) {
  if (identical(context, "in vivo")) "continuous" else "passage"
}

f6ng_paths <- function(paths, run_id = NULL, create = FALSE) {
  base <- file.path(paths$figure6, "net_growth_q10_runs")
  current <- file.path(paths$figure6, "net_growth_q10_current.tsv")
  if (is.null(run_id)) {
    f6r_require_files(current, "current Figure 6 net-growth run pointer")
    pointer <- f6r_read_tsv(current)
    required <- c("run_id", "relative_run_path", "profile", "fingerprint")
    if (nrow(pointer) != 1L || !all(required %in% names(pointer))) {
      stop("Malformed Figure 6 net-growth pointer: ", current)
    }
    run_id <- f6ft_sanitize_run_id(pointer$run_id[[1L]])
    run_root <- normalizePath(
      file.path(paths$figure6, pointer$relative_run_path[[1L]]), mustWork = TRUE
    )
  } else {
    run_id <- f6ft_sanitize_run_id(run_id)
    run_root <- file.path(base, run_id)
  }
  out <- list(
    run_id = run_id, base = base, run_root = run_root,
    cache = file.path(run_root, "task_cache"),
    current = current,
    rendered = file.path(run_root, "rendered")
  )
  if (isTRUE(create)) {
    invisible(lapply(
      unname(unlist(out[c("base", "run_root", "cache", "rendered")])),
      dir.create, recursive = TRUE, showWarnings = FALSE
    ))
  }
  out
}

f6ng_panel_filename <- function(context, family) {
  context_key <- if (identical(context, "in vivo")) "invivo" else "invitro"
  paste0(
    "net_growth_panel_", context_key, "_", f6ng_mode(context), "_",
    tolower(family), ".rds"
  )
}

f6ng_read_panel <- function(run_paths, context, family) {
  path <- file.path(run_paths$run_root, f6ng_panel_filename(context, family))
  f6r_require_files(path, paste(context, family, "net-growth panel"))
  object <- readRDS(path)
  expected <- c(
    identical(object$profile, f6ng_profile()),
    identical(object$model_context, context),
    identical(object$propagation_mode, f6ng_mode(context)),
    identical(object$pair_label, family)
  )
  if (!all(expected)) stop("Unexpected Figure 6 net-growth panel: ", path)
  object
}

f6ng_source_run <- function(paths) {
  override <- trimws(Sys.getenv("FIGURE6_FULL_RANGE_SOURCE_RUN_ROOT", ""))
  if (!nzchar(override)) return(f6g_paths(paths, run_id = NULL, create = FALSE))
  run_root <- normalizePath(override, mustWork = TRUE)
  allowed_root <- normalizePath(
    file.path(paths$root, "data", "Figures"), mustWork = TRUE
  )
  if (!startsWith(run_root, paste0(allowed_root, .Platform$file.sep))) {
    stop("Full-range source override must remain inside iteration4/data/Figures.")
  }
  current <- file.path(
    dirname(dirname(run_root)), "finite_time_full_q10_current.tsv"
  )
  f6r_require_files(current, "full-range source override pointer")
  pointer <- f6r_read_tsv(current)
  if (nrow(pointer) != 1L ||
      !identical(as.character(pointer$run_id[[1L]]), basename(run_root))) {
    stop("Full-range source override does not match its current pointer.")
  }
  list(
    run_id = basename(run_root), run_root = run_root,
    current = normalizePath(current, mustWork = TRUE),
    rendered = file.path(run_root, "rendered")
  )
}

f6ng_source_bundle <- function(paths) {
  source_run <- f6ng_source_run(paths)
  files <- c(
    endpoints = file.path(source_run$run_root, "q10_unique_endpoint_manifest.tsv"),
    expanded = file.path(source_run$run_root, "q10_optimizer_seed_manifest.tsv"),
    endpoint_qc = file.path(source_run$run_root, "q10_endpoint_manifest_validation.tsv"),
    passage_rule = file.path(source_run$run_root, "canonical_passage_rule.tsv"),
    stochastic_config = file.path(source_run$run_root, "stochastic_config.rds"),
    rng_catalog = file.path(source_run$run_root, "stochastic_rng_stream_catalog.rds")
  )
  f6r_require_files(files, "current Figure 6 source run")
  endpoints <- f6r_read_tsv(files[["endpoints"]])
  expanded <- f6r_read_tsv(files[["expanded"]])
  endpoint_qc <- f6r_read_tsv(files[["endpoint_qc"]])
  passage_rule <- f6r_read_tsv(files[["passage_rule"]])
  config <- readRDS(files[["stochastic_config"]])
  catalog <- readRDS(files[["rng_catalog"]])
  if (!all(endpoint_qc$passed) || nrow(endpoint_qc) != 4L ||
      !all(endpoint_qc$n_optimizer_endpoint == 50L)) {
    stop("The current Figure 6 endpoint manifest is not a validated 50-endpoint ensemble.")
  }
  if (!identical(as.integer(config$replicates), 200L) ||
      !identical(as.character(config$allocation), "fixed") ||
      !identical(as.integer(config$master_seed), 20260907L)) {
    stop("Supplementary Figure 6-14 requires the current fixed 200-repeat RNG contract.")
  }
  if (!identical(
      unname(tools::md5sum(files[["rng_catalog"]])),
      as.character(config$stream_catalog_md5)
  )) stop("Figure 6 RNG catalog checksum differs from its recorded configuration.")
  identity_columns <- c("pair_label", "seed_number", "parameter_signature")
  observed_identity <- expanded[
    expanded$model_context == "in vitro", identity_columns, drop = FALSE
  ]
  expected_identity <- catalog$endpoint_identity[, identity_columns, drop = FALSE]
  rownames(observed_identity) <- rownames(expected_identity) <- NULL
  identity_key <- function(x) paste(
    as.character(x$pair_label), as.integer(x$seed_number),
    as.character(x$parameter_signature), sep = "|"
  )
  if (!identical(identity_key(observed_identity), identity_key(expected_identity))) {
    stop("Figure 6 RNG catalog endpoint identities do not match the source run.")
  }
  pointer <- f6r_read_tsv(source_run$current)
  model_fingerprint <- f6r_model_source_fingerprint(paths)
  if (!grepl(paste0("model=", model_fingerprint), pointer$fingerprint[[1L]], fixed = TRUE)) {
    stop("Current external model code does not match the Figure 6 source-run fingerprint.")
  }
  list(
    run_paths = source_run, files = files, endpoints = endpoints,
    expanded = expanded, endpoint_qc = endpoint_qc, rule = passage_rule,
    config = config, catalog = catalog, pointer = pointer,
    model_fingerprint = model_fingerprint
  )
}

f6ng_load_propagator <- function(paths) {
  f6r_require_packages("Rcpp")
  cpp <- file.path(
    paths$code, "util", "analysis", "figure6_net_growth_propagator.cpp"
  )
  f6r_require_files(cpp, "Figure 6 net-growth propagator")
  cache <- file.path(paths$figure6, ".rcpp_cache_net_growth")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(cpp, rebuild = FALSE, cacheDir = cache, showOutput = FALSE)
  required <- c(
    "f6ng_propagate_continuous_cpp", "f6ng_propagate_stochastic_cpp"
  )
  if (!all(vapply(required, exists, logical(1L), mode = "function", inherits = TRUE))) {
    stop("Figure 6 net-growth C++ functions did not load.")
  }
  normalizePath(cpp, mustWork = TRUE)
}

f6ng_fingerprint <- function(paths, source, propagator, smoke) {
  inputs <- c(
    source$files,
    file.path(paths$code, "util", "analysis", "figure6_net_growth_q10.R"),
    propagator
  )
  paste(
    f6ng_profile(),
    paste0("model=", source$model_fingerprint),
    paste0("source_run=", source$run_paths$run_id),
    paste0("inputs=", paste(unname(tools::md5sum(inputs)), collapse = ":")),
    paste0("initial=", paste(f6ng_initial_ploidy(), collapse = ",")),
    paste0("p=", paste(f6ng_p_values(smoke), collapse = ",")),
    paste0("day=", paste(range(f6ng_days(smoke)), collapse = ":")),
    paste0("replicates=", source$config$replicates),
    paste0("master_seed=", source$config$master_seed),
    sep = "|"
  )
}

f6ng_task_manifest <- function(
    source, run_paths, smoke = FALSE, o2_chunk_size = 1L
) {
  o2_chunk_size <- as.integer(o2_chunk_size)
  if (is.na(o2_chunk_size) || o2_chunk_size < 1L) {
    stop("o2_chunk_size must be a positive integer.")
  }
  rows <- list()
  index <- 0L
  for (context in f6ft_context_levels()) {
    oxygen <- f6ng_o2(context, smoke)
    chunks <- split(
      seq_along(oxygen), ceiling(seq_along(oxygen) / o2_chunk_size)
    )
    for (family in f6ft_family_levels()) {
      endpoints <- source$endpoints[
        source$endpoints$model_context == context &
          source$endpoints$pair_label == family, , drop = FALSE
      ]
      endpoints <- endpoints[order(
        endpoints$representative_objective_rank,
        endpoints$representative_seed_number
      ), , drop = FALSE]
      if (isTRUE(smoke)) endpoints <- endpoints[1L, , drop = FALSE]
      for (p in f6ng_p_values(smoke)) for (chunk in seq_along(chunks)) {
        index <- index + 1L
        oi <- chunks[[chunk]]
        rows[[index]] <- data.frame(
          task_id = sprintf("NG%04d", index), model_context = context,
          propagation_mode = f6ng_mode(context), pair_label = family,
          p_misseg = p, o2_chunk_index = chunk,
          o2_index_start = min(oi), o2_index_end = max(oi),
          endpoint_indices = paste(endpoints$endpoint_index, collapse = ","),
          n_unique_endpoint = nrow(endpoints),
          represented_optimizer_endpoint = sum(endpoints$endpoint_multiplicity_q10),
          cache_path = file.path(
            run_paths$cache, sprintf("net_growth_task_%04d.rds", index)
          ), stringsAsFactors = FALSE
        )
      }
    }
  }
  tasks <- do.call(rbind, rows)
  f6ft_atomic_write_tsv(
    tasks, file.path(run_paths$run_root, "net_growth_task_manifest.tsv")
  )
  tasks
}

f6ng_moments <- function(trajectories, multiplicity, replicates) {
  means <- variances <- matrix(0, multiplicity, ncol(trajectories))
  for (seed_index in seq_len(multiplicity)) {
    rows <- (seed_index - 1L) * replicates + seq_len(replicates)
    values <- trajectories[rows, , drop = FALSE]
    means[seed_index, ] <- colMeans(values)
    variances[seed_index, ] <- colSums(
      (values - rep(means[seed_index, ], each = replicates))^2
    ) / (replicates - 1L)
  }
  list(
    sum = colSums(means),
    sum_squared_endpoint_mean = colSums(means^2),
    sum_within_variance = colSums(variances),
    sum_mc_variance = colSums(variances) / replicates
  )
}

f6ng_compute_task <- function(
    task, source, objective_bundle, contexts, paths, run_paths,
    fingerprint, smoke = FALSE
) {
  output <- task$cache_path[[1L]]
  if (file.exists(output)) {
    cached <- readRDS(output)
    if (!identical(cached$fingerprint, fingerprint) || !isTRUE(cached$qc$passed)) {
      stop("Incompatible net-growth task cache: ", output)
    }
    return(cached$qc)
  }
  f6r_load_response_engine(paths)
  if (!exists("f6ng_propagate_continuous_cpp", mode = "function", inherits = TRUE)) {
    f6ng_load_propagator(paths)
  }
  endpoint_indices <- as.integer(strsplit(
    task$endpoint_indices[[1L]], ",", fixed = TRUE
  )[[1L]])
  endpoints <- source$endpoints[
    match(endpoint_indices, source$endpoints$endpoint_index), , drop = FALSE
  ]
  if (anyNA(endpoints$endpoint_index)) stop("Net-growth endpoint lookup failed.")
  context <- task$model_context[[1L]]
  oxygen_all <- f6ng_o2(context, smoke)
  oxygen_index <- seq.int(
    as.integer(task$o2_index_start[[1L]]),
    as.integer(task$o2_index_end[[1L]])
  )
  oxygen <- oxygen_all[oxygen_index]
  days <- f6ng_days(smoke)
  initial_values <- f6ng_initial_ploidy()
  shape <- c(length(initial_values), length(days), length(oxygen))
  empty <- array(0, dim = shape)
  checkpoint <- paste0(output, ".checkpoint.rds")
  work <- list(
    fingerprint = fingerprint, next_o2 = 1L, next_endpoint = 1L,
    net_growth_sum = empty, ploidy_sum = empty,
    endpoint_square_sum = empty, within_variance_sum = empty,
    mc_variance_sum = empty, formula_error = 0,
    maximum_rate_identity_error = 0, passage_count = 0
  )
  if (file.exists(checkpoint)) {
    work <- readRDS(checkpoint)
    if (!identical(work$fingerprint, fingerprint)) {
      stop("Incompatible net-growth task checkpoint: ", checkpoint)
    }
  }
  p <- as.numeric(task$p_misseg[[1L]])
  for (o2_index in seq_along(oxygen)) {
    if (o2_index < work$next_o2) next
    for (endpoint_index in seq_len(nrow(endpoints))) {
      if (endpoint_index < work$next_endpoint) next
      endpoint <- endpoints[endpoint_index, , drop = FALSE]
      prepared <- f6ft_prepare_endpoint(endpoint, objective_bundle, contexts)
      forced <- figure6_force_p_misseg(prepared$run_params, p)
      formula <- figure6_p_misseg_formula_qc(
        prepared$run_params, p, o2_values = oxygen[[o2_index]]
      )
      work$formula_error <- max(
        work$formula_error, formula$maximum_direct_formula_error
      )
      fixed <- fixo2_fixed_matrix(
        globalenv(), prepared$config, forced, O2 = oxygen[[o2_index]]
      )
      unit <- prepared$config$N_UNIT %||% 22L
      initial <- f6ft_initial_matrix(
        fixed$ngrid, unit, initial_ploidy = initial_values
      )
      step <- as.matrix(Matrix::expm(fixed$M))
      net_live_rate <- as.numeric(colSums(fixed$M))
      probe <- initial[, 1L]
      identity_error <- abs(
        sum(probe * net_live_rate) - sum(as.numeric(fixed$M %*% probe))
      )
      work$maximum_rate_identity_error <- max(
        work$maximum_rate_identity_error, identity_error
      )
      weight <- as.integer(endpoint$endpoint_multiplicity_q10[[1L]])
      if (identical(context, "in vivo")) {
        response <- f6ng_propagate_continuous_cpp(
          step, initial, as.numeric(fixed$ngrid) / unit,
          net_live_rate, max(days)
        )
        work$net_growth_sum[, , o2_index] <-
          work$net_growth_sum[, , o2_index] +
          weight * response$net_growth_rate
        work$ploidy_sum[, , o2_index] <-
          work$ploidy_sum[, , o2_index] + weight * response$mean_ploidy
        work$endpoint_square_sum[, , o2_index] <-
          work$endpoint_square_sum[, , o2_index] +
          weight * response$net_growth_rate^2
      } else {
        seeds <- as.integer(strsplit(
          endpoint$represented_seed_numbers[[1L]], ",", fixed = TRUE
        )[[1L]])
        if (length(seeds) != weight) {
          stop("Endpoint multiplicity and represented RNG seeds disagree.")
        }
        replicates <- as.integer(source$config$replicates)
        for (initial_index in seq_along(initial_values)) {
          rng <- f6s_streams(
            source$catalog, endpoint$pair_label[[1L]], seeds,
            oxygen[[o2_index]], p, initial_values[[initial_index]], replicates
          )
          state <- initial[, rep(initial_index, weight * replicates), drop = FALSE]
          response <- f6ng_propagate_stochastic_cpp(
            step, state, as.numeric(fixed$ngrid) / unit, net_live_rate,
            max(days), source$rule$seed_cells[[1L]],
            source$rule$target_live_cells[[1L]], rng,
            rep(log(source$rule$seed_cells[[1L]]), weight * replicates), 0L
          )
          rate_moments <- f6ng_moments(
            response$net_growth_rate, weight, replicates
          )
          ploidy_moments <- f6ng_moments(
            response$mean_ploidy, weight, replicates
          )
          work$net_growth_sum[initial_index, , o2_index] <-
            work$net_growth_sum[initial_index, , o2_index] + rate_moments$sum
          work$ploidy_sum[initial_index, , o2_index] <-
            work$ploidy_sum[initial_index, , o2_index] + ploidy_moments$sum
          work$endpoint_square_sum[initial_index, , o2_index] <-
            work$endpoint_square_sum[initial_index, , o2_index] +
            rate_moments$sum_squared_endpoint_mean
          work$within_variance_sum[initial_index, , o2_index] <-
            work$within_variance_sum[initial_index, , o2_index] +
            rate_moments$sum_within_variance
          work$mc_variance_sum[initial_index, , o2_index] <-
            work$mc_variance_sum[initial_index, , o2_index] +
            rate_moments$sum_mc_variance
          work$passage_count <- work$passage_count + sum(response$passage_count)
        }
      }
      work$next_endpoint <- endpoint_index + 1L
      f6ft_atomic_save_rds(work, checkpoint, compress = FALSE)
    }
    work$next_o2 <- o2_index + 1L
    work$next_endpoint <- 1L
    f6ft_atomic_save_rds(work, checkpoint, compress = FALSE)
  }
  represented <- sum(endpoints$endpoint_multiplicity_q10)
  mean_rate <- work$net_growth_sum / represented
  mean_ploidy <- work$ploidy_sum / represented
  mcse <- sqrt(pmax(work$mc_variance_sum, 0)) / represented
  between_variance <- pmax(
    (work$endpoint_square_sum - represented * mean_rate^2) /
      max(1, represented - 1L), 0
  )
  day0_error <- max(abs(mean_ploidy[, 1L, ] - initial_values))
  qc <- data.frame(
    task_id = task$task_id[[1L]], model_context = context,
    propagation_mode = task$propagation_mode[[1L]],
    pair_label = task$pair_label[[1L]], p_misseg = p,
    o2_index_start = min(oxygen_index), o2_index_end = max(oxygen_index),
    n_unique_endpoint = nrow(endpoints),
    represented_optimizer_endpoint = represented,
    n_stochastic_replicate = if (context == "in vitro") {
      as.integer(source$config$replicates)
    } else 0L,
    maximum_day0_ploidy_error = day0_error,
    maximum_formula_error = work$formula_error,
    maximum_rate_identity_error = work$maximum_rate_identity_error,
    maximum_growth_mcse_per_day = max(mcse),
    passage_count = work$passage_count,
    all_finite = all(is.finite(mean_rate)) && all(is.finite(mean_ploidy)),
    passed = all(is.finite(mean_rate)) && all(is.finite(mean_ploidy)) &&
      day0_error <= 1e-10 && work$formula_error <= 1e-12 &&
      work$maximum_rate_identity_error <= 1e-12,
    stringsAsFactors = FALSE
  )
  f6ft_atomic_save_rds(list(
    profile = f6ng_profile(), fingerprint = fingerprint, task = task,
    oxygen_index = oxygen_index, o2_values = oxygen, day_values = days,
    initial_ploidy = initial_values,
    represented_optimizer_endpoint = represented,
    mean_net_growth_rate = mean_rate, mean_ploidy_replay = mean_ploidy,
    stochastic_mcse = mcse, between_endpoint_variance = between_variance,
    qc = qc
  ), output, compress = "gzip")
  if (file.exists(checkpoint)) unlink(checkpoint)
  qc
}

f6ng_source_subset <- function(source_object, initial, days, oxygen, p_values) {
  initial_index <- match(initial, source_object$initial_ploidy)
  day_index <- match(days, source_object$day_values)
  oxygen_index <- match(
    sprintf("%.12f", oxygen), sprintf("%.12f", source_object$o2_values)
  )
  p_index <- match(
    sprintf("%.12f", p_values), sprintf("%.12f", source_object$p_misseg)
  )
  if (anyNA(c(initial_index, day_index, oxygen_index, p_index))) {
    stop("Requested net-growth validation subset is absent from Figure 6 data.")
  }
  source_object$mean_ploidy[
    initial_index, day_index, oxygen_index, p_index, drop = FALSE
  ]
}

f6ng_panel_qc <- function(object, source_object, source_config) {
  source_ploidy <- f6ng_source_subset(
    source_object, object$initial_ploidy, object$day_values,
    object$o2_values, object$p_misseg
  )
  replay_error <- max(
    abs(object$mean_ploidy_replay - source_ploidy), na.rm = TRUE
  )
  source_mask_match <- identical(
    is.na(object$mean_ploidy_replay), is.na(source_ploidy)
  )
  deterministic <- identical(object$model_context, "in vivo")
  replay_tolerance <- if (deterministic) 1e-10 else 0.05
  mcse_target <- if (deterministic) 0 else
    as.numeric(source_config$mcse_target_N %||% 0.01)
  maximum_mcse <- max(object$stochastic_mcse, na.rm = TRUE)
  data.frame(
    model_context = object$model_context,
    propagation_mode = object$propagation_mode,
    pair_label = object$pair_label,
    n_initial = length(object$initial_ploidy),
    n_day = length(object$day_values), n_o2 = length(object$o2_values),
    n_p = length(object$p_misseg),
    optimizer_endpoint_weight = unique(object$optimizer_endpoint_weight)[[1L]],
    minimum_net_growth_per_day = min(object$mean_net_growth_rate),
    maximum_net_growth_per_day = max(object$mean_net_growth_rate),
    maximum_growth_mcse_per_day = maximum_mcse,
    growth_mcse_target_per_day = mcse_target,
    maximum_ploidy_replay_error = replay_error,
    ploidy_replay_tolerance = replay_tolerance,
    source_missing_mask_match = source_mask_match,
    validation_policy = if (deterministic) {
      "exact deterministic replay"
    } else {
      "fixed-200 stochastic estimate versus prior variable-repeat estimate"
    },
    passed = replay_error <= replay_tolerance && source_mask_match &&
      (deterministic || maximum_mcse <= mcse_target) &&
      all(is.finite(object$mean_net_growth_rate)),
    stringsAsFactors = FALSE
  )
}

f6ng_aggregate <- function(tasks, source, run_paths, fingerprint, smoke = FALSE) {
  task_qc <- do.call(rbind, lapply(tasks$cache_path, function(path) readRDS(path)$qc))
  task_qc_path <- f6ft_atomic_write_tsv(
    task_qc, file.path(run_paths$run_root, "net_growth_task_validation.tsv")
  )
  if (!all(task_qc$passed)) stop("One or more net-growth tasks failed validation.")
  days <- f6ng_days(smoke)
  p_values <- f6ng_p_values(smoke)
  panel_qc <- list()
  panel_paths <- character()
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
    rate <- ploidy <- mcse <- between <- array(
      NA_real_, dim = shape, dimnames = dimension_names
    )
    selected <- tasks[
      tasks$model_context == context & tasks$pair_label == family, , drop = FALSE
    ]
    for (task_index in seq_len(nrow(selected))) {
      task <- selected[task_index, , drop = FALSE]
      cached <- readRDS(task$cache_path[[1L]])
      if (!identical(cached$fingerprint, fingerprint)) {
        stop("Net-growth aggregate encountered an incompatible task cache.")
      }
      pi <- match(as.numeric(task$p_misseg[[1L]]), p_values)
      rate[, , cached$oxygen_index, pi] <- cached$mean_net_growth_rate
      ploidy[, , cached$oxygen_index, pi] <- cached$mean_ploidy_replay
      mcse[, , cached$oxygen_index, pi] <- cached$stochastic_mcse
      between[, , cached$oxygen_index, pi] <- cached$between_endpoint_variance
    }
    if (any(!is.finite(rate)) || any(!is.finite(ploidy))) {
      stop("Incomplete net-growth aggregate for ", context, " ", family)
    }
    if (isTRUE(smoke)) {
      # Smoke runs replay one unique endpoint, whereas the published source
      # panel contains all 50 optimizer endpoints. Exact kernel identity is
      # tested independently in f6ng_kernel_validation().
      replay_error <- 0
      source_mask_match <- TRUE
    } else {
      source_object <- f6g_read_panel(
        source$run_paths, context, f6ng_mode(context), family
      )
      source_ploidy <- f6ng_source_subset(
        source_object, f6ng_initial_ploidy(), days, oxygen, p_values
      )
      replay_error <- max(abs(ploidy - source_ploidy), na.rm = TRUE)
      source_mask_match <- identical(is.na(ploidy), is.na(source_ploidy))
    }
    represented <- if (isTRUE(smoke)) {
      unique(selected$represented_optimizer_endpoint)[[1L]]
    } else 50L
    object <- list(
      profile = f6ng_profile(), fingerprint = fingerprint,
      source_full_range_run_id = source$run_paths$run_id,
      model_context = context, propagation_mode = f6ng_mode(context),
      pair_label = family, initial_ploidy = f6ng_initial_ploidy(),
      day_values = days, o2_values = oxygen, p_misseg = p_values,
      optimizer_endpoint_weight = rep(represented, length(p_values)),
      mean_net_growth_rate = rate, mean_ploidy_replay = ploidy,
      stochastic_mcse = mcse, between_endpoint_variance = between,
      metric_definition = paste0(
        "sum_N f_N(t)*colSums(M)_N = 1^T M x(t)/1^T x(t); day^-1; ",
        "passage-day state is post-sampling and dilution is excluded"
      )
    )
    panel_path <- file.path(
      run_paths$run_root, f6ng_panel_filename(context, family)
    )
    f6ft_atomic_save_rds(object, panel_path, compress = "gzip")
    panel_paths <- c(panel_paths, panel_path)
    panel_index <- panel_index + 1L
    panel_qc[[panel_index]] <- if (isTRUE(smoke)) {
      data.frame(
        model_context = context, propagation_mode = f6ng_mode(context),
        pair_label = family, n_initial = length(f6ng_initial_ploidy()),
        n_day = length(days), n_o2 = length(oxygen), n_p = length(p_values),
        optimizer_endpoint_weight = represented,
        minimum_net_growth_per_day = min(rate),
        maximum_net_growth_per_day = max(rate),
        maximum_growth_mcse_per_day = max(mcse),
        growth_mcse_target_per_day = source$config$mcse_target_N %||% 0.01,
        maximum_ploidy_replay_error = 0,
        ploidy_replay_tolerance = NA_real_, source_missing_mask_match = TRUE,
        validation_policy = "smoke kernel identity",
        passed = all(is.finite(rate)), stringsAsFactors = FALSE
      )
    } else {
      f6ng_panel_qc(object, source_object, source$config)
    }
  }
  panel_qc <- do.call(rbind, panel_qc)
  panel_qc_path <- f6ft_atomic_write_tsv(
    panel_qc, file.path(run_paths$run_root, "net_growth_panel_validation.tsv")
  )
  if (!all(panel_qc$passed)) stop("Net-growth panel validation failed.")
  list(
    task_qc = task_qc_path, panel_qc = panel_qc_path,
    panels = normalizePath(panel_paths, mustWork = TRUE)
  )
}

f6ng_kernel_validation <- function(
    source, objective_bundle, contexts, paths, run_paths
) {
  endpoint <- source$endpoints[
    source$endpoints$model_context == "in vitro" &
      source$endpoints$pair_label == "C01", , drop = FALSE
  ][1L, , drop = FALSE]
  prepared <- f6ft_prepare_endpoint(endpoint, objective_bundle, contexts)
  # Use a source-grid condition known to exercise many real passage events;
  # this makes the stochastic replay test cover the sampling branch.
  p <- 0.005
  oxygen <- 0
  fixed <- fixo2_fixed_matrix(
    globalenv(), prepared$config,
    figure6_force_p_misseg(prepared$run_params, p), O2 = oxygen
  )
  unit <- prepared$config$N_UNIT %||% 22L
  initial <- f6ft_initial_matrix(fixed$ngrid, unit, 2)
  step <- as.matrix(Matrix::expm(fixed$M))
  net_rate <- as.numeric(colSums(fixed$M))
  seed_number <- as.integer(strsplit(
    endpoint$represented_seed_numbers[[1L]], ",", fixed = TRUE
  )[[1L]][[1L]])
  replicates <- 3L
  streams <- f6s_streams(
    source$catalog, "C01", seed_number, oxygen, p, 2, replicates
  )
  state <- initial[, rep(1L, replicates), drop = FALSE]
  arguments <- list(
    step, state, as.numeric(fixed$ngrid) / unit, net_rate, 1000L,
    source$rule$seed_cells[[1L]], source$rule$target_live_cells[[1L]],
    streams, rep(log(source$rule$seed_cells[[1L]]), replicates), 0L
  )
  first <- do.call(f6ng_propagate_stochastic_cpp, arguments)
  second <- do.call(f6ng_propagate_stochastic_cpp, arguments)
  legacy <- f6s_propagate_cpp(
    step, state, as.numeric(fixed$ngrid) / unit, 1000L,
    source$rule$seed_cells[[1L]], source$rule$target_live_cells[[1L]],
    streams, rep(log(source$rule$seed_cells[[1L]]), replicates), 0L, TRUE, 0L
  )
  random_composition <- seq_along(net_rate)
  random_composition <- random_composition / sum(random_composition)
  identity_error <- abs(
    sum(random_composition * net_rate) -
      sum(as.numeric(fixed$M %*% random_composition))
  )
  replay_identical <- identical(first, second)
  ploidy_error <- max(abs(first$mean_ploidy - legacy$mean_ploidy))
  passage_identical <- identical(first$passage_count, legacy$passage_count)
  passage_exercised <- sum(first$passage_count) > 0
  day0_rate_error <- max(abs(
    first$net_growth_rate[, 1L] - net_rate[which(initial[, 1L] == 1)]
  ))
  validation <- data.frame(
    check = c(
      "fixed_rng_replay_identical", "legacy_ploidy_replay",
      "legacy_passage_count_replay", "rate_identity",
      "actual_passage_exercised", "day0_rate_direct"
    ),
    observed = c(
      as.character(replay_identical), sprintf("%.17g", ploidy_error),
      as.character(passage_identical), sprintf("%.17g", identity_error),
      as.character(passage_exercised), sprintf("%.17g", day0_rate_error)
    ),
    expected = c("TRUE", "0", "TRUE", "0", "TRUE", "0"),
    tolerance = c(NA, 1e-12, NA, 1e-12, NA, 1e-12),
    passed = c(
      replay_identical,
      is.finite(ploidy_error) && ploidy_error <= 1e-12,
      passage_identical,
      is.finite(identity_error) && identity_error <= 1e-12,
      passage_exercised,
      is.finite(day0_rate_error) && day0_rate_error <= 1e-12
    ),
    stringsAsFactors = FALSE
  )
  path <- f6ft_atomic_write_tsv(
    validation, file.path(run_paths$run_root, "net_growth_kernel_validation.tsv")
  )
  if (!all(validation$passed)) stop("Net-growth kernel validation failed.")
  path
}

f6ng_publish_current <- function(paths, run_paths, fingerprint) {
  pointer <- data.frame(
    run_id = run_paths$run_id,
    relative_run_path = file.path(
      "net_growth_q10_runs", run_paths$run_id
    ),
    profile = f6ng_profile(), fingerprint = fingerprint,
    published_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  f6ft_atomic_write_tsv(pointer, run_paths$current)
}

f6ng_finalize_existing <- function(
    workspace_root = f6r_find_workspace_root(), run_id
) {
  paths <- f6r_paths(workspace_root)
  source <- f6ng_source_bundle(paths)
  run_paths <- f6ng_paths(paths, run_id = run_id, create = FALSE)
  task_qc_path <- file.path(run_paths$run_root, "net_growth_task_validation.tsv")
  f6r_require_files(task_qc_path, "completed net-growth task validation")
  task_qc <- f6r_read_tsv(task_qc_path)
  if (nrow(task_qc) != 1020L || !all(task_qc$passed)) {
    stop("Existing net-growth run does not contain 1020 validated tasks.")
  }
  rows <- list()
  fingerprints <- character()
  index <- 0L
  for (context in f6ft_context_levels()) for (family in f6ft_family_levels()) {
    object <- f6ng_read_panel(run_paths, context, family)
    source_object <- f6g_read_panel(
      source$run_paths, context, f6ng_mode(context), family
    )
    index <- index + 1L
    rows[[index]] <- f6ng_panel_qc(object, source_object, source$config)
    fingerprints <- c(fingerprints, object$fingerprint)
  }
  if (length(unique(fingerprints)) != 1L) {
    stop("Existing net-growth panels do not share one computation fingerprint.")
  }
  validation <- do.call(rbind, rows)
  validation_path <- f6ft_atomic_write_tsv(
    validation, file.path(run_paths$run_root, "net_growth_panel_validation.tsv")
  )
  if (!all(validation$passed)) stop("Existing net-growth panel validation failed.")
  provenance <- data.frame(
    key = c(
      "profile", "run_id", "computation_fingerprint",
      "source_full_range_run_id", "metric", "unit",
      "optimizer_endpoint_aggregation", "stochastic_replicates",
      "stochastic_master_seed", "invitro_replay_policy", "finalized_at"
    ),
    value = c(
      f6ng_profile(), run_paths$run_id, unique(fingerprints),
      source$run_paths$run_id,
      "population-weighted instantaneous net-live growth rate = 1^T M f(t)",
      "day^-1", "arithmetic mean with q10 multiplicity restored",
      as.character(source$config$replicates),
      as.character(source$config$master_seed),
      "fixed-200 estimate; prior Figure 6 used condition-specific variable repeats",
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ), stringsAsFactors = FALSE
  )
  provenance_path <- f6ft_atomic_write_tsv(
    provenance, file.path(run_paths$run_root, "net_growth_run_provenance.tsv")
  )
  f6ng_publish_current(paths, run_paths, unique(fingerprints))
  invisible(list(
    validation = validation_path, provenance = provenance_path,
    pointer = run_paths$current
  ))
}

f6ng_data <- function(
    workspace_root = f6r_find_workspace_root(), n_core = 1L,
    run_id = f6ft_resolve_run_id(), smoke = FALSE,
    publish_current = !isTRUE(smoke), o2_chunk_size = 1L
) {
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  f6r_require_packages(c("Matrix", "Rcpp", "future", "future.apply"))
  paths <- f6r_paths(workspace_root)
  run_paths <- f6ng_paths(paths, run_id = run_id, create = TRUE)
  f6r_load_response_engine(paths)
  f6g_load_propagator(paths)
  propagator <- f6ng_load_propagator(paths)
  source <- f6ng_source_bundle(paths)
  objective_bundle <- f6r_objective_selection(paths)
  contexts <- lapply(
    unique(source$endpoints$pair_id), f6r_pair_model_context,
    selected = objective_bundle$selected, paths = paths
  )
  names(contexts) <- unique(source$endpoints$pair_id)
  fingerprint <- f6ng_fingerprint(paths, source, propagator, smoke)
  kernel_validation <- f6ng_kernel_validation(
    source, objective_bundle, contexts, paths, run_paths
  )
  tasks <- f6ng_task_manifest(
    source, run_paths, smoke = smoke, o2_chunk_size = o2_chunk_size
  )
  task_list <- split(tasks, seq_len(nrow(tasks)))
  compute_one <- function(task) tryCatch(
    f6ng_compute_task(
      task, source, objective_bundle, contexts, paths, run_paths,
      fingerprint, smoke = smoke
    ), error = function(error) structure(list(
      task_id = task$task_id[[1L]], message = conditionMessage(error)
    ), class = "f6ng_error")
  )
  message(
    "Supplementary Figure 6-14 net growth: ", nrow(tasks),
    " tasks, workers=", min(as.integer(n_core), nrow(tasks)),
    ", run_id=", run_paths$run_id
  )
  results <- f6ft_parallel_lapply(task_list, compute_one, n_core = n_core)
  failed <- vapply(results, inherits, logical(1L), "f6ng_error")
  if (any(failed)) stop(
    "Net-growth task failures: ",
    paste(vapply(results[failed], `[[`, character(1L), "message"), collapse = "; ")
  )
  aggregation <- f6ng_aggregate(
    tasks, source, run_paths, fingerprint, smoke = smoke
  )
  if (!identical(fingerprint, f6ng_fingerprint(paths, source, propagator, smoke))) {
    stop("Net-growth inputs changed while the run was active.")
  }
  provenance <- data.frame(
    key = c(
      "profile", "run_id", "fingerprint", "workspace_root",
      "model_code_root", "model_source_fingerprint",
      "source_full_range_run_id", "source_full_range_fingerprint",
      "metric", "unit", "initial_ploidy", "day_grid",
      "invivo_o2_grid", "invitro_o2_grid", "p_misseg_grid",
      "optimizer_endpoint_aggregation", "stochastic_replicates",
      "stochastic_master_seed", "passage_day_state", "passage_dilution"
    ),
    value = c(
      f6ng_profile(), run_paths$run_id, fingerprint,
      normalizePath(paths$root, mustWork = TRUE),
      normalizePath(paths$oxygen_code, mustWork = TRUE),
      source$model_fingerprint, source$run_paths$run_id,
      source$pointer$fingerprint[[1L]],
      "population-weighted instantaneous net-live growth rate = 1^T M f(t)",
      "day^-1", paste(f6ng_initial_ploidy(), collapse = ","),
      paste(range(f6ng_days(smoke)), collapse = ":"),
      paste0(paste(range(f6ng_o2("in vivo", smoke)), collapse = ":"),
             ":", length(f6ng_o2("in vivo", smoke))),
      paste0(paste(range(f6ng_o2("in vitro", smoke)), collapse = ":"),
             ":", length(f6ng_o2("in vitro", smoke))),
      paste(f6ng_p_values(smoke), collapse = ","),
      "arithmetic mean with q10 optimizer-endpoint multiplicity restored",
      as.character(source$config$replicates),
      as.character(source$config$master_seed),
      "post-sampling right-continuous composition on each passage day",
      "excluded from biological net growth"
    ), stringsAsFactors = FALSE
  )
  provenance_path <- f6ft_atomic_write_tsv(
    provenance, file.path(run_paths$run_root, "net_growth_run_provenance.tsv")
  )
  if (isTRUE(publish_current)) {
    f6ng_publish_current(paths, run_paths, fingerprint)
  }
  invisible(list(
    paths = run_paths, source = source$run_paths,
    tasks = tasks, kernel_validation = kernel_validation,
    aggregation = aggregation, provenance = provenance_path
  ))
}
