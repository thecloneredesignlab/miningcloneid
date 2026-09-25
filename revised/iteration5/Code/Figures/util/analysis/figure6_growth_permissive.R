#!/usr/bin/env Rscript

# Supplementary Figures 6-20/21: condition on the Figure 6B *instantaneous*
# population-weighted net-live growth rate at day 1000. In vitro eligibility
# is per stochastic trajectory (optimizer seed x passage replicate), not per
# dominant eigenvalue or per already-averaged heatmap cell.

f6gp_profile <- function() "figure6_day1000_positive_growth_q10_v6"
f6gp_p_values <- function() f6ft_p_values()
f6gp_o2_values <- function() seq(0, 5, by = 0.1)
f6gp_initial_values <- function() c(2, 4)
f6gp_day <- function(smoke = FALSE) if (isTRUE(smoke)) 20L else 1000L

f6gp_paths <- function(paths, smoke = FALSE, create = FALSE) {
  root <- file.path(paths$root, "data", "Figures", if (isTRUE(smoke)) {
    "Supp_Figure6_20_21_day1000_positive_growth_v6_smoke"
  } else "Supp_Figure6_20_21_day1000_positive_growth_v6")
  result <- list(
    root = root,
    cache = file.path(root, "day1000_positive_growth_task_cache"),
    task = file.path(root, "day1000_positive_growth_tasks.tsv"),
    curve = file.path(root, "day1000_positive_growth_curves.tsv"),
    provenance = file.path(root, "day1000_positive_growth_provenance.tsv"),
    source_validation = file.path(root, "day1000_positive_growth_validation.tsv"),
    monotonicity = file.path(root, "day1000_positive_growth_monotonicity.tsv")
  )
  if (isTRUE(create)) {
    dir.create(result$cache, recursive = TRUE, showWarnings = FALSE)
  }
  result
}

f6gp_fingerprint <- function(paths, source, smoke = FALSE) {
  inputs <- c(
    unname(source$files),
    file.path(paths$code, "util", "analysis", "figure6_growth_permissive.R"),
    file.path(paths$code, "util", "analysis", "figure6_net_growth_q10.R"),
    file.path(paths$code, "util", "analysis", "figure6_net_growth_propagator.cpp")
  )
  f6r_require_files(inputs, "day-1000 positive-growth fingerprint inputs")
  paste(
    f6gp_profile(), source$run_paths$run_id, source$model_fingerprint,
    paste(unname(tools::md5sum(inputs)), collapse = ":"),
    paste(f6gp_p_values(), collapse = ","),
    paste(f6gp_o2_values(), collapse = ","),
    f6gp_day(smoke), source$config$replicates,
    if (isTRUE(smoke)) "smoke" else "full", sep = "|"
  )
}

f6gp_prepare <- function(workspace_root, smoke = FALSE) {
  paths <- f6r_paths(workspace_root)
  output <- f6gp_paths(paths, smoke = smoke, create = TRUE)
  source <- f6ng_source_bundle(paths)
  fingerprint <- f6gp_fingerprint(paths, source, smoke)
  oxygen <- if (isTRUE(smoke)) 0.5 else f6gp_o2_values()
  p_values <- if (isTRUE(smoke)) 0.1 else f6gp_p_values()
  grid <- expand.grid(
    model_context = f6ft_context_levels(),
    pair_label = f6ft_family_levels(),
    p_misseg = p_values, O2_pct = oxygen,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid$task_index <- seq_len(nrow(grid))
  grid$cache_path <- file.path(
    output$cache, sprintf("day1000_positive_growth_%04d.rds", grid$task_index)
  )
  if (file.exists(output$task)) {
    old <- f6r_read_tsv(output$task)
    if (!isTRUE(all.equal(old[, names(grid), drop = FALSE], grid,
                          tolerance = 1e-12, check.attributes = FALSE))) {
      stop("Existing positive-growth task manifest differs; use a fresh output root.")
    }
  } else f6ft_atomic_write_tsv(grid, output$task)
  provenance <- data.frame(
    key = c("profile", "source_full_range_run_id", "source_fingerprint",
            "model_code_root", "model_source_fingerprint", "fingerprint",
            "condition", "growth_metric", "growth_unit", "invitro_rule",
            "invitro_replicates", "invitro_master_seed", "initial_ploidy",
            "day", "oxygen_grid", "p_misseg_grid"),
    value = c(
      f6gp_profile(), source$run_paths$run_id, source$pointer$fingerprint[[1L]],
      normalizePath(paths$oxygen_code, mustWork = TRUE),
      source$model_fingerprint, fingerprint,
      "per simulated trajectory at the selected day: growth > 0",
      "population-weighted instantaneous net-live growth = 1^T M f(t)",
      "day^-1", "threshold-triggered stochastic passage; post-sampling state",
      source$config$replicates, source$config$master_seed,
      paste(f6gp_initial_values(), collapse = ","), f6gp_day(smoke),
      paste(range(oxygen), collapse = ":"), paste(p_values, collapse = ",")
    ), stringsAsFactors = FALSE
  )
  if (file.exists(output$provenance)) {
    old <- f6r_read_tsv(output$provenance)
    if (!identical(old$key, provenance$key) ||
        !identical(as.character(old$value), as.character(provenance$value))) {
      stop("Existing positive-growth provenance differs; refuse cache reuse.")
    }
  } else f6ft_atomic_write_tsv(provenance, output$provenance)
  list(paths = paths, output = output, source = source,
       fingerprint = fingerprint, tasks = grid)
}

f6gp_compute_task <- function(task, source, objective_bundle, contexts,
                              paths, fingerprint, smoke = FALSE) {
  output <- as.character(task$cache_path[[1L]])
  if (file.exists(output)) {
    cached <- readRDS(output)
    if (!identical(cached$fingerprint, fingerprint) ||
        !identical(cached$task$task_index, task$task_index) ||
        !isTRUE(cached$qc$passed)) {
      stop("Incompatible positive-growth task cache: ", output)
    }
    return(cached$qc)
  }
  context <- as.character(task$model_context[[1L]])
  family <- as.character(task$pair_label[[1L]])
  oxygen <- as.numeric(task$O2_pct[[1L]])
  p <- as.numeric(task$p_misseg[[1L]])
  day <- f6gp_day(smoke)
  endpoints <- source$endpoints[
    source$endpoints$model_context == context &
      source$endpoints$pair_label == family, , drop = FALSE
  ]
  endpoints <- endpoints[order(endpoints$representative_objective_rank,
                               endpoints$representative_seed_number), , drop = FALSE]
  if (isTRUE(smoke)) endpoints <- endpoints[1L, , drop = FALSE]
  if (!isTRUE(smoke) && sum(endpoints$endpoint_multiplicity_q10) != 50L) {
    stop("The q10 group does not represent exactly 50 optimizer endpoints.")
  }
  if (!exists("f6ng_propagate_continuous_cpp", mode = "function",
              inherits = TRUE)) f6ng_load_propagator(paths)
  initial_values <- f6gp_initial_values()
  n_all <- n_positive <- sum_all_ploidy <- sum_positive_ploidy <-
    sum_all_growth <- numeric(length(initial_values))
  passage_count <- 0L
  formula_error <- 0
  for (index in seq_len(nrow(endpoints))) {
    endpoint <- endpoints[index, , drop = FALSE]
    prepared <- f6ft_prepare_endpoint(endpoint, objective_bundle, contexts)
    forced <- figure6_force_p_misseg(prepared$run_params, p)
    formula <- figure6_p_misseg_formula_qc(
      prepared$run_params, p, o2_values = oxygen
    )
    formula_error <- max(formula_error, formula$maximum_direct_formula_error)
    fixed <- fixo2_fixed_matrix(
      globalenv(), prepared$config, forced, O2 = oxygen
    )
    unit <- prepared$config$N_UNIT %||% 22L
    initial <- f6ft_initial_matrix(fixed$ngrid, unit, initial_values)
    step <- as.matrix(Matrix::expm(fixed$M))
    net_live_rate <- as.numeric(colSums(fixed$M))
    ploidy_grid <- as.numeric(fixed$ngrid) / unit
    weight <- as.integer(endpoint$endpoint_multiplicity_q10[[1L]])
    if (identical(context, "in vivo")) {
      response <- f6ng_propagate_continuous_cpp(
        step, initial, ploidy_grid, net_live_rate, day
      )
      for (j in seq_along(initial_values)) {
        mean <- as.numeric(response$mean_ploidy[j, day + 1L])
        growth <- as.numeric(response$net_growth_rate[j, day + 1L])
        if (!is.finite(mean) || !is.finite(growth)) {
          stop("Non-finite in-vivo day-1000 trajectory.")
        }
        n_all[j] <- n_all[j] + weight
        sum_all_ploidy[j] <- sum_all_ploidy[j] + weight * mean
        sum_all_growth[j] <- sum_all_growth[j] + weight * growth
        if (growth > 0) {
          n_positive[j] <- n_positive[j] + weight
          sum_positive_ploidy[j] <- sum_positive_ploidy[j] + weight * mean
        }
      }
    } else {
      seeds <- as.integer(strsplit(endpoint$represented_seed_numbers[[1L]],
                                   ",", fixed = TRUE)[[1L]])
      replicates <- as.integer(source$config$replicates)
      if (length(seeds) != weight) {
        stop("Stochastic seed count differs from endpoint multiplicity.")
      }
      for (j in seq_along(initial_values)) {
        streams <- f6s_streams(
          source$catalog, family, seeds, oxygen, p,
          initial_values[[j]], replicates
        )
        state <- initial[, rep(j, weight * replicates), drop = FALSE]
        response <- f6ng_propagate_stochastic_cpp(
          step, state, ploidy_grid, net_live_rate, day,
          source$rule$seed_cells[[1L]],
          source$rule$target_live_cells[[1L]], streams,
          rep(log(source$rule$seed_cells[[1L]]), weight * replicates), 0L
        )
        mean <- as.numeric(response$mean_ploidy[, day + 1L])
        growth <- as.numeric(response$net_growth_rate[, day + 1L])
        if (length(mean) != weight * replicates ||
            any(!is.finite(mean)) || any(!is.finite(growth))) {
          stop("Invalid in-vitro day-1000 stochastic trajectories.")
        }
        positive <- growth > 0
        n_all[j] <- n_all[j] + length(mean)
        n_positive[j] <- n_positive[j] + sum(positive)
        sum_all_ploidy[j] <- sum_all_ploidy[j] + sum(mean)
        sum_positive_ploidy[j] <- sum_positive_ploidy[j] + sum(mean[positive])
        sum_all_growth[j] <- sum_all_growth[j] + sum(growth)
        passage_count <- passage_count + sum(response$passage_count)
      }
    }
  }
  expected <- sum(endpoints$endpoint_multiplicity_q10) * if (
    identical(context, "in vitro")) as.integer(source$config$replicates) else 1L
  if (any(n_all != expected) || any(n_positive < 0) ||
      any(n_positive > n_all) || formula_error > 1e-12) {
    stop("Positive-growth task denominator or model-formula validation failed.")
  }
  rows <- data.frame(
    model_context = context, pair_label = family,
    initial_ploidy = initial_values, p_misseg = p, O2_pct = oxygen,
    day = day, unfiltered_mean_ploidy = sum_all_ploidy / n_all,
    growth_positive_mean_ploidy = ifelse(
      n_positive > 0, sum_positive_ploidy / pmax(n_positive, 1), NA_real_
    ),
    unfiltered_mean_growth_per_day = sum_all_growth / n_all,
    growth_positive_count = as.integer(n_positive),
    trajectory_count = as.integer(n_all),
    growth_positive_fraction = n_positive / n_all,
    stringsAsFactors = FALSE
  )
  qc <- data.frame(
    task_index = as.integer(task$task_index[[1L]]),
    model_context = context, pair_label = family, p_misseg = p,
    O2_pct = oxygen, n_trajectory_per_initial = expected,
    n_positive_2N = as.integer(n_positive[[1L]]),
    n_positive_4N = as.integer(n_positive[[2L]]),
    n_passages = as.integer(passage_count),
    maximum_formula_error = formula_error,
    passed = TRUE, stringsAsFactors = FALSE
  )
  f6ft_atomic_save_rds(list(
    profile = f6gp_profile(), fingerprint = fingerprint,
    task = task, rows = rows, qc = qc
  ), output, compress = "gzip")
  qc
}

f6gp_monotonicity <- function(curves) {
  groups <- split(curves, interaction(
    curves$model_context, curves$pair_label, curves$initial_ploidy,
    curves$p_misseg, drop = TRUE
  ))
  result <- lapply(groups, function(x) {
    x <- x[order(x$O2_pct), , drop = FALSE]
    kept <- x[is.finite(x$growth_positive_mean_ploidy), , drop = FALSE]
    adjacent <- which(abs(diff(kept$O2_pct) - 0.1) < 1e-10)
    delta <- if (length(adjacent)) diff(kept$growth_positive_mean_ploidy)[
      adjacent] else numeric()
    data.frame(
      model_context = x$model_context[[1L]],
      pair_label = x$pair_label[[1L]],
      initial_ploidy = x$initial_ploidy[[1L]],
      p_misseg = x$p_misseg[[1L]],
      n_positive_o2 = nrow(kept), n_adjacent_pairs = length(delta),
      n_adjacent_increase = sum(delta > 1e-8),
      n_adjacent_decrease = sum(delta < -1e-8),
      maximum_adjacent_increase_N = if (length(delta)) max(delta) else NA_real_,
      minimum_adjacent_change_N = if (length(delta)) min(delta) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

f6gp_finalize <- function(workspace_root, smoke = FALSE) {
  setup <- f6gp_prepare(workspace_root, smoke = smoke)
  rows <- qc <- vector("list", nrow(setup$tasks))
  for (i in seq_len(nrow(setup$tasks))) {
    file <- setup$tasks$cache_path[[i]]
    f6r_require_files(file, "day-1000 positive-growth task cache")
    object <- readRDS(file)
    if (!identical(object$profile, f6gp_profile()) ||
        !identical(object$fingerprint, setup$fingerprint) ||
        !identical(as.integer(object$task$task_index[[1L]]), i) ||
        !isTRUE(object$qc$passed)) {
      stop("Incompatible positive-growth task at index ", i)
    }
    rows[[i]] <- object$rows
    qc[[i]] <- object$qc
  }
  curves <- do.call(rbind, rows)
  rownames(curves) <- NULL
  if (nrow(curves) != 2L * nrow(setup$tasks) ||
      any(!is.finite(curves$unfiltered_mean_ploidy)) ||
      any(!is.finite(curves$unfiltered_mean_growth_per_day)) ||
      any(curves$growth_positive_count > curves$trajectory_count) ||
      any(is.na(curves$growth_positive_mean_ploidy) !=
            (curves$growth_positive_count == 0L))) {
    stop("Incomplete positive-growth summary.")
  }
  validation <- list()
  if (!isTRUE(smoke)) {
    growth_run <- f6ng_paths(setup$paths)
    growth_manifest <- f6r_read_tsv(file.path(
      growth_run$run_root, "net_growth_task_manifest.tsv"
    ))
    for (context in f6ft_context_levels()) for (family in f6ft_family_levels()) {
      mode <- f6ng_mode(context)
      original <- f6g_read_panel(setup$source$run_paths, context, mode, family)
      growth <- f6ng_read_panel(growth_run, context, family)
      selected <- curves[curves$model_context == context &
                           curves$pair_label == family, , drop = FALSE]
      delta_source_ploidy <- delta_replay_ploidy <- delta_growth <-
        numeric(nrow(selected))
      last_cache_file <- ""
      cached_replay <- NULL
      for (j in seq_len(nrow(selected))) {
        row <- selected[j, , drop = FALSE]
        indices <- c(
          match(row$initial_ploidy, original$initial_ploidy),
          match(row$day, original$day_values),
          match(sprintf("%.3f", row$O2_pct),
                sprintf("%.3f", original$o2_values)),
          match(sprintf("%.3f", row$p_misseg),
                sprintf("%.3f", original$p_misseg))
        )
        growth_indices <- c(
          match(row$initial_ploidy, growth$initial_ploidy),
          match(row$day, growth$day_values),
          match(sprintf("%.3f", row$O2_pct),
                sprintf("%.3f", growth$o2_values)),
          match(sprintf("%.3f", row$p_misseg),
                sprintf("%.3f", growth$p_misseg))
        )
        if (anyNA(c(indices, growth_indices))) {
          stop("Original Figure 6B grid does not cover replayed cell.")
        }
        delta_source_ploidy[[j]] <- abs(
          row$unfiltered_mean_ploidy -
            original$mean_ploidy[indices[1], indices[2], indices[3], indices[4]]
        )
        replay_task <- growth_manifest[
          growth_manifest$model_context == context &
            growth_manifest$pair_label == family &
            abs(growth_manifest$p_misseg - row$p_misseg) < 1e-12 &
            growth_manifest$o2_index_start <= growth_indices[3] &
            growth_manifest$o2_index_end >= growth_indices[3],
          , drop = FALSE
        ]
        if (nrow(replay_task) != 1L) {
          stop("Cannot identify unique fixed-200 net-growth replay task.")
        }
        replay_file <- file.path(
          growth_run$cache, basename(replay_task$cache_path[[1L]])
        )
        if (!identical(replay_file, last_cache_file)) {
          f6r_require_files(replay_file, "fixed-200 net-growth replay cache")
          cached_replay <- readRDS(replay_file)
          if (!identical(cached_replay$profile, f6ng_profile()) ||
              !isTRUE(cached_replay$qc$passed)) {
            stop("Invalid fixed-200 net-growth replay cache.")
          }
          last_cache_file <- replay_file
        }
        replay_indices <- c(
          match(row$initial_ploidy, cached_replay$initial_ploidy),
          match(row$day, cached_replay$day_values),
          match(sprintf("%.3f", row$O2_pct),
                sprintf("%.3f", cached_replay$o2_values))
        )
        if (anyNA(replay_indices)) {
          stop("Fixed-200 replay cache does not cover requested cell.")
        }
        delta_replay_ploidy[[j]] <- abs(
          row$unfiltered_mean_ploidy -
            cached_replay$mean_ploidy_replay[
              replay_indices[1], replay_indices[2], replay_indices[3]
            ]
        )
        delta_growth[[j]] <- abs(
          row$unfiltered_mean_growth_per_day -
            growth$mean_net_growth_rate[
              growth_indices[1], growth_indices[2],
              growth_indices[3], growth_indices[4]
            ]
        )
      }
      validation[[length(validation) + 1L]] <- data.frame(
        model_context = context, pair_label = family,
        source_figure6_run_id = setup$source$run_paths$run_id,
        source_net_growth_run_id = growth_run$run_id,
        maximum_ploidy_replay_error = max(delta_replay_ploidy),
        maximum_original_figure6b_ploidy_difference =
          max(delta_source_ploidy),
        original_figure6b_ploidy_tolerance = if (
          identical(context, "in vivo")) 1e-8 else 0.05,
        maximum_growth_replay_error_per_day = max(delta_growth),
        passed = max(delta_replay_ploidy) <= 1e-8 &&
          max(delta_growth) <= 1e-8 &&
          max(delta_source_ploidy) <= if (
            identical(context, "in vivo")) 1e-8 else 0.05,
        stringsAsFactors = FALSE
      )
      rm(original, growth)
      invisible(gc(FALSE))
    }
  }
  validation <- if (length(validation)) do.call(rbind, validation) else {
    data.frame(model_context = "smoke", pair_label = "smoke",
               passed = TRUE)
  }
  if (!all(validation$passed)) {
    stop("All-data replay differs from Figure 6B or its growth companion.")
  }
  f6ft_atomic_write_tsv(curves, setup$output$curve)
  f6ft_atomic_write_tsv(f6gp_monotonicity(curves),
                        setup$output$monotonicity)
  f6ft_atomic_write_tsv(validation, setup$output$source_validation)
  invisible(list(curves = curves, validation = validation,
                 output = setup$output))
}

f6gp_run <- function(workspace_root, n_core = 1L, smoke = FALSE) {
  setup <- f6gp_prepare(workspace_root, smoke = smoke)
  f6r_require_packages(c("Matrix", "Rcpp", "future", "future.apply"))
  f6r_load_response_engine(setup$paths)
  f6ng_load_propagator(setup$paths)
  objective_bundle <- f6r_objective_selection(setup$paths)
  contexts <- lapply(
    unique(setup$source$endpoints$pair_id), f6r_pair_model_context,
    selected = objective_bundle$selected, paths = setup$paths
  )
  names(contexts) <- unique(setup$source$endpoints$pair_id)
  tasks <- split(setup$tasks, seq_len(nrow(setup$tasks)))
  message("Day-1000 positive-growth tasks=", length(tasks),
          "; workers=", min(as.integer(n_core), length(tasks)))
  result <- f6ft_parallel_lapply(tasks, function(task) tryCatch(
    f6gp_compute_task(task, setup$source, objective_bundle, contexts,
                      setup$paths, setup$fingerprint, smoke = smoke),
    error = function(e) list(error = conditionMessage(e),
                             task_index = task$task_index[[1L]])
  ), n_core = n_core)
  failed <- vapply(result, function(x) !is.null(x$error), logical(1L))
  if (any(failed)) {
    detail <- vapply(result[failed], function(x) paste0(
      x$task_index, ":", x$error
    ), character(1L))
    stop("Positive-growth task failures: ", paste(detail, collapse = "; "))
  }
  f6gp_finalize(workspace_root, smoke = smoke)
}
