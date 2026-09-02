#!/usr/bin/env Rscript

# Finite-time q10 endpoint analysis for Figure 7C-F and Supplementary
# Figures 7-5 through 7-7.
#
# Scientific boundary: this file contains only figure-workspace orchestration.
# Model construction and fixed-O2 numerical helpers are loaded from the
# read-only O2_supply_demand_MAP root resolved by f7r_paths().

options(stringsAsFactors = FALSE, warn = 1)

f7ft_profile <- function() "figure7_fixed_pmisseg_finite_time_q10_v1"
f7ft_family_levels <- function() c("C01", "C02")
f7ft_context_levels <- function() c("in vivo", "in vitro")
f7ft_initial_ploidy <- function() 2:6
f7ft_p_values <- function() c(0.005, 0.01, 0.10, 0.20, 0.30)
f7ft_o2_values <- function(smoke = FALSE) {
  if (isTRUE(smoke)) c(0, 1, 5) else seq(0, 5, length.out = 201L)
}
f7ft_day_values <- function(smoke = FALSE) {
  if (isTRUE(smoke)) 0:10 else 0:1000
}
f7ft_diagnostic_o2 <- function() c(0, 0.5, 1, 2, 5)
f7ft_diagnostic_days <- function() c(1, 10, 25, 50, 100, 250, 500, 1000)

f7ft_sanitize_run_id <- function(run_id) {
  run_id <- trimws(as.character(run_id[[1L]]))
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$", run_id)) {
    stop("Invalid Figure 7 finite-time run id: ", run_id, call. = FALSE)
  }
  run_id
}

f7ft_resolve_run_id <- function(run_id = Sys.getenv(
    "FIGURE7_FINITE_TIME_RUN_ID", unset = ""
)) {
  if (!nzchar(trimws(run_id))) {
    run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  }
  f7ft_sanitize_run_id(run_id)
}

f7ft_paths <- function(paths, run_id = NULL, create = FALSE) {
  base <- file.path(paths$figure7, "finite_time_q10_runs")
  current <- file.path(paths$figure7, "finite_time_q10_current.tsv")
  if (is.null(run_id)) {
    if (!file.exists(current)) {
      stop("Missing current Figure 7 finite-time run pointer: ", current)
    }
    pointer <- f7r_read_tsv(current)
    if (nrow(pointer) != 1L || !all(c("run_id", "relative_run_path") %in% names(pointer))) {
      stop("Malformed Figure 7 finite-time run pointer: ", current)
    }
    run_id <- f7ft_sanitize_run_id(pointer$run_id[[1L]])
    run_root <- normalizePath(
      file.path(paths$figure7, pointer$relative_run_path[[1L]]),
      mustWork = TRUE
    )
  } else {
    run_id <- f7ft_sanitize_run_id(run_id)
    run_root <- file.path(base, run_id)
  }
  out <- list(
    run_id = run_id,
    base = base,
    run_root = run_root,
    cache = file.path(run_root, "task_cache"),
    diagnostics_cache = file.path(run_root, "diagnostic_cache"),
    panels = file.path(run_root, "panels"),
    rendered = file.path(run_root, "rendered"),
    current = current,
    supp7_5 = file.path(paths$root, "data", "Figures", "Supp_Figure7_5"),
    supp7_6 = file.path(paths$root, "data", "Figures", "Supp_Figure7_6"),
    supp7_7 = file.path(paths$root, "data", "Figures", "Supp_Figure7_7")
  )
  if (isTRUE(create)) {
    dirs <- unname(unlist(out[c(
      "base", "run_root", "cache", "diagnostics_cache", "panels", "rendered",
      "supp7_5", "supp7_6", "supp7_7"
    )]))
    invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  }
  out
}

f7ft_atomic_save_rds <- function(object, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = compress)
  if (!file.rename(temporary, path)) {
    if (!file.copy(temporary, path, overwrite = TRUE, copy.mode = TRUE)) {
      stop("Failed to publish RDS cache: ", path)
    }
    unlink(temporary)
  }
  normalizePath(path, mustWork = TRUE)
}

f7ft_atomic_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.table(
    x, file = temporary, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
  if (!file.rename(temporary, path)) {
    if (!file.copy(temporary, path, overwrite = TRUE, copy.mode = TRUE)) {
      stop("Failed to publish TSV: ", path)
    }
    unlink(temporary)
  }
  normalizePath(path, mustWork = TRUE)
}

f7ft_context_spec <- function() {
  data.frame(
    model_context = f7ft_context_levels(),
    parameter_value_column = c("vivo_natural", "vitro_natural"),
    signature_column = c("parameter_signature", "parameter_signature_invitro"),
    endpoint_group_column = c(
      "parameter_endpoint_group", "parameter_endpoint_group_invitro"
    ),
    simulation_mode = c("joint", "invitro"),
    panel_letter = c("C", "E"),
    stringsAsFactors = FALSE
  )
}

f7ft_build_endpoint_manifest <- function(paths, objective_bundle, run_paths) {
  acceptance <- objective_bundle$objectives
  required <- c(
    "pair_id", "pair_label", "seed_number", "seed", "objective",
    "objective_rank", "hard_qc_pass", "eligible_q10",
    "parameter_signature", "parameter_endpoint_group",
    "parameter_signature_invitro", "parameter_endpoint_group_invitro"
  )
  if (!all(required %in% names(acceptance))) {
    stop("Figure 7 q10 acceptance table lacks required endpoint fields.")
  }
  selected <- acceptance[
    acceptance$hard_qc_pass & acceptance$eligible_q10 &
      acceptance$pair_label %in% f7ft_family_levels(),
    , drop = FALSE
  ]
  if (!identical(sort(unique(selected$pair_label)), f7ft_family_levels()) ||
      any(table(selected$pair_label) != 50L)) {
    stop("Figure 7 finite-time analysis requires exactly 50 q10 seeds for C01 and C02.")
  }

  context_spec <- f7ft_context_spec()
  endpoint_rows <- list()
  expanded_rows <- list()
  endpoint_index <- 0L
  expanded_index <- 0L
  for (context_index in seq_len(nrow(context_spec))) {
    spec <- context_spec[context_index, , drop = FALSE]
    signature_column <- spec$signature_column[[1L]]
    group_column <- spec$endpoint_group_column[[1L]]
    for (family in f7ft_family_levels()) {
      z <- selected[selected$pair_label == family, , drop = FALSE]
      z <- z[order(z$objective_rank, z$seed_number), , drop = FALSE]
      groups <- split(z, z[[group_column]], drop = TRUE)
      groups <- groups[order(vapply(
        groups, function(g) min(g$objective_rank), numeric(1L)
      ))]
      for (group in groups) {
        group <- group[order(group$objective_rank, group$seed_number), , drop = FALSE]
        endpoint_index <- endpoint_index + 1L
        endpoint_rows[[endpoint_index]] <- data.frame(
          model_context = spec$model_context[[1L]],
          pair_id = group$pair_id[[1L]],
          pair_label = family,
          panel_letter = if (spec$model_context[[1L]] == "in vivo") {
            if (family == "C01") "C" else "D"
          } else {
            if (family == "C01") "E" else "F"
          },
          endpoint_group = group[[group_column]][[1L]],
          parameter_signature = group[[signature_column]][[1L]],
          representative_seed_number = group$seed_number[[1L]],
          representative_seed = group$seed[[1L]],
          representative_objective = group$objective[[1L]],
          representative_objective_rank = group$objective_rank[[1L]],
          endpoint_multiplicity_q10 = nrow(group),
          represented_seed_numbers = paste(group$seed_number, collapse = ","),
          parameter_value_column = spec$parameter_value_column[[1L]],
          simulation_mode = spec$simulation_mode[[1L]],
          stringsAsFactors = FALSE
        )
        for (j in seq_len(nrow(group))) {
          expanded_index <- expanded_index + 1L
          expanded_rows[[expanded_index]] <- data.frame(
            model_context = spec$model_context[[1L]],
            pair_id = group$pair_id[[j]],
            pair_label = family,
            endpoint_group = group[[group_column]][[j]],
            seed = group$seed[[j]],
            seed_number = group$seed_number[[j]],
            objective = group$objective[[j]],
            objective_rank = group$objective_rank[[j]],
            parameter_signature = group[[signature_column]][[j]],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  endpoints <- do.call(rbind, endpoint_rows)
  expanded <- do.call(rbind, expanded_rows)
  rownames(endpoints) <- NULL
  rownames(expanded) <- NULL
  endpoints$endpoint_index <- seq_len(nrow(endpoints))
  endpoints <- endpoints[, c("endpoint_index", setdiff(names(endpoints), "endpoint_index"))]

  represented <- stats::aggregate(
    endpoint_multiplicity_q10 ~ model_context + pair_label,
    endpoints, sum
  )
  unique_counts <- stats::aggregate(
    endpoint_group ~ model_context + pair_label,
    endpoints, length
  )
  names(unique_counts)[names(unique_counts) == "endpoint_group"] <-
    "n_unique_parameter_endpoint"
  qc <- merge(represented, unique_counts,
              by = c("model_context", "pair_label"), sort = FALSE)
  names(qc)[names(qc) == "endpoint_multiplicity_q10"] <- "n_optimizer_endpoint"
  qc$expected_optimizer_endpoint <- 50L
  qc$passed <- qc$n_optimizer_endpoint == qc$expected_optimizer_endpoint
  qc <- qc[order(match(qc$model_context, f7ft_context_levels()),
                   match(qc$pair_label, f7ft_family_levels())), ]
  if (nrow(qc) != 4L || !all(qc$passed)) {
    stop("Figure 7 finite-time q10 endpoint manifest validation failed.")
  }

  paths_out <- c(
    endpoints = f7ft_atomic_write_tsv(
      endpoints, file.path(run_paths$run_root, "q10_unique_endpoint_manifest.tsv")
    ),
    expanded = f7ft_atomic_write_tsv(
      expanded, file.path(run_paths$run_root, "q10_optimizer_seed_manifest.tsv")
    ),
    qc = f7ft_atomic_write_tsv(
      qc, file.path(run_paths$run_root, "q10_endpoint_manifest_validation.tsv")
    )
  )
  list(endpoints = endpoints, expanded = expanded, qc = qc, paths = paths_out)
}

f7ft_model_fingerprint <- function(paths, manifest_bundle) {
  model <- f7r_model_source_fingerprint(paths)
  inputs <- c(
    manifest_bundle$paths[["endpoints"]],
    manifest_bundle$paths[["expanded"]]
  )
  paste(
    f7ft_profile(),
    paste0("model=", model),
    paste0("manifest=", paste(unname(tools::md5sum(inputs)), collapse = ":")),
    paste0("initial=", paste(f7ft_initial_ploidy(), collapse = ",")),
    paste0("p=", paste(format(f7ft_p_values(), scientific = FALSE), collapse = ",")),
    paste0("o2=", paste(range(f7ft_o2_values()), collapse = ":"), ":", length(f7ft_o2_values())),
    paste0("day=", paste(range(f7ft_day_values()), collapse = ":"), ":", length(f7ft_day_values())),
    sep = "|"
  )
}

f7ft_prepare_endpoint <- function(endpoint, objective_bundle, contexts) {
  pair_id <- endpoint$pair_id[[1L]]
  seed_number <- as.integer(endpoint$representative_seed_number[[1L]])
  rows <- objective_bundle$master[
    objective_bundle$master$pair_id == pair_id &
      objective_bundle$master$seed_number == seed_number,
    , drop = FALSE
  ]
  shared <- f7r_shared_parameters()
  value_column <- endpoint$parameter_value_column[[1L]]
  if (nrow(rows) != length(shared) || !setequal(rows$parameter, shared) ||
      !value_column %in% names(rows)) {
    stop("Incomplete endpoint parameters for ", endpoint$endpoint_group[[1L]])
  }
  params <- stats::setNames(as.numeric(rows[[value_column]]), rows$parameter)
  context <- contexts[[pair_id]]
  simulation_mode <- endpoint$simulation_mode[[1L]]
  if (identical(simulation_mode, "joint")) {
    params[["rho_2N"]] <- context$rho_2N
  }
  run_params <- prepare_run_params(
    param_values = params, simulation = simulation_mode,
    cfg = context$config, fixed_o2 = 5
  )
  config <- prepare_sim_cfg(
    context$config, argv = list(), fixed_o2 = 5, run_params = run_params
  )
  run_params$O2_growth <- isTRUE(config$O2_growth)
  run_params$ploidy_O2_death <- config$ploidy_O2_death
  list(run_params = run_params, config = config, context = context)
}

f7ft_initial_matrix <- function(ngrid, n_unit, initial_ploidy = f7ft_initial_ploidy()) {
  target_n <- as.numeric(initial_ploidy) * as.numeric(n_unit)
  index <- match(target_n, as.numeric(ngrid))
  if (anyNA(index)) {
    stop("The model state grid does not contain all requested initial ploidies.")
  }
  x <- matrix(0, nrow = length(ngrid), ncol = length(initial_ploidy))
  x[cbind(index, seq_along(index))] <- 1
  colnames(x) <- paste0("init_", initial_ploidy, "N")
  x
}

f7ft_normalize_columns <- function(x) {
  x <- Re(as.matrix(x))
  x[!is.finite(x)] <- NA_real_
  x[x < 0 & x > -1e-10] <- 0
  x <- pmax(x, 0)
  totals <- colSums(x, na.rm = TRUE)
  bad <- !is.finite(totals) | totals <= 0
  if (any(bad)) x[, bad] <- NA_real_
  good <- !bad
  if (any(good)) x[, good] <- sweep(x[, good, drop = FALSE], 2L, totals[good], "/")
  x
}

f7ft_rescale_columns <- function(x) {
  x <- Re(as.matrix(x))
  scales <- apply(abs(x), 2L, max, na.rm = TRUE)
  bad <- !is.finite(scales) | scales <= 0
  if (any(bad)) x[, bad] <- NA_real_
  good <- !bad
  if (any(good)) x[, good] <- sweep(x[, good, drop = FALSE], 2L, scales[good], "/")
  x
}

f7ft_expm_daily_mean_ploidy <- function(
    M, ngrid, n_unit, day_values = f7ft_day_values()
) {
  if (!identical(as.integer(day_values), seq.int(0L, max(day_values)))) {
    stop("Daily Figure 7 propagation requires a complete integer day grid from zero.")
  }
  state <- f7ft_initial_matrix(ngrid, n_unit)
  ploidy_grid <- as.numeric(ngrid) / as.numeric(n_unit)
  out <- matrix(
    NA_real_, nrow = length(f7ft_initial_ploidy()), ncol = length(day_values),
    dimnames = list(paste0(f7ft_initial_ploidy(), "N"), paste0("day", day_values))
  )
  out[, 1L] <- as.numeric(crossprod(ploidy_grid, state))
  step <- as.matrix(Matrix::expm(M))
  for (day_index in 2:length(day_values)) {
    state <- f7ft_normalize_columns(step %*% state)
    out[, day_index] <- as.numeric(crossprod(ploidy_grid, state))
  }
  out
}

f7ft_task_manifest <- function(
    endpoint_manifest, run_paths, smoke = FALSE, chunk_size = 2L
) {
  endpoints <- endpoint_manifest$endpoints
  if (isTRUE(smoke)) {
    endpoints <- do.call(rbind, lapply(
      split(endpoints, interaction(endpoints$model_context, endpoints$pair_label)),
      function(z) z[1L, , drop = FALSE]
    ))
    p_values <- f7ft_p_values()[1:2]
  } else {
    p_values <- f7ft_p_values()
  }
  rows <- list()
  index <- 0L
  strata <- split(
    endpoints,
    interaction(endpoints$model_context, endpoints$pair_label,
                drop = TRUE, lex.order = TRUE)
  )
  for (z in strata) {
    z <- z[order(z$representative_objective_rank,
                 z$representative_seed_number), , drop = FALSE]
    chunks <- split(seq_len(nrow(z)), ceiling(seq_len(nrow(z)) / chunk_size))
    for (p_value in p_values) {
      for (chunk_index in seq_along(chunks)) {
        index <- index + 1L
        selected <- z[chunks[[chunk_index]], , drop = FALSE]
        rows[[index]] <- data.frame(
          task_id = sprintf("T%04d", index),
          model_context = z$model_context[[1L]],
          pair_id = z$pair_id[[1L]],
          pair_label = z$pair_label[[1L]],
          panel_letter = z$panel_letter[[1L]],
          p_misseg = p_value,
          chunk_index = chunk_index,
          endpoint_indices = paste(selected$endpoint_index, collapse = ","),
          n_unique_endpoint = nrow(selected),
          represented_optimizer_endpoint = sum(selected$endpoint_multiplicity_q10),
          cache_path = file.path(run_paths$cache, sprintf("task_%04d.rds", index)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  tasks <- do.call(rbind, rows)
  rownames(tasks) <- NULL
  f7ft_atomic_write_tsv(tasks, file.path(run_paths$run_root, "finite_time_task_manifest.tsv"))
  tasks
}

f7ft_compute_task <- function(
    task, endpoint_manifest, objective_bundle, contexts, paths, run_paths,
    fingerprint, smoke = FALSE
) {
  f7r_load_response_engine(paths)
  endpoint_indices <- as.integer(strsplit(task$endpoint_indices[[1L]], ",", fixed = TRUE)[[1L]])
  endpoints <- endpoint_manifest$endpoints[
    match(endpoint_indices, endpoint_manifest$endpoints$endpoint_index),
    , drop = FALSE
  ]
  if (anyNA(endpoints$endpoint_index)) stop("Task endpoint lookup failed.")
  o2_values <- f7ft_o2_values(smoke)
  day_values <- f7ft_day_values(smoke)
  sum_array <- array(
    0,
    dim = c(length(f7ft_initial_ploidy()), length(day_values), length(o2_values)),
    dimnames = list(
      initial_ploidy = paste0(f7ft_initial_ploidy(), "N"),
      day = as.character(day_values), O2_pct = as.character(o2_values)
    )
  )
  p_value <- as.numeric(task$p_misseg[[1L]])
  operator_count <- 0L
  maximum_nonfinite <- 0L
  maximum_initial_error <- 0
  maximum_override_error <- 0
  maximum_formula_error <- 0
  for (i in seq_len(nrow(endpoints))) {
    endpoint <- endpoints[i, , drop = FALSE]
    prepared <- f7ft_prepare_endpoint(endpoint, objective_bundle, contexts)
    forced <- figure7_force_p_misseg(prepared$run_params, p_value)
    formula_qc <- figure7_p_misseg_formula_qc(
      prepared$run_params, p_value
    )
    maximum_override_error <- max(
      maximum_override_error, abs(as.numeric(forced$p_misseg) - p_value)
    )
    maximum_formula_error <- max(
      maximum_formula_error, formula_qc$maximum_direct_formula_error
    )
    weight <- as.integer(endpoint$endpoint_multiplicity_q10[[1L]])
    for (o2_index in seq_along(o2_values)) {
      fm <- fixo2_fixed_matrix(
        model_env = globalenv(), cfg = prepared$config,
        run_params = forced, O2 = o2_values[[o2_index]]
      )
      trajectory <- f7ft_expm_daily_mean_ploidy(
        fm$M, fm$ngrid, prepared$config$N_UNIT %||% 22L, day_values
      )
      maximum_nonfinite <- max(maximum_nonfinite, sum(!is.finite(trajectory)))
      maximum_initial_error <- max(
        maximum_initial_error,
        max(abs(trajectory[, 1L] - f7ft_initial_ploidy()), na.rm = TRUE)
      )
      sum_array[, , o2_index] <- sum_array[, , o2_index] + weight * trajectory
      operator_count <- operator_count + 1L
    }
  }
  qc <- data.frame(
    task_id = task$task_id[[1L]],
    model_context = task$model_context[[1L]],
    pair_label = task$pair_label[[1L]],
    p_misseg = p_value,
    n_unique_endpoint = nrow(endpoints),
    represented_optimizer_endpoint = sum(endpoints$endpoint_multiplicity_q10),
    n_operator = operator_count,
    n_o2 = length(o2_values),
    n_day = length(day_values),
    n_initial_ploidy = length(f7ft_initial_ploidy()),
    maximum_nonfinite_per_operator = maximum_nonfinite,
    maximum_day0_abs_error = maximum_initial_error,
    maximum_p_misseg_override_error = maximum_override_error,
    maximum_direct_formula_error = maximum_formula_error,
    passed = maximum_nonfinite == 0L && maximum_initial_error <= 1e-10 &&
      maximum_override_error <= 1e-12 && maximum_formula_error <= 1e-12,
    cache_path = task$cache_path[[1L]],
    stringsAsFactors = FALSE
  )
  f7ft_atomic_save_rds(
    list(
      profile = f7ft_profile(), fingerprint = fingerprint,
      task = task, endpoints = endpoints, o2_values = o2_values,
      day_values = day_values, weighted_sum = sum_array, qc = qc
    ),
    task$cache_path[[1L]], compress = FALSE
  )
  qc
}

f7ft_parallel_lapply <- function(X, FUN, n_core = 1L) {
  n_core <- max(1L, min(as.integer(n_core), length(X)))
  if (n_core == 1L) return(lapply(X, FUN))
  f7r_require_packages(c("future", "future.apply"))
  requested_plan <- tolower(trimws(Sys.getenv(
    "FIGURE7_FINITE_TIME_FUTURE_PLAN", "auto"
  )))
  if (!requested_plan %in% c("auto", "multicore", "multisession")) {
    stop(
      "FIGURE7_FINITE_TIME_FUTURE_PLAN must be auto, multicore, or multisession."
    )
  }
  multicore_supported <- isTRUE(future::supportsMulticore())
  plan_name <- requested_plan
  if (identical(plan_name, "auto")) {
    plan_name <- if (
      identical(.Platform$OS.type, "unix") &&
      !identical(unname(Sys.info()[["sysname"]]), "Darwin") &&
      multicore_supported
    ) {
      "multicore"
    } else {
      "multisession"
    }
  }
  if (identical(plan_name, "multicore") && !multicore_supported) {
    stop("The requested multicore future backend is not supported in this runtime.")
  }
  previous_plan <- future::plan()
  previous_limit <- getOption("parallelly.maxWorkers.localhost")
  previous_size <- getOption("future.globals.maxSize")
  on.exit({
    future::plan(previous_plan)
    options(
      parallelly.maxWorkers.localhost = previous_limit,
      future.globals.maxSize = previous_size
    )
  }, add = TRUE)
  options(
    parallelly.maxWorkers.localhost = max(3L, n_core),
    future.globals.maxSize = 4 * 1024^3
  )
  message(
    "Figure 7 finite-time parallel backend: ", plan_name,
    ", workers=", n_core
  )
  if (identical(plan_name, "multicore")) {
    future::plan(future::multicore, workers = n_core)
  } else {
    future::plan(future::multisession, workers = n_core)
  }
  future.apply::future_lapply(
    X, FUN, future.seed = TRUE, future.scheduling = 1
  )
}

f7ft_aggregate_tasks <- function(
    tasks, endpoint_manifest, run_paths, fingerprint, smoke = FALSE
) {
  objects <- lapply(tasks$cache_path, readRDS)
  if (any(!vapply(objects, function(x) identical(x$fingerprint, fingerprint), logical(1L)))) {
    stop("Finite-time task caches do not share the current input fingerprint.")
  }
  task_qc <- do.call(rbind, lapply(objects, `[[`, "qc"))
  if (!all(task_qc$passed)) stop("One or more finite-time tasks failed QC.")
  f7ft_atomic_write_tsv(
    task_qc, file.path(run_paths$run_root, "finite_time_task_validation.tsv")
  )

  panel_objects <- list()
  panel_rows <- unique(tasks[, c(
    "model_context", "pair_id", "pair_label", "panel_letter"
  )])
  panel_rows <- panel_rows[order(panel_rows$panel_letter), , drop = FALSE]
  for (panel_index in seq_len(nrow(panel_rows))) {
    panel <- panel_rows[panel_index, , drop = FALSE]
    p_values <- sort(unique(tasks$p_misseg[
      tasks$model_context == panel$model_context[[1L]] &
        tasks$pair_label == panel$pair_label[[1L]]
    ]))
    arrays <- vector("list", length(p_values))
    total_weights <- integer(length(p_values))
    for (p_index in seq_along(p_values)) {
      keep <- which(
        tasks$model_context == panel$model_context[[1L]] &
          tasks$pair_label == panel$pair_label[[1L]] &
          abs(tasks$p_misseg - p_values[[p_index]]) < 1e-12
      )
      selected <- objects[keep]
      weighted_sum <- Reduce(`+`, lapply(selected, `[[`, "weighted_sum"))
      total_weight <- sum(vapply(
        selected,
        function(x) sum(x$endpoints$endpoint_multiplicity_q10), integer(1L)
      ))
      expected_weight <- if (isTRUE(smoke)) {
        sum(endpoint_manifest$endpoints$endpoint_multiplicity_q10[
          endpoint_manifest$endpoints$model_context == panel$model_context[[1L]] &
            endpoint_manifest$endpoints$pair_label == panel$pair_label[[1L]]
        ][1L])
      } else 50L
      if (total_weight != expected_weight) {
        stop("Unexpected optimizer-endpoint weight for panel ", panel$panel_letter[[1L]])
      }
      arrays[[p_index]] <- weighted_sum / total_weight
      total_weights[[p_index]] <- total_weight
    }
    values <- simplify2array(arrays)
    dimnames(values)[[4L]] <- format(p_values, scientific = FALSE, trim = TRUE)
    panel_object <- list(
      profile = f7ft_profile(), fingerprint = fingerprint,
      panel_letter = panel$panel_letter[[1L]],
      pair_id = panel$pair_id[[1L]], pair_label = panel$pair_label[[1L]],
      model_context = panel$model_context[[1L]],
      initial_ploidy = f7ft_initial_ploidy(),
      day_values = f7ft_day_values(smoke),
      o2_values = f7ft_o2_values(smoke),
      p_misseg = p_values,
      optimizer_endpoint_weight = total_weights,
      mean_ploidy = values
    )
    panel_path <- file.path(
      run_paths$run_root,
      paste0("finite_time_panel_", tolower(panel$panel_letter[[1L]]), ".rds")
    )
    f7ft_atomic_save_rds(panel_object, panel_path, compress = "gzip")
    panel_objects[[panel$panel_letter[[1L]]]] <- panel_object
  }

  validation_rows <- lapply(panel_objects, function(object) {
    values <- object$mean_ploidy
    day0 <- values[, 1L, , , drop = FALSE]
    expected_n_p <- if (isTRUE(smoke)) 2L else 5L
    expected_n_o2 <- if (isTRUE(smoke)) 3L else 201L
    expected_n_day <- if (isTRUE(smoke)) 11L else 1001L
    expected <- array(
      object$initial_ploidy,
      dim = dim(day0), dimnames = dimnames(day0)
    )
    data.frame(
      panel_letter = object$panel_letter,
      model_context = object$model_context,
      pair_label = object$pair_label,
      n_initial_ploidy = length(object$initial_ploidy),
      n_day = length(object$day_values),
      n_o2 = length(object$o2_values),
      n_o2_0_to_0p5 = sum(
        object$o2_values >= 0 & object$o2_values <= 0.5 + 1e-12
      ),
      oxygen_grid_exact = isTRUE(all.equal(
        as.numeric(object$o2_values), as.numeric(f7ft_o2_values(smoke))
      )),
      n_p_misseg = length(object$p_misseg),
      minimum_mean_ploidy = min(values, na.rm = TRUE),
      maximum_mean_ploidy = max(values, na.rm = TRUE),
      maximum_day0_abs_error = max(abs(day0 - expected), na.rm = TRUE),
      all_finite = all(is.finite(values)),
      weights_match_expected = all(object$optimizer_endpoint_weight ==
        if (isTRUE(smoke)) object$optimizer_endpoint_weight[[1L]] else 50L),
      passed = all(is.finite(values)) && min(values, na.rm = TRUE) >= 1 - 1e-8 &&
        max(values, na.rm = TRUE) <= 7 + 1e-8 &&
        max(abs(day0 - expected), na.rm = TRUE) <= 1e-10 &&
        length(object$initial_ploidy) == 5L &&
        length(object$p_misseg) == expected_n_p &&
        length(object$o2_values) == expected_n_o2 &&
        isTRUE(all.equal(
          as.numeric(object$o2_values), as.numeric(f7ft_o2_values(smoke))
        )) &&
        length(object$day_values) == expected_n_day,
      stringsAsFactors = FALSE
    )
  })
  validation <- do.call(rbind, validation_rows)
  f7ft_atomic_write_tsv(
    validation, file.path(run_paths$run_root, "finite_time_panel_validation.tsv")
  )
  if (nrow(validation) != 4L || !all(validation$passed)) {
    stop("Figure 7 finite-time panel validation failed.")
  }
  invisible(list(panels = panel_objects, validation = validation, task_qc = task_qc))
}

f7ft_matrix_power_apply <- function(matrix, state, exponent) {
  exponent <- as.integer(exponent)
  if (exponent < 0L) stop("Matrix exponent must be non-negative.")
  if (exponent == 0L) return(f7ft_normalize_columns(state))
  power <- as.matrix(matrix)
  result <- as.matrix(state)
  n <- exponent
  while (n > 0L) {
    if (n %% 2L == 1L) {
      result <- power %*% result
      result <- f7ft_rescale_columns(result)
    }
    n <- n %/% 2L
    if (n > 0L) {
      power <- power %*% power
      scale <- max(abs(power), na.rm = TRUE)
      if (!is.finite(scale) || scale <= 0) stop("Non-finite matrix power.")
      power <- power / scale
    }
  }
  f7ft_normalize_columns(result)
}

f7ft_expm_anchor_states <- function(M, initial, time_values) {
  time_values <- sort(unique(as.numeric(time_values)))
  state <- initial
  current <- 0
  output <- vector("list", length(time_values))
  for (i in seq_along(time_values)) {
    delta <- time_values[[i]] - current
    if (delta < -1e-12) stop("Anchor times must be sorted.")
    if (delta > 1e-12) {
      state <- f7ft_normalize_columns(as.matrix(Matrix::expm(M * delta)) %*% state)
    }
    output[[i]] <- state
    current <- time_values[[i]]
  }
  output
}

f7ft_eigen_anchor_states <- function(M, initial, time_values) {
  eig <- eigen(as.matrix(M), only.values = FALSE)
  coefficients <- solve(eig$vectors, initial)
  lambda_ref <- max(Re(eig$values))
  states <- lapply(time_values, function(tt) {
    weights <- exp((eig$values - lambda_ref) * tt)
    f7ft_normalize_columns(eig$vectors %*% (weights * coefficients))
  })
  order_real <- order(Re(eig$values), decreasing = TRUE)
  lambda1 <- Re(eig$values[order_real[[1L]]])
  lambda2 <- Re(eig$values[order_real[[min(2L, length(order_real))]]])
  dominant <- Re(eig$vectors[, order_real[[1L]], drop = FALSE])
  if (sum(dominant, na.rm = TRUE) < 0) dominant <- -dominant
  list(
    states = states,
    dominant = f7ft_normalize_columns(dominant),
    spectral_gap = lambda1 - lambda2,
    eigenvector_condition_number = tryCatch(kappa(eig$vectors), error = function(e) NA_real_)
  )
}

f7ft_euler_anchor_states <- function(M, initial, time_values, dt = 0.1) {
  steps <- as.integer(round(as.numeric(time_values) / dt))
  if (any(abs(steps * dt - time_values) > 1e-8)) {
    stop("Euler anchor times must be exact multiples of dt.")
  }
  transition <- diag(nrow(M)) + dt * as.matrix(M)
  lapply(steps, function(step) f7ft_matrix_power_apply(transition, initial, step))
}

f7ft_state_mean_ploidy <- function(state, ngrid, n_unit) {
  as.numeric(crossprod(as.numeric(ngrid) / as.numeric(n_unit), state))
}

f7ft_compute_diagnostic_endpoint <- function(
    endpoint, objective_bundle, contexts, paths, run_paths, fingerprint
) {
  f7r_load_response_engine(paths)
  prepared <- f7ft_prepare_endpoint(endpoint, objective_bundle, contexts)
  time_values <- f7ft_diagnostic_days()
  rows <- list()
  row_index <- 0L
  for (p_value in f7ft_p_values()) {
    forced <- figure7_force_p_misseg(prepared$run_params, p_value)
    for (o2 in f7ft_diagnostic_o2()) {
      fm <- fixo2_fixed_matrix(
        model_env = globalenv(), cfg = prepared$config,
        run_params = forced, O2 = o2
      )
      n_unit <- prepared$config$N_UNIT %||% 22L
      initial <- f7ft_initial_matrix(fm$ngrid, n_unit)
      expm_states <- f7ft_expm_anchor_states(fm$M, initial, time_values)
      eigen_bundle <- f7ft_eigen_anchor_states(fm$M, initial, time_values)
      euler_states <- f7ft_euler_anchor_states(fm$M, initial, time_values, dt = 0.1)
      steady <- f7ft_state_mean_ploidy(eigen_bundle$dominant, fm$ngrid, n_unit)[[1L]]
      for (time_index in seq_along(time_values)) {
        expm_mean <- f7ft_state_mean_ploidy(
          expm_states[[time_index]], fm$ngrid, n_unit
        )
        eigen_mean <- f7ft_state_mean_ploidy(
          eigen_bundle$states[[time_index]], fm$ngrid, n_unit
        )
        euler_mean <- f7ft_state_mean_ploidy(
          euler_states[[time_index]], fm$ngrid, n_unit
        )
        for (initial_index in seq_along(f7ft_initial_ploidy())) {
          row_index <- row_index + 1L
          rows[[row_index]] <- data.frame(
            model_context = endpoint$model_context[[1L]],
            pair_id = endpoint$pair_id[[1L]],
            pair_label = endpoint$pair_label[[1L]],
            panel_letter = endpoint$panel_letter[[1L]],
            endpoint_group = endpoint$endpoint_group[[1L]],
            representative_seed_number = endpoint$representative_seed_number[[1L]],
            endpoint_multiplicity_q10 = endpoint$endpoint_multiplicity_q10[[1L]],
            p_misseg = p_value,
            O2_pct = o2,
            initial_ploidy = f7ft_initial_ploidy()[[initial_index]],
            day = time_values[[time_index]],
            eigen_mean_ploidy = eigen_mean[[initial_index]],
            expm_mean_ploidy = expm_mean[[initial_index]],
            euler_mean_ploidy = euler_mean[[initial_index]],
            steady_mean_ploidy = steady,
            spectral_gap = eigen_bundle$spectral_gap,
            eigenvector_condition_number = eigen_bundle$eigenvector_condition_number,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  diagnostics <- do.call(rbind, rows)
  context_tag <- if (endpoint$model_context[[1L]] == "in vitro") "vitro" else "vivo"
  cache_path <- file.path(
    run_paths$diagnostics_cache,
    paste0("diagnostic_", context_tag, "_", endpoint$endpoint_group[[1L]], ".rds")
  )
  qc <- data.frame(
    model_context = endpoint$model_context[[1L]],
    pair_label = endpoint$pair_label[[1L]],
    endpoint_group = endpoint$endpoint_group[[1L]],
    represented_optimizer_endpoint = endpoint$endpoint_multiplicity_q10[[1L]],
    n_row = nrow(diagnostics),
    core_all_finite = all(vapply(
      diagnostics[c("expm_mean_ploidy", "euler_mean_ploidy",
                    "steady_mean_ploidy", "spectral_gap",
                    "eigenvector_condition_number")],
      function(x) all(is.finite(x)), logical(1L)
    )),
    eigen_finite_rows = sum(is.finite(diagnostics$eigen_mean_ploidy)),
    eigen_nonfinite_rows = sum(!is.finite(diagnostics$eigen_mean_ploidy)),
    eigen_finite_fraction = mean(is.finite(diagnostics$eigen_mean_ploidy)),
    minimum_required_eigen_finite_fraction = 0.98,
    maximum_abs_eigen_expm = max(abs(
      diagnostics$eigen_mean_ploidy - diagnostics$expm_mean_ploidy
    ), na.rm = TRUE),
    cache_path = cache_path,
    stringsAsFactors = FALSE
  )
  qc$passed <- qc$core_all_finite &&
    qc$eigen_finite_fraction >= qc$minimum_required_eigen_finite_fraction
  f7ft_atomic_save_rds(
    list(
      profile = f7ft_profile(), fingerprint = fingerprint,
      endpoint = endpoint, diagnostics = diagnostics, qc = qc
    ),
    cache_path, compress = "gzip"
  )
  qc
}

f7ft_weighted_quantile <- function(x, w, probability) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]; w <- w[keep]
  if (!length(x)) return(NA_real_)
  order_index <- order(x)
  x <- x[order_index]; w <- w[order_index]
  x[[which(cumsum(w) / sum(w) >= probability)[[1L]]]]
}

f7ft_error_metrics <- function(data, reference, estimate, group_columns) {
  groups <- split(data, interaction(data[group_columns], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(z) {
    total_rows <- nrow(z)
    total_weight_all <- sum(
      z$endpoint_multiplicity_q10[is.finite(z$endpoint_multiplicity_q10)],
      na.rm = TRUE
    )
    keep <- is.finite(z[[reference]]) & is.finite(z[[estimate]]) &
      is.finite(z$endpoint_multiplicity_q10) & z$endpoint_multiplicity_q10 > 0
    z <- z[keep, , drop = FALSE]
    if (!nrow(z)) stop("No finite matched diagnostic rows for error metrics.")
    residual <- z[[estimate]] - z[[reference]]
    weight <- z$endpoint_multiplicity_q10
    total_weight <- sum(weight)
    data.frame(
      z[1L, group_columns, drop = FALSE],
      comparison = paste0(estimate, "_vs_", reference),
      n_unique_row = nrow(z),
      n_total_row = total_rows,
      finite_row_fraction = nrow(z) / total_rows,
      weighted_n = total_weight,
      weighted_n_total = total_weight_all,
      weighted_coverage = total_weight / total_weight_all,
      bias = sum(weight * residual) / total_weight,
      rmse = sqrt(sum(weight * residual^2) / total_weight),
      q95_absolute_error = f7ft_weighted_quantile(abs(residual), weight, 0.95),
      maximum_absolute_error = max(abs(residual), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

f7ft_euler_step_sensitivity <- function(
    endpoint_manifest, objective_bundle, contexts, paths, run_paths, fingerprint
) {
  representatives <- do.call(rbind, lapply(
    split(
      endpoint_manifest$endpoints,
      interaction(
        endpoint_manifest$endpoints$model_context,
        endpoint_manifest$endpoints$pair_label,
        drop = TRUE, lex.order = TRUE
      )
    ),
    function(z) z[which.min(z$representative_objective_rank), , drop = FALSE]
  ))
  rows <- list()
  row_index <- 0L
  for (endpoint_index in seq_len(nrow(representatives))) {
    endpoint <- representatives[endpoint_index, , drop = FALSE]
    prepared <- f7ft_prepare_endpoint(endpoint, objective_bundle, contexts)
    for (p_value in c(0.01, 0.30)) {
      forced <- figure7_force_p_misseg(prepared$run_params, p_value)
      for (o2 in c(0, 1, 5)) {
        fm <- fixo2_fixed_matrix(
          model_env = globalenv(), cfg = prepared$config,
          run_params = forced, O2 = o2
        )
        n_unit <- prepared$config$N_UNIT %||% 22L
        initial <- f7ft_initial_matrix(fm$ngrid, n_unit)
        time_values <- f7ft_diagnostic_days()
        expm_states <- f7ft_expm_anchor_states(fm$M, initial, time_values)
        euler_0p1 <- f7ft_euler_anchor_states(
          fm$M, initial, time_values, dt = 0.1
        )
        euler_0p05 <- f7ft_euler_anchor_states(
          fm$M, initial, time_values, dt = 0.05
        )
        for (time_index in seq_along(time_values)) {
          expm_mean <- f7ft_state_mean_ploidy(
            expm_states[[time_index]], fm$ngrid, n_unit
          )
          mean_0p1 <- f7ft_state_mean_ploidy(
            euler_0p1[[time_index]], fm$ngrid, n_unit
          )
          mean_0p05 <- f7ft_state_mean_ploidy(
            euler_0p05[[time_index]], fm$ngrid, n_unit
          )
          for (initial_index in seq_along(f7ft_initial_ploidy())) {
            row_index <- row_index + 1L
            rows[[row_index]] <- data.frame(
              panel_letter = endpoint$panel_letter[[1L]],
              model_context = endpoint$model_context[[1L]],
              pair_label = endpoint$pair_label[[1L]],
              endpoint_group = endpoint$endpoint_group[[1L]],
              p_misseg = p_value,
              O2_pct = o2,
              initial_ploidy = f7ft_initial_ploidy()[[initial_index]],
              day = time_values[[time_index]],
              expm_mean_ploidy = expm_mean[[initial_index]],
              euler_dt_0p1_mean_ploidy = mean_0p1[[initial_index]],
              euler_dt_0p05_mean_ploidy = mean_0p05[[initial_index]],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  sensitivity <- do.call(rbind, rows)
  sensitivity$abs_euler_step_difference <- abs(
    sensitivity$euler_dt_0p1_mean_ploidy -
      sensitivity$euler_dt_0p05_mean_ploidy
  )
  sensitivity$abs_euler_0p05_expm_difference <- abs(
    sensitivity$euler_dt_0p05_mean_ploidy - sensitivity$expm_mean_ploidy
  )
  summary <- do.call(rbind, lapply(
    split(sensitivity, sensitivity$panel_letter),
    function(z) data.frame(
      panel_letter = z$panel_letter[[1L]],
      model_context = z$model_context[[1L]],
      pair_label = z$pair_label[[1L]],
      n_row = nrow(z),
      maximum_abs_dt_0p1_minus_0p05 = max(
        z$abs_euler_step_difference, na.rm = TRUE
      ),
      q95_abs_dt_0p1_minus_0p05 = unname(stats::quantile(
        z$abs_euler_step_difference, 0.95, na.rm = TRUE
      )),
      maximum_abs_euler_0p05_minus_expm = max(
        z$abs_euler_0p05_expm_difference, na.rm = TRUE
      ),
      all_finite = all(vapply(
        z[c(
          "expm_mean_ploidy", "euler_dt_0p1_mean_ploidy",
          "euler_dt_0p05_mean_ploidy"
        )],
        function(x) all(is.finite(x)), logical(1L)
      )),
      stringsAsFactors = FALSE
    )
  ))
  summary$passed <- summary$all_finite
  paths_out <- c(
    rows = f7ft_atomic_save_rds(
      sensitivity,
      file.path(run_paths$run_root, "finite_time_euler_step_sensitivity_rows.rds"),
      compress = "gzip"
    ),
    summary = f7ft_atomic_write_tsv(
      summary,
      file.path(run_paths$run_root, "finite_time_euler_step_sensitivity.tsv")
    )
  )
  if (nrow(summary) != 4L || !all(summary$passed)) {
    stop("Finite-time Euler step-size sensitivity validation failed.")
  }
  list(data = sensitivity, summary = summary, paths = paths_out)
}

f7ft_compute_diagnostics <- function(
    endpoint_manifest, objective_bundle, contexts, paths, run_paths,
    fingerprint, n_core = 1L
) {
  endpoints <- endpoint_manifest$endpoints
  endpoint_list <- split(endpoints, seq_len(nrow(endpoints)))
  compute_one <- function(endpoint) {
    tryCatch(
      f7ft_compute_diagnostic_endpoint(
        endpoint, objective_bundle, contexts, paths, run_paths, fingerprint
      ),
      error = function(e) structure(
        list(
          endpoint_group = endpoint$endpoint_group[[1L]],
          message = conditionMessage(e)
        ),
        class = "f7ft_error"
      )
    )
  }
  results <- f7ft_parallel_lapply(endpoint_list, compute_one, n_core = n_core)
  failed <- vapply(results, inherits, logical(1L), "f7ft_error")
  if (any(failed)) {
    stop(
      "Figure 7 finite-time diagnostic failures: ",
      paste(vapply(results[failed], `[[`, character(1L), "message"), collapse = "; ")
    )
  }
  qc <- do.call(rbind, results)
  f7ft_atomic_write_tsv(
    qc, file.path(run_paths$run_root, "finite_time_diagnostic_endpoint_validation.tsv")
  )
  if (!all(qc$passed)) stop("Figure 7 finite-time diagnostic endpoint QC failed.")
  objects <- lapply(qc$cache_path, readRDS)
  diagnostics <- do.call(rbind, lapply(objects, `[[`, "diagnostics"))
  diagnostics <- diagnostics[order(
    diagnostics$panel_letter, diagnostics$endpoint_group,
    diagnostics$p_misseg, diagnostics$O2_pct,
    diagnostics$initial_ploidy, diagnostics$day
  ), , drop = FALSE]
  raw_path <- f7ft_atomic_save_rds(
    diagnostics,
    file.path(run_paths$run_root, "finite_time_diagnostic_rows.rds"),
    compress = "gzip"
  )
  metrics <- rbind(
    f7ft_error_metrics(
      diagnostics, "euler_mean_ploidy", "eigen_mean_ploidy",
      c("panel_letter", "model_context", "pair_label")
    ),
    f7ft_error_metrics(
      diagnostics, "euler_mean_ploidy", "expm_mean_ploidy",
      c("panel_letter", "model_context", "pair_label")
    ),
    f7ft_error_metrics(
      diagnostics, "expm_mean_ploidy", "eigen_mean_ploidy",
      c("panel_letter", "model_context", "pair_label")
    )
  )
  metrics_path <- f7ft_atomic_write_tsv(
    metrics, file.path(run_paths$run_root, "finite_time_diagnostic_error_summary.tsv")
  )
  euler_sensitivity <- f7ft_euler_step_sensitivity(
    endpoint_manifest, objective_bundle, contexts, paths, run_paths, fingerprint
  )
  list(data = diagnostics, qc = qc, metrics = metrics,
       euler_sensitivity = euler_sensitivity,
       paths = c(raw = raw_path, metrics = metrics_path,
                 euler_sensitivity$paths))
}

f7ft_chart_contract <- function(run_paths) {
  contract <- data.frame(
    artifact = c(
      "Figure 7C-F", "Supplementary Figure 7-5",
      "Supplementary Figure 7-6", "Supplementary Figure 7-7"
    ),
    analytical_question = c(
      "How do finite-time mean-ploidy trajectories vary with initial ploidy, fixed p_misseg, oxygen, family, and model context?",
      "Do full-eigen and expm finite-time solutions agree with an independent Euler integration on the q10 endpoint ensemble?",
      "How closely do the full-eigen and expm finite-time analytical solutions agree?",
      "At what time scale does the expm finite-time composition approach the dominant-eigenvector attractor?"
    ),
    chart_family = c(
      "five-by-five faceted heatmap", "calibration scatter matrix",
      "calibration scatter plus Bland-Altman residual", "time-faceted calibration scatter"
    ),
    data_grain = c(
      "4 panels x 5 initial ploidies x 5 fixed probabilities x 201 oxygen values x 1001 days; arithmetic mean over 50 q10 optimizer endpoints",
      "q10 unique parameter endpoints with optimizer-seed multiplicity weights across anchor oxygen, probability, initial-ploidy, and time grids",
      "same matched rows as Supplementary Figure 7-5",
      "same matched rows at days 25, 100, 500, and 1000"
    ),
    primary_encoding = c(
      "day on x, fixed oxygen on y, mean finite-time ploidy as Figure 7A-matched log color",
      "analytical method on x, Euler numerical integration on y, q10 endpoint density and identity line",
      "eigen on x, expm on y; residual against their mean in the lower block",
      "dominant-attractor ploidy on x and expm finite-time ploidy on y"
    ),
    caveat = c(
      "Optimizer endpoints are numerical solutions rather than biological replicates; oxygen and p_misseg are fixed within each trajectory",
      "Euler is a deterministic numerical-integration reference, not a stochastic biological simulation",
      "Near-degenerate or poorly conditioned eigenvectors can affect the full-eigen reconstruction",
      "Agreement at long time is conditional on a fixed environment and the fitted operator"
    ),
    stringsAsFactors = FALSE
  )
  f7ft_atomic_write_tsv(
    contract, file.path(run_paths$run_root, "figure7_finite_time_chart_contract.tsv")
  )
}

f7ft_publish_current <- function(paths, run_paths, fingerprint) {
  relative <- file.path("finite_time_q10_runs", run_paths$run_id)
  pointer <- data.frame(
    run_id = run_paths$run_id,
    relative_run_path = relative,
    profile = f7ft_profile(),
    fingerprint = fingerprint,
    published_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  f7ft_atomic_write_tsv(pointer, run_paths$current)
}

f7ft_data <- function(
    workspace_root = f7r_find_workspace_root(), n_core = 1L,
    run_id = f7ft_resolve_run_id(), smoke = FALSE,
    publish_current = !isTRUE(smoke), compute_diagnostics = !isTRUE(smoke)
) {
  Sys.setenv(
    KMP_USE_SHM = "0", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  f7r_require_packages(c("Matrix", "Rcpp", "data.table"))
  paths <- f7r_paths(workspace_root)
  run_paths <- f7ft_paths(paths, run_id = run_id, create = TRUE)
  if (length(list.files(run_paths$cache, all.files = TRUE, no.. = TRUE)) ||
      length(list.files(run_paths$diagnostics_cache, all.files = TRUE, no.. = TRUE))) {
    stop(
      "Fresh Figure 7 finite-time run directory is not empty: ",
      run_paths$run_root,
      ". Choose a new run id rather than reusing old caches."
    )
  }
  f7r_load_response_engine(paths)
  objective_bundle <- f7r_objective_selection(paths)
  endpoint_manifest <- f7ft_build_endpoint_manifest(
    paths, objective_bundle, run_paths
  )
  fingerprint <- f7ft_model_fingerprint(paths, endpoint_manifest)
  f7ft_chart_contract(run_paths)
  contexts <- lapply(
    unique(endpoint_manifest$endpoints$pair_id),
    f7r_pair_model_context,
    selected = objective_bundle$selected, paths = paths
  )
  names(contexts) <- unique(endpoint_manifest$endpoints$pair_id)
  tasks <- f7ft_task_manifest(
    endpoint_manifest, run_paths, smoke = smoke, chunk_size = 2L
  )
  task_list <- split(tasks, seq_len(nrow(tasks)))
  compute_one <- function(task) {
    tryCatch(
      f7ft_compute_task(
        task, endpoint_manifest, objective_bundle, contexts,
        paths, run_paths, fingerprint, smoke = smoke
      ),
      error = function(e) structure(
        list(task_id = task$task_id[[1L]], message = conditionMessage(e)),
        class = "f7ft_error"
      )
    )
  }
  message(
    "Figure 7 finite-time q10: ", nrow(tasks), " tasks, workers=",
    min(as.integer(n_core), nrow(tasks)), ", run_id=", run_paths$run_id
  )
  results <- f7ft_parallel_lapply(task_list, compute_one, n_core = n_core)
  failed <- vapply(results, inherits, logical(1L), "f7ft_error")
  if (any(failed)) {
    stop(
      "Figure 7 finite-time task failures: ",
      paste(vapply(results[failed], `[[`, character(1L), "message"), collapse = "; ")
    )
  }
  aggregation <- f7ft_aggregate_tasks(
    tasks, endpoint_manifest, run_paths, fingerprint, smoke = smoke
  )
  diagnostics <- NULL
  if (isTRUE(compute_diagnostics)) {
    diagnostics <- f7ft_compute_diagnostics(
      endpoint_manifest, objective_bundle, contexts, paths, run_paths,
      fingerprint, n_core = n_core
    )
  }
  provenance <- data.frame(
    key = c(
      "profile", "run_id", "fingerprint", "workspace_root", "model_code_root",
      "model_source_fingerprint", "joint_result_root", "n_core",
      "smoke", "optimizer_endpoint_per_panel", "aggregation"
    ),
    value = c(
      f7ft_profile(), run_paths$run_id, fingerprint,
      normalizePath(paths$root, mustWork = TRUE),
      normalizePath(paths$oxygen_code, mustWork = TRUE),
      f7r_model_source_fingerprint(paths),
      Sys.getenv("FIGURE_JOINT_RESULT_ROOT", unset = "not-recorded"),
      as.character(n_core), as.character(smoke),
      if (isTRUE(smoke)) "one unique endpoint per panel" else "50",
      "arithmetic mean with optimizer-seed multiplicity restored"
    ),
    stringsAsFactors = FALSE
  )
  provenance_path <- f7ft_atomic_write_tsv(
    provenance, file.path(run_paths$run_root, "finite_time_run_provenance.tsv")
  )
  if (isTRUE(publish_current)) {
    f7ft_publish_current(paths, run_paths, fingerprint)
  }
  invisible(list(
    paths = run_paths, endpoint_manifest = endpoint_manifest,
    tasks = tasks, aggregation = aggregation, diagnostics = diagnostics,
    fingerprint = fingerprint, provenance = provenance_path
  ))
}
