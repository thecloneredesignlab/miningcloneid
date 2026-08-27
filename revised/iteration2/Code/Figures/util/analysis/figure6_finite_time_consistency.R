#!/usr/bin/env Rscript

# Supplementary Figure 6-5: numerical consistency between finite-time
# propagation and the fixed-environment dominant eigenvector.  This is a
# post-fit diagnostic.  It does not refit parameters or reinterpret optimizer
# endpoints as biological replicates.

options(stringsAsFactors = FALSE, warn = 1)

s65_o2_values <- function() c(0, 1, 5, 20)
s65_separate_o2_values <- function() c(s65_o2_values(), 20.5)
s65_p_values <- function() c(0.01, 0.10, 0.20, 0.30)
s65_u_values <- function() {
  c(
    0, 0.1, 0.25, 0.5, 1, 2, 3, 5, 7.5,
    10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 30, 40
  )
}
s65_solver_u_values <- function() c(0, 1, 5, 20, 40)
s65_terminal_u <- function() 40
s65_seeds <- function() c(10L, 132L)
s65_shared_parameters <- function() f6r_shared_parameters()

s65_thresholds <- function() {
  c(
    requested_cin_error = 1e-8,
    solver_mean_ploidy_error = 1e-7,
    solver_total_variation = 1e-7,
    terminal_mean_ploidy_error = 1e-4,
    terminal_total_variation = 1e-4,
    initial_state_terminal_difference = 1e-4,
    normalized_mass_error = 1e-12,
    dominant_eigenpair_residual = 1e-10
  )
}

s65_paths <- function(workspace_root = f6r_find_workspace_root()) {
  base <- f6r_paths(workspace_root)
  data <- file.path(base$root, "data", "Figures", "Supp_Figure6_5")
  list(
    base = base,
    data = data,
    panels = file.path(data, "panels"),
    output_base = "supp_fig6-5_finite_time_eigen_consistency",
    output_png = file.path(
      base$figures, "supp_fig6-5_finite_time_eigen_consistency.png"
    ),
    output_pdf = file.path(
      base$figures, "supp_fig6-5_finite_time_eigen_consistency.pdf"
    ),
    manuscript_png = file.path(
      base$manuscript_figures,
      "supp_fig6-5_finite_time_eigen_consistency.png"
    ),
    manuscript_pdf = file.path(
      base$manuscript_figures,
      "supp_fig6-5_finite_time_eigen_consistency.pdf"
    )
  )
}

s65_require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(label, " is missing column(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

s65_normalize_state <- function(x) {
  value <- if (is.complex(x)) Re(x) else as.numeric(x)
  value <- as.numeric(value)
  value[!is.finite(value)] <- NA_real_
  value[value < 0 & value >= -1e-10] <- 0
  if (any(!is.finite(value)) || any(value < -1e-10)) {
    return(rep(NA_real_, length(value)))
  }
  value <- pmax(value, 0)
  total <- sum(value)
  if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(value)))
  value / total
}

s65_mean_ploidy <- function(state, ngrid, n_unit) {
  sum(as.numeric(state) * as.numeric(ngrid)) / as.numeric(n_unit)
}

s65_total_variation <- function(a, b) {
  0.5 * sum(abs(as.numeric(a) - as.numeric(b)))
}

s65_path_tokens <- function(x) strsplit(as.character(x), "_", fixed = TRUE)

s65_path_length <- function(x) lengths(s65_path_tokens(x))

s65_is_deprived_path <- function(x) {
  vapply(s65_path_tokens(x), function(tokens) {
    oxygen <- suppressWarnings(as.numeric(tokens))
    any(is.finite(oxygen) & oxygen < 20.5)
  }, logical(1L))
}

s65_experiment_duration <- function(lineage, cohort) {
  s65_require_columns(
    lineage,
    c("segment_id", "cohort", "selected_day"),
    "in-vitro lineage summary"
  )
  candidate <- lineage[
    lineage$cohort == cohort & s65_is_deprived_path(lineage$segment_id),
    , drop = FALSE
  ]
  if (!nrow(candidate)) stop("No oxygen-deprived lineage for cohort ", cohort)
  lengths <- s65_path_length(candidate$segment_id)
  terminal <- candidate$segment_id[[which.max(lengths)]]
  tokens <- s65_path_tokens(terminal)[[1L]]
  prefixes <- vapply(seq_along(tokens), function(index) {
    paste(tokens[seq_len(index)], collapse = "_")
  }, character(1L))
  selected_days <- vapply(prefixes, function(prefix) {
    values <- unique(lineage$selected_day[
      lineage$cohort == cohort & lineage$segment_id == prefix
    ])
    values <- as.numeric(values[is.finite(values)])
    if (length(values) != 1L) {
      stop(
        "Expected one selected day for ", cohort,
        " lineage prefix ", prefix
      )
    }
    values[[1L]]
  }, numeric(1L))
  data.frame(
    cohort = cohort,
    terminal_passage_count = length(tokens),
    experimental_duration_days = sum(selected_days),
    stringsAsFactors = FALSE
  )
}

s65_read_objective <- function(path) {
  x <- f6r_read_tsv(path)
  s65_require_columns(x, c("metric", "value"), "fit summary")
  keys <- c("objective_total", "optimizer_local_objective", "objective")
  for (key in keys) {
    value <- suppressWarnings(as.numeric(x$value[x$metric == key]))
    if (length(value) == 1L && is.finite(value)) return(value)
  }
  stop("Could not recover a finite fitted objective from: ", path)
}

s65_load_separate_seed <- function(seed_number, result_root) {
  seed_id <- paste0("seed", as.integer(seed_number))
  seed_dir <- file.path(result_root, seed_id)
  files <- c(
    fit_result = file.path(seed_dir, "fit_result.rds"),
    fit_summary = file.path(seed_dir, "fit_summary.tsv"),
    distribution = file.path(seed_dir, "invitro_distribution_summary.tsv"),
    lineage = file.path(seed_dir, "invitro_lineage_summary.tsv"),
    best_params = file.path(seed_dir, "best_params.tsv")
  )
  f6r_require_files(files, paste0(seed_id, " separate in-vitro fit"))
  fit <- readRDS(files[["fit_result"]])
  if (!is.list(fit) || !is.list(fit$cfg) ||
      !is.list(fit$best_components) || !is.list(fit$best_params)) {
    stop("Unexpected fit_result structure for ", seed_id)
  }
  shared <- s65_shared_parameters()
  values <- suppressWarnings(as.numeric(unlist(
    fit$best_params[shared], use.names = FALSE
  )))
  names(values) <- shared
  if (length(values) != length(shared) || any(!is.finite(values))) {
    stop("Incomplete shared fitted parameters for ", seed_id)
  }
  model_core <- fit$best_components$run_2N$model_core
  ngrid <- as.integer(model_core$grid_pre)
  init_2n <- s65_normalize_state(model_core$init_state_2N)
  init_4n <- s65_normalize_state(model_core$init_state_4N)
  if (length(ngrid) != length(init_2n) || length(ngrid) != length(init_4n) ||
      any(!is.finite(c(init_2n, init_4n)))) {
    stop("Invalid fitted initial chromosome-state distributions for ", seed_id)
  }
  distribution <- f6r_read_tsv(files[["distribution"]])
  lineage <- f6r_read_tsv(files[["lineage"]])
  list(
    seed_id = seed_id,
    seed_number = as.integer(seed_number),
    objective = s65_read_objective(files[["fit_summary"]]),
    params = values,
    config = fit$cfg,
    ngrid = ngrid,
    initial_states = list(`2N start` = init_2n, `4N start` = init_4n),
    distribution = distribution,
    lineage = lineage,
    files = normalizePath(files, mustWork = TRUE)
  )
}

s65_experimental_distribution <- function(seed_bundle) {
  x <- seed_bundle$distribution
  s65_require_columns(
    x,
    c("segment_id", "cohort", "N", "fraction"),
    paste0(seed_bundle$seed_id, " distribution summary")
  )
  x <- x[
    x$cohort == "4N" & s65_is_deprived_path(x$segment_id) &
      is.finite(x$N) & is.finite(x$fraction) & x$fraction >= 0,
    , drop = FALSE
  ]
  x$passage_sequence <- s65_path_length(x$segment_id)
  x <- stats::aggregate(
    fraction ~ passage_sequence + N,
    data = x,
    FUN = mean,
    na.rm = TRUE
  )
  mass <- stats::aggregate(
    fraction ~ passage_sequence, data = x, FUN = sum, na.rm = TRUE
  )
  x$fraction <- x$fraction /
    mass$fraction[match(x$passage_sequence, mass$passage_sequence)]
  summary <- stats::aggregate(
    cbind(
      weighted_N = x$N * x$fraction,
      mass = x$fraction,
      high_weighted_N = ifelse(x$N >= 66, x$N * x$fraction, 0),
      high_mass = ifelse(x$N >= 66, x$fraction, 0)
    ) ~ passage_sequence,
    data = x,
    FUN = sum,
    na.rm = TRUE
  )
  summary$mean_ploidy <- summary$weighted_N / summary$mass / 22
  summary$high_component_mean_ploidy <-
    summary$high_weighted_N / summary$high_mass / 22
  summary$high_component_fraction <- summary$high_mass / summary$mass
  summary$seed <- seed_bundle$seed_id
  summary$objective_rank <- match(seed_bundle$seed_number, s65_seeds())
  x$ploidy <- x$N / 22
  x$seed <- seed_bundle$seed_id
  x$objective_rank <- match(seed_bundle$seed_number, s65_seeds())
  list(distribution = x, summary = summary)
}

s65_prepare_fixed_model <- function(params, config, max_o2 = 20.5) {
  run_params <- prepare_run_params(
    param_values = params,
    simulation = "invitro",
    cfg = config,
    fixed_o2 = max_o2
  )
  sim_cfg <- prepare_sim_cfg(
    config,
    argv = list(),
    fixed_o2 = max_o2,
    run_params = run_params
  )
  run_params$O2_growth <- isTRUE(sim_cfg$O2_growth)
  run_params$ploidy_O2_death <- sim_cfg$ploidy_O2_death
  list(run_params = run_params, config = sim_cfg)
}

s65_dominant_bundle <- function(M, ngrid) {
  eig <- tryCatch(eigen(as.matrix(M), only.values = FALSE), error = function(e) NULL)
  if (is.null(eig)) stop("Fixed-environment eigen decomposition failed.")
  order_index <- order(Re(eig$values), decreasing = TRUE)
  eig$values <- eig$values[order_index]
  eig$vectors <- eig$vectors[, order_index, drop = FALSE]
  raw <- Re(eig$vectors[, 1L])
  if (sum(raw, na.rm = TRUE) < 0) {
    raw <- -raw
    eig$vectors[, 1L] <- -eig$vectors[, 1L]
  }
  nonnegative <- all(raw >= -1e-8)
  dominant <- s65_normalize_state(raw)
  if (any(!is.finite(dominant))) {
    stop("Dominant eigenvector could not be normalized.")
  }
  lambda1 <- Re(eig$values[[1L]])
  lambda2 <- Re(eig$values[[2L]])
  gap <- lambda1 - lambda2
  if (!is.finite(gap) || gap <= 0) {
    stop("Fixed-environment operator does not have a positive dominant spectral gap.")
  }
  residual <- max(abs(
    as.numeric(M %*% dominant) - lambda1 * dominant
  ))
  list(
    eig = eig,
    dominant = dominant,
    lambda1 = lambda1,
    lambda2 = lambda2,
    gap = gap,
    tau_days = 1 / gap,
    eigenpair_residual = residual,
    eigenvector_nonnegative = nonnegative,
    minimum_oriented_component = min(raw)
  )
}

s65_eigen_state <- function(dominant_bundle, initial_state, day) {
  eig <- dominant_bundle$eig
  coefficients <- tryCatch(
    solve(eig$vectors, as.numeric(initial_state)),
    error = function(e) NULL
  )
  if (is.null(coefficients)) return(rep(NA_real_, length(initial_state)))
  weights <- exp((eig$values - dominant_bundle$lambda1) * day) * coefficients
  s65_normalize_state(eig$vectors %*% weights)
}

s65_expm_state <- function(
    M, dominant_bundle, initial_state, day,
    method = "Higham08.b"
) {
  if (abs(day) <= 1e-15) return(s65_normalize_state(initial_state))
  shifted <- as.matrix(M)
  diag(shifted) <- diag(shifted) - dominant_bundle$lambda1
  propagated <- tryCatch(
    expm::expm(shifted * day, method = method) %*% as.numeric(initial_state),
    error = function(e) NULL
  )
  if (is.null(propagated)) return(rep(NA_real_, length(initial_state)))
  s65_normalize_state(propagated)
}

s65_fixed_operator_bundle <- function(params, config, oxygen, fixed_p) {
  prepared <- s65_prepare_fixed_model(params, config, max_o2 = max(20.5, oxygen))
  forced <- response_force_effective_p_misseg(
    prepared$run_params, fixed_p
  )
  fm <- fixo2_fixed_matrix(
    model_env = globalenv(), cfg = prepared$config,
    run_params = forced, O2 = oxygen
  )
  dominant <- s65_dominant_bundle(fm$M, fm$ngrid)
  dominant_mean <- s65_mean_ploidy(
    dominant$dominant, fm$ngrid, prepared$config$N_UNIT
  )
  actual_p <- fixo2_population_average_p_misseg(
    state = dominant$dominant, ngrid = fm$ngrid, O2 = oxygen,
    run_params = forced, model_env = globalenv()
  )
  list(
    prepared = prepared,
    forced = forced,
    fixed_matrix = fm,
    dominant = dominant,
    dominant_mean_ploidy = dominant_mean,
    actual_p_miss_eff = actual_p
  )
}

s65_transition_matrix <- function(M, dominant_bundle, day, method) {
  if (abs(day) <= 1e-15) return(diag(nrow(M)))
  shifted <- as.matrix(M)
  diag(shifted) <- diag(shifted) - dominant_bundle$lambda1
  transition <- tryCatch(
    expm::expm(shifted * day, method = method),
    error = function(e) NULL
  )
  if (is.null(transition) || any(!is.finite(transition))) {
    stop("Non-finite matrix exponential from method ", method, ".")
  }
  transition
}

s65_operator_diagnostic <- function(
    endpoint_scope, endpoint_label, family, endpoint_multiplicity,
    params, config, oxygen, fixed_p
) {
  bundle <- s65_fixed_operator_bundle(params, config, oxygen, fixed_p)
  dominant <- bundle$dominant
  data.frame(
    endpoint_scope = endpoint_scope,
    endpoint_label = endpoint_label,
    family = family,
    endpoint_multiplicity = as.integer(endpoint_multiplicity),
    oxygen_pct = as.numeric(oxygen),
    fixed_p_miss_eff = as.numeric(fixed_p),
    actual_p_miss_eff = as.numeric(bundle$actual_p_miss_eff),
    requested_cin_error = abs(bundle$actual_p_miss_eff - fixed_p),
    spectral_gap = dominant$gap,
    relaxation_time_days = dominant$tau_days,
    dominant_mean_ploidy = bundle$dominant_mean_ploidy,
    dominant_eigenpair_residual = dominant$eigenpair_residual,
    eigenvector_nonnegative = dominant$eigenvector_nonnegative,
    minimum_oriented_eigenvector_component =
      dominant$minimum_oriented_component,
    stringsAsFactors = FALSE
  )
}

s65_operator_trajectory <- function(
    endpoint_scope, endpoint_label, family, endpoint_multiplicity,
    params, config, initial_states, oxygen, fixed_p,
    u_values = s65_u_values(), solver_u_values = s65_solver_u_values(),
    run_solver_audit = TRUE
) {
  bundle <- s65_fixed_operator_bundle(params, config, oxygen, fixed_p)
  prepared <- bundle$prepared
  fm <- bundle$fixed_matrix
  dominant <- bundle$dominant
  dominant_mean <- bundle$dominant_mean_ploidy
  actual_p <- bundle$actual_p_miss_eff
  days <- as.numeric(u_values) / dominant$gap
  primary_transition <- lapply(days, function(day) {
    s65_transition_matrix(fm$M, dominant, day, "Higham08.b")
  })
  solver_requested <- vapply(as.numeric(u_values), function(u) {
    isTRUE(run_solver_audit) && any(abs(u - solver_u_values) <= 1e-12)
  }, logical(1L))
  secondary_transition <- lapply(seq_along(days), function(index) {
    if (solver_requested[[index]]) {
      s65_transition_matrix(fm$M, dominant, days[[index]], "Ward77")
    } else {
      NULL
    }
  })
  rows <- list()
  states <- list()
  for (initial_label in names(initial_states)) {
    initial <- s65_normalize_state(initial_states[[initial_label]])
    if (length(initial) != length(fm$ngrid) || any(!is.finite(initial))) {
      stop("Initial state does not match fixed-operator grid for ", endpoint_label)
    }
    for (u_index in seq_along(u_values)) {
      u <- as.numeric(u_values[[u_index]])
      day <- days[[u_index]]
      primary_state <- s65_normalize_state(
        primary_transition[[u_index]] %*% initial
      )
      if (any(!is.finite(primary_state))) {
        stop(
          "Non-finite Higham matrix-exponential propagation for ",
          endpoint_label, ", O2=", oxygen, ", p=", fixed_p,
          ", initial=", initial_label, ", gap-scaled time=", u
        )
      }
      secondary_state <- if (solver_requested[[u_index]]) {
        s65_normalize_state(
          secondary_transition[[u_index]] %*% initial
        )
      } else {
        rep(NA_real_, length(initial))
      }
      if (solver_requested[[u_index]] &&
          any(!is.finite(secondary_state))) {
        stop("Non-finite Ward matrix-exponential propagation for ", endpoint_label)
      }
      mean_ploidy <- s65_mean_ploidy(
        primary_state, fm$ngrid, prepared$config$N_UNIT
      )
      solver_mean <- if (solver_requested[[u_index]]) {
        s65_mean_ploidy(
          secondary_state, fm$ngrid, prepared$config$N_UNIT
        )
      } else {
        NA_real_
      }
      row_id <- length(rows) + 1L
      rows[[row_id]] <- data.frame(
        endpoint_scope = endpoint_scope,
        endpoint_label = endpoint_label,
        family = family,
        endpoint_multiplicity = as.integer(endpoint_multiplicity),
        initial_state = initial_label,
        oxygen_pct = as.numeric(oxygen),
        fixed_p_miss_eff = as.numeric(fixed_p),
        actual_p_miss_eff = as.numeric(actual_p),
        requested_cin_error = abs(actual_p - fixed_p),
        spectral_gap = dominant$gap,
        relaxation_time_days = dominant$tau_days,
        gap_scaled_time = as.numeric(u),
        day = day,
        mean_ploidy = mean_ploidy,
        dominant_mean_ploidy = dominant_mean,
        mean_ploidy_error = abs(mean_ploidy - dominant_mean),
        total_variation_to_dominant = s65_total_variation(
          primary_state, dominant$dominant
        ),
        normalized_mass_error = abs(sum(primary_state) - 1),
        solver_audited = solver_requested[[u_index]],
        solver_mean_ploidy_error = if (solver_requested[[u_index]]) {
          abs(mean_ploidy - solver_mean)
        } else {
          NA_real_
        },
        solver_total_variation = if (solver_requested[[u_index]]) {
          s65_total_variation(primary_state, secondary_state)
        } else {
          NA_real_
        },
        dominant_eigenpair_residual = dominant$eigenpair_residual,
        eigenvector_nonnegative = dominant$eigenvector_nonnegative,
        minimum_oriented_eigenvector_component =
          dominant$minimum_oriented_component,
        stringsAsFactors = FALSE
      )
      states[[row_id]] <- list(
        ngrid = fm$ngrid,
        primary_state = primary_state,
        secondary_state = secondary_state,
        dominant_state = dominant$dominant,
        row = rows[[row_id]]
      )
    }
  }
  list(
    trajectory = do.call(rbind, rows),
    states = states,
    ngrid = fm$ngrid,
    dominant = dominant$dominant,
    dominant_mean_ploidy = dominant_mean,
    spectral_gap = dominant$gap
  )
}

s65_separate_audit <- function(seed_bundles, n_core = 1L) {
  tasks <- expand.grid(
    seed_number = vapply(seed_bundles, `[[`, integer(1L), "seed_number"),
    oxygen_pct = s65_separate_o2_values(),
    fixed_p_miss_eff = s65_p_values(),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  compute <- function(index) {
    task <- tasks[index, , drop = FALSE]
    seed <- seed_bundles[[match(
      task$seed_number[[1L]],
      vapply(seed_bundles, `[[`, integer(1L), "seed_number")
    )]]
    full_curve <- abs(task$fixed_p_miss_eff[[1L]] - 0.01) <= 1e-12 &&
      task$oxygen_pct[[1L]] %in% c(0, 5)
    s65_operator_trajectory(
      endpoint_scope = "Separate in vitro",
      endpoint_label = seed$seed_id,
      family = seed$seed_id,
      endpoint_multiplicity = 1L,
      params = seed$params,
      config = seed$config,
      initial_states = seed$initial_states,
      oxygen = task$oxygen_pct[[1L]],
      fixed_p = task$fixed_p_miss_eff[[1L]],
      u_values = if (full_curve) s65_u_values() else s65_terminal_u(),
      solver_u_values = if (full_curve) {
        s65_solver_u_values()
      } else {
        s65_terminal_u()
      },
      run_solver_audit = TRUE
    )
  }
  result <- f6r_resilient_lapply(seq_len(nrow(tasks)), compute, n_core = n_core)
  list(
    trajectory = do.call(rbind, lapply(result, `[[`, "trajectory")),
    operator_diagnostic = unique(do.call(rbind, lapply(
      result, function(x) x$trajectory[, c(
        "endpoint_scope", "endpoint_label", "family",
        "endpoint_multiplicity", "oxygen_pct", "fixed_p_miss_eff",
        "actual_p_miss_eff", "requested_cin_error", "spectral_gap",
        "relaxation_time_days", "dominant_mean_ploidy",
        "dominant_eigenpair_residual", "eigenvector_nonnegative",
        "minimum_oriented_eigenvector_component"
      )]
    ))),
    result = result,
    tasks = tasks
  )
}

s65_joint_endpoint_bundle <- function(paths) {
  displayed <- f6r_read_tsv(file.path(
    paths$base$figure6, "figure6a_displayed_pairs.tsv"
  ))
  manifest <- f6r_read_tsv(file.path(
    paths$base$figure6, "joint_invitro_q20_endpoint_manifest.tsv"
  ))
  parameters <- f6r_read_tsv(file.path(
    paths$base$figure6, "joint_acceptable_seed_invitro_parameters.tsv"
  ))
  s65_require_columns(
    displayed, c("display_label", "pair_id"), "Figure 6A displayed-pair table"
  )
  s65_require_columns(
    manifest,
    c(
      "display_label", "pair_id", "parameter_endpoint_group",
      "representative_seed_number", "representative_objective_rank",
      "endpoint_multiplicity_q10"
    ),
    "joint in-vitro endpoint manifest"
  )
  manifest <- manifest[
    manifest$pair_id %in% displayed$pair_id &
      manifest$endpoint_multiplicity_q10 > 0,
    , drop = FALSE
  ]
  manifest$family <- displayed$display_label[
    match(manifest$pair_id, displayed$pair_id)
  ]
  represented <- tapply(
    manifest$endpoint_multiplicity_q10, manifest$family, sum
  )
  if (nrow(displayed) != 3L || anyNA(manifest$family) ||
      !identical(sort(unique(manifest$family)), c("C01", "C02", "C03")) ||
      any(represented != 50L)) {
    stop("Expected three displayed C families with 50 q10 endpoints each.")
  }
  selected <- f6r_selected_results(paths$base)
  contexts <- lapply(displayed$pair_id, f6r_pair_model_context,
    selected = selected, paths = paths$base
  )
  names(contexts) <- displayed$pair_id
  list(
    displayed = displayed,
    manifest = manifest,
    parameters = parameters,
    contexts = contexts,
    paths = c(
      displayed = normalizePath(file.path(
        paths$base$figure6, "figure6a_displayed_pairs.tsv"
      ), mustWork = TRUE),
      manifest = normalizePath(file.path(
        paths$base$figure6, "joint_invitro_q20_endpoint_manifest.tsv"
      ), mustWork = TRUE),
      parameters = normalizePath(file.path(
        paths$base$figure6,
        "joint_acceptable_seed_invitro_parameters.tsv"
      ), mustWork = TRUE)
    )
  )
}

s65_joint_audit <- function(
    paths, initial_states, n_core = 1L
) {
  bundle <- s65_joint_endpoint_bundle(paths)
  manifest <- bundle$manifest
  tasks <- merge(
    manifest[, c(
      "family", "pair_id", "parameter_endpoint_group",
      "representative_seed_number", "representative_objective_rank",
      "endpoint_multiplicity_q10"
    )],
    expand.grid(
      oxygen_pct = s65_o2_values(),
      fixed_p_miss_eff = s65_p_values(),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  )
  tasks <- tasks[order(
    tasks$family, tasks$representative_objective_rank,
    tasks$oxygen_pct, tasks$fixed_p_miss_eff
  ), , drop = FALSE]
  compute <- function(index) {
    task <- tasks[index, , drop = FALSE]
    z <- bundle$parameters[
      bundle$parameters$pair_id == task$pair_id[[1L]] &
        bundle$parameters$seed_number ==
          task$representative_seed_number[[1L]],
      , drop = FALSE
    ]
    shared <- s65_shared_parameters()
    if (nrow(z) != length(shared) || !setequal(z$parameter, shared)) {
      stop(
        "Incomplete in-vitro parameters for ", task$family[[1L]],
        " endpoint rank ", task$representative_objective_rank[[1L]]
      )
    }
    params <- stats::setNames(as.numeric(z$value), z$parameter)
    context <- bundle$contexts[[task$pair_id[[1L]]]]
    diagnostic <- s65_operator_diagnostic(
      endpoint_scope = "Figure 6A joint ensemble",
      endpoint_label = paste0(
        task$family[[1L]], " endpoint ",
        task$representative_objective_rank[[1L]]
      ),
      family = task$family[[1L]],
      endpoint_multiplicity = task$endpoint_multiplicity_q10[[1L]],
      params = params,
      config = context$config,
      oxygen = task$oxygen_pct[[1L]],
      fixed_p = task$fixed_p_miss_eff[[1L]]
    )
    trajectory <- if (task$representative_objective_rank[[1L]] == 1L) {
      s65_operator_trajectory(
        endpoint_scope = "Figure 6A joint ensemble",
        endpoint_label = paste0(
          task$family[[1L]], " endpoint ",
          task$representative_objective_rank[[1L]]
        ),
        family = task$family[[1L]],
        endpoint_multiplicity = task$endpoint_multiplicity_q10[[1L]],
        params = params,
        config = context$config,
        initial_states = initial_states,
        oxygen = task$oxygen_pct[[1L]],
        fixed_p = task$fixed_p_miss_eff[[1L]],
        u_values = s65_terminal_u(),
        solver_u_values = s65_terminal_u(),
        run_solver_audit = TRUE
      )$trajectory
    } else {
      NULL
    }
    list(diagnostic = diagnostic, trajectory = trajectory)
  }
  result <- f6r_resilient_lapply(seq_len(nrow(tasks)), compute, n_core = n_core)
  list(
    trajectory = do.call(rbind, Filter(
      Negate(is.null), lapply(result, `[[`, "trajectory")
    )),
    operator_diagnostic = do.call(rbind, lapply(
      result, `[[`, "diagnostic"
    )),
    tasks = tasks,
    bundle = bundle
  )
}

s65_snapshot_data <- function(separate_result, seed_bundle) {
  target_index <- which(
    separate_result$tasks$seed_number == 10L &
      abs(separate_result$tasks$oxygen_pct - 0) <= 1e-12 &
      abs(separate_result$tasks$fixed_p_miss_eff - 0.01) <= 1e-12
  )
  if (length(target_index) != 1L) {
    stop("Could not resolve the seed10 fixed 0% O2, p=0.01 snapshot operator.")
  }
  operator <- separate_result$result[[target_index]]
  experiment_duration <- do.call(rbind, lapply(
    c("2N", "4N"),
    function(cohort) s65_experiment_duration(seed_bundle$lineage, cohort)
  ))
  state_rows <- list()
  row_id <- 0L
  for (initial_label in names(seed_bundle$initial_states)) {
    cohort <- sub(" start$", "", initial_label)
    initial <- seed_bundle$initial_states[[initial_label]]
    gap <- operator$spectral_gap
    duration <- experiment_duration$experimental_duration_days[
      experiment_duration$cohort == cohort
    ]
    times <- c(
      `Initial` = 0,
      `Experimental duration` = duration,
      `5 relaxation times` = 5 / gap,
      `20 relaxation times` = 20 / gap,
      `40 relaxation times` = 40 / gap
    )
    target_state <- operator$dominant
    for (label in names(times)) {
      row_id <- row_id + 1L
      state <- if (identical(label, "Initial")) {
        s65_normalize_state(initial)
      } else {
        prepared <- s65_prepare_fixed_model(
          seed_bundle$params, seed_bundle$config, 20.5
        )
        fm <- fixo2_fixed_matrix(
          globalenv(), prepared$config,
          response_force_effective_p_misseg(prepared$run_params, 0.01), 0
        )
        dominant <- s65_dominant_bundle(fm$M, fm$ngrid)
        s65_expm_state(
          fm$M, dominant, initial, times[[label]], method = "Higham08.b"
        )
      }
      state_rows[[row_id]] <- data.frame(
        initial_state = initial_label,
        timepoint = label,
        day = as.numeric(times[[label]]),
        N = operator$ngrid,
        ploidy = operator$ngrid / 22,
        fraction = state,
        stringsAsFactors = FALSE
      )
    }
    row_id <- row_id + 1L
    state_rows[[row_id]] <- data.frame(
      initial_state = initial_label,
      timepoint = "Dominant eigenvector",
      day = Inf,
      N = operator$ngrid,
      ploidy = operator$ngrid / 22,
      fraction = target_state,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, state_rows)
  out$timepoint <- factor(
    out$timepoint,
    levels = c(
      "Initial", "Experimental duration", "5 relaxation times",
      "20 relaxation times", "40 relaxation times",
      "Dominant eigenvector"
    )
  )
  out
}

s65_group_maxima <- function(trajectory, operator_diagnostic) {
  terminal <- trajectory[
    abs(trajectory$gap_scaled_time - s65_terminal_u()) <= 1e-12,
  ]
  solver <- trajectory[trajectory$solver_audited, ]
  groups <- unique(trajectory[, c("endpoint_scope", "family")])
  rows <- lapply(seq_len(nrow(groups)), function(index) {
    group <- groups[index, , drop = FALSE]
    t <- terminal[
      terminal$endpoint_scope == group$endpoint_scope &
        terminal$family == group$family,
      , drop = FALSE
    ]
    s <- solver[
      solver$endpoint_scope == group$endpoint_scope &
        solver$family == group$family,
      , drop = FALSE
    ]
    d <- operator_diagnostic[
      operator_diagnostic$endpoint_scope == group$endpoint_scope &
        operator_diagnostic$family == group$family,
      , drop = FALSE
    ]
    represented <- unique(d[, c(
      "endpoint_label", "endpoint_multiplicity"
    )])
    data.frame(
      endpoint_scope = group$endpoint_scope,
      family = group$family,
      represented_endpoint_count = sum(represented$endpoint_multiplicity),
      unique_parameter_endpoint_count = nrow(represented),
      audited_fixed_operator_count = nrow(unique(d[, c(
        "endpoint_label", "oxygen_pct", "fixed_p_miss_eff"
      )])),
      finite_time_endpoint_count = length(unique(t$endpoint_label)),
      operator_initial_condition_count = nrow(t),
      maximum_requested_cin_error = max(d$requested_cin_error),
      maximum_dominant_eigenpair_residual =
        max(d$dominant_eigenpair_residual),
      maximum_terminal_mean_ploidy_error = max(t$mean_ploidy_error),
      maximum_terminal_total_variation =
        max(t$total_variation_to_dominant),
      maximum_normalized_mass_error = max(t$normalized_mass_error),
      maximum_solver_mean_ploidy_error = if (nrow(s)) {
        max(s$solver_mean_ploidy_error, na.rm = TRUE)
      } else {
        NA_real_
      },
      maximum_solver_total_variation = if (nrow(s)) {
        max(s$solver_total_variation, na.rm = TRUE)
      } else {
        NA_real_
      },
      all_dominant_eigenvectors_nonnegative = all(d$eigenvector_nonnegative),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

s65_initial_state_terminal_difference <- function(trajectory) {
  terminal <- trajectory[
    abs(trajectory$gap_scaled_time - s65_terminal_u()) <= 1e-12,
  ]
  keys <- unique(terminal[, c(
    "endpoint_scope", "endpoint_label", "family", "oxygen_pct",
    "fixed_p_miss_eff", "endpoint_multiplicity"
  )])
  rows <- lapply(seq_len(nrow(keys)), function(index) {
    key <- keys[index, , drop = FALSE]
    z <- terminal[
      terminal$endpoint_scope == key$endpoint_scope &
        terminal$endpoint_label == key$endpoint_label &
        terminal$family == key$family &
        abs(terminal$oxygen_pct - key$oxygen_pct) <= 1e-12 &
        abs(terminal$fixed_p_miss_eff - key$fixed_p_miss_eff) <= 1e-12,
      , drop = FALSE
    ]
    if (nrow(z) != 2L || !setequal(z$initial_state, c("2N start", "4N start"))) {
      stop("Expected two terminal initial-state rows for ", key$endpoint_label)
    }
    data.frame(
      endpoint_scope = key$endpoint_scope,
      endpoint_label = key$endpoint_label,
      family = key$family,
      endpoint_multiplicity = key$endpoint_multiplicity,
      oxygen_pct = key$oxygen_pct,
      fixed_p_miss_eff = key$fixed_p_miss_eff,
      terminal_mean_ploidy_difference = diff(range(z$mean_ploidy)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

s65_validation_table <- function(
    separate_bundles, experimental_summary, separate_trajectory,
    joint_trajectory, separate_diagnostic, joint_diagnostic,
    initial_difference
) {
  thresholds <- s65_thresholds()
  combined <- rbind(separate_trajectory, joint_trajectory)
  diagnostic <- rbind(separate_diagnostic, joint_diagnostic)
  terminal <- combined[
    abs(combined$gap_scaled_time - s65_terminal_u()) <= 1e-12,
  ]
  solver <- combined[combined$solver_audited, ]
  exp_delta <- do.call(rbind, lapply(split(
    experimental_summary, experimental_summary$seed
  ), function(z) {
    z <- z[order(z$passage_sequence), ]
    data.frame(
      seed = z$seed[[1L]],
      start_mean_ploidy = z$mean_ploidy[[1L]],
      end_mean_ploidy = z$mean_ploidy[[nrow(z)]],
      change_mean_ploidy = z$mean_ploidy[[nrow(z)]] - z$mean_ploidy[[1L]],
      start_high_component_mean_ploidy =
        z$high_component_mean_ploidy[[1L]],
      end_high_component_mean_ploidy =
        z$high_component_mean_ploidy[[nrow(z)]],
      change_high_component_mean_ploidy =
        z$high_component_mean_ploidy[[nrow(z)]] -
          z$high_component_mean_ploidy[[1L]],
      stringsAsFactors = FALSE
    )
  }))
  joint_unique_endpoints <- unique(joint_diagnostic[, c(
    "family", "endpoint_label", "endpoint_multiplicity"
  )])
  joint_represented_counts <- tapply(
    joint_unique_endpoints$endpoint_multiplicity,
    joint_unique_endpoints$family,
    sum
  )
  joint_represented_counts <- joint_represented_counts[
    c("C01", "C02", "C03")
  ]
  checks <- data.frame(
    check = c(
      "predeclared_separate_seeds_present",
      "experimental_4N_mean_declines_in_both_seeds",
      "experimental_high_component_mean_declines_in_both_seeds",
      "all_requested_cin_values_reproduced",
      "all_dominant_eigenpair_residuals_pass",
      "all_solver_mean_ploidy_checks_pass",
      "all_solver_distribution_checks_pass",
      "all_terminal_mean_ploidy_checks_pass",
      "all_terminal_distribution_checks_pass",
      "all_initial_state_terminal_checks_pass",
      "all_normalized_state_mass_checks_pass",
      "all_dominant_eigenvectors_nonnegative",
      "three_joint_families_represent_50_endpoints_each"
    ),
    observed = c(
      paste(sort(vapply(separate_bundles, `[[`, character(1L), "seed_id")), collapse = ","),
      paste(signif(exp_delta$change_mean_ploidy, 8), collapse = ","),
      paste(signif(exp_delta$change_high_component_mean_ploidy, 8), collapse = ","),
      format(max(diagnostic$requested_cin_error), scientific = TRUE),
      format(
        max(diagnostic$dominant_eigenpair_residual), scientific = TRUE
      ),
      format(max(solver$solver_mean_ploidy_error, na.rm = TRUE), scientific = TRUE),
      format(max(solver$solver_total_variation, na.rm = TRUE), scientific = TRUE),
      format(max(terminal$mean_ploidy_error), scientific = TRUE),
      format(max(terminal$total_variation_to_dominant), scientific = TRUE),
      format(max(initial_difference$terminal_mean_ploidy_difference), scientific = TRUE),
      format(max(terminal$normalized_mass_error), scientific = TRUE),
      as.character(all(terminal$eigenvector_nonnegative)),
      paste(joint_represented_counts, collapse = ",")
    ),
    expected = c(
      "seed10,seed132", "both < 0", "both < 0",
      paste0("<=", thresholds[["requested_cin_error"]]),
      paste0("<=", thresholds[["dominant_eigenpair_residual"]]),
      paste0("<=", thresholds[["solver_mean_ploidy_error"]]),
      paste0("<=", thresholds[["solver_total_variation"]]),
      paste0("<=", thresholds[["terminal_mean_ploidy_error"]]),
      paste0("<=", thresholds[["terminal_total_variation"]]),
      paste0("<=", thresholds[["initial_state_terminal_difference"]]),
      paste0("<=", thresholds[["normalized_mass_error"]]),
      "TRUE", "50,50,50"
    ),
    passed = c(
      identical(
        sort(vapply(separate_bundles, `[[`, character(1L), "seed_id")),
        c("seed10", "seed132")
      ),
      all(exp_delta$change_mean_ploidy < 0),
      all(exp_delta$change_high_component_mean_ploidy < 0),
      max(diagnostic$requested_cin_error) <=
        thresholds[["requested_cin_error"]],
      max(diagnostic$dominant_eigenpair_residual) <=
        thresholds[["dominant_eigenpair_residual"]],
      max(solver$solver_mean_ploidy_error, na.rm = TRUE) <=
        thresholds[["solver_mean_ploidy_error"]],
      max(solver$solver_total_variation, na.rm = TRUE) <=
        thresholds[["solver_total_variation"]],
      max(terminal$mean_ploidy_error) <=
        thresholds[["terminal_mean_ploidy_error"]],
      max(terminal$total_variation_to_dominant) <=
        thresholds[["terminal_total_variation"]],
      max(initial_difference$terminal_mean_ploidy_difference) <=
        thresholds[["initial_state_terminal_difference"]],
      max(terminal$normalized_mass_error) <=
        thresholds[["normalized_mass_error"]],
      all(diagnostic$eigenvector_nonnegative),
      length(joint_represented_counts) == 3L &&
        all(joint_represented_counts == 50L)
    ),
    stringsAsFactors = FALSE
  )
  list(checks = checks, experimental_delta = exp_delta)
}

s65_write_provenance <- function(paths, separate_bundles, joint_bundle) {
  model_files <- f6x_model_source_files(paths$base)
  source_files <- c(
    unlist(lapply(separate_bundles, `[[`, "files"), use.names = FALSE),
    unname(joint_bundle$paths),
    model_files
  )
  source_files <- unique(normalizePath(source_files, mustWork = TRUE))
  role <- ifelse(
    source_files %in% normalizePath(model_files, mustWork = TRUE),
    "current-branch fixed-environment model source",
    ifelse(
      grepl("seed10|seed132", source_files),
      "predeclared separate in-vitro fit input",
      "frozen Figure 6A joint-ensemble input"
    )
  )
  provenance <- data.frame(
    role = role,
    source_file = source_files,
    md5 = vapply(source_files, f6r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(
    provenance,
    file.path(paths$data, "supp_figure6-5_source_provenance.tsv")
  )
}

s65_data <- function(
    workspace_root = f6r_find_workspace_root(),
    invitro_result_root,
    n_core = 4L
) {
  f6r_require_packages(c("Matrix", "expm", "Rcpp", "data.table"))
  paths <- s65_paths(workspace_root)
  dir.create(paths$data, recursive = TRUE, showWarnings = FALSE)
  f6r_load_response_engine(paths$base)
  if (as.integer(n_core) != 1L) {
    message(
      "Supplementary Figure 6-5 uses sequential linear algebra to avoid ",
      "forking an OpenMP-backed matrix exponential."
    )
  }
  n_core <- 1L

  separate_bundles <- lapply(
    s65_seeds(), s65_load_separate_seed,
    result_root = invitro_result_root
  )
  experimental <- lapply(separate_bundles, s65_experimental_distribution)
  experimental_distribution <- do.call(rbind, lapply(
    experimental, `[[`, "distribution"
  ))
  experimental_summary <- do.call(rbind, lapply(
    experimental, `[[`, "summary"
  ))
  durations <- do.call(rbind, lapply(separate_bundles, function(seed) {
    x <- do.call(rbind, lapply(
      c("2N", "4N"),
      function(cohort) s65_experiment_duration(seed$lineage, cohort)
    ))
    x$seed <- seed$seed_id
    x
  }))

  message("Supplementary Figure 6-5: separate-fit fixed-operator audit")
  separate <- s65_separate_audit(
    separate_bundles, n_core = as.integer(n_core)
  )
  message("Supplementary Figure 6-5: Figure 6A joint-ensemble audit")
  joint <- s65_joint_audit(
    paths,
    initial_states = separate_bundles[[1L]]$initial_states,
    n_core = as.integer(n_core)
  )
  combined <- rbind(separate$trajectory, joint$trajectory)
  combined_diagnostic <- rbind(
    separate$operator_diagnostic, joint$operator_diagnostic
  )
  initial_difference <- s65_initial_state_terminal_difference(combined)
  group_maxima <- s65_group_maxima(combined, combined_diagnostic)
  snapshots <- s65_snapshot_data(separate, separate_bundles[[1L]])
  validation <- s65_validation_table(
    separate_bundles, experimental_summary,
    separate$trajectory, joint$trajectory,
    separate$operator_diagnostic, joint$operator_diagnostic,
    initial_difference
  )

  output_paths <- c(
    experimental_distribution = f6r_write_tsv(
      experimental_distribution,
      file.path(paths$data, "experimental_4N_deprived_distribution.tsv")
    ),
    experimental_summary = f6r_write_tsv(
      experimental_summary,
      file.path(paths$data, "experimental_4N_deprived_summary.tsv")
    ),
    experimental_duration = f6r_write_tsv(
      durations,
      file.path(paths$data, "experimental_duration_summary.tsv")
    ),
    separate_trajectory = f6r_write_tsv(
      separate$trajectory,
      file.path(paths$data, "separate_seed_fixed_operator_trajectory.tsv")
    ),
    joint_trajectory = f6r_write_tsv(
      joint$trajectory,
      file.path(paths$data, "joint_ensemble_fixed_operator_terminal.tsv")
    ),
    separate_operator_diagnostic = f6r_write_tsv(
      separate$operator_diagnostic,
      file.path(paths$data, "separate_seed_fixed_operator_diagnostic.tsv")
    ),
    joint_operator_diagnostic = f6r_write_tsv(
      joint$operator_diagnostic,
      file.path(paths$data, "joint_ensemble_fixed_operator_diagnostic.tsv")
    ),
    initial_difference = f6r_write_tsv(
      initial_difference,
      file.path(paths$data, "initial_state_terminal_difference.tsv")
    ),
    group_maxima = f6r_write_tsv(
      group_maxima,
      file.path(paths$data, "numerical_consistency_group_maxima.tsv")
    ),
    snapshots = f6r_write_tsv(
      snapshots,
      file.path(paths$data, "seed10_distribution_snapshots.tsv")
    ),
    experimental_delta = f6r_write_tsv(
      validation$experimental_delta,
      file.path(paths$data, "experimental_4N_change_summary.tsv")
    ),
    validation = f6r_write_tsv(
      validation$checks,
      file.path(paths$data, "supp_figure6-5_validation.tsv")
    )
  )
  provenance <- s65_write_provenance(
    paths, separate_bundles, joint$bundle
  )
  output_paths <- c(output_paths, provenance = provenance)
  manifest <- data.frame(
    artifact = names(output_paths),
    path = unname(output_paths),
    md5 = vapply(unname(output_paths), f6r_md5, character(1L)),
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(
    manifest,
    file.path(paths$data, "supp_figure6-5_data_manifest.tsv")
  )
  if (!all(validation$checks$passed)) {
    failed <- validation$checks$check[!validation$checks$passed]
    stop(
      "Supplementary Figure 6-5 numerical validation failed: ",
      paste(failed, collapse = ", ")
    )
  }
  invisible(list(
    paths = paths,
    experimental_distribution = experimental_distribution,
    experimental_summary = experimental_summary,
    durations = durations,
    separate_trajectory = separate$trajectory,
    joint_trajectory = joint$trajectory,
    separate_operator_diagnostic = separate$operator_diagnostic,
    joint_operator_diagnostic = joint$operator_diagnostic,
    initial_difference = initial_difference,
    group_maxima = group_maxima,
    snapshots = snapshots,
    validation = validation$checks,
    output_paths = output_paths
  ))
}

s65_theme <- function(base_size = 9.2) {
  ggplot2::theme_bw(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = base_size + 0.8, hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        color = "#4B5055", size = base_size - 0.5, hjust = 0
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "#202020"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.3),
      strip.background = ggplot2::element_rect(
        fill = "#F1F2F4", color = "#B8BDC3", linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "#E4E6E8", linewidth = 0.25
      ),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

s65_panel_a <- function(paths) {
  distribution <- f6r_read_tsv(file.path(
    paths$data, "experimental_4N_deprived_distribution.tsv"
  ))
  summary <- f6r_read_tsv(file.path(
    paths$data, "experimental_4N_deprived_summary.tsv"
  ))
  distribution$seed <- factor(
    distribution$seed, levels = c("seed10", "seed132"),
    labels = c("seed10, objective rank 1", "seed132, objective rank 2")
  )
  summary$seed <- factor(
    summary$seed, levels = c("seed10", "seed132"),
    labels = c("seed10, objective rank 1", "seed132, objective rank 2")
  )
  ggplot2::ggplot(
    distribution,
    ggplot2::aes(passage_sequence, ploidy, fill = fraction)
  ) +
    ggplot2::geom_tile(width = 0.98, height = 1 / 22) +
    ggplot2::geom_line(
      data = summary,
      ggplot2::aes(passage_sequence, mean_ploidy),
      inherit.aes = FALSE,
      color = "#D55E00", linewidth = 0.72
    ) +
    ggplot2::geom_line(
      data = summary,
      ggplot2::aes(passage_sequence, high_component_mean_ploidy),
      inherit.aes = FALSE,
      color = "white", linewidth = 0.85, linetype = 2
    ) +
    ggplot2::facet_wrap(~seed, ncol = 2) +
    ggplot2::scale_fill_gradientn(
      colors = c("white", "#56B4E9", "#0072B2", "#F0E442"),
      trans = scales::pseudo_log_trans(sigma = 0.002),
      name = "Fraction"
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(2.6, 4.5), expand = FALSE) +
    ggplot2::labs(
      title = "A  Finite-time redistribution in 4N-start hypoxic cultures",
      subtitle = paste0(
        "Orange, population mean; white dashed, mean within states at or above 3N"
      ),
      x = "Oxygen-deprived passage sequence",
      y = "Modeled ploidy"
    ) +
    s65_theme(9.2) +
    ggplot2::theme(legend.position = "right")
}

s65_panel_b <- function(paths) {
  trajectory <- f6r_read_tsv(file.path(
    paths$data, "separate_seed_fixed_operator_trajectory.tsv"
  ))
  trajectory <- trajectory[
    trajectory$oxygen_pct %in% c(0, 5) &
      abs(trajectory$fixed_p_miss_eff - 0.01) <= 1e-12,
    , drop = FALSE
  ]
  trajectory$seed <- factor(
    trajectory$endpoint_label,
    levels = c("seed10", "seed132"),
    labels = c("seed10, rank 1", "seed132, rank 2")
  )
  trajectory$condition <- factor(
    trajectory$oxygen_pct,
    levels = c(0, 5),
    labels = c("0% O2", "5% O2")
  )
  trajectory$initial_state <- factor(
    trajectory$initial_state,
    levels = c("2N start", "4N start")
  )
  dominant <- unique(trajectory[, c(
    "seed", "condition", "dominant_mean_ploidy"
  )])
  ggplot2::ggplot(
    trajectory,
    ggplot2::aes(
      gap_scaled_time, mean_ploidy,
      color = initial_state, linetype = initial_state
    )
  ) +
    ggplot2::geom_hline(
      data = dominant,
      ggplot2::aes(yintercept = dominant_mean_ploidy),
      color = "#222222", linetype = 3, linewidth = 0.55
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_grid(seed ~ condition) +
    ggplot2::scale_color_manual(
      values = c("2N start" = "#0072B2", "4N start" = "#D55E00")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("2N start" = 1, "4N start" = 2)
    ) +
    ggplot2::labs(
      title = "B  Matched fixed-environment trajectories approach one dominant mean",
      subtitle = paste0(
        "Fixed effective missegregation probability = 0.01; dotted lines are dominant-eigenvector means"
      ),
      x = "Time scaled by the dominant spectral gap",
      y = "Mean ploidy",
      color = "Initial distribution",
      linetype = "Initial distribution"
    ) +
    s65_theme(9.2) +
    ggplot2::theme(legend.position = "top")
}

s65_panel_c <- function(paths) {
  snapshots <- f6r_read_tsv(file.path(
    paths$data, "seed10_distribution_snapshots.tsv"
  ))
  snapshots$initial_state <- factor(
    snapshots$initial_state, levels = c("2N start", "4N start")
  )
  snapshots$timepoint <- factor(
    snapshots$timepoint,
    levels = c(
      "Initial", "Experimental duration", "5 relaxation times",
      "20 relaxation times", "40 relaxation times",
      "Dominant eigenvector"
    ),
    labels = c(
      "Initial", "Experiment length", "5 tau", "20 tau", "40 tau",
      "Dominant"
    )
  )
  ggplot2::ggplot(
    snapshots,
    ggplot2::aes(ploidy, fraction, color = initial_state)
  ) +
    ggplot2::geom_line(linewidth = 0.62) +
    ggplot2::facet_grid(initial_state ~ timepoint) +
    ggplot2::scale_color_manual(
      values = c("2N start" = "#0072B2", "4N start" = "#D55E00")
    ) +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 1e-4),
      breaks = c(0, 0.001, 0.01, 0.1)
    ) +
    ggplot2::coord_cartesian(xlim = c(1, 7)) +
    ggplot2::labs(
      title = "C  Full chromosome-state distributions converge, not only their means",
      subtitle = paste0(
        "seed10 at fixed 0% O2 and effective missegregation probability = 0.01"
      ),
      x = "Modeled ploidy",
      y = "Normalized fraction"
    ) +
    s65_theme(8.5) +
    ggplot2::theme(legend.position = "none")
}

s65_panel_d <- function(paths) {
  maxima <- f6r_read_tsv(file.path(
    paths$data, "numerical_consistency_group_maxima.tsv"
  ))
  thresholds <- s65_thresholds()
  maxima$group <- maxima$family
  maxima$group <- factor(
    maxima$group,
    levels = c("seed10", "seed132", "C01", "C02", "C03")
  )
  metric_map <- data.frame(
    column = c(
      "maximum_terminal_total_variation",
      "maximum_terminal_mean_ploidy_error",
      "maximum_solver_total_variation",
      "maximum_dominant_eigenpair_residual",
      "maximum_requested_cin_error"
    ),
    metric = c(
      "State vs dominant at 40 tau",
      "Mean vs dominant at 40 tau",
      "Ward vs Higham state",
      "Dominant eigenpair residual",
      "Requested p_miss,eff error"
    ),
    threshold = c(
      thresholds[["terminal_total_variation"]],
      thresholds[["terminal_mean_ploidy_error"]],
      thresholds[["solver_total_variation"]],
      thresholds[["dominant_eigenpair_residual"]],
      thresholds[["requested_cin_error"]]
    ),
    stringsAsFactors = FALSE
  )
  long <- do.call(rbind, lapply(seq_len(nrow(metric_map)), function(index) {
    data.frame(
      group = maxima$group,
      endpoint_scope = maxima$endpoint_scope,
      metric = metric_map$metric[[index]],
      value = as.numeric(maxima[[metric_map$column[[index]]]]),
      threshold = metric_map$threshold[[index]],
      stringsAsFactors = FALSE
    )
  }))
  long$plot_value <- pmax(long$value, 1e-16)
  long$metric <- factor(long$metric, levels = metric_map$metric)
  ggplot2::ggplot(
    long,
    ggplot2::aes(group, plot_value, shape = endpoint_scope)
  ) +
    ggplot2::geom_hline(
      data = unique(long[, c("metric", "threshold")]),
      ggplot2::aes(yintercept = threshold),
      color = "#8B1A1A", linetype = 2, linewidth = 0.55
    ) +
    ggplot2::geom_point(
      size = 2.4, stroke = 0.55, color = "#1F4E79", fill = "white"
    ) +
    ggplot2::facet_wrap(~metric, nrow = 1, scales = "free_y") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_shape_manual(
      values = c("Separate in vitro" = 21, "Figure 6A joint ensemble" = 24)
    ) +
    ggplot2::labs(
      title = "D  Numerical error bounds are satisfied across fitted regimes",
      subtitle = "Points show group maxima; red dashed lines are acceptance thresholds",
      x = NULL,
      y = "Maximum absolute error",
      shape = "Audit scope"
    ) +
    s65_theme(8.3) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
    )
}

s65_draw <- function(workspace_root = f6r_find_workspace_root()) {
  f6r_require_packages(c("ggplot2", "patchwork", "scales"))
  Sys.setenv(XDG_CACHE_HOME = tempdir())
  paths <- s65_paths(workspace_root)
  required <- file.path(paths$data, c(
    "experimental_4N_deprived_distribution.tsv",
    "experimental_4N_deprived_summary.tsv",
    "separate_seed_fixed_operator_trajectory.tsv",
    "seed10_distribution_snapshots.tsv",
    "numerical_consistency_group_maxima.tsv",
    "supp_figure6-5_validation.tsv"
  ))
  f6r_require_files(required, "Supplementary Figure 6-5 drawing input")
  validation <- f6r_read_tsv(file.path(
    paths$data, "supp_figure6-5_validation.tsv"
  ))
  if (!all(validation$passed)) {
    stop("Cannot draw Supplementary Figure 6-5 from failed validation data.")
  }
  dir.create(paths$panels, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(paths$output_png), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(paths$manuscript_png), recursive = TRUE, showWarnings = FALSE)

  panels <- list(
    A = s65_panel_a(paths),
    B = s65_panel_b(paths),
    C = s65_panel_c(paths),
    D = s65_panel_d(paths)
  )
  panel_dimensions <- list(
    A = c(10.8, 3.25), B = c(10.8, 4.4),
    C = c(10.8, 3.2), D = c(10.8, 3.1)
  )
  panel_manifest <- lapply(names(panels), function(label) {
    png <- file.path(paths$panels, paste0("supp_fig6-5", tolower(label), ".png"))
    pdf <- file.path(paths$panels, paste0("supp_fig6-5", tolower(label), ".pdf"))
    dims <- panel_dimensions[[label]]
    ggplot2::ggsave(
      png, panels[[label]], width = dims[[1L]], height = dims[[2L]],
      units = "in", dpi = 300, bg = "white"
    )
    ggplot2::ggsave(
      pdf, panels[[label]], width = dims[[1L]], height = dims[[2L]],
      units = "in", device = grDevices::cairo_pdf, bg = "white"
    )
    data.frame(
      panel = label, png = normalizePath(png), pdf = normalizePath(pdf),
      stringsAsFactors = FALSE
    )
  })
  panel_manifest <- do.call(rbind, panel_manifest)
  assembled <-
    panels$A / panels$B / panels$C / panels$D +
    patchwork::plot_layout(heights = c(0.88, 1.16, 0.88, 0.86))
  ggplot2::ggsave(
    paths$output_png, assembled,
    width = 11.2, height = 14.8, units = "in", dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    paths$output_pdf, assembled,
    width = 11.2, height = 14.8, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
  if (!file.copy(paths$output_png, paths$manuscript_png, overwrite = TRUE) ||
      !file.copy(paths$output_pdf, paths$manuscript_pdf, overwrite = TRUE)) {
    stop("Failed to publish Supplementary Figure 6-5 to manuscript figures.")
  }
  output_files <- c(
    output_png = paths$output_png,
    output_pdf = paths$output_pdf,
    manuscript_png = paths$manuscript_png,
    manuscript_pdf = paths$manuscript_pdf
  )
  hashes <- vapply(output_files, f6r_md5, character(1L))
  validation_out <- data.frame(
    check = c(
      "panel_count", "published_png_identity", "published_pdf_identity",
      "all_outputs_nonempty"
    ),
    observed = c(
      nrow(panel_manifest),
      identical(hashes[["output_png"]], hashes[["manuscript_png"]]),
      identical(hashes[["output_pdf"]], hashes[["manuscript_pdf"]]),
      all(file.info(output_files)$size > 0)
    ),
    expected = c(4, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  validation_out$passed <- as.character(validation_out$observed) ==
    as.character(validation_out$expected)
  if (!all(validation_out$passed)) {
    stop("Supplementary Figure 6-5 drawing validation failed.")
  }
  f6r_write_tsv(
    panel_manifest,
    file.path(paths$data, "supp_figure6-5_panel_manifest.tsv")
  )
  f6r_write_tsv(
    validation_out,
    file.path(paths$data, "supp_figure6-5_draw_validation.tsv")
  )
  output_manifest <- data.frame(
    artifact = names(output_files),
    path = normalizePath(output_files, mustWork = TRUE),
    md5 = hashes,
    stringsAsFactors = FALSE
  )
  f6r_write_tsv(
    output_manifest,
    file.path(paths$data, "supp_figure6-5_output_manifest.tsv")
  )
  invisible(list(
    paths = paths,
    panels = panels,
    assembled = assembled,
    validation = validation_out,
    output_manifest = output_manifest
  ))
}
