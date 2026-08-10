#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DEoptim))
suppressPackageStartupMessages(library(dplyr))

.o2_bootstrap_script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) > 0L) {
    return(dirname(frame_files[[length(frame_files)]]))
  }
  getwd()
})
SCRIPT_DIR <- normalizePath(.o2_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
OXYGEN_ROOT <- normalizePath(file.path(WORKFLOW_ROOT, "..", ".."), mustWork = FALSE)

source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(WORKFLOW_ROOT, "model", "model_O2_supply_demand_MAP.R"), local = environment())
source(
  file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_invitro_utils.R"),
  local = environment(),
  chdir = TRUE
)

`%||%` <- if (exists("%||%", inherits = TRUE)) get("%||%", inherits = TRUE) else function(x, y) {
  if (is.null(x) || !length(x)) y else x
}
rm(.o2_bootstrap_script_dir)

parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool
.first_non_null_local <- o2sd_first_non_null

default_parameter_table_path <- function(script_dir = SCRIPT_DIR,
                                         must_exist = FALSE) {
  path <- ivt_parameter_table_path(
    repo_root = OXYGEN_ROOT
  )
  if (isTRUE(must_exist) && !file.exists(path)) {
    stop("Default in vitro parameter table not found: ", path)
  }
  normalizePath(path, mustWork = FALSE)
}

default_fit_objects_dir <- function(script_dir = SCRIPT_DIR, must_exist = FALSE) {
  path <- normalizePath(
    file.path(OXYGEN_ROOT, "ploidyOxygen", "data", "fit_objects"),
    mustWork = FALSE
  )
  if (isTRUE(must_exist) && !dir.exists(path)) {
    stop("Default in vitro fit_objects directory not found: ", path)
  }
  path
}

default_flow_density_path <- function(script_dir = SCRIPT_DIR) {
  normalizePath(
    file.path(OXYGEN_ROOT, "data", "g0g1_ploidy_density_grid.csv"),
    mustWork = FALSE
  )
}

default_death_data_path <- function(script_dir = SCRIPT_DIR) {
  normalizePath(
    file.path(
      OXYGEN_ROOT,
      "..",
      "data",
      "InVitroData",
      "sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"
    ),
    mustWork = FALSE
  )
}

default_out_dir <- function(script_dir = SCRIPT_DIR) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  normalizePath(
    file.path(OXYGEN_ROOT, "results", paste0("fit_model_O2_supply_demand_MAP_invitro_", stamp)),
    mustWork = FALSE
  )
}

normalize_invitro_n_cores <- o2sd_normalize_n_cores

build_invitro_de_initial_population <- function(NP, lower, upper, init) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  init <- as.numeric(init)
  if (length(lower) == 0L || length(upper) != length(lower) ||
      length(init) != length(lower)) {
    stop("DE initial-population vectors must have the same non-zero length.")
  }
  if (any(!is.finite(lower)) || any(!is.finite(upper)) ||
      any(!is.finite(init)) || any(lower > upper)) {
    stop("DE initial-population bounds and init values must be finite and ordered.")
  }
  NP <- as.integer(NP)
  if (length(NP) != 1L || is.na(NP) || NP < 1L) {
    stop("DE initial-population NP must be one positive integer.")
  }
  population <- matrix(
    stats::runif(NP * length(lower)),
    nrow = NP,
    ncol = length(lower)
  )
  population <- sweep(population, 2L, upper - lower, `*`)
  population <- sweep(population, 2L, lower, `+`)
  population[1L, ] <- pmin(pmax(init, lower), upper)
  population
}

start_invitro_deoptim_cluster <- function(n_cores) {
  n_use <- normalize_invitro_n_cores(n_cores)
  if (n_use <= 1L) return(NULL)
  if (.Platform$OS.type == "unix" && exists("makeForkCluster", envir = asNamespace("parallel"), inherits = FALSE)) {
    return(parallel::makeForkCluster(n_use))
  }
  parallel::makePSOCKcluster(n_use)
}

INVITRO_DE_PREFLIGHT_MAX_RETRIES <- 5L
INVITRO_DE_PENALTY_OBJECTIVE <- 1e9

invitro_de_preflight_text <- function(x) {
  value <- if (is.null(x) || !length(x)) NA_character_ else as.character(x[[1L]])
  if (is.na(value)) return(NA_character_)
  trimws(gsub("[\t\r\n]+", " ", value))
}

write_invitro_de_preflight_audit <- function(audit, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    audit,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
  invisible(path)
}

initialize_invitro_deoptim_workers <- function(cluster, cpp_info) {
  wrapper_path <- as.character(cpp_info$wrapper_path %||% "")
  dll_path <- as.character(cpp_info$path %||% "")
  required_cpp <- c(
    "cpp_o2simps_build_G_for_o2_triplet",
    "cpp_o2simps_simulate_one",
    "cpp_o2simps_objective_components_map"
  )
  worker_results <- parallel::clusterCall(
    cluster,
    function(wrapper_path, dll_path, required_cpp) {
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1"
      )
      tryCatch({
        if (!nzchar(wrapper_path) || !file.exists(wrapper_path)) {
          stop("C++ wrapper does not exist: ", wrapper_path)
        }
        if (!nzchar(dll_path) || !file.exists(dll_path)) {
          stop("C++ DLL does not exist: ", dll_path)
        }

        # The generated sourceCpp wrapper begins with dyn.load(dll_path), then
        # reconstructs process-local Rcpp functions bound to that DLL.
        source(wrapper_path, local = .GlobalEnv)
        missing_cpp <- required_cpp[!vapply(
          required_cpp,
          exists,
          logical(1),
          envir = .GlobalEnv,
          mode = "function",
          inherits = TRUE
        )]
        if (length(missing_cpp)) {
          stop(
            "Worker C++ wrapper is missing required functions: ",
            paste(missing_cpp, collapse = ", ")
          )
        }
        build_formals <- names(formals(get(
          "cpp_o2simps_build_G_for_o2_triplet",
          envir = .GlobalEnv,
          inherits = TRUE
        )))
        sim_formals <- names(formals(get(
          "cpp_o2simps_simulate_one",
          envir = .GlobalEnv,
          inherits = TRUE
        )))
        objective_formals <- names(formals(get(
          "cpp_o2simps_objective_components_map",
          envir = .GlobalEnv,
          inherits = TRUE
        )))
        if (!("p_wgd" %in% build_formals) ||
            !("sim_args" %in% sim_formals) ||
            !all(c("scenario_data", "objective_data", "state_data", "sim_args") %in%
              objective_formals)) {
          stop("Worker C++ wrapper has stale or incompatible function signatures.")
        }
        list(
          ok = TRUE,
          status = "PASS_WRAPPER_DLL_LOADED",
          error = NA_character_
        )
      }, error = function(e) {
        list(
          ok = FALSE,
          status = "WORKER_CPP_INIT_ERROR",
          error = conditionMessage(e)
        )
      })
    },
    wrapper_path,
    dll_path,
    required_cpp
  )

  do.call(rbind, lapply(seq_along(worker_results), function(i) {
    result <- worker_results[[i]]
    data.frame(
      worker = as.integer(i),
      role = "worker_init",
      ok = isTRUE(result$ok),
      status = invitro_de_preflight_text(result$status),
      objective = NA_real_,
      penalty_reason = NA_character_,
      error = invitro_de_preflight_text(result$error),
      stringsAsFactors = FALSE
    )
  }))
}

invitro_de_preflight_evaluate <- function(fn,
                                          par,
                                          penalty_value = INVITRO_DE_PENALTY_OBJECTIVE) {
  tryCatch({
    comp <- fn(par)
    objective <- suppressWarnings(as.numeric(comp$objective))
    penalty_reason <- if (!is.null(comp$penalty_reason) && length(comp$penalty_reason)) {
      as.character(comp$penalty_reason[[1L]])
    } else {
      NA_character_
    }
    protocol_penalty <- !is.na(penalty_reason) &&
      grepl("^protocol_infeasible:", trimws(penalty_reason))
    objective_numeric_ok <- length(objective) == 1L && is.finite(objective)
    objective_fit_ok <- objective_numeric_ok && objective < penalty_value
    objective_environment_ok <- objective_fit_ok ||
      (objective_numeric_ok && protocol_penalty && objective >= penalty_value)
    reason_environment_ok <- is.na(penalty_reason) ||
      !nzchar(trimws(penalty_reason)) || protocol_penalty
    status <- if (protocol_penalty && objective_numeric_ok) {
      "MODEL_PROTOCOL_INFEASIBLE"
    } else if (!objective_fit_ok) {
      if (length(objective) == 1L && is.finite(objective) && objective >= penalty_value) {
        "PENALTY_OBJECTIVE"
      } else {
        "INVALID_OBJECTIVE"
      }
    } else if (!reason_environment_ok) {
      "PENALTY_REASON"
    } else {
      "PASS"
    }
    list(
      ok = isTRUE(objective_environment_ok) && isTRUE(reason_environment_ok),
      status = status,
      objective = if (length(objective) == 1L) objective else NA_real_,
      penalty_reason = penalty_reason,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      ok = FALSE,
      status = "OBJECTIVE_ERROR",
      objective = NA_real_,
      penalty_reason = NA_character_,
      error = conditionMessage(e)
    )
  })
}

run_invitro_de_worker_preflight <- function(cluster,
                                            objective_from_free,
                                            objective_value,
                                            init_free,
                                            penalty_objective = INVITRO_DE_PENALTY_OBJECTIVE,
                                            objective_tolerance = 1e-8) {
  master_probe <- invitro_de_preflight_evaluate(
    objective_from_free,
    init_free,
    penalty_objective
  )
  worker_probes <- tryCatch(
    parallel::clusterCall(
      cluster,
      function(fn, par, penalty_value, evaluator) {
        evaluator(fn, par, penalty_value)
      },
      objective_from_free,
      init_free,
      penalty_objective,
      invitro_de_preflight_evaluate
    ),
    error = function(e) {
      list(list(
        ok = FALSE,
        status = "CLUSTER_CALL_ERROR",
        objective = NA_real_,
        penalty_reason = NA_character_,
        error = conditionMessage(e)
      ))
    }
  )

  deoptim_path_error <- NULL
  deoptim_path_values <- tryCatch({
    parallel::clusterExport(
      cluster,
      varlist = "objective_value",
      envir = environment()
    )
    probe_matrix <- matrix(
      rep(as.numeric(init_free), times = length(cluster)),
      nrow = length(cluster),
      byrow = TRUE
    )
    colnames(probe_matrix) <- names(init_free)
    as.numeric(parallel::parApply(
      cluster,
      probe_matrix,
      MARGIN = 1L,
      FUN = objective_value
    ))
  }, error = function(e) {
    deoptim_path_error <<- conditionMessage(e)
    numeric(0)
  })

  probes <- c(list(master_probe), worker_probes)
  worker_ids <- c(0L, seq_along(worker_probes))
  rows <- do.call(rbind, lapply(seq_along(probes), function(i) {
    probe <- probes[[i]]
    data.frame(
      worker = worker_ids[[i]],
      role = if (worker_ids[[i]] == 0L) "master" else "worker",
      ok = isTRUE(probe$ok),
      status = invitro_de_preflight_text(probe$status),
      objective = suppressWarnings(as.numeric(probe$objective %||% NA_real_)),
      penalty_reason = invitro_de_preflight_text(probe$penalty_reason),
      error = invitro_de_preflight_text(probe$error),
      stringsAsFactors = FALSE
    )
  }))
  deoptim_rows <- if (length(deoptim_path_values) == length(cluster)) {
    protocol_penalty_ok <- identical(
      as.character(master_probe$status),
      "MODEL_PROTOCOL_INFEASIBLE"
    )
    deoptim_ok <- is.finite(deoptim_path_values) &
      (deoptim_path_values < penalty_objective |
         (protocol_penalty_ok &
            abs(deoptim_path_values - master_probe$objective) <= objective_tolerance *
              max(1, abs(master_probe$objective))))
    data.frame(
      worker = seq_along(deoptim_path_values),
      role = "deoptim_path",
      ok = deoptim_ok,
      status = ifelse(
        deoptim_ok,
        if (protocol_penalty_ok) "MODEL_PROTOCOL_INFEASIBLE" else "PASS",
        "PENALTY_OR_INVALID_OBJECTIVE"
      ),
      objective = deoptim_path_values,
      penalty_reason = NA_character_,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      worker = NA_integer_,
      role = "deoptim_path",
      ok = FALSE,
      status = "DEOPTIM_PATH_ERROR",
      objective = NA_real_,
      penalty_reason = NA_character_,
      error = deoptim_path_error %||% "DEoptim-style preflight returned the wrong number of values",
      stringsAsFactors = FALSE
    )
  }
  rows <- rbind(rows, deoptim_rows)

  if (isTRUE(master_probe$ok) && all(rows$ok)) {
    comparable_rows <- which(rows$role %in% c("worker", "deoptim_path"))
    worker_values <- rows$objective[comparable_rows]
    scale <- max(1, abs(master_probe$objective))
    mismatch <- !is.finite(worker_values) |
      abs(worker_values - master_probe$objective) > objective_tolerance * scale
    if (any(mismatch)) {
      bad_rows <- comparable_rows[mismatch]
      rows$ok[bad_rows] <- FALSE
      rows$status[bad_rows] <- "OBJECTIVE_MISMATCH"
      rows$error[bad_rows] <- paste0(
        "worker objective differs from master objective ",
        format(master_probe$objective, digits = 17)
      )
    }
  }

  rows
}

start_invitro_deoptim_cluster_with_preflight <- function(
    n_cores,
    objective_from_free,
    objective_value,
    init_free,
    cpp_info,
    audit_path,
    max_retries = INVITRO_DE_PREFLIGHT_MAX_RETRIES,
    cluster_factory = start_invitro_deoptim_cluster,
    cluster_stopper = function(cluster) parallel::stopCluster(cluster),
    worker_initializer = initialize_invitro_deoptim_workers,
    preflight_runner = run_invitro_de_worker_preflight,
    sleep_fn = Sys.sleep) {
  retries_use <- suppressWarnings(as.integer(max_retries))
  if (!is.finite(retries_use) || is.na(retries_use) || retries_use < 0L) {
    retries_use <- INVITRO_DE_PREFLIGHT_MAX_RETRIES
  }
  max_attempts <- retries_use + 1L
  audit_rows <- list()

  append_audit <- function(rows, attempt) {
    for (column in c("status", "penalty_reason", "error")) {
      rows[[column]] <- vapply(
        rows[[column]],
        invitro_de_preflight_text,
        character(1)
      )
    }
    rows$attempt <- as.integer(attempt)
    rows$retry_number <- as.integer(attempt - 1L)
    rows$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    rows <- rows[, c(
      "timestamp", "attempt", "retry_number", "worker", "role", "ok",
      "status", "objective", "penalty_reason", "error"
    ), drop = FALSE]
    audit_rows[[length(audit_rows) + 1L]] <<- rows
    audit <- do.call(rbind, audit_rows)
    write_invitro_de_preflight_audit(audit, audit_path)
    audit
  }

  last_error <- "unknown preflight failure"
  for (attempt in seq_len(max_attempts)) {
    cluster <- NULL
    start_error <- NULL
    cluster <- tryCatch(
      cluster_factory(n_cores),
      error = function(e) {
        start_error <<- conditionMessage(e)
        NULL
      }
    )

    if (is.null(cluster)) {
      last_error <- if (!is.null(start_error)) start_error else "cluster factory returned NULL"
      attempt_rows <- data.frame(
        worker = NA_integer_,
        role = "cluster",
        ok = FALSE,
        status = "CLUSTER_START_ERROR",
        objective = NA_real_,
        penalty_reason = NA_character_,
        error = last_error,
        stringsAsFactors = FALSE
      )
    } else {
      init_rows <- tryCatch(
        worker_initializer(
          cluster = cluster,
          cpp_info = cpp_info
        ),
        error = function(e) {
          data.frame(
            worker = NA_integer_,
            role = "worker_init",
            ok = FALSE,
            status = "WORKER_CPP_INIT_ERROR",
            objective = NA_real_,
            penalty_reason = NA_character_,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          )
        }
      )
      if (is.data.frame(init_rows) && nrow(init_rows) &&
          all(c("ok", "status", "error") %in% names(init_rows)) &&
          isTRUE(all(init_rows$ok))) {
        preflight_rows <- tryCatch(
          preflight_runner(
            cluster = cluster,
            objective_from_free = objective_from_free,
            objective_value = objective_value,
            init_free = init_free
          ),
          error = function(e) {
            data.frame(
              worker = NA_integer_,
              role = "cluster",
              ok = FALSE,
              status = "PREFLIGHT_ERROR",
              objective = NA_real_,
              penalty_reason = NA_character_,
              error = conditionMessage(e),
              stringsAsFactors = FALSE
            )
          }
        )
        attempt_rows <- rbind(init_rows, preflight_rows)
      } else {
        attempt_rows <- init_rows
      }
      required_preflight_columns <- c(
        "worker", "role", "ok", "status", "objective",
        "penalty_reason", "error"
      )
      if (!is.data.frame(attempt_rows) || !nrow(attempt_rows) ||
          !all(required_preflight_columns %in% names(attempt_rows))) {
        attempt_rows <- data.frame(
          worker = NA_integer_,
          role = "cluster",
          ok = FALSE,
          status = "INVALID_PREFLIGHT_RESULT",
          objective = NA_real_,
          penalty_reason = NA_character_,
          error = "preflight runner returned an invalid audit table",
          stringsAsFactors = FALSE
        )
      }
      failed_rows <- attempt_rows[!attempt_rows$ok, , drop = FALSE]
      if (nrow(failed_rows)) {
        detail <- failed_rows$error
        detail[is.na(detail) | !nzchar(detail)] <- failed_rows$penalty_reason[
          is.na(detail) | !nzchar(detail)
        ]
        detail[is.na(detail) | !nzchar(detail)] <- failed_rows$status[
          is.na(detail) | !nzchar(detail)
        ]
        last_error <- paste(unique(detail), collapse = " | ")
      }
    }

    audit <- append_audit(attempt_rows, attempt)
    if (!is.null(cluster) && isTRUE(all(attempt_rows$ok))) {
      message(
        "[fit_invitro] DEoptim worker preflight passed on attempt ", attempt,
        " (retries=", attempt - 1L, ")."
      )
      return(list(
        cluster = cluster,
        active_cores = length(cluster),
        attempts = as.integer(attempt),
        retries = as.integer(attempt - 1L),
        audit = audit
      ))
    }

    if (!is.null(cluster)) {
      try(cluster_stopper(cluster), silent = TRUE)
    }
    if (attempt < max_attempts) {
      message(
        "[fit_invitro] DEoptim worker preflight failed on attempt ", attempt,
        "; retrying with a fresh cluster (retry ", attempt, "/", retries_use,
        "). Reason: ", last_error
      )
      sleep_fn(min(0.25 * attempt, 1))
    }
  }

  stop(
    "[fit_invitro] DEoptim worker preflight failed after ", max_attempts,
    " attempts (", retries_use, " retries). See ", audit_path,
    ". Last error: ", last_error,
    call. = FALSE
  )
}

resolve_optional_flow_density_path <- function(raw_path = NULL) {
  if (!is.null(raw_path)) {
    path <- normalizePath(raw_path, mustWork = FALSE)
    if (!file.exists(path)) {
      stop("flow_density_path not found: ", path)
    }
    return(path)
  }
  default_path <- default_flow_density_path()
  if (file.exists(default_path)) default_path else default_path
}

resolve_death_data_path <- function(raw_path = NULL) {
  path <- if (is.null(raw_path) || !length(raw_path) ||
              is.na(raw_path[[1]]) || !nzchar(trimws(raw_path[[1]]))) {
    default_death_data_path()
  } else {
    normalizePath(raw_path[[1]], mustWork = FALSE)
  }
  if (!file.exists(path)) {
    stop("death_data_path not found: ", path)
  }
  path
}

ivt_load_fit_objects_compat <- function(fit_objects_dir,
                                        flow_density_path = NULL,
                                        death_data_path = NULL) {
  load_formals <- names(formals(ivt_load_fit_objects))
  call_args <- list(repo_root = OXYGEN_ROOT)
  if ("fit_objects_dir" %in% load_formals) {
    call_args$fit_objects_dir <- fit_objects_dir
  }
  if ("flow_csv_path" %in% load_formals) {
    call_args$flow_csv_path <- flow_density_path
  }
  if ("death_data_path" %in% load_formals) {
    call_args$death_data_path <- death_data_path
  }
  do.call(ivt_load_fit_objects, call_args)
}

build_invitro_cfg <- function(parameter_table,
                              dt = 0.05,
                              init_total_size = 1e6,
                              o2_upper_bound = 21,
                              fixed_oxygen = TRUE) {
  cfg <- ivt_build_default_cfg(
    repo_root = OXYGEN_ROOT,
    dt = dt,
    init_total_size = init_total_size,
    o2_upper_bound = o2_upper_bound,
    fixed_oxygen = fixed_oxygen
  )
  cfg$parameter_table <- normalizePath(parameter_table, mustWork = FALSE)
  cfg <- normalize_sim_cfg_common(cfg, context = "fit")
  cfg
}

validate_invitro_parameter_table <- function(parameter_table,
                                             dt = 0.05,
                                             init_total_size = 1e6,
                                             o2_upper_bound = 21,
                                             fixed_oxygen = TRUE) {
  if (!file.exists(parameter_table)) {
    stop("In vitro parameter table not found: ", parameter_table)
  }
  cfg <- build_invitro_cfg(
    parameter_table = parameter_table,
    dt = dt,
    init_total_size = init_total_size,
    o2_upper_bound = o2_upper_bound,
    fixed_oxygen = fixed_oxygen
  )
  ivt_optimizer_spec(cfg)
  invisible(cfg)
}

validate_invitro_fit_objects <- function(fit_objects_dir,
                                         flow_density_path = NULL,
                                         death_data_path = NULL) {
  ivt_load_fit_objects_compat(
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path,
    death_data_path = death_data_path
  )
  invisible(TRUE)
}

make_penalty_components <- function(objective = 1e9, reason = "penalty") {
  empty_summary <- data.frame()
  empty_result <- list(
    adapter = NULL,
    model_core = NULL,
    grid_pre = integer(0),
    segment_results = list()
  )
  list(
    objective = as.numeric(objective),
    likelihood_objective = as.numeric(objective),
    total_loglik = -as.numeric(objective),
    growth_loglik = -as.numeric(objective),
    ploidy_loglik = 0.0,
    flow_loglik = 0.0,
    death_loglik = 0.0,
    passage_time_loglik = 0.0,
    growth_loglik_sum = -as.numeric(objective),
    growth_passage_loglik_sum = -as.numeric(objective),
    ploidy_loglik_sum = 0.0,
    flow_loglik_sum = 0.0,
    death_loglik_sum = 0.0,
    passage_time_loglik_sum = 0.0,
    buffer_prior_weight = 0.0,
    buffer_prior_raw_penalty = 0.0,
    buffer_prior_penalty = 0.0,
    buffer_prior_terms = data.frame(),
    sigma_growth = NA_real_,
    sigma_kary = NA_real_,
    sigma_flow_ploidy = NA_real_,
    sigma_death_logit = NA_real_,
    death_fraction_eps = NA_real_,
    death_weight = NA_real_,
    passage_time_weight = NA_real_,
    passage_time_tolerance_days = NA_real_,
    passage_time_sigma_days = NA_real_,
    passage_time_df_value = NA_real_,
    n_growth = 0L,
    n_growth_timepoints = 0L,
    n_growth_passages = 0L,
    n_growth_observed = 0L,
    n_growth_missing_pred = 0L,
    n_growth_negative_pred = 0L,
    n_ploidy_passages = 0L,
    n_kary_cells = 0L,
    n_flow_passages = 0L,
    n_flow_samples = 0L,
    n_death_passages = 0L,
    n_passage_time_observations = 0L,
    n_scenarios = 0L,
    n_insufficient_boundaries = 0L,
    all_passage_boundaries_feasible = FALSE,
    protocol_feasibility_status = "UNKNOWN",
    max_boundary_scale = NA_real_,
    summary = empty_summary,
    growth_df = data.frame(),
    growth_count_diagnostics_df = data.frame(),
    growth_passage_df = data.frame(),
    ploidy_df = data.frame(),
    flow_df = data.frame(),
    death_df = data.frame(),
    passage_time_df = data.frame(),
    flow_overlay_df = data.frame(),
    objective_hierarchy = data.frame(),
    growth_lineage_loglik = data.frame(),
    growth_cohort_loglik = data.frame(),
    ploidy_lineage_loglik = data.frame(),
    ploidy_cohort_loglik = data.frame(),
    flow_lineage_loglik = data.frame(),
    flow_cohort_loglik = data.frame(),
    death_lineage_loglik = data.frame(),
    death_cohort_loglik = data.frame(),
    passage_time_lineage_loglik = data.frame(),
    passage_time_cohort_loglik = data.frame(),
    run_2N = empty_result,
    run_4N = empty_result,
    penalty_reason = as.character(reason)
  )
}

invitro_protocol_penalty_objective <- function(condition,
                                                base_penalty = 1e6,
                                                remaining_segment_penalty = 1e4,
                                                relative_shortfall_penalty = 1e3) {
  if (!inherits(condition, "invitro_protocol_infeasible")) {
    return(INVITRO_DE_PENALTY_OBJECTIVE)
  }
  ordinal <- suppressWarnings(as.integer(condition$segment_ordinal))
  cohort_count <- suppressWarnings(as.integer(condition$segment_count))
  cohort <- as.character(condition$segment$cohort %||% "")
  if (length(ordinal) != 1L || !is.finite(ordinal) || ordinal < 1L ||
      length(cohort_count) != 1L || !is.finite(cohort_count) || cohort_count < 1L) {
    return(INVITRO_DE_PENALTY_OBJECTIVE - 1)
  }
  prior_completed <- if (identical(cohort, "4N")) cohort_count else 0L
  completed_segments <- prior_completed + ordinal - 1L
  total_segments <- 2L * cohort_count
  remaining_segments <- max(total_segments - completed_segments, 1L)
  threshold <- suppressWarnings(as.numeric(condition$selection$threshold_target_cells))
  max_live <- suppressWarnings(as.numeric(condition$selection$max_live_cells_in_search))
  relative_shortfall <- if (length(threshold) == 1L && is.finite(threshold) && threshold > 0 &&
                            length(max_live) == 1L && is.finite(max_live)) {
    max((threshold - max_live) / threshold, 0)
  } else {
    1
  }
  score <- as.numeric(base_penalty) +
    as.numeric(remaining_segment_penalty) * remaining_segments +
    as.numeric(relative_shortfall_penalty) * min(relative_shortfall, 10)
  min(score, INVITRO_DE_PENALTY_OBJECTIVE - 1)
}

write_tsv_if_nonempty <- o2sd_write_tsv_if_nonempty

run_invitro_rscript_to_log <- function(label, script_path, args, log_path) {
  if (!file.exists(script_path)) {
    stop("Missing ", label, " script: ", script_path, call. = FALSE)
  }
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "Rscript",
    args = c(normalizePath(script_path, mustWork = TRUE), args),
    stdout = log_path,
    stderr = log_path
  )
  if (!identical(status, 0L)) {
    tail_txt <- if (file.exists(log_path)) {
      paste(utils::tail(readLines(log_path, warn = FALSE), 25L), collapse = "\n")
    } else {
      ""
    }
    stop(
      label,
      " failed with exit status ",
      status,
      ". See ",
      log_path,
      if (nzchar(tail_txt)) paste0("\nLast log lines:\n", tail_txt) else "",
      call. = FALSE
    )
  }
  invisible(status)
}

run_invitro_auto_viz_report <- function(out_dir) {
  postfit_script <- file.path(
    WORKFLOW_ROOT,
    "runner",
    "run_postfit_pipeline.R"
  )
  run_invitro_rscript_to_log(
    label = "invitro post-fit pipeline",
    script_path = postfit_script,
    args = c(
      paste0("--fit_dir=", normalizePath(out_dir, mustWork = FALSE)),
      "--scope=invitro"
    ),
    log_path = file.path(out_dir, "postfit_status.log")
  )
  invisible(TRUE)
}

invitro_prov_cell <- o2sd_provenance_cell
invitro_command_text <- o2sd_command_text

invitro_parse_effective_args <- function(args, source = "fit_command") {
  rows <- list()
  if (length(args) && !startsWith(args[[1]], "--")) {
    rows[[length(rows) + 1L]] <- data.frame(source = source, key = "script", value = args[[1]], stringsAsFactors = FALSE)
    args <- args[-1L]
  }
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      rows[[length(rows) + 1L]] <- data.frame(source = source, key = paste0("positional_", i), value = arg, stringsAsFactors = FALSE)
      i <- i + 1L
      next
    }
    stripped <- sub("^--", "", arg)
    if (grepl("=", stripped, fixed = TRUE)) {
      rows[[length(rows) + 1L]] <- data.frame(
        source = source,
        key = sub("=.*$", "", stripped),
        value = sub("^[^=]*=", "", stripped),
        stringsAsFactors = FALSE
      )
      i <- i + 1L
    } else {
      value <- "TRUE"
      if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
        value <- args[[i + 1L]]
        i <- i + 1L
      }
      rows[[length(rows) + 1L]] <- data.frame(source = source, key = stripped, value = value, stringsAsFactors = FALSE)
      i <- i + 1L
    }
  }
  if (!length(rows)) return(data.frame(source = character(0), key = character(0), value = character(0)))
  out <- do.call(rbind, rows)
  out[] <- lapply(out, invitro_prov_cell)
  out
}

write_invitro_run_provenance <- function(out_dir, argv, parameter_table,
                                         fit_objects_dir, flow_density_path,
                                         death_data_path, death_weight,
                                         ploidy_weight,
                                         sigma_death_logit,
                                         death_fraction_eps,
                                         buffer_prior_weight,
                                         buffer_prior_center_smax,
                                         buffer_prior_sd_smax,
                                         buffer_prior_center_log10_beta,
                                         buffer_prior_sd_log10_beta,
                                         buffer_prior_center_log10_n_exp,
                                         buffer_prior_sd_log10_n_exp,
                                         passage_time_weight,
                                         passage_time_tolerance_days,
                                         passage_time_sigma_days,
                                         passage_time_df,
                                         seed, itermax, NP,
                                         de_reltol, de_steptol, n_cores,
                                         de_include_parameter_init) {
  command_text <- Sys.getenv("O2SD_RUN_COMMAND", unset = NA_character_)
  if (is.na(command_text) || !nzchar(command_text)) {
    command_text <- invitro_command_text("Rscript", commandArgs(trailingOnly = FALSE))
  }
  writeLines(command_text, file.path(out_dir, "fit_command.txt"), useBytes = TRUE)
  args <- c(
    "--fit_invitro",
    paste0("--seed=", seed),
    paste0("--out_dir=", out_dir),
    paste0("--parameter_table=", parameter_table),
    paste0("--fit_objects_dir=", fit_objects_dir),
    paste0("--itermax=", itermax),
    paste0("--NP=", NP),
    paste0("--de_reltol=", de_reltol),
    paste0("--de_steptol=", de_steptol),
    paste0("--de_include_parameter_init=", de_include_parameter_init),
    paste0("--n_cores=", n_cores),
    paste0("--death_data_path=", death_data_path),
    paste0("--death_weight=", death_weight),
    paste0("--ploidy_weight=", ploidy_weight),
    paste0("--sigma_death_logit=", sigma_death_logit),
    paste0("--death_fraction_eps=", death_fraction_eps),
    paste0("--buffer_prior_weight=", buffer_prior_weight),
    paste0("--buffer_prior_center_smax=", buffer_prior_center_smax),
    paste0("--buffer_prior_sd_smax=", buffer_prior_sd_smax),
    paste0("--buffer_prior_center_log10_beta=", buffer_prior_center_log10_beta),
    paste0("--buffer_prior_sd_log10_beta=", buffer_prior_sd_log10_beta),
    paste0("--buffer_prior_center_log10_n_exp=", buffer_prior_center_log10_n_exp),
    paste0("--buffer_prior_sd_log10_n_exp=", buffer_prior_sd_log10_n_exp),
    paste0("--passage_time_weight=", passage_time_weight),
    paste0("--passage_time_tolerance_days=", passage_time_tolerance_days),
    paste0("--passage_time_sigma_days=", passage_time_sigma_days),
    paste0("--passage_time_df=", passage_time_df)
  )
  if (!is.null(flow_density_path) && nzchar(flow_density_path)) {
    args <- c(args, paste0("--flow_density_path=", flow_density_path))
  }
  utils::write.table(
    invitro_parse_effective_args(args, source = "fit_command"),
    file = file.path(out_dir, "run_effective_args.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  prov <- data.frame(
    section = c(
      "execution", "execution", "execution", "execution",
      "scripts", "input_config", "input_config", "input_config", "input_config",
      "fit", "fit", "fit", "fit", "fit", "fit", "fit", "fit", "fit",
      "prior", "prior", "prior", "prior", "prior", "prior", "prior",
      "optimizer", "optimizer", "optimizer", "optimizer", "optimizer", "optimizer",
      "slurm", "slurm"
    ),
    key = c(
      "timestamp", "hostname", "user", "fit_command_file",
      "array_script", "parameter_table", "fit_objects_dir", "flow_density_path", "death_data_path",
      "seed", "death_weight", "ploidy_weight", "sigma_death_logit", "death_fraction_eps",
      "passage_time_weight", "passage_time_tolerance_days",
      "passage_time_sigma_days", "passage_time_df",
      "buffer_prior_weight", "buffer_prior_center_smax", "buffer_prior_sd_smax",
      "buffer_prior_center_log10_beta", "buffer_prior_sd_log10_beta",
      "buffer_prior_center_log10_n_exp", "buffer_prior_sd_log10_n_exp",
      "itermax", "NP", "de_reltol", "de_steptol", "de_include_parameter_init", "n_cores",
      "array_job_id", "array_task_id"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      Sys.info()[["nodename"]],
      Sys.info()[["user"]],
      file.path(out_dir, "fit_command.txt"),
      Sys.getenv("O2SD_ARRAY_SCRIPT", unset = NA_character_),
      parameter_table,
      fit_objects_dir,
      flow_density_path %||% "",
      death_data_path,
      seed,
      death_weight,
      ploidy_weight,
      sigma_death_logit,
      death_fraction_eps,
      passage_time_weight,
      passage_time_tolerance_days,
      passage_time_sigma_days,
      passage_time_df,
      buffer_prior_weight,
      buffer_prior_center_smax,
      buffer_prior_sd_smax,
      buffer_prior_center_log10_beta,
      buffer_prior_sd_log10_beta,
      buffer_prior_center_log10_n_exp,
      buffer_prior_sd_log10_n_exp,
      itermax,
      NP,
      de_reltol,
      de_steptol,
      de_include_parameter_init,
      n_cores,
      Sys.getenv("O2SD_SLURM_ARRAY_JOB_ID", unset = NA_character_),
      Sys.getenv("O2SD_SLURM_ARRAY_TASK_ID", unset = NA_character_)
    ),
    stringsAsFactors = FALSE
  )
  prov[] <- lapply(prov, invitro_prov_cell)
  utils::write.table(prov, file.path(out_dir, "run_provenance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(TRUE)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  parameter_table <- if (!is.null(argv$parameter_table)) {
    argv$parameter_table
  } else {
    default_parameter_table_path(must_exist = TRUE)
  }
  fit_objects_dir <- if (!is.null(argv$fit_objects_dir)) {
    argv$fit_objects_dir
  } else {
    default_fit_objects_dir(must_exist = TRUE)
  }
  flow_density_path <- resolve_optional_flow_density_path(argv$flow_density_path)
  death_data_path <- resolve_death_data_path(argv$death_data_path)
  out_dir <- if (!is.null(argv$out_dir)) argv$out_dir else default_out_dir()

  seed <- as.integer(.first_non_null_local(argv$seed, 1L))
  itermax_requested <- as.integer(.first_non_null_local(argv$itermax, 500L))
  itermax_max <- as.integer(.first_non_null_local(argv$itermax_max, 500L))
  if (!is.finite(itermax_requested) || is.na(itermax_requested)) itermax_requested <- 500L
  if (!is.finite(itermax_max) || is.na(itermax_max) || itermax_max < 1L) itermax_max <- 500L
  itermax <- min(max(itermax_requested, 1L), itermax_max)
  NP_requested <- as.integer(.first_non_null_local(argv$NP, 80L))
  n_cores_requested <- normalize_invitro_n_cores(.first_non_null_local(argv$n_cores, 1L))
  de_reltol <- as.numeric(.first_non_null_local(argv$de_reltol, 1e-4))
  de_steptol <- as.integer(.first_non_null_local(argv$de_steptol, 25L))
  de_include_parameter_init <- as_bool(
    .first_non_null_local(argv$de_include_parameter_init, FALSE),
    FALSE
  )
  if (!is.finite(de_reltol) || de_reltol <= 0) de_reltol <- 1e-4
  if (!is.finite(de_steptol) || is.na(de_steptol) || de_steptol < 1L) de_steptol <- 25L
  dt_use <- as.numeric(.first_non_null_local(argv$dt, 0.05))
  init_total_size_use <- as.numeric(.first_non_null_local(argv$init_total_size, 1e6))
  o2_upper_bound_use <- as.numeric(.first_non_null_local(argv$o2_upper_bound, 21))
  fixed_oxygen_use <- TRUE
  auto_viz <- as_bool(.first_non_null_local(argv$auto_viz, TRUE), TRUE)
  death_weight <- as.numeric(.first_non_null_local(argv$death_weight, 1))
  ploidy_weight <- as.numeric(.first_non_null_local(argv$ploidy_weight, 1))
  sigma_death_logit <- as.numeric(.first_non_null_local(argv$sigma_death_logit, 0.75))
  death_fraction_eps <- as.numeric(.first_non_null_local(argv$death_fraction_eps, 1e-4))
  buffer_prior_weight <- as.numeric(.first_non_null_local(argv$buffer_prior_weight, 0))
  buffer_prior_center_smax <- as.numeric(.first_non_null_local(argv$buffer_prior_center_smax, 0.98))
  buffer_prior_sd_smax <- as.numeric(.first_non_null_local(argv$buffer_prior_sd_smax, 0.10))
  buffer_prior_center_log10_beta <- as.numeric(.first_non_null_local(
    argv$buffer_prior_center_log10_beta,
    log10(0.7)
  ))
  buffer_prior_sd_log10_beta <- as.numeric(.first_non_null_local(argv$buffer_prior_sd_log10_beta, 0.30))
  buffer_prior_center_log10_n_exp <- as.numeric(.first_non_null_local(
    argv$buffer_prior_center_log10_n_exp,
    log10(5)
  ))
  buffer_prior_sd_log10_n_exp <- as.numeric(.first_non_null_local(argv$buffer_prior_sd_log10_n_exp, 0.30))
  passage_time_weight <- as.numeric(.first_non_null_local(argv$passage_time_weight, 0.25))
  passage_time_tolerance_days <- as.numeric(.first_non_null_local(argv$passage_time_tolerance_days, 1))
  passage_time_sigma_days <- as.numeric(.first_non_null_local(argv$passage_time_sigma_days, 1))
  passage_time_df <- as.numeric(.first_non_null_local(argv$passage_time_df, 4))
  if (length(death_weight) != 1L || !is.finite(death_weight) || death_weight < 0) {
    stop("death_weight must be one finite non-negative value.")
  }
  if (length(ploidy_weight) != 1L || !is.finite(ploidy_weight) || ploidy_weight < 0) {
    stop("ploidy_weight must be one finite non-negative value.")
  }
  if (length(sigma_death_logit) != 1L || !is.finite(sigma_death_logit) || sigma_death_logit <= 0) {
    stop("sigma_death_logit must be one finite strictly positive value.")
  }
  if (length(death_fraction_eps) != 1L || !is.finite(death_fraction_eps) ||
      death_fraction_eps <= 0 || death_fraction_eps >= 0.5) {
    stop("death_fraction_eps must be one finite value strictly between 0 and 0.5.")
  }
  if (length(buffer_prior_weight) != 1L || !is.finite(buffer_prior_weight) ||
      buffer_prior_weight < 0) {
    stop("buffer_prior_weight must be one finite non-negative value.")
  }
  if (length(buffer_prior_center_smax) != 1L || !is.finite(buffer_prior_center_smax) ||
      buffer_prior_center_smax < 0 || buffer_prior_center_smax > 1) {
    stop("buffer_prior_center_smax must be one finite value in [0, 1].")
  }
  prior_sd_values <- c(
    buffer_prior_sd_smax,
    buffer_prior_sd_log10_beta,
    buffer_prior_sd_log10_n_exp
  )
  if (any(!is.finite(prior_sd_values)) || any(prior_sd_values <= 0)) {
    stop("All buffer prior SD values must be finite and strictly positive.")
  }
  prior_center_values <- c(
    buffer_prior_center_log10_beta,
    buffer_prior_center_log10_n_exp
  )
  if (any(!is.finite(prior_center_values))) {
    stop("All transformed buffer prior centers must be finite.")
  }
  if (length(passage_time_weight) != 1L ||
      !is.finite(passage_time_weight) || passage_time_weight < 0) {
    stop("passage_time_weight must be one finite non-negative value.")
  }
  if (length(passage_time_tolerance_days) != 1L ||
      !is.finite(passage_time_tolerance_days) || passage_time_tolerance_days < 0) {
    stop("passage_time_tolerance_days must be one finite non-negative value.")
  }
  if (length(passage_time_sigma_days) != 1L ||
      !is.finite(passage_time_sigma_days) || passage_time_sigma_days <= 0) {
    stop("passage_time_sigma_days must be one finite strictly positive value.")
  }
  if (length(passage_time_df) != 1L || !is.finite(passage_time_df) || passage_time_df <= 0) {
    stop("passage_time_df must be one finite strictly positive value.")
  }

  validate_invitro_parameter_table(
    parameter_table = parameter_table,
    dt = dt_use,
    init_total_size = init_total_size_use,
    o2_upper_bound = o2_upper_bound_use,
    fixed_oxygen = fixed_oxygen_use
  )
  validate_invitro_fit_objects(
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path,
    death_data_path = death_data_path
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(parameter_table, file.path(out_dir, "parameter_table_input.csv"), overwrite = TRUE)
  file.copy(death_data_path, file.path(out_dir, "invitro_death_likelihood_input.tsv"), overwrite = TRUE)
  write_invitro_run_provenance(
    out_dir = out_dir,
    argv = argv,
    parameter_table = parameter_table,
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path,
    death_data_path = death_data_path,
    death_weight = death_weight,
    ploidy_weight = ploidy_weight,
    sigma_death_logit = sigma_death_logit,
    death_fraction_eps = death_fraction_eps,
    buffer_prior_weight = buffer_prior_weight,
    buffer_prior_center_smax = buffer_prior_center_smax,
    buffer_prior_sd_smax = buffer_prior_sd_smax,
    buffer_prior_center_log10_beta = buffer_prior_center_log10_beta,
    buffer_prior_sd_log10_beta = buffer_prior_sd_log10_beta,
    buffer_prior_center_log10_n_exp = buffer_prior_center_log10_n_exp,
    buffer_prior_sd_log10_n_exp = buffer_prior_sd_log10_n_exp,
    passage_time_weight = passage_time_weight,
    passage_time_tolerance_days = passage_time_tolerance_days,
    passage_time_sigma_days = passage_time_sigma_days,
    passage_time_df = passage_time_df,
    seed = seed,
    itermax = itermax,
    NP = NP_requested,
    de_reltol = de_reltol,
    de_steptol = de_steptol,
    n_cores = n_cores_requested,
    de_include_parameter_init = de_include_parameter_init
  )
  set.seed(seed)

  cfg_local <- build_invitro_cfg(
    parameter_table = parameter_table,
    dt = dt_use,
    init_total_size = init_total_size_use,
    o2_upper_bound = o2_upper_bound_use,
    fixed_oxygen = fixed_oxygen_use
  )
  fit_objects <- ivt_load_fit_objects_compat(
    fit_objects_dir = fit_objects_dir,
    flow_density_path = flow_density_path,
    death_data_path = death_data_path
  )

  optim_spec <- ivt_optimizer_spec(cfg_local)
  free_spec <- optim_spec[optim_spec$estimate, , drop = FALSE]
  if (nrow(free_spec) == 0L) {
    stop("In vitro parameter table must estimate at least one parameter.")
  }
  free_names <- free_spec$param_name
  init_all <- setNames(as.numeric(optim_spec$init), optim_spec$param_name)
  init_free <- setNames(as.numeric(free_spec$init), free_names)
  lower_free <- setNames(as.numeric(free_spec$lower), free_names)
  upper_free <- setNames(as.numeric(free_spec$upper), free_names)

  objective_from_free <- function(par_free_t) {
    par_t <- init_all
    par_t[free_names] <- as.numeric(par_free_t)
    run_params <- ivt_optim_par_to_run_params(par_t = par_t, cfg = cfg_local)
    comp <- tryCatch(
      ivt_objective_components(
        run_params = run_params,
        fit_objects = fit_objects,
        cfg = cfg_local,
        fallback_max_passage_days = 14,
        growth_weight = 1,
        ploidy_weight = ploidy_weight,
        flow_weight = 1,
        death_weight = death_weight,
        passage_time_weight = passage_time_weight,
        passage_time_tolerance_days = passage_time_tolerance_days,
        passage_time_sigma_days = passage_time_sigma_days,
        passage_time_df = passage_time_df,
        sigma_death_logit = sigma_death_logit,
        death_fraction_eps = death_fraction_eps
      ),
      error = function(e) {
        reason <- if (inherits(e, "invitro_protocol_infeasible")) {
          conditionMessage(e)
        } else {
          paste0("simulation_error: ", conditionMessage(e))
        }
        objective <- if (inherits(e, "invitro_protocol_infeasible")) {
          invitro_protocol_penalty_objective(e)
        } else {
          INVITRO_DE_PENALTY_OBJECTIVE
        }
        make_penalty_components(objective = objective, reason = reason)
      }
    )
    if (is.null(comp$penalty_reason)) {
      buffer_prior <- .ivt_buffer_soft_prior(
        run_params = run_params,
        weight = buffer_prior_weight,
        center_smax = buffer_prior_center_smax,
        sd_smax = buffer_prior_sd_smax,
        center_log10_beta = buffer_prior_center_log10_beta,
        sd_log10_beta = buffer_prior_sd_log10_beta,
        center_log10_n_exp = buffer_prior_center_log10_n_exp,
        sd_log10_n_exp = buffer_prior_sd_log10_n_exp
      )
      comp$likelihood_objective <- comp$objective
      comp$buffer_prior_weight <- buffer_prior$weight
      comp$buffer_prior_raw_penalty <- buffer_prior$raw_penalty
      comp$buffer_prior_penalty <- buffer_prior$weighted_penalty
      comp$buffer_prior_terms <- buffer_prior$terms
      comp$objective <- comp$likelihood_objective + comp$buffer_prior_penalty
    }
    comp$run_params <- run_params
    comp$full_t <- par_t
    comp
  }

  objective_value <- function(par_free_t) {
    objective_from_free(par_free_t)$objective
  }

  cpp_info <- tryCatch(
    o2simps_cpp_dll_info(),
    error = function(e) {
      stop(
        "[fit_invitro] Failed to resolve compiled C++ wrapper/DLL: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  NP_use <- max(NP_requested, 10L * length(free_names))
  de_ctrl <- list(
    trace = TRUE,
    itermax = max(itermax, 1L),
    NP = NP_use,
    reltol = de_reltol,
    steptol = de_steptol
  )
  if (isTRUE(de_include_parameter_init)) {
    de_ctrl$initialpop <- build_invitro_de_initial_population(
      NP = NP_use,
      lower = lower_free,
      upper = upper_free,
      init = init_free
    )
    colnames(de_ctrl$initialpop) <- free_names
    utils::write.table(
      data.frame(
        parameter = free_names,
        transformed_init = as.numeric(init_free),
        transformed_lower = as.numeric(lower_free),
        transformed_upper = as.numeric(upper_free),
        stringsAsFactors = FALSE
      ),
      file.path(out_dir, "de_parameter_init.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message("[fit_invitro] Parameter-table init inserted as DE population member 1.")
  }
  de_cluster <- NULL
  de_active_cores <- 1L
  if (n_cores_requested > 1L) {
    message("[fit_invitro] DEoptim parallel requested with n_cores=", n_cores_requested, ".")
    preflight <- start_invitro_deoptim_cluster_with_preflight(
      n_cores = n_cores_requested,
      objective_from_free = objective_from_free,
      objective_value = objective_value,
      init_free = init_free,
      cpp_info = cpp_info,
      audit_path = file.path(out_dir, "de_worker_preflight.tsv")
    )
    de_cluster <- preflight$cluster
    on.exit(try(parallel::stopCluster(de_cluster), silent = TRUE), add = TRUE)
    parallel::clusterExport(
      de_cluster,
      varlist = c("objective_value"),
      envir = environment()
    )
    de_ctrl$cluster <- de_cluster
    de_active_cores <- preflight$active_cores
    message(
      "[fit_invitro] DEoptim parallel enabled: workers=", de_active_cores,
      "; preflight_retries=", preflight$retries, "."
    )
  } else {
    de_ctrl$parallelType <- "none"
    message("[fit_invitro] DEoptim running in serial mode (n_cores=1).")
  }
  de_fit <- DEoptim::DEoptim(
    fn = objective_value,
    lower = lower_free,
    upper = upper_free,
    control = de_ctrl
  )
  de_best_free_t <- as.numeric(de_fit$optim$bestmem)
  names(de_best_free_t) <- free_names
  de_best_objective <- suppressWarnings(as.numeric(de_fit$optim$bestval))
  best_free_t <- de_best_free_t

  local_maxit <- as_int(.first_non_null_local(argv$local_optim_maxit, argv$optim_maxit, 200L), 200L)
  if (!is.finite(local_maxit) || is.na(local_maxit) || local_maxit < 1L) local_maxit <- 200L
  local_attempted <- FALSE
  local_accepted <- FALSE
  local_fit <- NULL
  local_best_objective <- NA_real_
  local_convergence <- NA_integer_
  local_message <- NA_character_
  if (is.finite(de_best_objective)) {
    local_attempted <- TRUE
    message("[fit_invitro] Starting L-BFGS-B local refinement from DEoptim best; maxit=", local_maxit, ".")
    local_fit <- tryCatch(
      suppressWarnings(
        optim(
          par = best_free_t,
          fn = objective_value,
          method = "L-BFGS-B",
          lower = lower_free,
          upper = upper_free,
          control = list(maxit = local_maxit)
        )
      ),
      error = function(e) {
        warning("[fit_invitro] L-BFGS-B local refinement failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.list(local_fit)) {
      local_best_objective <- suppressWarnings(as.numeric(local_fit$value))
      local_convergence <- suppressWarnings(as.integer(local_fit$convergence))
      local_message <- as.character(.first_non_null_local(local_fit$message, NA_character_))
      if (is.finite(local_best_objective) && local_best_objective < de_best_objective) {
        best_free_t <- as.numeric(local_fit$par)
        names(best_free_t) <- free_names
        best_free_t <- o2sd_clip(best_free_t, lower_free, upper_free)
        local_accepted <- TRUE
        message(
          "[fit_invitro] L-BFGS-B local refinement improved objective: ",
          signif(de_best_objective, 8), " -> ", signif(local_best_objective, 8), "."
        )
      } else {
        message("[fit_invitro] L-BFGS-B local refinement did not improve objective; keeping DEoptim best.")
      }
    }
  }
  de_method <- if (de_active_cores > 1L) "DEoptim_parallel" else "DEoptim_serial"
  optimizer_method <- if (isTRUE(local_attempted)) paste0(de_method, "_plus_LBFGSB_serial") else de_method
  optimizer_trace <- list(
    method = optimizer_method,
    deoptim_objective = as.numeric(de_best_objective),
    local_objective = as.numeric(local_best_objective),
    local_attempted = isTRUE(local_attempted),
    local_accepted = isTRUE(local_accepted),
    local_convergence = as.integer(local_convergence),
    local_maxit = as.integer(local_maxit),
    local_message = local_message
  )

  best_comp <- objective_from_free(best_free_t)
  best_run_params <- best_comp$run_params
  best_full_t <- best_comp$full_t
  best_penalty_reason <- invitro_de_preflight_text(best_comp$penalty_reason)
  if (!is.na(best_penalty_reason) && nzchar(best_penalty_reason)) {
    failure <- data.frame(
      status = if (grepl("^protocol_infeasible:", best_penalty_reason)) {
        "PROTOCOL_INFEASIBLE"
      } else {
        "FIT_PENALTY"
      },
      objective = suppressWarnings(as.numeric(best_comp$objective)),
      reason = best_penalty_reason,
      seed = seed,
      itermax = itermax,
      NP_used = NP_use,
      n_cores_used = de_active_cores,
      stringsAsFactors = FALSE
    )
    utils::write.table(
      failure,
      file = file.path(out_dir, "fit_failure.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )
    stop(
      "[fit_invitro] No feasible best solution; failure recorded in ",
      file.path(out_dir, "fit_failure.tsv"),
      ". Reason: ", best_penalty_reason,
      call. = FALSE
    )
  }

  best_numeric_params <- best_run_params[vapply(best_run_params, is.numeric, logical(1))]
  best_numeric_params <- filter_family_specific_run_params_for_output_common(best_numeric_params)
  best_numeric_params <- best_numeric_params[!vapply(best_numeric_params, is.null, logical(1))]

  best_params_df <- data.frame(
    parameter = names(best_numeric_params),
    value = as.numeric(unlist(best_numeric_params)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  write.table(best_params_df, file = file.path(out_dir, "best_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  best_transformed_df <- data.frame(
    transformed_parameter = names(best_full_t),
    transformed_value = as.numeric(best_full_t),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  write.table(best_transformed_df, file = file.path(out_dir, "best_params_transformed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  postfit_tables <- ivt_collect_postfit_tables(best_comp)
  for (table_name in names(postfit_tables)) {
    write_tsv_if_nonempty(
      postfit_tables[[table_name]],
      file.path(out_dir, paste0(table_name, ".tsv"))
    )
  }
  write_tsv_if_nonempty(best_comp$ploidy_df, file.path(out_dir, "invitro_ploidy_loglik.tsv"))
  write_tsv_if_nonempty(best_comp$flow_df, file.path(out_dir, "invitro_flow_loglik.tsv"))
  write_tsv_if_nonempty(best_comp$flow_overlay_df, file.path(out_dir, "invitro_flow_overlay.tsv"))
  write_tsv_if_nonempty(best_comp$objective_hierarchy, file.path(out_dir, "invitro_objective_hierarchy.tsv"))
  write_tsv_if_nonempty(best_comp$passage_time_df, file.path(out_dir, "invitro_passage_time_loglik.tsv"))
  write_tsv_if_nonempty(best_comp$buffer_prior_terms, file.path(out_dir, "invitro_buffer_prior_terms.tsv"))

  dist_summary <- dplyr::bind_rows(
    ivt_collect_distribution_summary(best_comp$run_2N),
    ivt_collect_distribution_summary(best_comp$run_4N)
  )
  write_tsv_if_nonempty(dist_summary, file.path(out_dir, "invitro_distribution_summary.tsv"))

  ploidy_quantile_probs <- seq(0.01, 0.99, length.out = 50L)
  dist_quantiles <- dplyr::bind_rows(
    ivt_collect_distribution_quantiles(best_comp$run_2N, probs = ploidy_quantile_probs),
    ivt_collect_distribution_quantiles(best_comp$run_4N, probs = ploidy_quantile_probs)
  )
  write_tsv_if_nonempty(dist_quantiles, file.path(out_dir, "invitro_distribution_quantiles.tsv"))

  observed_kary <- dplyr::bind_rows(
    ivt_collect_observed_kary_summary(best_comp$run_2N, fit_objects$fit_data),
    ivt_collect_observed_kary_summary(best_comp$run_4N, fit_objects$fit_data)
  )
  write_tsv_if_nonempty(observed_kary, file.path(out_dir, "invitro_observed_kary.tsv"))

  observed_flow <- dplyr::bind_rows(
    ivt_collect_observed_flow_summary(best_comp$run_2N, fit_objects$fit_data),
    ivt_collect_observed_flow_summary(best_comp$run_4N, fit_objects$fit_data)
  )
  write_tsv_if_nonempty(observed_flow, file.path(out_dir, "invitro_observed_flow.tsv"))

  summary_df <- data.frame(
    metric = c(
      "fit_mode",
      "optimizer_method",
      "optimizer_deoptim_objective",
      "optimizer_local_objective",
      "optimizer_local_attempted",
      "optimizer_local_accepted",
      "optimizer_local_convergence",
      "optimizer_local_maxit",
      "objective_total",
      "objective_likelihood",
      "total_loglik",
      "growth_loglik",
      "ploidy_loglik",
      "flow_loglik",
      "death_loglik",
      "passage_time_loglik",
      "growth_weight",
      "ploidy_weight",
      "flow_weight",
      "death_weight",
      "passage_time_weight",
      "buffer_prior_weight",
      "buffer_prior_raw_penalty",
      "buffer_prior_penalty",
      "buffer_prior_center_smax",
      "buffer_prior_sd_smax",
      "buffer_prior_center_log10_beta",
      "buffer_prior_sd_log10_beta",
      "buffer_prior_center_log10_n_exp",
      "buffer_prior_sd_log10_n_exp",
      "growth_loglik_sum",
      "growth_passage_loglik_sum",
      "ploidy_loglik_sum",
      "flow_loglik_sum",
      "death_loglik_sum",
      "passage_time_loglik_sum",
      "sigma_growth",
      "sigma_kary",
      "sigma_flow_ploidy",
      "sigma_death_logit",
      "death_fraction_eps",
      "passage_time_tolerance_days",
      "passage_time_sigma_days",
      "passage_time_df",
      "n_growth",
      "n_growth_timepoints",
      "n_growth_passages",
      "n_growth_observed",
      "n_growth_missing_pred",
      "n_growth_negative_pred",
      "n_ploidy_passages",
      "n_kary_cells",
      "n_flow_passages",
      "n_flow_samples",
      "n_death_passages",
      "n_passage_time_observations",
      "n_scenarios",
      "n_insufficient_boundaries",
      "all_passage_boundaries_feasible",
      "protocol_feasibility_status",
      "max_boundary_scale",
      "seed",
      "itermax",
      "itermax_requested",
      "itermax_max",
      "de_reltol",
      "de_steptol",
      "de_include_parameter_init",
      "NP_requested",
      "NP_used",
      "n_cores_requested",
      "n_cores_used",
      "dt",
      "init_total_size",
      "parameter_table",
      "fit_objects_dir",
      "flow_density_path",
      "death_data_path"
    ),
    value = c(
      "fit_invitro",
      optimizer_method,
      as.character(de_best_objective),
      as.character(local_best_objective),
      as.character(local_attempted),
      as.character(local_accepted),
      as.character(local_convergence),
      as.character(local_maxit),
      as.character(best_comp$objective),
      as.character(best_comp$likelihood_objective),
      as.character(best_comp$total_loglik),
      as.character(best_comp$growth_loglik),
      as.character(best_comp$ploidy_loglik),
      as.character(best_comp$flow_loglik),
      as.character(best_comp$death_loglik),
      as.character(best_comp$passage_time_loglik),
      "1",
      as.character(ploidy_weight),
      "1",
      as.character(death_weight),
      as.character(passage_time_weight),
      as.character(best_comp$buffer_prior_weight),
      as.character(best_comp$buffer_prior_raw_penalty),
      as.character(best_comp$buffer_prior_penalty),
      as.character(buffer_prior_center_smax),
      as.character(buffer_prior_sd_smax),
      as.character(buffer_prior_center_log10_beta),
      as.character(buffer_prior_sd_log10_beta),
      as.character(buffer_prior_center_log10_n_exp),
      as.character(buffer_prior_sd_log10_n_exp),
      as.character(best_comp$growth_loglik_sum),
      as.character(best_comp$growth_passage_loglik_sum),
      as.character(best_comp$ploidy_loglik_sum),
      as.character(best_comp$flow_loglik_sum),
      as.character(best_comp$death_loglik_sum),
      as.character(best_comp$passage_time_loglik_sum),
      as.character(best_comp$sigma_growth),
      as.character(best_comp$sigma_kary),
      as.character(best_comp$sigma_flow_ploidy),
      as.character(best_comp$sigma_death_logit),
      as.character(best_comp$death_fraction_eps),
      as.character(passage_time_tolerance_days),
      as.character(passage_time_sigma_days),
      as.character(passage_time_df),
      as.character(best_comp$n_growth),
      as.character(best_comp$n_growth_timepoints),
      as.character(best_comp$n_growth_passages),
      as.character(best_comp$n_growth_observed),
      as.character(best_comp$n_growth_missing_pred),
      as.character(best_comp$n_growth_negative_pred),
      as.character(best_comp$n_ploidy_passages),
      as.character(best_comp$n_kary_cells),
      as.character(best_comp$n_flow_passages),
      as.character(best_comp$n_flow_samples),
      as.character(best_comp$n_death_passages),
      as.character(best_comp$n_passage_time_observations),
      as.character(best_comp$n_scenarios),
      as.character(best_comp$n_insufficient_boundaries),
      as.character(best_comp$all_passage_boundaries_feasible),
      as.character(best_comp$protocol_feasibility_status),
      as.character(best_comp$max_boundary_scale),
      as.character(seed),
      as.character(itermax),
      as.character(itermax_requested),
      as.character(itermax_max),
      as.character(de_reltol),
      as.character(de_steptol),
      as.character(de_include_parameter_init),
      as.character(NP_requested),
      as.character(NP_use),
      as.character(n_cores_requested),
      as.character(de_active_cores),
      as.character(dt_use),
      as.character(init_total_size_use),
      normalizePath(parameter_table, mustWork = FALSE),
      normalizePath(fit_objects_dir, mustWork = FALSE),
      normalizePath(flow_density_path, mustWork = FALSE),
      normalizePath(death_data_path, mustWork = FALSE)
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  summary_df <- filter_fit_summary_metrics_for_output_common(summary_df)
  write.table(summary_df, file = file.path(out_dir, "fit_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  saveRDS(
    list(
      cfg = cfg_local,
      deoptim = de_fit,
      local_optim = local_fit,
      optimizer_trace = optimizer_trace,
      best_components = best_comp,
      best_params = best_run_params,
      fit_objects_dir = normalizePath(fit_objects_dir, mustWork = FALSE),
      flow_density_path = normalizePath(flow_density_path, mustWork = FALSE),
      death_data_path = normalizePath(death_data_path, mustWork = FALSE),
      death_weight = death_weight,
      ploidy_weight = ploidy_weight,
      sigma_death_logit = sigma_death_logit,
      death_fraction_eps = death_fraction_eps,
      buffer_prior_weight = buffer_prior_weight,
      buffer_prior_center_smax = buffer_prior_center_smax,
      buffer_prior_sd_smax = buffer_prior_sd_smax,
      buffer_prior_center_log10_beta = buffer_prior_center_log10_beta,
      buffer_prior_sd_log10_beta = buffer_prior_sd_log10_beta,
      buffer_prior_center_log10_n_exp = buffer_prior_center_log10_n_exp,
      buffer_prior_sd_log10_n_exp = buffer_prior_sd_log10_n_exp,
      de_include_parameter_init = de_include_parameter_init,
      passage_time_weight = passage_time_weight,
      passage_time_tolerance_days = passage_time_tolerance_days,
      passage_time_sigma_days = passage_time_sigma_days,
      passage_time_df = passage_time_df
    ),
    file = file.path(out_dir, "fit_result.rds")
  )

  if (isTRUE(auto_viz)) {
    run_invitro_auto_viz_report(out_dir)
  }

  message("Done. Results written to: ", normalizePath(out_dir, mustWork = FALSE))
  invisible(normalizePath(out_dir, mustWork = FALSE))
}

if (sys.nframe() == 0) {
  main()
}
