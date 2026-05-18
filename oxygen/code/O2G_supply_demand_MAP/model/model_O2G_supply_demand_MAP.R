suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

# ----------------------------------------------------------------------------
# Align miningcloneid oxygen model to Richard's buffering.R logic.
# model_O2G_supply_demand_MAP extension:
# - Keep karyotype dynamics identical to Richard-aligned model.
# - Use a logarithmic oxygen supply-demand target with floor:
#   O2_target = max(o2_min, o2_S0 - kappa_O * log(1 + N_eff / o2_Nref)).
# - Self-contained model file (no dependency on model_functions_ploidy_buffer.R)
# - Core dynamics (lambda, misseg delta, B/G construction) are defined here
# - Enforce defaults requested:
#   * WGD branch offspring weight = +1 (not +2)
#   * boundary default = "drop"
# ----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Function: .resolve_align_script_dir
# Purpose: Resolve script directory path for robust relative file loading.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.resolve_align_script_dir <- function() {
  env_dir <- Sys.getenv("MININGCLONEID_OXYGEN_CODE_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = FALSE))
  }

  frs <- sys.frames()
  for (i in rev(seq_along(frs))) {
    ofile <- frs[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
  }

  args <- commandArgs(trailingOnly = FALSE)
  farg <- args[grepl("^--file=", args)]
  if (length(farg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", farg[[1]]), mustWork = FALSE)))
  }

  getwd()
}

.ALIGN_MODEL_DIR <- .resolve_align_script_dir()
.ALIGN_WORKFLOW_ROOT <- normalizePath(file.path(.ALIGN_MODEL_DIR, ".."), mustWork = FALSE)
source(file.path(.ALIGN_WORKFLOW_ROOT, "util", "o2g_supply_demand_map_common_semantics.R"), local = FALSE)

.init_cpp_o2simps_backend <- local({
  initialized <- FALSE
  available <- FALSE

# -----------------------------------------------------------------------------
# Function: acquire_dir_lock
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - lock_dir: Function-specific input argument.
#   - timeout_sec: Function-specific input argument.
#   - poll_sec: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
  acquire_dir_lock <- function(lock_dir, timeout_sec = 300, poll_sec = 0.1) {
    start <- Sys.time()
    repeat {
      if (dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)) {
        return(TRUE)
      }

      # Recover from stale lock directories left by crashed workers.
      if (dir.exists(lock_dir)) {
        info <- suppressWarnings(file.info(lock_dir))
        mtime <- info$mtime[[1]]
        if (is.finite(as.numeric(mtime))) {
          age <- as.numeric(difftime(Sys.time(), mtime, units = "secs"))
          if (is.finite(age) && age > timeout_sec) {
            unlink(lock_dir, recursive = TRUE, force = TRUE)
            Sys.sleep(0.02)
            next
          }
        }
      }

      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      if (is.finite(elapsed) && elapsed >= timeout_sec) return(FALSE)
      Sys.sleep(poll_sec + stats::runif(1L, min = 0, max = 0.05))
    }
  }

  function() {
    if (initialized) return(available)
    initialized <<- TRUE

    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      stop("Rcpp package is required for model_O2G_supply_demand_MAP but is not installed.")
    }

    # Dedicated O2 invivo backend (do not use Richard/shared cpp here).
    cpp_path <- file.path(.ALIGN_MODEL_DIR, "model_O2G_supply_demand_MAP.cpp")
    if (!file.exists(cpp_path)) {
      stop("Cannot find required C++ backend file: ", cpp_path)
    }

    cache_root <- file.path(.ALIGN_MODEL_DIR, ".rcpp_cache_o2g_supply_demand_map")
    cache_dir <- file.path(cache_root, "shared")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    rebuild_cpp <- tolower(trimws(Sys.getenv("MININGCLONEID_RCPP_REBUILD", unset = "FALSE"))) %in% c("1", "true", "t", "yes", "y")
    lock_timeout_sec <- suppressWarnings(as.numeric(Sys.getenv("MININGCLONEID_RCPP_LOCK_TIMEOUT_SEC", unset = "300")))
    if (!is.finite(lock_timeout_sec) || lock_timeout_sec <= 0) lock_timeout_sec <- 300

    lock_dir <- file.path(cache_root, ".sourcecpp_lock")
    lock_ok <- acquire_dir_lock(lock_dir, timeout_sec = lock_timeout_sec, poll_sec = 0.1)
    if (!isTRUE(lock_ok)) {
      stop("Timed out waiting for sourceCpp lock: ", lock_dir)
    }
    on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)

    tryCatch({
      Rcpp::sourceCpp(
        file = cpp_path,
        rebuild = rebuild_cpp,
        showOutput = FALSE,
        verbose = FALSE,
        cacheDir = cache_dir
      )
    }, error = function(e) {
      stop("Failed to compile/load model_O2G_supply_demand_MAP.cpp: ", conditionMessage(e))
    })

    required_fns <- c(
      "cpp_o2simps_pr_delta_vec",
      "cpp_o2simps_build_B_total_triplet",
      "cpp_o2simps_build_B_WGD_triplet",
      "cpp_o2simps_o2_window_supply",
      "cpp_o2simps_build_G_for_o2_triplet",
      "cpp_o2simps_simulate_one",
      "cpp_o2simps_objective_components_map"
    )
    missing_fns <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function", inherits = TRUE)]
    if (length(missing_fns) > 0L) {
      stop(
        "model_O2G_supply_demand_MAP C++ backend loaded but required symbols are missing: ",
        paste(missing_fns, collapse = ", ")
      )
    }

    wrappers_need_rebuild <- FALSE
    wrapper_mismatch_reason <- character(0)
    check_wrapper_formals <- function(fn_name, must_have = character(0), must_absent = character(0)) {
      if (!exists(fn_name, mode = "function", inherits = TRUE)) {
        return(FALSE)
      }
      f <- get(fn_name, mode = "function", inherits = TRUE)
      nms <- names(formals(f))
      miss <- setdiff(must_have, nms)
      bad <- intersect(must_absent, nms)
      if (length(miss) > 0L || length(bad) > 0L) {
        wrappers_need_rebuild <<- TRUE
        wrapper_mismatch_reason <<- c(
          wrapper_mismatch_reason,
          paste0(
            fn_name, " formal mismatch",
            if (length(miss) > 0L) paste0(" missing{", paste(miss, collapse = ","), "}") else "",
            if (length(bad) > 0L) paste0(" forbidden{", paste(bad, collapse = ","), "}") else ""
          )
        )
      }
      TRUE
    }
    check_wrapper_formals(
      "cpp_o2simps_build_G_for_o2_triplet",
      must_have = c("O2_crit", "glucose", "p_wgd", "p_wgd_max", "O2_wgd", "O2_growth", "n_O", "ploidy_O2_death"),
      must_absent = c("o2_ref_pct")
    )
    check_wrapper_formals(
      "cpp_o2simps_simulate_one",
      must_have = c("sim_args"),
      must_absent = c("o2_ref_pct")
    )
    check_wrapper_formals(
      "cpp_o2simps_objective_components_map",
      must_have = c("scenario_data", "objective_data", "state_data", "sim_args"),
      must_absent = c("o2_ref_pct", "O2_growth", "crowding_enabled", "burden_log_eps")
    )

    if (isTRUE(wrappers_need_rebuild) && !isTRUE(rebuild_cpp)) {
      # Stale sourceCpp wrapper cache can keep outdated formals; force rebuild once.
      tryCatch({
        Rcpp::sourceCpp(
          file = cpp_path,
          rebuild = TRUE,
          showOutput = FALSE,
          verbose = FALSE,
          cacheDir = cache_dir
        )
      }, error = function(e) {
        stop(
          "Failed forced rebuild for model_O2G_supply_demand_MAP.cpp after wrapper mismatch [",
          paste(wrapper_mismatch_reason, collapse = "; "),
          "]: ", conditionMessage(e)
        )
      })

      wrapper_mismatch_reason <- character(0)
      wrappers_need_rebuild <- FALSE
      check_wrapper_formals(
        "cpp_o2simps_build_G_for_o2_triplet",
        must_have = c("O2_crit", "glucose", "p_wgd", "p_wgd_max", "O2_wgd", "O2_growth", "n_O", "ploidy_O2_death"),
        must_absent = c("o2_ref_pct")
      )
      check_wrapper_formals(
        "cpp_o2simps_simulate_one",
        must_have = c("sim_args"),
        must_absent = c("o2_ref_pct")
      )
      check_wrapper_formals(
        "cpp_o2simps_objective_components_map",
        must_have = c("scenario_data", "objective_data", "state_data", "sim_args"),
        must_absent = c("o2_ref_pct", "O2_growth", "crowding_enabled", "burden_log_eps")
      )
      if (isTRUE(wrappers_need_rebuild)) {
        stop(
          "model_O2G_supply_demand_MAP wrapper signatures are inconsistent after forced rebuild: ",
          paste(wrapper_mismatch_reason, collapse = "; ")
        )
      }
    }

    available <<- TRUE

    available
  }
})

.USE_CPP_O2SIMPS_BACKEND <- .init_cpp_o2simps_backend()

# -----------------------------------------------------------------------------
# Function: o2simps_cpp_dll_info
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
o2simps_cpp_dll_info <- function() {
  if (!isTRUE(.USE_CPP_O2SIMPS_BACKEND)) {
    stop("C++ backend is not initialized.")
  }
  if (!exists("cpp_o2simps_simulate_one", mode = "function", inherits = TRUE)) {
    stop("Required wrapper function missing: cpp_o2simps_simulate_one")
  }

  loaded <- getLoadedDLLs()
  if (length(loaded) == 0L) stop("No loaded DLLs found after sourceCpp initialization.")
  dll_names <- names(loaded)
  dll_paths <- vapply(loaded, function(x) as.character(x[["path"]]), FUN.VALUE = character(1), USE.NAMES = FALSE)
  valid <- nzchar(dll_paths) & file.exists(dll_paths)

  # Prefer DLLs from this model's dedicated sourceCpp cache.
  cache_pat <- ".rcpp_cache_o2g_supply_demand_map"
  in_cache <- valid & grepl(cache_pat, dll_paths, fixed = TRUE)
  candidate <- if (any(in_cache)) in_cache else (valid & grepl("sourceCpp", dll_names, fixed = TRUE))
  if (!any(candidate)) {
    stop("Unable to resolve sourceCpp DLL path from loaded DLL list.")
  }

  cand_idx <- which(candidate)
  mt <- suppressWarnings(file.info(dll_paths[cand_idx])$mtime)
  best <- cand_idx[[1]]
  if (length(cand_idx) > 1L && any(is.finite(as.numeric(mt)))) {
    ord <- order(mt, decreasing = TRUE, na.last = TRUE)
    best <- cand_idx[[ord[[1]]]]
  }
  dll_path <- normalizePath(dll_paths[[best]], mustWork = TRUE)
  wrapper_candidates <- list.files(
    dirname(dll_path),
    pattern = "\\.cpp\\.R$",
    full.names = TRUE
  )
  wrapper_path <- ""
  if (length(wrapper_candidates) > 0L) {
    if (length(wrapper_candidates) == 1L) {
      wrapper_path <- wrapper_candidates[[1]]
    } else {
      hit <- vapply(
        wrapper_candidates,
        function(p) {
          txt <- tryCatch(readLines(p, warn = FALSE), error = function(e) character(0))
          any(grepl("cpp_o2simps_simulate_one", txt, fixed = TRUE))
        },
        logical(1)
      )
      if (any(hit)) wrapper_path <- wrapper_candidates[which(hit)[1]]
    }
  }
  if (!nzchar(wrapper_path) || !file.exists(wrapper_path)) {
    stop("Unable to resolve sourceCpp wrapper file (*.cpp.R) for O2G_supply_demand_MAP backend.")
  }

  list(
    name = as.character(dll_names[[best]]),
    path = dll_path,
    wrapper_path = normalizePath(wrapper_path, mustWork = TRUE)
  )
}

# -----------------------------------------------------------------------------
# Function: .first_non_null
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - ...: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.first_non_null <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v)) return(v)
  }
  NULL
}

# Euler stepper for generator matrix dynamics.
# -----------------------------------------------------------------------------
# Function: step_dt
# Purpose: Advance state vector by repeated generator-matrix time steps.
# Parameters:
#   - G: Generator or transition matrix used for state propagation.
#   - x: Input value or vector to process.
#   - dt: Time-step length in days.
#   - steps: Number of repeated matrix-step applications.
#   - normalize: Logical flag to normalize resulting state distribution.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
step_dt <- function(G, x, dt, steps = 1L, normalize = FALSE) {
  I <- Diagonal(n = nrow(G))
  A <- I + dt * G
  v <- as.numeric(x)
  for (k in seq_len(steps)) {
    v <- as.numeric(A %*% v)
    if (normalize) v <- v / sum(v)
  }
  v
}

# Chromosome-length-weighted ploidy mapping on fixed autosomes 1..22.
# -----------------------------------------------------------------------------
# Function: default_chr_lengths_bp_1to22
# Purpose: Return fixed chromosome lengths (bp) for autosomes 1..22.
# Parameters:
#   - (none): This helper consumes surrounding scope or global options.
# Returns:
#   Numeric vector of chromosome lengths in base pairs.
# -----------------------------------------------------------------------------
default_chr_lengths_bp_1to22 <- function() {
  c(
    `1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555,
    `5` = 181538259, `6` = 170805979, `7` = 159345973, `8` = 145138636,
    `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
    `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345,
    `17` = 83257441, `18` = 80373285, `19` = 58617616, `20` = 64444167,
    `21` = 46709983, `22` = 50818468
  )
}

# Precompute immutable defaults once; these are reused throughout the fit/eval path.
.chr_lengths_default_bp_1to22 <- default_chr_lengths_bp_1to22()
.chr_lengths_default_ord_desc <- order(.chr_lengths_default_bp_1to22, decreasing = TRUE)
.chr_lengths_default_denom <- sum(.chr_lengths_default_bp_1to22)

# -----------------------------------------------------------------------------
# Function: normalize_chr_lengths_bp_1to22
# Purpose: Validate and normalize chromosome-length vector for autosomes 1..22.
# Parameters:
#   - chr_lengths_bp: Optional chromosome-length vector.
# Returns:
#   Numeric vector of length 22 with positive finite values.
# -----------------------------------------------------------------------------
normalize_chr_lengths_bp_1to22 <- function(chr_lengths_bp = NULL) {
  w <- if (is.null(chr_lengths_bp)) default_chr_lengths_bp_1to22() else as.numeric(chr_lengths_bp)
  if (length(w) != 22L) stop("chr_lengths_bp must have length 22 (autosomes 1..22).")
  if (any(!is.finite(w) | w <= 0)) stop("chr_lengths_bp must be all positive finite values.")
  names(w) <- as.character(seq_len(22L))
  w
}

# -----------------------------------------------------------------------------
# Function: weighted_ploidy_from_total_N
# Purpose: Map total copy-number state N to chromosome-length-weighted ploidy.
# Parameters:
#   - N_total: Total copy number state(s) on the integer grid.
#   - chr_lengths_bp: Optional chromosome-length vector.
# Returns:
#   Numeric weighted ploidy values.
# -----------------------------------------------------------------------------
weighted_ploidy_from_total_N <- function(N_total, chr_lengths_bp = NULL) {
  if (is.null(chr_lengths_bp)) {
    w <- .chr_lengths_default_bp_1to22
    ord <- .chr_lengths_default_ord_desc
    n_chr <- 22L
    denom <- .chr_lengths_default_denom
  } else {
    w <- normalize_chr_lengths_bp_1to22(chr_lengths_bp)
    ord <- order(w, decreasing = TRUE)
    n_chr <- length(w)
    denom <- sum(w)
  }
  Nv <- as.numeric(N_total)
  vapply(Nv, function(nn) {
    if (!is.finite(nn)) return(NA_real_)
    n_int <- as.integer(round(nn))
    if (n_int < 0L) n_int <- 0L
    base <- n_int %/% n_chr
    rem <- n_int %% n_chr
    cn <- rep.int(base, n_chr)
    if (rem > 0L) cn[ord[seq_len(rem)]] <- cn[ord[seq_len(rem)]] + 1L
    sum(cn * w) / denom
  }, numeric(1))
}

# -----------------------------------------------------------------------------
# Function: map_ploidy_to_N_by_chrlen
# Purpose: Map ploidy value(s) to nearest integer N state under weighted mapping.
# Parameters:
#   - ploidy_values: Observed ploidy value(s).
#   - N_grid: Integer N grid.
#   - chr_lengths_bp: Optional chromosome-length vector.
# Returns:
#   Integer N states on N_grid.
# -----------------------------------------------------------------------------
map_ploidy_to_N_by_chrlen <- function(ploidy_values, N_grid, chr_lengths_bp = NULL) {
  grid <- as.integer(sort(unique(N_grid)))
  p_grid <- weighted_ploidy_from_total_N(grid, chr_lengths_bp = chr_lengths_bp)
  pv <- as.numeric(ploidy_values)
  vapply(pv, function(p) {
    if (!is.finite(p)) return(NA_integer_)
    d <- abs(p_grid - p)
    k <- which.min(d)
    as.integer(grid[[k]])
  }, integer(1))
}

# Build a normalized initial distribution on an integer N grid from ploidy values.
# -----------------------------------------------------------------------------
# Function: create_initial_dist
# Purpose: Construct initial ploidy-state distribution/state vector for simulation.
# Parameters:
#   - ploidy_values: Function-specific input argument.
#   - N_grid: Function-specific input argument.
#   - N_unit: Number of modeled chromosome classes.
#   - chr_lengths_bp: Optional chromosome-length vector for weighted ploidy mapping.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
create_initial_dist <- function(ploidy_values, N_grid, N_unit = 22L, chr_lengths_bp = NULL) {
  N_values <- map_ploidy_to_N_by_chrlen(
    ploidy_values = ploidy_values,
    N_grid = N_grid,
    chr_lengths_bp = chr_lengths_bp
  )
  N_counts <- table(N_values)
  N_fracs <- as.numeric(N_counts) / sum(N_counts)
  names(N_fracs) <- names(N_counts)
  x_vec <- rep(0, length(N_grid))
  names(x_vec) <- N_grid
  valid_names <- names(N_fracs)[names(N_fracs) %in% names(x_vec)]
  x_vec[valid_names] <- N_fracs[valid_names]
  x_vec
}

# Build explicit initial state vectors on a single chromosome-count grid.
# - Either pass ploidy or a custom initial fraction vector.
# - Returns absolute cell counts with total mass = total_size.
# -----------------------------------------------------------------------------
# Function: make_init_state
# Purpose: Construct initial ploidy-state distribution/state vector for simulation.
# Parameters:
#   - grid_pre: Chromosome-count grid.
#   - ploidy: Initial ploidy label used to choose nearest N state.
#   - init_dist: Optional named initial fractions on grid_pre.
#   - N_UNIT: Function-specific input argument.
#   - total_size: Function-specific input argument.
#   - chr_lengths_bp: Optional chromosome-length vector for weighted ploidy mapping.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
make_init_state <- function(grid_pre,
                            ploidy = c(2, 4),
                            init_dist = NULL,
                            N_UNIT = 22L, total_size = 1e6, chr_lengths_bp = NULL) {
  # Single-layer initialization: all cohorts are initialized on one N grid.
  ploidy <- match.arg(as.character(ploidy), choices = c("2", "4"))
  Pnum <- as.numeric(ploidy)
  grid_N <- as.integer(grid_pre)
  x <- rep(0, length(grid_N))
  names(x) <- grid_N

  if (!is.null(init_dist)) {
    nm <- intersect(names(init_dist), names(x))
    x[nm] <- x[nm] + as.numeric(init_dist[nm])
  } else {
    N_delta <- as.integer(map_ploidy_to_N_by_chrlen(
      ploidy_values = Pnum,
      N_grid = grid_N,
      chr_lengths_bp = chr_lengths_bp
    ))
    stopifnot(N_delta %in% grid_N)
    x[as.character(N_delta)] <- 1
  }

  s <- sum(x)
  if (s <= 0) stop("Init mass is zero.")
  x / s * total_size
}

# -----------------------------------------------------------------------------
# Function: .clip01
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.clip01 <- function(x) pmin(pmax(x, 0), 1)
# -----------------------------------------------------------------------------
# Function: .clip_o2pct
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.clip_o2pct <- function(x) pmin(pmax(x, 0), 100)
# -----------------------------------------------------------------------------
# Function: .assert_o2_pct
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - x: Input value or vector to process.
#   - label: Text label used for logging and progress messages.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.assert_o2_pct <- function(x, label = "O2") {
  x_num <- as.numeric(x)
  if (any(!is.finite(x_num))) stop(label, " must be finite.")
  if (any(x_num < 0 | x_num > 100)) stop(label, " must be in percent scale [0, 100].")
  x_num
}

# -----------------------------------------------------------------------------
# Function: .require_cpp_o2simps_fn
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - fn_name: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.require_cpp_o2simps_fn <- function(fn_name) {
  if (!(isTRUE(.USE_CPP_O2SIMPS_BACKEND) &&
        exists(fn_name, mode = "function", inherits = TRUE))) {
    stop("Required C++ backend function is unavailable: ", fn_name)
  }
}

# O2 target from a logarithmic supply-demand model on effective demand:
#   O2_target = max(o2_min, o2_S0 - kappa_O * log(1 + N_eff / o2_Nref)),
# where N_eff can encode ploidy-weighted demand.
# -----------------------------------------------------------------------------
# Function: .o2_supply_demand_from_burden
# Purpose: Compute oxygen target from viable burden under the supply-demand model.
# Parameters:
#   - Ntot: Effective oxygen-demand proxy at current time.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - o2_Nref: Fixed viable-cell scaling constant for demand normalization.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.o2_supply_demand_from_burden <- function(Ntot,
                                          run_params,
                                          o2_Nref = 1e6) {
  Ntot_use <- pmax(as.numeric(Ntot), 0)
  o2_s0_upper_use <- as.numeric(.first_non_null(run_params$o2_S0_upper_bound, run_params$o2_S0_max, 5.0))
  if (!is.finite(o2_s0_upper_use) || o2_s0_upper_use <= 0) o2_s0_upper_use <- 5.0
  o2_S0 <- as.numeric(.first_non_null(run_params$o2_S0, 0.5))
  kappa_O <- as.numeric(.first_non_null(run_params$kappa_O, 1.0))
  o2_min <- as.numeric(.first_non_null(run_params$o2_min, 0.5))
  if (!is.finite(kappa_O) || kappa_O <= 0) kappa_O <- 1.0
  if (!is.finite(o2_min) || o2_min < 0) o2_min <- 0.5
  if (!is.finite(o2_Nref) || o2_Nref <= 0) o2_Nref <- 1e6
  o2_S0 <- max(0, min(o2_s0_upper_use, o2_S0))
  o2_min <- max(0, min(o2_s0_upper_use, o2_min))

  .require_cpp_o2simps_fn("cpp_o2simps_o2_window_supply")
  return(as.numeric(cpp_o2simps_o2_window_supply(
    Ntot = as.numeric(Ntot_use),
    o2_S0 = as.numeric(o2_S0),
    kappa_O = as.numeric(kappa_O),
    o2_Nref = as.numeric(o2_Nref),
    o2_min = as.numeric(o2_min)
  )))
}

# Richard buffering.R-style lambda: O2-only saturating rate.
# -----------------------------------------------------------------------------
# Function: growth_lambda
# Purpose: Compute oxygen-dependent proliferation rate for a given ploidy state.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - N: Ploidy state value or chromosome-copy count.
#   - lam_min: Lower asymptote of proliferation rate.
#   - lam_max: Upper asymptote of proliferation rate.
#   - k_o: Oxygen-sensitivity parameter for proliferation rate.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
growth_lambda <- function(O2, N, lam_min, lam_max, k_o) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  lam_min_use <- as.numeric(lam_min)
  lam_max_use <- as.numeric(lam_max)
  k_o_use <- as.numeric(k_o)
  if (!is.finite(lam_min_use)) stop("lam_min must be finite.")
  if (!is.finite(lam_max_use)) stop("lam_max must be finite.")
  if (!is.finite(k_o_use) || k_o_use <= 0) stop("k_o must be > 0.")
  k_o_use <- max(k_o_use, 1e-12)

  frac <- O2_use / (O2_use + k_o_use)
  lam <- lam_min_use + (lam_max_use - lam_min_use) * frac
  rep(pmax(lam, 0), length(N))
}

# Main-path baseline-plus-increment missegregation helper (aligned with C++).
# -----------------------------------------------------------------------------
# Function: .mu_eff_of_O2
# Purpose: Compute state-specific hypoxia death rate under the main model.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - N: Ploidy state value or chromosome-copy count.
#   - O2_crit: Optional hypoxia Hill critical scale override.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.resource_stress_of_O2 <- function(O2, run_params, O2_crit = NULL, G = NULL) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  G_use <- if (is.null(G)) O2_use else .assert_o2_pct(G, label = "G")
  n_out <- max(length(O2_use), length(G_use))
  if (!(length(O2_use) %in% c(1L, n_out) && length(G_use) %in% c(1L, n_out))) {
    stop("O2 and G must have compatible lengths.")
  }
  O2_vec <- rep_len(as.numeric(O2_use), n_out)
  G_vec <- rep_len(as.numeric(G_use), n_out)
  O2_crit_use <- as.numeric(.first_non_null(O2_crit, run_params$O2_crit, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  n_O <- as.numeric(.first_non_null(run_params$n_O, 1.0))
  if (!is.finite(n_O) || n_O < 0) {
    stop("run_params$n_O must be finite and >= 0.")
  }
  o2_c <- pmax(O2_crit_use, 1e-12)
  h_o2 <- (o2_c^n_O) / ((o2_c^n_O) + (pmax(O2_vec, 0)^n_O))
  h_o2 <- .clip01(h_o2)
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null(run_params$glucose, TRUE),
    default = TRUE
  ))
  if (!isTRUE(glucose_use)) {
    return(h_o2)
  }
  h_g <- (o2_c^n_O) / ((o2_c^n_O) + (pmax(G_vec, 0)^n_O))
  h_g <- .clip01(h_g)
  .clip01(1 - (1 - h_o2) * (1 - h_g))
}

# Main-path proliferation helper aligned with the current C++ runtime dispatch.
# -----------------------------------------------------------------------------
# Function: .lambda_eff_of_O2
# Purpose: Compute state-specific effective proliferation under the active model.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - N: Ploidy state value or chromosome-copy count.
#   - O2_crit: Optional hypoxia Hill critical scale override.
#   - O2_growth: Optional runtime override for the growth-penalty switch.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.lambda_eff_of_O2 <- function(O2, run_params, N = 44, O2_crit = NULL, O2_growth = TRUE, G = NULL) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  N_use <- as.numeric(N)
  if (any(!is.finite(N_use))) stop("N must be finite.")
  G_use <- if (is.null(G)) O2_use else .assert_o2_pct(G, label = "G")
  n_out <- max(length(O2_use), length(N_use), length(G_use))
  if (!(length(O2_use) %in% c(1L, n_out) &&
        length(N_use) %in% c(1L, n_out) &&
        length(G_use) %in% c(1L, n_out))) {
    stop("O2, G, and N must have compatible lengths.")
  }
  O2_vec <- rep_len(as.numeric(O2_use), n_out)
  N_vec <- rep_len(N_use, n_out)
  G_vec <- rep_len(as.numeric(G_use), n_out)
  O2_crit_use <- as.numeric(.first_non_null(O2_crit, run_params$O2_crit, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0
  n_O <- as.numeric(.first_non_null(run_params$n_O, 1.0))
  if (!is.finite(n_O) || n_O < 0) stop("run_params$n_O must be finite and >= 0.")
  lam_min_use <- as.numeric(.first_non_null(run_params$lam_min, 0.0))
  lam_max_use <- as.numeric(.first_non_null(run_params$lam_max, lam_min_use))
  k_o_use <- as.numeric(.first_non_null(run_params$k_o, 1.0))
  if (!is.finite(k_o_use) || k_o_use <= 0) k_o_use <- 1e-12
  alpha_o2_use <- pmax(as.numeric(.first_non_null(run_params$alpha_o2, 0.0)), 0)
  gamma_growth_use <- pmax(as.numeric(.first_non_null(run_params$gamma_growth, 1.0)), 1e-12)
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null(run_params$glucose, TRUE),
    default = TRUE
  ))
  o2_c <- pmax(O2_crit_use, 1e-12)
  h_o2 <- (o2_c^n_O) / ((o2_c^n_O) + (pmax(O2_vec, 0)^n_O))
  h_o2 <- .clip01(h_o2)

  if (!isTRUE(glucose_use)) {
    frac <- O2_vec / (O2_vec + pmax(k_o_use, 1e-12))
    lam_base <- lam_min_use + (lam_max_use - lam_min_use) * frac
    if (!isTRUE(O2_growth)) return(pmax(lam_base, 0))
    denom <- 1 + alpha_o2_use * h_o2 * ((pmax(N_vec, 0) / 44)^gamma_growth_use)
    return(pmax(lam_base / pmax(denom, 1e-12), 0))
  }

  h_g <- (o2_c^n_O) / ((o2_c^n_O) + (pmax(G_vec, 0)^n_O))
  h_g <- .clip01(h_g)
  R_resource <- .clip01((1 - h_o2) * (1 - h_g))
  lam_base <- lam_min_use + (lam_max_use - lam_min_use) * R_resource
  if (!isTRUE(O2_growth)) return(pmax(lam_base, 0))
  h_resource <- .clip01(1 - (1 - h_o2) * (1 - h_g))
  denom <- 1 + alpha_o2_use * h_resource * ((pmax(N_vec, 0) / 44)^gamma_growth_use)
  pmax(lam_base / pmax(denom, 1e-12), 0)
}

# Main-path baseline-plus-increment missegregation helper (aligned with C++).
# -----------------------------------------------------------------------------
# Function: .mu_eff_of_O2
# Purpose: Compute state-specific hypoxia death rate under the main model.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - N: Ploidy state value or chromosome-copy count.
#   - O2_crit: Optional hypoxia Hill critical scale override.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.mu_eff_of_O2 <- function(O2, run_params, N = 44, O2_crit = NULL, G = NULL) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  N_use <- as.numeric(N)
  if (any(!is.finite(N_use))) stop("N must be finite.")
  G_use <- if (is.null(G)) O2_use else .assert_o2_pct(G, label = "G")
  n_out <- max(length(O2_use), length(N_use), length(G_use))
  if (!(length(O2_use) %in% c(1L, n_out) &&
        length(N_use) %in% c(1L, n_out) &&
        length(G_use) %in% c(1L, n_out))) {
    stop("O2, G, and N must have compatible lengths.")
  }
  O2_vec <- rep_len(O2_use, n_out)
  N_vec <- rep_len(N_use, n_out)
  G_vec <- rep_len(G_use, n_out)

  h_resource <- .resource_stress_of_O2(
    O2 = O2_vec,
    run_params = run_params,
    O2_crit = O2_crit,
    G = G_vec
  )

  mu_hp_use <- as.numeric(.first_non_null(run_params$mu_hp, 0.0))
  if (!is.finite(mu_hp_use) || mu_hp_use < 0) mu_hp_use <- 0.0
  gamma_mu_use <- as.numeric(.first_non_null(run_params$gamma_mu, 1.0))
  if (!is.finite(gamma_mu_use) || gamma_mu_use <= 0) gamma_mu_use <- 1.0
  ploidy_O2_death_mode <- assert_canonical_ploidy_o2_death_mode(
    .first_non_null(run_params$ploidy_O2_death, "diploid_NULL")
  )

  if (identical(ploidy_O2_death_mode, "uniform")) {
    mu_eff <- mu_hp_use * h_resource
  } else if (identical(ploidy_O2_death_mode, "diploid_NULL")) {
    above_dip <- pmax(N_vec / 44.0 - 1.0, 0.0)
    mu_eff <- mu_hp_use * h_resource * (1.0 + (above_dip^gamma_mu_use))
  } else {
    ratio <- pmax(N_vec / 44.0, 0.0)
    mu_eff <- mu_hp_use * h_resource * (ratio^gamma_mu_use)
  }
  pmax(mu_eff, 0.0)
}

# Main-path death-linked O2/ploidy missegregation helper (aligned with C++).
# -----------------------------------------------------------------------------
# Function: .pmisseg_of_O2
# Purpose: Compute state-specific missegregation probability (baseline + death-linked increment).
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - N: Ploidy state value or chromosome-copy count.
#   - O2_crit: Optional hypoxia Hill critical scale override.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.pmisseg_of_O2 <- function(O2, run_params, N = 44, O2_crit = NULL, G = NULL) {
  mu_eff <- .mu_eff_of_O2(
    O2 = O2,
    run_params = run_params,
    N = N,
    O2_crit = O2_crit,
    G = G
  )

  p_base <- as.numeric(.first_non_null(run_params$p_mis_base, 1e-5))
  if (!is.finite(p_base) || p_base < 0) p_base <- 1e-5
  p_amp <- as.numeric(.first_non_null(run_params$p_misseg, 0.0))
  if (!is.finite(p_amp) || p_amp < 0) p_amp <- 0.0
  k_o_mis <- as.numeric(.first_non_null(run_params$k_o_mis, 50.0))
  if (!is.finite(k_o_mis) || k_o_mis <= 0) k_o_mis <- 1e-12
  frac <- mu_eff / (mu_eff + k_o_mis)
  delta_p <- p_amp * frac
  .clip01(p_base + delta_p)
}

# -----------------------------------------------------------------------------
# Function: .p_wgd_of_O2
# Purpose: Return the constant per-division WGD probability under the main model.
# Parameters:
#   - O2: Oxygen level used only to match diagnostic vector lengths.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
# Returns:
#   Numeric vector used by downstream diagnostics and plotting helpers.
# -----------------------------------------------------------------------------
.p_wgd_of_O2 <- function(O2, run_params) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  p_wgd_use <- as.numeric(.first_non_null(run_params$p_wgd, 0.0))
  if (!is.finite(p_wgd_use) || p_wgd_use < 0) p_wgd_use <- 0.0
  rep(.clip01(p_wgd_use), length(O2_use))
}

# Intrinsic-buffer delta weight formula (aligned with asymmetric_intrinsic_buffer).
# -----------------------------------------------------------------------------
# Function: .pr_delta_vec
# Purpose: Compute missegregation delta-kernel probabilities over ploidy shifts.
# Parameters:
#   - N: Ploidy state value or chromosome-copy count.
#   - p: Missegregation probability parameter.
#   - eps_tail: Small truncation threshold for tail probabilities.
#   - buffer_smax: Maximum per-copy survival factor.
#   - buffer_beta: Ploidy-buffering strength.
#   - buffer_n_exp: Ploidy-buffering exponent.
#   - N_unit: Number of modeled chromosome classes for buffering scale.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.pr_delta_vec <- function(N, p, eps_tail = 1e-8, mr_lethality = 0.9,
                          N_unit = 22L,
                          buffer_smax = 1.0, buffer_beta = 0.0,
                          buffer_n_exp = 1.0) {
  .require_cpp_o2simps_fn("cpp_o2simps_pr_delta_vec")
  res <- cpp_o2simps_pr_delta_vec(
    as.integer(N),
    as.numeric(p),
    eps_tail = as.numeric(eps_tail),
    buffer_smax = as.numeric(buffer_smax),
    buffer_beta = as.numeric(buffer_beta),
    buffer_n_exp = as.numeric(buffer_n_exp),
    N_unit = as.integer(N_unit)
  )
  out <- as.numeric(res$prob)
  names(out) <- as.character(res$ts)
  attr(out, "mass_dropped") <- as.numeric(res$mass_dropped)
  return(out)
}

# -----------------------------------------------------------------------------
# Function: .build_B_total
# Purpose: Build total missegregation transition operator on ploidy grid.
# Parameters:
#   - Nmin: Minimum ploidy state on source grid.
#   - Nmax: Maximum ploidy state on source grid.
#   - p_vec: State-specific missegregation probability vector.
#   - mr_lethality: Probability of lethal outcome after severe missegregation.
#   - boundary: Boundary handling mode when transitions leave the ploidy grid.
#   - eps_tail: Small truncation threshold for tail probabilities.
#   - return_sparse: Function-specific input argument.
#   - buffer_smax: Maximum per-copy survival factor.
#   - buffer_beta: Ploidy-buffering strength.
#   - buffer_n_exp: Ploidy-buffering exponent.
#   - N_unit: Number of modeled chromosome classes for buffering scale.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_B_total <- function(Nmin, Nmax, p_vec, mr_lethality = 0.9,
                           boundary = c("drop", "absorb_minmax"),
                           eps_tail = 1e-8, return_sparse = TRUE,
                           N_unit = 22L,
                           buffer_smax = 1.0, buffer_beta = 0.0,
                           buffer_n_exp = 1.0) {
  boundary <- match.arg(boundary)
  R <- Nmax - Nmin + 1L
  if (length(p_vec) == 1L) p_vec <- rep(p_vec, R)

  .require_cpp_o2simps_fn("cpp_o2simps_build_B_total_triplet")
  tri <- cpp_o2simps_build_B_total_triplet(
    as.integer(Nmin),
    as.integer(Nmax),
    as.numeric(p_vec),
    boundary = boundary,
    eps_tail = as.numeric(eps_tail),
    buffer_smax = as.numeric(buffer_smax),
    buffer_beta = as.numeric(buffer_beta),
    buffer_n_exp = as.numeric(buffer_n_exp),
    N_unit = as.integer(N_unit)
  )
  B <- sparseMatrix(
    i = as.integer(tri$i),
    j = as.integer(tri$j),
    x = as.numeric(tri$x),
    dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
    repr = "C"
  )
  return(if (isTRUE(return_sparse)) B else as.matrix(B))
}

# -----------------------------------------------------------------------------
# Function: .build_B_WGD
# Purpose: Build WGD transition operator between source and doubled-ploidy grids.
# Parameters:
#   - N0min: Minimum ploidy state on source grid.
#   - N0max: Maximum ploidy state on source grid.
#   - N1min: Minimum ploidy state on target grid for doubled states.
#   - N1max: Maximum ploidy state on target grid for doubled states.
#   - boundary: Boundary handling mode when transitions leave the ploidy grid.
#   - return_sparse: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_B_WGD <- function(N0min, N0max, N1min, N1max,
                         boundary = c("drop", "absorb_minmax"),
                         return_sparse = TRUE) {
  boundary <- match.arg(boundary)
  R0 <- N0max - N0min + 1L
  R1 <- N1max - N1min + 1L

  .require_cpp_o2simps_fn("cpp_o2simps_build_B_WGD_triplet")
  tri <- cpp_o2simps_build_B_WGD_triplet(
    as.integer(N0min),
    as.integer(N0max),
    as.integer(N1min),
    as.integer(N1max),
    boundary = boundary,
    wgd_value = 1.0
  )
  B <- sparseMatrix(
    i = as.integer(tri$i),
    j = as.integer(tri$j),
    x = as.numeric(tri$x),
    dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
    repr = "C"
  )
  return(if (isTRUE(return_sparse)) B else as.matrix(B))
}

# -----------------------------------------------------------------------------
# Function: .build_G_with_WGD
# Purpose: Build generator matrix at the current oxygen/burden condition.
# Parameters:
#   - N0min: Minimum ploidy state on source grid.
#   - N0max: Maximum ploidy state on source grid.
#   - lambda0_vec: Function-specific input argument.
#   - p0_vec: Function-specific input argument.
#   - wgd_prob_vec: Function-specific input argument.
#   - mr_lethality0: Function-specific input argument.
#   - mr_buffer_by_ploidy: Function-specific input argument.
#   - N_unit: Ploidy scaling unit used to map integer states to N values.
#   - P_low: Function-specific input argument.
#   - P_high: Function-specific input argument.
#   - boundary: Boundary handling mode when transitions leave the ploidy grid.
#   - eps_tail: Small truncation threshold for tail probabilities.
#   - buffer_smax: Maximum per-copy survival factor.
#   - buffer_beta: Ploidy-buffering strength.
#   - buffer_n_exp: Ploidy-buffering exponent.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_G_with_WGD <- function(
    N0min, N0max, lambda0_vec, p0_vec, wgd_prob_vec,
    mr_lethality0 = 0.9, mr_lethality1 = 0.9,
    mr_buffer_by_ploidy = TRUE, N_unit = 22L, P_low = 2.0, P_high = 4.0,
    boundary = "drop", eps_tail = 1e-8,
    buffer_smax = 1.0, buffer_beta = 0.0,
    buffer_n_exp = 1.0
) {
  R0 <- N0max - N0min + 1L
  if (length(lambda0_vec) == 1L) lambda0_vec <- rep(lambda0_vec, R0)
  if (length(p0_vec) == 1L) p0_vec <- rep(p0_vec, R0)
  if (length(wgd_prob_vec) == 1L) wgd_prob_vec <- rep(wgd_prob_vec, R0)
  wgd_prob_vec <- .clip01(wgd_prob_vec)

  B0 <- .build_B_total(
    N0min, N0max, p_vec = p0_vec,
    boundary = boundary, eps_tail = eps_tail,
    N_unit = N_unit,
    buffer_smax = buffer_smax,
    buffer_beta = buffer_beta,
    buffer_n_exp = buffer_n_exp
  )
  BW <- .build_B_WGD(N0min, N0max, N0min, N0max, boundary = boundary)
  L0 <- Diagonal(x = lambda0_vec)
  S0 <- Diagonal(x = (1 - wgd_prob_vec))
  SW <- Diagonal(x = wgd_prob_vec)
  ((B0 %*% S0) %*% L0) + ((BW %*% SW) %*% L0) - L0
}

# -----------------------------------------------------------------------------
# Function: run_all_sims
# Purpose: Run forward simulations for all configured scenarios and collect outputs.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
run_all_sims <- function(run_params) {
  all_results_list <- list()
  passage_times <- list()

  init_P_2N <- x$P[x$passage == 0 & x$ploidy == "2N"]
  init_P_4N <- x$P[x$passage == 0 & x$ploidy == "4N"]

  buffer_smax <- as.numeric(.first_non_null(run_params$buffer_smax, 1.0))
  buffer_beta <- as.numeric(.first_non_null(run_params$buffer_beta, 0.0))
  buffer_n_exp <- as.numeric(.first_non_null(run_params$buffer_n_exp, 1.0))
  boundary_mode <- as.character(.first_non_null(run_params$boundary, "drop"))
  O2_crit_use <- as.numeric(.first_non_null(run_params$O2_crit, 1.0))
  if (!is.finite(O2_crit_use) || O2_crit_use < 0) O2_crit_use <- 1.0

  .require_cpp_o2simps_fn("cpp_o2simps_build_G_for_o2_triplet")
  lam_min_use <- as.numeric(.first_non_null(run_params$lam_min, 1.0))
  lam_max_use <- as.numeric(.first_non_null(run_params$lam_max, lam_min_use))
  k_o_use <- as.numeric(.first_non_null(run_params$k_o, 50.0))
  has_p_misseg <- !is.null(run_params$p_misseg)
  mu_hp_use <- as.numeric(.first_non_null(run_params$mu_hp, 0.0))
  gamma_mu_use <- as.numeric(.first_non_null(run_params$gamma_mu, 1.0))
  n_O_use <- as.numeric(.first_non_null(run_params$n_O, 1.0))
  o2_Nref_use <- as.numeric(.first_non_null(run_params$o2_Nref, if (exists("cfg", inherits = TRUE)) get("cfg", inherits = TRUE)$o2_Nref else NULL, 1e6))
  cfg_o2_growth <- if (exists("cfg", inherits = TRUE)) get("cfg", inherits = TRUE)$O2_growth else NULL
  o2_growth_use <- isTRUE(.first_non_null(run_params$O2_growth, cfg_o2_growth, TRUE))
  ploidy_O2_death_mode_use <- assert_canonical_ploidy_o2_death_mode(
    .first_non_null(run_params$ploidy_O2_death, "diploid_NULL")
  )
  glucose_use <- isTRUE(canonical_glucose_enabled(
    .first_non_null(run_params$glucose, TRUE),
    default = TRUE
  ))
  if (!is.finite(mu_hp_use) || mu_hp_use < 0) mu_hp_use <- 0.0
  if (!is.finite(gamma_mu_use) || gamma_mu_use <= 0) gamma_mu_use <- 1.0
  if (!is.finite(n_O_use) || n_O_use < 0) stop("run_params$n_O must be finite and >= 0.")
  if (!is.finite(o2_Nref_use) || o2_Nref_use <= 0) o2_Nref_use <- 1e6

  G_cache <- new.env(parent = emptyenv())
  build_G_for_resources <- function(O2) {
    O2_use <- .assert_o2_pct(as.numeric(O2), label = "O2")
    key <- sprintf("%.3f", O2_use)
    if (!exists(key, envir = G_cache, inherits = FALSE)) {
      tri <- cpp_o2simps_build_G_for_o2_triplet(
        O2 = as.numeric(O2_use),
        O2_crit = as.numeric(O2_crit_use),
        N0min = as.integer(N_MIN),
        N0max = as.integer(N_MAX),
        N1min = as.integer(N_MIN),
        N1max = as.integer(N_MAX),
        lam_min = as.numeric(lam_min_use),
        lam_max = as.numeric(lam_max_use),
        k_o = as.numeric(k_o_use),
        has_p_misseg = isTRUE(has_p_misseg),
        p_mis_base = as.numeric(.first_non_null(run_params$p_mis_base, 1e-5)),
        p_misseg = as.numeric(.first_non_null(run_params$p_misseg, 0.0)),
        k_o_mis = as.numeric(.first_non_null(run_params$k_o_mis, 50.0)),
        has_pmis_endpoints = FALSE,
        pmis_O2_0 = 0.0,
        pmis_O2_1 = 0.0,
        p_const = 0.0,
        glucose = isTRUE(glucose_use),
        p_wgd = as.numeric(.first_non_null(run_params$p_wgd, 0.0)),
        p_wgd_max = as.numeric(.first_non_null(run_params$p_wgd_max, 0.0)),
        O2_wgd = as.numeric(.first_non_null(run_params$O2_wgd, 0.1)),
        boundary = as.character(boundary_mode),
        eps_tail = as.numeric(1e-8),
        buffer_smax = as.numeric(buffer_smax),
        buffer_beta = as.numeric(buffer_beta),
        buffer_n_exp = as.numeric(buffer_n_exp),
        N_unit = as.integer(N_UNIT),
        beta_size = 0.0,
        O2_growth = isTRUE(o2_growth_use),
        alpha_o2 = as.numeric(.first_non_null(run_params$alpha_o2, 0.0)),
        gamma_growth = as.numeric(.first_non_null(run_params$gamma_growth, 1.0)),
        mu_hp = as.numeric(mu_hp_use),
        gamma_mu = as.numeric(gamma_mu_use),
        n_O = as.numeric(n_O_use),
        ploidy_O2_death = as.character(ploidy_O2_death_mode_use)
      )
      G <- sparseMatrix(
        i = as.integer(tri$i),
        j = as.integer(tri$j),
        x = as.numeric(tri$x),
        dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
        repr = "C"
      )
      assign(key, G, envir = G_cache)
    }
    get(key, envir = G_cache, inherits = FALSE)
  }

  for (sim in sim_configs) {
    O2_LEVEL <- .assert_o2_pct(as.numeric(sim$O2), label = "sim$O2")

    init_P_values <- if (sim$init_ploidy == "2N") init_P_2N else init_P_4N
    x_current <- create_initial_dist(
      init_P_values,
      grid_pre,
      N_UNIT,
      chr_lengths_bp = default_chr_lengths_bp_1to22()
    )
    x_current <- x_current / sum(x_current)

    sim_passage_times <- numeric(PASSAGES_TO_RUN)

    for (p in 1:PASSAGES_TO_RUN) {
      pop_start <- sum(x_current)
      pop_target <- pop_start * POP_GROWTH_FACTOR
      time_in_passage <- 0.0

      while (sum(x_current) < pop_target) {
        x_prev <- as.numeric(x_current)
        mu_vec_step <- as.numeric(.mu_eff_of_O2(
          O2 = O2_LEVEL,
          run_params = run_params,
          N = grid_pre,
          O2_crit = O2_crit_use
        ))
        G <- build_G_for_resources(O2_LEVEL)
        x_div <- step_dt(G, x_prev, DT, 1L)
        x_next <- x_div - DT * mu_vec_step * x_prev
        x_next[!is.finite(x_next) | x_next < 0] <- 0
        x_current <- x_next
        time_in_passage <- time_in_passage + DT
        if (sum(x_current) < pop_start * 1e-3 || time_in_passage > 1000) {
          break
        }
      }
      sim_passage_times[p] <- time_in_passage

      if (p %in% REPORT_PASSAGES) {
        pop_total <- sum(x_current)
        dist_df <- data.frame(
          sim_id = sim$id,
          passage = p,
          layer = "single",
          N = grid_pre,
          fraction = x_current / pop_total
        )
        all_results_list[[length(all_results_list) + 1]] <- dist_df
      }
      x_current <- x_current / sum(x_current) * pop_start
    }

    passage_times[[sim$id]] <- data.frame(
      sim_id = sim$id,
      passage = 1:PASSAGES_TO_RUN,
      duration = sim_passage_times
    )
  }

  all_dists <- do.call(rbind, all_results_list)
  all_passage_times <- do.call(rbind, passage_times)

  list(all_dists = all_dists, all_passage_times = all_passage_times)
}

# -----------------------------------------------------------------------------
# Function: plot_misseg_interp
# Purpose: Plot missegregation response curve under current parameterization.
# Parameters:
#   - par: Function-specific input argument.
#   - o2_ref: Reference oxygen level used for plotting interpolation curves.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
plot_misseg_interp <- function(par, o2_ref = 20.5) {
  O2 <- seq(0, 100, length.out = 401)
  p <- .pmisseg_of_O2(
    O2 = O2,
    run_params = par,
    N = 44,
    O2_crit = as.numeric(.first_non_null(par$O2_crit, 1.0))
  )
  df <- data.frame(O2_pct = O2, p = as.numeric(p))
  ggplot(df, aes(O2_pct, p)) +
    geom_line(linewidth = 1, color = "black") +
    geom_point(data = df[c(1, nrow(df)), ], size = 2, color = "red") +
    labs(x = "Oxygen (%)", y = "Missegregation rate") +
    theme_bw()
}
