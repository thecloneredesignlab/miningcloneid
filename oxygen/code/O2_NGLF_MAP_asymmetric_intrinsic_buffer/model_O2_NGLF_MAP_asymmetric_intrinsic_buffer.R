suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

# ----------------------------------------------------------------------------
# Align miningcloneid oxygen model to Richard's buffering.R logic.
# model_O2_NGLF_MAP_asymmetric_intrinsic_buffer extension:
# - Keep karyotype dynamics identical to Richard-aligned model.
# - Replace burden->O2 feedback with asymmetric sigmoid rise:
#   gompertz or generalized logistic (selected by one mode flag).
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
      stop("Rcpp package is required for model_O2_NGLF_MAP_asymmetric_intrinsic_buffer but is not installed.")
    }

    # Dedicated O2 invivo backend (do not use Richard/shared cpp here).
    cpp_path <- file.path(.ALIGN_MODEL_DIR, "model_O2_NGLF_MAP_asymmetric_intrinsic_buffer.cpp")
    if (!file.exists(cpp_path)) {
      stop("Cannot find required C++ backend file: ", cpp_path)
    }

    cache_root <- file.path(.ALIGN_MODEL_DIR, ".rcpp_cache_o2_nglf_map_asymmetric_intrinsic_buffer")
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
      stop("Failed to compile/load model_O2_NGLF_MAP_asymmetric_intrinsic_buffer.cpp: ", conditionMessage(e))
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
        "model_O2_NGLF_MAP_asymmetric_intrinsic_buffer C++ backend loaded but required symbols are missing: ",
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
      must_have = c("O2_cap", "beta_loss"),
      must_absent = c("o2_ref_pct")
    )
    check_wrapper_formals(
      "cpp_o2simps_simulate_one",
      must_have = c("O2_cap", "beta_loss"),
      must_absent = c("o2_ref_pct")
    )
    check_wrapper_formals(
      "cpp_o2simps_objective_components_map",
      must_have = c("O2_cap", "beta_loss"),
      must_absent = c("o2_ref_pct")
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
          "Failed forced rebuild for model_O2_NGLF_MAP_asymmetric_intrinsic_buffer.cpp after wrapper mismatch [",
          paste(wrapper_mismatch_reason, collapse = "; "),
          "]: ", conditionMessage(e)
        )
      })

      wrapper_mismatch_reason <- character(0)
      wrappers_need_rebuild <- FALSE
      check_wrapper_formals(
        "cpp_o2simps_build_G_for_o2_triplet",
        must_have = c("O2_cap", "beta_loss"),
        must_absent = c("o2_ref_pct")
      )
      check_wrapper_formals(
        "cpp_o2simps_simulate_one",
        must_have = c("O2_cap", "beta_loss"),
        must_absent = c("o2_ref_pct")
      )
      check_wrapper_formals(
        "cpp_o2simps_objective_components_map",
        must_have = c("O2_cap", "beta_loss"),
        must_absent = c("o2_ref_pct")
      )
      if (isTRUE(wrappers_need_rebuild)) {
        stop(
          "model_O2_NGLF_MAP_asymmetric_intrinsic_buffer wrapper signatures are inconsistent after forced rebuild: ",
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
  cache_pat <- ".rcpp_cache_o2_nglf_map_asymmetric_intrinsic_buffer"
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
    stop("Unable to resolve sourceCpp wrapper file (*.cpp.R) for O2_NGLF_MAP_asymmetric_intrinsic_buffer backend.")
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
# Purpose: Map total ploidy state N to observed ploidy scale under the
#   single-variable total-ploidy model.
# Parameters:
#   - N_total: Total ploidy state(s) on the integer grid.
#   - chr_lengths_bp: Deprecated/unused in this model.
# Returns:
#   Numeric ploidy values P = N / N_unit (with N_unit fixed at 22 autosomes).
# -----------------------------------------------------------------------------
weighted_ploidy_from_total_N <- function(N_total, chr_lengths_bp = NULL) {
  N_unit <- 22.0
  Nv <- as.numeric(N_total)
  vapply(Nv, function(nn) {
    if (!is.finite(nn)) return(NA_real_)
    as.numeric(nn) / N_unit
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
  p_grid <- weighted_ploidy_from_total_N(grid, chr_lengths_bp = NULL)
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
#   - N_unit: Ploidy scaling unit used to map integer states to N values.
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

# Build explicit initial state vector on a single total-ploidy axis.
# - Returns absolute cell counts with total mass = total_size.
# -----------------------------------------------------------------------------
# Function: make_init_state
# Purpose: Construct initial ploidy-state distribution/state vector for simulation.
# Parameters:
#   - N_grid: Integer ploidy grid for the single N-axis.
#   - ploidy: Initial ploidy mode ("2" or "4").
#   - N_UNIT: Ploidy scaling unit used to map N to P=N/N_UNIT.
#   - total_size: Total initial mass.
#   - chr_lengths_bp: Optional chromosome-length vector for weighted ploidy mapping.
# Returns:
#   Numeric state vector aligned to N_grid.
# -----------------------------------------------------------------------------
make_init_state <- function(N_grid,
                            ploidy = c(2, 4),
                            N_UNIT = 22L,
                            total_size = 1e6,
                            chr_lengths_bp = NULL) {
  ploidy <- match.arg(as.character(ploidy), choices = c("2", "4"))
  Pnum <- as.numeric(ploidy)

  grid_use <- as.integer(sort(unique(N_grid)))
  x_vec <- rep(0, length(grid_use))
  names(x_vec) <- grid_use

  N_delta <- as.integer(map_ploidy_to_N_by_chrlen(
    ploidy_values = Pnum,
    N_grid = grid_use,
    chr_lengths_bp = chr_lengths_bp
  ))
  stopifnot(N_delta %in% grid_use)
  x_vec[as.character(N_delta)] <- 1

  x_vec * as.numeric(total_size)
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
# Function: .sigmoid01
# Purpose: Internal helper used by the model fitting and simulation pipeline.
# Parameters:
#   - z: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.sigmoid01 <- function(z) 1 / (1 + exp(-z))
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

# O2(N) monotone-decreasing, size-aware NGLF:
# - xi_size = (1/3) * [log10(N + eps) - log10(N_anchor + eps)]
# - glogistic: O2 = O2_min + (O2_cap - O2_min) / (1 + exp(rate * (xi_size - xi0)))^(1/v)
# - gompertz:  O2 = O2_min + (O2_cap - O2_min) * exp(-exp(rate * (xi_size - xi0_g)))
# where xi0/xi0_g are solved from the anchor condition O2(N_anchor) = o2_init_pct.
# -----------------------------------------------------------------------------
# Function: .o2_sigmoid_supply_from_burden
# Purpose: Compute oxygen supply fraction/level from burden under the selected O2 model.
# Parameters:
#   - Ntot: Total predicted cell count (or burden proxy) at current time.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - O2_cap: Function-specific input argument.
#   - o2_logN_eps: Function-specific input argument.
#   - o2_anchor_N: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.o2_sigmoid_supply_from_burden <- function(Ntot,
                                           run_params,
                                           O2_cap = 5.0,
                                           o2_logN_eps = 1.0,
                                           o2_anchor_N = 1e6) {
  Ntot_use <- pmax(as.numeric(Ntot), 0)
  curve_type <- as.character(.first_non_null(run_params$o2_curve_type, "gompertz"))
  O2_cap_use <- .assert_o2_pct(as.numeric(.first_non_null(run_params$o2_cap, O2_cap)), label = "o2_cap")
  o2_init <- as.numeric(.first_non_null(run_params$o2_init_pct, 0.5))
  o2_rate <- as.numeric(.first_non_null(run_params$o2_rate, 1.0))
  o2_shape_v <- as.numeric(.first_non_null(run_params$o2_shape_v, 1.0))
  if (!is.finite(o2_rate) || o2_rate <= 0) o2_rate <- 1.0
  if (!is.finite(o2_shape_v) || o2_shape_v <= 0) o2_shape_v <- 1.0
  if (!is.finite(o2_anchor_N) || o2_anchor_N < 0) o2_anchor_N <- 1e6
  eps_use <- as.numeric(o2_logN_eps)
  if (!is.finite(eps_use) || eps_use <= 0) eps_use <- 1.0
  o2_eps <- 1e-9
  o2_init <- max(o2_eps, min(O2_cap_use - o2_eps, o2_init))

  .require_cpp_o2simps_fn("cpp_o2simps_o2_window_supply")
  return(as.numeric(cpp_o2simps_o2_window_supply(
    Ntot = as.numeric(Ntot_use),
    curve_type = as.character(curve_type),
    O2_cap = as.numeric(O2_cap_use),
    o2_init = as.numeric(o2_init),
    o2_rate = as.numeric(o2_rate),
    o2_shape_v = as.numeric(o2_shape_v),
    o2_anchor_N = as.numeric(o2_anchor_N),
    o2_logN_eps = as.numeric(eps_use)
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

# Richard buffering.R-style O2-dependent missegregation.
# Endpoint interpolation branch uses O2 in percent scale [0, 100].
# -----------------------------------------------------------------------------
# Function: .pmisseg_of_O2
# Purpose: Compute oxygen-dependent missegregation probability.
# Parameters:
#   - O2: Oxygen level used by model rate functions.
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.pmisseg_of_O2 <- function(O2, run_params) {
  O2_use <- .assert_o2_pct(O2, label = "O2")
  p0 <- as.numeric(run_params$p_misseg)
  k_o_mis <- as.numeric(run_params$k_o_mis)
  if (!is.finite(p0) || p0 < 0) stop("run_params$p_misseg must be finite and >= 0.")
  if (!is.finite(k_o_mis) || k_o_mis <= 0) stop("run_params$k_o_mis must be > 0.")
  p <- p0 * (1 - (O2_use / (O2_use + k_o_mis)))
  .clip01(p)
}

# Richard buffering.R delta weight formula.
# -----------------------------------------------------------------------------
# Function: .pr_delta_vec
# Purpose: Compute missegregation delta-kernel probabilities over ploidy shifts.
# Parameters:
#   - N: Ploidy state value or chromosome-copy count.
#   - p: Missegregation probability parameter.
#   - eps_tail: Small truncation threshold for tail probabilities.
#   - mr_lethality: Probability of lethal outcome after severe missegregation.
#   - beta_buffer: Legacy placeholder; inactive in intrinsic-buffer variant.
#   - n_exp: Legacy placeholder; inactive in intrinsic-buffer variant.
#   - smax: Legacy placeholder; inactive in intrinsic-buffer variant.
#   - N_unit: Ploidy scaling unit used to map integer states to N values.
#   - beta_loss: Loss-burden penalty coefficient.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.pr_delta_vec <- function(N, p, eps_tail = 1e-8, mr_lethality = 0.9,
                          beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L,
                          beta_loss = 0.25) {
  .require_cpp_o2simps_fn("cpp_o2simps_pr_delta_vec")
  res <- cpp_o2simps_pr_delta_vec(
    as.integer(N),
    as.numeric(p),
    eps_tail = as.numeric(eps_tail),
    beta_buffer = as.numeric(beta_buffer),
    n_exp = as.numeric(n_exp),
    smax = as.numeric(smax),
    N_unit = as.integer(N_unit),
    beta_loss = as.numeric(beta_loss)
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
#   - beta_buffer: Legacy placeholder; inactive in intrinsic-buffer variant.
#   - n_exp: Legacy placeholder; inactive in intrinsic-buffer variant.
#   - smax: Legacy placeholder; inactive in intrinsic-buffer variant.
#   - N_unit: Ploidy scaling unit used to map integer states to N values.
#   - beta_loss: Loss-burden penalty coefficient.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_B_total <- function(Nmin, Nmax, p_vec, mr_lethality = 0.9,
                           boundary = c("drop", "absorb_minmax"),
                           eps_tail = 1e-8, return_sparse = TRUE,
                           beta_buffer = 0.0, n_exp = 1.0, smax = 1.0, N_unit = 22L,
                           beta_loss = 0.25) {
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
    beta_buffer = as.numeric(beta_buffer),
    n_exp = as.numeric(n_exp),
    smax = as.numeric(smax),
    N_unit = as.integer(N_unit),
    beta_loss = as.numeric(beta_loss)
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
# Purpose: Build direct WGD transition on single N-axis: N -> 2N.
# Parameters:
#   - N0min: Minimum ploidy state on source grid.
#   - N0max: Maximum ploidy state on source grid.
#   - boundary: Boundary handling mode when transitions leave the ploidy grid.
#   - return_sparse: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
.build_B_WGD <- function(N0min, N0max,
                         boundary = c("drop", "absorb_minmax"),
                         return_sparse = TRUE) {
  boundary <- match.arg(boundary)

  .require_cpp_o2simps_fn("cpp_o2simps_build_B_WGD_triplet")
  tri <- cpp_o2simps_build_B_WGD_triplet(
    as.integer(N0min),
    as.integer(N0max),
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
.build_G_with_WGD <- function(
    N0min, N0max, lambda0_vec, p0_vec, wgd_prob_vec,
    N_unit = 22L,
    boundary = "drop", eps_tail = 1e-8,
    beta_buffer = 0.0, n_exp = 1.0, smax = 1.0,
    beta_loss = 0.25
) {
  R <- N0max - N0min + 1L
  if (length(lambda0_vec) == 1L) lambda0_vec <- rep(lambda0_vec, R)
  if (length(p0_vec) == 1L) p0_vec <- rep(p0_vec, R)
  if (length(wgd_prob_vec) == 1L) wgd_prob_vec <- rep(wgd_prob_vec, R)
  stopifnot(length(lambda0_vec) == R, length(p0_vec) == R, length(wgd_prob_vec) == R)
  wgd_prob_vec <- .clip01(wgd_prob_vec)

  B <- .build_B_total(
    N0min, N0max, p_vec = p0_vec,
    boundary = boundary, eps_tail = eps_tail,
    beta_buffer = beta_buffer, n_exp = n_exp, smax = smax, N_unit = N_unit,
    beta_loss = beta_loss
  )
  BW <- .build_B_WGD(N0min, N0max, boundary = boundary)

  L <- Diagonal(x = lambda0_vec)
  S0 <- Diagonal(x = (1 - wgd_prob_vec))
  SW <- Diagonal(x = wgd_prob_vec)

  (B %*% S0) %*% L + (BW %*% SW) %*% L - L
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
  stop(
    "run_all_sims() has been retired in this single-variable total-ploidy variant. ",
    "Use cpp_o2simps_simulate_one() via the fit/prediction pipeline."
  )
}

# -----------------------------------------------------------------------------
# Function: run_in_vivo_crowd
# Purpose: Run in vivo simulation pipeline and return aggregated trajectory outputs.
# Parameters:
#   - run_params: Model parameters on natural scale used by simulation and loss evaluation.
#   - O2_schedule: Function-specific input argument.
#   - T_end: Function-specific input argument.
#   - sample_days: Function-specific input argument.
#   - N_UNIT: Function-specific input argument.
#   - DT: Function-specific input argument.
#   - K: Function-specific input argument.
#   - crowding: Function-specific input argument.
#   - N_grid: Integer ploidy grid for the single N-axis.
#   - init_state: Function-specific input argument.
# Returns:
#   Object used by downstream model fitting/simulation steps.
# -----------------------------------------------------------------------------
run_in_vivo_crowd <- function(run_params,
                              O2_schedule = list(c(t0 = 0, t1 = Inf, O2 = 5.0)),
                              T_end = 28, sample_days = c(0, 7, 14, 21, 28),
                              N_UNIT = 22L, DT = 0.1,
                              K = 1e9, crowding = c("logistic", "gompertz"),
                              N_grid,
                              init_state,
                              chr_lengths_bp = default_chr_lengths_bp_1to22()) {
  crowding <- match.arg(crowding)
  grid_use <- as.integer(sort(unique(N_grid)))
  if (length(grid_use) == 0L) stop("run_in_vivo_crowd: empty ploidy grid.")
  N0min <- min(grid_use)
  N0max <- max(grid_use)
  R <- length(grid_use)

  init_use <- as.numeric(init_state)
  if (length(init_use) != R) {
    stop("run_in_vivo_crowd: init_state length must equal length(N_grid).")
  }

  beta_buffer <- as.numeric(.first_non_null(run_params$beta_buffer, 0.0))
  n_exp <- as.numeric(.first_non_null(run_params$n_exp, 1.0))
  smax <- as.numeric(.first_non_null(run_params$smax, 1.0))
  beta_loss <- as.numeric(.first_non_null(run_params$beta_loss, 0.25))
  boundary_mode <- as.character(.first_non_null(run_params$boundary, "drop"))
  pwgd_val <- as.numeric(.first_non_null(run_params$p_wgd, 0))
  o2_burden_feedback <- isTRUE(.first_non_null(run_params$o2_burden_feedback, TRUE))
  o2_cap <- .assert_o2_pct(as.numeric(.first_non_null(run_params$o2_cap, 5.0)), label = "o2_cap")
  tau_O2_use <- as.numeric(.first_non_null(run_params$tau_O2, 2.0))
  if (!is.finite(tau_O2_use) || tau_O2_use <= 0) tau_O2_use <- 2.0
  alpha_tau <- 1 - exp(-DT / tau_O2_use)
  rho_2N <- as.numeric(.first_non_null(run_params$rho_2N, 4e4))
  if (!is.finite(rho_2N) || rho_2N <= 0) rho_2N <- 4e4
  beta_size_vol <- as.numeric(.first_non_null(run_params$beta_size, 1.0))
  if (!is.finite(beta_size_vol)) beta_size_vol <- 1.0
  p_weighted <- weighted_ploidy_from_total_N(grid_use, chr_lengths_bp = chr_lengths_bp)
  vol_by_N <- (1 / rho_2N) * (p_weighted / 2)^beta_size_vol

  get_O2 <- function(t) {
    for (seg in O2_schedule) {
      if (t >= seg["t0"] && t < seg["t1"]) return(as.numeric(seg["O2"]))
    }
    as.numeric(O2_schedule[[length(O2_schedule)]]["O2"])
  }

  live_burden <- function(v_state) {
    sum(as.numeric(v_state) * as.numeric(vol_by_N), na.rm = TRUE)
  }

  apply_O2_feedback <- function(O2_base, burden_live) {
    O2_base <- .assert_o2_pct(as.numeric(O2_base), label = "O2_schedule value")
    if (!o2_burden_feedback) return(O2_base)
    .o2_sigmoid_supply_from_burden(
      Ntot = burden_live,
      run_params = run_params,
      O2_cap = min(o2_cap, O2_base),
      o2_logN_eps = 1.0,
      o2_anchor_N = as.numeric(.first_non_null(run_params$o2_anchor_N, sum(init_use), 1e6))
    )
  }

  G_cache <- new.env(parent = emptyenv())
  build_G_for_O2 <- function(O2) {
    O2_use <- .assert_o2_pct(as.numeric(O2), label = "O2")
    key <- sprintf("%.3f", O2_use)
    if (!exists(key, envir = G_cache, inherits = FALSE)) {
      .require_cpp_o2simps_fn("cpp_o2simps_build_G_for_o2_triplet")

      lam_min_use <- as.numeric(.first_non_null(run_params$lam_min, 1.0))
      lam_max_use <- as.numeric(.first_non_null(run_params$lam_max, lam_min_use))
      k_o_use <- as.numeric(.first_non_null(run_params$k_o, 50.0))
      has_p_misseg <- !is.null(run_params$p_misseg)
      death_on <- isTRUE(.first_non_null(run_params$death, TRUE))
      mu_hp_use <- if (death_on) as.numeric(.first_non_null(run_params$mu_hp, 0.0)) else 0.0
      if (!is.finite(mu_hp_use) || mu_hp_use < 0) mu_hp_use <- 0.0

      tri <- cpp_o2simps_build_G_for_o2_triplet(
        O2 = as.numeric(O2_use),
        O2_cap = as.numeric(o2_cap),
        N0min = as.integer(N0min),
        N0max = as.integer(N0max),
        lam_min = as.numeric(lam_min_use),
        lam_max = as.numeric(lam_max_use),
        k_o = as.numeric(k_o_use),
        has_p_misseg = isTRUE(has_p_misseg),
        p_misseg = as.numeric(.first_non_null(run_params$p_misseg, 0.0)),
        k_o_mis = as.numeric(.first_non_null(run_params$k_o_mis, 50.0)),
        has_pmis_endpoints = FALSE,
        pmis_O2_0 = 0.0,
        pmis_O2_1 = 0.0,
        p_const = 0.0,
        p_wgd = as.numeric(pwgd_val),
        boundary = as.character(boundary_mode),
        eps_tail = as.numeric(1e-8),
        beta_buffer = as.numeric(beta_buffer),
        n_exp = as.numeric(n_exp),
        smax = as.numeric(smax),
        N_unit = as.integer(N_UNIT),
        beta_loss = as.numeric(beta_loss),
        beta_size = as.numeric(.first_non_null(run_params$beta_size, 0.0)),
        alpha_o2 = as.numeric(.first_non_null(run_params$alpha_o2, 0.0)),
        gamma_growth = as.numeric(.first_non_null(run_params$gamma_growth, 1.0)),
        growth_penalty_ploidy = isTRUE(.first_non_null(run_params$growth_penalty_ploidy, FALSE)),
        growth_penalty_hypoxia = isTRUE(.first_non_null(run_params$growth_penalty_hypoxia, FALSE)),
        mu_hp = as.numeric(mu_hp_use)
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

  crowd <- function(Ntot) {
    if (crowding == "logistic") return(max(0, 1 - Ntot / K))
    exp(-Ntot / K)
  }

  I <- Diagonal(n = length(init_use))
  v <- init_use
  times <- seq(0, T_end, by = DT)
  snapshots <- list()
  size_trace <- data.frame(day = 0, Ntot = sum(v))
  O2_state <- apply_O2_feedback(get_O2(0), live_burden(v))

  for (t in times) {
    if (t %in% sample_days) {
      snapshots[[as.character(t)]] <- data.frame(
        day = t,
        layer = rep("live", R),
        N = grid_use,
        fraction = v / max(sum(v), 1e-300),
        pop = sum(v)
      )
    }
    if (t >= T_end) break
    Ntot <- sum(v)
    O2_target <- apply_O2_feedback(get_O2(t), live_burden(v))
    O2_state <- O2_state + alpha_tau * (O2_target - O2_state)
    O2t <- .assert_o2_pct(as.numeric(O2_state), label = "O2_eff")
    G <- build_G_for_O2(O2t)
    cfac <- crowd(Ntot)
    v <- as.numeric((I + DT * (cfac * G)) %*% v)
    v[!is.finite(v) | v < 0] <- 0
    size_trace <- rbind(size_trace, data.frame(day = t + DT, Ntot = sum(v)))
    if (sum(v) <= 1e-9) break
  }

  list(all_dists = do.call(rbind, snapshots), tumor_size = size_trace)
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
  p0 <- as.numeric(par$p_misseg)
  k_o_mis <- as.numeric(.first_non_null(par$k_o_mis, 50.0))
  if (!is.finite(p0) || p0 < 0) stop("par$p_misseg must be finite and >= 0.")
  if (!is.finite(k_o_mis) || k_o_mis <= 0) stop("par$k_o_mis must be > 0.")
  p <- p0 * (1 - O2 / (O2 + k_o_mis))
  df <- data.frame(O2_pct = O2, p = pmax(p, 0))
  ggplot(df, aes(O2_pct, p)) +
    geom_line(linewidth = 1, color = "black") +
    geom_point(data = df[c(1, nrow(df)), ], size = 2, color = "red") +
    labs(x = "Oxygen (%)", y = "Missegregation rate") +
    theme_bw()
}
