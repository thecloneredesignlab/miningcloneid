required_pkgs <- c("testthat", "Matrix", "Rcpp", "dplyr", "ggplot2", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing required package(s) for O2 supply-demand tests: ", paste(missing_pkgs, collapse = ", "))
}

find_repo_root <- function(start_dir = getwd()) {
  cur <- normalizePath(start_dir, mustWork = FALSE)
  up_path <- function(base, levels) {
    if (levels <= 0L) return(base)
    normalizePath(file.path(base, paste(rep("..", levels), collapse = "/")), mustWork = FALSE)
  }
  for (i in 0:8) {
    base_dir <- up_path(cur, i)
    candidate <- normalizePath(
      file.path(
        base_dir,
        "oxygen",
        "code",
        "O2G_supply_demand_MAP",
        "model",
        "model_O2G_supply_demand_MAP.R"
      ),
      mustWork = FALSE
    )
    if (file.exists(candidate)) {
      return(list(root = base_dir, model = candidate))
    }
  }
  stop("Cannot locate repo root/model_O2G_supply_demand_MAP.R from start dir: ", start_dir)
}

repo_info <- find_repo_root(getwd())
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(repo_info$model))
source(
  file.path(dirname(dirname(repo_info$model)), "util", "o2g_supply_demand_map_shared.R"),
  local = .GlobalEnv
)
source(repo_info$model, local = .GlobalEnv)

oracle_balanced_hidden_copies <- function(N, N_unit = 22L) {
  N_int <- as.integer(N)
  U <- as.integer(N_unit)
  if (!is.finite(N_int) || N_int < 0L) stop("N must be a nonnegative integer.")
  if (!is.finite(U) || U <= 0L) stop("N_unit must be a positive integer.")
  q <- N_int %/% U
  r <- N_int %% U
  out <- rep.int(q, U)
  if (r > 0L) out[seq_len(r)] <- out[seq_len(r)] + 1L
  out
}

oracle_prob_no_null_after_loss <- function(N, m_loss, N_unit = 22L) {
  N_int <- as.integer(N)
  m <- as.integer(m_loss)
  if (m <= 0L) return(1.0)
  if (m > N_int) return(0.0)

  copies <- oracle_balanced_hidden_copies(N_int, N_unit = N_unit)
  if (any(copies <= 0L)) return(0.0)

  dp <- numeric(m + 1L)
  dp[1L] <- 1.0
  for (a in copies) {
    max_take <- min(a - 1L, m)
    ways <- choose(a, 0:max_take)
    next_dp <- numeric(m + 1L)
    for (k in 0:m) {
      base <- dp[k + 1L]
      if (base == 0.0) next
      jmax <- min(max_take, m - k)
      idx <- 0:jmax
      next_dp[k + idx + 1L] <- next_dp[k + idx + 1L] + base * ways[idx + 1L]
    }
    dp <- next_dp
  }

  safe_ways <- dp[m + 1L]
  total_ways <- choose(N_int, m)
  if (!is.finite(total_ways) || total_ways <= 0.0) return(0.0)
  safe <- safe_ways / total_ways
  max(0.0, min(1.0, safe))
}

oracle_Rnull <- function(N, m_loss, N_unit = 22L) {
  if (m_loss <= 0L) return(0.0)
  1.0 - oracle_prob_no_null_after_loss(N, m_loss, N_unit = N_unit)
}

oracle_Sloss <- function(N, m_loss, gamma_loss = 0.1, N_unit = 22L) {
  if (m_loss <= 0L) return(1.0)
  gamma_use <- as.numeric(gamma_loss)
  if (!is.finite(gamma_use) || gamma_use < 0.0) stop("gamma_loss must be finite and >= 0.")
  if (gamma_use == 0.0) return(1.0)
  safe <- 1.0 - oracle_Rnull(N, m_loss, N_unit = N_unit)
  safe <- max(0.0, min(1.0, safe))
  safe^gamma_use
}

triplet_to_sparse <- function(tri) {
  Matrix::sparseMatrix(
    i = as.integer(tri$i),
    j = as.integer(tri$j),
    x = as.numeric(tri$x),
    dims = c(as.integer(tri$nrow), as.integer(tri$ncol)),
    repr = "C"
  )
}

shift_weight <- function(delta_res, t) {
  idx <- which(as.integer(delta_res$ts) == as.integer(t))
  if (length(idx) == 0L) return(0.0)
  sum(as.numeric(delta_res$prob[idx]))
}

oracle_live_mass_per_division <- function(N, p, gamma_loss = 0.1, N_unit = 22L, eps_tail = 0.0) {
  N_int <- as.integer(N)
  n_vals <- 0:N_int
  pn <- stats::dbinom(n_vals, size = N_int, prob = as.numeric(p))
  eps_use <- as.numeric(eps_tail)
  if (is.finite(eps_use) && eps_use > 0.0) {
    pn[pn < eps_use] <- 0.0
  }
  loss_survival <- vapply(
    n_vals,
    function(n) oracle_Sloss(N_int, n, gamma_loss = gamma_loss, N_unit = N_unit),
    numeric(1)
  )
  sum(pn * (1.0 + loss_survival))
}
