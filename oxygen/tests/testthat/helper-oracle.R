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
        "O2_supply_demand_MAP",
        "model",
        "model_O2_supply_demand_MAP.R"
      ),
      mustWork = FALSE
    )
    if (file.exists(candidate)) {
      return(list(root = base_dir, model = candidate))
    }
  }
  stop("Cannot locate repo root/model_O2_supply_demand_MAP.R from start dir: ", start_dir)
}

repo_info <- find_repo_root(getwd())
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(repo_info$model))
source(
  file.path(dirname(dirname(repo_info$model)), "util", "o2_supply_demand_map_shared.R"),
  local = .GlobalEnv
)
source(repo_info$model, local = .GlobalEnv)

oracle_buffering_survival <- function(N, m_misseg, buffer_smax = 1.0,
                                      buffer_beta = 0.0, buffer_n_exp = 1.0,
                                      N_unit = 22L) {
  if (m_misseg <= 0L) return(1.0)
  if (N <= 0L) return(1.0)
  smax <- max(0.0, min(1.0, as.numeric(buffer_smax)))
  beta <- max(0.0, as.numeric(buffer_beta))
  n_exp <- max(0.0, as.numeric(buffer_n_exp))
  n_chr <- if (N_unit > 0L) as.numeric(N_unit) else 22.0
  ratio <- (2.0 * n_chr) / as.numeric(N)
  sN <- smax * exp(-beta * max(ratio, 0.0)^n_exp)
  sN <- max(0.0, min(1.0, sN))
  if (sN <= 0.0) return(0.0)
  max(0.0, min(1.0, sN^as.integer(m_misseg)))
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

oracle_live_mass_per_division <- function(N, p, N_unit = 22L, eps_tail = 0.0,
                                          buffer_smax = 1.0,
                                          buffer_beta = 0.0,
                                          buffer_n_exp = 1.0) {
  N_int <- as.integer(N)
  n_vals <- 0:N_int
  pn <- stats::dbinom(n_vals, size = N_int, prob = as.numeric(p))
  eps_use <- as.numeric(eps_tail)
  if (is.finite(eps_use) && eps_use > 0.0) {
    pn[pn < eps_use] <- 0.0
  }
  survival <- vapply(
    n_vals,
    function(n) oracle_buffering_survival(
      N_int,
      n,
      buffer_smax = buffer_smax,
      buffer_beta = buffer_beta,
      buffer_n_exp = buffer_n_exp,
      N_unit = N_unit
    ),
    numeric(1)
  )
  sum(pn * (2.0 * survival))
}
