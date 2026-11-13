# Shared natural/free parameter mapping and validation for diagnostic profiles.
# This module does not source a simulator or implement an optimizer.

o2prof_symbols <- c(
  "lam_max", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta", "buffer_n_exp",
  "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu", "O2_crit", "n_O",
  "p_mis_base", "sigma_growth", "sigma_kary", "init_mean_2N", "init_sd_2N",
  "init_mean_4N", "init_sd_4N")

o2prof_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

o2prof_hash <- function(path) unname(digest::digest(file = path, algo = "sha256"))

o2prof_map <- function(parameter_table, base_params, domain = c("production", "zero_extended"),
                       fixed = list(), penalty_active = character(),
                       zero_scales = c(p_misseg = 0.015, buffer_beta = 1)) {
  domain <- match.arg(domain)
  tab <- if (is.data.frame(parameter_table)) parameter_table else read.csv(parameter_table, stringsAsFactors = FALSE)
  required <- c("param_symbol", "init_value", "lower_bound", "upper_bound", "use_invitro_fit")
  if (!all(required %in% names(tab))) stop("Incomplete production parameter table")
  flagged <- tolower(trimws(as.character(tab$use_invitro_fit))) %in% c("true", "t", "1", "yes", "y")
  if (anyDuplicated(tab$param_symbol) || !setequal(tab$param_symbol[flagged], o2prof_symbols)) stop("Unexpected production fitted-parameter specification")
  if (!all(o2prof_symbols %in% names(base_params))) stop("Missing full-precision baseline parameters")
  if (length(fixed) && (is.null(names(fixed)) || anyDuplicated(names(fixed)) || any(!names(fixed) %in% o2prof_symbols))) stop("Unknown or unnamed fixed parameter")
  if (any(!penalty_active %in% o2prof_symbols)) stop("Unknown penalty-active parameter")
  if (!setequal(names(zero_scales), c("p_misseg", "buffer_beta")) || any(!is.finite(zero_scales) | zero_scales <= 0)) stop("Zero-transform scales must be fixed positive constants")
  p <- tab[match(o2prof_symbols, tab$param_symbol), , drop = FALSE]
  rownames(p) <- NULL
  linear <- c("buffer_smax", "gamma_growth", "gamma_mu", "n_O", "init_mean_2N", "init_mean_4N")
  p$parameter <- p$param_symbol
  p$lower <- as.numeric(p$lower_bound)
  p$upper <- as.numeric(p$upper_bound)
  p$transform <- ifelse(p$parameter %in% linear, "identity", "log10")
  p$scale <- 1
  if (domain == "zero_extended") {
    hit <- match(names(zero_scales), p$parameter)
    p$lower[hit] <- 0
    p$transform[hit] <- "log1p"
    p$scale[hit] <- zero_scales
  }
  if (any(!is.finite(p$lower) | !is.finite(p$upper) | p$lower >= p$upper) || any(p$lower[p$transform == "log10"] <= 0)) stop("Invalid parameter domain")
  p$status <- "active"
  p$reason <- "fitted nuisance"
  p$data_active <- TRUE
  p$penalty_active <- p$parameter %in% penalty_active
  p$fixed_value <- NA_real_
  for (n in names(fixed)) {
    i <- match(n, p$parameter)
    value <- as.numeric(fixed[[n]])
    if (length(value) != 1L || !is.finite(value) || value < p$lower[i] || value > p$upper[i]) stop("Fixed value outside declared domain: ", n)
    p$status[i] <- "fixed"
    p$fixed_value[i] <- value
    p$reason[i] <- "exact target constraint"
  }
  inactive <- character()
  if ("p_misseg" %in% names(fixed) && fixed$p_misseg == 0) inactive <- c(inactive, k_o_mis = "zero inducible amplitude")
  if ("buffer_beta" %in% names(fixed) && fixed$buffer_beta == 0) inactive <- c(inactive, buffer_n_exp = "zero ploidy buffering strength")
  for (n in names(inactive)) {
    i <- match(n, p$parameter)
    p$data_active[i] <- FALSE
    if (p$status[i] == "active" && !p$penalty_active[i]) {
      p$status[i] <- "inactive"
      p$fixed_value[i] <- as.numeric(base_params[[n]])
      p$reason[i] <- paste(inactive[[n]], "and no complete-objective contribution")
    } else if (p$penalty_active[i]) p$reason[i] <- "data-inactive but retained for objective penalty"
  }
  prefix <- ifelse(p$transform == "identity", "", paste0(p$transform, "_"))
  p$coordinate <- paste0(prefix, p$parameter)
  map <- list(table = p[, c("parameter", "status", "reason", "data_active", "penalty_active", "lower", "upper", "transform", "scale", "coordinate", "fixed_value")],
              base_params = base_params, domain = domain, fixed = fixed, zero_scales = zero_scales)
  class(map) <- "o2prof_parameter_map"
  map
}

o2prof_transform <- function(x, table, inverse = FALSE) {
  out <- as.numeric(x)
  log <- table$transform == "log10"
  zero <- table$transform == "log1p"
  if (inverse) {
    out[log] <- 10^out[log]
    out[zero] <- table$scale[zero] * expm1(out[zero])
  } else {
    if (any(out[log] <= 0) || any(out[zero] < 0)) stop("Invalid natural value for transform")
    out[log] <- log10(out[log])
    out[zero] <- log1p(out[zero] / table$scale[zero])
  }
  out
}

o2prof_project_start <- function(map, candidate, clip = FALSE) {
  p <- map$table
  if (!all(p$parameter %in% names(candidate))) stop("Start must contain all fitted natural parameters")
  values <- vapply(p$parameter, function(n) as.numeric(candidate[[n]]), numeric(1))
  if (any(!is.finite(values))) stop("Nonfinite start")
  before <- values
  constrained <- p$status != "active"
  values[constrained] <- p$fixed_value[constrained]
  outside <- values < p$lower | values > p$upper
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(p$lower), abs(p$upper))
  materially_outside <- values < p$lower - tolerance | values > p$upper + tolerance
  if (any(materially_outside) && !clip) stop("Start outside declared bounds")
  if (any(outside)) values <- pmin(pmax(values, p$lower), p$upper)
  full <- map$base_params
  for (i in seq_along(values)) full[[p$parameter[i]]] <- unname(values[i])
  reason <- ifelse(constrained, p$status, ifelse(materially_outside, "projected_to_bound", ifelse(outside, "roundoff_to_bound", "unchanged")))
  list(params = full, projection = data.frame(parameter = p$parameter, provided_value = before, actual_value = values, reason = reason))
}

o2prof_to_free <- function(map, params) {
  p <- map$table
  v <- vapply(p$parameter, function(n) as.numeric(params[[n]]), numeric(1))
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(p$lower), abs(p$upper))
  if (length(v) != nrow(p) || any(!is.finite(v)) || any(v < p$lower - tolerance | v > p$upper + tolerance)) stop("Infeasible full parameter vector")
  fixed <- p$status != "active"
  if (any(v[fixed] != p$fixed_value[fixed])) stop("Full vector disagrees with fixed/inactive values")
  active <- !fixed
  v[active] <- pmin(pmax(v[active], p$lower[active]), p$upper[active])
  setNames(o2prof_transform(v[active], p[active, , drop = FALSE]), p$coordinate[active])
}

o2prof_from_free <- function(map, free) {
  p <- map$table
  active <- p$status == "active"
  wanted <- p$coordinate[active]
  if (length(free) != length(wanted) || (length(names(free)) && !identical(names(free), wanted))) stop("Wrong free-vector coordinates")
  v <- p$fixed_value
  v[active] <- o2prof_transform(free, p[active, , drop = FALSE], inverse = TRUE)
  tolerance <- 64 * .Machine$double.eps * pmax(1, abs(p$lower), abs(p$upper))
  if (any(!is.finite(v)) || any(v < p$lower - tolerance | v > p$upper + tolerance)) stop("Free vector outside declared bounds")
  # Roundoff at transformed endpoints is projected only onto those exact bounds.
  v[active] <- pmin(pmax(v[active], p$lower[active]), p$upper[active])
  full <- map$base_params
  for (i in seq_along(v)) full[[p$parameter[i]]] <- unname(v[i])
  full
}

o2prof_spec <- function(map, initial_params) {
  p <- map$table[map$table$status == "active", , drop = FALSE]
  list(param_name = p$coordinate, lower = o2prof_transform(p$lower, p),
       upper = o2prof_transform(p$upper, p), init = unname(o2prof_to_free(map, initial_params)))
}

o2prof_initial_population <- function(map, NP, starts, provided = NULL) {
  if (!length(starts) || is.null(names(starts)) || any(!nzchar(names(starts))) || anyDuplicated(names(starts))) stop("Provide uniquely named start vectors")
  projected <- lapply(starts, function(x) o2prof_project_start(map, x))
  vectors <- lapply(projected, function(x) o2prof_to_free(map, x$params))
  spec <- o2prof_spec(map, projected[[1]]$params)
  d <- length(spec$param_name)
  if (NP < length(starts) || NP < 5L || d < 1L) stop("Invalid initial population size")
  if (is.null(provided)) {
    population <- matrix(stats::runif(NP * d), nrow = NP, ncol = d)
    population <- sweep(sweep(population, 2, spec$upper - spec$lower, "*"), 2, spec$lower, "+")
    for (i in seq_along(vectors)) population[i, ] <- vectors[[i]]
  } else {
    population <- as.matrix(provided)
    if (!identical(dim(population), c(as.integer(NP), as.integer(d)))) stop("Wrong provided population dimensions")
    if (!identical(colnames(population), spec$param_name)) stop("Wrong provided population coordinates")
    for (i in seq_along(vectors)) if (!identical(as.numeric(population[i, ]), unname(vectors[[i]]))) stop("Provided population disagrees with declared injected starts")
  }
  colnames(population) <- spec$param_name
  if (any(!is.finite(population)) || any(sweep(population, 2, spec$lower, "<")) || any(sweep(population, 2, spec$upper, ">"))) stop("Infeasible initial population")
  list(population = population, starts = projected, labels = names(starts))
}

o2prof_check_samples <- function(comp, manifest) {
  expected <- list(growth_df = manifest$observation_id[manifest$growth_scored],
                   ploidy_df = manifest$observation_id[manifest$karyotype_scored],
                   flow_df = manifest$observation_id[manifest$flow_scored])
  for (n in names(expected)) {
    tab <- comp[[n]]
    column <- if (n == "growth_df") "loglik" else "mean_loglik"
    if (is.null(tab) || !all(c("passage_id", column) %in% names(tab)) ||
        anyDuplicated(tab$passage_id) || !setequal(as.character(tab$passage_id), expected[[n]]) ||
        any(!is.finite(tab[[column]]))) stop("Retained sample contract violated: ", n)
  }
  if (!is.finite(comp$objective) || comp$n_growth_missing_pred != 0L) stop("Nonfinite objective or missing prediction")
  invisible(TRUE)
}

o2prof_accept_candidate <- function(candidate, expected_fingerprint, map, initial_population_hash,
                                    score_tolerance = 1e-8) {
  if (!identical(candidate$fingerprint, expected_fingerprint)) stop("Stale or wrong-run fingerprint")
  if (!isTRUE(candidate$launch_completed) || !identical(candidate$launch_exit_code, 0L)) stop("Interrupted or failed launch")
  if (!identical(candidate$initial_population_hash, initial_population_hash)) stop("Wrong actual starting population")
  o2prof_to_free(map, candidate$params)
  if (!isTRUE(candidate$samples_valid) || !is.finite(candidate$objective) || !is.finite(candidate$reevaluated_objective) ||
      abs(candidate$objective - candidate$reevaluated_objective) > score_tolerance) stop("Failed independent candidate validation")
  list(accepted = TRUE, search_status = if (isTRUE(candidate$budget_limited)) "budget_limited_provisional" else "completed_search_candidate")
}

o2prof_assert_resume <- function(directory, fingerprint) {
  path <- file.path(directory, "run_manifest.json")
  if (!file.exists(path)) stop("No compatible checkpoint manifest")
  old <- jsonlite::read_json(path, simplifyVector = TRUE)
  if (!identical(old$fingerprint, fingerprint)) stop("Resume fingerprint mismatch")
  if (file.exists(file.path(directory, "validated_result.json"))) stop("Completed result must not be overwritten")
  invisible(TRUE)
}
