#!/usr/bin/env Rscript

o2ipa_null_coalesce <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

o2ipa_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!grepl("^--", arg)) {
      i <- i + 1L
      next
    }
    kv <- sub("^--", "", arg)
    eq <- regexpr("=", kv, fixed = TRUE)
    if (eq > 0L) {
      key <- substr(kv, 1L, eq - 1L)
      val <- substr(kv, eq + 1L, nchar(kv))
      out[[key]] <- val
      i <- i + 1L
    } else {
      key <- kv
      if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      } else {
        out[[key]] <- TRUE
        i <- i + 1L
      }
    }
  }
  out
}

o2ipa_as_chr <- function(x, default = "") {
  val <- o2ipa_null_coalesce(x, default)
  val <- as.character(val[[1]])
  if (!nzchar(val)) default else val
}

o2ipa_as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric(o2ipa_null_coalesce(x, default)[[1]]))
  if (is.finite(val)) val else default
}

o2ipa_as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer(o2ipa_null_coalesce(x, default)[[1]]))
  if (!is.na(val)) val else default
}

o2ipa_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(isTRUE(default))
  if (is.logical(x[[1]])) return(isTRUE(x[[1]]))
  tolower(trimws(as.character(x[[1]]))) %in% c("1", "true", "t", "yes", "y", "on")
}

o2ipa_split_csv <- function(x, default = character()) {
  txt <- trimws(o2ipa_as_chr(x, paste(default, collapse = ",")))
  if (!nzchar(txt)) return(default)
  vals <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}

o2ipa_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
  }
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
}

o2ipa_find_workflow_root <- function(path = o2ipa_script_dir()) {
  cur <- normalizePath(path, mustWork = FALSE)
  if (file.exists(cur) && !dir.exists(cur)) cur <- dirname(cur)
  for (i in seq_len(8L)) {
    if (
      file.exists(file.path(cur, "util", "o2_supply_demand_map_shared.R")) &&
        file.exists(file.path(cur, "model", "model_O2_supply_demand_MAP.R"))
    ) {
      return(normalizePath(cur, mustWork = FALSE))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(file.path(path, "..", ".."), mustWork = FALSE)
}

o2ipa_workflow_root <- function(script_dir = o2ipa_script_dir()) {
  o2ipa_find_workflow_root(script_dir)
}

o2ipa_repo_root <- function(script_dir = o2ipa_script_dir()) {
  normalizePath(file.path(o2ipa_workflow_root(script_dir), "..", "..", ".."), mustWork = FALSE)
}

o2ipa_default_out_dir <- function(script_dir = o2ipa_script_dir()) {
  file.path(o2ipa_repo_root(script_dir), "oxygen", "results", "analysis")
}

o2ipa_mkdirs <- function(out_dir) {
  dirs <- file.path(out_dir, c("tables", "figures", "cache", "logs", "report"))
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

o2ipa_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

o2ipa_read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

o2ipa_read_csv_or_tsv <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  if (length(first) && grepl(",", first, fixed = TRUE) && !grepl("\t", first, fixed = TRUE)) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    o2ipa_read_tsv(path)
  }
}

o2ipa_seed_number <- function(seed_id) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

o2ipa_norm_seed <- function(x) {
  s <- as.character(x)
  s <- trimws(s)
  out <- ifelse(grepl("^seed[0-9]+$", s), s, ifelse(grepl("^[0-9]+$", s), paste0("seed", s), s))
  out
}

o2ipa_order_seeds <- function(seed_id) {
  n <- o2ipa_seed_number(seed_id)
  order(ifelse(is.na(n), Inf, n), seed_id)
}

o2ipa_target_params <- function() {
  c(
    "O2_crit", "alpha_o2", "gamma_growth", "mu_hp", "p_misseg", "k_o_mis",
    "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O", "lam_max",
    "p_mis_base", "p_wgd", "gamma_mu", "o2_S0", "kappa_O", "eta_o2", "rho_2N"
  )
}

o2ipa_param_module <- function(parameter) {
  map <- c(
    O2_crit = "hypoxia_sensing", alpha_o2 = "proliferation", gamma_growth = "proliferation",
    mu_hp = "death", p_misseg = "CIN_missegregation", k_o_mis = "CIN_missegregation",
    buffer_smax = "aneuploidy_buffering", buffer_beta = "aneuploidy_buffering",
    buffer_n_exp = "aneuploidy_buffering", n_O = "hypoxia_sensing", lam_max = "proliferation",
    p_mis_base = "CIN_missegregation", p_wgd = "WGD", gamma_mu = "death",
    o2_S0 = "O2_supply_demand", kappa_O = "O2_supply_demand",
    eta_o2 = "O2_supply_demand", rho_2N = "O2_supply_demand"
  )
  unname(map[parameter])
}

o2ipa_param_aliases <- function() {
  list(
    O2_crit = c("O2_crit", "o2_crit", "O2crit", "o2crit"),
    alpha_o2 = c("alpha_o2", "alpha_O2"),
    gamma_growth = c("gamma_growth"),
    mu_hp = c("mu_hp"),
    p_misseg = c("p_misseg", "p_mis", "p_missegregation"),
    k_o_mis = c("k_o_mis", "ko_mis", "k_O_mis"),
    buffer_smax = c("buffer_smax", "buffer_s_max"),
    buffer_beta = c("buffer_beta"),
    buffer_n_exp = c("buffer_n_exp", "buffer_n"),
    n_O = c("n_O", "n_o", "nO"),
    lam_max = c("lam_max", "lambda_max"),
    p_mis_base = c("p_mis_base", "p_misseg_base"),
    p_wgd = c("p_wgd", "wgd_prob"),
    gamma_mu = c("gamma_mu"),
    o2_S0 = c("o2_S0", "O2_S0", "S0_o2"),
    kappa_O = c("kappa_O", "kappa_o"),
    eta_o2 = c("eta_o2", "eta_O2"),
    rho_2N = c("rho_2N", "rho2N")
  )
}

o2ipa_optimizer_transform <- function(parameter) {
  log10_params <- c(
    "O2_crit", "alpha_o2", "mu_hp", "p_misseg", "k_o_mis",
    "buffer_beta", "buffer_n_exp", "lam_max", "p_mis_base", "p_wgd",
    "o2_S0", "kappa_O", "eta_o2", "rho_2N"
  )
  if (parameter %in% log10_params) "log10" else "identity"
}

o2ipa_transform_parameter_value <- function(parameter, value, epsilon = 1e-12) {
  v <- as.numeric(value)
  if (!is.finite(v)) return(NA_real_)
  tr <- o2ipa_optimizer_transform(parameter)
  if (identical(tr, "log10")) {
    if (v <= 0) return(NA_real_)
    return(log10(v))
  }
  if (identical(tr, "logit")) {
    vv <- min(max(v, epsilon), 1 - epsilon)
    return(stats::qlogis(vv))
  }
  v
}

o2ipa_transform_metadata <- function(parameter_tables = list()) {
  params <- o2ipa_target_params()
  out <- data.frame(
    parameter = params,
    module = vapply(params, o2ipa_param_module, character(1)),
    raw_scale = ifelse(params %in% c("p_misseg", "p_mis_base", "p_wgd", "buffer_smax"), "bounded_or_probability", "positive_or_identity"),
    transform = vapply(params, o2ipa_optimizer_transform, character(1)),
    epsilon = ifelse(params %in% c("p_misseg", "p_mis_base", "p_wgd", "buffer_smax"), 1e-12, NA_real_),
    optimizer_scale = vapply(params, o2ipa_optimizer_transform, character(1)),
    source_file = NA_character_,
    stringsAsFactors = FALSE
  )
  for (pt in parameter_tables) {
    if (!is.data.frame(pt) || !all(c("param_prototype", "param_name") %in% names(pt))) next
    idx <- match(out$parameter, pt$param_prototype)
    hit <- which(!is.na(idx))
    if (length(hit)) {
      out$source_file[hit] <- pt$source_file[idx[hit]]
      out$optimizer_scale[hit] <- ifelse(grepl("^log10_", pt$param_name[idx[hit]]), "log10", out$optimizer_scale[hit])
    }
  }
  out
}

o2ipa_find_extra <- function(run_dir, file) {
  path <- file.path(run_dir, "extra_results", file)
  if (file.exists(path)) path else NA_character_
}

o2ipa_read_extra_summary <- function(run_dir) {
  path <- o2ipa_find_extra(run_dir, "seed_summary.tsv")
  if (is.na(path)) return(NULL)
  tab <- o2ipa_read_tsv(path)
  if (!"seed" %in% names(tab)) return(NULL)
  tab$seed_id <- o2ipa_norm_seed(tab$seed)
  tab
}

o2ipa_read_param_matrix <- function(run_dir) {
  parent <- dirname(normalizePath(run_dir, mustWork = FALSE))
  base <- basename(normalizePath(run_dir, mustWork = FALSE))
  candidates <- c(
    file.path(run_dir, "param_matrix_with_seed.tsv"),
    file.path(run_dir, "parameter_matrix_with_seed.tsv"),
    file.path(parent, paste0(base, "_param_matrix_with_seed.tsv")),
    file.path(parent, paste0(base, "_param_matrix.tsv"))
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(NULL)
  tab <- o2ipa_read_tsv(candidates[[1]])
  if (!"seed" %in% names(tab)) tab$seed <- seq_len(nrow(tab))
  tab$seed_id <- o2ipa_norm_seed(tab$seed)
  attr(tab, "source_file") <- candidates[[1]]
  tab
}

o2ipa_discover_seeds <- function(run_dir) {
  dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("^seed[^/]*[0-9]+$", basename(dirs))]
  seed_ids <- o2ipa_norm_seed(basename(dirs))
  seed_dir_map <- setNames(normalizePath(dirs, mustWork = FALSE), seed_ids)

  summary_tab <- o2ipa_read_extra_summary(run_dir)
  matrix_tab <- o2ipa_read_param_matrix(run_dir)
  all_ids <- unique(c(seed_ids, if (!is.null(summary_tab)) summary_tab$seed_id, if (!is.null(matrix_tab)) matrix_tab$seed_id))
  all_ids <- all_ids[o2ipa_order_seeds(all_ids)]
  data.frame(
    seed_id = all_ids,
    seed_dir = unname(seed_dir_map[all_ids]),
    stringsAsFactors = FALSE
  )
}

o2ipa_metric_map <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(list())
  tab <- o2ipa_read_tsv(path)
  if (!all(c("metric", "value") %in% names(tab))) return(list())
  vals <- as.list(tab$value)
  names(vals) <- tab$metric
  vals
}

o2ipa_num_from_map <- function(map, key) {
  if (is.null(map[[key]])) return(NA_real_)
  suppressWarnings(as.numeric(map[[key]]))
}

o2ipa_chr_from_map <- function(map, key) {
  if (is.null(map[[key]])) return(NA_character_)
  as.character(map[[key]])
}

o2ipa_read_best_params <- function(path) {
  if (!file.exists(path)) return(list(values = setNames(numeric(0), character(0)), aliases = data.frame()))
  tab <- o2ipa_read_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain parameter and value columns: ", path)
  }
  vals <- suppressWarnings(as.numeric(tab$value))
  names(vals) <- trimws(as.character(tab$parameter))
  list(values = vals, aliases = data.frame())
}

o2ipa_extract_param <- function(seed_id, parameter, best_vals, summary_row, matrix_row) {
  aliases <- o2ipa_param_aliases()[[parameter]]
  sources <- list()
  if (length(best_vals)) {
    sources$best_params <- best_vals
  }
  if (!is.null(summary_row) && nrow(summary_row) == 1L) {
    vals <- numeric(0)
    for (a in aliases) {
      col <- paste0("value__", a)
      if (col %in% names(summary_row)) vals[a] <- suppressWarnings(as.numeric(summary_row[[col]][[1]]))
    }
    sources$seed_summary_value_cols <- vals
  }
  if (!is.null(matrix_row) && nrow(matrix_row) == 1L) {
    vals <- numeric(0)
    for (a in aliases) {
      for (col in c(paste0("final_", a), a)) {
        if (col %in% names(matrix_row)) vals[a] <- suppressWarnings(as.numeric(matrix_row[[col]][[1]]))
      }
    }
    sources$param_matrix <- vals
  }

  for (src in names(sources)) {
    vals <- sources[[src]]
    for (a in aliases) {
      if (a %in% names(vals)) {
        v <- suppressWarnings(as.numeric(vals[[a]]))
        if (is.finite(v)) {
          return(list(value = v, source = src, alias = a))
        }
      }
    }
  }
  list(value = NA_real_, source = NA_character_, alias = NA_character_)
}

o2ipa_extract_all_params <- function(seed_id, seed_dir, summary_tab, matrix_tab) {
  best_path <- if (!is.na(seed_dir) && nzchar(seed_dir)) file.path(seed_dir, "best_params.tsv") else NA_character_
  best_vals <- if (!is.na(best_path) && file.exists(best_path)) o2ipa_read_best_params(best_path)$values else setNames(numeric(0), character(0))
  summary_row <- if (!is.null(summary_tab)) summary_tab[summary_tab$seed_id == seed_id, , drop = FALSE] else NULL
  if (!is.null(summary_row) && nrow(summary_row) > 1L) summary_row <- summary_row[1, , drop = FALSE]
  matrix_row <- if (!is.null(matrix_tab)) matrix_tab[matrix_tab$seed_id == seed_id, , drop = FALSE] else NULL
  if (!is.null(matrix_row) && nrow(matrix_row) > 1L) matrix_row <- matrix_row[1, , drop = FALSE]
  params <- o2ipa_target_params()
  rows <- lapply(params, function(p) {
    hit <- o2ipa_extract_param(seed_id, p, best_vals, summary_row, matrix_row)
    data.frame(
      seed_id = seed_id,
      parameter = p,
      value = hit$value,
      parameter_source = hit$source,
      matched_alias = hit$alias,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

o2ipa_choose_objective <- function(summary_vals, objective_source = "auto") {
  raw <- NA_real_
  b <- o2ipa_num_from_map(summary_vals, "objective_burden_neg2loglik_raw")
  p <- o2ipa_num_from_map(summary_vals, "objective_ploidy_neg2loglik_raw")
  if (is.finite(b) && is.finite(p)) raw <- b + p
  data <- o2ipa_num_from_map(summary_vals, "objective_data")
  total <- o2ipa_num_from_map(summary_vals, "objective_total")
  if (!is.finite(total)) total <- o2ipa_num_from_map(summary_vals, "objective")
  source <- match.arg(objective_source, c("auto", "raw_likelihood", "data", "total"))
  if (identical(source, "raw_likelihood")) return(list(value = raw, source = "raw_likelihood"))
  if (identical(source, "data")) return(list(value = data, source = "objective_data"))
  if (identical(source, "total")) return(list(value = total, source = if (is.finite(o2ipa_num_from_map(summary_vals, "objective_total"))) "objective_total" else "objective"))
  if (is.finite(raw)) return(list(value = raw, source = "raw_likelihood"))
  if (is.finite(data)) return(list(value = data, source = "objective_data"))
  list(value = total, source = if (is.finite(total)) "objective_total_or_objective" else NA_character_)
}

o2ipa_collect_seed_inputs <- function(run_dir, objective_source = "auto") {
  manifest0 <- o2ipa_discover_seeds(run_dir)
  summary_tab <- o2ipa_read_extra_summary(run_dir)
  matrix_tab <- o2ipa_read_param_matrix(run_dir)
  boundary_path <- o2ipa_find_extra(run_dir, "parameter_boundary_long.tsv")
  boundary_long <- if (!is.na(boundary_path)) o2ipa_read_tsv(boundary_path) else NULL
  if (!is.null(boundary_long) && "seed" %in% names(boundary_long)) {
    boundary_long$seed_id <- o2ipa_norm_seed(boundary_long$seed)
  }

  param_rows <- list()
  manifest_rows <- vector("list", nrow(manifest0))
  for (i in seq_len(nrow(manifest0))) {
    seed_id <- manifest0$seed_id[[i]]
    seed_dir <- manifest0$seed_dir[[i]]
    if (is.na(seed_dir)) seed_dir <- ""
    fit_summary_path <- if (nzchar(seed_dir)) file.path(seed_dir, "fit_summary.tsv") else NA_character_
    summary_vals <- if (!is.na(fit_summary_path) && file.exists(fit_summary_path)) {
      o2ipa_metric_map(fit_summary_path)
    } else {
      row <- if (!is.null(summary_tab)) summary_tab[summary_tab$seed_id == seed_id, , drop = FALSE] else NULL
      if (!is.null(row) && nrow(row)) as.list(row[1, , drop = TRUE]) else list()
    }
    obj <- o2ipa_choose_objective(summary_vals, objective_source = objective_source)
    params_long <- o2ipa_extract_all_params(seed_id, seed_dir, summary_tab, matrix_tab)
    param_rows[[i]] <- params_long
    missing_params <- params_long$parameter[!is.finite(params_long$value)]

    bseed <- if (!is.null(boundary_long)) boundary_long[boundary_long$seed_id == seed_id, , drop = FALSE] else NULL
    target_boundary <- if (!is.null(bseed) && nrow(bseed)) {
      pname <- if ("param_prototype" %in% names(bseed)) bseed$param_prototype else if ("parameter" %in% names(bseed)) bseed$parameter else bseed$param_name
      bseed[pname %in% o2ipa_target_params(), , drop = FALSE]
    } else {
      data.frame()
    }
    boundary_risk <- FALSE
    n_near <- 0L
    if (nrow(target_boundary)) {
      status_col <- if ("bound_status" %in% names(target_boundary)) "bound_status" else if ("status" %in% names(target_boundary)) "status" else NA_character_
      if (!is.na(status_col)) {
        boundary_risk <- any(!is.na(target_boundary[[status_col]]) & !target_boundary[[status_col]] %in% c("interior", ""))
        n_near <- sum(!is.na(target_boundary[[status_col]]) & !target_boundary[[status_col]] %in% c("interior", ""))
      }
    }
    fit_success <- length(missing_params) == 0L && is.finite(obj$value)
    convergence <- o2ipa_chr_from_map(summary_vals, "deoptim_stop_reason")
    if (is.na(convergence)) convergence <- o2ipa_chr_from_map(summary_vals, "optimizer_local_convergence")
    failure_parts <- character(0)
    if (!is.finite(obj$value)) failure_parts <- c(failure_parts, "missing_objective")
    if (length(missing_params)) failure_parts <- c(failure_parts, paste0("missing_params:", paste(missing_params, collapse = ",")))
    manifest_rows[[i]] <- data.frame(
      seed_id = seed_id,
      seed_dir = if (nzchar(seed_dir)) normalizePath(seed_dir, mustWork = FALSE) else NA_character_,
      fit_success = fit_success,
      convergence_status = convergence,
      objective = obj$value,
      objective_source = obj$source,
      objective_total = o2ipa_num_from_map(summary_vals, "objective_total"),
      objective_data = o2ipa_num_from_map(summary_vals, "objective_data"),
      objective_burden = o2ipa_num_from_map(summary_vals, "objective_burden"),
      objective_ploidy = o2ipa_num_from_map(summary_vals, "objective_ploidy"),
      burden_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_burden_neg2loglik_raw"),
      ploidy_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_ploidy_neg2loglik_raw"),
      runtime = o2ipa_num_from_map(summary_vals, "runtime"),
      parameter_file = if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir, "best_params.tsv"))) file.path(seed_dir, "best_params.tsv") else NA_character_,
      config_file = if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir, "fit_config.rds"))) file.path(seed_dir, "fit_config.rds") else NA_character_,
      visualization_available = !is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir, "viz_status.log")),
      failure_reason = if (fit_success) NA_character_ else paste(failure_parts, collapse = ";"),
      boundary_risk = boundary_risk,
      number_of_parameters_near_boundary = n_near,
      stringsAsFactors = FALSE
    )
  }
  manifest <- do.call(rbind, manifest_rows)
  params_long <- do.call(rbind, param_rows)
  list(manifest = manifest, params_long = params_long, boundary_long = boundary_long, summary_tab = summary_tab, matrix_tab = matrix_tab)
}

o2ipa_params_wide <- function(params_long, value_col = "value") {
  seeds <- unique(params_long$seed_id)
  params <- o2ipa_target_params()
  mat <- matrix(NA_real_, nrow = length(seeds), ncol = length(params), dimnames = list(seeds, params))
  for (i in seq_len(nrow(params_long))) {
    mat[params_long$seed_id[[i]], params_long$parameter[[i]]] <- params_long[[value_col]][[i]]
  }
  as.data.frame(mat, check.names = FALSE)
}

o2ipa_source_model <- function(script_dir = o2ipa_script_dir()) {
  workflow_root <- o2ipa_workflow_root(script_dir)
  model_file <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")
  if (!file.exists(model_file)) stop("Model file not found: ", model_file)
  shared_file <- file.path(workflow_root, "util", "o2_supply_demand_map_shared.R")
  common_file <- file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R")
  if (file.exists(shared_file)) source(shared_file, local = globalenv())
  if (file.exists(common_file)) source(common_file, local = globalenv())
  env <- new.env(parent = globalenv())
  old_dir <- Sys.getenv("MININGCLONEID_OXYGEN_CODE_DIR", unset = NA_character_)
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_file))
  on.exit({
    if (is.na(old_dir)) Sys.unsetenv("MININGCLONEID_OXYGEN_CODE_DIR") else Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = old_dir)
  }, add = TRUE)
  source(model_file, local = env)
  env
}

o2ipa_call_model <- function(model_env, fn, ...) {
  if (!exists(fn, envir = model_env, inherits = TRUE)) stop("Model helper missing: ", fn)
  get(fn, envir = model_env, inherits = TRUE)(...)
}

o2ipa_auc <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

o2ipa_crossing <- function(x, y, target, decreasing = FALSE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(list(value = NA_real_, status = "insufficient_grid"))
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  if (target < min(y) || target > max(y)) return(list(value = NA_real_, status = "no_crossing_in_grid"))
  if (isTRUE(decreasing)) {
    ord2 <- order(y)
    return(list(value = as.numeric(stats::approx(y[ord2], x[ord2], xout = target, ties = "ordered")$y), status = "crossed"))
  }
  list(value = as.numeric(stats::approx(y, x, xout = target, ties = "ordered")$y), status = "crossed")
}

o2ipa_feature_row <- function(seed_id, scope, module, feature, value, feature_type = "positive", source = "model_probe") {
  data.frame(
    seed_id = seed_id,
    fingerprint_scope = scope,
    module = module,
    feature = feature,
    raw_value = as.numeric(value),
    feature_type = feature_type,
    source = source,
    stringsAsFactors = FALSE
  )
}

o2ipa_feature_name_num <- function(prefix, value) {
  paste0(prefix, gsub("[^A-Za-z0-9]+", "p", format(value, trim = TRUE, scientific = FALSE)))
}

o2ipa_build_static_one <- function(seed_id, run_params, model_env, scope = "full") {
  rp <- as.list(run_params)
  names(rp) <- names(run_params)
  rp$ploidy_O2_death <- o2ipa_null_coalesce(rp$ploidy_O2_death, "ploidy_related")
  rp$O2_growth <- o2ipa_null_coalesce(rp$O2_growth, TRUE)
  rp$o2_min <- as.numeric(o2ipa_null_coalesce(rp$o2_min, 0))
  rp$o2_S0_upper_bound <- as.numeric(o2ipa_null_coalesce(rp$o2_S0_upper_bound, 5))
  rp$o2_Nref <- as.numeric(o2ipa_null_coalesce(rp$o2_Nref, 1e6))
  if (!is.finite(rp$o2_S0_upper_bound) || rp$o2_S0_upper_bound <= 0) rp$o2_S0_upper_bound <- 5
  o2_grid <- sort(unique(pmin(c(0.05, 0.1, 0.25, 0.5, 1, 2, 5), rp$o2_S0_upper_bound)))
  state_grid <- c(44, 66, 88, 132)
  state_grid <- state_grid[state_grid >= 22 & state_grid <= 154]
  p_mis_grid <- c(1e-4, 1e-3, 1e-2, 5e-2)
  demand_grid <- c(0, 1e3, 1e4, 1e5, 1e6, 1e7)
  burden_grid <- c(0.01, 0.1, 1, 10, 100)
  rows <- list()
  add <- function(module, feature, value, feature_type = "positive", source = "model_probe") {
    rows[[length(rows) + 1L]] <<- o2ipa_feature_row(seed_id, scope, module, feature, value, feature_type, source)
  }

  stress <- o2ipa_call_model(model_env, ".resource_stress_of_O2", O2 = o2_grid, run_params = rp)
  for (i in seq_along(o2_grid)) add("hypoxia_sensing", o2ipa_feature_name_num("stress_at_O2_", o2_grid[[i]]), stress[[i]], "fraction")
  add("hypoxia_sensing", "stress_AUC", o2ipa_auc(o2_grid, stress), "positive")
  dense_o2 <- seq(min(o2_grid), max(o2_grid), length.out = 400L)
  dense_stress <- o2ipa_call_model(model_env, ".resource_stress_of_O2", O2 = dense_o2, run_params = rp)
  for (thr in c(0.1, 0.5, 0.9)) {
    cr <- o2ipa_crossing(dense_o2, dense_stress, thr, decreasing = TRUE)
    add("hypoxia_sensing", paste0("O2_at_stress_", as.integer(thr * 100), "pct"), cr$value, "threshold")
    add("hypoxia_sensing", paste0("crossing_status_stress_", as.integer(thr * 100), "pct_crossed"), ifelse(cr$status == "crossed", 1, 0), "binary")
  }
  cr10 <- o2ipa_crossing(dense_o2, dense_stress, 0.1, decreasing = TRUE)
  cr90 <- o2ipa_crossing(dense_o2, dense_stress, 0.9, decreasing = TRUE)
  add("hypoxia_sensing", "stress_transition_width", if (is.finite(cr10$value) && is.finite(cr90$value)) abs(cr10$value - cr90$value) else NA_real_, "positive")

  grid <- expand.grid(O2 = o2_grid, N = state_grid)
  lam <- o2ipa_call_model(model_env, ".lambda_eff_of_O2", O2 = grid$O2, run_params = rp, N = grid$N, O2_growth = isTRUE(rp$O2_growth))
  mu <- o2ipa_call_model(model_env, ".mu_eff_of_O2", O2 = grid$O2, run_params = rp, N = grid$N)
  pm <- o2ipa_call_model(model_env, ".pmisseg_of_O2", O2 = grid$O2, run_params = rp, N = grid$N)
  net <- lam - mu
  for (i in seq_len(nrow(grid))) {
    suffix <- paste0("_O2_", gsub("[^0-9]+", "p", format(grid$O2[[i]], trim = TRUE, scientific = FALSE)), "_N_", grid$N[[i]])
    add("proliferation", paste0("lambda", suffix), lam[[i]], "positive")
    add("death", paste0("mu", suffix), mu[[i]], "positive")
    add("net_growth", paste0("net_growth", suffix), net[[i]], "signed")
    add("CIN_missegregation", paste0("p_misseg", suffix), pm[[i]], "probability")
  }
  for (n in state_grid) {
    idx <- grid$N == n
    add("proliferation", paste0("proliferation_AUC_over_O2_N_", n), o2ipa_auc(grid$O2[idx], lam[idx]), "positive")
    add("death", paste0("death_AUC_over_O2_N_", n), o2ipa_auc(grid$O2[idx], mu[idx]), "positive")
    add("net_growth", paste0("positive_net_growth_AUC_N_", n), o2ipa_auc(grid$O2[idx], pmax(net[idx], 0)), "positive")
    add("net_growth", paste0("negative_net_growth_AUC_N_", n), o2ipa_auc(grid$O2[idx], pmax(-net[idx], 0)), "positive")
    cr <- o2ipa_crossing(grid$O2[idx], net[idx], 0, decreasing = FALSE)
    add("net_growth", paste0("critical_O2_net_growth_zero_N_", n), cr$value, "threshold")
  }
  dip_hi <- lam[grid$N == 44 & grid$O2 == max(o2_grid)][[1]]
  dip_lo <- lam[grid$N == 44 & grid$O2 == min(o2_grid)][[1]]
  add("proliferation", "low_O2_growth_suppression_ratio_diploid", if (is.finite(dip_hi) && dip_hi != 0) dip_lo / dip_hi else NA_real_, "ratio")
  high_hi <- lam[grid$N == max(state_grid) & grid$O2 == max(o2_grid)][[1]]
  add("proliferation", "high_ploidy_vs_diploid_growth_ratio_high_O2", if (is.finite(dip_hi) && dip_hi != 0) high_hi / dip_hi else NA_real_, "ratio")
  mu_hi <- mu[grid$N == 44 & grid$O2 == max(o2_grid)][[1]]
  mu_lo <- mu[grid$N == 44 & grid$O2 == min(o2_grid)][[1]]
  add("death", "low_O2_death_induction_ratio_diploid", if (is.finite(mu_hi) && mu_hi != 0) mu_lo / mu_hi else NA_real_, "ratio")
  high_mu <- mu[grid$N == max(state_grid) & grid$O2 == max(o2_grid)][[1]]
  add("death", "high_ploidy_vs_diploid_death_ratio_high_O2", if (is.finite(mu_hi) && mu_hi != 0) high_mu / mu_hi else NA_real_, "ratio")
  pm_hi <- pm[grid$N == 44 & grid$O2 == max(o2_grid)][[1]]
  pm_lo <- pm[grid$N == 44 & grid$O2 == min(o2_grid)][[1]]
  add("CIN_missegregation", "hypoxia_induced_misseg_increment_diploid", pm_lo - pm_hi, "probability")
  add("CIN_missegregation", "low_O2_to_high_O2_p_misseg_ratio_diploid", if (is.finite(pm_hi) && pm_hi != 0) pm_lo / pm_hi else NA_real_, "ratio")
  add("CIN_missegregation", "expected_missegregation_events_per_100_divisions_low_O2_diploid", 100 * pm_lo, "positive")

  for (n in state_grid) {
    for (p in p_mis_grid) {
      delta <- o2ipa_call_model(model_env, ".pr_delta_vec",
        N = n, p = p, eps_tail = 1e-8, N_unit = 22L,
        buffer_smax = as.numeric(rp$buffer_smax),
        buffer_beta = as.numeric(rp$buffer_beta),
        buffer_n_exp = as.numeric(rp$buffer_n_exp)
      )
      shifts <- suppressWarnings(as.numeric(names(delta)))
      prob <- as.numeric(delta)
      md <- as.numeric(attr(delta, "mass_dropped"))
      suffix <- paste0("_N_", n, "_p_", gsub("[^0-9]+", "p", format(p, trim = TRUE, scientific = FALSE)))
      add("aneuploidy_buffering", paste0("surviving_transition_mass", suffix), sum(prob, na.rm = TRUE), "positive")
      add("aneuploidy_buffering", paste0("mass_dropped", suffix), md, "fraction")
      add("aneuploidy_buffering", paste0("mean_absolute_chromosome_shift_after_buffering", suffix), sum(abs(shifts) * prob, na.rm = TRUE), "positive")
      mu_shift <- sum(shifts * prob, na.rm = TRUE)
      add("aneuploidy_buffering", paste0("transition_variance", suffix), sum(((shifts - mu_shift)^2) * prob, na.rm = TRUE), "positive")
      add("aneuploidy_buffering", paste0("large_shift_survival_fraction", suffix), sum(prob[abs(shifts) >= 3], na.rm = TRUE), "fraction")
    }
  }

  pw <- as.numeric(rp$p_wgd)
  add("WGD", "p_wgd_per_division", pw, "probability")
  add("WGD", "expected_WGD_events_per_100_divisions", 100 * pw, "positive")
  add("WGD", "probability_at_least_one_WGD_in_100_divisions", 1 - (1 - min(max(pw, 0), 1))^100, "probability")

  o2_eff <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = demand_grid, run_params = rp, o2_Nref = as.numeric(rp$o2_Nref))
  for (i in seq_along(demand_grid)) add("O2_supply_demand", o2ipa_feature_name_num("O2_at_effective_demand_", demand_grid[[i]]), o2_eff[[i]], "oxygen")
  for (n in state_grid) {
    demand <- burden_grid * as.numeric(rp$rho_2N) * ((n / 44)^as.numeric(rp$eta_o2))
    o2_obs <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = demand, run_params = rp, o2_Nref = as.numeric(rp$o2_Nref))
    for (i in seq_along(burden_grid)) {
      add("O2_supply_demand", paste0(o2ipa_feature_name_num("O2_at_observed_burden_", burden_grid[[i]]), "_N_", n), o2_obs[[i]], "oxygen")
      add("O2_supply_demand", paste0(o2ipa_feature_name_num("effective_demand_at_observed_burden_", burden_grid[[i]]), "_N_", n), demand[[i]], "positive")
    }
    for (target_o2 in c(2, 1, 0.5)) {
      dense_b <- 10^seq(log10(min(burden_grid)), log10(max(burden_grid)), length.out = 200L)
      dense_d <- dense_b * as.numeric(rp$rho_2N) * ((n / 44)^as.numeric(rp$eta_o2))
      dense_o <- o2ipa_call_model(model_env, ".o2_supply_demand_from_burden", Ntot = dense_d, run_params = rp, o2_Nref = as.numeric(rp$o2_Nref))
      cr <- o2ipa_crossing(log10(dense_b), dense_o, target_o2, decreasing = TRUE)
      add("O2_supply_demand", paste0("log_burden_at_O2_", gsub("[^0-9]+", "p", target_o2), "pct_N_", n), cr$value, "signed")
    }
  }
  add("O2_supply_demand", "minimum_O2_floor", min(o2_eff, na.rm = TRUE), "oxygen")
  add("O2_supply_demand", "oxygen_depletion_AUC_effective_demand", o2ipa_auc(log10(demand_grid + 1), pmax(max(o2_eff, na.rm = TRUE) - o2_eff, 0)), "positive")
  do.call(rbind, rows)
}

o2ipa_feature_transform <- function(value, feature_type, oxygen_floor = 1e-6) {
  v <- as.numeric(value)
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v)
  ft <- as.character(feature_type)
  for (tp in unique(ft)) {
    idx <- ok & ft == tp
    if (!any(idx)) next
    if (tp %in% c("probability", "fraction", "binary")) {
      vv <- pmin(pmax(v[idx], 1e-12), 1 - 1e-12)
      out[idx] <- stats::qlogis(vv)
    } else if (tp %in% c("positive", "ratio")) {
      out[idx] <- ifelse(v[idx] > 0, log10(v[idx]), NA_real_)
    } else if (tp %in% c("time")) {
      out[idx] <- ifelse(v[idx] >= 0, log1p(v[idx]), NA_real_)
    } else if (tp %in% c("oxygen")) {
      out[idx] <- ifelse(v[idx] >= 0, log10(pmax(v[idx], oxygen_floor)), NA_real_)
    } else if (tp %in% c("threshold")) {
      out[idx] <- ifelse(v[idx] > 0, log10(v[idx]), NA_real_)
    } else {
      out[idx] <- v[idx]
    }
  }
  out
}

o2ipa_scale_long_features <- function(long_df, missing_feature_max = 0.5, missing_seed_max = 0.05, scale_floor = sqrt(.Machine$double.eps)) {
  if (!nrow(long_df)) return(list(long = long_df, wide = data.frame(), metadata = data.frame(), missingness = data.frame(), excluded_seeds = character()))
  long_df$transformed_value <- o2ipa_feature_transform(long_df$raw_value, long_df$feature_type)
  long_df$feature_key <- paste(long_df$fingerprint_scope, long_df$module, long_df$feature, sep = "||")
  meta_base <- unique(long_df[, c("feature_key", "fingerprint_scope", "module", "feature", "feature_type")])
  seeds <- unique(long_df$seed_id)
  trans_wide <- o2ipa_long_to_wide(long_df[, c("seed_id", "feature_key", "transformed_value")], "feature_key", "transformed_value")
  trans_mat <- as.matrix(trans_wide[, setdiff(names(trans_wide), "seed_id"), drop = FALSE])
  rownames(trans_mat) <- trans_wide$seed_id
  trans_mat <- trans_mat[seeds, meta_base$feature_key, drop = FALSE]
  stats_rows <- lapply(meta_base$feature_key, function(k) {
    vals <- trans_mat[, k]
    miss <- mean(!is.finite(vals))
    finite_vals <- vals[is.finite(vals)]
    med <- stats::median(finite_vals, na.rm = TRUE)
    mad <- stats::mad(finite_vals, center = med, constant = 1, na.rm = TRUE)
    if (!is.finite(med)) med <- NA_real_
    if (!is.finite(mad)) mad <- NA_real_
    data.frame(
      feature_key = k,
      n_seed = length(seeds),
      n_observed = sum(is.finite(vals)),
      n_missing = sum(!is.finite(vals)),
      missing_fraction = miss,
      center = med,
      scale = mad,
      zero_variance = !is.finite(mad) || mad <= scale_floor,
      stringsAsFactors = FALSE
    )
  })
  stat <- do.call(rbind, stats_rows)
  meta <- merge(meta_base, stat, by = "feature_key", all.x = TRUE)
  meta$scale_floor <- scale_floor
  meta$retained_for_clustering <- meta$missing_fraction <= missing_feature_max & !meta$zero_variance
  long_df <- merge(long_df, meta[, c("feature_key", "center", "scale", "retained_for_clustering")], by = "feature_key", all.x = TRUE)
  retained_keys <- meta$feature_key[meta$retained_for_clustering %in% TRUE]
  long_df$missing_before_imputation <- !is.finite(long_df$transformed_value)
  long_df$imputed_for_clustering <- long_df$retained_for_clustering %in% TRUE & long_df$missing_before_imputation
  long_df$clustering_transformed_value <- long_df$transformed_value
  imp_idx <- long_df$imputed_for_clustering %in% TRUE
  long_df$clustering_transformed_value[imp_idx] <- long_df$center[imp_idx]
  long_df$scaled_value <- with(long_df, ifelse(retained_for_clustering %in% TRUE, (clustering_transformed_value - center) / scale, NA_real_))
  retained <- long_df[long_df$feature_key %in% retained_keys, , drop = FALSE]
  if (length(retained_keys)) {
    retained_raw_mat <- trans_mat[, retained_keys, drop = FALSE]
    seed_miss <- data.frame(
      seed_id = rownames(retained_raw_mat),
      missing_fraction_retained_features = rowMeans(!is.finite(retained_raw_mat)),
      stringsAsFactors = FALSE
    )
  } else {
    seed_miss <- data.frame(seed_id = seeds, missing_fraction_retained_features = NA_real_, stringsAsFactors = FALSE)
  }
  excluded <- seed_miss$seed_id[is.finite(seed_miss$missing_fraction_retained_features) & seed_miss$missing_fraction_retained_features > missing_seed_max]
  module_counts <- table(meta$module[meta$retained_for_clustering %in% TRUE])
  weight_map <- setNames(1 / sqrt(as.numeric(module_counts)), names(module_counts))
  long_df$module_weight <- unname(weight_map[long_df$module])
  long_df$module_weight[is.na(long_df$module_weight)] <- NA_real_
  long_df$clustering_value <- long_df$scaled_value * long_df$module_weight
  if (length(retained_keys)) {
    center_map <- setNames(meta$center, meta$feature_key)
    scale_map <- setNames(meta$scale, meta$feature_key)
    module_map <- setNames(meta$module, meta$feature_key)
    weight_key_map <- setNames(as.numeric(weight_map[module_map[retained_keys]]), retained_keys)
    retained_mat <- trans_mat[, retained_keys, drop = FALSE]
    for (k in retained_keys) {
      miss <- !is.finite(retained_mat[, k])
      retained_mat[miss, k] <- center_map[[k]]
      retained_mat[, k] <- (retained_mat[, k] - center_map[[k]]) / scale_map[[k]]
      retained_mat[, k] <- retained_mat[, k] * as.numeric(weight_key_map[[k]])
    }
    keep_seed <- !(rownames(retained_mat) %in% excluded)
    retained_mat <- retained_mat[keep_seed, , drop = FALSE]
    wide <- data.frame(seed_id = rownames(retained_mat), retained_mat, check.names = FALSE)
    rownames(wide) <- NULL
  } else {
    wide <- data.frame()
  }
  list(long = long_df, wide = wide, metadata = meta, missingness = seed_miss, excluded_seeds = excluded)
}

o2ipa_long_to_wide <- function(df, feature_col, value_col) {
  if (!nrow(df)) return(data.frame())
  seeds <- unique(df$seed_id)
  feats <- unique(df[[feature_col]])
  mat <- matrix(NA_real_, nrow = length(seeds), ncol = length(feats), dimnames = list(seeds, feats))
  for (i in seq_len(nrow(df))) {
    mat[df$seed_id[[i]], df[[feature_col]][[i]]] <- as.numeric(df[[value_col]][[i]])
  }
  out <- data.frame(seed_id = rownames(mat), mat, check.names = FALSE)
  rownames(out) <- NULL
  out
}

o2ipa_distance <- function(wide) {
  if (!is.data.frame(wide) || nrow(wide) < 2L || ncol(wide) < 3L) return(NULL)
  mat <- as.matrix(wide[, setdiff(names(wide), "seed_id"), drop = FALSE])
  rownames(mat) <- wide$seed_id
  keep <- rowSums(!is.finite(mat)) == 0
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 1L) return(NULL)
  d <- as.matrix(stats::dist(mat, method = "euclidean"))
  diag(d) <- 0
  d
}

o2ipa_distance_long <- function(dmat) {
  if (is.null(dmat)) return(data.frame())
  idx <- which(upper.tri(dmat), arr.ind = TRUE)
  data.frame(
    seed_i = rownames(dmat)[idx[, 1]],
    seed_j = colnames(dmat)[idx[, 2]],
    distance = dmat[idx],
    stringsAsFactors = FALSE
  )
}

o2ipa_write_distance <- function(dmat, name, out_dir) {
  if (is.null(dmat)) {
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", paste0("distance_", name, "_matrix.tsv")))
    o2ipa_write_tsv(data.frame(), file.path(out_dir, "tables", paste0("distance_", name, "_long.tsv")))
    return(invisible(NULL))
  }
  mat_df <- data.frame(seed_id = rownames(dmat), dmat, check.names = FALSE)
  o2ipa_write_tsv(mat_df, file.path(out_dir, "tables", paste0("distance_", name, "_matrix.tsv")))
  o2ipa_write_tsv(o2ipa_distance_long(dmat), file.path(out_dir, "tables", paste0("distance_", name, "_long.tsv")))
}

o2ipa_upper_vec <- function(dmat) {
  if (is.null(dmat)) return(numeric(0))
  dmat[upper.tri(dmat)]
}

o2ipa_distance_correlations <- function(dlist) {
  pairs <- combn(names(dlist), 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    a <- dlist[[p[[1]]]]
    b <- dlist[[p[[2]]]]
    common <- intersect(rownames(a), rownames(b))
    if (length(common) < 3L) {
      rho <- NA_real_
      n <- 0L
    } else {
      aa <- a[common, common, drop = FALSE]
      bb <- b[common, common, drop = FALSE]
      av <- o2ipa_upper_vec(aa)
      bv <- o2ipa_upper_vec(bb)
      ok <- is.finite(av) & is.finite(bv)
      n <- sum(ok)
      rho <- if (n >= 3L) suppressWarnings(stats::cor(av[ok], bv[ok], method = "spearman")) else NA_real_
    }
    data.frame(distance_a = p[[1]], distance_b = p[[2]], n_pairs = n, spearman_rho = rho, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

o2ipa_silhouette_mean <- function(dmat, clusters) {
  n <- length(clusters)
  if (n < 3L || length(unique(clusters)) < 2L) return(NA_real_)
  vals <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    same <- clusters == clusters[[i]]
    other_clusters <- setdiff(unique(clusters), clusters[[i]])
    a <- if (sum(same) > 1L) mean(dmat[i, same & (seq_len(n) != i)]) else 0
    b <- min(vapply(other_clusters, function(cl) mean(dmat[i, clusters == cl]), numeric(1)))
    vals[[i]] <- if (max(a, b) > 0) (b - a) / max(a, b) else 0
  }
  mean(vals, na.rm = TRUE)
}

o2ipa_cluster_distance <- function(dmat, n_bootstrap = 200L, random_seed = 1L, feature_wide = NULL, feature_meta = NULL) {
  if (is.null(dmat) || nrow(dmat) < 3L) {
    return(list(diagnostics = data.frame(), membership = data.frame(), medoids = data.frame(), consensus = NULL, stability = data.frame(), recommended_k = NA_integer_, reason = "fewer_than_3_seeds"))
  }
  seeds <- rownames(dmat)
  hc <- stats::hclust(stats::as.dist(dmat), method = "ward.D2")
  kmax <- min(10L, nrow(dmat) - 1L)
  ks <- 2L:kmax
  diag_rows <- lapply(ks, function(k) {
    cl <- stats::cutree(hc, k = k)
    sizes <- as.integer(table(cl))
    data.frame(k = k, mean_silhouette = o2ipa_silhouette_mean(dmat, cl), min_cluster_size = min(sizes), max_cluster_size = max(sizes), cluster_size_imbalance = max(sizes) / min(sizes), stringsAsFactors = FALSE)
  })
  diag <- do.call(rbind, diag_rows)
  eligible <- diag[diag$min_cluster_size >= 2L & is.finite(diag$mean_silhouette), , drop = FALSE]
  membership <- do.call(rbind, lapply(ks, function(k) data.frame(seed_id = seeds, k = k, cluster = stats::cutree(hc, k = k), stringsAsFactors = FALSE)))
  if (!nrow(eligible)) {
    return(list(
      diagnostics = diag,
      membership = membership,
      medoids = data.frame(),
      consensus = NULL,
      stability = data.frame(seed_id = seeds, assignment_stability = NA_real_, stringsAsFactors = FALSE),
      recommended_k = NA_integer_,
      reason = "no_candidate_without_singleton_clusters",
      hclust = hc
    ))
  }
  recommended_k <- eligible$k[which.max(eligible$mean_silhouette)]
  cl_rec <- stats::cutree(hc, k = recommended_k)
  medoids <- o2ipa_medoids(dmat, cl_rec, recommended_k)

  consensus <- matrix(NA_real_, nrow = length(seeds), ncol = length(seeds), dimnames = list(seeds, seeds))
  stability <- data.frame(seed_id = seeds, assignment_stability = NA_real_, stringsAsFactors = FALSE)
  if (!is.null(feature_wide) && n_bootstrap > 0L && ncol(feature_wide) > 2L) {
    set.seed(random_seed)
    feature_mat <- as.matrix(feature_wide[match(seeds, feature_wide$seed_id), setdiff(names(feature_wide), "seed_id"), drop = FALSE])
    co <- matrix(0, nrow = length(seeds), ncol = length(seeds), dimnames = list(seeds, seeds))
    used <- 0L
    for (b in seq_len(n_bootstrap)) {
      cols <- sample(seq_len(ncol(feature_mat)), size = ncol(feature_mat), replace = TRUE)
      boot_mat <- feature_mat[, cols, drop = FALSE]
      if (ncol(boot_mat) < 1L || any(!is.finite(boot_mat))) next
      bd <- as.matrix(stats::dist(boot_mat))
      rownames(bd) <- colnames(bd) <- seeds
      bhc <- stats::hclust(stats::as.dist(bd), method = "ward.D2")
      bcl <- stats::cutree(bhc, k = recommended_k)
      co <- co + outer(bcl, bcl, "==")
      used <- used + 1L
    }
    if (used > 0L) {
      consensus <- co / used
      diag(consensus) <- 1
      stability$assignment_stability <- vapply(seq_along(seeds), function(i) {
        same <- cl_rec == cl_rec[[i]]
        if (sum(same) <= 1L) return(1)
        mean(consensus[i, same], na.rm = TRUE)
      }, numeric(1))
      diag$mean_within_cluster_consensus <- vapply(diag$k, function(k) {
        cl <- stats::cutree(hc, k = k)
        idx <- which(upper.tri(consensus), arr.ind = TRUE)
        same <- cl[idx[, 1]] == cl[idx[, 2]]
        if (!any(same)) NA_real_ else mean(consensus[cbind(idx[same, 1], idx[same, 2])], na.rm = TRUE)
      }, numeric(1))
    }
  }
  list(diagnostics = diag, membership = membership, medoids = medoids, consensus = consensus, stability = stability, recommended_k = recommended_k, reason = "ok", hclust = hc)
}

o2ipa_medoids <- function(dmat, clusters, k) {
  rows <- lapply(sort(unique(clusters)), function(cl) {
    idx <- which(clusters == cl)
    sub <- dmat[idx, idx, drop = FALSE]
    means <- if (length(idx) == 1L) 0 else rowMeans(sub)
    seed <- rownames(sub)[which.min(means)]
    data.frame(k = k, cluster = cl, medoid_seed = seed, mean_within_cluster_distance = min(means), cluster_size = length(idx), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

o2ipa_pairwise_diagnostics <- function(d_parameter, d_static, d_static18, d_dynamic, manifest, static_membership = NULL, dynamic_membership = NULL) {
  if (is.null(d_parameter) || is.null(d_static)) return(data.frame())
  common <- Reduce(intersect, Filter(Negate(is.null), list(rownames(d_parameter), rownames(d_static), if (!is.null(d_static18)) rownames(d_static18), if (!is.null(d_dynamic)) rownames(d_dynamic))))
  if (length(common) < 2L) return(data.frame())
  dp <- d_parameter[common, common, drop = FALSE]
  ds <- d_static[common, common, drop = FALSE]
  ds18 <- if (!is.null(d_static18)) d_static18[common, common, drop = FALSE] else matrix(NA_real_, length(common), length(common), dimnames = list(common, common))
  dd <- if (!is.null(d_dynamic)) d_dynamic[common, common, drop = FALSE] else matrix(NA_real_, length(common), length(common), dimnames = list(common, common))
  pv <- o2ipa_upper_vec(dp)
  sv <- o2ipa_upper_vec(ds)
  dv <- o2ipa_upper_vec(dd)
  th <- data.frame(
    distance = c("parameter", "static_process", "dynamic_phenotype"),
    close_threshold = c(stats::quantile(pv, 0.25, na.rm = TRUE), stats::quantile(sv, 0.25, na.rm = TRUE), stats::quantile(dv, 0.25, na.rm = TRUE)),
    far_threshold = c(stats::quantile(pv, 0.75, na.rm = TRUE), stats::quantile(sv, 0.75, na.rm = TRUE), stats::quantile(dv, 0.75, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  idx <- which(upper.tri(dp), arr.ind = TRUE)
  rows <- lapply(seq_len(nrow(idx)), function(ii) {
    i <- idx[ii, 1]
    j <- idx[ii, 2]
    si <- rownames(dp)[i]
    sj <- colnames(dp)[j]
    cls <- "intermediate"
    if (dp[i, j] >= th$far_threshold[th$distance == "parameter"] && ds[i, j] <= th$close_threshold[th$distance == "static_process"]) cls <- "parameter_far_process_close"
    if (ds[i, j] >= th$far_threshold[th$distance == "static_process"] && is.finite(dd[i, j]) && dd[i, j] <= th$close_threshold[th$distance == "dynamic_phenotype"]) cls <- "process_far_dynamic_close"
    if (ds[i, j] >= th$far_threshold[th$distance == "static_process"] && is.finite(dd[i, j]) && dd[i, j] >= th$far_threshold[th$distance == "dynamic_phenotype"]) cls <- "process_far_dynamic_far"
    if (dp[i, j] <= th$close_threshold[th$distance == "parameter"] && ds[i, j] >= th$far_threshold[th$distance == "static_process"]) cls <- "parameter_close_process_far"
    oi <- manifest$objective[match(si, manifest$seed_id)]
    oj <- manifest$objective[match(sj, manifest$seed_id)]
    data.frame(seed_i = si, seed_j = sj, parameter_distance = dp[i, j], static_process_distance = ds[i, j], static_18only_distance = ds18[i, j], dynamic_distance = dd[i, j], objective_i = oi, objective_j = oj, max_delta_objective = max(oi, oj, na.rm = TRUE) - min(manifest$objective, na.rm = TRUE), same_static_cluster = NA, same_dynamic_cluster = NA, diagnostic_class = cls, stringsAsFactors = FALSE)
  })
  list(pairs = do.call(rbind, rows), thresholds = th)
}

o2ipa_dynamic_fingerprints <- function(run_dir, seed_ids) {
  burden_path <- o2ipa_find_extra(run_dir, "predict1000_burden_total_seed_day_mean.tsv")
  ploidy_path <- o2ipa_find_extra(run_dir, "predict1000_ploidy_seed_day_mean.tsv")
  rows <- list()
  unavailable <- c("O2_min", "O2_terminal", "O2_AUC", "time_below_O2_2pct", "time_below_O2_1pct", "time_below_O2_0.5pct", "terminal_high_ploidy_fraction", "maximum_high_ploidy_fraction", "ploidy_entropy_AUC", "terminal_ploidy_entropy", "WGD_or_high_ploidy_lineage_fraction")
  if (!is.na(burden_path)) {
    tab <- o2ipa_read_tsv(burden_path)
    tab$seed_id <- o2ipa_norm_seed(tab$seed)
    tab <- tab[tab$seed_id %in% seed_ids, , drop = FALSE]
    for (key in split(seq_len(nrow(tab)), paste(tab$seed_id, tab$cohort, sep = "||"))) {
      d <- tab[key, , drop = FALSE]
      d <- d[order(d$day), , drop = FALSE]
      seed <- d$seed_id[[1]]
      cohort <- d$cohort[[1]]
      y <- as.numeric(d$burden_value)
      x <- as.numeric(d$day)
      prefix <- paste0("cohort_", cohort, "_")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "terminal_burden"), tail(y, 1), "positive", "existing_prediction")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "burden_AUC"), o2ipa_auc(x, y), "positive", "existing_prediction")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "maximum_burden"), max(y, na.rm = TRUE), "positive", "existing_prediction")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_burden", paste0(prefix, "time_to_burden_double"), o2ipa_crossing(x, y, 2 * y[[1]])$value, "time", "existing_prediction")
    }
  }
  if (!is.na(ploidy_path)) {
    tab <- o2ipa_read_tsv(ploidy_path)
    tab$seed_id <- o2ipa_norm_seed(tab$seed)
    tab <- tab[tab$seed_id %in% seed_ids, , drop = FALSE]
    for (key in split(seq_len(nrow(tab)), paste(tab$seed_id, tab$cohort, sep = "||"))) {
      d <- tab[key, , drop = FALSE]
      d <- d[order(d$day), , drop = FALSE]
      seed <- d$seed_id[[1]]
      cohort <- d$cohort[[1]]
      y <- as.numeric(d$ploidy_value)
      x <- as.numeric(d$day)
      prefix <- paste0("cohort_", cohort, "_")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_ploidy", paste0(prefix, "terminal_mean_chromosome"), tail(y, 1), "positive", "existing_prediction")
      rows[[length(rows) + 1L]] <- o2ipa_feature_row(seed, "dynamic", "dynamic_ploidy", paste0(prefix, "maximum_mean_chromosome"), max(y, na.rm = TRUE), "positive", "existing_prediction")
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(long = out, unavailable = data.frame(feature = unavailable, reason = "not_available_in_existing_prediction_outputs", stringsAsFactors = FALSE))
}

o2ipa_plot_outputs <- function(out_dir, d_parameter, d_static, d_static18, d_dynamic, static_scaled, dynamic_scaled, static_cluster, manifest, static_long) {
  fig <- function(name) file.path(out_dir, "figures", paste0(name, ".pdf"))
  scatter <- function(a, b, xlab, ylab, name) {
    if (is.null(a) || is.null(b)) return()
    common <- intersect(rownames(a), rownames(b))
    if (length(common) < 3L) return()
    av <- o2ipa_upper_vec(a[common, common, drop = FALSE])
    bv <- o2ipa_upper_vec(b[common, common, drop = FALSE])
    grDevices::pdf(fig(name), width = 6, height = 5)
    plot(av, bv, pch = 16, xlab = xlab, ylab = ylab, main = name)
    grDevices::dev.off()
  }
  scatter(d_parameter, d_static, "Parameter distance", "Static process distance", "parameter_vs_static_process_distance")
  scatter(d_parameter, d_static18, "Parameter distance", "Static 18-only distance", "parameter_vs_static_18only_distance")
  scatter(d_static, d_dynamic, "Static process distance", "Dynamic phenotype distance", "static_process_vs_dynamic_phenotype_distance")
  scatter(d_static, d_static18, "Static full distance", "Static 18-only distance", "static_full_vs_static_18only_distance")

  if (is.data.frame(static_scaled) && nrow(static_scaled) >= 2L && ncol(static_scaled) >= 3L) {
    mat <- as.matrix(static_scaled[, setdiff(names(static_scaled), "seed_id"), drop = FALSE])
    rownames(mat) <- static_scaled$seed_id
    mat <- mat[rowSums(!is.finite(mat)) == 0, colSums(!is.finite(mat)) == 0, drop = FALSE]
    if (nrow(mat) >= 2L && ncol(mat) >= 2L) {
      grDevices::pdf(fig("static_process_feature_heatmap"), width = 9, height = 7)
      heatmap(mat, scale = "none", margins = c(5, 5), main = "Static process features")
      grDevices::dev.off()
      grDevices::pdf(fig("static_process_hierarchical_dendrogram"), width = 7, height = 5)
      plot(stats::hclust(stats::dist(mat), method = "ward.D2"), main = "Static process dendrogram", xlab = "")
      grDevices::dev.off()
      pc <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
      grDevices::pdf(fig("static_process_pca"), width = 6, height = 5)
      plot(pc$x[, 1], pc$x[, 2], pch = 16, xlab = "PC1", ylab = "PC2", main = "Static process PCA")
      text(pc$x[, 1], pc$x[, 2], labels = rownames(mat), pos = 3, cex = 0.6)
      grDevices::dev.off()
    }
  }
  if (is.data.frame(dynamic_scaled) && nrow(dynamic_scaled) >= 2L && ncol(dynamic_scaled) >= 3L) {
    mat <- as.matrix(dynamic_scaled[, setdiff(names(dynamic_scaled), "seed_id"), drop = FALSE])
    rownames(mat) <- dynamic_scaled$seed_id
    mat <- mat[rowSums(!is.finite(mat)) == 0, colSums(!is.finite(mat)) == 0, drop = FALSE]
    if (nrow(mat) >= 2L && ncol(mat) >= 2L) {
      grDevices::pdf(fig("dynamic_phenotype_heatmap"), width = 8, height = 6)
      heatmap(mat, scale = "none", margins = c(5, 5), main = "Dynamic phenotype features")
      grDevices::dev.off()
      pc <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
      grDevices::pdf(fig("dynamic_phenotype_pca"), width = 6, height = 5)
      plot(pc$x[, 1], pc$x[, 2], pch = 16, xlab = "PC1", ylab = "PC2", main = "Dynamic phenotype PCA")
      grDevices::dev.off()
    }
  }
  if (!is.null(static_cluster$consensus) && all(is.finite(static_cluster$consensus))) {
    grDevices::pdf(fig("bootstrap_consensus_matrix"), width = 7, height = 6)
    image(seq_len(nrow(static_cluster$consensus)), seq_len(ncol(static_cluster$consensus)), static_cluster$consensus, xlab = "", ylab = "", main = "Bootstrap consensus")
    grDevices::dev.off()
  }
  if (nrow(static_cluster$membership)) {
    rec <- static_cluster$recommended_k
    memb <- static_cluster$membership[static_cluster$membership$k == rec, , drop = FALSE]
    mf <- merge(manifest, memb, by = "seed_id", all.x = FALSE)
    if (nrow(mf)) {
      grDevices::pdf(fig("objective_by_cluster"), width = 6, height = 5)
      boxplot(objective ~ cluster, data = mf, xlab = "Cluster", ylab = "Objective", main = "Objective by cluster")
      grDevices::dev.off()
      grDevices::pdf(fig("boundary_risk_by_cluster"), width = 6, height = 5)
      barplot(tapply(as.numeric(mf$boundary_risk), mf$cluster, mean, na.rm = TRUE), ylab = "Boundary risk fraction", xlab = "Cluster", main = "Boundary risk by cluster")
      grDevices::dev.off()
    }
  }
  if (nrow(static_long) && nrow(static_cluster$medoids)) {
    med <- static_cluster$medoids$medoid_seed
    curve_df <- static_long[static_long$seed_id %in% med & static_long$module %in% c("hypoxia_sensing", "proliferation", "death", "CIN_missegregation") & grepl("_O2_|stress_at_O2", static_long$feature), , drop = FALSE]
    if (nrow(curve_df)) {
      grDevices::pdf(fig("cluster_medoids_response_curves"), width = 8, height = 6)
      par(mfrow = c(2, 2))
      for (mod in unique(curve_df$module)[seq_len(min(4L, length(unique(curve_df$module))))]) {
        dd <- curve_df[curve_df$module == mod, , drop = FALSE]
        stripchart(raw_value ~ seed_id, data = dd, vertical = TRUE, method = "jitter", pch = 16, main = mod, ylab = "value")
      }
      grDevices::dev.off()
    }
  }
}

o2ipa_report <- function(out_dir, run_dir, manifest, param_meta, feature_meta, dist_cor, cluster, unavailable, limitations) {
  git_branch <- tryCatch(system("git branch --show-current", intern = TRUE), error = function(e) NA_character_)
  git_hash <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA_character_)
  lines <- c(
    "# In vivo process analysis summary",
    "",
    paste0("- Input directory: `", run_dir, "`"),
    paste0("- Output directory: `", out_dir, "`"),
    paste0("- Git branch: `", paste(git_branch, collapse = " "), "`"),
    paste0("- Git commit: `", paste(git_hash, collapse = " "), "`"),
    paste0("- Analysis time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- Discovered seeds: ", nrow(manifest)),
    paste0("- Valid seeds: ", sum(manifest$fit_success %in% TRUE)),
    paste0("- Objective source(s): ", paste(unique(manifest$objective_source), collapse = ", ")),
    paste0("- Final static feature count: ", sum(feature_meta$retained_for_clustering %in% TRUE, na.rm = TRUE)),
    paste0("- Recommended static cluster k: ", ifelse(is.na(cluster$recommended_k), "not available", cluster$recommended_k)),
    "",
    "## Scientific interpretation boundary",
    "",
    "Process clustering represents data-compatible mechanistic regimes; it does not automatically equal true biological heterogeneity. Only differences that exceed fitting uncertainty and predict independent data should be interpreted as real biological mechanisms.",
    "",
    "## Distance correlations",
    ""
  )
  if (nrow(dist_cor)) {
    lines <- c(lines, utils::capture.output(print(dist_cor, row.names = FALSE)))
  } else {
    lines <- c(lines, "No distance correlations available.")
  }
  lines <- c(lines, "", "## Cluster medoids", "")
  if (nrow(cluster$medoids)) {
    lines <- c(lines, utils::capture.output(print(cluster$medoids, row.names = FALSE)))
  } else {
    lines <- c(lines, "No stable cluster medoids available.")
  }
  lines <- c(lines, "", "## Parameter transforms", "")
  lines <- c(lines, utils::capture.output(print(param_meta[, c("parameter", "module", "transform", "optimizer_scale")], row.names = FALSE)))
  lines <- c(lines, "", "## Unavailable dynamic features", "")
  lines <- c(lines, if (nrow(unavailable)) paste0("- ", unavailable$feature, ": ", unavailable$reason) else "None recorded.")
  lines <- c(lines, "", "## Limitations", "", paste0("- ", limitations))
  dir.create(file.path(out_dir, "report"), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(out_dir, "report", "analysis_summary.md"))
}

o2ipa_make_mock_run <- function(run_dir, n_seed = 5L) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  params <- o2ipa_target_params()
  base <- c(
    O2_crit = 0.2, alpha_o2 = 1.2, gamma_growth = 3, mu_hp = 0.05, p_misseg = 0.1,
    k_o_mis = 0.003, buffer_smax = 0.9, buffer_beta = 0.2, buffer_n_exp = 4,
    n_O = 1.5, lam_max = 0.22, p_mis_base = 1e-4, p_wgd = 1e-4, gamma_mu = 2.5,
    o2_S0 = 3.5, kappa_O = 0.8, eta_o2 = 1.1, rho_2N = 1e5
  )
  set.seed(123)
  for (i in seq_len(n_seed)) {
    sdir <- file.path(run_dir, paste0("seed", i))
    dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
    mult <- exp(stats::rnorm(length(params), 0, 0.12))
    vals <- base[params] * mult
    vals[c("gamma_growth", "gamma_mu", "n_O", "buffer_smax")] <- base[c("gamma_growth", "gamma_mu", "n_O", "buffer_smax")] + stats::rnorm(4, 0, 0.05)
    vals["buffer_smax"] <- min(max(vals["buffer_smax"], 0.05), 0.99)
    best <- data.frame(parameter = names(vals), value = as.numeric(vals), stringsAsFactors = FALSE)
    best <- rbind(best, data.frame(parameter = c("o2_min", "o2_S0_upper_bound", "o2_Nref"), value = c(0, 5, 1e6)))
    o2ipa_write_tsv(best, file.path(sdir, "best_params.tsv"))
    fs <- data.frame(metric = c("objective", "objective_data", "objective_burden", "objective_ploidy", "objective_burden_neg2loglik_raw", "objective_ploidy_neg2loglik_raw", "deoptim_stop_reason"), value = c(10 + i / 10, 9 + i / 10, 4 + i / 20, 5 + i / 20, 1000 + i, 2000 + i, "mock_converged"))
    o2ipa_write_tsv(fs, file.path(sdir, "fit_summary.tsv"))
    pt <- data.frame(param_name = paste0(ifelse(params %in% c("gamma_growth", "gamma_mu", "n_O", "buffer_smax"), "", "log10_"), params), estimate = TRUE, init_value = 0, lower_bound = -10, upper_bound = 10, param_prototype = params, prototype_init_value = base[params], prototype_lower_bound = base[params] / 10, prototype_upper_bound = base[params] * 10, source = "mock", note = "mock", stringsAsFactors = FALSE)
    utils::write.csv(pt, file.path(sdir, "parameter_table.csv"), row.names = FALSE)
  }
  extra <- file.path(run_dir, "extra_results")
  dir.create(extra, recursive = TRUE, showWarnings = FALSE)
  days <- 0:20
  burden <- do.call(rbind, lapply(seq_len(n_seed), function(i) {
    do.call(rbind, lapply(c("2N", "4N"), function(cohort) data.frame(seed = paste0("seed", i), cohort = cohort, day = days, burden_value = (10 + i) * exp(0.05 * days) * ifelse(cohort == "4N", 1.1, 1), stringsAsFactors = FALSE)))
  }))
  ploidy <- do.call(rbind, lapply(seq_len(n_seed), function(i) {
    do.call(rbind, lapply(c("2N", "4N"), function(cohort) data.frame(seed = paste0("seed", i), cohort = cohort, day = days, ploidy_value = ifelse(cohort == "4N", 88, 44) + i * 0.2 + sin(days / 5), stringsAsFactors = FALSE)))
  }))
  o2ipa_write_tsv(burden, file.path(extra, "predict1000_burden_total_seed_day_mean.tsv"))
  o2ipa_write_tsv(ploidy, file.path(extra, "predict1000_ploidy_seed_day_mean.tsv"))
  invisible(run_dir)
}
